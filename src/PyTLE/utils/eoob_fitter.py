from datetime import timedelta, datetime
import numpy as np
import pandas as pd
import scipy.optimize

from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs72
from sgp4.propagation import sgp4 as sgprop

import PyTLE
from PyTLE import tle_fitter
from PyTLE.utils import julian

import astropy.coordinates
import astropy.units as u
import astropy.time

from PyTLE.utils.sgp4_helpers import propTLE, propTLEtoJD
'''
code here will take EO observations that have the following fields:
    - `teme_lv` : TEME look vector of the observation
    - `teme_p`  : TEME position vector of the sensor
'''
# ---------------------------------------------------------------------------------------
def spherical_to_xyz( ra, dec, radius=1 ):
    tx = radius * np.cos( dec ) * np.cos( ra )
    ty = radius * np.cos( dec ) * np.sin( ra )
    tz = radius * np.sin( dec )
    return np.hstack( ( tx[:,np.newaxis], ty[:,np.newaxis], tz[:,np.newaxis] ) )

# ---------------------------------------------------------------------------------------
def setup_udleoob( data : pd.DataFrame ):
    '''
    this takes in unclass UDL EO observations and sets them up with 
    TEME position and look vectors for solving; since we iterate within
    the TLE, it makes more sense to conver the truth data to that format
    so that we aren't constantly converting (and costing cycles)
    '''
    obtime         = pd.to_datetime( data['obTime'] )
    data['obTime'] = obtime.dt.tz_convert(None)
    data           = data.sort_values(by='obTime')
    atime          = astropy.time.Time( obtime )
    sen_lla        = astropy.coordinates.EarthLocation.from_geodetic(
                    lat    = data['senlat'] * u.deg,
                    lon    = data['senlon'] * u.deg,
                    height = data['senalt'] * u.km
    )
    sen_gcrs       = sen_lla.get_gcrs( obstime=atime )
    # the obs from UDL are in J2K; consider those to be GCRS as the root frame
    gcrs_ap = astropy.coordinates.GCRS( 
                                      ra        = data['ra'].values * u.deg,
                                      dec       = data['declination'].values * u.deg ,
                                      obsgeoloc = sen_gcrs.cartesian )
    # convert those ob-vectors into TEME and store the look vector
    teme_ap = gcrs_ap.transform_to( astropy.coordinates.TEME( obstime=atime ) )
    teme_lv = teme_ap.cartesian.xyz.to_value(u.one).T
    # get the sensor location in GCRS so we can convert to TEME (we need this later)
    gcrs    = sen_lla.get_gcrs( obstime=atime )
    teme    = gcrs.transform_to( astropy.coordinates.TEME( obstime=atime ) )
    teme_p  = teme.cartesian.xyz.to_value(u.km).T
    data['teme_p'] = teme_p.tolist()
    
    # ///////////////////////////////////////////////////////////////////////////////
    # version 1 of conversion
    # calculate the look vector in j2k
    gcrslv = spherical_to_xyz( np.radians( data['ra'].values ),
                               np.radians( data['declination'].values ),
                               radius = 1e6)
    
    # # and convert this to TEME
    gcrs = astropy.coordinates.GCRS( 
                                    obstime = atime,
                                    x = gcrslv[:,0],
                                    y = gcrslv[:,1],
                                    z = gcrslv[:,2],
                                    representation_type='cartesian')
    teme = gcrs.transform_to( astropy.coordinates.TEME( obstime=atime ) )
    Zteme_lv = teme.cartesian.xyz.to_value(u.one).T
    # make sure that these vectors are normalized post-transform
    Zteme_lv = Zteme_lv / np.linalg.norm( Zteme_lv, axis=1)[:,np.newaxis]
    for A,B in zip( teme_lv, Zteme_lv):
        print(A)
        print(B)
        print(A-B)
        print(np.arccos( np.sum(A*B) ) * 3600 )
        print()
    assert 1 == 2
    # //////////////////////////////////////////////////////////////////////////////////
    data['teme_lv'] = teme_lv.tolist()
    return data
    
# ---------------------------------------------------------------------------------------
def fit_TLE_to_obs( init_tle : tle_fitter,
                    obs_df   : pd.DataFrame,
                    fields   : list[str] = None):
    sensor_p  = np.vstack( obs_df['teme_p'] )
    sensor_lv = np.vstack( obs_df['teme_lv'] )
    obs_jd    = np.array( [julian.datetime_to_julian(X) for X in obs_df['obTime'] ])
    # assume the TLE is already epoch'd to the proper place (we could force last ob here)
    offset_m  = 1440 * ( obs_jd - julian.datetime_to_julian(init_tle.epoch) )
    if fields is None: fields = init_tle.get_available_field_names()
    initX  = init_tle.get_fields( fields ) 
    initS  = np.ones( shape = (len(initX) + 1, len(initX) ) ) * initX 
    for i in range(1, len(initX) ):
        initS[i,i-1] += 0.08
        initS[i,i]   -= 0.08
    
    ans = scipy.optimize.minimize( opt,
                initX,
                args    = (init_tle, offset_m, sensor_p, sensor_lv, fields ),
                method  = 'Nelder-Mead',
                options = {'fatol' : 5, 'xatol' : 10 }) # ,'initial_simplex' : initS })
    if ans.success : 
        return init_tle.from_fields( fields, ans.x ) 
    return None
    
# ----------------------------------------------------------------------------------------------------- 
def opt( X, fitR, offsets_m, senpos, actlv, fields, verbose=True ):
    fitR.from_fields( fields, X )
    try : 
        L1, L2 = fitR.generateLines()
        eph = propTLE( L1, L2, offsets_m )
    except: return np.inf
    # this return does NOT return JD as the first column... built for speed so 
    # we don't re-calculate TLE offset every call
    test_lv = eph[:,0:3] - senpos
    test_lv = test_lv / np.linalg.norm( test_lv, axis=1 )[:,np.newaxis]
    dotp    = np.sum( test_lv * actlv, axis=1 )
    ang     = np.degrees( np.arccos( dotp ) ) * 3600
    rms     = np.sqrt( np.sum( ang**2 ) / len(ang) )
    if verbose:
        print('Points : {}  RMS: {:015.8f} RMS-per-point: {:015.8f}          '.format( 
                len(ang), rms, rms/len(ang) ) ,end='\r') 
    return rms


def test_goes():
    import copy
    from PyTLE.utils.guess_tools import move_epoch_rough

    obs_df = pd.read_json('~/Downloads/60133.json')
    obs_df = setup_udleoob( obs_df )
    
    L1 = '1 60133U 24119A   26291.41659517 -.00000246  00000-0  00000-0 0  9991'
    L2 = '2 60133   0.0121  93.5376 0000410 180.3882 187.8035  1.00273432  8201'
    TLE = PyTLE.tle_fitter().parseLines( L1 , L2 )
    # set the new epoch position
    old_epoch = copy.copy( TLE.epoch )
    new_epoch = max( obs_df['obTime'] )
    print('Setting epoch from : {} to : {}'.format( old_epoch, new_epoch ) )
    print(TLE)
    TLE = move_epoch_rough( TLE, new_epoch )
    print()
    TLE = fit_TLE_to_obs( TLE, obs_df, fields=['mean_motion', 'eccentricity', 'inclination', 'argp', 'raan', 'mean_anomaly'] )
    print()
    print()
    TLE.satno = 99999
    print(TLE)
    
# =================================================================================================
if __name__ == '__main__':
    test_goes()
    # from PyTLE.utils.guess_tools import move_epoch_rough
    # # load up some UDL obs and prepare them
    # obs_df = pd.read_json('~/Downloads/39070.json')
    # obs_df = setup_udleoob( obs_df )
    
    # # use this as our initial TLE
    # L1  = '1 39070U 13004A   26291.18317172  .00000084  00000-0  00000-0 0  9999'
    # L2  = '2 39070   3.7116  33.0351 0002359  48.3828 197.0572  1.00264806 47926'
    # TLE = PyTLE.tle_fitter().parseLines( L1 , L2 )
    
    # new_epoch = max( obs_df['obTime'] )
    # TLE = move_epoch_rough( TLE, new_epoch )
    
    # print(TLE)
    # TLE = fit_TLE_to_obs( TLE, obs_df)
    # print()
    # TLE.satno = 99999
    # print(TLE)
    