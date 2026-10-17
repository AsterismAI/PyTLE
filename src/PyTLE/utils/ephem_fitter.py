from datetime import timedelta, datetime
import numpy as np
import scipy.optimize

from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs72
from sgp4.propagation import sgp4 as sgprop

import PyTLE
from PyTLE import tle_fitter
from PyTLE.utils import julian

# ----------------------------------------------------------------------------------------------------- 
def propTLE( L1, L2, offsets_m ):
    # tle object
    tleo = twoline2rv( L1, L2, wgs72 )
    eph = np.vstack( [ np.hstack( sgprop(tleo,X) ) for X in offsets_m ] )
    return eph

# ----------------------------------------------------------------------------------------------------- 
def propTLEtoJD( L1, L2, JD ):
    # tle object
    tleo  = twoline2rv( L1, L2, wgs72 )
    off_m = (JD - tleo.jdsatepoch) * 1440 
    eph   = np.vstack( [ np.hstack( sgprop(tleo,X) ) for X in off_m ] )
    eph   = np.hstack( (JD[:,np.newaxis], eph ) )
    return eph

# ----------------------------------------------------------------------------------------------------- 
def opt( X, fitR, offsets_m, truth_eph, fields, verbose=True ):
    fitR.from_fields( fields, X )
    try : 
        L1, L2 = fitR.generateLines()
        eph = propTLE( L1, L2, offsets_m )
    except: return np.inf
    diff = eph[:,0:3] - truth_eph[:,0:3]
    diff = np.linalg.norm( diff, axis=1 )
    rms  = np.sqrt( np.sum( diff**2 ) / len(diff) )
    if verbose:
        print('Points : {}  RMS: {:015.8f} RMS per point: {:015.8f}          '.format( 
                len(diff), rms, rms/len(diff) ) ,end='\r') 
    return rms

# ----------------------------------------------------------------------------------------------------- 
def fit_ephemeris( ephemeris : np.array, init_tle : tle_fitter, fields = None ):
    L1,L2  = init_tle.generateLines()
    new    = tle_fitter( PyTLE.TLE.parseLines( L1, L2 ) )
    tle_jd = julian.datetime_to_julian( new.epoch )
    eph_jd = ephemeris[:,0]
    eph_pv = ephemeris[:,1:]
    off_m  = 1440 * ( eph_jd - tle_jd )
    if fields is None:
        fields = new.get_available_field_names()
    initX  = new.get_fields( fields ) 
    initS  = initX * np.random.normal( 1, 0.05, size=( len(initX) + 1, len(initX) ) )
    initS[0,:] = initX
    ans = scipy.optimize.minimize( opt,
                initX,
                args = (new, off_m, eph_pv, fields ),
                method = 'Nelder-Mead',
                options = {'fatol' : 0.1, 'xatol' : 10 }) # 'initial_simplex' : initS })
    if ans.success : 
        return new.from_fields( fields, ans.x ) 
    return None

# ----------------------------------------------------------------------------------------------------- 
def reepochTLE( init_tle    : tle_fitter, 
                new_epoch   : datetime,
                eph_spacing : int = 5,
                days_around : int = 7): 
    assert isinstance( init_tle, tle_fitter )
    new_jd    = julian.datetime_to_julian( new_epoch ) + np.arange( -days_around, days_around+0.0001, eph_spacing/1440 )
    L1, L2    = init_tle.generateLines()
    new_eph   = propTLEtoJD( L1, L2, new_jd)
    middle    = new_eph[ new_eph.shape[0] // 2, :]
    jd, P, V  = middle[0], middle[1:4], middle[4:]
    temp_tle  = init_tle.fromPV( julian.julian_to_datetime( jd ), P, V )
    print(init_tle.generateLines())
    init_tle  = init_tle.from_array( temp_tle.to_array() )
    print(init_tle.generateLines())
    # this is a bad idea; we're just setting the new epoch and not modifying any other orbit elements
    # for now though, just brute-force it
    
    return fit_ephemeris( new_eph, init_tle, None )
  
# -----------------------------------------------------------------------------------------------------   
def test_reepoch():
    L1 = '1 43556U 18046C   22321.55519027  .00025005  00000+0  49749-3 0  9993'
    L2 = '2 43556  51.6329 154.1269 0008144 222.8163 137.2191 15.46745497242947'
    init_tle  = tle_fitter( PyTLE.TLE.parseLines( L1, L2 ) )
    newepoch  = datetime(year=2023, month=1, day=20 )
    newtle = reepochTLE( init_tle, newepoch, eph_spacing=5, days_around=10 )
    newtle.satno = 99999
    print('Re-epoch test, re-epoching :')
    print('\n'.join( init_tle.generateLines() ) )
    print('To : {}'.format( newepoch ))
    print('\n'.join( newtle.generateLines() ) )
    print()

# -----------------------------------------------------------------------------------------------------   
def test_fitephem():
    L1 = '1 43556U 18046C   22321.55519027  .00025005  00000+0  49749-3 0  9993'
    L2 = '2 43556  51.6329 154.1269 0008144 222.8163 137.2191 15.46745497242947'
    init_tle  = PyTLE.TLE.parseLines( L1, L2 )
    new_jd    = julian.datetime_to_julian( init_tle.epoch ) + np.arange( 90, 90+14, 5/1440 )
    new_eph   = propTLEtoJD( L1, L2, new_jd)
    new_tle   = PyTLE.TLE.fromPV( julian.julian_to_datetime(new_eph[0,0]),
                                 new_eph[0,1:4], new_eph[0,4:])
    new_tle = fit_ephemeris( new_eph, new_tle,['mean_motion','mean_anomaly'] )
    print('Fit 14 day span over 5 day intervals; moved epoch 90 days')
    print(L1)
    print(L2)
    print(new_tle)
    

# ----------------------------------------------------------------------------------------------------- 
if __name__ == '__main__':
    test_reepoch()
    # test_fitephem()

