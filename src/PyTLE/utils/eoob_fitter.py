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


# ---------------------------------------------------------------------------------------
def spherical_to_xyz( ra, dec ):
    tx = np.cos( dec ) * np.cos( ra )
    ty = np.cos( dec ) * np.sin( ra )
    tz = np.sin( dec )
    return np.hstack( ( tx[:,np.newaxis], ty[:,np.newaxis], tz[:,np.newaxis] ) )

# ---------------------------------------------------------------------------------------
def udleoob_setup( data : pd.DataFrame ):
    obtime = pd.to_datetime( data['obTime'] )
    atime  = astropy.time.Time( obtime )
    sen_lla = astropy.coordinates.EarthLocation.from_geodetic(
        lat    = data['senlat'] * u.deg,
        lon    = data['senlon'] * u.deg,
        height = data['senalt'] * u.km
    )
    itrs = sen_lla.get_itrs( obstime=atime )
    gcrs = itrs.transform_to( astropy.coordinates.TEME( obstime=atime ) )
    gcrs_p = gcrs.cartesian.xyz.to_value(u.km).T
    data['teme_x'] = gcrs_p[:,0]
    data['teme_y'] = gcrs_p[:,1]
    data['teme_z'] = gcrs_p[:,2]
    
    # calculate the look vector in j2k
    gcrslv = spherical_to_xyz( np.radians( data['ra'].values ),
                           np.radians( data['declination'].values ) )
    # and convert this to TEME
    gcrs = astropy.coordinates.GCRS( 
                                    obstime = atime,
                                    x = gcrslv[:,0],
                                    y = gcrslv[:,1],
                                    z = gcrslv[:,2],
                                    representation_type='cartesian')
    teme = gcrs.transform_to( astropy.coordinates.TEME( obstime=atime ) )
    teme_lv = teme.cartesian.xyz.to_value(u.one).T
    data['teme_lv'] = teme_lv.tolist()
    
    return data
    
# =================================================================================================
if __name__ == '__main__':
    D = pd.read_json('~/Downloads/39070.json')
    D = udleoob_setup( D )
    print(D)