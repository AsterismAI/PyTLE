import numpy as np

from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs72
from sgp4.propagation import sgp4 as sgprop

from PyTLE import tle_fitter
from PyTLE.utils import julian

# ----------------------------------------------------------------------------------------------------- 
def propTLE( L1, L2, offsets_m ):
    '''
    this outputs:
    x,y,z,dx,dy,z
    '''
    tleo = twoline2rv( L1, L2, wgs72 )
    eph = np.vstack( [ np.hstack( sgprop(tleo,X) ) for X in offsets_m ] )
    return eph

# ----------------------------------------------------------------------------------------------------- 
def propTLEtoJD( L1, L2, JD ):
    '''
    this outputs:
    jd,x,y,z,dx,dy,z
    '''
    tleo  = twoline2rv( L1, L2, wgs72 )
    off_m = (JD - tleo.jdsatepoch) * 1440 
    eph   = np.vstack( [ np.hstack( sgprop(tleo,X) ) for X in off_m ] )
    eph   = np.hstack( (JD[:,np.newaxis], eph ) )
    return eph