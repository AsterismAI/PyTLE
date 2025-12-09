from datetime import timedelta, datetime
import numpy as np
from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs72
from sgp4.propagation import sgp4 as sgprop
import scipy.optimize

from PyTLE import tle_fitter
from PyTLE import TLE


# ----------------------------------------------------------------------------------------------------- 
def propTLE( L1, L2, offsets_m ):
    # tle object
    tleo = twoline2rv( L1, L2, wgs72 )
    eph = np.vstack( [ np.hstack( sgprop(tleo,X) ) for X in offsets_m ] )
    return eph

# ----------------------------------------------------------------------------------------------------- 
def opt( X, fitR, offsets_m, truth_eph ):
    fitR.from_array( X )
    try : 
        L1, L2 = fitR.generateLines()
        eph = propTLE( L1, L2, offsets_m )
    except: return np.inf
    diff = eph[:,0:3] - truth_eph[:,0:3]
    diff = np.linalg.norm( diff, axis=1 )
    rms  = np.sqrt( np.sum( diff**2 ) / len(diff) )
    print('{:015.10f}                        '.format( rms ) ,end='\r') 
    return rms

# ----------------------------------------------------------------------------------------------------- 
def reepochTLE( L1 : str, L2 : str, 
               newepoch : datetime,
               SPAN = 10, 
               SPACING = 5 ):
    # parse the TLE 
    old = tle_fitter( TLE.parseLines( L1, L2 ) )
    new = tle_fitter( TLE.parseLines( L1, L2 ) )

    # NOTE: we have two different epochs here
    #   (1) is from the original TLE to the interval that we will fit
    #   (2) is the offset from the NEW TLE to those data that we will fit
    #   we precompute them for efficiency (or we can stuff dates in here)

    # with the new epoch set, find the minutes offset so that we can propagate
    offset_mid   = ( newepoch - old.epoch ).total_seconds() / 60
    offset_m     = np.arange( offset_mid - (1440*(SPACING/2)),
                              offset_mid + (1440*(SPACING/2)),
                              SPACING )

    # this is for the *new* TLE
    new_offset_m = offset_m - offset_mid

    # truth ephemeris for the new interval
    truth_eph = propTLE( L1, L2, offset_m )
    
    # set new values
    new.epoch = newepoch
    new.satno = 99998
    ans = scipy.optimize.minimize( opt,
                                  new.to_array(),
                                  args = (new, new_offset_m, truth_eph ),
                                  method = 'Nelder-Mead' )
    if ans.success : 
        new = new.from_array( ans.x ) 
        return new

    return None


# ----------------------------------------------------------------------------------------------------- 
if __name__ == '__main__':
    L1 = '1 43556U 18046C   22321.55519027  .00025005  00000+0  49749-3 0  9993'
    L2 = '2 43556  51.6329 154.1269 0008144 222.8163 137.2191 15.46745497242947'

    #L1 = '1 25544U 98067A   25343.50984768  .00014141  00000-0  26026-3 0  9994'
    #L2 = '2 25544  51.6310 158.4724 0003400 223.8185 136.2534 15.49447168542357'

    X = reepochTLE( L1, L2, datetime(year=2023, month=2, day=1 ) )
    print('OLD')
    print('\n'.join([L1,L2]))
    print('NEW')
    print('\n'.join( X.generateLines() ) )

