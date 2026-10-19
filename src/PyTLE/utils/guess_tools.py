from datetime import timedelta, datetime
import numpy as np

from PyTLE import tle_fitter
from PyTLE.utils import julian


from PyTLE.utils.sgp4_helpers import propTLE, propTLEtoJD

# ----------------------------------------------------------------------------------------------------- 
def move_epoch_rough( init_tle : tle_fitter, new_epoch : datetime ):
    # roughly move this TLE to that new time
    off_m     = (new_epoch - init_tle.epoch).total_seconds() / 60
    eph       = propTLE( *init_tle.generateLines(), [off_m])[0]
    return init_tle.fromPV( new_epoch, eph[0:3], eph[3:] )