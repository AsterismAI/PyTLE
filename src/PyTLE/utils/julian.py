import datetime
import numpy as np

_EPOCH_JD = datetime.datetime(2000, 1, 1, 12, 0, 0)
# -----------------------------------------------------------------------------------------------------
def datetime_to_julian(dt: datetime.datetime) -> float:
    """Convert a timezone-naive UTC datetime object to a Julian Date."""
    # Epoch anchor: JD 2451545.0 corresponds to 2000-01-01 12:00:00 UTC
    delta = dt - _EPOCH_JD
    return 2451545.0 + (delta.total_seconds() / 86400.0)

# -----------------------------------------------------------------------------------------------------
def julian_to_datetime(jd: float) -> datetime.datetime:
    """Convert a Julian Date float to a timezone-naive UTC datetime object."""
    delta_days = jd - 2451545.0
    return _EPOCH_JD + datetime.timedelta(days=delta_days)

# =====================================================================================================
if __name__ == '__main__':
    import pandas as pd
    N           = 10000
    NOW         = datetime.datetime.now()
    in_dates    = np.array( [ NOW + datetime.timedelta( days=X ) for X in np.random.uniform( 0, 365*100, N ) ] )
    print('Actual, Computed Julian, Pandas Julian, Re-convert')
    for D in in_dates:
        computed = datetime_to_julian(D)
        print('{} {} {} {} {}'.format(
            D,
            computed,
            pd.to_datetime(D).to_julian_date(),
            julian_to_datetime(computed),
            (D - julian_to_datetime(computed) ).total_seconds()
            )
        )

