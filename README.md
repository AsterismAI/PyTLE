# PyTLE : convenience routines for TLE parsing and fitting

Kerry N. Wood (kerry.n.wood@gmail.com)

#### PyTLE.TLE
<<<<<<< HEAD
- routines to parse, store, and re-generate two line TLE 
=======
- routines to parse, store, and re-generate two line TLE  (type 0 and 4)
>>>>>>> packaging
- data fields are stored in their native units (e.g. degrees for RAAN, inclination, etc)
- convenience routines for initializing *new* TLE
	* `fromCOE` : from classical osculating elements
	* `fromPV`  : from state vectors (in native frame and units / TEME / km / km/s)
<<<<<<< HEAD
	* `fromPoliastro` : from Poliastro orbit object
=======
>>>>>>> packaging

#### PyTLE.tle_fitter
- wraps `PyTLE.TLE` and maps TLE fields to ranges useful for optimization 
- example `ephem_fit` function will fit a TLE (depends on Brandon Rhode's SGP4 code) to an ephemeris frame (in TEME)

## Credits:
- alpha routines borrowed and modified from Brandon Rhodes SGP4 library
-`julian.py` is taken from Daniel Zawada's `julian.py` code; it was a complete implementation of the well-known open source algorithm (in most cases, Astropy.Time works, but this helps make the code standalone 
- TODO: will probably be replaced with python DateTime routines for true independence


