# ###############################################################################
# MIT License
# 
# Copyright (c) 2023 Kerry Wood
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ###############################################################################

from datetime import datetime, timedelta
import numpy as np
import PyTLE
from PyTLE.mapper import mapper

# -------------------------------------------------------------
#def generate_fake_GEO( ):
#    YEAR_DAY = datetime.utcnow().strftime('%y%j')
#    GEOL1 = '1 00001U SYNTGEO  {}.00000000  .00000000  00000+0  00000+0 0  9999'.format( YEAR_DAY )
#    GEOL2 = '2 00001   0.0000 000.0000 0002168 000.0000 000.0000  1.00270000  9999'
#    GEO_epoch = datetime.strptime( YEAR_DAY,'%y%j' )
#    return GEOL1, GEOL2

from PyTLE.utils import julian

# common fields
# TODO
# NOTE: in from_array below, we assume that the new range is 0-1 (modulo operator), this should be fixed if we want this to be variable

MAP = [
        #( human name, data struct name, range, map-to-range)
        ('mean_motion',     '_mm'       , mapper( [0,20]    ,[0,1] ) ),
        ('eccentricity',    '_ecc'      , mapper( [1e-15,1] ,[0,1] ) ),
        ('inclination',     '_incl'     , mapper( [0,180]   ,[0,1] ) ),
        ('argp',            '_argp'     , mapper( [0,360]   ,[0,1], range1_period=360, range2_period=1 ) ),
        ('raan',            '_raan'     , mapper( [0,360]   ,[0,1], range1_period=360, range2_period=1 ) ),
        ('mean_anomaly',    '_ma'       , mapper( [0,360]   ,[0,1], range1_period=360, range2_period=1 ) ),
        ('n_dot',           '_ndot'     , mapper( [-1,1]    ,[0,1] ) ),
        ('n_dot_dot',       '_ndotdot'  , mapper( [-1,1]    ,[0,1] ) )
        ]

#MAP = [
#        #( human name, data struct name, range, map-to-range)
#        ('mean_motion',     '_mm'       ,[0,20]     ,[0,1] ),
#        ('eccentricity',    '_ecc'      ,[1e-15,1]  ,[0,1] ),
#        ('inclination',     '_incl'     ,[0,180]    ,[0,1] ),
#        ('argp',            '_argp'     ,[0,360]    ,[0,1] ),
#        ('raan',            '_raan'     ,[0,360]    ,[0,1] ),
#        ('mean_anomaly',    '_ma'       ,[0,360]    ,[0,1] ),
#        ('n_dot',           '_ndot'     ,[-1,1]     ,[0,1] ),
#        ('n_dot_dot',       '_ndotdot'  ,[-1,1]     ,[0,1] )
#        ]

# for type 0/2
MAP_T0  = MAP + [ ('Bstar','_bstar',mapper( [-1,1],[0,1] )) ]
NAME_T0 = { X[0] : i for i,X in enumerate(MAP_T0) }
POS_T0  = { i : X[0] for i,X in enumerate(MAP_T0) }
N_T0    = len(MAP_T0)

# for type 4
MAP_T4 = MAP + [ 
        ('Bterm',           '_B',       mapper([-1,1]      ,[0,1] ) ),
        ('AGOM',            '_agom',    mapper([1e-15,100] ,[0,1] ) )
        ]
NAME_T4 = { X[0] : i for i,X in enumerate(MAP_T4) }
POS_T4  = { i : X[0] for i,X in enumerate(MAP_T4) }
N_T4    = len(MAP_T4)



# -----------------------------------------------------------------------------------------------------
class tle_fitter( PyTLE.TLE ):
    '''
    convenience routines to map internal TLE fields to a range an optimizer can use (generally 0--1)
    '''
    def __init__(self, inTLE : PyTLE.TLE = None , tletype : int = 0):
        if inTLE is None:
            if tletype == 0 : self._tle = PyTLE.TLE.get_type0()
            if tletype == 4 : self._tle = PyTLE.TLE.get_type4()
        else: self._tle = inTLE

    def get_type4( self ):
        self._tle = PyTLE.TLE.get_type4()
        return self

    @property 
    def epoch( self ): return self._tle.epoch
    
    @epoch.setter
    def epoch( self, epoch ): self._tle.epoch = epoch

    def set_epoch( self, epoch ):
        self._tle.epoch = epoch
        return self

    @property
    def satno( self ): return self._tle.satno

    @satno.setter
    def satno( self, newno ) : 
        # TODO : support for alpha
        assert newno > 0 and newno < 99999
        self._tle.satno = newno

    def set_satno( self, satno ):
        self._tle.satno = satno
        return self

    def set_note( self, note ):
        self._tle.set_note( note )
        return self

    def get_map( self ):
        if self._tle._type == 0 or self._tle._type == 2: return MAP_T0
        if self._tle._type == 4 : return MAP_T4

    def _val_to_mapval( self, M ):
        human, field, fmap = M
        # use the mapper to get the value
        return fmap.forward( getattr( self._tle,field) )
        #return np.interp( getattr(self._tle,field), orig_range, new_range )

    def to_array( self ):
        return np.array( [ self._val_to_mapval(X) for X in self.get_map() ] )

    def name_to_pos( self ):
        if self._tle._type == 0 or self._tle._type == 2: return NAME_T0
        if self._tle._type == 4 : return NAME_T4 
    
    def pos_to_name( self ):
        if self._tle._type == 0 or self._tle._type == 2: return POS_T0
        if self._tle._type == 4 : return POS_T4

    def get_fields( self, fields ):
        A      = self.to_array()
        N2P    = self.name_to_pos()
        return np.array( [ A[ N2P[F] ] for F in fields ] )

    def from_fields( self, fields, arr ):
        assert len(arr) == len(fields)
        A      = self.to_array()
        lookup = self.name_to_pos()
        for fname, val in zip( fields, arr ):
            A[ lookup[fname] ] = val
        self.from_array( A )
        return self
            
    def from_dict( self,  arr_dict ):
        ''' arr_dict = {'mean_motion' : 0.3, 'inclination' : 0.1...} '''
        A      = self.to_array()
        lookup = self.name_to_pos()
        for k,v in arr_dict.items():
            A[ lookup[k] ] = v
        self.from_array( A )
        return self

    def to_dict( self ):
        A = self.to_array()
        return { X[0] : A[i] for i,X in enumerate( self.get_map() ) }

    def from_array( self, array, note=None, satno=None, epoch=None ):
        # assume that order is preserved
        MAP = self.get_map()
        for i,M in enumerate(MAP): 
            human, field, fmap= M
            # map it back to the range
            #val = (array[i] + 1) % 1   # <---- TODO: this should auto-range, assumed 1 for now 
            #setattr(self._tle,field,np.interp( val, new_range, orig_range ) )
            setattr( self._tle, field, fmap.backward( array[i] ) )
        if satno : self._tle.satno = satno
        if note  : self._tle.set_note( note )
        if epoch : self._tle.epoch = epoch
        return self

    def testme( self, **kwargs ):
        print(kwargs)

    @staticmethod
    def parseLines( L1, L2 ):
        return tle_fitter( PyTLE.TLE.parseLines(L1,L2) )

    def __str__( self ):
        return "\n".join( self._tle.generateLines() )
    
    def generateLine1( self ):
        return self._tle.generateLine1() 

    def generateLine2( self ):
        return self._tle.generateLine2() 
    
    def generateLines( self ):
        return self._tle.generateLines() 


def test() :
    from sgp4.earth_gravity import wgs72
    from sgp4.io import twoline2rv
    from sgp4.propagation import sgp4 as sgprop

    # Aerocube 12A
    L1 = '1 43556U 18046C   22321.55519027  .00025005  00000+0  49749-3 0  9993'
    L2 = '2 43556  51.6329 154.1269 0008144 222.8163 137.2191 15.46745497242947'

    # parse the TLE and wrap a fitter
    TEST = tle_fitter( PyTLE.TLE.parseLines( L1, L2 ) )
    print('Str of fitter')
    print(str(TEST))

    # get array and then parse back into structure
    arr = TEST.to_array()
    print('Array: {}'.format(arr))
    print('Re-injecting...')
    TEST.from_array( arr )
    print(str(TEST))

    print('Modifying entries and re-injecting (wrap test)')
    N = len(arr)
    mod =arr + np.ones(N) * 0.7
    print( TEST.from_array( mod ) ) 



    # test ephemeris
    tle = twoline2rv( L1, L2, wgs72 )
    print()
    print('tle epoch is ', julian.from_jd( tle.jdsatepoch ) )
    tledate = julian.from_jd( tle.jdsatepoch )
    
    mins   = np.arange(0,1440*5,10)
    jdates = tle.jdsatepoch + mins/1440
    eph = np.vstack( [np.hstack( sgprop(tle,D)) for D in mins ] )

        # init as if we only had a state-vector
    FIT = tle_fitter(PyTLE.TLE.fromPV( epoch=tledate, P=eph[0,0:3], V=eph[0,3:] ) )
    print()
    print( str(FIT) )


# =====================================================================================================
if __name__ == '__main__':
    test()
