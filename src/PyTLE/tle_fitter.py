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
from .mapper import mapper

# -------------------------------------------------------------
#def generate_fake_GEO( ):
#    YEAR_DAY = datetime.utcnow().strftime('%y%j')
#    GEOL1 = '1 00001U SYNTGEO  {}.00000000  .00000000  00000+0  00000+0 0  9999'.format( YEAR_DAY )
#    GEOL2 = '2 00001   0.0000 000.0000 0002168 000.0000 000.0000  1.00270000  9999'
#    GEO_epoch = datetime.strptime( YEAR_DAY,'%y%j' )
#    return GEOL1, GEOL2

from .utils import julian

# this version uses the mapper to move values into and out of the proper range
# mostly this is done to AVOID wrapping inclination (not much use to a solver)
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
    def __init__(self, in_arg  = None ):
        super().__init__()
              
    def get_available_field_names( self ):
        return list( self.name_to_pos().keys() )
    
    def _val_to_mapval( self, M ):
        human, field, fmap = M
        # use the mapper to get the value
        return fmap.forward( getattr( self,field) )
        
    def get_map( self ):
        if self._type == 0 or self_type == 2: return MAP_T0
        if self._type == 4 : return MAP_T4

    def name_to_pos( self ):
        if self._type == 0 or self._type == 2: return NAME_T0
        if self._type == 4 : return NAME_T4 
    
    def pos_to_name( self ):
        if self._type == 0 or self._type == 2: return POS_T0
        if self._type == 4 : return POS_T4
        
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
    
    def to_array( self ):
        return np.array( [ self._val_to_mapval(X) for X in self.get_map() ] )
    
    def from_array( self, array ):
        # assume that order is preserved
        for i,M in enumerate( self.get_map() ): 
            human, field, fmap= M
            setattr( self, field, fmap.backward( array[i] ) )
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
# ------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    L1 = '1 43556U 18046C   22321.55519027  .00025005  00000+0  49749-3 0  9993'
    L2 = '2 43556  51.6329 154.1269 0008144 222.8163 137.2191 15.46745497242947'
    X = tle_fitter().parseLines( L1, L2)
    print(X.get_map())
    print(X.get_available_field_names())
    print(X.to_array())
    print( X.get_fields( X.get_available_field_names() ))
    print( X.get_fields( X.get_available_field_names()[:3] ))
    print(X.to_dict())
    X.fromPV( datetime(year=2020,month=1,day=1), np.array([7000,0,0]), np.array([0,9,2] ))
    print(X)