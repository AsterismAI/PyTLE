import numpy as np


class mapper:
    '''
    a mapper class for orbital elements; we want to wrap some values and not wrap others
        eg. wrap anomaly
            do not wrap inclination
    '''
    def __init__( self, range1 : list, range2 : list, range1_period = None, range2_period = None ):
        '''
        if range_period is not input, it defaults to clamping to range
        '''
        assert len(range1) == 2
        assert len(range2) == 2
        # set the ranges
        self._range1 = range1
        self._range2 = range2

        if range1_period:
            self.forward = lambda X: np.interp( X % range1_period, self._range1, self._range2 )
        else:
            self.forward = lambda X: np.interp( X, self._range1, self._range2 )
        
        if range2_period:
            self.backward = lambda X: np.interp( X % range2_period, self._range2, self._range1 )
        else:
            self.backward = lambda X: np.interp( X, self._range2, self._range1 )

    def forward( self, value ):
        return self.forward( value )

    def backward( self, value ):
        return self.backward( value )

# =====================================================================================================
if __name__ == "__main__":
    incl = mapper( [0,180], [0,1] )
    def print_incl( val ):
        print('Inclination test forward  for {}   : {}'.format( val, incl.forward(val) ) )
        print('Inclination test backward for {}   : {}'.format( val, incl.backward( incl.forward(val) ) ) )
    for v in [0,90,180,181,270] : print_incl(v)

    anom = mapper( [0,360], [0,1], range1_period=360, range2_period=1)
    def print_anom( val ):
        print('Anom test forward  for {}   : {}'.format( val, anom.forward(val) ) )
        print('Anom test backward for {}   : {}'.format( val, anom.backward( anom.forward(val)  ) ))
    
    print()
    for v in [0,90,180,181,270,360,390] : print_anom(v)
