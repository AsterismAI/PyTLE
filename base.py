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
from sgp4.ext import rv2coe
import numpy as np

from .alpha import alpha_to_integer, integer_to_alpha

from .formatters import generate_expo_format, process_expo_format
from .formatters import epoch_str_todatetime, datetime_to_epochstr

WGS84  = 398600.5
gL1    = '1 25544U 98067A   23137.83559306  .00011914  00000-0  21418-3 0  9990'
gL2    = '2 25544  51.6409 118.9691 0006630 359.0829  72.4864 15.50282135397083'
gEPOCH = datetime.fromisoformat('2000-01-01T12:00:00.000') 

# -----------------------------------------------------------------------------------------------------
def format_ecc( ecc ):
    return '{:9.7f}'.format(ecc)[2:].ljust(7,'0')

# -----------------------------------------------------------------------------------------------------
class TLE:
    def __init__(self, L1=None, L2=None):
        self.clear()

    def clear( self ):
        self._epoch = gEPOCH
        self._line1 = None
        self._line2 = None
        self._satno = 99999
        self._class = 'U'
        self._intld = ''
        self._elset = 9999
        self._incl    = 0.
        self._raan    = 0.
        self._ecc     = 0.
        self._argp    = 0.
        self._ma      = 0.
        self._mm      = 0. 
        # type 2
        self._ndot    = 0.
        self._ndotdot = 0.
        self._type    = 0.
        self._bstar   = 0.
        # type 4
        self._B       = 0.
        self._agom    = 0.

        # human readable / helper
        self._perigee = None
        self._apogee  = None

        # old fields
        self._elset   = 0

    @property
    def epoch( self ): return self._epoch

    @epoch.setter
    def epoch( self, epoch ):
        if isinstance(epoch,str)      : self._epoch = datetime.fromisoformat( epoch )
        if isinstance(epoch,datetime) : self._epoch = epoch

    @property
    def satno( self ): return self._satno

    @satno.setter
    def satno( self, newno ) : 
        # TODO : support for alpha
        assert newno > 0 and newno < 99999
        self._satno = newno

    @property
    def inclination( self ): return self._incl

    @inclination.setter
    def inclination( self, incl ):
        self._incl = (incl + 180) % 180

    @property
    def RAAN( self ) : return self._raan

    @RAAN.setter
    def RAAN( self, raan ):
        self._raan = (raan + 360) % 360

    @property 
    def eccentricity( self ): return self._ecc

    @eccentricity.setter
    def eccentricity( self, ecc ):
        self._ecc = np.clip( ecc, 0, 1 )

    @property 
    def arg_perigee( self ): return self._argp
    @arg_perigee.setter
    def arg_perigee( self, argp ):
        self._argp = (argp+ 360) % 360

    @property 
    def mean_anomaly( self ) : return self._ma
    @mean_anomaly.setter
    def mean_anomaly( self, ma ):
        self._ma = (ma + 360) % 360

    @property 
    def mean_motion( self ) : return self._mm
    # TODO: add some error checking here
    @mean_motion.setter
    def mean_motion( self, mm ): self._mm = mm

    @property
    def B( self ): return self._B
    @B.setter
    def B(self, newB ) : 
        self._B = newB

    @property
    def AGOM( self ) : return self._agom
    @AGOM.setter
    def AGOM( self, agom ):
        self._agom = agom

    def set_note( self, note : str ):
        self._intld = note[:8]
        return self
        
    @property
    def perigee( self ):
        if self._perigee is None: self._calculate_apogee_perigee()
        return self._perigee
    
    @property
    def apogee( self ):
        if self._apogee is None: self._calculate_apogee_perigee()
        return self._apogee

    def parseDate( self, S ):
        ''' assume this is just the date string '''
        self._epoch = epoch_str_todatetime( S )

    @staticmethod 
    def parseLines( L1, L2 ):
        if L1[62] == '0' or L1[62] == '2' : return TLE_2( L1, L2 )
        if L1[62] == '4' : return TLE_4( L1, L2 )

    @staticmethod
    def get_type0( ):
        return TLE_2()

    @staticmethod
    def get_type4( ):
        return TLE_4()


    def _calculate_apogee_perigee( self, earth_rad = 6378.135):
        ''' Default value for earth_rad is taken from space-track.
        space-track : https://www.space-track.org/documentation#/faq
        Additional references: http://www.satobs.org/seesat/Dec-2002/0197.html
        [I192G'''
        semi_major = (8681663.653 / self.mean_motion) ** (2.0/3.0)
        self.perigee = ( semi_major * (1 - self.eccentricity) ) - earth_rad
        self.apogee =  ( semi_major * (1 + self.eccentricity) ) - earth_rad

    def _calculate_apogee_perigee( self, earth_rad = 6378.135):
        ''' Default value for earth_rad is taken from space-track.
        space-track : https://www.space-track.org/documentation#/faq
        Additional references: http://www.satobs.org/seesat/Dec-2002/0197.html
        '''
        semi_major = (8681663.653 / self.mean_motion) ** (2.0/3.0)
        self._perigee = ( semi_major * (1 - self.eccentricity) ) - earth_rad
        self._apogee =  ( semi_major * (1 + self.eccentricity) ) - earth_rad


    @staticmethod
    def fromCOE(epoch : datetime,
                type  : int = 0,
                satno : int = 99999,
                a=7000, ecc=1e-10, incl=1e-3, argp=0, raan=0, mean_anomaly=0,   # COE elements
                bstar : float = 0,                                              # TLE type 0/2 only
                bterm : float  = 0,
                agom  : float  = 0,
                EARTHMU : float = WGS84,
                **kwargs):
        '''
        degrees and km
        '''
        if type == 0 or type == 2: 
            tle = TLE_2()
            tle._bstar = bstar
        if type == 4 : 
            tle = TLE_4()
            tle._B = bterm
            tle._agom = agom
        if a < 0:
            print('cannot init from COE with semi-major axis (a) < 0')
            return tle 
        tle._satno = satno
        tle._incl = incl
        tle._ecc  = ecc
        tle._argp = argp
        tle._raan = raan 
        tle._ma   = mean_anomaly
        # calculate the mean motion (these are km values, so we get rads/s, convert to TLE units)
        mm = np.sqrt(EARTHMU / a ** 3)
        tle._mm = (mm * 86400) / (2 * np.pi)
        tle._epoch = epoch
        return tle


    @staticmethod
    def fromPV(epoch : datetime, 
               P : np.array,
               V : np.array,
               type : int = 0,
               satno : int = 99999,
               bstar : float = 0,                                              # TLE type 0/2 only
               bterm : float  = 0,
               agom  : float  = 0,
               EARTHMU : float = WGS84,
               **kwargs):
        '''
        fromPV : given state position and velocity (in TEME), build an initial TLE
        note   : this is *not* going to build mean elements
        '''
        # return p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper
        if V[2] == 0:
            raise Exception('cannot init an orbit with perfectly zero inclination (velocity[Z] ~ 1e-5km/s minimum)')
            return TLE.fromCOE( epoch, type=type, satno=satno )

        try: p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe(P, V, EARTHMU )
        except Exception as e:
            print('could not init TLE from P: {} V: {}'.format( P,V ) )
            return TLE.fromCOE( epoch, type=type, satno=satno )

        # rv2coe in sgp4 returns > 999999 to indicate undefined or infinite... (see code)
        if any( (X > 999999. for X in [p,a,ecc,incl,omega,argp,nu,m] ) ):
            print('rv2coe returned undefined (> 999999) value')
            print([p,a,ecc,incl,omega,argp,nu,m] )
            return TLE.fromCOE( epoch , type=type)

        # a = 7000, ecc = 1e-10, incl = 1e-3, argp = 0, raan = 0, mean_anomaly = 0,
        return TLE.fromCOE( epoch, 
                type=type, 
                satno=satno, 
                a=a, 
                ecc=ecc, 
                incl=np.degrees(incl), 
                argp=np.degrees(argp), 
                raan=np.degrees(omega), 
                mean_anomaly=np.degrees(m),
                bstar = bstar,
                bterm = bterm,
                agom  = agom,
                EARTHMU = EARTHMU)


    def __str__( self ): return '\n'.join( self.generateLines() )
    
    def __repr__( self ): return str(self)




# -----------------------------------------------------------------------------------------------------
class TLE_2( TLE ):
    def __init__(self, L1=None, L2=None):
        self.clear()
        if L1 and L2: self.parseLines(L1,L2)
        self._type = 0

    def parseLine1( self, S ):
        if S[0] != '1' : raise Exception('LINE1 must begin with 1')
        self._satno = alpha_to_integer( S[2:7] )
        self._class = S[7]
        self._intld = S[9:17]
        self.parseDate( S[18:32] ) 
        self._ndot    = float( S[33:43] )
        self._ndotdot = process_expo_format( S[44:52] )
        self._bstar  = process_expo_format( S[53:61] )
        self._type   = int(S[62])
        self._elset = int(S[64:68])
    
    def parseLine2( self, S ):
        if S[0] != '2' : raise Exception('LINE2 must begin with 2')
        if alpha_to_integer( S[2:7] ) != self._satno : raise Exception('satno does not match')
        self._incl = float(S[8:16])
        self._raan = float(S[17:25])
        self._ecc  = float( '0.{}'.format( S[26:33] ) )
        self._argp = float(S[34:42])
        self._ma   = float(S[43:51])
        self._mm   = float(S[52:63])

    def parseLines( self, L1, L2 ):
        self.parseLine1( L1 )
        self.parseLine2( L2 )

    def generateLine1( self ):
        #1 25544U 98067A   23137.83559306  .00011914  00000-0  21418-3 0  9990
        L1 = '1 {:5}{:1} {:8} {:14} {} {} {} {} {}0'.format(
                integer_to_alpha( self._satno ).rjust(5,'0'),
                self._class[0],
                self._intld[:8].rjust(8,' '),
                datetime_to_epochstr( self._epoch ),
                "{:09.8f}".format( self._ndot ).ljust(9,' '),
                generate_expo_format( self._ndotdot),
                generate_expo_format( self._bstar ),
                0,  # <--------- TYPE FLAG
                "{}".format( self._elset ).rjust(4,' ')
                )

                #self._elset )
        return L1

    def generateLine2( self ):
        #L2='2 12345   9.7332 113.4837 7006332 206.5371  38.9576 01.00149480000003'
        L2 = '2 {:5} {:8} {:8} {:7} {:8} {:8} {} {}3 '.format(
                integer_to_alpha( self._satno ).rjust(5,'0'),
                "{:>8.4f}".format( self._incl )[:8], 
                "{:>8.4f}".format( self._raan)[:8], 
                #"{}".format( self._ecc )[2:10],
                format_ecc( self._ecc ),
                "{:>8.4f}".format( self._argp)[:8], 
                "{:>8.4f}".format( self._ma)[:8], 
                "{:>011.8f}".format( self._mm)[:12], 
                "{:04d}".format( self._elset)[:4]
                )
        return L2
    
    def generateLines( self ):
        return ( self.generateLine1(), self.generateLine2() )



# -----------------------------------------------------------------------------------------------------
class TLE_4( TLE ):
    def __init__(self, L1=None, L2=None):
        self.clear()
        if L1 and L2: self.parseLines(L1,L2)
        self._type = 4

    def parseLine1( self, S ):
        if S[0] != '1' : raise Exception('LINE1 must begin with 1')
        self._satno = alpha_to_integer( S[2:7] )
        self._class = S[7]
        self._intld = S[9:17]
        self.parseDate( S[18:32] ) 
        self._agom = process_expo_format( S[44:52] )
        self._B    = process_expo_format( S[53:61] )
        self._type = int(S[62])
        self._elset = int(S[64:68])
    
    def parseLine2( self, S ):
        if S[0] != '2' : raise Exception('LINE2 must begin with 2')
        if alpha_to_integer( S[2:7] ) != self._satno : raise Exception('satno does not match')
        self._incl = float(S[8:16])
        self._raan = float(S[17:25])
        self._ecc  = float( '0.{}'.format( S[26:33] ) )
        self._argp = float(S[34:42])
        self._ma   = float(S[43:51])
        self._mm   = float(S[52:63])

    def parseLines( self, L1, L2 ):
        self.parseLine1( L1 )
        self.parseLine2( L2 )

    def generateLine1( self ):
        #L1='1 12345U xyzzyz   23038.45547454 +.00000000 +46171+0 +33000-1 4 99992'
        L1 = '1 {:5}{:1} {:8} {:14} +.00000000 {} {} 4 {}'.format(
                integer_to_alpha( self._satno ).rjust(5,'0'),
                self._class[0],
                self._intld[:8].rjust(8,' '),
                datetime_to_epochstr( self._epoch ),
                generate_expo_format( self._agom ),
                generate_expo_format( self._B ),
                "{:d}".format( self._elset).rjust(4,'0')[-4:]
                 )
        return L1

    def generateLine2( self ):
        #L2='2 12345   9.7332 113.4837 7006332 206.5371  38.9576 01.00149480000003'
        L2 = '2 {:5} {:8} {:8} {:7} {:8} {:8} {} {}'.format(
                integer_to_alpha( self._satno ).rjust(5,'0'),
                "{:>8.4f}".format( self._incl )[:8], 
                "{:>8.4f}".format( self._raan)[:8], 
                #"{}".format( self._ecc )[2:10],
                format_ecc( self._ecc ),
                "{:>8.4f}".format( self._argp)[:8], 
                "{:>8.4f}".format( self._ma)[:8], 
                "{:>011.8f}".format( self._mm)[:12], 
                "{:d}".format( self._elset).rjust(4,'0')[-4:]
                )
        return L2
    
    def generateLines( self ):
        return ( self.generateLine1(), self.generateLine2() )

# -----------------------------------------------------------------------------------------------------
def demo():
    from sgp4.io import twoline2rv
    from sgp4.earth_gravity import wgs72
    from sgp4.propagation import sgp4 as sgprop
    L1='1 12345U xyzzyz   23038.45547454 +.00000000 +46171+0 +33000-1 4 99992'
    L2='2 12345   9.7332 113.4837 7006332 206.5371  38.9576 01.00149480000003'

    print('-----------------------------------------------------------------------------------------------------')
    print('Original TLE')
    print(L1)
    print(L2)
    print('-----------------------------------------------------------------------------------------------------')

    print()
    print('-----------------------------------------------------------------------------------------------------')
    print('Type 4 testing')
    print('CODE: test = TLE.parseLines( L1, L2 )')
    test = TLE.parseLines( L1, L2 )
    print('original   ', L1)
    print('generated  ',test.generateLine1())
    print('original   ',L2)
    print('generated  ',test.generateLine2())

    print()
    print('Setting international designator to {}'.format( 'ABC123' ))
    test._intld = 'ABC123'
    print('generated  ',test.generateLine1())
    print('generated  ',test.generateLine2())
    print('Setting international designator to {}'.format('ABCDEFGHIJKLMNO'))
    test._intld = 'ABCDEFGHIJKLMNO'
    print('generated  ',test.generateLine1())
    print('generated  ',test.generateLine2())

    print()
    RAAN = 15000000
    print('Setting RAAN to {}'.format( RAAN ))
    test.raan = RAAN
    print('generated  ',test.generateLine1())
    print('generated  ',test.generateLine2())

    print()
    AGOM = 10e6
    print('Setting AGOM to {}'.format(AGOM))
    test._agom = AGOM
    print('generated  ',test.generateLine1())
    print('generated  ',test.generateLine2())


    print('-----------------------------------------------------------------------------------------------------')
    print('Type 0 testing')
    print('CODE: test = TLE.parseLines( L1, L2 )')
    L1 = '1 25544U 98067A   23137.83559306  .00011914  00000-0  21418-3 0  9990'
    L2 = '2 25544  51.6409 118.9691 0006630 359.0829  72.4864 15.50282135397083'
    test = TLE.parseLines( L1, L2 )
    print('original   ', L1)
    print('generated  ',test.generateLine1())
    print('original   ',L2)
    print('generated  ',test.generateLine2())

    print()
    print('Setting RAAN to 999999999')
    test._raan = 999999999
    print('original   ', L1)
    print('generated  ',test.generateLine1())
    print('original   ',L2)
    print('generated  ',test.generateLine2())

    print()
    print('Setting inclination to 999999999')
    test._incl = 999999999
    print('original   ', L1)
    print('generated  ',test.generateLine1())
    print('original   ',L2)
    print('generated  ',test.generateLine2())


    print()
    print('Reparsing ORIGINAL and outputting dict')
    Y = TLE.parseLines( L1, L2 )
    print(Y.__dict__)

    X = TLE.fromCOE( datetime.utcnow() ).set_note('fromCOE')
    print('From COE:')
    print('\n'.join( X.generateLines())) 

    P = [9000,0,0]
    V = [0,8,0.001]
    print('\nFrom P: {} V: {}'.format( P, V ))
    X = TLE.fromPV( datetime.utcnow(), P, V,
                        type = 0 )
                        
    L1,L2 = X.generateLines()
    print(L1)
    print(L2)
    T = twoline2rv(L1,L2,wgs72)
    print('Repropagated to epoch: {}'.format( sgprop(T,0)) )
    
    

# =====================================================================================================
if __name__ == '__main__':
    demo()
