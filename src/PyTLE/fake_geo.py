from datetime import datetime
# -------------------------------------------------------------
def generate_fake_GEO( ):
    YEAR_DAY = datetime.utcnow().strftime('%y%j')
    GEOL1 = '1 00001U SYNTGEO  {}.00000000  .00000000  00000+0  00000+0 0  9999'.format( YEAR_DAY )
    GEOL2 = '2 00001   0.0000 000.0000 0002168 000.0000 000.0000  1.00270000  9999'
    GEO_epoch = datetime.strptime( YEAR_DAY,'%y%j' )
    return GEOL1, GEOL2
