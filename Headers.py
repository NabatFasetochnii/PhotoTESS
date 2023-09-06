from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time, TimeDelta
from astropy.utils import iers
import warnings
warnings.simplefilter("ignore")
iers.conf.auto_download = False


def Get_Times(Header):
    Exp = TimeDelta(Header['EXPTIME'], format='sec')
    t = Time(Header['DATE-OBS'], format='fits')
    t = t + Exp / 2.

    Obj = SkyCoord([Header['ALPHA'] + ' ' + Header['DELTA']],
                   unit=(u.hourangle, u.deg), frame='icrs')

    Site = EarthLocation.from_geodetic(lon=Header['LONGITUD'],
                                       lat=Header['LATITUDE'],
                                       height=Header['ALTITUDE'])

    helio = t.light_travel_time(Obj, 'heliocentric', location=Site)
    hjd = t + helio

    bary = t.light_travel_time(Obj, location=Site)
    bjd = t + bary

    return t, hjd.jd[0], bjd.jd[0]
