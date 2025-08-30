from photutils.centroids import centroid_sources, centroid_com
from photutils.utils import CutoutImage
from scipy import stats
import glob
import gzip
import os
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time, TimeDelta
from astropy.utils import iers
import warnings
warnings.simplefilter("ignore")
iers.conf.auto_download = False


def get_fits_list(path2data):
    file_list = []
    for f in os.listdir(path2data):
        if f.count('.fits') or f.count('.fts') or f.count('.fit'):
            file_list.append(path2data + '/' + f)
    return file_list


def get_com(Data, Cat, Bbox):
    x, y = centroid_sources(Data, Cat['X'], Cat['Y'], box_size=Bbox,
                            centroid_func=centroid_com)
    Cat['X'] = x
    Cat['Y'] = y
    new_Max = []
    for Obj in Cat:
        new_Max.append(CutoutImage(Data, (Obj['Y'], Obj['X']),
                                   (5, 5), mode='partial').data.max())
    Cat.add_column(new_Max, name='Max')
    return Cat


def Get_Times(Header, is_master):
    Exp = TimeDelta(Header['EXPTIME'], format='sec')
    t = Time(Header['DATE-OBS'], format='fits')
    t = t + Exp / 2.

    Obj = SkyCoord([Header['ALPHA'] + ' ' + Header['DELTA']],
                   unit=(u.hourangle, u.deg), frame='icrs')

    Site = EarthLocation.from_geodetic(lon=Header['LONGITUD'] if is_master else Header['LONGDEG'],
                                       lat=Header['LATITUDE'] if is_master else Header['LATDEG'],
                                       height=Header['ALTITUDE'])

    helio = t.light_travel_time(Obj, 'heliocentric', location=Site)
    hjd = t + helio

    bary = t.light_travel_time(Obj, location=Site)
    bjd = t + bary

    return t, hjd.jd[0], bjd.jd[0]


def unzip(path2data, is_del_zips=True):
    print(f'Unziping fits-files in {path2data}:')
    a = glob.glob(path2data + '/' + '*.fits.gz')
    for i in range(0, len(a)):
        print('unpacking', a[i])
        tar = gzip.open(a[i], 'rb')
        outF = open(a[i][0:-3], 'wb')
        outF.write(tar.read())
        tar.close()
        outF.close()
        if is_del_zips:
            os.remove(a[i])
    print('Unziping done')


def print_info(TTF):
    print('Target: ', TTF['Name'][0])
    print('V mag: ', TTF['V'][0])
    print('Start-end: ', TTF['start time'][0], '-', TTF['end time'][0])
    print('Depth: ', TTF['depth(mmag)'][0])
    print('Priority:', TTF['priority'][0])
    print('Comments: ', TTF['comments'][0])
    return TTF['Name'][0]


def target_to_tic_number(target: str):
    if 'TIC' in target:
        target = target[3:]
    if ' ' in target:
        return target.split(' ')
    if '.' in target:
        return target.split('.')
    if '_' in target:
        return target.split('_')
    return target


def draw_my_annotate(ax, t):
    ax.annotate('{:.3f}'.format(t), xy=(t, 1.006), xytext=(t, 0.973), va='center', ha='center', color='y',
                arrowprops={'width': 0.9, 'headwidth': 0, 'linestyle': '-', 'color': 'y'},
                fontsize=12, zorder=5)


def what_transit(time_series, t0, t1):
    deltaTime = [time_series[i] - time_series[i - 1] for i in range(1, len(time_series + 1))]
    gap = max(deltaTime) > stats.mode(deltaTime)[0]
    ingress = time_series[0] < t0
    egress = time_series[-1] > t1
    if ingress and egress and not gap:
        return 'a full transit'
    elif ingress and egress and gap:
        return 'a full with gapped'
    elif ingress and not egress and not gap:
        return 'an ingress '
    elif ingress and not egress and gap:
        return 'a gapped ingress'
    elif not ingress and egress and not gap:
        return 'an egress'
    elif not ingress and egress and gap:
        return 'a gapped egress'
