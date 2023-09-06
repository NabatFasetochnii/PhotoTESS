from photutils.centroids import centroid_sources, centroid_com
from photutils.utils import CutoutImage
from scipy import stats
import glob
import gzip
import os



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
                arrowprops={'width': 0.1, 'headwidth': 0, 'linestyle': '-', 'color': 'y'},
                fontsize=4, zorder=5)


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
