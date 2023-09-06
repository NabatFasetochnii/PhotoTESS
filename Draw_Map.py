from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np
from photutils.aperture import CircularAperture


def Draw_Map(Image, Size, Header, Cat, Name, RAper, Object, DT):
    wcs = WCS(Header)
    Image = np.log10(Image)
    X = Image.shape[1]
    Y = Image.shape[0]
    _mean, _median, _std = sigma_clipped_stats(Image[Y - 50:Y + 50, X - 50:X + 50])
    _max = _median + 10. * _std
    _min = _median - 1. * _std

    plt.switch_backend('pdf')  # для фикса какой-то тупой ошибки в нарнии pyplot
    fig = plt.figure(figsize=(7, 7))

    ax = plt.subplot(projection=wcs, position=[0.1, 0.1, 0.8, 0.8])
    plt.imshow(Image, vmin=_min, vmax=_max, cmap='gray_r')

    ax.set_xlim(int((Header['Naxis2'] - Size) / 2),
                int((Header['Naxis2'] + Size) / 2))
    ax.set_ylim(int((Header['Naxis1'] - Size) / 2),
                int((Header['Naxis1'] + Size) / 2))

    XY = wcs.all_world2pix(Cat['Ra'], Cat['Dec'], 0)
    XY = np.vstack((XY[0], XY[1])).T
    aper = CircularAperture(XY, r=RAper)
    aper.plot(color='blue', lw=1.5, alpha=0.5)
    aper = CircularAperture(XY[0], r=RAper)
    aper.plot(color='red', lw=1.5, alpha=0.8)
    for i, txt in enumerate(Cat['ID']):
        plt.annotate(txt, (XY[i, 0], XY[i, 1]), color='blue', alpha=0.8)

    if Header['CD2_2'] < 0:
        plt.gca().invert_xaxis()
    if Header['CD1_1'] > 0:
        plt.gca().invert_yaxis()

    Title = Object + ', ' + DT + '\n'
    Title += 'Filter=' + Header['FILTER']
    Title += ', aperture radius =' + '{:.1f}'.format(3600 * RAper * abs(Header['CD1_1'])) + '"'
    plt.title(Title)
    ax.coords[1].set_ticklabel(rotation=90)
    ax.coords[0].set_major_formatter('hh:mm:ss')
    ax.coords[1].set_major_formatter('dd:mm:ss')
    ax.coords[0].set_axislabel('RA')
    ax.coords[1].set_axislabel('Dec')
    ax.coords.grid(color='blue', ls='--', alpha=0.7)
    # plt.show()
    fig.savefig(Name)
