from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np
from photutils.aperture import CircularAperture


def draw_map(image, size, header, cat, name, r_aper, tess_object, dt, scale):
    wcs = WCS(header)
    image = np.log10(image)
    X = image.shape[1]
    Y = image.shape[0]
    _mean, _median, _std = sigma_clipped_stats(image[Y - 50:Y + 50, X - 50:X + 50])
    _max = _median + 10. * _std
    _min = _median - 1. * _std

    plt.switch_backend('pdf')  # для фикса какой-то тупой ошибки в нарнии pyplot
    fig = plt.figure(figsize=(7, 7))

    ax = plt.subplot(projection=wcs, position=[0.1, 0.1, 0.8, 0.8])
    plt.imshow(image, vmin=_min, vmax=_max, cmap='gray_r')

    ax.set_xlim(int((header['Naxis2'] - size) / 2),
                int((header['Naxis2'] + size) / 2))
    ax.set_ylim(int((header['Naxis1'] - size) / 2),
                int((header['Naxis1'] + size) / 2))

    XY = wcs.all_world2pix(cat['Ra'], cat['Dec'], 0)
    XY = np.vstack((XY[0], XY[1])).T
    aper = CircularAperture(XY, r=r_aper)
    aper.plot(color='blue', lw=1.5, alpha=0.5)
    aper = CircularAperture(XY[0], r=r_aper)
    aper.plot(color='red', lw=1.5, alpha=0.8)
    for i, txt in enumerate(cat['ID']):
        plt.annotate(txt, (XY[i, 0], XY[i, 1]), color='blue', alpha=0.8)

    # scale = 3600 * np.sqrt(header['CD1_1'] ** 2 + header['CD1_2'] ** 2)
    r_aper_2_5 = 2.5*60 / scale

    aper = CircularAperture(XY[0], r=r_aper_2_5)  # only for target
    aper.plot(color='y', lw=1.2, alpha=0.7)

    plt.text(XY[0, 0], XY[0, 1] - 84, s='Radius = 2.5\'', color='y', alpha=0.7)

    if header['CD2_2'] < 0:
        plt.gca().invert_xaxis()
    if header['CD1_1'] > 0:
        plt.gca().invert_yaxis()

    title = tess_object + ', ' + dt + '\n'
    title += 'Filter=' + header['FILTER']
    title += ', aperture radius =' + '{:.1f}'.format(r_aper * scale) + '"'
    plt.title(title)
    ax.coords[1].set_ticklabel(rotation=90)
    ax.coords[0].set_major_formatter('hh:mm:ss')
    ax.coords[1].set_major_formatter('dd:mm:ss')
    ax.coords[0].set_axislabel('RA')
    ax.coords[1].set_axislabel('Dec')
    ax.coords.grid(color='blue', ls='--', alpha=0.7)
    # plt.show()
    fig.savefig(name)
