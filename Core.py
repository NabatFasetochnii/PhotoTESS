import os
# disable warnings
import warnings

# import astropy modules
import astropy.io.fits as fits
import astropy.wcs as wcs
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.stats import SigmaClip
from photutils.aperture import CircularAperture, aperture_photometry
from photutils.background import Background2D, MedianBackground
from unidecode import unidecode

# import my modules
from Draw_Map import draw_map
from Get_Cat import Get_GAIA
from Postproc import Post_Proc
from TTF import TTF_Query
from Utils import get_com, unzip, print_info, get_fits_list, Get_Times

warnings.simplefilter("ignore")


def photometry(path2data, RAper=None, observer='chazov'):

    tess_object = None
    priority = 5  # set 5 for old data
    is_master = None

    ######################################################################
    Cat_R = 0.1  # catalog each cone radius
    Cat_len = 400  # catalog max length
    V_lim = 15  # catalog max depth
    Bbox = 7  # centering/morphology box (pixels)
    Pic_size = 256  # picture size (pixels)

    ######################################################################
    # unzip and del all .gz files
    unzip(path2data)

    ######################################################################
    # read directory and create list of fits-files
    file_list = get_fits_list(path2data)

    # set paths and suffix for output files
    Path2Save = path2data + '/Photometry'
    if not os.path.exists(Path2Save):
        os.makedirs(Path2Save)

    # clean old photometry
    for f in os.listdir(Path2Save):
        if f.count('Phot'):
            os.remove(Path2Save + '/' + f)

    # read first frame for object name and other information
    print('First frame: ' + file_list[0].split('/')[-1])
    hduList = fits.open(file_list[0])

    Header = hduList[0].header
    try:
        scope = Header['PROJECT']  # MASTER-II-Ural
        is_master = True
    except:
        scope = Header['TELESCOP']  # APM-RoboPhot
        is_master = False
    scale = 3600 * np.sqrt(Header['CD1_1'] ** 2 + Header['CD1_2'] ** 2)
    Data = hduList[0].data.copy()
    hduList.verify('fix')
    hduList.close()

    if RAper is None:
        if is_master:
            RAper = [3, 4, 5]  # apertures list
        else:
            RAper = [5, 6, 7]  # apertures list

    # get times: jd, hjd, bjd
    Times = Get_Times(Header, is_master)
    DT = Times[0].datetime.strftime('%m-%d-%Y')

    # try:
    # get info from TESS database
    if os.path.isfile(Path2Save + '/TTF_Info.txt'):
        print('Local TTF')
        TTF = ascii.read(Path2Save + '/TTF_Info.txt',
                         delimiter='\t', format='commented_header')
    else:
        print('Download TTF')
        if tess_object is None:
            tess_object = Header['OBJNAME']
            tess_object = tess_object[0:3] + '+' + tess_object[3:].replace('_', '.')
        TTF = TTF_Query(DT, priority, 0.5, 16, tess_object)
        # replace unicode symbols to ascii
        TTF['comments'][0] = unidecode(TTF['comments'][0])
        ascii.write(TTF, Path2Save + '/TTF_Info.txt',
                    fast_writer=False,
                    overwrite=True, delimiter='\t', format='commented_header', fill_values=[(ascii.masked, '0')])

    tess_object = print_info(TTF)

    # make coordinate object
    C = SkyCoord(Header['ALPHA'] + ' ' + Header['DELTA'],
                 unit=(u.hourangle, u.deg), frame='icrs')

    # read local catalog or create from GAIA data
    if os.path.isfile(Path2Save + '/Cat.txt'):
        print('Local catalog')
        Catalog = ascii.read(Path2Save + '/Cat.txt')
    else:
        print('Download GAIA')
        Catalog = Get_GAIA(C.ra.degree, C.dec.degree, Cat_R, V_lim, Cat_len)
        ascii.write(Catalog, Path2Save + '/Cat.txt', overwrite=True, delimiter='\t')

    # # set prefix and suffix for output files

    # draw and save picture
    pic_name = (Path2Save + '/' + tess_object.replace(' ', '').replace('.', '-') + '_' +
                Times[0].datetime.strftime('%Y%m%d') + ('_kourovka0.4_' if is_master else '_kourovka0.6_') +
                Header['FILTER'] + '_field.pdf')
    draw_map(Data, Pic_size, Header, Catalog, pic_name, RAper[0], tess_object, DT, scale)

    # open log file
    df = open(Path2Save + '/Time.txt', 'w')
    if is_master:
        df.write('DATE-OBS\tJD\tHJD\tBJD\t' +
                 'EXPTIME\tFILTER\tAIRMASS\tEXTINCT\t' +
                 'SKY-TEMP\tSEXFWHM\tSEXELL \t' +
                 'Sky\tX\tY\tMax\n')
    else:
        df.write('DATE-OBS\tJD\tHJD\tBJD\t' +
                 'EXPTIME\tFILTER\tEXTINCT\t' +
                 'SKY-TEMP\tSEXFWHM\tSEXELL\t' +
                 'Sky\tX\tY\tMax\n')

    Outputs = [open(Path2Save + '/Phot' + str(r) + '.txt', 'a') for r in RAper]

    counter = len(file_list)
    for f in file_list:
        print('----------<>-----------')
        counter = counter - 1
        print('Frames: ' + str(counter))
        print(f.split('/')[-1])

        hduList = fits.open(f)
        Header = hduList[0].header
        Data = hduList[0].data
        hduList.verify('fix')
        hduList.close()
        w = wcs.WCS(Header)

        Times = Get_Times(Header, is_master)
        df.write(Times[0].datetime.isoformat(timespec='milliseconds') + '\t')  # DATE-OBS
        df.write('{:.7f}'.format(Times[0].jd) + '\t')  # JD
        df.write('{:.7f}'.format(Times[1]) + '\t')  # HJD
        df.write('{:.7f}'.format(Times[2]) + '\t')  # BJD
        df.write('{:.1f}'.format(Header['EXPTIME']) + '\t')  # EXPTIME
        df.write(Header['FILTER'] + '\t')  # FILTER

        if is_master:
            df.write('{:.2f}'.format(Header['AIRMASS']) + '\t')
        df.write('{:.2f}'.format(Header['EXTINCT']) + '\t')
        df.write('{:.2f}'.format(Header['SKY-TEMP']) + '\t')
        df.write('{:.2f}'.format(Header['SEXFWHM']) + '\t')
        df.write('{:.2f}'.format(Header['SEXELL']) + '\t')
        # df.write('{:.2f}'.format(Header['ZEROPOI']) + '\t')

        X, Y = w.all_world2pix(Catalog['Ra'], Catalog['Dec'], 0)
        Catalog.add_columns([X, Y], names=['X', 'Y'])

        # delete background
        sigmaclip = SigmaClip(sigma=3.)
        bkg_estimator = MedianBackground()
        bkg = Background2D(Data, (100, 100), filter_size=(9, 9),
                           sigma_clip=sigmaclip, bkg_estimator=bkg_estimator)
        Data = Data - bkg.background
        Sky = np.median(bkg.background)
        df.write('{:.1f}'.format(Sky) + '\t')  # Sky
        print('Sky={0:.1f}'.format(Sky))

        # centroiding and image properties
        Catalog = get_com(Data, Catalog, Bbox)  # center of mass / fast and better!
        #     Catalog = get_1dg(Data, Catalog, Bbox) # Gauss 1d fit

        if is_master:
            df.write('{0:.3f}\t{1:.3f}\t{2:.1f}\n'.format(Header['XSTART'] + Catalog['X'][0],
                                                          Header['XSTART'] + Catalog['Y'][0],
                                                          Catalog['Max'][0]))
        else:
            df.write('{0:.3f}\t{1:.3f}\t{2:.1f}\n'.format(Catalog['X'][0],
                                                          Catalog['Y'][0],
                                                          Catalog['Max'][0]))
        print('Max count={0:.1f}'.format(Catalog['Max'][0]))

        Positions = np.transpose([Catalog['X'], Catalog['Y']])
        Stellar_aper = [CircularAperture(Positions, r=r) for r in RAper]
        Stellar_phot = aperture_photometry(Data, Stellar_aper, method='exact')

        # write to log
        for count, value in enumerate(RAper):
            Flux = Stellar_phot['aperture_sum_' + str(count)].value
            np.savetxt(Outputs[count], Flux, fmt='%1.3f', delimiter='\t', newline='\t')
            Outputs[count].write('\n')

        Catalog.remove_columns(names=('X', 'Y', 'Max'))

    for out in Outputs:
        out.close()

    df.close()

    Post_Proc.post_proc(apers=RAper, path2data=Path2Save, is_master=is_master, observer=observer, scale=scale)
