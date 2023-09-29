import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.io import ascii
from astropy.stats import sigma_clipped_stats
from astropy.time import Time as aTime
from astropy.timeseries import TimeSeries, aggregate_downsample
from pytransit import TransitAnalysis  # может быть проблема совместимости numpy 1.24 и numba (17.02.23)
from pytransit.orbits import orbits_py

from Postproc.Condition_Report import Condition_Report
from Postproc.aligner import aligner
from Utils import print_info, draw_my_annotate, what_transit
from Postproc.Check_best_aperture import check_best_aperture
from corner import corner

# install arviz, celerite, corner

warnings.simplefilter("ignore")


def post_proc(apers, observer, path2data, is_master, scale, save_figs=False):
    # tess_object = 'Obj'  # set object name if TTF data not available
    # T0 = 2459851.3069  # transit start jd time
    # T1 = 2459851.3572  # transit end jd time

    best_aperture = check_best_aperture(apers=apers, path2_data=path2data, save_fig=save_figs)
    print(f'best aperture = {best_aperture}')
    min_bin_size = 3  # bin factor for lightcurve
    max_extinction = 0.35

    observatory = 'kourovka0.4' if is_master else 'kourovka0.6'

    tess_object = None
    T0 = None
    T1 = None
    flux = None
    TTF = None

    #######################################################################
    # read data
    catalog = ascii.read(path2data + '/Cat.txt')
    time_list = ascii.read(path2data + '/Time.txt')

    for f in os.listdir(path2data):
        if f.count('Phot' + str(best_aperture)):
            flux = np.genfromtxt(path2data + '/' + f)

    if os.path.isfile(path2data + '/TTF_Info.txt'):
        TTF = ascii.read(path2data + '/TTF_Info.txt',
                         delimiter='\t', format='commented_header')
        tess_object = print_info(TTF)
        # transit start-end
        T0 = aTime(TTF['jd_start'][0] + 2450000, format='jd')
        T1 = aTime(TTF['jd_end'][0] + 2450000, format='jd')
    else:
        print('No TTF data')

    ########################################################################
    # set prefix
    transit_date = aTime(time_list['JD'][0], format='jd')
    # prefix = transit_date.datetime.strftime('%Y%m%d') + '_Kourovka0.4_' + time_list['FILTER'][0]

    exp = time_list['EXPTIME'][0]
    fil = time_list['FILTER'][0]

    save_string_for_additions = f'{path2data}/{tess_object}_{observatory}_{time_list["FILTER"][0]}_' \
                                f'{time_list["DATE-OBS"][0].split("T")[0]}'
    report_string = f"{path2data}/{tess_object.replace(' ', '').replace('.', '-')}_" \
                    f"{transit_date.datetime.strftime('%Y%m%d')}_{observatory}_{fil}"

    ########################################################################
    # condition report
    # if is_master:
    Condition_Report(tess_object, T0, T1, time_list, report_string, observatory)

    ########################################################################
    # plot raw data
    if save_figs:
        plt.figure(figsize=(7, 7))
        plt.plot(time_list['JD'], 20 - 2.5 * np.log10(flux))
        plt.ylabel('instrumental magnitude', fontsize=6)
        plt.title('Raw data', fontsize=8)
        plt.axvspan(T0.jd, T1.jd, facecolor='k', alpha=0.2)
        locs, labels = plt.xticks()
        t = aTime(locs, format='jd')
        x_ticks_labels = []
        for x in t:
            x_ticks_labels.append(str(x.iso).split(' ')[1].split('.')[0])
        plt.xticks(locs, x_ticks_labels, rotation='vertical', fontsize=6)
        plt.xlabel('transit_date-time_list (UTC), ' + time_list['DATE-OBS'][0].split('T')[0], fontsize=6)
        plt.xlim(locs[0], locs[-2])
        plt.tick_params(axis='both', labelsize=6, direction='in')
        plt.gca().invert_yaxis()
        plt.grid()
        plt.show()

        plt.savefig(save_string_for_additions + '_Raw data.pdf')

    ########################################################################
    # choose ensemble by colors and magnitudes
    # Ensemble = get_ensemble(catalog)
    # Flux_e = flux[:, Ensemble['ID']]
    # Sxc, Trend, Catalog_Corr = aligner(Flux_e, max_extinction, catalog)
    # Flux_Corr = flux / Trend[:, np.newaxis]

    # all stars (better at first glance)
    Sxc, Trend, Catalog_Corr = aligner(flux, max_extinction, catalog)
    if Sxc == 1:
        os._exit(0)
    Flux_Corr = flux / Trend[:, np.newaxis]

    ########################################################################
    # check result
    # plot std vs mag
    if save_figs:
        M = 20. - 2.5 * np.log10(Flux_Corr)
        S = np.nanstd(M, axis=0)
        plt.plot(catalog['B'], S, 'b.', label='exluded')
        plt.plot(catalog['B'][Catalog_Corr['ID']], S[Catalog_Corr['ID']], 'g.', label='in ensemble')
        plt.plot(catalog['B'][0], S[0], 'r*', label='target')
        plt.ylabel('std(mag)', fontsize=6)
        plt.xlabel('GAIA Bmag', fontsize=6)
        plt.tick_params(axis='both', labelsize=6, direction='in')
        plt.legend()
        plt.title('STD vs Mag')
        plt.grid()
        plt.show()
        plt.savefig(save_string_for_additions + '_std vs mag.pdf')

    # plot reference stars
    # normalize flux
    S_Flux = Flux_Corr[:, 1:] / np.nanmedian(Flux_Corr[:, 1:], axis=0)
    NEB_Flux = S_Flux
    IDs = catalog['ID'][1:]
    STD = np.nanstd(S_Flux, 0)

    # sort by std
    Sort = np.argsort(STD)
    S_Flux = S_Flux[:, Sort]
    IDs = IDs[Sort]
    plt.plot(time_list['JD'], Flux_Corr[:, 0] / np.median(Flux_Corr[:, 0]), 'r*',
             label='tess_object, BP=' + '{:.2f}'.format(catalog['B'][0]) +
                   ', STD = ' + '{:.4f}'.format(np.std(Flux_Corr[:, 0] / np.median(Flux_Corr[:, 0]))))

    # plot 6 best stars
    for i, Val in enumerate(IDs[0:7]):
        plt.plot(time_list['JD'], S_Flux[:, i] - 0.05 * (i + 1), '.',
                 label='Ref#' + str(Val) + ', BP=' + '{:.2f}'.format(catalog['B'][Val]) +
                       ', STD = ' + '{:.4f}'.format(np.std(S_Flux[:, i])))

    if save_figs:
        plt.ylabel('Normalized flux', fontsize=6)
        plt.axvspan(T0.jd, T1.jd, facecolor='k', alpha=0.2)
        locs, labels = plt.xticks()
        x_ticks_labels = []
        t = aTime(locs, format='jd')
        for x in t:
            x_ticks_labels.append(str(x.iso).split(' ')[1].split('.')[0])
        plt.xticks(locs[0:-1], x_ticks_labels[0:-1], rotation='vertical', fontsize=6)
        plt.xlabel('transit_date-time_list (UTC), ' + time_list['DATE-OBS'][0].split('T')[0], fontsize=6)
        plt.tick_params(axis='both', labelsize=6, direction='in')
        plt.gcf().set_size_inches(10, 4)
        plt.gca().set_position([0.1, 0.15, 0.65, 0.75])
        plt.legend(fontsize=6, bbox_to_anchor=(1.01, 1.))
        plt.title('Best reference stars')
        plt.grid()
        plt.show()

        plt.savefig(save_string_for_additions + '_Best reference stars.pdf')

    ########################################################################
    # start reports
    Title = tess_object
    Title += '\n'
    Title += f'{observatory}, ' + time_list['DATE-OBS'][0].split('T')[0]
    Title += ', filter=' + fil
    Title += ', Extime=' + '{:.1f}'.format(exp)
    Title += '\n'
    Title += 'Observation start=' + time_list['DATE-OBS'][0].split('.')[0]
    Title += ', '
    Title += 'end=' + time_list['DATE-OBS'][-1].split('.')[0]
    Title += ', ' + str(np.round((time_list['JD'][-1] - time_list['JD'][0]) * 24 * 60, 0)) + ' min.'
    Title += ', ' + str(len(time_list['JD'])) + ' frames'

    # zero level and dip
    Index = np.where((time_list['JD'] < T0.jd) | (time_list['JD'] > T1.jd))[0]
    Zero_Flux = sigma_clipped_stats(Flux_Corr[Index], sigma=3, maxiters=5,
                                    cenfunc=np.nanmedian, stdfunc=np.nanstd,
                                    axis=0)[1]
    # Depth = 1 / 10**(TTF['depth(mmag)'][0]/2500)

    # start output table
    Target_Report_Tbl = time_list.copy()
    Target_Report_Tbl.remove_columns(['EXPTIME', 'FILTER'])
    Target_Report_Tbl['X'] = np.round(Target_Report_Tbl['X'] -
                                      np.mean(Target_Report_Tbl['X']), 2)
    Target_Report_Tbl['Y'] = np.round(Target_Report_Tbl['Y'] -
                                      np.mean(Target_Report_Tbl['Y']), 2)
    Target_Report_Tbl['Rel_Flux_Object'] = np.round(Flux_Corr[:, 0] / Zero_Flux[0], 5)

    # calc error
    err = np.sqrt(flux[:, 0] + 3.14 * best_aperture ** 2 * (time_list['Sky'])) / (flux[:, 0])
    Target_Report_Tbl['Rel_Flux_Err_Object'] = np.round(err, 5)

    ZERO = int(Target_Report_Tbl['BJD'][0])

    ########################################################################
    print('start transit model fitting')
    transit_analysis = TransitAnalysis(name='tess_transit',
                                       # passbands='TESS',
                                       passbands=fil,
                                       times=Target_Report_Tbl['BJD'],
                                       fluxes=Target_Report_Tbl['Rel_Flux_Object'],
                                       # errors=Target_Report_Tbl['Rel_Flux_Err_Object'],
                                       # tm=RoadRunnerModel(),
                                       exptimes=exp / (60. * 60. * 24.)
                                       )
    transit_analysis.set_prior('p_1', 'NP', TTF['period(days)'], 1e-3)

    mid_time = (T1.jd + T0.jd) / 2.
    sigma = (T1.jd - T0.jd) / 4.
    transit_analysis.set_prior('tc_1', 'NP', mid_time, sigma)

    transit_analysis.optimize_global(niter=500, npop=100, use_tqdm=False)  # find posterior distribution
    transit_analysis.sample_mcmc(niter=2000, thin=20, repeats=3,
                                 save=False, use_tqdm=False)  # find posterior estimation

    samples = transit_analysis.sampler.flatchain
    max_prob = np.argmax(transit_analysis.sampler.flatlnprobability)  # индекс в цепи с наибольшей вероятностью
    popt = samples[max_prob]  # элемент цепи с наибольшей вероятностью
    model_flux = transit_analysis.flux_model(popt)  # модельный поток

    # print(popt)

    # параметры орбиты
    rho, tc, period, impact_parameter, radius_ratio2, seqw, sesw, q1, q2, wn_loge_0 = popt
    radius_ratio = np.sqrt(radius_ratio2)
    eccentricity = seqw ** 2 + sesw ** 2
    argument_of_periastron = np.arctan(sesw / seqw)
    smaxis = orbits_py.as_from_rhop(rho, period)
    inclination = orbits_py.i_from_ba(impact_parameter, smaxis)

    # время начала/конца транзита и дна
    t14 = orbits_py.d_from_pkaiews(period, radius_ratio, smaxis, inclination,
                                   eccentricity, argument_of_periastron, 1, kind=14)
    t23 = orbits_py.d_from_pkaiews(period, radius_ratio, smaxis, inclination,
                                   eccentricity, argument_of_periastron, 1, kind=23)
    t23 = np.float32(t23)
    # print('t23:', t23, 'type', type(t23))
    df = transit_analysis.posterior_samples()
    corner(df.posterior).savefig(report_string + '_corner_posterior.pdf')
    try:
        corner(df.derived_parameters).savefig(report_string + '_corner_derived_parameters.pdf')
    except:
        pass
    ########################################################################
    # BINNING

    actual_bin_size = int(len(Target_Report_Tbl['BJD']) / 10)
    if actual_bin_size > 20:
        actual_bin_size = 20
    if actual_bin_size >= min_bin_size:
        print('start binning')
        time_series = TimeSeries(time=Target_Report_Tbl['DATE-OBS'])
        time_series['flux'] = Target_Report_Tbl['Rel_Flux_Object']
        time_series_binned = aggregate_downsample(time_series=time_series, time_bin_size=actual_bin_size * u.min,
                                                  aggregate_func=np.nanmean)
        time_series_binned_std = aggregate_downsample(time_series=time_series, time_bin_size=actual_bin_size * u.min,
                                                      aggregate_func=np.nanstd)
        if not np.isnan(t23):  # если транзит v образный, то t23 не существует, значит глубину надо считать иначе
            indxs = np.where((tc + t23 * 0.5 > time_series_binned.time_bin_start.jd) &
                             (time_series_binned.time_bin_start.jd > tc - t23 * 0.5))[0]
        else:
            indxs = np.where((tc + t14 * 0.05 > time_series_binned.time_bin_start.jd) &
                             (time_series_binned.time_bin_start.jd > tc - t14 * 0.05))[0]
        transit_flux_depth = time_series_binned['flux'][indxs]

    else:
        if t23 is not np.nan:
            indxs = np.where((tc + t23 * 0.5 > Target_Report_Tbl['BJD']) &
                             (Target_Report_Tbl['BJD'] > tc - t23 * 0.5))[0]
        else:
            indxs = np.where((tc + t14 * 0.05 > Target_Report_Tbl['BJD']) &
                             (Target_Report_Tbl['BJD'] > tc - t14 * 0.05))[0]
        transit_flux_depth = Target_Report_Tbl['Rel_Flux_Object'][indxs]
        time_series_binned = None
        time_series_binned_std = None

    transit_depth_sigma_clip = sigma_clipped_stats(transit_flux_depth, sigma=3,
                                                   cenfunc=np.nanmedian, stdfunc=np.nanstd, axis=0)
    transit_depth_median = transit_depth_sigma_clip[1]
    transit_depth_sigma = transit_depth_sigma_clip[2]
    ########################################################################
    # start plots
    print("start plots")
    fig, axs = plt.subplots(3, 1, figsize=(6, 7), dpi=125)
    fig.suptitle(Title, fontsize=8)

    # flux vs JD
    pos = [0.125, 0.73, 0.8, 0.15]
    axs[0].set_position(pos)
    axs[0].plot(Target_Report_Tbl['BJD'] - ZERO, model_flux, 'y-', markersize=3, zorder=4,
                linewidth=1.2, label='model',
                linestyle='solid')  # рисуем модель
    axs[0].errorbar(Target_Report_Tbl['BJD'] - ZERO, Target_Report_Tbl['Rel_Flux_Object'],
                    Target_Report_Tbl['Rel_Flux_Err_Object'], fmt='b.',
                    label=r'1$\sigma$ errorbar')  # рисуем данные
    if time_series_binned is not None:
        axs[0].errorbar(time_series_binned.time_bin_start.jd - ZERO, time_series_binned['flux'],
                        time_series_binned_std['flux'], fmt='g.', markersize=9, zorder=4,
                        label=f'x{actual_bin_size} binned' + r', 1$\sigma$ errorbar')  # рисуем бины
    axs[0].axvspan(T0.tdb.jd - ZERO, T1.tdb.jd - ZERO, facecolor='k', alpha=0.2)  # продолжительность транзита
    # по изначальной информации
    axs[0].axhline(transit_depth_median, linewidth=0.5, color='r', linestyle='dashed')  # глубина транзита
    draw_my_annotate(axs[0], tc - t14 * 0.5 - ZERO)  # начало транзита по наблюдению
    if t23 is not np.nan:
        draw_my_annotate(axs[0], tc - t23 * 0.5 - ZERO)  # дно1
        draw_my_annotate(axs[0], tc + t23 * 0.5 - ZERO)  # дно2
    draw_my_annotate(axs[0], tc + t14 * 0.5 - ZERO)  # конец транзита по наблюдению
    axs[0].legend(loc=0, fontsize=6)
    axs[0].set_ylabel('Normalized flux', fontsize=6)
    locs = axs[0].get_xticks()
    t = aTime(locs, format='jd')
    x_ticks_labels = []
    for x in t:
        # x_ticks_labels.append(str(x.iso).split(' ')[1].split('.')[0])
        x_ticks_labels.append(str('{:.2f}'.format(x.tdb.jd)))
    axs[0].set_xticklabels(x_ticks_labels, fontsize=5)  # rotation='vertical',
    axs[0].set_xlabel(f'BJD - {ZERO}', fontsize=6)

    depth_ppt = np.round((1 - transit_depth_median) * 1000, 1)
    depth_ppt_sigma = np.round(transit_depth_sigma, 1)

    axs[0].tick_params(axis='both', labelsize=6, direction='in')
    axs[0].set_title(f'Target. Depth: {depth_ppt}ppt, '
                     f'std: {depth_ppt_sigma}. Model: impact_parameter: {np.round(impact_parameter, 1)},'
                     f'radius_ratio: {np.round(radius_ratio, 2)}, eccentricity: {np.round(eccentricity, 2)}, '
                     f'\nargument_of_periastron: {np.round(argument_of_periastron, 2)}, smaxis: {np.round(smaxis, 1)}, '
                     f'inclination: {np.round(inclination, 2)}',
                     loc='left', fontsize=6)
    axs[0].grid()

    # plot reference stars
    pos = [0.125, 0.35, 0.53, 0.31]
    axs[1].set_position(pos)
    axs[1].plot(time_list['BJD'] - ZERO, Flux_Corr[:, 0] / np.median(Flux_Corr[:, 0]), 'r*',
                label='tess_object, BP=' + '{:.2f}'.format(catalog['B'][0]) +
                      r', $\sigma$=' + '{:.4f}'.format(np.std(Flux_Corr[:, 0] / np.median(Flux_Corr[:, 0]))))
    # plot 6 best stars
    for i, Val in enumerate(IDs[0:7]):
        axs[1].plot(time_list['BJD'] - ZERO, S_Flux[:, i] - 0.05 * (i + 1), '.',
                    label='Ref#' + str(Val) + ', BP=' + '{:.2f}'.format(catalog['B'][Val]) +
                          r', $\sigma$=' + '{:.4f}'.format(np.std(S_Flux[:, i])))
    axs[1].legend(fontsize=6, bbox_to_anchor=(1.005, 1.))
    axs[1].set_ylabel('Normalized flux', fontsize=5)
    axs[1].axvspan(T0.tdb.jd - ZERO, T1.tdb.jd - ZERO, facecolor='k', alpha=0.2)
    locs = axs[1].get_xticks()
    t = aTime(locs, format='jd')
    x_ticks_labels = []
    for x in t:
        # x_ticks_labels.append(str(x.iso).split(' ')[1].split('.')[0])
        x_ticks_labels.append(str('{:.2f}'.format(x.tdb.jd)))
    axs[1].set_xticklabels(x_ticks_labels, fontsize=5)  # rotation='vertical',
    axs[1].set_xlabel(f'BJD-{ZERO}', fontsize=6)
    axs[1].tick_params(axis='both', labelsize=6, direction='in')
    axs[1].set_title('Reference stars', loc='left', fontsize=6)
    axs[1].grid()

    # plot errors
    pos = [0.125, 0.05, 0.8, 0.23]
    axs[2].set_position(pos)
    M = 20. - 2.5 * np.log10(Flux_Corr)
    S = np.nanstd(M, axis=0)
    axs[2].plot(catalog['B'], S, 'b.', label='exluded')
    axs[2].plot(catalog['B'][Catalog_Corr['ID']], S[Catalog_Corr['ID']], 'g.',
                label='in ensemble')
    axs[2].plot(catalog['B'][0], S[0], 'r*', label='tess_object')
    axs[2].set_ylabel('std(mag)', fontsize=6)
    axs[2].set_xlabel('GAIA Bmag', fontsize=6)
    axs[2].tick_params(axis='both', labelsize=6, direction='in')
    axs[2].legend(loc=2, fontsize=6)
    axs[2].set_title('std vs magnitudes', loc='left', fontsize=6)
    axs[2].set_ylim(0, 0.05)
    axs[2].grid()

    # plt.show()
    plt.savefig(report_string + '_report.pdf')

    ########################################################################
    # NEB
    fig1, ax = plt.subplots(1, 1, figsize=(6, 7), dpi=125)
    fig1.suptitle(Title, fontsize=8)
    # plot stars
    pos = [0.125, 0.07, 0.8, 0.8]
    ax.set_position(pos)
    k = 0

    NEB_Flux = NEB_Flux / np.median(NEB_Flux, axis=0)
    ax.plot(time_list['BJD'] - ZERO, NEB_Flux[:, 0],
            'r*', label='Target, RMS=' +
                        str(np.round(np.nanstd(NEB_Flux[:, 0]), 3)))
    for i in range(1, NEB_Flux.shape[1]):
        if k < 10:
            ax.plot(time_list['BJD'] - ZERO, NEB_Flux[:, i] - 0.05 * i,
                    '.', label='Ref#' + str(i) + ',' + \
                               #                 'Dist='+str(np.round(catalog['Dist'][i]*3600., 1))+'", '
                               'RMS=' + str(np.round(np.nanstd(NEB_Flux[:, i]), 3)))

            k = k + 1
    ax.legend(loc=2, fontsize=6)
    ax.set_ylabel('Normalized flux', fontsize=6)
    # ax.axvspan(T0.jd, T1.jd, facecolor='k', alpha=0.2)
    locs = ax.get_xticks()
    t = aTime(locs, format='jd')
    x_ticks_labels = []
    for x in t:
        # x_ticks_labels.append(str(x.iso).split(' ')[1].split('.')[0])
        x_ticks_labels.append(str('{:.2f}'.format(x.tdb.jd)))
    ax.set_xticklabels(x_ticks_labels, fontsize=5)  # rotation='vertical',
    ax.set_xlabel(f'BJD-{ZERO}', fontsize=6)
    ax.tick_params(axis='both', labelsize=6, direction='in')
    ax.set_title(str(k) + ' nearest stars, sorted by distance from target', loc='left', fontsize=6)
    ax.grid()
    plt.show()
    # if save_figs:
    plt.savefig(report_string + '_NEBs_checking.pdf')

    # report
    k = 0
    for i in range(0, len(Sort)):
        if k < 3:
            id = IDs[i]
            F = 'Rel_Flux_Ref#' + str(id)
            E = 'Rel_Flux_Err_Ref#' + str(id)
            Target_Report_Tbl[F] = np.round(S_Flux[:, i], 5)
            Target_Report_Tbl[E] = np.round(np.sqrt(flux[:, id] +
                                                    3.14 * best_aperture ** 2 * (time_list['Sky'])) / (flux[:, id]), 5)
        k = k + 1

    ascii.write(Target_Report_Tbl, report_string + '_lightcurve.dat',
                overwrite=True, delimiter='\t', fast_writer=False,
                format='commented_header')
    # scale = 3600 * np.sqrt(Header['CD1_1'] ** 2 + header['CD1_2'] ** 2)
    FS = np.argmax(catalog['R'])
    dt = time_list["DATE-OBS"][0].split("T")[0]
    dtr = dt.replace("-", "")
    raper = np.round(best_aperture * scale, 1)
    mean_fwhm = np.round(np.mean(time_list["SEXFWHM"]) * scale, 1)
    transit = what_transit(Target_Report_Tbl['JD'], T0.jd, T1.jd)
    eps = 5
    deltaT = np.round((T0.tdb.jd - (tc - t14 * 0.5)) / 24 / 60, 0)
    inTime = f"{np.abs(deltaT)} min. "
    # {"early" if  else "late"}

    if np.abs(deltaT) < eps:
        inTime = 'approx on time'
    elif deltaT >= eps:
        inTime = f"{np.abs(deltaT)} min. early"
    elif deltaT <= -eps:
        inTime = f"{np.abs(deltaT)} min. late"

    Target_Report_Txt = \
        f'''{"MASTER - Ural, 2x0.4m, 2xApogeeAltaU16m" if is_master else "RoboPhot, 0.6m, MicroLine ML4240"} cameras, observed a {tess_object} at {dt} in {fil} filter with exptime={np.round(exp, 1)} s.

Pixel scale = {np.round(scale, 2)}"/pix, mean FWHM of PSF is {mean_fwhm}". Photometric aperture radius is {best_aperture} pix or {raper}".
Duration of time series is {np.round((time_list["JD"][-1] - time_list["JD"][0]) * 24 * 60, 0)} and number of frames is {len(time_list["JD"])}.

Faintest star in FoV is #{catalog["ID"][FS]}. It is {np.round(catalog["R"][FS] - catalog["R"][0], 1)} mag(GAIA RP) weaker than the target.
Previous TTF comments: {TTF["comments"][0]}

Tag: {dtr}_{observer}_{observatory}_

Model parameters: 
* The model is built in the PyTransit 
* Radius ratio: {np.round(radius_ratio, 4)}
* Orbital semi-major axis divided by the stellar radius: {np.round(smaxis, 4)}
* Orbital inclination: {np.round(inclination, 4)}
* Orbital eccentricity: {np.round(eccentricity, 4)}
* Argument of periastron: {np.round(argument_of_periastron, 4)}

{tess_object} at UT{dtr} from Kourovka Observatory {"0.4m (x2)" if is_master else "0.6m"} in {fil}.

Nikita Chazov/Kourovka Observatory {"0.4m (x2)" if is_master else "0.6m"} observed {transit} on {dtr} in {fil} and detected a {inTime} {depth_ppt}ppt using {raper}" target aperture. 
1. Aperture radius = {raper}"
2. Typical FWHM = {mean_fwhm}"
3. Predicted Tc = {np.round((T1.tdb.jd + T0.tdb.jd) / 2 - 2450000, 4)}
'''
    if transit == 'a full transit' or transit == 'a full with gapped':
        Target_Report_Txt += f'4. Measured Tc = {np.round(tc - 2450000, 4)}'

    df = open(report_string + '_report_note.txt', 'w')
    df.write(Target_Report_Txt)
    df.close()
