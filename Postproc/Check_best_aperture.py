import numpy as np
import os

from astropy.io import ascii
from Postproc.aligner import get_ensemble, aligner
import matplotlib.pyplot as plt

# disable warnings
import warnings

warnings.simplefilter("ignore")

# #####################################################################3897


def check_best_aperture(apers, Path2Data=r'G:\tess\good\sent\TIC22951086601\3988\EAST_V\Photometry', save_fig=False):
    Max_extinction = 0.35

    ########################################################################
    # read data
    Catalog = ascii.read(Path2Data + '/Cat.txt')
    # Time = ascii.read(Path2Data + '/Time.txt')

    Phot_list = []
    for f in os.listdir(Path2Data):
        if f.count('Phot'):
            Phot_list.append(Path2Data + '/' + f)

    ########################################################################
    Ensemble = get_ensemble(Catalog)

    plt.switch_backend('pdf')  # для фикса какой-то тупой ошибки в нарнии pyplot
    fig, axs = plt.subplots(3, 1, figsize=(8, 6), dpi=125)
    CStar = np.zeros(3, dtype='int')

    check_list = []
    for f in Phot_list:
        print(f.split('/')[-1].split('.')[0])
        Flux = np.genfromtxt(f)
        Flux_e = Flux[:, Ensemble['ID']]
        Sxc, Trend, Catalog_Corr = aligner(Flux_e, Max_extinction, Catalog)
        Flux = Flux / Trend[:, np.newaxis]
        NFlux = Flux / np.mean(Flux, axis=0)
        line, = axs[0].plot(NFlux[:, 0],
                            label=f.split('/')[-1].split('.')[0] +
                                  ', std = ' + '{:.4f}'.format(np.std(Flux[:, 0] / np.mean(Flux[:, 0]))))
        if np.sum(CStar) == 0:
            CStar = np.argsort(np.std(NFlux, axis=0))[:3]
        for i, S in enumerate(CStar):
            line, = axs[1].plot(NFlux[:, S] - i * 0.02, c=line.get_color())

        line.set_label(f.split('/')[-1].split('.')[0] +
                       ', std: ' + '{:.4f}'.format(np.std(NFlux[:, CStar[0]])) +
                       ', {:.4f}'.format(np.std(NFlux[:, CStar[1]])) +
                       ', {:.4f}'.format(np.std(NFlux[:, CStar[2]])))

        Mag = 20 - 2.5 * np.log10(np.mean(Flux, axis=0))
        Err = np.std(Flux, axis=0) / np.mean(Flux, axis=0)
        axs[2].scatter(Mag, Err, marker='.', color=line.get_color(),
                       label=f.split('/')[-1].split('.')[0])
        axs[2].scatter(Mag[0], Err[0], marker='*', color=line.get_color(),
                       label=f.split('/')[-1].split('.')[0])
        check_list.append(np.sum(Err[CStar]**2) + Err[0]**2)
        print()

    axs[0].set_position([0.07, 0.7, 0.65, 0.25])
    axs[0].set_title('Object', fontsize=6, loc='right')
    axs[0].set_ylabel('normalized flux', fontsize=6)
    axs[0].set_xlabel('frame', fontsize=6)
    axs[0].tick_params(axis='both', labelsize=6, direction='in')
    axs[0].grid()
    axs[0].legend(fontsize=6, bbox_to_anchor=(1.01, 1.))

    axs[1].set_position([0.07, 0.38, 0.65, 0.25])
    axs[1].set_title('Check stars: ' + str(CStar), fontsize=6, loc='right')
    axs[1].set_ylabel('normalized flux', fontsize=6)
    axs[1].set_xlabel('frame', fontsize=6)
    axs[1].tick_params(axis='both', labelsize=6, direction='in')
    axs[1].grid()
    axs[1].legend(fontsize=6, bbox_to_anchor=(1.01, 1.))

    axs[2].set_position([0.07, 0.05, 0.65, 0.25])
    axs[2].set_title('error vs mag', fontsize=6, loc='right')
    axs[2].set_ylabel('error', fontsize=6)
    axs[2].set_xlabel('instrumental mag', fontsize=6)
    axs[2].set_ylim(0, 0.07)
    axs[2].tick_params(axis='both', labelsize=6, direction='in')
    axs[2].grid()
    axs[2].legend(fontsize=6, bbox_to_anchor=(1.01, 1.))

    plt.show()
    if save_fig:
        plt.savefig(Path2Data + '/check_best_aperture.pdf')
    return apers[np.argmin(check_list)]
