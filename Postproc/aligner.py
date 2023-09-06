import numpy as np


# return 50 stars with best colors and magnitudes
def get_ensemble(Cat):
    Cat['BR'] = Cat['B'] - Cat['R']
    Gmag = Cat['B'][0]
    BR = Cat['BR'][0]
    Cat['B'] = abs(Cat['B'] - Gmag)
    Cat['BR'] = abs(Cat['BR'] - BR)
    Cat.sort(['BR', 'B'])
    if len(Cat) > 50:
        Cat = Cat[0:50]
    Cat.sort('ID')
    return Cat


def aligner(Flux, Max_extinction, Catalog):
    Catalog_Corr = Catalog.copy()
    #    ##Search bad frames
    #    Max_extinction = 1/(10**(Max_extinction/2.5))
    #    Sum_Flux = np.sum(Flux, axis=1)
    #    Sum_Flux = Sum_Flux/np.median(Sum_Flux)
    #    ##delete bad frames
    #    Index = np.where(Sum_Flux>Max_extinction)[0]
    #      plt.plot(Sum_Flux, 'r')
    #      plt.plot(Sum_Flux[Index], 'g')
    #      plt.hlines(Max_extinction, 0, len(Sum_Flux), 'r', label = 'Max extinction')
    #      plt.xlabel('Frame')
    #      plt.ylabel('Normalized flux')
    #      plt.legend()
    #      plt.grid()
    #      plt.show()
    #    Flux = Flux[Index]
    #    Time.remove_rows(Index)

    # delete bad stars
    Index = np.where(np.isnan(Flux)[0])
    if 0 in Index:
        print('Object has NaNs!')
        return 1, 0, 0
    Flux = np.delete(Flux, Index, axis=1)
    Catalog_Corr.remove_rows(Index)

    # main circle
    flag = 0
    while flag != 1:
        Ensemble = Flux[:, 1:]
        Trend = np.sum(Ensemble, axis=1)
        Trend = Trend / np.mean(Trend)
        Ensemble = Ensemble / Trend[:, np.newaxis]

        # find and remove the worst star
        Std = np.std(Ensemble, 0) / np.sqrt(np.mean(Ensemble))
        if np.max(Std) > 3:
            Index = np.argmax(Std) + 1
            print('Delete star #', Catalog_Corr['ID'][Index],
                  ' with STD=', np.max(Std))
            #            plt.plot(Ensemble)
            #            plt.plot(Ensemble[:, Index-1], 'r--', label = 'worst star')
            #            plt.plot(Flux[:, 0]/ Trend, 'r.', label = 'object')
            #            plt.legend()
            #            plt.show()

            Flux = np.delete(Flux, Index, axis=1)
            Catalog_Corr.remove_rows(Index)
        else:
            print('Stars in ensemble: ', Flux.shape[1])
            flag = 1
    return 0, Trend, Catalog_Corr
