##from astropy.time import Time
import requests
from astropy.io import ascii


def create_TTF_entry_url(date_obs, priority, min_depth, max_mag, target):
    url = 'https://astro.swarthmore.edu/telescope/tess-secure/' + \
          'print_eclipses.cgi?observatory_string=57.036537%3B' + \
          '59.545735%3BAsia%2FYekaterinburg%3BKourovka+' + \
          'Observatory+0.4m+%28x2%29%3BKourovka+Obs+0.4m+' + \
          '%28x2%29&use_utc=1&observatory_latitude=57.036537' + \
          '&observatory_longitude=59.545735&timezone=UTC' + \
          '&start_date={}&days_to_print=1&days_in_past=0' + \
          '&minimum_start_elevation=30' + \
          '&and_vs_or=or&minimum_end_elevation=30' + \
          '&minimum_ha=-12&maximum_ha=12&baseline_hrs=0' + \
          '&show_unc=1&maximum_priority={}&minimum_depth={}' + \
          '&maximum_V_mag={}&target_string={}' + \
          '&lco_only=0&single_object=0&ra=&dec=&epoch=' + \
          '&period=&duration=&target=&show_ephemeris=0' + \
          '&print_html=2&twilight=-12&max_airmass=2.4'

    return (url.format(date_obs, priority, min_depth, max_mag, target))


def TTF_Query(DateObs, priority, min_depth, max_mag, Target):
    TTF_URL = create_TTF_entry_url(DateObs, priority, min_depth, max_mag, Target)
    ssn = requests.session()
    print(TTF_URL)
    ssn.auth = ('tess_nda_observer', 'F1nd_TE$S_PlaNets!')
    req = ssn.get(TTF_URL)
    if req.status_code != 200:
        message = None
    else:
        ##        print(req.text)
        ##        Cleared = req.text
        ##        lines = req.text.splitlines()
        ##        for line in lines:
        ##            print(line)
        ##            print()
        ##            Index = [i for i, ltr in enumerate(line) if ltr == '"']
        ##            if len(Index)>0:
        ##                Comment = '"'+line[Index[0]+1:Index[-1]].replace('"', 'arcsec')+'"'
        ##                line = line[:Index[0]]+Comment
        ##            Cleared =  Cleared + line + '\n'
        message = ascii.read(req.text, delimiter=',',
                             format='commented_header', quotechar='"')  # , \
    # delimiter = ',', format = 'commented_header')
    ##    print(message.info)

    return message

##print(TTF_Query('2022-08-28', 5, 0.5, 16, 'TIC+427153307.01'))
