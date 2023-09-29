import urllib3
import requests
from bs4 import BeautifulSoup
from time import sleep
import os
from Utils import target_to_tic_number

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


def download_file(url, path2save):
    http = urllib3.PoolManager(cert_reqs='CERT_NONE')
    header = urllib3.util.make_headers(basic_auth='googlov:qwertygooglov!')
    r = http.request('GET', url, headers=header)
    data = r.data
    print("Ok")
    print("downloading " + url)
    with open(path2save + '/' + os.path.basename(url), "wb") as local_file:
        local_file.write(data)
        print('save to', path2save)


def get_links(soup: BeautifulSoup):
    list_urls = []
    for i in soup.findAll("a"):
        link = str(i.get("href"))
        list_urls.append(link)
    return list_urls[1:]


def create_entry_url(target: str):
    url = f'https://master.kourovka.ru/data/master/imdata/procdata/trimmed/{target}/'
    return url.format()


def obs_query_from_master(target: str, path2save=r'E:/tess/master/'):
    target = ''.join(['TIC', *target_to_tic_number(target)])
    url = create_entry_url(target)
    ssn = requests.session()
    ssn.auth = ('googlov', 'qwertygooglov!')
    req = ssn.get(url)
    sleep(0.1)

    if req.status_code != 200:
        return None

    print('Successful request to observatory server')
    # print(req.text)
    soup = BeautifulSoup(req.text, "lxml")
    target_obs = get_links(soup)[-1]
    req = ssn.get(url + target_obs)
    sleep(0.1)
    soup = BeautifulSoup(req.text, "lxml")
    print('Get filters')
    fil = get_links(soup)
    print(f'Filters: {fil}')
    paths = []
    for f in fil:
        req = ssn.get(url + target_obs + f)
        print('Get frames from filter ', f[:-1])
        sleep(0.1)
        soup = BeautifulSoup(req.text, "lxml")
        frames = (get_links(soup))
        if not os.path.exists(path2save + target + '/' + target_obs + '/' + f):
            os.makedirs(path2save + target + '/' + target_obs + '/' + f)
        for frame in frames:
            download_file(url=url + target_obs + f + frame, path2save=path2save + target + '/' + target_obs + '/' + f)
        paths.append(path2save + target + '/' + target_obs + '/' + f[:-1])
    return paths
