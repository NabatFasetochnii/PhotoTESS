from Core import photometry
from download import obs_query
from Postproc import Post_Proc

observer = 'chazov'
# observer = 'sysolina'
# observer = 'krushinsky'

# paths = obs_query('TIC417129824_01')
# print('paths - ', paths)
# paths = [r'G:\tess\TIC14278542101\4280\EAST_V', r'G:\tess\TIC14278542101\4280\WEST_R']
paths = [r'E:\tess\TIC41712982401\4302\WEST_R']

for path in paths:
    photometry(path2data=path, observer=observer)

# Post_Proc.post_proc(apers=[3, 4, 5], observer=observer, path2data=r'E:\tess\TIC41712982401\4302\EAST_V\Photometry')
