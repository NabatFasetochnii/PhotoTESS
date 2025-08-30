from Core import photometry
from download import obs_query_from_master
from Postproc import Post_Proc

observer = 'chazov'
# observer = 'sysolina'
# observer = 'krushinsky'

# paths = obs_query_from_master('TIC250803548_01')
paths = [r'E:\tess\robophot\2023_09_08_TIC233047097_01\r\Photometry',
         r'E:\tess\robophot\2023_09_08_TIC233047097_01\i\Photometry']
print('paths - ', paths)

# paths = [r'E:\tess\master\good\sent\TIC32249810101\4325\EAST_I\Photometry']

for path in paths:
    # photometry(path2data=path, RAper=[3, 4, 5, 6, 7])
# photometry(r'E:\tess\TIC32681580411\4306\WEST_R', observer=observer)
    Post_Proc.post_proc(apers=[3, 4, 5, 6, 7], observer=observer,
                        path2data=path,
                        is_master=False, scale=1.3)
