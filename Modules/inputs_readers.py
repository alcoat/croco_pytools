import numpy as np
import sys

def lookvar(input):

    if input == 'mercator':
        dico={ 'depth':'depth',\
               'lonr':'longitude','lonu':'longitude','lonv':'longitude',\
               'latr':'latitude','latu':'latitude','latv':'latitude',\
               'ssh':'zos',\
               'temp':'thetao',\
               'salt':'so',\
               'u': 'uo',\
               'v': 'vo',\
               'time': 'time',\
               'time_dim':'time'\
             }

    elif input == 'eccov4':
        dico={ 'depth':'dep',\
               'lonr':'lon','lonu':'lon','lonv':'lon',\
               'latr':'lat','latu':'lat','latv':'lat',\
               'ssh':'ETAN',\
               'temp': 'THETA',\
               'salt': 'SALT',\
               'u': 'EVEL',\
               'v': 'NVEL',\
               'time':'tim',\
               'time_dim':'i1'\
             }

    elif input == 'soda':
        dico={ 'depth':'st_ocean',\
               'lonr':'xt_ocean','lonu':'xu_ocean','lonv':'xu_ocean',\
               'latr':'yt_ocean','latu':'yu_ocean','latv':'yu_ocean',\
               'ssh':'ssh',\
               'temp': 'temp',\
               'salt': 'salt',\
               'u': 'u',\
               'v': 'v',\
               'time':'time',\
               'time_dim':'time'\
         }
#    elif input == 'croco':
#        dico={ 'ssh':'ssh',\
#               'temp': 'temp',\
#               'salt': 'salt',\
#               'u': 'u',\
#               'v': 'v'\
#         }
    else:
        sys.exit(print('No \'%s\' dico available in Modules/dico.py. Please add it an run the script' % input))

    return dico



