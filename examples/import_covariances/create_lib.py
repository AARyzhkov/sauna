# This example generates .xlsx covariances in the 28-group ABBN and
# 56-group SCALE approximations, and it saves it into the coviariances
# folder in the root folder. The data was already saved in the folder.

import sauna

libraries = ['ENDF-B-VII.1']
exts = ['.dat']

for lib, ext in zip(libraries, exts):

    covariances = sauna.Covariances()
    covariances.library = lib
    covariances.group_structure = sauna.GROUP_STRUCTURES['SCALE-56']

    if (lib == 'ENDF-B-VII.1') | (lib == 'ENDF-B-VIII.0')| (lib== 'ENDF-B-VIII.1')  | (lib == 'JEFF-3.3'):
        folder = f'../../../NuclearData/ENDF_Libraries/{lib}/n/'
    else:
        folder = f'../../../NuclearData/ENDF_Libraries/{lib}/'
    
    covariances.from_endfs(folder, extension=ext, parallel=True)
    covariances.to_excels(f'../../covariances/{lib}-56/')

    covariances = sauna.Covariances()
    covariances.library = lib
    covariances.group_structure = sauna.GROUP_STRUCTURES['ABBN-28']

    if (lib == 'ENDF-B-VII.1') | (lib == 'ENDF-B-VIII.0')| (lib== 'ENDF-B-VIII.1')  | (lib == 'JEFF-3.3'):
        folder = f'../../../NuclearData/ENDF_Libraries/{lib}/n/'
    else:
        folder = f'../../../NuclearData/ENDF_Libraries/{lib}/'
    
    covariances.from_endfs(folder, extension=ext, parallel=True)
    covariances.to_excels(f'../../covariances/{lib}-28/')