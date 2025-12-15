import sauna

# From Excel
covariances = sauna.Covariances()
covariances.library = 'ENDF-B-VII.1'
covariances.group_structure = sauna.GROUP_STRUCTURES['SCALE-56']
covariances.from_excels(f'../../covariances/ENDF-B-VII.1-56/')

# From SCALE
sensitivities_1 = sauna.Sensitivities()
sensitivities_1.group_structure = sauna.GROUP_STRUCTURES['SCALE-56']
sensitivities_1.from_scale('../../models/MET1000_FC_NFP.2025.02.28T10.09.44.sdf', 'B', 'Eigenvalue')

sensitivities_2 = sauna.Sensitivities()
sensitivities_2.group_structure = sauna.GROUP_STRUCTURES['SCALE-56']
sensitivities_2.from_scale('../../models/Jezebel_NFP.2025.03.07T19.12.19.sdf', 'B', 'Eigenvalue')


E11    = sauna.Analysis.compare(sensitivities_1, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'E')
E12    = sauna.Analysis.compare(sensitivities_1, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'E')
E21    = sauna.Analysis.compare(sensitivities_2, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'E')
E22    = sauna.Analysis.compare(sensitivities_2, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'E')
G11    = sauna.Analysis.compare(sensitivities_1, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G')
G12    = sauna.Analysis.compare(sensitivities_1, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G')
G21    = sauna.Analysis.compare(sensitivities_2, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G')
G22    = sauna.Analysis.compare(sensitivities_2, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G')
G11_x  = sauna.Analysis.compare(sensitivities_1, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G', reactions = [2,4,18,102,103,104,105,106,107])
G12_x  = sauna.Analysis.compare(sensitivities_1, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G', reactions = [2,4,18,102,103,104,105,106,107])
G21_x  = sauna.Analysis.compare(sensitivities_2, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G', reactions = [2,4,18,102,103,104,105,106,107])
G22_x  = sauna.Analysis.compare(sensitivities_2, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G', reactions = [2,4,18,102,103,104,105,106,107])
c_k11  = sauna.Analysis.compare(sensitivities_1, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'c_k', covariances = covariances)
c_k12  = sauna.Analysis.compare(sensitivities_1, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'c_k', covariances = covariances)
c_k21  = sauna.Analysis.compare(sensitivities_2, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'c_k', covariances = covariances)
c_k22  = sauna.Analysis.compare(sensitivities_2, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'c_k', covariances = covariances)

print('----ck----------------E---------------------G---------')
print('----m1------m2--------m1---------m2---------m1---------m2')
print(round(c_k11,3),     round(c_k12,3),     round(E11,3),       round(E12,3),     round(G11_x,3), round(G12_x,3))
print(round(c_k21,3),     round(c_k22,3),     round(E21,3),       round(E22,3),     round(G21_x,3), round(G22_x,3))
