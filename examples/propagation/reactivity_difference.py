import sauna

covariances = sauna.Covariances()
covariances.group_structure = sauna.GROUP_STRUCTURES['SCALE-56']
covariances.from_excels(f'../../covariances/ENDF-B-VII.1-56/')

sensitivities_1 = sauna.Sensitivities()
sensitivities_1.group_structure = sauna.GROUP_STRUCTURES['SCALE-56']
sensitivities_1.from_scale('../../models/MET1000_FC_NFP.2025.02.28T10.09.44.sdf', 'B', 'Eigenvalue')

sensitivities_2 = sauna.Sensitivities()
sensitivities_2.group_structure = sauna.GROUP_STRUCTURES['SCALE-56']
sensitivities_2.from_scale('../../models/MET1000_FC_NFP.2025.03.02T00.03.35.sdf', 'B', 'Eigenvalue')


# 1.016566, 1.023644, 0.000058, 0.000058 - SCALE 
# 1.01657, 1.02378, 0.00006, 0.00006 - Serpent 
sensitivities_3 = sauna.Sensitivities.reactivity_difference(sensitivities_1, sensitivities_2, 1.016566, 1.023644, 0.000058, 0.000058, functional = 'Void reactivity')

dfs = sauna.Analysis.get_breakdown(sensitivities_3, covariances)
print(dfs)