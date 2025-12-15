import sauna

# From Serpent
sensitivities = sauna.Sensitivities()
sensitivities.group_structure = sauna.GROUP_STRUCTURES['ABBN-28']
sensitivities.from_serpent('../../models/MET1000_FC_28_sens0.m')

# Nuclides in the fuel
fuel_nuclides = [922340, 922350, 922360, 922380, 932370, 942360, 942380,  
                 942400, 942410, 942420, 952410, 952421, 952430, 962420,
                 962430, 962440, 962450, 962460,
                 400900, 400910, 400920, 400940, 400960,
                 420920, 420940, 420950, 420960, 420970, 420980, 421000]

uncertainty_den = sauna.Analysis.get_density_uncertainty(sensitivities    = sensitivities,
                                                        uncertainty      = 0.01,
                                                        targets          = fuel_nuclides,
                                                        uncertainty_type = 'normal')

print(f'The uncertainty influence of density is: {uncertainty_den}')

