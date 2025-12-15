######################################################################################
# This example intends to reproduce the results that were obtained via the
# SCALE code system via a Perfetti-like approach, presented in the following work:
# Ryzhkov, A.A., Tikhomirov, G.V., Ternovykh, M.Yu., Gerasimov, A.S., 2023.
# Evaluation of technological uncertainties using the sensitivity to nuclear data.
# Atomic Energy 135, 162–165. https://doi.org/10.1007/s10512-024-01103-w,
# but with the Serpent sensitivity results. Athough the scheme here is simplified 
# due to having only total sensitivities over the whole core model,
# the results are in excellent agreement.
#######################################################################################

import sauna

# From Serpent
sensitivities = sauna.Sensitivities()
sensitivities.group_structure = sauna.GROUP_STRUCTURES['ABBN-28']
sensitivities.from_serpent('../../models/MET1000_FC_28_sens0.m')

# Uncertainty in atomic fraction of Pu239 in Pu
pu_nuclides = [942360, 942380, 942400, 942410, 942420]
a_1 = 0.5390 # inner core
a_2 = 0.5133 # outer core
ao_239 = (a_1*78 + a_2*102)/180 # volume-averaged ao over the whole core
uncertainty_ao = sauna.Analysis.get_concentration_uncertainty(sensitivities    = sensitivities,
                                                              uncertainty      = 0.01,
                                                              targets          = [942390],
                                                              background_zams  = pu_nuclides,
                                                              fraction         = ao_239,
                                                              fraction_type    = 'ao',
                                                              uncertainty_type = 'normal')

print(f'The influence of uncertainty in atomic fraction of Pu239 in Pu is: {uncertainty_ao}')

# Uncertainty in weight fraction of Pu in the heavy components
heavy_nuclides = [922340, 922350, 922360, 922380, 932370, 942360, 942380,  
                 942400, 942410, 942420, 952410, 952421, 952430, 962420,
                 962430, 962440, 962450, 962460]
target_nuclides = pu_nuclides
w_1    = 0.09311 # inner core
w_2    = 0.11330 # outer core
wo_239 = (w_1*78 + w_2*102)/180 # volume-averaged wo over the whole core
uncertainty_wo = sauna.Analysis.get_concentration_uncertainty(sensitivities    = sensitivities,
                                                              uncertainty      = 0.01,
                                                              targets          = target_nuclides,
                                                              background_zams  = heavy_nuclides,
                                                              fraction         = wo_239,
                                                              fraction_type    = 'wo',
                                                              uncertainty_type = 'normal')

print(f'The influence of uncertainty in weight fraction of Pu in the heavy component is: {uncertainty_wo}')

