---
icon: percent
---

# Density

The uncertainty propagation in the material density assumes that all the components are perturbed. The first step after importing the package and sensitivities is to provide the composition of the material. This example considers the fuel region. Since the fuel is metallic, the alloying material, $$\text{Zr}$$ has to be accounted. It should be noted that the nuclides have to correspond to only one material of interest to get a correct solution.

```python
import sauna

sensitivities = sauna.Sensitivities()
sensitivities.from_serpent('../../models/MET1000_FC_28_sens0.m')

fuel_nuclides = [922340, 922350, 922360, 922380, 932370, 942360, 942380,  
                 942400, 942410, 942420, 952410, 952421, 952430, 962420,
                 962430, 962440, 962450, 962460,
                 400900, 400910, 400920, 400940, 400960,
                 420920, 420940, 420950, 420960, 420970, 420980, 421000]
```

The following function takes the list of nuclides in the material of interest and the uncertainty of the density, assumed to be equal to 1%. One may also provide the type of the uncertainty distribution {'normal', 'uniform', 'triangular'} with 'normal' as a default.

```python
uncertainty_den = sauna.Analysis.get_density_uncertainty(sensitivities    = sensitivities,
                                                        uncertainty      = 0.01,
                                                        targets          = fuel_nuclides,
                                                        uncertainty_type = 'normal')
```

The influence is moderate and equals 0.026%.
