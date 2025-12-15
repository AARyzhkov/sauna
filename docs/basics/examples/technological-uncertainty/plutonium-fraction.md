---
icon: weight-hanging
---

# Plutonium Fraction

To assess the effect of the $$\text{Pu}$$ fraction among the heavy component influences on the functional a similar approach to the previous one can be applied. The same may also be done for the $$^{235}\text{U}$$ uranium enrichment. Firstly, a part of the previous example is reproduced.

```python
import sauna

sensitivities = sauna.Sensitivities()
sensitivities.from_serpent('../../models/MET1000_FC_28_sens0.m')
pu_nuclides = [942360, 942380, 942400, 942410, 942420]
```

Next, heavy nuclides are provided as only they are accounted when fissile fraction is considered. It is also widely admitted that this fraction is a weight one, and the value has to be calculated beforehand.

```python
heavy_nuclides = [922340, 922350, 922360, 922380, 932370, 942360, 942380,  
                 942400, 942410, 942420, 952410, 952421, 952430, 962420,
                 962430, 962440, 962450, 962460]
wo_239 = 0.10455
```

It is useful to note that here the list of heavy nuclides contains the nuclides of interest, but it should not as it is mentioned in [the previous example](atomic-fraction-of-pu239-in-pu.md). As a result, SAUNA automatically removes redundant nuclides. Having the inputs prepare, one calls a function, explicitly mentioning that the fraction type is `'wo'`. In this example, the input uncertainty is assumed 1%, the same as in the previous example.

```python
uncertainty_wo = sauna.Analysis.get_concentration_uncertainty(sensitivities    = sensitivities,
                                                              uncertainty      = 0.01,
                                                              targets          = target_nuclides,
                                                              background_zams  = heavy_nuclides,
                                                              fraction         = wo_239,
                                                              fraction_type    = 'wo',
                                                              uncertainty_type = 'normal')
```

The expected output is about 0.26% for the eigenvalue.
