---
icon: cookie
---

# Atomic Fraction of Pu239 in Pu

This page is aimed to demonstrate how to compute the uncertainty influence in the atomic fraction of $$^{239}\text{Pu}$$ in the $$\text{Pu}$$ composition.

```python
import sauna
```

The implementation relies upon the sensitivities of a functional to the total cross section. Consequently, it is necessary to get them from the corresponding neutron transport simulation results. In this example, it is the sensitivities of MET1000 obtained via Serpent.&#x20;

```python
sensitivities = sauna.Sensitivities()
sensitivities.from_serpent('../../models/MET1000_FC_28_sens0.m')
```

After importing the data, one has to state the nuclides that are used in the normalization of the total atomic fraction to unity. In this case, they are the $$\text{Pu}$$ isotopes, but the ones, we are interested in, have to be excluded, that is, $$^{239}\text{Pu}$$.

```python
pu_nuclides = [942360, 942380, 942400, 942410, 942420]
```

Then, the atomic fraction of $$^{239}\text{Pu}$$ has to be provided for perturbing. In the example the value is about 0.5244:

```python
ao_239 = 0.5244
```

Having the values, it is possible now to call the function to assess the impact of the uncertainty, assuming its value is equal to 1%.

```python
uncertainty_ao = sauna.Analysis.get_concentration_uncertainty(sensitivities    = sensitivities,
                                                              uncertainty      = 0.01,
                                                              targets          = [942390],
                                                              background_zams  = pu_nuclides,
                                                              fraction         = ao_239,
                                                              fraction_type    = 'ao',
                                                              uncertainty_type = 'normal')
```

The expected output is about 0.31% for the eigenvalue.
