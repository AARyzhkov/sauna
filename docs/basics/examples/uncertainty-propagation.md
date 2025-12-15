---
icon: waves-sine
---

# Uncertainty Propagation

This subsection shows a minimum example how to propagate the nuclear data uncertainty via SAUNA. The example is provided in `./examples/propagation` from the SAUNA root.

Firstly, one shall import the package.

```python
import sauna
```

Next, one has to get sensitivities. These are taken from Serpent results for MET1000.

```python
serp_sensitivities = sauna.Sensitivities()
serp_sensitivities.group_structure = sauna.GROUP_STRUCTURES['SCALE-56']
serp_sensitivities.from_serpent('../../models/MET1000_FC_56_sens0.m')
```

After that, covariances are necessary to propagate them on the functionals.

```python
covariances = sauna.Covariances()
covariances.library = 'ENDF-B-VII.1'
covariances.group_structure = sauna.GROUP_STRUCTURES['SCALE-56']
covariances.from_excels(f'../../covariances/{covariances.library}-56/')
```

Finally, combining the sensitivities and covariances yields the functional uncertainties and their sources as a dataframe in the same way as it is presented in [a previous section](../user-guide/uncertainty-analysis/uncertainty-propagation.md).

```python
serp_dfs  = sauna.Analysis.get_breakdown(serp_sensitivities, covariances)
```

