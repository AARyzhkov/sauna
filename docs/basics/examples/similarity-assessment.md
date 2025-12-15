---
icon: reflect-both
---

# Similarity Assessment

This subsection presents a similarity quantification via similarity indices. The examples are provided in `'./examples/similarity'` from the root directory of SAUNA.

Quantifying how two different models are similar via similarity indices is possible to do via the SAUNA capabilities. Firstly, SAUNA has to be imported.&#x20;

```python
import sauna
```

For similarity assessment, one needs two sensitivities instances, containing sensitivities for two different models though the same models can be used to verify whether the capability works correctly expecting unity for an index. This example uses the eigenvalue sensitivities of MET1000 and Jezebel obtained via the SCALE code system.

```python
sensitivities_1 = sauna.Sensitivities()
sensitivities_1.from_scale('../../models/MET1000_FC_NFP.2025.02.28T10.09.44.sdf', 'B', 'Eigenvalue')

sensitivities_2 = sauna.Sensitivities()
sensitivities_2.from_scale('../../models/Jezebel_NFP.2025.03.07T19.12.19.sdf', 'B', 'Eigenvalue')
```

Now, it is possible to compute the $$E$$ and $$G$$ indices in a mass manner with different combinations.

```python
E11    = sauna.Analysis.compare(sensitivities_1, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'E')
E12    = sauna.Analysis.compare(sensitivities_1, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'E')
E21    = sauna.Analysis.compare(sensitivities_2, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'E')
E22    = sauna.Analysis.compare(sensitivities_2, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'E')
G11    = sauna.Analysis.compare(sensitivities_1, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G')
G12    = sauna.Analysis.compare(sensitivities_1, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G')
G21    = sauna.Analysis.compare(sensitivities_2, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G')
G22    = sauna.Analysis.compare(sensitivities_2, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G')
```

The $$G$$ index may be calculated with a number of chosen reactions. For instance, accounting the same reactions as ones in SCALE can be done in the following way, by providing a list of the reactions of interest as an argument.

```python
G11_x  = sauna.Analysis.compare(sensitivities_1, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G', reactions = [2,4,18,102,103,104,105,106,107])
G12_x  = sauna.Analysis.compare(sensitivities_1, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G', reactions = [2,4,18,102,103,104,105,106,107])
G21_x  = sauna.Analysis.compare(sensitivities_2, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G', reactions = [2,4,18,102,103,104,105,106,107])
G22_x  = sauna.Analysis.compare(sensitivities_2, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'G', reactions = [2,4,18,102,103,104,105,106,107])
```

The final index, $$c_k$$ requires an additional argument — covariances. Hence, one has to provide them to calculate this index.

```python
covariances = sauna.Covariances()
covariances.from_excels(f'../../covariances/ENDF-B-VII.1-56/')
```

Finally, it is possible now to assess how much uncertainties are shared between two models.

```python
c_k11  = sauna.Analysis.compare(sensitivities_1, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'c_k', covariances = covariances)
c_k12  = sauna.Analysis.compare(sensitivities_1, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'c_k', covariances = covariances)
c_k21  = sauna.Analysis.compare(sensitivities_2, sensitivities_1, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'c_k', covariances = covariances)
c_k22  = sauna.Analysis.compare(sensitivities_2, sensitivities_2, functional_1 = 'Eigenvalue', functional_2 = 'Eigenvalue', type = 'c_k', covariances = covariances)

```
