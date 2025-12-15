---
icon: gear-complex
---

# Complex Reactions

SAUNA allows one computing complex sensitivities via the $$\texttt{Sensitivities}$$ class and its methods.&#x20;

## Reactivity difference

Computing the reactivity difference sensitivity $$S(\Delta\rho,\alpha)$$ requires two $$\texttt{Sensitivities}$$ instances, two eigenvalues, corresponding to each instance, with their statistical uncertainties. Here is how to get a $$\texttt{Sensitivities}$$ instance for the reactivity difference. As the first six arguments have been mentioned, the last argument, $$\texttt{functional}$$ is provided to define the functional name to the new $$\texttt{Sensitivities}$$ instance.

```python
reactivity_sensitivities = sauna.Sensitivities.reactivity_difference(sensitivities_1, sensitivities_2, 1.016566, 1.023644, 0.000058, 0.000058, functional = 'Void reactivity')

```

This instance can be further applied in a conventional manner as other sensitivities.

