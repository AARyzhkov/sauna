---
icon: wrench
---

# Technological Uncertainty

The propagation of technological (modeling) uncertainties in SAUNA is based upon the equivalence of a change in the total cross section $$\sigma_{t,i}$$ and concentration $$C_i$$ of the corresponding nuclide $$i\in mat$$. That is, the sensitivities of a functional $$R$$ to the parameters are the same: $$S(R,\sigma_{t,i})=S(R,C_{i})$$. In the same manner, the sensitivity to the density of material $$\rho_{mat}$$ is the sum of the changes:

$$
S(R,\rho_{mat}) = \sum_{i\in mat}S(R,\sigma_{t,i})
$$

This can be expanded onto other parameters such as atomic fraction $$a_i$$ and weight fraction $$w_i$$.  However, these parameters are required to be renormalized to ensure the sums of the corresponding values are equal to unity. The renormalization is conducted in a Perfetti-like approach:

$$
S'(R,a_i)=S(R,\sigma_{t,i})-\frac{a_i}{1-a_i}\sum_{j\neq i}S(R,\sigma_{t,j})
$$

$$
S'(R,w_i)=\frac{S(R,\sigma_{t,i})}{1-w_i}-\frac{w_i}{1-w_i}\sum_{j\neq i}\frac{S(R,\sigma_{t,j})}{1-w_j}
$$

This approach does not pose any approximations beyond first-order perturbation theory since it does not assume perturbing the other technological parameters such as temperature and geometry that require special techniques to address.
