---
icon: underline
---

# Uncertainties of uncertainties

In many cases, propagating uncertainties relies on neutron transport simulations via stochastic codes. This is currently the case for SAUNA since it uses the sensitivity results of KENO-VI and Serpent. Consequently, the sensitivity uncertainties introduce uncertainties into results of uncertainty propagations. To be aware of statistical accuracy of calculated uncertainty, one should also propagate the statistical sensitivity uncertainty into the uncertainty of the functional. Besides some sensitivities are calculated through other functionals, and this also shall be taken into account.

## Accounting uncertainties of complex functionals

A calculation of complex functionals requires at least two functionals, and this results into a different uncertainty of the complex functional. The simplest case, considered here, is the reactivity difference $$\Delta\rho$$:

$$
\text{d}(\Delta\rho)=-\frac{\text{d}k}{k^2}+\frac{\text{d}k_0}{k_0^2}
$$

$$
\Delta(\Delta\rho)=\sqrt{\left(\frac{\Delta k}{k^2}\right)^2 + \left(\frac{\Delta k_0}{k_0^2}\right)^2}
$$

In addition, the statistical uncertainties of the complex functional sensitivities shall be computed:

$$
\Delta S (\Delta \rho, \alpha)=|S(\Delta\rho,\alpha)| \sqrt{\left[\frac{\Delta (S(k,\alpha)/k-S(k_0,\alpha)/k_0)}{S(k,\alpha)/k-S(k_0,\alpha)/k_0}\right]^2 + \left[\frac{\Delta(\Delta\rho)}{\Delta\rho}\right]^2 }
$$

or in terms of the eigenvalues:

$$
\Delta S (\Delta \rho, \alpha)=|S(\Delta\rho,\alpha)| \sqrt{\left[\frac{\Delta (S(k,\alpha)k_0-S(k_0,\alpha)k)}{S(k,\alpha)k_0-S(k_0,\alpha)k}\right]^2 + \left[\frac{\Delta(k_0-k)}{k_0-k}\right]^2 }
$$

## Statistical uncertainties of uncertainties

Provided the uncertainties for the inputs, one is able to compute their influence on the functional uncertainties. The first-order approximation is:

$$
\sigma^2=SCS^T
$$

and one may differentiate that:

$$
2\sigma\text{d}(\sigma)=2SC\text{d}S^T
$$

$$
\text{d}\sigma=\frac{SC\text{d}S^T}{\sigma}=\frac{SC\text{d}S^T}{\sqrt{SCS^T}}
$$

This equation permits obtaining how a functional uncertainty $$\sigma$$ is impacted by the statistical uncertainties of the sensitivities:

$$
\Delta\sigma=\sqrt{\frac{S'C'\text{d}S'^T}{SCS^T}}
$$

where $$[S']_i \equiv [S]_i^2$$ and $$[C']_{i,j} \equiv [C]_{i,j}^2$$.

However, sensitivities are often not equal, and have to be calculated taking this into account:

$$
\Delta\sigma=\frac{1}{2}\sqrt{\frac{\text{d}S'_1C'S'^T_2 + S'_1C'\text{d}S_2'^T}{S_1CS_2^T}}
$$

One may notice that using this equation with $$S_1=S_2$$ provides a different result from previous equation by a factor of $$\sqrt{2}$$. The latter formula assumes that the uncertainties of $$S_1$$ and $$S_2$$ are independent, while the former one assumes that they are fully correlated, because the values are the same. The difference can be seen using a formula for a covariance of a function of two arguments when $$S_1=S_2$$:

$$
\text{var}(\sigma)=\text{var}(AS_1)+2\text{cov}(AS_1,BS_2)+\text{var}(BS_2)=2\text{var}(AS)+2\text{cov}(AS,AS)=4\text{var}(AS)
$$

where A and B are the influence coefficients between $$\sigma$$ and the corresponding sensitivities $$S_1$$ and $$S_2$$, respectively.

According to this, the second formula assumes that $$\text{cov}(AS_1,BS_2)=0$$.
