---
icon: gear-complex
---

# Complex Sensitivities

Besides the eigenvalue sensitivities, there is a number of functionals of interest for analysis such as the reactivity coefficients $$\Delta\rho$$, effective delayed neutron fraction $$\beta_{eff}$$, effective neutron generation time $$\Lambda_{eff}$$, effective neutron lifetime $$\ell_{eff}$$, etc. However, not many stochastic neutron transport simulation tools are capable of computing these sensitivities directly, requiring further data curation.

## Reactivity difference&#x20;

Reactivity coefficients $$\Delta\rho$$ are a crucial part in a safety analysis. Neutron transport codes do not usually compute sensitivities to them, and obtaining them requires some additional efforts to get the uncertainty influence on them. Since this kind functionals are comprised of the eigenvalues in two states, the reactivity effect sensitivity can also be constructed from the corresponding two eigenvalue sensitivities.

$$
S(\Delta\rho,\alpha)=\frac{\alpha}{\Delta\rho}\frac{\text{d}(\Delta\rho)}{\text{d}\alpha}
$$

where $$\Delta\rho\equiv \frac{1}{k_0} - \frac{1}{k}$$ and $$k_0$$ is the reference state eigenvalue.

Expanding $$\Delta\rho$$ yields the following widely used expression:

$$
S(\Delta\rho,\alpha)=\frac{S(k,\alpha)/k-S(k_0,\alpha)/k_0}{\Delta\rho}=\frac{S(k,\alpha)k_0-S(k_0,\alpha)k}{k-k_0}
$$

## Effective delayed neutron fraction

In the case of the effective delayed neutron fraction $$\beta_{eff}$$, it can be approximated as:

$$
\beta_{eff}\approx1-\frac{k_p}{k}
$$

Through some algebraic actions one can obtain the following widely used equation for quantifying the sensitivity of the effective delayed neutron fraction via the results of two simulation:

$$
S(\beta_{eff},\alpha)=\frac{\alpha}{\beta_{eff}}\frac{\textnormal{d}\beta_{eff}}{\textnormal{d}\alpha}=\frac{\alpha}{\beta_{eff}}\frac{\textnormal{d}}{\textnormal{d}\alpha}\left(1-\frac{k_p}{k}\right)=\frac{k_p}{k-k_p}\left[S(k,\alpha)-S(k_p,\alpha)\right]
$$

where $$k_p$$ is the prompt eigenvalue, i.e. with no delayed neutrons.

## Prompt decay constant

Another functional is the prompt decay constant $$\alpha_p$$ (prompt alpha eigenvalue, Rossia-alpha), which is of interest from the point of validation and was firstly applied in the production of ENDF/B-VII.1. The sensitivity can be obtained via the following equation:

$$
S(\alpha_p,\alpha)=\frac{\alpha}{\alpha_p}\frac{\textnormal{d}\alpha_p}{\textnormal{d}\alpha}=\frac{\alpha}{\alpha_p}\frac{\textnormal{d}}{\textnormal{d}\alpha}\left(\frac{k_p-1}{\ell_{eff}}\right)=\frac{k_p}{k_p-1}S(k_p,\alpha)-S(\ell_{eff},\alpha)
$$

## Breeding ratio

The breeding ratio ($$BR$$) is another functional for consideration defining a ratio of fissile nuclide production rate, e.g. $$^{233}\textnormal{U}$$ and $$^{239}\textnormal{Pu}$$, and fissile nuclide removal rate, e.g. $$^{232}\textnormal{Th}$$ and $$^{238}\textnormal{U}$$. Modern stochastic tools do now allow one to compute the sensitivity to it directly and are limited to ratios of only one score in both numerator and denominator. This significantly limits the capability of computing sensitivities of more complex ratios such as the breeding ratio and, consequently, their analysis. To circumvent this limitation, it is suggested to compute the breeding ratio sensitivities through two different sensitivities:

$$
S(BR,\alpha)=\frac{\alpha}{BR}\frac{\textnormal{d}BR}{\textnormal{d}\alpha}=\frac{\alpha}{BR}\frac{\textnormal{d}BR}{\textnormal{d}\alpha}=-\frac{1}{1+R_f/R_{\gamma}}S(R_{\gamma},α)-\frac{1}{1+R_{\gamma}/R_f}S(R_f,α)
$$

where $$R_f$$ is the ratio of the $$(n,f)$$ reaction on a fissile nuclide and the $$(n,\gamma)$$ reaction on a fertile nuclide; $$R_{\gamma}$$ is the ratio of the $$(n,\gamma)$$ reaction on a fissile nuclide and the $$(n,\gamma)$$ reaction on a fertile nuclide.

It is possible to account any number of reactions using general formulation:

$$
S(BR,\alpha)=-\frac{1}{BR}\sum_i BR^2_i\left( \sum_j[R_{\gamma,j,i}S(R_{\gamma,j,i},\alpha)+R_{f,j,i}S(R_{f,j,i},\alpha)]\right)
$$

where $$BR_i$$ is the production rate on the $$i$$-nuclide, i.e. $$BR=\sum_i BR_i$$; $$R_{x,j,i}$$ is the ratio of the $$x$$ reaction on a fissile $$j$$-nuclide and the $$(n,\gamma)$$ reaction on the fertile nuclide.

Although the approach is general, it would require more quantifications of different ratios making it rather cumbersome. For example, accounting $$^{240}\textnormal{Pu}$$ and $$^{241}\textnormal{Pu}$$ would require calculating eight ratios. Thus, only previous formulation is considered reasonable to use in SAUNA.

## Ratio of functionals

There are some functionals that can be represented as a ratio of them. The basic idea for these functionals is that the sensitivity of them is merely a difference two sensitivities representing the numerator and denominator. For a critical system the prompt decay constant is $$\alpha_p=-\beta_{eff}/\Lambda_{eff}$$ and $$\Lambda_{eff}=\ell_{eff}$$, and the sensitivity equation becomes simpler:

$$
S(\alpha_p,\alpha)=S(\beta_{eff},\alpha)-S(\Lambda_{eff},\alpha)=S(\beta_{eff},\alpha)-S(\ell_{eff},\alpha)
$$

Another example is the simplified representation of the neutron generation time $$\Lambda_{eff}$$ if the effective neutron lifetime $$\ell_{eff}$$ is known (the case of using Serpent):
