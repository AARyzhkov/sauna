---
icon: sensor
---

# Sensitivity Analysis

## Definition

Sensitivity analysis is intended to assess the influence of perturbing parameters on a functional (output parameter, system response, etc.) allowing one to determine the parameters which require more attention while analyzing a system. To quantify the influence, the sensitivities (coefficients) are used. The sensitivity $$S(R,\alpha)$$ of an arbitrary functional $$R$$ that depends on a parameter $$\alpha$$ is defined as (fractional/relative) derivative:

$$
S(R,\alpha)\equiv\frac{\alpha}{R}\frac{\textnormal{d}R}{\textnormal{d}\alpha}
$$

Using this value one can assess the change of the functional $$R$$ from the value of $$R_0$$ in a linear approximation. If $$I$$ parameters $$\alpha_i$$ are perturbed by a value of $$\delta\alpha_i$$, $$i\in I$$, it can be represented as first-order Taylor series.

$$
R=R_0+\sum_{i=1}^I{S_i \delta\alpha_i}
$$

where $$S_i\equiv S(R_i,\alpha_i)$$ is used for brevity.

## Sensitivity assessment

Sensitivities are quantified via different approaches that can be divided into two groups relied on the direct and perturbation/variational methods. The first method is usually not considered due because it requires at least one additional simulation per parameter to perturb. The perturbation-based approach generally uses first-order perturbation theory. For the eigenvalue $$k$$ the equation has the following form:

$$
S(k,\alpha)=\frac{\alpha}{k}\frac{\textnormal{d} k}{\textnormal{d} \alpha} = -\alpha \frac{\langle\psi^*,\left(\frac{\partial \hat{A}}{\partial \alpha} - \frac{1}{k}\frac{\partial \hat{F}}{\partial \alpha}\right)\rangle}{\langle \psi^*,\frac{1}{k}\hat{F}\psi \rangle} +\mathcal{O}(\delta\psi)
$$

where $$\langle f,g \rangle \equiv \int f(\vec{\xi})g(\vec{\xi}) \textnormal{d}\vec{\xi}$$ is the scalar product of arbitrary functions $$f(\vec{\xi})$$ and $$g(\vec{\xi})$$; $$\vec{\xi}$$ is the phase-space vector; $$\hat{A}$$ is the transport operator; $$k$$ is the eigenvalue (multiplication factor); $$\hat{F}$$ is the fission operator; $$\delta\psi$$ is the perturbation in $$\psi$$ due to the perturbation in $$\alpha$$.&#x20;

This formula clearly shows the advantage of the perturbation-based approach: it requires only the knowledge of the forward and adjoint fluxes, i.e. requires two neutron transport simulations to get the sensitivities to all the $$I$$ parameters while the direct approach requires at least $$I$$ perturbed simulations besides one unperturbed forward simulation, making this approach rather inefficient.

This is the classical approach for the eigenvalue $$k$$, and, to calculate the sensitivity of other functionals, one shall turn to Generalized perturbation theory also known as GPT introduced by Usachev. GPT has a similar advantage, but it might require two addition simulations if the functional is bi-linear (or their ratio) such as the effective delayed neutron fraction $$\beta_{eff}$$ or the other kinetic parameters. GPT has two limitations: 1) the functional $$R$$ has to satisfy the condition $$\langle \psi, \frac{\partial R}{\partial \psi} \rangle= 0$$ (e.g. $$^{135}\textnormal{Xe}$$ equilibrium concentration); 2) if the number of functionals is large enough (e.g. power distribution), the advantage of the perturbation-based approach diminishes.&#x20;

{% hint style="warning" %}
SAUNA does not compute sensitivities: it relies on the tools, having these techniques implemented such as Serpent and KENO-VI of SCALE. Still, they may be created manually, but the approach is rather inefficient.
{% endhint %}

The formulations were developed long ago though they are still widely used, and their development is rather concerned with how a stochastic neutron transport code can compute sensitivities in a single simulation. These techniques are implemented in a number of well-known tools such as MCNP, KENO-VI, Serpent, TRIPOLI-4, etc.
