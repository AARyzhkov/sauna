---
icon: bullseye-arrow
---

# Target Accuracy Requirements

The uncertainties of functionals are often large and should be decreased. The allowed maximum uncertainty of a functional are defined through target values called target accuracy requirements (TARs). TARs are used to determine whether the uncertainties satisfy the target values or not. The widely accepted target for the eigenvalue $$k$$ is 0.3%, and it usually shows the largest difficulty to achieve among other functionals. When too large uncertainties are revealed, it is possible to pose an optimization problem, which is solved based upon sensitivity and uncertainty analysis results. The problem permits one to find the nuclides, reactions, and energy ranges (groups) to reduce the uncertainty to the level of the target accuracy requirements corresponding to minimal cost. The optimization problem is usually defined as:

$$
\begin{cases} \displaystyle\min_{\delta\alpha_i}\sum_{i}\frac{\lambda_i}{\delta\alpha_i^2} \\ S_jCS^T_j\leq\delta R_{j,TAR}^2 \\ (\delta \alpha_i)_{\min}\leq\delta\alpha_i\leq(\delta\alpha_i)_0 \end{cases}
$$

where $$\delta\alpha_i$$ is the uncertainty of parameter $$\alpha$$; $$\lambda_i$$ is the cost parameter; $$S_j$$ is the sensitivity vector to $$j$$-functional $$R_j$$; $$C$$ is the covariance matrix; $$\delta R_{j,TAR}$$ is the target uncertainty for functional $$R_j$$; $$(\delta\alpha_i)_{\min}$$ is the minimal uncertainty; $$(\delta\alpha_i)_0$$ is the base uncertainty.

The first term in the problem is the objective function, defining the cost function. This function makes an assumption that the cost for each parameter $$\alpha_i$$ is proportional to its weight $$w_i$$, which in statistics corresponds to the inverse dispersion:

$$
w_i\equiv\frac{1}{\delta\alpha_i^2}
$$

This shows that the cost function is the sum of the costs per unit of weight. However, each reaction and energy range require different efforts. Thus, a proportionality coefficient $$\lambda_i$$ is introduced. That is, the objective function is the sum of $$\lambda_i$$-weighted inverse dispersion.

The second expression is the non-linear inequality constrain defining that the result of the uncertainty propagation is lower than the corresponding TAR. The number of functionals is arbitrary, but a larger number affects the convergence speed.

The third expression is the linear inequality constraint for the uncertainties of the input parameters. The first term, $$(\delta\alpha_i)_{\min}$$ represents the minimum uncertainty, which may be set to zero in the simplest case. Besides, the achievable value is limited by a measurement instrument. The value is accepted at a level of 0.5%; still, the value may be too optimistic in some cases, and 1% or more are more reasonable values, especially if the resonance region is considered.\
It is important to note that tight conditions might make the achievement of TAR impossible that is easy to demonstrate via the simplest eigenvalue $$k$$ representation:

$$
k=\frac{\nu\Sigma_f}{\Sigma_a}
$$

where $$\nu$$ is the neutron multiplicity; $$\Sigma_f$$ is the fission cross section; $$\Sigma_a$$ is the absorption cross section.\
It is obvious that if the neutron multiplicity uncertainty is larger than 0.3%, an eigenvalue uncertainty of 0.3% cannot be achieved, because the sensitivity $$S(k,\nu)$$ is unity.



