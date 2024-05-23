# pyORFs

pyORFs is the code for arXiv:240X.XXXXX by Keisuke Inomata, Marc Kamionkowski, Celia M. Toral, and Stephen R. Taylor.
We numerically calculate the overlap reduction functions (ORFs) for the redshift modification in pulsar timing arrays (PTAs) and the deflections in astrometry:

- [auto_pta.ipynb](auto_pta.ipynb): auto-correlation of the PTA redshifts, $\Gamma^{p,X,zz}\_{LM}$,
- [auto_astrometry.ipynb](auto_astrometry.ipynb): auto-correlation of the astrometry deflections, $\Gamma^{p,X,ST}\_{LM,\alpha\beta}$,
- [cross.ipynb](cross.ipynb): cross-correlation between the PTA redshift and the astrometry deflection, $\Gamma^{p,X,Sz}\_{LM,\alpha}$ or $\Gamma^{p,X,zS}\_{LM,\alpha}$,  

where $S,T \in \\{E,B\\}$ for the astrometry modes.

We calculate the ORF $\Gamma^{p,X,ST}\_{LM}(\hat n_a, \hat n_b)$ (see the paper for its definition) in the following coordinates:  

$$
  \hat n_a = (0,0,1),\  \hat n_b = (\sin \theta, 0, \cos \theta).
$$  

Specifically, the following functions calculate the ORFs with overall imaginary numbers neglected:  
  
**Gamma_X_VW($L, M, \theta \\,\text{\[rad\]}, \ell_\text{max}, \alpha, \beta$) for $p = t$,**  
**Gamma_X_VW_v($L, M, \theta \\,\text{\[rad\]}, \ell_\text{max}, \alpha, \beta$) for $p = v$**,  
  
where $V,W \in \\{z,E,B\\}$ and $\alpha$ and/or $\beta$ may not be there if $V$ and/or $W$ is $z$.
$\ell_\text{max}$ is the maximum value of $\ell$ that we consider (see the paper). $\alpha,\beta = 0$ corresponds to $\theta$ and $\alpha,\beta = 1$ to $\phi$. For example, $\Gamma^{v,I,EB}\_{LM,\theta \phi}(\theta)$ can be calculated with Gamma_I_EB_v($L, M, \theta \\,\text{\[rad\]}, \ell_\text{max}, 0, 1$) except for the overall imaginary numbers.


When we plot the results, we recover the imaginary numbers that we neglect in the functions, which change the overall factor of the y-axis in some figures.


We also implement the exact analytic results for the redshift response due to the $I$-mode GW anisotropies, obtained in [Gair et al. 2014](https://arxiv.org/abs/1406.4664) for spin-2 GWs and [Gair et al. 2015](https://arxiv.org/abs/1506.08668) for spin-1 GWs:

**Gamma_I_zz_exact($L, M, \theta \\,\text{\[rad\]}, \ell_\text{max}$) for $p = t$,**  
**Gamma_I_zz_v_exact($L, M, \theta \\,\text{\[rad\]}, \ell_\text{max}$) for $p = v$**.



## Author
- [Keisuke Inomata](mailto:kinomat1@jhu.edu) (Johns Hopkins University)

## Dependencies
- Python
- matplotlib
- numpy, math, scipy, sympy


