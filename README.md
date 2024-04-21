# Notebooks for overlap reduction functions

Jupyter Notebooks associated with arXiv:240X.XXXXX by Keisuke Inomata, Marc Kamionkowski, Stephen Taylor, and Celia Toral.
These numerically calculate the overlap reduction functions (ORFs) for the redshift modification in pulsar timing array (PTA) and the deflections in astrometry.

- [auto_pta.ipynb](auto_pta.ipynb): auto-correlation of PTA redshifts. 
- [auto_astrometry.ipynb](auto_astrometry.ipynb): auto-correlation of astrometry deflections. 
- [cross.ipynb](cross.ipynb): cross-correlation between PTA redshifts and astrometry deflections. 

The codes calculate the ORF $\Gamma^{p,X,ST}\_{LM,\alpha\beta}(\hat n_a, \hat n_b)$ in the following coordinates:
$$
  \hat n_a = (0,0,1),\  \hat n_b = (\sin \theta, 0, \cos \theta).
$$
In this coordinate choice, $\Gamma^{p,X,ST}\_{LM,\alpha\beta}$ only depends on $\theta$, which are calculated by the following functions:
Gamma_X_ST$(L, M, \theta \text{\[rad\]}, l_\text{max}, alpha, beta) \text{for} p = t$  
Gamma_X_ST_v$(L, M, \theta \text{\[rad\]}, l_\text{max}, alpha, beta) \text{for} p = v$,  
where l_max is the maximum value of $\ell$ that we take into account (see the paper) and $\alpha$ and/or $\beta$ may not be there if $S$ and/or $T$ are $z$. $\alpha,\beta = 0$ corresponds to $\theta$ and $\alpha,\beta = 1$ to $\phi$.  
For example, Gamma_I_EB_v(L, M, theta\[rad\], l_max, 0, 1) $= \Gamma^{v,I,EB}\_{LM,\theta \phi}(\theta)$. 

## Author
- [Keisuke Inomata](mailto:kinomat1@jhu.edu) (Johns Hopkins University)

## Dependencies
- Python
- matplotlib
- numpy, math, scipy, sympy

