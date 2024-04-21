# Notebooks for overlap reduction functions

Jupyter Notebooks associated with arXiv:240X.XXXXX by Keisuke Inomata, Marc Kamionkowski, Stephen Taylor, and Celia Toral.
These numerically calculate the overlap reduction functions (ORFs) for the redshift modification in pulsar timing array (PTA) and the deflections in astrometry.

- [auto_pta.ipynb](auto_pta.ipynb): auto-correlation of PTA redshifts. 
- [auto_astrometry.ipynb](auto_astrometry.ipynb): auto-correlation of astrometry deflections. 
- [cross.ipynb](cross.ipynb): cross-correlation between PTA redshifts and astrometry deflections. 

The ORF $\Gamma^{p,X,ST}_{\alpha \beta}$ can be calculated by  
Gamma_X_ST(L, M, theta\[rad\], l_max, alpha, beta) for $p = t$  
Gamma_X_ST_v(L, M, theta\[rad\], l_max, alpha, beta) for $p = v$,  
where $\alpha$ and/or $\beta$ may not be there if $S$ and/or $T$ are $z$. 
$\alpha,\beta = 0$ corresponds to $\theta$ and $\alpha,\beta = 1$ to $\phi$.

## Author
- [Keisuke Inomata](mailto:kinomat1@jhu.edu) (Johns Hopkins University)

## Dependencies
- Python
- matplotlib
- numpy, math, scipy, sympy

