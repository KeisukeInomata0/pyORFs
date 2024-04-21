# Notebooks for overlap reduction functions

Jupyter Notebooks associated with arXiv:240X.XXXXX by Keisuke Inomata, Marc Kamionkowski, Celia Toral, and Stephen Taylor.
These numerically calculate the overlap reduction functions (ORFs) for the redshift modification in pulsar timing array (PTA) and the deflections in astrometry.

- [auto_pta.ipynb](auto_pta.ipynb): auto-correlation of PTA redshifts. 
- [auto_astrometry.ipynb](auto_astrometry.ipynb): auto-correlation of astrometry deflections. 
- [cross.ipynb](cross.ipynb): cross-correlation between PTA redshifts and astrometry deflections. 

The codes calculate the ORF $\Gamma^{p,X,ST}\_{LM,\alpha\beta}(\hat n_a, \hat n_b)$ (see the paper for its definition) in the following coordinates:  

$$
  \hat n_a = (0,0,1),\  \hat n_b = (\sin \theta, 0, \cos \theta).
$$  

In this coordinate choice, $\Gamma^{p,X,ST}\_{LM,\alpha\beta}$ is calculated by the following Python functions:  
**Gamma_X_ST($L, M, \theta \\,\text{\[rad\]}, \ell_\text{max}, \alpha, \beta$) for $p = t$,**  
**Gamma_X_ST_v($L, M, \theta \\,\text{\[rad\]}, \ell_\text{max}, \alpha, \beta$) for $p = v$**,  
where $\ell_\text{max}$ is the maximum value of $\ell$ that we take into account (see the paper) and $\alpha$ and/or $\beta$ may not be there if $S$ and/or $T$ are $z$. $\alpha,\beta = 0$ corresponds to $\theta$ and $\alpha,\beta = 1$ to $\phi$. For example, $\Gamma^{v,I,EB}\_{LM,\theta \phi}(\theta)$ is calculated with Gamma_I_EB_v($L, M, \theta \\,\text{\[rad\]}, \ell_\text{max}, 0, 1$) except for the imaginary values. 


The Python functions neglect the imaginary numbers. So, when we plot the results, we recover the imaginary values that we neglect in the functions, which change the expression of the y-axis of the figures.  



## Author
- [Keisuke Inomata](mailto:kinomat1@jhu.edu) (Johns Hopkins University)

## Dependencies
- Python
- matplotlib
- numpy, math, scipy, sympy

