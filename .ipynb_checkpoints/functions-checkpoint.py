## module for functions 

# Import Python packages
import scipy.special as sc
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_3j
import numpy as np
import math



####### scalar and vector spherical harmonics 

def Ysph(l,m,theta): # Y_{lm}(theta,0): spherical harmonics with phi = 0
    yy=0
    if(abs(m)<=l):
        yy = sc.sph_harm(m,l,0,theta)
    #NOTE: sc.sph_harm(m,l,phi,theta) = Y_{lm}(theta,phi), https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm.html
    return yy

def Ysph_v(l,m,theta,S,alpha): # Y^S_{lm,alpha}(theta,0): vector spherical harmonics with phi = 0 and with an imaginary number i omitted, S=0 for Y^E, S=1 for Y^B, alpha = 0 for Y_theta, alpha = 1 for Y_phi. The expressions are from Eqs.(D1) and (D2) in arXiv:1810.02369.
    yy = 0
    if S==0: # Y^E
        if alpha == 0: # Y_theta
            yy = -1/2/np.sqrt(l*(l+1))*( np.sqrt((l-m)*(l+m+1))*Ysph(l,m+1,theta) - np.sqrt((l+m)*(l-m+1))*Ysph(l,m-1,theta) )
        elif alpha == 1: # Y_phi/i, 
            yy = -1/np.sqrt(l*(l+1))*m/np.sin(theta)*Ysph(l,m,theta)
        else:
            print('error! alpha is neither 0 nor 1.')
    elif S==1: # Y^B 
        if alpha == 0: # Y_theta/i
            yy = 1/np.sqrt(l*(l+1))*m/np.sin(theta)*Ysph(l,m,theta)
        elif alpha == 1: # Y_phi
            yy = -1/2/np.sqrt(l*(l+1))*( np.sqrt((l-m)*(l+m+1))*Ysph(l,m+1,theta) - np.sqrt((l+m)*(l-m+1))*Ysph(l,m-1,theta) )
        else:
            print('error! alpha is neither 0 nor 1.')
    else: 
        print('error! S is neither 0 nor 1.')
    return yy



####### bipolar spherical harmonics

def YY(l,ld,L,M,theta): # {Y_l(\hat n_a) x Y_ld(\hat n_b)}_{LM} with theta = arccos(\hat n_a \cdot \hat n_b)
    yy=0
    for m in range(-l,l+1):
        md = M-m # used the fact <l,m,ld,md|L,M> is nonzero only when M = m + md
        if md >= -ld and md <= ld: # restric the sum ranges 
            yy = yy+ CG(l,m,ld,md,L,M).doit()*Ysph(l,m,0)*Ysph(ld,md,theta)
            #NOTE: CG(l1,m1,l2,m2,L,M) = <l1,m1,l2,m2|L,M>, https://docs.sympy.org/latest/modules/physics/quantum/cg.html#sympy.physics.quantum.cg.cg_simp
    return yy

def YS_Y(l,ld,L,M,theta,S,alpha): # {Y^S_{l,alpha} x Y_ld} in S={E,B} and alpha,beta={theta,phi}
    yy=0
    for m in range(-l,l+1):
        md = M-m
        if md >= -ld and md <= ld:
            yy = yy+ CG(l,m,ld,md,L,M).doit()*Ysph_v(l,m,0.001,S,alpha)*Ysph(ld,md,theta)
            # NOTE: 0.001 is substituted to avoid the coordinate singularity 
    return yy

def Y_YS(l,ld,L,M,theta,S,alpha): # Y_l x Y^S_{ld,alpha} in S={E,B} and alpha,beta={theta,phi}
    yy=0
    for m in range(-l,l+1):
        md = M-m
        if md >= -ld and md <= ld:
            yy = yy+ CG(l,m,ld,md,L,M).doit()*Ysph(l,m,0)*Ysph_v(ld,md,theta,S,alpha)
    return yy

def YS_YT(l,ld,L,M,theta,S,T,alpha,beta): # Y^S_{l,alpha} x Y^T_{ld,beta} in S,T={E,B} and alpha,beta={theta,phi}
    yy=0
    for m in range(-l,l+1):
        md = M-m
        if md >= -ld and md <= ld: 
            yy = yy+ CG(l,m,ld,md,L,M).doit()*Ysph_v(l,m,0.001,S,alpha)*Ysph_v(ld,md,theta,T,beta) 
            # NOTE: 0.001 is substituted to avoid the coordinate singularity 
    return yy



####### z, E, and B coefficients with i omitted

def zlt(l): # z^t_l
    return pow(-1,l)*np.sqrt(4*np.pi*(2*l + 1)*math.gamma(l-2+1)/math.gamma(l+2+1))
def zlv(l): # z^v_l
    zlv0 = pow(-1,l+1)*np.sqrt(4*np.pi*(2*l + 1))/np.sqrt(l*(l+1))
    if l==1:
        zlv0 = zlv0 - pow(-1,l+1)*np.sqrt(4*np.pi*(2*l + 1))*np.sqrt(2)/3
    return zlv0

def Elt(l): # E^t_l
    return 2/np.sqrt(l*(l+1))*zlt(l)
def Elv(l): # E^v_l
    elv0 = pow(-1,l+1)*np.sqrt(4*np.pi*(2*l+1)/l/(l+1))/np.sqrt(l*(l+1))
    if l==1: 
        elv0 = elv0 - pow(-1,l+1)*np.sqrt(4*np.pi*(2*l+1)/l/(l+1))*2*np.sqrt(2)/3
    return elv0

def Blt(l): # B^t_l/i, NOTE: i omitted. Strictly, B^t_l = (-i) E^t_l
    return -Elt(l)
def Blv(l): # B^v_l/i, NOTE: i omitted.
    return -zlv(l)/np.sqrt(l*(l+1))



####### F coefficients with i omitted

def xLll(l,ld,L): # X^L_{l ld} = 1 for l+ld+L = even, = 0 for l+ld+L = odd
    return (1+pow(-1,l+ld+L))/2


## for zz with spin-2 GWs
def FI_zz(l,ld,L,M):
    return 2*zlt(l)*zlt(ld)*wigner_3j(l,ld,L,-2,2,0)*xLll(l,ld,L) 
    #wigner_3j(l1,l2,l3,m1,m2,m3), https://docs.sympy.org/latest/modules/physics/wigner.html#sympy.physics.wigner.wigner_3j
def FV_zz(l,ld,L,M):
    return -2*zlt(l)*zlt(ld)*wigner_3j(l,ld,L,-2,2,0)*(1-xLll(l,ld,L))
def FE_zz(l,ld,L,M):
    return 2*zlt(l)*zlt(ld)*wigner_3j(l,ld,L,2,2,-4)*xLll(l,ld,L)
def FB_zz(l,ld,L,M):
    return 2*zlt(l)*zlt(ld)*wigner_3j(l,ld,L,2,2,-4)*(1-xLll(l,ld,L))

## for zz with spin-1 GWs
def FI_zz_v(l,ld,L,M):
    return -2*zlv(l)*zlv(ld)*wigner_3j(l,ld,L,-1,1,0)*xLll(l,ld,L)
def FV_zz_v(l,ld,L,M):
    return 2*zlv(l)*zlv(ld)*wigner_3j(l,ld,L,-1,1,0)*(1-xLll(l,ld,L))
def FE_zz_v(l,ld,L,M):
    return 2*zlv(l)*zlv(ld)*wigner_3j(l,ld,L,1,1,-2)*xLll(l,ld,L)
def FB_zz_v(l,ld,L,M):
    return 2*zlv(l)*zlv(ld)*wigner_3j(l,ld,L,1,1,-2)*(1-xLll(l,ld,L))


## for EB with spin-2 GWs, NOTE: the coefficients below are actually F^{p,X,EB}/i
def FI_eb(l,ld,L,M):
    return -2*Elt(l)*Blt(ld)*wigner_3j(l,ld,L,-2,2,0)*(1-xLll(l,ld,L))
def FV_eb(l,ld,L,M):
    return 2*Elt(l)*Blt(ld)*wigner_3j(l,ld,L,-2,2,0)*xLll(l,ld,L)
def FE_eb(l,ld,L,M):
    return -2*Elt(l)*Blt(ld)*wigner_3j(l,ld,L,2,2,-4)*(1-xLll(l,ld,L))
def FB_eb(l,ld,L,M):
    return -2*Elt(l)*Blt(ld)*wigner_3j(l,ld,L,2,2,-4)*xLll(l,ld,L)

## for EB with spin-1 GWs
def FI_eb_v(l,ld,L,M):
    return 2*Elv(l)*Blv(ld)*wigner_3j(l,ld,L,-1,1,0)*(1-xLll(l,ld,L))
def FV_eb_v(l,ld,L,M):
    return -2*Elv(l)*Blv(ld)*wigner_3j(l,ld,L,-1,1,0)*xLll(l,ld,L)
def FE_eb_v(l,ld,L,M):
    return -2*Elv(l)*Blv(ld)*wigner_3j(l,ld,L,1,1,-2)*(1-xLll(l,ld,L))
def FB_eb_v(l,ld,L,M):
    return -2*Elv(l)*Blv(ld)*wigner_3j(l,ld,L,1,1,-2)*xLll(l,ld,L)


## for BE with spin-2 GWs, NOTE: the coefficients below are actually F^{p,X,BE}/i
### Added additional (-1) in the end to take into account B^* = -B^* because B \propto i
def FI_be(l,ld,L,M):
    return -2*Blt(l)*Elt(ld)*wigner_3j(l,ld,L,-2,2,0)*(1-xLll(l,ld,L))*(-1)
def FV_be(l,ld,L,M):
    return 2*Blt(l)*Elt(ld)*wigner_3j(l,ld,L,-2,2,0)*xLll(l,ld,L)*(-1)
def FE_be(l,ld,L,M):
    return 2*Blt(l)*Elt(ld)*wigner_3j(l,ld,L,2,2,-4)*(1-xLll(l,ld,L))*(-1)
def FB_be(l,ld,L,M):
    return 2*Blt(l)*Elt(ld)*wigner_3j(l,ld,L,2,2,-4)*xLll(l,ld,L)*(-1)
    
## for BE with spin-1 GWs
def FI_be_v(l,ld,L,M):
    return 2*Blv(l)*Elv(ld)*wigner_3j(l,ld,L,-1,1,0)*(1-xLll(l,ld,L))*(-1)
def FV_be_v(l,ld,L,M):
    return -2*Blv(l)*Elv(ld)*wigner_3j(l,ld,L,-1,1,0)*xLll(l,ld,L)*(-1)
def FE_be_v(l,ld,L,M):
    return 2*Blv(l)*Elv(ld)*wigner_3j(l,ld,L,1,1,-2)*(1-xLll(l,ld,L))*(-1)
def FB_be_v(l,ld,L,M):
    return 2*Blv(l)*Elv(ld)*wigner_3j(l,ld,L,1,1,-2)*xLll(l,ld,L)*(-1)

## NOTE: above expressions for FX_be() are used for Gamma^{Sz} and Gamma^{zS}


## for EE with spin-2 GWs
def FI_ee(l,ld,L,M):
    return 4/np.sqrt(l*(l+1)*ld*(ld+1))*FI_zz(l,ld,L,M)
def FV_ee(l,ld,L,M):
    return 4/np.sqrt(l*(l+1)*ld*(ld+1))*FV_zz(l,ld,L,M)
def FE_ee(l,ld,L,M):
    return 4/np.sqrt(l*(l+1)*ld*(ld+1))*FE_zz(l,ld,L,M)    
def FB_ee(l,ld,L,M):
    return 4/np.sqrt(l*(l+1)*ld*(ld+1))*FB_zz(l,ld,L,M)

## for EE with spin-1 GWs
def FI_ee_v(l,ld,L,M):
    return Elv(l)*Elv(ld)/zlv(l)/zlv(ld)*FI_zz_v(l,ld,L,M)
def FV_ee_v(l,ld,L,M):
    return Elv(l)*Elv(ld)/zlv(l)/zlv(ld)*FV_zz_v(l,ld,L,M)
def FE_ee_v(l,ld,L,M):
    return Elv(l)*Elv(ld)/zlv(l)/zlv(ld)*FE_zz_v(l,ld,L,M)
def FB_ee_v(l,ld,L,M):
    return Elv(l)*Elv(ld)/zlv(l)/zlv(ld)*FB_zz_v(l,ld,L,M)


## for BB with spin-2 GWs
def FI_bb(l,ld,L,M):
    return 2*Blt(l)*Blt(ld)*wigner_3j(l,ld,L,-2,2,0)*xLll(l,ld,L)
def FV_bb(l,ld,L,M):
    return -2*Blt(l)*Blt(ld)*wigner_3j(l,ld,L,-2,2,0)*(1-xLll(l,ld,L))
def FE_bb(l,ld,L,M):
    return -2*Blt(l)*Blt(ld)*wigner_3j(l,ld,L,2,2,-4)*xLll(l,ld,L)
def FB_bb(l,ld,L,M):
    return -2*Blt(l)*Blt(ld)*wigner_3j(l,ld,L,2,2,-4)*(1-xLll(l,ld,L))

## for BB with spin-1 GWs
def FI_bb_v(l,ld,L,M):
    return -2*Blv(l)*Blv(ld)*wigner_3j(l,ld,L,-1,1,0)*xLll(l,ld,L)
def FV_bb_v(l,ld,L,M):
    return 2*Blv(l)*Blv(ld)*wigner_3j(l,ld,L,-1,1,0)*(1-xLll(l,ld,L))
def FE_bb_v(l,ld,L,M):
    return -2*Blv(l)*Blv(ld)*wigner_3j(l,ld,L,1,1,-2)*xLll(l,ld,L)
def FB_bb_v(l,ld,L,M):
    return -2*Blv(l)*Blv(ld)*wigner_3j(l,ld,L,1,1,-2)*(1-xLll(l,ld,L))


## for zE with spin-2 GWs
def FI_ze(l,ld,L,M):
    return zlt(l)/Elt(l)*FI_ee(l,ld,L,M)
def FV_ze(l,ld,L,M):
    return zlt(l)/Elt(l)*FV_ee(l,ld,L,M)
def FE_ze(l,ld,L,M):
    return zlt(l)/Elt(l)*FE_ee(l,ld,L,M)
def FB_ze(l,ld,L,M):
    return zlt(l)/Elt(l)*FB_ee(l,ld,L,M)

## for zE with spin-1 GWs
def FI_ze_v(l,ld,L,M):
    return zlv(l)/Elv(l)*FI_ee_v(l,ld,L,M)
def FV_ze_v(l,ld,L,M):
    return zlv(l)/Elv(l)*FV_ee_v(l,ld,L,M)
def FE_ze_v(l,ld,L,M):
    return zlv(l)/Elv(l)*FE_ee_v(l,ld,L,M)
def FB_ze_v(l,ld,L,M):
    return zlv(l)/Elv(l)*FB_ee_v(l,ld,L,M)


## for Ez with spin-2 GWs
def FI_ez(l,ld,L,M):
    return zlt(ld)/Elt(ld)*FI_ee(l,ld,L,M)
def FV_ez(l,ld,L,M):
    return zlt(ld)/Elt(ld)*FV_ee(l,ld,L,M)
def FE_ez(l,ld,L,M):
    return zlt(ld)/Elt(ld)*FE_ee(l,ld,L,M)
def FB_ez(l,ld,L,M):
    return zlt(ld)/Elt(ld)*FB_ee(l,ld,L,M)

## for Ez with spin-1 GWs
def FI_ez_v(l,ld,L,M):
    return zlv(ld)/Elv(ld)*FI_ee_v(l,ld,L,M)
def FV_ez_v(l,ld,L,M):
    return zlv(ld)/Elv(ld)*FV_ee_v(l,ld,L,M)
def FE_ez_v(l,ld,L,M):
    return zlv(ld)/Elv(ld)*FE_ee_v(l,ld,L,M)
def FB_ez_v(l,ld,L,M):
    return zlv(ld)/Elv(ld)*FB_ee_v(l,ld,L,M)


## for zB with spin-2 GWs
def FI_zb(l,ld,L,M):
    return zlt(l)/Elt(l)*FI_eb(l,ld,L,M)
def FV_zb(l,ld,L,M):
    return zlt(l)/Elt(l)*FV_eb(l,ld,L,M)
def FE_zb(l,ld,L,M):
    return zlt(l)/Elt(l)*FE_eb(l,ld,L,M)
def FB_zb(l,ld,L,M):
    return zlt(l)/Elt(l)*FB_eb(l,ld,L,M)

## for zB with spin-1 GWs
def FI_zb_v(l,ld,L,M):
    return zlv(l)/Elv(l)*FI_eb_v(l,ld,L,M)
def FV_zb_v(l,ld,L,M):
    return zlv(l)/Elv(l)*FV_eb_v(l,ld,L,M)
def FE_zb_v(l,ld,L,M):
    return zlv(l)/Elv(l)*FE_eb_v(l,ld,L,M)
def FB_zb_v(l,ld,L,M):
    return zlv(l)/Elv(l)*FB_eb_v(l,ld,L,M)


## for Bz with spin-2 GWs
def FI_bz(l,ld,L,M):
    return zlt(ld)/Elt(ld)*FI_be(l,ld,L,M)
def FV_bz(l,ld,L,M):
    return zlt(ld)/Elt(ld)*FV_be(l,ld,L,M)
def FE_bz(l,ld,L,M):
    return zlt(ld)/Elt(ld)*FE_be(l,ld,L,M)
def FB_bz(l,ld,L,M):
    return zlt(ld)/Elt(ld)*FB_be(l,ld,L,M)

## for Bz with spin-1 GWs
def FI_bz_v(l,ld,L,M):
    return zlv(ld)/Elv(ld)*FI_be_v(l,ld,L,M)
def FV_bz_v(l,ld,L,M):
    return zlv(ld)/Elv(ld)*FV_be_v(l,ld,L,M)
def FE_bz_v(l,ld,L,M):
    return zlv(ld)/Elv(ld)*FE_be_v(l,ld,L,M)
def FB_bz_v(l,ld,L,M):
    return zlv(ld)/Elv(ld)*FB_be_v(l,ld,L,M)

####### ORFs 

## for zz with spin-2 GWs
def Gamma_I_zz(L,M,theta,lmax):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld: # used the fact <l m ld md | L M > is nonzero only when abs(l-ld)<=L and L<= l+ld
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_zz(l,ld,L,M)*YY(l,ld,L,M,theta)
    return gam

def Gamma_V_zz(L,M,theta,lmax):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_zz(l,ld,L,M)*YY(l,ld,L,M,theta)
    return gam

def Gamma_E_zz(L,M,theta,lmax):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_zz(l,ld,L,M)*YY(l,ld,L,M,theta)
    return gam

def Gamma_B_zz(L,M,theta,lmax):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_zz(l,ld,L,M)*YY(l,ld,L,M,theta)
    return gam

## for zz with spin-1 GWs   
def Gamma_I_zz_v(L,M,theta,lmax):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_zz_v(l,ld,L,M)*YY(l,ld,L,M,theta)
    return gam

def Gamma_V_zz_v(L,M,theta,lmax):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_zz_v(l,ld,L,M)*YY(l,ld,L,M,theta)
    return gam

def Gamma_E_zz_v(L,M,theta,lmax):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_zz_v(l,ld,L,M)*YY(l,ld,L,M,theta)
    return gam

def Gamma_B_zz_v(L,M,theta,lmax):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_zz_v(l,ld,L,M)*YY(l,ld,L,M,theta)
    return gam



## for EB with spin-2 GWs
def Gamma_I_EB(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_eb(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,1,alpha,beta)
    return gam
    
def Gamma_V_EB(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_eb(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,1,alpha,beta)
    return gam  
    
def Gamma_E_EB(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_eb(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,1,alpha,beta)
    return gam  
    
def Gamma_B_EB(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_eb(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,1,alpha,beta)
    return gam    

## for EB with spin-1 GWs
def Gamma_I_EB_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_eb_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,1,alpha,beta)
    return gam
    
def Gamma_V_EB_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_eb_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,1,alpha,beta)
    return gam
    
def Gamma_E_EB_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_eb_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,1,alpha,beta)
    return gam
    
def Gamma_B_EB_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_eb_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,1,alpha,beta)
    return gam



## for BE with spin-2 GWs
def Gamma_I_BE(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_be(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,0,alpha,beta)
    return gam
    
def Gamma_V_BE(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_be(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,0,alpha,beta)
    return gam   
    
def Gamma_E_BE(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_be(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,0,alpha,beta)
    return gam  
    
def Gamma_B_BE(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_be(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,0,alpha,beta)
    return gam    

## for BE with spin-1 GWs
def Gamma_I_BE_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_be_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,0,alpha,beta)
    return gam

def Gamma_V_BE_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_be_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,0,alpha,beta)
    return gam
    
def Gamma_E_BE_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_be_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,0,alpha,beta)
    return gam
    
def Gamma_B_BE_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_be_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,0,alpha,beta)
    return gam



## for EE with spin-2 GWs
def Gamma_I_EE(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_ee(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,0,alpha,beta)
    return gam
    
def Gamma_V_EE(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_ee(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,0,alpha,beta)
    return gam   
    
def Gamma_E_EE(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_ee(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,0,alpha,beta)
    return gam  
    
def Gamma_B_EE(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_ee(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,0,alpha,beta)
    return gam    

## for EE with spin-1 GWs
def Gamma_I_EE_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_ee_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,0,alpha,beta)
    return gam
    
def Gamma_V_EE_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_ee_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,0,alpha,beta)
    return gam
    
def Gamma_E_EE_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_ee_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,0,alpha,beta)
    return gam
    
def Gamma_B_EE_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_ee_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,0,0,alpha,beta)
    return gam



## for BB with spin-2 GWs
def Gamma_I_BB(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_bb(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,1,alpha,beta)
    return gam
    
def Gamma_V_BB(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_bb(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,1,alpha,beta)
    return gam   
    
def Gamma_E_BB(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_bb(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,1,alpha,beta)
    return gam  
    
def Gamma_B_BB(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_bb(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,1,alpha,beta)
    return gam    

## for BB with spin-1 GWs
def Gamma_I_BB_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_bb_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,1,alpha,beta)
    return gam
    
def Gamma_V_BB_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_bb_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,1,alpha,beta)
    return gam
    
def Gamma_E_BB_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_bb_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,1,alpha,beta)
    return gam
    
def Gamma_B_BB_v(L,M,theta,lmax,alpha,beta):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_bb_v(l,ld,L,M)*YS_YT(l,ld,L,M,theta,1,1,alpha,beta)
    return gam



## for Ez with spin-2 GWs
def Gamma_I_Ez(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_ez(l,ld,L,M)*YS_Y(l,ld,L,M,theta,0,alpha)
    return gam
    
def Gamma_V_Ez(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_ez(l,ld,L,M)*YS_Y(l,ld,L,M,theta,0,alpha)
    return gam
    
def Gamma_E_Ez(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_ez(l,ld,L,M)*YS_Y(l,ld,L,M,theta,0,alpha)
    return gam
    
def Gamma_B_Ez(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_ez(l,ld,L,M)*YS_Y(l,ld,L,M,theta,0,alpha)
    return gam

    
## for Ez with spin-1 GWs
def Gamma_I_Ez_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_ez_v(l,ld,L,M)*YS_Y(l,ld,L,M,theta,0,alpha)
    return gam
    
def Gamma_V_Ez_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_ez_v(l,ld,L,M)*YS_Y(l,ld,L,M,theta,0,alpha)
    return gam
    
def Gamma_E_Ez_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_ez_v(l,ld,L,M)*YS_Y(l,ld,L,M,theta,0,alpha)
    return gam
    
def Gamma_B_Ez_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_ez_v(l,ld,L,M)*YS_Y(l,ld,L,M,theta,0,alpha)
    return gam
    


## for Bz with spin-2 GWs
def Gamma_I_Bz(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_bz(l,ld,L,M)*YS_Y(l,ld,L,M,theta,1,alpha)
    return gam
    
def Gamma_V_Bz(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_bz(l,ld,L,M)*YS_Y(l,ld,L,M,theta,1,alpha)
    return gam
    
def Gamma_E_Bz(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_bz(l,ld,L,M)*YS_Y(l,ld,L,M,theta,1,alpha)
    return gam
    
def Gamma_B_Bz(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_bz(l,ld,L,M)*YS_Y(l,ld,L,M,theta,1,alpha)
    return gam
    

## for Bz with spin-1 GWs
def Gamma_I_Bz_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_bz_v(l,ld,L,M)*YS_Y(l,ld,L,M,theta,1,alpha)
    return gam
    
def Gamma_V_Bz_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_bz_v(l,ld,L,M)*YS_Y(l,ld,L,M,theta,1,alpha)
    return gam
    
def Gamma_E_Bz_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_bz_v(l,ld,L,M)*YS_Y(l,ld,L,M,theta,1,alpha)
    return gam
    
def Gamma_B_Bz_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_bz_v(l,ld,L,M)*YS_Y(l,ld,L,M,theta,1,alpha)
    return gam
    


## for zE with spin-2 GWs
def Gamma_I_zE(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_ze(l,ld,L,M)*Y_YS(l,ld,L,M,theta,0,alpha)
    return gam
    
def Gamma_V_zE(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_ze(l,ld,L,M)*Y_YS(l,ld,L,M,theta,0,alpha)
    return gam
    
def Gamma_E_zE(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_ze(l,ld,L,M)*Y_YS(l,ld,L,M,theta,0,alpha)
    return gam

def Gamma_B_zE(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_ze(l,ld,L,M)*Y_YS(l,ld,L,M,theta,0,alpha)
    return gam
    

## for zE with spin-1 GWs
def Gamma_I_zE_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_ze_v(l,ld,L,M)*Y_YS(l,ld,L,M,theta,0,alpha)
    return gam

def Gamma_V_zE_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_ze_v(l,ld,L,M)*Y_YS(l,ld,L,M,theta,0,alpha)
    return gam

def Gamma_E_zE_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_ze_v(l,ld,L,M)*Y_YS(l,ld,L,M,theta,0,alpha)
    return gam

def Gamma_B_zE_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_ze_v(l,ld,L,M)*Y_YS(l,ld,L,M,theta,0,alpha)
    return gam



## for zB with spin-2 GWs
def Gamma_I_zB(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_zb(l,ld,L,M)*Y_YS(l,ld,L,M,theta,1,alpha)
    return gam

def Gamma_V_zB(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_zb(l,ld,L,M)*Y_YS(l,ld,L,M,theta,1,alpha)
    return gam

def Gamma_E_zB(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_zb(l,ld,L,M)*Y_YS(l,ld,L,M,theta,1,alpha)
    return gam

def Gamma_B_zB(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(2,lmax+1):
        for ld in range(2,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_zb(l,ld,L,M)*Y_YS(l,ld,L,M,theta,1,alpha)
    return gam


## for zB with spin-1 GWs
def Gamma_I_zB_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FI_zb_v(l,ld,L,M)*Y_YS(l,ld,L,M,theta,1,alpha)
    return gam
    
def Gamma_V_zB_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FV_zb_v(l,ld,L,M)*Y_YS(l,ld,L,M,theta,1,alpha)
    return gam

def Gamma_E_zB_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FE_zb_v(l,ld,L,M)*Y_YS(l,ld,L,M,theta,1,alpha)
    return gam

def Gamma_B_zB_v(L,M,theta,lmax,alpha):
    gam = 0
    for l in range(1,lmax+1):
        for ld in range(1,lmax+1):
            if abs(l-ld) <= L and L <= l+ld:
               gam = gam + pow(-1,L)*np.sqrt(np.pi)*FB_zb_v(l,ld,L,M)*Y_YS(l,ld,L,M,theta,1,alpha)
    return gam



