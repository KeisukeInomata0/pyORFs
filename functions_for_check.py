## module for functions 

# Import Python packages
import numpy as np
import math

####### functions for cross-checks: exact ORFs for PTA redshift response due to tensor and vector GWs (see Appendix E of arXiv:1406.4664 (Gair, Romano, Taylor, Mingarelli) for tensor GWs and Appendix J of arXiv:1506.08668  (Gair, Romano, Taylor) for vector GWs)

def facto(n): # shorthand function for factorial
    return math.factorial(n)

def Nlm(l,m): # N^m_l
    return np.sqrt((2*l+1)/(4*np.pi)*facto(l-m)/facto(l+m))

def Fm_qrlm(q,r,L,m,theta): # F^-_{q,r,L,m}, from Eqs.(E19) and (E23) in arXiv:1406.4664
    F = 0
    if r==0:
        for i in range(0,q+1):
            for j in range(m,L+1):
                F = F + pow(2,i-j)*pow(-1,q-i+j+m)*facto(q)*facto(L+j)*(pow(2,q-i+j-m+1) - pow(1+np.cos(theta),q-i+j-m+1))/(facto(i)*facto(q-i)*facto(j)*facto(L-j)*facto(j-m)*(q-i+j-m+1))
    elif r==1:
        for i in range(0,q):
            for j in range(m,L+1):
                F = F + pow(2,i-j)*pow(-1,q-i+j+m)*facto(q)*facto(L+j)*(pow(2,q-i+j-m) - pow(1+np.cos(theta),q-i+j-m))/(facto(i)*facto(q-i)*facto(j)*facto(L-j)*facto(j-m)*(q-i+j-m))
        for j in range(m+1,L+1):
            F = F + pow(2,q-j)*pow(-1,j+m)*facto(L+j)*(pow(2,j-m)-pow(1+np.cos(theta),j-m))/(facto(j)*facto(L-j)*facto(j-m)*(j-m))
        F = F + pow(2,q-m)*facto(L+m)/(facto(m)*facto(L-m))*np.log(2/(1+np.cos(theta)))
    elif r==-1:
        for i in range(0,q+1):
            for j in range(m,L+1):
                F = F + pow(2,i-j)*pow(-1,q-i+j+m)*facto(q)*facto(L+j)*(pow(2,q-i+j-m+2) - pow(1+np.cos(theta),q-i+j-m+2))/(facto(i)*facto(q-i)*facto(j)*facto(L-j)*facto(j-m)*(q-i+j-m+2))
    else: 
        print('error! r is neither 0 nor 1.')
    return F

def Fp_qrlm(q,r,L,m,theta): # F^-_{q,r,L,m}, from Eqs.(E21) in arXiv:1406.4664
    F = 0
    if r==0:
        for i in range(0,q+1):
            for j in range(m,L+1):
                F = F + pow(2,i-j)*pow(-1,L+q-i+j)*facto(q)*facto(L+j)*(pow(2,q-i+j-m+1) - pow(1-np.cos(theta),q-i+j-m+1))/(facto(i)*facto(q-i)*facto(j)*facto(L-j)*facto(j-m)*(q-i+j-m+1))
    elif r==1:
        for i in range(0,q):
            for j in range(m,L+1):
                F = F + pow(2,i-j)*pow(-1,L+q-i+j)*facto(q)*facto(L+j)*(pow(2,q-i+j-m) - pow(1-np.cos(theta),q-i+j-m))/(facto(i)*facto(q-i)*facto(j)*facto(L-j)*facto(j-m)*(q-i+j-m))
        for j in range(m+1,L+1):
            F = F + pow(2,q-j)*pow(-1,L+j)*facto(L+j)*(pow(2,j-m)-pow(1-np.cos(theta),j-m))/(facto(j)*facto(L-j)*facto(j-m)*(j-m))
        F = F + pow(-1,L+m)*pow(2,q-m)*facto(L+m)/(facto(m)*facto(L-m))*np.log(2/(1-np.cos(theta)))
    else: 
        print('error! r is neither 0 nor 1.')
    return F

def Gamma_I_zz_exact_posi(l,m,theta): # Gamma^{t,I}_{LM}(theta) for M>=0, from Eq.(E12) in arXiv:1406.4664
    yy=0
    if m==0:
        yy= -(1+np.cos(theta))*Fm_qrlm(0,0,l,0,theta) - (1-np.cos(theta))*Fp_qrlm(1,1,l,0,theta)
        if l==0:
            yy = yy + 1+np.cos(theta)/3
        elif l==1:
            yy = yy - (1+np.cos(theta))/3
        elif l==2:
            yy = yy + 2/15*np.cos(theta)
        yy = yy*1/2*np.sqrt((2*l+1)*np.pi)
    elif m==1:
        yy = -pow(1+np.cos(theta),3/2)/pow(1-np.cos(theta),1/2)*Fm_qrlm(1,0,l,1,theta) - pow(1-np.cos(theta),3/2)/pow(1+np.cos(theta),1/2)*Fp_qrlm(2,1,l,1,theta)
        if l==1:
            yy = yy + 2*np.sin(theta)/3
        elif l==2:
            yy = yy - 2*np.sin(theta)/5
        yy = yy*1/4*np.sqrt((2*l+1)*np.pi)*np.sqrt(facto(l-1)/facto(l+1))
    elif m>=2:
        yy = -1/4*np.sqrt((2*l+1)*np.pi)*np.sqrt(facto(l-m)/facto(l+m))* (
            pow(1+np.cos(theta),m/2+1)/pow(1-np.cos(theta),m/2)*Fm_qrlm(m,0,l,m,theta) 
            - pow(1+np.cos(theta),m/2)/pow(1-np.cos(theta),m/2-1)*Fm_qrlm(m-1,-1,l,m,theta)
            + pow(1-np.cos(theta),m/2+1)/pow(1+np.cos(theta),m/2)*Fp_qrlm(m+1,1,l,m,theta)
            - pow(1-np.cos(theta),m/2)/pow(1+np.cos(theta),m/2-1)*Fp_qrlm(m,0,l,m,theta)
            )
    return yy
    
def Gamma_I_zz_exact(l,m,theta):# Gamma^{t,I,zz}_{LM}(theta) for any M, from Eqs.(E12) and (E13) in arXiv:1406.4664
    yy = Gamma_I_zz_exact_posi(l,abs(m),theta)
    if m<0:
        yy = yy*pow(-1,m)
    return yy

def I_lm(l,m,theta): # I_lm, from Eq.(J4) in arXiv:1506.08668
    yy=0
    if m==1:
        yy = 2*pow(-1,l+1)
        if l==0:
            yy = yy + 2
        elif l==1:
            yy = yy - 4/3
        elif l==2:
            yy = yy + 4/5
        yy = yy*np.sin(theta)/2
    if m== -1:
        yy = 2*pow(-1,l+1)
        if l==0:
            yy = yy + 2
        elif l==1:
            yy = yy - 4/3
        elif l==2:
            yy = yy + 4/5
        yy = yy*(-Nlm(l,1)/Nlm(l,-1))*np.sin(theta)/2       
    return yy

def J_lm_posi(l,m,theta): # J_lm for m>=0, from Eqs.(J6) and (J7) in arXiv:1506.08668
    yy = 0
    if m==0:
        yy = Fm_qrlm(1,0,l,0,theta) + 2*Fp_qrlm(0,1,l,0,theta)-Fp_qrlm(1,0,l,0,theta)
        if l==0:
            yy = yy - 2/3*np.cos(theta) -2
        elif l==1:
            yy = yy + 2/3*np.cos(theta)
        elif l==2:
            yy = yy - 2/3*2/5*np.cos(theta)
    else:
        yy = pow(1+np.cos(theta),m/2)/pow(1-np.cos(theta),m/2)*(Fm_qrlm(m,0,l,m,theta) - Fm_qrlm(m-1,0,l,m,theta)) - pow(1-np.cos(theta),m/2)/pow(1+np.cos(theta),m/2)*(Fp_qrlm(m,0,l,m,theta) - Fp_qrlm(m,1,l,m,theta))
    return yy

def J_lm(l,m,theta): # J_lm for any m, from Eqs.(J6) and (J7) in arXiv:1506.08668
    yy = J_lm_posi(l,abs(m),theta)
    if m<0:
        yy = yy*pow(-1,-m)*Nlm(l,-m)/Nlm(l,m)
    return yy

def Gamma_I_zz_v_exact(l,m,theta):# Gamma^{v,I,zz}_{LM}(theta) for any M, from Eqs.(E12) and (E13) in arXiv:1406.4664
    return 2*np.pi*Nlm(l,m)*(I_lm(l,m,theta) + J_lm(l,m,theta))

