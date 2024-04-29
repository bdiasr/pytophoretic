import time

import numpy as np


import scipy as sp 
import scipy.integrate as integrate
import scipy.special as special



from mpmath import fac, exp, sin, cos, besselj, legenp, sqrt, conj, pi

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

from legendre import *
from specials import * 
from cte import * 
from Aq import * 
from beamShape import *


def cs_n(M, mi, n, x):
    dem = (M*mi*((RiccatiBessel_2Type(n, x)*derivativePsiBesseln(n, x)) - derivativeSphericalHankel_n(n, x)*psiBessel(n, x)))
    num = (mi*RiccatiBessel_2Type(n, x)*derivativePsiBesseln(n, x*M) - M*derivativeSphericalHankel_n(n, x)*psiBessel(n, x*M))
    return dem/num

#coefficient D_n
def ds_n(M, mi, n, x):
    dem = (M*M)*(RiccatiBessel_2Type(n, x)*derivativePsiBesseln(n, x) - derivativeSphericalHankel_n(n, x)*psiBessel(n, x))
    num = (M*RiccatiBessel_2Type(n, x)*derivativePsiBesseln(n, M*x)) - (mi*derivativeSphericalHankel_n(n, x)*psiBessel(n, M*x))
    return dem/num


def r_n(m, n, x):
    psiConj = np.conjugate(psiBessel(n, m*x))
    dem = (m*psiBessel(n+1, m*x)*psiConj)
    num = m*m
    return (dem.imag)/(num.imag)

def S_n(M, n, x):
    
    m = M
    a = ((1j)/(2*((m**2).imag)))
    b = (m*((abs(psiBessel(n, m*x))**2)))
    c = (np.conjugate(m)) * (abs(psiBessel(n+1, m*x)))**2
    d = (m +((((2*(n + 1)*((m*m).real)/m)))))*r_n(m, n, x)
    e = (2*n + 1)*(np.conjugate(m)*r_n(m, n+1, x))

    return -a*(x*(b+c) - d + e)

def A_mn(n, M, x, v,alpha1, rho_0, phi_0, x_0, y_0, z_0):
    
    soma1 = []
    soma2 = []
    soma3 = []
    soma4 = []
    total1 = 0
    total2= 0
    total3=0
    total4=0
    
    a = 1/((n + 1)**2)
    
    for m in range(0, (n - 1 + 1)):
        
        c1 = cs_n(M, mi, n, x)
        c1conj = conj(cs_n(M, mi, n+1, x))
        c1_gTM = g_mnTM(m+1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c1_gTM_conj = conj(g_mnTM(m, n+1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        c1_fact_num = fac(n + m + 1)
        c1_fact_dem = fac(n - m - 1)
        
        soma1 = c1*c1conj*c1_gTM*c1_gTM_conj*(c1_fact_num/c1_fact_dem)
        total1 += soma1
    
    for m in range(0, n + 1):
        
        c2_conj = conj(cs_n(M, mi, n, x))
        c2 = cs_n(M, mi, n+1, x)
        c2_gTM_conj = conj(g_mnTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        c2_gTM = g_mnTM(m+1, n+1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c2_fact_num = fac(n + m + 2)
        c2_fact_dem = fac(n - m)
        
        soma2 = c2_conj*c2*c2_gTM_conj*c2_gTM*(c2_fact_num/c2_fact_dem)
        total2 += soma2
        
    for m in range(-n-1, -1 + 1):
        c3 = cs_n(M, mi, n, x)
        c3_conj = conj(cs_n(M, mi, n+1, x))
        c3_gTM = g_mnTM(m+1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c3_gTM_conj = conj(g_mnTM(m, n+1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        c3_fact_num = fac(n + np.fabs(m) + 1)
        c3_fact_dem = fac(n - np.fabs(m) + 1)
        
        soma3 = c3*c3_conj*c3_gTM*c3_gTM_conj*(c3_fact_num/c3_fact_dem)
        total3 += soma3
        
    for m in range(-n, -1 +1 ):
        c4_conj = conj(cs_n(M, mi, n, x))
        c4 = cs_n(M, mi, n+1, x)
        c4_gTM_conj = conj(g_mnTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        c4_gTM =g_mnTM(m+1, n+1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c4_fact_num = fac(n + np.fabs(m))
        c4_fact_dem = fac(n - np.fabs(m))
        
        soma4 = c4_conj*c4*c4_gTM_conj*c4_gTM*(c4_fact_num/c4_fact_dem)
        total4 += soma4
    
    resul = a*(total1 + total2 - total3 - total4)
    
    return resul

def B_mn(n, M, x, v,alpha1, rho_0, phi_0, x_0, y_0, z_0):
    
    
    soma1 = []
    soma2 = []
    soma3 = []
    soma4 = []
    total1 = 0
    total2= 0
    total3= 0
    total4= 0
    
    eta_r = Eta_r(M, 1)
    const = (abs(eta_r)**2)/(n + 1)**2
    
    for m in range(0, (n - 1 + 1)):
        
        c1 = ds_n(M, mi, n, x)
        c1conj = conj(ds_n(M, mi, n+1, x))
        c1_gTE = g_mnTE(m+1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c1_gTE_conj = conj(g_mnTE(m, n+1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        c1_fact_num = fac((n + m + 1))
        c1_fact_dem = fac((n - m - 1))
        
        soma1 = c1*c1conj*c1_gTE*c1_gTE_conj*(c1_fact_num/c1_fact_dem)
        total1 += soma1
    
    for m in range(0, n+1):
        
        c2_conj = conj(ds_n(M, mi, n, x))
        c2 = ds_n(M, mi, n+1, x)
        c2_gTE_conj = conj(g_mnTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        c2_gTE =g_mnTE(m+1, n+1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c2_fact_num = fac(n + m + 2)
        c2_fact_dem = fac(n - m)
        
        soma2 = c2_conj*c2*c2_gTE_conj*c2_gTE*(c2_fact_num/c2_fact_dem)
        total2 += soma2
        
    for m in range(-n-1, -1+1):
        
        c3 = ds_n(M, mi, n, x)
        c3_conj = conj(ds_n(M, mi, n+1, x))
        c3_gTE = g_mnTE(m+1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c3_gTE_conj = conj(g_mnTE(m, n+1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        c3_fact_num = fac(n + abs(m) + 1)
        c3_fact_dem = fac(n - abs(m) + 1)
        
        soma3 = c3*c3_conj*c3_gTE*c3_gTE_conj*(c3_fact_num/c3_fact_dem)
        total3 += soma3
        
    for m in range(-n, -1 + 1):
        
        c4_conj = conj(ds_n(M, mi, n, x))
        c4 = ds_n(M, mi, n+1, x)
        c4_gTE_conj = conj(g_mnTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        c4_gTE = g_mnTE(m+1, n+1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c4_fact_num = fac(n + abs(m))
        c4_fact_dem = fac(n - abs(m))
        
        soma4 = c4_conj*c4*c4_gTE_conj*c4_gTE*(c4_fact_num/c4_fact_dem)
        total4 += soma4
    
    resul = const*(total1 + total2 - total3 - total4)
    
    return resul

def C_mnx(n, M, x, v,alpha1, rho_0, phi_0, x_0, y_0, z_0):
    
    soma1 = []
    soma2 = []
    total1 = 0
    total2 = 0
    
    eta_r = Eta_r(M, 1)
    frac = ((2*n) + 1)/((n**2)*((n+1)**2))
    conj_eta = conj(eta_r)
    cn = cs_n(M, mi, n, x)
    dn = conj(ds_n(M, mi, n, x))
    
    for m in range(0, (n)):
        
        c1_gTM = g_mnTM(m+1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c1_conj_gTE = conj(g_mnTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        c1_gn_TM2 =  g_mnTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c1_conj_gTE2 = conj(g_mnTE(m+1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        
        c1_fact_num = fac(n + m + 1)
        c1_fact_dem = fac(n - m - 1)
        
        soma1 = (c1_gTM*c1_conj_gTE + c1_gn_TM2*c1_conj_gTE2)*(c1_fact_num/c1_fact_dem )
        total1 += soma1
    
    for m in range(-n, -1 + 1):
        
        c2_gTM = g_mnTM(m+1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c2_conj_gTE = conj(g_mnTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        c2_gn_TM2 =  g_mnTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c2_conj_gTE2 = conj(g_mnTE(m+1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        
        c2_fact_num = fac(n + abs(m))
        c2_fact_dem = fac(n - abs(m))
        
        soma2 = (c2_gTM*c2_conj_gTE + c2_gn_TM2*c2_conj_gTE2)*(c2_fact_num/c2_fact_dem )
        total2 += soma2
        
    resul = frac*conj_eta*cn*dn*(total1 - total2)
    
    return resul

def C_mny(n, M, x, v,alpha1, rho_0, phi_0, x_0, y_0, z_0):
    
    soma1 = []
    soma2 = []
    total1 = 0
    total2= 0
    
    eta_r = Eta_r(M, 1)
    frac = (2*n + 1)/((n*(n+1))**2)
    conj_eta = conj(eta_r)
    cn = cs_n(M, mi, n, x)
    dn = conj(ds_n(M, mi, n, x))
    
    for m in range(0, (n - 1 + 1)):
        
        c1_gTM = g_mnTM(m+1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c1_conj_gTE = np.conjugate(g_mnTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        c1_gn_TM2 =  g_mnTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c1_conj_gTE2 = np.conjugate(g_mnTE(m+1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        
        c1_fact_num = special.factorial(n + m + 1)
        c1_fact_dem = special.factorial(n - m - 1)
        
        soma1 = (c1_gTM*c1_conj_gTE - c1_gn_TM2*c1_conj_gTE2)*(c1_fact_num/c1_fact_dem )
        total1 += soma1
    
    for m in range(-n, -1 + 1):
        
        c2_gTM = g_mnTM(m+1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c2_conj_gTE = np.conjugate(g_mnTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        c2_gn_TM2 =  g_mnTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        c2_conj_gTE2 = np.conjugate(g_mnTE(m+1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        
        c2_fact_num = special.factorial(n + abs(m))
        c2_fact_dem = special.factorial(n - abs(m))
        
        soma2 = (c2_gTM*c2_conj_gTE - c2_gn_TM2*c2_conj_gTE2)*(c2_fact_num/c2_fact_dem )
        total2 += soma2
        
    resul = frac*conj_eta*cn*dn*(total1 - total2)
    
    return resul

def D_mn(m, n, M, x, v,alpha, rho_0, phi_0, x_0, y_0, z_0):
    
    alpha1 = alpha
    a1 = 1/(n + 1)**2
    
    bcsn = cs_n(M, mi, n, x)
    bc = np.conjugate(bcsn)
    
    b = bc*cs_n(M, mi, n+1, x) 
    
    cgs = g_mnTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
    c = np.conjugate(cgs)
    
    d = g_mnTM(m, n + 1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
    
    e_num = special.factorial(n + np.abs(m) + 1, True)
    e_dem = special.factorial(n - np.abs(m), True)
    e =  e_num/e_dem

    return (a1*b*c*d*e)

def E_mn(m, n, M, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0):
    
    alpha1 = alpha
    eta_r = Eta_r(M, 1)

    a1 = (abs(eta_r)**2) / ((n + 1)**2)
    
    b = np.conjugate(ds_n(M, mi, n, x))*ds_n(M, mi, n+1, x)
    c = np.conjugate(g_mnTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
    d = g_mnTE(m, n+1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
    
    e = special.factorial(n + np.abs(m) + 1, True)/(special.factorial(n - np.abs(m), True))

    return (a1*b*c*d*e)


def F_mn(m, n, M, x, v,alpha, rho_0, phi_0, x_0, y_0, z_0):
    
    alpha1 = alpha
    eta_r = Eta_r(M, 1)
    a1 = (m*np.conjugate(eta_r)) * ((2*n + 1)/((n*(n + 1))**2))
    
    b = (cs_n(M, mi, n, x))*(np.conjugate(ds_n(M, mi, n, x)))
    
    c = g_mnTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
    d = np.conjugate(g_mnTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
    
    e = special.factorial(n + np.abs(m), True)/special.factorial(n - np.abs(m), True)

    return (a1 * b * c * d * e)

'''


print(D_mn(1, 2, M, x_tam, 0, ax, 0.1*spoti, 0, 0, 0, 0))
print(A_mn( 2, M, x_tam, 0, ax, 0.1*spoti, 0, 0, 0, 0))
'''
def teste_coeff_serial():

    ax = pi/360
    x_tam = 1
    spoti = ((405/100)/(k *sin(ax)))
    inicio = time.perf_counter()

    A_mn( 2, M, x_tam, 0, ax, 0.1*spoti, 0, 0, 0, 0)
    A_mn( 2, M, x_tam, 0, ax, 0.1*spoti, 0, 0, 0, 0)
    A_mn( 2, M, x_tam, 0, ax, 0.1*spoti, 0, 0, 0, 0)
    A_mn( 2, M, x_tam, 0, ax, 0.1*spoti, 0, 0, 0, 0)
    A_mn( 2, M, x_tam, 0, ax, 0.1*spoti, 0, 0, 0, 0)
    A_mn( 2, M, x_tam, 0, ax, 0.1*spoti, 0, 0, 0, 0)
    B_mn( 2, M, x_tam, 0, ax, 0.1*spoti, 0, 0, 0, 0)

    fim = time.perf_counter()

    total = fim - inicio
    print(f"tempo total serial: {total} segundos")

#teste_coeff_serial()