import os
import time

import numpy as np
import sympy as sym

import scipy as sp 
import scipy.integrate as integrate
import scipy.special as s

import cmath
import math

from mpmath import fac, exp, sin, cos, besselj, fabs

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

from legendre import *
from specials import * 
from cte import * 
from Aq import * 

import warnings
warnings.filterwarnings(action='ignore', category=UserWarning)

'''
def g_mnTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
    
    cosAxicon =  cos(alpha1)
    sinAxicon =  sin(alpha1)
    
    alpha = alpha1

    a1 = 1j*(-1/2)*(1j**(m+1))*(-1)**((m - abs(m)) * 1/2)   
    b = (fac(n - m))/(fac(n + abs(m)))
    c = exp(1j*k*cosAxicon*z_0)
    d = besselj(m-v-1, k*sinAxicon*rho_0)*exp(-1j*(m-v-1)*phi_0)
    e = (tau_mn(m,n, cosAxicon)/cosAxicon) + (m*pi_mn(m, n, alpha))
    f = besselj(m-v+1, k*sinAxicon*rho_0)*exp(-1j*(m-v+1)*phi_0)
    g = (-tau_mn(m, n, cosAxicon)/cosAxicon) + (m*pi_mn(m, n, alpha))
    #print(e)
    
    result = ((a1*b*c)*((d*e)-(f * g)))
    #print(a1, b, ((a1*b*c)*((d*e)-(f * g))))
    
    return result

def g_mnTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
    
    cosAxicon =  cos(alpha1)
    sinAxicon =  sin(alpha1)
    
    alpha = alpha1

    a1 = (1/2) * (1j**(m+1)) * (-1)**((m - abs(m)) * (1/2))
    b = (fac(n - m))/(fac(n + abs(m))) 
    c = exp(1j*k*cosAxicon*z_0) 
    #print(c)
    d = (besselj(m - v - 1,  k*sinAxicon*rho_0))*exp(-1j*(m - v - 1)*phi_0)
    e = (tau_mn(m, n, cosAxicon)/cosAxicon) + (m*pi_mn(m, n, alpha))
    f = (besselj(m - v + 1,  k*sinAxicon*rho_0))*exp(-1j*(m - v + 1)*phi_0) 
    g = (-tau_mn(m, n, cosAxicon)/cosAxicon) + (m*pi_mn(m, n, alpha)) 
    
    result =  (a1*b*c)*((d*e) + f*g)
    
    return result
'''
def g_mnTE_len(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
    
    soma = []
    total = 0

    a1 = 1j*(-1/2) * (1j**(m + 1 )) * (-1)**((m - abs(m))/2)  
    b = fac(n - m)/(fac(n + abs(m)))
    #qs = np.linspace(-1, 1, 3)
    print(b)
    for q in qs:

        c = (Aq(q, L)*exp(1j*beta_q(q)*z_0))
        d = s.jv(m - v - 1, h_q(q)*rho_0) * exp(-1j*(m-v-1)*phi_0) # ta sempre dando 1, ok
        e = (m*(pi_mn(m, n, beta_q(q)/k)))*(k/beta_q(q)) + (tau_mn(m, n, beta_q(q)/k))
        f = s.jv(m-v+1, h_q(q)*rho_0)*exp(-1j*(m-v+1)*phi_0)
        g = (m*pi_mn(m, n, beta_q(q)/k)) - (tau_mn(m, n, beta_q(q)/k)*(k/beta_q(q)))
        #print(c, d, e, f, g)
        #print(c)
        soma = c*((d*e) - (f*g))
        #print("q:", q, 'soma:', soma, '\n')
        total += soma

    #print(a1, b, a1*b*total)
    return a1*b*total

def g_mnTM_len(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
    
    soma = []
    total = 0
    
    a1 = (1/2) * (1j**(m+1)) * (-1)**((m - abs(m)) * (1/2))
    b = (fac(n - m))/(fac(n + abs(m))) 
    for q in qs:
        c = (Aq(q, L)*exp(1j*beta_q(q)*z_0))
        
        d = (besselj(m - v - 1,   h_q(q)*rho_0))*exp(-1j*(m - v - 1)*phi_0)
        e = (k/beta_q(q))*(tau_mn(m, n, beta_q(q)/k)) + (m*pi_mn(m, n, beta_q(q)/k))

        f = (besselj(m - v + 1,  h_q(q)*rho_0))*exp(-1j*(m - v + 1)*phi_0) 
        g = (k/beta_q(q))*(-tau_mn(m, n, beta_q(q)/k)) + (m*pi_mn(m, n, beta_q(q)/k)) 
        
        soma = c*((d*e) + (f*g))
        total+=soma
    
    return a1*b*total


def g_mnTE_FW(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):

    soma = []
    total = 0
    cosAxicon =  cos(alpha1)
    sinAxicon =  sin(alpha1)
    
    alpha = alpha1

    a1 = 1j*(-1/2)*(1j**(m+1))*(-1)**((m - abs(m)) * 1/2)   
    b = (fac(n - m))/(fac(n + abs(m)))
    print(b)

    for q in qs:

        c = (Aq(q, L)*cmath.exp(1j*beta_q(q)*z_0))
        d = besselj(m-v-1, k*sinAxicon*rho_0)*exp(-1j*(m-v-1)*phi_0)
        e = (tau_mn(m,n, cosAxicon)/cosAxicon) + (m*pi_mn(m, n, alpha))
        f = besselj(m-v+1, k*sinAxicon*rho_0)*exp(-1j*(m-v+1)*phi_0)
        g = (-tau_mn(m, n, cosAxicon)/cosAxicon) + (m*pi_mn(m, n, alpha))
        soma = c*((d*e) - (f*g))
        total += soma
    
    return a1*b*total

def g_mnTM_FW(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):

    soma = []
    total = 0

    cosAxicon =  cos(alpha1)
    sinAxicon =  sin(alpha1)
    
    alpha = alpha1

    a1 = (1/2) * (1j**(m+1)) * (-1)**((m - abs(m)) * (1/2))
    b = (fac(n - m))/(fac(n + abs(m))) 
    print(b)
    for q in qs:
        c = Aq(q, L)*exp(1j*k*cosAxicon*z_0) 
        d = (besselj(m - v - 1,  k*sinAxicon*rho_0))*exp(-1j*(m - v - 1)*phi_0)
        e = (tau_mn(m, n, cosAxicon)/cosAxicon) + (m*pi_mn(m, n, alpha))
        f = (besselj(m - v + 1,  k*sinAxicon*rho_0))*exp(-1j*(m - v + 1)*phi_0) 
        g = (-tau_mn(m, n, cosAxicon)/cosAxicon) + (m*pi_mn(m, n, alpha)) 
        
        soma = c*((d*e) - (f*g))
        total += soma

    return a1*b*total 

def g_mnTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
    soma = []
    total = 0
    
    a1 = (1/2) * (1j**(m+1)) * (-1)**((m - abs(m)) * (1/2))
    b = (fac(n - m))/(fac(n + abs(m))) 
    for q in qs:
        c = (Aq(q, L)*exp(1j*beta_q(q)*z_0))
        
        d = (besselj(m - v - 1,   h_q(q)*rho_0))*exp(-1j*(m - v - 1)*phi_0)
        e = (k/beta_q(q))*(tau_mn(m, n, beta_q(q)/k)) + (m*pi_mn(m, n, beta_q(q)/k))

        f = (besselj(m - v + 1,  h_q(q)*rho_0))*exp(-1j*(m - v + 1)*phi_0) 
        g = (k/beta_q(q))*(-tau_mn(m, n, beta_q(q)/k)) + (m*pi_mn(m, n, beta_q(q)/k)) 
        
        soma = c*((d*e) + (f*g))
        total+=soma
    
    return a1*b*total

def g_mnTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
    soma = []
    total = 0

    a1 = 1j*(-1/2) * (1j**(m + 1 )) * (-1)**((m - abs(m))/2)  
    b = fac(n - m)/(fac(n + abs(m)))
    for q in qs:

        c = (Aq(q, L)*exp(1j*beta_q(q)*z_0))
        d = s.jv(m - v - 1, h_q(q)*rho_0) * exp(-1j*(m-v-1)*phi_0) 
        e = (m*(pi_mn(m, n, beta_q(q)/k)))*(k/beta_q(q)) + (tau_mn(m, n, beta_q(q)/k))
        f = s.jv(m-v+1, h_q(q)*rho_0)*exp(-1j*(m-v+1)*phi_0)
        g = (m*pi_mn(m, n, beta_q(q)/k)) - (tau_mn(m, n, beta_q(q)/k)*(k/beta_q(q)))
        soma = c*((d*e) - (f*g))
        total += soma

    return a1*b*total

def gn_FrozenWave(n, N, k, L, Q, z0):

    soma = []
    total = 0

    primeiroTermoDem = n*(n + 1)
    qs = np.arange(N, (-N + 1))

    for q in qs:

        k_zq = Q + ((2*np.pi*q)/L)
        k_termo = k_zq/k
        primeiroTermoSum = Aq(q, L)/(1 + k_termo)
        primeiroTermoMul = (pi_mn(1, n, k_termo)) + (tau_mn(1, n, k_termo))
        exponencial = np.exp(1j * k_zq * z0)

        soma = primeiroTermoSum * primeiroTermoMul * exponencial
        total += soma
  
    gn = (-2/primeiroTermoDem)*(total)

    return gn

'''
axicon = np.pi/360
spoti = ((405/100)/(k * math.sin(axicon)))
'''

def gn_testes():
    '''
    VALORES PARA TESTAR DE AXICON
    0.1**(-10)
    np.pi/360
    np.pi/180
    np.pi/20
    np.pi/6
    np.pi/3
    '''
    '''
    print("Valores de gmn_TE para n = 1 \n")
    print("args = 0.1**(-10)")
    print("g_mnTE lentes: ", (g_mnTE_FW(1, 1, 0,  0.1**(-10), lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTE: ", (g_mnTE(1, 1, 0,  0.1**(-10), lamb, 0,0, 0, 0, 0)), "")
    print("\n args = np.pi/360")
    print("g_mnTE lentes: ", (g_mnTE_FW(1, 1, 0,  np.pi/360, lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTE: ", (g_mnTE(1, 1, 0,  np.pi/360, lamb, 0,0, 0, 0, 0)), "")
    print("\n args = np.pi/3")
    print("g_mnTE lentes: ", (g_mnTE_FW(1, 1, 0,  np.pi/3, lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTE: ", (g_mnTE(1, 1, 0,  np.pi/3, lamb, 0,0, 0, 0, 0)), "\n")

    print("Valores de gmn_TM para n = 1 \n")
    print("args = 0.1**(-10)")
    print("g_mnTM lentes: ", (g_mnTM_FW(1, 1, 0,  0.1**(-10), lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTM: ", (g_mnTM(1, 1, 0,  0.1**(-10), lamb, 0,0, 0, 0, 0)), "")
    print("\n args = np.pi/360")
    print("g_mnTM lentes: ", (g_mnTM_FW(1, 1, 0,  np.pi/360, lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTM: ", (g_mnTM(1, 1, 0,  np.pi/360, lamb, 0,0, 0, 0, 0)), "")
    print("\n args = np.pi/3")
    print("g_mnTM lentes: ", (g_mnTM_FW(1, 1, 0,  np.pi/3, lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTM: ", (g_mnTM(1, 1, 0,  np.pi/3, lamb, 0,0, 0, 0, 0)), "\n")
    '''
    
    print("Valores de gmn_TE para n = 2 \n")
    print("args = 0.1**(-10)")
    print("g_mnTE lentes: ", (g_mnTE_len(1, 2, 0,  0.1**(-10), lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTE: ", (g_mnTE(1, 2, 0,  0.1**(-10), lamb, 0,0, 0, 0, 0)), "")
    print("\n args = np.pi/360")
    print("g_mnTE lentes: ", (g_mnTE_len(1, 2, 0,  np.pi/360, lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTE: ", (g_mnTE(1, 2, 0,  np.pi/360, lamb, 0,0, 0, 0, 0)), "")
    print("\n args = np.pi/3")
    print("g_mnTE lentes: ", (g_mnTE_len(1, 2, 0,  np.pi/3, lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTE: ", (g_mnTE(1, 2, 0,  np.pi/3, lamb, 0,0, 0, 0, 0)), "")

    print("Valores de gmn_TM para n = 2 \n")
    print("args = 0.1**(-10)")
    print("g_mnTM lentes: ", (g_mnTM_len(1, 2, 0,  0.1**(-10), lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTM: ", (g_mnTM(1, 2, 0,  0.1**(-10), lamb, 0,0, 0, 0, 0)), "")
    print("\n args = np.pi/360")
    print("g_mnTM lentes: ", (g_mnTM_len(1, 2, 0,  np.pi/360, lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTM: ", (g_mnTM(1, 2, 0,  np.pi/360, lamb, 0,0, 0, 0, 0)), "")
    print("\n args = np.pi/3")
    print("g_mnTM lentes: ", (g_mnTM_len(1, 2, 0,  np.pi/3, lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTM: ", (g_mnTM(1, 2, 0,  np.pi/3, lamb, 0,0, 0, 0, 0)), "\n")

    
    print("Valores de gmn_TE para n = 3 \n")
    print("args = 0.1**(-10)")
    print("g_mnTE lentes: ", (g_mnTE_len(1, 3, 0,  0.1**(-10), lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTE: ", (g_mnTE(1, 3, 0,  0.1**(-10), lamb, 0,0, 0, 0, 0)), "")
    print("\n args = np.pi/360")
    print("g_mnTE lentes: ", (g_mnTE_len(1, 3, 0,  np.pi/360, lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTE: ", (g_mnTE(1, 3, 0,  np.pi/360, lamb, 0,0, 0, 0, 0)), "")
    print("\n args = np.pi/3")
    print("g_mnTE lentes: ", (g_mnTE_len(1, 3, 0,  np.pi/3, lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTE: ", (g_mnTE(1, 3, 0,  np.pi/3, lamb, 0,0, 0, 0, 0)), "")

    print("Valores de gmn_TM para n = 3 \n")
    print("args = 0.1**(-10)")
    print("g_mnTM lentes: ", (g_mnTM_len(1, 3, 0,  0.1**(-10), lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTM: ", (g_mnTM(1, 3, 0,  0.1**(-10), lamb, 0,0, 0, 0, 0)), "")
    print("\n args = np.pi/360")
    print("g_mnTM lentes: ", (g_mnTM_len(1, 3, 0,  np.pi/360, lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTM: ", (g_mnTM(1, 3, 0,  np.pi/360, lamb, 0,0, 0, 0, 0)), "")
    print("\n args = np.pi/3")
    print("g_mnTM lentes: ", (g_mnTM_len(1, 3, 0,  np.pi/3, lamb, 0,0, 0, 0, 0)), "")
    print("g_mnTM: ", (g_mnTM(1, 3, 0,  np.pi/3, lamb, 0,0, 0, 0, 0)), "\n")
    
    

#gn_testes()

#print("para n = 1, g_mnTE: ", abs(g_mnTE(1, 1, 0,  0.001, lamb, 0,0, 0, 0, 1)), "\n g_mnTE lentes: ", abs(g_mnTE_len(1, 1, 0,  0.001, lamb, 0,0, 0, 0, 1)))
#print("\n para n = 2, g_mnTE: ", abs(g_mnTE(1, 2, 0,  0.001, lamb, 0,0, 0, 0, 1)), "\n g_mnTE lentes:", abs(g_mnTE_len(1, 2, 0,  0.001, lamb, 0,0, 0, 0, 1)))
#print("\n para n = 3, g_mnTE: ", abs(g_mnTE(1, 3, 0,  0.001, lamb, 0,0, 0, 0, 1)), "\n g_mnTE lentes:", abs(g_mnTE_len(1, 3, 0,  0.001, lamb, 0,0, 0, 0, 1)))
#print("\n g_mnTE: ", abs(-0.5j* g_mnTE_len(1, 1, 0,  0.001, lamb, 0,0, 0, 0, 1)), "n = 1")
#print("\n g_mnTE: ", abs(-0.5j*g_mnTE_len(1, 2, 0,  0.001, lamb, 0,0, 0, 0,1)), "n = 2")
#print("\n g_mnTE: ", abs(-0.5j*g_mnTE_len(1, 3, 0,  0.001, lamb, 0,0, 0, 0, 1)), "n = 3")
#print("g_mnTM: ", abs(g_mnTM(1, 2, 0,  np.pi/360, lamb, 0,0, 0, 0, 1)))

#print("\n g_mnTE: ", abs(g_mnTE_len(1, 2, 0,  0.001, lamb, 0,0, 0, 0,1)), "n = 2")
#print("\n g_mnTE: ", abs(g_mnTE_len(1, 3, 0,  0.001, lamb, 0,0, 0, 0, 1)), "n = 3")

#print(" g_mnTE: ", abs(g_mnTE(1, 1, 0,  0.001, lamb, 0,0, 0, 0, 1)), "n = 1 \n ")
#print("\n g_mnTE: ", abs(g_mnTE(1, 2, 0,  0.001, lamb, 0,0, 0, 0,1)), "n = 2")
#print("\n g_mnTE: ", abs(g_mnTE(1, 3, 0,  0.001, lamb, 0,0, 0, 0, 1)), "n = 3")