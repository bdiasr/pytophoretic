import numpy as np
import cmath

from mpmath import mp, fp

# Definindo a precisão da aritmética de ponto flutuante
fp.prec = 30 # Número de bits significativos em ponto flutuante

# Definindo a precisão decimal
mp.dps = 30  # Número de dígitos decimais

lamb = (514)*10**(-9)
k = (2 * np.pi)/lamb
L = 0.35
Q = 0.999965*k
qs = np.arange(-23, 24)    
spot = 0.0000235133  

M = (1.57 - 0.038j)
M1 = (1.57 - 0.38j)
M2 = 1.57 - 0.01*1j
M3 = 1.57 - 1j

M2_laop = (1.57 - 0.19j)
M3_laop = (1.57 - 0.95j)

mi = 1
eta_r = 1/M

# artigo de lentes - valores usados experimentalmente 
lamb1_len = (450)*10**(-9)
lamb2_len = (700)*10**(-9)
lamb3_len = (950)*10**(-9)
M1_lentes = 1.64 - 0.002205j # para lamb 1 
M2_lentes = 1.588 - 0.003453j
M3_lentes = 1.564 - 0.004224j

f1 = 150*10**(-3)
f2 = 50*10**(-3)
f3 = 150*10**(-3)
f4 = 25*10**(-3)

def beta_q(q):
    return Q + 2*np.pi*q/L

def h_q(q):
  return cmath.sqrt(2)*k*cmath.sqrt(1 - beta_q(q)/k)

def k_rho(k, alpha):
    return (k*np.sin(alpha))

def k_z(k, alpha):
    return (k*np.cos(alpha))

def Z0(k, alpha, x_0, y_0):
    rho_0 = (x_0**2 + y_0**2)**(1/2)
    return np.float64(rho_0 * k_rho(k, alpha))

def epsilonR2(m, ur):
    epsilon = (m**2) / ur
    return - epsilon.imag

def Eta_r(m, ur):
    epsilon = (m**2) / ur
    divUrEpsilon = ur/epsilon 
    return divUrEpsilon ** (1/2)  

def tam_x(a, k):
    return k*a

def tam_a(x, k):
    return (x/k)
