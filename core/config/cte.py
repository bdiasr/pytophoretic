import numpy as np
import cmath

from mpmath import mp, fp

# Definindo a precisão da aritmética de ponto flutuante
fp.prec = 30 # Número de bits significativos em ponto flutuante

# Definindo a precisão decimal
mp.dps = 30  # Número de dígitos decimais

points = 100

lamb = (514)*10**(-9)
k = (2 * np.pi)/lamb
L = 0.5
Q = 0.99997*k

# -----------------------------------------------------------------------------
# N = 29: truncamento da Frozen Wave (2N+1 = 59 feixes de Bessel)
# -----------------------------------------------------------------------------
# Limite FISICO, nao arbitrario. O artigo (Sec. 2) exige 0 <= beta_n <= w0/c = k
# para que so existam ondas propagantes em +z, sem evanescentes. Com
# beta_n = Q + 2 pi n / L:
#
#     n <= (k - Q) L / (2 pi) = (1.2224e7)(3e-5)(0.5)/(2 pi) = 29.18
#
# Portanto N = 29 e' o MAIOR truncamento admissivel. Acima disso beta_n > k,
# k_rho_n = sqrt(k^2 - beta_n^2) fica imaginario e J_v vira Bessel modificada
# I_v, que cresce com rho -- exatamente a onda evanescente que a condicao exclui.
#
# (Antes: N = 28, sem justificativa registrada. A mudanca desloca I_x/I_z em
#  ~0.2%. `fw_field.n_max_propagante()` recalcula este limite a partir de k,Q,L.)
# -----------------------------------------------------------------------------
N_FW = 29
qs = np.arange(-N_FW, N_FW + 1)
Qs = qs          # alias historico

spot = 0.0000235133

M =  1.64 - 0.002205j


mi = 1
eta_r = 1/M

# artigo de lentes - valores usados experimentalmente 
M1_lentes = 1.64 - 0.002205j # para lamb 1 


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
