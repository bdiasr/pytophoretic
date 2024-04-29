import numpy as np
import scipy.special as s
import cmath
from numba import njit

from cte import *
from coefficients import * 
from Aq import * 
from legendre import *

import warnings
warnings.filterwarnings(action='ignore', category=UserWarning)

#@njit
def Psi(rho, z):

    soma = []
    total = 0
    
    for q in qs:

        a = Aq(q, L)
        j0 = s.jv(0, (h_q(q) * rho))
        exponencial = np.exp(-1j * beta_q(q) * z)
        soma = a * j0 * exponencial
      
        total += soma

    return total

def Aq_val(L):
    aq = []
    for q in qs:
        aq.append(Aq(q, L))
    
    return aq, len(aq)
  
#aq_vals = Aq_val(L)[0]
#print(aq_vals, Aq_val(L)[1])



'''
# usar como referencia 
epslon1 = (M**2).real
epslon2 = -(M**2).imag

epslon22 = -(M2**2).imag
epslon23 = -(M3**2).imag
'''


