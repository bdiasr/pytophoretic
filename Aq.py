
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

import torch
from torchquad import Simpson, set_up_backend

from cte import *
from Fz import *



device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
set_up_backend("torch", data_type="float64")
num_gpus = torch.cuda.device_count()
torch.set_printoptions(precision=20)
dimension = 1
simp = Simpson()
integration_domain = [[0, L]] * dimension


def Aq(q, L):

    frac = 1 / L 
    integration_domain = [[0, L]]
    
    integrand = lambda z: (F_z6(z) * 
                           torch.exp(torch.as_tensor(1j * z * (2 * torch.pi * q) / L)))
    
    real_integrand = lambda x: integrand(x).real
    complex_integrand = lambda x: integrand(x).imag
    
    real_integral = simp.integrate((real_integrand), dim=1, N=301, integration_domain=integration_domain).item()
    complex_integral = simp.integrate((complex_integrand), dim=1, N=301, integration_domain=integration_domain).item()
    
    if complex_integral == 0.0:
        integral_val = real_integral * frac
    
    else:
        integral_val = (real_integral + 1j*complex_integral) * frac
        
    return integral_val

def Aq_1(q, L):

    frac = 1 / L 
    integration_domain = [[0, L]]
    
    integrand = lambda z: (F_z5(z) * 
                           torch.exp(torch.as_tensor(1j * z * (2 * torch.pi * q) / L)))
    
    real_integrand = lambda x: integrand(x).real
    complex_integrand = lambda x: integrand(x).imag
    
    real_integral = simp.integrate((real_integrand), dim=1, N=301, integration_domain=integration_domain).item()
    complex_integral = simp.integrate((complex_integrand), dim=1, N=301, integration_domain=integration_domain).item()
    
    if complex_integral == 0.0:
        integral_val = real_integral * frac
    
    else:
        integral_val = (real_integral + 1j*complex_integral) * frac
        
    return integral_val

def Aq_2(q):

    frac = 1 / L 
    integration_domain = [[0, L]]
    
    integrand = lambda z: (Fz(z) * 
                           torch.exp(torch.as_tensor(1j * z * (2 * torch.pi * q) / L)))
    
    real_integrand = lambda x: integrand(x).real
    complex_integrand = lambda x: integrand(x).imag
    
    real_integral = simp.integrate((real_integrand), dim=1, N=301, integration_domain=integration_domain).item()
    complex_integral = simp.integrate((complex_integrand), dim=1, N=301, integration_domain=integration_domain).item()
    
    if complex_integral == 0.0:
        integral_val = real_integral * frac
    
    else:
        integral_val = (real_integral + 1j*complex_integral) * frac
        
    return integral_val
