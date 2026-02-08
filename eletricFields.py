import os
import time

import numpy as np
import sympy as sym

import scipy as sp 
import scipy.integrate as integrate
import scipy.special as s

import cmath
import math

from mpmath import fac, exp, sin, cos, besselj, fabs, legenp

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

from legendre import *
from specials import * 
from cte import * 
from Aq import * 
from coefficients import *

import warnings
warnings.filterwarnings(action='ignore', category=UserWarning)


def E_r_TM(E_0, k, r, theta, phi):

    ns = np.linspace(1, 100, 100)
    soma = []
    total = 0 
    for n in ns:
        for m in range(-n, n):
            c_pw_n = cpw_n(k, n)
            g_mnTM = 0
            psi_2_n = 0 
            psi_n = 0 
            p_mn = legenp(n, fabs(m), cos(theta))
            soma = c_pw_n*g_mnTM*(psi_2_n + psi_n)*p_mn*(exp(np.imag(phi)))
            total += soma

    return k*E_0*total

            

