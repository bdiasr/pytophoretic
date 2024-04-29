

import scipy 
import numpy as np
import sympy as sym
import math


import scipy.integrate as integrate
import cmath
import scipy.special as s
import matplotlib.pyplot as plt
import os
import mpmath
from mpmath import fac, exp, sin, cos, besselj, legenp, sqrt


#spherical bessel function 1 type 
# j_n = spherical_jn(n, z, derivative = bool) = j_n
def js_n(n, x):
    return s.spherical_jn(n, x, False)

def js_nz(n, z):
    
    if(n<0):
        return 0
    else:
        return np.where(not isinstance(z, complex) and z<0, (-1)**n*s.spherical_jn(n, -z, False), s.spherical_jn(n, z, False))
        #return s.spherical_jn(n, z, False)
        
def j0_nz(n, z):
    return s.jv(n, z, out=None)

#spherical bessel function 2 type
#y_n = spherical_jn(n, z, derivative = bool)
def ys_n(n, x):
    return s.spherical_yn(n, x, False)

#Ricatti Bessel 1type given by = x*(j_n(x))
def psiBessel(n, x):
    a = x*(s.spherical_jn(n, x, False))
    return a

#Ricatti Bessel 1type derivative
derivativePsiBesseln = (lambda n, x: ((1 + n)*s.spherical_jn(n, x, False) - (x*s.spherical_jn((n + 1), x, False))))

#Spherical Hankel H2 function 
def sphericalHankel_n(n, x):
    js_n = s.spherical_jn(n, x, False)
    ys_n = s.spherical_yn(n, x, False)
    return js_n - (ys_n*1j)

#Ricatti Bessel 2 type 
def RiccatiBessel_2Type(n, x):
    return x*(sphericalHankel_n(n, x))

#spherical hankel H2 derivative
def derivativeSphericalHankel_n(n, x):
    return (1+n)*sphericalHankel_n(n, x) - x*sphericalHankel_n((n+1), x)


