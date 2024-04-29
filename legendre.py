import numpy as np
from mpmath import fac, exp, sin, cos, besselj, legenp, sqrt
import cmath
import scipy.special as s
import time
from numpy import pi 

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True


from mpmath import mp, fp

# Definindo a precisão da aritmética de ponto flutuante
#fp.prec = 50 # Número de bits significativos em ponto flutuante
# Definindo a precisão decimal


def pi_n(m, n, x):

    if x == 1:
        fator =  cmath.sqrt(1 - 0.99999999999999)
    else:
        fator = cmath.sqrt(1 - x**2)
    
    #fator = cmath.sqrt(1 - x**2)
    frac =  1/fator
    aux = s.lpmn(m, n, x)[0]
    pi_val = aux[m][n]*frac
    
    return pi_val

def tau_n(m, n, x):

    fator = cmath.sqrt(1 - x**2)
    aux = s.lpmn(m, n, x)[1]
    tau_val = -(fator)*aux[m][n]
    return tau_val

def pi_mn(m, n, x):
    """
    sen(alpha) nao pode ser 0
    """
    
    num = legenp(n, m, cos(x))
    if(sin(x) == 0):
        return 0
    else:
        return num/sin(x)

def tau_mn(m, n, x):
    
    """
    Returns generalized Legendre function tau

    Derivative is calculated based on relation 14.10.5:
    http://dlmf.nist.gov/14.10
    funcionando corretamente

    x = cos de alpha 
    para frozen wave: x = beta_q/k
    x nao pode ser 0
    ou seja, cos (x) nao pode ser 1
    """

    pmn = legenp(n, m, x)
    pmn1 =  legenp(n+1, m, x)

    fristTerm = -(n + 1) * x * pmn
    secondTerm = (n - m + 1) * pmn1

    num = fristTerm + secondTerm
    
    if x == 1:
        den = sqrt(1 - 0.999999999)
    else:
        den = sqrt(1 - (x**2))
    
    return (num/den)

'''
print('Versão Mpmath')
start_time = time.time()
print(pi_mn(1, 1, cos(pi/360)))
print(tau_mn(1, 1, cos(pi/360)))
end_time = time.time()
print('Tempo de execução (Mpmath): {:.6f} segundos'.format(end_time - start_time))
'''
def fig_tau():

    n1 = []
    n2 = []
    n3 = []

    X = np.linspace(0, 1, 100)

    for x in X:
        n1.append(tau_mn(1, 1, x))

    for x in X:
        n2.append(tau_mn(1, 2, x))
  
    for x in X:
        n3.append(tau_mn(1, 3, x))

    plt.figure(figsize=[7,5])
    plt.plot(X, n1, 'b', label=r'$n_1$')
    plt.plot(X, n2, 'r', label=r'$n_3$')
    plt.plot(X, n3, 'g', label=r'$n_3$')

    plt.ylim(-6, 6)
    plt.xlabel(r'$\cos{\alpha}$')
    #plt.ylabel('')
    plt.title(r'$\tau^{m}_{n}(\cos{\alpha})$')
    plt.legend()
    plt.grid()
    plt.show()

def fig_pi():

    n1 = []
    n2 = []
    n3 = []

    X = np.linspace(0, 1, 100)

    for x in X:
        n1.append(pi_mn(1, 1, x))

    for x in X:
        n2.append(pi_mn(1, 2, x))
  
    for x in X:
        n3.append(pi_mn(1, 3, x))

    plt.plot(X, n1, 'b', label ='n1')
    plt.plot(X, n2, 'r', label='n2')
    plt.plot(X, n3, 'g', label='n3')
 
    plt.ylim(-6, 2)
    plt.xlabel('x')
    #plt.ylabel('')
    plt.legend()
    plt.grid()
    plt.show()

#fig_pi()
#fig_tau()  
