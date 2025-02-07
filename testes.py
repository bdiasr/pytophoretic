import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'

from cte import *
from coefficients import * 
from frozenWave import * 
from asymetryVector import *

import os.path
path = r'D:/Documentos/USP/pytophoretic/cache/'

'''
mp.dps = 80 # Número de dígitos decimais
print(mp)
eps = epsilonR2(M1_lentes, 1)
print(I_z(M1_lentes, 0.1, eps, 0, 0, 0, L/2))
#fp.prec = 30 # Número de bits significativos em ponto flutuante
'''

ceillingX = lambda x : math.ceil(x + (4.05 * (x ** (1/3))) + 2)
print(ceillingX(1))

