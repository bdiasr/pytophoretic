import os
import json 

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

path = r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/'

def cache_FW1_Psi():

    Z = np.linspace(0, L, 500)
    cnt = 0

    FW_Psi = {"FW_Psi":[]}
    
    print('Calculando Psi de FW1...\n')
    for zi in Z:
        print("i:", cnt)

        FW_Psi['FW_Psi'].append({
            "Psi": np.float64(abs(Psi_2(0, zi))**2),
            "z": zi
        })
        json.dump(FW_Psi, open(os.path.abspath(path + 'FW_Psi.json'), 'w'))
        cnt = cnt + 1

def cache_FW5_Psi():

    Z = np.linspace(0, L, 500)
    cnt = 0

    FW5_Psi = {"FW5_Psi":[]}
    
    print('Calculando Psi de FW5...\n')
    for zi in Z:
        print("i:", cnt)

        FW5_Psi['FW5_Psi'].append({
            "Psi": np.float64(abs(Psi_1(0, zi))**2),
            "z": zi
        })
        json.dump(FW5_Psi, open(os.path.abspath(path + 'FW5_Psi.json'), 'w'))
        cnt = cnt + 1


def cache_FW5_Jx():

    Z = np.linspace(0, L, 100)
    cnt = 0
    eps = epsilonR2(M1_lentes, 1)
    rho = np.linspace(-float(spot), float(spot), 150)

    jx_FW5_01_L3 = {"jx_FW5_01_L3":[]}
    jx_FW5_05_L3 = {"jx_FW5_05_L3":[]}

    jx_FW5_01_L2 = {"jx_FW5_01_L2":[]}
    jx_FW5_05_L2 = {"jx_FW5_05_L2":[]}

    jx_FW5_01_L5 = {"jx_FW5_01_L5":[]}
    jx_FW5_05_L5 = {"jx_FW5_05_L5":[]}

    print('Calculando Jx de FW5...\n')
    for rhoi in rho:
        print("i:", cnt)

        jx_FW5_01_L3['jx_FW5_01_L3'].append({
            "jx": np.float64(I_x(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.1)),
            "rho": rhoi
        })
        json.dump(jx_FW5_01_L3, open(os.path.abspath(path + 'jx_FW5_01_L3.json'), 'w'))

        jx_FW5_05_L3['jx_FW5_05_L3'].append({
            "jx": np.float64(I_x(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.1)),
            "rho": rhoi
        })
        json.dump(jx_FW5_05_L3, open(os.path.abspath(path + 'jx_FW5_05_L3.json'), 'w'))

        jx_FW5_01_L2['jx_FW5_01_L2'].append({
            "jx": np.float64(I_x(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.25)),
            "rho": rhoi
        })
        json.dump(jx_FW5_01_L2, open(os.path.abspath(path + 'jx_FW5_01_L2.json'), 'w'))
        jx_FW5_05_L2['jx_FW5_05_L2'].append({
            "jx": np.float64(I_x(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.25)),
            "rho": rhoi
        })
        json.dump(jx_FW5_05_L2, open(os.path.abspath(path + 'jx_FW5_05_L2.json'), 'w'))

        
        jx_FW5_01_L5['jx_FW5_01_L5'].append({
            "jx": np.float64(I_x(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.4)),
            "rho": rhoi
        })
        json.dump(jx_FW5_01_L5, open(os.path.abspath(path + 'jx_FW5_01_L5.json'), 'w'))
        jx_FW5_05_L5['jx_FW5_05_L5'].append({
            "jx": np.float64(I_x(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.4)),
            "rho": rhoi
        })
        json.dump(jx_FW5_05_L5, open(os.path.abspath(path + 'jx_FW5_05_L5.json'), 'w'))
        

        cnt = cnt + 1

def cache_FW5_Jz():

    eps = epsilonR2(M1_lentes, 1)
    zval = np.linspace(0, L, 100)
    cnt = 0

    jz_FW5_01 = {"jz_FW5_01":[]}
    jz_FW5_1 = {"jz_FW5_1":[]}
    jz_FW5_3 = {"jz_FW5_3":[]}

    #adicionar verificação se o arquivo ja existe 
    print('Calculando Jz de FW5...\n')
    for zi in zval:
        print("i:", cnt)

        jz_FW5_01['jz_FW5_01'].append({
            "jz": np.float64(I_z(M1_lentes, 0.1, eps, 0, 0, 0, zi)),
            "z": zi
        })
        json.dump(jz_FW5_01, open(os.path.abspath(path + 'jz_FW5_01.json'), 'w'))

        jz_FW5_1['jz_FW5_1'].append({
            "jz": np.float64(I_z(M1_lentes, 1, eps, 0, 0, 0, zi)),
            "z": zi
        })
        json.dump(jz_FW5_1, open(os.path.abspath(path + 'jz_FW5_1.json'), 'w'))

        jz_FW5_3['jz_FW5_3'].append({
            "jz": np.float64(I_z(M1_lentes, 3, eps, 0, 0, 0, zi)),
            "z": zi
        })
        json.dump(jz_FW5_3, open(os.path.abspath(path + 'jz_FW5_3.json'), 'w'))


        cnt = cnt + 1


def cache_FW5_Iy():
    """
    Plot de Iy para:
    z = 0.1
    z = 0.25
    z = 0.4
    """

    cnt = 0
    eps = epsilonR2(M1_lentes, 1)
    rho = np.linspace(-3*spot, 3*spot, 140)

    jy_FW5_01_L3 = {"jy_FW5_01_L3":[]}
    jy_FW5_05_L3 = {"jy_FW5_05_L3":[]}
    jy_FW5_1_L3 = {"jy_FW5_1_L3":[]}

    jy_FW5_01_L2 = {"jy_FW5_01_L2":[]}
    jy_FW5_05_L2 = {"jy_FW5_05_L2":[]}
    jy_FW5_1_L2 = {"jy_FW5_1_L2":[]}

    jy_FW5_01_L5 = {"jy_FW5_01_L5":[]}
    jy_FW5_05_L5 = {"jy_FW5_05_L5":[]}
    jy_FW5_1_L5 = {"jy_FW5_1_L5":[]}
    
    print(L)
    print("Iy para x = 0.1 , 0.5 e 1 , configuracao off axis, x-displacement (-3spot a 3 spot), z=0.125")
    for rhoi in rho:
        print("i:", cnt)
        
        #primeiro pico 1
        jy_FW5_01_L3['jy_FW5_01_L3'].append({
            "jy": np.float64(I_y(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.1)),
            "rho": rhoi
        })
        json.dump(jy_FW5_01_L3, open(os.path.abspath(path + 'jy_FW5_01_L3.json'), 'w'))

        jy_FW5_05_L3['jy_FW5_05_L3'].append({
            "jy": np.float64(I_y(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.1)),
            "rho": rhoi
        })
        json.dump(jy_FW5_05_L3, open(os.path.abspath(path + 'jy_FW5_05_L3.json'), 'w'))

        jy_FW5_1_L3['jy_FW5_1_L3'].append({
            "jy": np.float64(I_y(M1_lentes, 1, eps, 0, 0, rhoi, 0.1)),
            "rho": rhoi
        })
        json.dump(jy_FW5_1_L3, open(os.path.abspath(path + 'jy_FW5_1_L3.json'), 'w'))
        
        #vale 1 
        jy_FW5_01_L2['jy_FW5_01_L2'].append({
            "jy": np.float64(I_y(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.25)),
            "rho": rhoi
        })
        json.dump(jy_FW5_01_L2, open(os.path.abspath(path + 'jy_FW5_01_L2.json'), 'w'))

        jy_FW5_05_L2['jy_FW5_05_L2'].append({
            "jy": np.float64(I_y(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.25)),
            "rho": rhoi
        })
        json.dump(jy_FW5_05_L2, open(os.path.abspath(path + 'jy_FW5_05_L2.json'), 'w'))

        jy_FW5_1_L2['jy_FW5_1_L2'].append({
            "jy": np.float64(I_y(M1_lentes, 1, eps, 0, 0, rhoi, 0.25)),
            "rho": rhoi
        })
        json.dump(jy_FW5_1_L2, open(os.path.abspath(path + 'jy_FW5_1_L2.json'), 'w'))

        #segundo pico 0.1 0.5 e 1 
        jy_FW5_01_L5['jy_FW5_01_L5'].append({
            "jy": np.float64(I_y(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.4)),
            "rho": rhoi
        })
        json.dump(jy_FW5_01_L5, open(os.path.abspath(path + 'jy_FW5_01_L5.json'), 'w'))

        jy_FW5_05_L5['jy_FW5_05_L5'].append({
            "jy": np.float64(I_y(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.4)),
            "rho": rhoi
        })
        json.dump(jy_FW5_05_L5, open(os.path.abspath(path + 'jy_FW5_05_L5.json'), 'w'))

        jy_FW5_1_L5['jy_FW5_1_L5'].append({
            "jx": np.float64(I_y(M1_lentes, 1, eps, 0, 0, rhoi, 0.4)),
            "rho": rhoi
        })
        json.dump(jy_FW5_1_L5, open(os.path.abspath(path + 'jy_FW5_1_L5.json'), 'w'))

        cnt = cnt + 1



#cache_FW1_Psi()
#cache_FW5_Psi()
#cache_FW5_Jz()
#cache_FW5_Jx()