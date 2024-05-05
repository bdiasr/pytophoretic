aq_vals = [(-0.007630218676331162-1.0755285551056204e-16j), (0.012395268320895927-2.3332030751888055e-16j), (-0.009920225811417304-3.903127820947816e-17j), (0.0011471509682414492+2.6020852139652106e-17j), (0.009072280814612406+3.469446951953614e-18j), (-0.014597419014510766-1.9255430583342559e-16j), (0.011609993409216084-5.898059818321144e-17j), (-0.0011359339942944718+1.2576745200831851e-16j), (-0.011123305069012181+3.2959746043559335e-17j), (0.01778903804097489-2.949029909160572e-17j), (-0.014108430239921126+5.204170427930421e-17j), (0.0011268920658265916-5.724587470723463e-17j), (0.014271689008734465+4.85722573273506e-17j), (-0.02281807384565021-6.765421556309548e-17j), (0.018160381582904533-1.1102230246251565e-16j), (-0.001119942189787988+1.734723475976807e-17j), (-0.0197161404085515-1.214306433183765e-16j), (0.03188924077835213+3.2959746043559335e-17j), (-0.025835025583664163-2.7755575615628914e-17j), 0.0011150216489739188, (0.031395178845923286+2.0816681711721685e-17j), (-0.053086570366362114-3.122502256758253e-17j), (0.045823143709341554+5.551115123125783e-17j), (-0.0011120866000039614+1.3877787807814457e-17j), (-0.07425264036694881-2.0816681711721685e-17j), (0.1591665792157262+4.163336342344337e-17j), (-0.2258689101551024-2.42861286636753e-17j), 0.2511111111111111, (-0.2258689101551024+2.42861286636753e-17j), (0.1591665792157262-4.163336342344337e-17j), (-0.07425264036694881+2.0816681711721685e-17j), (-0.0011120866000039614-1.3877787807814457e-17j), (0.045823143709341554-5.551115123125783e-17j), (-0.053086570366362114+3.122502256758253e-17j), (0.031395178845923286-2.0816681711721685e-17j), 0.0011150216489739188, (-0.025835025583664163+2.7755575615628914e-17j), (0.03188924077835213-3.2959746043559335e-17j), (-0.0197161404085515+1.214306433183765e-16j), (-0.001119942189787988-1.734723475976807e-17j), (0.018160381582904533+1.1102230246251565e-16j), (-0.02281807384565021+6.765421556309548e-17j), (0.014271689008734465-4.85722573273506e-17j), (0.0011268920658265916+5.724587470723463e-17j), (-0.014108430239921126-5.204170427930421e-17j), (0.01778903804097489+2.949029909160572e-17j), (-0.011123305069012181-3.2959746043559335e-17j), (-0.0011359339942944718-1.2576745200831851e-16j), (0.011609993409216084+5.898059818321144e-17j), (-0.014597419014510766+1.9255430583342559e-16j), (0.009072280814612406-3.469446951953614e-18j), (0.0011471509682414492-2.6020852139652106e-17j), (-0.009920225811417304+3.903127820947816e-17j), (0.012395268320895927+2.3332030751888055e-16j), (-0.007630218676331162+1.0755285551056204e-16j)]


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
path = r'D:/Documentos/USP/pytophoretic/cache/'

def cache_FW6_Psi():

    Z = np.linspace(0, L, 500)
    cnt = 0

    FW_Psi = {"FW_Psi":[]}
    
    for zi in Z:
        print("i:", cnt)

        FW_Psi['FW_Psi'].append({
            "Psi": np.float64(abs(Psi(0, zi))**2),
            "z": zi
        })
        json.dump(FW_Psi, open(os.path.abspath(path + 'FW_Psi.json'), 'w'))
        cnt = cnt + 1

def cache_FW4_Psi():

    Z = np.linspace(0, L, 500)
    cnt = 0

    FW_Psi = {"FW_Psi":[]}
    
    for zi in Z:
        print("i:", cnt)

        FW_Psi['FW_Psi'].append({
            "Psi": np.float64(abs(Psi(0, zi))**2),
            "z": zi
        })
        json.dump(FW_Psi, open(os.path.abspath(path + 'FW4_Psi.json'), 'w'))
        cnt = cnt + 1

def cache_FW3_Psi():

    Z = np.linspace(0, L, 500)
    cnt = 0

    FW_Psi = {"FW_Psi":[]}
    
    for zi in Z:
        print("i:", cnt)

        FW_Psi['FW_Psi'].append({
            "Psi": np.float64(abs(Psi(0, zi))**2),
            "z": zi
        })
        json.dump(FW_Psi, open(os.path.abspath(path + 'FW3_Psi.json'), 'w'))
        cnt = cnt + 1


def cache_FW6_Jx():

    Z = np.linspace(0, L, 70)
    cnt = 0
    eps = epsilonR2(M1_lentes, 1)
    rho = np.linspace(-float(spot), float(spot), 70)

    jx_FW6_01_L3 = {"jx_FW6_01_L3":[]}
    jx_FW6_05_L3 = {"jx_FW6_05_L3":[]}

    jx_FW6_01_L2 = {"jx_FW6_01_L2":[]}
    jx_FW6_05_L2 = {"jx_FW6_05_L2":[]}

    jx_FW6_01_L5 = {"jx_FW6_01_L5":[]}
    jx_FW6_05_L5 = {"jx_FW6_05_L5":[]}
    
    for rhoi in rho:
        print("i:", cnt)

        jx_FW6_01_L3['jx_FW6_01_L3'].append({
            "jx": np.float64(I_x(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.125)),
            "rho": rhoi
        })
        json.dump(jx_FW6_01_L3, open(os.path.abspath(path + 'jx_FW6_01_L3.json'), 'w'))

        jx_FW6_05_L3['jx_FW6_05_L3'].append({
            "jx": np.float64(I_x(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.125)),
            "rho": rhoi
        })
        json.dump(jx_FW6_05_L3, open(os.path.abspath(path + 'jx_FW6_05_L3.json'), 'w'))

        jx_FW6_01_L2['jx_FW6_01_L2'].append({
            "jx": np.float64(I_x(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.175)),
            "rho": rhoi
        })
        json.dump(jx_FW6_01_L2, open(os.path.abspath(path + 'jx_FW6_01_L2.json'), 'w'))
        jx_FW6_05_L2['jx_FW6_05_L2'].append({
            "jx": np.float64(I_x(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.175)),
            "rho": rhoi
        })
        json.dump(jx_FW6_05_L2, open(os.path.abspath(path + 'jx_FW6_05_L2.json'), 'w'))

        '''
        jx_FW6_01_L5['jx_FW6_01_L5'].append({
            "jx": np.float64(I_x(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.225)),
            "rho": rhoi
        })
        json.dump(jx_FW6_01_L5, open(os.path.abspath(path + 'jx_FW6_01_L5.json'), 'w'))
        jx_FW6_05_L5['jx_FW6_05_L5'].append({
            "jx": np.float64(I_x(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.225)),
            "rho": rhoi
        })
        json.dump(jx_FW6_05_L5, open(os.path.abspath(path + 'jx_FW6_05_L5.json'), 'w'))
        '''

        cnt = cnt + 1

def cache_FW6_Jz():

    eps = epsilonR2(M1_lentes, 1)
    zval = np.linspace(0, L, 70)
    cnt = 0

    jz_FW6_01 = {"jz_FW6_01":[]}
    jz_FW6_1 = {"jz_FW6_1":[]}
    jz_FW6_3 = {"jz_FW6_1":[]}

    #adicionar verificação se o arquivo ja existe 
    #parou no 85
    #print(zval[41:])
    
    for zi in zval:
        print("i:", cnt)

        jz_FW6_01['jz_FW6_01'].append({
            "jz": np.float64(I_z(M1_lentes, 0.1, eps, 0, 0, 0, zi)),
            "z": zi
        })
        json.dump(jz_FW6_01, open(os.path.abspath(path + 'jz_FW6_01.json'), 'w'))

        jz_FW6_1['jz_FW6_1'].append({
            "jz": np.float64(I_z(M1_lentes, 1, eps, 0, 0, 0, zi)),
            "z": zi
        })
        json.dump(jz_FW6_1, open(os.path.abspath(path + 'jz_FW6_1.json'), 'w'))

        jz_FW6_3['jz_FW6_3'].append({
            "jz": np.float64(I_z(M1_lentes, 3, eps, 0, 0, 0, zi)),
            "z": zi
        })
        json.dump(jz_FW6_3, open(os.path.abspath(path + 'jz_FW6_3.json'), 'w'))


        cnt = cnt + 1



def cache_FW_Jz_x():

    eps = epsilonR2(M1_lentes, 1)
    cnt = 0
    xval = np.linspace(0.1, 20, 70)
    '''
    jz_BB_x_M = {"jz_BB_x_M_5":[]}
    jz_BB_x_M1 = {"jz_BB_x_M1_5":[]}
    jz_BB_x_M2 = {"jz_BB_x_M2_5":[]}
    jz_BB_x_M3 = {"jz_BB_x_M3_5":[]}
    '''

    jz_FW_x_M1 = {"jz_FW_x_M1":[]}

    '''
    with open(os.path.abspath(path + 'jz_FW_x_M1.json'), 'r') as jz_FW:
        jz_FW = json.load(jz_FW)
    '''

    #calculo do espectro de Iz, ao longo de diferentes tamanhos de particulas, para diferentes indices de refração 
    for xi in xval:
        print("i:", cnt)

        jz_FW_x_M1['jz_FW_x_M1'].append({
            "jz": np.float64(I_z(M, xi, eps, 0, math.radians(5), 0, 0)),
            "x": xi
        })
        json.dump(jz_FW_x_M1, open(os.path.abspath(path + 'jz_FW_x_M1.json'), 'w'))
        
        cnt = cnt + 1

def cache_FW1_Jz():

    eps = epsilonR2(M1_lentes, 1)
    zval = np.linspace(0, 400*10**(-6), 70)
    cnt = 0

    jz1_FW1_01 = {"jz1_FW1_01":[]}

    '''
    with open(os.path.abspath(path + 'jz1_FW1_01.json'), 'r') as jz_FW:
        jz_FW = json.load(jz_FW)
    '''

    #adicionar verificação se o arquivo ja existe 
    #parou no 85
    #print(zval[41:])
    
    for zi in zval:
        print("i:", cnt)

        jz1_FW1_01['jz1_FW1_01'].append({
            "jz": np.float64(I_z(M1_lentes, 0.1, eps, 0, 0, 0, zi)),
            "z": zi
        })
        json.dump(jz1_FW1_01, open(os.path.abspath(path + 'jz1_FW1_01.json'), 'w'))

        cnt = cnt + 1

def cache_load_FW6_Psi():
    
    Psi = []
    FW_Fz6 = []
    Ztorch = torch.linspace(0, L, 500)
    z_axis = np.linspace(0, 35, 500)
    count = np.linspace(0, 499, 500, dtype=int)

    with open(os.path.abspath(path + 'FW_Psi.json'), 'r') as FW_Psi:
        FW_Psi = json.load(FW_Psi)

    for zt in Ztorch:
        FW_Fz6.append(abs(F_z6(zt))**2)

    
    for i in count:
        Psi.append(FW_Psi['FW_Psi'][i]['Psi'])


    plt.plot(z_axis, Psi,"r", label = r'$|\Psi (0, z)|^2$')
    plt.plot(z_axis, torch.as_tensor(FW_Fz6).cpu(),"k--", label = r'$|F(z)|^2$')

    plt.title(r'$|\Psi (0, z)|^2$')
    plt.xlabel(r'$z (cm)$')
    plt.grid()
    plt.legend()
    #plt.savefig(path, dpi=500)
    plt.show()

def cache_load_FW4_Psi():
    
    Psi = []
    FW_Fz6 = []
    Ztorch = torch.linspace(0, L, 500)
    z_axis = np.linspace(0, 35, 500)
    count = np.linspace(0, 499, 500, dtype=int)

    with open(os.path.abspath(path + 'FW4_Psi.json'), 'r') as FW_Psi:
        FW_Psi = json.load(FW_Psi)

    for zt in Ztorch:
        FW_Fz6.append(abs(F_z4(zt))**2)

    
    for i in count:
        Psi.append(FW_Psi['FW_Psi'][i]['Psi'])


    plt.plot(z_axis, Psi,"r", label = r'$|\Psi (0, z)|^2$')
    plt.plot(z_axis, torch.as_tensor(FW_Fz6).cpu(),"k--", label = r'$|F(z)|^2$')

    plt.title(r'$|\Psi (0, z)|^2$')
    plt.xlabel(r'$z (cm)$')
    plt.grid()
    plt.legend()
    #plt.savefig(path, dpi=500)
    plt.show()

def cache_load_FW3_Psi():
    
    Psi = []
    FW_Fz6 = []
    Ztorch = torch.linspace(0, L, 500)
    z_axis = np.linspace(0, 35, 500)
    count = np.linspace(0, 499, 500, dtype=int)

    with open(os.path.abspath(path + 'FW3_Psi.json'), 'r') as FW_Psi:
        FW_Psi = json.load(FW_Psi)

    for zt in Ztorch:
        FW_Fz6.append(abs(F_z3(zt))**2)

    
    for i in count:
        Psi.append(FW_Psi['FW_Psi'][i]['Psi'])


    plt.plot(z_axis, Psi,"r", label = r'$|\Psi (0, z)|^2$')
    plt.plot(z_axis, torch.as_tensor(FW_Fz6).cpu(),"k--", label = r'$|F(z)|^2$')

    plt.title(r'$|\Psi (0, z)|^2$')
    plt.xlabel(r'$z (cm)$')
    plt.grid()
    plt.legend()
    #plt.savefig(path, dpi=500)
    plt.show()
    
def cache_load_FW1_JZ():
    
    with open(os.path.abspath(path + 'jz1_FW1_01.json'), 'r') as jz_FW:
        jz_FW = json.load(jz_FW)

    jz = []
    #z_axis = np.linspace(0, 400, 100)
    count = np.linspace(0, 69, 100, dtype=int)
    z_axis = []
    print(jz_FW)
    print(len(jz_FW))
    #M = 1.57 - 0.038j
    for i in count:
        jz.append(jz_FW['jz1_FW1_01'][i]['jz'])
        z_axis.append(jz_FW['jz1_FW1_01'][i]['z'])


    #plt.figure(figsize=[10,8])
    plt.plot(z_axis, jz,"b")

    plt.title(r'$J_z$')
    plt.xlabel(r'$z (\mu m)$')
    plt.grid()
    #plt.savefig(path, dpi=500)
    plt.show()

def cache_FW6_Jx_2():

    cnt = 0
    eps = epsilonR2(M1_lentes, 1)
    rho = np.linspace(-2*spot, 2*spot, 70)

    jx_FW6_01_L3 = {"jx_FW6_01_L3":[]}
    jx_FW6_05_L3 = {"jx_FW6_05_L3":[]}
    jx_FW6_1_L3 = {"jx_FW6_1_L3":[]}

    jx_FW6_01_L2 = {"jx_FW6_01_L2":[]}
    jx_FW6_05_L2 = {"jx_FW6_05_L2":[]}
    jx_FW6_1_L2 = {"jx_FW6_1_L2":[]}

    jx_FW6_01_L5 = {"jx_FW6_01_L5":[]}
    jx_FW6_05_L5 = {"jx_FW6_05_L5":[]}
    jx_FW6_1_L5 = {"jx_FW6_1_L5":[]}
    
    for rhoi in rho:
        print(L)
        print("i:", cnt)
        
        
        #x = 0.1
        jx_FW6_01_L3['jx_FW6_01_L3'].append({
            "jx": np.float64(I_x(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.125)),
            "rho": rhoi
        })
        json.dump(jx_FW6_01_L3, open(os.path.abspath(path + 'jx_FW6_01_L3_2.json'), 'w'))

         jx_FW6_01_L5['jx_FW6_01_L5'].append({
            "jx": np.float64(I_x(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.225)),
            "rho": rhoi
        })
        json.dump(jx_FW6_01_L5, open(os.path.abspath(path + 'jx_FW6_01_L5.json'), 'w'))

        jx_FW6_01_L2['jx_FW6_01_L2'].append({
            "jx": np.float64(I_x(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.175)),
            "rho": rhoi
        })
        json.dump(jx_FW6_01_L2, open(os.path.abspath(path + 'jx_FW6_01_L2_2.json'), 'w'))


        '''
        # x = 0.5
         jx_FW6_05_L3['jx_FW6_05_L3'].append({
            "jx": np.float64(I_x(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.125)),
            "rho": rhoi
        })
        json.dump(jx_FW6_05_L3, open(os.path.abspath(path + 'jx_FW6_05_L3_2.json'), 'w'))
        jx_FW6_05_L5['jx_FW6_05_L5'].append({
            "jx": np.float64(I_x(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.225)),
            "rho": rhoi
        })
        jx_FW6_05_L2['jx_FW6_05_L2'].append({
            "jx": np.float64(I_x(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.175)),
            "rho": rhoi
        })
        json.dump(jx_FW6_05_L2, open(os.path.abspath(path + 'jx_FW6_05_L2_2.json'), 'w'))
        '''
        
      '''
        #x = 1 
        jx_FW6_1_L2['jx_FW6_1_L2'].append({
            "jx": np.float64(I_x(M1_lentes, 1, eps, 0, 0, rhoi, 0.175)),
            "rho": rhoi
        })
        json.dump(jx_FW6_1_L2, open(os.path.abspath(path + 'jx_FW6_1_L2.json'), 'w'))
        

        jx_FW6_1_L3['jx_FW6_1_L3'].append({
            "jx": np.float64(I_x(M1_lentes, 1, eps, 0, 0, rhoi, 0.125)),
            "rho": rhoi
        })
        json.dump(jx_FW6_1_L3, open(os.path.abspath(path + 'jx_FW6_1_L3.json'), 'w'))
        
        json.dump(jx_FW6_05_L5, open(os.path.abspath(path + 'jx_FW6_05_L5.json'), 'w'))
        jx_FW6_1_L5['jx_FW6_1_L5'].append({
            "jx": np.float64(I_x(M1_lentes, 1, eps, 0, 0, rhoi, 0.225)),
            "rho": rhoi
        })
        json.dump(jx_FW6_1_L5, open(os.path.abspath(path + 'jx_FW6_1_L5.json'), 'w'))
        
        '''
        cnt = cnt + 1

def cache_FW6_Jx_2_pc2():

    cnt = 0
    eps = epsilonR2(M1_lentes, 1)
    rho = np.linspace(-2*spot, 2*spot, 70)

    jx_FW6_01_L3 = {"jx_FW6_01_L3":[]}
    jx_FW6_05_L3 = {"jx_FW6_05_L3":[]}
    jx_FW6_1_L3 = {"jx_FW6_1_L3":[]}

    jx_FW6_01_L2 = {"jx_FW6_01_L2":[]}
    jx_FW6_05_L2 = {"jx_FW6_05_L2":[]}
    jx_FW6_1_L2 = {"jx_FW6_1_L2":[]}

    jx_FW6_01_L5 = {"jx_FW6_01_L5":[]}
    jx_FW6_05_L5 = {"jx_FW6_05_L5":[]}
    jx_FW6_1_L5 = {"jx_FW6_1_L5":[]}
    
    for rhoi in rho:
        print(L)
        print("i:", cnt)
        
        '''
        #x = 0.1
        jx_FW6_01_L3['jx_FW6_01_L3'].append({
            "jx": np.float64(I_x(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.125)),
            "rho": rhoi
        })
        json.dump(jx_FW6_01_L3, open(os.path.abspath(path + 'jx_FW6_01_L3_2.json'), 'w'))

         jx_FW6_01_L5['jx_FW6_01_L5'].append({
            "jx": np.float64(I_x(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.225)),
            "rho": rhoi
        })
        json.dump(jx_FW6_01_L5, open(os.path.abspath(path + 'jx_FW6_01_L5.json'), 'w'))

        jx_FW6_01_L2['jx_FW6_01_L2'].append({
            "jx": np.float64(I_x(M1_lentes, 0.1, eps, 0, 0, rhoi, 0.175)),
            "rho": rhoi
        })
        json.dump(jx_FW6_01_L2, open(os.path.abspath(path + 'jx_FW6_01_L2_2.json'), 'w'))
        '''

        
        # x = 0.5
         jx_FW6_05_L3['jx_FW6_05_L3'].append({
            "jx": np.float64(I_x(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.125)),
            "rho": rhoi
        })
        json.dump(jx_FW6_05_L3, open(os.path.abspath(path + 'jx_FW6_05_L3_2.json'), 'w'))
        jx_FW6_05_L5['jx_FW6_05_L5'].append({
            "jx": np.float64(I_x(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.225)),
            "rho": rhoi
        })
        jx_FW6_05_L2['jx_FW6_05_L2'].append({
            "jx": np.float64(I_x(M1_lentes, 0.5, eps, 0, 0, rhoi, 0.175)),
            "rho": rhoi
        })
        json.dump(jx_FW6_05_L2, open(os.path.abspath(path + 'jx_FW6_05_L2_2.json'), 'w'))
        
        
      '''
        #x = 1 
        jx_FW6_1_L2['jx_FW6_1_L2'].append({
            "jx": np.float64(I_x(M1_lentes, 1, eps, 0, 0, rhoi, 0.175)),
            "rho": rhoi
        })
        json.dump(jx_FW6_1_L2, open(os.path.abspath(path + 'jx_FW6_1_L2.json'), 'w'))
        

        jx_FW6_1_L3['jx_FW6_1_L3'].append({
            "jx": np.float64(I_x(M1_lentes, 1, eps, 0, 0, rhoi, 0.125)),
            "rho": rhoi
        })
        json.dump(jx_FW6_1_L3, open(os.path.abspath(path + 'jx_FW6_1_L3.json'), 'w'))
        
        json.dump(jx_FW6_05_L5, open(os.path.abspath(path + 'jx_FW6_05_L5.json'), 'w'))
        jx_FW6_1_L5['jx_FW6_1_L5'].append({
            "jx": np.float64(I_x(M1_lentes, 1, eps, 0, 0, rhoi, 0.225)),
            "rho": rhoi
        })
        json.dump(jx_FW6_1_L5, open(os.path.abspath(path + 'jx_FW6_1_L5.json'), 'w'))
        
        '''
        cnt = cnt + 1
    
#cache_FW1_Jz() 
#cache_FW_Jz_x()
#cache_load_FW1_JZ()
#cache_FW_Psi()
#cache_load_FW6_Psi()

#print(I_x(M1_lentes, 0.1, epsilonR2(M1_lentes, 1), 0, 0, spot/2, 0.15)) #tem que dar 0 
#print(I_x(M1_lentes, 0.1, epsilonR2(M1_lentes, 1), 0, 0, spot/2, 0.125)) #tem que dar alguma coisa 
#print(I_x(M1_lentes, 0.1, epsilonR2(M1_lentes, 1), 0, 0, spot/2, 0.23))  #tem que dar alguma coisa 

#cache_FW6_Jx()
cache_FW6_Jx_2()
cache_FW6_Jx_2_pc2()

