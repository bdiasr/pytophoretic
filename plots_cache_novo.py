import os
import json 

import numpy as np



import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'
import matplotlib.ticker as ticker

from cte import *
from coefficients import * 
from frozenWave import * 
from asymetryVector import *

import os.path

path = r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/'

def Iz_FW5():

    with open(os.path.abspath(path + 'jz_FW5_01.json'), 'r') as jz1_FW5_01:
        jz1_FW5_01 = json.load(jz1_FW5_01)
    with open(os.path.abspath(path + 'jz_FW5_3.json'), 'r') as jz1_FW5_3:
        jz1_FW5_3 = json.load(jz1_FW5_3)
    with open(os.path.abspath(path + 'jz_FW5_1.json'), 'r') as jz1_FW5_1:
        jz1_FW5_1 = json.load(jz1_FW5_1)

    jz01 = []
    jz3 = []
    jz1 = []
    z_axis = []

    for i in range(len(jz1_FW5_1['jz_FW5_1'])):
        jz01.append((jz1_FW5_01['jz_FW5_01'][i]['jz'])*100000)
        jz3.append((jz1_FW5_3['jz_FW5_3'][i]['jz']))
        jz1.append(jz1_FW5_1['jz_FW5_1'][i]['jz'])
        z_axis.append(jz1_FW5_1['jz_FW5_1'][i]['z'])

    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(18)
    ax.xaxis.get_offset_text().set_fontsize(18)
    plt.plot(z_axis, jz01,"b", label=r'$x = 0.1 \text{ }(\times 10) $')
    plt.plot(z_axis, jz1,"r--", label=r'$x = 1$')
    plt.plot(z_axis, jz3,"k--", label=r'$x = 3$')
    plt.legend(prop = { "size": 16 })
    plt.xlabel(r'$z \text{ }(\text{m})$', fontsize=18)
    plt.ylabel(r'$I_z \text{ }(z\text{ = }0\text{ m}) $', fontsize=18)
    plt.grid()
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='minor', labelsize=18)
    # Ajuste o formato dos ticks do eixo y para função de 10^-3
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: r'${:.1f}$'.format(y * 1e4)))
    ax.yaxis.set_label_coords(-0.1, 0.5)  # ajusta a posição do label do eixo y
    # Adiciona a notação de potência no canto superior esquerdo
    plt.text(-0.0, 1.02, r'$\times 10^{-4}$', fontsize=18, transform=ax.transAxes)
    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/Iz_FW5.svg', format='svg', dpi=1000)
    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/Iz_FW5.png', format='png', dpi=1000)
    plt.show()


def Ix_FW5_L5():
    with open(os.path.abspath(path + 'jx_FW5_01_L5.json'), 'r') as jx_FW5_01_L5:
        jx_FW5_01_L5 = json.load(jx_FW5_01_L5)

    with open(os.path.abspath(path + 'jx_FW5_05_L5.json'), 'r') as jx_FW5_05_L5:
        jx_FW5_05_L5 = json.load(jx_FW5_05_L5)

    #with open(os.path.abspath(path + 'jx_FW5_1_L5.json'), 'r') as jx_FW5_1_L5:
        #jx_FW5_1_L5 = json.load(jx_FW5_1_L5)

    jz01 = []
    jz05 = []
    jz1 = []
    x_axis = []
    


    for i in range(len(jx_FW5_01_L5['jx_FW5_01_L5'])):
        jz01.append(jx_FW5_01_L5['jx_FW5_01_L5'][i]['jx']*10)
        jz05.append(jx_FW5_05_L5['jx_FW5_05_L5'][i]['jx'])
        #jz1.append(jx_FW5_1_L5['jx_FW5_1_L5'][i]['jx'])
        x_axis.append(jx_FW5_05_L5['jx_FW5_05_L5'][i]['rho'])

    Psi_values =[]
    for xi in x_axis: 
        Psi_values.append((abs(Psi(xi, 0.4))**2)*10**(-4))



    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(16)
    ax.xaxis.get_offset_text().set_fontsize(16)
    plt.plot(x_axis, jz01,"b", label=r'$x = 0.1 \text{ }(\times 10) $')
    plt.plot(x_axis, jz05,"k--", label=r'$x = 0.5$')
    #plt.plot(x_axis, jz1,"r--", label=r'$x = 1$')
    plt.plot(x_axis, torch.as_tensor(Psi_values).cpu(), 'k', linewidth=2)
    plt.legend(prop = { "size": 18 })
    plt.xlabel(r'$\rho_0 \text{ }(\mu \text{m})$', fontsize=18)
    plt.ylabel(r'$I_x \text{ }(z\text{ = }0.4\text{ m}) $', fontsize=18)
    plt.grid()
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='minor', labelsize=18)
    # Ajuste o formato dos ticks do eixo y para função de 10^-3
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: r'${:.1f}$'.format(y * 1e4)))
    ax.yaxis.set_label_coords(-0.1, 0.5)  # ajusta a posição do label do eixo y
    # Adiciona a notação de potência no canto superior esquerdo
    plt.text(-0.0, 1.02, r'$\times 10^{-4}$', fontsize=18, transform=ax.transAxes)
    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/Ix_FW5_L5.svg', format='svg', dpi=1000)
    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/Ix_FW5_L5.png', format='png', dpi=1000)
    plt.show()

def Ix_FW5_L2():
    with open(os.path.abspath(path + 'jx_FW5_01_L2.json'), 'r') as jx_FW5_01_L2:
        jx_FW5_01_L2 = json.load(jx_FW5_01_L2)

    with open(os.path.abspath(path + 'jx_FW5_05_L2.json'), 'r') as jx_FW5_05_L2:
        jx_FW5_05_L2 = json.load(jx_FW5_05_L2)

    #with open(os.path.abspath(path + 'jx_FW5_1_L5.json'), 'r') as jx_FW5_1_L5:
        #jx_FW5_1_L5 = json.load(jx_FW5_1_L5)

    jz01 = []
    jz05 = []
    jz1 = []
    x_axis = []
    


    for i in range(len(jx_FW5_01_L2['jx_FW5_01_L2'])):
        jz01.append(jx_FW5_01_L2['jx_FW5_01_L2'][i]['jx']*10)
        jz05.append(jx_FW5_05_L2['jx_FW5_05_L2'][i]['jx'])
        #jz1.append(jx_FW5_1_L5['jx_FW5_1_L5'][i]['jx'])
        x_axis.append(jx_FW5_05_L2['jx_FW5_05_L2'][i]['rho'])

    Psi_values =[]
    for xi in x_axis: 
        Psi_values.append((abs(Psi(xi, 0.25))**2)*10**(-4))



    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(16)
    ax.xaxis.get_offset_text().set_fontsize(16)
    plt.plot(x_axis, jz01,"b", label=r'$x = 0.1 \text{ }(\times 10) $')
    plt.plot(x_axis, jz05,"k--", label=r'$x = 0.25$')
    #plt.plot(x_axis, jz1,"r--", label=r'$x = 1$')
    plt.plot(x_axis, torch.as_tensor(Psi_values).cpu(), 'k', linewidth=2)
    plt.legend(prop = { "size": 18 })
    plt.xlabel(r'$\rho_0 \text{ }(\mu \text{m})$', fontsize=18)
    plt.ylabel(r'$I_x \text{ }(z\text{ = }0.25\text{ m}) $', fontsize=18)
    plt.grid()
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='minor', labelsize=18)
    # Ajuste o formato dos ticks do eixo y para função de 10^-3
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: r'${:.1f}$'.format(y * 1e4)))
    ax.yaxis.set_label_coords(-0.1, 0.5)  # ajusta a posição do label do eixo y
    # Adiciona a notação de potência no canto superior esquerdo
    plt.text(-0.0, 1.02, r'$\times 10^{-4}$', fontsize=18, transform=ax.transAxes)
    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/Ix_FW5_L2.svg', format='svg', dpi=1000)
    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/Ix_FW5_L2.png', format='png', dpi=1000)
    plt.show()

def Ix_FW5_L3():
    with open(os.path.abspath(path + 'jx_FW5_01_L3.json'), 'r') as jx_FW5_01_L3:
        jx_FW5_01_L3 = json.load(jx_FW5_01_L3)

    with open(os.path.abspath(path + 'jx_FW5_05_L3.json'), 'r') as jx_FW5_05_L3:
        jx_FW5_05_L3 = json.load(jx_FW5_05_L3)

    #with open(os.path.abspath(path + 'jx_FW5_1_L5.json'), 'r') as jx_FW5_1_L5:
        #jx_FW5_1_L5 = json.load(jx_FW5_1_L5)

    jz01 = []
    jz05 = []
    jz1 = []
    x_axis = []
    


    for i in range(len(jx_FW5_01_L3['jx_FW5_01_L3'])):
        jz01.append(jx_FW5_01_L3['jx_FW5_01_L3'][i]['jx']*10)
        jz05.append(jx_FW5_05_L3['jx_FW5_05_L3'][i]['jx'])
        #jz1.append(jx_FW5_1_L5['jx_FW5_1_L5'][i]['jx'])
        x_axis.append(jx_FW5_05_L3['jx_FW5_05_L3'][i]['rho'])

    Psi_values =[]
    for xi in x_axis: 
        Psi_values.append((abs(Psi(xi, 0.1))**2)*10**(-4))



    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(16)
    ax.xaxis.get_offset_text().set_fontsize(16)
    plt.plot(x_axis, jz01,"b", label=r'$x = 0.1 \text{ }(\times 10) $')
    plt.plot(x_axis, jz05,"k--", label=r'$x = 0.5$')
    #plt.plot(x_axis, jz1,"r--", label=r'$x = 1$')
    plt.plot(x_axis, torch.as_tensor(Psi_values).cpu(), 'k', linewidth=2)
    plt.legend(prop = { "size": 18 })
    plt.xlabel(r'$\rho_0 \text{ }(\mu \text{m})$', fontsize=18)
    plt.ylabel(r'$I_x \text{ }(z\text{ = }0.1\text{ m}) $', fontsize=18)
    plt.grid()
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='minor', labelsize=18)
    # Ajuste o formato dos ticks do eixo y para função de 10^-3
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: r'${:.1f}$'.format(y * 1e4)))
    ax.yaxis.set_label_coords(-0.1, 0.5)  # ajusta a posição do label do eixo y
    # Adiciona a notação de potência no canto superior esquerdo
    plt.text(-0.0, 1.02, r'$\times 10^{-4}$', fontsize=18, transform=ax.transAxes)
    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/Ix_FW5_L3.svg', format='svg', dpi=1000)
    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/Ix_FW5_L3.png', format='png', dpi=1000)
    plt.show()


#Iz_FW5()
Ix_FW5_L5()
Ix_FW5_L2()
Ix_FW5_L3()