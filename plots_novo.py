import torch, time
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import matplotlib as mpl
from mpl_toolkits import mplot3d
import matplotlib.ticker as ticker

from numba import njit, jit

from Fz import *
from cte import * 
from Aq import * 
from frozenWave import *
from legendre import * 
from beamShape import g_mnTE, g_mnTM, gn_FrozenWave

from matplotlib import cbook, cm
from matplotlib.colors import LightSource

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.rcParams['text.usetex'] = True
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'

path = r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/'


def fig_Psi_1_3d():
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    x = np.linspace(-spot*5, spot*5, 5000)
    z = np.linspace(0, L, 5000)
    X, Z = np.meshgrid(x, z)

    Y = abs(Psi_1(X, Z))**2

    surf = ax.plot_surface(X, Z, Y, cmap='gnuplot', label = r'$|\Psi (0, z)|^2$')

    ax.set_xlabel(r'$\rho \text{ } (\mu \text{m})$',  fontsize=18)
    ax.set_ylabel(r'$z \text{ (m)} $',  fontsize=18)
    ax.set_zlabel(r'$|\Psi (\rho, z)|^2$',  fontsize=18)
    fig.colorbar(surf, shrink=0.3, aspect=5)

    # Ajuste o formato dos ticks dos eixos x, y e z para função de 10^-3
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: r'${:.1f}$'.format(x * 1e6)))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: r'${:.1f}$'.format(y * 1e0)))
    ax.yaxis.get_offset_text().set_fontsize(14)
    ax.xaxis.get_offset_text().set_fontsize(14)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='minor', labelsize=18)
    plt.grid()
    #plt.text(-0.0, 1.02, r'$\times 10^{-4}$', fontsize=18, transform=ax.transAxes)
    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/FW5_3D.svg', format='svg', dpi=1000)
    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/FW5_3D.png', format='png', dpi=1000)

    plt.show()

def fig_Psi_2_3d():
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    x = np.linspace(-spot*5, spot*5, 5000)
    z = np.linspace(0, L, 5000)
    X, Z = np.meshgrid(x, z)

    Y = abs(Psi_2(X, Z))**2

    surf = ax.plot_surface(X, Z, Y, cmap='gnuplot', label = r'$|\Psi (0, z)|^2$')

    ax.set_xlabel(r'$\rho \text{ } (\mu \text{m})$',  fontsize=18)
    ax.set_ylabel(r'$z \text{ (m)} $',  fontsize=18)
    ax.set_zlabel(r'$|\Psi (\rho, z)|^2$',  fontsize=18)
    fig.colorbar(surf, shrink=0.3, aspect=5)

    # Ajuste o formato dos ticks dos eixos x, y e z para função de 10^-3
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: r'${:.1f}$'.format(x * 1e6)))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: r'${:.1f}$'.format(y * 1e0)))
    ax.yaxis.get_offset_text().set_fontsize(14)
    ax.xaxis.get_offset_text().set_fontsize(14)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='minor', labelsize=18)
    plt.grid()
    #plt.text(-0.0, 1.02, r'$\times 10^{-4}$', fontsize=18, transform=ax.transAxes)
    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/FW_3D.svg', format='svg', dpi=1000)
    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/FW_3D.png', format='png', dpi=1000)

    plt.show()


def plot_Psi_1():

    Psi_values = []
    F_values = []

    Z = np.linspace(0, L, 500)
    Ztorch = torch.linspace(0, L, 500)

    for z in Z: 
        Psi_values.append(abs(Psi_1(0, z))**2)
    
    #for zt in Ztorch:
        #F_values.append(abs(F_z5(zt))**2)


    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(18)
    ax.xaxis.get_offset_text().set_fontsize(18)
    #plt.plot(Ztorch.cpu(), torch.as_tensor(F_values).cpu(), 'k--', label = r'$F(z)$')
    plt.plot(Ztorch.cpu(), torch.as_tensor(Psi_values).cpu(), 'r', label = r'$|\Psi (0, z)|^2$')
    plt.title(r'$|\Psi (0, z)|^2$', fontsize=26)
    plt.legend(prop = { "size": 12 })
    plt.xlabel(r'$z \text{ }(\text{m})$', fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='minor', labelsize=18)
    plt.grid()

    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/PSI_1.svg', format='svg', dpi=1000)
    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/PSI_1.png', format='png', dpi=1000)
    plt.show()


def plot_Psi_2():

    Psi_values = []
    F_values = []

    Z = np.linspace(0, L, 500)
    Ztorch = torch.linspace(0, L, 500)

    for z in Z: 
        Psi_values.append(abs(Psi_2(0, z))**2)
    
    #for zt in Ztorch:
        #F_values.append(abs(Fz(zt))**2)


    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(18)
    ax.xaxis.get_offset_text().set_fontsize(18)
    #plt.plot(Z, torch.as_tensor(F_values).cpu(), 'k--', label = r'$F(z)$')
    plt.plot(Z, torch.as_tensor(Psi_values).cpu(), 'r', label = r'$|\Psi (0, z)|^2$')
    plt.title(r'$|\Psi (0, z)|^2$', fontsize=26)
    plt.legend(prop = { "size": 12 })
    plt.xlabel(r'$z \text{ }(\text{m})$', fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='minor', labelsize=18)
    plt.grid()

    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/PSI_2.svg', format='svg', dpi=1000)
    plt.savefig(r'C:\Users\Administrador\Documents\USP\pytophoretic\cache_novo/PSI_2.png', format='png', dpi=1000)
    plt.show()

fig_Psi_1_3d()
fig_Psi_2_3d()

#plot_Psi_1()
#plot_Psi_2()