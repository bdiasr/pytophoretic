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






def plot_Fz():
    
    """ Figure DPI """
    TEST_DPI = 100
    mpl.rcParams["lines.linewidth"] = 2
    mpl.rcParams["text.usetex"] = True
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["axes.titlesize"] = 32
    mpl.rcParams["axes.labelsize"] = 24
    mpl.rcParams["axes.titlepad"] = 24
    mpl.rcParams["xtick.labelsize"] = "xx-large"
    mpl.rcParams["ytick.labelsize"] = "xx-large"
    mpl.rcParams["legend.fontsize"] = 24
    mpl.rcParams['figure.figsize'] = [12, 8]
    mpl.rcParams['figure.dpi'] = TEST_DPI # Use 200 for fine images - note this is slow

    
    F_values = []
    #L = (1000*10**(-6))
    L = 400*10**(-6)

    Z = np.linspace(0, L, 100)

    print('Plot de F_z com Pytorch')
    start_time = time.time()

    for z in Z: 
        F_values.append(abs(F_z4(z, L))**2)

    end_time = time.time()
    print('Tempo de execução (versão pytorch): {} segundos'.format(end_time - start_time))

    plt.plot(Z, torch.as_tensor(F_values).cpu(), 'r')
    plt.xlabel(r'$z$')
    plt.title(r'$F(z)$')
    plt.grid()
    plt.show()


def plot_gnTE():

    print("plot_gnTE iniciado")

    n1 = []
    n2 = []
    n3 = []

    Z = np.linspace(0, 0.35, 100)

    for zi in Z:
        #print(zi)
        n1.append(
            cmath.phase(g_mnTE(1, 1, 0, 0, lamb, -0.00005, 0, 0, 0, zi)))
        n2.append(
            cmath.phase(g_mnTE(1, 1, 0, 0, lamb, 0.00005, np.pi, 0, 0, zi)))

    plt.plot(Z, n1, 'b', label = 'rho < 0')
    plt.plot(Z, n2, 'k--', label = 'rho > 0')
    #plt.plot(Z, n3, 'g', label = 'n = 3')
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel(r'$z_0(\mu m)$')
    plt.title(r'$ARG(g_{TE, n}^m), \rho_0 = -1$')
    plt.show()

def plot_gnTE1():

    n1 = []
    n2 = []
    n3 = []
    print("plot_gnTE1 iniciado")

    Z = np.linspace(0, 0.35, 100)

    for zi in Z:
        #print(z)
        n1.append(
            cmath.phase(g_mnTE(1, 1, 0, 0, lamb, 0.00005, np.pi, 0, 0, zi)))
        n2.append(
            cmath.phase(g_mnTE(1, 2, 0, 0, lamb, 0.00005, np.pi, 0, 0, zi)))
        n3.append(
            cmath.phase(g_mnTE(1, 3, 0, 0, lamb, 0.00005, np.pi, 0, 0, zi)))

    plt.plot(Z, n1, 'b', label = 'n = 1')
    plt.plot(Z, n2, 'r', label = 'n = 2')
    plt.plot(Z, n3, 'g', label = 'n = 3')
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel(r'$z_0(\mu m)$')
    plt.title(r'$ARG(g_{TE, n}^m), \rho_0 = 1$')
    plt.show()

def plot_gnTM():

    n1 = []
    n2 = []
    n3 = []

    Z = np.linspace(0, L, 100)

    for z in Z:
        #print(z)
        n1.append(
            abs(g_mnTM(1, 1, 0, 0, lamb, 0, 0, 0, 0, z)))
        n2.append(
            abs(g_mnTM(1, 2, 0, 0, lamb, 0, 0, 0, 0, z)))
        n3.append(
            abs(g_mnTM(1, 3, 0, 0, lamb, 0, 0, 0, 0, z)))

    plt.plot(Z, n1, 'b', label = 'n = 1')
    plt.plot(Z, n2, 'r', label = 'n = 2')
    plt.plot(Z, n3, 'g', label = 'n = 3')
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel(r'$z_0(\mu m)$')
    plt.title(r'$|g_{TM, n}^m|$')
    plt.show()



def plot_gnOn():

    n1 = []
    n2 = []
    n3 = []

    Z = np.linspace(0, L, 100)

    for zi in Z:
        n1.append(abs(gn_FrozenWave(1, -27, k, L, Q, zi)))
        n2.append(abs(gn_FrozenWave(2, -27, k, L, Q, zi)))
        n3.append(abs(gn_FrozenWave(3, -27, k, L, Q, zi)))

    print(n1)
    print(n2)
    print(n3)


    plt.plot(Z, n1, 'b')
    plt.plot(Z, n2, 'r')
    plt.plot(Z, n3, 'g')
    plt.xlabel(r'$z_0(\mu m)$')
    plt.title(r'$|g_n|$')
    plt.show()



def fig_tau():

    n1 = []
    n2 = []
    n3 = []

    X = np.linspace(0.0001, 0.99, 100)

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

    X = np.linspace(0.0001, 1, 100)

    for x in X:
        n1.append(pi_mn(1, 1, x))

    for x in X:
        n2.append(pi_mn(1, 2, x))
  
    for x in X:
        n3.append(pi_mn(1, 3, x))

    plt.plot(X, n1, 'b', label ='n1')
    plt.plot(X, n2, 'r', label='n2')
    plt.plot(X, n3, 'g', label='n3')

    #plt.ylim(-6, 2)
    plt.xlabel('x')
    #plt.ylabel('')
    plt.legend()
    plt.grid()
    plt.show()


"""
Funcoes de teste
"""
#plot_Psi()
#plot_gnOn()
#plot_gnTM()


#fig_pi()
#fig_tau()

#plot_gnTE()
#plot_gnTE1()

#F_z5 e F_z2
def plot_Psi():

    Psi_values = []
    F_values = []

    Z = np.linspace(0, L, 500)
    Ztorch = torch.linspace(0, L, 500)

    for z in Z: 
        Psi_values.append(abs(Psi(0, z))**2)
    
    for zt in Ztorch:
        F_values.append(abs(F_z2(zt))**2)


    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(26)
    ax.xaxis.get_offset_text().set_fontsize(26)
    plt.plot(Ztorch.cpu(), torch.as_tensor(F_values).cpu(), 'k--', label = r'$F(z)$')
    plt.plot(Ztorch.cpu(), torch.as_tensor(Psi_values).cpu(), 'r', label = r'$|\Psi (0, z)|^2$')
    plt.title(r'$|\Psi (0, z)|^2$', fontsize=26)
    plt.xlabel(r'$z \text{ }(\text{m})$', fontsize=26)
    plt.legend(prop = { "size": 18 })
    plt.tick_params(axis='both', which='major', labelsize=26)
    plt.tick_params(axis='both', which='minor', labelsize=26)
    plt.grid()
    plt.savefig(r'C:Users/Administrador/Documents/USP/pytophoretic/figs/PSI.svg', format='svg', dpi=1000)
    plt.savefig(r'C:Users/Administrador/Documents/USP/pytophoretic/figs/PSI.png', format='png', dpi=1000)
    plt.show()


'''funcoes artigo
#F_z5 e F_z2
'''
#ta certo 
def plot_Psi_1():

    Psi_values = []
    F_values = []

    Z = np.linspace(0, L, 200)
    Ztorch = torch.linspace(0, L, 200)

    for z in Z: 
        Psi_values.append(abs(Psi_1(0, z))**2)
    
    for zt in Ztorch:
        F_values.append(abs(F_z5(zt))**2)


    ax = plt.gca()
    #ax.yaxis.get_offset_text().set_fontsize(26)
    #ax.xaxis.get_offset_text().set_fontsize(26)
    plt.plot(Ztorch.cpu(), torch.as_tensor(F_values).cpu(), 'k--', label = r'$F(z)$')
    plt.plot(Ztorch.cpu(), torch.as_tensor(Psi_values).cpu(), 'r', label = r'$|\Psi (0, z)|^2$')
    plt.title(r'$|\Psi (0, z)|^2$', fontsize=26)
    plt.xlabel(r'$z \text{ }(\text{m})$', fontsize=26)
    #plt.legend(prop = { "size": 18 })
    #plt.tick_params(axis='both', which='major', labelsize=26)
    #plt.tick_params(axis='both', which='minor', labelsize=26)
    plt.grid()

    # Criar o diretório se ele não existir
    output_dir = r'figs/artigo/'
    os.makedirs(output_dir, exist_ok=True)

    plt.savefig(r'C:Users/Administrador/Documents/USP/pytophoretic/figs/artigo/PSI_1.svg', format='svg', dpi=1000)
    plt.savefig(r'C:Users/Administrador/Documents/USP/pytophoretic/figs/artigo/PSI_1.png', format='png', dpi=1000)
    plt.show()

#ta errado 
def plot_Psi_2():

    Psi_values = []
    F_values = []

    Z = np.linspace(0, L, 200)
    Ztorch = torch.linspace(0, L, 200)

    for z in Z: 
        Psi_values.append(abs(Psi_2(0, z))**2)
    
    for zt in Ztorch:
        F_values.append(abs(Fz(zt))**2)


    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(26)
    ax.xaxis.get_offset_text().set_fontsize(26)
    plt.plot(Z, torch.as_tensor(F_values).cpu(), 'k--', label = r'$F(z)$')
    plt.plot(Z, torch.as_tensor(Psi_values).cpu(), 'r', label = r'$|\Psi (0, z)|^2$')
    plt.title(r'$|\Psi (0, z)|^2$', fontsize=26)
    plt.xlabel(r'$z \text{ }(\text{m})$', fontsize=26)
    #plt.legend(prop = { "size": 18 })
    #plt.tick_params(axis='both', which='major', labelsize=26)
    #plt.tick_params(axis='both', which='minor', labelsize=26)
    plt.grid()

    # Criar o diretório se ele não existir
    output_dir = r'figs/artigo/'
    os.makedirs(output_dir, exist_ok=True)

    plt.savefig(r'C:Users/Administrador/Documents/USP/pytophoretic/figs/artigo/PSI_2.svg', format='svg', dpi=1000)
    plt.savefig(r'C:Users/Administrador/Documents/USP/pytophoretic/figs/artigo/PSI_2.png', format='png', dpi=1000)
    plt.show()

def fig_Psi_3d():
    ax = plt.axes(projection="3d")

    x = np.linspace(-spot*20, spot*20, 500)
    z = np.linspace(0, L, 500)
    X, Z = np.meshgrid(x, z)

    Y = abs(Psi(X, Z))**2

    ax.plot_surface(X, Z, Y, cmap='gnuplot', label = r'$|\Psi (0, z)|^2$')

    ax.set_xlabel(r'$\rho \text{ } (\mu \text{m})$',  fontsize=18)
    ax.set_ylabel(r'$z \text{ (m)} $',  fontsize=18)
    ax.set_zlabel(r'$|\Psi (\rho, z)|^2$',  fontsize=18)

    # Ajuste o formato dos ticks dos eixos x, y e z para função de 10^-3
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: r'${:.1f}$'.format(x * 1e6)))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: r'${:.1f}$'.format(y * 1e0)))
    ax.yaxis.get_offset_text().set_fontsize(14)
    ax.xaxis.get_offset_text().set_fontsize(14)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='minor', labelsize=18)
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/pytophoretic/figs/PSI_3D.svg', format='svg', dpi=1000)
    plt.savefig(r'D:/Documentos/USP/pytophoretic/figs/PSI_3D.png', format='png', dpi=1000)

    plt.show()

def plot_Psi_rho1():

    Psi_values = []
    
    x = np.linspace(-spot*3, spot*3, 500)

    for xi in x: 
        Psi_values.append(abs(Psi(xi, 0.125))**2)


    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(26)
    ax.xaxis.get_offset_text().set_fontsize(26)
    #plt.plot(Ztorch.cpu(), torch.as_tensor(F_values).cpu(), 'k--', label = r'$F(z)$')
    plt.plot(x, torch.as_tensor(Psi_values).cpu(), 'r', label = r'$|\Psi (\rho, z = 0.125)|^2$')
    plt.title(r'$|\Psi (\rho, z = 0.125)|^2$', fontsize=26)
    plt.xlabel(r'$z \text{ }(\text{m})$', fontsize=26)
    plt.legend(prop = { "size": 18 })
    plt.tick_params(axis='both', which='major', labelsize=26)
    plt.tick_params(axis='both', which='minor', labelsize=26)
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/pytophoretic/figs/PSI_rho1.svg', format='svg', dpi=1000)
    plt.savefig(r'D:/Documentos/USP/pytophoretic/figs/PSI_rho1.png', format='png', dpi=1000)
    plt.show()

def plot_Psi_rho2():

    Psi_values = []

    x = np.linspace(-spot*3, spot*3, 500)

    for xi in x: 
        Psi_values.append(abs(Psi(xi, 0.175))**2)


    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(26)
    ax.xaxis.get_offset_text().set_fontsize(26)
    plt.plot(x, torch.as_tensor(Psi_values).cpu(), 'r', label = r'$|\Psi (\rho, z = 0.175)|^2$')
    plt.title(r'$|\Psi (\rho, z = 0.175)|^2$', fontsize=26)
    plt.xlabel(r'$z \text{ }(\text{m})$', fontsize=26)
    plt.legend(prop = { "size": 18 })
    plt.tick_params(axis='both', which='major', labelsize=26)
    plt.tick_params(axis='both', which='minor', labelsize=26)
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/pytophoretic/figs/PSI_rho2.svg', format='svg', dpi=1000)
    plt.savefig(r'D:/Documentos/USP/pytophoretic/figs/PSI_rho2.png', format='png', dpi=1000)
    plt.show()

def plot_Psi_rho3():

    Psi_values = []

    x = np.linspace(-spot*3, spot*3, 500)

    for xi in x: 
        Psi_values.append(abs(Psi(xi, 0.225))**2)


    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(26)
    ax.xaxis.get_offset_text().set_fontsize(26)
    plt.plot(x, torch.as_tensor(Psi_values).cpu(), 'r', label = r'$|\Psi (\rho, z = 0.225)|^2$')
    plt.title(r'$|\Psi (\rho, z = 0.225)|^2$', fontsize=26)
    plt.xlabel(r'$z \text{ }(\text{m})$', fontsize=26)
    plt.legend(prop = { "size": 18 })
    plt.tick_params(axis='both', which='major', labelsize=26)
    plt.tick_params(axis='both', which='minor', labelsize=26)
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/pytophoretic/figs/PSI_rho3.svg', format='svg', dpi=1000)
    plt.savefig(r'D:/Documentos/USP/pytophoretic/figs/PSI_rho3.png', format='png', dpi=1000)
    plt.show()


#plot_Psi_rho1()
#plot_Psi_rho2()
#plot_Psi_rho3()

#plot_Psi()
#fig_Psi_3d()
plot_Psi_1()
plot_Psi_2()
