import torch, time
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import matplotlib as mpl

from numba import njit, jit

from Fz import *
from cte import * 
from Aq import * 
from frozenWave import *
from legendre import * 
from beamShape import g_mnTE, g_mnTM, gn_FrozenWave, g_mnTE_FW, g_mnTM_FW



def plot_Psi():

    TEST_DPI = 50
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
    mpl.rcParams['figure.dpi'] = TEST_DPI

    Psi_values = []
    F_values = []

    Z = np.linspace(0, L, 500)
    Ztorch = torch.linspace(0, L, 500)

    for z in Z: 
        Psi_values.append(abs(Psi(0, z))**2)
    
    for zt in Ztorch:
        F_values.append(abs(F_z6(zt))**2)


    plt.plot(Ztorch.cpu(), torch.as_tensor(F_values).cpu(), 'k--')
    plt.plot(Ztorch.cpu(), torch.as_tensor(Psi_values).cpu(), 'r')
    plt.title(r'$|\Psi (0, z)|^2$')
    plt.xlabel(r'$z (cm)$')
    plt.grid()
    plt.show()


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

    n1 = []
    n2 = []
    n3 = []

    Z = np.linspace(0, L, 100)
    rho = np.linspace(-float(spot), float(spot), 100)

    for rhoi in rho:
        #print(z)
        n1.append(
            abs(g_mnTE(1, 1, 0, 0, lamb, rhoi, 0, 0, 0, L/3)))
        n2.append(
            abs(g_mnTE(1, 2, 0, 0, lamb, rhoi, 0, 0, 0, L/3)))
        n3.append(
            abs(g_mnTE(1, 3, 0, 0, lamb, rhoi, 0, 0, 0, L/3)))

    plt.plot(rho, n1, 'b', label = 'n = 1')
    plt.plot(rho, n2, 'r', label = 'n = 2')
    plt.plot(rho, n3, 'g', label = 'n = 3')
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel(r'$z_0(\mu m)$')
    plt.title(r'$|g_{TE, n}^m|$')
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
#plot_gnTE()
#plot_gnTM()


#fig_pi()
#fig_tau()
