import numpy as np
import math

# I_x e I_y pelas expressoes de D. Pereira
# (asymetryVector_pereira); I_z inalterado. Antigo: asymetryVector_antigo.py

# NOTA (2026-09): argumento opcional beam= escolhe o feixe incidente (fonte dos
# BSCs). beam=None (default) = Frozen Wave linear em x de beamShape.py, pelo
# MESMO caminho de antes. beam=BesselBeam(alpha, v, lamb) = feixe de Bessel
# circularmente simetrico do TCC. Ver beams.py / besselBeam.py.
#
# NOTA (2026-09): I_x/I_y/I_z aceitam phi_0 como ultimo argumento (default 0.0).
# Antes estava fixo em 0 dentro de cada funcao. Chamadas posicionais existentes
# continuam validas.

from coefficients import *
from cte import k


def I_z(M, x, epislon, v, alpha, rho_0, z_0, phi_0=0.0, beam=None):

    # phi_0 e' agora parametro (default 0.0, o comportamento anterior).
    # Para v >= 1 o feixe tem fase azimutal e^{i v phi}, entao a dependencia em
    # phi_0 deixa de ser trivial: deslocar em phi_0 = 0 e em pi/2 ve campos
    # genuinamente diferentes.
    x_0 = 0.0
    y_0 = 0.0

    const = (6 * epislon) / ((abs(M) ** 2) * (x ** 3))

    ceillingX = (lambda x: math.ceil(x + 4.05 * x ** (1 / 3) + 2))
    nMax = ceillingX(x)

    total = 0.0 + 0.0j

    for n in range(1, nMax + 1):
        # These are constant for fixed n
        sn = S_n(M, n, x)
        sn_conj = np.conjugate(sn)
        r_n_np1 = r_n(M, n + 1, x)
        r_n_n = r_n(M, n, x)

        b_const = sn_conj + (((n + 1) / M) * r_n_np1)
        d_const = (-sn) + (((n + 1) / M) * r_n_n)

        for m in range(-n, n + 1):
            a1 = D_mn(m, n, M, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0, beam=beam)
            c = E_mn(m, n, M, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0, beam=beam)
            e = 1j * F_mn(m, n, M, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0, beam=beam)

            total += (a1 * b_const) + (c * d_const) + (e * sn)

    return const * total.imag


def I_x(M, x, epislon, v, alpha, rho_0, z_0, phi_0=0.0, beam=None):
    """J_x pelas expressoes CORRIGIDAS de D. Pereira.
    Formulacao antiga congelada em asymetryVector_antigo.py."""
    from asymetryVector_pereira import I_xy  # import local: evita ciclo de importacao
    return I_xy(M, x, epislon, v, alpha, rho_0, z_0, phi_0, beam)[0]


def I_y(M, x, epislon, v, alpha, rho_0, z_0, phi_0=0.0, beam=None):
    """J_y pelas expressoes CORRIGIDAS de D. Pereira."""
    from asymetryVector_pereira import I_xy
    return I_xy(M, x, epislon, v, alpha, rho_0, z_0, phi_0, beam)[1]
