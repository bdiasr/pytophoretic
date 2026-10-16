# -*- coding: utf-8 -*-
# FORMULACAO ANTIGA (Ambrosio 2022 / Sun 2022) de I_x, I_y -- congelada em 01/10/2026,
# quando asymetryVector.py passou a usar as expressoes de Pereira como padrao.
# Usar SO para comparacao/reproducao de resultados historicos. I_z e' identico.
import numpy as np
import math


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

    # phi_0 e' agora parametro (default 0.0, o comportamento anterior).
    # Para v >= 1 o feixe tem fase azimutal e^{i v phi}, entao a dependencia em
    # phi_0 deixa de ser trivial: deslocar em phi_0 = 0 e em pi/2 ve campos
    # genuinamente diferentes.
    x_0 = 0.0
    y_0 = 0.0

    const = (3 * epislon) / ((abs(M) ** 2) * (x ** 3))

    ceillingX = (lambda x: math.ceil(x + 4.05 * x ** (1 / 3) + 2))
    nMax = ceillingX(x)

    total = 0.0 + 0.0j

    for n in range(1, nMax + 1):

        a1 = A_mn(n, M, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0, beam=beam)
        b = B_mn(n, M, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0, beam=beam)
        c = 1j * C_mnx(n, M, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0, beam=beam)

        sn = S_n(M, n, x)
        sn_conj = np.conjugate(sn)
        r_np1 = r_n(M, n + 1, x)
        r_n_n = r_n(M, n, x)

        soma = (
            a1 * (sn_conj + ((n + 1) / M) * r_np1) +
            b * (-sn + ((n + 1) / M) * r_n_n) +
            c * sn
        )
        total += soma

    return const * total.imag


def I_y(M, x, epislon, v, alpha, rho_0, z_0, phi_0=0.0, beam=None):

    # phi_0 e' agora parametro (default 0.0, o comportamento anterior).
    # Para v >= 1 o feixe tem fase azimutal e^{i v phi}, entao a dependencia em
    # phi_0 deixa de ser trivial: deslocar em phi_0 = 0 e em pi/2 ve campos
    # genuinamente diferentes.
    x_0 = 0.0
    y_0 = 0.0

    const = (3 * epislon) / ((abs(M) ** 2) * (x ** 3))

    ceillingX = (lambda x: math.ceil(x + 4.05 * x ** (1 / 3) + 2))
    nMax = ceillingX(x)

    total = 0.0 + 0.0j

    for n in range(1, nMax + 1):

        a1 = A_mn(n, M, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0, beam=beam)
        b = B_mn(n, M, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0, beam=beam)
        c = 1j * C_mny(n, M, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0, beam=beam)

        sn = S_n(M, n, x)
        sn_conj = np.conjugate(sn)
        r_np1 = r_n(M, n + 1, x)
        r_n_n = r_n(M, n, x)

        soma = (
            a1 * (sn_conj + ((n + 1) / M) * r_np1) +
            b * (-sn + ((n + 1) / M) * r_n_n) +
            c * sn
        )
        total += soma

    return const * total.real
