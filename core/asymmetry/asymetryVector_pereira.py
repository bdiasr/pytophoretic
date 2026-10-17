# -*- coding: utf-8 -*-
"""
Transverse components J_x, J_y of the asymmetry vector (D. Pereira's correction)
================================================================================
References:
  D. Pereira, "Expressoes Jx e Jy do artigo 'Towards Photophoresis...'", technical
  note (Oct. 2025), Eqs. (8)-(11);
  D. J. S. Pereira, "3D Photophoresis in Slip-flow", PhD thesis, Univ. Coimbra (2026),
  Eqs. (6.10)-(6.13) and Appendix H.

Same normalization as asymetryVector (returns J_i = 3 eps'' I_i / (8 pi a^3)).
J_z is re-exported unchanged from asymetryVector.

    J_{x;y} = -3 eps'' / (|M|^2 x^3) {Im; Re} sum_n [ A_n (S_n^* + (n+1) R_{n+1} / M)
                                                   + B_n (-S_n + (n+1) R_n / M)
                                                   - i C_n S_n ]          Eq. (6.10)

One function per term, upper sign (x) and lower sign (y) written separately:
    A_n_x, A_n_y  Eq. (6.11)      B_n_x, B_n_y  Eq. (6.12)      C_n_x, C_n_y  Eq. (6.13)

Verification: tests/diagnostico_pereira2025.py (ratio 1 +/- 3e-11 against the
defining volume integral).
"""
import math
from functools import lru_cache

import numpy as np

from coefficients import cs_n, ds_n, S_n, r_n, _fact_ratio
from cte import Eta_r, lamb as _LAMBDA_FW
from beams import resolver_bsc
from asymetryVector import I_z  # noqa: F401  (longitudinal component, unchanged)

MU_R = 1   # relative permeability used in c_n, d_n


def n_max(x):
    """Truncation of the multipole series (Wiscombe)."""
    return math.ceil(x + 4.05 * x ** (1 / 3) + 2)


@lru_cache(maxsize=None)
def fact_ratio(a, b):
    """a!/b!; zero when b < 0 (the term does not exist)."""
    return _fact_ratio(a, b) if b >= 0 else 0.0


def bsc_functions(beam, v, alpha, rho_0, phi_0, z_0, n_top):
    """g_tm(m, n) and g_te(m, n) of the beam at this position. Each value is computed
    once and kept, because the same BSC appears in several sums."""
    gTM, gTE = resolver_bsc(beam)
    lamb = getattr(beam, 'lamb', _LAMBDA_FW) if beam is not None else _LAMBDA_FW
    args = (v, alpha, lamb, rho_0, phi_0, 0.0, 0.0, z_0)
    cache_tm, cache_te = {}, {}

    def g_tm(m, n):
        if (m, n) not in cache_tm:
            cache_tm[m, n] = gTM(m, n, *args)
        return cache_tm[m, n]

    def g_te(m, n):
        if (m, n) not in cache_te:
            cache_te[m, n] = gTE(m, n, *args)
        return cache_te[m, n]

    g_tm(0, n_top)   # Frozen Wave: build its internal table for the largest n only once
    return g_tm, g_te


# ----------------------------------------------------------------------------- A_n
def A_n_x(n, M, x, g_tm):
    """Eq. (6.11), upper sign (J_x)."""
    cn = cs_n(M, MU_R, n, x)
    cn1 = cs_n(M, MU_R, n + 1, x)
    coef = np.conj(cn) * cn1 / (n + 1) ** 2

    total1 = total2 = total3 = total4 = 0j

    for m in range(0, n):                       # 1st sum: m = 0 .. n-1
        g1 = np.conj(g_tm(m + 1, n))
        g2 = g_tm(m, n + 1)
        total1 += g1 * g2 * fact_ratio(n + m + 1, n - m - 1)

    for m in range(0, n + 1):                   # 2nd sum: m = 0 .. n
        g1 = np.conj(g_tm(m, n))
        g2 = g_tm(m + 1, n + 1)
        total2 += g1 * g2 * fact_ratio(n + m + 2, n - m)

    for m in range(-n - 1, 0):                  # 3rd sum: m = -n-1 .. -1
        am = abs(m)
        g1 = np.conj(g_tm(m + 1, n))
        g2 = g_tm(m, n + 1)
        total3 += g1 * g2 * fact_ratio(n + am + 1, n - am + 1)

    for m in range(-n, 0):                      # 4th sum: m = -n .. -1
        am = abs(m)
        g1 = np.conj(g_tm(m, n))
        g2 = g_tm(m + 1, n + 1)
        total4 += g1 * g2 * fact_ratio(n + am, n - am)

    return coef * (-total1 + total2 + total3 - total4)


def A_n_y(n, M, x, g_tm):
    """Eq. (6.11), lower sign (J_y)."""
    cn = cs_n(M, MU_R, n, x)
    cn1 = cs_n(M, MU_R, n + 1, x)
    coef = np.conj(cn) * cn1 / (n + 1) ** 2

    total1 = total2 = total3 = total4 = 0j

    for m in range(0, n):                       # 1st sum: m = 0 .. n-1
        g1 = np.conj(g_tm(m + 1, n))
        g2 = g_tm(m, n + 1)
        total1 += g1 * g2 * fact_ratio(n + m + 1, n - m - 1)

    for m in range(0, n + 1):                   # 2nd sum: m = 0 .. n
        g1 = np.conj(g_tm(m, n))
        g2 = g_tm(m + 1, n + 1)
        total2 += g1 * g2 * fact_ratio(n + m + 2, n - m)

    for m in range(-n - 1, 0):                  # 3rd sum: m = -n-1 .. -1
        am = abs(m)
        g1 = np.conj(g_tm(m + 1, n))
        g2 = g_tm(m, n + 1)
        total3 += g1 * g2 * fact_ratio(n + am + 1, n - am + 1)

    for m in range(-n, 0):                      # 4th sum: m = -n .. -1
        am = abs(m)
        g1 = np.conj(g_tm(m, n))
        g2 = g_tm(m + 1, n + 1)
        total4 += g1 * g2 * fact_ratio(n + am, n - am)

    return coef * (total1 + total2 - total3 - total4)


# ----------------------------------------------------------------------------- B_n
def B_n_x(n, M, x, g_te, eta):
    """Eq. (6.12), upper sign (J_x)."""
    dn = ds_n(M, MU_R, n, x)
    dn1 = ds_n(M, MU_R, n + 1, x)
    coef = abs(eta) ** 2 * np.conj(dn) * dn1 / (n + 1) ** 2

    total1 = total2 = total3 = total4 = 0j

    for m in range(0, n):                       # 1st sum: m = 0 .. n-1
        g1 = np.conj(g_te(m + 1, n))
        g2 = g_te(m, n + 1)
        total1 += g1 * g2 * fact_ratio(n + m + 1, n - m - 1)

    for m in range(0, n + 1):                   # 2nd sum: m = 0 .. n
        g1 = np.conj(g_te(m, n))
        g2 = g_te(m + 1, n + 1)
        total2 += g1 * g2 * fact_ratio(n + m + 2, n - m)

    for m in range(-n - 1, 0):                  # 3rd sum: m = -n-1 .. -1
        am = abs(m)
        g1 = np.conj(g_te(m + 1, n))
        g2 = g_te(m, n + 1)
        total3 += g1 * g2 * fact_ratio(n + am + 1, n - am + 1)

    for m in range(-n, 0):                      # 4th sum: m = -n .. -1
        am = abs(m)
        g1 = np.conj(g_te(m, n))
        g2 = g_te(m + 1, n + 1)
        total4 += g1 * g2 * fact_ratio(n + am, n - am)

    return coef * (-total1 + total2 + total3 - total4)


def B_n_y(n, M, x, g_te, eta):
    """Eq. (6.12), lower sign (J_y)."""
    dn = ds_n(M, MU_R, n, x)
    dn1 = ds_n(M, MU_R, n + 1, x)
    coef = abs(eta) ** 2 * np.conj(dn) * dn1 / (n + 1) ** 2

    total1 = total2 = total3 = total4 = 0j

    for m in range(0, n):                       # 1st sum: m = 0 .. n-1
        g1 = np.conj(g_te(m + 1, n))
        g2 = g_te(m, n + 1)
        total1 += g1 * g2 * fact_ratio(n + m + 1, n - m - 1)

    for m in range(0, n + 1):                   # 2nd sum: m = 0 .. n
        g1 = np.conj(g_te(m, n))
        g2 = g_te(m + 1, n + 1)
        total2 += g1 * g2 * fact_ratio(n + m + 2, n - m)

    for m in range(-n - 1, 0):                  # 3rd sum: m = -n-1 .. -1
        am = abs(m)
        g1 = np.conj(g_te(m + 1, n))
        g2 = g_te(m, n + 1)
        total3 += g1 * g2 * fact_ratio(n + am + 1, n - am + 1)

    for m in range(-n, 0):                      # 4th sum: m = -n .. -1
        am = abs(m)
        g1 = np.conj(g_te(m, n))
        g2 = g_te(m + 1, n + 1)
        total4 += g1 * g2 * fact_ratio(n + am, n - am)

    return coef * (total1 + total2 - total3 - total4)


# ----------------------------------------------------------------------------- C_n
def C_n_x(n, M, x, g_tm, g_te, eta):
    """Eq. (6.13), upper sign (J_x)."""
    cn = cs_n(M, MU_R, n, x)
    dn = ds_n(M, MU_R, n, x)
    coef = (2 * n + 1) / (n * (n + 1)) ** 2 * np.conj(eta) * cn * np.conj(dn)

    total1 = total2 = 0j

    for m in range(0, n):                       # 1st sum: m = 0 .. n-1
        term1 = g_tm(m + 1, n) * np.conj(g_te(m, n))
        term2 = g_tm(m, n) * np.conj(g_te(m + 1, n))
        total1 += (term1 + term2) * fact_ratio(n + m + 1, n - m - 1)

    for m in range(-n, 0):                      # 2nd sum: m = -n .. -1
        am = abs(m)
        term1 = g_tm(m + 1, n) * np.conj(g_te(m, n))
        term2 = g_tm(m, n) * np.conj(g_te(m + 1, n))
        total2 += (term1 + term2) * fact_ratio(n + am, n - am)

    return coef * (total1 - total2)


def C_n_y(n, M, x, g_tm, g_te, eta):
    """Eq. (6.13), lower sign (J_y)."""
    cn = cs_n(M, MU_R, n, x)
    dn = ds_n(M, MU_R, n, x)
    coef = (2 * n + 1) / (n * (n + 1)) ** 2 * np.conj(eta) * cn * np.conj(dn)

    total1 = total2 = 0j

    for m in range(0, n):                       # 1st sum: m = 0 .. n-1
        term1 = g_tm(m + 1, n) * np.conj(g_te(m, n))
        term2 = g_tm(m, n) * np.conj(g_te(m + 1, n))
        total1 += (term1 - term2) * fact_ratio(n + m + 1, n - m - 1)

    for m in range(-n, 0):                      # 2nd sum: m = -n .. -1
        am = abs(m)
        term1 = g_tm(m + 1, n) * np.conj(g_te(m, n))
        term2 = g_tm(m, n) * np.conj(g_te(m + 1, n))
        total2 += (term1 - term2) * fact_ratio(n + am, n - am)

    return coef * (total1 - total2)


# ----------------------------------------------------------------------------- J_x, J_y
def I_xy(M, x, eps_imag, v, alpha, rho_0, z_0, phi_0=0.0, beam=None):
    """(J_x, J_y) with Pereira's expressions, Eq. (6.10).

    M: relative refractive index (n - i kappa); x: size parameter; eps_imag: eps'';
    v, alpha: beam order and half-cone angle (legacy signature);
    (rho_0, phi_0, z_0): position of the beam axis in the particle frame;
    beam: beam object with g_mnTM/g_mnTE (None = Frozen Wave of beamShape).
    """
    N = n_max(x)
    g_tm, g_te = bsc_functions(beam, v, alpha, rho_0, phi_0, z_0, N + 1)
    eta = Eta_r(M, 1)

    total_x = total_y = 0j
    for n in range(1, N + 1):
        Sn = S_n(M, n, x)
        radial_A = np.conj(Sn) + (n + 1) / M * r_n(M, n + 1, x)
        radial_B = -Sn + (n + 1) / M * r_n(M, n, x)

        total_x += (A_n_x(n, M, x, g_tm) * radial_A
                    + B_n_x(n, M, x, g_te, eta) * radial_B
                    - 1j * C_n_x(n, M, x, g_tm, g_te, eta) * Sn)
        total_y += (A_n_y(n, M, x, g_tm) * radial_A
                    + B_n_y(n, M, x, g_te, eta) * radial_B
                    - 1j * C_n_y(n, M, x, g_tm, g_te, eta) * Sn)

    prefactor = -3.0 * eps_imag / (abs(M) ** 2 * x ** 3)
    return prefactor * total_x.imag, prefactor * total_y.real


def I_x(M, x, eps_imag, v, alpha, rho_0, z_0, phi_0=0.0, beam=None):
    return I_xy(M, x, eps_imag, v, alpha, rho_0, z_0, phi_0, beam)[0]


def I_y(M, x, eps_imag, v, alpha, rho_0, z_0, phi_0=0.0, beam=None):
    return I_xy(M, x, eps_imag, v, alpha, rho_0, z_0, phi_0, beam)[1]


# ----------------------------------------------------------------------------- legacy helpers
# Old private interface, kept only for tests/comparacao_ABC_antigo_novo.py and
# tests/verificacao_forma_antiga.py: the four sums of A_n/B_n (times coef_n^* coef_{n+1},
# without 1/(n+1)^2) and the two cross sums of C_n, with callable BSCs bsc(m, n).
_nmax = n_max


def _AB(n, coef_n, coef_n1, bsc):
    cc = np.conj(coef_n) * coef_n1
    t1 = sum(np.conj(bsc(m + 1, n)) * bsc(m, n + 1) * fact_ratio(n + m + 1, n - m - 1) for m in range(0, n))
    t2 = sum(np.conj(bsc(m, n)) * bsc(m + 1, n + 1) * fact_ratio(n + m + 2, n - m) for m in range(0, n + 1))
    t3 = sum(np.conj(bsc(m + 1, n)) * bsc(m, n + 1) * fact_ratio(n - m + 1, n + m + 1) for m in range(-n - 1, 0))
    t4 = sum(np.conj(bsc(m, n)) * bsc(m + 1, n + 1) * fact_ratio(n - m, n + m) for m in range(-n, 0))
    return cc * t1, cc * t2, cc * t3, cc * t4


def _C(n, g_tm, g_te):
    pa = sum(g_tm(m + 1, n) * np.conj(g_te(m, n)) * fact_ratio(n + m + 1, n - m - 1) for m in range(0, n))
    pb = sum(g_tm(m, n) * np.conj(g_te(m + 1, n)) * fact_ratio(n + m + 1, n - m - 1) for m in range(0, n))
    qa = sum(g_tm(m + 1, n) * np.conj(g_te(m, n)) * fact_ratio(n - m, n + m) for m in range(-n, 0))
    qb = sum(g_tm(m, n) * np.conj(g_te(m + 1, n)) * fact_ratio(n - m, n + m) for m in range(-n, 0))
    return pa - qa, pb - qb
