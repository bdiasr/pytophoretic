import numpy as np
import scipy.special as special
from functools import lru_cache

from specials import *
from cte import *
from beamShape import g_mnTM, g_mnTE   # mantido por compatibilidade
from beams import resolver_bsc        # beam=None -> FW linear-x (beamShape)


def _fact_ratio(a: int, b: int) -> float:
    """a!/b! as float (fast, stable)."""
    return float(np.exp(special.gammaln(a + 1.0) - special.gammaln(b + 1.0)))


def cpw_n(k_in, n):
    a = 1 / (1j * k_in)
    b = (-1j) ** n
    c = (2 * n + 1) / (n * (n + 1))
    return a * b * c


@lru_cache(maxsize=None)
def cs_n(M_in, mi_in, n, x):
    """
    Cached coefficient c_n.

    M_in and mi_in are converted to builtin complex to keep hashing stable.
    """
    M0 = complex(M_in)
    mi0 = complex(mi_in)
    n0 = int(n)
    x0 = float(x)

    dem = (M0 * mi0 * ((RiccatiBessel_2Type(n0, x0) * derivativePsiBesseln(n0, x0)) - derivativeSphericalHankel_n(n0, x0) * psiBessel(n0, x0)))
    num = (mi0 * RiccatiBessel_2Type(n0, x0) * derivativePsiBesseln(n0, x0 * M0) - M0 * derivativeSphericalHankel_n(n0, x0) * psiBessel(n0, x0 * M0))
    return dem / num


@lru_cache(maxsize=None)
def ds_n(M_in, mi_in, n, x):
    """
    Cached coefficient d_n.
    """
    M0 = complex(M_in)
    mi0 = complex(mi_in)
    n0 = int(n)
    x0 = float(x)

    dem = (M0 * M0) * (RiccatiBessel_2Type(n0, x0) * derivativePsiBesseln(n0, x0) - derivativeSphericalHankel_n(n0, x0) * psiBessel(n0, x0))
    num = (M0 * RiccatiBessel_2Type(n0, x0) * derivativePsiBesseln(n0, M0 * x0)) - (mi0 * derivativeSphericalHankel_n(n0, x0) * psiBessel(n0, M0 * x0))
    return dem / num


@lru_cache(maxsize=None)
def r_n(m_in, n, x):
    m0 = complex(m_in)
    n0 = int(n)
    x0 = float(x)

    psiConj = np.conjugate(psiBessel(n0, m0 * x0))
    dem = (m0 * psiBessel(n0 + 1, m0 * x0) * psiConj)
    num = m0 * m0
    return (dem.imag) / (num.imag)


@lru_cache(maxsize=None)
def S_n(M_in, n, x):
    m0 = complex(M_in)
    n0 = int(n)
    x0 = float(x)

    a = (1j) / (2 * ((m0 ** 2).imag))
    b = (m0 * ((abs(psiBessel(n0, m0 * x0)) ** 2)))
    c = (np.conjugate(m0)) * (abs(psiBessel(n0 + 1, m0 * x0))) ** 2
    d = (m0 + (((2 * (n0 + 1) * ((m0 * m0).real) / m0)))) * r_n(m0, n0, x0)
    e = (2 * n0 + 1) * (np.conjugate(m0) * r_n(m0, n0 + 1, x0))
    return -a * (x0 * (b + c) - d + e)


def A_mn(n, M_in, x, v, alpha1, rho_0, phi_0, x_0, y_0, z_0, beam=None):
    gTM, gTE = resolver_bsc(beam)
    total1 = 0.0 + 0.0j
    total2 = 0.0 + 0.0j
    total3 = 0.0 + 0.0j
    total4 = 0.0 + 0.0j

    a = 1.0 / ((n + 1) ** 2)

    cn = cs_n(M_in, mi, n, x)
    cn1 = cs_n(M_in, mi, n + 1, x)
    cn_conj = np.conjugate(cn)
    cn1_conj = np.conjugate(cn1)

    for m in range(0, n):
        g1 = gTM(m + 1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        g2 = np.conjugate(gTM(m, n + 1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        total1 += (cn * cn1_conj * g1 * g2 * _fact_ratio(n + m + 1, n - m - 1))

    for m in range(0, n + 1):
        g1 = np.conjugate(gTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        g2 = gTM(m + 1, n + 1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        total2 += (cn_conj * cn1 * g1 * g2 * _fact_ratio(n + m + 2, n - m))

    for m in range(-n - 1, 0):
        am = abs(m)
        g1 = gTM(m + 1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        g2 = np.conjugate(gTM(m, n + 1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        total3 += (cn * cn1_conj * g1 * g2 * _fact_ratio(n + am + 1, n - am + 1))

    for m in range(-n, 0):
        am = abs(m)
        g1 = np.conjugate(gTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        g2 = gTM(m + 1, n + 1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        total4 += (cn_conj * cn1 * g1 * g2 * _fact_ratio(n + am, n - am))

    return a * (total1 + total2 - total3 - total4)


def B_mn(n, M_in, x, v, alpha1, rho_0, phi_0, x_0, y_0, z_0, beam=None):
    gTM, gTE = resolver_bsc(beam)
    total1 = 0.0 + 0.0j
    total2 = 0.0 + 0.0j
    total3 = 0.0 + 0.0j
    total4 = 0.0 + 0.0j

    eta_r = Eta_r(M_in, 1)
    const = (abs(eta_r) ** 2) / ((n + 1) ** 2)

    dn = ds_n(M_in, mi, n, x)
    dn1 = ds_n(M_in, mi, n + 1, x)
    dn_conj = np.conjugate(dn)
    dn1_conj = np.conjugate(dn1)

    for m in range(0, n):
        g1 = gTE(m + 1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        g2 = np.conjugate(gTE(m, n + 1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        total1 += (dn * dn1_conj * g1 * g2 * _fact_ratio(n + m + 1, n - m - 1))

    for m in range(0, n + 1):
        g1 = np.conjugate(gTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        g2 = gTE(m + 1, n + 1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        total2 += (dn_conj * dn1 * g1 * g2 * _fact_ratio(n + m + 2, n - m))

    for m in range(-n - 1, 0):
        am = abs(m)
        g1 = gTE(m + 1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        g2 = np.conjugate(gTE(m, n + 1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        total3 += (dn * dn1_conj * g1 * g2 * _fact_ratio(n + am + 1, n - am + 1))

    for m in range(-n, 0):
        am = abs(m)
        g1 = np.conjugate(gTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
        g2 = gTE(m + 1, n + 1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        total4 += (dn_conj * dn1 * g1 * g2 * _fact_ratio(n + am, n - am))

    return const * (total1 + total2 - total3 - total4)


def C_mnx(n, M_in, x, v, alpha1, rho_0, phi_0, x_0, y_0, z_0, beam=None):
    gTM, gTE = resolver_bsc(beam)
    total1 = 0.0 + 0.0j
    total2 = 0.0 + 0.0j

    eta_r = Eta_r(M_in, 1)
    frac = ((2 * n) + 1) / ((n ** 2) * ((n + 1) ** 2))
    conj_eta = np.conjugate(eta_r)

    cn = cs_n(M_in, mi, n, x)
    dn = np.conjugate(ds_n(M_in, mi, n, x))

    for m in range(0, n + 1):
        gTM_1 = gTM(m + 1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        gTE_1 = np.conjugate(gTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))

        gTM_2 = gTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        gTE_2 = np.conjugate(gTE(m + 1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))

        total1 += (gTM_1 * gTE_1 + gTM_2 * gTE_2) * _fact_ratio(n + m + 1, n - m - 1)

    for m in range(-n, 0):
        am = abs(m)
        gTM_1 = gTM(m + 1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        gTE_1 = np.conjugate(gTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))

        gTM_2 = gTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        gTE_2 = np.conjugate(gTE(m + 1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))

        total2 += (gTM_1 * gTE_1 + gTM_2 * gTE_2) * _fact_ratio(n + am, n - am)

    return frac * conj_eta * cn * dn * (total1 - total2)


def C_mny(n, M_in, x, v, alpha1, rho_0, phi_0, x_0, y_0, z_0, beam=None):
    gTM, gTE = resolver_bsc(beam)
    total1 = 0.0 + 0.0j
    total2 = 0.0 + 0.0j

    eta_r = Eta_r(M_in, 1)
    frac = (2 * n + 1) / ((n * (n + 1)) ** 2)
    conj_eta = np.conjugate(eta_r)

    cn = cs_n(M_in, mi, n, x)
    dn = np.conjugate(ds_n(M_in, mi, n, x))

    for m in range(0, n + 1):
        gTM_1 = gTM(m + 1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        gTE_1 = np.conjugate(gTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))

        gTM_2 = gTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        gTE_2 = np.conjugate(gTE(m + 1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))

        total1 += (gTM_1 * gTE_1 - gTM_2 * gTE_2) * _fact_ratio(n + m + 1, n - m - 1)

    for m in range(-n, 0):
        am = abs(m)
        gTM_1 = gTM(m + 1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        gTE_1 = np.conjugate(gTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))

        gTM_2 = gTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
        gTE_2 = np.conjugate(gTE(m + 1, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))

        total2 += (gTM_1 * gTE_1 - gTM_2 * gTE_2) * _fact_ratio(n + am, n - am)

    return frac * conj_eta * cn * dn * (total1 - total2)


def D_mn(m, n, M_in, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0, beam=None):
    gTM, gTE = resolver_bsc(beam)
    alpha1 = alpha
    a1 = 1.0 / ((n + 1) ** 2)

    cn = cs_n(M_in, mi, n, x)
    b = np.conjugate(cn) * cs_n(M_in, mi, n + 1, x)

    c = np.conjugate(gTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
    d = gTM(m, n + 1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)

    am = abs(int(m))
    e = _fact_ratio(n + am + 1, n - am)

    return a1 * b * c * d * e


def E_mn(m, n, M_in, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0, beam=None):
    gTM, gTE = resolver_bsc(beam)
    alpha1 = alpha
    eta_r = Eta_r(M_in, 1)

    a1 = (abs(eta_r) ** 2) / ((n + 1) ** 2)
    b = np.conjugate(ds_n(M_in, mi, n, x)) * ds_n(M_in, mi, n + 1, x)

    c = np.conjugate(gTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))
    d = gTE(m, n + 1, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)

    am = abs(int(m))
    e = _fact_ratio(n + am + 1, n - am)

    return a1 * b * c * d * e


def F_mn(m, n, M_in, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0, beam=None):
    gTM, gTE = resolver_bsc(beam)
    alpha1 = alpha
    eta_r = Eta_r(M_in, 1)

    a1 = (m * np.conjugate(eta_r)) * ((2 * n + 1) / ((n * (n + 1)) ** 2))
    b = cs_n(M_in, mi, n, x) * np.conjugate(ds_n(M_in, mi, n, x))

    c = gTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
    d = np.conjugate(gTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0))

    am = abs(int(m))
    e = _fact_ratio(n + am, n - am)

    return a1 * b * c * d * e
