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
                                                   - i C_n S_n ]

    A_n = c_n^* c_{n+1} / (n+1)^2 * ( -/+ S1 + S2 +/- S3 - S4 )   (TM BSCs)
    B_n = |eta_r|^2 d_n^* d_{n+1} / (n+1)^2 * ( same with TE BSCs )
    C_n = (2n+1) / [n(n+1)]^2 * eta_r^* c_n d_n^* * ( C_plus +/- C_minus )

with the partial sums S1..S4, C_plus, C_minus defined in _ab_partial_sums and
_c_partial_sums. Upper sign: x; lower sign: y.

Verification: tests/diagnostico_pereira2025.py (ratio 1 +/- 3e-11 against the
defining volume integral).
"""
import math
import numpy as np

from coefficients import cs_n, ds_n, S_n, r_n, _fact_ratio
from cte import Eta_r, lamb as _LAMBDA_FW
from beams import resolver_bsc
from asymetryVector import I_z  # noqa: F401  (longitudinal component, unchanged)

MU_R = 1   # relative permeability used in c_n, d_n


def n_max(x):
    """Truncation of the multipole series (Wiscombe)."""
    return math.ceil(x + 4.05 * x ** (1 / 3) + 2)


def _factorial_ratio(a, b):
    """a!/b!; zero when b < 0 (the term does not exist)."""
    return _fact_ratio(a, b) if b >= 0 else 0.0


def _ab_partial_sums(n, coef_n, coef_n1, bsc):
    """Partial sums of A_n (coef = c, bsc = TM) or B_n (coef = d, bsc = TE),
    already multiplied by coef_n^* coef_{n+1} but without the 1/(n+1)^2 prefactor.

        S1 = sum_{m=0}^{n-1}    g^{m+1*}_n g^m_{n+1}     (n+m+1)!/(n-m-1)!
        S2 = sum_{m=0}^{n}      g^{m*}_n   g^{m+1}_{n+1} (n+m+2)!/(n-m)!
        S3 = sum_{m=-n-1}^{-1}  g^{m+1*}_n g^m_{n+1}     (n+|m|+1)!/(n-|m|+1)!
        S4 = sum_{m=-n}^{-1}    g^{m*}_n   g^{m+1}_{n+1} (n+|m|)!/(n-|m|)!
    """
    s1 = s2 = s3 = s4 = 0j
    for m in range(0, n):
        s1 += np.conj(bsc(m + 1, n)) * bsc(m, n + 1) * _factorial_ratio(n + m + 1, n - m - 1)
    for m in range(0, n + 1):
        s2 += np.conj(bsc(m, n)) * bsc(m + 1, n + 1) * _factorial_ratio(n + m + 2, n - m)
    for m in range(-n - 1, 0):
        s3 += np.conj(bsc(m + 1, n)) * bsc(m, n + 1) * _factorial_ratio(n - m + 1, n + m + 1)
    for m in range(-n, 0):
        s4 += np.conj(bsc(m, n)) * bsc(m + 1, n + 1) * _factorial_ratio(n - m, n + m)
    coef_product = np.conj(coef_n) * coef_n1
    return coef_product * s1, coef_product * s2, coef_product * s3, coef_product * s4


def _c_partial_sums(n, g_tm, g_te):
    """The two TM-TE cross sums of C_n, without prefactor:

        C_plus  = sum_{m>=0} g^{m+1}_TM g^{m*}_TE (n+m+1)!/(n-m-1)!  -  sum_{m<0} (same) (n+|m|)!/(n-|m|)!
        C_minus = sum_{m>=0} g^m_TM g^{m+1*}_TE   (n+m+1)!/(n-m-1)!  -  sum_{m<0} (same) (n+|m|)!/(n-|m|)!

    C_n is proportional to C_plus + C_minus (x) or C_plus - C_minus (y).
    """
    plus_pos = minus_pos = plus_neg = minus_neg = 0j
    for m in range(0, n):
        weight = _factorial_ratio(n + m + 1, n - m - 1)
        plus_pos += g_tm(m + 1, n) * np.conj(g_te(m, n)) * weight
        minus_pos += g_tm(m, n) * np.conj(g_te(m + 1, n)) * weight
    for m in range(-n, 0):
        weight = _factorial_ratio(n - m, n + m)
        plus_neg += g_tm(m + 1, n) * np.conj(g_te(m, n)) * weight
        minus_neg += g_tm(m, n) * np.conj(g_te(m + 1, n)) * weight
    return plus_pos - plus_neg, minus_pos - minus_neg


def I_xy(M, x, eps_imag, v, alpha, rho_0, z_0, phi_0=0.0, beam=None):
    """(J_x, J_y) with Pereira's expressions.

    M: relative refractive index (n - i kappa); x: size parameter; eps_imag: eps'';
    v, alpha: beam order and half-cone angle (legacy signature);
    (rho_0, phi_0, z_0): position of the beam axis in the particle frame;
    beam: beam object with g_mnTM/g_mnTE (None = Frozen Wave of beamShape).
    """
    g_tm_func, g_te_func = resolver_bsc(beam)
    lamb = getattr(beam, 'lamb', _LAMBDA_FW) if beam is not None else _LAMBDA_FW
    position = (v, alpha, lamb, rho_0, phi_0, 0.0, 0.0, z_0)

    def g_tm(m, n):
        return g_tm_func(m, n, *position)

    def g_te(m, n):
        return g_te_func(m, n, *position)

    eta = Eta_r(M, 1)
    sum_x = sum_y = 0j
    for n in range(1, n_max(x) + 1):
        c_n, c_n1 = cs_n(M, MU_R, n, x), cs_n(M, MU_R, n + 1, x)
        d_n, d_n1 = ds_n(M, MU_R, n, x), ds_n(M, MU_R, n + 1, x)

        # A_n (TM) and B_n (TE): same partial sums, combined with x/y-dependent signs
        a1, a2, a3, a4 = _ab_partial_sums(n, c_n, c_n1, g_tm)
        b1, b2, b3, b4 = _ab_partial_sums(n, d_n, d_n1, g_te)
        pref_a = 1.0 / (n + 1) ** 2
        pref_b = abs(eta) ** 2 / (n + 1) ** 2
        A_x = pref_a * (-a1 + a2 + a3 - a4)
        A_y = pref_a * (a1 + a2 - a3 - a4)
        B_x = pref_b * (-b1 + b2 + b3 - b4)
        B_y = pref_b * (b1 + b2 - b3 - b4)

        # C_n: TM-TE cross term
        c_plus, c_minus = _c_partial_sums(n, g_tm, g_te)
        pref_c = (2 * n + 1) / (n * (n + 1)) ** 2 * np.conj(eta) * c_n * np.conj(d_n)
        C_x = pref_c * (c_plus + c_minus)
        C_y = pref_c * (c_plus - c_minus)

        # radial integrals of the internal field (the same factors as in J_z)
        s_n = S_n(M, n, x)
        radial_a = np.conj(s_n) + (n + 1) / M * r_n(M, n + 1, x)
        radial_b = -s_n + (n + 1) / M * r_n(M, n, x)

        sum_x += A_x * radial_a + B_x * radial_b - 1j * C_x * s_n
        sum_y += A_y * radial_a + B_y * radial_b - 1j * C_y * s_n

    prefactor = -3.0 * eps_imag / (abs(M) ** 2 * x ** 3)
    return prefactor * sum_x.imag, prefactor * sum_y.real


def I_x(M, x, eps_imag, v, alpha, rho_0, z_0, phi_0=0.0, beam=None):
    return I_xy(M, x, eps_imag, v, alpha, rho_0, z_0, phi_0, beam)[0]


def I_y(M, x, eps_imag, v, alpha, rho_0, z_0, phi_0=0.0, beam=None):
    return I_xy(M, x, eps_imag, v, alpha, rho_0, z_0, phi_0, beam)[1]


# Old private names, still used by tests/comparacao_ABC_antigo_novo.py and
# tests/verificacao_forma_antiga.py.
_nmax = n_max
_fr = _factorial_ratio
_AB = _ab_partial_sums
_C = _c_partial_sums
