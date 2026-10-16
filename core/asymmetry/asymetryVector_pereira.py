# -*- coding: utf-8 -*-
"""
Vetor de assimetria com as componentes TRANSVERSAIS corrigidas por D. Pereira
=============================================================================
Referencias:
  D. Pereira, "Expressoes Jx e Jy do artigo 'Towards Photophoresis...'", nota
  tecnica (out. 2025), Eqs. (8)-(11);
  D. J. S. Pereira, "3D Photophoresis in Slip-flow", tese, Univ. Coimbra (2026),
  Eqs. (6.10)-(6.13) e Apendice H.

Modulo PARALELO: asymetryVector.py e coefficients.py nao sao alterados.
Mesma assinatura e mesma normalizacao de asymetryVector (retorna J_i, com
J_i = 3 eps'' I_i / (8 pi a^3)):

    I_x(M, x, epislon, v, alpha, rho_0, z_0, phi_0=0.0, beam=None)
    I_y(...), I_z(...)  -- I_z e' reexportado de asymetryVector (inalterado).

Forma implementada (Pereira, Eq. 8):
    I_{x;y} = -(8 pi a^3 / |M|^2 x^3) {Im;Re} sum_n [ A_n (S_n^* + (n+1)R_{n+1}/M)
                                                 + B_n (-S_n + (n+1)R_n/M) - i C_n S_n ]
    A_n = 1/(n+1)^2 { -/+ sum_{m=0}^{n-1}   c_n^* c_{n+1} g^{m+1*}_{n,TM} g^m_{n+1,TM} (n+m+1)!/(n-m-1)!
                      +   sum_{m=0}^{n}     c_n^* c_{n+1} g^{m*}_{n,TM}   g^{m+1}_{n+1,TM} (n+m+2)!/(n-m)!
                      +/- sum_{m=-n-1}^{-1} c_n^* c_{n+1} g^{m+1*}_{n,TM} g^m_{n+1,TM} (n+|m|+1)!/(n-|m|+1)!
                      -   sum_{m=-n}^{-1}   c_n^* c_{n+1} g^{m*}_{n,TM}   g^{m+1}_{n+1,TM} (n+|m|)!/(n-|m|)! }
    (sinal superior: x; inferior: y).  B_n: idem com d_n, TE e prefator |eta_r|^2/(n+1)^2.
    C_n = (2n+1)/[n(n+1)]^2 eta_r^* c_n d_n^* { sum_{m=0}^{n-1} (g^{m+1}_TM g^{m*}_TE +/- g^m_TM g^{m+1*}_TE) (n+m+1)!/(n-m-1)!
                                              - sum_{m=-n}^{-1}  (idem) (n+|m|)!/(n-|m|)! }
Verificacao: tests/diagnostico_pereira2025.py (R = 1 +/- 3e-11 contra a integral
de definicao, casos B1-B4 de tests/verificacao_ras.py).
"""
import math
import numpy as np

from coefficients import cs_n, ds_n, S_n, r_n, _fact_ratio
from cte import Eta_r, lamb as _LAMB_FW
from beams import resolver_bsc
from asymetryVector import I_z  # noqa: F401  (componente longitudinal: inalterada)

MI = 1


def _fr(a, b):
    return _fact_ratio(a, b) if b >= 0 else 0.0


def _nmax(x):
    return math.ceil(x + 4.05 * x ** (1 / 3) + 2)


def _AB(n, s_n, s_n1, G):
    """Retorna (termos 1 e 3, termos 2 e 4) de A_n/B_n SEM prefator e SEM sinal
    dependente da componente: t13 = t1 - t3... separados para montar x e y."""
    cc = np.conj(s_n) * s_n1
    t1 = t2 = t3 = t4 = 0j
    for m in range(0, n):
        t1 += np.conj(G(m + 1, n)) * G(m, n + 1) * _fr(n + m + 1, n - m - 1)
    for m in range(0, n + 1):
        t2 += np.conj(G(m, n)) * G(m + 1, n + 1) * _fr(n + m + 2, n - m)
    for m in range(-n - 1, 0):
        am = -m
        t3 += np.conj(G(m + 1, n)) * G(m, n + 1) * _fr(n + am + 1, n - am + 1)
    for m in range(-n, 0):
        am = -m
        t4 += np.conj(G(m, n)) * G(m + 1, n + 1) * _fr(n + am, n - am)
    return cc * t1, cc * t2, cc * t3, cc * t4


def _C(n, GTM, GTE):
    """(soma com +, soma com -) dos dois produtos de C_n, sem prefator."""
    pa = pb = qa = qb = 0j
    for m in range(0, n):
        f = _fr(n + m + 1, n - m - 1)
        pa += GTM(m + 1, n) * np.conj(GTE(m, n)) * f
        pb += GTM(m, n) * np.conj(GTE(m + 1, n)) * f
    for m in range(-n, 0):
        f = _fr(n - m, n + m)
        qa += GTM(m + 1, n) * np.conj(GTE(m, n)) * f
        qb += GTM(m, n) * np.conj(GTE(m + 1, n)) * f
    return (pa - qa), (pb - qb)          # C ~ (pa-qa) +/- (pb-qb)


def I_xy(M, x, epislon, v, alpha, rho_0, z_0, phi_0=0.0, beam=None):
    """(J_x, J_y) pela correcao de Pereira. Mesma convencao de asymetryVector."""
    gTM, gTE = resolver_bsc(beam)
    lamb = getattr(beam, 'lamb', _LAMB_FW) if beam is not None else _LAMB_FW
    args = (v, alpha, lamb, rho_0, phi_0, 0.0, 0.0, z_0)
    GTM = lambda m, nn: gTM(m, nn, *args)
    GTE = lambda m, nn: gTE(m, nn, *args)
    eta = Eta_r(M, 1)
    Zx = Zy = 0j
    for n in range(1, _nmax(x) + 1):
        cn, cn1 = cs_n(M, MI, n, x), cs_n(M, MI, n + 1, x)
        dn, dn1 = ds_n(M, MI, n, x), ds_n(M, MI, n + 1, x)
        a1, a2, a3, a4 = _AB(n, cn, cn1, GTM)
        b1, b2, b3, b4 = _AB(n, dn, dn1, GTE)
        pA = 1.0 / (n + 1) ** 2
        pB = abs(eta) ** 2 / (n + 1) ** 2
        Ax = pA * (-a1 + a2 + a3 - a4); Ay = pA * (a1 + a2 - a3 - a4)
        Bx = pB * (-b1 + b2 + b3 - b4); By = pB * (b1 + b2 - b3 - b4)
        ca, cb = _C(n, GTM, GTE)
        pC = (2 * n + 1) / (n * (n + 1)) ** 2 * np.conj(eta) * cn * np.conj(dn)
        Cx = pC * (ca + cb); Cy = pC * (ca - cb)
        Sn = S_n(M, n, x)
        FA = np.conj(Sn) + (n + 1) / M * r_n(M, n + 1, x)
        FB = -Sn + (n + 1) / M * r_n(M, n, x)
        Zx += Ax * FA + Bx * FB - 1j * Cx * Sn
        Zy += Ay * FA + By * FB - 1j * Cy * Sn
    const = -3.0 * epislon / (abs(M) ** 2 * x ** 3)
    return const * Zx.imag, const * Zy.real


def I_x(M, x, epislon, v, alpha, rho_0, z_0, phi_0=0.0, beam=None):
    return I_xy(M, x, epislon, v, alpha, rho_0, z_0, phi_0, beam)[0]


def I_y(M, x, epislon, v, alpha, rho_0, z_0, phi_0=0.0, beam=None):
    return I_xy(M, x, epislon, v, alpha, rho_0, z_0, phi_0, beam)[1]
