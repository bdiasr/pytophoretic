# -*- coding: utf-8 -*-
"""
Forca fotoforetica -- Modelo de Loesche-Husmann Modificado GENERALIZADO
======================================================================
Referencia unica: D. J. S. Pereira, "3D Photophoresis in Slip-flow", tese de
doutorado, Univ. Coimbra (2026):
  Eq. (6.14)  vetor forca fotoforetica (feixes arbitrarios);
  Eq. (5.3)   T~ = T_pS (temperatura media da superficie da particula);
  Eq. (5.5)   k_g medio (Zhang & Bar-Ziv);
  Eq. (5.6)   T_gS (temperatura media do gas na superficie);
  Eq. (5.10)  f_pc(T*), T* = T_gS/T_inf  (valida para ar atmosferico);
  Eq. (6.15)  J0 I0 = C_abs I_lambda / (4 pi Rp^2);
  Eqs. (6.16)-(6.17) C_ext, C_sca da GLMT (Gouesbet & Grehan);
  Tabela 2.1  coeficientes de deslizamento (acomodacao completa);
  Apendice E  propriedades do ar (Eqs. E.1, E.2, Tabela E.1).

  F_ph = -4 pi f_pc(T*) mu_g(T_inf) [ mu_g(T_gS)/(rho_g(T_gS) T_gS)
                                      (K_S + 3 K_h K_m Kn)/(1 + 3 K_m Kn) ]
         x ( 1/(1 + 2 K_T Kn) * I_lambda J1 /
             ( k_p/Rp + (k_g_bar/Rp) 2/(1 + 2 K_T Kn) + 4 sigma eps T~^3 ) )

J1 = (J_x, J_y, J_z): vetor de assimetria com a normalizacao da Eq. (6.8) --
a mesma de asymetryVector / asymetryVector_pereira (J_x, J_y corrigidos por
Pereira; J_z inalterado).

PREFATOR 1/(1 + 2 K_T Kn) -- as Eqs. (5.7)/(6.14)
impressas trazem 1/(2 + 2 K_T Kn), mas com ela as forcas dos casos de teste
1-3 (Figs. 6.8, 6.10, 6.11) saem 1,894x menores que as da tese, enquanto
T_gS bate exatamente. Com 1/(1 + 2 K_T Kn) (prefator do L-H original, Eq. 2.12)
Fz = 0,404/43,45/-0,169 e |Fx| = 0,165/14,58/0,0353 vs tese 0,4/43/-0,17 e
0,16/14,5/0,035. Adota-se 1/(1 + 2 K_T Kn); a forma impressa fica disponivel
via prefator_KT='impresso_6.14' (apenas para comparacao).

GRANDEZAS QUE A TESE NAO FIXA NUMERICAMENTE e que, por isso, sao ARGUMENTOS
OBRIGATORIOS (sem valor padrao): T_inf, T_sur, livre caminho medio ell
(Kn = ell/Rp, Eq. 1.1) e aceleracao da gravidade (so para F/mg).

Convencao de posicao: (rho_0, phi_0) = posicao do EIXO DO FEIXE no referencial
da particula (a do projeto). Para a particula em x_P sobre Ox (Fig. 6.7 da
tese): rho_0 = |x_P|, phi_0 = pi se x_P > 0.

Modulo novo; nao altera nenhum outro arquivo.
"""
import math
import numpy as np
from scipy.integrate import quad

from specials import psiBessel, derivativePsiBesseln, RiccatiBessel_2Type, \
    derivativeSphericalHankel_n
from coefficients import _fact_ratio
from beams import resolver_bsc
from cte import epsilonR2
import asymetryVector_pereira as AVP

SIGMA_SB = 5.670374419e-8  # W m^-2 K^-4 (CODATA)

# ---------------------------------------------------------------- Apendice E
_POLI = {  # a4, a3, a2, a1, a0  (Tabela E.1)
    'mu': (-3.967e-18, 2.1509e-14, -4.7024e-11, 6.9898e-8, 1.0925e-6),
    'cp': (1.1229e-10, -5.3562e-7, 8.2717e-4, -2.9542e-1, 1.0321e3),
    'k':  (-3.2253e-15, 1.9649e-11, -5.1570e-8, 1.0334e-4, -6.4045e-4),
}


def _poli(nome, T):
    a4, a3, a2, a1, a0 = _POLI[nome]
    return (((a4 * T + a3) * T + a2) * T + a1) * T + a0


def rho_ar(T):
    """Eq. (E.1): rho = e^{5.86697} T^{-1}  [kg/m^3], T em K."""
    return math.exp(5.86697) / T


def mu_ar(T):
    """Eq. (E.2), Tabela E.1: viscosidade dinamica [Pa s]."""
    return _poli('mu', T)


def cp_ar(T):
    return _poli('cp', T)


def k_ar(T):
    """Eq. (E.2), Tabela E.1: condutividade termica [W/(m K)]."""
    return _poli('k', T)


# ---------------------------------------------------------------- Tabela 2.1
TABELA_2_1 = dict(K_S=1.152, K_T=2.179, K_m=1.131)   # K_h: 1-3 (ajustavel)


def f_pc(Tstar):
    """Eq. (5.10) -- valida para ar atmosferico."""
    return -0.0171 * Tstar ** 2 + 0.2851 * Tstar + 0.7327


def k_g_medio(T_inf, T_gS):
    """Eq. (5.5)."""
    if abs(T_gS - T_inf) < 1e-9:
        return k_ar(T_inf)
    return quad(k_ar, T_inf, T_gS)[0] / (T_gS - T_inf)


# ------------------------------------------------- GLMT: a_n, b_n, C_ext, C_sca
def ab_n(M, n, x):
    """Coeficientes de espalhamento externos (GLMT, exp(+i w t), mu_r = 1):
    a_n = [psi(x) psi'(Mx) - M psi'(x) psi(Mx)] / [xi(x) psi'(Mx) - M xi'(x) psi(Mx)]
    b_n = [M psi(x) psi'(Mx) - psi'(x) psi(Mx)] / [M xi(x) psi'(Mx) - xi'(x) psi(Mx)]
    xi = x h^(2)."""
    ps, dps = psiBessel(n, x), derivativePsiBesseln(n, x)
    psM, dpsM = psiBessel(n, M * x), derivativePsiBesseln(n, M * x)
    xi, dxi = RiccatiBessel_2Type(n, x), derivativeSphericalHankel_n(n, x)
    a = (ps * dpsM - M * dps * psM) / (xi * dpsM - M * dxi * psM)
    b = (M * ps * dpsM - dps * psM) / (M * xi * dpsM - dxi * psM)
    return a, b


def _nmax(x):
    return math.ceil(x + 4.05 * x ** (1 / 3) + 2)


def secoes_choque(M, x, beam, v, rho_0, phi_0, z_0, lamb, n_max=None):
    """Eqs. (6.16)-(6.17). Retorna (C_ext, C_sca, C_abs) em m^2."""
    gTM, gTE = resolver_bsc(beam)
    args = (v, 0.0, lamb, rho_0, phi_0, 0.0, 0.0, z_0)
    se = ss = 0.0
    for n in range(1, (n_max or _nmax(x)) + 1):
        a, b = ab_n(M, n, x)
        pref = (2 * n + 1) / (n * (n + 1))
        for m in range(-n, n + 1):
            f = pref * _fact_ratio(n + abs(m), n - abs(m))
            t2 = abs(gTM(m, n, *args)) ** 2
            e2 = abs(gTE(m, n, *args)) ** 2
            se += f * (a * t2 + b * e2).real
            ss += f * (abs(a) ** 2 * t2 + abs(b) ** 2 * e2)
    c = lamb ** 2 / math.pi
    return c * se, c * ss, c * (se - ss)


# ------------------------------------------------------------ temperaturas
def temperaturas(I0J0, Rp, emiss, Kn, K_T, T_inf, T_sur, tol=1e-10, itmax=500):
    """Eqs. (5.3), (5.6), (5.5) em laco, iniciando com k_g_bar = k_g(T_inf).
    Retorna dict(T_pS=T~, T_gS, k_g_bar, iter)."""
    kb = k_ar(T_inf)
    Tt = T_inf
    TgS = T_inf
    s = SIGMA_SB * emiss
    for it in range(1, itmax + 1):
        a1 = 1.0 / (1.0 + K_T * Kn)
        # Eq. (5.3): ponto fixo em T~ (T~ nos dois lados)
        for _ in range(200):
            num = I0J0 + a1 * kb / Rp * T_inf + s * (3 * Tt ** 4 + T_sur ** 4)
            den = kb / Rp * a1 + 4 * s * Tt ** 3
            Tn = num / den
            if abs(Tn - Tt) < tol * Tt:
                Tt = Tn
                break
            Tt = Tn
        # Eq. (5.6)
        num = I0J0 - 4 * s * Tt ** 3 * T_inf + s * (3 * Tt ** 4 + T_sur ** 4)
        TgS_n = T_inf + a1 * num / (kb / Rp * a1 + 4 * s * Tt ** 3)
        kb_n = k_g_medio(T_inf, TgS_n)                    # Eq. (5.5)
        conv = abs(TgS_n - TgS) < tol * TgS_n and abs(kb_n - kb) < tol * kb_n
        TgS, kb = TgS_n, kb_n
        if conv:
            break
    return dict(T_pS=Tt, T_gS=TgS, k_g_bar=kb, iter=it)


# ------------------------------------------------------------- forca (6.14)
def forca(J1, I_lambda, C_abs, Rp, k_p, emiss, ell, T_inf, T_sur,
          K_h, coef=TABELA_2_1, tol=1e-10, prefator_KT='LH_2.12', rho_g_fixo=None):
    """Eq. (6.14). J1 = (Jx, Jy, Jz) [adim.], I_lambda [W/m^2], C_abs [m^2].
    Retorna (F [N, vetor 3], dict de temperaturas/intermediarios)."""
    K_S, K_T, K_m = coef['K_S'], coef['K_T'], coef['K_m']
    Kn = ell / Rp                                          # Eq. (1.1)
    I0J0 = C_abs * I_lambda / (4 * math.pi * Rp ** 2)      # Eq. (6.15)
    T = temperaturas(I0J0, Rp, emiss, Kn, K_T, T_inf, T_sur, tol)
    TgS, Tt, kb = T['T_gS'], T['T_pS'], T['k_g_bar']
    fpc = f_pc(TgS / T_inf)
    rho = rho_ar(TgS) if rho_g_fixo is None else float(rho_g_fixo)   # rho_g_fixo: p.ex. Tab. 1 de Sun 2022 (30 Torr)
    colch = mu_ar(TgS) / (rho * TgS) * (K_S + 3 * K_h * K_m * Kn) / (1 + 3 * K_m * Kn)
    den = k_p / Rp + kb / Rp * 2.0 / (1 + 2 * K_T * Kn) + 4 * SIGMA_SB * emiss * Tt ** 3
    if prefator_KT == 'LH_2.12':
        dKT = 1 + 2 * K_T * Kn          # adotado (ver cabecalho)
    elif prefator_KT == 'impresso_6.14':
        dKT = 2 + 2 * K_T * Kn          # forma impressa na tese
    else:
        raise ValueError("prefator_KT: 'LH_2.12' ou 'impresso_6.14'")
    fator = -4 * math.pi * fpc * mu_ar(T_inf) * colch / dKT * I_lambda / den
    F = fator * np.asarray(J1, float)
    T.update(Kn=Kn, I0J0=I0J0, f_pc=fpc, fator=fator)
    return F, T


def forca_feixe(M, Rp, beam, rho_0, phi_0, z_0, I_lambda, k_p, emiss, ell,
                T_inf, T_sur, K_h, coef=TABELA_2_1, prefator_KT='LH_2.12', rho_g_fixo=None):
    """Encadeia: J1 (asymetryVector_pereira) + C_abs (6.16-6.17) + Eq. (6.14).
    beam: objeto de feixe do projeto (com .lamb e .v)."""
    lamb, v = beam.lamb, getattr(beam, 'v', 0)
    x = 2 * math.pi / lamb * Rp
    eps2 = epsilonR2(M, 1)
    jx, jy = AVP.I_xy(M, x, eps2, v, 0.0, rho_0, z_0, phi_0, beam=beam)
    jz = AVP.I_z(M, x, eps2, v, 0.0, rho_0, z_0, phi_0, beam=beam)
    Ce, Cs, Ca = secoes_choque(M, x, beam, v, rho_0, phi_0, z_0, lamb)
    F, info = forca((jx, jy, jz), I_lambda, Ca, Rp, k_p, emiss, ell, T_inf, T_sur, K_h, coef,
                    prefator_KT=prefator_KT, rho_g_fixo=rho_g_fixo)
    info.update(J=(jx, jy, jz), C_ext=Ce, C_sca=Cs, C_abs=Ca)
    return F, info


def peso(Rp, rho_p, g):
    """m g, com m = rho_p 4/3 pi Rp^3."""
    return rho_p * 4.0 / 3.0 * math.pi * Rp ** 3 * g
