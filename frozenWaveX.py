# -*- coding: utf-8 -*-
"""
Frozen Wave POLARIZADA EM x construida a partir do feixe de Bessel de Wang 2017
===============================================================================
Referencia da deducao: relatorios/deducao_bscs_coeficientes.tex, Secao 3.6,
Eqs. (24)-(25). A FW e' a superposicao de 2N+1 feixes de Bessel polarizados em x
(Wang et al., Opt. Commun. 387, 102 (2017), Eqs. 15/17) com

  beta_q = Q + 2 pi q / L,   cos(alpha_q) = beta_q / k,   k_rho_q = sqrt(k^2 - beta_q^2)

e pesos A_q (Aq.py, perfil ativo). Por linearidade da projecao (Eq. 5 da deducao):

  g^m_{n,X}(FW) = sum_q A_q g^m_{n,X}(Bessel x; alpha_q, ordem v) com e^{i k_z z0} -> e^{i beta_q z0}.

Diferencas em relacao a FW de beamShape.py (NAO alterada): modelo vetorial de
Bessel de Wang (colchetes [tau +/- m pi] no TM e no TE), k_rho_q EXATO (nao o
paraxial h_q), amplitude g(alpha_q) (1/2 = aplanatico).

Interface identica aos demais feixes: g_mnTM / g_mnTE(m, n, v, alpha1, lamb,
rho_0, phi_0, x_0, y_0, z_0); alpha1 e lamb sao ignorados.
Validacao: testes/validacao_frozenWaveX.py.
"""
import numpy as np
import scipy.special as sp

import cte
from Aq import Aq_2_vec
from bsc_common import legendre_table, P_nm_from_table, tau_from_table, fator_nm


class FrozenWaveX:
    nome = 'Frozen Wave pol. x (Bessel de Wang 2017)'

    def __init__(self, v=0, lamb=None, Q=None, L=None, N=None, g=0.5, Aq=None, qs=None):
        self.lamb = float(cte.lamb if lamb is None else lamb)
        self.k = 2 * np.pi / self.lamb
        self.Q = float(cte.Q if Q is None else Q)
        self.L = float(cte.L if L is None else L)
        N = int(cte.N_FW if N is None else N)
        self.qs = np.arange(-N, N + 1) if qs is None else np.asarray(qs, dtype=int)
        self.v = int(v); self.l = self.v; self.g = float(g)
        self.beta = self.Q + 2 * np.pi * self.qs / self.L
        if np.any(self.beta <= 0) or np.any(self.beta >= self.k):
            raise ValueError('beta_q fora de (0, k): componentes evanescentes ou contrapropagantes')
        self.A = (Aq_2_vec(self.qs) if Aq is None else np.asarray(Aq)).astype(complex)
        self.cos_a = self.beta / self.k
        self.k_rho = np.sqrt(self.k ** 2 - self.beta ** 2)
        self.sin_a = self.k_rho / self.k
        self._n_tab = -1; self._ponto = None; self._memo = {}; self._nJ = -1

    def _tab(self, n):
        if n > self._n_tab:
            nn = max(n, 2 * max(self._n_tab, 1))
            self._P = legendre_table(self.cos_a, nn, nn); self._n_tab = nn

    def _ctx(self, rho_0, phi_0, z_0, n):
        chave = (float(rho_0), float(phi_0), float(z_0))
        if chave != self._ponto or n > self._nJ:
            nn = max(n, 2 * self._nJ) if chave == self._ponto else n
            ords = np.arange(-nn - abs(self.v) - 3, nn + abs(self.v) + 4)
            self._J = sp.jv(ords[:, None], (self.k_rho * chave[0])[None, :])
            self._o0 = int(-ords[0])
            self._ker = self.A * self.g * np.exp(1j * self.beta * chave[2])
            self._phi0 = chave[1]; self._ponto = chave; self._nJ = nn; self._memo = {}

    def _g(self, modo, m, n, rho_0, phi_0, z_0):
        m = int(m); n = int(n)
        if abs(m) > n:
            return 0j
        self._tab(n + 1); self._ctx(rho_0, phi_0, z_0, n + 1)
        key = (modo, m, n)
        if key in self._memo:
            return self._memo[key]
        pi = P_nm_from_table(self._P, m, n) / self.sin_a
        tau = tau_from_table(self._P, self.cos_a, m, n)
        l = self.v
        u = (1j ** (l - m + 1)) * np.exp(1j * (l - m + 1) * self._phi0) * self._J[l - m + 1 + self._o0] * (tau + m * pi)
        w = (1j ** (l - m - 1)) * np.exp(1j * (l - m - 1) * self._phi0) * self._J[l - m - 1 + self._o0] * (tau - m * pi)
        pref = ((-1) ** ((m - abs(m)) // 2)) * fator_nm(m, n)
        s = np.sum(self._ker * ((u + w) if modo == 'TM' else (u - w)))
        res = -pref * s if modo == 'TM' else 1j * pref * s
        self._memo[key] = res
        return res

    def _conferir(self, v):
        if int(v) != self.v:
            raise ValueError(f'ordem v={v} difere da ordem do feixe {self.v}')

    def g_mnTM(self, m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
        self._conferir(v); return self._g('TM', m, n, rho_0, phi_0, z_0)

    def g_mnTE(self, m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
        self._conferir(v); return self._g('TE', m, n, rho_0, phi_0, z_0)

    def campo(self, x, y, z, x0, y0, z0=0.0):
        """Campo analitico: soma ponderada das Eqs. (1) de Wang 2017 (validacao)."""
        rhoG = np.hypot(x - x0, y - y0); phG = np.arctan2(y - y0, x - x0); l = self.v
        Ex = Ey = Ez = 0j
        for Aq, b, ca, sa, kr in zip(self.A, self.beta, self.cos_a, self.sin_a, self.k_rho):
            s = kr * rhoG
            e = lambda o: (-1j) ** o * np.exp(1j * o * phG) * sp.jv(o, s)
            ph = Aq * self.g * np.exp(-1j * b * (z - z0))
            Ex = Ex + ph * ((1 + ca) * e(l) - 0.5 * (1 - ca) * (e(l - 2) + e(l + 2)))
            Ey = Ey + ph * (1 - ca) / (2j) * (e(l - 2) - e(l + 2))
            Ez = Ez - ph * sa * (e(l - 1) + e(l + 1))
        return Ex, Ey, Ez
