# -*- coding: utf-8 -*-
"""
Feixe de Bessel circularmente simetrico POLARIZADO EM x -- BSCs (GLMT)
======================================================================
Referencia: J. J. Wang, T. Wriedt, L. Madler, Y. P. Han, P. Hartmann,
Opt. Commun. 387, 102-109 (2017), Eqs. (15) e (17); convencao exp(+i w t)
(a mesma do projeto). Sun et al., Appl. Opt. 61, 7632 (2022), Eq. (24) traz a
mesma expressao (com e^{-i k_z z0}).

  g^m_{n,TM} = -g(a0) (-1)^{(m-|m|)/2} (n-m)!/(n+|m|)! e^{i k_z z0}
               { i^{l-m+1} e^{i(l-m+1)phi0} J_{l-m+1}(s0) [tau_n^m + m pi_n^m]
               + i^{l-m-1} e^{i(l-m-1)phi0} J_{l-m-1}(s0) [tau_n^m - m pi_n^m] }
  g^m_{n,TE} =  i g(a0) (-1)^{(m-|m|)/2} (n-m)!/(n+|m|)! e^{i k_z z0}
               { i^{l-m+1} e^{i(l-m+1)phi0} J_{l-m+1}(s0) [tau_n^m + m pi_n^m]
               - i^{l-m-1} e^{i(l-m-1)phi0} J_{l-m-1}(s0) [tau_n^m - m pi_n^m] }

  s0 = k sin(a0) rho0, (rho0, phi0) = posicao do centro do feixe no referencial
  da particula; pi_n^m = P_n^m(cos a0)/sin a0, tau_n^m = dP_n^m/d a0; P_n^m com
  fase de Condon-Shortley (Eq. 3 de Wang 2017).
  g(a0) = 1/2: feixe aplanatico (espectro angular); (1+cos a0)/4: Davis.

Ordens negativas: `neg` define P_n^m para m<0:
  'rel' -> P_n^{-|m|} = (-1)^{|m|}(n-|m|)!/(n+|m|)! P_n^{|m|}  (padrao do projeto)
  'abs' -> P_n^{|m|}
A escolha foi fixada pela reconstrucao do campo incidente contra a Eq. (1) de
Wang 2017 (tests/validacao_besselX.py).
"""
import numpy as np
import scipy.special as sp

from bsc_common import legendre_table, P_nm_from_table, tau_from_table, fator_nm

ALPHA_MIN = 1e-4


class BesselBeamX:
    nome = 'Bessel circularmente simetrico, polarizacao x (Wang 2017)'

    def __init__(self, alpha, l=0, lamb=9.58e-6, g=0.5, neg='rel'):
        alpha = float(alpha)
        if not (ALPHA_MIN <= alpha < np.pi / 2):
            raise ValueError('alpha fora de [1e-4, pi/2)')
        self.alpha = alpha; self.l = int(l); self.v = int(l)
        self.lamb = float(lamb); self.k = 2 * np.pi / self.lamb
        self.g = float(g); self.neg = neg
        self.cos_a = float(np.cos(alpha)); self.sin_a = float(np.sin(alpha))
        self.k_z = self.k * self.cos_a; self.k_rho = self.k * self.sin_a
        self._x = np.array([self.cos_a]); self._n_tab = -1; self._P = None
        self._ponto = None; self._memo = {}

    @property
    def spot(self):
        return 4.05 / self.k_rho

    def _leg(self, m, n):
        if n > self._n_tab:
            nn = max(n, 2 * max(self._n_tab, 1))
            self._P = legendre_table(self._x, nn, nn); self._n_tab = nn
        mm = m if (m >= 0 or self.neg == 'rel') else -m
        pi = P_nm_from_table(self._P, mm, n)[0] / self.sin_a
        tau = tau_from_table(self._P, self._x, mm, n)[0]
        return pi, tau

    def _ctx(self, rho_0, phi_0, z_0, n):
        chave = (float(rho_0), float(phi_0), float(z_0))
        if chave != self._ponto or n > self._nJ:
            nn = max(n, 2 * getattr(self, '_nJ', 1)) if chave == self._ponto else n
            ords = np.arange(-nn - abs(self.l) - 3, nn + abs(self.l) + 4)
            self._J = sp.jv(ords, self.k_rho * chave[0]); self._o0 = int(-ords[0])
            self._fz = np.exp(1j * self.k_z * chave[2]); self._phi0 = chave[1]
            self._ponto = chave; self._nJ = nn; self._memo = {}

    def _g(self, modo, m, n, rho_0, phi_0, z_0):
        m = int(m); n = int(n)
        if abs(m) > n:
            return 0j
        self._ctx(rho_0, phi_0, z_0, n + 1)
        key = (modo, m, n)
        if key in self._memo:
            return self._memo[key]
        pi, tau = self._leg(m, n)
        l = self.l
        u = (1j ** (l - m + 1)) * np.exp(1j * (l - m + 1) * self._phi0) * self._J[l - m + 1 + self._o0] * (tau + m * pi)
        w = (1j ** (l - m - 1)) * np.exp(1j * (l - m - 1) * self._phi0) * self._J[l - m - 1 + self._o0] * (tau - m * pi)
        pref = ((-1) ** ((m - abs(m)) // 2)) * fator_nm(m, n) * self._fz * self.g
        res = -pref * (u + w) if modo == 'TM' else 1j * pref * (u - w)
        self._memo[key] = res
        return res

    def _conferir(self, v):
        if int(v) != self.l:
            raise ValueError(f'ordem v={v} difere da ordem do feixe l={self.l}')

    def g_mnTM(self, m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
        self._conferir(v); return self._g('TM', m, n, rho_0, phi_0, z_0)

    def g_mnTE(self, m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
        self._conferir(v); return self._g('TE', m, n, rho_0, phi_0, z_0)

    # campo analitico (Wang 2017, Eq. 1) para validacao
    def campo(self, x, y, z, x0, y0, z0=0.0):
        rhoG = np.hypot(x - x0, y - y0); phG = np.arctan2(y - y0, x - x0)
        s = self.k_rho * rhoG; l = self.l; ca = self.cos_a
        J = lambda o: sp.jv(o, s)
        e = lambda o: (-1j) ** o * np.exp(1j * o * phG) * J(o)
        ph = self.g * np.exp(-1j * self.k_z * (z - z0))
        Ex = ph * ((1 + ca) * e(l) - 0.5 * (1 - ca) * (e(l - 2) + e(l + 2)))
        Ey = ph * (1 - ca) / (2j) * (e(l - 2) - e(l + 2))
        Ez = -ph * self.sin_a * (e(l - 1) + e(l + 1))
        return Ex, Ey, Ez
