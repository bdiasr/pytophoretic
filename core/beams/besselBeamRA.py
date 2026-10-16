# -*- coding: utf-8 -*-
"""
Feixe de Bessel circularmente simetrico nos modos TM (polarizacao radial, l=0)
e TE (polarizacao azimutal, l=0) -- BSCs (GLMT)
============================================================================
Referencia: Sun et al., Appl. Opt. 61, 7632 (2022), Eqs. (29), (30) e (32),
obtidas pela abordagem de Wang et al. 2017. O fator de fase em z0 foi escrito
como e^{+i k_z z0}, como na Eq. (15) de Wang 2017 (convencao exp(+i w t));
nas figuras usadas z0 = 0.

  modo TM (forma VALIDADA, convencao exp(+i w t)):
            g^m_{n,TM} = -g (-1)^{(m-|m|)/2} (n-m)!/(n+|m|)! e^{i k_z z0}
                          J_{m-l}(k_rho rho0) e^{i(l-m)phi0} i^{m-l}   tau_n^m
            g^m_{n,TE} =  g (-1)^{(m-|m|)/2} (n-m)!/(n+|m|)! e^{i k_z z0}
                          i^{m-l+1} J_{m-l}(k_rho rho0) e^{i(l-m)phi0} m pi_n^m
  Diferenca em relacao a Sun Eqs. (29)-(30), escritas em exp(-i w t): as potencias
  i^{l-m} e i^{l-m-1} viram i^{m-l} e i^{m-l+1} (i -> -i). Com a forma impressa,
  o campo reconstruido falha na translacao off-axis (erro ~0,6-0,9) e, para l != 0,
  E_z nao e' proporcional a J_l; com a forma acima ambos os testes dao ~1e-14.
  modo TE (Eq. 32):  g_TM^{(TE)} = g_TE^{(TM)},  g_TE^{(TE)} = -g_TM^{(TM)}.
Validacao: tests/validacao_besselRA.py.
"""
import numpy as np
import scipy.special as sp
from bsc_common import legendre_table, P_nm_from_table, tau_from_table, fator_nm


class BesselBeamRA:
    def __init__(self, alpha, l=0, lamb=9.58e-6, g=0.5, modo='TM'):
        self.alpha = float(alpha); self.l = int(l); self.v = int(l)
        self.lamb = float(lamb); self.k = 2 * np.pi / self.lamb; self.g = float(g)
        self.modo = modo.upper()
        self.nome = 'Bessel %s (%s)' % ('radial' if self.modo == 'TM' else 'azimutal', self.modo)
        self.cos_a = float(np.cos(alpha)); self.sin_a = float(np.sin(alpha))
        self.k_z = self.k * self.cos_a; self.k_rho = self.k * self.sin_a
        self._x = np.array([self.cos_a]); self._n_tab = -1; self._ponto = None; self._memo = {}; self._nJ = -1

    @property
    def spot(self):
        return 4.05 / self.k_rho

    def _leg(self, m, n):
        if n > self._n_tab:
            nn = max(n, 2 * max(self._n_tab, 1)); self._P = legendre_table(self._x, nn, nn); self._n_tab = nn
        return P_nm_from_table(self._P, m, n)[0] / self.sin_a, tau_from_table(self._P, self._x, m, n)[0]

    def _ctx(self, rho_0, phi_0, z_0, n):
        chave = (float(rho_0), float(phi_0), float(z_0))
        if chave != self._ponto or n > self._nJ:
            nn = max(n, 2 * self._nJ) if chave == self._ponto else n
            ords = np.arange(-nn - abs(self.l) - 3, nn + abs(self.l) + 4)
            self._J = sp.jv(ords, self.k_rho * chave[0]); self._o0 = int(-ords[0])
            self._fz = np.exp(1j * self.k_z * chave[2]); self._phi0 = chave[1]
            self._ponto = chave; self._nJ = nn; self._memo = {}

    def _gTMmode(self, qual, m, n, rho_0, phi_0, z_0):
        m = int(m); n = int(n)
        if abs(m) > n:
            return 0j
        self._ctx(rho_0, phi_0, z_0, n + 1)
        key = (qual, m, n)
        if key in self._memo:
            return self._memo[key]
        pi, tau = self._leg(m, n); l = self.l
        base = ((-1) ** ((m - abs(m)) // 2)) * fator_nm(m, n) * self._fz * self.g \
            * self._J[m - l + self._o0] * np.exp(1j * (l - m) * self._phi0)
        # Potencias de i na convencao exp(+i w t): i^{m-l} e i^{m-l+1}. As Eqs. (29)-(30)
        # de Sun (exp(-i w t)) trazem i^{l-m} e i^{l-m-1}; a forma impressa NAO passa na
        # validacao (feixe on-axis de ordem l != 0 e translacao off-axis). 
        if qual == 'TM':
            res = -base * (1j ** (m - l)) * tau
        else:
            res = base * (1j ** (m - l + 1)) * m * pi
        self._memo[key] = res
        return res

    def g_mnTM(self, m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
        if self.modo == 'TM':
            return self._gTMmode('TM', m, n, rho_0, phi_0, z_0)
        return self._gTMmode('TE', m, n, rho_0, phi_0, z_0)

    def g_mnTE(self, m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
        if self.modo == 'TM':
            return self._gTMmode('TE', m, n, rho_0, phi_0, z_0)
        return -self._gTMmode('TM', m, n, rho_0, phi_0, z_0)
