# -*- coding: utf-8 -*-
"""
Feixe de Bessel CIRCULARMENTE SIMETRICO -- Beam Shape Coefficients (GLMT)
=========================================================================
Referencia: TCC de B. H. D. Rocha (EESC/USP, 2023), Eqs. (3.34) e (3.35),
que seguem Wang et al., Opt. Commun. 387 (2017), Eqs. 7-9:

  g^m_{n,TM} = (1/2) i^{m+1} (-1)^{(m-|m|)/2} (n-m)!/(n+|m|)! e^{i k_z z0}
               x { J_{m-v-1}(Z0) e^{-i(m-v-1)phi0} [ m pi_n^m(cos a) + tau_n^m(cos a)/cos a ]
                 + J_{m-v+1}(Z0) e^{-i(m-v+1)phi0} [ m pi_n^m(cos a) - tau_n^m(cos a)/cos a ] }

  g^m_{n,TE} = -(i/2) i^{m+1} (-1)^{(m-|m|)/2} (n-m)!/(n+|m|)! e^{i k_z z0}
               x { J_{m-v-1}(Z0) e^{-i(m-v-1)phi0} [ m pi_n^m(cos a) + tau_n^m(cos a)/cos a ]
                 - J_{m-v+1}(Z0) e^{-i(m-v+1)phi0} [ m pi_n^m(cos a) - tau_n^m(cos a)/cos a ] }

  com  k_z = k cos a,  k_rho = k sin a  (EXATO),  Z0 = k_rho * rho0,
       pi_n^m(cos a) = P_n^m(cos a)/sin a,  tau_n^m = dP_n^m(cos t)/dt |_{t=a}.

NAO confundir com a Frozen Wave de beamShape.py: aquela e' polarizada
LINEARMENTE EM X (outro modelo), com o k/beta no termo m*pi do colchete TE e
k_rho paraxial. Com um termo so', a FW NAO reduz a este feixe -- e' esperado.

Uso
---
    from besselBeam import BesselBeam
    bb = BesselBeam(alpha=np.radians(5), v=0, lamb=1064e-9)
    I_z(M, x, eps, 0, 0.0, rho0, z0, beam=bb)

Os argumentos `alpha1` e `lamb` da assinatura legada g_mn(...) sao ignorados:
angulo de axicon e comprimento de onda vem do objeto. O argumento `v` e'
conferido contra bb.v (erro se divergir).

Precisao: float64. Para alpha < ALPHA_MIN = 1e-4 rad, pi_n^m = P/sin(alpha)
perde precisao (erro relativo medido: 1e-4 rad -> 3e-9; 1e-6 rad -> 9e-5), por
isso o construtor recusa angulos abaixo disso. Nao ha limite de onda plana.
"""
import numpy as np
import scipy.special as sp

from bsc_common import (legendre_table, P_nm_from_table, tau_from_table,
                        prefator_TM, prefator_TE, fator_nm)

ALPHA_MIN = 1e-4          # rad; abaixo disso float64 perde precisao em pi_n^m
SPOT_FATOR = 4.05         # spot = 4,05 / (k sin alpha)   (TCC, Tabela 2)


class BesselBeam:
    """Feixe de Bessel circularmente simetrico de ordem v (TCC, Eqs. 3.34/3.35)."""

    nome = 'Bessel circularmente simetrico'

    def __init__(self, alpha, v=0, lamb=1064e-9):
        alpha = float(alpha)
        if not (ALPHA_MIN <= alpha < np.pi / 2):
            raise ValueError(
                f'angulo de axicon alpha = {alpha:.3e} rad fora do intervalo suportado '
                f'[{ALPHA_MIN:g}, pi/2). Abaixo de {ALPHA_MIN:g} rad o calculo em float64 '
                f'de pi_n^m = P_n^m(cos a)/sin a perde precisao (limite de onda plana nao '
                f'implementado).')
        self.alpha = alpha
        self.v = int(v)
        self.lamb = float(lamb)
        self.k = 2.0 * np.pi / self.lamb
        self.cos_a = float(np.cos(alpha))
        self.sin_a = float(np.sin(alpha))
        self.k_z = self.k * self.cos_a
        self.k_rho = self.k * self.sin_a

        # tabela de Legendre em cos(alpha): depende so' do feixe, nao do ponto
        self._n_tab = -1
        self._P = None
        self._x = np.array([self.cos_a])
        # contexto do ponto corrente
        self._ponto = None
        self._J = None
        self._ord0 = 0
        self._n_ponto = -1
        self._memo = {}

    # ------------------------------------------------------------ grandezas
    @property
    def spot(self):
        """Delta rho = 4,05 / (k sin alpha)  (TCC, Tabela 2)."""
        return SPOT_FATOR / self.k_rho

    def perfil_escalar(self, rho):
        """|J_v(k sin(alpha) rho)|^2 -- perfil transversal escalar, p/ figuras."""
        return np.abs(sp.jv(self.v, self.k_rho * np.asarray(rho, dtype=float))) ** 2

    def __repr__(self):
        return (f'BesselBeam(alpha={np.degrees(self.alpha):.4g} deg, v={self.v}, '
                f'lamb={self.lamb*1e9:.1f} nm)')

    # -------------------------------------------------------------- tabelas
    def _garantir_legendre(self, n_max):
        if n_max > self._n_tab:
            n_novo = max(n_max, 2 * self._n_tab if self._n_tab > 0 else n_max)
            self._P = legendre_table(self._x, n_novo, n_novo)
            self._n_tab = n_novo

    def _garantir_ponto(self, n_max, rho_0, phi_0, z_0):
        chave = (float(rho_0), float(phi_0), float(z_0))
        if chave != self._ponto or n_max > self._n_ponto:
            n_novo = max(n_max, self._n_ponto) if chave == self._ponto else n_max
            ords = np.arange(-n_novo - abs(self.v) - 2, n_novo + abs(self.v) + 3)
            self._J = sp.jv(ords, self.k_rho * chave[0]).astype(np.float64)
            self._ord0 = int(-ords[0])
            self._fase_z = np.exp(1j * self.k_z * chave[2])
            self._phi0 = chave[1]
            self._ponto = chave
            self._n_ponto = n_novo
            self._memo = {}

    def _conferir_v(self, v):
        if int(v) != self.v:
            raise ValueError(f'ordem v = {v} passada ao calculo difere da ordem do feixe '
                             f'({self.v}). Construa BesselBeam(..., v={v}).')

    # ------------------------------------------------------------- nucleo
    def _g(self, modo, m, n, rho_0, phi_0, z_0):
        m = int(m)
        n = int(n)
        if abs(m) > n:
            return 0.0 + 0.0j              # modo inexistente
        self._garantir_legendre(n + 1)
        self._garantir_ponto(n + 1, rho_0, phi_0, z_0)
        chave = (modo, m, n)
        if chave in self._memo:
            return self._memo[chave]

        # pi e tau exatamente como no TCC (Listings 4.1 e 4.2):
        #   pi  = P_n^m(cos a) / sin(a)
        #   tau = [-(n+1) x P_n^m + (n-m+1) P_{n+1}^m] / sqrt(1-x^2),  x = cos a
        pi = P_nm_from_table(self._P, m, n)[0] / self.sin_a
        tau = tau_from_table(self._P, self._x, m, n)[0]

        v = self.v
        d = self._J[m - v - 1 + self._ord0] * np.exp(-1j * (m - v - 1) * self._phi0)
        f = self._J[m - v + 1 + self._ord0] * np.exp(-1j * (m - v + 1) * self._phi0)
        e = (tau / self.cos_a) + (m * pi)
        g = (-tau / self.cos_a) + (m * pi)

        if modo == 'TM':
            res = prefator_TM(m) * fator_nm(m, n) * self._fase_z * ((d * e) + (f * g))
        else:
            res = prefator_TE(m) * fator_nm(m, n) * self._fase_z * ((d * e) - (f * g))
        self._memo[chave] = res
        return res

    # -------------------------------------------- assinatura legada g_mn(...)
    def g_mnTM(self, m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
        """g^m_{n,TM} -- TCC Eq. (3.34). alpha1/lamb ignorados (vem do objeto)."""
        self._conferir_v(v)
        return self._g('TM', m, n, rho_0, phi_0, z_0)

    def g_mnTE(self, m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
        """g^m_{n,TE} -- TCC Eq. (3.35). alpha1/lamb ignorados (vem do objeto)."""
        self._conferir_v(v)
        return self._g('TE', m, n, rho_0, phi_0, z_0)
