# -*- coding: utf-8 -*-
"""
Onda plana x-polarizada, propagando em +z -- BSCs exatos da GLMT.

Ambrosio, Wang & Gouesbet, JQSRT 288, 108266 (2022), Eq. (40), convencao
exp(+i w t):
    g^{+1}_{n,TM} = g^{-1}_{n,TM} = i g^{+1}_{n,TE} = -i g^{-1}_{n,TE} = g_n/2
    g^m_{n,TM} = g^m_{n,TE} = 0  para |m| != 1,      g_n = exp(i k z0).
(Valores identicos ao limite on-axis de besselBeam.BesselBeam.)

Mesma interface legada g_mn(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
que os demais feixes; serve como beam= em asymetryVector.I_x/I_y/I_z.
"""
import numpy as np


class PlaneWave:
    nome = 'onda plana x-polarizada'

    def __init__(self, lamb=1064e-9):
        self.lamb = float(lamb)
        self.k = 2.0 * np.pi / self.lamb
        self.v = 0

    def __repr__(self):
        return f'PlaneWave(lamb={self.lamb*1e9:.1f} nm)'

    def _gn(self, z_0):
        return np.exp(1j * self.k * float(z_0))

    def g_mnTM(self, m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
        if int(v) != 0:
            raise ValueError('onda plana: v deve ser 0')
        return 0.5 * self._gn(z_0) if abs(int(m)) == 1 and abs(int(m)) <= int(n) else 0.0 + 0.0j

    def g_mnTE(self, m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
        if int(v) != 0:
            raise ValueError('onda plana: v deve ser 0')
        m = int(m)
        if m == 1:
            return -0.5j * self._gn(z_0)
        if m == -1:
            return 0.5j * self._gn(z_0)
        return 0.0 + 0.0j
