# -*- coding: utf-8 -*-
"""
Campo eletrico INTERNO de uma esfera homogenea (GLMT, exp(+i w t)), com
permeabilidade relativa mu_r arbitraria.

Generaliza tests/verificacao_ras.campo_interno (validado com mu_r = 1:
integral de definicao de J e C_abs) seguindo Ambrosio, JQSRT 255, 107276
(2020), Eqs. (6a)-(6c) e (10a)-(10b):
  - c_n, d_n = coefficients.cs_n / ds_n com mi = mu_r  (Eqs. 10a/10b)
  - termo TE do campo multiplicado por mu_r/M  (Eqs. 6b/6c: mu_sp k/(mu k_sp))
Para mu_r = 1 reduz-se exatamente a rotina validada.

  E_r, E_theta, E_phi (/E0) numa grade produto r x mu(=cos theta) x phi.
"""
import numpy as np
import scipy.special as sp

from bsc_common import legendre_table
from coefficients import cs_n, ds_n
from beams import resolver_bsc


def nmax_wiscombe(x):
    return int(np.ceil(x + 4.05 * x ** (1 / 3) + 2))


def leg(mu, nmax):
    """P[i,m,n], tau[i,m,n], pi[i,m,n] (Condon-Shortley), mu = cos(theta) != +-1."""
    mu = np.atleast_1d(mu).astype(float)
    P = legendre_table(mu, nmax, nmax)
    s = np.sqrt(np.maximum(1 - mu ** 2, 0.0))
    tau = np.zeros((mu.size, nmax + 1, nmax + 1)); pi = np.zeros_like(tau)
    for m in range(nmax + 1):
        for n in range(max(m, 1), nmax + 1):
            tau[:, m, n] = (-(n + 1) * mu * P[:, m, n] + (n - m + 1) * P[:, m, n + 1]) / s
            pi[:, m, n] = P[:, m, n] / s
    return P[:, :, :nmax + 1].copy(), tau, pi


def bscs(beam, nmax, v, rho0, phi0, z0=0.0, lamb=None):
    """Tabelas g[n, m+nmax]. beam=None -> Frozen Wave (beamShape)."""
    gTM, gTE = resolver_bsc(beam)
    if beam is None:
        import beamShape
        beamShape._ACTIVE_CTX = None
    lam = lamb if lamb is not None else getattr(beam, 'lamb', None)
    gtm = np.zeros((nmax + 1, 2 * nmax + 1), complex); gte = np.zeros_like(gtm)
    for n in range(1, nmax + 1):
        for p in range(-n, n + 1):
            gtm[n, p + nmax] = gTM(p, n, v, 0.0, lam, rho0, phi0, 0.0, 0.0, z0)
            gte[n, p + nmax] = gTE(p, n, v, 0.0, lam, rho0, phi0, 0.0, 0.0, z0)
    return gtm, gte


def campo(gtm, gte, M, mu_r, k, nmax, x, r, mu, ph):
    """E_r, E_theta, E_phi / E0 em r (m) x mu x ph."""
    P, tau, pi = leg(mu, nmax)
    z = M * k * np.asarray(r, float)
    Er = np.zeros((z.size, np.size(mu), np.size(ph)), complex)
    Et = np.zeros_like(Er); Ep = np.zeros_like(Er)
    fte = mu_r / M
    for n in range(1, nmax + 1):
        cn = cs_n(M, mu_r, n, x); dn = ds_n(M, mu_r, n, x)
        En = (2 * n + 1) / (n * (n + 1))
        jn = sp.spherical_jn(n, z); jn1 = sp.spherical_jn(n + 1, z)
        psi = z * jn; dpsi = (n + 1) * jn - z * jn1
        f = (-1j) ** (n + 1)
        for p in range(-n, n + 1):
            a = abs(p); gt, ge = gtm[n, p + nmax], gte[n, p + nmax]
            if gt == 0 and ge == 0:
                continue
            e = np.exp(1j * p * np.asarray(ph))[None, None, :]
            Er += (f * (2 * n + 1) * cn * gt * psi / z ** 2)[:, None, None] * P[:, a, n][None, :, None] * e
            Et += (f * En / z)[:, None, None] * (
                cn * gt * dpsi[:, None, None] * tau[:, a, n][None, :, None]
                + p * fte * dn * ge * psi[:, None, None] * pi[:, a, n][None, :, None]) * e
            Ep += (1j * f * En / z)[:, None, None] * (
                p * cn * gt * dpsi[:, None, None] * pi[:, a, n][None, :, None]
                + fte * dn * ge * psi[:, None, None] * tau[:, a, n][None, :, None]) * e
    return Er, Et, Ep


def intensidade(*args, **kw):
    Er, Et, Ep = campo(*args, **kw)
    return np.abs(Er) ** 2 + np.abs(Et) ** 2 + np.abs(Ep) ** 2
