# -*- coding: utf-8 -*-
"""
fw_field.py -- Campo escalar de Frozen Waves de ordem arbitraria v.

Implementa as duas equacoes do artigo de forma GENERICA no perfil F(z), na
ordem v e no truncamento N -- ao contrario de frozenWave.py/Aq.py, onde cada
perfil exige uma Psi_N/Aq_N dedicada e a ordem esta fixa em J_0.

    Psi(rho, phi, z) = sum_{n=-N}^{N} A_n J_v(k_rho_n rho) e^{i v phi} e^{-i beta_n z}   (1)

    A_n = (1/L) int_0^L F(z) e^{i 2 pi n z / L} dz                                       (2)

com beta_n = Q + 2 pi n / L  e  k_rho_n = sqrt(k^2 - beta_n^2).

Convencoes e escolhas (registradas para rastreabilidade):

* k_rho_n -- o nucleo de producao (cte.h_q) usa a forma PARAXIAL
  h_q = sqrt(2) k sqrt(1 - beta/k), que aproxima sqrt(k^2-beta^2) quando
  beta ~ k. Erro medido no grid deste projeto: 1.6e-06 a 1.3e-05 ao longo dos
  47 modos. O default aqui e' `paraxial=True` para que Psi fique CONSISTENTE
  com os BSCs que ele serve para validar; `paraxial=False` usa a forma exata
  da Eq. (1).

* A quadratura de A_n reproduz a de Aq.py: Simpson composto com 301 amostras
  uniformes em [0, L], mesmo sinal (+i) no expoente.

Este modulo NAO altera frozenWave.py nem Aq.py -- o caminho critico dos BSCs
ja validados permanece intocado.
"""
import numpy as np
import scipy.special as sp

from cte import k as K_WAVE, L as L_DEFAULT, Q as Q_DEFAULT, N_FW as N_DEFAULT


# -----------------------------------------------------------------------------
# Eq. (2): coeficientes A_n
# -----------------------------------------------------------------------------
_SIMPSON_N = 301          # mesmo numero de amostras de Aq.py
_AN_CACHE = {}            # (id_perfil, N, L) -> array de A_n


def _simpson_uniform(y, dz):
    """Simpson composto em grid uniforme com numero IMPAR de amostras."""
    n = y.shape[-1]
    if n < 3 or (n % 2) == 0:
        return np.trapezoid(y, dx=dz, axis=-1) if hasattr(np, 'trapezoid') else np.trapz(y, dx=dz, axis=-1)
    s = y[..., 0] + y[..., -1] + 4.0 * y[..., 1:-1:2].sum(-1) + 2.0 * y[..., 2:-2:2].sum(-1)
    return (dz / 3.0) * s


def An(F, N=N_DEFAULT, L=None, cache_key=None):
    """
    Eq. (2): A_n = (1/L) int_0^L F(z) exp(i 2 pi n z / L) dz, para n = -N..N.

    F         : callable F(z) vetorizado em z (array -> array)
    N         : truncamento; devolve 2N+1 coeficientes
    L         : comprimento; default cte.L
    cache_key : chave opcional para memoizar (ex.: nome do perfil)

    Retorna array complexo de shape (2N+1,), indexado n = -N ... +N.
    """
    Lf = float(L_DEFAULT if L is None else L)
    key = (cache_key if cache_key is not None else id(F), int(N), Lf)
    if key in _AN_CACHE:
        return _AN_CACHE[key]

    z = np.linspace(0.0, Lf, _SIMPSON_N, dtype=np.float64)
    dz = z[1] - z[0]
    Fz_vals = np.asarray(F(z), dtype=np.complex128)

    n_arr = np.arange(-int(N), int(N) + 1, dtype=np.int64)
    phase = np.exp(1j * 2.0 * np.pi * n_arr[:, None] * z[None, :] / Lf)

    out = (1.0 / Lf) * _simpson_uniform(phase * Fz_vals[None, :], dz)
    out = np.asarray(out, dtype=np.complex128)
    _AN_CACHE[key] = out
    return out


# -----------------------------------------------------------------------------
# Numeros de onda
# -----------------------------------------------------------------------------
def beta_n(N=N_DEFAULT, L=None, Q=None):
    """beta_n = Q + 2 pi n / L, para n = -N..N."""
    Lf = float(L_DEFAULT if L is None else L)
    Qf = float(Q_DEFAULT if Q is None else Q)
    n_arr = np.arange(-int(N), int(N) + 1, dtype=np.float64)
    return Qf + 2.0 * np.pi * n_arr / Lf


def k_rho_n(N=N_DEFAULT, L=None, Q=None, k=None, paraxial=True):
    """
    k_rho_n. paraxial=True  -> sqrt(2) k sqrt(1 - beta/k)   (forma de cte.h_q)
             paraxial=False -> sqrt(k^2 - beta^2)           (forma da Eq. 1)
    """
    kk = float(K_WAVE if k is None else k)
    b = beta_n(N, L, Q)
    if paraxial:
        return np.sqrt(2.0) * kk * np.sqrt(1.0 - b / kk + 0j)
    return np.sqrt(kk ** 2 - b ** 2 + 0j)


_warned = {}


def n_max_propagante(L=None, Q=None, k=None):
    """
    Maior N com todos os modos propagantes: beta_n = Q + 2 pi n/L <= k
    =>  n <= (k - Q) L / (2 pi).   (condicao 0 <= beta_n <= w0/c do artigo)
    """
    Lf = float(L_DEFAULT if L is None else L)
    Qf = float(Q_DEFAULT if Q is None else Q)
    kk = float(K_WAVE if k is None else k)
    return int(np.floor((kk - Qf) * Lf / (2.0 * np.pi)))


# -----------------------------------------------------------------------------
# Eq. (1): campo Psi
# -----------------------------------------------------------------------------
def Psi(F, rho, z, v=0, phi=0.0, N=N_DEFAULT, L=None, Q=None, k=None,
        paraxial=True, cache_key=None):
    """
    Eq. (1): Psi(rho,phi,z) = sum_n A_n J_v(k_rho_n rho) e^{i v phi} e^{-i beta_n z}

    F    : callable F(z) do perfil longitudinal desejado
    rho  : escalar ou array (m). Aceita valores negativos (usa |rho|; para v
           impar J_v e' impar, entao o sinal e' propagado).
    z    : escalar ou array (m)
    v    : ordem do feixe (carga topologica)
    phi  : angulo azimutal (rad)

    rho e z sao "broadcast" um contra o outro: passe rho[:,None] e z[None,:]
    para obter a malha 2D de uma vez.

    Retorna array complexo com o shape do broadcast de (rho, z).
    """
    v = int(v)
    A = An(F, N=N, L=L, cache_key=cache_key)             # (2N+1,)
    b = beta_n(N, L, Q)                                   # (2N+1,)
    kr = k_rho_n(N, L, Q, k, paraxial)                    # (2N+1,) complexo

    rho = np.asarray(rho, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)

    # Modos evanescentes: a condicao 0 <= beta_n <= k (artigo, Sec. 2) delimita
    # N. Acima disso beta_n > k, k_rho_n vira imaginario puro e J_v(k_rho rho)
    # deixa de ser uma onda transversal propagante -- vira I_v (Bessel
    # modificada), que CRESCE com rho. Tomar np.real() ali silenciaria o
    # problema tratando o modo como propagante com k_rho = 0.
    #
    # Aqui: detectamos, avisamos uma vez, e avaliamos J_v no argumento
    # complexo (scipy.special.jv aceita), que e' o prolongamento analitico
    # correto -- assim o resultado e' honesto sobre o que esses modos fazem.
    n_evan = int(np.count_nonzero(np.abs(np.imag(kr)) > 0))
    if n_evan and not _warned.get('evan'):
        import warnings as _w
        _w.warn(
            f"N={N} inclui {n_evan} modo(s) evanescente(s): beta_n > k viola "
            f"0 <= beta_n <= k (artigo, Sec. 2). N maximo propagante = "
            f"{n_max_propagante(L, Q, k)}.",
            RuntimeWarning, stacklevel=2)
        _warned['evan'] = True

    # J_v(k_rho |rho|), com o sinal de rho propagado quando v e' impar
    arg = kr[:, None] * np.abs(rho).ravel()[None, :]            # (2N+1, Nrho)
    Jv = sp.jv(v, arg)
    if v % 2 == 1:
        Jv = Jv * np.sign(rho).ravel()[None, :]

    # e^{-i beta_n z}
    expz = np.exp(-1j * b[:, None] * z.ravel()[None, :])        # (2N+1, Nz)

    # soma_n A_n Jv[n, i] expz[n, j]  -> (Nrho, Nz)
    W = A[:, None] * Jv                                          # (2N+1, Nrho)
    out = np.einsum('ni,nj->ij', W, expz)

    out = out * np.exp(1j * v * float(phi))

    # Shape de saida: sempre a malha (Nrho, Nz) do produto cartesiano.
    # Escalares sao reduzidos, para que Psi(F, 0.0, 0.175) devolva um escalar.
    if rho.ndim == 0 and z.ndim == 0:
        return out[0, 0]
    if rho.ndim == 0:
        return out[0, :]
    if z.ndim == 0:
        return out[:, 0]
    return out


def intensity(F, rho, z, **kw):
    """|Psi|^2 para o perfil F."""
    return np.abs(Psi(F, rho, z, **kw)) ** 2
