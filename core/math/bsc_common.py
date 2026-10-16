# -*- coding: utf-8 -*-
"""
Utilitarios COMUNS aos Beam Shape Coefficients (BSCs) de todos os feixes.

Extraidos de beamShape.py (2026-09) sem alterar nenhuma operacao: as funcoes
abaixo reproduzem, na mesma ordem de operacoes em ponto flutuante, o que
beamShape fazia internamente. A regressao tests/regressao_fw.py verifica
igualdade BIT A BIT da Frozen Wave antes/depois da extracao.

Usado por:
  beamShape.py   Frozen Wave, polarizacao linear em x   (Eqs. 18/19 do artigo)
  besselBeam.py  feixe de Bessel circularmente simetrico (TCC, Eqs. 3.34/3.35)
"""
import numpy as np
import scipy.special as sp


def _fact_ratio(a: int, b: int) -> float:
    """a!/b! as float (fast, stable)."""
    return float(np.exp(sp.gammaln(a + 1.0) - sp.gammaln(b + 1.0)))


def _P_neg_m_from_pos(P_pos, n: int, m_abs: int):
    """P_n^{-m_abs} from P_n^{m_abs} using the standard relation."""
    return ((-1) ** m_abs) * _fact_ratio(n - m_abs, n + m_abs) * P_pos


def _lpmn(m_max: int, n_max: int, x: float):
    """
    Drop-in replacement for the removed scipy.special.lpmn(m, n, x).
    Returns an array P of shape (m_max+1, n_max+1) with P[m, n] = P_n^m(x)
    for m = 0..m_max, n = 0..n_max (non-negative orders only).
    Inclui a fase de Condon-Shortley (mesma convencao de lpmn/mpmath.legenp).
    """
    table = sp.assoc_legendre_p_all(n_max, m_max, x)[0]
    return table[:, :m_max + 1].T


def legendre_table(x_arr, m_max: int, n_max: int):
    """
    Tabela P[i, m, n] = P_n^m(x_arr[i]), m = 0..m_max, n = 0..n_max+1.
    (n_max+1 porque tau_n^m usa P_{n+1}^m.)
    """
    x_arr = np.atleast_1d(np.asarray(x_arr, dtype=np.float64))
    P = np.empty((x_arr.size, m_max + 1, n_max + 2), dtype=np.float64)
    for i in range(x_arr.size):
        P[i] = _lpmn(m_max, n_max + 1, float(x_arr[i]))
    return P


def P_nm_from_table(P, m: int, n: int):
    """P_n^m para todos os angulos da tabela, aceitando m negativo."""
    m_int = int(m)
    n_int = int(n)
    m_abs = abs(m_int)
    P_pos = P[:, m_abs, n_int]
    if m_int >= 0:
        return P_pos
    return _P_neg_m_from_pos(P_pos, n_int, m_abs)


def tau_from_table(P, x, m: int, n: int):
    """
    tau_n^m(x) = dP_n^m(cos t)/dt  para todos os angulos, x = cos(t):

        tau = [ -(n+1) x P_n^m + (n-m+1) P_{n+1}^m ] / sqrt(1 - x^2)

    (identico a beamShape._FWContext.tau_array e ao tau_mn do TCC, Listing 4.2)
    """
    m_int = int(m)
    n_int = int(n)
    m_abs = abs(m_int)

    Pn_pos = P[:, m_abs, n_int]
    Pn1_pos = P[:, m_abs, n_int + 1]
    if m_int < 0:
        Pn = _P_neg_m_from_pos(Pn_pos, n_int, m_abs)
        Pn1 = _P_neg_m_from_pos(Pn1_pos, n_int + 1, m_abs)
    else:
        Pn = Pn_pos
        Pn1 = Pn1_pos

    num = (-(n_int + 1) * x * Pn) + ((n_int - m_int + 1) * Pn1)

    den = np.sqrt(np.maximum(1.0 - x * x, 0.0))
    den = np.where(den == 0.0, np.sqrt(1.0 - 0.999999999), den)
    return num / den


def pi_from_table(P, x, m: int, n: int):
    """
    pi_n^m(x) = P_n^m(x) / sqrt(1 - x^2)  para todos os angulos, x = cos(t).
    (identico a beamShape._FWContext.pi_array, apos a correcao B1)
    """
    Pn = P_nm_from_table(P, m, n)
    den = np.sqrt(np.maximum(1.0 - x * x, 0.0))
    den = np.where(den == 0.0, np.sqrt(1.0 - 0.999999999), den)
    return Pn / den


def prefator_TM(m: int) -> complex:
    """a1 do modo TM: (1/2) i^(m+1) (-1)^((m-|m|)/2)."""
    m_int = int(m)
    return 0.5 * (1j ** (m_int + 1)) * ((-1) ** ((m_int - abs(m_int)) // 2))


def prefator_TE(m: int) -> complex:
    """a1 do modo TE: -(i/2) i^(m+1) (-1)^((m-|m|)/2)."""
    m_int = int(m)
    return (1j * (-0.5)) * (1j ** (m_int + 1)) * ((-1) ** ((m_int - abs(m_int)) // 2))


def fator_nm(m: int, n: int) -> float:
    """(n-m)! / (n+|m|)!  -- fator comum as Eqs. 3.34/3.35 do TCC e 18/19 da FW."""
    m_int = int(m)
    n_int = int(n)
    return _fact_ratio(n_int - m_int, n_int + abs(m_int))
