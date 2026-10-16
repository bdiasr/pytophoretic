import numpy as np
import scipy.special as sp


def _fact_ratio(a: int, b: int) -> float:
    """
    Returns a!/b! as float using gammaln (fast, stable).
    """
    return float(np.exp(sp.gammaln(a + 1.0) - sp.gammaln(b + 1.0)))


def _Pnm(n: int, m: int, x):
    """
    Associated Legendre P_n^m(x) with support for negative m via:
      P_n^{-m}(x) = (-1)^m * (n-m)!/(n+m)! * P_n^m(x), for m>0
    """
    m_int = int(m)
    n_int = int(n)
    x = np.asarray(x)
    if m_int >= 0:
        return sp.lpmv(m_int, n_int, x)
    m_abs = -m_int
    P_pos = sp.lpmv(m_abs, n_int, x)
    # (-1)^m_abs * (n - m_abs)!/(n + m_abs)! * P_n^{m_abs}
    return ((-1) ** m_abs) * _fact_ratio(n_int - m_abs, n_int + m_abs) * P_pos


def pi_n(m, n, x):
    """
    Fast version of legacy pi_n using lpmn output.
    Kept for compatibility (not the main hot path).
    """
    x = float(x)
    if x == 1.0:
        denom = np.sqrt(1.0 - 0.99999999999999)
    else:
        denom = np.sqrt(1.0 - x * x)
    if denom == 0:
        return 0.0
    # This matches the intent of the original: P_n^m(x)/sqrt(1-x^2)
    return _Pnm(int(n), int(m), x) / denom


def tau_n(m, n, x):
    x = float(x)
    denom = np.sqrt(1.0 - x * x)
    if denom == 0:
        return 0.0
    # d/dx P_n^m(x) is not used in the main code, so keep a simple stable definition:
    # tau_n = -sqrt(1-x^2) * dP/dx
    # Use scipy.special.lpmn to get derivative cheaply.
    m_abs = abs(int(m))
    n_int = int(n)
    P, dP = sp.lpmn(m_abs, n_int, x)
    dPmn = dP[m_abs, n_int]
    return -(np.sqrt(1.0 - x * x)) * dPmn


def pi_mn(m, n, x):
    """
    IMPORTANT: This preserves the ORIGINAL semantics in your repo:
      pi_mn(m,n,x) = P_n^m(cos(x)) / sin(x)
    even though for Frozen Waves you sometimes pass x = beta/k (already a cos).
    """
    x = float(x)
    s = np.sin(x)
    if s == 0.0:
        return 0.0
    return _Pnm(int(n), int(m), np.cos(x)) / s


def tau_mn(m, n, x):
    """
    IMPORTANT: This preserves the ORIGINAL semantics in your repo:
      tau_mn uses x directly in P_n^m(x) and P_{n+1}^m(x) with the standard relation.
    """
    x = float(x)
    pmn = _Pnm(int(n), int(m), x)
    pmn1 = _Pnm(int(n + 1), int(m), x)

    first = -(int(n) + 1) * x * pmn
    second = (int(n) - int(m) + 1) * pmn1
    num = first + second

    if x == 1.0:
        den = np.sqrt(1.0 - 0.999999999)
    else:
        den = np.sqrt(1.0 - x * x)

    if den == 0.0:
        return 0.0

    return num / den
