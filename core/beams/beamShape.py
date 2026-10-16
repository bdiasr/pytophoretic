import numpy as np
import scipy.special as sp

from cte import qs, Q, L, k
from Aq import Aq_2_vec

# -----------------------------------------------------------------------------
# Fast Beam Shape Coefficients for Frozen Waves -- POLARIZACAO LINEAR EM X
# -----------------------------------------------------------------------------
# Formulacao validada (Eqs. 18/19 do artigo): NAO ALTERAR. Refatoracoes devem
# passar em tests/regressao_fw.py (igualdade bit a bit).
#
# Difere do feixe de Bessel do TCC (besselBeam.py) por ser OUTRO modelo de
# polarizacao: aqui o colchete TE tem k/beta no termo m*pi; no Bessel
# circularmente simetrico o 1/cos(alpha) fica no tau, como no TM. Tambem usa
# k_rho paraxial h_q = sqrt(2) k sqrt(1 - beta_q/k).
# -----------------------------------------------------------------------------
# Hot path: g_mnTM / g_mnTE.
#
# Main optimisations:
#   1) Aq_2(q) is cached and vectorised (see Aq.py)
#   2) For a fixed point (rho_0, phi_0, z_0, v), we precompute:
#        - beta_q, h_q, kernel Aq_2(q) * exp(i beta_q z0)
#        - Legendre tables needed by pi_mn/tau_mn (preserving original semantics)
#        - Bessel J_v(order, h_q * rho_0) for all orders in a single broadcast call
#   3) g_mn values are memoised per (kind, m, n) for the active point.
#
# This keeps the original formulas but avoids millions of repeated special
# function calls.
# -----------------------------------------------------------------------------


# Utilitarios comuns a todos os BSCs (extraidos para bsc_common.py em 2026-09,
# sem alterar nenhuma operacao). Os nomes continuam importaveis daqui por
# compatibilidade (ex.: testePolarizacao.py importa _P_neg_m_from_pos).
from bsc_common import (_fact_ratio, _P_neg_m_from_pos, _lpmn,  # noqa: F401
                        legendre_table, tau_from_table, pi_from_table,
                        prefator_TM, prefator_TE, fator_nm)


class _FWContext:
    __slots__ = (
        "v", "rho_0", "phi_0", "z_0",
        "n_max", "m_max",
        "q_arr", "beta", "x1", "x2", "sin_x1", "k_over_beta",
        "kernel", "h", "arg",
        "P1",
        "orders", "order0", "J",
        "cache_tm", "cache_te",
    )

    def __init__(self, n_max: int, v, rho_0, phi_0, z_0):
        self.v = int(v)
        self.rho_0 = float(rho_0)
        self.phi_0 = float(phi_0)
        self.z_0 = float(z_0)

        self.n_max = int(n_max)
        self.m_max = int(n_max)  # |m| <= n

        self.q_arr = np.asarray(qs, dtype=np.int64)
        self.beta = (float(Q) + (2.0 * np.pi * self.q_arr / float(L))).astype(np.float64)
        self.x1 = (self.beta / float(k)).astype(np.float64)   # passed to tau_mn (original semantics)
        self.x2 = np.cos(self.x1).astype(np.float64)           # passed to pi_mn via cos(x)
        self.sin_x1 = np.sin(self.x1).astype(np.float64)

        self.k_over_beta = (float(k) / self.beta).astype(np.float64)

        # h_q = sqrt(2)*k*sqrt(1 - beta/k)   (original cte.h_q)
        self.h = (np.sqrt(2.0) * float(k) * np.sqrt(1.0 - self.x1 + 0j)).astype(np.complex128)
        self.arg = (self.h * self.rho_0).astype(np.complex128)

        # kernel(q) = Aq_2(q) * exp(i * beta_q * z0)
        Aq_vals = Aq_2_vec(self.q_arr).astype(np.complex128)
        self.kernel = Aq_vals * np.exp(1j * self.beta * self.z_0)

        # Tabela de Legendre P_n^m(x1), x1 = beta_q/k = cos(alpha_q), usada por
        # tau_n^m e pi_n^m. (A antiga tabela P2, em cos(x1), so' servia ao pi
        # anterior a correcao B1 e deixou de ser usada.)
        self.P1 = legendre_table(self.x1, self.m_max, self.n_max)

        # Precompute Bessel J for all needed orders at once:
        # Orders needed: m - v ± 1 for m in [-n_max, n_max]
        min_order = (-self.n_max - abs(self.v) - 2)
        max_order = (self.n_max + abs(self.v) + 2)
        self.orders = np.arange(min_order, max_order + 1, dtype=np.int64)
        self.order0 = int(-min_order)  # index shift
        # Broadcast: (Norders, 1) vs (1, Nq) => (Norders, Nq)
        self.J = sp.jv(self.orders[:, None], self.arg[None, :]).astype(np.complex128)

        self.cache_tm = {}
        self.cache_te = {}

    def tau_array(self, m: int, n: int):
        """tau_n^m(x1) para todo q (ver bsc_common.tau_from_table)."""
        return tau_from_table(self.P1, self.x1, m, n)

    def pi_array(self, m: int, n: int):
        """
        pi_n^m(x1) = P_n^m(x1) / sqrt(1 - x1^2) para todo q.

        CORRECAO B1 (2026-09): a versao anterior avaliava P_n^m(cos(x1))/sin(x1),
        tratando x1 = beta_q/k como se fosse o ANGULO (cos/sin aplicados duas
        vezes; erro de ~539x a partir de n = 2). Ver bsc_common.pi_from_table.
        """
        return pi_from_table(self.P1, self.x1, m, n)


# Global active context (kept simple and deterministic)
_ACTIVE_CTX = None


def _get_ctx(n: int, v, rho_0, phi_0, z_0) -> _FWContext:
    """
    Get a context for the current evaluation point; rebuild only when needed.
    """
    global _ACTIVE_CTX
    n_needed = int(n) + 1  # some formulas require n+1
    if _ACTIVE_CTX is None:
        _ACTIVE_CTX = _FWContext(n_needed, v, rho_0, phi_0, z_0)
        return _ACTIVE_CTX

    same_point = (
        _ACTIVE_CTX.v == int(v) and
        _ACTIVE_CTX.rho_0 == float(rho_0) and
        _ACTIVE_CTX.phi_0 == float(phi_0) and
        _ACTIVE_CTX.z_0 == float(z_0)
    )
    if (not same_point) or (_ACTIVE_CTX.n_max < n_needed):
        _ACTIVE_CTX = _FWContext(max(_ACTIVE_CTX.n_max, n_needed), v, rho_0, phi_0, z_0)
    return _ACTIVE_CTX


def g_mnTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
    """
    TM beam-shape coefficient g_{mn}^{TM} for Frozen Waves.

    Signature is preserved; alpha1, lamb, x_0, y_0 are unused in this FW version.
    """
    if abs(int(m)) > int(n):
        return 0.0 + 0.0j  # modo inexistente (|m| > n): BSC = 0
    ctx = _get_ctx(n, v, rho_0, phi_0, z_0)
    key = (int(m), int(n))
    if key in ctx.cache_tm:
        return ctx.cache_tm[key]

    m_int = int(m)
    n_int = int(n)

    a1 = prefator_TM(m_int)      # (1/2) i^(m+1) (-1)^((m-|m|)/2)
    b = fator_nm(m_int, n_int)   # (n-m)! / (n+|m|)!

    # Precompute arrays over q
    beta_fac = ctx.k_over_beta  # k/beta_q

    tau = ctx.tau_array(m_int, n_int)
    pi = ctx.pi_array(m_int, n_int)

    # Bessel parts
    order_d = m_int - ctx.v - 1
    order_f = m_int - ctx.v + 1

    d = ctx.J[order_d + ctx.order0] * np.exp(-1j * (m_int - ctx.v - 1) * ctx.phi_0)
    f = ctx.J[order_f + ctx.order0] * np.exp(-1j * (m_int - ctx.v + 1) * ctx.phi_0)

    e = (beta_fac * tau) + (m_int * pi)
    g = (beta_fac * (-tau)) + (m_int * pi)

    total = np.sum(ctx.kernel * ((d * e) + (f * g)))

    res = (a1 * b) * total
    ctx.cache_tm[key] = res
    return res


def g_mnTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
    """
    TE beam-shape coefficient g_{mn}^{TE} for Frozen Waves.

    Signature is preserved; alpha1, lamb, x_0, y_0 are unused in this FW version.
    """
    if abs(int(m)) > int(n):
        return 0.0 + 0.0j  # modo inexistente (|m| > n): BSC = 0
    ctx = _get_ctx(n, v, rho_0, phi_0, z_0)
    key = (int(m), int(n))
    if key in ctx.cache_te:
        return ctx.cache_te[key]

    m_int = int(m)
    n_int = int(n)

    a1 = prefator_TE(m_int)      # -(i/2) i^(m+1) (-1)^((m-|m|)/2)
    b = fator_nm(m_int, n_int)   # (n-m)! / (n+|m|)!

    beta_fac = ctx.k_over_beta

    tau = ctx.tau_array(m_int, n_int)
    pi = ctx.pi_array(m_int, n_int)

    order_d = m_int - ctx.v - 1
    order_f = m_int - ctx.v + 1

    d = ctx.J[order_d + ctx.order0] * np.exp(-1j * (m_int - ctx.v - 1) * ctx.phi_0)
    f = ctx.J[order_f + ctx.order0] * np.exp(-1j * (m_int - ctx.v + 1) * ctx.phi_0)

    # Eq. (19): colchetes [ m*(k/beta_q)*pi + tau ] e [ m*(k/beta_q)*pi - tau ].
    # CORRECAO B4 (2026-09): o termo 'g' tinha o fator k/beta_q aplicado a tau
    # em vez de a m*pi -- isto e, usava (m*pi - tau*k/beta) no lugar de
    # (m*(k/beta)*pi - tau). Note que TE difere mesmo de TM (Eq. 18), onde o
    # k/beta multiplica tau; a assimetria entre as duas e real, mas dentro do
    # TE os dois colchetes devem ter o k/beta no MESMO termo (m*pi).
    e = beta_fac * (m_int * pi) + tau
    g = beta_fac * (m_int * pi) - tau

    total = np.sum(ctx.kernel * ((d * e) - (f * g)))

    res = (a1 * b) * total
    ctx.cache_te[key] = res
    return res


def gn_FrozenWave(n, N, k_in, L_in, Q_in, z0):
    """
    Kept for compatibility with plotting scripts.
    Uses the same definition as the original implementation.
    """
    n = int(n)
    N = int(N)

    primeiroTermoDem = n * (n + 1)
    q_arr = np.arange(N, (-N + 1), dtype=np.int64)

    beta = float(Q_in) + (2.0 * np.pi * q_arr / float(L_in))
    k_termo = beta / float(k_in)

    # Aq for this helper should follow the "Aq(q, L)" definition (legacy profile),
    # but it is rarely used in the hot path. Keep a direct import to avoid cycles.
    from Aq import Aq
    Aq_vals = np.asarray([Aq(int(q), L_in) for q in q_arr], dtype=np.complex128)

    # pi_mn and tau_mn are imported from legendre.py, which preserves original semantics.
    from legendre import pi_mn, tau_mn
    termoMul = np.asarray([pi_mn(1, n, kt) + tau_mn(1, n, kt) for kt in k_termo], dtype=np.complex128)

    total = np.sum((Aq_vals / (1.0 + k_termo)) * termoMul * np.exp(1j * beta * float(z0)))

    return (-2.0 / primeiroTermoDem) * total
