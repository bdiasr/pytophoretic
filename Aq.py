import numpy as np

from cte import L, qs
from Fz import Fz, F_z5

# -----------------------------------------------------------------------------
# Fast Aq computation
# -----------------------------------------------------------------------------
# The original implementation computed Aq via torchquad Simpson integration at
# every call (very expensive). Here we keep the same mathematical definition,
# but compute Aq on a fixed, uniform grid and cache results for all q in cte.qs.
#
# NOTE: The codebase mainly uses Aq_2(q) (profile Fz). We keep Aq and Aq_1 for
# compatibility.
# -----------------------------------------------------------------------------

# Match the original Simpson setup: N=301 samples on [0, L]
_SIMpson_N = 301
_Z = np.linspace(0.0, float(L), _SIMpson_N, dtype=np.float64)
_DZ = _Z[1] - _Z[0]


def _simpson_uniform(y: np.ndarray, dz: float) -> complex:
    """
    Composite Simpson rule for an odd number of samples (N=2k+1) on a uniform grid.
    Returns a float/complex scalar.
    """
    n = y.shape[0]
    if n < 3 or (n % 2) == 0:
        # Fallback to trapezoid if something unexpected happens
        return np.trapz(y, dx=dz)
    s = y[0] + y[-1] + 4.0 * np.sum(y[1:-1:2]) + 2.0 * np.sum(y[2:-2:2])
    return (dz / 3.0) * s


# Cached tables (lazy-initialised)
_AQ2_TABLE = None  # profile Fz
_AQ1_TABLE = None  # profile F_z5
_AQ_TABLE = None   # Aq(q, L) legacy -- agora aponta para Fz (Funcao 2)


def _build_table(profile: str) -> dict:
    """
    Build Aq table for all integer q in cte.qs.
    """
    frac = 1.0 / float(L)
    table = {}

    if profile == "Fz":
        F = Fz(_Z).astype(np.complex128, copy=False)
    elif profile == "F_z5":
        F = F_z5(_Z).astype(np.complex128, copy=False)
    else:
        raise ValueError(f"Unknown profile: {profile}")

    # Precompute exp(i 2π q z / L) for all q with broadcasting
    q_arr = np.asarray(qs, dtype=np.int64)
    phase = np.exp(1j * (2.0 * np.pi) * (q_arr[:, None] * _Z[None, :] / float(L))).astype(np.complex128)

    # Integrate for each q using Simpson (same N as original)
    # Integral = (1/L) ∫ F(z) exp(i 2π q z / L) dz
    integrand = phase * F[None, :]
    for i, q in enumerate(q_arr):
        table[int(q)] = frac * _simpson_uniform(integrand[i], _DZ)

    return table


def Aq(q, L_in=None):
    """
    Aq legada (assinatura Aq(q, L); L_in e' ignorado, usa-se cte.L).

    MUDANCA (2026-09): antes usava o perfil F_z6 (dois degraus binarios), que
    foi removido do projeto. Agora usa `Fz` (Funcao 2) -- o MESMO perfil de
    Aq_2 e dos BSCs. Isso tambem corrige a inconsistencia em que
    frozenWave.Psi() (que chama esta funcao) desenhava um feixe diferente
    daquele que gerava as forcas.
    """
    global _AQ_TABLE
    if _AQ_TABLE is None:
        _AQ_TABLE = _build_table("Fz")
    return _AQ_TABLE[int(q)]


def Aq_1(q):
    """
    Aq_1 uses profile F_z5.
    """
    global _AQ1_TABLE
    if _AQ1_TABLE is None:
        _AQ1_TABLE = _build_table("F_z5")
    return _AQ1_TABLE[int(q)]


# -----------------------------------------------------------------------------
# Perfil ATIVO dos BSCs
# -----------------------------------------------------------------------------
# Aq_2 alimenta g_mnTM/g_mnTE (beamShape) e portanto todo o pipeline de forcas.
# Antes o perfil estava fixo em Fz (Funcao 2). Agora e' selecionavel, para que
# se possa simular tambem com F_z5 (Funcao 1) sem editar o arquivo.
#
#     set_perfil_ativo('F1')  ou  set_perfil_ativo('F2')
#
# Trocar o perfil invalida os caches de Aq_2 E o contexto de beamShape (que
# guarda o kernel A_q * exp(i beta z0)) -- ambos sao limpos aqui.
# -----------------------------------------------------------------------------
_PERFIL_ATIVO = 'F2'          # 'F1' -> F_z5 ; 'F2' -> Fz
_PERFIL_MAP = {'F1': 'F_z5', 'F2': 'Fz'}


def set_perfil_ativo(nome):
    """Define o perfil que alimenta os BSCs ('F1' ou 'F2')."""
    global _PERFIL_ATIVO, _AQ2_TABLE
    nome = str(nome).upper()
    if nome not in _PERFIL_MAP:
        raise ValueError(f"perfil '{nome}' invalido; use 'F1' ou 'F2'")
    if nome != _PERFIL_ATIVO:
        _PERFIL_ATIVO = nome
        _AQ2_TABLE = None                     # invalida a tabela de A_q
        import beamShape as _bs                # invalida o contexto dos BSCs
        _bs._ACTIVE_CTX = None
    return _PERFIL_ATIVO


def get_perfil_ativo():
    return _PERFIL_ATIVO


def Aq_2(q):
    """
    A_q do perfil ATIVO (default F2 = Fz). Alimenta g_mnTM/g_mnTE.
    """
    global _AQ2_TABLE
    if _AQ2_TABLE is None:
        _AQ2_TABLE = _build_table(_PERFIL_MAP[_PERFIL_ATIVO])
    return _AQ2_TABLE[int(q)]


def Aq_2_vec(q_values):
    """
    Vectorised lookup for Aq_2 on an iterable/array of integer q.
    """
    global _AQ2_TABLE
    if _AQ2_TABLE is None:
        _AQ2_TABLE = _build_table(_PERFIL_MAP[_PERFIL_ATIVO])
    q_arr = np.asarray(q_values, dtype=np.int64)
    return np.asarray([_AQ2_TABLE[int(q)] for q in q_arr], dtype=np.complex128)
