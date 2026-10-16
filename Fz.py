import numpy as np
import scipy.special as sp

from cte import L


def _asarray(z):
    return np.asarray(z, dtype=np.float64)


def Fz(z):
    """
    Default Frozen Wave longitudinal profile used by Aq_2.

    Original implementation used torch operations. This numpy version is
    equivalent and much cheaper to import/execute.
    """
    z = _asarray(z)
    Lf = float(L)

    condition = (z >= Lf / 12.0) & (z <= 11.0 * Lf / 12.0)
    expoente = -5.0 * ((z - 0.5 * Lf) ** 2) / (Lf ** 2)

    out = np.zeros_like(z, dtype=np.float64)
    out[condition] = np.exp(expoente[condition]) * np.cos(6.0 * np.pi * z[condition] / Lf)
    return out


def F_z1(z, L_in):
    z = _asarray(z)
    Lf = float(L_in)
    condition = (z >= 3.0 * Lf / 8.0) & (z <= 5.0 * Lf / 8.0)
    out = np.zeros_like(z, dtype=np.float64)
    out[condition] = 1.0
    return out


def F_z1_L(z):
    """F_z1 com L fixado em cte.L (assinatura de 1 argumento, para o registro)."""
    return F_z1(z, L)


def F_z3(z):
    # Same form as Fz in the original file
    return Fz(z)


def F_z4(z):
    z = _asarray(z)
    Lf = float(L)
    condition = (z >= 0.25 * Lf) & (z <= 0.75 * Lf)
    out = np.zeros_like(z, dtype=np.float64)
    out[condition] = sp.j0((1.6e-6) * z[condition])
    return out


def F_z5(z):
    """
    Piecewise polynomial profile used by Aq_1.
    """
    z = _asarray(z)
    Lf = float(L)

    l1 = Lf / 10.0
    l2 = 3.0 * Lf / 10.0
    l3 = 4.0 * Lf / 10.0
    l4 = 6.0 * Lf / 10.0
    l5 = 7.0 * Lf / 10.0
    l6 = 9.0 * Lf / 10.0

    out = np.zeros_like(z, dtype=np.float64)

    cond1 = (z >= l1) & (z <= l2)
    cond2 = (z >= l3) & (z <= l4)
    cond3 = (z >= l5) & (z <= l6)

    out[cond1] = (-4.0 * (z[cond1] - l1) * (z[cond1] - l2)) / ((l2 - l1) ** 2)
    out[cond2] = 1.0
    out[cond3] = (-4.0 * (z[cond3] - l5) * (z[cond3] - l6)) / ((l6 - l5) ** 2)

    return out


# -----------------------------------------------------------------------------
# Registro de perfis longitudinais F(z)
# -----------------------------------------------------------------------------
# Mapeia um nome estavel -> (funcao, rotulo LaTeX para figuras).
#
# Usado por fw_field.py para calcular A_n (Eq. 2) e Psi (Eq. 1) para qualquer
# perfil, sem precisar de uma Aq_N/Psi_N dedicada por funcao.
#
# NOTA: F_z3 nao entra no registro -- e' `return Fz(z)`, duplicata exata da
# Funcao 2 (F2_gauss_cos).
#
# REMOVIDAS em 2026-09 (decisao do projeto, nao serao implementadas):
#   F_z2 / 'F_exp'    -- rampa exponencial
#   F_z6 / 'F_step2'  -- dois degraus binarios
# Backup do arquivo anterior: _old/Fz.py.bak-pre-limpeza
#
# Padrao do projeto (definido em 2026-09):
#   Funcao 1 = F1_poly        (polinomial por partes)
#   Funcao 2 = F2_gauss_cos   (gaussiana modulada por cosseno)
#   Funcao 3 = a definir
# -----------------------------------------------------------------------------

PROFILES = {
    # --- padroes do projeto (2026-09) ---
    'F1_poly':      (F_z5, 'Funcao 1: polinomial por partes'),
    'F2_gauss_cos': (Fz,   r'Funcao 2: $e^{-5(z-L/2)^2/L^2}\cos(6\pi z/L)$'),
    # --- auxiliares, nao sao padrao do projeto ---
    'F_box':        (F_z1_L, r'Degrau unico: $[3L/8,5L/8]$'),
    'F_bessel':     (F_z4, r'$J_0(1.6\times10^{-6}z)$ em $[L/4,3L/4]$'),
}

# Apelidos curtos: "F1" e "F2" referem-se sempre aos dois padroes acima.
F1 = 'F1_poly'
F2 = 'F2_gauss_cos'

def get_profile(name):
    """Devolve (funcao, rotulo) do perfil registrado."""
    if name not in PROFILES:
        raise KeyError(f"perfil '{name}' nao registrado. Disponiveis: {list(PROFILES)}")
    return PROFILES[name]


def list_profiles():
    """Nomes dos perfis registrados, em ordem."""
    return list(PROFILES)
