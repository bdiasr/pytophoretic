# -*- coding: utf-8 -*-
"""
Selecao do feixe incidente (fonte dos BSCs) para o vetor de assimetria.

Qualquer objeto com os metodos g_mnTM / g_mnTE na assinatura legada
    g_mn(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0)
serve como `beam=` em asymetryVector.I_x / I_y / I_z.

  beam=None               -> Frozen Wave linear em x (beamShape.py), EXATAMENTE
                             o caminho de sempre -- scripts antigos inalterados
  beam=BesselBeam(...)    -> feixe de Bessel circularmente simetrico (TCC)
"""
import beamShape
from besselBeam import BesselBeam  # noqa: F401  (reexportado por conveniencia)


def resolver_bsc(beam):
    """Devolve o par (g_TM, g_TE) de funcoes a usar."""
    if beam is None:
        return beamShape.g_mnTM, beamShape.g_mnTE
    return beam.g_mnTM, beam.g_mnTE
