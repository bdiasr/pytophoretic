# -*- coding: utf-8 -*-
"""Core modules of pytophoretic, grouped by stage of the calculation.

    config/     cte, Fz, Aq                       (constants, F(z) profiles, A_q)
    math/       specials, bsc_common, legendre    (Riccati-Bessel, P_n^m, pi, tau)
    beams/      beamShape, fw_field, frozenWaveX, besselBeam, besselBeamX,
                besselBeamRA, planeWave, beams    (BSCs of each beam)
    glmt/       coefficients, campo_interno       (Mie coefficients, internal field)
    asymmetry/  asymetryVector, asymetryVector_pereira, asymetryVector_antigo
    force/      fotoforese_pereira                (photophoretic force)

`import core` only adds these sub-folders to sys.path, so every module keeps its
short import name (import cte, import beamShape, ...). It does NOT import any
module: scripts that redefine cte before importing the beam modules
(simulations/fw1064_comum.py, simulations/ipc2024.py) keep working.
"""
import os as _os
import sys as _sys

_BASE = _os.path.dirname(_os.path.abspath(__file__))
SUBFOLDERS = ('config', 'math', 'beams', 'glmt', 'asymmetry', 'force')
for _sub in reversed(SUBFOLDERS):
    _p = _os.path.join(_BASE, _sub)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)
