import numpy as np
import scipy.special as s


def js_n(n, x):
    """
    Spherical Bessel function of the first kind j_n(x).
    """
    return s.spherical_jn(n, x, False)


def js_nz(n, z):
    if n < 0:
        return 0
    # Keep the original sign-handling for negative real z (rare in this project)
    if (not isinstance(z, complex)) and (z < 0):
        return (-1) ** n * s.spherical_jn(n, -z, False)
    return s.spherical_jn(n, z, False)


def j0_nz(n, z):
    return s.jv(n, z, out=None)


def ys_n(n, x):
    """
    Spherical Bessel function of the second kind y_n(x).
    """
    return s.spherical_yn(n, x, False)


def psiBessel(n, x):
    """
    Riccati-Bessel ψ_n(x) = x * j_n(x)
    """
    return x * s.spherical_jn(n, x, False)


# Riccati-Bessel derivative: ψ'_n(x) = (n+1) j_n(x) - x j_{n+1}(x)
derivativePsiBesseln = (lambda n, x: ((1 + n) * s.spherical_jn(n, x, False) - (x * s.spherical_jn((n + 1), x, False))))


def sphericalHankel_n(n, x):
    """
    Spherical Hankel of the second kind: h_n^(2)(x) = j_n(x) - i y_n(x)
    """
    return s.spherical_jn(n, x, False) - 1j * s.spherical_yn(n, x, False)


def RiccatiBessel_2Type(n, x):
    """
    Riccati-Bessel of the second kind: ξ_n(x) = x * h_n^(2)(x)
    """
    return x * sphericalHankel_n(n, x)


def derivativeSphericalHankel_n(n, x):
    """
    Derivative of x*h_n^(2)(x): ξ'_n(x) = (n+1)h_n^(2)(x) - x h_{n+1}^(2)(x)
    """
    return (1 + n) * sphericalHankel_n(n, x) - x * sphericalHankel_n((n + 1), x)
