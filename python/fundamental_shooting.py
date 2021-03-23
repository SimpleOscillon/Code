import numpy as np
from scipy.special import j0 as BesselJ0, j1 as BesselJ1
from scipy.optimize import root


def shoot_S1(central_value: float, w: float, R: np.ndarray,
             coeffs: np.ndarray) -> np.ndarray:
    """
    Parameters
    ----------
    central_value : the value of S1 at r=0
    w : the frequency of the oscillon
    R : the grid of radii
    coeffs : the Fourier coefficients of the potential, normalized

    Returns
    -------
    S1 : the values of S1 over the grid when shooting from the center
    """
    def f_(s):
        return s * w**2 - 2 * (coeffs * BesselJ1(s * np.arange(
            1,
            len(coeffs) + 1)) / np.arange(1,
                                          len(coeffs) + 1)).sum()

    dr = R[1] - R[0]
    S1 = np.empty_like(R)
    S1[0], S1[1] = central_value, central_value
    for i in range(2, len(S1)):
        S1[i] = (
            2 * S1[i-1] -
            dr**2 * f_(S1[i-1]) +
            S1[i-2] * (2 * dr/(2*R[i-1]) - 1)) \
            / (2 * dr/(2*R[i-1]) + 1)
    return S1


def initial_S1(w: float, R: np.ndarray, coeffs: np.ndarray) -> np.ndarray:

    c = root(
        lambda x: 0.5 * x**2 * w**2 + 2 * (coeffs * (BesselJ0(x * np.arange(
            1,
            len(coeffs) + 1)) - 1) / np.arange(1,
                                               len(coeffs) + 1)**2).sum(),
        10).x[0]
    l, r = c, c
    left_condition = (shoot_S1(l, w, R, coeffs)[-1] >= 0)
    while not left_condition:
        l = 0.95 * l
        left_condition = (shoot_S1(l, w, R, coeffs)[-1] >= 0)

    right_condition = shoot_S1(r, w, R, coeffs)[-1] < 0
    while not right_condition:
        r = 1.5 * r
        right_condition = shoot_S1(r, w, R, coeffs)[-1] < 0

    for _ in range(60):
        m = (l + r) / 2
        S1 = shoot_S1(m, w, R, coeffs)
        if S1[-1] >= 0:
            l = m
        else:
            r = m
    S1[np.abs(S1).argmin():] = 0.0
    return S1
