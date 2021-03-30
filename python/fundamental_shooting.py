import numpy as np
from scipy.special import j0 as BesselJ0, j1 as BesselJ1, jn as BesselJ
from scipy.optimize import root


def shoot_S1(central_value: float,
             w: float,
             R: np.ndarray,
             coeffs: np.ndarray,
             S_harmonics: np.ndarray = None) -> np.ndarray:
    """
    Shoots S1 from the center, starting from a central_value and
    zero-drivative.

    Parameters
    ----------
    central_value : the value of S1 at r=0
    w : the frequency of the oscillon
    R : the grid of radii
    coeffs : the Fourier coefficients of the potential, normalized
    S_harmonics : the in-phase perturbative radiative harmonics

    Returns
    -------
    S1 : the values of S1 over the grid when shooting from the center
    """

    S1 = np.empty_like(R)
    S1[0], S1[1] = central_value, central_value
    dr = R[1] - R[0]

    if S_harmonics is None:

        def f_(i):
            return (
                S1[i - 1] * w**2 - 2 *
                (coeffs * BesselJ1(S1[i - 1] * np.arange(1,
                                                         len(coeffs) + 1)) /
                 np.arange(1,
                           len(coeffs) + 1)).sum())
    else:
        N_harmonics = S_harmonics.shape[0]

        def f_(i):
            return (
                S1[i - 1] * w**2 - 2 *
                (coeffs * BesselJ1(S1[i - 1] * np.arange(1,
                                                         len(coeffs) + 1)) /
                 np.arange(1,
                           len(coeffs) + 1)).sum() -
                ((coeffs[:, np.newaxis] * BesselJ(
                    2 * np.arange(0, N_harmonics, 1) + 2, S1[i - 1] *
                    np.arange(1,
                              len(coeffs) + 1, 1)[:, np.newaxis])).sum(axis=0)
                 * S_harmonics[:, i - 1]).sum() +
                ((coeffs[:, np.newaxis] * BesselJ(
                    2 * np.arange(0, N_harmonics, 1) + 4, S1[i - 1] *
                    np.arange(1,
                              len(coeffs) + 1, 1)[:, np.newaxis])).sum(axis=0)
                 * S_harmonics[:, i - 1]).sum())

    for i in range(2, len(S1)):
        S1[i] = (
            2 * S1[i-1] -
            dr**2 * f_(i) +
            S1[i-2] * (2 * dr/(2*R[i-1]) - 1)) \
            / (2 * dr/(2*R[i-1]) + 1)
    return S1


def initial_S1(w: float,
               R: np.ndarray,
               coeffs: np.ndarray,
               S_harmonics: np.ndarray = None) -> np.ndarray:
    """
    Defines the binary search procedure to find the initial condition which
    shoots to zero at infinity.
    """
    # find the value at the same potential energy as the zero-field. We know
    # the true value will be slightly higher due to friction.
    c = root(
        lambda x: 0.5 * x**2 * w**2 + 2 * (coeffs * (BesselJ0(x * np.arange(
            1,
            len(coeffs) + 1)) - 1) / np.arange(1,
                                               len(coeffs) + 1)**2).sum(),
        10).x[0]

    # define the left- and right-boundaries of the search and push
    # these values apart until they have the appropriate signs:
    left, right = c, c
    left_condition = (shoot_S1(left, w, R, coeffs, S_harmonics)[-1] >= 0)
    while not left_condition:
        left = 0.95 * left
        left_condition = (shoot_S1(left, w, R, coeffs, S_harmonics)[-1] >= 0)

    right_condition = shoot_S1(right, w, R, coeffs, S_harmonics)[-1] < 0
    while not right_condition:
        right = 1.1 * right
        right_condition = shoot_S1(right, w, R, coeffs, S_harmonics)[-1] < 0

    # perform the binary search for 60 steps:
    for _ in range(60):
        m = (left + right) / 2
        S1 = shoot_S1(m, w, R, coeffs, S_harmonics)
        if S1[-1] >= 0:
            left = m
        else:
            right = m

    # zero-out the far-field:
    S1[np.abs(S1).argmin():] = 0.0
    return S1
