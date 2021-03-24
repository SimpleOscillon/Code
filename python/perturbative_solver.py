import numpy as np
from scipy import sparse
from scipy.special import jn as BesselJ

from fundamental_shooting import initial_S1


def bessel_wrap(coeffs: np.ndarray, indices: np.ndarray, arguments: np.ndarray,
                indices_signs: np.ndarray) -> np.ndarray:
    """
    Calculates sums of Bessel functions of given indices, over potential
    coefficients.
    """
    return ((coeffs / np.arange(1,
                                len(coeffs) + 1))[:, np.newaxis] *
            (indices_signs[:, np.newaxis, np.newaxis] * BesselJ(
                np.abs(indices)[:, np.newaxis, np.newaxis],
                np.arange(1,
                          len(coeffs) + 1)[np.newaxis, :, np.newaxis] *
                arguments[np.newaxis, np.newaxis, :])).sum(axis=0)).sum(axis=0)


def make_c_block(k: int, l: int, coeffs: np.ndarray, w: float, R: np.ndarray,
                 S1: np.ndarray) -> sparse.spmatrix:
    """
    Calculate the single-harmonic block between harmonics c_k and c_l.
    """
    N = len(R)
    dr = R[1] - R[0]
    if k == l:
        main_diagonal = -2 / dr**2 + (k * w)**2 - bessel_wrap(
            coeffs=coeffs,
            indices=np.array([0, 2 * k]),
            arguments=S1,
            indices_signs=np.array([1.0, 1.0]))
        main_diagonal[-1] += 1 / dr**2
        return sparse.diags([
            np.ones((N - 1)) / dr**2, main_diagonal,
            np.ones((N - 1)) / dr**2
        ], [-1, 0, 1])
    else:
        return sparse.diags([
            -bessel_wrap(coeffs=coeffs,
                         indices=np.array([k + l, np.abs(k - l)]),
                         arguments=S1,
                         indices_signs=np.array([1.0, 1.0]))
        ], [0])


def make_S_block(k: int, l: int, coeffs: np.ndarray, w: float, R: np.ndarray,
                 S1: np.ndarray) -> sparse.spmatrix:
    """
    Calculate the single-harmonic block between harmonics S_k and S_l.
    """
    N = len(R)
    dr = R[1] - R[0]
    if k == l:
        main_diagonal = -2 / dr**2 + (k * w)**2 - bessel_wrap(
            coeffs=coeffs,
            indices=np.array([0, 2 * k]),
            arguments=S1,
            indices_signs=np.array([1.0, -1.0]))
        main_diagonal[-1] += 1 / dr**2
        return sparse.diags([
            np.ones((N - 1)) / dr**2, main_diagonal,
            np.ones((N - 1)) / dr**2
        ], [-1, 0, 1])
    else:
        return sparse.diags([
            +bessel_wrap(coeffs=coeffs,
                         indices=np.array([k + l, np.abs(k - l)]),
                         arguments=S1,
                         indices_signs=np.array([1.0, -1.0]))
        ], [0])


def make_boundary_block(k: int, l: int, w: float,
                        R: np.ndarray) -> sparse.spmatrix:
    """
    Connection block between c_k and S_l. Its negative is the connection block
    between S_k and c_l.
    """
    N = len(R)
    dr = R[1] - R[0]
    if k == l:
        return sparse.vstack([
            sparse.csr_matrix((N - 1, N), dtype=np.float),
            np.concatenate([np.zeros(N - 1), [-np.sqrt((k * w)**2 - 1) / dr]])
        ])
    else:
        return sparse.csr_matrix((N, N), dtype=np.float)


def solve_oscillon(w: float,
                   coeffs: np.ndarray = np.array([1.0]),
                   N_harmonics: int = 3,
                   dr: float = 0.01,
                   L: float = 20.0):
    """
    Main function analyzing the oscillon at a given frequency.

    Parameters
    ----------
    w : float
        Frequency of the oscillon.
    coeffs : np.ndarray
        Fourier coefficients of the potential. Must sum to one.
    N_harmonics : int, default value 3
        Number of perturbative harmonics to consider.
    dr : float, default value 0.01
        Resolution of the spatial grid.
    L : float, default value 20.0
        The size of the box.

    Returns
    -------
    R : np.ndarray of size (N,)
        The equally-spaced grid of spatial points between 0 and L.
    S1 : np.ndarray of size (N,)
        The fundamental, non-perturbative harmonic obtained by
        shooting.
    c_harmonics : np.ndarray of size (N_harmonics, N)
        c_harmonics[i, :] is the orthogonal deformation for harmonic
        at frequency (2 * i + 3) * w.
    S_harmonics : np.ndarray of size (N_harmonics, N)
        S_harmonics[i, :] is the quasibreather harmonic
        at frequency (2 * i + 3) * w.
    power : float
        The radiated power in the leading N_harmonics radiative harmonics.
    energy : float
        The energy in the non-perturbative fundamental harmonic.
    """

    R = np.arange(dr, L, dr)
    N = len(R)

    S1 = initial_S1(w, R, coeffs)

    c_megablock = sparse.vstack([
        sparse.hstack([
            make_c_block(2 * k + 3, 2 * l + 3, coeffs, w, R, S1)
            for l in range(N_harmonics)
        ]) for k in range(N_harmonics)
    ])

    cS_megablock = sparse.vstack([
        sparse.hstack([
            make_boundary_block(2 * k + 3, 2 * l + 3, w, R)
            for l in range(N_harmonics)
        ]) for k in range(N_harmonics)
    ])

    Sc_megablock = sparse.vstack([
        sparse.hstack([
            -make_boundary_block(2 * k + 3, 2 * l + 3, w, R)
            for l in range(N_harmonics)
        ]) for k in range(N_harmonics)
    ])

    S_megablock = sparse.vstack([
        sparse.hstack([
            make_S_block(2 * k + 3, 2 * l + 3, coeffs, w, R, S1)
            for l in range(N_harmonics)
        ]) for k in range(N_harmonics)
    ])

    A = sparse.csr_matrix(
        sparse.vstack([
            sparse.hstack([c_megablock, cS_megablock]),
            sparse.hstack([Sc_megablock, S_megablock])
        ]))
    b = np.hstack([np.zeros(N * N_harmonics)] + [
        2 * R * bessel_wrap(coeffs=coeffs,
                            indices=np.array([2 * k + 3]),
                            arguments=S1,
                            indices_signs=np.array([1.0]))
        for k in range(N_harmonics)
    ])

    x = sparse.linalg.spsolve(A, b).reshape([2 * N_harmonics, N])
    x = x / R

    c_harmonics = x[:N_harmonics]
    S_harmonics = x[N_harmonics:]

    power = (2 * np.pi * L**2 *
             (c_harmonics[:, -1]**2 + S_harmonics[:, -1]**2) *
             np.array([(2 * k + 3) * w * np.sqrt(((2 * k + 3) * w)**2 - 1)
                       for k in range(N_harmonics)])).sum()

    energy = 4 * np.pi * (
        R[1:-1]**2 *
        (0.5 * w**2 * S1[1:-1]**2 + 0.5 * (S1[2:] - S1[:-2])**2 /
         (2 * dr)**2 - bessel_wrap(coeffs=coeffs,
                                   indices=np.array([0]),
                                   arguments=S1[1:-1],
                                   indices_signs=np.array([1.0])) +
         (coeffs / np.arange(1,
                             len(coeffs) + 1)**2).sum()) * dr).sum()

    return R, S1, c_harmonics, S_harmonics, power, energy
