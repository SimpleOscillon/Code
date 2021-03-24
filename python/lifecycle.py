import numpy as np
from perturbative_solver import solve_oscillon
from matplotlib import pyplot as plt
from progress.bar import Bar

############################################################################
# Edit these parameters:
############################################################################

# the values of the frequency to consider:
w_range = np.linspace(0.7, 0.99, 100)
# the Fourier coefficients of the potential. If they do not sum to one,
# another one will be added to satisfy the sum:
coeffs = np.array([0.7, -0.5, 0.5])
# the size of the spatial box:
L = 20.0
# the spatial step size:
dr = 0.01
# number of perturbative harmonics:
N_harmonics = 3

############################################################################
# Compute power curve and lifetime:
############################################################################


def calculate_lifecycle(w_range, coeffs, N_harmonics=3):
    """
    Auxiliary function to compute lifetime over a range of frequencies.
    """
    power_range = np.empty_like(w_range)
    energy_range = np.empty_like(w_range)
    # iterate through frequencies and collect power and energy information:
    with Bar('Processing', max=len(w_range)) as bar:
        for i, w in enumerate(w_range):
            R, S1, c_harmonics, S_harmonics, power, energy = solve_oscillon(
                w, coeffs=coeffs, N_harmonics=N_harmonics, dr=dr, L=L)
            power_range[i] = power
            energy_range[i] = energy
            bar.next()
    bar.finish()
    # lifetime is only integrated over segments of decreasing energy:
    lifetime = -(np.diff(energy_range)[np.diff(energy_range) < 0] /
                 power_range[1:][np.diff(energy_range) < 0]).sum()
    return np.log10(lifetime), power_range, energy_range


if __name__ == '__main__':
    # add the coefficient to satisfy the sum-to-one criterion, if needed:
    if coeffs.sum() != 1.0:
        coeffs = np.hstack((coeffs, [1.0 - coeffs.sum()]))

    log10lifetime, power_curve, energy_curve = calculate_lifecycle(
        w_range, coeffs)

    print('log10(lifetime)=', log10lifetime)

    # plot decreasing-energy and increasing-energy segments separately:
    for i in range(len(power_curve) - 1):
        if energy_curve[i + 1] - energy_curve[i] <= 0:
            plt.plot(w_range[[i, i + 1]],
                     power_curve[[i, i + 1]],
                     'b-',
                     lw=2.0)
        else:
            plt.plot(w_range[[i, i + 1]],
                     power_curve[[i, i + 1]],
                     'r--',
                     lw=1.0,
                     alpha=0.5)
    plt.xlabel('Frequency (m)', fontsize=14)
    plt.ylabel(r'Power ($f^2$)', fontsize=14)
    plt.yscale('log')
    plt.show()
