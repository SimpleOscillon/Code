import numpy as np
from perturbative_solver import solve_oscillon
from matplotlib import pyplot as plt
from progress.bar import Bar

# edit these parameters:
w_range = np.linspace(0.5, 0.99, 50)
coeffs = np.array([1.25, -0.25])


def calculate_lifecycle(w_range, coeffs, N_harmonics=3, verbose=False):
    power_range = np.empty_like(w_range)
    energy_range = np.empty_like(w_range)
    with Bar('Processing', max=len(w_range)) as bar:
        for i, w in enumerate(w_range):
            R, S1, c_harmonics, S_harmonics, power, energy = solve_oscillon(
                w, coeffs=coeffs, N_harmonics=N_harmonics)
            power_range[i] = power
            energy_range[i] = energy
            bar.next()
    bar.finish()
    lifetime = -(np.diff(energy_range)[np.diff(energy_range) < 0] /
                 power_range[1:][np.diff(energy_range) < 0]).sum()
    return np.log10(lifetime), power_range, energy_range


if __name__ == '__main__':
    if coeffs.sum() != 1.0:
        coeffs = np.hstack((coeffs, [1.0 - coeffs.sum()]))
    log10lifetime, power_curve, energy_curve = calculate_lifecycle(
        w_range, coeffs)
    print('log10(lifetime)=', log10lifetime)
    plt.plot(w_range, power_curve, 'b-', lw=2.0)
    plt.xlabel('Frequency (m)', fontsize=14)
    plt.ylabel('Power (m^2)', fontsize=14)
    plt.yscale('log')
    plt.show()
