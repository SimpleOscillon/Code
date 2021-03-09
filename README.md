# Simple Oscillon

Software to calculate oscillon properties in the PQB framework. Requires Matlab to run.

The high-level interface can be accessed by running the code in `PublicPerturbativePeriodic.m`, where an oscillon model is defined by the choice of:
- `Vcoeff` as the list of Fourier coefficients of the potential. _Additional term will be automatically added to enforce sum to be `omegaMax&2`._ Do not separate inputs by commas as per Matlab syntax.
- `omegaMax` (default value 1) is the fundamental periodicity of the potential. Higher values means higher resolution.
- `Radius` (default value 15) is the maximum radius of the simulation.
- `dr` (default value 0.01) is the spatial resolution.
- `NHarmonics` (default value 2) is the number of perturbative harmonics to consider. _In this code, only the fundamental is treated non-perturbatively._
- `OmegaList` (default range `0.80:0.01:0.94`) determines the values of the frequency at which the oscillon is computed.

This code outputs:
- oscillon lifetime in the specificed frequency range, expressed as a logarithm in base-10.
- the plot of power radiated versus frequency.
