# Simple Oscillon

Software to calculate oscillon properties in the PQB framework. Requires Matlab to run.

The high-level interface can be accessed by running the code in `PublicPerturbativePeriodic.m`, where an oscillon model is defined by the choice of:
- `Vcoeff` is the list of Fourier coefficients of the potential. _Additional term will be automatically added to enforce sum to be `omegaMax^2`._ Do not separate inputs by commas as per Matlab syntax.
- `thetaMax` (default value 1) is the fundamental periodicity of the potential.
- `Radius` (default value 15) is the maximum radius of the simulation.
- `dr` (default value 0.01) is the spatial resolution.
- `NHarmonics` (default value 2) is the number of perturbative harmonics to consider. _In this code, only the fundamental is treated non-perturbatively._
- `OmegaList` (default range `0.80:0.01:0.94`) determines the values of the frequency at which the oscillon is computed.

This code outputs:
- oscillon lifetime estimate in the specificed frequency range, expressed as a logarithm in base-10.
- the plot of log-base-10 power radiated versus frequency.
- other outputs can be accessed through the function `PublicPowerCurve.m` in order 
  - `[PowerVsOmegaList,EnergyVsOmegaList,Lifetime,PowerInHarmonics,SList,CList,r]`
  - The list of radiated power and frequency, labeled `PowerVsOmegaList`
  - The list of bound oscillon energies and frequency, labeled `EnergyVsOmegaList`
  - The estimated lifetime `Lifetime`
  - The power stored in each harmonic `PowerInHarmonics`
  - The harmonic decomposition of the PQB `SList`
  - The harmonic decomposition of the OD `CList`
  - The position coordinates `r` corresponding to the entries of `SList` and `CList`
