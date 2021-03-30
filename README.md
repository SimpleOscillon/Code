# Physical Quasibreather Oscillon Calculator

Code to compute the properties of the oscillon in the PQB formalism. 
The folders Matlab and Python contain equivalent versions of the procedure outlined in the paper arXiv:.

The procedure allows you to obtain the oscillon's power curves, energy, and integrated lifetime at desired resolution, for a scalar field potential provided in terms of a Fourier expansion. The calculations are restricted to work with one non-perturbative harmonic in order to maximize ease-of-use and speed. We provide a parameter which controls the number of times the code is iterated to account for linear back-reaction: if one needs to increase this number beyond 3 for good convergence, it is highly likely that the higher harmonics have become non-perturbatively large. In order to do more precise calculations in the PQB formalism taking into account the non-perturbativity of higher harmonics, the user is encouraged to read appendix A of arXiv:.
