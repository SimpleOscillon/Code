# Physical Quasibreather Oscillon Calculator

Code to compute the properties of the oscillon in the PQB formalism. 
The folders Matlab and Python contain equivalent versions of the procedure outlined in the paper _The Structure of the Oscillon: The Dynamics of Attractive Self-Interaction_, available at [arXiv:2104.02069](https://arxiv.org/abs/2104.02069).

The procedure allows you to obtain the oscillon's power curves, energy, and integrated lifetime at desired resolution for a scalar field potential provided in terms of a Fourier expansion:

<img width="336" alt="Screen Shot 2021-04-01 at 9 40 35 AM" src="https://user-images.githubusercontent.com/53380799/113302674-514f4580-92ce-11eb-8f35-ce6b13e2553f.png">

The calculations are restricted to work with one non-perturbative harmonic in order to maximize ease-of-use and speed. We provide a parameter which controls the number of times the code is iterated to account for linear back-reaction: if one needs to increase this number beyond 3 for good convergence, it is highly likely that the higher harmonics have become non-perturbatively large. In order to do more precise calculations in the PQB formalism taking into account the non-perturbativity of higher harmonics, the user is encouraged to read appendix A of [arXiv:2104.02069](https://arxiv.org/abs/2104.02069).
