# Harmonic-decomposition-approach-to-dynamical-friction-for-eccentric-orbits
This repository contains code and results used in the paper "Harmonic-decomposition approach to dynamical friction for eccentric orbits"

The repository contains the simulations' measurements for the angular momentum and energy dissipation rates used in the paper, specifically for figures 4 and 5. The results are found in the files "angular_dissipation.csv" and "energy_dissipation.csv".

The code used to calculate the friction coefficients is also available and follows the analytic results of the paper. We provided the Julia file "hansen_calculator.jl" to calculate recursively the Hansen coefficients used to calculate the g^{(j)}_{lm} function. It creates a cache that contains the coefficients and is used by other files and another file with some Newcomb operators that are calculated along the way. 
