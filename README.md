# LAPIS
LAPIS is an approach that reconstructs landscapes by integrating single-snapshot data with stochastic dynamical modeling.
The specific steps for running LAPIS are as follows:
(1) Run first_screening.m to identify better candidate initial parameters for subsequent optimization.
(2) Run second_refine.m to determine the optimal model parameters and diffusion coefficient, aiming to achieve the best possible fit to the data. LAPIS-optimized parameters will then serve as the basis for subsequent landscape reconstruction and sensitivity analysis.
