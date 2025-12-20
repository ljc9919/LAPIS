# Landscape reconstruction And Parameter Inference from Single-snapshot data (LAPIS)
<img width="11426" height="7460" alt="Image" src="https://github.com/user-attachments/assets/89ec1978-26d9-46c0-a238-3e12d32f73de" />
LAPIS is a landscape reconstruction approach by integrating single-snapshot data with stochastic dynamical modeling.

The specific steps for running LAPIS are as follows:

(1) Run first_screening.m to identify better candidate initial parameters for subsequent optimization.

(2) Run second_refining.m to determine the optimal model parameters and diffusion coefficient, aiming to achieve the best possible fit to the snapshot data. LAPIS-optimized parameters will then serve as the basis for subsequent landscape reconstruction, sensitivity analysis, and landscape control.

