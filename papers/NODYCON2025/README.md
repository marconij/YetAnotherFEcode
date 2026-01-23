# Backbone optimization of a MEMS gyro using an invariant manifold-based reduced order model

This code provides a MATLAB implementation of the numerical example presented at the Fourth International Nonlinear Dynamics Conference (NODYCON 2025).

## Installation

To run the code, clone the [YetAnotherFEcode](https://github.com/jain-shobhit/YetAnotherFEcode) repository and run the `startup.m` script to add the necessary paths.

## Examples

In this example, we tune the nonlinear coefficient of the reduced dynamics on a 3rd-order expansion of the Lyapunov Subcenter Manifold (LSM) using parametric optimization (fmincon). Sensitivities are computed efficiently using the adjoint method.

The system we consider is a MEMS gyroscope prototype (named TS4), which is modelled using the Multi-Point Constraint (MPC) method: masses are represented by rigid bodies, connected to beams with a master-slave approach. In this demo, we also tune the drive and sense vibration modes to be at specific target frequencies in order to test the mode-tracking algorithm.

To validate the results, the package [SSMTool 2.6](https://github.com/jain-shobhit/SSMTool) is required.

Two detailed demos on the optimization of the $\gamma$ coefficient are provided as MATLAB live scripts:

1. `demo_gamma_popt.mlx`: constraint on $\gamma \leq -1e10$.
2. `demo_gamma_popt_v2.mlx`: constraint on $\gamma \leq -1e11$.

## Citation
To cite this package, please cite the following:

M. Pozzi, J. Marconi, S. Jain, M. Li, and F. Braghin. "Backbone optimization of a MEMS gyro using an invariant manifold-based reduced order model". *Proceedings of the Fourth International Nonlinear Dynamics Conference* (2025, to appear).

## Contact
If you have any questions, create an issue or email <jacopo.marconi@polimi.it> or <matteo1.pozzi@polimi.it>.
