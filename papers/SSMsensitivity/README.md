# Sensitivity analysis of arbitrary order Spectral Submanifolds (SSM) for geometrically nonlinear mechanical systems

This repository provides a MATLAB implementation for computing the damped backbone curve of geometrically nonlinear mechanical systems using the Spectral Submanifold (SSM) method.
The normal-form style parametrization is used to compute the coefficients of the reduced dynamics.
This allows for the analytical computation of the backbone curve.

Additionally, the code provides the computation of the SSM and backbone curve sensitivities with respect to the system parameters.
Two approaches are implemented: the direct method and the adjoint method.

The code follows the multi-index nomenclature and the methodology described in the following papers:

T. Thurnher, G. Haller, and S. Jain. "Nonautonomous spectral submanifolds for model reduction of nonlinear mechanical systems under parametric resonance". *Chaos: An Interdisciplinary Journal of Nonlinear Science* (2024). DOI: [10.1063/5.0168431](https://doi.org/10.1063/5.0168431).

M. Pozzi, J. Marconi, S. Jain, M. Li, and F. Braghin. "Topology optimization of nonlinear structural dynamics with invariant manifold-based reduced order models". *Structural and Multidisciplinary Optimization* (2025). DOI: [10.1007/s00158-025-04010-1](https://doi.org/10.1007/s00158-025-04010-1).

The [SSMTool](https://github.com/jain-shobhit/SSMTool), a package for computing SSMs for generic nonlinear systems, was used to validate the code implementation.

## Installation

To run the code, clone the [YetAnotherFEcode](https://github.com/jain-shobhit/YetAnotherFEcode) repository and run the `startup.m` script to add the necessary paths.

## Examples

The following examples are provided in the `paper_examples` folder:

- Time comparison of the SSM and backbone curve sensitivities computed with the direct and adjoint methods.

- Perturbed backbone of two coupled nonlinear oscillators.

- Optimization of a von Kármán beam.

- Optimization of a MEMS gyroscope.

## Citation
To cite this package, please cite the following:

Matteo Pozzi, Jacopo Marconi, Shobhit Jain, Mingwu Li, Francesco Braghin; Adjoint sensitivities for the optimization of nonlinear structural dynamics via spectral submanifolds. Proc. A 1 December 2025; 481 (2328): 20250244. https://doi.org/10.1098/rspa.2025.0244

## Contact
If you have any questions, create an issue or email <jacopo.marconi@polimi.it> or <matteo1.pozzi@polimi.it>.
