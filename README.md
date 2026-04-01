# Optical Coefficient Inversion

## Overview

This project focuses on the numerical inversion of optical coefficients from measured or simulated data. The objective is to recover physical parameters such as absorption and scattering coefficients by solving an inverse problem based on a forward physical model.

## Physical Model

The system is modeled using [describe your model here, e.g. radiative transfer equation or diffusion approximation].

Key assumptions:

* Medium properties: [homogeneous / heterogeneous]
* Boundary conditions: [specify]
* Source/illumination: [specify]

## Mathematical Formulation

The problem is formulated as an inverse problem.

Given measured data:

```
y = F(x)
```

Recover:

```
x
```

Where:

* F is the forward model
* x represents the optical coefficients
* y corresponds to the observed data

Example governing equation (if applicable):

```
∇ · (D ∇ϕ) - μ_a ϕ = -S
```

## Numerical Method

The inversion is performed using the Pointwise Unconstrained Minimization Approach (PUMA), which formulates the inverse problem as a pixel-wise (or pointwise) optimization problem.

In this framework, the optical coefficients are estimated by minimizing a cost function that measures the discrepancy between observed data and model predictions.

The forward model is embedded within an analysis-by-synthesis scheme:

* A trial set of optical parameters is proposed
* The forward model is solved to generate synthetic data
* The mismatch with measured data is evaluated
* Parameters are iteratively updated to reduce this mismatch

The resulting optimization problem is solved using the Spectral Gradient Method (SGM), an efficient first-order iterative method that improves convergence over standard gradient descent by incorporating spectral step-length estimation.

Key features of the method:

* Fully iterative and gradient-based
* Does not require explicit Hessian computation
* Suitable for large-scale inverse problems

Important considerations:

* The inverse problem is inherently ill-posed and may require regularization
* Convergence depends strongly on the initial guess
* Sensitivity to noise in the observed data can affect stability


## Results

Representative results include:

* Reconstructed optical coefficients
* Error analysis
* Comparison with ground truth (if available)

(Insert plots or figures here)

## Project Structure

```
src/        # core algorithms and solvers
notebooks/  # exploratory analysis and testing
results/    # generated figures and outputs
```

## Installation

Using Conda (recommended):

```
conda env create -f environment.yml
conda activate [your-environment-name]
```

## Usage

```
python main.py
```

## Future Work

* Improve robustness of the inversion algorithm
* Extend to higher-dimensional models (2D/3D)
* Incorporate experimental data
* Explore more advanced regularization techniques

## Author

Gabriel Alvarez Castrillo
