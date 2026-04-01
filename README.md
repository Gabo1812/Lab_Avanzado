# Optical Coefficient Inversion in Thin Films

## Overview

This project focuses on the numerical inversion of optical properties of thin films from transmittance data. The goal is to recover the refractive index \( n(\lambda) \), extinction coefficient \( k(\lambda) \), and film thickness \( d \) by solving a nonlinear inverse problem.

The implementation is based on the **PUMA (Pointwise Unconstrained Minimization Approach)** framework, applied to experimental and simulated UV-Vis transmittance spectra of layered systems such as:

- ITO (Indium Tin Oxide)
- WO₃ (Tungsten Oxide)
- Soda–lime glass (SLG) substrates

---

## Physical Model

The forward model describes the optical transmittance of a thin absorbing film deposited on a transparent substrate, accounting for:

- Multiple internal reflections (interference effects)
- Absorption within the film
- Wavelength-dependent refractive indices

The transmittance is modeled as:

\[
T(\lambda) = \frac{A x}{B - Cx + Dx^2}
\]

with:

- \( x = e^{-\alpha d} \), \( \alpha = \frac{4\pi k}{\lambda} \)
- \( \phi = \frac{4\pi n d}{\lambda} \)

This model is **nonlinear, highly coupled, and oscillatory**, making inversion nontrivial.

### Key assumptions

- Planar, homogeneous thin films
- Known substrate refractive index \( s(\lambda) \)
- Normal incidence
- Coherent multiple reflections
- No scattering (pure absorption via \( k(\lambda) \))

---

## Inverse Problem Formulation

Given measured transmittance data:

\[
T_{\text{meas}}(\lambda_i)
\]

we solve:

\[
\min_{d,\,n(\lambda),\,k(\lambda)} \sum_i \left[T_{\text{model}}(\lambda_i) - T_{\text{meas}}(\lambda_i)\right]^2
\]

### Ill-posedness

This problem is fundamentally underdetermined:

- For each \( \lambda \): 1 equation, 2 unknowns \( (n, k) \)
- Infinite solution manifold without additional constraints

To regularize the problem, **physical constraints are imposed**:

- \( n(\lambda) \geq 1 \), \( k(\lambda) \geq 0 \)
- Monotonic behavior: \( n'(\lambda) \leq 0 \), \( k'(\lambda) \leq 0 \)
- Convexity conditions on \( n(\lambda) \) and \( k(\lambda) \)
- Inflection point in \( k(\lambda) \)

These constraints encode **physical dispersion relations** near the absorption edge.

---

## Numerical Method

### PUMA Framework

Instead of directly solving a constrained optimization problem, PUMA reformulates it as an **unconstrained problem** via a change of variables:

- Positivity enforced via squares (e.g. \( n = 1 + u^2 \))
- Convexity enforced through second derivatives

This transforms the problem into a high-dimensional nonlinear minimization:

\[
\min f(x)
\]

where \( x \) includes:

- Film thickness \( d \)
- Inflection point \( \lambda_{\text{infl}} \)
- Discretized representations of \( n(\lambda) \), \( k(\lambda) \)

---

### Optimization Algorithm

The minimization is performed using the **Spectral Gradient Method (SGM)**:

- First-order method (gradient-based)
- No Hessian required
- Adaptive step size (Barzilai–Borwein type)
- Suitable for large-scale problems (\( \mathcal{O}(10^4) \) variables)

Key properties:

- Fast per iteration
- Sensitive to scaling and initialization
- Can converge to local minima

---

## Practical Implementation Notes

From actual simulations:

- Each run can take **3–4 hours** (cluster execution)
- Strong sensitivity to:
  - Initial parameter ranges
  - Inflection point \( \lambda_{\text{infl}} \)
- Noise in experimental data significantly affects stability

### Observed issues

- Non-physical solutions (e.g. \( k < 0 \)) when constraints are poorly enforced
- Slow convergence due to ill-conditioning
- Multiple local minima

### Mitigation strategies

- Savitzky–Golay filtering of experimental data
- Careful tuning of parameter ranges
- Using PUMA results as **initial seed** for more robust forward models (e.g. Fortran implementations)

---

## Results

Preliminary results show:

- Good agreement between simulated and experimental transmittance
- RMSE as low as ~0.06% in favorable cases
- Physically consistent \( n(\lambda) \), \( k(\lambda) \) after parameter tuning

However:

- Full convergence is not always achieved
- Trade-off between fit quality and physical plausibility

*(Insert plots: transmittance fit, \( n(\lambda) \), \( k(\lambda \))*

---

## Project Structure
- src/ # core numerical routines and PUMA interface
- notebooks/ # data analysis and visualization
- data/ # experimental transmittance spectra
- results/ # output coefficients and plots


---

## Installation

```bash
conda env create -f environment.yml
conda activate optical-inversion
```

## Future Work
- Improve regularization (e.g. Tikhonov, Bayesian approaches)
- Extend to multilayer systems (ITO–WO₃ stacks)
- Incorporate reflectance data (better conditioning)
- Replace PUMA with hybrid methods (global + local optimization)
- Parallelize parameter sweeps

## Author

Gabriel Alvarez Castrillo
Physics Student, Universidad de Costa Rica