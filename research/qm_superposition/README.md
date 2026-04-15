# Superposition from Linearity of Weak-Field Phi

## Problem

Explain quantum superposition as a consequence of TGP field
equations, rather than as a fundamental postulate.

## Connection to TGP

The soliton ODE linearizes near vacuum where g(Phi) ~ 1.
In this regime the radial equation becomes:

    f'' + (2/r) f' + f = 0

This is a LINEAR equation, so solutions superpose: any linear
combination of solutions is also a solution. This is the origin
of the superposition principle in TGP.

However, nonlinearity of g(Phi) becomes important for large
field amplitudes (dense Phi regions). This means superposition
is approximate -- it holds in the weak-field (dilute) regime
but breaks down when soliton cores overlap.

## Key Questions

1. What is the precise validity range of linearization?
   At what Phi/Phi_0 ratio does nonlinearity become significant?
2. How does the breakdown of superposition manifest physically?
   Is it related to measurement or decoherence?
3. Can we quantify corrections to superposition from nonlinear
   terms as a perturbation series?
4. Does the nonlinear regime connect to the classical limit?

## Status

The linearization argument is straightforward. Main open work is
quantifying the nonlinear corrections and connecting the breakdown
of superposition to decoherence (see qm_decoherence).

## Dependencies

- Q1: Soliton ODE and g(Phi) profile
- qm_decoherence: Nonlinear regime and classical limit
