# Spin 1/2 from Soliton Topology

## Problem

Derive half-integer spin from the topology of 3D solitons in TGP,
without postulating spinor fields.

## Connection to TGP

A 3D soliton is a map from S^3 (compactified space) to R (field
values). The relevant homotopy group is pi_3(S^3) = Z, which
classifies solitons by an integer winding number n.

Under a 2pi spatial rotation, a soliton with winding number n
picks up a sign factor (-1)^n. This means:

- Odd winding number n: rotation by 2pi gives a sign change (-1).
  This is exactly spin-1/2 behavior.
- Even winding number n: rotation by 2pi gives no sign change (+1).
  This is integer spin behavior.

The connection to the SU(2) double cover of SO(3) is natural:
the soliton's topological structure automatically implements
the double cover. A 4pi rotation is needed to return an odd-n
soliton to its original state.

## Key Questions

1. Does the winding number n map exactly to 2*spin?
   (n=1 -> spin 1/2, n=2 -> spin 1, etc.)
2. How does the spin-statistics connection emerge?
   (See qm_statistics for the statistics side.)
3. Can we recover the full SU(2) representation theory
   (Clebsch-Gordan coefficients, addition of angular momenta)?
4. What determines which winding numbers are stable?

## Status

The topological argument is clean and compelling. Main open
work is connecting it to the full angular momentum algebra
and verifying stability of different winding sectors.

## Dependencies

- Q1: Soliton existence in 3D
- qm_statistics: Spin-statistics theorem
