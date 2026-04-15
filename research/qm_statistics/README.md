# Fermi-Dirac vs Bose-Einstein from Topology

## Problem

Derive quantum statistics (fermions vs bosons) from the
topological properties of TGP solitons, without postulating
the symmetrization postulate.

## Connection to TGP

Solitons in TGP carry a topological charge given by their
winding number n. When two identical solitons are exchanged
(swapped in position), the total wave function picks up a
phase factor (-1)^n.

This gives:
- Odd winding number (n = 1, 3, ...): exchange phase = -1.
  These are fermions, obeying Fermi-Dirac statistics and
  the Pauli exclusion principle.
- Even winding number (n = 0, 2, ...): exchange phase = +1.
  These are bosons, obeying Bose-Einstein statistics.

The spin-statistics connection follows naturally: odd winding
gives both spin-1/2 (from qm_spin) and Fermi statistics.
This is the spin-statistics theorem derived from topology
rather than from relativistic QFT axioms.

## Key Questions

1. Is the exchange phase exactly (-1)^n, or are there
   corrections from soliton interactions during exchange?
2. Does the exclusion principle follow rigorously -- do two
   n=1 solitons at the same location have zero amplitude?
3. Can we recover parastatistics or anyonic behavior for
   restricted geometries (2D substrate)?
4. How robust is the topological argument against
   perturbations of the substrate?

## Status

The topological argument parallels the Skyrme model result
for baryons. Formal derivation of the exchange phase from
the path integral over soliton configurations is in progress.

## Dependencies

- qm_spin: Winding number and rotation properties
- Q1: Soliton stability under exchange paths
