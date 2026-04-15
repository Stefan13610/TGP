# Entanglement from Substrate Graph Correlations

## Problem

Derive quantum entanglement from the structure of the TGP
substrate, without postulating nonlocal state vectors.

## Connection to TGP

The substrate Gamma = (V, E) is a graph with vertices V and
edges E. The scalar field Phi is defined on this graph. Two
solitons that share substrate edges have correlated Phi
fluctuations -- the field value at shared vertices couples
their dynamics.

This provides a mechanism for entanglement: correlations are
mediated by the substrate graph topology, not by signals
traveling through space. Two solitons can be "entangled" if
they were formed from the same substrate region and still
share graph edges even after spatial separation.

## Key Questions

1. Can substrate-mediated correlations violate Bell inequalities?
   This is the critical test for genuine entanglement.
2. How does the graph structure preserve correlations over
   large spatial separations?
3. What is the entanglement entropy in terms of shared edges?
4. How does decoherence (loss of shared edges) destroy
   entanglement?
5. Does this mechanism reproduce the correct tensor product
   structure of quantum Hilbert spaces?

## Status

Most speculative folder in the QM research program. The basic
idea (shared substrate edges = correlations) is stated, but no
quantitative results yet. Bell inequality analysis is the
highest priority next step.

## Dependencies

- Substrate graph model (continuum_limit)
- qm_decoherence: For understanding loss of entanglement
