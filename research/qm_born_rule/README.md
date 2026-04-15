# Born Rule from Soliton Tail Interference

## Problem

Derive the Born rule P = |psi|^2 from soliton dynamics in TGP,
without postulating it as a separate axiom.

## Connection to TGP

Solitons in TGP have exponentially decaying tails beyond their core.
When two soliton tails overlap in the substrate, the interaction
energy depends on |A_tail|^2, where A_tail is the tail amplitude.
This quadratic dependence on amplitude is exactly the Born rule:
probability of detection (interaction) is proportional to the
square of the wave amplitude.

The key mechanism is:
- Soliton tail amplitude A_tail plays the role of psi.
- Overlap integral of tails gives interaction strength ~ |A_tail|^2.
- Detection probability = interaction probability ~ |psi|^2.

## Key Questions

1. Can we derive the exact |psi|^2 form from the soliton ODE,
   or only the quadratic scaling?
2. How does normalization emerge (probabilities summing to 1)?
3. Does the interference pattern (cross terms) come out correctly
   for superpositions of soliton tails?
4. What sets the proportionality constant?

## Status

Depends on Q1 (soliton existence) results. Preliminary analysis
suggests the quadratic dependence is robust, but the full
derivation requires confirmed soliton profiles from numerics.

## Dependencies

- Q1: Soliton existence and tail profile
- qm_superposition: Linearity of weak-field regime
