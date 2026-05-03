---
title: "Born Rule from Soliton Tail Interference"
date: 2026-05-03
tgp_status:
  folder_status: "needs-bridge"
  level: L1
  kind: derivation
  core_compatibility: "unknown"
  last_reviewed_against_core: "unknown"
  may_edit_core: false
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges: []
  depends_on: []
  impacts: []
  source_of_status:
    - "README.md H1: 'Born Rule from Soliton Tail Interference'"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: "2026-05-03"
---

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
