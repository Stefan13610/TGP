"""
TGP n-body interaction package
===============================
Core modules for investigating n-body physics in the Theory of Generated Space.

PHYSICAL GLOSSARY
-----------------
C (source strength):
    Dimensionless scalar characterizing an EFFECTIVE point source in the
    N-body package.
    Enters the screened Yukawa profile:
        delta(r) = C * exp(-m_sp*r) / r.
    Canonical interpretation in `nbody`: C is the post-projection coupling
    C_eff obtained from the classical defect sector.

m_i (inertial mass):
    TGP AXIOM: inertial mass equals source strength, m_i = C_i.
    This is the equivalence principle analog in TGP — the source of the
    scalar field IS the inertia. All dynamics modules enforce this convention.

beta, gamma (self-interaction couplings):
    Appear in the TGP energy functional:
        E[Phi] = ... + (beta/3*Phi_0)*Phi^3 - (gamma/4*Phi_0^2)*Phi^4
    Vacuum condition: beta = gamma (ensures Phi = Phi_0 is a stable vacuum).

m_sp (screening mass):
    Effective Yukawa screening mass used by the N-body layer:
        m_sp = sqrt(3*gamma - 2*beta)
    At vacuum condition (beta=gamma): m_sp = sqrt(gamma).
    This is not the raw classical defect tail scale: the unprojected defect
    equation has an oscillatory tail ~ sin(r)/r.

softening (epsilon):
    NUMERICAL REGULATOR, not a physical cutoff. Replaces r → sqrt(r^2 + eps^2)
    to prevent 1/r divergences in numerical integration. Results should be
    checked for convergence as epsilon → 0. See dynamics_v2.py docstring.

THEORY-LAYER LABELS
-------------------
To avoid mixing the classical and effective descriptions, it is useful to
distinguish:
  CLASSICAL:    soliton/defect ODE with oscillatory tail ~ sin(r)/r
  EFT-DERIVED:  bridge from defect to effective source (C_eff, m_sp)
  N-BODY:       equations of motion built from Yukawa-source V_2 and V_3

STATUS LEVELS
-------------
Throughout this package, results are classified by rigor:
  EXACT:        Closed-form analytical result from the theory.
  APPROXIMATE:  Analytical approximation (e.g., saddle-point, perturbative).
  NUMERICAL:    Computed by numerical methods (integration, root-finding).
  HEURISTIC:    Threshold-based or phenomenological (legacy, being replaced).

MODULES
-------
tgp_field              -- Field profiles, gradients, energy density
bridge_nbody           -- Canonical bridge: defect -> C_eff -> Yukawa source -> EOM inputs
pairwise               -- 2-body effective potential and forces (EXACT)
three_body_terms       -- Irreducible 3-body interactions (APPROXIMATE/NUMERICAL)
three_body_force_exact -- Exact 3-body forces via Feynman 2D integral (EXACT)
equilibria             -- Equilibrium finder: pairwise (EXACT) + corrected (APPROX)
stability              -- Hessian eigenvalue analysis with zero-mode projection
dynamics_v2            -- Forces + integrators (forces EXACT/APPROX, dynamics NUMERICAL)
eom_tgp                -- Jawna postać EOM: V_2 + V_3(I); siły z ∂I/∂d (warstwa analityczna)
dynamics_backends      -- integration pair; COM / P / L / Hamiltonian helpers; -grad V (FD); TGP_INTEGRATION_BACKENDS
lyapunov               -- Benettin (RK4 / leapfrog + styczna), spektrum QR; P1 ex148–ex194 + verify_nbody_lyapunov_quick; TeX: tgp_lyapunov_benettin.tex
configurations         -- Standard n-body geometries (equilateral, collinear, n-gon)
nbody_energy           -- Total n-body energy assembler

ROADMAP: PROVEN vs CONJECTURED
------------------------------
  PROVEN (exact pairwise sector):
    - Equilateral/collinear equilibria exist for beta > 9*C/2
    - Equilibrium polynomial: d^2 - 4*beta*d + 18*gamma*C = 0
    - TGP violates Earnshaw's theorem (positive-definite pairwise Hessian
      at d_well, verified by rigorous zero-mode projection)

  NUMERICAL EVIDENCE (strong):
    - Bounded oscillations at d_well under perturbation (leapfrog + RK853)
    - |V_3/V_2| ~ O(C) at equilibrium (3-body is a small correction for C<<1)
    - Perturbative correction d* = d0 + delta_d converges for C < 1

  CONJECTURED (requires further investigation):
    - Stability persists for full V_2 + V_3 Hessian (not just pairwise)
    - Accumulation of 3-body corrections in N >> 1 systems (galaxy rotation)
    - TGP pairwise makes rotation curves steeper (opposite of dark matter)
"""

from . import tgp_field
from . import bridge_nbody
from . import pairwise
from . import three_body_terms
from . import three_body_force_exact
from . import nbody_energy
