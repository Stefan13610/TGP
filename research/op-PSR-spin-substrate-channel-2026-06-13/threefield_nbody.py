# -*- coding: utf-8 -*-
# op-PSR-spin-substrate-channel-2026-06-13 / 3-field N-body (EXPLORATORY)
# Question (user 2026-06-13): treat it 3-field -- P1 = scalar of body 1,
# P2 = scalar of body 2, P3 = environment (cosmic background) -- and do an
# N-body calc of the superpositions P3-P1, P3-P2, P3-P1-P2.
#
# READING A (within TGP, ONE field): P1,P2,P3 are contributions to a single Phi:
#   delta-Phi_1 (from m1) + delta-Phi_2 (from m2) + Phi_bar (background P3).
#   Then "N-body P3-P1-P2" = nonlinear superposition of 2 sources on the cosmo
#   background. We compute the radiated power per decomposition and test R=P_phi/P_GR.
# READING B (genuinely 3 independent fields) handled analytically at the end.
# Run: python threefield_nbody.py

import sympy as sp
import numpy as np
line = "=" * 78
def hdr(s): print(line); print(s); print(line)

# ----------------------------------------------------------------------
hdr("PART 1 (SYMBOLIC) - R = P_phi/P_GR for an ARBITRARY quadrupole moment")
# scalar breathing power  P_phi = (G/30 c^5) <Idddot:Idddot>   (Pdot T3.3)
# tensor (GR) power        P_GR  = (G/ 5 c^5) <Idddot:Idddot>   (Pdot T3.4, same I)
# => ratio is a FIXED geometric constant, independent of WHICH masses / how many
#    bodies / how the source is decomposed. Show it for a symbolic I.
Iddd = sp.MatrixSymbol('Iddd', 3, 3)      # any 3rd-derivative reduced quadrupole
contraction = sum(Iddd[i, j]**2 for i in range(3) for j in range(3))
G, c5 = sp.symbols('G c5', positive=True)
P_phi = G/(30*c5) * contraction
P_GR  = G/(5*c5)  * contraction
R = sp.simplify(P_phi/P_GR)
print(f"P_phi/P_GR = {R}   for ANY quadrupole Iddot_ij (decomposition-independent)")
print("=> No regrouping of sources (P1, P2, P3, cross terms) can move this ratio:")
print("   it is the fixed scalar/tensor geometric factor of the SAME contraction.")

# ----------------------------------------------------------------------
hdr("PART 2 (N-BODY) - explicit P3-P1, P3-P2, P3-P1-P2 decomposition")
t, w = sp.symbols('t omega', positive=True)
m1, m2, d = sp.symbols('m1 m2 d', positive=True)
M = m1 + m2
# CM-frame positions (P3 background = constant Phi_bar; sets asymptotics only)
x1 = (m2/M)*d*sp.Matrix([sp.cos(w*t), sp.sin(w*t), 0])
x2 = -(m1/M)*d*sp.Matrix([sp.cos(w*t), sp.sin(w*t), 0])

def reduced_quad(mass, x):
    Mij = mass * (x * x.T)
    return Mij - sp.eye(3)*sp.trace(Mij)/3

def avg_contraction(I):
    I3 = sp.diff(I, t, 3)
    s = sp.simplify(sum(I3[i, j]**2 for i in range(3) for j in range(3)))
    T = 2*sp.pi/w
    return sp.simplify(sp.integrate(s, (t, 0, T))/T)

I1 = reduced_quad(m1, x1)         # P3-P1 : body 1 on background
I2 = reduced_quad(m2, x2)         # P3-P2 : body 2 on background
Itot = sp.simplify(I1 + I2)       # P3-P1-P2 : full system

A11 = avg_contraction(I1)
A22 = avg_contraction(I2)
Afull = avg_contraction(Itot)
# cross term: <Idddot_1 : Idddot_2>
I1d3, I2d3 = sp.diff(I1, t, 3), sp.diff(I2, t, 3)
cross = sp.simplify(sp.integrate(
    sum(I1d3[i, j]*I2d3[i, j] for i in range(3) for j in range(3)),
    (t, 0, 2*sp.pi/w))/(2*sp.pi/w))
print(f"<Iddd^2>  P3-P1 (body1 only) = {A11}")
print(f"<Iddd^2>  P3-P2 (body2 only) = {A22}")
print(f"cross  2*<Iddd_1:Iddd_2>     = {sp.simplify(2*cross)}")
print(f"<Iddd^2>  P3-P1-P2 (full)    = {Afull}")
print(f"identity  full == self1 + 2*cross + self2 ? "
      f"{sp.simplify(Afull - (A11 + 2*cross + A22)) == 0}")
# express full in terms of reduced mass to confirm standard quadrupole
mu = m1*m2/M
print(f"full in reduced mass: {sp.simplify(Afull - 32*mu**2*d**4*w**6)} (==0 means "
      f"<Iddd^2>=32 mu^2 d^4 w^6)")
print("\nR per decomposition (P_phi/P_GR = (1/6) for each, by Part 1):")
print("  P3-P1 : 1/6     P3-P2 : 1/6     P3-P1-P2 : 1/6")
print("=> the N-body cross term BUILDS the system quadrupole (it is already in")
print("   the falsified calc); it does NOT add a new channel or shift R.")

# ----------------------------------------------------------------------
hdr("PART 3 - what P3 (the environment) contributes radiatively")
print("P3 = cosmological background Phi_bar(t). It evolves on the Hubble time")
print("1/H0 ~ 4e17 s, vs orbital period Pb ~ 8.8e3 s: d/dt Phi_bar at orbital")
print("scale is ~ Pb*H0 ~ 2e-14 of any orbital variation -> radiatively inert.")
print("Its ONLY role is the dispersion gap omega_gap=H0 (transparent, prior calc).")
print("=> P3 carries no quadrupole; P3-P1 and P3-P2 reduce to P1 and P2 alone.")

# ----------------------------------------------------------------------
hdr("PART 4 (READING B) - genuinely THREE INDEPENDENT fields: can it rescue?")
print("If P1,P2,P3 are independent DYNAMICAL fields (a new multi-field theory):")
print(" (i) If all three are SCALAR (helicity-0): each radiates a quadrupole with")
print("     the SAME 1/6 geometric ratio. N scalar channels ADD:")
print("       R_total = (n_active_scalars)/6  >= 1/6  -> verdict only WORSENS.")
print("     Three scalars cannot interfere to GR's value without fine-tuned")
print("     destructive cancellation = Branch D (screens statics too) = forbidden.")
print(" (ii) The observed Pdot_b matches a HELICITY-2 (spin-2) quadrupole to 6e-5.")
print("     No sum of helicity-0 scalar fields can synthesize helicity-2 radiation")
print("     (Lorentz representation theory). To match GR, P3 must be a genuine")
print("     SPIN-2 tensor field carrying ~P_GR -> that is a GRAVITON.")
print(" (iii) TGP FOUNDATIONS rule 6 explicitly FORBIDS a graviton (fundamental or")
print("     composite sigma_munu=d-phi(x)d-phi). So the only multi-field structure")
print("     that could rescue is OUTSIDE TGP's axioms = a new theory, requiring its")
print("     own Newton/1PN/GWTC-3/NICER re-derivation and Phase-0 scoping.")

hdr("VERDICT (conditional)")
print("Reading A (3-field = 3 contributions to one Phi):")
print("  - R = 1/6 is a fixed geometric ratio, independent of source decomposition")
print("    (Part 1). The P1-P2 cross term is the system quadrupole, already in the")
print("    falsified calc (Part 2). P3 is radiatively inert (Part 3). NO rescue.")
print("Reading B (3 genuinely independent fields):")
print("  - 3 scalars -> R >= 1/6 (adds channels, worsens). Matching GR needs a spin-2")
print("    field = graviton = forbidden by TGP axiom. That is the structural exit,")
print("    a NEW theory, not an N-body recombination within TGP.")
print("=> Either way the orbital verdict stands. The 3-field idea, taken to its")
print("   honest conclusion, RE-DERIVES why the only escape is abandoning the")
print("   single-(scalar)-field, no-graviton core of TGP itself.")
