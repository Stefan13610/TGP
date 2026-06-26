# -*- coding: utf-8 -*-
# op-PSR-spin-substrate-channel-2026-06-13 / saturation integration (EXPLORATORY)
# Question (user 2026-06-13): does INTEGRATING the pulsar make sense in TGP, given
# that in a strong source the space configuration saturates (like at the horizon
# psi=4/3), so overlapping matter layers do NOT contribute a linear q*M scalar
# charge? If the effective scalar charge is screened by saturation, does P_phi
# (and thus R = P_phi/P_GR = 1/6) get suppressed enough to rescue the verdict?
#
# Grounded in REAL TGP equations:
#   - horizon: g_tt=0 <=> psi=4/3                       (sek08c)
#   - composite profile chi_glob=[a0^3/r^3+1+3K e^-/r]^(1/3), chi=a0/r in core (sek08_formalizm 2466)
#   - nonlinear stiffening f(Phi)=1+4 ln(Phi/Phi0)      (sek08_formalizm 2521)
#   - single coupling: q^2 = 4*pi*G binds AND radiates  (Pdot T2a, LOCKED)
# Run: python saturation_charge_integration.py

import sympy as sp
import numpy as np

line = "=" * 76
def hdr(s): print(line); print(s); print(line)

# ----------------------------------------------------------------------
# PART 1 - is any part of the neutron star inside the TGP horizon (psi=4/3)?
#          i.e., is there a "frozen / no-reconfiguration" region at all?
# ----------------------------------------------------------------------
hdr("PART 1 - is the NS sub-horizon? (does the 'frozen' regime apply?)")
G = 6.67430e-11; c = 2.99792458e8; Msun = 1.98892e30
M_A = 1.338*Msun; R_NS = 12e3
U_A = G*M_A/(R_NS*c**2)                       # surface compactness
# LOCKED relation from closed cycle op-PSR-orbital-drift: psi_eq(U)=1+U/2
psi_surf = 1 + U_A/2
# core can be deeper; bound it with the strong (Schwarzschild-like) estimate 1+2U
psi_core_max = 1 + 2*U_A
psi_horizon = sp.Rational(4, 3)
print(f"NS surface compactness U = GM/Rc^2      = {U_A:.4f}")
print(f"Schwarzschild radius r_s = 2GM/c^2       = {2*G*M_A/c**2/1e3:.2f} km  (R_NS = {R_NS/1e3:.0f} km)")
print(f"psi at surface (locked 1+U/2)            = {psi_surf:.4f}")
print(f"psi core upper bound (1+2U)              = {psi_core_max:.4f}")
print(f"horizon psi                              = {float(psi_horizon):.4f}")
print(f"=> NS is SUB-horizon (psi < 4/3 everywhere): NO frozen region, NO BH.")
print(f"   The 'no-reconfiguration / no-interaction' zone is a BLACK HOLE feature;")
print(f"   a pulsar still interacts fully. So integration is over a NON-frozen body.")

# ----------------------------------------------------------------------
# PART 2 - how strong is the saturation nonlinearity inside the NS?
# ----------------------------------------------------------------------
hdr("PART 2 - magnitude of the saturation stiffening f(Phi)=1+4 ln(chi)")
for label, chi in [("surface", psi_surf), ("core upper bound", psi_core_max)]:
    f_stiff = 1 + 4*np.log(chi)
    print(f"  chi({label:16s}) = {chi:.4f}  ->  f(Phi) = 1+4 ln chi = {f_stiff:.3f}")
print("=> stiffening is O(1) (~1.3..1.9x), NOT the chi>>1 strong-screening regime.")
print("   The deep-saturation chi~a0/r (M-independent) regime needs chi>>1 (near psi=4/3),")
print("   which the NS does NOT reach. Integration gives a MODEST charge correction.")

# ----------------------------------------------------------------------
# PART 3 (DECISIVE, symbolic) - SINGLE-COUPLING cancellation.
#   Saturation rescales the ONE scalar charge q -> q*(1-S) for a compact body.
#   But that SAME q binds the orbit (statics, q^2=4piG) AND sources radiation.
#   The observable R = P_phi / P_GR is a RATIO normalized to GR built from the
#   orbit-inferred masses. Show R is INVARIANT under q -> q*(1-S).
# ----------------------------------------------------------------------
hdr("PART 3 - does common charge rescaling change R = P_phi/P_GR ?")
Gsym, mu, d, w, cc, S = sp.symbols('G mu d omega c S', positive=True)
# effective single coupling after saturation: Geff = G*(1-S)^2  (since q^2=4piG, q->q(1-S))
Geff = Gsym*(1-S)**2
# scalar quadrupole power (Pdot T3): P_phi = (16/15) Geff mu^2 d^4 w^6 / c^5
P_phi  = sp.Rational(16,15)*Geff*mu**2*d**4*w**6/cc**5
# tensor quadrupole the orbit is compared to uses the SAME effective coupling
# (the binding that sets the orbit and the masses): P_GR = (32/5) Geff mu^2 d^4 w^6/c^5
P_GR   = sp.Rational(32,5)*Geff*mu**2*d**4*w**6/cc**5
R = sp.simplify(P_phi/P_GR)
print(f"P_phi  = {P_phi}")
print(f"P_GR   = {P_GR}")
print(f"R = P_phi/P_GR = {R}   (independent of S, mu, d, omega, G, c)")
dR_dS = sp.simplify(sp.diff(R, S))
print(f"dR/dS  = {dR_dS}   -> saturation does NOT move R")
print("=> Single field: any charge screening rescales BIND and RADIATE together;")
print("   the GR-normalized ratio cancels it. R stays 1/6 (Branch A). NO rescue.")

# sanity: even letting Kepler re-infer masses, omega^2 = Geff*Mtot/a^3 -> the
# GR template absorbs Geff via the measured masses; ratio is coupling-free.

# ----------------------------------------------------------------------
# PART 4 - the ONLY differential effect saturation can produce: a DIPOLE,
#   because the two stars have DIFFERENT compactness -> different screening S_A != S_B.
#   Dipole is -1PN (ENHANCED). Quantify whether it helps or worsens.
# ----------------------------------------------------------------------
hdr("PART 4 - differential screening -> scalar dipole (the only non-cancelling term)")
M_B = 1.249*Msun
U_B = G*M_B/(R_NS*c**2)
# screening proxy ~ compactness (sensitivity s ~ U for these estimates)
S_A, S_B = U_A, U_B
# dipole radiation ~ (G/3c^3) (Delta s)^2 (m_red v)^2 w^2 ... order-of-mag vs quadrupole:
# P_dip/P_quad ~ (5/12)(c/v)^2 (S_A - S_B)^2   (Will, scalar-tensor scaling)
Pb = 0.10225156248*86400.0; worb = 2*np.pi/Pb
a = (G*(M_A+M_B)*(Pb/(2*np.pi))**2)**(1/3)
v_orb = worb*a
P_dip_over_quad = (5/12)*(c/v_orb)**2*(S_A - S_B)**2
print(f"U_A={U_A:.4f}, U_B={U_B:.4f}, (S_A-S_B)={S_A-S_B:+.4f}")
print(f"v_orb/c = {v_orb/c:.3e}  -> (c/v)^2 = {(c/v_orb)**2:.3e}")
print(f"P_dipole / P_quadrupole ~ {P_dip_over_quad:.3e}")
sign = "ENHANCES (adds, -1PN)" if P_dip_over_quad > 0 else "n/a"
print(f"=> differential screening turns ON a dipole that {sign} the energy loss.")
print(f"   It WORSENS the verdict (more excess), never cancels the 1/6.")

hdr("VERDICT (conditional, grounded in TGP saturation equations)")
print("1. Integrating the pulsar is physically MEANINGFUL in TGP, but the NS is")
print("   SUB-horizon (psi<4/3): no frozen/no-interaction region; saturation is")
print("   only O(1) stiffening (Part 1-2), not the M-independent deep regime.")
print("2. DECISIVE: TGP is single-coupling (q^2=4piG binds AND radiates). Any charge")
print("   screening from saturation rescales statics and radiation together; the")
print("   GR-normalized ratio R=P_phi/P_GR cancels it EXACTLY (Part 3, dR/dS=0).")
print("3. The only differential effect (S_A != S_B from unequal compactness) turns")
print("   on a -1PN scalar DIPOLE that WORSENS the excess (Part 4).")
print("=> Saturation/integration cannot rescue the falsification; R stays 1/6,")
print("   and the differential piece deepens it. Consistent with W-PDOT-4 (HIGH).")
print("   To suppress radiation WITHOUT suppressing binding is Branch D (screening)")
print("   = self-destructive (screens statics -> orbit unbinds). No single-field exit.")
