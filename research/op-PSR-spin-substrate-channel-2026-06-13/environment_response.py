# -*- coding: utf-8 -*-
# op-PSR-spin-substrate-channel-2026-06-13 / environment response (EXPLORATORY)
# Question (user 2026-06-13): the previous calcs treated only the SOURCE's own
# system (bare 2-body emission). But in TGP the mass forces a space configuration
# which forces the next, etc. Did we include the RESONANT RESPONSE OF THE
# SURROUNDING SUBSTRATE (the medium), or only the isolated binary?
#
# Here we treat the substrate as a MEDIUM with its own response and ask whether
# it (a) is transparent/dispersive, (b) resonates at orbital frequency, (c) emits
# a motion wake (Cherenkov), (d) can pump energy BACK into the orbit (the only way
# the environment could rescue Pdot).
#
# TGP medium: delta-Phi obeys  (1/c^2) d_tt - Lap + gamma) dPhi = source
#   => dispersion  omega^2 = c^2 k^2 + c^2 gamma,  gap omega_gap = c*sqrt(gamma)
#   with m_sp ~ H0/c (cosmological), gamma = (H0/c)^2  (Appendix E / G-substrate)
# Run: python environment_response.py

import numpy as np
line = "=" * 76
def hdr(s): print(line); print(s); print(line)

c   = 2.99792458e8
H0  = 70.0*1000/(3.0857e22)            # s^-1  (70 km/s/Mpc)
gamma = (H0/c)**2                      # m^-2   (substrate mass term m_sp^2)
w_gap = c*np.sqrt(gamma)               # = H0 : the ONLY intrinsic substrate freq
Pb  = 0.10225156248*86400.0
w_orb = 2*np.pi/Pb
w_gw  = 2*w_orb                        # quadrupole emission frequency

# orbit kinematics (J0737)
G=6.6743e-11; Msun=1.98892e30; M=(1.338+1.249)*Msun
a = (G*M*(Pb/(2*np.pi))**2)**(1/3); v_orb = w_orb*a

# ----------------------------------------------------------------------
hdr("PART A - is the substrate medium TRANSPARENT at orbital frequency?")
print(f"substrate mass gap  omega_gap = c*sqrt(gamma) = H0 = {w_gap:.3e} rad/s")
print(f"emission frequency  omega_GW  = 2*omega_orb     = {w_gw:.3e} rad/s")
ratio = w_gw/w_gap
print(f"omega_GW / omega_gap                            = {ratio:.3e}")
vg_over_c = np.sqrt(1 - (w_gap/w_gw)**2)
med_corr = 1 - vg_over_c                # fractional deviation of transport from c
print(f"group velocity v_g/c = sqrt(1-(gap/w)^2)        = {vg_over_c:.15f}")
print(f"fractional medium correction to transport       = {med_corr:.3e}")
print("=> emission is ~14 orders ABOVE the gap: medium is transparent, v_g=c,")
print("   no dispersive trapping, no absorption band. Flux at infinity = emitted.")

# ----------------------------------------------------------------------
hdr("PART B - does the orbit's motion drive a WAKE (Cherenkov/superluminal)?")
# substrate signal speed c_s = c0/sqrt(psi) ~ c0 in the weak orbital region
c_s = c                                # psi ~ 1 between the stars
mach = v_orb/c_s
print(f"orbital velocity v_orb       = {v_orb:.3e} m/s")
print(f"substrate signal speed c_s   ~ {c_s:.3e} m/s  (psi~1 at orbital separation)")
print(f"Mach number v_orb/c_s        = {mach:.3e}")
print("=> strongly SUB-luminal (Mach<<1): NO Cherenkov wake, NO motion-resonance")
print("   emission beyond the quadrupole. The 'wake' channel (v>c_s) is closed.")

# ----------------------------------------------------------------------
hdr("PART C - can the environment RESONATE and pump energy BACK to the orbit?")
# A medium can only ANTI-damp (return energy to source) if it has a mode at the
# driving frequency. The substrate's only intrinsic mode is omega_gap = H0.
detuning = w_orb/w_gap
print(f"substrate intrinsic mode (only one): omega_gap = H0 = {w_gap:.3e} rad/s")
print(f"drive (orbital) frequency:           omega_orb     = {w_orb:.3e} rad/s")
print(f"detuning omega_orb/omega_gap                       = {detuning:.3e}")
print("=> driven ~14 orders off the only resonance: zero resonant back-transfer.")
print("   The environment cannot anti-damp the orbit; no rescue channel exists.")

# ----------------------------------------------------------------------
hdr("PART D - the DECISIVE bookkeeping: what does Pdot_b actually measure?")
print("Pdot_b measures energy LEAVING THE ORBIT (near-zone radiation reaction),")
print("not energy arriving at infinity. So what the environment does downstream")
print("(absorb, scatter, re-emit) is IRRELEVANT to Pdot_b:")
print("  - if the substrate ABSORBS the scalar wave -> orbit still lost it -> R unchanged")
print("  - only ENVIRONMENT->ORBIT back-pumping could shrink R, and Part C kills that")
print("Near-zone radiation reaction uses the RETARDED Green fn of the medium;")
print("for omega>>gap it equals the massless one to O((gap/omega)^2) ~ %.0e." % ((w_gap/w_gw)**2))
print("=> environmental response changes the near-zone balance by ~1e-29. Negligible.")

# ----------------------------------------------------------------------
hdr("PART E - background DRESSING: does propagation on configured Phi-bar")
hdr("         change the RATIO R = P_phi/P_GR ?")
print("The cascade 'mass->config->next config' is the STATIC self-consistent")
print("profile chi_glob(r); a stationary configuration radiates nothing. Its only")
print("radiative role is to be the curved BACKGROUND the waves propagate on")
print("(redshift, backscatter tails). But in TGP BOTH the scalar signal AND the")
print("'tensor' GW are perturbations of the SAME field Phi on the SAME background")
print("=> any dressing multiplies numerator and denominator equally -> cancels in R.")
print("   (Same single-coupling cancellation as the saturation calc: dR=0.)")

hdr("VERDICT (conditional)")
print("We had computed the bare 2-body emission. Adding the ENVIRONMENT's response:")
print(f" A. medium transparent at omega_GW (14 orders above gap; correction ~{med_corr:.0e})")
print(f" B. sub-luminal orbit (Mach {mach:.0e}) -> no Cherenkov/wake channel")
print(f" C. no substrate mode near omega_orb (detuning {detuning:.0e}) -> no back-pumping")
print(" D. Pdot_b = energy leaving the ORBIT; downstream absorption cannot rescue it")
print(" E. background dressing acts on scalar & tensor alike -> cancels in R")
print("=> Including the resonant response of the surroundings does NOT change the")
print("   verdict. R stays 1/6. The environment is transparent, non-resonant, and")
print("   (decisively) downstream of the orbital energy balance that Pdot_b probes.")
