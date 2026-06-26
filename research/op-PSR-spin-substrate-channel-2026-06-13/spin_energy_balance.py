# -*- coding: utf-8 -*-
# op-PSR-spin-substrate-channel-2026-06-13  (EXPLORATORY, STRUCTURAL_CONDITIONAL)
# Question (user 2026-06-13): in TGP momentum/rotation is a resonance with the
# substrate. Does INCLUDING the pulsar SPIN open an energy channel that changes
# the orbital energy balance R = P_phi/P_GR of J0737-3039 (the falsified cycle)?
#
# This does NOT reopen the locked op-PSR-Pdot cycle. It tests whether spin
# rescues, worsens, or leaves unchanged the verdict, using TGP's OWN coupling:
#   - scalar source q*rho, q^2 = 4*pi*G (LOCKED in Pdot T2a)
#   - radiation only from ACCELERATION -> Phi fluctuations (FOUNDATIONS 6.2)
#   - motion-resonance dPhi from moving source (STRUCTURAL_CONDITIONAL)
# All radiated-power numbers conditional on the resonance being a true radiative
# (not near-field) channel. Run: python spin_energy_balance.py

import sympy as sp
import numpy as np

line = "=" * 74
def hdr(s): print(line); print(s); print(line)

# ----------------------------------------------------------------------
# PART 1 (SYMBOLIC) — does STEADY AXISYMMETRIC spin source a time-varying
# mass quadrupole?  Radiation (scalar OR tensor) needs M3_ij = d^3/dt^3 != 0.
# ----------------------------------------------------------------------
hdr("PART 1 - symmetry: lab-frame quadrupole of a rigidly spinning body")
t, Om = sp.symbols('t Omega', real=True, positive=True)
Ix, Iy, Iz = sp.symbols('I_x I_y I_z', positive=True)   # principal moments
# Rotation about z at angular velocity Omega
R = sp.Matrix([[sp.cos(Om*t), -sp.sin(Om*t), 0],
               [sp.sin(Om*t),  sp.cos(Om*t), 0],
               [0, 0, 1]])
# (a) AXISYMMETRIC body about spin axis: Ix == Iy
I_body_axi = sp.diag(Ix, Ix, Iz)
I_lab_axi  = sp.simplify(R * I_body_axi * R.T)
M3_axi = sp.simplify(sp.diff(I_lab_axi, t, 3))
print("axisymmetric (I_x=I_y): lab inertia tensor time-dependent? ",
      sp.simplify(sp.diff(I_lab_axi, t)) != sp.zeros(3, 3))
print("axisymmetric: third derivative M3_ij =\n", M3_axi)
axi_radiates = (M3_axi != sp.zeros(3, 3))
print("=> AXISYMMETRIC SPIN RADIATES (scalar or tensor)?  ", bool(axi_radiates))

# (b) NON-axisymmetric (ellipticity eps): Ix = Iz(1+eps), Iy = Iz(1-eps) approx
eps = sp.symbols('epsilon', positive=True)
I_body_def = sp.diag(Iz*(1+eps), Iz*(1-eps), Iz)
I_lab_def = sp.simplify(R * I_body_def * R.T)
M3_def = sp.simplify(sp.diff(I_lab_def, t, 3))
# dominant time dependence: components oscillate at 2*Omega
print("\nnon-axisymmetric: M3_xx =", sp.simplify(M3_def[0, 0]))
print("=> radiates at frequency 2*Omega (the 2*Omega_spin CW channel)")

# ----------------------------------------------------------------------
# PART 2 (NUMERIC) — J0737-3039 numbers
# ----------------------------------------------------------------------
hdr("PART 2 - J0737-3039 energy reservoirs (SI)")
G   = 6.67430e-11
c   = 2.99792458e8
Msun = 1.98892e30
# masses (Kramer et al. 2021)
mA = 1.338 * Msun
mB = 1.249 * Msun
M  = mA + mB
mu = mA * mB / M
# orbit
Pb   = 0.10225156248 * 86400.0          # s
worb = 2*np.pi / Pb                      # rad/s
e    = 0.0877775
a    = (G*M*(Pb/(2*np.pi))**2)**(1/3)    # relative semi-major axis
# eccentricity enhancement (Peters 1964)
fe = (1 + 73/24*e**2 + 37/96*e**4) / (1 - e**2)**3.5
# GR orbital GW luminosity (tensor quadrupole, circular form * f(e))
P_GR = (32/5) * (G/c**5) * mu**2 * a**4 * worb**6 * fe
P_phi_orb = P_GR / 6.0                   # TGP scalar excess (LOCKED 1/6 result)
print(f"a (relative)          = {a:.4e} m")
print(f"omega_orb             = {worb:.4e} rad/s   (GW at 2*omega = {2*worb:.3e})")
print(f"P_GR (orbital tensor) = {P_GR:.4e} W")
print(f"P_phi orbital = P_GR/6 = {P_phi_orb:.4e} W   <-- the excess that falsified TGP")

# spin sector
I_NS = 1.3e38                            # kg m^2 (1.3e45 g cm^2, standard NS)
# pulsar A
P_A, Pdot_A = 0.022699378, 1.7596e-18
Om_A  = 2*np.pi / P_A
Omd_A = -2*np.pi * Pdot_A / P_A**2
Edot_spin_A = abs(I_NS * Om_A * Omd_A)   # spin-down luminosity
# pulsar B
P_B = 2.7734607
Om_B = 2*np.pi / P_B
print(f"\nOmega_spin,A          = {Om_A:.3e} rad/s  (nu_A = {1/P_A:.2f} Hz)")
print(f"Omega_spin,B          = {Om_B:.3e} rad/s")
print(f"E_rot,A = 0.5 I Om^2  = {0.5*I_NS*Om_A**2:.3e} J")
print(f"spin-down luminosity A= {Edot_spin_A:.3e} W   (-> magnetosphere/EM, NOT orbit)")

# ----------------------------------------------------------------------
# PART 3 — could spin energy enter the ORBITAL balance and mask the 1/6?
# ----------------------------------------------------------------------
hdr("PART 3 - can spin rescue the orbital verdict?")
# 3a. axisymmetric spin substrate radiation = 0 (Part 1). Conservative case.
print("3a. AXISYMMETRIC spin: M3_ij=0 (Part 1) => zero substrate radiation,")
print("    zero contribution to dE/dt. Spin is radiatively silent in TGP too.")

# 3b. spin-orbit (1.5PN) contribution to orbital Pdot, order of magnitude.
#     spin-orbit is CONSERVATIVE at leading order; its dissipative imprint on
#     Pdot is suppressed by (v/c)^3 * (S/(mu*a*v)) relative to quadrupole.
v_orb = worb * a
S_A = I_NS * Om_A                        # spin angular momentum of A
L_orb = mu * np.sqrt(G*M*a*(1-e**2))     # orbital angular momentum
so_suppression = (v_orb/c)**3 * (S_A / L_orb)
print(f"\n3b. v_orb/c            = {v_orb/c:.3e}")
print(f"    S_A / L_orb        = {S_A/L_orb:.3e}")
print(f"    spin-orbit imprint on Pdot ~ {so_suppression:.3e} x quadrupole")
print(f"    => fractional shift in R far below the 1/6 = 0.167 gap; cannot rescue.")

# 3c. energy-transfer ceiling: even if ALL spin-down power were (impossibly)
#     funneled into the orbit, how does it compare to the needed offset?
needed_offset = P_phi_orb                # must cancel ~ +1/6 P_GR
print(f"\n3c. offset needed to mask TGP excess = {needed_offset:.3e} W")
print(f"    spin-down reservoir A            = {Edot_spin_A:.3e} W (ratio "
      f"{Edot_spin_A/needed_offset:.1f}x)")
print("    BUT no TGP mechanism couples spin-down (EM/magnetospheric) into the")
print("    companion orbit; spin-orbit is conservative (3b). Reservoir is")
print("    large but DISCONNECTED. No rescue path in LOCKED+conditional inventory.")

# ----------------------------------------------------------------------
# PART 4 — the genuine NEW TGP prediction the spin question surfaces:
#          a SCALAR component of spin-down CW emission at 2*Omega_spin.
# ----------------------------------------------------------------------
hdr("PART 4 - new TGP-native channel: scalar spin-down at 2*Omega (eps != 0)")
# tensor CW power for ellipticity eps:  P = (32/5)(G/c^5) I^2 eps^2 Om^6
# TGP adds scalar channel = (1/6) of that (same I3 structure, P_phi=P_GR/6)
def P_tensor_cw(I, eps, Om): return (32/5)*(G/c**5)*I**2*eps**2*Om**6
P_t_unit = P_tensor_cw(I_NS, 1.0, Om_A)          # per eps^2
# eps that would saturate the observed spin-down budget (tensor+scalar):
P_total_unit = P_t_unit * (1 + 1/6)
eps_sat = np.sqrt(Edot_spin_A / P_total_unit)
print(f"tensor CW power per eps^2 (A)      = {P_t_unit:.3e} W")
print(f"epsilon saturating spin-down budget= {eps_sat:.2e}")
print(f"  (cf. LVK NS ellipticity limits ~1e-6..1e-8 -> physically allowed band)")
print(f"emission frequency                 = 2*nu_spin = {2/P_A:.1f} Hz")
print(f"  -> a CONTINUOUS-WAVE band, 6 orders above the orbital 2*f_orb="
      f"{2*worb/(2*np.pi)*1e3:.3f} mHz")
print("  -> separate observable (CW / spin-down budget); does NOT enter Pdot_b.")

hdr("VERDICT (conditional)")
print("- Axisymmetric steady spin: radiates ZERO into substrate (Part 1 symmetry).")
print("- Spin-orbit imprint on orbital Pdot: ~%.0e x quadrupole (Part 3b) - negligible." % so_suppression)
print("- Spin-down reservoir is large but DISCONNECTED from the orbit (Part 3c).")
print("- Spin DOES give a new TGP scalar CW channel at 2*Omega_spin (Part 4),")
print("  but it is a DIFFERENT observable and does not touch R=P_phi/P_GR.")
print("=> Including spin does NOT rescue the orbital verdict; R stays 1/6 (Branch A).")
print("   The picture is unchanged, exactly as the user permitted it might be.")
