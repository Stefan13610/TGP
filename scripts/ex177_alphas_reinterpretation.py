#!/usr/bin/env python3
"""
ex177_alphas_reinterpretation.py
Sesja v45, 2026-04-05

Reinterpretacja alpha_s z ODE substratowym (alpha=1).

PROBLEM:
  Formula alpha_s = N_c^2 * g0* / (4*Phi_0) uzywala g0* z warunku
  B_tail(g0*) = 0 (H1). Z alpha=1 H1 NIE ISTNIEJE (ex176).

  Pytanie: Czy formula ZMIENIA SIE z alpha, czy tylko g0* wymaga
  nowej definicji?

PODEJSCIA:
  1. Reverse-engineering: jaki g0* daje alpha_s = 0.1179?
  2. Czy running (ex144: mu_TGP ~ 120 GeV) rozwiazuje problem?
  3. Nowe definicje g0*: phi-FP, self-consistency, topologiczne
  4. Czy formuła alpha_s ma postac phi-bazowa (bez g0*)?

LANCUCH DERYWACJI (z dodatekV):
  g_s^2 = 2 * J_c * e^2 * a_sub^2 / hbar_0^2
  J_c = N_c * g0* / Phi_0
  alpha_s = g_s^2/(4*pi) = N_c * g0* * kappa
  kappa = N_c / (4*Phi_0)   [soliton density parameter]
  => alpha_s = N_c^2 * g0* / (4*Phi_0)
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

PHI = (1 + np.sqrt(5)) / 2
PHI0 = 25.0
PHI0_B = 24.783  # Brannen
N_c = 3
ALPHA_S_PDG = 0.1179

print("=" * 72)
print("ex177: Reinterpretacja alpha_s z ODE alpha=1")
print("=" * 72)

# ---- ODE solver ----
def solve_ode(g0, r_max=150, n_points=40000, alpha=1):
    def rhs(r, y):
        g, gp = y
        if g < 1e-12: g = 1e-12
        if alpha == 1:
            src = 1 - g
        else:
            src = g**2 * (1 - g)
        gpp = src - (alpha/g)*gp**2 - (2.0/r)*gp if r > 1e-12 else src/3.0
        return [gp, gpp]
    r_eval = np.linspace(1e-6, r_max, n_points)
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0], method='RK45',
                    t_eval=r_eval, rtol=1e-10, atol=1e-12, max_step=0.04)
    return sol.t, sol.y[0]

def extract_AB(g0, alpha=1, window=(60, 120)):
    r, g = solve_ode(g0, r_max=150, alpha=alpha)
    delta = g - 1.0
    mask = (r > window[0]) & (r < window[1])
    r_w, d_w = r[mask], delta[mask]
    if len(r_w) < 300:
        return None, None
    y = d_w * r_w
    M = np.column_stack([np.cos(r_w), np.sin(r_w)])
    result = np.linalg.lstsq(M, y, rcond=None)
    return result[0][0], result[0][1]

def get_Atail(g0, alpha=1):
    A, B = extract_AB(g0, alpha=alpha)
    if A is None: return None
    return np.sqrt(A**2 + B**2)

# ===== 1. REVERSE ENGINEERING =====
print("\n--- 1. Reverse engineering: jaki g0* daje alpha_s = 0.1179? ---\n")

g0_needed_25 = 4 * PHI0 * ALPHA_S_PDG / N_c**2
g0_needed_B = 4 * PHI0_B * ALPHA_S_PDG / N_c**2

print(f"  Formula: alpha_s = N_c^2 * g0* / (4 * Phi_0)")
print(f"  => g0* = 4 * Phi_0 * alpha_s / N_c^2")
print()
print(f"  Phi_0 = 25.000: g0* potrzebne = {g0_needed_25:.5f}")
print(f"  Phi_0 = 24.783: g0* potrzebne = {g0_needed_B:.5f}")
print(f"  Stare g0* (H1, alpha=2) =    1.24915")
print()
print(f"  Roznica wzgledem starego g0*:")
print(f"    Phi_0=25: {(g0_needed_25 - 1.24915)/1.24915*100:+.2f}%")
print(f"    Phi_0=B:  {(g0_needed_B - 1.24915)/1.24915*100:+.2f}%")

# ===== 2. RUNNING ARGUMENT (ex144) =====
print("\n--- 2. Running argument: alpha_s^bare -> alpha_s(M_Z) ---\n")

# If alpha_s^TGP is the BARE value at soliton scale mu_TGP,
# then 1-loop running gives alpha_s(M_Z).
# ex144: alpha_s^bare = 0.1134, mu_TGP = 120 GeV -> alpha_s(M_Z) = 0.1179

M_Z = 91.1876
N_f = 5
b0 = (11*N_c - 2*N_f) / 3

# For each Phi_0, what g0* is needed to give bare alpha_s that RUNS to 0.1179?
# alpha_s(M_Z) = alpha_bare / [1 + b0*alpha_bare/(2pi)*ln(M_Z^2/mu^2)]
# We need: alpha_bare such that after running FROM mu_TGP TO M_Z, we get 0.1179

# Test range of mu_TGP
print(f"  b0 = {b0:.1f} (N_f={N_f})")
print()
print(f"  {'mu_TGP':>8s}  {'alpha_bare':>12s}  {'g0*_needed':>12s}  {'g0*/phi':>8s}  {'g0*/phi^2':>9s}")
print("  " + "-" * 58)

for mu_TGP in [80, 90, 100, 110, 120, 130, 150, 200, 300]:
    # Solve: 0.1179 = alpha_bare / [1 + b0*alpha_bare/(2pi)*ln(M_Z^2/mu^2)]
    # Let x = alpha_bare
    # 0.1179*(1 + b0*x/(2pi)*ln(M_Z^2/mu^2)) = x
    # 0.1179 + 0.1179*b0/(2pi)*ln(M_Z^2/mu^2)*x = x
    # 0.1179 = x*(1 - 0.1179*b0/(2pi)*ln(M_Z^2/mu^2))
    # x = 0.1179 / (1 - 0.1179*b0/(2pi)*ln(M_Z^2/mu^2))

    L = np.log(M_Z**2 / mu_TGP**2)
    alpha_bare = ALPHA_S_PDG / (1 - ALPHA_S_PDG * b0 / (2*np.pi) * L)
    g0_need = 4 * PHI0_B * alpha_bare / N_c**2
    print(f"  {mu_TGP:8.0f}  {alpha_bare:12.5f}  {g0_need:12.5f}  {g0_need/PHI:8.4f}  {g0_need/PHI**2:9.4f}")

# ===== 3. TEST: g0* = phi (golden ratio) =====
print("\n--- 3. Test: g0* = phi? ---\n")

alpha_s_phi_25 = N_c**2 * PHI / (4 * PHI0)
alpha_s_phi_B = N_c**2 * PHI / (4 * PHI0_B)

print(f"  alpha_s(g0*=phi) = N_c^2 * phi / (4 * Phi_0)")
print(f"    Phi_0 = 25:     alpha_s = {alpha_s_phi_25:.5f} ({(alpha_s_phi_25/ALPHA_S_PDG - 1)*100:+.2f}%)")
print(f"    Phi_0 = 24.783: alpha_s = {alpha_s_phi_B:.5f} ({(alpha_s_phi_B/ALPHA_S_PDG - 1)*100:+.2f}%)")
print()

# With running: what mu_TGP would make this work?
for label, Phi0, a_bare in [("Phi0=25", PHI0, alpha_s_phi_25), ("Phi0=B", PHI0_B, alpha_s_phi_B)]:
    # a_bare = a(M_Z)/(1 + b0*a(M_Z)/(2pi)*ln(M_Z^2/mu^2))
    # Actually a(M_Z) = a_bare/(1 + b0*a_bare/(2pi)*ln(M_Z^2/mu^2))
    # If a_bare > a(M_Z): mu > M_Z (TGP scale above M_Z)
    # ln(mu^2/M_Z^2) = (a_bare/a(M_Z) - 1) * 2*pi / (b0*a_bare)
    if a_bare > ALPHA_S_PDG:
        ln_ratio = (a_bare/ALPHA_S_PDG - 1) * 2*np.pi / (b0 * a_bare)
        mu_eff = M_Z * np.exp(ln_ratio/2)
        # Run to verify
        alpha_check = a_bare / (1 + b0*a_bare/(2*np.pi)*np.log(M_Z**2/mu_eff**2))
        print(f"  {label}: Jesli alpha_bare = {a_bare:.5f} przy mu = {mu_eff:.1f} GeV")
        print(f"           -> alpha_s(M_Z) = {alpha_check:.5f}")

# ===== 4. TEST: g0* jako wartosc uniwersalna z phi-FP =====
print("\n--- 4. Uniwersalne staleCCCC z phi-FP ---\n")

# In alpha=1 ODE, the phi-FP electron is at g0^e = 0.86901
# The phi-FP gives: g0^mu = phi * g0^e = 1.4054
# g0^tau = phi^2 * g0^e = 2.274 (approximately)
g0_e = 0.86901

# Key: g0^e depends on the ODE. With alpha=1: g0^e = 0.869
# With alpha=2 (old): g0^e = 1.2301
# The RATIO g0^mu/g0^e = phi is UNIVERSAL

# Question: is there a UNIVERSAL g0 (independent of particle) in alpha=1?
# Old g0* = 1.24915 was close to g0^e_old = 1.2301 (2% off)
# If g0* ~ g0^e, then with alpha=1: g0* ~ 0.869

# But alpha_s(g0^e=0.869) = 9*0.869/(4*24.783) = 0.0789. Too small!

# Alternative: maybe the formula changes. In alpha=1:
# K(g) = g^alpha = g (for alpha=1)  vs  g^2 (for alpha=2)
# The effective coupling might scale differently.

print(f"  g0^e(alpha=1) = {g0_e:.5f}")
print(f"  g0^mu(alpha=1) = {PHI*g0_e:.5f}")
print()
print(f"  alpha_s candidaci:")

candidates = [
    ("g0^e", g0_e),
    ("g0^mu = phi*g0^e", PHI*g0_e),
    ("phi", PHI),
    ("phi^2", PHI**2),
    ("1/phi", 1/PHI),
    ("phi*g0^e*phi (=phi^2*g0^e)", PHI**2*g0_e),
    ("sqrt(g0^e*g0^mu)", np.sqrt(g0_e * PHI*g0_e)),
    ("(g0^e + g0^mu)/2", (g0_e + PHI*g0_e)/2),
    ("g0^e * N_c/2", g0_e * N_c/2),
]

print(f"  {'kandydat':>30s}  {'g0*':>8s}  {'alpha_s':>8s}  {'dev%':>7s}")
print("  " + "-" * 60)
for name, val in candidates:
    a_s = N_c**2 * val / (4 * PHI0_B)
    dev = (a_s/ALPHA_S_PDG - 1)*100
    marker = " <--" if abs(dev) < 5 else ""
    print(f"  {name:>30s}  {val:8.5f}  {a_s:8.5f}  {dev:+7.2f}%{marker}")

# ===== 5. KLUCZOWY TEST: phi-formula bez g0* =====
print("\n--- 5. Formuly phi-bazowe (bez g0*) ---\n")

# Maybe alpha_s has a pure phi-formula?
# alpha_s = N_c^2 / (4 * Phi_0 * F(phi))

formulas = [
    ("N_c^2 * phi / (4*Phi0)", N_c**2 * PHI / (4*PHI0_B)),
    ("N_c^2 / (4*Phi0*phi)", N_c**2 / (4*PHI0_B*PHI)),
    ("N_c / (4*Phi0)", N_c / (4*PHI0_B)),  # just N_c instead of N_c^2
    ("N_c / (Phi0*phi^2)", N_c / (PHI0_B*PHI**2)),
    ("N_c^2 / (4*Phi0*phi^(1/2))", N_c**2 / (4*PHI0_B*PHI**0.5)),
    ("N_c^2 / (Phi0 * 4*phi^2)", N_c**2 / (PHI0_B * 4 * PHI**2)),
    ("1/(N_c * Phi0 * phi / pi)", 1/(N_c * PHI0_B * PHI / np.pi)),
    ("(N_c/Phi0)^2 * pi", (N_c/PHI0_B)**2 * np.pi),
    ("N_c^2 * phi^2 / (4*Phi0^(3/2))", N_c**2 * PHI**2 / (4*PHI0_B**1.5)),
    ("N_c^2/(4*Phi0) * phi/(phi+1)", N_c**2/(4*PHI0_B) * PHI/(PHI+1)),
    ("N_c^2/(4*Phi0) * 2/(phi+2)", N_c**2/(4*PHI0_B) * 2/(PHI+2)),
    ("N_c^2/(4*Phi0) * (phi-1/phi)", N_c**2/(4*PHI0_B) * (PHI - 1/PHI)),
]

print(f"  {'formula':>40s}  {'alpha_s':>8s}  {'dev%':>7s}")
print("  " + "-" * 62)
for name, val in formulas:
    dev = (val/ALPHA_S_PDG - 1)*100
    marker = " <--" if abs(dev) < 5 else " !!!" if abs(dev) < 1 else ""
    print(f"  {name:>40s}  {val:8.5f}  {dev:+7.2f}%{marker}")

# ===== 6. SOLITON COUPLING INTEGRAL =====
print("\n--- 6. Coupling integral z solitonu alpha=1 ---\n")

# In the derivation, J_c depends on the soliton SHAPE via an integral.
# J_c = N_c * <g> / Phi_0 where <g> is an effective average
# With different ODE, the soliton shape changes -> <g> changes
#
# Let's compute the "effective g" for the soliton:
# <g> = integral g(r) * weight(r) dr / integral weight(r) dr
# Natural weight: r^2 (volume element) or |delta|*r^2

print("  Soliton shape integrals dla g0^e = 0.86901, alpha=1:")
r, g = solve_ode(g0_e, r_max=80, alpha=1)
delta = g - 1.0

# Various effective <g> definitions
r2 = r**2
mask = r > 0.1  # avoid origin

# Volume-weighted average
g_vol = np.trapezoid(g[mask] * r2[mask], r[mask]) / np.trapezoid(r2[mask], r[mask])

# Delta-weighted (deviation from 1)
w_delta = np.abs(delta[mask]) * r2[mask]
g_delta = np.trapezoid(g[mask] * w_delta, r[mask]) / np.trapezoid(w_delta, r[mask])

# Energy-weighted
gp = np.gradient(g, r)
T = 0.5 * g**2 * gp**2
V = 0.5 * (1-g)**2
E_dens = (T + V) * r2
g_energy = np.trapezoid(g[mask] * E_dens[mask], r[mask]) / np.trapezoid(E_dens[mask], r[mask])

# Central value
g_center = g[0]

# Soliton "charge" Q = integral (1-g) * 4*pi*r^2 dr (topological)
Q = np.trapezoid((1-g[mask]) * 4*np.pi*r2[mask], r[mask])

# Coupling integral I = integral g^2 * (1-g) * 4*pi*r^2 dr
I_coupling = np.trapezoid(g[mask]**2 * (1-g[mask]) * 4*np.pi*r2[mask], r[mask])

print(f"  g0 (central)     = {g_center:.5f}")
print(f"  <g> (vol-avg)    = {g_vol:.5f}")
print(f"  <g> (delta-wgt)  = {g_delta:.5f}")
print(f"  <g> (energy-wgt) = {g_energy:.5f}")
print(f"  Q (topological)  = {Q:.4f}")
print(f"  I_coupling       = {I_coupling:.4f}")
print()

# What alpha_s do these give?
print(f"  alpha_s z roznych <g>:")
for name, geff in [("g0", g_center), ("<g>_vol", g_vol), ("<g>_delta", g_delta),
                    ("<g>_energy", g_energy), ("Q/(4pi)", Q/(4*np.pi)),
                    ("I_c/(4pi)", I_coupling/(4*np.pi))]:
    a_s = N_c**2 * geff / (4 * PHI0_B)
    dev = (a_s/ALPHA_S_PDG - 1)*100
    marker = " <--" if abs(dev) < 10 else ""
    print(f"    {name:>15s}: g_eff = {geff:8.5f}, alpha_s = {a_s:8.5f} ({dev:+.1f}%){marker}")

# ===== 7. NOWA HIPOTEZA: alpha_s z relacji phi-FP =====
print("\n--- 7. Nowa hipoteza: alpha_s z phi-FP relacji ---\n")

# The phi-FP gives r21 = (A_mu/A_e)^4.
# The amplitude A_e encodes the coupling strength.
# Maybe alpha_s is related to A_e itself?

A_e = get_Atail(g0_e, alpha=1)
A_mu = get_Atail(PHI * g0_e, alpha=1)

print(f"  A_tail(g0^e)  = {A_e:.6f}")
print(f"  A_tail(g0^mu) = {A_mu:.6f}")
print(f"  r21 = {(A_mu/A_e)**4:.2f}")
print()

# Try: alpha_s = A_e * something?
# Or: alpha_s = A_e^2 * something?
# Dimensional analysis: alpha_s is dimensionless, A_e has dimension of mass^(1/4) in TGP

# Actually: m_e = c_M * A_e^4
# alpha_em = e^2/(4pi) ~ 1/137
# alpha_s = N_c * g_s^2/(4pi)
# Maybe: alpha_s / alpha_em = function(A_e, phi)?

print(f"  Relacje z A_tail:")
alpha_em = 1/137.036

print(f"    A_e^2 = {A_e**2:.6f}")
print(f"    A_e^4 = {A_e**4:.8f}")
print(f"    A_e^2 * N_c = {A_e**2 * N_c:.6f}")
print(f"    alpha_s/alpha_em = {ALPHA_S_PDG/alpha_em:.4f}")
print(f"    A_mu/A_e = {A_mu/A_e:.6f}")
print(f"    (A_mu/A_e)^2 = {(A_mu/A_e)**2:.4f}")
print()

# ===== 8. KLUCZOWA OBSERWACJA: g0* z starego ODE w nowym ODE =====
print("\n--- 8. Stare g0* w nowym kontekscie ---\n")

# The OLD g0* = 1.24915 was from the OLD ODE (parametric, f(g)=1+2*alpha*ln(g))
# NOT from alpha=2 canonical ODE!
#
# In the old parametric ODE:
#   g'' + f(g)/g * (g')^2 + (2/r)*g' = g*f(g)*(1-g)
# where f(g) = 1 + 2*alpha_kin*ln(g)
#
# The B_tail=0 condition selected g0* in THAT ODE.
# In alpha=1 substrate ODE, the soliton is different.
#
# Question: Is g0*=1.249 a solution property of the OLD ODE only,
# or does it have a GEOMETRIC meaning independent of ODE?

# If g0* ~ 5/4 = 1.25, this suggests a simple fraction.
# Or g0* ~ phi/(phi-1+1/phi) = phi/phi = 1. No.
# g0* ~ (3+phi)/phi^2 = 4.618/2.618 = 1.764. No.
# g0* ~ phi^(1/phi) = 1.618^0.618 = ?
phi_to_1_over_phi = PHI**(1/PHI)
# g0* ~ phi^(phi-1) = phi^(1/phi)
phi_to_phi_minus_1 = PHI**(PHI-1)

print(f"  g0*(stare) = 1.24915")
print(f"  5/4 = {5/4:.5f} (rozn: {(5/4 - 1.24915)/1.24915*100:+.3f}%)")
print(f"  phi^(1/phi) = {phi_to_1_over_phi:.5f} (rozn: {(phi_to_1_over_phi - 1.24915)/1.24915*100:+.3f}%)")
print(f"  phi^(phi-1) = {phi_to_phi_minus_1:.5f} (rozn: {(phi_to_phi_minus_1 - 1.24915)/1.24915*100:+.3f}%)")
print(f"  phi/(1+1/phi) = {PHI/(1+1/PHI):.5f}")
print(f"  2*phi-phi^2 = {2*PHI-PHI**2:.5f}")  # = 2phi - phi-1 = phi-1 = 1/phi

# ===== 9. ROZWIAZANIE: running z g0^mu =====
print("\n--- 9. Najlepszy kandydat: g0^mu z running ---\n")

# g0^mu = phi * g0^e = 1.4054
# alpha_s(g0^mu) = 9 * 1.4054 / (4*24.783) = 0.1276
# This is 8% above PDG.
# With 1-loop running from mu ~ 200 GeV to M_Z:

g0_mu = PHI * g0_e
alpha_bare_mu = N_c**2 * g0_mu / (4 * PHI0_B)

print(f"  g0^mu = phi * g0^e = {g0_mu:.5f}")
print(f"  alpha_bare(g0^mu) = {alpha_bare_mu:.5f}")
print()

# Find mu_TGP where this bare value runs to 0.1179:
if alpha_bare_mu > ALPHA_S_PDG:
    ln_ratio = (alpha_bare_mu/ALPHA_S_PDG - 1) * 2*np.pi / (b0 * alpha_bare_mu)
    mu_eff = M_Z * np.exp(ln_ratio/2)
    alpha_check = alpha_bare_mu / (1 + b0*alpha_bare_mu/(2*np.pi)*np.log(M_Z**2/mu_eff**2))
    print(f"  mu_TGP = {mu_eff:.1f} GeV (gdyby alpha_bare = alpha_s(mu_TGP))")
    print(f"  alpha_s(M_Z) z running = {alpha_check:.5f}")
    print(f"  Weryfikacja: {abs(alpha_check - ALPHA_S_PDG) < 1e-4}")
else:
    print(f"  alpha_bare < alpha_PDG -> mu_TGP < M_Z")

# ===== 10. SYSTEMATYCZNY SKAN: g0_eff z roznych mu_TGP =====
print("\n--- 10. Jaki g0_eff daje alpha_s(M_Z)=0.1179 po running? ---\n")

print("  Logika: alpha_bare = N_c^2*g0_eff/(4*Phi0)")
print("          alpha_s(M_Z) = running(alpha_bare, mu_TGP -> M_Z)")
print()
print(f"  {'mu_TGP':>8s}  {'alpha_bare':>12s}  {'g0_eff':>10s}  {'g0_eff/phi':>10s}  {'g0_eff/g0e':>10s}")
print("  " + "-" * 60)

for mu in [M_Z, 100, 120, 150, 200, 300, 500, 1000]:
    L = np.log(M_Z**2 / mu**2)
    if abs(L) < 1e-10:
        a_bare = ALPHA_S_PDG
    else:
        a_bare = ALPHA_S_PDG / (1 - ALPHA_S_PDG * b0 / (2*np.pi) * L)
    g0_eff = 4 * PHI0_B * a_bare / N_c**2
    print(f"  {mu:8.1f}  {a_bare:12.5f}  {g0_eff:10.5f}  {g0_eff/PHI:10.5f}  {g0_eff/g0_e:10.5f}")

# ===== WNIOSKI =====
print("\n" + "=" * 72)
print("WNIOSKI ex177")
print("=" * 72)
print(f"""
  1. PROBLEM: B_tail=0 (H1) nie istnieje w alpha=1 ODE.
     g0* = 1.24915 bylo specyficzne dla starego parametrycznego ODE.

  2. REVERSE ENGINEERING:
     Aby alpha_s = 0.1179, potrzebne g0* = {g0_needed_B:.4f} (Phi0=24.783)
     Stare g0* = 1.2492 (roznica {(g0_needed_B - 1.24915)/1.24915*100:+.2f}%)

  3. RUNNING ARGUMENT:
     Jesli alpha_s^TGP jest wartoscia BARE przy skali mu_TGP,
     to 1-loop running do M_Z daje alpha_s(M_Z) = 0.1179.
     - g0^mu = phi*g0^e = {g0_mu:.4f}: bare = {alpha_bare_mu:.4f},
       mu_TGP ~ {mu_eff:.0f} GeV
     - Stare g0* = 1.249: bare = 0.1134, mu_TGP ~ 120 GeV (ex144)

  4. INTERPRETACJA:
     Formula alpha_s moze uzyc g0_eff ktore ZALEZY od skali mu_TGP.
     Problem G2 rozwiazuje sie jesli:
     a) g0_eff = g0^mu (muon phi-FP) przy mu_TGP ~ {mu_eff:.0f} GeV, LUB
     b) g0_eff = 1.249 (stara wartosc) przy mu_TGP ~ 120 GeV, LUB
     c) Nowa derywacja z alpha=1 daje inna formule

  5. STATUS G2:
     CZESCIOWO ROZWIAZANY: running argument zachowuje spojnosc.
     Otwarte: FIZYCZNA motywacja g0_eff i mu_TGP z alpha=1 ODE.
     Stare g0* = 1.249 moze byc ARTEFAKTEM starego ODE.

  6. REKOMENDACJA:
     Przeformulowac alpha_s jako alpha_s = f(N_c, Phi_0, phi, mu_TGP/M_Z).
     g0* moze nie byc fundamentalne -- moze byc jedynie parametrem
     pasujacym do konkretnej formy ODE.
""")
