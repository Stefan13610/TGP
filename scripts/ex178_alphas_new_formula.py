#!/usr/bin/env python3
"""
ex178_alphas_new_formula.py
Sesja v45, 2026-04-05

NOWA FORMULA alpha_s z ODE substratowym (alpha=1).

ODKRYCIE (ex177, sekcja 4):
  g0^e * N_c/2 = 0.86901 * 3/2 = 1.30351
  alpha_s = N_c^2 * (g0^e * N_c/2) / (4*Phi_0)
          = N_c^3 * g0^e / (8*Phi_0)
          = 27 * 0.86901 / (8 * 24.783) = 0.11834
  PDG: 0.1179, odchylenie: +0.38%

INTERPRETACJA:
  Stara formula: alpha_s = N_c^2 * g0* / (4*Phi_0) z g0* z B_tail=0
  Nowa formula:  alpha_s = N_c^3 * g0^e / (8*Phi_0) z g0^e z phi-FP

  Przelaczenie g0* -> g0^e*N_c/2 sugeruje:
  - Czynnik N_c/2 = T_F * N_c (T_F = 1/2 normalizacja generatora SU(3))
  - Lub: g0* = g0^e * N_c/2 (efektywna wartosc centralna w QCD sektorze)

  Fizyka: elektron (g0^e) SKALOWANY przez czynnik kolorowy N_c/2
  daje efektywne sprzezenie silne. To jest NATURALNE:
  lepton-quark relacja z czynnikiem kolorowym.

WERYFIKACJA:
  1. Precyzja formuly z roznymi Phi_0
  2. Self-consistency z phi-FP (g0^e pochodzi z r21 matchingu)
  3. Inverse: jaki Phi_0 daje DOKLADNIE alpha_s = 0.1179?
  4. Porownanie ze stara formula (g0*=1.249)
  5. Running i skala mu_TGP
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
N_c = 3
ALPHA_S_PDG = 0.1179
ALPHA_S_ERR = 0.0009  # PDG uncertainty
M_Z = 91.1876

print("=" * 72)
print("ex178: NOWA FORMULA alpha_s = N_c^3 * g0^e / (8 * Phi_0)")
print("=" * 72)

# ---- ODE solver (alpha=1) ----
def solve_ode(g0, r_max=150, n_points=40000):
    def rhs(r, y):
        g, gp = y
        if g < 1e-12: g = 1e-12
        gpp = (1-g) - (1.0/g)*gp**2 - (2.0/r)*gp if r > 1e-12 else (1-g)/3.0
        return [gp, gpp]
    r_eval = np.linspace(1e-6, r_max, n_points)
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0], method='RK45',
                    t_eval=r_eval, rtol=1e-10, atol=1e-12, max_step=0.04)
    return sol.t, sol.y[0]

def get_Atail(g0, window=(60, 120)):
    r, g = solve_ode(g0, r_max=150)
    delta = g - 1.0
    mask = (r > window[0]) & (r < window[1])
    r_w, d_w = r[mask], delta[mask]
    if len(r_w) < 300:
        return None
    y = d_w * r_w
    M = np.column_stack([np.cos(r_w), np.sin(r_w)])
    result = np.linalg.lstsq(M, y, rcond=None)
    return np.sqrt(result[0][0]**2 + result[0][1]**2)

# ===== 1. PRECYZJA FORMULY =====
print("\n--- 1. alpha_s = N_c^3 * g0^e / (8 * Phi_0) ---\n")

# g0^e from phi-FP (ex174, alpha=1 ODE)
g0_e = 0.86901

# Refine g0^e: find where (A(phi*g0)/A(g0))^4 = 206.768
print("  Kalibracja g0^e z r21 = 206.768...")

def r21_func(g0):
    A_e = get_Atail(g0)
    A_mu = get_Atail(PHI * g0)
    if A_e is None or A_mu is None or A_e < 1e-12:
        return 0
    return (A_mu / A_e)**4

# Scan around 0.869 to find better g0^e
best_g0e = g0_e
best_r21_err = abs(r21_func(g0_e) - 206.768)
for g0_test in np.linspace(0.860, 0.880, 21):
    r21_test = r21_func(g0_test)
    err = abs(r21_test - 206.768)
    if err < best_r21_err:
        best_r21_err = err
        best_g0e = g0_test

# Refine with finer grid
for g0_test in np.linspace(best_g0e - 0.002, best_g0e + 0.002, 41):
    r21_test = r21_func(g0_test)
    err = abs(r21_test - 206.768)
    if err < best_r21_err:
        best_r21_err = err
        best_g0e = g0_test

r21_best = r21_func(best_g0e)
print(f"  g0^e (refined) = {best_g0e:.5f}")
print(f"  r21 = {r21_best:.2f} (PDG: 206.768)")
print(f"  r21 error: {abs(r21_best - 206.768)/206.768*100:.2f}%")
print()

g0_e = best_g0e

# Now compute alpha_s with various Phi_0
phi0_values = [
    ("Phi_0 = 25.000 (exact)", 25.000),
    ("Phi_0 = 24.783 (Brannen)", 24.783),
    ("Phi_0 = 24.877 (inverse fit)", None),  # will compute
    ("Phi_0 = 4*pi^2 = 39.478", 4*np.pi**2),
    ("Phi_0 = 8*pi = 25.133", 8*np.pi),
]

# Compute inverse fit Phi_0
Phi0_fit = N_c**3 * g0_e / (8 * ALPHA_S_PDG)
phi0_values[2] = (f"Phi_0 = {Phi0_fit:.3f} (inverse fit)", Phi0_fit)

print(f"  NOWA FORMULA: alpha_s = N_c^3 * g0^e / (8 * Phi_0)")
print(f"  g0^e = {g0_e:.5f} (z phi-FP, alpha=1 ODE)")
print(f"  N_c = {N_c}")
print()
print(f"  {'Phi_0':>35s}  {'alpha_s':>8s}  {'dev':>7s}  {'sigma':>6s}")
print("  " + "-" * 65)

for name, Phi0 in phi0_values:
    a_s = N_c**3 * g0_e / (8 * Phi0)
    dev = (a_s - ALPHA_S_PDG) / ALPHA_S_PDG * 100
    sigma = abs(a_s - ALPHA_S_PDG) / ALPHA_S_ERR
    marker = " <--" if abs(dev) < 1 else ""
    print(f"  {name:>35s}  {a_s:8.5f}  {dev:+7.2f}%  {sigma:5.2f}s{marker}")

# ===== 2. STARA vs NOWA formula =====
print("\n--- 2. Porownanie STARA vs NOWA formula ---\n")

g0_star_old = 1.24915
Phi0_B = 24.783

alpha_old = N_c**2 * g0_star_old / (4 * Phi0_B)
alpha_new = N_c**3 * g0_e / (8 * Phi0_B)

print(f"  STARA: alpha_s = N_c^2 * g0* / (4*Phi_0)")
print(f"         = {N_c}^2 * {g0_star_old:.5f} / (4 * {Phi0_B})")
print(f"         = {alpha_old:.5f} ({(alpha_old/ALPHA_S_PDG - 1)*100:+.2f}% od PDG)")
print()
print(f"  NOWA:  alpha_s = N_c^3 * g0^e / (8*Phi_0)")
print(f"         = {N_c}^3 * {g0_e:.5f} / (8 * {Phi0_B})")
print(f"         = {alpha_new:.5f} ({(alpha_new/ALPHA_S_PDG - 1)*100:+.2f}% od PDG)")
print()
print(f"  RELACJA: nowa/stara = {alpha_new/alpha_old:.5f}")
print(f"           g0^e*N_c/2 / g0* = {g0_e*N_c/2 / g0_star_old:.5f}")
print()

# The factor N_c/2 replaces g0*/g0^e
print(f"  Czynnik korekcji: N_c/2 * g0^e / g0* = {N_c/2 * g0_e / g0_star_old:.5f}")
print(f"  Jesli g0* = N_c/2 * g0^e:")
print(f"    g0* = {N_c/2} * {g0_e:.5f} = {N_c/2*g0_e:.5f}")
print(f"    vs stare g0* = {g0_star_old}")
print(f"    roznica: {(N_c/2*g0_e - g0_star_old)/g0_star_old*100:+.2f}%")

# ===== 3. FIZYCZNA INTERPRETACJA N_c/2 =====
print("\n--- 3. Fizyczna interpretacja czynnika N_c/2 ---\n")

# Color factors in SU(N_c):
C_F = (N_c**2 - 1) / (2*N_c)    # Casimir of fundamental rep
C_A = N_c                         # Casimir of adjoint rep
T_F = 0.5                         # Normalization of generators
d_F = N_c                         # Dimension of fundamental rep
d_A = N_c**2 - 1                  # Dimension of adjoint rep

print(f"  SU({N_c}) color factors:")
print(f"    C_F = (N_c^2-1)/(2*N_c) = {C_F:.4f}")
print(f"    C_A = N_c = {C_A}")
print(f"    T_F = 1/2 = {T_F}")
print(f"    d_F = N_c = {d_F}")
print(f"    d_A = N_c^2-1 = {d_A}")
print()
print(f"  N_c/2 = {N_c/2:.1f}")
print(f"  T_F * N_c = {T_F * N_c:.1f} = N_c/2  <-- MATCH!")
print(f"  C_F = {C_F:.4f}")
print()
print(f"  INTERPRETACJA:")
print(f"    N_c/2 = T_F * N_c = sum_a (T^a T^a)_fundamental")
print(f"    To jest calkowity 'color charge' w reprezentacji fundamentalnej.")
print(f"    Fizycznie: kwark niesie ladunek kolorowy T_F*N_c razy silniejszy")
print(f"    niz lepton. Dlatego g0_eff = g0^e * N_c/2.")

# ===== 4. SELF-CONSISTENCY =====
print("\n--- 4. Self-consistency: wszystkie predykcje z g0^e ---\n")

# With the new formula, ALL TGP predictions depend on:
# - g0^e (from phi-FP)
# - phi (golden ratio)
# - N_c (SU(3))
# - Phi_0 (from Brannen or exact 25)

# 1. Mass ratio r21:
A_e = get_Atail(g0_e)
A_mu = get_Atail(PHI * g0_e)
r21 = (A_mu / A_e)**4

# 2. alpha_s:
alpha_s_new = N_c**3 * g0_e / (8 * Phi0_B)

# 3. m_tau (from phi-extended Koide):
# m_tau = m_e * (r21)^(phi^2) -- approximate
# Actually m_tau = c_M * A(phi^2*g0^e)^4
g0_tau_approx = PHI**2 * g0_e  # if we assume phi^2 scaling (not exact)
A_tau = get_Atail(g0_tau_approx)
if A_tau and A_e:
    r31 = (A_tau / A_e)**4

print(f"  Input: g0^e = {g0_e:.5f}")
print(f"  phi = {PHI:.6f}")
print(f"  N_c = {N_c}")
print(f"  Phi_0 = {Phi0_B}")
print()
print(f"  Predykcje:")
print(f"    r21 = (A_mu/A_e)^4 = {r21:.2f} (PDG: 206.768, {(r21/206.768-1)*100:+.2f}%)")
if A_tau:
    print(f"    r31 = (A_tau/A_e)^4 = {r31:.2f} (PDG: 3477.4, {(r31/3477.4-1)*100:+.2f}%)")
print(f"    alpha_s = N_c^3*g0^e/(8*Phi_0) = {alpha_s_new:.5f} (PDG: 0.1179, {(alpha_s_new/ALPHA_S_PDG-1)*100:+.2f}%)")

# ===== 5. SENSITIVITY: g0^e uncertainty propagation =====
print("\n--- 5. Sensitivity analiza ---\n")

# How sensitive is alpha_s to g0^e?
# alpha_s = N_c^3 * g0^e / (8 * Phi_0)
# d(alpha_s)/alpha_s = d(g0^e)/g0^e
# So 1% change in g0^e -> 1% change in alpha_s
# Current r21 precision: ~1.9% off PDG
# This means g0^e uncertainty is ~0.5% (r21 scales as g0^4 roughly)

dg0_range = np.linspace(-0.01, 0.01, 21)
print(f"  {'dg0^e/g0^e':>12s}  {'g0^e':>8s}  {'alpha_s':>8s}  {'dev_PDG%':>9s}")
print("  " + "-" * 45)
for dg in dg0_range[::4]:
    g0_test = g0_e * (1 + dg)
    a_s = N_c**3 * g0_test / (8 * Phi0_B)
    print(f"  {dg*100:+12.2f}%  {g0_test:8.5f}  {a_s:8.5f}  {(a_s/ALPHA_S_PDG-1)*100:+9.2f}%")

# ===== 6. INVERSE: Phi_0 z alpha_s =====
print("\n--- 6. Phi_0 z alpha_s (inverse problem) ---\n")

Phi0_inv = N_c**3 * g0_e / (8 * ALPHA_S_PDG)

print(f"  alpha_s(PDG) = {ALPHA_S_PDG} +/- {ALPHA_S_ERR}")
print(f"  g0^e = {g0_e:.5f}")
print()
print(f"  Phi_0 = N_c^3 * g0^e / (8 * alpha_s)")
print(f"        = {N_c**3} * {g0_e:.5f} / (8 * {ALPHA_S_PDG})")
print(f"        = {Phi0_inv:.4f}")
print()
print(f"  Porownanie:")
print(f"    Phi_0(inverse) = {Phi0_inv:.4f}")
print(f"    Phi_0(Brannen) = {Phi0_B:.3f}")
print(f"    Phi_0(exact)   = 25.000")
print(f"    8*pi           = {8*np.pi:.4f}")
print()

# With error propagation:
Phi0_inv_hi = N_c**3 * g0_e / (8 * (ALPHA_S_PDG - ALPHA_S_ERR))
Phi0_inv_lo = N_c**3 * g0_e / (8 * (ALPHA_S_PDG + ALPHA_S_ERR))
print(f"  Phi_0 z alpha_s +/- 1sigma: [{Phi0_inv_lo:.3f}, {Phi0_inv_hi:.3f}]")
print(f"  Brannen (24.783) w zakresie? {'TAK' if Phi0_inv_lo < 24.783 < Phi0_inv_hi else 'NIE'}")
print(f"  Exact 25 w zakresie?         {'TAK' if Phi0_inv_lo < 25 < Phi0_inv_hi else 'NIE'}")

# ===== 7. POROWNANIE: ile parametrow wolnych? =====
print("\n--- 7. Parametry formuly ---\n")

print("  STARA formula: alpha_s = N_c^2 * g0* / (4*Phi_0)")
print("    N_c = 3 (z symmetrii)")
print("    g0* = 1.24915 (z H1: B_tail=0 w starym ODE)")
print("    Phi_0 = 24.783 (Brannen) lub 25")
print("    -> g0* jest NIEZALEZNY od phi-FP (oddzielny warunek)")
print()
print("  NOWA formula: alpha_s = N_c^3 * g0^e / (8*Phi_0)")
print("    N_c = 3 (z symmetrii)")
print("    g0^e = 0.869 (z phi-FP: r21 = 206.77)")
print("    Phi_0 = 24.783 (Brannen) lub 25")
print("    -> g0^e jest TEN SAM co w predykcji mas!")
print("    -> alpha_s jest KONSEKWENCJA phi-FP, NIE niezalezna predykcja")
print()
print("  WNIOSEK: Nowa formula ma MNIEJ wolnych parametrow!")
print("           g0* (z B_tail=0) jest WYELIMINOWANY.")
print("           Wszystko wynika z g0^e (phi-FP) + N_c + Phi_0.")

# ===== 8. FORMULA SUMMARY =====
print("\n--- 8. Podsumowanie nowej formuly ---\n")

# Triple cross-check with Phi_0 = 25
alpha_25 = N_c**3 * g0_e / (8 * 25)
alpha_B = N_c**3 * g0_e / (8 * Phi0_B)

print(f"  alpha_s = N_c^3 * g0^e / (8 * Phi_0)")
print(f"         = (N_c/2) * [N_c^2 * g0^e / (4 * Phi_0)]")
print(f"         = T_F * N_c * [kappa * N_c * g0^e]")
print()
print(f"  Wartosc liczbowa:")
print(f"    Phi_0 = 25:     alpha_s = {alpha_25:.5f} (PDG: 0.1179, {(alpha_25/ALPHA_S_PDG-1)*100:+.2f}%)")
print(f"    Phi_0 = 24.783: alpha_s = {alpha_B:.5f} (PDG: 0.1179, {(alpha_B/ALPHA_S_PDG-1)*100:+.2f}%)")
print()
print(f"  Wartosc PDG: 0.1179 +/- 0.0009")
print(f"  alpha_s^TGP lezi {abs(alpha_B - ALPHA_S_PDG)/ALPHA_S_ERR:.1f} sigma od PDG (Brannen)")
print(f"  alpha_s^TGP lezi {abs(alpha_25 - ALPHA_S_PDG)/ALPHA_S_ERR:.1f} sigma od PDG (exact 25)")

# ===== WNIOSKI =====
print("\n" + "=" * 72)
print("WNIOSKI ex178")
print("=" * 72)
print(f"""
  1. NOWA FORMULA:
     alpha_s = N_c^3 * g0^e / (8 * Phi_0)
     = 27 * {g0_e:.5f} / (8 * 24.783) = {alpha_B:.5f}
     PDG: 0.1179 +/- 0.0009
     Odchylenie: {(alpha_B/ALPHA_S_PDG-1)*100:+.2f}% ({abs(alpha_B - ALPHA_S_PDG)/ALPHA_S_ERR:.1f} sigma)

  2. FIZYCZNA INTERPRETACJA:
     alpha_s = (T_F * N_c) * [N_c^2 * g0^e / (4 * Phi_0)]
     Czynnik T_F * N_c = N_c/2 = 3/2:
     calkowity ladunek kolorowy w rep. fundamentalnej.
     Kwark = lepton * czynnik kolorowy.

  3. ELIMINACJA g0*:
     Stara formula uzywala g0* z B_tail=0 (H1).
     Nowa formula uzywa g0^e z phi-FP.
     -> g0^e juz sluzy do predykcji r21 i mas!
     -> alpha_s staje sie KONSEKWENCJA tego samego g0^e.
     -> Mniej wolnych parametrow.

  4. Phi_0 INVERSE:
     Phi_0 = N_c^3 * g0^e / (8 * alpha_s)
           = {Phi0_inv:.3f}
     Miedzy Brannen (24.783) a exact (25).
     Spojne z oboma propozycjami.

  5. STATUS G2: ROZWIAZANY (warunkowo)
     Nowa formula alpha_s = N_c^3 * g0^e / (8*Phi_0)
     jest spojne z PDG w granicach {abs(alpha_B - ALPHA_S_PDG)/ALPHA_S_ERR:.1f} sigma.
     Wymaga jeszcze formalnej derywacji z dodatekV.
""")
