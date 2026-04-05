#!/usr/bin/env python3
"""
ex181_alpha_em_connection.py
Sesja v45, 2026-04-05

Badanie polaczenia alpha_em z ramami TGP.

MOTYWACJA:
  alpha_s = N_c^3 * g0^e / (8*Phi_0) laczy lepton z QCD.
  Czy alpha_em = f(g0^e, Phi_0, phi, N_c) jest rowniez predykowalne?

  alpha_em = 1/137.036 = 0.007297 (PDG, najdokladniej zmierzona stala)
  alpha_s = 0.1179 (PDG)

OBSERWACJA (ex179):
  alpha_s / alpha_em = 16.157
  Bliskie 10*phi = 16.180 (0.14% roznica)

HIPOTEZY:
  H1: alpha_em = g0^e^2 / (4*Phi_0)  [kwadrat g0^e]
  H2: alpha_em = alpha_s / (10*phi)   [relacja phi-bazowa]
  H3: alpha_em = g0^e / (Phi_0 * C)   [nowy czynnik C]
  H4: alpha_em = 1/(4*pi*Phi_0*F)     [relacja geometryczna]
"""
import numpy as np
from scipy.integrate import solve_ivp

PHI = (1 + np.sqrt(5)) / 2
N_c = 3
ALPHA_S_PDG = 0.1179
ALPHA_EM_PDG = 1/137.036
PHI0 = 25.0

print("=" * 72)
print("ex181: Polaczenie alpha_em z TGP")
print("=" * 72)

# ODE solver
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

g0_e = 0.86941  # from ex180 refined

# ===== 1. RELACJA alpha_s / alpha_em =====
print("\n--- 1. Stosunek alpha_s / alpha_em ---\n")

ratio = ALPHA_S_PDG / ALPHA_EM_PDG
print(f"  alpha_s = {ALPHA_S_PDG}")
print(f"  alpha_em = 1/137.036 = {ALPHA_EM_PDG:.7f}")
print(f"  alpha_s / alpha_em = {ratio:.4f}")
print()

# Test closed forms
ratio_forms = [
    ("10*phi", 10*PHI),
    ("N_c^3 / (2*g0^e)", N_c**3 / (2*g0_e)),
    ("2*Phi_0*phi/5", 2*PHI0*PHI/5),  # = 10*phi if Phi_0=25
    ("4*pi", 4*np.pi),
    ("N_c^2*phi", N_c**2*PHI),
    ("Phi_0*phi/pi", PHI0*PHI/np.pi),
    ("phi^6", PHI**6),
    ("2*N_c^2", 2*N_c**2),
    ("N_c^3", N_c**3),
    ("8*phi^2", 8*PHI**2),
    ("phi^3 * N_c", PHI**3 * N_c),
    ("5*phi^3", 5*PHI**3),
]

ratio_forms_sorted = sorted(ratio_forms, key=lambda x: abs(x[1] - ratio))
print(f"  {'formula':>25s}  {'value':>10s}  {'dev%':>8s}")
print("  " + "-" * 50)
for name, val in ratio_forms_sorted[:8]:
    dev = (val - ratio)/ratio * 100
    marker = " <--" if abs(dev) < 1 else ""
    print(f"  {name:>25s}  {val:10.4f}  {dev:+8.3f}%{marker}")

print()
print(f"  NAJLEPSZA: alpha_s/alpha_em ~ 10*phi = {10*PHI:.4f} (dev: {(10*PHI - ratio)/ratio*100:+.3f}%)")

# ===== 2. HIPOTEZA H1: alpha_em = g0^e^2 / (4*Phi_0) =====
print("\n--- 2. H1: alpha_em = g0^e^2 / (4*Phi_0) ---\n")

alpha_em_H1 = g0_e**2 / (4*PHI0)
print(f"  alpha_em(H1) = {g0_e:.5f}^2 / (4*{PHI0:.0f}) = {alpha_em_H1:.7f}")
print(f"  alpha_em(PDG) = {ALPHA_EM_PDG:.7f}")
print(f"  Odchylenie: {(alpha_em_H1/ALPHA_EM_PDG - 1)*100:+.2f}%")
print(f"  1/alpha_em(H1) = {1/alpha_em_H1:.2f} (PDG: 137.036)")

# If H1 true, then alpha_s/alpha_em = N_c^3*g0^e/(8*Phi_0) / (g0^e^2/(4*Phi_0))
#                                     = N_c^3 / (2*g0^e)
r_H1 = N_c**3 / (2*g0_e)
print(f"\n  Z H1: alpha_s/alpha_em = N_c^3/(2*g0^e) = {r_H1:.4f}")
print(f"  Obserwacja:               = {ratio:.4f}")
print(f"  Odchylenie: {(r_H1/ratio - 1)*100:+.2f}%")

# ===== 3. HIPOTEZA H2: alpha_em = alpha_s / (10*phi) =====
print("\n--- 3. H2: alpha_em = alpha_s / (10*phi) ---\n")

alpha_em_H2 = ALPHA_S_PDG / (10*PHI)
print(f"  alpha_em(H2) = {ALPHA_S_PDG} / (10*phi) = {alpha_em_H2:.7f}")
print(f"  alpha_em(PDG) = {ALPHA_EM_PDG:.7f}")
print(f"  Odchylenie: {(alpha_em_H2/ALPHA_EM_PDG - 1)*100:+.3f}%")
print(f"  1/alpha_em(H2) = {1/alpha_em_H2:.2f} (PDG: 137.036)")

# If H2 true, combined with alpha_s formula:
# alpha_em = N_c^3 * g0^e / (80*Phi_0*phi)
alpha_em_H2_full = N_c**3 * g0_e / (80*PHI0*PHI)
print(f"\n  Polaczenie z alpha_s formula:")
print(f"  alpha_em = N_c^3 * g0^e / (80*Phi_0*phi)")
print(f"           = 27*{g0_e:.5f} / (80*{PHI0:.0f}*{PHI:.4f})")
print(f"           = {alpha_em_H2_full:.7f}")
print(f"  1/alpha_em = {1/alpha_em_H2_full:.2f}")
print(f"  Odch od PDG: {(alpha_em_H2_full/ALPHA_EM_PDG - 1)*100:+.3f}%")

# ===== 4. HIPOTEZA H3: alpha_em z A_tail =====
print("\n--- 4. H3: alpha_em z amplitudy ogona ---\n")

A_e = get_Atail(g0_e)
print(f"  A_tail(g0^e) = {A_e:.6f}")
print(f"  A_e^2 = {A_e**2:.6f}")
print(f"  A_e^4 = {A_e**4:.8f}")
print()

# Test: alpha_em = A_e^2 * C for some C?
C_needed = ALPHA_EM_PDG / A_e**2
print(f"  alpha_em = A_e^2 * C => C = {C_needed:.5f}")
print(f"  C ~ {C_needed:.3f}")

# Test: alpha_em = A_e^2 / (2*phi)?
alpha_em_Ae = A_e**2 / (2*PHI)
print(f"\n  alpha_em = A_e^2/(2*phi) = {alpha_em_Ae:.7f}")
print(f"  Odch: {(alpha_em_Ae/ALPHA_EM_PDG - 1)*100:+.2f}%")

# alpha_em = A_e^2 / 2?
alpha_em_Ae2 = A_e**2 / 2
print(f"  alpha_em = A_e^2/2 = {alpha_em_Ae2:.7f}")
print(f"  Odch: {(alpha_em_Ae2/ALPHA_EM_PDG - 1)*100:+.2f}%")

# alpha_em = A_e^2 * pi / 3?
alpha_em_Ae3 = A_e**2 * np.pi / 3
print(f"  alpha_em = A_e^2*pi/3 = {alpha_em_Ae3:.7f}")
print(f"  Odch: {(alpha_em_Ae3/ALPHA_EM_PDG - 1)*100:+.2f}%")

# ===== 5. SYSTEMATYCZNE PRZESZUKANIE =====
print("\n--- 5. Systematyczne przeszukanie formul alpha_em ---\n")

# alpha_em = f(g0^e, Phi_0, phi, N_c, A_e)
formulas = [
    ("g0^e^2 / (4*Phi_0)", g0_e**2 / (4*PHI0)),
    ("g0^e^2 / (4*Phi_0*phi^(1/3))", g0_e**2 / (4*PHI0*PHI**(1/3))),
    ("g0^e / (Phi_0*phi^3)", g0_e / (PHI0*PHI**3)),
    ("g0^e / (N_c * Phi_0 * phi)", g0_e / (N_c * PHI0 * PHI)),
    ("g0^e^3 / (4*Phi_0*phi)", g0_e**3 / (4*PHI0*PHI)),
    ("1/(4*pi*Phi_0/phi)", 1/(4*np.pi*PHI0/PHI)),
    ("1/(Phi_0 * N_c * phi^2)", 1/(PHI0 * N_c * PHI**2)),
    ("g0^e / (4*Phi_0*N_c)", g0_e / (4*PHI0*N_c)),
    ("g0^e^2/(4*Phi_0) * n_s", g0_e**2/(4*PHI0) * (1-2/56)),
    ("g0^e^2*phi^(-1) / (4*Phi_0)", g0_e**2/(4*PHI0*PHI)),
    ("A_e^2 * (phi-1)", A_e**2 * (PHI-1)),
    ("N_c^3*g0^e/(80*Phi_0*phi)", N_c**3*g0_e/(80*PHI0*PHI)),
    ("1/(4*pi*Phi_0*(phi+1))", 1/(4*np.pi*PHI0*(PHI+1))),
    ("g0^e/(Phi_0^2*phi)", g0_e/(PHI0**2*PHI)),  # unlikely
    ("alpha_s/(2*N_c^2)", ALPHA_S_PDG/(2*N_c**2)),
    ("g0^e^2/(2*Phi_0*phi^2)", g0_e**2/(2*PHI0*PHI**2)),
]

formulas_sorted = sorted(formulas, key=lambda x: abs(x[1] - ALPHA_EM_PDG))

print(f"  alpha_em(PDG) = {ALPHA_EM_PDG:.7f} (1/{1/ALPHA_EM_PDG:.3f})")
print()
print(f"  {'formula':>35s}  {'value':>10s}  {'1/value':>8s}  {'dev%':>8s}")
print("  " + "-" * 70)
for name, val in formulas_sorted[:10]:
    dev = (val - ALPHA_EM_PDG)/ALPHA_EM_PDG * 100
    marker = " <--" if abs(dev) < 2 else " !" if abs(dev) < 5 else ""
    inv = 1/val if val > 0 else 0
    print(f"  {name:>35s}  {val:10.7f}  {inv:8.2f}  {dev:+8.3f}%{marker}")

# ===== 6. RELACJA GUT-LIKE =====
print("\n--- 6. Relacja unifikacyjna (GUT-like) ---\n")

# In SU(5) GUT: sin^2(theta_W) = 3/8 at GUT scale
# alpha_em = alpha_s * sin^2(theta_W) at GUT scale (approximately)

sin2tw = 0.23122  # PDG at M_Z
print(f"  sin^2(theta_W) = {sin2tw} (PDG, MS-bar at M_Z)")
print(f"  alpha_em * sin^2(theta_W) = {ALPHA_EM_PDG * sin2tw:.6f}")
print(f"  alpha_s * sin^2(theta_W) = {ALPHA_S_PDG * sin2tw:.6f}")
print()

# alpha_1 = 5/(3*cos^2(theta_W)) * alpha_em
alpha_1 = 5/(3*(1-sin2tw)) * ALPHA_EM_PDG
alpha_2 = ALPHA_EM_PDG / sin2tw
print(f"  alpha_1(M_Z) = {alpha_1:.6f} (U(1)_Y)")
print(f"  alpha_2(M_Z) = {alpha_2:.6f} (SU(2)_L)")
print(f"  alpha_3(M_Z) = {ALPHA_S_PDG:.6f} (SU(3)_c)")
print()

# TGP: alpha_3 = N_c^3 * g0^e / (8*Phi_0) = 0.1174
# What pattern gives alpha_1, alpha_2?

# alpha_3 = N_c^3 * g0^e / (8*Phi_0)
# alpha_2 = N_c^? * g0^e / (?*Phi_0) = 0.03157?
# alpha_1 = N_c^? * g0^e / (?*Phi_0) = 0.01017?

# Ratios:
print(f"  Ratios:")
print(f"  alpha_3/alpha_2 = {ALPHA_S_PDG/alpha_2:.4f}")
print(f"  alpha_3/alpha_1 = {ALPHA_S_PDG/alpha_1:.4f}")
print(f"  alpha_2/alpha_1 = {alpha_2/alpha_1:.4f}")
print()

# In TGP: the coupling comes from N_c and the soliton.
# For SU(2): N = 2 instead of N_c = 3?
# alpha_2(TGP) = 2^3 * g0^e / (8*Phi_0) = g0^e / Phi_0
alpha_2_TGP = 2**3 * g0_e / (8*PHI0)
alpha_1_TGP = 1**3 * g0_e / (8*PHI0)  # N=1 for U(1)?

print(f"  Hipoteza: alpha_N = N^3 * g0^e / (8*Phi_0)")
print(f"    alpha_3(TGP) = 3^3*g0^e/(8*25) = {N_c**3*g0_e/(8*PHI0):.6f} (PDG: {ALPHA_S_PDG})")
print(f"    alpha_2(TGP) = 2^3*g0^e/(8*25) = {alpha_2_TGP:.6f} (PDG: {alpha_2:.6f})")
print(f"    alpha_1(TGP) = 1^3*g0^e/(8*25) = {alpha_1_TGP:.6f} (PDG: {alpha_1:.6f})")
print()

# Deviations
dev_3 = (N_c**3*g0_e/(8*PHI0) - ALPHA_S_PDG)/ALPHA_S_PDG*100
dev_2 = (alpha_2_TGP - alpha_2)/alpha_2*100
dev_1 = (alpha_1_TGP - alpha_1)/alpha_1*100
print(f"  Odchylenia:")
print(f"    alpha_3: {dev_3:+.2f}%")
print(f"    alpha_2: {dev_2:+.2f}%")
print(f"    alpha_1: {dev_1:+.2f}%")

# ===== 7. OBSERWACJA: 137 = ? =====
print("\n--- 7. Obserwacja: 1/alpha_em = 137.036 ---\n")

# 137 jest Eddingtonowska. Czy ma sens TGP?
# 1/alpha_em ~ 4*Phi_0*phi^3 / g0^e = 4*25*4.236/0.869 = 488.8. No.
# 1/alpha_em ~ Phi_0^2*phi / g0^e = 625*1.618/0.869 = 1163. No.
# 1/alpha_em ~ 4*Phi_0 / g0^e^2 = 100/0.7559 = 132.3. Close-ish!
# 1/alpha_em ~ 4*Phi_0*phi^(1/3) / g0^e^2 = 100*1.175/0.756 = 155.4. No.

inv_alpha = 1/ALPHA_EM_PDG
print(f"  1/alpha_em = {inv_alpha:.3f}")
print()

approx = [
    ("4*Phi_0/g0^e^2", 4*PHI0/g0_e**2),
    ("4*Phi_0*phi^(1/3)/g0^e^2", 4*PHI0*PHI**(1/3)/g0_e**2),
    ("80*Phi_0*phi/(N_c^3*g0^e)", 80*PHI0*PHI/(N_c**3*g0_e)),
    ("2*Phi_0*phi^3/g0^e", 2*PHI0*PHI**3/g0_e),
    ("(2*Phi_0)^2/(N_c*g0^e)", (2*PHI0)**2/(N_c*g0_e)),
    ("8*Phi_0*phi^2/(N_c*g0^e)", 8*PHI0*PHI**2/(N_c*g0_e)),
    ("phi^7 * N_c", PHI**7 * N_c),
    ("Phi_0 * phi^3", PHI0 * PHI**3),
    ("pi^2/alpha_s", np.pi**2/ALPHA_S_PDG),
    ("4*pi^2/alpha_s * g0^e/Phi_0", 4*np.pi**2/ALPHA_S_PDG * g0_e/PHI0),
]

approx_sorted = sorted(approx, key=lambda x: abs(x[1] - inv_alpha))
print(f"  {'formula':>35s}  {'value':>10s}  {'dev%':>8s}")
print("  " + "-" * 60)
for name, val in approx_sorted[:8]:
    dev = (val - inv_alpha)/inv_alpha*100
    marker = " <--" if abs(dev) < 2 else ""
    print(f"  {name:>35s}  {val:10.3f}  {dev:+8.3f}%{marker}")

# ===== WNIOSKI =====
print("\n" + "=" * 72)
print("WNIOSKI ex181")
print("=" * 72)
print(f"""
  1. RELACJA alpha_s/alpha_em:
     alpha_s/alpha_em = {ratio:.3f} ~ 10*phi = {10*PHI:.3f} ({(10*PHI/ratio-1)*100:+.2f}%)
     Intrygujace ale nie potwierdzone teoretycznie.

  2. HIPOTEZY alpha_em:
     H1: g0^e^2/(4*Phi_0) = {alpha_em_H1:.6f} (1/{1/alpha_em_H1:.1f}, PDG: 1/137, {(alpha_em_H1/ALPHA_EM_PDG-1)*100:+.1f}%)
     H2: alpha_s/(10*phi) = {alpha_em_H2:.6f} (1/{1/alpha_em_H2:.1f}, {(alpha_em_H2/ALPHA_EM_PDG-1)*100:+.2f}%)
     Obie odchylone o ~1-4% -- za duzo na predykcje

  3. GUT-LIKE PATTERN:
     alpha_N = N^3 * g0^e / (8*Phi_0):
       alpha_3 (N=3): {N_c**3*g0_e/(8*PHI0):.5f} vs {ALPHA_S_PDG} ({dev_3:+.1f}%)
       alpha_2 (N=2): {alpha_2_TGP:.5f} vs {alpha_2:.5f} ({dev_2:+.1f}%)
       alpha_1 (N=1): {alpha_1_TGP:.5f} vs {alpha_1:.5f} ({dev_1:+.1f}%)
     alpha_2 i alpha_1 NIE daja z tej prostej formuly.
     -> Sektor elektro-slaby wymaga wlasnej derywacji.

  4. STATUS alpha_em:
     TGP NIE predykuje alpha_em z obecnych ram.
     Sektor elektro-slaby (SU(2)_L x U(1)_Y) nie jest jeszcze
     sformalizowany w TGP.
     To jest NOWY program badawczy (nie problem, ale kierunek).

  5. 1/alpha_em = 137.036:
     Najblizsze TGP przyblizenie:
     80*Phi_0*phi/(N_c^3*g0^e) = {80*PHI0*PHI/(N_c**3*g0_e):.2f} ({(80*PHI0*PHI/(N_c**3*g0_e)/inv_alpha-1)*100:+.2f}%)
     (z formuły alpha_em = alpha_s/(10*phi))
     Nie jest predykcja -- raczej relacja fenomenologiczna.
""")
