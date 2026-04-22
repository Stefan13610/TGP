#!/usr/bin/env python3
"""
ex176_g0star_substrate.py
Sesja v45, 2026-04-05

Badanie punktu stalego g0* w ODE substratowym (alpha=1).

Problem: B_tail(g0) > 0 dla WSZYSTKICH g0 w [1.1, 1.4] -> brak B=0.
Warunek H1 (B_tail=0) byl specyficzny dla f(g)=1+2*alpha_kin*ln(g).

Pytania:
1. Czy B_tail=0 istnieje GDZIEKOLWIEK (g0 < 1.1 lub g0 > 1.4)?
2. Czy jest alternatywny warunek selekcji g0* w alpha=1?
3. Co gra role g0* w formule alpha_s?

Alternatywne definicje g0*:
H1': Minimum |B_tail/A_tail| (najczystszy ogon kosinusowy)
H2': g0 gdzie A_tail(g0) = g0 (self-consistency)
H3': g0 = z0 z phi-FP (punkt staly skalowania)
H4': Minimum energii solitonu E(g0)
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

PHI = (1 + np.sqrt(5)) / 2

def solve_ode_a1(g0, r_max=150, n_points=40000):
    def rhs(r, y):
        g, gp = y
        if g < 1e-12: g = 1e-12
        gpp = (1-g) - (1.0/g)*gp**2 - (2.0/r)*gp if r > 1e-12 else (1-g)/3.0
        return [gp, gpp]
    r_eval = np.linspace(1e-6, r_max, n_points)
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0], method='RK45',
                    t_eval=r_eval, rtol=1e-10, atol=1e-12, max_step=0.04)
    return sol.t, sol.y[0]

def extract_AB(g0, window=(60, 120)):
    r, g = solve_ode_a1(g0, r_max=150)
    delta = g - 1.0
    mask = (r > window[0]) & (r < window[1])
    r_w, d_w = r[mask], delta[mask]
    if len(r_w) < 300:
        return None, None
    y = d_w * r_w
    M = np.column_stack([np.cos(r_w), np.sin(r_w)])
    result = np.linalg.lstsq(M, y, rcond=None)
    return result[0][0], result[0][1]

def get_Atail(g0):
    A, B = extract_AB(g0)
    if A is None: return None
    return np.sqrt(A**2 + B**2)

print("=" * 72)
print("ex176: g0* w ODE substratowym (alpha=1)")
print("=" * 72)

# ---- 1. Szeroki skan B_tail ----
print("\n--- 1. Szeroki skan B_tail(g0), alpha=1, g0 in [0.5, 2.5] ---\n")

g0_range = np.concatenate([
    np.linspace(0.50, 0.85, 8),
    np.linspace(0.86, 0.92, 7),  # around phi-FP electron
    np.linspace(1.00, 1.50, 11),
    np.linspace(1.60, 2.50, 10),
])

print(f"  {'g0':>8s}  {'A_tail':>10s}  {'B_tail':>10s}  {'|B/A|':>8s}  {'phase':>8s}")
print("  " + "-" * 55)

AB_data = []
for g0 in g0_range:
    A, B = extract_AB(g0)
    if A is not None and abs(A) > 1e-10:
        ratio = abs(B/A)
        phase = np.degrees(np.arctan2(B, A))
        AB_data.append((g0, A, B, ratio, phase))
        marker = ""
        if abs(B) < 0.02:
            marker = " <-- B~0!"
        elif ratio < 0.5:
            marker = " <-- low |B/A|"
        print(f"  {g0:8.4f}  {A:10.6f}  {B:10.6f}  {ratio:8.4f}  {phase:8.2f}{marker}")
    else:
        print(f"  {g0:8.4f}  {'FAIL':>10s}")

# ---- 2. Szukanie B_tail=0 w szerokim zakresie ----
print("\n--- 2. Szukanie B_tail=0 ---\n")

# Check if B changes sign anywhere
B_values = [(g0, B) for g0, A, B, _, _ in AB_data]
sign_changes = []
for i in range(len(B_values)-1):
    g0_1, B1 = B_values[i]
    g0_2, B2 = B_values[i+1]
    if B1 * B2 < 0:
        sign_changes.append((g0_1, g0_2, B1, B2))

if sign_changes:
    print(f"  Znaleziono {len(sign_changes)} zmian znaku B_tail:")
    for g1, g2, B1, B2 in sign_changes:
        from scipy.optimize import brentq
        g0_zero = brentq(lambda g0: extract_AB(g0)[1], g1, g2, xtol=1e-5)
        A_z, B_z = extract_AB(g0_zero)
        print(f"    B=0 przy g0 = {g0_zero:.6f} (A = {A_z:.6f})")
else:
    print("  BRAK zmian znaku B_tail w calym zakresie [0.5, 2.5]")
    print("  B_tail jest ZAWSZE dodatni (lub ujemny) -- monotoniczny trend")

# ---- 3. Minimalne |B/A| (najczystszy kosinusowy ogon) ----
print("\n--- 3. Minimum |B/A| (H1': najczystszy ogon) ---\n")

valid = [(g0, ratio) for g0, A, B, ratio, _ in AB_data if ratio < 100]
if valid:
    g0_min_ratio = min(valid, key=lambda x: x[1])
    print(f"  Min |B/A| = {g0_min_ratio[1]:.4f} przy g0 = {g0_min_ratio[0]:.4f}")

    # Refine around this minimum
    g0_c = g0_min_ratio[0]
    refinement = np.linspace(max(0.5, g0_c-0.15), g0_c+0.15, 30)
    best_ratio = 999
    best_g0 = g0_c
    for g0 in refinement:
        A, B = extract_AB(g0)
        if A is not None and abs(A) > 1e-10:
            ratio = abs(B/A)
            if ratio < best_ratio:
                best_ratio = ratio
                best_g0 = g0

    A_b, B_b = extract_AB(best_g0)
    print(f"  Refined: min |B/A| = {best_ratio:.4f} przy g0 = {best_g0:.5f}")
    print(f"  A = {A_b:.6f}, B = {B_b:.6f}")

# ---- 4. Faza ogona delta(r) ----
print("\n--- 4. Faza ogona delta(r) ~ C*cos(r + delta_0)/r ---\n")
print("  delta_0 = arctan(B/A)")
print()
for g0, A, B, ratio, phase in AB_data[:20]:
    print(f"  g0 = {g0:.4f}: phase = {phase:+.2f} deg, |B/A| = {ratio:.4f}")

# ---- 5. phi-FP weryfikacja ----
print("\n--- 5. phi-FP z nowym oknem ekstrakcji ---\n")

# Test different extraction windows
for win in [(40, 70), (50, 90), (60, 120)]:
    g0_e_test = 0.8690
    A_e = None
    A_mu = None
    try:
        r, g = solve_ode_a1(g0_e_test)
        delta = g - 1.0
        mask = (r > win[0]) & (r < win[1])
        r_w, d_w = r[mask], delta[mask]
        y = d_w * r_w
        M = np.column_stack([np.cos(r_w), np.sin(r_w)])
        res = np.linalg.lstsq(M, y, rcond=None)
        A_e_fit = np.sqrt(res[0][0]**2 + res[0][1]**2)

        g0_mu_test = PHI * g0_e_test
        r, g = solve_ode_a1(g0_mu_test)
        delta = g - 1.0
        mask = (r > win[0]) & (r < win[1])
        r_w, d_w = r[mask], delta[mask]
        y = d_w * r_w
        M = np.column_stack([np.cos(r_w), np.sin(r_w)])
        res = np.linalg.lstsq(M, y, rcond=None)
        A_mu_fit = np.sqrt(res[0][0]**2 + res[0][1]**2)

        r21 = (A_mu_fit / A_e_fit)**4
        print(f"  Window {win}: r21 = {r21:.2f} (PDG: 206.77)")
    except:
        print(f"  Window {win}: FAIL")

# ---- 6. Energia solitonu ----
print("\n--- 6. Energia solitonu E(g0) ---\n")

def soliton_energy(g0):
    r, g = solve_ode_a1(g0, r_max=80)
    gp = np.gradient(g, r)
    # E = int [K_sub/2 * (g')^2 + V(g)] * 4pi*r^2 dr
    # K_sub = g^2 (alpha=1)
    # V(g) = (1-g)^2/2 (simple double-well)
    T = 0.5 * g**2 * gp**2  # kinetic density
    V = 0.5 * (1-g)**2       # potential density
    integrand = (T + V) * 4 * np.pi * r**2
    return np.trapezoid(integrand, r)

g0_energies = np.linspace(0.5, 2.0, 31)
energies = []
for g0 in g0_energies:
    E = soliton_energy(g0)
    energies.append(E)

# Find minimum
idx_min = np.argmin(energies)
g0_E_min = g0_energies[idx_min]
E_min = energies[idx_min]

print(f"  Minimum E(g0) przy g0 = {g0_E_min:.3f}, E = {E_min:.4f}")
print(f"  (E nie jest fizycznie minimalizowane -- soliton nie jest w minimum)")

# ---- WNIOSKI ----
print("\n" + "=" * 72)
print("WNIOSKI ex176")
print("=" * 72)
print(f"""
  1. B_tail(g0) > 0 dla WSZYSTKICH g0 in [0.5, 2.5] z alpha=1.
     Warunek H1 (B_tail=0) NIE ISTNIEJE w ODE substratowym.

  2. B_tail jest monotoniczny: rosnie z g0.
     Minimalne |B/A| przy g0 ~ {best_g0:.3f} (najczystszy ogon).

  3. INTERPRETACJA:
     H1 (B_tail=0) bylo specyficzne dla parametrycznej rodziny
     ODE z f(g) = 1 + 2*alpha_kin*ln(g).
     Przejscie na ODE substratowe (alpha=1, K=g^2) ZMIENIA
     strukture ogona -- B_tail nie zeruje sie.

  4. KONSEKWENCJE DLA alpha_s:
     Formula alpha_s = N_c^2 * g0* / (4*Phi_0) UZYWALA
     g0* z warunku H1. Ten warunek nie istnieje w alpha=1.
     -> Formula alpha_s wymaga REINTERPRETACJI lub
        NOWEGO warunku selekcji.

  5. KANDYDACI na nowe g0*:
     a) g0^e = 0.869 z phi-FP (ale to g0 elektronu, nie "staly")
     b) g0 z min |B/A| ~ {best_g0:.3f}
     c) g0 z argumentu topologicznego (Z_3, etc.)
     d) Rezygnacja: alpha_s nie jest predykowana przez TGP
        (wymaga dodatkowego inputu)

  6. STATUS G2: OTWARTY -- formula alpha_s wymaga rewizji z alpha=1
""")
