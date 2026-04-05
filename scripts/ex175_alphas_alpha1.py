#!/usr/bin/env python3
"""
ex175_alphas_alpha1.py
Sesja v45, 2026-04-05

Weryfikacja alpha_s(M_Z) z ODE substratowym (alpha=1).

Pytanie: formula alpha_s = N_c^2 * g0* / (4*Phi_0) uzywa g0*
z warunku B_tail(g0*) = 0. Czy g0* zmienia sie z alpha?

g0* zostal wyznaczony z ODE z f(g) = 1+2*alpha_kin*ln(g).
Teraz uzywamy ODE substratowego (alpha=1):
  g'' + (1/g)(g')^2 + (2/r)g' = 1-g

Porownanie:
- ODE kanoniczne (alpha=2): g'' + (2/g)(g')^2 + (2/r)g' = g^2(1-g)
- ODE substratowe (alpha=1): g'' + (1/g)(g')^2 + (2/r)g' = 1-g
- Stare obliczenia (ex58, ex88): alpha_kin parametryczne

Szukamy g0* (B_tail=0) i z0 (phi-FP) z ODE alpha=1.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
PHI0 = 25.0
PHI0_BRANNEN = 24.783

print("=" * 72)
print("ex175: alpha_s z ODE substratowym (alpha=1)")
print("=" * 72)

# ---- ODE substratowe (alpha=1) ----
def solve_ode(g0, r_max=120, n_points=30000, alpha=1):
    """Solve g'' + (alpha/g)(g')^2 + (2/r)g' = source(g)
    alpha=1: source = 1-g  (substrate)
    alpha=2: source = g^2*(1-g)  (canonical)
    """
    def rhs(r, y):
        g, gp = y
        if g < 1e-12: g = 1e-12
        if alpha == 1:
            src = 1 - g
        else:  # alpha=2 canonical
            src = g**2 * (1 - g)
        gpp = src - (alpha/g)*gp**2 - (2.0/r)*gp if r > 1e-12 else src/3.0
        return [gp, gpp]

    r_eval = np.linspace(1e-6, r_max, n_points)
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0], method='RK45',
                    t_eval=r_eval, rtol=1e-10, atol=1e-12, max_step=0.05)
    return sol.t, sol.y[0]


def extract_tail(g0, r_max=120, alpha=1, window=(50, 100)):
    """Extract A_tail and B_tail (cos and sin coefficients in tail)."""
    r, g = solve_ode(g0, r_max, alpha=alpha)
    delta = g - 1.0
    mask = (r > window[0]) & (r < window[1])
    r_w = r[mask]
    d_w = delta[mask]
    if len(r_w) < 200:
        return None, None

    # Fit: delta(r) ~ [A*cos(r) + B*sin(r)] / r
    # d_w * r_w ~ A*cos(r_w) + B*sin(r_w)
    y = d_w * r_w
    cos_r = np.cos(r_w)
    sin_r = np.sin(r_w)
    # Linear regression
    M = np.column_stack([cos_r, sin_r])
    result = np.linalg.lstsq(M, y, rcond=None)
    A, B = result[0]
    return A, B


def get_Atail(g0, alpha=1):
    A, B = extract_tail(g0, alpha=alpha)
    if A is None:
        return None
    return np.sqrt(A**2 + B**2)  # amplitude


# ---- 1. Find g0* (B_tail = 0) for alpha=1 ----
print("\n--- 1. Punkt staly g0* (B_tail=0) z alpha=1 ---\n")

def B_tail_func(g0, alpha=1):
    A, B = extract_tail(g0, alpha=alpha)
    if B is None:
        return 1e10
    return B

# Scan B_tail as function of g0
print("  Skan B_tail(g0) dla alpha=1:")
g0_scan = np.linspace(1.10, 1.35, 26)
for g0 in g0_scan:
    A, B = extract_tail(g0, alpha=1)
    if B is not None:
        marker = " <-- B~0" if abs(B) < 0.3 else ""
        print(f"    g0 = {g0:.3f}: A = {A:+.4f}, B = {B:+.4f}{marker}")

# Find exact g0* where B=0
try:
    g0_star_a1 = brentq(lambda g0: B_tail_func(g0, alpha=1), 1.20, 1.30, xtol=1e-6)
    A_star, B_star = extract_tail(g0_star_a1, alpha=1)
    print(f"\n  g0*(alpha=1) = {g0_star_a1:.6f}")
    print(f"  A_tail = {A_star:.6f}, B_tail = {B_star:.2e}")
except ValueError:
    print("\n  Nie znaleziono zera B_tail w [1.20, 1.30]")
    # Try wider range
    try:
        g0_star_a1 = brentq(lambda g0: B_tail_func(g0, alpha=1), 1.10, 1.40, xtol=1e-6)
        A_star, B_star = extract_tail(g0_star_a1, alpha=1)
        print(f"\n  g0*(alpha=1) = {g0_star_a1:.6f}")
        print(f"  A_tail = {A_star:.6f}, B_tail = {B_star:.2e}")
    except ValueError:
        g0_star_a1 = None
        print("\n  Brak zera B_tail w [1.10, 1.40]")

# ---- 2. phi-FP point z0 for alpha=1 ----
print("\n--- 2. phi-FP punkt z0 (phi-FP fixed point) z alpha=1 ---\n")

# z0 is defined as g0 where phi-FP gives EXACT r21:
# (A_tail(phi*g0)/A_tail(g0))^4 = r21_lepton
# For ex174: g0^e = 0.86901
g0_e = 0.86901
g0_mu = PHI * g0_e

print(f"  g0^e = {g0_e:.5f} (z phi-FP, ex174)")
print(f"  g0^mu = phi * g0^e = {g0_mu:.5f}")

A_e = get_Atail(g0_e, alpha=1)
A_mu = get_Atail(g0_mu, alpha=1)
if A_e and A_mu:
    r21 = (A_mu/A_e)**4
    print(f"  r21 = {r21:.2f} (PDG: 206.77)")

# ---- 3. alpha_s formula with alpha=1 values ----
print("\n--- 3. alpha_s z ODE alpha=1 ---\n")

N_c = 3

if g0_star_a1 is not None:
    # Formula: alpha_s = N_c^2 * g0* / (4 * Phi_0)
    alpha_s_TGP_25 = N_c**2 * g0_star_a1 / (4 * PHI0)
    alpha_s_TGP_brannen = N_c**2 * g0_star_a1 / (4 * PHI0_BRANNEN)

    print(f"  g0*(alpha=1) = {g0_star_a1:.6f}")
    print(f"  g0*(stare, alpha=2) = 1.24915 (Tw. J2-FP)")
    print()
    print(f"  alpha_s = N_c^2 * g0* / (4*Phi_0):")
    print(f"    Phi_0 = 25:     alpha_s = {alpha_s_TGP_25:.4f}")
    print(f"    Phi_0 = 24.783: alpha_s = {alpha_s_TGP_brannen:.4f}")
    print(f"    PDG:             alpha_s = 0.1179")
    print()
    print(f"  Odchylenie:")
    print(f"    Phi_0=25:     {(alpha_s_TGP_25 - 0.1179)/0.1179*100:+.2f}%")
    print(f"    Phi_0=24.783: {(alpha_s_TGP_brannen - 0.1179)/0.1179*100:+.2f}%")

    # ---- 4. Running scale mu_TGP ----
    print("\n--- 4. Skala mu_TGP (1-loop running) ---\n")

    alpha_s_PDG = 0.1179
    M_Z = 91.1876

    # 1-loop: alpha_s(mu) = alpha_s(M_Z) / [1 + b0*alpha_s(M_Z)/(2pi)*ln(mu^2/M_Z^2)]
    # -> alpha_s^TGP = alpha_s(mu_TGP)
    # Solve for mu_TGP:
    N_f = 5  # active flavors at M_Z
    b0 = (11*N_c - 2*N_f) / 3

    alpha_TGP = N_c**2 * g0_star_a1 / (4 * PHI0_BRANNEN)

    # alpha_TGP = alpha_PDG / [1 + b0*alpha_PDG/(2pi)*ln(mu_TGP^2/M_Z^2)]
    # 1 + b0*alpha_PDG/(2pi)*ln(mu^2/M_Z^2) = alpha_PDG/alpha_TGP
    # ln(mu^2/M_Z^2) = (alpha_PDG/alpha_TGP - 1) * 2*pi / (b0 * alpha_PDG)

    ratio = alpha_s_PDG / alpha_TGP
    ln_mu2 = (ratio - 1) * 2 * np.pi / (b0 * alpha_s_PDG)
    mu_TGP = M_Z * np.exp(ln_mu2 / 2)

    print(f"  alpha_s^TGP = {alpha_TGP:.4f}")
    print(f"  b0 = (11*3 - 2*5)/3 = {b0:.1f}")
    print(f"  mu_TGP = {mu_TGP:.1f} GeV")
    print(f"  mu_TGP / M_Z = {mu_TGP/M_Z:.3f}")

    # Verify: run alpha_s from mu_TGP to M_Z
    alpha_at_MZ = alpha_TGP / (1 + b0*alpha_TGP/(2*np.pi) * np.log(M_Z**2/mu_TGP**2))
    print(f"\n  Weryfikacja: alpha_s(M_Z) z running = {alpha_at_MZ:.4f} (PDG: 0.1179)")

    # ---- 5. Porownanie g0* z alpha=1 i alpha=2 ----
    print("\n--- 5. Porownanie g0* ---\n")
    g0_star_old = 1.24915  # from ex58 with alpha=2 parametrization

    # Note: the old g0* was found with a DIFFERENT ODE!
    # alpha=2 canonical: g'' + (2/g)(g')^2 + (2/r)g' = g^2(1-g)
    # Let's find g0* for alpha=2 too
    try:
        g0_star_a2 = brentq(lambda g0: B_tail_func(g0, alpha=2), 1.10, 1.40, xtol=1e-6)
        print(f"  g0*(alpha=1, substrate) = {g0_star_a1:.6f}")
        print(f"  g0*(alpha=2, canonical) = {g0_star_a2:.6f}")
        print(f"  g0*(stare, ex58)        = {g0_star_old:.5f}")
        print(f"  Roznica alpha1 vs alpha2: {(g0_star_a1-g0_star_a2)/g0_star_a2*100:+.2f}%")
    except:
        print(f"  g0*(alpha=1) = {g0_star_a1:.6f}")
        print(f"  g0*(alpha=2): bariera duchowa - moze nie istniec B_tail=0")

else:
    print("  g0* nie znaleziony - kontynuacja niemozliwa")

# ---- WNIOSKI ----
print("\n" + "=" * 72)
print("WNIOSKI ex175")
print("=" * 72)

if g0_star_a1 is not None:
    print(f"""
  1. g0*(alpha=1) = {g0_star_a1:.5f}
     g0*(alpha=2) = stare wartosci (1.24915)
     Zmiana ODE ZMIENIA g0* -> ZMIENIA alpha_s!

  2. alpha_s^TGP(alpha=1) = N_c^2 * g0* / (4*Phi_0)
     = 9 * {g0_star_a1:.5f} / (4 * 24.783) = {alpha_TGP:.4f}
     PDG: 0.1179
     Odchylenie: {(alpha_TGP - 0.1179)/0.1179*100:+.2f}%

  3. mu_TGP = {mu_TGP:.1f} GeV (skala efektywna TGP)
     1-loop running do M_Z daje alpha_s(M_Z) = {alpha_at_MZ:.4f}

  4. WNIOSEK: przejscie na alpha=1 ZMIENIA predykcje alpha_s!
     Nowa wartosc moze byc blizej lub dalej od PDG.
""")
