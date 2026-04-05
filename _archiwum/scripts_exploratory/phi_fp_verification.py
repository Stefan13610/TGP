#!/usr/bin/env python3
"""
TGP phi-FP (Sciezka 9) -- Weryfikacja numeryczna
==================================================
Implementuje pelny lancuch phi-FP z dodatekJ/J2:
1. ODE solitonu TGP z funkcja kinetyczna f(g) = 1 + 4*ln(g)
2. Ekstrakcja amplitudy ogona oscylacyjnego A_tail
3. Wyznaczenie g_0* z warunku phi-FP: (A(phi*g_0*)/A(g_0*))^4 = r_21^PDG
4. Weryfikacja r_31 z phi^2 i z Koide

Referencja: dodatekJ_ogon_masy.tex, dodatekJ2_sciezka9_formalizacja.tex
TGP v1 -- 2026-03-31
"""

import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.signal import argrelextrema
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# Stale
# ============================================================
PHI = (1 + np.sqrt(5)) / 2        # zlota proporcja = 1.6180339...
R21_PDG = 206.768                  # m_mu / m_e (PDG 2024)
R31_PDG = 3477.48                  # m_tau / m_e (PDG 2024)
ALPHA_TGP = 2                      # alpha z thm:D-uniqueness

# Punkt duchowy (ghost): f(g*) = 0 => 1 + 2*alpha*ln(g*) = 0
# g* = exp(-1/(2*alpha)) = exp(-1/4)
G_GHOST = np.exp(-1 / (2 * ALPHA_TGP))  # ~ 0.7788

print("=" * 70)
print("TGP phi-FP (Sciezka 9) -- Weryfikacja numeryczna")
print("=" * 70)
print(f"  phi = {PHI:.10f}")
print(f"  r_21^PDG = {R21_PDG}")
print(f"  r_31^PDG = {R31_PDG}")
print(f"  alpha_TGP = {ALPHA_TGP}")
print(f"  g* (ghost) = {G_GHOST:.6f}")
print()

# ============================================================
# 1. ODE solitonu TGP z pelna funkcja kinetyczna
# ============================================================
# Rownanie z dodatekJ:
#   f(g) * g'' + (2/r) * g' = V'(g)
# gdzie:
#   f(g) = 1 + 2*alpha*ln(g) = 1 + 4*ln(g)
#   V(g) = g^3/3 - g^4/4  (beta = gamma)
#   V'(g) = g^2 - g^3 = g^2*(1 - g)
#
# Uklad 1-go rzedu: y = [g, g']
#   g' = p
#   p' = [V'(g) - (2/r)*p] / f(g)   -- uwaga na f(g) -> 0

def f_kin(g):
    """Funkcja kinetyczna f(g) = 1 + 4*ln(g)."""
    if g <= 0:
        return -np.inf
    return 1 + 2 * ALPHA_TGP * np.log(g)

def Vprime(g):
    """V'(g) = g^2 - g^3 = g^2*(1 - g)."""
    return g**2 * (1 - g)

def solve_soliton(g0, r_max=300, n_points=30000):
    """
    Rozwiazanie ODE solitonu TGP z f(g)*g'' + (2/r)*g' = V'(g).
    Warunki: g(0) = g_0, g'(0) = 0, g(inf) -> 1.

    Start: rozwiniecie Taylora g = g0 + c2*r^2 + ...
    c2 = V'(g0) / (3*f(g0))  [z warunku regularnosci w r=0]
    """
    fg0 = f_kin(g0)
    if abs(fg0) < 1e-12:
        raise ValueError(f"g0={g0:.6f} jest zbyt blisko punktu duchowego g*={G_GHOST:.6f}")

    c2 = Vprime(g0) / (3 * fg0)

    r_start = 0.01
    g_init = g0 + c2 * r_start**2
    p_init = 2 * c2 * r_start

    def rhs(r, y):
        g, p = y
        if g <= 1e-15:
            return [p, 0.0]
        fg = f_kin(g)
        # Zabezpieczenie: jesli f(g) ~ 0, ograniczamy p''
        if abs(fg) < 1e-10:
            return [p, 0.0]
        if r < 1e-10:
            gpp = Vprime(g) / fg / 3  # limit r->0
        else:
            gpp = (Vprime(g) - (2/r) * p) / fg
        return [p, gpp]

    # Event: stop jesli g diverguje
    def event_diverge(r, y):
        return 100.0 - abs(y[0])
    event_diverge.terminal = True

    r_eval = np.linspace(r_start, r_max, n_points)
    sol = solve_ivp(rhs, [r_start, r_max], [g_init, p_init],
                    t_eval=r_eval, method='RK45', rtol=1e-11, atol=1e-13,
                    max_step=0.05, events=[event_diverge])

    return sol.t, sol.y[0], sol.y[1], sol.status

# ============================================================
# 2. Ekstrakcja amplitudy ogona A_tail
# ============================================================
def extract_Atail(r, g, r_min_fit=100, r_max_fit=250):
    """
    Ekstrakcja A_tail z ogona oscylacyjnego.
    Dla r >> 1: g(r) - 1 ~ (B*cos(r) + C*sin(r)) / r
    wiec (g(r) - 1) * r ~ B*cos(r) + C*sin(r)
    A_tail = sqrt(B^2 + C^2)

    Metoda: fit sinusoidalny do (g-1)*r w przedziale [r_min, r_max].
    """
    mask = (r >= r_min_fit) & (r <= r_max_fit)
    r_fit = r[mask]
    tail = (g[mask] - 1) * r_fit

    if len(r_fit) < 10:
        return np.nan, np.nan, np.nan

    # Fit: tail = B*cos(r) + C*sin(r)
    # Macierz designu
    A_mat = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    result = np.linalg.lstsq(A_mat, tail, rcond=None)
    B, C = result[0]

    A_tail = np.sqrt(B**2 + C**2)
    residual = np.std(tail - A_mat @ np.array([B, C]))

    return A_tail, B, C

# ============================================================
# 3. Skan A_tail(g0)
# ============================================================
print("--- Faza 1: Skan A_tail(g_0) ---")
print()

# Skanujemy g0 w przedziale (g* + eps, 3.5)
g0_values = np.linspace(G_GHOST + 0.05, 3.5, 80)
Atail_values = []
valid_g0 = []

for g0 in g0_values:
    try:
        r, g, p, status = solve_soliton(g0, r_max=300, n_points=30000)
        # Sprawdz czy rozwiazanie doszlo do konca
        if r[-1] < 250:
            Atail_values.append(np.nan)
            valid_g0.append(g0)
            continue
        At, B, C = extract_Atail(r, g, r_min_fit=120, r_max_fit=260)
        Atail_values.append(At)
        valid_g0.append(g0)
    except Exception as e:
        Atail_values.append(np.nan)
        valid_g0.append(g0)

valid_g0 = np.array(valid_g0)
Atail_values = np.array(Atail_values)

# Filtruj NaN
mask_valid = np.isfinite(Atail_values) & (Atail_values > 0)
g0_ok = valid_g0[mask_valid]
At_ok = Atail_values[mask_valid]

print(f"  Rozwiazano {np.sum(mask_valid)}/{len(g0_values)} profilow")
if len(g0_ok) > 0:
    print(f"  g0 range: [{g0_ok[0]:.4f}, {g0_ok[-1]:.4f}]")
    print(f"  A_tail range: [{At_ok.min():.6f}, {At_ok.max():.6f}]")
print()

# ============================================================
# 4. Stosunek mas R(g0) = (A_tail(phi*g0) / A_tail(g0))^4
# ============================================================
print("--- Faza 2: Wyznaczanie phi-FP ---")
print()

def Atail_interpolated(g0_target, g0_arr, At_arr):
    """Interpolacja A_tail dla dowolnego g0."""
    if g0_target < g0_arr[0] or g0_target > g0_arr[-1]:
        return np.nan
    return np.interp(g0_target, g0_arr, At_arr)

def ratio_func(g0):
    """R(g0) = (A_tail(phi*g0) / A_tail(g0))^4 - r21_PDG"""
    A_e = Atail_interpolated(g0, g0_ok, At_ok)
    A_mu = Atail_interpolated(PHI * g0, g0_ok, At_ok)
    if np.isnan(A_e) or np.isnan(A_mu) or A_e < 1e-15:
        return np.nan
    return (A_mu / A_e)**4 - R21_PDG

# Skan ratio
g0_scan = g0_ok[(g0_ok > G_GHOST + 0.1) & (PHI * g0_ok < g0_ok[-1])]
ratio_scan = []
for g0 in g0_scan:
    r = ratio_func(g0)
    ratio_scan.append(r)
ratio_scan = np.array(ratio_scan)

# Szukamy zmiany znaku
sign_changes = []
for i in range(len(ratio_scan) - 1):
    if np.isfinite(ratio_scan[i]) and np.isfinite(ratio_scan[i+1]):
        if ratio_scan[i] * ratio_scan[i+1] < 0:
            sign_changes.append(i)

if len(sign_changes) > 0:
    print(f"  Znaleziono {len(sign_changes)} zmian(e) znaku w R(g0) - r21_PDG")

    for idx in sign_changes:
        g0_lo, g0_hi = g0_scan[idx], g0_scan[idx + 1]
        print(f"  Przedzial: [{g0_lo:.4f}, {g0_hi:.4f}]")

        # Uzyj gęstszego skanu z bezposrednim obliczaniem ODE
        # zamiast interpolacji — dla dokladnosci
        def ratio_direct(g0_val):
            """Bezposrednie obliczenie R(g0) z ODE."""
            try:
                r1, g1, _, s1 = solve_soliton(g0_val, r_max=300, n_points=30000)
                if r1[-1] < 250:
                    return np.nan
                A_e, _, _ = extract_Atail(r1, g1, 120, 260)

                r2, g2, _, s2 = solve_soliton(PHI * g0_val, r_max=300, n_points=30000)
                if r2[-1] < 250:
                    return np.nan
                A_mu, _, _ = extract_Atail(r2, g2, 120, 260)

                if A_e < 1e-15 or np.isnan(A_e) or np.isnan(A_mu):
                    return np.nan
                return (A_mu / A_e)**4 - R21_PDG
            except:
                return np.nan

        # Bisekcja reczna (brentq wymaga ciaglosci)
        g_lo, g_hi = g0_lo, g0_hi
        for _ in range(30):  # ~30 iteracji daje ~10^-9 dokladnosc
            g_mid = (g_lo + g_hi) / 2
            val = ratio_direct(g_mid)
            if np.isnan(val):
                break
            val_lo = ratio_direct(g_lo)
            if np.isnan(val_lo):
                break
            if val * val_lo < 0:
                g_hi = g_mid
            else:
                g_lo = g_mid

        g0_star = (g_lo + g_hi) / 2

        # Weryfikacja koncowa
        r_e, g_e, _, _ = solve_soliton(g0_star, r_max=300, n_points=30000)
        A_e, _, _ = extract_Atail(r_e, g_e, 120, 260)

        g0_mu = PHI * g0_star
        r_mu, g_mu, _, _ = solve_soliton(g0_mu, r_max=300, n_points=30000)
        A_mu, _, _ = extract_Atail(r_mu, g_mu, 120, 260)

        r21_calc = (A_mu / A_e)**4

        print()
        print("  === WYNIK phi-FP ===")
        print(f"  g_0* (elektron)  = {g0_star:.6f}")
        print(f"  g_0^mu (mion)    = {g0_mu:.6f}  [= phi * g_0*]")
        print(f"  A_tail(e)        = {A_e:.6f}")
        print(f"  A_tail(mu)       = {A_mu:.6f}")
        print(f"  r_21 = (A_mu/A_e)^4 = {r21_calc:.4f}")
        print(f"  r_21^PDG             = {R21_PDG}")
        print(f"  Odchylenie           = {abs(r21_calc - R21_PDG)/R21_PDG*100:.4f}%")
        print()

        # --- Tau: phi^2 * g0* ---
        g0_tau_phi2 = PHI**2 * g0_star
        print(f"  --- Tau (naiwne phi^2) ---")
        print(f"  g_0^tau = phi^2 * g_0* = {g0_tau_phi2:.6f}")

        if g0_tau_phi2 < g0_ok[-1]:
            try:
                r_tau, g_tau, _, _ = solve_soliton(g0_tau_phi2, r_max=300, n_points=30000)
                if r_tau[-1] >= 250:
                    A_tau, _, _ = extract_Atail(r_tau, g_tau, 120, 260)
                    r31_phi2 = (A_tau / A_e)**4
                    print(f"  A_tail(tau)      = {A_tau:.6f}")
                    print(f"  r_31 (phi^2)     = {r31_phi2:.2f}")
                    print(f"  r_31^PDG         = {R31_PDG}")
                    print(f"  Odchylenie       = {abs(r31_phi2 - R31_PDG)/R31_PDG*100:.2f}%")
                else:
                    print(f"  [ODE diverguje przy r={r_tau[-1]:.1f}]")
            except Exception as e:
                print(f"  [Blad: {e}]")
        else:
            print(f"  [g0_tau poza zakresem skanu]")

        # --- Tau: z Koide Q=3/2 ---
        print()
        print(f"  --- Tau (Koide Q = 3/2) ---")
        # Analitycznie: 2*(1 + sqrt(r21) + sqrt(r31))^2 = 3*(1 + r21 + r31)
        # Rozwiazanie: a = 1 + sqrt(r21), r31 = (1/4)*(2a + sqrt(6a^2 - 3 - 3*r21))^2
        a = 1 + np.sqrt(r21_calc)
        discriminant = 6 * a**2 - 3 - 3 * r21_calc
        if discriminant > 0:
            r31_koide = 0.25 * (2*a + np.sqrt(discriminant))**2
            print(f"  r_31 (Koide)     = {r31_koide:.2f}")
            print(f"  r_31^PDG         = {R31_PDG}")
            print(f"  Odchylenie       = {abs(r31_koide - R31_PDG)/R31_PDG*100:.4f}%")

            # Masa tau z Koide
            m_e = 0.51099895  # MeV
            m_tau_koide = m_e * r31_koide
            print(f"  m_tau (Koide)    = {m_tau_koide:.2f} MeV")
            print(f"  m_tau (PDG)      = 1776.86 MeV")
            print(f"  Odchylenie       = {abs(m_tau_koide - 1776.86)/1776.86*100:.4f}%")

            # Q Koidego
            Q = (1 + np.sqrt(r21_calc) + np.sqrt(r31_koide))**2 / (1 + r21_calc + r31_koide)
            print(f"  Q (Koide)        = {Q:.8f}")
        print()

        # --- Tau: harmoniczna 2*g0* ---
        g0_tau_harm = 2 * g0_star
        print(f"  --- Tau (harmoniczna 2*g_0*) ---")
        print(f"  g_0^tau = 2 * g_0* = {g0_tau_harm:.6f}")

        if g0_tau_harm < g0_ok[-1]:
            try:
                r_tau2, g_tau2, _, _ = solve_soliton(g0_tau_harm, r_max=300, n_points=30000)
                if r_tau2[-1] >= 250:
                    A_tau2, _, _ = extract_Atail(r_tau2, g_tau2, 120, 260)
                    r31_harm = (A_tau2 / A_e)**4
                    print(f"  A_tail(tau)      = {A_tau2:.6f}")
                    print(f"  r_31 (harm.)     = {r31_harm:.2f}")
                    print(f"  Odchylenie       = {abs(r31_harm - R31_PDG)/R31_PDG*100:.2f}%")
                else:
                    print(f"  [ODE diverguje przy r={r_tau2[-1]:.1f}]")
            except Exception as e:
                print(f"  [Blad: {e}]")
        print()

else:
    print("  BRAK zmiany znaku -- phi-FP nie znaleziony w tym przedziale")
    print("  Skan R(g0):")
    for i in range(0, len(g0_scan), max(1, len(g0_scan)//20)):
        if np.isfinite(ratio_scan[i]):
            r_val = ratio_scan[i] + R21_PDG
            print(f"    g0={g0_scan[i]:.4f}: R={(r_val):.2f}")
    print()
    print("  Profilowanie: wyswietlam A_tail(g0) dla diagnostyki")

# ============================================================
# 5. Wykresy
# ============================================================
script_dir = os.path.dirname(os.path.abspath(__file__))

# Wykres 1: A_tail(g0)
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

ax1 = axes[0]
ax1.semilogy(g0_ok, At_ok, 'b-', linewidth=1.5)
ax1.axvline(G_GHOST, color='red', linestyle='--', alpha=0.7, label=f'g* = {G_GHOST:.4f} (ghost)')
ax1.set_xlabel(r'$g_0$', fontsize=13)
ax1.set_ylabel(r'$A_{\mathrm{tail}}(g_0)$', fontsize=13)
ax1.set_title('Amplituda ogona solitonu TGP')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Wykres 2: R(g0) = (A(phi*g0)/A(g0))^4
ax2 = axes[1]
valid_ratio = np.isfinite(ratio_scan)
if np.any(valid_ratio):
    ax2.semilogy(g0_scan[valid_ratio], ratio_scan[valid_ratio] + R21_PDG, 'g-', linewidth=1.5)
    ax2.axhline(R21_PDG, color='red', linestyle='--', alpha=0.7, label=f'r_21^PDG = {R21_PDG}')
    ax2.set_xlabel(r'$g_0$', fontsize=13)
    ax2.set_ylabel(r'$R(g_0) = (A_{\mathrm{tail}}(\varphi g_0)/A_{\mathrm{tail}}(g_0))^4$', fontsize=13)
    ax2.set_title('Stosunek mas phi-FP')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(os.path.join(script_dir, 'phi_fp_atail.png'), dpi=150, bbox_inches='tight')
print(f"  Wykres zapisany: scripts/phi_fp_atail.png")

# Wykres 3: Profile solitonow (e, mu, tau)
fig2, ax3 = plt.subplots(figsize=(10, 6))
try:
    # Ponownie oblicz najlepsze profile
    if len(sign_changes) > 0:
        g0_best = g0_star
    else:
        g0_best = 1.25  # fallback

    profiles = {}
    for label, g0_val in [('e', g0_best), ('mu', PHI * g0_best), ('tau (phi^2)', PHI**2 * g0_best)]:
        try:
            r_p, g_p, _, _ = solve_soliton(g0_val, r_max=200, n_points=20000)
            profiles[label] = (r_p, g_p, g0_val)
            ax3.plot(r_p, g_p, linewidth=1.5, label=f'{label}: $g_0$={g0_val:.4f}')
        except:
            pass

    ax3.axhline(1.0, color='gray', linestyle=':', alpha=0.5, label='vacuum g=1')
    ax3.axhline(G_GHOST, color='red', linestyle='--', alpha=0.4, label=f'ghost g*={G_GHOST:.4f}')
    ax3.set_xlabel(r'$r$', fontsize=13)
    ax3.set_ylabel(r'$g(r) = \Phi(r)/\Phi_0$', fontsize=13)
    ax3.set_title('Profile solitonowe TGP (e, mu, tau)')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, 100)
except Exception as e:
    print(f"  [Blad profili: {e}]")

fig2.savefig(os.path.join(script_dir, 'phi_fp_profiles.png'), dpi=150, bbox_inches='tight')
print(f"  Wykres zapisany: scripts/phi_fp_profiles.png")

# ============================================================
# 6. Podsumowanie
# ============================================================
print()
print("=" * 70)
print("PODSUMOWANIE phi-FP")
print("=" * 70)
print(f"  ODE: f(g)*g'' + (2/r)*g' = V'(g)")
print(f"  f(g) = 1 + 4*ln(g),  V'(g) = g^2*(1-g)")
print(f"  Ghost: g* = exp(-1/4) = {G_GHOST:.6f}")
print(f"  phi = (1+sqrt(5))/2  = {PHI:.10f}")
print()
if len(sign_changes) > 0:
    print(f"  phi-FP ZNALEZIONY:")
    print(f"    g_0* = {g0_star:.6f}")
    print(f"    r_21 = {r21_calc:.4f}  (PDG: {R21_PDG})")
    if discriminant > 0:
        print(f"    r_31 (Koide) = {r31_koide:.2f}  (PDG: {R31_PDG})")
else:
    print(f"  phi-FP NIE ZNALEZIONY w przedziale [{g0_values[0]:.2f}, {g0_values[-1]:.2f}]")
    print(f"  Diagnostyka: sprawdz A_tail(g0) na wykresie phi_fp_atail.png")

print()
print("GOTOWE.")
