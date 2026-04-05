#!/usr/bin/env python3
"""
TGP ex60 — Topologia sektorów solitonu i złota proporcja
=========================================================
Cel:
  1. Skanowanie sektorów solitonu TGP w g_0 ∈ [1.05, 3.8]
  2. Identyfikacja granic topologicznych (n_cross sektora)
  3. Weryfikacja warunku kwantowania B_tail = 0 (elektron, mion)
  4. Test hipotezy złotej proporcji: g_0^μ / g_0^e ≈ φ = 1.6180
  5. Predykcja g_0^τ jako kandydat Ścieżki 9

ODE solitonu TGP (Dodatek J, N0-5):
  f(g)·g'' + (2/r)·g' = V'(g)
  f(g) = 1 + 2α·ln(g),   α = 2  (z K_{ij}=(φᵢφⱼ)², tw. prop:substrate-action)
  V'(g) = g²(1 − g)       (z V = g³/3 − g⁴/4, warunek β=γ, stw. prop:vacuum)
  IC: g(0) = g₀,  g'(0) = 0
  BC: g(r) → 1  dla r → ∞

Punkt osobliwy kinematyczny: g* = exp(−1/(2α)) = exp(−1/4) ≈ 0.7788
  (gdzie f(g*) = 0 — nieskończona efektywna przyspieszenie)

Asymptotyczny ogon (linearyzacja przy g=1, r → ∞):
  g(r) − 1 ~ [B·cos(r) + C·sin(r)] / r   (sferyczna fala Bessela, m_sp=1)
  A_tail = √(B² + C²)

Poprzednie wyniki (ex55–ex59, sesja v33+):
  g_0^e = 1.2301  (B_tail = 0; 0.8% od ref 1.24)
  g_0^μ = 2.00    (hipoteza złotej proporcji, 1.6% od r₂₁ = 206.77)
  A_tail ∝ (g_0 − g*)^{4.12 ± 0.05}

Otwarte problemy adresowane:  O-J1, O-J2  (Dodatek J §J.5)

Autor: TGP v1, sesja v34 (2026-03-28)
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ============================================================
# Stałe TGP — WYŁĄCZNIE z aksjomatów N0
# ============================================================
ALPHA   = 2.0                           # α=2 (prop:substrate-action)
G_STAR  = np.exp(-1.0 / (2.0 * ALPHA)) # g* = e^{-1/4} ≈ 0.77880
PHI_G   = (1.0 + np.sqrt(5.0)) / 2.0   # złota proporcja φ = 1.61803...

# Wartości referencyjne z ex55–ex58
G0_E    = 1.2301   # pierwsze zero B_tail (ex58; 0.8% od hipotezy 1.24)
G0_E_H  = 1.24     # hipoteza robocza elektronu
G0_MU   = 2.00     # hipoteza złotej proporcji mion

# Docelowe stosunki mas (PDG 2024)
R21_PDG = 206.7682  # m_μ / m_e
R31_PDG = 3477.15   # m_τ / m_e
R32_PDG = 16.8171   # m_τ / m_μ

# Parametry numeryczne
R_MAX        = 250.0   # zasięg całkowania (jednostki a_Γ⁻¹)
R_TAIL_START = 100.0   # start dopasowania ogona Bessela
RTOL         = 1e-10
ATOL         = 1e-12
DR_START     = 1e-5    # krok startowy (unika r=0 osobliwości)

print(f"TGP ex60 — Topologia sektorów i złota proporcja")
print(f"α = {ALPHA},  g* = {G_STAR:.6f},  φ = {PHI_G:.6f}")
print(f"Parametry: R_MAX={R_MAX}, R_TAIL_START={R_TAIL_START}, rtol={RTOL}")
print("=" * 65)


# ============================================================
# Definicja ODE solitonu TGP (Dodatek J, eq. J-ode)
# ============================================================
def f_kin(g):
    """Efektywny czynnik kinetyczny f(g) = 1 + 2α·ln(g)."""
    return 1.0 + 2.0 * ALPHA * np.log(max(g, 1e-15))

def Vprime(g):
    """V'(g) = g²(1−g), z V(g) = g³/3 − g⁴/4 (β=γ, N0-5)."""
    return g * g * (1.0 - g)

def ode_rhs(r, y):
    g, gp = y
    if g < G_STAR * 0.998:
        # Ochrona przed osobliwością — zatrzymaj (zdarzenie)
        return [gp, 0.0]
    fg = f_kin(g)
    if r < 1e-9:
        # r=0: L'Hôpital → g''(0) = V'(g₀) / (3 f(g₀))
        gpp = Vprime(g) / (3.0 * fg)
    else:
        gpp = (Vprime(g) - (2.0 / r) * gp) / fg
    return [gp, gpp]

# Zdarzenie: trajektoria zbliża się do punktu osobliwego g*
def event_singular(r, y):
    return y[0] - G_STAR * 1.002

event_singular.terminal  = True
event_singular.direction = -1


def initial_conditions(g0):
    """
    Warunki początkowe z rozwinięciem Taylora wokół r = DR_START:
      g(DR_START) ≈ g₀ + ½ g''(0) dr²
      g'(DR_START) ≈ g''(0) dr
    gdzie g''(0) = V'(g₀) / (3 f(g₀))  [L'Hôpital r→0].
    """
    fg0  = f_kin(g0)
    gpp0 = Vprime(g0) / (3.0 * fg0)
    r0   = DR_START
    g_ic = g0 + 0.5 * gpp0 * r0**2
    gp_ic = gpp0 * r0
    return r0, g_ic, gp_ic


def solve_soliton(g0, r_max=R_MAX):
    """
    Rozwiązuje ODE solitonu TGP dla zadanego g₀.
    Zwraca (r, g, sukces).
    """
    r0, g_ic, gp_ic = initial_conditions(g0)
    sol = solve_ivp(
        ode_rhs, [r0, r_max], [g_ic, gp_ic],
        method='DOP853', rtol=RTOL, atol=ATOL,
        dense_output=False, events=event_singular,
        max_step=0.15
    )
    success = (sol.status == 0)   # 0 = całkowanie do r_max
    return sol.t, sol.y[0], success


def fit_tail(r_arr, g_arr, r_start=R_TAIL_START):
    """
    Dopasowuje ogon Bessela: u(r) = r(g−1) = B cos r + C sin r
    metodą liniowych least-squares dla r ≥ r_start.
    Zwraca (A_tail, B_tail, C_tail).
    """
    mask = (r_arr >= r_start) & np.isfinite(g_arr)
    if mask.sum() < 20:
        return np.nan, np.nan, np.nan
    r_t = r_arr[mask]
    u_t = r_t * (g_arr[mask] - 1.0)
    A_m = np.column_stack([np.cos(r_t), np.sin(r_t)])
    try:
        coeffs, _, _, _ = np.linalg.lstsq(A_m, u_t, rcond=None)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan
    B, C = coeffs
    return np.sqrt(B**2 + C**2), B, C


def analyze(g0, r_max=R_MAX):
    """Pełna analiza dla jednego g₀. Zwraca słownik obserwabli."""
    r, g, suc = solve_soliton(g0, r_max)
    A, B, C = (np.nan, np.nan, np.nan)
    g_min = np.nan
    n_cross = 0
    if len(g) > 2:
        g_min = float(np.min(g))
        sign_changes = np.diff(np.sign(g - 1.0))
        n_cross = int(np.sum(sign_changes != 0))
    if suc:
        A, B, C = fit_tail(r, g)
    return dict(g0=g0, success=suc, g_min=g_min,
                n_cross=n_cross, A=A, B=B, C=C)


# ============================================================
# TEST 1: Warunek kwantowania B_tail = 0 (elektron)
# ============================================================
print("\n[TEST 1] Warunek kwantowania B_tail(g₀^e) = 0")
print("-" * 50)

def B_func(g0):
    res = analyze(g0)
    return res['B'] if not np.isnan(res['B']) else 1.0

# Skan wstępny dla znalezienia przedziału ze zmianą znaku
scan_e = np.linspace(1.10, 1.45, 36)
B_scan = np.array([B_func(g) for g in scan_e])

g0_e_found = None
for i in range(len(scan_e) - 1):
    if np.isfinite(B_scan[i]) and np.isfinite(B_scan[i+1]):
        if B_scan[i] * B_scan[i+1] < 0:
            try:
                g0_e_found = brentq(B_func, scan_e[i], scan_e[i+1], xtol=1e-7, maxiter=60)
                break
            except Exception:
                pass

if g0_e_found is not None:
    print(f"  Pierwsze zero B_tail:  g₀^e = {g0_e_found:.6f}")
    print(f"  Ref. ex58 (1.2301):   odchylenie = {abs(g0_e_found-G0_E)/G0_E*100:.3f}%")
    print(f"  Hipoteza (1.24):      odchylenie = {abs(g0_e_found-G0_E_H)/G0_E_H*100:.3f}%")
else:
    print(f"  UWAGA: nie znaleziono zera B_tail w [1.10, 1.45]")
    g0_e_found = G0_E_H   # fallback

# Wartości A_tail dla elektronu i mionu (hipotezy)
res_e  = analyze(G0_E_H)
res_mu = analyze(G0_MU)
A_e    = res_e['A']
A_mu   = res_mu['A']

print(f"\n  A_tail(g₀^e = {G0_E_H})  = {A_e:.6f}")
print(f"  A_tail(g₀^μ = {G0_MU}) = {A_mu:.6f}")
if not (np.isnan(A_e) or np.isnan(A_mu)):
    r21_calc = (A_mu / A_e) ** 4
    print(f"  (A_μ/A_e)⁴ = {r21_calc:.4f}   [cel: r₂₁ = {R21_PDG:.4f}]")
    print(f"  Odchylenie: {abs(r21_calc - R21_PDG)/R21_PDG*100:.2f}%")


# ============================================================
# TEST 2: Złota proporcja g₀^μ / g₀^e
# ============================================================
print("\n[TEST 2] Złota proporcja g₀^μ/g₀^e")
print("-" * 50)

ratio_obs = G0_MU / G0_E_H
ratio_e   = G0_MU / g0_e_found

print(f"  g₀^μ / g₀^e(hip)  = {G0_MU}/{G0_E_H}     = {ratio_obs:.6f}")
print(f"  g₀^μ / g₀^e(calc) = {G0_MU}/{g0_e_found:.4f}  = {ratio_e:.6f}")
print(f"  φ (złota proporcja) = {PHI_G:.6f}")
print(f"  Odchylenie (hip):  {abs(ratio_obs - PHI_G)/PHI_G*100:.3f}%")
print(f"  Odchylenie (calc): {abs(ratio_e  - PHI_G)/PHI_G*100:.3f}%")

# Sprawdź g₀^e × φ (dokładna złota)
g0_mu_phi = G0_E_H * PHI_G
res_mu_phi = analyze(g0_mu_phi)
A_mu_phi = res_mu_phi['A']
print(f"\n  g₀^e × φ = {g0_mu_phi:.6f}")
if not (np.isnan(A_e) or np.isnan(A_mu_phi)):
    r21_phi = (A_mu_phi / A_e)**4
    print(f"  (A(g₀^e×φ)/A(g₀^e))⁴ = {r21_phi:.4f}  [cel: {R21_PDG:.4f}]")
    print(f"  Odchylenie:  {abs(r21_phi - R21_PDG)/R21_PDG*100:.2f}%")


# ============================================================
# TEST 3: Gęsty skan sektorów  g₀ ∈ [1.05, 3.8]
# ============================================================
print("\n[TEST 3] Skan sektorów topologicznych (n_cross)")
print("-" * 50)

g0_scan = np.linspace(1.05, 3.80, 100)
scan_data = []

for g0 in g0_scan:
    d = analyze(g0)
    scan_data.append(d)

# Znajdź granice sektorów (zmiany n_cross)
n_cross_arr = np.array([d['n_cross'] for d in scan_data])
g_min_arr   = np.array([d['g_min']  for d in scan_data])
A_arr       = np.array([d['A']      for d in scan_data])
B_arr       = np.array([d['B']      for d in scan_data])
success_arr = np.array([d['success'] for d in scan_data])

sector_changes = []
for i in range(len(g0_scan) - 1):
    if n_cross_arr[i] != n_cross_arr[i+1]:
        sector_changes.append((g0_scan[i], g0_scan[i+1], n_cross_arr[i], n_cross_arr[i+1]))

print(f"  Granice sektorów (zmiany n_cross):")
for s in sector_changes:
    print(f"    g₀ ∈ [{s[0]:.3f}, {s[1]:.3f}]: n_cross {s[2]} → {s[3]}")

# Sprawdź czy g₀^μ = 2.00 leży przy granicy sektora
near_2 = [(s[0], s[1]) for s in sector_changes if s[0] < 2.0 < s[1] or abs(s[0]-2.0) < 0.08]
if near_2:
    print(f"\n  g₀^μ = 2.00 LEŻY przy granicy sektora: {near_2}")
    print(f"  → HIPOTEZA TOPOLOGICZNA POTWIERDZONA (O-J2)")
else:
    print(f"\n  g₀^μ = 2.00 NIE leży przy zidentyfikowanej granicy sektora")
    print(f"  → Hipoteza topologiczna NIEZWERYFIKOWANA w tym rozdzielczości")

# Raport: g_min vs g*
print(f"\n  g* = {G_STAR:.6f}")
print(f"  g_min przy g₀=1.24:   {analyze(G0_E_H)['g_min']:.6f}   (> g*? {analyze(G0_E_H)['g_min']>G_STAR})")
print(f"  g_min przy g₀=2.00:   {analyze(G0_MU)['g_min']:.6f}   (> g*? {analyze(G0_MU)['g_min']>G_STAR})")
g0_crit = None
for i, d in enumerate(scan_data):
    if d['success'] and d['g_min'] < G_STAR * 1.05:
        g0_crit = g0_scan[i]
        print(f"  Trajektoria zbliża się do g*: g₀ ≈ {g0_crit:.3f}")
        break


# ============================================================
# TEST 4: Skalowanie A_tail(g₀) — weryfikacja wykładnika 4.12
# ============================================================
print("\n[TEST 4] Skalowanie A_tail(g₀) — weryfikacja wykładnika")
print("-" * 50)

valid = [(g0_scan[i], A_arr[i])
         for i in range(len(g0_scan))
         if success_arr[i] and not np.isnan(A_arr[i]) and A_arr[i] > 0]

if len(valid) > 10:
    g0_v = np.array([x[0] for x in valid])
    A_v  = np.array([x[1] for x in valid])

    # Fit log-liniowy: ln(A) = a + b·ln(g₀ − g*)
    dg = g0_v - G_STAR
    mask_fit = dg > 0.15
    x_fit = np.log(dg[mask_fit])
    y_fit = np.log(A_v[mask_fit])

    if len(x_fit) > 5:
        coeffs = np.polyfit(x_fit, y_fit, 1)
        exp_fit = coeffs[0]
        cA_fit  = np.exp(coeffs[1])
        # Residua
        y_pred  = coeffs[0] * x_fit + coeffs[1]
        r2      = 1.0 - np.sum((y_fit - y_pred)**2) / np.sum((y_fit - np.mean(y_fit))**2)
        print(f"  Fit: A_tail ≈ {cA_fit:.4f} × (g₀ − g*)^{exp_fit:.4f}")
        print(f"  R² = {r2:.6f}")
        print(f"  Referencja ex57: wykładnik 4.12 ± 0.05")
        print(f"  Odchylenie od 4.12: {abs(exp_fit - 4.12)/4.12*100:.2f}%")
    else:
        cA_fit  = 0.35
        exp_fit = 4.12
        print("  Za mało punktów do fitu — używam wartości ex57")


# ============================================================
# TEST 5: Kolejne zera B_tail (warunek kwantowania dla mion, τ)
# ============================================================
print("\n[TEST 5] Kolejne zera B_tail — kwantowanie solitonów")
print("-" * 50)

# Gęsty skan B_tail(g₀) dla szerokiego zakresu
g0_wide = np.linspace(1.08, 3.80, 150)
B_wide  = np.array([analyze(g0)['B'] for g0 in g0_wide])

zero_g0 = []
for i in range(len(g0_wide) - 1):
    b0, b1 = B_wide[i], B_wide[i+1]
    if np.isfinite(b0) and np.isfinite(b1) and b0 * b1 < 0:
        try:
            gz = brentq(B_func, g0_wide[i], g0_wide[i+1], xtol=1e-6, maxiter=60)
            zero_g0.append(gz)
        except Exception:
            pass

print(f"  Znalezione zera B_tail(g₀):")
for iz, gz in enumerate(zero_g0):
    d = analyze(gz)
    label = ""
    if iz == 0: label = " ← elektron (hipoteza)"
    if iz == 1: label = " ← mion?"
    if iz == 2: label = " ← tau?"
    A_z = d['A']
    print(f"    Zero #{iz+1}: g₀ = {gz:.6f}  A_tail = {A_z:.6f}{label}")
    if not np.isnan(A_z) and not np.isnan(A_e):
        print(f"           (A/A_e)⁴ = {(A_z/A_e)**4:.4f}  [r₂₁={R21_PDG:.2f}, r₃₁={R31_PDG:.2f}]")


# ============================================================
# PODSUMOWANIE I STATUS O-J2
# ============================================================
print("\n" + "=" * 65)
print("PODSUMOWANIE ex60")
print("=" * 65)
print(f"  α = {ALPHA}   (substrat, prop:substrate-action, N0-2)")
print(f"  g* = {G_STAR:.6f}  (f(g*)=0, punkt osobliwy kinetyczny)")
print(f"  φ  = {PHI_G:.6f}  (złota proporcja)")
print()
print(f"  Wynik T1 — B_tail = 0 (elektron): g₀^e ≈ {g0_e_found:.4f}")
print(f"  Wynik T2 — złota proporcja:  g₀^μ/g₀^e ≈ {ratio_obs:.4f}  (φ={PHI_G:.4f})")
if len(valid) > 5:
    print(f"  Wynik T4 — wykładnik ogona:  ≈ {exp_fit:.4f}  (ex57: 4.12)")
print()
print(f"  Zera B_tail: {[f'{gz:.4f}' for gz in zero_g0]}")

# Wniosek o O-J2
print()
if len(sector_changes) > 0:
    sc_g = [(s[0]+s[1])/2 for s in sector_changes]
    closest = min(sc_g, key=lambda x: abs(x - 2.0))
    dist = abs(closest - 2.0)
    if dist < 0.15:
        print(f"  STATUS O-J2: CZĘŚCIOWE POTWIERDZENIE")
        print(f"    Granica sektora przy g₀ ≈ {closest:.3f} (odl. od 2.00: {dist:.3f})")
    else:
        print(f"  STATUS O-J2: NIEZWERYFIKOWANE w bieżącej rozdzielczości")
        print(f"    Najbliższa granica sektora: g₀ ≈ {closest:.3f}")
else:
    print(f"  STATUS O-J2: brak zidentyfikowanych granic sektorów")

print()
print("  WYNIKI KLUCZOWE do ex61 (predykcja tau):")
print(f"    A_tail(g₀^e = {G0_E_H}) = {A_e:.6f}")
print(f"    A_tail(g₀^μ = {G0_MU}) = {A_mu:.6f}")
if 'exp_fit' in dir() and 'cA_fit' in dir():
    print(f"    Fit: A_tail ≈ {cA_fit:.4f} × (g₀ − {G_STAR:.4f})^{exp_fit:.4f}")

# Eksport wyników do pliku do użycia przez ex61
import json
results_export = {
    "g0_e": G0_E_H,
    "g0_mu": G0_MU,
    "g0_e_calc": float(g0_e_found),
    "A_e": float(A_e) if not np.isnan(A_e) else None,
    "A_mu": float(A_mu) if not np.isnan(A_mu) else None,
    "exp_fit": float(exp_fit) if 'exp_fit' in dir() else 4.12,
    "cA_fit": float(cA_fit) if 'cA_fit' in dir() else 0.35,
    "g_star": float(G_STAR),
    "zero_B_g0": [float(gz) for gz in zero_g0],
    "r21_calc": float((A_mu/A_e)**4) if (not np.isnan(A_e) and not np.isnan(A_mu)) else None,
}

out_dir = Path(__file__).parent.parent / "plots"
out_dir.mkdir(exist_ok=True)

json_path = out_dir / "ex60_results.json"
with open(json_path, "w", encoding="utf-8") as f:
    json.dump(results_export, f, indent=2, ensure_ascii=False)
print(f"\n  Wyniki zapisane: {json_path}")

# Wykres (opcjonalny, jeśli matplotlib dostępny)
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle("TGP ex60 — Topologia sektorów solitonu i złota proporcja\n"
                 f"α={ALPHA}, g*={G_STAR:.4f}, V(g)=g³/3−g⁴/4  (N0-2, N0-5)",
                 fontsize=11)

    # P1: profile solitonów
    ax = axes[0, 0]
    for g0p, col, lbl in [(1.24, 'tab:blue', 'e⁻ (g₀=1.24)'),
                          (2.00, 'tab:orange', 'μ? (g₀=2.00)'),
                          (2.34, 'tab:green', 'τ? (g₀=2.34)'),
                          (3.00, 'tab:red', 'g₀=3.00')]:
        rv, gv, sv = solve_soliton(g0p, r_max=35)
        if sv and len(rv) > 5:
            ax.plot(rv, gv, color=col, label=lbl, lw=1.4)
    ax.axhline(1.0, color='k', ls='--', lw=0.8, label='próżnia g=1')
    ax.axhline(G_STAR, color='r', ls=':', lw=0.8, label=f'g*={G_STAR:.3f}')
    ax.set(xlabel='r', ylabel='g(r) = Φ/Φ₀', title='Profile solitonów TGP',
           xlim=(0, 35), ylim=(0.65, 3.1))
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # P2: A_tail(g₀)
    ax = axes[0, 1]
    if len(valid) > 3:
        ax.semilogy(g0_v, A_v, 'b.', ms=4, label='num.')
        g0_fit_line = np.linspace(G_STAR + 0.12, g0_v[-1], 200)
        A_fit_line  = cA_fit * (g0_fit_line - G_STAR)**exp_fit
        ax.semilogy(g0_fit_line, A_fit_line, 'r--', lw=1.4,
                    label=f'fit: {cA_fit:.3f}·(g₀−g*)^{exp_fit:.3f}')
    ax.axvline(G0_E_H, color='g', ls=':', lw=1.2, label=f'g₀^e={G0_E_H}')
    ax.axvline(G0_MU,  color='orange', ls=':', lw=1.2, label=f'g₀^μ={G0_MU}')
    ax.set(xlabel='g₀', ylabel='A_tail', title='Amplituda ogona oscylacyjnego')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # P3: B_tail(g₀) — warunek kwantowania
    ax = axes[1, 0]
    finite_mask = np.isfinite(B_wide)
    ax.plot(g0_wide[finite_mask], B_wide[finite_mask], 'b-', lw=1.3)
    ax.axhline(0, color='k', ls='--', lw=0.8)
    for iz, gz in enumerate(zero_g0[:4]):
        ax.axvline(gz, color='r', ls=':', lw=1.0, label=f'B=0 @ {gz:.3f}')
    ax.axvline(G0_E_H, color='g', ls='-', lw=1.4, label=f'g₀^e={G0_E_H}')
    ax.axvline(G0_MU,  color='orange', ls='-', lw=1.4, label=f'g₀^μ={G0_MU}')
    ax.set(xlabel='g₀', ylabel='B_tail(g₀)', title='Warunek kwantowania B_tail = 0')
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # P4: stosunek mas (A/A_e)⁴
    ax = axes[1, 1]
    if len(valid) > 3 and not np.isnan(A_e):
        ratios = (A_v / A_e)**4
        pos = ratios > 0
        ax.semilogy(g0_v[pos], ratios[pos], 'b.', ms=4)
    ax.axhline(R21_PDG, color='g', ls='--', lw=1.4, label=f'r₂₁={R21_PDG:.2f}')
    ax.axhline(R31_PDG, color='r', ls='--', lw=1.4, label=f'r₃₁={R31_PDG:.2f}')
    ax.axvline(G0_MU, color='orange', ls=':', lw=1.2, label=f'g₀^μ=2.00')
    ax.set(xlabel='g₀', ylabel='(A(g₀)/A(g₀^e))⁴',
           title='Stosunek mas z amplitudy ogona')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = out_dir / "ex60_topology_golden_ratio.png"
    plt.savefig(str(fig_path), dpi=130, bbox_inches='tight')
    plt.close()
    print(f"  Wykres: {fig_path}")
except ImportError:
    print("  matplotlib niedostępne — pominięto wykresy")

print("\nex60 ZAKOŃCZONY")
