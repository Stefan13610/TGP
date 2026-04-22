"""
ex61_atail_substrate_coupling.py
==================================
STATUS: LEGACY-TRANSLATIONAL

This script captures an older analytic route to `A_tail^4` and the historical
bridge from substrate language to lepton ratios. It remains useful as context,
but it is not the canonical synchronized bridge used by current `nbody`.

Analityczne i numeryczne wyprowadzenie A_tail^4 = r_21 z pierwszych zasad TGP.

KONTEKST (dod. J, hip:J-mass-Atail4):
  Hipoteza: m ∝ A_tail(g0)^4
  Wynik ex57: (A(2.0)/A(1.24))^4 = 210 ≈ 207 (1.6%)
  ODE TGP: f(g)·g'' + (2/r)·g' = V'(g)
  f(g) = 1 + 2α·ln(g), V(g) = g³/3 - g⁴/4 (β=γ, α=2)

CEL EX61:
  1. Weryfikacja skalowania A_tail ∝ (g0-g*)^p przez precyzyjny skan
  2. Wyznaczenie analitycznego wykładnika p (heurystyka WKB)
  3. Obliczenie "krytycznej" proporcji g0_mu/g0_e = r_21^(1/4p)
     dając predykcję analityczną (niezależną od ex57/ex58)
  4. Porównanie predykcji r_21^(1/4p) z φ i z wartością empiryczną 2.0/1.24
  5. WKB całka przy barierze kinetycznej: estymacja A_tail(g0) → f(g0)
  6. Test: czy A_tail(g0)^(1/4p) ∝ g0 (liniowe skalowanie masy z g0)?

ANALITYCZNE SZACOWANIE A_tail (ODE WKB, dod. J, sekcja app:J-analityczne):
  W pobliżu bariery f(g)=0 (g=g*), trajetoria g(r) przebywa czas T_bar ∝ 1/√ε
  gdzie ε = g0 - g* jest miarą "nadwyżki" ponad barierę.
  Amplituda oscylacyjna ogona ~ exp(-S_WKB), gdzie S_WKB jest akcją WKB.
  Dla V(g) ~ g^4 (dominujący człon): S_WKB ∝ (g0/g*)^4 → A_tail ∝ exp(-c·(g*/g0)^4)
  Ale to jest zbyt ogólne. Numerycznie A_tail ∝ g0^4.12.

DEFINICJE:
  g* = exp(-1/(2α)) ≈ 0.7788   (bariera kinetyczna)
  A_tail(g0) = sqrt(B^2 + C^2) z g(r)-1 ~ (B·cos(r)+C·sin(r))/r dla r→∞

TESTY (5):
  T1: Fit log-log A_tail ~ g0^p: p ∈ [3.5, 5.0] (z ex57: p=4.12)
  T2: r_21^(1/(4p)) ∈ [1.55, 1.70] → blisko φ=1.618 (test złotej proporcji)
  T3: g0_e·r_21^(1/4p) ∈ [1.90, 2.15] → blisko g0_mu=2.0
  T4: (A(g0_e·φ)/A(g0_e))^4 ∈ [180, 240] (odchylenie max 20% od 207)
  T5: Fit jakości R² > 0.99 (liniowe skalowanie w log-log)

Sesja: TGP v33+, 2026-03-28
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import os, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# Stałe TGP
# ─────────────────────────────────────────────────────────────────────────────
ALPHA    = 2.0
G_GHOST  = np.exp(-1.0 / (2.0 * ALPHA))   # ~0.7788
PHI_GOLD = (1.0 + np.sqrt(5.0)) / 2.0
R_21     = 206.768

# Parametry referencyjne
G0_E   = 1.24;  A_E   = 0.287
G0_MU  = 2.00;  A_MU  = 1.09

# Parametry ODE
G_BOUNCE = G_GHOST + 0.005
R_MAX    = 40.0
R_START  = 1e-4
R_TAIL_L = 20.0
R_TAIL_R = 36.0
RTOL     = 1e-10
ATOL     = 1e-13
MAX_STEP = 0.015

print("=" * 70)
print("EX61: ANALITYCZNE WYPROWADZENIE A_tail^4 = r21 W TGP")
print("=" * 70)
print(f"  g* = {G_GHOST:.5f},  φ = {PHI_GOLD:.6f}")
print(f"  Cel: r21 = {R_21}")
print(f"  g0_e = {G0_E}, g0_mu = {G0_MU}")
print(f"  Predykcja złotej serii: g0_e·φ = {G0_E*PHI_GOLD:.4f}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# ODE i integracja (precyzyjna)
# ─────────────────────────────────────────────────────────────────────────────
def rhs_tgp(r, y):
    g, gp = y
    g  = max(g, G_BOUNCE + 1e-8)
    fg = 1.0 + 2.0 * ALPHA * np.log(g)
    fg_eff = max(abs(fg), 1e-9)
    Vp = g**2 - g**3
    if r < 1e-4:
        gpp = Vp / (3.0 * fg_eff)
    else:
        gpp = (Vp - (2.0/r)*gp*fg_eff) / fg_eff
    return [gp, gpp]

def measure_tail(g0, r_max=R_MAX, r_tail_l=R_TAIL_L, r_tail_r=R_TAIL_R):
    """Precyzyjny pomiar A_tail, B_tail, C_tail."""
    y0  = [g0, 0.0]
    sol = solve_ivp(rhs_tgp, [R_START, r_max], y0,
                    dense_output=True, rtol=RTOL, atol=ATOL,
                    max_step=MAX_STEP, method='RK45')
    if not sol.success:
        return np.nan, np.nan, np.nan
    r_pts = np.linspace(r_tail_l, r_tail_r, 400)
    g_pts = sol.sol(r_pts)[0]
    dg    = g_pts - 1.0
    B_mat = np.column_stack([np.cos(r_pts)/r_pts, np.sin(r_pts)/r_pts])
    try:
        coeffs, _, _, _ = np.linalg.lstsq(B_mat, dg, rcond=None)
        B, C = coeffs
    except Exception:
        return np.nan, np.nan, np.nan
    A = np.sqrt(B**2 + C**2)
    return B, C, A

# ─────────────────────────────────────────────────────────────────────────────
# 1. Precyzyjny skan g0 ∈ [1.1, 2.5], N=80 (sektor n=0)
# ─────────────────────────────────────────────────────────────────────────────
print("─" * 60)
print("CZĘŚĆ 1: Precyzyjny skan A_tail(g0) w sektorze n=0")
print("─" * 60)

N_SCAN = 80
g0_arr = np.linspace(1.10, 2.50, N_SCAN)
A_arr  = np.zeros(N_SCAN)
B_arr  = np.zeros(N_SCAN)

for i, g0 in enumerate(g0_arr):
    B, C, A = measure_tail(g0)
    A_arr[i] = A
    B_arr[i] = B
    if i % 16 == 0:
        print(f"  g0={g0:.3f}  A={A:.4f}  B={B:.4f}")

# Sprawdź wartości referencyjne
B_e, C_e, A_e_meas = measure_tail(G0_E)
B_mu, C_mu, A_mu_meas = measure_tail(G0_MU)
B_phi, C_phi, A_phi_meas = measure_tail(G0_E * PHI_GOLD)

print()
print(f"  Zmierzone wartości referencyjne:")
print(f"    A(g0=1.24) = {A_e_meas:.4f}  (ex57: {A_E})")
print(f"    A(g0=2.00) = {A_mu_meas:.4f}  (ex57: {A_MU})")
print(f"    A(g0·φ={G0_E*PHI_GOLD:.3f}) = {A_phi_meas:.4f}")
if not np.isnan(A_e_meas) and A_e_meas > 0:
    r21_ex61 = (A_mu_meas/A_e_meas)**4
    r21_phi  = (A_phi_meas/A_e_meas)**4
    print(f"    (A_mu/A_e)^4       = {r21_ex61:.2f}  (cel: {R_21})")
    print(f"    (A(g0·φ)/A_e)^4   = {r21_phi:.2f}")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Fit log-log: A_tail ∝ g0^p
# ─────────────────────────────────────────────────────────────────────────────
print()
print("─" * 60)
print("CZĘŚĆ 2: Fit log-log A_tail ~ c·g0^p")
print("─" * 60)

valid = (~np.isnan(A_arr)) & (A_arr > 0.02) & (A_arr < 5.0)
g0_fit = g0_arr[valid]
A_fit  = A_arr[valid]

def power_law(g0, logc, p):
    return logc + p * np.log(g0)

try:
    popt_raw, pcov_raw = curve_fit(power_law, g0_fit, np.log(A_fit),
                                    p0=[0.0, 4.0], maxfev=5000)
    logc_fit, p_fit = popt_raw
    c_fit = np.exp(logc_fit)
    # R²
    residuals = np.log(A_fit) - power_law(g0_fit, *popt_raw)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((np.log(A_fit) - np.mean(np.log(A_fit)))**2)
    R2_raw = 1.0 - ss_res/ss_tot if ss_tot > 0 else 0

    print(f"  Fit (g0 całość): A ~ {c_fit:.4f}·g0^{p_fit:.3f}  (R²={R2_raw:.5f})")
except Exception as e:
    p_fit = 4.12  # fallback z ex57
    c_fit = 0.35
    R2_raw = 0.0
    print(f"  Fit nie powiódł się ({e}), używam p=4.12 (ex57)")

# Fit na węższym zakresie [1.1, 2.0] (sektor n=0)
mask_n0 = (g0_fit >= 1.1) & (g0_fit <= 2.0)
if np.sum(mask_n0) >= 5:
    try:
        popt_n0, _ = curve_fit(power_law, g0_fit[mask_n0], np.log(A_fit[mask_n0]),
                                p0=[0.0, 4.0], maxfev=5000)
        logc_n0, p_n0 = popt_n0
        c_n0 = np.exp(logc_n0)
        residuals_n0 = np.log(A_fit[mask_n0]) - power_law(g0_fit[mask_n0], *popt_n0)
        ss_res_n0 = np.sum(residuals_n0**2)
        ss_tot_n0 = np.sum((np.log(A_fit[mask_n0]) - np.mean(np.log(A_fit[mask_n0])))**2)
        R2_n0 = 1.0 - ss_res_n0/ss_tot_n0 if ss_tot_n0 > 0 else 0
        print(f"  Fit (n=0, [1.1,2.0]): A ~ {c_n0:.4f}·g0^{p_n0:.3f}  (R²={R2_n0:.5f})")
        p_use = p_n0
    except Exception:
        p_use = p_fit
        R2_n0 = R2_raw
else:
    p_use = p_fit
    R2_n0 = R2_raw

# ─────────────────────────────────────────────────────────────────────────────
# 3. Analityczna predykcja: g0_mu/g0_e z r_21^{1/(4p)}
# ─────────────────────────────────────────────────────────────────────────────
print()
print("─" * 60)
print("CZĘŚĆ 3: Analityczna predykcja g0_mu/g0_e z r_21^{1/(4p)}")
print("─" * 60)

# Z A_tail ∝ g0^p → (A_mu/A_e)^4 = (g0_mu/g0_e)^{4p} = r_21
# Stąd: g0_mu/g0_e = r_21^{1/(4p)}
for p_val in [p_use, 4.12, 4.0, 4.5]:
    if p_val > 0:
        ratio_pred = R_21**(1.0/(4.0*p_val))
        g0_mu_pred = G0_E * ratio_pred
        delta_phi  = abs(ratio_pred - PHI_GOLD) / PHI_GOLD * 100
        delta_mu   = abs(g0_mu_pred - G0_MU) / G0_MU * 100
        print(f"  p = {p_val:.3f}:  r_21^(1/4p) = {ratio_pred:.4f}  "
              f"(φ={PHI_GOLD:.4f}, Δφ={delta_phi:.1f}%),  "
              f"g0_mu_pred = {g0_mu_pred:.4f}  (g0_mu=2.00, Δ={delta_mu:.1f}%)")

print()
print(f"  Pytanie kluczowe: czy p ≈ ln(r_21)/(4·ln(φ))?")
p_golden = np.log(R_21) / (4.0 * np.log(PHI_GOLD))
print(f"  p_złota = ln({R_21:.1f}) / (4·ln(φ)) = {p_golden:.4f}")
print(f"  Zmierzony p = {p_use:.4f}  (Δ = {abs(p_use-p_golden)/p_golden*100:.1f}%)")
print()
print(f"  Interpretacja: jeśli p = p_złota = {p_golden:.4f},")
print(f"  to r_21^(1/(4p)) = φ DOKŁADNIE → złota proporcja wynikałaby")
print(f"  z wykładnika skalowania A_tail.")

# ─────────────────────────────────────────────────────────────────────────────
# 4. WKB szacowanie wykładnika A_tail
# ─────────────────────────────────────────────────────────────────────────────
print()
print("─" * 60)
print("CZĘŚĆ 4: Szacowanie WKB przy barierze kinetycznej")
print("─" * 60)

# Bariera: f(g) = 1 + 2α·ln(g) = 0 przy g = g* = exp(-1/(2α))
# V(g) = g³/3 - g⁴/4, V'(g) = g² - g³ = g²(1-g)
# W pobliżu bariery: dla g ~ g* + δg (δg mała)
#   f(g) ≈ 2α·δg/g* (Taylorek)
#   V(g) ≈ V(g*) + V'(g*)·δg
#
# Całka WKB (schematyczna) dla przejścia przez barierę:
# A_tail ~ exp(-∫ |κ(g)| dg), κ(g) ~ sqrt(|V_eff(g)|/|f(g)|)
# Dla V_eff ~ V'(g*)·δg i f ~ 2α·δg/g*:
# κ ~ sqrt(|V'(g*)| / (2α/g* · δg)) ~ (V'(g*)/g*)^{1/2} / sqrt(2α·δg)
# Całka od g* do g0: ∫_{g*}^{g0} dδg / sqrt(δg) ~ 2√(g0-g*)

g_star = G_GHOST
V_prime_star = g_star**2 * (1.0 - g_star)
print(f"  g* = {g_star:.5f}")
print(f"  V'(g*) = g*²(1-g*) = {V_prime_star:.5f}")
print(f"  2α/g* = {2*ALPHA/g_star:.5f}")
print()

# Dla A_tail ∝ (g0-g*)^{1/2} (WKB w pobliżu bariery):
# Jest to tylko przybliżenie "lokalne" przy barierze.
# Globalne skalowanie (g0 od g* do ~2): bardziej skomplikowane.
# Wynik numeryczny: A_tail ∝ g0^{4.12} ≈ g0^4.

# Szacowanie: dla potencjału V(g) ~ -g⁴/4 (dominujący przy małych g>g*):
# Akcja WKB: S = ∫_g*^g0 sqrt(|V(g)|/|f(g)|) dg
# V(g) ~ g³/3 - g⁴/4 ~ g³(1/3 - g/4) → dąży do 0 przy g=g*
# f(g) = 1 + 2α·ln(g) → dąży do 0 przy g=g*
# Stosunek: |V(g)|/|f(g)| ~ g³/(2α·|ln(g/g*)|) ≈ const dla g ~ g*

# Numeryczne szacowanie całki WKB
g_vals = np.linspace(g_star + 1e-5, 2.5, 1000)
V_vals = g_vals**3/3 - g_vals**4/4
f_vals = 1.0 + 2*ALPHA*np.log(g_vals)
integrand = np.abs(V_vals / (f_vals + 1e-10))**0.5

print("  Numeryczna całka WKB: S(g0) = ∫_{g*}^{g0} sqrt(|V/f|) dg")
g0_wkb_pts = [1.24, 1.50, 2.00, 2.50]
S_wkb = []
A_wkb = []
for g0_wkb in g0_wkb_pts:
    mask_wkb = g_vals <= g0_wkb
    if np.sum(mask_wkb) > 1:
        S = np.trapz(integrand[mask_wkb], g_vals[mask_wkb])
    else:
        S = 0.0
    S_wkb.append(S)
    A_wkb.append(np.exp(-S))

print("  g0     S_WKB   A_WKB_pred  A_tail_num")
for i, g0_wkb in enumerate(g0_wkb_pts):
    # Znormalizuj do g0_e=1.24
    a_pred_norm = A_wkb[i] / A_wkb[0] * A_e_meas if not np.isnan(A_e_meas) else np.nan
    B_tmp, C_tmp, A_num = measure_tail(g0_wkb)
    print(f"  {g0_wkb:.2f}  {S_wkb[i]:.4f}  {a_pred_norm:.4f}       {A_num:.4f}")

if len(S_wkb) >= 2 and S_wkb[0] > 0:
    S_ratio = S_wkb[-2] / S_wkb[0]  # S(2.0)/S(1.24)
    A_ratio_WKB = np.exp(-(S_wkb[-2] - S_wkb[0]))
    A_ratio_num = A_mu_meas / A_e_meas if not np.isnan(A_mu_meas) else np.nan
    print(f"\n  S(2.0)/S(1.24) = {S_ratio:.3f}")
    print(f"  A_WKB(2.0)/A_WKB(1.24) = exp(-(S_2-S_e)) = {A_ratio_WKB:.4f}")
    print(f"  A_num(2.0)/A_num(1.24)  = {A_ratio_num:.4f}")
    print(f"  (A_WKB_ratio)^4 = {A_ratio_WKB**4:.2f}  (cel: {R_21})")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Testy formalne
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 60)
print("TESTY FORMALNE")
print("=" * 60)

# T1: p ∈ [3.5, 5.0]
t1 = 3.5 <= p_use <= 5.0
print(f"T1 [{'PASS' if t1 else 'FAIL'}]: Wykładnik p = {p_use:.3f} ∈ [3.5, 5.0]")

# T2: r_21^(1/(4p)) ∈ [1.55, 1.70] → blisko φ
ratio_from_p = R_21**(1.0/(4.0*p_use))
t2 = 1.55 <= ratio_from_p <= 1.70
print(f"T2 [{'PASS' if t2 else 'FAIL'}]: r_21^(1/(4p)) = {ratio_from_p:.4f} ∈ [1.55, 1.70]  (φ={PHI_GOLD:.4f})")

# T3: g0_e · r_21^(1/(4p)) ∈ [1.90, 2.15]
g0_mu_pred_val = G0_E * ratio_from_p
t3 = 1.90 <= g0_mu_pred_val <= 2.15
print(f"T3 [{'PASS' if t3 else 'FAIL'}]: g0_e·r_21^(1/4p) = {g0_mu_pred_val:.4f} ∈ [1.90, 2.15]  (g0_mu=2.00)")

# T4: (A(g0_e·φ)/A(g0_e))^4 ∈ [180, 240]
t4 = False
r21_phi_val = np.nan
if not np.isnan(A_phi_meas) and not np.isnan(A_e_meas) and A_e_meas > 0:
    r21_phi_val = (A_phi_meas / A_e_meas)**4
    t4 = 180 <= r21_phi_val <= 240
print(f"T4 [{'PASS' if t4 else 'FAIL'}]: (A(g0_e·φ)/A(g0_e))^4 = {r21_phi_val:.1f} ∈ [180, 240]")

# T5: R² > 0.99
t5 = R2_raw > 0.99 or R2_n0 > 0.99
print(f"T5 [{'PASS' if t5 else 'FAIL'}]: R² = {max(R2_raw, R2_n0):.5f} > 0.99")

n_pass = sum([t1, t2, t3, t4, t5])
print()
print(f"WYNIK: {n_pass}/5 testów zdanych")

# ─────────────────────────────────────────────────────────────────────────────
# 6. Wykres
# ─────────────────────────────────────────────────────────────────────────────
print()
print("Generowanie wykresu...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle(
    f"EX61: Analityczne A_tail^4 = r21 w TGP\n"
    f"p = {p_use:.3f},  r21^(1/4p) = {ratio_from_p:.4f},  φ = {PHI_GOLD:.4f},  "
    f"Wyniki: {n_pass}/5",
    fontsize=11, fontweight='bold'
)

# Panel 1: A_tail(g0) log-log
ax = axes[0]
valid_p = (~np.isnan(A_arr)) & (A_arr > 0.01)
ax.plot(np.log(g0_arr[valid_p]), np.log(A_arr[valid_p]), 'k.', ms=3, label="dane")
g0_plot = np.linspace(1.1, 2.5, 100)
A_pred  = np.exp(logc_fit) * g0_plot**p_fit if 'logc_fit' in dir() else None
if A_pred is not None:
    ax.plot(np.log(g0_plot), np.log(A_pred), 'r-', lw=1.5,
            label=f"fit: $A \\propto g_0^{{{p_fit:.2f}}}$  R²={R2_raw:.3f}")
ax.axvline(np.log(G0_E),  color='blue',  lw=1.5, ls='--', label=f"ln g₀ᵉ={np.log(G0_E):.3f}")
ax.axvline(np.log(G0_MU), color='red',   lw=1.5, ls='--', label=f"ln g₀ᵘ={np.log(G0_MU):.3f}")
ax.axvline(np.log(G0_E*PHI_GOLD), color='green', lw=1.5, ls=':', label=f"ln g₀ᵉφ")
ax.set_xlabel("ln g₀"); ax.set_ylabel("ln A_tail")
ax.set_title("Log-log: A_tail vs g₀")
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

# Panel 2: (A(g0)/A_e)^4 vs g0 — predykcja masy
ax = axes[1]
A4_ratio = (A_arr / (A_e_meas if not np.isnan(A_e_meas) else A_E))**4
valid2 = valid_p & (A4_ratio < 500)
ax.plot(g0_arr[valid2], A4_ratio[valid2], 'k-', lw=1.0, label="$(A/A_e)^4$")
ax.axhline(R_21, color='orange', lw=2, ls='--', label=f"r_21={R_21:.1f}")
ax.axvline(G0_E,  color='blue', lw=1.5, ls='--')
ax.axvline(G0_MU, color='red',  lw=1.5, ls='--', label=f"g₀ᵘ=2.0")
ax.axvline(G0_E*PHI_GOLD, color='green', lw=1.5, ls=':', label=f"g₀ᵉφ={G0_E*PHI_GOLD:.3f}")
# Adnotacja wartości
for g0v, colv, lblv in [(G0_MU,'red',f'{(A_mu_meas/A_e_meas)**4:.1f}'),
                         (G0_E*PHI_GOLD,'green',f'{r21_phi_val:.1f}' if not np.isnan(r21_phi_val) else '?')]:
    idx_v = np.argmin(np.abs(g0_arr - g0v))
    if valid2[idx_v]:
        ax.annotate(lblv, xy=(g0v, A4_ratio[idx_v]),
                    xytext=(g0v+0.05, A4_ratio[idx_v]+5),
                    fontsize=8, color=colv)
ax.set_xlabel("g₀"); ax.set_ylabel("$(A_{\\rm tail}/A_e)^4$")
ax.set_title("Predykcja masy $(A/A_e)^4$ vs g₀")
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)
ax.set_ylim(0, min(A4_ratio[valid2].max()*1.1, R_21*2))

# Panel 3: Predykcja analityczna g0_mu vs p
ax = axes[2]
p_range = np.linspace(3.5, 5.5, 200)
ratio_pred_range = R_21**(1.0/(4.0*p_range))
g0_mu_range = G0_E * ratio_pred_range
ax.plot(p_range, ratio_pred_range, 'b-', lw=2, label="$r_{21}^{1/(4p)}$")
ax.axhline(PHI_GOLD, color='purple', lw=1.5, ls='--', label=f"φ={PHI_GOLD:.3f}")
ax.axhline(G0_MU/G0_E, color='red', lw=1.5, ls='--', label=f"g₀ᵘ/g₀ᵉ={G0_MU/G0_E:.3f}")
ax.axvline(p_use, color='orange', lw=1.5, ls='--', label=f"p_fit={p_use:.3f}")
ax.axvline(p_golden, color='green', lw=1.5, ls=':', label=f"p_φ={p_golden:.3f}")
# Zaznacz punkt przecięcia z φ
ax.plot(p_golden, PHI_GOLD, 'g*', ms=12, zorder=5, label=f"(p_φ, φ)")
ax.set_xlabel("p (wykładnik skalowania)"); ax.set_ylabel("$r_{21}^{1/(4p)}$ lub $g_0^\\mu/g_0^e$")
ax.set_title(f"Analityczna predykcja g₀ᵘ/g₀ᵉ\nφ = r21^(1/4p) przy p={p_golden:.3f}")
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)
ax.set_ylim(1.3, 2.1)

plt.tight_layout()
outdir = os.path.join(os.path.dirname(__file__), '..', 'plots')
os.makedirs(outdir, exist_ok=True)
outpath = os.path.join(outdir, 'ex61_atail_substrate_coupling.png')
plt.savefig(outpath, dpi=120, bbox_inches='tight')
print(f"  Wykres zapisany: {outpath}")
plt.close()

# ─────────────────────────────────────────────────────────────────────────────
# 7. Wnioski analityczne
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("WNIOSKI ANALITYCZNE (dod. J, sekcja app:J-analityczne)")
print("=" * 70)
print()
print(f"1. WYKŁADNIK SKALOWANIA A_tail:")
print(f"   Numeryczny p = {p_use:.3f} ± ~ 0.1")
print(f"   Heurystyka WKB: p ~ 4 (z V ~ g⁴)")
print(f"   p_złota = ln(r21)/(4·ln(φ)) = {p_golden:.4f}")
print(f"   Czy p = p_złota: Δ = {abs(p_use-p_golden)/p_golden*100:.1f}%")
print()
print(f"2. PREDYKCJA r21^(1/(4p)) PRZY p={p_use:.3f}:")
print(f"   r21^(1/(4·{p_use:.3f})) = {ratio_from_p:.4f}")
print(f"   Złota proporcja φ         = {PHI_GOLD:.4f}")
print(f"   Odchylenie od φ           = {abs(ratio_from_p-PHI_GOLD)/PHI_GOLD*100:.1f}%")
print()
print(f"3. KLUCZOWY WNIOSEK:")
print(f"   Jeśli p = ln(r21)/(4·ln(φ)) = {p_golden:.4f}, to złota proporcja")
print(f"   φ = r21^(1/(4p)) jest DOKŁADNIE spełniona.")
print(f"   Innymi słowy: r21 = φ^(4p)")
print(f"   Weryfikacja: φ^(4·{p_golden:.4f}) = {PHI_GOLD**(4*p_golden):.2f}  (cel: {R_21})")
print()
print(f"4. OTWARTE PYTANIE O-J1 (dod. J):")
print(f"   Czy p = p_złota jest wyprowadzalne z TGP?")
print(f"   Hipoteza: V(g) = g³/3 - g⁴/4 z V''(1) = -1 (faza metryczna)")
print(f"   → szczegółowa analiza WKB akcji przy g* jest potrzebna.")
print()
print(f"WYNIK KOŃCOWY: {n_pass}/5 testów zdanych")
print(f"Patrz: dodatekJ_ogon_masy.tex, problem O-J1, hipoteza hip:J-mass-Atail4")
