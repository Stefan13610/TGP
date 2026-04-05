#!/usr/bin/env python3
"""
TGP ex61 — Predykcja masy tau z amplitudy ogona solitonu
=========================================================
Cel:
  Wyznaczyć g₀^τ (parametr strzału tau leptonu) i stosunek
  r₃₁ = m_τ/m_e z hipotezy O-J3:

      m ∝ A_tail(g₀)⁴   (hip:J-mass-Atail4)

  Kandydaci na g₀^τ (do zbadania):
    C1: Zero B_tail(g₀) — trzecie zero warunku kwantowania
    C2: Złota proporcja: g₀^τ = g₀^μ × φ
    C3: Złota proporcja: g₀^τ = g₀^e × φ²
    C4: Numeryczne rozwiązanie (A(g₀^τ)/A(g₀^e))⁴ = r₃₁^{PDG}
    C5: Skalowanie z ex60: g₀ taki że sektor n_bounce = 2

ODE solitonu TGP (Dodatek J, N0-5):
  f(g)·g'' + (2/r)·g' = V'(g)
  f(g) = 1 + 2α·ln(g),   α = 2
  V'(g) = g²(1 − g)

Dane wejściowe z ex60:
  g₀^e = 1.24,  A_tail(g₀^e) ≈ 0.287
  g₀^μ = 2.00,  A_tail(g₀^μ) ≈ 1.09
  (A_μ/A_e)⁴ ≈ 206.77 (1.6% od r₂₁)

Dane PDG 2024:
  m_e  = 0.51099895 MeV
  m_μ  = 105.6583755 MeV
  m_τ  = 1776.86    MeV
  r₂₁  = m_μ/m_e = 206.7682
  r₃₁  = m_τ/m_e = 3477.15
  r₃₂  = m_τ/m_μ = 16.8171

Otwarte problemy adresowane: O-J3  (Dodatek J §J.5)

Autor: TGP v1, sesja v34 (2026-03-28)
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from pathlib import Path
import json, warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ============================================================
# Stałe TGP (z aksjomatów N0)
# ============================================================
ALPHA   = 2.0
G_STAR  = np.exp(-1.0 / (2.0 * ALPHA))   # e^{−1/4} ≈ 0.77880
PHI_G   = (1.0 + np.sqrt(5.0)) / 2.0     # φ = 1.61803...

# Wartości referencyjne
G0_E    = 1.24      # hipoteza elektronu
G0_MU   = 2.00      # hipoteza mion / złota proporcja

# Stosunki mas (PDG 2024)
R21_PDG = 206.7682
R31_PDG = 3477.15
R32_PDG = 16.8171

# Numeryczne
R_MAX        = 250.0
R_TAIL_START = 100.0
RTOL         = 1e-10
ATOL         = 1e-12
DR_START     = 1e-5

print("TGP ex61 — Predykcja masy tau z amplitudy ogona solitonu")
print(f"  α={ALPHA}, g*={G_STAR:.6f}, φ={PHI_G:.6f}")
print(f"  g₀^e={G0_E}, g₀^μ={G0_MU}")
print(f"  R₃₁(PDG)={R31_PDG}, R₃₂(PDG)={R32_PDG}")
print("=" * 65)


# ============================================================
# Solver (identyczny jak ex60)
# ============================================================
def f_kin(g):
    return 1.0 + 2.0 * ALPHA * np.log(max(g, 1e-15))

def Vprime(g):
    return g * g * (1.0 - g)

def ode_rhs(r, y):
    g, gp = y
    if g < G_STAR * 0.998:
        return [gp, 0.0]
    fg = f_kin(g)
    if r < 1e-9:
        gpp = Vprime(g) / (3.0 * fg)
    else:
        gpp = (Vprime(g) - (2.0 / r) * gp) / fg
    return [gp, gpp]

def event_singular(r, y):
    return y[0] - G_STAR * 1.002
event_singular.terminal  = True
event_singular.direction = -1

def initial_conditions(g0):
    fg0  = f_kin(g0)
    gpp0 = Vprime(g0) / (3.0 * fg0)
    r0   = DR_START
    return r0, g0 + 0.5*gpp0*r0**2, gpp0*r0

def solve_soliton(g0, r_max=R_MAX):
    r0, g_ic, gp_ic = initial_conditions(g0)
    sol = solve_ivp(ode_rhs, [r0, r_max], [g_ic, gp_ic],
                    method='DOP853', rtol=RTOL, atol=ATOL,
                    dense_output=False, events=event_singular,
                    max_step=0.15)
    return sol.t, sol.y[0], (sol.status == 0)

def fit_tail(r_arr, g_arr):
    mask = (r_arr >= R_TAIL_START) & np.isfinite(g_arr)
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

def A_tail_of(g0):
    r, g, suc = solve_soliton(g0)
    if not suc:
        return np.nan
    A, B, C = fit_tail(r, g)
    return A

def B_tail_of(g0):
    r, g, suc = solve_soliton(g0)
    if not suc:
        return np.nan
    A, B, C = fit_tail(r, g)
    return B if not np.isnan(B) else np.nan


# ============================================================
# Wyznaczenie wartości bazowych
# ============================================================
print("\n[Krok 0] Wyznaczanie wartości bazowych")
print("-" * 50)

A_e  = A_tail_of(G0_E)
A_mu = A_tail_of(G0_MU)

print(f"  A_tail(g₀^e = {G0_E}) = {A_e:.6f}")
print(f"  A_tail(g₀^μ = {G0_MU}) = {A_mu:.6f}")

if not (np.isnan(A_e) or np.isnan(A_mu)):
    r21_check = (A_mu / A_e)**4
    print(f"  (A_μ/A_e)⁴ = {r21_check:.4f}   [PDG: {R21_PDG:.4f}]  Δ={abs(r21_check-R21_PDG)/R21_PDG*100:.2f}%")

# A_tail wymagane dla tau
if not np.isnan(A_e):
    A_tau_target = A_e * R31_PDG**0.25
    print(f"\n  Wymagane A_tail(g₀^τ) = A_e × r₃₁^(1/4) = {A_e:.4f} × {R31_PDG**0.25:.4f} = {A_tau_target:.4f}")
else:
    A_tau_target = 2.20  # szacowanie z ex57


# ============================================================
# KANDYDAT C1: trzecie zero B_tail
# ============================================================
print("\n[C1] Trzecie zero B_tail (warunek kwantowania tau)")
print("-" * 50)

# Skanowanie B w szerokim zakresie — szukamy 3+ zer
g0_scan = np.linspace(1.10, 4.50, 250)
B_scan  = np.array([B_tail_of(g0) for g0 in g0_scan])

zero_indices = []
for i in range(len(g0_scan) - 1):
    b0, b1 = B_scan[i], B_scan[i+1]
    if np.isfinite(b0) and np.isfinite(b1) and b0 * b1 < 0:
        zero_indices.append(i)

C1_zeros = []
for i in zero_indices:
    try:
        gz = brentq(
            lambda g: B_tail_of(g) if not np.isnan(B_tail_of(g)) else 1.0,
            g0_scan[i], g0_scan[i+1], xtol=1e-6, maxiter=60
        )
        C1_zeros.append(gz)
    except Exception:
        pass

print(f"  Znalezione zera B_tail: {len(C1_zeros)}")
for iz, gz in enumerate(C1_zeros):
    Azz = A_tail_of(gz)
    labels = ["elektron", "?", "tau?"]
    lbl = labels[iz] if iz < len(labels) else ""
    print(f"    #{iz+1}: g₀ = {gz:.6f}  A_tail={Azz:.6f}  ({lbl})")
    if not np.isnan(Azz) and not np.isnan(A_e):
        ratio31 = (Azz/A_e)**4
        print(f"         (A/A_e)⁴ = {ratio31:.4f}  [r₃₁={R31_PDG:.2f}]  Δ={abs(ratio31-R31_PDG)/R31_PDG*100:.1f}%")
        ratio32 = (Azz/A_mu)**4 if not np.isnan(A_mu) else np.nan
        if not np.isnan(ratio32):
            print(f"         (A/A_μ)⁴ = {ratio32:.4f}  [r₃₂={R32_PDG:.2f}]  Δ={abs(ratio32-R32_PDG)/R32_PDG*100:.1f}%")

g0_tau_C1 = C1_zeros[2] if len(C1_zeros) >= 3 else None


# ============================================================
# KANDYDAT C2: złota proporcja g₀^τ = g₀^μ × φ
# ============================================================
print("\n[C2] Złota proporcja: g₀^τ = g₀^μ × φ")
print("-" * 50)

g0_tau_C2 = G0_MU * PHI_G
print(f"  g₀^τ = {G0_MU} × {PHI_G:.4f} = {g0_tau_C2:.6f}")
A_C2 = A_tail_of(g0_tau_C2)
print(f"  A_tail(g₀^τ) = {A_C2:.6f}")
if not np.isnan(A_C2) and not np.isnan(A_e):
    r31_C2 = (A_C2/A_e)**4
    r32_C2 = (A_C2/A_mu)**4 if not np.isnan(A_mu) else np.nan
    print(f"  (A_τ/A_e)⁴ = {r31_C2:.4f}  [r₃₁={R31_PDG:.2f}]  Δ={abs(r31_C2-R31_PDG)/R31_PDG*100:.1f}%")
    if not np.isnan(r32_C2):
        print(f"  (A_τ/A_μ)⁴ = {r32_C2:.4f}  [r₃₂={R32_PDG:.2f}]  Δ={abs(r32_C2-R32_PDG)/R32_PDG*100:.1f}%")


# ============================================================
# KANDYDAT C3: złota proporcja g₀^τ = g₀^e × φ²
# ============================================================
print("\n[C3] Złota proporcja: g₀^τ = g₀^e × φ²")
print("-" * 50)

g0_tau_C3 = G0_E * PHI_G**2
print(f"  g₀^τ = {G0_E} × φ² = {G0_E} × {PHI_G**2:.4f} = {g0_tau_C3:.6f}")
A_C3 = A_tail_of(g0_tau_C3)
print(f"  A_tail(g₀^τ) = {A_C3:.6f}")
if not np.isnan(A_C3) and not np.isnan(A_e):
    r31_C3 = (A_C3/A_e)**4
    r32_C3 = (A_C3/A_mu)**4 if not np.isnan(A_mu) else np.nan
    print(f"  (A_τ/A_e)⁴ = {r31_C3:.4f}  [r₃₁={R31_PDG:.2f}]  Δ={abs(r31_C3-R31_PDG)/R31_PDG*100:.1f}%")
    if not np.isnan(r32_C3):
        print(f"  (A_τ/A_μ)⁴ = {r32_C3:.4f}  [r₃₂={R32_PDG:.2f}]  Δ={abs(r32_C3-R32_PDG)/R32_PDG*100:.1f}%")


# ============================================================
# KANDYDAT C4: rozwiązanie odwrotne (A(g₀^τ)/A(g₀^e))⁴ = r₃₁
# ============================================================
print("\n[C4] Rozwiązanie odwrotne: (A(g₀)/A(g₀^e))⁴ = r₃₁^{PDG}")
print("-" * 50)

def objective_r31(g0):
    A = A_tail_of(g0)
    if np.isnan(A) or np.isnan(A_e) or A_e == 0:
        return np.nan
    return (A / A_e)**4 - R31_PDG

# Szukamy g₀ ∈ [2.0, 4.0]
g0_lo, g0_hi = 2.05, 4.0
f_lo = objective_r31(g0_lo)
f_hi = objective_r31(g0_hi)
print(f"  Sprawdzenie zakresu:")
print(f"    objective({g0_lo:.2f}) = {f_lo:.2f}")
print(f"    objective({g0_hi:.2f}) = {f_hi:.2f}")

g0_tau_C4 = None
if np.isfinite(f_lo) and np.isfinite(f_hi) and f_lo * f_hi < 0:
    try:
        g0_tau_C4 = brentq(objective_r31, g0_lo, g0_hi, xtol=1e-5, maxiter=80)
        A_C4 = A_tail_of(g0_tau_C4)
        r31_C4 = (A_C4 / A_e)**4
        r32_C4 = (A_C4 / A_mu)**4 if not np.isnan(A_mu) else np.nan
        print(f"\n  → g₀^τ (odwrotne) = {g0_tau_C4:.6f}")
        print(f"    A_tail(g₀^τ) = {A_C4:.6f}  [target: {A_tau_target:.6f}]")
        print(f"    (A_τ/A_e)⁴ = {r31_C4:.4f}  [r₃₁={R31_PDG:.4f}]")
        if not np.isnan(r32_C4):
            print(f"    (A_τ/A_μ)⁴ = {r32_C4:.4f}  [r₃₂={R32_PDG:.4f}]")
    except Exception as e:
        print(f"  Błąd brentq: {e}")
else:
    # Skan by znaleźć przybliżony g₀
    print("  Nie znaleziono zmiany znaku — skan gęsty:")
    g0_skan = np.linspace(1.8, 4.2, 50)
    f_skan = np.array([objective_r31(g) for g in g0_skan])
    sign_chg = np.where(np.diff(np.sign([f for f in f_skan if np.isfinite(f)])) != 0)[0]
    if len(sign_chg) > 0:
        g0_tau_C4 = g0_skan[sign_chg[0]]
        print(f"  Przybliżona wartość: g₀^τ ≈ {g0_tau_C4:.3f}")
    else:
        # Interpolacja z fit
        valid_f = [(g0_skan[i], f_skan[i]) for i in range(len(g0_skan))
                   if np.isfinite(f_skan[i])]
        if valid_f:
            g0_vals = np.array([x[0] for x in valid_f])
            f_vals  = np.array([x[1] for x in valid_f])
            closest = g0_vals[np.argmin(np.abs(f_vals))]
            A_closest = A_tail_of(closest)
            r31_closest = (A_closest/A_e)**4 if not np.isnan(A_closest) else np.nan
            print(f"  Minimalna odległość od r₃₁: g₀ ≈ {closest:.3f}, (A/A_e)⁴ ≈ {r31_closest:.1f}")
            g0_tau_C4 = closest


# ============================================================
# KANDYDAT C5: estymacja z fitowania A_tail ∝ (g₀ − g*)^n
# ============================================================
print("\n[C5] Estymacja analityczna z fitu potęgowego")
print("-" * 50)

# Wczytaj dane z ex60 jeśli dostępne
exp_val = 4.12    # domyślnie z ex57
cA_val  = 0.35

json_path = Path(__file__).parent.parent / "plots" / "ex60_results.json"
if json_path.exists():
    with open(json_path) as fp:
        ex60 = json.load(fp)
    A_e_ex60 = ex60.get("A_e") or A_e
    exp_val  = ex60.get("exp_fit") or 4.12
    cA_val   = ex60.get("cA_fit") or 0.35
    print(f"  Wczytano ex60: A_e={A_e_ex60:.6f}, wykładnik={exp_val:.4f}, cA={cA_val:.4f}")
else:
    print(f"  Używam wartości ex57: wykładnik={exp_val}, cA={cA_val}")

# Cel: cA × (g₀ − g*)^exp_val = A_tau_target
dg_tau_target = (A_tau_target / cA_val) ** (1.0 / exp_val)
g0_tau_C5 = G_STAR + dg_tau_target
print(f"  A_tau_target = {A_tau_target:.4f}")
print(f"  Wymagane (g₀ − g*) = {dg_tau_target:.4f}")
print(f"  → g₀^τ (analityczne) = {g0_tau_C5:.6f}")
A_C5 = A_tail_of(g0_tau_C5)
if not np.isnan(A_C5) and not np.isnan(A_e):
    r31_C5 = (A_C5/A_e)**4
    print(f"  Weryfikacja: (A({g0_tau_C5:.3f})/A_e)⁴ = {r31_C5:.4f}  [r₃₁={R31_PDG:.4f}]")
    print(f"  Odchylenie: {abs(r31_C5-R31_PDG)/R31_PDG*100:.2f}%")


# ============================================================
# TABLICA PORÓWNAWCZA
# ============================================================
print("\n" + "=" * 65)
print("TABLICA PORÓWNAWCZA KANDYDATÓW g₀^τ")
print("=" * 65)
print(f"{'Kandydat':<12} {'g₀^τ':>9} {'A_tail':>10} {'(A/Ae)⁴':>12} {'Δr₃₁%':>8} {'(A/Aμ)⁴':>12} {'Δr₃₂%':>8}")
print("-" * 65)

candidates = [
    ("C2 g₀^μ×φ", g0_tau_C2),
    ("C3 g₀^e×φ²", g0_tau_C3),
    ("C5 fit inv.", g0_tau_C5),
]
if g0_tau_C4 is not None:
    candidates.append(("C4 odwrotne", g0_tau_C4))
if g0_tau_C1 is not None:
    candidates.append(("C1 B=0#3", g0_tau_C1))

best_cand = None
best_err31 = np.inf

for name, g0t in candidates:
    if g0t is None:
        continue
    At = A_tail_of(g0t)
    if np.isnan(At) or np.isnan(A_e):
        print(f"  {name:<12} {g0t:>9.5f}  {'---':>10}")
        continue
    r31 = (At / A_e)**4
    r32 = (At / A_mu)**4 if not np.isnan(A_mu) else np.nan
    d31 = abs(r31 - R31_PDG)/R31_PDG*100
    d32 = abs(r32 - R32_PDG)/R32_PDG*100 if not np.isnan(r32) else np.nan
    d32_str = f"{d32:>8.1f}" if not np.isnan(d32) else "    ---"
    r32_str = f"{r32:>12.4f}" if not np.isnan(r32) else "         ---"
    print(f"  {name:<12} {g0t:>9.5f} {At:>10.6f} {r31:>12.4f} {d31:>8.1f} {r32_str} {d32_str}")
    if d31 < best_err31:
        best_err31 = d31
        best_cand  = (name, g0t, At, r31)

print("-" * 65)
print(f"  PDG cel:                              {R31_PDG:>12.4f}          {R32_PDG:>12.4f}")


# ============================================================
# WERYFIKACJA KOIDEGO — sprawdzenie formuly Koidego dla (e,μ,τ)
# ============================================================
print("\n[Bonus] Weryfikacja formuły Koidego")
print("-" * 50)
# Formuła Koidego: (√me + √mμ + √mτ)² = (2/3)(me + mμ + mτ)
# Używamy mas z PDG
m_e_pdg  = 0.51099895   # MeV
m_mu_pdg = 105.6583755
m_tau_pdg = 1776.86

LHS = (np.sqrt(m_e_pdg) + np.sqrt(m_mu_pdg) + np.sqrt(m_tau_pdg))**2
RHS = (2.0/3.0) * (m_e_pdg + m_mu_pdg + m_tau_pdg)
print(f"  PDG: (√me + √mμ + √mτ)² = {LHS:.6f}")
print(f"       (2/3)(me + mμ + mτ) = {RHS:.6f}")
print(f"  LHS/RHS = {LHS/RHS:.8f}  (idealne Koide: 1.000000)")

# Sprawdź predykcję tau z najlepszego kandydata
if best_cand is not None:
    name_b, g0t_b, At_b, r31_b = best_cand
    m_tau_pred = m_e_pdg * r31_b
    LHS_pred = (np.sqrt(m_e_pdg) + np.sqrt(m_mu_pdg) + np.sqrt(m_tau_pred))**2
    RHS_pred = (2.0/3.0) * (m_e_pdg + m_mu_pdg + m_tau_pred)
    koide_ratio = LHS_pred / RHS_pred
    print(f"\n  Predykcja ({name_b}): m_τ = {m_tau_pred:.2f} MeV")
    print(f"  Koide ratio (predykcja): {koide_ratio:.8f}")
    print(f"  Odchylenie od Koide=1: {abs(koide_ratio-1)*100:.4f}%")
    print(f"  Odchylenie m_τ od PDG: {abs(m_tau_pred-m_tau_pdg)/m_tau_pdg*100:.2f}%")


# ============================================================
# PODSUMOWANIE I STATUS O-J3
# ============================================================
print("\n" + "=" * 65)
print("PODSUMOWANIE ex61 — STATUS O-J3")
print("=" * 65)

if best_cand is not None:
    name_b, g0t_b, At_b, r31_b = best_cand
    print(f"  NAJLEPSZY KANDYDAT: {name_b}")
    print(f"    g₀^τ = {g0t_b:.6f}")
    print(f"    A_tail(g₀^τ) = {At_b:.6f}")
    print(f"    r₃₁ = (A_τ/A_e)⁴ = {r31_b:.4f}   [PDG: {R31_PDG:.4f}]")
    print(f"    Odchylenie od PDG: {best_err31:.2f}%")

print()
print(f"  Ciąg g₀: e={G0_E} → μ={G0_MU} → τ≈{g0t_b:.3f}" if best_cand else "")
print(f"  Stosunek g₀^μ/g₀^e = {G0_MU/G0_E:.4f}  (φ = {PHI_G:.4f})")
if best_cand:
    print(f"  Stosunek g₀^τ/g₀^μ = {g0t_b/G0_MU:.4f}")
    print(f"  Stosunek g₀^τ/g₀^e = {g0t_b/G0_E:.4f}  (φ² = {PHI_G**2:.4f})")

print()
print("  STATUS:")
if best_err31 < 5.0:
    print(f"  [PROPOZYCJA] g₀^τ ≈ {g0t_b:.4f} reprodukuje r₃₁ z dokładnością {best_err31:.1f}%")
    print(f"  → O-J3 CZĘŚCIOWO ZAMKNIĘTY: predykcja ilościowa istnieje,")
    print(f"    ale zasada selekcji g₀^τ pozostaje otwarta.")
elif best_err31 < 20.0:
    print(f"  [SZKIC] Najlepsza predykcja odchyla się o {best_err31:.1f}%")
    print(f"  → O-J3 OTWARTY: potrzebna dodatkowa zasada selekcji")
else:
    print(f"  [PROGRAM] Żaden kandydat nie reprodukuje r₃₁ z akceptowalną dokładnością")
    print(f"  → O-J3 OTWARTY: mechanizm wyznaczania g₀^τ nieznany")

# Eksport wyników
out_dir = Path(__file__).parent.parent / "plots"
out_dir.mkdir(exist_ok=True)
results = {
    "g0_e": G0_E, "g0_mu": G0_MU,
    "A_e": float(A_e) if not np.isnan(A_e) else None,
    "A_mu": float(A_mu) if not np.isnan(A_mu) else None,
    "candidates": {
        name: {
            "g0_tau": float(g0t),
            "A_tau": float(A_tail_of(g0t)) if not np.isnan(A_tail_of(g0t)) else None,
            "r31": float((A_tail_of(g0t)/A_e)**4)
                   if (not np.isnan(A_tail_of(g0t)) and not np.isnan(A_e)) else None
        }
        for name, g0t in candidates if g0t is not None
    },
    "best_candidate": best_cand[0] if best_cand else None,
    "g0_tau_best": float(best_cand[1]) if best_cand else None,
    "r31_best": float(best_cand[3]) if best_cand else None,
    "r31_pdg": R31_PDG,
    "r21_check": float(r21_check) if 'r21_check' in dir() and not np.isnan(r21_check) else None,
}

with open(out_dir / "ex61_results.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)
print(f"\n  Wyniki zapisane: {out_dir / 'ex61_results.json'}")
print("\nex61 ZAKOŃCZONY")
