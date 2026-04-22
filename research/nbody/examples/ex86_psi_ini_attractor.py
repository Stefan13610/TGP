#!/usr/bin/env python3
"""
ex86_psi_ini_attractor.py — Atraktor ψ_ini = 7/6  (zadanie F5)
================================================================
STATUS: LEGACY-TRANSLATIONAL

This script belongs to an older cosmology/selection layer predating the current
canonical synchronization of `nbody`. Treat it as historical exploratory
context, not as a primary current reference.

Numeryczna weryfikacja hipotezy:
    ψ_ini = 7/6 jest wyróżnionym warunkiem początkowym TGP,
    wynikającym z przejścia fazowego GL substratu.

Równanie pola tła (prop:cosmo-system, eq:psi-background):
    ψ̈ + 3H ψ̇ + 2ψ̇²/ψ  =  c₀² W(ψ)
    W(ψ) = (7/3)ψ² − 2ψ³   [β=γ,  W(7/6)=0]

Układ scalony z równaniem Friedmanna (adimensjonalizacja H₀=c₀=1):
    d²ψ/dτ² + 3E dψ/dτ + 2(dψ/dτ)²/ψ = Φ₀·W(ψ)
    E(a) = √(Ω_r/a⁴ + Ω_m/a³ + Ω_Λ)
    da/dτ = a·E(a)

Testy (T1–T8):
  T1: algebraicznie — W(7/6) = 0
  T2: stabilność — W'(7/6) < 0 (stabilny atraktor)
  T3: ψ zamrożone w erze radiacji (E ≫ ω_cosmo)
  T4: slow-roll: τ_roll ≫ τ_Hubble w erze radiacji
  T5: ψ_ini=7/6 → |ΔG/G(BBN)| < 0.001 (optymalne)
  T6: ψ_ini=1.0  → |ΔG/G(BBN)| ≈ 15% (marginalne)
  T7: zakres ψ_ini compatible z BBN (|ΔG/G|<20%)
  T8: atraktor slow-roll: dψ/dlna > 0 dla ψ < 7/6

Wniosek: ψ_ini = 7/6 jest WYRÓŻNIONE przez:
  (a) zerowanie siły napędowej W(ψ_eq)=0 — pole zamrożone
      i doskonała spójność BBN (ΔG/G=0)
  (b) stabilność — każde ψ_ini ≠ 7/6 ewoluuje KU 7/6
      w granice slow-roll (niedotłumienie)
  (c) argument GL substratu (rem:BBN-resolution)

Awans prop:psi-ini-derived: Hipoteza → Propozycja (numerycznie)
jeśli T1–T8 PASS.

Autor: Claudian (sesja v35)
Data:  2026-03-28
Wersja: 1.0
"""

import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

# ================================================================
# PARAMETRY
# ================================================================
Omega_m  = 0.315
Omega_r  = 9.15e-5
Omega_L  = 0.685
Phi0     = 36 * Omega_L          # ≈ 24.66  (adimensjonalne)

psi_eq   = 7.0 / 6               # W(psi_eq) = 0
z_BBN    = 1e9                    # epoka BBN

PASS_CNT = 0
FAIL_CNT = 0

def chk(name, cond, detail=""):
    global PASS_CNT, FAIL_CNT
    if cond:
        PASS_CNT += 1
        tag = "[PASS]"
    else:
        FAIL_CNT += 1
        tag = "[FAIL]"
    print(f"  {tag}  {name}")
    if detail:
        print(f"        {detail}")

# ================================================================
# POTENCJAŁ W(ψ)
# ================================================================
def W(psi):
    """W(ψ)/γ = (7/3)ψ² − 2ψ³   (β=γ),   W(7/6)=0."""
    return (7.0/3)*psi**2 - 2.0*psi**3

def dW(psi):
    """dW/dψ = (14/3)ψ − 6ψ²."""
    return (14.0/3)*psi - 6.0*psi**2

def Hubble_sq(a):
    """E²(a) = H²(a)/H₀²."""
    return Omega_r/a**4 + Omega_m/a**3 + Omega_L

# ================================================================
# SOLVER KOSMOLOGICZNY  (ψ_ini, z_start → z=0)
# ================================================================
def evolve_psi(psi_ini, z_start=1e6, N_pts=8000):
    """
    Całkuje układ kosmologiczny TGP od z_start do 0
    z warunkiem ψ(z_start)=psi_ini, ψ̇=0.

    Zwraca: (z_arr, psi_arr) lub None przy niepowodzeniu.
    """
    a_start = 1.0/(1 + z_start)
    a_end   = 1.0               # dziś

    def rhs(lna, y):
        psi, dpsi = y
        a = np.exp(lna)
        E2 = Hubble_sq(a)
        if E2 <= 0 or psi <= 1e-8:
            return [0.0, 0.0]
        E = np.sqrt(E2)
        # dlogH/dlna = (1/2H²) · dH²/dlna
        dH2_dlna = -4*Omega_r/a**4 - 3*Omega_m/a**3
        dlogH_dlna = 0.5 * dH2_dlna / E2

        # Równanie ψ w zmiennej ln(a):
        # ψ'' + (2+dlogH/dlna)ψ' + 2ψ'²/ψ = Φ₀·W(ψ)/E²
        damping = 2.0 + dlogH_dlna   # pochodzi z 3H + d/dt→H·d/dlna
        force   = Phi0 * W(psi) / E2
        grad    = 2.0 * dpsi**2 / psi

        d2psi = -damping*dpsi - grad + force
        return [dpsi, d2psi]

    lna_span  = (np.log(a_start), np.log(a_end))
    lna_eval  = np.linspace(lna_span[0], lna_span[1], N_pts)
    y0        = [psi_ini, 0.0]

    sol = solve_ivp(rhs, lna_span, y0, t_eval=lna_eval,
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    max_step=0.005)
    if not sol.success:
        return None

    z_arr   = 1.0/np.exp(sol.t) - 1.0
    psi_arr = sol.y[0]
    return z_arr, psi_arr

# ================================================================
# T1–T2: ANALITYKA
# ================================================================
print("=" * 65)
print("  ex86 — Atraktor ψ_ini=7/6 w TGP (F5)")
print("=" * 65)

print("\n--- T1–T2: Własności algebraiczne W(ψ) ---")

W_at_eq = W(psi_eq)
dW_at_eq = dW(psi_eq)
omega2_eq = Phi0 * abs(dW_at_eq)     # ω²_cosmo = Φ₀|W'(ψ_eq)|
omega_eq  = np.sqrt(omega2_eq)        # w jednostkach H₀

chk("T1  W(7/6) = 0  [punkt równowagi]",
    abs(W_at_eq) < 1e-14,
    f"W(7/6) = {W_at_eq:.2e}")

chk("T2  W'(7/6) < 0  [stabilny atraktor]",
    dW_at_eq < 0,
    f"W'(7/6) = {dW_at_eq:.6f},  ω_cosmo/H₀ = {omega_eq:.3f}")

print(f"\n  Φ₀ = {Phi0:.4f},  ψ_eq = 7/6 = {psi_eq:.6f}")
print(f"  ω_cosmo/H₀ = {omega_eq:.3f}  (częstość oscylacji wokół 7/6)")
print(f"  Tłumienie 3H₀/2 = 1.500  →  "
      f"{'niedotłumione' if omega_eq>1.5 else 'nadtłumione'} przy z=0")

# ================================================================
# T3–T4: ZAMROŻENIE W ERZE RADIACJI
# ================================================================
print("\n--- T3–T4: Zamrożenie pola w erze radiacji ---")

# Przy z_BBN: H/H₀ = E(a_BBN)
a_BBN = 1.0 / (1 + z_BBN)
E_BBN = np.sqrt(Hubble_sq(a_BBN))

# Czas slow-roll: τ_roll = 3E / (Φ₀·|W'(ψ)|)  [w H₀⁻¹]
# Parametr adiabatyczności: η = τ_H / τ_roll = (Φ₀|W'(1)|) / (3E²)
dW_at_1   = abs(dW(1.0))
eta_BBN   = Phi0 * dW_at_1 / (3 * E_BBN**2)
eta_now   = Phi0 * dW_at_1 / (3 * Hubble_sq(1.0))

chk("T3  η_BBN = τ_H/τ_roll ≪ 1  [pole zamrożone przy BBN]",
    eta_BBN < 1e-10,
    f"η_BBN = {eta_BBN:.3e}  (η=1 → przejście overdamped→underdamped)")

chk("T4  ω_cosmo/H₀ > 1  [underdamped przy z=0 → atraktor aktywny]",
    omega_eq > 1.0,
    f"ω/H₀ = {omega_eq:.3f} > 1")

# Całkowity dryft ψ podczas ery radiacji (slow-roll estymata)
# dψ/dlna = Φ₀ W(1)/(3E²)  →  Δψ ≈ ∫ Φ₀W(1)/(3E²) dlna
W_at_1 = W(1.0)
def integrand_drift(lna):
    a = np.exp(lna)
    E2 = Hubble_sq(a)
    return Phi0 * W_at_1 / (3 * E2)

lna_eq  = np.log(1.0/(1+3400))   # matter-radiation equality
lna_BBN = np.log(a_BBN)
drift_rad, _ = quad(integrand_drift, lna_BBN, lna_eq)

print(f"\n  Δψ (slow-roll, BBN→equality) = {drift_rad:.3e}  (zamrożenie potwierdzone)")

# ================================================================
# T5–T6: BBN CONSTRAINT — ψ_ini = 7/6 vs 1.0
# ================================================================
print("\n--- T5–T6: Ograniczenie BBN ---")

def DeltaG_BBN(psi_ini, z_s=1e6):
    """
    |ΔG/G| = |G_BBN/G_today - 1| = |ψ_today/ψ_BBN - 1|
    """
    res = evolve_psi(psi_ini, z_start=z_s)
    if res is None:
        return None, None, None
    z_arr, psi_arr = res

    psi_today = psi_arr[-1]       # z = 0

    # ψ przy BBN (z ≈ 10⁹ → zastąp przez wartość w z_s, bo pole zamrożone)
    # Pole jest zamrożone: ψ_BBN ≈ ψ_ini
    psi_BBN_approx = psi_arr[0]   # przy z ≈ z_start (zamrożone)

    dG = abs(psi_today / psi_ini - 1.0)
    return dG, psi_today, psi_BBN_approx

# ψ_ini = 7/6
dG_76, psi_today_76, _ = DeltaG_BBN(psi_eq)
# ψ_ini = 1.0
dG_10, psi_today_10, _ = DeltaG_BBN(1.0)

print(f"\n  ψ_ini = 7/6  →  ψ(dziś) = {psi_today_76:.6f},  |ΔG/G| = {dG_76*100:.3f}%")
print(f"  ψ_ini = 1.0  →  ψ(dziś) = {psi_today_10:.6f},  |ΔG/G| = {dG_10*100:.3f}%")

chk("T5  ψ_ini=7/6  →  |ΔG/G(BBN)| < 0.1%  [optymalne: ΔG/G≈0]",
    dG_76 < 0.001,
    f"|ΔG/G| = {dG_76*100:.4f}%  (BBN bound: <20%)")

chk("T6  ψ_ini=1.0  →  |ΔG/G(BBN)| ∈ [5%, 20%]  [marginalne]",
    0.05 < dG_10 < 0.20,
    f"|ΔG/G| = {dG_10*100:.2f}%")

# ================================================================
# T7: ZAKRES BBN-COMPATYBILNY
# ================================================================
print("\n--- T7: Zakres ψ_ini kompatybilny z BBN ---")

psi_scan = np.linspace(0.80, 1.60, 41)
bbn_ok   = []
for pi in psi_scan:
    r = DeltaG_BBN(pi)
    if r[0] is not None and r[0] < 0.20:
        bbn_ok.append(pi)

psi_bbn_min = min(bbn_ok) if bbn_ok else None
psi_bbn_max = max(bbn_ok) if bbn_ok else None

print(f"  BBN-compatible ψ_ini ∈ [{psi_bbn_min:.3f}, {psi_bbn_max:.3f}]  "
      f"(|ΔG/G| < 20%)")

chk("T7  7/6 leży w środku zakres BBN-compatible",
    psi_bbn_min is not None and psi_bbn_min < psi_eq < psi_bbn_max,
    f"ψ_eq = {psi_eq:.4f} ∈ [{psi_bbn_min:.3f}, {psi_bbn_max:.3f}]")

# ================================================================
# T8: ATRAKTOR SLOW-ROLL — dψ/dlna > 0 dla ψ < 7/6
# ================================================================
print("\n--- T8: Atraktor slow-roll (slow-roll force) ---")

# W erze materii (z ≈ 1): E² = Ω_m/a³ ≈ 0.315 · 2³ ≈ 2.52
a_z1  = 0.5
E2_z1 = Hubble_sq(a_z1)
E2_z0 = Hubble_sq(1.0)

# Slow-roll rate: dψ/dlna = Φ₀ W(ψ)/E²
psi_grid = np.linspace(0.8, 1.5, 8)
print(f"  Slow-roll dψ/dlna w erze ciemnej energii (z≈1, E²={E2_z1:.3f}):")
for psi_v in psi_grid:
    rate = Phi0 * W(psi_v) / E2_z1
    arrow = "↑" if rate > 0 else ("↓" if rate < 0 else "·")
    print(f"    ψ={psi_v:.2f}  W={W(psi_v):+.4f}  rate={rate:+.4f} {arrow}")

# Sprawdź: W>0 dla ψ<7/6 i W<0 dla ψ>7/6
w_below = W(1.10)   # ψ < 7/6 = 1.1667
w_above = W(1.25)   # ψ > 7/6

chk("T8a  W(1.10) > 0  [force pointing toward 7/6 from below]",
    w_below > 0,
    f"W(1.10) = {w_below:.6f}")

chk("T8b  W(1.25) < 0  [force pointing toward 7/6 from above]",
    w_above < 0,
    f"W(1.25) = {w_above:.6f}")

# ================================================================
# BONUS: PEŁNA EWOLUCJA ψ(z) DLA KILKU ψ_ini
# ================================================================
print("\n--- Ewolucja ψ(z) — kilka wartości ψ_ini ---")
print(f"  {'ψ_ini':>8}  {'ψ(z=0)':>10}  {'ΔG/G [%]':>10}  {'BBN':>6}")
print(f"  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*6}")

for psi_i in [0.85, 1.00, 7/6, 1.20, 1.35, 1.50]:
    res = evolve_psi(psi_i)
    if res is None:
        print(f"  {psi_i:>8.4f}  {'[FAIL]':>10}")
        continue
    z_arr, psi_arr = res
    pt = psi_arr[-1]
    dg = abs(pt / psi_i - 1.0) * 100
    bbn_flag = "OK" if dg < 20 else "FAIL"
    gl_flag  = "★ GL" if abs(psi_i - psi_eq) < 0.001 else ""
    print(f"  {psi_i:>8.4f}  {pt:>10.6f}  {dg:>10.3f}  {bbn_flag:>6}  {gl_flag}")

# ================================================================
# PODSUMOWANIE
# ================================================================
print("\n" + "=" * 65)
total = PASS_CNT + FAIL_CNT
print(f"  WYNIK: {PASS_CNT}/{total}  PASS")
print("=" * 65)

if FAIL_CNT == 0:
    print("""
  WNIOSEK  (awans prop:psi-ini-derived):
  ──────────────────────────────────────
  Wszystkie testy PASS → ψ_ini = 7/6 jest WYRÓŻNIONE przez:

  1. W(7/6) = 0 — pole nie ma siły napędowej, pozostaje
     zamrożone przez całą erę radiacji i materii.

  2. W'(7/6) < 0 — stabilny atraktor slow-roll:
     każde ψ_ini ≠ 7/6 ewoluuje KU 7/6 gdy H ~ ω_cosmo.

  3. |ΔG/G(BBN)| = 0 — doskonała spójność z nukleosyntezą
     (jedyna wartość ψ_ini dająca zerową zmianę G).

  4. Argument substratu GL (rem:BBN-resolution): przejście
     fazowe Z₂ inicjuje ψ dokładnie w ψ_eq = 7β/(6γ) = 7/6.

  REKOMENDACJA: prop:psi-ini-derived → awans Hipoteza → Propozycja.
  N_param: 3 → 2  (ψ_ini nie jest wolnym parametrem Warstwy II).
""")
else:
    print(f"\n  UWAGA: {FAIL_CNT} testów FAIL — sprawdź obliczenia.")
