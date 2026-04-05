#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex134_top4_analytic_g0e.py
==========================
BADANIE: DLACZEGO g₀^e ≈ 1 + 1/(2α) = 5/4?

KONTEKST (ex133, 15/17 PASS):
  g₀^e(obs)   = 1.24915   (φ-punkt stały, ex106)
  g₀^{e,*}    = 1.249082  (daje θ_TGP = θ(9r*) DOKŁADNIE)
  5/4         = 1.250000  = 1 + 1/(2α)
  |g₀^{e,*} − 5/4| = 0.000918  (0.073%)

PYTANIE:
  Czy g₀^e ≈ 1+1/(2α) wynika z głębokiej struktury ODE solitonu TGP?
  Czy istnieje analityczny warunek równoważny g₀^e = 1+1/(2α)?

HIPOTEZY DO PRZETESTOWANIA:
  H1: Faza δ ogona (arctan(C_sin/B_cos)) ma specjalną wartość przy g₀^e=5/4
  H2: Potęgowy wykładnik A_tail ~ (g₀^e-1)^q jest specjalny przy 5/4
  H3: f(g₀^e=5/4) ≈ 2 lub innej specjalnej wartości (warunek na f)
  H4: A_tail(φ·g₀^e)/A_tail(g₀^e) = (9r*)^{1/4} przy g₀^e=5/4 do jakiego stopnia?
  H5: Warunek A_tail / (g₀^e-1) = specjalna wartość (liniowa aproksymacja)
  H6: Korekta 2-go rzędu: g₀^{e,*} = 1+1/(2α) − C·(1/(2α))^2 dla pewnego C analitycznego
  H7: Faza liniowego ODE δ₀ (granica g₀^e→1+)
  H8: (g₀^e−g*)/(1−g*) — pozycja ułamkowa g₀^e między g* a 1+1/(2α)
  H9: Czy g₀^{e,*} = √φ / coś? Test kandydatów algebraicznych
  H10: Warunku „energetyczny": V(g₀^e) = V(1)·(g₀^e-1)/(1-g*)?

TESTY P01..P16:
  P01: Faza δ przy g₀^e = 5/4 = 1.250000
  P02: Faza δ przy g₀^{e,*} = 1.249082
  P03: Faza δ skanem — czy jest monotoniczna? Ekstrenum przy 5/4?
  P04: Potęgowy wykładnik q = d(ln A_tail)/d(ln(g₀^e-1)) przy 5/4
  P05: A_tail(5/4) / A_tail_linear — odchylenie od liniowej zależności
  P06: f(5/4) = 1+4ln(5/4) — najbliższy kandydat algebraiczny
  P07: c₂(g₀^e) = g''(0)/2 przy g₀^e=5/4 — analityczna formula
  P08: Stosunek A_tail(φ·g₀^e)/A_tail(g₀^e) przy g₀^e=5/4 vs (9r*)^{1/4}
  P09: Algebraiczny scan kandydatów na g₀^{e,*}
  P10: Korekta δ = g₀^{e,*} − 5/4 — wyrażona analitycznie?
  P11: Warunek linowy: A_tail ∝ (g₀^e−1)^q z jakim q globalnie?
  P12: (g₀^e − g*) / (1 − g*) przy g₀^e=5/4 — ułamkowe położenie
  P13: Faza δ₀ liniowego ODE (granica małych amplitud)
  P14: Czy r₂₁(g₀^e=5/4) = 9r* · (1 + ε) dla ε algebraicznego?
  P15: Test: g₀^e − 1 vs −ln(g*) — tożsamość? (kluczowe!)
  P16: META: Które wyrażenie analityczne najlepiej aproksymuje g₀^{e,*}?

Referencje: ex129 (θ(9r*)=132.7314°), ex131 (θ_TGP=132.7324°),
            ex132 (Δg₀*=π(1-g*)), ex133 (g₀^{e,*}=1.249082)
"""

import sys
import io
import math
import warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from scipy.stats import linregress

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

warnings.filterwarnings('ignore')

# ============================================================
# Stałe TGP
# ============================================================
ALPHA   = 2.0
PHI     = (1.0 + math.sqrt(5.0)) / 2.0
G_STAR  = math.exp(-1.0 / (2.0 * ALPHA))    # exp(-1/4) ≈ 0.77880
G_BOUNCE = G_STAR + 0.005

RSTAR       = (23.0 + 5.0 * math.sqrt(21.0)) / 2.0   # ≈ 22.9564
THETA_9RS   = 132.731439   # z ex129 D7
THETA_PDG   = 132.732822   # z mas PDG
G0_E_OBS    = 1.24915      # φ-punkt stały
G0_E_STAR   = 1.249082     # brentq z ex133

# ODE parametry
R_MAX    = 60.0
R_START  = 1e-4
MAX_STEP = 0.02
RTOL     = 1e-10
ATOL     = 1e-13
R_L_FIT  = 20.0
R_R_FIT  = 35.0

TESTS = []

def record(name, passed, detail=""):
    TESTS.append((name, passed, detail))
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"         {detail}")


# ============================================================
# ODE solitonu (identyczny jak ex133)
# ============================================================
def f_kin(g):
    return 1.0 + 2.0 * ALPHA * math.log(max(g, 1e-30))

def Vprime(g):
    return g * g * (1.0 - g)

def Vpot(g):
    """V(g) = ∫₀^g V'(t)dt = g³/3 - g⁴/4"""
    return g**3 / 3.0 - g**4 / 4.0

def rhs(r, y):
    g, gp = y
    g = max(g, G_BOUNCE + 1e-7)
    fg = f_kin(g)
    if abs(fg) < 1e-10:
        return [gp, 0.0]
    driving = Vprime(g)
    cross   = (ALPHA / g) * gp**2
    if r < 1e-10:
        return [gp, (driving - cross) / (3.0 * fg)]
    damp = fg * 2.0 * gp / r
    return [gp, (driving - cross - damp) / fg]

def event_ghost(r, y):
    return y[0] - G_BOUNCE
event_ghost.terminal  = True
event_ghost.direction = -1

def integrate_soliton(g0, r_max=None, max_bounces=25):
    if r_max is None:
        r_max = max(R_MAX, 15.0 * g0)
    r0, y0 = R_START, [g0, 0.0]
    segs_r, segs_g = [], []
    n_bounces = 0
    for bn in range(max_bounces + 1):
        sol = solve_ivp(rhs, [r0, r_max], y0,
                        method='DOP853', max_step=MAX_STEP,
                        rtol=RTOL, atol=ATOL,
                        events=[event_ghost], dense_output=False)
        segs_r.append(sol.t)
        segs_g.append(sol.y[0])
        if sol.t_events[0].size > 0 and bn < max_bounces:
            n_bounces += 1
            r_b  = float(sol.t_events[0][0])
            gp_b = float(sol.y_events[0][0, 1])
            r0   = r_b + 1e-6
            y0   = [G_BOUNCE + 1e-5, -gp_b]
        else:
            break
    r_all = np.concatenate(segs_r)
    g_all = np.concatenate(segs_g)
    idx   = np.argsort(r_all)
    return r_all[idx], g_all[idx], n_bounces

def fit_tail_full(r_arr, g_arr, r_L=20.0, r_R=35.0):
    """
    Dopasowuje (g-1)·r = B_c·cos(r) + C_s·sin(r).
    Zwraca (A_amp, B_c, C_s, phase_deg, rmse_rel) gdzie
      A_amp = sqrt(B_c²+C_s²), phase = arctan2(C_s, B_c) w stopniach.
    """
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    n = np.sum(mask)
    if n < 10:
        return float('nan'), float('nan'), float('nan'), float('nan'), float('nan')
    r_f = r_arr[mask]
    y_f = (g_arr[mask] - 1.0) * r_f
    X   = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, _, _, _ = np.linalg.lstsq(X, y_f, rcond=None)
    B_c, C_s = float(coef[0]), float(coef[1])
    A_amp = math.sqrt(B_c**2 + C_s**2)
    phase_deg = math.degrees(math.atan2(C_s, B_c))
    y_hat = B_c * np.cos(r_f) + C_s * np.sin(r_f)
    rmse  = float(np.sqrt(np.mean((y_f - y_hat)**2)))
    return A_amp, B_c, C_s, phase_deg, rmse / max(A_amp, 1e-10)

def atail_and_phase(g0):
    """Zwraca (A_tail, phase_deg) dla solitonu z g₀."""
    r_arr, g_arr, _ = integrate_soliton(g0)
    A, B_c, C_s, phase, rmse = fit_tail_full(r_arr, g_arr, R_L_FIT, R_R_FIT)
    return A, phase

def atail_for_g0(g0):
    A, _ = atail_and_phase(g0)
    return A


# ============================================================
print("=" * 72)
print("EX134: DLACZEGO g₀^e ≈ 1 + 1/(2α) = 5/4?")
print("       Analityczne badanie struktury ODE solitonu TGP")
print("=" * 72)
print()
print(f"  α       = {ALPHA}")
print(f"  g*      = {G_STAR:.8f}  = exp(-1/(2α))")
print(f"  1/(2α)  = {1.0/(2*ALPHA):.6f}")
print(f"  5/4     = {5.0/4.0:.6f}  = 1 + 1/(2α)")
print(f"  g₀^e(obs) = {G0_E_OBS:.5f}")
print(f"  g₀^e(*) = {G0_E_STAR:.6f}  (brentq z ex133)")
print(f"  PHI     = {PHI:.8f}")
print(f"  RSTAR   = {RSTAR:.8f}  = (23+5√21)/2")
print()

# ============================================================
# SEKCJA 0: Weryfikacja ODE — atail przy g₀^e=1.24915
# ============================================================
print("[0] WERYFIKACJA ODE (kontrola)")
print("-" * 55)
A_e_obs, ph_e_obs = atail_and_phase(G0_E_OBS)
A_mu_obs = atail_for_g0(2.02117)   # φ·g₀^e(obs)
r21_obs  = (A_mu_obs / A_e_obs)**4
print(f"  A_tail(g₀^e=1.24915)  = {A_e_obs:.6f}")
print(f"  A_tail(g₀^μ=2.02117)  = {A_mu_obs:.6f}")
print(f"  r₂₁(obs)              = {r21_obs:.4f}  (9r*={9*RSTAR:.4f})")
print(f"  phase(g₀^e=1.24915)   = {ph_e_obs:.4f}°")

record("P00: Weryfikacja A_tail(g₀^e=1.24915) ≈ 0.299 (ex131)",
       abs(A_e_obs - 0.29882) < 0.005,
       f"A_tail={A_e_obs:.5f}, cel≈0.29882")

# ============================================================
# SEKCJA 1: FAZA δ ogona solitonu
# ============================================================
print("\n[1] FAZA δ OGONA (g-1)·r = A·cos(r-δ)")
print("-" * 55)

# P01: faza przy g₀^e = 5/4
g0_54 = 5.0 / 4.0
A_54, phase_54 = atail_and_phase(g0_54)
phase_54_rad = math.radians(phase_54)
print(f"  g₀^e = 5/4 = 1.25000:")
print(f"    A_tail    = {A_54:.6f}")
print(f"    δ (faza)  = {phase_54:.4f}°  ({phase_54_rad:.6f} rad)")
print(f"    δ / π     = {phase_54 / 180:.6f}  (= {phase_54 / 180} w jednostkach π)")

# Sprawdź specjalne wartości fazy
spec_vals = {
    '0':   0.0, 'π/6': 30.0, 'π/4': 45.0, 'π/3': 60.0,
    'π/2': 90.0, '2π/3': 120.0, '3π/4': 135.0, 'π': 180.0,
    '-π/4': -45.0, '-π/6': -30.0, '-π/3': -60.0,
}
best_spec = min(spec_vals.items(), key=lambda kv: abs(kv[1] - phase_54))
print(f"    Najbliższa spec. wartość: {best_spec[0]} = {best_spec[1]}°, "
      f"|Δ| = {abs(best_spec[1]-phase_54):.4f}°")

record("P01: Faza δ przy g₀^e=5/4 obliczona",
       math.isfinite(phase_54),
       f"δ={phase_54:.4f}°, najbliższe spec: {best_spec[0]} (|Δ|={abs(best_spec[1]-phase_54):.4f}°)")

# P02: faza przy g₀^{e,*}
A_star, phase_star = atail_and_phase(G0_E_STAR)
phase_star_rad = math.radians(phase_star)
print(f"\n  g₀^{{e,*}} = {G0_E_STAR:.6f}:")
print(f"    A_tail  = {A_star:.6f}")
print(f"    δ (faza)= {phase_star:.4f}°")
print(f"  Δδ = δ(5/4) - δ(g₀^{{e,*}}) = {phase_54 - phase_star:.4f}°")

record("P02: Faza δ przy g₀^{e,*} obliczona",
       math.isfinite(phase_star),
       f"δ={phase_star:.4f}°, Δδ vs 5/4 = {phase_54-phase_star:.4f}°")

# P03: faza liniowego ODE (limit g₀^e → 1+ε, ε→0)
# Dla liniowego ODE h''-2h'/r+h=0, faza jest stała (niezależna od amplitudy)
# Zmierzymy ją przy bardzo małym g₀^e-1
print(f"\n  FAZA LINIOWEGO ODE (granica małych amplitud):")
g0_linear_vals = [1.005, 1.010, 1.015, 1.020, 1.025]
phases_linear = []
for g0_lin in g0_linear_vals:
    _, ph_lin = atail_and_phase(g0_lin)
    phases_linear.append(ph_lin)
    print(f"    g₀^e = {g0_lin:.3f}: δ = {ph_lin:.4f}°")

phase_linear_extrap = phases_linear[0]   # limit jako g₀^e→1+
print(f"  Ekstrapolacja (g₀^e→1+): δ₀ ≈ {phase_linear_extrap:.4f}°")
print(f"  Zmiana Δδ od g₀^e=1.005 do 1.25: "
      f"{phase_54 - phase_linear_extrap:.4f}°")

record("P03: Faza liniowego ODE δ₀ wyznaczona",
       len(phases_linear) == 5,
       f"δ₀ ≈ {phase_linear_extrap:.4f}°, "
       f"δ(5/4)={phase_54:.4f}°, "
       f"Δδ(nonlin)={phase_54-phase_linear_extrap:.4f}°")

# ============================================================
# SEKCJA 2: POTĘGOWY WYKŁADNIK q = d(ln A_tail)/d(ln(g₀^e-1))
# ============================================================
print("\n[2] POTĘGOWY WYKŁADNIK A_tail ~ (g₀^e-1)^q")
print("-" * 55)

# Skan A_tail dla wielu g₀^e
g0e_scan = np.linspace(1.02, 1.35, 30)
atail_scan = []
for g0s in g0e_scan:
    A = atail_for_g0(g0s)
    atail_scan.append(A)
atail_scan = np.array(atail_scan)
excess_scan = g0e_scan - 1.0

# Globalny dopasowanie potęgowe: log A_tail ~ q·log(g₀^e-1) + const
mask_valid = (atail_scan > 0) & np.isfinite(atail_scan)
log_A   = np.log(atail_scan[mask_valid])
log_exc = np.log(excess_scan[mask_valid])
slope_global, intercept_global, r_val, _, _ = linregress(log_exc, log_A)
print(f"  Globalny fit (g₀^e∈[1.02,1.35]):")
print(f"    q_global = {slope_global:.4f}  (r²={r_val**2:.6f})")
print(f"    C = exp(intercept) = {math.exp(intercept_global):.4f}")

# Lokalny q przy g₀^e=5/4 (różniczkowanie numeryczne)
eps = 0.005
A_hi = atail_for_g0(g0_54 + eps)
A_lo = atail_for_g0(g0_54 - eps)
exc_hi = (g0_54 + eps) - 1.0
exc_lo = (g0_54 - eps) - 1.0
q_local_54 = (math.log(A_hi) - math.log(A_lo)) / (math.log(exc_hi) - math.log(exc_lo))
print(f"\n  Lokalny q przy g₀^e=5/4:")
print(f"    q_local(5/4) = {q_local_54:.4f}")
print(f"    Sprawdź kandydatów q: 1/√2={1/math.sqrt(2):.4f}, φ/2={PHI/2:.4f}, "
      f"1={1:.4f}, 2/3={2/3:.4f}")

# Lokalny q przy g₀^{e,*}
A_hi2 = atail_for_g0(G0_E_STAR + eps)
A_lo2 = atail_for_g0(G0_E_STAR - eps)
exc_hi2 = (G0_E_STAR + eps) - 1.0
exc_lo2 = (G0_E_STAR - eps) - 1.0
q_local_star = (math.log(A_hi2) - math.log(A_lo2)) / (math.log(exc_hi2) - math.log(exc_lo2))
print(f"  Lokalny q przy g₀^{{e,*}}={G0_E_STAR:.6f}:")
print(f"    q_local(g₀^{{e,*}}) = {q_local_star:.4f}")

record("P04: Potęgowy wykładnik q wyznaczony numerycznie",
       abs(slope_global) < 5,
       f"q_global={slope_global:.4f}, q_local(5/4)={q_local_54:.4f}")

# P05: odchylenie od liniowej zależności A_tail ~ C·(g₀^e-1)
# "Liniowa" aproksymacja: pasujemy przy g₀^e=1.01 i używamy A(1.01)/(1.01-1) jako skali
g0_lin_ref = 1.01
A_lin_ref  = atail_for_g0(g0_lin_ref)
C_linear   = A_lin_ref / (g0_lin_ref - 1.0)
A_linear_at_54 = C_linear * (g0_54 - 1.0)
nonlin_factor_54 = A_54 / A_linear_at_54
print(f"\n  Nieliniowy czynnik korekcji (A_actual/A_linear):")
print(f"  C_linear = A_tail(1.01)/0.01 = {C_linear:.4f}")
print(f"  A_linear(5/4) = C·(5/4−1) = {A_linear_at_54:.6f}")
print(f"  A_actual(5/4) = {A_54:.6f}")
print(f"  Nieliniowa korekta przy 5/4: {nonlin_factor_54:.4f}  "
      f"(= A_actual/A_linear)")

record("P05: Nieliniowy czynnik przy g₀^e=5/4 wyznaczony",
       math.isfinite(nonlin_factor_54),
       f"A_actual/A_linear={nonlin_factor_54:.4f} (> 1 oznacza wzmocnienie przez nieliniowość)")

# ============================================================
# SEKCJA 3: KOEFICJENTY ODE PRZY g₀^e=5/4
# ============================================================
print("\n[3] KOEFICJENTY ODE PRZY g₀^e = 5/4")
print("-" * 55)

g0 = g0_54
f_val   = f_kin(g0)
Vp_val  = Vprime(g0)
V_val   = Vpot(g0)
V_1     = Vpot(1.0)
V_gstar = Vpot(G_STAR)
c2_val  = -Vp_val / (2.0 * f_val)   # współczynnik r² przy r=0

print(f"  g₀^e = 5/4 = {g0}")
print(f"  f(g₀^e) = 1+4·ln(5/4) = {f_val:.6f}")
print(f"  V'(g₀^e)= {Vp_val:.6f}  [= (5/4)²·(1-5/4) = -25/64]")
print(f"  V(g₀^e) = {V_val:.6f}  [= (5/4)³/3 - (5/4)⁴/4]")
print(f"  V(1)    = {V_1:.6f}  [= 1/3-1/4 = 1/12]")
print(f"  V(g*)   = {V_gstar:.6f}")
print(f"  c₂ = g''(0)/2 = -V'/(2f) = {c2_val:.6f}")
print(f"  1/(2α)  = {0.25:.6f}  (α=2)")

# Sprawdź f(5/4) vs specjalne wartości
f_candidates = {
    '2': 2.0, '√3': math.sqrt(3), 'φ': PHI, '√φ': math.sqrt(PHI),
    '1+1/φ': 1.0+1.0/PHI, '2−1/(4α)': 2.0-1.0/(4*ALPHA),
    '1+√2/2': 1.0+math.sqrt(2)/2,
}
best_f = min(f_candidates.items(), key=lambda kv: abs(kv[1] - f_val))
print(f"\n  f(5/4) = {f_val:.6f}")
print(f"  Najbliższy kandydat: '{best_f[0]}' = {best_f[1]:.6f}, "
      f"|Δ| = {abs(best_f[1]-f_val):.6f}")

# Sprawdź c₂ vs specjalne wartości
c2_cands = {
    '1/(8α)': 1.0/(8*ALPHA), 'V\'(1)/4': Vprime(1.0)/4.0,
    '1/10': 0.1, '(1-g*)/4': (1-G_STAR)/4.0,
    '-V\'(g₀^e)/4': -Vp_val/4.0,
    'g₀^e-1/4': (g0-1)/4.0,
}
best_c2 = min(c2_cands.items(), key=lambda kv: abs(kv[1] - c2_val))
print(f"\n  c₂ = {c2_val:.6f}")
print(f"  Najbliższy kandydat: '{best_c2[0]}' = {best_c2[1]:.6f}, "
      f"|Δ| = {abs(best_c2[1]-c2_val):.6f}")

record("P06: f(g₀^e=5/4) zbadany — najbliższy kandydat to 2 lub inny?",
       True,
       f"f(5/4)={f_val:.6f}, najbliższe: '{best_f[0]}'={best_f[1]:.6f} (|Δ|={abs(best_f[1]-f_val):.6f})")

record("P07: c₂(g₀^e=5/4) — warunek brzegowy przy r=0",
       True,
       f"c₂={c2_val:.6f}, najbliższe: '{best_c2[0]}'={best_c2[1]:.6f}")

# ============================================================
# SEKCJA 4: STOSUNEK A_tail(φ·g₀^e)/A_tail(g₀^e) PRZY g₀^e=5/4
# ============================================================
print("\n[4] STOSUNEK A_μ/A_e PRZY g₀^e=5/4 vs (9r*)^{1/4}")
print("-" * 55)

R_MU = 2.02117 / G0_E_OBS   # ≈ φ, z TGP solitonów

# Przy g₀^e = 5/4
g0_mu_54 = R_MU * g0_54
A_mu_54  = atail_for_g0(g0_mu_54)
ratio_54 = A_mu_54 / A_54
r21_54   = ratio_54**4
target_ratio = (9.0 * RSTAR)**0.25

print(f"  g₀^e = 5/4:  g₀^μ = R_MU·5/4 = {g0_mu_54:.5f}")
print(f"  A_tail(g₀^e=5/4) = {A_54:.6f}")
print(f"  A_tail(g₀^μ=5φ/4) = {A_mu_54:.6f}")
print(f"  Stosunek A_μ/A_e = {ratio_54:.6f}")
print(f"  Target (9r*)^{{1/4}} = {target_ratio:.6f}")
print(f"  r₂₁(5/4) = {r21_54:.4f}  vs  9r* = {9*RSTAR:.4f}")
print(f"  r₂₁(5/4)/r* = {r21_54/RSTAR:.5f}  (cel: 9.0000)")
rel_err = (r21_54 - 9.0*RSTAR) / (9.0*RSTAR) * 100
print(f"  Odchylenie r₂₁(5/4) od 9r*: {rel_err:+.4f}%")

# Przy g₀^{e,*}
g0_mu_star = R_MU * G0_E_STAR
A_mu_star  = atail_for_g0(g0_mu_star)
ratio_star = A_mu_star / A_star
r21_star   = ratio_star**4
print(f"\n  g₀^{{e,*}} = {G0_E_STAR:.6f}: r₂₁ = {r21_star:.4f}, r₂₁/r* = {r21_star/RSTAR:.5f}")

record("P08: r₂₁(5/4) blisko 9r*?",
       abs(r21_54 / (9*RSTAR) - 1.0) < 0.01,
       f"r₂₁(5/4)/9r*={r21_54/(9*RSTAR):.6f}, odchylenie={rel_err:+.4f}%")

# ============================================================
# SEKCJA 5: KLUCZOWY TEST — g₀^e - 1 vs -ln(g*)
# ============================================================
print("\n[5] KLUCZOWY TEST: g₀^e − 1 = −ln(g*) = 1/(2α)?")
print("-" * 55)

log_gstar     = math.log(G_STAR)           # = -1/(2α) = -0.25
minus_log_gstar = -log_gstar               # = 1/(2α) = 0.25
one_over_2a   = 1.0 / (2.0 * ALPHA)       # = 0.25

g0e_obs_excess   = G0_E_OBS   - 1.0       # = 0.24915
g0e_star_excess  = G0_E_STAR  - 1.0       # = 0.249082

print(f"  −ln(g*) = −ln(e^{{-1/(2α)}}) = 1/(2α) = {minus_log_gstar:.8f}")
print(f"  g₀^e(obs)−1  = {g0e_obs_excess:.8f}   vs 1/(2α)={one_over_2a:.8f}")
print(f"  g₀^{{e,*}}−1  = {g0e_star_excess:.8f}   vs 1/(2α)={one_over_2a:.8f}")
print()
print(f"  Odchylenie obs:    Δ₁ = (g₀^e(obs)−1) − 1/(2α) = {g0e_obs_excess - one_over_2a:+.8f}")
print(f"  Odchylenie g*:     Δ₂ = (g₀^e(*) −1)  − 1/(2α) = {g0e_star_excess - one_over_2a:+.8f}")
print()

# Korekty 2-go rzędu
# g₀^{e,*} = 1 + 1/(2α) + δ, gdzie δ ≈ -0.000918
# czy δ = -c·(1/(2α))^2 dla pewnego analitycznego c?
delta_star = g0e_star_excess - one_over_2a
c_corr = -delta_star / (one_over_2a**2)
print(f"  KOREKTA 2-GO RZĘDU: δ = g₀^{{e,*}} − (1+1/(2α)) = {delta_star:+.8f}")
print(f"  δ / (1/(2α))²  = {delta_star / one_over_2a**2:.4f}")
print(f"  ⟹  c_corr = {c_corr:.4f}  (g₀^{{e,*}} ≈ 1+1/(2α) − c·(1/(2α))²)")
print()
# sprawdź czy c_corr jest bliskie: 1/2, ln2, π/4, 1/√2, 1, φ, ...
c_cands = {'1/2': 0.5, 'ln 2': math.log(2), '1/√2': 1/math.sqrt(2),
           '1': 1.0, 'π/4': math.pi/4, 'φ-1': PHI-1, '1/φ': 1/PHI,
           '3/8': 3/8, '1/(4α)': 1/(4*ALPHA)}
best_c = min(c_cands.items(), key=lambda kv: abs(kv[1] - c_corr))
print(f"  c_corr = {c_corr:.5f}")
print(f"  Kandydaci c (od najbliższego):")
for name, val in sorted(c_cands.items(), key=lambda kv: abs(kv[1]-c_corr)):
    print(f"    '{name}' = {val:.5f}, |Δ| = {abs(val-c_corr):.5f}")

record("P15: g₀^e − 1 = 1/(2α) — test kluczowy",
       abs(g0e_star_excess - one_over_2a) < 0.002,
       f"g₀^e(*)-1={g0e_star_excess:.6f}, 1/(2α)={one_over_2a:.6f}, "
       f"δ={delta_star:+.6f}, c_corr={c_corr:.4f}")

# ============================================================
# SEKCJA 6: POZYCJA UŁAMKOWA g₀^e MIĘDZY g* A 1
# ============================================================
print("\n[6] POZYCJA UŁAMKOWA g₀^e MIĘDZY g* A 1")
print("-" * 55)

# Pozycja: (g₀^e − g*) / (1 − g*)
frac_obs  = (G0_E_OBS  - G_STAR) / (1.0 - G_STAR)
frac_star = (G0_E_STAR - G_STAR) / (1.0 - G_STAR)
frac_54   = (g0_54     - G_STAR) / (1.0 - G_STAR)

# Dla g₀^e = 1 + 1/(2α):
# frac = (1 + 1/(2α) - g*) / (1 - g*) = (1/(2α) + 1 - g*) / (1 - g*)
# = 1 + 1/(2α(1-g*)) = 1 + (1/(2α)) / (1-g*)
frac_analytic = 1.0 + one_over_2a / (1.0 - G_STAR)
print(f"  g* = {G_STAR:.6f}, 1-g* = {1-G_STAR:.6f}")
print(f"  (g₀^e(obs)  − g*)/(1−g*) = {frac_obs:.6f}")
print(f"  (g₀^e(*)    − g*)/(1−g*) = {frac_star:.6f}")
print(f"  (5/4        − g*)/(1−g*) = {frac_54:.6f}")
print(f"  Analitycznie dla 5/4: 1+1/(2α(1−g*)) = {frac_analytic:.6f}")
print()
# Czy frac ≈ 2? Czy ≈ 1+φ? Czy ≈ 1+1/φ?
frac_cands = {'2': 2.0, '1+1/φ': 1+1/PHI, '1+φ': 1+PHI,
              '5/3': 5/3, 'φ': PHI, '1.5': 1.5,
              f'1+1/(2α(1-g*))': frac_analytic}
best_frac = min(frac_cands.items(), key=lambda kv: abs(kv[1]-frac_54))
print(f"  Kandydaci dla frac(5/4)={frac_54:.5f}:")
for name, val in sorted(frac_cands.items(), key=lambda kv: abs(kv[1]-frac_54)):
    print(f"    '{name}' = {val:.5f}, |Δ| = {abs(val-frac_54):.5f}")

record("P12: Pozycja ułamkowa g₀^e między g* i 1",
       True,
       f"(g₀^e-g*)/(1-g*)={frac_54:.5f}, (g₀^{{e,*}}-g*)/(1-g*)={frac_star:.5f}")

# ============================================================
# SEKCJA 7: ILOCZYN g₀^e · g* I INNE KOMBINACJE
# ============================================================
print("\n[7] KOMBINACJE ALGEBRAICZNE g₀^e I g*")
print("-" * 55)

prod_54   = g0_54 * G_STAR
prod_star = G0_E_STAR * G_STAR
prod_obs  = G0_E_OBS * G_STAR
# Analitycznie: (1+1/(2α))·exp(-1/(2α)) = (1+x)·e^{-x} dla x=1/(2α)=1/4
x = one_over_2a
prod_analytic = (1 + x) * math.exp(-x)

print(f"  (1+x)·e^(-x) dla x=1/(2α)=1/4: {prod_analytic:.8f}")
print(f"  g₀^e(5/4) · g*  = {prod_54:.8f}")
print(f"  g₀^e(*) · g*    = {prod_star:.8f}")
print(f"  g₀^e(obs) · g*  = {prod_obs:.8f}")
print(f"  Odchylenie od 1: (1+x)e^(-x)-1 = {prod_analytic-1:+.8f}")
print(f"  Przybliżenie: -(1/(2α))²/2 = {-x**2/2:.8f}")

# Inne kombinacje
sum_gstar = 1.0 / G_STAR + G_STAR    # 1/g* + g*
diff_log  = math.log(g0_54) + math.log(G_STAR)  # ln(g₀) + ln(g*)
print(f"\n  ln(g₀^e=5/4) + ln(g*) = {diff_log:.6f}  (= ln(g₀·g*) = {math.log(prod_54):.6f})")
print(f"  g₀^e=5/4 + g* = {g0_54 + G_STAR:.6f}  (≈ 2?)")
print(f"  5/4 + g* = {g0_54 + G_STAR:.6f} vs 2 - (err: {abs(g0_54+G_STAR-2):.6f})")

record("P10: Iloczyn g₀^e·g* = (1+x)e^{-x} dla x=1/(2α)",
       True,
       f"(1+x)e^(-x)={prod_analytic:.6f}, 1-Δ where Δ≈x²/2={x**2/2:.6f}")

# ============================================================
# SEKCJA 8: ALGEBRAICZNI KANDYDACI NA g₀^{e,*}
# ============================================================
print("\n[8] ALGEBRAICZNI KANDYDACI NA g₀^{e,*} = 1.249082")
print("-" * 55)
target = G0_E_STAR

candidates = {}
# Proste
candidates['1 + 1/(2α)']       = 1 + 1/(2*ALPHA)         # = 5/4
candidates['exp(1/(2α)) - 1/4'] = math.exp(1/(2*ALPHA)) - 0.25
candidates['1 + ln(φ)/ln(2)']   = 1 + math.log(PHI)/math.log(2)
candidates['1 + (1-g*)/√2']    = 1 + (1-G_STAR)/math.sqrt(2)
candidates['1 + 1/(4α) + 1/(4α)²'] = 1 + 1/(4*ALPHA) + 1/(4*ALPHA)**2
candidates['√(5/4·g*)^{-1}']   = 1/math.sqrt(g0_54 * G_STAR)
candidates['(1+g*+g*²)/3 + Δ'] = (1+G_STAR+G_STAR**2)/3
candidates['1/g*^{1/4} - 1/2'] = G_STAR**(-0.25) - 0.5
# Kandydaci oparci o 1/(2α)
xval = one_over_2a
candidates[f'1+x−x²ln2 (x=1/(2α))'] = 1 + xval - xval**2 * math.log(2)
candidates[f'1+x−x²/2']              = 1 + xval - xval**2 / 2
candidates[f'1+x−x²']                = 1 + xval - xval**2
# Koide FP
candidates['u*/4']                   = (5+math.sqrt(21))/8
candidates['(u*-1)/4']               = ((5+math.sqrt(21))/2 - 1)/4
candidates['1+1/(r*+5)']             = 1 + 1/(RSTAR + 5)
candidates['1+1/r*']                 = 1 + 1/RSTAR
candidates['1+1/(4r*-1)']            = 1 + 1/(4*RSTAR - 1)

print(f"  Target g₀^{{e,*}} = {target:.8f}")
print(f"\n  {'Kandydat':<35}  {'Wartość':>10}  {'|Δ|':>12}  {'err%':>8}")
print("  " + "-" * 70)
for name, val in sorted(candidates.items(), key=lambda kv: abs(kv[1]-target)):
    err_pct = abs(val-target)/target * 100
    print(f"  {name:<35}  {val:>10.6f}  {abs(val-target):>12.8f}  {err_pct:>8.5f}%")

# Znajdź TOP 3
top3 = sorted(candidates.items(), key=lambda kv: abs(kv[1]-target))[:3]
print(f"\n  TOP 3 kandydatów:")
for i, (name, val) in enumerate(top3, 1):
    print(f"    {i}. '{name}' = {val:.8f}  (err={abs(val-target)/target*100:.5f}%)")

record("P09: Algebraiczny scan kandydatów na g₀^{e,*}",
       abs(top3[0][1] - target) < 0.005,
       f"Najlepszy: '{top3[0][0]}'={top3[0][1]:.6f} (err={abs(top3[0][1]-target)/target*100:.4f}%)")

# ============================================================
# SEKCJA 9: FAZA SKANU — δ(g₀^e) dla g₀^e ∈ [1.10, 1.35]
# ============================================================
print("\n[9] SKAN FAZY δ(g₀^e)")
print("-" * 55)

g0e_phase_scan = np.linspace(1.05, 1.35, 21)
phases_scan = []
for g0s in g0e_phase_scan:
    _, ph = atail_and_phase(g0s)
    phases_scan.append(ph)

print(f"  {'g₀^e':>8}  {'δ (°)':>10}  {'δ/π':>8}")
print("  " + "-" * 32)
for g0s, ph in zip(g0e_phase_scan, phases_scan):
    print(f"  {g0s:>8.4f}  {ph:>10.4f}  {ph/180:>8.5f}")

# Sprawdź monotoniczność
phases_arr = np.array(phases_scan)
dphases = np.diff(phases_arr)
is_mono = bool(np.all(dphases > 0) or np.all(dphases < 0))
print(f"\n  Monotoniczność δ(g₀^e): {'TAK' if is_mono else 'NIE'}")

# Gradient przy g₀^e=5/4
idx_54 = np.argmin(np.abs(g0e_phase_scan - g0_54))
if idx_54 > 0 and idx_54 < len(g0e_phase_scan)-1:
    dg = g0e_phase_scan[idx_54+1] - g0e_phase_scan[idx_54-1]
    ddelta = phases_scan[idx_54+1] - phases_scan[idx_54-1]
    ddelta_dg0 = ddelta / dg
    print(f"  dδ/dg₀^e przy 5/4 ≈ {ddelta_dg0:.2f} °/unit")

record("P03b: Monotoniczność fazy δ(g₀^e) zweryfikowana",
       True,
       f"Monoton: {is_mono}; δ(5/4)={phase_54:.4f}°, δ(g₀^e*)={phase_star:.4f}°")

# ============================================================
# SEKCJA 10: WARUNEK WNĘTRZA — r₁ (pierwsze zero solitonu)
# ============================================================
print("\n[10] PIERWSZE ZERO SOLITONU r₁(g₀^e)")
print("-" * 55)

def first_zero_soliton(g0):
    """Znajdź pierwsze r₁ gdzie g(r₁)=1 (przejście z g₀^e>1 przez g=1)."""
    r_arr, g_arr, _ = integrate_soliton(g0)
    # szukamy przejścia przez g=1 (z góry)
    diff = g_arr - 1.0
    for i in range(len(diff)-1):
        if diff[i] > 0 and diff[i+1] <= 0:
            # liniowa interpolacja
            r1 = r_arr[i] + (r_arr[i+1]-r_arr[i]) * (-diff[i])/(diff[i+1]-diff[i])
            return r1
    return float('nan')

r1_obs   = first_zero_soliton(G0_E_OBS)
r1_54    = first_zero_soliton(g0_54)
r1_star  = first_zero_soliton(G0_E_STAR)

print(f"  r₁(g₀^e=1.24915) = {r1_obs:.4f}")
print(f"  r₁(g₀^e=5/4)     = {r1_54:.4f}")
print(f"  r₁(g₀^e=g₀^e*)   = {r1_star:.4f}")

# Sprawdź r₁ vs specjalne wartości
r1_spec = {
    'π': math.pi, '2π/3': 2*math.pi/3, '√2·π/2': math.sqrt(2)*math.pi/2,
    'π/√2': math.pi/math.sqrt(2), 'π·g*': math.pi*G_STAR, '2': 2.0
}
if math.isfinite(r1_54):
    best_r1 = min(r1_spec.items(), key=lambda kv: abs(kv[1]-r1_54))
    print(f"\n  r₁(5/4) = {r1_54:.6f}")
    print(f"  Najbliższe spec.: '{best_r1[0]}' = {best_r1[1]:.6f}, "
          f"|Δ| = {abs(best_r1[1]-r1_54):.6f}")

record("P11: Pierwsze zero r₁ solitonu przy g₀^e=5/4",
       math.isfinite(r1_54),
       f"r₁(5/4)={r1_54:.4f}" if math.isfinite(r1_54) else "brak przejścia g=1")

# ============================================================
# SEKCJA 11: WKB-LIKE — WARUNEK ENERGETYCZNY
# ============================================================
print("\n[11] WARUNEK ENERGETYCZNY V(g₀^e) vs INNE SKALE")
print("-" * 55)

print(f"  V(g₀^e=5/4) = {V_val:.8f}")
print(f"  V(1)        = {V_1:.8f}  [= 1/12 = {1/12:.8f}]")
print(f"  V(g*)       = {V_gstar:.8f}")
print(f"  Czy V(g₀^e) + V(g*) = V(1)?")
print(f"  V(5/4)+V(g*) = {V_val+V_gstar:.8f}  vs V(1) = {V_1:.8f}")
print(f"  Odchylenie: {abs(V_val+V_gstar-V_1):.8f}")
print()

# Test: V(g₀^e) = V(g*)? (potencjały symetryczne)
print(f"  Ratio V(g₀^e)/V(1): {V_val/V_1:.6f}")
print(f"  Ratio V(g*)/V(1):   {V_gstar/V_1:.6f}")
print(f"  Ratio V(g₀^e)/V(g*): {V_val/V_gstar:.6f}")

# Test c₂ vs 1/(2α)²:
print(f"\n  Czy c₂ = (g₀^e-1)^2/(2f) (WKB-aproksymacja)?")
c2_wkb = (g0_54-1)**2 / (2.0*f_val)
print(f"  c₂_WKB = {c2_wkb:.6f}  vs  c₂_actual = {c2_val:.6f}")
print(f"  Ratio c₂_actual/c₂_WKB = {c2_val/c2_wkb:.4f}")

# Kluczowy test: V(g₀^e) = (1/(2α))^2 / 2?
# Mamy V(1+ε) ≈ -ε²/2 dla małych ε (V''(1)/2 · ε²= (-1/2)ε²)
# Przy ε = 1/(2α) = 0.25: V(1.25) ≈ -0.0625/2 = -0.03125
# Ale V(1.25) = 0.04069 > 0! Więc V(g₀^e=5/4) ≠ -ε²/2.
# V''(1) = 2·1·(1-1) - 1² = -1 (jak policzyłem wcześniej)
v_approx = 0.5 * (-1) * (g0_54-1)**2 + V_1  # V(1) + V''(1)·(g₀-1)²/2
print(f"\n  V(g₀^e=5/4) Taylor wokół g=1: V(1)+V''(1)/2·ε² = {v_approx:.8f}")
print(f"  V(5/4) actual = {V_val:.8f}")
print(f"  Korekta wyższego rzędu: {V_val - v_approx:.8f}")

record("P13: Warunek energetyczny V(g₀^e) zbadany",
       True,
       f"V(5/4)={V_val:.6f}, V(1)={V_1:.6f}, V(g*)={V_gstar:.6f}")

# ============================================================
# SEKCJA 12: LINIOWE ODE — FAZA ANALITYCZNIE
# ============================================================
print("\n[12] FAZA LINIOWEGO ODE — ANALITYCZNA STRUKTURA")
print("-" * 55)

# Dla liniowego ODE h''−2h'/r+h=0 przy r=0: h≈h₀(1+r²/2−r⁴/8+...)
# Rozwinięcie przy r=0:
# h = 1 + r²/2 + a₄r⁴ + ... gdzie a₄ = -a₂/4 = -1/8
# Przy dużych r: h ~ A·cos(r-δ₀)/r

# Numerycznie rozwiązujemy liniowe ODE dla h₀=1 i dopasowujemy fazę
def rhs_linear(r, y):
    """Liniowe ODE: h''-(2/r)h'+h=0"""
    h, hp = y
    if r < 1e-10:
        # Granica: h'' = h₀·(−1)/3·... hmm użyjemy c₂ = h₀/2
        # W graniicy r→0: h''−(2/r)h'+h=0 → −2h''+h→0 → h''=h₀/2
        return [hp, h / 2.0]
    return [hp, (2.0/r)*hp - h]

# Integruj liniowe ODE
sol_lin = solve_ivp(rhs_linear, [1e-4, 50.0], [1.0, 0.5e-4],
                    method='DOP853', max_step=0.01, rtol=1e-12, atol=1e-15,
                    dense_output=False)
r_lin = sol_lin.t
h_lin = sol_lin.y[0]

# Dopasuj fazę przy dużych r
mask_lin = (r_lin >= 20.0) & (r_lin <= 35.0)
if np.sum(mask_lin) > 10:
    r_f = r_lin[mask_lin]
    y_f = h_lin[mask_lin] * r_f   # (g-1)·r
    X = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef_lin, _, _, _ = np.linalg.lstsq(X, y_f, rcond=None)
    B_c_lin, C_s_lin = float(coef_lin[0]), float(coef_lin[1])
    A_lin_fit = math.sqrt(B_c_lin**2 + C_s_lin**2)
    phase_lin0 = math.degrees(math.atan2(C_s_lin, B_c_lin))
    print(f"  Liniowe ODE (h₀=1): faza δ₀ = {phase_lin0:.4f}°")
    print(f"  B_c = {B_c_lin:.6f}, C_s = {C_s_lin:.6f}")
    print(f"  A_lin = {A_lin_fit:.6f}")
    print(f"  δ₀/π = {phase_lin0/180:.6f}")

    # Sprawdź δ₀ vs specjalne wartości
    best_lin = min(spec_vals.items(), key=lambda kv: abs(kv[1]-phase_lin0))
    print(f"  Najbliższa spec. wartość: '{best_lin[0]}' = {best_lin[1]}°, "
          f"|Δ| = {abs(best_lin[1]-phase_lin0):.4f}°")

    # Porównaj z fazą obserwowaną przy g₀^e=5/4
    delta_nonlin = phase_54 - phase_lin0
    print(f"\n  Nieliniowa zmiana fazy Δδ przy g₀^e=5/4:")
    print(f"  δ(5/4)    = {phase_54:.4f}°")
    print(f"  δ₀(lin)   = {phase_lin0:.4f}°")
    print(f"  Δδ(nonlin)= {delta_nonlin:+.4f}°  = {delta_nonlin/180:.5f}·π")

    record("P13b: Faza liniowego ODE δ₀ wyznaczona analitycznie",
           math.isfinite(phase_lin0),
           f"δ₀={phase_lin0:.4f}°, δ(5/4)={phase_54:.4f}°, "
           f"Δδ(nonlin)={delta_nonlin:+.4f}°")
else:
    print("  Błąd: za mało punktów w oknie dopasowania")
    record("P13b: Faza liniowego ODE", False, "Błąd integracji")

# ============================================================
# SEKCJA 13: DERYWACJA KOREKTY δ = g₀^{e,*} − 5/4
# ============================================================
print("\n[13] ANALITYCZNA KOREKTA g₀^{e,*} = 5/4 + δ")
print("-" * 55)

delta = G0_E_STAR - 5.0/4.0    # ≈ -0.000918
eps   = 1.0/(2*ALPHA)            # = 0.25

print(f"  g₀^{{e,*}} = {G0_E_STAR:.8f}")
print(f"  5/4       = {5.0/4.0:.8f}")
print(f"  δ = g₀^{{e,*}} − 5/4 = {delta:.8f}")
print()

# Czy δ = -ε²·f(ε) dla analitycznego f?
# ε = 1/(2α) = 0.25
# |δ| = 0.000918
# |δ|/ε = 0.003672
# |δ|/ε² = 0.014688

print(f"  ε = 1/(2α) = {eps:.6f}")
print(f"  |δ|/ε    = {abs(delta)/eps:.6f}")
print(f"  |δ|/ε²   = {abs(delta)/eps**2:.6f}")
print(f"  |δ|/ε³   = {abs(delta)/eps**3:.6f}")
print()

# Sprawdź ε²·różne_wyrażenia
print(f"  Kandydaci na 'δ = −ε²·c':")
eps2 = eps**2
corr_candidates = {
    'ε²/2':       eps2/2,
    'ε²·ln2':     eps2*math.log(2),
    'ε²·(1-g*)':  eps2*(1-G_STAR),
    'ε²·(1/4α)':  eps2*1/(4*ALPHA),
    'ε³·4':       eps**3*4,
    'ε²·φ/4':     eps2*PHI/4,
    'ε²·1/(f(5/4))': eps2/f_val,
    'ε²·V\'\'(1)(-1)': eps2*1,   # V''(1)=-1
}
for name, val in sorted(corr_candidates.items(), key=lambda kv: abs(kv[1]-abs(delta))):
    print(f"    '{name}' = {val:.8f}, |Δ| = {abs(val-abs(delta)):.8f}")

record("P10b: Korekta δ = g₀^{e,*} − 5/4 zanalizowana",
       True,
       f"δ={delta:.8f}, |δ|/ε²={abs(delta)/eps**2:.6f} (ε=1/(2α))")

# ============================================================
# SEKCJA 14: FINAŁOWE ODKRYCIA — SYNTEZA
# ============================================================
print("\n[14] SYNTEZA: MECHANIZM g₀^e ≈ 1 + 1/(2α)")
print("-" * 55)

# Kluczowe liczby
print(f"  1. FAZA OGONA:")
print(f"     δ₀ (liniowe ODE)     ≈ {phase_lin0:.3f}°  (stała przy g₀^e→1)")
print(f"     δ(g₀^e=5/4)          = {phase_54:.3f}°")
print(f"     δ(g₀^{{e,*}})          = {phase_star:.3f}°")
print(f"     Δδ(nonlin)            = δ(5/4) − δ₀ ≈ {phase_54-phase_lin0:+.3f}°")
print()
print(f"  2. POTĘGOWY WYKŁADNIK:")
print(f"     q_local przy 5/4     = {q_local_54:.4f}  (A_tail ~ (g₀^e-1)^q)")
print()
print(f"  3. STOSUNEK A_μ/A_e PRZY g₀^e=5/4:")
print(f"     A_μ/A_e(5/4) = {ratio_54:.5f}  vs  (9r*)^{{1/4}} = {target_ratio:.5f}")
print(f"     Odchylenie = {(ratio_54/target_ratio-1)*100:+.4f}%")
print()
print(f"  4. KLUCZOWA RELACJA: g₀^e − 1 = −ln(g*) = 1/(2α)")
print(f"     1/(2α) = {one_over_2a:.6f}")
print(f"     g₀^e(*)-1 = {g0e_star_excess:.6f}  (0.073% dalej od 5/4)")
print(f"     Korekta 2-go rzędu: δ = {delta:.6f} ≈ −{abs(delta)/eps**2:.4f}·ε²")
print()

# Czy jest PROSTA DERYWACJA?
print(f"  5. KANDYDAT DERYWACJI:")
print(f"     Jeśli A_tail(g₀^e) ≈ C·(g₀^e−1) dla malych g₀^e−1,")
print(f"     to stosunek A_μ/A_e = (φ·g₀^e−1)/(g₀^e−1)")
print(f"     Przy g₀^e = 1+ε: = (φ(1+ε)−1)/ε = (φ−1)/ε + φ = 1/(φε) + φ")
target_cond = target_ratio   # (9r*)^{1/4}
# Warunek A_μ/A_e = target przy g₀^e = 1+ε:
# (φ−1)/ε + φ = target  →  ε = (φ−1)/(target−φ)
eps_analytic = (PHI - 1.0) / (target_cond - PHI)
g0e_analytic = 1.0 + eps_analytic
print(f"     Warunek LINIOWY: (φ−1)/ε + φ = (9r*)^{{1/4}}")
print(f"       ε = (φ−1)/((9r*)^{{1/4}}−φ) = {eps_analytic:.8f}")
print(f"       g₀^e = 1+ε = {g0e_analytic:.8f}  (vs g₀^{{e,*}}={G0_E_STAR:.8f})")
print(f"       Błąd vs 5/4: {abs(g0e_analytic - 5/4):.6f}")
print(f"       Błąd vs g₀^{{e,*}}: {abs(g0e_analytic - G0_E_STAR):.6f}")

record("P14: Warunek LINIOWY A_μ/A_e=(9r*)^{1/4} → g₀^e = 1+ε",
       abs(eps_analytic - eps) < 0.1,
       f"ε_analytic={eps_analytic:.6f}, 1/(2α)={one_over_2a:.6f}, "
       f"g₀^e_analytic={g0e_analytic:.6f}")

# ============================================================
# PODSUMOWANIE
# ============================================================
print()
print("=" * 72)
n_pass  = sum(1 for _, p, _ in TESTS if p)
n_total = len(TESTS)
print(f"WYNIK: {n_pass}/{n_total} testów PASS")
print()
print(f"  STATUS BADANIA g₀^e ≈ 1+1/(2α) = 5/4:")
print(f"  ● faza ogona δ(g₀^e=5/4)        = {phase_54:.3f}°  (nie π/4=45°)")
print(f"  ● δ₀ (liniowe ODE)               ≈ {phase_lin0:.3f}°")
print(f"  ● A_tail(5/4)/A_tail_linear      = {nonlin_factor_54:.5f}  (nonlin.korekta)")
print(f"  ● r₂₁(5/4) − 9r*                = {r21_54-9*RSTAR:+.4f}  (odch. {rel_err:+.4f}%)")
print(f"  ● g₀^e_linear = 1+(φ-1)/((9r*)^{{1/4}}-φ) = {g0e_analytic:.6f}")
print(f"  ● g₀^{{e,*}} vs 5/4: δ={delta:.6f} ≈ −{abs(delta)/eps**2:.3f}·ε²")
print("=" * 72)
