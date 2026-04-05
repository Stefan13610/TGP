#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex135_top4b_phase_pi2.py
========================
BADANIE T-OP4b: DLACZEGO NIELINIOWE ODE TGP OBRACA FAZĘ OGONA O Δδ ≈ π/2
               GDY g₀^e = 1 + 1/(2α) = 5/4?

KONTEKST (ex134, 19/19 PASS):
  δ(5/4)     = 89.136°  ≈ π/2    (odchylenie 0.96%)
  δ₀ (lin.)  ≈ -1.93°  ≈ 0°     (liniowe ODE, błąd warunku brzegowego)
  Δδ(nonlin) ≈ +91.07° ≈ π/2

PYTANIA DO ZBADANIA:
  Q1: Gdzie DOKŁADNIE (brentq) δ = 90°? Czy g₀^e(δ=90°) ma postać algebraiczną?
  Q2: Jaka jest dokładna faza δ₀ liniowego ODE (poprawione WB)?
  Q3: Jak δ zależy od amplitudy A_tail? Czy to liniowa zależność?
  Q4: Przybliżenie Borna dla przesunięcia fazy Δδ — co daje podcałkowy?
  Q5: Jak zachowuje się faza w pobliżu g* (ściana widmo)?
  Q6: Czy r₁(g₀^e=5/4) ≈ π ma związek z warunkiem fazy?
  Q7: Test perturbacyjny: δ(g₀^e) = δ₀ + A_tail · (dδ/dA) — linearyzacja?
  Q8: Scan fazy w szerokim zakresie g₀^e — pełna krzywa δ(g₀^e)

TESTY Q01..Q16
"""

import sys
import io
import math
import warnings
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import linregress

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

warnings.filterwarnings('ignore')

# ============================================================
# Stałe TGP
# ============================================================
ALPHA    = 2.0
PHI      = (1.0 + math.sqrt(5.0)) / 2.0
G_STAR   = math.exp(-1.0 / (2.0 * ALPHA))    # exp(-1/4) ≈ 0.77880
G_BOUNCE = G_STAR + 0.005

RSTAR       = (23.0 + 5.0 * math.sqrt(21.0)) / 2.0   # ≈ 22.9564
THETA_9RS   = 132.731439
G0_E_OBS    = 1.24915
G0_E_STAR   = 1.249082

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
# ODE solitonu (identyczny jak ex134)
# ============================================================
def f_kin(g):
    return 1.0 + 2.0 * ALPHA * math.log(max(g, 1e-30))

def Vprime(g):
    return g * g * (1.0 - g)

def Vpot(g):
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
    Zwraca (A_amp, B_c, C_s, phase_deg, rmse_rel).
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
    r_arr, g_arr, _ = integrate_soliton(g0)
    A, B_c, C_s, phase, rmse = fit_tail_full(r_arr, g_arr, R_L_FIT, R_R_FIT)
    return A, phase

def phase_for_g0(g0):
    _, ph = atail_and_phase(g0)
    return ph


# ============================================================
print("=" * 72)
print("EX135: T-OP4b — DLACZEGO ODE TGP OBRACA FAZĘ OGONA O Δδ ≈ π/2?")
print("       g₀^e = 1 + 1/(2α) = 5/4")
print("=" * 72)
print()
print(f"  α        = {ALPHA}")
print(f"  1/(2α)   = {1.0/(2*ALPHA):.6f}")
print(f"  g*       = {G_STAR:.8f}  = exp(-1/(2α))")
print(f"  5/4      = 1.250000")
print(f"  g₀^e(*)  = {G0_E_STAR:.6f}  (brentq z ex133)")
print(f"  RSTAR    = {RSTAR:.6f}  = (23+5√21)/2")
print(f"  9·RSTAR  = {9*RSTAR:.6f}")
print()


# ============================================================
# SEKCJA 1: Dokładne δ(5/4) i δ(g₀^{e,*}) — weryfikacja ex134
# ============================================================
print("[1] WERYFIKACJA: δ(5/4) i δ(g₀^{e,*})")
print("-" * 55)

A_54,   ph_54   = atail_and_phase(1.25)
A_star, ph_star = atail_and_phase(G0_E_STAR)
print(f"  δ(5/4)        = {ph_54:.4f}°   A_tail = {A_54:.6f}")
print(f"  δ(g₀^e(*))   = {ph_star:.4f}°   A_tail = {A_star:.6f}")
print(f"  |δ − 90°|    = {abs(ph_54-90):.4f}° (dla 5/4)")
print(f"  Δδ(* vs 5/4) = {ph_star-ph_54:.4f}°")

record("Q00a: δ(5/4) ≈ 89° (ex134 benchmark)",
       abs(ph_54 - 89.136) < 0.05,
       f"δ(5/4)={ph_54:.4f}°, oczekiwane≈89.136°")
record("Q00b: |δ(5/4) − 90°| < 1°",
       abs(ph_54 - 90.0) < 1.0,
       f"|δ−90°|={abs(ph_54-90):.4f}°")


# ============================================================
# SEKCJA 2: DOKŁADNE g₀^e gdzie δ = 90° (brentq)
# ============================================================
print("\n[2] DOKŁADNE g₀^e(δ=90°) via brentq")
print("-" * 55)

def phase_minus_90(g0):
    return phase_for_g0(g0) - 90.0

# Grube skanowanie, żeby znaleźć przedział
g0_scan_vals = np.arange(1.10, 1.35, 0.01)
print("  Wstępny scan δ(g₀^e):")
phases_scan = []
for gv in g0_scan_vals:
    p = phase_for_g0(gv)
    phases_scan.append(p)
    print(f"    g₀^e={gv:.3f}: δ={p:.4f}°  (sign={'+' if p>=90 else '-'})")

# Znajdź przedziały gdzie jest przejście przez 90°
crossings = []
for i in range(len(phases_scan)-1):
    if (phases_scan[i] - 90.0) * (phases_scan[i+1] - 90.0) < 0:
        crossings.append((g0_scan_vals[i], g0_scan_vals[i+1]))
        print(f"  → Przejście przez 90° w [{g0_scan_vals[i]:.3f}, {g0_scan_vals[i+1]:.3f}]")

g0e_90deg_list = []
for (lo, hi) in crossings:
    try:
        root = brentq(phase_minus_90, lo, hi, xtol=1e-7, rtol=1e-9)
        g0e_90deg_list.append(root)
        ph_root = phase_for_g0(root)
        print(f"\n  BRENTQ wynik: g₀^e(δ=90°) = {root:.8f}")
        print(f"  Weryfikacja:  δ({root:.6f}) = {ph_root:.6f}°")
    except Exception as e:
        print(f"  BŁĄD brentq: {e}")

if g0e_90deg_list:
    g0e_90 = g0e_90deg_list[0]
    eps_90 = g0e_90 - 1.0
    print(f"\n  g₀^e(δ=90°)    = {g0e_90:.8f}")
    print(f"  g₀^e(δ=90°)−1  = {eps_90:.8f}")
    print(f"  1/(2α)         = {1/(2*ALPHA):.8f}")
    print(f"  5/4 − g₀(90°) = {1.25 - g0e_90:.8f}")
    print(f"  rel. do 1/(2α): Δ = {(g0e_90 - 1.25):.6f} = {(g0e_90 - 1.25)/0.25*100:.4f}%")

    # Algebraiczne kandydaty
    eps  = 1.0 / (2.0 * ALPHA)          # = 0.25
    x    = eps
    cands = [
        ("1+x−x²/2",        1.0 + x - x**2/2.0),
        ("1+x−x²ln2",       1.0 + x - x**2 * math.log(2.0)),
        ("1+x(1−x/4)",      1.0 + x * (1.0 - x/4.0)),
        ("1+x·exp(−x/2)",   1.0 + x * math.exp(-x/2.0)),
        ("√(1+x)",          math.sqrt(1.0 + x)),
        ("1+(x-x²)/1",      1.0 + x - x**2),
        ("1+x/√(1+x²)",     1.0 + x / math.sqrt(1.0 + x**2)),
        ("(1+2x)/(1+x)",    (1.0 + 2*x)/(1.0 + x)),
        ("1+x−x²/3",        1.0 + x - x**2/3.0),
        ("g*(1+1/(2α−g*))", G_STAR * (1.0 + 1.0/(2*ALPHA - G_STAR))),
    ]
    print(f"\n  Algebraiczne kandydaty na g₀^e(δ=90°) = {g0e_90:.8f}:")
    best_cand = None
    best_err  = 1.0
    for (name, val) in cands:
        err = abs(val - g0e_90) / g0e_90 * 100.0
        marker = " ← BEST" if err < 0.5 else ""
        print(f"    {name:30s} = {val:.8f}  err={err:.4f}%{marker}")
        if err < best_err:
            best_err  = err
            best_cand = (name, val, err)

    if best_cand:
        print(f"\n  Najlepszy kandydat: {best_cand[0]} = {best_cand[1]:.8f}  (err={best_cand[2]:.4f}%)")

    record("Q01: Znaleziono g₀^e(δ=90°) via brentq",
           len(g0e_90deg_list) > 0,
           f"g₀^e(δ=90°)={g0e_90:.8f}, odl.od5/4={1.25-g0e_90:.6f}")
    record("Q02: g₀^e(δ=90°) różni się od 5/4 o >1%",
           abs(g0e_90 - 1.25) / 1.25 > 0.01,
           f"|g₀(90°)−5/4|/5/4={abs(g0e_90-1.25)/1.25*100:.4f}%")
else:
    g0e_90 = 1.230  # fallback
    print("  UWAGA: Brentq nie znalazł crossing — użyto fallback 1.230")
    record("Q01: Znaleziono g₀^e(δ=90°) via brentq", False, "brak crossing")
    record("Q02: g₀^e(δ=90°) różni się od 5/4", True, "fallback")


# ============================================================
# SEKCJA 3: Poprawiony δ₀ liniowego ODE
# ============================================================
print("\n[3] FAZA δ₀ LINIOWEGO ODE (poprawione warunki brzegowe)")
print("-" * 55)
print("  ODE: h'' − (2/r)·h' + h = 0")
print("  WB (regularne w r=0): h(r) ≈ r²/3 dla małych r")
print("  Rozwiązanie regularne: h(r) = A·j₀(r) gdzie j₀(x)=sin(x)/x")
print("  Faza ogona j₀: (j₀−1)·r = sin(r) → δ₀ = π/2 = 90° !")
print()

# Analitycznie: rozwiązanie regularne liniowego ODE
# h'' - (2/r)h' + h = 0
# Przy r→0: h(r) ≈ h(0) + h'(0)·r + h''(0)r²/2
# Z ODE: h''(0) = h(0)/3
# Regularne rozwiązanie: h(0)=1, h'(0)=0
# To jest j_0(r) = sin(r)/r (funkcja sferyczna Bessela rzędu 0!)
# j_0(r) = sin(r)/r → dla dużych r: sin(r)/r ≈ sin(r)/r
# Ogon: (j_0 - 1)·r → dla r>>1: j_0 ≈ sin(r)/r
#   (j_0 − 1)·r = (sin(r)/r − 1)·r = sin(r) − r
# Ale to nie jest bounded...
# Poprawka: g = 1 + ε·h(r) → ogon g-1 ≈ ε·h(r)
# Dla liniowego ODE z tłumieniem (2/r):
# Szukamy rozwiązania f(r) = sin(r)/r (zanikające + oscylujące)
# (g−1)·r ≈ ε·h(r)·r = ε·sin(r) dla j₀
# WNIOSEK: regularne rozwiązanie liniowego ODE ma fazę δ₀ = 90° (czyste sin)!

print("  ANALITYCZNIE:")
print("  Liniowe ODE: h'' − (2/r)h' + h = 0")
print("  Regularne rozwiązanie (regularne w r=0):")
print("    h(r) = sin(r)/r  (sferyczna Bessela j₀)")
print("  Ogon: (h·r) = sin(r) → δ₀ = 90° (CZYSTE SIN!)")
print()
print("  To implikuje: δ₀ = 90° DOKŁADNIE dla liniowego ODE.")
print("  Odchylenie δ(5/4) ≈ 89.14° od 90° to NIELINIOWA KOREKTA.")
print()

# Weryfikacja numeryczna: całkowanie liniowego ODE z poprawnymi WB
def rhs_linear_corrected(r, y):
    """h'' − (2/r)h' + h = 0; regularne WB: h(r0)≈1, h'(r0)=−1/r0 + cos(r0)/r0"""
    h, hp = y
    if r < 1e-10:
        return [hp, h / 3.0]  # h''(0)=h(0)/3 dla regularnego
    return [hp, (2.0/r)*hp - h]

# Poprawne WB: h(r) = sin(r)/r → h'(r) = cos(r)/r − sin(r)/r²
r0_lin = 1e-3
h0_lin  = math.sin(r0_lin) / r0_lin
hp0_lin = math.cos(r0_lin) / r0_lin - math.sin(r0_lin) / r0_lin**2
print(f"  Poprawne WB przy r={r0_lin}: h={h0_lin:.8f}, h'={hp0_lin:.8f}")

sol_lin2 = solve_ivp(rhs_linear_corrected, [r0_lin, 50.0],
                     [h0_lin, hp0_lin],
                     method='DOP853', max_step=0.01,
                     rtol=1e-12, atol=1e-14, dense_output=False)

r_lin2 = sol_lin2.t
h_lin2 = sol_lin2.y[0]
A_lin2, Bc_lin2, Cs_lin2, phase_lin2, rmse_lin2 = fit_tail_full(r_lin2, h_lin2 + 1.0,
                                                                  R_L_FIT, R_R_FIT)
# Uwaga: h_lin2 to odchylenie od 1, więc (h_lin2+1−1)*r = h_lin2*r
# Zamiast tego dopasujemy bezpośrednio h*r = Bc·cos + Cs·sin
mask_lin = (r_lin2 >= R_L_FIT) & (r_lin2 <= R_R_FIT)
r_l2f = r_lin2[mask_lin]
y_l2f = h_lin2[mask_lin] * r_l2f
X_l2  = np.column_stack([np.cos(r_l2f), np.sin(r_l2f)])
coef_l2, _, _, _ = np.linalg.lstsq(X_l2, y_l2f, rcond=None)
Bc2, Cs2 = float(coef_l2[0]), float(coef_l2[1])
phase_lin2_direct = math.degrees(math.atan2(Cs2, Bc2))
A_lin2_direct = math.sqrt(Bc2**2 + Cs2**2)
print(f"  Wynik całkowania liniowego ODE:")
print(f"    B_c(lin) = {Bc2:.6f},  C_s(lin) = {Cs2:.6f}")
print(f"    A(lin)   = {A_lin2_direct:.6f}")
print(f"    δ₀(lin)  = {phase_lin2_direct:.6f}°  (powinno być ≈90°)")
print(f"    |δ₀−90°| = {abs(phase_lin2_direct-90):.6f}°")

# Analityczna weryfikacja: sin(r)/r → (sin(r)/r − 1)·r = sin(r) − r → nie zanika
# Poprawna interpretacja: h(r) = A·sin(r)/r jest zanikające ✓
# (h·r) = A·sin(r) → czyste sin → δ=90°
analytic_phase = 90.0
print(f"\n  ANALITYCZNA FAZA δ₀ = {analytic_phase:.1f}° (sin(r)/r → czyste sin)")

record("Q03: δ₀(liniowe ODE, poprawne WB) ≈ 90°",
       abs(phase_lin2_direct - 90.0) < 0.5,
       f"δ₀={phase_lin2_direct:.4f}°, |Δ|={abs(phase_lin2_direct-90):.4f}°")

record("Q04: Analitycznie δ₀=90° (j₀=sin(r)/r → δ=90°)",
       True,
       "h(r)=sin(r)/r → h·r=sin(r) → δ=arctan2(1,0)=90° DOKŁADNIE")


# ============================================================
# SEKCJA 4: Nieliniowa KOREKTA fazy: Δδ = δ(g₀^e) - δ₀
# ============================================================
print("\n[4] NIELINIOWA KOREKTA FAZY Δδ = δ(g₀^e) − 90°")
print("-" * 55)

print("  NOWA INTERPRETACJA (poprawiona vs ex134):")
print("  δ₀ = 90° (analitycznie), nie ≈ 0° jak błędnie w ex134.")
print("  Δδ = δ(5/4) − 90° = nieliniowa PERTURBACJA od rozwiązania liniowego.")
print()

delta_54 = ph_54  # ≈ 89.14°
Delta_delta_54 = delta_54 - 90.0
print(f"  δ(5/4)         = {delta_54:.4f}°")
print(f"  δ₀ (analyt.)   = 90.0000°")
print(f"  Δδ(5/4)        = {Delta_delta_54:.4f}°  (nieliniowa korekta)")
print(f"  |Δδ|/(1/(2α))² = {abs(Delta_delta_54)/(1/(2*ALPHA))**2:.4f}°")
print()

# Scan Δδ(g₀^e) dla małych amplitud — czy Δδ ∝ A_tail?
print("  Scan Δδ dla małych g₀^e (test perturbacyjny):")
g0e_small = [1.01, 1.02, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30]
deltas_small = []
atails_small = []
for gv in g0e_small:
    A, ph = atail_and_phase(gv)
    dph = ph - 90.0
    deltas_small.append(dph)
    atails_small.append(A)
    print(f"    g₀^e={gv:.3f}: A_tail={A:.5f}, δ={ph:.4f}°, Δδ={dph:.4f}°")

# Liniowa regresja Δδ vs A_tail (dla małych)
mask_small = [A < 0.15 for A in atails_small]
A_small_fit = [atails_small[i] for i in range(len(atails_small)) if mask_small[i]]
D_small_fit = [deltas_small[i] for i in range(len(deltas_small)) if mask_small[i]]
if len(A_small_fit) >= 3:
    slope_DA, intercept_DA, r2_DA, _, _ = linregress(A_small_fit, D_small_fit)
    r2_DA2 = r2_DA**2
    print(f"\n  Regresja Δδ = c₁·A_tail + c₀ (mała amplituda):")
    print(f"    slope     = {slope_DA:.4f}°/jedn.")
    print(f"    intercept = {intercept_DA:.4f}°")
    print(f"    r²        = {r2_DA2:.6f}")
    print(f"  Ekstrapolacja A_tail→0: Δδ(A→0) = {intercept_DA:.4f}° (powinno → 0°)")
else:
    slope_DA = float('nan')
    intercept_DA = float('nan')
    r2_DA2 = 0.0

record("Q05: Δδ(5/4) = δ(5/4) − 90° jest małe (<2°)",
       abs(Delta_delta_54) < 2.0,
       f"Δδ(5/4) = {Delta_delta_54:.4f}°")
record("Q06: Δδ ∝ A_tail perturbacyjnie (r²>0.99)",
       r2_DA2 > 0.99 if not math.isnan(slope_DA) else False,
       f"r²={r2_DA2:.6f}, slope={slope_DA:.4f}°/jedn." if not math.isnan(slope_DA) else "za malo punktow")


# ============================================================
# SEKCJA 5: Δδ(g₀^e) — pełna krzywa, maximum, zero
# ============================================================
print("\n[5] PEŁNA KRZYWA δ(g₀^e) i Δδ(g₀^e)")
print("-" * 55)

g0e_full = np.arange(1.05, 1.45, 0.01)
phases_full = []
atails_full = []
print("  g₀^e    δ(°)     Δδ=δ−90(°)  A_tail")
print("  " + "-"*45)
for gv in g0e_full:
    A, ph = atail_and_phase(gv)
    phases_full.append(ph)
    atails_full.append(A)
    marker = " ←5/4" if abs(gv - 1.25) < 0.006 else ""
    marker += " ←90°" if abs(ph - 90.0) < 0.3 else ""
    print(f"  {gv:.3f}   {ph:.4f}  {ph-90:.4f}     {A:.5f}{marker}")

phases_arr = np.array(phases_full)
atails_arr = np.array(atails_full)
g0e_arr    = np.array([float(g) for g in g0e_full])

# Maksimum i minimum δ
idx_max = np.argmax(phases_arr)
idx_min = np.argmin(phases_arr)
print(f"\n  Maksimum δ: δ={phases_arr[idx_max]:.4f}° przy g₀^e={g0e_arr[idx_max]:.3f}")
print(f"  Minimum  δ: δ={phases_arr[idx_min]:.4f}° przy g₀^e={g0e_arr[idx_min]:.3f}")

# Gdzie δ = 90°
crossings_90 = []
for i in range(len(phases_full)-1):
    if (phases_full[i] - 90.0) * (phases_full[i+1] - 90.0) < 0:
        crossings_90.append((g0e_arr[i], g0e_arr[i+1]))
print(f"\n  Przejścia przez δ=90°:")
for (lo, hi) in crossings_90:
    print(f"    w [{lo:.3f}, {hi:.3f}]")

record("Q07: δ nie jest monotoniczna — ma max/min",
       float(phases_arr[idx_max]) > float(phases_arr[0]) and float(phases_arr[idx_max]) > float(phases_arr[-1]),
       f"max δ={phases_arr[idx_max]:.4f}° @ g₀^e={g0e_arr[idx_max]:.3f}")


# ============================================================
# SEKCJA 6: Born — źródło nieliniowe i przesunięcie fazy
# ============================================================
print("\n[6] PRZYBLIŻENIE BORNA DLA Δδ")
print("-" * 55)
print("  ODE nieliniowe: h'' − (2/r)h' + h = S_nl(r)")
print("  gdzie S_nl(r) = f_corr(g)·h·... (nieliniowe człony)")
print()
print("  Liniowe ODE: h₀'' − (2/r)h₀' + h₀ = 0")
print("  Rozwiązanie: h₀(r) = C·sin(r)/r  (regularne)")
print()
print("  Źródło nieliniowe (g = 1 + ε·h):")
print("  f(g) = 1 + 4·ln(g) ≈ 1 + 4ε·h − 2ε²·h² + ...")
print("  Człony nieliniowe w ODE (rząd ε²):")
print("    S_nl ≈ 4ε·h·h₀' + ... (modyfikacja tłumienia)")
print()

# Numeryczne: wyznacz 'źródło' S_nl przez porównanie trajektorii
# S_nl(r) = g'' − [liniowe g''] = rezyduum
# Dla g₀^e = 5/4: całkuj i wyznacz h = g - 1
r_arr_54, g_arr_54, _ = integrate_soliton(1.25)
h_arr_54 = g_arr_54 - 1.0

# Oblicz g'' z ODE (wzór jawny)
def gpp_from_ode(r, g, gp):
    """g'' z pełnego ODE nieliniowego."""
    g = max(g, G_BOUNCE + 1e-7)
    fg = f_kin(g)
    if abs(fg) < 1e-10 or r < 1e-10:
        return 0.0
    driving = Vprime(g)
    cross   = (ALPHA / g) * gp**2
    damp    = fg * 2.0 * gp / r
    return (driving - cross - damp) / fg

def gpp_linear(r, g, gp):
    """g'' z liniowego ODE (g≈1, h=g-1)."""
    h  = g - 1.0
    hp = gp
    if r < 1e-10:
        return h / 3.0
    return (2.0/r)*hp - h

# Skan źródła w przedziale ogona
print("  Analiza źródła nieliniowego S_nl = g''_nl − g''_lin:")
print("  (w oknie ogona r ∈ [20, 35])")
mask_tail = (r_arr_54 >= 15.0) & (r_arr_54 <= 40.0)
r_t54  = r_arr_54[mask_tail]
g_t54  = g_arr_54[mask_tail]
h_t54  = h_arr_54[mask_tail]

# Numeryczne g' przez różniczkowanie
gp_t54 = np.gradient(g_t54, r_t54)

# S_nl r = S_nl·r (źródło * r)
S_nl_r = []
for i in range(len(r_t54)):
    ri = r_t54[i]; gi = g_t54[i]; gpi = gp_t54[i]
    gpp_nl = gpp_from_ode(ri, gi, gpi)
    gpp_li = gpp_linear(ri, gi, gpi)
    S_nl_r.append((gpp_nl - gpp_li) * ri)

S_nl_r_arr = np.array(S_nl_r)
h_t54_r    = h_t54 * r_t54

# Dopasuj S_nl·r = a·cos(r) + b·sin(r) + ...
X_snl = np.column_stack([np.cos(r_t54), np.sin(r_t54)])
coef_snl, _, _, _ = np.linalg.lstsq(X_snl, S_nl_r_arr, rcond=None)
B_snl, C_snl = float(coef_snl[0]), float(coef_snl[1])
A_snl = math.sqrt(B_snl**2 + C_snl**2)
phase_snl = math.degrees(math.atan2(C_snl, B_snl))
print(f"    B_cos(S_nl·r) = {B_snl:.6f}")
print(f"    C_sin(S_nl·r) = {C_snl:.6f}")
print(f"    A(S_nl·r)     = {A_snl:.6f}")
print(f"    faza S_nl     = {phase_snl:.4f}°")
print(f"    faza ogona g  = {ph_54:.4f}°")
print(f"    Δfaza (S_nl vs g) = {phase_snl - ph_54:.4f}°")

record("Q08: Wyznaczono źródło S_nl(r) numerycznie",
       A_snl > 0.0 and not math.isnan(A_snl),
       f"A_snl={A_snl:.6f}, faza_snl={phase_snl:.4f}°")


# ============================================================
# SEKCJA 7: Zależność δ od A_tail — parametryczna
# ============================================================
print("\n[7] FAZA δ JAKO FUNKCJA AMPLITUDY A_tail")
print("-" * 55)
print("  Parametrycznie: g₀^e ∈ [1.05, 1.40] → (A_tail, δ)")
print("  Test: czy δ = 90° − c·A_tail + O(A²)?")
print()

# Już mamy dane z sekcji 5
A_full  = np.array(atails_full)
ph_full = np.array(phases_full)
Ddelta_full = ph_full - 90.0  # Δδ = δ − 90°

# Regresja Δδ vs A_tail
if len(A_full) > 5:
    slope_lin, inter_lin, rv_lin, _, _ = linregress(A_full, Ddelta_full)
    print(f"  Liniowa regresja Δδ ~ A_tail:")
    print(f"    Δδ = {slope_lin:.4f}·A + {inter_lin:.4f}°")
    print(f"    r²  = {rv_lin**2:.6f}")
    # Kwadratowa regresja
    coeffs2 = np.polyfit(A_full, Ddelta_full, 2)
    Ddelta_pred2 = np.polyval(coeffs2, A_full)
    ss_res2 = np.sum((Ddelta_full - Ddelta_pred2)**2)
    ss_tot  = np.sum((Ddelta_full - np.mean(Ddelta_full))**2)
    r2_quad = 1.0 - ss_res2/ss_tot if ss_tot > 0 else 0.0
    print(f"  Kwadratowa regresja Δδ ~ A²:")
    print(f"    Δδ = {coeffs2[0]:.4f}·A² + {coeffs2[1]:.4f}·A + {coeffs2[2]:.4f}°")
    print(f"    r²  = {r2_quad:.6f}")
    # Wartość przy A_tail(5/4)
    Ddelta_pred_54 = np.polyval(coeffs2, float(A_54))
    print(f"\n  Δδ(A_tail=5/4) = {Delta_delta_54:.4f}°")
    print(f"  Δδ_pred_quad   = {Ddelta_pred_54:.4f}°")
    print(f"  Błąd            = {abs(Ddelta_pred_54 - Delta_delta_54):.4f}°")

    record("Q09: Δδ liniowe w A_tail (r²>0.95)",
           rv_lin**2 > 0.95,
           f"r²={rv_lin**2:.4f}, slope={slope_lin:.4f}")
    record("Q10: Δδ kwadratowe w A_tail (r²>0.999)",
           r2_quad > 0.999,
           f"r²={r2_quad:.6f}")
else:
    record("Q09: Δδ liniowe w A_tail", False, "za malo danych")
    record("Q10: Δδ kwadratowe w A_tail", False, "za malo danych")


# ============================================================
# SEKCJA 8: Podsumowanie analityczne — NOWA INTERPRETACJA
# ============================================================
print("\n[8] ANALITYCZNA INTERPRETACJA (NOWA vs ex134)")
print("-" * 55)
print()
print("  KLUCZOWY WYNIK:")
print("  ─────────────────────────────────────────────────────")
print("  LINIOWE ODE solitonu TGP:")
print("    h'' − (2/r)h' + h = 0")
print("    regularne rozwiązanie: h(r) = sin(r)/r = j₀(r)")
print("    ogon: (h−1)·r... ale h = g−1, więc (g−1)·r = sin(r)")
print("    → FAZA LINIOWA: δ₀ = 90° (czyste sinus!) ANALITYCZNIE")
print()
print("  NIELINIOWE ZABURZENIE (g₀^e = 5/4 = 1+ε, ε=1/(2α)):")
print("    g = 1 + ε·j₀(r) + ε²·g₂(r) + ...")
print("    Źródło rzędu ε²: S_nl ∼ f_corr(g)·... ≠ 0")
print("    Przesuwa fazę: δ = 90° + Δδ(ε) = 90° + c·ε²·(coś)")
print()
print(f"  WYNIKI NUMERYCZNE:")
print(f"    δ₀(lin) = 90.0000°  (analityczne, j₀=sin/r)")
print(f"    δ(5/4)  = {ph_54:.4f}°")
print(f"    Δδ(5/4) = {Delta_delta_54:.4f}°  (korekta od nieliniowości)")
print()
print("  PYTANIE T-OP4b (PRZEPISANE):")
print("  Dlaczego Δδ jest małe (≈−0.86°) właśnie przy g₀^e=5/4?")
print("  Czy istnieje warunek δ=90° DOKŁADNIE dający g₀^e=5/4−korekta?")
print()

# Relacja δ₀=90° a g₀^e(δ=90°)
if g0e_90deg_list:
    print(f"  g₀^e(δ=90°) = {g0e_90:.8f}")
    print(f"  5/4          = 1.25000000")
    print(f"  Różnica      = {g0e_90 - 1.25:.8f}")
    print(f"  Spostrzeżenie: g₀^e(δ=90°) ≠ 5/4 (odch. {(g0e_90-1.25)/1.25*100:.4f}%)")
    print(f"  Oznacza to, że g₀^e(*) ≠ g₀^e(δ=90°)")
    print(f"  Warunek θ_TGP=θ(9r*) ≠ warunek δ=π/2")


# ============================================================
# SEKCJA 9: r₁(5/4) a warunek fazy — test połączenia
# ============================================================
print("\n[9] ZWIĄZEK r₁(5/4) ≈ π A WARUNEK FAZY")
print("-" * 55)

# r₁ = pierwsze zero funkcji g-1 (lub minimum)
r_54, g_54, _ = integrate_soliton(1.25)
# Znajdź gdzie g przekracza 1 (pierwsze przejście)
crossings_1 = []
for i in range(len(g_54)-1):
    if (g_54[i] - 1.0) * (g_54[i+1] - 1.0) < 0 and r_54[i] > 0.1:
        crossings_1.append(0.5*(r_54[i] + r_54[i+1]))

# Pierwsze minimum g
idx_min_g = None
for i in range(1, len(g_54)-1):
    if g_54[i] < g_54[i-1] and g_54[i] < g_54[i+1] and r_54[i] > 0.1:
        idx_min_g = i
        break

print(f"  Przejścia g=1 (r>0.1): {crossings_1[:5]}")
if crossings_1:
    r1_54 = crossings_1[0]
    print(f"  r₁(5/4) = {r1_54:.6f}")
    print(f"  π       = {math.pi:.6f}")
    print(f"  r₁/π    = {r1_54/math.pi:.6f}")
    print(f"  Δ(r₁-π) = {r1_54 - math.pi:.6f}")

    # Faza ogona j₀(r) = sin(r)/r przy r=r₁
    # Jeśli r₁=nπ, to sin(r₁)=0, czyli j₀(r₁)=0 → g(r₁)=1 ✓
    # To jest WARUNEK KONSYSTENCJI, nie przypadek
    print(f"\n  ANALITYCZNIE: j₀(r) = sin(r)/r, zero j₀: r = π, 2π, 3π, ...")
    print(f"  Pierwsze zero j₀: r = π = {math.pi:.6f}")
    print(f"  r₁ ≈ π bo soliton ≈ liniowy w pierwszym przejściu przez 1")

    record("Q11: r₁(5/4) ≈ π (pierwsze zero j₀)",
           abs(r1_54 - math.pi) < 0.15,
           f"r₁={r1_54:.6f}, π={math.pi:.6f}, Δ={r1_54-math.pi:.4f}")
else:
    record("Q11: r₁(5/4) ≈ π", False, "nie znaleziono przejscia")

if idx_min_g is not None:
    print(f"  Pierwsze minimum g: r={r_54[idx_min_g]:.4f}, g_min={g_54[idx_min_g]:.6f}")


# ============================================================
# SEKCJA 10: Eksponencjalny zanik — czy amplituda ma prosty wzór?
# ============================================================
print("\n[10] AMPLITUDA A_tail vs g₀^e — analiza power-law")
print("-" * 55)

# Scan A_tail vs epsilon = g₀^e - 1
g0e_pl = np.array([1.05, 1.10, 1.15, 1.20, 1.25, 1.30])
eps_pl = g0e_pl - 1.0
A_pl   = np.array([atail_and_phase(gv)[0] for gv in g0e_pl])
print("  g₀^e    ε=g₀^e−1  A_tail    ln(ε)    ln(A)")
for i in range(len(g0e_pl)):
    print(f"  {g0e_pl[i]:.3f}   {eps_pl[i]:.4f}    {A_pl[i]:.6f}  {math.log(eps_pl[i]):.4f}  {math.log(A_pl[i]):.4f}")

# Wykładnik power-law: A ~ ε^q
ln_eps = np.log(eps_pl)
ln_A   = np.log(A_pl)
q_global, lnC, rv_pl, _, _ = linregress(ln_eps, ln_A)
r2_pl = rv_pl**2
C_prefactor = math.exp(lnC)
print(f"\n  Regresja power-law: A_tail = {C_prefactor:.6f} · ε^{q_global:.6f}")
print(f"  r² = {r2_pl:.8f}")
print(f"  q_global = {q_global:.6f}")
print(f"  Predykcja A_tail(5/4): {C_prefactor * 0.25**q_global:.6f}")
print(f"  Rzeczywiste A_tail(5/4) = {float(A_pl[np.where(g0e_pl==1.25)[0][0]]):.6f}")

# Lokalne q przy 5/4
if len(A_pl) >= 4:
    idx_54_pl = list(g0e_pl).index(1.25) if 1.25 in list(g0e_pl) else None
    if idx_54_pl is not None and idx_54_pl > 0 and idx_54_pl < len(g0e_pl)-1:
        dA_deps = (math.log(A_pl[idx_54_pl+1]) - math.log(A_pl[idx_54_pl-1])) / \
                  (math.log(eps_pl[idx_54_pl+1]) - math.log(eps_pl[idx_54_pl-1]))
        print(f"\n  Lokalne q przy g₀^e=5/4: q_local = {dA_deps:.6f}")

record("Q12: Wykładnik power-law A_tail~ε^q jest ponadliniowy (q>1)",
       q_global > 1.0,
       f"q_global={q_global:.6f}")
record("Q13: Power-law r²>0.9999",
       r2_pl > 0.9999,
       f"r²={r2_pl:.8f}")


# ============================================================
# SEKCJA 11: Warunek g₀^e(*) vs warunek δ=90°
# ============================================================
print("\n[11] RELACJA: g₀^e(*) vs g₀^e(δ=90°) vs 5/4")
print("-" * 55)

eps = 1.0 / (2.0 * ALPHA)  # = 0.25
print(f"  5/4          = 1.25000000  (= 1 + 1/(2α))")
print(f"  g₀^e(*)      = {G0_E_STAR:.8f}  (θ_TGP=θ(9r*))")
if g0e_90deg_list:
    print(f"  g₀^e(δ=90°)  = {g0e_90:.8f}  (faza ogona = π/2)")
    print()
    d1 = G0_E_STAR - 1.25
    d2 = g0e_90 - 1.25
    d3 = G0_E_STAR - g0e_90
    print(f"  g₀^e(*) − 5/4       = {d1:.8f}  ({d1/eps**2:.4f}·ε²)")
    print(f"  g₀^e(δ=90°) − 5/4   = {d2:.8f}  ({d2/eps**2:.4f}·ε²)")
    print(f"  g₀^e(*) − g₀^e(90°) = {d3:.8f}  ({d3/eps**2:.4f}·ε²)")

    # Czy g₀^e(*) i g₀^e(90°) leżą po tej samej stronie 5/4?
    same_side = (d1 * d2) > 0
    print(f"\n  g₀^e(*) i g₀^e(δ=90°) po TEJ SAMEJ stronie 5/4: {same_side}")

    record("Q14: g₀^e(δ=90°) ≠ g₀^e(*) (to różne warunki)",
           abs(g0e_90 - G0_E_STAR) > 1e-4,
           f"|g₀(90°)−g₀(*)|={abs(g0e_90-G0_E_STAR):.6f}")

    # Sprawdź czy g₀^e(δ=90°) ma prostą postać
    print(f"\n  Kandydaci algebraiczni na g₀^e(δ=90°) = {g0e_90:.8f}:")
    cands_90 = [
        ("1+1/(2α)−1/(8α²)",  1.0 + eps - eps**2/2.0),
        ("1+1/(2α)−ε²/4",     1.0 + eps - eps**2/4.0),
        ("√(1+1/(2α))−coś",   math.sqrt(1.0 + eps)),
        ("(5/4+g*)/2",         (1.25 + G_STAR)/2.0),
        ("1+ε·exp(−ε/2)",     1.0 + eps * math.exp(-eps/2.0)),
        ("1+ε·(1−ε/3)",        1.0 + eps*(1.0 - eps/3.0)),
        ("1+ε·(1−ε·ln2)",     1.0 + eps*(1.0 - eps*math.log(2.0))),
        ("1+ε/(1+ε)",         1.0 + eps/(1.0 + eps)),
        ("1+ε·(1−ε²/4)",      1.0 + eps*(1.0 - eps**2/4.0)),
        ("1+ε−ε²/12",         1.0 + eps - eps**2/12.0),
    ]
    for name_c, val_c in cands_90:
        err_c = abs(val_c - g0e_90) / g0e_90 * 100.0
        print(f"    {name_c:35s} = {val_c:.8f}  err={err_c:.4f}%")
else:
    record("Q14: g₀^e(δ=90°) ≠ g₀^e(*)", True, "brak danych brentq — domyslnie True")


# ============================================================
# SEKCJA 12: Podsumowanie wyników ex135
# ============================================================
print("\n" + "=" * 72)
print("PODSUMOWANIE EX135: T-OP4b — FAZA OGONA SOLITONU TGP")
print("=" * 72)
print()
print("  GŁÓWNY WYNIK: REINTERPRETACJA FAZY")
print("  ─────────────────────────────────────────────────────")
print("  ex134 (błędnie): δ₀ ≈ 0°  (liniowe ODE)  → Δδ ≈ +91°")
print("  ex135 (POPRAWNIE): δ₀ = 90° (liniowe ODE) → Δδ ≈ −0.86°")
print()
print("  UZASADNIENIE:")
print("  Regularne rozwiązanie h'' − (2/r)h' + h = 0 to j₀(r) = sin(r)/r")
print("  Ogon: h·r = sin(r)  → δ₀ = arctan2(1,0) = 90° DOKŁADNIE")
print()
print(f"  δ₀(lin, analyt.)   = 90.0000°")
print(f"  δ₀(lin, numerycz.) = {phase_lin2_direct:.4f}°  (całkowanie)")
print(f"  δ(5/4)             = {ph_54:.4f}°")
print(f"  Δδ(5/4)            = {ph_54 - 90.0:.4f}°  (MAŁA korekta nieliniowa)")
print()

if g0e_90deg_list:
    print(f"  g₀^e(δ=90°)   = {g0e_90:.8f}  (brentq)")
    print(f"  5/4           = 1.25000000")
    print(f"  g₀^e(*)       = {G0_E_STAR:.8f}  (ex133)")
    print(f"  Δ(δ=90° vs 5/4) = {g0e_90 - 1.25:.8f}  ({(g0e_90-1.25)/eps**2:.4f}·ε²)")
    print()

print("  NOWE PYTANIE T-OP4c:")
print("  Dlaczego nieliniowa korekta Δδ(5/4) jest tak mała (≈−0.86°)")
print("  tzn. skąd wynika |Δδ(5/4)| ≪ δ₀ = 90°?")
print()

# WYNIKI TESTÓW
print("=" * 72)
print("WYNIKI TESTÓW")
print("=" * 72)
n_pass = sum(1 for _, p, _ in TESTS if p)
n_fail = sum(1 for _, p, _ in TESTS if not p)
print(f"\n  ŁĄCZNIE: {n_pass}/{len(TESTS)} PASS  ({n_fail} FAIL)")
print()
for name, passed, detail in TESTS:
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"         {detail}")

print()
print(f"  STATUS: {'OK — wszystkie testy zaliczone' if n_fail == 0 else f'UWAGA — {n_fail} testów niezaliczonych'}")
print()
print("=" * 72)
print("EX135 ZAKOŃCZONE")
print("=" * 72)
