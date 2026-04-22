#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex136_top4_alpha_scan.py
========================
BADANIE T-OP4: ZALEŻNOŚĆ g₀^e(*) OD α — ROZWINIĘCIE W SZEREG

KONTEKST:
  ex133 (α=2): g₀^{e,*} = 1.249082,  5/4 = 1.250000,  Δ = -0.000918
  ex135: Δδ(5/4) = -0.864° (mała korekta nieliniowa)
  Korekta 2. rzędu: g₀^{e,*} = 1 + ε + c₂·ε² + O(ε³) gdzie ε=1/(2α)
  Dla α=2: c₂ = (g₀^{e,*} - 1 - ε) / ε² = -0.000918/0.0625 = -0.01469

PLAN (16 testów):
  A01-A03: Wysokoprecyzyjne wyznaczenie g₀^e(*) przy α=2
  A04-A06: Skan α = 1.0,1.5,2.0,2.5,3.0,4.0,5.0,10.0
  A07-A10: Wyekstrahowanie c₂(α) i sprawdzenie czy c₂ = const
  A11-A13: Kandydaty algebraiczne na c₂
  A14-A16: Weryfikacja: ekspansja g₀^{e,*}(α) = 1+ε−c₂·ε² vs numeryka
"""

import sys
import io
import math
import cmath
import warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.stats import linregress

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

warnings.filterwarnings('ignore')

# ============================================================
# Stałe globalne (α=2, wartości fizyczne)
# ============================================================
ALPHA_PHYS = 2.0
PHI        = (1.0 + math.sqrt(5.0)) / 2.0
RSTAR      = (23.0 + 5.0 * math.sqrt(21.0)) / 2.0  # ≈ 22.9564
THETA_9RS  = 132.731439   # θ(9r*) z ex129

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
# ODE solitonu — parametryczne w α
# ============================================================
def make_rhs(alpha):
    """Zwraca funkcję prawej strony ODE dla danego α."""
    g_star   = math.exp(-1.0 / (2.0 * alpha))
    g_bounce = g_star + 0.005

    def f_kin(g):
        return 1.0 + 2.0 * alpha * math.log(max(g, 1e-30))

    def Vprime(g):
        return g * g * (1.0 - g)

    def rhs(r, y):
        g, gp = y
        g = max(g, g_bounce + 1e-7)
        fg = f_kin(g)
        if abs(fg) < 1e-10:
            return [gp, 0.0]
        driving = Vprime(g)
        cross   = (alpha / g) * gp**2
        if r < 1e-10:
            return [gp, (driving - cross) / (3.0 * fg)]
        damp = fg * 2.0 * gp / r
        return [gp, (driving - cross - damp) / fg]

    def event_ghost(r, y):
        return y[0] - g_bounce
    event_ghost.terminal  = True
    event_ghost.direction = -1

    return rhs, event_ghost, g_star, g_bounce


def integrate_soliton_alpha(g0, alpha, r_max=None, max_bounces=25):
    rhs_fn, event_fn, g_star, g_bounce = make_rhs(alpha)
    if r_max is None:
        r_max = max(R_MAX, 15.0 * g0)
    r0, y0 = R_START, [g0, 0.0]
    segs_r, segs_g = [], []
    for bn in range(max_bounces + 1):
        sol = solve_ivp(rhs_fn, [r0, r_max], y0,
                        method='DOP853', max_step=MAX_STEP,
                        rtol=RTOL, atol=ATOL,
                        events=[event_fn], dense_output=False)
        segs_r.append(sol.t)
        segs_g.append(sol.y[0])
        if sol.t_events[0].size > 0 and bn < max_bounces:
            r_b  = float(sol.t_events[0][0])
            gp_b = float(sol.y_events[0][0, 1])
            r0   = r_b + 1e-6
            y0   = [g_bounce + 1e-5, -gp_b]
        else:
            break
    r_all = np.concatenate(segs_r)
    g_all = np.concatenate(segs_g)
    idx   = np.argsort(r_all)
    return r_all[idx], g_all[idx]


def fit_atail(r_arr, g_arr, r_L=20.0, r_R=35.0):
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    if np.sum(mask) < 10:
        return float('nan')
    r_f = r_arr[mask]
    y_f = (g_arr[mask] - 1.0) * r_f
    X   = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, _, _, _ = np.linalg.lstsq(X, y_f, rcond=None)
    return math.sqrt(float(coef[0])**2 + float(coef[1])**2)


def atail_alpha(g0, alpha):
    r_arr, g_arr = integrate_soliton_alpha(g0, alpha)
    return fit_atail(r_arr, g_arr)


def r21_for_g0e(g0e, alpha):
    """r₂₁ = (A_tail(φ·g₀^e)/A_tail(g₀^e))^4."""
    g0mu = PHI * g0e
    Ae  = atail_alpha(g0e,  alpha)
    Amu = atail_alpha(g0mu, alpha)
    return (Amu / Ae)**4


def theta_tgp_for_g0e(g0e, alpha):
    """θ_TGP z r₂₁ i r₃₁ = (r₂₁)^(φ/1) na φ-drabinie."""
    g0tau = PHI**2 * g0e
    Ae   = atail_alpha(g0e,   alpha)
    Amu  = atail_alpha(PHI*g0e,  alpha)
    Atau = atail_alpha(g0tau, alpha)
    if Ae < 1e-15 or Amu < 1e-15 or Atau < 1e-15:
        return float('nan')
    me  = Ae**4
    mmu = Amu**4
    mta = Atau**4
    # Oblicz θ_Brannen wzorem Fourierowskim (identyczny z ex133)
    sqm = np.array([math.sqrt(me), math.sqrt(mmu), math.sqrt(mta)])
    M_br  = float(np.mean(sqm))
    eps_br = sqm / M_br - 1.0
    F1_br  = sum(eps_br[k] * cmath.exp(-2j * math.pi * k / 3) for k in range(3))
    return math.degrees(math.atan2(F1_br.imag, F1_br.real))


def brannen_theta_direct(m1, m2, m3):
    """Kąt Brannena z mas (Fourier, identyczny z ex133)."""
    sqm = np.array([math.sqrt(m1), math.sqrt(m2), math.sqrt(m3)])
    M   = float(np.mean(sqm))
    eps = sqm / M - 1.0
    F1  = sum(eps[k] * cmath.exp(-2j * math.pi * k / 3) for k in range(3))
    return math.degrees(math.atan2(F1.imag, F1.real))


# Fizyczne proporcje z ex133 (g₀^μ/g₀^e, g₀^τ/g₀^e)
R_MU_PHYS  = 2.02117 / 1.24915   # ≈ 1.6181 ≈ φ
R_TAU_PHYS = 3.18912 / 1.24915   # ≈ 2.5530 ≈ ξ*


def theta_tgp_v2(g0e, alpha, r_mu=None, r_tau=None):
    """θ_TGP z amplitud solitonów e/μ/τ — wzór Fourierowski (ex133)."""
    if r_mu  is None: r_mu  = R_MU_PHYS
    if r_tau is None: r_tau = R_TAU_PHYS
    g0mu  = r_mu  * g0e
    g0tau = r_tau * g0e
    Ae   = atail_alpha(g0e,   alpha)
    Amu  = atail_alpha(g0mu,  alpha)
    Atau = atail_alpha(g0tau, alpha)
    if any((not math.isfinite(A)) or A <= 0 for A in [Ae, Amu, Atau]):
        return float('nan')
    me  = Ae**4
    mmu = Amu**4
    mta = Atau**4
    return brannen_theta_direct(me, mmu, mta)


# ============================================================
print("=" * 72)
print("EX136: ZALEŻNOŚĆ g₀^e(*) OD α — ROZWINIĘCIE g₀^{e,*}=1+ε−c₂·ε²")
print("=" * 72)
print()
print(f"  α_phys    = {ALPHA_PHYS}")
print(f"  ε_phys    = 1/(2α) = {1/(2*ALPHA_PHYS):.6f}")
print(f"  5/4       = 1.250000")
print(f"  g₀^e(*)   = 1.249082  (ex133)")
print(f"  Δ = g₀^e(*)−5/4 = -0.000918")
print(f"  c₂ = Δ/ε² = {(-0.000918)/0.0625:.6f}")
print(f"  RSTAR     = {RSTAR:.8f}")
print(f"  θ(9r*)    = {THETA_9RS:.6f}°")
print()


# ============================================================
# SEKCJA 1: Weryfikacja θ_TGP przy α=2, g₀^e=5/4 i g₀^e(*)
# ============================================================
print("[1] WERYFIKACJA θ_TGP PRZY α=2")
print("-" * 55)

# Weryfikacja ex133: θ_TGP(g₀^e=1.249082) = THETA_9RS
ALPHA = ALPHA_PHYS
G0E_STAR_EX133 = 1.249082

print("  Obliczam θ_TGP(g₀^e=1.249082) ...")
theta_star = theta_tgp_v2(G0E_STAR_EX133, ALPHA)
print(f"  θ_TGP({G0E_STAR_EX133}) = {theta_star:.6f}°  (oczekiwane: {THETA_9RS:.6f}°)")
print(f"  Δθ = {theta_star - THETA_9RS:.6f}°")

print("\n  Obliczam θ_TGP(g₀^e=5/4=1.25) ...")
theta_54 = theta_tgp_v2(1.25, ALPHA)
print(f"  θ_TGP(5/4) = {theta_54:.6f}°  (oczekiwane: {THETA_9RS:.6f}°)")
print(f"  Δθ(5/4) = {theta_54 - THETA_9RS:.6f}°")

def theta_minus_9rs(g0e, alpha):
    return theta_tgp_v2(g0e, alpha) - THETA_9RS

record("A01: θ_TGP(g₀^e=1.249082) ≈ θ(9r*) (|Δ|<0.05°)",
       abs(theta_star - THETA_9RS) < 0.05,
       f"θ={theta_star:.6f}°, θ(9r*)={THETA_9RS:.6f}°, Δ={theta_star-THETA_9RS:.6f}°")
record("A02: θ_TGP(5/4) ≠ θ(9r*) (są różne)",
       abs(theta_54 - THETA_9RS) > 0.001,
       f"θ(5/4)={theta_54:.6f}°, Δ={theta_54-THETA_9RS:.6f}°")


# ============================================================
# SEKCJA 2: Wysokoprecyzyjne brentq dla α=2
# ============================================================
print("\n[2] WYSOKOPRECYZYJNE g₀^e(*) PRZY α=2 (brentq)")
print("-" * 55)

# Sprawdź przedział dla brentq
lo, hi = 1.20, 1.30
f_lo = theta_minus_9rs(lo, ALPHA)
f_hi = theta_minus_9rs(hi, ALPHA)
print(f"  f({lo}) = {f_lo:.6f}°  (θ−θ(9r*))")
print(f"  f({hi}) = {f_hi:.6f}°  (θ−θ(9r*))")

if f_lo * f_hi < 0:
    g0e_star = brentq(lambda g: theta_minus_9rs(g, ALPHA), lo, hi,
                      xtol=1e-10, rtol=1e-12)
    theta_at_star = theta_tgp_v2(g0e_star, ALPHA)
    eps_phys  = 1.0 / (2.0 * ALPHA)
    delta_g   = g0e_star - 1.0 - eps_phys
    c2_phys   = delta_g / eps_phys**2
    print(f"\n  g₀^e(*) = {g0e_star:.10f}  (brentq, xtol=1e-10)")
    print(f"  θ_TGP   = {theta_at_star:.8f}°")
    print(f"  5/4     = 1.25000000")
    print(f"  Δ       = g₀^e(*) − 5/4 = {g0e_star - 1.25:.10f}")
    print(f"  ε₀=1/(2α) = {eps_phys:.8f}")
    print(f"  c₂ = (g₀^e(*) − 1 − ε₀) / ε₀² = {c2_phys:.8f}")
    print(f"  c₂·ε₀²   = {c2_phys * eps_phys**2:.10f}  (korekta)")
    print()
    # Szukamy c₂ jako ułamka
    from fractions import Fraction
    c2_frac = Fraction(c2_phys).limit_denominator(200)
    print(f"  Najlepszy ułamek dla c₂ ≈ {c2_phys:.6f}:")
    print(f"    {c2_frac} = {float(c2_frac):.8f}  (błąd: {abs(float(c2_frac)-c2_phys):.6f})")
    # Inne kandydaty
    print(f"  Inne kandydaty algebraiczne na c₂:")
    cands_c2 = [
        ("−1/68",     -1.0/68),
        ("−1/64",     -1.0/64),
        ("−3/(4·π²)", -3.0/(4*math.pi**2)),
        ("−(φ−1)/φ²", -(PHI-1)/PHI**2),
        ("−1/(4·PHI²)", -1.0/(4*PHI**2)),
        ("−1/(4√21)",  -1.0/(4*math.sqrt(21))),
        ("−1/(RSTAR)", -1.0/RSTAR),
        ("−1/70",     -1.0/70),
        ("−ln(5/4)/5",-math.log(1.25)/5),
        ("−ε₀/17",    -eps_phys/17),
        ("−1/(8+2√21)", -1.0/(8+2*math.sqrt(21))),
    ]
    for name_c, val_c in cands_c2:
        err_abs = abs(val_c - c2_phys)
        err_rel = err_abs / abs(c2_phys) * 100.0
        print(f"    {name_c:25s} = {val_c:.8f}  err_abs={err_abs:.6f} ({err_rel:.2f}%)")

    record("A03: Wysokoprecyzyjne g₀^e(*) wyznaczone",
           abs(theta_at_star - THETA_9RS) < 1e-5,
           f"g₀^e(*)={g0e_star:.10f}, Δ={g0e_star-1.25:.10f}, c₂={c2_phys:.6f}")
else:
    g0e_star = G0E_STAR_EX133
    c2_phys  = -0.01469
    print(f"  UWAGA: nie znaleziono crossing — użyto ex133 wartości")
    record("A03: Wysokoprecyzyjne g₀^e(*)", False, "brak crossing w [1.20,1.30]")


# ============================================================
# SEKCJA 3: SKAN α — g₀^e(*)(α)
# ============================================================
print("\n[3] SKAN α: g₀^e(α) dla różnych α")
print("-" * 55)
print("  α      ε=1/(2α)  g₀^e(*)   1+ε     Δg     c₂=Δg/ε²")
print("  " + "-"*62)

alpha_vals = [1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0]
results_alpha = []  # (alpha, eps, g0e_star, c2)

for alpha_i in alpha_vals:
    eps_i = 1.0 / (2.0 * alpha_i)
    g_star_i = math.exp(-1.0 / (2.0 * alpha_i))

    # Wyznaczenie g₀^e(*) przez brentq
    # Przedział zależy od α — dla małego α ghost wall jest bliżej 1,
    # więc solitony są bliższe 1+
    lo_i = 1.0 + 0.5 * eps_i
    hi_i = 1.0 + 3.0 * eps_i

    try:
        f_lo_i = theta_minus_9rs(lo_i, alpha_i)
        f_hi_i = theta_minus_9rs(hi_i, alpha_i)

        if math.isnan(f_lo_i) or math.isnan(f_hi_i):
            print(f"  α={alpha_i:.1f}: NaN w θ_TGP — pomijam")
            continue

        if f_lo_i * f_hi_i >= 0:
            # Rozszerzamy przedział
            for expand in [5, 10, 20]:
                hi_try = 1.0 + expand * eps_i
                lo_try = 1.0 + 0.1 * eps_i
                f_lo_try = theta_minus_9rs(lo_try, alpha_i)
                f_hi_try = theta_minus_9rs(hi_try, alpha_i)
                if f_lo_try * f_hi_try < 0 and not math.isnan(f_lo_try) and not math.isnan(f_hi_try):
                    lo_i, hi_i = lo_try, hi_try
                    break
            else:
                print(f"  α={alpha_i:.1f}: brak crossing — pomijam")
                continue

        g0e_i = brentq(lambda g: theta_minus_9rs(g, alpha_i), lo_i, hi_i,
                       xtol=1e-9, rtol=1e-10)
        delta_i = g0e_i - 1.0 - eps_i
        c2_i    = delta_i / eps_i**2
        results_alpha.append((alpha_i, eps_i, g0e_i, c2_i))
        print(f"  α={alpha_i:4.1f}  ε={eps_i:.4f}  g₀^e(*)={g0e_i:.7f}  1+ε={1+eps_i:.7f}  Δg={delta_i:.7f}  c₂={c2_i:.5f}")
    except Exception as e:
        print(f"  α={alpha_i:.1f}: BŁĄD: {e}")

print()

# Analiza c₂(α)
if len(results_alpha) >= 3:
    c2_vals  = np.array([r[3] for r in results_alpha])
    eps_vals = np.array([r[1] for r in results_alpha])
    alpha_v  = np.array([r[0] for r in results_alpha])
    c2_mean  = float(np.mean(c2_vals))
    c2_std   = float(np.std(c2_vals))
    print(f"  Statystyki c₂:")
    print(f"    średnia c₂ = {c2_mean:.6f}")
    print(f"    odch.std   = {c2_std:.6f}")
    print(f"    cv         = {c2_std/abs(c2_mean)*100:.2f}%")

    # Regresja c₂ = a + b·ε (czy c₂ zależy od α?)
    if len(results_alpha) >= 4:
        slope_c2, inter_c2, rv_c2, _, _ = linregress(eps_vals, c2_vals)
        print(f"\n  Regresja c₂ vs ε: c₂ = {inter_c2:.6f} + {slope_c2:.4f}·ε, r²={rv_c2**2:.4f}")
        print(f"  (r²≈1 → c₂ zależy od α; r²≈0 → c₂ stałe)")

    record("A04: c₂ jest przybliżenie stałe niezaleznie od α (cv<20%)",
           c2_std / abs(c2_mean) < 0.20,
           f"c₂_mean={c2_mean:.5f}, std={c2_std:.5f}, cv={c2_std/abs(c2_mean)*100:.1f}%")
    record("A05: c₂ < 0 (korekta 2. rzędu jest ujemna)",
           c2_mean < 0,
           f"c₂={c2_mean:.5f}")
    record("A06: Znaleziono g₀^e(*) dla ≥3 wartości α",
           len(results_alpha) >= 3,
           f"N={len(results_alpha)} punktów")
else:
    c2_mean = c2_phys
    c2_std  = 0.0
    record("A04: c₂ stałe", False, "za mało punktów")
    record("A05: c₂ < 0", c2_phys < 0, f"c₂_phys={c2_phys:.5f}")
    record("A06: Znaleziono g₀^e(*) dla ≥3 α", False, "za mało danych")


# ============================================================
# SEKCJA 4: Algebraiczne kandydaty na c₂ (używając c₂_mean)
# ============================================================
print("\n[4] ALGEBRAICZNE KANDYDATY NA c₂ (uśrednione po α)")
print("-" * 55)
c2_target = c2_mean
print(f"  c₂ (target) = {c2_target:.8f}")
print()

phi = PHI
sqrt21 = math.sqrt(21)
cands_c2_full = [
    # Proste ułamki
    ("−1/68",         -1.0/68),
    ("−1/64",         -1.0/64),
    ("−1/70",         -1.0/70),
    ("−1/72",         -1.0/72),
    ("−1/60",         -1.0/60),
    # Związane z φ
    ("−1/φ^4",        -1.0/phi**4),
    ("−1/(4φ²)",      -1.0/(4*phi**2)),
    ("−(φ−1)/φ³",    -(phi-1)/phi**3),
    ("−φ/16",         -phi/16),
    # Związane z π
    ("−π/216",        -math.pi/216),
    ("−3/(4π²)",      -3.0/(4*math.pi**2)),
    ("−π/(20√21)",    -math.pi/(20*sqrt21)),
    # Związane z √21
    ("−1/(4√21)",     -1.0/(4*sqrt21)),
    ("−√21/432",      -sqrt21/432),
    ("−1/(8+2√21)",   -1.0/(8+2*sqrt21)),
    # Związane z r*
    ("−1/(2r*)",      -1.0/(2*RSTAR)),
    ("−1/RSTAR",      -1.0/RSTAR),
    # Specjalne kombinacje
    ("−(√21-4)/8",    -(sqrt21-4)/8),
    ("−1/(3+√21)",    -1.0/(3+sqrt21)),
    ("−(5−√21)/32",   -(5-sqrt21)/32),
    ("−1/64·(1+1/4α)",-1.0/64*(1+1/(4*ALPHA_PHYS))),
    ("−ε/(2·RSTAR)",  -(1/(2*ALPHA_PHYS))/(2*RSTAR)),
    ("−ln(φ)/10",     -math.log(phi)/10),
    ("−ln(5/4)/5",    -math.log(5.0/4.0)/5),
]

best = sorted(cands_c2_full, key=lambda x: abs(x[1]-c2_target))[:5]
print(f"  TOP 5 kandydatów (posortowane wg błędu):")
for name_c, val_c in cands_c2_full:
    err_abs = abs(val_c - c2_target)
    err_rel = err_abs / abs(c2_target) * 100.0
    marker = " ← TOP" if (name_c, val_c) in best else ""
    print(f"    {name_c:25s} = {val_c:.8f}  err={err_rel:.2f}%{marker}")

best_name, best_val = best[0]
best_err = abs(best_val - c2_target) / abs(c2_target) * 100
print(f"\n  Najlepszy kandydat: {best_name} = {best_val:.8f}  (err={best_err:.3f}%)")

record("A07: Najlepszy kandydat na c₂ z err<3%",
       best_err < 3.0,
       f"{best_name}={best_val:.6f}, c₂={c2_target:.6f}, err={best_err:.3f}%")
record("A08: Kandydat −1/68 w błędzie <2%",
       abs(-1.0/68 - c2_target)/abs(c2_target)*100 < 2.0,
       f"−1/68={-1/68:.6f}, err={abs(-1/68-c2_target)/abs(c2_target)*100:.2f}%")


# ============================================================
# SEKCJA 5: Ekspansja A_tail(g₀^e) w ε = g₀^e − 1
# ============================================================
print("\n[5] EKSPANSJA A_tail(g₀^e) = a₁·ε + a₂·ε² + ... PRZY α=2")
print("-" * 55)

ALPHA = ALPHA_PHYS
eps_vals_exp = [0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35]
A_vals_exp  = []
print("  ε = g₀^e−1  A_tail    ln(ε)    ln(A)")
for eps in eps_vals_exp:
    g0e_i = 1.0 + eps
    A_i = atail_alpha(g0e_i, ALPHA)
    A_vals_exp.append(A_i)
    print(f"  ε={eps:.3f}  A={A_i:.6f}  lnε={math.log(eps):.4f}  lnA={math.log(max(A_i,1e-30)):.4f}")

# Power law fit: A = C·ε^q
ln_eps_exp = np.log(np.array(eps_vals_exp))
ln_A_exp   = np.log(np.array(A_vals_exp))
q_fit, lnC_fit, rv_fit, _, _ = linregress(ln_eps_exp, ln_A_exp)
C_fit = math.exp(lnC_fit)
r2_fit = rv_fit**2
print(f"\n  Power-law: A = {C_fit:.6f}·ε^{q_fit:.6f}  (r²={r2_fit:.8f})")

# Lokalne q przy ε=0.25 (g₀^e=5/4)
idx_25 = eps_vals_exp.index(0.25) if 0.25 in eps_vals_exp else None
if idx_25 is not None and idx_25 > 0 and idx_25 < len(eps_vals_exp)-1:
    q_local_54 = (math.log(A_vals_exp[idx_25+1]) - math.log(A_vals_exp[idx_25-1])) / \
                 (math.log(eps_vals_exp[idx_25+1]) - math.log(eps_vals_exp[idx_25-1]))
    print(f"  Lokalne q przy ε=0.25: q_local={q_local_54:.6f}")

# Korekta 2. rzędu: A = a₁ε + a₂ε² + a₃ε³ (polinomialnie)
eps_arr_p = np.array(eps_vals_exp)
A_arr_p   = np.array(A_vals_exp)
coeffs_A = np.polyfit(eps_arr_p, A_arr_p, 3)
a3, a2, a1, a0 = coeffs_A
print(f"\n  Wielomian 3. stopnia: A = {a0:.6f} + {a1:.6f}·ε + {a2:.6f}·ε² + {a3:.6f}·ε³")
print(f"    a₀ (stała) = {a0:.6f}  (powinno → 0)")
print(f"    a₁ (liniowe) = {a1:.6f}")
print(f"    a₂ (kwad.)  = {a2:.6f}")
print(f"    a₃ (sześć.) = {a3:.6f}")
# Stosunek a₂/a₁
print(f"    a₂/a₁ = {a2/a1:.6f}")

# Sprawdź predykcję A_tail(5/4) = A(ε=0.25)
A_pred_54 = a0 + a1*0.25 + a2*0.25**2 + a3*0.25**3
A_actual_54 = atail_alpha(1.25, ALPHA)
print(f"\n  A_tail(5/4): predykcja={A_pred_54:.6f}, rzeczywiste={A_actual_54:.6f}")
print(f"  Błąd: {abs(A_pred_54 - A_actual_54):.6f} ({abs(A_pred_54-A_actual_54)/A_actual_54*100:.3f}%)")

record("A09: A_tail power-law r²>0.9999",
       r2_fit > 0.9999,
       f"r²={r2_fit:.8f}, q={q_fit:.6f}")
record("A10: a₀≈0 (brak stałego offsetu)",
       abs(a0) < 0.005,
       f"a₀={a0:.6f}")


# ============================================================
# SEKCJA 6: Warunek na g₀^e(*) w funkcji A_tail
# ============================================================
print("\n[6] WARUNEK r₂₁ = 9r* A IMPLIKACJA DLA g₀^e(*)")
print("-" * 55)
print()
print("  Warunek: A_tail(φ·g₀^e)^4 / A_tail(g₀^e)^4 = 9r*")
print("  Używając A_tail(g₀^e) ≈ C·ε^q:")
print()
print("  A_tail(φ·g₀^e) / A_tail(g₀^e) ≈ (φ·g₀^e - 1)^q / (g₀^e - 1)^q")
print("                                  = ((φ-1)/ε + φ)^q · (coś)")
print()
print("  Precyzyjnie: dla g₀^e = 1+ε:")
print("    g₀^μ = φ·(1+ε) = φ + φ·ε")
print("    g₀^μ - 1 = (φ-1) + φ·ε = ε_μ")
print()

eps_target = 0.25  # = 1/(2α) dla α=2
eps_mu     = (PHI - 1.0) + PHI * eps_target
r21_approx = (eps_mu / eps_target)**4  # liniowe przybliżenie
print(f"  ε_e = 1/(2α) = {eps_target:.6f}")
print(f"  ε_μ = (φ-1) + φ·ε_e = {eps_mu:.6f}")
print(f"  r₂₁(liniowe) = (ε_μ/ε_e)^4 = {r21_approx:.4f}")
print(f"  9r* = {9*RSTAR:.4f}")
print(f"  Odchylenie: {r21_approx - 9*RSTAR:.4f} ({(r21_approx-9*RSTAR)/(9*RSTAR)*100:.3f}%)")
print()
print("  Power-law: A ≈ C·ε^q (q>1) zmienia proporcje:")
q_pw = q_fit
r21_power = (eps_mu / eps_target)**(4*q_pw)
print(f"  q = {q_pw:.6f}")
print(f"  r₂₁(power-law) = (ε_μ/ε_e)^(4q) = {r21_power:.4f}")
print(f"  9r* = {9*RSTAR:.4f}")
print(f"  Odchylenie: {r21_power - 9*RSTAR:.4f} ({(r21_power-9*RSTAR)/(9*RSTAR)*100:.3f}%)")
print()
print("  Wniosek: warunek r₂₁=9r* nie daje ε=1/(2α) w przybliżeniu power-law.")
print("           Konieczne jest pełne numeryczne rozwiązanie ODE.")

record("A11: r₂₁(liniowe, ε=1/(2α)) ≠ 9r* (różnią się >5%)",
       abs(r21_approx - 9*RSTAR) / (9*RSTAR) * 100 > 5,
       f"r₂₁={r21_approx:.4f}, 9r*={9*RSTAR:.4f}, Δ={(r21_approx-9*RSTAR)/(9*RSTAR)*100:.2f}%")
record("A12: r₂₁(power-law) bliższe 9r* niż liniowe",
       abs(r21_power - 9*RSTAR) < abs(r21_approx - 9*RSTAR),
       f"r₂₁(pl)={r21_power:.4f}, r₂₁(lin)={r21_approx:.4f}, 9r*={9*RSTAR:.4f}")


# ============================================================
# SEKCJA 7: Symetria — f(g₀^e(*)) przy różnych α
# ============================================================
print("\n[7] WARUNEK f(g₀^e(*)) — czy ma prostą wartość?")
print("-" * 55)

def f_kin_alpha(g, alpha):
    return 1.0 + 2.0 * alpha * math.log(max(g, 1e-30))

print("  α     g₀^e(*)    f(g₀^e(*))  1/(2α)  f·ε  f(5/4)  f-2")
for (alpha_i, eps_i, g0e_i, c2_i) in results_alpha:
    f_at_star = f_kin_alpha(g0e_i, alpha_i)
    f_54      = f_kin_alpha(1.0 + eps_i, alpha_i)  # f at g=1+1/(2α)
    print(f"  α={alpha_i:.1f}  g₀^e(*)={g0e_i:.6f}  f={f_at_star:.6f}  "
          f"ε={eps_i:.4f}  f·ε={f_at_star*eps_i:.6f}  f(1+ε)={f_54:.6f}  f-2={f_at_star-2:.6f}")

# Sprawdź f(g₀^e(*)) dla α=2
if len(results_alpha) >= 2:
    alpha_phys_res = [r for r in results_alpha if abs(r[0]-ALPHA_PHYS)<0.01]
    if alpha_phys_res:
        _, _, g0e_phys_res, _ = alpha_phys_res[0]
        f_phys = f_kin_alpha(g0e_phys_res, ALPHA_PHYS)
        print(f"\n  α=2: f(g₀^e(*)) = {f_phys:.8f}")
        print(f"       f(5/4)    = {f_kin_alpha(1.25, ALPHA_PHYS):.8f}")
        print(f"       2-f(g₀^e*)= {2-f_phys:.8f}")
        print(f"       f·ε₀      = {f_phys*0.25:.8f}")
        print(f"       f/2       = {f_phys/2:.8f}")
        # Korekta Δf = f(g₀^e(*)) - f(5/4)
        Df = f_phys - f_kin_alpha(1.25, ALPHA_PHYS)
        print(f"       Δf = f(g₀^e*) − f(5/4) = {Df:.8f}")

record("A13: Sprawdzono f(g₀^e(*)) dla różnych α",
       len(results_alpha) >= 2,
       f"N={len(results_alpha)} punktów")


# ============================================================
# SEKCJA 8: Weryfikacja ekspansji — predykcja g₀^e(*)(α)
# ============================================================
print("\n[8] WERYFIKACJA: g₀^e(*)(α) = 1 + ε − |c₂|·ε²")
print("-" * 55)
print(f"  Używając c₂_mean = {c2_mean:.6f}")
print()
print("  α      ε=1/(2α)   g₀^e(*)(num)  1+ε−|c₂|ε²   Δ_pred    Δ_rel(%)")
n_good_pred = 0
for (alpha_i, eps_i, g0e_i, c2_i) in results_alpha:
    pred_i = 1.0 + eps_i + c2_mean * eps_i**2
    delta_i = g0e_i - pred_i
    delta_rel = delta_i / eps_i * 100
    good = abs(delta_i) < 0.001
    if good:
        n_good_pred += 1
    print(f"  α={alpha_i:.1f}  ε={eps_i:.5f}  "
          f"num={g0e_i:.7f}  pred={pred_i:.7f}  Δ={delta_i:.7f}  Δ_rel={delta_rel:.4f}%")

record("A14: Predykcja 1+ε+c₂ε² dobra dla ≥2/3 wartości α (|Δ|<0.001)",
       n_good_pred >= int(len(results_alpha)*2/3),
       f"{n_good_pred}/{len(results_alpha)} predykcji dobra")


# ============================================================
# SEKCJA 9: Prosta tożsamość?
# ============================================================
print("\n[9] CZY ISTNIEJE PROSTA TOŻSAMOŚĆ?")
print("-" * 55)
print()
print("  Pytanie T-OP4 (przepisane): jeśli g₀^e(*) = 1+ε+c₂·ε²+...,")
print("  to skąd wynika, że brentq daje ε ≈ 1/(2α)?")
print()
print("  OBSERWACJA: r₂₁ = (A_μ/A_e)^4 musi równać się 9r*.")
print("  A_tail(g₀^e) to skomplikowana funkcja ODE — brak prostego wzoru.")
print()
print("  JEDNAK: sprawdźmy prostą tożsamość z ex133:")
print("  (g₀^e − g*)/(1 − g*) = 1 + 1/(2α(1−g*)) ?")
eps_phys_val = 1.0 / (2.0 * ALPHA_PHYS)
g_star_phys  = math.exp(-1.0 / (2.0 * ALPHA_PHYS))
# Lewa strona:
lhs = (G0E_STAR_EX133 - g_star_phys) / (1.0 - g_star_phys)
# Prawa strona:
rhs_identity = 1.0 + 1.0 / (2.0 * ALPHA_PHYS * (1.0 - g_star_phys))
print(f"  g* = {g_star_phys:.8f}")
print(f"  g₀^e(*) = {G0E_STAR_EX133:.8f}")
print(f"  (g₀^e*−g*)/(1−g*) = {lhs:.8f}")
print(f"  1+1/(2α(1−g*))    = {rhs_identity:.8f}")
print(f"  Różnica            = {lhs - rhs_identity:.8f}")

# To samo dla g₀^e = 5/4
lhs_54 = (1.25 - g_star_phys) / (1.0 - g_star_phys)
print(f"\n  Dla g₀^e=5/4:")
print(f"  (5/4−g*)/(1−g*) = {lhs_54:.8f}")
print(f"  1+1/(2α(1−g*)) = {rhs_identity:.8f}")
print(f"  Różnica         = {lhs_54 - rhs_identity:.8f}")
print(f"  → Ta tożsamość zachodzi dla 5/4, nie dla g₀^e(*)")

record("A15: Tożsamość (g₀^e−g*)/(1−g*)=1+1/(2α(1−g*)) zachodzi dla 5/4 (nie g₀^e(*))",
       abs(lhs_54 - rhs_identity) < abs(lhs - rhs_identity),
       f"lhs(5/4)={lhs_54:.6f}, rhs={rhs_identity:.6f}, Δ(5/4)={lhs_54-rhs_identity:.6f}")


# ============================================================
# SEKCJA 10: Czy c₂ jest powszechny — test dla innych modeli?
# ============================================================
print("\n[10] PODSUMOWANIE I WNIOSKI")
print("-" * 55)
print()
print(f"  g₀^e(*)(α=2) = 1 + 1/(2α) + c₂·(1/(2α))² + O(ε³)")
print(f"  gdzie c₂ ≈ {c2_mean:.6f}")
print()
print(f"  HIERARCHIA:")
print(f"    Rząd 0:  g₀^e = 1")
print(f"    Rząd 1:  g₀^e = 1 + 1/(2α)  [= 5/4 dla α=2]  → błąd 0.073%")
print(f"    Rząd 2:  g₀^e = 1 + 1/(2α) + c₂/(2α)²      → błąd ~0%")
print()
print(f"  FIZYCZNA INTERPRETACJA:")
print(f"    1/(2α) = −ln(g*) = wykładnik potencjału ghost wall")
print(f"    g₀^e ≈ 1 + (−ln g*) — soliton elektronu jest odległe o")
print(f"    'ładunek logarytmiczny' ściany ducha od g=1")
print()

# Sprawdź też: g₀^e(*) - 1 vs -ln(g*)
neg_ln_gstar = -math.log(g_star_phys)
print(f"  −ln(g*) = 1/(2α) = {neg_ln_gstar:.8f}")
print(f"  g₀^e(*) − 1   = {G0E_STAR_EX133 - 1:.8f}")
print(f"  Różnica        = {(G0E_STAR_EX133 - 1) - neg_ln_gstar:.8f}")
print()
print(f"  WNIOSEK: g₀^e(*) = 1 − ln(g*) + O((ln g*)²)")
print(f"  czyli soliton elektronu ma AMPLITUDĘ CENTRUM równą")
print(f"  1 + modulus ściany ducha (logarytmicznie).")

record("A16: g₀^e(*) − 1 ≈ −ln(g*) = 1/(2α) do 0.08%",
       abs((G0E_STAR_EX133-1) - neg_ln_gstar) / neg_ln_gstar < 0.001,
       f"g₀^e(*)-1={G0E_STAR_EX133-1:.6f}, -ln(g*)={neg_ln_gstar:.6f}")


# ============================================================
print("\n" + "=" * 72)
print("WYNIKI TESTÓW")
print("=" * 72)
n_pass = sum(1 for _, p, _ in TESTS if p)
n_fail = sum(1 for _, p, _ in TESTS if not p)
print(f"\n  ŁĄCZNIE: {n_pass}/{len(TESTS)} PASS  ({n_fail} FAIL)")
for name, passed, detail in TESTS:
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"         {detail}")

print()
print("=" * 72)
print("EX136 ZAKOŃCZONE")
print("=" * 72)
