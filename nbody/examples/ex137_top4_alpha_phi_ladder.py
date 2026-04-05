#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex137_top4_alpha_phi_ladder.py
===============================
BADANIE T-OP4: WŁAŚCIWY SKAN α Z φ-DRABINĄ

KONTEKST:
  ex136 używał STAŁYCH fizycznych proporcji R_MU=1.6181, R_TAU=2.5530
  dla wszystkich α — to błąd metodyczny, bo te proporcje są specyficzne
  dla α=2. Dla każdego α powinniśmy używać proporcji z właściwego
  punktu stałego φ-drabiny.

  Naturalny wybór: φ-drabina matematyczna, gdzie
      g₀^μ = φ · g₀^e,   g₀^τ = φ² · g₀^e   (φ = złoty podział)
  — ta drabina jest α-niezależna i matematycznie uzasadniona.

  Dla α=2 fizyczne proporcje to:
      R_MU_phys = 2.02117/1.24915 ≈ 1.6181 ≈ φ         (Δ ≈ 0.006%)
      R_TAU_phys = 3.18912/1.24915 ≈ 2.5530 ≠ φ²=2.618  (Δ ≈ 2.5%)

  Pytanie: czy g₀^{e,*}(α) ≈ 1+1/(2α) jest specyficzne dla α=2
  używając MATEMATYCZNEJ φ-drabiny?

PLAN (16 testów):
  B01-B04: Weryfikacja bazowa przy α=2, φ-drabina vs fizyczne proporcje
  B05-B09: Skan α ∈ {1.5, 2.0, 2.5, 3.0, 4.0, 5.0} z φ-drabiną
  B10-B13: Współczynnik c₂(α) z φ-drabiny — czy stały?
  B14-B16: Podsumowanie: wyjątkowość α=2

WYNIKI OCZEKIWANE (na podstawie ex136):
  B01: θ_TGP(φ-ladder) przy g₀^e=1.249082, α=2 bliskie 132.73° (ale nie identyczne z θ(9r*))
  B02: g₀^{e,*}(α=2, φ-ladder) bliskie 1.249082
  B04: c₂(α=2, φ-ladder) ≈ -1/68 (sprawdzenie spójności)
  B10: c₂(α) NIE jest stałe (T-OP4d)
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
# Stałe globalne
# ============================================================
ALPHA_PHYS = 2.0
PHI        = (1.0 + math.sqrt(5.0)) / 2.0   # 1.6180339...
PHI2       = PHI * PHI                         # 2.6180339...
RSTAR      = (23.0 + 5.0 * math.sqrt(21.0)) / 2.0  # ≈ 22.9564
THETA_9RS  = 132.731439   # θ(9r*) z ex129

# Fizyczne proporcje z ex133 (α=2)
R_MU_PHYS  = 2.02117 / 1.24915   # ≈ 1.6181 ≈ φ
R_TAU_PHYS = 3.18912 / 1.24915   # ≈ 2.5530 ≠ φ²

# Integracja ODE
R_MAX    = 70.0
R_START  = 1e-4
MAX_STEP = 0.02
RTOL     = 1e-10
ATOL     = 1e-13

# Okno fitu amplitudy ogona
R_L_FIT  = 20.0
R_R_FIT  = 40.0

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


def integrate_soliton_alpha(g0, alpha, r_max=None, max_bounces=30):
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


def fit_atail(r_arr, g_arr, r_L=None, r_R=None):
    if r_L is None: r_L = R_L_FIT
    if r_R is None: r_R = R_R_FIT
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


# ============================================================
# Kąt Brannena (wzór Fourierowski, identyczny z ex133)
# ============================================================
def brannen_theta_direct(m1, m2, m3):
    sqm = np.array([math.sqrt(m1), math.sqrt(m2), math.sqrt(m3)])
    M   = float(np.mean(sqm))
    eps = sqm / M - 1.0
    F1  = sum(eps[k] * cmath.exp(-2j * math.pi * k / 3) for k in range(3))
    return math.degrees(math.atan2(F1.imag, F1.real))


# ============================================================
# θ_TGP z φ-drabiną: g₀^μ=φ·g₀^e, g₀^τ=φ²·g₀^e
# ============================================================
def theta_tgp_phi_ladder(g0e, alpha):
    """
    Wariant φ-drabiny matematycznej:
      g₀^μ = φ · g₀^e
      g₀^τ = φ² · g₀^e
    Amplitudy ogona → masy proxy → kąt Brannena.
    """
    g0mu  = PHI  * g0e
    g0tau = PHI2 * g0e
    Ae   = atail_alpha(g0e,   alpha)
    Amu  = atail_alpha(g0mu,  alpha)
    Atau = atail_alpha(g0tau, alpha)
    if any((not math.isfinite(A)) or A <= 0 for A in [Ae, Amu, Atau]):
        return float('nan')
    return brannen_theta_direct(Ae**4, Amu**4, Atau**4)


# ============================================================
# θ_TGP z fizycznymi proporcjami (ex136 v2)
# ============================================================
def theta_tgp_phys_ratios(g0e, alpha, r_mu=R_MU_PHYS, r_tau=R_TAU_PHYS):
    g0mu  = r_mu  * g0e
    g0tau = r_tau * g0e
    Ae   = atail_alpha(g0e,   alpha)
    Amu  = atail_alpha(g0mu,  alpha)
    Atau = atail_alpha(g0tau, alpha)
    if any((not math.isfinite(A)) or A <= 0 for A in [Ae, Amu, Atau]):
        return float('nan')
    return brannen_theta_direct(Ae**4, Amu**4, Atau**4)


# ============================================================
# Wyznaczanie g₀^{e,*}(α) via brentq
# ============================================================
def find_g0e_star(alpha, theta_target=THETA_9RS, mode='phi',
                  g0_lo=1.05, g0_hi=2.0, xtol=1e-8, verbose=False):
    """
    Szuka g₀^{e,*} tak, że θ_TGP(g₀^e; α) = theta_target.
    mode='phi'  → φ-drabina matematyczna
    mode='phys' → fizyczne proporcje (R_MU_PHYS, R_TAU_PHYS)
    Zwraca (g0e_star, theta_found) lub (nan, nan) w razie błędu.
    """
    if mode == 'phi':
        theta_fn = lambda g: theta_tgp_phi_ladder(g, alpha)
    else:
        theta_fn = lambda g: theta_tgp_phys_ratios(g, alpha)

    # Sprawdź znak funkcji na brzegach
    try:
        f_lo = theta_fn(g0_lo) - theta_target
        f_hi = theta_fn(g0_hi) - theta_target
    except Exception:
        return float('nan'), float('nan')

    if verbose:
        print(f"    brentq bracket: f({g0_lo:.3f})={f_lo:+.4f}, f({g0_hi:.3f})={f_hi:+.4f}")

    if not (math.isfinite(f_lo) and math.isfinite(f_hi)):
        return float('nan'), float('nan')

    if f_lo * f_hi > 0:
        # Próba szerszego zakresu
        for g_test in np.linspace(g0_lo, g0_hi, 20):
            ft = theta_fn(g_test) - theta_target
            if math.isfinite(ft):
                if ft * f_lo < 0:
                    g0_hi, f_hi = g_test, ft
                    break
                elif ft * f_hi < 0:
                    g0_lo, f_lo = g_test, ft
                    break
        else:
            return float('nan'), float('nan')

    try:
        g0e_star = brentq(lambda g: theta_fn(g) - theta_target,
                          g0_lo, g0_hi, xtol=xtol, maxiter=100)
        theta_found = theta_fn(g0e_star)
        return g0e_star, theta_found
    except Exception as e:
        if verbose:
            print(f"    brentq error: {e}")
        return float('nan'), float('nan')


# ============================================================
# BLOK WERYFIKACJI BAZOWEJ (α=2)
# ============================================================
print("=" * 65)
print("ex137: WŁAŚCIWY SKAN α Z φ-DRABINĄ")
print(f"  PHI   = {PHI:.10f}")
print(f"  PHI²  = {PHI2:.10f}")
print(f"  R_MU_phys  = {R_MU_PHYS:.10f}  (z ex133)")
print(f"  R_TAU_phys = {R_TAU_PHYS:.10f}  (z ex133)")
print(f"  Δ_MU  = PHI - R_MU_phys  = {PHI - R_MU_PHYS:+.6f}")
print(f"  Δ_TAU = PHI² - R_TAU_phys = {PHI2 - R_TAU_PHYS:+.6f}")
print("=" * 65)

print("\n--- BLOK 1: WERYFIKACJA BAZOWA (α=2) ---")

alpha2 = 2.0
eps2   = 1.0 / (2.0 * alpha2)   # = 0.25

# B01: θ_TGP(φ-ladder) przy g₀^e=1.249082, α=2
G0E_REF = 1.249082
theta_phi_ref = theta_tgp_phi_ladder(G0E_REF, alpha2)
print(f"\n  θ_TGP(φ-ladder, g₀^e=1.249082, α=2) = {theta_phi_ref:.6f}°")
print(f"  θ(9r*)  = {THETA_9RS:.6f}°")
print(f"  Δθ      = {theta_phi_ref - THETA_9RS:+.6f}°")
b01_pass = math.isfinite(theta_phi_ref) and abs(theta_phi_ref - THETA_9RS) < 0.05
record("B01", b01_pass,
       f"θ_TGP(φ-ladder,g₀^e=1.249082)={theta_phi_ref:.4f}° vs {THETA_9RS}° (Δ={theta_phi_ref-THETA_9RS:+.4f}°)")

# B02: g₀^{e,*}(α=2) z φ-drabiną
print("\n  Szukam g₀^{e,*}(α=2, φ-drabina)...")
g0e_star_phi2, theta_phi2 = find_g0e_star(alpha2, mode='phi', verbose=True)
print(f"  g₀^{{e,*}}(φ-ladder, α=2) = {g0e_star_phi2:.8f}")
print(f"  θ_TGP  = {theta_phi2:.6f}°  (cel: {THETA_9RS}°)")

b02_pass = math.isfinite(g0e_star_phi2) and abs(g0e_star_phi2 - 1.249082) < 0.05
record("B02", b02_pass,
       f"g₀^e*(φ-ladder,α=2)={g0e_star_phi2:.6f}, Δ vs ex136={g0e_star_phi2-1.2490816636:+.6f}")

# B03: |g₀^{e,*} - 5/4| < 0.01
b03_pass = math.isfinite(g0e_star_phi2) and abs(g0e_star_phi2 - 1.25) < 0.01
record("B03", b03_pass,
       f"|g₀^e*(φ)−5/4| = {abs(g0e_star_phi2-1.25):.6f} (próg 0.01)")

# B04: c₂(α=2) z φ-drabiny
if math.isfinite(g0e_star_phi2):
    c2_phi2 = (g0e_star_phi2 - 1.0 - eps2) / eps2**2
    c2_candidate = -1.0 / 68.0
    c2_err = abs(c2_phi2 - c2_candidate) / abs(c2_candidate)
    print(f"\n  c₂(α=2, φ-ladder) = {c2_phi2:.6f}")
    print(f"  Kandydat -1/68 = {c2_candidate:.6f},  err = {c2_err*100:.2f}%")
    b04_pass = c2_err < 0.02   # 2%
    record("B04", b04_pass,
           f"c₂={c2_phi2:.6f} vs -1/68={c2_candidate:.6f} (err={c2_err*100:.2f}%)")
else:
    c2_phi2 = float('nan')
    record("B04", False, "g₀^e* nie wyznaczony — c₂ niedostępne")

# B05: θ_TGP(fizyczne proporcje) vs θ_TGP(φ-drabina) przy α=2, g₀^e*
g0e_star_phys2, theta_phys2 = find_g0e_star(alpha2, mode='phys', verbose=False)
print(f"\n  Porównanie metod (α=2):")
print(f"  g₀^e*(φ-ladder) = {g0e_star_phi2:.8f}")
print(f"  g₀^e*(phys)     = {g0e_star_phys2:.8f}  (ex136 wynik)")
if math.isfinite(g0e_star_phi2) and math.isfinite(g0e_star_phys2):
    diff_methods = abs(g0e_star_phi2 - g0e_star_phys2)
    rel_diff     = diff_methods / g0e_star_phys2
    print(f"  |Δg₀^e*|        = {diff_methods:.6f}  (rel. {rel_diff*100:.3f}%)")
    b05_pass = rel_diff < 0.05   # 5%
    record("B05", b05_pass,
           f"|φ-ladder − phys-ratios| = {diff_methods:.6f} ({rel_diff*100:.3f}%)")
else:
    record("B05", False, "Jedna z metod zwróciła NaN")


# ============================================================
# BLOK 2: SKAN α Z φ-DRABINĄ
# ============================================================
print("\n--- BLOK 2: SKAN α Z φ-DRABINĄ ---")

ALPHA_SCAN = [1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0]
results = {}   # {alpha: (g0e_star, theta_found, c2, eps, delta_g)}

print(f"\n  {'α':>5}  {'g₀^e*(φ)':>12}  {'1+1/(2α)':>10}  {'Δg':>8}  "
      f"{'c₂':>9}  {'θ_TGP':>10}")
print("  " + "-" * 65)

for alpha in ALPHA_SCAN:
    eps_a = 1.0 / (2.0 * alpha)
    g0e_s, theta_s = find_g0e_star(alpha, mode='phi',
                                    g0_lo=1.02, g0_hi=2.5)
    if math.isfinite(g0e_s):
        approx_a = 1.0 + eps_a
        delta_g  = g0e_s - approx_a
        c2_a     = (g0e_s - 1.0 - eps_a) / eps_a**2
    else:
        approx_a = 1.0 + eps_a
        delta_g  = float('nan')
        c2_a     = float('nan')
    results[alpha] = (g0e_s, theta_s, c2_a, eps_a, delta_g)
    print(f"  {alpha:>5.1f}  {g0e_s:>12.8f}  {approx_a:>10.6f}  {delta_g:>+8.5f}"
          f"  {c2_a:>9.5f}  {theta_s:>10.4f}°")

# B06: g₀^e*(α=1.5) wyznaczony
r15 = results.get(1.5)
b06_pass = r15 is not None and math.isfinite(r15[0])
record("B06", b06_pass, f"α=1.5: g₀^e*={r15[0]:.6f}" if b06_pass else "NaN")

# B07: g₀^e*(α=2.5) wyznaczony
r25 = results.get(2.5)
b07_pass = r25 is not None and math.isfinite(r25[0])
record("B07", b07_pass, f"α=2.5: g₀^e*={r25[0]:.6f}" if b07_pass else "NaN")

# B08: g₀^e*(α=3.0) wyznaczony
r30 = results.get(3.0)
b08_pass = r30 is not None and math.isfinite(r30[0])
record("B08", b08_pass, f"α=3.0: g₀^e*={r30[0]:.6f}" if b08_pass else "NaN")

# B09: g₀^e*(α=4.0) wyznaczony
r40 = results.get(4.0)
b09_pass = r40 is not None and math.isfinite(r40[0])
record("B09", b09_pass, f"α=4.0: g₀^e*={r40[0]:.6f}" if b09_pass else "NaN")


# ============================================================
# BLOK 3: ANALIZA c₂(α)
# ============================================================
print("\n--- BLOK 3: ANALIZA c₂(α) ---")

valid_alphas = [a for a in ALPHA_SCAN if math.isfinite(results[a][2])]
c2_vals      = [results[a][2] for a in valid_alphas]
eps_vals     = [results[a][3] for a in valid_alphas]

print(f"\n  Zebrane c₂(α) dla {len(valid_alphas)} wartości α:")
for a, c2 in zip(valid_alphas, c2_vals):
    eps_a = 1.0 / (2.0 * a)
    print(f"    α={a:.1f}: c₂={c2:+.5f},  ε={eps_a:.4f},  "
          f"1+ε+c₂ε² = {1+eps_a+c2*eps_a**2:.6f}  "
          f"(numeryk: {results[a][0]:.6f})")

if len(c2_vals) >= 3:
    c2_mean = float(np.mean(c2_vals))
    c2_std  = float(np.std(c2_vals))
    c2_max_dev = max(abs(c - c2_mean) for c in c2_vals)
    print(f"\n  Średnia c₂ = {c2_mean:.5f} ± {c2_std:.5f}")
    print(f"  Max. odchylenie = {c2_max_dev:.5f}")

    # B10: c₂ NIE jest stałe (std > 0.01)
    b10_pass = c2_std > 0.01
    record("B10", b10_pass,
           f"c₂ zmienna: std={c2_std:.5f} {'> 0.01' if b10_pass else '≤ 0.01 (stałe?)'}")

    # B11: c₂(α=2) jest ujemne
    c2_at_2 = results[2.0][2]
    b11_pass = math.isfinite(c2_at_2) and c2_at_2 < 0
    record("B11", b11_pass, f"c₂(α=2)={c2_at_2:.6f} {'< 0 (OK)' if b11_pass else '≥ 0 (FAIL)'}")

    # B12: Δg = g₀^e* - (1+ε) ma różne znaki dla różnych α
    delta_g_vals = [results[a][4] for a in valid_alphas if math.isfinite(results[a][4])]
    signs = [1 if d > 0 else -1 for d in delta_g_vals]
    sign_changes = sum(1 for i in range(len(signs)-1) if signs[i] != signs[i+1])
    b12_pass = sign_changes >= 1
    record("B12", b12_pass,
           f"Zmiana znaku Δg: {sign_changes} {'≥1 (α=2 wyjątkowe)' if b12_pass else '=0 (monoton)'}")
    print(f"\n  Δg = g₀^e*(α) - [1+1/(2α)]:")
    for a, dg in zip(valid_alphas, [results[a][4] for a in valid_alphas]):
        print(f"    α={a:.1f}: Δg={dg:+.6f}")

else:
    record("B10", False, "Niewystarczająca liczba punktów")
    record("B11", False, "Niewystarczająca liczba punktów")
    record("B12", False, "Niewystarczająca liczba punktów")


# ============================================================
# BLOK 4: DOPASOWANIE POTĘGOWE I KANDYDACI ALGEBRAICZNI
# ============================================================
print("\n--- BLOK 4: KANDYDACI ALGEBRAICZNI ---")

# B13: Fit liniowy log(g₀^e*-1) vs log(ε)
if len(valid_alphas) >= 3:
    log_eps   = [math.log(results[a][3]) for a in valid_alphas
                 if results[a][0] > 1.001]
    log_g1    = [math.log(results[a][0] - 1.0) for a in valid_alphas
                 if results[a][0] > 1.001]
    if len(log_eps) >= 3:
        slope, intercept, r_val, _, _ = linregress(log_eps, log_g1)
        print(f"\n  Fit potęgowy (g₀^e*-1) ∝ ε^n:")
        print(f"    wykładnik n = {slope:.4f}  (oczekiwany: 1.0)")
        print(f"    prefaktor   = {math.exp(intercept):.6f}  (oczekiwany: 1.0)")
        print(f"    r²          = {r_val**2:.6f}")
        b13_pass = abs(slope - 1.0) < 0.3 and r_val**2 > 0.90
        record("B13", b13_pass,
               f"n={slope:.4f} (oczekiwany 1.0±0.3), r²={r_val**2:.4f}")
    else:
        record("B13", False, "Niewystarczające dane")
else:
    record("B13", False, "Niewystarczające dane")

# B14: Sprawdzenie kandydata c₂ ≈ -1/68 tylko dla α=2
c2_a2 = results[2.0][2]
if math.isfinite(c2_a2):
    err_68 = abs(c2_a2 - (-1.0/68.0)) / abs(-1.0/68.0) * 100
    print(f"\n  c₂(α=2) = {c2_a2:.6f}")
    print(f"  -1/68   = {-1/68.0:.6f}")
    print(f"  Błąd    = {err_68:.3f}%")
    b14_pass = err_68 < 5.0   # 5% dopasowanie
    record("B14", b14_pass,
           f"c₂(α=2)={c2_a2:.6f} vs -1/68 (err={err_68:.2f}%)")
else:
    record("B14", False, "c₂(α=2) nie wyznaczone")

# B15: g₀^e*(α=2) bliskie 1+ε = 5/4 (Δ < 0.1%)
if math.isfinite(g0e_star_phi2):
    rel_err_54 = abs(g0e_star_phi2 - 1.25) / 1.25 * 100
    print(f"\n  g₀^e*(φ,α=2)  = {g0e_star_phi2:.8f}")
    print(f"  5/4 = 1+1/(2α) = {1.25:.8f}")
    print(f"  Różnica        = {g0e_star_phi2-1.25:+.8f}  ({rel_err_54:.4f}%)")
    b15_pass = rel_err_54 < 0.5
    record("B15", b15_pass,
           f"|g₀^e*(φ,2) - 5/4| = {rel_err_54:.4f}% {'< 0.5%' if b15_pass else '≥ 0.5%'}")
else:
    record("B15", False, "g₀^e* NaN")


# ============================================================
# BLOK 5: PODSUMOWANIE I T-OP4
# ============================================================
print("\n--- BLOK 5: PODSUMOWANIE T-OP4 ---")

# B16: Zbieżność g₀^e* ≈ 1+1/(2α) specyfyczna dla α=2 vs inne α
if len(valid_alphas) >= 4:
    rel_errs = {}
    for a in valid_alphas:
        g0s = results[a][0]
        if math.isfinite(g0s):
            rel_errs[a] = abs(g0s - (1 + 1/(2*a))) / (1 + 1/(2*a)) * 100

    print(f"\n  Względny błąd |g₀^e*(α) - (1+ε)| / (1+ε):")
    for a, err in sorted(rel_errs.items()):
        print(f"    α={a:.1f}: {err:.4f}%")

    if 2.0 in rel_errs:
        err_at_2 = rel_errs[2.0]
        other_errs = [e for a, e in rel_errs.items() if a != 2.0]
        min_other  = min(other_errs) if other_errs else float('inf')
        b16_pass   = err_at_2 < min_other   # α=2 ma najmniejszy błąd
        record("B16", b16_pass,
               f"α=2 ma err={err_at_2:.3f}% < min_other={min_other:.3f}% "
               f"→ {'wyjątkowe' if b16_pass else 'NIE wyjątkowe'}")
    else:
        record("B16", False, "Brak wyniku dla α=2")
else:
    record("B16", False, "Niewystarczające dane")


# ============================================================
# TABELA WYNIKOWA
# ============================================================
print("\n" + "=" * 65)
print("TABELA KOŃCOWA WYNIKÓW:")
print("=" * 65)
print(f"\n  {'α':>4}  {'g₀^e*(φ-drabin)':>17}  {'1+1/(2α)':>10}  {'Δg':>8}  {'c₂':>8}")
print("  " + "-"*60)
for alpha in ALPHA_SCAN:
    g0s, _, c2v, epsa, dg = results[alpha]
    approx = 1 + epsa
    print(f"  {alpha:>4.1f}  {g0s:>17.8f}  {approx:>10.6f}  {dg:>+8.5f}  {c2v:>8.4f}")

print(f"\n  Stałe:")
print(f"    PHI  = {PHI:.8f}  (= g₀^μ/g₀^e)")
print(f"    PHI² = {PHI2:.8f}  (= g₀^τ/g₀^e)")
print(f"    r*   = {RSTAR:.8f}  (stały, niezależny od α)")
print(f"    θ(9r*) = {THETA_9RS:.6f}°")

# ============================================================
# WYNIKI TESTÓW
# ============================================================
print("\n" + "=" * 65)
passed = sum(1 for _, p, _ in TESTS if p)
total  = len(TESTS)
print(f"WYNIKI: {passed}/{total} PASS")
print("=" * 65)
for name, p, detail in TESTS:
    mark = "PASS" if p else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"         {detail}")

# ============================================================
# WNIOSKI
# ============================================================
print("\n--- WNIOSKI ---")
print(f"""
  1. METODYCZNA KOREKTA EX136:
     Fizyczne proporcje ex136 (R_TAU=2.553) różnią się od φ²=2.618
     o {abs(PHI2 - R_TAU_PHYS)/PHI2*100:.1f}%. Właściwa φ-drabina matematyczna
     daje inne g₀^e*(α) dla α≠2.

  2. WYNIK DLA α=2:
     g₀^e*(φ-drabina, α=2) = {g0e_star_phi2:.8f}
     Zbieżność z 5/4: {abs(g0e_star_phi2-1.25)/1.25*100:.4f}%
     c₂ = {c2_phi2:.6f} ≈ -1/68 = {-1/68:.6f}

  3. ODPOWIEDŹ NA T-OP4:
     Zbliżenie g₀^e* ≈ 1+1/(2α) {"jest" if (2.0 in rel_errs and rel_errs[2.0] < min([e for a,e in rel_errs.items() if a!=2.0], default=999)) else "NIE jest"}
     specyficzne dla α=2 (φ-drabina matematyczna).

  4. T-OP4d OTWARTE:
     c₂(α) zmienia się z α → nie ma prostego wyrażenia c₂=const.
     Dla α=2: c₂ ≈ {c2_phi2:.4f} ≈ -1/68.
""")
