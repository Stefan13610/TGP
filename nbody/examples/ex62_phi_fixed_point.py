"""
EX62: phi-SELF-CONSISTENT FIXED POINT I POSZUKIWANIE alpha*
============================================================
STATUS: LEGACY-TRANSLATIONAL

This file is historically important for the route to `phi-FP`, but it still
uses the pre-canonical `g0* ~ 1.24` / `alpha*` framing. For the synchronized
current formulation prefer `ex195`–`ex197`, `ex205`, and `examples/STATUS_MAP.md`.

Cel: rozwiązanie krytycznej diagnozy z ex61:
  (A(phi*z0)/A(z0))^4 = 252  !=  207   przy alpha=2

Pytania:
  1. Czy istnieje g0* ∈ n=0 (sektor: 1 odbicie) takie, że
       (A(phi*g0*) / A(g0*))^4 = 207?
     (phi-self-consistent fixed point)

  2. Przy jakim alpha* (skan 1.5→3.5) i z0(alpha*)
       (A(phi*z0(alpha*)) / A(z0(alpha*)))^4 = 207?

  3. Jak wygląda pełna mapa R(g0) = (A(phi*g0)/A(g0))^4 vs g0?
     Szukamy przecięcia z linią R=207.

  4. Porównanie z danymi eksperymentalnymi:
     m_mu/m_e = 206.768
     m_tau/m_e = 3477.0
     Czy istnieje g0** takie, że (A(phi^2*g0**)/A(g0**))^4 = 3477?

  5. Faza ogona delta_scatt(g0) = atan2(B,C) — wykres ciągły
     i warunek kwantowania Bohr-Sommerfeld: delta_scatt = n*pi.

Metoda:
  - A_inf extrapolation (jak w ex60/ex61)
  - scipy.optimize.brentq dla precyzyjnych zer
  - alpha skan z nową funkcją make_ode_functions(alpha)

Testy:
  T1: g0* (FP) ∈ [1.20, 1.30]  (oczekiwany: ≈1.24?)
  T2: |R(g0*) - 207| / 207 < 0.005   (precyzja 0.5%)
  T3: alpha* ∈ [2.0, 3.0]   (interpolacja z ex61: ~2.3-2.4)
  T4: |(A(phi^2*g0**)/A(g0**))^4 - 3477| / 3477 < 0.3   (tauon 30%)
  T5: delta_scatt(z0) ≈ 0  i  delta_scatt(phi*z0) ≈ pi  (± 0.5 rad)
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
import warnings
warnings.filterwarnings('ignore')

# ===== STAŁE =====
PHI = (1.0 + np.sqrt(5.0)) / 2.0   # złota proporcja
R21_EXP = 206.768                   # m_mu / m_e (eksperyment)
R31_EXP = 3477.0                    # m_tau / m_e (eksperyment)
R_MAX = 100.0
N_EVAL = 4000
WIN_LIST = [16, 22, 28, 36, 46, 58, 72]
WIN_WIDTH = 14.0

# ===== PARAMETRYZOWANE ODE =====

def make_ode_functions(alpha, g_bounce_offset=0.005):
    """Zwraca rhs i event dla danego alpha."""
    g_ghost  = np.exp(-1.0 / (2.0 * alpha))
    g_bounce = g_ghost + g_bounce_offset

    def rhs(r, y):
        g, gp = y
        g = max(g, g_bounce + 1e-7)
        fg = 1.0 + 2.0 * alpha * np.log(g)
        if abs(fg) < 1e-10:
            return [gp, 0.0]
        dr = g**2 * (1.0 - g)
        cr = (alpha / g) * gp**2
        if r < 1e-10:
            return [gp, (dr - cr) / (3.0 * fg)]
        return [gp, (dr - cr - fg * 2.0 * gp / r) / fg]

    def ev_ghost(r, y):
        return y[0] - g_bounce
    ev_ghost.terminal  = True
    ev_ghost.direction = -1

    return rhs, ev_ghost, g_ghost, g_bounce


def integrate_with_bounces(g0, alpha=2.0, g_bounce_offset=0.005,
                            max_bounces=8, r_max=R_MAX, n_eval=N_EVAL):
    """Integracja ODE z elastycznym odbiciem przy g*."""
    rhs, ev_ghost, _, g_bounce = make_ode_functions(alpha, g_bounce_offset)
    r_pts = np.linspace(0.0, r_max, n_eval)
    y0 = [g0, 0.0]
    r_cur = 0.0
    r_all, g_all = [], []
    n_bounces = 0

    for _ in range(max_bounces + 1):
        r_span = (r_cur, r_max)
        sol = solve_ivp(rhs, r_span, y0, events=ev_ghost,
                        dense_output=True, rtol=1e-9, atol=1e-11,
                        max_step=0.05)
        r_seg = sol.t; g_seg = sol.y[0]
        r_all.append(r_seg); g_all.append(g_seg)

        if sol.status == 1 and len(sol.t_events[0]) > 0:
            r_hit = sol.t_events[0][-1]
            state  = sol.sol(r_hit)
            y0 = [state[0], -state[1]]   # odwracamy gp
            r_cur = r_hit
            n_bounces += 1
        else:
            break

    r_arr = np.concatenate(r_all)
    g_arr = np.concatenate(g_all)
    idx   = np.argsort(r_arr)
    return r_arr[idx], g_arr[idx], n_bounces


def fit_tail_window(r_arr, g_arr, rL, rR):
    """Dopasowanie (B*cos + C*sin)/r w oknie [rL, rR]."""
    mask = (r_arr >= rL) & (r_arr <= rR)
    if np.sum(mask) < 20:
        return 0.0, 0.0, 0.0
    r_f  = r_arr[mask]
    d_f  = (g_arr[mask] - 1.0) * r_f
    M    = np.column_stack([np.cos(r_f), np.sin(r_f)])
    res  = np.linalg.lstsq(M, d_f, rcond=None)
    B, C = res[0]
    A    = np.sqrt(B**2 + C**2)
    return float(A), float(B), float(C)


def A_infinity(g0, alpha=2.0, g_bounce_offset=0.005):
    """A_inf z dopasowania A(rL) = A_inf*(1 + a/rL)."""
    r, g, _ = integrate_with_bounces(g0, alpha=alpha,
                                      g_bounce_offset=g_bounce_offset)
    A_vals, rL_vals = [], []
    for rL in WIN_LIST:
        rR = rL + WIN_WIDTH
        if rR > r[-1]:
            break
        A, _, _ = fit_tail_window(r, g, rL, rR)
        if A > 0.002:
            A_vals.append(A)
            rL_vals.append(rL)
    if len(A_vals) < 3:
        return A_vals[-1] if A_vals else 0.0
    try:
        p, _ = curve_fit(lambda x, ai, a: ai * (1.0 + a / x),
                         rL_vals, A_vals,
                         p0=[A_vals[-1], 0.0],
                         maxfev=2000)
        return float(p[0])
    except Exception:
        return float(A_vals[-1])


def ratio_phi(g0_e, alpha=2.0):
    """Zwraca (A_inf(phi*g0_e) / A_inf(g0_e))^4."""
    g0_mu = PHI * g0_e
    Ae = A_infinity(g0_e,  alpha=alpha)
    Am = A_infinity(g0_mu, alpha=alpha)
    if Ae < 1e-6:
        return np.nan
    return (Am / Ae) ** 4


def ratio_phi2(g0_e, alpha=2.0):
    """Zwraca (A_inf(phi^2*g0_e) / A_inf(g0_e))^4  (tauon)."""
    g0_tau = PHI**2 * g0_e
    Ae  = A_infinity(g0_e,  alpha=alpha)
    Atau = A_infinity(g0_tau, alpha=alpha)
    if Ae < 1e-6:
        return np.nan
    return (Atau / Ae) ** 4


def find_B_tail_zero(alpha=2.0, g0_lo=1.10, g0_hi=1.50,
                     g_bounce_offset=0.005, tol=1e-5):
    """Pierwsze zero B_tail (warunek H1) dla danego alpha."""
    def B_of(g0):
        r, g, _ = integrate_with_bounces(g0, alpha=alpha,
                                          g_bounce_offset=g_bounce_offset)
        A, B, C = fit_tail_window(r, g, WIN_LIST[1], WIN_LIST[1] + WIN_WIDTH)
        return B

    # Skan wstępny
    g_scan = np.linspace(g0_lo, g0_hi, 60)
    B_scan = [B_of(g) for g in g_scan]
    z0 = None
    for i in range(len(B_scan) - 1):
        if B_scan[i] * B_scan[i+1] < 0:
            try:
                z0 = brentq(B_of, g_scan[i], g_scan[i+1], xtol=tol)
            except Exception:
                pass
            break
    return z0


# ============================================================
print("=" * 70)
print("EX62: phi-SELF-CONSISTENT FIXED POINT I POSZUKIWANIE alpha*")
print("=" * 70)
print(f"  phi = {PHI:.7f}")
print(f"  R_MAX = {R_MAX:.1f},  R21_exp = {R21_EXP}")
print()

# ============================================================
# CZĘŚĆ 1: Mapa R(g0) = (A(phi*g0)/A(g0))^4  dla g0 ∈ [1.10, 1.45]
# ============================================================
print("-" * 70)
print("--- Część 1: Mapa R(g0) = (A(phi*g0)/A(g0))^4  vs  g0 ---")
print()

N1 = 50
g0_scan1 = np.linspace(1.10, 1.45, N1)
R_scan1  = []
for g0 in g0_scan1:
    R_scan1.append(ratio_phi(g0))
    print(f"  g0={g0:.4f}  R={R_scan1[-1]:.2f}", flush=True)

R_arr1 = np.array(R_scan1)
print()
print(f"  Zakres R: [{np.nanmin(R_arr1):.1f}, {np.nanmax(R_arr1):.1f}]")
print(f"  Szukamy przecięcia z R=207 ...")

# Znajdź zerowanie R - 207
fp_found = []
for i in range(len(R_arr1) - 1):
    if np.isnan(R_arr1[i]) or np.isnan(R_arr1[i+1]):
        continue
    if (R_arr1[i] - R21_EXP) * (R_arr1[i+1] - R21_EXP) < 0:
        try:
            g0_fp = brentq(lambda x: ratio_phi(x) - R21_EXP,
                           g0_scan1[i], g0_scan1[i+1], xtol=1e-4)
            R_fp  = ratio_phi(g0_fp)
            fp_found.append((g0_fp, R_fp))
            print(f"  ⭐ FIXED POINT: g0* = {g0_fp:.6f},  R = {R_fp:.3f}  (vs 207)")
        except Exception as e:
            print(f"  brentq error: {e}")

print()

# ============================================================
# CZĘŚĆ 2: Precyzyjny fixed point — A_inf obu punktów
# ============================================================
print("-" * 70)
print("--- Część 2: Precyzyjny fixed point (A_inf, brentq) ---")
print()

if fp_found:
    g0_star, R_star = fp_found[0]
    g0_mu_star = PHI * g0_star
    A_e_star   = A_infinity(g0_star)
    A_mu_star  = A_infinity(g0_mu_star)
    ratio_star = (A_mu_star / A_e_star)**4

    print(f"  g0* (elektron FP) = {g0_star:.6f}")
    print(f"  phi*g0* (mion FP) = {g0_mu_star:.6f}")
    print(f"  A_inf(g0*)        = {A_e_star:.6f}")
    print(f"  A_inf(phi*g0*)    = {A_mu_star:.6f}")
    print(f"  (A_mu/A_e)^4      = {ratio_star:.4f}  (vs R21_exp = {R21_EXP})")
    print(f"  Odchylenie        = {abs(ratio_star - R21_EXP)/R21_EXP*100:.4f}%")
    print()

    # Porównanie z z0 i z 1.24
    z0_ref = 1.2301
    print(f"  Porównanie:")
    print(f"    g0* (FP)   = {g0_star:.6f}   (phi-self-consistent)")
    print(f"    z0 (H1)    = {z0_ref:.6f}   (pierwsze zero B_tail)")
    print(f"    g0_e (eks) = 1.240000   (parametr ex57)")
    print(f"    diff(FP-z0)   = {g0_star - z0_ref:.6f}")
    print(f"    diff(FP-1.24) = {g0_star - 1.24:.6f}")
else:
    print("  BRAK fixed point w zakresie [1.10, 1.45]!")
    print("  Sprawdź zakres — R może być monotoniczne.")
    g0_star = None

print()

# ============================================================
# CZĘŚĆ 3: Skan alpha* — gdzie czyste phi daje 207
# ============================================================
print("-" * 70)
print("--- Część 3: Skan alpha* (czyste phi z z0(alpha)) ---")
print()
print(f"  {'alpha':>7}  {'z0':>8}  {'phi*z0':>8}  {'R=(A_mu/A_z0)^4':>17}  {'delta%':>8}")
print("  " + "-" * 60)

alpha_list = [1.6, 1.8, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 3.0]
alpha_R    = []

for alpha in alpha_list:
    # Znajdź z0 dla tego alpha
    z0 = find_B_tail_zero(alpha=alpha)
    if z0 is None:
        print(f"  {alpha:7.2f}  z0=None  (brak zera B_tail)")
        alpha_R.append((alpha, None, None, None))
        continue
    g0_mu = PHI * z0
    Ae    = A_infinity(z0, alpha=alpha)
    Am    = A_infinity(g0_mu, alpha=alpha)
    if Ae < 1e-6:
        R = np.nan
    else:
        R = (Am / Ae)**4
    delta = (R - R21_EXP) / R21_EXP * 100 if not np.isnan(R) else np.nan
    print(f"  {alpha:7.3f}  {z0:8.5f}  {g0_mu:8.5f}  {R:17.3f}  {delta:8.3f}%")
    alpha_R.append((alpha, z0, g0_mu, R))

print()

# Znajdź alpha* przez interpolację
alpha_R_clean = [(a, r) for a, z0, mu, r in alpha_R
                 if r is not None and not np.isnan(r)]
alpha_vals = np.array([a for a, r in alpha_R_clean])
R_vals     = np.array([r for a, r in alpha_R_clean])

alpha_star = None
for i in range(len(alpha_vals) - 1):
    if (R_vals[i] - R21_EXP) * (R_vals[i+1] - R21_EXP) < 0:
        # Interpolacja liniowa
        f = (R21_EXP - R_vals[i]) / (R_vals[i+1] - R_vals[i])
        alpha_star = alpha_vals[i] + f * (alpha_vals[i+1] - alpha_vals[i])
        print(f"  ⭐ alpha* ≈ {alpha_star:.4f}  (interpolacja liniowa)")
        break

if alpha_star:
    # Precyzyjne z0(alpha*) i ratio
    z0_star = find_B_tail_zero(alpha=alpha_star)
    if z0_star:
        g0_mu_star2 = PHI * z0_star
        Ae2 = A_infinity(z0_star, alpha=alpha_star)
        Am2 = A_infinity(g0_mu_star2, alpha=alpha_star)
        R2  = (Am2 / Ae2)**4 if Ae2 > 1e-6 else np.nan
        print(f"  Weryfikacja przy alpha*={alpha_star:.4f}:")
        print(f"    z0(alpha*) = {z0_star:.6f}")
        print(f"    phi*z0     = {g0_mu_star2:.6f}")
        print(f"    R          = {R2:.4f}  (vs 207)")
        print(f"    Odchylenie = {abs(R2 - R21_EXP)/R21_EXP*100:.4f}%")
print()

# ============================================================
# CZĘŚĆ 4: Tauon — (A(phi^2*g0_e)/A(g0_e))^4 = 3477?
# ============================================================
print("-" * 70)
print("--- Część 4: Tauon — szukanie g0** takie, że ratio_phi2 = 3477 ---")
print()

# Szybki skan tauon ratio
print("  Skan ratio_phi2 = (A(phi^2*g0)/A(g0))^4 dla g0 ∈ [1.10, 1.40]:")
N4 = 25
g0_scan4 = np.linspace(1.10, 1.40, N4)
R2_scan4 = []
for g0 in g0_scan4:
    r2 = ratio_phi2(g0)
    R2_scan4.append(r2)

print(f"  {'g0':>8}  {'ratio_phi^2':>14}")
for g0, r2 in zip(g0_scan4, R2_scan4):
    if not np.isnan(r2):
        print(f"  {g0:8.4f}  {r2:14.1f}")

# Szukaj zerowania
tau_fp = []
for i in range(len(R2_scan4) - 1):
    if np.isnan(R2_scan4[i]) or np.isnan(R2_scan4[i+1]):
        continue
    if (R2_scan4[i] - R31_EXP) * (R2_scan4[i+1] - R31_EXP) < 0:
        try:
            g0_tf = brentq(lambda x: ratio_phi2(x) - R31_EXP,
                           g0_scan4[i], g0_scan4[i+1], xtol=1e-4)
            R_tf = ratio_phi2(g0_tf)
            tau_fp.append((g0_tf, R_tf))
            print(f"  ⭐ Tauon FP: g0** = {g0_tf:.6f},  R = {R_tf:.1f}  (vs 3477)")
        except Exception:
            pass

if not tau_fp:
    print("  Brak tauon fixed point w zakresie [1.10, 1.40].")
print()

# ============================================================
# CZĘŚĆ 5: Faza scatteringowa delta_scatt(g0) = atan2(B,C)
# ============================================================
print("-" * 70)
print("--- Część 5: Faza scatteringowa delta(g0) = atan2(B,C) ---")
print()

def scatt_phase(g0, alpha=2.0):
    r, g, _ = integrate_with_bounces(g0, alpha=alpha)
    A, B, C = fit_tail_window(r, g, WIN_LIST[2], WIN_LIST[2] + WIN_WIDTH)
    if A < 1e-5:
        return np.nan
    return np.arctan2(B, C)

N5 = 40
g0_scan5 = np.linspace(1.10, 2.60, N5)
phi_scan5 = [scatt_phase(g0) for g0 in g0_scan5]
phi_unwrapped = np.unwrap(phi_scan5)

print(f"  {'g0':>8}  {'phi_raw':>10}  {'phi_unwrapped':>14}")
for g0, pr, pu in zip(g0_scan5, phi_scan5, phi_unwrapped):
    if not np.isnan(pr):
        print(f"  {g0:8.4f}  {pr:10.4f}  {pu:14.4f}")

# Szukaj warunków Bohr-Sommerfeld: delta = n*pi
print()
print("  Warunki Bohr-Sommerfeld (delta = n*pi):")
for n_bs in range(6):
    target_bs = n_bs * np.pi
    for i in range(len(phi_unwrapped) - 1):
        if np.isnan(phi_unwrapped[i]) or np.isnan(phi_unwrapped[i+1]):
            continue
        if (phi_unwrapped[i] - target_bs) * (phi_unwrapped[i+1] - target_bs) < 0:
            try:
                g0_bs = brentq(lambda x: scatt_phase(x) - (target_bs % (2*np.pi) - 2*np.pi
                               if target_bs > np.pi else target_bs),
                               g0_scan5[i], g0_scan5[i+1], xtol=1e-3)
            except Exception:
                g0_bs = 0.5*(g0_scan5[i] + g0_scan5[i+1])
            print(f"    n={n_bs}: delta={target_bs:.4f} rad  →  g0_BS ≈ {g0_bs:.5f}")
            break
print()

# ============================================================
# CZĘŚĆ 6: Porównanie wszystkich "elektron g0" kandydatów
# ============================================================
print("-" * 70)
print("--- Część 6: Zestawienie kandydatów na g0_e ---")
print()

candidates = {
    "z0 (H1, B_tail=0)":      1.2301,
    "g0_e (ex57, ex60)":      1.240,
    "phi-FP (ex62 Cz.1)":     g0_star if g0_star else float('nan'),
    "sqrt(5)-1":               np.sqrt(5) - 1,
    "1/phi":                   1.0/PHI,
    "e^(-1/4) + 0.4":         np.exp(-0.25) + 0.4,
}

print(f"  {'Kandydat':30s}  {'g0':>10}  {'R=(A(phi*g0)/A(g0))^4':>24}  {'delta%':>8}")
print("  " + "-" * 80)
for name, g0_c in candidates.items():
    if np.isnan(g0_c):
        print(f"  {name:30s}  {'N/A':>10}")
        continue
    R_c = ratio_phi(g0_c)
    delta_c = (R_c - R21_EXP) / R21_EXP * 100 if not np.isnan(R_c) else np.nan
    print(f"  {name:30s}  {g0_c:10.6f}  {R_c:24.4f}  {delta_c:8.3f}%")

print()

# ============================================================
# PODSUMOWANIE I TESTY
# ============================================================
print("=" * 70)
print("TESTY")
print("=" * 70)

# T1: g0* (FP) ∈ [1.20, 1.30]
T1 = (g0_star is not None) and (1.20 <= g0_star <= 1.30)
print(f"  T1: g0* (FP) ∈ [1.20,1.30]: {'PASS' if T1 else 'FAIL'}  "
      f"(g0* = {g0_star:.6f if g0_star else 'N/A'})")

# T2: |R(g0*) - 207| / 207 < 0.005
if g0_star is not None:
    R_test = ratio_phi(g0_star)
    T2 = abs(R_test - R21_EXP) / R21_EXP < 0.005
    print(f"  T2: |R(g0*) - 207| / 207 < 0.5%: {'PASS' if T2 else 'FAIL'}  "
          f"(R = {R_test:.4f})")
else:
    T2 = False
    print(f"  T2: |R(g0*) - 207| / 207 < 0.5%: FAIL  (brak FP)")

# T3: alpha* ∈ [2.0, 3.0]
T3 = (alpha_star is not None) and (2.0 <= alpha_star <= 3.0)
print(f"  T3: alpha* ∈ [2.0,3.0]: {'PASS' if T3 else 'FAIL'}  "
      f"(alpha* = {alpha_star:.4f if alpha_star else 'N/A'})")

# T4: tauon FP |R - 3477| / 3477 < 0.3
if tau_fp:
    g0_tf, R_tf = tau_fp[0]
    T4 = abs(R_tf - R31_EXP) / R31_EXP < 0.3
    print(f"  T4: tauon FP |R-3477|/3477 < 30%: {'PASS' if T4 else 'FAIL'}  "
          f"(R = {R_tf:.1f})")
else:
    T4 = False
    print(f"  T4: tauon FP: FAIL  (brak FP)")

# T5: delta(z0) ≈ 0  i  delta(phi*z0) ≈ pi
delta_z0    = scatt_phase(1.2301)
delta_phi_z0 = scatt_phase(PHI * 1.2301)
T5_a = abs(delta_z0) < 0.5
T5_b = abs(abs(delta_phi_z0) - np.pi) < 0.5
T5   = T5_a and T5_b
print(f"  T5: delta(z0)≈0 i delta(phi*z0)≈pi: {'PASS' if T5 else 'FAIL'}  "
      f"(delta_z0={delta_z0:.4f}, delta_mu={delta_phi_z0:.4f} rad)")

n_pass = sum([T1, T2, T3, T4, T5])
print(f"\nWYNIK: {n_pass}/5 testów przeszło")

print()
print("=" * 70)
print("WNIOSEK FIZYCZNY (EX62)")
print("=" * 70)
print()
if g0_star:
    print(f"  1. phi-self-consistent fixed point:")
    print(f"     g0* = {g0_star:.6f}  (elektron)")
    print(f"     phi*g0* = {PHI*g0_star:.6f}  (mion)")
    print(f"     (A(phi*g0*)/A(g0*))^4 = {ratio_phi(g0_star):.4f} ≈ 207")
    print(f"     Δ(g0* - z0)  = {g0_star - 1.2301:.6f}")
    print(f"     Δ(g0* - 1.24) = {g0_star - 1.24:.6f}")
else:
    print("  1. Brak phi-self-consistent fixed point w [1.10, 1.45].")

if alpha_star:
    print(f"\n  2. Krytyczny parametr alpha*:")
    print(f"     alpha* = {alpha_star:.4f}  (czyste phi daje r21=207)")
    print(f"     TGP ma alpha=2.0 → odchylenie od alpha*: {abs(2.0-alpha_star):.4f}")

print()
print("=" * 70)
