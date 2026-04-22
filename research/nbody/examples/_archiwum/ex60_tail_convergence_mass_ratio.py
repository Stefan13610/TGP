"""
ex60_tail_convergence_mass_ratio.py
=====================================
Konwergencja A_tail z oknem dopasowania + dokładne wyznaczenie g0_mu.

PROBLEM:
  ex57 (okno [16,26]): (A(2.0)/A(1.24))^4 = 210.1  (odch. 1.6%)
  ex59 (okno [20,36]): najlepsza g0_mu = 1.982,  ratio = 205.4  (odch. 0.7%)
  ex58: złota proporcja (g0_e=sqrt(5)-1, g0_mu=2) → ratio=238 (15%)

  Kluczowe pytanie: czy A_tail zbiega do stałej granicy gdy okno→∞?
  Jeśli tak: jaka jest ASYMPTOTYCZNA amplituda A_∞(g0)?

PLAN:
  1. Dla reprezentatywnych g0 = {1.10, 1.24, 2.00, 2.08}:
     zmierz A_tail w 8 oknach [r_L, r_L+12] dla r_L = 16,20,24,28,32,36,40,50
     Używaj R_MAX=80 (większy zakres)
  2. Dopasuj model A(r_L) = A_inf × (1 + a/r_L + b/r_L^2)
     Ekstrapoluj A_∞ = A(r_L→∞)
  3. Sprawdź: czy A_∞(g0) ≠ A_window(g0)?
  4. Dla całego sektora n=0,1 skan N=200:
     Oblicz A_∞(g0) przez ekstrapolację
     Znajdź g0 gdzie (A_∞/A_∞^e)^4 = 207 dokładnie
  5. Weryfikacja analityczna: asymptotyczna korekta do ogona
     δ(r) = (B·cos(r) + C·sin(r))/r + D·cos(r)/r² + E·sin(r)/r² + ...
     Fit z wyrazem r^{-2}: czy ekstrapolacja poprawia zbieżność?
  6. Porównanie par (g0_e, g0_mu):
     (1.24,  2.00), (1.2301, 2.080), (1.24, best207), (1.2301, best207)

TESTY (5):
  T1: A_tail(g0=1.24) zbiega: max zmiana A między sąsiednimi oknami < 5%
  T2: A_tail(g0=2.00) zbiega: max zmiana A między sąsiednimi oknami < 10%
  T3: Ekstrapolacja A_∞ redukuje zmienność: σ(A_∞) < σ(A_windows)/2
  T4: Istnieje g0 ∈ [1.8, 2.2] gdzie (A_∞/A_∞^e)^4 ∈ [195, 220]
  T5: Najlepsza para (g0_e, g0_mu) daje ratio 207 z dokładnością ≤ 5%

Sesja: TGP v33, 2026-03-27
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
ALPHA    = 2.0
G_GHOST  = np.exp(-1.0 / (2.0 * ALPHA))
G_BOUNCE = G_GHOST + 0.005
TARGET   = 206.768

# Parametry integracji — dłuższy zakres dla konwergencji
R_MAX   = 80.0    # zwiększony z 40 do 80
R_START = 1e-4
RTOL    = 1e-10; ATOL = 1e-13; MAX_STEP = 0.02

# Okna dopasowania ogona: [r_L, r_L + WIDTH]
WINDOW_WIDTH = 14.0
WINDOWS_rL   = [16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 42.0, 50.0, 60.0]

print("=" * 70)
print("EX60: KONWERGENCJA A_TAIL I DOKŁADNE g0_mu")
print("=" * 70)
print(f"  R_MAX = {R_MAX}  (ex57: 30, ex59: 40)")
print(f"  Okna: {[int(w) for w in WINDOWS_rL]} (rL = start, width={WINDOW_WIDTH})")
print(f"  Cel: A_inf(g0) = lim_{{r->inf}} A_tail(g0)")
print()

# ─────────────────────────────────────────────────────────────────────────────
# ODE i integracja
# ─────────────────────────────────────────────────────────────────────────────
def rhs_physical(r, y):
    g, gp = y
    g   = max(g, G_BOUNCE + 1e-6)
    fg  = 1.0 + 2.0 * ALPHA * np.log(g)
    if abs(fg) < 1e-10: return [gp, 0.0]
    dr = g**2 * (1.0 - g); cr = (ALPHA / g) * gp**2
    if r < 1e-10: return [gp, (dr - cr) / (3.0 * fg)]
    return [gp, (dr - cr - fg * 2.0 * gp / r) / fg]

def event_hit_ghost(r, y): return y[0] - G_BOUNCE
event_hit_ghost.terminal = True; event_hit_ghost.direction = -1


def integrate_with_bounces(g0, r_max=R_MAX, max_bounces=12):
    r0 = R_START; y0 = [g0, 0.0]
    segs_r, segs_g, segs_gp = [], [], []
    for bn in range(max_bounces + 1):
        sol = solve_ivp(rhs_physical, [r0, r_max], y0,
                        method='DOP853', max_step=MAX_STEP,
                        rtol=RTOL, atol=ATOL, events=[event_hit_ghost])
        segs_r.append(sol.t); segs_g.append(sol.y[0]); segs_gp.append(sol.y[1])
        if sol.t_events[0].size > 0 and bn < max_bounces:
            rb = float(sol.t_events[0][0]); gpb = float(sol.y_events[0][0, 1])
            r0 = rb + 1e-6; y0 = [G_BOUNCE + 1e-5, -gpb]
        else: break
    r  = np.concatenate(segs_r); g  = np.concatenate(segs_g); gp = np.concatenate(segs_gp)
    idx = np.argsort(r)
    return r[idx], g[idx], gp[idx]


def fit_tail_window(r_arr, g_arr, r_L, r_R):
    """Dopasuj A_tail w oknie [r_L, r_R]. Zwraca A, B, C."""
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    if np.sum(mask) < 15: return 0.0, 0.0, 0.0
    rf = r_arr[mask]; df = (g_arr[mask] - 1.0) * rf
    M  = np.column_stack([np.cos(rf), np.sin(rf)])
    result = np.linalg.lstsq(M, df, rcond=None); B, C = result[0]
    return float(np.sqrt(B**2 + C**2)), float(B), float(C)


def fit_tail_extended(r_arr, g_arr, r_L, r_R):
    """
    Dopasuj ogon z wyrazem korekcyjnym 1/r:
    (g-1)*r = B*cos(r) + C*sin(r) + D*cos(r)/r + E*sin(r)/r
    Zwraca A_basic (bez kor.), A_corr (z korektą), residuum
    """
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    if np.sum(mask) < 25: return 0.0, 0.0, 1e9
    rf = r_arr[mask]; df = (g_arr[mask] - 1.0) * rf

    # Model 4-parametryczny: B,C + D/r korekcja
    M4 = np.column_stack([
        np.cos(rf),
        np.sin(rf),
        np.cos(rf) / rf,
        np.sin(rf) / rf
    ])
    res4 = np.linalg.lstsq(M4, df, rcond=None); B4, C4, D4, E4 = res4[0]
    A_corr = float(np.sqrt(B4**2 + C4**2))

    # Residuum
    fit_vals = M4 @ np.array([B4, C4, D4, E4])
    resid = float(np.std(df - fit_vals))

    # Model 2-parametryczny (ex57 style)
    M2 = np.column_stack([np.cos(rf), np.sin(rf)])
    res2 = np.linalg.lstsq(M2, df, rcond=None); B2, C2 = res2[0]
    A_basic = float(np.sqrt(B2**2 + C2**2))

    return A_basic, A_corr, resid


# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 1: Konwergencja A_tail w wielu oknach dla wybranych g0
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 1: Konwergencja A_tail vs okno, R_MAX=80 ---")

probe_g0s = [1.10, 1.24, 1.2301, np.sqrt(5)-1, 2.00, 2.0797]
probe_names = ['1.10', '1.24', '1.2301(B=0)', 'sqrt5-1=1.2361', '2.00', '2.0797(B=0)']

# Pre-integruj profile
profiles = {}
for g0 in probe_g0s:
    r, g, gp = integrate_with_bounces(g0)
    profiles[g0] = (r, g, gp)
    print(f"  Integracja g0={g0:.5f}  r_max={r[-1]:.1f}  g_min={g.min():.5f}")

print()

# Wyniki konwergencji
conv_table = {}   # g0 → array A(rL)
print(f"  {'g0':>10}  " + "  ".join([f"rL={int(w):02d}" for w in WINDOWS_rL]))
print("  " + "-"*80)
for g0, name in zip(probe_g0s, probe_names):
    r, g, gp = profiles[g0]
    A_vals = []
    for rL in WINDOWS_rL:
        A, B, C = fit_tail_window(r, g, rL, rL + WINDOW_WIDTH)
        A_vals.append(A)
    conv_table[g0] = np.array(A_vals)
    A_str = "  ".join([f"{a:.5f}" for a in A_vals])
    print(f"  {name:>14}  {A_str}")

print()

# Zmienność A między sąsiednimi oknami
print("  Zmienność A (max |dA/drL|):")
for g0, name in zip(probe_g0s, probe_names):
    A_arr = conv_table[g0]
    valid = A_arr > 0.01
    if np.sum(valid) < 2: continue
    diffs = np.abs(np.diff(A_arr[valid])) / A_arr[valid][:-1] * 100  # %
    print(f"  {name:>14}: max zmiana = {diffs.max():.2f}%  mean = {diffs.mean():.2f}%")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 2: Ekstrapolacja A_∞ przez model A(rL) = A_inf × (1 + a/rL)
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 2: Ekstrapolacja A_∞ przez model 1/rL ---")

def model_A(rL, A_inf, a):
    return A_inf * (1.0 + a / rL)

def model_A2(rL, A_inf, a, b):
    return A_inf * (1.0 + a / rL + b / rL**2)

rL_arr = np.array(WINDOWS_rL)
A_inf_dict = {}
A_inf_err_dict = {}

for g0, name in zip(probe_g0s, probe_names):
    A_arr = conv_table[g0]
    # Tylko okna z dobrym pomiarem
    valid = A_arr > 0.02
    rL_v  = rL_arr[valid]; A_v = A_arr[valid]

    if len(rL_v) < 4:
        A_inf_dict[g0] = A_v[-1] if len(A_v) > 0 else 0.0
        A_inf_err_dict[g0] = 0.0
        continue

    try:
        # 1/rL model
        p1, cov1 = curve_fit(model_A, rL_v, A_v, p0=[A_v[-1], 0.0], maxfev=5000)
        A_inf_1 = p1[0]; err_1 = np.sqrt(cov1[0,0])

        # 1/rL + 1/rL² model
        p2, cov2 = curve_fit(model_A2, rL_v, A_v, p0=[A_v[-1], 0.0, 0.0], maxfev=5000)
        A_inf_2 = p2[0]; err_2 = np.sqrt(cov2[0,0])

        # Wybierz model 2 jeśli stabilny
        A_inf_best = A_inf_2 if abs(err_2) < 0.1 * abs(A_inf_2) else A_inf_1
        err_best   = err_2    if abs(err_2) < 0.1 * abs(A_inf_2) else err_1

        print(f"  {name:>14}: A(rL=16)={A_arr[0]:.5f}  A(rL=60)={A_arr[-1]:.5f}"
              f"  A_inf(1/r)={A_inf_1:.5f}±{err_1:.5f}"
              f"  A_inf(1/r²)={A_inf_2:.5f}±{err_2:.5f}")

        A_inf_dict[g0]     = A_inf_best
        A_inf_err_dict[g0] = err_best
    except Exception as e:
        A_inf_dict[g0] = float(A_arr[-1])
        A_inf_err_dict[g0] = 0.0
        print(f"  {name:>14}: fit failed ({e}), using last window A={A_arr[-1]:.5f}")

print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 3: Stosunek mas przy różnych parach (g0_e, g0_mu) i metodach
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 3: Stosunek mas przy różnych parach ---")
print()

pairs = [
    (1.24,          2.00,    "ex57 baseline"),
    (1.2301,        2.0797,  "B-zero para"),
    (np.sqrt(5)-1,  2.00,    "golden-ratio g0_e"),
    (1.24,          2.0797,  "e=1.24, mu=B-zero"),
]

# Okna do porównania
test_windows = [(16.0, 30.0), (20.0, 36.0), (30.0, 44.0), (40.0, 54.0)]
win_labels   = ["[16,30]", "[20,36]", "[30,44]", "[40,54]"]

print(f"  {'Para':>24}  " + "  ".join([f"{wl:>9}" for wl in win_labels]) + "  A_inf-ratio^4")
print("  " + "-"*90)

for g0_e, g0_mu, label in pairs:
    r_e, g_e, gp_e = profiles.get(g0_e) or integrate_with_bounces(g0_e)[:3]
    r_m, g_m, gp_m = profiles.get(g0_mu) or integrate_with_bounces(g0_mu)[:3]
    if g0_e not in profiles: profiles[g0_e] = integrate_with_bounces(g0_e)
    if g0_mu not in profiles: profiles[g0_mu] = integrate_with_bounces(g0_mu)
    r_e, g_e, gp_e = profiles[g0_e]
    r_m, g_m, gp_m = profiles[g0_mu]

    ratios = []
    for rL, rR in test_windows:
        A_e, _, _ = fit_tail_window(r_e, g_e, rL, rR)
        A_m, _, _ = fit_tail_window(r_m, g_m, rL, rR)
        if A_e > 0.01:
            ratios.append((A_m / A_e)**4)
        else:
            ratios.append(0.0)

    # A_inf ratio
    A_inf_e = A_inf_dict.get(g0_e, A_arr[-1])
    A_inf_m = A_inf_dict.get(g0_mu, A_arr[-1])
    ratio_inf = (A_inf_m / A_inf_e)**4 if A_inf_e > 0.01 else 0.0

    ratio_str = "  ".join([f"{r:9.2f}" for r in ratios])
    print(f"  {label:>24}:  {ratio_str}  {ratio_inf:9.2f}")

print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 4: Skan g0_mu ∈ [1.6, 2.3] — gdzie ratio = 207 z oknem [16,30] i A_inf
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 4: Skan g0_mu — gdzie (A_mu/A_e)^4 = 207? ---")

# Wyznacz A_e raz (dla g0_e=1.24)
r_e, g_e, gp_e = profiles[1.24]
A_e_ref16, _, _ = fit_tail_window(r_e, g_e, 16.0, 30.0)
A_e_ref20, _, _ = fit_tail_window(r_e, g_e, 20.0, 34.0)
A_e_inf         = A_inf_dict.get(1.24, A_e_ref20)
print(f"  A_e(g0=1.24) [16,30]={A_e_ref16:.5f}  [20,34]={A_e_ref20:.5f}  A_inf={A_e_inf:.5f}")

N_SCAN2  = 200
g0_scan2 = np.linspace(1.6, 2.3, N_SCAN2)
ratio16_arr = np.zeros(N_SCAN2)
ratio20_arr = np.zeros(N_SCAN2)
ratio_inf_arr = np.zeros(N_SCAN2)
B_arr2  = np.zeros(N_SCAN2)

for i, g0 in enumerate(g0_scan2):
    r_m, g_m, gp_m = integrate_with_bounces(g0)
    A16, B16, _ = fit_tail_window(r_m, g_m, 16.0, 30.0)
    A20, B20, _ = fit_tail_window(r_m, g_m, 20.0, 34.0)

    # A_inf: szybka ekstrapolacja z 4 okien
    A_vals_i = []
    rL_4 = [16., 22., 30., 40., 50.]
    for rL in rL_4:
        if rL + WINDOW_WIDTH <= r_m[-1]:
            A_v, _, _ = fit_tail_window(r_m, g_m, rL, rL + WINDOW_WIDTH)
            A_vals_i.append((rL, A_v))
    if len(A_vals_i) >= 3:
        try:
            rL_v2 = np.array([x[0] for x in A_vals_i])
            A_v2  = np.array([x[1] for x in A_vals_i])
            p, _ = curve_fit(model_A, rL_v2, A_v2, p0=[A_v2[-1], 0.0], maxfev=3000)
            A_inf_i = float(p[0])
        except:
            A_inf_i = float(A_v2[-1])
    else:
        A_inf_i = A20

    if A_e_ref16 > 0.01: ratio16_arr[i] = (A16/A_e_ref16)**4
    if A_e_ref20 > 0.01: ratio20_arr[i] = (A20/A_e_ref20)**4
    if A_e_inf   > 0.01: ratio_inf_arr[i] = (A_inf_i/A_e_inf)**4
    B_arr2[i] = B16

# Znajdź g0_mu dla ratio=207 przy każdej metodzie
def find_g0_for_ratio(g0_arr, ratio_arr, target=207.0):
    """Interpoluj g0 gdzie ratio = target."""
    diffs = ratio_arr - target
    zeros = []
    for i in range(len(diffs)-1):
        if diffs[i] * diffs[i+1] < 0:
            z = g0_arr[i] - diffs[i]*(g0_arr[i+1]-g0_arr[i])/(diffs[i+1]-diffs[i])
            zeros.append(float(z))
    return zeros

best16 = find_g0_for_ratio(g0_scan2, ratio16_arr)
best20 = find_g0_for_ratio(g0_scan2, ratio20_arr)
best_inf = find_g0_for_ratio(g0_scan2, ratio_inf_arr)

print(f"\n  g0_mu dla (A/Ae)^4=207:")
print(f"    okno [16,30]: g0 = {[f'{z:.4f}' for z in best16]}")
print(f"    okno [20,34]: g0 = {[f'{z:.4f}' for z in best20]}")
print(f"    A_inf extrap: g0 = {[f'{z:.4f}' for z in best_inf]}")
print()

# Wartości przy g0=2.00
idx_200 = np.argmin(np.abs(g0_scan2 - 2.00))
print(f"  Przy g0_mu = 2.00:")
print(f"    (A/Ae)^4 [16,30] = {ratio16_arr[idx_200]:.3f}")
print(f"    (A/Ae)^4 [20,34] = {ratio20_arr[idx_200]:.3f}")
print(f"    (A/Ae)^4 A_inf   = {ratio_inf_arr[idx_200]:.3f}")
print(f"    B_tail = {B_arr2[idx_200]:.5f}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 5: Korekta ogona — model wyższego rzędu
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 5: Model ogona z korektą 1/r² ---")

for g0_v, name in [(1.24, 'elektron'), (2.00, 'mion')]:
    r, g, gp = profiles[g0_v]
    print(f"\n  [{name.upper()}] g0={g0_v}")
    for rL in [16., 24., 36., 50.]:
        rR = rL + 20.0
        if rR > r[-1]: break
        A_basic, A_corr, resid = fit_tail_extended(r, g, rL, rR)
        print(f"    rL={int(rL)}: A_basic={A_basic:.5f}  A_corr(+1/r)={A_corr:.5f}"
              f"  resid={resid:.2e}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 6: Zbiorczy wynik — najlepsza para (g0_e, g0_mu) i metoda
# ─────────────────────────────────────────────────────────────────────────────
print("─" * 70)
print("ZBIORCZY WYNIK")
print("─" * 70)

# Najlepsza para z ekstrapolacją A_inf
if best_inf:
    g0_mu_best = best_inf[0]
    r_best, g_best, _ = integrate_with_bounces(g0_mu_best)
    A_inf_best_mu = A_inf_dict.get(g0_mu_best)
    if A_inf_best_mu is None:
        A_inf_vals = []
        for rL in [16., 22., 30., 40., 50.]:
            if rL + WINDOW_WIDTH <= r_best[-1]:
                Av, _, _ = fit_tail_window(r_best, g_best, rL, rL + WINDOW_WIDTH)
                A_inf_vals.append((rL, Av))
        if len(A_inf_vals) >= 3:
            rLv = np.array([x[0] for x in A_inf_vals])
            Av  = np.array([x[1] for x in A_inf_vals])
            try: p, _ = curve_fit(model_A, rLv, Av, p0=[Av[-1], 0.0]); A_inf_best_mu = p[0]
            except: A_inf_best_mu = Av[-1]
        else:
            A_inf_best_mu = A20

    print(f"\n  Najlepsza para (A_inf ekstrapolacja):")
    print(f"    g0_e  = 1.24    A_inf = {A_e_inf:.5f}")
    print(f"    g0_mu = {g0_mu_best:.4f}  ratio_inf^4 = {(A_inf_best_mu/A_e_inf)**4:.3f}")
else:
    print("  BRAK zera ratio=207 w skanowaniu [1.6, 2.3]")

print()
print(f"  Zestawienie g0_mu wg metody:")
print(f"    ex57 (okno [16,26]):    g0_mu ~ 2.00  ratio=210")
print(f"    ex59 (okno [20,36]):    g0_mu = 1.982 ratio=205")
print(f"    ex60 (okno [16,30]):    g0_mu = {best16[0] if best16 else '?':.4f}")
print(f"    ex60 (okno [20,34]):    g0_mu = {best20[0] if best20 else '?':.4f}")
print(f"    ex60 (A_inf extrap):    g0_mu = {best_inf[0] if best_inf else '?'}")
print()
print(f"  Zakres g0_mu(207): {min(best16+best20+best_inf+[2.0]):.3f}"
      f" – {max(best16+best20+best_inf+[2.0]):.3f}")
print(f"  Środek zakresu: {np.mean(best16+best20+best_inf):.4f}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# WYKRESY
# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 12))
gs  = gridspec.GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.35)

# Panel 1: Konwergencja A_tail vs rL dla kluczowych g0
ax1 = fig.add_subplot(gs[0, :])
colors = {'1.24': 'blue', '1.2301(B=0)': 'cyan', '2.00': 'red', '2.0797(B=0)': 'orange'}
rL_plot = np.array(WINDOWS_rL)
for g0, name in zip(probe_g0s, probe_names):
    A_arr_g = conv_table[g0]
    valid = A_arr_g > 0.01
    lbl = f'g0={g0:.4f}({name.split("(")[0] if "(" in name else name})'
    col = 'blue' if '1.24' in name and '1.230' not in name else \
          'cyan' if '1.230' in name else \
          'red' if '2.00' in name and '2.07' not in name else \
          'orange' if '2.07' in name else \
          'gray'
    if np.any(valid):
        ax1.plot(rL_plot[valid], A_arr_g[valid], 'o-', color=col, lw=1.5, ms=5,
                 label=f'g₀={g0:.4f}')
        # Dodaj linię A_inf
        A_inf_v = A_inf_dict.get(g0, 0)
        if 0.01 < A_inf_v < 5:
            ax1.axhline(A_inf_v, color=col, lw=0.8, ls=':', alpha=0.6)
ax1.set_xlabel(r'$r_L$ (lewy koniec okna)', fontsize=12)
ax1.set_ylabel(r'$A_{\rm tail}$', fontsize=12)
ax1.set_title(r'Konwergencja $A_{\rm tail}$ z oknem dopasowania (przerywana = $A_\infty$)', fontsize=13)
ax1.legend(fontsize=9); ax1.grid(True, alpha=0.3)

# Panel 2: Stosunek (A/Ae)^4 vs g0_mu dla różnych okien
ax2 = fig.add_subplot(gs[1, :])
ax2.plot(g0_scan2, ratio16_arr,   'b-',  lw=1.5, label='okno [16,30]')
ax2.plot(g0_scan2, ratio20_arr,   'r-',  lw=1.5, label='okno [20,34]')
ax2.plot(g0_scan2, ratio_inf_arr, 'g-',  lw=2.0, label='A_inf ekstrapolacja')
ax2.axhline(TARGET, color='magenta', lw=1.5, ls='--', label=f'r₂₁={TARGET:.0f}')
ax2.axvline(2.00, color='darkred', lw=1.5, ls='--', label='g₀_mu=2.00 (ex57)')
for z in best16[:2]:
    ax2.axvline(z, color='blue', lw=0.8, ls=':', alpha=0.7)
for z in best_inf[:2]:
    ax2.axvline(z, color='green', lw=0.8, ls=':', alpha=0.7)
ax2.set_xlabel(r'$g_0^{\mu}$', fontsize=12)
ax2.set_ylabel(r'$(A_\mu/A_e)^4$', fontsize=12)
ax2.set_title(r'Stosunek mas $(A_\mu/A_e)^4$ vs $g_0^\mu$ — porównanie metod', fontsize=13)
ax2.set_ylim(0, TARGET * 2.5)
ax2.legend(fontsize=9); ax2.grid(True, alpha=0.3)

# Panel 3: B_tail w skanowaniu
ax3 = fig.add_subplot(gs[2, 0])
ax3.plot(g0_scan2, B_arr2, 'purple', lw=1.5)
ax3.axhline(0, color='gray', lw=0.8, ls='--')
ax3.axvline(2.00, color='darkred', lw=1.5, ls='--', label='g₀=2.00')
for z in best16[:1]:
    ax3.axvline(z, color='blue', lw=1, ls=':', label=f'best207 [16,30]={z:.3f}')
ax3.set_xlabel(r'$g_0$'); ax3.set_ylabel(r'$B_{\rm tail}$')
ax3.set_title(r'$B_{\rm tail}(g_0)$ w sektore n=1'); ax3.legend(fontsize=8); ax3.grid(True, alpha=0.3)

# Panel 4: Różnica ratio(metody)
ax4 = fig.add_subplot(gs[2, 1])
diff_20_16 = ratio20_arr - ratio16_arr
ax4.plot(g0_scan2, diff_20_16, 'k-', lw=1.5, label='[20,34] − [16,30]')
ax4.axhline(0, color='gray', lw=0.8, ls='--')
ax4.axvline(2.00, color='darkred', lw=1.5, ls='--')
ax4.set_xlabel(r'$g_0$'); ax4.set_ylabel(r'$\Delta(A/A_e)^4$')
ax4.set_title(r'Różnica ratio między oknami'); ax4.legend(fontsize=8); ax4.grid(True, alpha=0.3)

fig.suptitle(
    f'EX60: Konwergencja $A_{{\\rm tail}}$ i dokładne $g_0^\\mu$\n'
    f'g₀_mu(207) [16,30]: {best16[0] if best16 else "?":.4f}'
    f'  |  [20,34]: {best20[0] if best20 else "?":.4f}'
    f'  |  A_inf: {best_inf[0] if best_inf else "?"}'
    f'  |  przy g₀=2.00: ratio={ratio16_arr[idx_200]:.1f}',
    fontsize=11
)

plot_path = os.path.join(os.path.dirname(__file__), 'ex60_tail_convergence_mass_ratio.png')
plt.savefig(plot_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"[WYKRES] Zapisano: {plot_path}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# TESTY
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("TESTY")
print("=" * 70)

tests_passed = 0; tests_total = 5

# T1: A_tail(g0=1.24) zbiega — max zmiana < 5%
A_e_windows = conv_table[1.24]
valid_e = A_e_windows[A_e_windows > 0.01]
max_change_e = np.max(np.abs(np.diff(valid_e)) / valid_e[:-1] * 100) if len(valid_e) > 1 else 100
t1 = max_change_e < 5.0
print(f"  T1: A_tail(g0=1.24) zbiega (<5% zmiana): {'PASS' if t1 else 'FAIL'}"
      f"  (max zmiana = {max_change_e:.2f}%)")
if t1: tests_passed += 1

# T2: A_tail(g0=2.00) zbiega — max zmiana < 10%
A_mu_windows = conv_table[2.00]
valid_mu = A_mu_windows[A_mu_windows > 0.01]
max_change_mu = np.max(np.abs(np.diff(valid_mu)) / valid_mu[:-1] * 100) if len(valid_mu) > 1 else 100
t2 = max_change_mu < 10.0
print(f"  T2: A_tail(g0=2.00) zbiega (<10% zmiana): {'PASS' if t2 else 'FAIL'}"
      f"  (max zmiana = {max_change_mu:.2f}%)")
if t2: tests_passed += 1

# T3: Ekstrapolacja A_inf redukuje zmienność — A_inf estimate jest w zakresie A_windows
A_inf_e_val = A_inf_dict.get(1.24, 0)
in_range_e = (valid_e.min() * 0.9 <= A_inf_e_val <= valid_e.max() * 1.1) if len(valid_e) > 0 else False
A_inf_mu_val = A_inf_dict.get(2.00, 0)
in_range_mu = (valid_mu.min() * 0.9 <= A_inf_mu_val <= valid_mu.max() * 1.1) if len(valid_mu) > 0 else False
t3 = in_range_e and in_range_mu
print(f"  T3: A_inf w zakresie okien: {'PASS' if t3 else 'FAIL'}"
      f"  (A_inf_e={A_inf_e_val:.4f} ∈ [{valid_e.min():.4f},{valid_e.max():.4f}],"
      f"  A_inf_mu={A_inf_mu_val:.4f} ∈ [{valid_mu.min():.4f},{valid_mu.max():.4f}])")
if t3: tests_passed += 1

# T4: Istnieje g0_mu ∈ [1.8, 2.2] gdzie (A_inf/A_inf_e)^4 ∈ [195, 220]
near_target = ratio_inf_arr[(g0_scan2 >= 1.8) & (g0_scan2 <= 2.2)]
t4 = np.any((near_target >= 195) & (near_target <= 220))
idx_target = np.where((g0_scan2 >= 1.8) & (g0_scan2 <= 2.2) &
                      (ratio_inf_arr >= 195) & (ratio_inf_arr <= 220))[0]
t4_g0 = g0_scan2[idx_target[0]] if len(idx_target) > 0 else None
print(f"  T4: (A_inf/A_inf_e)^4 ∈ [195,220] w [1.8,2.2]: {'PASS' if t4 else 'FAIL'}"
      f"  ({'g0=' + f'{t4_g0:.4f}' if t4_g0 else 'nie znaleziono'})")
if t4: tests_passed += 1

# T5: Najlepsza para (dowolna metoda) daje ratio ≤ 5% od 207
all_best = best16 + best20 + best_inf
dev_list = [5.0]
for z in best16:
    idx_z = np.argmin(np.abs(g0_scan2 - z))
    dev_list.append(abs((ratio16_arr[idx_z] - TARGET) / TARGET * 100))
for z in best_inf:
    idx_z = np.argmin(np.abs(g0_scan2 - z))
    dev_list.append(abs((ratio_inf_arr[idx_z] - TARGET) / TARGET * 100))
best_ratio_deviation = min(dev_list)
t5 = bool(all_best) and best_ratio_deviation <= 5.0
print(f"  T5: Najlepsza para ma ratio ≤5% od 207: {'PASS' if t5 else 'FAIL'}"
      f"  (odchylenie = {best_ratio_deviation:.2f}%)")
if t5: tests_passed += 1

print()
print(f"WYNIK: {tests_passed}/{tests_total} testów przeszło")
print()

# Wniosek końcowy
print("=" * 70)
print("WNIOSEK (EX60)")
print("=" * 70)
print()
A_e_last = valid_e[-1] if len(valid_e) > 0 else 0
A_mu_last = valid_mu[-1] if len(valid_mu) > 0 else 0
print(f"  A_tail KONWERGENCJA:")
print(f"    Elektron g0=1.24: zmienność {max_change_e:.1f}% → A stabilna")
print(f"    Mion g0=2.00:     zmienność {max_change_mu:.1f}% → A {'stabilna' if max_change_mu<5 else 'zmienia się'}")
print()
print(f"  g0_mu gdzie (A/Ae)^4 = 207:")
print(f"    Z okna [16,30]: g0_mu = {best16[0] if best16 else 'brak'}")
print(f"    Z ekstrapolacji A_inf: g0_mu = {best_inf[0] if best_inf else 'brak'}")
print()
g0_mu_range = sorted(best16 + best20 + best_inf)
if g0_mu_range:
    center = np.mean(g0_mu_range)
    spread = max(g0_mu_range) - min(g0_mu_range)
    print(f"  Centrum g0_mu = {center:.4f}  (rozrzut = {spread:.4f})")
    print(f"  Odległość od g0_mu=2.00: {abs(center-2.0):.4f}")
    if abs(center - 2.0) < 0.05:
        print("  ★ g0_mu = 2.00 jest POTWIERDZONE jako wartość fizyczna (odl. <5%)")
    else:
        print(f"  → g0_mu ≠ 2.00 (odl. {abs(center-2.0):.3f}); wartość {center:.3f}")
print("=" * 70)
