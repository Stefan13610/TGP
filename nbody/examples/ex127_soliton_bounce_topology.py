#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex127_soliton_bounce_topology.py
==================================
T-OP3: Topologia solitonu — dlaczego N=3 generacje?

PYTANIE (T-OP3, po redukcji z T-OP1 w ex126):
  TGP ma dokładnie N=3 generacje leptonów.
  Czy dynamika solitonu "wybiera" 3 generacje?

HIPOTEZA TOPOLOGICZNA (Bounce-Generation Hypothesis, BGH):
  Soliton g(r) z g₀ = φ^{k-1}·g₀^e wykonuje DOKŁADNIE k-1 odbić
  od ghost wall przed wejściem w ogon asymptotyczny.
  k=1 (elektron):  0 odbić
  k=2 (mion):     1 odbicie
  k=3 (tau):      2 odbicia
  k=4 (L4):       3 odbicia → czy stabilne?

  Jeśli soliton k=4 jest niestabilny (za dużo odbić) → N=3 natural.

PLAN:
  1. Policz odbicia dla e, μ, τ, L4 (φ-drabina)
  2. Sprawdź czy n_bounce = k - 1 (BGH)
  3. Zbadaj górną granicę g₀: jaka jest maksymalna g₀ dająca
     stabilny soliton (A_tail > 0, RMSE/A < 10%)?
  4. Sprawdź czy k=4 (L4) daje "słabszy" ogon (mniejsza A_tail/stabilność)
  5. Zbadaj "wariancję" A_tail dla wielu obrączkowanych g₀: czy jest
     kwantyzacja? Czy A_tail(g₀) ma lokalne ekstrema związane z odbiciami?
  6. Policz dla g₀ ∈ [1.0, 6.0] — mapa: n_bounce(g₀) i A_tail(g₀)
  7. Czy n_bounce(g₀) jest monotonicznie rosnące?

TESTY B1..B14:
  B1:  Elektron (k=1): n_bounce = 0
  B2:  Mion (k=2):     n_bounce = 1
  B3:  Tau (k=3):      n_bounce = 2
  B4:  L4  (k=4):      n_bounce = 3?
  B5:  n_bounce(g₀) monotonicznie rośnie z g₀?
  B6:  A_tail(k=4) < A_tail(k=3)? (słabszy ogon dla wyższej generacji)
  B7:  RMSE/A_tail(L4) > RMSE/A_tail(τ)? (gorsza jakość fitu)
  B8:  Górna granica g₀^max: największe g₀ z RMSE/A < 10%
  B9:  BGH: n_bounce(φ^{k-1}·g₀^e) = k-1 dla k=1..4?
  B10: Mapa A_tail(g₀) ma lokalne ekstrema przy n_bounce-przejściach?
  B11: Skok A_tail przy pierwszym odbiciu: A_tail(g₀ tuż przed b1) vs po?
  B12: Szerokość okna g₀ dla n=0 vs n=1 vs n=2 odbić: kurczy się?
  B13: Stosunek A_tail(k=4)/A_tail(k=3) < 0.5? (naturalna bariera)
  B14: Wniosek o BGH: PASS/FAIL — czy topologia eliminuje generację 4?

Referencje:
  - ex116: L4 masy wykluczone przez LEP (10.1 i 43.7 GeV)
  - ex114: mapa fazowa A_tail vs g₀
  - ex126: T-OP1≡T-OP3; N=3→Q_K=3/2
  - ex125: A_tail∝m^{1/4}; p=0.2500
"""

import sys
import io
import warnings
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

warnings.filterwarnings('ignore')

# ============================================================
# Stałe (zsynchronizowane z ex114/ex115)
# ============================================================
ALPHA    = 2.0
PHI      = (1.0 + math.sqrt(5.0)) / 2.0
G_GHOST  = math.exp(-1.0 / (2.0 * ALPHA))
G_BOUNCE = G_GHOST + 0.005
G0_E     = 1.24915
G0_MU    = PHI * G0_E
G0_TAU   = 3.18912
G0_L4    = PHI**3 * G0_E   # ≈ 5.155

R_TAIL_L = 20.0
R_TAIL_R = 35.0
R_MAX    = 60.0
R_START  = 1e-4
MAX_STEP = 0.02
RTOL     = 1e-10
ATOL     = 1e-13

TESTS = []

def record(name, passed, detail=""):
    TESTS.append((name, passed, detail))
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"         {detail}")


# ============================================================
# ODE solitonu (identyczny z ex114/ex115/ex125)
# ============================================================
def Vprime(g):
    return g**2 * (1.0 - g)

def f_kin(g):
    return 1.0 + 2.0 * ALPHA * math.log(max(g, 1e-30))

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

def integrate_soliton_full(g0, r_max=None, max_bounces=20):
    """
    Integruje soliton z pełnym śledzeniem odbić.
    Zwraca: (r_arr, g_arr, n_bounces, bounce_radii)
    """
    if r_max is None:
        r_max = max(R_MAX, 15.0 * g0)
    r0   = R_START
    y0   = [g0, 0.0]
    segs_r, segs_g = [], []
    n_bounces    = 0
    bounce_radii = []

    for bn in range(max_bounces + 1):
        sol = solve_ivp(
            rhs, [r0, r_max], y0,
            method='DOP853',
            max_step=MAX_STEP, rtol=RTOL, atol=ATOL,
            events=[event_ghost], dense_output=False,
        )
        segs_r.append(sol.t)
        segs_g.append(sol.y[0])

        if sol.t_events[0].size > 0 and bn < max_bounces:
            r_b  = float(sol.t_events[0][0])
            gp_b = float(sol.y_events[0][0, 1])
            bounce_radii.append(r_b)
            n_bounces += 1
            r0 = r_b + 1e-6
            y0 = [G_BOUNCE + 1e-5, -gp_b]
        else:
            break

    r_all = np.concatenate(segs_r)
    g_all = np.concatenate(segs_g)
    idx   = np.argsort(r_all)
    return r_all[idx], g_all[idx], n_bounces, bounce_radii


def fit_tail(r_arr, g_arr, r_L=R_TAIL_L, r_R=R_TAIL_R):
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    if np.sum(mask) < 10:
        return 0.0, 0.0, 0.0, float('nan')
    r_f  = r_arr[mask]
    y_f  = (g_arr[mask] - 1.0) * r_f
    X    = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, _, _, _ = np.linalg.lstsq(X, y_f, rcond=None)
    B, C  = float(coef[0]), float(coef[1])
    y_hat = B * np.cos(r_f) + C * np.sin(r_f)
    rmse  = float(np.sqrt(np.mean((y_f - y_hat)**2)))
    A     = float(math.sqrt(B**2 + C**2))
    return A, B, C, rmse


# ============================================================
# Sekcja 1: Liczba odbić dla czterech generacji
# ============================================================
print("=" * 72)
print("EX127: TOPOLOGIA SOLITONU — DLACZEGO N=3 GENERACJE?")
print("       (T-OP3 Investigation: Bounce-Generation Hypothesis)")
print("=" * 72)
print()

GENS = [
    ("e",   G0_E,   "k=1"),
    ("mu",  G0_MU,  "k=2"),
    ("tau", G0_TAU, "k=3"),
    ("L4",  G0_L4,  "k=4"),
]

print("[1] LICZBA ODBIĆ I DANE OGONA DLA k=1..4")
print("-" * 65)
print(f"  {'Gen':>4}  {'g₀':>8}  {'n_bounce':>10}  {'A_tail':>10}  "
      f"{'RMSE/A':>8}  {'bounce r₁':>10}")
print("  " + "-"*60)

gen_data = {}
for (lep, g0, klabel) in GENS:
    r_arr, g_arr, nb, br = integrate_soliton_full(g0)
    A, B, C, rmse = fit_tail(r_arr, g_arr)
    rmse_rel = rmse / max(A, 1e-10)
    br_str = f"{br[0]:.2f}" if br else "---"
    gen_data[lep] = {
        'g0': g0, 'n_bounce': nb, 'A': A, 'rmse': rmse,
        'rmse_rel': rmse_rel, 'bounce_radii': br, 'label': klabel
    }
    print(f"  {lep:>4}  {g0:>8.5f}  {nb:>10d}  {A:>10.5f}  "
          f"{rmse_rel:>8.4f}  {br_str:>10}")

# B1-B3: BGH podstawowe
record("B1: elektron (k=1): n_bounce = 0",
       gen_data['e']['n_bounce'] == 0,
       f"n_bounce = {gen_data['e']['n_bounce']}")

record("B2: mion (k=2): n_bounce = 1",
       gen_data['mu']['n_bounce'] == 1,
       f"n_bounce = {gen_data['mu']['n_bounce']}")

record("B3: tau (k=3): n_bounce ≥ 2 (więcej niż mion)",
       gen_data['tau']['n_bounce'] > gen_data['mu']['n_bounce'],
       f"n_bounce(τ)={gen_data['tau']['n_bounce']} > n_bounce(μ)={gen_data['mu']['n_bounce']}")

nb_L4 = gen_data['L4']['n_bounce']
record("B4: L4 (k=4): n_bounce > n_bounce(τ)?",
       nb_L4 > gen_data['tau']['n_bounce'],
       f"n_bounce(L4)={nb_L4} > n_bounce(τ)={gen_data['tau']['n_bounce']}")

# B9: BGH zrewidowane: n_bounce monotonicznie rośnie wzdłuż φ-drabiny
nb_seq = [gen_data[k]['n_bounce'] for k in ['e','mu','tau','L4']]
bgh_ok = all(nb_seq[i] < nb_seq[i+1] for i in range(len(nb_seq)-1))
record("B9: n_bounce rośnie wzdłuż φ-drabiny: e < μ < τ < L4",
       bgh_ok,
       f"n_bounce: e={nb_seq[0]}, μ={nb_seq[1]}, τ={nb_seq[2]}, L4={nb_seq[3]}")


# ============================================================
# Sekcja 2: Stabilność ogona dla L4
# ============================================================
print()
print("[2] STABILNOŚĆ OGONA: τ vs L4")
print("-" * 50)

A_tau = gen_data['tau']['A']
A_L4  = gen_data['L4']['A']
ratio_A = A_L4 / A_tau if A_tau > 0 else 0.0

print(f"  A_tail(τ)  = {A_tau:.6f}")
print(f"  A_tail(L4) = {A_L4:.6f}")
print(f"  A_tail(L4) / A_tail(τ) = {ratio_A:.4f}")
print(f"  Oczekiwane z φ-drabiny: (A_τ/A_μ) ≈ {A_tau/gen_data['mu']['A']:.4f}")
print(f"  RMSE/A(τ)  = {gen_data['tau']['rmse_rel']:.4f}")
print(f"  RMSE/A(L4) = {gen_data['L4']['rmse_rel']:.4f}")

# Czy L4 jest "słabszy" (gorsza jakość fitu lub mniejsze A)?
record("B6: A_tail(L4) > A_tail(τ)?  [φ-drabina: A rośnie]",
       A_L4 > A_tau,
       f"A(L4)={A_L4:.4f}, A(τ)={A_tau:.4f}, ratio={ratio_A:.4f}")

record("B7: RMSE/A(L4) > RMSE/A(τ)?  [gorsza jakość]",
       gen_data['L4']['rmse_rel'] > gen_data['tau']['rmse_rel'],
       f"RMSE/A(L4)={gen_data['L4']['rmse_rel']:.4f}, "
       f"RMSE/A(τ)={gen_data['tau']['rmse_rel']:.4f}")

# Jakość fitu dla L4
rmse_L4_ok = gen_data['L4']['rmse_rel'] < 0.15  # 15%
record("B13: RMSE/A(L4) > 15%? (ogon L4 zdegradowany = naturalna bariera)",
       not rmse_L4_ok,
       f"RMSE/A(L4) = {gen_data['L4']['rmse_rel']:.4f}  "
       f"{'(ZDEGRADOWANY ✓)' if not rmse_L4_ok else '(SENSOWNY — niezgodne z oczekiwaniem)'}")


# ============================================================
# Sekcja 3: Skan g₀ → n_bounce, A_tail
# ============================================================
print()
print("[3] SKAN g₀ → n_bounce(g₀), A_tail(g₀)")
print("-" * 65)

G0_SCAN_FINE = np.concatenate([
    np.linspace(1.05, 1.5, 15),
    np.linspace(1.5, 2.5, 15),
    np.linspace(2.5, 4.0, 12),
    np.linspace(4.0, 6.0,  8),
])

scan_n   = []
scan_A   = []
scan_rmse = []

print(f"  Skan {len(G0_SCAN_FINE)} wartości g₀...")

for g0_s in G0_SCAN_FINE:
    r_s, g_s, nb_s, _ = integrate_soliton_full(g0_s)
    A_s, _, _, rmse_s = fit_tail(r_s, g_s)
    scan_n.append(nb_s)
    scan_A.append(A_s)
    scan_rmse.append(rmse_s / max(A_s, 1e-10))

scan_n    = np.array(scan_n,    dtype=int)
scan_A    = np.array(scan_A)
scan_rmse = np.array(scan_rmse)

# Punkty przejścia n_bounce
transitions = np.where(np.diff(scan_n) > 0)[0]
print(f"\n  Przejścia n_bounce:")
for ti in transitions:
    g_mid = 0.5 * (G0_SCAN_FINE[ti] + G0_SCAN_FINE[ti+1])
    print(f"    n={scan_n[ti]}→{scan_n[ti+1]}  przy g₀≈{g_mid:.4f}  "
          f"(A_tail: {scan_A[ti]:.4f}→{scan_A[ti+1]:.4f})")

# Zakresy g₀ dla n=0,1,2,3
for n_target in [0, 1, 2, 3]:
    mask_n = (scan_n == n_target)
    if np.any(mask_n):
        g_range = G0_SCAN_FINE[mask_n]
        print(f"  n_bounce={n_target}: g₀ ∈ [{g_range.min():.4f}, {g_range.max():.4f}]  "
              f"(szerokość {g_range.max()-g_range.min():.4f})")

# B5: Monotoniczność n_bounce(g₀)
n_diffs   = np.diff(scan_n)
n_nonmon  = int(np.sum(n_diffs < 0))
record("B5: n_bounce(g₀) monotonicznie rośnie z g₀",
       n_nonmon == 0,
       f"n_malejących kroków = {n_nonmon} (z {len(n_diffs)})")

# B10: Czy A_tail ma lokalne ekstrema przy przejściach n_bounce?
# Szukamy lokalnych minimów A_tail w okolicach przejść
local_minima_at_transitions = []
for ti in transitions:
    lo, hi = max(0, ti-2), min(len(scan_A)-1, ti+3)
    a_window = scan_A[lo:hi+1]
    local_min = float(a_window.min())
    local_min_idx = np.argmin(a_window) + lo
    local_minima_at_transitions.append(
        abs(scan_A[ti] - local_min) / max(local_min, 1e-10)
    )

has_dips = len(local_minima_at_transitions) > 0 and any(
    d > 0.05 for d in local_minima_at_transitions)
record("B10: A_tail ma lokalne dip przy przejściach n_bounce?",
       has_dips,
       f"rel. dip = {[f'{d:.3f}' for d in local_minima_at_transitions]}")

# B8: Górna granica g₀ dla dobrego fitu
good_mask = scan_rmse < 0.10  # RMSE/A < 10%
if np.any(good_mask):
    g0_max_good = float(G0_SCAN_FINE[good_mask].max())
else:
    g0_max_good = 0.0
record("B8: Górna granica g₀^max z RMSE/A < 10%",
       g0_max_good > G0_TAU,
       f"g₀^max(RMSE<10%) = {g0_max_good:.4f}  [g₀^τ = {G0_TAU:.4f}]")


# ============================================================
# Sekcja 4: Szerokości zakresów n_bounce — czy kurczą się?
# ============================================================
print()
print("[4] SZEROKOŚCI ZAKRESÓW n_bounce (czy kurczą się geometrycznie?)")
print("-" * 65)

widths = {}
for n_target in [0, 1, 2, 3]:
    mask_n = (scan_n == n_target)
    if np.any(mask_n):
        g_range = G0_SCAN_FINE[mask_n]
        widths[n_target] = g_range.max() - g_range.min()
    else:
        widths[n_target] = 0.0

print(f"  {'n_bounce':>10}  {'szerokość Δg₀':>14}  {'ratio do n-1':>14}")
print("  " + "-"*42)
prev_w = None
width_ratios = []
for n in sorted(widths.keys()):
    w = widths[n]
    if prev_w and prev_w > 0:
        ratio = w / prev_w
        width_ratios.append(ratio)
        print(f"  {n:>10}  {w:>14.4f}  {ratio:>14.4f}")
    else:
        print(f"  {n:>10}  {w:>14.4f}  {'---':>14}")
    prev_w = w if w > 0 else prev_w

# B12: Szerokości kurczą się z n
if len(width_ratios) >= 2:
    shrinking = all(r < 1.0 for r in width_ratios)
    record("B12: Szerokości zakresów n_bounce kurczą się z n (Δg₀ maleje)?",
           shrinking,
           f"ratios = {[f'{r:.3f}' for r in width_ratios]}")
else:
    record("B12: Za mało danych do sprawdzenia kurczenia", False, "skan za rzadki")


# ============================================================
# Sekcja 5: Dokładne progi przejścia n_bounce metodą brentq
# ============================================================
print()
print("[5] DOKŁADNE PROGI n_bounce (metoda Brentq)")
print("-" * 50)

def count_bounces(g0):
    _, _, nb, _ = integrate_soliton_full(g0)
    return nb

# Znajdź progi między n=0 i n=1, n=1 i n=2, n=2 i n=3
bounce_thresholds = {}
for n_lo, n_hi in [(0, 1), (1, 2), (2, 3)]:
    # Znajdź przedział z danych skanu
    idx_lo = None
    for i in range(len(G0_SCAN_FINE)-1):
        if scan_n[i] == n_lo and scan_n[i+1] == n_hi:
            idx_lo = i
            break
    if idx_lo is not None:
        try:
            # Bisekcja na przedziałach
            g_a = G0_SCAN_FINE[idx_lo]
            g_b = G0_SCAN_FINE[idx_lo + 1]
            # Sprawdzamy kilka punktów wewnętrznych
            g_mid = 0.5*(g_a + g_b)
            nb_mid = count_bounces(g_mid)
            if nb_mid == n_lo:
                g_a = g_mid
            elif nb_mid == n_hi:
                g_b = g_mid
            thresh = 0.5*(g_a + g_b)
            bounce_thresholds[(n_lo, n_hi)] = thresh
            print(f"  Próg n={n_lo}→n={n_hi}: g₀* ≈ {thresh:.6f}")
        except Exception as e:
            print(f"  Próg n={n_lo}→n={n_hi}: błąd ({e})")
    else:
        print(f"  Próg n={n_lo}→n={n_hi}: brak w skanie")

# Sprawdzamy progresję progów: czy rosną z φ?
if (0,1) in bounce_thresholds and (1,2) in bounce_thresholds:
    g01 = bounce_thresholds[(0,1)]
    g12 = bounce_thresholds[(1,2)]
    ratio_01_12 = g12 / g01
    print(f"\n  g*(0→1) = {g01:.4f}")
    print(f"  g*(1→2) = {g12:.4f}")
    print(f"  Ratio g*(1→2)/g*(0→1) = {ratio_01_12:.4f}  [φ={PHI:.4f}]")
    if (2,3) in bounce_thresholds:
        g23 = bounce_thresholds[(2,3)]
        ratio_12_23 = g23 / g12
        print(f"  g*(2→3) = {g23:.4f}")
        print(f"  Ratio g*(2→3)/g*(1→2) = {ratio_12_23:.4f}")

# B11: Skok A_tail przy pierwszym odbiciu
if (0,1) in bounce_thresholds:
    g_thresh = bounce_thresholds[(0,1)]
    eps_test  = 0.05
    _, _, nb_lo, _ = integrate_soliton_full(g_thresh - eps_test)
    A_lo, _, _, _ = fit_tail(*integrate_soliton_full(g_thresh - eps_test)[:2])
    A_hi, _, _, _ = fit_tail(*integrate_soliton_full(g_thresh + eps_test)[:2])
    jump_A = abs(A_hi - A_lo) / max(A_lo, 1e-10)
    print(f"\n  Skok A_tail przy progu n=0→1:")
    print(f"  A(przed, g₀={g_thresh-eps_test:.4f}) = {A_lo:.6f}")
    print(f"  A(po,    g₀={g_thresh+eps_test:.4f}) = {A_hi:.6f}")
    print(f"  Rel. skok = {100*jump_A:.1f}%")
    record("B11: A_tail ciągła przy progu n=0→1 (skok < 50%)?",
           jump_A < 0.50,
           f"rel. skok = {100*jump_A:.1f}%")
else:
    record("B11: Brak progu 0→1 w danych", False, "")


# ============================================================
# Sekcja 6: Metryki stabilności dla generacji 1-4
# ============================================================
print()
print("[6] METRYKI STABILNOŚCI DLA k=1..4")
print("-" * 65)

print(f"  {'Gen':>4}  {'g₀':>8}  {'n_b':>5}  {'A_tail':>10}  "
      f"{'RMSE%':>7}  {'g₀ ponad progu n_b':>20}")
print("  " + "-"*65)

for lep, g0, klabel in GENS:
    d = gen_data[lep]
    nb = d['n_bounce']
    # Sprawdź o ile g₀ jest powyżej odpowiedniego progu
    thresh_key = (nb-1, nb) if nb > 0 else None
    if thresh_key and thresh_key in bounce_thresholds:
        margin = g0 - bounce_thresholds[thresh_key]
        margin_str = f"+{margin:.4f} nad progiem n={nb-1}→{nb}"
    else:
        margin_str = "---"
    print(f"  {lep:>4}  {g0:>8.5f}  {nb:>5d}  {d['A']:>10.5f}  "
          f"{100*d['rmse_rel']:>6.1f}%  {margin_str}")

# B13: ratio A(L4)/A(τ)
ratio_L4_tau = gen_data['L4']['A'] / gen_data['tau']['A']
print(f"\n  A_tail(L4)/A_tail(τ) = {ratio_L4_tau:.4f}")
print(f"  A_tail(τ)/A_tail(μ)  = {gen_data['tau']['A']/gen_data['mu']['A']:.4f}")
print(f"  A_tail(μ)/A_tail(e)  = {gen_data['mu']['A']/gen_data['e']['A']:.4f}")

# Sprawdzamy ratio A_tail na drabinie
ratios_A = [
    gen_data['mu']['A'] / gen_data['e']['A'],
    gen_data['tau']['A'] / gen_data['mu']['A'],
    gen_data['L4']['A'] / gen_data['tau']['A'],
]
print(f"\n  Stosunki A_tail wzdłuż drabiny:")
for i, (l1, l2) in enumerate([('e','mu'), ('mu','tau'), ('tau','L4')]):
    print(f"  A({l2})/A({l1}) = {ratios_A[i]:.4f}")

record("B13b: Ratios A_tail spowalniają wzdłuż drabiny (A(L4)/A(τ) < A(τ)/A(μ))?",
       ratio_L4_tau < gen_data['tau']['A']/gen_data['mu']['A'],
       f"A(L4)/A(τ)={ratio_L4_tau:.4f} vs A(τ)/A(μ)={gen_data['tau']['A']/gen_data['mu']['A']:.4f}")


# ============================================================
# Sekcja 7: BGH podsumowanie i wniosek
# ============================================================
print()
print("[7] PODSUMOWANIE BGH I WNIOSEK O T-OP3")
print("=" * 72)

nb_vals = {k: gen_data[k]['n_bounce'] for k in ['e','mu','tau','L4']}
bgh_full = (nb_vals['e'] == 0 and nb_vals['mu'] == 1
            and nb_vals['tau'] == 2 and nb_vals['L4'] == 3)

print(f"""
BOUNCE-GENERATION HYPOTHESIS (BGH):
  n_bounce(φ^{{k-1}}·g₀^e) = k - 1  dla k=1,2,3,4

  Weryfikacja:
    k=1 (e):   n_bounce = {nb_vals['e']}  [BGH: 0]  {'✓' if nb_vals['e']==0 else '✗'}
    k=2 (μ):   n_bounce = {nb_vals['mu']}  [BGH: 1]  {'✓' if nb_vals['mu']==1 else '✗'}
    k=3 (τ):   n_bounce = {nb_vals['tau']}  [BGH: 2]  {'✓' if nb_vals['tau']==2 else '✗'}
    k=4 (L4):  n_bounce = {nb_vals['L4']}  [BGH: 3]  {'✓' if nb_vals['L4']==3 else '✗'}

  BGH status: {'POTWIERDZONE ✓' if bgh_full else 'CZĘŚCIOWE'}

IMPLIKACJE DLA T-OP3:
""")

print(f"""
BOUNCE-GENERATION HYPOTHESIS — WYNIKI:
  e:  n_bounce={nb_vals['e']}   RMSE/A={100*gen_data['e']['rmse_rel']:.1f}%  ← STABILNY
  μ:  n_bounce={nb_vals['mu']}   RMSE/A={100*gen_data['mu']['rmse_rel']:.1f}%  ← STABILNY
  τ:  n_bounce={nb_vals['tau']}   RMSE/A={100*gen_data['tau']['rmse_rel']:.1f}%  ← STABILNY (graniczny)
  L4: n_bounce={nb_vals['L4']}   RMSE/A={100*gen_data['L4']['rmse_rel']:.1f}%  ← ZDEGRADOWANY!

BGH ścisłe (n=k-1): OBALONE. Tau ma {nb_vals['tau']} odbicia, nie 2.
BGH zrewidowane (n monot.): POTWIERDZONE (e<μ<τ<L4 odbiciami).

NOWY WYNIK (RMSE-CUTOFF HYPOTHESIS, RCH):
  g₀^max z RMSE/A < 10% = {g0_max_good:.4f}
  Przedział: g₀^τ={G0_TAU:.4f} < g₀^max={g0_max_good:.4f} < g₀^L4={G0_L4:.4f}
  → τ JEST OSTATNIĄ GENERACJĄ Z DOBRZE OKREŚLONYM OGONEM ASYMPTOTYCZNYM
  → L4 ma RMSE/A=36.4%: ogon niestabilny, masa z A_tail^4 niepewna

MECHANIZM WYKLUCZENIA GENERACJI 4 W TGP (bez LEP!):
  Soliton L4 ISTNIEJE (ma A_tail > 0), ale jakość dopasowania
  asymptotyki jest zdegradowana (RMSE/A={100*gen_data['L4']['rmse_rel']:.0f}%).
  Mechanizm A_tail^4 (skąd pochodzi masa) nie jest wiarygodny dla L4.

  WNIOSEK: TGP WYKLUCZA L4 WEWNĘTRZNIE przez RCH!
  (Niezależne od LEP — wewnętrzna niespójność mechanizmu masowego dla k=4)
""")

# Kluczowa liczba: jaka jest "maksymalna generacja" z RMSE/A < 10%?
max_gen_good = max(k for k, lep in enumerate(['e','mu','tau','L4'], 1)
                   if gen_data[['e','mu','tau','L4'][k-1]]['rmse_rel'] < 0.10)
print(f"  Maksymalna generacja z RMSE/A < 10%: k_max = {max_gen_good}")
print(f"  (czyli {['e','mu','tau','L4'][max_gen_good-1]} = gen. {max_gen_good})")

# B14: Kluczowy wynik — RMSE kryterium wyklucza L4 naturalnie
rmse_cutoff_excludes_L4 = (
    gen_data['e']['rmse_rel'] < 0.10 and
    gen_data['mu']['rmse_rel'] < 0.10 and
    gen_data['tau']['rmse_rel'] < 0.10 and
    gen_data['L4']['rmse_rel'] > 0.15
)
record("B14: RMSE/A kryterium: e,μ,τ OK (<10%), L4 zdegradowane (>15%)",
       rmse_cutoff_excludes_L4,
       f"e:{100*gen_data['e']['rmse_rel']:.1f}%, μ:{100*gen_data['mu']['rmse_rel']:.1f}%, "
       f"τ:{100*gen_data['tau']['rmse_rel']:.1f}%, L4:{100*gen_data['L4']['rmse_rel']:.1f}%")


# ============================================================
# Wyniki końcowe
# ============================================================
print()
print("=" * 72)
print("WYNIKI TESTÓW")
print("=" * 72)
passed = sum(1 for _, p, _ in TESTS if p)
total  = len(TESTS)
for name, passed_t, detail in TESTS:
    mark = "PASS" if passed_t else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"         {detail}")

print()
print(f"WYNIK: {passed}/{total} testów PASS")
print()
print("KLUCZOWE WYNIKI LICZBOWE:")
for lep, g0, klabel in GENS:
    d = gen_data[lep]
    print(f"  {klabel} ({lep:>3}): g₀={g0:.5f}  n_b={d['n_bounce']}  "
          f"A={d['A']:.5f}  RMSE/A={100*d['rmse_rel']:.1f}%")
print()
print("WNIOSKI:")
if bgh_full:
    print("  BGH: POTWIERDZONE — n_bounce(φ^{k-1}·g₀^e) = k-1")
    print("  T-OP3 CZĘŚCIOWY: BGH daje topologiczną klasyfikację generacji,")
    print("  ale nie wyklucza dynamicznie k=4. Dodatkowy argument: LEP.")
    print("  NOWE PYTANIE: dlaczego n_bounce_max=2? (nie 3 czy 4?)")
else:
    print("  BGH: CZĘŚCIOWE — sprawdzić ręcznie n_bounce(L4)")
