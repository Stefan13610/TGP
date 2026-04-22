"""
EX69: ZASADA SELEKCJI g0* — zaleznosc od alpha i warunki kwantowania
STATUS: LEGACY-TRANSLATIONAL

This file is part of the older `g0* ~ 1.24`, `alpha*`, and quantization
selection program. It remains useful as a record of historical hypotheses, but
it is not aligned with the canonical synchronized `nbody` framing.

Pytanie: Co wyznacza g0* = 1.24771 z aksjomatow TGP (bez wejscia r21)?

Strategia: analizujemy f(g0, alpha) = (A(phi*g0, alpha) / A(g0, alpha))^4
           i szukamy warunku, ktory przy alpha=2 daje g0* = 1.24771.

Badane hipotezy:
  H1: g0* = z0(alpha) — punkt B-S zero (B_tail=0)
      → przy alpha=2: z0=1.2301 ≠ g0*=1.2477 (FAIL)
      → czy przy alpha=alpha* z0(alpha*) = g0*(alpha*)?

  H2: g0*(alpha) = z0(alpha) + delta(alpha)
      gdzie delta jest znana funkcja alpha

  H3: g0* jest minimum/maximum/punklem przegięcia f(g0, alpha)

  H4: g0*(alpha) * alpha = stala (kwantowanie)

  H5: g0*(alpha=2) = z0(alpha=2) + (alpha*-alpha) * dg0/dalpha
      (perturbacyjna poprawka)

Tests T1-T5:
  T1: g0*(alpha) vs z0(alpha) dla alpha in {2.0, 2.2, 2.44} — roznica
  T2: f'(g0*, alpha=2) vs f'(z0, alpha=2) — gdzie f szybciej sie zmienia?
  T3: g0*(alpha) * alpha = stala? (sprawdz alpha*g0* dla kilku alpha)
  T4: Ekstrapolacja g0*(alpha) -> g0*(alpha=2) z danych alpha in {2.2, 2.44}
  T5: Czy g0*(alpha) = z0(alpha) dla alpha=alpha*? (unifikacja)
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from scipy.stats import linregress
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# STALE
# ============================================================
PHI     = (1.0 + np.sqrt(5.0)) / 2.0
R21_EXP = 206.768
G_OFF   = 0.005
R_MAX   = 100.0
WIN_LIST  = [16, 22, 28, 36, 46, 58, 72]
WIN_WIDTH = 14.0

# ============================================================
# ODE i solver (z ex62b — z odbiciami)
# ============================================================
def make_ode(alpha, g_off=G_OFF):
    g_bounce = np.exp(-1.0/(2.0*alpha)) + g_off
    def rhs(r, y):
        g, gp = y
        g = max(g, g_bounce + 1e-7)
        fg = 1.0 + 2.0*alpha*np.log(g)
        if abs(fg) < 1e-10: return [gp, 0.0]
        dr = g**2*(1.0-g); cr = (alpha/g)*gp**2
        if r < 1e-10: return [gp, (dr-cr)/(3.0*fg)]
        return [gp, (dr-cr-fg*2.0*gp/r)/fg]
    def ev(r, y): return y[0] - (np.exp(-1.0/(2.0*alpha)) + g_off)
    ev.terminal = True; ev.direction = -1
    return rhs, ev

def integrate(g0, alpha=2.0, g_off=G_OFF, maxb=8):
    rhs, ev = make_ode(alpha, g_off)
    y0 = [g0, 0.0]; r0 = 1e-10; ra, ga = [], []
    for _ in range(maxb+1):
        sol = solve_ivp(rhs, (r0, R_MAX), y0, events=ev,
                        dense_output=True, rtol=1e-9, atol=1e-11, max_step=0.05)
        ra.append(sol.t); ga.append(sol.y[0])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]; st = sol.sol(rh)
            y0 = [st[0], -st[1]]; r0 = rh
        else: break
    r = np.concatenate(ra); g = np.concatenate(ga)
    idx = np.argsort(r)
    return r[idx], g[idx]

def fit_win(r, g, rL, rR):
    mask = (r >= rL) & (r <= rR)
    if np.sum(mask) < 20: return 0., 0., 0.
    rf = r[mask]; df = (g[mask]-1.0)*rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    B, C = bc
    return float(np.sqrt(B**2+C**2)), float(B), float(C)

def A_inf(g0, alpha=2.0):
    r, g = integrate(g0, alpha)
    Av, rv = [], []
    for rL in WIN_LIST:
        if rL + WIN_WIDTH > r[-1]: break
        A, _, _ = fit_win(r, g, rL, rL+WIN_WIDTH)
        if A > 0.002: Av.append(A); rv.append(float(rL))
    if not Av: return 0.
    if len(Av) < 3: return Av[-1]
    try:
        p, _ = curve_fit(lambda x, ai, a: ai*(1+a/x), rv, Av, p0=[Av[-1], 0.], maxfev=2000)
        return float(p[0])
    except: return float(Av[-1])

def B_coeff(g0, alpha=2.0):
    r, g = integrate(g0, alpha)
    _, B, _ = fit_win(r, g, WIN_LIST[2], WIN_LIST[2]+WIN_WIDTH)
    return B

def find_z0(alpha=2.0, lo=1.05, hi=1.60):
    gs = np.linspace(lo, hi, 60)
    Bs = [B_coeff(g, alpha) for g in gs]
    for i in range(len(Bs)-1):
        if Bs[i]*Bs[i+1] < 0:
            try: return brentq(lambda x: B_coeff(x, alpha), gs[i], gs[i+1], xtol=1e-7)
            except: pass
    return None

def f_ratio(g0, alpha=2.0):
    Ae = A_inf(g0, alpha)
    Am = A_inf(PHI*g0, alpha)
    return (Am/Ae)**4 if Ae > 1e-6 else np.nan

def find_g0star(alpha=2.0, target=R21_EXP, lo=1.10, hi=1.60):
    """Znajdz g0* takie ze f_ratio(g0*, alpha) = target."""
    gs = np.linspace(lo, hi, 30)
    fs = [f_ratio(g, alpha) for g in gs]
    # Szukaj przejscia przez target
    for i in range(len(fs)-1):
        if not (np.isnan(fs[i]) or np.isnan(fs[i+1])):
            if (fs[i]-target)*(fs[i+1]-target) < 0:
                try:
                    return brentq(lambda x: f_ratio(x, alpha) - target, gs[i], gs[i+1], xtol=1e-7)
                except: pass
    return None

# ============================================================
# CZESC 1: f(g0, alpha) dla alpha = 2.0, 2.2, 2.44
# ============================================================
print("=" * 68)
print("EX69: ZASADA SELEKCJI g0* — ZALEZNOSC OD ALPHA")
print("=" * 68)

ALPHAS = [2.0, 2.20, 2.44]

print("\n" + "-" * 68)
print("--- Czesc 1: z0(alpha), g0*(alpha), f(g0*, alpha) ---")
print()

z0_tab    = {}
g0st_tab  = {}
f_at_z0   = {}
f_at_g0st = {}

print(f"  {'alpha':>6}  {'z0':>10}  {'g0*':>10}  {'g0*-z0':>8}  "
      f"{'f(z0)':>10}  {'f(g0*)':>10}  {'g0*-z0 as frac z0':>18}")
print(f"  {'-'*6}  {'-'*10}  {'-'*10}  {'-'*8}  "
      f"{'-'*10}  {'-'*10}  {'-'*18}")

for alpha in ALPHAS:
    z0   = find_z0(alpha)
    g0st = find_g0star(alpha)
    fz0  = f_ratio(z0, alpha) if z0 else np.nan
    fg0s = f_ratio(g0st, alpha) if g0st else np.nan
    z0_tab[alpha]    = z0
    g0st_tab[alpha]  = g0st
    f_at_z0[alpha]   = fz0
    f_at_g0st[alpha] = fg0s
    diff  = (g0st - z0) if (g0st and z0) else np.nan
    frac  = diff / z0 if z0 else np.nan
    g0st_s = f"{g0st:.7f}" if g0st else "N/A"
    z0_s   = f"{z0:.7f}"   if z0   else "N/A"
    fz0_s  = f"{fz0:.3f}"  if not np.isnan(fz0) else "N/A"
    fg0s_s = f"{fg0s:.3f}" if not np.isnan(fg0s) else "N/A"
    diff_s = f"{diff:.6f}" if not np.isnan(diff) else "N/A"
    frac_s = f"{frac:.5f}" if not np.isnan(frac) else "N/A"
    print(f"  {alpha:6.2f}  {z0_s:>10}  {g0st_s:>10}  {diff_s:>8}  "
          f"{fz0_s:>10}  {fg0s_s:>10}  {frac_s:>18}")

# ============================================================
# CZESC 2: Mapa f(g0) dla alpha=2: szukaj specjalnych wlasnosci w g0*
# ============================================================
print("\n" + "-" * 68)
print("--- Czesc 2: f(g0) i df/dg0 w okolicach g0* (alpha=2) ---")
print()

alpha_ref = 2.0
z0_ref   = z0_tab.get(alpha_ref, 1.2301)
g0st_ref = g0st_tab.get(alpha_ref, 1.24771)

g0_scan = np.linspace(1.18, 1.32, 29)
f_scan  = [f_ratio(g, alpha_ref) for g in g0_scan]

print(f"  {'g0':>8}  {'f(g0)':>10}  {'f-r21':>10}")
print(f"  {'-'*8}  {'-'*10}  {'-'*10}")
for g0, f in zip(g0_scan, f_scan):
    marker = ""
    if z0_ref and abs(g0-z0_ref) < 0.01: marker = " <z0"
    if g0st_ref and abs(g0-g0st_ref) < 0.01: marker = " <g0*"
    f_str = f"{f:.4f}" if not np.isnan(f) else "N/A"
    fr_str = f"{f-R21_EXP:.4f}" if not np.isnan(f) else "N/A"
    print(f"  {g0:8.4f}  {f_str:>10}  {fr_str:>10}{marker}")

# Numeryczne pochodne f'(g0) i f''(g0)
f_arr = np.array(f_scan)
dg    = g0_scan[1] - g0_scan[0]
valid = ~np.isnan(f_arr)
if valid.sum() >= 5:
    # Numeryczna pierwsza pochodna
    df_arr = np.gradient(f_arr, dg)
    ddf_arr = np.gradient(df_arr, dg)

    # Znajdz punkty gdzie df/dg0 zmienia znak lub jest ekstremalny
    sign_changes = []
    for i in range(len(df_arr)-1):
        if df_arr[i]*df_arr[i+1] < 0:
            sign_changes.append((g0_scan[i]+g0_scan[i+1])/2)

    print(f"\n  Zmiany znaku df/dg0 (ekstrema f): {sign_changes}")

    # Wartosc df/dg0 w z0 i g0*
    idx_z0  = np.argmin(abs(g0_scan - z0_ref))  if z0_ref  else None
    idx_g0s = np.argmin(abs(g0_scan - g0st_ref)) if g0st_ref else None

    if idx_z0 is not None:
        print(f"  df/dg0 @ z0={z0_ref:.4f}:  {df_arr[idx_z0]:.4f}")
        print(f"  d2f/dg0^2 @ z0:           {ddf_arr[idx_z0]:.4f}")
    if idx_g0s is not None:
        print(f"  df/dg0 @ g0*={g0st_ref:.5f}: {df_arr[idx_g0s]:.4f}")
        print(f"  d2f/dg0^2 @ g0*:           {ddf_arr[idx_g0s]:.4f}")

# ============================================================
# CZESC 3: Hipoteza H4: g0*(alpha) * alpha = stala?
# ============================================================
print("\n" + "-" * 68)
print("--- Czesc 3: Hipoteza g0*(alpha)*alpha = const? ---")
print()

products = {}
for alpha in ALPHAS:
    g0st = g0st_tab.get(alpha)
    if g0st:
        products[alpha] = g0st * alpha
        print(f"  alpha = {alpha:.2f}: g0* = {g0st:.7f}, g0**alpha = {products[alpha]:.6f}")

if len(products) >= 2:
    vals = list(products.values())
    rng  = max(vals) - min(vals)
    print(f"\n  Rozrzut g0**alpha: {rng:.6f} (wzgledny: {rng/np.mean(vals)*100:.2f}%)")
    H4 = rng / np.mean(vals) < 0.02
    print(f"  H4 (rng/mean < 2%): {'PASS' if H4 else 'FAIL'}")

# ============================================================
# CZESC 4: Hipoteza H2: delta(alpha) = g0*(alpha) - z0(alpha) vs alpha
# ============================================================
print("\n" + "-" * 68)
print("--- Czesc 4: delta(alpha) = g0*(alpha) - z0(alpha) ---")
print()

deltas = {}
for alpha in ALPHAS:
    z0   = z0_tab.get(alpha)
    g0st = g0st_tab.get(alpha)
    if z0 and g0st:
        d = g0st - z0
        deltas[alpha] = d
        print(f"  alpha = {alpha:.2f}: z0 = {z0:.6f}, g0* = {g0st:.6f}, delta = {d:.6f}")

if len(deltas) >= 2:
    al_arr = np.array(list(deltas.keys()))
    d_arr  = np.array(list(deltas.values()))
    slope_d, ic_d, r2_d, _, _ = linregress(al_arr, d_arr)
    print(f"\n  Fit delta(alpha) = a*alpha + b: a = {slope_d:.6f}, b = {ic_d:.6f}, R^2 = {r2_d**2:.5f}")
    alpha_pred = -ic_d / slope_d if slope_d != 0 else np.nan
    print(f"  Zero delta @ alpha = {alpha_pred:.4f}  (oczekiwane: alpha* ≈ 2.44)")

# ============================================================
# CZESC 5: g0*(alpha=2) obliczone roznymi metodami
# ============================================================
print("\n" + "-" * 68)
print("--- Czesc 5: Rne wyrazenia dla g0* ---")
print()

g0s_alpha2 = g0st_tab.get(2.0, 1.24771102)
z0_alpha2  = z0_tab.get(2.0,  1.2301)

print(f"  g0* (alpha=2)     = {g0s_alpha2:.7f}  (z ex62b)")
print(f"  z0  (alpha=2)     = {z0_alpha2:.7f}  (B=0 warunek)")
print(f"  g0* - z0          = {g0s_alpha2 - z0_alpha2:.7f}")
print(f"  g0* / z0          = {g0s_alpha2/z0_alpha2:.7f}")
print(f"  ln(g0*/z0)        = {np.log(g0s_alpha2/z0_alpha2):.7f}")
print(f"  g0* / phi         = {g0s_alpha2/PHI:.7f}")
print(f"  g0* * phi         = {g0s_alpha2*PHI:.7f}")
print(f"  g0* * phi^2       = {g0s_alpha2*PHI**2:.7f}")
print(f"  1/g0*             = {1/g0s_alpha2:.7f}")
print(f"  1/g0* - 1         = {1/g0s_alpha2 - 1:.7f}")
print(f"  exp(1/g0*)        = {np.exp(1/g0s_alpha2):.7f}")
print(f"  g0*^2             = {g0s_alpha2**2:.7f}")
print(f"  g0*^pi            = {g0s_alpha2**np.pi:.7f}")
print(f"  g0* * 2 * pi      = {g0s_alpha2 * 2 * np.pi:.7f}")
print(f"  g0* / z0 - 1      = {g0s_alpha2/z0_alpha2 - 1:.7f}")

# Sprawdz A_inf(g0*)
A_e = A_inf(g0s_alpha2)
A_z0 = A_inf(z0_alpha2)
print(f"\n  A_inf(g0*)        = {A_e:.7f}")
print(f"  A_inf(z0)         = {A_z0:.7f}")
print(f"  A_inf(g0*)/A_inf(z0) = {A_e/A_z0:.7f}")
print(f"  A_inf(g0*)/A_inf(z0) - 1 = {A_e/A_z0-1:.7f}")

# ============================================================
# TESTY
# ============================================================
print("\n" + "=" * 68)
print("TESTY")
print("=" * 68)

# T1: g0*(alpha) vs z0(alpha) — roznica vs alpha
deltas_v = list(deltas.values())
T1_varz0g0 = any(abs(d) > 0.01 for d in deltas_v)   # czy g0* ≠ z0 dla alpha < alpha*
T1_z0_at_alpha_star = False
if len(deltas) >= 2:
    # Czy delta(alpha*) ≈ 0?
    alpha_star_pred = -ic_d/slope_d if slope_d != 0 else np.nan
    T1_z0_at_alpha_star = not np.isnan(alpha_star_pred) and 2.3 < alpha_star_pred < 2.6
    T1 = T1_varz0g0 and T1_z0_at_alpha_star
else:
    T1 = False
print(f"  T1: g0* ≠ z0 dla alpha<alpha*, delta->0 @ alpha*:  {'PASS' if T1 else 'FAIL'}  "
      f"(delta(alpha=2) = {deltas.get(2.0, np.nan):.5f})")

# T2: f(z0) vs f(g0*) — jaka roznica?
fz0_2  = f_at_z0.get(2.0, np.nan)
fg0s_2 = f_at_g0st.get(2.0, np.nan)
T2 = not np.isnan(fz0_2) and not np.isnan(fg0s_2) and abs(fg0s_2 - R21_EXP)/R21_EXP < 0.001
print(f"  T2: f(g0*) = r21 (0.1%):  {'PASS' if T2 else 'FAIL'}  "
      f"(f(z0)={fz0_2:.2f}, f(g0*)={fg0s_2:.3f})")

# T3: g0*alpha = const?
T3 = H4 if 'H4' in dir() else False
print(f"  T3: g0*(alpha) * alpha = const:  {'PASS' if T3 else 'FAIL'}  "
      f"({list(products.values()) if products else 'N/A'})")

# T4: delta(alpha) -> 0 @ alpha* ≈ 2.44
T4 = T1_z0_at_alpha_star
alpha_pred_str = f"{-ic_d/slope_d:.4f}" if (slope_d != 0 and len(deltas) >= 2) else "N/A"
print(f"  T4: delta=0 @ alpha* in [2.3,2.6]:  {'PASS' if T4 else 'FAIL'}  "
      f"(alpha_pred = {alpha_pred_str})")

# T5: czy g0*(alpha*) ≈ z0(alpha*)
alpha_star_alpha = 2.44
z0_astar  = z0_tab.get(alpha_star_alpha,  find_z0(alpha_star_alpha))
g0s_astar = g0st_tab.get(alpha_star_alpha, find_g0star(alpha_star_alpha))
T5 = False
if z0_astar and g0s_astar:
    diff5 = abs(g0s_astar - z0_astar)
    T5 = diff5 < 0.01
    print(f"  T5: g0*(alpha*)  ≈ z0(alpha*) (diff<0.01):  {'PASS' if T5 else 'FAIL'}  "
          f"(z0={z0_astar:.5f}, g0*={g0s_astar:.5f}, diff={diff5:.5f})")
else:
    print(f"  T5: g0*(alpha*)  ≈ z0(alpha*):  FAIL  (brak danych)")

n_pass = sum([T1, T2, T3, T4, T5])
print(f"\nWYNIK: {n_pass}/5 testow przeszlo")

print("\n" + "=" * 68)
print("WNIOSEK (EX69)")
print("=" * 68)
print(f"  z0(alpha=2)   = {z0_tab.get(2.0, 'N/A')}")
print(f"  g0*(alpha=2)  = {g0st_tab.get(2.0, 'N/A')}")
print(f"  delta = g0*-z0 = {deltas.get(2.0, np.nan):.6f}")
print(f"  f(z0) = {f_at_z0.get(2.0, np.nan):.4f}  vs  f(g0*) = {f_at_g0st.get(2.0, np.nan):.4f}")
if len(deltas) >= 2:
    print(f"  delta->0 @ alpha = {-ic_d/slope_d:.4f}  (alpha* = 2.441)")
print(f"  Konkluzja: g0*(alpha*)={g0s_astar}, z0(alpha*)={z0_astar}")
print("=" * 68)
