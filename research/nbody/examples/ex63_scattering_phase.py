"""
EX63: FAZA SCATTERINGOWA delta(g0) I WARUNKI KWANTOWANIA
=========================================================
Cel: zrozumiec, jakie kryterium fizyczne wyznacza g0* = 1.24771.

Hipoteza: warunek kwantowania TGP to para faz:
  delta_scatt(g0_e) = phi_e  (wartsc "podstawowa")
  delta_scatt(g0_mu) = phi_mu
  gdzie relacja phi_e, phi_mu wynika z aksjomatu TGP.

Pytania:
  1. Jak wyglada delta_scatt(g0) w zakresie n=0 i n=1 sektorow?
     (ciagla mapa, np.unwrap)

  2. Przy jakim g0 jest delta_scatt = 0, pi/4, pi/2, pi, 3pi/2, 2pi?
     (Bohr-Sommerfeld dla roznych n)

  3. Czy para (g0_e, g0_mu) spelnia warunek:
     delta_scatt(g0_mu) - delta_scatt(g0_e) = pi?
     Lub: delta_scatt(g0_mu) = n*pi?

  4. Jaka jest faza scatteringowa przy g0* = 1.24771?
     Czy to liczba wyrozniajaca sie fizycznie?

  5. "Rezonansowe" kryterium: g0 gdzie A_tail jest ekstremalny lokalnie?
     (czy to pokrywa sie z g0*?)

  6. Czy istnieje warunek "samodualna" faza:
     delta(g0_mu) = pi - delta(g0_e)?

Testy:
  T1: |delta_scatt(z0)| < 0.05 rad  (H1: B=0 => delta=0)
  T2: |delta_scatt(phi*z0) - pi| < 0.5 rad  (B-S muon?)
  T3: |delta_scatt(g0*) - delta_BS| < 0.1 rad  (g0* ma sens fazowy)
  T4: delta_scatt ciagla (max jump < 0.3 po unwrap) w [1.10, 2.20]
  T5: (A(g0_mu_BS)/A(g0_e_BS))^4 in [150, 300]  (B-S para daje blisko 207?)
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
import warnings
warnings.filterwarnings('ignore')

# ============================================================
PHI      = (1.0 + np.sqrt(5.0)) / 2.0
R21_EXP  = 206.768
R_MAX    = 100.0
WIN_LIST  = [16, 22, 28, 36, 46, 58, 72]
WIN_WIDTH = 14.0
G0_STAR  = 1.24771102   # phi-FP z ex62

# ============================================================
# ODE helpers
# ============================================================

ALPHA = 2.0
G_GHOST  = np.exp(-1.0 / (2.0 * ALPHA))   # = 0.77880
G_BOUNCE = G_GHOST + 0.005                 # = 0.78380

def rhs(r, y):
    g, gp = y
    g = max(g, G_BOUNCE + 1e-7)
    fg = 1.0 + 2.0 * ALPHA * np.log(g)
    if abs(fg) < 1e-10: return [gp, 0.0]
    dr = g**2 * (1.0 - g)
    cr = (ALPHA / g) * gp**2
    if r < 1e-10: return [gp, (dr - cr) / (3.0 * fg)]
    return [gp, (dr - cr - fg * 2.0 * gp / r) / fg]

def ev_ghost(r, y): return y[0] - G_BOUNCE
ev_ghost.terminal  = True
ev_ghost.direction = -1


def integrate(g0, max_bounces=8):
    y0 = [g0, 0.0]; r0 = 0.0; ra, ga = [], []
    for _ in range(max_bounces + 1):
        sol = solve_ivp(rhs, (r0, R_MAX), y0, events=ev_ghost,
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
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M  = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    B, C = bc
    return float(np.sqrt(B**2+C**2)), float(B), float(C)


def A_inf(g0):
    r, g = integrate(g0)
    Av, rv = [], []
    for rL in WIN_LIST:
        if rL + WIN_WIDTH > r[-1]: break
        A, _, _ = fit_win(r, g, rL, rL+WIN_WIDTH)
        if A > 0.002: Av.append(A); rv.append(rL)
    if len(Av) < 3: return Av[-1] if Av else 0.
    try:
        p, _ = curve_fit(lambda x, ai, a: ai*(1+a/x), rv, Av,
                         p0=[Av[-1], 0.], maxfev=2000)
        return float(p[0])
    except: return float(Av[-1])


def delta_scatt(g0, rL=WIN_LIST[2], rW=WIN_WIDTH):
    """Faza scatteringowa delta = atan2(B, C) w oknie [rL, rL+rW]."""
    r, g = integrate(g0)
    A, B, C = fit_win(r, g, rL, rL+rW)
    if A < 1e-5: return np.nan
    return np.arctan2(B, C)


# ============================================================
print("="*68)
print("EX63: FAZA SCATTERINGOWA delta(g0) I WARUNKI KWANTOWANIA")
print("="*68)
print(f"  phi = {PHI:.7f},  g0* = {G0_STAR:.7f}")
print(f"  G_GHOST = {G_GHOST:.5f},  G_BOUNCE = {G_BOUNCE:.5f}")
print()

# ============================================================
# CZESC 1: Mapa delta_scatt(g0) — gestskan
# ============================================================
print("-"*68)
print("--- Czesc 1: Mapa delta_scatt(g0) i A_inf(g0) ---")
print()

N1   = 80
g0_1 = np.linspace(1.10, 2.30, N1)
delta_1 = []
A_1     = []
n_bounce_1 = []

for g0 in g0_1:
    r, g = integrate(g0)
    # Policz n_bounce: ile razy przerywalo
    # (prosta heurystyka przez liczbe minimow g < G_BOUNCE+0.01)
    below = np.where(g < G_BOUNCE + 0.03)[0]
    nb = 0
    if len(below) > 0:
        # policz segmenty
        prev = below[0]
        nb = 1
        for idx in below[1:]:
            if idx > prev + 20:
                nb += 1
            prev = idx
    n_bounce_1.append(nb)
    A_1.append(A_inf(g0))
    delta_1.append(delta_scatt(g0))
    print(f"  g0={g0:.4f}  nb={nb}  A={A_1[-1]:.5f}  delta={delta_1[-1]:.4f} rad  ({delta_1[-1]/np.pi:.4f}*pi)")

delta_arr = np.array(delta_1)
delta_uw  = np.unwrap(delta_arr)
A_arr     = np.array(A_1)
n_b_arr   = np.array(n_bounce_1)

print()
print(f"  Zakres delta_raw: [{np.nanmin(delta_arr):.4f}, {np.nanmax(delta_arr):.4f}] rad")
print(f"  Zakres delta_uw:  [{np.nanmin(delta_uw):.4f}, {np.nanmax(delta_uw):.4f}] rad")
print()

# ============================================================
# CZESC 2: Faza przy kluczowych punktach
# ============================================================
print("-"*68)
print("--- Czesc 2: Faza przy kluczowych g0 ---")
print()

key_pts = [
    ("z0 (H1)",         1.2301),
    ("g0_e (ex57)",     1.2400),
    ("g0* (phi-FP)",    G0_STAR),
    ("phi*z0",          PHI*1.2301),
    ("phi*g0*",         PHI*G0_STAR),
    ("z0_mu (H1 n=1)",  2.0797),
    ("g0_e_BS?",        1.9726),
]

print(f"  {'Nazwa':20s}  {'g0':10s}  {'delta(rad)':12s}  {'delta/pi':10s}  {'A_inf':10s}")
print("  " + "-"*68)
for name, g0 in key_pts:
    d  = delta_scatt(g0)
    Ai = A_inf(g0)
    print(f"  {name:20s}  {g0:10.6f}  {d:12.6f}  {d/np.pi:10.6f}  {Ai:10.6f}")
print()

# ============================================================
# CZESC 3: Znajdz g0 gdzie delta = n*pi dla n=0,1,2,...
# ============================================================
print("-"*68)
print("--- Czesc 3: Bohr-Sommerfeld: g0 gdzie delta(g0) = n*pi ---")
print()

# Uzywamy delta_uw (unwrapped) dla ciaglej mapy
# Znajdz zerowania delta_uw - n*pi
print(f"  {'n':>3}  {'g0_BS':>10}  {'delta_BS(rad)':>14}  {'A_inf':>10}  {'sector n_b':>12}")
print("  " + "-"*58)
BS_points = {}
for n_target in range(7):
    target = n_target * np.pi
    for i in range(len(delta_uw)-1):
        if np.isnan(delta_uw[i]) or np.isnan(delta_uw[i+1]): continue
        if (delta_uw[i] - target) * (delta_uw[i+1] - target) < 0:
            # Interpolacja liniowa
            f = (target - delta_uw[i]) / (delta_uw[i+1] - delta_uw[i])
            g0_bs = g0_1[i] + f * (g0_1[i+1] - g0_1[i])
            # Weryfikacja
            d_bs  = delta_scatt(g0_bs)
            A_bs  = A_inf(g0_bs)
            nb_bs = n_b_arr[i]
            print(f"  {n_target:>3}  {g0_bs:10.6f}  {d_bs:14.6f}  {A_bs:10.6f}  {nb_bs:12d}")
            BS_points[n_target] = (g0_bs, d_bs, A_bs, nb_bs)
            break

print()

# ============================================================
# CZESC 4: Ratio (A_mu/A_e)^4 dla par B-S
# ============================================================
print("-"*68)
print("--- Czesc 4: Ratio A-ratio^4 dla par Bohr-Sommerfeld ---")
print()

print(f"  Para (n_e, n_mu): r21 = (A(g0_mu_BS)/A(g0_e_BS))^4")
print(f"  {'(ne,nmu)':10s}  {'g0_e':9s}  {'g0_mu':9s}  {'ratio^4':10s}  {'delta%':9s}")
print("  " + "-"*55)
for ne in range(4):
    for nmu in range(ne+1, 6):
        if ne in BS_points and nmu in BS_points:
            g0_e_bs, _, A_e_bs, _ = BS_points[ne]
            g0_mu_bs, _, A_mu_bs, _ = BS_points[nmu]
            if A_e_bs > 1e-6:
                ratio = (A_mu_bs / A_e_bs)**4
                delta = (ratio - R21_EXP) / R21_EXP * 100
                mark = " <===  !!!" if abs(delta) < 10 else ""
                print(f"  ({ne},{nmu}):        {g0_e_bs:9.5f}  {g0_mu_bs:9.5f}  {ratio:10.3f}  {delta:9.3f}%{mark}")

print()

# ============================================================
# CZESC 5: Faza przy g0* — interpretacja
# ============================================================
print("-"*68)
print("--- Czesc 5: Analiza fazy przy phi-FP g0* = 1.24771 ---")
print()

d_star     = delta_scatt(G0_STAR)
d_star_mu  = delta_scatt(PHI * G0_STAR)
A_star     = A_inf(G0_STAR)
A_star_mu  = A_inf(PHI * G0_STAR)

print(f"  delta_scatt(g0*)            = {d_star:.6f} rad = {d_star/np.pi:.6f} pi")
print(f"  delta_scatt(phi*g0*)        = {d_star_mu:.6f} rad = {d_star_mu/np.pi:.6f} pi")
print(f"  delta_mu - delta_e          = {d_star_mu - d_star:.6f} rad = {(d_star_mu-d_star)/np.pi:.6f} pi")
print(f"  |delta_mu - delta_e - pi|   = {abs(d_star_mu - d_star - np.pi):.6f} rad")
print(f"  |delta_mu - delta_e - pi/2| = {abs(d_star_mu - d_star - np.pi/2):.6f} rad")
print(f"  |delta_e + delta_mu - pi|   = {abs(d_star + d_star_mu - np.pi):.6f} rad")
print()

# A-ratio rekonstrukcja
ratio_FP = (A_star_mu / A_star)**4
print(f"  A_inf(g0*)    = {A_star:.7f}")
print(f"  A_inf(phi*g0*) = {A_star_mu:.7f}")
print(f"  ratio^4        = {ratio_FP:.6f}  (vs {R21_EXP})")
print()

# Dodatkowa kandydatura: delta_e = pi/6 (30 deg)?
print(f"  Specjalne wartosci delta_e:")
specials = [0, np.pi/8, np.pi/6, np.pi/4, np.pi/3, np.pi/2]
for s in specials:
    print(f"    delta={s/np.pi:.4f}*pi:  |delta(g0*) - {s:.4f}| = {abs(d_star - s):.6f}")
print()

# ============================================================
# CZESC 6: Ciagla relacja delta(phi*g0) - delta(g0) vs g0
# ============================================================
print("-"*68)
print("--- Czesc 6: delta(phi*g0) - delta(g0) vs g0 ---")
print()

g0_6 = np.linspace(1.10, 1.45, 30)
print(f"  {'g0_e':>8}  {'delta_e':>10}  {'g0_mu=phi*g0_e':>16}  {'delta_mu':>10}  {'Ddelta':>10}  {'Ddelta/pi':>10}")
print("  " + "-"*75)
for g0e in g0_6:
    g0m = PHI * g0e
    de  = delta_scatt(g0e)
    dm  = delta_scatt(g0m)
    if np.isnan(de) or np.isnan(dm): continue
    Dd  = dm - de
    print(f"  {g0e:8.5f}  {de:10.6f}  {g0m:16.5f}  {dm:10.6f}  {Dd:10.6f}  {Dd/np.pi:10.6f}")
print()

# ============================================================
# CZESC 7: Czy istnieje g0 gdzie delta(phi*g0) - delta(g0) = pi?
# ============================================================
print("-"*68)
print("--- Czesc 7: delta(phi*g0) - delta(g0) = pi --- (oczekiwane g0 ~?) ---")
print()

def dDelta(g0):
    de = delta_scatt(g0)
    dm = delta_scatt(PHI * g0)
    if np.isnan(de) or np.isnan(dm): return np.nan
    return dm - de - np.pi

g0_7 = np.linspace(1.10, 1.45, 40)
dD_7 = [dDelta(g0) for g0 in g0_7]

pi_zeros = []
for i in range(len(dD_7)-1):
    if np.isnan(dD_7[i]) or np.isnan(dD_7[i+1]): continue
    if dD_7[i] * dD_7[i+1] < 0:
        try:
            g0z = brentq(dDelta, g0_7[i], g0_7[i+1], xtol=1e-6)
            de_z  = delta_scatt(g0z)
            dm_z  = delta_scatt(PHI*g0z)
            Ai_z  = A_inf(g0z)
            Aim_z = A_inf(PHI*g0z)
            r_z   = (Aim_z/Ai_z)**4 if Ai_z > 1e-6 else np.nan
            pi_zeros.append((g0z, de_z, dm_z, r_z))
            print(f"  ZERO: g0 = {g0z:.7f}  delta_e={de_z:.4f}  delta_mu={dm_z:.4f}")
            print(f"        ratio^4 = {r_z:.3f}  (vs {R21_EXP})")
        except Exception as e:
            print(f"  brentq err: {e}")

if not pi_zeros:
    print("  Brak zerowania dDelta - pi w zakresie [1.10, 1.45]")
    print("  Wartosc dDelta w kluczowych punktach:")
    for g0 in [1.2301, G0_STAR, 1.24, 1.25]:
        dD = dDelta(g0)
        print(f"    g0={g0:.5f}: dDelta = {dD:.5f} rad = {dD/np.pi:.5f}*pi")
print()

# ============================================================
# TESTY
# ============================================================
print("="*68)
print("TESTY")
print("="*68)

# T1: |delta_scatt(z0)| < 0.05 rad
d_z0 = delta_scatt(1.2301)
T1   = abs(d_z0) < 0.05
print(f"  T1: |delta(z0)| < 0.05 rad:            {'PASS' if T1 else 'FAIL'}  (delta={d_z0:.6f})")

# T2: |delta_scatt(phi*z0) - pi| < 0.5 rad
d_phi_z0 = delta_scatt(PHI * 1.2301)
T2 = abs(d_phi_z0 - np.pi) < 0.5 or abs(abs(d_phi_z0) - np.pi) < 0.5
print(f"  T2: |delta(phi*z0) - pi| < 0.5 rad:    {'PASS' if T2 else 'FAIL'}  (delta={d_phi_z0:.6f}, |d-pi|={abs(abs(d_phi_z0)-np.pi):.4f})")

# T3: Sprawdz czy g0* ma sens fazowy — porownaj delta(g0*) z najblizszym n*pi/4
d_star_2 = delta_scatt(G0_STAR)
best_frac = min([abs(d_star_2 - n*np.pi/4) for n in range(-4, 12)])
T3 = best_frac < 0.1
print(f"  T3: delta(g0*) blisko n*pi/4 (< 0.1):  {'PASS' if T3 else 'FAIL'}  (delta={d_star_2:.6f}, best_frac={best_frac:.4f})")

# T4: ciaglosc delta_uw — max jump < 0.3 rad (po unwrap)
# Sprawdzamy segment n=0 (g0 < 1.61)
n0_mask = [i for i,g0 in enumerate(g0_1) if g0 <= 1.62]
if len(n0_mask) > 3:
    seg = delta_uw[n0_mask]
    seg_valid = seg[~np.isnan(seg)]
    jumps = np.abs(np.diff(seg_valid))
    max_jump = np.max(jumps) if len(jumps) > 0 else 0.0
    T4 = max_jump < 0.3
    print(f"  T4: max jump delta_uw (n=0) < 0.3:     {'PASS' if T4 else 'FAIL'}  (max_jump={max_jump:.5f})")
else:
    T4 = False
    print(f"  T4: max jump < 0.3: FAIL  (brak danych)")

# T5: (A(g0_mu_BS)/A(g0_e_BS))^4 in [150, 300] dla jakiejs pary B-S
T5 = False
if BS_points:
    for ne in [0]:
        for nmu in [1,2]:
            if ne in BS_points and nmu in BS_points:
                _, _, Ae_bs, _ = BS_points[ne]
                _, _, Am_bs, _ = BS_points[nmu]
                if Ae_bs > 1e-6:
                    r_bs = (Am_bs/Ae_bs)**4
                    if 150 <= r_bs <= 300:
                        T5 = True
                        print(f"  T5: A-ratio^4 in [150,300] dla BS ({ne},{nmu}): PASS  (r={r_bs:.1f})")
                        break
        if T5: break
if not T5:
    print(f"  T5: A-ratio^4 in [150,300] dla BS pary:FAIL")

n_pass = sum([T1,T2,T3,T4,T5])
print(f"\nWYNIK: {n_pass}/5 testow przeszlo")

print()
print("="*68)
print("WNIOSEK FIZYCZNY (EX63)")
print("="*68)
print()
print(f"  delta_scatt(g0*)     = {d_star:.6f} rad = {d_star/np.pi:.6f} pi")
print(f"  delta_scatt(phi*g0*) = {d_star_mu:.6f} rad = {d_star_mu/np.pi:.6f} pi")
print(f"  Delta_delta          = {d_star_mu - d_star:.6f} rad = {(d_star_mu-d_star)/np.pi:.6f} pi")
print()
if pi_zeros:
    g0_pi, de_pi, dm_pi, r_pi = pi_zeros[0]
    print(f"  Punkt gdzie Delta_delta = pi: g0 = {g0_pi:.6f}  ratio = {r_pi:.2f}")
print()
print("="*68)
