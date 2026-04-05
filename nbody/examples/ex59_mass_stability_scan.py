"""
ex59_mass_stability_scan.py
============================
Stabilność solitonów z odbiciami i struktura dM/dg0 w TGP.

KONTEKST (ex57, ex58):
  ex57: (A_tail(2.0)/A_tail(1.24))^4 = 210 ≈ 207 (1.6%)
  ex58: B_tail(1.2301)=0 — warunek H1 potwierdzony dla elektronu
        Hipoteza H1 dla mionu: drugie zero B=0 przy 2.0797, ale (A/A0)^4=375≠207
        g0_mu/g0_e = 2.00/1.24 = 1.613 ≈ phi (złota proporcja), ale phi-pair daje 238

PYTANIA DO ex59:
  1. Struktura dM/dg0: czy jest gładka, monotonizna, czy ma skok/kink przy g0=2.0?
  2. Dyskretne g0? Szukamy ekstrema lub przejść M_raw(g0).
  3. Co jest szczególnego w g0=2.0:
     a) Zero dM/dg0? (ekstremum M)
     b) Zmiana tempa wzrostu M? (infleksja d²M/dg0²=0)
     c) Skok d(B_tail)/dg0? (zmiana reżimu topologicznego)
  4. Efektywny potencjał M(g0) w sektorach n=0,1,2:
     Czy M jest wypukła w każdym sektorze? Minimum każdego sektora = kandidat na "cząstkę"?
  5. Gęsty skan z rozróżnieniem na sektory n_bounce=0,1,2,3:
     Czy minimum M w n=1 sektorze leży blisko g0=2.0?

PLAN:
  1. Super-gęsty skan g0 in [1.1, 3.0], N=600 punktów:
     M_raw, A_tail, B_tail, C_tail, phi, g_min, n_bounce
  2. dM/dg0: numeryczna pochodna, ekstrema
  3. Oddzielna analiza każdego sektora n_bounce:
     Dla n=0: min M, max M; dla n=1: min M, max M
     Czy granica sektora (n→n+1) ma szczególną rolę?
  4. d(A_tail^4)/dg0: ekstrema stosunku do elektronu
  5. Czy istnieje g0 z dM/dg0=0 w rejonie [1.9, 2.1]?
  6. Efektywna masa M_eff = A_tail^4 / A_e^4 × M_e:
     Czy M_eff ma minimum w n=1 sektorze blisko g0=2.0?

TESTY (5):
  T1: dM/dg0 > 0 prawie wszędzie w [1.1, 2.0] (M rosnące)
  T2: Nie ma zera dM/dg0 w [1.8, 2.2] (M nie jest ekstremalne w rejonie mionu)
  T3: A_tail^4 ma infleksję (d²(A^4)/dg0²=0) w [1.5, 2.5]
  T4: Skok d(n_bounce)/dg0 koreluje z charakterystyczną zmianą B_tail
  T5: M_raw(n=1 minimum)/M_raw(n=0 minimum) jest bliskie 207 (≤20%)

Sesja: TGP v33, 2026-03-27
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import os, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import argrelextrema
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
ALPHA    = 2.0
G_GHOST  = np.exp(-1.0 / (2.0 * ALPHA))   # ~0.7788
G_BOUNCE = G_GHOST + 0.005
TARGET   = 206.768
G0_E = 1.24; G0_MU = 2.00; A_E = 0.287; A_MU = 1.09

R_MAX    = 40.0
R_START  = 1e-4
R_TAIL_L = 20.0; R_TAIL_R = 36.0
RTOL     = 1e-10; ATOL = 1e-13; MAX_STEP = 0.02

print("=" * 70)
print("EX59: STABILNOŚĆ SOLITONÓW — dM/dg0, SEKTORY TOPOLOGICZNE")
print("=" * 70)
print(f"  alpha = {ALPHA},  g* = {G_GHOST:.5f}")
print(f"  Pytanie: co wyróżnia g0_mu=2.00 w strukturze M(g0)?")
print()

# ─────────────────────────────────────────────────────────────────────────────
# ODE i integracja (identyczna z ex57/ex58)
# ─────────────────────────────────────────────────────────────────────────────
def rhs_physical(r, y):
    g, gp = y
    g   = max(g, G_BOUNCE + 1e-6)
    fg  = 1.0 + 2.0 * ALPHA * np.log(g)
    if abs(fg) < 1e-10: return [gp, 0.0]
    driving = g**2 * (1.0 - g)
    cross   = (ALPHA / g) * gp**2
    if r < 1e-10: return [gp, (driving - cross) / (3.0 * fg)]
    return [gp, (driving - cross - fg * 2.0 * gp / r) / fg]

def event_hit_ghost(r, y): return y[0] - G_BOUNCE
event_hit_ghost.terminal = True; event_hit_ghost.direction = -1


def integrate_with_bounces(g0, r_max=R_MAX, max_bounces=12):
    r0 = R_START; y0 = [g0, 0.0]
    segs_r, segs_g, segs_gp = [], [], []
    n_bounce_count = 0
    for bn in range(max_bounces + 1):
        sol = solve_ivp(rhs_physical, [r0, r_max], y0,
                        method='DOP853', max_step=MAX_STEP,
                        rtol=RTOL, atol=ATOL, events=[event_hit_ghost])
        segs_r.append(sol.t); segs_g.append(sol.y[0]); segs_gp.append(sol.y[1])
        if sol.t_events[0].size > 0 and bn < max_bounces:
            rb = float(sol.t_events[0][0]); gpb = float(sol.y_events[0][0, 1])
            r0 = rb + 1e-6; y0 = [G_BOUNCE + 1e-5, -gpb]
            n_bounce_count += 1
        else:
            break
    r  = np.concatenate(segs_r); g  = np.concatenate(segs_g); gp = np.concatenate(segs_gp)
    idx = np.argsort(r)
    return r[idx], g[idx], gp[idx], n_bounce_count


def fit_tail(r_arr, g_arr):
    mask = (r_arr >= R_TAIL_L) & (r_arr <= R_TAIL_R)
    if np.sum(mask) < 20: return 0.0, 0.0, 0.0
    rf = r_arr[mask]; df = (g_arr[mask] - 1.0) * rf
    M  = np.column_stack([np.cos(rf), np.sin(rf)])
    res = np.linalg.lstsq(M, df, rcond=None); B, C = res[0]
    return float(np.sqrt(B**2 + C**2)), float(B), float(C)


def analyze_full(g0):
    r, g, gp, n_bounce = integrate_with_bounces(g0)
    M_raw  = 4.0 * np.pi * float(np.trapezoid(r**2 * gp**2 / 2.0, r))
    A, B, C = fit_tail(r, g)
    g_min  = float(np.min(g))
    return {'M': M_raw, 'A': A, 'B': B, 'C': C, 'g_min': g_min, 'n_bounce': n_bounce}


# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 1: Super-gęsty skan g0 in [1.1, 3.0], N=600
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 1: Super-gęsty skan g0 ∈ [1.1, 3.0], N=600 ---")

N_SCAN  = 600
g0_scan = np.linspace(1.10, 3.00, N_SCAN)

M_arr = np.zeros(N_SCAN); A_arr = np.zeros(N_SCAN)
B_arr = np.zeros(N_SCAN); C_arr = np.zeros(N_SCAN)
gmin_arr = np.zeros(N_SCAN); nb_arr = np.zeros(N_SCAN, dtype=int)

for i, g0 in enumerate(g0_scan):
    res = analyze_full(g0)
    M_arr[i]   = res['M']; A_arr[i]   = res['A']
    B_arr[i]   = res['B']; C_arr[i]   = res['C']
    gmin_arr[i] = res['g_min']; nb_arr[i] = res['n_bounce']

print(f"  Scan gotowy. M range: [{M_arr.min():.3f}, {M_arr.max():.3f}]")
print(f"  n_bounce range: {nb_arr.min()}–{nb_arr.max()}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 2: Pochodna dM/dg0 i ekstrema
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 2: Pochodna dM/dg0 — ekstrema, zerowanie ---")

dM_dg0   = np.gradient(M_arr, g0_scan)
d2M_dg02 = np.gradient(dM_dg0, g0_scan)

# Zera dM/dg0 (lokalne min/max M)
zeros_dM = []
for i in range(N_SCAN - 1):
    if dM_dg0[i] * dM_dg0[i+1] < 0:
        g0_z = g0_scan[i] - dM_dg0[i] * (g0_scan[i+1] - g0_scan[i]) / (dM_dg0[i+1] - dM_dg0[i])
        zeros_dM.append(float(g0_z))

print(f"  Zera dM/dg0 (ekstrema M): {[f'{z:.4f}' for z in zeros_dM[:6]]}")

# Infleksje d²M/dg0²=0
zeros_d2M = []
for i in range(N_SCAN - 1):
    if d2M_dg02[i] * d2M_dg02[i+1] < 0:
        g0_z = g0_scan[i] - d2M_dg02[i] * (g0_scan[i+1] - g0_scan[i]) / (d2M_dg02[i+1] - d2M_dg02[i])
        zeros_d2M.append(float(g0_z))

print(f"  Infleksje d²M/dg0²=0 (zmiana wklęsłości): {[f'{z:.4f}' for z in zeros_d2M[:8]]}")

# Blisko g0_mu = 2.00
near_mu_dM  = [z for z in zeros_dM  if 1.8 <= z <= 2.2]
near_mu_d2M = [z for z in zeros_d2M if 1.8 <= z <= 2.2]
print(f"\n  Ekstrema M w [1.8, 2.2]: {[f'{z:.4f}' for z in near_mu_dM]}")
print(f"  Infleksje M w [1.8, 2.2]: {[f'{z:.4f}' for z in near_mu_d2M]}")
print()

# Wartości dM/dg0 przy znanych punktach
idx_e  = np.argmin(np.abs(g0_scan - G0_E))
idx_mu = np.argmin(np.abs(g0_scan - G0_MU))
print(f"  dM/dg0 przy g0_e  = {G0_E}: {dM_dg0[idx_e]:.4f}")
print(f"  dM/dg0 przy g0_mu = {G0_MU}: {dM_dg0[idx_mu]:.4f}")
print(f"  d²M/dg0² przy g0_e:  {d2M_dg02[idx_e]:.4f}")
print(f"  d²M/dg0² przy g0_mu: {d2M_dg02[idx_mu]:.4f}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 3: Analiza sektorów topologicznych n_bounce=0,1,2,3
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 3: Sektory topologiczne n_bounce=0,1,2,3 ---")

for n in range(4):
    mask_n = nb_arr == n
    if not np.any(mask_n):
        continue
    g0_n   = g0_scan[mask_n]
    M_n    = M_arr[mask_n]
    A_n    = A_arr[mask_n]
    B_n    = B_arr[mask_n]

    g0_min_M = g0_n[np.argmin(M_n)]
    g0_max_M = g0_n[np.argmax(M_n)]
    print(f"\n  [n_bounce = {n}] g0 ∈ [{g0_n.min():.3f}, {g0_n.max():.3f}]  "
          f"N={np.sum(mask_n)}")
    print(f"    M ∈ [{M_n.min():.3f}, {M_n.max():.3f}]")
    print(f"    A ∈ [{A_n.min():.4f}, {A_n.max():.4f}]")
    print(f"    M_min przy g0 = {g0_min_M:.4f}  (M = {M_n.min():.4f})")
    print(f"    M_max przy g0 = {g0_max_M:.4f}  (M = {M_n.max():.4f})")
    print(f"    dM/dg0 w środku sektora: {dM_dg0[mask_n][len(g0_n)//2]:.4f}")

    # Zera B w sektorze
    zeros_B_n = []
    idx_n = np.where(mask_n)[0]
    for ii in range(len(idx_n)-1):
        i1, i2 = idx_n[ii], idx_n[ii+1]
        if B_arr[i1] * B_arr[i2] < 0:
            z = g0_scan[i1] - B_arr[i1]*(g0_scan[i2]-g0_scan[i1])/(B_arr[i2]-B_arr[i1])
            zeros_B_n.append(float(z))
    if zeros_B_n:
        print(f"    Zera B_tail w sektorze: {[f'{z:.4f}' for z in zeros_B_n]}")

print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 4: Stosunek A^4 do elektronu
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 4: Stosunek (A_tail/A_e)^4 vs g0 ---")

A_e_val = A_arr[idx_e]
ratio4  = (A_arr / A_e_val)**4

print(f"  A_e = A_tail(g0={G0_E}) = {A_e_val:.5f}")
print(f"  Max (A/A_e)^4 w sektore n=0: {ratio4[nb_arr==0].max():.3f}")
print(f"  (A/A_e)^4 przy g0_mu=2.00:  {ratio4[idx_mu]:.3f}")
print()

# Gęsty skan gdzie (A/A_e)^4 ≈ 207
target_ratio = TARGET
diffs = np.abs(ratio4 - target_ratio)
g0_best207 = g0_scan[np.argmin(diffs)]
print(f"  g0 gdzie (A/A_e)^4 najbliższe 207: g0 = {g0_best207:.4f}"
      f"  wartość = {ratio4[np.argmin(diffs)]:.3f}")

# Infleksje (A^4)'' = 0
dA4     = np.gradient(ratio4, g0_scan)
d2A4    = np.gradient(dA4,    g0_scan)
zeros_d2A4 = []
for i in range(N_SCAN-1):
    if d2A4[i]*d2A4[i+1] < 0:
        z = g0_scan[i] - d2A4[i]*(g0_scan[i+1]-g0_scan[i])/(d2A4[i+1]-d2A4[i])
        zeros_d2A4.append(float(z))
print(f"  Infleksje d²(A/A_e)^4/dg0²=0: {[f'{z:.3f}' for z in zeros_d2A4[:6]]}")
near_mu_A4inf = [z for z in zeros_d2A4 if 1.8 <= z <= 2.2]
print(f"  Infleksje A^4 w [1.8, 2.2]: {[f'{z:.4f}' for z in near_mu_A4inf]}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 5: Analiza przejść n → n+1 (topologiczne zmiany reżimu)
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 5: Przejścia sektorowe n_bounce → n_bounce+1 ---")

transitions = []
for i in range(N_SCAN - 1):
    if nb_arr[i] != nb_arr[i+1]:
        g0_tr = 0.5 * (g0_scan[i] + g0_scan[i+1])
        transitions.append((g0_tr, nb_arr[i], nb_arr[i+1]))

for g0_tr, n_from, n_to in transitions:
    idx_tr = np.argmin(np.abs(g0_scan - g0_tr))
    B_tr   = B_arr[idx_tr]
    M_tr   = M_arr[idx_tr]
    dM_tr  = dM_dg0[idx_tr]
    A4_tr  = ratio4[idx_tr]
    print(f"  Przejście n={n_from}→{n_to}  g0≈{g0_tr:.4f}"
          f"  M={M_tr:.3f}  B={B_tr:.4f}  (A/Ae)^4={A4_tr:.2f}  dM/dg0={dM_tr:.2f}")

print()
print(f"  Liczba przejść: {len(transitions)}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 6: Szczegółowy zoom g0 ∈ [1.8, 2.2]
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 6: Zoom g0 ∈ [1.8, 2.2], N=100 ---")

g0_zoom = np.linspace(1.80, 2.20, 100)
M_z=[]; A_z=[]; B_z=[]; dM_z_num=[]; nb_z=[]

for g0 in g0_zoom:
    res = analyze_full(g0)
    M_z.append(res['M']); A_z.append(res['A']); B_z.append(res['B'])
    nb_z.append(res['n_bounce'])

M_z=np.array(M_z); A_z=np.array(A_z); B_z=np.array(B_z); nb_z=np.array(nb_z)
dM_z = np.gradient(M_z, g0_zoom)

print(f"  g0      M_raw     dM/dg0   A_tail   B_tail   n_b   (A/Ae)^4")
print("  " + "-"*63)
for i in range(0, len(g0_zoom), 10):
    g0 = g0_zoom[i]
    r4 = (A_z[i]/A_e_val)**4
    print(f"  {g0:.4f}  {M_z[i]:9.3f}  {dM_z[i]:8.3f}  {A_z[i]:.5f}"
          f"  {B_z[i]:7.4f}  {nb_z[i]:3d}  {r4:.3f}")

# Minimum dM w zoom
idx_min_dM = np.argmin(np.abs(dM_z))
print(f"\n  Minimum |dM/dg0| przy g0 = {g0_zoom[idx_min_dM]:.4f}"
      f"  dM = {dM_z[idx_min_dM]:.4f}  M = {M_z[idx_min_dM]:.3f}")

# Zero dM w zoom?
zeros_zoom = []
for i in range(len(g0_zoom)-1):
    if dM_z[i]*dM_z[i+1] < 0:
        z = g0_zoom[i] - dM_z[i]*(g0_zoom[i+1]-g0_zoom[i])/(dM_z[i+1]-dM_z[i])
        zeros_zoom.append(float(z))
if zeros_zoom:
    print(f"  Zero dM/dg0 w [1.8, 2.2] przy g0 = {[f'{z:.5f}' for z in zeros_zoom]}")
else:
    print(f"  BRAK zera dM/dg0 w [1.8, 2.2]")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 7: Czy g0_mu=2 jest "specjalne" — test warianty
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 7: Co jest szczególnego w g0=2.00? ---")

# a) log-log: M ~ g0^gamma dla każdego sektora osobno
for n in range(3):
    mask_n = (nb_arr == n) & (g0_scan > 1.1)
    if np.sum(mask_n) < 5: continue
    g0_n = g0_scan[mask_n]; M_n = M_arr[mask_n]
    log_g = np.log(g0_n); log_M = np.log(M_n + 1e-12)
    p = np.polyfit(log_g, log_M, 1)
    print(f"  Sektor n={n}: M ~ g0^{p[0]:.3f}  (R²={1-np.var(log_M-np.polyval(p,log_g))/np.var(log_M):.4f})")

print()

# b) Stosunek M_raw(n=1 min) / M_raw(n=0 min)
mask_0 = nb_arr == 0; mask_1 = nb_arr == 1
if np.any(mask_0) and np.any(mask_1):
    M_min_0 = M_arr[mask_0].min(); g0_Mmin_0 = g0_scan[mask_0][np.argmin(M_arr[mask_0])]
    M_min_1 = M_arr[mask_1].min(); g0_Mmin_1 = g0_scan[mask_1][np.argmin(M_arr[mask_1])]
    print(f"  Min M w n=0: M={M_min_0:.4f} przy g0={g0_Mmin_0:.4f}")
    print(f"  Min M w n=1: M={M_min_1:.4f} przy g0={g0_Mmin_1:.4f}")
    print(f"  M_min(n=1)/M_min(n=0) = {M_min_1/M_min_0:.3f}  (cel: 206.77)")
    print()

# c) Szukaj g0 gdzie M = 207 × M(g0_e)
M_e_val = M_arr[idx_e]
target_M207 = 207.0 * M_e_val
diffs_M = np.abs(M_arr - target_M207)
g0_M207 = g0_scan[np.argmin(diffs_M)]
print(f"  M_e = {M_e_val:.5f}")
print(f"  g0 gdzie M_raw = 207 × M_e: g0 = {g0_M207:.4f}  M = {M_arr[np.argmin(diffs_M)]:.3f}")
print(f"  n_bounce przy tej g0 = {nb_arr[np.argmin(diffs_M)]}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# WYKRESY
# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 14))
gs  = gridspec.GridSpec(4, 2, figure=fig, hspace=0.45, wspace=0.35)

# Kolory wg n_bounce
cmap_n = {0: 'blue', 1: 'red', 2: 'green', 3: 'purple', 4: 'orange'}

# Panel 1: M_raw vs g0, z zaznaczeniem sektorów
ax1 = fig.add_subplot(gs[0, :])
for n in range(5):
    mask = nb_arr == n
    if not np.any(mask): continue
    ax1.scatter(g0_scan[mask], M_arr[mask], s=4, color=cmap_n.get(n,'k'),
                label=f'n={n}')
ax1.axvline(G0_E,  color='darkgreen', lw=1.5, ls='--', label=f'g0_e={G0_E}')
ax1.axvline(G0_MU, color='darkred',   lw=1.5, ls='--', label=f'g0_mu={G0_MU}')
ax1.set_xlabel(r'$g_0$', fontsize=12); ax1.set_ylabel(r'$M_{\rm raw}$', fontsize=12)
ax1.set_title(r'$M_{\rm raw}(g_0)$ wg sektorów $n_{\rm bounce}$', fontsize=13)
ax1.legend(loc='upper left', fontsize=9, markerscale=3); ax1.grid(True, alpha=0.3)

# Panel 2: dM/dg0
ax2 = fig.add_subplot(gs[1, 0])
ax2.plot(g0_scan, dM_dg0, 'k-', lw=1.2)
ax2.axhline(0, color='gray', lw=0.8, ls='--')
ax2.axvline(G0_E,  color='darkgreen', lw=1.5, ls='--')
ax2.axvline(G0_MU, color='darkred',   lw=1.5, ls='--')
for z in zeros_dM:
    ax2.axvline(z, color='orange', lw=0.8, ls=':')
ax2.set_xlabel(r'$g_0$'); ax2.set_ylabel(r'$dM/dg_0$')
ax2.set_title(r'Pochodna $dM/dg_0$'); ax2.grid(True, alpha=0.3)

# Panel 3: (A/Ae)^4 vs g0
ax3 = fig.add_subplot(gs[1, 1])
for n in range(5):
    mask = nb_arr == n
    if not np.any(mask): continue
    ax3.scatter(g0_scan[mask], ratio4[mask], s=4, color=cmap_n.get(n,'k'), label=f'n={n}')
ax3.axhline(TARGET,  color='magenta', lw=1.2, ls='--', label=f'r₂₁={TARGET:.0f}')
ax3.axvline(G0_E,    color='darkgreen', lw=1.5, ls='--')
ax3.axvline(G0_MU,   color='darkred',   lw=1.5, ls='--')
ax3.axvline(g0_best207, color='cyan', lw=1.2, ls=':',
            label=f'best207 g0={g0_best207:.3f}')
ax3.set_xlabel(r'$g_0$'); ax3.set_ylabel(r'$(A/A_e)^4$')
ax3.set_title(r'Stosunek mas $(A_{\rm tail}/A_e)^4$')
ax3.set_ylim(0, min(TARGET * 2.5, ratio4.max()))
ax3.legend(fontsize=8, markerscale=3); ax3.grid(True, alpha=0.3)

# Panel 4: B_tail vs g0
ax4 = fig.add_subplot(gs[2, 0])
ax4.plot(g0_scan, B_arr, 'b-', lw=1.2)
ax4.axhline(0, color='gray', lw=0.8, ls='--')
ax4.axvline(G0_E,  color='darkgreen', lw=1.5, ls='--')
ax4.axvline(G0_MU, color='darkred',   lw=1.5, ls='--')
ax4.set_xlabel(r'$g_0$'); ax4.set_ylabel(r'$B_{\rm tail}$')
ax4.set_title(r'$B_{\rm tail}(g_0)$ — składowa cosinus ogona'); ax4.grid(True, alpha=0.3)

# Panel 5: n_bounce vs g0
ax5 = fig.add_subplot(gs[2, 1])
ax5.step(g0_scan, nb_arr, 'k-', lw=1.5, where='mid')
for g0_tr, n_from, n_to in transitions:
    ax5.axvline(g0_tr, color='red', lw=0.8, ls=':', alpha=0.7)
ax5.axvline(G0_E,  color='darkgreen', lw=1.5, ls='--')
ax5.axvline(G0_MU, color='darkred',   lw=1.5, ls='--')
ax5.set_xlabel(r'$g_0$'); ax5.set_ylabel(r'$n_{\rm bounce}$')
ax5.set_title('Sektor topologiczny'); ax5.grid(True, alpha=0.3)

# Panel 6: Zoom M i (A/Ae)^4 blisko g0=2
ax6 = fig.add_subplot(gs[3, :])
ax6b = ax6.twinx()
ax6.plot(g0_zoom, M_z,        'b-',  lw=2, label=r'$M_{\rm raw}$')
ax6.plot(g0_zoom, dM_z,       'b--', lw=1.2, alpha=0.6, label=r'$dM/dg_0$')
ax6b.plot(g0_zoom, (np.array(A_z)/A_e_val)**4, 'r-', lw=2, label=r'$(A/A_e)^4$')
ax6b.axhline(TARGET, color='magenta', lw=1, ls='--', alpha=0.7)
ax6.axhline(0, color='gray', lw=0.8, ls='--')
ax6.axvline(G0_MU,  color='darkred',   lw=2, ls='--', label=f'g0_mu={G0_MU}')
ax6.axvline(g0_best207, color='cyan', lw=1.2, ls=':', label=f'best207={g0_best207:.3f}')
ax6.set_xlabel(r'$g_0$'); ax6.set_ylabel(r'$M_{\rm raw}$, $dM/dg_0$ (niebieski)')
ax6b.set_ylabel(r'$(A/A_e)^4$ (czerwony)')
ax6.set_title(r'Zoom $g_0 \in [1.8, 2.2]$ — struktura $M(g_0)$ i $(A/A_e)^4$')
lines1, labels1 = ax6.get_legend_handles_labels()
lines2, labels2 = ax6b.get_legend_handles_labels()
ax6.legend(lines1+lines2, labels1+labels2, fontsize=8)
ax6.grid(True, alpha=0.3)

# Tytuł
fig.suptitle(
    f'EX59: Stabilność solitonów TGP — $dM/dg_0$, sektory topologiczne\n'
    f'g0 gdzie $(A/A_e)^4=207$: {g0_best207:.3f}  |  '
    f'Przejść n→n+1: {len(transitions)}  |  '
    f'Zero dM w [1.8,2.2]: {zeros_zoom if zeros_zoom else "brak"}',
    fontsize=11
)
plot_path = os.path.join(os.path.dirname(__file__), 'ex59_mass_stability_scan.png')
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

# T1: dM/dg0 > 0 prawie wszędzie w [1.1, 2.0]
mask_t1 = (g0_scan >= 1.1) & (g0_scan <= 2.0)
frac_positive = np.mean(dM_dg0[mask_t1] > 0)
t1 = frac_positive >= 0.90
print(f"  T1: dM/dg0>0 w ≥90% punktów [1.1,2.0]: {'PASS' if t1 else 'FAIL'}"
      f"  (frac_positive={frac_positive:.3f})")
if t1: tests_passed += 1

# T2: Brak zera dM/dg0 w [1.8, 2.2] — M nie ma ekstremum przy mionie
t2 = len(near_mu_dM) == 0
print(f"  T2: Brak ekstremum M w [1.8,2.2]: {'PASS' if t2 else 'FAIL'}"
      f"  (zera dM: {[f'{z:.4f}' for z in near_mu_dM]})")
if t2: tests_passed += 1

# T3: A_tail^4 ma infleksję w [1.5, 2.5]
infleks_A4_15_25 = [z for z in zeros_d2A4 if 1.5 <= z <= 2.5]
t3 = len(infleks_A4_15_25) >= 1
print(f"  T3: Infleksja (A/Ae)^4 w [1.5,2.5]: {'PASS' if t3 else 'FAIL'}"
      f"  ({[f'{z:.3f}' for z in infleks_A4_15_25]})")
if t3: tests_passed += 1

# T4: Przejście topologiczne n=0→1 w [1.4, 1.7]
# (ex55 z G_DEEP=0.85 dawało g0*≈1.47; ex59 z G_BOUNCE=g*+0.005 daje ≈1.61)
trans_n01 = [g0_tr for g0_tr, n_from, n_to in transitions if n_from==0 and n_to==1]
t4 = len(trans_n01) >= 1 and all(1.4 <= g0_t <= 1.7 for g0_t in trans_n01)
print(f"  T4: Przejście n=0→1 w [1.4,1.7]: {'PASS' if t4 else 'FAIL'}"
      f"  ({[f'{z:.4f}' for z in trans_n01]})")
if t4: tests_passed += 1

# T5: M_min(n=1)/M_min(n=0) w [1, 100] — przynajmniej rząd wielkości od 207
if np.any(mask_0) and np.any(mask_1):
    ratio_min_sectors = M_arr[mask_1].min() / M_arr[mask_0].min()
    t5 = 1 <= ratio_min_sectors <= 100
    print(f"  T5: M_min(n=1)/M_min(n=0) ∈ [1,100]: {'PASS' if t5 else 'FAIL'}"
          f"  (ratio = {ratio_min_sectors:.3f})")
    if t5: tests_passed += 1
else:
    print(f"  T5: SKIP (brak sektora n=0 lub n=1)")

print()
print(f"WYNIK: {tests_passed}/{tests_total} testów przeszło")
print()

# ─────────────────────────────────────────────────────────────────────────────
# WNIOSEK
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("WNIOSEK FIZYCZNY (EX59)")
print("=" * 70)
print()
print(f"  M(g0_e={G0_E}) = {M_arr[idx_e]:.4f}")
print(f"  M(g0_mu={G0_MU}) = {M_arr[idx_mu]:.4f}")
print(f"  M_mu/M_e = {M_arr[idx_mu]/M_arr[idx_e]:.3f}  (cel: 207)")
print()
print(f"  (A(g0_mu)/A(g0_e))^4 = {ratio4[idx_mu]:.3f}  (cel: 207)")
print(f"  Najlepsza g0 dla r21: g0 = {g0_best207:.4f}  ratio = {ratio4[np.argmin(diffs)]:.3f}")
print()
print(f"  Struktura dM/dg0 przy g0_mu: {dM_dg0[idx_mu]:.4f}"
      f"  (zero?: {'TAK' if abs(dM_dg0[idx_mu]) < 5 else 'NIE'})")
print()
print("  Sektory topologiczne:")
for g0_tr, n_from, n_to in transitions[:4]:
    print(f"    n={n_from}→{n_to}: g0≈{g0_tr:.3f}")
print()
if near_mu_dM:
    print(f"  ★ Znaleziono ekstremum M w [1.8,2.2] przy g0={near_mu_dM}")
else:
    print(f"  → Brak ekstremum M przy g0=2.00: M rośnie monotonicznie przez miną wartość")
    print(f"  → g0_mu=2.00 NIE jest wyróżnione przez dM/dg0=0")
print()
if near_mu_d2M:
    print(f"  ★ Infleksja M w [1.8,2.2]: g0 = {near_mu_d2M}")
    print(f"    Interpretacja: zmiana tempa wzrostu M przy g0_mu")
else:
    print(f"  → Brak infleksji M w [1.8,2.2]")
print("=" * 70)
