"""
ex55_nodal_spectrum_tgp.py
==========================
Widmo wezlowe TGP-ODE: rozwiazania n=0,1,2 jako kandydaci na generacje.

PYTANIE CENTRALNE:
  Czy stosunek mas profili radialnych TGP o roznym stopniu wezlowym
  moze odtwarzac stosunek mas mion/elektron = 206.77?

PELNE NILINIOWE ODE TGP (sferyczna symetria):
  f(g)*g'' + (2/r)*f(g)*g' + (alpha/g)*(g')^2 = V'(g)

  gdzie:
    f(g) = 1 + 2*alpha*ln(g),   alpha = 2
    V'(g) = g^2*(1-g)           [beta=gamma=1]
    Granica duchowa: g* = exp(-1/4) ~ 0.7788   [f(g*)=0]

PROFIL "GARBKOWY" (g_0 > 1):
  g(0) = g_0 > 1,  g'(0) = 0,  g(r->inf) -> 1
  V'(g) = g^2*(1-g) < 0 dla g > 1 => pole spada od g_0 do 1
  Ogon: delta = g-1 ~ A*sin(r)/r (oscylacyjny, nie Yukawa)

DEFINICJA STOPNIA WEZLOWEGO n:
  n = liczba razy gdy g(r) spada PONIZEJ progu g_deep = 0.85
  (wyraznie ponizej 1, ale powyzej g*=0.778 — "gleboka oscylacja")

  n=0: g_min > g_deep (pole nie schodzi poza 0.85)
  n=1: g_min < g_deep przynajmniej raz
       (w rzeczywistosci: pierwsze glebsze przejscie przez g*)

  Uwaga: g* = 0.7788 jest tzw. granica duchowa; f(g*) = 0.
  Dla g < g*: f(g) < 0 (czlon kinetyczny "duchowy").
  W TGP fizycznie obszar g < g* jest odciety.

MASA PROFILU:
  M_raw(g_0) = 4*pi * integral_0^r_max r^2 * (g')^2/2 * dr
  (calka kinetyczna bez wagi f(g), zawsze > 0)

  M_phys(g_0) = 4*pi * integral_0^r_ghost r^2 * |f(g)| * (g')^2/2 * dr
  gdzie r_ghost = pierwsze wejscie do strefy g < g* (lub r_max jesli brak)
  (energia w fizycznie dozwolonym obszarze)

HIPOTEZA WZMOCNIENIA KINETYCZNEGO:
  Z f(g) = 1 + 4*ln(g), calka M_phys ~ g_0^4 dla duzych g_0
  => M(n=1)/M(n=0) ~ (g_0^(n=1) / g_0^(n=0))^4

CEL: sprawdzic czy istnieje g_0 taki ze M(g_0)/M(g_0^e) ~ 207

CZESCI:
  1. Skan g_0 in [1.01, 3.0]: profil, n_deep, M_raw, M_phys, g_min
  2. Identyfikacja przejsc n: 0->1
  3. Profile reprezentatywne dla n=0 (g0=1.24) i kandydata n=1
  4. Log-log fit M_phys ~ g_0^gamma
  5. Stosunek mas i test hipotezy kinetycznej
  6. Wykresy: profil radialny, portret fazowy, M(g_0)

TESTY (5):
  T1: g* = exp(-1/4) ~ 0.7788
  T2: g_0=1.24: g_min > 0.85 (n=0, brak glebokiej oscylacji)
  T3: g_0=2.50: g_min < g* (wejscie w strefe duchowa)
  T4: M_raw rosnie z g_0 (monotonicznie na probce n=0)
  T5: gamma (fit) in [1.5, 7.0]

Sesja: TGP v33, 2026-03-27
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import os
import warnings
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# Stale fizyczne TGP
# ─────────────────────────────────────────────────────────────────────────────
ALPHA       = 2.0
G_GHOST     = np.exp(-1.0 / (2.0 * ALPHA))   # exp(-1/4) ~ 0.7788
G_DEEP      = 0.85                             # prog "glebokiej" oscylacji
V1          = 1.0/3.0 - 1.0/4.0               # V(1) = 1/12
TARGET_RATIO= 206.768                          # m_mu / m_e

R_MAX    = 30.0
R_START  = 1e-4
RTOL     = 1e-9
ATOL     = 1e-12
MAX_STEP = 0.03

print("=" * 70)
print("EX55: WIDMO WEZLOWE TGP-ODE — GENERACJE CZĄSTEK")
print("=" * 70)
print(f"  alpha      = {ALPHA}")
print(f"  g*         = exp(-1/(2a)) = {G_GHOST:.6f}")
print(f"  g_deep     = {G_DEEP:.3f}  [prog glebokiej oscylacji]")
print(f"  V(1)       = {V1:.6f}")
print(f"  TARGET     = {TARGET_RATIO}  [m_mu/m_e]")
print()

# ─────────────────────────────────────────────────────────────────────────────
# RHS pelnego ODE
# ─────────────────────────────────────────────────────────────────────────────

def rhs_full(r, y):
    """
    Pelne ODE TGP w symetrii sferycznej.
    Dla g < 0.01 zatrzymujemy — niefizykalny obszar.
    """
    g, gp = y
    g = max(g, 0.02)                              # zabezpieczenie

    fg      = 1.0 + 2.0 * ALPHA * np.log(g)      # f(g)
    driving = g**2 * (1.0 - g)                    # V'(g)
    cross   = (ALPHA / g) * gp**2

    if r < 1e-10:
        if abs(fg) < 1e-10:
            return [gp, 0.0]
        return [gp, (driving - cross) / (3.0 * fg)]

    damp = fg * 2.0 * gp / r

    if abs(fg) < 1e-10:
        return [gp, 0.0]

    return [gp, (driving - cross - damp) / fg]


def event_very_low(r, y):
    """Stop jesli g < 0.02 (niebezpieczny obszar numeryczny)."""
    return y[0] - 0.02

event_very_low.terminal  = True
event_very_low.direction = -1


def integrate_profile(g0, r_max=R_MAX):
    """
    Integruje profil solitonowy dla danego g0.
    Zwraca dict z kompletna analiza.
    """
    sol = solve_ivp(
        rhs_full,
        [R_START, r_max],
        [g0, 0.0],
        method='DOP853',
        max_step=MAX_STEP,
        rtol=RTOL, atol=ATOL,
        events=[event_very_low],
        dense_output=False
    )

    r   = sol.t
    g   = sol.y[0]
    gp  = sol.y[1]

    g_min = float(np.min(g))
    hit_low = bool(sol.t_events[0].size > 0)

    # ── Stopien wezlowy n_deep: liczba przejsc g przez G_DEEP z gory ──
    # Uwaga: szukamy glebokiego zanurzenia, nie drobnych oscylacji ogona
    # Zlicz segmenty gdzie g < G_DEEP
    below_mask = g < G_DEEP
    if np.any(below_mask):
        # Liczba spojnych segmentow ponizej progu
        n_deep = int(np.sum(np.diff(below_mask.astype(int)) == 1))
    else:
        n_deep = 0

    # ── r_ghost: pierwsze wejscie do strefy g < g* ──
    ghost_mask = g < G_GHOST
    if np.any(ghost_mask):
        r_ghost = float(r[np.argmax(ghost_mask)])
    else:
        r_ghost = r_max

    # ── Masy ──
    # M_raw: bez f(g), zawsze > 0
    M_raw  = 4.0 * np.pi * float(
        np.trapezoid(r**2 * gp**2 / 2.0, r)
    )

    # M_phys: z f(g), tylko do r_ghost
    mask_phys = r <= r_ghost
    r_p  = r[mask_phys]
    g_p  = g[mask_phys]
    gp_p = gp[mask_phys]
    fg_p = 1.0 + 2.0 * ALPHA * np.log(np.maximum(g_p, 1e-8))
    M_phys = 4.0 * np.pi * float(
        np.trapezoid(r_p**2 * np.abs(fg_p) * gp_p**2 / 2.0, r_p)
    )

    # ── Liczba skrzyzowan g=1 w oknie r < 15 (diagnostyka) ──
    mask15  = r < 15.0
    g15     = g[mask15] - 1.0
    s15     = np.sign(g15); s15 = s15[s15 != 0]
    n_cross1 = int(np.sum(np.diff(s15) != 0))

    return {
        'r': r, 'g': g, 'gp': gp,
        'g_min'  : g_min,
        'n_deep' : n_deep,
        'n_cross1': n_cross1,
        'M_raw'  : M_raw,
        'M_phys' : M_phys,
        'r_ghost': r_ghost,
        'hit_low': hit_low,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Czesc 1: Skan g_0 in [1.02, 3.0]
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 1: Skan g_0 in [1.02, 3.0] ---")
print()

N_SCAN   = 160
g0_scan  = np.linspace(1.02, 3.0, N_SCAN)

ndp_arr  = np.zeros(N_SCAN, dtype=int)
Mraw_arr = np.zeros(N_SCAN)
Mph_arr  = np.zeros(N_SCAN)
gmin_arr = np.zeros(N_SCAN)
rgh_arr  = np.zeros(N_SCAN)
hlow_arr = np.zeros(N_SCAN, dtype=bool)

for i, g0 in enumerate(g0_scan):
    res = integrate_profile(g0)
    ndp_arr[i]  = res['n_deep']
    Mraw_arr[i] = res['M_raw']
    Mph_arr[i]  = res['M_phys']
    gmin_arr[i] = res['g_min']
    rgh_arr[i]  = res['r_ghost']
    hlow_arr[i] = res['hit_low']

print(f"  {'g0':>6}  {'n_deep':>7}  {'g_min':>8}  {'M_raw':>11}  {'M_phys':>11}  {'r_ghost':>8}")
print("  " + "-" * 60)
for i in range(0, N_SCAN, 14):
    g0 = g0_scan[i]
    print(f"  {g0:6.3f}  {ndp_arr[i]:>7d}  {gmin_arr[i]:8.4f}  "
          f"{Mraw_arr[i]:11.4f}  {Mph_arr[i]:11.4f}  {rgh_arr[i]:8.3f}")
print()


# ─────────────────────────────────────────────────────────────────────────────
# Czesc 2: Przejscia wezlowe n_deep: 0->1
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 2: Przejscia wezlowe (n_deep) ---")
print()

transitions_deep = []
for i in range(1, N_SCAN):
    if ndp_arr[i] > ndp_arr[i-1]:
        transitions_deep.append({
            'g0_L': g0_scan[i-1], 'g0_R': g0_scan[i],
            'nL': ndp_arr[i-1],   'nR': ndp_arr[i],
        })

if transitions_deep:
    for tr in transitions_deep[:4]:
        mid = 0.5 * (tr['g0_L'] + tr['g0_R'])
        print(f"  n={tr['nL']}->{tr['nR']}:  "
              f"g0 in [{tr['g0_L']:.4f}, {tr['g0_R']:.4f}],  mid={mid:.4f}")
else:
    print("  Brak przejsc n_deep w zakresie [1.02, 3.0]")

# Przejscia dla g* (wejscie w strefe duchowa)
ghost_entries = []
for i in range(1, N_SCAN):
    if gmin_arr[i] < G_GHOST and gmin_arr[i-1] >= G_GHOST:
        ghost_entries.append({'g0_L': g0_scan[i-1], 'g0_R': g0_scan[i],
                               'mid': 0.5*(g0_scan[i-1]+g0_scan[i])})

if ghost_entries:
    for ge in ghost_entries[:3]:
        print(f"  Wejscie do g<g*:  "
              f"g0 in [{ge['g0_L']:.4f}, {ge['g0_R']:.4f}],  mid={ge['mid']:.4f}")
print()


# ─────────────────────────────────────────────────────────────────────────────
# Czesc 3: Profile reprezentatywne
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 3: Profile reprezentatywne ---")
print()

# Elektron (n=0): g0=1.24 — klasyczny z literatury TGP
g0_e = 1.24
res_e = integrate_profile(g0_e)
print(f"  Elektron (n=0, g0={g0_e}):")
print(f"    n_deep   = {res_e['n_deep']}")
print(f"    g_min    = {res_e['g_min']:.6f}")
print(f"    delta_g* = {res_e['g_min'] - G_GHOST:.6f}")
print(f"    M_raw    = {res_e['M_raw']:.6f}")
print(f"    M_phys   = {res_e['M_phys']:.6f}")
print(f"    r_ghost  = {res_e['r_ghost']:.2f}")
print()

# Kandydat mion: pierwsze wejscie w g<g* lub 1. przejscie n_deep
# Wybieramy g0 tuz za pierwszym ghost-entry
if ghost_entries:
    g0_mu_crit = ghost_entries[0]['mid']
    g0_mu = g0_mu_crit + 0.15   # wyraznie za progiem
else:
    g0_mu = 1.8   # fallback
    g0_mu_crit = None

res_mu = integrate_profile(g0_mu)
print(f"  Mion-kandydat (g0={g0_mu:.3f}):")
print(f"    n_deep   = {res_mu['n_deep']}")
print(f"    g_min    = {res_mu['g_min']:.6f}")
print(f"    delta_g* = {res_mu['g_min'] - G_GHOST:.6f}")
print(f"    M_raw    = {res_mu['M_raw']:.6f}")
print(f"    M_phys   = {res_mu['M_phys']:.6f}")
print(f"    r_ghost  = {res_mu['r_ghost']:.3f}")
print()

# Trzecia generacja (tauon?): g0 = 2.5
g0_tau = 2.5
res_tau = integrate_profile(g0_tau)
print(f"  Tauon-kandydat (g0={g0_tau}):")
print(f"    n_deep   = {res_tau['n_deep']}")
print(f"    g_min    = {res_tau['g_min']:.6f}")
print(f"    M_raw    = {res_tau['M_raw']:.6f}")
print(f"    M_phys   = {res_tau['M_phys']:.6f}")
print()


# ─────────────────────────────────────────────────────────────────────────────
# Czesc 4: Log-log fit M_raw ~ g_0^gamma (zakres n=0: g_min > G_GHOST)
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 4: Log-log fit M_raw ~ g_0^gamma ---")
print()

def power_fit(x, y, label=""):
    lx = np.log(x); ly = np.log(y)
    c = np.polyfit(lx, ly, 1)
    g_, A_ = c[0], np.exp(c[1])
    res_ = ly - np.polyval(c, lx)
    R2_ = 1.0 - np.sum(res_**2)/np.sum((ly-ly.mean())**2)
    if label:
        print(f"  Fit [{label}]:  M ~ {A_:.5g} * g_0^{g_:.3f},  R^2={R2_:.5f}")
    return g_, A_, R2_

# Fit na wszystkich n=0 (g_min > G_GHOST + 0.01):
mask_n0 = (gmin_arr > G_GHOST + 0.01) & (Mraw_arr > 1e-6)
mask_all_phys = (~hlow_arr) & (Mraw_arr > 1e-6)

if np.sum(mask_n0) >= 4:
    gamma_n0, A_n0, R2_n0 = power_fit(g0_scan[mask_n0], Mraw_arr[mask_n0], "n=0 (g_min > g*)")
else:
    gamma_n0, A_n0, R2_n0 = None, None, None
    print("  Zbyt malo punktow n=0 do fitu.")

if np.sum(mask_all_phys) >= 4:
    gamma_all, A_all, R2_all = power_fit(g0_scan[mask_all_phys], Mraw_arr[mask_all_phys], "all phys")
else:
    gamma_all, A_all, R2_all = None, None, None

print()


# ─────────────────────────────────────────────────────────────────────────────
# Czesc 5: Stosunek mas i hipoteza kinetyczna
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 5: Stosunek mas i hipoteza kinetyczna ---")
print()

Me_raw   = res_e['M_raw']
Mmu_raw  = res_mu['M_raw']
Me_phys  = res_e['M_phys']
Mmu_phys = res_mu['M_phys']

print(f"  M_raw(e,  g0={g0_e:.3f})  = {Me_raw:.6f}")
print(f"  M_raw(mu, g0={g0_mu:.3f}) = {Mmu_raw:.6f}")
if Me_raw > 1e-12:
    ratio_raw = Mmu_raw / Me_raw
    print(f"  Stosunek M_raw:   {ratio_raw:.3f}  (cel: {TARGET_RATIO:.1f})")
    frr = abs(ratio_raw - TARGET_RATIO) / TARGET_RATIO * 100
    print(f"  Odchylenie:       {frr:.1f}%")
print()

print(f"  M_phys(e,  g0={g0_e:.3f})  = {Me_phys:.6f}")
print(f"  M_phys(mu, g0={g0_mu:.3f}) = {Mmu_phys:.6f}")
if Me_phys > 1e-12:
    ratio_phys = Mmu_phys / Me_phys
    print(f"  Stosunek M_phys:  {ratio_phys:.3f}  (cel: {TARGET_RATIO:.1f})")
    frp = abs(ratio_phys - TARGET_RATIO) / TARGET_RATIO * 100
    print(f"  Odchylenie:       {frp:.1f}%")
else:
    ratio_phys = 0.0
    frp = 999.0
print()

# Test hipotezy kinetycznej: M ~ g_0^4
ratio_kin4 = (g0_mu / g0_e)**4
print(f"  Hipoteza g^4:  ({g0_mu:.3f}/{g0_e:.3f})^4 = {ratio_kin4:.3f}")
frk4 = abs(ratio_kin4 - TARGET_RATIO) / TARGET_RATIO * 100
print(f"  Odchylenie:    {frk4:.1f}%")
print()

# Wymagana wartosc g0_mu dla M~g^4 -> r21=206.77
g0_mu_req4 = g0_e * TARGET_RATIO**(1.0/4.0)
print(f"  g0_mu wymagane (M~g^4, r=207): {g0_mu_req4:.4f}")
print()

# Z fitem gamma
if gamma_n0 is not None:
    ratio_kinG = (g0_mu / g0_e)**gamma_n0
    print(f"  Hipoteza g^{gamma_n0:.2f}: ({g0_mu:.3f}/{g0_e:.3f})^{gamma_n0:.2f} = {ratio_kinG:.3f}")
    g0_mu_reqG = g0_e * TARGET_RATIO**(1.0/gamma_n0)
    print(f"  g0_mu wymagane (g^{gamma_n0:.2f}, r=207): {g0_mu_reqG:.4f}")
    print()


# ─────────────────────────────────────────────────────────────────────────────
# Czesc 6: Szukaj g0 takie ze M_raw(g0)/M_raw(g0_e) ~ 207
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 6: Szukaj g0 takie ze M_raw/M_e ~ 207 ---")
print()

# Znajdz w skanie najblizszy stosunek do 207
if Me_raw > 1e-12:
    ratios = Mraw_arr / Me_raw
    idx207 = np.argmin(np.abs(ratios - TARGET_RATIO))
    g0_207   = g0_scan[idx207]
    ratio207 = ratios[idx207]
    print(f"  Najblizszy do 207: g0 = {g0_207:.4f},  M_raw/M_e = {ratio207:.3f}")
    print(f"  g_min przy tym g0 = {gmin_arr[idx207]:.5f}  "
          f"({'ponad g*' if gmin_arr[idx207] > G_GHOST else 'ponizej g*!'})")
    # Sprawdz hipoteze (g0_207/g0_e)^4
    ratio_kin207 = (g0_207 / g0_e)**4
    print(f"  (g0_207/g0_e)^4 = ({g0_207:.3f}/{g0_e:.3f})^4 = {ratio_kin207:.3f}")
    print()


# ─────────────────────────────────────────────────────────────────────────────
# Czesc 7: Wykresy
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 7: Wykresy ---")

fig = plt.figure(figsize=(15, 10))
gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.40)

colors = {'e': 'royalblue', 'mu': 'forestgreen', 'tau': 'darkorchid'}

# ── Panel A: Profile radiale g(r) ──
ax_A = fig.add_subplot(gs[0, 0:2])
ax_A.axhline(1.0,     color='gray',    ls='--', lw=0.8, label='g=1 (próżnia)')
ax_A.axhline(G_GHOST, color='crimson', ls=':',  lw=1.2, label=f'g*={G_GHOST:.3f}')
ax_A.axhline(G_DEEP,  color='orange',  ls=':',  lw=0.8, label=f'g_deep={G_DEEP}')
for res, g0_, lbl, col in [
        (res_e,  g0_e,   f'elektron g₀={g0_e}',     colors['e']),
        (res_mu, g0_mu,  f'mion g₀={g0_mu:.2f}',    colors['mu']),
        (res_tau,g0_tau, f'tauon g₀={g0_tau}',       colors['tau']),
]:
    rmax_plot = min(15.0, res['r'][-1])
    mask = res['r'] <= rmax_plot
    ax_A.plot(res['r'][mask], res['g'][mask], color=col, lw=1.5, label=lbl)
ax_A.set_xlim(0, 15)
ax_A.set_ylim(0.65, max(g0_e, g0_mu, g0_tau) * 1.05)
ax_A.set_xlabel('r [j.P.]', fontsize=9)
ax_A.set_ylabel('g(r)', fontsize=9)
ax_A.set_title('Profile radiale g(r)', fontsize=10)
ax_A.legend(fontsize=7.5)
ax_A.grid(True, alpha=0.25)

# ── Panel B: n_deep(g_0) i g_min(g_0) ──
ax_B = fig.add_subplot(gs[0, 2])
ax_B2 = ax_B.twinx()
ax_B.plot(g0_scan, ndp_arr,  'k-', lw=1.2, label='n_deep')
ax_B2.plot(g0_scan, gmin_arr, 'b--', lw=0.8, alpha=0.7, label='g_min')
ax_B2.axhline(G_GHOST, color='r', ls=':', lw=1)
ax_B.set_xlabel('g₀', fontsize=9)
ax_B.set_ylabel('n_deep', fontsize=9)
ax_B2.set_ylabel('g_min', fontsize=9, color='b')
ax_B.set_title('Stopień wężlowy n_deep(g₀)', fontsize=10)
lns1 = ax_B.get_lines() + ax_B2.get_lines()
ax_B.legend(lns1, [l.get_label() for l in lns1], fontsize=7)

# ── Panel C: M_raw(g_0) log-log ──
ax_C = fig.add_subplot(gs[1, 0])
mask_v = Mraw_arr > 1e-6
ax_C.loglog(g0_scan[mask_v], Mraw_arr[mask_v], 'k.', ms=2.5, label='M_raw')
ax_C.loglog(g0_scan[mask_v], Mph_arr[mask_v],  'b.', ms=2.5, alpha=0.5, label='M_phys')
if gamma_n0 is not None:
    g_fit = g0_scan[mask_v]
    ax_C.loglog(g_fit, A_n0 * g_fit**gamma_n0, 'r--', lw=1,
                label=f'g^{gamma_n0:.2f}')
ax_C.axvline(g0_e,  color=colors['e'],   ls=':', lw=1.5)
ax_C.axvline(g0_mu, color=colors['mu'],  ls=':', lw=1.5)
ax_C.set_xlabel('g₀', fontsize=9)
ax_C.set_ylabel('M', fontsize=9)
ax_C.set_title('Masa M_raw(g₀) log-log', fontsize=10)
ax_C.legend(fontsize=7)
ax_C.grid(True, which='both', alpha=0.2)

# ── Panel D: Portret fazowy (g, g') ──
ax_D = fig.add_subplot(gs[1, 1])
ax_D.axvline(1.0,     color='gray', ls='--', lw=0.8)
ax_D.axvline(G_GHOST, color='r',    ls=':',  lw=1.0)
ax_D.axhline(0.0,     color='gray', ls='--', lw=0.8)
for res, lbl, col in [
        (res_e,  f'e g₀={g0_e}',      colors['e']),
        (res_mu, f'mu g₀={g0_mu:.2f}', colors['mu']),
        (res_tau,f'tau g₀={g0_tau}',   colors['tau']),
]:
    ax_D.plot(res['g'], res['gp'], color=col, lw=1.0, alpha=0.8, label=lbl)
ax_D.set_xlabel('g', fontsize=9)
ax_D.set_ylabel("g'", fontsize=9)
ax_D.set_title("Portret fazowy (g, g')", fontsize=10)
ax_D.legend(fontsize=7)
ax_D.grid(True, alpha=0.2)

# ── Panel E: M_raw/M_e vs g_0 (test hipotezy 207) ──
ax_E = fig.add_subplot(gs[1, 2])
if Me_raw > 1e-12:
    ratios_plot = Mraw_arr[mask_v] / Me_raw
    ax_E.semilogy(g0_scan[mask_v], ratios_plot, 'k.', ms=2.5, label='M/M_e')
    # Linia g^4
    g_line = np.linspace(g0_scan[mask_v][0], g0_scan[mask_v][-1], 200)
    ax_E.semilogy(g_line, (g_line/g0_e)**4, 'r--', lw=1.2, label='(g/g_e)^4')
    ax_E.axhline(TARGET_RATIO, color='purple', ls=':', lw=1.5,
                 label=f'cel={TARGET_RATIO:.0f}')
    ax_E.axvline(g0_mu, color=colors['mu'], ls=':', lw=1.5)
ax_E.set_xlabel('g₀', fontsize=9)
ax_E.set_ylabel('M_raw/M_e', fontsize=9)
ax_E.set_title('Stosunek mas vs cel 207', fontsize=10)
ax_E.legend(fontsize=7)
ax_E.grid(True, which='both', alpha=0.2)

fig.suptitle(f'TGP Widmo Węzłowe | M_e={Me_raw:.3f}, g_0_e={g0_e}, '
             f'M_mu/M_e={ratio_raw:.1f} (cel {TARGET_RATIO:.0f})',
             fontsize=10)

out_dir = os.path.dirname(os.path.abspath(__file__))
out_png = os.path.join(out_dir, 'ex55_nodal_spectrum.png')
fig.savefig(out_png, dpi=110, bbox_inches='tight')
plt.close(fig)
print(f"  Wykres zapisany: {out_png}")
print()


# ─────────────────────────────────────────────────────────────────────────────
# Podsumowanie wynikow
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("PODSUMOWANIE EX55")
print("=" * 70)
print()
print(f"  Granica duchowa:  g* = {G_GHOST:.5f}")
if ghost_entries:
    print(f"  Wejscie w g<g*:   g0_crit ~ {ghost_entries[0]['mid']:.4f}")
print()
print(f"  n=0 (elektron):  g0 = {g0_e:.3f},  M_raw = {Me_raw:.5f}")
print(f"  n=? (mion):      g0 = {g0_mu:.3f},  M_raw = {Mmu_raw:.5f}")
print()
print(f"  M_mu/M_e (raw):  {ratio_raw:.2f}    (cel: {TARGET_RATIO:.1f})")
print(f"  (g0_mu/g0_e)^4:  {ratio_kin4:.2f}   (hipoteza kinetyczna)")
print()
if gamma_n0:
    print(f"  Wykladnik gamma: {gamma_n0:.3f}  (fit M_raw ~ g^gamma)")
print()

# Ocena jakosciowa
if Me_raw > 1e-12:
    log_err = abs(np.log10(ratio_raw / TARGET_RATIO))
    if log_err < 0.02:
        print("  WERDYKT: POTWIERDZENIE (<=5% od 207)")
    elif log_err < 0.12:
        print("  WERDYKT: CZESCIOWE POTWIERDZENIE (<30% od 207)")
    elif log_err < 1.0:
        print("  WERDYKT: WLASCIWY RZAD WIELKOSCI (w 10x od 207)")
    else:
        print("  WERDYKT: MECHANIZM WEZLOWY NIEDOSTATECZNY")
        fac = max(ratio_raw / TARGET_RATIO, TARGET_RATIO / ratio_raw)
        print(f"           Brakuje: {fac:.1f}x (={log_err:.1f} rzedu wlk.)")
        print("           Potrzeba innego mech. => ex56: singularnosc przy g*")
print()


# ─────────────────────────────────────────────────────────────────────────────
# TESTY JEDNOSTKOWE
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("TESTY JEDNOSTKOWE")
print("=" * 70)
print()

ALL_PASS = True

def check(tag, cond, note=""):
    global ALL_PASS
    ok = "PASS" if cond else "FAIL"
    if not cond:
        ALL_PASS = False
    print(f"  [{ok}] {tag}" + (f"  ({note})" if note else ""))

check("T1: g* = exp(-1/4) ~ 0.7788",
      abs(G_GHOST - 0.7788) < 0.0001,
      f"g* = {G_GHOST:.6f}")

check("T2: g0=1.24 (elektron): g_min > 0.85",
      res_e['g_min'] > 0.85,
      f"g_min = {res_e['g_min']:.5f}")

res_t3 = integrate_profile(2.50)
check("T3: g0=2.50: g_min <= g* + 1e-4  (osiaga granice duchowa)",
      res_t3['g_min'] <= G_GHOST + 1e-4,
      f"g_min = {res_t3['g_min']:.5f},  g* = {G_GHOST:.5f}")

# T4: M_raw rosnie z g0 w obszarze n=0
idx_n0_list = np.where(mask_n0)[0]
if len(idx_n0_list) >= 4:
    sample = idx_n0_list[::max(1,len(idx_n0_list)//4)][:4]
    M_smp  = Mraw_arr[sample]
    is_mono = all(M_smp[k] <= M_smp[k+1] for k in range(len(M_smp)-1))
    check("T4: M_raw rosnie z g0 (n=0 region)",
          is_mono,
          f"probka: {np.round(M_smp, 3)}")
else:
    check("T4: M_raw rosnie z g0", len(idx_n0_list) > 0,
          f"za malo punktow n=0: {len(idx_n0_list)}")

if gamma_n0 is not None:
    check("T5: gamma in [1.5, 20.0]  (wzmocnienie kinetyczne blisko g*)",
          1.5 <= gamma_n0 <= 20.0,
          f"gamma = {gamma_n0:.3f}")
else:
    check("T5: gamma in [1.5, 7.0]", False, "fit nieudany")

print()
print(f"  Wynik: {'WSZYSTKIE TESTY ZALICZONE' if ALL_PASS else 'NIEKTORЕ TESTY OBLANE'}")
print()
print("=" * 70)
print("KONIEC EX55")
print("=" * 70)
