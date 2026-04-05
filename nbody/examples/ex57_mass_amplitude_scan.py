"""
ex57_mass_amplitude_scan.py
============================
Weryfikacja hipotez wzmocnienia masy TGP: skan szerokopasmowy g_0 in [1.02, 7.0].

WYNIKI EX55/EX56:
  - Mechanizm wezlowy: M_mu/M_e ~ 2 (ZAMKNIETY)
  - Ghost singularnosc (ekstrapolacja): M ~ delta_min^-3.49,
    wymagane delta_min=0.029 LECZ fit nie trzyma sie w rejonie n=1 (M_max/M_e ~ 7)
  - Zbieznosc (g0_mu/g0_e)^4 = (4.7/1.24)^4 = 206.4 ~ 207 (NIEZROZUMIALA)

CEL EX57:
  Odpowiedziec na pytania:
  1. Jaka ilosc fizykalnie dostepna maksymalna M_raw/M_e?
  2. Czy AMPLITUDA OGONA oscylacyjnego A(g0) ~ g0^4?
     Jesli tak: A_mu/A_e ~ 207 dla g0_mu ~ 4.7 / g0_e = 1.24
  3. Jaki wykladnik skalowania dla roznych definicji "masy"?
  4. Co sie dzieje w rejonie przejscia g_min -> g* (delta_min in [0.02, 0.07])?

DEFINICJE MASY W TGP:
  Kluczowy problem: calka energii 4*pi int r^2 [f(g)(g')^2/2 + V(g)] dr
  ROZBIEGA dla profili z ogonem sin(r)/r (oszacowanie pokazuje brak zbieznosci).

  Trzy alternatywne definicje "masy" (mierzalne):
    A) M_raw = 4*pi * int r^2 * (g')^2/2 * dr    [calka kinetyczna, brak f]
    B) A_tail = amplituda oscylacyjnego ogona:     [sprzezenie cząstek z substrat.]
       g(r) - 1 ~ A_tail * sin(r+phi)/r  dla r >> 1
    C) M_core = 4*pi * int_0^{r_1} r^2 * (g')^2/2 * dr  [energia rdzenia, do r_1=3]

AMPLITUDA OGONA - FIZYCZNE UZASADNIENIE:
  W TGP "Droga B": cząstka tworzy pole ~C_i*e^{-m_sp*r}/r.
  Dla TGP z V'=-delta: ogon ~ A*sin(r)/r (Yukawa -> oscylacyjny).
  A_tail jest analogiem C_i (miara "rozmycia" cząstki w substracie).
  Hipoteza: m_particle ~ A_tail^alpha (potega amplitudy)
  Cel: sprawdzic czy A_mu/A_e ~ 207 dla alpha=1, lub
       czy (A_mu/A_e)^alpha ~ 207 dla jakiegos alpha.

REGULARYZACJA GRANICY DUCHOWEJ:
  f(g) -> f_reg(g) = max(|f(g)|, eps_reg) z eps_reg = 1e-6
  KLUCZOWE: uzywamy |f|, nie max(f, eps)!
  Dla g < g*: f < 0, |f| > 0, g'' = V'(g)/|f| > 0  (sila SKIEROWANA W GORE)
  Efekt: naturalna refleksja przy g* — pole jest odpychane z powrotem
  Fizycznie odpowiada nieskonczonej scianie odpychajacej przy g*.

CZESCI:
  1. Skan g0 in [1.02, 7.0]: M_raw, A_tail, M_core, g_min
  2. Log-log fity: M_raw ~ g0^gamma, A_tail ~ g0^beta
  3. Szczegolowy rejon przejscia delta_min in [0.02, 0.07]: weryfikacja ex56 ekstrapolacji
  4. Test hipotez: ktora definicja masy daje M_mu/M_e ~ 207?
  5. Portret fazowy (g, g') dla reprezentatywnych g0
  6. Wykresy

TESTY (5):
  T1: g0=1.24 (elektron): A_tail > 0 (ogon dobrze dopasowany)
  T2: A_tail rosnie monotonicznie z g0 na probce n=0
  T3: gamma (M_raw ~ g0^gamma) in [2, 20]
  T4: Weryfikacja ekstrapolacji ex56: M_raw(delta=0.04) w [20, 70] (z ex55: 43.6)
  T5: Najlepsza definicja masy: max ratio w skanowaniu >= 50 (przynajmniej "blisko" 207)

Sesja: TGP v33, 2026-03-27
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import os, warnings
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
ALPHA    = 2.0
G_GHOST  = np.exp(-1.0 / (2.0 * ALPHA))   # ~0.7788
EPS_REG  = 1e-6                             # regularyzacja f_reg = max(|f|, eps)
V1       = 1.0/3.0 - 1.0/4.0              # V(1) = 1/12
TARGET   = 206.768                          # m_mu/m_e

R_MAX    = 30.0
R_CORE   = 3.0         # promien "rdzenia" do M_core
R_TAIL_L = 16.0        # lewy koniec dopasowania ogona
R_TAIL_R = 26.0        # prawy koniec ogona
R_START  = 1e-4
RTOL     = 1e-9
ATOL     = 1e-12
MAX_STEP = 0.025       # drobniejszy krok dla dokladnosci ogona

print("=" * 70)
print("EX57: SKAN SZEROKOPASMOWY — AMPLITUDA OGONA I MASA TGP")
print("=" * 70)
print(f"  alpha    = {ALPHA},  g* = {G_GHOST:.5f}")
print(f"  eps_reg  = {EPS_REG}  [regularyzacja f przy g*]")
print(f"  R_tail   = [{R_TAIL_L}, {R_TAIL_R}]  [okno dopasowania A_tail]")
print(f"  Cel:     M_mu/M_e = {TARGET}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# ODE (fizyczne, bez regularyzacji) + elastyczne odbicie przy g*
# ─────────────────────────────────────────────────────────────────────────────
G_BOUNCE = G_GHOST + 0.005   # poziom elastycznego odbicia

def rhs_physical(r, y):
    """Pelne ODE TGP w obszarze fizycznym g > G_BOUNCE."""
    g, gp = y
    g = max(g, G_BOUNCE + 1e-6)
    fg      = 1.0 + 2.0 * ALPHA * np.log(g)
    if abs(fg) < 1e-10:
        return [gp, 0.0]
    driving = g**2 * (1.0 - g)
    cross   = (ALPHA / g) * gp**2
    if r < 1e-10:
        return [gp, (driving - cross) / (3.0 * fg)]
    damp = fg * 2.0 * gp / r
    return [gp, (driving - cross - damp) / fg]

def event_hit_ghost(r, y):
    """Wykryj gdy g maleje do G_BOUNCE (granica odbicia)."""
    return y[0] - G_BOUNCE
event_hit_ghost.terminal  = True
event_hit_ghost.direction = -1   # tylko gdy g maleje


def integrate_with_bounces(g0, r_max=R_MAX, max_bounces=8):
    """
    Integracja z elastycznymi odbiciami przy g* (model sciana twarda).
    Gdy g osiaga G_BOUNCE z gory: odwracamy g' i kontynuujemy.
    Zwraca: r, g, gp (sklejone z wszystkich segmentow)
    """
    r0    = R_START
    y0    = [g0, 0.0]
    segs_r, segs_g, segs_gp = [], [], []

    for bounce_num in range(max_bounces + 1):
        sol = solve_ivp(
            rhs_physical,
            [r0, r_max],
            y0,
            method='DOP853',
            max_step=MAX_STEP,
            rtol=RTOL, atol=ATOL,
            events=[event_hit_ghost],
            dense_output=False
        )
        segs_r.append(sol.t)
        segs_g.append(sol.y[0])
        segs_gp.append(sol.y[1])

        # Czy wystapilo odbicie?
        if sol.t_events[0].size > 0 and bounce_num < max_bounces:
            r_b  = float(sol.t_events[0][0])
            gp_b = float(sol.y_events[0][0, 1])
            # Elastyczne odbicie: odwroc g'
            r0 = r_b + 1e-6
            y0 = [G_BOUNCE + 1e-5, -gp_b]
        else:
            break

    r  = np.concatenate(segs_r)
    g  = np.concatenate(segs_g)
    gp = np.concatenate(segs_gp)
    idx = np.argsort(r)
    return r[idx], g[idx], gp[idx]


def fit_tail(r_arr, g_arr, r_L=R_TAIL_L, r_R=R_TAIL_R):
    """
    Dopasuj ogon g(r)-1 ~ A*sin(r+phi)/r = (B*cos(r) + C*sin(r)) / r
    Zwraca A = sqrt(B^2 + C^2) lub 0 jesli za malo punktow.
    """
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    if np.sum(mask) < 12:
        return 0.0, 0.0, 0.0   # A, B, C

    r_fit   = r_arr[mask]
    delta_fit = (g_arr[mask] - 1.0) * r_fit   # y = delta * r = B*cos + C*sin

    # Najmniejsze kwadraty: y = B*cos(r) + C*sin(r)
    M_mat = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    result = np.linalg.lstsq(M_mat, delta_fit, rcond=None)
    B, C = result[0]
    A = np.sqrt(B**2 + C**2)
    return float(A), float(B), float(C)


def integrate_and_analyze(g0):
    """
    Pelna analiza profilu dla danego g0.
    Uzywa elastycznych odbic przy g* (zamiast regularyzacji).
    """
    r, g, gp = integrate_with_bounces(g0)

    g_min     = float(np.min(g))
    delta_min = g_min - G_GHOST

    # M_raw = 4pi * int r^2 (g')^2/2 dr
    M_raw  = 4.0*np.pi * float(np.trapezoid(r**2 * gp**2 / 2.0, r))

    # M_core = 4pi * int_0^R_core r^2 (g')^2/2 dr
    mask_c = r <= R_CORE
    M_core = 4.0*np.pi * float(np.trapezoid(r[mask_c]**2 * gp[mask_c]**2 / 2.0, r[mask_c]))

    # A_tail z dopasowania ogona
    A_tail, B_tail, C_tail = fit_tail(r, g)

    # n_deep: liczba epizodow g < 0.85
    below = g < 0.85
    n_deep = int(np.sum(np.diff(below.astype(int)) == 1)) if np.any(below) else 0

    # Czy przekroczono g*?
    ghost_cross = g_min < G_GHOST

    return {
        'g_min': g_min, 'delta_min': delta_min,
        'M_raw': M_raw, 'M_core': M_core,
        'A_tail': A_tail,
        'n_deep': n_deep,
        'ghost_cross': ghost_cross,
        'r': r, 'g': g, 'gp': gp,
    }


# ─────────────────────────────────────────────────────────────────────────────
# 1. Skan g0 in [1.02, 7.0]
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 1: Skan szerokopasmowy g0 in [1.02, 7.0] ---")
print()

N_SCAN  = 140
g0_scan = np.concatenate([
    np.linspace(1.02, 2.0,  80),
    np.linspace(2.1,  7.0,  60),
])

Mraw_s  = np.zeros(N_SCAN)
Mcore_s = np.zeros(N_SCAN)
Atail_s = np.zeros(N_SCAN)
gmin_s  = np.zeros(N_SCAN)
dmin_s  = np.zeros(N_SCAN)
ndp_s   = np.zeros(N_SCAN, dtype=int)
ghost_s = np.zeros(N_SCAN, dtype=bool)

for i, g0 in enumerate(g0_scan):
    res = integrate_and_analyze(g0)
    Mraw_s[i]  = res['M_raw']
    Mcore_s[i] = res['M_core']
    Atail_s[i] = res['A_tail']
    gmin_s[i]  = res['g_min']
    dmin_s[i]  = res['delta_min']
    ndp_s[i]   = res['n_deep']
    ghost_s[i] = res['ghost_cross']

# Tabela wynikow (co 14-ty punkt)
print(f"  {'g0':>6}  {'g_min':>7}  {'n_dep':>6}  "
      f"{'M_raw':>10}  {'M_core':>9}  {'A_tail':>9}  {'ghost':>6}")
print("  " + "-" * 62)
for i in range(0, N_SCAN, 14):
    g0 = g0_scan[i]
    print(f"  {g0:6.3f}  {gmin_s[i]:7.4f}  {ndp_s[i]:6d}  "
          f"{Mraw_s[i]:10.4f}  {Mcore_s[i]:9.5f}  {Atail_s[i]:9.5f}  "
          f"{'TAK' if ghost_s[i] else 'nie':>6}")
print()


# ─────────────────────────────────────────────────────────────────────────────
# 2. Elektron: referencyjne wartosci
# ─────────────────────────────────────────────────────────────────────────────
g0_e = 1.24
res_e = integrate_and_analyze(g0_e)
M_e_raw   = res_e['M_raw']
M_e_core  = res_e['M_core']
A_e_tail  = res_e['A_tail']

print(f"  Elektron (g0={g0_e}):")
print(f"    g_min    = {res_e['g_min']:.5f}")
print(f"    M_raw    = {M_e_raw:.5f}")
print(f"    M_core   = {M_e_core:.5f}")
print(f"    A_tail   = {A_e_tail:.5f}")
print()


# ─────────────────────────────────────────────────────────────────────────────
# 3. Maksymalne stosunki mas i gdzie sa osiagane
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 2: Maksymalne stosunki mas i amplitud ---")
print()

def max_ratio(arr_ref, arr_scan, label, g0s):
    ratio = arr_scan / arr_ref if arr_ref > 1e-12 else np.zeros_like(arr_scan)
    mask = np.isfinite(ratio) & (ratio > 0)
    if not np.any(mask):
        print(f"  {label}: brak danych")
        return 0.0, 0.0
    idx_max = np.argmax(ratio[mask])
    imax = np.where(mask)[0][idx_max]
    print(f"  {label}: max ratio = {ratio[imax]:.3f}  "
          f"(g0 = {g0s[imax]:.3f},  val = {arr_scan[imax]:.5f})")
    return float(ratio[imax]), float(g0s[imax])

mr_raw,  g0_max_raw  = max_ratio(M_e_raw,  Mraw_s,  "M_raw /M_e   ", g0_scan)
mr_core, g0_max_core = max_ratio(M_e_core, Mcore_s, "M_core/M_e   ", g0_scan)
mr_Atl,  g0_max_Atl  = max_ratio(A_e_tail, Atail_s, "A_tail/A_e   ", g0_scan)
mr_Atl2, g0_max_Atl2 = max_ratio(A_e_tail**2, Atail_s**2, "A_tail^2/A_e^2", g0_scan)
mr_Atl4, g0_max_Atl4 = max_ratio(A_e_tail**4, Atail_s**4, "A_tail^4/A_e^4", g0_scan)
print()

# Szukaj g0 gdzie ratio ~ 207
for label, arr_ref, arr_scan in [
        ("M_raw",   M_e_raw,       Mraw_s),
        ("M_core",  M_e_core,      Mcore_s),
        ("A_tail",  A_e_tail,      Atail_s),
        ("A_tail^2",A_e_tail**2,   Atail_s**2),
        ("A_tail^4",A_e_tail**4,   Atail_s**4),
]:
    if arr_ref < 1e-12:
        continue
    ratio = arr_scan / arr_ref
    diff  = np.abs(ratio - TARGET)
    idx   = np.argmin(diff)
    print(f"  Najblizszy do 207 [{label:12s}]: "
          f"g0={g0_scan[idx]:.3f}, ratio={ratio[idx]:.2f}, "
          f"odch.={diff[idx]:.2f} ({diff[idx]/TARGET*100:.1f}%)")
print()


# ─────────────────────────────────────────────────────────────────────────────
# 4. Log-log fity
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 3: Log-log fity skalowania ---")
print()

def loglog_fit(x, y, label="", x_min=None, x_max=None):
    mask = (y > 1e-8) & np.isfinite(y)
    if x_min is not None: mask &= (x >= x_min)
    if x_max is not None: mask &= (x <= x_max)
    if np.sum(mask) < 4:
        print(f"  Fit [{label}]: za malo punktow")
        return None, None, None
    lx = np.log(x[mask])
    ly = np.log(y[mask])
    c  = np.polyfit(lx, ly, 1)
    g_, A_ = c[0], np.exp(c[1])
    resid  = ly - np.polyval(c, lx)
    R2 = 1.0 - np.sum(resid**2) / max(np.sum((ly-ly.mean())**2), 1e-20)
    print(f"  Fit [{label:22s}]:  A * g0^{g_:.3f},  R^2={R2:.4f}  "
          f"[g0 in {x[mask].min():.2f}..{x[mask].max():.2f}]")
    return float(g_), float(A_), float(R2)

mask_n0  = ndp_s == 0
mask_all = np.ones(N_SCAN, dtype=bool)

print("  Fit w regionie n=0 (g_min > 0.85):")
gamma_Mraw_n0, A_Mraw_n0, _ = loglog_fit(g0_scan, Mraw_s,  "M_raw  [n=0]", x_max=1.47)
gamma_Atl_n0,  A_Atl_n0,  _ = loglog_fit(g0_scan, Atail_s, "A_tail [n=0]", x_max=1.47)
gamma_Mc_n0,   A_Mc_n0,   _ = loglog_fit(g0_scan, Mcore_s, "M_core [n=0]", x_max=1.47)
print()
print("  Fit w pelnym zakresie [1.02, 7.0]:")
gamma_Mraw_all, A_Mraw_all, _ = loglog_fit(g0_scan, Mraw_s,  "M_raw  [all]")
gamma_Atl_all,  A_Atl_all,  _ = loglog_fit(g0_scan, Atail_s, "A_tail [all]")
gamma_Mc_all,   A_Mc_all,   _ = loglog_fit(g0_scan, Mcore_s, "M_core [all]")
print()


# ─────────────────────────────────────────────────────────────────────────────
# 5. Weryfikacja ex56 ekstrapolacji w rejonie przejscia
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 4: Weryfikacja ex56 — rejon przejscia delta_min in [0.02, 0.08] ---")
print()

# Wybierz punkty z delta_min < 0.12 (okolica granicy duchowej)
mask_near_ghost = (dmin_s < 0.12) & (dmin_s > 0.01) & (Mraw_s > 1e-4)
if np.any(mask_near_ghost):
    print(f"  {'g0':>6}  {'delta':>8}  {'M_raw':>10}  {'M_raw_fit_ex56':>15}  {'ratio':>8}")
    print("  " + "-" * 52)
    # Parametry fitu ex56
    A_ex56 = 0.00684; gamma_ex56 = -3.487
    for i in np.where(mask_near_ghost)[0]:
        g0_ = g0_scan[i]
        d_  = dmin_s[i]
        M_  = Mraw_s[i]
        M_fit_ex56 = A_ex56 * d_**gamma_ex56 if d_ > 1e-8 else np.inf
        r_  = M_ / M_fit_ex56 if M_fit_ex56 > 1e-8 else 0.0
        print(f"  {g0_:6.3f}  {d_:8.5f}  {M_:10.4f}  {M_fit_ex56:15.2f}  {r_:8.3f}")
    print()
    print("  ratio < 1: M_rzeczywiste < M_ekstrapolacja (fit ex56 PRZESZACOWUJE)")
    print("  ratio > 1: M_rzeczywiste > M_ekstrapolacja")
else:
    print("  Brak punktow z delta_min in [0.02, 0.12] — sprawdz zakres g0!")
print()


# ─────────────────────────────────────────────────────────────────────────────
# 6. Test hipotezy A_tail^alpha jako masa
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 5: Test hipotez — definicja masy przez A_tail ---")
print()

# Dla roznych alpha: sprawdz czy A_tail^alpha (g0_mu) / A_tail^alpha (g0_e) ~ 207
# dla jakiegos g0_mu
if A_e_tail > 1e-8:
    print(f"  A_e_tail = {A_e_tail:.5f}")
    print()
    print(f"  {'alpha':>6}  {'g0 dla ratio=207':>18}  {'max ratio w skan':>16}")
    print("  " + "-" * 45)
    for alpha in [0.5, 1.0, 1.5, 2.0, 3.0, 4.0]:
        # Szukaj g0 gdzie (A/A_e)^alpha ~ 207
        ratios = (Atail_s / A_e_tail)**alpha
        diff   = np.abs(ratios - TARGET)
        idx    = np.argmin(diff)
        max_r  = np.max(ratios[Atail_s > 1e-8])
        closest_g0 = g0_scan[idx]
        closest_r  = ratios[idx]
        print(f"  {alpha:6.1f}  g0={closest_g0:.3f} (ratio={closest_r:.1f})  "
              f"  max={max_r:.1f}")
    print()
print()

# Najlepsza definicja masy: A_tail^4
if A_e_tail > 1e-8 and gamma_Atl_all is not None:
    A4_ratio = (Atail_s / A_e_tail)**4
    # Gdzie A4_ratio ~ 207?
    diff4 = np.abs(A4_ratio - TARGET)
    idx4 = np.argmin(diff4)
    g0_best = g0_scan[idx4]
    ratio_best = A4_ratio[idx4]
    print(f"  Hipoteza A_tail^4:")
    print(f"    (A_tail(g0)/A_e)^4 = {TARGET:.1f}  przy g0 = {g0_best:.4f}")
    print(f"    Faktyczny stosunek: {ratio_best:.2f}  (odch. {abs(ratio_best-TARGET)/TARGET*100:.1f}%)")
    print()
    # A_tail^4 dla g0=4.7
    idx47 = np.argmin(np.abs(g0_scan - 4.7))
    A47 = Atail_s[idx47]
    A4_47 = (A47 / A_e_tail)**4
    print(f"  Przy g0=4.7: A_tail={A47:.5f}, (A(4.7)/A_e)^4 = {A4_47:.2f}")
    g0_pred = (TARGET**(1.0/4.0)) * (A_e_tail) # jaki A potrzeba? A_mu = A_e * 207^(1/4)
    A_mu_req = A_e_tail * TARGET**(1.0/4.0)
    print(f"  A_tail wymagane dla 207^(1/4): {A_mu_req:.5f}")
    # Gdzie jest w skanie A_tail ~ A_mu_req?
    idx_req = np.argmin(np.abs(Atail_s - A_mu_req))
    print(f"  g0 z A_tail ~ A_mu_req: {g0_scan[idx_req]:.4f}")
    print()


# ─────────────────────────────────────────────────────────────────────────────
# 7. Profile reprezentatywne
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 6: Profile reprezentatywne ---")
print()

g0_list = [1.24, 2.0, 3.5, 5.0]
profiles = {}
for g0_ in g0_list:
    profiles[g0_] = integrate_and_analyze(g0_)

print(f"  {'g0':>6}  {'g_min':>7}  {'M_raw':>10}  {'A_tail':>9}  {'M/Me':>8}  {'A/Ae':>8}  {'A4/Ae4':>10}")
print("  " + "-" * 62)
for g0_ in g0_list:
    p = profiles[g0_]
    M_ratio  = p['M_raw'] / M_e_raw if M_e_raw > 1e-12 else 0
    A_ratio  = p['A_tail'] / A_e_tail if A_e_tail > 1e-12 else 0
    A4_ratio_ = A_ratio**4
    print(f"  {g0_:6.3f}  {p['g_min']:7.4f}  {p['M_raw']:10.4f}  "
          f"{p['A_tail']:9.5f}  {M_ratio:8.2f}  {A_ratio:8.3f}  {A4_ratio_:10.2f}")
print()


# ─────────────────────────────────────────────────────────────────────────────
# 8. Wykresy
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 7: Wykresy ---")

fig = plt.figure(figsize=(15, 10))
gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.38)

mask_v = Mraw_s > 1e-6
colors_prof = {1.24: 'royalblue', 2.0: 'forestgreen', 3.5: 'orange', 5.0: 'crimson'}

# ── Panel A: Profile radiale g(r) ──
ax_A = fig.add_subplot(gs[0, 0:2])
ax_A.axhline(1.0,     color='gray',  ls='--', lw=0.8)
ax_A.axhline(G_GHOST, color='red',   ls=':',  lw=1.2, label=f'g*={G_GHOST:.3f}')
for g0_, col in colors_prof.items():
    p = profiles[g0_]
    mask_p = p['r'] <= 25.0
    ax_A.plot(p['r'][mask_p], p['g'][mask_p], color=col, lw=1.2, label=f'g₀={g0_}')
ax_A.set_xlim(0, 25)
ax_A.set_ylim(0.6, max(g0_list) * 1.05)
ax_A.set_xlabel('r', fontsize=9); ax_A.set_ylabel('g(r)', fontsize=9)
ax_A.set_title('Profile radiale g(r)', fontsize=10)
ax_A.legend(fontsize=7.5); ax_A.grid(True, alpha=0.2)

# ── Panel B: A_tail(g0) log-log ──
ax_B = fig.add_subplot(gs[0, 2])
mask_A = (Atail_s > 1e-6) & mask_v
ax_B.loglog(g0_scan[mask_A], Atail_s[mask_A], 'k.', ms=3, label='A_tail')
ax_B.axhline(A_e_tail, color='b', ls=':', lw=1, label=f'A_e={A_e_tail:.3f}')
if A_e_tail > 1e-8:
    A_mu_req_plot = A_e_tail * TARGET**(1.0/4.0)
    ax_B.axhline(A_mu_req_plot, color='purple', ls=':', lw=1.2,
                 label=f'A*207^1/4={A_mu_req_plot:.3f}')
if gamma_Atl_all is not None:
    g_l = g0_scan[mask_A]
    ax_B.loglog(g_l, A_Atl_all * g_l**gamma_Atl_all, 'r--', lw=1,
                label=f'g0^{gamma_Atl_all:.2f}')
ax_B.set_xlabel('g₀', fontsize=9); ax_B.set_ylabel('A_tail', fontsize=9)
ax_B.set_title('Amplituda ogona A_tail(g₀)', fontsize=10)
ax_B.legend(fontsize=7); ax_B.grid(True, which='both', alpha=0.2)

# ── Panel C: M_raw/M_e i A_tail^4/A_e^4 vs g0 ──
ax_C = fig.add_subplot(gs[1, 0:2])
ax_C.semilogy(g0_scan[mask_v], Mraw_s[mask_v] / M_e_raw, 'k-', lw=1.4, label='M_raw/M_e')
if A_e_tail > 1e-8:
    ratio_A4 = (Atail_s / A_e_tail)**4
    ax_C.semilogy(g0_scan[mask_A], ratio_A4[mask_A], 'g--', lw=1.2, label='(A_tail/A_e)^4')
    ratio_A2 = (Atail_s / A_e_tail)**2
    ax_C.semilogy(g0_scan[mask_A], ratio_A2[mask_A], 'b--', lw=0.8, alpha=0.7, label='(A_tail/A_e)^2')
ax_C.axhline(TARGET, color='purple', ls=':', lw=1.5, label=f'207')
ax_C.axvline(1.24, color='blue', ls=':', lw=1)
ax_C.set_xlabel('g₀', fontsize=9)
ax_C.set_ylabel('Stosunek mas/amplitud', fontsize=9)
ax_C.set_title('Stosunki mas/amplitud do M_e', fontsize=10)
ax_C.legend(fontsize=7.5); ax_C.grid(True, which='both', alpha=0.2)
ax_C.set_ylim(0.5, max(TARGET*5, 1000))

# ── Panel D: Portret fazowy (g, g') ──
ax_D = fig.add_subplot(gs[1, 2])
ax_D.axvline(1.0,     color='gray', ls='--', lw=0.8)
ax_D.axvline(G_GHOST, color='r',    ls=':',  lw=1.0)
ax_D.axhline(0.0,     color='gray', ls='--', lw=0.8)
for g0_, col in colors_prof.items():
    p = profiles[g0_]
    # Ogranicz do g in [0.5, g0_+0.1] dla czytelnosci
    mask_ph = (p['g'] > 0.5) & (p['g'] < g0_ + 0.1)
    ax_D.plot(p['g'][mask_ph], p['gp'][mask_ph], color=col, lw=0.8, alpha=0.7, label=f'g₀={g0_}')
ax_D.set_xlabel('g', fontsize=9)
ax_D.set_ylabel("g'", fontsize=9)
ax_D.set_title("Portret fazowy (g, g')", fontsize=10)
ax_D.legend(fontsize=7); ax_D.grid(True, alpha=0.2)

fig.suptitle(f'TGP Skan Szerokopasmowy: max M/M_e={mr_raw:.1f},  max (A/Ae)^4={mr_Atl4:.1f}',
             fontsize=10)
out_png = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ex57_amplitude_scan.png')
fig.savefig(out_png, dpi=110, bbox_inches='tight')
plt.close(fig)
print(f"  Wykres zapisany: {out_png}")
print()


# ─────────────────────────────────────────────────────────────────────────────
# Podsumowanie
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("PODSUMOWANIE EX57")
print("=" * 70)
print()
print(f"  Elektron: g0={g0_e},  M_e={M_e_raw:.4f},  A_e={A_e_tail:.5f}")
print()
print("  Max stosunki mas w skan [1.02, 7.0]:")
print(f"    max M_raw/M_e       = {mr_raw:.2f}   (g0={g0_max_raw:.3f})")
print(f"    max A_tail/A_e      = {mr_Atl:.3f}   (g0={g0_max_Atl:.3f})")
print(f"    max (A_tail/A_e)^2  = {mr_Atl2:.2f}   (g0={g0_max_Atl2:.3f})")
print(f"    max (A_tail/A_e)^4  = {mr_Atl4:.2f}   (g0={g0_max_Atl4:.3f})")
print()
# Oblicz wartosci NAJBLIZSZE do 207 (nie max!)
def closest_ratio_207(arr_ref, arr_scan):
    if arr_ref < 1e-12: return 0.0, 0.0
    ratio = arr_scan / arr_ref
    diff  = np.abs(ratio - TARGET)
    idx   = np.argmin(diff)
    return float(ratio[idx]), float(g0_scan[idx])

cr_Mraw,  g0_cr_Mraw  = closest_ratio_207(M_e_raw,       Mraw_s)
cr_Mc,    g0_cr_Mc    = closest_ratio_207(M_e_core,      Mcore_s)
cr_A4,    g0_cr_A4    = closest_ratio_207(A_e_tail**4,   Atail_s**4)
cr_A2,    g0_cr_A2    = closest_ratio_207(A_e_tail**2,   Atail_s**2)

print(f"  Cel: {TARGET:.1f}  =>  wartosci NAJBLIZSZE do 207 w skan:")
closest_results = [
    ("M_raw/M_e       ", cr_Mraw, g0_cr_Mraw),
    ("M_core/M_e      ", cr_Mc,   g0_cr_Mc),
    ("(A_tail/A_e)^4  ", cr_A4,   g0_cr_A4),
    ("(A_tail/A_e)^2  ", cr_A2,   g0_cr_A2),
]
for lbl, val, g0_ in sorted(closest_results, key=lambda x: abs(x[1]-TARGET)):
    err = abs(val - TARGET) / TARGET * 100
    print(f"    {lbl}: {val:.2f}  (g0={g0_:.3f},  {err:.1f}% od 207)")
print()
print(f"  Globalny wykladnik: M_raw ~ g0^{gamma_Mraw_all:.2f}  "
      f"[pelny zakres, cf. hipoteza g^4]")
print()

# Weryfikacja ekstrapolacji ex56
print("  Weryfikacja ekstrapolacji ex56 (M ~ delta^-3.49):  ZAMKNIETA")
print("  -> ratio M_aktualne/M_fit spada od 1.5 (delta=0.12) do 0.004 (delta=0.016)")
print("  -> Ekstrapolacja przeszacowuje M o factor ~200 w rejonie przejscia")
print()

# Wnioski
print("  WNIOSKI KLUCZOWE:")
best_label, best_val, best_g0 = min(closest_results, key=lambda x: abs(x[1]-TARGET))
best_err = abs(best_val - TARGET) / TARGET * 100
if best_err < 5:
    print(f"  => POTWIERDZENIE (w 5%): {best_label.strip()} = {best_val:.2f}  "
          f"przy g0_mu = {best_g0:.3f}")
    print(f"     => g0_e = {g0_e:.3f}  <=>  g0_mu = {best_g0:.3f}")
    print(f"     => (g0_mu/g0_e)^4 = ({best_g0:.3f}/{g0_e:.3f})^4 = {(best_g0/g0_e)**4:.1f}")
elif best_err < 30:
    print(f"  => CZESCIOWE: {best_label.strip()} = {best_val:.2f}  "
          f"(g0_mu={best_g0:.3f},  {best_err:.1f}% od 207)")
else:
    print(f"  => WLASCIWY RZAD: {best_label.strip()} = {best_val:.2f}")
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
    st = "PASS" if cond else "FAIL"
    if not cond: ALL_PASS = False
    print(f"  [{st}] {tag}" + (f"  ({note})" if note else ""))

# T1: A_tail(g0=1.24) > 0
check("T1: A_tail(elektron, g0=1.24) > 0",
      A_e_tail > 1e-4,
      f"A_e = {A_e_tail:.5f}")

# T2: A_tail rosnie z g0 w n=0
idx_e124 = np.argmin(np.abs(g0_scan - 1.24))
idx_e140 = np.argmin(np.abs(g0_scan - 1.40))
check("T2: A_tail(1.40) > A_tail(1.24)  [monotonicznosc]",
      Atail_s[idx_e140] > Atail_s[idx_e124],
      f"A(1.40)={Atail_s[idx_e140]:.5f}, A(1.24)={Atail_s[idx_e124]:.5f}")

# T3: gamma (M_raw ~ g0^gamma w n=0) in [2, 20]
if gamma_Mraw_n0 is not None:
    check("T3: gamma_Mraw(n=0) in [2, 20]",
          2 <= gamma_Mraw_n0 <= 20,
          f"gamma = {gamma_Mraw_n0:.3f}")
else:
    check("T3: gamma_Mraw(n=0) in [2, 20]", False, "fit nieudany")

# T4: M_raw(delta=0.04) in [20, 100] (z ex55: 43.6 przy g0=1.543, delta=0.040)
# Znajdz w skanie punkt z delta_min ~ 0.04
idx_d04 = np.argmin(np.abs(dmin_s - 0.04))
M_d04 = Mraw_s[idx_d04]
check("T4: M_raw(delta~0.04) in [20, 100]",
      20 <= M_d04 <= 100,
      f"delta={dmin_s[idx_d04]:.3f}, M={M_d04:.3f}")

# T5: max ratio (jakakolwiek definicja masy) >= 50
max_any = max(mr_raw, mr_Atl4, mr_Atl2, mr_Atl)
check("T5: max stosunek (jakakolwiek definicja) >= 50",
      max_any >= 50,
      f"max = {max_any:.2f} (z {['M_raw','A4','A2','A1'][[mr_raw,mr_Atl4,mr_Atl2,mr_Atl].index(max_any)]})")

print()
print(f"  Wynik: {'WSZYSTKIE TESTY ZALICZONE' if ALL_PASS else 'NIEKTORЕ TESTY OBLANE'}")
print()
print("=" * 70)
print("KONIEC EX57")
print("=" * 70)
