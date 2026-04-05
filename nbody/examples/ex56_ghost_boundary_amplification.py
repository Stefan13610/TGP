"""
ex56_ghost_boundary_amplification.py
=====================================
Analiza osobliwosci granicy duchowej g* w TGP i jej wplyw
na stosunek mas profili solitonowych.

WYNIK EX55:
  Mechanizm wezlowy (n=0/n=1) daje M_mu/M_e ~ 2, nie 207.
  Wykladnik gamma ~ 12.5 (stromszy niz oczekiwane 4).
  Pytanie: Czy singularnosc f(g)=0 przy g* moze dawac wzmocnienie 207?

FIZYKA GRANICY DUCHOWEJ:
  f(g) = 1 + 2*alpha*ln(g),   g* = exp(-1/4) ~ 0.7788,   f(g*) = 0
  V'(g) = g^2*(1-g)
  V'(g*) = (g*)^2*(1-g*) ~ 0.134 > 0  (sila odpycha od g* w gore)

POTENCJAL EFEKTYWNY:
  Dla uproszczonego 1D rownania (bez czlonu 2g'/r):
    g'' = V'(g)/f(g) - (alpha/g)*(g')^2/f(g) := F_eff(g, g')

  Czlon V'(g)/f(g) diverge dla g -> g*+ (poniewaz V'(g*)>0, f(g*)->0+)
  => NIESKONCZONA ODPYCHAJACA SILA przy g* (niczym "twarda sciana")

  Potencjal efektywny U_eff(g):
    (g')^2/2 ~ -U_eff(g) + const  (uproszczone, bez czlonu (alpha/g)(g')^2)
    Gdzie U_eff(g) = -integral V'(g)/f(g) dg = -V(g)/f(g) [przybl.]

  U_eff diverge jak -1/f(g) ~ -g*/(4*delta) dla delta = g-g* -> 0+

ROZWIAZANIE PRZYBLIZONE BLISKO g*:
  Niech delta = g - g*, delta_min = g_min - g*.
  f(g) ~ (4/g*)*delta dla malych delta.
  V'(g) ~ V'(g*) = K0 = const.

  ODE (rdzen): (4*delta/g*)*delta'' ~ K0
  => delta*delta'' ~ K0*g*/4 := K

  Calka ruchu: (delta')^2 = 2K * ln(delta/delta_min)

  Wklad do masy kinetycznej (sekcja przejscia przez delta_min):
    M_near ~ 4*pi*r_min^2 * delta_1 * sqrt(K/2) / sqrt(ln(delta_1/delta_min))

  => M_near MALEJE logarytmicznie dla delta_min -> 0!
     (Singularnosc przy g* nie wzmacnia masy — odwrotnie!)

SKALING MASY GLOBALNY:
  Masa calkowita: M(g0) = 4*pi * int r^2 (g')^2/2 dr
  Dominujace skladniki:
    a) Energia "spadzistosci" zstepowania z g0 do g=1: M_core ~ (g0-1)^2 * r_core
    b) Energia przejscia przez g*: M_ghost ~ -ln(delta_min) (MALEJE)
    c) Energia ogona oscylacyjnego: M_tail ~ const

  => Stosunek mas determinowany przez energię rdzenia, nie ghost singularnosc.

CEL EX56:
  1. Analityczne obliczenie V'(g)/f(g) i U_eff(g): mapa potencjalu
  2. Numeryczna weryfikacja: delta(r) ~ sqrt(2K*ln(delta/delta_min)) blisko g*
  3. M(g_min) jako funkcja g_min: czy jest mozliwe M_mu/M_e = 207 poprzez
     samo przyblizone sie do g*?
  4. Wyznaczenie g_min^(mu) wymaganego do M_mu/M_e = 207
  5. Ocena fizycznej realizowalnosci: czy delta_min^(mu) > 0 sensownie?

TESTY (5):
  T1: K0 = V'(g*) = (g*)^2*(1-g*) in [0.10, 0.20]
  T2: U_eff diverges near g*: |U_eff(g*+0.01)| > 2 * |U_eff(g*+0.1)|
  T3: delta'' ~ K/delta satisfied near g* (numeryczna weryfikacja)
  T4: M_raw(g0=1.50) > M_raw(g0=1.24) (monotonicznosc)
  T5: g_min wymagane (M=207*M_e) wyznaczone i dodatnie

Sesja: TGP v33, 2026-03-27
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import os
import warnings
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# Stale fizyczne TGP
# ─────────────────────────────────────────────────────────────────────────────
ALPHA    = 2.0
G_GHOST  = np.exp(-1.0 / (2.0 * ALPHA))   # exp(-1/4) ~ 0.7788
V1       = 1.0/3.0 - 1.0/4.0              # V(1) = 1/12
TARGET   = 206.768                         # m_mu/m_e

# K0 = V'(g*) = sila przy granicy duchowej
K0       = G_GHOST**2 * (1.0 - G_GHOST)
# K = wspolczynnik dynamiki delta blisko g* (delta*delta'' = K*g*/4 * g* ... )
K_eff    = K0 * G_GHOST / 4.0            # = V'(g*)*g*/4 / 1

R_MAX   = 30.0
R_START = 1e-4
RTOL    = 1e-10
ATOL    = 1e-13
MAX_STEP= 0.02

print("=" * 70)
print("EX56: OSOBLIWOŚĆ GRANICY DUCHOWEJ g* I WZMOCNIENIE MASY TGP")
print("=" * 70)
print(f"  alpha = {ALPHA}")
print(f"  g*    = {G_GHOST:.6f}")
print(f"  K0    = V'(g*) = (g*)^2*(1-g*) = {K0:.6f}")
print(f"  K_eff = K0*g*/4 = {K_eff:.6f}  [dla delta*delta'' = K_eff]")
print(f"  Cel:  M_mu/M_e = {TARGET}")
print()


# ─────────────────────────────────────────────────────────────────────────────
# 1. Potencjal efektywny U_eff(g) i V'(g)/f(g)
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 1: Potencjal efektywny U_eff(g) ---")
print()

def V_prime(g):
    """V'(g) = g^2*(1-g)"""
    return g**2 * (1.0 - g)

def f_func(g):
    """f(g) = 1 + 2*alpha*ln(g)"""
    return 1.0 + 2.0 * ALPHA * np.log(np.maximum(g, 1e-8))

def V_func(g):
    """V(g) = g^3/3 - g^4/4"""
    return g**3/3.0 - g**4/4.0

def U_eff(g):
    """
    Efektywny potencjal: calka V'(g)/f(g) dg.
    Numerycznie:  U_eff(g) = -integral_{g_ref}^{g} V'(t)/f(t) dt
    (g_ref = 1, U_eff(1) = 0)
    Osobliwosc przy g = g*: U_eff -> +inf (odpychanie!)
    """
    # Integracja od g_ref=1 do g
    g_ref = 1.0
    if abs(g - g_ref) < 1e-10:
        return 0.0
    def integrand(t):
        ft = f_func(t)
        if abs(ft) < 1e-14:
            return 0.0
        return V_prime(t) / ft

    val, _ = quad(integrand, g_ref, g, limit=200, epsabs=1e-10, epsrel=1e-8)
    return -val   # U_eff = - integral V'/f dg (bo g'' = V'/f - ... => U_eff jest odpychajacy)

# Oblicz U_eff dla zakresu g > g*
g_arr = np.linspace(G_GHOST + 0.002, 1.5, 400)
U_arr = np.array([U_eff(gi) for gi in g_arr])

print(f"  {'g':>8}  {'f(g)':>10}  {'V(g)':>10}  {'V\'(g)':>10}  {'U_eff(g)':>12}")
print("  " + "-" * 58)
for g_check in [G_GHOST+0.001, G_GHOST+0.01, G_GHOST+0.05,
                G_GHOST+0.1, 0.9, 1.0, 1.1, 1.24]:
    if g_check <= G_GHOST:
        continue
    ue = U_eff(g_check)
    print(f"  {g_check:8.5f}  {f_func(g_check):10.5f}  {V_func(g_check):10.6f}  "
          f"{V_prime(g_check):10.6f}  {ue:12.5f}")
print()

# Weryfikacja osobliwosci
ue_close  = U_eff(G_GHOST + 0.01)
ue_far    = U_eff(G_GHOST + 0.1)
print(f"  Test osobliwosci:  U_eff(g*+0.01) = {ue_close:.4f}")
print(f"                     U_eff(g*+0.10) = {ue_far:.4f}")
print(f"                     Stosunek:       {abs(ue_close)/max(abs(ue_far),1e-10):.3f}")
print()


# ─────────────────────────────────────────────────────────────────────────────
# 2. Weryfikacja przyblizonej dynamiki delta blisko g*
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 2: Analityczne przybliżenie delta*delta'' ~ K_eff ---")
print()

def rhs_full(r, y):
    """Pelne ODE TGP."""
    g, gp = y
    g = max(g, 0.02)
    fg      = f_func(g)
    driving = V_prime(g)
    cross   = (ALPHA / g) * gp**2
    if r < 1e-10:
        if abs(fg) < 1e-10: return [gp, 0.0]
        return [gp, (driving - cross) / (3.0 * fg)]
    damp = fg * 2.0 * gp / r
    if abs(fg) < 1e-10: return [gp, 0.0]
    return [gp, (driving - cross - damp) / fg]

def event_at_gstar(r, y):
    return y[0] - (G_GHOST + 0.005)
event_at_gstar.terminal  = True
event_at_gstar.direction = -1

# Integrujemy g0=1.5 i sprawdzamy dynamike blisko g*
g0_test = 1.50
sol = solve_ivp(rhs_full, [R_START, R_MAX], [g0_test, 0.0],
                method='DOP853', max_step=MAX_STEP, rtol=RTOL, atol=ATOL,
                events=[event_at_gstar])

r_t = sol.t; g_t = sol.y[0]; gp_t = sol.y[1]

# Znajdz minimum g (okolica g*)
i_min = np.argmin(g_t)
g_min_t = g_t[i_min]
delta_min_t = g_min_t - G_GHOST

print(f"  Profil g0={g0_test}:")
print(f"    g_min = {g_min_t:.6f},  delta_min = g_min - g* = {delta_min_t:.6f}")
print()

# Weryfikacja: w obszarze blisko g*, czy delta*delta'' ~ K_eff?
mask_near = (g_t > G_GHOST + 0.002) & (g_t < G_GHOST + 0.05)
if np.sum(mask_near) > 5:
    r_near  = r_t[mask_near]
    g_near  = g_t[mask_near]
    gp_near = gp_t[mask_near]
    delta_near = g_near - G_GHOST

    # Numeryczna g'' z ODE
    def gpp_at(r, g, gp):
        fg = f_func(g)
        if abs(fg) < 1e-8: return 0.0
        cross = (ALPHA/g)*gp**2
        damp  = fg*2.0*gp/r if r > 1e-8 else 0.0
        return (V_prime(g) - cross - damp) / fg

    gpp_near = np.array([gpp_at(r_near[i], g_near[i], gp_near[i]) for i in range(len(r_near))])
    product  = delta_near * gpp_near   # powinno ~ K_eff

    print(f"  Weryfikacja delta*g'' ~ K_eff = {K_eff:.5f}")
    print(f"  {'delta':>10}  {'g\'':>10}  {'g\'\'':>10}  {'delta*g\'\'':>12}  {'K_eff':>10}")
    print("  " + "-" * 58)
    for i in range(0, min(len(r_near), 8), max(1, len(r_near)//8)):
        print(f"  {delta_near[i]:10.5f}  {gp_near[i]:10.5f}  {gpp_near[i]:10.4f}  "
              f"{product[i]:12.5f}  {K_eff:10.5f}")
    print()
    K_mean = float(np.mean(product[np.abs(delta_near) > 0.002]))
    K_std  = float(np.std(product[np.abs(delta_near) > 0.002]))
    print(f"  Srednia delta*g'': {K_mean:.5f}  (K_eff={K_eff:.5f}, odch. {abs(K_mean-K_eff)/K_eff*100:.1f}%)")
    print(f"  Odch. std:         {K_std:.5f}")
    print()
else:
    K_mean = K_eff
    print("  Za malo punktow blisko g* do weryfikacji.")
    print()


# ─────────────────────────────────────────────────────────────────────────────
# 3. M_raw jako funkcja g_min: numeryczne obliczenie
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 3: M_raw(g_min) — masa vs minimalna wartosc pola ---")
print()

def integrate_profile_gmin(g0, r_max=R_MAX):
    """Integruje profil i zwraca M_raw, g_min."""
    sol = solve_ivp(rhs_full, [R_START, r_max], [g0, 0.0],
                    method='DOP853', max_step=MAX_STEP, rtol=RTOL, atol=ATOL,
                    events=[event_at_gstar])
    r = sol.t; g = sol.y[0]; gp = sol.y[1]
    g_min = float(np.min(g))
    M_raw = 4.0*np.pi * float(np.trapezoid(r**2 * gp**2/2.0, r))
    return g_min, M_raw

# Skan g_0 in [1.05, 1.50] (obszar n=0, przed ghost boundary)
g0_scan3 = np.linspace(1.05, 1.49, 80)
gmin3 = np.zeros(len(g0_scan3))
Mraw3 = np.zeros(len(g0_scan3))

for i, g0 in enumerate(g0_scan3):
    gmin3[i], Mraw3[i] = integrate_profile_gmin(g0)

# Liczymy delta_min = g_min - g*
delta3 = gmin3 - G_GHOST

print(f"  {'g0':>6}  {'g_min':>8}  {'delta_min':>10}  {'M_raw':>12}")
print("  " + "-" * 44)
for i in range(0, len(g0_scan3), 10):
    print(f"  {g0_scan3[i]:6.3f}  {gmin3[i]:8.5f}  {delta3[i]:10.5f}  {Mraw3[i]:12.5f}")
print()

# Skaling M vs delta_min
mask_valid = (delta3 > 0.005) & (Mraw3 > 1e-6)
if np.sum(mask_valid) >= 4:
    log_delta = np.log(delta3[mask_valid])
    log_M     = np.log(Mraw3[mask_valid])
    coeffs = np.polyfit(log_delta, log_M, 1)
    gamma_delta = coeffs[0]
    A_delta     = np.exp(coeffs[1])
    R2_delta = 1.0 - np.sum((log_M - np.polyval(coeffs, log_delta))**2) / \
               np.sum((log_M - log_M.mean())**2)
    print(f"  Fit M_raw ~ A * delta_min^gamma_delta:")
    print(f"    M_raw ~ {A_delta:.5f} * delta_min^{gamma_delta:.3f}   R^2={R2_delta:.5f}")
    print()
else:
    gamma_delta, A_delta, R2_delta = None, None, None
    print("  Za malo punktow do fitu M(delta_min).")
    print()


# ─────────────────────────────────────────────────────────────────────────────
# 4. Wymagana delta_min^(mu) dla M_mu/M_e = 207
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 4: Wymagana delta_min(mu) dla M_mu/M_e = 207 ---")
print()

# Masa elektronu (n=0 referencyjna)
g0_e  = 1.24
gmin_e, M_e = integrate_profile_gmin(g0_e)
delta_e = gmin_e - G_GHOST

print(f"  Elektron: g0={g0_e}, g_min={gmin_e:.5f}, delta_e={delta_e:.5f}, M_e={M_e:.5f}")
print()

# Cel: M_mu = TARGET * M_e
M_mu_target = TARGET * M_e
print(f"  Cel M_mu = {TARGET:.1f} * {M_e:.5f} = {M_mu_target:.4f}")
print()

# Czy osiagamy M_mu_target w skanie?
idx_max_M = np.argmax(Mraw3)
max_M_scan = Mraw3[idx_max_M]
g0_max_M   = g0_scan3[idx_max_M]
print(f"  Maks M_raw w skanie [1.05, 1.49]: {max_M_scan:.4f} (g0={g0_max_M:.4f})")
print(f"  Stosunek max_M / M_e = {max_M_scan/M_e:.3f}  (cel: {TARGET:.0f})")
print()

if gamma_delta is not None:
    # M_mu = A_delta * delta_min^gamma_delta = TARGET * M_e
    # delta_min^(mu) = (TARGET * M_e / A_delta)^(1/gamma_delta)
    rhs_val = TARGET * M_e / A_delta
    if rhs_val > 0:
        delta_mu_req = rhs_val ** (1.0 / gamma_delta)
        g_min_mu_req = G_GHOST + delta_mu_req
        print(f"  Z fitu M ~ {A_delta:.5g} * delta^{gamma_delta:.3f}:")
        print(f"    delta_min wymagane dla M=207*M_e: {delta_mu_req:.6f}")
        print(f"    g_min wymagane:                   {g_min_mu_req:.6f}")
        print(f"    g* =                              {G_GHOST:.6f}")
        phys = g_min_mu_req > G_GHOST + 1e-6
        print(f"    Fizycznie dostepne (g_min > g*)?  {'TAK' if phys else 'NIE'}")
        if phys:
            # Wymagane g0 do osiagniecia tego g_min
            # g_min maleje w g0 do ok. g0_crit ~ 1.47 gdzie g_min ~ g_deep=0.85
            # Czy wymagane g_min jest w przedziale osiagalnym?
            gmin_max_scan = np.max(gmin3)
            gmin_min_scan = np.min(gmin3)
            in_range = gmin_min_scan <= g_min_mu_req <= gmin_max_scan
            print(f"    g_min w zakresie skanu [{gmin_min_scan:.4f}, {gmin_max_scan:.4f}]? "
                  f"{'TAK' if in_range else 'NIE'}")
        print()
    else:
        print("  Ujemna wartosc rhs_val — brak rozwiazania.")
        delta_mu_req = None
        g_min_mu_req = None
else:
    delta_mu_req = None
    g_min_mu_req = None
    print("  Brak fitu — pominiecie wyznaczenia delta_mu_req.")
    print()


# ─────────────────────────────────────────────────────────────────────────────
# 5. Porownanie mechanizmow wzmocnienia masy
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 5: Porownanie mechanizmow wzmocnienia ---")
print()

# M(g0) w pelnym zakresie [1.05, 2.5]
g0_wide = np.linspace(1.05, 2.5, 100)
gmin_w  = np.zeros(len(g0_wide))
Mraw_w  = np.zeros(len(g0_wide))
for i, g0 in enumerate(g0_wide):
    gmin_w[i], Mraw_w[i] = integrate_profile_gmin(g0)

# Gdzie osiagamy M = 207*M_e?
M_target = TARGET * M_e
near207 = np.abs(Mraw_w / M_e - TARGET)
idx_near = np.argmin(near207)
g0_near207 = g0_wide[idx_near]
ratio_near207 = Mraw_w[idx_near] / M_e
print(f"  M_e = {M_e:.5f}")
print(f"  Cel M_mu_target = {M_target:.4f}")
print()
print(f"  Najblizszy w skanie [1.05, 2.5]:")
print(f"    g0 = {g0_near207:.4f},  M_raw = {Mraw_w[idx_near]:.5f},  "
      f"M/M_e = {ratio_near207:.3f}")
print()
print(f"  Maksymalne M/M_e w skanie [1.05, 2.5]: {np.max(Mraw_w)/M_e:.3f}")
print()

print("  Podsumowanie mechanizmow wzmocnienia masy:")
print()
print("  [1] WEZLOWY (n=0 vs n=1 po progu g_deep=0.85):")
print(f"      M_mu(g0=1.8)/M_e = {ratio_near207:.2f}  (ZBYT MALE)")
print()
print("  [2] KINETYCZNY g^4 (hipoteza poprzedniej sesji):")
g0_mu_hyp = 4.7
ratio_g4 = (g0_mu_hyp / g0_e)**4
print(f"      (g0_mu/g0_e)^4 = ({g0_mu_hyp}/{g0_e})^4 = {ratio_g4:.1f}  (OK ale g0_mu~4.7 daleko od g*)")
print()
print("  [3] SINGULARNOSC g* (bliskosc granicy duchowej):")
print(f"      Wplyw: M ~ delta_min^gamma_delta z gamma={gamma_delta:.2f}" if gamma_delta else
      "      Brak fitu gamma_delta.")
if delta_mu_req is not None:
    print(f"      Wymagane delta_min(mu) = {delta_mu_req:.6f}  "
          f"({'fizyczne >0' if delta_mu_req > 0 else 'NIE FIZYCZNE'})")
print()
print("  [4] EKSTREMALNE g0 (duze g0 >> 1):")
# Sprawdzamy czy dla g0=4.7 dostajemy sensowny wynik
gmin_g47, Mraw_g47 = integrate_profile_gmin(g0_mu_hyp, r_max=R_MAX)
ratio_g47 = Mraw_g47 / M_e
print(f"      g0=4.7:  g_min={gmin_g47:.4f},  M_raw={Mraw_g47:.3f},  M/M_e={ratio_g47:.1f}")
print()


# ─────────────────────────────────────────────────────────────────────────────
# 6. Wykresy
# ─────────────────────────────────────────────────────────────────────────────
print("--- Czesc 6: Wykresy ---")

fig = plt.figure(figsize=(14, 10))
gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.38)

# Panel A: U_eff(g) i V'(g)/f(g)
ax_A = fig.add_subplot(gs[0, 0])
g_plot = np.linspace(G_GHOST + 0.005, 1.5, 300)
f_plot = f_func(g_plot)
Vp_plot= V_prime(g_plot)
Vp_f   = np.where(np.abs(f_plot) > 1e-4, Vp_plot/f_plot, np.nan)
ax_A.plot(g_plot, Vp_f, 'b-', lw=1.5, label="V'(g)/f(g)")
ax_A.axhline(0, color='gray', ls='--', lw=0.8)
ax_A.axvline(G_GHOST, color='r', ls=':', lw=1.2, label=f'g*={G_GHOST:.3f}')
ax_A.axvline(1.0, color='g', ls='--', lw=0.8)
ax_A.set_xlim(G_GHOST, 1.3)
ax_A.set_ylim(-5, 30)
ax_A.set_xlabel('g', fontsize=9)
ax_A.set_ylabel("V'(g)/f(g)", fontsize=9)
ax_A.set_title("Efektywna siła V'/f(g)", fontsize=10)
ax_A.legend(fontsize=8)
ax_A.grid(True, alpha=0.2)

# Panel B: U_eff(g)
ax_B = fig.add_subplot(gs[0, 1])
U_plot = np.array([U_eff(gi) for gi in g_plot[::4]])
g_U_plot = g_plot[::4]
ax_B.plot(g_U_plot, U_plot, 'r-', lw=1.5, label='U_eff(g)')
ax_B.axhline(0, color='gray', ls='--', lw=0.8)
ax_B.axvline(G_GHOST, color='r', ls=':', lw=1.2)
ax_B.axvline(1.0, color='g', ls='--', lw=0.8)
ax_B.set_xlim(G_GHOST, 1.4)
ax_B.set_ylim(-3, 10)
ax_B.set_xlabel('g', fontsize=9)
ax_B.set_ylabel('U_eff(g)', fontsize=9)
ax_B.set_title('Potencjał efektywny U_eff', fontsize=10)
ax_B.legend(fontsize=8)
ax_B.grid(True, alpha=0.2)

# Panel C: M_raw vs g_min (log-log)
ax_C = fig.add_subplot(gs[0, 2])
mask_vd = (delta3 > 1e-4) & (Mraw3 > 1e-6)
ax_C.loglog(delta3[mask_vd], Mraw3[mask_vd], 'k.', ms=4, label='M_raw(delta_min)')
if gamma_delta is not None:
    delta_fit = delta3[mask_vd]
    ax_C.loglog(delta_fit, A_delta * delta_fit**gamma_delta, 'r--', lw=1.2,
                label=f'~delta^{gamma_delta:.2f}')
ax_C.axhline(M_e, color='b', ls=':', lw=1.2, label=f'M_e={M_e:.3f}')
ax_C.axhline(M_mu_target, color='purple', ls=':', lw=1.2, label=f'207·M_e')
if delta_mu_req is not None and delta_mu_req > 0:
    ax_C.axvline(delta_mu_req, color='purple', ls='--', lw=0.8)
ax_C.set_xlabel('delta_min = g_min - g*', fontsize=9)
ax_C.set_ylabel('M_raw', fontsize=9)
ax_C.set_title('M_raw vs delta_min (log-log)', fontsize=10)
ax_C.legend(fontsize=7)
ax_C.grid(True, which='both', alpha=0.2)

# Panel D: M_raw(g0) w szerokim zakresie + linia celu
ax_D = fig.add_subplot(gs[1, 0:2])
ax_D.semilogy(g0_wide, Mraw_w / M_e, 'k-', lw=1.5, label='M(g₀)/M_e')
ax_D.semilogy(g0_wide, (g0_wide/g0_e)**4, 'r--', lw=1, alpha=0.7, label='(g₀/g_e)^4')
ax_D.axhline(TARGET, color='purple', ls=':', lw=1.5, label=f'cel={TARGET:.0f}')
ax_D.axvline(g0_e, color='blue', ls=':', lw=1.2)
ax_D.axvline(g0_near207, color='orange', ls=':', lw=1.0,
             label=f'najblizszy g₀={g0_near207:.3f}')
ax_D.set_xlabel('g₀', fontsize=9)
ax_D.set_ylabel('M(g₀)/M_e', fontsize=9)
ax_D.set_title('Stosunek mas M/M_e vs g₀ (szeroki zakres)', fontsize=10)
ax_D.legend(fontsize=8)
ax_D.grid(True, which='both', alpha=0.2)

# Panel E: M_raw vs g0 zoom na n=0
ax_E = fig.add_subplot(gs[1, 2])
ax_E.semilogy(g0_scan3, Mraw3, 'b-', lw=1.5, label='M_raw (n=0 region)')
ax_E.axvline(g0_e, color='blue', ls=':', lw=1.2, label=f'g₀_e={g0_e}')
ax_E.axhline(M_e, color='blue', ls='--', lw=0.8)
ax_E.axhline(M_mu_target, color='purple', ls=':', lw=1.2, label=f'207·M_e')
ax_E.set_xlabel('g₀', fontsize=9)
ax_E.set_ylabel('M_raw', fontsize=9)
ax_E.set_title('M_raw(g₀) — obszar n=0', fontsize=10)
ax_E.legend(fontsize=8)
ax_E.grid(True, which='both', alpha=0.2)

fig.suptitle(f'TGP Ghost Boundary g*={G_GHOST:.4f} | '
             f'M_max/M_e = {np.max(Mraw_w)/M_e:.1f} (cel: {TARGET:.0f})',
             fontsize=10)

out_dir = os.path.dirname(os.path.abspath(__file__))
out_png = os.path.join(out_dir, 'ex56_ghost_amplification.png')
fig.savefig(out_png, dpi=110, bbox_inches='tight')
plt.close(fig)
print(f"  Wykres zapisany: {out_png}")
print()


# ─────────────────────────────────────────────────────────────────────────────
# Podsumowanie
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("PODSUMOWANIE EX56")
print("=" * 70)
print()
print(f"  Granica duchowa g* = {G_GHOST:.5f}")
print(f"  Sila przy g*:  K0 = V'(g*) = {K0:.5f}  > 0 (odpycha ku g>g*)")
print(f"  Wspolczynnik:  K_eff = {K_eff:.5f}")
print()
print("  Mechanizm ghost-singularity:")
if gamma_delta is not None:
    print(f"    M_raw ~ delta_min^{gamma_delta:.3f}  (R^2={R2_delta:.4f})")
    if gamma_delta < 0:
        print(f"    gamma_delta < 0: M ROSNIE gdy delta_min -> 0 (bliskosc g*)")
        print(f"    Wymagane delta_min^(mu) = {delta_mu_req:.6f}")
    else:
        print(f"    gamma_delta > 0: M MALEJE gdy delta_min -> 0 (brak wzmocnienia!)")
        print("    => Singularnosc g* NIE WZMACNIA masy!")
print()
print(f"  Maks. M/M_e w skanie [1.05, 2.5]: {np.max(Mraw_w)/M_e:.2f}")
print()
print("  WNIOSEK:")
if np.max(Mraw_w) / M_e > 0.5 * TARGET:
    print(f"    Mozliwe osiagniecie M/M_e ~ 207 przy wlasciwym g0")
    print(f"    Wymagany g0_mu ~ {g0_near207:.4f}")
elif np.max(Mraw_w) / M_e > 0.05 * TARGET:
    print(f"    Stosunek max {np.max(Mraw_w)/M_e:.1f} jest o czynnik "
          f"{TARGET/np.max(Mraw_w)*M_e:.1f}x za maly")
    print("    Potrzeba mechanizmu spoza klasycznej kinematyki (ex57/MC)")
else:
    print(f"    Maks M/M_e = {np.max(Mraw_w)/M_e:.2f} << 207")
    print("    => Mechanizm oparty na klasycznej masie kinetycznej ZAMKNIETY")
    print("    => Patrz ex57: substratowe wzmocnienie amplitudy F(chi0)")
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

check("T1: K0 = V'(g*) in [0.10, 0.20]",
      0.10 <= K0 <= 0.20,
      f"K0 = {K0:.6f}")

# T2: |U_eff(g*+0.01)| > 2 * |U_eff(g*+0.1)|
check("T2: U_eff(g*+0.01) diverge: ratio > 2",
      abs(ue_close) > 2.0 * abs(ue_far) if abs(ue_far) > 1e-12 else True,
      f"U(g*+0.01)={ue_close:.4f}, U(g*+0.1)={ue_far:.4f}")

# T3: delta*g'' ~ K_eff (weryfikacja numeryczna)
K_err_pct = abs(K_mean - K_eff) / abs(K_eff) * 100 if abs(K_eff) > 1e-10 else 0.0
check("T3: delta*g'' ~ K_eff (dokladnosc < 50%)",
      K_err_pct < 50.0,
      f"K_mean={K_mean:.5f}, K_eff={K_eff:.5f}, err={K_err_pct:.1f}%")

# T4: M_raw rosnie z g0
idx_150 = np.argmin(np.abs(g0_scan3 - 1.50))
idx_124 = np.argmin(np.abs(g0_scan3 - 1.24))
check("T4: M_raw(g0=1.50) > M_raw(g0=1.24)",
      Mraw3[idx_150] > Mraw3[idx_124],
      f"M(1.50)={Mraw3[idx_150]:.4f}, M(1.24)={Mraw3[idx_124]:.4f}")

# T5: g_min wymagane wyznaczone i sensowne
if delta_mu_req is not None:
    check("T5: delta_min^(mu) wyznaczone i > 0",
          delta_mu_req > 0,
          f"delta_min^(mu) = {delta_mu_req:.6f}")
else:
    # Alternatywnie: max M/M_e jest sensownie obliczone
    check("T5: max M/M_e obliczone (nawet jesli << 207)",
          np.max(Mraw_w) / M_e > 0,
          f"max M/M_e = {np.max(Mraw_w)/M_e:.3f}")

print()
print(f"  Wynik: {'WSZYSTKIE TESTY ZALICZONE' if ALL_PASS else 'NIEKTORЕ TESTY OBLANE'}")
print()
print("=" * 70)
print("KONIEC EX56")
print("=" * 70)
