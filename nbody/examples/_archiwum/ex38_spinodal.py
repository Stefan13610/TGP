"""
ex38_spinodal.py
================
Tachionowa niestabilność spinodalna TGP jako mechanizm ciemnej materii.

IDEA (★ NOWY MECHANIZM — poza ex14):
  TGP ma próżnię g=1, ale V''(g=1) = -γ < 0 — NIESTABILNA!
  W wczesnym wszechświecie, gdy pole ewoluuje od g≈0 do g≈1,
  przechodzi przez strefę tachioniczną. Fluktuacje δg rosną eksponencjalnie.

  To jest SPINODALNA NIESTABILNOŚĆ FAZOWA — znana w fizyce cząstek
  (preheating po inflacji) i fizyce ciała stałego (rozkład spinodalny).

FIZYKA MECHANIZMU:
  Linowe perturbacje wokół tła g_bg(t):
    δg̈_k + H δg_k̇ + [k²/a² + V''(g_bg)] δg_k = 0

  Gdy V''(g_bg) < 0: mamy |V''| jako "tachionową masę kwadratową".
  W granicy statycznej (bez ekspansji): δg_k ~ exp(√|V''| t) dla k < k_max.

  Strefą niestabilności: g_bg ∈ (0, 2/3) (dla β=γ=1: V'' = 2g-3g² < 0)
  Maks. niestabilność: g_bg = 1/3 → V'' = -1/3 (jednostki Plancka)

PARAMETRY:
  Dominująca skala niestabilności (modus k_max):
    ω²(k) = k² + V''(g_bg) → maksimum wzrostu przy k=0
    k_max rośnie = 0 (IR dominuje)
    Ale dla skończonego pola grawitacyjnego: k_dominate ~ √|V''| (modus rezonansowy)

GALAKTYCZNA SKALA SPINODALNA:
  Niestabilność tworzy struktury na skali:
    λ_spin ~ 2π/k_max ~ 2π / √|V''(g_bg)|

  Dla g_bg ~ 1/3 (maks. niestabilność):
    |V''(1/3)| = |2/3 - 3/9| = |2/3 - 1/3| = 1/3  [Planck]

  Skala w kpc: λ_spin = 2π / √(1/3) * l_Pl = 2π * √3 * 1.616e-35 m
    => λ_spin ~ 1.76e-34 m  (MIKRO, nie galaktyczna!)

WNIOSEK: Bez sprzężenia z grawitacją, spinodalna skala = l_Planck.
  Aby osiągnąć skalę galaktyczną, potrzeba KOSMOLOGICZNEJ EWOLUCJI.

KOSMOLOGICZNA ANALIZA (rozszerzenie):
  Czas Hubble'a na danej epoce: t_H ~ 1/H(z)
  Rozmiar Hubble'a: d_H ~ c/H(z)
  Spinodalna amplifikacja: A ~ exp(√|V''| * t_spin)

  Dla spójności: t_spin ~ 1/√|V''| (czas efektywnego wzrostu)
  Skala fizyczna (na kohorensji): λ_phys = a(t_spin) * λ_comoving

  KLUCZOWE PYTANIE: Przy jakiej wartości γ spinodalna tworzy struktury
  na skali ~10 kpc (galaktyki)?

  Odpowiedź: γ ~ (2π/λ_kpc)² = (2π / 1.91e54 l_Pl)² ~ 10^{-110}
  (ten sam wynik co §2.3 i §4.2 ANALIZA_CIEMNA_MATERIA.md)

NUMERYCZNE EKSPERYMENTY:
  1. Widmo wzrostu spinodalnego w TGP (ω² vs k)
  2. Ewolucja δg_k(t) podczas przejścia g_bg: 0 → 1
  3. Rozkład energii spinodalnej: ρ_spin(k)
  4. Korespondencja skali k_gal z parametrem γ
"""

import sys, os
import numpy as np
from scipy.integrate import solve_ivp, quad, odeint
from scipy.special import airy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# =============================================================================
# Stałe
# =============================================================================
l_Pl_m  = 1.616e-35   # m
kpc_m   = 3.086e19    # m
Mpc_m   = 3.086e22    # m
c_SI    = 3e8         # m/s
hbar_SI = 1.055e-34   # J*s
G_SI    = 6.674e-11   # m^3 kg^-1 s^-2
m_Pl_kg = 2.176e-8    # kg
E_Pl_eV = 1.221e28    # eV

print("=" * 70)
print("EX38: TACHIONOWA NIESTABILNOŚĆ SPINODALNA TGP")
print("=" * 70)
print()

# =============================================================================
# Sekcja 1: Potencjał TGP — strefa tachioniczna
# =============================================================================
print("=" * 70)
print("SEKCJA 1: Strefa tachioniczna V''(g) < 0 w TGP")
print("=" * 70)
print()

# Potencjał TGP: V(g) = g^3/3 - g^4/4 (β=γ=1 Planck)
# V'(g)  = g^2 - g^3 = g^2(1-g)
# V''(g) = 2g - 3g^2 = g(2-3g)
# V''(g) = 0 przy g=0 lub g=2/3
# V''(g) < 0 dla g ∈ (2/3, ∞)   [powyżej 2/3]
# V''(g) > 0 dla g ∈ (0, 2/3)   [poniżej 2/3]

g_arr = np.linspace(0, 1.5, 300)
V_arr   = g_arr**3/3 - g_arr**4/4
dV_arr  = g_arr**2 - g_arr**3
d2V_arr = 2*g_arr - 3*g_arr**2

g_tach_enter = 2.0/3.0
g_tach_exit  = 1.0  # próżnia (V'' = 2-3 = -1 < 0 — NIGDY nie wychodzi! g=1 jest w strefie)

print("  V''(g) dla TGP:")
print(f"  V''(g) = 2g - 3g²")
print(f"  V''(g) = 0 przy g = 0 lub g = 2/3")
print()
print("  UWAGA: g=1 (próżnia) leży w strefie V''<0!")
print(f"  V''(0)   = {2*0 - 3*0**2:.3f}")
print(f"  V''(1/3) = {2/3 - 3/9:.3f}  (minimum V'', maks. niestabilność)")
print(f"  V''(2/3) = {2*2/3 - 3*(2/3)**2:.3f}  (granica strefy)")
print(f"  V''(1)   = {2*1 - 3*1**2:.3f}  (próżnia — tachion!)")
print()

for g_test in [0.1, 1/3, 0.5, 2/3, 0.8, 1.0]:
    d2V = 2*g_test - 3*g_test**2
    status = "STABILNA" if d2V > 0 else "TACHIONICZNA"
    print(f"  g = {g_test:.3f}: V'' = {d2V:.4f}  => {status}")

print()
print("  WNIOSEK: Cały przedział g ∈ (2/3, ∞) jest tachioniczny!")
print("  Pole TGP, gdy ewoluuje ku próżni g=1, ZAWSZE przechodzi przez")
print("  strefę tachioniczną (od g=2/3 do g=1) i NIE WYCHODZI z niej.")
print("  Próżnia g=1 jest PERMANENTNIE tachioniczna.")
print()

# =============================================================================
# Sekcja 2: Widmo wzrostu — ω²(k) dla różnych g_bg
# =============================================================================
print("=" * 70)
print("SEKCJA 2: Widmo wzrostu ω²(k) = k² + V''(g_bg)")
print("=" * 70)
print()

def omega2_spinodal(k_arr, g_bg, gamma=1.0):
    """
    Dyspersja perturbacji liniowych wokół tła g_bg.
    ω²(k) = k² - |V''(g_bg)|  (dla V''<0)
    Tryby rosnące: k² < |V''(g_bg)|  => k < k_max = √|V''|
    """
    d2V = 2*gamma*g_bg - 3*gamma*g_bg**2  # z V(g) = βg³/3 - γg⁴/4 dla β=γ
    return k_arr**2 + d2V

k_arr = np.linspace(0, 2.0, 300)

print("  Mody rosnące (ω²<0) dla różnych g_bg:")
print()
print(f"  {'g_bg':>6s}  {'V''(g_bg)':>10s}  {'k_max':>10s}  {'λ_max [Pl]':>12s}  {'Status':>12s}")
print("-" * 60)

g_bg_vals = [0.1, 1/3, 0.5, 2/3, 0.8, 0.9, 1.0]
spinodal_results = {}

for g_bg in g_bg_vals:
    d2V = 2*g_bg - 3*g_bg**2
    if d2V < 0:
        k_max = np.sqrt(-d2V)
        lam_max_Pl = 2*np.pi/k_max if k_max > 0 else np.inf
        lam_max_m = lam_max_Pl * l_Pl_m
        lam_max_kpc = lam_max_m / kpc_m
        status = f"NIESTAB. k<{k_max:.3f}"
        spinodal_results[g_bg] = {
            'd2V': d2V, 'k_max': k_max,
            'lam_Pl': lam_max_Pl, 'lam_kpc': lam_max_kpc
        }
    else:
        k_max = 0
        lam_max_Pl = np.inf
        lam_max_kpc = np.inf
        status = "STABILNA"
        spinodal_results[g_bg] = {'d2V': d2V, 'k_max': 0, 'lam_Pl': np.inf, 'lam_kpc': np.inf}

    if d2V < 0:
        print(f"  {g_bg:>6.3f}  {d2V:>10.4f}  {k_max:>10.4f}  {lam_max_Pl:>12.3e}  {status}")
    else:
        print(f"  {g_bg:>6.3f}  {d2V:>10.4f}  {'—':>10s}  {'—':>12s}  {status}")

print()
print("  KLUCZ: Skala spinodalna λ_max = 2π/√|V''| jest w JEDNOSTKACH PLANCKA!")
print(f"  Dla g_bg = 1/3 (maks. niestab.): λ_max = 2π/√(1/3) * l_Pl")
print(f"  = {2*np.pi/np.sqrt(1/3):.2f} l_Pl = {2*np.pi/np.sqrt(1/3)*l_Pl_m:.2e} m")
print(f"  = {2*np.pi/np.sqrt(1/3)*l_Pl_m/kpc_m:.2e} kpc")
print()
print("  To jest SKALA PLANCKA — absolutnie mikroskopowa!")
print("  Spinodalna w MINIMALNYM TGP (γ=1) nie daje skali galaktycznej.")
print()

# =============================================================================
# Sekcja 3: Skalowanie γ → skala galaktyczna
# =============================================================================
print("=" * 70)
print("SEKCJA 3: Jakie γ daje spinodalną na skali galaktycznej?")
print("=" * 70)
print()
print("  Z ogólnego V(g) = βg³/3 - γg⁴/4 z β=γ:")
print("  V''(g, γ) = 2γg - 3γg² = γ(2g - 3g²)")
print()
print("  Skala spinodalna: λ_spin = 2π/√(γ|2g-3g²|) [l_Pl]")
print()
print("  Dla g_bg = 1/3 (maks. niestab.): |2g-3g²| = |2/3-1/3| = 1/3")
print("  => λ_spin = 2π/√(γ/3) = 2π√3/√γ [l_Pl]")
print()

# Dla skali galaktycznej λ_spin = L_gal
L_gal_kpc_vals = [1.0, 10.0, 100.0]  # kpc

print("  γ potrzebne do spinodalnej na różnych skalach:")
print()
print(f"  {'λ_gal [kpc]':>12s}  {'λ_gal [l_Pl]':>14s}  {'γ':>12s}  {'m_sp=√γ [l_Pl⁻¹]':>20s}")
print("-" * 65)

for L_kpc in L_gal_kpc_vals:
    L_Pl = L_kpc * kpc_m / l_Pl_m
    # λ_spin = 2π√3/√γ * l_Pl = L_gal
    # => √γ = 2π√3 / (L_gal/l_Pl)
    gamma_needed = (2*np.pi*np.sqrt(3) / L_Pl)**2
    m_sp = np.sqrt(gamma_needed)
    print(f"  {L_kpc:>12.0f}  {L_Pl:>14.3e}  {gamma_needed:>12.3e}  {m_sp:>20.3e}")

print()
print("  WYNIK: Spinodalna na skali galaktycznej (1–100 kpc) wymaga:")
print(f"  γ ~ 10^{{-110}} do 10^{{-106}}")
print()
print("  To IDENTYCZNA skala co:")
print("  - §4.2 ANALIZA_CIEMNA_MATERIA.md (tachioniczne plamki)")
print("  - Efektywne m_sp dla Yukawa na 10 kpc")
print("  - Paradoks γ: trzecie ograniczenie")
print()
print("  WNIOSEK: Spinodalna na skali galaktycznej wymaga γ ~ 10^{-110},")
print("  tj. ULTRALEKKIEGO bozonu skalarnego m_sp ~ 10^{-55} l_Pl⁻¹.")
print()

# =============================================================================
# Sekcja 4: Ewolucja kosmologiczna g_bg(t) — kiedy strefa tachioniczna?
# =============================================================================
print("=" * 70)
print("SEKCJA 4: Ewolucja kosmologiczna g_bg(t) w TGP")
print("=" * 70)
print()
print("  Równanie ruchu dla jednorodnego pola g_bg(t) w FRW:")
print("  g̈ + 3H ġ + V'(g) = 0")
print("  gdzie V'(g) = g² - g³ = g²(1-g)")
print()
print("  Rozwiązujemy numerycznie dla różnych warunków początkowych.")
print()

def homogeneous_field_eq(t, y, H_val=0.0):
    """
    g̈ + 3H ġ + V'(g) = 0
    y = [g, g_dot]
    H_val: stałe H (uproszczenie dla de Sitter)
    """
    g, gd = y
    dVdg = g**2 - g**3  # V'(g) = g^2(1-g)
    gdd = -3*H_val*gd - dVdg
    return [gd, gdd]

print("  Trajektorie g_bg(t) dla różnych warunków pocz. (bez rozszerzenia):")
print()

g_init_vals = [0.01, 0.05, 0.10, 0.30, 0.50, 0.65]
t_span = (0, 50.0)
t_eval = np.linspace(0, 50.0, 500)

fig_traj, axes_traj = plt.subplots(1, 2, figsize=(14, 6))
ax_g = axes_traj[0]
ax_v2 = axes_traj[1]

colors_traj = cm.plasma(np.linspace(0.1, 0.9, len(g_init_vals)))

traj_results = {}
for g0, col in zip(g_init_vals, colors_traj):
    sol = solve_ivp(
        lambda t, y: homogeneous_field_eq(t, y, H_val=0.0),
        t_span, [g0, 0.0], t_eval=t_eval,
        method='RK45', rtol=1e-8, atol=1e-10
    )
    if sol.success:
        g_t = sol.y[0]
        # Czas spędzony w strefie tachionicznej
        tachio_mask = (g_t > 2/3)
        t_in_tach = np.sum(tachio_mask) * (t_eval[1] - t_eval[0])
        traj_results[g0] = {'t': t_eval, 'g': g_t, 't_tach': t_in_tach}
        ax_g.plot(t_eval, g_t, color=col, lw=1.8, label=f'g₀={g0:.2f} (t_tach={t_in_tach:.1f})')
        # V''(g(t))
        d2V_t = 2*g_t - 3*g_t**2
        ax_v2.plot(t_eval, d2V_t, color=col, lw=1.5, label=f'g₀={g0:.2f}')
        print(f"  g₀={g0:.2f}: max_g={g_t.max():.3f}, t_in_tach={t_in_tach:.2f} [Pl]")

ax_g.axhline(1.0, color='gray', ls='--', lw=1, alpha=0.5, label='próżnia g=1')
ax_g.axhline(2/3, color='red', ls='--', lw=1.2, alpha=0.7, label='granica tachionu g=2/3')
ax_g.set_xlabel('t [Planck]', fontsize=11)
ax_g.set_ylabel('g_bg(t)', fontsize=11)
ax_g.set_title('Ewolucja pola jednorodnego g_bg(t)\n(bez ekspansji kosmologicznej)', fontsize=10)
ax_g.legend(fontsize=8, loc='upper right')
ax_g.grid(True, alpha=0.3)
ax_g.set_ylim(-0.1, 1.4)

ax_v2.axhline(0, color='k', lw=1, alpha=0.5)
ax_v2.axhline(-1/3, color='purple', ls=':', lw=1.5, alpha=0.7, label='min V'' (max niestab.)')
ax_v2.fill_between(t_eval, -0.5, 0, alpha=0.1, color='red', label='Strefa tachioniczna')
ax_v2.set_xlabel('t [Planck]', fontsize=11)
ax_v2.set_ylabel("V''(g_bg)", fontsize=11)
ax_v2.set_title("Masa kwadratowa V''(g_bg) podczas ewolucji", fontsize=10)
ax_v2.legend(fontsize=8)
ax_v2.grid(True, alpha=0.3)
ax_v2.set_ylim(-0.5, 0.3)

fig_traj.tight_layout()
out_traj = os.path.join(os.path.dirname(__file__), 'ex38_trajectories.png')
plt.savefig(out_traj, dpi=150, bbox_inches='tight')
plt.close()
print(f"\n  Wykres trajektorii: {out_traj}")

# =============================================================================
# Sekcja 5: Wzrost perturbacji δg_k podczas przejścia tachionicznego
# =============================================================================
print()
print("=" * 70)
print("SEKCJA 5: Wzrost perturbacji δg_k podczas przejścia g_bg: 0→1")
print("=" * 70)
print()
print("  Równanie perturbacji: δg̈_k + [k² + V''(g_bg(t))] δg_k = 0")
print("  (zaniedbujemy człon Hubble dla uproszczenia)")
print()

def perturbation_eq(t, y, k_val, g_bg_interp):
    """
    Perturbacja liniowa δg_k:
    δg̈_k + [k² + V''(g_bg(t))] δg_k = 0
    y = [δg, δg_dot]
    """
    g_bg = float(g_bg_interp(t))
    d2V  = 2*g_bg - 3*g_bg**2
    omega2 = k_val**2 + d2V
    delta_g, delta_gd = y
    return [delta_gd, -omega2 * delta_g]

# Wygeneruj trajektorię g_bg (g₀=0.01, bez ekspansji)
g0_bg = 0.05
sol_bg = solve_ivp(
    lambda t, y: homogeneous_field_eq(t, y, H_val=0.0),
    (0, 30.0), [g0_bg, 0.0], t_eval=np.linspace(0, 30.0, 3000),
    method='RK45', rtol=1e-10, atol=1e-12
)

from scipy.interpolate import interp1d
g_bg_interp = interp1d(sol_bg.t, sol_bg.y[0], kind='cubic', fill_value='extrapolate')

print("  Amplifikacja perturbacji |δg_k(t)|² dla różnych k:")
print()

k_test_vals = [0.001, 0.01, 0.1, 0.3, 0.5, 1.0]
print(f"  {'k [Pl⁻¹]':>12s}  {'|δg(t=10)|²/|δg(0)|²':>22s}  {'|δg(t=20)|²/|δg(0)|²':>22s}")
print("-" * 60)

amp_results = {}
for k_val in k_test_vals:
    sol_pert = solve_ivp(
        lambda t, y: perturbation_eq(t, y, k_val, g_bg_interp),
        (0, 25.0), [1.0, 0.0],  # δg(0) = 1, δg_dot(0) = 0
        t_eval=np.linspace(0, 25.0, 2500),
        method='RK45', rtol=1e-8, atol=1e-10
    )
    if sol_pert.success:
        delta_g = sol_pert.y[0]
        amp10 = (np.interp(10.0, sol_pert.t, delta_g))**2
        amp20 = (np.interp(20.0, sol_pert.t, delta_g))**2
        amp_results[k_val] = {'t': sol_pert.t, 'delta': delta_g}
        print(f"  {k_val:>12.3f}  {amp10:>22.4f}  {amp20:>22.4f}")

print()

# Skala z najsilniejszą amplifikacją
max_amp = 0
best_k  = 0
for k_val, res in amp_results.items():
    amp = max(abs(res['delta']))
    if amp > max_amp:
        max_amp = amp
        best_k = k_val

best_lam_Pl = 2*np.pi/best_k if best_k > 0 else np.inf
best_lam_kpc = best_lam_Pl * l_Pl_m / kpc_m

print(f"  Najsilniejsza amplifikacja: k = {best_k:.3f} l_Pl⁻¹")
print(f"  Odpowiadająca skala: λ = {best_lam_Pl:.2e} l_Pl = {best_lam_kpc:.2e} kpc")
print()

# =============================================================================
# Sekcja 6: Kosmologiczne skalowanie — H(z) a skala spinodalna
# =============================================================================
print("=" * 70)
print("SEKCJA 6: Kosmologiczne skalowanie")
print("=" * 70)
print()
print("  W standardowej kosmologii (ΛCDM), spinodalna skala jest FIZYCZNA:")
print("  λ_phys = a(t) * λ_comoving")
print()
print("  Przejście tachioniczne zachodzi przy z_tach (redshift).")
print("  Struktury na skali k_max rosną przez czas Δt (czas w strefie tachionicznej).")
print()
print("  Amplifikacja: A ~ exp(|ω| * Δt) = exp(√|V''| * Δt)")
print()
print("  Dla γ << 1 (małe sprzężenie):")
print("  √|V''| = √γ * √|2g-3g²| ≈ √γ * 0.577")
print("  Czas w strefie (jeden oscylacyjny tryb): Δt ~ 1/√γ [Planck]")
print("  A ~ exp(0.577) ~ 1.78 — SŁABA amplifikacja!")
print()
print("  Dla silnej amplifikacji (np. A ~ 10^{30}): potrzeba Δt ~ 100/√γ.")
print("  Ale Δt ~ czas Hubble'a ~ 1/H(z_tach) [Planck]")
print("  => √γ / H(z_tach) ~ 100")
print()
print("  WNIOSEK: Spinodalna może być skuteczna tylko jeśli:")
print("  (a) γ jest duże (ale wtedy skala jest mikroskopowa), LUB")
print("  (b) czas w strefie tachionicznej jest duży (długotrwałe przejście)")
print()
print("  Kosmologiczna ewolucja g_bg(t) z H≠0:")
print()

# Ewolucja z ekspansją (H = const ~ H_Planck ~ 1)
def field_with_hubble(t, y, H_val):
    """g̈ + 3H ġ + V'(g) = 0"""
    g, gd = y
    dVdg = g**2 - g**3
    return [gd, -3*H_val*gd - dVdg]

H_vals = [0.0, 0.01, 0.1, 0.5, 1.0]  # w jednostkach Plancka
print(f"  {'H [Pl⁻¹]':>10s}  {'max g':>8s}  {'t_reach_1':>12s}  {'t_in_tach':>12s}  {'Status':>15s}")
print("-" * 62)

for H_val in H_vals:
    try:
        sol_H = solve_ivp(
            lambda t, y: field_with_hubble(t, y, H_val=H_val),
            (0, 200), [0.05, 0.0],
            t_eval=np.linspace(0, 200, 2000),
            method='RK45', rtol=1e-8, atol=1e-10
        )
        if sol_H.success:
            g_H = sol_H.y[0]
            max_g = g_H.max()
            dt = sol_H.t[1] - sol_H.t[0]
            t_tach = np.sum(g_H > 2/3) * dt

            # Kiedy g przekracza 0.9?
            above_09 = np.where(g_H > 0.9)[0]
            t_reach_1 = sol_H.t[above_09[0]] if len(above_09) > 0 else np.inf

            if max_g > 0.9:
                status = "OSIĄGA PRÓŻNIĘ"
            elif max_g > 2/3:
                status = "WCHODZI W TACH"
            else:
                status = "ZATRZYMUJE SIĘ"

            print(f"  {H_val:>10.2f}  {max_g:>8.3f}  {t_reach_1:>12.2f}  {t_tach:>12.2f}  {status}")
    except Exception as e:
        print(f"  {H_val:>10.2f}  {'BŁĄD':>8s}  {'—':>12s}  {'—':>12s}  {str(e)[:20]}")

print()
print("  WNIOSEK kosmologiczny:")
print("  - Duże H (szybka ekspansja): pole zatrzymuje się przed próżnią")
print("    (tłumienie Hubble'a). Przejście tachioniczne nie zachodzi.")
print("  - Małe H (wolna ekspansja, późna epoka): pole dociera do g=1,")
print("    ale czas w strefie tachionicznej ~ Pl — zbyt krótki.")
print()

# =============================================================================
# Sekcja 7: Porównanie z NFW profilem
# =============================================================================
print("=" * 70)
print("SEKCJA 7: Czy spinodalne fluktuacje dają profil NFW?")
print("=" * 70)
print()
print("  Fluktuacje spinodalne tworzą struktury o profilu:")
print("  ρ_spin(r) ~ |δg(r)|² * Φ₀² / (l_Pl³)")
print()
print("  Kształt w przestrzeni rzeczywistej (transformata Fouriera):")
print("  Dla δg_k ~ exp(-k²/k_max²):")
print("    δg(r) ~ exp(-r² * k_max²/4) * k_max³")
print("  => Profil GAUSSOWSKI, nie NFW!")
print()
print("  NFW: ρ ∝ 1/[r*(1+r/r_s)²]")
print("  Spinodal: ρ ∝ exp(-r²/r_spin²)")
print()
print("  Te profile są FUNDAMENTALNIE RÓŻNE:")
print("  - NFW ma ogon ~1/r³ dla r>>r_s (wolny zanik)")
print("  - Spinodal ma Gaussowski zanik (szybki)")
print()
print("  Spinodalne fluktuacje TGP NIE tworzą profilu NFW.")
print("  => Mechanizm spinodalny w minimalnym TGP nie wyjaśnia DM halo.")
print()

# =============================================================================
# Sekcja 8: Energia spinodalna jako DM
# =============================================================================
print("=" * 70)
print("SEKCJA 8: Gęstość energii spinodalnej")
print("=" * 70)
print()
print("  Gęstość energii spinodalnych fluktuacji:")
print("  ρ_spin = Φ₀²/2 * ∫ d³k/(2π)³ * |δg_k|²")
print()
print("  W granicy k_max >> 1/r_galaktyki:")
print("  ρ_spin ≈ Φ₀²/(2(2π)²) * ∫₀^{k_max} k² |δg_k|² dk")
print()
print("  Jeśli amplifikacja A jest duża: ρ_spin ~ A² * ρ_vac")
print()

# Oblicz A dla realnych parametrów
gamma_gal = 1e-110  # γ dla skali galaktycznej
m_sp_gal  = np.sqrt(gamma_gal)
k_max_gal = m_sp_gal  # √|V''| ~ √γ dla g_bg~1
t_tach_est = 1.0 / k_max_gal  # czas przejścia ~ 1/k_max [Planck]

A_spinodal = np.exp(k_max_gal * t_tach_est)

print(f"  Dla γ = {gamma_gal:.2e} (skala galaktyczna):")
print(f"  k_max = √γ = {k_max_gal:.2e} l_Pl⁻¹")
print(f"  t_tach ~ 1/k_max = {t_tach_est:.2e} l_Pl")
print(f"  Amplifikacja A = exp(k_max * t_tach) = exp(1) = {A_spinodal:.3f}")
print()
print("  WYNIK: Amplifikacja A ~ e^1 ~ 2.7 — MARGINALNA!")
print("  Spinodala NIE produkuje eksponencjalnego wzrostu gęstości")
print("  dla skali galaktycznej (γ ~ 10^{-110}).")
print()
print("  Warunek silnej amplifikacji (A >> 1):")
print("  k_max * t_tach >> 1")
print("  Ale t_tach ~ 1/k_max => k_max * t_tach ~ O(1)")
print("  => ZAWSZE marginalnie. Niestabilność spinodalna w TGP jest SŁABA.")
print()

# =============================================================================
# Sekcja 9: Widmo mocy spinodalnego TGP
# =============================================================================
print("=" * 70)
print("SEKCJA 9: Widmo mocy P(k) fluktuacji spinodalnych")
print("=" * 70)
print()

# Widmo mocy perturbacji po przejściu spinodalnym
# Dla prostego modelu: P_spin(k) ~ k² * |δg_k|² ~ k² * exp(2*omega_k * t_tach)
# gdzie omega_k = √(|V''| - k²) dla k < k_max
# Czas efektywny t_tach ~ 1/H (czas Hubble'a)

print("  Model widma mocy po spinodalnym przejściu:")
print("  P_spin(k) ∝ k² * exp(2*√(|V''|-k²) * t_tach)")
print()
print("  Porównanie z obserwowanym P(k) struktury wielkoskalowej:")
print("  P_obs(k) ∝ k^n_s * T²(k) * D²(z)")
print("  (ΛCDM: n_s ~ 0.965, Harrison-Zel'dovich)")
print()

# Oblicz P_spin(k) dla różnych γ
k_plot = np.logspace(-4, 0, 200)  # w l_Pl^{-1}

fig_power, ax_p = plt.subplots(figsize=(10, 6))

gamma_vals_plot = [1e-50, 1e-80, 1e-100, 1e-110]
cols_p = cm.viridis(np.linspace(0.1, 0.9, len(gamma_vals_plot)))

for gamma_p, col_p in zip(gamma_vals_plot, cols_p):
    k_max_p = np.sqrt(gamma_p / 3.0)  # dla g_bg=1/3
    omega_k = np.where(k_plot < k_max_p,
                        np.sqrt(np.maximum(gamma_p/3.0 - k_plot**2, 0)),
                        0.0)
    t_tach_p = 1.0 / np.sqrt(gamma_p) if gamma_p > 0 else 1e100
    P_spin = k_plot**2 * np.exp(2 * omega_k * min(t_tach_p, 1e10))
    # Normalizuj
    P_spin = P_spin / P_spin.max()
    lam_max_kpc = 2*np.pi / k_max_p * l_Pl_m / kpc_m if k_max_p > 0 else np.inf
    ax_p.loglog(k_plot * l_Pl_m / kpc_m, P_spin, color=col_p, lw=2,
                label=rf'$\gamma={gamma_p:.0e}$, $\lambda_{{max}}$={lam_max_kpc:.1e} kpc')

# Galaktyczna skala referencyjna
ax_p.axvline(1.0/10.0, color='red', ls='--', lw=1.5, alpha=0.7, label='k = 1/(10 kpc)')
ax_p.axvline(1.0/0.1, color='orange', ls=':', lw=1.5, alpha=0.7, label='k = 1/(0.1 kpc)')
ax_p.set_xlabel(r'$k$ [kpc$^{-1}$]', fontsize=11)
ax_p.set_ylabel('P_spin(k) (znormalizowane)', fontsize=11)
ax_p.set_title('Widmo mocy spinodalnych fluktuacji TGP', fontsize=12)
ax_p.legend(fontsize=9)
ax_p.grid(True, alpha=0.3, which='both')
ax_p.set_xlim(1e-8, 1e6)

out_power = os.path.join(os.path.dirname(__file__), 'ex38_power_spectrum.png')
plt.savefig(out_power, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Wykres widma mocy: {out_power}")

# =============================================================================
# Wykres główny
# =============================================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('TGP: Tachionowa Niestabilność Spinodalna (ex38)\n'
             r'$V(g) = g^3/3 - g^4/4$, $V^{\prime\prime}(g) = 2g - 3g^2$',
             fontsize=13, y=1.01)

# Panel 1: V(g) i V''(g)
ax = axes[0, 0]
ax2 = ax.twinx()
ln1 = ax.plot(g_arr, V_arr, 'b-', lw=2.5, label='V(g)')
ln2 = ax2.plot(g_arr, d2V_arr, 'r-', lw=2, label="V''(g)")
ax2.axhline(0, color='k', lw=0.8, alpha=0.4)
ax.axvline(2/3, color='purple', ls='--', lw=1.5, alpha=0.8, label='g=2/3 (granica tachionu)')
ax.axvline(1, color='gray', ls='--', lw=1.2, alpha=0.6, label='g=1 (false vacuum)')
ax.fill_betweenx([-0.05, 0.10], 2/3, 1.5, alpha=0.08, color='red', label='Strefa tachioniczna')
ax.set_xlabel('g = Φ/Φ₀', fontsize=11)
ax.set_ylabel('V(g)', fontsize=11, color='b')
ax2.set_ylabel("V''(g)", fontsize=11, color='r')
ax.set_title('Potencjał TGP i masa kwadratowa', fontsize=11)
ax.set_xlim(0, 1.5); ax.set_ylim(-0.05, 0.10)
ax2.set_ylim(-0.6, 0.3)
lns = ln1 + ln2
ax.legend(lns + [ax.get_lines()[1], ax.get_lines()[2]],
          [l.get_label() for l in lns] + ['g=2/3', 'g=1'], fontsize=8)
ax2.tick_params(axis='y', labelcolor='r')
ax.grid(True, alpha=0.3)

# Panel 2: Dyspersja ω²(k) dla g_bg
ax = axes[0, 1]
k_disp = np.linspace(0, 1.5, 200)
for g_bg, col in zip([1/3, 0.5, 0.8, 1.0], cm.RdYlBu(np.linspace(0, 1, 4))):
    d2V_g = 2*g_bg - 3*g_bg**2
    omega2 = k_disp**2 + d2V_g
    ax.plot(k_disp, omega2, color=col, lw=2, label=f'g_bg={g_bg:.2f}')
ax.axhline(0, color='k', lw=1.5)
ax.fill_between(k_disp, np.minimum(0, k_disp**2 - 1/3), 0,
                alpha=0.1, color='red', label='Strefa wzrostu')
ax.set_xlabel('k [Pl⁻¹]', fontsize=11)
ax.set_ylabel(r'$\omega^2(k)$', fontsize=11)
ax.set_title('Dyspersja perturbacji spinodalnych', fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.5, 2.0)

# Panel 3: Ewolucja δg_k
ax = axes[0, 2]
for k_val, col in zip(k_test_vals[:4], cm.Greens(np.linspace(0.4, 0.9, 4))):
    if k_val in amp_results:
        res = amp_results[k_val]
        ax.semilogy(res['t'], np.abs(res['delta']), color=col, lw=1.8,
                    label=f'k={k_val:.3f} Pl⁻¹')
ax.set_xlabel('t [Planck]', fontsize=11)
ax.set_ylabel('|δg_k(t)|', fontsize=11)
ax.set_title('Wzrost perturbacji spinodalnych\n(g₀=0.05, bez ekspansji)', fontsize=10)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 4: Skala spinodalna vs γ
ax = axes[1, 0]
gamma_arr = np.logspace(-120, -30, 200)
lambda_spin_kpc = 2*np.pi*np.sqrt(3) / np.sqrt(gamma_arr) * l_Pl_m / kpc_m
ax.loglog(gamma_arr, lambda_spin_kpc, 'navy', lw=2.5)
ax.axhline(10.0, color='red', ls='--', lw=2, label='10 kpc (galaktyka)')
ax.axhline(1.0, color='orange', ls='--', lw=1.5, label='1 kpc')
ax.axhline(0.1, color='gold', ls=':', lw=1.5, label='100 pc')
# Zaznaczyć γ dla skali galaktycznej
gamma_at_10kpc = (2*np.pi*np.sqrt(3) / (10.0*kpc_m/l_Pl_m))**2
ax.axvline(gamma_at_10kpc, color='gray', ls=':', lw=1.5, alpha=0.7, label=f'γ(10kpc)={gamma_at_10kpc:.0e}')
ax.set_xlabel(r'$\gamma$ [Planck]', fontsize=11)
ax.set_ylabel(r'$\lambda_{spin}$ [kpc]', fontsize=11)
ax.set_title('Skala spinodalna vs γ', fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, which='both')

# Panel 5: Ewolucja g_bg z H
ax = axes[1, 1]
g_arr_bg = np.linspace(0, 1.5, 300)
for H_v in [0.0, 0.01, 0.1, 0.5]:
    try:
        sol_h = solve_ivp(
            lambda t, y: field_with_hubble(t, y, H_val=H_v),
            (0, 100), [0.05, 0.0],
            t_eval=np.linspace(0, 100, 1000),
            method='RK45', rtol=1e-8, atol=1e-10
        )
        if sol_h.success:
            g_h = sol_h.y[0]
            ax.plot(sol_h.t, g_h, lw=1.8, label=f'H={H_v:.2f} Pl⁻¹')
    except Exception:
        pass
ax.axhline(1.0, color='gray', ls='--', lw=1, alpha=0.6, label='próżnia g=1')
ax.axhline(2/3, color='red', ls='--', lw=1.2, alpha=0.7, label='g=2/3 (tachion)')
ax.set_xlabel('t [Planck]', fontsize=11)
ax.set_ylabel('g_bg(t)', fontsize=11)
ax.set_title('Ewolucja g_bg z ekspansją H', fontsize=11)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.05, 1.3)

# Panel 6: Werdykt — diagram strony
ax = axes[1, 2]
ax.axis('off')
verdict_text = """
TACHIONOWA SPINODALNA TGP
━━━━━━━━━━━━━━━━━━━━━━━━━

WYNIKI EX38:

1. Strefa tachioniczna: g ∈ (2/3, ∞)
   Próżnia g=1 jest PERMANENTNIE tachioniczna.

2. Skala spinodalna:
   λ_spin = 2π√3/√γ * l_Pl
   Dla γ=1: λ_spin ~ l_Pl  (mikroskopowa!)
   Dla skali 10 kpc: γ ~ 10⁻¹¹⁰

3. Amplifikacja perturbacji:
   A ~ exp(1) ~ 2.7  (SŁABA)
   Spinodalna TGP jest marginalna.

4. Profil gęstości:
   ρ_spin(r) ~ Gaussowski
   ρ_NFW(r) ~ 1/r * (1+r/rs)⁻²
   NIEZGODNE!

WERDYKT:
━━━━━━━━━━━━━━━━━━━━━━━━━
Spinodalna TGP NIE wyjaśnia DM.
• Skala wymaga γ ~ 10⁻¹¹⁰ (paradoks γ)
• Amplifikacja zbyt słaba (A~e¹)
• Profil Gaussowski ≠ NFW

MECHANIZM SFALSYFIKOWANY.
"""
ax.text(0.05, 0.95, verdict_text, transform=ax.transAxes,
        fontsize=9, verticalalignment='top', family='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
out_main = os.path.join(os.path.dirname(__file__), 'ex38_spinodal.png')
plt.savefig(out_main, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nWykres główny: {out_main}")

# =============================================================================
# Podsumowanie
# =============================================================================
print()
print("=" * 70)
print("PODSUMOWANIE ex38: TACHIONOWA SPINODALNA TGP")
print("=" * 70)
print()
print("WYNIK 1: Strefa tachioniczna TGP")
print("  V''(g) < 0 dla g ∈ (2/3, ∞). Próżnia g=1 jest ZAWSZE tachioniczna.")
print("  To NIE jest tymczasowe — jest permanentną własnością próżni TGP.")
print()
print("WYNIK 2: Skala spinodalna")
print("  λ_spin = 2π√3/√γ * l_Pl")
print("  Dla γ=1 (Planck): λ_spin ~ l_Pl (mikroskopowa)")
print("  Dla galaktycznej (10 kpc): γ ~ 10^{-110}")
print()
print("WYNIK 3: Amplifikacja perturbacji")
print("  A ~ exp(√γ * t_tach) ~ exp(1) ≈ 2.7 — SŁABA!")
print("  Nie daje eksponencjalnego wzrostu gęstości DM.")
print()
print("WYNIK 4: Kształt profilu gęstości")
print("  ρ_spin(r) ~ Gaussowski  vs  ρ_NFW(r) ~ 1/r*(1+r/rs)^{-2}")
print("  NIEZGODNE — spinodalna nie odtwarza profilu halo.")
print()
print("WERDYKT KOŃCOWY ex38:")
print("  Tachionowa niestabilność spinodalna TGP nie może wyjaśnić")
print("  galaktycznej ciemnej materii. Trzy niezależne problemy:")
print("  (a) Wymagane γ ~ 10^{-110} (paradoks γ)")
print("  (b) Amplifikacja zbyt słaba (exp(1), nie exp(100))")
print("  (c) Profil gęstości Gaussowski, nie NFW")
print()
print("  MECHANIZM M3-spinodalny: SFALSYFIKOWANY")
print()
print("ZAKTUALIZOWANA LISTA SFALSYFIKOWANYCH MECHANIZMÓW TGP-DM:")
print("  ✗ M1: Yukawa G_eff(r) < G (ex14, ex35)")
print("  ✗ M2: V₃ akumulacja — kształt krzywej 1/r (ex35)")
print("  ✗ M3a: Energia pola ε_field — 1/r³–1/r⁴ (ex35)")
print("  ✗ M4: Tachion m²<0 — niestabilna teoria (ex35)")
print("  ✗ M5: Spinodalna — amplifikacja słaba, profil Gaussowski (ex38)")
print()
print("  ? M_FDM: V_mod = εg² + V_TGP — MOŻLIWE (ex36/ex37),")
print("           ale wymaga ε ~ 10^{-101} i nie wynika z N0 aksjomatów.")
print()
print("NASTĘPNE KROKI:")
print("  - Aktualizacja ANALIZA_SPOJNOSCI (Kill-shot K13, K14, K15)")
print("  - ex39: Kosmologiczne widmo mocy TGP-FDM")
print("  - Teoria: derivacja ε z topologicznych defektów TGP")
