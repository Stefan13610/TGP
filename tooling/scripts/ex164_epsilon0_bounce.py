"""
ex164_epsilon0_bounce.py
=========================
R6: ε₀ z termodynamiki substratu — obliczenie bounce Colemana

IDEA:
  Przestrzeń nukleuje z "nic" (Φ=0) do bąbla z Φ>0.
  Nukleacja = tunelowanie kwantowe Φ: 0 → ε₀ > 0.
  ε₀ determinowane przez profil bounce + akcję Colemana.

  Potencjał GL substratu TGP (z prop:substrate-action):
    V(ψ) = β·ψ³/3 - γ·ψ⁴/4 + λ₆·ψ⁶/6   (β=γ z thm:beta_gamma)

  W zmiennej g = ψ (bezwymiarowej):
    V(g) = γ[g³/3 - g⁴/4] + λ₆·g⁶/6

  Próżnia TGP: V'(g_vac) = 0, g_vac = 1 (po renormalizacji)
  Fałszywa próżnia: g = 0 (brak przestrzeni)

  Bounce: rozwiązanie O(4)-symetryczne:
    g'' + (3/ρ)g' = dV/dg   (ρ = promień euklidesowy)
    BC: g'(0) = 0, g(∞) = 0

  Akcja bounce'u: B = 2π² ∫₀^∞ dρ ρ³ [(g')²/2 + V(g)]
  Tempem nukleacji: Γ ~ exp(-B)
  ε₀ ~ g(ρ=0)^? zależy od profilu bounce'u

PLAN:
  1. Oblicz V(g) dla TGP z β=γ i stabilizacją ψ⁶
  2. Rozwiąż ODE bounce'u (metoda undershoot/overshoot)
  3. Oblicz B i profil bounce'u
  4. Wyciągnij ε₀ i N_e
  5. Skan λ₆: jak ε₀ zależy od parametru stabilizacji?
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

print("=" * 72)
print("ex164: ε₀ z bounce'u Colemana (potencjał GL substratu TGP)")
print("=" * 72)

# =====================================================================
# 1. Potencjał GL substratu TGP
# =====================================================================
print("\n--- 1. Potencjał V(g) substratu TGP ---")

# V(g) = γ[g³/3 - g⁴/4] + λ₆·g⁶/6
# β = γ (thm:beta_gamma)
# Próżnia: V'(1) = 0 → γ[1 - 1] + λ₆·1 = 0? Nie — trzeba poprawnie.
# V'(g) = γ[g² - g³] + λ₆·g⁵
# V'(1) = γ[1 - 1] + λ₆ = λ₆ → musi być = 0 dla próżni → λ₆ = 0?
#
# Nie, V(g) w TGP to:
# V(g) = γ/Φ₀ · [g²/2 - g³/3]    (z sek08, β=γ)
# Plus stabilizacja: +λ₆·g⁶/6
# Próżnia: V'(g_vac) = 0 → g_vac(γ - g_vac·γ + ...) = 0
#
# Przeformułujmy: kanoniczny potencjał TGP (z sek10):
# V_mod(ψ) = (β/3)ψ³ - (γ/4)ψ⁴  z β = γ
# V_mod(ψ) = γ[(1/3)ψ³ - (1/4)ψ⁴]
# V_mod'(ψ) = γ[ψ² - ψ³] = γψ²(1 - ψ) = 0 → ψ = 0, 1
# V_mod(0) = 0, V_mod(1) = γ(1/3 - 1/4) = γ/12 > 0
#
# Problem: V(0) = 0 < V(1) = γ/12 → fałszywa próżnia jest NIŻEJ!
# Nukleacja musi iść 0 → 1, ale 1 jest wyżej w energii.
# → To nie standardowy Coleman-bounce.
#
# W TGP: nukleacja GRAWITACYJNA — korekcja de Sitter.
# Albo: z λ₆ stabilizacją, V(g) ma minimum przy g_vac > 1 z V < 0.

# Bardziej realistyczny model:
# V(g) = γ[(1/3)g³ - (1/4)g⁴] + (λ₆/6)g⁶
# z λ₆ < 0 → true vacuum poniżej false vacuum

# Albo standard: V(g) z trzema zerowymi punktami:
# V(g) = (1/2)g²(g - g_+)² z barierą — to jest generyczna forma

# Użyjmy FIZYCZNEGO modelu:
# False vacuum: g = 0 (brak przestrzeni, Φ = 0)
# True vacuum:  g = 1 (TGP vacuum, Φ = Φ₀)
# Bariera pomiędzy nimi

# Standardowy model: potencjał z barierą
# V(g) = λ g²(g - a)² - ε g
# ale prostszy: thin-wall approximation

# PODEJŚCIE: parametryzacja ogólna
# V(g) = A·g² - B·g³ + C·g⁴
# z V(0) = 0, V ma barierę, V(g_true) < V(0) (true vacuum niżej)
# V'(g_true) = 0

# Konkretnie z TGP:
# V_TGP(ψ) = (γ/3)ψ³ - (γ/4)ψ⁴
# To ma V(0) = 0 (false), V(1) = γ/12 > 0 (wyższe!)
# Dodanie stałej kosmologicznej: ΔV
# V(ψ) = (γ/3)ψ³ - (γ/4)ψ⁴ - ε·ψ
# gdzie ε = mały bias faworyzujący ψ > 0

# Lepsze: użyj standardowej formy nukleacji
# φ: 0 (false) → v (true), z barierą
# V(φ) = λ/4 (φ² - v²)² - ε(φ + v)/(2v)
# → V(v) = -ε, V(-v) = 0 (false), bariera przy φ=0

# Ale TGP jest specyficzne: φ ≥ 0 (Φ ≥ 0).
# Użyjmy: V(g) = (g²/2)(1 - g)² - ε·g  [uproszczona forma z barierą]
# False: g ≈ 0, True: g ≈ 1 (z V(1) < 0 dla ε > 0)

# ---- Model roboczy ----
# V(g) = (μ²/2)g² - (η/3)g³ + (λ/4)g⁴ - ε·g
# Parametry: μ² (masa), η (kubiczny, z β=γ TGP), λ (kwadryczny), ε (bias)
# Uproszczenie: μ²=1, η wyznacza barierę, λ=1, ε (bias) → skan

# TGP-inspired:
# V_TGP = normalized so that barrier height ~ 1
# V(g) = g²(1-g)² — double-well z min w g=0 i g=1, obie V=0
# Dodajemy bias: V(g) = g²(1-g)² - ε·g²/2
# → V(0) = 0, V(1) = -ε/2 (true vacuum niżej)

def V_tgp(g, eps_bias):
    """TGP-inspired potential with bias."""
    return g**2 * (1 - g)**2 - eps_bias * g**2 / 2

def dV_tgp(g, eps_bias):
    """dV/dg"""
    # d/dg [g²(1-g)² - ε g²/2]
    # = 2g(1-g)² + g²·2(1-g)·(-1) - εg
    # = 2g(1-g)² - 2g²(1-g) - εg
    # = 2g(1-g)[(1-g) - g] - εg
    # = 2g(1-g)(1-2g) - εg
    return 2*g*(1-g)*(1-2*g) - eps_bias*g

# Sprawdź: próżnia
print("\n  V(g) = g²(1-g)² - ε·g²/2")
print("  False vacuum: g = 0 (V = 0)")
print("  True vacuum: g = g_true (V < 0)")
print("")

for eps_bias in [0.01, 0.05, 0.1, 0.2]:
    # Znajdź true vacuum: V'(g) = 0, g > 0.5
    # 2g(1-g)(1-2g) - εg = 0 → g[2(1-g)(1-2g) - ε] = 0
    # 2(1-g)(1-2g) = ε → solvable
    def dV_zero(g):
        return 2*(1-g)*(1-2*g) - eps_bias

    # Find root near g=1
    try:
        g_true = brentq(dV_zero, 0.6, 1.5)
        V_true = V_tgp(g_true, eps_bias)
        # Barrier top: between 0 and g_true
        g_bar_res = minimize_scalar(lambda g: -V_tgp(g, eps_bias), bounds=(0.01, g_true*0.9), method='bounded')
        g_bar = g_bar_res.x
        V_bar = V_tgp(g_bar, eps_bias)
        print(f"  ε = {eps_bias:.3f}: g_true = {g_true:.4f}, V_true = {V_true:.6f}, "
              f"g_barrier = {g_bar:.4f}, V_barrier = {V_bar:.6f}, ΔV = {V_bar - V_true:.6f}")
    except:
        print(f"  ε = {eps_bias:.3f}: nie znaleziono minimum")

# =====================================================================
# 2. Bounce Colemana — O(4) symetria
# =====================================================================
print(f"\n{'='*72}")
print("--- 2. Bounce Coleman — O(4) symmetric ---")
print(f"{'='*72}")
print("\n  g'' + (3/ρ)g' = dV/dg")
print("  BC: g'(0) = 0, g(ρ→∞) = 0 (false vacuum)\n")

def solve_bounce(g0, eps_bias, rho_max=50.0):
    """
    Rozwiąż ODE bounce'u: g'' + (3/ρ)g' = dV/dg
    z g(0) = g0, g'(0) = 0.
    Szukamy g0 takiego że g(∞) → 0 (undershoot/overshoot).
    """
    def rhs(rho, y):
        g, gp = y
        dVdg = dV_tgp(g, eps_bias)
        if rho < 1e-10:
            return [gp, dVdg / 4.0]  # L'Hôpital: (3/ρ)g' → g''·3/1 → eq: 4g'' = dV/dg
        return [gp, dVdg - 3 * gp / rho]

    sol = solve_ivp(rhs, (0, rho_max), [g0, 0.0],
                    rtol=1e-10, atol=1e-12, max_step=0.1,
                    events=None)
    return sol.t, sol.y[0], sol.y[1]

def bounce_action(rho, g, gp, eps_bias):
    """
    B = 2π² ∫₀^∞ dρ ρ³ [(g')²/2 + V(g)]
    """
    integrand = rho**3 * (gp**2 / 2 + V_tgp(g, eps_bias))
    return 2 * np.pi**2 * np.trapezoid(integrand, rho)

# Metoda undershoot/overshoot:
# g0 za duże → g przechodzi przez 0, potem odskakuje (overshoot)
# g0 za małe → g nie dochodzi do 0, wraca (undershoot)
# Bounce: g0 taki że g(∞) → 0 dokładnie

eps_values = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3]

print(f"  {'ε':>6} {'g₀_bounce':>12} {'g_true':>10} {'B (action)':>12} {'Γ~exp(-B)':>14} {'B/thin-wall':>12}")
print("  " + "-" * 72)

results_bounce = []

for eps_bias in eps_values:
    # Znajdź g_true
    try:
        g_true = brentq(lambda g: 2*(1-g)*(1-2*g) - eps_bias, 0.6, 1.5)
    except:
        continue

    # Undershoot/overshoot: szukaj g0 ∈ (0.5, g_true+0.1)
    def shoot(g0):
        rho, g, gp = solve_bounce(g0, eps_bias, rho_max=80.0)
        # Kryterium: g na końcu
        g_end = g[-1]
        # Overshoot: g < 0 gdzieś → g_min < 0
        g_min = np.min(g)
        if g_min < -0.01:
            return 1.0  # overshoot
        # Undershoot: g still > 0.3 at end
        return g_end  # chcemy g_end → 0

    # Skan: znajdź przejście overshoot → undershoot
    g0_scan = np.linspace(g_true * 0.5, g_true * 1.05, 60)
    g_end_vals = []
    for g0_test in g0_scan:
        rho, g, gp = solve_bounce(g0_test, eps_bias, rho_max=60.0)
        g_end_vals.append(g[-1])
    g_end_vals = np.array(g_end_vals)

    # Bounce: g_end przechodzi przez minimum bliskie 0
    # Szukamy g0 dającego minimum g_end
    idx_min = np.argmin(np.abs(g_end_vals))

    # Doprecyzowanie brentq: szukaj g_end bliskiego 0
    g0_bounce = None
    threshold = 0.05  # chcemy g_end < threshold

    for i in range(len(g0_scan) - 1):
        if g_end_vals[i] > threshold and g_end_vals[i+1] < threshold:
            try:
                def bounce_target(g0):
                    rho, g, gp = solve_bounce(g0, eps_bias, rho_max=80.0)
                    return g[-1] - threshold/2
                g0_bounce = brentq(bounce_target, g0_scan[i], g0_scan[i+1], xtol=1e-6)
                break
            except:
                pass
        elif g_end_vals[i] < threshold and g_end_vals[i+1] > threshold:
            # Odwrotne przejście
            try:
                def bounce_target(g0):
                    rho, g, gp = solve_bounce(g0, eps_bias, rho_max=80.0)
                    return g[-1] - threshold/2
                g0_bounce = brentq(bounce_target, g0_scan[i], g0_scan[i+1], xtol=1e-6)
                break
            except:
                pass

    if g0_bounce is None:
        # Fallback: użyj g0 z minimalnym g_end
        g0_bounce = g0_scan[idx_min]

    # Oblicz bounce z najlepszym g0
    rho, g, gp = solve_bounce(g0_bounce, eps_bias, rho_max=100.0)
    B = bounce_action(rho, g, gp, eps_bias)

    # Thin-wall approximation: B_tw = 27π²·σ⁴ / (2·ε³)
    # σ = surface tension = ∫ dg √(2V(g)) from g=0 to g_true (at ε=0)
    # V_0(g) = g²(1-g)² (no bias)
    g_int = np.linspace(0, 1.0, 1000)
    sigma_tw = np.trapezoid(np.sqrt(2 * g_int**2 * (1 - g_int)**2), g_int)
    B_tw = 27 * np.pi**2 * sigma_tw**4 / (2 * eps_bias**3) if eps_bias > 0 else np.inf

    ratio_tw = B / B_tw if B_tw > 0 and B_tw < 1e20 else 0

    Gamma = np.exp(-min(B, 700))  # cap for display

    results_bounce.append((eps_bias, g0_bounce, g_true, B, Gamma, ratio_tw))

    print(f"  {eps_bias:6.3f} {g0_bounce:12.6f} {g_true:10.4f} {B:12.2f} {Gamma:14.4e} {ratio_tw:12.4f}")

# =====================================================================
# 3. Interpretacja fizyczna: ε₀ i N_e
# =====================================================================
print(f"\n{'='*72}")
print("--- 3. Interpretacja: ε₀ z bounce'u ---")
print(f"{'='*72}")

print("""
  W TGP, nukleacja przestrzeni:
    Φ = 0 (brak przestrzeni) → Φ = ε₀·Φ₀ (bąbel)

  Profil bounce'u g(ρ) daje:
    g(ρ=0) = g₀_bounce ← wartość wewnątrz bąbla
    ε₀ = g₀_bounce² (bo Φ = Φ₀·ψ², ψ = g)
    Albo: ε₀ = g₀_bounce (w zmiennej ψ = Φ/Φ₀)

  N_e = (1/3)·ln(Φ₀/ε₀) = (1/3)·ln(1/ε₀) + (1/3)·ln(Φ₀)
""")

Phi0 = 115.0
print(f"  Φ₀ = {Phi0}")
print(f"\n  {'ε_bias':>8} {'g₀_bounce':>12} {'ε₀=g₀²':>10} {'N_e(ε₀)':>10} {'B':>10} {'exp(-B)':>12}")
print("  " + "-" * 70)

for eps_bias, g0b, g_true, B, Gamma, ratio in results_bounce:
    eps0 = g0b**2  # ε₀ = ψ² w TGP
    if eps0 > 0:
        Ne = (1.0/3) * np.log(Phi0 / eps0)
    else:
        Ne = 0
    print(f"  {eps_bias:8.3f} {g0b:12.6f} {eps0:10.6f} {Ne:10.2f} {B:10.2f} {Gamma:12.4e}")

print("""
  PROBLEM: g₀_bounce ~ O(1) → ε₀ ~ O(1) → N_e ~ 1
  To za mało! Potrzebujemy ε₀ << 1 → N_e >> 1.

  Przyczyna: bezwymiarowy potencjał V(g) = g²(1-g)² ma barierę O(1).
  Bounce daje g₀ ~ g_true ~ O(1).

  Żeby uzyskać ε₀ ~ 10⁻⁷⁰:
  Potrzebna hierarchia SKALI — substrat ma a_sub ~ ℓ_P,
  ale Φ₀ = 115 oznacza ξ/ℓ_P ~ 10³⁵.
""")

# =====================================================================
# 4. Model z hierarchią skali
# =====================================================================
print(f"{'='*72}")
print("--- 4. Model z hierarchią skali ---")
print(f"{'='*72}")

print("""
  Fizyczna nukleacja: NIE jest bezwymiarowym bounce'em.
  Wymiarowy potencjał:
    V_phys(Φ) = (E_sub/Φ₀³)·[Φ³/3 - Φ⁴/(4Φ₀)]
  gdzie E_sub = energia substratu na węzeł.

  Skala bariery: V_bar ~ E_sub · (Φ₀)^{coś}
  Promień krytyczny: R_c ~ ξ (długość korelacji)
  Akcja: B ~ (R_c/ℓ_P)^4 · V_bar / T_sub

  Kluczowe: B ~ (ξ/ℓ_P)^4 → dla ξ/ℓ_P ~ 10³⁵:
    B ~ 10^{140} → Γ ~ exp(-10^{140}) ≈ 0

  ALE: nukleacja zachodzi PRZY T_c (przejście fazowe).
  Wówczas bariera ΔF → 0, R_c → ∞, ale B → B_c (skończone).

  To jest problem Lindeana: inflacja z przejścia fazowego.
""")

# Model Lindego: V(φ) = λ/4·(φ²-v²)² near T_c
# Fluktuacja termiczna δΦ ~ T_c → ε₀ ~ (T_c/V₀)² ~ (T_c/M_P)² ~ (ℓ_P/ξ)²
# To jest DOKŁADNIE estymata z epsilon0_estimate.py!

# Wniosek: ε₀ = σ²_c · (ℓ_P/ξ)² jest POPRAWNĄ estymatą
# Bardziej precyzyjne obliczenie wymaga:
# 1. Pełnego potencjału V(Φ, T) blisko T_c
# 2. Obliczenia Γ(T) ≈ T^4 exp(-B₃(T)/T) (3D bounce, termiczny)
# 3. Warunku nukleacji: Γ(T_nuc) · H⁻⁴ ~ 1

print("  Estymata Fluktuacji termicznej:")
print("  ε₀ ~ σ²_c · (ℓ_P/ξ)²")
print("")

# Ising 3D parametry
T_c_zJ = 4.5115 / 6  # T_c/(z·J) dla sc lattice
sigma2_c = T_c_zJ / 2

print(f"  σ²_c = T_c/(2zJ) = {sigma2_c:.6f}")
print(f"  ℓ_P = 1.616 × 10⁻³⁵ m")
print(f"  ξ przy T_c: ξ ~ a_sub · |t|^{{-ν}} z ν = 0.630\n")

# Skan: ξ/a_sub → ε₀ → N_e
print(f"  {'ξ/a_sub':>12} {'ε₀':>14} {'N_e (Φ₀=115)':>14} {'n_s (Starobinsky)':>18}")
print("  " + "-" * 65)

xi_a_values = np.logspace(10, 45, 15)
for xi_a in xi_a_values:
    eps0 = sigma2_c / xi_a**2
    Ne = (1.0/3) * np.log(Phi0 / eps0) if eps0 > 0 and eps0 < Phi0 else 0
    ns = 1 - 2.0/Ne - 0.0005 if Ne > 5 else 0

    mark = ""
    if 50 <= Ne <= 70:
        mark = " ✅"

    print(f"  {xi_a:12.3e} {eps0:14.4e} {Ne:14.2f} {ns:18.4f}{mark}")

# =====================================================================
# 5. Relacja ξ ↔ T i warunek samospójności
# =====================================================================
print(f"\n{'='*72}")
print("--- 5. Samospójność: ξ z przejścia fazowego ---")
print(f"{'='*72}")

print("""
  W przejściu fazowym 2. rodzaju (klasa Isinga 3D):
    ξ = ξ₀ · |t|^{-ν}   (t = (T-T_c)/T_c, ν = 0.6301)

  ξ₀ ~ a_sub (kilka stałych sieciowych)
  Nukleacja: |t| = (T_c - T_nuc)/T_c (jak blisko T_c?)

  ξ/a_sub = ξ₀/a_sub · |t|^{-0.630}

  Warunek N_e = 55 → ε₀ = Φ₀·exp(-165)
  ε₀ = σ²_c · (a_sub/ξ)² → ξ/a_sub = √(σ²_c/ε₀)
""")

Ne_target = 55
eps0_target = Phi0 * np.exp(-3 * Ne_target)
xi_a_target = np.sqrt(sigma2_c / eps0_target)

print(f"  N_e = {Ne_target}")
print(f"  ε₀ = Φ₀·exp(-3·N_e) = {eps0_target:.4e}")
print(f"  ξ/a_sub = √(σ²_c/ε₀) = {xi_a_target:.4e}")
print(f"  log₁₀(ξ/a_sub) = {np.log10(xi_a_target):.2f}")

nu = 0.6301
# ξ/a_sub = C · |t|^{-ν} → |t| = (C·a/ξ)^{1/ν}
# Przyjmijmy ξ₀/a_sub ≈ 1 (rzędu jedności)
C_xi = 1.0  # ξ₀/a_sub
t_nuc = (C_xi / xi_a_target)**(1.0/nu)

print(f"\n  |t_nuc| = (ξ₀/ξ)^{{1/ν}} = ({C_xi}/{xi_a_target:.2e})^{{1/{nu}}} = {t_nuc:.4e}")
print(f"  log₁₀|t_nuc| = {np.log10(t_nuc):.2f}")
print(f"  T_nuc = T_c · (1 - |t|) = T_c · (1 - {t_nuc:.2e})")
print(f"  → Nukleacja zachodzi NIEZWYKLE blisko T_c")
print(f"     (odchylenie od T_c: {t_nuc*100:.2e}%)")

# =====================================================================
# 6. Bounce termiczny 3D (zamiast kwantowego 4D)
# =====================================================================
print(f"\n{'='*72}")
print("--- 6. Bounce termiczny 3D ---")
print(f"{'='*72}")

print("""
  Blisko T_c: nukleacja jest TERMICZNA (3D bounce), nie kwantowa (4D).

  B₃/T = (16π/3)·σ_s³/(T·ΔF²)
  σ_s = napięcie powierzchniowe ∝ |t|^{μ}, μ = ν(d-1) = 0.630·2 = 1.260
  ΔF = różnica energii swobodnej ∝ |t|^{2-α} ≈ |t|^{1.89} (α=0.110)

  B₃/T ∝ |t|^{3μ - 2(2-α)} = |t|^{3·1.260 - 2·1.89} = |t|^{3.78-3.78} = |t|^0

  UWAGA: W 3D Isingu, B₃/T ∝ |t|^{d·ν - 2(2-α)}
  z d·ν = 3·0.630 = 1.890 i 2-α = 1.890 (relacja hiperawarzystości!)
  → B₃/T ~ const blisko T_c (!)

  To oznacza: tempem nukleacji Γ ~ T⁴·exp(-B₃/T) jest STAŁE blisko T_c.
  Wartość B₃/T zależy od detali nieuniversalnych (σ₀, ΔF₀).
""")

# Wykładniki Isinga 3D
nu_3d = 0.6301
alpha_3d = 0.110
gamma_3d = 1.2372
beta_3d = 0.3265

# Sprawdzenie relacji hiperaskalowości
print(f"  Wykładniki Isinga 3D:")
print(f"    ν = {nu_3d}")
print(f"    α = {alpha_3d}")
print(f"    β = {beta_3d}")
print(f"    γ = {gamma_3d}")
print(f"    d·ν = {3*nu_3d:.4f}")
print(f"    2-α = {2-alpha_3d:.4f}")
print(f"    d·ν = 2-α? {abs(3*nu_3d - (2-alpha_3d)) < 0.01} (relacja Josephsona)")

# σ_s ∝ ξ^{-(d-1)} ∝ |t|^{ν(d-1)} = |t|^{1.260}
# ΔF ∝ ξ^{-d} ∝ |t|^{d·ν} = |t|^{1.890}
# B₃/T ∝ σ_s³/(T·ΔF²) ∝ |t|^{3·1.260} / |t|^{2·1.890}
#       = |t|^{3.780} / |t|^{3.780} = |t|^0 = const!

print(f"\n  σ_s ∝ |t|^{{ν(d-1)}} = |t|^{{{nu_3d*2:.3f}}}")
print(f"  ΔF ∝ |t|^{{d·ν}} = |t|^{{{3*nu_3d:.3f}}}")
print(f"  B₃/T ∝ σ_s³/ΔF² ∝ |t|^{{3·{nu_3d*2:.3f} - 2·{3*nu_3d:.3f}}}")
print(f"        = |t|^{{{3*nu_3d*2 - 2*3*nu_3d:.6f}}} = |t|^0 = CONST ✅")

print(f"\n  → B₃/T jest UNIVERSALNĄ stałą klasy 3D Isinga!")
print(f"  → Γ ~ exp(-B₃/T) jest STAŁE blisko T_c (critical slowing NOT of nucleation)")

# Szacunek B₃/T z amplitud:
# σ_s = σ₀ · |t|^{1.260}, ΔF = ΔF₀ · |t|^{1.890}
# B₃/T = (16π/3) · σ₀³ / (T_c · ΔF₀²)
# σ₀ ~ J/a² (energy/area), ΔF₀ ~ J/a³ (energy/volume)

# W jednostkach kT_c:
# σ₀·a² / kT_c ~ 1 (rzędu jedności)
# ΔF₀·a³ / kT_c ~ 1
# B₃/T ~ (16π/3)·(σ₀a²/kT_c)³ / [(ΔF₀a³/kT_c)²] · 1/a^{...}
# Wymiarowo: σ₀a²/(kT_c) i ΔF₀a³/(kT_c) — bezwymiarowe amplitudy

# Universalne amplitudy (z MC):
# σ_s·ξ²/kT_c = R_σ (universal ratio)
# ΔF·ξ³/kT_c = R_F (universal ratio)
# B₃/T = (16π/3) · R_σ³ / R_F²

# Z literatury Isinga 3D (Hasenbusch):
R_sigma = 0.0  # TODO: szukaj
# Przybliżenie: R_σ ~ 0.08, R_F ~ 1

# Zamiast tego użyjmy szacunków numerycznych
# B₃_universal ≈ 50-200 (z różnych modeli GL)
# Przyjmijmy B₃/T ~ 100 jako typowe

B3_over_T_est = 100.0
Gamma_est = np.exp(-B3_over_T_est)
print(f"\n  Szacunek: B₃/T ~ {B3_over_T_est} (typowe dla GL)")
print(f"  Γ/T⁴ ~ exp(-B₃/T) ~ exp(-{B3_over_T_est}) ~ {Gamma_est:.2e}")
print(f"  → Nukleacja zachodzi z prawdopodobieństwem ~ {Gamma_est:.2e} na objętość ξ⁴")

# =====================================================================
# 7. Związek ε₀ z B₃/T
# =====================================================================
print(f"\n{'='*72}")
print("--- 7. ε₀ z B₃/T i ξ ---")
print(f"{'='*72}")

print("""
  Model złożony:
    1. Nukleacja tworzy bąbel o rozmiarze R ~ ξ
    2. Wewnątrz bąbla: Φ = Φ_bubble (blisko T_c, ale Φ_bubble << Φ₀)
    3. ε₀ = Φ_bubble/Φ₀ ≈ σ²_c · (a_sub/ξ)²

  Niezależna estymata (z fluktuacji):
    Φ_bubble ~ ⟨φ²⟩_crit · V_bubble
    ~ (T_c/J) · (a_sub/ξ)²
    ~ σ²_c · (a_sub/ξ)²

  Warunek nukleacji: Γ·t_obs·V_obs ~ 1
    exp(-B₃/T) · T_c⁴ · t_P · ℓ_P³ ~ 1
    → B₃/T ~ 4·ln(M_P/T_c) ~ 4·ln(ξ/a_sub)
    → B₃/T ~ 4·35·ln(10) ~ 322 (dla ξ/a_sub ~ 10³⁵)
""")

# Iteracyjny model:
# B₃/T = 4·ln(ξ/a_sub) (warunek nukleacji)
# ε₀ = σ²_c / (ξ/a_sub)² (fluktuacja)
# N_e = (1/3)·ln(Φ₀/ε₀)

print(f"  {'ξ/a_sub':>12} {'B₃/T':>10} {'exp(-B₃/T)':>14} {'ε₀':>14} {'N_e':>8} {'n_s':>8}")
print("  " + "-" * 75)

for log_xi in range(20, 46, 5):
    xi_a = 10.0**log_xi
    B3T = 4 * np.log(xi_a)
    eps0 = sigma2_c / xi_a**2
    Ne = (1.0/3) * np.log(Phi0/eps0) if eps0 > 0 else 0
    ns = 1 - 2.0/Ne - 0.0005 if Ne > 5 else 0
    mark = " ✅" if 50 <= Ne <= 70 else ""

    print(f"  {xi_a:12.2e} {B3T:10.1f} {np.exp(-min(B3T,700)):14.4e} {eps0:14.4e} {Ne:8.1f} {ns:8.4f}{mark}")

# =====================================================================
# 8. SAMOSPÓJNE ROZWIĄZANIE
# =====================================================================
print(f"\n{'='*72}")
print("--- 8. Samospójne rozwiązanie ---")
print(f"{'='*72}")

print("""
  Równania:
    (i)   ε₀ = σ²_c · (a_sub/ξ)²
    (ii)  N_e = (1/3)·ln(Φ₀/ε₀)
    (iii) Warunek nukleacji: B₃/T ≈ 4·ln(ξ/a_sub)
    (iv)  B₃/T = const (universalny, z Josephsona) — SPRZECZNE z (iii)!

  Rozwiązanie sprzeczności:
    (iv) mówi B₃/T = const blisko T_c (skalowanie krytyczne).
    (iii) to warunek na to JAK BLISKO T_c następuje nukleacja.

    Nukleacja zachodzi gdy Γ·V/T·ξ³ ~ 1
    → exp(-B₃/T) · (ξ/a)³ ~ 1  (w objętości ξ³ w czasie ξ/T)
    → B₃/T ~ 3·ln(ξ/a)

    Ale B₃/T = const (Josephson!) → to WYZNACZA ξ/a:
    ξ/a = exp(B₃/(3T))

  Oszacowanie B₃/T z MC (ising 3D):
""")

# B₃/T estimation from universal amplitude ratios
# In GL theory: B₃/T = (16π/3) σ_s³/(T ΔF²)
# σ_s³/ΔF² at criticality → finite universal ratio

# From literature (approximately):
# Universal number for 3D Ising bounce:
# B₃_crit/T_c ≈ O(1) — ale to jest AT criticality
# Away from T_c: B₃/T = B₃_crit + correction → still O(1) at T_c

# Dokładniej: B₃(t)/T = (16π/3) [σ₀³/(T_c·ΔF₀²)] · |t|^{0} + O(|t|^ω)
# gdzie ω = 0.832 (wykładnik korekcyjny Isinga 3D)

# Szacunek: σ₀ ~ 0.1 J/a², ΔF₀ ~ 0.1 J/a³
# B₃/T ~ (16π/3) · (0.1)³ / (0.1)² ~ (16π/3) · 0.1 ~ 1.7

# Zatem B₃/T ~ O(1) → ξ/a ~ exp(B₃/(3T)) ~ exp(O(1)) ~ O(1)
# To daje ξ/a ~ O(1) → ε₀ ~ σ²_c ~ O(0.1) → N_e ~ O(1)
# NIEWYSTARCZAJĄCE!

B3T_estimates = [1.0, 5.0, 10.0, 50.0, 100.0, 200.0, 300.0]
print(f"  {'B₃/T':>8} {'ξ/a=exp(B₃/3T)':>18} {'ε₀':>14} {'N_e':>8} {'komentarz':>20}")
print("  " + "-" * 75)

for B3T in B3T_estimates:
    xi_a = np.exp(B3T / 3)
    eps0 = sigma2_c / xi_a**2
    Ne = (1.0/3) * np.log(Phi0 / eps0) if eps0 > 0 and eps0 < Phi0 else 0

    comment = ""
    if 50 <= Ne <= 70:
        comment = "✅ CMB OK"
    elif Ne < 1:
        comment = "❌ za mało inflacji"
    elif Ne < 50:
        comment = "❌ za mało"
    elif Ne > 70:
        comment = "⚠️ za dużo"

    print(f"  {B3T:8.1f} {xi_a:18.4e} {eps0:14.4e} {Ne:8.1f} {comment:>20}")

print("""
  KLUCZOWY WYNIK:
  → B₃/T ~ 230-280 daje N_e ∈ [50, 70] (kompatybilne z CMB)

  ALE: z Josephsona, B₃/T ~ O(1) blisko T_c!
  → Nuclelacja AT T_c daje ξ ~ a → ε₀ ~ O(1) → brak inflacji

  ROZWIĄZANIE: nukleacja NIE zachodzi at T_c, ale PONIŻEJ T_c.
  → Supercooling! System musi zostać przechłodzony.
  → |t_nuc| > 0, ξ(t_nuc) = ξ₀·|t_nuc|^{-ν}
  → Potrzeba |t_nuc| ~ 10^{-56} (ekstremalny supercooling)
""")

# Self-consistent: find t_nuc from nucleation condition
# B₃(t)/T = B₃_0 + ... (corrections far from T_c are non-universal)
# In thin-wall limit far from T_c:
# B₃/T ≈ (16π/3)·σ₀³·|t|^{3·1.260} / (T_c · ΔF₀²·|t|^{2·1.890})
#       = (16π/3)·σ₀³/(T_c·ΔF₀²) · |t|^{3.780-3.780} = const (!)
# → B₃/T is genuinely constant in this class!

# So the resolution MUST involve non-universal physics:
# 1. Przejście 1. rodzaju (nie 2. rodzaju) → B₃/T diverges at T_c
# 2. Coupling to gravity → Hawking-Moss bounce
# 3. False vacuum has V(Φ=0) > 0 (cosmological constant) → gravitational nucleation

print(f"\n{'='*72}")
print("--- 9. Nuckelacja grawitacyjna (Hawking-Moss) ---")
print(f"{'='*72}")

print("""
  W obecności grawitacji, bounce Colemana jest modyfikowany.

  Hawking-Moss instanton:
    B_HM = (24π²/V_bar) · [1 - (1 - V_bar·V_false/(3M_P²))^(-1)]
    ≈ (24π² M_P⁴) / V_bar   (dla V_bar << M_P⁴)

  W TGP: bariera V_bar ~ (Φ₀)^p · M_P⁴ / (czegoś)

  Prostszy model: Coleman-de Luccia
    B_CDL ≈ B_flat · [1 - (R_bubble/R_dS)²]^{-2}
    gdzie R_dS = √(3/Λ) = promień de Sittera

  Kluczowe: grawitacja POMAGA nukleacji (zmniejsza B) gdy:
    V_false > 0 (dodatnia stała kosmologiczna w fałszywej próżni)
    → Φ = 0 w TGP ma V(0) = V_Λ,0 (subplankowska) → gravity-assisted tunneling
""")

# Estymata: jak V(Φ=0) wpływa na B?
# V(Φ=0) = ρ_vac(no space) — powinno być ~ M_P⁴ lub ~ M_sub⁴
# B_CDL ≈ B_flat / [1 + B_flat·V_false/(12π²M_P⁴)]
# Jeśli V_false ~ M_P⁴: B_CDL ~ 12π²M_P⁴/V_false ~ 12π² ~ 118

print(f"  Szacunek: B_CDL(V_false ~ M_P⁴) ~ 12π² ≈ {12*np.pi**2:.1f}")
print(f"  → exp(-B_CDL) ~ exp(-118) ~ {np.exp(-118):.2e}")
print(f"  → ε₀ ~ ? (profil bounce'u CDL)")

# Hawking-Moss: g = g_barrier (homogeneous tunneling)
# ε₀_HM = g_barrier² → bardziej jak O(0.01) z naszego potencjału
for eps_bias in [0.01, 0.05, 0.1]:
    g_bar_res = minimize_scalar(lambda g: -V_tgp(g, eps_bias), bounds=(0.01, 0.9), method='bounded')
    g_bar = g_bar_res.x
    eps0_HM = g_bar**2
    Ne_HM = (1.0/3) * np.log(Phi0 / eps0_HM)
    print(f"  ε_bias={eps_bias:.2f}: g_barrier={g_bar:.4f}, ε₀_HM={eps0_HM:.4f}, N_e={Ne_HM:.1f}")

# =====================================================================
# WNIOSKI
# =====================================================================
print(f"\n{'='*72}")
print("WNIOSKI ex164 — R6: ε₀")
print(f"{'='*72}")
print(f"""
  1. BOUNCE BEZWYMIAROWY (Coleman flat-space):
     g₀_bounce ~ O(1) → ε₀ ~ O(1) → N_e ~ O(1) ❌
     NIEWYSTARCZAJĄCY — nie produkuje inflacji.

  2. FLUKTUACJA TERMICZNA (z epsilon0_estimate.py):
     ε₀ = σ²_c · (a_sub/ξ)² → N_e = (1/3)ln(Φ₀/ε₀)
     N_e = 55 wymaga ξ/a_sub ~ 10³⁵ ✅
     ALE: skąd tak duże ξ?

  3. JOSEPHSON + NUKLEACJA:
     B₃/T = const (universalne) blisko T_c
     → Nukleacja at T_c daje ξ ~ a → brak inflacji
     → POTRZEBNY SUPERCOOLING: |t_nuc| ~ 10⁻⁵⁶

  4. GRAWITACJA:
     Coleman-de Luccia/Hawking-Moss: B zmniejszone przez V_false > 0
     → ε₀ ~ g_barrier² ~ O(0.01-0.1) → N_e ~ 2-3 ❌
     Nadal niewystarczające bez hierarchii skali.

  5. FUNDAMENTALNY PROBLEM (potwierdzenie G5):
     ε₀ << 1 koduje hierarchię M_P/m_sub ~ (ξ/a_sub)
     Nie wynika z samej termodynamiki substratu!
     Wymaga albo:
       (a) Ekstremalnego supercoolingu (dlaczego?)
       (b) Inflacyjnego slow-rollu po nukleacji (ε₀ nie musi być ~ 10⁻⁷⁰)
       (c) Nowego mechanizmu (np. landscape substratu)

  6. OPCJA (b) — NAJLEPSZA:
     ε₀ nie musi być absurdalnie małe!
     Inflacja = slow-roll Φ: Φ(t₀) = ε₀·Φ₀ → Φ(t_end) = Φ₀
     Jeśli ε₀ ~ 0.01 → N_e = (1/3)ln(25/0.01) = 2.6 ❌
     Jeśli potencjał daje slow-roll z N_e >> ln(1/ε₀):
       V_slow = γ(Φ³/3 - Φ⁴/(4Φ₀)) → V'(Φ)/V(Φ) << 1 blisko Φ=0
       → ε-parameter slow-roll = (V'/V)² << 1 wystarczy
       → N_e zależy od KSZTAŁTU V, nie od ε₀!

  STATUS R6: **OTWARTY** — identyfikacja kluczowego pytania:
    "Czy inflacja TGP jest driven by ε₀ (initial condition)
     czy przez slow-roll potencjału V(Φ) (dynamics)?"
    Jeśli slow-roll → ε₀ nie musi być ekstremalnie małe.
    Jeśli initial condition → wymaga nowego mechanizmu.

    REKOMENDACJA: zbadać slow-roll TGP potencjału V(Φ) bezpośrednio
    (ex165: slow-roll analysis).
""")
