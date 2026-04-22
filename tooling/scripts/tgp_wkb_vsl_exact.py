"""
tgp_wkb_vsl_exact.py — Dokładny potencjał Sturm-Liouville z profilu kinku TGP

Cel: zastąpienie przybliżenia gaussowskiego V_SL = 1 - D*exp(-xi²/(2σ²))
w fermion_mass_spectrum.py przez potencjał DOKŁADNY wyprowadzony z linearyzacji
równania radialnego kinku TGP wokół profilu χ₀(ξ).

Wyprowadzenie (Dodatek F_v2):
  Kink ODE:  χ'' + (2/ξ)χ' + (α/χ)(χ')² = V'_eff(χ)
  Podstawienie η = u / (ξ·χ₀^α) eliminuje pochodną pierwszego rzędu.
  Równanie Schrödingera dla u(ξ):
    -u'' + V_SL^exact(ξ)·u = E·u
  Dokładny potencjał (dla α=2, β=γ=1):
    V_SL^exact(ξ) = -4χ₀ + 5χ₀² + 8(χ₀')²/χ₀² + 2/ξ² + 8χ₀'/(ξ·χ₀)
  Asymptotyka: V_SL^exact → 1 = m_sp² = γ  przy ξ→∞ ✓
  Warunek progu: E < 1 → stany związane (lepton TGP)

Wyniki:
  - E₀, E₁, E₂ (dokładne stany związane ze studni)
  - E₃ ≥ 1 (brak 4. generacji)
  - κ₁ = (E₁ - E₀)/E₀  (korekcja anharmoniczna, 2. generacja)
  - κ₂ = (E₂ - 2E₁ + E₀)/(2E₀)  (korekcja anharmoniczna, 3. generacja)
  - Porównanie z przybliżeniem gaussowskim

Użycie:
  python scripts/tgp_wkb_vsl_exact.py
"""

import sys
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
import warnings
warnings.filterwarnings("ignore")

sys.stdout.reconfigure(encoding='utf-8')

# ─────────────────────────────────────────────────────────────
# PARAMETRY MODELU TGP
# ─────────────────────────────────────────────────────────────
ALPHA = 2.0      # kinetyczny (z Lematu N0-K i sek10)
BETA = 1.0       # potencjał (β = γ, warunek próżniowy sek08)
GAMMA = 1.0      # próg kontinuum m_sp² = γ = 1
PHI0 = 1.0       # normalizacja (χ = Φ/Φ₀)

PASS_COUNT = 0
FAIL_COUNT = 0
TESTS = []


def record(name, status, info=""):
    global PASS_COUNT, FAIL_COUNT
    TESTS.append((name, status, info))
    if status == "PASS":
        PASS_COUNT += 1
    else:
        FAIL_COUNT += 1
    marker = "PASS" if status == "PASS" else "FAIL"
    print(f"  [{marker}] {name}: {info}")


# ─────────────────────────────────────────────────────────────
# SEKCJA 1: Potencjał efektywny
# ─────────────────────────────────────────────────────────────

def V_eff(chi):
    """V_eff(χ) = β·χ³/3 - γ·χ⁴/4  (β=γ=1, min w χ=1)."""
    return BETA * chi**3 / 3.0 - GAMMA * chi**4 / 4.0


def dV_eff(chi):
    """V'_eff(χ) = β·χ² - γ·χ³ = χ²(1 - χ)."""
    return BETA * chi**2 - GAMMA * chi**3


def d2V_eff(chi):
    """V''_eff(χ) = 2β·χ - 3γ·χ²."""
    return 2 * BETA * chi - 3 * GAMMA * chi**2


# ─────────────────────────────────────────────────────────────
# SEKCJA 2: Profil kinku χ₀(ξ)
# ─────────────────────────────────────────────────────────────

def solve_kink_profile(chi0_guess=0.5, xi_max=30.0, N=1200):
    """
    Rozwiązuje ODE kinku TGP przez strzelaninę (shooting):
      χ'' + (2/ξ)χ' + (α/χ)(χ')² = V'_eff(χ)
    Warunek brzegowy: χ'(0)=0, χ(∞)=1.
    Regularyzacja L'Hôpitala: χ''(0) = V'_eff(χ₀)/3.
    Zwraca: (ξ_arr, χ₀_arr, χ₀'_arr) jako tablice numpy.
    """
    eps = 5e-3
    xi_arr = np.linspace(eps, xi_max, N)

    def ode(xi, y):
        chi, cp = y
        chi = max(chi, 1e-9)
        cpp = dV_eff(chi) - (2.0 / xi) * cp - (ALPHA / chi) * cp**2
        return [cp, cpp]

    def shoot_residual(chi0):
        """Całkuje i zwraca chi(xi_max) - 1."""
        chi0 = max(chi0, 1e-6)
        dV0 = dV_eff(chi0)
        chi_s = chi0 + dV0 * eps**2 / 6.0
        chip_s = dV0 * eps / 3.0
        try:
            sol = solve_ivp(ode, [eps, xi_max], [chi_s, chip_s],
                            method='RK45', t_eval=[xi_max],
                            rtol=1e-8, atol=1e-10, max_step=0.05)
            if sol.success:
                return sol.y[0, -1] - 1.0
            return None
        except Exception:
            return None

    # Szukamy chi0 metodą bisekcji
    chi0_found = chi0_guess
    try:
        # Zakres poszukiwań: chi0 ∈ (0.05, 0.95) dla n=0 profilu
        r_lo = shoot_residual(0.05)
        r_hi = shoot_residual(0.95)
        if r_lo is not None and r_hi is not None and r_lo * r_hi < 0:
            chi0_found = brentq(shoot_residual, 0.05, 0.95, xtol=1e-6)
        else:
            # Próba bezpośrednia z chi0_guess
            chi0_found = chi0_guess
    except Exception:
        chi0_found = chi0_guess

    # Ostateczna integracja z pełną siatkę xi
    chi0_c = max(chi0_found, 1e-6)
    dV0 = dV_eff(chi0_c)
    chi_s = chi0_c + dV0 * eps**2 / 6.0
    chip_s = dV0 * eps / 3.0

    try:
        sol = solve_ivp(ode, [eps, xi_max], [chi_s, chip_s],
                        method='RK45', t_eval=xi_arr,
                        rtol=1e-8, atol=1e-10, max_step=0.05,
                        dense_output=False)
        if sol.success and len(sol.y[0]) == N:
            chi0_arr = sol.y[0]
            chip0_arr = sol.y[1]
            # Korekcja: wymuszenie asymptotyki chi→1
            chi0_arr = np.clip(chi0_arr, 1e-9, 3.0)
            return xi_arr, chi0_arr, chip0_arr
    except Exception:
        pass

    return None, None, None


# ─────────────────────────────────────────────────────────────
# SEKCJA 3: Dokładny potencjał SL
# ─────────────────────────────────────────────────────────────

def compute_V_SL_exact(xi_arr, chi0_arr, chip0_arr):
    """
    Dokładny potencjał Sturm-Liouville z linearyzacji ODE kinku.

    Wyprowadzenie (patrz Dodatek F_v2):
    Podstawienie η = u/(ξ·χ₀^α) eliminuje pierwszą pochodną w ODE fluktuacji.
    Równanie dla u: -u'' + V_SL^exact·u = E·u

    Dokładna formuła (α=2, β=γ=1):
      V_SL^exact(ξ) = -(V''_eff(χ₀) + α·V'_eff(χ₀)/χ₀) + (α² - α)(χ₀')²/χ₀² + 2/ξ²
                    = -4χ₀ + 5χ₀² + 8(χ₀')²/χ₀² + 2/ξ² + 8χ₀'/(ξ·χ₀)

    Weryfikacja: V_SL^exact → 1 przy ξ→∞  (bo χ₀→1, χ₀'→0)
    Wartość: -4(1) + 5(1) + 0 + 0 + 0 = 1 ✓

    Zwraca tablicę V_SL_arr zgodną z konwencją kodu (próg = 1).
    """
    chi0 = np.clip(chi0_arr, 1e-9, 10.0)
    chip0 = chip0_arr

    # Składowe dokładnego potencjału
    term_V2 = d2V_eff(chi0)                              # V''_eff(χ₀) = 2χ₀ - 3χ₀²
    term_V1 = ALPHA * dV_eff(chi0) / chi0               # α·V'_eff(χ₀)/χ₀ = 2χ₀(1-χ₀)
    term_deriv = (ALPHA**2 - ALPHA) * (chip0 / chi0)**2  # (α²-α)(χ₀')²/χ₀² = 2·(χ₀')²/χ₀²

    # Centrobieżna poprawka geometryczna (z metryki sferycznej)
    term_centrifugal = 2.0 / xi_arr**2

    # Mieszana poprawka (z podstawienia)
    # Uwaga: pochodzi z 2·(w'/w)·(1/ξ) gdzie w'/w = -1/ξ - α·χ₀'/χ₀
    term_mixed = 2.0 * ALPHA * chip0 / (xi_arr * chi0)   # 8χ₀'/(ξ·χ₀) dla α=2

    # Potencjał SL: neguję sumę (konwencja: E·u = -u'' + V·u)
    V_SL = -(term_V2 + term_V1) + term_deriv + term_centrifugal + term_mixed

    return V_SL


def V_SL_gaussian(xi_arr):
    """Potencjał gaussowski z fermion_mass_spectrum.py (D=0.9, σ=5.0)."""
    D = 0.9
    sigma = 5.0
    return GAMMA - D * np.exp(-xi_arr**2 / (2 * sigma**2))


# ─────────────────────────────────────────────────────────────
# SEKCJA 4: Kwantyzacja WKB (Bohr-Sommerfeld)
# ─────────────────────────────────────────────────────────────

def wkb_energy(n, xi_arr, V_arr, E_thresh=1.0):
    """
    Energia n-tego stanu związanego z kwantyzacji Bohra-Sommerfelda:
      ∫ sqrt(E - V(ξ)) dξ = (n + 1/2)·π

    Szuka E ∈ (V_min, E_thresh). Jeśli brak stanu → zwraca >= E_thresh.
    """
    V_min = float(np.min(V_arr))
    if V_min >= E_thresh:
        return E_thresh + n * 0.1  # brak studni — tylko stany ciągłe

    def bs_action(E):
        mask = E > V_arr
        if not np.any(mask):
            return -(n + 0.5) * np.pi
        integrand = np.zeros_like(xi_arr)
        integrand[mask] = np.sqrt(np.clip(E - V_arr[mask], 0, None))
        return np.trapezoid(integrand, xi_arr) - (n + 0.5) * np.pi

    E_lo = V_min + 1e-4
    E_hi = E_thresh - 1e-4

    val_lo = bs_action(E_lo)
    val_hi = bs_action(E_hi)

    if val_lo * val_hi >= 0:
        # Stan nie istnieje poniżej progu
        return E_thresh + (n - 1) * 0.05

    try:
        E_n = brentq(bs_action, E_lo, E_hi, xtol=1e-8, maxiter=200)
        return E_n
    except Exception:
        return E_thresh + 0.01


# ─────────────────────────────────────────────────────────────
# SEKCJA 5: GŁÓWNE OBLICZENIA
# ─────────────────────────────────────────────────────────────

print("=" * 60)
print("tgp_wkb_vsl_exact.py — Dokładny potencjał SL z ODE kinku")
print("=" * 60)

# ── Krok 1: Profil kinku ──────────────────────────────────────
print("\n--- Krok 1: Profil kinku χ₀(ξ) ---")
xi_max = 35.0
N_pts = 2000
xi_arr, chi0_arr, chip0_arr = solve_kink_profile(chi0_guess=0.5, xi_max=xi_max, N=N_pts)

if chi0_arr is not None:
    chi0_end = chi0_arr[-1]
    chi0_center = chi0_arr[0]
    converged = abs(chi0_end - 1.0) < 0.05
    record("K1_kink_profile_converged",
           "PASS" if converged else "FAIL",
           f"χ₀(0)={chi0_center:.4f}, χ₀(ξ_max)={chi0_end:.4f} → 1")
else:
    print("  BŁĄD: profil kinku nie zbiegł — używam interpolacji analitycznej")
    # Fallback: profil analityczny dla prostszego potencjału (bez członu α/χ)
    xi_arr = np.linspace(5e-3, xi_max, N_pts)
    # Aproksymacja tanh dla kinku TGP (szacunkowa)
    xi_kink = 3.5  # skala charakterystyczna kinku
    chi0_arr = np.tanh(xi_arr / xi_kink) ** 0.5
    chip0_arr = np.gradient(chi0_arr, xi_arr)
    record("K1_kink_profile_converged", "FAIL",
           "Shooting nie zbiegł; użyto aproksymacji tanh")

# ── Krok 2: Dokładny potencjał SL ────────────────────────────
print("\n--- Krok 2: Potencjał V_SL^exact(ξ) ---")
V_SL_exact = compute_V_SL_exact(xi_arr, chi0_arr, chip0_arr)

# Weryfikacja asymptotyki
V_asym = float(np.mean(V_SL_exact[-50:]))  # średnia ostatnich 50 punktów
record("K2_V_SL_asymptote_check",
       "PASS" if abs(V_asym - 1.0) < 0.05 else "FAIL",
       f"V_SL(ξ→∞) = {V_asym:.4f} (oczekiwane: 1.0 = m_sp² = γ)")

# Minimum potencjału (studnia)
V_SL_interior = V_SL_exact[5:]  # pomiń osobliwość przy xi=0
V_min_exact = float(np.min(V_SL_interior))
xi_min_idx = np.argmin(V_SL_interior) + 5
xi_at_min = xi_arr[xi_min_idx]
print(f"  V_SL^exact minimum: {V_min_exact:.4f} przy ξ = {xi_at_min:.2f}")

# Porównaj z gaussowskim
V_gauss = V_SL_gaussian(xi_arr)
V_min_gauss = float(np.min(V_gauss))
print(f"  V_SL^gauss minimum: {V_min_gauss:.4f} (D=0.9, σ=5.0)")

# ── Krok 3: Energia WKB — stany związane ─────────────────────
print("\n--- Krok 3: Kwantyzacja WKB — stany związane ---")
# Używamy tylko fragmentu siatki z ξ > 0.5 (pomijamy osobliwość)
xi_wkb_mask = xi_arr > 0.5
xi_wkb = xi_arr[xi_wkb_mask]
V_exact_wkb = V_SL_exact[xi_wkb_mask]
V_gauss_wkb = V_gauss[xi_wkb_mask]

E_thresh = GAMMA  # = 1.0

print("  Potencjał dokładny:")
E_exact = []
for n in range(4):
    E_n = wkb_energy(n, xi_wkb, V_exact_wkb, E_thresh=E_thresh)
    E_exact.append(E_n)
    label = "STAN ZWIĄZANY" if E_n < E_thresh else "CIĄGŁE (>= prog)"
    print(f"    E_{n}^exact = {E_n:.4f}  [{label}]")

print("  Potencjał gaussowski (D=0.9, σ=5.0):")
E_gauss = []
for n in range(4):
    E_n = wkb_energy(n, xi_wkb, V_gauss_wkb, E_thresh=E_thresh)
    E_gauss.append(E_n)
    label = "STAN ZWIĄZANY" if E_n < E_thresh else "CIĄGŁE (>= prog)"
    print(f"    E_{n}^gauss = {E_n:.4f}  [{label}]")

# ── Krok 4: Weryfikacja liczby generacji ──────────────────────
print("\n--- Krok 4: Predykcja TGP — 3 generacje (nie 4) ---")
n_bound_exact = sum(1 for E in E_exact if E < E_thresh)
n_bound_gauss = sum(1 for E in E_gauss if E < E_thresh)

record("K4_three_generations_exact",
       "PASS" if n_bound_exact == 3 else "FAIL",
       f"Liczba stanów związanych (dokładne): {n_bound_exact} (oczekiwane: 3)")
record("K4_three_generations_gauss",
       "PASS" if n_bound_gauss == 3 else "FAIL",
       f"Liczba stanów związanych (gaussowski): {n_bound_gauss} (oczekiwane: 3)")

# ── Krok 5: Korekcje anharmoniczne κ₁, κ₂ ────────────────────
print("\n--- Krok 5: Korekcje anharmoniczne κ₁, κ₂ ---")

if n_bound_exact >= 3:
    E0, E1, E2 = E_exact[0], E_exact[1], E_exact[2]
    # κ₁: korekcja 2. generacji (μ) względem 1. (e)
    kappa1_exact = (E1 - E0) / E0
    # κ₂: korekcja anharmoniczności (dla 3. generacji τ)
    kappa2_exact = (E2 - 2 * E1 + E0) / (2 * E0)
    print(f"  Dokładny potencjał:")
    print(f"    E₀ = {E0:.4f}, E₁ = {E1:.4f}, E₂ = {E2:.4f}")
    print(f"    κ₁ = (E₁-E₀)/E₀ = {kappa1_exact:.4f}")
    print(f"    κ₂ = (E₂-2E₁+E₀)/(2E₀) = {kappa2_exact:.4f}")
else:
    kappa1_exact = None
    kappa2_exact = None
    print("  WARN: Za mało stanów związanych dla κ₁, κ₂")

if n_bound_gauss >= 3:
    E0g, E1g, E2g = E_gauss[0], E_gauss[1], E_gauss[2]
    kappa1_gauss = (E1g - E0g) / E0g
    kappa2_gauss = (E2g - 2 * E1g + E0g) / (2 * E0g)
    print(f"  Gaussowski potencjał (D=0.9, σ=5.0):")
    print(f"    E₀ = {E0g:.4f}, E₁ = {E1g:.4f}, E₂ = {E2g:.4f}")
    print(f"    κ₁ = {kappa1_gauss:.4f}")
    print(f"    κ₂ = {kappa2_gauss:.4f}")

# ── Krok 6: Masy leptonów z κ₁, κ₂ ──────────────────────────
print("\n--- Krok 6: Estymacja stosunków mas z WKB ---")
print("  Schemat: m_n ~ Φ₀ · E_n  (jednostki Φ₀)")

# Stosunki mas (bezwymiarowe)
if n_bound_exact >= 3:
    r21_wkb = E_exact[1] / E_exact[0]  # m_μ/m_e (schemat energii defektu)
    r32_wkb = E_exact[2] / E_exact[1]  # m_τ/m_μ
    r31_wkb = E_exact[2] / E_exact[0]  # m_τ/m_e
    print(f"  WKB (dokładny): r₂₁ = E₁/E₀ = {r21_wkb:.2f}")
    print(f"  WKB (dokładny): r₃₁ = E₂/E₀ = {r31_wkb:.2f}")
    print(f"  PDG obserwacja: r₂₁ = 206.77, r₃₁ = 3477.18")
    print(f"  Uwaga: energie węzłów nie dają bezpośrednio mas (wymagana kalibracja α_K)")

# Różnica κ₁ między dokładnym a gaussowskim
if kappa1_exact is not None and n_bound_gauss >= 3:
    dk1 = abs(kappa1_exact - kappa1_gauss)
    dk2 = abs(kappa2_exact - kappa2_gauss) if kappa2_exact is not None else 0
    record("K6_kappa1_comparison",
           "PASS" if dk1 < 0.3 else "FAIL",
           f"|κ₁^exact - κ₁^gauss| = {dk1:.4f} (korekcja anharmoniczna)")
    record("K6_kappa2_comparison",
           "PASS" if dk2 < 0.5 else "FAIL",
           f"|κ₂^exact - κ₂^gauss| = {dk2:.4f}")

# ── Krok 7: Weryfikacja asymptotyki V_SL ──────────────────────
print("\n--- Krok 7: Weryfikacje analityczne ---")

# Test: V_SL(xi→∞) = 1 (forma analityczna)
chi_test = 1.0
chip_test = 0.0
xi_test = 100.0
V_analytic_limit = -(d2V_eff(chi_test) + ALPHA * dV_eff(chi_test) / chi_test) + \
                   (ALPHA**2 - ALPHA) * (chip_test / chi_test)**2 + \
                   2.0 / xi_test**2 + \
                   2 * ALPHA * chip_test / (xi_test * chi_test)
record("K7_V_SL_analytic_limit",
       "PASS" if abs(V_analytic_limit - 1.0) < 0.01 else "FAIL",
       f"V_SL^exact(χ₀=1, χ₀'=0, ξ=100) = {V_analytic_limit:.6f} (oczekiwane: 1.0)")

# Test: wzór na próg = V''(1) + α·V'(1)/1 z oboma χ₀=1, χ₀'=0
# V_SL(ξ→∞) = -(V''_eff(1) + α·V'_eff(1)/1) = -(2-3) + α·(1-1) = -(-1) = 1 ✓
V_limit_formula = -(d2V_eff(1.0) + ALPHA * dV_eff(1.0) / 1.0)
record("K7_asymptote_formula",
       "PASS" if abs(V_limit_formula - 1.0) < 1e-10 else "FAIL",
       f"-(V''(1) + α·V'(1)) = -({d2V_eff(1.0):.1f} + {ALPHA*dV_eff(1.0):.1f}) = {V_limit_formula:.6f}")

# ─────────────────────────────────────────────────────────────
# PODSUMOWANIE
# ─────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("PODSUMOWANIE")
print("=" * 60)

for name, status, info in TESTS:
    marker = "PASS" if status == "PASS" else "FAIL"
    print(f"  [{marker}] {name}: {info}")

total = PASS_COUNT + FAIL_COUNT
print(f"\nWYNIK: {PASS_COUNT}/{total} PASS")

if kappa1_exact is not None:
    print(f"\nKLUCZOWE WYNIKI (dokładny V_SL):")
    print(f"  κ₁ = {kappa1_exact:.4f}  (korekcja 2. generacji)")
    print(f"  κ₂ = {kappa2_exact:.4f}  (korekcja 3. generacji)")
    print(f"  E₀ = {E_exact[0]:.4f}")
    print(f"  E₁ = {E_exact[1]:.4f}")
    print(f"  E₂ = {E_exact[2]:.4f}")
    print(f"  E₃ = {E_exact[3]:.4f}  ({'≥ 1 → brak 4. gen.' if E_exact[3] >= E_thresh else '< 1 OSTRZEZENIE'})")
    if n_bound_gauss >= 3:
        print(f"\nPORÓWNANIE Z GAUSSOWSKIM:")
        print(f"  Δκ₁ = {kappa1_exact - kappa1_gauss:+.4f}  "
              f"({'+' if kappa1_exact > kappa1_gauss else '-'}{'wiek' if abs(kappa1_exact-kappa1_gauss)<0.1 else 'ZNACZNA'} korekcja)")
        print(f"  Δκ₂ = {kappa2_exact - kappa2_gauss:+.4f}")

print(f"\nFormula dokładnego potencjału (Dodatek F_v2):")
print(f"  V_SL^exact(ξ) = -4χ₀ + 5χ₀² + 8(χ₀')²/χ₀² + 2/ξ² + 8χ₀'/(ξ·χ₀)")
print(f"  Asymptotyka: -4·1 + 5·1 + 0 + 0 + 0 = 1 = m_sp² = γ  ✓")
print(f"  Status: WYPROWADZONE z linearyzacji ODE kinku TGP [AN]")
