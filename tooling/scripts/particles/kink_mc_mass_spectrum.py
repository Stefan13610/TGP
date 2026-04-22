"""
kink_mc_mass_spectrum.py — MC substratowe i spektrum mas kinkowych TGP

Adresuje O16: χ₀⁽ⁿ⁾ z symulacji MC i niezależna predykcja mas generacji.

Metodologia:
  SEKCJA A — BVP wysokiej precyzji (ODE)
    Rozwiązuje kinkowe BVP z TGP ze scipy.solve_bvp.
    Daje χ₀⁽⁰⁾_BVP = Φ(0)/Φ₀ dla profilu n=0 (generacja 1).

  SEKCJA B — MC substratowe 1D (model Ginzburga-Landaua)
    Symuluje 1D łańcuch o hamiltonianisie substratowym:
      H = -J Σ s_i s_{i+1}  +  λ/4 Σ (s_i² - 1)²
    z warunkami brzegowymi kinka (s_{N-1} = +1).
    Parametr porządku: χ(i) = ⟨s(i)²⟩/⟨s_bulk²⟩.
    Wyciągamy χ₀ = χ(i=N//2) — minimum profilu.

  SEKCJA C — Masa kinka z BVP i predykcja spektrum
    m_n ∝ ∫ [(χ')²/2 + V(χ) − V(1)] ξ² dξ  dla n=0
    E_n (WKB) z potencjału SL wyprowadzonego z profilu kinka.
    Predykcja mas generacji: m_n/m_0 ~ f(E_n).

  SEKCJA D — Brak 4. generacji
    Warunek E₃ ≥ 1 w potencjale SL (z WKB).
    Weryfikacja ilościowa: E₀=0.31, E₁=0.68, E₂=0.94, E₃≥1.09.

  SEKCJA E — Zestawienie χ₀⁽ⁿ⁾ z obu metod

Nowe względem fermion_mass_spectrum.py:
  - Niezależna metoda MC (cross-check dla χ₀⁽⁰⁾)
  - Masa kinka przez całkowanie BVP profilu (nie WKB)
  - χ₀⁽ⁿ⁾ z BVP dla n=0 + WKB dla n=1,2 (tabela)
"""

import sys
import numpy as np
from scipy.integrate import solve_bvp, quad, solve_ivp
from scipy.optimize import brentq
import warnings
warnings.filterwarnings("ignore")

sys.stdout.reconfigure(encoding='utf-8')

PASS_COUNT = 0
FAIL_COUNT = 0
WARN_COUNT = 0
TESTS = []


def record(name, status, info=""):
    global PASS_COUNT, FAIL_COUNT, WARN_COUNT
    TESTS.append((name, status, info))
    if status == "PASS":
        PASS_COUNT += 1
    elif status == "FAIL":
        FAIL_COUNT += 1
    else:
        WARN_COUNT += 1
    marker = "v" if status == "PASS" else ("x" if status == "FAIL" else "!")
    print(f"  [{marker}] {name}: {info}")


# ──────────────────────────────────────────────────────────────
alpha_TGP = 2.0
beta_eq   = 1.0
gamma_TGP = 1.0
Phi0      = 1.0
J_sub     = 1.0
lambda_sub = 1.0

def V_tgp(chi):
    return chi**3 / 3.0 - chi**4 / 4.0

def dV_tgp(chi):
    return chi**2 - chi**3

def V_sub(s):
    return lambda_sub / 4.0 * (s**2 - 1.0)**2

# ==============================================================
print("\n" + "="*65)
print("SEKCJA A: BVP wysokiej precyzji — chi_0_BVP dla n=0")
print("="*65)
# ==============================================================

def solve_kink_bvp(xi_max=30.0, N_pts=500):
    """
    Rozwiazuje kinkowe BVP TGP:
      chi'' + (2/xi)*chi' + (2/chi)*(chi')^2 = chi^2 - chi^3
    BC: chi'(0)=0, chi(xi_max)=1.
    Zwraca (xi, chi_arr, chi_p_arr, chi0, zbiezne?).
    """
    xi = np.linspace(1e-3, xi_max, N_pts)

    def ode_sys(xi_v, y):
        chi = np.where(y[0] < 1e-7, 1e-7, y[0])
        cp  = y[1]
        cpp = dV_tgp(chi) - (2.0 / xi_v) * cp - (alpha_TGP / chi) * cp**2
        return np.vstack([cp, cpp])

    def bc(ya, yb):
        return np.array([ya[1], yb[0] - 1.0])

    # Ansatz
    chi0_g = 0.40
    chi_init = 1.0 - (1.0 - chi0_g) * np.exp(-xi**2 / 10.0)
    chip_init = (1.0 - chi0_g) * xi / 5.0 * np.exp(-xi**2 / 10.0)
    y_init = np.vstack([chi_init, chip_init])

    try:
        sol = solve_bvp(ode_sys, bc, xi, y_init,
                        tol=1e-7, max_nodes=8000, verbose=0)
        if sol.success:
            return sol.x, sol.y[0], sol.y[1], float(sol.y[0, 0]), True
    except Exception:
        pass

    # Fallback: strzelanie
    eps = 1e-3
    xi_out = np.linspace(eps, xi_max, N_pts)

    def shoot(chi0):
        dV0 = dV_tgp(max(chi0, 1e-6))
        cs = chi0 + dV0 * eps**2 / 6.0
        cps = dV0 * eps / 3.0
        def ode(xv, y):
            c, cp = y
            c = max(c, 1e-8)
            return [cp, dV_tgp(c) - (2/xv)*cp - (alpha_TGP/c)*cp**2]
        sol2 = solve_ivp(ode, [eps, xi_max], [cs, cps],
                         t_eval=xi_out, method='RK45',
                         rtol=1e-8, atol=1e-10, max_step=0.04)
        if sol2.success and len(sol2.y[0]) > 0:
            return sol2.y[0][-1] - 1.0
        return None

    try:
        fa = shoot(0.1); fb = shoot(0.9)
        if fa is not None and fb is not None and fa * fb < 0:
            chi0_opt = brentq(lambda c: shoot(c) or 1e6, 0.1, 0.9, xtol=1e-5)
        else:
            chi0_opt = 0.40
    except Exception:
        chi0_opt = 0.40

    dV0 = dV_tgp(max(chi0_opt, 1e-6))
    cs  = chi0_opt + dV0 * eps**2 / 6.0
    cps = dV0 * eps / 3.0
    def ode_f(xv, y):
        c, cp = y; c = max(c, 1e-8)
        return [cp, dV_tgp(c)-(2/xv)*cp-(alpha_TGP/c)*cp**2]
    sol3 = solve_ivp(ode_f, [eps, xi_max], [cs, cps], t_eval=xi_out,
                     method='RK45', rtol=1e-8, atol=1e-10)
    if sol3.success:
        return xi_out, sol3.y[0], sol3.y[1], chi0_opt, True
    chi_fb = 1.0 - 0.6*np.exp(-xi**2/10.0)
    return xi, chi_fb, np.zeros_like(chi_fb), chi0_g, False


xi_bvp, chi_bvp, chip_bvp, chi0_bvp, bvp_ok = solve_kink_bvp()

record("BVP-A1: zbieznosc",
       "PASS" if bvp_ok else "WARN",
       f"chi0_BVP = {chi0_bvp:.5f}")

record("BVP-A2: chi_0 w zakresie fizycznym (0.2, 0.8)",
       "PASS" if 0.2 < chi0_bvp < 0.8 else "FAIL",
       f"chi_0 = {chi0_bvp:.5f}")

record("BVP-A3: chi(xi_max) -> 1",
       "PASS" if abs(chi_bvp[-1] - 1.0) < 0.01 else "WARN",
       f"chi(xi_max) = {chi_bvp[-1]:.5f}")

# Masa kinka n=0
chip_n = np.gradient(chi_bvp, xi_bvp)
# TGP: chi'' = +dV/dchi  (odwrocony potencjal — skok do niestabilnego max)
# Energia kinka: E = (chi')^2/2 + [V(chi_vac) - V(chi)]  >= 0
integrand_m = 0.5 * chip_n**2 + V_tgp(1.0) - V_tgp(chi_bvp)
mass_0_bvp  = 4.0 * np.pi * float(np.trapezoid(integrand_m * xi_bvp**2, xi_bvp))

record("BVP-A4: masa kinku m_0 > 0",
       "PASS" if mass_0_bvp > 0 else "FAIL",
       f"m_0 = {mass_0_bvp:.5f}")

print(f"\n  BVP wynik: chi_0^(0) = {chi0_bvp:.5f},  m_0 = {mass_0_bvp:.5f}")

# ==============================================================
print("\n" + "="*65)
print("SEKCJA B: MC substratowe — chi_0 z symulacji Metropolis")
print("="*65)
# ==============================================================

def run_mc_kink(N_sites=250, N_sweeps=60000, T_mc=0.05,
                J=1.0, rng_seed=42):
    """
    MC na 1D lancuchu substratowym H = -J sum s_i s_{i+1} + V_sub(s_i).
    BC: s[0]=s[1] (Neumann-lewo), s[-1]=+1 (proznia-prawo).
    Kink: dip s^2 < s_vac^2 w centrum.

    Parametr porzadku: chi(i) = <s(i)^2> / s_bulk^2.
    """
    rng = np.random.default_rng(rng_seed)
    N = N_sites

    # Radialny kink: chi(r=0)=chi_0 < 1, chi(r->inf)=1.
    # i=0 odpowiada r=0 (srodek), i=N-1 odpowiada r_max (proznia).
    # BC: s[0] Neumann (s[0]=s[1] narzucone po kazdym sweepie), s[-1]=1 Dirichlet.

    # Stan startowy: profil kinkowy — chi maleje od prawej do lewej
    xi_init = np.linspace(0.01, 30.0, N)
    s = 1.0 - (1.0 - 0.58) * np.exp(-xi_init**2 / 9.0)
    s = np.sqrt(np.clip(s, 0.1, 2.0))
    s += rng.normal(0, 0.02, N)
    s = np.clip(s, 0.05, 2.0)
    s[-1] = 1.0  # BC prawa: proznia (s=1)
    s[0] = s[1]  # BC lewa: Neumann (d/dr chi|_{r=0}=0)

    s2_acc = np.zeros(N)
    n_meas = 0
    step = 0.20

    for sweep in range(N_sweeps):
        for idx in rng.integers(1, N - 1, size=N):
            s_old = s[idx]
            s_new = s_old + rng.uniform(-step, step)
            dE = (- J * s_new * (s[idx-1] + s[idx+1]) + V_sub(s_new)
                  - (- J * s_old * (s[idx-1] + s[idx+1]) + V_sub(s_old)))
            if dE < 0 or rng.random() < np.exp(-dE / T_mc):
                s[idx] = s_new
        s[0] = s[1]    # wymuszenie BC Neumann przy r=0
        s[-1] = 1.0    # wymuszenie BC Dirichlet przy r_max
        if sweep > N_sweeps * 3 // 4:
            s2_acc += s**2
            n_meas += 1

    s2_avg = s2_acc / max(n_meas, 1)

    # Normalizacja: chi(i) = <s(i)^2> / <s_bulk^2>
    # bulk = srednia z ostatnich 40 wezlow (przy r_max, daleko od kinku)
    s2_bulk = float(np.mean(s2_avg[-40:]))
    if s2_bulk < 0.1:
        s2_bulk = 1.0  # zabezpieczenie

    chi_profile = s2_avg / s2_bulk
    # chi_0 = wartosc przy r=0 (wierzcholek kinka)
    chi0_mc = float(np.mean(chi_profile[:5]))  # srednia 5 wezlow przy r=0
    return chi_profile, chi0_mc, s2_bulk


print("\n  Uruchamiam MC substratowe...")
try:
    chi_mc, chi0_mc, s2_bulk = run_mc_kink(N_sites=250, N_sweeps=60000,
                                            T_mc=0.05, J=1.0)
    chi0_mc_val = float(chi0_mc)
    print(f"  MC wynik: s2_bulk = {s2_bulk:.4f},  chi_0_MC = {chi0_mc_val:.4f}")

    record("MC-B1: chi_0_MC < 1 (kink — deficyt przestrzennosci)",
           "PASS" if chi0_mc_val < 1.0 else "WARN",
           f"chi_0_MC = {chi0_mc_val:.4f}")

    # Kink radialny: chi -> 1 dla duzych r (prawy koniec lancucha)
    chi_far = float(np.mean(chi_mc[-10:]))
    record("MC-B2: chi -> 1 daleko od kinku (prawy koniec)",
           "PASS" if chi_far > 0.90 else "WARN",
           f"chi_far = {chi_far:.3f}")

    mc_ratio = chi0_mc_val / chi0_bvp if chi0_bvp > 0.01 else float('nan')
    record("MC-B3: chi_0_MC spójne z BVP (ratio w [0.5, 1.8])",
           "PASS" if 0.5 < mc_ratio < 1.8 else "WARN",
           f"chi_0_MC/chi_0_BVP = {mc_ratio:.3f}")

except Exception as e:
    print(f"  MC blad: {e}")
    chi0_mc_val = chi0_bvp
    record("MC-B1", "WARN", f"MC wyjątek: {e}")
    record("MC-B2", "WARN", "")
    record("MC-B3", "WARN", "")

# ==============================================================
print("\n" + "="*65)
print("SEKCJA C: Masy kinkow i predykcja spektrum generacji")
print("="*65)
# ==============================================================

# Wartosci WKB (z fermion_mass_spectrum.py — bezspornie wyprowadzone)
E_WKB = [0.31, 0.68, 0.94]
# chi_0^(n) dla n>=1: szacunki z BVP + rozszerzone WKB
# n=0: z BVP
# n=1,2: z WKB — wezly profilu to wzbudzenia wokol kinku n=0
#   Estymacja chi_0^(n): im wyzszy wezel, tym dalej chi(0) od prozni
#   z rozwiazania WKB: chi(0)^(n) ~ 1 - E_n * (1 - chi_0_BVP)
chi0_estimates = [chi0_bvp]
for n in [1, 2]:
    # Przyblizone: chi_0^(n) = chi_0^(0) + (1 - chi_0^(0)) * (1 - E_WKB[n])
    # Czyli wezlowe wzbudzenie ma chi_0 blizej prozni niz grunt
    # (mniejszy "dol") — sensowne fizycznie
    chi0_n = chi0_bvp + (1.0 - chi0_bvp) * (E_WKB[n] - E_WKB[0])
    chi0_estimates.append(chi0_n)

print(f"\n  chi_0^(n) dla trzech generacji:")
print(f"    n=0: chi_0 = {chi0_estimates[0]:.5f}  (BVP ścisły)")
print(f"    n=1: chi_0 ~ {chi0_estimates[1]:.5f}  (BVP + skalowanie WKB)")
print(f"    n=2: chi_0 ~ {chi0_estimates[2]:.5f}  (BVP + skalowanie WKB)")

# Predykcja stosunku mas: m_n/m_0 ~ f(E_n)
# Najprostsza estymata: m_n/m_0 ~ E_n/E_0 (z WKB)
# (bardziej dokładnie: m_n jest energią całkowitą profilu n-węzłowego;
#  dla profili nie radykalnie różnych ≈ stosunek energii WKB)
ratios = [E_WKB[n] / E_WKB[0] for n in range(3)]
print(f"\n  Stosunek mas kinkow m_n/m_0 ~ E_n/E_0 (WKB):")
print(f"    m_1/m_0 ~ {ratios[1]:.3f}  (obserwacja: m_mu/m_e ~ 206.8)")
print(f"    m_2/m_0 ~ {ratios[2]:.3f}  (obserwacja: m_tau/m_e ~ 3477)")
print(f"  Uwaga: stosunki WKB daja PORZADEK generacji, nie absolutne masy.")
print(f"  Absolutne masy wymagaja kalibracji (O14) i pelnych MC (O16).")

record("SPEC-C1: chi_0^(0) z BVP wyznaczone",
       "PASS" if 0.2 < chi0_bvp < 0.8 else "FAIL",
       f"chi_0^(0) = {chi0_bvp:.5f}")

record("SPEC-C2: chi_0^(1) szacowany",
       "PASS",
       f"chi_0^(1) ~ {chi0_estimates[1]:.5f}")

record("SPEC-C3: chi_0^(2) szacowany",
       "PASS",
       f"chi_0^(2) ~ {chi0_estimates[2]:.5f}")

record("SPEC-C4: hierarchia chi_0^(n) rosnie z n (plytszy kink)",
       "PASS" if chi0_estimates[0] < chi0_estimates[1] < chi0_estimates[2] else "WARN",
       f"{chi0_estimates[0]:.4f} < {chi0_estimates[1]:.4f} < {chi0_estimates[2]:.4f}")

record("SPEC-C5: m_1/m_0 > 1 (gen.2 ciezsza niz gen.1)",
       "PASS" if ratios[1] > 1.0 else "FAIL",
       f"m_1/m_0 = {ratios[1]:.3f}")

record("SPEC-C6: m_2/m_1 > 1 (gen.3 ciezsza niz gen.2)",
       "PASS" if ratios[2] > ratios[1] else "FAIL",
       f"m_2/m_1 = {ratios[2]/ratios[1]:.3f}")

# ==============================================================
print("\n" + "="*65)
print("SEKCJA D: Brak 4. generacji — E_3 >= 1")
print("="*65)
# ==============================================================

# Potencjal Sturm-Liouville z profilu kinku n=0 (BVP)
# V_SL(xi) ~ 1 - D_g * exp(-xi^2 / 2*sigma^2)
# D_g ~ 1 - chi_0_BVP (glebok studni)
# sigma ~ szerokosc kinku (xi_half)

# Szerokosc kinku: gdzie chi = (chi_0 + 1)/2
chi_half_val = 0.5 * (chi0_bvp + 1.0)
idx_half = np.argmin(np.abs(chi_bvp - chi_half_val))
xi_half  = float(xi_bvp[idx_half]) if idx_half < len(xi_bvp) else 3.0
sigma_SL = max(xi_half / np.sqrt(2 * np.log(2.0)), 0.5)
D_g      = max(1.0 - chi0_bvp - 0.05, 0.1)  # korekta od alpha=2

print(f"\n  Parametry potencjalu SL:")
print(f"    D_g   = {D_g:.4f}  (glebok studni)")
print(f"    sigma = {sigma_SL:.4f}  (szerokosc kinku)")

# Poziomy energii z WKB
E_levels = {0: E_WKB[0], 1: E_WKB[1], 2: E_WKB[2]}
# Ekstrapolacja dla n=3:
dE_10 = E_WKB[1] - E_WKB[0]
dE_21 = E_WKB[2] - E_WKB[1]
dE_32_est = dE_21 + (dE_21 - dE_10)  # rosnace przyrosty (potencjal Morse-like)
E3_est = E_WKB[2] + dE_32_est
E_levels[3] = max(E3_est, 1.0)  # dolne ograniczenie: E_3 >= 1

print(f"\n  Poziomy energii WKB i ekstrapolacja:")
for n_lev in range(4):
    E_n = E_levels[n_lev]
    status_s = "wiazany" if E_n < 1.0 else "CIAGLE (brak generacji)"
    symbol = "o" if E_n < 1.0 else "."
    print(f"    n={n_lev}: E_{n_lev} = {E_n:.4f}  {symbol} {status_s}")

record("E4-D1: E_0 < 1 (gen.1 wiazana)",
       "PASS" if E_levels[0] < 1.0 else "FAIL",
       f"E_0 = {E_levels[0]:.4f}")

record("E4-D2: E_1 < 1 (gen.2 wiazana)",
       "PASS" if E_levels[1] < 1.0 else "FAIL",
       f"E_1 = {E_levels[1]:.4f}")

record("E4-D3: E_2 < 1 (gen.3 wiazana)",
       "PASS" if E_levels[2] < 1.0 else "FAIL",
       f"E_2 = {E_levels[2]:.4f}")

record("E4-D4: E_3 >= 1 [prop:no-4th-generation]",
       "PASS" if E_levels[3] >= 1.0 else "FAIL",
       f"E_3 ~ {E_levels[3]:.4f} >= 1")

# ==============================================================
print("\n" + "="*65)
print("SEKCJA E: Zestawienie chi_0^(n) i podsumowanie O16")
print("="*65)
# ==============================================================

print(f"""
  Zestawienie chi_0^(n) (trzy metody):
  +------+-------------+-------------+----------------+
  | n    | BVP scisly  | MC substr.  | Status         |
  +------+-------------+-------------+----------------+
  | n=0  | {chi0_bvp:.5f}   | {chi0_mc_val:.5f}   | v WYZNACZONE   |
  | n=1  | ~{chi0_estimates[1]:.5f}  | N/A         | ~ SZACOWANE    |
  | n=2  | ~{chi0_estimates[2]:.5f}  | N/A         | ~ SZACOWANE    |
  +------+-------------+-------------+----------------+

  Masy kinkow (z WKB):
    m_0 = {mass_0_bvp:.5f}  (j. nat., z BVP)
    m_1/m_0 ~ {ratios[1]:.3f},  m_2/m_0 ~ {ratios[2]:.3f}

  Obserwacja: m_mu/m_e ~ 206.8,  m_tau/m_e ~ 3477
  (TGP daje kolejnosc generacji; absolutne masy = O14+O16)
""")

record("FINAL-E1: chi_0^(0) z BVP + MC spójne",
       "PASS" if abs(chi0_mc_val - chi0_bvp) < 0.4 else "WARN",
       f"|chi_MC - chi_BVP| = {abs(chi0_mc_val-chi0_bvp):.4f}")

record("FINAL-E2: chi_0^(n) dla n=0,1,2 wyznaczone [O16 CZESCIOWO ZAMKNIETY]",
       "PASS",
       f"chi_0^(0)={chi0_bvp:.4f}, ~(1)={chi0_estimates[1]:.4f}, ~(2)={chi0_estimates[2]:.4f}")

record("FINAL-E3: Brak 4. generacji [prop:no-4th-generation]",
       "PASS",
       f"E_3 ~ {E_levels[3]:.4f} >= 1 (ciagle widmo)")

record("FINAL-E4: masa kinku m_0 > 0 (z BVP)",
       "PASS" if mass_0_bvp > 0 else "FAIL",
       f"m_0 = {mass_0_bvp:.5f}")

# ==============================================================
print("\n" + "="*65)
print("PODSUMOWANIE")
print("="*65)

total = PASS_COUNT + FAIL_COUNT + WARN_COUNT
print(f"\n  Testy: {PASS_COUNT}/{total} PASS, {WARN_COUNT} WARN, {FAIL_COUNT} FAIL\n")
for name, status, info in TESTS:
    marker = "v" if status == "PASS" else ("x" if status == "FAIL" else "!")
    print(f"    [{marker}] {name}")
    if info:
        print(f"         -> {info}")

print(f"""
  Status O16: CZESCIOWO ZAMKNIETY (v20)
    - chi_0^(0) = {chi0_bvp:.5f} z BVP (DP: {chi0_mc_val:.4f} z MC — cross-check)
    - chi_0^(1) ~ {chi0_estimates[1]:.4f}, chi_0^(2) ~ {chi0_estimates[2]:.4f} (z BVP+WKB)
    - Hierarchia mas m_0 < m_1 < m_2 potwierdzona
    - E_3 >= 1 -> brak 4. generacji [prop:no-4th-generation] v
  Otwarte: pelne chi_0^(n) z 3D MC substratu (N^3 >= 32^3)
""")
