#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r3_p_alpha_analytical.py
=========================

PURPOSE
-------
Probowac wyprowadzic analitycznie wzor p(alpha) = 5 - alpha dla R3 mass formula
m_obs = c * A_tail^p(alpha).

KONTEKST
--------
Z r3_observable_vs_full_mass.py odkrylo sie empirycznie:
  alpha=1.00: p=4.001  (= 5-1 dokladnie)
  alpha=2.00: p=3.001  (= 5-2 dokladnie)
  alpha=1.50: p=3.428  (vs 5-1.5=3.5, diff 2.1%)

Liniowy wzor w polowie rozsadnym zakresie alpha. Pytanie: czy wynika z fizyki?

HIPOTEZA ANALITYCZNA
--------------------
Dla solitonu w d=3 z asymptotyka tail:
  g(r) - 1 ~= A * sin(r+delta) / r   (linearyzacja wokol vacuum g=1)

Kinetic density:
  T(r) = (1/2) g^{2*alpha} (g')^2

Asymptotic (g -> 1):
  g^{2*alpha} -> 1     (bezalfowe — nie powinno wnosic zaleznosci od alpha)
  (g')^2 ~= (A/r)^2    (dla cos~1 srednia)

ALE rdzen (r ~ A^{1/?}):
  g(0) = g0 (a nie 1), wiec g^{2*alpha} != 1 w rdzeniu
  Dla excess solitonow g0 > 1, g rozni sie od 1 nawet w rdzeniu

Wnioski wstepne (do testu numerycznego):

  K_total = K_tail + K_core
  K_tail ~ A^2 * R_max / 4         [universal, nie zalezy od alpha]
  K_core ~ A^? * <g^{2*alpha}>_core  [tu wchodzi alpha!]

  Dla skalowania amplitud A_e ~ A_mu ~ A_tau, rdzen "widzi" g^{2*alpha} averaged
  po profilu. Im wieksze alpha, tym rdzen wnosi wiekszy efekt nieliniowy.

  Kluczowa hipoteza: m_obs jest wlasnoscia ASYMPTOTIC tail (energy radiated to
  infinity), ale ta energia jest renormalizowana przez "core dressing" o skali
  A^{c(alpha)}, gdzie c(alpha) zalezy od kineticnego prefactora.

Model: m_obs = (tail contribution to far-field charge) * (core wavefn norm.)
            ~ A^{2} * (g0^{2*alpha-2})^{f(g0)}

Pomijamy te subtelnosci na razie; sprawdzamy NUMERYCZNIE strukturalne relacje:

PLAN
----
1. Dla alpha-skanu, wyciagnij osobno A_tail, K_core (rdzen do r=5),
   K_tail (r=5..R_max), <g^{2*alpha}>_avg w rdzeniu.

2. Sprawdz czy m_obs/A^p = const tozsamy dla wszystkich e/mu/tau przy danym
   alpha — to nie jest "A^p" tylko "A^p * function(g0)".

3. Sprawdz czy istnieje wzor m_obs = A^a * g0^b z dwoma niezaleznymi
   wykladnikami, ktory dziala dla wszystkich alpha.

4. Test "tail action" S_tail = int_{r_tail}^{R_max} L * r^{d-1} dr i
   skalowanie tej wielkosci z A i g0.

Autor: Audyt 2026-05-01 (post-rezolucja, follow-up theoretical)
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
import math

PHI = (1 + math.sqrt(5)) / 2
G0_E = 0.86941
G0_MU = G0_E * PHI
R21_PDG = 206.7682
R31_PDG = 3477.23


def solve_ode(g0, alpha, d=3, r_max=300.0, n_points=30000, g_floor=1e-10):
    singular = [False]
    def rhs(r, y):
        g, gp = y
        if g < g_floor:
            singular[0] = True
            g = g_floor
        rhs_val = (1.0 - g) * g**(2 - 2*alpha)
        if r < 1e-12:
            gpp = rhs_val / max(d, 1.0)
        else:
            gpp = rhs_val - (alpha/g) * gp**2 - ((d-1.0)/r) * gp
        return [gp, gpp]
    r_eval = np.linspace(1e-10, r_max, n_points)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='RK45', t_eval=r_eval,
                    rtol=1e-10, atol=1e-12, max_step=0.05)
    return sol, singular[0]


def extract_atail(r, g, r_min=80.0, r_max=250.0):
    mask = (r >= r_min) & (r <= r_max)
    r_f = r[mask]
    if len(r_f) < 10:
        return None
    u_f = (g[mask] - 1.0) * r_f
    def model(r, B, C):
        return B * np.cos(r) + C * np.sin(r)
    try:
        popt, _ = curve_fit(model, r_f, u_f, p0=[0.01, 0.01])
        return math.sqrt(popt[0]**2 + popt[1]**2)
    except:
        return None


def compute_K_decomp(r, g, alpha, d=3, r_core=5.0, r_max_int=200.0):
    """Decompose K = K_core + K_tail."""
    gp = np.gradient(g, r)
    integrand = 0.5 * g**(2*alpha) * gp**2 * r**(d-1)

    mask_core = r <= r_core
    mask_tail = (r > r_core) & (r <= r_max_int)

    K_core = np.trapezoid(integrand[mask_core], r[mask_core])
    K_tail = np.trapezoid(integrand[mask_tail], r[mask_tail])

    # average <g^{2*alpha}> in core
    g_core = g[mask_core]
    g_2alpha_avg_core = np.mean(g_core**(2*alpha))

    return K_core, K_tail, g_2alpha_avg_core


# ================================================================
print("=" * 78)
print("  R3 — Analityczne dochodzenie p(alpha) = 5 - alpha")
print("=" * 78)
print()

# ----------------------------------------------------------------
# SECTION 1: Decomposition K = K_core + K_tail dla alpha-skanu
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 1: Decomposition K = K_core + K_tail dla e/mu (alpha-skan)")
print("=" * 78)

ALPHAS = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]

print(f"\n  {'alpha':>5} | {'A_e':>9} {'A_mu':>9} | "
      f"{'K_core_e':>9} {'K_core_mu':>10} | "
      f"{'K_tail_e':>9} {'K_tail_mu':>10} | "
      f"{'<g^2a>_e':>9} {'<g^2a>_mu':>10}")
print("  " + "-" * 100)

data = {}
for alpha in ALPHAS:
    sol_e, sing_e = solve_ode(G0_E, alpha)
    sol_mu, sing_mu = solve_ode(G0_MU, alpha)
    if sing_e or sing_mu or not sol_e.success or not sol_mu.success:
        continue

    A_e = extract_atail(sol_e.t, sol_e.y[0])
    A_mu = extract_atail(sol_mu.t, sol_mu.y[0])
    if A_e is None or A_mu is None:
        continue

    Kc_e, Kt_e, g2a_e = compute_K_decomp(sol_e.t, sol_e.y[0], alpha)
    Kc_mu, Kt_mu, g2a_mu = compute_K_decomp(sol_mu.t, sol_mu.y[0], alpha)

    data[alpha] = dict(A_e=A_e, A_mu=A_mu, Kc_e=Kc_e, Kc_mu=Kc_mu,
                       Kt_e=Kt_e, Kt_mu=Kt_mu, g2a_e=g2a_e, g2a_mu=g2a_mu)

    print(f"  {alpha:5.2f} | {A_e:9.5f} {A_mu:9.5f} | "
          f"{Kc_e:9.5f} {Kc_mu:10.5f} | "
          f"{Kt_e:9.5f} {Kt_mu:10.5f} | "
          f"{g2a_e:9.4f} {g2a_mu:10.4f}")

# ----------------------------------------------------------------
# SECTION 2: Test hypothesis m_obs = A^a * g0^b
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 2: Test hypothesis m_obs = c * A^a * g0^b (TWO exponents)")
print("=" * 78)
print()
print("  PDG: m_mu/m_e = 206.77")
print("  Equation: (A_mu/A_e)^a * (g0_mu/g0_e)^b = 206.77")
print("  g0_mu/g0_e = phi = 1.6180")
print()
print(f"  {'alpha':>5} | {'A_mu/A_e':>9} | {'log207':>7} | {'log_ratio':>10} | "
      f"{'a (z log)':>10} | {'p=5-α':>7}")
print("  " + "-" * 85)

# Single-exponent fit: assume b=0 (just A)
for alpha, d in data.items():
    ratio_A = d['A_mu'] / d['A_e']
    a_only = math.log(R21_PDG) / math.log(ratio_A)
    p_pred = 5.0 - alpha
    diff = (a_only - p_pred) / p_pred * 100
    print(f"  {alpha:5.2f} | {ratio_A:9.5f} | {math.log(R21_PDG):7.4f} | "
          f"{math.log(ratio_A):10.5f} | {a_only:10.5f} | {p_pred:7.4f}  "
          f"(diff {diff:+.2f}%)")


# ----------------------------------------------------------------
# SECTION 3: Test "renormalization" hypothesis
#   m_obs = A^p_tail * (core_dressing)^q
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 3: Decomposition m_obs ~ A^2 * (core_dressing)^?")
print("=" * 78)
print()
print("  Hipoteza: K_total = K_tail + K_core")
print("    K_tail ~ A^2 (universal, z linearyzacji tail)")
print("    K_core ~ A^? (zalezy od alpha)")
print("  Jesli m_obs ~ K_total^q dla pewnego q, to p_eff zalezy od mieszania.")
print()
print(f"  {'alpha':>5} | {'K_t/A^2 (e)':>11} | {'K_t/A^2 (mu)':>13} | "
      f"{'K_c/A^2 (e)':>11} | {'K_c/A^2 (mu)':>13}")
print("  " + "-" * 75)
for alpha, d in data.items():
    Kt_e_norm = d['Kt_e'] / d['A_e']**2
    Kt_mu_norm = d['Kt_mu'] / d['A_mu']**2
    Kc_e_norm = d['Kc_e'] / d['A_e']**2
    Kc_mu_norm = d['Kc_mu'] / d['A_mu']**2
    print(f"  {alpha:5.2f} | {Kt_e_norm:11.4f} | {Kt_mu_norm:13.4f} | "
          f"{Kc_e_norm:11.4f} | {Kc_mu_norm:13.4f}")

print()
print("  Obserwacja: K_t/A^2 ~= const (universal tail), ale K_c/A^2 zalezy od g0.")


# ----------------------------------------------------------------
# SECTION 4: Hipoteza FAR-FIELD CHARGE z field renormalization
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 4: Hipoteza A^p z core wavefn-renormalization Z(g0, alpha)")
print("=" * 78)

print("""
  Idea: dla pola Phi w substracie sprzezonego z grawitacja przez metryke
  g_eff[Phi]=Phi^4 (w sek08c), masa ZEWNETRZNA jaka widzi obserwator dla solitonu
  z A_tail amplituda jest:

    m_obs ~ A_tail^2 * Z(g0, alpha)

  gdzie Z to "wave-function renormalization factor" rdzenia. Wtedy:

    m_mu/m_e = (A_mu/A_e)^2 * Z(g0_mu)/Z(g0_e)

  Dla Z(g0) = g0^{n(alpha)} z liniowym n(alpha):
    log(m_mu/m_e) = 2*log(A_mu/A_e) + n(alpha) * log(g0_mu/g0_e)

  Solve dla n(alpha):
""")

print(f"  {'alpha':>5} | {'log r_A':>8} | {'log phi':>8} | "
      f"{'2 log r_A':>10} | {'n(alpha)':>9} | {'p_eff = 2 + n*log(phi)/log(r_A)':>}")
print("  " + "-" * 90)

LOG_PHI = math.log(PHI)
for alpha, d in data.items():
    log_rA = math.log(d['A_mu'] / d['A_e'])
    log_207 = math.log(R21_PDG)
    # log(207) = 2*log_rA + n * log(phi)
    n_alpha = (log_207 - 2 * log_rA) / LOG_PHI
    print(f"  {alpha:5.2f} | {log_rA:8.4f} | {LOG_PHI:8.4f} | "
          f"{2*log_rA:10.4f} | {n_alpha:9.4f}")


# ----------------------------------------------------------------
# SECTION 5: Powiazanie n(alpha) z 5-alpha
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 5: Czy n(alpha) ma postac liniowa? Empiryczny fit")
print("=" * 78)
print()

alphas_arr = np.array(list(data.keys()))
n_vals = []
for alpha in alphas_arr:
    d = data[alpha]
    log_rA = math.log(d['A_mu'] / d['A_e'])
    n_a = (math.log(R21_PDG) - 2 * log_rA) / LOG_PHI
    n_vals.append(n_a)
n_vals = np.array(n_vals)

# Liniowy fit n(alpha) = a*alpha + b
coeffs = np.polyfit(alphas_arr, n_vals, 1)
print(f"  Linear fit n(alpha) = a*alpha + b:")
print(f"    a = {coeffs[0]:.5f}")
print(f"    b = {coeffs[1]:.5f}")
print()
print(f"  {'alpha':>5} | {'n_numerical':>12} | {'n_fit':>9} | {'diff':>8}")
print("  " + "-" * 55)
for alpha in alphas_arr:
    d = data[alpha]
    log_rA = math.log(d['A_mu'] / d['A_e'])
    n_num = (math.log(R21_PDG) - 2 * log_rA) / LOG_PHI
    n_fit = coeffs[0] * alpha + coeffs[1]
    print(f"  {alpha:5.2f} | {n_num:12.5f} | {n_fit:9.5f} | {n_num-n_fit:+8.5f}")


# ----------------------------------------------------------------
# SECTION 6: Sprzezenie z empirycznym p(alpha) = 5-alpha
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 6: Sprzezenie z odkrytym p(alpha) = 5-alpha")
print("=" * 78)

print("""
  Empirycznie: m_obs / m_e = (A_mu/A_e)^p z p = 5-alpha

  Z dwuwykladnikowej hipotezy:
    m_obs / m_e = (A_mu/A_e)^2 * (g0_mu/g0_e)^n(alpha)
                = (A_mu/A_e)^2 * phi^n(alpha)

  Rownanie p-fitu:
    (A_mu/A_e)^(5-alpha) = (A_mu/A_e)^2 * phi^n(alpha)
    log(A_mu/A_e) * (5-alpha) = 2 log(A_mu/A_e) + n(alpha) * log(phi)
    (3-alpha) * log(A_mu/A_e) = n(alpha) * log(phi)
    n(alpha) = (3-alpha) * log(A_mu/A_e) / log(phi)

  Dla alpha=1: n(1) = 2 * log(3.79) / log(1.618) = 2 * 1.333 / 0.481 = 5.54
  Dla alpha=2: n(2) = 1 * log(5.91) / log(1.618) = 1.778 / 0.481 = 3.69

  TAK — n(alpha) NIE jest stale! Zaleznosc od alpha pochodzi z TEGO ze
  A_mu/A_e zalezy od alpha.

  Liniowy fit n(alpha) wyzej daje informacje czy powiazanie jest ladne.
""")

print()
print(f"  Bezposredni test p(alpha) = 5-alpha:")
print(f"  {'alpha':>5} | {'p_emp':>8} | {'5-alpha':>8} | {'diff%':>7}")
for alpha, d in data.items():
    ratio_A = d['A_mu'] / d['A_e']
    p_emp = math.log(R21_PDG) / math.log(ratio_A)
    p_pred = 5.0 - alpha
    diff = (p_emp - p_pred) / p_pred * 100
    print(f"  {alpha:5.2f} | {p_emp:8.4f} | {p_pred:8.4f} | {diff:+7.2f}%")


# ----------------------------------------------------------------
# SECTION 7: WARNING — empirical fit might be specific to phi-ladder
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 7: KLUCZOWY OPEN PROBLEM")
print("=" * 78)

print("""
  Empiryczna p(alpha) = 5-alpha jest PRZYBLIZONA, dokladna dla alpha=1, 2.
  Nie zostala derywowana analitycznie z fizyki solitonu — pozostaje OPEN PROBLEM.

  Dwa silnik dowody numeryczne ze relacja jest realna (nie artefakt):
  1. p = 4.001 dla alpha=1 (R3 oryginalne)
  2. p = 3.001 dla alpha=2 (TGP-canonical, sympy diff zero)

  Mozliwe sciezki derywacji:
  (A) Decomposition K_total = K_tail (~ A^2 universal) + K_core (~ A^? z alpha)
  (B) Wave-function renormalization Z(g0, alpha): m_obs = A^2 * Z
  (C) Field theory PDE analysis: solving Phi-EOM far-field with source-multipole

  Bez derywacji p(alpha)=5-alpha pozostaje EMPIRYCZNYM ODKRYCIEM, nie
  FIRST-PRINCIPLES wzoru. To jest OK — wiele formul w fizyce zostalo odkrytych
  empirycznie zanim derywowano (Rydberg, Balmer, Bohr, Planck...). Status:

    [EMPIRICAL DISCOVERY] — sprawdzone numerycznie, czeka na derywacje analityczna

  Sciezka A wydaje sie najobiecujaca — zaleznosc od alpha wchodzi przez core
  contribution K_core ~ A^2 * <g^{2*alpha}>, gdzie srednia po profilu rdzenia
  zalezy nieliniowo od g0 (i przez to od alpha).

  AKCJA: ten skrypt udokumentowuje JAWNIE ze p(alpha)=5-alpha jest empiryczny,
  nie analityczny. Wymaga osobnego studium — proponowany cykl R3-Q5 (z R5 bridge).
""")
