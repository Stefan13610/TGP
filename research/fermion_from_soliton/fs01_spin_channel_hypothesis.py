"""
fs01_spin_channel_hypothesis.py — test hipotezy "Pauli-gap = spin-isospin channel gap"

Kontekst:
  nuclear_from_soliton dał residual overbinding 1.45× dla triton (po nfs05
  multi-Gauss convergence), 1.22-1.36× dla alpha. Residual jest fizyczny
  (nie wariacyjny artefakt). Zdiagnozowany jako "Pauli-gap".

Hipoteza:
  W S-wave jądrowych podstaw, spatial wavefunction JEST symmetryczna. Pauli
  nie antysymetryzuje spatial directly. Zamiast tego Pauli ZMUSZA pary do
  konkretnych kanałów spin×izospin:
    - (T=0, S=1) deuteron-like triplet: silnie wiązany (V_triplet)
    - (T=1, S=0) np singlet / nn,pp singlet: słabiej wiązany (V_singlet)
    - T=1, S=1 / T=0, S=0 w L=0 ZABRONIONE przez Pauli.

  TGP V_NN z nfs04 jest SKALARNY (niezależny od spin/izospin), zfitowany do
  deuteronu = pure T=0,S=1 triplet. Stosując go do triton/alpha z "all triplet"
  założeniem przeszacowujemy binding.

Prawidłowe channel weighting (DERYWACJA przez fs02 explicit Slater det):

  Explicit SU(2)_spin × SU(2)_iso analiza (fs02) daje per-pair:
    <P_{T=1}>_avg = <P_{T=0}>_avg = 1/2
    <P_{S=1}>_avg = <P_{S=0}>_avg = 1/2
    W L=0 allowed tylko T0S1 i T1S0, po 1/2 każdy.

  Triton (S=1/2, T=1/2, 3 pary):
    <Σ P_{T1S0}> = 1.5,  <Σ P_{T0S1}> = 1.5 (derywacja fs02)
    <V>_triton_eff = 1.5 V_{T1S0} + 1.5 V_{T0S1} = 1.5 V_t (1 + f_s)
    gdzie f_s = V_{T1S0}/V_{T0S1}, V_t = V_{T0S1}
    Naive: 3 V_t → reduction factor w = (1 + f_s) / 2

  Alpha (S=0, T=0, 6 par):
    <Σ P_{T1S0}> = 3.0,  <Σ P_{T0S1}> = 3.0 (derywacja fs02)
    <V>_alpha_eff = 3 V_{T1S0} + 3 V_{T0S1} = 3 V_t (1 + f_s)
    Naive: 6 V_t → reduction factor (1 + f_s) / 2 — IDENTYCZNY!

  Zatem hipoteza PRZEWIDUJE:
    Pauli-channel reduction factor w = (1+f_s)/2 taki sam dla triton i alpha
    → spójność systemu (zgodne z explicit Slater det fs02).

  UWAGA: wersja v1 tego pliku używała formuły w=(1+2f_s)/3 (heurystyka
  "nn-forced + np 50/50") — okazała się błędna po fs02. Poprawna jest
  symetryczna 50-50 per-pair, która daje w=(1+f_s)/2.

Test:
  T1: Extract T, V separately dla triton single-Gauss (baseline)
  T2: Znajdź channel factor w taki że E_eff = T + w·V = E_obs = -8.48 MeV
  T3: Wyliczyć implied f_s = 2w-1, sprawdzić że 0 < f_s < 1 (fizyczny)
  T4: Sprawdzić że ten sam f_s daje alpha E_eff bliżej E_obs = -28.3 MeV
  T5: Compare f_s z OPE-like expectations (V_singlet ~ 0.5-0.9 V_triplet)
"""

import math
import sys
import io
import numpy as np
from scipy.special import erfc, erfcx

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
np.random.seed(42)

print("=" * 78)
print("  fs01 — Spin-isospin channel hypothesis for nuclear Pauli-gap")
print("=" * 78)

PASS_COUNT = 0
FAIL_COUNT = 0

def check(cond, label, info=""):
    global PASS_COUNT, FAIL_COUNT
    status = "PASS" if cond else "FAIL"
    if cond: PASS_COUNT += 1
    else:    FAIL_COUNT += 1
    print(f"  [{status}] {label}" + (f"  ({info})" if info else ""))

# Physical constants
hbarc = 197.327
m_N = 938.918
m_pi = 138.04
m_rho = 770.0
a_a = hbarc / m_pi
a_r = hbarc / m_rho

# V_NN parameters (z nfs04, fit do deuteronu = pure T=0,S=1 triplet channel)
V_a = 86.734
V_r = 2400.0

# Observed binding energies
E_triton_obs = -8.48   # MeV
E_alpha_obs = -28.3    # MeV

print(f"\n  V_NN triplet channel (z nfs04, deuteron fit):")
print(f"    V_a = {V_a} MeV,  a_a = {a_a:.3f} fm  (pion attraction)")
print(f"    V_r = {V_r} MeV,  a_r = {a_r:.3f} fm  (rho repulsion)")
print(f"\n  Observed: E_triton = {E_triton_obs} MeV,  E_alpha = {E_alpha_obs} MeV")

# ---------------------------------------------------------------------------
# Analytical Yukawa-Gauss integrals (z nfs05)
# ---------------------------------------------------------------------------
def I1(A, gamma):
    s = gamma / (2.0 * math.sqrt(A))
    if s < 6.0:
        term = 0.5 * gamma * math.sqrt(math.pi / A) * math.exp(s*s) * erfc(s)
    else:
        term = 0.5 * gamma * math.sqrt(math.pi / A) * erfcx(s)
    return (1.0 / (2.0 * A)) * (1.0 - term)

def V_pair_triton(alpha, V_0, a_yuk):
    """<V_Y(r_12)>_{density ∝ exp(-α·Σr_ij²)} 3-body triton."""
    A_eff = 3.0 * alpha
    gamma = math.sqrt(2.0) / a_yuk
    prefactor = (4.0 * math.pi * V_0 * a_yuk / math.sqrt(2.0)) * (3.0*alpha/math.pi)**1.5
    return prefactor * I1(A_eff, gamma)

def V_NN_triton(alpha):
    """Per-pair V_NN_triplet for 3-body Jastrow, density α."""
    return V_pair_triton(alpha, -V_a, a_a) + V_pair_triton(alpha, V_r, a_r)

def T_triton(beta):
    """Total kinetic for 3-body Jastrow ψ(β) symmetric."""
    return 9.0 * beta * hbarc**2 / (2.0 * m_N)

# ---------------------------------------------------------------------------
# Alpha: 4-body kinetic + potential
# Kinetic scaling: for N-body Gauss, T_coeff = 3·C(N,2)·β·ℏ²/(2m)
# For alpha N=4: T_coeff = 3·6 = 18β·ℏ²/(2m)
# ---------------------------------------------------------------------------
def T_alpha(beta):
    """Kinetic for 4-body Jastrow."""
    return 18.0 * beta * hbarc**2 / (2.0 * m_N)

def V_pair_alpha(alpha, V_0, a_yuk):
    """Per-pair V w Jastrow 4-body. Σr_ij² = 4·(|ξ₁|²+|ξ₂|²+|ξ₃|²) dla 4 body.
    Pair r_12 = Jacobi coord. Używam heurystycznego scaling z nfs03:
    dla 4-body, per-pair geometry inna niż 3-body. Zwracam poprzez uśrednienie."""
    # Dla 4 body CM-frame, 3 Jacobi coords. Σr_ij² = 4·Σ|ξ_k|² (6 pairs × relation)
    # Sampling ξ_k z σ² = 1/(8β) daje density |ψ|²∝exp(-β·Σr²)=exp(-4β·Σ|ξ|²)
    # Ekvivalent: per-pair integral ma A_eff = 4α (dla 4-body density α)
    A_eff = 4.0 * alpha
    gamma = math.sqrt(2.0) / a_yuk
    prefactor = (4.0 * math.pi * V_0 * a_yuk / math.sqrt(2.0)) * (4.0*alpha/math.pi)**1.5
    return prefactor * I1(A_eff, gamma)

def V_NN_alpha(alpha):
    """<V_NN>_per pair alpha 4-body."""
    return V_pair_alpha(alpha, -V_a, a_a) + V_pair_alpha(alpha, V_r, a_r)

# ---------------------------------------------------------------------------
# T1: single-Gauss triton baseline, extract T i V separately
# ---------------------------------------------------------------------------
print("\n[T1] Triton single-Gauss: extract T, V separately")
betas_scan = np.logspace(-1.5, 0.5, 60)
best = None
for b in betas_scan:
    T = T_triton(b)
    V = 3.0 * V_NN_triton(b)  # 3 pary
    E = T + V
    if best is None or E < best[3]:
        best = (b, T, V, E)

b_t, T_t, V_t, E_t = best
print(f"  Optymalny β = {b_t:.4f}")
print(f"  T (3-body)  = {T_t:+.3f} MeV")
print(f"  V (naive)   = {V_t:+.3f} MeV")
print(f"  E = T + V   = {E_t:.3f} MeV  (nfs05 analityczny: -11.17)")
check(abs(E_t - (-11.17)) < 1.0,
      "T1: triton single-Gauss zgodny z nfs05",
      f"E = {E_t:.2f} MeV")

# ---------------------------------------------------------------------------
# T2: znajdź channel factor w takie że E_eff = -8.48
# E_eff(β, w) = T(β) + w · V_naive(β) → minimize w/r β, find w matching E_obs
# ---------------------------------------------------------------------------
def E_triton_channel(w, n_scan=60):
    """Znajdź min E dla danego channel factor w."""
    E_best = 1e10
    for b in np.logspace(-1.5, 0.5, n_scan):
        T = T_triton(b)
        V = 3.0 * V_NN_triton(b)
        E = T + w * V
        if E < E_best:
            E_best = E
            b_best = b
    return E_best, b_best

print("\n[T2] Scan channel factor w do match E_triton_obs:")
ws = np.linspace(0.5, 1.0, 51)
E_vals = []
for w in ws:
    E, b = E_triton_channel(w)
    E_vals.append((w, E, b))

# Find w such that E ≈ E_triton_obs (bisection-like)
w_match = None
for i in range(len(ws)-1):
    E1 = E_vals[i][1]
    E2 = E_vals[i+1][1]
    if (E1 - E_triton_obs) * (E2 - E_triton_obs) <= 0:
        # crossing
        frac = (E_triton_obs - E1) / (E2 - E1)
        w_match = ws[i] + frac * (ws[i+1] - ws[i])
        break

if w_match is None:
    print("  [ERROR] Nie znaleziono matching w w zakresie [0.5, 1.0]")
    # scan wider
    ws2 = np.linspace(0.3, 1.2, 91)
    for w in ws2:
        E, b = E_triton_channel(w)
        if E <= E_triton_obs + 0.5 and E >= E_triton_obs - 0.5:
            w_match = w
            break

print(f"  Channel factor w match: {w_match:.4f}")
E_check, b_check = E_triton_channel(w_match)
print(f"  Weryfikacja: E_triton(w={w_match:.3f}) = {E_check:.3f} MeV  (cel: {E_triton_obs})")

# Bazowe stwierdzenie: w <= 1 (Pauli zmniejsza binding, nie zwiększa)
check(w_match is not None and 0.5 < w_match < 1.0,
      "T2: channel factor w fizycznym zakresie (0.5, 1.0)",
      f"w = {w_match:.4f}")

# Robustność: też przetestuj z multi-Gauss baseline z nfs05
# Multi-Gauss converged triton E = -12.28 MeV (nfs05 N=9)
# Użyjemy T_single + V_multi = -12.28 jako baseline
print("\n[T2b] Robustność: multi-Gauss baseline (E=-12.28 MeV, nfs05):")
E_triton_multi_baseline = -12.28
# Single-Gauss T ~26 MeV (nfs05 dominant β), więc V_multi ≈ -12.28 - 26 = -38.28
# Proste przeskalowanie: jeśli V_multi jest więcej negative, ten sam w da
# inne E_effective. Scan:
V_multi_approx = E_triton_multi_baseline - T_t  # T_t z T1
w_ratio = (E_triton_obs - T_t) / V_multi_approx
print(f"  Approx V_multi ≈ {V_multi_approx:.2f} MeV (zakładając T ~ single)")
print(f"  Implied w = (E_obs - T)/V_multi = {w_ratio:.4f}")
f_s_multi = 2 * w_ratio - 1   # z fs02: w = (1+f_s)/2 → f_s = 2w - 1
print(f"  Implied f_s (multi) = {f_s_multi:.4f}  (wzór: f_s = 2w - 1, z fs02)")
check(0.3 < f_s_multi < 0.95,
      "T2b: f_s z multi-Gauss baseline też fizyczny",
      f"f_s_multi = {f_s_multi:.3f}")

# ---------------------------------------------------------------------------
# T3: implied f_s = V_singlet/V_triplet
# w = (1 + f_s)/2 → f_s = 2w - 1  (z fs02 explicit Slater derywacji)
# ---------------------------------------------------------------------------
if w_match is not None:
    f_s = 2 * w_match - 1
    print(f"\n[T3] Implied singlet/triplet ratio:")
    print(f"  f_s = V_{{T1S0}}/V_{{T0S1}} = 2w - 1 = {f_s:.4f}")
    print(f"  (Wzór z fs02 explicit Slater det: w = (1+f_s)/2)")
    print(f"  (Experimental estimate: f_s ~ 0.5-0.9 z analizy NN scattering)")
    check(0.0 < f_s < 1.0,
          "T3: f_s w fizycznym zakresie (0, 1)",
          f"f_s = {f_s:.3f}")
else:
    f_s = None
    check(False, "T3: nie wyznaczono f_s")

# ---------------------------------------------------------------------------
# T4: ten sam f_s dla alpha — sprawdź consistency
# Alpha channel factor = (1 + 2 f_s)/3 (IDENTYCZNE z triton, przewidywanie!)
# ---------------------------------------------------------------------------
if f_s is not None:
    print(f"\n[T4] Alpha prediction z tym samym f_s = {f_s:.3f}:")
    # w_alpha = (1 + 2 f_s)/3 = (2 + 4 f_s)/6
    # Ale można użyć równoważnej: w_alpha = w_triton (hipoteza)
    w_alpha = w_match

    # Scan alpha
    E_alpha_naive = 1e10
    E_alpha_eff = 1e10
    for b in np.logspace(-1.5, 0.5, 80):
        T_a = T_alpha(b)
        V_a_naive_val = 6.0 * V_NN_alpha(b)  # 6 par
        E_naive = T_a + V_a_naive_val
        E_eff = T_a + w_alpha * V_a_naive_val
        if E_naive < E_alpha_naive:
            E_alpha_naive = E_naive
            b_naive = b
        if E_eff < E_alpha_eff:
            E_alpha_eff = E_eff
            b_eff = b

    print(f"  Naive alpha (single-channel TGP, HC): E = {E_alpha_naive:.3f} MeV @ β={b_naive:.3f}")
    print(f"  Z channel weighting f_s={f_s:.3f}:  E = {E_alpha_eff:.3f} MeV @ β={b_eff:.3f}")
    print(f"  Observed alpha: {E_alpha_obs} MeV")
    print(f"  Residual po weighting: {E_alpha_eff - E_alpha_obs:+.3f} MeV")

    # Residual should be small if hypothesis works
    alpha_residual = abs(E_alpha_eff - E_alpha_obs)
    check(alpha_residual < 5.0,
          "T4: alpha E po channel weighting bliskie obs",
          f"residual = {alpha_residual:.2f} MeV vs {abs(E_alpha_obs - E_alpha_naive):.2f} MeV bez weighting")
else:
    check(False, "T4: skip (f_s not determined)")
    E_alpha_eff = None

# ---------------------------------------------------------------------------
# T5: porównanie z OPE expectations
# OPE central: V_{T1S0} / V_{T0S1} = 1 (same sign, same magnitude from σ·σ · τ·τ)
# Rzeczywistość: V_{T1S0} słabszy z powodu tensor + channel-dep HC
# Oczekiwana f_s: 0.5-0.9 z fits to NN scatterning
# ---------------------------------------------------------------------------
if f_s is not None:
    print(f"\n[T5] Porównanie z fizyką NN:")
    print(f"  OPE-only prediction: f_s = 1 (central V same, ignoruje tensor)")
    print(f"  Realistic NN fits: f_s ~ 0.6-0.9 (z tensor/HC channel dep)")
    print(f"  Nasz wynik: f_s = {f_s:.3f}")

    # Assessment
    if f_s < 0.3:
        interp = "EXTREME: V_singlet dużo słabsze — nietypowe"
    elif 0.3 <= f_s < 0.6:
        interp = "UMIARKOWANE: V_singlet ~0.5 V_triplet, zgodne z fitami NN"
    elif 0.6 <= f_s < 0.95:
        interp = "REALISTYCZNE: V_singlet ~0.7-0.9 V_triplet, typowe dla CD-Bonn/AV18"
    else:
        interp = "TRIVIAL: f_s≈1, channel weighting nie zmienia praktycznie nic"

    print(f"  Interpretacja: {interp}")
    check(0.3 < f_s < 0.95,
          "T5: f_s w zakresie eksperymentalnie rozsądnym",
          f"f_s = {f_s:.3f} → {interp}")
else:
    check(False, "T5: skip")

# ---------------------------------------------------------------------------
# Werdyk
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print(f"  fs01 — WYNIK: {PASS_COUNT}/{PASS_COUNT+FAIL_COUNT} PASS")
print("=" * 78)

if w_match is not None and f_s is not None:
    print(f"""
  HIPOTEZA CHANNEL-WEIGHTING:

  Dane wejściowe:
    E_triton_TGP_naive = {E_t:.2f} MeV  (single-Gauss, 3-pairs V_t)
    E_triton_obs       = {E_triton_obs} MeV

  Wynik:
    Channel factor w = {w_match:.4f}  (V_NN_effective = w · V_triplet)
    Implied f_s = V_s/V_t = {f_s:.4f}
    Consistent alpha prediction: {E_alpha_eff if E_alpha_eff else 'N/A':.2f} MeV (obs {E_alpha_obs})

  IMPLIKACJA:
    Jeśli f_s matching eksperymentalne wartości → Pauli-gap JEST "spin-channel
    gap" fenomenologicznie. TGP potrzebuje V_NN(T, S) zamiast scalar V_NN.

    Most strukturalny: SU(2) × SU(2) (spin × izospin) sektor TGP powinien
    generować channel dependence naturally.
""")
