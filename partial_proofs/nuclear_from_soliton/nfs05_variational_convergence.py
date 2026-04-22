"""
nfs05_variational_convergence.py — test czy residual overbinding to Pauli czy wariacyjne?

Kontekst: nfs04 z HC (V_a=86.7, V_r=2400) dawało tryton E = -11.27 MeV (obs -8.48),
residual overbinding 1.33×. Diagnozowane jako "Pauli-gap" (residual).

ALE: może residual to wariacyjny błąd single-Gauss ansatz? Gauss-Jastrow zbyt
sztywne — Jastrow exp(-β/2·Σr_ij²) clumpi wszystkie pary jednocześnie. Prawdziwy
triton może mieć bardziej rozproszoną strukturę.

Test: multi-Gauss basis ψ = Σ_k c_k · exp(-β_k/2·Σr_ij²), solve generalized
eigenvalue H·c = E·S·c.

KLUCZOWA ZMIANA vs poprzedniej wersji: ANALYTYCZNE Yukawa-Gauss całki (bez MC).
To eliminuje szum numeryczny który w poprzedniej próbie dawał niestabilne GEV
(N=5 dało E=-285 MeV — artefakt MC noise).

Analytyczne wyrażenia (3-body CM, Gaussian basis, Yukawa pair V):
  Σ_{i<j} r_ij² = 3·(|x|² + |y|²)  [x,y = Jacobi]
  r_12 = √2 |x| (i podobnie dla innych par przez symetrię)

  S_ij = (2π/(3·α))^3  where α = (β_i+β_j)/2   [overlap, up to const]
  T_ij/S_ij = 9·β_i·β_j/α · ℏ²/(2m)            [kinetic / overlap ratio]
  V_ij/S_ij = 3 × <V_pair>_{3α}   [3 pary × integral Gaussa-Yukawy]

  <V_pair>_{A} dla V_Y(r) = ±V_0·exp(-r/a)/(r/a) = ±V_0·a·exp(-r/a)/r:

    <V_Y>_{A} = (integral V_Y(√2r) · exp(-3A·r²)·4πr² dr) /
                (integral exp(-3A·r²)·4πr² dr)
              = (2·V_0·a/(√2)) · I_1(3A, √2/a) / (π/(3A))^(3/2)

  I_1(A, γ) = ∫_0^∞ r·exp(-A·r² - γ·r) dr
           = (1/(2A))·[1 - γ·√(π/A)/2 · exp(γ²/(4A)) · erfc(γ/(2√A))]

Testy:
  T1: single-Gauss baseline reproduces nfs04 (E ≈ -11.27 MeV)
  T2: multi-Gauss (N=3) daje E_min ≤ single-Gauss (variational principle)
  T3: multi-Gauss (N=5) mała zmiana vs N=3 (near convergence)
  T4: final E_multi vs E_obs: residual fizyczny? (Pauli pozostaje)
  T5: r_rms multi-Gauss w fizycznym zakresie
"""

import math
import sys
import io
import numpy as np
from scipy.linalg import eigh
from scipy.special import erfc

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
np.random.seed(42)

print("=" * 78)
print("  nfs05 — Multi-Gauss basis: test czy residual to Pauli czy wariacyjne?")
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

V_a = 86.734
V_r = 2400.0

print(f"\n  V_NN parametry (z nfs04, dual-gap calibrated):")
print(f"    V_a (attrakcyjne, pion-like)  = {V_a} MeV,  a_a = {a_a:.3f} fm")
print(f"    V_r (repulsyjne, rho-like)    = {V_r} MeV,  a_r = {a_r:.3f} fm")

# ---------------------------------------------------------------------------
# Analytical Gauss-Yukawa integral
# ---------------------------------------------------------------------------
def I1(A, gamma):
    """∫_0^∞ r·exp(-A·r² - γ·r) dr = (1/(2A))·[1 - g·√(π/A)/2·exp(g²/(4A))·erfc(g/(2√A))]"""
    s = gamma / (2.0 * math.sqrt(A))
    # Numerical stability: exp(s²)·erfc(s) = erfcx(s) for large s
    # For small s, direct; for large, use erfcx
    if s < 6.0:
        term = 0.5 * gamma * math.sqrt(math.pi / A) * math.exp(s*s) * erfc(s)
    else:
        # erfcx(s) ≈ 1/(s√π) for large s
        from scipy.special import erfcx
        term = 0.5 * gamma * math.sqrt(math.pi / A) * erfcx(s)
    return (1.0 / (2.0 * A)) * (1.0 - term)

def V_matrix_element_pair(alpha, V_0, a_yuk):
    """
    <V_Y(r_12)>_{ψ(α)} gdzie ψ(α) ∝ exp(-α·Σr_ij²) w 3-body CM.

    Zwraca <V_Y>_{density ∝ exp(-α·Σr²)} = numerator / S_pair_norm

    V_Y(r) = ±V_0·exp(-r/a)/(r/a) = ±V_0·a·exp(-r/a)/r
    r_12 = √2 |x|, more convenient to integrate in x-space.

    Full 6D integral:
      ∫ exp(-α·Σr²) V_Y(r_12) d^6ξ
      = ∫ dy exp(-3α|y|²) · ∫ dx exp(-3α|x|²) · V_0·a·exp(-√2|x|/a)/(√2|x|)·sign
      = (π/(3α))^(3/2) · (V_0·a/√2) · 4π · ∫_0^∞ r·exp(-3α·r² - √2/a·r) dr
      = (π/(3α))^(3/2) · (4π·V_0·a/√2) · I_1(3α, √2/a)

    Dividing by overlap S = (π/(3α))^3:
      <V_Y>_{per pair, density norm} = (4π·V_0·a/√2) / (π/(3α))^(3/2)) · I_1(3α, √2/a)
                                     = (4π·V_0·a/√2) · (3α/π)^(3/2) · I_1(3α, √2/a)
    """
    A_eff = 3.0 * alpha
    gamma = math.sqrt(2.0) / a_yuk
    prefactor = (4.0 * math.pi * V_0 * a_yuk / math.sqrt(2.0)) * (3.0*alpha/math.pi)**1.5
    return prefactor * I1(A_eff, gamma)

def V_total_triton_analytic(alpha):
    """
    <V_NN_total>_{density ∝ exp(-α·Σr²)} for triton (sum over 3 pairs).
    V_NN = -V_a·Yukawa(a_a) + V_r·Yukawa(a_r)
    """
    V_att_per_pair = V_matrix_element_pair(alpha, -V_a, a_a)
    V_rep_per_pair = V_matrix_element_pair(alpha, +V_r, a_r)
    return 3.0 * (V_att_per_pair + V_rep_per_pair)

def T_over_S_triton(alpha):
    """T_ii/S_ii = 9·α · ℏ²/(2m) for single Gauss with β=α (density α = β single-Gauss)."""
    return 9.0 * alpha * hbarc**2 / (2.0 * m_N)

def energy_single_gauss(beta):
    """Total E(β) = T(β) + V_total(β) gdzie β = density param, ψ(β)² density = exp(-β·Σr²)."""
    T = T_over_S_triton(beta)
    V = V_total_triton_analytic(beta)
    return T + V, T, V

# ---------------------------------------------------------------------------
# Sanity check: analytical vs MC-based nfs04 wynik
# ---------------------------------------------------------------------------
print("\n[T1] Single-Gauss baseline analitycznie (scan β):")
betas_scan = np.logspace(-1.5, 0.5, 60)
E_scan = []
for b in betas_scan:
    E, T, V = energy_single_gauss(b)
    E_scan.append((b, T, V, E))

E_scan = sorted(E_scan, key=lambda x: x[3])
b_single_opt, T_s, V_s, E_single = E_scan[0]
# r_rms: Σr_ij² density <Σr_ij²>_β = 3/β (3 pairs, each <r_ij²> = 1/β)
# So average pair <r_ij²> = 1/β, r_ij_rms = 1/√β
r_pair_rms = 1.0 / math.sqrt(b_single_opt)
# r_rms for nucleon = sqrt(<|r_i|²>) with CM frame: for 3 equal masses, Σ|r_i|² = (Σr_ij²)/3
# So <|r_i|²> = <Σr_ij²>/9 = (3/β)/9 = 1/(3β). r_rms = 1/√(3β)
r_single = 1.0 / math.sqrt(3.0 * b_single_opt)
print(f"  Optimalny single-Gauss β = {b_single_opt:.4f}")
print(f"  E_single = {E_single:.3f} MeV  (T={T_s:.2f}, V={V_s:.2f})")
print(f"  r_rms (particle) = {r_single:.3f} fm,  r_pair_rms = {r_pair_rms:.3f} fm")
print(f"  (nfs04 MC dawało E ≈ -11.27 MeV — sanity check)")
check(abs(E_single - (-11.27)) < 2.0,
      "T1: single-Gauss analytic zgodny z nfs04 MC (~ -11.27 MeV)",
      f"E_single = {E_single:.2f} MeV")

# ---------------------------------------------------------------------------
# Multi-Gauss matrix elements (analytical)
# ---------------------------------------------------------------------------
def S_element(b_i, b_j):
    """S_ij = (π/(3α))^3, α = (b_i+b_j)/2. Unit: powinno być skonsystentne."""
    alpha = (b_i + b_j) / 2.0
    return (math.pi / (3.0 * alpha))**3

def T_element(b_i, b_j):
    """T_ij = S_ij · 9·b_i·b_j/α · ℏ²/(2m), α = (b_i+b_j)/2"""
    alpha = (b_i + b_j) / 2.0
    return S_element(b_i, b_j) * 9.0 * b_i * b_j / alpha * hbarc**2 / (2.0 * m_N)

def V_element(b_i, b_j):
    """V_ij = S_ij · <V_NN_total>_α, α = (b_i+b_j)/2"""
    alpha = (b_i + b_j) / 2.0
    return S_element(b_i, b_j) * V_total_triton_analytic(alpha)

def build_matrices(betas):
    nB = len(betas)
    S = np.zeros((nB, nB))
    T = np.zeros((nB, nB))
    V = np.zeros((nB, nB))
    for i in range(nB):
        for j in range(nB):
            S[i,j] = S_element(betas[i], betas[j])
            T[i,j] = T_element(betas[i], betas[j])
            V[i,j] = V_element(betas[i], betas[j])
    H = T + V
    return H, S

def solve_gev(betas):
    H, S = build_matrices(betas)
    # Condition number of S
    eig_S = np.linalg.eigvalsh(S)
    cond_S = eig_S[-1] / max(eig_S[0], 1e-30)
    try:
        eigvals, eigvecs = eigh(H, S)
        return eigvals[0], eigvecs[:, 0], H, S, cond_S
    except np.linalg.LinAlgError as e:
        print(f"  EIGH failed: {e}")
        # Regularize: S + ε·I
        eps = 1e-8 * np.trace(S) / len(betas)
        eigvals, eigvecs = eigh(H, S + eps*np.eye(len(betas)))
        return eigvals[0], eigvecs[:, 0], H, S, cond_S

# ---------------------------------------------------------------------------
# T2: N=3 multi-Gauss
# ---------------------------------------------------------------------------
print("\n[T2] Multi-Gauss (N=3) test (analityczny):")
betas_3 = np.geomspace(0.3 * b_single_opt, 3.0 * b_single_opt, 3)
print(f"  Basis β: {betas_3}")
E_3, c_3, H_3, S_3, cS_3 = solve_gev(betas_3)
print(f"  E(N=3) = {E_3:.4f} MeV,  cond(S) = {cS_3:.2e}")
print(f"  Δ vs single-Gauss  = {E_3 - E_single:+.4f} MeV  (musi być ≤ 0)")
check(E_3 <= E_single + 0.01,
      "T2: multi-Gauss (N=3) ≤ single-Gauss (variational principle)",
      f"ΔE = {E_3 - E_single:+.4f} MeV")

# ---------------------------------------------------------------------------
# T3: N=5 multi-Gauss
# ---------------------------------------------------------------------------
print("\n[T3] Multi-Gauss (N=5) convergence test:")
betas_5 = np.geomspace(0.1 * b_single_opt, 10.0 * b_single_opt, 5)
print(f"  Basis β: {betas_5}")
E_5, c_5, H_5, S_5, cS_5 = solve_gev(betas_5)
print(f"  E(N=5) = {E_5:.4f} MeV,  cond(S) = {cS_5:.2e}")
print(f"  Δ vs N=3 = {E_5 - E_3:+.4f} MeV")
rel_change_35 = abs(E_5 - E_3) / max(abs(E_3), 0.1)
print(f"  Relative change = {rel_change_35*100:.1f}%")
check(rel_change_35 < 0.10,
      "T3: N=5 zbieżne vs N=3 (<10% zmiana)",
      f"{rel_change_35*100:.1f}%")

# ---------------------------------------------------------------------------
# T3b: N=7 dla extra convergence
# ---------------------------------------------------------------------------
print("\n[T3b] Multi-Gauss (N=7) extra convergence:")
betas_7 = np.geomspace(0.05 * b_single_opt, 20.0 * b_single_opt, 7)
E_7, c_7, _, _, cS_7 = solve_gev(betas_7)
print(f"  E(N=7) = {E_7:.4f} MeV,  cond(S) = {cS_7:.2e}")
print(f"  Δ vs N=5 = {E_7 - E_5:+.4f} MeV")

# N=9
print("\n[T3c] Multi-Gauss (N=9) super-convergence:")
betas_9 = np.geomspace(0.03 * b_single_opt, 30.0 * b_single_opt, 9)
E_9, c_9, _, _, cS_9 = solve_gev(betas_9)
print(f"  E(N=9) = {E_9:.4f} MeV,  cond(S) = {cS_9:.2e}")
print(f"  Δ vs N=7 = {E_9 - E_7:+.4f} MeV")

# ---------------------------------------------------------------------------
# T4: final residual vs observation
# ---------------------------------------------------------------------------
E_obs = -8.48
E_best = min(E_3, E_5, E_7, E_9)
residual = E_best - E_obs
residual_abs = abs(residual)
drop_total = E_single - E_best

print(f"\n[T4] Finalne porównanie:")
print(f"  E_single-Gauss  = {E_single:.3f} MeV")
print(f"  E_multi-Gauss   = {E_best:.3f} MeV  (best of N=3..9)")
print(f"  E_observed      = {E_obs} MeV")
print(f"  Residual (multi)   = {residual:+.3f} MeV")
print(f"  Residual (single)  = {E_single - E_obs:+.3f} MeV")
print(f"  Wariacyjny drop single → multi = {drop_total:.3f} MeV")

if residual_abs < 0.5:
    interp = "WARIACYJNY ARTEFAKT — single-Gauss był niezbieżny, Pauli-diagnosis słaba"
elif drop_total > 1.5 and residual_abs > 1.5:
    interp = "MIESZANY — wariacja miała wpływ, residual trzyma (Pauli fizyczny)"
elif drop_total < 0.5 and residual_abs > 1.0:
    interp = "FIZYCZNY RESIDUAL — Pauli przeżywa wariację"
elif drop_total < 1.0 and residual_abs > 1.5:
    interp = "FIZYCZNY RESIDUAL — wariacja mała, Pauli silny"
else:
    interp = "POŚREDNI"

# Success criterion: jeśli multi-Gauss NIE ZMIENIA dużo E single (wariacyjny poprawny)
# I residual pozostaje > 1 MeV, to Pauli-diagnosis trzyma
check(residual_abs > 1.0 and drop_total < 2.0,
      "T4: residual trzyma po konwergencji (Pauli-diagnosis PHYSICAL)",
      f"residual {residual:+.2f}, drop {drop_total:.2f} → {interp}")

# ---------------------------------------------------------------------------
# T5: r_rms po multi-Gauss
# ---------------------------------------------------------------------------
print(f"\n[T5] r_rms po multi-Gauss — dominant β analysis:")
# Dla eigvector c_9, weighted <β> = Σ_ij c_i·c_j·S_ij·β_avg_ij / (c^T·S·c)
# Use effective β = <Σr²>^(-1) through energy ratio
norm_9 = c_9.T @ S_9 if 'S_9' in dir() else 1.0  # reuse
# Simple: normalize c_9 and compute effective β via <1/r²> type estimator
# <Σr_ij²>_multi = Σ_ij c_i·c_j · (Σr_ij²)_{ψ_i ψ_j} / S_full
# For Gauss basis, <Σr²>_{ψ_i ψ_j} / S_ij = 3/α where α = (b_i+b_j)/2
# So <Σr²>_multi = Σ_ij c_i·c_j·S_ij·3/α_ij / (c^T·S·c)
# Then <r_pair²> = <Σr²>/3, <|r_i|²> = <Σr²>/9

H_9, S_9_fresh = build_matrices(betas_9)
norm_9 = float(c_9 @ S_9_fresh @ c_9)
c_n = c_9 / math.sqrt(norm_9)

sum_r2 = 0.0
for i in range(len(betas_9)):
    for j in range(len(betas_9)):
        alpha_ij = (betas_9[i] + betas_9[j]) / 2.0
        sum_r2 += c_n[i] * c_n[j] * S_9_fresh[i,j] * (3.0 / alpha_ij)
r_pair_multi = math.sqrt(sum_r2 / 3.0)
r_particle_multi = math.sqrt(sum_r2 / 9.0)
print(f"  <Σr_ij²>_multi = {sum_r2:.3f} fm²")
print(f"  <r_pair²>^(1/2) = {r_pair_multi:.3f} fm  (single: {r_pair_rms:.3f})")
print(f"  <|r_i|²>^(1/2)  = {r_particle_multi:.3f} fm  (single: {r_single:.3f})")
check(1.0 < r_particle_multi < 3.0,
      "T5: r_rms multi-Gauss w fizycznym zakresie",
      f"r_rms = {r_particle_multi:.3f} fm")

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print(f"  nfs05 — WYNIK: {PASS_COUNT}/{PASS_COUNT+FAIL_COUNT} PASS")
print("=" * 78)

print(f"""
  TEST KONWERGENCJI WARIACYJNEJ (V_NN z nfs04, ANALITYCZNE Yukawa-Gauss):

    Ansatz          |  E (MeV)  |  ΔE vs single  |  cond(S)
    ----------------+-----------+----------------+-----------
    Single-Gauss    |  {E_single:>6.3f}  |    (baseline)  |  --
    Multi-Gauss N=3 |  {E_3:>6.3f}  |  {E_3 - E_single:>+6.3f}       |  {cS_3:.2e}
    Multi-Gauss N=5 |  {E_5:>6.3f}  |  {E_5 - E_single:>+6.3f}       |  {cS_5:.2e}
    Multi-Gauss N=7 |  {E_7:>6.3f}  |  {E_7 - E_single:>+6.3f}       |  {cS_7:.2e}
    Multi-Gauss N=9 |  {E_9:>6.3f}  |  {E_9 - E_single:>+6.3f}       |  {cS_9:.2e}
    Observed triton |  {E_obs:>6.3f}  |       --       |  --

  INTERPRETACJA: {interp}

  WNIOSEK:""")

if drop_total < 0.5:
    print("    • Multi-Gauss prawie niezmienia E vs single → single już zbieżny")
elif drop_total < 1.5:
    print(f"    • Multi-Gauss E spada o {drop_total:.2f} MeV — wariacja umiarkowana")
else:
    print(f"    • Multi-Gauss E spada o {drop_total:.2f} MeV — wariacja duża!")

if residual_abs > 1.5:
    print(f"    • Residual {residual:.2f} MeV — pozostaje fizyczna luka (Pauli/inne)")
elif residual_abs < 0.5:
    print(f"    • Residual {residual:.2f} MeV — TGP V_NN + wariacja zbieżna do obs!")
else:
    print(f"    • Residual {residual:.2f} MeV — marginalny, niejasny wniosek")

print(f"""
  ZNACZENIE:
    Ten test dystynguje: 'TGP brakuje fermionów' (głębokie) vs 'Gauss-Jastrow
    zbyt sztywne' (łatwo fixable). Multi-Gauss daje elastyczniejszy ansatz,
    więc jeśli E zbiega do -8.48 MeV → wariacyjny artefakt.
    Jeśli E stabilne ~ -11 MeV → fizyczny gap (Pauli ?)
""")
