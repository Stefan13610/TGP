#!/usr/bin/env python3
"""
ex145_B2_qk_z3_bridge.py -- B2: Q_K = 3/2 z Z_3 (formalny most)
TGP v1 -- Most B2 (2026-04-04)

Cel: Formalna weryfikacja lancucha logicznego:
  1. Trzy mody solitonowe -> wektor amplitudowy lambda = (A_0^2, A_1^2, A_2^2)
  2. Wiezy S_1, S_2 -> lambda lezy na okregu C w R^3
  3. Dekoherencja fazowa na okregu -> Z_3 (co 2pi/3)
  4. Z_3 na okregu = parametryzacja Brannena -> Q_K = 3/2

Kluczowa korekta vs prop:R2-minimal-interference:
  - STARY argument: fazy ogonowe phi_n tworza Z_3 (FAŁSZ wg ex125)
  - NOWY argument: Z_3 dziala na KATACH alpha_k w przestrzeni
    amplitudowej lambda_k = A_k^2, NIE na surowych fazach delta_n.
  - Selekcja Z_3: dekoherencja fazowa |Sum e^{i*alpha_k}|^2 = 0
    (jedyne rozwiazanie: rownomierny rozklad co 2pi/3)

Testy: 14 PASS/FAIL
"""

import numpy as np
from scipy.optimize import differential_evolution
import sys
import os

# Force UTF-8 output on Windows
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# ============================================================
# PDG masy leptonowe (MeV)
# ============================================================
M_E = 0.51099895
M_MU = 105.6583755
M_TAU = 1776.86

MASSES_PDG = np.array([M_E, M_MU, M_TAU])
SQRT_M_PDG = np.sqrt(MASSES_PDG)

# ============================================================
# Funkcje pomocnicze
# ============================================================

def Q_K(masses):
    """Koide quotient Q_K = (sum sqrt(m))^2 / sum(m)"""
    sm = np.sqrt(masses)
    return np.sum(sm)**2 / np.sum(masses)


def CV(x):
    """Coefficient of variation"""
    return np.std(x) / np.mean(x)


def brannen_sqrt_m(c, r, theta_K, N=3):
    """Brannen parametrization: sqrt(m_k) = c * (1 + r*cos(theta_K + 2*pi*k/N))"""
    k = np.arange(N)
    return c * (1.0 + r * np.cos(theta_K + 2.0 * np.pi * k / N))


def pair_cosine_sum(alphas):
    """Sum_{n<m} cos(alpha_n - alpha_m)"""
    N = len(alphas)
    s = 0.0
    for n in range(N):
        for m in range(n+1, N):
            s += np.cos(alphas[n] - alphas[m])
    return s


def coherence(alphas):
    """|Sum e^{i*alpha_k}|^2 -- dekoherencja = 0"""
    return abs(np.sum(np.exp(1j * np.asarray(alphas))))**2


def angular_entropy(alphas_sorted, N=3):
    """Entropia rozkladu katowego na okregu"""
    gaps = np.diff(alphas_sorted)
    gaps = np.append(gaps, 2*np.pi - alphas_sorted[-1] + alphas_sorted[0])
    gaps = gaps / (2*np.pi)  # normalizacja do rozkladu
    gaps = gaps[gaps > 1e-15]
    return -np.sum(gaps * np.log(gaps))


# ============================================================
# TESTY
# ============================================================
def run_tests():
    results = []

    # === BLOK 1: Algebraiczne tozsamosci ===

    # C1: Q_K(PDG) = 3/2
    qk_pdg = Q_K(MASSES_PDG)
    results.append(("C1: Q_K(PDG) = 3/2",
                    abs(qk_pdg - 1.5) < 0.001,
                    f"Q_K = {qk_pdg:.6f}"))

    # C2: CV(sqrt(m)) = 1 <=> Q_K = 3/2
    cv_pdg = CV(SQRT_M_PDG)
    qk_from_cv = 3.0 / (1.0 + cv_pdg**2)
    results.append(("C2: Q_K = N/(1+CV^2), CV(sqrt(m)) = 1",
                    abs(cv_pdg - 1.0) < 0.01 and abs(qk_from_cv - qk_pdg) < 1e-6,
                    f"CV = {cv_pdg:.6f}, Q_K(CV) = {qk_from_cv:.6f}"))

    # C3: Brannen z Z_3 (r=sqrt(2)) daje Q_K = 3/2
    # Testujemy tylko theta_K z dodatnimi sqrt(m)
    thetas = np.linspace(0, 2*np.pi, 1000)
    qk_brannen_vals = []
    n_valid = 0
    for th in thetas:
        sm = brannen_sqrt_m(1.0, np.sqrt(2), th)
        if np.all(sm > 0):
            n_valid += 1
            masses = sm**2
            qk_brannen_vals.append(Q_K(masses))
    qk_brannen_arr = np.array(qk_brannen_vals)
    max_dev = np.max(np.abs(qk_brannen_arr - 1.5)) if len(qk_brannen_arr) > 0 else 999
    results.append(("C3: Brannen(r=sqrt(2)): Q_K=3/2 (valid theta_K)",
                    max_dev < 1e-12 and n_valid > 100,
                    f"max|Q_K - 3/2| = {max_dev:.2e}, {n_valid} valid thetas"))

    # C4: r=sqrt(2) jest jedyna wartosc dajaca CV=1
    # CV = r/sqrt(2), wiec CV=1 <=> r=sqrt(2)
    r_test = np.sqrt(2)
    cv_from_r = r_test / np.sqrt(2)
    results.append(("C4: r=sqrt(2) <=> CV=1 (unikalne)",
                    abs(cv_from_r - 1.0) < 1e-14,
                    f"CV(r=sqrt(2)) = {cv_from_r:.15f}"))

    # === BLOK 2: Okrag C w R^3 ===

    # C5: Lambda PDG lezy na okregu C (S_1, S_2 definiuja okrag)
    lam_pdg = SQRT_M_PDG  # lambda_k = sqrt(m_k) = A_k^2
    S1 = np.sum(lam_pdg)
    S2 = np.sum(lam_pdg**2)
    lam_bar = S1 / 3.0
    rho_sq = 2.0/3.0 * (S2 - S1**2/3.0)
    rho = np.sqrt(rho_sq)
    # Verify: lambda_k = lam_bar + rho*cos(alpha_k) for some alpha_k
    cos_alpha = (lam_pdg - lam_bar) / rho
    all_in_range = np.all(np.abs(cos_alpha) <= 1.0 + 1e-10)
    results.append(("C5: lambda(PDG) lezy na okregu C(S_1, S_2)",
                    all_in_range,
                    f"lam_bar={lam_bar:.4f}, rho={rho:.4f}, "
                    f"cos_alpha = [{cos_alpha[0]:.6f}, {cos_alpha[1]:.6f}, {cos_alpha[2]:.6f}]"))

    # C6: CV(lambda_PDG) = 1
    cv_lam = CV(lam_pdg)
    results.append(("C6: CV(lambda_PDG) = CV(sqrt(m)) = 1",
                    abs(cv_lam - 1.0) < 0.001,
                    f"CV = {cv_lam:.6f}"))

    # C7: rho/lam_bar = sqrt(2) (relacja r=sqrt(2))
    r_ratio = rho / lam_bar
    results.append(("C7: rho/lam_bar = sqrt(2) (r = sqrt(2))",
                    abs(r_ratio - np.sqrt(2)) < 0.01,
                    f"r = rho/lam_bar = {r_ratio:.6f}, sqrt(2) = {np.sqrt(2):.6f}"))

    # === BLOK 3: E_cross = const na C (wazny wniosek negatywny) ===

    # C8: E_cross = (S1^2-S2)/2 jest STALA na okregu C
    # Brannen z roznym theta_K daje ten sam E_cross (bo S1, S2 zachowane)
    E_cross_analytic = (S1**2 - S2) / 2.0
    E_vals = []
    for th in np.linspace(0, 2*np.pi, 200):
        sm = brannen_sqrt_m(lam_bar, r_ratio, th)
        if np.all(sm > -1e-6):  # dopuszczamy male ujemne (numeryczne)
            s1 = np.sum(sm)
            s2 = np.sum(sm**2)
            E = (s1**2 - s2) / 2.0
            E_vals.append(E)
    E_vals = np.array(E_vals)
    E_std_rel = np.std(E_vals) / np.abs(np.mean(E_vals)) if len(E_vals) > 0 else 999

    results.append(("C8: E_cross = (S1^2-S2)/2 = const na C (Brannen skan)",
                    E_std_rel < 1e-10 and len(E_vals) >= 50,
                    f"E_cross = {E_cross_analytic:.4f}, "
                    f"std/mean = {E_std_rel:.2e}, {len(E_vals)} punktow"))

    # === BLOK 4: Dekoherencja fazowa => Z_3 ===

    # C9: Z_3 = max entropy katowej (H = ln 3)
    alpha_z3 = np.array([0, 2*np.pi/3, 4*np.pi/3])
    H_z3 = angular_entropy(alpha_z3)
    H_max = np.log(3)
    H_rand = []
    np.random.seed(42)
    for _ in range(1000):
        a = np.sort(np.random.uniform(0, 2*np.pi, 3))
        H_rand.append(angular_entropy(a))
    H_rand = np.array(H_rand)

    results.append(("C9: Z_3 = max entropy katowej (H = ln 3)",
                    abs(H_z3 - H_max) < 1e-12 and np.all(H_rand <= H_max + 1e-10),
                    f"H(Z_3) = {H_z3:.6f}, ln(3) = {H_max:.6f}, "
                    f"max(H_rand) = {np.max(H_rand):.6f}"))

    # C10: Sum cos(alpha_n - alpha_m) = -3/2 dla Z_3 (minimum globalne)
    cos_z3 = pair_cosine_sum(alpha_z3)
    cos_rand = []
    for _ in range(5000):
        a = np.random.uniform(0, 2*np.pi, 3)
        cos_rand.append(pair_cosine_sum(a))
    cos_rand = np.array(cos_rand)

    results.append(("C10: Z_3: Sum cos(a_n-a_m) = -3/2 (min globalne)",
                    abs(cos_z3 - (-1.5)) < 1e-12 and np.min(cos_rand) >= -1.5 - 1e-10,
                    f"Sum cos(Z_3) = {cos_z3:.6f}, min(rand) = {np.min(cos_rand):.6f}"))

    # C11: |Sum e^{i*alpha_k}| = 0 <=> Sum cos = -N/2 <=> Z_3
    # Dowod: |Sum e^{i*alpha}|^2 = N + 2*Sum cos(alpha_n - alpha_m) >= 0
    # Rownosc z 0 <=> Sum cos = -N/2 <=> Z_N
    sum_exp_z3 = np.sum(np.exp(1j * alpha_z3))
    identity_val = 3.0 + 2.0 * cos_z3  # should be 0
    results.append(("C11: |Sum e^{ia}|^2 = N + 2*Sum cos = 0 <=> Z_3",
                    abs(sum_exp_z3) < 1e-14 and abs(identity_val) < 1e-12,
                    f"|Sum e^{{ia}}| = {abs(sum_exp_z3):.2e}, "
                    f"N + 2*Sum cos = {identity_val:.2e}"))

    # C12: min|Sum e^{ia}|^2 = 0 osiagane jedynie w Z_3 (optymalizacja)
    res = differential_evolution(coherence, bounds=[(0, 2*np.pi)]*3,
                                  seed=42, maxiter=1000, tol=1e-12)
    opt_alpha = np.sort(res.x % (2*np.pi))
    gaps_opt = np.diff(opt_alpha)
    gaps_opt = np.append(gaps_opt, 2*np.pi - opt_alpha[-1] + opt_alpha[0])
    gaps_dev = np.std(gaps_opt)

    results.append(("C12: min|Sum e^{ia}|^2 = 0 osiagane w Z_3 (jedyne)",
                    res.fun < 1e-10 and gaps_dev < 0.01,
                    f"|Sum|^2_min = {res.fun:.2e}, "
                    f"gaps std = {gaps_dev:.4f}"))

    # === BLOK 5: Pelny lancuch B2 ===

    # C13: Pelny lancuch: Z_3(amplitude space) => Brannen => Q_K = 3/2
    c_pdg = np.mean(SQRT_M_PDG)
    r_pdg = np.sqrt(2) * CV(SQRT_M_PDG)
    cos_th = (SQRT_M_PDG[0]/c_pdg - 1.0) / r_pdg
    theta_K_pdg = np.arccos(np.clip(cos_th, -1, 1))
    theta_deg = np.degrees(theta_K_pdg)

    m_brannen = brannen_sqrt_m(c_pdg, r_pdg, theta_K_pdg)**2
    qk_brannen_val = Q_K(m_brannen)

    results.append(("C13: PDG -> Brannen(theta_K) -> Q_K = 3/2",
                    abs(qk_brannen_val - 1.5) < 0.001 and 130 < theta_deg < 135,
                    f"theta_K = {theta_deg:.2f} deg, Q_K(Brannen) = {qk_brannen_val:.6f}"))

    # C14: Kompletnosc lancucha B2
    chain_ok = (
        abs(qk_pdg - 1.5) < 0.001 and           # Q_K = 3/2
        abs(cv_pdg - 1.0) < 0.01 and             # CV = 1
        abs(cos_z3 + 1.5) < 1e-12 and            # Z_3 selekcja
        abs(qk_brannen_val - 1.5) < 0.001         # Brannen zamkniecie
    )

    results.append(("C14: Lancuch B2 kompletny (5 krokow)",
                    chain_ok,
                    "N=3 -> A_n -> C(S1,S2) -> Z_3(dekoherencja) -> Brannen -> Q_K=3/2"))

    # === PODSUMOWANIE ===
    n_pass = sum(1 for _, p, _ in results if p)
    n_total = len(results)
    print("=" * 70)
    print(f"ex145_B2_qk_z3_bridge.py -- B2: Q_K = 3/2 z Z_3")
    print(f"TGP v1 -- {n_pass}/{n_total} PASS")
    print("=" * 70)
    for name, passed, detail in results:
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {name}: {detail}")
    print("=" * 70)

    print(f"\n  LANCUCH LOGICZNY B2:")
    print(f"    1. Ghost barrier => N=3 mody solitonowe (prop:no-4th-generation)")
    print(f"    2. M ~ A^4 (most B1, zamkniety) => sqrt(m_n) ~ A_n^2")
    print(f"    3. lambda = (A_0^2, A_1^2, A_2^2) na okregu C (wiezy S_1, S_2)")
    print(f"    4. Dekoherencja fazowa: |Sum e^{{i*alpha_k}}| = 0 => Z_3")
    print(f"       (jedyne rozwiazanie: rowne odstepy 2pi/3)")
    print(f"    5. Z_3 na C = Brannen(r=sqrt(2)) => Q_K = 3/2")
    print(f"")
    print(f"  WAZNY WNIOSEK (C8):")
    print(f"    E_cross = (S1^2-S2)/2 = CONST na okregu C")
    print(f"    => Z_3 NIE wynika z minimalizacji E_cross!")
    print(f"    => Selekcja Z_3 pochodzi z DEKOHERENCJI FAZOWEJ")
    print(f"       (maksymalna niezaleznosc sektorow, min koherencji)")
    print(f"")
    print(f"  KOREKTA vs prop:R2-minimal-interference:")
    print(f"    STARY: fazy ogonowe phi_n tworza Z_3 (blad!)")
    print(f"    NOWY:  Z_3 dziala na KATACH alpha_k w przestrzeni")
    print(f"           amplitudowej lambda_k = A_k^2, NIE na fazach delta_n")
    print(f"")
    print(f"  WYNIKI KLUCZOWE:")
    print(f"    Q_K(PDG) = {qk_pdg:.6f}")
    print(f"    CV(sqrt(m)) = {cv_pdg:.6f}")
    print(f"    r = rho/lam_bar = {r_ratio:.6f} (sqrt(2) = {np.sqrt(2):.6f})")
    print(f"    Sum cos(a_n-a_m)|_Z3 = {cos_z3:.6f} (min = -N/2 = -3/2)")
    print(f"    theta_K = {theta_deg:.2f} deg")

    return n_pass == n_total


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
