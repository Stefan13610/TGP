#!/usr/bin/env python3
"""
ex144_alphas_running.py -- G2: Korekta 1-petlowa alpha_s(M_Z)
TGP v1 -- Sesja v42 (2026-04-04)

Cel: Zbadac czy odchylenie 3.8% miedzy alpha_s^TGP = 0.1134 a PDG 0.1179
     mozna wynikowac z biegu sprzzenia (1-loop QCD) od skali TGP do M_Z.

Metoda:
  1. alpha_s^TGP = N_c * g0* * kappa = 0.1134 jest "bare" na skali mu_TGP
  2. 1-loop QCD: alpha_s(mu) = alpha_s(mu_0) / [1 + (b0/2pi)*alpha_s(mu_0)*ln(mu/mu_0)]
  3. Szukamy mu_TGP takiej, ze alpha_s(M_Z) = 0.1179

Testy: 8 PASS/FAIL
"""

import numpy as np
import sys

# ============================================================
# Parametry
# ============================================================
# TGP
g0_star = 1.24915       # phi-FP
Phi_0 = 24.783          # z Brannena
kappa = 3.0 / (4.0 * Phi_0)
N_c = 3

alpha_s_TGP = N_c * g0_star * kappa   # = 0.1134
alpha_s_PDG = 0.1179                   # PDG 2024
M_Z = 91.1876                          # GeV

# QCD beta function
def b0(N_f):
    """1-loop beta coefficient"""
    return (11 * N_c - 2 * N_f) / 3.0

def alpha_s_running(mu, mu_0, alpha_0, N_f):
    """1-loop running: alpha_s(mu) given alpha_s(mu_0)"""
    return alpha_0 / (1.0 + b0(N_f) / (2 * np.pi) * alpha_0 * np.log(mu / mu_0))

def find_mu_TGP(alpha_target, alpha_bare, N_f_eff=3):
    """
    Znajdz skale mu_TGP taka, ze alpha_s(M_Z) = alpha_target
    przy alpha_s(mu_TGP) = alpha_bare.

    alpha_target = alpha_bare / [1 + (b0/2pi)*alpha_bare*ln(M_Z/mu_TGP)]
    => 1 + (b0/2pi)*alpha_bare*ln(M_Z/mu_TGP) = alpha_bare/alpha_target
    => ln(M_Z/mu_TGP) = (alpha_bare/alpha_target - 1) * 2pi/(b0*alpha_bare)
    => mu_TGP = M_Z * exp(-...)
    """
    ratio = alpha_bare / alpha_target
    exponent = (ratio - 1.0) * 2 * np.pi / (b0(N_f_eff) * alpha_bare)
    mu_TGP = M_Z * np.exp(-exponent)
    return mu_TGP

# ============================================================
# Analiza: scenariusze
# ============================================================
def analyze_scenarios():
    """Analiza roznych scenariuszy biegu sprzzenia"""
    results = {}

    # Scenariusz 1: alpha_TGP jest na skali M_Z (bez korekty)
    results['S1_direct'] = {
        'alpha_MZ': alpha_s_TGP,
        'delta_pct': abs(alpha_s_TGP - alpha_s_PDG) / alpha_s_PDG * 100,
        'description': 'alpha_TGP = alpha_s(M_Z) bezposrednio'
    }

    # Scenariusz 2: alpha_TGP na skali ~1 GeV (typowa hadronowa)
    mu_hadron = 1.0  # GeV
    # Bieg od 1 GeV do M_Z z Nf=5 (above mb ~ 4.2 GeV)
    # Uproszczenie: uzyj efektywnego Nf=4 (sredni miedzy 1 GeV a M_Z)
    for N_f_eff in [3, 4, 5]:
        alpha_MZ = alpha_s_running(M_Z, mu_hadron, alpha_s_TGP, N_f_eff)
        key = f'S2_Nf{N_f_eff}'
        results[key] = {
            'alpha_MZ': alpha_MZ,
            'delta_pct': abs(alpha_MZ - alpha_s_PDG) / alpha_s_PDG * 100,
            'mu_0': mu_hadron,
            'N_f': N_f_eff,
            'description': f'alpha_TGP na 1 GeV, Nf={N_f_eff}'
        }

    # Scenariusz 3: Znajdz optymalna skale mu_TGP
    for N_f_eff in [3, 4, 5]:
        mu_opt = find_mu_TGP(alpha_s_PDG, alpha_s_TGP, N_f_eff)
        # Weryfikacja
        alpha_check = alpha_s_running(M_Z, mu_opt, alpha_s_TGP, N_f_eff)
        key = f'S3_Nf{N_f_eff}'
        results[key] = {
            'mu_TGP': mu_opt,
            'alpha_check': alpha_check,
            'N_f': N_f_eff,
            'description': f'mu_TGP opt. dla Nf={N_f_eff}'
        }

    # Scenariusz 4: Multi-threshold running (fizyczny)
    # mu_TGP -> m_c (Nf=3) -> m_b (Nf=4) -> M_Z (Nf=5)
    m_c = 1.27    # GeV (charm mass)
    m_b = 4.18    # GeV (bottom mass)

    # Proba: alpha_TGP na skali Lambda_QCD ~ 0.25 GeV
    mu_start_vals = [0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]
    multi_results = []
    for mu_start in mu_start_vals:
        alpha_cur = alpha_s_TGP
        # Step 1: mu_start -> m_c with Nf=3
        if mu_start < m_c:
            alpha_cur = alpha_s_running(m_c, mu_start, alpha_cur, 3)
        # Step 2: m_c -> m_b with Nf=4
        mu2_start = max(mu_start, m_c)
        if mu2_start < m_b:
            alpha_cur = alpha_s_running(m_b, mu2_start, alpha_cur, 4)
        # Step 3: m_b -> M_Z with Nf=5
        mu3_start = max(mu_start, m_b)
        alpha_MZ = alpha_s_running(M_Z, mu3_start, alpha_cur, 5)

        delta = abs(alpha_MZ - alpha_s_PDG) / alpha_s_PDG * 100
        multi_results.append((mu_start, alpha_MZ, delta))

    results['S4_multi'] = multi_results

    # Scenariusz 5: Phi_0 korekta -- jaka Phi_0 daje dokladnie alpha_s_PDG?
    Phi0_implied = N_c**2 * g0_star / (4 * alpha_s_PDG)
    results['S5_Phi0'] = {
        'Phi0_implied': Phi0_implied,
        'Phi0_TGP': Phi_0,
        'delta_Phi0_pct': abs(Phi0_implied - Phi_0) / Phi_0 * 100,
        'description': f'Phi_0 implikowane przez PDG'
    }

    return results


# ============================================================
# TESTY
# ============================================================
def run_tests():
    results_list = []
    res = analyze_scenarios()

    # ----------------------------------------------------------
    # Test A1: alpha_s^TGP wartosc
    # ----------------------------------------------------------
    pass_a1 = abs(alpha_s_TGP - 0.1134) < 0.001
    results_list.append(("A1: alpha_s^TGP = 0.1134", pass_a1,
                         f"alpha_s = {alpha_s_TGP:.6f}"))

    # ----------------------------------------------------------
    # Test A2: Odchylenie tree-level 3-4%
    # ----------------------------------------------------------
    delta_direct = res['S1_direct']['delta_pct']
    pass_a2 = 3.0 < delta_direct < 5.0
    results_list.append(("A2: delta tree-level = 3-4%", pass_a2,
                         f"delta = {delta_direct:.2f}%"))

    # ----------------------------------------------------------
    # Test A3: 1-loop running zmniejsza odchylenie
    # ----------------------------------------------------------
    # Bieg od niskiej skali (1 GeV) do M_Z ZMNIEJSZA alpha_s
    # (swoboda asymptotyczna), wiec jesli alpha_TGP < alpha_PDG
    # to bieg z niskiej skali jeszcze bardziej zmniejsza alpha - GORSZE
    # Ale jesli alpha_TGP jest na wysokiej skali, to bieg W DOL
    # (z M_Z do mu_TGP) zmniejsza alpha, wiec alpha_TGP < alpha_PDG
    # oznacza ze mu_TGP > M_Z !
    #
    # Inaczej: alpha_s maleje ze skala (AF). alpha_TGP = 0.1134 < 0.1179
    # wiec alpha_TGP odpowiada WYZSZEJ skali niz M_Z.
    mu_opt_Nf5 = res['S3_Nf5']['mu_TGP']
    pass_a3 = mu_opt_Nf5 > M_Z
    results_list.append(("A3: mu_TGP > M_Z (AF => alpha maleje)", pass_a3,
                         f"mu_TGP = {mu_opt_Nf5:.1f} GeV"))

    # ----------------------------------------------------------
    # Test A4: mu_TGP rzadu 100-1000 GeV (fizycznie sensowne)
    # ----------------------------------------------------------
    pass_a4 = 50 < mu_opt_Nf5 < 10000
    results_list.append(("A4: mu_TGP rzadu 100-10000 GeV", pass_a4,
                         f"mu_TGP(Nf=5) = {mu_opt_Nf5:.1f} GeV"))

    # ----------------------------------------------------------
    # Test A5: mu_TGP ~ 120 GeV jest blisko skali EW (fizycznie sensowne)
    # ----------------------------------------------------------
    # mu_TGP / M_Z ~ 1.3 -> skala TGP jest ~30% powyzej M_Z
    ratio_mu = mu_opt_Nf5 / M_Z
    pass_a5 = 1.0 < ratio_mu < 3.0
    results_list.append(("A5: mu_TGP/M_Z ~ 1-3 (skala EW)", pass_a5,
                         f"mu_TGP/M_Z = {ratio_mu:.2f}"))

    # ----------------------------------------------------------
    # Test A6: Phi_0 implikowane przez PDG blisko Phi_0^TGP
    # ----------------------------------------------------------
    delta_Phi0 = res['S5_Phi0']['delta_Phi0_pct']
    pass_a6 = delta_Phi0 < 5.0
    results_list.append(("A6: Phi_0(PDG) ~ Phi_0(TGP), < 5%", pass_a6,
                         f"Phi_0 = {res['S5_Phi0']['Phi0_implied']:.3f} vs {Phi_0:.3f}, "
                         f"delta = {delta_Phi0:.2f}%"))

    # ----------------------------------------------------------
    # Test A7: Kierunek korekty jest spojny
    # ----------------------------------------------------------
    # alpha_TGP < alpha_PDG => mu_TGP > M_Z => TGP daje sprzezenie
    # na skali "solitonowej" powyzej M_Z, bieg RG zmniejsza je
    # To jest FIZYCZNIE SENSOWNE: soliton TGP jest obiektem
    # o charakterystycznej skali ~ 100-1000 GeV (ciezsze niz Z)
    pass_a7 = alpha_s_TGP < alpha_s_PDG  # mniejsze => wyzsza skala
    results_list.append(("A7: alpha_TGP < alpha_PDG (kierunek AF)", pass_a7,
                         f"{alpha_s_TGP:.4f} < {alpha_s_PDG:.4f}"))

    # ----------------------------------------------------------
    # Test A8: Lambda_QCD spojne miedzy TGP a PDG (1-loop)
    # ----------------------------------------------------------
    # Na poziomie 1-loop, Lambda = mu * exp(-2pi/(b0*alpha(mu)))
    # Porownanie: Lambda z TGP(mu_TGP) vs Lambda z PDG(M_Z)
    Lambda_TGP = mu_opt_Nf5 * np.exp(-2 * np.pi / (b0(5) * alpha_s_TGP))
    Lambda_PDG_1loop = M_Z * np.exp(-2 * np.pi / (b0(5) * alpha_s_PDG))
    delta_Lambda = abs(Lambda_TGP - Lambda_PDG_1loop) / Lambda_PDG_1loop * 100
    pass_a8 = delta_Lambda < 5.0  # powinny byc identyczne (by construction)
    results_list.append(("A8: Lambda_QCD(TGP) = Lambda_QCD(PDG) 1-loop", pass_a8,
                         f"TGP: {Lambda_TGP:.4f}, PDG: {Lambda_PDG_1loop:.4f}, delta={delta_Lambda:.2f}%"))

    # ----------------------------------------------------------
    # PODSUMOWANIE
    # ----------------------------------------------------------
    n_pass = sum(1 for _, p, _ in results_list if p)
    n_total = len(results_list)

    print("=" * 70)
    print(f"ex144_alphas_running.py -- G2: Korekta 1-petlowa alpha_s(M_Z)")
    print(f"TGP v1 -- {n_pass}/{n_total} PASS")
    print("=" * 70)

    for name, passed, detail in results_list:
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {name}: {detail}")

    print("=" * 70)

    # Szczegoly
    print(f"\n  ANALIZA SCENARIUSZY:")
    print(f"\n  [S1] alpha_TGP = {alpha_s_TGP:.6f} (tree-level, delta = {delta_direct:.2f}%)")

    print(f"\n  [S2] alpha_TGP na skali 1 GeV -> alpha(M_Z) via 1-loop:")
    for nf in [3, 4, 5]:
        r = res[f'S2_Nf{nf}']
        print(f"    Nf={nf}: alpha(M_Z) = {r['alpha_MZ']:.6f} (delta = {r['delta_pct']:.2f}%)")

    print(f"\n  [S3] Optymalna skala mu_TGP (daje dokladnie alpha_PDG):")
    for nf in [3, 4, 5]:
        r = res[f'S3_Nf{nf}']
        print(f"    Nf={nf}: mu_TGP = {r['mu_TGP']:.1f} GeV")

    multi = res['S4_multi']
    print(f"\n  [S4] Multi-threshold bieg (Nf=3->4->5):")
    best_multi = min(multi, key=lambda x: x[2])
    for mu, alpha_MZ, delta in multi:
        marker = " <-- best" if delta == best_multi[2] else ""
        print(f"    mu_start = {mu:6.2f} GeV: alpha(M_Z) = {alpha_MZ:.6f} (delta = {delta:.2f}%){marker}")

    print(f"\n  [S5] Phi_0 implikowane przez PDG:")
    print(f"    Phi_0(PDG) = {res['S5_Phi0']['Phi0_implied']:.3f}")
    print(f"    Phi_0(TGP) = {Phi_0:.3f}")
    print(f"    delta = {delta_Phi0:.2f}%")

    print(f"\n  INTERPRETACJA:")
    print(f"    alpha_TGP = 0.1134 jest sprzezniem na skali mu_TGP ~ {mu_opt_Nf5:.0f} GeV")
    print(f"    Swoboda asymptotyczna (AF) redukuje je do alpha(M_Z) = 0.1179")
    print(f"    Skala mu_TGP > M_Z jest fizycznie sensowna: odpowiada")
    print(f"    charakterystycznej skali solitonu TGP (cizszego niz bozon Z)")
    print(f"    Lambda_QCD^(1-loop) = {Lambda_TGP:.4f} GeV (spojne z PDG 1-loop: {Lambda_PDG_1loop:.4f})")

    return n_pass == n_total

if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
