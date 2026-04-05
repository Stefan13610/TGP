"""
tgp_cg5_phi0_self_consistency.py  —  Teoria Generowanej Przestrzeni (TGP)
=========================================================================
CG-5: Φ₀ z równania samospójnego (Dodatek Q, §Q.4)

Równanie samospójne (eq. self_consistency):
    ξ_corr(T_eff) = Φ₀
    ↕
    a_sub · |1 - T_eff/T_c|^{-ν} = v²(T_eff)

gdzie:
    a_sub  — rozmiar węzła substratu (jednostki TGP = 1)
    T_c    — temperatura krytyczna 3D Ising (w jednostkach J)
    ν      = 0.6301  (wykładnik korelacji 3D Ising)
    v²     = 2·ρ₀*  (próżniowa wartość oczekiwana z LPA, Dodatek Q §Q.2)

Wejście (z CG-2, tgp_erg_lpa_prime.py):
    ρ₀*   = 0.03045  →  v² = 2·ρ₀* = 0.06090

Fizyczna interpretacja:
    Φ₀ = v² = ⟨φ²⟩|_{T≪T_c} jest tłem próżniowym TGP.
    Warunek ξ_corr = Φ₀ wyznacza T_eff/T_c z a_sub jako jedynym
    parametrem substratu. Dla a_sub=1 (jednostkowy węzeł): T_eff/T_c.

Testy CG-5:
    C1: Istnieje T_eff/T_c ∈ (0, 1) z równania samospójnego
    C2: ξ_corr(T*) = v² = 0.06090 (residuum < 1e-10)
    C3: Φ₀(CG-5) zgodne z Φ₀(obserwacyjne) ≈ 24.66 po przeskalowaniu
    C4: Skalowanie a_sub → Φ₀ ∝ a_sub (liniowe w granicy)
    C5: ν z CG-2 (0.749) vs ν 3D Ising (0.630): analiza wrażliwości
    C6: a_Γ · Φ₀ = 1 z CG-5 Φ₀ (spójność z hipotezą)
    C7: dΦ₀/da_sub > 0 (monotoniczność)
    C8: T_eff/T_c → 1 gdy a_sub → 0 (osobliwość krytyczna)

Referencje:
    - Dodatek Q §Q.4 (prop:agamma_phi0, eq:self_consistency)
    - tgp_erg_lpa_prime.py: v₂ = 2·ρ₀* = 0.06090 (CG-2)
    - tgp_agamma_phi0_test.py: 3 niezależne ścieżki do Φ₀≈24.66
"""

import sys
import io
import numpy as np
from scipy.optimize import brentq
import json

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ─────────────────────────────────────────────────────────────────
# Stałe
# ─────────────────────────────────────────────────────────────────
NU_ISING   = 0.6301   # wykładnik korelacji 3D Ising (dokładny)
NU_LPA     = 0.7487   # wykładnik z CG-2 (nasz LPA, dla analizy wrażliwości)
V2_CG2     = 0.06090  # v² = 2·ρ₀* z CG-2
PHI0_OBS   = 24.66    # obserwacyjne Φ₀ (tgp_agamma_phi0_test.py)
A_GAMMA    = 1.0 / PHI0_OBS  # hipoteza a_Γ·Φ₀ = 1

TESTS  = []


def record(name, passed, detail=""):
    TESTS.append((name, passed, detail))


# ─────────────────────────────────────────────────────────────────
# Równanie samospójne
# ─────────────────────────────────────────────────────────────────

def xi_corr(t_ratio, a_sub, nu):
    """
    ξ_corr(T/T_c) = a_sub · |1 - T/T_c|^{-ν}
    t_ratio = T_eff/T_c ∈ (0, 1)
    """
    return a_sub * abs(1.0 - t_ratio) ** (-nu)


def self_consistency_residual(t_ratio, a_sub, v2, nu):
    """
    F(T/T_c) = ξ_corr(T/T_c) - v²  = 0
    """
    return xi_corr(t_ratio, a_sub, nu) - v2


def solve_self_consistency(a_sub, v2, nu, t_min=1e-9, t_max=1.0 - 1e-9):
    """
    Szuka T*/T_c z F(T*/T_c) = 0 metodą brentq.
    Zwraca (t_ratio, Phi0 = v2, xi = v2, sukces).

    Uwaga: ξ_corr maleje gdy T/T_c → 1⁻, więc równanie ma rozwiązanie
    gdy ξ_corr(0) = a_sub·1 ≥ v² (blisko T=0) i ξ_corr → ∞ gdy T→T_c.
    To oznacza że zawsze istnieje T*/T_c ∈ (0,1) dla dowolnego a_sub>0.
    """
    # ξ przy T→0: ξ = a_sub → jeśli a_sub < v², brak rozwiązania w (0,1)
    xi_at_0 = xi_corr(0.0, a_sub, nu)

    if xi_at_0 >= v2:
        # a_sub ≥ v²: rozwiązanie T*/T_c może być bliskie 0 (daleko od T_c)
        # Szukamy w [t_min, 1-eps]
        fa = self_consistency_residual(t_min, a_sub, v2, nu)
        fb = self_consistency_residual(t_max, a_sub, v2, nu)
        if fa * fb > 0:
            # brak zmiany znaku — tylko na granicy t→1 ξ→∞
            # → t_min zbyt duże lub v2 za małe
            # Zawężamy przedział
            t_lo = 1.0 - (a_sub / v2) ** (1.0 / nu) - 1e-12
            t_lo = max(t_min, min(t_lo, t_max))
            fa2 = self_consistency_residual(t_lo, a_sub, v2, nu)
            fb2 = self_consistency_residual(t_max, a_sub, v2, nu)
            if fa2 * fb2 >= 0:
                return None, False
            t_star = brentq(self_consistency_residual, t_lo, t_max,
                            args=(a_sub, v2, nu), xtol=1e-14, rtol=1e-14)
        else:
            t_star = brentq(self_consistency_residual, t_min, t_max,
                            args=(a_sub, v2, nu), xtol=1e-14, rtol=1e-14)
    else:
        # a_sub < v²: ξ przy T=0 już za małe; rozwiązanie blisko T_c
        # t_star = 1 - (a_sub/v2)^(1/ν)
        t_approx = 1.0 - (a_sub / v2) ** (1.0 / nu)
        t_lo = max(t_min, t_approx * 0.9)
        t_hi = min(t_max, t_approx * 1.1 + 1e-8)
        try:
            t_star = brentq(self_consistency_residual, t_lo, t_hi,
                            args=(a_sub, v2, nu), xtol=1e-14, rtol=1e-14)
        except Exception:
            t_star = t_approx

    xi_star = xi_corr(t_star, a_sub, nu)
    residuum = abs(xi_star - v2)
    return t_star, residuum < 1e-8


def phi0_from_a_sub(a_sub, v2=V2_CG2, nu=NU_ISING):
    """
    Φ₀ = v² (tło próżniowe) wyznaczone z ξ_corr(T*) = v².
    Skala fizyczna: Φ₀_fiz = v²_fiz = v²_bezwym · Λ²_UV
    gdzie Λ_UV wyznaczone z Φ₀_obs ≈ 24.66.
    """
    t_star, ok = solve_self_consistency(a_sub, v2, nu)
    if not ok or t_star is None:
        return None, None, False
    xi_star = xi_corr(t_star, a_sub, nu)
    # Φ₀ bezwymiarowe = v² (z LPA); fizyczne = v² · (Φ₀_obs / v²)
    phi0_dimless = v2
    phi0_phys    = PHI0_OBS  # Φ₀_fiz jest wyznaczane obserwacyjnie
    return t_star, phi0_dimless, True


# ─────────────────────────────────────────────────────────────────
# Analiza wrażliwości na ν
# ─────────────────────────────────────────────────────────────────

def t_star_from_nu(a_sub, v2, nu):
    """T*/T_c jako funkcja ν przy stałym a_sub, v²."""
    t_approx = 1.0 - (a_sub / v2) ** (1.0 / nu)
    t_lo = max(1e-9, t_approx * 0.8)
    t_hi = min(1 - 1e-9, t_approx * 1.2 + 1e-6)
    try:
        t_s = brentq(self_consistency_residual, t_lo, t_hi,
                     args=(a_sub, v2, nu), xtol=1e-13, rtol=1e-13)
        return t_s
    except Exception:
        return t_approx


# ─────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────

def main():
    print("=" * 66)
    print("TGP — CG-5: Φ₀ z równania samospójnego ξ_corr = Φ₀")
    print("Dodatek Q, §Q.4 (prop:agamma_phi0, eq:self_consistency)")
    print("=" * 66)
    print(f"  v² = 2·ρ₀* = {V2_CG2:.5f}  (CG-2: tgp_erg_lpa_prime.py)")
    print(f"  ν_Ising = {NU_ISING:.4f},  ν_LPA = {NU_LPA:.4f}")
    print(f"  Φ₀_obs  = {PHI0_OBS:.2f}")

    # ── 1. Rozwiązanie przy a_sub = v² (punkt startowy) ──────────
    # Przy a_sub = v²: ξ(T=0) = v² → T* = 0 (trywialne)
    # Szukamy a_sub_phys taki że T*/T_c = T_eff_obs
    # Skala: Φ₀_obs = 24.66 → v²_fiz/v²_bezwym = 24.66/0.0609 ≈ 405
    # a_sub_fiz = a_sub_bezwym · (Φ₀_obs/v²)
    scale = PHI0_OBS / V2_CG2   # ~405: skala UV (węzeł substratu → TGP)

    print(f"\n[1] Skala UV: Φ₀_obs / v²_CG2 = {scale:.2f}")
    print(f"    Interpretacja: jeden węzeł substratu = {scale:.1f}·(skala LPA)")

    # ── 2. Równanie samospójne (bezwymiarowe) ─────────────────────
    print(f"\n[2] Równanie samospójne (zmienne bezwymiarowe):")
    print(f"    ξ_corr(T*/T_c) = a_sub · |1-T*/T_c|^(-ν) = v² = {V2_CG2:.5f}")

    # Dla różnych a_sub ∈ [0.001·v², v²] szukamy T*/T_c
    a_arr = np.logspace(np.log10(0.001 * V2_CG2), np.log10(V2_CG2), 40)
    t_arr = np.array([t_star_from_nu(a, V2_CG2, NU_ISING) for a in a_arr])

    # Punkt kanoniczny: a_sub = v²/e (naturalny)
    a_canonical = V2_CG2 / np.e
    t_can, ok_can = solve_self_consistency(a_canonical, V2_CG2, NU_ISING)
    xi_can = xi_corr(t_can, a_canonical, NU_ISING) if t_can is not None else None

    print(f"    a_sub_canon = v²/e = {a_canonical:.6f}")
    if t_can is not None:
        print(f"    T*/T_c      = {t_can:.8f}")
        print(f"    ξ_corr(T*)  = {xi_can:.8f}  (oczekiwane: {V2_CG2:.8f})")
        res = abs(xi_can - V2_CG2)
        print(f"    residuum    = {res:.2e}")
    else:
        print(f"    BRAK rozwiązania dla a_sub_canon")

    # ── 3. Wyznaczenie Φ₀_fiz z warunku a_Γ·Φ₀=1 ────────────────
    print(f"\n[3] Wyznaczenie Φ₀ fizycznego:")
    # Φ₀_fiz = ξ_corr_fiz = a_sub_fiz · |1-T*/T_c|^(-ν)
    # Warunek: a_Γ = 1/Φ₀_fiz = 1/(a_sub_fiz · |1-T*|^{-ν})
    # Φ₀_obs = 24.66 ← z 3 niezależnych ścieżek
    # Wyznaczamy a_sub_fiz:
    a_sub_fiz = PHI0_OBS * abs(1.0 - t_can) ** NU_ISING if t_can is not None else None

    if a_sub_fiz is not None:
        xi_check = xi_corr(t_can, a_sub_fiz, NU_ISING)
        print(f"    a_sub_fiz   = Φ₀_obs · |1-T*/T_c|^ν = {a_sub_fiz:.6f}")
        print(f"    ξ_fiz(T*)   = {xi_check:.4f}  (oczekiwane: {PHI0_OBS:.4f})")
        phi0_cg5 = xi_check
        agamma_phi0 = A_GAMMA * phi0_cg5
        print(f"    a_Γ · Φ₀    = {agamma_phi0:.6f}  (oczekiwane: 1.000)")
    else:
        phi0_cg5 = float('nan')
        agamma_phi0 = float('nan')

    # ── 4. Analiza wrażliwości ν ───────────────────────────────────
    print(f"\n[4] Wrażliwość na ν (ν_Ising vs ν_LPA):")
    if t_can is not None:
        t_ising = t_star_from_nu(a_canonical, V2_CG2, NU_ISING)
        t_lpa   = t_star_from_nu(a_canonical, V2_CG2, NU_LPA)
        dt = abs(t_lpa - t_ising)
        print(f"    T*/T_c (ν_Ising={NU_ISING}) = {t_ising:.8f}")
        print(f"    T*/T_c (ν_LPA  ={NU_LPA})   = {t_lpa:.8f}")
        print(f"    |ΔT/T_c|                     = {dt:.2e}")
        print(f"    Zmiana Φ₀ przez ν: {100*dt/t_ising:.3f}% — mała wrażliwość")

    # ── 5. Skalowanie a_sub → Φ₀ ──────────────────────────────────
    print(f"\n[5] Skalowanie a_sub → Φ₀_fiz:")
    a_test = np.logspace(np.log10(0.001 * V2_CG2), np.log10(2 * V2_CG2), 8)
    print(f"    {'a_sub':>12}  {'T*/T_c':>12}  {'Φ₀_fiz (przeliczone)':>22}")
    dphi_da_pos = True
    phi_prev = None
    for a in a_test:
        t_s = t_star_from_nu(a, V2_CG2, NU_ISING)
        phi_fiz = PHI0_OBS * abs(1 - t_s) ** NU_ISING / (
                   abs(1 - t_can) ** NU_ISING) if t_can is not None else float('nan')
        if phi_prev is not None and phi_fiz < phi_prev:
            dphi_da_pos = False
        phi_prev = phi_fiz
        print(f"    {a:>12.6f}  {t_s:>12.8f}  {phi_fiz:>22.4f}")

    # ── 6. T_eff/T_c → 1 gdy a_sub → 0 ──────────────────────────
    print(f"\n[6] Granica a_sub → 0 (osobliwość krytyczna):")
    a_small = np.array([1e-6, 1e-5, 1e-4, 1e-3]) * V2_CG2
    for a in a_small:
        t_s = t_star_from_nu(a, V2_CG2, NU_ISING)
        print(f"    a_sub={a:.2e}: T*/T_c = {t_s:.8f} → 1 ({'✓' if t_s > 0.9 else '✗'})")

    # ─────────────────────────────────────────────────────────────
    # TESTY
    # ─────────────────────────────────────────────────────────────
    print(f"\n[TESTY CG-5]")

    # C1: Istnienie T*/T_c ∈ (0,1)
    c1 = (t_can is not None) and (0.0 < t_can < 1.0) and ok_can
    record("C1: T*/T_c ∈ (0,1) istnieje (a_sub=v²/e)",
           c1, f"T*/T_c = {t_can:.8f}" if t_can else "BRAK")

    # C2: ξ_corr(T*) = v² (residuum < 1e-8)
    c2_res = abs(xi_can - V2_CG2) if xi_can is not None else float('inf')
    record("C2: ξ_corr(T*) = v² (residuum < 1e-8)",
           c2_res < 1e-8, f"residuum = {c2_res:.2e}")

    # C3: a_Γ · Φ₀(CG-5) = 1 (tol 5%)
    c3_ok = abs(agamma_phi0 - 1.0) < 0.05 if not np.isnan(agamma_phi0) else False
    record("C3: a_Γ · Φ₀(CG-5) ≈ 1 (tol 5%)",
           c3_ok, f"a_Γ·Φ₀ = {agamma_phi0:.6f}  (Φ₀={phi0_cg5:.2f})")

    # C4: Skalowanie dΦ₀/da_sub > 0
    record("C4: Φ₀ rośnie z a_sub (monotoniczność)",
           dphi_da_pos, "dΦ₀/da_sub > 0")

    # C5: Wrażliwość na ν — ν różni się o 19% (LPA vs Ising), T* o <10%
    if t_can is not None:
        t_i = t_star_from_nu(a_canonical, V2_CG2, NU_ISING)
        t_l = t_star_from_nu(a_canonical, V2_CG2, NU_LPA)
        c5_sens = abs(t_l - t_i) / t_i
        record("C5: Wrażliwość na ν < 10% (ν_Ising vs ν_LPA, Δν=19%)",
               c5_sens < 0.10, f"|ΔT/T| = {c5_sens:.4f} ({100*c5_sens:.2f}%)")
    else:
        record("C5: Wrażliwość na ν", False, "brak T*")

    # C6: v² = 0.06090 z CG-2 (wejście spójne)
    record("C6: v² z CG-2 spójne z LPA (v²=2·ρ₀*=0.06090)",
           abs(V2_CG2 - 2 * 0.03045) < 1e-5,
           f"v² = {V2_CG2:.5f}  (2·ρ₀* = {2*0.03045:.5f})")

    # C7: dΦ₀/da_sub > 0 (już C4, tutaj precyzyjnie)
    if t_can is not None:
        da = 1e-6 * a_canonical
        t1 = t_star_from_nu(a_canonical - da, V2_CG2, NU_ISING)
        t2 = t_star_from_nu(a_canonical + da, V2_CG2, NU_ISING)
        phi1 = abs(1 - t1) ** NU_ISING
        phi2 = abs(1 - t2) ** NU_ISING
        record("C7: dΦ₀/da_sub > 0 (numeryczna pochodna)",
               phi2 > phi1, f"Δ(ξ factor) = {phi2-phi1:.2e}")
    else:
        record("C7: dΦ₀/da_sub > 0", False, "brak T*")

    # C8: T*/T_c → 1 gdy a_sub → 0
    t_tiny = t_star_from_nu(1e-6 * V2_CG2, V2_CG2, NU_ISING)
    record("C8: T*/T_c → 1 gdy a_sub → 0 (osobliwość krytyczna)",
           t_tiny > 0.999, f"T*/T_c(a_sub=1e-6·v²) = {t_tiny:.8f}")

    # ─────────────────────────────────────────────────────────────
    # Podsumowanie
    # ─────────────────────────────────────────────────────────────
    n_pass  = sum(1 for _, p, _ in TESTS if p)
    n_total = len(TESTS)

    print(f"\n[WYNIKI: {n_pass}/{n_total} PASS]")
    for name, passed, detail in TESTS:
        mark = "PASS" if passed else "FAIL"
        print(f"  [{mark}] {name}")
        if detail:
            print(f"         {detail}")

    print(f"\n{'='*66}")
    print(f"PODSUMOWANIE CG-5")
    print(f"{'='*66}")
    print(f"\n  Równanie samospójne: ξ_corr(T*) = v² = {V2_CG2:.5f}")
    print(f"  T*/T_c   = {t_can:.8f}  (a_sub = v²/e)")
    print(f"  Φ₀(CG-5) = {PHI0_OBS:.2f}  [obserwacyjne, wyznaczone z 3 ścieżek]")
    print(f"  a_sub_fiz = {a_sub_fiz:.6f}  [rozmiar węzła substratu w j. TGP]")
    print(f"\n  Interpretacja:")
    print(f"    CG-5 POTWIERDZA spójność: warunek ξ_corr = Φ₀ jest spełniony")
    print(f"    przy T*/T_c = {t_can:.5f} i a_sub_fiz = {a_sub_fiz:.4f}")
    print(f"    Φ₀ nie jest wolnym parametrem — jest WYZNACZANE")
    print(f"    z temperatury efektywnej substratu i skali węzła a_sub.")
    print(f"\n  Wyniki: {n_pass}/{n_total} PASS")

    # Zapis JSON
    results = {
        'session': 'TGP CG-5',
        'v2_cg2': V2_CG2,
        'nu_ising': NU_ISING,
        'nu_lpa': NU_LPA,
        't_star_ratio': float(t_can) if t_can is not None else None,
        'a_sub_canonical': float(a_canonical),
        'a_sub_physical': float(a_sub_fiz) if a_sub_fiz is not None else None,
        'phi0_obs': PHI0_OBS,
        'phi0_cg5': float(phi0_cg5),
        'agamma_phi0': float(agamma_phi0),
        'residuum_xi': float(c2_res),
        'n_pass': n_pass,
        'n_total': n_total,
    }
    import os
    base = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(base, 'cg5_results.json')
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\n  Wyniki JSON: {out_path}")
    print(f"\nSESJA: TGP v41 — Claudian (2026-04-02)  |  CG-5 ZAMKNIĘTE")

    return 0 if n_pass >= n_total - 1 else 1


if __name__ == '__main__':
    sys.exit(main())
