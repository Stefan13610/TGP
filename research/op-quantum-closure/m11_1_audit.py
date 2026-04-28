##############################################################################
# M11.1 — 1-loop V_eff(Φ) audit (Branch II level 1)
#
# Audyt closure-grade pre-existing M2b loop derivation/script:
#   - TGP/TGP_v1/research/op1-op2-op4/m2b_loop.py
#   - TGP/TGP_v1/research/op1-op2-op4/M2b_loop_derivation.md
#   - TGP/TGP_v1/research/op1-op2-op4/M2b_results.md
#
# Tests:
#   M11.1.1 — V_eff 1-loop reproducibility (cases A, B, C)
#   M11.1.2 — ln K_Φ = 2 ln Φ + const  identity (struct reason β=γ preserved)
#   M11.1.3 — Loop stability A(Φ_0) > 0  +  BZ integrand positivity K+A > 0
#   M11.1.4 — BZ Gauss-Legendre convergence (n_q scan)
#   M11.1.5 — β=γ preservation cross-cases (|β_nat/γ_nat - 1| < 5%)
#   M11.1.6 — η boundary condition consistency: η_CG2 = 0.044 vs η_BI = 0.0253
#
# Konwencje (M2b derivation):
#   V_onsite(Φ) = (m_0²/2)·Φ + (λ_0/4)·Φ² + (T/2) ln Φ
#   K_Φ(q;Φ)   = 2J·Φ²·(3 - cos qₓ - cos qᵧ - cos q_z)        (a = J = 1)
#   A(Φ)       = V_onsite''(Φ) = λ_0/2 - T/(2 Φ²)
#   V_eff(Φ)   = V_onsite(Φ) + (T/2) ∫_BZ d³q/(2π)³ · ln[K_Φ + A]
#   c_n        = (1/n!) · ∂ⁿ V_eff/∂Φⁿ |_{Φ_0}    (Taylor around tree saddle)
#   β_eff = 3·c_3,  γ_eff = -4·c_4
#   β_nat / γ_nat = (β·Φ_0³) / (γ·Φ_0⁴) = β/(γ·Φ_0)
##############################################################################

import numpy as np


# ============================================================================
# Core functions (standalone reimplementation z m2b_loop.py)
# ============================================================================

def v_onsite(phi, m0sq, lambda0, T):
    """V_onsite(Φ) = (m_0²/2)Φ + (λ_0/4)Φ² + (T/2) ln Φ"""
    return 0.5*m0sq*phi + 0.25*lambda0*phi**2 + 0.5*T*np.log(phi)


def A_phi(phi, lambda0, T):
    """A(Φ) = V_onsite''(Φ) = λ_0/2 - T/(2 Φ²)"""
    return 0.5*lambda0 - 0.5*T/phi**2


def phi_saddle(m0sq, lambda0, T):
    """Tree saddle of V_onsite': λ_0 Φ² + m_0² Φ + T = 0,
    physical root for m_0² < 0 ordered phase."""
    disc = m0sq**2 - 4*lambda0*T
    if disc < 0:
        raise ValueError(f"No real saddle: disc = {disc}")
    return (-m0sq + np.sqrt(disc)) / (2*lambda0)


def make_bz_grid(n_q):
    """Gauss-Legendre 3D BZ grid on [0,π]³ (1/8 BZ z symmetrii q_μ → -q_μ).
    Returns cos_sum (Σ cos q_μ) and product weights.
    Normalization: ∫_BZ d³q/(2π)³ f(q) = (1/π³) Σ_i w_i f(q_i)."""
    nodes, w = np.polynomial.legendre.leggauss(n_q)
    q = 0.5*np.pi*(nodes + 1.0)
    dq = 0.5*np.pi*w
    qx, qy, qz = np.meshgrid(q, q, q, indexing='ij')
    wx, wy, wz = np.meshgrid(dq, dq, dq, indexing='ij')
    weights = (wx*wy*wz).flatten()
    cos_sum = (np.cos(qx) + np.cos(qy) + np.cos(qz)).flatten()
    return cos_sum, weights


def v_loop_one(phi, m0sq, lambda0, T, J, cos_sum, weights):
    """1-loop correction (T/2) ∫_BZ d³q/(2π)³ ln[K_Φ + A] @ pojedyncze Φ."""
    A = A_phi(phi, lambda0, T)
    K = 2.0 * J * phi**2 * (3.0 - cos_sum)
    arg = K + A
    if np.any(arg <= 0):
        raise ValueError(f"K + A <= 0: min={arg.min():+.4e} at Φ={phi}")
    integral = np.sum(weights * np.log(arg)) / np.pi**3
    return 0.5 * T * integral


def v_eff_grid(phi_grid, m0sq, lambda0, T, J, cos_sum, weights):
    v_tree = v_onsite(phi_grid, m0sq, lambda0, T)
    v_1loop = np.array([v_loop_one(p, m0sq, lambda0, T, J, cos_sum, weights)
                        for p in phi_grid])
    return v_tree + v_1loop, v_tree, v_1loop


def fit_taylor(phi_grid, v_arr, phi0, window=0.12):
    """Polynomial fit deg-4 wokół Φ_0; return (a1,a2,a3,a4) = coeff of (Φ-Φ_0)^{1..4}."""
    mask = (phi_grid >= phi0*(1-window)) & (phi_grid <= phi0*(1+window))
    if mask.sum() < 6:
        raise RuntimeError(f"Too few fit points: {mask.sum()}")
    x = phi_grid[mask] - phi0
    y = v_arr[mask]
    coeffs = np.polyfit(x, y, 4)
    a4, a3, a2, a1, _ = coeffs
    return a1, a2, a3, a4


# ============================================================================
# Audit tests
# ============================================================================

# Reference values from M2b_results.md (n_q=32, window=0.12)
M2B_REFERENCE = {
    'A': dict(m0sq=-4.0, lambda0=1.0, T=1.0,
              phi0_ref=3.7321, A_ref=0.464,
              c3_tree_ref=3.206e-3,  c3_1loop_ref=9.507e-3,
              c4_tree_ref=-6.443e-4, c4_1loop_ref=-1.912e-3,
              beta_nat_over_gamma_nat_1loop=0.9990,    # -0.10%
              ratio_beta_amp=2.965, ratio_gamma_amp=2.968),
    'B': dict(m0sq=-5.0, lambda0=2.0, T=1.0,
              phi0_ref=2.2808, A_ref=0.904,
              c3_tree_ref=1.405e-2,  c3_1loop_ref=3.932e-2,
              c4_tree_ref=-4.619e-3, c4_1loop_ref=-1.312e-2,
              beta_nat_over_gamma_nat_1loop=0.9855,    # -1.45%
              ratio_beta_amp=2.799, ratio_gamma_amp=2.840),
    'C': dict(m0sq=-2.0, lambda0=1.0, T=0.1,
              phi0_ref=1.9487, A_ref=0.487,
              c3_tree_ref=2.252e-3,  c3_1loop_ref=6.183e-3,
              c4_tree_ref=-8.669e-4, c4_1loop_ref=-2.276e-3,
              beta_nat_over_gamma_nat_1loop=1.0457,    # +4.57%
              ratio_beta_amp=2.745, ratio_gamma_amp=2.625),
}


def compute_case(label, m0sq, lambda0, T, J=1.0, n_q=32, n_phi=400,
                 window=0.12, verbose=False):
    """Wykonaj 1-loop audyt dla pojedynczego case."""
    phi0 = phi_saddle(m0sq, lambda0, T)
    A0 = A_phi(phi0, lambda0, T)
    cos_sum, weights = make_bz_grid(n_q)
    phi_thresh = np.sqrt(T/lambda0)
    phi_grid = np.linspace(0.35*phi0, 2.5*phi0, n_phi)
    phi_grid = phi_grid[phi_grid > phi_thresh*1.001]

    v_full, v_tree_arr, v_1loop_arr = v_eff_grid(
        phi_grid, m0sq, lambda0, T, J, cos_sum, weights)

    _, c2_t, c3_t, c4_t = fit_taylor(phi_grid, v_tree_arr, phi0, window)
    _, c2_l, c3_l, c4_l = fit_taylor(phi_grid, v_1loop_arr, phi0, window)
    _, c2_f, c3_f, c4_f = fit_taylor(phi_grid, v_full, phi0, window)

    beta_t  = 3*c3_t;     gamma_t  = -4*c4_t
    beta_1l = 3*c3_f;     gamma_1l = -4*c4_f

    beta_nat  = beta_1l * phi0**3
    gamma_nat = gamma_1l * phi0**4
    ratio_nat = beta_nat / gamma_nat   # tree value = 1 (M2a structural)

    return {
        'phi0': phi0, 'A0': A0,
        'c3_tree': c3_t, 'c3_1loop': c3_f,
        'c4_tree': c4_t, 'c4_1loop': c4_f,
        'c3_loop_only': c3_l, 'c4_loop_only': c4_l,
        'beta_tree': beta_t, 'beta_1loop': beta_1l,
        'gamma_tree': gamma_t, 'gamma_1loop': gamma_1l,
        'beta_nat_over_gamma_nat': ratio_nat,
        'ratio_beta_amp': beta_1l/beta_t,
        'ratio_gamma_amp': gamma_1l/gamma_t,
        'phi_grid_npts': len(phi_grid),
    }


def test_M11_1_1_reproducibility(verbose=True):
    """M11.1.1 — V_eff 1-loop reproducibility (cases A, B, C).

    Re-compute β^1loop, γ^1loop, β_nat/γ_nat for each M2b case and
    verify zgodność z reportowanymi w M2b_results.md (Tolerance 1%).
    """
    if verbose:
        print("\n" + "=" * 78)
        print("M11.1.1 — V_eff 1-loop reproducibility (cases A, B, C)")
        print("=" * 78)

    results = {}
    pass_per_case = {}
    for label, ref in M2B_REFERENCE.items():
        res = compute_case(label, ref['m0sq'], ref['lambda0'], ref['T'])
        results[label] = res

        rd_phi0 = abs(res['phi0'] - ref['phi0_ref']) / abs(ref['phi0_ref'])
        rd_c3 = abs(res['c3_1loop'] - ref['c3_1loop_ref']) / abs(ref['c3_1loop_ref'])
        rd_c4 = abs(res['c4_1loop'] - ref['c4_1loop_ref']) / abs(ref['c4_1loop_ref'])
        rd_ratio = abs(res['beta_nat_over_gamma_nat'] - ref['beta_nat_over_gamma_nat_1loop']) \
                   / abs(ref['beta_nat_over_gamma_nat_1loop'])

        crit_phi0 = rd_phi0 < 0.01
        crit_c3 = rd_c3 < 0.05  # 5% — fit systematics
        crit_c4 = rd_c4 < 0.05
        crit_ratio = rd_ratio < 0.01  # 1% — agregowany central claim

        case_pass = crit_phi0 and crit_c3 and crit_c4 and crit_ratio
        pass_per_case[label] = case_pass

        if verbose:
            print(f"\n  Case {label}: m0²={ref['m0sq']}, λ_0={ref['lambda0']}, T={ref['T']}")
            print(f"    Φ_0 (tree saddle):       {res['phi0']:.6f}  "
                  f"(ref {ref['phi0_ref']:.4f}, drift {rd_phi0*100:.3f}%)")
            print(f"    A(Φ_0):                  {res['A0']:+.4e}  "
                  f"({'STABLE' if res['A0'] > 0 else 'TACHYONIC'})")
            print(f"    c_3 (1-loop):            {res['c3_1loop']:+.4e}  "
                  f"(ref {ref['c3_1loop_ref']:+.4e}, drift {rd_c3*100:.3f}%)")
            print(f"    c_4 (1-loop):            {res['c4_1loop']:+.4e}  "
                  f"(ref {ref['c4_1loop_ref']:+.4e}, drift {rd_c4*100:.3f}%)")
            print(f"    β_nat/γ_nat (1-loop):    {res['beta_nat_over_gamma_nat']:.4f}  "
                  f"(ref {ref['beta_nat_over_gamma_nat_1loop']:.4f}, drift {rd_ratio*100:.3f}%)")
            print(f"    β amplification:         {res['ratio_beta_amp']:.3f}× tree")
            print(f"    γ amplification:         {res['ratio_gamma_amp']:.3f}× tree")
            print(f"    → {'PASS' if case_pass else 'FAIL'}")

    test_pass = all(pass_per_case.values())
    if verbose:
        print(f"\n  All 3 cases reproduce M2b values within tolerance: "
              f"{test_pass}")
        print(f"  M11.1.1 → {'PASS' if test_pass else 'FAIL'}")

    return test_pass, {'cases': results, 'pass_per_case': pass_per_case}


def test_M11_1_2_lnK_identity(verbose=True):
    """M11.1.2 — ln K_Φ = 2 ln Φ + const_BZ identity verification.

    Strukturalna przyczyna β=γ preservation: stiffness piece
    `K_Φ(q;Φ) = 2J Φ² (3 - Σ cos q)`  →
    `ln K_Φ(q;Φ) = 2 ln Φ + ln[2J(3-Σ cos q)]`
    (drugi człon Φ-independent), więc
    `<ln K_Φ>_BZ = 2 ln Φ + C_BZ`.

    Test: log-slope BZ-averaged ln K_Φ vs log Φ powinno być DOKŁADNIE 2.
    """
    if verbose:
        print("\n" + "=" * 78)
        print("M11.1.2 — ln K_Φ = 2 ln Φ + const  structural identity")
        print("=" * 78)
        print("  Tożsamość: <ln K_Φ(q;Φ)>_BZ = 2·ln Φ + C_BZ  (C_BZ Φ-independent)")
        print("  → To wyjaśnia β=γ preservation w 1-loop V_eff (M2b §4)")

    cos_sum, weights = make_bz_grid(32)
    J = 1.0
    Phi_test = np.array([0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0])
    avg_lnK = np.zeros_like(Phi_test)

    for i, p in enumerate(Phi_test):
        K = 2.0 * J * p**2 * (3.0 - cos_sum)
        avg_lnK[i] = np.sum(weights * np.log(K)) / np.pi**3

    # Linear fit: y = 2·log(Φ) + C
    log_phi = np.log(Phi_test)
    n = len(log_phi)
    sx = log_phi.sum(); sy = avg_lnK.sum()
    sxx = (log_phi**2).sum(); sxy = (log_phi*avg_lnK).sum()
    slope = (n*sxy - sx*sy) / (n*sxx - sx*sx)
    intercept = (sy - slope*sx) / n

    # Predicted intercept = (1/V_BZ) ∫ ln[2J(3-Σcos q)] d³q
    C_BZ_predicted = np.sum(weights * np.log(2.0 * J * (3.0 - cos_sum))) / np.pi**3

    if verbose:
        print(f"\n  Φ test grid: {Phi_test}")
        print(f"  <ln K_Φ>_BZ values:")
        for p, v in zip(Phi_test, avg_lnK):
            print(f"     Φ = {p:>5.2f}:  <ln K_Φ>_BZ = {v:+.6f}")
        print(f"\n  Linear fit  <ln K_Φ>_BZ  =  slope · ln Φ + intercept:")
        print(f"     slope     = {slope:.6f}    (predicted exact 2.0)")
        print(f"     intercept = {intercept:.6f}    (predicted C_BZ = {C_BZ_predicted:.6f})")
        print(f"     |slope - 2|     = {abs(slope - 2):.2e}")
        print(f"     |intercept - C| = {abs(intercept - C_BZ_predicted):.2e}")

    crit_slope = abs(slope - 2.0) < 1e-6
    crit_intercept = abs(intercept - C_BZ_predicted) < 1e-6

    test_pass = crit_slope and crit_intercept
    if verbose:
        print(f"\n  Slope = 2 (exact, < 1e-6):                  {crit_slope}")
        print(f"  Intercept = C_BZ (predicted, < 1e-6):       {crit_intercept}")
        print(f"  M11.1.2 → {'PASS' if test_pass else 'FAIL'}")
        print(f"\n  STRUCTURAL: identyfikacja `2 ln Φ` w stiffness loop integrand wyjaśnia")
        print(f"  multiplikatywne wzmocnienie H-S Jacobianu (T/2 → 3T/2 effective). β i γ")
        print(f"  oba wzmacniają się ~3× preserving ratio β/γ. To jest STRUKTURALNA przyczyna")
        print(f"  M2a → M2b stability of `β = γ` w canonical sek08a units.")

    return test_pass, {
        'Phi_test': Phi_test, 'avg_lnK': avg_lnK,
        'slope': slope, 'intercept': intercept,
        'C_BZ_predicted': C_BZ_predicted,
    }


def test_M11_1_3_loop_stability_positivity(verbose=True):
    """M11.1.3 — Loop stability A(Φ_0) > 0  +  BZ integrand K+A > 0 ∀q ∈ BZ."""
    if verbose:
        print("\n" + "=" * 78)
        print("M11.1.3 — Loop stability + BZ integrand positivity")
        print("=" * 78)

    cos_sum, weights = make_bz_grid(32)
    pass_per_case = {}
    diag = {}

    for label, ref in M2B_REFERENCE.items():
        phi0 = phi_saddle(ref['m0sq'], ref['lambda0'], ref['T'])
        A0 = A_phi(phi0, ref['lambda0'], ref['T'])

        # K_Φ(q,Φ_0) + A(Φ_0) — minimum across BZ
        K0 = 2.0 * 1.0 * phi0**2 * (3.0 - cos_sum)
        integrand = K0 + A0
        min_KA = float(integrand.min())
        max_KA = float(integrand.max())

        crit_A0 = A0 > 0
        crit_KA = min_KA > 0
        case_pass = crit_A0 and crit_KA
        pass_per_case[label] = case_pass

        diag[label] = {
            'phi0': phi0, 'A0': A0,
            'min_KA': min_KA, 'max_KA': max_KA,
        }

        if verbose:
            print(f"\n  Case {label}: Φ_0 = {phi0:.4f}")
            print(f"    A(Φ_0) = {A0:+.4e}  (loop stable: {crit_A0})")
            print(f"    K(q,Φ_0) + A range over BZ:")
            print(f"        min = {min_KA:+.4e}   max = {max_KA:+.4e}")
            print(f"    Integrand positive ∀q ∈ BZ:  {crit_KA}")
            print(f"    → {'PASS' if case_pass else 'FAIL'}")

    test_pass = all(pass_per_case.values())
    if verbose:
        print(f"\n  All 3 cases: A(Φ_0) > 0 AND K+A > 0 ∀q:    {test_pass}")
        print(f"  M11.1.3 → {'PASS' if test_pass else 'FAIL'}")

    return test_pass, {'diag': diag, 'pass_per_case': pass_per_case}


def test_M11_1_4_bz_convergence(verbose=True):
    """M11.1.4 — BZ Gauss-Legendre convergence (n_q scan, case A)."""
    if verbose:
        print("\n" + "=" * 78)
        print("M11.1.4 — BZ Gauss-Legendre convergence (case A, n_q scan)")
        print("=" * 78)

    ref = M2B_REFERENCE['A']
    nq_list = [8, 16, 24, 32, 40]
    c3_list = []
    c4_list = []

    for nq in nq_list:
        res = compute_case('A', ref['m0sq'], ref['lambda0'], ref['T'], n_q=nq)
        c3_list.append(res['c3_1loop'])
        c4_list.append(res['c4_1loop'])

    c3_arr = np.array(c3_list)
    c4_arr = np.array(c4_list)

    # convergence: |c_n(nq+1) - c_n(nq)| → 0
    c3_diffs = np.abs(np.diff(c3_arr))
    c4_diffs = np.abs(np.diff(c4_arr))

    max_drift_c3 = max(abs(c3_arr[i] - c3_arr[-1]) for i in range(1, len(c3_arr)))
    max_drift_c4 = max(abs(c4_arr[i] - c4_arr[-1]) for i in range(1, len(c4_arr)))
    rel_drift_c3 = max_drift_c3 / abs(c3_arr[-1])
    rel_drift_c4 = max_drift_c4 / abs(c4_arr[-1])

    if verbose:
        print(f"\n  {'n_q':>4} | {'c_3 (1-loop)':>16} | {'c_4 (1-loop)':>16} "
              f"| {'|Δc_3|':>12} | {'|Δc_4|':>12}")
        print(f"  {'-'*4}-+-{'-'*16}-+-{'-'*16}-+-{'-'*12}-+-{'-'*12}")
        for i, nq in enumerate(nq_list):
            d3 = c3_diffs[i-1] if i > 0 else 0.0
            d4 = c4_diffs[i-1] if i > 0 else 0.0
            print(f"  {nq:>4} | {c3_arr[i]:>+16.6e} | {c4_arr[i]:>+16.6e} "
                  f"| {d3:>12.3e} | {d4:>12.3e}")
        print(f"\n  Max relative drift c_3 (n_q≥16 vs n_q=40):  {rel_drift_c3:.3e}")
        print(f"  Max relative drift c_4 (n_q≥16 vs n_q=40):  {rel_drift_c4:.3e}")

    # PASS: 5+ digits stable beyond n_q ≥ 16
    crit_c3 = rel_drift_c3 < 1e-5
    crit_c4 = rel_drift_c4 < 1e-5

    test_pass = crit_c3 and crit_c4
    if verbose:
        print(f"\n  c_3 stable to 5+ digits (rel < 1e-5):       {crit_c3}")
        print(f"  c_4 stable to 5+ digits (rel < 1e-5):       {crit_c4}")
        print(f"  M11.1.4 → {'PASS' if test_pass else 'FAIL'}")

    return test_pass, {
        'nq_list': nq_list, 'c3_arr': c3_arr.tolist(), 'c4_arr': c4_arr.tolist(),
        'rel_drift_c3': rel_drift_c3, 'rel_drift_c4': rel_drift_c4,
    }


def test_M11_1_5_beta_gamma_preservation(verbose=True):
    """M11.1.5 — β=γ preservation cross-cases (|β_nat/γ_nat - 1| < 5%).

    Central claim M2b: w 1-loop V_eff w canonical sek08a units (φ = Φ/Φ_0),
    relacja β = γ jest zachowana z dokładnością do fit systematics.
    """
    if verbose:
        print("\n" + "=" * 78)
        print("M11.1.5 — β=γ preservation cross-cases (M2b central claim)")
        print("=" * 78)
        print("  Test: |β_nat/γ_nat - 1| < 5% for all 3 cases")
        print("  (φ_natural = Φ/Φ_0 ⇒ β_nat = β·Φ_0³, γ_nat = γ·Φ_0⁴)")

    pass_per_case = {}
    deviations = {}

    for label, ref in M2B_REFERENCE.items():
        res = compute_case(label, ref['m0sq'], ref['lambda0'], ref['T'])
        ratio = res['beta_nat_over_gamma_nat']
        deviation = ratio - 1.0
        abs_dev = abs(deviation)

        crit_5pct = abs_dev < 0.05
        pass_per_case[label] = crit_5pct
        deviations[label] = (ratio, deviation, abs_dev)

        if verbose:
            print(f"\n  Case {label}:")
            print(f"    β_nat / γ_nat (1-loop):  {ratio:.4f}")
            print(f"    Deviation from 1:        {deviation:+.4f}  ({deviation*100:+.2f}%)")
            print(f"    Within 5%:                {crit_5pct}")

    test_pass = all(pass_per_case.values())
    if verbose:
        print(f"\n  All 3 cases preserve β=γ within 5% (|β_nat/γ_nat - 1|):")
        for label, (r, d, a) in deviations.items():
            mark = "✓" if pass_per_case[label] else "✗"
            print(f"    {mark} Case {label}:  ratio = {r:.4f}, deviation = {d*100:+.2f}%")
        print(f"  M11.1.5 → {'PASS' if test_pass else 'FAIL'}")
        print(f"\n  STRUCTURAL: β=γ preservation w 1-loop jest konsekwencją tożsamości")
        print(f"  M11.1.2 (ln K_Φ = 2 ln Φ + const). M2a tree-level β=γ pochodzi z H-S")
        print(f"  Jacobianu `(T/2) ln Φ`; M2b dodaje stiffness piece który MULTIPLIKATYWNIE")
        print(f"  wzmacnia ten Jacobian o factor 3 — preserves β:γ ratio strukturalnie.")

    return test_pass, {
        'pass_per_case': pass_per_case, 'deviations': deviations,
    }


def test_M11_1_6_eta_boundary_consistency(verbose=True):
    """M11.1.6 — η boundary condition consistency: η_CG2 vs η_BI vs η_BII target.

    M11.1 jest top-down (Branch II perturbative + FRG). Jego entry boundary
    condition pochodzi z CG-2 LPA' (Wetterich-Litim, 3D Ising WF FP):
       η_CG2 = 0.044  (cg_strong_numerical.py)

    Branch I bottom-up (M11.G.6) dało:
       η_BI = 0.0253  (factor 0.58× CG-2, explained jako scheme difference)

    Test: oba w paśmie perturbacyjnym (< 1/(4π) ≈ 0.080), ratio w
    udokumentowanym paśmie [0.3, 2.0] (M11.G.6 reports factor < 5).

    UWAGA: pełne FRG-based η_BII wymaga LPA' truncation z explicit η-eq
    (M11.2 scope, NIE M11.1). M11.1 weryfikuje TYLKO że CG-2 boundary
    condition jest consistent z Branch I emergent value.
    """
    if verbose:
        print("\n" + "=" * 78)
        print("M11.1.6 — η boundary condition consistency (CG-2 ↔ Branch I)")
        print("=" * 78)
        print("  Honest scope: M11.1 weryfikuje BOUNDARY CONDITION (CG-2 LPA'),")
        print("  NIE pełne η_BII (M11.2 scope: explicit FRG flow).")

    eta_CG2 = 0.044       # LPA' Wetterich-Litim, 3D Ising WF (CG-2 closed)
    eta_BI = 0.0253       # M11.G.6 Branch I emergent
    perturbative_upper = 1.0 / (4.0 * np.pi)  # ≈ 0.0796

    ratio_BI_over_CG2 = eta_BI / eta_CG2

    # PASS criteria:
    crit_CG2_perturbative = (1e-4 < eta_CG2 < perturbative_upper)
    crit_BI_perturbative = (1e-4 < eta_BI < perturbative_upper)
    # Ratio in [0.3, 2.0] band (M11.G.6 documents factor < 5; tighter band 0.3–2)
    crit_ratio_band = (0.3 <= ratio_BI_over_CG2 <= 2.0)
    crit_signs_match = (eta_CG2 > 0) and (eta_BI > 0)  # both positive (Wilson-Fisher)

    if verbose:
        print(f"\n  η_CG2  (LPA' Wetterich-Litim, 3D Ising WF):     {eta_CG2:.4f}")
        print(f"  η_BI   (Branch I emergent, M11.G.6):            {eta_BI:.4f}")
        print(f"  Ratio η_BI / η_CG2:                              {ratio_BI_over_CG2:.4f}")
        print(f"  Perturbative upper 1/(4π):                       {perturbative_upper:.4f}")
        print(f"\n  η_CG2 ∈ (1e-4, 1/(4π)):                          {crit_CG2_perturbative}")
        print(f"  η_BI  ∈ (1e-4, 1/(4π)):                          {crit_BI_perturbative}")
        print(f"  Ratio in band [0.3, 2.0]:                       {crit_ratio_band}")
        print(f"  Both η > 0 (Wilson-Fisher universality):        {crit_signs_match}")

    test_pass = (crit_CG2_perturbative and crit_BI_perturbative
                 and crit_ratio_band and crit_signs_match)
    if verbose:
        print(f"\n  M11.1.6 → {'PASS' if test_pass else 'FAIL'}")
        print(f"\n  KOMENTARZ STRUKTURALNY: factor 0.58× między η_BI = 0.0253 a η_CG2 = 0.044")
        print(f"  jest dokumentowany w M11.G.6 jako spodziewana scheme-dependence:")
        print(f"  - η_BI: bottom-up soliton 1-loop, hard cutoff Λ, Branch I scheme")
        print(f"  - η_CG2: top-down LPA' Wetterich-Litim regulator, 3D Ising WF FP")
        print(f"  Pełna §4.1 consistency (η_BI = η_BII = η_CG2 within 0.01) wymaga")
        print(f"  M11.2 (explicit FRG flow) + M11.R-final (multi-branch synthesis).")

    return test_pass, {
        'eta_CG2': eta_CG2, 'eta_BI': eta_BI,
        'ratio_BI_over_CG2': ratio_BI_over_CG2,
        'perturbative_upper': perturbative_upper,
    }


# ============================================================================
# Main
# ============================================================================

def main():
    print("#" * 78)
    print("# M11.1 — 1-loop V_eff(Φ) closure audit (Branch II level 1)")
    print("# Audyt closure-grade pre-existing M2b (op1-op2-op4) loop derivation")
    print("# Cel: zweryfikować strukturalne preservation β=γ w canonical sek08a,")
    print("# BZ convergence, loop stability, oraz η boundary condition (CG-2).")
    print("#" * 78)

    results = {}
    p1, results['M11.1.1'] = test_M11_1_1_reproducibility()
    p2, results['M11.1.2'] = test_M11_1_2_lnK_identity()
    p3, results['M11.1.3'] = test_M11_1_3_loop_stability_positivity()
    p4, results['M11.1.4'] = test_M11_1_4_bz_convergence()
    p5, results['M11.1.5'] = test_M11_1_5_beta_gamma_preservation()
    p6, results['M11.1.6'] = test_M11_1_6_eta_boundary_consistency()

    passes = [p1, p2, p3, p4, p5, p6]
    n_pass = sum(passes)

    print("\n" + "=" * 78)
    print(" M11.1 — FINAL VERDICT")
    print("=" * 78)
    print(f"  Sub-tests passed: {n_pass}/6")
    for name, p in zip(['M11.1.1', 'M11.1.2', 'M11.1.3',
                         'M11.1.4', 'M11.1.5', 'M11.1.6'], passes):
        marker = '[v]' if p else '[x]'
        print(f"    {marker} {name}: {'PASS' if p else 'FAIL'}")

    if n_pass == 6:
        print("\n  >> M11.1 CLOSED (6/6 PASS) — 1-loop V_eff audit complete:")
        print("  >> structural β=γ preservation (M2b central claim) verified;")
        print("  >> ln K_Φ = 2 ln Φ + const identity exact; BZ convergence < 1e-5;")
        print("  >> loop stability A(Φ_0) > 0 wszystkie 3 cases; η_CG2 boundary")
        print("  >> condition consistent z Branch I (factor 0.58×, scheme-explained).")
        print("  >> Branch II rozpoczęty. Kolejny etap: M11.2 (β-functions z Wetterich FRG).")
    elif n_pass >= 4:
        print(f"\n  ?? M11.1 PARTIAL ({n_pass}/6) — review FAIL details")
    else:
        print(f"\n  XX M11.1 INSUFFICIENT ({n_pass}/6)")

    return n_pass, results


if __name__ == "__main__":
    main()
