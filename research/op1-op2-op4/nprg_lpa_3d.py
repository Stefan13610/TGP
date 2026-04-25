"""
M8 — NPRG (Wetterich) cross-check at the 3D Ising WF FP.

Polynomial-truncation LPA fixed-point solver with Litim regulator,
d=3, Z_2 symmetric scalar.

Setup (see M8_NPRG_setup.md):

    Flow eq (LPA, η=0, d=3, Litim):
        ∂_t v(ρ̃) = -3 v(ρ̃) + ρ̃ v'(ρ̃) + c · 1/(1 + m²(ρ̃))
        c = 1/(6π²)
        m²(ρ̃) = v'(ρ̃) + 2 ρ̃ v''(ρ̃)

    FP equation:
        3 v_*(ρ̃) - ρ̃ v_*'(ρ̃) = c · 1/(1 + m_*²(ρ̃))

    Polynomial truncation v_*(ρ̃) = Σ_{k=1}^N a_k ρ̃^k yields
    N coupled equations:
        (3-k) a_k = c · [coeff of ρ̃^k in 1/(1+m_*²)],   k=1..N

    Mapping to MK-RG conventions (V(ŝ) = Σ c_{2k}/(2k) ŝ^{2k}, ŝ²=2ρ̃):
        r = a_1, u = a_2, B = (3/4) a_3, Γ = a_4 / 2
        B/Γ = (3/2) (a_3 / a_4)

WF identification: at the Wilson-Fisher FP in d=3, the linearised
flow has EXACTLY ONE positive (relevant) eigenvalue y_t = 1/ν,
with ν ≈ 0.6496 in LPA (Litim 2001). Other FPs (Gaussian,
multicritical, regulator-artefact) have different eigenvalue
spectra and are rejected.

Strategy: at each truncation N, scan initial guesses in a grid,
run fsolve on each, accept only WF-validated solutions
(nu ∈ [0.55, 0.70], exactly one positive eigenvalue).

Goal: compute B*/Γ* and compare with MK-RG result B*/Γ* = -0.5687
(M3 at N_ops=8). If both negative ⇒ OP-2b confirmed open at the
universality-class level (M3-M7 verdict universal, not MK-specific).

Usage:
    python nprg_lpa_3d.py
"""

import itertools

import numpy as np
from scipy.optimize import fsolve

# Loop integral prefactor for d=3, Litim regulator, LPA (η=0).
# I_d = 1/(6π²) for d=3.
C_LITIM_3D = 1.0 / (6.0 * np.pi ** 2)

# WF validation gates
NU_LO, NU_HI = 0.55, 0.70   # acceptable LPA ν range (literature N→∞: 0.6496)
N_POS_EIG = 1               # WF FP: exactly one positive eigenvalue


def m2_coeffs(a, N):
    """Coefficients of m²(ρ̃) = v'(ρ̃) + 2ρ̃ v''(ρ̃) as polynomial in ρ̃.

    Given v(ρ̃) = Σ_{k=1}^N a_k ρ̃^k, formula:
        [m²]_j = (j+1)(2j+1) a_{j+1},   j = 0, 1, ..., N-1.

    Returns array of length N.
    """
    return np.array([(j + 1) * (2 * j + 1) * a[j] for j in range(N)])


def invert_polynomial(A, N):
    """Formal series inversion of polynomial A(ρ̃) up to order N.

    Given A(ρ̃) = A[0] + A[1] ρ̃ + ... + A[m] ρ̃^m (m = len(A)-1),
    returns B such that A·B = 1 + O(ρ̃^{N+1}).

    Output B has length N+1.
    """
    if A[0] == 0:
        raise ValueError("Cannot invert: A[0] = 0")
    B = np.zeros(N + 1)
    B[0] = 1.0 / A[0]
    m = len(A) - 1
    for j in range(1, N + 1):
        s = 0.0
        for i in range(1, min(j, m) + 1):
            s += A[i] * B[j - i]
        B[j] = -s / A[0]
    return B


def fp_residuals(a, N, c=C_LITIM_3D):
    """Residuals F_k = (3-k) a_k - c · B[k] for k = 1..N."""
    m2 = m2_coeffs(a, N)
    A = np.zeros(N)
    A[0] = 1.0 + m2[0]
    if N >= 2:
        A[1:N] = m2[1:N]
    B = invert_polynomial(A, N)
    F = np.empty(N)
    for k in range(1, N + 1):
        F[k - 1] = (3 - k) * a[k - 1] - c * B[k]
    return F


def evaluate_v_derivs(a, rho, k_max=4):
    """Evaluate v_*(rho), v_*'(rho), ..., up to v_*^(k_max).
    Returns array of length k_max+1."""
    N = len(a)
    derivs = np.zeros(k_max + 1)
    # v(rho) = sum_{k=1}^N a_k rho^k
    # v^(d)(rho) = sum_{k=d}^N a_k * (k!/(k-d)!) * rho^{k-d}
    for d in range(k_max + 1):
        s = 0.0
        for k in range(max(d, 1), N + 1):
            coef = 1.0
            for j in range(d):
                coef *= (k - j)
            s += a[k - 1] * coef * (rho ** (k - d) if k > d else 1.0)
        derivs[d] = s
    return derivs


def find_rho_min(a, rho_init=0.03):
    """Solve v_*'(rho_0) = 0 by Newton iteration starting from rho_init."""
    from scipy.optimize import brentq
    # Bracket: v'(0) = a[0] < 0, v'(rho_far) > 0 if FP function is well-behaved
    f = lambda r: evaluate_v_derivs(a, r, k_max=1)[1]
    if f(0.0) >= 0:
        return None  # No minimum at positive rho
    # Find a positive value where f > 0
    rho_far = 1e-3
    for _ in range(50):
        if f(rho_far) > 0:
            break
        rho_far *= 1.4
    else:
        return None
    try:
        rho_min = brentq(f, 0.0, rho_far, xtol=1e-12, rtol=1e-10)
    except (ValueError, RuntimeError):
        return None
    return rho_min


def beta_gamma_at_fp_min(a):
    """Compute the v2 source-shifted (β, γ) at the FP minimum ρ̃_0.

    Setup: Φ_0² = 2ρ̃_0; expand v(Φ) = v(Φ_0 + φ) around φ = 0:
        v(Φ_0 + φ) = const + (1/2) m² φ² + (1/6) v'''(Φ_0) φ³
                            + (1/24) v''''(Φ_0) φ⁴ + ...
    Identification with U(φ) = (β/3) φ³ - (γ/4) φ⁴:
        β = v'''(Φ_0) / 2,    γ = -v''''(Φ_0) / 6
        β/γ = -3 v'''(Φ_0) / v''''(Φ_0)

    Compute v''', v'''' from v_*(ρ̃) via chain rule:
        Let w(ρ̃) = v_*(ρ̃) and Φ² = 2ρ̃. Then
        v_Φ           = Φ w'
        v_ΦΦ          = w' + 2ρ̃ w''
        v_ΦΦΦ         = Φ (3 w'' + 2ρ̃ w''')
        v_ΦΦΦΦ        = 3 w'' + 12 ρ̃ w''' + 4 ρ̃² w''''.
    At ρ̃ = ρ̃_0 (where w'(ρ̃_0) = 0):
        v_ΦΦ(Φ_0)   = 2 ρ̃_0 w''(ρ̃_0) = m²(ρ̃_0)
        v_ΦΦΦ(Φ_0)  = Φ_0 (3 w''(ρ̃_0) + 2 ρ̃_0 w'''(ρ̃_0))
        v_ΦΦΦΦ(Φ_0) = 3 w''(ρ̃_0) + 12 ρ̃_0 w'''(ρ̃_0) + 4 ρ̃_0² w''''(ρ̃_0).

    Returns dict with rho_min, Phi_0, beta, gamma, beta_over_gamma, m_sq.
    """
    rho_0 = find_rho_min(a)
    if rho_0 is None:
        return None
    derivs = evaluate_v_derivs(a, rho_0, k_max=4)
    w_p, w_pp, w_ppp, w_pppp = derivs[1], derivs[2], derivs[3], derivs[4]
    # w_p ≈ 0 by construction
    Phi_0 = np.sqrt(2.0 * rho_0)
    v_PhiPhi = w_p + 2.0 * rho_0 * w_pp        # ≈ 2 rho_0 w''
    v_PhiPhiPhi = Phi_0 * (3.0 * w_pp + 2.0 * rho_0 * w_ppp)
    v_PhiPhiPhiPhi = 3.0 * w_pp + 12.0 * rho_0 * w_ppp + 4.0 * rho_0 ** 2 * w_pppp
    beta = v_PhiPhiPhi / 2.0
    gamma = -v_PhiPhiPhiPhi / 6.0
    return {
        "rho_min": rho_0,
        "Phi_0": Phi_0,
        "v_sq": 2.0 * rho_0,
        "m_sq": v_PhiPhi,
        "v_PhiPhiPhi": v_PhiPhiPhi,
        "v_PhiPhiPhiPhi": v_PhiPhiPhiPhi,
        "beta": beta,
        "gamma": gamma,
        "beta_over_gamma": beta / gamma if gamma != 0 else None,
    }


def jacobian_numerical(a, N, eps=1e-8):
    """Numerical Jacobian of fp_residuals w.r.t. a, central differences."""
    J = np.empty((N, N))
    for j in range(N):
        a_p = a.copy()
        a_m = a.copy()
        a_p[j] += eps
        a_m[j] -= eps
        F_p = fp_residuals(a_p, N)
        F_m = fp_residuals(a_m, N)
        J[:, j] = (F_p - F_m) / (2 * eps)
    return J


def is_wf_solution(a_star, N, tol_res=1e-6):
    """Check whether (a_star) satisfies WF criteria:
        - low residual
        - exactly N_POS_EIG positive eigenvalues
        - ν ∈ [NU_LO, NU_HI]
    Returns (is_wf, info_dict).
    """
    res = fp_residuals(a_star, N)
    res_max = np.max(np.abs(res))
    info = {"residual_max": res_max, "rejected": None}
    if res_max > tol_res:
        info["rejected"] = f"residual {res_max:.2e}"
        return False, info
    if abs(a_star[1]) < 1e-3:
        info["rejected"] = "trivial Gaussian (a_2 ≈ 0)"
        return False, info
    J = jacobian_numerical(a_star, N)
    eigvals = np.linalg.eigvals(J)
    eigvals_real = np.sort(eigvals.real)[::-1]
    info["eigvals_real"] = eigvals_real
    pos_eigs = eigvals_real[eigvals_real > 1e-6]
    info["n_pos_eigs"] = len(pos_eigs)
    if len(pos_eigs) != N_POS_EIG:
        info["rejected"] = f"n_pos_eigs={len(pos_eigs)} ≠ {N_POS_EIG}"
        return False, info
    nu = 1.0 / pos_eigs[0]
    info["nu_lpa"] = nu
    if not (NU_LO <= nu <= NU_HI):
        info["rejected"] = f"ν={nu:.4f} outside [{NU_LO}, {NU_HI}]"
        return False, info
    return True, info


def solve_fp_single(N, a_init, tol_res=1e-6):
    """Run fsolve from a single initial guess. Return (a_star, info) or None."""
    if len(a_init) < N:
        a_init = np.concatenate([a_init, np.zeros(N - len(a_init))])
    elif len(a_init) > N:
        a_init = a_init[:N]
    try:
        a_star, _, ier, _ = fsolve(fp_residuals, a_init, args=(N,),
                                    full_output=True, xtol=1e-13,
                                    maxfev=20000)
    except (FloatingPointError, ValueError, ZeroDivisionError):
        return None
    if ier != 1:
        return None
    is_wf, info = is_wf_solution(a_star, N, tol_res=tol_res)
    if is_wf:
        info["a_star"] = a_star
        return info
    return None


def make_initial_grid(N):
    """Build a grid of initial guesses for the WF FP at truncation N.

    Strategy: alternating-sign, MK-RG-inspired magnitudes, plus the
    N=2 analytic result and a Gaussian-like seed.
    """
    guesses = []

    # MK-RG-inspired (alternating signs, increasing magnitudes — matches
    # the M3 result rescaled to NPRG conventions: a_3 < 0, a_4 > 0, etc.)
    for r0 in (-0.03, -0.05, -0.08, -0.12, -0.20):
        for u0 in (0.5, 1.0, 1.5, 2.0):
            for s_a3 in (-1.0, -3.0, -5.0):  # a_3 negative
                for s_a4 in (1.0, 5.0, 15.0, 30.0):  # a_4 positive
                    g = np.zeros(N)
                    g[0] = r0
                    g[1] = u0
                    if N >= 3:
                        g[2] = s_a3
                    if N >= 4:
                        g[3] = s_a4
                    # Higher: alternating signs, decaying
                    for k in range(4, N):
                        sign = (-1) ** (k - 2)
                        g[k] = sign * (5.0 ** (k - 3)) * 1.5
                    guesses.append(g)

    # N=2 analytic seed (no higher coeffs)
    g0 = np.zeros(N)
    g0[0] = -1.0 / 13.0
    g0[1] = (12.0 / 13.0) ** 3 * np.pi ** 2 / 6.0
    guesses.append(g0)

    # Random fallback
    rng = np.random.default_rng(seed=42)
    for _ in range(20):
        g = np.zeros(N)
        g[0] = -rng.uniform(0.01, 0.30)
        g[1] = rng.uniform(0.3, 3.0)
        for k in range(2, N):
            sign = (-1) ** k
            g[k] = sign * rng.uniform(0.5, 50.0) * (1.5 ** (k - 2))
        guesses.append(g)

    return guesses


def solve_fp(N, warm_start=None, verbose=False):
    """Robust WF-FP solver: scan initial guesses, return WF-validated result."""
    candidates = []

    # Warm start first (if provided)
    if warm_start is not None:
        result = solve_fp_single(N, warm_start)
        if result is not None:
            candidates.append(result)

    # Then scan grid
    grid = make_initial_grid(N)
    for g in grid:
        result = solve_fp_single(N, g)
        if result is not None:
            candidates.append(result)

    if not candidates:
        if verbose:
            print(f"  N={N}: no WF candidate found")
        return None

    # Pick smallest residual; among tied, pick ν closest to literature 0.6496
    candidates.sort(key=lambda c: (c["residual_max"], abs(c["nu_lpa"] - 0.6496)))
    best = candidates[0]

    # Cluster-deduplicate (count distinct WF candidates by rounded a_3/a_4)
    distinct_a3a4 = set()
    for c in candidates:
        a = c["a_star"]
        if N >= 4:
            distinct_a3a4.add(round(a[2] / a[3], 4))
    best["distinct_a3a4_count"] = len(distinct_a3a4)
    best["n_candidates"] = len(candidates)

    # Mapping to MK-RG observables
    a = best["a_star"]
    if N >= 4:
        best["a3_over_a4"] = a[2] / a[3]
        best["B_over_Gamma"] = 1.5 * (a[2] / a[3])
    else:
        best["a3_over_a4"] = None
        best["B_over_Gamma"] = None
    best["v_sq_proxy"] = abs(a[0]) / a[1]
    best["N"] = N

    if verbose:
        print(f"  N={N}: WF found from {best['n_candidates']} candidates, "
              f"{best['distinct_a3a4_count']} distinct a_3/a_4 values, "
              f"ν={best['nu_lpa']:.4f}")
    return best


def analytic_check_N2():
    """Closed-form N=2 truncation values.

        a_1* = -1/13,   a_2* = (12/13)^3 · π²/6,   ν(N=2) ≈ 0.5429
    """
    a1 = -1.0 / 13.0
    a2 = (12.0 / 13.0) ** 3 * np.pi ** 2 / 6.0
    return a1, a2


def main():
    print("=" * 72)
    print("M8 — NPRG-LPA at 3D Ising WF FP, polynomial truncation around ρ̃=0")
    print(f"     Litim regulator, d=3, Z2.  c_loop = 1/(6π²) = {C_LITIM_3D:.6f}")
    print(f"     WF validation: ν ∈ [{NU_LO}, {NU_HI}], n_pos_eigs = {N_POS_EIG}")
    print("=" * 72)

    # ---- Analytic N=2 sanity check ----
    a1_exact, a2_exact = analytic_check_N2()
    print(f"\n[Analytic N=2] a_1* = {a1_exact:.8f}    a_2* = {a2_exact:.8f}")

    # N=2 has only one WF FP and ν=0.543 < NU_LO=0.55, so we accept it as
    # a special case (literature confirms N=2 LPA ν ≈ 0.543).
    res2_a = np.array([a1_exact, a2_exact])
    res2_info = is_wf_solution(res2_a, 2)
    # Force-acceptance for N=2 baseline reporting:
    print(f"[N=2 analytic] residual_max = {is_wf_solution(res2_a, 2)[1]['residual_max']:.2e}")
    eigvals_N2 = np.linalg.eigvals(jacobian_numerical(res2_a, 2))
    eigs_real = np.sort(eigvals_N2.real)[::-1]
    nu_N2 = 1.0 / eigs_real[0] if eigs_real[0] > 0 else None
    print(f"[N=2 analytic] eigvals = {eigs_real}, ν = {nu_N2:.6f}  "
          "(literature N=2: 0.543; gate widens for higher N)")

    # ---- Scan higher N ----
    print("\n" + "=" * 72)
    print("Robust WF-FP scan over truncation order N")
    print("=" * 72)
    print(f"{'N':>3}  {'a_1*':>10}  {'a_2*':>10}  {'a_3*':>10}  {'a_4*':>10}  "
          f"{'a_3/a_4':>10}  {'B*/Γ*':>10}  {'ν_LPA':>8}  {'cands':>5}")
    print("-" * 100)

    results = {}
    a_warm = None
    for N in [4, 5, 6, 7, 8, 10, 12, 14, 16]:
        res = solve_fp(N, warm_start=a_warm, verbose=True)
        if res is None:
            print(f"{N:>3}  --- WF FP not found at this truncation ---")
            continue
        results[N] = res
        a_warm = res["a_star"]

        a = res["a_star"]
        ag = lambda i: f"{a[i]:>+10.6f}" if i < N else " " * 10
        a3a4_s = (f"{res['a3_over_a4']:>+10.6f}"
                  if res["a3_over_a4"] is not None else " " * 10)
        BoG_s = (f"{res['B_over_Gamma']:>+10.6f}"
                 if res["B_over_Gamma"] is not None else " " * 10)
        nu_s = f"{res['nu_lpa']:>8.5f}" if res["nu_lpa"] is not None else " " * 8
        print(f"{N:>3}  {ag(0)}  {ag(1)}  {ag(2)}  {ag(3)}  "
              f"{a3a4_s}  {BoG_s}  {nu_s}  {res['n_candidates']:>5}")

        # Compute scheme-independent (β, γ) at FP minimum
        bg = beta_gamma_at_fp_min(a)
        if bg is not None:
            res["fp_min"] = bg

    # ---- Detailed per-N report ----
    print("\n" + "=" * 72)
    print("Detailed FP coefficients (selected truncations)")
    print("=" * 72)
    for N in sorted(results):
        if N not in (4, 6, 8, 10, 12, 14, 16):
            continue
        r = results[N]
        print(f"\n--- N = {N} ---")
        for i, ai in enumerate(r["a_star"], start=1):
            print(f"  a_{i:>2} = {ai:+.10f}")
        print(f"  ν_LPA       = {r['nu_lpa']:.6f}")
        print(f"  a_3/a_4     = {r['a3_over_a4']:+.6f}")
        print(f"  B*/Γ*       = (3/2)(a_3/a_4) = {r['B_over_Gamma']:+.6f}")
        print(f"  |a_1|/a_2  = {r['v_sq_proxy']:.6f}  "
              "(analogue of v*² in MK-RG)")
        print(f"  residual    = {r['residual_max']:.2e}")
        evs = r["eigvals_real"][:min(6, N)]
        ev_str = ", ".join(f"{e:+.4f}" for e in evs)
        print(f"  eigenvalues (top {min(6, N)}): [{ev_str}]")
        print(f"  candidates  = {r['n_candidates']}, "
              f"distinct a_3/a_4 = {r['distinct_a3a4_count']}")

    # ---- Scheme-independent β/γ at FP minimum ----
    print("\n" + "=" * 72)
    print("Scheme-independent β/γ at the FP minimum ρ̃_0")
    print("=" * 72)
    print("(NPRG: from v_*(ρ̃) Taylor expansion, derivatives at ρ̃_0)")
    print(f"{'N':>3}  {'ρ̃_0':>10}  {'Φ_0':>10}  {'m²(ρ̃_0)':>10}  "
          f"{'β':>12}  {'γ':>12}  {'β/γ':>10}")
    print("-" * 90)
    for N in sorted(results):
        r = results[N]
        if "fp_min" not in r:
            continue
        b = r["fp_min"]
        bog = b["beta_over_gamma"]
        bog_s = f"{bog:>+10.5f}" if bog is not None else " " * 10
        print(f"{N:>3}  {b['rho_min']:>10.6f}  {b['Phi_0']:>10.6f}  "
              f"{b['m_sq']:>10.5f}  {b['beta']:>+12.5f}  {b['gamma']:>+12.5f}  "
              f"{bog_s}")

    # ---- Same observable for MK-RG FP (as a sanity-cross-check) ----
    print("\n--- Same observable evaluated at the MK-RG FP polynomial ---")
    # MK FP at N_ops=8: r* = -2.45694, u* = +3.01611, B* = -5.12163, Γ* = +9.00600
    # → NPRG variables: a = (r, u, (4/3)B, 2Γ, ...) up to a_4
    a_MK = np.array([-2.45694, 3.01611, (4.0 / 3.0) * (-5.12163),
                     2.0 * 9.00600])
    bg_MK = beta_gamma_at_fp_min(a_MK)
    if bg_MK is not None:
        print(f"  MK FP (a = {a_MK}):")
        print(f"    ρ̃_0 = {bg_MK['rho_min']:.6f}, Φ_0 = {bg_MK['Phi_0']:.6f}")
        print(f"    m² = {bg_MK['m_sq']:.5f}")
        print(f"    β  = {bg_MK['beta']:+.5f}")
        print(f"    γ  = {bg_MK['gamma']:+.5f}")
        print(f"    β/γ = {bg_MK['beta_over_gamma']:+.5f}")
    else:
        print("  MK FP: no FP minimum found in (a_1, a_2, a_3, a_4) "
              "polynomial truncation.")

    # ---- Comparison with M3 (polynomial-coefficient ratio) ----
    print("\n" + "=" * 72)
    print("Comparison with MK-RG (polynomial-coefficient ratio B*/Γ*)")
    print("=" * 72)
    print("Note: B*/Γ* is a SCHEME-DEPENDENT polynomial ratio.")
    print("M3 (MK, N_ops=8): B*/Γ* = -0.5687  ⟺  a_3/a_4 = -0.3791")
    print()
    for N in sorted(results):
        r = results[N]
        if r["B_over_Gamma"] is None:
            continue
        BG = r["B_over_Gamma"]
        a34 = r["a3_over_a4"]
        delta = BG - (-0.5687)
        sign_match = "MATCH" if BG < 0 else "OPPOSITE"
        print(f"  N={N:>3}: B*/Γ* = {BG:+8.4f}   a_3/a_4 = {a34:+8.4f}   "
              f"Δ = {delta:+7.4f}  [sign: {sign_match}]")

    # ---- Verdict ----
    print("\n" + "=" * 72)
    print("Verdict")
    print("=" * 72)
    if not results:
        print("→ NO WF FP found at any truncation. Polynomial-around-ρ̃=0 may")
        print("  not converge for this regulator. Need broken-phase scheme or")
        print("  ODE shooting.")
        return
    best_N = max(N for N in results if results[N]["B_over_Gamma"] is not None)
    BG = results[best_N]["B_over_Gamma"]
    a34 = results[best_N]["a3_over_a4"]
    nu_best = results[best_N]["nu_lpa"]
    fp_min_best = results[best_N].get("fp_min", None)
    print(f"At highest converged N = {best_N}:")
    print(f"   ν_LPA    = {nu_best:.6f}  (literature ν_LPA = 0.6496)")
    print(f"   Polynomial-ratio observable B*/Γ* = (3/2)(a_3/a_4) = {BG:+.6f}")
    print(f"   MK-RG    : B*/Γ* = -0.5687  (scheme-dependent)")
    if fp_min_best is not None:
        bog = fp_min_best["beta_over_gamma"]
        print(f"   Scheme-independent β/γ at ρ̃_0 = {fp_min_best['rho_min']:.5f}: "
              f"{bog:+.5f}")
    print()
    print("-- Sign analysis --")
    if BG < -0.2:
        print("  Polynomial ratio: NPRG agrees with MK-RG (negative).")
    elif BG > 0.2:
        print("  Polynomial ratio: NPRG OPPOSITE SIGN to MK-RG.")
        print("  ⟹ MK-RG's B*/Γ* < 0 is a SCHEME ARTEFACT, not universal.")
    else:
        print("  Polynomial ratio: near zero (inconclusive).")
    print()
    if fp_min_best is not None and fp_min_best["beta_over_gamma"] is not None:
        bog = fp_min_best["beta_over_gamma"]
        print("-- Scheme-independent verdict for OP-2b (β/γ at FP minimum) --")
        if abs(bog - 1.0) < 0.1:
            print(f"  β/γ ≈ {bog:.3f}: OP-2b CLOSED at NPRG-LPA (β = γ within 10%).")
        elif 0.5 <= bog <= 1.5:
            print(f"  β/γ ≈ {bog:.3f}: OP-2b PARTIAL CLOSURE")
            print("  (right sign, magnitude within factor 2 of v1 prediction β=γ).")
        elif 0 < bog < 0.5 or bog > 1.5:
            print(f"  β/γ ≈ {bog:.3f}: WRONG MAGNITUDE but RIGHT SIGN.")
            print("  Higher truncation / LPA' needed for sharper statement.")
        elif bog < 0:
            print(f"  β/γ ≈ {bog:.3f}: WRONG SIGN. OP-2b open at NPRG level.")
        else:
            print(f"  β/γ ≈ {bog:.3f}: degenerate case (γ ≈ 0).")
    print()


if __name__ == "__main__":
    main()
