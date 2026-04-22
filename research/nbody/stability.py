"""
stability.py  --  TGP n-body: stability analysis
==================================================
Compute the Hessian matrix of the total potential energy at equilibrium
configurations and classify stability via eigenvalue analysis.

METHOD:
  1. Compute the 3n x 3n Hessian H = d^2V / dx_a dx_b numerically.
  2. Construct an orthonormal basis Z for the zero-mode subspace
     (3 translations + up to 3 rotations).
  3. Project: H_phys = P^T H P, where P = I - Z Z^T.
  4. Eigenvalues of H_phys on the physical subspace determine stability.

  This is rigorous: zero modes are explicitly constructed from the
  geometry (not identified by a numerical threshold).

STATUS:
  - compute_hessian_numerical(): NUMERICAL (finite differences of V)
  - compute_hessian_generic(): accepts any V(pos) callable (P2 extension)
  - normal_mode_analysis(): mass-weighted Hessian -> omega^2 (P2 extension)
  - stability_comparison(): TGP vs Newton at same config (P2 extension)
  - stability_bifurcation_scan(): (beta,C) parameter scan (P2 extension)
  - classify_stability():
      with positions: RIGOROUS (explicit zero-mode projection)
      without positions: HEURISTIC (threshold-based, legacy mode)
  - earnshaw_test(): uses PAIRWISE energy only (no 3-body in Hessian)

KEY PHYSICS:
In Newtonian gravity, Earnshaw's theorem states that no stable static
equilibrium exists for point masses in a 1/r^2 force law.
TGP breaks this because the repulsive 1/d^2 term creates potential wells.

IMPORTANT: stability results from compute_hessian_numerical() depend on
which energy function is passed. By default, PAIRWISE ONLY energy is used.
For full 2b+3b stability, pass full_potential_tgp from dynamics_v2.py.
"""

import numpy as np


def compute_hessian_numerical(positions, C_values, beta, gamma=None,
                              energy_func=None, dx=1e-5):
    """
    Compute the 3n × 3n Hessian matrix of V_total by finite differences.

    H_{αβ} = ∂²V / ∂x_α ∂x_β

    where α, β run over all 3n Cartesian coordinates of n bodies.

    Parameters
    ----------
    positions : (n, 3) ndarray
    C_values : (n,) ndarray
    beta : float
    gamma : float or None (defaults to beta)
    energy_func : callable
        V(positions, C_values, beta, gamma) → float
    dx : float
        Finite difference step.

    Returns
    -------
    H : (3n, 3n) ndarray
        Hessian matrix.
    """
    if gamma is None:
        gamma = beta

    if energy_func is None:
        # Default pairwise energy
        def energy_func(pos, Cv, b, g):
            V = 0.0
            for i in range(len(Cv)):
                for j in range(i+1, len(Cv)):
                    d = np.linalg.norm(pos[i] - pos[j])
                    d = max(d, 1e-8)
                    Ci, Cj = Cv[i], Cv[j]
                    V += -4*np.pi*Ci*Cj/d \
                         + 8*np.pi*b*Ci*Cj/d**2 \
                         - 24*np.pi*g*Ci*Cj*(Ci+Cj)/(2*d**3)
            return V

    n = len(C_values)
    ndof = 3 * n
    pos_flat = positions.flatten().copy()

    H = np.zeros((ndof, ndof))

    V0 = energy_func(positions, C_values, beta, gamma)

    for a in range(ndof):
        for b in range(a, ndof):
            # H_ab = (V(+dx_a,+dx_b) - V(+dx_a,-dx_b) - V(-dx_a,+dx_b) + V(-dx_a,-dx_b)) / (4*dx^2)
            pos_pp = pos_flat.copy()
            pos_pm = pos_flat.copy()
            pos_mp = pos_flat.copy()
            pos_mm = pos_flat.copy()

            pos_pp[a] += dx; pos_pp[b] += dx
            pos_pm[a] += dx; pos_pm[b] -= dx
            pos_mp[a] -= dx; pos_mp[b] += dx
            pos_mm[a] -= dx; pos_mm[b] -= dx

            Vpp = energy_func(pos_pp.reshape(n, 3), C_values, beta, gamma)
            Vpm = energy_func(pos_pm.reshape(n, 3), C_values, beta, gamma)
            Vmp = energy_func(pos_mp.reshape(n, 3), C_values, beta, gamma)
            Vmm = energy_func(pos_mm.reshape(n, 3), C_values, beta, gamma)

            H[a, b] = (Vpp - Vpm - Vmp + Vmm) / (4 * dx**2)
            H[b, a] = H[a, b]

    return H


def _build_zero_mode_basis(positions):
    """
    Construct an orthonormal basis for the null space of translations
    and infinitesimal rotations.

    For n bodies with coordinates (x1,y1,z1, x2,y2,z2, ...):
      - 3 translation modes: uniform shift in x, y, z
      - up to 3 rotation modes: infinitesimal rotations about CoM

    For collinear configurations (all on one axis): only 2 rotations
    are independent (rotation about the collinear axis is trivial).

    Returns
    -------
    Z : (3n, n_zero) ndarray
        Orthonormal columns spanning the zero-mode subspace.
    n_trans : int
    n_rot : int
    is_collinear : bool
    """
    n = len(positions)
    ndof = 3 * n

    modes = []

    # --- Translation modes ---
    for axis in range(3):
        v = np.zeros(ndof)
        for i in range(n):
            v[3 * i + axis] = 1.0
        v /= np.linalg.norm(v)
        modes.append(v)

    # --- Rotation modes ---
    # CoM
    com = np.mean(positions, axis=0)
    rel = positions - com  # relative positions

    # Check collinearity: SVD of relative position matrix
    _, svals, _ = np.linalg.svd(rel, full_matrices=False)
    # Collinear if only 1 singular value is significant
    sig_svals = np.sum(svals > svals[0] * 1e-8)
    is_collinear = (sig_svals <= 1)

    # Infinitesimal rotation about axis e_a:
    # delta_r_i = e_a x (r_i - com)
    axes_rot = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]

    for e_a in axes_rot:
        v = np.zeros(ndof)
        for i in range(n):
            cross = np.cross(e_a, rel[i])
            v[3 * i: 3 * i + 3] = cross
        norm = np.linalg.norm(v)
        if norm > 1e-10:
            v /= norm
            modes.append(v)

    # Orthonormalize via Gram-Schmidt
    ortho = []
    for v in modes:
        for u in ortho:
            v = v - np.dot(v, u) * u
        norm = np.linalg.norm(v)
        if norm > 1e-10:
            ortho.append(v / norm)

    Z = np.column_stack(ortho) if ortho else np.zeros((ndof, 0))
    n_trans = 3
    n_rot = Z.shape[1] - 3

    return Z, n_trans, n_rot, is_collinear


def classify_stability(H, positions=None, project_zero_modes=True):
    """
    Classify stability of an equilibrium from its Hessian.

    If positions are provided and project_zero_modes=True, the Hessian
    is projected onto the physical subspace (orthogonal to translations
    and rotations) before eigenvalue analysis. This is the rigorous
    method — it explicitly constructs the zero-mode basis and projects
    it out, rather than relying on a numerical threshold.

    Parameters
    ----------
    H : (3n, 3n) ndarray
        Hessian matrix at equilibrium.
    positions : (n, 3) ndarray or None
        Equilibrium positions. Required for proper zero-mode projection.
        If None, falls back to threshold-based classification (heuristic).
    project_zero_modes : bool
        If True and positions given, project out translations + rotations.

    Returns
    -------
    dict with:
        eigenvalues_full : all eigenvalues of H (sorted)
        eigenvalues_physical : eigenvalues of projected H_phys (sorted)
        n_positive : positive physical eigenvalues
        n_negative : negative physical eigenvalues
        n_zero_modes : number of projected-out zero modes
        stable : bool (all physical eigenvalues > 0)
        classification : 'stable', 'unstable', 'saddle', 'degenerate'
        frequencies : sqrt(positive eigenvalues) — oscillation frequencies
        is_collinear : bool
        method : 'projection' or 'threshold' (which method was used)
    """
    eigenvalues_full = np.linalg.eigvalsh(H)

    if positions is not None and project_zero_modes:
        # --- Rigorous projection method ---
        n = len(positions)
        Z, n_trans, n_rot, is_collinear = _build_zero_mode_basis(positions)
        n_zero_modes = Z.shape[1]

        # Projector onto physical subspace: P = I - Z @ Z^T
        ndof = 3 * n
        P = np.eye(ndof) - Z @ Z.T

        # Projected Hessian: H_phys = P^T @ H @ P
        H_phys = P.T @ H @ P

        # Eigenvalues of projected Hessian
        evals_proj = np.linalg.eigvalsh(H_phys)

        # The projected-out modes give ~0 eigenvalues; keep only
        # the (3n - n_zero_modes) physical ones
        evals_sorted = np.sort(evals_proj)
        evals_physical = evals_sorted[n_zero_modes:]

        method = 'projection'

    else:
        # --- Fallback: threshold-based (heuristic) ---
        n = H.shape[0] // 3
        max_abs = max(np.abs(eigenvalues_full).max(), 1e-10)
        threshold = max_abs * 1e-6

        evals_physical = eigenvalues_full[np.abs(eigenvalues_full) >= threshold]
        n_zero_modes = int(np.sum(np.abs(eigenvalues_full) < threshold))
        is_collinear = None
        method = 'threshold'

    # Threshold for classifying physical eigenvalues
    if len(evals_physical) > 0:
        phys_max = max(np.abs(evals_physical).max(), 1e-10)
        phys_thresh = phys_max * 1e-8
    else:
        phys_thresh = 1e-10

    n_positive = int(np.sum(evals_physical > phys_thresh))
    n_negative = int(np.sum(evals_physical < -phys_thresh))

    if len(evals_physical) == 0:
        classification = "degenerate"
    elif np.all(evals_physical > -phys_thresh):
        classification = "stable"
    elif np.all(evals_physical < phys_thresh):
        classification = "unstable"
    else:
        classification = "saddle"

    # Oscillation frequencies for stable modes
    pos_evals = evals_physical[evals_physical > phys_thresh]
    frequencies = np.sqrt(pos_evals) if len(pos_evals) > 0 else np.array([])

    return {
        'eigenvalues_full': np.sort(eigenvalues_full),
        'eigenvalues_physical': np.sort(evals_physical),
        'n_positive': n_positive,
        'n_negative': n_negative,
        'n_zero_modes': n_zero_modes,
        'stable': classification == "stable",
        'classification': classification,
        'frequencies': np.sort(frequencies),
        'is_collinear': is_collinear,
        'method': method,
    }


# ── P2 Extensions: Full Hessian Stability Mapping ─────────────────────────


def compute_hessian_generic(positions, potential_fn, dx=1e-5):
    """
    Compute the 3n x 3n Hessian of any potential V(positions).

    Unlike compute_hessian_numerical(), this accepts a GENERIC callable
    with signature potential_fn(positions) -> float, making it usable
    for Newton, TGP pairwise, TGP full (V2+V3), or any other potential.

    Parameters
    ----------
    positions : (n, 3) ndarray
    potential_fn : callable
        V(positions) -> float. Must accept (n, 3) ndarray.
    dx : float
        Finite difference step.

    Returns
    -------
    H : (3n, 3n) ndarray
        Symmetric Hessian matrix.
    """
    n = len(positions)
    ndof = 3 * n
    pos_flat = positions.flatten().copy()

    H = np.zeros((ndof, ndof))

    for a in range(ndof):
        for b in range(a, ndof):
            pos_pp = pos_flat.copy()
            pos_pm = pos_flat.copy()
            pos_mp = pos_flat.copy()
            pos_mm = pos_flat.copy()

            pos_pp[a] += dx; pos_pp[b] += dx
            pos_pm[a] += dx; pos_pm[b] -= dx
            pos_mp[a] -= dx; pos_mp[b] += dx
            pos_mm[a] -= dx; pos_mm[b] -= dx

            Vpp = potential_fn(pos_pp.reshape(n, 3))
            Vpm = potential_fn(pos_pm.reshape(n, 3))
            Vmp = potential_fn(pos_mp.reshape(n, 3))
            Vmm = potential_fn(pos_mm.reshape(n, 3))

            H[a, b] = (Vpp - Vpm - Vmp + Vmm) / (4 * dx**2)
            H[b, a] = H[a, b]

    return H


def normal_mode_analysis(positions, C_values, potential_fn, dx=1e-5):
    """
    Full normal mode analysis: mass-weighted Hessian eigenvalues and eigenvectors.

    Computes the dynamical matrix D = M^{-1/2} H M^{-1/2} where
    M = diag(C_1, C_1, C_1, C_2, ..., C_n, C_n, C_n) is the mass matrix
    (TGP axiom: m_i = C_i).

    Physical frequencies: omega_k = sqrt(lambda_k) for lambda_k > 0.

    Zero modes (translations + rotations) are projected out using
    the rigorous geometric construction from _build_zero_mode_basis().

    Parameters
    ----------
    positions : (n, 3) ndarray
        Equilibrium positions.
    C_values : (n,) ndarray
        Source strengths (= inertial masses).
    potential_fn : callable
        V(positions) -> float.
    dx : float
        Finite difference step for Hessian.

    Returns
    -------
    dict with:
        omega2_physical : (n_phys,) eigenvalues of projected dynamical matrix
        omega_physical : sqrt of positive eigenvalues (oscillation frequencies)
        eigenvectors_physical : (3n, n_phys) mass-weighted eigenvectors
        H : (3n, 3n) raw Hessian
        D : (3n, 3n) dynamical matrix M^{-1/2} H M^{-1/2}
        n_stable : number of positive omega^2 modes
        n_unstable : number of negative omega^2 modes
        n_zero_modes : number of projected-out zero modes
        classification : 'stable', 'unstable', 'saddle', 'degenerate'
        is_collinear : bool
        mode_characters : list of dicts describing each physical mode
    """
    n = len(C_values)
    ndof = 3 * n

    # Step 1: Hessian
    H = compute_hessian_generic(positions, potential_fn, dx)

    # Step 2: Mass matrix M^{-1/2}
    M_inv_sqrt = np.zeros(ndof)
    for i in range(n):
        m_inv_sqrt_i = 1.0 / np.sqrt(max(C_values[i], 1e-30))
        M_inv_sqrt[3*i:3*i+3] = m_inv_sqrt_i
    M_inv_sqrt_diag = np.diag(M_inv_sqrt)

    # Dynamical matrix: D = M^{-1/2} H M^{-1/2}
    D = M_inv_sqrt_diag @ H @ M_inv_sqrt_diag

    # Step 3: Zero-mode projection
    Z, n_trans, n_rot, is_collinear = _build_zero_mode_basis(positions)
    n_zero_modes = Z.shape[1]

    # Mass-weighted zero modes: Z_mw = M^{1/2} Z, then re-orthonormalize
    M_sqrt = np.zeros(ndof)
    for i in range(n):
        M_sqrt[3*i:3*i+3] = np.sqrt(max(C_values[i], 1e-30))
    Z_mw_raw = np.diag(M_sqrt) @ Z
    # Orthonormalize
    Q, R_qr = np.linalg.qr(Z_mw_raw)
    # Keep only columns with nonzero norm
    norms = np.abs(np.diag(R_qr))
    keep = norms > 1e-10
    Z_mw = Q[:, keep]
    n_zero_mw = Z_mw.shape[1]

    # Projector in mass-weighted space
    P = np.eye(ndof) - Z_mw @ Z_mw.T
    D_phys = P.T @ D @ P

    # Step 4: Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(D_phys)

    # Extract physical modes (skip the n_zero_mw smallest)
    idx_sorted = np.argsort(eigenvalues)
    eigenvalues_sorted = eigenvalues[idx_sorted]
    eigenvectors_sorted = eigenvectors[:, idx_sorted]

    omega2_physical = eigenvalues_sorted[n_zero_mw:]
    evecs_physical = eigenvectors_sorted[:, n_zero_mw:]

    # Convert eigenvectors back from mass-weighted to Cartesian
    evecs_cartesian = M_inv_sqrt_diag @ evecs_physical

    # Classify physical modes
    # Use a minimum absolute threshold to distinguish genuine eigenvalues
    # from numerical noise in the finite-difference Hessian.
    # With dx=1e-5, the Hessian precision is ~V/dx^2 * machine_eps ~ 1e-6,
    # so eigenvalues below ~1e-6 are unreliable.
    ABS_THRESH = 1e-6
    phys_max = max(np.abs(omega2_physical).max(), 1e-10) if len(omega2_physical) > 0 else 1e-10
    phys_thresh = max(phys_max * 1e-8, ABS_THRESH)

    n_stable = int(np.sum(omega2_physical > phys_thresh))
    n_unstable = int(np.sum(omega2_physical < -phys_thresh))
    n_marginal = len(omega2_physical) - n_stable - n_unstable

    if len(omega2_physical) == 0:
        classification = "degenerate"
    elif n_stable > 0 and n_unstable == 0:
        classification = "stable"
    elif n_unstable > 0 and n_stable == 0:
        classification = "unstable"
    elif n_stable > 0 and n_unstable > 0:
        classification = "saddle"
    else:
        # All eigenvalues are below threshold — flat/degenerate potential
        classification = "marginal"

    # Oscillation frequencies
    pos_omega2 = omega2_physical[omega2_physical > phys_thresh]
    omega_physical = np.sqrt(pos_omega2)

    # Mode character analysis
    mode_characters = []
    for k in range(len(omega2_physical)):
        evec = evecs_cartesian[:, k]
        # Decompose into per-body contributions
        body_amps = np.array([np.linalg.norm(evec[3*i:3*i+3]) for i in range(n)])
        # Radial vs tangential (relative to CoM)
        com = np.mean(positions, axis=0)
        radial_frac = 0.0
        total_amp = 0.0
        for i in range(n):
            r_hat = positions[i] - com
            r_norm = np.linalg.norm(r_hat)
            if r_norm > 1e-12:
                r_hat = r_hat / r_norm
                radial_comp = abs(np.dot(evec[3*i:3*i+3], r_hat))
            else:
                radial_comp = 0.0
            radial_frac += radial_comp**2
            total_amp += np.dot(evec[3*i:3*i+3], evec[3*i:3*i+3])

        radial_frac = radial_frac / max(total_amp, 1e-30)

        if radial_frac > 0.8:
            mode_type = "breathing"
        elif radial_frac < 0.2:
            mode_type = "tangential"
        else:
            mode_type = "mixed"

        mode_characters.append({
            'omega2': omega2_physical[k],
            'omega': np.sqrt(abs(omega2_physical[k])),
            'stable': omega2_physical[k] > phys_thresh,
            'type': mode_type,
            'radial_fraction': radial_frac,
            'body_amplitudes': body_amps,
        })

    return {
        'omega2_physical': omega2_physical,
        'omega_physical': np.sort(omega_physical),
        'eigenvectors_physical': evecs_cartesian,
        'H': H,
        'D': D,
        'n_stable': n_stable,
        'n_unstable': n_unstable,
        'n_zero_modes': n_zero_mw,
        'classification': classification,
        'is_collinear': is_collinear,
        'mode_characters': mode_characters,
    }


def stability_comparison(positions, C_values, beta, gamma=None,
                         softening=1e-6, G_newton=None,
                         include_3body=False, dx=1e-5):
    """
    Compare TGP vs Newton stability at the same configuration.

    For fair comparison, G_newton = 4*pi by default (matching TGP leading term).

    Parameters
    ----------
    positions : (n, 3) ndarray
    C_values : (n,) ndarray
    beta, gamma : float
    softening : float
    G_newton : float or None
        If None, uses 4*pi (TGP-matched).
    include_3body : bool
        If True, include V3 in TGP potential.
    dx : float
        Finite difference step.

    Returns
    -------
    dict with:
        newton : normal_mode_analysis result for Newton
        tgp : normal_mode_analysis result for TGP
        omega2_ratio : element-wise ratio of physical eigenvalues
        stability_differs : bool (different classification)
        n_modes_stabilized : count of modes stable in TGP but not Newton
        n_modes_destabilized : count of modes stable in Newton but not TGP
    """
    if gamma is None:
        gamma = beta
    if G_newton is None:
        G_newton = 4.0 * np.pi

    from .dynamics_v2 import potential_newton, potential_tgp, full_potential_tgp

    # Newton potential wrapper
    def pot_newton(pos):
        return potential_newton(pos, C_values, G=G_newton, softening=softening)

    # TGP potential wrapper
    if include_3body:
        def pot_tgp(pos):
            return full_potential_tgp(pos, C_values, beta, gamma,
                                      softening=softening, include_3body=True,
                                      use_yukawa=False)
    else:
        def pot_tgp(pos):
            return potential_tgp(pos, C_values, beta, gamma, softening=softening)

    # Normal mode analysis for both
    nma_newton = normal_mode_analysis(positions, C_values, pot_newton, dx)
    nma_tgp = normal_mode_analysis(positions, C_values, pot_tgp, dx)

    # Compare eigenvalues
    o2_N = nma_newton['omega2_physical']
    o2_T = nma_tgp['omega2_physical']

    n_phys = min(len(o2_N), len(o2_T))
    omega2_ratio = np.full(n_phys, np.nan)
    for k in range(n_phys):
        if abs(o2_N[k]) > 1e-10:
            omega2_ratio[k] = o2_T[k] / o2_N[k]

    # Count stabilized/destabilized modes
    thresh_N = max(np.abs(o2_N).max(), 1e-10) * 1e-8 if len(o2_N) > 0 else 1e-10
    thresh_T = max(np.abs(o2_T).max(), 1e-10) * 1e-8 if len(o2_T) > 0 else 1e-10

    n_stabilized = 0
    n_destabilized = 0
    for k in range(n_phys):
        newton_stable = o2_N[k] > thresh_N
        tgp_stable = o2_T[k] > thresh_T
        if tgp_stable and not newton_stable:
            n_stabilized += 1
        if newton_stable and not tgp_stable:
            n_destabilized += 1

    return {
        'newton': nma_newton,
        'tgp': nma_tgp,
        'omega2_ratio': omega2_ratio,
        'stability_differs': nma_newton['classification'] != nma_tgp['classification'],
        'n_modes_stabilized': n_stabilized,
        'n_modes_destabilized': n_destabilized,
    }


def stability_bifurcation_scan(config_fn, betas, C_values_or_Cs,
                                gamma_fn=None, softening=1e-6,
                                include_3body=False, dx=1e-5):
    """
    Scan stability classification over a grid of (beta, C) values.

    For each (beta, C) point:
    1. Find equilibrium configuration using config_fn
    2. Compute normal modes for TGP and Newton
    3. Record classification and key eigenvalues

    Parameters
    ----------
    config_fn : callable
        config_fn(C, beta, gamma) -> (positions, C_values) or None.
        Returns equilibrium positions or None if no equilibrium exists.
    betas : array-like
        Beta values to scan.
    C_values_or_Cs : array-like
        C values to scan. Each C is used as equal mass for all bodies.
    gamma_fn : callable or None
        gamma(beta) -> float. If None, gamma = beta.
    softening : float
    include_3body : bool
    dx : float

    Returns
    -------
    list of dicts, one per (beta, C) grid point, with:
        beta, C, gamma, has_equilibrium,
        tgp_classification, newton_classification,
        tgp_n_stable, newton_n_stable,
        tgp_omega2_min, newton_omega2_min,
        n_modes_stabilized, n_modes_destabilized
    """
    results = []

    for beta in betas:
        gamma = gamma_fn(beta) if gamma_fn else beta

        for C in C_values_or_Cs:
            config = config_fn(C, beta, gamma)
            if config is None:
                results.append({
                    'beta': beta, 'C': C, 'gamma': gamma,
                    'has_equilibrium': False,
                })
                continue

            positions, C_vals = config

            try:
                comp = stability_comparison(
                    positions, C_vals, beta, gamma,
                    softening=softening, include_3body=include_3body, dx=dx,
                )
            except Exception:
                results.append({
                    'beta': beta, 'C': C, 'gamma': gamma,
                    'has_equilibrium': True, 'error': True,
                })
                continue

            tgp_nma = comp['tgp']
            newt_nma = comp['newton']

            tgp_o2 = tgp_nma['omega2_physical']
            newt_o2 = newt_nma['omega2_physical']

            results.append({
                'beta': beta, 'C': C, 'gamma': gamma,
                'has_equilibrium': True,
                'tgp_classification': tgp_nma['classification'],
                'newton_classification': newt_nma['classification'],
                'tgp_n_stable': tgp_nma['n_stable'],
                'tgp_n_unstable': tgp_nma['n_unstable'],
                'newton_n_stable': newt_nma['n_stable'],
                'newton_n_unstable': newt_nma['n_unstable'],
                'tgp_omega2_min': float(tgp_o2.min()) if len(tgp_o2) > 0 else np.nan,
                'tgp_omega2_max': float(tgp_o2.max()) if len(tgp_o2) > 0 else np.nan,
                'newton_omega2_min': float(newt_o2.min()) if len(newt_o2) > 0 else np.nan,
                'newton_omega2_max': float(newt_o2.max()) if len(newt_o2) > 0 else np.nan,
                'n_modes_stabilized': comp['n_modes_stabilized'],
                'n_modes_destabilized': comp['n_modes_destabilized'],
                'stability_differs': comp['stability_differs'],
                'tgp_mode_characters': tgp_nma['mode_characters'],
            })

    return results


def earnshaw_test(C, beta, gamma=None):
    """
    Test whether TGP violates Earnshaw's theorem for specific configurations.

    STATUS: uses PAIRWISE energy only (no 3-body terms in Hessian).
    For full 2b+3b Earnshaw test, pass a full energy function to
    compute_hessian_numerical() directly.

    Uses rigorous zero-mode projection (not threshold heuristic).

    Parameters
    ----------
    C : float — source strength
    beta : float

    Returns
    -------
    dict with Earnshaw analysis for equilateral configurations.
    """
    if gamma is None:
        gamma = beta

    from . import equilibria

    results = {}

    # Equilateral triangle
    eq_tri = equilibria.equilateral_pairwise_equilibrium(C, beta, gamma)
    if eq_tri['exists']:
        for label, d in [('d_rep', eq_tri['d_rep']), ('d_well', eq_tri['d_well'])]:
            if d is None or d <= 0:
                continue
            # Build positions
            from . import configurations
            pos, Cv, _ = configurations.equilateral_triangle(d, C)
            H = compute_hessian_numerical(pos, Cv, beta, gamma)
            trace = np.trace(H)
            stab = classify_stability(H, positions=pos)

            results[f'equilateral_{label}'] = {
                'd': d,
                'trace_H': trace,
                'earnshaw_violated': trace > 0,  # positive trace = at least one stable direction
                'stability': stab['classification'],
                'eigenvalues_physical': stab['eigenvalues_physical'],
                'eigenvalues_full': stab['eigenvalues_full'],
            }

    return results
