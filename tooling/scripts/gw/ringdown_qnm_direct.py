#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ringdown_qnm_direct.py -- TGP QNM via Time-Domain Integration
==============================================================
Computes QNM spectrum for l=0 and l=2 perturbations using
direct time-domain evolution of the wave equation:

    d^2 Psi/dt^2 = d^2 Psi/dr*^2 - V_eff(r*) Psi

Method:
  1. Solve background profile f(r) via BVP
  2. Build V_eff(r*) on tortoise grid
  3. Evolve Gaussian initial pulse via leapfrog finite differences
  4. Extract QNM frequencies via:
     a) FFT of late-time signal
     b) Prony analysis (exponential fit)
  5. Compare l=0 and l=2

Reference: dodatekC_ringdown.tex, thm:ringdown
"""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import os
import numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ===================================================================
# PARAMETERS
# ===================================================================
ALPHA = 2.0
SOURCE_SIGN = -1
X_MAX = 60.0
N_MESH = 800

# Background profiles to compute
STRONG_PARAMS = [
    {"bg": 0.01, "S":  5.0,  "sigma": 0.5, "label": "S=5"},
    {"bg": 0.01, "S": 10.0,  "sigma": 0.5, "label": "S=10"},
    {"bg": 0.01, "S": 20.0,  "sigma": 0.5, "label": "S=20"},
]

# Time-domain evolution parameters
N_RSTAR = 3000       # spatial grid points in r*
T_MAX = 600.0        # total evolution time (long enough for damping)
CFL = 0.5            # Courant factor (optimal for leapfrog + Mur ABC)
R_OBS_FRAC = 0.5     # observation point as fraction of r*_max

# ===================================================================
# 1. RADIAL PROFILE SOLVER (same as ringdown_qnm.py)
# ===================================================================

def source_gaussian(x, S, sigma):
    return S * np.exp(-x**2 / (2 * sigma**2)) / (2 * np.pi * sigma**2)**1.5

def ode_fun(x, y, bg, S, sigma):
    f = np.maximum(y[0], 1e-15)
    fp = y[1]
    x_safe = np.maximum(x, 1e-10)
    src = SOURCE_SIGN * source_gaussian(x, S, sigma)
    fpp = (-(2.0 / x_safe) * fp
           - ALPHA * fp**2 / f
           - bg * f**2 + bg * f**3
           + src)
    return np.vstack([fp, fpp])

def bc_fun(ya, yb):
    return np.array([ya[1], yb[0] - 1.0])

def solve_profile(bg, S, sigma):
    x_min = 0.01
    x_mesh = np.linspace(x_min, X_MAX, N_MESH)
    peak = S / (4 * np.pi * sigma)
    f_guess = 1.0 + peak * np.exp(-x_mesh**2 / (2 * sigma**2))
    fp_guess = -peak * x_mesh / sigma**2 * np.exp(-x_mesh**2 / (2 * sigma**2))
    y_guess = np.vstack([f_guess, fp_guess])

    for tol, nodes in [(1e-6, 50000), (1e-4, 100000), (1e-3, 150000)]:
        sol = solve_bvp(
            fun=lambda x, y: ode_fun(x, y, bg, S, sigma),
            bc=bc_fun, x=x_mesh, y=y_guess, tol=tol,
            max_nodes=nodes, verbose=0,
        )
        if sol.success:
            break
        x_mesh = np.linspace(x_min, X_MAX, 2 * len(x_mesh))
        f_guess = 1.0 + peak * np.exp(-x_mesh**2 / (2 * sigma**2))
        fp_guess = -peak * x_mesh / sigma**2 * np.exp(-x_mesh**2 / (2 * sigma**2))
        y_guess = np.vstack([f_guess, fp_guess])

    if not sol.success:
        print(f"    [FAIL] BVP did not converge: {sol.message}")
        return None

    x = sol.x
    f = sol.y[0]
    fp = sol.y[1]
    rhs = ode_fun(x, sol.y, bg, S, sigma)
    fpp = rhs[1]
    return x, f, fp, fpp

# ===================================================================
# 2. EFFECTIVE POTENTIAL & TORTOISE COORDINATE
# ===================================================================

def compute_V_eff(x, f, fp, fpp, bg, ell=0):
    """V_eff from thm:ringdown (eq:V-eff in dodatekC)."""
    f_safe = np.maximum(f, 1e-15)
    x_safe = np.maximum(x, 1e-10)
    chi = f_safe

    bracket = (ell * (ell + 1) / x_safe**2
               + (9.0/4.0) * fpp / f_safe
               + (59.0/16.0) * (fp / f_safe)**2
               + 4.0 * fp / (x_safe * f_safe)
               - 2.0 * bg * chi
               + 3.0 * bg * chi**2)

    V = bracket / f_safe
    return V

def compute_tortoise(x, f):
    """r* = integral sqrt(f) dx."""
    integrand = np.sqrt(np.maximum(f, 1e-15))
    r_star = np.zeros_like(x)
    for i in range(1, len(x)):
        r_star[i] = r_star[i-1] + 0.5 * (integrand[i] + integrand[i-1]) * (x[i] - x[i-1])
    return r_star

def build_potential_on_tortoise(x, f, fp, fpp, bg, ell, n_rstar):
    """
    Build V_eff on a uniform r* grid using interpolation.
    Returns: r_star_uniform, V_eff_uniform
    """
    V = compute_V_eff(x, f, fp, fpp, bg, ell=ell)
    r_star = compute_tortoise(x, f)

    # Clip extreme values for numerical stability
    V_clipped = np.clip(V, -100.0, 100.0)

    # Remove any non-finite values
    mask = np.isfinite(V_clipped) & np.isfinite(r_star)
    r_star_clean = r_star[mask]
    V_clean = V_clipped[mask]

    # Ensure monotonic r*
    dr = np.diff(r_star_clean)
    mono_mask = np.concatenate([[True], dr > 0])
    r_star_clean = r_star_clean[mono_mask]
    V_clean = V_clean[mono_mask]

    if len(r_star_clean) < 10:
        return None, None

    # Interpolate onto uniform r* grid
    r_star_min = r_star_clean[2]   # skip very first points (numerically noisy)
    r_star_max = r_star_clean[-1]

    cs = CubicSpline(r_star_clean, V_clean, extrapolate=True)
    r_star_uni = np.linspace(r_star_min, r_star_max, n_rstar)
    V_uni = cs(r_star_uni)

    # Smooth boundaries to avoid reflections
    # Apply absorbing taper at both ends
    taper_width = int(n_rstar * 0.05)
    taper_left = np.linspace(0, 1, taper_width)**2
    taper_right = np.linspace(1, 0, taper_width)**2
    # Don't taper V itself -- taper is for the wave function (handled by ABCs)

    return r_star_uni, V_uni

# ===================================================================
# 3. TIME-DOMAIN EVOLUTION (Leapfrog / Verlet scheme)
# ===================================================================

def evolve_wave(r_star, V_eff, t_max, cfl=0.5, r0_pulse=None, sigma_pulse=None):
    """
    Evolve the wave equation:
        d^2 Psi/dt^2 = d^2 Psi/dr*^2 - V_eff(r*) Psi

    using leapfrog finite differences with:
    - Inner boundary: hard wall (reflecting) at V<0/V>0 transition
      (physically: boundary of linearization-valid region)
    - Outer boundary: Sommerfeld absorbing BC (outgoing wave condition)
      dPsi/dt + dPsi/dr* = 0   =>   Psi_N^{n+1} from one-way wave eq.
    - Outer sponge layer: gentle damping zone as backup absorption

    Returns: t_array, Psi_obs(t) at observation point,
             also returns auxiliary info dict.
    """
    # --- 1. Truncate inner V < 0 region ---
    V_work = V_eff.copy()
    i_inner_wall = 0
    for i in range(len(V_work) - 1):
        if V_work[i] < 0 and V_work[i+1] >= 0:
            i_inner_wall = i + 1
    if i_inner_wall > 0:
        r_star = r_star[i_inner_wall:]
        V_work = V_work[i_inner_wall:]

    N = len(r_star)
    if N < 100:
        return np.array([]), np.array([]), {}

    dr = r_star[1] - r_star[0]
    dt = cfl * dr
    Nt = int(t_max / dt)

    # --- 2. No sponge layer (pure Sommerfeld ABC at outer boundary) ---
    sponge_width = 0  # no sponge, purely Sommerfeld

    # --- 3. Observation point ---
    i_obs_target = int(N * 0.4)
    i_obs = max(10, min(i_obs_target, N - 20))

    # Also observe at a second point for consistency check
    i_obs2 = max(10, min(int(N * 0.3), N - 20))

    # --- 4. Initial pulse ---
    if r0_pulse is None:
        # Place pulse between inner wall and observation point
        r0_pulse = r_star[max(5, N // 5)]
    if sigma_pulse is None:
        sigma_pulse = (r_star[-1] - r_star[0]) * 0.02

    Psi = np.exp(-(r_star - r0_pulse)**2 / (2 * sigma_pulse**2))
    Psi_prev = Psi.copy()

    # Init for zero velocity
    d2Psi = np.zeros(N)
    d2Psi[1:-1] = (Psi[2:] - 2*Psi[1:-1] + Psi[:-2]) / dr**2
    accel = d2Psi - V_work * Psi
    Psi_prev = Psi - 0.5 * dt**2 * accel

    # Storage
    t_store = []
    psi_store = []
    psi_store2 = []
    store_every = max(1, Nt // 20000)

    ratio = (dt / dr)**2
    c_ratio = cfl  # dt/dr = CFL, used in Sommerfeld BC

    for n in range(Nt):
        Psi_next = np.zeros(N)

        # Interior: standard leapfrog (no sponge)
        Psi_next[1:-1] = (2 * Psi[1:-1] - Psi_prev[1:-1]
                          + ratio * (Psi[2:] - 2*Psi[1:-1] + Psi[:-2])
                          - dt**2 * V_work[1:-1] * Psi[1:-1])

        # Inner BC: hard wall (reflecting -- physically correct for rising V)
        Psi_next[0] = 0.0

        # Outer BC: Sommerfeld absorbing condition
        # One-way wave equation: dPsi/dt + dPsi/dr* = 0
        # Discretized: Psi_N^{n+1} = Psi_{N-1}^{n}
        #              + (CFL - 1)/(CFL + 1) * (Psi_{N-1}^{n+1} - Psi_N^{n})
        Psi_next[-1] = (Psi[-2]
                        + (c_ratio - 1) / (c_ratio + 1) * (Psi_next[-2] - Psi[-1]))

        if n % store_every == 0:
            t_store.append(n * dt)
            psi_store.append(Psi_next[i_obs])
            psi_store2.append(Psi_next[i_obs2])

        Psi_prev = Psi.copy()
        Psi = Psi_next.copy()

        if np.max(np.abs(Psi)) > 1e10:
            print(f"    [WARN] Blow-up at t = {n*dt:.1f}, truncating")
            break

    info = {
        "i_inner_wall": i_inner_wall,
        "i_obs": i_obs,
        "i_obs2": i_obs2,
        "r_obs": r_star[i_obs],
        "r_obs2": r_star[i_obs2],
        "r_pulse": r0_pulse,
        "N_eff": N,
        "dt": dt,
        "dr": dr,
        "psi2": np.array(psi_store2),
    }
    return np.array(t_store), np.array(psi_store), info

# ===================================================================
# 4. FREQUENCY EXTRACTION
# ===================================================================

def extract_qnm_fft(t, psi, t_start_frac=0.05, t_end_frac=0.6):
    """
    Extract dominant frequency from ringdown signal via FFT.
    Uses early-to-mid portion of signal (after initial transient,
    before signal decays too much or boundary effects dominate).
    """
    N = len(t)
    i_start = int(N * t_start_frac)
    i_end = int(N * t_end_frac)
    t_ring = t[i_start:i_end]
    psi_ring = psi[i_start:i_end]

    if len(psi_ring) < 64:
        return None, None, None

    # Apply Blackman window (better sidelobe suppression than Hanning)
    window = np.blackman(len(psi_ring))
    psi_windowed = psi_ring * window

    dt = t_ring[1] - t_ring[0]

    # Zero-pad for better frequency resolution
    n_fft = max(len(psi_windowed), 2**14)
    freqs = np.fft.rfftfreq(n_fft, d=dt)
    spectrum = np.abs(np.fft.rfft(psi_windowed, n=n_fft))

    if len(spectrum) < 2:
        return None, None, None

    # Skip DC and very low frequencies (below mass gap)
    omega_arr = 2 * np.pi * freqs
    omega_min_gap = 0.05  # slightly below mass gap
    mask_low = omega_arr > omega_min_gap
    spectrum_masked = spectrum.copy()
    spectrum_masked[~mask_low] = 0

    i_peak = np.argmax(spectrum_masked)
    f_peak = freqs[i_peak]
    omega_re = 2 * np.pi * f_peak

    return omega_re, freqs, spectrum

def extract_qnm_prony(t, psi, t_start_frac=0.05, t_end_frac=0.5, n_modes=6):
    """
    Extract QNM frequencies using Prony's method (Matrix Pencil variant).
    Fits sum of damped exponentials: Psi(t) = sum A_k exp(s_k t)
    where s_k = -gamma_k + i*omega_k

    Uses early ringdown phase (after initial transient, before decay).
    Returns list of (omega_re, omega_im, amplitude) tuples, deduplicated.
    """
    N = len(t)
    i_start = int(N * t_start_frac)
    i_end = int(N * t_end_frac)
    t_ring = t[i_start:i_end]
    psi_ring = psi[i_start:i_end]

    if len(psi_ring) < 4 * n_modes:
        return []

    # Subsample for efficiency, but keep enough points
    target_points = 600
    step = max(1, len(psi_ring) // target_points)
    t_s = t_ring[::step]
    y_s = psi_ring[::step]
    M = len(y_s)

    if M < 2 * n_modes + 2:
        return []

    dt_s = t_s[1] - t_s[0]

    # Matrix Pencil Method (more robust than standard Prony)
    p = n_modes
    L = M - p

    # Build data Hankel matrices Y0 and Y1
    Y0 = np.zeros((L, p))
    Y1 = np.zeros((L, p))
    for i in range(L):
        for j in range(p):
            Y0[i, j] = y_s[i + j]
            Y1[i, j] = y_s[i + j + 1]

    # Solve generalized eigenvalue problem: Y1 v = z Y0 v
    # via pseudoinverse: z = eigenvalues of pinv(Y0) @ Y1
    try:
        Y0_pinv = np.linalg.pinv(Y0, rcond=1e-8)
        pencil = Y0_pinv @ Y1
        eigenvalues = np.linalg.eigvals(pencil)
    except np.linalg.LinAlgError:
        return []

    # Convert eigenvalues z_k to complex frequencies s_k = ln(z_k)/dt
    results = []
    seen = set()
    for z in eigenvalues:
        if np.abs(z) < 1e-12 or np.abs(z) > 2.0:
            continue  # unphysical
        s = np.log(z + 0j) / dt_s
        omega_re = np.abs(s.imag)
        omega_im = s.real  # negative = damped

        if omega_re < 0.01:  # skip near-DC
            continue
        if omega_im > 0:  # skip growing modes (unphysical)
            continue

        # Deduplicate conjugate pairs (within tolerance)
        key = (round(omega_re, 3), round(omega_im, 3))
        if key in seen:
            continue
        seen.add(key)

        # Estimate amplitude from Vandermonde fit
        amp = np.abs(z)

        results.append((omega_re, omega_im, amp))

    # Sort by damping rate (least damped first = most visible)
    results.sort(key=lambda x: x[1], reverse=True)
    return results

def extract_damping_rate(t, psi, t_start_frac=0.05, t_end_frac=0.5):
    """
    Extract damping rate from envelope of |Psi(t)| during ringdown phase.
    Fits log|Psi_peaks| = -gamma*t + const.
    Uses only the ringdown phase (not late-time noise).
    """
    N = len(t)
    i_start = int(N * t_start_frac)
    i_end = int(N * t_end_frac)
    t_ring = t[i_start:i_end]
    psi_ring = np.abs(psi[i_start:i_end])

    # Find envelope via peaks of |Psi|
    min_dist = max(3, len(psi_ring) // 200)
    peaks, _ = find_peaks(psi_ring, distance=min_dist)
    if len(peaks) < 4:
        peaks, _ = find_peaks(psi_ring, distance=2)
    if len(peaks) < 4:
        return None

    t_peaks = t_ring[peaks]
    amp_peaks = psi_ring[peaks]

    # Remove near-zero values
    mask = amp_peaks > np.max(amp_peaks) * 1e-6
    t_peaks = t_peaks[mask]
    amp_peaks = amp_peaks[mask]

    if len(t_peaks) < 4:
        return None

    # Linear fit to log(amplitude)
    log_amp = np.log(amp_peaks)
    try:
        coeffs = np.polyfit(t_peaks, log_amp, 1)
        gamma = -coeffs[0]  # damping rate (positive = damped)
        if gamma <= 0:  # must be damped
            return None
        return gamma
    except Exception:
        return None

# ===================================================================
# 5. MAIN
# ===================================================================

def main():
    print("=" * 70)
    print("  TGP RINGDOWN: DIRECT TIME-DOMAIN INTEGRATION")
    print("  QNM for l=0 and l=2 (bypassing WKB limitation)")
    print("=" * 70)

    save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(save_dir, exist_ok=True)

    all_results = []

    for ps in STRONG_PARAMS:
        bg, S, sigma = ps["bg"], ps["S"], ps["sigma"]
        label = ps["label"]
        print(f"\n{'='*60}")
        print(f"  Profile: {label}")
        print(f"{'='*60}")

        # Solve background
        out = solve_profile(bg, S, sigma)
        if out is None:
            print("    SKIPPED (solver failed)")
            continue
        x, f, fp, fpp = out
        print(f"    f_max = {f.max():.4f} at x = {x[np.argmax(f)]:.3f}")

        for ell in [0, 2]:
            print(f"\n  --- l = {ell} ---")

            # Build potential on tortoise grid
            r_star_uni, V_uni = build_potential_on_tortoise(
                x, f, fp, fpp, bg, ell, N_RSTAR)

            if r_star_uni is None:
                print("    SKIPPED (potential construction failed)")
                continue

            # Clip potential for stability (very large V near origin)
            # For l=2 the centrifugal barrier is huge near r=0;
            # clip it to allow stable evolution
            V_max_clip = 5.0
            V_uni_clip = np.clip(V_uni, -5.0, V_max_clip)

            dr_star = r_star_uni[1] - r_star_uni[0]
            print(f"    r* range: [{r_star_uni[0]:.2f}, {r_star_uni[-1]:.2f}]")
            print(f"    dr* = {dr_star:.4f}")
            print(f"    V range: [{V_uni_clip.min():.4f}, {V_uni_clip.max():.4f}]")

            # Evolve
            print(f"    Evolving wave equation (T_max = {T_MAX}) ...")
            t_arr, psi_arr, evo_info = evolve_wave(
                r_star_uni, V_uni_clip, T_MAX, cfl=CFL)
            print(f"    Evolution complete: {len(t_arr)} time samples")
            if len(psi_arr) > 0:
                print(f"    max |Psi_obs| = {np.max(np.abs(psi_arr)):.6e}")
                print(f"    r*_obs = {evo_info.get('r_obs', '?'):.2f}, "
                      f"r*_pulse = {evo_info.get('r_pulse', '?'):.2f}, "
                      f"N_eff = {evo_info.get('N_eff', '?')}")

            if np.max(np.abs(psi_arr)) < 1e-30:
                print("    Signal too weak -- no QNM detected")
                continue

            # Extract frequencies via FFT
            omega_fft, freqs, spectrum = extract_qnm_fft(t_arr, psi_arr)
            if omega_fft is not None:
                print(f"    FFT: omega_re = {omega_fft:.6f}")
            else:
                print(f"    FFT: no clear peak")

            # Extract via Prony
            prony_modes = extract_qnm_prony(t_arr, psi_arr, n_modes=4)
            if prony_modes:
                print(f"    Prony analysis ({len(prony_modes)} damped modes found):")
                for k, (w_re, w_im, amp) in enumerate(prony_modes[:3]):
                    Q = -w_re / (2 * w_im) if w_im != 0 else float('inf')
                    print(f"      mode {k}: omega = {w_re:.5f} {w_im:+.5f}i  "
                          f"Q = {Q:.2f}  amp = {amp:.4f}")
            else:
                print(f"    Prony: no damped modes found")

            # Extract damping from envelope
            gamma_env = extract_damping_rate(t_arr, psi_arr)
            if gamma_env is not None:
                print(f"    Envelope damping: gamma = {gamma_env:.6f}")
                if omega_fft is not None:
                    omega_complex = complex(omega_fft, -gamma_env)
                    Q_env = omega_fft / (2 * gamma_env)
                    print(f"    Combined: omega = {omega_fft:.5f} {-gamma_env:+.5f}i"
                          f"  Q = {Q_env:.2f}")

            # Cross-check: Prony on second observation point
            psi2 = evo_info.get("psi2", np.array([]))
            if len(psi2) == len(t_arr):
                prony2 = extract_qnm_prony(t_arr, psi2, n_modes=4)
                if prony2:
                    pm2 = prony2[0]
                    Q2 = -pm2[0] / (2 * pm2[1]) if pm2[1] != 0 else float('inf')
                    print(f"    Cross-check (r*={evo_info.get('r_obs2',0):.1f}): "
                          f"omega = {pm2[0]:.5f} {pm2[1]:+.5f}i  Q = {Q2:.2f}")
                    # Check consistency
                    if prony_modes:
                        pm1 = prony_modes[0]
                        delta_w = abs(pm1[0] - pm2[0]) / max(pm1[0], 0.01)
                        status = "CONSISTENT" if delta_w < 0.15 else "INCONSISTENT"
                        print(f"    Frequency match: delta_omega/omega = {delta_w:.3f} ({status})")

            all_results.append({
                "label": label, "S": S, "ell": ell,
                "r_star": r_star_uni, "V_eff": V_uni_clip,
                "t": t_arr, "psi": psi_arr,
                "omega_fft": omega_fft,
                "freqs": freqs, "spectrum": spectrum,
                "prony_modes": prony_modes,
                "gamma_env": gamma_env,
                "evo_info": evo_info,
            })

    if not all_results:
        print("\nNo results.")
        return

    # =================================================================
    # PLOTS
    # =================================================================
    print(f"\n  Generating plots ...")

    # --- Fig 1: Time-domain signal for all cases ---
    fig1, axes1 = plt.subplots(len(STRONG_PARAMS), 2, figsize=(16, 4*len(STRONG_PARAMS)),
                                squeeze=False)
    fig1.suptitle("TGP Ringdown: Time-Domain Signals", fontsize=15, y=1.01)

    for res in all_results:
        row = [i for i, ps in enumerate(STRONG_PARAMS) if ps["label"] == res["label"]]
        if not row:
            continue
        row = row[0]
        col = res["ell"] // 2  # 0 for l=0, 1 for l=2

        ax = axes1[row, col]
        t_plot = res["t"]
        psi_plot = res["psi"]

        ax.plot(t_plot, psi_plot, lw=0.5, color="C0")
        ax.set_xlabel(r"$t$", fontsize=11)
        ax.set_ylabel(r"$\Psi(t, r_{\rm obs})$", fontsize=11)
        ax.set_title(f"{res['label']}, $\\ell = {res['ell']}$", fontsize=12)
        ax.grid(True, ls=":", alpha=0.4)

        # Mark damping envelope
        if res["gamma_env"] is not None and res["omega_fft"] is not None:
            A0 = np.max(np.abs(psi_plot[:len(psi_plot)//3]))
            t_env = np.linspace(t_plot[0], t_plot[-1], 200)
            env = A0 * np.exp(-res["gamma_env"] * (t_env - t_plot[0]))
            ax.plot(t_env, env, "r--", lw=1.2, alpha=0.7,
                    label=rf"$\gamma = {res['gamma_env']:.4f}$")
            ax.plot(t_env, -env, "r--", lw=1.2, alpha=0.7)
            ax.legend(fontsize=9)

    fig1.tight_layout()
    path1 = os.path.join(save_dir, "ringdown_timedomain.png")
    fig1.savefig(path1, dpi=150, bbox_inches='tight')
    print(f"    Saved {path1}")

    # --- Fig 2: FFT spectra ---
    fig2, axes2 = plt.subplots(1, 2, figsize=(14, 5))
    fig2.suptitle("TGP Ringdown: Frequency Spectra", fontsize=14)

    for ell_idx, ell in enumerate([0, 2]):
        ax = axes2[ell_idx]
        for res in all_results:
            if res["ell"] != ell or res["spectrum"] is None:
                continue
            omega_arr = 2 * np.pi * res["freqs"]
            ax.plot(omega_arr, res["spectrum"] / np.max(res["spectrum"]),
                    lw=1.5, label=res["label"])
            if res["omega_fft"] is not None:
                ax.axvline(res["omega_fft"], ls="--", lw=0.8, color="gray")

        # Mass gap
        ax.axvline(np.sqrt(0.01), color="red", ls=":", lw=1.0,
                   label=r"$\omega_{\rm min} = \sqrt{\gamma}$")
        ax.set_xlabel(r"$\omega$", fontsize=12)
        ax.set_ylabel("Normalized amplitude", fontsize=12)
        ax.set_title(rf"$\ell = {ell}$", fontsize=13)
        ax.set_xlim(0, 2.0)
        ax.legend(fontsize=9)
        ax.grid(True, ls=":", alpha=0.4)

    fig2.tight_layout()
    path2 = os.path.join(save_dir, "ringdown_spectra.png")
    fig2.savefig(path2, dpi=150, bbox_inches='tight')
    print(f"    Saved {path2}")

    # --- Fig 3: V_eff comparison l=0 vs l=2 ---
    fig3, axes3 = plt.subplots(1, len(STRONG_PARAMS), figsize=(6*len(STRONG_PARAMS), 5),
                                squeeze=False)
    fig3.suptitle(r"$V_{\rm eff}(r_*)$: $\ell=0$ vs $\ell=2$", fontsize=14)

    for i, ps in enumerate(STRONG_PARAMS):
        ax = axes3[0, i]
        for res in all_results:
            if res["label"] != ps["label"]:
                continue
            V_plot = np.clip(res["V_eff"], -2, 10)
            ax.plot(res["r_star"], V_plot, lw=1.5,
                    label=rf"$\ell={res['ell']}$")
        ax.axhline(0.01, color="gray", ls="--", lw=0.8, label=r"$\gamma$")
        ax.set_xlabel(r"$r_*$", fontsize=12)
        ax.set_ylabel(r"$V_{\rm eff}$", fontsize=12)
        ax.set_title(ps["label"], fontsize=13)
        ax.set_ylim(-1, 5)
        ax.legend(fontsize=10)
        ax.grid(True, ls=":", alpha=0.4)

    fig3.tight_layout()
    path3 = os.path.join(save_dir, "ringdown_V_comparison.png")
    fig3.savefig(path3, dpi=150, bbox_inches='tight')
    print(f"    Saved {path3}")

    # --- Fig 4: log|Psi| (damping visualization) ---
    fig4, axes4 = plt.subplots(1, 2, figsize=(14, 5))
    fig4.suptitle("TGP Ringdown: Damping (log scale)", fontsize=14)

    for ell_idx, ell in enumerate([0, 2]):
        ax = axes4[ell_idx]
        for res in all_results:
            if res["ell"] != ell:
                continue
            mask = np.abs(res["psi"]) > 1e-30
            if np.sum(mask) < 10:
                continue
            ax.semilogy(res["t"][mask], np.abs(res["psi"][mask]),
                        lw=0.8, label=res["label"])
        ax.set_xlabel(r"$t$", fontsize=12)
        ax.set_ylabel(r"$|\Psi(t)|$", fontsize=12)
        ax.set_title(rf"$\ell = {ell}$", fontsize=13)
        ax.legend(fontsize=10)
        ax.grid(True, ls=":", alpha=0.4)

    fig4.tight_layout()
    path4 = os.path.join(save_dir, "ringdown_damping.png")
    fig4.savefig(path4, dpi=150, bbox_inches='tight')
    print(f"    Saved {path4}")

    plt.close("all")

    # =================================================================
    # SUMMARY TABLE
    # =================================================================
    print("\n" + "=" * 70)
    print("  SUMMARY: QNM DIRECT INTEGRATION RESULTS")
    print("=" * 70)

    print(f"\n  {'Source':<10} {'l':<4} {'omega_re(FFT)':<14} "
          f"{'gamma(env)':<12} {'omega(Prony)':<22} {'Q':>6}")
    print("  " + "-" * 68)

    for res in all_results:
        w_fft = f"{res['omega_fft']:.5f}" if res['omega_fft'] else "---"
        g_env = f"{res['gamma_env']:.5f}" if res['gamma_env'] else "---"

        if res["prony_modes"]:
            pm = res["prony_modes"][0]
            w_prony = f"{pm[0]:.4f}{pm[1]:+.4f}i"
            Q_p = -pm[0] / (2 * pm[1]) if pm[1] != 0 else float('inf')
            q_str = f"{Q_p:.2f}"
        else:
            w_prony = "---"
            q_str = "---"

        print(f"  {res['label']:<10} {res['ell']:<4} {w_fft:<14} "
              f"{g_env:<12} {w_prony:<22} {q_str:>6}")

    # Mass gap reminder
    print(f"\n  Mass gap: omega_min = sqrt(gamma) = {np.sqrt(0.01):.4f}")

    # l=2 specific analysis
    print(f"\n  --- l=2 ANALYSIS ---")
    l2_results = [r for r in all_results if r["ell"] == 2]
    if l2_results:
        has_signal = any(np.max(np.abs(r["psi"])) > 1e-20 for r in l2_results)
        if has_signal:
            print("  l=2 modes DETECTED in time-domain evolution!")
            print("  (WKB missed these because potential lacks standard barrier shape)")
        else:
            print("  l=2 signal too weak -- modes may be strongly suppressed")
    else:
        print("  No l=2 data.")

    print("\n  INTERPRETATION:")
    print("  - Time-domain integration does not require barrier structure")
    print("  - Any trapped/resonant oscillation will appear in the signal")
    print("  - Prony method extracts complex frequencies directly")
    print("  - Compare Q factors between l=0 and l=2")
    print()

if __name__ == "__main__":
    main()
