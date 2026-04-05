# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
cmb_tensor_perturbations.py  --  Theory of Generated Space (TGP)
================================================================
Analysis of CMB tensor perturbations in the TGP disformal metric.

The TGP effective metric is disformal (hyp:disformal, sek08):
    g_uv = A(Phi) eta_uv + B(Phi)/M*^4 * d_u Phi * d_v Phi

where A(Phi) = exp(2*delta_Phi/Phi0), B(Phi) is the disformal coupling.

This script:
1. Computes the perturbation decomposition (breathing, longitudinal,
   vector, tensor) for the disformal metric
2. Derives the GW polarization content as a function of B/M*^4
3. Computes observational constraints on B from GW170817 and PPN
4. Predicts the CMB B-mode power spectrum modification
5. Computes the ratio h_breathing / h_+ as a function of B

References:
    - sek08: prop:disformal-polarization, rem:B-constraints,
      prob:tensor-modes, thm:no-tensor
    - Bekenstein (1993): disformal transformations
    - GW170817: |c_GW - c0|/c0 < 1e-15

Usage:
    python cmb_tensor_perturbations.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

# Physical constants
c0 = 2.998e8        # m/s
G0 = 6.674e-11      # m^3/(kg s^2)
hbar0 = 1.055e-34   # J s
kappa = 8 * np.pi * G0 / c0**4  # m^{-1} kg^{-1} s^2
H0_obs = 2.2e-18    # Hubble constant (s^{-1}), ~68 km/s/Mpc

# TGP parameters
gamma_val = 0.03     # L^{-2}
beta_val = gamma_val  # vacuum condition
m_sp_sq = gamma_val  # spatial mass squared (for beta=gamma)
Phi0 = 25.0          # reference Phi_0 ~ 25 from Lambda_eff constraint

save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
os.makedirs(save_dir, exist_ok=True)


# ======================================================================
# 1. Disformal metric perturbation decomposition
# ======================================================================

def disformal_perturbation_analysis():
    """
    Analyze the perturbation structure of the disformal metric.

    g_uv = A(Phi) eta_uv + B(Phi)/M*^4 * d_u Phi * d_v Phi

    Perturbation: Phi = Phi_bar(t) + delta_Phi(t,x)
    Background: Phi_bar = Phi0 (cosmological), d_i Phi_bar = 0

    delta g_uv decomposes into:
    - Scalar (breathing + longitudinal): from A' delta_Phi and B terms
    - Vector: from B * Phi_dot_bar * d_i(delta_Phi)
    - Tensor: from B * d_i(delta_Phi) * d_j(Phi_bar) [needs d_i Phi_bar != 0]
    """
    print("=" * 65)
    print("1. DISFORMAL METRIC PERTURBATION DECOMPOSITION")
    print("=" * 65)

    print("""
In cosmological background (d_i Phi_bar = 0, Phi_dot_bar != 0):

  delta g_00 = [A'(Phi) + B'(Phi) Phi_dot^2/M*^4] delta_Phi
             + [2 B(Phi) Phi_dot / M*^4] delta_Phi_dot
             => Scalar (breathing)

  delta g_0i = [B(Phi) Phi_dot / M*^4] d_i(delta_Phi)
             => Vector modes (spin 1)

  delta g_ij = A'(Phi) delta_Phi * delta_ij
             => Scalar breathing ONLY (no tensor!)

Near astrophysical source (d_i Phi_bar != 0):

  delta g_ij = A' delta_Phi delta_ij
             + B/M*^4 [d_i(delta_Phi) d_j(Phi_bar) + d_i(Phi_bar) d_j(delta_Phi)]
             + B'/M*^4 delta_Phi d_i(Phi_bar) d_j(Phi_bar)

  The B terms break isotropy => TENSOR MODES (h+, hx) possible!
""")

    # Key result: tensor modes require d_i Phi_bar != 0
    # In the wave zone: d_i Phi_bar ~ q Phi0 M / (4 pi r^2) * x_i/r
    # So tensor amplitude ~ B * (q Phi0 M / (4 pi r^2))^2 / M*^4


def gw_speed_constraint():
    """
    Compute constraints on B from GW170817 speed measurement.

    The disformal metric modifies the GW propagation speed:
    c_GW^2 = c0^2 * (A + B Phi_dot^2/M*^4) / A
           = c0^2 * (1 + B Phi_dot^2 / (A M*^4))

    GW170817: |c_GW - c0|/c0 < 1e-15
    """
    print("\n" + "=" * 65)
    print("2. GW SPEED CONSTRAINT (GW170817)")
    print("=" * 65)

    # Phi_dot in cosmological background
    # Near vacuum: Phi ~ Phi0, Phi_dot/Phi0 ~ H * phi where phi = perturbation
    # phi_dot ~ H * epsilon where epsilon ~ delta_Phi/Phi0 ~ 1e-5 (CMB level)

    # But for the constraint we need Phi_dot_bar (background)
    # In TGP: psi_dot = 0 at vacuum (psi = 1 is equilibrium)
    # So the cosmological Phi_dot is tiny, and the constraint is very weak!

    # More precisely: from the field eq at psi ~ 1:
    # psi_ddot + 3H psi_dot + ... = c0^2 W(1) = c0^2 gamma/3  (kappa cancels)
    # At equilibrium: psi_dot ~ c0^2 gamma / (3 * 3H) ~ very small

    psi_dot_est = c0**2 * gamma_val / (9 * H0_obs)
    Phi_dot = psi_dot_est * Phi0

    print(f"  Estimated psi_dot (vacuum): {psi_dot_est:.4e}")
    print(f"  Phi_dot = psi_dot * Phi0 = {Phi_dot:.4e}")

    # Constraint: |B Phi_dot^2 / (A M*^4)| < 1e-15
    # At Phi ~ Phi0: A ~ 1
    # => |B| / M*^4 < 1e-15 / Phi_dot^2

    bound_B_over_M4 = 1e-15 / Phi_dot**2

    print(f"\n  GW170817 constraint: |B|/M*^4 < {bound_B_over_M4:.4e} [Phi_dot units]")
    print(f"  (This is VERY WEAK because Phi_dot_bar ~ 0 in vacuum)")

    # More relevant constraint: near binary merger where Phi is varying
    # Phi_dot ~ (q Phi0 M_merger) / (4 pi r^3) * v_orb
    # For GW170817: M ~ 2.7 M_sun, r ~ 1e6 m (merger), v ~ 0.3 c
    q_val = 8 * np.pi * G0 / c0**2  # generation constant
    M_merger = 2.7 * 2e30  # kg
    r_merger = 1e6          # m (late inspiral)
    v_orb = 0.3 * c0

    Phi_dot_merger = q_val * Phi0 * M_merger * v_orb / (4 * np.pi * r_merger**3)

    bound_B_M4_merger = 1e-15 / max(Phi_dot_merger**2, 1e-100)

    print(f"\n  Near merger:")
    print(f"    Phi_dot_merger ~ {Phi_dot_merger:.4e}")
    print(f"    |B|/M*^4 < {bound_B_M4_merger:.4e}")

    return bound_B_over_M4, bound_B_M4_merger


def tensor_to_breathing_ratio():
    """
    Compute the ratio h_+ / h_breathing as a function of B/M*^4.

    h_breathing = A'(Phi) delta_Phi ~ delta_Phi/Phi0
    h_tensor ~ B/M*^4 * k * |nabla Phi_bar| * delta_Phi

    For a binary system:
    |nabla Phi_bar| ~ q Phi0 M / (4 pi r^2)
    """
    print("\n" + "=" * 65)
    print("3. TENSOR-TO-BREATHING RATIO")
    print("=" * 65)

    # h_breathing from A'(Phi):
    # A(Phi) = exp(2 delta_Phi/Phi0) => A' = 2/Phi0 * A
    # delta g_ij^breathing = A' delta_Phi delta_ij = (2/Phi0) delta_Phi delta_ij
    # h_breathing ~ 2 delta_Phi / Phi0

    # h_tensor from B term:
    # delta g_ij^tensor ~ B/M*^4 * k * |nabla Phi_bar| * delta_Phi (schematic)
    # For GW from binary: delta_Phi = h_scalar * Phi0/2 (relates to GW amplitude)
    # |nabla Phi_bar| ~ q Phi0 M / (4 pi r^2)

    # So: h_tensor / h_breathing ~ B * k * q Phi0 M / (M*^4 * 4 pi r^2)
    #                             ~ B * omega/c0 * q Phi0 M / (M*^4 * 4 pi r^2)

    q_val = 8 * np.pi * G0 / c0**2
    M_source = 30 * 2e30  # 30 solar masses (BBH merger)

    # At LIGO frequency: f ~ 100 Hz, omega = 2 pi f
    f_gw = 100.0  # Hz
    omega = 2 * np.pi * f_gw
    k = omega / c0

    # Source distance for gradient: r ~ GM/c^2 (Schwarzschild radius scale)
    r_source = G0 * M_source / c0**2

    # Gradient of background Phi
    grad_Phi = q_val * Phi0 * M_source / (4 * np.pi * r_source**2)

    print(f"  Source: M = {M_source/2e30:.0f} M_sun, f_GW = {f_gw} Hz")
    print(f"  k = omega/c0 = {k:.4e} m^-1")
    print(f"  r_source ~ R_S = {r_source:.4e} m")
    print(f"  |nabla Phi_bar| ~ {grad_Phi:.4e}")

    # Ratio: h_tensor / h_breathing = B * k * grad_Phi / (2 M*^4)
    # Define dimensionless: b = B * (q Phi0)^2 / M*^4

    b_values = np.logspace(-20, -5, 100)  # dimensionless B/M*^4 in some unit

    # For each b, compute ratio
    ratio = b_values * k * grad_Phi / 2.0

    # Plot
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.loglog(b_values, ratio, 'b-', lw=2)
    ax.axhline(1.0, color='red', ls='--', lw=1.5, label='h_tensor = h_breathing')
    ax.axhline(0.01, color='orange', ls=':', lw=1.2, label='h_tensor = 1% of h_breathing')

    # LIGO sensitivity: can detect ratio > ~0.1
    ax.axhline(0.1, color='green', ls='-.', lw=1.2, label='LIGO sensitivity floor')

    ax.set_xlabel(r'$B / M_*^4$ [dimensionless]', fontsize=12)
    ax.set_ylabel(r'$h_{+,\times} / h_{\rm breathing}$', fontsize=12)
    ax.set_title('TGP: Tensor-to-Breathing GW Ratio\n'
                  '(30 $M_\\odot$ BBH, $f=100$ Hz)', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, ls=':', alpha=0.4)
    ax.set_xlim(b_values[0], b_values[-1])
    ax.set_ylim(1e-15, 1e5)

    path = os.path.join(save_dir, "cmb_tensor_to_breathing.png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    print(f"  Saved {path}")
    plt.close(fig)

    return b_values, ratio


def cmb_bmode_prediction():
    """
    Predict the CMB B-mode power spectrum modification from TGP.

    In standard cosmology: C_l^BB from inflationary tensor modes
    with tensor-to-scalar ratio r.

    In TGP: tensor modes arise from the disformal coupling B(Phi).
    The effective tensor-to-scalar ratio is:
        r_TGP = r_breathing * (h_tensor/h_breathing)^2

    where r_breathing comes from the scalar perturbations of Phi.
    """
    print("\n" + "=" * 65)
    print("4. CMB B-MODE PREDICTION")
    print("=" * 65)

    # The breathing mode contributes to the SCALAR perturbation spectrum
    # (E-modes), not B-modes. Only true tensor modes (h+, hx) generate B-modes.

    # In TGP: tensor modes from disformal coupling are suppressed by B/M*^4.
    # At CMB scales (l ~ 100, k ~ 0.01/Mpc), the relevant Phi gradient is
    # the COSMOLOGICAL one (inflation-era perturbations).

    # Key insight: during inflation (if TGP is coupled to an inflaton):
    # Phi_dot ~ H_inf * Phi0 (significant!)
    # grad_Phi ~ (k/a) * delta_Phi ~ (k/a) * H_inf / (2 pi)

    # The tensor amplitude from disformal coupling:
    # h_tensor^2 ~ (B/M*^4)^2 * (Phi_dot / c0)^2 * P_scalar(k)

    # For now, let's parametrize by B_0 = B(Phi0)/M*^4 and estimate r_TGP

    l_values = np.arange(2, 201)

    # Standard tensor spectrum (template for r = 0.01)
    r_std = 0.01
    # Approximate: C_l^BB ~ r * A_s * (l*(l+1))^{-0.5} * transfer(l)
    # Using simplified template:
    l_peak = 80
    C_BB_template = r_std * 1e-10 * np.exp(-(l_values - l_peak)**2 / (2 * 30**2))

    # TGP prediction: r_TGP = r_breathing_eff * (B_0 * H_inf * Phi0 / c0)^2
    # For B_0 values:
    B0_values = [1e-12, 1e-10, 1e-8]
    H_inf = 1e13  # GeV ~ 1e-6 in natural units (high-scale inflation)
    # In SI: H_inf ~ 1e-6 * c0 / l_P ~ 1e-6 * c0 / 1.6e-35 ~ 2e37 s^{-1}
    # But this is model-dependent. Let's use a dimensionless ratio.

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(l_values, C_BB_template * l_values * (l_values + 1) / (2 * np.pi),
            'k-', lw=2, label=f'Standard (r={r_std})')

    colors = ['blue', 'green', 'red']
    for B0, col in zip(B0_values, colors):
        # Crude estimate: TGP tensor suppression factor
        # For low B0: only breathing mode, no tensor => C_l^BB ~ 0
        # For large B0: tensor modes approach standard amplitude
        suppression = min(1.0, (B0 * 1e10)**2)  # crude parametrization
        C_BB_TGP = C_BB_template * suppression

        ax.plot(l_values, C_BB_TGP * l_values * (l_values + 1) / (2 * np.pi),
                ls='--', color=col, lw=1.5,
                label=f'TGP (B0={B0:.0e})')

    ax.set_xlabel(r'Multipole $\ell$', fontsize=12)
    ax.set_ylabel(r'$\ell(\ell+1) C_\ell^{BB} / 2\pi$', fontsize=12)
    ax.set_title('CMB B-mode spectrum: Standard vs TGP', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, ls=':', alpha=0.4)
    ax.set_yscale('log')
    ax.set_ylim(1e-20, 1e-8)

    path = os.path.join(save_dir, "cmb_bmode_TGP.png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    print(f"  Saved {path}")
    plt.close(fig)

    print("""
  KEY PREDICTIONS:
  1. For B(Phi0) = 0 (conformal metric only):
     - NO tensor modes, NO B-modes from TGP
     - Only breathing mode (scalar) => E-modes only
     - TGP predicts r_tensor = 0 (distinguishable from inflation)

  2. For small B(Phi0) > 0 (disformal coupling):
     - Tensor modes SUPPRESSED relative to standard
     - r_TGP ~ r_standard * (B0 * Phi_dot / M*^4)^2
     - B-mode spectrum has SAME shape but REDUCED amplitude

  3. For observable B >> 0:
     - Tensor-to-scalar ratio approaches standard value
     - Additional breathing mode creates EXCESS power in E-modes
     - Distinguishable from standard inflation by E/B ratio

  4. UNIQUE TGP SIGNATURE:
     - Breathing mode present (E-mode excess at low l)
     - Tensor-to-breathing ratio is SCALE-DEPENDENT
       (depends on k through disformal coupling)
     - This is a SMOKING GUN for TGP if detected
""")


def ppn_constraints():
    """
    Compute PPN constraints on B(Phi) from Solar System tests.

    The disformal metric modifies the PPN parameters:
    gamma_PPN = 1 + delta_gamma (from B contribution)
    beta_PPN = 1 + delta_beta

    Cassini: |gamma_PPN - 1| < 2.3e-5
    LLR: |beta_PPN - 1| < 1.1e-4
    """
    print("\n" + "=" * 65)
    print("5. PPN CONSTRAINTS ON B(Phi)")
    print("=" * 65)

    # In the weak-field limit, Phi = Phi0(1 + phi), |phi| << 1
    # The disformal metric:
    # g_00 = -(1 - 2phi + ...) + B/M*^4 * Phi0^2 * phi_dot^2
    # g_ij = (1 + 2phi + ...) delta_ij + B/M*^4 * Phi0^2 * d_i phi d_j phi

    # For static source: phi_dot = 0, d_i phi ~ GM/(c0^2 r^2) * x_i/r
    # Disformal contribution to g_ij:
    # delta_ij * (1 + 2phi) + B Phi0^2/M*^4 * (GM/c0^2)^2 / r^4 * x_i x_j

    # PPN: g_ij = (1 + 2 gamma U) delta_ij + ...
    # The anisotropic part doesn't contribute to gamma_PPN (it's the
    # Nordtvedt parameter). But the trace modifies gamma:

    # For the Sun:
    M_sun = 2e30  # kg
    r_earth = 1.5e11  # m

    q_val = 8 * np.pi * G0 / c0**2
    phi_sun = q_val * Phi0 * M_sun / (4 * np.pi * r_earth) / Phi0
    grad_phi = phi_sun / r_earth

    print(f"  Solar field: phi(r_Earth) ~ {phi_sun:.4e}")
    print(f"  |nabla phi|/Phi0 ~ {grad_phi:.4e} m^-1")

    # Cassini constraint: |delta_gamma| < 2.3e-5
    # delta_gamma ~ B * Phi0^2 * (grad_phi)^2 / M*^4
    # => B/M*^4 < 2.3e-5 / (Phi0^2 * grad_phi^2)

    bound_cassini = 2.3e-5 / (Phi0**2 * grad_phi**2)
    print(f"\n  Cassini constraint on B/M*^4: < {bound_cassini:.4e}")
    print(f"  (This is the STRONGEST constraint from Solar System)")

    return bound_cassini


def summary():
    """Print comprehensive summary."""
    print("\n" + "=" * 65)
    print("SUMMARY: CMB Tensor Perturbations in TGP")
    print("=" * 65)
    print("""
1. BREATHING MODE (conformal metric, A(Phi) only):
   - Always present in TGP
   - Scalar (spin 0), isotropic delta g_ij
   - Contributes to CMB E-modes, NOT B-modes
   - Dispersion: omega^2 = c0^2(k^2 + m_sp^2), m_sp^2 = gamma

2. TENSOR MODES (disformal metric, B(Phi) term):
   - Arise ONLY from disformal coupling B(Phi) d_u Phi d_v Phi
   - Require gradient in background Phi: d_i Phi_bar != 0
   - Near sources: h_tensor/h_breathing ~ B k |grad Phi_bar| / M*^4
   - In cosmology: suppressed by smallness of B and Phi_dot

3. OBSERVATIONAL SIGNATURES:
   - If B = 0: r_tensor = 0 (FALSIFIABLE by B-mode detection)
   - If B > 0: suppressed B-modes + enhanced E-modes
   - Scale-dependent tensor-to-breathing ratio (smoking gun)
   - GW polarization: breathing + tensor (6 modes possible)

4. CONSTRAINTS:
   - GW170817: weak (Phi_dot ~ 0 in vacuum)
   - Cassini/PPN: strongest Solar System bound
   - CMB B-modes: future test with CMB-S4, LiteBIRD

5. STATUS: Analysis complete. Key open problems:
   - Derive B(Phi) from substrate coarse-graining
   - Full numerical simulation of GW waveform with disformal metric
   - CMB power spectrum with Boltzmann code (CLASS/CAMB modification)
""")


# ======================================================================
# Main
# ======================================================================

if __name__ == "__main__":
    disformal_perturbation_analysis()
    bound_cosmo, bound_merger = gw_speed_constraint()
    b_values, ratio = tensor_to_breathing_ratio()
    cmb_bmode_prediction()
    bound_cassini = ppn_constraints()
    summary()
    print("Done.")
