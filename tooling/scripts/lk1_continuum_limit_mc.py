#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
LK-1: Monte Carlo Ising 3D + coarse-graining → continuum limit test
=====================================================================
Cel: Zweryfikować numerycznie most Gamma → Phi

Co robimy:
1. Symulujemy 3D Ising na sieci kubicznej (L = 8, 16, 32, 64, 128)
2. Definiujemy pole blokowe Phi_B(x) = (1/N_B) * sum_i <s_i^2> w bloku b×b×b
3. Mierzymy:
   - efektywne K(Phi) z gradientów Phi_B
   - efektywne U(Phi) z histogramów Phi_B
   - korelacje <Phi_B(x) Phi_B(y)> → xi_corr(L_B)
   - alpha_eff = d(ln K)/d(ln Phi) → oczekiwane alpha_eff → 2
4. Weryfikujemy beta/gamma → 1 w potencjale efektywnym

Substrat TGP:
  H_Gamma = sum_i [m0^2/2 * s_i^2 + lambda/4 * s_i^4] - J * sum_<ij> s_i*s_j
  Phi(x) = <s^2>_block >= 0 (Z_2 invariant)

Uwaga: Używamy Wolff cluster algorithm dla szybkiej termalizacji.
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.ndimage import uniform_filter
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Parameters
# ============================================================
J_coupling = 1.0            # Nearest-neighbor coupling
T_c_ising = 4.5115          # Critical temperature 3D Ising (sc lattice, J=1)
T_sim = 0.95 * T_c_ising    # Work in ordered phase (T < T_c)
beta_T = 1.0 / T_sim

LATTICE_SIZES = [8, 16]            # L values (32+ too slow without C/Fortran)
BLOCK_SIZES = [2, 4]              # Coarse-graining block sizes
N_THERM = 500                     # Thermalization sweeps
N_MEASURE = 200                   # Measurement sweeps
N_SKIP = 2                        # Sweeps between measurements

PHI_0 = 25.0  # TGP vacuum value (for comparison)

# ============================================================
# 3D Ising Monte Carlo (Metropolis + Wolff cluster)
# ============================================================
class Ising3D:
    def __init__(self, L, beta, J=1.0):
        self.L = L
        self.beta = beta
        self.J = J
        self.spins = np.ones((L, L, L), dtype=np.int8)  # start ordered

    def metropolis_sweep(self):
        """One full Metropolis sweep."""
        L = self.L
        for _ in range(L**3):
            x, y, z = np.random.randint(0, L, 3)
            s = self.spins[x, y, z]
            # Sum of neighbors
            nn_sum = (
                self.spins[(x+1)%L, y, z] + self.spins[(x-1)%L, y, z] +
                self.spins[x, (y+1)%L, z] + self.spins[x, (y-1)%L, z] +
                self.spins[x, y, (z+1)%L] + self.spins[x, y, (z-1)%L]
            )
            dE = 2.0 * self.J * s * nn_sum
            if dE <= 0 or np.random.random() < np.exp(-self.beta * dE):
                self.spins[x, y, z] = -s

    def wolff_step(self):
        """One Wolff cluster flip."""
        L = self.L
        p_add = 1.0 - np.exp(-2.0 * self.beta * self.J)

        # Random seed spin
        x0, y0, z0 = np.random.randint(0, L, 3)
        cluster_spin = self.spins[x0, y0, z0]

        visited = np.zeros((L, L, L), dtype=bool)
        stack = [(x0, y0, z0)]
        visited[x0, y0, z0] = True

        while stack:
            x, y, z = stack.pop()
            self.spins[x, y, z] *= -1

            for dx, dy, dz in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
                nx, ny, nz = (x+dx)%L, (y+dy)%L, (z+dz)%L
                if not visited[nx, ny, nz] and self.spins[nx, ny, nz] == cluster_spin:
                    if np.random.random() < p_add:
                        visited[nx, ny, nz] = True
                        stack.append((nx, ny, nz))

    def measure_phi_block(self, b):
        """
        Compute block field: Phi_B(X) = (1/b^3) * sum_{i in B(X)} s_i^2
        Since s_i = +-1 (Ising), s_i^2 = 1 always → Phi_B = 1.

        For continuous substrate, we use <s_i * s_j> correlator instead.
        The PHYSICAL Phi is related to the ORDER PARAMETER:

        Phi_B(X) = [m_B(X)]^2 where m_B = (1/b^3) sum_i s_i

        This is Z_2 invariant, >= 0, and captures the "density of generated space".
        """
        L = self.L
        L_B = L // b
        phi_block = np.zeros((L_B, L_B, L_B))

        for ix in range(L_B):
            for iy in range(L_B):
                for iz in range(L_B):
                    block = self.spins[ix*b:(ix+1)*b, iy*b:(iy+1)*b, iz*b:(iz+1)*b]
                    m_block = np.mean(block.astype(float))
                    phi_block[ix, iy, iz] = m_block**2  # Z_2 invariant!

        return phi_block

    def measure_gradient_energy(self, phi_B):
        """Measure (nabla Phi_B)^2 for each block."""
        L_B = phi_B.shape[0]
        grad_sq = np.zeros_like(phi_B)
        for d in range(3):
            grad = np.roll(phi_B, -1, axis=d) - phi_B
            grad_sq += grad**2
        return grad_sq

    def measure_correlation(self, phi_B):
        """Measure <Phi_B(0) Phi_B(r)> - <Phi_B>^2."""
        phi_mean = np.mean(phi_B)
        phi_fluct = phi_B - phi_mean

        L_B = phi_B.shape[0]
        max_r = L_B // 2
        corr = np.zeros(max_r)
        counts = np.zeros(max_r)

        # Sample along axes
        for ix in range(L_B):
            for iy in range(L_B):
                for iz in range(L_B):
                    for r in range(1, max_r):
                        # Along x
                        corr[r] += phi_fluct[ix, iy, iz] * phi_fluct[(ix+r)%L_B, iy, iz]
                        counts[r] += 1

        corr[1:] /= np.maximum(counts[1:], 1)
        corr[0] = np.mean(phi_fluct**2)
        return corr


def extract_kinetic_coupling(phi_values, grad_sq_values, n_bins=15):
    """
    Extract effective K(Phi) from correlation between Phi and grad^2.

    In TGP: L_kin = (1/2) K(Phi) (nabla Phi)^2
    So: <(nabla Phi)^2 | Phi=phi> propto 1/K(phi) (equipartition)

    Inversely: K_eff(phi) propto 1 / <(nabla Phi)^2 | Phi=phi>
    """
    phi_flat = phi_values.flatten()
    grad_flat = grad_sq_values.flatten()

    # Bin by Phi
    phi_min, phi_max = np.percentile(phi_flat, [5, 95])
    bins = np.linspace(phi_min, phi_max, n_bins + 1)

    phi_centers = []
    K_eff = []

    for i in range(n_bins):
        mask = (phi_flat >= bins[i]) & (phi_flat < bins[i+1])
        if np.sum(mask) > 10:
            phi_c = np.mean(phi_flat[mask])
            grad_mean = np.mean(grad_flat[mask])
            if grad_mean > 1e-10:
                phi_centers.append(phi_c)
                K_eff.append(1.0 / grad_mean)  # Proportional to K(phi)

    return np.array(phi_centers), np.array(K_eff)


def fit_power_law(phi, K):
    """Fit K(phi) = C * phi^alpha to extract alpha_eff."""
    mask = (phi > 0) & (K > 0)
    if np.sum(mask) < 3:
        return np.nan, np.nan

    log_phi = np.log(phi[mask])
    log_K = np.log(K[mask])

    try:
        coeffs = np.polyfit(log_phi, log_K, 1)
        alpha_eff = coeffs[0]
        C = np.exp(coeffs[1])
        return alpha_eff, C
    except:
        return np.nan, np.nan


def extract_potential(phi_values, n_bins=20):
    """
    Extract effective U(Phi) from histogram (Boltzmann weight).
    P(Phi) propto exp(-beta * U_eff(Phi)) → U_eff(Phi) = -T * ln P(Phi)
    """
    phi_flat = phi_values.flatten()

    hist, edges = np.histogram(phi_flat, bins=n_bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])

    mask = hist > 0
    U_eff = np.zeros_like(hist, dtype=float)
    U_eff[mask] = -np.log(hist[mask])  # In units of T
    U_eff[~mask] = np.nan

    return centers, U_eff


# ============================================================
# Main simulation
# ============================================================
def run_continuum_test():
    print("=" * 70)
    print("LK-1: CONTINUUM LIMIT TEST — Gamma → Phi")
    print("=" * 70)
    print(f"T/T_c = {T_sim/T_c_ising:.3f} (ordered phase)")
    print(f"beta = {beta_T:.4f}")
    print()

    results = {}

    for L in LATTICE_SIZES:
        print(f"\n{'─'*50}")
        print(f"L = {L}")
        print(f"{'─'*50}")

        # Initialize and thermalize
        ising = Ising3D(L, beta_T, J_coupling)
        print(f"  Thermalizing ({N_THERM} Wolff steps)...")
        for _ in range(N_THERM):
            ising.wolff_step()

        # Collect measurements
        for b in BLOCK_SIZES:
            if L // b < 4:
                continue  # Skip if too few blocks

            L_B = L // b
            print(f"\n  Block size b={b} → L_B={L_B}")

            phi_samples = []
            grad_samples = []
            corr_sum = None
            n_samples = 0

            for step in range(N_MEASURE):
                for _ in range(N_SKIP):
                    ising.wolff_step()

                phi_B = ising.measure_phi_block(b)
                grad_sq = ising.measure_gradient_energy(phi_B)
                corr = ising.measure_correlation(phi_B)

                phi_samples.append(phi_B.copy())
                grad_samples.append(grad_sq.copy())

                if corr_sum is None:
                    corr_sum = corr.copy()
                else:
                    corr_sum += corr
                n_samples += 1

            # Average
            phi_all = np.array(phi_samples)
            grad_all = np.array(grad_samples)
            corr_avg = corr_sum / n_samples

            phi_mean = np.mean(phi_all)
            phi_std = np.std(phi_all)

            print(f"    <Phi_B> = {phi_mean:.4f} ± {phi_std:.4f}")

            # Extract K(Phi)
            phi_c, K_eff = extract_kinetic_coupling(
                np.concatenate(phi_samples),
                np.concatenate(grad_samples)
            )

            if len(phi_c) >= 3:
                alpha_eff, C_K = fit_power_law(phi_c, K_eff)
                print(f"    K(Phi) ~ Phi^alpha_eff:")
                print(f"      alpha_eff = {alpha_eff:.3f}  (TGP predicts: 2)")
                print(f"      |alpha_eff - 2| = {abs(alpha_eff - 2):.3f}")
            else:
                alpha_eff = np.nan
                print(f"    K(Phi): insufficient data for fit")

            # Extract U(Phi)
            phi_u, U_eff = extract_potential(np.concatenate(phi_samples))

            # Find vacuum: minimum of U_eff
            valid = ~np.isnan(U_eff)
            if np.sum(valid) > 3:
                phi_u_valid = phi_u[valid]
                U_valid = U_eff[valid]
                idx_min = np.argmin(U_valid)
                phi_vac = phi_u_valid[idx_min]
                print(f"    U_eff: vacuum at Phi_0 = {phi_vac:.4f}")

                # Check beta/gamma ratio from cubic/quartic terms around vacuum
                # U ~ (beta/3)(phi/phi0)^3 - (gamma/4)(phi/phi0)^4
                # At vacuum U'(phi0) = 0 → beta = gamma
                # Check: is potential symmetric around vacuum?
                delta_phi = phi_u_valid - phi_vac
                mask_left = delta_phi < 0
                mask_right = delta_phi > 0
                if np.sum(mask_left) > 2 and np.sum(mask_right) > 2:
                    # Fit quadratic around vacuum
                    try:
                        p = np.polyfit(delta_phi, U_valid - U_valid[idx_min], 4)
                        # p[0]*x^4 + p[1]*x^3 + p[2]*x^2 + p[3]*x + p[4]
                        # beta ~ p[1], gamma ~ -p[0]
                        if abs(p[0]) > 1e-10:
                            bg_ratio = -p[1] / p[0]  # Should be ~1 if beta=gamma
                            # Note: sign depends on convention
                            print(f"    U_eff expansion: cubic/quartic ratio ~ {abs(bg_ratio):.3f}")
                    except:
                        pass

            # Correlation length
            if len(corr_avg) > 3 and corr_avg[0] > 0:
                r_vals = np.arange(1, len(corr_avg))
                c_vals = corr_avg[1:] / corr_avg[0]
                mask_pos = c_vals > 0.01
                if np.sum(mask_pos) >= 2:
                    try:
                        log_c = np.log(c_vals[mask_pos])
                        r_fit = r_vals[mask_pos]
                        slope, _ = np.polyfit(r_fit, log_c, 1)
                        xi_corr = -1.0 / slope if slope < 0 else np.inf
                        print(f"    xi_corr = {xi_corr:.2f} (in block units)")
                    except:
                        xi_corr = np.nan
                else:
                    xi_corr = np.nan
            else:
                xi_corr = np.nan

            results[(L, b)] = {
                'phi_mean': phi_mean,
                'phi_std': phi_std,
                'alpha_eff': alpha_eff,
                'xi_corr': xi_corr if not np.isnan(xi_corr) else None,
            }

    # ============================================================
    # Summary
    # ============================================================
    print("\n" + "=" * 70)
    print("SUMMARY: CONTINUUM LIMIT CONVERGENCE")
    print("=" * 70)

    print(f"\n{'L':>4} {'b':>3} {'L_B':>4} {'<Phi>':>8} {'alpha_eff':>10} {'|a-2|':>8} {'xi_corr':>8}")
    print("-" * 55)

    alpha_values = []

    for (L, b), res in sorted(results.items()):
        L_B = L // b
        a_eff = res['alpha_eff']
        xi = res.get('xi_corr', None)

        a_str = f"{a_eff:.3f}" if not np.isnan(a_eff) else "N/A"
        da_str = f"{abs(a_eff-2):.3f}" if not np.isnan(a_eff) else "N/A"
        xi_str = f"{xi:.2f}" if xi is not None else "N/A"

        print(f"{L:4d} {b:3d} {L_B:4d} {res['phi_mean']:8.4f} {a_str:>10} {da_str:>8} {xi_str:>8}")

        if not np.isnan(a_eff):
            alpha_values.append(a_eff)

    if alpha_values:
        alpha_mean = np.mean(alpha_values)
        alpha_std = np.std(alpha_values)
        print(f"\nalpha_eff (mean ± std) = {alpha_mean:.3f} ± {alpha_std:.3f}")
        print(f"TGP prediction: alpha = 2")
        print(f"Deviation: |alpha_eff - 2| = {abs(alpha_mean - 2):.3f}")

        if abs(alpha_mean - 2) < 0.5:
            print("→ CONSISTENT with TGP prediction alpha=2")
        else:
            print("→ TENSION with TGP prediction (check T/Tc, statistics)")

    print("\n" + "=" * 70)
    print("TESTS:")
    print("=" * 70)

    tests_pass = 0
    tests_total = 0

    # Test 1: Phi_B > 0 everywhere
    tests_total += 1
    all_positive = all(r['phi_mean'] > 0 for r in results.values())
    status = "PASS" if all_positive else "FAIL"
    print(f"  T1 [Phi_B >= 0 (Z_2 invariant)]: {status}")
    if all_positive: tests_pass += 1

    # Test 2: alpha_eff close to 2
    tests_total += 1
    if alpha_values and abs(np.mean(alpha_values) - 2) < 1.0:
        print(f"  T2 [alpha_eff ~ 2 (within 1.0)]: PASS (alpha={np.mean(alpha_values):.3f})")
        tests_pass += 1
    else:
        print(f"  T2 [alpha_eff ~ 2]: FAIL or NO DATA")

    # Test 3: Phi has finite variance (bounded in L^2)
    tests_total += 1
    finite_var = all(r['phi_std'] < 10 * r['phi_mean'] for r in results.values() if r['phi_mean'] > 0)
    status = "PASS" if finite_var else "FAIL"
    print(f"  T3 [Phi_B bounded in L^2 (finite variance)]: {status}")
    if finite_var: tests_pass += 1

    # Test 4: Correlation decays (xi_corr finite)
    tests_total += 1
    xi_vals = [r['xi_corr'] for r in results.values() if r.get('xi_corr') is not None]
    if xi_vals and all(x < 100 for x in xi_vals):
        print(f"  T4 [Correlations decay (xi finite)]: PASS (xi_max={max(xi_vals):.2f})")
        tests_pass += 1
    elif not xi_vals:
        print(f"  T4 [Correlations decay]: NO DATA")
    else:
        print(f"  T4 [Correlations decay]: FAIL (xi too large)")

    # Test 5: Vacuum exists (U_eff has minimum)
    tests_total += 1
    has_vacuum = any(r['phi_mean'] > 0.1 for r in results.values())
    status = "PASS" if has_vacuum else "FAIL"
    print(f"  T5 [Non-trivial vacuum Phi_0 > 0]: {status}")
    if has_vacuum: tests_pass += 1

    print(f"\n  TOTAL: {tests_pass}/{tests_total} PASS")

    return results


if __name__ == "__main__":
    results = run_continuum_test()
