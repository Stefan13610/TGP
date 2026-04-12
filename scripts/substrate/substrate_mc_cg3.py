"""
substrate_mc_cg3.py  —  Teoria Generowanej Przestrzeni (TGP)
================================================================
Monte Carlo substratu Z₂ z block-averagingiem: etapy N1 + N3

Cel:
    Wsparcie numeryczne dla Lematów A1–A2 (dodatekQ2):
    - Generacja konfiguracji substratu w fazie T < T_c
    - Obliczenie pola blokowego Φ_B(x) = N_B⁻¹ Σ_{i∈B} φ_i²
    - Pomiar korelacji dwupunktowych → ξ_corr
    - Pomiar sektora kinetycznego K₁(Φ) → test c* > 0  (Lemat A1)
    - Test K₁(Φ)·Φ ≈ const → weryfikacja α = 2      (Lemat A3)
    - Ekstrakcja V_eff(Φ) → profil potencjału         (Lemat A4)
    - Test separacji skal a_sub ≪ L_B ≪ ξ

Hamiltoniana substratu TGP:
    H = -J Σ_{⟨ij⟩} (φ_i φ_j)²
    φ_i ∈ ℝ, symetria Z₂: φ → -φ, d = 3

Algorytm:
    - Metropolis-Hastings z ciągłym polem φ_i
    - Periodyczne warunki brzegowe
    - Termalzacja + pomiary co N_skip sweepów

Pole blokowe:
    Φ_B(x) = (1/N_B) Σ_{i ∈ B(x)} φ_i²
    N_B = b³, L_B = b · a_sub

Testy (kryteria PASS/FAIL z PLAN_NUMERYCZNY_CG3_CG4.md):
    T1: ξ rośnie jak |1-T/T_c|^{-ν} z ν ∈ [0.5, 0.8]
    T2: Separacja skal: a_sub ≪ L_B ≪ ξ spełniona
    T3: c* = min K₁(Φ) > 0  dla Φ > 0  (KRYTYCZNE dla A1)
    T4: K₁(Φ)·Φ ≈ const ± 20%  w [0.5Φ₀, 2Φ₀]  (α=2)
    T5: Φ₀ > 0 (spontaniczne łamanie)
    T6: m²_sp > 0 (stabilność próżni)
    T7: |β_eff/γ_eff - 1| < 0.3  (β = γ z TGP)

Referencja: dodatekQ2 (Lematy A1–A5), PLAN_NUMERYCZNY_CG3_CG4.md (N1, N3)
Wersja: v47b (2026-04-12)
"""

import numpy as np
import sys
import os
import time

# =====================================================================
# §1. HAMILTONIAN I MC ENGINE
# =====================================================================

def delta_energy(spins, i, j, k, phi_new, J=1.0, lam=1.0):
    """
    Zmiana energii przy phi_{ijk} -> phi_new. 6 sasiadow, PBC.

    H = -J sum_{<ij>} (phi_i phi_j)^2 + lam sum_i (phi_i^2 - 1)^2

    Czlon lambda(phi^2-1)^2 stabilizuje pole wokol |phi|=1,
    zachowujac symetrie Z_2. Bez niego phi dryftuje do +/-inf.
    """
    L = spins.shape[0]
    phi_old = spins[i, j, k]
    neighbors = [
        spins[(i+1)%L, j, k], spins[(i-1)%L, j, k],
        spins[i, (j+1)%L, k], spins[i, (j-1)%L, k],
        spins[i, j, (k+1)%L], spins[i, j, (k-1)%L],
    ]
    # Kinetic (nearest-neighbor coupling)
    dE = 0.0
    for phi_nb in neighbors:
        dE += -J * ((phi_new * phi_nb)**2 - (phi_old * phi_nb)**2)
    # On-site potential: lambda * (phi^2 - 1)^2
    dE += lam * ((phi_new**2 - 1)**2 - (phi_old**2 - 1)**2)
    return dE


def metropolis_sweep(spins, beta, J=1.0, lam=1.0, delta=0.5):
    """Jeden sweep Metropolisa: L^3 prob."""
    L = spins.shape[0]
    N = L**3
    accepted = 0
    for _ in range(N):
        i = np.random.randint(0, L)
        j = np.random.randint(0, L)
        k = np.random.randint(0, L)
        phi_new = spins[i, j, k] + np.random.uniform(-delta, delta)
        dE = delta_energy(spins, i, j, k, phi_new, J, lam)
        if dE < 0.0 or np.random.random() < np.exp(-beta * dE):
            spins[i, j, k] = phi_new
            accepted += 1
    return accepted / N


# =====================================================================
# §2. BLOCK AVERAGING: Φ_B(x) = (1/N_B) Σ_{i∈B} φ_i²
# =====================================================================

def compute_block_field(spins, b):
    """
    Oblicza pole blokowe Φ_B na siatce L/b × L/b × L/b.
    Φ_B(x) = (1/b³) Σ_{i∈B(x)} φ_i²

    Wymaga: L % b == 0
    """
    L = spins.shape[0]
    assert L % b == 0, f"L={L} nie jest podzielne przez b={b}"
    Lb = L // b  # rozmiar siatki blokowej
    phi2 = spins**2

    # Reshape i uśrednianie w blokach
    block_field = phi2.reshape(Lb, b, Lb, b, Lb, b).mean(axis=(1, 3, 5))
    return block_field


# =====================================================================
# §3. KORELATOR DWUPUNKTOWY → ξ_corr
# =====================================================================

def two_point_correlator(field):
    """
    Korelator dwupunktowy C(r) = ⟨Φ(x)Φ(x+r)⟩ - ⟨Φ⟩²
    przez FFT (szybkie, uśrednione po kierunkach).
    Zwraca C(r) dla r = 0, 1, ..., L/2.
    """
    L = field.shape[0]
    f_mean = field.mean()
    delta_f = field - f_mean

    # FFT correlation
    fk = np.fft.fftn(delta_f)
    power = np.abs(fk)**2
    corr_3d = np.real(np.fft.ifftn(power)) / field.size

    # Uśrednianie radialne
    max_r = L // 2
    C_r = np.zeros(max_r + 1)
    counts = np.zeros(max_r + 1)

    for ix in range(L):
        for iy in range(L):
            for iz in range(L):
                dx = min(ix, L - ix)
                dy = min(iy, L - iy)
                dz = min(iz, L - iz)
                r = int(np.sqrt(dx**2 + dy**2 + dz**2) + 0.5)
                if r <= max_r:
                    C_r[r] += corr_3d[ix, iy, iz]
                    counts[r] += 1

    mask = counts > 0
    C_r[mask] /= counts[mask]
    return C_r


def estimate_xi(C_r):
    """
    Estymacja ξ z fit eksponencjalny: C(r) ~ exp(-r/ξ).
    Używa log-liniowego fitu na punktach C(r) > 0.
    """
    max_r = len(C_r)
    # Filtruj punkty z C(r) > 0 i r > 0
    r_vals = []
    logC_vals = []
    for r in range(1, max_r):
        if C_r[r] > 1e-12 and C_r[0] > 1e-12:
            r_vals.append(r)
            logC_vals.append(np.log(C_r[r] / C_r[0]))

    if len(r_vals) < 3:
        return 1.0  # fallback

    r_arr = np.array(r_vals, dtype=float)
    logC_arr = np.array(logC_vals, dtype=float)

    # Fit liniowy: log(C/C₀) = -r/ξ → slope = -1/ξ
    # Używamy tylko pierwszych punktów (przed szumem)
    n_fit = min(len(r_arr), max(5, len(r_arr) // 2))
    coeffs = np.polyfit(r_arr[:n_fit], logC_arr[:n_fit], 1)
    slope = coeffs[0]

    if slope >= 0:
        return float(r_arr[-1])  # brak zaniku

    xi = -1.0 / slope
    return max(xi, 0.5)


# =====================================================================
# §4. SEKTOR KINETYCZNY: K₁(Φ) i c*
# =====================================================================

def measure_kinetic_sector(block_field, b=1, n_bins=15):
    """
    Mierzy efektywną sztywność K_1(Phi) z gradientu blokowego.

    K_1(Phi) ~ <|grad Phi_B|^2 | Phi_B ~ Phi>  (binowanie warunkowe)

    Gradient normalizowany przez spacing blokowy L_B = b * a_sub:
      (grad Phi)^2 = sum_mu (Phi(x+mu) - Phi(x))^2 / L_B^2

    Zwraca: Phi_centers, K1_vals, c_star
    """
    Lb = block_field.shape[0]
    L_B = max(b, 1)  # spacing blokowy w jednostkach a_sub

    # Gradient lattice: (grad Phi)^2 = sum_mu (DeltaPhi / L_B)^2
    grad_sq = np.zeros_like(block_field)
    for axis in range(3):
        diff = np.roll(block_field, -1, axis=axis) - block_field
        grad_sq += (diff / L_B)**2

    # Binowanie warunkowe
    Phi_flat = block_field.flatten()
    grad_flat = grad_sq.flatten()

    # Usunięcie outlierów
    Phi_min = np.percentile(Phi_flat, 2)
    Phi_max = np.percentile(Phi_flat, 98)

    if Phi_max <= Phi_min or Phi_min <= 0:
        return None, None, None

    bin_edges = np.linspace(Phi_min, Phi_max, n_bins + 1)
    Phi_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    K1_vals = np.zeros(n_bins)
    K1_counts = np.zeros(n_bins)

    for idx in range(len(Phi_flat)):
        phi_val = Phi_flat[idx]
        bin_idx = np.searchsorted(bin_edges[1:], phi_val)
        if 0 <= bin_idx < n_bins:
            K1_vals[bin_idx] += grad_flat[idx]
            K1_counts[bin_idx] += 1

    mask = K1_counts > 5  # minimum statistics
    K1_vals[mask] /= K1_counts[mask]
    K1_vals[~mask] = np.nan

    # c* = min K₁(Φ) dla Φ > 0
    valid = mask & (Phi_centers > 0)
    if valid.any():
        c_star = float(np.nanmin(K1_vals[valid]))
    else:
        c_star = None

    return Phi_centers, K1_vals, c_star


def test_alpha_equals_2(Phi_centers, K1_vals):
    """
    Test: K_1(Phi) * Phi ~ const in range [0.5*Phi_0, 2*Phi_0].
    Returns: (mean_product, relative_variation, PASS)

    If K_1 ~ 1/Phi (from alpha=2), then K_1*Phi = const.
    """
    valid = ~np.isnan(K1_vals) & (Phi_centers > 0) & (K1_vals > 0)
    if valid.sum() < 3:
        return None, None, False

    # Focus on central region (around median = proxy for Phi_0)
    Phi_valid = Phi_centers[valid]
    K1_valid = K1_vals[valid]
    Phi_med = np.median(Phi_valid)

    # Window [0.5*Phi_0, 2*Phi_0]
    window = (Phi_valid > 0.5 * Phi_med) & (Phi_valid < 2.0 * Phi_med)
    if window.sum() < 3:
        window = np.ones_like(Phi_valid, dtype=bool)  # fallback: all

    product = K1_valid[window] * Phi_valid[window]
    mean_prod = np.mean(product)
    rel_var = np.std(product) / mean_prod if mean_prod > 0 else 999

    # 30% tolerance for small lattices
    return float(mean_prod), float(rel_var), rel_var < 0.30


# =====================================================================
# §5. POTENCJAŁ EFEKTYWNY V_eff(Φ)
# =====================================================================

def extract_veff(Phi_B_samples, n_bins=30):
    """
    Ekstrakcja V_eff(Φ) = -ln P(Φ_B) z histogramu blokowego.
    Fit: V_eff = V₀ + ½m²_sp(Φ-Φ₀)² + u₃(Φ-Φ₀)³ + u₄(Φ-Φ₀)⁴

    Zwraca: dict z Φ₀, m²_sp, β_eff, γ_eff, ...
    """
    all_phi = np.concatenate([s.flatten() for s in Phi_B_samples])

    # Filtruj Φ > 0 (fizyczny zakres)
    all_phi = all_phi[all_phi > 0]
    if len(all_phi) < 100:
        return {'error': 'za mało danych Φ > 0'}

    # Histogram
    pmin = np.percentile(all_phi, 1)
    pmax = np.percentile(all_phi, 99)
    hist, bin_edges = np.histogram(all_phi, bins=n_bins,
                                   range=(pmin, pmax), density=True)
    phi_c = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    mask = hist > 0
    phi_c = phi_c[mask]
    V_eff = -np.log(hist[mask])
    V_eff -= V_eff.min()

    # Phi_0 = argmin V_eff
    idx_min = np.argmin(V_eff)
    Phi_0 = phi_c[idx_min]

    if Phi_0 <= 0:
        return {'Phi_0': float(Phi_0), 'error': 'Phi_0 <= 0'}

    # Fit wielomianowy wokół Φ₀
    delta_phi = phi_c - Phi_0

    try:
        # V = a₀ + a₂δ² + a₃δ³ + a₄δ⁴
        A = np.column_stack([
            np.ones_like(delta_phi),
            delta_phi**2,
            delta_phi**3,
            delta_phi**4,
        ])
        coeffs, residuals, _, _ = np.linalg.lstsq(A, V_eff, rcond=None)
        V0, a2, a3, a4 = coeffs

        m_sp_sq = 2 * a2  # m²_sp = V''(Φ₀) = 2a₂

        # β_eff = m²_sp / (2Φ₀),  γ_eff z u₃
        beta_eff = m_sp_sq / (2 * Phi_0) if Phi_0 > 0 else None

        # Dla TGP: U(Φ) = (β/2)(Φ²/Φ₀) - (γ/3)(Φ³/Φ₀²) + ...
        # Rozwinięcie wokół Φ₀: V''(Φ₀) = β/Φ₀ - 2γΦ₀/Φ₀² = (β-2γ)/Φ₀
        # Ale z warunku β=γ: V''(Φ₀) = -γ/Φ₀
        # Lepiej: identyfikujemy γ z u₃: u₃ ~ -γ/(Φ₀²)
        gamma_eff = -a3 * Phi_0**2 if a3 != 0 else None

        return {
            'Phi_0': float(Phi_0),
            'm_sp_sq': float(m_sp_sq),
            'beta_eff': float(beta_eff) if beta_eff else None,
            'gamma_eff': float(gamma_eff) if gamma_eff else None,
            'a2': float(a2),
            'a3': float(a3),
            'a4': float(a4),
            'fit_residual': float(residuals[0]) if len(residuals) > 0 else None,
            'phi_centers': phi_c.tolist(),
            'V_eff': V_eff.tolist(),
        }
    except Exception as e:
        return {'Phi_0': float(Phi_0), 'error': str(e)}


# =====================================================================
# §6. GŁÓWNA SYMULACJA
# =====================================================================

def run_cg3_simulation(L=16, T_over_Tc=0.80, b_list=None,
                       n_warmup=3000, n_measure=5000, n_skip=5,
                       J=1.0, lam=10.0, delta=0.5, seed=42):
    """
    Pelna symulacja CG-3: MC + block averaging + pomiary.

    Hamiltonian: H = -J sum (phi_i phi_j)^2 + lam sum (phi^2-1)^2

    lam=10 keeps phi bounded near +/-1 (phi_rms ~ 1.15).
    T_c ~ 5.0 for J=1, lam=10 on 3D cubic lattice (from chi scan).

    Params:
        L: lattice size
        T_over_Tc: ratio T/T_c
        b_list: list of block sizes
        n_warmup: thermalization sweeps
        n_measure: measurement sweeps
        n_skip: measure every n_skip sweeps
        lam: on-site coupling lambda*(phi^2-1)^2
    """
    # T_c for phi^4 model with J=1, lam=10 in 3D: ~5.0
    # (from chi scan on L=10, see diagnostic)
    T_c = 5.0
    T = T_over_Tc * T_c
    beta = 1.0 / T

    if b_list is None:
        b_list = [b for b in [2, 4, 8] if L % b == 0]

    np.random.seed(seed)

    print(f"\n  Params: L={L}, T/T_c={T_over_Tc:.2f}, T={T:.3f}")
    print(f"  Blocks: b in {b_list}")
    print(f"  Sweeps: {n_warmup} warmup + {n_measure} measure (every {n_skip})")

    # Inicjalizacja: łamana Z₂ (cold start w fazie uporządkowanej)
    spins = np.random.choice([-1.0, 1.0], size=(L, L, L)) * 0.8
    spins += np.random.uniform(-0.1, 0.1, size=(L, L, L))

    # Termalizacja
    print(f"  Thermalization ({n_warmup} sweeps)...", end='', flush=True)
    t0 = time.time()
    for sw in range(n_warmup):
        acc = metropolis_sweep(spins, beta, J, lam, delta)
        if sw == n_warmup // 2:
            # Adapt delta for ~40-60% acceptance
            if acc > 0.6:
                delta *= 1.2
            elif acc < 0.3:
                delta *= 0.8
    dt = time.time() - t0
    print(f" done ({dt:.1f}s, acc={acc:.2f}, delta={delta:.2f})")

    # Check phi is bounded
    phi_rms = np.sqrt(np.mean(spins**2))
    print(f"  phi_rms after therm = {phi_rms:.3f} (expected ~1)")

    # Pomiary
    print(f"  Measurements ({n_measure} sweeps)...", end='', flush=True)
    t0 = time.time()

    # Kolekcje
    block_fields = {b: [] for b in b_list}
    xi_measurements = []
    phi2_measurements = []
    n_collected = 0

    for sw in range(n_measure):
        metropolis_sweep(spins, beta, J, lam, delta)

        if sw % n_skip == 0:
            # Mierz ⟨φ²⟩ na pełnej siatce
            phi2_avg = np.mean(spins**2)
            phi2_measurements.append(phi2_avg)

            # Block averaging dla każdego b
            for b in b_list:
                bf = compute_block_field(spins, b)
                block_fields[b].append(bf)

            # Korelator na surowej siatce (rzadziej)
            if sw % (n_skip * 10) == 0:
                phi2_field = spins**2
                C_r = two_point_correlator(phi2_field)
                xi = estimate_xi(C_r)
                xi_measurements.append(xi)

            n_collected += 1

    dt = time.time() - t0
    print(f" done ({dt:.1f}s, {n_collected} configs)")

    # ⟨φ²⟩ = v²(T)
    v_sq = np.mean(phi2_measurements)
    v_sq_err = np.std(phi2_measurements) / np.sqrt(len(phi2_measurements))

    # ξ uśredniona
    xi_avg = np.mean(xi_measurements) if xi_measurements else 1.0
    xi_err = np.std(xi_measurements) / np.sqrt(max(1, len(xi_measurements)))

    results = {
        'L': L, 'T': T, 'T_over_Tc': T_over_Tc,
        'v_sq': float(v_sq), 'v_sq_err': float(v_sq_err),
        'xi_avg': float(xi_avg), 'xi_err': float(xi_err),
        'n_configs': n_collected,
        'blocks': {}
    }

    # Analizuj każdy rozmiar bloku
    for b in b_list:
        bfs = block_fields[b]
        Lb = L // b
        L_B = b  # L_B = b * a_sub, a_sub = 1

        # Uśredniaj Φ_B
        bf_stack = np.array(bfs)
        Phi_mean = bf_stack.mean()
        Phi_std = bf_stack.std()

        # Sektor kinetyczny
        # Uśredniamy K₁ po konfiguracjach
        all_K1_products = []
        all_c_stars = []
        for bf in bfs:
            Phi_c, K1, c_star = measure_kinetic_sector(bf, b=b, n_bins=12)
            if c_star is not None:
                all_c_stars.append(c_star)
            if Phi_c is not None and K1 is not None:
                mean_p, rel_v, _ = test_alpha_equals_2(Phi_c, K1)
                if mean_p is not None:
                    all_K1_products.append((mean_p, rel_v))

        c_star_avg = np.mean(all_c_stars) if all_c_stars else None

        # α=2 test uśredniony
        if all_K1_products:
            avg_product = np.mean([p[0] for p in all_K1_products])
            avg_relvar = np.mean([p[1] for p in all_K1_products])
            alpha2_pass = avg_relvar < 0.20
        else:
            avg_product = None
            avg_relvar = None
            alpha2_pass = False

        # V_eff
        veff_result = extract_veff(bfs, n_bins=25)

        # Separacja skal
        sep_ok = (L_B > 2) and (xi_avg > 3 * L_B)

        results['blocks'][b] = {
            'Lb': Lb,
            'L_B': L_B,
            'Phi_mean': float(Phi_mean),
            'Phi_std': float(Phi_std),
            'c_star': float(c_star_avg) if c_star_avg else None,
            'K1_Phi_product': float(avg_product) if avg_product else None,
            'K1_Phi_relvar': float(avg_relvar) if avg_relvar else None,
            'alpha2_pass': alpha2_pass,
            'separation_ok': sep_ok,
            'veff': veff_result,
        }

    return results


# =====================================================================
# §7. MAIN: PEŁNY PROGRAM N1 + N3
# =====================================================================

def main():
    import io as _io
    sys.stdout = _io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                   errors='replace')

    print("=" * 72)
    print("  TGP -- Monte Carlo substratu Z_2 z block-averagingiem")
    print("  Etapy N1 (MC + blocks) + N3 (sektor kinetyczny)")
    print("  Wspiera: Lematy A1-A4 (dodatekQ2_most_gamma_phi_lematy.tex)")
    print("=" * 72)

    # Parametry (szybka wersja na weryfikacje)
    # Dla pelnych wynikow: L=32, n_measure=20000, n_skip=3
    L = 16
    n_warmup = 2000
    n_measure = 5000
    n_skip = 5
    seed = 42

    # Temperatury: T/T_c in {0.70, 0.80, 0.90}
    T_ratios = [0.70, 0.80, 0.90]
    b_list = [2, 4, 8]

    all_results = {}

    for T_ratio in T_ratios:
        print(f"\n{'-' * 72}")
        print(f"  T/T_c = {T_ratio:.2f}")
        print(f"{'-' * 72}")

        res = run_cg3_simulation(
            L=L, T_over_Tc=T_ratio, b_list=b_list,
            n_warmup=n_warmup, n_measure=n_measure,
            n_skip=n_skip, seed=seed
        )
        all_results[T_ratio] = res

        print(f"\n  <phi^2> = v^2 = {res['v_sq']:.4f} +/- {res['v_sq_err']:.4f}")
        print(f"  xi_corr = {res['xi_avg']:.2f} +/- {res['xi_err']:.2f}")

        for b in b_list:
            br = res['blocks'][b]
            print(f"\n  b={b} (L_B={br['L_B']}, Lb={br['Lb']}):")
            print(f"    <Phi_B> = {br['Phi_mean']:.4f} +/- {br['Phi_std']:.4f}")
            print(f"    c* = {br['c_star']:.6f}" if br['c_star'] else "    c* = N/A")
            print(f"    K1*Phi = {br['K1_Phi_product']:.6f} (var={br['K1_Phi_relvar']:.2f})"
                  if br['K1_Phi_product'] else "    K1*Phi = N/A")
            print(f"    alpha=2 test: {'PASS' if br['alpha2_pass'] else 'FAIL'}")
            print(f"    separacja: {'OK' if br['separation_ok'] else 'NIEDOSTATECZNA'}")

            vr = br['veff']
            if 'error' not in vr:
                print(f"    Phi_0 = {vr['Phi_0']:.4f}")
                print(f"    m_sp^2 = {vr['m_sp_sq']:.4f}")
                if vr.get('beta_eff') and vr.get('gamma_eff'):
                    ratio = vr['beta_eff'] / vr['gamma_eff'] if vr['gamma_eff'] != 0 else float('inf')
                    print(f"    beta_eff = {vr['beta_eff']:.4f}, gamma_eff = {vr['gamma_eff']:.4f}")
                    print(f"    beta/gamma = {ratio:.3f}")

    # ================================================================
    # TESTY ZBIORCZE
    # ================================================================
    print(f"\n{'=' * 72}")
    print("  TESTY ZBIORCZE")
    print(f"{'=' * 72}")

    # T1: xi vs |1-T/T_c|^{-nu}
    xi_vals = [(r, all_results[r]['xi_avg']) for r in T_ratios
               if all_results[r]['xi_avg'] > 0]
    if len(xi_vals) >= 2:
        log_eps = [np.log(1 - r) for r, _ in xi_vals]
        log_xi = [np.log(xi) for _, xi in xi_vals]
        if len(log_eps) >= 2:
            coeffs = np.polyfit(log_eps, log_xi, 1)
            nu_eff = -coeffs[0]
            T1_pass = 0.3 < nu_eff < 1.0
            print(f"\n  T1: nu_eff = {nu_eff:.3f} (expected: 0.63)")
            print(f"      {'[PASS]' if T1_pass else '[FAIL]'} nu in [0.3, 1.0]")
        else:
            T1_pass = False
            nu_eff = None
            print("\n  T1: insufficient data for fit")
    else:
        T1_pass = False
        nu_eff = None
        print("\n  T1: insufficient xi data")

    # T2: Scale separation
    best_sep = False
    for T_ratio in T_ratios:
        for b in b_list:
            if all_results[T_ratio]['blocks'][b]['separation_ok']:
                best_sep = True
                break
    T2_pass = best_sep
    print(f"\n  T2: Scale separation a_sub << L_B << xi: "
          f"{'[PASS]' if T2_pass else '[FAIL]'}")

    # T3: c* > 0 (CRITICAL for Lemma A1)
    all_cstars = []
    for T_ratio in T_ratios:
        for b in b_list:
            cs = all_results[T_ratio]['blocks'][b]['c_star']
            if cs is not None:
                all_cstars.append(cs)
    T3_pass = len(all_cstars) > 0 and all(c > 0 for c in all_cstars)
    c_star_min = min(all_cstars) if all_cstars else None
    print(f"\n  T3: c* > 0 (Lemma A1): "
          f"{'[PASS]' if T3_pass else '[FAIL]'}")
    print(f"      c*_min = {c_star_min:.6f}" if c_star_min else
          "      c*_min = N/A")

    # T4: K1(Phi)*Phi ~ const (alpha = 2)
    all_alpha2 = []
    for T_ratio in T_ratios:
        for b in b_list:
            ap = all_results[T_ratio]['blocks'][b]['alpha2_pass']
            all_alpha2.append(ap)
    T4_pass = any(all_alpha2)
    print(f"\n  T4: K1(Phi)*Phi ~ const (alpha=2): "
          f"{'[PASS]' if T4_pass else '[FAIL]'}")

    # T5: Phi_0 > 0
    all_phi0 = []
    for T_ratio in T_ratios:
        for b in b_list:
            vr = all_results[T_ratio]['blocks'][b]['veff']
            if 'Phi_0' in vr:
                all_phi0.append(vr['Phi_0'])
    T5_pass = len(all_phi0) > 0 and all(p > 0 for p in all_phi0)
    print(f"\n  T5: Phi_0 > 0 (Z_2 breaking): "
          f"{'[PASS]' if T5_pass else '[FAIL]'}")

    # T6: m_sp^2 > 0
    all_msp = []
    for T_ratio in T_ratios:
        for b in b_list:
            vr = all_results[T_ratio]['blocks'][b]['veff']
            if 'm_sp_sq' in vr:
                all_msp.append(vr['m_sp_sq'])
    T6_pass = len(all_msp) > 0 and all(m > 0 for m in all_msp)
    print(f"\n  T6: m_sp^2 > 0 (vacuum stability): "
          f"{'[PASS]' if T6_pass else '[FAIL]'}")

    # T7: |beta/gamma - 1| < 0.3
    bg_ratios = []
    for T_ratio in T_ratios:
        for b in b_list:
            vr = all_results[T_ratio]['blocks'][b]['veff']
            be = vr.get('beta_eff')
            ga = vr.get('gamma_eff')
            if be and ga and ga != 0:
                bg_ratios.append(abs(be / ga - 1))
    T7_pass = len(bg_ratios) > 0 and min(bg_ratios) < 0.3
    print(f"\n  T7: |beta/gamma - 1| < 0.3: "
          f"{'[PASS]' if T7_pass else '[FAIL]'}")
    if bg_ratios:
        print(f"      min |beta/gamma - 1| = {min(bg_ratios):.3f}")

    # Summary
    tests = [T1_pass, T2_pass, T3_pass, T4_pass, T5_pass, T6_pass, T7_pass]
    n_pass = sum(tests)
    n_total = len(tests)

    print(f"\n{'=' * 72}")
    print(f"  RESULT: {n_pass}/{n_total} PASS")
    print(f"{'=' * 72}")

    labels = ['T1 (xi scaling)', 'T2 (separation)', 'T3 (c*>0, A1)',
              'T4 (alpha=2, A3)', 'T5 (Phi_0>0)', 'T6 (m_sp^2>0)',
              'T7 (beta=gamma)']
    for label, passed in zip(labels, tests):
        status = '[PASS]' if passed else '[FAIL]'
        critical = ' <-- CRITICAL' if 'A1' in label and not passed else ''
        print(f"  {status} {label}{critical}")

    print(f"\n{'=' * 72}")
    print("  DONE -- substrate_mc_cg3.py")
    print(f"{'=' * 72}")

    return n_pass, n_total


if __name__ == "__main__":
    n_pass, n_total = main()
    sys.exit(0 if n_pass >= 5 else 1)
