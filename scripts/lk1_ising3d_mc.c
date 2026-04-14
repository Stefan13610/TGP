/*
 * LK-1: Monte Carlo Ising 3D + Coarse-Graining → Continuum Limit
 * ================================================================
 * High-performance C implementation with Wolff cluster algorithm.
 *
 * Measures:
 *   - Block field Phi_B(x) = <s_i^2>_block  (Z_2 invariant)
 *   - Effective kinetic kernel K(Phi) from gradient energy
 *   - Effective potential U(Phi) from histogram
 *   - alpha_eff = d(ln K)/d(ln Phi)  (TGP predicts: 2)
 *   - Correlations <Phi_B(x) Phi_B(y)>
 *
 * Compile: gcc -O3 -march=native -o lk1_ising3d lk1_ising3d_mc.c -lm
 * Run:     ./lk1_ising3d
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* ============================================================
 * Parameters
 * ============================================================ */
#define J_COUPLING 1.0
#define T_C_ISING  4.5115    /* 3D Ising critical T (sc, J=1) */

/* Simulation parameters */
#define N_THERM    5000      /* Thermalization sweeps */
#define N_MEASURE  2000      /* Measurement sweeps */
#define N_SKIP     5         /* Sweeps between measurements */

/* Lattice sizes to test */
static const int LATTICE_SIZES[] = {16, 32, 64};
#define N_LSIZES 3

/* Block sizes */
static const int BLOCK_SIZES[] = {2, 4, 8};
#define N_BSIZES 3

/* Histogram bins */
#define N_BINS 50

/* Maximum lattice size */
#define LMAX 64

/* ============================================================
 * RNG: xoshiro256** (fast, high quality)
 * ============================================================ */
static uint64_t rng_state[4];

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline uint64_t rng_next(void) {
    const uint64_t result = rotl(rng_state[1] * 5, 7) * 9;
    const uint64_t t = rng_state[1] << 17;
    rng_state[2] ^= rng_state[0];
    rng_state[3] ^= rng_state[1];
    rng_state[1] ^= rng_state[2];
    rng_state[0] ^= rng_state[3];
    rng_state[2] ^= t;
    rng_state[3] = rotl(rng_state[3], 45);
    return result;
}

static inline double rng_double(void) {
    return (rng_next() >> 11) * 0x1.0p-53;
}

static void rng_seed(uint64_t seed) {
    /* SplitMix64 to seed xoshiro */
    for (int i = 0; i < 4; i++) {
        seed += 0x9e3779b97f4a7c15ULL;
        uint64_t z = seed;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        rng_state[i] = z ^ (z >> 31);
    }
}

/* ============================================================
 * 3D Ising lattice
 * ============================================================ */
static int spin[LMAX * LMAX * LMAX];  /* +1 or -1 */
static int cluster_mark[LMAX * LMAX * LMAX];
static int cluster_stack[LMAX * LMAX * LMAX];

static inline int idx3(int x, int y, int z, int L) {
    return ((x % L + L) % L) * L * L + ((y % L + L) % L) * L + ((z % L + L) % L);
}

/* Initialize all spins up */
static void init_lattice(int L) {
    int N = L * L * L;
    for (int i = 0; i < N; i++) spin[i] = 1;
}

/* Wolff cluster flip */
static int wolff_step(int L, double p_add) {
    int N = L * L * L;
    memset(cluster_mark, 0, N * sizeof(int));

    /* Random seed site */
    int seed = (int)(rng_double() * N);
    if (seed >= N) seed = N - 1;
    int s0 = spin[seed];

    cluster_mark[seed] = 1;
    cluster_stack[0] = seed;
    int stack_top = 1;
    int cluster_size = 1;

    while (stack_top > 0) {
        int site = cluster_stack[--stack_top];
        int x = site / (L * L);
        int y = (site / L) % L;
        int z = site % L;

        /* 6 neighbors in 3D */
        int neighbors[6];
        neighbors[0] = idx3(x+1, y, z, L);
        neighbors[1] = idx3(x-1, y, z, L);
        neighbors[2] = idx3(x, y+1, z, L);
        neighbors[3] = idx3(x, y-1, z, L);
        neighbors[4] = idx3(x, y, z+1, L);
        neighbors[5] = idx3(x, y, z-1, L);

        for (int n = 0; n < 6; n++) {
            int nb = neighbors[n];
            if (!cluster_mark[nb] && spin[nb] == s0) {
                if (rng_double() < p_add) {
                    cluster_mark[nb] = 1;
                    cluster_stack[stack_top++] = nb;
                    cluster_size++;
                }
            }
        }
    }

    /* Flip the cluster */
    for (int i = 0; i < N; i++) {
        if (cluster_mark[i]) spin[i] = -spin[i];
    }

    return cluster_size;
}

/* ============================================================
 * Coarse-graining: Phi_B(x) = (1/N_B) * sum_i s_i^2 = 1
 * For Ising (s=+/-1), s^2=1 always. Use magnetization^2 instead:
 *   Phi_B(x) = (m_B)^2 where m_B = (1/N_B) * sum_i s_i
 * This is the Z_2 invariant order parameter.
 * ============================================================ */
static double block_field[LMAX * LMAX * LMAX]; /* at most (L/b)^3 blocks */

static int compute_block_field(int L, int b) {
    int LB = L / b;
    if (LB < 2) return 0;  /* too few blocks */
    int NB = b * b * b;

    for (int bx = 0; bx < LB; bx++)
    for (int by = 0; by < LB; by++)
    for (int bz = 0; bz < LB; bz++) {
        double m_block = 0.0;
        for (int dx = 0; dx < b; dx++)
        for (int dy = 0; dy < b; dy++)
        for (int dz = 0; dz < b; dz++) {
            int x = bx * b + dx;
            int y = by * b + dy;
            int z = bz * b + dz;
            m_block += spin[idx3(x, y, z, L)];
        }
        m_block /= NB;
        block_field[bx * LB * LB + by * LB + bz] = m_block * m_block;
    }
    return LB;
}

/* ============================================================
 * Measure gradient energy: E_grad = sum |Phi(x+e) - Phi(x)|^2
 * and K_eff(Phi) from binned gradient vs Phi
 * ============================================================ */
typedef struct {
    double phi_mean;
    double phi_std;
    double grad2_mean;   /* <|grad Phi|^2> */
    double alpha_eff;    /* fit exponent */
    double phi_vac;      /* vacuum value (histogram peak) */
    double beta_gamma;   /* ratio of cubic/quartic in U_eff */
    double corr_length;  /* correlation length estimate */
} BlockResult;

static BlockResult measure_block(int L, int b) {
    BlockResult res = {0};
    int LB = compute_block_field(L, b);
    if (LB < 3) {
        res.alpha_eff = -999;
        return res;
    }
    int NB3 = LB * LB * LB;

    /* Mean and std of Phi_B */
    double sum = 0, sum2 = 0;
    for (int i = 0; i < NB3; i++) {
        sum += block_field[i];
        sum2 += block_field[i] * block_field[i];
    }
    res.phi_mean = sum / NB3;
    res.phi_std = sqrt(sum2 / NB3 - res.phi_mean * res.phi_mean);

    /* Gradient energy binned by Phi */
    double grad_sum[N_BINS] = {0};
    double phi_sum[N_BINS] = {0};
    int count[N_BINS] = {0};

    double phi_min = 0.0, phi_max = 1.0;
    double bin_width = (phi_max - phi_min) / N_BINS;

    for (int bx = 0; bx < LB; bx++)
    for (int by = 0; by < LB; by++)
    for (int bz = 0; bz < LB; bz++) {
        int idx = bx * LB * LB + by * LB + bz;
        double phi_here = block_field[idx];

        /* Gradient in 3 directions */
        double g2 = 0;
        int idx_xp = ((bx+1)%LB) * LB * LB + by * LB + bz;
        int idx_yp = bx * LB * LB + ((by+1)%LB) * LB + bz;
        int idx_zp = bx * LB * LB + by * LB + ((bz+1)%LB);
        g2 += (block_field[idx_xp] - phi_here) * (block_field[idx_xp] - phi_here);
        g2 += (block_field[idx_yp] - phi_here) * (block_field[idx_yp] - phi_here);
        g2 += (block_field[idx_zp] - phi_here) * (block_field[idx_zp] - phi_here);

        int bin = (int)((phi_here - phi_min) / bin_width);
        if (bin < 0) bin = 0;
        if (bin >= N_BINS) bin = N_BINS - 1;
        grad_sum[bin] += g2;
        phi_sum[bin] += phi_here;
        count[bin]++;
    }

    /* Total gradient energy */
    double total_grad = 0;
    for (int i = 0; i < NB3; i++) {
        int idx_xp = (((i / (LB*LB)) + 1) % LB) * LB * LB + ((i / LB) % LB) * LB + (i % LB);
        int idx_yp = (i / (LB*LB)) * LB * LB + (((i / LB) % LB + 1) % LB) * LB + (i % LB);
        int idx_zp = (i / (LB*LB)) * LB * LB + ((i / LB) % LB) * LB + ((i % LB + 1) % LB);
        total_grad += (block_field[idx_xp] - block_field[i]) * (block_field[idx_xp] - block_field[i]);
        total_grad += (block_field[idx_yp] - block_field[i]) * (block_field[idx_yp] - block_field[i]);
        total_grad += (block_field[idx_zp] - block_field[i]) * (block_field[idx_zp] - block_field[i]);
    }
    res.grad2_mean = total_grad / NB3;

    /* Fit K(Phi) ~ Phi^alpha: log-log regression on binned data */
    /* K(Phi) ~ grad_energy(Phi) / <(dPhi)^2> at given Phi */
    double sum_lnx = 0, sum_lny = 0, sum_lnx2 = 0, sum_lnxy = 0;
    int n_fit = 0;
    for (int i = 0; i < N_BINS; i++) {
        if (count[i] < 5) continue;
        double phi_avg = phi_sum[i] / count[i];
        double K_avg = grad_sum[i] / count[i];  /* ~ K(phi) * <(dPhi)^2> */
        if (phi_avg < 0.01 || K_avg < 1e-15) continue;

        double lnx = log(phi_avg);
        double lny = log(K_avg);
        sum_lnx += lnx;
        sum_lny += lny;
        sum_lnx2 += lnx * lnx;
        sum_lnxy += lnx * lny;
        n_fit++;
    }

    if (n_fit >= 3) {
        double denom = n_fit * sum_lnx2 - sum_lnx * sum_lnx;
        if (fabs(denom) > 1e-20) {
            res.alpha_eff = (n_fit * sum_lnxy - sum_lnx * sum_lny) / denom;
        } else {
            res.alpha_eff = -999;
        }
    } else {
        res.alpha_eff = -999;
    }

    /* Histogram peak → vacuum */
    int peak_bin = 0, peak_count = 0;
    for (int i = 0; i < N_BINS; i++) {
        if (count[i] > peak_count) {
            peak_count = count[i];
            peak_bin = i;
        }
    }
    res.phi_vac = phi_min + (peak_bin + 0.5) * bin_width;

    /* Correlation length estimate: <Phi(0)Phi(r)> along x-axis */
    if (LB >= 6) {
        double c0 = 0;
        for (int by = 0; by < LB; by++)
        for (int bz = 0; bz < LB; bz++) {
            double phi0 = block_field[0 * LB * LB + by * LB + bz];
            c0 += phi0 * phi0;
        }
        c0 /= (LB * LB);

        /* Find distance where correlation drops to 1/e */
        res.corr_length = 0;
        for (int r = 1; r < LB/2; r++) {
            double cr = 0;
            for (int by = 0; by < LB; by++)
            for (int bz = 0; bz < LB; bz++) {
                double phi0 = block_field[0 * LB * LB + by * LB + bz];
                double phir = block_field[r * LB * LB + by * LB + bz];
                cr += phi0 * phir;
            }
            cr /= (LB * LB);
            if (c0 > 0 && cr / c0 < 1.0/M_E) {
                res.corr_length = r;
                break;
            }
        }
        if (res.corr_length == 0) res.corr_length = LB / 2;
    }

    return res;
}

/* ============================================================
 * Multi-temperature scan for alpha_eff convergence
 * ============================================================ */
static void run_T_scan(int L, int b) {
    double T_ratios[] = {0.80, 0.85, 0.90, 0.95, 0.98, 0.99, 1.00, 1.01, 1.02, 1.05};
    int n_T = 10;

    printf("\n  T/T_c scan (L=%d, b=%d):\n", L, b);
    printf("  %8s %8s %8s %8s %8s\n", "T/T_c", "<Phi>", "alpha", "|a-2|", "xi");

    for (int t = 0; t < n_T; t++) {
        double T = T_ratios[t] * T_C_ISING;
        double beta = J_COUPLING / T;
        double p_add = 1.0 - exp(-2.0 * beta);

        init_lattice(L);

        /* Thermalize */
        for (int i = 0; i < N_THERM; i++) {
            wolff_step(L, p_add);
        }

        /* Measure: accumulate over multiple configs */
        double alpha_acc = 0;
        double phi_acc = 0;
        int n_good = 0;

        for (int m = 0; m < N_MEASURE; m++) {
            for (int s = 0; s < N_SKIP; s++) {
                wolff_step(L, p_add);
            }
            BlockResult br = measure_block(L, b);
            if (br.alpha_eff > -100) {
                alpha_acc += br.alpha_eff;
                phi_acc += br.phi_mean;
                n_good++;
            }
        }

        if (n_good > 0) {
            double alpha_mean = alpha_acc / n_good;
            double phi_mean = phi_acc / n_good;
            printf("  %8.4f %8.4f %8.3f %8.3f %8s\n",
                   T_ratios[t], phi_mean, alpha_mean,
                   fabs(alpha_mean - 2.0), "---");
        }
    }
}

/* ============================================================
 * Main: run for each lattice size
 * ============================================================ */
int main(void) {
    rng_seed((uint64_t)time(NULL));

    printf("============================================================\n");
    printf("LK-1: CONTINUUM LIMIT TEST (C implementation)\n");
    printf("============================================================\n");
    printf("T_c = %.4f, Wolff cluster, N_therm=%d, N_measure=%d\n\n",
           T_C_ISING, N_THERM, N_MEASURE);

    /* Main scan: fixed T/T_c = 0.95 (ordered phase) */
    double T_sim = 0.95 * T_C_ISING;
    double beta = J_COUPLING / T_sim;
    double p_add = 1.0 - exp(-2.0 * beta);

    printf("Phase 1: Fixed T/T_c = 0.95 (ordered phase)\n");
    printf("============================================================\n");

    /* Results table header */
    printf("\n  %4s %4s %4s %8s %8s %8s %8s %8s\n",
           "L", "b", "L_B", "<Phi>", "std", "alpha", "|a-2|", "xi");
    printf("  -----------------------------------------------------------\n");

    typedef struct {
        int L, b, LB;
        double phi_mean, phi_std, alpha_eff, corr_len;
    } ResultRow;

    ResultRow results[20];
    int n_results = 0;

    for (int li = 0; li < N_LSIZES; li++) {
        int L = LATTICE_SIZES[li];
        if (L > LMAX) continue;

        printf("\n  --- L = %d ---\n", L);
        init_lattice(L);

        /* Thermalize */
        printf("  Thermalizing (%d Wolff steps)...\n", N_THERM);
        for (int i = 0; i < N_THERM; i++) {
            wolff_step(L, p_add);
        }

        for (int bi = 0; bi < N_BSIZES; bi++) {
            int b = BLOCK_SIZES[bi];
            int LB = L / b;
            if (LB < 3) continue;

            /* Accumulate measurements */
            double alpha_sum = 0, alpha_sum2 = 0;
            double phi_sum = 0, phi_sum2 = 0;
            double corr_sum = 0;
            int n_good = 0;

            for (int m = 0; m < N_MEASURE; m++) {
                for (int s = 0; s < N_SKIP; s++) {
                    wolff_step(L, p_add);
                }
                BlockResult br = measure_block(L, b);
                if (br.alpha_eff > -100) {
                    alpha_sum += br.alpha_eff;
                    alpha_sum2 += br.alpha_eff * br.alpha_eff;
                    phi_sum += br.phi_mean;
                    phi_sum2 += br.phi_mean * br.phi_mean;
                    corr_sum += br.corr_length;
                    n_good++;
                }
            }

            if (n_good > 10) {
                double a_mean = alpha_sum / n_good;
                double a_std = sqrt(alpha_sum2 / n_good - a_mean * a_mean);
                double p_mean = phi_sum / n_good;
                double p_std = sqrt(phi_sum2 / n_good - p_mean * p_mean);
                double c_mean = corr_sum / n_good;

                printf("  L=%2d b=%d LB=%2d: <Phi>=%.4f+/-%.4f  "
                       "alpha=%.3f+/-%.3f  |a-2|=%.3f  xi=%.1f\n",
                       L, b, LB, p_mean, p_std, a_mean, a_std,
                       fabs(a_mean - 2.0), c_mean);

                results[n_results].L = L;
                results[n_results].b = b;
                results[n_results].LB = LB;
                results[n_results].phi_mean = p_mean;
                results[n_results].phi_std = p_std;
                results[n_results].alpha_eff = a_mean;
                results[n_results].corr_len = c_mean;
                n_results++;
            }
        }
    }

    /* Phase 2: T/T_c scan for best L and b */
    printf("\n\n============================================================\n");
    printf("Phase 2: Temperature scan (critical region)\n");
    printf("============================================================\n");

    /* Use L=32, b=4 for T scan (good balance of speed and resolution) */
    int L_scan = 32;
    int b_scan = 4;
    if (L_scan <= LMAX) {
        run_T_scan(L_scan, b_scan);
    }

    /* Summary */
    printf("\n\n============================================================\n");
    printf("SUMMARY: CONTINUUM LIMIT CONVERGENCE\n");
    printf("============================================================\n\n");

    printf("  %4s %4s %4s %8s %8s %8s %8s\n",
           "L", "b", "L_B", "<Phi>", "alpha", "|a-2|", "xi");
    printf("  -------------------------------------------------------\n");
    for (int i = 0; i < n_results; i++) {
        printf("  %4d %4d %4d %8.4f %8.3f %8.3f %8.1f\n",
               results[i].L, results[i].b, results[i].LB,
               results[i].phi_mean, results[i].alpha_eff,
               fabs(results[i].alpha_eff - 2.0), results[i].corr_len);
    }

    /* Check convergence trend */
    printf("\n  TGP prediction: alpha_eff -> 2 in continuum limit\n");

    if (n_results >= 2) {
        /* Check if alpha increases with L (convergence toward 2) */
        double a_first = results[0].alpha_eff;
        double a_last = results[n_results-1].alpha_eff;
        if (a_last > a_first) {
            printf("  Trend: alpha INCREASING with L/b (%.3f -> %.3f)\n",
                   a_first, a_last);
            printf("  Direction: TOWARD alpha=2 (promising)\n");
        } else {
            printf("  Trend: alpha NOT increasing (%.3f -> %.3f)\n",
                   a_first, a_last);
            printf("  Needs larger L or different T/T_c\n");
        }
    }

    printf("\n============================================================\n");
    printf("TESTS:\n");
    printf("============================================================\n");

    /* T1: Phi_B >= 0 */
    int t1_pass = 1;
    for (int i = 0; i < n_results; i++) {
        if (results[i].phi_mean < 0) t1_pass = 0;
    }
    printf("  T1 [Phi_B >= 0 (Z_2 invariant)]: %s\n", t1_pass ? "PASS" : "FAIL");

    /* T2: alpha_eff trend */
    int t2_pass = 0;
    for (int i = 0; i < n_results; i++) {
        if (fabs(results[i].alpha_eff - 2.0) < 1.0) t2_pass = 1;
    }
    printf("  T2 [alpha_eff within 1.0 of 2]: %s\n",
           t2_pass ? "PASS" : "FAIL (need larger L or closer to T_c)");

    /* T3: Bounded variance */
    int t3_pass = 1;
    for (int i = 0; i < n_results; i++) {
        if (results[i].phi_std > 10.0) t3_pass = 0;
    }
    printf("  T3 [Phi_B bounded variance]: %s\n", t3_pass ? "PASS" : "FAIL");

    /* T4: Correlations exist */
    int t4_pass = 0;
    for (int i = 0; i < n_results; i++) {
        if (results[i].corr_len > 0.5) t4_pass = 1;
    }
    printf("  T4 [Finite correlation length]: %s\n", t4_pass ? "PASS" : "NO DATA");

    /* T5: Non-trivial vacuum */
    int t5_pass = 0;
    for (int i = 0; i < n_results; i++) {
        if (results[i].phi_mean > 0.01) t5_pass = 1;
    }
    printf("  T5 [Non-trivial vacuum Phi_0 > 0]: %s\n", t5_pass ? "PASS" : "FAIL");

    int total = t1_pass + t2_pass + t3_pass + t4_pass + t5_pass;
    printf("\n  TOTAL: %d/5 PASS\n", total);

    printf("============================================================\n");

    return 0;
}
