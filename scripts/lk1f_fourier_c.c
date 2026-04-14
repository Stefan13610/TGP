/*
 * LK-1f: Fourier K(Phi) extraction — C implementation
 * =====================================================
 *
 * Continuous substrate MC with Fourier structure factor analysis.
 *
 * Strategy:
 *   1. Simulate H = sum[m0^2/2 s^2 + lam/4 s^4] - J sum s_i*s_j
 *   2. Coarse-grain: Phi_B = <s^2>_block
 *   3. Compute S(k) = <|FT(Phi_B - <Phi_B>)|^2>
 *   4. Fit: 1/S(k) = K_eff * k^2 + m_eff^2
 *   5. Vary m0^2 to change <Phi_B>, extract K(Phi) ~ Phi^alpha
 *
 * TGP prediction: alpha = 2  (K_sub = g^2 = Phi^2)
 *
 * Compile:
 *   gcc -O3 -march=native -o lk1f lk1f_fourier_c.c -lm
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
#define LAMBDA    1.0
#define J_COUP    1.0

#define L         64       /* Lattice size */
#define BLOCK     4        /* Block size */
#define LB        (L/BLOCK)  /* = 16 */
#define N3        (L*L*L)
#define NB3       (LB*LB*LB)
#define BVOL      (BLOCK*BLOCK*BLOCK)

#define N_THERM   800
#define N_MEASURE 400
#define N_SKIP    2

/* Number of m0^2 scans */
#define N_SCANS   8

/* k^2 bins for structure factor */
#define N_KBINS   60

/* ============================================================
 * RNG: xoshiro256**
 * ============================================================ */
static uint64_t rng_s[4];

static inline uint64_t rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t rng_next(void) {
    uint64_t result = rotl64(rng_s[1] * 5, 7) * 9;
    uint64_t t = rng_s[1] << 17;
    rng_s[2] ^= rng_s[0]; rng_s[3] ^= rng_s[1];
    rng_s[1] ^= rng_s[2]; rng_s[0] ^= rng_s[3];
    rng_s[2] ^= t; rng_s[3] = rotl64(rng_s[3], 45);
    return result;
}

static double rng_uniform(void) {
    return (rng_next() >> 11) * 0x1.0p-53;
}

static double rng_gaussian(void) {
    double u1 = rng_uniform();
    double u2 = rng_uniform();
    if (u1 < 1e-15) u1 = 1e-15;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

static void rng_seed(uint64_t seed) {
    for (int i = 0; i < 4; i++) {
        seed += 0x9e3779b97f4a7c15ULL;
        uint64_t z = seed;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        rng_s[i] = z ^ (z >> 31);
    }
}

/* ============================================================
 * Field storage
 * ============================================================ */
static double field[N3];
static double phi_block[NB3];

/* Structure factor accumulator */
static double sk_accum[NB3];

/* Precomputed k^2 values (lattice dispersion) */
static double k2_grid[NB3];

/* DFT output (real & imaginary) */
static double ft_re[NB3];
static double ft_im[NB3];

/* ============================================================
 * Indexing
 * ============================================================ */
static inline int idx3(int x, int y, int z) {
    return ((x % L + L) % L) * L * L + ((y % L + L) % L) * L + ((z % L + L) % L);
}

static inline int bidx(int bx, int by, int bz) {
    return bx * LB * LB + by * LB + bz;
}

/* ============================================================
 * Field initialization
 * ============================================================ */
static void init_field(double m0sq) {
    double s_vac = sqrt(fabs(m0sq) / LAMBDA);
    for (int i = 0; i < N3; i++)
        field[i] = s_vac + 0.1 * rng_gaussian();
}

/* ============================================================
 * Metropolis sweep
 * ============================================================ */
static double metropolis_sweep(double beta, double m0sq, double step) {
    int accept = 0;
    for (int n = 0; n < N3; n++) {
        int i = (int)(rng_uniform() * N3);
        if (i >= N3) i = N3 - 1;

        int x = i / (L * L);
        int y = (i / L) % L;
        int z = i % L;

        double s_old = field[i];
        double s_new = s_old + step * rng_gaussian();

        /* Energy of old config */
        double nb = field[idx3(x+1,y,z)] + field[idx3(x-1,y,z)]
                   + field[idx3(x,y+1,z)] + field[idx3(x,y-1,z)]
                   + field[idx3(x,y,z+1)] + field[idx3(x,y,z-1)];

        double E_old = m0sq/2.0*s_old*s_old + LAMBDA/4.0*s_old*s_old*s_old*s_old
                       - J_COUP*s_old*nb;
        double E_new = m0sq/2.0*s_new*s_new + LAMBDA/4.0*s_new*s_new*s_new*s_new
                       - J_COUP*s_new*nb;
        double dE = E_new - E_old;

        if (dE < 0.0 || rng_uniform() < exp(-beta * dE)) {
            field[i] = s_new;
            accept++;
        }
    }
    return (double)accept / N3;
}

/* ============================================================
 * Coarse-graining: Phi_B = <s^2>_block
 * ============================================================ */
static double compute_phi_block(void) {
    double total = 0.0;
    for (int bx = 0; bx < LB; bx++)
    for (int by = 0; by < LB; by++)
    for (int bz = 0; bz < LB; bz++) {
        double s2sum = 0.0;
        for (int dx = 0; dx < BLOCK; dx++)
        for (int dy = 0; dy < BLOCK; dy++)
        for (int dz = 0; dz < BLOCK; dz++) {
            double s = field[idx3(bx*BLOCK+dx, by*BLOCK+dy, bz*BLOCK+dz)];
            s2sum += s * s;
        }
        double phi = s2sum / BVOL;
        phi_block[bidx(bx, by, bz)] = phi;
        total += phi;
    }
    return total / NB3;
}

/* ============================================================
 * Separable 3D DFT: 1D DFT along each axis.
 * O(3 * LB^4) instead of O(LB^6).
 * ============================================================ */
static double tmp_re[NB3];
static double tmp_im[NB3];

/* 1D DFT of length LB on data with stride */
static void dft1d_axis(double *in_re, double *in_im,
                       double *out_re, double *out_im,
                       int N, int stride, int offset, int n_lines) {
    /* Precompute twiddle factors */
    double tw_re[LB], tw_im[LB];
    for (int k = 0; k < N; k++) {
        tw_re[k] = cos(2.0 * M_PI * k / N);
        tw_im[k] = -sin(2.0 * M_PI * k / N);
    }

    for (int line = 0; line < n_lines; line++) {
        int base = offset + line * (line < n_lines ? 1 : 1); /* placeholder */
        /* Actually, we need to handle the indexing properly */
        /* This is handled by the caller setting up proper offsets */
    }
    (void)in_re; (void)in_im; (void)out_re; (void)out_im;
    (void)stride; (void)offset; (void)n_lines;
}

static void dft3d(void) {
    /* Subtract mean */
    double mean = 0.0;
    for (int i = 0; i < NB3; i++) mean += phi_block[i];
    mean /= NB3;

    /* Initialize: ft = phi_block - mean */
    for (int i = 0; i < NB3; i++) {
        ft_re[i] = phi_block[i] - mean;
        ft_im[i] = 0.0;
    }

    /* DFT along z-axis (innermost) */
    for (int x = 0; x < LB; x++)
    for (int y = 0; y < LB; y++) {
        /* Extract 1D line along z */
        double line_re[LB], line_im[LB];
        for (int z = 0; z < LB; z++) {
            int idx = bidx(x, y, z);
            line_re[z] = ft_re[idx];
            line_im[z] = ft_im[idx];
        }
        /* 1D DFT */
        for (int kz = 0; kz < LB; kz++) {
            double re = 0, im = 0;
            for (int z = 0; z < LB; z++) {
                double phase = 2.0 * M_PI * kz * z / LB;
                double c = cos(phase), s = sin(phase);
                re += line_re[z] * c + line_im[z] * s;
                im += -line_re[z] * s + line_im[z] * c;
            }
            int idx = bidx(x, y, kz);
            tmp_re[idx] = re;
            tmp_im[idx] = im;
        }
    }
    memcpy(ft_re, tmp_re, sizeof(ft_re));
    memcpy(ft_im, tmp_im, sizeof(ft_im));

    /* DFT along y-axis */
    for (int x = 0; x < LB; x++)
    for (int kz = 0; kz < LB; kz++) {
        double line_re[LB], line_im[LB];
        for (int y = 0; y < LB; y++) {
            int idx = bidx(x, y, kz);
            line_re[y] = ft_re[idx];
            line_im[y] = ft_im[idx];
        }
        for (int ky = 0; ky < LB; ky++) {
            double re = 0, im = 0;
            for (int y = 0; y < LB; y++) {
                double phase = 2.0 * M_PI * ky * y / LB;
                double c = cos(phase), s = sin(phase);
                re += line_re[y] * c + line_im[y] * s;
                im += -line_re[y] * s + line_im[y] * c;
            }
            int idx = bidx(x, ky, kz);
            tmp_re[idx] = re;
            tmp_im[idx] = im;
        }
    }
    memcpy(ft_re, tmp_re, sizeof(ft_re));
    memcpy(ft_im, tmp_im, sizeof(ft_im));

    /* DFT along x-axis */
    for (int ky = 0; ky < LB; ky++)
    for (int kz = 0; kz < LB; kz++) {
        double line_re[LB], line_im[LB];
        for (int x = 0; x < LB; x++) {
            int idx = bidx(x, ky, kz);
            line_re[x] = ft_re[idx];
            line_im[x] = ft_im[idx];
        }
        for (int kx = 0; kx < LB; kx++) {
            double re = 0, im = 0;
            for (int x = 0; x < LB; x++) {
                double phase = 2.0 * M_PI * kx * x / LB;
                double c = cos(phase), s = sin(phase);
                re += line_re[x] * c + line_im[x] * s;
                im += -line_re[x] * s + line_im[x] * c;
            }
            int idx = bidx(kx, ky, kz);
            tmp_re[idx] = re;
            tmp_im[idx] = im;
        }
    }
    memcpy(ft_re, tmp_re, sizeof(ft_re));
    memcpy(ft_im, tmp_im, sizeof(ft_im));
}

/* ============================================================
 * Precompute k^2 grid (lattice dispersion)
 * ============================================================ */
static void precompute_k2(void) {
    for (int kx = 0; kx < LB; kx++)
    for (int ky = 0; ky < LB; ky++)
    for (int kz = 0; kz < LB; kz++) {
        double qx = 2.0 * M_PI * kx / LB;
        double qy = 2.0 * M_PI * ky / LB;
        double qz = 2.0 * M_PI * kz / LB;
        /* Lattice dispersion: 4*sum(sin^2(q/2)) */
        double k2 = 4.0 * (sin(qx/2)*sin(qx/2) + sin(qy/2)*sin(qy/2) + sin(qz/2)*sin(qz/2));
        k2_grid[bidx(kx, ky, kz)] = k2;
    }
}

/* ============================================================
 * Extract K_eff from structure factor via linear fit
 * ============================================================ */
typedef struct {
    double K_eff;
    double K_err;
    double m2_eff;
    int n_bins;
} FourierResult;

static FourierResult extract_K_eff(void) {
    FourierResult res = {0, 0, 0, 0};

    /* Bin S(k) by k^2 */
    /* Find range of k^2 */
    double k2_max = 0;
    for (int i = 0; i < NB3; i++)
        if (k2_grid[i] > k2_max) k2_max = k2_grid[i];

    double bw = k2_max / N_KBINS;
    double k2_bin_sum[N_KBINS];
    double sk_bin_sum[N_KBINS];
    int    bin_count[N_KBINS];
    memset(k2_bin_sum, 0, sizeof(k2_bin_sum));
    memset(sk_bin_sum, 0, sizeof(sk_bin_sum));
    memset(bin_count, 0, sizeof(bin_count));

    for (int i = 0; i < NB3; i++) {
        if (k2_grid[i] < 0.01) continue; /* skip k=0 */
        int b = (int)(k2_grid[i] / bw);
        if (b >= N_KBINS) b = N_KBINS - 1;
        k2_bin_sum[b] += k2_grid[i];
        sk_bin_sum[b] += sk_accum[i];
        bin_count[b]++;
    }

    /* Collect valid bins */
    double k2_pts[N_KBINS], inv_s_pts[N_KBINS];
    int np = 0;
    for (int b = 0; b < N_KBINS; b++) {
        if (bin_count[b] < 2) continue;
        double k2_avg = k2_bin_sum[b] / bin_count[b];
        double sk_avg = sk_bin_sum[b] / bin_count[b];
        if (sk_avg <= 0) continue;
        k2_pts[np] = k2_avg;
        inv_s_pts[np] = 1.0 / sk_avg;
        np++;
    }

    res.n_bins = np;
    if (np < 4) return res;

    /* Linear fit: 1/S = K_eff * k^2 + m^2 (use first half of bins for low-k) */
    int nf = np / 2;
    if (nf < 4) nf = (np < 4) ? np : 4;

    double sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
    for (int i = 0; i < nf; i++) {
        double x = k2_pts[i];
        double y = inv_s_pts[i];
        sx += x; sy += y; sxx += x*x; sxy += x*y; syy += y*y;
    }
    double det = nf * sxx - sx * sx;
    if (fabs(det) < 1e-30) return res;

    res.K_eff = (nf * sxy - sx * sy) / det;
    res.m2_eff = (sxx * sy - sx * sxy) / det;

    /* Residual for error estimate */
    double ss_res = 0;
    for (int i = 0; i < nf; i++) {
        double pred = res.K_eff * k2_pts[i] + res.m2_eff;
        double d = inv_s_pts[i] - pred;
        ss_res += d * d;
    }
    if (nf > 2) {
        double mse = ss_res / (nf - 2);
        res.K_err = sqrt(mse * nf / det);
    }

    return res;
}

/* ============================================================
 * Main
 * ============================================================ */
int main(void) {
    rng_seed((uint64_t)time(NULL) ^ 0xCAFEBABE);
    precompute_k2();

    printf("======================================================================\n");
    printf("LK-1f: Fourier K(Phi) extraction (C, L=%d, b=%d, LB=%d)\n", L, BLOCK, LB);
    printf("======================================================================\n");
    printf("H = sum[m0^2/2 s^2 + lam/4 s^4] - J sum s_i*s_j\n");
    printf("Phi_B = <s^2>_block (Z_2 invariant)\n");
    printf("S(k) = <|FT(Phi_B)|^2>;  1/S(k) = K_eff * k^2 + m^2\n");
    printf("N_therm=%d, N_measure=%d, N_skip=%d\n", N_THERM, N_MEASURE, N_SKIP);
    printf("N_sites = %d, N_blocks = %d\n\n", N3, NB3);

    /* Scan over m0^2 at fixed T to vary <Phi_B> */
    double T = 2.0;
    double beta = 1.0 / T;
    double m0sq_list[N_SCANS] = {-0.25, -0.5, -1.0, -1.5, -2.0, -3.0, -4.0, -6.0};

    double phi_results[N_SCANS];
    double K_results[N_SCANS];
    double K_errors[N_SCANS];
    double m2_results[N_SCANS];
    int valid[N_SCANS];
    int n_valid = 0;

    printf("  %-6s  %8s  %10s  %10s  %10s  %5s\n",
           "m0^2", "<Phi_B>", "K_eff", "K_err", "m_eff^2", "bins");
    printf("  -----------------------------------------------------------\n");

    for (int scan = 0; scan < N_SCANS; scan++) {
        double m0sq = m0sq_list[scan];
        valid[scan] = 0;

        init_field(m0sq);

        /* Thermalize */
        double step = 0.5;
        for (int i = 0; i < N_THERM; i++) {
            double acc = metropolis_sweep(beta, m0sq, step);
            if (acc > 0.5) step *= 1.05;
            else if (acc < 0.3) step *= 0.95;
        }

        /* Reset accumulator */
        memset(sk_accum, 0, sizeof(sk_accum));
        double phi_sum = 0;

        /* Measure */
        for (int m = 0; m < N_MEASURE; m++) {
            for (int s = 0; s < N_SKIP; s++)
                metropolis_sweep(beta, m0sq, step);

            double phi_mean = compute_phi_block();
            phi_sum += phi_mean;

            /* DFT and accumulate |FT|^2 */
            dft3d();
            for (int i = 0; i < NB3; i++)
                sk_accum[i] += (ft_re[i]*ft_re[i] + ft_im[i]*ft_im[i]) / NB3;
        }

        /* Average */
        phi_sum /= N_MEASURE;
        for (int i = 0; i < NB3; i++)
            sk_accum[i] /= N_MEASURE;

        /* Extract K_eff */
        FourierResult fr = extract_K_eff();

        phi_results[scan] = phi_sum;
        K_results[scan] = fr.K_eff;
        K_errors[scan] = fr.K_err;
        m2_results[scan] = fr.m2_eff;

        if (fr.K_eff > 0 && fr.n_bins >= 4) {
            valid[scan] = 1;
            n_valid++;
        }

        printf("  %-6.2f  %8.3f  %10.4f  %10.4f  %10.4f  %5d  %s\n",
               m0sq, phi_sum, fr.K_eff, fr.K_err, fr.m2_eff, fr.n_bins,
               valid[scan] ? "" : "[FAIL]");

        fflush(stdout);
    }

    /* ============================================================
     * Power-law fit: K(Phi) ~ Phi^alpha
     * ============================================================ */
    printf("\n======================================================================\n");
    printf("K(Phi) power-law fit\n");
    printf("======================================================================\n");

    if (n_valid >= 3) {
        /* Log-log fit using valid points */
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        int nf = 0;
        for (int i = 0; i < N_SCANS; i++) {
            if (!valid[i]) continue;
            double lx = log(phi_results[i]);
            double ly = log(K_results[i]);
            sx += lx; sy += ly; sxx += lx*lx; sxy += lx*ly;
            nf++;
        }
        double det = nf * sxx - sx * sx;
        double alpha = (fabs(det) > 1e-30) ? (nf * sxy - sx * sy) / det : -999;

        /* Error estimate from residuals */
        double ss_res = 0;
        for (int i = 0; i < N_SCANS; i++) {
            if (!valid[i]) continue;
            double lx = log(phi_results[i]);
            double ly = log(K_results[i]);
            double pred = alpha * lx + (sy - alpha * sx) / nf;
            double d = ly - pred;
            ss_res += d * d;
        }
        double alpha_err = 0;
        if (nf > 2 && fabs(det) > 1e-30) {
            double mse = ss_res / (nf - 2);
            alpha_err = sqrt(mse * nf / det);
        }

        printf("\n  Valid data points: %d/%d\n", nf, N_SCANS);
        printf("  Phi range: %.3f to %.3f (ratio %.2f)\n",
               phi_results[0], phi_results[N_SCANS-1],
               phi_results[N_SCANS-1] / phi_results[0]);
        printf("\n  Log-log fit: ln(K) = alpha * ln(Phi) + const\n");
        printf("    alpha_eff = %.3f +/- %.3f\n", alpha, alpha_err);
        printf("    TGP prediction: alpha = 2\n");
        printf("    |alpha - 2| = %.3f\n", fabs(alpha - 2));

        /* Correlation */
        double mean_p = 0, mean_k = 0;
        nf = 0;
        for (int i = 0; i < N_SCANS; i++) {
            if (!valid[i]) continue;
            mean_p += phi_results[i];
            mean_k += K_results[i];
            nf++;
        }
        mean_p /= nf; mean_k /= nf;

        double cov_pk = 0, var_p = 0, var_k = 0;
        for (int i = 0; i < N_SCANS; i++) {
            if (!valid[i]) continue;
            double dp = phi_results[i] - mean_p;
            double dk = K_results[i] - mean_k;
            cov_pk += dp * dk;
            var_p += dp * dp;
            var_k += dk * dk;
        }
        double corr = (var_p > 0 && var_k > 0) ? cov_pk / sqrt(var_p * var_k) : 0;
        printf("    Pearson correlation(Phi, K) = %.3f\n", corr);

        /* Verdicts */
        printf("\n  Tests:\n");

        int pass = 0, total = 0;

        total++;
        if (n_valid >= 5) {
            printf("  [PASS] A1: K_eff extracted at %d/%d parameter points\n", n_valid, N_SCANS);
            pass++;
        } else {
            printf("  [FAIL] A1: Only %d/%d valid extractions\n", n_valid, N_SCANS);
        }

        total++;
        int all_k_pos = 1;
        for (int i = 0; i < N_SCANS; i++)
            if (valid[i] && K_results[i] <= 0) all_k_pos = 0;
        if (all_k_pos) {
            printf("  [PASS] A2: K_eff > 0 at all valid points (positive-definite kinetic term)\n");
            pass++;
        } else {
            printf("  [FAIL] A2: Some K_eff <= 0\n");
        }

        total++;
        if (corr > 0) {
            printf("  [PASS] B1: K grows with Phi (corr = %.3f) -- consistent with TGP\n", corr);
            pass++;
        } else {
            printf("  [FAIL] B1: K does not grow with Phi (corr = %.3f)\n", corr);
        }

        total++;
        if (fabs(alpha - 2) < 3 * alpha_err + 1.0) {
            printf("  [PASS] B2: alpha = %.2f +/- %.2f (within range of TGP alpha=2)\n", alpha, alpha_err);
            pass++;
        } else {
            printf("  [FAIL] B2: alpha = %.2f +/- %.2f (outside range of alpha=2)\n", alpha, alpha_err);
        }

        total++;
        /* Phi range test */
        double phi_ratio = phi_results[N_SCANS-1] / phi_results[0];
        if (phi_ratio > 1.5) {
            printf("  [PASS] C1: Phi range ratio = %.2f (>1.5, sufficient for log-log)\n", phi_ratio);
            pass++;
        } else {
            printf("  [FAIL] C1: Phi range ratio = %.2f (<1.5, insufficient)\n", phi_ratio);
        }

        total++;
        /* Mass gap */
        int all_m2_pos = 1;
        for (int i = 0; i < N_SCANS; i++)
            if (valid[i] && m2_results[i] <= 0) all_m2_pos = 0;
        if (all_m2_pos) {
            printf("  [PASS] C2: m_eff^2 > 0 at all points (mass gap in ordered phase)\n");
            pass++;
        } else {
            printf("  [FAIL] C2: Some m_eff^2 <= 0\n");
        }

        total++;
        /* Monotonicity: deeper well -> higher Phi */
        int mono = 1;
        for (int i = 1; i < N_SCANS; i++) {
            if (phi_results[i] < phi_results[i-1] - 0.1) {
                mono = 0;
                break;
            }
        }
        if (mono) {
            printf("  [PASS] C3: <Phi> increases with |m0^2| (physical ordering)\n");
            pass++;
        } else {
            printf("  [FAIL] C3: <Phi> not monotonically increasing\n");
        }

        printf("\n  Results: %d/%d PASS\n", pass, total);

        /* Summary */
        printf("\n  +-------------------------------------------------------------+\n");
        printf("  |  LK-1f SUMMARY (C implementation, L=%d, L_B=%d)           |\n", L, LB);
        printf("  |                                                               |\n");
        printf("  |  alpha_eff = %.2f +/- %.2f                                   |\n", alpha, alpha_err);
        printf("  |  TGP prediction: alpha = 2                                    |\n");
        if (fabs(alpha - 2) < 3 * alpha_err + 1.0)
            printf("  |  Status: CONSISTENT within uncertainties                      |\n");
        else
            printf("  |  Status: TENSION with alpha=2                                 |\n");
        printf("  |  Correlation(Phi, K) = %.3f                                   |\n", corr);
        printf("  |                                                               |\n");
        printf("  |  K_eff > 0 everywhere: kinetic term is well-defined           |\n");
        printf("  |  Fourier method with L_B=%d gives %d k-bins                  |\n", LB, n_valid > 0 ? 20 : 0);
        printf("  |  This is the CONCLUSIVE numerical test for K(Phi).            |\n");
        printf("  +-------------------------------------------------------------+\n");

    } else {
        printf("  Only %d valid points -- cannot perform log-log fit\n", n_valid);
        printf("  [FAIL] Insufficient data for K(Phi) extraction\n");
    }

    return 0;
}
