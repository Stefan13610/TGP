/*
 * LK-1b: Continuous Substrate MC — TGP coarse-graining test
 * ===========================================================
 * Full TGP substrate with continuous field s_i in R:
 *
 *   H_Gamma = sum_i [m0^2/2 * s_i^2 + lambda/4 * s_i^4] - J * sum_<ij> s_i * s_j
 *
 * Block field (Z_2 invariant): Phi_B(x) = (1/N_B) * sum_{i in block} s_i^2
 *
 * This is the CORRECT substrate for TGP — unlike binary Ising,
 * the continuous field with lambda*s^4 self-interaction generates
 * a non-trivial kinetic kernel K(Phi) in the effective theory.
 *
 * TGP prediction: K(Phi) ~ Phi^alpha with alpha = 2 (canonical)
 *
 * Compile: gcc -O3 -march=native -o lk1_continuous lk1_continuous_substrate_mc.c -lm
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
/* Substrate Hamiltonian parameters */
#define M0_SQ    -1.0     /* m_0^2: NEGATIVE for broken Z_2 phase */
#define LAMBDA    1.0     /* lambda: self-interaction coupling */
#define J_COUP    1.0     /* J: nearest-neighbor coupling */

/* Simulation */
#define N_THERM   1000    /* Thermalization sweeps */
#define N_MEASURE 500     /* Measurement sweeps */
#define N_SKIP    1       /* Sweeps between measurements */

/* Lattice */
static const int LATTICE_SIZES[] = {32, 64};
#define N_LSIZES 2
static const int BLOCK_SIZES[] = {4, 8, 16};
#define N_BSIZES 3
#define LMAX 64

/* Temperature scan */
#define N_TEMPS 8

/* Histogram */
#define N_BINS 40

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

/* Box-Muller for Gaussian */
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
 * Lattice
 * ============================================================ */
static double field[LMAX * LMAX * LMAX];  /* continuous s_i */
static double block_phi[LMAX * LMAX * LMAX];  /* Phi_B */

static inline int idx3(int x, int y, int z, int L) {
    return ((x % L + L) % L) * L * L + ((y % L + L) % L) * L + ((z % L + L) % L);
}

/* Local energy at site i: m0^2/2 * s^2 + lambda/4 * s^4 - J*s*sum(neighbors) */
static double local_energy(int x, int y, int z, int L, double s) {
    double E = M0_SQ / 2.0 * s * s + LAMBDA / 4.0 * s * s * s * s;
    /* Neighbor sum */
    double nb_sum = field[idx3(x+1,y,z,L)] + field[idx3(x-1,y,z,L)]
                  + field[idx3(x,y+1,z,L)] + field[idx3(x,y-1,z,L)]
                  + field[idx3(x,y,z+1,L)] + field[idx3(x,y,z-1,L)];
    E -= J_COUP * s * nb_sum;
    return E;
}

/* Initialize: Gaussian around one of the minima */
static void init_field(int L) {
    /* V(s) = m0^2/2 s^2 + lambda/4 s^4
       With m0^2 < 0: V'(s) = m0^2 s + lambda s^3 = s(m0^2 + lambda s^2) = 0
       Minimum at s = +/- sqrt(|m0^2|/lambda)  (broken Z_2)
    */
    double s_vac = sqrt(fabs(M0_SQ) / LAMBDA);  /* vacuum value */
    int N = L * L * L;
    for (int i = 0; i < N; i++) {
        field[i] = s_vac + 0.1 * rng_gaussian();  /* start near +vacuum */
    }
}

/* Metropolis sweep with adaptive step size */
static double metropolis_sweep(int L, double beta, double step) {
    int N = L * L * L;
    int accept = 0;

    for (int n = 0; n < N; n++) {
        int i = (int)(rng_uniform() * N);
        if (i >= N) i = N - 1;

        int x = i / (L * L);
        int y = (i / L) % L;
        int z = i % L;

        double s_old = field[i];
        double s_new = s_old + step * rng_gaussian();

        double E_old = local_energy(x, y, z, L, s_old);
        double E_new = local_energy(x, y, z, L, s_new);
        double dE = E_new - E_old;

        if (dE < 0 || rng_uniform() < exp(-beta * dE)) {
            field[i] = s_new;
            accept++;
        }
    }
    return (double)accept / N;
}

/* Over-relaxation sweep (for continuous fields) */
static void overrelax_sweep(int L, double beta) {
    int N = L * L * L;
    for (int i = 0; i < N; i++) {
        int x = i / (L * L);
        int y = (i / L) % L;
        int z = i % L;

        double nb_sum = field[idx3(x+1,y,z,L)] + field[idx3(x-1,y,z,L)]
                      + field[idx3(x,y+1,z,L)] + field[idx3(x,y-1,z,L)]
                      + field[idx3(x,y,z+1,L)] + field[idx3(x,y,z-1,L)];

        /* For quadratic part: reflect around conditional mean */
        /* s_new = 2*<s>_cond - s_old, where <s>_cond = J*nb_sum / (m0^2 + lambda*s^2) */
        /* This is exact for Gaussian; approximate for s^4 but helps mixing */
        double s = field[i];
        double denom = M0_SQ + LAMBDA * s * s;
        if (fabs(denom) > 1e-10) {
            double s_mean = J_COUP * nb_sum / denom;
            field[i] = 2.0 * s_mean - s;
        }
    }
}

/* ============================================================
 * Compute block field Phi_B = <s^2>_block
 * ============================================================ */
static int compute_block_phi(int L, int b) {
    int LB = L / b;
    if (LB < 3) return 0;
    int NB = b * b * b;

    for (int bx = 0; bx < LB; bx++)
    for (int by = 0; by < LB; by++)
    for (int bz = 0; bz < LB; bz++) {
        double s2_sum = 0.0;
        for (int dx = 0; dx < b; dx++)
        for (int dy = 0; dy < b; dy++)
        for (int dz = 0; dz < b; dz++) {
            double s = field[idx3(bx*b+dx, by*b+dy, bz*b+dz, L)];
            s2_sum += s * s;
        }
        block_phi[bx * LB * LB + by * LB + bz] = s2_sum / NB;
    }
    return LB;
}

/* ============================================================
 * Measure alpha_eff from gradient energy
 * ============================================================ */
typedef struct {
    double phi_mean, phi_std;
    double alpha_eff;
    double xi_corr;
    int n_fit_pts;
} MeasResult;

static MeasResult measure_alpha(int L, int b) {
    MeasResult res = {0, 0, -999, 0, 0};
    int LB = compute_block_phi(L, b);
    if (LB < 3) return res;
    int NB3 = LB * LB * LB;

    /* Mean and std */
    double sum = 0, sum2 = 0;
    for (int i = 0; i < NB3; i++) {
        sum += block_phi[i];
        sum2 += block_phi[i] * block_phi[i];
    }
    res.phi_mean = sum / NB3;
    res.phi_std = sqrt(fabs(sum2 / NB3 - res.phi_mean * res.phi_mean));

    /* Bin gradient energy by Phi value */
    double grad_sum[N_BINS] = {0};
    double phi_bsum[N_BINS] = {0};
    int count[N_BINS] = {0};

    /* Auto-range from data */
    double phi_min = 1e30, phi_max = -1e30;
    for (int i = 0; i < NB3; i++) {
        if (block_phi[i] < phi_min) phi_min = block_phi[i];
        if (block_phi[i] > phi_max) phi_max = block_phi[i];
    }
    if (phi_max - phi_min < 1e-10) {
        phi_min -= 0.1;
        phi_max += 0.1;
    }
    double bw = (phi_max - phi_min) / N_BINS;

    for (int bx = 0; bx < LB; bx++)
    for (int by = 0; by < LB; by++)
    for (int bz = 0; bz < LB; bz++) {
        int idx = bx * LB * LB + by * LB + bz;
        double phi_here = block_phi[idx];

        /* Gradient in 3 directions */
        double g2 = 0;
        int inx = ((bx+1)%LB)*LB*LB + by*LB + bz;
        int iny = bx*LB*LB + ((by+1)%LB)*LB + bz;
        int inz = bx*LB*LB + by*LB + ((bz+1)%LB);
        g2 += (block_phi[inx] - phi_here) * (block_phi[inx] - phi_here);
        g2 += (block_phi[iny] - phi_here) * (block_phi[iny] - phi_here);
        g2 += (block_phi[inz] - phi_here) * (block_phi[inz] - phi_here);

        int bin = (int)((phi_here - phi_min) / bw);
        if (bin < 0) bin = 0;
        if (bin >= N_BINS) bin = N_BINS - 1;
        grad_sum[bin] += g2;
        phi_bsum[bin] += phi_here;
        count[bin]++;
    }

    /* Log-log fit: K(Phi) ~ Phi^alpha */
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int nf = 0;
    for (int i = 0; i < N_BINS; i++) {
        if (count[i] < 3) continue;
        double p = phi_bsum[i] / count[i];
        double k = grad_sum[i] / count[i];
        if (p < 1e-6 || k < 1e-20) continue;

        double lx = log(p), ly = log(k);
        sx += lx; sy += ly; sxx += lx*lx; sxy += lx*ly;
        nf++;
    }
    res.n_fit_pts = nf;

    if (nf >= 2) {
        double det = nf * sxx - sx * sx;
        if (fabs(det) > 1e-20)
            res.alpha_eff = (nf * sxy - sx * sy) / det;
    }

    /* Correlation length from <Phi(0)Phi(r)> */
    if (LB >= 6) {
        double c0 = 0;
        int cnt0 = 0;
        for (int by = 0; by < LB; by++)
        for (int bz = 0; bz < LB; bz++) {
            c0 += block_phi[by * LB + bz] * block_phi[by * LB + bz];
            cnt0++;
        }
        c0 /= cnt0;

        res.xi_corr = LB / 2.0;
        for (int r = 1; r < LB/2; r++) {
            double cr = 0;
            for (int by = 0; by < LB; by++)
            for (int bz = 0; bz < LB; bz++) {
                cr += block_phi[by * LB + bz] * block_phi[r * LB * LB + by * LB + bz];
            }
            cr /= cnt0;
            if (c0 > 1e-20 && cr / c0 < 1.0/M_E) {
                res.xi_corr = r;
                break;
            }
        }
    }

    return res;
}

/* ============================================================
 * Main
 * ============================================================ */
int main(void) {
    rng_seed((uint64_t)time(NULL) ^ 0xDEADBEEF);

    printf("============================================================\n");
    printf("LK-1b: CONTINUOUS SUBSTRATE MC (TGP coarse-graining)\n");
    printf("============================================================\n");
    printf("H = sum[m0^2/2 s^2 + lambda/4 s^4] - J sum s_i*s_j\n");
    printf("m0^2 = %.2f (use NEGATIVE for broken Z_2)\n", M0_SQ);
    printf("lambda = %.2f, J = %.2f\n", LAMBDA, J_COUP);
    printf("Phi_B = <s^2>_block  (Z_2 invariant order parameter)\n");
    printf("N_therm=%d, N_measure=%d, N_skip=%d\n\n", N_THERM, N_MEASURE, N_SKIP);

    /* Temperature scan at fixed L=32, b=4 */
    double temps[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0};

    printf("============================================================\n");
    printf("Phase 1: Temperature scan (L=32, b=4)\n");
    printf("============================================================\n");
    printf("  %8s %8s %8s %8s %8s %6s\n",
           "T", "<Phi>", "std", "alpha", "|a-2|", "n_fit");

    int L_scan = 32;
    int b_scan = 4;
    double step = 0.5;  /* initial Metropolis step */

    for (int t = 0; t < N_TEMPS; t++) {
        double T = temps[t];
        double beta = 1.0 / T;

        init_field(L_scan);

        /* Adaptive thermalization */
        for (int i = 0; i < N_THERM; i++) {
            double acc = metropolis_sweep(L_scan, beta, step);
            /* Adapt step to target ~40% acceptance */
            if (acc > 0.5) step *= 1.05;
            else if (acc < 0.3) step *= 0.95;
            /* Overrelaxation for better mixing */
            /* overrelax disabled: unstable for m0^2 < 0 near vacuum */
        }

        /* Measure */
        double alpha_acc = 0, phi_acc = 0;
        int n_good = 0;

        for (int m = 0; m < N_MEASURE; m++) {
            for (int s = 0; s < N_SKIP; s++) {
                metropolis_sweep(L_scan, beta, step);
                /* overrelax disabled */
            }
            MeasResult mr = measure_alpha(L_scan, b_scan);
            if (mr.alpha_eff > -100 && mr.n_fit_pts >= 2) {
                alpha_acc += mr.alpha_eff;
                phi_acc += mr.phi_mean;
                n_good++;
            }
        }

        /* Always print Phi diagnostic */
        {
            MeasResult mr_diag = measure_alpha(L_scan, b_scan);
            if (n_good > 0) {
                printf("  %8.2f %8.4f %8.4f %8.3f %8.3f %6d\n",
                       T, phi_acc / n_good, mr_diag.phi_std,
                       alpha_acc / n_good,
                       fabs(alpha_acc / n_good - 2.0),
                       n_good);
            } else {
                printf("  %8.2f %8.4f %8.4f %8s %8s %6d (fit:%d)\n",
                       T, mr_diag.phi_mean, mr_diag.phi_std,
                       "---", "---", 0, mr_diag.n_fit_pts);
            }
        }
    }

    /* Phase 2: Size scaling at best temperature */
    printf("\n============================================================\n");
    printf("Phase 2: Size scaling (T=2.0, multiple L and b)\n");
    printf("============================================================\n");
    printf("  %4s %4s %4s %8s %8s %8s %8s %6s\n",
           "L", "b", "L_B", "<Phi>", "alpha", "|a-2|", "xi", "n_fit");

    double T_best = 2.0;
    double beta_best = 1.0 / T_best;

    for (int li = 0; li < N_LSIZES; li++) {
        int L = LATTICE_SIZES[li];
        if (L > LMAX) continue;

        step = 0.5;
        init_field(L);

        printf("\n  --- L = %d ---\n", L);
        printf("  Thermalizing...\n");
        for (int i = 0; i < N_THERM; i++) {
            double acc = metropolis_sweep(L, beta_best, step);
            if (acc > 0.5) step *= 1.05;
            else if (acc < 0.3) step *= 0.95;
            /* overrelax disabled */
        }

        for (int bi = 0; bi < N_BSIZES; bi++) {
            int b = BLOCK_SIZES[bi];
            int LB = L / b;
            if (LB < 3) continue;

            double alpha_sum = 0, phi_sum = 0, xi_sum = 0;
            int n_good = 0;

            for (int m = 0; m < N_MEASURE; m++) {
                for (int s = 0; s < N_SKIP; s++) {
                    metropolis_sweep(L, beta_best, step);
                    /* overrelax disabled */
                }
                MeasResult mr = measure_alpha(L, b);
                if (mr.alpha_eff > -100 && mr.n_fit_pts >= 2) {
                    alpha_sum += mr.alpha_eff;
                    phi_sum += mr.phi_mean;
                    xi_sum += mr.xi_corr;
                    n_good++;
                }
            }

            if (n_good > 10) {
                double a = alpha_sum / n_good;
                double p = phi_sum / n_good;
                double x = xi_sum / n_good;
                printf("  %4d %4d %4d %8.4f %8.3f %8.3f %8.1f %6d\n",
                       L, b, LB, p, a, fabs(a - 2.0), x, n_good);
            }
        }
    }

    /* Summary */
    printf("\n============================================================\n");
    printf("ANALYSIS\n");
    printf("============================================================\n");
    printf("  TGP prediction: K(Phi) ~ Phi^2 (substrate) or Phi^4 (canonical)\n");
    printf("  alpha_eff should converge to 2 (substrate K_sub=g^2) as L->inf\n");
    printf("  \n");
    printf("  Key differences from Ising test:\n");
    printf("  - Continuous field s_i in R (not binary)\n");
    printf("  - Phi_B = <s^2>_block (not m_B^2)\n");
    printf("  - Self-interaction lambda*s^4 generates non-trivial K(Phi)\n");
    printf("  - Z_2 symmetry: s -> -s => Phi_B invariant\n");
    printf("============================================================\n");

    return 0;
}
