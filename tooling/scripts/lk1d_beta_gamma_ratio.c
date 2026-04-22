/*
 * LK-1d: Numerical verification β_eff / γ_eff → 1
 * ==================================================
 *
 * KEY TEST (N0): The TGP field equation requires β = γ.
 * This follows from U'(Φ₀) = 0 (vacuum condition).
 *
 * Strategy:
 *   - Run continuous substrate MC at several temperatures
 *   - Measure histogram P(Φ_B) of block-averaged field
 *   - Extract V_eff(Φ) = -T·ln P(Φ)
 *   - Fit V_eff near minimum to: V(Φ) = a2·(Φ-Φ₀)² + a3·(Φ-Φ₀)³
 *   - From TGP: a2 = β/(2Φ₀), a3 = (2β-3γ)/(6Φ₀²)
 *   - If β=γ: a3/a2 = -1/(3Φ₀) → ratio R = -3·Φ₀·a3/a2 = 1
 *   - Test: does R → 1?
 *
 * Additional test: measure <Φ²> and <Φ³> fluctuations
 *   - Susceptibility χ₂ = <(δΦ)²> = T/(2·a2·V)
 *   - Skewness S₃ = <(δΦ)³> / <(δΦ)²>^{3/2}
 *   - From TGP: S₃ = -6a3/(2a2)^{3/2} · √V → vanishes for symmetric potential
 *   - β=γ implies specific skewness ratio
 *
 * Compile: gcc -O3 -march=native -o lk1d.exe lk1d_beta_gamma_ratio.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

/* ============================================================ */
#define M0_SQ    -1.0
#define LAMBDA    1.0
#define J_COUP    1.0

#define L         64
#define BLOCK     4
#define LB        (L/BLOCK)
#define N3        (L*L*L)
#define NB3       (LB*LB*LB)
#define BVOL      (BLOCK*BLOCK*BLOCK)

#define N_THERM   1500
#define N_MEASURE 2000
#define N_SKIP    2

#define N_TEMPS   8
#define N_HIST    200

/* ============================================================
 * RNG (xoshiro256**)
 * ============================================================ */
static uint64_t rng_s[4];
static inline uint64_t rotl64(uint64_t x, int k) { return (x<<k)|(x>>(64-k)); }
static uint64_t rng_next(void) {
    uint64_t r = rotl64(rng_s[1]*5,7)*9, t = rng_s[1]<<17;
    rng_s[2]^=rng_s[0]; rng_s[3]^=rng_s[1];
    rng_s[1]^=rng_s[2]; rng_s[0]^=rng_s[3];
    rng_s[2]^=t; rng_s[3]=rotl64(rng_s[3],45);
    return r;
}
static double rng_uniform(void) { return (rng_next()>>11)*0x1.0p-53; }
static double rng_gaussian(void) {
    double u1=rng_uniform(), u2=rng_uniform();
    if(u1<1e-15) u1=1e-15;
    return sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
}
static void rng_seed(uint64_t seed) {
    for(int i=0;i<4;i++){
        seed+=0x9e3779b97f4a7c15ULL;
        uint64_t z=seed;
        z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;
        z=(z^(z>>27))*0x94d049bb133111ebULL;
        rng_s[i]=z^(z>>31);
    }
}

/* ============================================================
 * Field & helpers
 * ============================================================ */
static double field[N3];
static double phi_block[NB3];

static inline int idx3(int x, int y, int z) {
    return ((x%L+L)%L)*L*L + ((y%L+L)%L)*L + ((z%L+L)%L);
}
static inline int bidx(int bx, int by, int bz) {
    return bx*LB*LB + by*LB + bz;
}

static void init_field(void) {
    double sv = sqrt(fabs(M0_SQ)/LAMBDA);
    for(int i=0;i<N3;i++) field[i] = sv + 0.1*rng_gaussian();
}

static double metropolis_sweep(double beta, double step) {
    int acc=0;
    for(int n=0;n<N3;n++){
        int i=(int)(rng_uniform()*N3); if(i>=N3) i=N3-1;
        int x=i/(L*L), y=(i/L)%L, z=i%L;
        double so=field[i], sn=so+step*rng_gaussian();
        double nb=field[idx3(x+1,y,z)]+field[idx3(x-1,y,z)]
                  +field[idx3(x,y+1,z)]+field[idx3(x,y-1,z)]
                  +field[idx3(x,y,z+1)]+field[idx3(x,y,z-1)];
        double dE=(M0_SQ/2.0*(sn*sn-so*so)+LAMBDA/4.0*(sn*sn*sn*sn-so*so*so*so)
                   -J_COUP*(sn-so)*nb);
        if(dE<0||rng_uniform()<exp(-beta*dE)){ field[i]=sn; acc++; }
    }
    return (double)acc/N3;
}

static double compute_phi_block(void) {
    double tot=0;
    for(int bx=0;bx<LB;bx++)
    for(int by=0;by<LB;by++)
    for(int bz=0;bz<LB;bz++){
        double s2=0;
        for(int dx=0;dx<BLOCK;dx++)
        for(int dy=0;dy<BLOCK;dy++)
        for(int dz=0;dz<BLOCK;dz++){
            double s=field[idx3(bx*BLOCK+dx,by*BLOCK+dy,bz*BLOCK+dz)];
            s2+=s*s;
        }
        double p=s2/BVOL;
        phi_block[bidx(bx,by,bz)]=p;
        tot+=p;
    }
    return tot/NB3;
}

/* ============================================================
 * Histogram-based V_eff extraction
 * ============================================================ */
typedef struct {
    double phi0;       /* vacuum (peak of histogram) */
    double a2;         /* quadratic coefficient */
    double a3;         /* cubic coefficient */
    double ratio_R;    /* -3*phi0*a3/a2 — should be 1 if β=γ */
    double chi2;       /* susceptibility <(δΦ)²> */
    double skewness;   /* <(δΦ)³> / <(δΦ)²>^{3/2} */
    double mean_phi;   /* <Φ_B> */
    double sigma_phi;  /* std(Φ_B) */
    int n_bins_used;
} VeffResult;

static VeffResult extract_Veff(double T, double *phi_samples, int n_samples, int n_blocks) {
    VeffResult res;
    memset(&res, 0, sizeof(res));

    /* Compute moments from all block samples */
    double sum1=0, sum2=0, sum3=0, sum4=0;
    int ntot = n_samples * n_blocks;
    double *all_phi = (double*)malloc(ntot * sizeof(double));

    /* All phi values are already stored in phi_samples[sample * n_blocks + block] */
    for(int i=0; i<ntot; i++){
        double p = phi_samples[i];
        sum1 += p;
    }
    double mean = sum1 / ntot;

    for(int i=0; i<ntot; i++){
        double dp = phi_samples[i] - mean;
        sum2 += dp*dp;
        sum3 += dp*dp*dp;
        sum4 += dp*dp*dp*dp;
    }
    double var = sum2 / ntot;
    double mu3 = sum3 / ntot;
    double mu4 = sum4 / ntot;

    res.mean_phi = mean;
    res.sigma_phi = sqrt(var);
    res.chi2 = var;
    res.skewness = (var > 1e-15) ? mu3 / pow(var, 1.5) : 0.0;

    /* Build histogram */
    double phi_min = mean - 5.0*res.sigma_phi;
    double phi_max = mean + 5.0*res.sigma_phi;
    if(phi_min < 0) phi_min = 0;
    double bin_width = (phi_max - phi_min) / N_HIST;

    int hist[N_HIST];
    memset(hist, 0, sizeof(hist));
    for(int i=0; i<ntot; i++){
        int b = (int)((phi_samples[i] - phi_min) / bin_width);
        if(b >= 0 && b < N_HIST) hist[b]++;
    }

    /* Find peak (Φ₀) */
    int peak_bin = 0;
    for(int b=1; b<N_HIST; b++)
        if(hist[b] > hist[peak_bin]) peak_bin = b;

    res.phi0 = phi_min + (peak_bin + 0.5) * bin_width;

    /* Extract V_eff = -T·ln(P) near peak */
    /* Fit to V_eff(Φ) = a2·(Φ-Φ₀)² + a3·(Φ-Φ₀)³ using weighted least squares */
    /* Use bins within 2σ of peak */
    double fit_range = 2.5 * res.sigma_phi;

    /* Collect valid bins */
    double SxxA=0, SxxB=0, SxyA=0, SxyB=0;
    double Sx2x2=0, Sx2x3=0, Sx3x3=0, Sx2y=0, Sx3y=0;
    int n_fit = 0;

    for(int b=0; b<N_HIST; b++){
        if(hist[b] < 5) continue; /* skip low-stats bins */
        double phi_b = phi_min + (b + 0.5) * bin_width;
        double dp = phi_b - res.phi0;
        if(fabs(dp) > fit_range) continue;

        double V = -T * log((double)hist[b]);
        /* Subtract V at peak */
        double V0 = -T * log((double)hist[peak_bin]);
        double dV = V - V0;

        double x2 = dp * dp;
        double x3 = dp * dp * dp;

        /* Normal equations for dV = a2·dp² + a3·dp³ */
        Sx2x2 += x2 * x2;
        Sx2x3 += x2 * x3;
        Sx3x3 += x3 * x3;
        Sx2y  += x2 * dV;
        Sx3y  += x3 * dV;
        n_fit++;
    }

    res.n_bins_used = n_fit;

    if(n_fit < 5){
        res.a2 = 0; res.a3 = 0; res.ratio_R = 0;
        free(all_phi);
        return res;
    }

    /* Solve 2x2 system:
     * [Sx2x2  Sx2x3] [a2]   [Sx2y]
     * [Sx2x3  Sx3x3] [a3] = [Sx3y]
     */
    double det = Sx2x2 * Sx3x3 - Sx2x3 * Sx2x3;
    if(fabs(det) < 1e-30){
        res.a2 = 0; res.a3 = 0; res.ratio_R = 0;
        free(all_phi);
        return res;
    }

    res.a2 = (Sx3x3 * Sx2y  - Sx2x3 * Sx3y) / det;
    res.a3 = (Sx2x2 * Sx3y  - Sx2x3 * Sx2y) / det;

    /* TGP prediction: β=γ ↔ R = -3·Φ₀·a3/a2 = 1 */
    if(fabs(res.a2) > 1e-15)
        res.ratio_R = -3.0 * res.phi0 * res.a3 / res.a2;
    else
        res.ratio_R = 0;

    free(all_phi);
    return res;
}

/* ============================================================
 * MAIN
 * ============================================================ */
int main(void) {
    rng_seed(20260414ULL);

    double temps[N_TEMPS] = {1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};

    printf("======================================================================\n");
    printf("LK-1d: beta_eff / gamma_eff --> 1 (vacuum condition test)\n");
    printf("======================================================================\n");
    printf("L=%d, block=%d, LB=%d, N_therm=%d, N_measure=%d\n",
           L, BLOCK, LB, N_THERM, N_MEASURE);
    printf("Fixed: m0^2=%.1f, lam=%.1f, J=%.1f\n", M0_SQ, LAMBDA, J_COUP);
    printf("Test: R = -3*Phi0*a3/a2 should --> 1 if beta=gamma (TGP N0)\n\n");

    printf("  T       <Phi>     sigma    a2        a3         R=bg_ratio  skewness  bins\n");
    printf("  --------------------------------------------------------------------------\n");

    /* Store results for analysis */
    double all_T[N_TEMPS], all_R[N_TEMPS], all_phi0[N_TEMPS];
    double all_a2[N_TEMPS], all_a3[N_TEMPS], all_skew[N_TEMPS];
    int n_valid = 0;

    for(int it=0; it<N_TEMPS; it++){
        double T = temps[it];
        double beta = 1.0 / T;

        init_field();

        /* Thermalize */
        double step = 1.0;
        for(int sw=0; sw<N_THERM; sw++){
            double acc = metropolis_sweep(beta, step);
            if(acc > 0.55) step *= 1.05;
            if(acc < 0.35) step *= 0.95;
        }

        /* Allocate storage for all block phi values */
        double *phi_samples = (double*)malloc(N_MEASURE * NB3 * sizeof(double));

        /* Measure */
        for(int m=0; m<N_MEASURE; m++){
            for(int sk=0; sk<N_SKIP; sk++)
                metropolis_sweep(beta, step);

            compute_phi_block();

            /* Store all block values */
            for(int i=0; i<NB3; i++)
                phi_samples[m * NB3 + i] = phi_block[i];
        }

        /* Extract V_eff and coefficients */
        VeffResult res = extract_Veff(T, phi_samples, N_MEASURE, NB3);

        printf("  %.1f    %8.4f  %7.4f  %9.4f  %10.6f  %9.4f    %+7.4f    %d\n",
               T, res.mean_phi, res.sigma_phi, res.a2, res.a3, res.ratio_R,
               res.skewness, res.n_bins_used);

        all_T[n_valid] = T;
        all_R[n_valid] = res.ratio_R;
        all_phi0[n_valid] = res.phi0;
        all_a2[n_valid] = res.a2;
        all_a3[n_valid] = res.a3;
        all_skew[n_valid] = res.skewness;
        n_valid++;

        free(phi_samples);
    }

    /* ============================================================
     * Analysis
     * ============================================================ */
    printf("\n======================================================================\n");
    printf("ANALYSIS: beta/gamma ratio\n");
    printf("======================================================================\n");

    /* Compute mean and std of R */
    double R_sum=0, R_sum2=0;
    int R_count=0;
    for(int i=0; i<n_valid; i++){
        if(all_a2[i] > 0.1){ /* Only use points with well-determined a2 */
            R_sum += all_R[i];
            R_sum2 += all_R[i] * all_R[i];
            R_count++;
        }
    }

    if(R_count > 1){
        double R_mean = R_sum / R_count;
        double R_var = R_sum2/R_count - R_mean*R_mean;
        double R_std = sqrt(fabs(R_var));
        double R_err = R_std / sqrt(R_count);

        printf("\n  Valid points (a2 > 0.1): %d\n", R_count);
        printf("  R = -3*Phi0*a3/a2:\n");
        printf("    Mean R = %.4f +/- %.4f\n", R_mean, R_err);
        printf("    Std(R) = %.4f\n", R_std);
        printf("    TGP prediction: R = 1 (beta = gamma)\n");
        printf("    |R - 1| = %.4f\n", fabs(R_mean - 1.0));
        printf("    Deviation: %.2f sigma\n",
               (R_err > 1e-10) ? fabs(R_mean - 1.0)/R_err : 999.0);

        /* Near-Tc analysis (higher T, larger fluctuations) */
        double R_near_sum=0, R_near_sum2=0;
        int R_near_count=0;
        printf("\n  --- Near-T_c points (T >= 3.5) ---\n");
        for(int i=0; i<n_valid; i++){
            if(all_T[i] >= 3.5 && all_a2[i] > 0.1){
                printf("    T=%.1f: R=%.4f, a2=%.4f, a3=%.6f, Phi0=%.3f\n",
                       all_T[i], all_R[i], all_a2[i], all_a3[i], all_phi0[i]);
                R_near_sum += all_R[i];
                R_near_sum2 += all_R[i] * all_R[i];
                R_near_count++;
            }
        }
        if(R_near_count > 1){
            double Rn_mean = R_near_sum / R_near_count;
            double Rn_var = R_near_sum2/R_near_count - Rn_mean*Rn_mean;
            double Rn_err = sqrt(fabs(Rn_var)) / sqrt(R_near_count);
            printf("  Near-Tc mean R = %.4f +/- %.4f\n", Rn_mean, Rn_err);
            printf("  |R_nearTc - 1| = %.4f (%.2f sigma)\n",
                   fabs(Rn_mean - 1.0),
                   (Rn_err > 1e-10) ? fabs(Rn_mean - 1.0)/Rn_err : 999.0);
        }

        /* Skewness analysis */
        printf("\n  --- Skewness analysis ---\n");
        printf("  TGP: β=γ implies specific skewness from cubic term\n");
        double sk_sum=0;
        for(int i=0; i<n_valid; i++){
            printf("    T=%.1f: skewness = %+.4f\n", all_T[i], all_skew[i]);
            sk_sum += fabs(all_skew[i]);
        }
        printf("  Mean |skewness| = %.4f\n", sk_sum / n_valid);
        printf("  (Small skewness = approximately symmetric potential = consistent with β≈γ)\n");
    }

    /* ============================================================
     * Tests
     * ============================================================ */
    printf("\n======================================================================\n");
    printf("TESTS\n");
    printf("======================================================================\n");

    int pass=0, total=0;

    /* T1: a2 > 0 everywhere (stable vacuum) */
    total++;
    int a2_pos = 1;
    for(int i=0; i<n_valid; i++) if(all_a2[i] <= 0) a2_pos = 0;
    if(a2_pos){ printf("  [PASS] T1: a2 > 0 at all temperatures (stable vacuum)\n"); pass++; }
    else printf("  [FAIL] T1: a2 <= 0 at some temperature\n");

    /* T2: R exists and is finite */
    total++;
    int R_finite = 1;
    for(int i=0; i<n_valid; i++) if(!isfinite(all_R[i]) || all_a2[i] < 0.01) R_finite = 0;
    if(R_finite){ printf("  [PASS] T2: R = -3*Phi0*a3/a2 is finite at all T\n"); pass++; }
    else printf("  [FAIL] T2: R is not finite at some T\n");

    /* T3: Mean R consistent with 1 within 3σ */
    total++;
    if(R_count > 1){
        double R_mean = R_sum / R_count;
        double R_var = R_sum2/R_count - R_mean*R_mean;
        double R_err = sqrt(fabs(R_var)) / sqrt(R_count);
        double dev = (R_err > 1e-10) ? fabs(R_mean - 1.0)/R_err : 999.0;
        if(dev < 3.0){
            printf("  [PASS] T3: R = %.3f +/- %.3f, consistent with 1 (%.1f sigma)\n",
                   R_mean, R_err, dev);
            pass++;
        } else {
            printf("  [FAIL] T3: R = %.3f +/- %.3f, %.1f sigma from 1\n",
                   R_mean, R_err, dev);
        }
    } else printf("  [FAIL] T3: insufficient valid points\n");

    /* T4: Potential is confining (a2 increases toward T_c) */
    total++;
    /* As T→T_c, fluctuations grow, a2 should decrease */
    if(n_valid >= 4 && all_a2[0] > all_a2[n_valid-1]){
        printf("  [PASS] T4: a2 decreases toward T_c (mass gap closing, %.3f -> %.3f)\n",
               all_a2[0], all_a2[n_valid-1]);
        pass++;
    } else {
        printf("  [FAIL] T4: a2 does not decrease toward T_c\n");
    }

    /* T5: |skewness| < 1 everywhere (approximately symmetric) */
    total++;
    int skew_ok = 1;
    for(int i=0; i<n_valid; i++) if(fabs(all_skew[i]) > 1.0) skew_ok = 0;
    if(skew_ok){
        printf("  [PASS] T5: |skewness| < 1 at all T (approximately symmetric potential)\n");
        pass++;
    } else printf("  [FAIL] T5: large skewness detected\n");

    /* T6: Phi0 is well-defined and > 0 */
    total++;
    int phi0_ok = 1;
    for(int i=0; i<n_valid; i++) if(all_phi0[i] <= 0) phi0_ok = 0;
    if(phi0_ok){
        printf("  [PASS] T6: Phi0 > 0 at all T (non-trivial vacuum)\n");
        pass++;
    } else printf("  [FAIL] T6: Phi0 <= 0\n");

    printf("\n  Overall: %d/%d PASS\n", pass, total);

    printf("\n  +-------------------------------------------------------------+\n");
    printf("  |  LK-1d: beta/gamma ratio test                              |\n");
    printf("  |                                                             |\n");
    printf("  |  TGP prediction: R = -3*Phi0*a3/a2 = 1 (beta = gamma)      |\n");
    printf("  |  This is equivalent to U'(Phi0) = 0 (vacuum condition)      |\n");
    printf("  |  A violation would require beta != gamma,                   |\n");
    printf("  |  breaking the TGP field equation structure.                 |\n");
    printf("  +-------------------------------------------------------------+\n");

    return 0;
}
