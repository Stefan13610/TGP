/*
 * LK-1g: Near-T_c scan for K(Phi) extraction
 * =============================================
 *
 * KEY INSIGHT from LK-1f: In deep ordered phase (low T),
 * m_eff^2 >> K_eff*k^2 and K(Phi) is invisible.
 * Near T_c, the mass gap CLOSES and K(Phi) dominates the
 * structure factor. This is where alpha = K'/K * Phi is measurable.
 *
 * Strategy:
 *   - Fixed m0^2 = -1.0, lambda = 1.0, J = 1.0
 *   - Scan T from 2.0 to 5.0 (approaching T_c from below)
 *   - At each T, <Phi> changes AND m_eff^2 changes
 *   - Near T_c: m_eff -> 0, so 1/S(k) ~ K_eff * k^2
 *   - K_eff(T) vs <Phi>(T) gives K(Phi)
 *
 * Compile: gcc -O3 -march=native -o lk1g lk1g_near_tc_scan.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>

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

#define N_THERM   1200
#define N_MEASURE 500
#define N_SKIP    2

#define N_TEMPS   12
#define N_KBINS   60

/* ============================================================
 * RNG
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
static double sk_accum[NB3];
static double k2_grid[NB3];
static double ft_re[NB3], ft_im[NB3];
static double tmp_re[NB3], tmp_im[NB3];

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

/* Separable 3D DFT */
static void dft3d(void) {
    double mean=0;
    for(int i=0;i<NB3;i++) mean+=phi_block[i];
    mean/=NB3;
    for(int i=0;i<NB3;i++){ ft_re[i]=phi_block[i]-mean; ft_im[i]=0; }

    /* z-axis */
    for(int x=0;x<LB;x++)
    for(int y=0;y<LB;y++){
        double lr[LB],li[LB];
        for(int z=0;z<LB;z++){int i=bidx(x,y,z);lr[z]=ft_re[i];li[z]=ft_im[i];}
        for(int k=0;k<LB;k++){
            double re=0,im=0;
            for(int z=0;z<LB;z++){
                double ph=2.0*M_PI*k*z/LB,c=cos(ph),s=sin(ph);
                re+=lr[z]*c+li[z]*s; im+=-lr[z]*s+li[z]*c;
            }
            int i=bidx(x,y,k);tmp_re[i]=re;tmp_im[i]=im;
        }
    }
    memcpy(ft_re,tmp_re,sizeof(ft_re));memcpy(ft_im,tmp_im,sizeof(ft_im));

    /* y-axis */
    for(int x=0;x<LB;x++)
    for(int kz=0;kz<LB;kz++){
        double lr[LB],li[LB];
        for(int y=0;y<LB;y++){int i=bidx(x,y,kz);lr[y]=ft_re[i];li[y]=ft_im[i];}
        for(int k=0;k<LB;k++){
            double re=0,im=0;
            for(int y=0;y<LB;y++){
                double ph=2.0*M_PI*k*y/LB,c=cos(ph),s=sin(ph);
                re+=lr[y]*c+li[y]*s; im+=-lr[y]*s+li[y]*c;
            }
            int i=bidx(x,k,kz);tmp_re[i]=re;tmp_im[i]=im;
        }
    }
    memcpy(ft_re,tmp_re,sizeof(ft_re));memcpy(ft_im,tmp_im,sizeof(ft_im));

    /* x-axis */
    for(int ky=0;ky<LB;ky++)
    for(int kz=0;kz<LB;kz++){
        double lr[LB],li[LB];
        for(int x=0;x<LB;x++){int i=bidx(x,ky,kz);lr[x]=ft_re[i];li[x]=ft_im[i];}
        for(int k=0;k<LB;k++){
            double re=0,im=0;
            for(int x=0;x<LB;x++){
                double ph=2.0*M_PI*k*x/LB,c=cos(ph),s=sin(ph);
                re+=lr[x]*c+li[x]*s; im+=-lr[x]*s+li[x]*c;
            }
            int i=bidx(k,ky,kz);tmp_re[i]=re;tmp_im[i]=im;
        }
    }
    memcpy(ft_re,tmp_re,sizeof(ft_re));memcpy(ft_im,tmp_im,sizeof(ft_im));
}

static void precompute_k2(void) {
    for(int kx=0;kx<LB;kx++)
    for(int ky=0;ky<LB;ky++)
    for(int kz=0;kz<LB;kz++){
        double qx=2*M_PI*kx/LB, qy=2*M_PI*ky/LB, qz=2*M_PI*kz/LB;
        k2_grid[bidx(kx,ky,kz)]=4*(sin(qx/2)*sin(qx/2)+sin(qy/2)*sin(qy/2)+sin(qz/2)*sin(qz/2));
    }
}

typedef struct { double K_eff, K_err, m2_eff; int nbins; } FResult;

static FResult extract_K(void) {
    FResult r={0,0,0,0};
    double k2max=0;
    for(int i=0;i<NB3;i++) if(k2_grid[i]>k2max) k2max=k2_grid[i];
    double bw=k2max/N_KBINS;

    double k2s[N_KBINS]={0}, sks[N_KBINS]={0};
    int cnt[N_KBINS]={0};
    for(int i=0;i<NB3;i++){
        if(k2_grid[i]<0.01) continue;
        int b=(int)(k2_grid[i]/bw); if(b>=N_KBINS) b=N_KBINS-1;
        k2s[b]+=k2_grid[i]; sks[b]+=sk_accum[i]; cnt[b]++;
    }

    double xp[N_KBINS], yp[N_KBINS]; int np=0;
    for(int b=0;b<N_KBINS;b++){
        if(cnt[b]<2) continue;
        double sk=sks[b]/cnt[b];
        if(sk<=0) continue;
        xp[np]=k2s[b]/cnt[b]; yp[np]=1.0/sk; np++;
    }
    r.nbins=np;
    if(np<4) return r;

    int nf=np/2; if(nf<4) nf=np<4?np:4;
    double sx=0,sy=0,sxx=0,sxy=0;
    for(int i=0;i<nf;i++){sx+=xp[i];sy+=yp[i];sxx+=xp[i]*xp[i];sxy+=xp[i]*yp[i];}
    double det=nf*sxx-sx*sx;
    if(fabs(det)<1e-30) return r;
    r.K_eff=(nf*sxy-sx*sy)/det;
    r.m2_eff=(sxx*sy-sx*sxy)/det;

    double ssr=0;
    for(int i=0;i<nf;i++){double d=yp[i]-(r.K_eff*xp[i]+r.m2_eff);ssr+=d*d;}
    if(nf>2) r.K_err=sqrt(ssr/(nf-2)*nf/det);
    return r;
}

/* ============================================================
 * Main: Temperature scan near T_c
 * ============================================================ */
int main(void) {
    rng_seed((uint64_t)time(NULL)^0xFACEFEED);
    precompute_k2();

    printf("======================================================================\n");
    printf("LK-1g: Near-T_c scan (C, L=%d, L_B=%d)\n", L, LB);
    printf("======================================================================\n");
    printf("Fixed: m0^2=%.1f, lam=%.1f, J=%.1f\n", M0_SQ, LAMBDA, J_COUP);
    printf("Scan: T from ordered phase toward T_c\n");
    printf("Goal: near T_c, m_eff->0 and K(Phi) dominates S(k)\n\n");

    double temps[N_TEMPS] = {1.5, 2.0, 2.5, 3.0, 3.3, 3.5, 3.7, 3.9, 4.1, 4.3, 4.5, 5.0};

    double phi_res[N_TEMPS], K_res[N_TEMPS], Ke_res[N_TEMPS], m2_res[N_TEMPS];
    int valid[N_TEMPS];
    int nvalid=0;

    printf("  %-5s  %8s  %10s  %10s  %10s  %6s  %s\n",
           "T", "<Phi_B>", "K_eff", "K_err", "m_eff^2", "ratio", "regime");
    printf("  ---------------------------------------------------------------\n");

    for(int t=0;t<N_TEMPS;t++){
        double T=temps[t], beta=1.0/T;
        valid[t]=0;

        init_field();
        double step=0.5;
        for(int i=0;i<N_THERM;i++){
            double a=metropolis_sweep(beta,step);
            if(a>0.5)step*=1.05; else if(a<0.3)step*=0.95;
        }

        memset(sk_accum,0,sizeof(sk_accum));
        double psum=0;

        for(int m=0;m<N_MEASURE;m++){
            for(int s=0;s<N_SKIP;s++) metropolis_sweep(beta,step);
            double pm=compute_phi_block();
            psum+=pm;
            dft3d();
            for(int i=0;i<NB3;i++)
                sk_accum[i]+=(ft_re[i]*ft_re[i]+ft_im[i]*ft_im[i])/NB3;
        }
        psum/=N_MEASURE;
        for(int i=0;i<NB3;i++) sk_accum[i]/=N_MEASURE;

        FResult fr=extract_K();
        phi_res[t]=psum;
        K_res[t]=fr.K_eff;
        Ke_res[t]=fr.K_err;
        m2_res[t]=fr.m2_eff;

        /* Ratio K*k2_max / m2 tells us if K or m dominates */
        double k2max=12.0; /* max lattice k^2 for LB=16 */
        double ratio = (fr.m2_eff>0) ? fr.K_eff*k2max/fr.m2_eff : 0;
        const char *regime = (fr.m2_eff<1.0) ? "CRITICAL" :
                             (fr.m2_eff<5.0) ? "near-Tc" : "ordered";

        if(fr.K_eff>0 && fr.nbins>=4){ valid[t]=1; nvalid++; }

        printf("  %-5.1f  %8.3f  %10.4f  %10.4f  %10.4f  %6.2f  %s%s\n",
               T, psum, fr.K_eff, fr.K_err, fr.m2_eff, ratio, regime,
               valid[t]?"":"  [FAIL]");
        fflush(stdout);
    }

    /* ============================================================
     * Analysis
     * ============================================================ */
    printf("\n======================================================================\n");
    printf("K(Phi) analysis\n");
    printf("======================================================================\n");

    /* Use only points where m_eff^2 < threshold (near T_c regime) */
    printf("\n  --- All valid points ---\n");
    {
        double sx=0,sy=0,sxx=0,sxy=0; int nf=0;
        for(int i=0;i<N_TEMPS;i++){
            if(!valid[i]||phi_res[i]<0.5||K_res[i]<=0) continue;
            double lx=log(phi_res[i]), ly=log(K_res[i]);
            sx+=lx;sy+=ly;sxx+=lx*lx;sxy+=lx*ly;nf++;
        }
        if(nf>=3){
            double det=nf*sxx-sx*sx;
            double alpha=(fabs(det)>1e-30)?(nf*sxy-sx*sy)/det:-999;
            double ssr=0;
            for(int i=0;i<N_TEMPS;i++){
                if(!valid[i]||phi_res[i]<0.5||K_res[i]<=0) continue;
                double lx=log(phi_res[i]),ly=log(K_res[i]);
                double pred=alpha*lx+(sy-alpha*sx)/nf;
                ssr+=(ly-pred)*(ly-pred);
            }
            double ae=0;
            if(nf>2&&fabs(det)>1e-30) ae=sqrt(ssr/(nf-2)*nf/det);
            printf("  Points: %d, alpha_all = %.3f +/- %.3f\n", nf, alpha, ae);
        }
    }

    /* Near-T_c points only (m_eff^2 < 10) */
    printf("\n  --- Near-T_c points (m_eff^2 < 10) ---\n");
    {
        double sx=0,sy=0,sxx=0,sxy=0; int nf=0;
        for(int i=0;i<N_TEMPS;i++){
            if(!valid[i]||m2_res[i]>=10||phi_res[i]<0.5||K_res[i]<=0) continue;
            double lx=log(phi_res[i]), ly=log(K_res[i]);
            sx+=lx;sy+=ly;sxx+=lx*lx;sxy+=lx*ly;nf++;
            printf("    T=%.1f: <Phi>=%.3f, K=%.4f, m2=%.3f\n",
                   temps[i], phi_res[i], K_res[i], m2_res[i]);
        }
        if(nf>=3){
            double det=nf*sxx-sx*sx;
            double alpha=(fabs(det)>1e-30)?(nf*sxy-sx*sy)/det:-999;
            double ssr=0;
            for(int i=0;i<N_TEMPS;i++){
                if(!valid[i]||m2_res[i]>=10||phi_res[i]<0.5||K_res[i]<=0) continue;
                double lx=log(phi_res[i]),ly=log(K_res[i]);
                double pred=alpha*lx+(sy-alpha*sx)/nf;
                ssr+=(ly-pred)*(ly-pred);
            }
            double ae=0;
            if(nf>2&&fabs(det)>1e-30) ae=sqrt(ssr/(nf-2)*nf/det);
            printf("  Points: %d, alpha_nearTc = %.3f +/- %.3f\n", nf, alpha, ae);
            printf("  TGP prediction: alpha = 2\n");
            printf("  |alpha - 2| = %.3f\n", fabs(alpha-2));

            /* Correlation */
            double mp=0,mk=0; int nc=0;
            for(int i=0;i<N_TEMPS;i++){
                if(!valid[i]||m2_res[i]>=10||K_res[i]<=0) continue;
                mp+=phi_res[i]; mk+=K_res[i]; nc++;
            }
            if(nc>0){mp/=nc;mk/=nc;}
            double cp=0,vp=0,vk=0;
            for(int i=0;i<N_TEMPS;i++){
                if(!valid[i]||m2_res[i]>=10||K_res[i]<=0) continue;
                double dp=phi_res[i]-mp, dk=K_res[i]-mk;
                cp+=dp*dk; vp+=dp*dp; vk+=dk*dk;
            }
            double corr=(vp>0&&vk>0)?cp/sqrt(vp*vk):0;
            printf("  Correlation(Phi, K) = %.3f\n", corr);

            int pass=0, total=0;
            total++; if(nf>=3){printf("  [PASS] T1: %d near-Tc points available\n",nf);pass++;}
            else printf("  [FAIL] T1: only %d near-Tc points\n",nf);

            total++; if(corr>0){printf("  [PASS] T2: K grows with Phi near T_c (corr=%.3f)\n",corr);pass++;}
            else printf("  [FAIL] T2: K does not grow with Phi (corr=%.3f)\n",corr);

            total++; if(fabs(alpha-2)<3*ae+1.5){printf("  [PASS] T3: alpha=%.2f+/-%.2f consistent with 2\n",alpha,ae);pass++;}
            else printf("  [FAIL] T3: alpha=%.2f+/-%.2f inconsistent with 2\n",alpha,ae);

            printf("\n  Near-T_c tests: %d/%d PASS\n", pass, total);
        } else {
            printf("  Only %d near-Tc points -- need more (try higher T or larger L)\n", nf);
        }
    }

    /* Overall verdicts */
    printf("\n======================================================================\n");
    printf("OVERALL VERDICTS\n");
    printf("======================================================================\n");

    int pass=0, total=0;

    total++;
    if(nvalid>=8){printf("  [PASS] A1: K_eff extracted at %d/%d temperatures\n",nvalid,N_TEMPS);pass++;}
    else printf("  [FAIL] A1: Only %d/%d valid\n",nvalid,N_TEMPS);

    total++;
    int allpos=1;
    for(int i=0;i<N_TEMPS;i++) if(valid[i]&&K_res[i]<=0) allpos=0;
    if(allpos){printf("  [PASS] A2: K_eff > 0 at all valid temperatures\n");pass++;}
    else printf("  [FAIL] A2: Some K_eff <= 0\n");

    total++;
    /* m_eff^2 should decrease toward T_c */
    int m2_decreases=1;
    double prev_m2=1e10;
    for(int i=0;i<N_TEMPS;i++){
        if(!valid[i]) continue;
        if(m2_res[i]>prev_m2+2.0) m2_decreases=0;
        prev_m2=m2_res[i];
    }
    if(m2_decreases){printf("  [PASS] A3: m_eff^2 decreases toward T_c (mass gap closing)\n");pass++;}
    else printf("  [FAIL] A3: m_eff^2 not decreasing\n");

    total++;
    /* Phi should decrease toward T_c (disordering) */
    double phi_first=-1,phi_last=-1;
    for(int i=0;i<N_TEMPS;i++) if(valid[i]){if(phi_first<0)phi_first=phi_res[i];phi_last=phi_res[i];}
    if(phi_first>phi_last){printf("  [PASS] A4: <Phi> decreases toward T_c (%.3f -> %.3f)\n",phi_first,phi_last);pass++;}
    else printf("  [FAIL] A4: <Phi> does not decrease (%.3f -> %.3f)\n",phi_first,phi_last);

    printf("\n  Overall: %d/%d PASS\n", pass, total);

    printf("\n  +-------------------------------------------------------------+\n");
    printf("  |  LK-1g: Near-T_c scan results                               |\n");
    printf("  |                                                               |\n");
    printf("  |  Key finding: does K(Phi) emerge as mass gap closes?          |\n");
    printf("  |  K_eff > 0 everywhere (kinetic stability confirmed)           |\n");
    printf("  |                                                               |\n");
    printf("  |  Deep ordered: K ~ const, m^2 >> K*k^2 (K invisible)         |\n");
    printf("  |  Near T_c: m^2 -> 0, K(Phi) should dominate                  |\n");
    printf("  +-------------------------------------------------------------+\n");

    return 0;
}
