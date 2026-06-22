#!/usr/bin/env python3
"""
Phase 2b — czysta ekstrakcja sztywnosci pola Z_R z structure factor (M0, kappa_c).
Cel: cross-check c*>0 (rozwiazanie red-flag CG-34; tam Z~2.46-2.61) na substracie M0.
1/S(k) = const + Z_R * khat^2 ; khat^2 = 2(1-cos(2 pi n / L)).  Z_R = slope.
Anti-Lakatos: metoda structure-factor (forbidden #1: NIE <|grad phi|^2>).
"""
import numpy as np, sys, json, time
from Phase2_mc import metropolis_sweep, make_mask, neighbor_sum

def stiffness(L, kappa, lam, nther, nmeas, nskip, delta, seed, nk=6):
    rng=np.random.default_rng(seed); phi=rng.standard_normal((L,L,L))*0.3
    mask=make_mask(L); V=L**3
    for t in range(nther): metropolis_sweep(phi,kappa,lam,delta,mask,rng)
    Sacc=np.zeros(nk+1)
    for t in range(nmeas):
        for _ in range(nskip): metropolis_sweep(phi,kappa,lam,delta,mask,rng)
        ft=np.fft.fftn(phi)
        # usrednij po 3 osiach dla n=0..nk
        for n in range(nk+1):
            Sacc[n]+=(np.abs(ft[n,0,0])**2+np.abs(ft[0,n,0])**2+np.abs(ft[0,0,n])**2)/(3*V)
    S=Sacc/nmeas
    n=np.arange(1,nk+1)
    khat2=2*(1-np.cos(2*np.pi*n/L))
    invS=1.0/S[1:]
    # fit liniowy 1/S = a + Z * khat2
    A=np.vstack([np.ones_like(khat2),khat2]).T
    coef,res,_,_=np.linalg.lstsq(A,invS,rcond=None)
    a,Z=coef
    pred=A@coef; ss=1-np.sum((invS-pred)**2)/np.sum((invS-invS.mean())**2)
    return dict(L=L,kappa=kappa,Z_R=float(Z),intercept=float(a),R2=float(ss),
                S=[float(x) for x in S],khat2=[float(x) for x in khat2])

if __name__=="__main__":
    lam=1.0; kappa=float(sys.argv[1]) if len(sys.argv)>1 else 0.190
    t0=time.time(); res=[]
    for L in [16,24,32]:
        r=stiffness(L,kappa,lam,1200,1200,2,0.8,seed=7000+L)
        res.append(r)
        print(f"L={L} k={kappa}: Z_R={r['Z_R']:.3f}  intercept={r['intercept']:.4f}  R2={r['R2']:.4f}")
    json.dump(res,open(f"Phase2_stiffness_k{kappa}.json","w"),indent=1)
    Zs=[r['Z_R'] for r in res]
    print(f"\nZ_R(L): {[round(z,3) for z in Zs]}  -> spread {max(Zs)/min(Zs):.2f}x")
    print(f"CG-34 ref: Z~2.46-2.61.  c*>0 {'POTWIERDZONE' if min(Zs)>0 else 'NIE'}.  czas={time.time()-t0:.1f}s")
