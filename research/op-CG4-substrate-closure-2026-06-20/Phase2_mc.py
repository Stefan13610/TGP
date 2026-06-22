#!/usr/bin/env python3
"""
Phase 2 — MC kampania op-CG4-substrate-closure (2026-06-20). Model M0 (phi^4 GL, Z2/d3).
Anti-Lakatos: forbidden #4 -> walidacja silnika PRZED pomiarem fizycznym.

Akcja sieciowa (standardowy phi^4):
  S = sum_x [ phi_x^2 + lam*(phi_x^2 - 1)^2 - 2*kappa*sum_{mu} phi_x*phi_{x+mu} ]
Aktualizacja: checkerboard Metropolis (nn-only -> dwie podsieci niezalezne).

Tryby:
  validate : skan kappa (pik chi, Binder U4, <phi^2>, xi) -> lokalizacja kappa_c, walidacja Z2
  measure  : przy kappa_c pomiar structure factor pola + operatora zlozonego O=(grad phi)^2
"""
import numpy as np, sys, json, time

def neighbor_sum(phi):
    return (np.roll(phi,1,0)+np.roll(phi,-1,0)
           +np.roll(phi,1,1)+np.roll(phi,-1,1)
           +np.roll(phi,1,2)+np.roll(phi,-1,2))

def make_mask(L):
    idx=np.indices((L,L,L)).sum(0)
    return (idx%2==0)

def metropolis_sweep(phi, kappa, lam, delta, mask, rng):
    for sub in (mask, ~mask):
        h=neighbor_sum(phi)
        prop=phi+delta*(rng.random(phi.shape)*2-1)
        # S_local(phi) = phi^2 + lam(phi^2-1)^2 - 2k phi h
        def Sloc(x): return x*x+lam*(x*x-1)**2-2*kappa*x*h
        dS=Sloc(prop)-Sloc(phi)
        acc=(rng.random(phi.shape)<np.exp(-np.clip(dS,-50,50)))&sub
        phi[acc]=prop[acc]
    return phi

def run(L, kappa, lam, nther, nmeas, nskip, delta, seed, measure=False):
    rng=np.random.default_rng(seed)
    phi=rng.standard_normal((L,L,L))*0.3
    mask=make_mask(L)
    # termalizacja + auto-tune delta do ~50% akceptacji
    for t in range(nther):
        metropolis_sweep(phi,kappa,lam,delta,mask,rng)
        if t<nther//2 and t%20==19:
            # zgrubny tuning delta
            h=neighbor_sum(phi); prop=phi+delta*(rng.random(phi.shape)*2-1)
            dS=(prop*prop+lam*(prop*prop-1)**2-2*kappa*prop*h)-(phi*phi+lam*(phi*phi-1)**2-2*kappa*phi*h)
            a=np.mean(rng.random(phi.shape)<np.exp(-np.clip(dS,-50,50)))
            if a<0.4: delta*=0.9
            elif a>0.6: delta*=1.1
    V=L**3
    mags=[]; m2=[]; m4=[]; phi2=[]
    Skp=[]  # structure factor pola w najmniejszym k
    OkdotO=[]  # operator zlozony O=(grad phi)^2 structure factor w 0 i k_min
    kmin=2*np.pi/L
    for t in range(nmeas):
        for _ in range(nskip):
            metropolis_sweep(phi,kappa,lam,delta,mask,rng)
        m=phi.mean(); mags.append(m); m2.append(m*m); m4.append(m**4)
        phi2.append((phi*phi).mean())
        if measure:
            # structure factor pola: |phi_tilde(k)|^2 dla k wzdluz x
            ft=np.fft.fftn(phi)
            Sk0=(np.abs(ft[0,0,0])**2)/V
            Sk1=(np.abs(ft[1,0,0])**2)/V
            Skp.append((Sk0,Sk1))
            # operator zlozony O_x = sum_mu (phi_{x+mu}-phi_x)^2  (~ (grad phi)^2, gestosc kinetyczna)
            grad2=((np.roll(phi,-1,0)-phi)**2+(np.roll(phi,-1,1)-phi)**2+(np.roll(phi,-1,2)-phi)**2)
            Of=np.fft.fftn(grad2-grad2.mean())  # odjeta srednia (subtrakcja zerowej mody)
            O0=(np.abs(Of[0,0,0])**2)/V
            O1=(np.abs(Of[1,0,0])**2)/V
            O2=(np.abs(Of[2,0,0])**2)/V
            OkdotO.append((O0,O1,O2))
    mags=np.array(mags); m2=np.array(m2); m4=np.array(m4); phi2=np.array(phi2)
    chi=V*(np.mean(mags**2)-np.mean(np.abs(mags))**2)
    U4=1-np.mean(m4)/(3*np.mean(m2)**2)
    out=dict(L=L,kappa=kappa,lam=lam,
             abs_m=float(np.mean(np.abs(mags))),
             chi=float(chi),U4=float(U4),
             phi2=float(np.mean(phi2)),phi2_err=float(np.std(phi2)/np.sqrt(len(phi2))),
             delta=float(delta),nmeas=nmeas)
    if measure:
        Skp=np.array(Skp); Ok=np.array(OkdotO)
        out['Sk_field']=[float(Skp[:,0].mean()),float(Skp[:,1].mean())]
        out['O_struct']=[float(Ok[:,0].mean()),float(Ok[:,1].mean()),float(Ok[:,2].mean())]
        out['kmin']=float(kmin)
    return out

if __name__=="__main__":
    mode=sys.argv[1] if len(sys.argv)>1 else "validate"
    lam=1.0; t0=time.time()
    if mode=="validate":
        L=int(sys.argv[2]) if len(sys.argv)>2 else 12
        nth=int(sys.argv[3]) if len(sys.argv)>3 else 400
        nme=int(sys.argv[4]) if len(sys.argv)>4 else 400
        res=[]
        for kappa in [0.16,0.17,0.18,0.185,0.19,0.195,0.20,0.21]:
            r=run(L,kappa,lam,nth,nme,2,0.8,seed=12345+int(kappa*1000))
            res.append(r)
            print(f"L={L} k={kappa:.3f}  |m|={r['abs_m']:.4f}  chi={r['chi']:.3f}  "
                  f"U4={r['U4']:.4f}  <phi2>={r['phi2']:.4f}  d={r['delta']:.3f}")
        json.dump(res,open(f"Phase2_validate_L{L}.json","w"),indent=1)
        kc=max(res,key=lambda r:r['chi'])['kappa']
        print(f"\nPik chi przy kappa~{kc} (L={L}).  czas={time.time()-t0:.1f}s")
    elif mode=="measure":
        kappa=float(sys.argv[2]);
        nth=int(sys.argv[3]) if len(sys.argv)>3 else 800
        nme=int(sys.argv[4]) if len(sys.argv)>4 else 800
        res=[]
        for L in [16,24,32]:
            r=run(L,kappa,lam,nth,nme,2,0.8,seed=999+L,measure=True)
            res.append(r);
            print(f"L={L} k={kappa}: <phi2>={r['phi2']:.4f}+-{r['phi2_err']:.4f} chi={r['chi']:.2f} "
                  f"U4={r['U4']:.4f} Sk={r['Sk_field']} O={r['O_struct']}")
        json.dump(res,open(f"Phase2_measure_k{kappa}.json","w"),indent=1)
        print(f"czas={time.time()-t0:.1f}s")
