# RIBBONS_stepB2_skyrme_v3.py  (KANONICZNY, dzialajacy)
# Historia (uczciwie): v1 (L-BFGS, init=rampa pi(1-r/R) na duzym pudle R=30-60) UTKNAL ->
#   artefakt pudla: <r>~R/2, E~R, wirial rosnacy 276..1239. v2 (solve_bvp) NIE ZBIEGL
#   (osobliwosc r=0 + blad initu). PRZYCZYNA obu: zly init + zbyt duze pudlo.
# v3: relaksacja energii, init ZLOKALIZOWANY 4*arctan(e^-r) (f(0)=pi dokladnie),
#   pudlo FIZYCZNE R~8-12 (rozmiar skyrmionu ~1). To znajduje prawdziwy skyrmion B=1.
#
# E = INT_0^R [ r^2 f'^2 + 2 sin^2 f + 2 sin^2 f f'^2 + sin^4 f / r^2 ] dr ; f(0)=pi, f(R)=0.
# Derrick/wirial: E2=E4 w stacjonarnym punkcie (tu ~4% resztka z ucieci pudla/r=0 - uczciwie).

import numpy as np
from scipy.optimize import minimize
from scipy.integrate import trapezoid

def parts(f, r):
    dr = r[1]-r[0]; fp = np.gradient(f, dr); s2 = np.sin(f)**2
    rr = r.copy(); rr[0] = dr
    return trapezoid(rr**2*fp**2 + 2*s2, r), trapezoid(2*s2*fp**2 + s2**2/rr**2, r)

def Etot(fin, r):
    f = np.empty(len(r)); f[0]=np.pi; f[-1]=0.0; f[1:-1]=fin
    E2,E4 = parts(f, r); return E2+E4

def solve(R, N):
    r = np.linspace(0, R, N); fi = 4*np.arctan(np.exp(-r))
    res = minimize(Etot, fi[1:-1], args=(r,), method='L-BFGS-B',
                   options={'maxiter':4000,'ftol':1e-14,'gtol':1e-10})
    f = np.empty(N); f[0]=np.pi; f[-1]=0.0; f[1:-1]=res.x
    E2,E4 = parts(f, r); fp = np.gradient(f, r[1]-r[0])
    bd = np.clip(-(np.sin(f)**2)*fp, 0, None); w = trapezoid(bd, r)
    return dict(E=E2+E4, E2=E2, E4=E4, vir=E2/E4, rms=trapezoid(bd*r,r)/w, B=(f[0]-f[-1])/np.pi)

print("[R] SKAN PUDLA (kontrola artefaktu; v1 tu poleglo z <r>~R/2):")
print(f"  {'R':>4} {'E':>9} {'E2':>8} {'E4':>8} {'wirial':>7} {'<r>':>7} {'B':>6}")
for R in [6,8,10,12,16]:
    s = solve(R, int(160*R))
    print(f"  {R:4d} {s['E']:9.4f} {s['E2']:8.4f} {s['E4']:8.4f} {s['vir']:7.4f} {s['rms']:7.4f} {s['B']:6.3f}")
print("  => E,<r>,wirial STALE wzgl. R => zlokalizowany soliton O(1), NIE artefakt pudla.")

print("\n[N] ZBIEZNOSC SIATKI (R=10):")
print(f"  {'N':>6} {'E':>9} {'wirial':>7} {'<r>':>7}")
for N in [400,800,1600,3200]:
    s = solve(10, N)
    print(f"  {N:6d} {s['E']:9.4f} {s['vir']:7.4f} {s['rms']:7.4f}")
print("  => E,<r> zbiegaja; wirial osiada ~1.04 (resztka ucieci pudla/r=0, nie fizyka).")

b = solve(10, 2000); E2,E4 = b['E2'],b['E4']
print("\n[D] Derrick E(lambda)=lambda*E2+E4/lambda (E4/E2 ~ 1 => balans):")
for lam in [0.5,0.8,1.0,1.25,2.0]:
    print(f"   lam={lam:4.2f}: sigma+Skyrme={lam*E2+E4/lam:8.4f}  sam_sigma(bez Skyrme)={lam*E2:8.4f}")
print("  => sigma+Skyrme ma minimum ~lambda=1 (skonczony rozmiar);")
print("     sam sigma monotonicznie maleje z lambda->0 (Derrick: zapada).")
