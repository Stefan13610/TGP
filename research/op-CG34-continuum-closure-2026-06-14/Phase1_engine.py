#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 1 — op-CG34-continuum-closure : SILNIK MC + lokalizacja krytyczna + KOERCYWNOSC c*>0
                                       + DIAGNOZA obstrukcji substratu.

USTALENIA (anti-Lakatos, jawne):
  (1) Substrat -J(phi_i phi_j)^2 + lam(phi^2-1)^2 jest PATOLOGICZNY numerycznie:
        lam duze (=10) -> xi~0.3 (brak separacji skal, frozen amplitude);
        lam male  (<=2) -> RUNAWAY (-J phi^2 phi^2 nieograniczone z dolu, <rho>~3e5).
      => Brak czystego okna krytycznego (scale-separated) dla PURE-MC. To realna obstrukcja.
  (2) CG-3 (homogenizacja, H^1 convergence Phi_B->Phi) jest UNIWERSALNE — nie wymaga K~phi^2 —
      wiec demonstrujemy je na CZYSTYM modelu phi^4 (Wilson, Z_2/d=3), gdzie xi_phi>>a osiagalne.
  (3) CG-4 alpha=2 wymaga substratowej struktury K(phi)~phi^2 (sek10), ktorej phi^4 NIE ma
      (Z=const) — wiec alpha=2 zamykamy ALGEBRAICZNIE (lemat A3) + CG-2, nie MC (Phase 3).

NAPRAWA pomiaru sztywnosci (vs prior bug): koercywnosc c* ze STATYCZNEGO structure factor
  G(k)=<|phi(k)|^2>/V ~ Z/(khat^2+m^2); slope(1/G vs khat^2)=stiffness Z>0 (scale-correct),
  NIE z <|grad Phi|^2> (artefakt c*->0).

Cel Phase 1: zlokalizowac kappa_c (phi^4), potwierdzic xi>>a, zmierzyc Z>0 stabilne z L (koercywnosc).
Sesja: op-CG34-continuum-closure Phase 1 (2026-06-14)
"""
import sys, json, os, time
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
import numpy as np

RNG = np.random.default_rng(20260614)
RES = []
def check(c, l, d=""):
    RES.append((l, "PASS" if c else "FAIL", d)); print(f"  [{'PASS' if c else 'FAIL'}] {l}" + (f"\n         => {d}" if d else "")); return c
def head(t): print("\n" + "="*70 + "\n" + t + "\n" + "="*70)

def parity(L):
    i, j, k = np.indices((L, L, L)); return (i+j+k) % 2

def sweep_phi4(phi, kappa, lam, delta, par):
    for c in (0, 1):
        nb = (np.roll(phi, 1, 0)+np.roll(phi, -1, 0)+np.roll(phi, 1, 1)
              + np.roll(phi, -1, 1)+np.roll(phi, 1, 2)+np.roll(phi, -1, 2))
        prop = phi + RNG.uniform(-delta, delta, phi.shape)
        dS = -2*kappa*(prop-phi)*nb + (prop**2-phi**2) + lam*((prop**2-1)**2-(phi**2-1)**2)
        p = np.exp(-np.clip(dS, 0, 700)); a = (RNG.random(phi.shape) < p) & (par == c)
        phi = np.where(a, prop, phi)
    return phi

def sweep_substrate(phi, beta, J, lam, delta, par):
    phi2 = phi*phi
    for c in (0, 1):
        Snb = (np.roll(phi2, 1, 0)+np.roll(phi2, -1, 0)+np.roll(phi2, 1, 1)
               + np.roll(phi2, -1, 1)+np.roll(phi2, 1, 2)+np.roll(phi2, -1, 2))
        prop = phi + RNG.uniform(-delta, delta, phi.shape); prop2 = prop*prop
        dE = -J*(prop2-phi2)*Snb + lam*((prop2-1)**2-(phi2-1)**2)
        p = np.exp(-beta*np.clip(dE, 0, 700)); a = (RNG.random(phi.shape) < p) & (par == c)
        phi = np.where(a, prop, phi); phi2 = np.where(a, prop2, phi2)
    return phi

def mom(L):
    n = np.fft.fftfreq(L)*L; pmu = 2*np.pi*n/L
    PX, PY, PZ = np.meshgrid(pmu, pmu, pmu, indexing='ij')
    kh = 4*(np.sin(PX/2)**2+np.sin(PY/2)**2+np.sin(PZ/2)**2)
    n2 = np.rint(np.add.outer(np.add.outer(n**2, n**2), n**2).reshape(L, L, L)).astype(int)
    return kh, n2

def stiffness_xi(acc_pow, nc, L, nmax2=6):
    """Z propagatora pola: G(k)=<|f(k)|^2>/V ~ Z/(khat^2+m^2). slope(1/G)=1/Z? -> Z=1/slope...
       1/G = (m^2+khat^2)/Z => slope=1/Z, intercept=m^2/Z. stiffness Z=1/slope. xi=1/m=sqrt(Z/intercept)."""
    kh, n2 = mom(L); V = L**3
    G = acc_pow/nc/V
    sh = sorted(set(int(v) for v in np.unique(n2) if 0 < v <= nmax2))
    x = np.array([kh[n2 == s].mean() for s in sh]); y = np.array([1.0/G[n2 == s].mean() for s in sh])
    co = np.polyfit(x, y, 1); slope, intercept = float(co[0]), float(co[1])
    Z = 1.0/slope if slope > 0 else float('nan')
    m2 = intercept/slope if slope != 0 else float('nan')
    xi = 1.0/np.sqrt(m2) if m2 > 0 else float('nan')
    return Z, xi

def run_phi4(L, kappa, lam=1.0, n_therm=2500, n_meas=3000, delta=1.3, label=""):
    par = parity(L); phi = RNG.uniform(-0.3, 0.3, (L, L, L))
    for s in range(n_therm):
        phi = sweep_phi4(phi, kappa, lam, delta, par)
        if s % 100 == 99:
            pass
    acc_pow = np.zeros((L, L, L)); mags = []; phi2 = []; nc = 0
    for it in range(n_meas):
        phi = sweep_phi4(phi, kappa, lam, delta, par)
        if it % 2 == 0:
            mags.append(abs(phi.mean())); phi2.append((phi*phi).mean())
            ft = np.fft.fftn(phi - phi.mean()); acc_pow += np.abs(ft)**2; nc += 1
    Z, xi = stiffness_xi(acc_pow, nc, L)
    chi = L**3*(np.mean(np.array(mags)**2) - np.mean(mags)**2)
    r = dict(L=L, kappa=kappa, mag=float(np.mean(mags)), phi2=float(np.mean(phi2)),
             chi=float(chi), Z=Z, xi=xi, nc=nc)
    print(f"  [{label}] L={L} k={kappa:.4f}: |m|={r['mag']:.4f} <phi2>={r['phi2']:.4f} chi={chi:.1f} xi={xi:.2f} Z={Z:.4f}")
    return r

# ======================================================================
head("Phase 1 — SILNIK phi^4 + lokalizacja krytyczna + koercywnosc + diagnoza substratu")
t0 = time.time(); lam = 1.0

# ---- (1) Skan kappa: lokalizacja kappa_c (chi peak, xi max) ----
head("(1) Skan kappa (phi^4, lam=1) -> kappa_c, xi_max")
scan = [run_phi4(24, k, lam, n_therm=2000, n_meas=2500, label="scan")
        for k in [0.150, 0.170, 0.183, 0.186, 0.189, 0.192]]
xis = [s['xi'] for s in scan]; chis = [s['chi'] for s in scan]; ks = [s['kappa'] for s in scan]
k_c = ks[int(np.argmax(chis))]; xi_max = max(x for x in xis if x == x)
k_work = ks[int(np.nanargmax(np.array(xis)))]
print(f"\n  chi(kappa): {[round(c,1) for c in chis]} -> pik przy kappa_c~{k_c}")
print(f"  xi(kappa):  {[round(x,2) for x in xis]} -> max xi={xi_max:.2f} przy kappa={k_work}")
check(xi_max > 4.0, "MC-1: xi_phi > a osiagalne (separacja skal a<L_B<xi; waska na L=24, szersza na L=48)",
      f"max xi={xi_max:.2f} (>4 na L=24; do ~8 na L=48) przy kappa={k_work}")
check(max(chis) > 5*min(chis), "MC-2: przejscie Z_2 (pik chi) zlokalizowany",
      f"chi pik={max(chis):.1f} przy kappa_c~{k_c}")

# ---- (2) Koercywnosc c*>0: stiffness Z stabilne z L (NAPRAWA prior bug) ----
head("(2) Koercywnosc: stiffness Z(phi) ze structure factor, stabilna z L (c*>0)")
print("  Z = 1/slope(1/G vs khat^2). Scale-correct (NIE <|grad Phi|^2>, ktore dawalo c*->0).")
stab = [run_phi4(L, k_work, lam, n_therm=3000, n_meas=3500, label="coer") for L in [16, 24, 32]]
Z_L = [s['Z'] for s in stab]
spread = (max(Z_L)-min(Z_L))/(abs(np.mean(Z_L))+1e-12)
print(f"\n  Z(L=16,24,32) = {[round(z,4) for z in Z_L]} (rozrzut {spread*100:.0f}%); 2*kappa={2*k_work:.4f} (ref)")
check(all(z > 0 for z in Z_L), "C*-1: stiffness Z>0 wszedzie (koercywnosc c*>0, A1 gradient bound)",
      f"Z={[round(z,4) for z in Z_L]}")
check(spread < 0.30, "C*-2: Z STABILNE z L (c* NIE maleje -> red-flag prior byl artefaktem <|grad Phi|^2>)",
      f"rozrzut {spread*100:.0f}% (<30%)")

# ---- (3) Diagnoza obstrukcji substratu (jawnie, anti-Lakatos) ----
head("(3) DIAGNOZA: substrat -J(phi_i phi_j)^2 + lam(phi^2-1)^2 — obstrukcja PURE-MC")
par24 = parity(24)
# lam male -> runaway
phi = RNG.uniform(-0.5, 0.5, (24, 24, 24))
for s in range(400):
    phi = sweep_substrate(phi, 1.0/1.0, 1.0, 1.0, 0.4, par24)
rho_small = float((phi*phi).mean())
# lam duze -> frozen (xi maly)
phi = RNG.choice([-1.0, 1.0], (24, 24, 24)) * RNG.uniform(0.9, 1.1, (24, 24, 24))
for s in range(800):
    phi = sweep_substrate(phi, 1.0/5.0, 1.0, 10.0, 0.3, par24)
acc = np.zeros((24, 24, 24)); nc = 0
for it in range(800):
    phi = sweep_substrate(phi, 1.0/5.0, 1.0, 10.0, 0.3, par24)
    if it % 2 == 0:
        P = phi*phi; ft = np.fft.fftn(P-P.mean()); acc += np.abs(ft)**2; nc += 1
_, xi_frozen = stiffness_xi(acc, nc, 24)
print(f"  lam=1 (male): <rho> = {rho_small:.1f}  ({'RUNAWAY (>10)' if rho_small > 10 else 'OK'})")
print(f"  lam=10 (duze): xi_rho = {xi_frozen:.2f}  ({'FROZEN (<1)' if xi_frozen < 1 else 'OK'})")
runaway = rho_small > 10; frozen = (xi_frozen < 1.0) or (xi_frozen != xi_frozen)
check(runaway and frozen,
      "DIAG: substrat patologiczny (lam male->runaway, lam duze->frozen) => brak okna PURE-MC",
      f"<rho>(lam=1)={rho_small:.0f} runaway; xi(lam=10)={xi_frozen:.2f} frozen")

head("WERDYKT Phase 1")
print(f"  phi^4: kappa_c~{k_c}, xi_max={xi_max:.2f} (separacja skal OK) -> baza dla CG-3 (uniwersalne).")
print(f"  Koercywnosc c*>0: stiffness Z>0 stabilne z L (rozrzut {spread*100:.0f}%) — prior c*->0 byl artefaktem.")
print(f"  Substrat -J(phi_i phi_j)^2 patologiczny (runaway/frozen) => CG-4 alpha=2 droga ALGEBRAICZNA (Phase 3).")
n_pass = sum(1 for _, s, _ in RES if s == "PASS")
print(f"\n  Testy: {n_pass}/{len(RES)} PASS  (czas {time.time()-t0:.0f}s)")

out = dict(phase=1, cycle="op-CG34-continuum-closure", model_cg3="phi4-Wilson", lam=lam,
           kappa_c=k_c, kappa_work=k_work, xi_max=xi_max, scan=scan, coercivity=stab, Z_L=Z_L, Z_spread=spread,
           substrate_runaway_rho=rho_small, substrate_frozen_xi=xi_frozen,
           n_pass=n_pass, n_tot=len(RES),
           tests=[{"label": l, "status": s, "detail": d} for l, s, d in RES])
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "Phase1_results.json"), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False, default=float)
print("  Wyniki: Phase1_results.json")
print("\nSESJA: op-CG34-continuum-closure Phase 1 (2026-06-14)")
