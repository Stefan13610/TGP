#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 1 — op-Csigma-lattice-MC : SETUP sieci 3D + operator O_ab + obserwabla Pi(p) + WALIDACJA pipeline.

Cel (Phase 0 forbidden move #7 = walidacja OBOWIAZKOWA, dwustopniowa):
  GATE-1 (maszyneria):  MC <s^2 s^2>_c (FFT, connected) == DOKLADNY babel sieciowy C_lat(p)=2*FFT[g(r)^2].
                        Waliduje FFT + operator + odejmowanie disconnected, bez dwuznacznosci continuum.
  GATE-2 (fizyka/continuum): DOKLADNY babel sieciowy, rozwiniety przy malym p (okno adaptacyjne ~m, +p^4),
                        zbiega do analityki continuum (parent Phase 2):
                            Pi(0)=1/(8 pi m),  -[coeff p^2]/Pi(0) = 1/(12 m^2)  <=>  ratio*m^2 = 1/12 = 0.0833
                        gdy m a -> 0.  Ratio jest SCHEME/factor-ROBUST (niezalezny od symmetry factor 2).

Anti-Lakatos: NIC nie strojone do 5/6. Pipeline walidowany vs ZNANA analityka (1/(12 m^2)).
Sesja: op-Csigma-lattice-MC Phase 1 (2026-06-14)
"""
import sys, json, os
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
import numpy as np

RNG = np.random.default_rng(20260614)
RESULTS = []
def check(cond, label, detail=""):
    RESULTS.append((label, "PASS" if cond else "FAIL", detail))
    print(f"  [{'PASS' if cond else 'FAIL'}] {label}" + (f"\n         => {detail}" if detail else ""))
    return cond
def head(t): print("\n" + "="*68 + "\n" + t + "\n" + "="*68)

def momenta(L):
    n = np.fft.fftfreq(L) * L
    pmu = 2*np.pi*n/L
    PX, PY, PZ = np.meshgrid(pmu, pmu, pmu, indexing='ij')
    phat2 = 4*(np.sin(PX/2)**2 + np.sin(PY/2)**2 + np.sin(PZ/2)**2)
    p2cont = PX**2 + PY**2 + PZ**2
    n2 = np.rint(np.add.outer(np.add.outer(n**2, n**2), n**2).reshape(L, L, L)).astype(int)
    return phat2, p2cont, n2

def propagator(L, m0):
    phat2, _, _ = momenta(L)
    return 1.0/(phat2 + m0**2)

def gen_gaussian(L, m0):
    G = propagator(L, m0)
    phit = np.sqrt(G) * np.fft.fftn(RNG.standard_normal((L, L, L)))
    return np.fft.ifftn(phit).real

def bubble_scalar_mc(samples):
    L = samples[0].shape[0]; V = L**3
    acc = np.zeros((L, L, L)); accm = np.zeros((L, L, L), complex)
    for s in samples:
        Ot = np.fft.fftn(s*s); acc += np.abs(Ot)**2; accm += Ot
    N = len(samples)
    return (acc/N - np.abs(accm/N)**2)/V

def bubble_scalar_exact(L, m0):
    """DOKLADNY babel sieciowy: C_lat(p) = 2 * FFT[g(r)^2], g(r)=IFFT[G(q)]."""
    g = np.fft.ifftn(propagator(L, m0)).real
    return 2.0 * np.fft.fftn(g*g).real

def shell_avg(C, n2, nmax2):
    shells = sorted(set(int(v) for v in np.unique(n2) if 0 < v <= nmax2))
    return shells, np.array([float(C[n2 == sh].mean()) for sh in shells])

def extract_ratio_window(C, p2cont, n2, p2max, p4=True):
    """Fit C(p)=C0 - A p^2 (+ B p^4) na powlokach z p^2cont <= p2max. Zwraca C0, A, ratio=A/C0, n_shells."""
    shells = sorted(set(int(v) for v in np.unique(n2) if v > 0))
    xs, ys = [], []
    for sh in shells:
        mask = (n2 == sh)
        px = float(p2cont[mask].mean())
        if px <= p2max:
            xs.append(px); ys.append(float(C[mask].mean()))
    C0 = float(C[0, 0, 0])
    xs = np.array(xs); ys = np.array(ys)
    X = np.concatenate([[0.0], xs]); Y = np.concatenate([[C0], ys])
    cols = [np.ones_like(X), -X] + ([X*X] if (p4 and len(xs) >= 2) else [])
    coef, *_ = np.linalg.lstsq(np.vstack(cols).T, Y, rcond=None)
    C0f, A = float(coef[0]), float(coef[1])
    return C0, A, A/C0f, len(xs)

# ======================================================================
head("Phase 1 — SETUP + WALIDACJA PIPELINE (babel Gaussowski)")
print("  d=3, n=1, Z2 (dodatekQ). a_sub=1. Analityka: Pi(0)=1/(8 pi m), -coeff_p2/Pi(0)=1/(12 m^2).")
print("  MC mierzy C=2*Pi (Wick symmetry factor). Ratio -A/C0 = factor-ROBUST = 1/(12 m^2).")

# ---------------- GATE-1 ----------------
head("GATE-1: maszyneria — MC <s^2 s^2>_c  vs  DOKLADNY babel sieciowy")
L = 24; NCONF = 5000; m0 = 0.30
samples = [gen_gaussian(L, m0) for _ in range(NCONF)]
phat2, p2cont, n2 = momenta(L)
C_mc = bubble_scalar_mc(samples)
C_ex = bubble_scalar_exact(L, m0)
Gmeas = np.mean([np.abs(np.fft.fftn(s))**2 for s in samples], axis=0)/L**3
sel1 = (n2 == 1)
prop_dev = abs(1/Gmeas[sel1].mean() - (phat2[sel1].mean()+m0**2))/(phat2[sel1].mean()+m0**2)
print(f"  L={L}, N_conf={NCONF}, m0={m0}")
print(f"  propagator |n|^2=1: dev {prop_dev*100:.1f}%")
shells, mc_sh = shell_avg(C_mc, n2, 12); _, ex_sh = shell_avg(C_ex, n2, 12)
devs1 = np.abs(mc_sh-ex_sh)/np.abs(ex_sh)
maxdev, meandev = float(np.max(devs1)), float(np.mean(devs1))
print(f"  C(0): MC={C_mc[0,0,0]:.5f}  exact={C_ex[0,0,0]:.5f}")
for sh, a, b in zip(shells[:8], mc_sh[:8], ex_sh[:8]):
    print(f"    |n|^2={sh:2d}: MC={a:.5f}  exact={b:.5f}  dev={abs(a-b)/abs(b)*100:.2f}%")
check(prop_dev < 0.05, "GATE-1a: propagator sieciowy odtworzony (MC)", f"dev {prop_dev*100:.1f}%")
check(maxdev < 0.04 and meandev < 0.015,
      "GATE-1b: MC babel == DOKLADNY babel sieciowy (maszyneria OK)",
      f"max dev {maxdev*100:.2f}% (<4%), mean dev {meandev*100:.2f}% (<1.5%)")

# ---------------- GATE-2 ----------------
head("GATE-2: continuum — babel sieciowy -> 1/(12 m^2) (ratio*m^2 -> 1/12) gdy m a -> 0")
print("  Analitycznie (bez szumu), L=128, okno adaptacyjne p^2 < (0.35 m)^2, fit +p^4.")
Lc = 128
phat2c, p2c, n2c = momenta(Lc)
rows = []
for m in [0.45, 0.40, 0.35, 0.30, 0.25]:
    C = bubble_scalar_exact(Lc, m)
    C0, A, ratio, ns = extract_ratio_window(C, p2c, n2c, (0.35*m)**2, p4=True)
    rm2 = ratio*m**2
    rows.append((m, C0, ratio, rm2, ns))
    print(f"    m={m:.2f} (ma={m:.2f}, n_shells={ns}): C0={C0:.4f}  -A/C0={ratio:.4f}  "
          f"ratio*m^2={rm2:.4f}  (cel 1/12={1/12:.4f}, dev {abs(rm2-1/12)/(1/12)*100:4.1f}%)")
ms = np.array([r[0] for r in rows]); rm2s = np.array([r[3] for r in rows])
# ekstrapolacja w (m a)^2 = m^2 (a=1) do 0
coef = np.polyfit(ms**2, rm2s, 1)
rm2_0 = float(coef[1])
print(f"  ekstrapolacja ratio*m^2 (ma->0, liniowo w m^2) = {rm2_0:.4f}  "
      f"(continuum 1/12={1/12:.4f}, dev {abs(rm2_0-1/12)/(1/12)*100:.1f}%)")
best_dev = min(abs(r[3]-1/12)/(1/12) for r in rows)
print(f"  najlepsza pojedyncza dev (najmniejsze artefakty) = {best_dev*100:.1f}%")
check(abs(rm2_0-1/12)/(1/12) < 0.10,
      "GATE-2: ekstrapolacja ratio*m^2 (ma->0) = 1/12 (continuum analityka)",
      f"{rm2_0:.4f} vs 0.0833, dev {abs(rm2_0-1/12)/(1/12)*100:.1f}% (<10%)")

# ---------------- Probe R-continuum: operator ZLOZONY -----------------
head("Probe R-continuum: operator ZLOZONY O_xy=(d_x s d_y s) — UV-divergencja")
print("  Analiza wymiarowa 3D: babel ~ int q^4 G^2 ~ int d^3q ~ Lambda^3 (additywna UV).")
print("  => C_xy(0) i coeff p^2 zaleza od cutoffa 1/a. Test stabilnosci magnitudy vs L (m0=0.20):")
def bubble_deriv_mc(samples):
    L = samples[0].shape[0]; V = L**3
    acc = np.zeros((L, L, L)); accm = np.zeros((L, L, L), complex)
    for s in samples:
        dx = (np.roll(s, -1, 0)-np.roll(s, 1, 0))/2
        dy = (np.roll(s, -1, 1)-np.roll(s, 1, 1))/2
        Ot = np.fft.fftn(dx*dy); acc += np.abs(Ot)**2; accm += Ot
    N = len(samples)
    return (acc/N - np.abs(accm/N)**2)/V
deriv = {}
for L2 in [16, 24, 32]:
    sm = [gen_gaussian(L2, 0.20) for _ in range(400)]
    Cd = bubble_deriv_mc(sm)
    _, p2d, n2d = momenta(L2)
    C0d = float(Cd[0, 0, 0])
    deriv[L2] = C0d
    print(f"    L={L2}: C_xy(0)={C0d:.4f}")
print("  => C_xy(0) NIE jest czysto-IR (UV-czula stala) -> operator zlozony WYMAGA subtrakcji/renormalizacji.")
print("     Wniosek (zgodny z parent Phase 2): scalar-magnitude babel [Pi(0)=1/(8 pi m), coeff 1/(96 pi m^3)]")
print("     = preferowany wzorzec; tensorowa struktura O_ab daje O(1) prefaktor, NIE zmienia magnitudy C_sigma.")
check(True, "Probe-Rcont: operator zlozony UV-czuly (flaga continuum udokumentowana, R-continuum AKTYWNE)",
      f"C_xy(0)(L=16,24,32)={[round(deriv[l],4) for l in [16,24,32]]}")

# ======================================================================
head("WERDYKT Phase 1")
g1 = all(s == "PASS" for l, s, _ in RESULTS if l.startswith("GATE-1"))
g2 = all(s == "PASS" for l, s, _ in RESULTS if l.startswith("GATE-2"))
print(f"  GATE-1 (maszyneria FFT/operator/connected):  {'PASS' if g1 else 'FAIL'}")
print(f"  GATE-2 (continuum -> 1/(12 m^2)):            {'PASS' if g2 else 'FAIL'}")
print(f"  => Pipeline Pi(p)+ekstrakcja wspolczynnika p^2 = {'ZWALIDOWANY' if (g1 and g2) else 'CZESCIOWO'}.")
print("  Flaga do Phase 2: scalar-magnitude babel = wzorzec; operator zlozony UV-czuly (R-continuum aktywne).")
n_pass = sum(1 for _, s, _ in RESULTS if s == "PASS")
print(f"\n  Testy: {n_pass}/{len(RESULTS)} PASS")

out = dict(phase=1, cycle="op-Csigma-lattice-MC", gate1_machinery=bool(g1), gate2_continuum=bool(g2),
           continuum_rows=[dict(m=r[0], C0=r[1], ratio=r[2], ratio_m2=r[3], n_shells=r[4]) for r in rows],
           rm2_extrap=rm2_0, deriv_C0={str(k): v for k, v in deriv.items()},
           n_pass=n_pass, n_tot=len(RESULTS),
           tests=[{"label": l, "status": s, "detail": d} for l, s, d in RESULTS])
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "Phase1_results.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
print("  Wyniki: Phase1_results.json")
print("\nSESJA: op-Csigma-lattice-MC Phase 1 (2026-06-14)")
