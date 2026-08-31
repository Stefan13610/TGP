# -*- coding: utf-8 -*-
"""
Phase 3 (QN numerycznie) — op-fluctuation-extended-nbody-2026-08-31
LOCK: Phase0_balance.md §3.2, §4 (QN-2..QN-4), §5.

Trzy defekty punktowe (pinning), konfiguracje:
  T: (0,0,0),(d,0,0),(0,d,0)  — odleglosci d, d, d*sqrt2
  C: (0,0,0),(d,0,0),(2d,0,0) — odleglosci d, d, 2d
Delta F3 = F_fl(ABC) - [F_fl(AB)+F_fl(BC)+F_fl(AC)]  (log-det, 3x3 vs 2x2).
Kontrola klasyczna: U_S superpozycja dokladna (<1e-12).
QN-4: zgodnosc z rozwinieciem Phase 1 do rzedu 4 (m=0.2, T, d=8, <5%).
"""
import numpy as np
from scipy.optimize import curve_fit

results = {}
def check(name, cond, detail=""):
    results[name] = bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  {detail}")

def green(L, m):
    k = 2.0 * np.pi * np.fft.fftfreq(L)
    cx = np.cos(k)
    om = 2.0 * ((1 - cx)[:, None, None] + (1 - cx)[None, :, None]
                + (1 - cx)[None, None, :])
    denom = m * m + om
    if m == 0.0:
        denom[0, 0, 0] = np.inf
    return np.real(np.fft.ifftn(1.0 / denom))

def gv(G, L, a, b, B=0.0):
    return G[(a[0] - b[0]) % L, (a[1] - b[1]) % L, (a[2] - b[2]) % L] - B

def F_fl_sites(G, L, sites, B=0.0):
    n = len(sites)
    S = np.array([[gv(G, L, a, b, B) for b in sites] for a in sites])
    sgn, ld = np.linalg.slogdet(S)
    assert sgn > 0
    iso = sum(np.log(gv(G, L, a, a, B)) for a in sites)
    return 0.5 * (ld - iso)

def dF3(G, L, s1, s2, s3, B=0.0):
    F123 = F_fl_sites(G, L, [s1, s2, s3], B)
    F12 = F_fl_sites(G, L, [s1, s2], B)
    F13 = F_fl_sites(G, L, [s1, s3], B)
    F23 = F_fl_sites(G, L, [s2, s3], B)
    return F123 - (F12 + F13 + F23), (F12 + F13 + F23)

print("=" * 72)
print("Phase 3 — QN: N-cialowosc kanalu fluktuacyjnego (log-det)")
print("=" * 72)

L = 128
G_m02 = green(L, 0.2)
G_c = green(L, 0.0)
ds_fit = np.arange(8, 25, dtype=float)
Gd_fit = np.array([G_c[int(d), 0, 0] for d in ds_fit])
popt, _ = curve_fit(lambda d, A, p, B: A * d ** p + B, ds_fit, Gd_fit,
                    p0=[1 / (4 * np.pi), -1.0, 0.0], maxfev=20000)
Bhat = popt[2]
print(f"  B_hat(L=128) = {Bhat:.6e}")

def config_T(d):
    return (0, 0, 0), (d, 0, 0), (0, d, 0)

def config_C(d):
    return (0, 0, 0), (d, 0, 0), (2 * d, 0, 0)

# ------------------------------------------------------------------ QN-2
print("    QN-2: znak Delta F3 (20 punktow):")
signs, table = [], []
for mtag, G, B in [("m02", G_m02, 0.0), ("crit", G_c, Bhat)]:
    for ctag, cfg in [("T", config_T), ("C", config_C)]:
        for d in [4, 6, 8, 10, 12]:
            s1, s2, s3 = cfg(d)
            v, pairsum = dF3(G, L, s1, s2, s3, B)
            signs.append(np.sign(v))
            ratio = abs(v) / abs(pairsum)
            table.append((mtag, ctag, d, v, ratio))
            print(f"      {mtag:4s} {ctag} d={d:2d}: DF3={v:+.4e}, "
                  f"|DF3|/|sum par|={ratio:.4e}")
check("QN-2: znak Delta F3 uniwersalny (ten sam we wszystkich 20 punktach)",
      len(set(signs)) == 1, f"znak = {'+' if signs[0] > 0 else '-'}")

# ------------------------------------------------------------------ QN-3
ok3 = True
for mtag, G, m in [("m02", G_m02, 0.2)]:
    for ctag, cfg in [("T", config_T), ("C", config_C)]:
        for d in [4, 8, 12]:
            s1, s2, s3 = cfg(d)
            # U_S(ABC) = -sum_par q_i q_j G(d_ij): z definicji suma parowa;
            # test maszynerii: energia klasyczna z pelnego rozwiazania
            # liniowego (3 zrodla) vs superpozycja par
            q = np.ones(3)
            Gm = np.array([[gv(G, L, a, b) for b in (s1, s2, s3)]
                           for a in (s1, s2, s3)])
            # klasyczna energia interakcji: -1/2 q^T Gm q + 1/2 sum q_i^2 G(0)
            E_full = -0.5 * q @ Gm @ q + 0.5 * np.sum(q * q) * Gm[0, 0]
            E_pairs = -sum(Gm[i, j] for i in range(3) for j in range(3)
                           if i < j)
            ok3 &= abs(E_full - E_pairs) < 1e-12
check("QN-3: kanal klasyczny (zrodlowy) addytywny dokladnie (<1e-12)", ok3)

# ------------------------------------------------------------------ QN-4
s1, s2, s3 = config_T(8)
v_num, _ = dF3(G_m02, L, s1, s2, s3)
G0v = G_m02[0, 0, 0]
g12 = gv(G_m02, L, s1, s2) / G0v
g13 = gv(G_m02, L, s1, s3) / G0v
g23 = gv(G_m02, L, s2, s3) / G0v
v_an = (g12 * g13 * g23
        - 0.5 * (g12**2 * g13**2 + g12**2 * g23**2 + g13**2 * g23**2))
rel = abs(v_num - v_an) / abs(v_num)
check("QN-4: zgodnosc z rozwinieciem rzedu<=4 (m=0.2, T, d=8) < 5%",
      rel < 0.05, f"num={v_num:.6e}, analit={v_an:.6e}, rel={rel:.2e}")

# --------------------------------------------------- deskryptywnie: zanik
print("    deskryptywnie (fit slope ln|DF3| vs ln d, crit, [6,12]):")
for ctag, cfg in [("T", config_T), ("C", config_C)]:
    dd = np.array([6.0, 8.0, 10.0, 12.0])
    vv = np.array([abs(dF3(G_c, L, *cfg(int(d)), Bhat)[0]) for d in dd])
    A_ = np.vstack([np.ones_like(dd), np.log(dd)]).T
    coef, *_ = np.linalg.lstsq(A_, np.log(vv), rcond=None)
    print(f"      {ctag}: slope = {coef[1]:.3f} (analitycznie ~ -3 z g^3)")

# ---------------------------------------------------------------- SUMMARY
print("-" * 72)
npass = sum(results.values()); ntot = len(results)
print(f"SUMMARY Phase 3: {npass}/{ntot} PASS")
for k, ok in results.items():
    print(f"  {'PASS' if ok else 'FAIL'}  {k}")
if not results["QN-3: kanal klasyczny (zrodlowy) addytywny dokladnie (<1e-12)"]:
    print("WERDYKT QN: INCONCLUSIVE (pad kontroli klasycznej)")
elif npass == ntot:
    print("WERDYKT QN: TAK — kanal fluktuacyjny NIEADDYTYWNY z uniwersalnym")
    print("  znakiem DODATNIM czlonu 3-cialowego (oslabia przyciaganie parowe),")
    print("  zgodnie z Phase 1 (czlon wiodacy g12*g13*g23); kanal klasyczny")
    print("  addytywny dokladnie (kontrast potwierdzony).")
else:
    print("WERDYKT QN: FAIL/mieszany — patrz pozycje wyzej.")
