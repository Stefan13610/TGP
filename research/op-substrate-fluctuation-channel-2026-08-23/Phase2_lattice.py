# -*- coding: utf-8 -*-
"""
Phase 2 (QF numerycznie) — op-substrate-fluctuation-channel-2026-08-23
LOCK: Phase0_balance.md §3.1, §4 (QF-1..QF-4), §5 (okna) + Amendment A1.

Siatka kubiczna L^3 periodyczna, funkcja Greena przez FFT:
  G_m(r) = (1/L^3) sum_k e^{ikr} / (m^2 + sum_a 2(1-cos k_a));
  m=0: mod zerowy k=0 usuniety (projekcja, LOCK §3.1).

T0 (maszyneria): (-Delta_lat + m^2) G = delta_0 (rezidua < 1e-10).
QF-1: znaki + monotonicznosc F_fl w oknach; tabela znakow kanalow.
QF-2: zasiegi kappa_S = mu(m), kappa_D = 2 mu(m) +/- 10%.
QF-3 (A1): m=0 fit G = A d^p + B, p = -1 +/- 0.1; F_fl (connected) slope -2 +/- 0.2.
QF-4: kontrole negatywne (m=2 zanik; obrazy periodyczne; G_m > 0).
"""
import numpy as np

try:
    from scipy.optimize import curve_fit
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

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
        denom[0, 0, 0] = np.inf   # projekcja modu zerowego
    return np.real(np.fft.ifftn(1.0 / denom))

def mu_lat(m):
    return np.arccosh(1.0 + m * m / 2.0)

def linfit(x, y):
    A = np.vstack([np.ones_like(x), x]).T
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return coef[0], coef[1], r2

def F_fl(G0, Gd):
    return 0.5 * np.log(1.0 - (Gd / G0) ** 2)

print("=" * 72)
print("Phase 2 — QF numerycznie (siatka, FFT)")
print("=" * 72)

cases = {}
for L, ms in [(64, [0.0, 0.2, 0.5, 2.0]), (128, [0.0, 0.05, 0.2])]:
    for m in ms:
        G = green(L, m)
        cases[(L, m)] = G
        print(f"  policzono G: L={L}, m={m}, G(0)={G[0,0,0]:.6f}")

# ------------------------------------------------------------------- T0
ok_T0 = True
for (L, m), G in cases.items():
    if m == 0.0:
        continue  # projekcja psuje delte o stala 1/L^3 — test dla masywnych
    def lapG(i):
        s = (G[(i[0]+1) % L, i[1], i[2]] + G[(i[0]-1) % L, i[1], i[2]]
             + G[i[0], (i[1]+1) % L, i[2]] + G[i[0], (i[1]-1) % L, i[2]]
             + G[i[0], i[1], (i[2]+1) % L] + G[i[0], i[1], (i[2]-1) % L]
             - 6 * G[i[0], i[1], i[2]])
        return s
    r0 = -(lapG((0, 0, 0))) + m * m * G[0, 0, 0] - 1.0
    r5 = -(lapG((5, 0, 0))) + m * m * G[5, 0, 0]
    ok_T0 &= (abs(r0) < 1e-10) and (abs(r5) < 1e-10)
check("T0: (-Delta+m^2)G = delta_0 (masywne, rezidua <1e-10)", ok_T0)

# ------------------------------------------------------------------ QF-1
windows_QF1 = [(128, 0.05, 10, 30), (64, 0.2, 5, 15), (64, 0.5, 2, 8),
               (64, 0.0, 6, 16), (128, 0.0, 8, 24)]
ok_sign, ok_mono = True, True
for (L, m, d1, d2) in windows_QF1:
    G = cases[(L, m)]
    ds = np.arange(d1, d2 + 1)
    Fv = np.array([F_fl(G[0, 0, 0], G[d, 0, 0]) for d in ds])
    s_ok = np.all(Fv < 0)
    m_ok = np.all(np.diff(Fv) > 0)
    ok_sign &= s_ok; ok_mono &= m_ok
    print(f"    QF-1 L={L} m={m} okno[{d1},{d2}]: F_fl<0: {s_ok}, "
          f"monotoniczne: {m_ok}, F_fl({d1})={Fv[0]:.3e}, F_fl({d2})={Fv[-1]:.3e}")
check("QF-1a: F_fl < 0 we wszystkich oknach", ok_sign)
check("QF-1b: F_fl rosnace z d (przyciaganie) we wszystkich oknach", ok_mono)

# tabela znakow: m=0.2, L=64, d=8
G = cases[(64, 0.2)]; G0v, Gdv = G[0, 0, 0], G[8, 0, 0]
def F_class(v1, v2, G0v, Gdv):
    Sig = np.array([[G0v, Gdv], [Gdv, G0v]])
    v = np.array([v1, v2])
    return 0.5 * v @ np.linalg.solve(Sig, v) - (v1**2 + v2**2) / (2 * G0v)
a = 0.3
tab = {
    "M-D klas. (+a,+a)": F_class(a, a, G0v, Gdv),
    "M-D klas. (+a,-a)": F_class(a, -a, G0v, Gdv),
    "M-D klas. (+a, 0)": F_class(a, 0.0, G0v, Gdv),
    "M-S (q=+1,+1)": -1.0 * 1.0 * Gdv,
    "M-S (q=+1,-1)": -1.0 * (-1.0) * Gdv,
    "M-D flukt. (kazde v)": F_fl(G0v, Gdv),
}
print("    tabela znakow (m=0.2, L=64, d=8):")
for k2, val in tab.items():
    print(f"      {k2:24s} = {val:+.4e}  ({'przyciaganie' if val < 0 else 'odpychanie'})")
flip_ok = (tab["M-D klas. (+a,+a)"] < 0 and tab["M-D klas. (+a,-a)"] > 0
           and tab["M-S (q=+1,+1)"] < 0 and tab["M-S (q=+1,-1)"] > 0
           and tab["M-D flukt. (kazde v)"] < 0)
check("QF-1c: kanaly klasyczne flipuja znak, fluktuacyjny NIE", flip_ok)

# ------------------------------------------------------------------ QF-2
ok_QF2 = True
for (L, m, d1, d2) in [(64, 0.5, 2, 8), (64, 0.2, 5, 15)]:
    G = cases[(L, m)]
    ds = np.arange(d1, d2 + 1, dtype=float)
    Gd_ = np.array([G[int(d), 0, 0] for d in ds])
    Fd_ = np.array([abs(F_fl(G[0, 0, 0], G[int(d), 0, 0])) for d in ds])
    # ln G + ln d = a - kappa_S d ; ln|F| + 2 ln d = a - kappa_D d
    _, s_S, r2S = linfit(ds, np.log(Gd_) + np.log(ds))
    _, s_D, r2D = linfit(ds, np.log(Fd_) + 2 * np.log(ds))
    kS, kD = -s_S, -s_D
    mu = mu_lat(m)
    okS = abs(kS - mu) / mu < 0.10 and r2S > 0.99
    okD = abs(kD - 2 * mu) / (2 * mu) < 0.10 and r2D > 0.99
    ok_QF2 &= okS and okD
    print(f"    QF-2 m={m}: kappa_S={kS:.4f} vs mu={mu:.4f} "
          f"({100*abs(kS-mu)/mu:.1f}%, R2={r2S:.5f}) | "
          f"kappa_D={kD:.4f} vs 2mu={2*mu:.4f} "
          f"({100*abs(kD-2*mu)/(2*mu):.1f}%, R2={r2D:.5f})")
check("QF-2: zasiegi kappa_S=mu, kappa_D=2mu (+/-10%, R2>0.99)", ok_QF2)

# ------------------------------------------------------------- QF-3 (A1)
ok_QF3 = True
p_fit = {}
for (L, d1, d2) in [(64, 6, 16), (128, 8, 24)]:
    G = cases[(L, 0.0)]
    ds = np.arange(d1, d2 + 1, dtype=float)
    Gd_ = np.array([G[int(d), 0, 0] for d in ds])
    if HAVE_SCIPY:
        def model(d, A, p, B):
            return A * d ** p + B
        popt, _ = curve_fit(model, ds, Gd_, p0=[1 / (4 * np.pi), -1.0, 0.0],
                            maxfev=20000)
        A_, p_, B_ = popt
        yhat = model(ds, *popt)
        r2 = 1 - np.sum((Gd_ - yhat) ** 2) / np.sum((Gd_ - Gd_.mean()) ** 2)
    else:
        raise SystemExit("scipy wymagane dla QF-3 (A1)")
    p_fit[L] = p_
    ok_p = abs(p_ - (-1.0)) < 0.10 and r2 > 0.99
    # connected F_fl:
    G0c = G[0, 0, 0] - B_
    Fc = np.array([abs(0.5 * np.log(1 - ((G[int(d), 0, 0] - B_) / G0c) ** 2))
                   for d in ds])
    _, slope, r2F = linfit(np.log(ds), np.log(Fc))
    ok_F = abs(slope - (-2.0)) < 0.20 and r2F > 0.99
    ok_QF3 &= ok_p and ok_F
    print(f"    QF-3 L={L}: p={p_:.4f} (cel -1.0+/-0.1, R2={r2:.6f}), "
          f"B={B_:.5f}, A={A_:.5f} (1/4pi={1/(4*np.pi):.5f}); "
          f"slope F_fl={slope:.4f} (cel -2.0+/-0.2, R2={r2F:.5f})")
drift = abs(p_fit[64] - p_fit[128])
ok_QF3 &= drift < 0.10
print(f"    QF-3 dryf p(64)-p(128) = {drift:.4f} (cel <0.1)")
check("QF-3: krytyczne p=-1, slope F_fl=-2, dryf<0.1 (Amendment A1)", ok_QF3)

# ------------------------------------------------------------------ QF-4
G2 = cases[(64, 2.0)]
ratio = abs(F_fl(G2[0, 0, 0], G2[6, 0, 0])) / abs(F_fl(G2[0, 0, 0], G2[2, 0, 0]))
check("QF-4a: m=2 zanik |F_fl(6)|/|F_fl(2)| < 1e-3", ratio < 1e-3,
      f"ratio = {ratio:.3e}")
G64, G128 = cases[(64, 0.2)], cases[(128, 0.2)]
rels = [abs(G64[d, 0, 0] - G128[d, 0, 0]) / abs(G128[d, 0, 0])
        for d in range(1, 17)]
check("QF-4b: obrazy periodyczne — |G64-G128|/G128 < 1e-2 (d<=16, m=0.2)",
      max(rels) < 1e-2, f"max rel = {max(rels):.3e}")
ok_pos = True
for (L, m), G in cases.items():
    if m > 0:
        ok_pos &= (G.min() > 0)
check("QF-4c: G_m > 0 wszedzie (przypadki masywne)", ok_pos)

# ---------------------------------------------------------------- SUMMARY
print("-" * 72)
npass = sum(results.values()); ntot = len(results)
print(f"SUMMARY Phase 2: {npass}/{ntot} PASS")
for k2, ok in results.items():
    print(f"  {'PASS' if ok else 'FAIL'}  {k2}")
qf4 = results["QF-4a: m=2 zanik |F_fl(6)|/|F_fl(2)| < 1e-3"] and \
      results["QF-4b: obrazy periodyczne — |G64-G128|/G128 < 1e-2 (d<=16, m=0.2)"] and \
      results["QF-4c: G_m > 0 wszedzie (przypadki masywne)"]
qf123 = all(results[k2] for k2 in results if k2.startswith(("QF-1", "QF-2", "QF-3")))
if not qf4:
    print("WERDYKT QF: INCONCLUSIVE (pad kontroli QF-4 — maszyneria)")
elif qf123 and results["T0: (-Delta+m^2)G = delta_0 (masywne, rezidua <1e-10)"]:
    print("WERDYKT QF: PASS (QF-1 ∧ QF-2 ∧ QF-3 ∧ QF-4)")
else:
    print("WERDYKT QF: FAIL (co najmniej jedno z QF-1..QF-3 przy zdrowym QF-4)")
