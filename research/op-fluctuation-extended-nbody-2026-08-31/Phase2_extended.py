# -*- coding: utf-8 -*-
"""
Phase 2 (QE numerycznie) — op-fluctuation-extended-nbody-2026-08-31
LOCK: Phase0_balance.md §3.1, §4 (QE-0..QE-5, kontrole), §5 (okna).

Klastry = kule dyskretne {x: |x|<=R}, R in {1,2,3} (7/33/123 wezlow),
pinning; F_fl(d;R) = 1/2 [ln det Sigma(d) - ln det Sigma_A - ln det Sigma_B].
Krytycznie (m=0, projekcja modu zerowego): wpisy z propagatora connected
G_c(r) = G(r) - B_hat, B_hat z fitu S-typu na defekcie punktowym
(okno [8,24], dziedziczone z Amendment A1 rodzica; jedno B_hat na L).
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

def mu_lat(m):
    return np.arccosh(1.0 + m * m / 2.0)

def linfit(x, y):
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return coef[0], coef[1], r2

def ball_sites(R):
    if R == 0:
        return [(0, 0, 0)]
    rng = range(-R, R + 1)
    return [(x, y, z) for x in rng for y in rng for z in rng
            if x * x + y * y + z * z <= R * R]

def gval(G, L, dx, dy, dz, B=0.0):
    return G[dx % L, dy % L, dz % L] - B

def cov_matrix(G, L, sites, B=0.0):
    n = len(sites)
    S = np.empty((n, n))
    for i, a in enumerate(sites):
        for j, b in enumerate(sites):
            S[i, j] = gval(G, L, a[0] - b[0], a[1] - b[1], a[2] - b[2], B)
    return S

def F_fl_clusters(G, L, R, d, B=0.0):
    """F_fl dwoch klastrow promienia R, srodki w odleglosci d (os x)."""
    sA = ball_sites(R)
    sB = [(x + d, y, z) for (x, y, z) in sA]
    full = cov_matrix(G, L, sA + sB, B)
    SA = cov_matrix(G, L, sA, B)
    n = len(sA)
    sgn_f, ld_f = np.linalg.slogdet(full)
    sgn_a, ld_a = np.linalg.slogdet(SA)
    SB = full[n:, n:]
    sgn_b, ld_b = np.linalg.slogdet(SB)
    assert sgn_f > 0 and sgn_a > 0 and sgn_b > 0, "det <= 0 (maszyneria)"
    return 0.5 * (ld_f - ld_a - ld_b)

def local_slopes(ds, Fabs):
    """3-punktowe lokalne nachylenia log-log (centralne)."""
    ln_d, ln_F = np.log(ds), np.log(Fabs)
    return ds[1:-1], (ln_F[2:] - ln_F[:-2]) / (ln_d[2:] - ln_d[:-2])

print("=" * 72)
print("Phase 2 — QE: obiekty rozciagle (siatka, FFT, log-det)")
print("=" * 72)

G128_m02 = green(128, 0.2)
G128_c = green(128, 0.0)
G96_c = green(96, 0.0)
print(f"  G policzone: L=128 m=0.2, L=128 m=0, L=96 m=0")

# B_hat (Amendment A1 rodzica): fit G = A d^p + B, defekt punktowy, [8,24]
Bhat = {}
for L, G in [(128, G128_c), (96, G96_c)]:
    ds = np.arange(8, 25, dtype=float)
    Gd = np.array([G[int(d), 0, 0] for d in ds])
    popt, _ = curve_fit(lambda d, A, p, B: A * d ** p + B, ds, Gd,
                        p0=[1 / (4 * np.pi), -1.0, 0.0], maxfev=20000)
    A_, p_, B_ = popt
    Bhat[L] = B_
    print(f"  B_hat(L={L}) = {B_:.6e} (fit: A={A_:.5f}, p={p_:.4f})")

# ------------------------------------------------------------------ QE-0
ok0 = True
G0v = G128_m02[0, 0, 0]
for d in [4, 8, 12, 16]:
    Gd = G128_m02[d, 0, 0]
    closed = 0.5 * np.log(1.0 - (Gd / G0v) ** 2)
    viadet = F_fl_clusters(G128_m02, 128, 0, d)
    ok0 &= abs(closed - viadet) < 1e-12
    print(f"    QE-0 d={d}: zamkniete={closed:.15e}, log-det={viadet:.15e}, "
          f"|diff|={abs(closed-viadet):.2e}")
check("QE-0: R=0 log-det == forma zamknieta (<1e-12)", ok0)

# ------------------------ obliczenie F_fl dla wszystkich potrzebnych par
F = {}   # (tag, R) -> (ds, F_fl array)
cases = [("m02", G128_m02, 128, 0.0, 1, np.arange(4, 29)),
         ("m02", G128_m02, 128, 0.0, 2, np.arange(6, 29)),
         ("crit128", G128_c, 128, Bhat[128], 1, np.arange(4, 29)),
         ("crit128", G128_c, 128, Bhat[128], 2, np.arange(6, 29)),
         ("crit128", G128_c, 128, Bhat[128], 3, np.arange(8, 29)),
         ("crit96", G96_c, 96, Bhat[96], 2, np.arange(6, 25))]
for tag, G, L, B, R, ds in cases:
    F[(tag, R)] = (ds, np.array([F_fl_clusters(G, L, R, int(d), B)
                                 for d in ds]))
    print(f"  F_fl policzone: {tag} R={R} d in [{ds[0]},{ds[-1]}]")

# ------------------------------------------------------------------ QE-1
ok_sgn, ok_mono = True, True
qe1_windows = [("m02", 1, 8, 18), ("m02", 2, 10, 20),
               ("crit128", 1, 4, 28), ("crit128", 2, 6, 28),
               ("crit128", 3, 8, 28)]
for tag, R, d1, d2 in qe1_windows:
    ds, Fv = F[(tag, R)]
    sel = (ds >= d1) & (ds <= d2)
    s_ok = np.all(Fv[sel] < 0)
    m_ok = np.all(np.diff(Fv[sel]) > 0)
    ok_sgn &= s_ok; ok_mono &= m_ok
    print(f"    QE-1 {tag} R={R} [{d1},{d2}]: F<0: {s_ok}, mono: {m_ok}, "
          f"F({d1})={Fv[sel][0]:.3e}, F({d2})={Fv[sel][-1]:.3e}")
check("QE-1a: F_fl(d;R) < 0 wszedzie (uniwersalnosc przezywa rozciaglosc)",
      ok_sgn)
check("QE-1b: F_fl rosnace z d (przyciaganie) wszedzie", ok_mono)

# ------------------------------------------------------------------ QE-2
ok2 = True
mu = mu_lat(0.2)
for R, d1, d2 in [(1, 8, 18), (2, 10, 20)]:
    ds, Fv = F[("m02", R)]
    sel = (ds >= d1) & (ds <= d2)
    dsw = ds[sel].astype(float)
    _, s, r2 = linfit(dsw, np.log(np.abs(Fv[sel])) + 2 * np.log(dsw))
    kD = -s
    ok = abs(kD - 2 * mu) / (2 * mu) < 0.10 and r2 > 0.99
    ok2 &= ok
    print(f"    QE-2 R={R}: kappa_D={kD:.4f} vs 2mu={2*mu:.4f} "
          f"({100*abs(kD-2*mu)/(2*mu):.1f}%, R2={r2:.5f})")
check("QE-2: zasieg 2mu przezywa rozciaglosc (+/-10%, R2>0.99)", ok2)

# ------------------------------------------------------------------ QE-3
slopes = {}
amp = {}
for R, d1, d2 in [(1, 10, 28), (2, 14, 28), (3, 16, 28)]:
    ds, Fv = F[("crit128", R)]
    sel = (ds >= d1) & (ds <= d2)
    dsw = ds[sel].astype(float)
    _, s, r2 = linfit(np.log(dsw), np.log(np.abs(Fv[sel])))
    slopes[R] = (s, r2)
    # QE-5: amplituda przy narzuconym -2
    amp[R] = float(np.exp(np.mean(np.log(np.abs(Fv[sel])) + 2 * np.log(dsw))))
    print(f"    QE-3 R={R} okno[{d1},{d2}]: slope={s:.4f} (R2={r2:.6f})")
all_m2 = all(abs(s + 2.0) < 0.25 and r2 > 0.99 for s, r2 in slopes.values())
any_m1 = any(abs(s + 1.0) < 0.15 and r2 > 0.99 for s, r2 in slopes.values())
print(f"    QE-3: wszystkie slope=-2+/-0.25: {all_m2}; "
      f"ktorykolwiek slope=-1+/-0.15: {any_m1}")
check("QE-3: fit daleki wykonany, R2>0.99 wszedzie",
      all(r2 > 0.99 for _, r2 in slopes.values()),
      f"slopes = {dict((R, round(s, 3)) for R, (s, _) in slopes.items())}")

# ------------------------------------------------------------------ QE-4
run_found = {}
for R in [1, 2, 3]:
    ds, Fv = F[("crit128", R)]
    dm, pl = local_slopes(ds.astype(float), np.abs(Fv))
    hits = np.abs(pl - (-1.0)) <= 0.15
    best = cur = 0
    for h in hits:
        cur = cur + 1 if h else 0
        best = max(best, cur)
    run_found[R] = best
    prof = ", ".join(f"{int(d)}:{p:.2f}" for d, p in zip(dm, pl))
    print(f"    QE-4 R={R} profil p_loc(d): {prof}")
    print(f"    QE-4 R={R}: najdluzszy przebieg |p+1|<=0.15: {best}")
qe4_yes = any(b >= 3 for b in run_found.values())
print(f"    QE-4: przebieg >=3 kolejnych d z p_loc=-1+/-0.15: {qe4_yes}")
check("QE-4: profil p_loc policzony na pelnym zakresie (raport bez selekcji)",
      all((("crit128", R) in F) for R in [1, 2, 3]))

# ------------------------------------------------------------------ QE-5
print("    QE-5 (deskryptywnie): amplituda A(R) przy |F|~A/d^2 i pojemnosci:")
for R in [1, 2, 3]:
    sites = ball_sites(R)
    S = cov_matrix(G128_c, 128, sites, Bhat[128])
    CR = float(np.ones(len(sites)) @ np.linalg.solve(S, np.ones(len(sites))))
    pred = 0.5 * CR * CR / (4 * np.pi) ** 2
    print(f"      R={R}: A(R)={amp[R]:.5f}, C_R={CR:.4f}, "
          f"monopol 0.5*C^2/(4pi)^2={pred:.5f}, A/pred={amp[R]/pred:.3f}")

# -------------------------------------------------------------- kontrole
ds128, F128 = F[("crit128", 2)]
ds96, F96 = F[("crit96", 2)]
sel128 = (ds128 >= 14) & (ds128 <= 24)
sel96 = (ds96 >= 14) & (ds96 <= 24)
_, s128, _ = linfit(np.log(ds128[sel128].astype(float)),
                    np.log(np.abs(F128[sel128])))
_, s96, _ = linfit(np.log(ds96[sel96].astype(float)),
                   np.log(np.abs(F96[sel96])))
drift = abs(s128 - s96)
check("QE-kontrola(c): dryf slope L=96 vs L=128 (R=2, [14,24]) < 0.1",
      drift < 0.1, f"s128={s128:.4f}, s96={s96:.4f}, dryf={drift:.4f}")
allF = np.concatenate([Fv for _, Fv in F.values()])
check("QE-kontrola(b): F_fl <= 0 wszedzie (Fischer)", np.all(allF < 0),
      f"max F = {allF.max():.3e}")
check("QE-kontrola(d): G_m>0 (m=0.2)", G128_m02.min() > 0)

# ---------------------------------------------------------------- SUMMARY
print("-" * 72)
npass = sum(results.values()); ntot = len(results)
print(f"SUMMARY Phase 2: {npass}/{ntot} PASS")
for k, ok in results.items():
    print(f"  {'PASS' if ok else 'FAIL'}  {k}")
ctrl = all(results[k] for k in results if k.startswith(("QE-0", "QE-kontrola")))
if not ctrl:
    print("WERDYKT QE: INCONCLUSIVE (pad kontroli maszynerii)")
else:
    if any_m1 or qe4_yes:
        print("WERDYKT QE (pytanie N1): TAK — istnieje rezim wykladnika -1"
              f" (QE-3 any={any_m1}, QE-4 przebieg={qe4_yes})")
    else:
        print("WERDYKT QE (pytanie N1): NIE — wykladnik dalekiego pola "
              f"pozostaje -2 (slopes: "
              f"{dict((R, round(s, 3)) for R, (s, _) in slopes.items())}), "
              "R wchodzi tylko w amplitude; brak przebiegu p_loc=-1.")
