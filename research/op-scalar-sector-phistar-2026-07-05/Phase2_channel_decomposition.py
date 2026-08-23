# Phase 2: dekompozycja kanalow kwadrupolem S/O — G3, G4, G5, G6 (+P5)
# Implementacja SCISLE wg Phase0_balance.md (LOCKED, sekcje 2-4, P1-P5).

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt = 0.02
S_GUARD = 10.0
H_TOL = 1e-10
TAU_RELAX = 200
CHECK_AT = 150
TWO_PI = 2.0 * np.pi
SIGNS = {"S": (+1, +1, -1, -1), "O": (+1, -1, -1, +1)}
D_GRID = [5.0, 6.0, 8.0, 10.0, 12.0]
D_GRID_CTRL = [5.0, 6.0, 8.0]
B_PRED = 13.5460          # Phase0_output (glowna); kontrola: 13.4861
B_PRED_CTRL = 13.4861

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)
Phi_star = s_star ** 2
Vpp_star = a - 2.0 * b_par * s_star + 3.0 * c * s_star ** 2
xi2 = kappa / Vpp_star
xi = np.sqrt(xi2)
J = kappa * Phi_star
C2 = TWO_PI * J
m_a = np.sqrt(Vpp_star / kappa)
M_BAND = (0.8 * m_a, 1.2 * m_a)


def Vpot(m):
    return 0.5 * a * m**2 - (b_par / 3.0) * m**3 + 0.25 * c * m**4


def make_grid(N, L):
    dxg = L / N
    coords = (np.arange(N) - N // 2) * dxg
    X, Y = np.meshgrid(coords, coords, indexing="ij")
    # sieciowa funkcja Greena: Lap G = delta - mean (5-pkt)
    rho = np.zeros((N, N))
    rho[0, 0] = 1.0 / (dxg * dxg)
    kx = TWO_PI * np.fft.fftfreq(N, d=dxg)
    KX, KY = np.meshgrid(kx, kx, indexing="ij")
    lam = -((2 - 2 * np.cos(KX * dxg)) + (2 - 2 * np.cos(KY * dxg))) / (dxg * dxg)
    lam[0, 0] = 1.0
    Gh = np.fft.fft2(rho) / lam
    Gh[0, 0] = 0.0
    G_lat = np.real(np.fft.ifft2(Gh))
    return {"N": N, "L": L, "dx": dxg, "X": X, "Y": Y, "Z": X + 1j * Y,
            "coords": coords, "G": G_lat}


def energy(psi, g):
    dxg = g["dx"]
    gx = (np.roll(psi, -1, 0) - psi) / dxg
    gy = (np.roll(psi, -1, 1) - psi) / dxg
    return float(np.sum(0.5 * kappa * (np.abs(gx)**2 + np.abs(gy)**2)
                        + Vpot(np.abs(psi))) * dxg * dxg)


def euler_step(psi, g):
    dxg = g["dx"]
    m = np.abs(psi)
    lp = (np.roll(psi, 1, 0) + np.roll(psi, -1, 0)
          + np.roll(psi, 1, 1) + np.roll(psi, -1, 1) - 4.0 * psi) / (dxg * dxg)
    return psi + dt * (kappa * lp - psi * (a - b_par * m + c * m * m))


def make_cfg(g, vortices):
    L = g["L"]
    th = np.zeros_like(g["X"])
    A = np.full_like(g["X"], s_star)
    for (x0, y0, n) in vortices:
        for k in range(-2, 3):
            th = th + n * np.angle(
                np.sinh(np.pi * (g["Z"] - (x0 + k * L) - 1j * y0) / L))
        ddx = (g["X"] - x0 + 0.5 * L) % L - 0.5 * L
        ddy = (g["Y"] - y0 + 0.5 * L) % L - 0.5 * L
        rr = np.sqrt(ddx * ddx + ddy * ddy)
        A = A * rr / np.sqrt(rr * rr + 2.0 * xi2)
    return A * np.exp(1j * th)


def wrap(x):
    return (x + np.pi) % TWO_PI - np.pi


def vortex_clusters(psi, g):
    N, dxg, coords = g["N"], g["dx"], g["coords"]
    th = np.angle(psi)
    dthx = wrap(np.roll(th, -1, 0) - th)
    dthy = wrap(np.roll(th, -1, 1) - th)
    W = (dthx + np.roll(dthy, -1, 0) - np.roll(dthx, -1, 1) - dthy) / TWO_PI
    cells = [tuple(cc) for cc in np.argwhere(np.abs(W) > 0.5)]
    Phi = np.abs(psi) ** 2
    out = []
    unassigned = set(cells)
    while unassigned:
        seed_cell = unassigned.pop()
        group = [seed_cell]
        grew = True
        while grew:
            grew = False
            for cand in list(unassigned):
                for gc in group:
                    di = abs((cand[0] - gc[0] + N // 2) % N - N // 2)
                    dj = abs((cand[1] - gc[1] + N // 2) % N - N // 2)
                    if max(di, dj) <= 3:
                        group.append(cand)
                        unassigned.discard(cand)
                        grew = True
                        break
        n_tot = float(sum(round(W[cc]) for cc in group))
        i0, j0 = group[0]
        ri = np.mean([((cc[0] - i0 + N // 2) % N - N // 2) for cc in group])
        rj = np.mean([((cc[1] - j0 + N // 2) % N - N // 2) for cc in group])
        ic, jc = int(round(i0 + ri)) % N, int(round(j0 + rj)) % N
        off = np.arange(-6, 7)
        ii, jj = (ic + off) % N, (jc + off) % N
        Pw = Phi[np.ix_(ii, jj)]
        wgt = np.clip(0.5 * Phi_star - Pw, 0.0, None)
        OX, OY = np.meshgrid(off * dxg, off * dxg, indexing="ij")
        sx = float(np.sum(wgt * OX) / wgt.sum()) if wgt.sum() > 0 else 0.0
        sy = float(np.sum(wgt * OY) / wgt.sum()) if wgt.sum() > 0 else 0.0
        out.append((coords[ic] + 0.5 * dxg + sx, coords[jc] + 0.5 * dxg + sy, n_tot))
    return out


def G_at(g, ddx, ddy):
    N, L, dxg = g["N"], g["L"], g["dx"]
    fx, fy = (ddx % L) / dxg, (ddy % L) / dxg
    i0, j0 = int(np.floor(fx)) % N, int(np.floor(fy)) % N
    tx, ty = fx - np.floor(fx), fy - np.floor(fy)
    i1, j1 = (i0 + 1) % N, (j0 + 1) % N
    GL = g["G"]
    return ((1 - tx) * (1 - ty) * GL[i0, j0] + tx * (1 - ty) * GL[i1, j0]
            + (1 - tx) * ty * GL[i0, j1] + tx * ty * GL[i1, j1])


def gph_of(g, cls):
    s = 0.0
    for i in range(len(cls)):
        for j in range(i + 1, len(cls)):
            xi_, yi, ni = cls[i]
            xj, yj, nj = cls[j]
            s += ni * nj * G_at(g, xi_ - xj, yi - yj)
    return -TWO_PI * s


def mindist(g, p, q):
    L = g["L"]
    ddx = (p[0] - q[0] + 0.5 * L) % L - 0.5 * L
    ddy = (p[1] - q[1] + 0.5 * L) % L - 0.5 * L
    return float(np.hypot(ddx, ddy))


def measure(psi, g):
    cls = vortex_clusters(psi, g)
    ok = (len(cls) == 4 and all(abs(n) == 1 for (_, _, n) in cls)
          and sum(n for (_, _, n) in cls) == 0)
    left = [cc for cc in cls if cc[0] < 0]
    right = [cc for cc in cls if cc[0] >= 0]
    if len(left) == 2 and len(right) == 2:
        dL = mindist(g, left[0][:2], left[1][:2])
        dR = mindist(g, right[0][:2], right[1][:2])
        d_cfg = 0.5 * (dL + dR)
    else:
        d_cfg = np.nan
        ok = False
    return {"E": energy(psi, g), "gph": gph_of(g, cls), "d": d_cfg, "ok": ok}


def relax_measure(g, X0, d0, sg):
    pos = [(-X0, +0.5 * d0), (-X0, -0.5 * d0), (+X0, +0.5 * d0), (+X0, -0.5 * d0)]
    psi = make_cfg(g, [(x, y, n) for (x, y), n in zip(pos, sg)])
    H_prev, H_inc = energy(psi, g), 0.0
    fl = {"nonfinite": False, "runaway": False}
    snaps = {}
    for stp in range(1, TAU_RELAX + 1):
        psi = euler_step(psi, g)
        if stp % 10 == 0:
            if not np.all(np.isfinite(psi)):
                fl["nonfinite"] = True
                break
            if np.max(np.abs(psi)) > S_GUARD:
                fl["runaway"] = True
                break
            H = energy(psi, g)
            if H > H_prev:
                H_inc = max(H_inc, H - H_prev)
            H_prev = H
            if stp in (CHECK_AT, TAU_RELAX):
                snaps[stp] = measure(psi, g)
    fl["H_mono"] = H_inc <= H_TOL
    return snaps, fl


def joint_fit(rows, form, m_grid=None):
    """rows: list of dict(E, gph, d). form: 'P' albo 'Y'."""
    E = np.array([r["E"] for r in rows])
    gph = np.array([r["gph"] for r in rows])
    d = np.array([r["d"] for r in rows])
    n = len(E)
    if form == "P":
        A_ = np.vstack([np.ones(n), gph, -np.log(d / xi) / d**2]).T
        coef, *_ = np.linalg.lstsq(A_, E, rcond=None)
        rss = float(np.sum((E - A_ @ coef) ** 2))
        return {"K": coef[0], "Cphi": coef[1], "B": coef[2], "rss": rss,
                "aic": n * np.log(max(rss, 1e-30) / n) + 2 * 3, "k": 3}
    best = None
    for m in m_grid:
        A_ = np.vstack([np.ones(n), gph, -np.exp(-m * d)]).T
        coef, *_ = np.linalg.lstsq(A_, E, rcond=None)
        rss = float(np.sum((E - A_ @ coef) ** 2))
        if best is None or rss < best["rss"]:
            best = {"K": coef[0], "Cphi": coef[1], "A": coef[2], "m": m, "rss": rss}
    best["aic"] = n * np.log(max(best["rss"], 1e-30) / n) + 2 * 4
    best["k"] = 4
    return best


def res_sym(rows, fit):
    """P2: E_res_sym(d0) = srednia rezyduow S,O (bez skladnika amp)."""
    out = {}
    for r in rows:
        res = r["E"] - fit["K"] - fit["Cphi"] * r["gph"]
        out.setdefault(r["d0"], []).append(res)
    return {d0: float(np.mean(v)) for d0, v in sorted(out.items())}


def run_scan(g, X0, d_list, tag):
    print(f"--- skan S/O [{tag}]: X={X0}, d0 in {d_list} ---")
    rows200, rows150, fls = [], [], []
    for d0 in d_list:
        for cfg, sg in SIGNS.items():
            snaps, fl = relax_measure(g, X0, d0, sg)
            fls.append(fl)
            m2, m1 = snaps[TAU_RELAX], snaps[CHECK_AT]
            rows200.append({**m2, "d0": d0, "cfg": cfg})
            rows150.append({**m1, "d0": d0, "cfg": cfg})
            print(f"  d0={d0:4.1f} {cfg}: d_meas={m2['d']:7.4f}  "
                  f"E={m2['E']:.5f}  gph={m2['gph']:+9.4f}  kret_ok={m2['ok']}  "
                  f"H_mono={fl['H_mono']}")
    return rows200, rows150, fls


print("=== Phase 2: dekompozycja kanalow (LOCK: Phase0_balance.md) ===")
print(f"s*={s_star:.6f} xi={xi:.4f} 2piJ={C2:.4f} m_a={m_a:.4f} "
      f"B_pred={B_PRED} (Phase0)")
print()

G_MAIN = make_grid(256, 128.0)
rows200, rows150, fls = run_scan(G_MAIN, 32.0, D_GRID, "GLOWNA dx=0.5")
print()

fP = joint_fit(rows200, "P")
fY = joint_fit(rows200, "Y", np.linspace(M_BAND[0], M_BAND[1], 801))
fP150 = joint_fit(rows150, "P")
dAIC = fY["aic"] - fP["aic"]
print("--- fit laczny (P1): E = K + Cphi*gph + amp(d) ---")
print(f"  (P): K={fP['K']:.4f}  Cphi={fP['Cphi']:.4f}  B={fP['B']:.4f}  "
      f"RSS={fP['rss']:.3e}  AIC={fP['aic']:+.2f}")
print(f"  (Y): K={fY['K']:.4f}  Cphi={fY['Cphi']:.4f}  A={fY['A']:.4f}  "
      f"m={fY['m']:.4f}  RSS={fY['rss']:.3e}  AIC={fY['aic']:+.2f}")
print(f"  Delta_AIC (Y - P) = {dAIC:+.2f}  (>= +10: wygrywa P; <= -10: Y)")
print()

rs200 = res_sym(rows200, fP)
rs150 = res_sym(rows150, fP150)
floor = max(abs(rs200[d0] - rs150[d0]) for d0 in rs200)
vals = [rs200[d0] for d0 in sorted(rs200)]
rng = max(vals) - min(vals)
diffs = np.diff(vals)
spearman_full = bool(np.all(diffs > 0))
print("--- P2: przebieg kanalu slepego na znak (rezyduum bezmodelowe) ---")
print("   d0     E_res_sym(200)   E_res_sym(150)   |roznica|")
for d0 in sorted(rs200):
    print(f"  {d0:4.1f}   {rs200[d0]:+.5f}       {rs150[d0]:+.5f}       "
          f"{abs(rs200[d0]-rs150[d0]):.5f}")
print(f"  floor = {floor:.5f}   rozstep = {rng:.5f}   "
      f"rozstep/floor = {rng/max(floor,1e-12):.1f}")
print()

# ---------------- G6: kontrola dx=0.25 ----------------
G_CTRL = make_grid(256, 64.0)
rows200c, rows150c, flsc = run_scan(G_CTRL, 16.0, D_GRID_CTRL, "KONTROLA dx=0.25")
fPc = joint_fit(rows200c, "P")
rsc = res_sym(rows200c, fPc)
valsc = [rsc[d0] for d0 in sorted(rsc)]
mono_c = bool(np.all(np.diff(valsc) > 0))
print(f"  kontrola: Cphi={fPc['Cphi']:.4f}  B={fPc['B']:.4f} "
      f"(B_pred_ctrl={B_PRED_CTRL})  E_res_sym rosnace: {mono_c}")
print()

# ---------------- P5: obserwacja dynamiczna ----------------
print("--- P5 (obserwacja): kwadrupol S, d0=4, dynamika tau<=100 ---")
pos = [(-32.0, +2.0), (-32.0, -2.0), (+32.0, +2.0), (+32.0, -2.0)]
psi = make_cfg(G_MAIN, [(x, y, n) for (x, y), n in zip(pos, SIGNS["S"])])
recs = []
fl5 = {"nonfinite": False, "runaway": False}
for stp in range(1, 5001):
    psi = euler_step(psi, G_MAIN)
    if stp % 50 == 0:
        m = measure(psi, G_MAIN)
        recs.append((stp * dt, m["d"]))
for (t, dv) in recs[::20] + [recs[-1]]:
    print(f"  tau={t:6.1f}  d_ss={dv:7.4f}")
ds = np.array([dv for (_, dv) in recs])
ts = np.array([t for (t, _) in recs])
vd = np.gradient(ds, ts) * ds
print(f"  v*d: poczatek {vd[2]:.4f} -> koniec {vd[-3]:.4f} "
      f"(czysta faza: ~const; spadek na malych d = slad przyciagajacej "
      f"poprawki dylatacyjnej — raport jakosciowy)")
print()

# ---------------- EWALUACJA ----------------
print("=== EWALUACJA (LOCKED) ===")
g1 = (all(not (f["nonfinite"] or f["runaway"]) and f["H_mono"]
          for f in fls + flsc)
      and all(r["ok"] for r in rows200 + rows200c))
print(f"G1: flagi + H_mono + kret (4x |n|=1, n_tot=0) wszedzie: "
      f"{'PASS' if g1 else 'FAIL'}")
g3 = 0.8 <= fP["Cphi"] / C2 <= 1.2
print(f"G3 (dekompozycja dziala): Cphi/(2piJ) = {fP['Cphi']/C2:.4f} "
      f"in [0.8, 1.2]: {g3}  -> {'PASS' if g3 else 'FAIL'}")
g4 = spearman_full and rng > 10 * floor
print(f"G4 (RDZEN — uniwersalne przyciaganie): E_res_sym scisle rosnace: "
      f"{spearman_full}; rozstep > 10*floor: {rng > 10*floor} "
      f"({rng:.4f} vs {10*floor:.4f})  -> {'PASS' if g4 else 'FAIL'}")
if dAIC >= 10:
    winner = "P"
    band_ok = 0.5 <= fP["B"] / B_PRED <= 2.0
    print(f"G5: wygrywa (P); B/B_pred = {fP['B']/B_PRED:.4f} in [0.5, 2.0]: "
          f"{band_ok}  -> {'PASS' if band_ok else 'FAIL'}")
    g5 = band_ok
elif dAIC <= -10:
    winner = "Y"
    band_ok = M_BAND[0] <= fY["m"] <= M_BAND[1]
    print(f"G5: wygrywa (Y); m = {fY['m']:.4f} in {M_BAND}: {band_ok}  "
          f"-> {'PASS' if band_ok else 'FAIL'}")
    g5 = band_ok
else:
    winner = "BRAK"
    g5 = False
    print(f"G5: BRAK rozstrzygniecia formy (|Delta_AIC| < 10) — raport; "
          f"(P): B/B_pred={fP['B']/B_PRED:.4f}; (Y): m={fY['m']:.4f}")
g6 = mono_c and np.sign(valsc[-1] - valsc[0]) == np.sign(vals[-1] - vals[0])
print(f"G6 (kontrola dx=0.25): monotonia+znak zgodne: {'PASS' if g6 else 'FAIL'}")
print()
print(f"PODSUMOWANIE PHASE 2: G1={'PASS' if g1 else 'FAIL'} "
      f"G3={'PASS' if g3 else 'FAIL'} G4={'PASS' if g4 else 'FAIL'} "
      f"G5={'PASS' if g5 else ('BRAK' if winner == 'BRAK' else 'FAIL')} "
      f"G6={'PASS' if g6 else 'FAIL'}  [forma: {winner}]")
