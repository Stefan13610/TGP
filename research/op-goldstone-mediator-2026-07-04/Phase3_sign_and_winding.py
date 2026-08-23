# Phase 3: dyskryminator znaku (G5) + skalowanie z kretem (G6)
# Implementacja SCISLE wg Phase0_balance.md (LOCKED, sekcje 3-5 + D7, D8).
#
# sign_test (D7): konfiguracja neutralna 4-wirowa (+1,+1,-1,-1);
#   pary rownoimienne d_ss=8, grupy przeciwne w L/2 (sily miedzygrupowe
#   znosza sie z symetrii) -> czy rownoimienne sie ODPYCHAJA?
# winding (D8): scan (+2,-2) vs baseline (+1,-1) TA SAMA procedura ->
#   ratio zboczy C2, predykcja n1*n2 = 4, pasmo G6 [3.2, 4.8].

import numpy as np

# ---------------- parametry LOCKED (identyczne z Phase 2) ----------------
kappa = 0.50
a, b, c = 0.50, 1.60, 1.00
dt = 0.02
EPS = 0.30
S_GUARD = 10.0
REC_EVERY = 10
H_TOL = 1e-10
TAU_RELAX = 200
CHECK_AT = 150
D_GRID = [6.0, 8.0, 12.0, 16.0, 24.0]
TWO_PI = 2.0 * np.pi
N, L = 256, 128.0
dx = L / N

disc = np.sqrt(b * b - 4.0 * a * c)
s_star = (b + disc) / (2.0 * c)
Phi_star = s_star ** 2
Vpp_star = a - 2.0 * b * s_star + 3.0 * c * s_star ** 2
xi = np.sqrt(kappa / Vpp_star)
J = kappa * Phi_star
C2_PRED = TWO_PI * J

coords = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
Z = X + 1j * Y


def Vpot(m):
    return 0.5 * a * m**2 - (b / 3.0) * m**3 + 0.25 * c * m**4


def energy(psi):
    gx = (np.roll(psi, -1, 0) - psi) / dx
    gy = (np.roll(psi, -1, 1) - psi) / dx
    m = np.abs(psi)
    return float(np.sum(0.5 * kappa * (np.abs(gx)**2 + np.abs(gy)**2)
                        + Vpot(m)) * dx * dx)


def euler_step(psi):
    m = np.abs(psi)
    lap = (np.roll(psi, 1, 0) + np.roll(psi, -1, 0)
           + np.roll(psi, 1, 1) + np.roll(psi, -1, 1) - 4.0 * psi) / (dx * dx)
    return psi + dt * (kappa * lap - psi * (a - b * m + c * m * m))


def make_config(vortices):
    th = np.zeros_like(X)
    A = np.full_like(X, s_star)
    for (x0, y0, n) in vortices:
        for k in range(-2, 3):
            th = th + n * np.angle(np.sinh(np.pi * (Z - (x0 + k * L) - 1j * y0) / L))
        ddx = (X - x0 + 0.5 * L) % L - 0.5 * L
        ddy = (Y - y0 + 0.5 * L) % L - 0.5 * L
        rr = np.sqrt(ddx * ddx + ddy * ddy)
        A = A * rr / np.sqrt(rr * rr + 2.0 * xi * xi)
    return A * np.exp(1j * th)


def wrap(x):
    return (x + np.pi) % TWO_PI - np.pi


def vortex_clusters(psi):
    th = np.angle(psi)
    dthx = wrap(np.roll(th, -1, 0) - th)
    dthy = wrap(np.roll(th, -1, 1) - th)
    W = (dthx + np.roll(dthy, -1, 0) - np.roll(dthx, -1, 1) - dthy) / TWO_PI
    cells = [tuple(cc) for cc in np.argwhere(np.abs(W) > 0.5)]
    if not cells:
        return []
    clusters = []
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
        clusters.append(group)
    Phi = np.abs(psi) ** 2
    out = []
    for group in clusters:
        n_tot = float(sum(round(W[cc]) for cc in group))
        i0, j0 = group[0]
        ri = np.mean([((cc[0] - i0 + N // 2) % N - N // 2) for cc in group])
        rj = np.mean([((cc[1] - j0 + N // 2) % N - N // 2) for cc in group])
        ic, jc = int(round(i0 + ri)) % N, int(round(j0 + rj)) % N
        off = np.arange(-6, 7)
        ii = (ic + off) % N
        jj = (jc + off) % N
        Pw = Phi[np.ix_(ii, jj)]
        wgt = np.clip(0.5 * Phi_star - Pw, 0.0, None)
        OX, OY = np.meshgrid(off * dx, off * dx, indexing="ij")
        if wgt.sum() > 0:
            sx = float(np.sum(wgt * OX) / wgt.sum())
            sy = float(np.sum(wgt * OY) / wgt.sum())
        else:
            sx = sy = 0.0
        out.append((coords[ic] + 0.5 * dx + sx, coords[jc] + 0.5 * dx + sy, n_tot))
    return out


def min_image_dist(p, q):
    ddx = (p[0] - q[0] + 0.5 * L) % L - 0.5 * L
    ddy = (p[1] - q[1] + 0.5 * L) % L - 0.5 * L
    return float(np.hypot(ddx, ddy))


def sign_group_positions(cls):
    """Zwraca (pozycje_plus, pozycje_minus) jako listy (x, y)."""
    pos_p = [(x, y) for (x, y, n) in cls if n > 0]
    pos_m = [(x, y) for (x, y, n) in cls if n < 0]
    return pos_p, pos_m


def pair_distance(cls):
    """d_meas miedzy centroidami znakow (waga |n|), min-image."""
    out = {}
    for sign in (+1, -1):
        sel = [(x, y, n) for (x, y, n) in cls if np.sign(n) == sign]
        if not sel:
            return np.nan
        x0, y0, _ = sel[0]
        wsum = xacc = yacc = 0.0
        for (x, y, n) in sel:
            rx = (x - x0 + 0.5 * L) % L - 0.5 * L
            ry = (y - y0 + 0.5 * L) % L - 0.5 * L
            xacc += abs(n) * rx
            yacc += abs(n) * ry
            wsum += abs(n)
        out[sign] = (x0 + xacc / wsum, y0 + yacc / wsum)
    return min_image_dist(out[+1], out[-1])


def evolve(psi0, max_steps, snap_at=(), track=None, stop_d=None):
    psi = psi0.copy()
    H0 = energy(psi)
    H_prev, H_inc_max = H0, 0.0
    flags = {"nonfinite": False, "runaway": False}
    snaps, tracks = {}, []
    for stp in range(1, max_steps + 1):
        psi = euler_step(psi)
        if stp % REC_EVERY:
            continue
        if not np.all(np.isfinite(psi)):
            flags["nonfinite"] = True
            break
        if np.max(np.abs(psi)) > S_GUARD:
            flags["runaway"] = True
            break
        H = energy(psi)
        if H > H_prev:
            H_inc_max = max(H_inc_max, H - H_prev)
        H_prev = H
        cls = vortex_clusters(psi)
        if track is not None:
            tracks.append((stp * dt, track(cls)))
        if stp in snap_at:
            snaps[stp] = {"E": H, "d": pair_distance(cls), "cls": cls}
        if stop_d is not None:
            d = pair_distance(cls)
            if not np.isfinite(d) or d < stop_d or not cls:
                break
    flags["H_inc_max"] = H_inc_max
    flags["H_mono"] = H_inc_max <= H_TOL * max(1.0, abs(H0))
    return psi, flags, snaps, tracks


def fit_L(dv, Ev):
    A_ = np.vstack([np.ones(len(dv)), np.log(dv)]).T
    coef, *_ = np.linalg.lstsq(A_, Ev, rcond=None)
    rss = float(np.sum((Ev - A_ @ coef) ** 2))
    return {"C1": coef[0], "C2": coef[1], "rss": rss}


def static_scan(n_wind, label):
    rows, fls = [], []
    for d0 in D_GRID:
        psi0 = make_config([(-0.5 * d0, 0.0, +n_wind), (+0.5 * d0, 0.0, -n_wind)])
        _, fl, snaps, _ = evolve(psi0, TAU_RELAX, snap_at=(CHECK_AT, TAU_RELAX))
        s200 = snaps[TAU_RELAX]
        n_tot = sum(n for (_, _, n) in s200["cls"])
        n_plus = sum(n for (_, _, n) in s200["cls"] if n > 0)
        rows.append({"d0": d0, "d": s200["d"], "E": s200["E"],
                     "n_tot": n_tot, "n_plus": n_plus})
        fls.append(fl)
        print(f"  [{label}] d0={d0:5.1f}  d_meas={s200['d']:7.4f}  "
              f"E={s200['E']:.4f}  kret(+)={n_plus:+.0f} n_tot={n_tot:+.0f}  "
              f"H_mono={fl['H_mono']}")
    return rows, fls


print("=== Phase 3: znak i kret (LOCK: Phase0_balance.md) ===")
print(f"N={N} L={L} dx={dx} dt={dt}  s*={s_star:.6f} xi={xi:.4f} "
      f"J={J:.4f}  C2_pred={C2_PRED:.4f}")
print()

# ---------------- A. sign_test (G5, D7) ----------------
print("--- sign_test (D7): (+1,+1) w x=-32, (-1,-1) w x=+32, d_ss=8, tau=200 ---")
VORT4 = [(-32.0, -4.0, +1), (-32.0, +4.0, +1),
         (+32.0, -4.0, -1), (+32.0, +4.0, -1)]


def track_ss(cls):
    pos_p, pos_m = sign_group_positions(cls)
    dp = min_image_dist(pos_p[0], pos_p[1]) if len(pos_p) == 2 else np.nan
    dm = min_image_dist(pos_m[0], pos_m[1]) if len(pos_m) == 2 else np.nan
    return dp, dm


psi0 = make_config(VORT4)
_, fl_sign, _, tracks = evolve(psi0, 10000, track=track_ss)
taus = np.array([t for (t, _) in tracks])
dss_p = np.array([v[0] for (_, v) in tracks])
dss_m = np.array([v[1] for (_, v) in tracks])
for tq in (0.2, 25.0, 50.0, 100.0, 150.0, 200.0):
    k = int(np.argmin(np.abs(taus - tq)))
    print(f"  tau={taus[k]:6.1f}   d_ss(+,+)={dss_p[k]:7.4f}   "
          f"d_ss(-,-)={dss_m[k]:7.4f}")
delta_p = float(dss_p[-1] - dss_p[0])
delta_m = float(dss_m[-1] - dss_m[0])
print(f"  Delta d_ss(+,+) = {delta_p:+.4f}   Delta d_ss(-,-) = {delta_m:+.4f}   "
      f"(prog D7: +/-1.0)")
if delta_p >= 1.0 and delta_m >= 1.0:
    g5 = "LADUNKOWY (rownoimienne sie ODPYCHAJA)"
elif delta_p <= -1.0 and delta_m <= -1.0:
    g5 = "UNIWERSALNY (rownoimienne sie przyciagaja — gravity-like)"
else:
    g5 = "NIEROZSTRZYGNIETE"
print(f"  G5 klasyfikator: {g5}")
print(f"  flagi: nf={fl_sign['nonfinite']} run={fl_sign['runaway']} "
      f"H_mono={fl_sign['H_mono']}")
print()

# ---------------- B. winding (G6, D8) ----------------
print("--- winding (D8): scan (+2,-2) vs baseline (+1,-1), ta sama procedura ---")
rows1, fls1 = static_scan(1, "(+1,-1)")
rows2, fls2 = static_scan(2, "(+2,-2)")
d1 = np.array([r["d"] for r in rows1]); E1 = np.array([r["E"] for r in rows1])
d2 = np.array([r["d"] for r in rows2]); E2 = np.array([r["E"] for r in rows2])
f1 = fit_L(d1, E1)
f2 = fit_L(d2, E2)
ratio = f2["C2"] / f1["C2"]
print(f"  C2(+1,-1) = {f1['C2']:.4f}   C2(+2,-2) = {f2['C2']:.4f}   "
      f"ratio = {ratio:.4f}  (predykcja |n1*n2| = 4)")
wind_ok1 = all(r["n_tot"] == 0 and r["n_plus"] == +1 for r in rows1)
wind_ok2 = all(r["n_tot"] == 0 and r["n_plus"] == +2 for r in rows2)
print(f"  kret zachowany: (+1,-1): {wind_ok1};  (+2,-2): n_plus=+2 wszedzie: "
      f"{wind_ok2}")
print()

# ---------------- C. obserwacja: rozszczepienie rdzenia n=2 ----------------
print("--- obserwacja (poza kryterium): stabilnosc rdzenia n=2, d0=16, tau=100 ---")
psi0 = make_config([(-8.0, 0.0, +2), (+8.0, 0.0, -2)])


def track_split(cls):
    pos_p, _ = sign_group_positions(cls)
    if len(pos_p) < 2:
        return len(pos_p), 0.0
    dmax = max(min_image_dist(p, q)
               for i, p in enumerate(pos_p) for q in pos_p[i + 1:])
    return len(pos_p), dmax


_, fl_split, _, tr_split = evolve(psi0, 5000, track=track_split)
for tq in (0.2, 10.0, 25.0, 50.0, 75.0, 100.0):
    tt = np.array([t for (t, _) in tr_split])
    k = int(np.argmin(np.abs(tt - tq)))
    ncl, dmax = tr_split[k][1]
    print(f"  tau={tt[k]:6.1f}   klastrow(+): {ncl}   "
          f"max separacja wewnetrzna: {dmax:.3f}")
print(f"  flagi: nf={fl_split['nonfinite']} run={fl_split['runaway']} "
      f"H_mono={fl_split['H_mono']}")
print()

# ---------------- EWALUACJA ----------------
print("=== EWALUACJA (LOCKED) ===")
all_fl = [fl_sign] + fls1 + fls2 + [fl_split]
g1_ok = all(not (f["nonfinite"] or f["runaway"]) and f["H_mono"] for f in all_fl)
print(f"G1 (czesc Phase3): finite/no-runaway + H_Gamma nierosnace: "
      f"{'PASS' if g1_ok else 'FAIL'}")
print(f"G5 (klasyfikator znaku, bez PASS/FAIL): {g5}")
g6 = (3.2 <= ratio <= 4.8) and wind_ok1 and wind_ok2
print(f"G6 (skalowanie z kretem): ratio={ratio:.4f} in [3.2, 4.8]: "
      f"{3.2 <= ratio <= 4.8}; kret zachowany: {wind_ok1 and wind_ok2}  "
      f"-> {'PASS' if g6 else 'FAIL'}")
print()
print(f"PODSUMOWANIE PHASE 3: G1(czesc)={'PASS' if g1_ok else 'FAIL'} "
      f"G5={g5.split(' ')[0]} G6={'PASS' if g6 else 'FAIL'}")
