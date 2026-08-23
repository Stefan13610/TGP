# Phase 2: prawo oddzialywania pary wirow — G3 (istnienie/stabilnosc), G4 (RDZEN)
# Implementacja SCISLE wg Phase0_balance.md (LOCKED, sekcje 1-5 + D1-D10).
#
# Scenariusze: vortex_core (profil #1), pair_static(d) (fit L vs Y, #2),
# pair_dynamic(d0) (v*d ~ const, C2 z drogi dynamicznej D6),
# kontrola dx=0.25 + test pinningu (D9).

import numpy as np

# ---------------- parametry LOCKED ----------------
kappa = 0.50
a, b, c = 0.50, 1.60, 1.00
dt = 0.02
EPS = 0.30
S_GUARD = 10.0
REC_EVERY = 10
VBASE = 2.0
KOFF = int(round(VBASE / (REC_EVERY * dt)))     # D5: baza Delta_tau = 2.0
H_TOL = 1e-10
TAU_RELAX = 200                                  # D3: dokladnie 200 krokow
CHECK_AT = 150                                   # D3: punkt kontrolny
D_GRID = [6.0, 8.0, 12.0, 16.0, 24.0]
D0_DYN = [12.0, 16.0, 24.0]
D_GRID_CTRL = [6.0, 8.0, 12.0]
TWO_PI = 2.0 * np.pi

disc = np.sqrt(b * b - 4.0 * a * c)
s_star = (b + disc) / (2.0 * c)
Phi_star = s_star ** 2
Vpp_star = a - 2.0 * b * s_star + 3.0 * c * s_star ** 2
xi = np.sqrt(kappa / Vpp_star)
J = kappa * Phi_star
C2_PRED = TWO_PI * J
BAND = (0.8 * C2_PRED, 1.2 * C2_PRED)
M_BAND = (0.3, 2.0)                              # D4: kryterialne pasmo m dla (Y)


def Vpot(m):
    return 0.5 * a * m**2 - (b / 3.0) * m**3 + 0.25 * c * m**4


def Vprime(m):
    return m * (a - b * np.abs(m) + c * m * m)


def make_grid(N, L):
    dxg = L / N
    coords = (np.arange(N) - N // 2) * dxg
    X, Y = np.meshgrid(coords, coords, indexing="ij")
    return {"N": N, "L": L, "dx": dxg, "X": X, "Y": Y, "coords": coords}


def energy(psi, g):
    dxg = g["dx"]
    gx = (np.roll(psi, -1, 0) - psi) / dxg
    gy = (np.roll(psi, -1, 1) - psi) / dxg
    m = np.abs(psi)
    return float(np.sum(0.5 * kappa * (np.abs(gx)**2 + np.abs(gy)**2)
                        + Vpot(m)) * dxg * dxg)


def euler_step(psi, g):
    dxg = g["dx"]
    m = np.abs(psi)
    lap = (np.roll(psi, 1, 0) + np.roll(psi, -1, 0)
           + np.roll(psi, 1, 1) + np.roll(psi, -1, 1) - 4.0 * psi) / (dxg * dxg)
    return psi + dt * (kappa * lap - psi * (a - b * m + c * m * m))


def theta_ansatz(g, vortices):
    """D1: theta = sum_i n_i sum_{k=-2..2} arg sinh(pi*(z - z_i - kL)/L)."""
    Z = g["X"] + 1j * g["Y"]
    L = g["L"]
    th = np.zeros_like(g["X"])
    for (x0, y0, n) in vortices:
        for k in range(-2, 3):
            th = th + n * np.angle(np.sinh(np.pi * (Z - (x0 + k * L) - 1j * y0) / L))
    return th


def amp_ansatz(g, vortices):
    A = np.full_like(g["X"], s_star)
    L = g["L"]
    for (x0, y0, n) in vortices:
        ddx = (g["X"] - x0 + 0.5 * L) % L - 0.5 * L
        ddy = (g["Y"] - y0 + 0.5 * L) % L - 0.5 * L
        rr = np.sqrt(ddx * ddx + ddy * ddy)
        A = A * rr / np.sqrt(rr * rr + 2.0 * xi * xi)
    return A


def make_pair(g, vortices):
    return amp_ansatz(g, vortices) * np.exp(1j * theta_ansatz(g, vortices))


def wrap(x):
    return (x + np.pi) % TWO_PI - np.pi


def vortex_clusters(psi, g):
    """D2: winding plakietowy -> klastry rdzeni -> pozycje (centroid deficytu Phi).
    Zwraca liste (x, y, n_tot)."""
    N, dxg, coords = g["N"], g["dx"], g["coords"]
    th = np.angle(psi)
    dthx = wrap(np.roll(th, -1, 0) - th)
    dthy = wrap(np.roll(th, -1, 1) - th)
    W = (dthx + np.roll(dthy, -1, 0) - np.roll(dthx, -1, 1) - dthy) / TWO_PI
    cells = [tuple(cc) for cc in np.argwhere(np.abs(W) > 0.5)]
    if not cells:
        return []
    # klastrowanie po sasiedztwie (Chebyshev <= 3 wezly, torus)
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
        # sredni indeks wzgledem pierwszej komorki (torus)
        ri = np.mean([((cc[0] - i0 + N // 2) % N - N // 2) for cc in group])
        rj = np.mean([((cc[1] - j0 + N // 2) % N - N // 2) for cc in group])
        ic, jc = int(round(i0 + ri)) % N, int(round(j0 + rj)) % N
        # centroid deficytu Phi w oknie +/-6 wezlow (D2)
        off = np.arange(-6, 7)
        ii = (ic + off) % N
        jj = (jc + off) % N
        Pw = Phi[np.ix_(ii, jj)]
        wgt = np.clip(0.5 * Phi_star - Pw, 0.0, None)
        OX, OY = np.meshgrid(off * dxg, off * dxg, indexing="ij")
        if wgt.sum() > 0:
            sx = float(np.sum(wgt * OX) / wgt.sum())
            sy = float(np.sum(wgt * OY) / wgt.sum())
        else:
            sx = sy = 0.0
        # srodek plakietki (i,j) lezy w coords + dx/2
        x = coords[ic] + 0.5 * dxg + sx
        y = coords[jc] + 0.5 * dxg + sy
        out.append((x, y, n_tot))
    return out


def pair_distance(cls, g):
    """d_meas = odleglosc min-image miedzy centroidami znakow (waga |n|)."""
    L = g["L"]
    pos = {+1: None, -1: None}
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
        pos[sign] = (x0 + xacc / wsum, y0 + yacc / wsum)
    ddx = (pos[+1][0] - pos[-1][0] + 0.5 * L) % L - 0.5 * L
    ddy = (pos[+1][1] - pos[-1][1] + 0.5 * L) % L - 0.5 * L
    return float(np.hypot(ddx, ddy))


def evolve(psi0, g, max_steps, stop_d=None, snap_at=()):
    """Ewolucja z metrologia D10; zapisy co REC_EVERY: (tau, d_meas, H).
    snap_at: kroki, w ktorych zapisywany pelny pomiar (E, d, klastry)."""
    psi = psi0.copy()
    H0 = energy(psi, g)
    H_prev, H_inc_max = H0, 0.0
    flags = {"nonfinite": False, "runaway": False}
    taus, dmeas, Hs = [], [], []
    snaps = {}
    for stp in range(1, max_steps + 1):
        psi = euler_step(psi, g)
        if stp % REC_EVERY:
            continue
        if not np.all(np.isfinite(psi)):
            flags["nonfinite"] = True
            break
        if np.max(np.abs(psi)) > S_GUARD:
            flags["runaway"] = True
            break
        H = energy(psi, g)
        if H > H_prev:
            H_inc_max = max(H_inc_max, H - H_prev)
        H_prev = H
        cls = vortex_clusters(psi, g)
        d = pair_distance(cls, g)
        taus.append(stp * dt); dmeas.append(d); Hs.append(H)
        if stp in snap_at:
            snaps[stp] = {"E": H, "d": d, "cls": cls}
        if stop_d is not None and (not np.isfinite(d) or d < stop_d or not cls):
            break
    flags["H_inc_max"] = H_inc_max
    flags["H_mono"] = H_inc_max <= H_TOL * max(1.0, abs(H0))
    return psi, np.array(taus), np.array(dmeas), np.array(Hs), flags, snaps


def fit_L(dv, Ev):
    A_ = np.vstack([np.ones(len(dv)), np.log(dv)]).T
    coef, *_ = np.linalg.lstsq(A_, Ev, rcond=None)
    rss = float(np.sum((Ev - A_ @ coef) ** 2))
    n = len(dv)
    aic = n * np.log(max(rss, 1e-30) / n) + 2 * 2
    return {"C1": coef[0], "C2": coef[1], "rss": rss, "aic": aic}


def fit_Y(dv, Ev, m_grid):
    best = None
    n = len(dv)
    for m in m_grid:
        A_ = np.vstack([np.ones(n), -np.exp(-m * dv)]).T
        coef, *_ = np.linalg.lstsq(A_, Ev, rcond=None)
        rss = float(np.sum((Ev - A_ @ coef) ** 2))
        if best is None or rss < best["rss"]:
            best = {"C1": coef[0], "A": coef[1], "m": m, "rss": rss}
    best["aic"] = n * np.log(max(best["rss"], 1e-30) / n) + 2 * 3
    return best


def solve_core_ode():
    dr = 0.05
    r = np.arange(1, int(25.0 / dr)) * dr
    f = s_star * r / np.sqrt(r * r + 2.0 * xi * xi)
    dt_ode = 0.4 * dr * dr / (2.0 * kappa)
    for it in range(400000):
        fl = np.concatenate(([0.0], f[:-1]))
        fr = np.concatenate((f[1:], [s_star]))
        d2 = (fl + fr - 2.0 * f) / (dr * dr)
        d1 = (fr - fl) / (2.0 * dr)
        rhs = kappa * (d2 + d1 / r - f / (r * r)) - Vprime(f)
        f = f + dt_ode * rhs
        if it % 2000 == 0 and np.max(np.abs(rhs)) < 1e-9:
            break
    return r, f


def radial_profile(psi, g, pos, rmax=10.0, dbin=0.5):
    L = g["L"]
    ddx = (g["X"] - pos[0] + 0.5 * L) % L - 0.5 * L
    ddy = (g["Y"] - pos[1] + 0.5 * L) % L - 0.5 * L
    rr = np.hypot(ddx, ddy).ravel()
    mm = np.abs(psi).ravel()
    edges = np.arange(0.0, rmax + dbin, dbin)
    cent, prof = [], []
    for k in range(len(edges) - 1):
        sel = (rr >= edges[k]) & (rr < edges[k + 1])
        if sel.sum():
            cent.append(0.5 * (edges[k] + edges[k + 1]))
            prof.append(float(mm[sel].mean()))
    return np.array(cent), np.array(prof)


def core_check(g, d_core, label, ode_r, ode_f):
    psi0 = make_pair(g, [(-0.5 * d_core, 0.0, +1), (+0.5 * d_core, 0.0, -1)])
    psi, _, _, _, fl, snaps = evolve(psi0, g, TAU_RELAX,
                                     snap_at=(CHECK_AT, TAU_RELAX))
    cls = snaps[TAU_RELAX]["cls"]
    ns = sorted(n for (_, _, n) in cls)
    pos_p = [(x, y) for (x, y, n) in cls if n > 0][0]
    rc, pc = radial_profile(psi, g, pos_p)
    print(f"  [{label}] para referencyjna d={d_core}: klastry n={ns} "
          f"(n_tot={sum(ns):.0f})  flagi: nf={fl['nonfinite']} "
          f"run={fl['runaway']} H_mono={fl['H_mono']}")
    print("     r     |psi|_zmierzone   f_ODE(r)    roznica")
    devs = []
    for rq in (0.75, 1.25, 2.25, 3.25, 4.25, 6.25, 8.25):
        mq = float(np.interp(rq, rc, pc))
        fq = float(np.interp(rq, ode_r, ode_f))
        devs.append(abs(mq - fq))
        print(f"   {rq:5.2f}     {mq:8.5f}      {fq:8.5f}    {mq-fq:+.5f}")
    m_far = float(np.interp(8.25, rc, pc))
    ok_star = abs(m_far / s_star - 1.0) <= 0.03
    print(f"     |psi|(r=8.25)/s* = {m_far/s_star:.4f} (prog: 1 +/- 0.03: {ok_star}); "
          f"max |roznica| vs ODE = {max(devs):.4f}")
    return ns, ok_star, max(devs), fl


def static_scan(g, d_list, label):
    print(f"--- pair_static [{label}]: relaksacja {TAU_RELAX} krokow (D3) ---")
    rows = []
    all_fl = []
    for d0 in d_list:
        psi0 = make_pair(g, [(-0.5 * d0, 0.0, +1), (+0.5 * d0, 0.0, -1)])
        _, _, _, _, fl, snaps = evolve(psi0, g, TAU_RELAX,
                                       snap_at=(CHECK_AT, TAU_RELAX))
        s150, s200 = snaps[CHECK_AT], snaps[TAU_RELAX]
        ns = sorted(n for (_, _, n) in s200["cls"])
        rows.append({"d0": d0, "d150": s150["d"], "E150": s150["E"],
                     "d200": s200["d"], "E200": s200["E"], "ns": ns})
        all_fl.append(fl)
        print(f"  d0={d0:5.1f}  d_meas(150)={s150['d']:7.4f} E(150)={s150['E']:.4f}"
              f"   d_meas(200)={s200['d']:7.4f} E(200)={s200['E']:.4f}"
              f"   kret={ns}  H_mono={fl['H_mono']}")
    return rows, all_fl


print("=== Phase 2: prawo pary wirow (LOCK: Phase0_balance.md) ===")
print(f"kappa={kappa} a={a} b={b} c={c} dt={dt} eps={EPS}  s*={s_star:.6f} "
      f"Phi*={Phi_star:.6f} xi={xi:.4f} J={J:.4f}")
print(f"predykcja: C2 = 2*pi*J = {C2_PRED:.4f}, pasmo G4 = "
      f"[{BAND[0]:.3f}, {BAND[1]:.3f}]; (Y) kryterialnie m in {M_BAND} (D4)")
print()

ode_r, ode_f = solve_core_ode()
G_MAIN = make_grid(256, 128.0)

# ---------------- A. vortex_core (G3) ----------------
print("--- vortex_core: profil rdzenia po relaksacji vs rownanie #1 ---")
ns_core, okstar_main, dev_main, fl_core = core_check(
    G_MAIN, 64.0, "GLOWNA dx=0.5", ode_r, ode_f)
print()

# ---------------- B. pair_static (G4 RDZEN) ----------------
rows, fl_stat = static_scan(G_MAIN, D_GRID, "GLOWNA dx=0.5")
d200 = np.array([r["d200"] for r in rows])
E200 = np.array([r["E200"] for r in rows])
d150 = np.array([r["d150"] for r in rows])
E150 = np.array([r["E150"] for r in rows])
wind_ok_stat = all(r["ns"] == [-1.0, 1.0] for r in rows)
print()

fL = fit_L(d200, E200)
fY = fit_Y(d200, E200, np.linspace(M_BAND[0], M_BAND[1], 1701))
fYfree = fit_Y(d200, E200, np.geomspace(0.005, 3.0, 2500))
fL150 = fit_L(d150, E150)
dAIC = fY["aic"] - fL["aic"]
print("--- fit (L) vs (Y) na E(d_meas) po 200 krokach (D4) ---")
print(f"  (L) log:     C1={fL['C1']:.4f}  C2={fL['C2']:.4f}  "
      f"RSS={fL['rss']:.3e}  AIC={fL['aic']:+.2f}")
print(f"  (Y) Yukawa (m in {M_BAND}): C1={fY['C1']:.4f}  A={fY['A']:.4f}  "
      f"m={fY['m']:.4f}  RSS={fY['rss']:.3e}  AIC={fY['aic']:+.2f}")
print(f"  Delta_AIC (Y - L) = {dAIC:+.2f}   (>= +10: wygrywa L; <= -10: wygrywa Y)")
print(f"  [czulosc D4] (Y) z m swobodnym [0.005, 3]: m={fYfree['m']:.4f} "
      f"RSS={fYfree['rss']:.3e} AIC={fYfree['aic']:+.2f} "
      f"(degeneracja m->0 zapisana z gory w D4, poza kryterium)")
print(f"  [czulosc D3] fit (L) na punktach z kroku 150: C2={fL150['C2']:.4f} "
      f"(stabilnosc transjentu: Delta C2 = {fL150['C2']-fL['C2']:+.4f})")
print()

# ---------------- C. pair_dynamic (D5, D6) ----------------
print("--- pair_dynamic: d(tau), test v*d ~ const (D5), C2_dyn (D6) ---")
dyn_results = {}
fl_dyn = []
for d0 in D0_DYN:
    psi0 = make_pair(G_MAIN, [(-0.5 * d0, 0.0, +1), (+0.5 * d0, 0.0, -1)])
    _, taus, dm, Hs, fl, _ = evolve(psi0, G_MAIN, 120000, stop_d=3.0)
    fl_dyn.append(fl)
    fin = np.isfinite(dm)
    # znak: monotoniczne zblizanie
    incs = np.diff(dm[fin])
    frac_up = float(np.mean(incs > 0.0)) if len(incs) else np.nan
    # D5: v*d na oknie d in [4, d0-2]
    idx = np.arange(KOFF, len(taus) - KOFF)
    vd_med = vd_lo = vd_hi = np.nan
    if len(idx) > 3:
        v = -(dm[idx + KOFF] - dm[idx - KOFF]) / (taus[idx + KOFF] - taus[idx - KOFF])
        dmid = dm[idx]
        sel = np.isfinite(v) & np.isfinite(dmid) & (dmid >= 4.0) & (dmid <= d0 - 2.0)
        if sel.sum() > 3:
            vd = v[sel] * dmid[sel]
            vd_med = float(np.median(vd))
            vd_lo, vd_hi = map(float, np.percentile(vd, [10, 90]))
    # D6: C2_dyn = dH/dln(d) na oknie d in [5, d0-3]
    selE = fin & (dm >= 5.0) & (dm <= d0 - 3.0)
    C2_dyn = np.nan
    if selE.sum() > 5:
        C2_dyn = fit_L(dm[selE], Hs[selE])["C2"]
    dyn_results[d0] = {"C2_dyn": C2_dyn, "vd": (vd_med, vd_lo, vd_hi),
                       "frac_up": frac_up, "tau_end": taus[-1] if len(taus) else 0,
                       "d_end": dm[fin][-1] if fin.any() else np.nan}
    print(f"  d0={d0:5.1f}: tau_end={dyn_results[d0]['tau_end']:8.1f}  "
          f"d: {d0:.1f} -> {dyn_results[d0]['d_end']:.3f}  "
          f"frakcja krokow z dd>0: {frac_up:.4f}")
    print(f"           v*d (okno d in [4, {d0-2:.0f}]): mediana={vd_med:.4f} "
          f"[p10={vd_lo:.4f}, p90={vd_hi:.4f}]  "
          f"C2_dyn (D6, okno d in [5, {d0-3:.0f}]) = {C2_dyn:.4f}")
print()

# ---------------- D. pinning na siatce glownej (D9) ----------------
print("--- test pinningu (D9): przesuniecie rdzeni o pol wezla ---")


def pin_test(g, d0):
    h = 0.5 * g["dx"]
    E = {}
    for tag, off in (("baza", 0.0), ("shift", h)):
        psi0 = make_pair(g, [(-0.5 * d0 + off, off, +1),
                             (+0.5 * d0 + off, off, -1)])
        _, _, _, _, _, snaps = evolve(psi0, g, TAU_RELAX, snap_at=(TAU_RELAX,))
        E[tag] = snaps[TAU_RELAX]["E"]
    return abs(E["shift"] - E["baza"])


E_range_main = float(E200.max() - E200.min())
dE_pin_main = pin_test(G_MAIN, 12.0)
print(f"  GLOWNA dx=0.5:  Delta_E_pin={dE_pin_main:.3e}  "
      f"rozstep E(d)={E_range_main:.4f}  "
      f"pin/rozstep={100*dE_pin_main/E_range_main:.3f}%  (prog 2%)")

# ---------------- E. kontrola dx=0.25 (D9) ----------------
print()
print("--- kontrola dx=0.25 (D9): L=64, N=256 ---")
G_CTRL = make_grid(256, 64.0)
ns_ctrl, okstar_ctrl, dev_ctrl, fl_cc = core_check(
    G_CTRL, 32.0, "KONTROLA dx=0.25", ode_r, ode_f)
rows_c, fl_statc = static_scan(G_CTRL, D_GRID_CTRL, "KONTROLA dx=0.25")
d200c = np.array([r["d200"] for r in rows_c])
E200c = np.array([r["E200"] for r in rows_c])
wind_ok_ctrl = all(r["ns"] == [-1.0, 1.0] for r in rows_c)
fLc = fit_L(d200c, E200c)
E_range_ctrl = float(E200c.max() - E200c.min())
dE_pin_ctrl = pin_test(G_CTRL, 8.0)
print(f"  C2_ctrl (3 punkty, raport) = {fLc['C2']:.4f}")
print(f"  Delta_E_pin={dE_pin_ctrl:.3e}  rozstep E(d)={E_range_ctrl:.4f}  "
      f"pin/rozstep={100*dE_pin_ctrl/E_range_ctrl:.3f}%  (prog 2%)")
# mobilnosc rdzeni na drobnej siatce (krotki run dynamiczny)
psi0 = make_pair(G_CTRL, [(-4.0, 0.0, +1), (+4.0, 0.0, -1)])
_, tausc, dmc, _, flc2, _ = evolve(psi0, G_CTRL, 4000, stop_d=3.0)
finc = np.isfinite(dmc)
print(f"  mobilnosc (d0=8): d {dmc[finc][0]:.3f} -> {dmc[finc][-1]:.3f} "
      f"w tau={tausc[-1]:.1f} (rdzenie ruchome: {dmc[finc][-1] < dmc[finc][0] - 0.5})")
print()

# ---------------- EWALUACJA ----------------
print("=== EWALUACJA (LOCKED) ===")
all_flags = ([fl_core] + fl_stat + fl_dyn + [fl_cc] + fl_statc + [flc2])
g1_ok = all(not (f["nonfinite"] or f["runaway"]) and f["H_mono"]
            for f in all_flags)
print(f"G1 (czesc Phase2): finite/no-runaway + H_Gamma nierosnace "
      f"na wszystkich runach: {'PASS' if g1_ok else 'FAIL'}")

pin_ok = (dE_pin_main <= 0.02 * E_range_main) and (dE_pin_ctrl <= 0.02 * E_range_ctrl)
mobile_ok = bool(dmc[finc][-1] < dmc[finc][0] - 0.5) and \
    all(np.isfinite(dyn_results[d0]["d_end"]) and
        dyn_results[d0]["d_end"] < d0 - 2.0 for d0 in D0_DYN)
g3 = (ns_core == [-1.0, 1.0]) and wind_ok_stat and wind_ok_ctrl \
    and okstar_main and okstar_ctrl and pin_ok and mobile_ok
print(f"G3 (istnienie i stabilnosc wirow):")
print(f"  kret zachowany po relaksacji (core/statyka/kontrola): "
      f"{ns_core == [-1.0, 1.0]}/{wind_ok_stat}/{wind_ok_ctrl}")
print(f"  |psi| -> s* dla r >> xi (dx=0.5 / dx=0.25): {okstar_main}/{okstar_ctrl}")
print(f"  pinning nieistotny (D9, prog 2%): {pin_ok};  rdzenie ruchome: {mobile_ok}")
print(f"  -> {'PASS' if g3 else 'FAIL'}")

winner = "L" if dAIC >= 10 else ("Y" if dAIC <= -10 else "BRAK")
c2_ok = BAND[0] <= fL["C2"] <= BAND[1]
g4 = (winner == "L") and c2_ok
print(f"G4 (RDZEN, prawo dlugozasiegowe): Delta_AIC(Y-L)={dAIC:+.2f} -> "
      f"wygrywa: {winner}")
print(f"  C2={fL['C2']:.4f} in [{BAND[0]:.3f}, {BAND[1]:.3f}]: {c2_ok} "
      f"(predykcja 2*pi*kappa*Phi* = {C2_PRED:.4f}; "
      f"C2/C2_pred = {fL['C2']/C2_PRED:.4f})")
print(f"  -> {'PASS' if g4 else 'FAIL'}")
print()
print(f"PODSUMOWANIE PHASE 2: G1(czesc)={'PASS' if g1_ok else 'FAIL'} "
      f"G3={'PASS' if g3 else 'FAIL'} G4={'PASS' if g4 else 'FAIL'}"
      f"  [forma: {'logarytmiczna (2D Newton)' if winner == 'L' else winner}]")
