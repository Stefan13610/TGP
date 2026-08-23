# Phase 2: skan ugiecia pakietu na wirze — G3, G4, G6 + obserwacje
# Implementacja SCISLE wg Phase0_balance.md (LOCKED; sekcje 2-4, D1-D10).
#
# Kazdy run: psi_w (tlo+pakiet) i psi_ref (tlo bez pakietu) w LOCKSTEP
# tym samym integratorem (D3); delta_psi = psi_w - psi_ref;
# pozycje rdzeni z psi_ref; b_eff = min odleglosc centroidu |dpsi|^2
# od rdzenia rozpraszajacego (D5); alpha_pred = ray tracing przy b_eff
# na tym samym zrelaksowanym tle (D6).
#
# Runy: kryterialne omega {1.1,1.3} x b {8,12,16};
#       G6: (1.1, 8, u0/2); obserwacje: (1.1, 6), (1.0, 8).

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt_flow = 0.02
dt = 0.05
N, L = 256, 128.0
dx = L / N
TWO_PI = 2.0 * np.pi
S_GUARD = 10.0
H_TOL = 1e-10
SEED = 20260704

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)
Phi_star = s_star ** 2
Vpp_star = a - 2.0 * b_par * s_star + 3.0 * c * s_star ** 2
xi2 = kappa / Vpp_star

coords = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
Z = X + 1j * Y
V1 = (-32.0, -32.0)   # rdzen rozpraszajacy (+1 na tle zwyklym)
V2 = (+32.0, +32.0)


def gfun(m):
    return a - b_par * np.abs(m) + c * m * m


def Vpp(m):
    return a - 2.0 * b_par * np.abs(m) + 3.0 * c * m * m


def Vpot(m):
    return 0.5 * a * m**2 - (b_par / 3.0) * m**3 + 0.25 * c * m**4


def lap(u):
    return (np.roll(u, 1, 0) + np.roll(u, -1, 0)
            + np.roll(u, 1, 1) + np.roll(u, -1, 1) - 4.0 * u) / (dx * dx)


def rhs(psi):
    return kappa * lap(psi) - psi * gfun(np.abs(psi))


def energy(psi, psit):
    gx = (np.roll(psi, -1, 0) - psi) / dx
    gy = (np.roll(psi, -1, 1) - psi) / dx
    return float(np.sum(0.5 * np.abs(psit) ** 2
                        + 0.5 * kappa * (np.abs(gx) ** 2 + np.abs(gy) ** 2)
                        + Vpot(np.abs(psi))) * dx * dx)


def wrap(x):
    return (x + np.pi) % TWO_PI - np.pi


def theta_h_pair(z1, z2):
    th = np.zeros_like(X)
    for (x0, y0), n in ((z1, +1), (z2, -1)):
        for k in range(-2, 3):
            th = th + n * np.angle(np.sinh(np.pi * (Z - (x0 + k * L) - 1j * y0) / L))
    return th


def theta_v_pair(z1, z2):
    Wt = Y + 1j * X
    th = np.zeros_like(X)
    for (x0, y0), n in ((z1, +1), (z2, -1)):
        for k in range(-2, 3):
            th = th - n * np.angle(np.sinh(np.pi * (Wt - (y0 + k * L) - 1j * x0) / L))
    return th


def make_pair(mirror=False):
    A_, B_ = V1, V2
    C_ = (B_[0], A_[1])
    th = theta_h_pair(A_, C_) + theta_v_pair(C_, B_)
    if mirror:
        th = -th
    Amp = np.full_like(X, s_star)
    for (x0, y0) in (V1, V2):
        ddx = (X - x0 + 0.5 * L) % L - 0.5 * L
        ddy = (Y - y0 + 0.5 * L) % L - 0.5 * L
        rr = np.sqrt(ddx * ddx + ddy * ddy)
        Amp = Amp * rr / np.sqrt(rr * rr + 2.0 * xi2)
    return Amp * np.exp(1j * th)


def core_positions(psi):
    th = np.angle(psi)
    dthx = wrap(np.roll(th, -1, 0) - th)
    dthy = wrap(np.roll(th, -1, 1) - th)
    W = (dthx + np.roll(dthy, -1, 0) - np.roll(dthx, -1, 1) - dthy) / TWO_PI
    Phi = np.abs(psi) ** 2
    out = {}
    for sign in (+1, -1):
        cells = np.argwhere(W * sign > 0.5)
        if len(cells) == 0:
            return None
        ic, jc = cells[0]
        off = np.arange(-6, 7)
        ii, jj = (ic + off) % N, (jc + off) % N
        Pw = Phi[np.ix_(ii, jj)]
        wgt = np.clip(0.5 * Phi_star - Pw, 0.0, None)
        OX, OY = np.meshgrid(off * dx, off * dx, indexing="ij")
        sx = float(np.sum(wgt * OX) / wgt.sum())
        sy = float(np.sum(wgt * OY) / wgt.sum())
        out[sign] = (coords[ic] + 0.5 * dx + sx, coords[jc] + 0.5 * dx + sy)
    ntot = float(np.rint(W).sum())
    return out, ntot


def pdist(p, q):
    ddx = (p[0] - q[0] + 0.5 * L) % L - 0.5 * L
    ddy = (p[1] - q[1] + 0.5 * L) % L - 0.5 * L
    return float(np.hypot(ddx, ddy))


# ---------------- eikonal na rzeczywistym tle (D6) ----------------
def build_n2(psi_bg, om):
    f_bg = np.abs(psi_bg)
    n2 = (om * om - Vpp(f_bg)) / (om * om - Vpp_star)
    gx = (np.roll(n2, -1, 0) - np.roll(n2, 1, 0)) / (2 * dx)
    gy = (np.roll(n2, -1, 1) - np.roll(n2, 1, 1)) / (2 * dx)
    return n2, gx, gy


def bilin(G, x, y):
    fx = ((x - coords[0]) / dx) % N
    fy = ((y - coords[0]) / dx) % N
    i0, j0 = int(fx) % N, int(fy) % N
    i1, j1 = (i0 + 1) % N, (j0 + 1) % N
    tx, ty = fx - int(fx), fy - int(fy)
    return ((1 - tx) * (1 - ty) * G[i0, j0] + tx * (1 - ty) * G[i1, j0]
            + (1 - tx) * ty * G[i0, j1] + tx * ty * G[i1, j1])


def trace_bg(n2g, b_signed, core_y=V1[1]):
    n2, gxg, gyg = n2g
    x, y = -60.0, core_y + b_signed
    px, py = np.sqrt(bilin(n2, x, y)), 0.0
    ds = 0.02

    def F(st):
        return np.array([st[2], st[3],
                         0.5 * bilin(gxg, st[0], st[1]),
                         0.5 * bilin(gyg, st[0], st[1])])

    st = np.array([x, y, px, py])
    for _ in range(500000):
        k1 = F(st)
        k2 = F(st + 0.5 * ds * k1)
        k3 = F(st + 0.5 * ds * k2)
        k4 = F(st + ds * k3)
        st = st + ds / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        if st[0] > 40.0:
            sgn = -1.0 if b_signed > 0 else +1.0
            return sgn * np.degrees(np.arctan2(st[3], st[2]))
        if np.hypot(st[0] - V1[0], st[1] - V1[1]) < 0.6:
            return np.nan
    return np.nan


# ---------------- run ugiecia (D1-D5, D8, D9) ----------------
def deflection_run(psi_bg, n2cache, om, b_nom, u0=1e-3, scatter_sign=+1,
                   tag=""):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    max_steps = 8000 if om < 1.05 else 6000
    th_bg = np.angle(psi_bg)
    eith = np.exp(1j * th_bg)
    x0, y0 = -60.0, V1[1] + b_nom
    env = u0 * np.exp(-((X - x0) ** 2) / (2 * 8.0 ** 2)
                      - ((Y - y0) ** 2) / (2 * 8.0 ** 2))
    psi_w = psi_bg + eith * env * np.cos(k0 * (X - x0))
    psi_wp = psi_bg + eith * env * np.cos(k0 * (X - x0) + om * dt)  # D2
    psi_r = psi_bg.copy()
    psi_rp = psi_bg.copy()

    flags = {"nonfinite": False, "runaway": False, "wrap": False,
             "no5": False, "core_lost": False}
    Es, taus_E = [], []
    core_series = []          # (tau, pos_scatter, pos_other)
    lncore_series = []        # (tau, ln RMS dpsi w oknach rdzeni)
    b_eff = np.inf
    triggered, x_thr, tau_trig = False, None, None
    samples, next_s = [], None
    n_tot_end = None

    for stp in range(1, max_steps + 1):
        tau = stp * dt
        psi_wn = 2.0 * psi_w - psi_wp + dt * dt * rhs(psi_w)
        psi_rn = 2.0 * psi_r - psi_rp + dt * dt * rhs(psi_r)
        if stp % 10 == 0:
            if not (np.all(np.isfinite(psi_wn)) and np.all(np.isfinite(psi_rn))):
                flags["nonfinite"] = True
                break
            if np.max(np.abs(psi_wn)) > S_GUARD:
                flags["runaway"] = True
                break
            dpsi = psi_w - psi_r
            w2 = np.abs(dpsi) ** 2
            Wtot = float(w2.sum())
            xc = float((w2 * X).sum() / Wtot)
            yc = float((w2 * Y).sum() / Wtot)
            if np.hypot(xc - x0, yc - y0) > 0.9 * L:
                flags["wrap"] = True
                break
            cp = core_positions(psi_r)
            if cp is None:
                flags["core_lost"] = True
                break
            pos_sc = cp[0][scatter_sign]
            core_series.append((tau, pos_sc, cp[0][-scatter_sign]))
            n_tot_end = cp[1]
            # gamma_core: RMS dpsi w oknach r<5 wokol obu rdzeni
            m_core = np.zeros((N, N), bool)
            for s_ in (+1, -1):
                px_, py_ = cp[0][s_]
                rr = np.hypot((X - px_ + L / 2) % L - L / 2,
                              (Y - py_ + L / 2) % L - L / 2)
                m_core |= rr < 5.0
            lncore_series.append(
                (tau, float(np.log(np.sqrt(np.mean(w2[m_core])) + 1e-300))))
            # energia (dryft sekularny, D8)
            ut_w = (psi_wn - psi_wp) / (2 * dt)
            Es.append(energy(psi_w, ut_w))
            taus_E.append(tau)
            if not triggered:
                b_eff = min(b_eff, pdist((xc, yc), pos_sc))
                thr = pos_sc[0] + 10.0
                m_tr = X > thr
                if float(w2[m_tr].sum()) > 0.5 * Wtot:
                    triggered, x_thr, tau_trig = True, thr, tau
                    next_s = tau
            if triggered and next_s is not None and tau >= next_s \
                    and len(samples) < 5:
                m_tr = X > x_thr
                ut_d = ((psi_wn - psi_rn) - (psi_wp - psi_rp)) / (2 * dt)
                gx = (np.roll(dpsi, -1, 0) - np.roll(dpsi, 1, 0)) / (2 * dx)
                gy = (np.roll(dpsi, -1, 1) - np.roll(dpsi, 1, 1)) / (2 * dx)
                Px = float(-np.real(np.sum(np.conj(ut_d)[m_tr] * gx[m_tr])) * dx * dx)
                Py = float(-np.real(np.sum(np.conj(ut_d)[m_tr] * gy[m_tr])) * dx * dx)
                samples.append((tau, Px, Py))
                next_s += 4.0
            if len(samples) >= 5:
                psi_wp, psi_w = psi_w, psi_wn
                psi_rp, psi_r = psi_r, psi_rn
                break
        psi_wp, psi_w = psi_w, psi_wn
        psi_rp, psi_r = psi_r, psi_rn

    out = {"om": om, "b_nom": b_nom, "u0": u0, "tag": tag, "flags": flags,
           "b_eff": b_eff if np.isfinite(b_eff) else np.nan,
           "tau_trig": tau_trig, "n_tot_end": n_tot_end}
    sgn = -1.0 if b_nom > 0 else +1.0
    if len(samples) == 5:
        Pxm = float(np.mean([s[1] for s in samples]))
        Pym = float(np.mean([s[2] for s in samples]))
        out["alpha"] = sgn * float(np.degrees(np.arctan2(Pym, Pxm)))
        al_s = [sgn * float(np.degrees(np.arctan2(s[2], s[1])))
                for s in samples]
        out["alpha_samples"] = al_s
        out["sigma"] = float(np.std(al_s))
    else:
        flags["no5"] = True
        out["alpha"], out["sigma"], out["alpha_samples"] = np.nan, np.nan, []
    # dryft sekularny energii
    q = max(len(Es) // 4, 1)
    out["E_drift"] = abs(np.mean(Es[-q:]) - np.mean(Es[:q])) / abs(Es[0])
    # gamma_core: fit po tranzycie (od wyzwolenia; D8) + druga polowa
    tt = np.array([p[0] for p in lncore_series])
    yy = np.array([p[1] for p in lncore_series])

    def slope(mask):
        if mask.sum() < 5:
            return np.nan
        A_ = np.vstack([np.ones(mask.sum()), tt[mask]]).T
        coef, *_ = np.linalg.lstsq(A_, yy[mask], rcond=None)
        return float(coef[1])

    out["gamma_core_post"] = slope(tt >= (tau_trig if tau_trig else tt[-1]))
    out["gamma_core_2half"] = slope(tt >= 0.5 * tt[-1])
    # dryf rdzenia w oknie runu
    out["core_start"] = core_series[0][1]
    out["core_end"] = core_series[-1][1]
    out["core_disp"] = pdist(core_series[0][1], core_series[-1][1])
    # alpha_pred przy b_eff (D5): ray tracing na tym samym tle
    if np.isfinite(out["b_eff"]):
        out["alpha_pred"] = trace_bg(n2cache, np.sign(b_nom) * out["b_eff"])
    else:
        out["alpha_pred"] = np.nan
    return out


def report(r):
    fl = r["flags"]
    print(f"  run [{r['tag']}] omega={r['om']:.1f} b_nom={r['b_nom']:+.0f} "
          f"u0={r['u0']:.1e}:")
    cs, ce = r['core_start'], r['core_end']
    print(f"    b_eff = {r['b_eff']:.3f}  tau_wyzw = {r['tau_trig']}  "
          f"rdzen: ({cs[0]:.2f},{cs[1]:.2f}) -> ({ce[0]:.2f},{ce[1]:.2f}) "
          f"(dryf {r['core_disp']:.2f})")
    if r["alpha_samples"]:
        print(f"    probki alpha [deg]: "
              f"{[f'{x:+.2f}' for x in r['alpha_samples']]}")
    print(f"    alpha = {r['alpha']:+.3f} +/- {r['sigma']:.3f} deg   "
          f"alpha_pred(b_eff) = {r['alpha_pred']:+.3f} deg   "
          f"alpha/alpha_pred = "
          f"{r['alpha']/r['alpha_pred'] if np.isfinite(r['alpha_pred']) and r['alpha_pred'] != 0 else np.nan:.3f}")
    print(f"    G1: E_dryft = {r['E_drift']:.2e}; gamma_core(po tranzycie) = "
          f"{r['gamma_core_post']:+.5f}, (2. polowa) = "
          f"{r['gamma_core_2half']:+.5f}; n_tot = {r['n_tot_end']:.0f}; "
          f"flagi: {[k for k, v in fl.items() if v] or 'brak'}")


print("=== Phase 2: skan ugiecia (LOCK: Phase0_balance.md, D1-D10) ===")
print(f"N={N} L={L} dx={dx} dt={dt}; tlo: para (+1,-1); "
      f"pakiet u0=1e-3 sx=sy=8, start x0=-60, y0=-32+b")
print()

# ---------------- tlo ----------------
psi = make_pair()
for stp in range(1, 2001):
    psi = psi + dt_flow * rhs(psi)
psi_bg = psi
print(f"tlo zrelaksowane (2000 krokow dt={dt_flow}); "
      f"residuum max|RHS| = {float(np.max(np.abs(rhs(psi_bg)))):.3e}")
n2c = {om: build_n2(psi_bg, om) for om in (1.0, 1.1, 1.3)}
print()

# ---------------- skan glowny (G3, G4) ----------------
print("--- skan glowny: omega {1.1, 1.3} x b {8, 12, 16} (KRYTERIALNE) ---")
runs = {}
for om in (1.1, 1.3):
    for bb in (8.0, 12.0, 16.0):
        r = deflection_run(psi_bg, n2c[om], om, bb, tag="kryt")
        runs[(om, bb)] = r
        report(r)
print()

# ---------------- G6: liniowosc ----------------
print("--- G6: liniowosc, (omega=1.1, b=8), u0 -> u0/2 ---")
r_half = deflection_run(psi_bg, n2c[1.1], 1.1, 8.0, u0=5e-4, tag="G6")
report(r_half)
print()

# ---------------- obserwacje poza kryteriami ----------------
print("--- obserwacje POZA kryteriami: b=6 (silne pole), omega=1.0 ---")
r_b6 = deflection_run(psi_bg, n2c[1.1], 1.1, 6.0, tag="obs-b6")
report(r_b6)
r_o10 = deflection_run(psi_bg, n2c[1.0], 1.0, 8.0, tag="obs-om1.0")
report(r_o10)
print()

# ---------------- EWALUACJA ----------------
print("=== EWALUACJA (LOCKED) ===")
crit = [runs[(om, bb)] for om in (1.1, 1.3) for bb in (8.0, 12.0, 16.0)]
g1_flags = all(not any(r["flags"].values()) for r in crit + [r_half])
g1_E = all(r["E_drift"] <= 1e-4 for r in crit + [r_half])
g1_gc = all(r["gamma_core_post"] <= 0.01 for r in crit + [r_half])
g1_n = all(r["n_tot_end"] == 0 for r in crit + [r_half])
print(f"G1 (czesc Phase2): flagi czyste: {g1_flags}; E_dryft <= 1e-4: "
      f"{g1_E}; gamma_core <= 0.01: {g1_gc}; kret zachowany: {g1_n}  "
      f"-> {'PASS' if all([g1_flags, g1_E, g1_gc, g1_n]) else 'FAIL'}")

r0 = runs[(1.1, 8.0)]
g3_a = np.isfinite(r0["alpha"]) and abs(r0["alpha"]) > 3 * r0["sigma"]
g3_b = (abs(runs[(1.1, 8.0)]["alpha"]) > abs(runs[(1.1, 12.0)]["alpha"])
        > abs(runs[(1.1, 16.0)]["alpha"]))
g3_c = abs(runs[(1.1, 8.0)]["alpha"]) > abs(runs[(1.3, 8.0)]["alpha"])
g3 = g3_a and g3_b and g3_c
print(f"G3 (istnienie): |alpha(8,1.1)|={abs(r0['alpha']):.2f} > "
      f"3*sigma={3*r0['sigma']:.2f}: {g3_a}; maleje z b: {g3_b}; "
      f"maleje z omega: {g3_c}  -> {'PASS' if g3 else 'FAIL'}")

g4_ok = True
print("G4 (RDZEN, policzalnosc): znak KU (alpha>0) ORAZ "
      "alpha/alpha_pred(b_eff) in [0.5, 2.0]:")
for om in (1.1, 1.3):
    for bb in (8.0, 12.0, 16.0):
        r = runs[(om, bb)]
        ratio = (r["alpha"] / r["alpha_pred"]
                 if np.isfinite(r["alpha_pred"]) and r["alpha_pred"] != 0
                 else np.nan)
        ok = (np.isfinite(r["alpha"]) and r["alpha"] > 0
              and np.isfinite(ratio) and 0.5 <= ratio <= 2.0)
        g4_ok = g4_ok and ok
        print(f"    omega={om:.1f} b={bb:2.0f}: alpha={r['alpha']:+7.3f} "
              f"pred(b_eff={r['b_eff']:.2f})={r['alpha_pred']:+7.3f} "
              f"ratio={ratio:5.3f}  -> {'OK' if ok else 'POZA'}")
print(f"  -> G4: {'PASS' if g4_ok else 'FAIL'}")

d6 = abs(r0["alpha"] - r_half["alpha"])
thr6 = 0.1 * max(abs(r0["alpha"]), abs(r_half["alpha"]))
g6 = d6 <= thr6
print(f"G6 (liniowosc): |{r0['alpha']:+.3f} - {r_half['alpha']:+.3f}| = "
      f"{d6:.3f} <= {thr6:.3f}: {g6}  -> {'PASS' if g6 else 'FAIL'}")
print()
print(f"PODSUMOWANIE PHASE 2: "
      f"G1(czesc)={'PASS' if all([g1_flags, g1_E, g1_gc, g1_n]) else 'FAIL'} "
      f"G3={'PASS' if g3 else 'FAIL'} G4={'PASS' if g4_ok else 'FAIL'} "
      f"G6={'PASS' if g6 else 'FAIL'}")
