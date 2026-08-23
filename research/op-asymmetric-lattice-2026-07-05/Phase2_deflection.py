# Phase 2: skan ugiecia na siatce SKOSNEJ — G4, G5a, G5b (RDZEN), G6b
#          + obserwacje poza kryteriami: (1.1, +28), (1.1, -28)
# Implementacja SCISLE wg Phase0_balance.md proby #4 (LOCKED; sekcje
# 2-4, D1-D6, D8-D10, D12-D14, D18).
# WARUNEK WSTEPNY (protokol, D19): bramka G2 PASS w Phase 1.
#
# Skan kryterialny (ZAMKNIETY, 5 par +/-b):
#   omega=1.1 x b in {+14,-14,+20,-20,+24,-24}
#   omega=1.5 x b in {+8,-8,+12,-12}
# G5a: WSZYSTKIE 10 punktow: sign(alpha)=sign(alpha_pred) ORAZ
#      alpha/alpha_pred(b_eff) in [0.5, 2.0]
# G5b (D18): per para: Delta_meas = alpha(+b)-alpha(-b);
#      Delta_eik = alpha_pred(b_eff+)-alpha_pred(b_eff-) (kazda przy
#      WLASNYM b_eff i wlasnej plaszczyznie); sigma_Delta =
#      sqrt(sigma+^2+sigma-^2); |Delta_meas - Delta_eik| <= 3*sigma_Delta
#      + regula wzmacniajaca (sekcja 4). Delta_meas ZAWSZE obok
#      Delta_eik (forbidden move #5).
# alpha > 0 = KU rdzeniowi A dla OBU znakow b (D4: -sign(b)*atan2).

import numpy as np
from itertools import permutations

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt_flow = 0.02
dt = 0.05
N, L = 512, 256.0
dx = L / N
TWO_PI = 2.0 * np.pi
S_GUARD = 10.0
SEED = 20260704

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)
Phi_star = s_star ** 2
Vpp_star = a - 2.0 * b_par * s_star + 3.0 * c * s_star ** 2
xi2 = kappa / Vpp_star

coords = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
Z = X + 1j * Y

VORTS = [((-64.0, -64.0), +1), ((+64.0, -64.0), -1),
         ((-128.0, 0.0), +1), ((0.0, 0.0), -1),
         ((+64.0, +64.0), +1), ((-64.0, +64.0), -1),
         ((0.0, -128.0), +1), ((+128.0, -128.0), -1)]
NOM = {+1: [(-64.0, -64.0), (-128.0, 0.0), (+64.0, +64.0), (0.0, -128.0)],
       -1: [(+64.0, -64.0), (0.0, 0.0), (-64.0, +64.0), (+128.0, -128.0)]}
Y_ROW = -64.0
X0 = -100.0


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


def pdist(p, q):
    ddx = (p[0] - q[0] + 0.5 * L) % L - 0.5 * L
    ddy = (p[1] - q[1] + 0.5 * L) % L - 0.5 * L
    return float(np.hypot(ddx, ddy))


def theta_h_pair(z1, z2):
    th = np.zeros_like(X)
    for (x0, y0), n in ((z1, +1), (z2, -1)):
        for k in range(-2, 3):
            th = th + n * np.angle(np.sinh(np.pi * (Z - (x0 + k * L) - 1j * y0) / L))
    return th


def make_lattice():
    th = theta_h_pair((-64.0, -64.0), (+64.0, -64.0)) \
        + theta_h_pair((-128.0, 0.0), (0.0, 0.0)) \
        + theta_h_pair((+64.0, +64.0), (-64.0, +64.0)) \
        + theta_h_pair((0.0, -128.0), (+128.0, -128.0))
    Amp = np.full_like(X, s_star)
    for (x0, y0), _n in VORTS:
        ddx = (X - x0 + 0.5 * L) % L - 0.5 * L
        ddy = (Y - y0 + 0.5 * L) % L - 0.5 * L
        rr = np.sqrt(ddx * ddx + ddy * ddy)
        Amp = Amp * rr / np.sqrt(rr * rr + 2.0 * xi2)
    return Amp * np.exp(1j * th)


def core_positions(psi, nom=NOM):
    """Detektor 4+4 (D14)."""
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
        clusters = []
        for (i, j) in cells:
            placed = False
            for cl in clusters:
                i0, j0 = cl[0]
                di = (i - i0 + N // 2) % N - N // 2
                dj = (j - j0 + N // 2) % N - N // 2
                if abs(di) <= 8 and abs(dj) <= 8:
                    cl.append((i, j))
                    placed = True
                    break
            if not placed:
                clusters.append([(i, j)])
        if len(clusters) != 4:
            return None
        pos = []
        for cl in clusters:
            ic, jc = cl[0]
            off = np.arange(-6, 7)
            ii, jj = (ic + off) % N, (jc + off) % N
            Pw = Phi[np.ix_(ii, jj)]
            wgt = np.clip(0.5 * Phi_star - Pw, 0.0, None)
            OX, OY = np.meshgrid(off * dx, off * dx, indexing="ij")
            sx = float(np.sum(wgt * OX) / wgt.sum())
            sy = float(np.sum(wgt * OY) / wgt.sum())
            pos.append((coords[ic] + 0.5 * dx + sx, coords[jc] + 0.5 * dx + sy))
        best = None
        for perm in permutations(range(4)):
            ssum = sum(pdist(pos[perm[t]], nom[sign][t]) for t in range(4))
            if best is None or ssum < best[0]:
                best = (ssum, perm)
        for t in range(4):
            out[(sign, t)] = pos[best[1][t]]
    ntot = float(np.rint(W).sum())
    return out, ntot


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


def trace_bg(n2g, b_signed, planes):
    """Ray na tle siatki skosnej; kat [deg, KU rdzeniowi A] w
    plaszczyznach x = planes (dict); None = przechwyt (r<0.6)."""
    n2, gxg, gyg = n2g
    x, y = X0, Y_ROW + b_signed
    px, py = np.sqrt(bilin(n2, x, y)), 0.0
    ds = 0.02
    sgn = -1.0 if b_signed > 0 else +1.0

    def F(st):
        return np.array([st[2], st[3],
                         0.5 * bilin(gxg, st[0], st[1]),
                         0.5 * bilin(gyg, st[0], st[1])])

    st = np.array([x, y, px, py])
    out = {}
    todo = sorted(planes)
    for _ in range(1000000):
        k1 = F(st)
        k2 = F(st + 0.5 * ds * k1)
        k3 = F(st + 0.5 * ds * k2)
        k4 = F(st + ds * k3)
        st_new = st + ds / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        while todo and st.item(0) <= todo[0] < st_new.item(0):
            f_ = (todo[0] - st[0]) / (st_new[0] - st[0])
            pxi = st[2] + f_ * (st_new[2] - st[2])
            pyi = st[3] + f_ * (st_new[3] - st[3])
            out[todo.pop(0)] = sgn * float(np.degrees(np.arctan2(pyi, pxi)))
        st = st_new
        if not todo:
            return out
        for (x0, y0), _n in VORTS:
            if np.hypot((st[0] - x0 + L/2) % L - L/2,
                        (st[1] - y0 + L/2) % L - L/2) < 0.6:
                return None
    return None


# ---------------- run ugiecia (D1-D6, D8, D9) ----------------
def deflection_run(psi_bg, n2cache, om, b_nom, u0=1e-3, tag=""):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    lam = TWO_PI / k0
    max_steps = 6000                              # D9
    th_bg = np.angle(psi_bg)
    eith = np.exp(1j * th_bg)
    x0, y0 = X0, Y_ROW + b_nom
    env = u0 * np.exp(-((X - x0) ** 2) / (2 * 8.0 ** 2)
                      - ((Y - y0) ** 2) / (2 * 8.0 ** 2))
    psi_w = psi_bg + eith * env * np.cos(k0 * (X - x0))
    psi_wp = psi_bg + eith * env * np.cos(k0 * (X - x0) + om * dt)  # D2
    psi_r = psi_bg.copy()
    psi_rp = psi_bg.copy()

    flags = {"nonfinite": False, "runaway": False, "wrap": False,
             "no5": False, "core_lost": False}
    Es = []
    lncore_series = []
    core_start, core_end = None, None
    b_eff = np.inf
    triggered, tau_trig = False, None
    x_thr1, x_thr2, m_win = None, None, None
    samples, next_s = [], None
    xc3 = None
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
            pos_sc = cp[0][(+1, 0)]           # rdzen A
            pos_x2 = cp[0][(-1, 0)]           # drugi wir rzedu y=-64
            if core_start is None:
                core_start = pos_sc
            core_end = pos_sc
            n_tot_end = cp[1]
            # gamma_core: RMS dpsi w oknach r<5 wokol WSZYSTKICH 8 rdzeni
            m_core = np.zeros((N, N), bool)
            for key in cp[0]:
                px_, py_ = cp[0][key]
                rr = np.hypot((X - px_ + L / 2) % L - L / 2,
                              (Y - py_ + L / 2) % L - L / 2)
                m_core |= rr < 5.0
            lncore_series.append(
                (tau, float(np.log(np.sqrt(np.mean(w2[m_core])) + 1e-300))))
            ut_w = (psi_wn - psi_wp) / (2 * dt)
            Es.append(energy(psi_w, ut_w))
            if not triggered:
                b_eff = min(b_eff, pdist((xc, yc), pos_sc))
                if xc > pos_sc[0] + 10.0:     # wyzwolenie CENTROIDEM (D4)
                    triggered, tau_trig = True, tau
                    x_thr1 = pos_sc[0] + 10.0             # ZAMROZONE
                    x_thr2 = pos_x2[0] - 10.0             # ZAMROZONE
                    m_win = (X > x_thr1) & (X < x_thr2)
                    next_s = tau
            if triggered and next_s is not None and tau >= next_s \
                    and len(samples) < 5:
                ut_d = ((psi_wn - psi_rn) - (psi_wp - psi_rp)) / (2 * dt)
                gx = (np.roll(dpsi, -1, 0) - np.roll(dpsi, 1, 0)) / (2 * dx)
                gy = (np.roll(dpsi, -1, 1) - np.roll(dpsi, 1, 1)) / (2 * dx)
                Px = float(-np.real(np.sum(np.conj(ut_d)[m_win] * gx[m_win])) * dx * dx)
                Py = float(-np.real(np.sum(np.conj(ut_d)[m_win] * gy[m_win])) * dx * dx)
                samples.append((tau, Px, Py))
                if len(samples) == 3:
                    xc3 = xc                  # plaszczyzna ewaluacji D6
                next_s += 4.0
            if len(samples) >= 5:
                psi_wp, psi_w = psi_w, psi_wn
                psi_rp, psi_r = psi_r, psi_rn
                break
        psi_wp, psi_w = psi_w, psi_wn
        psi_rp, psi_r = psi_r, psi_rn

    out = {"om": om, "b_nom": b_nom, "u0": u0, "tag": tag, "flags": flags,
           "b_eff": b_eff if np.isfinite(b_eff) else np.nan,
           "lam": lam, "tau_trig": tau_trig, "n_tot_end": n_tot_end,
           "x_thr1": x_thr1, "x_thr2": x_thr2, "xc3": xc3}
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
    q = max(len(Es) // 4, 1)
    out["E_drift"] = abs(np.mean(Es[-q:]) - np.mean(Es[:q])) / abs(Es[0])
    tt = np.array([p[0] for p in lncore_series])
    yy = np.array([p[1] for p in lncore_series])

    def slope(mask):
        if mask.sum() < 5:
            return np.nan
        A_ = np.vstack([np.ones(mask.sum()), tt[mask]]).T
        coef, *_ = np.linalg.lstsq(A_, yy[mask], rcond=None)
        return float(coef[1])

    out["gamma_core_post"] = slope(tt >= (tau_trig if tau_trig else tt[-1]))
    out["core_start"], out["core_end"] = core_start, core_end
    out["core_disp"] = (pdist(core_start, core_end)
                        if core_start is not None else np.nan)
    # alpha_pred przy WLASNYM b_eff (D5/D6/D18)
    if np.isfinite(out["b_eff"]) and xc3 is not None:
        planes = sorted({round(xc3, 6), round(x_thr1, 6), round(x_thr2, 6)})
        tr = trace_bg(n2cache, np.sign(b_nom) * out["b_eff"], planes)
        if tr is None:
            out["alpha_pred"] = np.nan
            out["pred_lo"], out["pred_hi"] = np.nan, np.nan
        else:
            out["alpha_pred"] = tr[round(xc3, 6)]
            out["pred_lo"] = tr[round(x_thr1, 6)]
            out["pred_hi"] = tr[round(x_thr2, 6)]
    else:
        out["alpha_pred"] = np.nan
        out["pred_lo"], out["pred_hi"] = np.nan, np.nan
    return out


def report(r):
    fl = r["flags"]
    print(f"  run [{r['tag']}] omega={r['om']:.1f} b_nom={r['b_nom']:+.0f} "
          f"u0={r['u0']:.1e}:")
    if r["core_start"] is not None:
        cs, ce = r["core_start"], r["core_end"]
        print(f"    b_eff = {r['b_eff']:.3f}  b_eff/lambda = "
              f"{r['b_eff'] / r['lam']:.3f}  tau_wyzw = {r['tau_trig']}  "
              f"okno = ({r['x_thr1']:.3f}, {r['x_thr2']:.3f})")
        print(f"    rdzen A: ({cs[0]:.2f},{cs[1]:.2f}) -> "
              f"({ce[0]:.2f},{ce[1]:.2f}) (dryf {r['core_disp']:.3f})")
    if r["alpha_samples"]:
        print(f"    probki alpha [deg]: "
              f"{[f'{x:+.2f}' for x in r['alpha_samples']]}")
    ratio = (r["alpha"] / r["alpha_pred"]
             if np.isfinite(r.get("alpha_pred", np.nan)) and r["alpha_pred"] != 0
             else np.nan)
    print(f"    alpha = {r['alpha']:+.3f} +/- {r['sigma']:.3f} deg   "
          f"alpha_pred(b_eff; x_eval="
          f"{r['xc3'] if r['xc3'] is None else round(r['xc3'], 2)}) = "
          f"{r['alpha_pred']:+.3f} deg   ratio = {ratio:.3f}")
    print(f"    czulosc plaszczyzny (krance okna): "
          f"[{r['pred_lo']:+.3f}, {r['pred_hi']:+.3f}] deg")
    print(f"    G1: E_dryft = {r['E_drift']:.2e}; gamma_core(po tranzycie) = "
          f"{r['gamma_core_post']:+.5f}; n_tot = "
          f"{r['n_tot_end'] if r['n_tot_end'] is not None else 'NA'}; "
          f"flagi: {[k for k, v in fl.items() if v] or 'brak'}")


print("=== Phase 2: skan ugiecia na siatce skosnej L=256 (LOCK proby #4: "
      "D1-D6/D8-D10/D12-D14/D18) ===")
print(f"N={N} L={L} dx={dx} dt={dt}; tlo: siatka skosna 8 wirow; "
      f"pakiet u0=1e-3 sx=sy=8, start x0={X0}, y0=-64+b (b ZE ZNAKIEM)")
print("WARUNEK WSTEPNY: bramka G2 PASS (Phase1_output.txt)")
print()

# ---------------- tlo ----------------
psi = make_lattice()
for stp in range(1, 2001):
    psi = psi + dt_flow * rhs(psi)
psi_bg = psi
print(f"tlo zrelaksowane (2000 krokow dt={dt_flow}); "
      f"residuum max|RHS| = {float(np.max(np.abs(rhs(psi_bg)))):.3e}")
n2c = {om: build_n2(psi_bg, om) for om in (1.1, 1.5)}
print()

# ---------------- skan glowny (G4, G5a, G5b) ----------------
print("--- skan glowny: 5 par +/-b (KRYTERIALNE, zestaw ZAMKNIETY) ---")
PAIRS = [(1.1, 14.0), (1.1, 20.0), (1.1, 24.0), (1.5, 8.0), (1.5, 12.0)]
CRIT = []
for om, bb in PAIRS:
    CRIT.append((om, +bb))
    CRIT.append((om, -bb))
runs = {}
for om, bb in CRIT:
    r = deflection_run(psi_bg, n2c[om], om, bb, tag="kryt")
    runs[(om, bb)] = r
    report(r)
print()

# ---------------- G6b: liniowosc ----------------
print("--- G6b: liniowosc, (omega=1.1, b=+14), u0 -> u0/2 ---")
r_half = deflection_run(psi_bg, n2c[1.1], 1.1, +14.0, u0=5e-4, tag="G6b")
report(r_half)
print()

# ---------------- obserwacje poza kryteriami ----------------
print("--- obserwacje POZA kryteriami (D12; niedolaczalne do G5) ---")
r_op = deflection_run(psi_bg, n2c[1.1], 1.1, +28.0, tag="obs+28")
report(r_op)
r_om = deflection_run(psi_bg, n2c[1.1], 1.1, -28.0, tag="obs-28")
report(r_om)
print()

# ---------------- EWALUACJA ----------------
print("=== EWALUACJA (LOCKED) ===")
crit = [runs[k] for k in CRIT]
g1_flags = all(not any(r["flags"].values()) for r in crit + [r_half])
g1_E = all(r["E_drift"] <= 1e-4 for r in crit + [r_half])
g1_gc = all(r["gamma_core_post"] <= 0.01 for r in crit + [r_half])
g1_n = all(r["n_tot_end"] == 0 for r in crit + [r_half])
g1 = all([g1_flags, g1_E, g1_gc, g1_n])
print(f"G1 (czesc Phase2): flagi czyste: {g1_flags}; E_dryft <= 1e-4: "
      f"{g1_E}; gamma_core <= 0.01: {g1_gc}; kret zachowany (n_tot=0, "
      f"8 rdzeni): {g1_n}  -> {'PASS' if g1 else 'FAIL'}")

r0 = runs[(1.1, +14.0)]
g4_a = np.isfinite(r0["alpha"]) and abs(r0["alpha"]) > 3 * r0["sigma"]
g4_bp = (abs(runs[(1.1, +14.0)]["alpha"]) > abs(runs[(1.1, +20.0)]["alpha"])
         > abs(runs[(1.1, +24.0)]["alpha"]))
g4_bm = (abs(runs[(1.1, -14.0)]["alpha"]) > abs(runs[(1.1, -20.0)]["alpha"])
         > abs(runs[(1.1, -24.0)]["alpha"]))
g4 = g4_a and g4_bp and g4_bm
print(f"G4 (istnienie): |alpha(1.1,+14)|={abs(r0['alpha']):.2f} > "
      f"3*sigma={3*r0['sigma']:.2f}: {g4_a}; maleje z |b| strona +: "
      f"{g4_bp}; strona -: {g4_bm}  -> {'PASS' if g4 else 'FAIL'}")

g5a_ok = True
print("G5a (policzalnosc bez oslony symetrii, P2-a): sign(alpha) = "
      "sign(alpha_pred) ORAZ ratio in [0.5, 2.0] dla WSZYSTKICH "
      "10 punktow:")
for om, bb in CRIT:
    r = runs[(om, bb)]
    ratio = (r["alpha"] / r["alpha_pred"]
             if np.isfinite(r["alpha_pred"]) and r["alpha_pred"] != 0
             else np.nan)
    ok = (np.isfinite(r["alpha"]) and np.isfinite(r["alpha_pred"])
          and np.sign(r["alpha"]) == np.sign(r["alpha_pred"])
          and np.isfinite(ratio) and 0.5 <= ratio <= 2.0)
    g5a_ok = g5a_ok and ok
    print(f"    omega={om:.1f} b={bb:+3.0f}: alpha={r['alpha']:+7.3f} "
          f"pred(b_eff={r['b_eff']:.2f})={r['alpha_pred']:+7.3f} "
          f"ratio={ratio:5.3f}  b_eff/lambda={r['b_eff']/r['lam']:.2f}"
          f"  -> {'OK' if ok else 'POZA'}")
print(f"  -> G5a: {'PASS' if g5a_ok else 'FAIL'}")

print("G5b (RDZEN — zero kanalu cyrkulacyjnego, P2-b; D18): per para:")
print("    |Delta_meas - Delta_eik| <= 3*sigma_Delta")
g5b_ok = True
fails = []
for om, bb in PAIRS:
    rp = runs[(om, +bb)]
    rm = runs[(om, -bb)]
    Dm = rp["alpha"] - rm["alpha"]
    De = rp["alpha_pred"] - rm["alpha_pred"]
    sD = float(np.sqrt(rp["sigma"] ** 2 + rm["sigma"] ** 2))
    diff = Dm - De
    ok = np.isfinite(diff) and abs(diff) <= 3 * sD
    g5b_ok = g5b_ok and ok
    if not ok:
        fails.append(((om, bb), diff))
    print(f"    omega={om:.1f} |b|={bb:2.0f}: Delta_meas={Dm:+7.3f}  "
          f"Delta_eik={De:+7.3f}  roznica={diff:+7.3f}  "
          f"3*sigma_Delta={3*sD:6.3f}  -> {'OK' if ok else 'FAIL'}")
print(f"  -> G5b: {'PASS' if g5b_ok else 'FAIL'}")
if not g5b_ok:
    if len(fails) >= 2:
        signs = {np.sign(d) for _, d in fails}
        if len(signs) == 1:
            print("  regula wzmacniajaca (sekcja 4): FAIL w >= 2 parach "
                  "ze SPOJNYM znakiem nadwyzki -> 'kanal cyrkulacyjny "
                  "zmierzony' (ODKRYCIE wg reguly decyzyjnej; "
                  "bez deklaracji GR — forbidden move #6)")
        else:
            print("  regula wzmacniajaca: FAIL w >= 2 parach, znaki "
                  "nadwyzki NIESPOJNE -> raport bez deklaracji odkrycia")
    else:
        print("  regula wzmacniajaca: FAIL w 1 parze -> 'anomalia "
              "pojedynczej pary' (raport, bez deklaracji odkrycia)")

d6 = abs(r0["alpha"] - r_half["alpha"])
thr6 = 0.1 * max(abs(r0["alpha"]), abs(r_half["alpha"]))
g6b = d6 <= thr6
print(f"G6b (liniowosc): |{r0['alpha']:+.3f} - {r_half['alpha']:+.3f}| = "
      f"{d6:.3f} <= {thr6:.3f}: {g6b}  -> {'PASS' if g6b else 'FAIL'}")
print()
print("obserwacje (POZA kryteriami, raport; Delta_meas ZAWSZE obok "
      "Delta_eik — forbidden move #5):")
for r, lbl in ((r_op, "1.1/+28"), (r_om, "1.1/-28")):
    ratio = (r["alpha"] / r["alpha_pred"]
             if np.isfinite(r.get("alpha_pred", np.nan)) and r["alpha_pred"] != 0
             else np.nan)
    print(f"    {lbl}: alpha={r['alpha']:+.3f} pred={r['alpha_pred']:+.3f} "
          f"ratio={ratio:.3f} b_eff/lambda={r['b_eff']/r['lam']:.2f}")
Dm_o = r_op["alpha"] - r_om["alpha"]
De_o = r_op["alpha_pred"] - r_om["alpha_pred"]
sD_o = float(np.sqrt(r_op["sigma"] ** 2 + r_om["sigma"] ** 2))
print(f"    para |b|=28 (obs): Delta_meas={Dm_o:+.3f}  "
      f"Delta_eik={De_o:+.3f}  sigma_Delta={sD_o:.3f} "
      f"(NIEDOLACZALNA do G5b — forbidden move #4)")
print()
print(f"PODSUMOWANIE PHASE 2: G1(czesc)={'PASS' if g1 else 'FAIL'} "
      f"G4={'PASS' if g4 else 'FAIL'} G5a={'PASS' if g5a_ok else 'FAIL'} "
      f"G5b={'PASS' if g5b_ok else 'FAIL'} "
      f"G6b={'PASS' if g6b else 'FAIL'}")
print("(G6a — symetrie: Phase3_symmetry.py)")
