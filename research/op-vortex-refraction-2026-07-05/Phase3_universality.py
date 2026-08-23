# Phase 3: uniwersalnosc w znaku kretu — G5 + alpha_odd
# Implementacja SCISLE wg Phase0_balance.md (LOCKED; #1.3, sekcja 2, D1-D5, D11).
#
# Cztery warianty przy (omega=1.1, |b|=8):
#   A: tlo zwykle  (+1 w (-32,-32)), b = +8
#   B: tlo zwykle,                    b = -8
#   C: tlo lustrzane (-1 w (-32,-32)), b = +8
#   D: tlo lustrzane,                  b = -8
# alpha_even = [alpha(+b) + alpha(-b)]/2 (konwencja: alpha>0 = KU wirowi)
# alpha_odd  = [alpha(+b) - alpha(-b)]/2 (cyrkulacyjna; RAPORT, poza kryterium)
# G5: |alpha_even(zwykle) - alpha_even(lustro)| <= 0.2 * |srednia|
# D11 (pre-rejestracja): warianty lustrzane = sprzezenie zespolone ->
#   oczekiwane alpha(C)=alpha(A), alpha(D)=alpha(B) z dokladnoscia zaokraglen.

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt_flow = 0.02
dt = 0.05
N, L = 256, 128.0
dx = L / N
TWO_PI = 2.0 * np.pi
S_GUARD = 10.0

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)
Phi_star = s_star ** 2
Vpp_star = a - 2.0 * b_par * s_star + 3.0 * c * s_star ** 2
xi2 = kappa / Vpp_star

coords = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
Z = X + 1j * Y
V1 = (-32.0, -32.0)
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
        th = -th                      # D11
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


def trace_bg(n2g, b_signed):
    n2, gxg, gyg = n2g
    x, y = -60.0, V1[1] + b_signed
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


def deflection_run(psi_bg, n2cache, om, b_nom, u0=1e-3, scatter_sign=+1,
                   tag=""):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    max_steps = 6000
    th_bg = np.angle(psi_bg)
    eith = np.exp(1j * th_bg)
    x0, y0 = -60.0, V1[1] + b_nom
    env = u0 * np.exp(-((X - x0) ** 2) / (2 * 8.0 ** 2)
                      - ((Y - y0) ** 2) / (2 * 8.0 ** 2))
    psi_w = psi_bg + eith * env * np.cos(k0 * (X - x0))
    psi_wp = psi_bg + eith * env * np.cos(k0 * (X - x0) + om * dt)
    psi_r = psi_bg.copy()
    psi_rp = psi_bg.copy()

    flags = {"nonfinite": False, "runaway": False, "wrap": False,
             "no5": False, "core_lost": False}
    Es = []
    lncore_series = []
    core_series = []
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
            core_series.append((tau, pos_sc))
            n_tot_end = cp[1]
            m_core = np.zeros((N, N), bool)
            for s_ in (+1, -1):
                px_, py_ = cp[0][s_]
                rr = np.hypot((X - px_ + L / 2) % L - L / 2,
                              (Y - py_ + L / 2) % L - L / 2)
                m_core |= rr < 5.0
            lncore_series.append(
                (tau, float(np.log(np.sqrt(np.mean(w2[m_core])) + 1e-300))))
            ut_w = (psi_wn - psi_wp) / (2 * dt)
            Es.append(energy(psi_w, ut_w))
            if not triggered:
                b_eff = min(b_eff, pdist((xc, yc), pos_sc))
                thr = pos_sc[0] + 10.0
                if float(w2[X > thr].sum()) > 0.5 * Wtot:
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
                break
        psi_wp, psi_w = psi_w, psi_wn
        psi_rp, psi_r = psi_r, psi_rn

    out = {"om": om, "b_nom": b_nom, "tag": tag, "flags": flags,
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
    out["core_disp"] = pdist(core_series[0][1], core_series[-1][1])
    if np.isfinite(out["b_eff"]):
        out["alpha_pred"] = trace_bg(n2cache, np.sign(b_nom) * out["b_eff"])
    else:
        out["alpha_pred"] = np.nan
    return out


def report(r):
    print(f"  wariant [{r['tag']}] b_nom={r['b_nom']:+.0f}:")
    print(f"    b_eff = {r['b_eff']:.3f}  tau_wyzw = {r['tau_trig']}  "
          f"dryf rdzenia {r['core_disp']:.2f}")
    if r["alpha_samples"]:
        print(f"    probki alpha [deg]: "
              f"{[f'{x:+.2f}' for x in r['alpha_samples']]}")
    print(f"    alpha = {r['alpha']:+.4f} +/- {r['sigma']:.4f} deg   "
          f"alpha_pred(b_eff) = {r['alpha_pred']:+.3f} deg")
    print(f"    G1: E_dryft = {r['E_drift']:.2e}; gamma_core(po tranzycie) = "
          f"{r['gamma_core_post']:+.5f}; n_tot = {r['n_tot_end']:.0f}; "
          f"flagi: {[k for k, v in r['flags'].items() if v] or 'brak'}")


print("=== Phase 3: uniwersalnosc G5 (LOCK: Phase0_balance.md, #1.3, D11) ===")
print("omega = 1.1, |b| = 8; cztery warianty A-D")
print()

OM = 1.1
results = {}
for mirror, lbl in ((False, "zwykle(+1)"), (True, "lustro(-1)")):
    psi = make_pair(mirror=mirror)
    for stp in range(1, 2001):
        psi = psi + dt_flow * rhs(psi)
    psi_bg = psi
    n2c = build_n2(psi_bg, OM)
    ssign = -1 if mirror else +1
    print(f"--- tlo {lbl}: zrelaksowane; residuum = "
          f"{float(np.max(np.abs(rhs(psi_bg)))):.3e} ---")
    for bb in (+8.0, -8.0):
        key = ("M" if mirror else "N", "+" if bb > 0 else "-")
        r = deflection_run(psi_bg, n2c, OM, bb, scatter_sign=ssign,
                           tag=f"{lbl} b={bb:+.0f}")
        results[key] = r
        report(r)
    print()

# ---------------- EWALUACJA ----------------
print("=== EWALUACJA (LOCKED) ===")
aN_p, aN_m = results[("N", "+")]["alpha"], results[("N", "-")]["alpha"]
aM_p, aM_m = results[("M", "+")]["alpha"], results[("M", "-")]["alpha"]
even_N, odd_N = 0.5 * (aN_p + aN_m), 0.5 * (aN_p - aN_m)
even_M, odd_M = 0.5 * (aM_p + aM_m), 0.5 * (aM_p - aM_m)
print(f"tlo zwykle (+1):  alpha(+8) = {aN_p:+.4f}, alpha(-8) = {aN_m:+.4f}")
print(f"                  alpha_even = {even_N:+.4f}, "
      f"alpha_odd = {odd_N:+.4f} (cyrkulacyjna, POZA kryterium)")
print(f"tlo lustrzane:    alpha(+8) = {aM_p:+.4f}, alpha(-8) = {aM_m:+.4f}")
print(f"                  alpha_even = {even_M:+.4f}, "
      f"alpha_odd = {odd_M:+.4f}")
print(f"D11 (pre-rejestracja, sprzezenie): |alpha_A - alpha_C| = "
      f"{abs(aN_p - aM_p):.2e}, |alpha_B - alpha_D| = {abs(aN_m - aM_m):.2e}")
mean_even = 0.5 * (abs(even_N) + abs(even_M))
diff_even = abs(even_N - even_M)
g5 = diff_even <= 0.2 * mean_even
print(f"G5: |alpha_even(zwykle) - alpha_even(lustro)| = {diff_even:.4f} "
      f"<= 0.2*|srednia| = {0.2 * mean_even:.4f}: {g5}  "
      f"-> {'PASS' if g5 else 'FAIL'}")
g1_ok = all(not any(r["flags"].values()) and r["E_drift"] <= 1e-4
            and r["gamma_core_post"] <= 0.01 and r["n_tot_end"] == 0
            for r in results.values())
print(f"G1 (czesc Phase3): flagi/E_dryft/gamma_core/kret: "
      f"{'PASS' if g1_ok else 'FAIL'}")
print()
print(f"PODSUMOWANIE PHASE 3: G5={'PASS' if g5 else 'FAIL'} "
      f"G1(czesc)={'PASS' if g1_ok else 'FAIL'}")
