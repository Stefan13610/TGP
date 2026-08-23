# Phase 3: symetrie — G6a (LOCK proby #5: sekcja 1.2/4, D11, D16)
# Pre-rejestracja (D11): testy IMPLEMENTACJI ramki skosnej, nie
# fizyki; wynik nie jest pomiarem P2 (forbidden move #7).
#
# Runy: baseline (omega=1.1, b=+14) — powtorka runu Phase 2
#       (pakiet standardowy D2, przod +u);
#       S-inv: (1.1, b=-14) pakiet WSTECZNY D16 (-u, start
#       r0L = A_nom + 36u - 14n, faza -omega*dt, wyzwolenie
#       xi_c < -10, okno {-130 < xi < -10} - dyski) — inwersja
#       wokol A jest dokladna symetria szachownicy;
#       S-conj: tlo lustrzane (theta -> -theta), pakiet standardowy
#       b=+14 — oczekiwanie bitowe (wzorzec D11 prob #2/#3).
# Kryterium G6a: |alpha(+14,+u) - alpha(-14,-u)| <= 0.05*|alpha(+14)|
#           ORAZ |alpha(lustro,+14) - alpha(+14)| <= 0.05*|alpha(+14)|.

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt_flow = 0.02
dt = 0.05
N, L = 512, 256.0
dx = L / N
TWO_PI = 2.0 * np.pi
S_GUARD = 10.0
SEED = 20260704
S5 = np.sqrt(5.0)

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)
Phi_star = s_star ** 2
Vpp_star = a - 2.0 * b_par * s_star + 3.0 * c * s_star ** 2
xi2 = kappa / Vpp_star

coords = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
Z = X + 1j * Y

VORTS = [((-64.0, -64.0), +1), ((+64.0, -64.0), -1),
         ((-64.0, +64.0), -1), ((+64.0, +64.0), +1)]
NOM = {+1: [(-64.0, -64.0), (+64.0, +64.0)],
       -1: [(+64.0, -64.0), (-64.0, +64.0)]}
A_NOM = (-64.0, -64.0)
UVEC = np.array([2.0, 1.0]) / S5
NVEC = np.array([-1.0, 2.0]) / S5
XI_START = -36.0


def gfun(m):
    return a - b_par * np.abs(m) + c * m * m


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


def xi_eta(x, y, A):
    dxp = (x - A[0] + 0.5 * L) % L - 0.5 * L
    dyp = (y - A[1] + 0.5 * L) % L - 0.5 * L
    return (2.0 * dxp + dyp) / S5, (-dxp + 2.0 * dyp) / S5


def theta_h_pair(z1, z2):
    th = np.zeros_like(X)
    for (x0, y0), n in ((z1, +1), (z2, -1)):
        for k in range(-2, 3):
            th = th + n * np.angle(np.sinh(np.pi * (Z - (x0 + k * L) - 1j * y0) / L))
    return th


def make_lattice(mirror=False):
    th = theta_h_pair((-64.0, -64.0), (+64.0, -64.0)) \
        + theta_h_pair((+64.0, +64.0), (-64.0, +64.0))
    if mirror:
        th = -th
    Amp = np.full_like(X, s_star)
    for (x0, y0), _n in VORTS:
        ddx = (X - x0 + 0.5 * L) % L - 0.5 * L
        ddy = (Y - y0 + 0.5 * L) % L - 0.5 * L
        rr = np.sqrt(ddx * ddx + ddy * ddy)
        Amp = Amp * rr / np.sqrt(rr * rr + 2.0 * xi2)
    return Amp * np.exp(1j * th)


def core_positions(psi, nom=NOM):
    """Detektor 2+2 (D14). Tlo lustrzane: znaki nominalow odwrocone."""
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
        if len(clusters) != 2:
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
        d00 = pdist(pos[0], nom[sign][0]) + pdist(pos[1], nom[sign][1])
        d01 = pdist(pos[0], nom[sign][1]) + pdist(pos[1], nom[sign][0])
        if d01 < d00:
            pos = [pos[1], pos[0]]
        out[(sign, 0)] = pos[0]
        out[(sign, 1)] = pos[1]
    ntot = float(np.rint(W).sum())
    return out, ntot


def deflection_run(psi_bg, om, b_nom, u0=1e-3, tag="", mode="fwd",
                   sc_sign=+1):
    """Run ugiecia w ramce skosnej (D1-D5, D8, D9; bez eikonalu —
    G6a porownuje alpha miedzy wariantami).
    mode='fwd': pakiet D2 (+u); mode='bwd': pakiet WSTECZNY D16 (-u).
    sc_sign: znak rdzenia A w polu (tlo lustrzane: -1)."""
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    max_steps = 6000
    th_bg = np.angle(psi_bg)
    eith = np.exp(1j * th_bg)
    if mode == "fwd":
        r0 = np.array(A_NOM) + XI_START * UVEC + b_nom * NVEC
        ph = +om * dt
    else:
        r0 = np.array(A_NOM) + 36.0 * UVEC + b_nom * NVEC       # D16
        ph = -om * dt
    env = u0 * np.exp(-((X - r0[0]) ** 2 + (Y - r0[1]) ** 2)
                      / (2 * 8.0 ** 2))
    s_par = (2.0 * (X - r0[0]) + (Y - r0[1])) / S5
    psi_w = psi_bg + eith * env * np.cos(k0 * s_par)
    psi_wp = psi_bg + eith * env * np.cos(k0 * s_par + ph)
    psi_r = psi_bg.copy()
    psi_rp = psi_bg.copy()
    nom = NOM if sc_sign == +1 else {s: NOM[-s] for s in NOM}

    flags = {"nonfinite": False, "runaway": False, "wrap": False,
             "no5": False, "core_lost": False}
    Es = []
    b_eff = np.inf
    triggered, tau_trig = False, None
    m_win = None
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
            if np.hypot(xc - r0[0], yc - r0[1]) > 0.9 * L:
                flags["wrap"] = True
                break
            cp = core_positions(psi_r, nom)
            if cp is None:
                flags["core_lost"] = True
                break
            pos_A = cp[0][(sc_sign, 0)]
            n_tot_end = cp[1]
            xi_c, eta_c = xi_eta(xc, yc, pos_A)
            ut_w = (psi_wn - psi_wp) / (2 * dt)
            Es.append(energy(psi_w, ut_w))
            if not triggered:
                b_eff = min(b_eff, pdist((xc, yc), pos_A))
                trig_now = (xi_c > 10.0) if mode == "fwd" else (xi_c < -10.0)
                if trig_now:
                    triggered, tau_trig = True, tau
                    DXp = (X - pos_A[0] + 0.5 * L) % L - 0.5 * L
                    DYp = (Y - pos_A[1] + 0.5 * L) % L - 0.5 * L
                    XIg = (2.0 * DXp + DYp) / S5
                    if mode == "fwd":
                        m_win = (XIg > 10.0) & (XIg < 130.0)
                    else:
                        m_win = (XIg > -130.0) & (XIg < -10.0)
                    for key in cp[0]:
                        px_, py_ = cp[0][key]
                        rr = np.hypot((X - px_ + L / 2) % L - L / 2,
                                      (Y - py_ + L / 2) % L - L / 2)
                        m_win &= rr >= 12.0
                    next_s = tau
            if triggered and next_s is not None and tau >= next_s \
                    and len(samples) < 5:
                ut_d = ((psi_wn - psi_rn) - (psi_wp - psi_rp)) / (2 * dt)
                gx = (np.roll(dpsi, -1, 0) - np.roll(dpsi, 1, 0)) / (2 * dx)
                gy = (np.roll(dpsi, -1, 1) - np.roll(dpsi, 1, 1)) / (2 * dx)
                Px = float(-np.real(np.sum(np.conj(ut_d)[m_win] * gx[m_win])) * dx * dx)
                Py = float(-np.real(np.sum(np.conj(ut_d)[m_win] * gy[m_win])) * dx * dx)
                samples.append((tau, Px, Py))
                next_s += 4.0
            if len(samples) >= 5:
                psi_wp, psi_w = psi_w, psi_wn
                psi_rp, psi_r = psi_r, psi_rn
                break
        psi_wp, psi_w = psi_w, psi_wn
        psi_rp, psi_r = psi_r, psi_rn

    out = {"tag": tag, "b_nom": b_nom, "mode": mode, "flags": flags,
           "b_eff": b_eff if np.isfinite(b_eff) else np.nan,
           "tau_trig": tau_trig, "n_tot_end": n_tot_end}
    sgn = -1.0 if b_nom > 0 else +1.0
    if len(samples) == 5:
        Pxm = float(np.mean([s[1] for s in samples]))
        Pym = float(np.mean([s[2] for s in samples]))
        Pu = (2.0 * Pxm + Pym) / S5
        Pn = (-Pxm + 2.0 * Pym) / S5
        if mode == "fwd":
            out["alpha"] = sgn * float(np.degrees(np.arctan2(Pn, Pu)))
        else:
            out["alpha"] = sgn * float(np.degrees(np.arctan2(Pn, -Pu)))
        al_s = []
        for _t, sPx, sPy in samples:
            su = (2.0 * sPx + sPy) / S5
            sn = (-sPx + 2.0 * sPy) / S5
            if mode == "fwd":
                al_s.append(sgn * float(np.degrees(np.arctan2(sn, su))))
            else:
                al_s.append(sgn * float(np.degrees(np.arctan2(sn, -su))))
        out["alpha_samples"] = al_s
        out["sigma"] = float(np.std(al_s))
    else:
        flags["no5"] = True
        out["alpha"], out["sigma"], out["alpha_samples"] = np.nan, np.nan, []
    q = max(len(Es) // 4, 1)
    out["E_drift"] = abs(np.mean(Es[-q:]) - np.mean(Es[:q])) / abs(Es[0])
    return out


def report(r):
    print(f"  run [{r['tag']}] b_nom={r['b_nom']:+.0f} ({r['mode']}):")
    print(f"    b_eff = {r['b_eff']:.3f}  tau_wyzw = {r['tau_trig']}")
    if r["alpha_samples"]:
        print(f"    probki alpha [deg]: "
              f"{[f'{x:+.4f}' for x in r['alpha_samples']]}")
    print(f"    alpha = {r['alpha']:+.5f} +/- {r['sigma']:.5f} deg;  "
          f"E_dryft = {r['E_drift']:.2e}; n_tot = {r['n_tot_end']}; "
          f"flagi: {[k for k, v in r['flags'].items() if v] or 'brak'}")


print("=== Phase 3: symetrie — G6a (LOCK proby #5: sekcja 1.2/4, "
      "D11, D16) ===")
print("(pre-rejestracja: testy IMPLEMENTACJI ramki skosnej; S-inv: "
      "zaokraglenia, S-conj: bitowo; NIE pomiar P2)")
print()

psi = make_lattice()
for stp in range(1, 2001):
    psi = psi + dt_flow * rhs(psi)
psi_bg = psi
print(f"tlo zrelaksowane; residuum max|RHS| = "
      f"{float(np.max(np.abs(rhs(psi_bg)))):.3e}")
psi_m = make_lattice(mirror=True)
for stp in range(1, 2001):
    psi_m = psi_m + dt_flow * rhs(psi_m)
psi_bg_m = psi_m
print(f"tlo lustrzane zrelaksowane; residuum = "
      f"{float(np.max(np.abs(rhs(psi_bg_m)))):.3e}; "
      f"kontrola D11b: max|psi_bg_m - conj(psi_bg)| = "
      f"{float(np.max(np.abs(psi_bg_m - np.conj(psi_bg)))):.3e}")
# diagnostyka D11a: inwersja indeksow wokol nominalnego A (-64,-64)
# (i,j) -> (256 - i) mod 512: I = roll(psi[::-1,::-1], 257, obie osie)
psi_inv = np.roll(psi_bg[::-1, ::-1], (257, 257), axis=(0, 1))
print(f"diagnostyka D11a: max|psi_bg - I[psi_bg]| (inwersja wokol "
      f"A_nom) = {float(np.max(np.abs(psi_bg - psi_inv))):.3e} "
      f"(ansatz z obcieciem obrazow: oczekiwane niezerowe male)")
print()

print("--- baseline: (omega=1.1, b=+14) pakiet standardowy D2 "
      "(powtorka runu Phase 2) ---")
r_base = deflection_run(psi_bg, 1.1, +14.0, tag="baseline")
report(r_base)
print()
print("--- S-inv: (omega=1.1, b=-14) pakiet WSTECZNY D16 (-u) ---")
r_inv = deflection_run(psi_bg, 1.1, -14.0, tag="S-inv bwd", mode="bwd")
report(r_inv)
print()
print("--- S-conj: tlo lustrzane (theta -> -theta), b=+14, pakiet "
      "standardowy ---")
r_mir = deflection_run(psi_bg_m, 1.1, +14.0, tag="S-conj lustro",
                       sc_sign=-1)
report(r_mir)
print()

print("=== EWALUACJA (LOCKED) ===")
a0 = r_base["alpha"]
d_si = abs(a0 - r_inv["alpha"])
d_sc = abs(a0 - r_mir["alpha"])
thr = 0.05 * abs(a0)
ok_si = d_si <= thr
ok_sc = d_sc <= thr
g1_p3 = (all(not any(r["flags"].values()) for r in (r_base, r_inv, r_mir))
         and all(r["E_drift"] <= 1e-4 for r in (r_base, r_inv, r_mir))
         and all(r["n_tot_end"] == 0 for r in (r_base, r_inv, r_mir)))
print(f"G1 (czesc Phase3): flagi/E_dryft/kret  -> "
      f"{'PASS' if g1_p3 else 'FAIL'}")
print(f"G6a-S-inv: |alpha(+14,+u) - alpha(-14,-u)| = |{a0:+.5f} - "
      f"{r_inv['alpha']:+.5f}| = {d_si:.2e} <= {thr:.4f}: {ok_si}")
print(f"G6a-S-conj: |alpha(+14) - alpha(lustro,+14)| = |{a0:+.5f} - "
      f"{r_mir['alpha']:+.5f}| = {d_sc:.2e} <= {thr:.4f}: {ok_sc}")
print(f"  -> G6a: {'PASS' if ok_si and ok_sc else 'FAIL'}")
print()
print("UWAGA (pre-rejestrowana, forbidden move #7): wynik S-inv/S-conj")
print("jest testem implementacji ramki skosnej i sprzezenia — NIE")
print("pomiarem P2 (jedyny pomiar P2 to G5b na parach +/-b, Phase 2)")
print()
print(f"PODSUMOWANIE PHASE 3: G6a={'PASS' if ok_si and ok_sc else 'FAIL'}")
