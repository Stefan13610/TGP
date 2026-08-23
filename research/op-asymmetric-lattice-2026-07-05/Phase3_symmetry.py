# Phase 3: symetrie tla — G6a (S-inv, S-conj)
# Implementacja SCISLE wg Phase0_balance.md proby #4 (LOCKED;
# sekcje 1.2, 2, 4; D11, D17).
# Pre-rejestracja (sekcja 1.2 + D11): obie rownosci to testy
# IMPLEMENTACJI (dokladne symetrie KONFIGURACJI tla), nie fizyki;
# wynik G6a NIE jest pomiarem P2 (forbidden move #7).
#
# Runy:
#   baseline: (omega=1.1, b=+14) prawobiezny — powtorka runu Phase 2;
#   S-inv:    (omega=1.1, b=-14) z pakietem LEWOBIEZNYM (D17):
#             start (-28, -78), kierunek -x, okno lustrzane wzgledem
#             inwersji wokol rdzenia A, centroid x periodyczny;
#   S-conj:   tlo sprzezone + pakiet przesuniety o (128,0),
#             zaimplementowane w ukladzie przesunietym (D11b):
#             tlo_c = np.roll(conj(psi_bg), N/2, axis=0), pakiet
#             standardowy D2 w (-100, -50) ukladu przesunietego.
# Kryterium G6a: |alpha(+14,prawy) - alpha(-14,lewy)| <= 0.05*|alpha(+14)|
#           ORAZ S-conj bitowo / w tej samej tolerancji.

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


def deflection_run_right(psi_bg, om, b_nom, u0=1e-3, tag=""):
    """Run PRAWOBIEZNY jak Phase 2 (D1-D5, D8, D9), bez eikonalu
    (G6a porownuje alpha miedzy wariantami)."""
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    max_steps = 6000
    th_bg = np.angle(psi_bg)
    eith = np.exp(1j * th_bg)
    x0, y0 = X0, Y_ROW + b_nom
    env = u0 * np.exp(-((X - x0) ** 2) / (2 * 8.0 ** 2)
                      - ((Y - y0) ** 2) / (2 * 8.0 ** 2))
    psi_w = psi_bg + eith * env * np.cos(k0 * (X - x0))
    psi_wp = psi_bg + eith * env * np.cos(k0 * (X - x0) + om * dt)
    psi_r = psi_bg.copy()
    psi_rp = psi_bg.copy()

    flags = {"nonfinite": False, "runaway": False, "wrap": False,
             "no5": False, "core_lost": False}
    Es = []
    b_eff = np.inf
    triggered, tau_trig = False, None
    x_thr1, x_thr2, m_win = None, None, None
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
            pos_sc = cp[0][(+1, 0)]
            pos_x2 = cp[0][(-1, 0)]
            n_tot_end = cp[1]
            ut_w = (psi_wn - psi_wp) / (2 * dt)
            Es.append(energy(psi_w, ut_w))
            if not triggered:
                b_eff = min(b_eff, pdist((xc, yc), pos_sc))
                if xc > pos_sc[0] + 10.0:
                    triggered, tau_trig = True, tau
                    x_thr1 = pos_sc[0] + 10.0
                    x_thr2 = pos_x2[0] - 10.0
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
                next_s += 4.0
            if len(samples) >= 5:
                psi_wp, psi_w = psi_w, psi_wn
                psi_rp, psi_r = psi_r, psi_rn
                break
        psi_wp, psi_w = psi_w, psi_wn
        psi_rp, psi_r = psi_r, psi_rn

    out = {"tag": tag, "b_nom": b_nom, "flags": flags,
           "b_eff": b_eff if np.isfinite(b_eff) else np.nan,
           "tau_trig": tau_trig, "n_tot_end": n_tot_end,
           "samples": samples}
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
    return out


def deflection_run_left(psi_bg, om, b_nom, u0=1e-3, tag=""):
    """Run LEWOBIEZNY (D17; wylacznie G6a-S-inv): start (-28, -64+b_nom)
    (dla b_nom=-14: (-28,-78)), kierunek -x; centroid x PERIODYCZNY;
    wyzwolenie s=(x1-x_c) mod L in (10,60); okno lustrzane wzgledem
    inwersji wokol A; alpha := -sign(b_nom)*atan2(Py, -Px)."""
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    max_steps = 6000
    th_bg = np.angle(psi_bg)
    eith = np.exp(1j * th_bg)
    x0L, y0L = -28.0, Y_ROW + b_nom
    env = u0 * np.exp(-((X - x0L) ** 2) / (2 * 8.0 ** 2)
                      - ((Y - y0L) ** 2) / (2 * 8.0 ** 2))
    psi_w = psi_bg + eith * env * np.cos(k0 * (x0L - X))
    psi_wp = psi_bg + eith * env * np.cos(k0 * (x0L - X) + om * dt)  # D17
    psi_r = psi_bg.copy()
    psi_rp = psi_bg.copy()

    flags = {"nonfinite": False, "runaway": False, "wrap": False,
             "no5": False, "core_lost": False}
    Es = []
    b_eff = np.inf
    triggered, tau_trig = False, None
    x1f, x2f, m_win = None, None, None
    samples, next_s = [], None
    n_tot_end = None
    EIX = np.exp(1j * TWO_PI * X / L)     # centroid kolowy w x (D17)

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
            xc = float(np.angle(np.sum(w2 * EIX)) * L / TWO_PI)
            yc = float((w2 * Y).sum() / Wtot)
            if pdist((xc, yc), (x0L, y0L)) > 0.9 * L:
                flags["wrap"] = True
                break
            cp = core_positions(psi_r)
            if cp is None:
                flags["core_lost"] = True
                break
            pos_sc = cp[0][(+1, 0)]           # rdzen A
            pos_x2 = cp[0][(-1, 0)]           # drugi wir rzedu (w -x)
            n_tot_end = cp[1]
            ut_w = (psi_wn - psi_wp) / (2 * dt)
            Es.append(energy(psi_w, ut_w))
            if not triggered:
                b_eff = min(b_eff, pdist((xc, yc), pos_sc))
                s_rel = (pos_sc[0] - xc) % L
                if 10.0 < s_rel < 60.0:       # wyzwolenie (D17)
                    triggered, tau_trig = True, tau
                    x1f = pos_sc[0]                       # ZAMROZONE
                    x2f = pos_x2[0]
                    dist_row = (x1f - x2f) % L
                    sarr = (x1f - X) % L
                    m_win = (sarr > 10.0) & (sarr < dist_row - 10.0)
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

    out = {"tag": tag, "b_nom": b_nom, "flags": flags,
           "b_eff": b_eff if np.isfinite(b_eff) else np.nan,
           "tau_trig": tau_trig, "n_tot_end": n_tot_end,
           "okno": (x1f, x2f)}
    sgn = -1.0 if b_nom > 0 else +1.0
    if len(samples) == 5:
        Pxm = float(np.mean([s[1] for s in samples]))
        Pym = float(np.mean([s[2] for s in samples]))
        out["alpha"] = sgn * float(np.degrees(np.arctan2(Pym, -Pxm)))  # D17
        al_s = [sgn * float(np.degrees(np.arctan2(s[2], -s[1])))
                for s in samples]
        out["alpha_samples"] = al_s
        out["sigma"] = float(np.std(al_s))
    else:
        flags["no5"] = True
        out["alpha"], out["sigma"], out["alpha_samples"] = np.nan, np.nan, []
    q = max(len(Es) // 4, 1)
    out["E_drift"] = abs(np.mean(Es[-q:]) - np.mean(Es[:q])) / abs(Es[0])
    return out


def report(r):
    print(f"  run [{r['tag']}] b_nom={r['b_nom']:+.0f}:")
    print(f"    b_eff = {r['b_eff']:.3f}  tau_wyzw = {r['tau_trig']}")
    if r["alpha_samples"]:
        print(f"    probki alpha [deg]: "
              f"{[f'{x:+.4f}' for x in r['alpha_samples']]}")
    print(f"    alpha = {r['alpha']:+.5f} +/- {r['sigma']:.5f} deg;  "
          f"E_dryft = {r['E_drift']:.2e}; n_tot = {r['n_tot_end']}; "
          f"flagi: {[k for k, v in r['flags'].items() if v] or 'brak'}")


print("=== Phase 3: symetrie tla — G6a (LOCK proby #4: sekcje 1.2/4, "
      "D11, D17) ===")
print("(pre-rejestracja: testy IMPLEMENTACJI — S-inv i S-conj sa "
      "dokladnymi symetriami KONFIGURACJI tla; nie pomiar P2)")
print()

psi = make_lattice()
for stp in range(1, 2001):
    psi = psi + dt_flow * rhs(psi)
psi_bg = psi
print(f"tlo zrelaksowane; residuum max|RHS| = "
      f"{float(np.max(np.abs(rhs(psi_bg)))):.3e}")
# diagnostyki symetrii tla (D11): surowe i gaugeo-swiadome
idxI = (256 - np.arange(N)) % N
inv_field = psi_bg[idxI][:, idxI]
d_inv_raw = float(np.max(np.abs(psi_bg - inv_field)))
d_inv_amp = float(np.max(np.abs(np.abs(psi_bg) - np.abs(inv_field))))
ph_i = np.angle(psi_bg * np.conj(inv_field))
d_inv_ph = float(np.std(wrap(ph_i - np.median(ph_i))))
tlo_c = np.roll(np.conj(psi_bg), N // 2, axis=0)
d_cj_raw = float(np.max(np.abs(tlo_c - psi_bg)))
d_cj_amp = float(np.max(np.abs(np.abs(tlo_c) - np.abs(psi_bg))))
ph_c = np.angle(psi_bg * np.conj(tlo_c))
d_cj_ph = float(np.std(wrap(ph_c - np.median(ph_c))))
print(f"  D11(a) inwersja wokol nominalnego A: max|psi - I[psi]| = "
      f"{d_inv_raw:.3e}; max||psi|-|I psi|| = {d_inv_amp:.3e}; "
      f"std fazy wzglednej = {d_inv_ph:.3e} rad")
print(f"  D11(b) tlo_c = roll(conj(psi_bg), N/2): max|tlo_c - psi_bg| = "
      f"{d_cj_raw:.3e}; max||tlo_c|-|psi_bg|| = {d_cj_amp:.3e}; "
      f"std fazy wzglednej = {d_cj_ph:.3e} rad")
print("  (roznica czysto fazowo-stala = ta sama fizyka: dynamika "
      "U(1)-niezmiennicza; pomiar P na delta_psi gaugeo-niezalezny)")
print()

print("--- baseline: (omega=1.1, b=+14) prawobiezny — powtorka "
      "runu Phase 2 ---")
r_base = deflection_run_right(psi_bg, 1.1, +14.0, tag="baseline")
report(r_base)
print()
print("--- S-inv: (omega=1.1, b=-14) pakiet LEWOBIEZNY (D17) ---")
r_inv = deflection_run_left(psi_bg, 1.1, -14.0, tag="S-inv lewy")
report(r_inv)
print()
print("--- S-conj: tlo sprzezone (= translacja o (128,0)), pakiet "
      "przesuniety; uklad przesuniety (D11b) ---")
r_cj = deflection_run_right(tlo_c, 1.1, +14.0, tag="S-conj")
report(r_cj)
# kontrola bitowa probek S-conj vs baseline
if r_base["samples"] and r_cj["samples"]:
    dPmax = max(max(abs(sb[1] - sc[1]), abs(sb[2] - sc[2]))
                for sb, sc in zip(r_base["samples"], r_cj["samples"]))
    print(f"    kontrola bitowa: max roznica probek (Px,Py) vs baseline "
          f"= {dPmax:.3e} (0.0 = bitowo)")
print()

print("=== EWALUACJA (LOCKED) ===")
a0 = r_base["alpha"]
d_si = abs(a0 - r_inv["alpha"])
d_sc = abs(a0 - r_cj["alpha"])
thr = 0.05 * abs(a0)
ok_si = d_si <= thr
ok_sc = d_sc <= thr
g1_p3 = (all(not any(r["flags"].values()) for r in (r_base, r_inv, r_cj))
         and all(r["E_drift"] <= 1e-4 for r in (r_base, r_inv, r_cj))
         and all(r["n_tot_end"] == 0 for r in (r_base, r_inv, r_cj)))
print(f"G1 (czesc Phase3): flagi/E_dryft/kret  -> "
      f"{'PASS' if g1_p3 else 'FAIL'}")
print(f"G6a-S-inv: |alpha(+14,prawy) - alpha(-14,lewy)| = |{a0:+.5f} - "
      f"{r_inv['alpha']:+.5f}| = {d_si:.2e} <= {thr:.4f}: {ok_si}")
print(f"G6a-S-conj: |alpha(+14) - alpha(conj+shift)| = |{a0:+.5f} - "
      f"{r_cj['alpha']:+.5f}| = {d_sc:.2e} <= {thr:.4f}: {ok_sc}")
print(f"  -> G6a: {'PASS' if ok_si and ok_sc else 'FAIL'}")
print()
print("UWAGA (pre-rejestrowana, forbidden move #7): G6a to testy "
      "implementacji na symetriach dokladnych tla; jedyny pomiar P2 "
      "to G5b na parach +/-b prawobieznych (Phase 2)")
print()
print(f"PODSUMOWANIE PHASE 3: G6a={'PASS' if ok_si and ok_sc else 'FAIL'}")
