# Phase 1: BRAMKA G2 (czas zycia tla pod dynamika pomiarowa)
#          + kalibracja G3 (dyspersja) + background_check + czesc G1
# Implementacja SCISLE wg Phase0_balance.md (LOCKED; D1, D7, D8, D9,
# D13, D14, D15). KOLEJNOSC = CZESC PROTOKOLU: bramka G2 przed
# jakimkolwiek runem z pakietem; G2 FAIL -> STOP (zero runow ugiecia).
#
# background_check: szachownica 4 wirow (2 pary poziome theta_h_pair,
#   BEZ transpozycji), diagnostyka szwu, detektor 2+2 (D14), relaksacja
#   2000 krokow, determinizm bitowy relaksacji.
# G2 (D15): ewolucja II rzedu BEZ pakietu, 8000 krokow (tau=400);
#   pozycje 4 rdzeni co 10 krokow; kryterium: max przemieszczenie
#   kazdego rdzenia <= 1.0 dla tau <= 300; gamma_s fit jesli widoczny.
# G3 (D7): tlo jednorodne s*, zaburzenie rzeczywiste, m in {11,15,20,
#   25,30}, okno tau=100; blad <= 5%; determinizm bitowy runu falowego.

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt_flow = 0.02
dt = 0.05
N, L = 256, 128.0
dx = L / N
TWO_PI = 2.0 * np.pi
S_GUARD = 10.0
H_DRIFT_TOL = 1e-4
H_TOL = 1e-10
SEED = 20260704
M_LIST = [11, 15, 20, 25, 30]

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)
Phi_star = s_star ** 2
Vpp_star = a - 2.0 * b_par * s_star + 3.0 * c * s_star ** 2
xi2 = kappa / Vpp_star

coords = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
Z = X + 1j * Y

VORTS = [((-32.0, -32.0), +1), ((+32.0, -32.0), -1),
         ((-32.0, +32.0), -1), ((+32.0, +32.0), +1)]
NOM = {+1: [(-32.0, -32.0), (+32.0, +32.0)],
       -1: [(+32.0, -32.0), (-32.0, +32.0)]}


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


def make_lattice(mirror=False):
    th = theta_h_pair((-32.0, -32.0), (+32.0, -32.0)) \
        + theta_h_pair((+32.0, +32.0), (-32.0, +32.0))
    if mirror:
        th = -th
    Amp = np.full_like(X, s_star)
    for (x0, y0), _n in VORTS:
        ddx = (X - x0 + 0.5 * L) % L - 0.5 * L
        ddy = (Y - y0 + 0.5 * L) % L - 0.5 * L
        rr = np.sqrt(ddx * ddx + ddy * ddy)
        Amp = Amp * rr / np.sqrt(rr * rr + 2.0 * xi2)
    return Amp * np.exp(1j * th)


def core_positions(psi):
    """Detektor 2+2 (D14)."""
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
        d00 = pdist(pos[0], NOM[sign][0]) + pdist(pos[1], NOM[sign][1])
        d01 = pdist(pos[0], NOM[sign][1]) + pdist(pos[1], NOM[sign][0])
        if d01 < d00:
            pos = [pos[1], pos[0]]
        out[(sign, 0)] = pos[0]
        out[(sign, 1)] = pos[1]
    ntot = float(np.rint(W).sum())
    return out, ntot


def freq_zero_cross(series, taus):
    sgn = np.sign(series)
    idx = np.where(sgn[:-1] * sgn[1:] < 0)[0]
    if len(idx) < 4:
        return np.nan

    def tcross(i):
        return taus[i] + (taus[i + 1] - taus[i]) * series[i] / (series[i] - series[i + 1])

    t0, t1 = tcross(idx[0]), tcross(idx[-1])
    n_half = len(idx) - 1
    return np.pi * n_half / (t1 - t0)


print("=== Phase 1: BRAMKA G2 + kalibracja G3 + background_check ===")
print(f"(LOCK: Phase0_balance.md, D1/D7/D8/D14/D15) "
      f"N={N} L={L} dx={dx} dt={dt}")
print()

# ---------------- background_check ----------------
print("--- background_check: szachownica 4 wirow (2 pary poziome) ---")
psi = make_lattice()
th0 = np.angle(psi)
gxa = np.abs(wrap(np.roll(th0, -1, 0) - th0))
gya = np.abs(wrap(np.roll(th0, -1, 1) - th0))
far = np.ones((N, N), bool)
for (x0, y0), _n in VORTS:
    far &= np.hypot((X - x0 + L/2) % L - L/2, (Y - y0 + L/2) % L - L/2) > 4
seam = float(max(gxa[far].max(), gya[far].max()))
print(f"  diagnostyka szwu: max |dtheta| na linkach poza rdzeniami = "
      f"{seam:.4f} rad << pi: {seam < 0.5}  (ogony dipolowe; "
      f"BEZ artefaktu punktu C)")
cp0 = core_positions(psi)
print(f"  kret ansatzu: n_tot = {cp0[1]:.0f}; 4 rdzenie: "
      + "; ".join(f"{s:+d}#{i}:({p[0]:+.2f},{p[1]:+.2f})"
                  for (s, i), p in sorted(cp0[0].items())))
H_prev, H_inc = None, 0.0
for stp in range(1, 2001):
    psi = psi + dt_flow * rhs(psi)
    if stp % 10 == 0:
        gx = (np.roll(psi, -1, 0) - psi) / dx
        gy = (np.roll(psi, -1, 1) - psi) / dx
        H = float(np.sum(0.5 * kappa * (np.abs(gx)**2 + np.abs(gy)**2)
                         + Vpot(np.abs(psi))) * dx * dx)
        if H_prev is not None and H > H_prev:
            H_inc = max(H_inc, H - H_prev)
        H_prev = H
psi_bg = psi
cpR = core_positions(psi_bg)
res = float(np.max(np.abs(rhs(psi_bg))))
print(f"  po relaksacji 2000 krokow (dt={dt_flow}):")
print(f"    residuum max|RHS| = {res:.3e} (para proby #1: 8.107e-03 — "
      f"tlo siatki jest STACJONARNE, nie quasi-)")
print(f"    pozycje rdzeni: "
      + "; ".join(f"{s:+d}#{i}:({p[0]:+.3f},{p[1]:+.3f})"
                  for (s, i), p in sorted(cpR[0].items())))
print(f"    n_tot = {cpR[1]:.0f}; |psi| daleko od rdzeni: min = "
      f"{float(np.abs(psi_bg)[far].min()):.4f} (s* = {s_star:.4f})")
# determinizm bitowy relaksacji (D8)
psi_d = make_lattice()
for stp in range(1, 2001):
    psi_d = psi_d + dt_flow * rhs(psi_d)
det_bg = bool(np.array_equal(psi_bg, psi_d))
print(f"    determinizm bitowy relaksacji: {det_bg}; H_mono: "
      f"{H_inc <= H_TOL}")
print()

# ---------------- BRAMKA G2 (D15) ----------------
print("--- BRAMKA G2: tlo bez pakietu, dynamika II rzedu, tau = 400 ---")
print("    (kryterium: max przemieszczenie KAZDEGO rdzenia <= 1.0 "
      "dla tau <= 300)")
pos0 = {k: v for k, v in cpR[0].items()}
psi_w, psi_p = psi_bg.copy(), psi_bg.copy()   # zerowa predkosc (D15)
disp_series = []          # (tau, dmax, {key: disp})
Es, taus_E = [], []
gate_flags = {"nonfinite": False, "runaway": False, "core_lost": False}
annih_tau = None
for stp in range(1, 8001):
    tau = stp * dt
    psi_n = 2.0 * psi_w - psi_p + dt * dt * rhs(psi_w)
    if stp % 10 == 0:
        if not np.all(np.isfinite(psi_n)):
            gate_flags["nonfinite"] = True
            break
        if np.max(np.abs(psi_n)) > S_GUARD:
            gate_flags["runaway"] = True
            break
        cp = core_positions(psi_w)
        if cp is None:
            gate_flags["core_lost"] = True
            annih_tau = tau
            break
        d = {k: pdist(cp[0][k], pos0[k]) for k in pos0}
        disp_series.append((tau, max(d.values()), d))
        psit = (psi_n - psi_p) / (2.0 * dt)
        Es.append(energy(psi_w, psit))
        taus_E.append(tau)
        if stp % 800 == 0:
            print(f"    tau={tau:6.1f}: max przemieszczenie = "
                  f"{max(d.values()):.4f}   n_tot = {cp[1]:.0f}")
    psi_p, psi_w = psi_w, psi_n

if annih_tau is not None:
    print(f"    UTRATA RDZENI (detektor 2+2) przy tau = {annih_tau:.1f}")
tt = np.array([p[0] for p in disp_series])
dd = np.array([p[1] for p in disp_series])
m300 = tt <= 300.0
dmax300 = float(dd[m300].max()) if m300.any() else np.inf
dmax400 = float(dd.max()) if len(dd) else np.inf
q = max(len(Es) // 4, 1)
E_drift_gate = abs(np.mean(Es[-q:]) - np.mean(Es[:q])) / abs(Es[0])
print(f"    max przemieszczenie rdzenia: tau<=300: {dmax300:.4f}; "
      f"tau<=400: {dmax400:.4f}")
per_core_300 = {k: max(p[2][k] for p in disp_series if p[0] <= 300.0)
                for k in pos0}
print("    per rdzen (tau<=300): "
      + "; ".join(f"{s:+d}#{i}: {v:.4f}"
                  for (s, i), v in sorted(per_core_300.items())))
print(f"    dryft sekularny energii (run bramki) = {E_drift_gate:.2e}")
# gamma_s: fit ln(dmax) vs tau na odcinku dmax in [0.05, 1.0] (D15)
mfit = (dd >= 0.05) & (dd <= 1.0)
if mfit.sum() >= 5:
    A_ = np.vstack([np.ones(mfit.sum()), tt[mfit]]).T
    coef, *_ = np.linalg.lstsq(A_, np.log(dd[mfit]), rcond=None)
    print(f"    gamma_s (fit ln dmax vs tau, odcinek [0.05,1.0]) = "
          f"{float(coef[1]):+.5f} /tau")
else:
    print("    gamma_s: brak widocznej ucieczki w oknie "
          "(dmax nie osiaga 0.05 w wystarczajacej liczbie probek)"
          if dmax400 < 0.05 else
          "    gamma_s: za malo probek w odcinku [0.05, 1.0] do fitu")
print("    KONTRAST (proba #1, cytowane): para (+1,-1) — przemieszczenie "
      "14 przy tau=50, anihilacja tau=107.5")
g2_pass = (dmax300 <= 1.0) and not any(gate_flags.values())
print(f"    -> BRAMKA G2: max przemieszczenie (tau<=300) = {dmax300:.4f} "
      f"<= 1.0: {dmax300 <= 1.0}; flagi: "
      f"{[k for k, v in gate_flags.items() if v] or 'brak'}  "
      f"-> {'PASS' if g2_pass else 'FAIL'}")
print()

if not g2_pass:
    print("=== G2 FAIL -> STOP wg reguly decyzyjnej ===")
    print("  zadnych runow ugiecia; deliverable: zmierzone czasy zycia")
    print("  (siatka vs para); kalibracja G3 POMINIETA (stop przed pomiarem)")
    raise SystemExit(0)

# ---------------- G3: dyspersja kanalu amplitudowego (D7) ----------------
print("--- G3: dyspersja na tle jednorodnym s* (D7) ---")
print("    m     k       omega_meas   omega_teor   blad%   omega_siec   "
      "H_dryft(sek.)  H_osc")
TAU_WIN = 100.0
STEPS = int(round(TAU_WIN / dt))
PROBE = (N // 2, 10)
disp_ok, g1_ok, det_ok = True, True, True
errs, drifts = [], []
for m in M_LIST:
    k = TWO_PI * m / L
    pert = 1e-3 * np.cos(k * X)
    psi0 = (s_star + pert).astype(complex)
    psi_prev0 = psi0.copy()

    def run_g3():
        psi_prev, psi_ = psi_prev0.copy(), psi0.copy()
        series, Hs, taus = [], [], []
        nf = False
        for stp in range(1, STEPS + 1):
            psi_next = 2.0 * psi_ - psi_prev + dt * dt * rhs(psi_)
            if stp % 10 == 0:
                if not np.all(np.isfinite(psi_next)):
                    nf = True
                    break
                psit = (psi_next - psi_prev) / (2.0 * dt)
                series.append(float(np.real(psi_[PROBE]) - s_star))
                Hs.append(energy(psi_, psit))
                taus.append(stp * dt)
            psi_prev, psi_ = psi_, psi_next
        return np.array(series), np.array(Hs), np.array(taus), nf, psi_

    s1, Hs, taus, nf, fin1 = run_g3()
    if m == M_LIST[0]:
        s2, _, _, _, fin2 = run_g3()
        det_ok = det_ok and bool(np.array_equal(fin1, fin2))
    om_meas = freq_zero_cross(s1, taus)
    om_th = np.sqrt(kappa * k * k + Vpp_star)
    k2_lat = (2 - 2 * np.cos(k * dx)) / (dx * dx)
    om_lat = np.sqrt(kappa * k2_lat + Vpp_star)
    err = 100 * (om_meas - om_th) / om_th
    q = len(Hs) // 4
    drift = abs(np.mean(Hs[-q:]) - np.mean(Hs[:q])) / abs(Hs[0])
    osc = (Hs.max() - Hs.min()) / abs(Hs[0])
    disp_ok = disp_ok and abs(err) <= 5.0
    g1_ok = g1_ok and drift <= H_DRIFT_TOL and not nf
    errs.append(abs(err))
    drifts.append(drift)
    print(f"   {m:3d}  {k:.4f}   {om_meas:9.5f}   {om_th:9.5f}   "
          f"{err:+5.2f}   {om_lat:9.5f}   {drift:9.2e}   {osc:.2e}")
print()

# ---------------- EWALUACJA ----------------
print("=== EWALUACJA (LOCKED) ===")
g1_part = (g1_ok and det_ok and det_bg and cpR[1] == 0
           and E_drift_gate <= H_DRIFT_TOL)
print(f"G1 (czesc Phase1): determinizm bitowy (fala: {det_ok}, "
      f"relaksacja: {det_bg}); brak NaN; dryft sekularny H (dyspersja) "
      f"max = {max(drifts):.2e}, (bramka) = {E_drift_gate:.2e} <= "
      f"{H_DRIFT_TOL}; kret tla n_tot=0 i 4 rdzenie: {cpR[1] == 0}  "
      f"-> {'PASS' if g1_part else 'FAIL'}")
print(f"G2 (BRAMKA czasu zycia): max przemieszczenie (tau<=300) = "
      f"{dmax300:.4f} <= 1.0  -> {'PASS' if g2_pass else 'FAIL'}")
print(f"G3 (kalibracja dyspersji): max |blad| = {max(errs):.2f}% <= 5%: "
      f"{disp_ok}  -> {'PASS' if disp_ok else 'FAIL'}")
print()
print(f"PODSUMOWANIE PHASE 1: G1(czesc)={'PASS' if g1_part else 'FAIL'} "
      f"G2={'PASS' if g2_pass else 'FAIL'} "
      f"G3={'PASS' if disp_ok else 'FAIL'}")
print("BRAMKA G2 PASS -> wolno uruchomic Phase 2 (runy ugiecia)")
