# Phase 1: kalibracja G2 (dyspersja kanalu amplitudowego) + background_check
# Implementacja SCISLE wg Phase0_balance.md (LOCKED; D1, D7, D8, D13).
#
# Sektor falowy: PELNY zespolony II rzedu:
#   d^2psi/dtau^2 = kappa*Lap(psi) - psi*(a - b|psi| + c|psi|^2)
# G2: tlo jednorodne s*; zaburzenie rzeczywiste u0*cos(kx), u_tau(0)=0
#     (pole pozostaje rzeczywiste -> czysty kanal amplitudowy);
#     m in {11,15,20,25,30} (D7), omega z przejsc przez zero, okno tau=100.
# background_check: para (+1,-1), ansatz zlozony, relaksacja 2000 krokow;
#     diagnostyka szwu, kret, residuum, dryf rdzeni, determinizm bitowy.

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt_flow = 0.02
dt = 0.05
N, L = 256, 128.0
dx = L / N
TWO_PI = 2.0 * np.pi
TAU_WIN = 100.0
STEPS = int(round(TAU_WIN / dt))
M_LIST = [11, 15, 20, 25, 30]          # D7: k in [0.54, 1.47]
H_DRIFT_TOL = 1e-4
H_TOL = 1e-10
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
        th = -th                      # D11: tlo lustrzane = negacja fazy
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


print("=== Phase 1: kalibracja + background_check "
      "(LOCK: Phase0_balance.md, D1/D7/D8) ===")
print(f"N={N} L={L} dx={dx} kappa={kappa} dt={dt}")
print()

# ---------------- G2: dyspersja kanalu amplitudowego ----------------
print("--- G2: dyspersja na tle jednorodnym s* (D7) ---")
print("    m     k       omega_meas   omega_teor   blad%   omega_siec   "
      "H_dryft(sek.)  H_osc")
PROBE = (N // 2, 10)                    # x=0: cos(k*0)=1 dla kazdego k
disp_ok, g1_ok, det_ok = True, True, True
errs, drifts = [], []
for m in M_LIST:
    k = TWO_PI * m / L
    pert = 1e-3 * np.cos(k * X)
    psi0 = (s_star + pert).astype(complex)
    psi_prev0 = psi0.copy()             # u_tau(0)=0 (fala stojaca)

    def run_g2():
        psi_prev, psi = psi_prev0.copy(), psi0.copy()
        series, Hs, taus = [], [], []
        nf = False
        for stp in range(1, STEPS + 1):
            psi_next = 2.0 * psi - psi_prev + dt * dt * rhs(psi)
            if stp % 10 == 0:
                if not np.all(np.isfinite(psi_next)):
                    nf = True
                    break
                psit = (psi_next - psi_prev) / (2.0 * dt)   # D1: centralna
                series.append(float(np.real(psi[PROBE]) - s_star))
                Hs.append(energy(psi, psit))
                taus.append(stp * dt)
            psi_prev, psi = psi, psi_next
        return np.array(series), np.array(Hs), np.array(taus), nf, psi

    s1, Hs, taus, nf, fin1 = run_g2()
    if m == M_LIST[0]:
        s2, _, _, _, fin2 = run_g2()
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

# ---------------- background_check ----------------
print("--- background_check: para (+1,-1), ansatz zlozony ---")
psi = make_pair()
th0 = np.angle(psi)
gxa = np.abs(wrap(np.roll(th0, -1, 0) - th0))
gya = np.abs(wrap(np.roll(th0, -1, 1) - th0))
far = (np.hypot((X - V1[0] + L/2) % L - L/2, (Y - V1[1] + L/2) % L - L/2) > 4) & \
      (np.hypot((X - V2[0] + L/2) % L - L/2, (Y - V2[1] + L/2) % L - L/2) > 4)
print(f"  diagnostyka szwu: max |dtheta| na linkach poza rdzeniami = "
      f"{max(gxa[far].max(), gya[far].max()):.4f} "
      f"(wadliwy ansatz: pi; artefakt C ~1.6 znany, relaksuje sie)")
pos0 = core_positions(psi)
print(f"  kret ansatzu: n_tot = {pos0[1]:.0f}; pozycje: "
      f"+1 {tuple(round(v, 2) for v in pos0[0][+1])}, "
      f"-1 {tuple(round(v, 2) for v in pos0[0][-1])}")
H_prev, H_inc = None, 0.0
pos1000 = None
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
    if stp == 1000:
        pos1000 = core_positions(psi)
psi_bg = psi
pos2000 = core_positions(psi_bg)
res = float(np.max(np.abs(rhs(psi_bg))))
drift_core = max(np.hypot(pos2000[0][s][0] - pos1000[0][s][0],
                          pos2000[0][s][1] - pos1000[0][s][1])
                 for s in (+1, -1))
# determinizm bitowy relaksacji (D8)
psi_d = make_pair()
for stp in range(1, 2001):
    psi_d = psi_d + dt_flow * rhs(psi_d)
det_bg = bool(np.array_equal(psi_bg, psi_d))
print(f"  po relaksacji 2000 krokow (dt={dt_flow}):")
print(f"    pozycje rdzeni: +1 {tuple(round(v, 3) for v in pos2000[0][+1])}, "
      f"-1 {tuple(round(v, 3) for v in pos2000[0][-1])}; "
      f"n_tot = {pos2000[1]:.0f}")
print(f"    |psi| daleko od rdzeni: min = "
      f"{float(np.abs(psi_bg)[far].min()):.4f} (s* = {s_star:.4f})")
print(f"    residuum max|RHS| = {res:.3e} (tlo QUASI-stacjonarne — "
      f"korekta przez b_eff/psi_ref, D3/D5; to NIE jest kryterium)")
print(f"    dryf rdzeni (krok 1000 -> 2000, tau=20) = {drift_core:.4f} "
      f"~ {drift_core/20:.5f}/tau (cytowane 0.015/tau)")
print(f"    determinizm bitowy relaksacji: {det_bg}; "
      f"H_mono: {H_inc <= H_TOL}")
print()

# ---------------- EWALUACJA ----------------
print("=== EWALUACJA (LOCKED) ===")
g1_part = g1_ok and det_ok and det_bg
print(f"G1 (czesc Phase1): determinizm bitowy (fala: {det_ok}, "
      f"relaksacja: {det_bg}); brak NaN; max dryft sekularny H = "
      f"{max(drifts):.2e} <= {H_DRIFT_TOL}: {g1_ok}; kret tla n_tot=0: "
      f"{pos2000[1] == 0}  -> {'PASS' if g1_part and pos2000[1] == 0 else 'FAIL'}")
g2 = disp_ok
print(f"G2 (kalibracja dyspersji): max |blad| = {max(errs):.2f}% <= 5%: "
      f"{disp_ok}  -> {'PASS' if g2 else 'FAIL'}")
print()
print(f"PODSUMOWANIE PHASE 1: G1(czesc)={'PASS' if g1_part and pos2000[1] == 0 else 'FAIL'} "
      f"G2={'PASS' if g2 else 'FAIL'}")
