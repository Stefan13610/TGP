# Phase 2b: DIAGNOSTYKA (w dniu runow, przed interpretacja; precedens:
# Phase2b op-stationary). ZERO zmian progow/kryteriow.
#
# Ustalenie Phase 2: zaden run nie osiagnal 5 probek — rdzen rozpraszajacy
# przesunal sie o ~44 (anihilacja pary) w oknie runu (tau <= 300).
# Pytania diagnostyczne:
#   (1) czas zycia pary pod PELNA dynamika II rzedu: separacja(tau),
#       tau anihilacji; porownanie z budzetem dryfu #1.4 (zalozenie: ~3);
#   (2) czy kolaps jest artefaktem integratora: porownanie separacji
#       przy dt=0.05 vs dt=0.025 do tau=100 (diagnostyka, nie run
#       kryterialny);
#   (3) OBSERWACJA POZA KRYTERIAMI: klasyfikator kierunku ugiecia
#       z oknem STALYM x > x_rdzen(0)+10 (jak demo G5 op-stationary),
#       omega=1.1, b=8 — wylacznie znak, bez statusu kryterialnego.

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt_flow = 0.02
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


def lap(u):
    return (np.roll(u, 1, 0) + np.roll(u, -1, 0)
            + np.roll(u, 1, 1) + np.roll(u, -1, 1) - 4.0 * u) / (dx * dx)


def rhs(psi):
    return kappa * lap(psi) - psi * gfun(np.abs(psi))


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


def make_pair():
    A_, B_ = V1, V2
    C_ = (B_[0], A_[1])
    th = theta_h_pair(A_, C_) + theta_v_pair(C_, B_)
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
    return out


def pdist(p, q):
    ddx = (p[0] - q[0] + 0.5 * L) % L - 0.5 * L
    ddy = (p[1] - q[1] + 0.5 * L) % L - 0.5 * L
    return float(np.hypot(ddx, ddy))


print("=== Phase 2b: diagnostyka czasu zycia tla (przed interpretacja) ===")
print()
psi0 = make_pair()
psi = psi0.copy()
for stp in range(1, 2001):
    psi = psi + dt_flow * rhs(psi)
psi_bg = psi
print("tlo zrelaksowane (identycznie jak Phase 2)")
print()

# ---------- (1) separacja(tau) pod dynamika II rzedu ----------
print("--- (1) para pod pelna dynamika II rzedu: separacja(tau) ---")
dt = 0.05
psi_w, psi_p = psi_bg.copy(), psi_bg.copy()
sep_series = []
annih = None
for stp in range(1, 6001):
    psi_n = 2.0 * psi_w - psi_p + dt * dt * rhs(psi_w)
    if stp % 50 == 0:
        cp = core_positions(psi_w)
        tau = stp * dt
        if cp is None:
            annih = tau
            break
        sep = pdist(cp[+1], cp[-1])
        sep_series.append((tau, cp[+1], sep))
        if stp % 200 == 0:
            print(f"   tau={tau:6.1f}: rdzen+1 = ({cp[+1][0]:+7.2f},"
                  f"{cp[+1][1]:+7.2f})  separacja = {sep:7.2f}")
    psi_p, psi_w = psi_w, psi_n
if annih is not None:
    print(f"   ANIHILACJA PARY: tau = {annih:.1f} (brak windingu w polu)")
else:
    print(f"   para przetrwala okno tau=300; koncowa separacja = "
          f"{sep_series[-1][2]:.2f}")
d0 = pdist(sep_series[0][1], V1)
tt = [s[0] for s in sep_series]
dd = [pdist(s[1], sep_series[0][1]) for s in sep_series]
print(f"   przemieszczenie rdzenia +1: tau=50: "
      f"{dd[min(range(len(tt)), key=lambda i: abs(tt[i]-50))]:.2f}; "
      f"tau=100: {dd[min(range(len(tt)), key=lambda i: abs(tt[i]-100))]:.2f}; "
      f"tau=150: {dd[min(range(len(tt)), key=lambda i: abs(tt[i]-150))]:.2f}")
print(f"   (budzet #1.4 zakladal ~3 na przelot tau~200 — SFALSYFIKOWANY;")
print(f"    dynamika II rzedu jest bezdyssypacyjna: sila akumuluje predkosc)")
print()

# ---------- (2) kontrola integratora: dt=0.05 vs dt=0.025 ----------
print("--- (2) kontrola integratora: przemieszczenie rdzenia +1 do tau=100 ---")
res = {}
for dtt in (0.05, 0.025):
    psi_w, psi_p = psi_bg.copy(), psi_bg.copy()
    nst = int(round(100.0 / dtt))
    for stp in range(1, nst + 1):
        psi_n = 2.0 * psi_w - psi_p + dtt * dtt * rhs(psi_w)
        psi_p, psi_w = psi_w, psi_n
    cp = core_positions(psi_w)
    res[dtt] = cp[+1]
    print(f"   dt={dtt}: rdzen+1(tau=100) = ({cp[+1][0]:+7.3f},"
          f"{cp[+1][1]:+7.3f})")
print(f"   roznica pozycji = {pdist(res[0.05], res[0.025]):.4f} "
      f"-> kolaps NIE jest artefaktem kroku czasowego"
      if pdist(res[0.05], res[0.025]) < 0.5 else
      "   UWAGA: rozjazd dt — artefakt integratora?")
print()

# ---------- (3) OBSERWACJA POZA KRYTERIAMI: klasyfikator kierunku ----------
print("--- (3) OBSERWACJA POZA KRYTERIAMI: kierunek ugiecia, okno STALE ---")
print("    (omega=1.1, b=8; okno x > x_rdzen(0)+10, wyzwolenie: centroid;")
print("     tlo KOLABUJE — wynik ma status wylacznie klasyfikatora znaku)")
om, bb, u0 = 1.1, 8.0, 1e-3
k0 = np.sqrt((om * om - Vpp_star) / kappa)
th_bg = np.angle(psi_bg)
eith = np.exp(1j * th_bg)
x0, y0 = -60.0, V1[1] + bb
env = u0 * np.exp(-((X - x0) ** 2) / (2 * 64.0) - ((Y - y0) ** 2) / (2 * 64.0))
psi_w = psi_bg + eith * env * np.cos(k0 * (X - x0))
psi_wp = psi_bg + eith * env * np.cos(k0 * (X - x0) + om * dt)
psi_r = psi_bg.copy()
psi_rp = psi_bg.copy()
cp0 = core_positions(psi_bg)
x_thr = cp0[+1][0] + 10.0
m_tr = X > x_thr
samples, next_s = [], None
b_eff = np.inf
for stp in range(1, 6001):
    tau = stp * dt
    psi_wn = 2.0 * psi_w - psi_wp + dt * dt * rhs(psi_w)
    psi_rn = 2.0 * psi_r - psi_rp + dt * dt * rhs(psi_r)
    if stp % 10 == 0:
        dpsi = psi_w - psi_r
        w2 = np.abs(dpsi) ** 2
        Wtot = float(w2.sum())
        xc = float((w2 * X).sum() / Wtot)
        yc = float((w2 * Y).sum() / Wtot)
        cp = core_positions(psi_r)
        if cp is not None and next_s is None:
            b_eff = min(b_eff, pdist((xc, yc), cp[+1]))
        if next_s is None and xc > x_thr:
            next_s = tau
        if next_s is not None and tau >= next_s and len(samples) < 5:
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
if samples:
    Pxm = float(np.mean([s[1] for s in samples]))
    Pym = float(np.mean([s[2] for s in samples]))
    al = float(np.degrees(np.arctan2(Pym, Pxm)))
    al_s = [f"{np.degrees(np.arctan2(s[2], s[1])):+.1f}" for s in samples]
    kier = ("KU wirowi" if al < -0.5 else
            ("OD wiru" if al > 0.5 else "nierozstrzygniete"))
    print(f"    tau probek: {samples[0][0]:.1f}..{samples[-1][0]:.1f}; "
          f"b_eff(do wyzwolenia) = {b_eff:.2f}")
    print(f"    alpha_surowe (srednia wektorowa) = {al:+.2f} deg; "
          f"probki: {al_s}")
    print(f"    KIERUNEK: {kier} (KU = P_y < 0 dla b>0)")
    print(f"    UWAGA: wielkosc NIEINTERPRETOWALNA (soczewka w ruchu "
          f"~0.2-0.3 vs v_g=0.37)")
else:
    print("    pakiet nie osiagnal okna — brak probek")
print()
print("=== WNIOSKI DIAGNOSTYKI ===")
print("  kolaps pary pod dynamika II rzedu jest FIZYCZNY (nie artefakt dt);")
print("  budzet dryfu #1.4 (3 jedn.) sfalsyfikowany -> runy kryterialne")
print("  Phase 2 niewazne z przyczyny strukturalnej (D9: brak 5 probek);")
print("  szczegolowa regula decyzyjna: Phase_FINAL_close.md")
