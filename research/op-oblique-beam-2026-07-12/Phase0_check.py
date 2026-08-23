# Phase 0 check: rachunki PRZED kodem produkcyjnym
# (LOCK: Phase0_balance.md proby #5, sekcje 1-2 + D2, D4, D6, D12, D13)
#   1) k0(omega), v_g, lambda — kanal amplitudowy + TABELA PROJEKTOWA
#      par kryterialnych +/-b: b_nom vs 1.5*lambda + 0.9 (D12)
#   2) GEOMETRIA PASA SKOSNEGO: obrazy soczewek w (xi, eta) wzgledem A
#      + odleglosci prostopadle od linii wiazki dla kazdego +/-b
#      (kryterialne i obserwacyjne) — weryfikacja rownania #2
#   3) dyspersja SKOSNA analitycznie: k = k0*(2,1)/sqrt(5), czynnik
#      anizotropii 17/25; mody G3 (2q,q) q in {10,13,20,25,30}
#   4) tlo rzeczywiste: szachownica 4 wirow L=256 (ansatz jak proba #3
#      BEZ ZMIAN), relaksacja 2000; szew, detektor 2+2, residuum
#      (oczekiwane 9.6e-07 — proba #3, ten sam kod)
#   5) ray tracing WIAZACY w ramce skosnej (D6): alpha_pred(+/-b)
#      OSOBNO dla kazdego punktu, plaszczyzny xi_eval in
#      {10, 40, 70, 100, 130}; Delta_eik per para (plaszczyzna po
#      plaszczyznie) + rozrzut plaszczyzn; porownanie z orientacyjna
#      formula ogonowa (sekcja 1.3)
#   6) budzet czasu (sekcja 1.4)

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt_flow = 0.02
N, L = 512, 256.0
dx = L / N
TWO_PI = 2.0 * np.pi
S5 = np.sqrt(5.0)

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)
Phi_star = s_star ** 2
Vpp_star = a - 2.0 * b_par * s_star + 3.0 * c * s_star ** 2
Vppp_star = -2.0 * b_par + 6.0 * c * s_star
xi2 = kappa / Vpp_star
m_TV = np.sqrt(Vpp_star)
A_t = Vppp_star * s_star * xi2

coords = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
Z = X + 1j * Y

# szachownica L=256 — DOKLADNIE jak proba #3
VORTS = [((-64.0, -64.0), +1), ((+64.0, -64.0), -1),
         ((-64.0, +64.0), -1), ((+64.0, +64.0), +1)]
NOM = {+1: [(-64.0, -64.0), (+64.0, +64.0)],
       -1: [(+64.0, -64.0), (-64.0, +64.0)]}
A_NOM = (-64.0, -64.0)

# wersory wiazki skosnej (D2)
UVEC = np.array([2.0, 1.0]) / S5
NVEC = np.array([-1.0, 2.0]) / S5
XI_START = -36.0

# pary kryterialne (ZAMKNIETE) i obserwacje (D12)
CRIT_PAIRS = [(1.1, 14.0), (1.1, 20.0), (1.5, 8.0), (1.5, 12.0)]
OBS_PAIRS = [(1.1, 24.0)]
# formula ogonowa orientacyjna (sekcja 1.3)
TAIL = {1.1: 698.0, 1.5: 169.0}


def gfun(m):
    return a - b_par * np.abs(m) + c * m * m


def Vpp(m):
    return a - 2.0 * b_par * np.abs(m) + 3.0 * c * m * m


def lap(u):
    return (np.roll(u, 1, 0) + np.roll(u, -1, 0)
            + np.roll(u, 1, 1) + np.roll(u, -1, 1) - 4.0 * u) / (dx * dx)


def rhs(psi):
    return kappa * lap(psi) - psi * gfun(np.abs(psi))


def wrap(x):
    return (x + np.pi) % TWO_PI - np.pi


def pdist(p, q):
    ddx = (p[0] - q[0] + 0.5 * L) % L - 0.5 * L
    ddy = (p[1] - q[1] + 0.5 * L) % L - 0.5 * L
    return float(np.hypot(ddx, ddy))


def xi_eta(x, y, A):
    """Ramka skosna D4: skladowe wrapowane do [-L/2, L/2) wzgledem A."""
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


def core_positions(psi):
    """Detektor 2+2 (D14) — identycznie jak proba #3."""
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


print("=== Phase 0 check: op-oblique-beam (LOCK: Phase0_balance.md, "
      "proba #5) ===")
print()
print("--- 1. kanal amplitudowy + TABELA PROJEKTOWA par +/-b "
      "(rownanie #1, D12) ---")
print(f"  s* = {s_star:.6f}  Phi* = {Phi_star:.6f}  V''(s*) = {Vpp_star:.6f}")
print(f"  m_TV = {m_TV:.4f}  xi^2 = {xi2:.4f}  A_t = {A_t:.4f}")
lam = {}
for om in (1.1, 1.5):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    lam[om] = TWO_PI / k0
    print(f"  omega={om:.1f}:  k0 = {k0:.4f}   lambda = {TWO_PI / k0:6.2f}   "
          f"v_g = {kappa * k0 / om:.4f}")
print()
print("  TABELA PROJEKTOWA (kryterium: |b_nom| >= 1.5*lambda + 0.9):")
for om, bb in CRIT_PAIRS + OBS_PAIRS:
    thr = 1.5 * lam[om] + 0.9
    tag = "KRYT para" if (om, bb) in CRIT_PAIRS else "obs para "
    print(f"    [{tag}] omega={om:.1f} b=+/-{bb:4.0f}: prog = {thr:5.2f}  "
          f"|b| - prog = {bb - thr:+5.2f}  |b|/lambda = {bb / lam[om]:.2f}"
          f"  -> {'W DOMENIE' if bb >= thr else 'POZA'}")
print()

# --- 2. geometria pasa skosnego ---
print("--- 2. GEOMETRIA PASA SKOSNEGO (rownanie #2): obrazy soczewek "
      "w (xi, eta) wzgledem A_nom ---")
print("    u = (2,1)/sqrt(5); n = (-1,2)/sqrt(5); phi = atan(1/2) = "
      f"{np.degrees(np.arctan(0.5)):.3f} deg")
print("    obrazy periodyczne 4 rdzeni z -50 < xi < 300 (pomijajac A):")
images = []
for (vx, vy), q in VORTS:
    for ki in range(-2, 3):
        for kj in range(-2, 3):
            px_, py_ = vx + ki * L, vy + kj * L
            ddx = px_ - A_NOM[0]
            ddy = py_ - A_NOM[1]
            xi_v = (2.0 * ddx + ddy) / S5
            eta_v = (-ddx + 2.0 * ddy) / S5
            if -50.0 < xi_v < 300.0 and abs(eta_v) < 180.0 \
                    and not (abs(ddx) < 1 and abs(ddy) < 1):
                images.append((xi_v, eta_v, q, (ddx / 128.0, ddy / 128.0)))
images.sort()
for xi_v, eta_v, q, (ii, jj) in images:
    in_strip = "W PASIE {10<xi<130}" if 10.0 < xi_v < 130.0 else ""
    print(f"      (i,j)=({ii:+.1f},{jj:+.1f}) [q={q:+d}]: "
          f"xi = {xi_v:7.2f}  eta = {eta_v:+8.2f}  {in_strip}")
print()
print("    odleglosci prostopadle |eta_soczewki - b| od linii wiazki "
      "(soczewki pasa i najblizsze):")
LENS = [(57.24, +114.49, -1, "(0,1)"), (114.49, -57.24, -1, "(1,0)"),
        (171.73, +57.24, +1, "(1,1)"), (286.22, 0.00, +1, "(2,1)")]
for om, bb in CRIT_PAIRS + OBS_PAIRS:
    for sgn in (+1.0, -1.0):
        b = sgn * bb
        row = f"      omega={om:.1f} b={b:+5.0f}: "
        row += "  ".join(f"{name}: {abs(eta_v - b):6.2f}"
                         for xi_v, eta_v, q, name in LENS)
        dmin = min(abs(eta_v - b) for xi_v, eta_v, q, name in LENS[:2])
        ok = dmin >= 33.0
        row += f"   min(pas) = {dmin:.2f} >= 33: {ok}"
        print(row)
print("    (obraz (2,1) lezy NA osi wiazki przy xi=286.2 — poza oknem "
      "130 i poza 5. probka; D12)")
print()

# --- 3. dyspersja skosna analitycznie ---
print("--- 3. dyspersja sieciowa SKOSNA analitycznie (dx=0.5; kierunek "
      "(2,1)/sqrt(5), anizotropia (kx^4+ky^4)/k^4 = 17/25) ---")
for om in (1.1, 1.5):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    kx, ky = 2.0 * k0 / S5, k0 / S5
    k2_lat = ((2 - 2 * np.cos(kx * dx)) + (2 - 2 * np.cos(ky * dx))) / (dx * dx)
    om_lat = np.sqrt(kappa * k2_lat + Vpp_star)
    err = 100.0 * (om_lat - om) / om
    # dla porownania: axialnie (proba #3)
    k2_ax = (2 - 2 * np.cos(k0 * dx)) / (dx * dx)
    om_ax = np.sqrt(kappa * k2_ax + Vpp_star)
    err_ax = 100.0 * (om_ax - om) / om
    print(f"    omega={om:.1f} (k0={k0:.4f}): skosnie blad = {err:+.2f}%  "
          f"(axialnie, proba #3: {err_ax:+.2f}%)")
print("    mody G3 (D7): k = (2 pi/L)(2q, q):")
for q in (10, 13, 20, 25, 30):
    kx = TWO_PI * 2 * q / L
    ky = TWO_PI * q / L
    kk = np.hypot(kx, ky)
    om_t = np.sqrt(kappa * kk * kk + Vpp_star)
    k2_lat = ((2 - 2 * np.cos(kx * dx)) + (2 - 2 * np.cos(ky * dx))) / (dx * dx)
    om_lat = np.sqrt(kappa * k2_lat + Vpp_star)
    err = 100.0 * (om_lat - om_t) / om_t
    print(f"      q={q:2d}: |k| = {kk:.4f}  omega_teor = {om_t:.5f}  "
          f"omega_siec(analit.) = {om_lat:.5f}  blad = {err:+.2f}%")
print("    prog G3: 5% (MIERZONE w Phase 1, nie zakladane)")
print()

# --- 4. tlo rzeczywiste ---
print("--- 4. tlo: szachownica 4 wirow L=256 (BEZ ZMIAN vs proba #3), "
      "relaksacja 2000 ---")
psi = make_lattice()
th0 = np.angle(psi)
gxa = np.abs(wrap(np.roll(th0, -1, 0) - th0))
gya = np.abs(wrap(np.roll(th0, -1, 1) - th0))
far = np.ones((N, N), bool)
for (x0, y0), _n in VORTS:
    far &= np.hypot((X - x0 + L/2) % L - L/2, (Y - y0 + L/2) % L - L/2) > 4
print(f"  diagnostyka szwu (ansatz): max |dtheta| poza rdzeniami = "
      f"{max(gxa[far].max(), gya[far].max()):.4f} rad (proba #3: 0.1244)")
for stp in range(1, 2001):
    psi = psi + dt_flow * rhs(psi)
psi_bg = psi
f_bg = np.abs(psi_bg)
res = float(np.max(np.abs(rhs(psi_bg))))
cpR = core_positions(psi_bg)
print(f"  po relaksacji 2000 krokow (dt={dt_flow}):")
print(f"    residuum max|RHS| = {res:.3e} (proba #3: 9.648e-07)")
if cpR is not None:
    print(f"    pozycje rdzeni: "
          + "; ".join(f"{s:+d}#{i}: ({p[0]:+.3f},{p[1]:+.3f})"
                      for (s, i), p in sorted(cpR[0].items())))
    print(f"    n_tot = {cpR[1]:.0f}")
A_REL = cpR[0][(+1, 0)]
print(f"    rdzen A (rozpraszajacy) po relaksacji: "
      f"({A_REL[0]:+.3f}, {A_REL[1]:+.3f})")
print()

# --- 5. ray tracing WIAZACY w ramce skosnej (D6) ---
print("--- 5. ray tracing WIAZACY (D6): alpha_pred [deg, KU A] czytany "
      "w plaszczyznach xi_eval in {10, 40, 70, 100, 130} ---")
print("    start: r = A + (-36)u + b*n, kierunek u; przechwyt r<0.6;")
print("    Delta_eik(para) = alpha_pred(+b) - alpha_pred(-b) "
      "plaszczyzna po plaszczyznie; rozrzut = max-min po plaszczyznach")

n2_grid, gx_grid, gy_grid = {}, {}, {}
for om in (1.1, 1.5):
    n2 = (om * om - Vpp(f_bg)) / (om * om - Vpp_star)
    n2_grid[om] = n2
    gx_grid[om] = (np.roll(n2, -1, 0) - np.roll(n2, 1, 0)) / (2 * dx)
    gy_grid[om] = (np.roll(n2, -1, 1) - np.roll(n2, 1, 1)) / (2 * dx)


def bilin_f(G, x, y):
    fx = ((x - coords[0]) / dx) % N
    fy = ((y - coords[0]) / dx) % N
    i0, j0 = int(fx) % N, int(fy) % N
    i1, j1 = (i0 + 1) % N, (j0 + 1) % N
    tx, ty = fx - int(fx), fy - int(fy)
    return ((1 - tx) * (1 - ty) * G[i0, j0] + tx * (1 - ty) * G[i1, j0]
            + (1 - tx) * ty * G[i0, j1] + tx * ty * G[i1, j1])


PLANES0 = (10.0, 40.0, 70.0, 100.0, 130.0)


def trace_oblique(om, b_signed, planes=PLANES0, A=A_REL):
    """Ray w ramce skosnej: dict {xi_plane: alpha_ku_deg} albo None
    (przechwyt r<0.6 od ktoregokolwiek rdzenia)."""
    n2, gxg, gyg = n2_grid[om], gx_grid[om], gy_grid[om]
    r0 = np.array(A) + XI_START * UVEC + b_signed * NVEC
    n_start = np.sqrt(bilin_f(n2, r0[0], r0[1]))
    p0 = n_start * UVEC
    ds = 0.02
    sgn = -1.0 if b_signed > 0 else +1.0

    def F(st):
        return np.array([st[2], st[3],
                         0.5 * bilin_f(gxg, st[0], st[1]),
                         0.5 * bilin_f(gyg, st[0], st[1])])

    st = np.array([r0[0], r0[1], p0[0], p0[1]])
    out = {}
    todo = sorted(planes)
    xi_old, _ = xi_eta(st[0], st[1], A)
    for _ in range(1000000):
        k1 = F(st)
        k2 = F(st + 0.5 * ds * k1)
        k3 = F(st + 0.5 * ds * k2)
        k4 = F(st + ds * k3)
        st_new = st + ds / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        xi_new, _ = xi_eta(st_new[0], st_new[1], A)
        while todo and xi_old <= todo[0] < xi_new:
            f_ = (todo[0] - xi_old) / (xi_new - xi_old)
            pxi = st[2] + f_ * (st_new[2] - st[2])
            pyi = st[3] + f_ * (st_new[3] - st[3])
            pu = (2.0 * pxi + pyi) / S5
            pn = (-pxi + 2.0 * pyi) / S5
            out[todo.pop(0)] = sgn * float(np.degrees(np.arctan2(pn, pu)))
        st = st_new
        xi_old = xi_new
        if not todo:
            return out
        for (x0, y0), _n in VORTS:
            if np.hypot((st[0] - x0 + L/2) % L - L/2,
                        (st[1] - y0 + L/2) % L - L/2) < 0.6:
                return None
    return None


def tail_orient(om, bb):
    """Formula ogonowa orientacyjna (sekcja 1.3): rozpraszacz +
    korekta sasiadow pasa (2x polowa przelotu, orientacyjnie)."""
    t = TAIL[om]
    base = t / bb**2
    d_plus = (t / (57.24 - (-bb))**2 - t / (114.49 - bb)**2)
    return base, d_plus


for om, bb in CRIT_PAIRS + OBS_PAIRS:
    crit_t = "[KRYT]" if (om, bb) in CRIT_PAIRS else "[obs] "
    res_p = trace_oblique(om, +bb)
    res_m = trace_oblique(om, -bb)
    print(f"    {crit_t} omega={om:.1f} para b=+/-{bb:.0f}:")
    for lbl, res_t in (("+b", res_p), ("-b", res_m)):
        if res_t is None:
            print(f"      {lbl}: PRZECHWYT (r<0.6)")
            continue
        row = "  ".join(f"xi={xp:3.0f}: {res_t[xp]:+7.3f}"
                        for xp in PLANES0 if xp in res_t)
        spread = max(res_t.values()) - min(res_t.values())
        print(f"      {lbl}: {row}   (rozrzut {spread:.3f})")
    if res_p is not None and res_m is not None:
        deltas = {xp: res_p[xp] - res_m[xp] for xp in PLANES0}
        row = "  ".join(f"xi={xp:3.0f}: {dv:+7.3f}"
                        for xp, dv in sorted(deltas.items()))
        spr = max(deltas.values()) - min(deltas.values())
        print(f"      Delta_eik(plaszczyzny): {row}")
        print(f"      Delta_eik rozrzut plaszczyzn = {spr:.3f} deg  "
              f"(odniesienie reguly wzmacniajacej G5b)")
    base, dor = tail_orient(om, bb)
    print(f"      [orient. sekcja 1.3: alpha ~ {base:+.2f}; "
          f"Delta ~ {2*dor:+.2f} (formula ogonowa, polowkowy przelot)]")
print()

# --- 6. budzet czasu ---
print("--- 6. budzet czasu (sekcja 1.4) ---")
for om in (1.1, 1.5):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    vg = kappa * k0 / om
    T_trig = (10.0 - XI_START) / vg
    T_last = T_trig + 16.0
    T_130 = (130.0 - XI_START) / vg
    print(f"  omega={om:.1f}: v_g={vg:.3f}; wyzwolenie (xi=+10) ~tau="
          f"{T_trig:5.1f}; 5. probka ~tau={T_last:5.1f}; "
          f"koniec okna (xi=130) ~tau={T_130:5.1f}")
print("  limit runow falowych tau=300 (6000 krokow) pokrywa 5. probke")
print("  obu serii z zapasem; bramka G2: okno tau=600, kryterium")
print("  tau<=450 >= 3x przelot (~140 przy omega=1.1)")
print("  koszt: jak proba #3 (bramka ~15 min; 11 runow Phase2 +")
print("  3 runy Phase3 ~ 1-1.5 h)")
print()
print("=== WNIOSKI PHASE 0 ===")
print("  wszystkie 4 pary kryterialne W DOMENIE eikonalu (|b| >=")
print("  1.5*lambda + 0.9) i w odleglosci >= 33 od soczewek pasa;")
print("  kierunek (2,1) poza wszystkimi lustrami D4 (obrazy soczewek")
print("  asymetryczne: eta=-57.2 vs +114.5) -> Delta_eik != 0,")
print("  policzalna, WIAZACO z ray tracingu OSOBNO dla +b i -b wraz")
print("  z rozrzutem plaszczyzn; dyspersja skosna w progu 5%")
print("  analitycznie (pomiar: G3 na modach (2q,q));")
print("  Phase 1 = bramka G2 (tau=600) PRZED jakimkolwiek pakietem")
