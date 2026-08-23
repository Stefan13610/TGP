# Phase 0 check: rachunki PRZED kodem produkcyjnym
# (LOCK: Phase0_balance.md, sekcje 1-2 + D6, D13, D14)
#   1) k0(omega), v_g, lambda — kanal amplitudowy (rownanie #1)
#   2) BILANS SIL na siatce szachownicowej (rownanie #2): suma
#      punktowo-wirowa po sasiadach i obrazach periodycznych,
#      powloki K=1..8; KONTRAST: para (+1,-1) proby #1
#   3) tlo rzeczywiste: szachownica 4 wirow (ansatz = suma dwoch par
#      poziomych theta_h_pair, BEZ transpozycji), relaksacja 2000;
#      diagnostyka szwu, detektor 2+2 (D14), residuum
#   4) n^2 na zrelaksowanym tle siatki vs ogon algebraiczny A_t/r^2
#   5) ray tracing na RZECZYWISTYM tle siatki: alpha_pred(b, omega)
#      z CZULOSCIA plaszczyzny ewaluacji (D6): x_eval in
#      {-22, -12, 0, +12, +22}; porownanie z tabela ORIENTACYJNA
#      (tlo pary, proba #1)
#   6) budzet czasu (sekcja 1.5): przelot vs okno bramki G2

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt_flow = 0.02
N, L = 256, 128.0
dx = L / N
TWO_PI = 2.0 * np.pi

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

# szachownica (sekcja 2): +1 (-32,-32), -1 (+32,-32), -1 (-32,+32), +1 (+32,+32)
VORTS = [((-32.0, -32.0), +1), ((+32.0, -32.0), -1),
         ((-32.0, +32.0), -1), ((+32.0, +32.0), +1)]
NOM = {+1: [(-32.0, -32.0), (+32.0, +32.0)],
       -1: [(+32.0, -32.0), (-32.0, +32.0)]}
V_SC = (-32.0, -32.0)   # rdzen rozpraszajacy (+1)

# tabela ORIENTACYJNA (tlo PARY, proba #1 Phase0_output.txt)
ORIENT = {(1.1, 8): +14.577, (1.1, 12): +4.993, (1.1, 16): +2.651,
          (1.3, 8): +5.808, (1.3, 12): +2.054, (1.3, 16): +1.090,
          (1.1, 6): +39.081, (1.0, 8): +50.330}


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


def theta_h_pair(z1, z2):
    """Para wspolliniowa pozioma: +1 w z1, -1 w z2 (konstrukcja D1 z B1)."""
    th = np.zeros_like(X)
    for (x0, y0), n in ((z1, +1), (z2, -1)):
        for k in range(-2, 3):
            th = th + n * np.angle(np.sinh(np.pi * (Z - (x0 + k * L) - 1j * y0) / L))
    return th


def make_lattice(mirror=False):
    """Szachownica: rzad y=-32: +1/-1; rzad y=+32: +1 w (+32), -1 w (-32)."""
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
    """Detektor 2+2 (D14): klastrowanie komorek windingu, centroid
    deficytu, przypisanie do nominalow. Zwraca (dict, n_tot) albo None."""
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


print("=== Phase 0 check: op-lattice-background (LOCK: Phase0_balance.md) ===")
print()
print("--- 1. kanal amplitudowy: k0, v_g, lambda (rownanie #1) ---")
print(f"  s* = {s_star:.6f}  Phi* = {Phi_star:.6f}  V''(s*) = {Vpp_star:.6f}")
print(f"  m_TV = {m_TV:.4f}  xi^2 = {xi2:.4f}  A_t = {A_t:.4f}")
for om in (1.0, 1.1, 1.3):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    print(f"  omega={om:.1f}:  k0 = {k0:.4f}   lambda = {TWO_PI / k0:6.2f}   "
          f"v_g = {kappa * k0 / om:.4f}")
print(f"  UWAGA dyfrakcyjna (D12): lambda(1.1)=7.72 ~ b — rezim falowy,")
print(f"  stad szerokie pasmo G5 [0.5, 2.0]")
print()

# --- 2. bilans sil na siatce (rownanie #2) ---
print("--- 2. bilans sil (punktowo-wirowy, obrazy periodyczne) ---")
print("    F_i = sum_j q_i q_j (r_i - r_j)/|r_i - r_j|^2; powloki obrazow")
print("    K = max(|kx|,|ky|); suma czesciowa po kazdej powloce")


def force_on(idx, vortlist, K):
    (xi_, yi_), qi = vortlist[idx]
    F = np.zeros(2)
    shells = []
    for k in range(0, K + 1):
        Fs = np.zeros(2)
        for j, ((xj, yj), qj) in enumerate(vortlist):
            for kx in range(-k, k + 1):
                for ky in range(-k, k + 1):
                    if max(abs(kx), abs(ky)) != k:
                        continue
                    if j == idx and kx == 0 and ky == 0:
                        continue
                    rx = xi_ - (xj + kx * L)
                    ry = yi_ - (yj + ky * L)
                    r2 = rx * rx + ry * ry
                    Fs += qi * qj * np.array([rx, ry]) / r2
        F += Fs
        shells.append((k, float(np.hypot(*F))))
    return F, shells


# SIATKA: na torusie obrazy tworza nieskonczona siec szachownicowa
# wokol KAZDEGO wiru: sasiedzi w (64i, 64j), ladunek (-1)^(i+j).
# Kasowanie jest DOKLADNE per powloka, jesli powloki grupowac
# SYMETRYCZNIE WOKOL WIRU (kazda powloka = unia orbit D4 o wspolnym
# ladunku). Obciecie pudelkowe po indeksie obrazow (nie-centrowane)
# zostawia czlon powierzchniowy ~K^-2 — artefakt obciecia, nie sila.
print("    SIATKA, wir +1 w (-32,-32): powloki centrowane na wirze")
print("    (sasiedzi w (64i,64j), ladunek (-1)^(i+j); powloka "
      "K = max(|i|,|j|))")
F_lat = np.zeros(2)
for k in range(1, 9):
    Fs = np.zeros(2)
    for i in range(-k, k + 1):
        for j in range(-k, k + 1):
            if max(abs(i), abs(j)) != k:
                continue
            qj = 1.0 if (i + j) % 2 == 0 else -1.0
            rx, ry = -64.0 * i, -64.0 * j
            r2 = rx * rx + ry * ry
            Fs += qj * np.array([rx, ry]) / r2
    F_lat += Fs
    print(f"      powloka K={k}: |F_powloki| = {float(np.hypot(*Fs)):.3e}"
          f"   |F_suma| = {float(np.hypot(*F_lat)):.3e}")
print(f"      -> F_netto = ({F_lat[0]:+.3e}, {F_lat[1]:+.3e})  "
      f"(kasowanie DOKLADNE z symetrii D4, do zaokraglen)")
F_box, sh_box = force_on(0, VORTS, 8)
print(f"      (kontrola: obciecie pudelkowe K=8 daje "
      f"|F| = {float(np.hypot(*F_box)):.3e} — czlon powierzchniowy "
      f"obciecia, maleje ~K^-2)")
PAIR = [((-32.0, -32.0), +1), ((+32.0, +32.0), -1)]
F_par, sh_par = force_on(0, PAIR, 8)
print("    KONTRAST — para proby #1, wir +1 w (-32,-32):")
print(f"      po powloce K=8: F = ({F_par[0]:+.3e}, {F_par[1]:+.3e}), "
      f"|F| = {float(np.hypot(*F_par)):.3e}  (!= 0)")
print("      pomiar proby #1 (Phase2b): przemieszczenie 14 przy tau=50,")
print("      anihilacja tau=107.5 — cytowane, bez powtorki")
print()

# --- 3. tlo rzeczywiste: szachownica ---
print("--- 3. tlo: szachownica 4 wirow, ansatz 2 pary poziome, "
      "relaksacja 2000 ---")
psi = make_lattice()
th0 = np.angle(psi)
gxa = np.abs(wrap(np.roll(th0, -1, 0) - th0))
gya = np.abs(wrap(np.roll(th0, -1, 1) - th0))
far = np.ones((N, N), bool)
for (x0, y0), _n in VORTS:
    far &= np.hypot((X - x0 + L/2) % L - L/2, (Y - y0 + L/2) % L - L/2) > 4
print(f"  diagnostyka szwu (ansatz): max |dtheta| na linkach poza rdzeniami"
      f" = {max(gxa[far].max(), gya[far].max()):.4f} rad"
      f"  (oczekiwane ogony dipolowe ~0.04; wadliwy ansatz: ~pi)")
cp0 = core_positions(psi)
if cp0 is None:
    print("  BLAD: detektor nie widzi 2+2 rdzeni w ansatzu")
else:
    print(f"  kret ansatzu: n_tot = {cp0[1]:.0f}; pozycje: "
          + "; ".join(f"{s:+d}#{i}: ({p[0]:+.2f},{p[1]:+.2f})"
                      for (s, i), p in sorted(cp0[0].items())))
for stp in range(1, 2001):
    psi = psi + dt_flow * rhs(psi)
psi_bg = psi
f_bg = np.abs(psi_bg)
res = float(np.max(np.abs(rhs(psi_bg))))
cpR = core_positions(psi_bg)
print(f"  po relaksacji 2000 krokow (dt={dt_flow}):")
print(f"    residuum max|RHS| = {res:.3e} "
      f"(proba #1, para: 8.107e-03)")
if cpR is not None:
    print(f"    pozycje rdzeni: "
          + "; ".join(f"{s:+d}#{i}: ({p[0]:+.3f},{p[1]:+.3f})"
                      for (s, i), p in sorted(cpR[0].items())))
    print(f"    n_tot = {cpR[1]:.0f}; max |poz - nominal| = "
          f"{max(pdist(p, NOM[s][i]) for (s, i), p in cpR[0].items()):.4f}")
print(f"    |psi| daleko od rdzeni: min = {float(f_bg[far].min()):.4f} "
      f"(s* = {s_star:.4f})")
print()

# --- 4. n^2 na zrelaksowanym tle siatki ---
print("--- 4. n^2 od rdzenia rozpraszajacego (srednia 4 kierunkow "
      "osiowych; kolumny omega 1.0/1.1/1.3) ---")


def bilin_f(G, x, y):
    fx = ((x - coords[0]) / dx) % N
    fy = ((y - coords[0]) / dx) % N
    i0, j0 = int(fx) % N, int(fy) % N
    i1, j1 = (i0 + 1) % N, (j0 + 1) % N
    tx, ty = fx - int(fx), fy - int(fy)
    return ((1 - tx) * (1 - ty) * G[i0, j0] + tx * (1 - ty) * G[i1, j0]
            + (1 - tx) * ty * G[i0, j1] + tx * ty * G[i1, j1])


sc = cpR[0][(+1, 0)] if cpR is not None else V_SC
for rq in (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0):
    fq = np.mean([bilin_f(f_bg, sc[0] + rq, sc[1]),
                  bilin_f(f_bg, sc[0] - rq, sc[1]),
                  bilin_f(f_bg, sc[0], sc[1] + rq),
                  bilin_f(f_bg, sc[0], sc[1] - rq)])
    row = f"   r={rq:5.1f}: "
    for om in (1.0, 1.1, 1.3):
        n2 = (om * om - Vpp(fq)) / (om * om - Vpp_star)
        row += f"  {n2:7.4f}"
    tail = 1.0 + A_t / (rq * rq * (1.21 - Vpp_star))
    row += f"   [ogon(1.1): {tail:.4f}]"
    print(row)
print("   (uwaga: siatka nie jest radialna — sasiedzi w +/-x, +/-y "
      "w odleglosci 64)")
print()

# --- 5. ray tracing na tle siatki + czulosc plaszczyzny ewaluacji ---
print("--- 5. ray tracing (D6): alpha_pred(b, omega) [deg], kat czytany")
print("    w plaszczyznach x_eval in {-22, -12, 0, +12, +22}")
print("    (nominalne krance okna: x1+10 = -22, x2-10 = +22);")
print("    w nawiasie [orient]: tabela orientacyjna tla PARY (proba #1)")

n2_grid, gx_grid, gy_grid = {}, {}, {}
for om in (1.0, 1.1, 1.3):
    n2 = (om * om - Vpp(f_bg)) / (om * om - Vpp_star)
    n2_grid[om] = n2
    gx_grid[om] = (np.roll(n2, -1, 0) - np.roll(n2, 1, 0)) / (2 * dx)
    gy_grid[om] = (np.roll(n2, -1, 1) - np.roll(n2, 1, 1)) / (2 * dx)

PLANES = (-22.0, -12.0, 0.0, +12.0, +22.0)


def trace_bg(om, b_signed, planes=PLANES, y_row=-32.0):
    """Ray na tle siatki; zwraca dict {x_plane: alpha_ku_deg} albo None
    (przechwyt: r<0.6 od ktoregokolwiek rdzenia)."""
    n2, gxg, gyg = n2_grid[om], gx_grid[om], gy_grid[om]
    x, y = -60.0, y_row + b_signed
    px, py = np.sqrt(bilin_f(n2, x, y)), 0.0
    ds = 0.02
    sgn = -1.0 if b_signed > 0 else +1.0

    def F(st):
        return np.array([st[2], st[3],
                         0.5 * bilin_f(gxg, st[0], st[1]),
                         0.5 * bilin_f(gyg, st[0], st[1])])

    st = np.array([x, y, px, py])
    out = {}
    todo = sorted(planes)
    for _ in range(500000):
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


for om in (1.0, 1.1, 1.3):
    print(f"    omega={om:.1f}:")
    for bb in (6.0, 8.0, 12.0, 16.0):
        res_t = trace_bg(om, bb)
        orient = ORIENT.get((om, int(bb)))
        otxt = f"  [orient para: {orient:+.2f}]" if orient else ""
        if res_t is None:
            print(f"      b={bb:4.0f}:  PRZECHWYT (r<0.6 od rdzenia){otxt}")
            continue
        row = "  ".join(f"x={xp:+.0f}: {res_t[xp]:+8.3f}"
                        for xp in PLANES if xp in res_t)
        mid = res_t.get(0.0, np.nan)
        spread = max(res_t.values()) - min(res_t.values())
        print(f"      b={bb:4.0f}:  {row}")
        print(f"                (srodek {mid:+.3f}; rozrzut plaszczyzn "
              f"{spread:.3f} deg){otxt}")
print()

# --- 6. budzet czasu (sekcja 1.5) ---
print("--- 6. budzet czasu ---")
for om in (1.0, 1.1, 1.3):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    vg = kappa * k0 / om
    T_trig = (-22.0 - (-60.0)) / vg           # start -> prog x1+10
    T_last = T_trig + 16.0                     # 5 probek co 4
    T_exit = (32.0 + 10.0 - (-60.0)) / vg      # za druga soczewke (D16)
    print(f"  omega={om:.1f}: v_g={vg:.3f}; wyzwolenie ~tau={T_trig:5.1f}; "
          f"5. probka ~tau={T_last:5.1f}; za 2. soczewka ~tau={T_exit:5.1f}")
print(f"  bramka G2: okno tau=400, kryterium tau<=300 "
      f">= 3x przelot (~120 przy omega=1.1) — spelnia definicje 1.5")
print()
print("=== WNIOSKI PHASE 0 ===")
print("  bilans sil siatki = 0 z symetrii (powloki kasuja sie do maszyn.);")
print("  zaklad eikonalny WIAZACY policzony na rzeczywistym tle siatki")
print("  wraz z czuloscia plaszczyzny ewaluacji (D6);")
print("  Phase 1 = bramka G2 (czas zycia) PRZED jakimkolwiek pakietem")
