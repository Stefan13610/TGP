# Phase 0 check: rachunki PRZED kodem produkcyjnym
# (LOCK: Phase0_balance.md, sekcje 1-2 + D6, D13, D14)
#   1) k0(omega), v_g, lambda — kanal amplitudowy (rownanie #1)
#      + TABELA PROJEKTOWA punktow kryterialnych (rownanie #2):
#      b_nom vs 1.5*lambda + 0.9; oczekiwane b_eff/lambda
#   2) BILANS SIL na siatce szachownicowej L=256: powloki D4
#      CENTROWANE NA WIRZE (wzorzec erraty #1 proby #2);
#      kontrola: obciecie pudelkowe (czlon powierzchniowy ~K^-2)
#   3) dyspersja sieciowa ANALITYCZNIE przy k0(1.1..1.7)
#      (pre-rejestracja 1.4: prog 5%; k=1.66 ~ -1.9%, k=2.01 ~ -2.9%)
#   4) tlo rzeczywiste: szachownica 4 wirow L=256 (ansatz = suma
#      dwoch par poziomych theta_h_pair, BEZ transpozycji),
#      relaksacja 2000; diagnostyka szwu, detektor 2+2 (D14), residuum
#   5) n^2 na zrelaksowanym tle siatki vs ogon algebraiczny A_t/r^2
#   6) ray tracing WIAZACY na RZECZYWISTYM tle siatki L=256:
#      alpha_pred(b, omega) z CZULOSCIA plaszczyzny ewaluacji (D6):
#      x_eval in {-54, -27, 0, +27, +54}; porownanie z ekstrapolacja
#      orientacyjna z tabeli L=128 (proba #2, sekcja 1.3)
#   7) budzet czasu i kosztu (sekcja 1.6)

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt_flow = 0.02
N, L = 512, 256.0
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

# szachownica L=256 (sekcja 2): +1 (-64,-64), -1 (+64,-64),
#                               -1 (-64,+64), +1 (+64,+64)
VORTS = [((-64.0, -64.0), +1), ((+64.0, -64.0), -1),
         ((-64.0, +64.0), -1), ((+64.0, +64.0), +1)]
NOM = {+1: [(-64.0, -64.0), (+64.0, +64.0)],
       -1: [(+64.0, -64.0), (-64.0, +64.0)]}
V_SC = (-64.0, -64.0)   # rdzen rozpraszajacy (+1)
Y_ROW = -64.0

# ekstrapolacje ORIENTACYJNE (sekcja 1.3; tabela L=128 proby #2
# + skalowanie ogona dla omega=1.5)
ORIENT = {(1.1, 14): +3.9, (1.1, 20): +1.8,
          (1.5, 8): +3.2, (1.5, 12): +1.3, (1.5, 16): +0.7}
# punkty kryterialne (rownanie #2, ZAMKNIETE)
CRIT = [(1.1, 14.0), (1.1, 16.0), (1.1, 20.0),
        (1.5, 8.0), (1.5, 12.0), (1.5, 16.0)]
# obserwacje poza kryteriami
OBS = [(1.1, 8.0), (1.3, 8.0), (1.7, 8.0)]


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
    """Para wspolliniowa pozioma: +1 w z1, -1 w z2 (wzorzec proby #2)."""
    th = np.zeros_like(X)
    for (x0, y0), n in ((z1, +1), (z2, -1)):
        for k in range(-2, 3):
            th = th + n * np.angle(np.sinh(np.pi * (Z - (x0 + k * L) - 1j * y0) / L))
    return th


def make_lattice(mirror=False):
    """Szachownica: rzad y=-64: +1/-1; rzad y=+64: +1 w (+64), -1 w (-64)."""
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
    """Detektor 2+2 (D14), nominaly +/-64."""
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


print("=== Phase 0 check: op-shortwave-lattice (LOCK: Phase0_balance.md) ===")
print()
print("--- 1. kanal amplitudowy + TABELA PROJEKTOWA (rownania #1, #2) ---")
print(f"  s* = {s_star:.6f}  Phi* = {Phi_star:.6f}  V''(s*) = {Vpp_star:.6f}")
print(f"  m_TV = {m_TV:.4f}  xi^2 = {xi2:.4f}  A_t = {A_t:.4f}")
lam = {}
for om in (1.1, 1.3, 1.5, 1.7):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    lam[om] = TWO_PI / k0
    print(f"  omega={om:.1f}:  k0 = {k0:.4f}   lambda = {TWO_PI / k0:6.2f}   "
          f"v_g = {kappa * k0 / om:.4f}")
print()
print("  TABELA PROJEKTOWA (kryterium: b_nom >= 1.5*lambda + 0.9;")
print("  0.9 = najwiekszy zmierzony dip b_eff w probie #2 [0.32, 0.84]):")
for om, bb in CRIT:
    thr = 1.5 * lam[om] + 0.9
    print(f"    omega={om:.1f} b={bb:4.0f}: prog = {thr:5.2f}  "
          f"b_nom - prog = {bb - thr:+5.2f}  "
          f"b_nom/lambda = {bb / lam[om]:.2f}  "
          f"(b_nom-0.9)/lambda = {(bb - 0.9) / lam[om]:.2f}"
          f"  -> {'W DOMENIE' if bb >= thr else 'POZA (BLAD PROJEKTU)'}")
print("  obserwacje POZA kryteriami:")
for om, bb in OBS:
    print(f"    omega={om:.1f} b={bb:4.0f}: b_nom/lambda = "
          f"{bb / lam[om]:.2f}  (rezim {'dyfrakcyjny' if bb/lam[om] < 1.5 else 'krotkofalowy'})")
print()

# --- 2. bilans sil na siatce L=256 ---
print("--- 2. bilans sil (punktowo-wirowy, powloki D4 CENTROWANE "
      "NA WIRZE) ---")
print("    siatka szachownicowa na torusie L=256: sasiedzi wiru +1")
print("    w (-64,-64) leza w (128i, 128j), ladunek (-1)^(i+j);")
print("    powloka K = max(|i|,|j|) = unia orbit D4 o wspolnym ladunku")
F_lat = np.zeros(2)
for k in range(1, 9):
    Fs = np.zeros(2)
    for i in range(-k, k + 1):
        for j in range(-k, k + 1):
            if max(abs(i), abs(j)) != k:
                continue
            qj = 1.0 if (i + j) % 2 == 0 else -1.0
            rx, ry = -128.0 * i, -128.0 * j
            r2 = rx * rx + ry * ry
            Fs += qj * np.array([rx, ry]) / r2
    F_lat += Fs
    print(f"      powloka K={k}: |F_powloki| = {float(np.hypot(*Fs)):.3e}"
          f"   |F_suma| = {float(np.hypot(*F_lat)):.3e}")
print(f"      -> F_netto = ({F_lat[0]:+.3e}, {F_lat[1]:+.3e})  "
      f"(kasowanie DOKLADNE z symetrii D4, do zaokraglen)")


def force_on_box(idx, vortlist, K):
    (xi_, yi_), qi = vortlist[idx]
    F = np.zeros(2)
    for kx in range(-K, K + 1):
        for ky in range(-K, K + 1):
            for j, ((xj, yj), qj) in enumerate(vortlist):
                if j == idx and kx == 0 and ky == 0:
                    continue
                rx = xi_ - (xj + kx * L)
                ry = yi_ - (yj + ky * L)
                r2 = rx * rx + ry * ry
                F += qi * qj * np.array([rx, ry]) / r2
    return F


F_box = force_on_box(0, VORTS, 8)
print(f"      (kontrola: obciecie pudelkowe K=8 daje "
      f"|F| = {float(np.hypot(*F_box)):.3e} — czlon powierzchniowy "
      f"obciecia ~K^-2, artefakt, nie sila; errata #1 proby #2)")
print("    KONTRAST (cytowane, bez powtorki): para proby #1 |F| ~ 1.7e-02,")
print("    anihilacja tau=107.5; siatka L=128 proby #2: 0.0000 przez tau=400")
print()

# --- 3. dyspersja sieciowa analitycznie (pre-rejestracja 1.4) ---
print("--- 3. dyspersja sieciowa ANALITYCZNIE (dx=0.5): "
      "omega_lat = sqrt(kappa*(2-2cos(k dx))/dx^2 + V''(s*)) ---")
for om in (1.1, 1.3, 1.5, 1.7):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    k2_lat = (2.0 - 2.0 * np.cos(k0 * dx)) / (dx * dx)
    om_lat = np.sqrt(kappa * k2_lat + Vpp_star)
    err = 100.0 * (om_lat - om) / om
    print(f"    omega={om:.1f} (k0={k0:.4f}, k0*dx={k0*dx:.3f}): "
          f"omega_lat = {om_lat:.5f}  blad = {err:+.2f}%"
          f"{'  [KRYTERIALNE pasmo G3]' if om in (1.1, 1.5) else '  [obserwacja]'}")
print("    prog G3: 5% (MIERZONE w Phase 1, nie zakladane)")
print()

# --- 4. tlo rzeczywiste: szachownica L=256 ---
print("--- 4. tlo: szachownica 4 wirow L=256, ansatz 2 pary poziome, "
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
      f"  (proba #2 przy L=128: 0.124; wadliwy ansatz: ~pi)")
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
      f"(proba #2, L=128: 1.6e-06)")
if cpR is not None:
    print(f"    pozycje rdzeni: "
          + "; ".join(f"{s:+d}#{i}: ({p[0]:+.3f},{p[1]:+.3f})"
                      for (s, i), p in sorted(cpR[0].items())))
    print(f"    n_tot = {cpR[1]:.0f}; max |poz - nominal| = "
          f"{max(pdist(p, NOM[s][i]) for (s, i), p in cpR[0].items()):.4f}")
print(f"    |psi| daleko od rdzeni: min = {float(f_bg[far].min()):.4f} "
      f"(s* = {s_star:.4f})")
print()

# --- 5. n^2 na zrelaksowanym tle siatki ---
print("--- 5. n^2 od rdzenia rozpraszajacego (srednia 4 kierunkow "
      "osiowych; kolumny omega 1.1/1.3/1.5/1.7) ---")


def bilin_f(G, x, y):
    fx = ((x - coords[0]) / dx) % N
    fy = ((y - coords[0]) / dx) % N
    i0, j0 = int(fx) % N, int(fy) % N
    i1, j1 = (i0 + 1) % N, (j0 + 1) % N
    tx, ty = fx - int(fx), fy - int(fy)
    return ((1 - tx) * (1 - ty) * G[i0, j0] + tx * (1 - ty) * G[i1, j0]
            + (1 - tx) * ty * G[i0, j1] + tx * ty * G[i1, j1])


sc = cpR[0][(+1, 0)] if cpR is not None else V_SC
for rq in (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 20.0, 24.0):
    fq = np.mean([bilin_f(f_bg, sc[0] + rq, sc[1]),
                  bilin_f(f_bg, sc[0] - rq, sc[1]),
                  bilin_f(f_bg, sc[0], sc[1] + rq),
                  bilin_f(f_bg, sc[0], sc[1] - rq)])
    row = f"   r={rq:5.1f}: "
    for om in (1.1, 1.3, 1.5, 1.7):
        n2 = (om * om - Vpp(fq)) / (om * om - Vpp_star)
        row += f"  {n2:7.4f}"
    tail = 1.0 + A_t / (rq * rq * (1.21 - Vpp_star))
    row += f"   [ogon(1.1): {tail:.4f}]"
    print(row)
print("   (uwaga: siatka nie jest radialna — sasiedzi w +/-x, +/-y "
      "w odleglosci 128)")
print()

# --- 6. ray tracing WIAZACY + czulosc plaszczyzny ewaluacji ---
print("--- 6. ray tracing WIAZACY (D6): alpha_pred(b, omega) [deg], kat")
print("    czytany w plaszczyznach x_eval in {-54, -27, 0, +27, +54}")
print("    (nominalne krance okna: x1+10 = -54, x2-10 = +54);")
print("    w nawiasie [orient]: ekstrapolacja z tabeli L=128 (sekcja 1.3)")

n2_grid, gx_grid, gy_grid = {}, {}, {}
for om in (1.1, 1.3, 1.5, 1.7):
    n2 = (om * om - Vpp(f_bg)) / (om * om - Vpp_star)
    n2_grid[om] = n2
    gx_grid[om] = (np.roll(n2, -1, 0) - np.roll(n2, 1, 0)) / (2 * dx)
    gy_grid[om] = (np.roll(n2, -1, 1) - np.roll(n2, 1, 1)) / (2 * dx)

PLANES = (-54.0, -27.0, 0.0, +27.0, +54.0)


def trace_bg(om, b_signed, planes=PLANES):
    """Ray na tle siatki L=256; dict {x_plane: alpha_ku_deg} albo None
    (przechwyt: r<0.6 od ktoregokolwiek rdzenia)."""
    n2, gxg, gyg = n2_grid[om], gx_grid[om], gy_grid[om]
    x, y = -100.0, Y_ROW + b_signed
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


for om in (1.1, 1.3, 1.5, 1.7):
    bset = {1.1: (8.0, 14.0, 16.0, 20.0), 1.3: (8.0,),
            1.5: (8.0, 12.0, 16.0), 1.7: (8.0,)}[om]
    print(f"    omega={om:.1f}:")
    for bb in bset:
        res_t = trace_bg(om, bb)
        orient = ORIENT.get((om, int(bb)))
        otxt = f"  [orient: {orient:+.2f}]" if orient else ""
        crit_t = "  [KRYT]" if (om, bb) in CRIT else "  [obs]"
        if res_t is None:
            print(f"      b={bb:4.0f}:{crit_t}  PRZECHWYT "
                  f"(r<0.6 od rdzenia){otxt}")
            continue
        row = "  ".join(f"x={xp:+.0f}: {res_t[xp]:+8.3f}"
                        for xp in PLANES if xp in res_t)
        mid = res_t.get(0.0, np.nan)
        spread = max(res_t.values()) - min(res_t.values())
        print(f"      b={bb:4.0f}:{crit_t}  {row}")
        print(f"                (srodek {mid:+.3f}; rozrzut plaszczyzn "
              f"{spread:.3f} deg){otxt}")
print()

# --- 7. budzet czasu i kosztu (sekcja 1.6) ---
print("--- 7. budzet czasu i kosztu ---")
for om in (1.1, 1.3, 1.5, 1.7):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    vg = kappa * k0 / om
    T_trig = (-54.0 - (-100.0)) / vg           # start -> prog x1+10
    T_last = T_trig + 16.0                     # 5 probek co 4
    T_exit = (64.0 + 10.0 - (-100.0)) / vg     # za druga soczewke (D16)
    print(f"  omega={om:.1f}: v_g={vg:.3f}; wyzwolenie ~tau={T_trig:5.1f}; "
          f"5. probka ~tau={T_last:5.1f}; za 2. soczewka ~tau={T_exit:5.1f}")
print(f"  bramka G2: okno tau=600, kryterium tau<=450 "
      f">= 3x przelot (~140 przy omega=1.1) — spelnia definicje 1.6")
print(f"  limit runow falowych tau=300 pokrywa 5. probke kazdej serii;")
print(f"  superpozycja (1.5,8): za 2. soczewke ~tau=315 -> dedykowany")
print(f"  run 10000 krokow (tau=500, wyjatek D9/D16)")
print(f"  koszt: N=512 -> ~4x koszt kroku proby #2 (informacyjnie)")
print()
print("=== WNIOSKI PHASE 0 ===")
print("  wszystkie 6 punktow kryterialnych W DOMENIE eikonalu")
print("  (b_nom >= 1.5*lambda + 0.9); bilans sil siatki L=256 = 0")
print("  z symetrii (powloki kasuja sie do maszynowych); dyspersja")
print("  sieciowa przy k0(1.5) w progu 5% analitycznie (pomiar: G3);")
print("  zaklad eikonalny WIAZACY policzony na rzeczywistym tle L=256")
print("  wraz z czuloscia plaszczyzny ewaluacji (D6);")
print("  Phase 1 = bramka G2 (czas zycia, tau=600) PRZED jakimkolwiek")
print("  pakietem")
