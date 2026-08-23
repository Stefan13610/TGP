# Phase 0 check: rachunki PRZED kodem produkcyjnym
# (LOCK: Phase0_balance.md proby #4, sekcje 1-2 + D6, D12-D14, D16)
#   1) k0(omega), v_g, lambda + TABELA PROJEKTOWA par kryterialnych:
#      b_nom vs 1.5*lambda + 0.9 ORAZ odleglosci wiazki od soczewek
#      rzedow posrednich (64-|b|, 64+|b|) vs 1.5*lambda (D12)
#   2) BILANS SIL na siatce SKOSNEJ: PARY INWERSYJNE v/-v (D16);
#      kontrola negatywna: obciecie polplaszczyznowe NIE kasuje sie;
#      kontrola pudelkowa (artefakt powierzchniowy) informacyjnie
#   3) weryfikacja LUSTRZANOSCI POZYCJI (y -> -128-y wzgl. y=-64)
#      i ZLAMANIA LUSTRA PRZEZ LADUNKI; weryfikacja S-inv (inwersja
#      wokol A) i S-conj (sprzezenie = translacja o (128,0))
#   4) dyspersja sieciowa ANALITYCZNIE (pre-rejestracja, prog G3 5%)
#   5) tlo rzeczywiste: siatka skosna 8 wirow (ansatz = suma CZTERECH
#      par poziomych theta_h_pair), relaksacja 2000; szew, detektor
#      4+4 (D14), residuum; diagnostyki bitowe symetrii tla (D11)
#   6) n^2 na zrelaksowanym tle vs ogon algebraiczny A_t/r^2
#   7) ray tracing WIAZACY: alpha_pred dla +b i -b OSOBNO na kazdym
#      punkcie kryterialnym, czulosc plaszczyzn {-54,-27,0,+27,+54},
#      Delta_eik per para (D18; projektowo przy b_nom)
#   8) budzet czasu

import numpy as np
from itertools import permutations

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

# SIATKA SKOSNA (sekcja 0): a1=(128,0), a2=(64,64), q=(-1)^(n1+n2)
VORTS = [((-64.0, -64.0), +1), ((+64.0, -64.0), -1),
         ((-128.0, 0.0), +1), ((0.0, 0.0), -1),
         ((+64.0, +64.0), +1), ((-64.0, +64.0), -1),
         ((0.0, -128.0), +1), ((+128.0, -128.0), -1)]
NOM = {+1: [(-64.0, -64.0), (-128.0, 0.0), (+64.0, +64.0), (0.0, -128.0)],
       -1: [(+64.0, -64.0), (0.0, 0.0), (-64.0, +64.0), (+128.0, -128.0)]}
V_SC = (-64.0, -64.0)   # rdzen A (rozpraszajacy, +1)
Y_ROW = -64.0

# pary kryterialne (ZAMKNIETE) i obserwacje (sekcja 2)
CRIT_PAIRS = [(1.1, 14.0), (1.1, 20.0), (1.1, 24.0),
              (1.5, 8.0), (1.5, 12.0)]
OBS_B = [(1.1, 28.0)]


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
    """Para wspolliniowa pozioma: +1 w z1, -1 w z2 (wzorzec prob #2/#3)."""
    th = np.zeros_like(X)
    for (x0, y0), n in ((z1, +1), (z2, -1)):
        for k in range(-2, 3):
            th = th + n * np.angle(np.sinh(np.pi * (Z - (x0 + k * L) - 1j * y0) / L))
    return th


def make_lattice():
    """Siatka skosna 8 wirow: 4 pary poziome (sekcja 2)."""
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
    """Detektor 4+4 (D14): DOKLADNIE 4 klastry na znak; przypisanie
    do nominalow po najmniejszej sumie odleglosci (24 permutacje)."""
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


print("=== Phase 0 check: op-asymmetric-lattice, proba #4 (LOCK: "
      "Phase0_balance.md) ===")
print()
print("--- 1. kanal amplitudowy + TABELA PROJEKTOWA (rownanie #1, D12) ---")
print(f"  s* = {s_star:.6f}  Phi* = {Phi_star:.6f}  V''(s*) = {Vpp_star:.6f}")
print(f"  m_TV = {m_TV:.4f}  xi^2 = {xi2:.4f}  A_t = {A_t:.4f}")
lam = {}
for om in (1.1, 1.5):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    lam[om] = TWO_PI / k0
    print(f"  omega={om:.1f}:  k0 = {k0:.4f}   lambda = {TWO_PI / k0:6.2f}   "
          f"v_g = {kappa * k0 / om:.4f}")
print()
print("  TABELA PROJEKTOWA (D12): b_nom >= 1.5*lambda + 0.9 ORAZ")
print("  odleglosc wiazki od soczewek rzedow posrednich >= 1.5*lambda")
print("  (soczewki posrednie w x=0: rzad y=0 [-1] i y=-128 [+1];")
print("   wiazka y = -64 +/- b: d_blizsza = 64-|b|, d_dalsza = 64+|b|)")
for om, bb in CRIT_PAIRS + OBS_B:
    thr = 1.5 * lam[om] + 0.9
    d_near, d_far = 64.0 - bb, 64.0 + bb
    ok_b = bb >= thr
    ok_d = d_near >= 1.5 * lam[om]
    tagc = "[KRYT +/-]" if (om, bb) in CRIT_PAIRS else "[obs +/-]"
    verdict = ("W DOMENIE" if (ok_b and ok_d) else "POZA (BLAD PROJEKTU)") \
        if (om, bb) in CRIT_PAIRS else "poza kryteriami"
    print(f"    omega={om:.1f} |b|={bb:4.0f} {tagc}: prog_b = {thr:5.2f} "
          f"(b-prog = {bb - thr:+5.2f}); b/lambda = {bb / lam[om]:.2f}; "
          f"d_posr = {d_near:.0f}/{d_far:.0f} "
          f"(prog {1.5 * lam[om]:.2f}: {'OK' if ok_d else 'ZA BLISKO'}) "
          f"-> {verdict}")
print()

# --- 2. bilans sil: PARY INWERSYJNE v/-v (D16) ---
print("--- 2. bilans sil (punktowo-wirowy): PARY INWERSYJNE v/-v (D16) ---")
print("    sasiedzi rdzenia: v = n1*a1 + n2*a2, a1=(128,0), a2=(64,64),")
print("    q(v) = (-1)^(n1+n2); q(-v) = q(v) -> pary v/-v kasuja sie")
print("    DOKLADNIE (F nieparzysta). UWAGA (README): NIE sumowac powlok")
print("    bez parowania — kasowanie jest wlasnoscia PAR, nie orbit D4.")
a1v = np.array([128.0, 0.0])
a2v = np.array([64.0, 64.0])
F_tot = np.zeros(2)
for K in range(1, 9):
    F_shell = np.zeros(2)
    max_pair = 0.0
    npairs = 0
    for n1 in range(-K, K + 1):
        for n2 in range(-K, K + 1):
            if max(abs(n1), abs(n2)) != K:
                continue
            if not (n1 > 0 or (n1 == 0 and n2 > 0)):
                continue          # reprezentant pary (D16)
            v = n1 * a1v + n2 * a2v
            q = 1.0 if (n1 + n2) % 2 == 0 else -1.0
            r2 = float(v @ v)
            F_pair = (-q * v / r2) + (-q * (-v) / r2)   # F(v) + F(-v)
            F_shell += F_pair
            max_pair = max(max_pair, float(np.hypot(*F_pair)))
            npairs += 1
    F_tot += F_shell
    print(f"      powloka indeksowa K={K}: par = {npairs:3d}; "
          f"max |F_pary| = {max_pair:.3e}; |F_powloki| = "
          f"{float(np.hypot(*F_shell)):.3e}; |F_suma| = "
          f"{float(np.hypot(*F_tot)):.3e}")
print(f"      -> F_netto (pary, K<=8) = ({F_tot[0]:+.3e}, {F_tot[1]:+.3e})"
      f"  (kasowanie DOKLADNE per para)")
# kontrola negatywna: polplaszczyzna (niezamknieta na inwersje)
F_half = np.zeros(2)
for n1 in range(-8, 9):
    for n2 in range(-8, 9):
        if (n1, n2) == (0, 0) or max(abs(n1), abs(n2)) > 8:
            continue
        if not (n1 > 0 or (n1 == 0 and n2 > 0)):
            continue
        v = n1 * a1v + n2 * a2v
        q = 1.0 if (n1 + n2) % 2 == 0 else -1.0
        F_half += -q * v / float(v @ v)
print(f"      kontrola NEGATYWNA (polplaszczyzna, bez partnerow -v): "
      f"|F| = {float(np.hypot(*F_half)):.3e} — NIE kasuje sie;")
print("      kasowanie pochodzi z par inwersyjnych, nie z D4")


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


Fb = force_on_box(0, VORTS, 8)
print(f"      kontrola pudelkowa (K=8, wszystkie 8 wirow + obrazy): "
      f"|F| = {float(np.hypot(*Fb)):.3e} — czlon powierzchniowy "
      f"obciecia ~K^-2 (artefakt, errata #1 proby #2)")
Fmax_all = 0.0
for idx in range(8):
    Fi = force_on_box(idx, VORTS, 8)
    Fmax_all = max(Fmax_all, float(np.hypot(*Fi)))
print(f"      max po WSZYSTKICH 8 rdzeniach (pudelko K=8): "
      f"{Fmax_all:.3e} (jednakowy artefakt — otoczenia identyczne "
      f"z inwersji)")
print()

# --- 3. lustro pozycji / zlamanie przez ladunki / S-inv / S-conj ---
print("--- 3. symetrie konfiguracji (weryfikacja konstrukcji, sekcja 1.2) ---")


def find_match(pos):
    for (p, q) in VORTS:
        if pdist(p, pos) < 1e-9:
            return (p, q)
    return None


# (a) lustro y -> -128 - y + sprzezenie (ladunki ZACHOWANE)
mism = []
for (p, q) in VORTS:
    img = (p[0], -128.0 - p[1])
    m = find_match(img)
    if m is None:
        mism.append((p, q, "brak pozycji"))
    elif m[1] != q:
        mism.append((p, q, f"ladunek {m[1]:+d} zamiast {q:+d} w "
                           f"({img[0]:+.0f},{img[1]:+.0f})"))
print("  (a) lustro osi wiazki (y -> -128-y) + sprzezenie "
      "[zachowuje ladunki]:")
print(f"      pozycje lustrzane: "
      f"{'TAK (8/8 pozycji przechodzi na pozycje)' if all(find_match((p[0], -128.0 - p[1])) is not None for p, _ in VORTS) else 'NIE'}")
print(f"      ladunki: NIEZGODNE w {len(mism)}/8 wirach -> lustro "
      f"ZLAMANE przez ladunki:")
for p, q, msg in mism:
    print(f"        wir {q:+d}@({p[0]:+.0f},{p[1]:+.0f}): {msg}")
# (b) S-inv: inwersja wokol A (zachowuje ladunki)
bad_inv = 0
for (p, q) in VORTS:
    img = (2 * V_SC[0] - p[0], 2 * V_SC[1] - p[1])
    m = find_match(img)
    if m is None or m[1] != q:
        bad_inv += 1
print(f"  (b) S-inv (r -> 2A - r wokol A=(-64,-64)): zgodnosc "
      f"pozycji+ladunkow: {8 - bad_inv}/8 "
      f"-> {'DOKLADNA symetria tla' if bad_inv == 0 else 'BLAD'}")
# (c) S-conj: sprzezenie (odwraca ladunki) = translacja o (128,0)
bad_conj = 0
for (p, q) in VORTS:
    img = (p[0] + 128.0, p[1])
    m = find_match((((img[0] + L / 2) % L) - L / 2, img[1]))
    if m is None or m[1] != -q:
        bad_conj += 1
print(f"  (c) S-conj (sprzezenie = translacja o (128,0)): translacja "
      f"kazdego wiru trafia w wir o PRZECIWNYM ladunku: "
      f"{8 - bad_conj}/8 -> "
      f"{'DOKLADNA rownowaznosc' if bad_conj == 0 else 'BLAD'}")
print()

# --- 4. dyspersja sieciowa analitycznie ---
print("--- 4. dyspersja sieciowa ANALITYCZNIE (dx=0.5): "
      "omega_lat = sqrt(kappa*(2-2cos(k dx))/dx^2 + V''(s*)) ---")
for om in (1.1, 1.5):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    k2_lat = (2.0 - 2.0 * np.cos(k0 * dx)) / (dx * dx)
    om_lat = np.sqrt(kappa * k2_lat + Vpp_star)
    err = 100.0 * (om_lat - om) / om
    print(f"    omega={om:.1f} (k0={k0:.4f}, k0*dx={k0*dx:.3f}): "
          f"omega_lat = {om_lat:.5f}  blad = {err:+.2f}%  "
          f"[KRYTERIALNE pasmo G3]")
print("    prog G3: 5% (MIERZONE w Phase 1, nie zakladane)")
print()

# --- 5. tlo rzeczywiste: siatka skosna 8 wirow ---
print("--- 5. tlo: siatka skosna 8 wirow L=256, ansatz 4 pary poziome, "
      "relaksacja 2000 ---")
psi = make_lattice()
th0 = np.angle(psi)
gxa = np.abs(wrap(np.roll(th0, -1, 0) - th0))
gya = np.abs(wrap(np.roll(th0, -1, 1) - th0))
far = np.ones((N, N), bool)
for (x0, y0), _n in VORTS:
    far &= np.hypot((X - x0 + L/2) % L - L/2, (Y - y0 + L/2) % L - L/2) > 4
print(f"  diagnostyka szwu (ansatz, 4 pary — nowa kombinacja): "
      f"max |dtheta| na linkach poza rdzeniami = "
      f"{max(gxa[far].max(), gya[far].max()):.4f} rad "
      f"(proba #3: 0.124; wadliwy ansatz: ~pi)")
cp0 = core_positions(psi)
if cp0 is None:
    print("  BLAD: detektor nie widzi 4+4 rdzeni w ansatzu")
else:
    print(f"  kret ansatzu: n_tot = {cp0[1]:.0f}; pozycje: "
          + "; ".join(f"{s:+d}#{i}:({p[0]:+.2f},{p[1]:+.2f})"
                      for (s, i), p in sorted(cp0[0].items())))
for stp in range(1, 2001):
    psi = psi + dt_flow * rhs(psi)
psi_bg = psi
f_bg = np.abs(psi_bg)
res = float(np.max(np.abs(rhs(psi_bg))))
cpR = core_positions(psi_bg)
print(f"  po relaksacji 2000 krokow (dt={dt_flow}):")
print(f"    residuum max|RHS| = {res:.3e} (proba #3: 9.6e-07)")
if cpR is not None:
    print(f"    pozycje rdzeni: "
          + "; ".join(f"{s:+d}#{i}:({p[0]:+.3f},{p[1]:+.3f})"
                      for (s, i), p in sorted(cpR[0].items())))
    print(f"    n_tot = {cpR[1]:.0f}; max |poz - nominal| = "
          f"{max(pdist(p, NOM[s][i]) for (s, i), p in cpR[0].items()):.4f}")
print(f"    |psi| daleko od rdzeni: min = {float(f_bg[far].min()):.4f} "
      f"(s* = {s_star:.4f})")
# diagnostyki bitowe symetrii tla (D11)
inv_field = psi_bg[(256 - np.arange(N)) % N][:, (256 - np.arange(N)) % N]
d_inv = float(np.max(np.abs(psi_bg - inv_field)))
d_conj = float(np.max(np.abs(np.roll(np.conj(psi_bg), N // 2, axis=0)
                             - psi_bg)))
print(f"    diagnostyka D11(a): max|psi_bg - I[psi_bg]| = {d_inv:.3e} "
      f"(inwersja wokol nominalnego A; 0 = bitowa)")
print(f"    diagnostyka D11(b): max|roll(conj(psi_bg), N/2) - psi_bg| "
      f"= {d_conj:.3e} (sprzezenie=translacja; 0 = bitowa)")
print()

# --- 6. n^2 na zrelaksowanym tle ---
print("--- 6. n^2 od rdzenia A (srednia 4 kierunkow osiowych; "
      "kolumny omega 1.1/1.5) ---")


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
    for om in (1.1, 1.5):
        n2 = (om * om - Vpp(fq)) / (om * om - Vpp_star)
        row += f"  {n2:7.4f}"
    tail = 1.0 + A_t / (rq * rq * (1.21 - Vpp_star))
    row += f"   [ogon(1.1): {tail:.4f}]"
    print(row)
print("   (gestosc 2x proby #3: najblizszy sasiad 90.5, nie 128)")
print()

# --- 7. ray tracing WIAZACY: +b i -b OSOBNO + Delta_eik per para ---
print("--- 7. ray tracing WIAZACY (D6/D18): alpha_pred(+/-b) [deg],")
print("    plaszczyzny x_eval in {-54, -27, 0, +27, +54}; projektowo")
print("    przy b_nom (WIAZACO per run: b_eff w Phase 2);")
print("    Delta_eik(para) = pred(+b) - pred(-b) w kazdej plaszczyznie")

n2_grid, gx_grid, gy_grid = {}, {}, {}
for om in (1.1, 1.5):
    n2 = (om * om - Vpp(f_bg)) / (om * om - Vpp_star)
    n2_grid[om] = n2
    gx_grid[om] = (np.roll(n2, -1, 0) - np.roll(n2, 1, 0)) / (2 * dx)
    gy_grid[om] = (np.roll(n2, -1, 1) - np.roll(n2, 1, 1)) / (2 * dx)

PLANES = (-54.0, -27.0, 0.0, +27.0, +54.0)


def trace_bg(om, b_signed, planes=PLANES):
    """Ray na tle siatki skosnej; dict {x_plane: alpha_ku_deg} albo
    None (przechwyt: r<0.6 od ktoregokolwiek z 8 rdzeni)."""
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


for om, bb in CRIT_PAIRS + OBS_B:
    tagc = "[KRYT]" if (om, bb) in CRIT_PAIRS else "[obs]"
    rp = trace_bg(om, +bb)
    rm = trace_bg(om, -bb)
    print(f"    omega={om:.1f} para |b|={bb:.0f} {tagc}:")
    for lbl, rr in ((f"+{bb:.0f}", rp), (f"-{bb:.0f}", rm)):
        if rr is None:
            print(f"      b={lbl:>4}: PRZECHWYT (r<0.6)")
            continue
        row = "  ".join(f"x={xp:+.0f}: {rr[xp]:+7.3f}"
                        for xp in PLANES if xp in rr)
        spread = max(rr.values()) - min(rr.values())
        print(f"      b={lbl:>4}: {row}  (rozrzut {spread:.3f})")
    if rp is not None and rm is not None:
        drow = "  ".join(f"x={xp:+.0f}: {rp[xp] - rm[xp]:+7.3f}"
                         for xp in PLANES)
        print(f"      Delta_eik(b_nom): {drow}")
print()

# --- 8. budzet czasu ---
print("--- 8. budzet czasu (sekcja 1.4) ---")
for om in (1.1, 1.5):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    vg = kappa * k0 / om
    T_trig = (-54.0 - (-100.0)) / vg
    T_last = T_trig + 16.0
    print(f"  omega={om:.1f}: v_g={vg:.3f}; wyzwolenie ~tau={T_trig:5.1f}; "
          f"5. probka ~tau={T_last:5.1f}")
print("  bramka G2: okno tau=600, kryterium tau<=450 >= 3x przelot "
      "(~140 przy omega=1.1)")
print("  limit runow falowych tau=300 pokrywa 5. probke kazdej serii")
print("  koszt: bramka ~10-15 min; run ~3-6 min; 13 runow Phase 2 "
      "+ 3 runy Phase 3 -> calosc ~1.5-2 h")
print()
print("=== WNIOSKI PHASE 0 ===")
print("  wszystkie 10 punktow kryterialnych W DOMENIE eikonalu")
print("  (b_nom >= 1.5*lambda+0.9; odleglosci od soczewek posrednich")
print("  >= 1.5*lambda); bilans sil = 0 DOKLADNIE parami inwersyjnymi")
print("  v/-v (kontrola negatywna polplaszczyzny potwierdza mechanizm);")
print("  lustro POZYCJI dokladne, ZLAMANE przez ladunki w 4/8 wirow")
print("  (oba rzedy posrednie) -> alpha_odd niewiazana zadna symetria;")
print("  S-inv i S-conj zweryfikowane na konfiguracji (testy G6a);")
print("  zaklad eikonalny WIAZACY policzony dla +b i -b OSOBNO")
print("  z Delta_eik per para; Phase 1 = bramka G2 PRZED pakietem")
