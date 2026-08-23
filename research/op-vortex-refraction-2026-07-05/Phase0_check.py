# Phase 0 check: rachunki PRZED kodem produkcyjnym
# (LOCK: Phase0_balance.md, sekcje 1-2 + D5, D6)
#   1) k0(omega), v_g, lambda — kanal amplitudowy (rownanie #1)
#   2) tlo rzeczywiste: para (+1,-1), ansatz zlozony, relaksacja 2000
#   3) n^2(r) na tle pary vs ogon algebraiczny A_t/r^2
#   4) ray tracing na RZECZYWISTYM tle: alpha_pred(b, omega) + b=16
#      + porownanie z tabela cytowana (izolowany profil ODE, ktora
#      przeliczamy tez lokalnie z dodatkiem b=16)
#   5) budzet dryfu tla (cytowany 0.015/tau) -> uzasadnienie b_eff

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
V1 = (-32.0, -32.0)   # +1, wir rozpraszajacy
V2 = (+32.0, +32.0)   # -1

# tabela CYTOWANA (op-stationary-background/Phase0_output.txt)
CITED = {(1.0, 6): -158.225, (1.0, 8): +57.871, (1.0, 12): +15.743,
         (1.1, 6): +43.015, (1.1, 8): +14.006, (1.1, 12): +5.299,
         (1.3, 6): +12.448, (1.3, 8): +5.236, (1.3, 12): +2.108}


def gfun(m):
    return a - b_par * np.abs(m) + c * m * m


def Vpp(m):
    return a - 2.0 * b_par * np.abs(m) + 3.0 * c * m * m


def Vprime(m):
    return m * (a - b_par * np.abs(m) + c * m * m)


def lap(u):
    return (np.roll(u, 1, 0) + np.roll(u, -1, 0)
            + np.roll(u, 1, 1) + np.roll(u, -1, 1) - 4.0 * u) / (dx * dx)


def rhs(psi):
    return kappa * lap(psi) - psi * gfun(np.abs(psi))


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


print("=== Phase 0 check: op-vortex-refraction (LOCK: Phase0_balance.md) ===")
print()
print("--- 1. kanal amplitudowy: k0, v_g, lambda (rownanie #1) ---")
print(f"  s* = {s_star:.6f}  Phi* = {Phi_star:.6f}  V''(s*) = {Vpp_star:.6f}")
print(f"  m_TV = {m_TV:.4f}  xi^2 = {xi2:.4f}  A_t = {A_t:.4f}")
for om in (1.0, 1.1, 1.3):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    print(f"  omega={om:.1f}:  k0 = {k0:.4f}   lambda = {TWO_PI / k0:6.2f}   "
          f"v_g = {kappa * k0 / om:.4f}")
print(f"  UWAGA dyfrakcyjna (D12): lambda(1.1)=7.72 ~ b — rezim falowy,")
print(f"  stad szerokie pasmo G4 [0.5, 2.0]")
print()

# --- 2. tlo rzeczywiste ---
print("--- 2. tlo: para (+1,-1), ansatz zlozony, relaksacja 2000 krokow ---")
psi = make_pair()
for stp in range(1, 2001):
    psi = psi + dt_flow * rhs(psi)
psi_bg = psi
f_bg = np.abs(psi_bg)
res = float(np.max(np.abs(rhs(psi_bg))))
print(f"  residuum max|RHS| po relaksacji = {res:.3e} "
      f"(op-stationary: 8.1e-03; tlo QUASI-stacjonarne)")
print(f"  |psi| daleko od rdzeni: min ~ {float(f_bg[(np.hypot((X-V1[0]+L/2)%L-L/2,(Y-V1[1]+L/2)%L-L/2)>4)&(np.hypot((X-V2[0]+L/2)%L-L/2,(Y-V2[1]+L/2)%L-L/2)>4)].min()):.4f}  (s*={s_star:.4f})")
print()

# --- 3. n^2(r) na tle pary ---
print("--- 3. n^2(r) na tle pary (rownanie #1) vs ogon algebraiczny ---")
print("  (srednia z 4 kierunkow osiowych od rdzenia +1; kolumny omega"
      " 1.0/1.1/1.3)")


def f_at(x, y):
    fx = ((x - coords[0]) / dx) % N
    fy = ((y - coords[0]) / dx) % N
    i0, j0 = int(fx) % N, int(fy) % N
    i1, j1 = (i0 + 1) % N, (j0 + 1) % N
    tx, ty = fx - int(fx), fy - int(fy)
    return ((1 - tx) * (1 - ty) * f_bg[i0, j0] + tx * (1 - ty) * f_bg[i1, j0]
            + (1 - tx) * ty * f_bg[i0, j1] + tx * ty * f_bg[i1, j1])


for rq in (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0):
    fq = np.mean([f_at(V1[0] + rq, V1[1]), f_at(V1[0] - rq, V1[1]),
                  f_at(V1[0], V1[1] + rq), f_at(V1[0], V1[1] - rq)])
    row = f"   r={rq:5.1f}: "
    for om in (1.0, 1.1, 1.3):
        n2 = (om * om - Vpp(fq)) / (om * om - Vpp_star)
        row += f"  {n2:7.4f}"
    tail = 1.0 + A_t / (rq * rq * (1.21 - Vpp_star))
    row += f"   [ogon(1.1): {tail:.4f}]"
    print(row)
print()

# --- 4. ray tracing na RZECZYWISTYM tle (D6) ---
print("--- 4. ray tracing na rzeczywistym tle: alpha_pred(b, omega) [deg] ---")
print("    (D6: start x=-60, wyjscie x>+40; alpha>0 = KU wirowi;")
print("     w nawiasach: tabela CYTOWANA — izolowany profil ODE, -60/+60)")

n2_grid = {}
gx_grid, gy_grid = {}, {}
for om in (1.0, 1.1, 1.3):
    n2 = (om * om - Vpp(f_bg)) / (om * om - Vpp_star)
    n2_grid[om] = n2
    gx_grid[om] = (np.roll(n2, -1, 0) - np.roll(n2, 1, 0)) / (2 * dx)
    gy_grid[om] = (np.roll(n2, -1, 1) - np.roll(n2, 1, 1)) / (2 * dx)


def bilin(G, x, y):
    fx = ((x - coords[0]) / dx) % N
    fy = ((y - coords[0]) / dx) % N
    i0, j0 = int(fx) % N, int(fy) % N
    i1, j1 = (i0 + 1) % N, (j0 + 1) % N
    tx, ty = fx - int(fx), fy - int(fy)
    return ((1 - tx) * (1 - ty) * G[i0, j0] + tx * (1 - ty) * G[i1, j0]
            + (1 - tx) * ty * G[i0, j1] + tx * ty * G[i1, j1])


def trace_bg(om, b, x_start=-60.0, x_exit=40.0, max_steps=500000):
    """Ray tracing na tle siatkowym; zwraca (alpha_ku_deg, status)."""
    n2, gxg, gyg = n2_grid[om], gx_grid[om], gy_grid[om]
    x, y = x_start, V1[1] + b
    px, py = np.sqrt(bilin(n2, x, y)), 0.0
    ds = 0.02

    def F(st):
        return np.array([st[2], st[3],
                         0.5 * bilin(gxg, st[0], st[1]),
                         0.5 * bilin(gyg, st[0], st[1])])

    st = np.array([x, y, px, py])
    for _ in range(max_steps):
        k1 = F(st)
        k2 = F(st + 0.5 * ds * k1)
        k3 = F(st + 0.5 * ds * k2)
        k4 = F(st + ds * k3)
        st = st + ds / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        if st[0] > x_exit:
            # alpha KU wirowi dla b>0: py<0 -> -atan2; ogolnie -sign(b)*atan2
            sgn = -1.0 if b > 0 else +1.0
            return sgn * np.degrees(np.arctan2(st[3], st[2])), "ok"
        rr = np.hypot(st[0] - V1[0], st[1] - V1[1])
        if rr < 0.6:
            return np.nan, "przechwyt(rdzen)"
    return np.nan, "przechwyt(petla)"


# 4a. profil izolowany ODE (metoda cytowana) — kontrola + b=16
dr = 0.05
r1d = np.arange(1, int(30.0 / dr)) * dr
fp = s_star * r1d / np.sqrt(r1d * r1d + 2.0 * xi2)
dt_ode = 0.4 * dr * dr / (2.0 * kappa)
for it in range(400000):
    fl = np.concatenate(([0.0], fp[:-1]))
    fr = np.concatenate((fp[1:], [s_star]))
    d2 = (fl + fr - 2.0 * fp) / (dr * dr)
    d1 = (fr - fl) / (2.0 * dr)
    rr_ = kappa * (d2 + d1 / r1d - fp / (r1d * r1d)) - Vprime(fp)
    fp = fp + dt_ode * rr_
    if it % 2000 == 0 and np.max(np.abs(rr_)) < 1e-9:
        break


def trace_ode(om, b):
    n2v = (om * om - Vpp(fp)) / (om * om - Vpp_star)
    dn2 = np.gradient(n2v, dr)

    def n2_at(rr):
        return np.interp(rr, r1d, n2v, left=n2v[0], right=1.0)

    def dn2_at(rr):
        return np.interp(rr, r1d, dn2, left=0.0, right=0.0)

    x, y = -60.0, b
    px, py = np.sqrt(n2_at(np.hypot(x, y))), 0.0
    ds = 0.02
    st = np.array([x, y, px, py])
    for _ in range(500000):
        def F(s_):
            rr = max(np.hypot(s_[0], s_[1]), 1e-9)
            g = 0.5 * dn2_at(rr)
            return np.array([s_[2], s_[3], g * s_[0] / rr, g * s_[1] / rr])
        k1 = F(st)
        k2 = F(st + 0.5 * ds * k1)
        k3 = F(st + 0.5 * ds * k2)
        k4 = F(st + ds * k3)
        st = st + ds / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        if st[0] > 60.0:
            return -np.degrees(np.arctan2(st[3], st[2]))
    return np.nan


B_LIST = (6.0, 8.0, 12.0, 16.0)
print("    b:            " + "".join(f"  {bb:8.0f}" for bb in B_LIST))
table_bg = {}
for om in (1.0, 1.1, 1.3):
    row = f"    omega={om:.1f}:  "
    row2 = "      (ODE):    "
    for bb in B_LIST:
        al, stt = trace_bg(om, bb)
        table_bg[(om, bb)] = (al, stt)
        row += f"  {al:+8.3f}" if np.isfinite(al) else f"  {'PRZECHW':>8}"
        alo = trace_ode(om, bb)
        row2 += f"  {alo:+8.3f}" if np.isfinite(alo) else f"  {'PRZECHW':>8}"
    print(row)
    print(row2)
print()
print("    roznice tlo-rzeczywiste vs cytowane (punkty wspolne):")
for (om, bb), cited in sorted(CITED.items()):
    al, stt = table_bg[(om, float(bb))]
    if np.isfinite(al) and abs(cited) < 90:
        print(f"      omega={om:.1f} b={bb:2d}: tlo {al:+8.3f} vs cyt "
              f"{cited:+8.3f}  ({100*(al-cited)/abs(cited):+5.1f}%)")
    else:
        print(f"      omega={om:.1f} b={bb:2d}: tlo "
              f"{al if np.isfinite(al) else stt} vs cyt {cited:+8.3f} "
              f"(strefa silna/przechwyt — poza kryterium)")
print()

# --- 5. budzet dryfu ---
print("--- 5. budzet dryfu tla (#1.4) ---")
print("  dryf pary (cytowany, op-stationary Phase2b): ~0.015/tau,"
      " przyspieszajacy")
for om in (1.0, 1.1, 1.3):
    k0 = np.sqrt((om * om - Vpp_star) / kappa)
    vg = kappa * k0 / om
    T = (abs(-60.0 - (V1[0] + 10.0))) / vg   # start -> okno x_rdzen+10
    print(f"  omega={om:.1f}: v_g={vg:.3f}, przelot do okna ~tau={T:5.1f} "
          f"-> przesuniecie rdzenia ~{0.015 * T:.2f}")
print("  wniosek: korekta b_eff (D5) + pozycje z pola (psi_ref, D3)"
      " OBOWIAZKOWE")
print()
print("=== WNIOSKI PHASE 0 ===")
print("  zaklad eikonalny na rzeczywistym tle policzony i zapisany;")
print("  Phase2 liczy alpha_pred PRZY b_eff na tym samym tle (D5);")
print("  kryterialne punkty: omega {1.1,1.3} x b {8,12,16}; pasmo [0.5,2.0]")
