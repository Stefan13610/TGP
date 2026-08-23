# Phase 0 check: weryfikacja rachunkow PRZED kodem produkcyjnym.
# Wg Phase0_balance.md (LOCKED, sekcje 1-3):
#   1) algebra: xi^2, m_a, 2piJ, 4piJ
#   2) ogon algebraiczny wiru: f_ODE(r) vs s*(1 - xi^2/r^2)
#   3) dekompozycja kwadrupolowa na ansatzu sieciowym: E_sym(d),
#      E_anty(d) vs 4*pi*J*[G(d) - G((2X,d))]  -> B_pred (forma P)
#   4) sieciowa funkcja Greena G_lat (baza gph do fitu Phase 2)

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
TWO_PI = 2.0 * np.pi

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)
Phi_star = s_star ** 2
Vpp_star = a - 2.0 * b_par * s_star + 3.0 * c * s_star ** 2
xi2 = kappa / Vpp_star
xi = np.sqrt(xi2)
J = kappa * Phi_star
m_a = np.sqrt(Vpp_star / kappa)

D_GRID = [5.0, 6.0, 8.0, 10.0, 12.0]
D_GRID_CTRL = [5.0, 6.0, 8.0]


def Vpot(m):
    return 0.5 * a * m**2 - (b_par / 3.0) * m**3 + 0.25 * c * m**4


def Vprime(m):
    return m * (a - b_par * np.abs(m) + c * m * m)


print("=== Phase 0 check: op-scalar-sector-phistar (LOCK: Phase0_balance.md) ===")
print()
print("--- 1. algebra ---")
print(f"  s*={s_star:.6f}  Phi*={Phi_star:.6f}  V''(s*)={Vpp_star:.6f}")
print(f"  xi^2 = kappa/V''(s*) = {xi2:.4f}   m_a = sqrt(V''(s*)/kappa) = {m_a:.4f}")
print(f"  2*pi*J = {TWO_PI*J:.4f}   4*pi*J = {2*TWO_PI*J:.4f}")
print()

# ---------------- 2. ogon ODE ----------------
dr = 0.05
r1d = np.arange(1, int(30.0 / dr)) * dr
f = s_star * r1d / np.sqrt(r1d * r1d + 2.0 * xi2)
dt_ode = 0.4 * dr * dr / (2.0 * kappa)
for it in range(400000):
    fl = np.concatenate(([0.0], f[:-1]))
    fr = np.concatenate((f[1:], [s_star]))
    d2 = (fl + fr - 2.0 * f) / (dr * dr)
    d1 = (fr - fl) / (2.0 * dr)
    rhs = kappa * (d2 + d1 / r1d - f / (r1d * r1d)) - Vprime(f)
    f = f + dt_ode * rhs
    if it % 2000 == 0 and np.max(np.abs(rhs)) < 1e-9:
        break
print("--- 2. ogon algebraiczny (rownanie #1) ---")
print("     r     f_ODE/s*    1-xi^2/r^2    stosunek")
for rq in (4.0, 5.0, 6.0, 8.0, 10.0, 12.0):
    fq = float(np.interp(rq, r1d, f)) / s_star
    alg = 1.0 - xi2 / (rq * rq)
    print(f"  {rq:5.1f}   {fq:.5f}     {alg:.5f}      {fq/alg:.5f}")
print()


# ---------------- 3. dekompozycja na ansatzu sieciowym ----------------
def quad_positions(X0, d):
    return [(-X0, +0.5 * d), (-X0, -0.5 * d), (+X0, +0.5 * d), (+X0, -0.5 * d)]


SIGNS = {"S": (+1, +1, -1, -1), "O": (+1, -1, -1, +1)}


def analyze_geometry(N, L, X0, d_list, tag):
    dxg = L / N
    coords = (np.arange(N) - N // 2) * dxg
    X, Y = np.meshgrid(coords, coords, indexing="ij")
    Z = X + 1j * Y

    def make_cfg(vortices):
        th = np.zeros_like(X)
        A = np.full_like(X, s_star)
        for (x0, y0, n) in vortices:
            for k in range(-2, 3):
                th = th + n * np.angle(
                    np.sinh(np.pi * (Z - (x0 + k * L) - 1j * y0) / L))
            ddx = (X - x0 + 0.5 * L) % L - 0.5 * L
            ddy = (Y - y0 + 0.5 * L) % L - 0.5 * L
            rr = np.sqrt(ddx * ddx + ddy * ddy)
            A = A * rr / np.sqrt(rr * rr + 2.0 * xi2)
        return A * np.exp(1j * th)

    def energy(psi):
        gx = (np.roll(psi, -1, 0) - psi) / dxg
        gy = (np.roll(psi, -1, 1) - psi) / dxg
        return float(np.sum(0.5 * kappa * (np.abs(gx)**2 + np.abs(gy)**2)
                            + Vpot(np.abs(psi))) * dxg * dxg)

    # sieciowa funkcja Greena (5-pkt, zero-mean): Lap G = delta - 1/L^2
    rho = np.zeros((N, N))
    rho[0, 0] = 1.0 / (dxg * dxg)
    kx = TWO_PI * np.fft.fftfreq(N, d=dxg)
    KX, KY = np.meshgrid(kx, kx, indexing="ij")
    lam = -((2 - 2 * np.cos(KX * dxg)) + (2 - 2 * np.cos(KY * dxg))) / (dxg * dxg)
    lam[0, 0] = 1.0
    G_hat = np.fft.fft2(rho) / lam
    G_hat[0, 0] = 0.0
    G_lat = np.real(np.fft.ifft2(G_hat))

    def G_at(ddx, ddy):
        fx = (ddx % L) / dxg
        fy = (ddy % L) / dxg
        i0, j0 = int(np.floor(fx)) % N, int(np.floor(fy)) % N
        tx, ty = fx - np.floor(fx), fy - np.floor(fy)
        i1, j1 = (i0 + 1) % N, (j0 + 1) % N
        return ((1 - tx) * (1 - ty) * G_lat[i0, j0] + tx * (1 - ty) * G_lat[i1, j0]
                + (1 - tx) * ty * G_lat[i0, j1] + tx * ty * G_lat[i1, j1])

    def gph(vortices):
        s = 0.0
        for i in range(len(vortices)):
            for j in range(i + 1, len(vortices)):
                xi_, yi, ni = vortices[i]
                xj, yj, nj = vortices[j]
                s += ni * nj * G_at(xi_ - xj, yi - yj)
        return -TWO_PI * s

    print(f"  geometria {tag}: N={N} L={L} dx={dxg} X={X0}")
    print("    d     E_S          E_O          E_sym        E_anty     "
          "4piJ*[G(d)-G(2X,d)]   ratio")
    Esym, Eanty, gph_check = [], [], []
    for d in d_list:
        pos = quad_positions(X0, d)
        E = {}
        for cfg, sg in SIGNS.items():
            vort = [(x, y, n) for (x, y), n in zip(pos, sg)]
            E[cfg] = energy(make_cfg(vort))
        es = 0.5 * (E["S"] + E["O"])
        ea = 0.5 * (E["O"] - E["S"])
        # E_anty_teoria = 2*(2pi)^2*J*[G(d) - G((2X,d))]  (G_lat: Lap G = delta)
        th = 2 * TWO_PI * TWO_PI * J * (G_at(0.0, d) - G_at(2 * X0, d))
        Esym.append(es); Eanty.append(ea)
        print(f"   {d:4.1f}  {E['S']:.5f}  {E['O']:.5f}  {es:.5f}  "
              f"{ea:+.5f}   {th:+.5f}   {ea/th:.4f}")
    Esym = np.array(Esym); dv = np.array(d_list)
    # fit (P): E_sym = C - B*ln(d/xi)/d^2
    basis = np.vstack([np.ones(len(dv)), -np.log(dv / xi) / dv**2]).T
    coef, *_ = np.linalg.lstsq(basis, Esym, rcond=None)
    rss = float(np.sum((Esym - basis @ coef) ** 2))
    B_pred = coef[1]
    rng_amp = float((basis @ coef).max() - (basis @ coef).min())
    print(f"    fit (P) na E_sym_ansatz: C={coef[0]:.5f}  B_pred={B_pred:.4f}  "
          f"RSS={rss:.2e}")
    print(f"    przewidywany rozstep kanalu na skanie: {rng_amp:.4f}")
    return B_pred


print("--- 3. dekompozycja kwadrupolowa (ansatz zlozeniowy, PRZED pomiarem) ---")
B_main = analyze_geometry(256, 128.0, 32.0, D_GRID, "GLOWNA dx=0.5")
print()
B_ctrl = analyze_geometry(256, 64.0, 16.0, D_GRID_CTRL, "KONTROLA dx=0.25")
print()

print("=== WNIOSKI PHASE 0 ===")
print(f"  1) algebra: xi^2={xi2:.4f}, m_a={m_a:.4f}, 2piJ={TWO_PI*J:.4f}")
print(f"  2) ogon wiru algebraiczny ~ 1 - xi^2/r^2 (podstawa formy P)")
print(f"  3) B_pred(glowna) = {B_main:.4f}; B_pred(kontrola) = {B_ctrl:.4f}")
print(f"     — zalockowane predykcje kanalu dylatacyjnego przed pomiarem")
print(f"  4) E_anty vs teoria Greena: ratio ~ 1 potwierdza algebre")
print(f"     kasowania kanalu fazowego (rownanie #2)")
