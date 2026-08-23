# Phase 2b: pre-rejestrowany wariant techniczny G3 — relaksacja x4
# (LOCK, regula decyzyjna: "G3 FAIL -> ... jeden dozwolony wariant:
#  dluzsze relaksowanie x4; jesli utrzymane -> stop")

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
    C_ = (V2[0], V1[1])
    th = theta_h_pair(V1, C_) + theta_v_pair(C_, V2)
    Amp = np.full_like(X, s_star)
    for (x0, y0) in (V1, V2):
        ddx = (X - x0 + 0.5 * L) % L - 0.5 * L
        ddy = (Y - y0 + 0.5 * L) % L - 0.5 * L
        rr = np.sqrt(ddx * ddx + ddy * ddy)
        Amp = Amp * rr / np.sqrt(rr * rr + 2.0 * xi2)
    return Amp * np.exp(1j * th)


def wrap(x):
    return (x + np.pi) % TWO_PI - np.pi


def core_positions(psi):
    th = np.angle(psi)
    dthx = wrap(np.roll(th, -1, 0) - th)
    dthy = wrap(np.roll(th, -1, 1) - th)
    W = (dthx + np.roll(dthy, -1, 0) - np.roll(dthx, -1, 1) - dthy) / TWO_PI
    Phi = np.abs(psi) ** 2
    out = {}
    for sign in (+1, -1):
        cells = np.argwhere(W * sign > 0.5)
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


print("=== Phase 2b: wariant x4 (LOCK, regula G3) ===")
psi = make_pair()
pos = {}
for stp in range(1, 8001):
    psi = psi + dt_flow * rhs(psi)
    if stp in (2000, 4000, 8000):
        pos[stp] = core_positions(psi)
        r = float(np.max(np.abs(rhs(psi))))
        p = pos[stp]
        print(f"  krok {stp}: +1 ({p[+1][0]:.3f},{p[+1][1]:.3f})  "
              f"-1 ({p[-1][0]:.3f},{p[-1][1]:.3f})  residuum={r:.3e}")
d48 = max(np.hypot(pos[8000][s][0] - pos[4000][s][0],
                   pos[8000][s][1] - pos[4000][s][1]) for s in (+1, -1))
print(f"  dryf pozycji (krok 4000 -> 8000) = {d48:.4f}  (prog 0.02)")
print(f"  -> wariant x4: {'PASS' if d48 <= 0.02 else 'FAIL — dryf STRUKTURALNY'}")
if d48 > 0.02:
    print("  interpretacja: konfiguracja szachownicowa = MAKSIMUM energii")
    print("  wzglednej pary na torusie; kazda asymetria inicjuje zjazd ->")
    print("  para (+1,-1) na torusie nie jest scisle stacjonarna. Wg reguly")
    print("  decyzyjnej LOCKa: STOP (G3 utrzymane FAIL).")
