# Phase 1: ogon amplitudowy wiru na siatce 2D — G2 (+ czesc G1)
# Wg Phase0_balance.md (LOCKED): |f(r)/[s*(1-xi^2/r^2)] - 1| <= 1% dla r in [6,10]

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt = 0.02
N, L = 256, 128.0
dx = L / N
TAU_RELAX = 200
S_GUARD = 10.0
H_TOL = 1e-10
TWO_PI = 2.0 * np.pi

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)
Vpp_star = a - 2.0 * b_par * s_star + 3.0 * c * s_star ** 2
xi2 = kappa / Vpp_star

coords = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
Z = X + 1j * Y


def Vpot(m):
    return 0.5 * a * m**2 - (b_par / 3.0) * m**3 + 0.25 * c * m**4


def energy(psi):
    gx = (np.roll(psi, -1, 0) - psi) / dx
    gy = (np.roll(psi, -1, 1) - psi) / dx
    return float(np.sum(0.5 * kappa * (np.abs(gx)**2 + np.abs(gy)**2)
                        + Vpot(np.abs(psi))) * dx * dx)


def step(psi):
    m = np.abs(psi)
    lp = (np.roll(psi, 1, 0) + np.roll(psi, -1, 0)
          + np.roll(psi, 1, 1) + np.roll(psi, -1, 1) - 4.0 * psi) / (dx * dx)
    return psi + dt * (kappa * lp - psi * (a - b_par * m + c * m * m))


def make_cfg(vortices):
    th = np.zeros_like(X)
    A = np.full_like(X, s_star)
    for (x0, y0, n) in vortices:
        for k in range(-2, 3):
            th = th + n * np.angle(np.sinh(np.pi * (Z - (x0 + k * L) - 1j * y0) / L))
        ddx = (X - x0 + 0.5 * L) % L - 0.5 * L
        ddy = (Y - y0 + 0.5 * L) % L - 0.5 * L
        rr = np.sqrt(ddx * ddx + ddy * ddy)
        A = A * rr / np.sqrt(rr * rr + 2.0 * xi2)
    return A * np.exp(1j * th)


def relax(psi0, steps):
    psi = psi0.copy()
    H_prev, H_inc = energy(psi), 0.0
    fl = {"nonfinite": False, "runaway": False}
    for stp in range(1, steps + 1):
        psi = step(psi)
        if stp % 10 == 0:
            if not np.all(np.isfinite(psi)):
                fl["nonfinite"] = True
                break
            if np.max(np.abs(psi)) > S_GUARD:
                fl["runaway"] = True
                break
            H = energy(psi)
            if H > H_prev:
                H_inc = max(H_inc, H - H_prev)
            H_prev = H
    fl["H_mono"] = H_inc <= H_TOL * max(1.0, abs(energy(psi0)))
    return psi, fl


print("=== Phase 1: ogon amplitudowy wiru (LOCK: Phase0_balance.md) ===")
print(f"N={N} L={L} dx={dx} dt={dt} tau_relax={TAU_RELAX*dt}")
print()
psi0 = make_cfg([(-32.0, 0.0, +1), (+32.0, 0.0, -1)])
psi, fl = relax(psi0, TAU_RELAX)
psi_b, _ = relax(psi0, TAU_RELAX)
det = bool(np.array_equal(psi, psi_b))

ddx = (X + 32.0 + 0.5 * L) % L - 0.5 * L
ddy = (Y + 0.5 * L) % L - 0.5 * L
rr = np.hypot(ddx, ddy).ravel()
mm = np.abs(psi).ravel()
print("     r     f_2D/s*     1-xi^2/r^2    stosunek   |odchylka|")
edges = np.arange(5.75, 10.30, 0.5)
devs = []
for k in range(len(edges) - 1):
    sel = (rr >= edges[k]) & (rr < edges[k + 1])
    rc = 0.5 * (edges[k] + edges[k + 1])
    fv = float(mm[sel].mean()) / s_star
    alg = 1.0 - xi2 / (rc * rc)
    ratio = fv / alg
    devs.append(abs(ratio - 1.0))
    print(f"  {rc:5.2f}   {fv:.5f}    {alg:.5f}      {ratio:.5f}   {abs(ratio-1)*100:.3f}%")
g2 = max(devs) <= 0.01
print()
print("=== EWALUACJA (LOCKED) ===")
print(f"G1 (czesc Phase1): determinizm bitowy: {det}; nf={fl['nonfinite']} "
      f"run={fl['runaway']} H_mono={fl['H_mono']}")
print(f"G2 (ogon algebraiczny): max |odchylka| = {100*max(devs):.3f}% <= 1%: "
      f"{g2}  -> {'PASS' if g2 else 'FAIL'}")
print()
print(f"PODSUMOWANIE PHASE 1: G1(czesc)={'PASS' if det and fl['H_mono'] else 'FAIL'} "
      f"G2={'PASS' if g2 else 'FAIL'}")
