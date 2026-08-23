# Phase 1: statyczna addytywnosc parowa — G1-G4
# Implementacja SCISLE wg Phase0_balance.md (LOCKED) tego opa.
# E_2 liczone dla DOKLADNIE tych samych ofsetow calkowitoliczbowych co w
# konfiguracji N -> delta izoluje czysto >=3-cialowe czlony potencjalu.

import numpy as np

N = 256
L = 128.0
dx = L / N
kappa = 0.50
a, b, c = 0.50, 1.60, 1.00
dt = 0.02
EPS = 0.30
S_GUARD = 10.0
w = 1.50
REC_EVERY = 10
H_TOL = 1e-10
R_OBJ = 8.0
G_SCAN = [2.0, 2.5, 3.0, 4.0]
G2_SCAN = [2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
FLOOR = 1e-4

Nc = N // 2
coords = (np.arange(N) - Nc) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")


def Vpot(s):
    return 0.5 * a * s**2 - (b / 3.0) * np.abs(s)**3 + 0.25 * c * s**4


def Vprime(s):
    return s * (a - b * np.abs(s) + c * s * s)


def lap(s):
    return (np.roll(s, 1, 0) + np.roll(s, -1, 0)
            + np.roll(s, 1, 1) + np.roll(s, -1, 1) - 4.0 * s) / (dx * dx)


def energy(s):
    gx = (np.roll(s, -1, 0) - s) / dx
    gy = (np.roll(s, -1, 1) - s) / dx
    return float(np.sum(0.5 * kappa * (gx * gx + gy * gy) + Vpot(s)) * dx * dx)


def gaussian(A0, cx, cy):
    return A0 * np.exp(-(((X - cx) ** 2 + (Y - cy) ** 2) / (2.0 * w * w)))


def make_cluster(A0, D):
    s = gaussian(A0, 0.0, 0.0)
    for k in range(6):
        ang = k * np.pi / 3.0
        s = s + gaussian(A0, D * np.cos(ang), D * np.sin(ang))
    return s


def half_width_x(s):
    phi = s[:, Nc] ** 2
    if phi[Nc] <= EPS:
        return np.nan
    i = Nc
    while i + 1 < N and phi[i + 1] > EPS:
        i += 1
    if i + 1 >= N:
        return np.nan
    xr = coords[i] + dx * (phi[i] - EPS) / (phi[i] - phi[i + 1])
    j = Nc
    while j - 1 >= 0 and phi[j - 1] > EPS:
        j -= 1
    if j - 1 < 0:
        return np.nan
    xl = coords[j] - dx * (phi[j] - EPS) / (phi[j] - phi[j - 1])
    return 0.5 * (xr - xl)


def gap_x(s):
    phi = s[:, Nc] ** 2
    if phi[Nc] > EPS:
        return 0.0
    i = Nc
    while i - 1 >= 0 and phi[i - 1] <= EPS:
        i -= 1
    if i - 1 < 0:
        return np.nan
    x_l = coords[i - 1] + dx * (phi[i - 1] - EPS) / (phi[i - 1] - phi[i])
    j = Nc
    while j + 1 < N and phi[j + 1] <= EPS:
        j += 1
    if j + 1 >= N:
        return np.nan
    x_r = coords[j + 1] - dx * (phi[j + 1] - EPS) / (phi[j + 1] - phi[j])
    return x_r - x_l


def prep_object(A0, r_target, max_steps=20000):
    s = make_cluster(A0, 3.0)
    H_prev = energy(s)
    H_inc_max = 0.0
    for step in range(1, max_steps + 1):
        s = s + dt * (kappa * lap(s) - Vprime(s))
        if step % REC_EVERY:
            continue
        H = energy(s)
        if H > H_prev:
            H_inc_max = max(H_inc_max, H - H_prev)
        H_prev = H
        R = half_width_x(s)
        if np.isfinite(R) and R >= r_target:
            return s, R, step * dt, H_inc_max
    raise RuntimeError("prep_object: nie osiagnieto r_target")


def roll2(s, ci, cj):
    return np.roll(np.roll(s, ci, 0), cj, 1)


print("=== Phase 1 (op-nbody-additivity): statyczna addytywnosc (LOCK) ===")
print(f"N={N} L={L} dx={dx} kappa={kappa} a={a} b={b} c={c} eps={EPS} "
      f"R_obj={R_OBJ} g_scan={G_SCAN} podloga={FLOOR}")
print()

# ---------------- G1: prep + determinizm ----------------
s_obj, R_act, tau_prep, hinc1 = prep_object(0.6, R_OBJ)
s_obj2, _, _, hinc2 = prep_object(0.6, R_OBJ)
deterministic = bool(np.array_equal(s_obj, s_obj2))
finite_ok = bool(np.all(np.isfinite(s_obj))) and float(np.max(np.abs(s_obj))) <= S_GUARD
H_obj = energy(s_obj)
print(f"obiekt: R_act={R_act:.4f} po tau={tau_prep:.1f}, H_obj={H_obj:.6f}; "
      f"determinizm bitowy prep x2: {deterministic}; "
      f"H_inc_max={max(hinc1, hinc2):.2e}; finite={finite_ok}")
print()

E2_cache = {}


def E2(delta_cells):
    ci, cj = delta_cells
    if ci < 0 or (ci == 0 and cj < 0):
        ci, cj = -ci, -cj
    key = (ci, cj)
    if key not in E2_cache:
        pair = s_obj + roll2(s_obj, ci, cj)
        E2_cache[key] = energy(pair) - 2.0 * H_obj
    return E2_cache[key]


def config_delta(centers_xy):
    cells = [(int(round(cx / dx)), int(round(cy / dx))) for cx, cy in centers_xy]
    s_cfg = np.zeros_like(s_obj)
    for ci, cj in cells:
        s_cfg = s_cfg + roll2(s_obj, ci, cj)
    E_N = energy(s_cfg) - len(cells) * H_obj
    psum = 0.0
    for i in range(len(cells)):
        for j in range(i + 1, len(cells)):
            psum += E2((cells[j][0] - cells[i][0], cells[j][1] - cells[i][1]))
    return E_N, psum, (E_N - psum) / abs(psum)


# ---------------- G2: ciaglosc prawa 2-cialowego ----------------
print("--- G2: pary osiowe, refit Yukawy (g_meas z kompozycji) ---")
g_meas_l, E2_l = [], []
for g in G2_SCAN:
    cell = int(round((R_act + 0.5 * (2 * R_act + g) - R_act) / dx))  # = d/2 zaokr.
    cell = int(round((2 * R_act + g) / (2 * dx)))
    pair = roll2(s_obj, -cell, 0) + roll2(s_obj, +cell, 0)
    gm = gap_x(pair)
    e2 = energy(pair) - 2.0 * H_obj
    g_meas_l.append(gm); E2_l.append(e2)
    print(f"  g0={g:3.1f}: g_meas={gm:7.4f}  E_2={e2:+.6e}")
g_meas_l, E2_l = np.array(g_meas_l), np.array(E2_l)
ok2 = np.isfinite(g_meas_l) & (E2_l < 0)
yv = np.log(-E2_l[ok2])
Amat = np.vstack([np.ones(ok2.sum()), -g_meas_l[ok2]]).T
coef, *_ = np.linalg.lstsq(Amat, yv, rcond=None)
pred = Amat @ coef
rss = float(np.sum((yv - pred) ** 2))
tss = float(np.sum((yv - yv.mean()) ** 2))
r2_g2 = 1.0 - rss / tss
m_g2, A_g2 = coef[1], np.exp(coef[0])
print(f"  fit: A={A_g2:.4f}  m={m_g2:.4f}  R^2(log)={r2_g2:.6f}")
print()

# ---------------- konfiguracje N-cialowe ----------------
print("--- delta_N(g): chain3 / triangle3 / hex6 ---")
res = {}
for g in G_SCAN:
    d = 2 * R_act + g
    rho = d / np.sqrt(3.0)
    chain = [(-d, 0.0), (0.0, 0.0), (d, 0.0)]
    tri = [(rho * np.cos(t), rho * np.sin(t))
           for t in (np.pi / 2, np.pi / 2 + 2 * np.pi / 3, np.pi / 2 + 4 * np.pi / 3)]
    hexa = [(d * np.cos(t), d * np.sin(t))
            for t in (np.pi / 6 + k * np.pi / 3 for k in range(6))]
    out = {}
    for name, cfg in (("chain3", chain), ("tri3", tri), ("hex6", hexa)):
        E_N, psum, delta = config_delta(cfg)
        out[name] = (E_N, psum, delta)
    res[g] = out
    print(f"  g={g:3.1f}:"
          + "".join(f"  {n}: E_N={v[0]:+.5e} sum_par={v[1]:+.5e} "
                    f"delta={v[2]:+.3e}" for n, v in out.items()))
print()

# ---------------- EWALUACJA ----------------
print("=== EWALUACJA (LOCKED) ===")
g1 = deterministic and finite_ok and max(hinc1, hinc2) <= H_TOL * max(1.0, abs(H_obj))
print(f"G1: determinizm={deterministic}, finite={finite_ok}, "
      f"H_nierosnace(prep)={max(hinc1, hinc2) <= H_TOL * max(1.0, abs(H_obj))} "
      f"-> {'PASS' if g1 else 'FAIL'}")

g2 = (0.8 <= m_g2 <= 1.2) and (r2_g2 >= 0.97)
print(f"G2 (ciaglosc): m={m_g2:.4f} in [0.8,1.2]: {0.8 <= m_g2 <= 1.2}; "
      f"R^2(log)={r2_g2:.4f} >= 0.97: {r2_g2 >= 0.97} -> {'PASS' if g2 else 'FAIL'}")

g3_items = [(g, n, res[g][n][2]) for g in G_SCAN if g >= 3.0
            for n in ("chain3", "tri3", "hex6")]
g3 = all(abs(dl) <= 0.02 for _, _, dl in g3_items)
worst = max(g3_items, key=lambda t: abs(t[2]))
print(f"G3 (rezim addytywny, g>=3): max|delta| = {abs(worst[2]):.3e} "
      f"({worst[1]}, g={worst[0]}) <= 0.02 -> {'PASS' if g3 else 'FAIL'}")

tri_g = np.array(G_SCAN)
tri_d = np.array([res[g]["tri3"][2] for g in G_SCAN])
above = np.abs(tri_d) > FLOOR
ga, da = tri_g[above], tri_d[above]
if len(ga) == 0:
    print("G4: NIEROZSTRZYGNIETE (wszystkie |delta_tri| ponizej podlogi 1e-4)")
    g4, g4_txt = None, "NIEROZSTRZYGNIETE"
else:
    sign_ok = bool(np.all(da < 0))
    order = np.argsort(ga)
    mono_ok = bool(np.all(np.diff(np.abs(da[order])) < 0))
    if len(ga) >= 3:
        Am = np.vstack([np.ones(len(ga)), -ga]).T
        cf, *_ = np.linalg.lstsq(Am, np.log(np.abs(da)), rcond=None)
        mu_fit, C_fit = cf[1], np.exp(cf[0])
        mu_ok = 0.4 <= mu_fit <= 1.1
        mu_txt = f"mu={mu_fit:.4f} in [0.4,1.1]: {mu_ok} (predykcja 0.73; model: 0.59)"
    else:
        mu_fit, mu_ok, mu_txt = np.nan, None, "mu: NIEROZSTRZYGALNE (<3 punktow)"
    geom_ok = all(abs(res[g]["chain3"][2]) < abs(res[g]["tri3"][2]) for g in G_SCAN
                  if abs(res[g]["tri3"][2]) > FLOOR)
    g4 = sign_ok and mono_ok and geom_ok and (mu_ok is not False)
    g4_txt = "PASS" if g4 else "FAIL"
    print(f"G4 (RDZEN, mechanizm Fermata): punkty nad podloga: {len(ga)}/{len(tri_g)}")
    print(f"  (a) znak delta_tri < 0: {sign_ok}")
    print(f"  (b) |delta_tri| rosnie przy malejacym g: {mono_ok}")
    print(f"  (c) {mu_txt}")
    print(f"  (d) |delta_chain| < |delta_tri| (nad podloga): {geom_ok}   "
          f"[chain: " + ", ".join(f"{res[g]['chain3'][2]:+.1e}" for g in G_SCAN) + "]")
    print(f"  -> {g4_txt}")

print()
print(f"PODSUMOWANIE PHASE 1: G1={'PASS' if g1 else 'FAIL'} "
      f"G2={'PASS' if g2 else 'FAIL'} G3={'PASS' if g3 else 'FAIL'} G4={g4_txt}")
