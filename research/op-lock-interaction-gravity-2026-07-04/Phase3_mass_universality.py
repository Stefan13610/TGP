# Phase 3: skalowanie masowe (G5) i uniwersalnosc stanu (G6)
# Implementacja SCISLE wg Phase0_balance.md (LOCKED, sekcje 3-5 + D1-D10).
#
# G5: pary {8x8, 12x12, 16x16} przy g0=4 -> ratio_mass = Dv(16)/Dv(8)
# G6: dwa obiekty R=8 przygotowane ROZNA historia (A0=0.6 dluzej vs A0=1.0
#     krocej) -> ratio_hist = Dv(hist2)/Dv(hist1); |ratio-1| <= 0.10
# Dv_bar (D5): srednia Delta_v w g in {2.0,2.5,3.0,3.5} na CZESCI WSPOLNEJ
# nosnikow porownywanych runow (raportowane, ktore punkty weszly).

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
VBASE = 2.0
KOFF = int(round(VBASE / (REC_EVERY * dt)))
H_TOL = 1e-10
G0 = 4.0
G_GRID = [2.0, 2.5, 3.0, 3.5]

Nc = N // 2
coords = (np.arange(N) - Nc) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
CHEB_ZONE = np.maximum(np.abs(X), np.abs(Y)) > 0.45 * L


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


def prep_multi(A0, r_targets, max_steps=30000):
    """Ewolucja cluster(A0, D=3, K=7); snapshoty przy kolejnych R >= target."""
    s = make_cluster(A0, 3.0)
    H_prev = energy(s)
    H_inc_max = 0.0
    snaps = {}
    todo = sorted(r_targets)
    for step in range(1, max_steps + 1):
        s = s + dt * (kappa * lap(s) - Vprime(s))
        if step % REC_EVERY:
            continue
        H = energy(s)
        if H > H_prev:
            H_inc_max = max(H_inc_max, H - H_prev)
        H_prev = H
        R = half_width_x(s)
        if todo and np.isfinite(R) and R >= todo[0]:
            snaps[todo.pop(0)] = (s.copy(), R, step * dt)
        if not todo:
            return snaps, H_inc_max
    raise RuntimeError("prep_multi: nie osiagnieto celow")


def evolve_recorded(s0, max_steps, measure):
    s = s0.copy()
    taus, vals = [0.0], [measure(s)]
    H_prev = energy(s)
    H_inc_max = 0.0
    flags = {"nonfinite": False, "runaway": False, "contact": False}
    for step in range(1, max_steps + 1):
        s = s + dt * (kappa * lap(s) - Vprime(s))
        if step % REC_EVERY:
            continue
        if not np.all(np.isfinite(s)):
            flags["nonfinite"] = True
            break
        if np.max(np.abs(s)) > S_GUARD:
            flags["runaway"] = True
            break
        H = energy(s)
        if H > H_prev:
            H_inc_max = max(H_inc_max, H - H_prev)
        H_prev = H
        m = measure(s)
        taus.append(step * dt); vals.append(m)
        if np.any((s * s)[CHEB_ZONE] > EPS):
            flags["contact"] = True
            break
        if not np.isfinite(m) or m < 1.0:
            break
    flags["H_inc_max"] = H_inc_max
    flags["H_monotone"] = H_inc_max <= H_TOL * max(1.0, abs(energy(s0)))
    return np.array(taus), np.array(vals), flags


def central_v(tau, val):
    n = len(val)
    idx = np.arange(KOFF, n - KOFF)
    v = (val[idx + KOFF] - val[idx - KOFF]) / (tau[idx + KOFF] - tau[idx - KOFF])
    return idx, v


def pair_deltav(s_obj, R_act, label):
    """Pair przy g0=4 + kontrola D4; zwraca funkcje interp Delta_v(g)."""
    cells = int(round((R_act + 0.5 * G0) / dx))
    pair0 = np.roll(s_obj, -cells, axis=0) + np.roll(s_obj, +cells, axis=0)
    g_init = gap_x(pair0)
    tau_c, R_c_, fl_c = evolve_recorded(s_obj, 2500, half_width_x)
    idx_c, v_c = central_v(tau_c, R_c_)
    tau_p, g_p, fl_p = evolve_recorded(pair0, 2500, gap_x)
    idx_p, vg = central_v(tau_p, g_p)
    npts = min(len(idx_p), len(v_c))
    dv = -vg[:npts] - 2.0 * v_c[:npts]
    gm = g_p[idx_p][:npts]
    gs, dvs = gm[::-1], dv[::-1]
    cover = (gs[0], gs[-1]) if len(gs) else (np.nan, np.nan)
    vals = {gq: (float(np.interp(gq, gs, dvs)) if gs[0] <= gq <= gs[-1] else None)
            for gq in G_GRID}
    print(f"  {label:<14} g_init={g_init:6.3f}  nosnik g=[{cover[0]:.3f},"
          f"{cover[1]:.3f}]  " + "  ".join(
              f"Dv({gq})={vals[gq]:+.5f}" if vals[gq] is not None else f"Dv({gq})=--"
              for gq in G_GRID)
          + f"  [nf={fl_p['nonfinite']} run={fl_p['runaway']} "
            f"bc={fl_p['contact']} H_mono={fl_p['H_monotone'] and fl_c['H_monotone']}]")
    return vals, (fl_c, fl_p)


def dv_bar(vals_list):
    """Srednia po czesci wspolnej pokrytych punktow G_GRID (D5)."""
    common = [gq for gq in G_GRID if all(v[gq] is not None for v in vals_list)]
    if not common:
        return [None] * len(vals_list), common
    return [float(np.mean([v[gq] for gq in common])) for v in vals_list], common


print("=== Phase 3: skalowanie masowe + uniwersalnosc (LOCK: Phase0_balance.md) ===")
print(f"N={N} L={L} dx={dx} kappa={kappa} a={a} b={b} c={c} dt={dt} eps={EPS} "
      f"g0={G0} G_GRID={G_GRID}")
print()

# ---------------- przygotowania ----------------
snaps_h1, hinc1 = prep_multi(0.6, [8.0, 12.0, 16.0])
snaps_h2, hinc2 = prep_multi(1.0, [8.0])
print("--- przygotowanie obiektow ---")
for R0, (s_, Ra, tp) in sorted(snaps_h1.items()):
    print(f"  hist1 (A0=0.6): R_target={R0:4.1f} -> R_act={Ra:.4f} po tau={tp:6.1f}")
for R0, (s_, Ra, tp) in sorted(snaps_h2.items()):
    print(f"  hist2 (A0=1.0): R_target={R0:4.1f} -> R_act={Ra:.4f} po tau={tp:6.1f}")
print(f"  H_inc_max prep: hist1={hinc1:.2e}, hist2={hinc2:.2e}")
print()

# ---------------- pary ----------------
print("--- pary przy g0=4 (kontrola D4 per snapshot) ---")
v8, f8 = pair_deltav(snaps_h1[8.0][0], snaps_h1[8.0][1], "8x8 (hist1)")
v12, f12 = pair_deltav(snaps_h1[12.0][0], snaps_h1[12.0][1], "12x12 (hist1)")
v16, f16 = pair_deltav(snaps_h1[16.0][0], snaps_h1[16.0][1], "16x16 (hist1)")
v8b, f8b = pair_deltav(snaps_h2[8.0][0], snaps_h2[8.0][1], "8x8 (hist2)")
print()

all_flags = [f for pair in (f8, f12, f16, f8b) for f in pair]
g1_ok = all(not (f["nonfinite"] or f["runaway"]) and f["H_monotone"]
            for f in all_flags) and hinc1 == 0.0 and hinc2 == 0.0

# ---------------- G5 ----------------
(b8, b12, b16), common5 = dv_bar([v8, v12, v16])
print("=== EWALUACJA (LOCKED) ===")
print(f"G1 (czesc Phase3): finite/no-runaway + H_nierosnace: "
      f"{'PASS' if g1_ok else 'FAIL'}")
if b8 and b16:
    ratio_mass = b16 / b8
    ratio_12 = b12 / b8
    cls = ("MASOWE (bulk, gravity-like)" if ratio_mass >= 2.0
           else "FRONTOWE (powierzchniowe)" if ratio_mass <= 1.5
           else "MIXED")
    print(f"G5 (dyskryminator; punkty wspolne g={common5}):")
    print(f"  Dv_bar: 8x8={b8:+.5f}  12x12={b12:+.5f}  16x16={b16:+.5f}")
    print(f"  ratio_mass = Dv(16x16)/Dv(8x8) = {ratio_mass:.4f}  "
          f"(12x12/8x8 = {ratio_12:.4f})")
    print(f"  klasyfikacja: {cls}")
else:
    ratio_mass, cls = np.nan, "NIEROZSTRZYGNIETE (brak wspolnego nosnika)"
    print(f"G5: {cls}")

# ---------------- G6 ----------------
(u1, u2), common6 = dv_bar([v8, v8b])
if u1 and u2:
    ratio_hist = u2 / u1
    g6 = abs(ratio_hist - 1.0) <= 0.10
    print(f"G6 (uniwersalnosc stanu; punkty wspolne g={common6}):")
    print(f"  Dv_bar: hist1={u1:+.5f}  hist2={u2:+.5f}  "
          f"ratio_hist={ratio_hist:.4f}  |ratio-1|={abs(ratio_hist-1):.4f} "
          f"(<=0.10) -> {'PASS' if g6 else 'FAIL'}")
else:
    g6, ratio_hist = False, np.nan
    print("G6: FAIL (brak wspolnego nosnika)")

print()
print(f"PODSUMOWANIE PHASE 3: G1(czesc)={'PASS' if g1_ok else 'FAIL'} "
      f"G5={cls} G6={'PASS' if g6 else 'FAIL'}")
