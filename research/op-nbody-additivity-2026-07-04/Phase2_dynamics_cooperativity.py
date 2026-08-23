# Phase 2: dynamika — klasyfikator kooperatywnosci G5 + obserwacja hex6
# Implementacja SCISLE wg Phase0_balance.md (LOCKED) tego opa.
# ratio_coop = Delta_v_bar(szczelina chain3, g0=4) / Delta_v_bar(para, g0=4)
# kontrola D4: pojedynczy obiekt z tego samego snapshotu.

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
R_OBJ = 8.0
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


def gap_1d(phi, k0, ax):
    """Szczelina wokol indeksu k0 w 1D profilu phi (os ax tylko do skali)."""
    if phi[k0] > EPS:
        return 0.0
    i = k0
    while i - 1 >= 0 and phi[i - 1] <= EPS:
        i -= 1
    if i - 1 < 0:
        return np.nan
    x_l = (i - 1) * dx + dx * (phi[i - 1] - EPS) / (phi[i - 1] - phi[i])
    j = k0
    while j + 1 < len(phi) and phi[j + 1] <= EPS:
        j += 1
    if j + 1 >= len(phi):
        return np.nan
    x_r = (j + 1) * dx - dx * (phi[j + 1] - EPS) / (phi[j + 1] - phi[j])
    return x_r - x_l


def prep_object(A0, r_target, max_steps=20000):
    s = make_cluster(A0, 3.0)
    H_prev = energy(s)
    for step in range(1, max_steps + 1):
        s = s + dt * (kappa * lap(s) - Vprime(s))
        if step % REC_EVERY:
            continue
        H_prev = energy(s)
        R = half_width_x(s)
        if np.isfinite(R) and R >= r_target:
            return s, R
    raise RuntimeError("prep")


def roll2(s, ci, cj):
    return np.roll(np.roll(s, ci, 0), cj, 1)


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


def deltav_interp(tau_p, g_p, v_c, label):
    idx_p, vg = central_v(tau_p, g_p)
    npts = min(len(idx_p), len(v_c))
    dv = -vg[:npts] - 2.0 * v_c[:npts]
    gm = g_p[idx_p][:npts]
    gs, dvs = gm[::-1], dv[::-1]
    vals = {gq: (float(np.interp(gq, gs, dvs)) if len(gs) and gs[0] <= gq <= gs[-1]
                 else None) for gq in G_GRID}
    print(f"  {label:<16} nosnik g=[{gs[0]:.3f},{gs[-1]:.3f}]  "
          + "  ".join(f"Dv({gq})={vals[gq]:+.5f}" if vals[gq] is not None
                      else f"Dv({gq})=--" for gq in G_GRID))
    return vals


print("=== Phase 2 (op-nbody-additivity): dynamika kooperatywnosci (LOCK) ===")
print(f"N={N} L={L} dx={dx} dt={dt} eps={EPS} g0={G0} G_GRID={G_GRID}")
print()

s_obj, R_act = prep_object(0.6, R_OBJ)
print(f"obiekt: R_act={R_act:.4f}")
d4 = 2 * R_act + G0
half = int(round(d4 / (2 * dx)))
D = int(round(d4 / dx))

# kontrola D4
tau_c, R_c_, fl_c = evolve_recorded(s_obj, 2500, half_width_x)
idx_c, v_c = central_v(tau_c, R_c_)
print(f"kontrola: R {R_c_[0]:.3f} -> {R_c_[-1]:.3f}  "
      f"[nf={fl_c['nonfinite']} run={fl_c['runaway']} bc={fl_c['contact']} "
      f"H_mono={fl_c['H_monotone']}]")
print()

print("--- pary/lancuch przy g0=4 ---")
pair0 = roll2(s_obj, -half, 0) + roll2(s_obj, +half, 0)
tau_p, g_p, fl_p = evolve_recorded(pair0, 2500, lambda s: gap_1d(s[:, Nc] ** 2, Nc, 0))
v_pair = deltav_interp(tau_p, g_p, v_c, "pair(g0=4)")

chain0 = roll2(s_obj, -D, 0) + s_obj + roll2(s_obj, +D, 0)
k_gap = Nc - D // 2                     # srodek lewej szczeliny
tau_h, g_h, fl_h = evolve_recorded(chain0, 2500,
                                   lambda s: gap_1d(s[:, Nc] ** 2, k_gap, 0))
v_chain = deltav_interp(tau_h, g_h, v_c, "chain3(g0=4)")
print()

# ---------------- obserwacja: hex6 przy g=3 vs para przy g=3 ----------------
print("--- OBSERWACJA (poza kryteriami): hex6(g=3) vs pair(g=3) ---")
g_obs = 3.0
d3 = 2 * R_act + g_obs
cells_hex = [(int(round(d3 * np.cos(t) / dx)), int(round(d3 * np.sin(t) / dx)))
             for t in (np.pi / 6 + k * np.pi / 3 for k in range(6))]
hex0 = np.zeros_like(s_obj)
for ci, cj in cells_hex:
    hex0 = hex0 + roll2(s_obj, ci, cj)
i_col = Nc + cells_hex[0][0]            # pionowa krawedz: wierzcholki -30/+30 st.
tau_x, g_x, fl_x = evolve_recorded(hex0, 2500,
                                   lambda s: gap_1d(s[i_col, :] ** 2, Nc, 1))
half3 = int(round(d3 / (2 * dx)))
pair3 = roll2(s_obj, -half3, 0) + roll2(s_obj, +half3, 0)
tau_q, g_q, fl_q = evolve_recorded(pair3, 2500,
                                   lambda s: gap_1d(s[:, Nc] ** 2, Nc, 0))


def t_reach(tau, g, gt):
    m = np.isfinite(g)
    tt, gg = tau[m], g[m]
    below = np.where(gg <= gt)[0]
    return float(tt[below[0]]) if len(below) else np.nan


print(f"  hex6:  g_init={g_x[0]:.3f} -> koniec g={g_x[-1]:.3f} przy tau={tau_x[-1]:.1f}; "
      f"tau(g<=1.5)={t_reach(tau_x, g_x, 1.5):.2f}  "
      f"[nf={fl_x['nonfinite']} run={fl_x['runaway']} bc={fl_x['contact']} "
      f"H_mono={fl_x['H_monotone']}]")
print(f"  pair3: g_init={g_q[0]:.3f} -> koniec g={g_q[-1]:.3f} przy tau={tau_q[-1]:.1f}; "
      f"tau(g<=1.5)={t_reach(tau_q, g_q, 1.5):.2f}  "
      f"[nf={fl_q['nonfinite']} run={fl_q['runaway']} bc={fl_q['contact']} "
      f"H_mono={fl_q['H_monotone']}]")
rr = t_reach(tau_x, g_x, 1.5) / t_reach(tau_q, g_q, 1.5)
print(f"  stosunek czasow domkniecia hex/para = {rr:.4f} (addytywnie ~1; "
      f"<1 = kooperatywne przyspieszenie)")
print()

# ---------------- EWALUACJA ----------------
print("=== EWALUACJA (LOCKED) ===")
all_fl = [fl_c, fl_p, fl_h, fl_x, fl_q]
g1 = all(not (f["nonfinite"] or f["runaway"]) and f["H_monotone"] for f in all_fl)
print(f"G1 (czesc Phase2): finite/no-runaway + H_nierosnace: "
      f"{'PASS' if g1 else 'FAIL'}")

common = [gq for gq in G_GRID if v_pair[gq] is not None and v_chain[gq] is not None]
if common:
    bp = float(np.mean([v_pair[gq] for gq in common]))
    bh = float(np.mean([v_chain[gq] for gq in common]))
    ratio = bh / bp
    cls = ("LOKALNY (sciana-ogon, addytywny)" if 0.9 <= ratio <= 1.1
           else "KOLEKTYWNY (wzmocnienie ponadparowe)" if ratio >= 1.3
           else "MIXED")
    print(f"G5 (klasyfikator; punkty wspolne g={common}):")
    print(f"  Dv_bar: para={bp:+.5f}  chain3={bh:+.5f}  ratio_coop={ratio:.4f}")
    print(f"  klasyfikacja: {cls}")
else:
    ratio, cls = np.nan, "NIEROZSTRZYGNIETE (brak wspolnego nosnika)"
    print(f"G5: {cls}")

print()
print(f"PODSUMOWANIE PHASE 2: G1(czesc)={'PASS' if g1 else 'FAIL'}  G5={cls}")
