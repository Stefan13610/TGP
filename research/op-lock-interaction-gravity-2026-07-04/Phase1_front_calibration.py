# Phase 1: kalibracja frontu — rownanie #3: v(R) = v_inf - c/R  (G2, elementy G1)
# Implementacja SCISLE wg Phase0_balance.md (LOCKED, sekcje 2-5 + D1-D10).
#
# Dynamika: ds/dtau = kappa*Lap(s) - V'(s), jawny Euler; ZERO clampow/pinningu.

import numpy as np

# ---------------- parametry LOCKED (sekcja 2) ----------------
N = 256
L = 128.0
dx = L / N                      # 0.5
kappa = 0.50
a, b, c = 0.50, 1.60, 1.00
dt = 0.02
EPS = 0.30
S_GUARD = 10.0
w = 1.50                        # profil gaussowski jak w op-bare-substrate-genesis
REC_EVERY = 10                  # D10: kadencja zapisow/proby H
VBASE = 2.0                     # D3: baza roznic centralnych
KOFF = int(round(VBASE / (REC_EVERY * dt)))   # offset rekordow dla roznicy centralnej
H_TOL = 1e-10
FIT_RMIN, FIT_RMAX = 6.0, 25.0  # D9

Nc = N // 2
coords = (np.arange(N) - Nc) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
CHEB_ZONE = np.maximum(np.abs(X), np.abs(Y)) > 0.45 * L   # D6


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


def make_cluster(A0, D):        # K = 7: centrum + 6 heks
    s = gaussian(A0, 0.0, 0.0)
    for k in range(6):
        ang = k * np.pi / 3.0
        s = s + gaussian(A0, D * np.cos(ang), D * np.sin(ang))
    return s


def half_width_x(s):
    """D8: polowa szerokosci Phi>eps wzdluz osi x (y=0), interpolacja subgridowa."""
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


def run_control(s0, max_steps, r_stop):
    """Ewolucja pojedynczego obiektu; zapis (tau, R, H) co REC_EVERY krokow."""
    s = s0.copy()
    taus, Rs, Hs = [0.0], [half_width_x(s)], [energy(s)]
    H_prev, H_inc_max = Hs[0], 0.0
    runaway = nonfinite = contact = False
    for step in range(1, max_steps + 1):
        s = s + dt * (kappa * lap(s) - Vprime(s))
        if step % REC_EVERY:
            continue
        if not np.all(np.isfinite(s)):
            nonfinite = True
            break
        if np.max(np.abs(s)) > S_GUARD:
            runaway = True
            break
        H = energy(s)
        if H > H_prev:
            H_inc_max = max(H_inc_max, H - H_prev)
        H_prev = H
        R = half_width_x(s)
        taus.append(step * dt); Rs.append(R); Hs.append(H)
        if np.any((s * s)[CHEB_ZONE] > EPS):
            contact = True
            break
        if np.isfinite(R) and R >= r_stop:
            break
    return {
        "s_final": s, "tau": np.array(taus), "R": np.array(Rs), "H": np.array(Hs),
        "H_inc_max": H_inc_max,
        "H_monotone": H_inc_max <= H_TOL * max(1.0, abs(Hs[0])),
        "runaway": runaway, "nonfinite": nonfinite, "boundary_contact": contact,
    }


def central_v(tau, val):
    """D3: roznice centralne, baza VBASE."""
    n = len(val)
    idx = np.arange(KOFF, n - KOFF)
    v = (val[idx + KOFF] - val[idx - KOFF]) / (tau[idx + KOFF] - tau[idx - KOFF])
    return idx, v


print("=== Phase 1: kalibracja frontu (LOCK: Phase0_balance.md) ===")
print(f"N={N} L={L} dx={dx} kappa={kappa} a={a} b={b} c={c} dt={dt} w={w} "
      f"eps={EPS} rec_every={REC_EVERY} vbase={VBASE} fit_R=[{FIT_RMIN},{FIT_RMAX}]")
print()

# ---------------- run kontrolny (x2 dla determinizmu, G1) ----------------
s0 = make_cluster(0.6, 3.0)
r1 = run_control(s0, 12000, 26.0)
r2 = run_control(s0, 12000, 26.0)
deterministic = bool(np.array_equal(r1["s_final"], r2["s_final"]))

print(f"run#1: {len(r1['tau'])} rekordow, tau_end={r1['tau'][-1]:.1f}, "
      f"R_end={r1['R'][-1]:.3f}, H: {r1['H'][0]:.4f} -> {r1['H'][-1]:.4f}")
print(f"flagi: nonfinite={r1['nonfinite']} runaway={r1['runaway']} "
      f"boundary_contact={r1['boundary_contact']} H_monotone={r1['H_monotone']} "
      f"(H_inc_max={r1['H_inc_max']:.2e})")
print(f"determinizm bitowy run#1 vs run#2: {deterministic}")
print()

# ---------------- v(R) i fit rownania #3 ----------------
idx, v = central_v(r1["tau"], r1["R"])
Rmid = r1["R"][idx]
ok = np.isfinite(v) & np.isfinite(Rmid) & (Rmid >= FIT_RMIN) & (Rmid <= FIT_RMAX)
Rf, vf = Rmid[ok], v[ok]
A_ = np.vstack([np.ones(len(Rf)), 1.0 / Rf]).T
coef, res_, rank_, sv_ = np.linalg.lstsq(A_, vf, rcond=None)
v_inf, c_fit = coef[0], -coef[1]
pred = A_ @ coef
rss = float(np.sum((vf - pred) ** 2))
tss = float(np.sum((vf - vf.mean()) ** 2))
r2fit = 1.0 - rss / tss

print("--- probki R(tau), v(R) (co 10. punkt fitu) ---")
for k in range(0, len(Rf), 10):
    print(f"  R={Rf[k]:7.3f}   v={vf[k]:.6f}   v_fit={pred[k]:.6f}")
print()
print(f"FIT #3:  v(R) = v_inf - c/R")
print(f"  n_punktow = {len(Rf)}")
print(f"  v_inf = {v_inf:.6f}   (Phase0 check, front 1D: 0.1599)")
print(f"  c     = {c_fit:.6f}   (predykcja LOCK: c ~ kappa = {kappa})")
print(f"  R^2   = {r2fit:.6f}")
print()

# ---------------- ewaluacja ----------------
g1_part = (deterministic and not r1["nonfinite"] and not r1["runaway"]
           and r1["H_monotone"])
g2 = (r2fit >= 0.98) and (0.2 <= c_fit <= 1.0)
print("=== EWALUACJA (LOCKED) ===")
print(f"G1 (czesc Phase1): determinizm={deterministic}, finite/no-runaway="
      f"{not (r1['nonfinite'] or r1['runaway'])}, H_nierosnace={r1['H_monotone']} "
      f"-> {'PASS' if g1_part else 'FAIL'}")
print(f"G2 (kalibracja frontu): R^2={r2fit:.4f} (>=0.98: {r2fit >= 0.98}), "
      f"c={c_fit:.4f} (in [0.2,1.0]: {0.2 <= c_fit <= 1.0}) "
      f"-> {'PASS' if g2 else 'FAIL'}")
print()
print(f"PODSUMOWANIE PHASE 1: G1(czesc)={'PASS' if g1_part else 'FAIL'} "
      f"G2={'PASS' if g2 else 'FAIL'}")
