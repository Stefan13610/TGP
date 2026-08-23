# Phase 1: genesis_u1 — ciaglosc fundamentu na polu zespolonym (G2 + czesc G1)
# Implementacja SCISLE wg Phase0_balance.md (LOCKED, sekcja 2-4 + D10, D11).
#
# Dynamika: dpsi/dtau = kappa*Lap(psi) - psi*(a - b*|psi| + c*|psi|^2)
# (rownowaznie V'(|psi|)*psi/|psi|, regularne w psi=0). Przeplyw gradientowy.
# Scenariusze D11: bare (szum * e^{i*0.7}) -> zanik; cluster -> lock Phi ~ Phi*.
# Dodatkowo: determinizm bitowy (bare x2) i kowariancja fazowa U(1).

import numpy as np

# ---------------- parametry LOCKED ----------------
N = 256
L = 128.0
dx = L / N
kappa = 0.50
a, b, c = 0.50, 1.60, 1.00
dt = 0.02
EPS = 0.30
S_GUARD = 10.0
A_MIN = 4.0 / (N * N)
NOISE_AMP = 0.05
RNG_SEED = 20260704
THETA0 = 0.7
w = 1.50
REC_EVERY = 10
H_TOL = 1e-10
STEPS_DECAY = 6000
STEPS_LOCK_MAX = 30000
R_TARGET = 8.0

disc = np.sqrt(b * b - 4.0 * a * c)
s_star = (b + disc) / (2.0 * c)
Phi_star = s_star ** 2

Nc = N // 2
coords = (np.arange(N) - Nc) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
FAR_ZONE = np.maximum(np.abs(X), np.abs(Y)) > 0.45 * L


def Vpot(m):
    return 0.5 * a * m**2 - (b / 3.0) * m**3 + 0.25 * c * m**4


def energy(psi):
    gx = (np.roll(psi, -1, 0) - psi) / dx
    gy = (np.roll(psi, -1, 1) - psi) / dx
    m = np.abs(psi)
    return float(np.sum(0.5 * kappa * (np.abs(gx)**2 + np.abs(gy)**2)
                        + Vpot(m)) * dx * dx)


def step(psi):
    m = np.abs(psi)
    lap = (np.roll(psi, 1, 0) + np.roll(psi, -1, 0)
           + np.roll(psi, 1, 1) + np.roll(psi, -1, 1) - 4.0 * psi) / (dx * dx)
    return psi + dt * (kappa * lap - psi * (a - b * m + c * m * m))


def gaussian(A0, cx, cy):
    return A0 * np.exp(-(((X - cx) ** 2 + (Y - cy) ** 2) / (2.0 * w * w)))


def make_bare():
    rng = np.random.default_rng(RNG_SEED)
    return rng.uniform(-NOISE_AMP, NOISE_AMP, size=(N, N)) * np.exp(1j * THETA0)


def make_cluster(A0, D):
    s = gaussian(A0, 0.0, 0.0)
    for k in range(6):
        ang = k * np.pi / 3.0
        s = s + gaussian(A0, D * np.cos(ang), D * np.sin(ang))
    return s * np.exp(1j * THETA0)


def half_width_x(psi):
    phi = np.abs(psi[:, Nc]) ** 2
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


def run(psi0, max_steps, label, stop_R=None):
    psi = psi0.copy()
    H0 = energy(psi)
    H_prev, H_inc_max = H0, 0.0
    flags = {"nonfinite": False, "runaway": False, "contact": False}
    rec = []
    stop_step = max_steps
    for stp in range(1, max_steps + 1):
        psi = step(psi)
        if stp % REC_EVERY:
            continue
        if not np.all(np.isfinite(psi)):
            flags["nonfinite"] = True
            stop_step = stp
            break
        if np.max(np.abs(psi)) > S_GUARD:
            flags["runaway"] = True
            stop_step = stp
            break
        H = energy(psi)
        if H > H_prev:
            H_inc_max = max(H_inc_max, H - H_prev)
        H_prev = H
        Phi = np.abs(psi) ** 2
        if np.any(Phi[FAR_ZONE] > EPS):
            flags["contact"] = True
        if stp % 1500 == 0 or stp == max_steps:
            rec.append((stp, float(Phi.max()), float(np.mean(Phi > EPS))))
        if stop_R is not None:
            R = half_width_x(psi)
            if np.isfinite(R) and R >= stop_R:
                stop_step = stp
                rec.append((stp, float(Phi.max()), float(np.mean(Phi > EPS))))
                break
    Phi = np.abs(psi) ** 2
    return {"psi": psi, "H0": H0, "H_inc_max": H_inc_max, "flags": flags,
            "rec": rec, "stop": stop_step,
            "Phi_max": float(Phi.max()), "area": float(np.mean(Phi > EPS)),
            "H_mono": H_inc_max <= H_TOL * max(1.0, abs(H0)), "label": label}


print("=== Phase 1: genesis_u1 (LOCK: Phase0_balance.md, D11) ===")
print(f"N={N} L={L} dx={dx} kappa={kappa} a={a} b={b} c={c} dt={dt} eps={EPS} "
      f"seed={RNG_SEED} theta0={THETA0}")
print()

# ---------------- bare: zanik + determinizm ----------------
r1 = run(make_bare(), STEPS_DECAY, "bare#1")
r2 = run(make_bare(), STEPS_DECAY, "bare#2")
det = bool(np.array_equal(r1["psi"], r2["psi"]))
print("--- bare (szum U(-0.05,0.05) * e^{i*0.7}), 6000 krokow ---")
for stp, pm, ar in r1["rec"]:
    print(f"  tau={stp*dt:7.1f}  Phi_max={pm:.3e}  metric_area={ar:.6f}")
print(f"  determinizm bitowy (run x2): {det}")
print(f"  flagi: nf={r1['flags']['nonfinite']} run={r1['flags']['runaway']} "
      f"bc={r1['flags']['contact']}  H_mono={r1['H_mono']} "
      f"(H_inc_max={r1['H_inc_max']:.2e})")
bare_dead = r1["area"] < A_MIN
print(f"  zanik: metric_area_final={r1['area']:.2e} < A_MIN={A_MIN:.2e}: {bare_dead}")
print()

# ---------------- kowariancja fazowa U(1) ----------------
rng = np.random.default_rng(RNG_SEED)
s_real = rng.uniform(-NOISE_AMP, NOISE_AMP, size=(N, N)).astype(np.complex128)
rr = run(s_real, STEPS_DECAY, "bare-real")
cov = float(np.max(np.abs(np.abs(r1["psi"]) - np.abs(rr["psi"]))))
print(f"--- kowariancja fazowa U(1): max| |psi|_theta0 - |psi|_real | = {cov:.3e} ---")
print()

# ---------------- cluster: lock ----------------
rc = run(make_cluster(0.6, 3.0), STEPS_LOCK_MAX, "cluster", stop_R=R_TARGET)
print("--- cluster(A0=0.6, D=3, K=7) * e^{i*0.7} ---")
for stp, pm, ar in rc["rec"]:
    print(f"  tau={stp*dt:7.1f}  Phi_max={pm:.5f}  metric_area={ar:.6f}")
Rfin = half_width_x(rc["psi"])
print(f"  R(Phi>eps) osiagniete: {Rfin:.3f} po tau={rc['stop']*dt:.1f}")
print(f"  flagi: nf={rc['flags']['nonfinite']} run={rc['flags']['runaway']} "
      f"bc={rc['flags']['contact']}  H_mono={rc['H_mono']} "
      f"(H_inc_max={rc['H_inc_max']:.2e})")
lock_ok = 0.8 * Phi_star <= rc["Phi_max"] <= 1.2 * Phi_star
print(f"  lock: Phi_max={rc['Phi_max']:.5f} in [0.8, 1.2]*Phi*="
      f"[{0.8*Phi_star:.4f}, {1.2*Phi_star:.4f}]: {lock_ok}")
theta_c = float(np.angle(rc["psi"][Nc, Nc]))
print(f"  faza w centrum po locku: {theta_c:.6f} (theta0={THETA0}) — "
      f"faza globalna zachowana (U(1))")
print()

# ---------------- EWALUACJA ----------------
print("=== EWALUACJA (LOCKED) ===")
g1_ok = det and all(not (r["flags"]["nonfinite"] or r["flags"]["runaway"])
                    and r["H_mono"] for r in (r1, r2, rr, rc))
print(f"G1 (czesc Phase1): determinizm + finite + H_Gamma nierosnace: "
      f"{'PASS' if g1_ok else 'FAIL'}")
g2 = bare_dead and lock_ok
print(f"G2 (ciaglosc fundamentu na U(1)): bare zanika: {bare_dead}, "
      f"cluster lockuje z Phi ~ Phi*: {lock_ok}  -> {'PASS' if g2 else 'FAIL'}")
print()
print(f"PODSUMOWANIE PHASE 1: G1(czesc)={'PASS' if g1_ok else 'FAIL'} "
      f"G2={'PASS' if g2 else 'FAIL'}")
