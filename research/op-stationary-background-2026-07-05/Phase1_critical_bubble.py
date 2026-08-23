# Phase 1: babel krytyczny (K1) — G2 + czesc G1 (LOCK: Phase0_balance.md)
# Bisekcja R_c (rosnie/kurczy sie) + lambda_min zamrozonego babla ~R_c.

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt = 0.02
N, L = 256, 128.0
dx = L / N
EPS = 0.30
S_GUARD = 10.0
R_C_PRED = 3.0626
W_TANH = 1.5

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)

coords = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
RR = np.hypot(X, Y)


def Vprime(s):
    return s * (a - b_par * np.abs(s) + c * s * s)


def Vpp(s):
    return a - 2.0 * b_par * np.abs(s) + 3.0 * c * s * s


def lap(s):
    return (np.roll(s, 1, 0) + np.roll(s, -1, 0)
            + np.roll(s, 1, 1) + np.roll(s, -1, 1) - 4.0 * s) / (dx * dx)


def bubble(R):
    return s_star * 0.5 * (1.0 - np.tanh((RR - R) / W_TANH))


def classify(R, max_steps=8000):
    s = bubble(R)
    area0 = float(np.mean(s * s > EPS))
    fl = {"nonfinite": False, "runaway": False}
    for stp in range(1, max_steps + 1):
        s = s + dt * (kappa * lap(s) - Vprime(s))
        if stp % 500 == 0:
            if not np.all(np.isfinite(s)):
                fl["nonfinite"] = True
                return "ERR", fl
            if np.max(np.abs(s)) > S_GUARD:
                fl["runaway"] = True
                return "ERR", fl
            Phi = s * s
            if Phi.max() < EPS:
                return "SHRINK", fl
            if float(np.mean(Phi > EPS)) > 1.5 * max(area0, 1e-9):
                return "GROW", fl
    # rozstrzygniecie po trendzie
    return ("GROW" if float(np.mean(s * s > EPS)) > area0 else "SHRINK"), fl


print("=== Phase 1: babel krytyczny K1 (LOCK: Phase0_balance.md) ===")
print(f"N={N} L={L} dx={dx} dt={dt}  R_c_pred={R_C_PRED} "
      f"pasmo=[{0.8*R_C_PRED:.3f}, {1.2*R_C_PRED:.3f}]")
print()
lo, hi = 1.0, 6.0
v_lo, fl1 = classify(lo)
v_hi, fl2 = classify(hi)
print(f"  R={lo}: {v_lo};  R={hi}: {v_hi}   (oczekiwane: SHRINK / GROW)")
flags = [fl1, fl2]
while hi - lo > 0.125:
    mid = 0.5 * (lo + hi)
    v, fl = classify(mid)
    flags.append(fl)
    print(f"  R={mid:.4f}: {v}")
    if v == "SHRINK":
        lo = mid
    else:
        hi = mid
R_c = 0.5 * (lo + hi)
print(f"  R_c = {R_c:.4f} +/- {0.5*(hi-lo):.4f}")
print()

# lambda_min zamrozonego babla ~R_c (Hessian pola rzeczywistego)
s_b = bubble(R_c)
Vb = Vpp(s_b)
rng = np.random.default_rng(20260704)
v = rng.standard_normal((N, N))
v /= np.linalg.norm(v)
SHIFT = 20.0
lam = np.nan
for it in range(4000):
    Lv = -kappa * lap(v) + Vb * v
    u_new = SHIFT * v - Lv
    v = u_new / np.linalg.norm(u_new)
    if it % 500 == 499:
        lam = float(np.sum(v * (-kappa * lap(v) + Vb * v)))
amp2 = v * v
r_mode = float(np.sqrt(np.sum(amp2 * RR * RR)))
print(f"--- lambda_min zamrozonego babla R={R_c:.3f} ---")
print(f"  lambda_min = {lam:+.5f}  "
      + (f"(TACHION/mod wzrostu, gamma={np.sqrt(-lam):.4f}; "
         f"por. B2: gamma sciany=0.123)" if lam < 0 else "(brak modu ujemnego?)"))
print(f"  lokalizacja modu: <r> = {r_mode:.2f} (R_c = {R_c:.2f})")
print()

print("=== EWALUACJA (LOCKED) ===")
g1 = all(not (f["nonfinite"] or f["runaway"]) for f in flags)
in_band = 0.8 * R_C_PRED <= R_c <= 1.2 * R_C_PRED
g2 = in_band and lam < 0 and (hi - lo) <= 0.25
print(f"G1 (czesc Phase1): flagi czyste: {g1}")
print(f"G2 (K1): R_c={R_c:.3f} w pasmie [{0.8*R_C_PRED:.3f}, "
      f"{1.2*R_C_PRED:.3f}]: {in_band}; Delta_R<=0.25: {(hi-lo)<=0.25}; "
      f"lambda_min<0: {lam<0}")
print(f"  -> {'PASS' if g2 else 'FAIL'}: K1 = konfiguracja stacjonarna, "
      f"ale SIODLO (mod wzrostu) -> odrzucony jako tlo (zgodnie z "
      f"oczekiwaniem pre-code)")
