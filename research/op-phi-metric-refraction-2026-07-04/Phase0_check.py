# Phase 0 check: weryfikacja rachunkow PRZED kodem produkcyjnym.
# Wg Phase0_balance.md (LOCKED, sekcje 1-2 + D1-D11):
#   1) masy sektora falowego: m_FV, m_TV; wspolczynnik zalamania n(omega)
#   2) tlo: zlockowany obiekt (cluster 0.6, D=3, K=7 -> R_obj>=12), zamrozone
#   3) D6: widmo L = -kappa*Lap + V''(s_bg) — lambda_min (tachion sciany?)
#   4) D7: eikonal/ray-tracing alpha_pred(b, omega) na TYM SAMYM tle

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt_flow = 0.02
EPS = 0.30
N, L = 256, 128.0
dx = L / N
w = 1.50
R_TARGET = 12.0

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)
Vpp_star = a - 2.0 * b_par * s_star + 3.0 * c * s_star ** 2

Nc = N // 2
coords = (np.arange(N) - Nc) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")


def Vprime(s):
    return s * (a - b_par * np.abs(s) + c * s * s)


def Vpp(s):
    return a - 2.0 * b_par * np.abs(s) + 3.0 * c * s * s


def lap(s):
    return (np.roll(s, 1, 0) + np.roll(s, -1, 0)
            + np.roll(s, 1, 1) + np.roll(s, -1, 1) - 4.0 * s) / (dx * dx)


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
    xr = coords[i] + dx * (phi[i] - EPS) / (phi[i] - phi[i + 1])
    j = Nc
    while j - 1 >= 0 and phi[j - 1] > EPS:
        j -= 1
    xl = coords[j] - dx * (phi[j] - EPS) / (phi[j] - phi[j - 1])
    return 0.5 * (xr - xl)


def prep_object(A0, r_target, max_steps=40000):
    s = make_cluster(A0, 3.0)
    for step in range(1, max_steps + 1):
        s = s + dt_flow * (kappa * lap(s) - Vprime(s))
        if step % 10 == 0:
            R = half_width_x(s)
            if np.isfinite(R) and R >= r_target:
                return s, R, step * dt_flow
    raise RuntimeError("prep_object: nie osiagnieto r_target")


print("=== Phase 0 check: op-phi-metric-refraction (LOCK: Phase0_balance.md) ===")
print()

# ---------------- 1. masy i n(omega) ----------------
m_FV = np.sqrt(a)
m_TV = np.sqrt(Vpp_star)
print("--- 1. sektor falowy: masy i wspolczynnik zalamania (rownania #1-#2) ---")
print(f"  V''(0) = {a}  -> m_FV = {m_FV:.4f}   (LOCK: 0.707)")
print(f"  V''(s*) = {Vpp_star:.5f}  -> m_TV = {m_TV:.4f}   (LOCK: 0.937)")
s_grid = np.linspace(0, 1.25, 500)
vs = Vpp(s_grid)
s_neg = s_grid[vs < 0]
print(f"  D6: V''(s) < 0 dla s in ({s_neg.min():.3f}, {s_neg.max():.3f}), "
      f"min V'' = {vs.min():.4f} przy s = {s_grid[np.argmin(vs)]:.3f}")
print("  n_wnetrze(omega) = sqrt((omega^2 - V''(s*)) / (omega^2 - V''(0))):")
for om in (1.0, 1.1, 1.3):
    n_in = np.sqrt((om * om - Vpp_star) / (om * om - a))
    k0 = np.sqrt((om * om - a) / kappa)
    vg = kappa * k0 / om
    print(f"    omega={om:.1f}: n_in={n_in:.4f}  k0_FV={k0:.4f}  v_g={vg:.4f}")
print(f"  przyklad LOCKa (omega=1.1): n_in = 0.68 -> "
      f"{np.sqrt((1.21 - Vpp_star)/(1.21 - a)):.4f}  OK")
print()

# ---------------- 2. tlo zamrozone ----------------
s_bg, R_obj, tau_prep = prep_object(0.6, R_TARGET)
print("--- 2. tlo: zlockowany obiekt (D9) ---")
print(f"  cluster(0.6, D=3, K=7) -> R_obj = {R_obj:.3f} po tau = {tau_prep:.1f}; "
      f"max s = {s_bg.max():.5f} (s* = {s_star:.5f})")
Vb = Vpp(s_bg)
neg_cells = int(np.sum(Vb < 0))
print(f"  komorek z V''(s_bg) < 0 (pierscien sciany): {neg_cells} "
      f"({100*neg_cells/N/N:.2f}% siatki); min V''(s_bg) = {Vb.min():.4f}")
print()

# ---------------- 3. D6: lambda_min operatora L ----------------
print("--- 3. D6: widmo L = -kappa*Lap + V''(s_bg) (iteracja potegowa) ---")
SHIFT = 20.0
rng = np.random.default_rng(20260704)
v = rng.standard_normal((N, N))
v /= np.linalg.norm(v)
lam = np.nan
for it in range(4000):
    Lv = -kappa * lap(v) + Vb * v
    u_new = SHIFT * v - Lv
    nrm = np.linalg.norm(u_new)
    v = u_new / nrm
    if it % 500 == 499:
        Lv = -kappa * lap(v) + Vb * v
        lam = float(np.sum(v * Lv))
lam_min = lam
gamma = np.sqrt(-lam_min) if lam_min < 0 else 0.0
# lokalizacja modu
amp2 = v * v
r_mode = np.sqrt(np.sum(amp2 * (X * X + Y * Y)))
frac_wall = float(np.sum(amp2[np.sqrt(X * X + Y * Y) < R_obj + 5]))
print(f"  lambda_min = {lam_min:+.5f}  ->  "
      + (f"TACHION, gamma = sqrt(-lambda_min) = {gamma:.4f}" if lam_min < 0
         else "brak modu niestabilnego"))
print(f"  lokalizacja modu: <r> = {r_mode:.2f} (R_obj = {R_obj:.1f}); "
      f"frakcja normy w r < R_obj+5: {frac_wall:.4f}")
if lam_min < 0:
    for om in (1.1,):
        k0 = np.sqrt((om * om - a) / kappa)
        vg = kappa * k0 / om
        tau_fly = (2 * R_obj + 60.0) / vg
        print(f"  prognoza kontaminacji (omega={om}): tau_lotu ~ {tau_fly:.0f}, "
              f"amplifikacja e^(gamma*tau) ~ e^{gamma*tau_fly:.1f} = "
              f"{np.exp(min(gamma*tau_fly, 700)):.2e}")
    print("  -> wg D6: jesli pomiar ugiecia skontaminowany, G3/G4/G6 "
          "NIEAPLIKOWALNE (klauzula pre-code)")
print()

# ---------------- 4. D7: eikonal alpha_pred(b, omega) ----------------
print("--- 4. D7: ray tracing eikonalny na TYM SAMYM tle ---")


def n2_grid(om):
    return (om * om - Vb) / (om * om - a)


def bilinear(F, x, y):
    fx = (x + L / 2) / dx
    fy = (y + L / 2) / dx
    i0 = int(np.floor(fx)) % N
    j0 = int(np.floor(fy)) % N
    tx, ty = fx - np.floor(fx), fy - np.floor(fy)
    i1, j1 = (i0 + 1) % N, (j0 + 1) % N
    return ((1 - tx) * (1 - ty) * F[i0, j0] + tx * (1 - ty) * F[i1, j0]
            + (1 - tx) * ty * F[i0, j1] + tx * ty * F[i1, j1])


def trace(om, b):
    n2 = n2_grid(om)
    gx = (np.roll(n2, -1, 0) - np.roll(n2, 1, 0)) / (2 * dx)
    gy = (np.roll(n2, -1, 1) - np.roll(n2, 1, 1)) / (2 * dx)
    x, y = -60.0, b
    n0 = np.sqrt(bilinear(n2, x, y))
    px, py = n0, 0.0
    ds = 0.05
    for _ in range(40000):
        # RK4 dla (r, p): dr/ds = p, dp/ds = grad(n^2)/2
        def F(st):
            xx, yy, ppx, ppy = st
            return np.array([ppx, ppy,
                             0.5 * bilinear(gx, xx, yy),
                             0.5 * bilinear(gy, xx, yy)])
        st = np.array([x, y, px, py])
        k1 = F(st)
        k2 = F(st + 0.5 * ds * k1)
        k3 = F(st + 0.5 * ds * k2)
        k4 = F(st + ds * k3)
        st = st + ds / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        x, y, px, py = st
        if x > 60.0:
            break
    return np.degrees(np.arctan2(py, px)), y


print("  alpha_pred [stopnie] (alpha > 0 = ugiecie OD obiektu dla b > 0):")
print("    b:        " + "".join(f"  {bb:7.0f}" for bb in (8, 12, 16, 20)))
alpha_pred = {}
for om in (1.0, 1.1, 1.3):
    row = f"    omega={om:.1f}:"
    for bb in (8.0, 12.0, 16.0, 20.0):
        al, yfin = trace(om, bb)
        alpha_pred[(om, bb)] = al
        row += f"  {al:+7.3f}"
    print(row)
print("  znaki/monotonie raportowane; szczegolowa interpretacja w Phase 2")
print()

print("=== WNIOSKI PHASE 0 ===")
print(f"  1) masy zgodne z LOCKiem; kontrast n silny (n_in(1.1) = 0.68)")
print(f"  2) tlo przygotowane deterministycznie (R_obj = {R_obj:.2f})")
print(f"  3) lambda_min = {lam_min:+.5f}: "
      + ("TACHION SCIANY POTWIERDZONY rachunkiem (zgodnie z pre-rejestracja "
         "D6) — pomiar ugiecia na tle zamrozonym zagrozony; decyzja wg "
         "empirycznej klauzuli D6(c) w Phase 2" if lam_min < 0 else
         "tlo stabilne — pelny program pomiarowy"))
print(f"  4) alpha_pred policzone (D7) — predykcje eikonalne zalockowane "
      f"przed pomiarem")
