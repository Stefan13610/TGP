# Phase 0 check: rachunki PRZED kodem produkcyjnym (LOCK: Phase0_balance.md)
#   1) R_c_pred z prawa frontu B (cytowane v_inf, c)
#   2) profil wiru (ODE) -> n^2(r) dla fal amplitudowych; ogon A_t/r^2
#   3) ray tracing: alpha_pred(b, omega) — zaklad na op pomiarowy
#   4) progi: lambda_thr widmowe i gamma empiryczne

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
V_INF, C_FRONT = 0.15993, 0.4898          # op-lock-interaction, G2 (CYTOWANE)

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)
Vpp_star = a - 2.0 * b_par * s_star + 3.0 * c * s_star ** 2
Vppp_star = -2.0 * b_par + 6.0 * c * s_star
xi2 = kappa / Vpp_star
m_TV = np.sqrt(Vpp_star)


def Vpp(m):
    return a - 2.0 * b_par * np.abs(m) + 3.0 * c * m * m


def Vprime(m):
    return m * (a - b_par * np.abs(m) + c * m * m)


print("=== Phase 0 check: op-stationary-background (LOCK: Phase0_balance.md) ===")
print()
print("--- 1. babel krytyczny (rownanie #1) ---")
Rc = C_FRONT / V_INF
print(f"  v_inf={V_INF} c={C_FRONT} (op-lock-interaction, G2)  ->  "
      f"R_c_pred = {Rc:.4f}; pasmo G2: [{0.8*Rc:.3f}, {1.2*Rc:.3f}]")
print()

print("--- 2. profil wiru i n^2(r) (rownania #2-#3) ---")
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
A_t = Vppp_star * s_star * xi2
print(f"  V'''(s*) = {Vppp_star:.4f}   A_t = V'''(s*)*s**xi^2 = {A_t:.4f}")
print(f"  m_TV = {m_TV:.4f}; pasmo V''(f)<0: f in (0.190, 0.874) "
      f"-> pierscien r in ~({np.interp(0.190, f, r1d):.2f}, "
      f"{np.interp(0.874, f, r1d):.2f})  [zrodlo oczekiwanego modu L_toy]")
print("  n^2(r) dla omega (kolumny 1.0/1.1/1.3) + ogon algebraiczny:")
for rq in (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0):
    fq = float(np.interp(rq, r1d, f))
    row = f"   r={rq:5.1f}: "
    for om in (1.0, 1.1, 1.3):
        n2 = (om * om - Vpp(fq)) / (om * om - Vpp_star)
        row += f"  {n2:7.4f}"
    tail = 1.0 + A_t / (rq * rq * (1.21 - Vpp_star))
    row += f"   [ogon(1.1): {tail:.4f}]"
    print(row)
print()

print("--- 3. ray tracing: alpha_pred(b, omega) [deg] ---")
print("    (alpha > 0 = ugiecie KU wirowi; profil radialny z ODE)")
n2r = {}
for om in (1.0, 1.1, 1.3):
    n2r[om] = (om * om - Vpp(f)) / (om * om - Vpp_star)


def trace(om, b):
    n2v = n2r[om]
    dn2 = np.gradient(n2v, dr)

    def n2_at(rr):
        return np.interp(rr, r1d, n2v, left=n2v[0], right=1.0)

    def dn2_at(rr):
        return np.interp(rr, r1d, dn2, left=0.0, right=0.0)

    x, y = -60.0, b
    px, py = np.sqrt(n2_at(np.hypot(x, y))), 0.0
    ds = 0.02
    for _ in range(300000):
        def F(st):
            xx, yy, _, _ = st
            rr = max(np.hypot(xx, yy), 1e-9)
            g = 0.5 * dn2_at(rr)
            return np.array([st[2], st[3], g * xx / rr, g * yy / rr])
        st = np.array([x, y, px, py])
        k1 = F(st)
        k2 = F(st + 0.5 * ds * k1)
        k3 = F(st + 0.5 * ds * k2)
        k4 = F(st + ds * k3)
        st = st + ds / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        x, y, px, py = st
        if x > 60.0:
            break
    # KU wirowi (wir w (0,0), b>0) = py < 0 -> alpha_ku = -atan2(py,px)
    return -np.degrees(np.arctan2(py, px))


print("    b:        " + "".join(f"  {bb:7.0f}" for bb in (4, 6, 8, 12)))
for om in (1.0, 1.1, 1.3):
    row = f"    omega={om:.1f}:"
    for bb in (4.0, 6.0, 8.0, 12.0):
        row += f"  {trace(om, bb):+7.3f}"
    print(row)
print()

print("--- 4. progi (LOCK, sekcja 3) ---")
print(f"  G4a: lambda_min(L_full) >= -1e-3 (rozdzielczosc iteracji przy "
      f"kontinuum Goldstone'a od 0)")
print(f"  G4b: gamma_fit <= 0.01  (amplifikacja <= e^{0.01*150:.1f} = "
      f"{np.exp(1.5):.1f} na oknie tau=150; kontrast B2: gamma=0.125 -> e^19)")
print()
print("=== WNIOSKI PHASE 0 ===")
print(f"  R_c_pred = {Rc:.3f}; n(r) > 1 wszedzie (deficyt = osrodek gestszy,")
print(f"  refrakcja KU wirowi, UNIWERSALNA — n nie zalezy od znaku kretu);")
print(f"  ogon n^2-1 ~ {A_t:.3f}/(r^2*(omega^2-V''(s*))) — algebraiczny;")
print(f"  alpha_pred zalockowane jako zaklad dla opa pomiarowego")
