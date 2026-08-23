# Phase 0 check: weryfikacja rownan #1 i #3 PRZED eksperymentem wlasciwym.
# LOCK: Phase0_balance.md (sekcje 1-2, doprecyzowania D1-D10).
#
# #1: masa ekranowania m0 = sqrt(a/kappa) = 1.0; ogon radialny ~ K0(m0 r),
#     czyli ln(u*sqrt(r)) liniowe w r ze zboczem -m0.
# #3: front planarny 1D: v_planar = DeltaV / int(f'^2 dx); przewidywanie
#     c ~ kappa = 0.5 w prawie v(R) = v_inf - c/R. Oszacowanie v_inf.

import numpy as np

kappa = 0.50
a, b, c = 0.50, 1.60, 1.00
dt = 0.02

print("=== Phase 0 check: rachunki pre-eksperyment (LOCK: Phase0_balance.md) ===")
print()

# ---------- stale potencjalu ----------
disc = np.sqrt(b * b - 4 * a * c)
s_star = (b + disc) / (2 * c)          # prawdziwa proznia
s_bar = (b - disc) / (2 * c)           # bariera
def V(s): return 0.5 * a * s**2 - (b / 3.0) * np.abs(s)**3 + 0.25 * c * s**4
def Vp(s): return s * (a - b * np.abs(s) + c * s * s)
DeltaV = V(0.0) - V(s_star)
Vpp0 = a
Vpps = a - 2 * b * s_star + 3 * c * s_star**2
m0 = np.sqrt(Vpp0 / kappa)
print(f"s* = {s_star:.5f}, s_bar = {s_bar:.5f}")
print(f"V(s*) = {V(s_star):+.6f}  -> DeltaV = V(0)-V(s*) = {DeltaV:.6f}")
print(f"V''(0) = {Vpp0:.4f}, V''(s*) = {Vpps:.5f}")
print(f"ROWNANIE #1: m0 = sqrt(V''(0)/kappa) = {m0:.4f}  (predykcja LOCK: 1.0)")
ok_m0 = abs(m0 - 1.0) < 1e-12
print()

# ---------- #1: ogon radialny K0 (linearyzacja wokol s=0) ----------
# kappa*(u'' + u'/r) = a*u ; relaksacja BVP na siatce radialnej,
# u(r0)=1, u(Rmax)=0; test: slope d/dr ln(u*sqrt(r)) -> -m0.
r0, Rmax, drr = 2.0, 22.0, 0.01
r = np.arange(r0, Rmax + drr / 2, drr)
u = np.exp(-(r - r0))  # start
u[0], u[-1] = 1.0, 0.0
for _ in range(400000):
    upp = (u[2:] - 2 * u[1:-1] + u[:-2]) / drr**2
    upr = (u[2:] - u[:-2]) / (2 * drr) / r[1:-1]
    u_new = u.copy()
    u_new[1:-1] = u[1:-1] + 0.2 * drr**2 / kappa * (kappa * (upp + upr) - a * u[1:-1])
    if np.max(np.abs(u_new - u)) < 1e-14:
        u = u_new; break
    u = u_new
# fit zbocza w oknie [6, 16] (z dala od obu brzegow)
msk = (r >= 6.0) & (r <= 16.0)
y = np.log(u[msk] * np.sqrt(r[msk]))
A_ = np.vstack([r[msk], np.ones(msk.sum())]).T
slope, icpt = np.linalg.lstsq(A_, y, rcond=None)[0]
ok_tail = abs(-slope - m0) <= 0.05
print(f"OGON RADIALNY: slope ln(u*sqrt(r)) = {slope:+.5f}  (K0 => -m0 = {-m0:.4f})")
print(f"  |slope+m0| = {abs(slope + m0):.5f} <= 0.05 : {'PASS' if ok_tail else 'FAIL'}")
print()

# ---------- #3: front planarny 1D ----------
Nx, dx = 1200, 0.5
x = np.arange(Nx) * dx
s = np.where(x < 150.0, s_star, 0.0) * 0.5 * (1 - np.tanh((x - 150.0) / 2.0)) \
    + s_star * 0.5 * (1 - np.tanh((x - 150.0) / 2.0)) * 0.0
s = s_star * 0.5 * (1.0 - np.tanh((x - 150.0) / 2.0))
def lap1(f):
    g = np.empty_like(f)
    g[1:-1] = (f[2:] - 2 * f[1:-1] + f[:-2]) / dx**2
    g[0] = g[1]; g[-1] = g[-2]
    return g
def front_pos(f):
    eps_s = np.sqrt(0.30)
    idx = np.where(f > eps_s)[0]
    if len(idx) == 0: return np.nan
    i = idx[-1]
    if i + 1 >= len(f): return x[i]
    return x[i] + dx * (f[i] - eps_s) / (f[i] - f[i + 1])
pos, taus = [], []
steps = 30000
for step in range(steps + 1):
    if step % 100 == 0:
        pos.append(front_pos(s)); taus.append(step * dt)
    s = s + dt * (kappa * lap1(s) - Vp(s))
pos, taus = np.array(pos), np.array(taus)
half = len(pos) // 2                     # fit na drugiej polowie (ustalony ksztalt)
A_ = np.vstack([taus[half:], np.ones(len(taus) - half)]).T
v_planar, _ = np.linalg.lstsq(A_, pos[half:], rcond=None)[0]
sp = (s[1:] - s[:-1]) / dx
sigma_int = float(np.sum(sp * sp) * dx)  # int f'^2 dx
v_theory = DeltaV / sigma_int
rel = abs(v_planar - v_theory) / v_theory
ok_v = rel <= 0.05
print(f"ROWNANIE #3 (front planarny 1D, {steps} krokow):")
print(f"  v_planar(zmierzone)  = {v_planar:.6f}")
print(f"  int f'^2 dx          = {sigma_int:.6f}")
print(f"  v_teoria = DeltaV/int = {v_theory:.6f}   (rozn. wzgledna {rel*100:.2f}% <= 5%: "
      f"{'PASS' if ok_v else 'FAIL'})")
print(f"  => oszacowanie v_inf ~ {v_planar:.4f}; przewidywanie c ~ kappa = {kappa}")
Rc = kappa / v_planar
print(f"  => promien krytyczny R_c = kappa/v_inf ~ {Rc:.2f}  "
      f"(obiekty R_obj in {{8,12,16}} > R_c: {'OK, rosna' if Rc < 8 else 'UWAGA: kolaps!'})")
print(f"  => szacowany czas domkniecia szczeliny g0=8: tau ~ {8.0/(2*v_planar):.0f}")
print()
allok = ok_m0 and ok_tail and ok_v
print(f"PODSUMOWANIE PHASE 0 CHECK: m0={'PASS' if ok_m0 else 'FAIL'} "
      f"ogon_K0={'PASS' if ok_tail else 'FAIL'} front_1D={'PASS' if ok_v else 'FAIL'} "
      f"-> {'PASS' if allok else 'FAIL'}")
