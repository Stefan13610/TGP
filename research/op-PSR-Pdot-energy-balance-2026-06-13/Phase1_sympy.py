# -*- coding: utf-8 -*-
# Phase 1 — op-PSR-Pdot-energy-balance-2026-06-13
# Energy balance of TGP radiative sector vs J0737-3039 orbital decay.
# All derivations ab initio in sympy; observational anchors as locked numerics.
# NO hardcoded T_pass. Run: python Phase1_sympy.py

import sympy as sp

PASS, FAIL = "PASS", "FAIL"
results = []

def check(name, cond, note=""):
    verdict = PASS if cond else FAIL
    results.append((name, verdict, note))
    print(f"[{verdict}] {name}" + (f"  -- {note}" if note else ""))

print("=" * 72)
print("T0 — dipole cancellation for C_i = q*m_i (circular two-body, CM frame)")
print("=" * 72)
t, q, m1, m2, d, w, G, c = sp.symbols('t q m1 m2 d omega G c', positive=True)
Mtot = m1 + m2
mu = m1 * m2 / Mtot
x1 = (m2 / Mtot) * d * sp.Matrix([sp.cos(w * t), sp.sin(w * t), 0])
x2 = -(m1 / Mtot) * d * sp.Matrix([sp.cos(w * t), sp.sin(w * t), 0])
D_scalar = q * m1 * x1 + q * m2 * x2
check("T0 scalar dipole D = q*Sum(m_i x_i) == 0 identically",
      sp.simplify(D_scalar) == sp.zeros(3, 1))

print("=" * 72)
print("T2a — static calibration: q^2 = 4*pi*G from Newtonian limit")
print("=" * 72)
# Field eq: Box psi = -q rho  (Box = d_t^2 - Lap). Static: Lap psi = q rho.
# Point mass m1: psi(r) = -q*m1/(4 pi r).  Verify Lap(psi)=0 for r>0:
r = sp.symbols('r', positive=True)
psi_pt = -q * m1 / (4 * sp.pi * r)
lap = sp.diff(r**2 * sp.diff(psi_pt, r), r) / r**2
check("T2a.1 Lap(-q m1/(4 pi r)) == 0 for r>0", sp.simplify(lap) == 0)
# Interaction energy of m2 in psi_1: L_int = -q rho psi -> U = +q m2 psi_1(d)
U_int = q * m2 * psi_pt.subs(r, d)
q2_solved = sp.solve(sp.Eq(U_int, -G * m1 * m2 / d), q**2)
check("T2a.2 U_int = -G m1 m2/d  <=>  q^2 = 4*pi*G",
      len(q2_solved) == 1 and sp.simplify(q2_solved[0] - 4 * sp.pi * G) == 0,
      f"q^2 = {q2_solved}")

print("=" * 72)
print("T1 — angular integrals over the unit sphere")
print("=" * 72)
th, ph = sp.symbols('theta phi', real=True)
n = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])

def ang_int(expr):
    return sp.integrate(sp.integrate(sp.expand_trig(expr) * sp.sin(th),
                                     (th, 0, sp.pi)), (ph, 0, 2 * sp.pi))

check("T1.1  Int n_x^4 dOmega = 4*pi/5",
      sp.simplify(ang_int(n[0]**4) - sp.Rational(4, 5) * sp.pi) == 0)
check("T1.2  Int n_x^2 n_y^2 dOmega = 4*pi/15",
      sp.simplify(ang_int(n[0]**2 * n[1]**2) - sp.Rational(4, 15) * sp.pi) == 0)
check("T1.3  Int n_x^3 n_y dOmega = 0",
      sp.simplify(ang_int(n[0]**3 * n[1])) == 0)

print("=" * 72)
print("T3 — scalar quadrupole power for circular binary, ab initio")
print("=" * 72)
# Far field (retarded, c=1 internal units; ratio to P_GR is c-independent):
#   psi = (q/(4 pi r)) [ M + n.D_dot + (1/2) n_i n_j Mddot_ij ](t_ret)
# monopole const, dipole zero (T0)  =>  psi_dot = (q/(8 pi r)) n_i n_j M3_ij
M_ij = m1 * x1 * x1.T + m2 * x2 * x2.T          # second mass moment
trM = sp.simplify(sp.trace(M_ij))
check("T3.1 trace M_ij = mu*d^2 = const (no breathing for circular orbit)",
      sp.simplify(sp.diff(trM, t)) == 0, f"trM = {trM}")
M3 = sp.simplify(sp.diff(M_ij, t, 3))            # third time derivative
nn_M3 = sp.expand_trig(sp.simplify((n.T * M3 * n)[0, 0]))
psidot_sq_ang = ang_int(nn_M3**2)                # angular integral of (n.M3.n)^2
# time average over one orbital period (integrand periodic with 2w):
T_orb = 2 * sp.pi / w
P_phi = sp.simplify((q**2 / (64 * sp.pi**2)) *
                    sp.integrate(psidot_sq_ang, (t, 0, T_orb)) / T_orb)
P_phi = sp.simplify(P_phi.subs(q**2, 4 * sp.pi * G))
P_phi_expected = sp.Rational(16, 15) * G * mu**2 * d**4 * w**6
check("T3.2 P_phi = (16/15) G mu^2 d^4 omega^6   [ab initio]",
      sp.simplify(P_phi - P_phi_expected) == 0, f"P_phi = {P_phi}")

# Cross-form: P_phi == (G/30) <I3:I3> with I = traceless part of M
I_ij = M_ij - sp.eye(3) * trM / 3
I3 = sp.diff(I_ij, t, 3)
I3_sq = sp.simplify(sum(I3[i, j]**2 for i in range(3) for j in range(3)))
I3_avg = sp.simplify(sp.integrate(I3_sq, (t, 0, T_orb)) / T_orb)
check("T3.3 P_phi == (G/30)*<I3:I3>  (user's formula reproduced)",
      sp.simplify(P_phi - G * I3_avg / 30) == 0, f"<I3:I3> = {I3_avg}")

# GR anchor (LOCKED textbook): P_GR = (G/5) <I3:I3>
P_GR = sp.simplify(G * I3_avg / 5)
check("T3.4 P_GR = (32/5) G mu^2 d^4 omega^6  (lock consistency)",
      sp.simplify(P_GR - sp.Rational(32, 5) * G * mu**2 * d**4 * w**6) == 0)
ratio_scalar = sp.simplify(P_phi / P_GR)
check("T3.5 P_phi / P_GR = 1/6 EXACT",
      sp.simplify(ratio_scalar - sp.Rational(1, 6)) == 0,
      f"ratio = {ratio_scalar}")

print("=" * 72)
print("T4 — sigma-channel propagation gate at J0737-3039 frequencies")
print("=" * 72)
hbar_eVs = 6.582119569e-16          # eV*s
c_SI = 2.99792458e8                 # m/s
m_sigma_eV = 0.71e-3                # eV  (LOCKED: sigma-yukawa-audit)
Pb_s = 0.10225156248 * 86400.0      # s   (LOCKED: Kramer 2021)
w_orb = 2 * 3.141592653589793 / Pb_s
w_gw = 2 * w_orb                    # quadrupole radiates at 2*omega_orb
E_gw_eV = hbar_eVs * w_gw
check("T4.1 hbar*omega_GW < m_sigma c^2  =>  sigma channel NON-PROPAGATING",
      E_gw_eV < m_sigma_eV,
      f"hbar*w_GW = {E_gw_eV:.3e} eV vs m_sigma = {m_sigma_eV:.2e} eV "
      f"(ratio {E_gw_eV/m_sigma_eV:.1e})")
yukawa_range_m = hbar_eVs * c_SI / m_sigma_eV
check("T4.2 static range of sigma = hbar*c/m_sigma ~ sub-mm (cannot carry Newton)",
      yukawa_range_m < 1e-2, f"range = {yukawa_range_m*1e3:.3f} mm")
# dispersion: k^2 = w^2 - m^2 < 0 -> evanescent (symbolic)
ww, mm = sp.symbols('w_ m_', positive=True)
k_sq = ww**2 - mm**2
check("T4.3 k^2 < 0 for omega < m (evanescent, zero flux at infinity)",
      sp.simplify(k_sq.subs([(ww, 1), (mm, 2)])) < 0)

print("=" * 72)
print("T5 — Branch B: amplitude lock does NOT pin the energy flux")
print("=" * 72)
# sigma EOM: Box sigma = -xi T_TT  ->  far field sigma = (xi/(4 pi r)) S
# GR:        Box h_TT = -16 pi G T_TT -> h_GR = (16 pi G/(4 pi r)) S
# detector metric: h_eff = lam * sigma
lam, xi, S = sp.symbols('lambda_ xi_ S', positive=True)
h_eff_over_h_GR = sp.simplify((lam * xi) / (16 * sp.pi * G))
amp_condition = sp.Eq(lam * xi, 16 * sp.pi * G)
# Energy flux: canonical sigma kinetic -(1/4)(d sigma)^2  -> flux = (1/2)<sdot^2>
# GR Isaacson: flux = (1/(32 pi G)) <hdot^2>;  with h = lam*sigma:
kappa_E = sp.simplify(((sp.Rational(1, 2)) / (lam**(-2) * lam**2)) /
                      (lam**2 / (32 * sp.pi * G)))   # = 16 pi G / lam^2
kappa_E = sp.simplify(16 * sp.pi * G / lam**2)
# under amplitude condition: kappa_E = xi/lam  (independent of amplitude match)
kappa_E_amp = sp.simplify(kappa_E.subs(lam, 16 * sp.pi * G / xi))
check("T5.1 amplitude match h_eff=h_GR  <=>  lam*xi = 16 pi G (1 condition)",
      sp.simplify(h_eff_over_h_GR.subs(lam, 16 * sp.pi * G / xi) - 1) == 0)
check("T5.2 flux ratio kappa_E = xi^2/(16 pi G) under amplitude lock "
      "-> NOT fixed to 1 (2nd independent condition xi=lam required)",
      sp.simplify(kappa_E_amp - xi**2 / (16 * sp.pi * G)) == 0,
      f"kappa_E = {kappa_E_amp}")
# canonical symmetric choice xi = lam = sqrt(16 pi G) -> kappa_E = 1:
check("T5.3 symmetric choice xi=lam gives kappa_E = 1 (GR flux) exactly",
      sp.simplify(kappa_E_amp.subs(xi, sp.sqrt(16 * sp.pi * G)) - 1) == 0)

print("=" * 72)
print("T6 — verdicts vs J0737-3039 (R_obs = 0.999963 +/- 0.000063)")
print("=" * 72)
R_obs, sig_obs = 0.999963, 0.000063
R_A = 1.0 / 6.0                      # Branch A: sigma dead, scalar only
R_B_canonical = 1.0 / 6.0 + 1.0      # Branch B with kappa_E = 1
sigmas_A = abs(R_A - R_obs) / sig_obs
sigmas_B = abs(R_B_canonical - R_obs) / sig_obs
kappa_E_required = R_obs - 1.0 / 6.0
kappa_E_window_5s = 5 * sig_obs
check("T6.1 Branch A: R = 1/6 deviates by > 5 sigma",
      sigmas_A > 5, f"R_A = {R_A:.6f}, deviation = {sigmas_A:.0f} sigma")
check("T6.2 Branch B (canonical kappa_E=1): R = 7/6 deviates by > 5 sigma",
      sigmas_B > 5, f"R_B = {R_B_canonical:.6f}, deviation = {sigmas_B:.0f} sigma")
print(f"      Branch B rescue would need kappa_E = {kappa_E_required:.6f}"
      f" +/- {kappa_E_window_5s:.6f} (rel. {kappa_E_window_5s/kappa_E_required:.1e})"
      f" -> NEW tuned anchor, forbidden by Phase 0 par.6.2")
# cross-check Hulse-Taylor (qualitative; e=0.617 so ratio approximate):
R_obs_HT, sig_HT = 0.9983, 0.0016
print(f"      B1913+16 cross-check: Branch A {abs(R_A-R_obs_HT)/sig_HT:.0f} sigma,"
      f" Branch B {abs(R_B_canonical-R_obs_HT)/sig_HT:.0f} sigma (same verdicts)")

print("=" * 72)
print("T7 — eccentricity robustness (e = 0.0877775)")
print("=" * 72)
e_val = 0.0877775
# Peters tensor enhancement f(e); scalar-channel enhancement differs at O(e^2).
f_peters = (1 + 73.0/24*e_val**2 + 37.0/96*e_val**4) / (1 - e_val**2)**3.5
# conservative bound on ratio shift: full size of either enhancement
ratio_shift_bound = abs(f_peters - 1.0)
gap_A = abs(R_obs - 1.0/6.0)
gap_B = abs(R_B_canonical - R_obs)
check("T7.1 O(e^2) channel-ratio shift cannot bridge Branch A gap",
      (1.0/6.0) * ratio_shift_bound * 3 < gap_A,
      f"bound on ratio shift ~ {ratio_shift_bound*100:.1f}% vs gap {gap_A:.3f}")
check("T7.2 O(e^2) channel-ratio shift cannot bridge Branch B gap",
      1.0 * ratio_shift_bound * 3 < gap_B,
      f"bound ~ {ratio_shift_bound*100:.1f}% vs gap {gap_B:.3f}")

print("=" * 72)
n_pass = sum(1 for _, v, _ in results if v == PASS)
print(f"TOTAL: {n_pass}/{len(results)} PASS")
for name, v, _ in results:
    if v == FAIL:
        print(f"  FAILED: {name}")
