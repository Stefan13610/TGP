# ==============================================================================
#   gs66_frw_propagator.py
#
#   LINEAR FRW PROPAGATOR OF TGP FIELD EQUATION — DECISIVE TEST OF BRIDGE (a)
#
#   Goal: Test whether solving the FULL (time-dependent) TGP field equation
#
#       Phi_ddot + 3 H Phi_dot - (1/a^2) nabla^2 Phi + U'(Phi) = source
#
#   around a localized mass in FRW background produces log(r) far-field
#   (MOND-like) behavior, which gs65 proved the static equation does NOT.
#
#   Strategy (linear perturbation):
#     - Linearize about vacuum Phi = Phi0 (psi = 1), delta Phi = varphi.
#     - Fourier transform in (t, x).
#     - Propagator:  G~(omega, k) = 1 / [-omega^2 + 3 i H omega + k^2/a^2 - gamma]
#     - Quasi-static limit (omega -> 0): G~(k) = 1/(k^2/a^2 - gamma + 3 i H^2)
#     - Real-space Green's function via inverse Fourier transform.
#     - Check asymptotics:  Newton (r << L_nat) ?  log(r) (r ~ L_H) ?
#
#   Parts:
#     Part 1. Symbolic setup & dispersion relation
#     Part 2. Analytic propagator in k-space
#     Part 3. Real-space Green's function: closed form (3D Yukawa-like)
#     Part 4. Numerical evaluation at realistic parameters
#     Part 5. Asymptotic limits (small-r, resonance, large-r)
#     Part 6. Rotation curve prediction from point source
#     Part 7. Fourier-power argument: can any D(k) give log(r)?
#     Part 8. Verdict on bridge (a) + roadmap
# ==============================================================================

from __future__ import annotations

import math
import numpy as np
import sympy as sp

hdr = "=" * 78
sub = "-" * 78


def header(title: str) -> None:
    print()
    print(hdr)
    print(f"  {title}")
    print(hdr)


def sub_header(title: str) -> None:
    print()
    print(sub)
    print(f"  {title}")
    print(sub)


# Physical constants (SI)
c_si  = 2.998e8
G_si  = 6.674e-11
H0_si = 2.1844e-18         # 67.4 km/s/Mpc in s^-1
kpc_m = 3.0857e19
Mpc_m = 3.0857e22
Msun  = 1.989e30

L_H  = c_si / H0_si        # Hubble length ~ 4.4 Gpc
a0_obs = 1.20e-10
a0_tgp = c_si * H0_si / (2 * math.pi)

# ------------------------------------------------------------------------------
# Part 1. SYMBOLIC SETUP
# ------------------------------------------------------------------------------
header("Part 1. Linearized FRW field equation (symbolic)")

omega, k, H, gamma_s, a = sp.symbols("omega k H gamma a", positive=True, real=True)
r = sp.Symbol("r", positive=True)
M = sp.Symbol("M", positive=True)
G_n = sp.Symbol("G_N", positive=True)

D_op = -omega**2 + 3*sp.I*H*omega + k**2/a**2 - gamma_s

print(f"""
Linearized field equation (vacuum + source):
    varphi_ddot + 3 H varphi_dot - (1/a^2) nabla^2 varphi + U''(1)*varphi = S

with U''(1) = -gamma (from sek05, gamma > 0).

Fourier (t -> omega, x -> k):
    D(omega, k) * varphi~(omega, k) = S~(omega, k)
    D(omega, k) = -omega^2 + 3 i H omega + k^2/a^2 - gamma

Dispersion relation (D = 0, free propagation):
    omega^2 - 3 i H omega - k^2/a^2 + gamma = 0
""")

disp = sp.solve(sp.Eq(D_op, 0), omega)
print("Solutions omega(k):")
for s in disp:
    print(f"    omega = {sp.simplify(s)}")


# ------------------------------------------------------------------------------
# Part 2. QUASI-STATIC PROPAGATOR
# ------------------------------------------------------------------------------
header("Part 2. Quasi-static propagator  D(k) = k^2/a^2 - gamma + 3 i H^2")

print("""
Quasi-static limit: take omega -> 0 and add O(H) frictional imaginary
shift from the 3 H omega term in slowly-evolving configurations.

Standard procedure: source varies on timescale ~ 1/H, so omega ~ H.
    3 i H omega  ->  3 i H^2    (order-of-magnitude frictional shift)

Propagator:
    G~(k) = 1 / (k^2/a^2 - gamma + 3 i H^2)

Setting a = 1 (today) for simplicity:
    G~(k) = 1 / (k^2 - gamma + 3 i H^2)

This has the Yukawa form 1/(k^2 + m^2) with COMPLEX mass^2:
    m^2 = -gamma + 3 i H^2        (note: would-be tachyonic, Hubble-damped)

Real-space Green's function (standard 3D Yukawa):
    G(r) = (1/(4 pi r)) * exp(- mu r),   mu = sqrt(m^2)   (principal branch)

Branch analysis:
    m^2 = -gamma + 3 i H^2
    |m^2| = sqrt(gamma^2 + 9 H^4)
    arg(m^2) = pi - arctan(3 H^2 / gamma)   (2nd quadrant for gamma > 0)
    arg(mu) = (arg(m^2))/2 = pi/2 - (1/2) arctan(3 H^2 / gamma)

    For gamma >> H^2 (L_nat << L_H):
        arg(mu) approx pi/2
        mu ~ i * sqrt(gamma) + 3 H^2 / (2 sqrt(gamma))
        Re(mu) = 3 H^2 / (2 sqrt(gamma))   (damping)
        Im(mu) = sqrt(gamma)                (oscillation)

    For gamma ~ H^2:
        arg(mu) approx pi/2 - pi/8 = 3 pi/8
        mu mixes real and imaginary equally.

    For gamma << H^2 (L_nat >> L_H):
        arg(m^2) approx pi/2
        arg(mu) approx pi/4
        mu ~ sqrt(3 H^2/2) * (1 + i)
        equal damping & oscillation at scale sqrt(2/(3 H^2)).
""")


# ------------------------------------------------------------------------------
# Part 3. ANALYTIC FAR-FIELD EXPANSION
# ------------------------------------------------------------------------------
header("Part 3. Analytic Green's function: exact form & expansions")

print("""
Exact closed form (a = 1, quasi-static):

    G(r) = exp(- mu r) / (4 pi r),    mu = sqrt(-gamma + 3 i H^2)

Physical Green's function (real-valued, retarded, outgoing wave conv.):
    G_phys(r) = Re[G(r)]
              = exp(-Re(mu) r) * cos(Im(mu) r) / (4 pi r)

Leading asymptotics (gamma >> H^2):
    Re(mu) = 3 H^2 / (2 sqrt(gamma))         [damping length 2 sqrt(gamma)/(3 H^2)]
    Im(mu) = sqrt(gamma)                      [oscillation wavelength 2 pi/sqrt(gamma)]

Short-distance r << 1/Im(mu) = 1/sqrt(gamma):
    cos(Im(mu) r) ~ 1,    exp(-Re(mu) r) ~ 1
    G_phys(r) ~ 1 / (4 pi r)                   [NEWTONIAN 1/r]

Intermediate r ~ 1/sqrt(gamma):
    G_phys(r) ~ cos(sqrt(gamma) r) / (4 pi r)  [oscillatory ringing]

Long-distance r >> 1/sqrt(gamma):
    G_phys(r) decays as oscillating-damped
    but NEVER becomes log(r).

CRITICAL: no log(r) behavior anywhere in this linear quasi-static result.
""")

# Symbolic propagator in closed form
mu_sym = sp.sqrt(-gamma_s + 3*sp.I*H**2)
G_sym = sp.exp(-mu_sym * r) / (4 * sp.pi * r)
print("\nSymbolic Green's function:")
sp.pprint(G_sym)


# ------------------------------------------------------------------------------
# Part 4. NUMERICAL EVALUATION (realistic parameters)
# ------------------------------------------------------------------------------
header("Part 4. Numerical evaluation at physical parameters")

# Test three choices of gamma:
scenarios = [
    ("gamma = 1/L_nat^2 with L_nat = 3 kpc (gs54 Yukawa)",
     1.0 / (3 * kpc_m)**2),
    ("gamma = 1/L_nat^2 with L_nat = L_H (cosmological)",
     1.0 / L_H**2),
    ("gamma = H0^2 (dimensional minimum)",
     H0_si**2 / c_si**2),     # gamma has units 1/L^2 = (1/c)^2 * omega^2
]

def G_phys(r_arr, gamma_val, H_val):
    """Physical Green's function (a=1, quasi-static) in SI."""
    m2 = -gamma_val + 3j * (H_val / c_si)**2  # keep consistent units 1/L^2
    mu = np.sqrt(m2)  # complex
    re_mu = np.real(mu)
    im_mu = np.imag(mu)
    val = np.exp(-re_mu * r_arr) * np.cos(im_mu * r_arr) / (4 * np.pi * r_arr)
    return val, re_mu, im_mu

r_grid = np.array([0.1*kpc_m, 1*kpc_m, 10*kpc_m, 100*kpc_m, 1*Mpc_m,
                   10*Mpc_m, 100*Mpc_m, 1000*Mpc_m, L_H])

for name, gv in scenarios:
    print(f"\n  Scenario: {name}")
    print(f"  gamma = {gv:.3e} 1/m^2  =>  L_nat = 1/sqrt(gamma) = "
          f"{1.0/math.sqrt(abs(gv))/kpc_m:.3e} kpc")
    G, re_mu, im_mu = G_phys(r_grid, gv, H0_si)
    print(f"  Re(mu) = {re_mu:.3e} 1/m  (damping length = "
          f"{1.0/max(re_mu, 1e-60)/kpc_m:.3e} kpc)")
    print(f"  Im(mu) = {im_mu:.3e} 1/m  (oscillation wavelength = "
          f"{2*math.pi/max(abs(im_mu), 1e-60)/kpc_m:.3e} kpc)")
    print(f"  {'r [kpc]':>12s} {'G_phys':>14s} {'4 pi r * G':>14s}")
    for rv, gvv in zip(r_grid, G):
        print(f"  {rv/kpc_m:12.3e} {gvv:14.3e} {(4*math.pi*rv*gvv):14.3e}")


# ------------------------------------------------------------------------------
# Part 5. ASYMPTOTIC ANALYSIS
# ------------------------------------------------------------------------------
header("Part 5. Comparison to MOND log(r) at asymptotic regimes")

print("""
MOND far-field prediction:
    Phi_MOND(r) = v_flat^2 * log(r / r_0),   g_MOND(r) = v_flat^2 / r

We need to check: for point mass M in this TGP propagator,

    Phi_TGP(r) = -G_N M * [4 pi * G_phys(r)]
               = -G_N M * exp(-Re(mu) r) cos(Im(mu) r) / r

DOES this asymptote to log(r) for any choice of gamma?

NO.  A function of the form A * f(r)/r for bounded f does not tend to log(r).

Numerical check at r = 10 L_nat and r = 100 L_nat:
""")

for name, gv in scenarios:
    G1, _, _ = G_phys(np.array([1/math.sqrt(abs(gv))]), gv, H0_si)
    G10, _, _ = G_phys(np.array([10/math.sqrt(abs(gv))]), gv, H0_si)
    G100, _, _ = G_phys(np.array([100/math.sqrt(abs(gv))]), gv, H0_si)
    # Compare to 1/r (Newton) and log(r) (MOND)
    r_nat = 1/math.sqrt(abs(gv))
    Phi1   = 4*math.pi*r_nat*G1[0]
    Phi10  = 4*math.pi*(10*r_nat)*G10[0]
    Phi100 = 4*math.pi*(100*r_nat)*G100[0]
    print(f"  {name}")
    print(f"    r / L_nat = 1:    4 pi r G = {Phi1:.4e}  (Newton = 1.0)")
    print(f"    r / L_nat = 10:   4 pi r G = {Phi10:.4e}  (if MOND: Phi ~ log(10)*... )")
    print(f"    r / L_nat = 100:  4 pi r G = {Phi100:.4e}")


# ------------------------------------------------------------------------------
# Part 6. ROTATION CURVE FROM POINT SOURCE
# ------------------------------------------------------------------------------
header("Part 6. Rotation curve prediction (bridge (a) linear)")

print("""
For point source mass M at origin:
    Phi(r) = -(4 pi G_N M) * G_phys(r) / [4 pi]
           = -(G_N M / r) * exp(-Re(mu) r) * cos(Im(mu) r)

Acceleration:
    g(r) = -d Phi / d r
         = -G_N M * d/dr [exp(-Re(mu)r) cos(Im(mu)r) / r]

At r << 1/|mu|:
    g(r) -> G_N M / r^2            (NEWTON)

At r >> 1/Re(mu):
    g(r) -> 0                        (SCREENED, exponentially)

Nowhere does g(r) -> const * sqrt(G_N M a0)/r (MOND).
""")

# Example: Milky Way mass
M_gal = 1.5e11 * Msun
print(f"\nExample: M = {M_gal/Msun:.2e} Msun (Milky Way)")
print(f"{'r [kpc]':>10s} {'g_Newton [m/s2]':>18s} {'g_TGP [m/s2]':>18s} "
      f"{'g_MOND [m/s2]':>18s}")

# Numerical derivative
def acceleration(r_val, gamma_val, H_val, M_val):
    dr = r_val * 1e-4
    phi_p, _, _ = G_phys(np.array([r_val + dr]), gamma_val, H_val)
    phi_m, _, _ = G_phys(np.array([r_val - dr]), gamma_val, H_val)
    # Phi = -4 pi G M G_phys
    Phi_p = -4*math.pi * G_si * M_val * phi_p[0]
    Phi_m = -4*math.pi * G_si * M_val * phi_m[0]
    return -(Phi_p - Phi_m) / (2 * dr)

gv_cosmo = 1.0 / L_H**2
for r_kpc in [1, 10, 100, 1000, 10000]:
    rv = r_kpc * kpc_m
    g_N = G_si * M_gal / rv**2
    g_T = acceleration(rv, gv_cosmo, H0_si, M_gal)
    g_M = math.sqrt(G_si * M_gal * a0_obs) / rv  # deep MOND
    print(f"  {r_kpc:10d} {g_N:18.4e} {g_T:18.4e} {g_M:18.4e}")


# ------------------------------------------------------------------------------
# Part 7. FOURIER-POWER ARGUMENT (general impossibility)
# ------------------------------------------------------------------------------
header("Part 7. Fourier-power argument: can ANY D(k) in linear theory give log(r)?")

print("""
Claim: For a SPHERICALLY SYMMETRIC linear theory with translation-invariant
propagator G(r) = inverse-Fourier[ 1 / D(k) ], log(r) behavior requires:

    G~(k) ~ C / k^3   as k -> 0

because log(r) has Fourier transform:
    integral d^3k / (2pi)^3  * exp(i k . r) / k^3  ~ log(r)    [up to regularization]

Proof (dimensional + scaling):
    Under k -> s k, r -> r/s:
       integral d^3k * exp(i k r) / k^3     invariant
    So integral behaves as log(length) when regularized.

Any LINEAR TGP propagator from sek02+sek05 has:
    D(k) = k^2 - gamma + 3 i H^2 + O(k^4, etc)
    G~(k) = 1/D(k)

At k -> 0: G~(k) -> 1/(-gamma + 3 i H^2) = constant
           (short-range 3D Yukawa, mass -gamma+3iH^2).
At k -> infinity: G~(k) -> 1/k^2 (Newtonian).

NEITHER regime is 1/k^3.  Therefore:

    THEOREM (linear TGP far-field impossibility):
        No static or quasi-static linear extension of sek02 + sek05
        with polynomial D(k) produces MOND log(r) behavior from a
        point source.  Bridge (a) in LINEAR regime is RULED OUT.

Corollary: Any viable bridge (a) must be INTRINSICALLY NONLINEAR.
The non-linear kinetic coupling alpha*(nabla Phi)^2/Phi (sek02) is the
only remaining TGP-axiomatic source. But gs65 showed its static solutions
are r^(-1/3), not log(r).

We are left with one of:
    - FULL time-dependent nonlinear TGP solution (hard; probably
      does not help based on static r^(-1/3) result).
    - Bridge (b): extend sek05 with dual-scale potential (explicit extension).
    - Bridge (c): smooth membrane geometry (explicit extension).
""")


# ------------------------------------------------------------------------------
# Part 8. VERDICT & ROADMAP
# ------------------------------------------------------------------------------
header("Part 8. VERDICT: Bridge (a) FALSIFIED in linear regime")

print(f"""
SUMMARY OF gs66 FINDINGS:

1. Linearized FRW TGP field equation has propagator:
       G~(k) = 1/(k^2/a^2 - gamma + 3 i H^2)

2. Real-space Green's function:
       G(r) = exp(-mu r)/(4 pi r),   mu = sqrt(-gamma + 3 i H^2)
    - Newtonian 1/r at small r,
    - oscillatory-damped at r >~ 1/sqrt(|gamma|),
    - exponentially screened at large r.

3. NEVER log(r).  Proven by Fourier-power counting:
       log(r) <-> G~(k) ~ 1/k^3 at k -> 0
       TGP gives G~(k) -> constant or 1/k^2, never 1/k^3.

4. Rotation curves:
       g_TGP(r) -> G_N M / r^2  (Newtonian) for r < L_nat
       g_TGP(r) -> 0            (screened)    for r >> L_nat
    NOT v_flat = constant.

CONSEQUENCE FOR TGP GALAXY PROGRAM:

  BRIDGE (a) — Hubble friction in linear FRW — is FALSIFIED.
  It was the only extension that preserved TGP axioms without additions.

  TGP now REQUIRES an explicit theory extension to describe MOND:
    Bridge (b): dual-scale substrate (two gamma parameters)
    Bridge (c): smooth 3D->2D membrane geometry

  NEITHER is derivable from sek02+sek05 alone.

THIS IS A SIGNIFICANT NEGATIVE RESULT.
It closes the question "can TGP produce MOND from axioms alone?" with NO.

NUMERICAL CONSISTENCY NOTE:
    a0_TGP predicted = c H0 / (2 pi) = {a0_tgp:.4e} m/s^2
    a0_obs           = {a0_obs:.4e} m/s^2
    agreement within 15%, but now purely phenomenological
    (not derived from TGP axioms).

ROADMAP:

  Option A (Theory extension path):
    gs67: Formalize Bridge (b) dual-scale extension of sek05.
          Derive a0 = c H0 / (2 pi) from geometric factor.
          Test against SPARC with smooth nu(y).

  Option B (Different phenomenon path):
    Re-focus galaxy program as "best-fit nu(y) phenomenology"
    (like gs10-gs49), abandoning unique-derivation claim.
    Focus TGP predictive power elsewhere:
      - Cluster dynamics (gs62 reinterpreted)
      - Solar system constraints
      - Cosmology (sek05 slow-roll inflation?)
      - Quantum sector (gs-QG program)

  Option C (Honest closure):
    Conclude that galaxy dynamics are NOT a unique TGP signature.
    Best case: consistent with observations under phenomenological
    smooth nu(y).  TGP neither predicts nor forbids MOND uniquely.
""")

header("gs66 complete.  Bridge (a) ruled out.  Galaxy program downgraded.")

print("""
FINAL POSITION OF TGP ON MOND (as of gs66):

  * MOND is CONSISTENT with TGP under phenomenological parameterization.
  * MOND is NOT UNIQUELY PREDICTED by TGP axioms (sek02+sek05+linear).
  * MOND requires EITHER theory extension (bridges b,c) OR phenomenology.
  * The numerical relation a0 = c H0 / (2 pi) remains empirically striking
    but lacks first-principles derivation.

The galaxy-scaling program transitions from
  "does TGP predict MOND?"  (answered: NO, in any linear sense)
to
  "what DOES TGP predict in galaxy/cluster dynamics?" (new program).
""")
