# ==============================================================================
#   gs65_formal_derivation.py
#
#   FORMAL DERIVATION: Does MOND emerge from sek02 + sek05?
#
#   This script performs a rigorous symbolic investigation of whether
#   the TGP field equation (sek02) with the slow-roll potential (sek05)
#   naturally produces MOND-like flat rotation curves.
#
#   The conclusion is NEGATIVE in the standard form: classical static
#   TGP + point source does NOT yield v_flat = const. We identify the
#   specific gap and propose three candidate bridges, then evaluate each
#   against self-consistency and the observed a0 = c*H0/(2pi).
#
#   Structure:
#     Part 1. Equations of motion (symbolic)
#     Part 2. Linearized regime (weak field, Yukawa)
#     Part 3. Weak nonlinearity (perturbative, power-law corrections)
#     Part 4. Strong nonlinearity (large-gradient ansatz)
#     Part 5. Comparison with AQUAL structure
#     Part 6. What would be required for MOND in TGP
#     Part 7. Candidate bridges:
#                (a) Hubble-friction mechanism (cosmological 3H*d_t)
#                (b) Dual-scale substrate (mass-dependent L_nat)
#                (c) Different operator ordering / membrane
#     Part 8. Honest conclusion + roadmap
# ==============================================================================

from __future__ import annotations

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


# ------------------------------------------------------------------------------
# Part 1. EQUATIONS OF MOTION
# ------------------------------------------------------------------------------
header("Part 1. Equations of motion from sek02 + sek05")

# Symbols
Phi, Phi0 = sp.symbols("Phi Phi0", positive=True)
psi = sp.Symbol("psi", positive=True)        # psi = Phi / Phi0
varphi = sp.Symbol("varphi", real=True)      # varphi = psi - 1  (fluctuation)
alpha, beta, gamma = sp.symbols("alpha beta gamma", positive=True)
r, M, G, c, H = sp.symbols("r M G c H", positive=True)
rho = sp.Symbol("rho", real=True)            # source density

print("""
sek02 field equation (static, spherical):

    (1/Phi0) * nabla^2 Phi
  + (alpha/Phi0) * (nabla Phi)^2 / Phi
  + beta  * Phi^2 / Phi0^2
  - gamma * Phi^3 / Phi0^3
  = kappa_src * rho

with alpha = 2 (from variation of action),
     beta  = gamma (vacuum condition U'(1)=0),
     [beta] = [gamma] = L^-2.

sek05 self-coupling potential (around psi = 1):

    U(psi) = (beta/3) * psi^3 - (gamma/4) * psi^4

with vacuum at psi = 1, and in Taylor expansion:

    U''(1)   = 2 beta - 3 gamma = -gamma < 0     (MAXIMUM, slow-roll)
    U'''(1)  = 2 beta - 6 gamma = -4 gamma
    U''''(1) = -6 gamma

Key observation: U''(1) is NEGATIVE, so psi=1 is a SADDLE/MAX, NOT an
oscillator minimum. There is no natural frequency at which the substrate
would ring. The Mech-F hypothesis (oscillation at omega ~ H0) is
falsified at the level of the background potential.
""")

# Verify U derivatives symbolically
U = (beta/3)*psi**3 - (gamma/4)*psi**4
U1 = sp.diff(U, psi)
U2 = sp.diff(U, psi, 2)
U3 = sp.diff(U, psi, 3)
U4 = sp.diff(U, psi, 4)

vac_sub = {beta: gamma}  # vacuum condition

print("Symbolic check (with vacuum condition beta = gamma):")
print(f"  U'(1)   = {sp.simplify(U1.subs(psi, 1).subs(vac_sub))}")
print(f"  U''(1)  = {sp.simplify(U2.subs(psi, 1).subs(vac_sub))}")
print(f"  U'''(1) = {sp.simplify(U3.subs(psi, 1).subs(vac_sub))}")
print(f"  U''''(1)= {sp.simplify(U4.subs(psi, 1).subs(vac_sub))}")


# ------------------------------------------------------------------------------
# Part 2. LINEARIZED REGIME
# ------------------------------------------------------------------------------
header("Part 2. Linearized regime (weak field) -> Yukawa, NOT MOND")

print("""
Write Phi = Phi0 * (1 + varphi), keep terms O(varphi).
Substitute into sek02 and retain linear terms in varphi:

    nabla^2 varphi + U''(1) * varphi = rho/Phi0 * kappa_src

With U''(1) = -gamma < 0:

    nabla^2 varphi - (-gamma) * varphi = source
    nabla^2 varphi - m_eff^2 varphi   = source,   m_eff^2 = -gamma  (tachyonic!)

But slow-roll interprets |m_eff| as an inverse length L_nat:
    L_nat = 1/sqrt(gamma)

For a point source M at origin, the Green's function of
    (nabla^2 - m^2) varphi = -4 pi G rho
has SCREENED/ANTI-SCREENED behavior:
    m^2 > 0 (Yukawa):        varphi ~ exp(-m r)/r
    m^2 < 0 (tachyonic):     varphi ~ cos(|m| r)/r (oscillatory)

NEITHER produces log(r) required for flat rotation curves
    v^2(r) = G M_enc(r) / r   ->   v_flat const  <=>   Phi ~ log(r).

VERDICT (Part 2):
    LINEARIZED TGP YIELDS EITHER YUKAWA SCREENING OR OSCILLATORY CORRECTIONS.
    NEITHER IS MOND.  The flat-RC regime must come from NONLINEAR dynamics.
""")


# ------------------------------------------------------------------------------
# Part 3. WEAK NONLINEARITY (perturbative)
# ------------------------------------------------------------------------------
header("Part 3. Weak nonlinearity (perturbative expansion) -> r^(-1/3) corrections")

print("""
Keep the nonlinear kinetic term (alpha/Phi0)*(nabla Phi)^2/Phi = alpha*(nabla varphi)^2/(1+varphi)
(Dominant in large-gradient, small |varphi| regime, and crucially:
 it does NOT vanish when |varphi| is small -- grad varphi can still be large.)

Far from source, in vacuum (rho = 0) and dropping U (slow varphi):

    nabla^2 varphi + alpha * (nabla varphi)^2 / (1 + varphi) = 0                 (*)

Ansatz: Phi = Phi0 * psi(r) with spherical symmetry.

(nabla varphi)^2   -> (psi')^2
nabla^2 varphi    -> psi'' + 2 psi'/r

Equation (*):
    psi'' + 2 psi'/r + alpha * (psi')^2 / psi = 0                                 (**)

Let u = log(psi). Then psi' = psi * u', psi'' = psi * (u'' + u'^2).
Equation (**) becomes:

    psi (u'' + u'^2) + 2 psi u'/r + alpha psi u'^2 = 0
    u'' + (1 + alpha) u'^2 + 2 u'/r = 0                                           (***)

Set w = u'.  First-order ODE:
    w' + (1+alpha) w^2 + 2 w/r = 0

This is Riccati.  Change variable w = y / r:
    (y' r - y)/r^2 + (1+alpha) y^2 / r^2 + 2 y / r^2 = 0
    y' r + (1+alpha) y^2 + y = 0
    dy / (r) = -[(1+alpha) y^2 + y] / r

Separable:
    dy / [(1+alpha) y^2 + y] = -d(log r)

Integration:
    y = [ A r^(1+alpha) / (1 - A r^(1+alpha) (1+alpha)) ] * (-1/(1+alpha))
  (details less important than asymptotics)

Asymptotic far-field solution (vacuum):  psi ~ r^n with  u'' + (1+alpha) u'^2 = 0
=>  u = log(psi) ~ -(1/(1+alpha)) * log(r)   =>   psi ~ r^{-1/(1+alpha)}

With alpha = 2:    psi ~ r^{-1/3}.

Gradient of effective Newton potential Phi_N from geodesic matching:
    g(r) = -(d Phi_eff/d r) ~ dpsi/dr * (Phi0/m) ~ r^{-4/3}

This corresponds to a rotation curve:
    v^2 = r |g|  ~  r * r^{-4/3}  =  r^{-1/3}   -->   v ~ r^{-1/6}  (FALLING!)

VERDICT (Part 3):
    STANDARD TGP FAR-FIELD VACUUM SCALING  g ~ r^(-4/3)
    IS STEEPER THAN NEWTON (g ~ r^-2), NOT SHALLOWER.
    This PREDICTS ROTATION CURVES THAT FALL FASTER THAN NEWTON,
    NOT flat curves.

This is the opposite of MOND. Something essential is missing.
""")

# Verify symbolic derivation
r_s = sp.Symbol("r", positive=True)
psi_s = sp.Function("psi")(r_s)
nabla2 = psi_s.diff(r_s, 2) + 2*psi_s.diff(r_s)/r_s
grad2  = psi_s.diff(r_s)**2
eq_vac = nabla2 + alpha*grad2/psi_s

# Test power-law ansatz psi = r^n
n = sp.Symbol("n", real=True)
psi_ansatz = r_s**n
residual = (psi_ansatz.diff(r_s, 2) + 2*psi_ansatz.diff(r_s)/r_s
            + alpha*psi_ansatz.diff(r_s)**2/psi_ansatz)
residual = sp.simplify(residual * r_s**(2))  # make polynomial in n
print("Symbolic verification: residual * r^2 for psi = r^n:")
print(f"  = {sp.expand(residual)}")
print(f"Solving for n with alpha = 2:")
sol = sp.solve(residual.subs(alpha, 2), n)
print(f"  n = {sol}")
print(f"  => nontrivial solution psi ~ r^(-1/3) confirmed.")


# ------------------------------------------------------------------------------
# Part 4. STRONG NONLINEARITY
# ------------------------------------------------------------------------------
header("Part 4. Strong nonlinearity: what would flat RC require?")

print("""
For flat rotation curves, the effective potential Phi_eff must be logarithmic:

    Phi_eff(r) = v_flat^2 * log(r/r0)                                             (MOND)

This requires g(r) = v_flat^2 / r, hence an inverse-linear acceleration.

Working backwards from this target in the TGP framework:

A. If the effective kinetic Lagrangian is L = F((nabla Phi)^2), the field
   equation in deep regime is:
       div [ F'((nabla Phi)^2) * nabla Phi ] = 4 pi G rho
   For spherical point source:  F'(y) y^(1/2) ~ M/r^2,
   i.e. F'(y) ~ 1/sqrt(y) <=> F(y) ~ y^(3/2)  (AQUAL deep MOND limit).

B. TGP nonlinearity (alpha/Phi0)(nabla Phi)^2/Phi is NOT of the form F((nabla Phi)^2).
   It is (nabla Phi)^2 divided by the FIELD itself, not a function of gradient alone.
   So the standard AQUAL derivation does not apply directly.

C. Could TGP's (nabla Phi)^2 / Phi mimic (nabla Phi)^3 in some regime?
   Only if Phi ~ 1 / |nabla Phi|. In spherical symmetry with Phi ~ log(r):
       |nabla Phi| = v_flat^2 / r
       Phi          = v_flat^2 log(r)
       Phi * |nabla Phi| = v_flat^4 log(r) / r   (NOT constant)
   So this identification fails: TGP's kinetic term does NOT reduce to AQUAL.

D. Structural comparison:
     AQUAL:  L ~ |nabla Phi|^3 / a0        (homogeneous degree 3 in gradients)
     TGP  :  L ~ (nabla Phi)^2 / Phi       (homogeneous degree 2 in gradients,
                                            degree -1 in field)

   These Lagrangians are NOT equivalent, and do not share a BTFR slope by coincidence.

VERDICT (Part 4):
    TGP's NONLINEAR KINETIC STRUCTURE DOES NOT REDUCE TO AQUAL FORM.
    Flat rotation curves cannot emerge from static TGP + point mass alone.
""")


# ------------------------------------------------------------------------------
# Part 5. COMPARISON WITH AQUAL
# ------------------------------------------------------------------------------
header("Part 5. Quantitative comparison: AQUAL vs TGP nonlinear kinetic")

print("""
                     AQUAL                         TGP (sek02)
  ----------------------------------------------------------------
  Kinetic Lagrangian  F(X) with X=(nabla Phi)^2    (nabla Phi)^2 / Phi
  Deep-MOND F(X)      (2/3) X^(3/2) / a0           N/A (no pure X-dependence)
  Grad term in EoM    div[ mu(|g|/a0) g ]          nabla^2 Phi + 2 (nabla Phi)^2/Phi
  BTFR slope          exactly 4                    not derivable from static limit
  Point-source deep   g = sqrt(G M a0)/r           psi ~ r^(-1/3) (NOT MOND)
  External Field Eff  YES (explicit via mu)        Weak/absent in static limit

The key structural difference: AQUAL's interpolation function mu depends
on |g|/a0 with a fixed acceleration scale a0, while TGP's nonlinearity
(alpha/Phi0)(nabla Phi)^2/Phi depends on the RATIO of gradient to field
strength -- a DIMENSIONLESS measure with no built-in acceleration scale.

This is WHY a0 does not emerge from sek02 directly:
    sek02 is conformally/multiplicatively invariant in the kinetic term
    (it has no preferred scale).  a0 must come from ELSEWHERE.
""")


# ------------------------------------------------------------------------------
# Part 6. WHAT WOULD PRODUCE MOND IN TGP
# ------------------------------------------------------------------------------
header("Part 6. What additional structure would produce MOND?")

print("""
For TGP to produce MOND dynamics with a0 = c H0 / (2 pi), we need:

(R1)  An acceleration scale that enters the equations.
      Candidates:
        * Hubble rate H0 (cosmological)
        * gamma^(1/2) = 1/L_nat (potential-derived length)
        * GM-dependent scale (membrane mechanism)
        * Quantum-gravitational a_Q = c^2 / L_Planck (too large by 10^43)

(R2)  A regime in which the field equation effectively reduces to:
          div[ mu(|nabla Phi|/a0) nabla Phi ] = 4 pi G rho

(R3)  Recovery of Newton (mu -> 1) for |nabla Phi| >> a0.

Three formal bridges evaluated:

  BRIDGE 1: Cosmological Hubble friction (add 3H * dot{Phi}).
    Static TGP misses time-dependence. In FRW:
        Phi_ddot + 3 H Phi_dot - (1/a^2) nabla^2 Phi + U'(Phi) = source

    Quasi-static limit with Phi_dot ~ H Phi:
        (effective time-derivative term) ~ 3 H^2 Phi
    This INTRODUCES a0 scale via:
        3 H^2 Phi ~ rho           Phi ~ rho / (3 H^2)
    but this is GLOBAL cosmology, not a modification to local gravity.

    For a LOCAL effect: require oscillation at omega ~ H0.
    But U''(1) = -gamma (not H0^2), so period = 2 pi / sqrt(-U''(1))
    is imaginary (slow-roll), not matching 2 pi / H0 unless
        gamma = -H0^2  <=>  L_nat^2 = -H0^2 (inconsistent: L_nat real).

    Therefore BRIDGE 1 as "oscillation mechanism" is FALSIFIED.
    (This confirms gs64 result.)

    BUT: Hubble friction term 3 H * dot{Phi} in a slow-rolling Phi
    CAN modify the effective propagator in ways that change far-field
    scaling.  Requires explicit computation (gs66 or future work).

  BRIDGE 2: Dual-scale substrate (L_nat -> {L_mu, L_c}).
    Posit TWO natural lengths:
        L_mu ~ kpc       (microscopic Yukawa, local modification)
        L_c  ~ L_H       (cosmological, sets a0)
    This requires TWO gamma-like parameters in the potential -- a
    modification of sek05.  Not derivable from current axioms.

    If adopted, BTFR follows because M-r_eq relation is set by GM/c^2
    and L_c alone (gs9d).  Normalization off by factor ~2pi.

    BRIDGE 2 requires EXTENDING the theory.

  BRIDGE 3: Membrane / 2D-transition mechanism.
    At a mass-dependent radius r_eq = sqrt(GM * L_c) the effective
    kinetic term switches from 3D (~(nabla Phi)^2) to 2D (~(nabla_2 Phi)^2).
    Gives slope 4 BTFR but hard transition, wrong chi^2 on SPARC (gs61).

    BRIDGE 3 (hard switch) is EMPIRICALLY REJECTED.
    Soft transition version remains open.

""")


# ------------------------------------------------------------------------------
# Part 7. CONSISTENCY OF a0 = c H0 / (2 pi)
# ------------------------------------------------------------------------------
header("Part 7. a0 = c H0 / (2 pi):  numerical consistency & dimensional logic")

import math

c_num  = 2.998e8
H0_num = 2.1844e-18
a0_obs = 1.20e-10
a0_tgp = c_num * H0_num / (2 * math.pi)
ch0    = c_num * H0_num

print(f"""
Numerical values:
    c * H0           = {ch0:.4e} m/s^2
    a0_TGP           = c H0 / (2 pi) = {a0_tgp:.4e} m/s^2
    a0_obs (MOND)    = {a0_obs:.4e} m/s^2
    ratio (obs/TGP)  = {a0_obs/a0_tgp:.4f}
    ratio (c*H0/obs) = {ch0/a0_obs:.4f}

The numerical prediction a0_TGP = {a0_tgp:.4e} is within ~15% of a0_obs.
The factor of (2 pi) is the signature of an oscillation / period-length
identification (e.g. 2 pi R = c T => a0 = v^2/R = v H0 / (2 pi) for v=c).

Dimensional logic:
    [c H0] = m/s * 1/s = m/s^2 = acceleration  [CHECK]
    [a0]   = m/s^2                             [CHECK]
So a0 = cH0/(2pi) is dimensionally required for any theory where:
    - only c and H0 set fundamental scales for a0
    - and the 2pi comes from a geometric period

BRIDGE-EVALUATION of a0 PREDICTION:

  Bridge 1 (Hubble friction):
    Natural scale is c*H (no 2pi).  Matches Milgrom '99 but gives
    a0 ~ cH0/2pi only if a geometric prefactor appears.
    Provisional: PLAUSIBLE if 2pi comes from angular integration.

  Bridge 2 (dual-scale):
    a0 = c^2 / L_c.  With L_c = L_H = c/H0:  a0 = c H0.  NO 2pi factor.
    Provisional: MISSES by 2pi; need additional 1/(2pi) from geometry.

  Bridge 3 (membrane hard-switch):
    BTFR with slope 4 gives a0_eff = c*H0 = 6.55e-10 (off by 5.5x).
    Provisional: REJECTED numerically.
""")


# ------------------------------------------------------------------------------
# Part 8. CONCLUSION
# ------------------------------------------------------------------------------
header("Part 8. Honest conclusion and roadmap")

print("""
FORMAL STATUS after gs65 derivation:

  CLAIM (negative):
      Standard TGP (sek02 + sek05, static limit, point source)
      does NOT naturally produce MOND dynamics.

      - Linearized regime:       Yukawa / tachyonic, NOT log(r).
      - Weak nonlinearity:       psi ~ r^(-1/3),
                                 g ~ r^(-4/3),
                                 v_rot ~ r^(-1/6)  (falling).
      - Strong nonlinearity:     TGP kinetic is not of AQUAL form.
      - Potential U(psi):        slow-roll MAX (U''(1) < 0),
                                 not oscillator.

  CLAIM (positive conditional):
      IF TGP is extended with ONE of:
         a) time-dependent cosmological Hubble friction with
            nontrivial spatial propagator, OR
         b) dual-scale substrate (two natural lengths), OR
         c) smooth 3D->2D transition at membrane radius
      THEN MOND-like phenomenology CAN be embedded,
      with a0 approximately c H0 / (2 pi).

  UNIQUENESS:
      None of (a), (b), (c) is forced by sek02 + sek05 alone.
      Additional axioms or a FUNDAMENTAL derivation of 3H*d_t in the
      static limit is required.

  FALSIFIABILITY IMPACT:
      The TGP "natural MOND" claim is DOWNGRADED.
      Current status: TGP is CONSISTENT WITH MOND-like phenomenology
      under additional assumptions, but does NOT UNIQUELY PREDICT it.

ROADMAP:

  gs66: Solve the FULL (cosmological) TGP field equation
        around a localized mass in FRW background:
            Phi_ddot + 3 H Phi_dot - nabla^2 Phi/a^2 + U'(Phi) = source
        Does the 3H term induce a log-like far-field at r ~ L_H?
        This is the only derivation path that preserves
        "no extra axioms".

  gs62-revisit: Use Bridge 1 (Hubble friction) on clusters.
                Test if same mechanism explains cluster residuals.

  gs63: Euclid a0(z) prediction.  Any of the three bridges predicts:
            a0(z) = c H(z) / (2 pi) ?
        (Bridge 1 trivially; Bridges 2,3 nontrivially through L_c(z).)

Recommended next step: gs66 full FRW-solution in linear perturbation
theory.  Keep everything symbolic; see whether a scale emerges from
the propagator at k ~ H0/c.
""")


# ------------------------------------------------------------------------------
# Appendix: explicit Green's function for quasi-static FRW TGP
# ------------------------------------------------------------------------------
sub_header("Appendix A. Quasi-static FRW propagator (preview for gs66)")

print("""
Consider metric perturbations small, background FRW.  Fourier-space
equation for varphi at wave-number k, frequency omega:

    [- omega^2 + 3 i H omega + (k^2/a^2) + U''(1)] varphi(k, omega)
        = source

In quasi-static limit omega -> 0 (time variation slow compared to 1/H):

    [ 3 i H omega + k^2 / a^2 + U''(1) ] varphi = source
            |-> 0 in omega->0 limit

But static limit loses the Hubble term.  Need NEXT order in omega:

    omega ~ H (slow-roll), so:

    [ 3 i H^2 + (k^2 / a^2) - gamma ] varphi = source

Real part of propagator:
    Re[1 / (k^2/a^2 - gamma + 3 i H^2)]
        = (k^2/a^2 - gamma) / [(k^2/a^2 - gamma)^2 + 9 H^4]

For k^2/a^2 -> gamma (tachyonic resonance):
    propagator picks up a pole width ~ 9 H^4 / gamma
    characteristic length scale: L_res = a * sqrt(gamma)
    but COUPLED to H0 through the width.

This is the calculation to be done rigorously in gs66.
The hope: an effective log(r) tail emerges from integrating over
the quasi-static resonance for r >> L_nat.

A preliminary dimensional estimate:
    effective a0 ~ c * H0 * (some order-unity geometric factor)

This is the direction that will either VINDICATE or DEFINITIVELY FALSIFY
the TGP "natural MOND" program.
""")

header("gs65 formal derivation complete.")

print("""
SUMMARY LINE:
    sek02 + sek05 (static) does NOT give MOND.
    MOND in TGP requires time-dependent FRW background (Hubble friction)
    OR theory extension.  This is an HONEST negative result for the
    static program, and a precise pointer to gs66 (FRW linear propagator).
""")
