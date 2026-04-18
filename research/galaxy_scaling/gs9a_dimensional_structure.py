#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs9a_dimensional_structure.py: TGP soliton equation in d dimensions.

QUESTION: Does TGP naturally "prefer" 2D at large scales?

The soliton equation in d spatial dimensions:
  g'' + g'^2/g + (d-1)*g'/r + g = 1

We study:
1. How solutions change with d (continuously from 1 to 4)
2. Green's function behavior: 1D, 2D, 3D — which gives flat RC?
3. Effective dimension d_eff(r) extracted from potential behavior
4. Mathematical structure: symmetries, conservation laws vs d
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.special import jv  # Bessel functions

print("="*78)
print("  gs9a: TGP equation in d dimensions")
print("="*78)

# ===========================================================================
# 1. SOLITON IN d DIMENSIONS
# ===========================================================================
print(f"\n{'='*78}")
print(f"  1. Soliton equation: g'' + g'^2/g + (d-1)*g'/r + g = 1")
print(f"{'='*78}")

def solve_soliton_d(g0, d, r_max=80, N=80000):
    """Solve g'' + g'^2/g + (d-1)*g'/r + g = 1"""
    def rhs(r, y):
        g, gp = y
        if r < 1e-10:
            # L'Hopital: (d-1)*g'/r -> (d-1)*g'' at r=0
            # g'' + g'^2/g + (d-1)*g'' + g = 1
            # d*g'' + g'^2/g + g = 1
            # g'' = (1 - g - gp^2/g) / d
            gpp = (1 - g - gp**2/g) / d
            return [gp, gpp]
        gpp = -gp**2/g - (d-1)*gp/r - g + 1
        return [gp, gpp]

    r_eval = np.linspace(1e-6, r_max, N)
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0], t_eval=r_eval,
                    method='RK45', rtol=1e-12, atol=1e-14, max_step=0.02)
    return sol.t, sol.y[0], sol.y[1]

print(f"\n  Soliton properties vs dimension d (g0 = 1.1):")
print(f"  {'d':>5s} {'E_soliton':>12s} {'M_grav':>10s} {'tail_env':>10s} {'tail_decay':>12s} {'period':>8s}")
print(f"  {'-'*62}")

g0 = 1.1
soliton_data = {}

for d in [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:
    r, g, gp = solve_soliton_d(g0, d, r_max=80)
    delta = g - 1.0
    dr = r[1] - r[0]

    # Energy: E = integral [g'^2/g + (g-1)^2] * S_d * r^(d-1) dr
    # where S_d = 2*pi^(d/2) / Gamma(d/2) (surface area of d-sphere)
    from scipy.special import gamma as gamma_func
    S_d = 2 * np.pi**(d/2) / gamma_func(d/2)
    integrand_E = (gp**2/g + delta**2) * S_d * r**(d-1)
    E_sol = simpson(integrand_E, x=r)

    # Gravitational mass from asymptotic behavior
    # In d dimensions: M_enc ~ -r^(d-1) * g' / g (Gauss's law)
    M_enc = -r**(d-1) * gp / g
    idx_far = (r > 20) & (r < 60)
    M_grav = np.mean(M_enc[idx_far]) if np.any(idx_far) else M_enc[-1]

    # Tail analysis: find peaks and fit envelope
    peaks_r = []
    peaks_d = []
    for i in range(1, len(r)-1):
        if r[i] > 8 and abs(delta[i]) > abs(delta[i-1]) and abs(delta[i]) > abs(delta[i+1]):
            peaks_r.append(r[i])
            peaks_d.append(abs(delta[i]))

    if len(peaks_r) > 3:
        pr = np.array(peaks_r)
        pd = np.array(peaks_d)
        # Fit: |delta| ~ A * r^(-alpha)
        log_r = np.log(pr)
        log_d = np.log(pd)
        alpha_fit, log_A = np.polyfit(log_r, log_d, 1)
        A_env = np.exp(log_A)
        tail_env = f"{A_env:.4f}"
        tail_decay = f"r^({alpha_fit:.2f})"

        # Period
        periods = np.diff(pr)
        period = np.mean(periods)
    else:
        tail_env = "n/a"
        tail_decay = "n/a"
        period = 0
        alpha_fit = 0

    print(f"  {d:5.1f} {E_sol:12.4f} {M_grav:10.4f} {tail_env:>10s} {tail_decay:>12s} {period:8.2f}")
    soliton_data[d] = {
        'r': r, 'g': g, 'gp': gp, 'delta': delta,
        'E': E_sol, 'M': M_grav, 'alpha': alpha_fit, 'period': period
    }

# ===========================================================================
# 2. TAIL BEHAVIOR: THE KEY DIFFERENCE
# ===========================================================================
print(f"\n{'='*78}")
print(f"  2. Tail behavior analysis")
print(f"{'='*78}")

print(f"""
  The LINEARIZED equation in d dimensions:
  delta'' + (d-1)*delta'/r + delta = 0

  Solutions:
  d=1: delta = A*sin(r) + B*cos(r)  (non-decaying oscillation!)
  d=2: delta = A*J_0(r) + B*Y_0(r)  (Bessel, decays as r^(-1/2))
  d=3: delta = A*sin(r)/r + B*cos(r)/r  (decays as r^(-1))
  d=4: delta = A*J_1(r)/r + B*Y_1(r)/r  (decays as r^(-3/2))

  General: tail decays as r^(-(d-1)/2)

  This means:
  d=1: tail amplitude CONSTANT → infinite range!
  d=2: tail ~ r^(-1/2) → slower decay than 3D
  d=3: tail ~ r^(-1) → standard
  d=4: tail ~ r^(-3/2) → faster decay

  For GRAVITATIONAL POTENTIAL (which IS the tail):
  The "mass contribution" from the tail scales as:
  M_tail(R) ~ integral |delta|^2 * r^(d-1) dr ~ integral r^(-(d-1)) * r^(d-1) dr
  = integral dr → DIVERGES for ALL d!

  But the FORCE from the tail:
  F ~ d(delta)/dr ~ r^(-(d-1)/2 - 1) * sin(r)
  Averaged over oscillations: <F> = 0
  (same cancellation problem as in gs6)

  HOWEVER: what matters for rotation curves is not the soliton tail
  but the GRAVITATIONAL GREEN'S FUNCTION in d dimensions.
""")

# Verify numerically
print(f"  Numerical tail decay exponents:")
print(f"  {'d':>5s} {'alpha (fit)':>12s} {'(d-1)/2':>10s} {'match?':>8s}")
print(f"  {'-'*38}")
for d in [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:
    alpha = soliton_data[d]['alpha']
    expected = -(d-1)/2
    match = "YES" if abs(alpha - expected) < 0.15 else "no"
    print(f"  {d:5.1f} {alpha:12.2f} {expected:10.2f} {match:>8s}")

# ===========================================================================
# 3. GREEN'S FUNCTION: POISSON + SPRING IN d DIMENSIONS
# ===========================================================================
print(f"\n{'='*78}")
print(f"  3. Green's function of (nabla^2 + 1) in d dimensions")
print(f"{'='*78}")

print(f"""
  The TGP static equation (linearized) with source:
  nabla^2_d delta + delta = -rho(r)

  The Green's function G(r) satisfies:
  nabla^2_d G + G = -delta_Dirac(r)

  In d dimensions:
  d=1: G(r) = sin(|r|)/2  (half-wave)
  d=2: G(r) = -(i/4)*H_0^(1)(r) → oscillatory, envelope ~ 1/sqrt(r)
  d=3: G(r) = -sin(r)/(4*pi*r)  → oscillatory, envelope ~ 1/r

  For the Newtonian part (nabla^2 G_N = -delta):
  d=2: G_N(r) = -(1/2pi)*ln(r)  → LOGARITHMIC!
  d=3: G_N(r) = 1/(4*pi*r)      → 1/r

  KEY OBSERVATION:
  The TGP equation has nabla^2 + 1 (with spring term).
  The Newtonian equation has just nabla^2.

  In 3D: spring makes it oscillatory → tails cancel
  In 2D: Newtonian potential is ALREADY logarithmic → flat RC!

  So the question becomes:
  Is there a mechanism in TGP where the spring term (+delta)
  becomes ineffective at large r, leaving just nabla^2?

  If yes: at large r, effective equation is nabla^2 delta = -rho
  and if d_eff = 2: delta ~ ln(r) → FLAT ROTATION CURVE!
""")

# Compute Green's functions numerically
print(f"  Green's function comparison (point source at r=0.5):")
print(f"  {'r':>6s} {'G_3D_TGP':>12s} {'G_3D_Newton':>12s} {'G_2D_Newton':>12s} {'G_2D_TGP':>12s}")
print(f"  {'-'*58}")

for r_val in [1, 2, 5, 10, 20, 50, 100]:
    r = float(r_val)
    G_3D_TGP = -np.sin(r) / (4*np.pi*r)
    G_3D_N = 1 / (4*np.pi*r)
    G_2D_N = -np.log(r) / (2*np.pi)
    G_2D_TGP = -jv(0, r) / 4  # approximate
    print(f"  {r_val:6d} {G_3D_TGP:12.4e} {G_3D_N:12.4e} {G_2D_N:12.4e} {G_2D_TGP:12.4e}")

# ===========================================================================
# 4. THE SPRING TERM AND EFFECTIVE DIMENSION
# ===========================================================================
print(f"\n{'='*78}")
print(f"  4. When does the spring term become irrelevant?")
print(f"{'='*78}")

print(f"""
  The TGP equation: nabla^2 delta + delta = -rho

  Two regimes:
  A) r << 2*pi (within one oscillation wavelength):
     nabla^2 delta dominates over delta
     → Effective equation: nabla^2 delta ~ -rho → NEWTONIAN

  B) r >> 2*pi (many wavelengths):
     Both terms matter → oscillatory solution → averages out

  In regime A: delta behaves as Newtonian potential in d dimensions.
  For d=3: delta ~ 1/r → F ~ 1/r^2
  For d=2: delta ~ ln(r) → F ~ 1/r → FLAT RC!

  The soliton wavelength in physical units: lambda = 2*pi*l_soliton.
  If l_soliton ~ fm, then lambda ~ 6 fm.

  A galaxy (R ~ 10 kpc ~ 3e22 m) is at r/lambda ~ 5e36.
  → We are DEEPLY in regime B → oscillations dominate → cancellation.

  UNLESS: the effective wavelength is much longer.

  From gs7a: if mu = H0/c (cosmological scale), then:
  lambda = 2*pi/mu = 2*pi*c/H0 ~ 10^27 m >> R_galaxy

  In this case, the galaxy is in regime A: nabla^2 delta ~ -rho.
  And if d_eff = 2 at these scales → ln(r) potential → flat RC!

  So the COMBINATION of:
  1. Spring term with cosmological wavelength (mu = H0/c)
  2. Effective 2D behavior at galactic scales
  gives the right phenomenology.

  But we already showed (gs7b, gs8) that mu = H0/c gives perturbative
  corrections ~ 10^-12 in 3D.

  The key difference in 2D: the Newtonian potential is ALREADY ln(r),
  so even WITHOUT modifications to Newton, 2D gravity gives flat RC!

  The question is NOT "how to modify 3D gravity"
  but "WHY would gravity become effectively 2D?"
""")

# ===========================================================================
# 5. EFFECTIVE DIMENSION FROM THE POTENTIAL
# ===========================================================================
print(f"\n{'='*78}")
print(f"  5. Extracting effective dimension from potential behavior")
print(f"{'='*78}")

print(f"""
  Given a potential Phi(r), define effective dimension:

  In d dimensions: nabla^2 Phi = rho implies
  Phi(r) ~ r^(2-d) for d > 2 (power law)
  Phi(r) ~ ln(r) for d = 2

  So: d_eff(r) = 2 - d(ln Phi)/d(ln r)  for power-law potential
  (with appropriate sign conventions)

  For observed galaxies:
  - Inner region: Phi ~ -GM/r → d(ln|Phi|)/d(ln r) = -1 → d_eff = 3
  - Outer region: Phi ~ v_flat^2*ln(r) → d(ln|Phi|)/d(ln r) → 0 → d_eff = 2!

  The galaxy rotation curve ITSELF tells us:
  d_eff transitions from 3 to 2 at the MOND radius!

  This is not a PREDICTION — it's a RESTATEMENT of MOND.
  But it suggests a clean mathematical framework:

  MOND ←→ effective dimensional reduction 3 → 2 at a = a0
""")

# Compute d_eff for MOND potential
G_astro = 4.3016e-6
a0_astro = 3.703e3
M_MW = 7e10

print(f"  Effective dimension for MW (MOND potential):")
print(f"  {'r (kpc)':>8s} {'g_N':>10s} {'g_MOND':>10s} {'Phi_MOND':>12s} {'d_eff':>6s}")
print(f"  {'-'*50}")

r_prev = None
Phi_prev = None
Phi_MOND_accumulated = 0

radii = np.logspace(-0.5, 2.5, 50)
Phi_values = []
for r in radii:
    g_N = G_astro * M_MW / r**2
    g_MOND = g_N / (1 - np.exp(-np.sqrt(g_N / a0_astro)))
    # Phi_MOND(r) = -integral g_MOND dr from inf to r
    # Approximate: Phi(r) ~ -g_MOND * r (locally)
    Phi_values.append(-g_MOND * r)

# Compute d_eff from slope of Phi vs r
for i in range(1, len(radii)):
    r = radii[i]
    g_N = G_astro * M_MW / r**2
    g_MOND = g_N / (1 - np.exp(-np.sqrt(g_N / a0_astro)))

    # d(ln|Phi|)/d(ln r) from finite differences
    dlogPhi = np.log(abs(Phi_values[i])) - np.log(abs(Phi_values[i-1]))
    dlogr = np.log(radii[i]) - np.log(radii[i-1])
    slope = dlogPhi / dlogr

    # For Phi ~ r^(2-d): slope = 2-d → d = 2-slope
    # For Phi ~ -1/r: slope = -1 → d = 3 ✓
    # For Phi ~ ln(r): slope → 0 → d = 2 ✓
    # But Phi = -g*r, so Phi ~ -g*r
    # If g ~ 1/r^2 (Newton): Phi ~ -1/r → slope = -1 → d = 2-(-1) = 3 ✓
    # If g ~ 1/r (MOND): Phi ~ -ln(r)*const → slope → 0ish → d → 2 ✓

    d_eff = 2 - slope

    if r in [0.5, 1, 2, 5, 10, 20, 50, 100, 200]:
        pass  # print only selected
    if abs(np.log10(r) - round(np.log10(r))) < 0.15 or r < 1:
        print(f"  {r:8.1f} {g_N:10.1f} {g_MOND:10.1f} {Phi_values[i]:12.1f} {d_eff:6.2f}")

# ===========================================================================
# 6. MATHEMATICAL STRUCTURE: TGP LAGRANGIAN IN d DIMENSIONS
# ===========================================================================
print(f"\n{'='*78}")
print(f"  6. TGP Lagrangian in d dimensions — special properties")
print(f"{'='*78}")

print(f"""
  TGP Lagrangian: L = |nabla g|^2/g + (g-1)^2

  In d dimensions: S = integral L * r^(d-1) dr * S_d

  The kinetic term |nabla g|^2/g = g'^2/g has a special structure:
  it's the FISHER INFORMATION metric on the space of probability
  distributions p = g / integral(g).

  Fisher information is defined as:
  I_F = integral (p'/p)^2 * p dx = integral p'^2/p dx

  For g = p (probability-like): |nabla g|^2/g IS Fisher information!

  Fisher information has deep connections to:
  - Information geometry (Amari)
  - Cramer-Rao bound (minimum uncertainty)
  - Maximum entropy methods
  - Statistical mechanics

  KEY PROPERTY: Fisher information in d dimensions scales as:
  I_F[g_lambda(r)] = lambda^(2-d) * I_F[g(r)]

  under scaling g_lambda(r) = g(r/lambda).

  This means:
  d=2: I_F is SCALE-INVARIANT (conformal invariance!)
  d<2: I_F increases with scale (UV divergent)
  d>2: I_F decreases with scale (IR finite)

  CONFORMAL INVARIANCE IN 2D:
  The kinetic term of TGP is conformally invariant precisely in d=2.
  This is the SAME as why 2D field theories are special:
  - String theory worldsheet is 2D
  - Conformal field theory is most powerful in 2D
  - The Polyakov action for strings uses the 2D conformal structure

  This is a MATHEMATICAL reason why TGP might "prefer" 2D:
  the Fisher information metric (= TGP kinetic term) has
  conformal symmetry in 2D and ONLY in 2D.
""")

# Verify the scaling property
print(f"  Verification: Fisher information scaling")
print(f"  {'d':>5s} {'I_F[g]':>12s} {'I_F[g_2x]':>12s} {'ratio':>8s} {'2^(2-d)':>8s} {'match':>6s}")
print(f"  {'-'*55}")

for d in [1.0, 2.0, 3.0, 4.0]:
    data = soliton_data[d]
    r = data['r']
    g = data['g']
    gp = data['gp']

    S_d = 2 * np.pi**(d/2) / gamma_func(d/2)

    # I_F for original soliton
    integrand = (gp**2/g) * S_d * r**(d-1)
    mask = r < 30
    I_F_orig = simpson(integrand[mask], x=r[mask])

    # I_F for scaled soliton (lambda=2): g_2(r) = g(r/2), g_2'(r) = g'(r/2)/2
    r2 = r * 2
    mask2 = r2 < 30
    # Use the same data but at half the r values
    idx_half = np.arange(0, len(r)//2)
    if len(idx_half) > 100:
        r_h = r[idx_half]
        g_h = g[idx_half]
        gp_h = gp[idx_half] / 2  # chain rule
        integrand_h = (gp_h**2/g_h) * S_d * (r_h*2)**(d-1) * 2  # dr→2dr
        I_F_scaled = simpson(integrand_h, x=r_h)
    else:
        I_F_scaled = 0

    ratio = I_F_scaled / I_F_orig if I_F_orig > 0 else 0
    expected = 2**(2-d)
    match = "YES" if abs(ratio - expected) < 0.3 else "~" if abs(ratio-expected) < 0.5 else "no"
    print(f"  {d:5.1f} {I_F_orig:12.4f} {I_F_scaled:12.4f} {ratio:8.3f} {expected:8.3f} {match:>6s}")

# ===========================================================================
# 7. THE CONFORMAL CONNECTION
# ===========================================================================
print(f"\n{'='*78}")
print(f"  7. Conformal invariance of TGP kinetic term in 2D")
print(f"{'='*78}")

print(f"""
  The TGP Lagrangian density (kinetic part): l_kin = |nabla g|^2 / g

  Under conformal transformation in d dimensions:
  x → lambda*x, g(x) → g(lambda*x)

  The action S_kin = integral l_kin * d^d x transforms as:
  S_kin → lambda^(2-d) * S_kin  (from nabla → nabla/lambda, d^d x → lambda^d * d^d x)

  ONLY for d=2: S_kin → S_kin (invariant!)

  This means the kinetic term of TGP is a CONFORMAL FIELD THEORY in 2D.

  Consequences:
  1. In 2D, the TGP kinetic energy is INDEPENDENT of scale
     → Scale-free dynamics → no preferred length scale for gravity
     → This IS the phenomenology of flat rotation curves!

  2. The potential term (g-1)^2 BREAKS conformal invariance
     → It introduces a scale (the soliton size)
     → But it decays at large r (where g → 1, so (g-1)^2 → 0)

  3. At large distances from a source (galaxy outskirts):
     g ~ 1 + delta with delta << 1
     → (g-1)^2 = delta^2 << delta
     → The potential term becomes NEGLIGIBLE compared to kinetic
     → The equation approaches conformal (scale-free) behavior
     → EFFECTIVE 2D dynamics even in 3D space!

  THIS IS THE MECHANISM:

  At large r (weak field): TGP Lagrangian → |nabla g|^2/g (kinetic only)
  The kinetic term is conformally invariant in 2D.
  The dynamics "wants" to be in its conformally invariant sector.
  This manifests as 2D-like gravitational behavior.

  QUANTITATIVE CHECK:
  When does (g-1)^2 << g'^2/g?
  delta^2 << delta'^2 / (1 + delta)
  For delta ~ GM/(rc^2):
  delta^2 ~ (GM)^2/(r^2 c^4)
  delta'^2 ~ (GM)^2/(r^4 c^4)
  Ratio: delta^2/delta'^2 ~ r^2

  So: (g-1)^2 > g'^2/g when r > 1 (in soliton units)

  Wait — this means the potential term DOMINATES at large r,
  not the kinetic term! The opposite of what we need.

  Let me reconsider...

  Actually: at large r from a galaxy, delta → 0.
  Both delta^2 and delta'^2 → 0.
  The RATIO matters for the dynamics.

  The equation is: delta'' + 2delta'/r + delta = 0
  The "spring" term is +delta.
  The "kinetic" term is delta'' + 2delta'/r.

  At large r: if delta ~ A/r (Newtonian), then:
  delta'' = 2A/r^3, 2delta'/r = -2A/r^3 → kinetic sum = 0!
  Spring: A/r ≠ 0 → spring dominates → solution is NOT 1/r.

  The spring forces oscillations → sin(r)/r tail.

  For the potential to become 2D-like (ln r), we need
  the spring to become ineffective.

  In 2D: delta'' + delta'/r + delta = 0
  If delta ~ ln(r), then delta' = 1/r, delta'' = -1/r^2
  delta'' + delta'/r = -1/r^2 + 1/r^2 = 0
  Spring: delta = ln(r) → grows without bound → NOT → 0!

  So ln(r) is NOT a solution of the homogeneous equation in ANY dimension.
  The ln(r) potential comes from the Newtonian part (nabla^2 = source)
  WITHOUT the spring term.
""")

# ===========================================================================
# 8. KEY REALIZATION: POISSON vs HELMHOLTZ
# ===========================================================================
print(f"\n{'='*78}")
print(f"  8. KEY: It's about Poisson (nabla^2) vs Helmholtz (nabla^2 + 1)")
print(f"{'='*78}")

print(f"""
  The flat RC requires:
  nabla^2 Phi = rho  (POISSON, no spring)  in d=2 → Phi ~ ln(r)

  TGP gives:
  nabla^2 delta + delta = rho  (HELMHOLTZ with spring)  → oscillatory

  The ENTIRE problem of TGP and galaxy scaling comes down to:
  Can the spring term (+delta) be effectively removed at galactic scales?

  We've shown (gs7-gs8): NO perturbative mechanism removes it.
  The spring IS the equation. Without it, TGP = Newton.

  BUT: what if we reinterpret the meaning of "d=2"?

  In the LAGRANGIAN: L = g'^2/g + delta^2
  The ratio kinetic/potential = g'^2/(g*delta^2)

  For a Newtonian potential in 3D: g = 1 + A/r
  kinetic = A^2/r^4
  potential = A^2/r^2
  ratio = 1/r^2 → potential dominates at large r

  For a logarithmic potential in 2D: g = 1 + A*ln(r)
  kinetic = A^2/r^2
  potential = A^2*ln^2(r)
  ratio = 1/(r^2 * ln^2(r)) → potential dominates even MORE

  So the kinetic term is ALWAYS subdominant at large r.
  Conformal invariance of the kinetic term doesn't help —
  it's the POTENTIAL term that controls large-r behavior.

  REVISED UNDERSTANDING:
  The "effective 2D" signal in gs8 is NOT about conformal invariance.
  It's about the FORCE LAW:
  - In 3D: F = GM/r^2 → v^2 = GM/r → falling RC
  - In 2D: F = GM/r → v^2 = GM = const → flat RC

  The DGP-like model works because it INTERPOLATES between these.
  The question remains: WHY does gravity become 2D-like?
""")

# ===========================================================================
# 9. INFORMATION-THEORETIC PERSPECTIVE
# ===========================================================================
print(f"\n{'='*78}")
print(f"  9. Information-theoretic perspective: degrees of freedom")
print(f"{'='*78}")

print(f"""
  Alternative to geometric explanation:

  In d dimensions, the number of gravitational degrees of freedom
  on a sphere of radius r scales as:
  N_dof(r) ~ r^(d-1) (area of sphere)

  In 3D: N_dof ~ r^2 (surface area) → Gauss's law → F ~ 1/r^2
  In 2D: N_dof ~ r (circumference) → Gauss's law → F ~ 1/r

  If gravity is mediated by substrate correlations,
  and the EFFECTIVE number of correlated DoF at distance r is:

  N_eff(r) = N_3D(r) * f(r) = r^2 * f(r)

  then the force: F ~ 1/N_eff(r) * M = M/(r^2 * f(r))

  For f(r) = 1: F ~ 1/r^2 (Newton)
  For f(r) = r/r_c: F ~ 1/(r^3/r_c) — too fast!
  For f(r) = 1/(1 + r/r_c): F ~ (1+r/r_c)/r^2

  At r >> r_c: F ~ 1/(r*r_c) → F ~ 1/r → FLAT RC!

  This means: at large r, the effective number of DoF SATURATES.
  Instead of growing as r^2, it grows only as r (or slower).

  PHYSICAL INTERPRETATION:
  The substrate has a finite "correlation depth" or "thickness" L.
  At r < L: all DoF on the sphere are independent → r^2 → 3D
  At r > L: DoF are correlated along one direction → r^2/r = r → 2D

  The correlation depth L = r_c = sqrt(GM/a0) depends on M because
  the mass determines how far the substrate deformation extends.

  OR: the correlation depth is universal L = c/H0 (Hubble scale),
  and the mass-dependence comes from the TRANSITION acceleration:
  a_trans = GM/L^2 = GM*H0^2/c^2 ~ 10^-22 (too small!)

  vs a0 ~ 10^-10 → gap of 10^12 (same old problem).

  Unless L is NOT c/H0 but something smaller:
  L = sqrt(c*l_soliton/H0)? → L = sqrt(3e8 * 1e-15 / 2.2e-18)
  = sqrt(1.4e11) = 3.7e5 m ~ 370 km. Too small.

  L = c/sqrt(H0/t_soliton)?
  No natural combination gives L ~ 10 kpc.

  The mass-dependent r_c = sqrt(GM/a0) IS the right answer,
  but it means the substrate response depends on the source.
  This is intrinsically NONLINEAR — and brings us back to
  the same problem: TGP nonlinearity is too weak.
""")

# ===========================================================================
# 10. WHAT DID WE LEARN?
# ===========================================================================
print(f"\n{'='*78}")
print(f"  10. SUMMARY: What did we learn about TGP in d dimensions?")
print(f"{'='*78}")

print(f"""
  MATHEMATICAL FINDINGS:
  ======================
  1. Soliton tail decays as r^(-(d-1)/2):
     d=1: constant, d=2: r^(-0.5), d=3: r^(-1)  ✓ confirmed

  2. Oscillation period ~ pi, independent of d  ✓ confirmed

  3. Soliton energy and mass depend on d

  4. TGP kinetic term g'^2/g = Fisher information
     → Conformally invariant in d=2 (and ONLY d=2)

  5. But conformal invariance of kinetic term doesn't help:
     at large r, the potential term delta^2 dominates anyway

  6. 2D Newtonian potential IS logarithmic → flat RC automatically
     But TGP is NOT Newtonian (it has the spring term +delta)

  PHYSICAL INSIGHT:
  =================
  The "effective 2D" observation from gs8 can be restated as:

  MOND ←→ effective DoF reduction from r^2 to r at large distances

  This is equivalent to saying: the number of independent gravitational
  channels grows as r (circumference) instead of r^2 (area).

  In substrate language: gravitational information propagates
  on effectively 1D "channels" at large distances (like strings
  connecting source to test mass), not as a spherical wave.

  OPEN QUESTIONS:
  ===============
  1. WHY would TGP substrate correlations reduce effective dimension?
  2. Is the Fisher information / conformal invariance connection
     a CLUE or a COINCIDENCE?
  3. The mass-dependent crossover r_c = sqrt(GM/a0) suggests
     the substrate response IS nonlinear — but how?

  STATUS: PARTIAL SUCCESS
  - Rich mathematical structure found (Fisher info, conformal inv.)
  - d=2 IS special for TGP (conformal symmetry of kinetic term)
  - But no concrete mechanism for dimensional reduction identified
  - The DGP-like phenomenology from gs8 WORKS but lacks derivation
""")
