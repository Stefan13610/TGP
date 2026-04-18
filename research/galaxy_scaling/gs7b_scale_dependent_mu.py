#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs7b_scale_dependent_mu.py: Option E — Scale-dependent mu in TGP equation.

FROM gs7a: perturbative mu*r correction fails (10^12 too small).
The mechanism must be NON-PERTURBATIVE.

THIS SCRIPT: Instead of mu = const (perturbative), explore the possibility
that mu depends on the LOCAL field strength: mu = mu(delta).

Key idea: the TGP soliton equation g'' + g'^2/g + 2g'/r + g = 1
has "+g" which forces g→1 at infinity. But what if the effective
spring constant depends on the field value itself?

  g'' + g'^2/g + 2g'/r + mu^2(g) * (g - 1) = 0

This is still a nonlinear equation, but now mu runs with g.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp

print("="*78)
print("  OPTION E: Scale-dependent mu — non-perturbative approach")
print("="*78)

# ===========================================================================
# 1. THE IDEA: mu DEPENDS ON FIELD STRENGTH
# ===========================================================================
print(f"\n{'='*78}")
print(f"  1. The idea: mu(delta) — field-dependent mass term")
print(f"{'='*78}")

print(f"""
  Standard TGP soliton: g'' + g'^2/g + 2g'/r + g = 1
  Define delta = g - 1: delta'' + (delta')^2/(1+delta) + 2delta'/r + delta = 0

  In weak field (delta << 1): delta'' + 2delta'/r + delta = 0
  Solution: delta = A*sin(r)/r  (oscillating tail, wavelength 2*pi)

  NOW: what if the "+delta" term (spring constant mu^2 = 1) is actually
  mu^2 = mu^2(|delta|)?

  Physical motivation in TGP:
  - The substrate has a natural oscillation frequency (the "+g" term)
  - In strong deformation (delta ~ 1): substrate responds normally, mu ~ 1
  - In weak deformation (delta << 1): substrate response changes
  - Perhaps: mu^2 = 1 when |delta| > delta_c, mu^2 << 1 when |delta| << delta_c

  This is like a DIELECTRIC analogy:
  - Strong E-field: material polarizes normally (epsilon = const)
  - Weak E-field: material can't be polarized below thermal noise
  - Result: effective epsilon changes at critical field strength

  For gravity: "strong field" = normal Newton, "weak field" = modified (MOND-like)
  Transition at: delta_c ~ a0*r^2/(c^2*...) → some critical deformation
""")

# ===========================================================================
# 2. SIMPLEST MODEL: mu^2 = min(1, |delta|/delta_c)
# ===========================================================================
print(f"\n{'='*78}")
print(f"  2. Model 1: mu^2 = min(1, |delta|/delta_c)")
print(f"{'='*78}")

def solve_tgp_variable_mu(g0, delta_c, r_max=200, N=20000):
    """Solve g'' + g'^2/g + 2g'/r + mu^2(g)*(g-1) = 0
    with mu^2 = min(1, |g-1|/delta_c)"""
    def rhs(r, y):
        g, gp = y
        if r < 1e-10:
            return [gp, 0]
        delta = g - 1.0
        abs_delta = abs(delta)
        mu2 = min(1.0, abs_delta / delta_c) if delta_c > 0 else 1.0
        # g'' = -g'^2/g - 2g'/r - mu^2*(g-1)
        gpp = -gp**2/g - 2*gp/r - mu2 * delta
        return [gp, gpp]

    r_span = (1e-6, r_max)
    r_eval = np.linspace(1e-6, r_max, N)
    sol = solve_ivp(rhs, r_span, [g0, 0], t_eval=r_eval,
                    method='RK45', rtol=1e-10, atol=1e-12, max_step=0.05)
    return sol.t, sol.y[0], sol.y[1]

# Test with standard soliton first (delta_c = 0 means mu=1 always)
print(f"  Standard soliton (mu=1, delta_c=0):")
r_std, g_std, gp_std = solve_tgp_variable_mu(1.1, delta_c=0, r_max=100)
delta_std = g_std - 1.0
# Find tail envelope
peaks = []
for i in range(1, len(r_std)-1):
    if abs(delta_std[i]) > abs(delta_std[i-1]) and abs(delta_std[i]) > abs(delta_std[i+1]):
        if r_std[i] > 5:
            peaks.append((r_std[i], abs(delta_std[i])))

if len(peaks) > 1:
    r_p = np.array([p[0] for p in peaks])
    d_p = np.array([p[1] for p in peaks])
    # Fit: |delta| ~ A/r
    A_std = np.mean(d_p * r_p)
    print(f"    g0=1.1: tail envelope |delta|*r = {A_std:.4f}")
    print(f"    Oscillation: {len(peaks)} peaks in r=[5,100]")
    if len(peaks) >= 2:
        periods = np.diff(r_p)
        print(f"    Mean period: {np.mean(periods):.2f} (expected: ~pi = {np.pi:.2f})")

# Now test with variable mu
print(f"\n  Variable mu models:")
print(f"  {'delta_c':>10s} {'behavior':>40s} {'tail at r=50':>15s} {'tail*r':>10s}")
print(f"  {'-'*78}")

results_mu = {}
for delta_c in [1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.01, 0.05, 0.1, 0.5, 1.0]:
    try:
        r, g, gp = solve_tgp_variable_mu(1.1, delta_c, r_max=100)
        delta = g - 1.0
        # Check tail at r=50
        idx50 = np.argmin(np.abs(r - 50))
        d50 = delta[idx50]

        # Determine behavior
        if np.any(np.abs(delta[len(delta)//2:]) > 10):
            behavior = "DIVERGES"
        elif np.max(np.abs(delta[len(delta)//2:])) < 1e-15:
            behavior = "decays to zero (no tail)"
        else:
            # Check oscillation vs monotonic
            sign_changes = np.sum(np.diff(np.sign(delta[len(delta)//4:])) != 0)
            if sign_changes > 3:
                behavior = f"oscillating ({sign_changes} crossings)"
            else:
                behavior = f"monotonic/slow ({sign_changes} crossings)"

        print(f"  {delta_c:10.1e} {behavior:>40s} {d50:15.3e} {d50*50:10.3e}")
        results_mu[delta_c] = (r, g, gp, delta, behavior)
    except Exception as e:
        print(f"  {delta_c:10.1e} {'ERROR: '+str(e)[:40]:>40s}")

# ===========================================================================
# 3. BETTER MODEL: mu^2 TRANSITIONS SMOOTHLY
# ===========================================================================
print(f"\n{'='*78}")
print(f"  3. Model 2: mu^2 = |delta|^alpha (power-law running)")
print(f"{'='*78}")

print(f"""
  If mu^2 = |delta|^alpha:
  - At delta~1 (strong field): mu~1 (standard soliton)
  - At delta~10^-6 (galaxy scale): mu ~ 10^(-6*alpha/2)

  For the tail: delta ~ A*sin(mu*r)/r
  But mu depends on delta itself → self-consistent equation!

  For small delta: delta'' + 2delta'/r + |delta|^alpha * delta = 0

  If delta = B/r^beta (power-law decay):
  Then |delta|^alpha * delta = B^(1+alpha) / r^(beta*(1+alpha))
  And delta'' ~ beta*(beta+1)*B/r^(beta+2), 2delta'/r ~ -2*beta*B/r^(beta+2)

  Balance: beta*(beta+1) - 2*beta = beta*(beta-1) vs r^(beta*(1+alpha) - beta - 2)
  Need: beta*(1+alpha) = beta + 2  →  beta*alpha = 2  →  beta = 2/alpha

  So: delta ~ B / r^(2/alpha)

  For different alpha:
""")

print(f"  {'alpha':>6s} {'beta=2/alpha':>12s} {'delta(r)':>20s} {'v^2(r)':>20s} {'Flat RC?':>10s}")
print(f"  {'-'*72}")
for alpha in [0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]:
    beta = 2.0/alpha
    # v^2 = r * d(Phi)/dr ~ r * d(delta)/dr ~ r * B / r^(beta+1) = B / r^beta
    # For flat RC: v^2 = const → beta = 0 → alpha = inf (not physical)
    # For MOND-like: v^2 ~ const at large r means Phi ~ ln(r) → delta ~ ln(r)/r??

    # Actually in MOND: g = sqrt(g_N * a0) → g ~ sqrt(GM*a0)/r
    # → Phi ~ sqrt(GM*a0) * ln(r)
    # → delta = Phi/c^2 ~ sqrt(GM*a0)*ln(r)/c^2
    # → |delta| decays slower than any power law!

    v2_scaling = f"v^2 ~ r^(-{beta:.1f}+1) = r^{1-beta:.1f}"
    flat = "YES!" if abs(beta - 1.0) < 0.1 else ("close" if abs(beta-1.0) < 0.3 else "no")
    print(f"  {alpha:6.2f} {beta:12.2f} {'delta ~ r^(-'+f'{beta:.1f}'+')':>20s} {v2_scaling:>20s} {flat:>10s}")

print(f"""
  For FLAT rotation curve: need v^2 = const
  → need Phi ~ ln(r) → need delta ~ ln(r)/r (NOT a power law!)

  With alpha = 2: beta = 1, v^2 ~ const → FLAT! (v^2 = B = const)

  Wait — let's check: if delta ~ B/r, then:
  v^2(r) = r * |dPhi/dr| = r * c^2 * |delta'/r + ...|

  Actually v^2 = G*M_enc(r)/r. The extra mass from TGP:
  M_extra(r) ~ integral rho_extra * 4pi*r^2 dr
  where rho_extra comes from the modified Poisson equation.

  Let me solve this properly numerically.
""")

# ===========================================================================
# 4. NUMERICAL SOLUTION: mu^2 = |delta|^alpha
# ===========================================================================
print(f"\n{'='*78}")
print(f"  4. Numerical solution of delta'' + 2delta'/r + |delta|^alpha * delta = -S(r)")
print(f"{'='*78}")

def solve_modified_poisson(alpha, M_galaxy_code=1.0, a_scale=5.0, r_max=200, N=10000):
    """
    Solve the modified Poisson equation in TGP with mu^2 = |delta|^alpha.

    Equation: delta'' + 2delta'/r + |delta|^alpha * delta = -rho(r)

    where rho is Hernquist profile: rho = M*a/(2pi*r*(r+a)^3)

    Using delta = g-1, with boundary conditions:
    delta(0) finite, delta(inf) → 0

    We solve as a boundary value problem by shooting from large r inward,
    or as an initial value problem from r=0 outward with proper BC.

    Actually, for comparison with Newton, let's solve:
    - Newton: delta_N'' + 2delta_N'/r = -rho  → delta_N = -M/(4pi*r) * ...
    - TGP:    delta_T'' + 2delta_T'/r + |delta_T|^alpha * delta_T = -rho

    The extra term |delta|^alpha * delta acts as a "source sink" —
    it removes energy from the radial mode and stores it in the spring.
    """
    # Hernquist source
    def rho(r):
        return M_galaxy_code * a_scale / (2*np.pi * r * (r + a_scale)**3) if r > 1e-10 else 0

    # Enclosed mass (Hernquist)
    def M_enc(r):
        return M_galaxy_code * r**2 / (r + a_scale)**2

    # Newtonian potential: delta_N(r) = -M_enc(r) / r (in code units where 4piG=1)

    # Solve TGP equation from r=0 outward
    # Near r=0: delta ~ delta_0 + delta_2 * r^2 + ...
    # delta'' + 2delta'/r + |delta|^alpha * delta = -rho
    # At r=0: 3*delta_2 + |delta_0|^alpha * delta_0 = -rho(0)

    def rhs(r, y):
        delta, delta_p = y
        if r < 1e-6:
            return [delta_p, 0]

        # Source
        source = rho(r)

        # mu^2 term
        abs_d = abs(delta)
        if abs_d > 1e-30:
            mu2_delta = abs_d**alpha * delta  # mu^2(delta) * delta
        else:
            mu2_delta = 0

        # delta'' = -2delta'/r - mu^2*delta - rho
        deltapp = -2*delta_p/r - mu2_delta - source
        return [delta_p, deltapp]

    # Initial conditions: need to find delta(0) that gives delta→0 at infinity
    # This is a shooting problem. For now, try delta(0) = -M/(4pi*a) (Newtonian center)
    delta0_newton = -M_galaxy_code / (4*np.pi * a_scale)

    # Try several initial values
    r_eval = np.linspace(1e-4, r_max, N)

    best_sol = None
    best_delta0 = None
    best_residual = 1e30

    for factor in np.linspace(0.5, 2.0, 30):
        d0 = delta0_newton * factor
        try:
            sol = solve_ivp(rhs, (1e-4, r_max), [d0, 0], t_eval=r_eval,
                           method='RK45', rtol=1e-8, atol=1e-10, max_step=0.5)
            if sol.success:
                # Check boundary: delta should → 0 at large r
                residual = abs(sol.y[0][-1])
                if residual < best_residual and not np.any(np.isnan(sol.y[0])):
                    best_residual = residual
                    best_sol = sol
                    best_delta0 = d0
        except:
            continue

    if best_sol is not None:
        return best_sol.t, best_sol.y[0], best_sol.y[1], best_delta0
    else:
        return None, None, None, None

print(f"  Solving for different alpha values...")
print(f"  (Hernquist galaxy, M=1, a=5, code units)")
print()

# Also compute Newtonian reference
r_ref = np.linspace(0.5, 200, 5000)
a_scale = 5.0
M_enc_ref = r_ref**2 / (r_ref + a_scale)**2  # Hernquist enclosed mass (M=1)
v2_newton = M_enc_ref / r_ref  # v^2 = GM_enc/r (with G=1/4pi absorbed)

print(f"  {'alpha':>6s} {'delta(0)':>10s} {'v2(r=20)':>10s} {'v2(r=50)':>10s} {'v2(r=100)':>10s} {'v2_ratio(50/20)':>15s} {'Flat?':>8s}")
print(f"  {'-'*72}")

alpha_results = {}
for alpha in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]:
    r, delta, deltap, d0 = solve_modified_poisson(alpha, r_max=200)
    if r is not None:
        # Compute v^2 from the potential
        # v^2 = r * |dPhi/dr| = r * |delta'| (approximately)
        v2 = r * np.abs(deltap)

        idx20 = np.argmin(np.abs(r - 20))
        idx50 = np.argmin(np.abs(r - 50))
        idx100 = np.argmin(np.abs(r - 100))

        v2_20 = v2[idx20]
        v2_50 = v2[idx50]
        v2_100 = v2[idx100]

        ratio = v2_50 / v2_20 if v2_20 > 0 else 0
        flat = "YES!" if 0.8 < ratio < 1.2 else ("~flat" if 0.6 < ratio < 1.4 else "no")

        print(f"  {alpha:6.1f} {d0:10.4f} {v2_20:10.4f} {v2_50:10.4f} {v2_100:10.4f} {ratio:15.3f} {flat:>8s}")
        alpha_results[alpha] = (r, delta, deltap, v2)
    else:
        print(f"  {alpha:6.1f} {'FAILED':>10s}")

# ===========================================================================
# 5. COMPARISON: v^2 PROFILES
# ===========================================================================
print(f"\n{'='*78}")
print(f"  5. Detailed v^2(r) profiles for promising alpha values")
print(f"{'='*78}")

for alpha in [0.0, 1.0, 2.0]:
    if alpha in alpha_results:
        r, delta, deltap, v2 = alpha_results[alpha]
        print(f"\n  alpha = {alpha:.1f}:")
        print(f"    {'r':>6s} {'delta':>12s} {'v^2':>10s} {'v^2/v^2_max':>12s}")
        print(f"    {'-'*44}")
        v2_max = np.max(v2[:len(v2)//2]) if np.max(v2[:len(v2)//2]) > 0 else 1
        for ri in [2, 5, 10, 15, 20, 30, 50, 70, 100, 150]:
            idx = np.argmin(np.abs(r - ri))
            print(f"    {ri:6d} {delta[idx]:12.4e} {v2[idx]:10.4e} {v2[idx]/v2_max:12.3f}")

# ===========================================================================
# 6. ALTERNATIVE: LOGARITHMIC POTENTIAL
# ===========================================================================
print(f"\n{'='*78}")
print(f"  6. What potential gives flat rotation curves?")
print(f"{'='*78}")

print(f"""
  Working backwards from observation:

  If v^2 = v_flat^2 = const at large r, then:
    F = v^2/r → F ~ 1/r (not 1/r^2!)
    Phi(r) = v_flat^2 * ln(r) + const

  In MOND (deep regime): v^4 = G*M*a0
    v_flat = (G*M*a0)^(1/4)
    Phi(r) = sqrt(G*M*a0) * ln(r)

  In terms of TGP delta = Phi/c^2:
    delta(r) = sqrt(G*M*a0) * ln(r) / c^2

  For MW (M=7e10 M_sun, a0=1.2e-10):
""")

G_val = 6.674e-11
c_val = 2.998e8
M_MW = 7e10 * 1.989e30
a0_val = 1.2e-10

v_flat = (G_val * M_MW * a0_val)**0.25
delta_10kpc = np.sqrt(G_val * M_MW * a0_val) * np.log(10 * 3.086e19) / c_val**2
delta_50kpc = np.sqrt(G_val * M_MW * a0_val) * np.log(50 * 3.086e19) / c_val**2

print(f"    v_flat = (G*M*a0)^(1/4) = {v_flat:.0f} m/s = {v_flat/1e3:.0f} km/s")
print(f"    delta(10 kpc) = {delta_10kpc:.3e}")
print(f"    delta(50 kpc) = {delta_50kpc:.3e}")

print(f"""
  So delta ~ 10^(-6) at galactic scales — VERY weak field.

  For the TGP equation with mu^2 = |delta|^alpha:
  mu^2 ~ (10^-6)^alpha

  Need: mu^2 * delta ~ delta''/r ~ delta/r^2 (balance of terms)
  (10^-6)^alpha * 10^-6 ~ 10^-6 / r^2
  (10^-6)^(alpha+1) ~ 10^-6 / r^2
  r^2 ~ 10^(-6) / 10^(-6*(alpha+1)) = 10^(6*alpha)
  r ~ 10^(3*alpha) (in natural units)

  For transition at r ~ 10 kpc and soliton radius ~ 10^-15 m:
  r_galaxy / l_soliton ~ 3e22 / 10^-15 ~ 3e37

  Need: 10^(3*alpha) ~ 3e37 → alpha ~ 37/3 ~ 12

  This is a VERY steep power law — essentially a step function.

  INTERPRETATION:
  mu^2 = |delta|^12 means:
  - At delta = 1 (soliton core): mu = 1 (normal oscillation)
  - At delta = 10^-6 (galaxy): mu = 10^-36 (effectively zero!)
  - The "+g" term simply TURNS OFF at weak fields

  If "+g" turns off, what's left?
  g'' + g'^2/g + 2g'/r = 0  (no spring term)
  delta'' + 2delta'/r ~ 0  (weak field)
  Solution: delta = A/r + B  (pure 1/r, Newton!)

  → Turning off the spring just gives back Newton — NOT MOND!
""")

# ===========================================================================
# 7. THE REAL QUESTION: WHAT EQUATION GIVES ln(r) POTENTIAL?
# ===========================================================================
print(f"\n{'='*78}")
print(f"  7. What modification gives Phi ~ ln(r) at large r?")
print(f"{'='*78}")

print(f"""
  Standard Poisson: nabla^2 Phi = 4*pi*G*rho
  Solution: Phi ~ -M/r (3D)

  For Phi ~ ln(r): need nabla^2 ln(r) = ?
  In 3D spherical: (1/r^2) d/dr (r^2 d(ln r)/dr) = (1/r^2) d/dr (r) = 1/r^2

  So: nabla^2 Phi = v_flat^2 / r^2

  This means: there's an EXTRA SOURCE proportional to 1/r^2!
  The "phantom dark matter" density: rho_phantom = v_flat^2 / (4*pi*G*r^2)
  = a0 * M_bar / (4*pi*G*r^2)^(1/2) ... no, let me redo:

  rho_phantom(r) = v_flat^2 / (4*pi*G*r^2)
  = sqrt(G*M*a0) / (4*pi*G*r^2)
  = sqrt(M*a0/G) / (4*pi*r^2)

  This is EXACTLY the isothermal DM halo! (rho ~ 1/r^2)
  And it's what gs2 already computed!

  HOW could TGP produce this extra source?

  In TGP: the field equation is nabla^2 delta + (nabla delta)^2/(1+delta) + delta = source

  Rearrange: nabla^2 delta = source - delta - (nabla delta)^2/(1+delta)

  Compare with Newton: nabla^2 delta_N = source

  Extra terms: -delta - (nabla delta)^2/(1+delta)

  At galactic scale (delta ~ 10^-6):
  - "-delta" ~ 10^-6 (negligible vs source ~ delta/r ~ 10^-6/r)
  - "-(nabla delta)^2" ~ delta^2/r^2 ~ 10^-12/r^2 (NEGLIGIBLE)

  → Standard TGP gives corrections of order delta ~ 10^-6.
  → Need corrections of order 1 (100% of Newton).
  → Gap: 10^6 (same old problem!)

  UNLESS: there's a term we're missing that's NOT suppressed by delta.
""")

# ===========================================================================
# 8. KEY REALIZATION: COSMOLOGICAL BACKGROUND
# ===========================================================================
print(f"\n{'='*78}")
print(f"  8. KEY REALIZATION: the cosmological background field")
print(f"{'='*78}")

print(f"""
  We've been treating delta as the TOTAL deviation from g=1.
  But in an expanding universe, the BACKGROUND g is NOT 1!

  In TGP cosmology, the Friedmann-like equation gives:
  g_background = g_cosmo(t) where g_cosmo evolves with expansion.

  The PERTURBATION on top of cosmological background:
  g = g_cosmo + delta_local

  The effective equation for delta_local:
  delta_local'' + ... + mu_eff^2 * delta_local = source

  where mu_eff^2 depends on g_cosmo and its derivatives!

  Specifically, expanding around g_cosmo:
  g'' + g'^2/g + 2g'/r + g = 1

  Let g = g_cosmo(1 + phi), where phi is the local perturbation:
  g_cosmo*phi'' + ... + g_cosmo*(1 + nonlinear terms)*phi = -source

  The key: g_cosmo is NOT static — it evolves on timescale 1/H0.
  For a static galaxy, we need the TIME-DEPENDENT equation:

  (1/c^2) * d^2g/dt^2 - nabla^2 g - g'^2/g + g = 1 + source

  (The spatial TGP equation is the static limit of a wave equation!)

  In an expanding background:
  g_cosmo(t) evolves, and the perturbation phi satisfies an equation
  where H0 appears naturally as a damping/frequency term.

  THIS IS WHERE a0 = c*H0/(2pi) comes from!
  Not from spatial structure, but from the TIME evolution of the background.
""")

# ===========================================================================
# 9. THE FREQUENCY ARGUMENT
# ===========================================================================
print(f"\n{'='*78}")
print(f"  9. The frequency argument: static limit of wave equation")
print(f"{'='*78}")

print(f"""
  Full TGP wave equation (1+1D for simplicity):
  (1/c^2) g_tt - g_rr - g_r^2/g - 2g_r/r + g = 1 + S

  In expanding universe, use conformal time eta:
  g_tt → a^(-2) (g_eta,eta + (a'/a) g_eta)
  where a is scale factor, a'/a = a*H (conformal Hubble)

  For static local structure (galaxy), g_eta = 0 in comoving frame.
  But NOT in the substrate frame!

  In the substrate frame, a galaxy at comoving distance r has:
  - Static internal structure
  - But the substrate around it is EXPANDING

  The boundary condition at large r:
  g → g_cosmo(t) (expanding background), NOT g → 1 (static)

  This changes the equation from:
  g'' + g'^2/g + 2g'/r + g = 1    (static, g→1 at infinity)
  to:
  g'' + g'^2/g + 2g'/r + g = g_cosmo(t)   (quasi-static, g→g_cosmo)

  If g_cosmo(t) = 1 + epsilon(t) with epsilon ~ H0*t:
  delta'' + 2delta'/r + delta = epsilon(t)

  The solution has a HOMOGENEOUS part (oscillating as before)
  and a PARTICULAR solution: delta_part = epsilon(t) = const at fixed t.

  This doesn't help — it's just a constant offset.

  BUT: if we consider the TEMPORAL evolution properly:

  In Fourier space: delta(r, omega) satisfies:
  delta_rr + 2delta_r/r + (1 - omega^2/c^2)*delta = source

  For omega = 0 (truly static): standard soliton equation
  For omega = H0: 1 - H0^2/c^2 ~ 1 (negligible change)

  → The frequency H0 is too small to change the spatial equation.
  → Same 10^-36 problem.

  CONCLUSION FOR OPTION E:
  Scale-dependent mu doesn't help PERTURBATIVELY.
  Turning off the spring term just gives Newton back.
  The equation structure doesn't naturally produce ln(r) potential.
""")

# ===========================================================================
# 10. ASSESSMENT
# ===========================================================================
print(f"\n{'='*78}")
print(f"  10. ASSESSMENT of Option E")
print(f"{'='*78}")

print(f"""
  TESTED:
  1. mu^2 = min(1, |delta|/delta_c)  → modifies soliton tail, but
     at galactic scales delta ~ 10^-6 so mu → 0 → Newton (not MOND)

  2. mu^2 = |delta|^alpha → need alpha ~ 12 for transition at galaxy scale
     → but turning off spring gives Newton, not MOND

  3. Cosmological background correction:
     g → g_cosmo ≠ 1 changes BC but not force law

  4. Wave equation with omega = H0:
     (1 - H0^2/c^2) ~ 1 → no significant change

  FUNDAMENTAL PROBLEM:
  The TGP equation delta'' + 2delta'/r + delta = source
  has TWO terms that could modify Newton:
  a) The "+delta" spring term → gives oscillating tail sin(r)/r
  b) The nonlinear (delta')^2/(1+delta) → suppressed by delta ~ 10^-6

  Neither can produce a 1/r force (needed for flat RC) because:
  - The spring term adds oscillations, not 1/r modification
  - The nonlinear term is too weak in the weak-field regime

  To get MOND-like behavior, we need a term that:
  - Is NOT proportional to delta (would be too small)
  - IS proportional to sqrt(delta * something) or delta'/r
  - DOESN'T vanish in the weak-field limit

  STATUS: FAILED
  - mu(delta) idea: doesn't produce flat RC
  - Cosmological background: negligible correction
  - No path to MOND-like transition found
""")
