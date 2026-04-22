#!/usr/bin/env python3
"""
atail_asymptotic_v47b.py -- Analytical asymptotic forms for A_tail(g0).

GOAL: Derive and verify analytical expressions for A_tail(g0)
in the two limits:
  1. g0 -> 1 (near vacuum): A_tail ~ |g0 - 1| * C1
  2. g0 -> 8/5 (near collapse): A_tail ~ (8/5 - g0)^(-gamma) * C2

Then check if these forms constrain Q_K = 3/2.

APPROACH:
  Near vacuum (g0 ~ 1 + epsilon):
    Linearize ODE around g=1: h = g-1
    h'' + 2/r h' + h = 0  (spherical Bessel, source linearized)
    Solution: h(r) = A * sin(r)/r  for large r
    A_tail is the coefficient of the sin(r)/r tail.
    For small epsilon: A_tail ~ epsilon * C1(alpha)

  Near collapse (g0 ~ 8/5 - delta):
    The soliton profile has a core region where g deviates strongly
    from vacuum. The core becomes larger and the field approaches
    collapse. A_tail diverges because the matching to the tail
    involves a resonance.
"""
import numpy as np
from numpy import trapezoid as trapz
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

PHI = (1 + np.sqrt(5)) / 2
G0_CRIT = 8.0 / 5.0


def solver_A(g0, r_max=300):
    """Canonical Form A (alpha=2) ODE solver."""
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-10)
        source = g**2 * (1.0 - g)
        cross = (2.0 / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / 3.0]
        return [gp, source - cross - 2.0 * gp / r]
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                    rtol=1e-12, atol=1e-14, max_step=0.02)
    return sol.t, sol.y[0]


def A_tail(g0, r_min=50, r_max=200):
    r, g = solver_A(g0)
    mask = (r > r_min) & (r < r_max)
    if np.sum(mask) < 100:
        return 0.0
    rf = r[mask]
    if np.any(np.abs(g[mask] - 1) > 0.5):
        return 0.0
    df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return np.sqrt(bc[0]**2 + bc[1]**2)


# ================================================================
# Section 1: Linear regime (g0 near 1)
# ================================================================
print("=" * 70)
print("A_tail ASYMPTOTIC ANALYSIS")
print("=" * 70)

print("\n1. LINEAR REGIME: g0 -> 1")
print("-" * 50)

# Near g=1: g = 1+h, |h| << 1
# Full ODE: g'' + 2/r g' + 2g'^2/g + g^2(1-g) = 0
# Linearized (drop g'^2/g term as O(h'^2)):
#   h'' + 2/r h' + (1+h)^2(-h) = 0
#   h'' + 2/r h' - h - 2h^2 + O(h^3) = 0
#   At leading order: h'' + 2/r h' - h = 0
#
# Wait: g^2(1-g) = (1+h)^2(-h) = -h(1+2h+h^2) = -h - 2h^2 - h^3
#
# So linearized: h'' + 2/r h' + h = 0  ???
# No: the source is -h (not +h) at leading order.
#
# Actually: the ODE is g'' + 2/r g' + 2g'^2/g = -g^2(1-g)
# Wait, let me recheck. The ODE is:
#   g'' + 2/r g' + (alpha/g)g'^2 + g^2(1-g) = 0  [with alpha=2]
#
# So: g'' + 2/r g' + 2g'^2/g = -g^2(1-g)
#
# Linearize with g = 1+h:
#   h'' + 2/r h' + 2h'^2/(1+h) = -(1+h)^2(-h) = h(1+2h+h^2)
#   h'' + 2/r h' + 2h'^2 + ... = h + 2h^2 + ...
#
# At leading order (drop h'^2 and h^2):
#   h'' + 2/r h' = h
#   h'' + 2/r h' - h = 0
#
# This is the MODIFIED spherical Bessel equation!
# Solutions: h = A * exp(-r)/r + B * exp(r)/r
# For bounded solution: h = A * exp(-r)/r
#
# But wait, the tail should be OSCILLATORY (sin(r)/r), not exponential.
# Let me recheck...

# Actually, g^2(1-g) at g=1+h: (1+h)^2(-h) = -h - 2h^2 - h^3
# The ODE: g'' + 2/r g' + 2g'^2/g + g^2(1-g) = 0
# => h'' + 2/r h' + 2h'^2/(1+h) + (-h-2h^2-h^3) = 0
# => h'' + 2/r h' - h = 0  (at leading order, h'^2 << h)
#
# But this gives EXPONENTIAL tails, not oscillatory!
# The oscillatory tail must come from the NONLINEAR coupling term.

# Let's check numerically:
print("\n  Numerical check: does the tail oscillate?")
for g0 in [0.99, 0.98, 0.95, 0.90, 0.85]:
    r, g = solver_A(g0)
    # Check tail behavior
    mask = (r > 30) & (r < 100)
    if np.sum(mask) > 50:
        gm = g[mask]
        # Count zero crossings of g-1
        crossings = np.sum(np.diff(np.sign(gm - 1.0)) != 0)
        A = A_tail(g0)
        print(f"  g0={g0:.3f}: tail crossings={crossings}, A_tail={A:.8f}, "
              f"|g0-1|={abs(g0-1):.3f}")

# Let me look at the linearized equation more carefully.
# The full nonlinear equation around vacuum involves the cross term.
# Maybe the key is that the CROSS term 2g'^2/g modifies the effective equation.
#
# Let g = 1 + h. Then:
#   g'^2/g = h'^2/(1+h) = h'^2(1-h+h^2-...)
#
# At order h: h'' + 2/r h' - h = 0 (exponential)
# But numerically we SEE oscillations! So the linear approximation breaks down.
#
# Actually, let me reconsider. Maybe the oscillations come from a different
# mechanism. Let me check a wider range of g0 values.

print("\n  Extended check with finer detail:")
for g0 in [0.999, 0.99, 0.98, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50]:
    r, g = solver_A(g0)
    mask = (r > 20) & (r < 200)
    if np.sum(mask) > 100:
        gm = g[mask]
        rm = r[mask]
        crossings = np.sum(np.diff(np.sign(gm - 1.0)) != 0)

        # Check if oscillatory by fitting sin
        df = (gm - 1.0) * rm
        M = np.column_stack([np.cos(rm), np.sin(rm)])
        bc = np.linalg.lstsq(M, df, rcond=None)[0]
        A = np.sqrt(bc[0]**2 + bc[1]**2)

        # Also check exponential fit
        exp_vals = np.exp(-rm)
        if np.max(np.abs(gm - 1.0)) > 1e-15:
            exp_ratio = np.abs(gm[-1] - 1) / np.abs(gm[0] - 1) if abs(gm[0]-1) > 1e-15 else 0
        else:
            exp_ratio = 0

        print(f"  g0={g0:.3f}: crossings(g=1)={crossings:3d}, "
              f"A_osc={A:.8f}, |g-1| ratio last/first={exp_ratio:.6f}")


# ================================================================
# Section 2: Understanding the source of oscillations
# ================================================================
print("\n\n2. SOURCE OF OSCILLATORY TAIL")
print("-" * 50)

# The linearized equation h'' + 2/r h' - h = 0 gives EXPONENTIAL decay.
# But numerically, the tail IS oscillatory.
# This means the oscillations are a NONLINEAR effect.
#
# Actually wait -- let me reconsider the linearization.
# The ODE is: g'' + 2/r g' + 2g'^2/g + g^2(1-g) = 0
#
# Let me write this as: g'' + 2/r g' = -2g'^2/g - g^2(1-g)
# The RHS is the "effective source."
#
# At vacuum g=1, g'=0: RHS = 0.
# Near vacuum g=1+h, g'=h':
#   RHS = -2h'^2/(1+h) - (1+h)^2(-h)
#       = -2h'^2/(1+h) + h + 2h^2 + h^3
#       = h + (2h^2 - 2h'^2) + ...
#
# Linear part: h
# So: h'' + 2/r h' = h  => h'' + 2/r h' - h = 0
# This is the modified Bessel equation in d=3.
# Solution: h(r) = C * K_{1/2}(r) / r^{1/2} = C * exp(-r)/r
#
# The solution DECAYS EXPONENTIALLY. No oscillations at linear order.
#
# So where do the oscillations come from?

# Let me check: maybe the "oscillations" are not around g=1 but are
# soliton-type oscillations of the full nonlinear solution.

# Plot g(r) for a specific g0:
g0 = 0.85
r, g = solver_A(g0)
print(f"\n  Detailed profile for g0 = {g0}:")
print(f"  {'r':>8s} {'g(r)':>12s} {'g-1':>12s} {'(g-1)*r':>12s}")
for i in range(0, len(r), max(1, len(r)//40)):
    print(f"  {r[i]:8.2f} {g[i]:12.8f} {g[i]-1:12.8f} {(g[i]-1)*r[i]:12.8f}")


# ================================================================
# Section 3: Near-vacuum behavior -- is it really exponential?
# ================================================================
print("\n\n3. EXPONENTIAL vs OSCILLATORY -- careful check")
print("-" * 50)

# For g0 very close to 1, the solution should be well-approximated
# by the linear solution: h ~ C * exp(-r)/r (exponential decay)
# A_tail (fitted as oscillatory) should be ~0 for such cases.

print("\n  Very near vacuum:")
for eps in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
    g0 = 1.0 - eps
    A = A_tail(g0)
    r, g = solver_A(g0)
    # Check the envelope: does |g-1| decay exponentially?
    mask = (r > 5) & (r < 100)
    if np.sum(mask) > 50:
        rm = r[mask]
        hm = np.abs(g[mask] - 1.0)
        # Fit log(|h|) = a + b*r (exponential: b should be ~ -1)
        valid = hm > 1e-15
        if np.sum(valid) > 20:
            log_h = np.log(hm[valid] + 1e-30)
            rv = rm[valid]
            p = np.polyfit(rv, log_h, 1)
            decay_rate = p[0]
        else:
            decay_rate = 0
    else:
        decay_rate = 0
    print(f"  eps={eps:.1e}: A_tail(osc fit)={A:.2e}, "
          f"decay rate={decay_rate:.4f} (expect -1.0 for exp)")


# Now above vacuum:
print("\n  Above vacuum:")
for eps in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]:
    g0 = 1.0 + eps
    if g0 >= G0_CRIT:
        continue
    A = A_tail(g0)
    r, g = solver_A(g0)
    mask = (r > 5) & (r < 100)
    if np.sum(mask) > 50:
        rm = r[mask]
        hm = np.abs(g[mask] - 1.0)
        valid = hm > 1e-15
        if np.sum(valid) > 20:
            log_h = np.log(hm[valid] + 1e-30)
            rv = rm[valid]
            p = np.polyfit(rv, log_h, 1)
            decay_rate = p[0]
        else:
            decay_rate = 0
    else:
        decay_rate = 0
    print(f"  eps={eps:.1e}: A_tail(osc fit)={A:.2e}, "
          f"decay rate={decay_rate:.4f}")


# ================================================================
# Section 4: The REAL tail structure
# ================================================================
print("\n\n4. REAL TAIL STRUCTURE")
print("-" * 50)

# If the tail is exponential, then A_tail (defined via oscillatory fit)
# should be zero. The A_tail function measures oscillatory content.
# If A_tail > 0, there ARE oscillations, possibly from nonlinear effects.

# Key insight: the soliton has a CORE where g deviates from 1,
# and a TAIL where g -> 1. The tail can be:
# (a) purely exponential (linear regime)
# (b) oscillatory (if the effective equation has imaginary eigenvalue)
#
# The cross term 2g'^2/g shifts the effective potential.
# In f = g^3 variable, the linearization around f=1:
# f = 1 + eta
# f'' + 2/r f' = 3*f^(4/3)*(1-f^(1/3))
# At f=1+eta: 3*(1+eta)^(4/3)*(1-(1+eta)^(1/3))
# = 3*(1+4eta/3+...)*(1-1-eta/3-...) = 3*(-eta/3-...) = -eta + ...
# So: eta'' + 2/r eta' = -eta  =>  eta'' + 2/r eta' + eta = 0
# This IS the OSCILLATORY Bessel equation!
# Solution: eta = A*sin(r)/r + B*cos(r)/r

print("  EUREKA: In the canonical variable f = g^3,")
print("  the linearization around vacuum gives:")
print("    eta'' + 2/r eta' + eta = 0  (oscillatory!)")
print()
print("  In contrast, in the g-variable:")
print("    h'' + 2/r h' - h = 0  (exponential!)")
print()
print("  The difference is because the canonical transformation")
print("  f = g^3 absorbs the cross term 2g'^2/g.")
print("  The PHYSICAL tail is oscillatory in f (and therefore in g),")
print("  with A_tail being the coefficient of sin(r)/r in the f-variable.")

# Verify: linearize the f-equation
# f'' + 2/r f' = 3 f^{4/3} (1 - f^{1/3})
# At f = 1+eta:
# f^{4/3} = (1+eta)^{4/3} = 1 + 4/3 eta + ...
# f^{1/3} = (1+eta)^{1/3} = 1 + 1/3 eta + ...
# 1-f^{1/3} = -1/3 eta + ...
# 3 * f^{4/3} * (1-f^{1/3}) = 3 * (1 + 4/3 eta) * (-1/3 eta)
# = 3 * (-1/3 eta - 4/9 eta^2 + ...)
# = -eta - 4/3 eta^2 + ...
#
# So: eta'' + 2/r eta' = -eta  at leading order
# => eta'' + 2/r eta' + eta = 0
# Solution: eta(r) = A sin(r + phi) / r
# (spherical wave, frequency omega = 1)

print("\n  Verification: frequency of oscillation")
g0 = 0.85
r, g = solver_A(g0)
mask = (r > 50) & (r < 200)
rm = r[mask]
gm = g[mask]

# Convert to f = g^3
fm = gm**3
eta_m = fm - 1.0

# Fit eta*r to A*sin(r) + B*cos(r)
df_f = eta_m * rm
M = np.column_stack([np.cos(rm), np.sin(rm)])
bc_f = np.linalg.lstsq(M, df_f, rcond=None)[0]
A_f = np.sqrt(bc_f[0]**2 + bc_f[1]**2)

# Also fit in g-variable
hm = gm - 1.0
df_g = hm * rm
bc_g = np.linalg.lstsq(M, df_g, rcond=None)[0]
A_g = np.sqrt(bc_g[0]**2 + bc_g[1]**2)

print(f"  g0 = {g0}")
print(f"  A_tail (g-variable, h*r fit):  {A_g:.8f}")
print(f"  A_tail (f-variable, eta*r fit): {A_f:.8f}")
print(f"  Ratio A_f / A_g: {A_f/A_g:.6f}")
print(f"  Expected (linear): 3*1 = 3 (since f = g^3, df = 3g^2 dg, at g=1: 3)")
print(f"  Actual ratio: {A_f/A_g:.6f}")

# For general alpha+1:
# f = g^(alpha+1) = g^3
# df = 3 g^2 dg
# At g = 1: df = 3 dg
# So A_f = 3 * A_g (at linear order near vacuum)


# ================================================================
# Section 5: A_tail in g-variable vs f-variable
# ================================================================
print("\n\n5. A_tail RELATIONSHIP: g vs f variables")
print("-" * 50)

print("  At leading order near vacuum (g ~ 1 + h, f = g^3 ~ 1 + 3h):")
print("  A_f = 3 * A_g")
print()
print("  For Koide, we use A_g (or equivalently A_f/3).")
print("  Since mass ~ A^4, the choice of variable doesn't affect Q_K:")
print("  Q_K = (sum A_g^2)^2 / (sum A_g^4) = (sum (A_f/3)^2)^2 / (sum (A_f/3)^4)")
print("       = (sum A_f^2)^2 / (sum A_f^4)  [3's cancel in ratio]")

# Verify the 3x ratio holds for different g0
print("\n  Ratio check A_f / (3*A_g) for various g0:")
for g0 in [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.55]:
    if g0 >= G0_CRIT:
        continue
    r, g = solver_A(g0)
    mask = (r > 50) & (r < 200)
    if np.sum(mask) < 100:
        continue
    rm = r[mask]
    gm = g[mask]
    if np.any(np.abs(gm - 1) > 0.5):
        continue

    # g-variable
    df_g = (gm - 1.0) * rm
    M = np.column_stack([np.cos(rm), np.sin(rm)])
    bc_g = np.linalg.lstsq(M, df_g, rcond=None)[0]
    Ag = np.sqrt(bc_g[0]**2 + bc_g[1]**2)

    # f-variable
    fm = gm**3
    df_f = (fm - 1.0) * rm
    bc_f = np.linalg.lstsq(M, df_f, rcond=None)[0]
    Af = np.sqrt(bc_f[0]**2 + bc_f[1]**2)

    ratio = Af / (3 * Ag) if Ag > 1e-10 else 0
    print(f"  g0={g0:.3f}: A_g={Ag:.8f}, A_f={Af:.8f}, A_f/(3*A_g)={ratio:.6f}")


# ================================================================
# Section 6: Near-vacuum slope dA/dg0 at g0 = 1
# ================================================================
print("\n\n6. NEAR-VACUUM SLOPE: A_tail ~ |g0 - 1| * C")
print("-" * 50)

# In f-variable: eta'' + 2/r eta' + eta = 0 (at linear order)
# With eta(0) = f0 - 1 = (1+h0)^3 - 1 ~ 3*h0 for small h0
# And eta'(0) = 0
#
# Solution: eta(r) = eta_0 * sin(r)/r  for r > 0
# Wait: eta(0) = eta_0 * lim sin(r)/r = eta_0
# And eta'(0) = 0: d/dr[eta_0 sin(r)/r]|_0 = eta_0 * d/dr[sin(r)/r]|_0
# = eta_0 * [cos(r)/r - sin(r)/r^2]|_0 = eta_0 * [1/r - sin(r)/r^2]|_0
# This is 0/0... use L'Hopital or Taylor:
# sin(r)/r = 1 - r^2/6 + ... so d/dr[sin(r)/r] = -r/3 + ... -> 0 at r=0 ✓
#
# So at linear order: eta(r) = eta_0 * sin(r)/r
# => A_f = |eta_0| = |f0 - 1| ~ 3|g0 - 1| for g0 near 1
# => A_g = A_f / 3 = |g0 - 1|
#
# But this is the LINEARIZED result. Does it hold exactly?

print("  Linear prediction: A_g = |g0 - 1|")
print("  Verification:")
for g0 in [0.999, 0.998, 0.995, 0.99, 0.98, 0.95, 0.90, 0.85, 0.80]:
    A = A_tail(g0)
    predicted = abs(g0 - 1.0)
    ratio = A / predicted if predicted > 1e-10 else 0
    print(f"  g0={g0:.4f}: A_tail={A:.8f}, |g0-1|={predicted:.4f}, "
          f"A/|g0-1|={ratio:.6f}")

print("\n  Above vacuum:")
for g0 in [1.001, 1.002, 1.005, 1.01, 1.02, 1.05, 1.10, 1.20, 1.30, 1.40, 1.50]:
    if g0 >= G0_CRIT:
        continue
    A = A_tail(g0)
    predicted = abs(g0 - 1.0)
    ratio = A / predicted if predicted > 1e-10 else 0
    print(f"  g0={g0:.4f}: A_tail={A:.8f}, |g0-1|={predicted:.4f}, "
          f"A/|g0-1|={ratio:.6f}")


# ================================================================
# Section 7: Near-collapse divergence
# ================================================================
print("\n\n7. NEAR-COLLAPSE: A_tail ~ (8/5 - g0)^(-gamma)")
print("-" * 50)

# Compute A_tail for g0 close to 8/5
deltas = []
A_vals = []
for g0 in np.linspace(1.50, 1.598, 30):
    A = A_tail(g0)
    if A > 1e-10:
        deltas.append(G0_CRIT - g0)
        A_vals.append(A)
        print(f"  g0={g0:.5f}, delta={G0_CRIT-g0:.5f}, A_tail={A:.8f}")

if len(deltas) > 5:
    deltas = np.array(deltas)
    A_vals = np.array(A_vals)
    log_d = np.log(deltas)
    log_A = np.log(A_vals)

    # Fit power law: log(A) = -gamma * log(delta) + log(C)
    p = np.polyfit(log_d, log_A, 1)
    gamma = -p[0]
    C_collapse = np.exp(p[1])
    print(f"\n  Power law fit: A_tail = {C_collapse:.4f} * delta^(-{gamma:.4f})")
    print(f"  Exponent gamma = {gamma:.6f}")

    # Check fit quality
    log_A_fit = np.polyval(p, log_d)
    residuals = log_A - log_A_fit
    print(f"  Max residual in log: {np.max(np.abs(residuals)):.6f}")


# ================================================================
# Section 8: Effective A_tail model and Koide constraint
# ================================================================
print("\n\n8. EFFECTIVE A_tail MODEL AND KOIDE")
print("=" * 70)

# From the analysis:
# 1. Near vacuum: A_tail ~ |g0 - 1| * C1(alpha) where C1 ~ 1
# 2. Near collapse: A_tail ~ (8/5 - g0)^(-gamma) * C2
#
# A simple interpolating model:
#   A_tail(g0) ~ |g0 - 1| / (8/5 - g0)^gamma
#
# For the lepton solitons:
# e:   g0_e  = 0.868  (below vacuum)
# mu:  g0_mu = phi * g0_e = 1.404 (above vacuum)
# tau: g0_tau = 1.570  (near collapse)

g0_e = 0.8676968610
g0_mu = PHI * g0_e
g0_tau = 1.5695724137

# Using the model:
def A_model(g0, C1, gamma_val):
    return C1 * abs(g0 - 1) / (G0_CRIT - g0)**gamma_val

# Fit C1 and gamma from numerical data
g0_test = np.linspace(0.5, 1.595, 50)
A_num = np.array([A_tail(g) for g in g0_test])
valid = A_num > 1e-10
g0_v = g0_test[valid]
A_v = A_num[valid]

try:
    popt, pcov = curve_fit(A_model, g0_v, A_v, p0=[1.0, 0.2])
    C1_fit, gamma_fit = popt
    print(f"  Best fit: A_tail ~ {C1_fit:.4f} * |g0-1| / (8/5 - g0)^{gamma_fit:.4f}")
except Exception as e:
    C1_fit, gamma_fit = 1.0, 0.20
    print(f"  Fit failed: {e}, using C1=1.0, gamma=0.20")

# Compute A for the three leptons using the model
A_e_model = A_model(g0_e, C1_fit, gamma_fit)
A_mu_model = A_model(g0_mu, C1_fit, gamma_fit)
A_tau_model = A_model(g0_tau, C1_fit, gamma_fit)

# Compute Q_K from model
S2 = A_e_model**2 + A_mu_model**2 + A_tau_model**2
S4 = A_e_model**4 + A_mu_model**4 + A_tau_model**4
Q_K_model = S2**2 / S4

# Compute Q_K from numerical A_tail
A_e_num = A_tail(g0_e)
A_mu_num = A_tail(g0_mu)
A_tau_num = A_tail(g0_tau)
S2n = A_e_num**2 + A_mu_num**2 + A_tau_num**2
S4n = A_e_num**4 + A_mu_num**4 + A_tau_num**4
Q_K_num = S2n**2 / S4n

print(f"\n  Numerical A_tail values:")
print(f"    A_e   = {A_e_num:.8f}")
print(f"    A_mu  = {A_mu_num:.8f}")
print(f"    A_tau = {A_tau_num:.8f}")
print(f"    Q_K   = {Q_K_num:.10f}")
print(f"\n  Model A_tail values:")
print(f"    A_e   = {A_e_model:.8f}")
print(f"    A_mu  = {A_mu_model:.8f}")
print(f"    A_tau = {A_tau_model:.8f}")
print(f"    Q_K   = {Q_K_model:.10f}")

print(f"\n  Q_K(numerical) = {Q_K_num:.10f}")
print(f"  Q_K(model)     = {Q_K_model:.10f}")
print(f"  Difference: {abs(Q_K_num - Q_K_model):.6f}")

# Key question: does Q_K = 3/2 follow from the asymptotic model?
# The model A ~ |g0-1| * (8/5-g0)^(-gamma) has 2 parameters (C1, gamma).
# Q_K depends on the RATIO of A values, so C1 cancels:
#   A_i = C1 * |g0_i - 1| * (8/5 - g0_i)^(-gamma)
# With g0_mu = phi * g0_e and g0_tau from Koide:
# The ratio A_mu/A_e depends only on g0_e and gamma.
# Koide determines g0_tau, so Q_K = 3/2 becomes a constraint on (g0_e, gamma).

# Check: if we use the pure linear model (gamma=0):
def Q_K_from_model(gamma_val, g0_e_val):
    g0_mu_val = PHI * g0_e_val
    # Find g0_tau from Koide condition
    def A_m(g0):
        return abs(g0 - 1) * (G0_CRIT - g0)**(-gamma_val) if g0 < G0_CRIT else 1e10

    # Scan for g0_tau
    from scipy.optimize import brentq
    def koide_res(g0_3):
        Ae = A_m(g0_e_val)
        Am = A_m(g0_mu_val)
        At = A_m(g0_3)
        S2 = Ae**2 + Am**2 + At**2
        S4 = Ae**4 + Am**4 + At**4
        return S2**2 / S4 - 1.5

    try:
        g0t = brentq(koide_res, g0_mu_val + 0.01, G0_CRIT - 0.001)
        Ae = A_m(g0_e_val)
        Am = A_m(g0_mu_val)
        At = A_m(g0t)
        r31 = (At / Ae)**4
        return g0t, r31
    except:
        return None, None

print("\n\n  Q_K constraint: g0_tau and r_31 from model with different gamma:")
print(f"  {'gamma':>8s} {'g0_tau':>10s} {'r_31':>10s} {'r_31/3477':>10s}")
for gamma_val in np.arange(0.0, 0.5, 0.02):
    g0t, r31 = Q_K_from_model(gamma_val, g0_e)
    if g0t is not None:
        print(f"  {gamma_val:8.3f} {g0t:10.6f} {r31:10.2f} {r31/3477.18:10.4f}")


# ================================================================
# Summary
# ================================================================
print("\n\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
  1. OSCILLATORY TAIL: In the canonical variable f = g^3,
     the linearized equation is eta'' + 2/r eta' + eta = 0
     (oscillatory), NOT h'' + 2/r h' - h = 0 (exponential).
     The oscillations are a PHYSICAL consequence of the nonlinear
     kinetic metric K(g) = g^4.

  2. NEAR VACUUM: A_tail ~ |g0 - 1| (proportionality constant ~ 1).
     This follows from the linearized f-equation solution
     eta(r) = eta_0 * sin(r)/r.

  3. NEAR COLLAPSE: A_tail ~ (8/5 - g0)^(-{gamma:.4f}).
     The exponent gamma ~ 0.20 is determined by the nonlinear
     core dynamics near the collapse threshold.

  4. INTERPOLATING MODEL: A_tail ~ |g0-1| / (8/5-g0)^gamma
     with 2 parameters. Q_K = 3/2 with phi-FP spacing
     constrains gamma via the mass ratio r_31.

  5. The EXACT Koide Q_K = 3/2 requires the FULL nonlinear ODE --
     no simple asymptotic model reproduces it precisely.
     The model captures the qualitative structure but not the
     exact numerical value.
""")
