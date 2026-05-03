#!/usr/bin/env python3
"""
omicron2_K_test.py — czy Phi tracks rho_bar(t) cosmologically?
================================================================

Test K kandydata z AUDIT_blind_spots.md (sek01/sek08a structural
inconsistency hypothesis):

  Sek01 ontology: "Phi_0 jest generowane przez kosmologiczna materie"
  Sek08a formalism: Phi_0 = vacuum minimum V'(Phi_0)=0, niezalezne od rho

  M10.5/de2 ignored matter source term in Phi-EOM. Quick estimate
  z AUDYTU sugerowal ze rzeczywiste cosmological evolution Phi_0(t)
  moze byc ~ Omega_m ~ 0.3, czyli HUGE (vs 1e-8 z Buchert variance).

Method:
  Solve coupled FRW Phi-EOM with matter source explicit:

    psi_dot = u
    u_dot = -3 H u - V_0 V'(psi) - s * Omega_m_a(a)
    a_dot = a H
    H = sqrt(Omega_r0/a^4 + Omega_m0/a^3 + rho_psi)

  where:
    rho_psi = (1/2) u^2 + V_0 V(psi)
    V(psi) = 4 psi^3 - 3 psi^4 (de2 normalized; V(1)=1, V'(1)=0, V''(1)=-12)
    V'(psi) = 12 psi^2 (1-psi)

  Matter source coupling 's':
    From sek08a action, source = q/Phi_0 * rho.
    Dimensional analysis (in dimensionless units) gives
    s_natural ~ 3/(2 Phi_0^2) for natural TGP coupling.

  For Phi_0 = 24.65 (post-UV.3): s_natural ~ 2.5e-3.

Three regimes scanned:
  s = 0       : pure de2 slow-roll (no source, baseline M10.1)
  s = small   : natural TGP (~2.5e-3)
  s = 0.01    : 4x natural
  s = 0.1     : 40x natural
  s = 0.3     : ~ Omega_m today (full tracker-like)

Verdict targets:
  (a) ψ -> 1 attractor wins for all s -> M10.5 confirmed
  (b) ψ tracks rho_bar for s_natural -> Hubble tension solver
  (c) ψ tracks for large s, frozen for small s -> shows where transition lies
"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
from scipy.integrate import solve_ivp


# ============================================================================
# Cosmological parameters (Planck 2018)
# ============================================================================
OMEGA_M0  = 0.315
OMEGA_R0  = 9.1e-5
OMEGA_DE0 = 1.0 - OMEGA_M0 - OMEGA_R0

PHI_0_TGP = 24.65  # post-UV.3 effective Phi_0 (IR scale)

# Hubble tension target
H0_SHOES = 73.04
H0_PLANCK = 67.40
DELTA_H = (H0_SHOES - H0_PLANCK) / H0_PLANCK  # ~0.0837

# Natural TGP matter coupling (dimensional analysis in dimensionless EOM)
S_NATURAL = 3.0 / (2.0 * PHI_0_TGP**2)  # ~ 2.47e-3


# ============================================================================
# Potential and EOM
# ============================================================================
def V(psi):
    """Normalized potential V(psi) = 4 psi^3 - 3 psi^4 (de2 convention)."""
    return 4.0*psi**3 - 3.0*psi**4

def dVdpsi(psi):
    """V'(psi) = 12 psi^2 (1 - psi)."""
    return 12.0 * psi**2 * (1.0 - psi)

def Omega_m_a(a):
    """Omega_m(z): rho_m / rho_crit_0 = Omega_m0 / a^3."""
    return OMEGA_M0 / a**3

def Omega_r_a(a):
    return OMEGA_R0 / a**4

def rho_psi(psi, u):
    """Dark energy density (dimensionless)."""
    return 0.5 * u**2 + V_0_global * V(psi)

def H_dimless(psi, u, a):
    """Hubble (dimensionless, H_0 = 1)."""
    rho_total = Omega_r_a(a) + Omega_m_a(a) + rho_psi(psi, u)
    return np.sqrt(max(rho_total, 1e-30))

def eom(t, y, s_coupling):
    """
    Coupled FRW Phi-EOM in cosmic time t (units 1/H_0).
    State y = [psi, u, a].
    """
    psi, u, a = y
    H = H_dimless(psi, u, a)

    # Phi-EOM with matter source: psi_dot = u, u_dot = -3 H u - V_0 V'(psi) - s rho_m
    psi_dot = u
    u_dot = -3.0 * H * u - V_0_global * dVdpsi(psi) - s_coupling * Omega_m_a(a)
    a_dot = a * H

    return [psi_dot, u_dot, a_dot]


# ============================================================================
# Solver
# ============================================================================
def solve_FRW(s_coupling, V_0=None, a_init=1e-6, delta_init=1e-3,
              method='RK45', rtol=1e-8):
    """
    Solve from a_init (z = 1/a_init - 1, very early) to a = 1 (today).
    Initial conditions: psi(t_init) = 1 - delta_init, u(t_init) = 0 (slow-roll).
    """
    global V_0_global
    if V_0 is None:
        V_0 = OMEGA_DE0  # initial guess
    V_0_global = V_0

    psi_init = 1.0 - delta_init
    u_init = 0.0

    # Initial cosmic time t_init:
    # In radiation/matter era, H ~ a^-2 / a^-3/2; integrate backward from t=now.
    # For solve_ivp we use cosmic time. Approximation: t_init small.
    # Let's parameterize by log(a) and convert.
    # Use solve_ivp with t span large enough to cover a_init -> 1.

    # We need t_today such that a(t_today) = 1.
    # Strategy: integrate forward in t from 0 to some t_max, find where a = 1.
    # For our purposes, integrate to t=2 (2 Hubble times) and interpolate.

    # Better: change variable to N = log(a). Then dN/dt = H, so
    #   dpsi/dN = u/H
    #   du/dN = (-3 H u - V_0 V'(psi) - s rho_m) / H
    #   N goes from N_init = log(a_init) to N_final = 0

    def eom_N(N, y, s_c=s_coupling):
        psi, u = y
        a = np.exp(N)
        H = H_dimless(psi, u, a)
        if H < 1e-30:
            return [0, 0]
        psi_dot = u / H
        u_dot = (-3.0*H*u - V_0_global*dVdpsi(psi) - s_c*Omega_m_a(a)) / H
        return [psi_dot, u_dot]

    N_init = np.log(a_init)
    N_final = 0.0
    N_eval = np.linspace(N_init, N_final, 5000)

    sol = solve_ivp(eom_N, [N_init, N_final], [psi_init, u_init],
                    args=(s_coupling,),
                    t_eval=N_eval, method=method, rtol=rtol, atol=1e-10)

    if not sol.success:
        return None

    return {
        'N': N_eval,
        'a': np.exp(N_eval),
        'z': 1.0/np.exp(N_eval) - 1.0,
        'psi': sol.y[0],
        'u': sol.y[1],
    }


def shoot_V0(s_coupling, target=OMEGA_DE0, delta_init=1e-3, a_init=1e-6,
             tol=1e-5, max_iter=50):
    """
    Find V_0 such that rho_psi(today) = target Omega_DE0.
    """
    V0_low = 0.1
    V0_high = 2.0

    for _ in range(max_iter):
        V0_mid = 0.5*(V0_low + V0_high)
        result = solve_FRW(s_coupling, V_0=V0_mid, delta_init=delta_init, a_init=a_init)
        if result is None:
            V0_high = V0_mid
            continue
        rho_today = rho_psi(result['psi'][-1], result['u'][-1])
        diff = rho_today - target
        if abs(diff) < tol:
            return V0_mid, result
        if diff > 0:
            V0_high = V0_mid
        else:
            V0_low = V0_mid
    return V0_mid, result


# ============================================================================
# Run scenarios
# ============================================================================
print("=" * 75)
print("omicron2 K-test: Czy Phi tracks rho_bar(t) cosmologically?")
print("=" * 75)
print(f"\n  Cosmology: Omega_m_0 = {OMEGA_M0}, Omega_DE_0 = {OMEGA_DE0:.4f}")
print(f"  Phi_0 (post-UV.3)    = {PHI_0_TGP}")
print(f"  Natural TGP coupling: s = 3/(2 Phi_0^2) = {S_NATURAL:.4e}")
print(f"  Hubble tension target: |Delta H_0/H_0| = {DELTA_H:.4f} ({DELTA_H*100:.2f}%)")

# Scenarios
scenarios = [
    ("s=0       (de2 baseline, no matter source)", 0.0),
    ("s=natural (TGP, ~2.5e-3)",                   S_NATURAL),
    ("s=0.01    (4x natural)",                     0.01),
    ("s=0.1     (40x natural)",                    0.1),
    ("s=0.3     (~Omega_m today, full tracker)",   0.3),
    ("s=1.0     (huge coupling, sanity check)",    1.0),
]

results_table = []

for label, s_value in scenarios:
    print(f"\n  Solving: {label}")
    V0_fitted, result = shoot_V0(s_value, delta_init=1e-3, a_init=1e-5)
    if result is None:
        print(f"    [FAILED to converge]")
        continue

    # Extract values at z=0 and z=1100
    z = result['z']
    psi = result['psi']
    u = result['u']

    psi_today = psi[-1]
    u_today = u[-1]

    # Find z=1100 index
    idx_recomb = np.argmin(np.abs(z - 1100.0))
    psi_recomb = psi[idx_recomb]

    # Find z=2 (just before DE dominance) for additional probe
    idx_z2 = np.argmin(np.abs(z - 2.0))
    psi_z2 = psi[idx_z2]

    # delta from 1
    d_today = psi_today - 1.0
    d_recomb = psi_recomb - 1.0
    d_z2 = psi_z2 - 1.0
    Delta_psi_cosmo = psi_today - psi_recomb

    # Implied Lambda variation (from V(psi) shift)
    # Lambda_eff ~ V_0 V(psi); V(1) = 1, V'(1) = 0, V''(1) = -12
    # V(psi) ≈ 1 - 6*delta^2 for small delta
    Lambda_today_norm = V(psi_today)  # / V(1) = V(psi_today)
    Lambda_recomb_norm = V(psi_recomb)
    dLL_cosmo = (Lambda_today_norm - Lambda_recomb_norm) / Lambda_today_norm

    # Implied Delta H_0 / H_0
    # ΔH²/H² = ΔΛ/Λ * Omega_L; ΔH/H = (1/2) ΔΛ/Λ * Omega_L
    delta_H_over_H = 0.5 * abs(dLL_cosmo) * OMEGA_DE0

    # Tracker ratio: how much does psi shift relative to required Omega_m
    if abs(d_recomb) > 1e-15:
        tracker_ratio = abs(d_recomb) / Omega_m_a(1/1101)  # at recomb a~1/1101
    else:
        tracker_ratio = 0

    print(f"    V_0 fitted             = {V0_fitted:.6f}")
    print(f"    psi(today)             = {psi_today:.6f}  (delta = {d_today:+.3e})")
    print(f"    psi(z=2)               = {psi_z2:.6f}  (delta = {d_z2:+.3e})")
    print(f"    psi(z=1100)            = {psi_recomb:.6f}  (delta = {d_recomb:+.3e})")
    print(f"    Delta psi (recomb→today) = {Delta_psi_cosmo:+.3e}")
    print(f"    dLambda/Lambda cosmo    = {dLL_cosmo:+.3e}")
    print(f"    |Delta H_0/H_0|         = {delta_H_over_H:.3e}")

    results_table.append({
        'label': label,
        's': s_value,
        'V_0': V0_fitted,
        'psi_today': psi_today,
        'psi_recomb': psi_recomb,
        'd_today': d_today,
        'd_recomb': d_recomb,
        'Delta_psi_cosmo': Delta_psi_cosmo,
        'dLL_cosmo': dLL_cosmo,
        'delta_H': delta_H_over_H,
    })


# ============================================================================
# Summary table
# ============================================================================
print("\n" + "=" * 75)
print("SUMMARY")
print("=" * 75)
print(f"\n  {'Scenario':<45} {'delta(today)':<15} {'delta(recomb)':<15} {'|dH/H|':<10}")
print(f"  {'-'*45} {'-'*15} {'-'*15} {'-'*10}")
for r in results_table:
    print(f"  {r['label']:<45} {r['d_today']:+.3e}   {r['d_recomb']:+.3e}   {r['delta_H']:.2e}")

print(f"\n  Required for Hubble tension:                |dH/H| = {DELTA_H:.4f}")

# Find scenario that crosses threshold
threshold_scenario = None
for r in results_table:
    if r['delta_H'] >= DELTA_H * 0.5:  # half of required is "interesting"
        threshold_scenario = r
        break

if threshold_scenario:
    print(f"\n  >>> First scenario with dH/H > 50% of required: s = {threshold_scenario['s']:.3e}")
    print(f"      Ratio s/s_natural = {threshold_scenario['s']/S_NATURAL:.1f}x")
else:
    print(f"\n  >>> NO scenario reaches 50% of required Hubble tension shift")


# ============================================================================
# Verdict
# ============================================================================
print("\n" + "=" * 75)
print("VERDICT")
print("=" * 75)

natural_result = next((r for r in results_table if abs(r['s'] - S_NATURAL) < 1e-5), None)
zero_result = next((r for r in results_table if r['s'] == 0.0), None)

if natural_result and zero_result:
    print(f"""
  Natural TGP coupling (s ~ 2.5e-3):
    Δψ cosmologically  = {natural_result['Delta_psi_cosmo']:+.3e}
    dLambda/Lambda      = {natural_result['dLL_cosmo']:+.3e}
    |dH_0/H_0|          = {natural_result['delta_H']:.3e}

  vs s=0 (no source, slow-roll only):
    Δψ cosmologically  = {zero_result['Delta_psi_cosmo']:+.3e}

  Natural matter coupling effect (vs s=0):
    Amplification factor = {abs(natural_result['Delta_psi_cosmo']/zero_result['Delta_psi_cosmo']) if abs(zero_result['Delta_psi_cosmo']) > 1e-30 else float('inf'):.1f}x
""")

# Determine which scenario applies
gap_natural = DELTA_H / max(natural_result['delta_H'], 1e-30) if natural_result else float('inf')
print(f"  Gap (required vs TGP natural): {gap_natural:.2e} = 10^{np.log10(gap_natural):.1f}")

if gap_natural > 1e3:
    print("""
  CONCLUSION: psi -> 1 attractor DOMINATES.
    Hubble friction (3H psi_dot) >> matter source (s rho_m).
    Cosmologically, psi stays locked near 1 — Phi_0 effectively constant.
    M10.5 verdict REINFORCED: TGP natural coupling cannot solve H_0 tension.

    User's K hypothesis: technically present in formalism but quantitatively
    suppressed by 1/Phi_0^2 ~ 1e-3 dimensionless coupling.
""")
elif gap_natural < 10:
    print("""
  CONCLUSION: TGP natural coupling COULD solve Hubble tension (within 1 dex).
    K hypothesis: matter coupling is sufficient to give cosmological Phi_0(t)
    evolution at level required for H_0 tension.
    M10.5 had blind spot — its calculation missed the matter source term.

    !!! POTENTIALLY EPOCHAL RESULT !!! Needs full audit.
""")
else:
    print(f"""
  CONCLUSION: PARTIAL match.
    TGP natural coupling gives effect within {gap_natural:.0f}x of required.
    Hubble tension partially explained by Phi_0(t) evolution, but not fully.
    Could combine with other effects (RG running, baryonic feedback, etc.).
""")

print("=" * 75)
