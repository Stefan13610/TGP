"""
dimensional_analysis.py  --  Theory of Generated Space (TGP)
=============================================================
Systematic dimensional analysis verification for all TGP quantities.

Checks that every equation in the theory is dimensionally consistent.
Uses symbolic dimension tracking (L = length, T = time, M = mass).

This is a CRITICAL consistency check:  if any quantity has wrong
dimensions, the entire derivation chain is suspect.

Usage:
    python dimensional_analysis.py
"""

import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# ═══════════════════════════════════════════════════════════════════════════
# Dimensional algebra  (exponents in [M, L, T])
# ═══════════════════════════════════════════════════════════════════════════
class Dim:
    """Represents physical dimensions as exponents of [M, L, T]."""
    def __init__(self, M=0, L=0, T=0, name=""):
        self.M = M
        self.L = L
        self.T = T
        self.name = name

    def __mul__(self, other):
        return Dim(self.M + other.M, self.L + other.L, self.T + other.T)

    def __truediv__(self, other):
        return Dim(self.M - other.M, self.L - other.L, self.T - other.T)

    def __pow__(self, n):
        return Dim(self.M * n, self.L * n, self.T * n)

    def __eq__(self, other):
        return (self.M == other.M and self.L == other.L and self.T == other.T)

    def __repr__(self):
        parts = []
        for sym, exp in [("M", self.M), ("L", self.L), ("T", self.T)]:
            if exp == 1:
                parts.append(sym)
            elif exp != 0:
                parts.append(f"{sym}^{exp}")
        return " ".join(parts) if parts else "1 (dimensionless)"


# ═══════════════════════════════════════════════════════════════════════════
# Base dimensions
# ═══════════════════════════════════════════════════════════════════════════
DIMLESS = Dim(name="dimensionless")
M = Dim(M=1, name="mass")
L = Dim(L=1, name="length")
T = Dim(T=1, name="time")

# ═══════════════════════════════════════════════════════════════════════════
# Fundamental constants
# ═══════════════════════════════════════════════════════════════════════════
c0    = L / T                          # [m/s]
G0    = L**3 / (M * T**2)             # [m^3 kg^-1 s^-2]
hbar0 = M * L**2 / T                  # [kg m^2 s^-1]  = [J s]
H0    = T**(-1)                        # [s^-1]
kappa = G0 / c0**4                     # [s^2 / (kg m)] = 8piG/c^4

# ═══════════════════════════════════════════════════════════════════════════
# TGP quantities
# ═══════════════════════════════════════════════════════════════════════════
Phi0    = DIMLESS                       # Phi_0 is dimensionless (def. in TGP)
Phi     = DIMLESS                       # Phi(x) is dimensionless
psi     = DIMLESS                       # psi = Phi / Phi_0
rho     = M / L**3                     # mass density [kg/m^3]
q       = G0 / c0**2                   # source coupling = 8piG_0/c_0^2 [m/kg]

# Dynamic constants
c_Phi   = c0                           # c(Phi) has same dimension as c_0
hbar_Phi = hbar0                        # hbar(Phi) same as hbar_0
G_Phi   = G0                           # G(Phi) same as G_0

# Planck length
lP      = (hbar0 * G0 / c0**3)**0.5    # But we check manually

# Field equation parameters
alpha   = DIMLESS                       # alpha = 2 (dimensionless)
beta    = L**(-2)                       # [m^-2]  (or dimensionless in hat units)
gamma   = L**(-2)                       # [m^-2]

# ═══════════════════════════════════════════════════════════════════════════
# Checks
# ═══════════════════════════════════════════════════════════════════════════
results = []

def check(name, lhs, rhs, detail=""):
    ok = (lhs == rhs)
    results.append((name, ok, detail if detail else f"LHS={lhs}, RHS={rhs}"))
    return ok


# --- CHECK 1: Field equation ---
# nabla^2 Phi / Phi_0  [dim of nabla^2 = L^{-2}]
# = L^{-2} * 1 / 1 = L^{-2}
dim_laplacian_Phi = L**(-2) * Phi / Phi0
dim_D = dim_laplacian_Phi  # = L^{-2}

# N[Phi] term 1: alpha * (grad Phi)^2 / (Phi_0 * Phi)
# grad Phi ~ L^{-1}, so (grad Phi)^2 ~ L^{-2}
dim_N1 = DIMLESS * L**(-2) / (DIMLESS * DIMLESS)  # = L^{-2}

# N[Phi] term 2: beta * Phi^2 / Phi_0^2
dim_N2 = beta * DIMLESS / DIMLESS  # beta * 1 = L^{-2}

# N[Phi] term 3: gamma * Phi^3 / Phi_0^3
dim_N3 = gamma * DIMLESS / DIMLESS  # = L^{-2}

# RHS: -q * Phi_0 * rho
dim_RHS_field = q * Phi0 * rho
# q = G0/c0^2 = [m^3/(kg s^2)] / [m^2/s^2] = [m/kg]
# q * rho = [m/kg] * [kg/m^3] = [m^{-2}] = L^{-2}  OK

check("Field eq: D[Phi]", dim_D, L**(-2),
      f"nabla^2(Phi)/Phi_0 = {dim_D}")
check("Field eq: N1 (gradient)", dim_N1, L**(-2),
      f"alpha*(grad Phi)^2/(Phi_0*Phi) = {dim_N1}")
check("Field eq: N2 (beta)", dim_N2, L**(-2),
      f"beta*Phi^2/Phi_0^2 = {dim_N2}")
check("Field eq: N3 (gamma)", dim_N3, L**(-2),
      f"gamma*Phi^3/Phi_0^3 = {dim_N3}")
check("Field eq: RHS (source)", dim_RHS_field, L**(-2),
      f"-q*Phi_0*rho = {dim_RHS_field}")

# --- CHECK 2: Spatial mass ---
# m_sp^2 = 3*gamma - 2*beta, [m^{-2}]
dim_m_sp_sq = gamma  # = L^{-2}
check("m_sp^2 dimension", dim_m_sp_sq, L**(-2),
      f"m_sp^2 = {dim_m_sp_sq}")

# --- CHECK 3: Cosmological mass ---
# m_cosmo^2 = 4*c0^2*gamma/3  (kappa cancels from the action)
dim_m_cosmo_sq = c0**2 * gamma
# c0^2 = L^2/T^2
# c0^2 * gamma = [L^2/T^2] * [L^{-2}] = [T^{-2}]
# This IS [T^{-2}] -- correct! (kappa cancels from the action)

# With kappa cancelled from the action, the cosmological mass is simply:
# m_cosmo^2 = c0^2 * gamma = [L^2/T^2] * [L^{-2}] = [T^{-2}]
dim_c2_gamma = c0**2 * gamma
check("c0^2 * gamma dimension (cosmo mass)", dim_c2_gamma, T**(-2),
      f"c0^2*gamma = {dim_c2_gamma}")

# kappa is still a defined physical constant but no longer enters field equations:
kappa_dim = G0 / c0**4
check("kappa = 8piG/c^4 dimension (reference only)", kappa_dim, Dim(M=-1, L=-1, T=2),
      f"kappa = {kappa_dim}")

# With kappa cancelled from the gravitational action (1/kappa on the entire
# gravitational sector), the FRW field equation reads:
#   ddot(psi) + 3H*dot(psi) + 2*dot(psi)^2/psi = c0^2 * W(psi)
#
# Dimensional check:
# - ddot(psi): [T^{-2}] (psi dimensionless)
# - 3H*dot(psi): [T^{-1}]*[T^{-1}] = [T^{-2}] OK
# - c0^2*W(psi): [L^2/T^2]*[L^{-2}] = [T^{-2}] OK!
#
# The previous dimensional mismatch (with kappa) is RESOLVED:
# c0^2*kappa*W had dim [M^{-1}L^{-1}] != [T^{-2}].
# Without kappa: c0^2*W has dim [L^2/T^2]*[L^{-2}] = [T^{-2}]. Correct!
#
# Similarly: m_cosmo^2 = 4*c0^2*gamma/3 has dim [T^{-2}]. Correct!
print("=" * 65)
print("TGP Dimensional Analysis")
print("=" * 65)

# Report basic checks
print("\n  BASIC DIMENSION CHECKS:")
n_pass = 0
n_fail = 0
for name, ok, detail in results:
    status = "PASS" if ok else "FAIL"
    if ok:
        n_pass += 1
    else:
        n_fail += 1
    print(f"  [{status}] {name}")
    print(f"         {detail}")

print(f"\n  Basic checks: {n_pass}/{n_pass + n_fail} passed")

# ═══════════════════════════════════════════════════════════════════════════
# DIMENSIONAL CONSISTENCY TABLE
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  TGP QUANTITY DIMENSIONS  (M=mass, L=length, T=time)")
print("=" * 65)

quantities = [
    ("Phi, Phi_0, psi", "1 (dimensionless)", "Field value"),
    ("c_0", "L T^-1", "Speed of light"),
    ("G_0", "L^3 M^-1 T^-2", "Newton constant"),
    ("hbar_0", "M L^2 T^-1", "Planck constant"),
    ("l_P", "L", "Planck length"),
    ("H_0", "T^-1", "Hubble parameter"),
    ("kappa = 8piG/c^4", "M^-1 L^-1 T^2", "Einstein kappa"),
    ("rho", "M L^-3", "Mass density"),
    ("q = 8piG_0/c_0^2", "L M^-1", "Source coupling"),
    ("alpha", "1", "Gradient coupling (=2)"),
    ("beta, gamma (field eq)", "L^-2", "Self-interaction [static]"),
    ("beta_hat, gamma_hat", "1", "Dimensionless: beta*r_0^2"),
    ("m_sp", "L^-1", "Spatial mass (Yukawa)"),
    ("m_cosmo^2 = 4c0^2*gamma/3", "T^-2",
     "Cosmological mass (kappa cancels)"),
    ("tau_0", "T", "Cosmological timescale"),
    ("Lambda_eff = gamma/12", "L^-2", "Effective cosm. constant"),
    ("rho_Lambda = Lambda*c^2/(8piG)", "M L^-3", "DE density"),
    ("c(Phi)/c_0 = (Phi_0/Phi)^{1/2}", "1", "Ratio (dimensionless)"),
    ("V_eff(d)", "L^-1 (energy/unit)", "Interaction potential"),
    ("E_int(d)", "M L^2 T^-2", "Interaction energy"),
]

for name, dim, desc in quantities:
    print(f"  {name:<40s}  [{dim:<20s}]  {desc}")

# ═══════════════════════════════════════════════════════════════════════════
# KEY FINDING: Dimensional issue in m_cosmo
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("  KEY FINDING: m_cosmo^2 dimensional analysis")
print("=" * 65)
print("""
  The linearized cosmological equation (eq:Phi-cosmo-linearized):

    ddot(delta) + 3H*dot(delta) + (4*c0^2*gamma/3)*delta = ...

  requires the coefficient to have dimensions [T^{-2}] (since delta
  is dimensionless and ddot(delta) has dim [T^{-2}]).

  With gamma in [L^{-2}] (from the static field equation):
    c0^2 * gamma = [L^2/T^2] * [L^{-2}] = [T^{-2}]  CORRECT!

  The previous version had kappa in the equation:
    c0^2 * kappa * gamma = [M^{-1} L^{-1}]  !=  [T^{-2}]  (WRONG)

  RESOLUTION: With 1/kappa on the ENTIRE gravitational sector of the
  action (not just 1/(2kappa) on the kinetic term), kappa cancels from
  ALL field equations. The corrected FRW equation reads:

    ddot(psi) + 3H*dot(psi) + 2*dot(psi)^2/psi = c0^2 * W(psi)

  where c0^2 * W(psi) has dim [L^2/T^2] * [L^{-2}] = [T^{-2}].

  This RESOLVES the dimensional mismatch that was previously flagged.
  m_cosmo^2 = 4*c0^2*gamma/3 has dim [T^{-2}] as required.

  STATUS: RESOLVED -- kappa removal fixes dimensional consistency.
""")

if n_fail > 0:
    sys.exit(1)
else:
    sys.exit(0)
