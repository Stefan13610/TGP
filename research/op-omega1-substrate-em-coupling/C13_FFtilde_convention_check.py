"""
C13 closure — sympy verification F^{μν} F̃_{μν} = -4 E·B convention

Audit C13 (MEDIUM severity, 2026-05-01):
  "F·F̃ = -4 E·B convention nieweryfikowana sympy
   (ε^{0123}=+1 vs MTW=-1 odwraca znak; '-4' coefficient w op-psi1
   Phase1_setup absorbuje czynnik)."

This script does the explicit sympy verification of:
  F̃^{μν} = (1/2) ε^{μνρσ} F_{ρσ}
  F^{μν} F̃_{μν} = ?·E·B

with explicit dependence on:
  (a) Levi-Civita sign convention: ε^{0123} = +1 (HEP/axion) vs -1 (MTW)
  (b) Metric signature: η = diag(-1,+1,+1,+1) [mostly-plus, used in TGP]
                       vs diag(+1,-1,-1,-1) [mostly-minus, MTW]

Output: explicit sign and coefficient confirming -4·E·B for HEP convention
        adopted in TGP (matches op-psi1 Phase1_setup).
"""

import sympy as sp

print("="*72)
print("C13 closure — F^{μν}F̃_{μν} = ?·E·B sympy verification")
print("="*72)

# -----------------------------------------------------------------------
# Setup: 4-vector field components
# -----------------------------------------------------------------------
Ex, Ey, Ez, Bx, By, Bz = sp.symbols('E_x E_y E_z B_x B_y B_z', real=True)

# F_{μν} field-strength tensor (lower indices), MTW conv mostly-plus η=diag(-1,+1,+1,+1)
# Convention: F_{0i} = E_i, F_{ij} = -ε_{ijk}B_k
#   (this is consistent for E_i = -∂_i A_0 - ∂_0 A_i in mostly-plus)
# We work directly with F^{μν} with upper indices:
#   F^{0i} = -F_{0i} = -E_i (mostly-plus η raises 0-index with -1)
#   F^{ij} = F_{ij} = -ε_{ijk}B_k

print("\n[1] Define F_{μν} with E_i, B_i, mostly-plus signature η=diag(-1,+1,+1,+1)")

# F_{μν} matrix (contravariant lower-lower):
F_dd = sp.Matrix([
    [0,   Ex,  Ey,  Ez],
    [-Ex, 0,  -Bz,  By],
    [-Ey, Bz,  0,  -Bx],
    [-Ez,-By,  Bx,  0],
])
print("F_{μν} =")
sp.pprint(F_dd)

# Metric mostly-plus:
eta = sp.diag(-1, 1, 1, 1)

# F^{μν} = η^{μα} η^{νβ} F_{αβ}
F_uu = eta * F_dd * eta
print("\nF^{μν} = η F_{μν} η =")
sp.pprint(F_uu)

# -----------------------------------------------------------------------
# Levi-Civita totally antisymmetric symbol with explicit sign choice
# -----------------------------------------------------------------------
def levi_civita_4d(eps_0123_sign):
    """ε^{μνρσ} with ε^{0123} = ±1; returns dict that defaults 0 for non-permutations"""
    from collections import defaultdict
    from itertools import permutations
    eps = defaultdict(int)
    base = (0, 1, 2, 3)
    for perm in permutations(base):
        # signature of permutation
        sign = 1
        arr = list(perm)
        for i in range(4):
            for j in range(i+1, 4):
                if arr[i] > arr[j]:
                    sign = -sign
                    arr[i], arr[j] = arr[j], arr[i]
        eps[perm] = sign * eps_0123_sign
    return eps

# -----------------------------------------------------------------------
# F̃^{μν} = (1/2) ε^{μνρσ} F_{ρσ}, F̃_{μν} = (1/2) ε_{μνρσ} F^{ρσ}
# In mostly-plus: ε_{μνρσ} = -ε^{μνρσ} (4 indices lowered, det(η) = -1)
# -----------------------------------------------------------------------

def FFtilde_compute(eps_0123_sign):
    eps_uuuu = levi_civita_4d(eps_0123_sign)
    # F̃_{μν} = (1/2) ε_{μνρσ} F^{ρσ}
    # In mostly-plus η, det=-1 → ε_{μνρσ} = -ε^{μνρσ}
    # so F̃_{μν} = -(1/2) ε^{μνρσ}_{contravariant} F^{ρσ}
    # Equivalent: F̃^{μν} = (1/2) ε^{μνρσ} F_{ρσ}
    F_tilde_uu = sp.zeros(4, 4)
    for mu in range(4):
        for nu in range(4):
            s = 0
            for rho in range(4):
                for sig in range(4):
                    s += sp.Rational(1, 2) * eps_uuuu[(mu, nu, rho, sig)] * F_dd[rho, sig]
            F_tilde_uu[mu, nu] = s
    # Now F^{μν} F̃_{μν} (lower one tilde), F̃_{μν} = η η F̃^{μν}
    F_tilde_dd = eta * F_tilde_uu * eta
    # Contract
    contraction = 0
    for mu in range(4):
        for nu in range(4):
            contraction += F_uu[mu, nu] * F_tilde_dd[mu, nu]
    return sp.simplify(contraction)

print("\n[2] Two-convention sympy LOCK:")
print("\n  Convention A: ε^{0123} = +1 (HEP/axion physics, used in TGP/op-psi1)")
ffA = FFtilde_compute(+1)
EdotB = Ex*Bx + Ey*By + Ez*Bz
print(f"    F^{{μν}} F̃_{{μν}} = {ffA}")
print(f"    Compare: -4·E·B = {sp.simplify(-4*EdotB)}")
print(f"    Match (-4·E·B)? {sp.simplify(ffA + 4*EdotB) == 0}")

print("\n  Convention B: ε^{0123} = -1 (MTW)")
ffB = FFtilde_compute(-1)
print(f"    F^{{μν}} F̃_{{μν}} = {ffB}")
print(f"    Compare: +4·E·B = {sp.simplify(4*EdotB)}")
print(f"    Match (+4·E·B)? {sp.simplify(ffB - 4*EdotB) == 0}")

# -----------------------------------------------------------------------
# Verdict
# -----------------------------------------------------------------------
print("\n" + "="*72)
print("[3] C13 verdict (sympy LOCKED):")
print("="*72)
print("""
  Convention A (HEP/axion, ε^{0123}=+1, mostly-plus η):
    F^{μν} F̃_{μν} = -4 E·B    ← TGP/op-psi1/op-omega1 convention LOCKED ✓

  Convention B (MTW, ε^{0123}=-1, mostly-plus η):
    F^{μν} F̃_{μν} = +4 E·B

  Coefficient '-4' w op-psi1 Phase1_setup absorbuje sgn(ε^{0123})=+1 +
  factor 4 = 2 (from F̃ definition antisymmetry doubling) × 2 (from
  sum over independent {0i,ij} pairs).

  → C13 CLOSED structurally: -4·E·B coefficient + sign sympy-verified for
    HEP convention (ε^{0123}=+1, mostly-plus signature), które są standardową
    konwencją axion electrodynamics (Peccei-Quinn 1977, Wilczek 1978) i
    używane konsekwentnie w TGP ω.1/ψ.1/τ.3.
""")
