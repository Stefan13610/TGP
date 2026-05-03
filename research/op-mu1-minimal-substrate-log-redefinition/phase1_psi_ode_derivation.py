#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
μ.1 P1.1 — Sympy derivation ψ-ODE z g-ODE pod redefinicją ψ ≡ log g

Cel: pokazać analitycznie, że pod podstawieniem g(r) = exp(ψ(r)),
g-ODE przechodzi w smooth, well-posed ψ-ODE bez patologii.

g-ODE (oryginalna R3):
    g'' + (α/g)(g')² + (2/r)g' = (1-g)·g^(2-2α)

Pod g = exp(ψ):
    g'  = exp(ψ)·ψ'
    g'' = exp(ψ)·(ψ'' + (ψ')²)

Spodziewany wynik:
    ψ'' + (1+α)(ψ')² + (2/r)ψ' = (e^(-ψ) - 1)·e^((1-2α)ψ)
"""

import sympy as sp


def derive_psi_ode():
    print("=" * 78)
    print("  μ.1 P1.1 — Sympy derivation ψ-ODE z g-ODE")
    print("=" * 78)
    print()

    # Symbols
    r, alpha = sp.symbols('r alpha', positive=True, real=True)
    psi = sp.Function('psi')(r)

    # Definitions
    g_psi = sp.exp(psi)
    gp = sp.diff(g_psi, r)        # g' = exp(ψ)·ψ'
    gpp = sp.diff(g_psi, r, 2)    # g'' = exp(ψ)·(ψ'' + ψ'²)

    print("Oznaczenia:")
    print(f"  g(r)  = exp(ψ(r))")
    print(f"  g'(r) = {sp.simplify(gp)}")
    print(f"  g''(r) = {sp.simplify(gpp)}")
    print()

    # Original g-ODE: g'' + (α/g)(g')² + (2/r)g' - (1-g)·g^(2-2α) = 0
    # Express g^(2-2α) explicit jako exp((2-2α)ψ) żeby sympy zwijał
    g_pow = sp.exp((2 - 2 * alpha) * psi)  # = g^(2-2α) = exp(ψ)^(2-2α)
    g_ode_lhs = gpp + (alpha / g_psi) * gp**2 + (2 / r) * gp \
                - (1 - g_psi) * g_pow

    print("g-ODE LHS (substituted g = exp(ψ), g^(2-2α) → exp((2-2α)ψ)):")
    print(f"  {sp.simplify(g_ode_lhs)}")
    print()

    # Divide by exp(ψ) to get ψ-ODE
    psi_ode_lhs = sp.simplify(g_ode_lhs / g_psi)
    psi_ode_lhs = sp.powsimp(sp.expand(psi_ode_lhs), force=True)

    print("ψ-ODE LHS (after dividing by exp(ψ)):")
    print(f"  {psi_ode_lhs}")
    print()

    # Try to express in canonical form: ψ'' + (1+α)(ψ')² + (2/r)ψ' = RHS
    psi_pp = sp.diff(psi, r, 2)
    psi_p = sp.diff(psi, r)

    # Expected ψ-ODE
    expected_lhs = psi_pp + (1 + alpha) * psi_p**2 + (2 / r) * psi_p \
                   - (1 - sp.exp(psi)) * sp.exp((1 - 2 * alpha) * psi)
    expected_lhs_simplified = sp.powsimp(sp.expand(expected_lhs), force=True)

    print("Expected ψ-ODE LHS:")
    print(f"  {expected_lhs_simplified}")
    print()

    # Verify they are identical (force exp combination)
    diff = sp.simplify(psi_ode_lhs - expected_lhs_simplified)
    diff = sp.powsimp(diff, force=True)
    diff = sp.simplify(diff)
    print(f"Difference (should be 0): {diff}")
    print()

    if diff == 0:
        print("  ✓ DERIVATION CONFIRMED: ψ-ODE is")
        print()
        print("    ψ'' + (1+α)(ψ')² + (2/r)ψ' = (1 - e^ψ)·e^((1-2α)ψ)")
        print()
    else:
        print(f"  ✗ DIFFERENCE NON-ZERO — algebra needs review")

    return diff == 0


def analyze_psi_ode_for_alpha_2():
    """For α=2 analyse RHS and fixed point structure."""
    print("=" * 78)
    print("  α=2 specific analysis")
    print("=" * 78)
    print()

    psi = sp.symbols('psi', real=True)

    # RHS dla α=2: (1 - e^ψ)·e^(-3ψ)
    alpha_val = 2
    rhs = (1 - sp.exp(psi)) * sp.exp((1 - 2 * alpha_val) * psi)
    rhs_simplified = sp.expand(rhs)
    print(f"RHS dla α=2: {rhs_simplified}")
    print(f"           = e^(-3ψ) - e^(-2ψ)")
    print()

    # Fixed point ψ=0
    rhs_at_0 = rhs.subs(psi, 0)
    print(f"RHS w ψ=0: {rhs_at_0}  (oczekiwane 0 — vacuum fixed point)")
    print()

    # Linearization wokół ψ=0
    rhs_linear = sp.series(rhs, psi, 0, 2).removeO()
    print(f"Linearyzacja wokół ψ=0: RHS ≈ {sp.simplify(rhs_linear)}")
    print(f"  (ψ ≈ 0, RHS ≈ -ψ — restoring force, vacuum stabilny)")
    print()

    # Asymptotic for large ψ (g → ∞)
    print("Asymptotyka dla α=2:")
    print(f"  ψ → +∞ (g → ∞):  RHS = e^(-3ψ) - e^(-2ψ) → 0⁻ (oba człony → 0)")
    print(f"  ψ → -∞ (g → 0):  RHS = e^(-3ψ) - e^(-2ψ) → +∞ (singular)")
    print()
    print("UWAGA: ψ → -∞ odpowiada g → 0 (zerowanie metryki).")
    print("Fizyczne soliton solutions z g₀=1.77 → ψ₀=0.57 → ψ ∈ [0, ψ₀] — bezpieczne.")
    print()


def check_alpha_1_special_case():
    """For α=1 (Bridge theorem case) sprawdzić czy ODE jest specjalna."""
    print("=" * 78)
    print("  α=1 special case (Bridge theorem)")
    print("=" * 78)
    print()

    psi = sp.symbols('psi', real=True)

    rhs_alpha1 = (1 - sp.exp(psi)) * sp.exp((1 - 2 * 1) * psi)
    rhs_alpha1_expanded = sp.expand(rhs_alpha1)
    print(f"RHS dla α=1: {rhs_alpha1_expanded}")
    print(f"           = e^(-ψ) - 1")
    print()

    # Linearization
    rhs_linear = sp.series(rhs_alpha1, psi, 0, 3).removeO()
    print(f"Linearyzacja: {sp.simplify(rhs_linear)}")
    print()

    # ODE ma czysto wykładniczą strukturę dla α=1 — może być integrable
    print("OBSERWACJA: dla α=1 RHS = e^(-ψ) - 1, czysta exp-funkcja (BEZ r-zależności!)")
    print("            ma najprostszą exp-strukturę z całej rodziny α — możliwe że")
    print("            ODE jest integrable. To zgodne z bridge theorem (α=1 critical).")
    print()


def main():
    success = derive_psi_ode()
    if not success:
        print("DERIVATION FAILED — przerywam")
        return

    analyze_psi_ode_for_alpha_2()
    check_alpha_1_special_case()

    print("=" * 78)
    print("  P1.1 verdict")
    print("=" * 78)
    print()
    print("  ψ-ODE derivation: ✓ PASS (sympy potwierdza algebrę)")
    print("  α=2 RHS analysis:  vacuum stabilny, soliton range ψ≥0 bezpieczny")
    print("  α=1 special case: czysta exp-struktura, possibly integrable")
    print()
    print("  → przechodzimy do P1.2 (numerical soliton)")
    print()
    print("=" * 78)


if __name__ == "__main__":
    main()
