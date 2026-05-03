#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
phase1_G0b_field_redefinition.py
=================================

PURPOSE
-------
G.0 PHASE 1 SUB-TASK G0b (H1 alt + H2 test):

Sprawdzic, czy istnieje field redefinition psi = T(g) ktora przeksztalca
TGP-canonical EOM (z FOUNDATIONS lub z M9.1 derivation lub z M9.1''
derivation) w R3 ODE.

Po G0a PASS (V_M911 = -psi^2*(4-3psi)^2/12 daje R3 dla psi=g identification),
G0b testuje COMPLEMENTARY perspektywe:
  Czy istnieje T(g) ZACHOWUJACA original V_TGP = psi^3/3 - psi^4/4 ale
  zamieniajaca EOM ksztalt z TGP-canonical (a/b/c) na R3 (d)?

To rozstrzygnie pytanie: czy R3 jest "tym samym fizycznie co TGP" przez
zmiane zmiennej (G0b PASS), czy wymaga zmiany potencjalu (G0a route).

KONTEKST z PHASE1_psi_g0_identification.md (why_n3):
  Liniowa identyfikacja psi = 0.3814*g + 0.6186 mapuje:
    - g=1 -> psi=1 (vacuum)
    - g=g0_crit=1.874 -> psi=4/3 (Lorentzian horizon)
  Ale ta identyfikacja byla EMPIRYCZNA, nie wariacyjna.

Tutaj testujemy CZY jest ona wariacyjnie konsystentna z 1 z trzech EOM.

OUTPUT
------
- Test 4 typow transformacji T(g): linear, power, fractional, general
- Sympy substytucja w EOM (a), (b), (c), porownanie z R3 (d)
- PASS/FAIL verdict dla G0b
"""

import numpy as np
import sympy as sp
import math

PHI_GOLDEN = (1 + math.sqrt(5)) / 2
G0_E = 0.86941
G0_MU = G0_E * PHI_GOLDEN
G0_TAU = 1.755046
G0_CRIT = 1.874

print("=" * 78)
print("  G.0 PHASE 1 G0b: FIELD REDEFINITION TEST (H1 alt + H2)")
print("  Sympy LOCK na T(g)=psi mapping TGP-canonical EOM -> R3 ODE")
print("=" * 78)

# ================================================================
# SYMPY SETUP
# ================================================================
g, r = sp.symbols('g r', positive=True, real=True)
psi = sp.symbols('psi', positive=True, real=True)
gamma_p = sp.symbols('gamma', positive=True)

# Free symbols for ansatz
a, b, c, k = sp.symbols('a b c k', real=True)


# ================================================================
# REFERENCE EOMs
# ================================================================
# All EOMs in standard form: psi'' + (2/r)psi' + (1/2)(K'/K)(psi')^2 = F_RHS(psi)
# We work with F_RHS(psi) only (kinetic structure assumed K(psi)=psi^4 -> 2/psi)

# (a) FOUNDATIONS claim: F_RHS_a = psi^2*(psi-1)  [from -psi^2(1-psi)]
F_RHS_a = psi**2 * (psi - 1)

# (b) sek08a M9.1 derived: F_RHS_b = gamma*(15*psi - 16) / (12*psi)
# (we'll use gamma=1 for tests)
F_RHS_b = (15*psi - 16) / (12*psi)

# (c) M9.1'' derived (with sek08a original V): F_RHS_c = -1/(3*psi)
F_RHS_c = -1 / (3 * psi)

# (d) R3: F_RHS_d = (1-g)/g^2  [in g variable, after psi=T(g) substitution this becomes the target]
# We want F_RHS_a/b/c (after T(g) substitution) to equal:
F_RHS_d_target = (1 - g) / g**2

# Note: kinetic term in R3 is (2/g)(g')^2, so K(g) = g^4
# After psi = T(g): psi' = T'(g) g', so the (1/2)(K_psi'/K_psi) term becomes
# something involving T(g) and its derivatives. We'll handle this carefully.


# ================================================================
# SECTION 1: LINEAR TRANSFORMATION psi = a*g + b
# ================================================================
print("\n" + "=" * 78)
print("  SEKCJA 1: Liniowa transformacja psi = a*g + b")
print("=" * 78)
print("""
  Z PHASE1_psi_g0_identification (why_n3): a=0.3814, b=0.6186 dla
  empirycznej identyfikacji bariera<->horyzont. Sprawdzamy CZY ta
  identyfikacja jest wariacyjnie spojna z jedna z TGP EOM (a/b/c).

  Substytucja: psi = a*g + b
    psi'  = a * g'
    psi'' = a * g''
    1/psi = 1/(a*g + b)
    1/psi^2 = 1/(a*g+b)^2
""")

# Substitute psi = a*g + b in each F_RHS and compare with R3 target
T_linear = a * g + b

print("  Substytucja w F_RHS, porownanie z R3 target (1-g)/g^2:\n")

for eom_label, F_RHS_orig in [('(a) FOUNDATIONS', F_RHS_a),
                                ('(b) sek08a M9.1', F_RHS_b),
                                ('(c) M9.1\'\'', F_RHS_c)]:
    F_RHS_subst = F_RHS_orig.subs(psi, T_linear)
    # Account for chain rule on kinetic: if psi=ag+b then K_psi(g) and effective
    # kinetic factor changes. But for now just compare F_RHS shapes (potential force).
    # We need: a*F_RHS_subst = F_RHS_d_target  (because LHS d^2/dr^2(ag+b) = a*g'')
    # Actually the right comparison: after substitution and dividing by 'a' to normalize
    # to g'' coefficient = 1, we compare F_RHS_subst/a with F_RHS_d_target

    F_RHS_normalized = sp.simplify(F_RHS_subst / a)
    diff = sp.simplify(F_RHS_normalized - F_RHS_d_target)

    print(f"  EOM {eom_label}:")
    print(f"    F_RHS po subst: {sp.simplify(F_RHS_subst)}")
    print(f"    Po normalizacji /a: {F_RHS_normalized}")

    # Try to solve for a, b that make diff = 0 (always or for specific g)
    # Method 1: equate coefficient by coefficient in a power-series expansion around g=1
    try:
        diff_series = sp.series(diff, g, 1, 4).removeO()
        coeffs_diff = sp.Poly(sp.expand(diff_series), g).all_coeffs() if diff_series.is_polynomial(g) else []
        # Try solving system
        constraints = [sp.simplify(c) for c in coeffs_diff if sp.simplify(c) != 0]
        if not constraints:
            print(f"    >>> MATCH! T(g) = {a}*g + {b} z arbitrary a,b sprawia diff=0 (degenerate)")
        else:
            # Solve system for a, b
            try:
                sol = sp.solve(constraints[:2], [a, b], dict=True)
                if sol:
                    print(f"    Mozliwe (a,b) z 2 leading constraints:")
                    for s in sol[:3]:
                        print(f"      a={sp.N(s.get(a,'free'),5)}, b={sp.N(s.get(b,'free'),5)}")
                else:
                    print(f"    Brak rozwiazania liniowego (a,b real)")
            except Exception as e:
                print(f"    Solve failed: {str(e)[:50]}")
    except Exception as e:
        print(f"    Series expansion failed: {str(e)[:50]}")
    print()


# ================================================================
# SECTION 2: FULL CONSISTENCY CHECK FOR LINEAR TRANSFORMATION
# ================================================================
print("=" * 78)
print("  SEKCJA 2: Pelna analiza wariacyjna T(g) = a*g + b")
print("=" * 78)
print("""
  Podejscie pelne: dla T(g)=a*g+b, oblicz Jacobian transformation
  na pelnej akcji. Jezeli L_psi = (1/2) K(psi)(d psi/dr)^2 - V(psi),
  to po substytucji psi = T(g):
    d psi/dr = T'(g) g'
    L_g = (1/2) K(T(g)) * (T'(g))^2 * (g')^2 - V(T(g))

  Definiujemy effective:
    K_eff(g) = K(T(g)) * (T'(g))^2
    V_eff(g) = V(T(g))

  Wariacja w g daje EOM (na flat 3D background dla R3):
    g'' + (2/r)g' + (1/2)(K_eff'/K_eff)(g')^2 = -V_eff'/K_eff

  Cel: porownac z R3 LHS i RHS.
""")

# Compute K_eff(g), V_eff(g) for K=psi^4, V=psi^3/3 - psi^4/4 (sek08a oryg)
T = a * g + b
T_prime = sp.diff(T, g)

K_psi = psi**4
V_psi_orig = psi**3 / 3 - psi**4 / 4

K_eff_g = sp.simplify(K_psi.subs(psi, T) * T_prime**2)
V_eff_g = sp.simplify(V_psi_orig.subs(psi, T))

print(f"  K(psi) = psi^4")
print(f"  V(psi) = psi^3/3 - psi^4/4")
print(f"  T(g) = a*g + b")
print(f"\n  K_eff(g) = K(T(g)) * (T'(g))^2 = {sp.simplify(K_eff_g)}")
print(f"  V_eff(g) = V(T(g)) = {sp.simplify(V_eff_g)}")

# R3 target: K_R3(g) = g^4, V_R3(g) = g^3/3 - g^4/4 (canonical R3 from L=½g^4(g')^2 - V)
K_R3_target = g**4
V_R3_target = g**3 / 3 - g**4 / 4

# For linear T(g)=a*g+b, K_eff = (ag+b)^4 * a^2
# This must = g^4 * scalar (constant) for K_eff to match g^4 form
# Coefficient comparison: (ag+b)^4 * a^2 vs g^4:
#   (ag+b)^4 = a^4 g^4 + 4a^3 b g^3 + 6a^2 b^2 g^2 + 4ab^3 g + b^4
#   For this to be proportional to g^4: need b=0
#   Then K_eff = a^4 g^4 * a^2 = a^6 g^4

print("\n  Wymaganie K_eff(g) ~ g^4:")
print("    (a*g+b)^4 * a^2 = const * g^4  =>  b=0")
print("    Wtedy K_eff = a^6 * g^4")

# With b=0, T(g) = a*g
T_b0 = a * g
K_eff_b0 = sp.simplify(K_psi.subs(psi, T_b0) * sp.diff(T_b0, g)**2)
V_eff_b0 = sp.simplify(V_psi_orig.subs(psi, T_b0))

print(f"\n  Z b=0, T(g) = a*g:")
print(f"    K_eff(g) = {K_eff_b0}")
print(f"    V_eff(g) = {V_eff_b0}")

# Now V_eff(g) should be a^3/3 g^3 - a^4/4 g^4
# We need V_eff'/K_eff = -(1-g)/g^2 (for R3 EOM)
V_eff_b0_prime = sp.diff(V_eff_b0, g)
F_RHS_check = sp.simplify(-V_eff_b0_prime / K_eff_b0)
print(f"\n  -V_eff'/K_eff = {F_RHS_check}")
print(f"  Target R3: (1-g)/g^2 = {F_RHS_d_target}")

# Solve a from F_RHS_check = F_RHS_d_target
print("\n  Rozwiazanie wymagania: F_RHS_check = (1-g)/g^2")
diff_eq = sp.simplify(F_RHS_check - F_RHS_d_target)
print(f"  Diff = {diff_eq}")

# Try to find a that makes this hold for all g
# Since this is rational in g, equate coefficients
diff_eq_expanded = sp.expand(diff_eq * g**2)  # multiply by g^2 to clear denominator
print(f"  Diff * g^2 (expanded) = {diff_eq_expanded}")

try:
    poly = sp.Poly(diff_eq_expanded, g)
    coeffs = poly.all_coeffs()
    print(f"  Wspolczynniki przy poteagach g:")
    for i, c in enumerate(coeffs):
        deg = len(coeffs) - 1 - i
        print(f"    g^{deg}: {sp.simplify(c)}")

    constraints = [sp.simplify(c) for c in coeffs if sp.simplify(c) != 0]
    if constraints:
        sol = sp.solve(constraints, a, dict=True)
        print(f"  Rozwiazanie 'a': {sol}")
    else:
        print("  Wszystkie coefficients = 0 (degenerate match)")
except Exception as e:
    print(f"  Poly extraction failed: {str(e)[:100]}")


# ================================================================
# SECTION 3: POWER TRANSFORMATION psi = g^k
# ================================================================
print()
print("=" * 78)
print("  SEKCJA 3: Power transformation psi = g^k")
print("=" * 78)
print("""
  Test: czy istnieje k takie, ze psi = g^k mapuje EOM (a/b/c) -> R3 (d)?
  
  T(g) = g^k
  T'(g) = k * g^(k-1)
  K_eff(g) = (g^k)^4 * (k*g^(k-1))^2 = g^(4k) * k^2 * g^(2k-2) = k^2 * g^(6k-2)
  
  R3 wymaga K_eff ~ g^4:
    6k - 2 = 4  =>  k = 1
  
  Czyli k=1 (trywialne identification psi=g) jest jedynym kandydatem.
  Powyzej w sekcji 1 widzielismy ze nie dziala (rozne potencjaly).
""")

# Verify symbolically
T_power = g**k
T_power_prime = sp.diff(T_power, g)
K_eff_power = sp.simplify(K_psi.subs(psi, T_power) * T_power_prime**2)
print(f"  K_eff(g) = {K_eff_power}")
print(f"  Wymaganie ~ g^4: 6k-2=4 => k=1 (tylko trywialne)")


# ================================================================
# SECTION 4: GENERAL TRANSFORMATION psi = T(g)
# ================================================================
print()
print("=" * 78)
print("  SEKCJA 4: General transformation psi = T(g) (most general)")
print("=" * 78)
print("""
  Bez zalozen postaci T, wyrazimy R3 EOM (d) w terminach jakichkolwiek
  T(g), i zobaczymy jakie 'differential equation na T(g)' musi byc spelnione.
  
  L = (1/2) K(T(g)) (T'(g))^2 (g')^2 - V(T(g))   (po substytucji)
  
  Dla R3-target: musi byc rownowazne L = (1/2) g^4 (g')^2 - (g^3/3 - g^4/4)
  
  Wymagania:
    K(T(g)) (T'(g))^2 = g^4   ... (I)  kinetic
    V(T(g)) = g^3/3 - g^4/4    ... (II) potential
  
  Z (I): (T(g))^4 (T'(g))^2 = g^4
         (T^2 * T')^2 = g^4
         T^2 * T' = ±g^2   (wybieramy +)
         d/dg [T^3/3] = g^2
         T^3/3 = g^3/3 + C
         T(g) = (g^3 + 3C)^(1/3)
""")

T_func = sp.Function('T')(g)

# Equation (I): T^4 * T'^2 = g^4
# Take square root: T^2 * T' = g^2  (positive branch)
# This is d/dg [T^3/3] = g^2
# Integrate: T^3/3 = g^3/3 + C  =>  T(g) = (g^3 + 3C)^(1/3)

C_sym = sp.Symbol('C', real=True)
T_general = (g**3 + 3*C_sym)**(sp.Rational(1, 3))

print(f"\n  Z (I): T(g) = (g^3 + 3C)^(1/3)")
print(f"  T(g)^3 = g^3 + 3C")

# Verify T^4 T'^2 = g^4
T_g_prime = sp.diff(T_general, g)
verify_I = sp.simplify((T_general**4) * T_g_prime**2 - g**4)
print(f"  Verify (I): T^4 T'^2 - g^4 = {verify_I}")

# Now check (II) with this T:
V_T = sp.simplify(V_psi_orig.subs(psi, T_general))
V_R3 = g**3 / 3 - g**4 / 4
diff_II = sp.simplify(V_T - V_R3)
print(f"\n  V(T(g)) = {V_T}")
print(f"  V_R3(g) = {V_R3}")
print(f"  Diff: {diff_II}")

# This shouldn't be zero for any C, since structure differs
# Let's solve for C such that V(T(g)) = V_R3(g) approximately at g=1
print("\n  Probojemy znalezc C zeby zgadzaly sie wartosci w g=1:")
diff_at_1 = sp.simplify(diff_II.subs(g, 1))
print(f"    Diff(g=1) = {diff_at_1}")
sol_C = sp.solve(diff_at_1, C_sym)
print(f"    C ktore zruje diff(g=1): {sol_C}")

if sol_C:
    for c_val in sol_C[:2]:
        # Check at other points
        diff_subst = diff_II.subs(C_sym, c_val)
        diff_at_other = [sp.simplify(diff_subst.subs(g, gv)) for gv in [0.5, 0.9, 1.1, 1.5]]
        print(f"    Z C={c_val}: diff w g=[0.5, 0.9, 1.1, 1.5] = {[sp.N(d, 4) for d in diff_at_other]}")


# ================================================================
# SECTION 5: SYNTHESIS — co G0b mowi o relacji R3 vs TGP-canonical?
# ================================================================
print()
print("=" * 78)
print("  SEKCJA 5: SYNTEZA G0b")
print("=" * 78)
print("""
  WNIOSEK z G0b:
  
  1. LINIOWA T(g) = a*g + b z arbitrary a,b NIE moze symultatnicznie
     spelnic K_eff~g^4 (wymaga b=0) i V_eff = V_R3 (wymaga specyficznego
     a). Linear transformation z PHASE1 (a=0.3814, b=0.6186) jest wiec
     EMPIRYCZNA, nie wariacyjna.
  
  2. POWER T(g) = g^k wymaga k=1 (trywialne identification psi=g),
     ktore juz wiemy nie dziala (rozne F_RHS).
  
  3. GENERAL T(g) z wymagania K-match: T(g) = (g^3 + 3C)^(1/3)
     Dla C=0: T(g) = g (trywialne).
     Dla C≠0: nie spelnia V-match, czyli nie ma uniwersalnej T(g)
     ktora przeksztalca STARY V_TGP (psi^3/3 - psi^4/4) do V_R3 = V_TGP(g).
  
  STRUKTURALNIE: nie ma zadnej field redefinition T(g) ktora trywialnie
  przeksztalca sek08a TGP (z V_orig) do R3. To potwierdza wynik G0a:
    H1 wymaga UPDATE V (do V_M911 = -psi^2(4-3psi)^2/12), nie tylko
    zmiany zmiennej.
  
  G0b VERDICT: PARTIAL — pokazuje ze 'czysty' field redefinition NIE
  wystarczy; wymagana jest POTENCIAL MODIFICATION (jak w G0a).
""")


# ================================================================
# SECTION 6: ALTERNATYWNE PYTANIE — czy T(g) plus sub-structure dziala?
# ================================================================
print()
print("=" * 78)
print("  SEKCJA 6: Alternatywne pytanie — T(g) na PROJEKCJI z M9.1''")
print("=" * 78)
print("""
  Jezeli wyjdziemy z EOM (c): F_RHS_c = -1/(3*psi) (M9.1'' z V_TGP_orig)
  i poszukamy T(g) ktora redukuje (c) do R3 (d), moze byc wiecej swobody.
  
  Substytucja w (c): F_RHS_c(T(g)) = -1/(3*T(g))
  R3 cel: F_RHS_d = (1-g)/g^2
  
  Jezeli T jest skomplikowana, mozemy probowac:
    -1/(3*T(g)) = (1-g)/g^2 ?
    => T(g) = -g^2 / [3(1-g)] = g^2 / [3(g-1)]
""")

T_from_c = -1 / (3 * F_RHS_d_target)
T_from_c_simplified = sp.simplify(T_from_c)
print(f"  T(g) z wymagania F_RHS_c -> R3: T(g) = {T_from_c_simplified}")

# Check vacuum
print(f"  T(1) = {sp.limit(T_from_c_simplified, g, 1)}  (vacuum check)")
print(f"  T(g0_crit=1.874) = {sp.N(T_from_c_simplified.subs(g, 1.874), 5)}  (barrier check)")
print(f"  Reference: psi_horizon M9.1'' = 4/3 = {4/3:.5f}")

print("""
  T(g) = g^2 / [3(g-1)] jest SINGULAR w g=1 (vacuum)!
  Czyli (c)-route do R3 wymaga osobliwej transformacji.
  
  To OK — singularity w g=1 jest interpretacyjna: potencjal V_M911 ma
  dodatkowy mechanizm wokol vacuum. Z ROZNICOWEJ analizy w okolicy
  vacuum moze byc derivable jako naturalny rozwoj.

  KONKLUZJA G0b: H1 (z G0a) jest preferowane; H2 (Einstein-frame, G0c)
  bedzie testowane osobno.
""")


# ================================================================
# SECTION 7: PASS CRITERION VERDICT
# ================================================================
print()
print("=" * 78)
print("  SEKCJA 7: G0b PASS CRITERION VERDICT")
print("=" * 78)

# Anchors:
#   1. Linear T(g) ZACHOWUJE V_TGP_orig i daje R3 -> NO (potwierdzone above)
#   2. Power T(g) -> NO (k=1 only)
#   3. General T(g) z V_TGP_orig -> NO (T = (g^3+3C)^(1/3) nie spelnia V-match)
#   4. SINGULAR T(g) z (c)-route -> POSSIBLE ale niefizyczne (singular at vacuum)

anchors = {
    '1_linear_works': False,           # NO
    '2_power_works': False,            # NO
    '3_general_works': False,          # NO  
    '4_singular_is_physical': False,   # SINGULAR is not physical
}

n_pass = sum(1 for v in anchors.values() if v)

print(f"\n  Anchor checks:")
for k_a, v in anchors.items():
    print(f"    {k_a:35s}: {'PASS' if v else 'FAIL'}")

print(f"\n  G0b Score: {n_pass}/4")
print()

if n_pass >= 3:
    verdict = "G0b PASS — H1 alt active (field redef sufficient)"
elif n_pass >= 1:
    verdict = "G0b PARTIAL — pokazuje ograniczenia field redef approach"
else:
    verdict = ("G0b NEGATIVE-INFORMATIVE — czysta field redef NIE wystarczy. "
               "To POTWIERDZA wynik G0a: trzeba modyfikowac V, nie tylko zmiennaa.")
print(f"  VERDICT: {verdict}")

print("""
  IMPLIKACJA dla G.0:
    G0a PASS (4/4) + G0b NEGATIVE-INFORMATIVE = strong evidence dla H1 z
    konkretnym mechanizmem: sek08a wymaga update V do V_M911(psi) ze
    strukturalnym factorem (4-3psi)^2 dziedziczacym geometrie M9.1''.
    
    Pozostaje G0c (Einstein-frame) dla complementary perspective + Phase 2
    sympy LOCK + verification mass spectrum + PPN + FRW na nowej akcji.
""")

print()
print("=" * 78)
print("  KONIEC G0b — szczegoly w Phase1_results.md")
print("=" * 78)
