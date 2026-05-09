"""
Phase 1 N2 - Darwin Lagrangian z TGP delta_Phi retarded propagator

Cel: sprawdzic czy TGP scalar field z M9.1'' background reprodukuje:
  (a) Coulomb-like static potential (gravity-like)
  (b) Darwin-like velocity correction (magnetism-like)
  (c) Identify gap (jesli istnieje) miedzy TGP scalar i standard QED gauge field

Test: pure scalar field MASSIVE Klein-Gordon na M9.1'' background
"""
import sympy as sp
from sympy import (symbols, Function, sqrt, exp, pi, integrate, oo,
                   simplify, expand, factor, Matrix, diff, cos, sin)

print("=" * 72)
print("Phase 1 N2 - Darwin Lagrangian z TGP delta_Phi")
print("=" * 72)

# Symbols
r, t, c, m, q1, q2, gamma_param = symbols('r t c m q_1 q_2 gamma', positive=True)
v1, v2 = symbols('v_1 v_2', real=True)
hbar = symbols('hbar', positive=True)

# =============================================================================
# F1: Linearized delta_Phi EOM
# =============================================================================
print("\n--- F1: Linearized delta_Phi EOM (z op-Phi-decomposition-photon) ---")

print("""
Z op-Phi-decomposition-photon Phase 1, linearized EOM dla delta_Phi
wokol tla Phi_bar (psi=1):

  Box delta_Phi + gamma * delta_Phi = -q * Phi_0 * delta_rho

Box = (1/c^2) d^2/dt^2 - Laplacian      (mostly-plus signature)

To jest MASSIVE Klein-Gordon equation z m_eff^2 c^4 / hbar^2 = gamma * c^2
m_eff = hbar * sqrt(gamma) / c
""")

# =============================================================================
# F2: Static Green's function (Yukawa)
# =============================================================================
print("\n--- F2: Static Green's function (Yukawa potential) ---")

print("""
Static limit (d/dt = 0):
  -Laplacian G_static + gamma * G_static = delta^3(r)

Solution (standard):
  G_static(r) = exp(-sqrt(gamma) * r) / (4*pi*r)

For gamma = 0: G_static = 1/(4*pi*r)  (pure Coulomb / Newton)
For gamma > 0: Yukawa-like, exponential cutoff at distance ~1/sqrt(gamma)

W TGP: gamma ~ H_0^2/c^2, sqrt(gamma) ~ H_0/c ~ 10^-26 /m
Na laboratoryjnych skalach (r ~ 1m): sqrt(gamma)*r << 1
=> G_static ≈ 1/(4*pi*r)  (Coulomb-like effectively)
""")

# Symboliczne G_static
m_eff = sqrt(gamma_param)
G_static = exp(-m_eff * r) / (4*pi*r)
print(f"G_static(r) = exp(-sqrt(gamma)*r) / (4*pi*r) = {G_static}")

# Limit gamma -> 0
G_static_massless = sp.limit(G_static, gamma_param, 0)
print(f"\nLimit gamma -> 0: G_static -> {G_static_massless}")
print(f"To jest Coulomb potential 1/(4*pi*r) ✓")

# =============================================================================
# F3: Static interaction energy (two solitons)
# =============================================================================
print("\n--- F3: Static interaction energy two solitons ---")

print("""
Energia interakcji 2 sol w static limit:
  E_int_static = q_1 * q_2 * G_static(|r_1 - r_2|)
              = q_1 q_2 / (4*pi*r) * exp(-sqrt(gamma)*r)

Dla TGP gamma ~ H_0^2: praktycznie identyczne z Coulomb
(albo Newton, w zaleznosci od interpretation of q jako charge / mass).

KLUCZOWY POINT:
  "Coulomb" w TGP scalar field daje attractive lub repulsive force
  zaleznie od signs q_1, q_2.

  W standard EM: q_1 q_2 same sign -> repulsion
  W gravity: m_1 m_2 always positive -> attraction

  W TGP single scalar Phi: efective "charge" wynika z Phi-coupling.
  Pytanie: jaki sign?

  Jesli q_1 q_2 zawsze positive (jak masa) -> attraction = gravity-like
  Jesli moze byc obu znakow -> EM-like
""")

E_int_static = q1 * q2 * G_static
print(f"\nE_int_static = q_1 q_2 G_static = {E_int_static}")

# =============================================================================
# F4: Retarded Green's function (full)
# =============================================================================
print("\n--- F4: Retarded Green's function (full) ---")

print("""
Pełna retarded Green's function dla MASSIVE Klein-Gordon:

  G_ret(r, t) = theta(t) * [delta(t - r/c)/(4*pi*r) -
                m * theta(t^2 - r^2/c^2) * J_1(m*sqrt(t^2 - r^2/c^2)) /
                (4*pi*sqrt(t^2 - r^2/c^2))]

gdzie m = sqrt(gamma)*c (effective mass), J_1 to Bessel function.

W praktyce: 'sharp' part (delta function) + 'tail' z J_1.

Static limit recovery: integral nad t daje Yukawa.

Dla MOVING SOURCE z velocity v:
  source: rho(x, t) = delta^3(x - v*t)

  Retarded delta_Phi at point r' (test position):
    delta_Phi(r', t) = q * G_ret(|r' - v*t_ret|, t - t_ret)

  z retarded condition:
    t_ret = t - |r' - v*t_ret|/c   (implicit)
""")

# =============================================================================
# F5: Slow-motion expansion (CRITICAL TEST)
# =============================================================================
print("\n--- F5: Slow-motion expansion v/c <<1 (KRYTYCZNY) ---")

print("""
Dla v << c, expansion w (v/c):

LEADING ORDER (static): V_0 = q_1 q_2 / (4*pi*r)
  -> CONFIRMED: TGP daje Coulomb / Newton-like potential ✓

NEXT-TO-LEADING (v/c):
  Standard EM (gauge field A_mu): Darwin term
    L_Darwin = q_1 q_2 / (8*pi*r) * [v_1 . v_2 + (v_1.r_hat)(v_2.r_hat)]

  Pochodzi z VECTOR vector potential A_mu coupling z J_mu.
  Specific dla SPIN-1 gauge field.

  TGP scalar field (no vector A): co dostajemy?
""")

# Sketch derivation dla scalar field
print("""
Dla scalar field, source jest:
  source = rho(x,t)  (only scalar component)

NIE MA J_mu vector source.

Retarded propagator daje:
  delta_Phi_at_2 from 1:
    delta_Phi ~ q_1 / |r_2 - r_1(t_ret)|
    z relativistic correction (1 - n.v/c) factor

Energia interakcji 2 sol z velocities:
  E_int(r, v_1, v_2) = q_1 q_2 / (4 pi r) * [1 + correction(v_1, v_2)]
""")

# =============================================================================
# F6: Konkretne sprawdzenie - czy emerguje Darwin-like term?
# =============================================================================
print("\n--- F6: Sprawdzenie konkretne - Darwin term emergence? ---")

print("""
Standard derivation (Jackson, Sec. 12.7):
  Lienard-Wiechert potentials dla scalar field:
    delta_Phi(r, t) = q / (4*pi |r - r_s(t_ret)| (1 - n.v_s/c))

  Expand do O((v/c)^2) z source moving slowly:
    delta_Phi ≈ q/(4*pi*r) * [1 + (n.v_s)/c + (1/2)((v_s)^2 + (n.v_s)^2)/c^2 + ...]

Energia interakcji z test charge q_2 z velocity v_2 (moving):
  E = q_2 * delta_Phi - to jest only POSITION coupling

  Plus: kinetic terms for both particles in their own frames

  Plus: pole gradient force on moving particle: F_2 = -nabla V

KEY OBSERVATION:
  Dla SCALAR field, V_eff zalezy od:
    - r (position separation)
    - v_s (source velocity, jak source generates field)
    - NIE bezposredni v_2 dependence! (test particle nie sprzega z phase)

  To jest INNE niz standard EM (vector field), gdzie:
    - V_eff zawiera (v_1 . v_2) / c^2  -> Darwin
    - Lorentz force F = q(E + v x B) - velocity-dependent BY DESIGN

WERDYKT F6:
  Pure scalar field NIE daje natywnie Darwin (v_1 . v_2)/c^2 term.
  Ten term jest specific dla gauge field A_mu structure.

  TGP z czysto skalarnym Phi NIE reprodukuje standard EM kinematics.
""")

# =============================================================================
# F7: Phase rotation rescue?
# =============================================================================
print("\n--- F7: Phase rotation rescue mechanism? ---")

print("""
User's claim: 'magnetism = phase-rotated gravity'

Mathematical interpretation candidates:

OPTION A: Phase from M9.1'' background
  M9.1'' gives g_tt = -c^2(4-3psi)/psi
  Time component depends on local psi value
  Moving observer sees different g_tt -> phase shift

  But this is COORDINATE effect, nie physical phase rotation.

OPTION B: Phase from coupled spinor (z op-SPIN-SU2)
  Soliton z spinor S w external delta_Phi field
  Larmor-like precession
  Phase = SU(2) rotation angle

  This IS physical mechanism, ale wymaga that B-field couples
  z S poprzez zewnetrzny mechanizm (gauge field A_mu z Stage 2!)

OPTION C: Phase z time-dependent delta_Phi background
  delta_Phi(r, t) ma natywna oscylacje
  Test particle sees alternating gradient
  Net effect: drift force perpendicular to gradient

  This MIGHT give effective B-like structure.
  Ale wymaga that delta_Phi background oscillates - what drives it?

WERDYKT F7:
  'Phase rotation' jako simple mechanism w czysto skalarnym TGP
  jest NIE OBVIOUS i wymaga formal derivation.

  Najobiecujsza Option B (przez spinor coupling) - ale to znaczy ze:
    - B-field ontology jest standard A_mu (z Stage 2)
    - delta_Phi-mediated 'magnetism' jest separate phenomenon
    - NIE 'unified' z gravity w sense literal single mechanism
    - Unification jest CONCEPTUAL (jeden Phi field), nie operational
""")

# =============================================================================
# F8: Honest verdict
# =============================================================================
print("\n--- F8: HONEST VERDICT N2 ---")

print("""
Werdykt na podstawie F1-F7:

WHAT TGP SCALAR DELIVERS:
  ✓ Static Coulomb / Newton-like potential (Yukawa for nonzero gamma)
  ✓ Gravity-like attraction (jezeli q masa-positive)
  ✓ Massive Klein-Gordon dynamics
  ✓ Retarded propagator standard

WHAT TGP SCALAR DOES NOT DELIVER:
  ✗ Darwin (v_1 . v_2)/c^2 term (specific dla vector A_mu)
  ✗ Magnetic field B = nabla x A (no vector field)
  ✗ Lorentz force F = qv x B (no velocity-dep coupling natively)

UNIFICATION CLAIM RE-EVALUATION:

Autor's claim: 'magnetism = phase-rotated gravity'

W swiatle F1-F7:
  - Pure scalar Phi NIE daje EM-magnetism kinematics natively
  - 'Phase rotation' jako rescue wymaga additional structure
  - Najobiecujsza ścieżka: SPINOR COUPLING (z op-SPIN-SU2 N18)
  - Ale wtedy: B-field wciaz jest A_mu (Stage 2), 'magnetism' to
    Larmor precession z Lorentz force jako standard QM result

CONCLUSION:
  TGP unifikuje gravity z magnetism CONCEPTUALLY (jeden Phi field
  jest fundamental dla obu), ale OPERATIONALLY:
    - Gravity emerguje z M9.1''(Phi_bar)
    - Magnetism emerguje z A_mu (Stage 2) PLUS coupling z TGP spinor

  Czyli unifikacja jest na poziomie ontology (single field substrate)
  ale NIE na poziomie kinematics (oddzielne mechanisms).

WERDYKT N2: PARTIAL FAIL dla literal kinematic unification.
            PARTIAL SUCCESS dla ontological unification.

PROBABILITY UPDATE:
  Pełen DERIVED kinematic unification: 5-10% (very low)
  Pełen DERIVED ontological unification: 30-40%

  Interpretacja autor's claim 'silniejsza wersja grawitacji' MOZE
  wciaz byc valid w sense gravitomagnetic effects (Lense-Thirring
  in TGP framework), ale to jest GR-type effect, nie EM-type.
""")

print("=" * 72)
print("N2 ANALYSIS COMPLETE - PARTIAL FAIL FOR LITERAL UNIFICATION")
print("Honest result: TGP scalar NIE reprodukuje EM kinematics natively")
print("=" * 72)
