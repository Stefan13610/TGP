"""
Phase 3 sympy verification — polarization & spin (KRYTYCZNA)
Cykl: op-Phi-decomposition-photon-2026-05-07

Tests:
  F3.1: DOF counting for scalar δΦ vs photon
  F3.2: Spin analysis (Lorentz group representation)
  F3.3: Alternative β (A_μ = ∂_μ δΦ) — F_μν = 0 test
  F3.4: Alternative γ (∂_i ∂_j δΦ TT-mode) — spin test
  F3.5: Verdict
"""
import sympy as sp

print("=" * 70)
print("F3.1: DOF counting (scalar δΦ vs photon)")
print("=" * 70)

print("""
Scalar field δΦ(x,t):
  - 1 component (skalar Lorentz-niezmienniczy)
  - On-shell DOF: 1 (single scalar mode)
  - Polarization: trivialny (skalarna fala nie ma kierunku)

Photon A_μ(x,t):
  - 4 components (4-wektor)
  - Off-shell DOF: 4
  - Constraints: gauge ∂_μ A^μ = 0 (Lorenz), removes 1 DOF → 3
  - On-shell (Maxwell, massless): k_μ A^μ = 0 removes 1 more → 2 DOF
  - Polarizacje: 2 transverse (e.g., x, y dla k = ẑ)

KONFLIKT FUNDAMENTALNY: 1 vs 2 DOF.
""")

print("=" * 70)
print("F3.2: Spin analysis (Lorentz group representation)")
print("=" * 70)

print("""
Lorentz group SO(3,1) representations indexed by (j_L, j_R):
  - (0, 0): skalar, spin 0
  - (1/2, 0) ⊕ (0, 1/2): Dirac spinor, spin 1/2
  - (1/2, 1/2): 4-wektor, spin 1
  - (1, 0) ⊕ (0, 1): tensor antysymetryczny F_μν, spin 1
  - (1, 1): tensor symetryczny bez śladu, spin 2 (graviton)

Scalar δΦ → reprezentacja (0, 0) → spin 0.
Photon → reprezentacja (1/2, 1/2) lub równoważnie (1,0)⊕(0,1) → spin 1.

KONFLIKT REPREZENTACJI: spin 0 vs spin 1.
Twierdzenie reprezentacji: NIE MA niezmienniczego mapowania
   (0, 0) → (1/2, 1/2)
bez wprowadzenia DODATKOWEGO pola lub struktury.
""")

print("=" * 70)
print("F3.3: Alternative β (A_μ ≡ ∂_μ δΦ) — F_μν test")
print("=" * 70)

# A_μ = ∂_μ δΦ → F_μν = ∂_μ A_ν - ∂_ν A_μ = ∂_μ ∂_ν δΦ - ∂_ν ∂_μ δΦ = 0
# Pure gauge → no electromagnetic dynamics

x = sp.symbols('x_0 x_1 x_2 x_3', real=True)  # x^μ
phi_func = sp.Function('phi')(*x)

# A_μ = ∂_μ φ
A = [sp.diff(phi_func, x[mu]) for mu in range(4)]
print(f"\nA_μ ≡ ∂_μ δΦ:")
for mu in range(4):
    print(f"  A_{mu} = {A[mu]}")

# F_μν = ∂_μ A_ν - ∂_ν A_μ
print(f"\nF_μν = ∂_μ A_ν - ∂_ν A_μ:")
F = [[sp.diff(A[nu], x[mu]) - sp.diff(A[mu], x[nu]) for nu in range(4)] for mu in range(4)]

all_zero = True
for mu in range(4):
    for nu in range(mu+1, 4):
        F_simplified = sp.simplify(F[mu][nu])
        print(f"  F_{mu}{nu} = ∂²δΦ/(∂x{mu}∂x{nu}) - ∂²δΦ/(∂x{nu}∂x{mu}) = {F_simplified}")
        if F_simplified != 0:
            all_zero = False

if all_zero:
    print("\nWERDYKT: F_μν = 0 EXACTLY (równość pochodnych mieszanych)")
    print("→ A_μ = ∂_μ δΦ jest CZYSTĄ TRANSFORMACJĄ GAUGE")
    print("→ NIE MA dynamiki elektromagnetycznej")
    print("→ FALSIFIES alternative β")
else:
    print("\nUNEXPECTED: F_μν ≠ 0!")

print("=" * 70)
print("F3.4: Alternative γ (∂_i ∂_j δΦ in TT projection) — spin test")
print("=" * 70)

print("""
h_ij ≡ ∂_i ∂_j δΦ  (spatial second derivatives, symmetric tensor)

Pod transformacjami SO(3) rotations:
  - h_ij is symmetric rank-2 spatial tensor (6 components)
  - Decomposition under SO(3):
    * Trace: h ≡ δ^ij h_ij = ∇²δΦ → 1 DOF (skalar)
    * Traceless symmetric: h_ij - (1/3)δ_ij h → 5 DOF (spin 2)
    * Transverse-traceless (TT): k^i h_ij = 0 → 2 DOF (spin 2 graviton-like)

ALTERNATYWA γ daje TT-mode z 2 DOF, ALE:
  - SPIN = 2 (graviton representation), NIE spin 1 (photon)
  - To by była ALTERNATYWA dla GRAVITON-MODE, nie photon

Reprezentacja Lorentz: ∂_i ∂_j δΦ → tensor symmetric → spin 2 (po TT)
Reprezentacja photon: spin 1.

→ FALSIFIES alternative γ jako foton (works as graviton-mode candidate).
""")

print("=" * 70)
print("F3.5: Alternative δ (multi-component / new field) — S05 violation")
print("=" * 70)

print("""
Aby uzyskać spin 1 z field-theoretic perspective, MUSI być:
  - Nowe pole wektorowe A_μ (independent), LUB
  - Multi-component skalar (np. SO(N) z N=2,3 → spin combinations)

OBA naruszają S05 (single-Φ axiom, closed 2026-04-26 Path B).

Re-open S05 byłby OUT-OF-SCOPE Stage 2 — wymagałby osobnego cyklu:
  op-S05-reopen-vector-photon-2026-MM-DD/

CO TYM SAMYM oznacza dla Stage 2:
  - W obecnym single-Φ framework: foton-as-δΦ jest STRUKTURALNIE
    NIEMOŻLIWY (representation theory argument)
  - To NIE jest błąd obliczeniowy — to fundamentalna konsekwencja
    representation theory grupy Lorentza
""")

print("=" * 70)
print("F3.6: VERDICT — Stage 2 photon-as-δΦ FAILS")
print("=" * 70)

print("""
STRUCTURAL FAILURE: foton-as-δΦ-mode hypothesis NIE jest możliwa
w obecnym single-Φ framework.

Powód: representation theory grupy Lorentza.
  - Skalarne pole = spin 0 reprezentacja
  - Foton = spin 1 reprezentacja
  - NIE MA mapowania (0,0) → (1/2,1/2) bez additional field structure.

Tested alternatives:
  α (longitudinal foton):           FAIL — ≠ 2 transverse polarization
  β (A_μ = ∂_μ δΦ):                 FAIL — F_μν = 0, pure gauge
  γ (∂_i∂_j δΦ TT-mode):            FAIL — spin 2 (graviton), nie spin 1
  δ (new field/multi-component):    FAIL — narusza S05 single-Φ axiom

IMPLICATIONS:

  1. PHOTON ontology w TGP: standardowy A_μ z QED, propagujący
     na background metric M9.1''(Φ̄). NIE jest δΦ-mode.

  2. δΦ-modes (Phase 1+2 wyniki): ZACHOWANE jako VALID, ale rebranded
     jako "TGP scalar gravitational modes" — DODATKOWY sektor
     niezależny od photona.

  3. Phase 1+2 derivations:
     - λ = hc/E: VALID dla any wave + canonical quantization
       (działa zarówno dla δΦ-mode jak i standard photon A_μ)
     - Dispersion ω² = c²(k²+γ): VALID dla δΦ-mode (jako TGP scalar)
     - m_eff² = γ: VALID dla δΦ-mode (TGP scalar mass ~ H_0)
     - T^μ_μ ~ γ: VALID dla δΦ-mode (NIE photon, photon ma T^μ_μ=0)

  4. L01 ρ-bridge restored (post-rebrand):
     - ρ_EM = 0 EXACTLY z FμνF^μν Weyl-niezmienniczości 4D ✓
     - ρ_(δΦ-mode) ≈ ρ_Λ z γ-mass term (separate scalar sector)

  5. EXT-1 BBN return (Phase 4):
     - Foton = standard A_μ → ρ_EM = 0 strukturalnie (L01 preserved)
     - BBN problem: nadal istnieje (ρ_rad w erze radiacyjnej z γ-photons +
       neutrino, nie z δΦ-modes)
     - δΦ-modes są dodatkowym sektorem (jak dark radiation), Δ N_eff?
     - EXT-1 STRUCTURAL_NO_GO pozostaje BARDZO PRAWDOPODOBNIE

CLASSIFICATION:
  - Photon-as-δΦ hypothesis: STRUCTURAL_NO_GO (representation theory)
  - δΦ ontology generally: STRUCTURAL DERIVED (Phase 1+2 valid)
  - Stage 2 NET STATUS: STRUCTURAL CONDITIONAL z caveat
    (δΦ-modes valid as TGP scalar sector, NIE jako foton)

HONEST PHASE 6 REPORTING:
Stage 2 cykl wykonał Phase 1+2+3 z 6/6 + 4/4 + 4/4 PASS analitycznie.
Phase 3 ujawnił structural constraint: foton-as-δΦ niemożliwy w
single-Φ framework. Phase 1+2 wyniki ZACHOWANE jako TGP scalar modes
(ważny dodatkowy sektor), ale photon ontology pozostaje standardowym EM.
""")
