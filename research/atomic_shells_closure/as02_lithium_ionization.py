"""
as02_lithium_ionization.py — Li IE₁/IE₂/IE₃ vs TGP

Trzy testy falsyfikacyjne:

TEST 1. Ciąg jonizacyjny Li jako 3 hydrogeniczne "stany wirtualne":
   IE₁ (odrywanie 2s¹) = 5.3917 eV
   IE₂ (odrywanie z 1s²) = 75.640 eV
   IE₃ (hydrogenic Li²⁺) = 122.454 eV

   Dla każdego wyliczamy Z_eff z E = -Z_eff²·Ry/n².
   Czy Z_eff's tworzą jakiś prosty wzór TGP?

TEST 2. Koide-like dla Li IE:
   K_Li = (IE₁+IE₂+IE₃) / (√IE₁+√IE₂+√IE₃)²
   Leptony naładowane: K = 2/3 dokładnie. Jeśli Li IE chain ma
   strukturę solitonowo-drabinową (g₀ w progression φ?), oczekujemy K=2/3.
   W przeciwnym razie: K ≠ 2/3 → IE Li NIE są "drabiną" jak leptony.

TEST 3. Czy istnieje g₀-drabina dla alkali IE₁ (Li, Na, K, Rb, Cs)?
   Każdy alkali ma 1 valence elektron w ns¹. IE₁ w TGP powinien skalować jak
   A_valence⁴. Porównujemy z obserwowanymi.

   Dodatkowo: testujemy WIELKOŚĆ i ZNAK A_valence (2s Li) w odniesieniu
   do empirycznego A_s = -0.111 z SC paper.

TEST 4. Hipoteza tail-phase:
   Hydrogenic ψ_2s(r) = (1/√8π)(1/a_0)^(3/2) · (2-r/a_0) · exp(-r/2a_0)
   Ma NODE w r=2a_0, ZNAK ogona jest UJEMNY (2-r/a_0 < 0 dla r>2a_0).
   Propozycja: A_s(TGP-atom) ~ -|A_e|·f(screen) < 0 z powodu tail-sign flip.
   Weryfikacja: czy numeryczna amplituda A_tail_2s(Li) ≈ -0.111?
"""

import math
import sys
import io
import numpy as np
from scipy.integrate import quad

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ---------------------------------------------------------------------------
# Stałe
# ---------------------------------------------------------------------------
Ry = 13.605693122994        # eV
alpha_fs = 7.2973525693e-3
a_0 = 5.29177210903e-11     # m
Hartree_eV = 27.21138624598

# Soliton tail amplitude dla wolnego elektronu (z mass_scaling_k4/R5)
A_e_free = 0.1246           # TGP reference

# Empirical A_orb z SC/ρ(T) (ps5_global_fit, r13)
A_s_SC   = -0.111
A_sp_SC  = +0.207
A_d_SC   = +0.310
A_f_SC   = +2.034

# ---------------------------------------------------------------------------
# Dane eksperymentalne
# ---------------------------------------------------------------------------
# Li ciąg jonizacyjny (NIST ASD)
IE_Li = np.array([5.39171495, 75.64018, 122.45437])  # eV

# Alkali IE₁ + hydrogenic quantum number valency
alkalis = [
    ("H",  1,  1, 13.59844),   # for reference (1s)
    ("Li", 3,  2,  5.39171),   # 2s
    ("Na", 11, 3,  5.13908),   # 3s
    ("K",  19, 4,  4.34066),   # 4s
    ("Rb", 37, 5,  4.17713),   # 5s
    ("Cs", 55, 6,  3.89390),   # 6s
    ("Fr", 87, 7,  4.07272),   # 7s (relativistic!)
]

print("=" * 70)
print("  as02 — Lithium ionization chain vs TGP")
print("=" * 70)

# ---------------------------------------------------------------------------
# TEST 1: Z_eff dla ciągu IE Li
# ---------------------------------------------------------------------------
print("\n[1] Li IE chain as hydrogenic with Z_eff:")
print(f"    {'stage':<8}{'IE [eV]':>12}{'n':>4}{'Z_eff':>12}{'δ_n (defect)':>16}")
for i, (IE, n) in enumerate(zip(IE_Li, [2, 1, 1])):  # IE1→2s, IE2→He-like 1s, IE3→H-like 1s
    Zeff_sq = IE / Ry * n**2
    Zeff = math.sqrt(Zeff_sq)
    # True Z felt at infinity for valence (not for inner shells!)
    if i == 0:  # 2s valence
        Z_asymptotic = 1
        n_star = Z_asymptotic / math.sqrt(IE/Ry)
        delta = n - n_star
    else:
        delta = None  # not defined for inner shells
    delta_str = f"{delta:.4f}" if delta is not None else "n/a"
    print(f"    IE_{i+1:<5}{IE:>12.5f}{n:>4}{Zeff:>12.5f}{delta_str:>16}")

print(f"    IE₃ = {IE_Li[2]:.3f} eV;  hydrogenic 13.6·3² = {Ry*9:.3f} eV;  diff {100*(IE_Li[2]-Ry*9)/(Ry*9):.3f}%")
print("    → IE₃ jest DOKŁADNIE hydrogenic (zgodnie z as01). Sanity: OK.")
print("    → IE₂ wymaga Z_eff ≈ {:.2f} (<3 ze względu na screening 1s-elektronu)".format(
    math.sqrt(IE_Li[1]/Ry)))
print("    → IE₁ wymaga Z_eff ≈ 1.26 (quantum defect δ_s ≈ 0.41)")

# ---------------------------------------------------------------------------
# TEST 2: Koide K dla Li IE
# ---------------------------------------------------------------------------
print("\n[2] Koide-style relation for Li IE:")
sum_IE = IE_Li.sum()
sum_sqrt_IE = np.sqrt(IE_Li).sum()
K_Li = sum_IE / sum_sqrt_IE**2
print(f"    K_Li = Σ IE / (Σ √IE)² = {sum_IE:.4f} / {sum_sqrt_IE**2:.4f} = {K_Li:.6f}")
print(f"    Leptony naładowane (e, μ, τ): K = 2/3 = 0.666667")
print(f"    Li IE: K = {K_Li:.6f}, różnica od 2/3: {K_Li - 2/3:+.4f} ({100*(K_Li-2/3)/(2/3):+.1f}%)")

# Test Brannen ansatz θ scan
# √IE_i = a(1 + √2·cos(θ + 2π·i/3))
print("\n    Brannen √2 ansatz test: Czy Li IE chain pasuje do √m_i = a(1+√2·cos(θ+2πi/3))?")
sqrt2 = math.sqrt(2)
sorted_sqIE = np.sort(np.sqrt(IE_Li))  # zrosnąco
best_res = np.inf; best_theta = None
for theta in np.linspace(0, 2*math.pi, 10000):
    cos_vals = [math.cos(theta + 2*math.pi*i/3) for i in range(3)]
    # skala: a = avg(√IE) / (1 + √2·avg(cos)) — ale trzeba dopasować
    # wymaga sort; dopasowujemy przez permutację i skalę
    for perm in [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]]:
        model = np.array([1 + sqrt2*cos_vals[perm[i]] for i in range(3)])
        if np.all(model > 0):
            a_fit = sorted_sqIE.sum() / model.sum()
            pred = a_fit * model
            res = np.sum((pred - sorted_sqIE)**2) / np.sum(sorted_sqIE**2)
            if res < best_res:
                best_res = res; best_theta = theta
print(f"    Best residual²/norm² = {best_res:.4e}  (leptons: ~10⁻¹⁵ dokładnie)")
print(f"    → Li IE chain NIE obeysuje Brannen ansatz. Nie jest drabiną solitonową.")

# ---------------------------------------------------------------------------
# TEST 3: Alkali series — czy g₀-drabina dla ns¹ valence?
# ---------------------------------------------------------------------------
print("\n[3] Alkali series: czy A_valence skalą czteropotęgowo?")
print(f"    {'atom':<4}{'Z':>4}{'n':>4}{'IE₁ [eV]':>12}{'Z_eff':>10}{'(IE/Ry)^(1/4)':>18}")
for name, Z, n, IE in alkalis:
    Zeff = math.sqrt(IE/Ry) * n
    A_proxy = (IE/Ry)**(0.25)   # hipoteza: IE ∝ A⁴ → A ∝ IE^(1/4)
    print(f"    {name:<4}{Z:>4}{n:>4}{IE:>12.5f}{Zeff:>10.5f}{A_proxy:>18.5f}")

# Hypothesis: if TGP IE ∝ A^4, then A ∝ IE^(1/4).
# For Li: A_Li = (5.39171/13.6057)^(1/4) = 0.794 → zupełnie inne niż |A_s|=0.111
# Chyba że A_proxy jest w JEDNOSTKACH A_e=0.1246:
# A_Li_in_Ae = A_proxy / 1 (for H at IE=Ry)... ale H IE=13.6, (13.6/13.6)^(1/4)=1
# Nie uzgadnia się z A_s = -0.111 prostą drogą.

print("\n    A_proxy(Li) = (IE_Li/Ry)^(1/4) = 0.7938")
print("    Porównanie do empirycznego |A_s|(SC) = 0.111 — RZĄD WIELKOŚCI INNY (7×)")
print("    → 'IE ∝ A⁴' dla alkali NIE daje wartości zgodnych z A_s(SC).")
print("    → A_s(SC) jest w innej skali — prawdopodobnie NIE jest A_tail elektronu izolowanego")

# ---------------------------------------------------------------------------
# TEST 4: Hipoteza tail-phase z hydrogenicznej 2s
# ---------------------------------------------------------------------------
print("\n[4] Tail-phase hypothesis: ψ_2s node flips sign at r = 2a_0")
print("    Hydrogenic 2s: ψ_2s(r) = (1/√(8π))(1/a_0)^(3/2)·(2 - r/a_0)·exp(-r/2a_0)")
print("    NODE at r = 2a_0 → tail (r→∞) has NEGATIVE sign")
print()

# Dla Li 2s z Z_eff=1.26:
# ψ_2s^Li(r) ≈ (Z/√(8π))^(3/2) · (2 - Zr/a_0) · exp(-Zr/2a_0), Z=Z_eff
Zeff_Li = 1.259
def psi_2s(r_in_a0, Z=1.0):
    """Hydrogenic 2s wavefunction (a_0 units). Returns r·R_2s (radial prob amp)."""
    rho = Z * r_in_a0
    R = (1/math.sqrt(2)) * Z**1.5 * (1 - rho/2) * math.exp(-rho/2)  # normalized r·R_2s/sqrt(r)
    return R

# Dla Li: wyliczamy amplitude przy r = 3a_0 (za nodem, w regime "tail")
r_test_vals = [3.0, 5.0, 10.0]  # a_0 units
print(f"    Li 2s (Z_eff={Zeff_Li:.3f}) radial amplitude R_2s(r) at tail points:")
for r in r_test_vals:
    psi = psi_2s(r, Z=Zeff_Li)
    print(f"      r = {r:>4.1f} a_0:  R_2s = {psi:>+.6e}  (sign: {'negative' if psi<0 else 'positive'})")

# Porównaj znak z A_s(SC) = -0.111
print(f"\n    Znak A_s(SC) = -0.111 (UJEMNY)")
print(f"    Znak tail R_2s(Li) = UJEMNY (za nodem)")
print(f"    → ZGODA JAKOŚCIOWA: oba znaki ujemne!")

# Ale magnituda: R_2s w a.u. przy r=5a_0 ≈ -0.04 — też inna niż 0.111
# Stąd: A_s(SC) NIE JEST dosłownie wartością R_2s(r) w atomie, ale:
# 1) znak jest zachowany, 2) jest ~O(0.1) skala.

# ---------------------------------------------------------------------------
# TEST 5: Alkali ns¹ tail-sign → prediction dla A_orb per klasa
# ---------------------------------------------------------------------------
print("\n[5] Alkali ns¹ tail-sign prediction:")
print("    n=2 (Li): 2s ma 1 radial node → tail sign NEGATIVE")
print("    n=3 (Na): 3s ma 2 radial nodes → tail sign POSITIVE (2 flips)")
print("    n=4 (K):  4s ma 3 radial nodes → NEGATIVE")
print("    n=5 (Rb): 5s ma 4 radial nodes → POSITIVE")
print("    n=6 (Cs): 6s ma 5 radial nodes → NEGATIVE")
print()
print("    Jeśli A_s(SC) = -0.111 byłoby tylko dla Li, przewidywalibyśmy")
print("    PRZESKOKI ZNAKU wzdłuż kolumny s metali... ale w SC paper A_s = -0.111")
print("    traktowane jest JAKO STAŁA dla wszystkich s-metali (Hg, Tl, Sn dowodne).")
print()
print("    → hipoteza 'A_s = R_2s tail sign' ZAWODZI: nie wyjaśnia stałości A_s")
print("    → A_s(SC) musi odnosić się do STANU FERMIEGO w metalu (delokalizowanego),")
print("      nie izolowanego atomu. Koherentny fermion sea ma WSPÓLNY fazowy ogon.")

# ---------------------------------------------------------------------------
# WERDYK
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("  AS02 VERDICT")
print("=" * 70)
print("""
  T1 (Li IE chain):      hydrogenic Z_eff - trywialne skalowanie ✓ OK
  T2 (Koide-like):       K_Li = {:.4f} ≠ 2/3. Brannen NIE pasuje.
                         → Li IE NIE JEST drabiną solitonową.
  T3 (A ∝ IE^(1/4)):    A_Li ≈ 0.79, ale A_s(SC) ≈ 0.11. Rząd wielkości niezgodny.
                         → Naiwne TGP przenoszenie A^4 na atomy ZAWODZI.
  T4 (tail-phase sign):  Li 2s tail ujemne ✓ zgodne znakiem A_s.
                         Magnituda niezgodna.
  T5 (alkali consistency): stałość A_s w SC wyklucza hypoteza per-atom tail-sign.
                         → A_s(SC) jest WŁASNOŚCIĄ MNOGOCIAŁOWĄ Fermi morza.

  WNIOSEK as02:
    A_s, A_sp, A_d, A_f z SC paper NIE SĄ A_tail izolowanych atomów.
    Są to EMPIRYCZNE markery Fermi-level stanów w solidach.

    To wyjaśnia pewien observ: w ps41 (d-class closure) Lu pozostało
    nierozwiązane. TGP jako "orbital-class theory" nie ma pierwszych zasad
    do derywacji A_orb z atomowej fizyki.

  Pierwsze zadokumentowane miejsce, gdzie TGP-core NIE LĄCZY się z TGP-stosowanym:
    • core: mass ∝ A_tail⁴ DLA IZOLOWANEGO solitonu
    • applied (SC/ρT): A_orb klas w SIATCE metalu (empiryczne)
    • brakuje: derywacja A_orb z wieloelektronowej fizyki atomowej w metalu

  Kierunek dla as03/04:
    (a) porzucić ambicję "A_s = A_tail(2s Li)" — to nieprawda.
    (b) zamiast tego: testować ZUPEŁNIE inny aspekt, np. polaryzowalność
        atomu α_pol(Li) = 164 a.u. i jej związek z substratem.
    (c) lub: energia spójności (cohesive energy) Li metal — mostek
        atom↔metal, gdzie A_orb faktycznie ma sens.
""".format(K_Li))
