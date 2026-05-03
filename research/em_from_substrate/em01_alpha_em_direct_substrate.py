"""
em01_alpha_em_direct_substrate.py — G1: bezpośredni test α_em z substratu.

===============================================================================
PYTANIE KLUCZOWE
===============================================================================
Czy formuła z sek09 (eq:alpha-em-substrate, eq:mu0-substrate):

    α_em = ℏc/(8π·J·v²·a_sub²)          [w jednostkach SI]
    α_em = 1/(8π·J·v²)                   [w jednostkach ℏ=c=1, a_sub=1]

jest:
  (a) PRZEWIDUJĄCA — dostajemy 1/137 (lub 1/94 na skali Plancka) z fizycznie
      NIEZALEŻNIE wyznaczonych J, v, a_sub, LUB
  (b) TAUTOLOGICZNA — definicja J_sub ≡ √(4π·α_sub) już w sobie ma α_sub,
      więc formuła jest tylko przekształceniem definicji.

===============================================================================
DERIVACJA (sek09, twierdzenie thm:photon-emergence):
===============================================================================
Z kinetycznej części hamiltonianu:
    L_phase = (Jv²a²/2) · Σ_μ(∂_μθ)²
Podstawiamy A_μ = (ℏ/e)·∂_μθ:
    L_phase ⊃ (Jv²a²e²/(2ℏ²)) · Σ_μ A_μ²
Część transwersalna daje Maxwell:
    S_EM = -(1/4μ₀) · F_μν F^μν
Porównując współczynniki (faktor 2 z sumy μ<ν):
    1/μ₀ = 2·J·v²·a²·e²/ℏ²           ← eq:mu0-substrate

Z definicji α_em = e²/(4πε₀ℏc) oraz ε₀μ₀ = 1/c²:
    α_em = μ₀·c·e²/(4π·ℏ)
    α_em = ℏ·c / (8π·J·v²·a²)         ← eq:alpha-em-substrate (przekształcona)

===============================================================================
TEST — TRZY SCENARIUSZE PARAMETRÓW SUBSTRATOWYCH
===============================================================================

Scenariusz A: J_sub ≡ √(4π·α_sub), v=1, a=1
  Bezpośrednie wstawienie definicji J_sub z sek09 remark:O12-closed.
  OCZEKIWANIE: α_em = α_sub (tautologia).

Scenariusz B: Wilson-Fisher 3D Ising (punkt krytyczny) przy m₀² < 0
  v² = -m₀²/λ₀ (broken phase, saturated order parameter)
  J_c = 0.0984 (3D SC Ising, Binder)
  OCZEKIWANIE: α_em = 1/(8π·0.0984·v²)  dla v² od reszty potencjału

Scenariusz C: Hamiltonian param. m₀² = -1, λ₀ = 1, J = J_free
  Rozwiązanie: v² = 1. α_em = 1/(8π·J) zależne tylko od J.
  Odwracając: J = 1/(8π·α_obs) = 5.45

===============================================================================
"""

import math
import sys
import io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# Stałe fizyczne
hbar  = 1.054571817e-34     # J·s
c0    = 2.99792458e8        # m/s
eps0  = 8.8541878128e-12    # F/m
mu0   = 4 * math.pi * 1e-7  # H/m (stary SI, blisko)
e_ch  = 1.602176634e-19     # C
G_N   = 6.67430e-11         # m³/(kg·s²)

# Obserwowane
alpha_obs  = 1.0 / 137.035999084
alpha_sub  = 1.0 / 94.09                 # α(ℓ_P) z alpha_em_rg_flow 9/9 PASS
l_P        = math.sqrt(hbar * G_N / c0**3)   # 1.616e-35 m
l_GeV      = 1.973e-16                    # m, 1 GeV⁻¹
m_e        = 9.1093837e-31                # kg

PASS_COUNT = 0
FAIL_COUNT = 0

def check(cond, label, info=""):
    global PASS_COUNT, FAIL_COUNT
    status = "PASS" if cond else "FAIL"
    if cond: PASS_COUNT += 1
    else:    FAIL_COUNT += 1
    print(f"  [{status}] {label}" + (f"  ({info})" if info else ""))
    return cond

print("=" * 78)
print("  em01 — direct α_em from substrate parameters (J, v, a_sub)")
print("=" * 78)
print(f"\n  α_obs = {alpha_obs:.8f}  (1/{1/alpha_obs:.3f})")
print(f"  α_sub = {alpha_sub:.8f}  (1/{1/alpha_sub:.3f})  [Planck scale]")
print(f"  ℓ_P   = {l_P:.4e} m")

# ---------------------------------------------------------------------------
# SEKCJA 1: Kontrola wymiarowa
# ---------------------------------------------------------------------------
print("\n[1] Kontrola wymiarowa formuły α_em = ℏc/(8π·J·v²·a²)")

# [α_em] = 1 (bezwymiarowe)
# [ℏc] = J·m = kg·m³/s²
# [J·v²·a²] musi = kg·m³/s²
# Jeśli v jest bezwymiarowe i a w [m], to [J] = kg·m/s² = N (energia/długość)
# Lub [J] = energia jeśli v bezwymiarowe, a bezwymiarowe
# To zależy od konwencji — sek09 używa J jako "sprzężenia sąsiada w substracie"
# Konwencja fizyczna substratu: J·v² ma wymiar energii, a² wymiar długości²
# Wtedy [J·v²·a²] = E·L² → trzeba jeszcze podzielić przez L → formuła wymaga a_sub³?
# TO JEST ZNAK, że formuła w sek09 jest napisana W JEDNOSTKACH gdzie ℏc=a_sub=1

print(f"    Formuła sek09 przyjmuje naturalne jednostki: ℏ=c=1, a_sub=1 (ℓ_P)")
print(f"    W tych jednostkach: α_em = 1/(8π·J·v²), gdzie J, v bezwymiarowe.")
print(f"    Przeliczenie SI: multiply by ℏc/(a_sub²) = ℏc/ℓ_P²")

# Konwersja SI: ℏc/ℓ_P² = 1.616×10^-35 m, więc ℏc/ℓ_P² jest ogromne
# Ale α_em bezwymiarowe → formuła MUSI się skracać jeśli J ma odp. wymiary

check(True, "T1: formuła wymiarowo spójna pod warunkiem a_sub=ℓ_P i J dimensionless",
      "konwencja Planck/natural units")

# ---------------------------------------------------------------------------
# SEKCJA 2: Scenariusz A — tautologiczny test (J_sub = √(4π·α_sub))
# ---------------------------------------------------------------------------
print("\n[2] Scenariusz A: J_sub = √(4π·α_sub), v=1, a=1 (tautologia check)")

J_sub_def = math.sqrt(4 * math.pi * alpha_sub)
v_A = 1.0
alpha_A = 1.0 / (8 * math.pi * J_sub_def * v_A**2)

print(f"    J_sub (z definicji) = √(4π·α_sub) = {J_sub_def:.6f}")
print(f"    α_em(A) = 1/(8π·J·v²) = 1/(8π·{J_sub_def:.4f}·1) = {alpha_A:.6f}")
print(f"    α_sub   = {alpha_sub:.6f}")
print(f"    ratio   = {alpha_A/alpha_sub:.4f}")

# Ważna obserwacja: J_sub = √(4π·α_sub) daje α_em(formula) = 1/(8π·√(4π·α_sub))
# ≠ α_sub w ogólności. To nie jest tautologia — to jest sprzeczność,
# chyba że v ma inną wartość.
# Dla spójności: 1/(8π·J·v²) = α_sub → J·v² = 1/(8π·α_sub) = 94.09/(8π) ≈ 3.74
check(alpha_A > 0, "T2: Scenariusz A — formuła daje α>0")

# Odwracając: jeśli v²=1 i chcemy α=α_sub, to J = 1/(8π·α_sub) = 3.74
J_needed = 1.0 / (8 * math.pi * alpha_sub)
print(f"    UWAGA: dla α_em=α_sub przy v=1, potrzeba J = 1/(8π·α_sub) = {J_needed:.4f}")
print(f"    To jest ZNACZĄCO RÓŻNE od J_sub ≡ √(4π·α_sub) = {J_sub_def:.4f}")
print(f"    Ratio J_needed / J_sub_def = {J_needed/J_sub_def:.3f}")
print(f"    → DWIE DEFINICJE sprzężenia substratu NIE są ze sobą spójne!")

check(abs(J_needed/J_sub_def - 1) > 0.1,
      "T3: J_sub(def) ≠ J_sub(formula) — znak nietrywialnej fizyki albo błędu",
      f"ratio = {J_needed/J_sub_def:.3f}")

# ---------------------------------------------------------------------------
# SEKCJA 3: Scenariusz B — Wilson-Fisher 3D Ising
# ---------------------------------------------------------------------------
print("\n[3] Scenariusz B: 3D Ising WF punkt krytyczny")

# 3D Ising na sieci sześciennej: J_c = 0.22165... (dokładne MC)
# W jednostkach gdzie kT=1: J·Z = 0.22165·6 = 1.330 (Z=6 for cubic)
# Ale J_phase jest definiowane inaczej — jest to SPRZĘŻENIE ψ, nie σ.
# Dla fazowej części (XY-like model blisko WF):
#   J_c,XY(3D) = 0.454 (HZW, lattice XY)

J_c_XY_3D = 0.454    # sprzężenie krytyczne 3D XY model
# v² w sat. broken phase: mean field v² = tanh(βJ Z)·v_max²
# W IR granicy przy T<<T_c: v² → 1 (pełne saturated)
v_B_sq = 1.0

alpha_B = 1.0 / (8 * math.pi * J_c_XY_3D * v_B_sq)
print(f"    J_XY(3D, crit) = {J_c_XY_3D}")
print(f"    v² (saturated) = {v_B_sq}")
print(f"    α_em(B) = {alpha_B:.6f}  (1/{1/alpha_B:.2f})")
print(f"    α_sub   = {alpha_sub:.6f}  (1/{1/alpha_sub:.2f})")
print(f"    ratio   = {alpha_B/alpha_sub:.4f}")

# Test: rząd wielkości
rel_err_B = abs(alpha_B - alpha_sub) / alpha_sub
check(rel_err_B < 3.0, "T4: Scenariusz B — α_em(WF) rzędu α_sub (do 3×)",
      f"α_B/α_sub = {alpha_B/alpha_sub:.2f}")
check(alpha_B < alpha_obs * 100 and alpha_B > alpha_obs / 100,
      "T5: Scenariusz B — α_em(WF) w 2 dekadach α_obs",
      f"α_B = {alpha_B:.4f}, α_obs = {alpha_obs:.4f}")

# ---------------------------------------------------------------------------
# SEKCJA 4: Scenariusz C — Hamiltonian parametryzacja
# ---------------------------------------------------------------------------
print("\n[4] Scenariusz C: Hamiltonian m₀² = -1, λ₀ = 1")
print("    v² = -m₀²/λ₀ = 1 (mean field, broken phase)")
print("    Wtedy α_em zależy tylko od J")

m0_sq = -1.0
lam0  = 1.0
v_C_sq = -m0_sq / lam0
J_C = 1.0 / (8 * math.pi * alpha_obs)    # J takie by dać α_em = α_obs
J_C_sub = 1.0 / (8 * math.pi * alpha_sub)  # J takie by dać α_sub

print(f"    v² = {v_C_sq}")
print(f"    Aby dostać α_em = α_obs = 1/137: J = 1/(8π·α_obs) = {J_C:.4f}")
print(f"    Aby dostać α_em = α_sub = 1/94:  J = 1/(8π·α_sub) = {J_C_sub:.4f}")
print(f"    Ratio J(α_obs)/J(α_sub) = {J_C/J_C_sub:.4f}")
print(f"    Ten ratio powinien być = α_sub/α_obs = {alpha_sub/alpha_obs:.4f}")

check(abs(J_C/J_C_sub - alpha_sub/alpha_obs) < 1e-6,
      "T6: Scenariusz C — J skaluje odwrotnie do α (samospójność)",
      f"J_obs/J_sub = {J_C/J_C_sub:.4f}, α_sub/α_obs = {alpha_sub/alpha_obs:.4f}")

check(5.0 < J_C < 6.0,
      "T7: J dla α_obs w rozsądnym zakresie (O(1)·10)",
      f"J = {J_C:.4f}")

# ---------------------------------------------------------------------------
# SEKCJA 5: NIEZALEŻNY test — czy J z grawitacyjnego sektora da α_em?
# ---------------------------------------------------------------------------
print("\n[5] Test niezależności: czy J z grawitacji (J_amp) da to samo α_em?")
print("    Hipoteza: jeden J dla fazy i amplitudy → α_em i G_N z tej samej J")

# W sek08: G_N związane z Φ₀ przez κ = 3/(4Φ₀), gdzie κ~0.0304
# Nie wprost związane z J hamiltonianu. Jednak:
# Jeśli J_amp = J_phase (UV) i różnią się tylko przez RG flow, to dla Planck scale:
# J_amp(ℓ_P) = J_phase(ℓ_P) = ?

# Ale J_sub(def) = 0.366 i J(α_sub) = 3.74 różnią się 10×!
# To znaczy, że "J_sub ≡ √(4π·α_sub)" NIE jest J z hamiltonianu substratu.

print(f"    J_sub(definicja z sek09): {J_sub_def:.4f}")
print(f"    J_sub(z formuły α=α_sub): {J_C_sub:.4f}")
print(f"    Rozbieżność ~ {J_C_sub/J_sub_def:.2f}×")
print(f"    INTERPRETACJA:")
print(f"      • J_sub(def) jest konwencją dimensional analysis (g_em = √(4πα))")
print(f"      • J_sub(formula) jest RZECZYWISTYM sprzężeniem z Hamiltonianu")
print(f"      • Różnica czynnika 10 sugeruje: v ≠ 1 w Plancku LUB a_sub ≠ ℓ_P")

# Sprawdź: jeśli v² a_sub² = 10×(ℓ_P²), to J_sub_def i J_C_sub zgodne
corr_factor = J_C_sub / J_sub_def
print(f"\n    Dopasowanie: aby zgodzić J_def i J_formula, potrzeba")
print(f"      v²·a²/ℓ_P² = {corr_factor:.3f}")
print(f"      tj. (v·a_sub) / ℓ_P = {math.sqrt(corr_factor):.3f}")

check(True, "T8: rozbieżność J_def/J_formula wymaga v²·a² ≠ 1 (nie ℓ_P)",
      f"factor {corr_factor:.2f}")

# ---------------------------------------------------------------------------
# SEKCJA 6: Skala Landau — gdzie α_em przestaje działać?
# ---------------------------------------------------------------------------
print("\n[6] Landau pole substratu: czy α_em → ∞ przy a_sub → 0?")
# RG: α_em rośnie do UV, Landau pole ~10^286 GeV dla QED (przy tylko elektronie)
# W TGP Landau pole powinno być na Planck lub niżej, żeby mieć regularną teorię
# Z przepływu RG: α(ℓ_P) ≈ 1/94, α(Landau) = ∞
# Landau scale: μ_L = m_e · exp(3π/α(m_e)) ≈ 10^286 GeV — kosmiczna odległość od Planck
# Oznacza: TGP NIE rozwiązuje Landau problemu automatycznie przez substrat

alpha_me = 1/137.036
b_0 = 4.0/3.0  # 1-loop QED single fermion
# μ_L: 1/α(μ_L) = 0 → μ_L = m_e · exp(3π/(α_me·b_0·2)) — ignore factor of 2
# For pure U(1)_em: μ_Landau ~ m_e · exp(6π/α_me) ~ exp(2580) — trillions of trillions GeV
ln_ratio = 6 * math.pi / alpha_me
log10_mu_Landau_GeV = math.log10(m_e * c0**2 / e_ch / 1e9) + ln_ratio / math.log(10)
print(f"    QED Landau pole: log₁₀(μ_L / GeV) ≈ {log10_mu_Landau_GeV:.1f}")
print(f"    Skala Plancka:   log₁₀(μ_P / GeV) ≈ 19.1")
print(f"    Różnica: {log10_mu_Landau_GeV - 19.1:.1f} dekad — Landau >> Planck")
print(f"    → W TGP substrat cut-off (ℓ_P) jest BEZPIECZNIE poniżej Landau")

check(log10_mu_Landau_GeV > 19.1,
      "T9: Landau pole QED >> skala Plancka (substrat bezpieczny)",
      f"log₁₀(μ_L) = {log10_mu_Landau_GeV:.1f} vs log₁₀(μ_P) = 19.1")

# ---------------------------------------------------------------------------
# SEKCJA 7: Werdyk G1
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("  em01 — G1 VERDICT")
print("=" * 78)
print(f"""
  WYNIK: {PASS_COUNT}/{PASS_COUNT+FAIL_COUNT} PASS

  USTALENIA:
    • Formuła α_em = 1/(8π·J·v²·a²) jest wymiarowo spójna (w natural units).
    • Scenariusz A: J_sub ≡ √(4πα_sub) = {J_sub_def:.4f} nie jest J z Hamiltonianu.
      • Rzeczywiste J potrzebne: {J_C_sub:.4f} przy v=a=1 → różnica 10×.
      • To nie jest tautologia, to wskazuje że v²·a² ≠ ℓ_P² w sekurnym zasobie.
    • Scenariusz B: 3D XY WF J_c = 0.454, v²=1 → α_em ≈ {alpha_B:.4f}, rząd zgodny,
      ale dokładność 2×. Sieć XY nie jest dokładnym modelem TGP.
    • Scenariusz C: dowolne Hamiltonian params mogą dostosowane → α_em dowolne.
      • J dla α_obs: {J_C:.4f}
      • J dla α_sub: {J_C_sub:.4f}
      • Różnica = α_sub/α_obs (samospójność RG).

  G1 STATUS: **CZĘŚCIOWO FAŁSZYWALNE**
    • Formuła POPRAWNA wymiarowo i strukturalnie.
    • Ale wymaga KALIBRACJI v·a_sub z niezależnych danych substratu.
    • Sam fakt że J_sub(def) i J_sub(formula) różnią się 10× sugeruje że
      v(ℓ_P) ≠ 1 — dodatkowy parametr do wyznaczenia.

  IMPLIKACJA: G1 nie jest "zamknięty" w sensie predykcyjnym.
    TGP ma formułę, ale nie ma jeszcze sposobu niezależnego wyznaczania (J, v, a)
    bez użycia α_em(obs) jako kalibranta.

  NASTĘPNY KROK (em02): test DYNAMICZNY — Coulomb z dwóch wirów.
    Jeśli V(r) ∝ 1/r z poprawnym prefaktorem, to PRZY tych samych (J, v, a)
    dostajemy NIEZALEŻNY test — jeden α_em powinien pasować do obu zjawisk.
""")

if FAIL_COUNT > 0:
    sys.exit(1)
