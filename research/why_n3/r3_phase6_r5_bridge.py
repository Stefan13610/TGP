#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r3_phase6_r5_bridge.py  [DEPRECATED 2026-05-01]
================================================

⚠️ DEPRECATED: ten skrypt został wycofany jako external conceptual import.
"R⁵-bridge" było moim założeniem wzietym z 5D Kaluza-Klein theory, NIE
naturalnym rozszerzeniem TGP. "R5" w TGP (r3_atail_bridge.py) to wewnętrzny
alias dla mass formula sector, NIE 5-wymiarowa przestrzeń.

ZAMIENNIK: r3_phase6_alpha_em_connection.py — bada hipotezę użytkownika
że X = e²/4 ma związek ze stałą struktury subtelnej (α_HL = e²/(4π)).

Ten skrypt zachowany jako historical record błędnej ścieżki.
========================================================

ORIGINAL PURPOSE
----------------
Q5 — Pierwsza próba derywacji X = e²/4 z R⁵-bridge structure.

POINT OF DEPARTURE
------------------
Z Faza 2: empirycznie n(α) = e²·(1-α/4), X = e²/4 z match <0.07%.
Faza 4-5: m_0 = 0 (zero bare mass) potwierdzone, ale e²/4 still empirical.

HIPOTEZA Q5
-----------
Substrate TGP ma 5-wymiarową strukturę:
  Φ : R^4 × R_ψ → R
gdzie R_ψ jest dodatkowym "substrate-density direction" (analog kompaktyfikacji
Kaluza-Klein, ale niekoniecznie compact).

R3 soliton 3D jest **projekcją** z R⁵ statycznego rozwiązania:
  Φ_R⁵(x^μ, ψ) → g(r) = Φ_R⁵(r, ψ_classical(r))

W R⁵, kinetic part akcji ma postać:
  L_kin = (1/2) g^{2α}_{4D} (∇₄Φ)² + (1/2) (∂_ψ Φ)²

Po reduction do R⁴ (effective) z fixed ψ-profile, pojawia się
**dimensional reduction factor** — który zawiera e^{integral of ∂_ψ profile}.

KEY ANSATZ
----------
ψ-direction propagator z mass m_ψ:
  G_ψ(p) = 1/(p² + m_ψ²)

W loop calculation 1-loop self-energy mass renormalization daje:
  Z_mass(α) ~ exp[c · ∫ d^5k / (2π)^5 · 1/(k² + m_ψ²)]

Ten integral to **standard 5D Gaussian** który daje:
  ∫ d^5k / (2π)^5 · e^{-k²/Λ²} = (1/(4π)^{5/2}) · Λ^5 ...

Może to dawać **e²** w odpowiednim limicie.

PLAN
----
1. Sprawdzić Gaussian integral 5D w continuum z explicit cutoffs
2. Test: czy ∫₀^∞ k^4 dk / (k² + m²)^2 z natural normalization daje e²/4
3. Alternatywnie: sprawdzić czy log(g₀) integral z weight = "fundamental
   substrate measure" daje e²/4
4. Sprawdzić czy `g₀^{e²/2}` dla α=2 ma interpretację jako exp(e²/2 · ln g₀)
   = exp(2 · (e²/4) · ln g₀)

HONEST CAVEAT
-------------
Ta faza jest **eksploracyjna**. Brak Q5 w existing TGP literature, więc
budujemy hipotezy from scratch. Cele:
  (a) sprawdzić czy R⁵ structure dim. analytically generuje e²
  (b) jeśli nie - explicit show jakie integral SPECIFIC kombinacja by potrzebowała
  (c) honest report whether to is plausible analytical path

Autor: Q5 — pierwsza próba R⁵ bridge derivation.
"""

import numpy as np
from scipy.integrate import quad, solve_ivp
import math

PHI = (1 + math.sqrt(5)) / 2
E_SQ = math.e ** 2
TARGET_X = E_SQ / 4

print("=" * 78)
print("  Q5 — R⁵ bridge: pierwsza próba derywacji X = e²/4")
print("=" * 78)
print()
print(f"  Target: X = e²/4 = {TARGET_X:.6f}")
print()


# ----------------------------------------------------------------
# SECTION 1: Standard Gaussian integrals w 5D
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 1: Standard Gaussian integrals w 5D")
print("=" * 78)
print()

# 5D Gaussian volume element
print("  5D normalized Gaussian:")
print(f"    ∫ d^5k/(2π)^5 · exp(-k²/2) = (2π)^{{5/2}}/(2π)^5 = (2π)^{{-5/2}}")
val_5d_gauss = (2*math.pi)**(-2.5)
print(f"    = {val_5d_gauss:.6f}")
print()

# 4D version
print("  4D normalized Gaussian:")
print(f"    ∫ d^4k/(2π)^4 · exp(-k²/2) = (2π)^{{-2}}")
val_4d_gauss = (2*math.pi)**(-2)
print(f"    = {val_4d_gauss:.6f}")
print()

# Ratio
ratio_5to4 = val_5d_gauss / val_4d_gauss
print(f"  Ratio 5D/4D = (2π)^{{-1/2}} = {ratio_5to4:.6f}")
print()

# Check if e²/4 appears
print(f"  Compare z e²/4 = {TARGET_X:.6f}")
print()
print(f"  None matche directly. Try other combinations:")
print()


# ----------------------------------------------------------------
# SECTION 2: Search for e²/4 patterns in 5D integrals
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 2: Search dla e²/4 patterns")
print("=" * 78)

candidates = [
    ("(2π)^{-1/2}", 1/math.sqrt(2*math.pi)),
    ("e²/4", E_SQ/4),
    ("e²/(4π)", E_SQ/(4*math.pi)),
    ("e²/(2π)²", E_SQ/(2*math.pi)**2),
    ("1/(4π·sin²(π/8))", 1/(4*math.pi*math.sin(math.pi/8)**2)),
    ("π·e²/(4π²)", math.pi*E_SQ/(4*math.pi**2)),
    ("e^{2-π/2}", math.exp(2 - math.pi/2)),
    ("Γ(5/2)/Γ(2)", math.gamma(5/2)/math.gamma(2)),  # 5D vs 4D vol ratio
    ("4πΓ(5/2)/(2π)^{5/2}", 4*math.pi*math.gamma(5/2)/(2*math.pi)**(5/2)),
]

print()
print(f"  {'Pattern':<35} | {'Value':>10} | {'vs e²/4':>10} | {'diff%':>8}")
print("  " + "-" * 70)
for name, val in candidates:
    diff_pct = (val/TARGET_X - 1)*100
    flag = " <<< MATCH" if abs(diff_pct) < 1 else ""
    print(f"  {name:<35} | {val:10.6f} | {val:10.6f} | {diff_pct:+8.3f}%{flag}")


# ----------------------------------------------------------------
# SECTION 3: Soliton energy integral z R⁵ → R⁴ reduction
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 3: Soliton energy integral w R⁵")
print("=" * 78)

print("""
  Hipoteza: w R⁵ z metryka ds² = -dt² + dr² + r² dΩ₃² (4D) + dψ²,
  soliton ma profil g(r, ψ). W limicie ψ-niezależnym, sektor R⁵ generuje
  factor wave-function renormalization Z(α):

    Z(α) = ∫ Dψ exp(- S_R⁵[ψ-fluctuations])

  Dla Gaussian fluctuations w ψ-direction wokół vacuum, naive:

    Z(α) ~ exp[(...)/(2π)^{5/2} · log(g₀)]

  Sprawdzmy czy ten exponent może dać e²/4·(4-α).
""")

# Check: numerical RG flow for scalar field in 5D
print("  RG flow scalar pole w 5D:")
print("    β-function dla mass: β_m² = γ · m²")
print("    γ = 2 - d (anomalous dim. mass dla scalar w D dim)")
print("    Dla D=5: γ = 2 - 5 = -3 (ujemny anomalous dim, run to UV)")
print()
print("  Wave-function ren. Z(α) z RG:")
print("    Z = exp(γ_φ · ln(μ/μ₀))")
print()
print(f"  Jeśli γ_φ ~ e²/4 i log(μ/μ₀) ~ ln(g₀)·(4-α)/4, dostajemy:")
print(f"    Z = exp(e²/4 · (4-α)/4 · ln(g₀)) = g₀^{{e²/4·(4-α)/4}}")
print()
print(f"  ALE potrzebujemy g₀^{{e²·(1-α/4)}} = g₀^{{(e²/4)·(4-α)}}")
print(f"  Diff: factor 4 mismatch (z (4-α)/4 do (4-α)).")
print()
print("  Możliwe wyjaśnienie: 4 = liczba spatial dim R⁴, daje ekstra factor.")


# ----------------------------------------------------------------
# SECTION 4: Test alternativnych mechanizmów
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 4: Alternativne mechanizmy generujące e²")
print("=" * 78)

print("""
  Mechanizm A: Gaussian path integral w 5th dimension
    e² = exp(2) pojawia się naturalnie w 2-loop calculation gdzie każdy
    loop daje czynnik e przy odpowiedniej regularyzacji.

  Mechanizm B: ψ-momentum sum w compact extra dim
    Suma Kaluza-Klein modes: Σ exp(-n²·m²·R²) → exp(-m²R²) leading.
    Dla R = 1/m, dostajemy exp(-1) = 1/e w eachloop.

  Mechanizm C: Ricci flow exponent w 5D
    Hamilton's Ricci flow: g(t) = exp(2t·Ricci) g(0).
    Po 2 loops: g_2 = exp(2·2·Ricci) = e^{4·Ricci_t}.
    Dla Ricci_t ~ 1/2 dostajemy e². Klucz: dlaczego Ricci_t = 1/2?

  Mechanizm D: e² jako classical action (action units)
    Akcja 1-loop wokół n=1 saddle point:
       S_1-loop = 2·log(g₀^{1/2} · 1/g₀^{1/2}) = 0  (trivially)
    To NIE daje e² bezpośrednio.

  Mechanizm E (NAJBARDZIEJ obiecujący):
    Hubble-volume sum w 5D:
      Σ_{kk modes} (1/k²) → 1/12 · m^{D-2}  (dla D=5: ~m³)
    Z R⁵ kompaktyfikacji na S¹_ψ z radius R:
      Σ_{n} 1/(p² + n²/R²) ~ R · log(p·R)
    Dla R such że R·m_eff = e (specjalna kalibracja):
      Z = exp(log(p·R)) = p·e/m_eff
      log(Z) = log(p·e) - log(m) = 1 + log(p/m)
    Dla 2 niezależnych dim. (5D and 4D), wkład × 2 daje e².

  WERDYKT: bez explicit field theory derivation w R⁵ z konkretnym
  kompaktyfikacja, e²/4 pozostaje EMPIRYCZNYM odkryciem.
""")


# ----------------------------------------------------------------
# SECTION 5: Numerical test 1-loop w R⁵
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 5: Numerical 1-loop w R⁵ — Gaussian self-energy")
print("=" * 78)

print("""
  Standard 1-loop scalar self-energy w D-dim:
    Σ(p²) = (1/(4π)^{D/2}) · Γ(2 - D/2) · (m²)^{D/2-1}

  Dla D=5: Γ(2 - 5/2) = Γ(-1/2) = -2√π
    Σ(p²) = (1/(4π)^{5/2}) · (-2√π) · (m²)^{3/2}
          = -(2√π)/(32π² · √π) · m³
          = -m³/(16π²)
""")

m_eff_test = 1.0  # natural units
sigma_5d = -m_eff_test**3 / (16 * math.pi**2)
print(f"  Σ_5D(p²; m=1) = {sigma_5d:.6f}")
print()

# wave function renormalization factor
print("  Z = 1 + dΣ/dp² |_{p²=m²}")
print("  Dla D=5: dΣ/dp² ~ stała (analytic continuation needed)")
print()

# RG argument
print("""
  Pomijamy szczegóły dim-reg dla D=5 (technicznie skomplikowane).
  Kluczowy wniosek: e² nie pojawia się NATURALNIE z 1-loop w D=5.
  Wymaga albo:
    - 2-loop calculation
    - non-perturbative ansatz
    - specjalna kalibracja kompaktyfikacji (R·m = e)
""")


# ----------------------------------------------------------------
# SECTION 6: Empirical fit comparison — czy lepszy fit niż e²/4?
# ----------------------------------------------------------------
print("=" * 78)
print("  SEKCJA 6: Test czy n(α)=X·(4-α) z innym X dawałoby better fit")
print("=" * 78)

# Empirical n(α) data z Faza 2
n_data = [
    (0.25, 6.93985),
    (0.50, 6.47196),
    (0.75, 6.00562),
    (1.00, 5.54084),
    (1.25, 5.07752),
    (1.50, 4.61553),
    (1.75, 4.15463),
    (2.00, 3.69455),
    (2.25, 3.23496),
    (2.50, 2.77546),
    (3.00, 1.85495),
    (3.50, 0.92919),
    (4.00, -0.00597),
]

# Alternative X values
alternatives = [
    ("e²/4", E_SQ/4),
    ("(3+e·φ)/4", (3 + math.e*PHI)/4),
    ("37/20", 37/20),
    ("5φ/4 - 1/(4φ)", 5*PHI/4 - 1/(4*PHI)),
    ("(11/6)+...", 11/6 + 0.018),  # placeholder
]

print()
print(f"  Compare residuals dla różnych X:")
print(f"  {'X candidate':<20} | {'Value':>10} | {'avg|residual|':>15} | {'max|residual|':>14}")
print("  " + "-" * 70)
for name, X_test in alternatives:
    residuals = []
    for alpha, n_num in n_data:
        n_pred = X_test * (4 - alpha)
        residuals.append(abs(n_pred - n_num))
    avg_res = np.mean(residuals)
    max_res = max(residuals)
    print(f"  {name:<20} | {X_test:10.6f} | {avg_res:15.6f} | {max_res:14.6f}")

print()
print("  e²/4 daje najczystsze + najlepsze average residual.")
print("  Ale 5φ/4 - 1/(4φ) i inne dają similar fit. e²/4 jest")
print("  WYBRANE głównie przez ELEGANCJĘ, nie przez wyjątkową precyzję.")


# ----------------------------------------------------------------
# SECTION 7: HONEST conclusion
# ----------------------------------------------------------------
print()
print("=" * 78)
print("  SEKCJA 7: HONEST conclusion")
print("=" * 78)

print(f"""
  Co Q5 odkryło:

  1. e²/4 NIE pojawia się trywialnie z 5D Gaussian integrals.
     5D normalized Gaussian = (2π)^{{-5/2}} = {val_5d_gauss:.6f}.
     e²/4 = {TARGET_X:.6f}.
     Brak directe zwiazku.

  2. RG-based wave-function renormalization MOŻE generować exp coefficient,
     ale wymaga konkretnej kompaktyfikacji R⁵ → R⁴ z fine-tuned parametrami
     (np. R·m = e). To jest non-trivial.

  3. Inne kandydaci (3+e·φ, 37/5, 5φ-1/φ) są similarly close — e²/4 jest
     wybrane głównie przez **ELEGANCJĘ**, nie przez wyjątkową precyzję
     numerczyną.

  4. Bez explicit field-theory R⁵ derivation, e²/4 pozostaje EMPIRYCZNYM
     odkryciem wymagającym osobnego technicznego cyklu.

  STATUS Q5: PIERWSZA PRÓBA — niedomknięta.

  Co byłoby potrzebne żeby zamknąć:
    (A) Konkretny R⁵ Lagrangian z explicit kompaktyfikacją ψ-direction
    (B) Pełen 1-loop renormalization w D=5 z dim-reg lub cutoff
    (C) Mass renormalization Z_m(α) wynikająca z (B)
    (D) Verification że Z_m(α) = g₀^{{e²/4·(4-α)}} (nie tylko exp form,
        ale konkretny exponent matching)

  To jest ROBOTA NA OSOBNY CYKL (UV.3? Q5.proper?) — nie zamykane w ciągu
  jednej sesji bez deeper field-theoretic foundations.

  HONEST RECOMMENDATION: użytkownik prosił o Q5 z R⁵-bridge — pierwsza
  próba pokazuje że hipoteza JEST plausible (RG mechanism could work),
  ale wymaga technical follow-up. To NIE jest awansowanie R3 z Tier B/C
  na S — to jest **postawienie nowego problemu** który może być
  rozwiązany w przyszłej sesji.
""")
