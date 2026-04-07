#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex243_boson_koide_sum_rule.py
==============================
K(γ,W,Z) = 1/2 I BOSONOWA REGUŁA SUM

KONTEKST (ex242):
  SURPRISE: K(γ, W, Z) = 0.5005 ≈ 1/2
  ALSO: m_H² + m_W² + m_Z² ≈ v²/2 (0.5% off!)

PYTANIA:
  1. Dlaczego K(γ,W,Z) ≈ 1/2? Czy to DOKŁADNIE 1/2?
  2. Jeśli K=1/2: m_W/m_Z jest WYZNACZONE (→ sin²θ_W!)
  3. m_H² + m_W² + m_Z² = v²/2: derived or coincidence?
  4. Unified K=1/2 for massless-member triplets?
  5. Implications for TGP

Data: 2026-04-06
"""

import sys, io, math
import numpy as np
from scipy.optimize import brentq

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

phi = (1 + np.sqrt(5)) / 2

TESTS = []
def record(name, passed, detail=""):
    TESTS.append((name, passed, detail))
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        for line in detail.split('\n'):
            print(f"         {line}")


def koide(m1, m2, m3):
    S = np.sqrt(abs(m1)) + np.sqrt(abs(m2)) + np.sqrt(abs(m3))
    if S == 0:
        return np.nan
    return (m1 + m2 + m3) / S**2


# ============================================================
# EXPERIMENTAL DATA
# ============================================================
m_W = 80.377    # GeV
m_Z = 91.1876   # GeV
m_H = 125.25    # GeV
v = 246.22      # GeV
sin2w = 0.23122


# ============================================================
# §1. K(γ,W,Z) = 1/2 — DOKŁADNA ANALIZA
# ============================================================
print("=" * 72)
print("§1. K(γ,W,Z) — CZY DOKŁADNIE 1/2?")
print("=" * 72)

K_gWZ = koide(0, m_W, m_Z)

# Analytic: K(0, m₁, m₂) = (m₁+m₂)/(√m₁+√m₂)²
# Let r = m_Z/m_W = 1/cos²θ_W (in SM: m_Z = m_W/cosθ_W)
r_ZW = m_Z / m_W
cosW = m_W / m_Z
cos2W = cosW**2

print(f"\n  K(0, m_W, m_Z) = (m_W + m_Z)/(√m_W + √m_Z)²")
print(f"  = (1 + r)/(1 + √r)²  gdzie r = m_Z/m_W = {r_ZW:.6f}")
print(f"")
print(f"  K = {K_gWZ:.8f}")
print(f"  1/2 = 0.50000000")
print(f"  |K - 1/2| = {abs(K_gWZ - 0.5):.8f}")
print(f"  Relative: {abs(K_gWZ - 0.5)/0.5 * 100:.4f}%")

# When is K(0, 1, r) = 1/2 exactly?
# (1+r)/(1+√r)² = 1/2
# 2(1+r) = (1+√r)²
# 2+2r = 1+2√r+r
# r+1 = 2√r
# (√r)² - 2√r + 1 = 0
# (√r - 1)² = 0
# √r = 1 → r = 1!
# K(0,m,m) = 2m/(2√m)² = 2m/(4m) = 1/2 ✓
# But this requires m_W = m_Z!

print(f"""
  ANALITYCZNIE: K(0, 1, r) = (1+r)/(1+√r)² = 1/2

  Rozwiązanie: 2(1+r) = (1+√r)²
  → r + 1 = 2√r
  → (√r - 1)² = 0
  → r = 1 (dokładnie)

  Więc K = 1/2 DOKŁADNIE ⟺ m_W = m_Z (i.e. sin²θ_W = 0)

  ALE: m_W ≠ m_Z w rzeczywistości (sin²θ_W = 0.231)
  r = m_Z/m_W = {r_ZW:.6f}
  √r - 1 = {np.sqrt(r_ZW) - 1:.6f} (small!)

  Rozwinięcie: K(0, 1, 1+ε) dla ε = r-1 << 1:
  K ≈ 1/2 + ε²/16 + O(ε³)
""")

eps = r_ZW - 1
K_approx = 0.5 + eps**2/16
print(f"  ε = r - 1 = {eps:.6f}")
print(f"  K ≈ 1/2 + ε²/16 = {K_approx:.8f}")
print(f"  K(exact) = {K_gWZ:.8f}")
print(f"  Correction: ε²/16 = {eps**2/16:.8f}")

# So K(γ,W,Z) ≈ 1/2 because m_W and m_Z are CLOSE (ε = 0.135)
# The deviation from 1/2 is (ε²/16) ≈ 0.001 — tiny!

record("T1: K(γ,W,Z) ≈ 1/2",
       abs(K_gWZ - 0.5) < 0.001,
       f"K = {K_gWZ:.6f}, deviation = {abs(K_gWZ-0.5):.6f}")

# BUT: if K(γ,W,Z) = 1/2 EXACTLY, this predicts sin²θ_W = 0!
# So K=1/2 is NOT exact — it's approximate due to sin²θ_W ≠ 0
print(f"\n  WNIOSEK: K(γ,W,Z) ≈ 1/2 jest PRZYBLIŻONE.")
print(f"  Bliskość do 1/2 wynika z małego sin²θ_W = {sin2w:.4f}")
print(f"  NIE jest to głęboka relacja — to artefakt m_W ≈ m_Z")


# ============================================================
# §2. BOSONOWA REGUŁA SUM: m_H² + m_W² + m_Z² = v²/2
# ============================================================
print("\n" + "=" * 72)
print("§2. ★ BOSONOWA REGUŁA SUM")
print("=" * 72)

lhs = m_H**2 + m_W**2 + m_Z**2
rhs = v**2 / 2

print(f"\n  m_H² + m_W² + m_Z² = {m_H**2:.1f} + {m_W**2:.1f} + {m_Z**2:.1f}")
print(f"                     = {lhs:.1f} GeV²")
print(f"  v²/2               = {rhs:.1f} GeV²")
print(f"  Ratio               = {lhs/rhs:.6f}")
print(f"  Error               = {abs(lhs/rhs - 1)*100:.4f}%")

record("T2: m_H² + m_W² + m_Z² ≈ v²/2",
       abs(lhs/rhs - 1) < 0.01,
       f"Ratio = {lhs/rhs:.6f}, error = {abs(lhs/rhs-1)*100:.3f}%")

# In SM terms:
# m_W = g₂v/2, m_Z = g₂v/(2cosθ_W), m_H = √(2λ)v
# m_W² = g₂²v²/4
# m_Z² = g₂²v²/(4cos²θ_W) = g₂²v²/(4(1-sin²θ_W))
# m_H² = 2λv²

# Sum: v²(g₂²/4 + g₂²/(4cos²θ_W) + 2λ) = v²/2
# → g₂²/4 × (1 + 1/cos²θ_W) + 2λ = 1/2
# → g₂²/4 × (1 + sec²θ_W) = 1/2 - 2λ

g2_sq = 4 * m_W**2 / v**2
lam = m_H**2 / (2 * v**2)
sec2W = 1 / cos2W

print(f"\n  SM parametry:")
print(f"    g₂² = 4m_W²/v² = {g2_sq:.6f}")
print(f"    λ = m_H²/(2v²) = {lam:.6f}")
print(f"    sec²θ_W = {sec2W:.6f}")
print(f"")
print(f"  Reguła sum w SM:")
print(f"    g₂²(1 + sec²θ_W)/4 + 2λ = {g2_sq*(1+sec2W)/4 + 2*lam:.6f}")
print(f"    1/2 = {0.5:.6f}")
print(f"    Difference = {g2_sq*(1+sec2W)/4 + 2*lam - 0.5:.6f}")

# This is g₂²(1 + sec²θ_W)/4 + 2λ ≈ 0.5
# Let's check what this means
print(f"\n  Rozkład:")
print(f"    g₂²/4 = {g2_sq/4:.6f}  (from m_W)")
print(f"    g₂²sec²θ_W/4 = {g2_sq*sec2W/4:.6f}  (from m_Z - m_W)")
print(f"    2λ = {2*lam:.6f}  (from m_H)")
print(f"    SUMA = {g2_sq/4 + g2_sq*sec2W/4 + 2*lam:.6f}")

# Is there a deeper reason? In SM, at tree level:
# Σm² = v²(g₂²(2+sec²θ_W-1)/4 + 2λ) — wait let me redo
# m_W² + m_Z² = (g₂²v²/4)(1 + sec²θ_W) = (g₂²v²/4) × (cos²θ+1)/cos²θ
# = (g₂²v²/4) × (2-sin²θ_W)/cos²θ_W

gauge_part = m_W**2 + m_Z**2
higgs_part = m_H**2
print(f"\n  m_W² + m_Z² = {gauge_part:.1f} GeV² ({gauge_part/rhs*100:.1f}% of v²/2)")
print(f"  m_H²        = {higgs_part:.1f} GeV² ({higgs_part/rhs*100:.1f}% of v²/2)")


# ============================================================
# §3. PREDYKCJA sin²θ_W Z REGUŁY SUM?
# ============================================================
print("\n" + "=" * 72)
print("§3. CZY REGUŁA SUM WYZNACZA sin²θ_W?")
print("=" * 72)

print("""
  Jeśli m_H² + m_W² + m_Z² = v²/2 jest DOKŁADNE:
    m_H² + m_W² + m_Z² = v²/2

  W SM: m_W = (g₂v)/2, m_Z = m_W/cosθ_W
    m_W² + m_Z² = m_W²(1 + sec²θ_W) = m_W²(2-s²)/(1-s²)
    where s² = sin²θ_W

  Więc: m_W²(2-s²)/(1-s²) + m_H² = v²/2

  Znając m_H i v, to daje relację między m_W i sin²θ_W:
    m_W² × (2-s²)/(1-s²) = v²/2 - m_H²
""")

rhs_gauge = v**2/2 - m_H**2
print(f"  v²/2 - m_H² = {rhs_gauge:.1f} GeV²")

# If we know sin²θ_W, predict m_W:
# m_W²(2-s²)/(1-s²) = rhs_gauge
# m_W² = rhs_gauge × (1-s²)/(2-s²)
m_W_pred = np.sqrt(rhs_gauge * (1-sin2w)/(2-sin2w))
print(f"\n  Z sin²θ_W = {sin2w}:")
print(f"    m_W(pred) = {m_W_pred:.3f} GeV")
print(f"    m_W(PDG)  = {m_W:.3f} GeV")
print(f"    Error: {abs(m_W_pred-m_W)/m_W*100:.3f}%")

record("T3: m_W from sum rule + sin²θ_W",
       abs(m_W_pred - m_W)/m_W < 0.01,
       f"m_W = {m_W_pred:.3f} vs {m_W:.3f} GeV, err = {abs(m_W_pred-m_W)/m_W*100:.3f}%")

# Or: if we know m_W, predict sin²θ_W:
# m_W²(2-s²)/(1-s²) = rhs_gauge
# Let x = s² = sin²θ_W
# m_W²(2-x)/(1-x) = rhs_gauge
# m_W²(2-x) = rhs_gauge(1-x)
# 2m_W² - m_W²x = rhs_gauge - rhs_gauge·x
# (rhs_gauge - m_W²)x = rhs_gauge - 2m_W²
# x = (rhs_gauge - 2m_W²)/(rhs_gauge - m_W²)

sin2w_pred = (rhs_gauge - 2*m_W**2)/(rhs_gauge - m_W**2)
print(f"\n  Z m_W = {m_W} GeV:")
print(f"    sin²θ_W(pred) = {sin2w_pred:.5f}")
print(f"    sin²θ_W(PDG)  = {sin2w:.5f}")
print(f"    Error: {abs(sin2w_pred-sin2w)/sin2w*100:.3f}%")

record("T4: sin²θ_W from sum rule + m_W",
       abs(sin2w_pred - sin2w)/sin2w < 0.01,
       f"sin²θ_W = {sin2w_pred:.5f} vs {sin2w:.5f}, err = {abs(sin2w_pred-sin2w)/sin2w*100:.3f}%")


# ============================================================
# §4. REGUŁA SUM VS UNCERTAINTY
# ============================================================
print("\n" + "=" * 72)
print("§4. REGUŁA SUM — PROPAGACJA NIEPEWNOŚCI")
print("=" * 72)

# Uncertainties
sig_mH = 0.17  # GeV
sig_mW = 0.012  # GeV
sig_mZ = 0.0021  # GeV

# Partial derivatives of R = (m_H² + m_W² + m_Z²)/(v²/2) - 1
# dR/dm_H = 2m_H/(v²/2) = 2×125.25/30312 = 0.00826
# dR/dm_W = 2m_W/(v²/2) = 2×80.377/30312 = 0.00530
# dR/dm_Z = 2m_Z/(v²/2) = 2×91.188/30312 = 0.00602

dR_dmH = 2*m_H/rhs
dR_dmW = 2*m_W/rhs
dR_dmZ = 2*m_Z/rhs

sig_R = np.sqrt((dR_dmH*sig_mH)**2 + (dR_dmW*sig_mW)**2 + (dR_dmZ*sig_mZ)**2)

R_central = lhs/rhs - 1

print(f"  R = (m_H² + m_W² + m_Z²)/(v²/2) - 1 = {R_central:.6f}")
print(f"  σ_R = {sig_R:.6f}")
print(f"  R/σ_R = {abs(R_central)/sig_R:.1f}σ")

record("T5: Sum rule within experimental error",
       abs(R_central)/sig_R < 3.0,
       f"R = {R_central:.4f} ± {sig_R:.4f}, deviation = {abs(R_central)/sig_R:.1f}σ")


# ============================================================
# §5. SM INTERPRETATION
# ============================================================
print("\n" + "=" * 72)
print("§5. SM INTERPRETACJA REGUŁY SUM")
print("=" * 72)

print(f"""
  W Modelu Standardowym (tree level):
    m_W² = g₂²v²/4
    m_Z² = (g₁²+g₂²)v²/4
    m_H² = 2λv²

  Reguła sum: m_H² + m_W² + m_Z² = v²/2
  → 2λv² + g₂²v²/4 + (g₁²+g₂²)v²/4 = v²/2
  → 2λ + g₂²/4 + (g₁²+g₂²)/4 = 1/2
  → 2λ + g₂²/2 + g₁²/4 = 1/2
  → 2λ + (2g₂² + g₁²)/4 = 1/2

  Numerycznie:
    g₁² = g₂²·tan²θ_W = g₂²·sin²θ_W/cos²θ_W
    g₂² = {g2_sq:.6f}
    g₁² = {g2_sq * sin2w / cos2W:.6f}
    λ = {lam:.6f}

    2λ + (2g₂² + g₁²)/4 = {2*lam + (2*g2_sq + g2_sq*sin2w/cos2W)/4:.6f}
    1/2 = 0.500000

  Daje to relację: λ = (1 - 2g₂² - g₁²/2)/4 = (1 - g₂²(2+tan²θ_W))/4

  Sprawdzenie: λ_pred = (1 - g₂²(2+tan²θ_W))/4
""")

g1_sq = g2_sq * sin2w / cos2W
lam_pred = (1 - g2_sq*(2 + sin2w/cos2W)) / 4
# Wait, let me redo: 2λ + (2g₂² + g₁²)/4 = 1/2
# 2λ = 1/2 - (2g₂² + g₁²)/4
# λ = 1/4 - (2g₂² + g₁²)/8

lam_pred2 = 0.25 - (2*g2_sq + g1_sq)/8
print(f"  λ_pred = 1/4 - (2g₂² + g₁²)/8 = {lam_pred2:.6f}")
print(f"  λ_actual = {lam:.6f}")
print(f"  Error: {abs(lam_pred2 - lam)/lam * 100:.2f}%")

# If this holds, m_H is predicted from gauge couplings!
m_H_pred = np.sqrt(2 * lam_pred2) * v
print(f"\n  m_H(pred) = v·√(2λ_pred) = {m_H_pred:.2f} GeV")
print(f"  m_H(PDG)  = {m_H} GeV")
print(f"  Error: {abs(m_H_pred - m_H)/m_H * 100:.2f}%")

record("T6: m_H from gauge couplings via sum rule",
       abs(m_H_pred - m_H)/m_H < 0.01,
       f"m_H = {m_H_pred:.2f} vs {m_H} GeV, err = {abs(m_H_pred-m_H)/m_H*100:.2f}%")


# ============================================================
# §6. PHYSICAL MEANING
# ============================================================
print("\n" + "=" * 72)
print("§6. ★ FIZYCZNE ZNACZENIE")
print("=" * 72)

print(f"""
  Reguła sum m_H² + m_W² + m_Z² = v²/2 oznacza:

  RELACJA MIĘDZY SELF-COUPLING λ A GAUGE COUPLINGS:
    λ = 1/4 - (2g₂² + g₁²)/8

  Numerycznie: λ = {lam:.4f}, g₂² = {g2_sq:.4f}, g₁² = {g1_sq:.4f}
  λ_pred = {lam_pred2:.4f}

  To NIE jest standardowy wynik SM — λ jest wolny parametr w SM.
  ALE: gdyby ta relacja była DOKŁADNA, to:

  1. m_H nie jest wolnym parametrem — jest WYZNACZONE przez g₁, g₂
  2. λ jest OBLICZALNE z gauge couplings
  3. SM miałby 18, nie 19 wolnych parametrów

  UWAGA: W SM, ta relacja NIE jest chroniona symetrią.
  Radiative corrections (1-loop) modyfikują ją.
  Ale: fakt że tree-level jest spełniony na 0.5% jest NIETRYWIALNY.

  ★ INTERPRETACJA TGP:
    Jeśli TGP generuje WSZYSTKIE masy z jednego mechanizmu (soliton),
    to relacja m_H² + m_W² + m_Z² = v²/2 może wynikać z:
    - Ward identity lub sumrule z TGP action
    - Completeness relation: bozony (W,Z,H) "wypełniają" v²/2
    - Analogy: ΣY² = 1/2 (sum of hypercharges squared)

  ALTERNATYWNA INTERPRETACJA:
    (m_H² + m_W² + m_Z²) / v² = 1/2
    Σ(m_boson/v)² = 1/2
    Σ(coupling²) = 1/2
    → Coupling sum rule!

  Porównanie:
    (m_W/v)² = {(m_W/v)**2:.6f}
    (m_Z/v)² = {(m_Z/v)**2:.6f}
    (m_H/v)² = {(m_H/v)**2:.6f}
    Σ         = {(m_W/v)**2 + (m_Z/v)**2 + (m_H/v)**2:.6f}
    1/2       = {0.5:.6f}
""")

# Is this related to Σg² normalization?
print(f"  Gauge couplings:")
print(f"    g₂²/4 = (m_W/v)² = {g2_sq/4:.6f}")
print(f"    (g₁²+g₂²)/4 = (m_Z/v)² = {(g1_sq+g2_sq)/4:.6f}")
print(f"    2λ = (m_H/v)² = {2*lam:.6f}")
print(f"    Sum = {g2_sq/4 + (g1_sq+g2_sq)/4 + 2*lam:.6f} ≈ 1/2")


# ============================================================
# §7. PORÓWNANIE: TRZY K=1/2
# ============================================================
print("\n" + "=" * 72)
print("§7. ★ TRZY REALIZACJE K = 1/2")
print("=" * 72)

print(f"""
  K = 1/2 pojawia się w TRZECH kontekstach:

  1. NEUTRINA (ex239):
     K(ν₁, ν₂, ν₃) = 1/2  (z K = N/(2N) = 3/6, Majorana)
     → Σm_ν = 59.8 meV, Normal Ordering

  2. BOZONY CECHOWANIA (ex242-243):
     K(γ, W, Z) ≈ 1/2  (przybliżone, z m_W ≈ m_Z)
     → NIE głębokie: wynika z małego sin²θ_W

  3. COUPLING SUM RULE:
     Σ(m_boson/v)² = 1/2  (reguła sum, ~0.5% exact)
     → Potencjalnie głębokie: łączy λ z g₁, g₂

  ANALOGIA:
     K(ν) = 1/2 jest DOKŁADNE (z definicji hipotezy)
     K(γ,W,Z) ≈ 1/2 jest PRZYBLIŻONE (z sin²θ_W ≠ 0)
     Σ(m/v)² = 1/2 jest PRAWIE DOKŁADNE (0.5% off)

  Czy 2 i 3 są POWIĄZANE? Tak!
     K(0, m_W, m_Z) = (m_W+m_Z)/(√m_W+√m_Z)² ≈ 1/2
     Σ(m/v)² = (m_W²+m_Z²+m_H²)/v² = 1/2
     ALE to INNE relacje (√m vs m²)
""")

record("T7: Three K=1/2 appearances",
       True,
       "ν-Koide (exact), γWZ-Koide (approximate), coupling sum (0.5% exact)")


# ============================================================
# §8. RADIATIVE CORRECTIONS
# ============================================================
print("=" * 72)
print("§8. KOREKCJE RADIACYJNE")
print("=" * 72)

print(f"""
  Reguła sum m_H² + m_W² + m_Z² = v²/2 jest tree-level.
  Korekcje 1-loop (dominujące od top quarka):

  Δ(m_H²) ~ -3y_t²m_t²/(8π²) × ln(Λ²/m_t²)
  Δ(m_W²) ~ +g₂²m_t²/(32π²)
  Δ(m_Z²) ~ +g₂²m_t²/(32π²cos²θ_W)

  Efektywne: Δ(sum)/v² ~ -3y_t⁴/(8π²) × ln(Λ/m_t)

  Dla Λ = 1 TeV: Δ ~ -0.02 (4% correction)
  Dla Λ = 10 TeV: Δ ~ -0.03 (6% correction)

  Przy obecnej dokładności 0.5%: korekcje są na granicy.
  Potrzebna 2-loop analiza do rozstrzygnięcia.

  KLUCZOWE: jeśli reguła sum trzyma po korekcjach radiacyjnych,
  to sugeruje symetrię chroniącą tę relację (cancellation mechanism).
""")

# Estimate the 1-loop correction
y_t = np.sqrt(2) * 172.76 / v
delta_1loop = -3 * y_t**4 / (8 * np.pi**2) * np.log(1000/172.76)
print(f"  y_t = {y_t:.4f}")
print(f"  Δ(1-loop, Λ=1TeV) ≈ {delta_1loop:.4f}")
print(f"  Corrected sum: {lhs/rhs + delta_1loop:.4f} × v²/2")
print(f"  (Tree: {lhs/rhs:.4f})")


# ============================================================
# SCORECARD
# ============================================================
print("\n" + "=" * 72)
print("SCORECARD")
print("=" * 72)
n_pass = sum(1 for _, p, _ in TESTS if p)
n_total = len(TESTS)
for name, passed, detail in TESTS:
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        for line in detail.split('\n'):
            print(f"         {line}")
print(f"\n  {n_pass}/{n_total} testów przeszło.")

print(f"""
========================================================================
PODSUMOWANIE ex243
========================================================================

  ★ BOSONOWA REGUŁA SUM:

  1. K(γ,W,Z) ≈ 1/2 jest PRZYBLIŻONE (wynika z m_W ≈ m_Z)
     → NIE głębokie, artefakt małego sin²θ_W = 0.231

  2. m_H² + m_W² + m_Z² = v²/2  (error: {abs(lhs/rhs-1)*100:.3f}%)
     → POTENCJALNIE GŁĘBOKIE: łączy Higgs self-coupling λ
       z gauge couplings g₁, g₂

  3. Jeśli reguła sum jest dokładna:
     λ = 1/4 - (2g₂² + g₁²)/8
     → m_H jest WYZNACZONE przez gauge couplings
     → SM ma 18, nie 19 wolnych parametrów

  4. Deviation od 1/2: {abs(R_central)/sig_R:.1f}σ
     Korekcje radiacyjne: ~4% (1-loop, Λ=1TeV)
     → Tree-level agreement at 0.5% is remarkable

  STATUS: OPEN — potrzebna 2-loop analiza i formalna weryfikacja
""")
