#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
ex273 — Analytical candidate g₀ᵉ = √(3/4) for the TGP soliton ODE
====================================================================

The published g₀ᵉ = 0.86941 is tantalizingly close to √(3/4) ≈ 0.86603.
This script investigates whether any simple algebraic expression reproduces
the numerical phi-fixed-point values from the soliton ODE.

Numerical references
--------------------
  published  g₀ᵉ       = 0.86941
  K = g²     phi-FP     = 0.86947     (ex271-class)
  K = g⁴     phi-FP     = 0.86778     (ex271-class)
"""

import math

# ── 1. Reference values ─────────────────────────────────────────────────
g0_pub    = 0.86941          # published value
g0_Kg2    = 0.86947          # numerical phi-FP, K = g²
g0_Kg4    = 0.86778          # numerical phi-FP, K = g⁴

# Physical / group-theory constants used in candidate expressions
Omega_L   = 0.6847           # Planck 2018  Ω_Λ
N         = 3                # SU(3) / generation count
GL3F2     = 168              # |GL(3, F₂)|
alpha_s   = 0.1190           # α_s(M_Z)

sqrt34    = math.sqrt(3.0/4.0)   # = √3 / 2 = cos 30°

# ── 2. The puzzle ────────────────────────────────────────────────────────
print("=" * 72)
print("  ex273 — Is g₀ᵉ = √(3/4) an analytical candidate?")
print("=" * 72)
print()
print(f"  √(3/4) = √3/2 = cos 30° = {sqrt34:.8f}")
print(f"  published g₀ᵉ             = {g0_pub:.8f}")
print(f"  K=g²  phi-FP              = {g0_Kg2:.8f}")
print(f"  K=g⁴  phi-FP              = {g0_Kg4:.8f}")
print()

# ── 3. Exact distances ──────────────────────────────────────────────────
print("─" * 72)
print("  Exact distances from √(3/4)")
print("─" * 72)
targets = [
    ("published  g₀ᵉ", g0_pub),
    ("K=g²  phi-FP  ", g0_Kg2),
    ("K=g⁴  phi-FP  ", g0_Kg4),
]
for label, val in targets:
    diff = val - sqrt34
    pct  = abs(diff) / val * 100
    print(f"  |{label} − √(3/4)| = {abs(diff):.6f}  ({pct:.4f} %)")
print()

# ── 4. Candidate algebraic expressions ──────────────────────────────────
print("─" * 72)
print("  Candidate algebraic expressions for g₀ᵉ")
print("─" * 72)

candidates = []

def add(name, val):
    candidates.append((name, val))

add("√(3/4)  = √3/2 = cos 30°",               sqrt34)
add("7/8",                                       7.0/8.0)
add("(7/8)(1 − 1/168)   GL(3,F₂) corr.",       (7.0/8.0) * (1.0 - 1.0/168.0))
add("1 − 1/N!  = 1 − 1/6",                     1.0 - 1.0/math.factorial(N))
add("1 − Ω_Λ/N!",                              1.0 - Omega_L / math.factorial(N))
add("Ω_Λ / N^(1/N)",                            Omega_L / N**(1.0/N))
add("(2/3)^(1/3)",                              (2.0/3.0)**(1.0/3.0))
add("(56/57) × 0.88  Higgs-sector",            (56.0/57.0) * 0.88)
add("cos(30°)  [= √3/2]",                      math.cos(math.radians(30)))
add("sin(arctan(√3/2))",                        math.sin(math.atan(math.sqrt(3)/2)))
add("3^(3/8) / 2",                              3.0**(3.0/8.0) / 2.0)
add("Ω_Λ^(1/N) × (N−1)",                       Omega_L**(1.0/N) * (N - 1))
add("(1 + √5)/2 − φ/2  (φ = golden)",         (1 + math.sqrt(5))/2 - ((1+math.sqrt(5))/2)/2)
# ── more targeted attempts ──
add("√3/2 × (1 + α_s/π)",                      sqrt34 * (1 + alpha_s / math.pi))
add("√3/2 × (1 + 1/GL(3,F₂)/4)",              sqrt34 * (1 + 1.0/(GL3F2 * 4)))
add("√3/2 × (1 + 1/(2π)²)",                    sqrt34 * (1 + 1.0/(2*math.pi)**2))
add("√(3/4 + α_s/(4π))",                        math.sqrt(3.0/4.0 + alpha_s/(4*math.pi)))
add("√(3/4 + 1/168)",                           math.sqrt(3.0/4.0 + 1.0/168.0))
add("1 − 1/(N! + GL(3,F₂))",                   1.0 - 1.0/(math.factorial(N) + GL3F2))
add("1 − 1/(N × GL(3,F₂)/4)",                  1.0 - 1.0/(N * GL3F2 / 4.0))
add("exp(−1/N!) = e^(−1/6)",                    math.exp(-1.0/math.factorial(N)))
add("(π/2)^(1/3)",                              (math.pi/2)**(1.0/3.0))
add("Ω_m^(1/3) = 0.3153^(1/3)",                0.3153**(1.0/3.0))
add("(2π − 2Ω_Λ) / (2π)",                      (2*math.pi - 2*Omega_L)/(2*math.pi))
add("1 − exp(−2)",                              1.0 - math.exp(-2))
add("(GL3F2 − 1)/GL3F2 × 0.876",              (GL3F2 - 1.0)/GL3F2 * 0.876)
add("tanh(√3/2)",                               math.tanh(math.sqrt(3)/2))
add("2 − π/√N",                                 2.0 - math.pi / math.sqrt(N))
add("N/(N + Ω_Λ/2)",                            N / (N + Omega_L/2.0))
add("(167/168)^(1/48)",                         (167.0/168.0)**(1.0/48.0))
add("Ω_Λ^(1/N) × √(Ω_Λ + Ω_m)",              Omega_L**(1.0/N) * math.sqrt(Omega_L + 0.3153))

# ── 5. Rank by closeness to EACH reference ──────────────────────────────

def rank_table(ref_label, ref_val):
    scored = []
    for name, val in candidates:
        diff = abs(val - ref_val)
        pct  = diff / ref_val * 100
        star = " ★" if pct < 0.10 else ""
        scored.append((pct, diff, val, name, star))
    scored.sort()
    print()
    print(f"  Ranked by closeness to {ref_label} = {ref_val:.5f}")
    print(f"  {'Expression':<46s} {'value':>10s} {'|Δ|':>10s} {'%':>8s}")
    print(f"  {'─'*46} {'─'*10} {'─'*10} {'─'*8}")
    for pct, diff, val, name, star in scored[:20]:
        print(f"  {name:<46s} {val:10.6f} {diff:10.6f} {pct:7.4f}%{star}")

rank_table("published g₀ᵉ",  g0_pub)
rank_table("K=g²  phi-FP",   g0_Kg2)
rank_table("K=g⁴  phi-FP",   g0_Kg4)

# ── 6. Radiative-correction hypothesis ──────────────────────────────────
print()
print("─" * 72)
print("  Radiative-correction hypothesis")
print("  g₀(phys) = g₀(bare) × (1 + C × α_s / π)")
print("─" * 72)
g0_bare = sqrt34
for ref_label, ref_val in [("published", g0_pub), ("K=g² FP", g0_Kg2), ("K=g⁴ FP", g0_Kg4)]:
    ratio = ref_val / g0_bare
    C = (ratio - 1.0) / (alpha_s / math.pi)
    print(f"  {ref_label:12s}: ratio = {ratio:.6f},  C = {C:.4f}")
    print(f"     → g₀(bare)=√(3/4), C α_s/π = {C * alpha_s / math.pi:.6f}")
print()
print("  Interpretation:")
print("    C ≈ 0.10  is a natural 1-loop coefficient.")
print("    If C = 1/(8π) ≈ 0.0398 → too small.")
c_pub = (g0_pub / g0_bare - 1.0) / (alpha_s / math.pi)
c_ideal_candidates = [
    ("1/(3!)",          1.0/6.0),
    ("1/(2π)",          1.0/(2*math.pi)),
    ("1/(4π)",          1.0/(4*math.pi)),
    ("C_F = 4/3",       4.0/3.0),
    ("C_A/4 = 3/4",     3.0/4.0),
    ("C_F − C_A/4 = 7/12", 7.0/12.0),
    ("1/8",             1.0/8.0),
    ("π/32",            math.pi/32.0),
    ("1/10",            0.1),
    ("3/(8π²)",         3.0/(8*math.pi**2)),
]
print(f"\n  C needed (published) = {c_pub:.6f}")
print(f"  Candidate C values and resulting g₀:")
for cname, cval in c_ideal_candidates:
    g0_pred = g0_bare * (1 + cval * alpha_s / math.pi)
    diff = abs(g0_pred - g0_pub)
    pct = diff / g0_pub * 100
    star = " ★" if pct < 0.10 else ""
    print(f"    C = {cname:<20s} = {cval:.6f}  →  g₀ = {g0_pred:.6f}  Δ = {pct:.4f}%{star}")

# ── 7. Summary verdict ──────────────────────────────────────────────────
print()
print("=" * 72)
print("  VERDICT")
print("=" * 72)
print("""
  • √(3/4) = 0.86603 is 0.39% below the published g₀ᵉ = 0.86941.
    Neither K=g² nor K=g⁴ soliton ODE formulation predicts √(3/4) exactly.

  • The radiative-correction route:
      g₀(phys) = √(3/4) × (1 + C × α_s/π)
    requires C ≈ {c_pub:.2f}, which is a plausible 1-loop Wilson coefficient.

  • Best simple-fraction matches to the published value are shown in
    the ranked tables above.  Any expression matching < 0.1% gets a ★.

  • Conclusion: √(3/4) as a *bare* coupling with a 1-loop QCD
    dressing is the most economical analytical hypothesis, but it
    is not yet derivable from the soliton ODE alone — it would need
    a separate RG argument linking g₀(bare) to g₀(phys).
""".format(c_pub=c_pub))
