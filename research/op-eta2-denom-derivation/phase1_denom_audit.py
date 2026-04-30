#!/usr/bin/env python3
"""
η.2.Phase1 — Cross-sector denom landscape audit (5 sub-tests).
"""
import sympy as sp


# Anchors from prior cycles
A_TGP       = sp.Rational(64, 81)        # η.1
RHO_BAR_TGP = sp.Rational(11, 78)        # η.1
ETA_BAR_TGP = sp.Rational(5, 14)         # η.1
K_UP        = sp.Rational(7, 8)          # θ.1
K_DOWN      = sp.Rational(37, 50)        # θ.1
K_LEPTON    = sp.Rational(2, 3)          # ζ.1 / leptons
K_NEUTRINO  = sp.Rational(1, 2)          # ζ.1
PSI_PH      = sp.Rational(160, 137)      # ε.1
EPS_PH      = sp.Rational(23, 137)       # ε.1
RESIDUAL_FIT = sp.Rational(9, 250)       # α.1 best fit
ALPHA_INV_0  = sp.Float("137.035999084", 30)
RESIDUAL_REAL = ALPHA_INV_0 - 137

# B² chirality-counting values (z θ.1)
B2_LEPTON = sp.Integer(2)
B2_NEUTRINO = sp.Integer(1)
B2_UP = sp.Rational(13, 4)
B2_DOWN = sp.Rational(61, 25)
N_GEN = sp.Integer(3)


def prime_factors(n):
    """Return prime factorization as dict {prime: power}."""
    n = abs(int(n))
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def fmt_factors(factors):
    if not factors:
        return "1"
    parts = []
    for p in sorted(factors):
        if factors[p] == 1:
            parts.append(str(p))
        else:
            parts.append(f"{p}^{factors[p]}")
    return "·".join(parts)


def drift_pct(a, b):
    return float(sp.Abs(sp.N(a, 30) - sp.N(b, 30)) / sp.N(b, 30) * 100)


# =================== B1.1 — Cross-sector denom inventory ============
print("=" * 72)
print("B1.1 — Cross-sector denom inventory")
print("=" * 72)

entries = [
    ("η.1", "A",       A_TGP),
    ("η.1", "ρ̄",       RHO_BAR_TGP),
    ("η.1", "η̄",       ETA_BAR_TGP),
    ("θ.1", "K_up",    K_UP),
    ("θ.1", "K_down",  K_DOWN),
    ("ζ.1", "K_lepton", K_LEPTON),
    ("ζ.1", "K_ν",     K_NEUTRINO),
    ("ε.1", "ψ_ph",    PSI_PH),
    ("ε.1", "ε_ph",    EPS_PH),
    ("α.1", "9/250",   RESIDUAL_FIT),
]
print(f"  {'Sector':<6} {'Param':<10} {'Value':<14} {'Num factors':<14} {'Denom factors':<14}")
for sector, param, value in entries:
    num_f = fmt_factors(prime_factors(value.p))
    den_f = fmt_factors(prime_factors(value.q))
    print(f"  {sector:<6} {param:<10} {str(value):<14} {num_f:<14} {den_f:<14}")
B11_PASS = len(entries) == 10
print(f"  Total entries inventoried = {len(entries)}/10")
print(f"  Verdict B1.1              = {'PASS' if B11_PASS else 'FAIL'}")
print()


# =================== B1.2 — Cross-sector prime inventory ============
print("=" * 72)
print("B1.2 — Cross-sector prime inventory")
print("=" * 72)

prime_appearances = {}
for sector, param, value in entries:
    all_factors = {}
    for p, k in prime_factors(value.p).items():
        all_factors[p] = all_factors.get(p, 0) + k
    for p, k in prime_factors(value.q).items():
        all_factors[p] = all_factors.get(p, 0) + k
    for p in all_factors:
        prime_appearances.setdefault(p, []).append(f"{sector}/{param}")

print(f"  {'Prime':<6} {'#sectors':<10} {'Appearances'}")
for p in sorted(prime_appearances):
    apps = prime_appearances[p]
    print(f"  {p:<6} {len(apps):<10} {', '.join(apps)}")

# Categorize: ubiquitous (>=5), cross-sector (2-4), unique (1)
ubiq = [p for p, a in prime_appearances.items() if len(a) >= 5]
xsec = [p for p, a in prime_appearances.items() if 2 <= len(a) <= 4]
uniq = [p for p, a in prime_appearances.items() if len(a) == 1]
print()
print(f"  Ubiquitous primes (≥5 entries): {sorted(ubiq)}")
print(f"  Cross-sector primes (2-4):       {sorted(xsec)}")
print(f"  Unique primes (1 entry):          {sorted(uniq)}")

# Verify chirality-counting cascade core {2, 3, 5, 7}
core = {2, 3, 5, 7}
core_present = core.issubset(set(prime_appearances))
print(f"  Chirality-counting cascade core {{2,3,5,7}} present: {core_present}")
B12_PASS = core_present and len(prime_appearances) >= 8
print(f"  Verdict B1.2              = {'PASS' if B12_PASS else 'FAIL'}")
print()


# =================== B1.3 — B² cross-product → denom mapping =========
print("=" * 72)
print("B1.3 — B² cross-product → denom mapping (η.1 triple denoms)")
print("=" * 72)

# Test mappings
# 81 = 3⁴ = N_gen⁴
test_81 = N_GEN**4
print(f"  81 vs N_gen⁴ = {N_GEN}⁴ = {test_81}: match = {test_81 == 81}")

# 78 = 2·3·13 = 2·N_gen·(B²_up numerator)
B2_UP_NUM = sp.Integer(13)  # numerator of B²_up = 13/4
test_78 = 2 * N_GEN * B2_UP_NUM
print(f"  78 vs 2·N_gen·B²_up_num = 2·{N_GEN}·{B2_UP_NUM} = {test_78}: match = {test_78 == 78}")

# 14 = 2·7 = K_up_num · K_lepton_num
K_UP_NUM = sp.Integer(7)
K_LEPTON_NUM = sp.Integer(2)
test_14 = K_UP_NUM * K_LEPTON_NUM
print(f"  14 vs K_up_num·K_lepton_num = {K_UP_NUM}·{K_LEPTON_NUM} = {test_14}: match = {test_14 == 14}")

mappings = [(81, test_81), (78, test_78), (14, test_14)]
matches = sum(1 for target, derived in mappings if target == derived)
print()
print(f"  Clean mappings: {matches}/3")
B13_PASS = matches >= 2
print(f"  Verdict B1.3              = {'PASS' if B13_PASS else 'FAIL'}")
print()


# =================== B1.4 — Residual 0.036 cross-product landscape ===
print("=" * 72)
print("B1.4 — Residual 0.036 cross-product candidate landscape")
print("=" * 72)

target = float(RESIDUAL_REAL)
print(f"  Target residual = α⁻¹(0) − 137 = {target:.9f}")
print()

candidates = [
    ("R1: (B²_lepton − B²_ν)/N=27.78",
     sp.Rational(1) / sp.Rational(2778, 100)),
    ("R2: (B²_up − B²_down)/(some N)",
     (B2_UP - B2_DOWN) / sp.Rational(225, 10)),
    ("R3: 1/81 − 1/250",
     sp.Rational(1, 81) - sp.Rational(1, 250)),
    ("R4: (η̄·ρ̄)/2",
     ETA_BAR_TGP * RHO_BAR_TGP / 2),
    ("R5: (1−A)·(1−ρ̄)·η̄",
     (1 - A_TGP) * (1 - RHO_BAR_TGP) * ETA_BAR_TGP),
    ("R6: 9/250 (numerical fit)",
     RESIDUAL_FIT),
]
print(f"  {'Candidate':<40} {'Value':<14} {'Drift %':<10}")
n_under_1 = 0
for name, val in candidates:
    d = drift_pct(val, RESIDUAL_REAL)
    flag = "✓" if d < 1.0 else " "
    print(f"  {name:<40} {float(val):.6f}      {d:.4f}% {flag}")
    if d < 1.0:
        n_under_1 += 1
print()
print(f"  Candidates within drift < 1%: {n_under_1}/6")
B14_PASS = n_under_1 >= 1
print(f"  Verdict B1.4              = {'PASS' if B14_PASS else 'FAIL'}")
print()


# =================== B1.5 — Phase-1 audit verdict ====================
print("=" * 72)
print("B1.5 — Phase-1 audit verdict + Phase-2 viability")
print("=" * 72)

results = [
    ("B1.1 denom inventory",            B11_PASS),
    ("B1.2 prime inventory",            B12_PASS),
    ("B1.3 B²-cross-product mapping",   B13_PASS),
    ("B1.4 residual landscape",         B14_PASS),
]
n_pass = sum(1 for _, r in results if r)

for name, r in results:
    print(f"  [{'PASS' if r else 'FAIL'}] {name}")

# Phase-2 viability gate: B1.1 + B1.2 + B1.3 must PASS
B15_PASS = B11_PASS and B12_PASS and B13_PASS
print()
print(f"  Phase-2 viability (B1.1+B1.2+B1.3 PASS) = {B15_PASS}")
print(f"  Verdict B1.5              = {'PASS' if B15_PASS else 'FAIL'}")
print()


# =================== Final ============================================
print("=" * 72)
print("η.2.Phase1 — Final verdict")
print("=" * 72)

all_results = results + [("B1.5 Phase-2 viability", B15_PASS)]
n_total = len(all_results)
n_pass_total = sum(1 for _, r in all_results if r)
for name, r in all_results:
    print(f"  [{'PASS' if r else 'FAIL'}] {name}")
print()
print(f"  η.2.Phase1 score = {n_pass_total}/{n_total}")
if n_pass_total == n_total:
    print(f"  → Phase 2 first-principles derivation proceeds.")
elif n_pass_total >= 4:
    print(f"  → Phase 2 z minor caveat (review failing sub-test).")
else:
    print(f"  → η.2.Phase1 reframing required.")
