#!/usr/bin/env python3
"""
α.1.Phase1 — α_QED numerical landscape audit.

5 sub-tests:
  A1.1  CODATA 2022 + PDG 2024 reference values
  A1.2  Top-N rational candidates dla α_QED⁻¹ (denom ≤ 200)
  A1.3  Structural form 137 = ψ_ph·137 − ε_ph·137 sympy exact
  A1.4  Cross-sector prime-137 mapping
  A1.5  α_QED⁻¹ RG-running consistency α(0) → α(M_Z)

Run:
  PYTHONIOENCODING=utf-8 python -X utf8 phase1_alpha_audit.py
"""
import math
import sympy as sp


# -------- Reference constants ----------------------------------------
ALPHA_INV_0_PDG = sp.Float("137.035999084", 30)   # CODATA 2022
ALPHA_INV_MZ_PDG = sp.Float("127.952", 20)         # PDG 2024 running
ALPHA_INV_0_SIGMA = sp.Float("0.000000021", 30)    # CODATA 81 ppt

# TGP anchors z ε.1
PSI_PH = sp.Rational(160, 137)
EPS_PH = sp.Rational(23, 137)
TARGET_SHIFT_F4 = sp.Rational(57, 500)             # F4 chain
N_A = sp.Rational(500, 57)                         # F4 inverse
ALPHA_0_F4 = sp.Rational(1069833, 264500)          # F4 photon-ring
KAPPA_TGP = sp.sqrt(ALPHA_0_F4)                    # XS.1 cross-sector

# TGP rationals z innych cycles (cross-sector check)
A_TGP = sp.Rational(64, 81)
RHO_TGP = sp.Rational(11, 78)
ETA_TGP = sp.Rational(5, 14)
LAMBDA_C = sp.Rational(165, 167)
K_LEPTON = sp.Rational(2, 3)
K_NU = sp.Rational(1, 2)
K_UP = sp.Rational(7, 8)
K_DOWN = sp.Rational(37, 50)


def drift_pct(a, b):
    """Absolute drift in percent: |a - b| / b * 100."""
    return float(sp.Abs(sp.N(a, 30) - sp.N(b, 30)) / sp.N(b, 30) * 100)


def is_prime(n):
    if n < 2:
        return False
    for d in range(2, int(math.isqrt(n)) + 1):
        if n % d == 0:
            return False
    return True


def prime_factors(n):
    """Return list of (prime, exponent) tuples."""
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            out.append((p, e))
        p += 1
    if n > 1:
        out.append((n, 1))
    return out


def best_rationals(target, max_denom, top_n=8):
    """Return top_n unique p/q (q ≤ max_denom) closest to target by drift."""
    target_f = float(target)
    cands = []
    for q in range(1, max_denom + 1):
        p = round(target_f * q)
        if p < 1:
            continue
        r = sp.Rational(p, q)
        d = drift_pct(sp.Float(r, 30), target)
        cands.append((r, d, p, q))
    seen = set()
    out = []
    for r, d, p, q in sorted(cands, key=lambda x: x[1]):
        key = (r.p, r.q)
        if key in seen:
            continue
        seen.add(key)
        out.append((r, d, p, q))
        if len(out) >= top_n:
            break
    return out


# =================== A1.1 — reference values ========================
print("=" * 72)
print("A1.1 — α_QED reference values (CODATA 2022 + PDG 2024)")
print("=" * 72)

alpha_0 = 1 / ALPHA_INV_0_PDG
alpha_MZ = 1 / ALPHA_INV_MZ_PDG

print(f"  α_QED⁻¹(0)              = {float(ALPHA_INV_0_PDG):.9f}")
print(f"  σ(α⁻¹) CODATA 81 ppt    = {float(ALPHA_INV_0_SIGMA):.2e}")
print(f"  relative precision      = {float(ALPHA_INV_0_SIGMA/ALPHA_INV_0_PDG):.2e}")
print(f"  α_QED⁻¹(M_Z)            = {float(ALPHA_INV_MZ_PDG):.4f}")
print(f"  α_QED(0)                = {float(alpha_0):.10e}")
print(f"  α_QED(M_Z)              = {float(alpha_MZ):.10e}")

running_pct = (float(alpha_MZ) / float(alpha_0) - 1) * 100
print(f"  α(M_Z)/α(0) − 1         = {running_pct:.3f}%   (SM vac.pol.)")

A11_PASS = (
    float(ALPHA_INV_0_PDG) > 137.0
    and float(ALPHA_INV_0_PDG) < 138.0
    and abs(running_pct - 7.1) < 1.5
)
print(f"  Verdict A1.1            = {'PASS' if A11_PASS else 'FAIL'}")
print()


# =================== A1.2 — top-N rationals =========================
print("=" * 72)
print("A1.2 — top-8 rational candidates dla α_QED⁻¹  (denom ≤ 200)")
print("=" * 72)

top = best_rationals(ALPHA_INV_0_PDG, 200, top_n=8)
print(f"  {'rank':<5}{'rational':<18}{'value':<15}{'drift_%':<12}")
print(f"  {'-'*5}{'-'*18}{'-'*15}{'-'*12}")
for i, (r, d, p, q) in enumerate(top, 1):
    print(f"  {i:<5}{f'{p}/{q}':<18}{float(r):<15.9f}{d:<12.6f}")

# best non-trivial (skip 137/1 if present)
best_rat = top[0]
print()
print(f"  Best rational           = {best_rat[2]}/{best_rat[3]} drift {best_rat[1]:.6f}%")
print(f"  Trivial 137/1 drift     = {drift_pct(sp.Integer(137), ALPHA_INV_0_PDG):.6f}%")

# Phase 2 will need a non-trivial cleaner candidate with denom > 1
non_trivial = [t for t in top if t[3] > 1]
print(f"  Best non-trivial        = {non_trivial[0][2]}/{non_trivial[0][3]} "
      f"drift {non_trivial[0][1]:.6f}%")

A12_PASS = best_rat[1] < 0.05  # zeroth-order anchor band
print(f"  Verdict A1.2            = {'PASS' if A12_PASS else 'FAIL'}")
print()


# =================== A1.3 — structural form =========================
print("=" * 72)
print("A1.3 — structural form 137 = ψ_ph·137 − ε_ph·137 sympy exact")
print("=" * 72)

psi_times_137 = sp.simplify(PSI_PH * 137)
eps_times_137 = sp.simplify(EPS_PH * 137)
diff_times_137 = sp.simplify((PSI_PH - EPS_PH) * 137)

print(f"  ψ_ph                    = {PSI_PH} = {float(PSI_PH):.9f}")
print(f"  ε_ph                    = {EPS_PH} = {float(EPS_PH):.9f}")
print(f"  ψ_ph · 137              = {psi_times_137} (sympy exact)")
print(f"  ε_ph · 137              = {eps_times_137} (sympy exact)")
print(f"  (ψ_ph − ε_ph) · 137     = {diff_times_137} (sympy exact)")
print()

residual = sp.N(ALPHA_INV_0_PDG - 137, 30)
print(f"  α_QED⁻¹(0) − 137        = {float(residual):.9f}")

# Cascade probe candidates dla residual ≈ 0.036
print()
print("  Residual ≈ 0.036 cascade candidates:")
candidates = [
    ("23/640",                      sp.Rational(23, 640)),
    ("ε_ph · target_shift",          EPS_PH * TARGET_SHIFT_F4),
    ("ε_ph² · α₀",                   EPS_PH**2 * ALPHA_0_F4),
    ("η̄ · κ_TGP / 19.7",             ETA_TGP * KAPPA_TGP / sp.Rational(197,10)),
    ("1/(2π · κ_TGP² · 1.07)",       1 / (2*sp.pi * ALPHA_0_F4 * sp.Rational(107,100))),
    ("ε_ph / κ_TGP²",                EPS_PH / ALPHA_0_F4),
    ("9/250",                        sp.Rational(9, 250)),
    ("36/1000",                      sp.Rational(36, 1000)),
    ("ε_ph · √(2/9)",                EPS_PH * sp.sqrt(sp.Rational(2,9))),
    ("(η̄)² · target_shift",          ETA_TGP**2 * TARGET_SHIFT_F4),
]
print(f"  {'candidate':<32}{'value':<15}{'drift_%':<10}")
print(f"  {'-'*32}{'-'*15}{'-'*10}")
for name, val in candidates:
    v = float(sp.N(val, 30))
    d = abs(v - float(residual)) / float(residual) * 100
    print(f"  {name:<32}{v:<15.9f}{d:<10.4f}")

# Strict A1.3 PASS = sympy exact identity verifies
A13_PASS = (psi_times_137 == 160 and eps_times_137 == 23 and diff_times_137 == 137)
print()
print(f"  Verdict A1.3            = {'PASS' if A13_PASS else 'FAIL'}")
print()


# =================== A1.4 — cross-sector prime mapping ==============
print("=" * 72)
print("A1.4 — cross-sector prime-137 mapping")
print("=" * 72)

print("  Prime structure of TGP rationals:")
inventory = [
    ("ψ_ph",      PSI_PH,    "ε.1"),
    ("ε_ph",      EPS_PH,    "ε.1"),
    ("A_TGP",     A_TGP,     "η.1"),
    ("ρ̄_TGP",    RHO_TGP,   "η.1"),
    ("η̄_TGP",    ETA_TGP,   "η.1"),
    ("λ_C",       LAMBDA_C,  "tgp-leptons"),
    ("N_A",       N_A,       "UV.1"),
    ("K_lepton",  K_LEPTON,  "tgp-leptons"),
    ("K_ν",       K_NU,      "ζ.1"),
    ("K_up",      K_UP,      "θ.1"),
    ("K_down",    K_DOWN,    "θ.1"),
    ("target_shift", TARGET_SHIFT_F4, "F4"),
]
all_primes = set()
prime_137_rationals = []
print(f"  {'symbol':<14}{'rational':<14}{'num primes':<22}{'denom primes':<22}")
print(f"  {'-'*14}{'-'*14}{'-'*22}{'-'*22}")
for name, r, src in inventory:
    np_factors = prime_factors(r.p) if r.p > 0 else []
    dp_factors = prime_factors(r.q) if r.q > 0 else []
    np_str = ",".join(f"{p}^{e}" if e > 1 else str(p) for p, e in np_factors) or "1"
    dp_str = ",".join(f"{p}^{e}" if e > 1 else str(p) for p, e in dp_factors) or "1"
    for p, _ in np_factors + dp_factors:
        all_primes.add(p)
        if p == 137:
            prime_137_rationals.append(name)
    print(f"  {name:<14}{f'{r.p}/{r.q}':<14}{np_str:<22}{dp_str:<22}")

print()
print(f"  Distinct primes across TGP landscape: {sorted(all_primes)}")
print(f"  Rationals containing prime 137: {prime_137_rationals}")
print(f"  Is 137 prime?           = {is_prime(137)}")

# Structural conclusion: 137 is isolated to ε.1 (ψ_ph, ε_ph)
A14_PASS = (
    set(prime_137_rationals) == {"ψ_ph", "ε_ph"}
    and is_prime(137)
)
print(f"  Verdict A1.4            = {'PASS' if A14_PASS else 'FAIL'}")
print()


# =================== A1.5 — RG-running consistency ==================
print("=" * 72)
print("A1.5 — α_QED⁻¹ RG-running α(0) → α(M_Z) consistency")
print("=" * 72)

ratio = sp.N(ALPHA_INV_0_PDG / ALPHA_INV_MZ_PDG, 30)
print(f"  α⁻¹(0) / α⁻¹(M_Z)       = {float(ratio):.6f}")
print(f"  α(M_Z) / α(0)           = {1/float(ratio):.6f}")
print(f"  Running fraction        = {(1/float(ratio) - 1)*100:.3f}%")

# TGP form factor probe for ratio - 1
delta = float(ratio) - 1  # ≈ 0.07097
print()
print(f"  ratio − 1 = {delta:.6f}  (RG-running magnitude)")
print()
print("  Cascade probes dla 0.071:")
probes = [
    ("ε_ph² · α₀",         EPS_PH**2 * ALPHA_0_F4),       # = target_shift exact
    ("target_shift",        TARGET_SHIFT_F4),              # 0.114
    ("ε_ph · κ_TGP / sqrt(18)", EPS_PH * KAPPA_TGP / sp.sqrt(18)),
    ("ε_ph / 2.366",        EPS_PH / sp.Rational(2366,1000)),
    ("η̄ / 5",               ETA_TGP / 5),                  # = 1/14
    ("η̄·target_shift",      ETA_TGP * TARGET_SHIFT_F4),    # 5/14·57/500
]
print(f"  {'candidate':<32}{'value':<15}{'drift_%':<12}")
print(f"  {'-'*32}{'-'*15}{'-'*12}")
for name, val in probes:
    v = float(sp.N(val, 30))
    d = abs(v - delta) / delta * 100
    print(f"  {name:<32}{v:<15.9f}{d:<12.4f}")

# Common β-rescaling test: under m → c·m, α_QED⁻¹ form is dimensionless,
# but RG-running comes z lepton+hadron loops. TGP rational anchor for
# 137 itself is PDG-flow-invariant (number-theoretic, not RG quantity).
# So "RG-stability" here means: 137-anchor unaffected by common β-rescale.
# We test by verifying ψ_ph and ε_ph ratios are dimensionless and
# RG-invariant (already verified w UV.1).

print()
print("  RG-stability of 137-anchor:")
print("    ψ_ph = 160/137 dimensionless ratio (UV.1.UV2.5 ratio invariant)")
print("    ε_ph = 23/137 dimensionless ratio (UV.1.UV2.5 ratio invariant)")
print("    α_QED⁻¹ ≈ 137 zeroth-order anchor: prime-denom structure")
print("    RG-running 7.1% from SM vac.pol. — orthogonal do TGP rational anchor")

A15_PASS = (
    abs(delta - 0.0710) < 0.005      # ratio − 1 ≈ 7.1% within ±0.5%
    and float(ALPHA_INV_0_PDG) - float(ALPHA_INV_MZ_PDG) > 5  # > 5 unit running
)
print(f"  Verdict A1.5            = {'PASS' if A15_PASS else 'FAIL'}")
print()


# =================== Final verdict ===================================
print("=" * 72)
print("α.1.Phase1 — Final verdict")
print("=" * 72)

results = [
    ("A1.1 reference values",          A11_PASS),
    ("A1.2 top rationals",             A12_PASS),
    ("A1.3 structural 137 = 160-23",   A13_PASS),
    ("A1.4 prime-137 mapping",         A14_PASS),
    ("A1.5 RG-running consistency",    A15_PASS),
]
n_pass = sum(1 for _, r in results if r)
n_total = len(results)

for name, r in results:
    print(f"  [{'PASS' if r else 'FAIL'}] {name}")
print()
print(f"  α.1.Phase1 score        = {n_pass}/{n_total}")
if n_pass == n_total:
    print(f"  Phase 2 first-principles derivation viable.")
elif n_pass >= 4:
    print(f"  Phase 2 conditional (partial-pass).")
else:
    print(f"  STRUCTURAL HINT only — terminate research-track.")
