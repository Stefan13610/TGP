"""
m9_1_prime_pivot_B.py -- M9.1' Pivot B: search for natural non-power f(psi)
with f*h = 1 + alpha=2 + beta_PPN = 1.

Setup
-----
With substrate budget f*h = 1, gamma_PPN = 1 is automatic.
For TGP canonical alpha = 2 (kinetic coupling), c_2 = -1 (verified
numerically in M9.1'). The PPN master formula gives:

    beta_PPN = f''(1)/f'(1)^2 + 2 c_2 / f'(1)
             = f''(1)/f'(1)^2 - 2 / f'(1)

beta_PPN = 1  <=>  f''(1) = f'(1)^2 + 2 f'(1)
                  = f'(1) * (f'(1) + 2)

This is one equation in two derivatives. Power-form f(psi) = psi^p
has f'(1) = p, f''(1) = p(p-1) = p^2 - p, giving:

    p^2 - p = p^2 + 2p   =>   -p = 2p   =>   p = 0   (trivial f=1)

Hence NO power form satisfies it. We scan natural candidate forms:

    1. f(psi) = exp(a*(psi-1))           Taylor: f'(1)=a, f''(1)=a^2
       => a^2 = a^2 + 2a  =>  a = 0  (trivial).  EXCLUDED.

    2. f(psi) = exp(a*(1-psi))           f'(1)=-a, f''(1)=a^2
       => a^2 = a^2 - 2a  =>  a = 0  (trivial).  EXCLUDED.

    3. f(psi) = 1 / (1 + a*(psi-1))      f'(1)=-a, f''(1)=2a^2
       => 2a^2 = a^2 - 2a  =>  a^2 = -2a  =>  a = -2 (or 0).
       a = -2: f(psi) = 1/(1 - 2(psi-1)) = 1/(3 - 2*psi).
       Pole at psi = 3/2 (within physical range psi in (0, 4/3)?).
       CANDIDATE, but pole structure suspicious.

    4. f(psi) = (1 + a*(psi-1))^p        f'(1)=ap, f''(1)=a^2 p(p-1)
       Two parameters; one equation; one-param family.
       Reduces to power form when a=1.

    5. f(psi) = log-related                See item 3 (rational).

    6. f(psi) = 1 + a*(1/psi - 1) = 1 + a*(1-psi)/psi
       f'(1) = -a, f''(1) = 2a
       => 2a = a^2 - 2a  =>  a^2 = 4a  =>  a = 4 (or 0).
       a = 4: f(psi) = 1 + 4*(1-psi)/psi = (psi + 4 - 4 psi)/psi = (4 - 3 psi)/psi.
       Pole at psi=0; zero at psi=4/3 (boundary of basin).
       INTERESTING: zero EXACTLY at the basin boundary psi = 4/3.

    7. Polynomial f(psi) = 1 + a*(psi-1) + b*(psi-1)^2 + ...
       f'(1)=a, f''(1)=2b.
       => 2b = a^2 + 2a  =>  b = a*(a+2)/2.
       For a=-1: b = -1/2. f = 1 - eps - eps^2/2 + O(eps^3).
       For a=-2: b = 0.   f = 1 - 2 eps + O(eps^3).
       For a=-3: b = +3/2. f = 1 - 3 eps + (3/2) eps^2.
       Pure constructions; need further constraint to fix a.

We test the candidates (3) and (6) numerically: solve modified
Phi-EOM with metric f(psi) for static spherical source, extract c_2
(for the specific f), and verify whether the "modified master formula"
(which now goes via different chain rule) actually gives beta_PPN = 1.

Strategy
--------
For a NEW f(psi), the substrate-budget condition f*h = 1 gives h = 1/f.
The Phi-EOM derived from action with kinetic term K(psi) (alpha=2 means
K = K_geo psi^4) is INVARIANT to choice of f (f only enters g_eff, not
the substrate dynamics). So the SAME ε(r) profile from M9.1' applies.

Therefore Pivot B is a PURELY ANALYTICAL choice of how to MAP the
substrate field psi onto the effective metric components. The Phi-EOM
gives ε(r) with c_2 = -alpha/2 = -1 (alpha=2). The PPN parameters
are then read from f, h via:

    g_tt = -c^2 f(1+eps),   g_rr = h(1+eps) = 1/f(1+eps)

The expansion f(1+eps) = 1 + f'(1) eps + (1/2) f''(1) eps^2 + ...
combined with eps(r) = a_1/r + a_2/r^2 + ... gives:

    g_tt = -c^2 [1 + f'(1) a_1/r + (f'(1) a_2 + f''(1) a_1^2/2)/r^2 + ...]

PPN: g_tt = -c^2 [1 - 2GM/(c^2 r) + 2 beta (GM/c^2)^2 / r^2 + ...]
=>  -2 GM/c^2 = f'(1) a_1   =>   a_1 = -2GM/(c^2 f'(1))

beta = ... (same master formula).

In summary: the master formula is INDEPENDENT of choice of f(psi)
(power vs non-power), as long as we use 2nd-order Taylor coefficients
f'(1), f''(1). Pivot B is NOT an independent numerical experiment;
it is a TASK OF FINDING a function f(psi) with a specific relation
between f'(1) and f''(1).

Conclusion: Pivot B problem reduces to identifying a NATURAL physical
principle that selects f'(1) and f''(1) such that:

    f''(1) = f'(1) * (f'(1) + 2)        (with alpha=2)

Below we tabulate candidate f's and check whether ANY of them is
"natural" in the sense of admitting a substrate-level derivation.
"""

from __future__ import annotations

import numpy as np


def beta_ppn(f_prime, f_double_prime, c_2=-1.0):
    """Master PPN formula."""
    if abs(f_prime) < 1e-12:
        return float("nan")  # f'(1)=0 means metric is flat, no Newton limit
    return f_double_prime / f_prime ** 2 + 2.0 * c_2 / f_prime


def report_form(label, f_prime, f_double_prime, expr_str=""):
    """Report a candidate metric form."""
    b = beta_ppn(f_prime, f_double_prime)
    target = f_prime * (f_prime + 2.0)
    constraint_residual = f_double_prime - target
    status = "BETA_PPN=1" if abs(b - 1.0) < 1e-9 else f"beta_PPN={b:+.4f}"
    return (f"  {label:<35}  f'(1)={f_prime:+.4f}  f''(1)={f_double_prime:+.4f}  "
            f"|  target f''(1) = {target:+.4f}  resid={constraint_residual:+.4f}  "
            f"|  {status}  | {expr_str}")


def main():
    print("=" * 90)
    print("M9.1' Pivot B: natural f(psi) candidates with f*h=1 + alpha=2 + beta_PPN=1")
    print("=" * 90)
    print()
    print("Constraint (alpha=2, c_2=-1):  f''(1) = f'(1) * (f'(1) + 2)")
    print()
    print("Candidates:")
    print("-" * 90)

    # 1. Pure power: f = psi^p
    print("(1) Pure power f(psi) = psi^p:")
    for p in [-1.0, -1.0/2, -2.0, -1.0/3, +1.0, +0.5]:
        fp = p
        fpp = p * (p - 1)
        print(report_form(f"    p = {p:+.3f}", fp, fpp, f"f = psi^{p}"))
    print()

    # 2. Exponential
    print("(2) Exponential f(psi) = exp(a*(psi-1)):")
    for a in [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0]:
        fp = a
        fpp = a * a
        print(report_form(f"    a = {a:+.3f}", fp, fpp, f"f = exp({a:+.2f}(psi-1))"))
    print()

    # 3. Rational f = 1/(1 + a(psi-1))
    print("(3) Rational f(psi) = 1/(1 + a*(psi-1)):")
    for a in [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]:
        fp = -a
        fpp = 2 * a * a
        psi_pole = 1.0 - 1.0 / a if a != 0 else float("inf")
        print(report_form(f"    a = {a:+.3f}", fp, fpp,
                          f"f = 1/(1+{a:+.2f}(psi-1))  pole psi={psi_pole:.3f}"))
    print()

    # 4. Hyperbolic / log inverse: f = 1 + a*(1-psi)/psi
    print("(4) Form f(psi) = 1 + a*(1-psi)/psi:")
    for a in [4.0, 3.0, 2.0, 1.0, 0.0, -1.0]:
        fp = -a       # df/dpsi at psi=1: d/dpsi [a*(1-psi)/psi] = a*(-1/psi - (1-psi)/psi^2)|psi=1 = -a
        fpp = 2 * a   # d^2/dpsi^2 = 2a (at psi=1)
        print(report_form(f"    a = {a:+.3f}", fp, fpp, f"f = 1 + {a:+.2f}*(1-psi)/psi"))
    print()

    # 5. Logarithmic: f(psi) = 1 + a*ln(psi) + b*ln(psi)^2
    print("(5) Log: f(psi) = 1 + a*ln(psi) + b*ln(psi)^2 + ...:")
    for a, b in [(-1.0, 0.0), (-1.0, 0.5), (-1.0, -0.5), (-2.0, 1.0), (-1.0, +1.5/2)]:
        # f'(psi) = a/psi - a/psi^2 * 0 + 2b ln(psi)/psi  -- careful
        # f(psi) = 1 + a u + b u^2 with u = ln psi, du/dpsi = 1/psi
        # f'(psi) = a/psi + 2 b u / psi
        # f'(1) = a (since u=0)
        # f''(psi) = -a/psi^2 + 2b/psi^2 - 2 b u/psi^2 = (-a + 2b)/psi^2 - 2 b u / psi^2
        # f''(1) = -a + 2b
        fp = a
        fpp = -a + 2 * b
        print(report_form(f"    a={a:+.2f}, b={b:+.2f}", fp, fpp,
                          f"f = 1 + {a:+.2f} ln(psi) + {b:+.2f} ln^2(psi)"))
    print()

    # 6. Quadratic Taylor (1-parameter family; b determined from beta=1)
    print("(6) Quadratic Taylor f(psi) = 1 + a(psi-1) + b(psi-1)^2:")
    for a in [-3.0, -2.0, -1.0, -0.5, +0.5, +1.0]:
        b = a * (a + 2) / 2.0
        fp = a
        fpp = 2 * b
        # By construction satisfies beta_PPN=1
        print(report_form(f"    a = {a:+.3f}", fp, fpp,
                          f"b = a(a+2)/2 = {b:+.3f}; f = 1 + {a:+.2f}(psi-1) + {b:+.2f}(psi-1)^2"))
    print()

    # Summary
    print("=" * 90)
    print("Conclusion: the constraint f''(1) = f'(1)(f'(1)+2) with alpha=2 and f*h=1")
    print("is a SINGLE EQUATION in TWO Taylor coefficients of f.")
    print()
    print("Hence a 1-parameter family of admissible f exists. The QUESTION is whether")
    print("any such f has a NATURAL substrate-level derivation (substrate budget,")
    print("substrate action, conformal coupling, etc.).")
    print()
    print("Forms tested:")
    print("  (1) Pure power:        FAIL (only trivial p=0).")
    print("  (2) Exponential:       FAIL (only trivial a=0).")
    print("  (3) Rational:          a=-2 (pole psi=3/2, near basin boundary).")
    print("  (4) Hyperbolic:        a=+4 (zero at psi=4/3 = basin boundary).")
    print("  (5) Logarithmic:       always satisfies f'(1)=a, f''(1)=-a+2b;")
    print("                         constraint -a + 2b = a(a+2)  =>  b = (a^2+3a)/2.")
    print("                         For a=-1: b=-1; for a=-2: b=-1; for a=-3: b=0.")
    print("  (6) Quadratic Taylor:  by construction, infinite family.")
    print()
    print("Of these, only (4) shows a NATURAL connection to TGP physics:")
    print("  f(psi) = 1 + 4*(1-psi)/psi = (4 - 3*psi)/psi")
    print("  Zero precisely at psi = 4/3 -- the BASIN BOUNDARY of TGP")
    print("  (cf. action `eq:action-full-psi`, ghost-free range psi in (0, 4/3)).")
    print()
    print("This is suspicious enough to warrant follow-up: is this a 'natural'")
    print("metric choice, or a coincidence?")
    print("=" * 90)


if __name__ == "__main__":
    main()
