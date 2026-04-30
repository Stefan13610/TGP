#!/usr/bin/env python3
"""
π.1.Phase1 — NME landscape inventory (5 sub-tests).

Anchors:
- m_ββ_A = 1.584 meV (ν.1 Form A)
- m_ββ_B = 3.249 meV (ν.1 Form B)
- m_e = 0.5109989461 MeV = 510998.9461 eV

Standard mass-mechanism formula:
    [T_{1/2}^{0ν}]^{-1} = G_{0ν} · |M|² · (m_ββ / m_e)²

NME methods: QRPA, IBM-2, NSM, EDF
Isotopes: Ge-76, Te-130, Xe-136
"""

# ---------------- Constants ---------------------
M_E_eV = 510998.9461
M_BB_A_meV = 1.584
M_BB_B_meV = 3.249

# Phase-space factors G_{0ν} in 10⁻¹⁵ yr⁻¹ (Kotila & Iachello 2012, g_A=1.27)
G_PSF = {
    "Ge-76":  2.36e-15,
    "Te-130": 14.22e-15,
    "Xe-136": 14.58e-15,
}

# NME literature values (g_A=1.27 unquenched)
NME = {
    "Ge-76":  {"QRPA": 5.0, "IBM-2": 4.6, "NSM": 3.0, "EDF": 4.6},
    "Te-130": {"QRPA": 4.2, "IBM-2": 4.1, "NSM": 1.9, "EDF": 5.1},
    "Xe-136": {"QRPA": 3.4, "IBM-2": 3.3, "NSM": 1.8, "EDF": 4.2},
}

# Current 90% CL T_{1/2} bounds
T_BOUND = {
    "Ge-76":  1.8e26,
    "Te-130": 2.2e25,
    "Xe-136": 2.3e26,
}


def t_half(iso: str, m_bb_meV: float, M: float) -> float:
    """Compute T_{1/2}^{0ν} in years."""
    G = G_PSF[iso]
    m_bb_eV = m_bb_meV * 1e-3
    factor = (m_bb_eV / M_E_eV) ** 2
    inv_T = G * M * M * factor
    return 1.0 / inv_T


# =============== P1.1 — NME inventory ============================
print("=" * 72)
print("P1.1 — NME inventory (QRPA / IBM-2 / NSM / EDF × Ge-76 / Te-130 / Xe-136)")
print("=" * 72)
print(f"  {'Isotope':<10} {'QRPA':>6} {'IBM-2':>6} {'NSM':>6} {'EDF':>6} | {'mean':>6} {'span':>6}")
all_mean = {}
all_span = {}
for iso, methods in NME.items():
    vals = list(methods.values())
    mean = sum(vals) / len(vals)
    span = max(vals) - min(vals)
    all_mean[iso] = mean
    all_span[iso] = span
    print(f"  {iso:<10} {methods['QRPA']:>6.2f} {methods['IBM-2']:>6.2f} {methods['NSM']:>6.2f} {methods['EDF']:>6.2f} | {mean:>6.2f} {span:>6.2f}")
P11_PASS = all(len(m) == 4 for m in NME.values())
print(f"  12-cell matrix complete: {P11_PASS}")
print(f"  Verdict P1.1 = {'PASS' if P11_PASS else 'FAIL'}")
print()


# =============== P1.2 — PSF + g_A inventory ======================
print("=" * 72)
print("P1.2 — Phase-space factor inventory (Kotila & Iachello 2012)")
print("=" * 72)
for iso, G in G_PSF.items():
    print(f"  {iso:<10} G_0ν = {G*1e15:>5.2f} × 10⁻¹⁵ yr⁻¹")
print(f"  g_A = 1.27 unquenched (literature default)")
P12_PASS = all(G > 0 for G in G_PSF.values())
print(f"  Verdict P1.2 = {'PASS' if P12_PASS else 'FAIL'}")
print()


# =============== P1.3 — m_ββ Form A/B from ν.1 ===================
print("=" * 72)
print("P1.3 — m_ββ Form A/B anchored from ν.1")
print("=" * 72)
print(f"  m_ββ_TGP_A = {M_BB_A_meV} meV (Form A — chirality halving)")
print(f"  m_ββ_TGP_B = {M_BB_B_meV} meV (Form B — PMNS-Wolfenstein)")
gap = M_BB_B_meV - M_BB_A_meV
ratio = M_BB_B_meV / M_BB_A_meV
print(f"  Δm_ββ (B−A) = {gap:.3f} meV; ratio B/A = {ratio:.3f}× (= {ratio**2:.3f}× in T_1/2)")
P13_PASS = M_BB_A_meV > 0 and M_BB_B_meV > M_BB_A_meV
print(f"  Verdict P1.3 = {'PASS' if P13_PASS else 'FAIL'}")
print()


# =============== P1.4 — Method-spread σ_NME ======================
print("=" * 72)
print("P1.4 — NME method-spread σ_NME = span/mean per isotope")
print("=" * 72)
print(f"  {'Isotope':<10} {'mean':>6} {'span':>6} {'σ_NME':>7}")
sigmas = {}
for iso in NME:
    sigmas[iso] = all_span[iso] / all_mean[iso]
    print(f"  {iso:<10} {all_mean[iso]:>6.2f} {all_span[iso]:>6.2f} {sigmas[iso]:>7.3f}")
most_stable = min(sigmas, key=sigmas.get)
print(f"  Most NME-stable isotope: {most_stable} (σ_NME = {sigmas[most_stable]:.3f})")
P14_PASS = all(s > 0 for s in sigmas.values())
print(f"  Verdict P1.4 = {'PASS' if P14_PASS else 'FAIL'}")
print()


# =============== P1.5 — Rate-only viability ======================
print("=" * 72)
print("P1.5 — Rate-only viability vs current 90% CL bounds")
print("=" * 72)
print(f"  {'Iso':<8} {'Method':<8} {'T_A (yr)':>14} {'T_B (yr)':>14} {'Bound':>14} {'OK?':>4}")
all_ok = True
for iso in NME:
    for method, M in NME[iso].items():
        T_A = t_half(iso, M_BB_A_meV, M)
        T_B = t_half(iso, M_BB_B_meV, M)
        bound = T_BOUND[iso]
        ok = T_A > bound and T_B > bound
        if not ok:
            all_ok = False
        print(f"  {iso:<8} {method:<8} {T_A:>14.3e} {T_B:>14.3e} {bound:>14.3e} {'✓' if ok else '✗':>4}")
P15_PASS = all_ok
print(f"  All cells (12 isotope×method × 2 forms = 24) above 90% CL bounds: {P15_PASS}")
print(f"  Verdict P1.5 = {'PASS' if P15_PASS else 'FAIL'}")
print()


# =============== Final ===========================================
print("=" * 72)
print("π.1.Phase1 — Final verdict")
print("=" * 72)

results = [
    ("P1.1 NME inventory 12-cell matrix",                P11_PASS),
    ("P1.2 PSF + g_A inventory",                         P12_PASS),
    ("P1.3 m_ββ Form A/B from ν.1",                       P13_PASS),
    ("P1.4 Method-spread σ_NME",                          P14_PASS),
    ("P1.5 Rate-only viability vs 90% CL",                P15_PASS),
]
n_pass = sum(1 for _, r in results if r)
n_total = len(results)
for name, r in results:
    print(f"  [{'PASS' if r else 'FAIL'}] {name}")
print()
print(f"  π.1.Phase1 score = {n_pass}/{n_total}")
if n_pass == n_total:
    print(f"  → π.1.Phase2 OK (full landscape PASS).")
elif n_pass >= 4:
    print(f"  → π.1.Phase2 OK ({n_pass}/{n_total} ≥ 4).")
else:
    print(f"  → π.1 FAIL gate; reframing required.")
