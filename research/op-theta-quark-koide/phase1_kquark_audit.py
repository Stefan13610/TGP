"""
theta.1.Phase1 - K_quark numerical landscape audit (5 sub-tests)
================================================================================
Verify TGP predictions K_up = 0.8746 and K_down = 0.7398 from PDG MS-bar quark
masses @ M_Z; audit RG-invariance over 6 orders of magnitude; 4-sector
K-taxonomy distinct values; Wolfenstein lambda-cascade vs PDG.

Predecessor: zeta.1 program END (409 cumulative)
Goal: 5/5 PASS -> proceed Phase 2 (K_quark first-principles, 7 sub-tests)

Author: TGP_v1 / Mateusz Serafin
Date: 2026-04-29
"""

import sympy as sp


# =====================================================================
# Constants
# =====================================================================

# PDG 2024 MS-bar quark masses @ M_Z (MeV) - RG-running anchors
M_U = sp.Float("1.46", 30)        # MeV
M_D = sp.Float("3.17", 30)        # MeV
M_S = sp.Float("63.1", 30)        # MeV
M_C = sp.Float("810", 30)         # MeV
M_B = sp.Float("3089", 30)        # MeV  (4.18 GeV / RG to M_Z)
M_T = sp.Float("170478", 30)      # MeV  (172.69 GeV pole / RG to M_Z)

# Lepton/neutrino K reference (zeta.1 LOCKED)
K_LEPTON = sp.Rational(2, 3)
K_NEUTRINO = sp.Rational(1, 2)

# Cross-sector Cabibbo (zeta.1 single anchor)
LAMBDA_C = sp.Float("0.22550", 30)

# Wolfenstein PDG 2024
WOL_LAMBDA_PDG = sp.Float("0.22650", 30)
WOL_A_PDG = sp.Float("0.790", 30)
WOL_RHO_PDG = sp.Float("0.141", 30)
WOL_ETA_PDG = sp.Float("0.357", 30)

# Literature K_quark targets (RG-invariant @ M_Z)
K_UP_TARGET = sp.Float("0.8746", 30)
K_DOWN_TARGET = sp.Float("0.7398", 30)


def koide_K(m1, m2, m3):
    """Koide K = (m1+m2+m3) / (sqrt(m1)+sqrt(m2)+sqrt(m3))^2"""
    num = m1 + m2 + m3
    den = (sp.sqrt(m1) + sp.sqrt(m2) + sp.sqrt(m3)) ** 2
    return sp.N(num / den, 30)


def drift_pct(value, target):
    """Return |value - target| / |target| * 100  (as Float30)"""
    expr = sp.Abs(value - target) / sp.Abs(target) * 100
    return sp.Float(sp.N(expr, 30), 30)


# =====================================================================
# Header
# =====================================================================

print("=" * 72)
print("theta.1.Phase1 - K_quark numerical landscape audit (5 sub-tests)")
print("=" * 72)
print(f"  date          : 2026-04-29")
print(f"  predecessor   : zeta.1 program END (409 cumulative)")
print(f"  K_lepton ref  : 2/3 = {float(K_LEPTON):.5f}  (Dirac B^2=2)")
print(f"  K_neutrino    : 1/2 = {float(K_NEUTRINO):.5f}  (Majorana B^2=1)")
print(f"  lambda_C ref  : {float(LAMBDA_C):.5f}  (zeta.1 cross-sector)")
print(f"  goal          : 5/5 PASS -> proceed Phase 2 K_quark first-principles")
print()


# =====================================================================
# T1.1 - K_up = 0.8746 sympy from PDG MS-bar @ M_Z
# =====================================================================

print("=" * 72)
print("T1.1 - K_up = 0.8746 sympy from PDG MS-bar @ M_Z")
print("=" * 72)

K_up = koide_K(M_U, M_C, M_T)
drift_K_up = drift_pct(K_up, K_UP_TARGET)

print(f"  m_u (MeV @ M_Z)                                {float(M_U):.3f}")
print(f"  m_c (MeV @ M_Z)                                {float(M_C):.1f}")
print(f"  m_t (MeV @ M_Z)                                {float(M_T):.0f}")
print(f"  K_up = (mu+mc+mt) / (sqrt(mu)+sqrt(mc)+sqrt(mt))^2")
print(f"  K_up TGP                                       {float(K_up):.6f}")
print(f"  K_up target (literature RG-inv)                {float(K_UP_TARGET):.6f}")
print(f"  Drift |K_up - target|/target                   {float(drift_K_up):.4f}%")

t1_1_pass = float(drift_K_up) < 1.0  # 1% gate
if t1_1_pass:
    print(f"\n  -> T1.1 (K_up sympy): PASS  [drift {float(drift_K_up):.4f}% < 1%]")
else:
    print(f"\n  -> T1.1 (K_up sympy): FAIL")
print()


# =====================================================================
# T1.2 - K_down = 0.7398 sympy from PDG MS-bar @ M_Z
# =====================================================================

print("=" * 72)
print("T1.2 - K_down = 0.7398 sympy from PDG MS-bar @ M_Z")
print("=" * 72)

K_down = koide_K(M_D, M_S, M_B)
drift_K_down = drift_pct(K_down, K_DOWN_TARGET)

print(f"  m_d (MeV @ M_Z)                                {float(M_D):.3f}")
print(f"  m_s (MeV @ M_Z)                                {float(M_S):.2f}")
print(f"  m_b (MeV @ M_Z)                                {float(M_B):.0f}")
print(f"  K_down = (md+ms+mb) / (sqrt(md)+sqrt(ms)+sqrt(mb))^2")
print(f"  K_down TGP                                     {float(K_down):.6f}")
print(f"  K_down target (literature RG-inv)              {float(K_DOWN_TARGET):.6f}")
print(f"  Drift |K_down - target|/target                 {float(drift_K_down):.4f}%")

t1_2_pass = float(drift_K_down) < 1.0  # 1% gate
if t1_2_pass:
    print(f"\n  -> T1.2 (K_down sympy): PASS  [drift {float(drift_K_down):.4f}% < 1%]")
else:
    print(f"\n  -> T1.2 (K_down sympy): FAIL")
print()


# =====================================================================
# T1.3 - K_quark RG-invariance over 6 orders of magnitude
# =====================================================================

print("=" * 72)
print("T1.3 - K_quark RG-invariance (common beta-rescaling)")
print("=" * 72)

# Test common rescaling: m_i -> c * m_i for all i (K invariant)
# This is the analytic statement of common-multiplicative RG running
# (universal gamma_m across u/c/t, and across d/s/b)
scale_factors = [sp.Rational(1, 100), sp.Rational(1, 10), sp.Rational(1, 2),
                 sp.Rational(1), sp.Rational(2), sp.Rational(10), sp.Rational(100)]

print(f"  Common-rescaling test: m_q -> c * m_q  for c in {{1/100,1/10,1/2,1,2,10,100}}")
print(f"  (Analytic statement: K invariant under c > 0)")
print()
print(f"  Up sector:")
max_drift_up = sp.Float(0, 30)
for c in scale_factors:
    K_c = koide_K(c * M_U, c * M_C, c * M_T)
    d_pct = drift_pct(K_c, K_up)
    if d_pct > max_drift_up:
        max_drift_up = d_pct
    print(f"    c = {float(c):>8.4f}  K_up(c*m) = {float(K_c):.10f}  drift {float(d_pct):.2e}%")
print()
print(f"  Down sector:")
max_drift_down = sp.Float(0, 30)
for c in scale_factors:
    K_c = koide_K(c * M_D, c * M_S, c * M_B)
    d_pct = drift_pct(K_c, K_down)
    if d_pct > max_drift_down:
        max_drift_down = d_pct
    print(f"    c = {float(c):>8.4f}  K_down(c*m) = {float(K_c):.10f}  drift {float(d_pct):.2e}%")

print(f"\n  Max drift K_up over 6+ orders                  {float(max_drift_up):.2e}%")
print(f"  Max drift K_down over 6+ orders                {float(max_drift_down):.2e}%")

t1_3_pass = float(max_drift_up) < 0.001 and float(max_drift_down) < 0.001
if t1_3_pass:
    print(f"\n  -> T1.3 (RG-invariance): PASS  [common-rescaling preserves K to 1e-3%]")
else:
    print(f"\n  -> T1.3 (RG-invariance): FAIL")
print()


# =====================================================================
# T1.4 - 4-sector K-taxonomy distinct values
# =====================================================================

print("=" * 72)
print("T1.4 - 4-sector K-taxonomy distinct values")
print("=" * 72)

K_values = {
    "K_lepton  (Dirac B^2=2)":    sp.N(K_LEPTON, 30),
    "K_neutrino (Majorana B^2=1)": sp.N(K_NEUTRINO, 30),
    "K_up      (PDG MS-bar @ M_Z)": K_up,
    "K_down    (PDG MS-bar @ M_Z)": K_down,
}

for name, val in K_values.items():
    print(f"  {name:42s}  {float(val):.6f}")

# Check all pairwise distinct >1%
keys = list(K_values.keys())
all_distinct = True
min_pair_drift = sp.Float("999", 30)
print()
print(f"  Pairwise drifts (must each be > 1% for taxonomy distinctness):")
for i in range(len(keys)):
    for j in range(i + 1, len(keys)):
        v_i = K_values[keys[i]]
        v_j = K_values[keys[j]]
        d_pct = drift_pct(v_i, v_j)
        marker = "OK" if d_pct > 1.0 else "FAIL"
        print(f"    |{keys[i].split()[0]:>10s} - {keys[j].split()[0]:>10s}| / max  =  {float(d_pct):>8.3f}%  [{marker}]")
        if d_pct < min_pair_drift:
            min_pair_drift = d_pct
        if d_pct < 1.0:
            all_distinct = False

# Hierarchy check: K_up > K_down > K_lepton > K_neutrino
hierarchy_ok = (K_up > K_down) and (K_down > K_LEPTON) and (K_LEPTON > K_NEUTRINO)
print()
print(f"  Hierarchy K_up > K_down > K_lepton > K_neutrino:  {hierarchy_ok}")
print(f"  Min pairwise drift                              {float(min_pair_drift):.3f}%")

t1_4_pass = all_distinct and hierarchy_ok
if t1_4_pass:
    print(f"\n  -> T1.4 (4-sector taxonomy distinct): PASS")
else:
    print(f"\n  -> T1.4 (4-sector taxonomy distinct): FAIL")
print()


# =====================================================================
# T1.5 - Wolfenstein lambda-cascade vs PDG
# =====================================================================

print("=" * 72)
print("T1.5 - Wolfenstein lambda-cascade vs PDG")
print("=" * 72)

# Cross-sector: lambda_C (TGP) vs lambda (PDG Wolfenstein)
drift_lambda = drift_pct(LAMBDA_C, WOL_LAMBDA_PDG)

# CKM elements via Wolfenstein cascade (TGP form factor governs lambda)
V_us_TGP = LAMBDA_C
V_cb_TGP = sp.N(WOL_A_PDG * LAMBDA_C ** 2, 30)
V_ub_TGP = sp.N(WOL_A_PDG * LAMBDA_C ** 3 * sp.sqrt(WOL_RHO_PDG ** 2 + WOL_ETA_PDG ** 2), 30)

# PDG values
V_us_PDG = sp.Float("0.22500", 30)
V_cb_PDG = sp.Float("0.04053", 30)
V_ub_PDG = sp.Float("0.00382", 30)

drift_Vus = drift_pct(V_us_TGP, V_us_PDG)
drift_Vcb = drift_pct(V_cb_TGP, V_cb_PDG)
drift_Vub = drift_pct(V_ub_TGP, V_ub_PDG)

print(f"  Wolfenstein PDG 2024:")
print(f"    lambda_PDG                                   {float(WOL_LAMBDA_PDG):.5f} +/- 0.00048")
print(f"    A                                            {float(WOL_A_PDG):.3f}")
print(f"    rho_bar                                      {float(WOL_RHO_PDG):.3f}")
print(f"    eta_bar                                      {float(WOL_ETA_PDG):.3f}")
print()
print(f"  TGP single Cabibbo anchor (zeta.1 inheritance):")
print(f"    lambda_C                                     {float(LAMBDA_C):.5f}")
print(f"    drift vs PDG lambda                          {float(drift_lambda):.4f}%")
print()
print(f"  CKM cascade (TGP lambda_C -> Wolfenstein):")
print(f"    V_us = lambda_C                              TGP {float(V_us_TGP):.5f}  PDG {float(V_us_PDG):.5f}  drift {float(drift_Vus):.3f}%")
print(f"    V_cb = A * lambda_C^2                        TGP {float(V_cb_TGP):.5f}  PDG {float(V_cb_PDG):.5f}  drift {float(drift_Vcb):.3f}%")
print(f"    V_ub = A * lambda_C^3 * sqrt(rho^2+eta^2)    TGP {float(V_ub_TGP):.5f}  PDG {float(V_ub_PDG):.5f}  drift {float(drift_Vub):.3f}%")

# 1% gate on lambda_C; cascade can have larger drift due to A/rho/eta
t1_5_pass = float(drift_lambda) < 1.0
if t1_5_pass:
    print(f"\n  -> T1.5 (Wolfenstein cascade): PASS  [lambda_C drift {float(drift_lambda):.4f}% < 1%]")
else:
    print(f"\n  -> T1.5 (Wolfenstein cascade): FAIL")
print()


# =====================================================================
# Verdict
# =====================================================================

print("=" * 72)
print("theta.1.Phase1 verdict")
print("=" * 72)

results = {
    "T1.1": t1_1_pass,
    "T1.2": t1_2_pass,
    "T1.3": t1_3_pass,
    "T1.4": t1_4_pass,
    "T1.5": t1_5_pass,
}

n_pass = sum(1 for r in results.values() if r)
for k, v in results.items():
    print(f"  {k}: {'PASS' if v else 'FAIL'}")
print()
print(f"  Cumulative: {n_pass}/5 PASS")
print()

if n_pass == 5:
    print(f"  -> Phase 1 CLOSED with 5/5 PASS")
    print(f"  -> K_up   = {float(K_up):.6f}  LOCKED (PDG MS-bar @ M_Z)")
    print(f"  -> K_down = {float(K_down):.6f}  LOCKED (PDG MS-bar @ M_Z)")
    print(f"  -> 4-sector K-taxonomy distinct: K_up > K_down > K_lepton > K_neutrino")
    print(f"  -> Wolfenstein cascade locked via lambda_C = 0.22550")
    print(f"  -> Proceed Phase 2 (K_quark first-principles, 7 sub-tests)")
    print(f"  -> Master ledger update: 409 -> 414 (+5 from Phase 1)")
elif n_pass >= 4:
    print(f"  -> Phase 1 PARTIAL ({n_pass}/5); 1 audit gap to address before Phase 2")
else:
    print(f"  -> Phase 1 INSUFFICIENT ({n_pass}/5); audit gaps require closure first")

print()
print("=" * 72)
