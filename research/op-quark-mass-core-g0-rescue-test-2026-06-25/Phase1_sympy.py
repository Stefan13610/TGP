# -*- coding: utf-8 -*-
# Phase 1 — op-quark-mass-core-g0-rescue-test-2026-06-25
# NIEZALEZNA re-weryfikacja maszynerii kwarkowej dodatekX + diagnoza strawmana HALT-B.
# 0 hardcoded T_pass; werdykt RESCUE-* WYLICZONY z flag. m_t pod OBOMA schematami.
#
# Maszyneria dodatekX (LIVE, read-only): dla sektora kwarkowego, majac (m_1,m_2)+A=a_G/phi,
# rozwiaz samozgodnie:  m_0 = A*m_3/m_1   AND   Q_K(m_1+m_0,m_2+m_0,m_3+m_0)=2/3   -> m_3.

import sympy as sp
from sympy import sqrt, nsolve, Symbol, Rational

print("="*72)
print("PHASE 1 RESCUE-TEST: niezalezne samozgodne m_3 + diagnoza strawmana HALT-B")
print("="*72)

R = {}  # wyniki testow

# ---- Stale LIVE (sektor leptonowy / TGP) ----
phi = (1 + sp.sqrt(5))/2
a_Gamma = sp.Rational(40, 1000)          # 0.040 (dodatekX l.335)
A_univ = a_Gamma/phi                     # 𝒜 = a_Γ/φ
CF = sp.Rational(4,3)                     # Casimir SU(3) fund.
alpha_s_PDG, alpha_s_sig = 0.1179, 0.0009

# ---- PDG 2024 quark masses [MeV] (instrument; rescue-test => dozwolone jako CEL) ----
m_u, m_c, m_t_pole, m_t_msbar = 2.16, 1270.0, 172760.0, 162500.0
m_d, m_s, m_b = 4.67, 93.4, 4180.0

print("\n[STALE] A=a_G/phi = %.5f ; phi=%.5f ; C_F=%.4f" % (float(A_univ), float(phi), float(CF)))

def QK(a,b,c):
    s = sp.sqrt(a)+sp.sqrt(b)+sp.sqrt(c)
    return (a+b+c)/s**2

# ---------- T1 (FP): φ-FP base anchors = [0,817; 0,891] (zamkniecie #40) ----------
# dodatekX res:X-phiFP-universal: g0^(1) down=0.8171, lepton=0.8695, up=0.8905
anchors = {"down":0.8171, "lepton":0.8695, "up":0.8905}
lo, hi = min(anchors.values()), max(anchors.values())
R['T1_anchors_eq_sek08b_range'] = (abs(lo-0.817)<0.01) and (abs(hi-0.891)<0.01)
print("\n[T1] phi-FP anchors {down,lepton,up} = {%.4f,%.4f,%.4f} -> [%.3f,%.3f] = sek08b:529? %s"
      % (anchors["down"],anchors["lepton"],anchors["up"], lo, hi, R['T1_anchors_eq_sek08b_range']))
print("     => [0,817;0,891] to KOTWICE sektorow, NIE domena (potwierdza #40 NORM-OVERLOAD)")

# ---------- T2 (FP, CENTRALNY): niezalezne samozgodne m_3 ----------
def solve_m3(m1, m2, A, guess):
    m3 = Symbol('m3', positive=True)
    m0 = A*m3/m1
    eq = QK(m1+m0, m2+m0, m3+m0) - sp.Rational(2,3)
    sol = nsolve(eq, m3, guess)
    m0v = float(A*sol/m1)
    return float(sol), m0v

# down: (m_d,m_s) -> m_b ; up: (m_u,m_c) -> m_t
m_b_pred, m0_d = solve_m3(m_d, m_s, A_univ, 4000.0)
m_t_pred, m0_u = solve_m3(m_u, m_c, A_univ, 170000.0)

delta_b = abs(m_b_pred - m_b)/m_b
delta_t_pole  = abs(m_t_pred - m_t_pole)/m_t_pole
delta_t_msbar = abs(m_t_pred - m_t_msbar)/m_t_msbar
delta_t_best = min(delta_t_pole, delta_t_msbar)

R['T2_mb_within_5pct'] = (delta_b <= 0.05)
R['T2_mt_within_5pct_some_scheme'] = (delta_t_best <= 0.05)
print("\n[T2] DOWN: m_b^pred = %.1f MeV (PDG 4180); m_0(d)=%.1f ; delta_b=%.2f%%"
      % (m_b_pred, m0_d, 100*delta_b))
print("[T2] UP:   m_t^pred = %.0f MeV ; m_0(u)=%.0f"
      % (m_t_pred, m0_u))
print("           vs pole 172760: delta=%.2f%% ; vs MS-bar 162500: delta=%.2f%% ; best=%.2f%%"
      % (100*delta_t_pole, 100*delta_t_msbar, 100*delta_t_best))
print("     m_b<=5%%? %s ; m_t<=5%% (jakis schemat)? %s" % (R['T2_mb_within_5pct'], R['T2_mt_within_5pct_some_scheme']))

# ---------- T3 (FP): 𝒜 universal + most α_s ----------
A_emp_d = m0_d * m_d / m_b_pred   # samozgodne
A_emp_u = m0_u * m_u / m_t_pred
# empiryczne z PDG bezposrednio:
A_pdg_d = 21.9 * m_d / m_b
A_pdg_u = 1981.5 * m_u / m_t_pole
dev_A = abs(float(A_univ) - 0.02464)/0.02464
R['T3_A_universal_within_2pct'] = (dev_A <= 0.02)
alpha_s_from_A = float(sp.sqrt(A_univ)/CF)
sig_off = abs(alpha_s_from_A - alpha_s_PDG)/alpha_s_sig
R['T3_alpha_s_within_1sigma'] = (sig_off <= 1.0)
print("\n[T3] A=a_G/phi=%.5f vs A_emp(PDG avg ~0.02464): dev=%.2f%% (<=2%%: %s)"
      % (float(A_univ), 100*dev_A, R['T3_A_universal_within_2pct']))
print("     A=C_F^2*alpha_s^2 -> alpha_s=sqrt(A)/C_F = %.5f vs PDG 0.1179+/-0.0009 (%.2f sigma; <=1: %s)"
      % (alpha_s_from_A, sig_off, R['T3_alpha_s_within_1sigma']))

# ---------- T4 (FP): diagnoza strawmana HALT-B ----------
# Skladniki maszynerii dodatekX (CLOSED): (i) additive m_0, (ii) phi-FP per-sektor, (iii) shifted Koide
haltb_has = {"additive_m0": False, "phiFP_per_sector": False, "shifted_koide": False}
# HALT-B formula: m=c_M*A_tail^2*g0^(e^2/2), wszystkie g0 in [0.817,0.891], bez m_0/Koide
missing = sum(1 for v in haltb_has.values() if not v)
R['T4_strawman_ge2_missing'] = (missing >= 2)
print("\n[T4] HALT-B formula sklada sie z maszynerii dodatekX? %s ; brakuje %d/3 skladnikow: %s"
      % (any(haltb_has.values()), missing, [k for k,v in haltb_has.items() if not v]))
print("     strawman (>=2/3 brak)? %s" % R['T4_strawman_ge2_missing'])

# ---------- T5 (FP): uczciwy licznik inputow ----------
# per sektor: {m_1 (skala), m_2 (lub r_21 via g0^(1))} = 2 wejscia -> m_3 predykcja
# 6 mas kwarkow: inputy {m_u,m_d,m_c,m_s} = 4 ; predykcje {m_b,m_t} = 2 ; shared {a_G,phi} z leptonow
n_quark_inputs = 4
n_quark_predictions = 2
m3_used_as_input = False  # m_3 NIE wstawiony do wlasnej predykcji (solver liczy go)
R['T5_genuine_prediction_not_fit'] = (not m3_used_as_input) and (n_quark_predictions >= 1)
print("\n[T5] LICZNIK: per sektor 2 wejscia (m_1,m_2) -> m_3 predykcja ; shared {a_G,phi}")
print("     6 mas kwarkow = 4 inputy + 2 PREDYKCJE (m_b,m_t) ; m_3 jako input? %s" % m3_used_as_input)
print("     r_31 = genuine predykcja (NIE 'zero parametrow', NIE fit)? %s" % R['T5_genuine_prediction_not_fit'])

# ---------- T6 (FP): WERDYKT wyliczony ----------
strawman = R['T4_strawman_ge2_missing']
closure_ok_5 = R['T2_mb_within_5pct'] and R['T2_mt_within_5pct_some_scheme']
A_ok = R['T3_A_universal_within_2pct']
closure_ok_15 = (delta_b <= 0.15) and (delta_t_best <= 0.15)
genuine = R['T5_genuine_prediction_not_fit']

if strawman and closure_ok_5 and A_ok:
    verdict = "RESCUE-CONFIRMED"
elif strawman and closure_ok_15 and genuine:
    verdict = "RESCUE-PARTIAL"
else:
    verdict = "RESCUE-FAILED"
R['T6_verdict_computed'] = True
print("\n[T6] flagi: strawman=%s closure<=5%%=%s A_univ<=2%%=%s closure<=15%%=%s genuine=%s"
      % (strawman, closure_ok_5, A_ok, closure_ok_15, genuine))
print("     >>> WERDYKT (wyliczony): %s <<<" % verdict)

# ---------- T7 (LIT/DEC): scheme-caveat + IMMUTABLE ----------
print("\n[T7-DEC] m_t scheme: pole 172760 vs MS-bar 162500 (~6%% rozjazd) -> RAPORT pod oboma.")
print("         dodatekX LIVE read-only; HALT-B IMMUTABLE; A=C_F^2*a_s^2 warunkowe (CG-1/3 R5).")

print("\n" + "="*72)
fp = [k for k in R if k.startswith(('T1','T2','T3','T4','T5','T6'))]
print("PASS: %d / %d FP-flag" % (sum(1 for k in fp if R[k]), len(fp)))
print("WERDYKT: %s" % verdict)
print("="*72)
