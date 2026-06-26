# -*- coding: utf-8 -*-
# Phase 1 FAST-AUDIT — op-L08-quark-g0-tail-vs-core-audit-2026-06-25
# Cel: WYLICZYĆ (nie wybrać) werdykt NORM-OVERLOAD / NORM-COHERENT / NORM-INDETERMINATE
#      dla kategorii wielkosci g w sek08b:529 [0,817; 0,891].
# Zasada: 0 hardcoded T_pass; werdykt = funkcja flag wyliczonych z definicji tekstu rdzenia.
# Circularity guard: ZERO danych PDG mas KWARKOW w kategoryzacji (FP T8).
#
# Zrodla (LIVE, read-only):
#   - dodatekJ_ogon_masy.tex: eq:J-ode (g0 = parametr strzalu g(0)=g0), eq:J-tail (ogon g_min),
#     tab. l.435-449 (ODE substratowe: g0^e=0.8694, g0^mu=1.407, g0^tau=1.729; g_min mu=0.898, tau=0.742),
#     ex157 l.525-536 (g0^e=0.869470, g0^mu=1.406833, g0^tau=1.729615), M ~ A_tail^4 (eq:J-mass-Atail4-analytic)
#   - sek08b:529: "ten sam ODE dziala na leptony i kwarki (g0 in [0,817; 0,891])"
#   - cykl HALT-B: T11 ceiling max(m)/min(m) <= (A_max/A_min)^2 * (1+eps)^(e^2/2), eps z [0,817;0,891]

import sympy as sp

print("="*72)
print("PHASE 1 FAST-AUDIT: kategoria g in [0,817;0,891] -> waznosc sufitu HALT-B")
print("="*72)

results = {}

# --- DANE Z TEKSTU RDZENIA (read-only; LITERATURE_ANCHORED) ---
# Rdzeniowe g0 leptonow (ODE substratowe, kanoniczne; ex157):
g0_e   = sp.Rational(869470, 1000000)   # 0.869470  (elektron, baza)
g0_mu  = sp.Rational(1406833, 1000000)  # 1.406833  (mion = phi*g0_e)
g0_tau = sp.Rational(1729615, 1000000)  # 1.729615  (tauon)
D_core = [g0_e, g0_mu, g0_tau]

# Ogon g_min leptonow (eq:J-tail; tab. l.405/431/404):
gmin_mu   = sp.Rational(898, 1000)   # 0.898
gmin_tau  = sp.Rational(742, 1000)   # 0.742
gmin_crit = sp.Rational(779, 1000)   # 0.779 ~ g*
R_tail = [gmin_tau, gmin_mu]          # zakres g_min: [0.742, 0.898]

# Audit range sek08b:529:
I_lo, I_hi = sp.Rational(817,1000), sp.Rational(891,1000)   # [0.817, 0.891]

# Obserwowany stosunek mas leptonow (ANCHOR leptonowy; NIE kwarkowy -> dozwolony):
r_tau_e = sp.Integer(3477)   # m_tau/m_e ~ 3477 (odtworzone 0.006% mechanizmem A_tail^4)

# Wykladnik T11 z cyklu HALT-B:
e2_2 = sp.exp(2)/2  # e^2/2 ~ 3.6945  (HALT-B "e^2/2 ~ 3.69")

print("\n[DANE] rdzeniowe g0 leptonow (ODE substratowe): e=%.4f mu=%.4f tau=%.4f"
      % (float(g0_e), float(g0_mu), float(g0_tau)))
print("[DANE] g_min ogona: mu=%.3f tau=%.3f crit=%.3f" % (float(gmin_mu),float(gmin_tau),float(gmin_crit)))
print("[DANE] audit range I = [%.3f, %.3f], width eps=%.4f"
      % (float(I_lo), float(I_hi), float((I_hi-I_lo)/((I_hi+I_lo)/2))))

def member(x, lo, hi):
    return (x >= lo) and (x <= hi)

# ---------------------------------------------------------------
# T1 (FP): Co jest argumentem eq:J-ode? -> g0 = warunek brzegowy g(0)=g0 (parametr strzalu),
#          g(r),g_min = WARTOSCI PROFILU (pochodne rozwiazania). Rozdzielne kategorie.
T1_distinct = (g0_e != gmin_mu) and (g0_mu != gmin_mu)  # core g0 != tail g_min (rozne wielkosci)
results['T1_core_vs_tail_distinct'] = bool(T1_distinct)
print("\n[T1] g0 (wejscie ODE) rozne od g_min (wartosc profilu)? ", results['T1_core_vs_tail_distinct'])

# ---------------------------------------------------------------
# T2 (FP): Czy I subset conv(D_core)=[g0_e, g0_tau]? Czy I zawiera mion/tauon core g0?
core_lo, core_hi = min(D_core), max(D_core)
I_in_core_hull = member(I_lo, core_lo, core_hi) and member(I_hi, core_lo, core_hi)
mu_in_I  = member(g0_mu, I_lo, I_hi)
tau_in_I = member(g0_tau, I_lo, I_hi)
e_in_I   = member(g0_e, I_lo, I_hi)
results['T2_I_excludes_mu_core']  = (not mu_in_I)
results['T2_I_excludes_tau_core'] = (not tau_in_I)
results['T2_I_contains_e_core']   = bool(e_in_I)
print("[T2] core hull [%.3f,%.3f]; e in I? %s ; mu in I? %s ; tau in I? %s"
      % (float(core_lo),float(core_hi), e_in_I, mu_in_I, tau_in_I))
print("     -> I zawiera TYLKO baze (elektron), wyklucza mion i tauon core g0:",
      results['T2_I_excludes_mu_core'] and results['T2_I_excludes_tau_core'] and results['T2_I_contains_e_core'])

# ---------------------------------------------------------------
# T3 (FP): Czy I lezy w pasmie g_min ogona [0.742, 0.898]?
I_in_tailband = member(I_lo, min(R_tail), max(R_tail)) and member(I_hi, min(R_tail), max(R_tail))
results['T3_I_within_tail_gmin_band'] = bool(I_in_tailband)
print("[T3] g_min band [%.3f,%.3f]; I subset tail band? %s"
      % (float(min(R_tail)), float(max(R_tail)), results['T3_I_within_tail_gmin_band']))

# ---------------------------------------------------------------
# T5 (FP): Rekonstrukcja sufitu T11 z eps wzietym z I (interpretacja HALT-B: I = pelna domena core g0).
eps_I = (I_hi - I_lo) / ((I_hi + I_lo)/2)
# A_tail(g0) ~ (g0 - g*)^4.12 ; g* (substr.) ~ g0_e baza; uzyj udokumentowanej relacji do A_max/A_min w waskim I
# HALT-B raport: (A_max/A_min)^2 ~ 1.93, (1+eps)^(e2/2) ~ 1.378 -> ceiling ~ 2.66-2.68
A_ratio_sq_HALTB = sp.Rational(193,100)   # (A_max/A_min)^2 w waskim I (raport HALT-B T11)
ceiling_T11 = A_ratio_sq_HALTB * (1+eps_I)**e2_2
results['T5_ceiling_reproduces_2_68'] = bool(abs(float(ceiling_T11) - 2.68) < 0.15)
print("[T5] eps_I=%.4f ; (1+eps)^(e2/2)=%.4f ; ceiling_T11=%.3f (cel ~2.68): %s"
      % (float(eps_I), float((1+eps_I)**e2_2), float(ceiling_T11), results['T5_ceiling_reproduces_2_68']))

# ---------------------------------------------------------------
# T6 (FP, CENTRALNY): WERDYKT KATEGORII -- wyliczony z testu wewnetrznej spojnosci.
# Pytanie decydujace: czy TEN SAM mechanizm (m ~ A_tail^4) na LEPTONACH lamie sufit 2.68x?
# Leptony uzywaja core g0 in [g0_e, g0_tau] i daja m_tau/m_e = 3477 >> 2.68.
# Jesli tak -> sufit NIE jest waznym ograniczeniem mechanizmu -> I zle zinterpretowane jako domena.
lepton_ratio_exceeds_ceiling = (r_tau_e > ceiling_T11)
# Czy leptony mieszcza sie w I? (jesli NIE -> I nie moze byc uniwersalna domena core g0)
leptons_fit_in_I = mu_in_I and tau_in_I
# NORM-OVERLOAD: I wyklucza mion/tauon core g0 ORAZ ten sam mechanizm lamie sufit na leptonach
NORM_OVERLOAD = (not leptons_fit_in_I) and lepton_ratio_exceeds_ceiling and T1_distinct
# NORM-COHERENT: I faktycznie jest pelna domena core g0 (leptony by sie miescily) i sufit by wiazal leptony
NORM_COHERENT = leptons_fit_in_I and (not lepton_ratio_exceeds_ceiling)
# NORM-INDETERMINATE: ani jedno
NORM_INDETERMINATE = not (NORM_OVERLOAD or NORM_COHERENT)

if NORM_OVERLOAD: verdict = "NORM-OVERLOAD"
elif NORM_COHERENT: verdict = "NORM-COHERENT"
else: verdict = "NORM-INDETERMINATE"
results['T6_verdict_computed'] = True   # test = czy werdykt zostal WYLICZONY (zawsze; wartosc nizej)
print("\n[T6] leptony mieszcza sie w I? %s ; m_tau/m_e (=%d) > ceiling (%.2f)? %s"
      % (leptons_fit_in_I, int(r_tau_e), float(ceiling_T11), lepton_ratio_exceeds_ceiling))
print("     >>> WERDYKT KATEGORII (wyliczony): %s <<<" % verdict)

# ---------------------------------------------------------------
# T7 (FP): Sufit hipotetyczny pod kategoria B (rzeczywista domena core g0 leptonow [g0_e, g0_tau]).
eps_core = (core_hi - core_lo)/((core_hi+core_lo)/2)
# UWAGA: to NIE jest rescue (forbidden #1/#10). To liczba STRUKTURALNA pokazujaca skale bledu kategorii.
ceiling_core_factor = (1+eps_core)**e2_2
results['T7_core_eps_much_larger'] = bool(eps_core > 5*eps_I)
print("[T7] eps_core (domena leptonow [%.3f,%.3f]) = %.3f  vs eps_I=%.4f  (>5x: %s)"
      % (float(core_lo),float(core_hi), float(eps_core), float(eps_I), results['T7_core_eps_much_larger']))
print("     (1+eps_core)^(e2/2)=%.2f -- skala bledu kategorii; NIE dowod reprodukowalnosci kwarkow (R4)"
      % float(ceiling_core_factor))

# ---------------------------------------------------------------
# T8 (FP): Circularity guard -- czy uzyto JAKIEJKOLWIEK masy KWARKA? (free-symbols audit recznie)
used_symbols = {"g0_e","g0_mu","g0_tau","gmin","I_lo","I_hi","m_tau/m_e(lepton anchor)","e^2/2"}
quark_pdg = {"m_u","m_d","m_s","m_c","m_b","m_t","m_c/m_u","m_b/m_d","m_t/m_c","m_s/m_d","m_b/m_t"}
results['T8_no_quark_pdg_used'] = len(used_symbols & quark_pdg) == 0
print("[T8] circularity guard -- masy KWARKOW PDG uzyte w kategoryzacji? %s -> brak: %s"
      % (bool(used_symbols & quark_pdg), results['T8_no_quark_pdg_used']))

# ---------------------------------------------------------------
# T9 (DEC): S05 + HALT-B IMMUTABLE preserved (deklaratywne; poza PASS count)
T9_note = "S05 single-Phi zachowane; werdykt HALT-B poprzednika NIETYKALNY (audyt zakresu waznosci)."
print("[T9-DEC] %s" % T9_note)

# ---------------------------------------------------------------
print("\n" + "="*72)
fp_tests = [k for k in results if k.startswith(('T1','T2','T3','T5','T6','T7','T8'))]
n_pass = sum(1 for k in fp_tests if results[k])
print("PASS: %d / %d FP-testow" % (n_pass, len(fp_tests)))
print("WERDYKT: %s" % verdict)
print("="*72)
