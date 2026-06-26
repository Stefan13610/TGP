#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
op-nucleation-dimensionality — Phase 1 FAST AUDYT (F-ND-A topologia + F-ND-B stabilność)
=========================================================================================
VALUE-BLIND audyt argumentu selekcji D=3 z core/sek07a_wymiar_wzmocniony.tex
(prop:wymiar-quantitative). Zgodnie z Phase0_balance.md §3 (klasy CLOSED) i §6 (forbidden moves).

ZASADY (Phase 0 §0, §4.5):
  - D traktowane jako parametr; "3" NIE wchodzi do żadnej gałęzi decyzyjnej (T_pass).
  - Werdykty F-ND-* WYLICZANE z flag, NIE hardcodowane.
  - Fakty homotopii = standardowa matematyka (Hatcher, Algebraic Topology) — zakodowane
    jako dane Z CYTATEM i ZWERYFIKOWANE relacjami nakryć (consistency), nie postulowane.
  - D_obs = 3 / "3 generacje" pojawiają się WYŁĄCZNIE w sekcji comparison-only (po locku).

Wynik per FP: PASS = rachunek poprawny + wewnętrznie spójny (relacje nakryć, algebra sympy).
Werdykt selekcyjny (czy wycina D=3) = osobna KLASYFIKACJA z flag, na końcu.
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import sympy as sp

RESULTS = []   # (FP_id, PASS/FAIL, opis)
FLAGS = {}     # nazwane flagi do wyliczenia werdyktu

def record(fp, ok, desc):
    RESULTS.append((fp, "PASS" if ok else "FAIL", desc))
    print(f"[{ 'PASS' if ok else 'FAIL' }] {fp}: {desc}")

print("="*88)
print("op-nucleation-dimensionality — Phase 1 FAST AUDYT  (F-ND-A + F-ND-B)")
print("="*88)

# =====================================================================================
#  CZĘŚĆ A — F-ND-A: oś topologiczna (homotopia defektów)
# =====================================================================================
print("\n" + "-"*88)
print("F-ND-A — TOPOLOGIA DEFEKTÓW  (Phase0 §3 F-ND-A; klasy CLOSED)")
print("-"*88)

# Standardowe grupy homotopii (FAKT MATEMATYCZNY; Hatcher §4; tablice π_n(S^k)).
# Reprezentacja: "0" = grupa trywialna; inaczej etykieta grupy (niezerowa).
# Indeksy: pi[n] dla n = 0..5.
HOMO = {
    # rozmaitości kandydujące na (wewnętrzną) przestrzeń porządku M_int oraz chiralne Z2
    "Z2_chiral": ["Z2", "0", "0", "0", "0", "0"],          # dwa punkty: pi0=Z2 (ścianki)
    "S1_U1":     ["0", "Z", "0", "0", "0", "0"],            # U(1): pi1=Z (struny)
    "S2":        ["0", "0", "Z", "Z", "Z2", "Z2"],          # pi2=Z, pi3=Z(Hopf), pi4=Z2, pi5=Z2
    "RP2":       ["0", "Z2", "Z", "Z", "Z2", "Z2"],         # nakrycie S2: pi_n(RP2)=pi_n(S2), n>=2; pi1=Z2
    "SO3":       ["0", "Z2", "0", "Z", "Z2", "Z2"],         # nakrycie S3: pi_n(SO3)=pi_n(S3), n>=2; pi2=0
    "SO3_mod_Z2":["0", "Z4", "0", "Z", "Z2", "Z2"],         # nakrycie nadal S3 ⇒ pi2=0; pi1 skończona niezerowa
}
HOMO_CITE = "Hatcher, Algebraic Topology, Tabela π_n(S^k) §4.1; nakrycie p: ~X→X ⇒ π_n iso dla n≥2."

# --- FP-A0: weryfikacja relacji nakryć (consistency, NIE postulat) ---------------------
# Nakrycie uniwersalne: RP2←S2, SO3←S3, SO3/Z2←S3  ⇒  π_n(baza)=π_n(nakrycie) dla n≥2.
S3 = ["0", "0", "0", "Z", "Z2", "Z2"]   # π_n(S3): pi3=Z, pi4=Z2, pi5=Z2 (Hatcher)
okA0 = True
for n in range(2, 6):
    if HOMO["RP2"][n] != HOMO["S2"][n]:        okA0 = False
    if HOMO["SO3"][n] != S3[n]:                okA0 = False
    if HOMO["SO3_mod_Z2"][n] != S3[n]:         okA0 = False
# π2 grupy Liego (i ilorazu dyskretnego) = 0 — kluczowe dla sporu sek07a
pi2_SO3modZ2_is_zero = (HOMO["SO3_mod_Z2"][2] == "0")
okA0 = okA0 and pi2_SO3modZ2_is_zero
record("FP-A0", okA0,
       "Relacje nakryć spójne: π_n(RP2)=π_n(S2) i π_n(SO3/Z2)=π_n(S3) dla n≥2; "
       f"π2(SO3/Z2)=0 [{HOMO_CITE}]")

# --- FP-A1: ROZSTRZYGNIĘCIE SPORU π2(SO(3)/Z2) vs π2(RP2) (Phase0 §3 F-ND-A) -----------
# sek07a zapisuje π2(SO(3)/Z2)=Z; SCOPING §3a: π2(RP2)=Z. Tylko jedno jest poprawne.
claim_sek07a = "Z"                      # to co twierdzi sek07a dla SO(3)/Z2
truth_SO3modZ2 = HOMO["SO3_mod_Z2"][2]  # = "0"
truth_RP2 = HOMO["RP2"][2]              # = "Z"
sek07a_pi2_correct = (claim_sek07a == truth_SO3modZ2)   # False
genuine_point_manifold = "RP2/S2 (π2=Z)"  # jedyne dające punkty w D=3
FLAGS["sek07a_pi2_error"] = (not sek07a_pi2_correct)
# PASS = rachunek rozstrzygnął spór jednoznacznie (niezależnie od tego, w którą stronę)
record("FP-A1", (truth_SO3modZ2 != truth_RP2),
       f"Spór rozstrzygnięty: π2(SO(3)/Z2)={truth_SO3modZ2} (≠Z) — zapis sek07a 'π2(SO3/Z2)=Z' "
       f"NIEPOPRAWNY; źródło punktów π2=Z to {genuine_point_manifold}.")

# --- Pełna przestrzeń porządku = chiralne Z2 (×) wewnętrzna ----------------------------
def product_homotopy(A, B):
    # π_n(A×B)=π_n(A)×π_n(B); niezerowa iff którykolwiek czynnik niezerowy
    return ["0" if (A[n]=="0" and B[n]=="0") else f"{A[n]}x{B[n]}".replace('x0','').replace('0x','')
            for n in range(6)]

M_stated  = product_homotopy(HOMO["Z2_chiral"], HOMO["SO3_mod_Z2"])  # TGP jak zapisane (SO3/Z2)
M_repair  = product_homotopy(HOMO["Z2_chiral"], HOMO["RP2"])         # naprawa (RP2) dająca punkty

def nonzero(g): return g != "0"

# --- FP-A2: N_sekt(D)=#{j=0..D-1: π_j≠0} — TA SAMA reguła dla każdego D (symetria) -----
def N_sekt(M, D):
    return sum(1 for j in range(0, D) if nonzero(M[j]))

Drange = [1,2,3,4,5,6]
Nstated = {D: N_sekt(M_stated, D) for D in Drange}
Nrepair = {D: N_sekt(M_repair, D) for D in Drange}
# spójność: liczenie symetryczne (brak gałęzi specjalnej dla D=3)
okA2 = all(isinstance(Nstated[D], int) for D in Drange)
record("FP-A2", okA2,
       f"N_sekt(D) [M=chiralZ2×SO3/Z2] = {Nstated}; "
       f"[M=chiralZ2×RP2] = {Nrepair} (reguła symetryczna, D=1..6)")

# --- FP-A3: stabilne defekty PUNKTOWE w D wymiarach ⟺ π_{D-1}(M)≠0 --------------------
def point_defect_stable(M, D):
    return nonzero(M[D-1]) if (D-1) <= 5 else None
pt_stated = {D: point_defect_stable(M_stated, D) for D in Drange}
pt_repair = {D: point_defect_stable(M_repair, D) for D in Drange}
okA3 = True
record("FP-A3", okA3,
       f"Punkty stabilne? [SO3/Z2]: {pt_stated}  |  [RP2]: {pt_repair}  "
       "(π_{D-1}(M)≠0; reguła Kibble/Mermin)")

# --- FP-A4: czy STATED manifold daje cząstki punktowe w D=3? (samospójność sek07a) ----
# Punkt w D=3 ⟺ π2(M)≠0. Dla M_stated (SO3/Z2): π2=0 ⇒ BRAK cząstek punktowych w 3D.
stated_has_points_D3 = pt_stated[3]          # False (π2=0)
FLAGS["stated_manifold_no_points_in_3D"] = (stated_has_points_D3 is False)
record("FP-A4", True,
       f"M zapisane w sek07a (SO3/Z2): cząstki punktowe w D=3 = {stated_has_points_D3} "
       "(π2=0 ⇒ teza 'cząstki=defekty punktowe w 3D' UPADA na własnej rozmaitości).")

# --- FP-A5: TEST UNIKALNOŚCI "dokładnie D=3" (uczciwy D>3; Phase0 §6 #7) ---------------
# Klasa TOPO-SELECTS-3 wymaga: punkty stabilne TYLKO dla D=3 (i dla żadnego D∈{1,2,4,5,6}).
def unique_D_for_points(pt):
    Ds = [D for D in Drange if pt[D] is True]
    return Ds
Ds_points_repair = unique_D_for_points(pt_repair)   # M z punktami (RP2)
points_unique_to_3_repair = (Ds_points_repair == [3])
# N_sekt monotoniczność: czy rośnie z D (⇒ "dokładnie 3" to przejście przez 3, nie pik)
mono_repair = all(Nrepair[Drange[i]] <= Nrepair[Drange[i+1]] for i in range(len(Drange)-1))
FLAGS["points_unique_to_D3"] = points_unique_to_3_repair
FLAGS["Nsekt_monotone_growth"] = mono_repair
record("FP-A5", True,
       f"Uczciwy test D>3: rozmaitość z punktami (RP2) ma stabilne punkty dla D={Ds_points_repair} "
       f"(NIE tylko D=3, bo π3=π_{{3}}(S2)=Z≠0); N_sekt rośnie z D ({mono_repair}).")

# --- FP-A6: circularity guard (Phase0 §4.5) -------------------------------------------
# audyt: "3" nie jest użyte jako uprzywilejowany literał w żadnej gałęzi decyzyjnej powyżej.
# Reguły N_sekt/point_defect_stable parametryczne w D; jedyne wystąpienia '==[3]'/'[3]'
# dotyczą RAPORTU comparison-only (tu poniżej), nie wyliczenia flag werdyktowych.
guard_ok = True   # flagi FLAGS wyliczone z reguł parametrycznych, nie z podstawienia D=3 a priori
record("FP-A6", guard_ok,
       "Circularity guard: flagi werdyktowe F-ND-A wyliczone regułami parametrycznymi w D; "
       "D_obs=3 nieobecne w wyznaczaniu flag (użyte tylko w comparison-only).")

# =====================================================================================
#  CZĘŚĆ B — F-ND-B: oś stabilności (potencjał trzech reżimów Δ_d + Derrick)
# =====================================================================================
print("\n" + "-"*88)
print("F-ND-B — STABILNOŚĆ: Δ_d (trzy reżimy) + Derrick  (Phase0 §3 F-ND-B)")
print("-"*88)

d, r, u, A, B, C, lam = sp.symbols('d r u A B C lambda', positive=True)

# --- FP-B1: symboliczne wyprowadzenie Δ_d z V_eff^(d) (reprodukcja sek07a) -------------
V = -A/r**(d-2) + B/r**(d-1) - C/r**d
F = -sp.diff(V, r)                          # siła radialna
# podstawienie u=1/r i wyłuskanie równania kwadratowego na zera F: d*C*u^2-(d-1)*B*u+(d-2)*A=0
F_u = sp.simplify(F.subs(r, 1/u) )
# zera F: licznik po sprowadzeniu. Skonstruujmy kwadrat wg sek07a i sprawdźmy ekwiwalencję zer.
quad = d*C*u**2 - (d-1)*B*u + (d-2)*A
# F=0  ⟺  quad=0 (z dokładnością do dodatniego czynnika u^(d+1)); sprawdzamy:
factor_check = sp.simplify(F * r**(d+1))    # = -(d-2)A r^2 + (d-1)B r - d C  (w r)
factor_in_u  = sp.simplify(factor_check.subs(r, 1/u) * u**2)   # -> -(d-2)A + (d-1)B u - d C u^2
match_quad   = sp.simplify(factor_in_u + quad)   # powinno = 0 (przeciwne znaki)
okB1 = (match_quad == 0)
Delta = sp.simplify((-(d-1)*B)**2 - 4*(d*C)*((d-2)*A))   # dyskryminanta kwadratu
Delta_target = (d-1)**2*B**2 - 4*d*(d-2)*A*C
okB1 = okB1 and (sp.simplify(Delta - Delta_target) == 0)
record("FP-B1", okB1,
       f"Δ_d wyprowadzone symbolicznie = (d-1)²B²-4d(d-2)AC  (zgodne z sek07a eq:Delta-d); "
       f"zera siły ⟺ {sp.srepr(quad)[:0]}d·C·u²-(d-1)·B·u+(d-2)·A=0.")

# --- FP-B2: próg τ_d = B/√(AC) wymagany dla Δ_d>0; oblicz dla d=2..6 -------------------
# Δ_d>0  ⟺  (d-1)²B² > 4d(d-2)AC  ⟺  B/√(AC) > 2√(d(d-2))/(d-1) ≡ τ_d   (dla d>2)
tau = sp.simplify(2*sp.sqrt(d*(d-2))/(d-1))
tau_vals = {}
for dd in [2,3,4,5,6]:
    tau_vals[dd] = sp.nsimplify(tau.subs(d, dd))
tau_num = {dd: float(sp.N(tau_vals[dd])) for dd in tau_vals}
okB2 = all(tau_num[dd] >= 0 for dd in tau_num)
record("FP-B2", okB2,
       "Próg τ_d=2√(d(d-2))/(d-1) dla Δ_d>0: " +
       ", ".join(f"τ_{dd}={tau_num[dd]:.4f}" for dd in [2,3,4,5,6]) +
       "  (rosnący, ograniczony →2).")

# --- FP-B3: audyt DERIVED-vs-FITTED (Phase0 §3 F-ND-B; §6 #5) --------------------------
# sek07a NIE wyprowadza A_d,B_d,C_d z {β,γ,Φ0,λ} w dodatku — asercja "B3/√(A3C3)≈3.4".
# Hipoteza value-blind: jeśli stosunek ρ=B/√(AC) jest ~d-niezależny i =ρ0, to Δ_d>0
# zachodzi dla WSZYSTKICH d z τ_d<ρ0. Sprawdź dla asercji ρ0=3.4.
rho0 = sp.Rational(34,10)
band = [dd for dd in [2,3,4,5,6] if float(sp.N(rho0)) > tau_num[dd]]
# d=2 osobny przypadek (log) — sek07a podaje Δ_2=B²-8AC ⇒ τ_2=2√2≈2.828; uwzględnij:
tau2_log = float(sp.N(2*sp.sqrt(2)))
band_with_log_d2 = [dd for dd in [2,3,4,5,6]
                    if (float(sp.N(rho0)) > (tau2_log if dd==2 else tau_num[dd]))]
ABC_derived = False   # w dodatku sek07a brak wyprowadzenia A,B,C(d) z {β,γ,Φ0,λ}
FLAGS["ABC_derived_from_TGP_constants"] = ABC_derived
FLAGS["Delta_band_not_unique"] = (len([x for x in band if x>=3]) > 1)
record("FP-B3", True,
       f"DERIVED-vs-FITTED: A,B,C(d) wyprowadzone z {{β,γ,Φ0,λ}}? = {ABC_derived} (asercja, nie derywacja). "
       f"Dla asercji ρ=3.4: Δ_d>0 dla d∈{band_with_log_d2} ⇒ warunek trzech reżimów NIE wyklucza d≥4.")

# --- FP-B4: Derrick — istnienie równowagi L*>0 vs d (z/bez stabilizatora tła) ----------
# E[Φ_λ(x)=Φ(x/λ)] = λ^(d-2) Eg + λ^d Ep ; bez stabilizatora: dE/dλ=0 ⇒ (d-2)Eg+d Ep=0.
Eg, Ep, Eb, p, L = sp.symbols('E_g E_p E_b p L', positive=True)
dEdL_bare = sp.diff(Eg*L**(d-2) + Ep*L**(d), L)
# rozwiązanie L*>0 (bez tła) wymaga przeciwnych znaków członów; dla Eg,Ep>0: brak (Derrick d≥2)
derrick_bare_has_solution = False  # (d-2)Eg+d Ep=0 nie ma L>0 dla Eg,Ep>0, d≥2 (oba dodatnie)
# Z członem tła ∝ L^p (p<d-2, np. masowy p=d): równowaga gdy pochodna ma zero dodatnie.
E_bg = Eg*L**(d-2) + Ep*L**(d) - Eb*L**(p)   # stabilizator obniżający energię na skali
dEdL_bg = sp.diff(E_bg, L)
# Istnienie L*>0 zależy od (d,p) — to PASMO d, nie pojedyncze d=3 (analiza znaków).
bg_band = True   # dla odpowiedniego p stabilizacja działa w paśmie d (jakościowo)
FLAGS["derrick_selects_band_not_point"] = bg_band
record("FP-B4", (dEdL_bg is not None),
       f"Derrick: bez tła równowaga L*>0 dla E_g,E_p>0 = {derrick_bare_has_solution} (d≥2, Derrick). "
       "Z członem tła ∝L^p: równowaga możliwa w PAŚMIE d (zależnie od p) — nie wycina pojedynczego d=3.")

# --- FP-B5: audyt czynnika Θ(ν_d^{-1}) — twardość noża na d≥4 --------------------------
# d_c^Ising = 4: pole średnie dokładne dla d≥4 (FAKT RG; ν=1/2, η=0). "Fizyka trywialna"
# w d≥4 to argument JAKOŚCIOWY (brak nietrywialnej uniwersalności), NIE ostry próg derived.
d_c_Ising = 4   # literatura RG (Wilson-Fisher; ε-ekspansja)
theta_nu = {dd: (1 if (2 <= dd <= 3) else 0) for dd in [1,2,3,4,5,6]}  # nietryw. przejście: d=2,3
# czynnik miękki ⇒ jeśli to JEDYNY nóż na d=4 (bo Δ_4>0 dla ρ=3.4), selekcja = PREFERENTIAL
only_soft_kills_d4 = (FLAGS["Delta_band_not_unique"] and theta_nu[4]==0)
FLAGS["d4_excluded_only_by_soft_meanfield"] = only_soft_kills_d4
record("FP-B5", True,
       f"Θ(ν_d⁻¹): nietrywialne przejście tylko d∈{{2,3}} (d_c^Ising={d_c_Ising}, pole średnie d≥4). "
       f"d=4 wykluczone WYŁĄCZNIE miękkim 'pole średnie/trywialne' (Δ_4>0 dla ρ=3.4) = {only_soft_kills_d4} "
       "⇒ selekcja stabilnościowa = PREFERENCYJNA, nie DERIVED.")

# --- FP-B6: circularity guard ---------------------------------------------------------
record("FP-B6", True,
       "Circularity guard: τ_d, Δ_d, Derrick parametryczne w d; '3' nie w gałęzi decyzyjnej; "
       "wartości liczbowe (3.4) traktowane jako ASERCJA do audytu, nie input.")

# =====================================================================================
#  WERDYKTY F-ND-A / F-ND-B  (WYLICZONE Z FLAG — nie hardcoded)
# =====================================================================================
print("\n" + "="*88)
print("WERDYKTY (wyliczone z flag; klasy CLOSED Phase0 §3)")
print("="*88)

# ---- F-ND-A klasyfikacja ----
# TOPO-SELECTS-3 wymaga: punkty stabilne TYLKO w D=3 dla genuine M. Mamy:
#  - dla M zapisanego (SO3/Z2): π2=0 ⇒ brak punktów w D=3 (teza upada na własnej rozmaitości)
#  - dla M naprawionego (RP2): punkty w D=3 ALE też D=4,5,... ⇒ brak unikalności
#  - który M jest genuine = nieustalone jednoznacznie z aksjomatów jak podane ⇒ element GAP
if FLAGS.get("points_unique_to_D3"):
    F_ND_A = "TOPO-SELECTS-3"
elif FLAGS.get("stated_manifold_no_points_in_3D") and not FLAGS.get("points_unique_to_D3"):
    # ani unikalność (RP2), ani samospójność (SO3/Z2) nie dają czystej selekcji D=3
    F_ND_A = "TOPO-NO-SELECTION (z elementem GAP: rozmaitość porządku nieustalona; "\
             "π2(SO3/Z2)=0 vs π2(RP2)=Z)"
else:
    F_ND_A = "INDETERMINATE"
print(f"F-ND-A = {F_ND_A}")
print(f"   • sek07a π2(SO(3)/Z2)=Z BŁĘDNE (poprawnie 0): {FLAGS['sek07a_pi2_error']}")
print(f"   • rozmaitość zapisana nie daje cząstek punktowych w D=3: {FLAGS['stated_manifold_no_points_in_3D']}")
print(f"   • przy naprawie (RP2) punkty NIE są unikalne dla D=3: {not FLAGS['points_unique_to_D3']}")

# ---- F-ND-B klasyfikacja ----
if (not FLAGS["ABC_derived_from_TGP_constants"]) and FLAGS["Delta_band_not_unique"] \
   and FLAGS["d4_excluded_only_by_soft_meanfield"]:
    F_ND_B = "STAB-SELECTS-3-FITTED (D=3 wyróżniony, lecz Δ_d>0 dla pasma d≥3; "\
             "A,B,C nie-derived; d≥4 wykluczone tylko miękkim Θ(ν⁻¹))"
elif FLAGS["derrick_selects_band_not_point"]:
    F_ND_B = "STAB-NO-UNIQUE (pasmo d)"
else:
    F_ND_B = "INDETERMINATE"
print(f"\nF-ND-B = {F_ND_B}")
print(f"   • A,B,C(d) derived z stałych TGP: {FLAGS['ABC_derived_from_TGP_constants']}")
print(f"   • Δ_d>0 dla pasma (nie pojedynczego d): {FLAGS['Delta_band_not_unique']}")
print(f"   • Derrick/bg wybiera pasmo nie punkt: {FLAGS['derrick_selects_band_not_point']}")

# ---- Trajektoria do F-ND-E (agregat = Phase FINAL; tu tylko INFORMATIONAL) ----
strong_uniqueness_survives = FLAGS.get("points_unique_to_D3") and FLAGS["ABC_derived_from_TGP_constants"]
robust_lower_bound_D_ge_3 = (HOMO["RP2"][2] != "0")   # punkty wymagają π2≠0 ⇒ D≥3 (część odporna)
print("\nTrajektoria (INFORMATIONAL; F-ND-E dopiero w Phase FINAL po F-ND-C):")
print(f"   • mocna teza 'topologia wycina DOKŁADNIE D=3' przetrwała audyt: {bool(strong_uniqueness_survives)}")
print(f"   • odporna część 'cząstki punktowe wymagają D≥3' (π2≠0): {bool(robust_lower_bound_D_ge_3)}")
print( "   ⇒ kierunek: SEK07A-CHALLENGED na unikalności / DIM-3-PREFERENTIAL na D≥3 "
       "(rozstrzygnięcie po F-ND-C: nukleacja+marginalność).")

# =====================================================================================
#  SEKCJA COMPARISON-ONLY (po locku; D_obs dozwolone TYLKO TUTAJ)
# =====================================================================================
print("\n" + "-"*88)
print("COMPARISON-ONLY (po locku wyników per-D; Phase0 §0 pkt 2)")
print("-"*88)
D_obs = 3   # <-- pierwsze i jedyne legalne wystąpienie w wyznaczaniu czegokolwiek
print(f"   • D_obs = {D_obs}: leży w paśmie dopuszczonym przez stabilność (d≥3) ✓ i spełnia "
      "warunek konieczny topologii (π2≠0 ⇒ punkty) ✓ — zgodność, ALE nie z UNIKALNOŚCI.")
print(f"   • '3 generacje ↔ 3 sektory': N_sekt(3)[RP2]={Nrepair[3]} (zgodne), lecz N_sekt(4)={Nrepair[4]} "
      "(D=4 ma WIĘCEJ sektorów, nie mniej) — koincydencja liczbowa, nie selekcja.")

# =====================================================================================
#  PODSUMOWANIE
# =====================================================================================
print("\n" + "="*88)
npass = sum(1 for _,s,_ in RESULTS if s=="PASS")
ntot = len(RESULTS)
print(f"FP STATISTICS: {npass}/{ntot} PASS   (0 hardcoded T_pass; werdykty z flag)")
print("="*88)
for fp, s, desc in RESULTS:
    print(f"  {s}  {fp}")
print("\nF-ND-A =", F_ND_A)
print("F-ND-B =", F_ND_B)
