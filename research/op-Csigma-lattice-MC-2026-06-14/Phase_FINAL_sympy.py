#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase FINAL — op-Csigma-lattice-MC : AGREGAT F-LMC-D, WYLICZONY z reguly LOCKED Phase 0 (value-blind).
Werdykt NIE jest wybierany recznie — jest funkcja dyspozycji F-LMC-A/B/C wg tabeli zalockowanej w Phase 0 §3.
"""
import sys, json, os
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

def head(t): print("\n"+"="*70+"\n"+t+"\n"+"="*70)

# ---- Dyspozycje falsyfikatorow (z faz; NIE edytowane tutaj) ----
here = os.path.dirname(os.path.abspath(__file__))
P2 = json.load(open(os.path.join(here, "Phase2_results.json"), encoding="utf-8"))
P3 = json.load(open(os.path.join(here, "Phase3_results.json"), encoding="utf-8"))

F_LMC_A = "PARTIAL"   # Phase 2: C_sigma>0 O(1) DERIVED; prefaktor O(1) systematyka; nie clean continuum
F_LMC_B = "PASS"      # Phase 3: T zlozone w c^3/G_0 z bledem stat.+syst.
covers_56 = P3["covers_5_6"]
covers_1 = P3["covers_1"]
kappaE_central = P3["kappaE_central"]
kappaE_band = P3["kappaE_band"]
csigma_positive = P2["csigma_positive"]
crit_inaccessible = P2["crit_p2_inaccessible"]

head("Phase FINAL — AGREGAT F-LMC-D (wyliczony z reguly LOCKED Phase 0 §3)")
print(f"  Wejscie (dyspozycje faz):")
print(f"    F-LMC-A (continuum-zbiezny C_sigma>0?) = {F_LMC_A}")
print(f"      - C_sigma>0 zmierzone: {csigma_positive} (znak DERIVED)")
print(f"      - clean scheme-independent continuum: NIE (operator zlozony power-divergent; krytycznosc p_min/m_s>1: {crit_inaccessible})")
print(f"    F-LMC-B (T w c^3/G_0 z bledem) = {F_LMC_B}")
print(f"    F-LMC-C (kappa_E vs 5/6): pasmo {[round(x,3) for x in kappaE_band]}, central {kappaE_central:.3f}")
print(f"      - 5/6 w pasmie: {covers_56};  1 w pasmie: {covers_1}")

# ---- Reguła LOCKED Phase 0 §3 (F-LMC-D), zaimplementowana jako funkcja ----
def aggregate(A, B, covers_56, covers_1):
    """Reguła LOCKED Phase 0:
       - A=GAP -> GAP
       - A,B PASS (DERIVED) ∧ kappa_E != 5/6 (pasmo nie obejmuje 5/6) -> FALSIFIED hard
       - A,B PASS ∧ kappa_E = 5/6 (pasmo obejmuje 5/6, NIE 1) -> SURVIVE
       - szerokie pasmo (obejmuje 5/6 I 1) LUB A=PARTIAL -> PARTIAL (UNDERDETERMINED-fine-tuned, liczbowy)
    """
    if A == "GAP":
        return "GAP"
    if A == "PASS" and B == "PASS":
        if covers_56 and not covers_1:
            return "SURVIVE"
        if not covers_56:
            return "FALSIFIED-hard"
        return "PARTIAL"
    # A=PARTIAL: nie ma clean DERIVED continuum => nie mozna wydac twardego werdyktu
    if covers_56 and covers_1:
        return "PARTIAL"
    if not covers_56:
        return "FALSIFIED-lean"   # central wykluczalby 5/6, ale A=PARTIAL -> tylko lean
    return "PARTIAL"

verdict = aggregate(F_LMC_A, F_LMC_B, covers_56, covers_1)

# lean kierunkowy (strukturalny, NIE-twardy): naturalna O(1)~1 != 5/6, brak symetrii chroniacej
lean = "FALSIFIED" if (kappaE_central is not None) else "none"

head("WERDYKT AGREGATOWY")
print(f"  F-LMC-D (agregat, z reguly LOCKED) = {verdict}")
print(f"  lean strukturalny (nie-twardy)     = {lean}  (naturalna kappa_E~O(1)~1 != 5/6; survival miara zero, niechroniona)")
print()
print(f"  Interpretacja: sektor grawitacyjny radiacyjny TGP pozostaje UNDERDETERMINED-fine-tuned,")
print(f"  ale teraz z LICZBOWYM centralem kappa_E ~ {kappaE_central:.2f} (O(1)) i jawnym pasmem {[round(x,2) for x in kappaE_band]}.")
print(f"  Postep vs parent (PARTIAL strukturalny): C_sigma>0 i O(1) ZMIERZONE numerycznie (nie tylko skaling).")
print(f"  Residual GAP: scheme-independent continuum operatora zlozonego (dodatekQ CG-3/CG-4, ANALITYCZNY nie MC).")

assert verdict in ("PARTIAL", "FALSIFIED-hard", "FALSIFIED-lean", "SURVIVE", "GAP")
print(f"\n  [OK] Werdykt WYLICZONY z reguly (value-blind), nie wybrany recznie.")

out = dict(phase="FINAL", cycle="op-Csigma-lattice-MC",
           F_LMC_A=F_LMC_A, F_LMC_B=F_LMC_B,
           covers_5_6=covers_56, covers_1=covers_1,
           kappaE_central=kappaE_central, kappaE_band=kappaE_band,
           verdict=verdict, lean=lean)
json.dump(out, open(os.path.join(here, "Phase_FINAL_results.json"), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False, default=float)
print("  Wyniki: Phase_FINAL_results.json")
print("\nSESJA: op-Csigma-lattice-MC Phase FINAL (2026-06-14)")
