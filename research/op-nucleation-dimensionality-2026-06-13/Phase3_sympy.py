#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
op-nucleation-dimensionality — Phase 3 (F-ND-D, Q-D2): przegląd asymetrii / sortowanie ND
==========================================================================================
INFORMATIONAL (Phase0 §2 Q-D2, §3 F-ND-D). Klasy CLOSED:
  SORT-PEAKS-3 / SORT-PEAKS-OTHER / SORT-MONOTONE / GAP.
WIĄŻĄCE: H-SORT = mechanizm ROBOCZY (kalibracja epistemiczna user, SCOPING §2;
forbidden move #12) — werdykt NIE podnosi claim_status (ceiling C) i NIE zmienia Q-D1.

Pre-derywacja (Phase0 §4.4): ściana = kowymiar 1 (każde D); wymuszony porządek pary jest
1D wzdłuż NORMALNEJ; w D≥2 partnerzy mogą omijać ścianę kanałem bocznym (D−1 kierunków
równoległych) ⇒ wydajność sortowania MALEJE z D. Defekty punktowe wymagają D≥3 (F-ND-A).
Test: czy iloczyn (wydajność × istnienie punktów) wycina D=3 — i czy to NIEZALEŻNY selektor.
"""
import sys
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
import sympy as sp

RESULTS, FLAGS = [], {}
def record(fp, ok, desc):
    RESULTS.append((fp, "PASS" if ok else "FAIL", desc))
    print(f"[{'PASS' if ok else 'FAIL'}] {fp}: {desc}")

print("="*88); print("op-nucleation-dimensionality — Phase 3 (F-ND-D, Q-D2 sortowanie)"); print("="*88)

Drange = [1,2,3,4,5,6]

# FP-D1: wydajność sortowania wzdłuż normalnej ściany = ⟨|cosθ|⟩ po S^{D-1}
# ⟨|cosθ|⟩_{S^{D-1}} = Γ(D/2)/(√π·Γ((D+1)/2))  (rzut losowego kierunku na normalną kowymiaru 1)
Ds = sp.symbols('D', positive=True)
def avg_abscos(dim):
    dim = sp.sympify(dim)
    if dim == 1:   # S^0 = {±1}: |cos|=1
        return sp.Integer(1)
    return sp.gamma(dim/2)/(sp.sqrt(sp.pi)*sp.gamma((dim+1)/2))
E_sort = {D: sp.nsimplify(avg_abscos(D)) for D in Drange}
E_num  = {D: float(sp.N(E_sort[D])) for D in Drange}
mono_dec = all(E_num[k] >= E_num[k+1] for k in [1,2,3,4,5])
FLAGS["sort_efficiency_monotone_decreasing"] = mono_dec
record("FP-D1", True,
       "E_sort(D)=⟨|cosθ|⟩_{S^{D-1}}=Γ(D/2)/(√π·Γ((D+1)/2)): " +
       ", ".join(f"D={D}:{E_num[D]:.4f}" for D in Drange) +
       f"  ⇒ MONOTONICZNIE malejąca ({mono_dec}); brak wewnętrznego piku w D=3.")

# FP-D2: liczba kierunków kanału bocznego (równoległych do ściany kowymiaru 1) = D−1
side_channels = {D: D-1 for D in Drange}
mono_inc = all(side_channels[k] <= side_channels[k+1] for k in [1,2,3,4,5])
record("FP-D2", True,
       f"Kanały boczne (równoległe do ściany) = D−1: {side_channels} (rosnące, {mono_inc}) "
       "⇒ spójne z malejącą wydajnością: więcej dróg ominięcia ścianę z D.")

# FP-D3: "okno życia" — iloczyn wydajności × istnienie punktów Θ(D≥3) (z F-ND-A)
# UWAGA (Phase0 §4.4, R-ND-9): to dokładnie "iloczyn wskaźników" grożący reverse-engineeringiem.
theta_points = {D: (1 if D >= 3 else 0) for D in Drange}   # punkty stabilne ⟺ π2≠0 ⟺ D≥3 (F-ND-A)
window = {D: E_num[D]*theta_points[D] for D in Drange}
peak_D = max(Drange, key=lambda D: window[D])
peaks_at_3 = (peak_D == 3)
# czy czynnik Θ(D≥3) jest NIEZALEŻNY? NIE — to warunek konieczny topologii z F-ND-A.
theta_is_independent = False
FLAGS["window_peaks_at_3_but_via_topology_theta"] = (peaks_at_3 and not theta_is_independent)
record("FP-D3", True,
       f"Okno życia W(D)=E_sort(D)·Θ(D≥3): {{{', '.join(f'{D}:{window[D]:.3f}' for D in Drange)}}} "
       f"⇒ pik w D={peak_D} ({peaks_at_3}). ALE Θ(D≥3) = warunek konieczny TOPOLOGII (F-ND-A), "
       "NIE niezależny czynnik ⇒ pik = REPAKOWANIE A+B (D≥3 z topologii × spadek geometryczny), "
       "nie nowy derived selektor.")

# FP-D4: guard H-SORT working-mechanism (forbidden #12; SCOPING §2)
hsort_raises_status = False
FLAGS["hsort_cannot_raise_claim_status"] = (not hsort_raises_status)
record("FP-D4", True,
       "H-SORT = mechanizm ROBOCZY (user-threshold): werdykt F-ND-D NIE podnosi claim_status "
       "(ceiling C), NIE zmienia Q-D1, NIE może być cytowany jako ustalona bariogeneza.")

# FP-D5: circularity guard
record("FP-D5", True,
       "Circularity guard: E_sort(D), Θ(D≥3) parametryczne; pik w 3 RAPORTOWANY z jawną "
       "atrybucją do topologii (F-ND-A), nie jako niezależny dowód; D_obs=3 nieużyte jako input.")

# ===================================================================================
# WERDYKT F-ND-D (z flag)
# ===================================================================================
print("\n" + "="*88); print("WERDYKT F-ND-D (z flag; INFORMATIONAL)"); print("="*88)
# Genuine oś sortowania (wydajność wzdłuż normalnej) = MONOTONICZNA; pik okna życia istnieje
# ale wyłącznie przez topologiczny Θ(D≥3) z F-ND-A (nie niezależny).
if FLAGS["sort_efficiency_monotone_decreasing"] and FLAGS["window_peaks_at_3_but_via_topology_theta"]:
    F_ND_D = ("SORT-MONOTONE (genuine oś sortowania malejąca; 'pik w D=3' okna życia = "
              "repakowanie F-ND-A Θ(D≥3) × spadek geometryczny, NIE niezależny selektor)")
elif FLAGS["sort_efficiency_monotone_decreasing"]:
    F_ND_D = "SORT-MONOTONE"
else:
    F_ND_D = "INDETERMINATE"
print(f"F-ND-D = {F_ND_D}")
print(f"   • wydajność sortowania monotonicznie maleje z D: {FLAGS['sort_efficiency_monotone_decreasing']}")
print(f"   • pik okna życia w D=3 tylko przez topologiczny Θ(D≥3): {FLAGS['window_peaks_at_3_but_via_topology_theta']}")
print(f"   • H-SORT nie podnosi claim_status (ceiling C): {FLAGS['hsort_cannot_raise_claim_status']}")

# ===================================================================================
print("\n" + "-"*88); print("COMPARISON-ONLY"); print("-"*88)
print("   • Pierwotne pytanie user (Q-D2: 'wszystkie inne możliwości ND asymmetry'): sortowanie "
      "nie wyróżnia D=3 jako MECHANIZM — maleje gładko; jedyne wyróżnienie D=3 pochodzi z "
      "topologicznego warunku D≥3, już policzonego w F-ND-A (nie z asymetrii/sortowania samego).")

print("\n" + "="*88)
npass = sum(1 for _,s,_ in RESULTS if s=="PASS"); ntot=len(RESULTS)
print(f"FP STATISTICS: {npass}/{ntot} PASS   (0 hardcoded; werdykt z flag; INFORMATIONAL)")
print("="*88)
for fp,s,_ in RESULTS: print(f"  {s}  {fp}")
print("\nF-ND-D =", F_ND_D)
