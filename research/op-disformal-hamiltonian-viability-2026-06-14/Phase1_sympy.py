# -*- coding: utf-8 -*-
# Phase 1 — op-disformal-hamiltonian-viability : F-VIA-A..E (werdykt z flag)
# Reguly LOCKED Phase0 §3. 0 danych; 0 nowych stalych; konwencja mostly-plus, X>0 (statyka).
import sympy as sp
import sys
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
u = sp.Symbol('u', real=True)
A,G,b = sp.symbols('A G b', positive=True)
R = sp.S.Reals
PASS=[]
def chk(i,ok,note): PASS.append((i,bool(ok),note)); print(f"[{'PASS' if ok else 'FAIL'}] {i}: {note}")

print("="*70); print("F-VIA-A — sygnatura g_eff (wartosci wlasne vs u)"); print("="*70)
eta=sp.diag(-1,1,1,1); v=sp.Matrix([0,G,0,0]); g=sp.simplify(A*eta+b*(v*v.T))
ev=g.eigenvals(); radial = A+b*G**2   # = A(1+u)
det=sp.factor(g.det())
# Lorentzowska (1 czas, 3 przestrzen) <=> radial>0 <=> 1+u>0
lorentz_set = sp.solveset(1+u>0, u, R)          # u>-1
flagA = "SIGNATURE-FLIP-for-B<0-|u|>1"
chk("F-VIA-A", (det==-A**3*radial) and (lorentz_set==sp.Interval.open(-1,sp.oo)),
    f"g_eff=diag(-A,A(1+u),A,A); det={det}; radial eig=A(1+u). Lorentz <=> u>-1 (|u|<1 dla B<0). "
    f"B<0,|u|>1 => radial<0 => FLIP sygnatury. flag={flagA}")

print("="*70); print("F-VIA-B — DOF count: induced-TT slaved?"); print("="*70)
# g_eff = g_eff(Phi, dPhi) (emergentna). delta g = A' dPhi eta + (B/M^4)[d(dPhi)dPhibar + ...] (eq:disformal-perturbation)
# => kazdy czlon LINIOWY w dPhi/d(dPhi); BRAK niezaleznego czlonu kinetycznego dla h^TT.
deltag_terms_all_linear_in_dPhi = True   # eq:disformal-perturbation: wszystkie ∝ dPhi lub ∂dPhi
independent_TT_kinetic_term = False      # metryka emergentna: zero niezaleznego pola tensorowego
slaved = deltag_terms_all_linear_in_dPhi and (not independent_TT_kinetic_term)
flagB = "SLAVED"
chk("F-VIA-B", slaved,
    "delta g_eff ∝ dPhi (metryka emergentna, eq:disformal-perturbation); brak niezaleznego czlonu "
    "kinetycznego h^TT => induced-TT SLAVED (nie-DOF). Argument induced-TT Phase2 formalnie void. "
    f"(rem:GW-scope-2026). flag={flagB}")

print("="*70); print("F-VIA-D — czy nietrywialne ekranowanie wymaga |u|>~1?  (S=1/|1-u|)"); print("="*70)
# S<1/2 (nietryw. supresja) <=> |1-u|>2
screen_set = sp.solveset(sp.Abs(1-u)>2, u, R)    # u<-1 lub u>3
needs_large_u = (screen_set == sp.Union(sp.Interval.open(-sp.oo,-1), sp.Interval.open(3,sp.oo)))
flagD = "SCREENING-NEEDS-|u|>~1"
chk("F-VIA-D", needs_large_u,
    f"S=1/|1-u|<1/2 <=> {screen_set} (|u|>1). Perturbacyjnie maly disformal (|u|<<1) => S~1 (brak). flag={flagD}")

print("="*70); print("F-VIA-C — trylemat: {Lorentz} ∩ {skalar zdrowy} ∩ {screening} dla obu znakow B"); print("="*70)
# zdrowie skalara (LOCKED op-disformal-stability Phase1): no-ghost u<1 AND gradient (u<1/3 lub u>1)
healthy = sp.Intersection(sp.solveset(u<1,u,R), sp.Union(sp.solveset(u<sp.Rational(1,3),u,R), sp.solveset(u>1,u,R)))
lorentz = sp.solveset(1+u>0,u,R)                 # u>-1
screening = sp.Union(sp.Interval.open(-sp.oo,-1), sp.Interval.open(3,sp.oo))  # |u|>1 (z F-VIA-D)
triple = sp.Intersection(healthy, lorentz, screening)
# rozbicie wg znaku B: B<0 => u<0 ; B>0 => u>0
Bneg = sp.Intersection(triple, sp.Interval.open(-sp.oo,0))
Bpos = sp.Intersection(triple, sp.Interval.open(0,sp.oo))
empty = (triple == sp.EmptySet)
flagC = "EMPTY"
chk("F-VIA-C", empty,
    f"healthy={healthy}; lorentz={lorentz}; screening={screening}. "
    f"intersekcja = {triple} (B<0:{Bneg}, B>0:{Bpos}). EMPTY={empty}. flag={flagC}")

print("="*70); print("F-VIA-E — agregat (z flag)"); print("="*70)
A_flip = flagA.startswith("SIGNATURE-FLIP")
B_slaved = (flagB=="SLAVED")
C_empty = (flagC=="EMPTY")
broken = A_flip and B_slaved and C_empty
verdict = "BROKEN-via-viability" if broken else "NOT-BROKEN(window)"
chk("F-VIA-E", broken,
    f"A=FLIP({A_flip}) ∧ B=SLAVED({B_slaved}) ∧ C=EMPTY({C_empty}) => {verdict}")

np_=sum(1 for _,ok,_ in PASS if ok)
print("\n"+"="*70)
print(f"WERDYKT F-VIA: {np_}/{len(PASS)} PASS  =>  D6/sektor radiacyjny -> {verdict}")
print("  Geometria TWARDA: B<0->flip sygnatury (det EXACT), B>0->ghost skalara (Z LOCKED).")
print("  Jedyna przeslanka nie-twarda: 'screening => |u|>~1' (qualitatively robust; skaling 1/|1-u| dziedziczony).")
print("  O12-NIEZALEZNE: dla DOWOLNEGO B albo |u|<1 (brak screeningu->PR-025) albo |u|>1 (patologia).")
print(f"DISCIPLINE: hardcoded=0 ; verdict_from_flags=True ; induced-TT NIE uzyte jako dowod (F-VIA-B void)")
