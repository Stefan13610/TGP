#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 1 - op-CG-Kij-from-Hgamma : F-CGK-A (baseline) + F-CGK-C1 (relewancja RG).
================================================================================
Reguly LOCKED w Phase0_LOCK.md (value-blind). Endpoint +4/alpha=2 NIE wbudowany.

F-CGK-A: chain-rule wykladnik K w obu ramach (amplituda phi=s, gestosc Phi=s^2),
         dla substratu (kanoniczny s) i manuskryptu (K(phi)=Kgeo phi^4).
F-CGK-C1: wymiar skalowania operatora kinetycznego O_n = Phi^n (grad Phi)^2 przy 3D
          Z2 (Ising) WF fixed point; relewantny (Delta<d=3) czy irrelewantny (Delta>3)?

Kotwice literaturowe 3D Ising CFT (conformal bootstrap, Kos-Poland-Simmons-Duffin-Vichi 2016;
do CYTOWANIA, nie zakladania jako wynik TGP):
  Delta_sigma = 0.5181489   (=> eta = 2 Delta_sigma - (d-2) = 0.0362978)
  Delta_eps   = 1.412625    (operator energii eps ~ s^2 = nasze Phi=<s^2>)
"""
import sympy as sp

# ---- kotwice literaturowe (LITERATURE, nie TGP-derived) ----
Delta_sigma = sp.Rational(5181489, 10000000)   # 0.5181489
d = sp.Integer(3)
eta_lit = sp.nsimplify(2*Delta_sigma - (d-2), rational=False)
Delta_eps = sp.Float('1.412625')               # energy operator = s^2 = Phi

print("="*70)
print("Phase 1 - F-CGK-A + F-CGK-C1 (value-blind)")
print("="*70)
print(f"  [lit 3D Ising] Delta_sigma={float(Delta_sigma)}, eta={float(eta_lit):.5f}, Delta_eps={float(Delta_eps)}")

# =====================================================================
# F-CGK-A : wykladnik K w obu ramach, substrat vs manuskrypt (chain rule)
# =====================================================================
print("\n[F-CGK-A] Wykladnik e w K~(pole)^e : chain rule, obie ramy")
phi, Phi, c = sp.symbols('phi Phi c', positive=True)

def exponent(K_expr, var):
    """wykladnik e w K ~ var^e przez log-pochodna"""
    return sp.simplify(sp.diff(sp.log(K_expr), var)*var)

# --- substrat: kanoniczny w amplitudzie (phi = s) ---
K_sub_phi = c                      # const w phi (kanoniczna)
e_sub_phi = exponent(K_sub_phi, phi) if K_sub_phi.has(phi) else sp.Integer(0)
# do gestosci Phi=phi^2 : K_Phi = K_phi * (dphi/dPhi)^2, phi=sqrt(Phi)
dphi_dPhi = sp.diff(sp.sqrt(Phi), Phi)
K_sub_Phi = sp.simplify(c*dphi_dPhi**2)
e_sub_Phi = exponent(K_sub_Phi, Phi)
print(f"  substrat (kanoniczny s): e_amplituda={e_sub_phi}, e_gestosc={e_sub_Phi}  (e_gestosc=-1 => #49 anchor)")

# --- manuskrypt: K(phi)=Kgeo phi^4 (amplituda) ---
Kgeo = sp.symbols('K_geo', positive=True)
K_man_phi = Kgeo*phi**4
e_man_phi = exponent(K_man_phi, phi)
# do gestosci: K_man_Phi = K_man_phi(phi=sqrt(Phi)) * (dphi/dPhi)^2
K_man_Phi = sp.simplify(K_man_phi.subs(phi, sp.sqrt(Phi))*dphi_dPhi**2)
e_man_Phi = exponent(K_man_Phi, Phi)
print(f"  manuskrypt (K=Kgeo phi^4): e_amplituda={e_man_phi}, e_gestosc={e_man_Phi}")

# --- luki w SPOJNEJ ramie ---
gap_phi  = sp.simplify(e_man_phi - e_sub_phi)
gap_Phi  = sp.simplify(e_man_Phi - e_sub_Phi)
print(f"  LUKA spojna: amplituda Delta_e={gap_phi} (Delta_alpha={gap_phi/2}); gestosc Delta_e={gap_Phi} (Delta_alpha={gap_Phi/2})")
# anchor F-CGK-A: e_sub_Phi == -1 ?
A_anchor = sp.simplify(e_sub_Phi - (-1)) == 0
print(f"  F-CGK-A anchor (e_sub_gestosc == -1): {A_anchor}")
# DOUBT: porownanie mieszane (#49 headline -1 -> +4)
print(f"  [W-CGK-1 DOUBT] '#49 headline' Delta_e=5 = e_sub_gestosc(-1) vs e_man_AMPLITUDA(+4) = MIESZANA rama.")
print(f"                  W spojnej ramie luka = {gap_phi} (ampl.) lub {gap_Phi} (gest.). (B) ROBUST (luka != 0 w obu).")

# =====================================================================
# F-CGK-C1 : relewancja RG operatora kinetycznego O_n = Phi^n (grad Phi)^2
# =====================================================================
print("\n[F-CGK-C1] Wymiar RG O_n = Phi^n (grad Phi)^2 przy 3D WF (Phi=s^2=eps)")
n = sp.symbols('n')
# Delta[Phi]=Delta_eps ; Delta[grad]=1 ; O_n ma (n+2) pol Phi i 2 pochodne
Delta_On = (n + 2)*Delta_eps + 2
print(f"  Delta[O_n] = (n+2)*Delta_eps + 2 = {sp.simplify(Delta_On)}")
rows = []
for nv, label in [(-1, "substrat K~Phi^-1"), (0, "kanon (grad Phi)^2"), (1, "manuskrypt(gest) K~Phi^+1"), (2, "K~Phi^+2 (sek08b g^2)")]:
    Dv = float(Delta_On.subs(n, nv))
    relevant = Dv < 3
    rows.append((nv, label, Dv, relevant))
    print(f"    n={nv:+d} ({label:28s}): Delta={Dv:.3f}  {'RELEWANTNY (<3)' if relevant else 'IRRELEWANTNY (>3)'}")

any_relevant = any(r[3] for r in rows)
print(f"  Czy ktorykolwiek operator kinetyczny RELEWANTNY? {any_relevant}")
# C1 verdict
if any_relevant:
    C1 = "kandydat-na-pinowanie (sprawdz wspolczynnik w C2)"
else:
    C1 = "C-AXIOM-contribution (wszystkie irrelewantne => exponent nie pinowany przez FP)"
print(f"  F-CGK-C1 => {C1}")

# =====================================================================
# F-CGK-B (implikacja natychmiastowa z C1): czy eta moze zmostkowac luke?
# =====================================================================
print("\n[F-CGK-B] (implikacja z C1) eta potrzebne vs bound unitarnosci")
# w spojnej ramie gestosci luka = gap_Phi (=2); ale przejscie -1->+1 to ZMIANA OPERATORA,
# nie shift eta. eta_lit maly; bound unitarnosci 3D scalar: eta >= 0, male (O(0.01-0.1)).
print(f"  luka gestosc = {gap_Phi}; to wymaga DOMINACJI operatora Phi^+1 nad Phi^-1 w IR,")
print(f"  tj. RELEWANCJI (C1), NIE shiftu eta. eta_lit={float(eta_lit):.4f} (male).")
print(f"  prog LOCKED: eta>1 = 'duze'. eta_lit << 1 => B-REFUTED (eta nie mostkuje).")
B = "B-REFUTED"

# =====================================================================
print("\n" + "="*70)
print("WYNIK Phase 1 (F-CGK-A, C1; B-implikacja):")
print(f"  A: anchor e_sub_gestosc=-1 -> {A_anchor} (#49 potwierdzone jako baseline).")
print(f"     luka spojna: {gap_phi} (ampl.) / {gap_Phi} (gest.) -- NIE 5 (W-CGK-1).")
print(f"  C1: O_n dla n in -1,0,1,2 -> wszystkie Delta>3 -> IRRELEWANTNE.")
print(f"      => caly sektor kinetyczny kompozytu Phi=eps jest RG-irrelewantny (zgodne z #39).")
print(f"  B: eta_lit<<1 -> B-REFUTED.")
print(f"  ----")
print(f"  WSTEPNIE (do agregatu F-CGK-D w Phase 2 po C2):")
print(f"  (B-REFUTED) AND (C1=wszystkie irrelewantne) => silnie ku NON-DERIVABLE")
print(f"  (alpha=2 = nieredukowalny aksjomat; exponent kinetyczny nie pinowany przez FP).")
print(f"  Pozostaje C2 (czteropolowy bond: wspolczynnik fiksowany czy wolny) -> Phase 2.")
print("="*70)
