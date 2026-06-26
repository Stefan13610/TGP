# -*- coding: utf-8 -*-
# Phase 1 — op-sigma-kinetic-Csigma : strukturalne pinowanie kappa_E (sektor sigma_ab)
# Obiekt (poprawnie zidentyfikowany): kappa_E ~ C_sigma*sigma_0^2 (strumien energii sigma_ab),
# niezalezne od xi_eff (amplituda, dopasowana do GR). C_sigma = PROBLEM OTWARTY (rem:sigma-params).
import sympy as sp
import sys
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

xi,Cs,s0 = sp.symbols('xi C_sigma sigma_0', positive=True)
# (1) Amplituda _|_ strumien (reuse op-disformal-radiation-resolution det J; tu w C_sigma,sigma_0)
Oamp  = xi/(Cs*s0)      # ~ h^TT (thm:amplitude-matching pinuje to do GR: xi=4pi G0 sigma0 Phi0)
Oflux = Cs*s0**2        # ~ strumien Isaacsona sigma_ab  (= kappa_E * GR)
J = sp.Matrix([[sp.diff(Oamp,Cs), sp.diff(Oamp,s0)],
               [sp.diff(Oflux,Cs), sp.diff(Oflux,s0)]])
detJ = sp.simplify(J.det())
print("(1) det J(amp,flux | C_sigma,sigma_0) =", detJ, " != 0  => amplituda _|_ strumien (kappa_E swobodne)")

# (2) Bilans energii: P_b = kappa_E*P_GR + (1/6)*P_GR (sigma_ab tensor + nieunikniony skalar konforemny)
kE = sp.Symbol('kappa_E', positive=True)
survive = sp.solve(sp.Eq(kE + sp.Rational(1,6), 1), kE)[0]
print(f"(2) Warunek PRZEZYCIA (P_b=P_GR): kappa_E = {survive} (=5/6). Pojedynczy punkt => MIARA ZERO.")

# (3) Wartosci charakterystyczne
for k,lab in [(sp.Integer(1),"naturalna (sigma GR-matched w strumieniu)"),
              (sp.Rational(5,6),"jedyna przezywajaca"),
              (sp.Rational(1,6),"sigma martwy radiacyjnie")]:
    tot = sp.nsimplify(k + sp.Rational(1,6))
    print(f"    kappa_E={k} ({lab}): total P_b/P_GR = {tot}",
          "=> GR (PASS)" if tot==1 else "=> != GR (deviation)")
print("    Uwaga: kappa_E=1 => 7/6 = DOKLADNIE galaz B PR-025 (2646 sigma FALSIFIED).")

# (4) Status C_sigma
print("(4) C_sigma = stala kinetyczna sigma_ab; kappa_E = C_sigma*sigma_0^2 / (wartosc GR-lock).")
print("    rem:sigma-params: 'wyznaczalny w zasadzie z H_Gamma, ale obecnie NIEZOBLICZONY; problem otwarty'.")
print("    Bez wyprowadzenia C_sigma: kappa_E nie jest przypiete => sektor UNDERDETERMINED.")
print("    ALE survival = miara zero (kappa_E=5/6); generyczne C_sigma => P_b != P_GR => FALSIFIED.")

print("\nWERDYKT STRUKTURALNY: UNDERDETERMINED (kappa_E swobodne), ale SURVIVAL = MIARA ZERO (kappa_E=5/6).")
print("Brama decydujaca: wyprowadzenie C_sigma z H_Gamma (coarse-graining sigma_ab=<s s>^TF). Problem otwarty.")
print("Gdyby C_sigma policzone: kappa_E != 5/6 (generycznie) => FALSIFIED; kappa_E=5/6 (spisek) => SURVIVE.")
