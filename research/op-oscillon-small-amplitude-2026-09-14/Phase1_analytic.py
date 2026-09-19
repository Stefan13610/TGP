#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-oscillon-small-amplitude -- Phase 1 (analityczna; MD sec. 7).
P1a': gate cytatow form sympy(bin double) vs float, psi in {0.9,1,1.1},
      prog |d| <= 1e-12*max(1,|v|); dodatkowo simplify=0 dla cytatu M'
      i tozsamosci energetycznej dU (correction note 2 poprzednika).
P1b': CYTAT omega_2 = -139/24 (poprzednik op-r3-stationary-states
      Phase1_output.txt -- BEZ ponownego wyprowadzania) + tabela
      omega(a) ~= 1 + omega_2 a^2 (MIEKKA, deskryptywna; MD sec. 6).
Output: Phase1_output.txt. Zero ewolucji numerycznej.
"""
import sympy as sp

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-oscillon-small-amplitude-2026-09-14/")

import sys
sys.path.insert(0, BASE)
import engine_core as ec

out = []
psi = sp.Symbol('psi', positive=True)

# ---- formy CYTAT (MD sec. 1) --------------------------------------
M_s = psi**6/(4 - 3*psi)**2
Mp_cyt = 12*psi**5*(2 - psi)/(4 - 3*psi)**3
K_s = psi**4
Kp_s = 4*psi**3
U_s = psi**4/4 - psi**3/3
Up_s = psi**2*(psi - 1)
U2_s = 3*psi**2 - 2*psi
dU_s = (psi - 1)**2*(3*psi**2 + 2*psi + 1)/12

out.append("="*78)
out.append("PHASE 1 -- op-oscillon-small-amplitude (MD sec. 7; zero ewolucji)")
out.append("Formy CYTAT (MD sec. 1; zrodlo: op-action-audit Phase1_output.txt")
out.append("  via op-r3-stationary-states Phase1_output.txt):")
out.append("  M = psi^6/(4-3psi)^2 ; M' = 12 psi^5(2-psi)/(4-3psi)^3")
out.append("  K = psi^4 ; K' = 4 psi^3 ; U = psi^4/4 - psi^3/3 ;")
out.append("  U' = psi^2(psi-1) ; U'' = 3psi^2-2psi ; K_geo=gamma=c0=1 [LOCK]")
out.append("")
out.append("-"*78)
out.append("P1a': tozsamosci symboliczne (simplify = 0 wymagane):")
chk1 = sp.simplify(sp.diff(M_s, psi) - Mp_cyt)
out.append("  simplify(dM/dpsi - M'_cytat) = %s" % chk1)
chk2 = sp.simplify(dU_s - (U_s - U_s.subs(psi, 1)))
out.append("  simplify(dU_tozsamosc - (U(psi)-U(1))) = %s   "
           "(correction note 2 poprzednika -- dziedziczona)" % chk2)
chk3 = sp.simplify(sp.diff(U_s, psi) - Up_s)
out.append("  simplify(dU/dpsi - U'_cytat) = %s" % chk3)
chk4 = sp.simplify(sp.diff(K_s, psi) - Kp_s)
out.append("  simplify(dK/dpsi - K'_cytat) = %s" % chk4)
sym_ok = (chk1 == 0) and (chk2 == 0) and (chk3 == 0) and (chk4 == 0)
out.append("  tozsamosci: %s" % ("PASS" if sym_ok else "FAIL"))
out.append("")

# ---- gate sympy(bin double) vs float (prog 1e-12) -----------------
out.append("P1a' gate: sympy(na binarnym double) vs float (engine_core),")
out.append("  psi in {0.9, 1, 1.1}; prog |d| <= 1e-12*max(1,|v|)")
forms = [
    ("M",  M_s,  ec.Mfun),
    ("Mp", Mp_cyt, ec.Mpfun),
    ("K",  K_s,  ec.Kfun),
    ("Kp", Kp_s, ec.Kpfun),
    ("U",  U_s,  ec.Ufun),
    ("Up", Up_s, ec.Upfun),
    ("U2", U2_s, None),
    ("dU", dU_s, ec.dUfun),
]


def U2float(p):
    return 3.0*p*p - 2.0*p


npass = 0
ntot = 0
allpass = True
for pv in (0.9, 1.0, 1.1):
    pbin = sp.Rational(pv)  # dokladna binarna reprezentacja double
    for name, s_expr, ffun in forms:
        sv = float(sp.N(s_expr.subs(psi, pbin), 30))
        fv = ffun(pv) if ffun is not None else U2float(pv)
        d = abs(sv - fv)
        tol = 1e-12*max(1.0, abs(sv))
        ok = d <= tol
        allpass = allpass and ok
        ntot += 1
        npass += int(ok)
        out.append("  psi=%.1f %-3s sympy(bin)=%+.16e float=%+.16e "
                   "|d|=%.2e %s" % (pv, name, sv, fv, d,
                                    "PASS" if ok else "FAIL"))
out.append("  P1a': %s (%d/%d; + tozsamosci %s)"
           % ("PASS" if (allpass and sym_ok) else "FAIL",
              npass, ntot, "PASS" if sym_ok else "FAIL"))
out.append("")
out.append("-"*78)

# ---- P1b': CYTAT omega_2 + tabela omega(a) ------------------------
w2 = sp.Rational(-139, 24)
out.append("P1b': omega_2 -- CYTAT (op-r3-stationary-states Phase1_output.txt,")
out.append("  bez ponownego wyprowadzania):")
out.append("  'omega_2 = -139/24' ; float: %.6f" % float(w2))
out.append("  (tam: omega_2 = M1^2/16 + M1*U3/8 - M2/8 - 5*U3^2/48 + U4/16,")
out.append("   M'(1)=12, M''(1)=156, U'''(1)=4, U''''(1)=6)")
out.append("")
out.append("Tabela przewidywanych omega(a) ~= 1 + omega_2 a^2 -- MIEKKA,")
out.append("deskryptywna (MD sec. 6; pre-rejestrowany ZNAK omega_peak<1")
out.append("i monotonia w a; BEZ progu PASS/FAIL) -- zapisana PRZED Phase 3:")
for a in (0.02, 0.05, 0.08, 0.10):
    wa = 1.0 + float(w2)*a*a
    out.append("  a=%.2f  omega(a) = %.6f" % (a, wa))
out.append("")
out.append("Uwaga deskryptywna (nie zmienia progow): dla a=0.02 mapa LP daje")
out.append("omega=0.9977 > 0.99 (prog detektora, warunek 3) -- prog pozostaje")
out.append("litera LOCKa.")
out.append("")
out.append("="*78)
out.append("PHASE 1 PODSUMOWANIE: P1a' %s; P1b' CYTAT omega_2 = -139/24"
           % ("PASS" if (allpass and sym_ok) else "FAIL"))
out.append("  (-5.791667) + tabela omega(a) zapisana PRZED Phase 3")
out.append("="*78)

with open(BASE + "Phase1_output.txt", "w") as f:
    f.write("\n".join(out) + "\n")
print("PHASE1 DONE; P1a' %s" % ("PASS" if (allpass and sym_ok) else "FAIL"))
