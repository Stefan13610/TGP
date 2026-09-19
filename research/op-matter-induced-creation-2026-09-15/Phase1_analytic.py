#!/usr/bin/env python3
# -*- coding: ascii -*-
"""op-matter-induced-creation -- Phase 1: analityka pre-rejestrowana
(LOCK sec. 2; MD sec. 4). ZERO ewolucji numerycznej.

P1-I1: gate cytatow form (prog 1e-12, psi in {0.9,1,1.1}) + tozsamosci
       symboliczne (simplify = 0) + gate rampy lam(t).
P1-I2: PRE-REJESTRACJA ILOSCIOWA -- 0D (rho_hat=1):
       f(psi) = U(psi) + lam*psi^2/(4-3psi); psi_min(lam) (f'=0, f''>0,
       galaz ciagla od psi_min(0)=1); lam_fold (f'=f''=0).
       TABELA psi_min dla lam z listy LOCKa {0.06 ... 0.18}
       (+ kotwice deskryptywne 0.01/0.05/0.20). ZAPISANA PRZED Phase 3.
P1-I3: cytat kontekstu (psi=1 jedyny znany stan trwaly bez zrodla).
"""
import os
import sys

import numpy as np
import sympy as sp
from scipy.optimize import brentq, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import engine_core as ec  # noqa: E402

OUTP = os.path.join(HERE, "Phase1_output.txt")
OUT = []


def w(s=""):
    OUT.append(str(s))
    with open(OUTP, "w") as fh:
        fh.write("\n".join(OUT) + "\n")


BAR = "=" * 78
SEP = "-" * 78
w(BAR)
w("PHASE 1 -- analityka pre-rejestrowana (LOCK sec. 2; MD sec. 4)")
w("REJESTR [INPUT]: formy CYTAT (MD sec.1-2); rho_hat=exp(-r^2/18);")
w("  rampa lam(t)=lam*S((600+100-t)/100) [LOCK]; prog gate'u 1e-12;")
w("  punkty kontrolne psi in {0.9, 1, 1.1} [LOCK sec.2];")
w("  lista lam LOCKa: {0.06,0.08,0.10,0.12,0.14,0.16,0.18}")
w(BAR)

# ======================================================= P1-I1
w("P1-I1: gate cytatow form (sympy vs float silnika) + tozsamosci")
psi = sp.symbols("psi", positive=True)
lr = sp.symbols("lr", positive=True)   # lr = lam*rho_hat

M_s = psi**6/(4 - 3*psi)**2
Mp_s = 12*psi**5*(2 - psi)/(4 - 3*psi)**3
K_s = psi**4
Kp_s = 4*psi**3
U_s = psi**4/4 - psi**3/3
Up_s = psi**2*(psi - 1)
dU_s = (psi - 1)**2*(3*psi**2 + 2*psi + 1)/12
Umat_s = lr*psi**2/(4 - 3*psi)
Umatp_s = lr*psi*(8 - 3*psi)/(4 - 3*psi)**2
dUmat_s = lr*(psi - 1)*(psi + 4)/(4 - 3*psi)

ident = [
    ("M' - dM/dpsi", Mp_s - sp.diff(M_s, psi)),
    ("K' - dK/dpsi", Kp_s - sp.diff(K_s, psi)),
    ("U' - dU/dpsi", Up_s - sp.diff(U_s, psi)),
    ("[U-U(1)] - (psi-1)^2(3psi^2+2psi+1)/12",
     (U_s - U_s.subs(psi, 1)) - dU_s),
    ("dU_mat/dpsi - lr*psi(8-3psi)/(4-3psi)^2",
     sp.diff(Umat_s, psi) - Umatp_s),
    ("[U_mat-U_mat(1)] - lr(psi-1)(psi+4)/(4-3psi)",
     (Umat_s - Umat_s.subs(psi, 1)) - dUmat_s),
]
id_ok = True
for name, expr in ident:
    val = sp.simplify(expr)
    ok = (val == 0)
    id_ok = id_ok and ok
    w("  simplify(%s) = %s  %s" % (name, val, "PASS" if ok else "FAIL"))

pairs = [
    ("M", M_s, lambda p, l: ec.Mfun(p)),
    ("M'", Mp_s, lambda p, l: ec.Mpfun(p)),
    ("K", K_s, lambda p, l: ec.Kfun(p)),
    ("K'", Kp_s, lambda p, l: ec.Kpfun(p)),
    ("U", U_s, lambda p, l: ec.Ufun(p)),
    ("U'", Up_s, lambda p, l: ec.Upfun(p)),
    ("U-U(1)", dU_s, lambda p, l: ec.dUfun(p)),
    ("U_mat", Umat_s, None),
    ("dU_mat/dpsi", Umatp_s, lambda p, l: ec.Umat_pfun(p, l)),
    ("U_mat-U_mat(1)", dUmat_s, lambda p, l: ec.dUmatfun(p, l)),
]
w("  kontrola numeryczna (lr = lam*rho_hat = 1): |d| <= 1e-12*max(1,|v|)")
num_ok = True
for p0 in (0.9, 1.0, 1.1):
    for nm, ex, fn in pairs:
        if fn is None:
            continue
        vs = float(ex.subs({psi: sp.Float(p0, 30), lr: sp.Integer(1)}))
        vf = float(fn(np.float64(p0), np.float64(1.0)))
        d = abs(vs - vf)
        ok = d <= 1e-12*max(1.0, abs(vs))
        num_ok = num_ok and ok
        w("    psi=%.1f %-16s sympy=%+.16e float=%+.16e |d|=%.2e %s"
          % (p0, nm, vs, vf, d, "PASS" if ok else "FAIL"))

w("  gate rampy lam(t) [LOCK sec.3]: lam=1, t_off=600, Delta=100")
ramp_ok = True
for t, exp_v in ((0.0, 1.0), (600.0, 1.0), (650.0, 0.5), (700.0, 0.0),
                 (1700.0, 0.0)):
    v = ec.lam_of_t(1.0, t)
    ok = abs(v - exp_v) <= 1e-15
    ramp_ok = ramp_ok and ok
    w("    t=%7.1f lam(t)/lam = %.17f (oczek. %.1f) %s"
      % (t, v, exp_v, "PASS" if ok else "FAIL"))
x = sp.symbols("x")
S = 6*x**5 - 15*x**4 + 10*x**3
d0 = sp.diff(S, x).subs(x, 0)
d1 = sp.diff(S, x).subs(x, 1)
okd = (d0 == 0 and d1 == 0)
ramp_ok = ramp_ok and okd
w("    S'(0)=%s S'(1)=%s (gladkosc C1 na obu koncach) %s"
  % (d0, d1, "PASS" if okd else "FAIL"))
w("    lam(t)=0 DOKLADNIE dla t>=700: %s"
  % ("True" if ec.lam_of_t(3.14159, 700.0) == 0.0 else "False"))

p1i1 = id_ok and num_ok and ramp_ok
w("  P1-I1: %s" % ("PASS" if p1i1 else "FAIL"))
w(SEP)

# ======================================================= P1-I2
w("P1-I2: PRE-REJESTRACJA ILOSCIOWA -- statyczna odpowiedz nieliniowa 0D")
w("  f(psi) = U(psi) + lam*psi^2/(4-3psi)  (rho_hat=1, rdzen zrodla)")
lam_s = sp.symbols("lam", positive=True)
f_s = U_s + lam_s*psi**2/(4 - 3*psi)
fp_s = sp.simplify(sp.diff(f_s, psi))
fpp_s = sp.simplify(sp.diff(f_s, psi, 2))
w("  f'(psi)  = %s" % sp.simplify(fp_s))
w("  f''(psi) = %s" % sp.simplify(fpp_s))
# galaz lam(psi) z f'=0
lam_branch = sp.simplify(sp.solve(sp.Eq(fp_s, 0), lam_s)[0])
w("  f'=0  <=>  lam = LAM(psi) = %s" % lam_branch)
chk = sp.simplify(lam_branch - psi*(1 - psi)*(4 - 3*psi)**2/(8 - 3*psi))
w("  gate: simplify(LAM(psi) - psi(1-psi)(4-3psi)^2/(8-3psi)) = %s  %s"
  % (chk, "PASS" if chk == 0 else "FAIL"))

# lam_fold symbolicznie: f'=f''=0  <=>  dLAM/dpsi = 0
dlam = sp.simplify(sp.diff(lam_branch, psi))
roots = sp.solve(sp.Eq(sp.numer(sp.together(dlam)), 0), psi)
cand = []
for rt in roots:
    try:
        rv = complex(sp.N(rt))
    except TypeError:
        continue
    if abs(rv.imag) < 1e-20 and 0.0 < rv.real < 1.0:
        cand.append(sp.re(sp.N(rt, 30)))
w("  lam_fold (saddle-node f'=f''=0): pierwiastki d LAM/d psi = 0 w (0,1):")
for c in cand:
    w("    psi_fold = %s" % sp.N(c, 18))


def LAM(p):
    return p*(1.0 - p)*(4.0 - 3.0*p)**2/(8.0 - 3.0*p)


res = minimize_scalar(lambda p: -LAM(p), bounds=(1e-9, 1.0 - 1e-12),
                      method="bounded",
                      options={"xatol": 1e-14})
psi_fold_num = float(res.x)
lam_fold_num = float(LAM(psi_fold_num))
if cand:
    psi_fold_sym = float(sp.N(cand[0], 30))
    lam_fold_sym = float(sp.N(lam_branch.subs(psi, cand[0]), 30))
else:
    psi_fold_sym, lam_fold_sym = float("nan"), float("nan")
w("  lam_fold (sympy)   : psi_fold = %.15f  lam_fold = %.15f"
  % (psi_fold_sym, lam_fold_sym))
w("  lam_fold (numeryka): psi_fold = %.15f  lam_fold = %.15f"
  % (psi_fold_num, lam_fold_num))
w("  zgodnosc sympy/numeryka: |d lam_fold| = %.3e"
  % abs(lam_fold_sym - lam_fold_num))
w("  => dla lam > lam_fold lokalne minimum 0D NIE ISTNIEJE")
w("     (przewidywany MECHANIZM lam_crit; konfrontacja DESKRYPTYWNA)")


def psi_min_of(lam):
    """galaz ciagla od psi_min(0)=1: pierwiastek LAM(psi)=lam
    w (psi_fold, 1)."""
    if lam <= 0:
        return 1.0
    if lam > lam_fold_num:
        return float("nan")
    return brentq(lambda p: LAM(p) - lam, psi_fold_num, 1.0 - 1e-15,
                  xtol=1e-15, rtol=8.9e-16, maxiter=200)


def fpp(p, lam):
    return float(sp.N(fpp_s.subs({psi: sp.Float(p, 30),
                                  lam_s: sp.Float(lam, 30)}), 25))


LIST_LOCK = [0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18]
ANCHORS = [0.01, 0.05, 0.20]
w("")
w("  TABELA psi_min(lam) -- LISTA LOCKa (PRE-REJESTRACJA, PRZED Phase 3)")
w("    lam     psi_min      f''(psi_min)   psi_min<5/6?  d_psi_min=psi_min-1")
tab = {}
for lam in LIST_LOCK:
    pm = psi_min_of(lam)
    tab[lam] = pm
    w("    %.3f   %.9f  %+.6e  %-5s        %+.9f"
      % (lam, pm, fpp(pm, lam), "TAK" if pm < 5.0/6.0 else "nie", pm - 1.0))
w("  [INPUT-MD deskryptywne, dopisane PRZED Phase 3] kotwice spoza listy:")
for lam in ANCHORS:
    pm = psi_min_of(lam)
    tab[lam] = pm
    w("    %.3f   %.9f  %+.6e  %-5s        %+.9f"
      % (lam, pm, fpp(pm, lam), "TAK" if pm < 5.0/6.0 else "nie", pm - 1.0))
w("  prog detektora M911: 5/6 = %.9f" % (5.0/6.0))
lam_sub = None
for lam in sorted(tab):
    if tab[lam] == tab[lam] and tab[lam] < 5.0/6.0:
        lam_sub = lam if lam_sub is None else min(lam_sub, lam)
w("  najnizsze lam w tabeli z psi_min < 5/6 : %s" % lam_sub)

# --- dodatek deskryptywny [INPUT-MD], PRZED Phase 3 ---------------
w("")
w("  [INPUT-MD deskryptywne, dopisane PRZED Phase 3; NIE bramkuje")
w("   zadnego werdyktu] 0D zachowawczy punkt zwrotu ze startu psi=1,")
w("   psidot=0: DV(psi) = [U(psi)-U(1)] + lam*[U_mat(psi)-U_mat(1)] = 0")
w("   => lam_turn(psi) = -[U-U(1)]/[U_mat-U_mat(1)] ; psi_turn(lam):")


def lam_turn(p):
    dU = (p - 1.0)**2*(3*p*p + 2*p + 1.0)/12.0
    dm = (p - 1.0)*(p + 4.0)/(4.0 - 3.0*p)
    return -dU/dm


def psi_turn_of(lam):
    lo, hi = 1e-9, 1.0 - 1e-9
    if lam >= lam_turn(lo):
        return 0.0
    return brentq(lambda p: lam_turn(p) - lam, lo, hi, xtol=1e-14)


w("    lam_turn(psi->0+) = %.9f  (powyzej tej wartosci 0D-trajektoria"
  % lam_turn(1e-12))
w("     ze startu psi=1 siega PODLOGI psi=0 -- 0D-owy prog dynamiczny)")
w("    lam     psi_turn(0D, zachowawczy)")
for lam in [0.01, 0.05] + LIST_LOCK + [0.20]:
    pt = psi_turn_of(lam)
    w("    %.3f   %s" % (lam, ("%.9f" % pt) if pt > 0 else
                         "0 (siega podlogi)"))
w(SEP)

# ======================================================= P1-I3
w("P1-I3: kontekst predykcji Q-I2 (CYTAT, bez obliczen)")
w("  ../op-metric-pair-M911-2026-09-02/NEEDS.md:")
w("    \"Werdykty cyklu: Q-A-PASS (sektor samodomkniety bez podlog")
w("     i barier) + Q-B-FAIL (wszystko relaksuje do jednorodnej prozni")
w("     psi=1; zero nukleacji, zero 'pole wybiera granice').\"")
w("  Phase0_balance.md sec.0 (LOCK tego cyklu):")
w("    \"w galezi zdrowej BEZ zrodla nie ma statycznych solitonow")
w("     (Q-B-FAIL) ani oscylonow (Q-E/Q-G) ==> oczekiwany wynik Q-I2:")
w("     RETURN-TO-VACUUM lub COLLAPSE\"")
w("  Phase0_balance.md sec.2 (P1-I3, litera): \"po wygaszeniu (lam=0)")
w("    jedyne znane stany trwale to psi=1 (Q-B-FAIL, Q-A-PASS)\"")
w("  => predykcja pre-rejestrowana Q-I2 NIENARUSZALNA (LOCK sec.6).")
w("  P1-I3: PASS (cytat kontekstu)")
w(BAR)
w("PHASE 1 PODSUMOWANIE: P1-I1 %s ; P1-I2 ZAPISANE (tabela + lam_fold)"
  % ("PASS" if p1i1 else "FAIL"))
w("  lam_fold = %.12f (psi_fold = %.12f)" % (lam_fold_num, psi_fold_num))
w("  psi_min(lam) dla listy LOCKa: " +
  " ".join("%.2f:%.6f" % (l, tab[l]) for l in LIST_LOCK))
w("  P1-I3 PASS (cytat)")
w(BAR)
print("PHASE1 DONE")
