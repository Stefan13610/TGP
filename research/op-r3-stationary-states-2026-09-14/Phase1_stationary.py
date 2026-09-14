#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-r3-stationary-states (Phase 1) -- analityka stanow stacjonarnych.
LOCK: Phase0_balance.md sec. 2 Phase 1; decyzje FROZEN:
Phase_method_decisions.md (MD) sec. 1, 8.

P1a: linearyzacja stacjonarna psi = 1 + eps e^{-i om t} f(r):
     nabla^2 f = -kappa^2 f, kappa^2 = (M(1)om^2 - U''(1))/K(1)
     = om^2 - 1 (wyprowadzenie sympy); klasy om<1 / om>1;
     mapowanie na linearyzacje R3 (kappa=1 <=> om^2=2).
P1b: Lindstedt-Poincare O(a^2) dla modu jednorodnego rdzenia --
     znak omega_2 = PREDYKCJA pre-rejestrowana (MD sec. 8).
P1c (gate): M,M',K,K',U,U',U'',U''',U'''' w psi in {0.9,1,1.1}:
     sympy(na binarnym double) vs float, |d| <= 1e-12*max(1,|v|).

Formy CYTAT (MD sec. 1, dziedziczone z op-action-audit
Phase1_output.txt): M = psi^6/(4-3psi)^2; K = psi^4;
U = psi^4/4 - psi^3/3; K_geo=gamma=c0=1 [LOCK].
ZERO ewolucji w tej fazie.
"""
import sys
import sympy as sp

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-r3-stationary-states-2026-09-14/")

out = []


def em(s=""):
    print(s, flush=True)
    out.append(s)


psi, r, t, eps, om, kap = sp.symbols('psi r t epsilon omega kappa',
                                     positive=False)

# ---- formy (CYTAT, MD sec. 1) ------------------------------------
M = psi**6 / (4 - 3*psi)**2
K = psi**4
U = psi**4/4 - psi**3/3

Mp = sp.diff(M, psi)
Kp = sp.diff(K, psi)
Up = sp.diff(U, psi)

em("=" * 78)
em("PHASE 1 -- analityka stanow stacjonarnych; K_geo=gamma=c0=1 [INPUT]")
em("Formy CYTAT (MD sec.1; dziedziczone op-action-audit Phase1_output):")
em("  M(psi)   = psi^6/(4-3psi)^2 = %s" % sp.simplify(M))
em("  K(psi)   = %s" % K)
em("  U(psi)   = %s" % U)
em("  weryfikacja M' vs cytat 12psi^5(2-psi)/(4-3psi)^3: simplify = %s"
   % sp.simplify(Mp - 12*psi**5*(2 - psi)/(4 - 3*psi)**3))
em("  U'  = %s  (= psi^2(psi-1): %s)"
   % (sp.expand(Up), sp.simplify(Up - psi**2*(psi - 1))))
em("  U'' = %s ; U''(1) = %s" % (sp.expand(sp.diff(U, psi, 2)),
                                 sp.diff(U, psi, 2).subs(psi, 1)))
em("")

# ---- rownowaznosc pi-formulacji (deskryptywnie, MD sec. 2) -------
pv, pd2, pdot = sp.symbols('v a pidot')  # v=psidot, a=psiddot
RHS = sp.symbols('RHS')  # sila przestrzenna zbiorczo
# EOM: M a + 1/2 M' v^2 = RHS  =>  a = (RHS - M'v^2/2)/M
a_eom = (RHS - Mp*pv**2/2)/M
# pi = M v ; pidot = M' v^2 + M a ; forma Hamiltona:
# pidot_H = pi^2 M'/(2M^2) + RHS  (pi = M v)
pidot_from_eom = Mp*pv**2 + M*a_eom
pidot_H = (M*pv)**2 * Mp/(2*M**2) + RHS
em("pi-formulacja (MD sec.2, deskryptywnie): pidot_EOM - pidot_H")
em("  simplify = %s  (0 wymagane)"
   % sp.simplify(pidot_from_eom - pidot_H))
em("")

# ================== P1a: linearyzacja stacjonarna =================
em("-" * 78)
em("P1a: linearyzacja stacjonarna psi = 1 + eps e^{-i om t} f(r)")
f = sp.Function('f')
Psi = 1 + eps*sp.exp(-sp.I*om*t)*f(r)
Msub = M.subs(psi, Psi)
Mpsub = Mp.subs(psi, Psi)
Ksub = K.subs(psi, Psi)
Kpsub = Kp.subs(psi, Psi)
Upsub = Up.subs(psi, Psi)
lap = sp.diff(r**2 * Ksub * sp.diff(Psi, r), r)/r**2
EOM = (Msub*sp.diff(Psi, t, 2) + sp.Rational(1, 2)*Mpsub
       * sp.diff(Psi, t)**2
       - lap + sp.Rational(1, 2)*Kpsub*sp.diff(Psi, r)**2 + Upsub)
lin = sp.expand(sp.diff(EOM, eps).subs(eps, 0))
lin = sp.simplify(lin/sp.exp(-sp.I*om*t))
em("  O(eps) EOM (po podzieleniu przez e^{-i om t}):")
em("    %s = 0" % lin)
# forma docelowa: -om^2 f - f'' - 2 f'/r + f = 0
target = -om**2*f(r) - sp.diff(f(r), r, 2) - 2*sp.diff(f(r), r)/r + f(r)
em("  simplify(O(eps) - [-om^2 f - f'' - 2f'/r + f]) = %s"
   % sp.simplify(lin - target))
# => f'' + 2f'/r = (1 - om^2) f = -kappa^2 f, kappa^2 = om^2 - 1
kappa2_formula = (M.subs(psi, 1)*om**2
                  - sp.diff(U, psi, 2).subs(psi, 1))/K.subs(psi, 1)
em("  kappa^2 = (M(1)om^2 - U''(1))/K(1) = %s" % sp.expand(kappa2_formula))
em("  simplify(kappa^2 - (om^2 - 1)) = %s"
   % sp.simplify(kappa2_formula - (om**2 - 1)))
em("  KLASY: om<1 => kappa^2<0: f ~ e^{-|kappa| r}/r (ZLOKALIZOWANE);")
em("         om>1 => kappa^2>0: f ~ sin(kappa r)/r (kontinuum oscyl.)")
om2_at_k1 = sp.solve(sp.Eq(kappa2_formula, 1), om**2)
em("  mapowanie na linearyzacje R3: kappa=1 <=> om^2 = %s"
   % om2_at_k1)
em("  zgodnosc z LOCK (kappa=1 <=> om^2=2): %s"
   % ("PASS" if om2_at_k1 == [2] else "FAIL"))
p1a_pass = (sp.simplify(kappa2_formula - (om**2 - 1)) == 0
            and om2_at_k1 == [2]
            and sp.simplify(lin - target) == 0)
em("  P1a: %s" % ("PASS" if p1a_pass else "FAIL"))
em("")

# ================== P1b: Lindstedt-Poincare O(a^2) ================
em("-" * 78)
em("P1b: Lindstedt-Poincare, mod jednorodny rdzenia (MD sec. 8)")
em("  ODE: M(1+u) om^2 u_thth + 1/2 M'(1+u) om^2 u_th^2 + U'(1+u) = 0")
a, th = sp.symbols('a theta', positive=True)
U3s, U4s, M1s, M2s = sp.symbols('U3 U4 M1 M2')  # pochodne w psi=1
u1 = sp.Function('u1')
u2 = sp.Function('u2')
u3 = sp.Function('u3')
w1s, w2s = sp.symbols('omega1 omega2')
u = a*u1(th) + a**2*u2(th) + a**3*u3(th)
omega = 1 + a*w1s + a**2*w2s
# rozwiniecia Taylora wokol psi=1 (M(1)=1, U'(1)=0, U''(1)=1):
Mu = 1 + M1s*u + M2s*u**2/2
Mpu = M1s + M2s*u
Upu = u + U3s*u**2/2 + U4s*u**3/6
eq = (Mu*omega**2*sp.diff(u, th, 2)
      + sp.Rational(1, 2)*Mpu*omega**2*sp.diff(u, th)**2 + Upu)
eq = sp.expand(eq)
orders = [sp.expand(eq.coeff(a, n)) for n in (1, 2, 3)]
em("  O(a):   %s = 0" % orders[0])
# u1 = cos(theta)
sub1 = {u1(th): sp.cos(th)}
o1 = sp.simplify(orders[0].subs(sub1).doit())
em("    u1 = cos(theta): residuum = %s" % o1)
# O(a^2): L u2 = -(reszta); warunek sekularny -> w1
o2 = sp.expand(orders[1].subs(sub1).doit())
o2 = sp.simplify(sp.expand_trig(o2))
# wspolczynnik przy cos(theta) (sekularny) bez u2:
o2_nou2 = o2.subs({u2(th): 0, sp.diff(u2(th), th, 2): 0})
sec1 = sp.integrate(o2_nou2*sp.cos(th), (th, 0, 2*sp.pi))/sp.pi
em("  O(a^2) czlon sekularny (proj. cos): %s" % sp.simplify(sec1))
w1_sol = sp.solve(sp.Eq(sec1, 0), w1s)
em("    => omega_1 = %s" % w1_sol)
w1v = w1_sol[0] if w1_sol else sp.Integer(0)
# rozwiaz u2: L u2 = rhs2 (harmoniki 0 i 2theta)
rhs2 = sp.simplify(-(o2_nou2.subs(w1s, w1v)))
rhs2 = sp.expand_trig(sp.expand(rhs2))
c0 = sp.integrate(rhs2, (th, 0, 2*sp.pi))/(2*sp.pi)
c2 = sp.integrate(rhs2*sp.cos(2*th), (th, 0, 2*sp.pi))/sp.pi
em("  rhs2 = %s + %s cos(2 theta)" % (sp.simplify(c0), sp.simplify(c2)))
u2sol = c0 - c2/3*sp.cos(2*th)   # (u2''+u2 = c0 + c2 cos2th)
chk2 = sp.simplify(sp.diff(u2sol, th, 2) + u2sol - rhs2)
em("    u2 = %s ; kontrola L u2 - rhs2 = %s" % (sp.simplify(u2sol), chk2))
# O(a^3): warunek sekularny -> w2
o3 = sp.expand(orders[2].subs(sub1).doit().subs(w1s, w1v))
o3 = o3.subs({u2(th): u2sol,
              sp.diff(u2(th), th): sp.diff(u2sol, th),
              sp.diff(u2(th), th, 2): sp.diff(u2sol, th, 2)}).doit()
o3_nou3 = o3.subs({u3(th): 0, sp.diff(u3(th), th, 2): 0})
o3_nou3 = sp.expand_trig(sp.expand(o3_nou3))
sec3 = sp.integrate(o3_nou3*sp.cos(th), (th, 0, 2*sp.pi))/sp.pi
w2_sol = sp.solve(sp.Eq(sp.simplify(sec3), 0), w2s)
em("  O(a^3) sekularny: %s = 0" % sp.simplify(sec3))
em("  => omega_2 (symbolicznie, U3=U'''(1), U4=U''''(1), M1=M'(1),")
em("     M2=M''(1)):")
w2_sym = sp.simplify(w2_sol[0])
em("     omega_2 = %s" % w2_sym)
# wartosci form:
U3v = sp.diff(U, psi, 3).subs(psi, 1)
U4v = sp.diff(U, psi, 4).subs(psi, 1)
M1v = Mp.subs(psi, 1)
M2v = sp.diff(M, psi, 2).subs(psi, 1)
em("  wartosci form: U'''(1)=%s, U''''(1)=%s, M'(1)=%s, M''(1)=%s"
   % (U3v, U4v, M1v, M2v))
w2_val = sp.simplify(w2_sym.subs({U3s: U3v, U4s: U4v, M1s: M1v,
                                  M2s: M2v}))
em("  omega_2 = %s = %s" % (w2_val, sp.nsimplify(w2_val)))
em("  omega(a) = 1 + omega_2 a^2 + O(a^4); float: omega_2 = %.6f"
   % float(w2_val))
soft = bool(w2_val < 0)
em("")
em("  PREDYKCJA PRE-REJESTROWANA (P1b, LOCK sec.2; zakaz")
em("  reinterpretacji po Phase 3): omega_2 %s 0 => nieliniowosc %s"
   % ("<" if soft else ">=", "MIEKKA" if soft else "TWARDA/NEUTRALNA"))
em("  => warunek istnienia oscylonu malej amplitudy (om(a) < m=1): %s"
   % ("SPELNIONY (oscylony malej amplitudy OCZEKIWANE)" if soft
      else "NIESPELNIONY (oscylonow malej amplitudy NIE oczekujemy;"
           " pytanie Q-E pozostaje otwarte dla duzych amplitud)"))
em("")

# ================== P1c: gate sympy vs float ======================
em("-" * 78)
em("P1c (gate): sympy(bin double) vs float w psi in {0.9, 1, 1.1};")
em("  prog |d| <= 1e-12*max(1,|v|) (MD sec. 8)")


def floats(p):
    Mf = p**6/(4.0 - 3.0*p)**2
    Mpf = 12.0*p**5*(2.0 - p)/(4.0 - 3.0*p)**3
    Kf = p**4
    Kpf = 4.0*p**3
    Uf = p**4/4.0 - p**3/3.0
    Upf = p*p*(p - 1.0)
    U2f = p*(3.0*p - 2.0)
    U3f = 6.0*p - 2.0
    U4f = 6.0
    return dict(M=Mf, Mp=Mpf, K=Kf, Kp=Kpf, U=Uf, Up=Upf,
                U2=U2f, U3=U3f, U4=U4f)


syms = dict(M=M, Mp=Mp, K=K, Kp=Kp, U=U, Up=Up,
            U2=sp.diff(U, psi, 2), U3=sp.diff(U, psi, 3),
            U4=sp.diff(U, psi, 4))
allpass = True
for pval in (0.9, 1.0, 1.1):
    pbin = sp.Rational(pval)  # dokladna binarna reprezentacja double
    fv = floats(pval)
    for name in ("M", "Mp", "K", "Kp", "U", "Up", "U2", "U3", "U4"):
        sv = float(syms[name].subs(psi, pbin))
        d = abs(sv - fv[name])
        tol = 1e-12*max(1.0, abs(sv))
        ok = d <= tol
        allpass = allpass and ok
        em("  psi=%.1f %-3s sympy(bin)=%+.16e float=%+.16e |d|=%.2e %s"
           % (pval, name, sv, fv[name], d, "PASS" if ok else "FAIL"))
em("  P1c: %s" % ("PASS" if allpass else "FAIL"))
em("")
em("=" * 78)
em("PHASE 1 PODSUMOWANIE: P1a %s; P1b PREDYKCJA: omega_2 = %s (%.6f;"
   % ("PASS" if p1a_pass else "FAIL", w2_val, float(w2_val)))
em("  nieliniowosc %s); P1c %s"
   % ("MIEKKA -- oscylony malej amplitudy oczekiwane" if soft
      else "TWARDA -- oscylonow malej amplitudy nie oczekujemy",
      "PASS" if allpass else "FAIL"))
em("  kappa^2 = om^2 - 1 [WYPROWADZONE]; kappa=1 <=> om^2=2 [zgodne R3]")
em("=" * 78)

with open(BASE + "Phase1_output.txt", "w", encoding="ascii") as fh:
    fh.write("\n".join(out) + "\n")
print("zapisano:", BASE + "Phase1_output.txt")
if not (p1a_pass and allpass):
    sys.exit(1)
