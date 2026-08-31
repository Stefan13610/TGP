#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-symplectic-Jspectrum (Phase 1) -- formalizacja analityczna (sympy,
ZERO numeryki siatkowej). LOCK: Phase0_balance.md sec. 3 Phase 1;
decyzje FROZEN: Phase_method_decisions.md sec. 1.

P1a: inwariantnosc widmowa zanurzen II rzedu: kazde zanurzenie
     F(g)g_tt + F'/2 g_t^2 = F g_xx + F'/2 g_x^2 - V'(g) odtwarzajace
     statyke zrodla S = beta g^2 - gamma g^3 (przy beta=gamma: prozni g=1)
     ma dyspersje omega^2 = k^2 + V''(1)/F(1); dla V' = F*S (dowolne F!)
     oraz dla obu nazwanych zanurzen (energia: F=K=g^4, V=U;
     zrodlo #63: F=f=1+4ln g, V'=W'=S) wychodzi omega^2 = k^2 - gamma.

P1b: gradient flow g_t = -dE/dg: operator linearyzacji == -(druga
     wariacja E) == -L_plus -- stabilnosc <=> Hessian >= 0.

P1c: derywacja L_plus i L_minus z E[u] = int[K(|u|)/2 |u_x|^2 + U(|u|)]
     (i u_t = dE/du*, u = g0 + a + ib); tozsamosc L_minus g0 = 0 on-shell
     EXACT (FAIL => STOP calego cyklu); rownowaznosc Q_plus z forma
     operatora Phase 3 cyklu bloch; czynnik 1/2 dynamiki (Wirtinger);
     dyspersja JL wokol prozni u=1: lambda^2 = -k^2(k^2-gamma)/4 (do C3).
"""
import sys
import sympy as sp
from sympy.calculus.euler import euler_equations

x, t, k, w, eps, lam = sp.symbols('x t k omega epsilon lambda', real=True)
beta, gamma = sp.symbols('beta gamma', positive=True)

FAILS = []


def check(name, expr_zero, note=""):
    ok = sp.simplify(sp.expand(expr_zero)) == 0
    print("  [%s] %s %s" % ("PASS" if ok else "FAIL", name,
                            ("-- " + note) if note else ""))
    if not ok:
        FAILS.append(name)
    return ok


print("=" * 78)
print("Phase 1 -- analityka sympy (LOCK sec. 3 Phase 1). beta=gamma=1")
print("w cyklu; tu symbolicznie (proznia g=1 wymaga beta=gamma -- tak")
print("substytuujemy). [INPUT-ONTO] zlozenie u = g_d + a + ib (modul")
print("pola = g; LOCK sec. 2 Rejestr WEJSC).")
print("=" * 78)

# =========================================================== P1a
print("\n[P1a] inwariantnosc widmowa zanurzen II rzedu (energia vs zrodlo)")
g = sp.Function('g')(x, t)
F = sp.Function('F')
V = sp.Function('V')


def pde_embedding(Ffun, Vprime_expr):
    """F(g) g_tt + F'/2 g_t^2 - F g_xx - F'/2 g_x^2 + V'(g) = 0."""
    Fg = Ffun(g)
    Fp = sp.diff(Ffun(sp.Symbol('_z')), sp.Symbol('_z')).subs(
        sp.Symbol('_z'), g)
    return (Fg * sp.diff(g, t, 2) + Fp / 2 * sp.diff(g, t) ** 2
            - Fg * sp.diff(g, x, 2) - Fp / 2 * sp.diff(g, x) ** 2
            + Vprime_expr)


def dispersion(pde):
    """Linearyzacja wokol g=1: g = 1 + eps e^{i(kx - w t)} -> omega^2."""
    pert = sp.exp(sp.I * (k * x - w * t))
    lin = pde.subs(g, 1 + eps * pert).doit()
    c1 = sp.expand(lin.series(eps, 0, 2).removeO()).coeff(eps, 1)
    c1 = sp.simplify(c1 / pert)
    sols = sp.solve(sp.Eq(c1, 0), w ** 2)
    assert len(sols) == 1
    return sp.simplify(sols[0].subs(beta, gamma))


zz = sp.Symbol('_z')
S_of = lambda q: beta * q ** 2 - gamma * q ** 3  # zrodlo M-C

# (i) OGOLNE F, V' = F*S (statyka: F g'' + F'/2 g'^2 = F S -- forma M-C
#     z dowolna waga F; obejmuje energie F=K):
w2_gen = dispersion(pde_embedding(F, F(g) * S_of(g)))
print("  ogolne zanurzenie V'=F*S, DOWOLNE F: omega^2 =", w2_gen)
check("P1a-ogolne: omega^2 == k^2 - gamma (niezaleznie od F)",
      w2_gen - (k ** 2 - gamma))

# (ii) zanurzenie ENERGETYCZNE: F = K = g^4, V = U = beta g^7/7 - gamma g^8/8
Kfun = lambda q: q ** 4
Ufun = lambda q: beta * q ** 7 / 7 - gamma * q ** 8 / 8
Uprime = sp.diff(Ufun(zz), zz)
w2_en = dispersion(pde_embedding(Kfun, Uprime.subs(zz, g)))
print("  zanurzenie energetyczne (K=g^4, U'=K*S): omega^2 =", w2_en)
check("P1a-energia: omega^2 == k^2 - gamma", w2_en - (k ** 2 - gamma))
check("P1a-energia: U' == K*S (spojnosc akcji kanonicznej)",
      Uprime - Kfun(zz) * S_of(zz))

# (iii) zanurzenie ZRODLOWE #63: F = f = 1 + 4 ln g, V' = W' = S (waga 1)
ffun = lambda q: 1 + 4 * sp.log(q)
w2_src = dispersion(pde_embedding(ffun, S_of(g)))
print("  zanurzenie zrodlowe #63 (f=1+4ln g, W'=S): omega^2 =", w2_src)
check("P1a-zrodlo: omega^2 == k^2 - gamma", w2_src - (k ** 2 - gamma))
print("  WNIOSEK P1a: obie konwencje (W_source/V_energy) i KAZDA waga F")
print("  daja te sama dyspersje prozni omega^2 = k^2 - gamma -- relabeling")
print("  widmowo inwariantny w calej klasie II rzedu.")

# =========================================================== przygotowanie
# pola statyczne 1D do P1b/P1c
# [correction note 2026-08-31] positive=True: bez tego sympy nie redukuje
# (g0^2)^(7/2) -> g0^7 (U ma nieparzyste potegi rho); fizycznie g_d > 0.
g0 = sp.Function('g_0', positive=True)(x)
A = sp.Function('a')(x)
B = sp.Function('b')(x)
g0p, g0pp = sp.diff(g0, x), sp.diff(g0, x, 2)
Kg, Kp, Kpp = [sp.diff(Kfun(zz), zz, n).subs(zz, g0) for n in (0, 1, 2)]
Ug, Up, Upp = [sp.diff(Ufun(zz), zz, n).subs(zz, g0) for n in (0, 1, 2)]

# statyczne EL (dE/dg = 0):  EL0 = -(K g0')' + K'/2 g0'^2 + U'
EL0 = sp.expand(-sp.diff(Kg * g0p, x) + Kp / 2 * g0p ** 2 + Up)
g0pp_onshell = sp.solve(sp.Eq(EL0, 0), g0pp)[0]

# operatory docelowe (Phase_method_decisions sec. 1):
Qplus = Upp - Kp * g0pp - sp.Rational(1, 2) * Kpp * g0p ** 2
Qminus = Kp * g0p ** 2 / (2 * g0) + Up / g0
Lplus_A = sp.expand(-sp.diff(Kg * sp.diff(A, x), x) + Qplus * A)
Lminus_B = sp.expand(-sp.diff(Kg * sp.diff(B, x), x) + Qminus * B)

# =========================================================== P1b
print("\n[P1b] gradient flow g_t = -dE/dg <=> Hessian (druga wariacja)")
# [correction note 2026-08-31] konwencja sympy: euler_equations lhs =
# dL/df - d/dx dL/df' = +dE/df dla E = int L dx (test: L=f'^2/2 -> -f'');
# pierwotny znak '-' byl bledem odczytu konwencji. Pole positive=True (j.w.).
gd = sp.Function('g', positive=True)(x)
dens_real = Kfun(gd) / 2 * sp.diff(gd, x) ** 2 + Ufun(gd)
dEdg = euler_equations(dens_real, [gd], [x])[0].lhs  # dE/dg
lin_flow = sp.expand(dEdg.subs(gd, g0 + eps * A).doit()
                     .series(eps, 0, 2).removeO()).coeff(eps, 1)
check("P1b: linearyzacja dE/dg wokol g0 == L_plus (Hessian)",
      lin_flow - Lplus_A,
      "a_t = -L_plus a => stabilnosc gradient flow <=> L_plus >= 0")
check("P1b: czlon rzedu 0 == statyczne EL (tlo stacjonarne w flow)",
      sp.expand(dEdg.subs(gd, g0 + eps * A).doit()
                .series(eps, 0, 2).removeO()).coeff(eps, 0) - EL0)
print("  Uwaga (Sylvester): problem II rzedu ma forme wazona L_plus phi =")
print("  omega^2 K phi z K>0 -- znaki widma pokrywaja sie ze znakami L_plus;")
print("  gradient flow i dynamika II rzedu rzadza sie TYM SAMYM Hessianem.")

# =========================================================== P1c
print("\n[P1c] derywacja L_plus / L_minus z E[u] (pole zespolone)")
p = sp.Function('p')(x)
q = sp.Function('q')(x)
rho = sp.sqrt(p ** 2 + q ** 2)
dens = (Kfun(rho) / 2 * (sp.diff(p, x) ** 2 + sp.diff(q, x) ** 2)
        + Ufun(rho))
eqs = euler_equations(dens, [p, q], [x])
# dE/dp = dL/dp - d/dx dL/dp' = lhs euler_equations (E = int dens dx);
# [correction note 2026-08-31] usunieta martwa podwojna przypiska ze
# zbednym znakiem '-'.
E_p = eqs[0].lhs
E_q = eqs[1].lhs

sub_ab = {p: g0 + eps * A, q: eps * B}
Ep_ser = sp.expand(E_p.subs(sub_ab).doit().series(eps, 0, 2).removeO())
Eq_ser = sp.expand(E_q.subs(sub_ab).doit().series(eps, 0, 2).removeO())

check("P1c: E_p|_(rzad 0) == statyczne EL (tlo = punkt stacjonarny)",
      Ep_ser.coeff(eps, 0) - EL0)
check("P1c: E_q|_(rzad 0) == 0 (brak sily w kierunku fazowym)",
      Eq_ser.coeff(eps, 0))
check("P1c: E_p|_(rzad 1) == L_plus a  (amplitudowy Hessian)",
      Ep_ser.coeff(eps, 1) - Lplus_A)
check("P1c: E_q|_(rzad 1) == L_minus b (operator fazowy)",
      Eq_ser.coeff(eps, 1) - Lminus_B)

# --- czynnik 1/2 dynamiki (Wirtinger): dE/du* = (E_p + i E_q)/2 ---
u_, v_ = sp.Function('u')(x), sp.Function('v')(x)  # v = u* NIEZALEZNE
rho_uv = sp.sqrt(u_ * v_)
dens_uv = (Kfun(rho_uv) / 2 * sp.diff(u_, x) * sp.diff(v_, x)
           + Ufun(rho_uv))
EL_v = euler_equations(dens_uv, [u_, v_], [x])[1].lhs  # dE/du*
wirt = sp.expand(EL_v.subs({u_: p + sp.I * q, v_: p - sp.I * q}).doit()
                 - (E_p + sp.I * E_q) / 2)
# uzgodnienie sqrt((p+iq)(p-iq)) -> sqrt(p^2+q^2): przepisujemy argumenty
wirt = wirt.replace(lambda e: e.func == sp.Pow and e.exp.is_Rational
                    and not e.base.is_Symbol,
                    lambda e: sp.Pow(sp.expand(e.base), e.exp))
check("P1c: dE/du* == (E_p + i E_q)/2 (Wirtinger => czynnik 1/2)",
      wirt,
      "i u_t = dE/du* => a_t = +L_minus b/2, b_t = -L_plus a/2;"
      " JL = [[0, L_minus/2], [-L_plus/2, 0]], lambda^2 = -nu/4,"
      " nu = spec(L_minus L_plus)")

# --- tozsamosc centralna: L_minus g0 = 0 on-shell EXACT ---
Lminus_g0 = sp.expand(Lminus_B.subs(B, g0).doit())
check("P1c: L_minus g0 == EL0 (tozsamosc strukturalna, off-shell)",
      Lminus_g0 - EL0)
onshell = sp.simplify(Lminus_g0.subs(g0pp, g0pp_onshell))
check("P1c-CENTRALNA: L_minus g_d == 0 ON-SHELL (exact; mod fazowy)",
      onshell, "FAIL => STOP calego cyklu (LOCK sec. 3 Phase 1)")
check("P1c: Q_minus == (K g0')'/g0 on-shell (dyskretyzacja on-shell"
      " z Phase_method_decisions sec. 2 jest tozsamosciowa)",
      sp.simplify((Qminus - sp.diff(Kg * g0p, x) / g0)
                  .subs(g0pp, g0pp_onshell)))

# --- rownowaznosc Q_plus z forma operatora Phase 3 cyklu bloch ---
Q_bloch = Kg * (2 * beta * g0 - 3 * gamma * g0 ** 2 + 2 * g0p ** 2 / g0 ** 2)
check("P1c: Q_plus == K(2 beta g - 3 gamma g^2 + 2 g'^2/g^2) on-shell"
      " (operator Phase 3 bloch)",
      sp.simplify((Qplus - Q_bloch).subs(g0pp, g0pp_onshell)))

# --- dyspersja JL wokol prozni u=1 (do C3) ---
print("\n[P1c-vac] dyspersja JL wokol prozni u=1 (beta=gamma):")
K1 = Kfun(1)
Upp1 = sp.diff(Ufun(zz), zz, 2).subs(zz, 1).subs(beta, gamma)
Up1 = sp.diff(Ufun(zz), zz).subs(zz, 1).subs(beta, gamma)
Lp_sym = K1 * k ** 2 + Upp1          # symbol L_plus na e^{ikx}
Lm_sym = K1 * k ** 2 + Up1 / 1       # symbol L_minus (Q_minus(1)=U'(1)/1)
lam2 = sp.simplify(-(Lp_sym * Lm_sym) / 4)
print("  symbol L_plus  =", sp.simplify(Lp_sym), " (k^2 - gamma)")
print("  symbol L_minus =", sp.simplify(Lm_sym), " (k^2; U'(1)=0)")
print("  lambda^2(k) = -L_plus*L_minus/4 =", lam2)
check("P1c-vac: lambda^2 == -k^2 (k^2 - gamma)/4",
      lam2 - (-(k ** 2) * (k ** 2 - gamma) / 4))
check("P1c-vac: U'(1) == 0 przy beta=gamma (proznia stacjonarna)", Up1)
lam_max = sp.sqrt(sp.Rational(1, 4) * (gamma / 2) * (gamma - gamma / 2))
print("  0 < k^2 < gamma: lambda^2 > 0 -- PROZNIA NIESTABILNA takze")
print("  symplektycznie; max Re lambda = gamma/4 przy k^2 = gamma/2")
check("P1c-vac: max Re lambda == gamma/4 (przy k^2=gamma/2)",
      sp.simplify(lam_max - gamma / 4))

# =========================================================== werdykt
print("\n" + "=" * 78)
if FAILS:
    print("PHASE 1: FAIL -- %s" % ", ".join(FAILS))
    print("P1c FAIL => STOP cyklu (LOCK sec. 3 Phase 1 / sec. 6).")
    sys.exit(1)
print("PHASE 1: PASS (wszystkie tozsamosci sympy-exact).")
print("Operatory do Phase 2/3 (FROZEN, Phase_method_decisions sec. 1-2):")
print("  L_plus  = -d/dx(K d/dx) + K(2 beta g - 3 gamma g^2 + 2 g'^2/g^2)")
print("  L_minus = -d/dx(K d/dx) + (K g')'/g   [L_minus g_d = 0 exact]")
print("  JL = [[0, L_minus/2], [-L_plus/2, 0]];  lambda^2 = -nu/4")
print("  C3 (analitycznie): lambda^2(kappa) = -kappa^2(kappa^2 - gamma)/4")
print("=" * 78)
