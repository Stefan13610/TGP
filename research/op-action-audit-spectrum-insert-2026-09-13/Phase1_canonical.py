#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-action-audit-spectrum-insert (Phase 1) -- kanonika z JEDNEJ akcji.

LOCK: Phase0_balance.md sec. 2 Phase 1; MD (FROZEN) sec. 1, 3, 5.
- P1a: z literalnej akcji (eq:S-TGP-unified-M911-canonical; formy
  CYTAT MD sec. 1): L = 1/2 M psidot^2 - 1/2 Keff |grad psi|^2 - Ueff,
  M(psi) = WYNIK: sqrt(-g) * K * |g^tt|  [czlon czasowy LOCK sec. 1],
  Keff = K_geo psi^4 (odczyt B DZIEDZICZONY, MD poprzednika sec. 2),
  Ueff = w * V_M911 (tozsamosc wielomianowa); pi, H, EOM z M', K'.
- P1b (gate, FAIL => STOP): dH/dpsi|_{pi=0} == dE_PRIMARY/dpsi
  (sympy simplify = 0); dodatkowo Ueff == w*V (simplify = 0).
- P1c (gate implementacji): M, Keff, Ueff w {0.5, 1, 7/6, 1.3}
  sympy vs float -- zgodnosc 1e-12; POST-CORRECTION-1
  (Phase_correction_note_P1c_gate.md): sympy ewaluowane w IDENTYCZNYM
  wejsciu binarnym double (izolacja bledu implementacji -- cel gate'u);
  referencja float M uzywa math.fma dla (4-3psi); wartosci w doklad-
  nych punktach wymiernych raportowane rownolegle (deskryptywnie).
- Deskryptywnie (MD sec. 3 kaweat, NIE werdyktotworcze): wariant
  ZNAKOWANEGO g^tt (sygnatura (-,+,+,+) literalnie): L_sgn, EOM,
  statyka vs R3 ODE.

REJESTR WEJSC [INPUT]: K_geo = gamma = c0 = 1 [LOCK sec. 1].
Formy CYTAT: w = psi/(4-3psi) [eq:vol-element-M911, sek08c ~510-519];
V = -gamma psi^2 (4-3psi)^2 / 12 [eq:V-M911, sek08a ~977-981];
K = K_geo psi^4 [eq:K-coupling-unified, sek08a ~169-172];
ds^2 = -c0^2(4-3psi)/psi dt^2 + psi/(4-3psi) dx^2
[eq:metric-M911-canonical, sek08c ~420-424] => g^tt = -psi/(c0^2(4-3psi)),
g^ij = (4-3psi)/psi delta^ij; sqrt(-g) = c0 psi/(4-3psi).
"""
import sympy as sp

OUT = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
       "op-action-audit-spectrum-insert-2026-09-13/Phase1_output.txt")

lines = []


def em(s=""):
    print(s, flush=True)
    lines.append(s)


# ---------------------------------------------------------------- symbole
psi = sp.Symbol('psi', positive=True)
gamma, Kgeo, c0 = sp.Integer(1), sp.Integer(1), sp.Integer(1)  # INPUT LOCK

em("=" * 78)
em("PHASE 1 -- kanonika z jednej akcji (sympy); K_geo=gamma=c0=1 [INPUT]")
em("Formy CYTAT (MD sec.1): w=psi/(4-3psi) [eq:vol-element-M911];")
em("  V=-psi^2(4-3psi)^2/12 [eq:V-M911]; K=psi^4 [eq:K-coupling-unified];")
em("  ds^2=-(4-3psi)/psi dt^2 + psi/(4-3psi) dx^2 [eq:metric-M911-canonical]")
em("Odczyt kinetyki przestrzennej: B DZIEDZICZONY (MD poprzednika sec.2):")
em("  w * g^ij = [psi/(4-3psi)]*[(4-3psi)/psi] = 1 => Keff = K_geo psi^4")
em("=" * 78)

# formy z korpusu (CYTATY)
w = psi / (4 - 3 * psi)                       # eq:vol-element-M911
V = -gamma * psi ** 2 * (4 - 3 * psi) ** 2 / 12   # eq:V-M911
K = Kgeo * psi ** 4                            # eq:K-coupling-unified
sqrtg = c0 * psi / (4 - 3 * psi)               # eq:vol-element-M911
gtt_signed = -psi / (c0 ** 2 * (4 - 3 * psi))  # z eq:metric-M911-canonical
gtt_abs = psi / (c0 ** 2 * (4 - 3 * psi))      # |g^tt| na (0,4/3) (LOCK s.1)
gij = (4 - 3 * psi) / psi                      # inwersja przestrzenna

# ------------------------------------------------------------- P1a: M, K, U
em()
em("P1a: WYPROWADZENIE M, Keff, Ueff (M jest WYNIKIEM, nie wejsciem)")
M = sp.simplify(sqrtg * K * gtt_abs)           # czlon czasowy LOCKa sec.1
M_expected = psi ** 6 / (4 - 3 * psi) ** 2
em("  M(psi) = sqrt(-g)*K*|g^tt| = %s" % M)
em("  simplify(M - psi^6/(4-3psi)^2) = %s"
   % sp.simplify(M - M_expected))
Keff = sp.simplify(sqrtg * K * gij * w / sqrtg * 1)  # zapis odczytu B:
# w*g^ij == 1 dokladnie -> Keff = K_geo psi^4; weryfikacja tozsamosci:
id_B = sp.simplify(w * gij - 1)
em("  odczyt B: simplify(w*g^ij - 1) = %s  (== 0 wymagane)" % id_B)
Keff = K                                        # = K_geo psi^4 (odczyt B)
em("  Keff(psi) = %s" % Keff)
Ueff = sp.expand(sp.simplify(w * V))
Ueff_poly = gamma * (psi ** 4 / 4 - psi ** 3 / 3)
id_U = sp.simplify(Ueff - Ueff_poly)
em("  Ueff = w*V = %s ; simplify(Ueff - gamma(psi^4/4-psi^3/3)) = %s"
   % (Ueff, id_U))
Mp = sp.simplify(sp.diff(M, psi))
Kp = sp.diff(Keff, psi)
Up = sp.expand(sp.diff(Ueff, psi))
em("  M'(psi)  = %s" % sp.factor(Mp))
em("  Keff'    = %s" % Kp)
em("  Ueff'    = %s  (= gamma psi^2 (psi-1): %s)"
   % (Up, sp.simplify(Up - gamma * psi ** 2 * (psi - 1))))

# ------------------------------------------- L, pi, H, EOM (pole psi(t,x))
em()
em("Gestosc Lagranzjanu (forma kanoniczna LOCKa sec.2 P1a):")
em("  L = 1/2 M(psi) psidot^2 - 1/2 Keff(psi) |grad psi|^2 - Ueff(psi)")
t, x = sp.symbols('t x')
f = sp.Function('psi')(t, x)
ft, fx = sp.diff(f, t), sp.diff(f, x)
Mf = M.subs(psi, f)
Kf = Keff.subs(psi, f)
Uf = Ueff.subs(psi, f)
L = sp.Rational(1, 2) * Mf * ft ** 2 - sp.Rational(1, 2) * Kf * fx ** 2 - Uf

pi_expr = sp.diff(L, ft)
em("  pi = dL/dpsidot = M(psi) psidot   [= %s]" % (M * sp.Symbol('psidot')))
# Hamiltonian (gestosc): H = pi*psidot - L, psidot = pi/M
piS = sp.Symbol('pi')
H_dens = sp.simplify((piS * (piS / Mf) - L.subs(ft, piS / Mf)))
H_expected = (piS ** 2 / (2 * Mf) + sp.Rational(1, 2) * Kf * fx ** 2 + Uf)
em("  H[pi,psi] (gestosc) = pi^2/(2M) + 1/2 Keff |grad psi|^2 + Ueff")
em("  simplify(H - powyzsze) = %s" % sp.simplify(H_dens - H_expected))

# EOM Euler-Lagrange
EL = (sp.diff(sp.diff(L, ft), t) + sp.diff(sp.diff(L, fx), x)
      - sp.diff(L, f))
EOM_target = (Mf * sp.diff(f, t, 2)
              + sp.Rational(1, 2) * Mp.subs(psi, f) * ft ** 2
              - sp.diff(Kf * fx, x)
              + sp.Rational(1, 2) * Kp.subs(psi, f) * fx ** 2
              + Up.subs(psi, f))
em()
em("EOM (Euler-Lagrange, czlony M' i K' jawnie):")
em("  M psiddot + 1/2 M' psidot^2 = div(Keff grad psi)")
em("    - 1/2 Keff' |grad psi|^2 - Ueff'   [= -dE_PRIMARY/dpsi]")
em("  weryfikacja: simplify(EL - forma_docelowa) = %s"
   % sp.simplify(EL - EOM_target))

# --------------------------------------------------------- P1b: gate STOP
em()
em("P1b (GATE, FAIL => STOP): dH/dpsi|_{pi=0} == dE_PRIMARY/dpsi")
# wariacja H przy pi=0 (gestosc h = 1/2 Keff fx^2 + Ueff):
h0 = sp.Rational(1, 2) * Kf * fx ** 2 + Uf
dH0 = sp.diff(h0, f) - sp.diff(sp.diff(h0, fx), x)
# wariacja E_PRIMARY poprzednika: E = int 1/2 psi^4 |grad|^2 + Ueff
e_prim = sp.Rational(1, 2) * (f ** 4) * fx ** 2 + Uf
dE = sp.diff(e_prim, f) - sp.diff(sp.diff(e_prim, fx), x)
resid = sp.simplify(dH0 - dE)
em("  dH/dpsi|_{pi=0} = -(Keff psi_x)_x + 1/2 Keff' psi_x^2 + Ueff'")
em("  simplify(dH/dpsi|_{pi=0} - dE_PRIMARY/dpsi) = %s" % resid)
p1b_a = (resid == 0)
p1b_b = (id_U == 0)
p1b_c = (id_B == 0)
em("  tozsamosc wielomianowa Ueff = w*V: %s" % ("PASS" if p1b_b else "FAIL"))
em("  tozsamosc odczytu B  w*g^ij = 1 : %s" % ("PASS" if p1b_c else "FAIL"))
P1B = p1b_a and p1b_b and p1b_c
em("  P1b: %s" % ("PASS" if P1B else "FAIL => STOP"))

# --------------------------------------------------------- P1c: sympy/float
em()
em("P1c (GATE): M, Keff, Ueff w {0.5, 1, 7/6, 1.3} sympy vs float 1e-12")
em("  [POST-CORRECTION-1: identyczne wejscie binarne double; referencja")
em("   float M z math.fma dla (4-3psi) -- correction note zapisana PRZED")
em("   uzyciem wyniku; pierwotny output: Phase1_output_pre_correction1.txt]")
import math


def M_float(p):
    q = math.fma(-3.0, p, 4.0)      # (4-3p) pojedyncze zaokraglenie
    return p ** 6 / (q * q)


def Keff_float(p):
    return p ** 4


def Ueff_float(p):
    return p ** 4 / 4.0 - p ** 3 / 3.0


pts_exact = [sp.Rational(1, 2), sp.Integer(1), sp.Rational(7, 6),
             sp.Rational(13, 10)]
P1C = True
for pt_ex in pts_exact:
    fl = float(pt_ex)
    pt_bin = sp.Rational(fl)        # DOKLADNA konwersja wejscia double
    for name, expr, impl in (("M", M, M_float), ("Keff", Keff, Keff_float),
                             ("Ueff", Ueff, Ueff_float)):
        ex = float(sp.N(expr.subs(psi, pt_bin), 30))
        dv = abs(ex - impl(fl))
        ok = dv <= 1e-12
        P1C = P1C and ok
        ex_rat = float(sp.N(expr.subs(psi, pt_ex), 30))
        em("  psi=%-6s %-4s sympy(bin)=%+.16e float=%+.16e |d|=%.2e %s"
           "  [deskr. sympy(exact %s)=%+.16e]"
           % (str(pt_ex), name, ex, impl(fl), dv,
              "PASS" if ok else "FAIL", str(pt_ex), ex_rat))
em("  P1c: %s" % ("PASS" if P1C else "FAIL"))

# ------------------------------------- warunek STOP LOCKa (procedura wariacji)
em()
em("Warunek STOP LOCKa sec.2 Phase 1 (g_eff jako zmienna niezalezna /")
em("  wariacja euklidesowa): NIE AKTYWOWANY -- akcja w formie zlozonej")
em("  (metryka podstawiona jako funkcja psi PRZED wariacja; wariacja")
em("  standardowa Lorentzowska pola psi; MD sec. 3).")

# ------------------- deskryptywnie: wariant ZNAKOWANEGO g^tt (MD sec.3 kaweat)
em()
em("-" * 78)
em("DESKRYPTYWNIE (MD sec.3, NIE werdyktotworcze): wariant ZNAKOWANEGO")
em("g^tt = -psi/(c0^2(4-3psi)) (sygnatura (-,+,+,+) literalnie w")
em("kontrakcji +1/2 K g^{mu nu} d_mu psi d_nu psi):")
M_sgn = sp.simplify(sqrtg * K * gtt_signed)
em("  wspolczynnik czasowy = sqrt(-g)*K*g^tt = %s  (< 0 na (0,4/3))"
   % M_sgn)
em("  => L_sgn = -1/2 M psidot^2 + 1/2 Keff |grad psi|^2 - Ueff")
em("     (M = psi^6/(4-3psi)^2 jak wyzej); duch w czlonie czasowym")
L_sgn = (-sp.Rational(1, 2) * Mf * ft ** 2
         + sp.Rational(1, 2) * Kf * fx ** 2 - Uf)
EL_sgn = (sp.diff(sp.diff(L_sgn, ft), t) + sp.diff(sp.diff(L_sgn, fx), x)
          - sp.diff(L_sgn, f))
em("  EOM_sgn: %s = 0" % sp.simplify(-EL_sgn))
# statyka obu wariantow radialnie vs R3 ODE
r = sp.Symbol('r', positive=True)
g_r = sp.Function('psi')(r)
gr, grr = sp.diff(g_r, r), sp.diff(g_r, r, 2)
# kanoniczny (LOCK): statyka = dE_PRIMARY/dpsi = 0
stat_can = (-(1 / r ** 2) * sp.diff(r ** 2 * Keff.subs(psi, g_r) * gr, r)
            + sp.Rational(1, 2) * Kp.subs(psi, g_r) * gr ** 2
            + Up.subs(psi, g_r))
lhs_R3 = grr + 2 / r * gr + 2 * gr ** 2 / g_r
rhs_can = sp.simplify(sp.expand(lhs_R3 - stat_can / Keff.subs(psi, g_r)
                                - lhs_R3) * (-1))
em()
em("  statyka radialna (po podzieleniu przez Keff):")
em("    kanoniczna (LOCK, dE/dpsi=0):  psi''+2/r psi'+2psi'^2/psi = RHS,")
resid_can = sp.simplify(stat_can / Keff.subs(psi, g_r)
                        - (lhs_R3 - (g_r - 1) / g_r ** 2) * (-1))
# sprawdz: stat_can/Keff == -(lhs_R3 - (psi-1)/psi^2) ?
chk_can = sp.simplify(stat_can / Keff.subs(psi, g_r)
                      + lhs_R3 - (g_r - 1) / g_r ** 2)
em("      RHS = (psi-1)/psi^2 ;  weryfikacja simplify = %s" % chk_can)
# znakowany: statyka = wariacja int[1/2 K psi'^2 - Ueff] = 0
stat_sgn = (-(1 / r ** 2) * sp.diff(r ** 2 * Keff.subs(psi, g_r) * gr, r)
            + sp.Rational(1, 2) * Kp.subs(psi, g_r) * gr ** 2
            - Up.subs(psi, g_r))
chk_sgn = sp.simplify(stat_sgn / Keff.subs(psi, g_r)
                      + lhs_R3 - (1 - g_r) / g_r ** 2)
em("    znakowana (literalna (-,+,+,+)): RHS = (1-psi)/psi^2 = R3 ODE")
em("      [eq:R3-ODE, dowod prop:V-M911-canonical];")
em("      weryfikacja simplify = %s" % chk_sgn)
em("  NAPIECIE ZNAKU (raportowane wprost, user-gate w NEEDS; zero")
em("  samowolnych napraw): statyka formy kanonicznej |g^tt| daje RHS")
em("  (psi-1)/psi^2 (zaniki Yukawy wokol prozni), a R3 ODE korpusu ma")
em("  (1-psi)/psi^2 (ogony oscylacyjne) -- odpowiada wariacji znakowanej,")
em("  ktora w sektorze czasowym daje ducha/tachion (audyt rozdz. 2).")

# ------------------------------------------------------------------ werdykt
em()
em("=" * 78)
em("PHASE 1 PODSUMOWANIE: P1a WYPROWADZONE (M = psi^6/(4-3psi)^2 WYNIK);")
em("  P1b %s; P1c %s" % ("PASS" if P1B else "FAIL => STOP",
                         "PASS" if P1C else "FAIL"))
em("  M(psi)   = K_geo psi^6/(c0 (4-3psi)^2)   [WYNIK, czlon |g^tt|]")
em("  Keff(psi)= K_geo psi^4                   [odczyt B DZIEDZICZONY]")
em("  Ueff(psi)= gamma(psi^4/4 - psi^3/3)      [= w*V, tozsamosc]")
em("  pi = M psidot;  H = int[pi^2/(2M) + 1/2 Keff|grad|^2 + Ueff]d3x")
em("  EOM: M psiddot + 1/2 M' psidot^2 = div(Keff grad psi)")
em("       - 1/2 Keff'|grad psi|^2 - Ueff'")
em("=" * 78)

with open(OUT, "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("zapisano:", OUT)
if not (P1B and P1C):
    raise SystemExit(1)
