# -*- coding: utf-8 -*-
"""
Phase 1 (QF analitycznie) — op-substrate-fluctuation-channel-2026-08-23
LOCK: Phase0_balance.md (kryteria zamkniete PRZED kodem).

Wyprowadzenia exact (sympy):
  T1  F_fl(g) = (1/2) ln(1-g^2) < 0 dla g in (0,1)  [uniwersalne przyciaganie]
  T2  czesc klasyczna pinningu: wyraz wiodacy = -(v1*v2/G0^2)*Gd  [znak ~ v1*v2]
  T3  M-S (zrodla): F_int = -q1*q2*Gd  (dokladnie, dopelnienie kwadratu)
  T4  kontinuum: G = e^{-m r}/(4 pi r) spelnia -(G'' + 2G'/r) + m^2 G = 0
  T5  kontrola negatywna maszynerii: g>1 (niedodatnia macierz kowariancji)
      musi zostac WYKRYTA przez test
  T6  szereg F_fl: -g^2/2 - g^4/4 - ...  => zasieg 2m, potega -2 (sygnatura
      kanalu fluktuacyjnego vs m, -1 kanalu klasycznego)

Kazdy test ma osiagalny FAIL (forbidden move #5).
"""
import sympy as sp

results = {}

def check(name, cond, detail=""):
    results[name] = bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  {detail}")

print("=" * 72)
print("Phase 1 — QF analitycznie (sympy exact)")
print("=" * 72)

# ---------------------------------------------------------------- T1/T2
# Marginal 2-punktowy pola Gaussa: Sigma = [[G0, Gd],[Gd, G0]].
# Pinning phi(x1)=v1, phi(x2)=v2:
#   F(d) = (1/2) ln det Sigma + (1/2) v^T Sigma^{-1} v   (+ const)
# G0 > 0 (wariancja), 0 < g < 1 (dodatnia okreslonosc) — zalozenia konieczne,
# by sympy poprawnie upraszczal logarytmy (poprawka implementacyjna po
# 1. przebiegu: bez 'positive' T1a nie upraszczal sie mimo rownosci wyrazen).
G0 = sp.Symbol('G0', positive=True)
g = sp.Symbol('g', positive=True)
Gd, v1, v2 = sp.symbols('G_d v1 v2', real=True)
Sigma = sp.Matrix([[G0, Gd], [Gd, G0]])
v = sp.Matrix([v1, v2])

detS = sp.simplify(Sigma.det())
F_det = sp.Rational(1, 2) * sp.log(detS)
F_det_inf = sp.Rational(1, 2) * sp.log(G0**2)          # Gd(d->inf)=0
F_fl = sp.simplify(F_det - F_det_inf)                   # czesc fluktuacyjna
F_fl_g = sp.logcombine(sp.expand_log(F_fl.subs(Gd, g * G0), force=False),
                       force=False)
F_fl_expected = sp.Rational(1, 2) * sp.log(1 - g**2)
diff_T1a = sp.simplify(F_fl_g - F_fl_expected)
# rownowaznosc logarytmow przy G0>0, 0<g<1:
ok_T1a = sp.simplify(sp.expand_log(diff_T1a, force=True)) == 0
check("T1a: F_fl = (1/2)ln(1-g^2) (exact)", ok_T1a, f"F_fl = {F_fl_g}")

# znak: 0<g<1 => 0<1-g^2<1 => ln<0; siatka numeryczna (osiagalny FAIL)
grid = [sp.Rational(i, 100) for i in range(1, 100)]
signs_neg = all(sp.N(F_fl_expected.subs(g, gv)) < 0 for gv in grid)
# monotonicznosc w g: dF/dg = -g/(1-g^2) < 0 => |F| rosnie z g;
# poniewaz g(d) maleje z d (Yukawa), F_fl(d) rosnie do 0 => przyciaganie
dFdg = sp.simplify(sp.diff(F_fl_expected, g))
mono = all(sp.N(dFdg.subs(g, gv)) < 0 for gv in grid)
check("T1b: F_fl < 0 dla g in (0,1) [przyciaganie uniwersalne]", signs_neg)
check("T1c: dF_fl/dg < 0 (|F| rosnie z g; F_fl(d) rosnie z d)", mono)

# czesc klasyczna: (1/2) v^T Sigma^{-1} v - granica d->inf
F_cl = sp.Rational(1, 2) * (v.T * Sigma.inv() * v)[0, 0]
F_cl_inf = (v1**2 + v2**2) / (2 * G0)
F_cl_int = sp.simplify(F_cl - F_cl_inf)
series_cl = sp.series(F_cl_int, Gd, 0, 3).removeO().expand()
lead = series_cl.coeff(Gd, 1)
check("T2: wyraz wiodacy klasyczny = -(v1*v2/G0^2)*Gd",
      sp.simplify(lead + v1 * v2 / G0**2) == 0,
      f"coeff(Gd) = {lead}")

# ---------------------------------------------------------------- T3
# M-S: H = (1/2) x^T P x - q^T x, P = Sigma^{-1}.
# Dopelnienie kwadratu => F = -(1/2) q^T Sigma q => interakcja = -q1*q2*Gd.
q1, q2, x1, x2 = sp.symbols('q1 q2 x1 x2', real=True)
x = sp.Matrix([x1, x2]); qv = sp.Matrix([q1, q2])
P = Sigma.inv()
H = sp.Rational(1, 2) * (x.T * P * x)[0, 0] - (qv.T * x)[0, 0]
xstar = Sigma * qv                                     # minimum: x* = Sigma q
H_shift = sp.expand(H.subs({x1: xstar[0] + x1, x2: xstar[1] + x2}))
# po przesunieciu: H = (1/2) x^T P x - (1/2) q^T Sigma q  (czlon liniowy = 0)
lin_terms = sp.simplify(H_shift.coeff(x1, 1).subs({x2: 0})) == 0 and \
            sp.simplify(H_shift.coeff(x2, 1).subs({x1: 0})) == 0
F_source = sp.simplify(H_shift.subs({x1: 0, x2: 0}))
F_source_expected = -sp.Rational(1, 2) * (qv.T * Sigma * qv)[0, 0]
ok3 = lin_terms and sp.simplify(F_source - sp.expand(F_source_expected)) == 0
cross = sp.expand(F_source_expected).coeff(q1 * q2)
check("T3: F_zrodla = -(1/2)q^T Sigma q; czlon krzyzowy = -q1*q2*Gd",
      ok3 and sp.simplify(cross + Gd) == 0, f"czlon krzyzowy/(q1q2) = {cross}")

# ---------------------------------------------------------------- T4
r, m = sp.symbols('r m', positive=True)
Gcont = sp.exp(-m * r) / (4 * sp.pi * r)
lhs = sp.simplify(-(sp.diff(Gcont, r, 2) + 2 / r * sp.diff(Gcont, r)) + m**2 * Gcont)
check("T4: (-(d2/dr2 + (2/r)d/dr) + m^2) e^{-mr}/(4 pi r) = 0 (r>0)",
      sp.simplify(lhs) == 0, f"residuum = {sp.simplify(lhs)}")

# ---------------------------------------------------------------- T5
# kontrola negatywna maszynerii: g=1.2 => 1-g^2<0 => log niereal.
bad = F_fl_expected.subs(g, sp.Rational(12, 10))
is_detected = not sp.N(bad).is_real
check("T5: kontrola negatywna — g>1 wykryte jako niefizyczne", is_detected,
      f"F_fl(g=1.2) = {sp.N(bad)}")

# ---------------------------------------------------------------- T6
ser = sp.series(F_fl_expected, g, 0, 6).removeO()
check("T6a: F_fl = -g^2/2 - g^4/4 + O(g^6)",
      sp.simplify(ser - (-g**2 / 2 - g**4 / 4)) == 0, f"szereg = {ser}")
# sygnatura: g = G(d)/G(0), G ~ e^{-m d}/d  =>  F_fl ~ -(1/2G0^2) e^{-2md}/d^2
d_ = sp.symbols('d', positive=True)
F_asym = sp.simplify((-sp.Rational(1, 2)) * (Gcont.subs(r, d_) / G0)**2)
expected_asym = -sp.exp(-2 * m * d_) / (32 * sp.pi**2 * d_**2 * G0**2)
check("T6b: F_fl ~ -e^{-2 m d}/(32 pi^2 d^2 G0^2) [zasieg 2m, potega -2]",
      sp.simplify(F_asym - expected_asym) == 0, f"F_asym = {F_asym}")

# ---------------------------------------------------------------- SUMMARY
print("-" * 72)
npass = sum(results.values()); ntot = len(results)
print(f"SUMMARY Phase 1: {npass}/{ntot} PASS")
for k, ok in results.items():
    print(f"  {'PASS' if ok else 'FAIL'}  {k}")
print()
print("Wnioski analityczne (exact):")
print(" * kanal M-S (zrodlowy):   U = -q1 q2 G_m(d)      -> znak ~ sgn(q1 q2),")
print("   zasieg m, na krytycznosci ~ -1/d  (kanal LADUNKOWY)")
print(" * kanal M-D klasyczny:    ~ -(v1 v2/G0^2) G_m(d)  -> znak ~ sgn(v1 v2)")
print(" * kanal M-D fluktuacyjny: F_fl = (1/2)ln(1-g^2) < 0 ZAWSZE")
print("   (niezaleznie od v1, v2) -> zasieg 2m, na krytycznosci ~ -1/d^2")
print("   JEDYNY kanal o uniwersalnym znaku przyciagania na poziomie 0.")
if npass != ntot:
    raise SystemExit(1)
