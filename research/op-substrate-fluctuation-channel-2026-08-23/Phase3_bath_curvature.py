# -*- coding: utf-8 -*-
"""
Phase 3 (QB) — op-substrate-fluctuation-channel-2026-08-23
LOCK: Phase0_balance.md §3.2, §4 (QB-1..QB-3).

MFT substratu v2 (dodatekB rem:B-v2-status): wezel s w kapieli zamrozonych
sasiadow s_b, koordynacja z=6:
  V_eff(s; s_b) = (m0^2/2) s^2 + (lam0/4) s^4 + z J s^2 s_b^2 (s^2 - s_b^2)^2
Badana wielkosc: C(s_b) = d^2 V_eff / ds^2 |_{s=s_b}  (lokalna stabilnosc
konfiguracji jednorodnej o gestosci Phi_b = s_b^2).

QB-1: znak wkladu wiazania Delta_C_bond (sympy exact + kontrola numeryczna
      + kontrola negatywna z odwroconym znakiem J).
QB-2: prog rozrzedzenia Phi_c przy punkcie WF (r*=-2.251, u*=3.917, z=6,
      eq:B-WF dodatekB) + mapa na skanie r in [-3,-1], u in [2,6] (41x41).
QB-3: odczyt deskryptywny (zero claimow poziomu 1).
"""
import numpy as np
import sympy as sp

results = {}
def check(name, cond, detail=""):
    results[name] = bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  {detail}")

print("=" * 72)
print("Phase 3 — QB: krzywizna MFT v2 w kapieli")
print("=" * 72)

s, sb = sp.symbols('s s_b', positive=True)
m02, lam0, J, z = sp.symbols('m0^2 lambda0 J z', real=True)

V_onsite = m02 / 2 * s**2 + lam0 / 4 * s**4
V_bond = z * J * s**2 * sb**2 * (s**2 - sb**2) ** 2
V_eff = V_onsite + V_bond

# ------------------------------------------------------------------ QB-1
dC_bond = sp.simplify(sp.diff(V_bond, s, 2).subs(s, sb))
expected_bond = 8 * z * J * sb**6
check("QB-1a: Delta_C_bond = 8 z J s_b^6 (sympy exact)",
      sp.simplify(dC_bond - expected_bond) == 0, f"Delta_C_bond = {dC_bond}")

# kontrola numeryczna: 2. roznica centralna w losowych punktach
rng = np.random.default_rng(20260823)
ok_num = True
f_num = sp.lambdify((s, sb), V_bond.subs({z: 6, J: 1.0}), 'numpy')
for _ in range(20):
    sbv = float(rng.uniform(0.05, 1.5)); h = 1e-5
    num2 = (f_num(sbv + h, sbv) - 2 * f_num(sbv, sbv) + f_num(sbv - h, sbv)) / h**2
    exact = 8 * 6 * 1.0 * sbv**6
    ok_num &= abs(num2 - exact) / max(abs(exact), 1e-12) < 1e-4
check("QB-1b: kontrola numeryczna 2. roznicy (20 punktow, tol 1e-4)", ok_num)

# znak dla wszystkich s_b >= 0 przy z,J>0:
ok_sign = sp.simplify(expected_bond.subs({z: 6, J: 1})) == 8 * 6 * sb**6
check("QB-1c: Delta_C_bond >= 0 dla wszystkich s_b (z,J>0) — wiazanie "
      "gradientowe STABILIZUJE punkt jednorodny", ok_sign,
      "(parzysta potega s_b, wsp. dodatni)")

# kontrola negatywna: bond z odwroconym znakiem musi dac ujemny wklad
dC_flip = sp.simplify(sp.diff(-V_bond, s, 2).subs(s, sb)).subs({z: 6, J: 1,
                                                                sb: 0.7})
check("QB-1d: kontrola negatywna — flip znaku J wykrywany jako "
      "destabilizujacy", float(dC_flip) < 0, f"Delta_C(-J, s_b=0.7) = {dC_flip}")

# ------------------------------------------------------------------ QB-2
# C(s_b)/J = r + 3 u sigma + 8 z sigma^3,  sigma = s_b^2, r=m0^2/J, u=lam0/J
r_, u_, sig = sp.symbols('r u sigma', real=True)
C_red = r_ + 3 * u_ * sig + 8 * 6 * sig**3          # z=6
dC_dsig = sp.diff(C_red, sig)
# korekta implementacyjna po 1. przebiegu: oczekiwana stala to
# d(48 sigma^3)/dsigma = 144 sigma^2 (w 1. przebiegu blednie wpisano 288 —
# literowka w TEScie, nie w pochodnej; wniosek o jedynosci korzenia bez zmian)
check("QB-2a: dC/dsigma = 3u + 144 sigma^2 > 0 (u>0) — korzen jedyny",
      sp.simplify(dC_dsig - (3 * u_ + 144 * sig**2)) == 0)

rWF, uWF = -2.251, 3.917
sig_vac = abs(rWF) / uWF
poly = np.array([48.0, 0.0, 3 * uWF, rWF])           # 48 s^3 + 3u s + r
roots = np.roots(poly)
real_roots = [rt.real for rt in roots if abs(rt.imag) < 1e-12 and rt.real > 0]
sig_c = min(real_roots)
resid = 48 * sig_c**3 + 3 * uWF * sig_c + rWF
ratio = sig_c / sig_vac
sig_c0 = abs(rWF) / (3 * uWF)                        # prog bez wiazania
C_at_vac = rWF + 3 * uWF * sig_vac + 48 * sig_vac**3
print(f"    WF: sigma_vac = |r|/u = {sig_vac:.4f}; sigma_c = {sig_c:.4f} "
      f"(residuum {resid:.2e}); sigma_c^0 (bez bondu) = {sig_c0:.4f}")
print(f"    Phi_c/Phi_vac = {ratio:.4f};  C(sigma_vac) = {C_at_vac:.4f}")
check("QB-2b: istnieje prog rozrzedzenia Phi_c in (0, Phi_vac) przy WF",
      0 < sig_c < sig_vac and abs(resid) < 1e-12,
      f"Phi_c/Phi_vac = {ratio:.4f}")
check("QB-2c: kapiel o gestosci prozniowej STABILNA (C(sigma_vac) > 0)",
      C_at_vac > 0, f"C = {C_at_vac:.4f} = 2|r| + 48 sigma_vac^3")
check("QB-2d: bond ZAWEZA obszar niestabilny (sigma_c < sigma_c^0)",
      sig_c < sig_c0, f"{sig_c:.4f} < {sig_c0:.4f}")

# skan r in [-3,-1], u in [2,6], 41x41
rs = np.linspace(-3, -1, 41); us = np.linspace(2, 6, 41)
ratios = np.zeros((41, 41))
for i, rv in enumerate(rs):
    for j, uv in enumerate(us):
        rr = np.roots([48.0, 0.0, 3 * uv, rv])
        pos = [x.real for x in rr if abs(x.imag) < 1e-12 and x.real > 0]
        ratios[i, j] = min(pos) / (abs(rv) / uv)
ok_scan = np.all((ratios > 0) & (ratios < 1))
print(f"    skan 41x41: Phi_c/Phi_vac in [{ratios.min():.4f}, {ratios.max():.4f}]")
check("QB-2e: prog istnieje i Phi_c < Phi_vac na CALYM skanie", ok_scan)

# ---------------------------------------------------------------- SUMMARY
print("-" * 72)
npass = sum(results.values()); ntot = len(results)
print(f"SUMMARY Phase 3: {npass}/{ntot} PASS")
for k2, ok in results.items():
    print(f"  {'PASS' if ok else 'FAIL'}  {k2}")
print()
print("WERDYKT QB-1: wklad wiazania gradientowego do krzywizny w punkcie")
print("  jednorodnym = +8zJ s_b^6 >= 0 ZAWSZE. Znak tachionowy NIE moze")
print("  emergowac z wiazania gradientowego substratu na poziomie MFT —")
print("  jesli istnieje, pochodzi z czesci ON-SITE (rozrzedzenie).")
print("WERDYKT QB-2: prog rozrzedzenia ISTNIEJE: przy punkcie WF")
print(f"  Phi_c/Phi_vac = {ratio:.3f}; kapiel pelna stabilna, kapiel")
print(f"  rozrzedzona ponizej ~{100*ratio:.0f}% gestosci prozniowej -> krzywizna")
print("  tachionowa (spinodala substratu). Na skanie (r,u) prog zawsze")
print(f"  w (0,1): [{ratios.min():.3f}, {ratios.max():.3f}].")
print("QB-3 (deskryptywnie, ZERO claimow poziomu 1): wzor {pelna kapiel ->")
print("  stabilnie; rozrzedzenie -> tachionowo} jest strukturalnie rownolegly")
print("  do hipotezy dwoch sektorow N2 (izolowany/stabilny vs kapiel/tachion),")
print("  ale z ODWROTNA rola gestosci: na poziomie 0 to NISKA gestosc tla")
print("  destabilizuje. Porownanie czysto deskryptywne; poziom 1 rozstrzyga")
print("  op-bath-two-sectors (Q2), nie ten cykl.")
if npass != ntot:
    raise SystemExit(1)
