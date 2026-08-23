# Phase 0 check: numeryczna weryfikacja rownan #1-#3 LOCKa PRZED eksperymentem.
# Model ogona: s_i(x) = eta*exp(-m0*max(0, |x-c_i|-R)); J3 = int s1*s2*s3 /eta^3.
# Predykcje do zweryfikowania:
#   mu_pred = zbocze ln|delta_pred| vs g ~ (sqrt(3)-1)*m0 = 0.732 (trojkat)
#   |delta_chain| << |delta_tri| (Fermat za obiektem srodkowym)
#   rzad wielkosci C (raportowany, poza kryterium)

import numpy as np

m0 = 1.0
eta = np.sqrt(0.30)      # ~0.548
R = 8.0
b = 1.60
A2 = 3.50                # amplituda Yukawy 2-cialowej (poprzedni cykl, fit czulosci)
G_SCAN = [2.0, 2.5, 3.0, 4.0]

# drobna siatka lokalna do calki (niezalezna od siatki eksperymentu)
h = 0.05
ext = 40.0
xs = np.arange(-ext, ext + h / 2, h)
Xg, Yg = np.meshgrid(xs, xs, indexing="ij")


def tail(cx, cy):
    r = np.sqrt((Xg - cx) ** 2 + (Yg - cy) ** 2)
    return np.exp(-m0 * np.maximum(0.0, r - R))


def J3(centers):
    p = tail(*centers[0]) * tail(*centers[1]) * tail(*centers[2])
    return float(np.sum(p) * h * h)


print("=== Phase 0 check (op-nbody-additivity): rownania #1-#3 ===")
print(f"m0={m0} eta={eta:.4f} R={R} b={b} A2={A2} (z op-lock-interaction)")
print()
rows = []
for g in G_SCAN:
    d = 2 * R + g
    rho = d / np.sqrt(3.0)
    tri = [(rho * np.cos(a), rho * np.sin(a))
           for a in (np.pi / 2, np.pi / 2 + 2 * np.pi / 3, np.pi / 2 + 4 * np.pi / 3)]
    cha = [(-d, 0.0), (0.0, 0.0), (d, 0.0)]
    j_tri, j_cha = J3(tri), J3(cha)
    pairsum = 3.0 * A2 * np.exp(-m0 * g)          # trojkat: 3 pary przy szczelinie g
    delta_tri = -2.0 * b * eta**3 * j_tri / pairsum
    # lancuch: 2 pary przy g + 1 para przy 2g+2R (zaniedbywalna)
    pairsum_c = 2.0 * A2 * np.exp(-m0 * g) + A2 * np.exp(-m0 * (2 * g + 2 * R))
    delta_cha = -2.0 * b * eta**3 * j_cha / pairsum_c
    rows.append((g, j_tri, j_cha, delta_tri, delta_cha))
    print(f"  g={g:3.1f}: J3_tri={j_tri:.4e}  J3_chain={j_cha:.4e}  "
          f"delta_pred_tri={delta_tri:+.3e}  delta_pred_chain={delta_cha:+.3e}")
print()

gv = np.array([r[0] for r in rows])
dv = np.array([abs(r[3]) for r in rows])
Amat = np.vstack([np.ones(len(gv)), -gv]).T
coef, *_ = np.linalg.lstsq(Amat, np.log(dv), rcond=None)
mu_pred, C_pred = coef[1], np.exp(coef[0])
mu_theory = (np.sqrt(3.0) - 1.0) * m0
ok_mu = abs(mu_pred - mu_theory) <= 0.15
ratio_geom = max(abs(r[4]) / abs(r[3]) for r in rows)
ok_geom = ratio_geom < 0.05

print(f"ROWNANIE #2: mu_pred(fit modelu) = {mu_pred:.4f}  vs  (sqrt(3)-1)*m0 = "
      f"{mu_theory:.4f}   |roznica| <= 0.15: {'PASS' if ok_mu else 'FAIL'}")
print(f"  prefaktor C ~ {C_pred:.3e} -> |delta_tri(2)| ~ {C_pred*np.exp(-2*mu_pred):.2e} "
      f"(rzad; poza kryterium)")
print(f"ROWNANIE #3: max |delta_chain/delta_tri| = {ratio_geom:.2e} < 0.05: "
      f"{'PASS' if ok_geom else 'FAIL'}  (dyskryminator geometryczny)")
print(f"ZNAK (rownanie #1): delta_tri < 0 wszedzie: "
      f"{'PASS' if all(r[3] < 0 for r in rows) else 'FAIL'}")
print()
print(f"UWAGA pre-code: przewidywane |delta_tri| ~ 1e-3..1e-2 przy g=2, "
      f"malejace ~e^(-0.73g); punkty g=4 moga zejsc pod podloge 1e-4 -> "
      f"fit (c) na dostepnych punktach (zapisane w LOCK G4).")
allok = ok_mu and ok_geom and all(r[3] < 0 for r in rows)
print(f"PODSUMOWANIE PHASE 0 CHECK: {'PASS' if allok else 'FAIL'}")
