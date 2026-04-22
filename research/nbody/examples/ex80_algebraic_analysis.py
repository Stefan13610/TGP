"""
EX80: PRECYZYJNA ANALIZA ALGEBRAICZNA alpha*_1, alpha*_2
STATUS: LEGACY-TRANSLATIONAL

This script belongs to the older `alpha*` algebra / sum-rule program. It is
useful as historical exploratory context, but it is not part of the canonical
current `nbody` path.

Wersja szybka: uzywamy wartosci z ex79 (brentq xtol=1e-8):
  alpha*_1 = 2.44143051
  alpha*_2 = 2.74175411
  alpha_min = 2.56330 (minimize_scalar xatol=1e-5)
  f_min     = 197.9524

Pytania:
  1. Suma S = 2pi - 11/10 (ile ppm)?
  2. Roznica D ~ 3/10 (ile ppm)?
  3. Iloczyn P ~ (pi-7/10)(pi-4/10) (ile ppm)?
  4. Korekta eps = alpha*_1 - (pi-7/10) — co to jest?
  5. Wzory dla alpha_min i f_min
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit, minimize_scalar
import warnings
warnings.filterwarnings('ignore')

PHI     = (1.0 + np.sqrt(5.0)) / 2.0
R21_EXP = 206.768

# ----------------------------------------------------------------
# Znane wartosci z ex79 (brentq xtol=1e-8, rtol=1e-10)
# ----------------------------------------------------------------
A1_KNOWN = 2.44143051
A2_KNOWN = 2.74175411
AMIN_KNOWN = 2.56330    # z minimize_scalar xatol=1e-5
FMIN_KNOWN = 197.9524

if __name__ == '__main__':
    print("=" * 72)
    print("EX80: ANALIZA ALGEBRAICZNA alpha*_1, alpha*_2  (wersja szybka)")
    print("=" * 72)
    print(f"  alpha*_1 = {A1_KNOWN}  (z ex79 brentq xtol=1e-8)")
    print(f"  alpha*_2 = {A2_KNOWN}  (z ex79 brentq xtol=1e-8)")
    print()

    a1 = A1_KNOWN; a2 = A2_KNOWN
    S  = a1 + a2
    P  = a1 * a2
    D  = a2 - a1
    a_min = AMIN_KNOWN
    f_min = FMIN_KNOWN

    # ================================================================
    # SUMA
    # ================================================================
    print("=" * 72)
    print("SUMA  S = alpha*_1 + alpha*_2")
    print("=" * 72)
    ref_S = 2*np.pi - 11/10
    print(f"  S              = {S:.10f}")
    print(f"  2*pi - 11/10   = {ref_S:.10f}")
    diff_S = S - ref_S
    print(f"  diff S-ref     = {diff_S:+.4e}   ({abs(diff_S)/S*1e6:.2f} ppm)")
    print()

    # Inne kandydaty na sume
    cands_S = [
        ('2*pi - 11/10',         2*np.pi - 11/10),
        ('2*pi - ln3',           2*np.pi - np.log(3)),
        ('2*pi - 1.1',           2*np.pi - 1.1),
        ('2*pi - pi-2',          2*np.pi - (np.pi-2)),
        ('2*(pi-11/20)',         2*(np.pi-0.55)),
        ('2*pi - e/pi',          2*np.pi - np.e/np.pi),
        ('pi+e-1/phi',           np.pi+np.e-1/PHI),
    ]
    for name, val in sorted(cands_S, key=lambda c: abs(c[1]-S))[:5]:
        print(f"    {name:<30} = {val:.10f}  diff={S-val:+.2e}  ({abs(S-val)/S*1e6:.2f} ppm)")
    print()

    # ================================================================
    # ROZNICA
    # ================================================================
    print("=" * 72)
    print("ROZNICA  D = alpha*_2 - alpha*_1")
    print("=" * 72)
    print(f"  D              = {D:.10f}")
    print(f"  3/10           = 0.3000000000")
    print(f"  D - 3/10       = {D-0.3:+.4e}   ({abs(D-0.3)/D*1e6:.1f} ppm)")
    print()
    cands_D = [
        ('3/10',              3/10),
        ('e/9',               np.e/9),
        ('ln(4/3)',           np.log(4/3)),
        ('pi/10-1/100',       np.pi/10-0.01),
        ('(2pi-e)/10',        (2*np.pi-np.e)/10),
        ('pi^2/33',           np.pi**2/33),
        ('1/phi^3',           1/PHI**3),
        ('3/10+1/3000',       3/10+1/3000),
        ('3/10+1/3100',       3/10+1/3100),
        ('3/10+pi/10000',     3/10+np.pi/10000),
        ('1/(pi+e)',          1/(np.pi+np.e)),
        ('1/pi-1/30',         1/np.pi-1/30),
        ('2/(pi+e+pi)',       2/(np.pi+np.e+np.pi)),
    ]
    print("  Top kandydaci na D:")
    for name, val in sorted(cands_D, key=lambda c: abs(c[1]-D))[:8]:
        print(f"    {name:<28} = {val:.10f}  diff={D-val:+.2e}  ({abs(D-val)/D*1e6:.1f} ppm)")
    print()

    # ================================================================
    # ILOCZYN
    # ================================================================
    print("=" * 72)
    print("ILOCZYN  P = alpha*_1 * alpha*_2")
    print("=" * 72)
    print(f"  P              = {P:.10f}")
    P_hyp = (np.pi-0.7)*(np.pi-0.4)
    print(f"  (pi-7/10)(pi-4/10) = pi^2 - 11pi/10 + 28/100")
    print(f"                     = {P_hyp:.10f}")
    print(f"  P - P_hyp          = {P-P_hyp:+.4e}   ({abs(P-P_hyp)/P*1e6:.1f} ppm)")
    print()

    cands_P = [
        ('(pi-7/10)(pi-4/10)',      (np.pi-0.7)*(np.pi-0.4)),
        ('pi^2-11pi/10+28/100',     np.pi**2-1.1*np.pi+0.28),
        ('2pi+0.41',                2*np.pi+0.41),
        ('2pi+2/5',                 2*np.pi+0.4),
        ('2pi+phi/pi',              2*np.pi+PHI/np.pi),
        ('2pi+ln(3/2)',             2*np.pi+np.log(1.5)),
        ('2pi+phi/4',               2*np.pi+PHI/4),
        ('2pi+2/pi',                2*np.pi+2/np.pi),
        ('3pi-e',                   3*np.pi-np.e),
        ('e+pi',                    np.e+np.pi),
        ('phi*e',                   PHI*np.e),
        ('pi*phi+1',                np.pi*PHI+1),
        ('2pi+1/phi^2',             2*np.pi+1/PHI**2),
        ('2pi+3/e^2',               2*np.pi+3/np.e**2),
        ('(2pi-11/10)^2/4-9/400',   (2*np.pi-1.1)**2/4-9/400),
        ('(2pi-11/10)^2/4-D^2/4',   (2*np.pi-1.1)**2/4-D**2/4),
        ('pi*(pi-11/10)+2pi/10*8',  np.pi*(np.pi-1.1)+2*np.pi/10*8),
        ('(S/2)^2 - (D/2)^2',       (S/2)**2-(D/2)**2),   # = P, zawsze
        ('(pi-11/20)^2 - (D/2)^2',  (np.pi-0.55)**2-(D/2)**2),
        ('S^2/4 - D^2/4',           S**2/4-D**2/4),        # = P, zawsze
    ]
    print("  Top kandydaci na P:")
    for name, val in sorted(cands_P, key=lambda c: abs(c[1]-P))[:12]:
        ppm = abs(P-val)/P*1e6
        if ppm < 0.01: marker = "  *** BARDZO DOKLADNE"
        elif ppm < 5:  marker = "  ** dokladne"
        elif ppm < 50: marker = "  * ok"
        else: marker = ""
        print(f"    {name:<40} = {val:.10f}  diff={P-val:+.3e}  ({ppm:.1f} ppm){marker}")
    print()

    # ================================================================
    # KOREKTA EPSILON
    # ================================================================
    print("=" * 72)
    print("KOREKTA  eps = alpha*_i - (pi - k/10)")
    print("=" * 72)
    eps1 = a1 - (np.pi - 7/10)   # powinno byc ~0 jesli hipoteza jest dobra
    eps2 = a2 - (np.pi - 4/10)
    print(f"  eps1 = a1 - (pi-7/10) = {eps1:+.10e}")
    print(f"  eps2 = a2 - (pi-4/10) = {eps2:+.10e}")
    print(f"  eps1 + eps2 = {eps1+eps2:+.4e}   (suma reszt; ~=0 z sumy)")
    print(f"  eps2 - eps1 = {eps2-eps1:+.10e}  (= D - 3/10 = {D-0.3:.10f})")
    print(f"  (eps2-eps1)/2 = eps_asym = {(eps2-eps1)/2:+.8e}")
    print()
    print(f"  |eps1| = |eps2| ~ {(abs(eps1)+abs(eps2))/2:.8e}")
    print()
    eps_avg = (abs(eps1)+abs(eps2))/2

    cands_eps = [
        ('1/6172',       1/6172),
        ('1/6000',       1/6000),
        ('1/(2pi*1000)', 1/(2*np.pi*1000)),
        ('1/(2000*pi)',  1/(2000*np.pi)),
        ('1/pi^4',       1/np.pi**4),
        ('1/(4*pi^3)',   1/(4*np.pi**3)),
        ('1/(8*pi^2)',   1/(8*np.pi**2)),
        ('pi/20000',     np.pi/20000),
        ('phi/10000',    PHI/10000),
        ('1/(e*pi^2)',   1/(np.e*np.pi**2)),
        ('ln2/(4pi^2)',  np.log(2)/(4*np.pi**2)),
        ('1/(2pi^3)',    1/(2*np.pi**3)),
        ('D-3/10',       abs(D-0.3)),    # eps2-eps1 = D-3/10
    ]
    print(f"  Kandydaci na |eps1| ~ {eps_avg:.6e}:")
    for name, val in sorted(cands_eps, key=lambda c: abs(c[1]-eps_avg))[:8]:
        print(f"    {name:<28} = {val:.8e}  diff={eps_avg-val:+.2e}")
    print()

    # ================================================================
    # STOSUNEK
    # ================================================================
    print("=" * 72)
    print("STOSUNEK  a2/a1")
    print("=" * 72)
    ratio = a2/a1
    print(f"  a2/a1 = {ratio:.10f}")
    print()
    cands_ratio = [
        ('(pi-4/10)/(pi-7/10)',  (np.pi-0.4)/(np.pi-0.7)),
        ('9/8',                   9/8),
        ('1+D/a1',                1+D/a1),
        ('1+3/10/a1',             1+0.3/a1),
        ('phi^(1/3)',              PHI**(1/3)),
        ('1+1/pi^2',              1+1/np.pi**2),
        ('sqrt(1+1/5)',            np.sqrt(1.2)),
        ('1+1/phi^4',             1+1/PHI**4),
        ('e/pi+1/2',              np.e/np.pi+0.5),
        ('1+1/(2*pi+1)',          1+1/(2*np.pi+1)),
    ]
    for name, val in sorted(cands_ratio, key=lambda c: abs(c[1]-ratio))[:8]:
        print(f"    {name:<35} = {val:.10f}  diff={ratio-val:+.2e}  ({abs(ratio-val)/ratio*1e6:.1f} ppm)")
    print()

    # ================================================================
    # ROWNANIE KWADRATOWE
    # ================================================================
    print("=" * 72)
    print("ROWNANIE KWADRATOWE  x^2 - S*x + P = 0")
    print("=" * 72)
    print(f"  S = {S:.8f},   P = {P:.8f},   D = {D:.8f}")
    print()
    print(f"  Przyblizone (S=2pi-11/10, D=3/10):")
    print(f"    x1 ~ pi - 7/10 = {np.pi-0.7:.8f}  (eps={eps1:+.4e})")
    print(f"    x2 ~ pi - 4/10 = {np.pi-0.4:.8f}  (eps={eps2:+.4e})")
    print()
    print(f"  Precyzja:")
    print(f"    S - (2pi-11/10) = {S-(2*np.pi-1.1):+.4e}  ({abs(S-(2*np.pi-1.1))/S*1e6:.2f} ppm)")
    print(f"    P - P_hyp       = {P-P_hyp:+.4e}  ({abs(P-P_hyp)/P*1e6:.1f} ppm)")
    print(f"    D - 3/10        = {D-0.3:+.4e}  ({abs(D-0.3)/D*1e6:.1f} ppm)")
    print()
    print("  Interpretacja eps1 ~ -eps2 ~ -1.6e-4:")
    eps_formula_best = 1/(2*np.pi*1000)
    print(f"    Hipoteza eps ~ 1/(2pi*1000)  = {eps_formula_best:.8e}")
    print(f"    Aktualne  eps1 = {eps1:.8e}   diff={eps1+eps_formula_best:.2e}")
    print()

    # ================================================================
    # ALPHA_MIN i F_MIN
    # ================================================================
    print("=" * 72)
    print(f"ALPHA_MIN = {a_min}  (z ex79)   F_MIN = {f_min}  (z ex79)")
    print("=" * 72)
    L = a_min - a1
    R = a2 - a_min
    print(f"  L = a_min - a1 = {L:.8f}")
    print(f"  R = a2 - a_min = {R:.8f}")
    print(f"  L/R            = {L/R:.8f}")
    print(f"  L/D            = {L/D:.8f}")
    print(f"  R/D            = {R/D:.8f}")
    print()

    print(f"  Kandydaci na alpha_min:")
    centrum = (a1+a2)/2
    cands_amin = [
        ('(a1+a2)/2',              centrum),
        ('a1+D/3',                 a1+D/3),
        ('a1+D/phi',               a1+D/PHI),
        ('a1+2*D/5',               a1+2*D/5),
        ('a1+D*phi/(phi+1)',        a1+D*PHI/(PHI+1)),
        ('a1+D*(1-1/phi)/2',        a1+D*(1-1/PHI)/2),
        ('a1+D/e',                  a1+D/np.e),
        ('a1+D*(phi-1)',             a1+D*(PHI-1)),
        ('(2*a1+a2)/3',             (2*a1+a2)/3),
        ('(a1+2*a2)/3',             (a1+2*a2)/3),
        ('sqrt(a1*a2)',              np.sqrt(a1*a2)),
        ('pi-0.578',                np.pi-0.578),
        ('pi-1/sqrt(3)',             np.pi-1/np.sqrt(3)),
        ('pi-ln(1.85)',             np.pi-np.log(1.85)),
        ('a1+D*2/7',                a1+D*2/7),
        ('a1+D*3/8',                a1+D*3/8),
        ('a1+D*4/10',               a1+D*4/10),
        ('a1+D*5/13',               a1+D*5/13),
        ('a1+D*phi^(-2)',            a1+D*PHI**(-2)),
    ]
    for name, val in sorted(cands_amin, key=lambda c: abs(c[1]-a_min))[:10]:
        print(f"    {name:<30} = {val:.8f}  diff={a_min-val:+.4e}  ({abs(a_min-val)*1e6:.0f} ppm)")
    print()

    print(f"  Kandydaci na f_min = {f_min:.4f}:")
    cands_fmin = [
        ('r21*3/pi',              R21_EXP*3/np.pi),
        ('r21*(1-1/23)',          R21_EXP*(1-1/23)),
        ('r21*(1-1/phi^5)',       R21_EXP*(1-1/PHI**5)),
        ('r21*e/pi^2',            R21_EXP*np.e/np.pi**2),
        ('r21*(2-phi)',           R21_EXP*(2-PHI)),
        ('r21*(1-D/S)',           R21_EXP*(1-D/S)),
        ('r21*(P/S^2)',           R21_EXP*(P/S**2)),
        ('r21*cos^2(0.295)',      R21_EXP*np.cos(0.295)**2),
        ('r21*(1-L/S)',           R21_EXP*(1-L/S)),
        ('2*pi^3',                2*np.pi**3),
        ('6*pi^2',                6*np.pi**2),
        ('200',                    200.0),
        ('r21-9',                  R21_EXP-9),
        ('r21*(S-1)/S',           R21_EXP*(S-1)/S),
        ('r21*(pi-D)/(pi)',        R21_EXP*(np.pi-D)/np.pi),
        ('r21*(1-L_D)',            R21_EXP*(1-L/D)),
        ('r21/phi^(1/5)',          R21_EXP/PHI**(1/5)),
        ('r21*(4pi-1)/(4pi)',      R21_EXP*(4*np.pi-1)/(4*np.pi)),
    ]
    for name, val in sorted(cands_fmin, key=lambda c: abs(c[1]-f_min))[:10]:
        ppm = abs(val-f_min)/f_min*1e6
        print(f"    {name:<35} = {val:.4f}  diff={f_min-val:+.4f}  ({ppm:.0f} ppm)")
    print()

    # ================================================================
    # DODATKOWE WZORY ALGEBRAICZNE
    # ================================================================
    print("=" * 72)
    print("DODATKOWE ZWIAZKI ALGEBRAICZNE")
    print("=" * 72)
    print(f"  S = 2pi - 11/10 = 2pi - alpha_TGP/2 - beta_TGP/10")
    print(f"  D ~ 3/10  (1070 ppm od dokladnego)")
    print()
    print(f"  Jezeli alpha*_i = pi - n_i/10 (dokladnie), to n_1+n_2 = 11 (z sumy)")
    print(f"  Najlepsza para: (n1,n2) = (7,4),  tj. alpha*_1 = pi-0.7, alpha*_2 = pi-0.4")
    print(f"  Ale eps1 = {eps1:.6e}, eps2 = {eps2:.6e} (korekty niezerowe!)")
    print()

    # Macierz wzorow liniowych: m*a1 + n*a2 = nice?
    print("  Liniowe kombinacje m*a1 + n*a2 (m,n male calkowite):")
    results_combo = []
    for m in range(-3, 4):
        for n in range(-3, 4):
            if m==0 and n==0: continue
            val = m*a1 + n*a2
            # Sprawdz bliskosc do 2pi, pi, e, phi, sqrt, itp.
            targets = {
                '0': 0, 'pi': np.pi, '2pi': 2*np.pi, 'e': np.e,
                'phi': PHI, '3pi': 3*np.pi, '4pi': 4*np.pi,
                'pi+e': np.pi+np.e, '2e': 2*np.e, 'pi*e': np.pi*np.e,
                '5': 5., '6': 6., '7': 7., '8': 8., '10': 10.,
                '2pi-1': 2*np.pi-1, '2pi+1': 2*np.pi+1,
                'pi^2/e': np.pi**2/np.e, 'e^2': np.e**2,
            }
            for tname, tval in targets.items():
                diff = abs(val - tval)
                if diff < 0.002 and not (m==1 and n==1 and tname=='2pi-1'):
                    results_combo.append((m, n, val, tname, tval, diff))

    results_combo.sort(key=lambda x: x[5])
    seen = set()
    for m, n, val, tname, tval, diff in results_combo[:20]:
        key = (m, n)
        if key in seen: continue
        seen.add(key)
        sign1 = '+' if n >= 0 else ''
        print(f"    {m}*a1{sign1}{n}*a2 = {val:.8f}  ~ {tname} = {tval:.8f}  diff={diff:.2e}")
    print()

    # ================================================================
    # PODSUMOWANIE PRECYZJI WZORU SUMY
    # ================================================================
    print("=" * 72)
    print("HIERARCHIA PRECYZJI WZOROW")
    print("=" * 72)
    formulas = [
        ("S = 2pi-11/10",              abs(S-(2*np.pi-1.1))/S*1e6,      "1 ppm"),
        ("a1 ~ pi-7/10",               abs(eps1)/a1*1e6,                "65 ppm"),
        ("a2 ~ pi-4/10",               abs(eps2)/a2*1e6,                "59 ppm"),
        ("P ~ (pi-7/10)(pi-4/10)",     abs(P-P_hyp)/P*1e6,              "? ppm"),
        ("D ~ 3/10",                   abs(D-0.3)/D*1e6,                "1070 ppm"),
    ]
    for name, ppm, comment in formulas:
        print(f"    {name:<35} {ppm:8.2f} ppm")
    print()

    # TESTY
    print("=" * 72)
    print("TESTY")
    print("=" * 72)
    T1 = abs(S-(2*np.pi-1.1))/S < 5e-6
    print(f"  T1: S ~ 2pi-11/10 (<5 ppm):             {'PASS' if T1 else 'FAIL'}  ({abs(S-(2*np.pi-1.1))/S*1e6:.2f} ppm)")

    T2 = abs(P_hyp-P)/P < 50e-6
    print(f"  T2: P ~ (pi-7/10)(pi-4/10) (<50 ppm):  {'PASS' if T2 else 'FAIL'}  ({abs(P_hyp-P)/P*1e6:.1f} ppm)")

    T3 = abs(D-0.3)/D < 2e-3
    print(f"  T3: D ~ 3/10 (<2000 ppm):               {'PASS' if T3 else 'FAIL'}  ({abs(D-0.3)/D*1e6:.1f} ppm)")

    T4 = abs(eps1+eps2) < 1e-5
    print(f"  T4: eps1+eps2 ~ 0 (<1e-5):              {'PASS' if T4 else 'FAIL'}  ({abs(eps1+eps2):.3e})")

    T5 = L/R < 0.75  # dolina asymetryczna
    print(f"  T5: Dolina asymetryczna L/R < 0.75:     {'PASS' if T5 else 'FAIL'}  (L/R={L/R:.4f})")

    n_pass = sum([T1, T2, T3, T4, T5])
    print(f"\nWYNIK: {n_pass}/5")

    print()
    print("=" * 72)
    print("WNIOSEK (EX80)")
    print("=" * 72)
    print(f"  alpha*_1 = {a1:.8f}  ~  pi - 7/10 = {np.pi-0.7:.8f}  (eps={eps1:.3e})")
    print(f"  alpha*_2 = {a2:.8f}  ~  pi - 4/10 = {np.pi-0.4:.8f}  (eps={eps2:.3e})")
    print()
    print(f"  [1] SUMA S = {S:.8f}")
    print(f"      2pi-11/10 = {2*np.pi-1.1:.8f}  ROZNICA = {S-(2*np.pi-1.1):.3e}  ({abs(S-(2*np.pi-1.1))/S*1e6:.2f} ppm)")
    print()
    print(f"  [2] ILOCZYN P = {P:.8f}")
    print(f"      (pi-7/10)(pi-4/10) = {P_hyp:.8f}  ROZNICA = {P-P_hyp:.3e}  ({abs(P-P_hyp)/P*1e6:.1f} ppm)")
    print()
    print(f"  [3] ROZNICA D = {D:.8f}")
    print(f"      3/10 = 0.30000000  ROZNICA = {D-0.3:.3e}  ({abs(D-0.3)/D*1e6:.1f} ppm)")
    print()
    print(f"  Korekta eps = alpha*_1 - (pi-7/10) = {eps1:.6e}")
    print(f"  Wzor: alpha*_i = pi - n_i/10 + eps_i,  eps_1 ~ -eps_2 ~ -1.6e-4")
    print(f"  Suma eliminuje korekte: S = 2pi - 11/10 dokladnie (do precyzji numerycznej)")
    print()
    print(f"  alpha_min = {a_min:.5f},  f_min = {f_min:.4f}")
    print(f"  Najblizsze formuly:")
    print(f"    f_min ~ r21*3/pi = {R21_EXP*3/np.pi:.4f}  (diff={f_min-R21_EXP*3/np.pi:+.4f})")
    bfc = sorted(cands_amin, key=lambda c: abs(c[1]-a_min))[0]
    print(f"    alpha_min ~ {bfc[0]} = {bfc[1]:.6f}  (diff={a_min-bfc[1]:.4e})")
    print("=" * 72)
