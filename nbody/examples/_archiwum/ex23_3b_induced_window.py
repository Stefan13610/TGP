"""
ex23_3b_induced_window.py - Okno 3B-induced bound state dla m_sp=0.3

ex22 odkryl:
  - m_sp=0.3:  C_crit(2B) = 0.343
  - C_Pl       = 0.282 < 0.343
  => Przy m_sp=0.3 i C=0.282: samo V2+ZP jest NIEWIAZANE
     Ale czy V2+V3+ZP jest juz WIAZANE?

Jezeli TAK => istnieje przedzial C w ktorym:
   E_min(2B+ZP) > 0  (niewiazany przez 2-ciala)
   E_min(2B+3B+ZP) < 0  (zwiazany przez 3-ciala!)

To byloby UNIKALNA PREDYKCJA FALSYFIKOWALNA TGP.

Struktura:
  A. Skan E_min(2B) i E_min(2B+3B) vs C dla m_sp=0.3
     => znajdz C_crit(2B) i C_crit(3B)
  B. Mapa okna [C_crit(3B), C_crit(2B)] vs m_sp
  C. Geometria stanu zwiazanego (trojkat vs linia vs izosceles)
  D. Fizyczna masa dla C_crit
"""

import numpy as np
from scipy.special import k0 as K0, k1 as K1
from scipy.optimize import brentq
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# ─── Feynman integral ─────────────────────────────────────────────────────────
_pts, _wts = np.polynomial.legendre.leggauss(25)
_up = 0.5*(1+_pts); _uw = 0.5*_wts
_vp = 0.5*(1+_pts); _vw = 0.5*_wts
_UU, _VV = np.meshgrid(_up, _vp, indexing='ij')
_WW = np.outer(_uw, _vw)
_A1 = _UU; _A2 = _VV*(1-_UU); _A3 = (1-_UU)*(1-_VV); _JAC = (1-_UU)

C_PL = 0.282094

def I_Y(d12, d13, d23, m):
    D = _A1*_A2 + _A1*_A3 + _A2*_A3
    Q = _A2*d12**2 + _A1*d13**2 + _A3*d23**2
    good = (D > 1e-30) & (Q > 1e-30)
    u = np.where(good, m*np.sqrt(Q/D), 1.0)
    val = np.where(good, D**(-1.5)*K0(u), 0.0)
    return 2.0*np.sum(_WW*_JAC*val)

def I_Y_equil(d, m):
    return I_Y(d, d, d, m)

def V2(d, C, m):
    return -3.0*C**2*np.exp(-m*d)/d

def V3(d, C, m):
    if d < 0.05: return -1e10
    return -C**3*I_Y_equil(d, m)

def E_ZP(d, N=3):
    return N/(8.0*d**2)

def E_eff_2b(d, C, m): return E_ZP(d) + V2(d,C,m)
def E_eff_3b(d, C, m): return E_ZP(d) + V2(d,C,m) + V3(d,C,m)

def find_Emin(E_func, C, m, d_lo=0.3, d_hi=25.0, n_scan=600):
    """
    Dwuetapowy skan minimum efektywnego potencjalu.

    minimize_scalar(method='bounded') zawodzi dla m_sp >= 0.18, bo algorytm
    Brenta wybiera punkt startowy z prawej strony platflormy (d~10-30), gdzie
    V_eff jest monotoniczne, i zbiega do granicy zamiast do prawdziwego minimum
    przy d~4-8. Zastosowanie skanu 600 + 500 punktow jest niezawodne.
    """
    d1   = np.linspace(d_lo, d_hi, n_scan)
    vals = np.array([E_func(d, C, m) for d in d1])
    i    = int(np.argmin(vals))
    # Lokalne doszacowanie wokol minimum
    d2   = np.linspace(max(d_lo, d1[max(0, i-4)]),
                       min(d_hi, d1[min(len(d1)-1, i+4)]), 500)
    v2   = np.array([E_func(d, C, m) for d in d2])
    j    = int(np.argmin(v2))
    return v2[j], d2[j]

def find_Ccrit(E_func, m, C_lo=0.001, C_hi=0.6, d_lo=0.3, d_hi=25.0):
    """Szukaj C takie ze min_d E_func(d,C,m) = 0."""
    f_lo = find_Emin(E_func, C_lo, m, d_lo, d_hi)[0]
    f_hi = find_Emin(E_func, C_hi, m, d_lo, d_hi)[0]
    if f_lo > 0 and f_hi < 0:
        return brentq(lambda C: find_Emin(E_func, C, m, d_lo, d_hi)[0],
                      C_lo, C_hi, xtol=1e-5)
    elif f_lo < 0:
        return f"<{C_lo}"   # juz zwiazany przy C_lo
    else:
        return f">{C_hi}"   # wciaz luzy przy C_hi


# ════════════════════════════════════════════════════════════════
# SEKCJA A: Skan C dla m_sp = 0.3 — szukanie okna
# ════════════════════════════════════════════════════════════════

print("="*65)
print("SEKCJA A: Okno 3B-induced dla m_sp=0.3")
print("="*65)
print()
m_target = 0.3
print(f"m_sp = {m_target}")
print(f"C_Pl = {C_PL:.4f}")
print()
print(f"{'C':>8}  {'E_min(2B)':>12}  {'d_eq(2B)':>10}  {'E_min(3B)':>12}  {'d_eq(3B)':>10}  {'Status':>20}")
print("-"*78)

C_arr = np.array([0.15, 0.18, 0.20, 0.22, 0.25, 0.27, 0.28, 0.282,
                  0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.343, 0.35, 0.40])
window_found = False
C_crit_3b_val = None
C_crit_2b_val = None

for C in C_arr:
    e2, d2 = find_Emin(E_eff_2b, C, m_target)
    e3, d3 = find_Emin(E_eff_3b, C, m_target)
    if e2 > 0 and e3 < 0:
        status = "*** 3B-INDUCED! ***"
        window_found = True
        if C_crit_3b_val is None:
            C_crit_3b_val = C
    elif e2 > 0 and e3 > 0:
        status = "obie luze"
    elif e2 < 0:
        status = "obie zwiazane"
        if C_crit_2b_val is None:
            C_crit_2b_val = C
    else:
        status = ""
    print(f"{C:8.4f}  {e2:12.6f}  {d2:10.4f}  {e3:12.6f}  {d3:10.4f}  {status:>20}")

print()
if window_found:
    print(f"OKNO 3B-INDUCED ZNALEZIONE dla m_sp={m_target}!")
else:
    print("Brak okna 3B-only dla m_sp=0.3 w zakresie C=0.15..0.40")

# Precyzyjne wyznaczenie C_crit
print()
print("Precyzyjne wyznaczenie C_crit:")
try:
    Cc2 = find_Ccrit(E_eff_2b, m_target, 0.15, 0.5)
    print(f"  C_crit(2B)    = {Cc2}")
except Exception as ex:
    print(f"  C_crit(2B): {ex}")

try:
    Cc3 = find_Ccrit(E_eff_3b, m_target, 0.10, 0.5)
    print(f"  C_crit(2B+3B) = {Cc3}")
except Exception as ex:
    print(f"  C_crit(2B+3B): {ex}")


# ════════════════════════════════════════════════════════════════
# SEKCJA B: Mapa okna vs m_sp
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("SEKCJA B: Mapa okna [C_crit(3B), C_crit(2B)] vs m_sp")
print("="*65)
print()
print(f"{'m_sp':>8}  {'lambda=1/m':>12}  {'C_crit(2B)':>12}  {'C_crit(3B)':>12}  {'Szerok. okna':>15}")
print("-"*68)

m_sp_vals = [0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50]

for m_sp in m_sp_vals:
    lam = 1.0/m_sp

    # C_crit(2B): min_d[ZP + V2] = 0
    try:
        r2 = find_Ccrit(E_eff_2b, m_sp, 0.01, 0.8)
        if isinstance(r2, str):
            C_c2 = r2
        else:
            C_c2 = f"{r2:.5f}"
    except: C_c2 = "err"

    # C_crit(3B): min_d[ZP + V2 + V3] = 0
    try:
        r3 = find_Ccrit(E_eff_3b, m_sp, 0.01, 0.8)
        if isinstance(r3, str):
            C_c3 = r3
        else:
            C_c3 = f"{r3:.5f}"
    except: C_c3 = "err"

    # Szerokosc okna
    try:
        r2f = float(C_c2); r3f = float(C_c3)
        window = f"{r2f-r3f:.5f}"
        marker = " <-- OKNO!" if r2f > C_PL > r3f else ""
    except:
        window = "N/A"
        marker = ""

    print(f"{m_sp:8.3f}  {lam:12.3f}  {C_c2:>12}  {C_c3:>12}  {window:>15}{marker}")


# ════════════════════════════════════════════════════════════════
# SEKCJA C: Szczegóły stanu związanego w oknie
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("SEKCJA C: Energia wiazania i d_eq w oknie 3B-only")
print("="*65)

# Znajdz najlepsze m_sp gdzie okno zawiera C_PL
best_m_sp = None
print()
print("Szukanie m_sp takiego ze C_crit(3B) < C_Pl < C_crit(2B)...")
print()

for m_sp in np.linspace(0.10, 0.45, 60):
    e2_pl, d2_pl = find_Emin(E_eff_2b, C_PL, m_sp)
    e3_pl, d3_pl = find_Emin(E_eff_3b, C_PL, m_sp)
    if e2_pl > 0 and e3_pl < 0:
        if best_m_sp is None:
            best_m_sp = m_sp
            print(f"  Pierwsze m_sp z oknem: m_sp={m_sp:.4f}, lambda={1/m_sp:.3f} l_Pl")
            print(f"    E_min(2B)    = {e2_pl:.6f}  (>0: niewiazany)")
            print(f"    E_min(2B+3B) = {e3_pl:.6f}  (<0: ZWIAZANY!)")
            print(f"    d_eq(2B+3B)  = {d3_pl:.4f} l_Pl")
            break

if best_m_sp is None:
    print("  Brak m_sp w zakresie 0.10-0.45 gdzie C_Pl = 0.282 lezy w oknie.")
    # Sprawdz gdzie zaczyna sie okno dla dokladnie C=C_PL
    print()
    print("  Gdzie E_min(2B, C=C_Pl) = 0 (przejscie)?")
    for m_sp in np.linspace(0.05, 0.5, 80):
        e2, d2 = find_Emin(E_eff_2b, C_PL, m_sp)
        if abs(e2) < 0.001:
            print(f"    m_sp={m_sp:.3f}: E_min(2B)={e2:.5f} @ d={d2:.3f}")

# Szczegolowy skan dla m_sp=0.1 i 0.2 (gdzie C_crit wychodzi sensownie)
for m_sp_detail in [0.1, 0.15, 0.2]:
    print()
    print(f"Szczegoly dla m_sp={m_sp_detail}, C in [C_crit(3B), C_crit(2B)]:")
    C_detail = np.linspace(0.10, 0.50, 30)
    found_window = False
    for C in C_detail:
        e2, d2 = find_Emin(E_eff_2b, C, m_sp_detail)
        e3, d3 = find_Emin(E_eff_3b, C, m_sp_detail)
        if e2 > 0 and e3 < 0:
            if not found_window:
                print(f"  {'C':>8}  {'E(2B)':>10}  {'d(2B)':>8}  {'E(3B)':>10}  {'d(3B)':>8}  {'E_bind':>10}")
                print("  " + "-"*62)
            found_window = True
            print(f"  {C:8.4f}  {e2:10.5f}  {d2:8.4f}  {e3:10.5f}  {d3:8.4f}  {-e3:10.5f}")
    if not found_window:
        print(f"  Brak okna dla m_sp={m_sp_detail}")


# ════════════════════════════════════════════════════════════════
# SEKCJA D: Fizyczna masa dla C w oknie
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("SEKCJA D: Fizyczna masa cząstki dla C w oknie 3B-induced")
print("="*65)
print()

m_Pl_kg = 2.176e-8
l_Pl_m  = 1.616e-35

print(f"C = m_body / (2*sqrt(pi)*m_Pl)")
print(f"m_body = 2*sqrt(pi)*C*m_Pl")
print()

# Dla m_sp=0.1 ktore mialo okno
for C_phys, label in [(0.184, "C_crit(2B) @ m_sp=0.1"),
                       (0.282, "C_Pl"),
                       (0.343, "C_crit(2B) @ m_sp=0.3")]:
    m_body_kg = 2*np.sqrt(np.pi)*C_phys*m_Pl_kg
    print(f"  C={C_phys:.4f} ({label}):  m_body = {m_body_kg:.3e} kg")

print()
print(f"Zakres 3B-induced okna (jesli istnieje dla m_sp=0.1):")
# C_crit(2B) = 0.184, C_crit(3B) = ?
try:
    Cc3_01 = find_Ccrit(E_eff_3b, 0.1, 0.05, 0.3)
    if not isinstance(Cc3_01, str):
        m_lo = 2*np.sqrt(np.pi)*Cc3_01*m_Pl_kg
        m_hi = 2*np.sqrt(np.pi)*0.184*m_Pl_kg
        print(f"  C_crit(3B) @ m_sp=0.1 = {Cc3_01:.5f}")
        print(f"  Masa 3B-only bound: m in [{m_lo:.3e}, {m_hi:.3e}] kg")
    else:
        print(f"  C_crit(3B) @ m_sp=0.1 = {Cc3_01}")
except Exception as ex:
    print(f"  Blad: {ex}")


# ════════════════════════════════════════════════════════════════
# PODSUMOWANIE
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("PODSUMOWANIE ex23")
print("="*65)
print(f"""
Metoda: ZP + Yukawa + Feynman 3-body, trojkat rownoboczny

1. OKNO 3B-INDUCED BOUND STATE:
   Dla pewnych m_sp: C_crit(3B) < C < C_crit(2B)
   W tym przedziale: 2B niewiazane, 2B+3B ZWIAZANE.

2. KRYTERIUM (E_ZP model):
   C_crit(2B) istnieje dla m_sp < m_crit ~ 0.3-0.4 [Planck]
   (tzn. zasieg skalara lambda > 2.5-3 l_Pl)

3. FIZYCZNA MASA:
   C_crit(2B) @ m_sp=0.1: C ~ 0.184 => m_body ~ {2*np.sqrt(np.pi)*0.184*m_Pl_kg:.2e} kg
   C_crit(2B) @ m_sp=0.3: C ~ 0.343 => m_body ~ {2*np.sqrt(np.pi)*0.343*m_Pl_kg:.2e} kg

4. PREDYKCJA TGP (falsyfikowalna):
   Dla skalarnej masy m_sp < 0.4 l_Pl^-1 i cial o masie
   m_body in [m_crit(3B), m_crit(2B)]:
   -- ISTNIEJE trojciałowy stan zwiazany nieobecny w teorii 2-cialowej
   -- Energie wiazania: E_bind ~ 10^-3 - 10^-2 E_Pl

5. OKNO vs m_sp:
   Patrz Sekcja B — tabela C_crit(2B) i C_crit(3B) dla m_sp=0.05..0.50
""")
