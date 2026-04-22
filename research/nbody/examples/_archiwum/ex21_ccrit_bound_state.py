"""
ex21_ccrit_bound_state.py - Krytyczne C dla trojcialowego stanu zwiazanego

Pytanie: dla jakiego C trojka obiektow ma E_total(2B+3B) < 0
         podczas gdy E_total(2B) > 0 ?

To jest predykcja falsyfikowalna TGP: istnieja trojkaty zwiazane
przez CZYSTO 3-cialowa interakcje.

Metoda:
  Twierdzenie wirialne dla mieszanego potencjalu V2+V3:
    2T = -(r dV2/dr) - (r dV3/dr)
    E_total = T + V2 + V3

  Dla potencjalu Yukawa:
    r dV2/dr = -C^2(1+mr)exp(-mr)/r = -(1+mr)V2_ij/pair
    Srednia: <r dV2/dr> = -(1+md_avg) * V2_total

  Dla calki Feynmana (skalowanie):
    V3 ~ d^alpha * I_Y(d): r dV3/dr ≈ -(alpha + m_eff*d)*V3
    Gdzie alpha ~ -1.5 z calki Feynmana (z asymptoty K_0)

  Warunek stanu zwiazanego: E_total < 0
    T + V2 + V3 < 0
    T = -(1+md)/2 * V2_total - (alpha+md)/2 * V3_total
    => (1-(1+md)/2)*V2 + (1-(alpha+md)/2)*V3 < 0

  Wariacja: szukamy C takie ze E_virial = 0:
    E_virial(C) = a(d)*C^2 + b(d)*C^3 = 0  (wyznacza C_crit)
"""

import numpy as np
from scipy.special import k0 as K0, k1 as K1
from scipy.optimize import brentq
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# ─── Feynman integral (wektoryzowany) ─────────────────────────────────────────
_pts, _wts = np.polynomial.legendre.leggauss(25)
_up = 0.5*(1+_pts); _uw = 0.5*_wts
_vp = 0.5*(1+_pts); _vw = 0.5*_wts
_UU, _VV = np.meshgrid(_up, _vp, indexing='ij')
_WW = np.outer(_uw, _vw)
_A1 = _UU; _A2 = _VV*(1-_UU); _A3 = (1-_UU)*(1-_VV); _JAC = (1-_UU)

M_SP = 1.0
C_PL = 0.282094

def I_Y(d12, d13, d23, m=M_SP):
    D = _A1*_A2 + _A1*_A3 + _A2*_A3
    Q = _A2*d12**2 + _A1*d13**2 + _A3*d23**2
    good = (D > 1e-30) & (Q > 1e-30)
    u = np.where(good, m*np.sqrt(Q/D), 1.0)
    val = np.where(good, D**(-1.5)*K0(u), 0.0)
    return 2.0*np.sum(_WW*_JAC*val)

def V2_equilat(d, C, m=M_SP):
    """V2 trojkata rownobocznego (3 pary)."""
    return -3.0 * C**2 * np.exp(-m*d)/d

def V3_equilat(d, C, m=M_SP):
    """V3 trojkata rownobocznego."""
    return -C**3 * I_Y(d, d, d, m)

def rdVdr_2b(d, C, m=M_SP):
    """<r dV2/dr> = sum_pairs C^2*(mr+1)*exp(-mr)/r * (-1) = -(1+md)*V2/V2_coeff."""
    # V2 = -3C^2*exp(-md)/d => dV2/dr = 3C^2*(m+1/d)*exp(-md)/d
    # r*dV2/dr = 3C^2*(md+1)*exp(-md)/d = -(1+md)*V2_total
    return -(1 + m*d) * V2_equilat(d, C, m)  # = +(1+md)*3C^2*exp(-md)/d > 0

def rdVdr_3b_approx(d, C, m=M_SP):
    """
    Przyblizenie <r dV3/dr> dla trojkata rownobocznego.
    Z calki Feynmana: V3 ~ -C^3 * f(md) gdzie f(t) = I_Y(t,t,t;1)/t^3 * t^3 ...
    Numeryczne: oblicz d*dV3/dd przez rozniczkowanie skonczone.
    """
    eps = d * 1e-4
    V3p = V3_equilat(d+eps, C, m)
    V3m = V3_equilat(d-eps, C, m)
    dV3_dd = (V3p - V3m)/(2*eps)
    return d * dV3_dd

def E_virial(d, C, m=M_SP):
    """
    Energia calkowita przez twierdzenie wirialne:
    2T = -(rdV2dr + rdV3dr)
    E = T + V2 + V3 = (V2+V3) - (rdV2dr+rdV3dr)/2
    """
    v2 = V2_equilat(d, C, m)
    v3 = V3_equilat(d, C, m)
    rV2 = rdVdr_2b(d, C, m)
    rV3 = rdVdr_3b_approx(d, C, m)
    return (v2 + v3) - (rV2 + rV3)/2.0

def E_virial_2b_only(d, C, m=M_SP):
    """Energia wirialna bez V3 (tylko 2B)."""
    v2 = V2_equilat(d, C, m)
    rV2 = rdVdr_2b(d, C, m)
    return v2 - rV2/2.0

# ════════════════════════════════════════════════════════════════
# SEKCJA A: E_virial vs d dla roznych C
# ════════════════════════════════════════════════════════════════

print("="*65)
print("SEKCJA A: Energia wirialna vs d dla roznych C")
print("="*65)
print("Trojkat rownoboczny, m_sp=1.0")
print()
print("E_virial < 0 = stan zwiazany, > 0 = niewiazany")
print()

C_test_vals = [0.05, 0.10, 0.15, 0.20, 0.25, 0.282]

for C in C_test_vals:
    d_scan = np.linspace(0.3, 4.0, 100)
    E_2b_arr  = np.array([E_virial_2b_only(d, C) for d in d_scan])
    E_tot_arr = np.array([E_virial(d, C)         for d in d_scan])

    d_min_2b  = d_scan[np.argmin(E_2b_arr)]
    E_min_2b  = np.min(E_2b_arr)
    d_min_tot = d_scan[np.argmin(E_tot_arr)]
    E_min_tot = np.min(E_tot_arr)

    bound_2b  = "ZWIAZ" if E_min_2b  < 0 else "luzy "
    bound_tot = "ZWIAZ" if E_min_tot < 0 else "luzy "

    print(f"C={C:.3f}:  E_min(2B)={E_min_2b:+.4f} [{bound_2b}] @ d={d_min_2b:.2f}   "
          f"E_min(2B+3B)={E_min_tot:+.4f} [{bound_tot}] @ d={d_min_tot:.2f}")

# ════════════════════════════════════════════════════════════════
# SEKCJA B: Szukanie C_crit
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("SEKCJA B: Szukanie C_crit - przejscie 2B+3B z luzy w zwiazany")
print("="*65)
print()

def E_min_total(C):
    """Minimalna energia wirialna (2B+3B) po d."""
    d_scan = np.linspace(0.3, 5.0, 80)
    E_arr = np.array([E_virial(d, C) for d in d_scan])
    return np.min(E_arr)

def E_min_2b(C):
    """Minimalna energia wirialna (tylko 2B) po d."""
    d_scan = np.linspace(0.3, 5.0, 80)
    E_arr = np.array([E_virial_2b_only(d, C) for d in d_scan])
    return np.min(E_arr)

print("Skan E_min vs C:")
print(f"{'C':>8}  {'E_min(2B)':>12}  {'E_min(2B+3B)':>14}  {'Delta E':>10}  {'Status':>20}")
print("-"*75)

C_scan = np.array([0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20, 0.22, 0.25, 0.282])
E2b_arr, E3b_arr = [], []

for C in C_scan:
    e2 = E_min_2b(C)
    e3 = E_min_total(C)
    E2b_arr.append(e2)
    E3b_arr.append(e3)
    status = ""
    if e2 > 0 and e3 < 0:
        status = "*** 3B-INDUCED BOUND! ***"
    elif e2 > 0:
        status = "obie luze"
    else:
        status = "obie zwiazane"
    print(f"{C:8.3f}  {e2:12.5f}  {e3:14.5f}  {e3-e2:10.5f}  {status}")

E2b_arr = np.array(E2b_arr)
E3b_arr = np.array(E3b_arr)

# Szukanie C_crit(2B): gdzie E_min_2b = 0
print()
try:
    # Sprawdz czy jest znak zmiana
    if np.any(E2b_arr < 0) and np.any(E2b_arr > 0):
        C_crit_2b = brentq(E_min_2b, C_scan[E2b_arr > 0][-1]*0.9,
                           C_scan[E2b_arr < 0][0]*1.1, xtol=1e-4)
        print(f"C_crit(2B) = {C_crit_2b:.5f}  (powyzej: zwiazany 2B)")
    else:
        C_crit_2b = None
        print(f"E_min(2B) zawsze {'<0' if E2b_arr[-1]<0 else '>0'} w zakresie C={C_scan[0]:.3f}-{C_scan[-1]:.3f}")
except Exception as ex:
    print(f"Szukanie C_crit(2B): {ex}")

# Szukanie C_crit(3B): gdzie E_min(2B+3B) = 0
try:
    if np.any(E3b_arr < 0) and np.any(E3b_arr > 0):
        C_crit_3b = brentq(E_min_total, C_scan[E3b_arr > 0][-1]*0.9,
                           C_scan[E3b_arr < 0][0]*1.1, xtol=1e-4)
        print(f"C_crit(3B) = {C_crit_3b:.5f}  (powyzej: zwiazany 2B+3B)")
    else:
        print(f"E_min(2B+3B) zawsze {'<0' if E3b_arr[-1]<0 else '>0'} w zakresie badanym")
        # sprawdz dla bardzo malych C
        E_small = E_min_total(0.001)
        print(f"  E_min(2B+3B) @ C=0.001: {E_small:.5f}")
except Exception as ex:
    print(f"Szukanie C_crit(3B): {ex}")

# ════════════════════════════════════════════════════════════════
# SEKCJA C: Rownianie na C_crit — analitycznie
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("SEKCJA C: Analityczne rownianie dla C_crit")
print("="*65)
print()

# Energia wirialna:
# E = V2 + V3 - (rdV2dr + rdV3dr)/2
# = -A*C^2 + B(d)*C^2 - D*C^3 + E(d)*C^3
# gdzie A,B,D,E sa funkcjami d

# Przy minimum po d (dd(E)/dd=0):
# E = a2(d)*C^2 + a3(d)*C^3 = 0  => C* = -a2(d)/a3(d)

print("Koeficjenty energii wirialnej dla trojkata rownobocznego:")
print("E_virial = a2(d)*C^2 + a3(d)*C^3")
print()

def a2_coeff(d, m=M_SP):
    """Koeficjent przy C^2: (V2 - rdV2dr/2)/C^2"""
    V2_unit = -3.0*np.exp(-m*d)/d
    rV2_unit = -(1+m*d)*V2_unit  # = +3(1+md)*exp(-md)/d
    return V2_unit - rV2_unit/2.0

def a3_coeff(d, m=M_SP):
    """Koeficjent przy C^3: (V3 - rdV3dr/2)/C^3"""
    C_unit = 1.0
    eps = d*1e-4
    V3u = -I_Y(d,d,d,m)          # V3 = C^3 * V3u
    V3p = -I_Y(d+eps,d+eps,d+eps,m)
    V3m = -I_Y(d-eps,d-eps,d-eps,m)
    dV3_dd = (V3p-V3m)/(2*eps)
    rV3 = d*dV3_dd
    return V3u - rV3/2.0

print(f"{'d':>6}  {'a2(d)':>12}  {'a3(d)':>12}  {'C*=-a2/a3':>12}  {'C_Pl?':>8}")
print("-"*60)

C_crit_vals = []
for d in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0]:
    a2 = a2_coeff(d)
    a3 = a3_coeff(d)
    if a3 != 0 and -a2/a3 > 0:
        C_star = -a2/a3
        below = "YES" if C_star < C_PL else ""
    else:
        C_star = float('nan')
        below = ""
    C_crit_vals.append((d, C_star))
    print(f"{d:6.2f}  {a2:12.5f}  {a3:12.5f}  {C_star:12.5f}  {below:>8}")

# Minimum C_crit
valid = [(d,c) for d,c in C_crit_vals if not np.isnan(c) and c > 0]
if valid:
    d_opt, C_opt = min(valid, key=lambda x: x[1])
    print()
    print(f"Minimalne C_crit = {C_opt:.5f} przy d = {d_opt:.2f} l_Pl")
    print(f"C_Planck = {C_PL:.5f}")
    if C_opt < C_PL:
        print("=> ISTNIEJE okno C_crit < C < C_Pl gdzie tylko 3B tworzy stan zwiazany!")
    else:
        print("=> C_crit > C_Pl: dla obiektow Plancka stan zwiazany istnieje juz przez 2B")

# ════════════════════════════════════════════════════════════════
# SEKCJA D: Konkretne wartosci fizyczne stanu zwiazanego
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("SEKCJA D: Energia i promien stanu zwiazanego N=3, Planck")
print("="*65)
print()

# Minimalizuj E_virial(d) dla C=C_PL
from scipy.optimize import minimize_scalar
res = minimize_scalar(lambda d: E_virial(d, C_PL),
                      bounds=(0.3, 5.0), method='bounded')
d_eq = res.x
E_eq = res.fun

v2_eq  = V2_equilat(d_eq, C_PL)
v3_eq  = V3_equilat(d_eq, C_PL)
rV2_eq = rdVdr_2b(d_eq, C_PL)
rV3_eq = rdVdr_3b_approx(d_eq, C_PL)
T_eq   = -(rV2_eq + rV3_eq)/2.0

print(f"Rownowaga wirialna dla trojkata C={C_PL:.4f}, m_sp={M_SP}:")
print(f"  d_eq    = {d_eq:.4f} l_Pl")
print(f"  V2      = {v2_eq:.5f} E_Pl")
print(f"  V3      = {v3_eq:.5f} E_Pl  ({100*v3_eq/v2_eq:.1f}% V2)")
print(f"  <r dV2/dr> = {rV2_eq:.5f}")
print(f"  <r dV3/dr> = {rV3_eq:.5f}")
print(f"  T (wirial)  = {T_eq:.5f} E_Pl")
print(f"  E_total     = {E_eq:.5f} E_Pl")
print()

# Fizyczne jednostki
E_Pl_J = 1.956e9   # J (energia Plancka)
l_Pl_m = 1.616e-35 # m
print(f"W jednostkach fizycznych:")
print(f"  d_eq    = {d_eq*l_Pl_m:.3e} m = {d_eq:.2f} l_Pl")
print(f"  E_total = {E_eq*E_Pl_J:.3e} J = {E_eq:.5f} E_Pl")
print(f"  |E_bind|= {abs(E_eq)*E_Pl_J:.3e} J")

# Stan 2B sam w sobie:
E_eq_2b = minimize_scalar(lambda d: E_virial_2b_only(d, C_PL),
                           bounds=(0.3,5.0), method='bounded').fun
print()
print(f"Dla porownania — stan 2B (bez 3-cialowej):")
print(f"  E_min(2B) = {E_eq_2b:.5f} E_Pl")
print(f"  Wzmocnienie przez V3: {(E_eq-E_eq_2b)/abs(E_eq_2b)*100:.1f}%")

# ════════════════════════════════════════════════════════════════
# SEKCJA E: Mapa (C, m_sp) — granica stanu zwiazanego
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("SEKCJA E: Granica zwiazana w przestrzeni (C, d_eq)")
print("="*65)
print()
print("Wyznaczanie d_rownowagi i E_total vs C:")
print()
print(f"{'C':>8}  {'d_eq':>8}  {'E_tot':>10}  {'E_2B':>10}  {'V3/V2':>8}  {'Stan':>15}")
print("-"*65)

for C in np.array([0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20, 0.22, 0.25, 0.282]):
    res2 = minimize_scalar(lambda d: E_virial(d, C), bounds=(0.3,5.0), method='bounded')
    res3 = minimize_scalar(lambda d: E_virial_2b_only(d, C), bounds=(0.3,5.0), method='bounded')
    d_e = res2.x; E_e = res2.fun; E_2b_e = res3.fun
    v2e = V2_equilat(d_e, C); v3e = V3_equilat(d_e, C)
    frac = v3e/v2e if v2e != 0 else 0

    if E_2b_e > 0 and E_e < 0:
        status = "3B-INDUCED!"
    elif E_e < 0:
        status = "zwiazany (2B+3B)"
    else:
        status = "luzy"

    print(f"{C:8.3f}  {d_e:8.3f}  {E_e:10.5f}  {E_2b_e:10.5f}  {frac:8.4f}  {status:>15}")

# ════════════════════════════════════════════════════════════════
# PODSUMOWANIE
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("PODSUMOWANIE ex21")
print("="*65)
print(f"""
Metoda: twierdzenie wirialne dla potencjalu V2+V3

1. ENERGIA WIRIALNA:
   E_virial = (V2+V3) - (r*dV2/dr + r*dV3/dr)/2
   Uwzglednia wklad kinetyczny — bardziej realistyczne niz samo V_pot.

2. STAN 2B-ONLY:
   Dla trojkata rownobocznego z C=C_Pl:
   E_min(2B) = ? (zbadaj wyniki Sekcji E)

3. STAN 2B+3B:
   E_min(2B+3B) < E_min(2B) — sila 3-cialowa zawsze glębiej zwiazuje.

4. OKNO 3B-INDUCED:
   Jezeli istnieje C takie ze E_min(2B)>0 ale E_min(2B+3B)<0 =>
   czysto 3-cialowy stan zwiazany — unikalna predykcja TGP.

5. OPTYMALNE C i d:
   d_eq (Planck) = wynik Sekcji D
   E_total = wynik Sekcji D
""")
