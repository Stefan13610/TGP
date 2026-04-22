"""
ex22_quantum_bound_state.py - Kwantowe zero-point energy i C_crit^QM

Problem klasyczny (ex21): Yukawa V ~ -1/r kolapuje do r=0 => zawsze zwiazany.
Problem kwantowy: zasada nieoznaczonosci Heisenberga:
  Delta_x * Delta_p >= hbar/2
  Dla trojkata o boku d: Delta_p >= hbar/(2d)
  E_kin^ZP = (Delta_p)^2 / (2m) = hbar^2 / (8 m d^2)  [na czastke]

W jednostkach Plancka hbar=1, m=1 => E_ZP = 1/(8d^2) na czastke.

Efektywna energia calkowita (przyblizenie):
  E_eff(d, C) = 3 * E_ZP(d) + V2(d,C) + V3(d,C)
             = 3/(8d^2) + V2 + V3

Minimum E_eff(d) wyznacza rownowagowe d_eq^QM i energie wiazania.

Warunek 3B-induced bound state:
  min_d [3/(8d^2) + V2(d,C)] > 0   (sam 2B + ZP: NIEWIAZANY)
  min_d [3/(8d^2) + V2(d,C) + V3(d,C)] < 0  (2B+3B + ZP: ZWIAZANY)

To wyznacza OKNO C w ktorym czysto 3-cialowy stan zwiazany istnieje!
"""

import numpy as np
from scipy.special import k0 as K0, k1 as K1
from scipy.optimize import minimize_scalar, brentq
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# ─── Feynman integral (n=25) ───────────────────────────────────────────────────
_pts, _wts = np.polynomial.legendre.leggauss(25)
_up = 0.5*(1+_pts); _uw = 0.5*_wts
_vp = 0.5*(1+_pts); _vw = 0.5*_wts
_UU, _VV = np.meshgrid(_up, _vp, indexing='ij')
_WW = np.outer(_uw, _vw)
_A1 = _UU; _A2 = _VV*(1-_UU); _A3 = (1-_UU)*(1-_VV); _JAC = (1-_UU)

M_SP = 1.0
C_PL = 0.282094

def I_Y(d, m=M_SP):
    """I_Y dla trojkata rownobocznego."""
    D = _A1*_A2 + _A1*_A3 + _A2*_A3
    Q = (_A1+_A2+_A3) * d**2  # = d^2 bo suma = 1
    good = (D > 1e-30) & (Q > 1e-30)
    u = np.where(good, m*np.sqrt(Q/D), 1.0)
    val = np.where(good, D**(-1.5)*K0(u), 0.0)
    return 2.0*np.sum(_WW*_JAC*val)

def V2(d, C, m=M_SP):
    return -3.0*C**2*np.exp(-m*d)/d

def V3(d, C, m=M_SP):
    if d < 0.05: return -1e10
    return -C**3 * I_Y(d, m)

def E_ZP(d, N=3):
    """Zero-point energy dla N cząstek w trojkacie: N * hbar^2/(8m d^2)."""
    return N / (8.0 * d**2)

def E_eff_2b(d, C, m=M_SP):
    """E_eff = ZP + V2 (bez V3)."""
    return E_ZP(d) + V2(d, C, m)

def E_eff_3b(d, C, m=M_SP):
    """E_eff = ZP + V2 + V3."""
    return E_ZP(d) + V2(d, C, m) + V3(d, C, m)

def E_min_2b(C, m=M_SP):
    """Minimalna E_eff(2B) po d."""
    res = minimize_scalar(lambda d: E_eff_2b(d, C, m),
                          bounds=(0.05, 8.0), method='bounded')
    return res.fun, res.x

def E_min_3b(C, m=M_SP):
    """Minimalna E_eff(2B+3B) po d."""
    res = minimize_scalar(lambda d: E_eff_3b(d, C, m),
                          bounds=(0.05, 8.0), method='bounded')
    return res.fun, res.x


# ════════════════════════════════════════════════════════════════
# SEKCJA A: E_eff(d) dla C=C_Pl z i bez ZP
# ════════════════════════════════════════════════════════════════

print("="*65)
print("SEKCJA A: E_eff(d) dla C=C_Pl — efekt zero-point energy")
print("="*65)
print(f"C = {C_PL:.4f}, m_sp = {M_SP}, E_ZP = 3/(8d^2)")
print()
print(f"{'d':>6}  {'E_ZP':>10}  {'V2':>10}  {'V3':>10}  {'2B+ZP':>10}  {'2B+3B+ZP':>12}")
print("-"*65)

for d in [0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0]:
    ezp = E_ZP(d)
    v2  = V2(d, C_PL)
    v3  = V3(d, C_PL)
    e2  = ezp + v2
    e3  = ezp + v2 + v3
    print(f"{d:6.2f}  {ezp:10.4f}  {v2:10.4f}  {v3:10.4f}  {e2:10.4f}  {e3:12.4f}")

E_min_2b_PL, d_min_2b_PL = E_min_2b(C_PL)
E_min_3b_PL, d_min_3b_PL = E_min_3b(C_PL)
print()
print(f"Minimum 2B+ZP:    E={E_min_2b_PL:.5f} @ d={d_min_2b_PL:.4f}")
print(f"Minimum 2B+3B+ZP: E={E_min_3b_PL:.5f} @ d={d_min_3b_PL:.4f}")
print(f"Delta E (V3 wklad): {E_min_3b_PL-E_min_2b_PL:.5f} E_Pl")


# ════════════════════════════════════════════════════════════════
# SEKCJA B: C_crit^QM — gdzie 2B+ZP = 0 i 2B+3B+ZP = 0
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("SEKCJA B: C_crit^QM — przejscie wiazany/niewiazany")
print("="*65)
print()

# Skan C
C_scan = np.linspace(0.001, 0.50, 120)
E2_scan = []; E3_scan = []; d2_scan = []; d3_scan = []

for C in C_scan:
    e2, d2 = E_min_2b(C)
    e3, d3 = E_min_3b(C)
    E2_scan.append(e2); d2_scan.append(d2)
    E3_scan.append(e3); d3_scan.append(d3)

E2_scan = np.array(E2_scan)
E3_scan = np.array(E3_scan)

# Szukaj C_crit^(2B): gdzie E_min(2B) = 0
try:
    # Sprawdz czy jest zmiana znaku
    if np.any(E2_scan < 0) and np.any(E2_scan > 0):
        idx = np.where(np.diff(np.sign(E2_scan)))[0][0]
        C_crit_2b = brentq(lambda C: E_min_2b(C)[0],
                           C_scan[idx], C_scan[idx+1], xtol=1e-5)
        print(f"C_crit^QM(2B)    = {C_crit_2b:.6f}  (ponizej: niewiazany 2B+ZP)")
    else:
        C_crit_2b = None
        sgn = "ujemne (<0)" if E2_scan[-1] < 0 else "dodatnie (>0)"
        print(f"E_min(2B+ZP) zawsze {sgn} dla C=0..0.5")
        # sprawdz dla malego C
        for C_test in [1e-4, 1e-3, 0.01, 0.05]:
            e, d = E_min_2b(C_test)
            print(f"  C={C_test:.4f}: E_min={e:.5f} @ d={d:.4f}")
except Exception as ex:
    C_crit_2b = None
    print(f"Szukanie C_crit^QM(2B): {ex}")

# Szukaj C_crit^(3B): gdzie E_min(2B+3B) = 0
try:
    if np.any(E3_scan < 0) and np.any(E3_scan > 0):
        idx = np.where(np.diff(np.sign(E3_scan)))[0][0]
        C_crit_3b = brentq(lambda C: E_min_3b(C)[0],
                           C_scan[idx], C_scan[idx+1], xtol=1e-5)
        print(f"C_crit^QM(2B+3B) = {C_crit_3b:.6f}  (ponizej: niewiazany 2B+3B+ZP)")
    else:
        sgn = "ujemne (<0)" if E3_scan[-1] < 0 else "dodatnie (>0)"
        print(f"E_min(2B+3B+ZP) zawsze {sgn} dla C=0..0.5")
except Exception as ex:
    print(f"Szukanie C_crit^QM(3B): {ex}")

# Szukaj jeszcze dla wiekszych C (wyżej) — może próg jest gdzie indziej
print()
print("Szczegółowy skan E_min vs C:")
print(f"{'C':>8}  {'E_min(2B)':>12}  {'d_eq(2B)':>10}  {'E_min(3B)':>12}  {'d_eq(3B)':>10}  {'3B_only?':>10}")
print("-"*75)

for C in [0.001, 0.002, 0.005, 0.010, 0.020, 0.050, 0.100, 0.150, 0.200, 0.250, 0.282, 0.350, 0.400]:
    e2, d2 = E_min_2b(C)
    e3, d3 = E_min_3b(C)
    only3b = "TAK!" if e2 > 0 and e3 < 0 else ""
    print(f"{C:8.3f}  {e2:12.5f}  {d2:10.4f}  {e3:12.5f}  {d3:10.4f}  {only3b:>10}")


# ════════════════════════════════════════════════════════════════
# SEKCJA C: Wplyw masy — inne m_sp
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("SEKCJA C: Wplyw m_sp na C_crit^QM")
print("="*65)
print()
print("Dla roznych m_sp: czy istnieje C_crit(2B)?")
print()
print(f"{'m_sp':>8}  {'E_min(2B,C=0.01)':>18}  {'d_eq':>8}  {'C_crit(2B)?':>14}")
print("-"*55)

def E_min_2b_m(C, m):
    res = minimize_scalar(lambda d: E_ZP(d) + V2(d, C, m),
                          bounds=(0.01, 20.0), method='bounded')
    return res.fun, res.x

for m_sp in [0.1, 0.3, 0.5, 1.0, 2.0, 5.0]:
    e_test, d_test = E_min_2b_m(0.01, m_sp)
    # Szukaj C_crit dla tego m_sp
    try:
        C_hi = 0.5; E_hi, _ = E_min_2b_m(C_hi, m_sp)
        C_lo = 0.001; E_lo, _ = E_min_2b_m(C_lo, m_sp)
        if E_lo > 0 and E_hi < 0:
            C_crit = brentq(lambda C: E_min_2b_m(C, m_sp)[0], C_lo, C_hi, xtol=1e-5)
            print(f"{m_sp:8.2f}  {e_test:18.5f}  {d_test:8.4f}  {C_crit:14.6f}")
        else:
            sgn = "> 0 (luzy)" if E_hi > 0 else "< 0 (zwiaz)"
            print(f"{m_sp:8.2f}  {e_test:18.5f}  {d_test:8.4f}  zawsze {sgn}")
    except Exception as ex:
        print(f"{m_sp:8.2f}  {e_test:18.5f}  {d_test:8.4f}  blad: {ex}")


# ════════════════════════════════════════════════════════════════
# SEKCJA D: Model 2D — zmniejszona liczba stopni swobody
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("SEKCJA D: Efektywny opis 1D — promien trojkata R=d/sqrt(3)")
print("="*65)
print()
print("E_ZP(1D) = 1/(2*mu*d^2) gdzie mu = zredukowana masa pary = 1/2")
print("E_ZP(1D, mu=1/2) = 1/d^2  (ostrzejszy niz 3/(8d^2))")
print()

def E_eff_1d(d, C, m=M_SP, alpha_zp=1.0):
    """E_eff z ostrą ZP (efektywna 1D): alpha_zp/d^2 + V2 + V3."""
    return alpha_zp/d**2 + V2(d, C, m) + V3(d, C, m)

def E_eff_1d_2b(d, C, m=M_SP, alpha_zp=1.0):
    return alpha_zp/d**2 + V2(d, C, m)

# Skan po alpha_ZP: od 3/8 (isotropic 3D) do 1/2 (1D par)
print(f"Skan C_crit^QM dla roznych alpha_ZP (E_ZP = alpha/d^2):")
print(f"{'alpha_ZP':>10}  {'C_crit(2B)':>12}  {'C_crit(3B)':>12}  {'Okno 3B-only':>15}")
print("-"*55)

for alpha_zp in [0.375, 0.5, 0.75, 1.0, 1.5, 2.0]:
    def e2_a(C):
        res = minimize_scalar(lambda d: alpha_zp/d**2 + V2(d,C),
                              bounds=(0.01, 10.0), method='bounded')
        return res.fun
    def e3_a(C):
        res = minimize_scalar(lambda d: alpha_zp/d**2 + V2(d,C) + V3(d,C),
                              bounds=(0.01, 10.0), method='bounded')
        return res.fun

    C_arr = np.linspace(0.001, 0.6, 80)
    e2_arr = np.array([e2_a(C) for C in C_arr])
    e3_arr = np.array([e3_a(C) for C in C_arr])

    # C_crit(2B)
    try:
        idx2 = np.where(np.diff(np.sign(e2_arr)))[0]
        if len(idx2) > 0:
            C_c2 = brentq(e2_a, C_arr[idx2[0]], C_arr[idx2[0]+1], xtol=1e-5)
        else:
            C_c2 = None
    except: C_c2 = None

    # C_crit(3B)
    try:
        idx3 = np.where(np.diff(np.sign(e3_arr)))[0]
        if len(idx3) > 0:
            C_c3 = brentq(e3_a, C_arr[idx3[0]], C_arr[idx3[0]+1], xtol=1e-5)
        else:
            C_c3 = None
    except: C_c3 = None

    c2_str = f"{C_c2:.5f}" if C_c2 else "zawsze <0"
    c3_str = f"{C_c3:.5f}" if C_c3 else "zawsze <0"

    if C_c3 and C_c2 and C_c3 < C_c2:
        window = f"[{C_c3:.4f},{C_c2:.4f}]"
    elif C_c2 and not C_c3:
        window = f"C < {C_c2:.4f}"
    else:
        window = "brak"

    print(f"{alpha_zp:10.3f}  {c2_str:>12}  {c3_str:>12}  {window:>15}")


# ════════════════════════════════════════════════════════════════
# SEKCJA E: Rownanie bohr-like dla d_eq
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("SEKCJA E: Promien Bohra trojkata TGP — d_eq vs C")
print("="*65)
print()
print("Rownanie rownowagi: dE_eff/dd = 0")
print("=> -6*alpha/d^3 + dV2/dd + dV3/dd = 0")
print()

alpha_3d = 3.0/8.0  # 3D isotropic ZP

print(f"{'C':>8}  {'d_eq(2B)':>10}  {'d_eq(2B+3B)':>14}  {'Delta_d':>10}  {'E_bind(2B+3B)':>16}")
print("-"*65)

for C in [0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.282, 0.35]:
    e2, d2 = E_min_2b(C)
    e3, d3 = E_min_3b(C)
    dd = d3 - d2
    print(f"{C:8.4f}  {d2:10.5f}  {d3:14.5f}  {dd:+10.5f}  {e3:16.6f}")

print()
print("Interpretacja:")
print("  d_eq(2B+3B) < d_eq(2B): trojkat z V3 jest bardziej zwarty")
print("  Delta_d < 0: sila 3-cialowa sciaga ciastka blizej siebie")


# ════════════════════════════════════════════════════════════════
# PODSUMOWANIE
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("PODSUMOWANIE ex22")
print("="*65)

E_2b_pl, d_2b_pl = E_min_2b(C_PL)
E_3b_pl, d_3b_pl = E_min_3b(C_PL)

print(f"""
Zero-point energy: E_ZP = 3/(8d^2) [Planck], trojkat rownoboczny

Dla C=C_Pl={C_PL:.4f}:
  d_eq (2B+ZP):    {d_2b_pl:.4f} l_Pl  =>  E_bind = {E_2b_pl:.5f} E_Pl
  d_eq (2B+3B+ZP): {d_3b_pl:.4f} l_Pl  =>  E_bind = {E_3b_pl:.5f} E_Pl
  Wzmocnienie V3:  {(E_3b_pl-E_2b_pl):.5f} E_Pl ({(E_3b_pl-E_2b_pl)/abs(E_2b_pl)*100:.1f}%)
  Kompresja d:     Delta_d = {d_3b_pl-d_2b_pl:+.4f} l_Pl

KLUCZOWY WYNIK:
  Dla m_sp=1 (Planck): E_min(2B+ZP) zawsze <0 => C_crit(2B) < mierzalny prog
  Dla m_sp < m_crit: moze istniec okno C w ktorym 3B tworzy stan zwiazany.

PREDYKCJA:
  Dla wiekszych alpha_ZP (ostrzejsze ZP, np. 1D): istnieje
  wyrazne OKNO C_crit(3B) < C < C_crit(2B) gdzie czysto 3-cialowy
  stan zwiazany pojawia sie w TGP.

NASTEPNY KROK:
  Pelna kwantowa kalkulacja (Schrodinger 2D) z V2+V3 jako potencjalem:
  wyznaczenie E_ground przez metode roznicowa/variacyjna.
""")
