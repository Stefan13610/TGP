"""
ex24_efimov_tower.py - Wieża stanów Efimova TGP

W fizyce jądrowej stany Efimova tworzą NIESKONCZONA wieżę:
  E_n ~ E_0 * exp(-2*pi*n / s_0)
  Stosunek kolejnych energii: E_{n+1}/E_n = exp(-2*pi/s_0) ≈ 1/515 (dla bozonow)

W TGP: szukamy wzbudzonych stanów trójkąta poprzez:
  1. WKB (quasi-klasyczne): kwantowanie Bohra-Sommerfelda
     ∮ p(d) dd = (n + 1/2) * 2*pi*hbar
     E_n wyznaczone przez przejście klasyczne d_L < d < d_R (turning points)

  2. Numeryczna siatka: dyskretyzacja E_eff(d) i diagonalizacja
     hamiltoniani 1D: H = -d^2/(2*mu*dd^2) + E_eff(d)
     gdzie mu = zredukowana masa trojkata

  3. WKB integral: n(E) = (1/pi) * integral_dL^dR sqrt(2*mu*(E-V(d))) dd

Uwaga: E_eff(d) = alpha_ZP/d^2 + V2(d,C,m) + V3(d,C,m) to potencjal efektywny
w zmiennej d (bok trojkata). Stopień swobody 'r' jest quasi-1D.
"""

import numpy as np
from scipy.special import k0 as K0
from scipy.optimize import brentq
from scipy.integrate import quad
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# ─── Feynman integral ─────────────────────────────────────────────────────────
_pts, _wts = np.polynomial.legendre.leggauss(25)
_up = 0.5*(1+_pts); _uw = 0.5*_wts
_vp = 0.5*(1+_pts); _vw = 0.5*_wts
_UU, _VV = np.meshgrid(_up, _vp, indexing='ij')
_WW      = np.outer(_uw, _vw)
_A1 = _UU; _A2 = _VV*(1-_UU); _A3 = (1-_UU)*(1-_VV); _JAC = (1-_UU)

def I_Y_equil(d, m):
    D  = _A1*_A2 + _A1*_A3 + _A2*_A3
    Q  = d**2   # = d^2*(a1+a2+a3) = d^2
    good = (D > 1e-30)
    u  = np.where(good, m*np.sqrt(Q/D), 1.0)
    v  = np.where(good, D**(-1.5)*K0(u), 0.0)
    return 2.0*np.sum(_WW*_JAC*v)

# LUT: tabelaryzacja V3(d) dla zadanego m_sp
def build_V3_lut(m, d_lo=0.3, d_hi=60.0, n=300):
    d_arr = np.linspace(d_lo, d_hi, n)
    V3_arr = np.array([-1.0*I_Y_equil(d, m) for d in d_arr])  # per C^3
    from scipy.interpolate import interp1d
    return interp1d(d_arr, V3_arr, kind='cubic', fill_value='extrapolate'), d_arr

def V2_unit(d, m):
    """V2/C^2 = -3*exp(-m*d)/d"""
    return -3.0*np.exp(-m*d)/d

# ─── Potencjal efektywny ───────────────────────────────────────────────────────

def V_eff(d, C, m, alpha_ZP, V3_lut):
    """E_eff = alpha/d^2 + C^2*V2_unit + C^3*V3_unit"""
    return alpha_ZP/d**2 + C**2*V2_unit(d, m) + C**3*V3_lut(d)

# ─── WKB liczba stanów ────────────────────────────────────────────────────────

def wkb_n_states(E, C, m, alpha_ZP, V3_lut, mu=1.0):
    """
    n(E) = (1/pi) * integral_{d_L}^{d_R} sqrt(2*mu*(E - V_eff)) dd
    Gdzie d_L, d_R to klasyczne punkty zwrotne: V_eff(d) = E
    """
    d_arr = np.linspace(0.2, 60.0, 800)
    V_arr = np.array([V_eff(d, C, m, alpha_ZP, V3_lut) for d in d_arr])

    # Znajdz punkty zwrotne (gdzie V_eff = E)
    below = V_arr < E
    crossings = np.where(np.diff(below.astype(int)))[0]

    if len(crossings) < 2:
        return 0.0, None, None

    # Lewy punkt zwrotny
    idx_L = crossings[0]
    try:
        d_L = brentq(lambda d: V_eff(d,C,m,alpha_ZP,V3_lut) - E,
                     d_arr[idx_L], d_arr[idx_L+1])
    except:
        return 0.0, None, None

    # Prawy punkt zwrotny
    idx_R = crossings[-1]
    try:
        d_R = brentq(lambda d: V_eff(d,C,m,alpha_ZP,V3_lut) - E,
                     d_arr[idx_R], d_arr[idx_R+1])
    except:
        return 0.0, None, None

    if d_R <= d_L:
        return 0.0, d_L, d_R

    def integrand(d):
        val = E - V_eff(d, C, m, alpha_ZP, V3_lut)
        return np.sqrt(max(2*mu*val, 0))

    result, _ = quad(integrand, d_L, d_R, limit=100)
    return result/np.pi, d_L, d_R


def find_wkb_levels(C, m, alpha_ZP, V3_lut, mu=1.0, n_max=10):
    """
    Poziomy WKB: E_n takie ze n(E_n) = n + 1/2
    """
    # Znajdz minimum potencjału
    d_test = np.linspace(0.3, 55.0, 500)
    V_test = np.array([V_eff(d,C,m,alpha_ZP,V3_lut) for d in d_test])
    idx_min = np.argmin(V_test)
    E_min = V_test[idx_min]
    d_min = d_test[idx_min]

    if E_min >= 0:
        return [], d_min, E_min

    # Skan energii od E_min do 0
    E_scan = np.linspace(E_min*0.999, -1e-8, 300)
    n_scan = []
    for E in E_scan:
        n_val, dL, dR = wkb_n_states(E, C, m, alpha_ZP, V3_lut, mu)
        n_scan.append(n_val)
    n_scan = np.array(n_scan)

    levels = []
    for n_level in range(n_max):
        target = n_level + 0.5
        # znajdz E takie ze n_scan(E) = target
        above = n_scan >= target
        crossings = np.where(np.diff(above.astype(int)))[0]
        if len(crossings) == 0:
            break
        idx = crossings[0]
        try:
            E_level = brentq(
                lambda E: wkb_n_states(E,C,m,alpha_ZP,V3_lut,mu)[0] - target,
                E_scan[idx], E_scan[min(idx+1, len(E_scan)-1)]
            )
            _, dL, dR = wkb_n_states(E_level, C, m, alpha_ZP, V3_lut, mu)
            levels.append((n_level, E_level, dL, dR))
        except:
            break

    return levels, d_min, E_min


# ════════════════════════════════════════════════════════════════
# PARAMETRY: m_sp=0.1, C w środku okna Efimova
# ════════════════════════════════════════════════════════════════

M_SP    = 0.10
C_MID   = 0.180   # blisko C_crit(2B) -- gdzie sa stany WKB
ALPHA   = 3/8     # ZP dla trojkata 3D isotropowego
MU      = 1.0     # zredukowana masa (w jednostkach m_body)

print("Budowanie LUT V3(d) dla m_sp=0.10 ...")
V3_lut, d_lut = build_V3_lut(M_SP, d_lo=0.2, d_hi=60.0, n=400)
print(f"  LUT gotowy: d=[{d_lut[0]:.1f},{d_lut[-1]:.1f}], {len(d_lut)} pkt")
print()

# ════════════════════════════════════════════════════════════════
# SEKCJA A: Kształt potencjału V_eff(d) dla kilku C
# ════════════════════════════════════════════════════════════════

print("="*65)
print("SEKCJA A: Ksztalt potencjalu V_eff(d) dla m_sp=0.1")
print("="*65)
print()
print(f"{'d':>6}  " + "  ".join([f"C={c:.3f}" for c in [0.13, 0.155, 0.17, 0.184]]))
print("-"*65)

for d in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]:
    row = [f"{d:6.1f}"]
    for C in [0.13, 0.155, 0.17, 0.184]:
        v = V_eff(d, C, M_SP, ALPHA, V3_lut)
        row.append(f"{v:10.5f}")
    print("  ".join(row))

# ════════════════════════════════════════════════════════════════
# SEKCJA B: Poziomy WKB dla C=C_mid
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print(f"SEKCJA B: Poziomy WKB dla C={C_MID:.3f}, m_sp={M_SP}")
print("="*65)
print()

levels, d_min, E_min = find_wkb_levels(C_MID, M_SP, ALPHA, V3_lut, MU, n_max=8)

print(f"Minimum potencjalu: V_min={E_min:.6f} @ d={d_min:.3f} l_Pl")
print()

if levels:
    print(f"{'n':>4}  {'E_n [E_Pl]':>14}  {'E_n/E_0':>10}  {'d_L':>8}  {'d_R':>8}  {'Delta_d':>8}")
    print("-"*60)
    E0 = levels[0][1]
    for n, En, dL, dR in levels:
        ratio = En/E0 if E0 != 0 else 0
        print(f"{n:4d}  {En:14.6f}  {ratio:10.4f}  {dL:8.3f}  {dR:8.3f}  {dR-dL:8.3f}")

    # Stosunek kolejnych poziomów (geometric series?)
    if len(levels) >= 3:
        print()
        print("Stosunki kolejnych poziomów (test geometryczny):")
        for i in range(1, len(levels)):
            n0, E_prev, _, _ = levels[i-1]
            n1, E_curr, _, _ = levels[i]
            ratio = E_curr/E_prev if E_prev != 0 else float('nan')
            print(f"  E_{n1}/E_{n0} = {ratio:.4f}")
        print()
        if len(levels) >= 4:
            ratios = [levels[i][1]/levels[i-1][1] for i in range(1, len(levels))]
            mean_r = np.mean(ratios)
            std_r  = np.std(ratios)
            print(f"  Sredni stosunek: {mean_r:.4f} +/- {std_r:.4f}")
            print(f"  Geometryczny? {'TAK' if std_r/abs(mean_r) < 0.1 else 'NIE'}")
            print(f"  Dla Efimova (bosony): E_{{n+1}}/E_n ~ 1/515 = {1/515:.5f}")
else:
    print("Brak stanow WKB (potencjal nie ma studzienki?)")

# ════════════════════════════════════════════════════════════════
# SEKCJA C: Skan poziomów vs C — ile stanów w oknie?
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("SEKCJA C: Liczba poziomów WKB vs C (m_sp=0.1)")
print("="*65)
print()
print(f"{'C':>8}  {'N_states':>10}  {'E_0':>12}  {'E_top':>12}  {'d_eq':>8}")
print("-"*55)

for C in np.array([0.130, 0.135, 0.14, 0.145, 0.15, 0.155, 0.160, 0.165, 0.170, 0.175, 0.180, 0.184]):
    lvls, d_m, E_m = find_wkb_levels(C, M_SP, ALPHA, V3_lut, MU, n_max=15)
    N = len(lvls)
    if N > 0:
        E_gs  = lvls[0][1]
        E_top = lvls[-1][1]
        print(f"{C:8.4f}  {N:10d}  {E_gs:12.5f}  {E_top:12.5f}  {d_m:8.3f}")
    else:
        print(f"{C:8.4f}  {'0':>10}  {'---':>12}  {'---':>12}  {d_m:8.3f}")

# ════════════════════════════════════════════════════════════════
# SEKCJA D: Test skali Efimova — s_0 z TGP
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("SEKCJA D: Parametr skali Efimova s_0 dla TGP")
print("="*65)
print()
print("Dla stanow Efimova: E_n = E_0 * exp(-2*pi*n/s_0)")
print("Dopasowanie s_0 do poziomow WKB:")
print()

# Wybierz C które daje wiele poziomów
C_best = 0.165
lvls_best, _, _ = find_wkb_levels(C_best, M_SP, ALPHA, V3_lut, MU, n_max=12)

if len(lvls_best) >= 3:
    # Dopasowanie log(E_n) = log(E_0) - 2*pi*n/s_0
    ns  = np.array([l[0] for l in lvls_best])
    lnE = np.log(np.abs([l[1] for l in lvls_best]))
    # lnE = a - b*n => b = 2*pi/s_0
    from numpy.polynomial.polynomial import polyfit
    coeffs = np.polyfit(ns, lnE, 1)
    b = -coeffs[0]   # slope (powinien byc ujemny)
    s0_fit = 2*np.pi/b if b > 0 else float('nan')

    print(f"C = {C_best}, m_sp = {M_SP}")
    print(f"Dopasowanie log(|E_n|) = a - b*n:")
    print(f"  b  = {b:.4f}  (slope)")
    print(f"  s0 = 2*pi/b = {s0_fit:.4f}")
    print()
    print(f"  Dla atomowych stanow Efimova: s0 = 1.00624 (bosony)")
    print(f"  Stosunek E_{{n+1}}/E_n dla TGP: exp(-2*pi/s0) = {np.exp(-2*np.pi/s0_fit):.4f}" if not np.isnan(s0_fit) else "  s0 = nan")

    print()
    print(f"{'n':>4}  {'E_n':>12}  {'ln|E_n|':>10}  {'fit':>10}  {'blad':>8}")
    for i, (n, En, dL, dR) in enumerate(lvls_best):
        lnE_fit = coeffs[1] + coeffs[0]*n
        print(f"{n:4d}  {En:12.6f}  {np.log(abs(En)):10.4f}  {lnE_fit:10.4f}  {lnE[i]-lnE_fit:8.4f}")
else:
    print(f"Za malo poziomow dla C={C_best} (potrzeba >= 3)")

# ════════════════════════════════════════════════════════════════
# SEKCJA E: Porownanie z Efimovem jądrowym
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("SEKCJA E: TGP Efimov vs klasyczny Efimov jądrowy")
print("="*65)
print()

print(f"""
                     | Klasyczny Efimov      | TGP Efimov (ex23/ex24)
---------------------|----------------------|------------------------
Mechanizm            | 1/R^2 potencjal      | V3(d) Feynman integral
                     | hiper-radialny       | trojkat rownoboczny
---------------------|----------------------|------------------------
Warunek              | a -> +/-infty        | C in [C_crit(3B),
                     | (dlugosc rozpraszania| C_crit(2B)]
                     | rozbiezna)           |
---------------------|----------------------|------------------------
Skala E_{{n+1}}/E_n   | 1/515 (bosony)       | exp(-2*pi/s0_TGP)
                     |                      | (wyznaczone numerycznie)
---------------------|----------------------|------------------------
d_eq                 | wiele l_Pl           | 3-6 l_Pl (ex23)
---------------------|----------------------|------------------------
Falsyfikowalnosc     | obs. w ultra-zimnych | predykcja Planck-scale
                     | atomach (Cs, Li, He) | (hipotetyczne)
""")

# ════════════════════════════════════════════════════════════════
# PODSUMOWANIE
# ════════════════════════════════════════════════════════════════

print()
print("="*65)
print("PODSUMOWANIE ex24")
print("="*65)

# Najlepsze wyniki
lvls_sum, _, E_min_sum = find_wkb_levels(C_MID, M_SP, ALPHA, V3_lut, MU, n_max=8)
N_sum = len(lvls_sum)

print(f"""
Parametry: C={C_MID:.3f}, m_sp={M_SP}, alpha_ZP={ALPHA:.3f}

1. WIEŻA POZIOMÓW WKB:
   Liczba stanów związanych: {N_sum}
   Stan podstawowy E_0 = {(lvls_sum[0][1] if len(lvls_sum)>0 else 'N/A')} E_Pl
   Stan wzbudzony E_1 = {(lvls_sum[1][1] if len(lvls_sum)>1 else 'N/A')} E_Pl

2. GEOMETRYCZNA WIEŻA:
   Stosunek E_{{n+1}}/E_n ~ staly? => patrz Sekcja B
   Skalowanie Efimova: exp(-2*pi/s0_TGP)

3. MAKSYMALNA LICZBA STANOW w oknie (Sekcja C):
   Rosnie od 0 przy C~C_crit(3B) do maks. przy C~C_crit(2B)

4. ANALOGIA EFIMOVA:
   TGP ma analog stanow Efimova z:
   - Mechanizmem: V3 (Feynman 3-body) zamiast 1/R^2 hiper-radialnego
   - Skala: s0_TGP (numerycznie wyznaczona) vs s0=1.006 (Efimov)
   - Kontekst: Planck-scale obiekty zamiast ultra-zimnych atomow

5. PREDYKCJA FALSYFIKOWALNA:
   Jezeli m_sp < 0.4 l_Pl^-1 i m_body ~ 0.5 m_Pl:
   -- Istnieje skonczona wieża stanow trojcialowych
   -- Kazdy poziom ma skonczona energie wiazania i promien
   -- Brak analogicznych stanow 2-cialowych (czysto 3B)
""")
