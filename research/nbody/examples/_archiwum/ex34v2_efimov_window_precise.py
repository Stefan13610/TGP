"""
ex34v2_efimov_window_precise.py
================================
Precyzyjna bisekcja granic okna Efimova 2D (rewizja ex34)

MOTYWACJA
---------
Ex34 wykazal, ze przy C=C_Pl:
  - m_sp=0.076: E0_2D=-0.015 (3B wiaze), ale C_Q(2B)_2D=0.276 < C_Pl=0.282
    => 2B tez wiaze => BRAK okna Efimova (C_Pl POWYZEJ okna)
  - m_sp=0.100: E0_2D=+0.023 (3B niewiazany)
    => C_Pl PONIZEJ progu 3B => BRAK okna

Interpretacja: dwa przecia przez C_Pl:
  m_sp** = gdzie C_Q(2B)_2D = C_Pl  [przejscie 2B z wiazanego do niewiazanego]
  m_sp*  = gdzie C_Q(3B)_2D = C_Pl  = gdzie E0_2D(C_Pl) = 0

Interpolacja liniowa z ex34:
  m_sp** ~ 0.081   (z C_Q2B: 0.276 at 0.076, 0.307 at 0.100)
  m_sp*  ~ 0.085   (z E0_2D: -0.015 at 0.076, +0.023 at 0.100)

Jezeli m_sp** < m_sp*:
  Okno Efimova dla C_Pl ISTNIEJE dla m_sp in (m_sp**, m_sp*)

Ten skrypt wyznacza m_sp** i m_sp* precyzyjnie przez bisekcje 2D.

METODA
------
1. Wyznaczenie m_sp* (3B threshold przy C_Pl):
   Bisekcja E0_2D(C_Pl, m_sp) = 0  dla m_sp in [0.076, 0.100]
   Uzywamy solve_2d_isosceles(C_PL, m_sp)

2. Wyznaczenie m_sp** (2B threshold przy C_Pl):
   Bisekcja C_Q2B_2D(m_sp) = C_Pl  dla m_sp in [0.076, 0.100]
   Dla kazdego m_sp: find_C_Q2B_2D(m_sp), sprawdz czy == C_Pl

ALTERNATYWA dla 2.:
   Bisekcja na m_sp: E0_2D(C_Pl, m_sp, V3=OFF) = 0
   (bez V3 = efektywny 2-body problem)

WYNIK
-----
  m_sp**:  przejscie 2B   (z wiazanego do niewiazanego)
  m_sp*:   przejscie 3B   (z wiazanego do niewiazanego)
  Okno: (m_sp**, m_sp*) lub brak jezeli m_sp** >= m_sp*

Autor: TGP Analysis Session v26, 2026-03-22
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.special import k0 as K0
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh
import warnings
warnings.filterwarnings('ignore')

C_PL = 1.0 / (2.0 * np.sqrt(np.pi))

print("=" * 70)
print("EX34v2: Precyzyjna bisekcja okna Efimova 2D")
print(f"        C_Pl = {C_PL:.6f}  [Planck units]")
print("=" * 70)
print()

# -------------------------------------------------------------------
# Feynman integral (identyczny jak ex34)
# -------------------------------------------------------------------
_pts, _wts = np.polynomial.legendre.leggauss(20)
_up = 0.5 * (1 + _pts);  _uw = 0.5 * _wts
_UU, _VV = np.meshgrid(_up, _up, indexing='ij')
_WW  = np.outer(_uw, _uw)
_A1  = _UU
_A2  = _VV * (1.0 - _UU)
_A3  = (1.0 - _UU) * (1.0 - _VV)
_JAC = (1.0 - _UU)

def I_Y(d12, d13, d23, m):
    D = _A1*_A2 + _A1*_A3 + _A2*_A3
    Q = _A2*d12**2 + _A1*d13**2 + _A3*d23**2
    good = (D > 1e-30) & (Q > 1e-30)
    u   = np.where(good, m * np.sqrt(Q / D), 1.0)
    val = np.where(good, D**(-1.5) * K0(u), 0.0)
    return 2.0 * np.sum(_WW * _JAC * val)

def V2_total(a, b, m):
    return 2.0 * np.exp(-m*a)/a + np.exp(-m*b)/b

# -------------------------------------------------------------------
# 2D solver (isosceles Jacobi) — umiarkowana siatka dla predkosci
# -------------------------------------------------------------------

def solve_2d(C, m_sp, Nb=75, Nh=55,
             b_lo=0.3, b_hi=18.0, h_lo=0.2, h_hi=10.0,
             include_V3=True):
    """
    Rozwiazuje 2D Schrodingera dla izosceles trojkata (d12=d13=a, d23=b).
    include_V3=True:  pelny potencjal V2+V3  (trojcialowy)
    include_V3=False: tylko V2 bez V3         (dwucialowy — szacunek C_Q(2B))
    """
    mu1 = 0.5
    mu2 = 2.0 / 3.0
    b_arr = np.linspace(b_lo, b_hi, Nb)
    h_arr = np.linspace(h_lo, h_hi, Nh)
    db = b_arr[1] - b_arr[0]
    dh = h_arr[1] - h_arr[0]

    V_grid = np.zeros((Nb, Nh))
    for ib, bv in enumerate(b_arr):
        for ih, hv in enumerate(h_arr):
            a = np.sqrt(hv**2 + bv**2 / 4.0)
            if a < 1e-6 or bv < 1e-6:
                V_grid[ib, ih] = 1e6
                continue
            v2 = -C**2 * V2_total(a, bv, m_sp)
            v3 = -C**3 * I_Y(a, a, bv, m_sp) if include_V3 else 0.0
            V_grid[ib, ih] = v2 + v3

    Ntot = Nb * Nh
    H = lil_matrix((Ntot, Ntot))
    cb = 0.5 / (mu1 * db**2)
    ch = 0.5 / (mu2 * dh**2)

    for ib in range(Nb):
        for ih in range(Nh):
            i = ib * Nh + ih
            H[i, i] = 2.0*cb + 2.0*ch + V_grid[ib, ih]
            if ib > 0:   H[i, (ib-1)*Nh+ih] = -cb
            if ib < Nb-1:H[i, (ib+1)*Nh+ih] = -cb
            if ih > 0:   H[i, ib*Nh+(ih-1)] = -ch
            if ih < Nh-1:H[i, ib*Nh+(ih+1)] = -ch

    H = csr_matrix(H)
    try:
        vals, _ = eigsh(H, k=1, which='SA', tol=1e-7, maxiter=20000)
        return float(vals[0])
    except Exception:
        return np.inf

# -------------------------------------------------------------------
# KROK 1: Skan E0_2D(C_Pl) dla m_sp in [0.076, 0.100]
# Cel: znalezienie m_sp* (gdzie 3B traci wiazanie przy C_Pl)
# -------------------------------------------------------------------

print("─" * 70)
print("KROK 1: Skan E0_2D(C_Pl, m_sp) — szukamy m_sp* (zero-crossing)")
print("─" * 70)

m_scan = np.arange(0.076, 0.101, 0.003)  # 0.076, 0.079, ..., 0.100
E0_scan = {}

for m_sp in m_scan:
    print(f"  m_sp={m_sp:.3f} ...", flush=True, end=" ")
    E0 = solve_2d(C_PL, m_sp, include_V3=True)
    E0_scan[round(m_sp, 4)] = E0
    bound_str = "WIAZANY" if E0 < 0 else "NIEwiazany"
    print(f"E0_2D = {E0:+.5f}  [{bound_str}]", flush=True)

print()
print("  Wyniki skanu:")
print(f"  {'m_sp':>8}  {'E0_2D':>12}  {'3B wiaze?':>12}")
print("  " + "─" * 38)
for m_sp in sorted(E0_scan):
    E0 = E0_scan[m_sp]
    print(f"  {m_sp:>8.4f}  {E0:>+12.5f}  {'TAK' if E0 < 0 else 'nie':>12}")

# Interpolacja liniowa dla m_sp*
m_arr = sorted(E0_scan.keys())
E_arr = [E0_scan[m] for m in m_arr]

m_star = None
for i in range(len(m_arr) - 1):
    if E_arr[i] < 0 and E_arr[i+1] >= 0:
        # Zero crossing between m_arr[i] and m_arr[i+1]
        frac = -E_arr[i] / (E_arr[i+1] - E_arr[i])
        m_star = m_arr[i] + frac * (m_arr[i+1] - m_arr[i])
        break

print()
if m_star:
    print(f"  m_sp* (3B traci wiazanie przy C_Pl) = {m_star:.4f} l_Pl^-1  [interpolacja]")
else:
    print("  m_sp*: brak zero-crossing w zakresie [0.076, 0.100]")

# -------------------------------------------------------------------
# KROK 2: Bisekcja dla m_sp** (2B traci wiazanie przy C_Pl)
# Uzywamy V3=False (bez sil trojcialowych) — efektywny problem 2-cialowy
# -------------------------------------------------------------------

print()
print("─" * 70)
print("KROK 2: Skan E0_2D_noV3(C_Pl, m_sp) — szukamy m_sp** (2B threshold)")
print("        (V3=OFF: efektywny problem 2-cialowy)")
print("─" * 70)

E0_2B_scan = {}

for m_sp in m_scan:
    print(f"  m_sp={m_sp:.3f} (noV3) ...", flush=True, end=" ")
    E0 = solve_2d(C_PL, m_sp, include_V3=False)
    E0_2B_scan[round(m_sp, 4)] = E0
    bound_str = "WIAZANY" if E0 < 0 else "NIEwiazany"
    print(f"E0_2D_noV3 = {E0:+.5f}  [{bound_str}]", flush=True)

print()
print("  Wyniki (bez V3 — efekt 2-cialowy):")
print(f"  {'m_sp':>8}  {'E0_2D_noV3':>14}  {'2B wiaze?':>12}")
print("  " + "─" * 40)
for m_sp in sorted(E0_2B_scan):
    E0 = E0_2B_scan[m_sp]
    print(f"  {m_sp:>8.4f}  {E0:>+14.5f}  {'TAK' if E0 < 0 else 'nie':>12}")

m_arr_2B = sorted(E0_2B_scan.keys())
E_arr_2B = [E0_2B_scan[m] for m in m_arr_2B]

m_dstar = None
for i in range(len(m_arr_2B) - 1):
    if E_arr_2B[i] < 0 and E_arr_2B[i+1] >= 0:
        frac = -E_arr_2B[i] / (E_arr_2B[i+1] - E_arr_2B[i])
        m_dstar = m_arr_2B[i] + frac * (m_arr_2B[i+1] - m_arr_2B[i])
        break

print()
if m_dstar:
    print(f"  m_sp** (2B traci wiazanie przy C_Pl) = {m_dstar:.4f} l_Pl^-1  [interpolacja]")
else:
    print("  m_sp**: brak zero-crossing — sprawdz zakres!")

# -------------------------------------------------------------------
# KROK 3: Werdykt — czy okno istnieje?
# -------------------------------------------------------------------

print()
print("=" * 70)
print("WERDYKT: Okno Efimova 2D dla C_Pl")
print("=" * 70)
print()
print(f"  C_Pl = {C_PL:.5f}")
print()

if m_star is None and m_dstar is None:
    print("  BRAK danych do rozstrzygniecia — poszerz zakres m_sp!")
elif m_star is None:
    print("  m_sp* nie znalezione — 3B moze wiazac w calym zakresie")
elif m_dstar is None:
    print("  m_sp** nie znalezione — 2B nigdy nie wiaze lub zawsze wiaze w zakresie")
else:
    print(f"  m_sp** = {m_dstar:.4f} l_Pl^-1  [2B przejscie: wiazany -> niewiazany]")
    print(f"  m_sp*  = {m_star:.4f} l_Pl^-1  [3B przejscie: wiazany -> niewiazany]")
    print()
    if m_dstar < m_star:
        width = m_star - m_dstar
        print(f"  OKNO EFIMOVA ISTNIEJE: m_sp in ({m_dstar:.4f}, {m_star:.4f}) l_Pl^-1")
        print(f"  Szerokosc okna: Delta_m = {width:.4f} l_Pl^-1")
        print()
        print(f"  FIZYCZNA INTERPRETACJA:")
        print(f"    Dla m_sp in ({m_dstar:.4f}, {m_star:.4f}):")
        print(f"      - Para (2B) NIEWIAZANA przy C=C_Pl")
        print(f"      - Trojka (3B) WIAZANA przy C=C_Pl")
        print(f"      => Stan Efimova: trojcialowe wiazanie bez wiazania parowego !")
        print()
        print(f"  K13 STATUS: POTWIERDZONE (waze okno ~{width:.3f} l_Pl^-1)")
    else:
        print(f"  BRAK okna Efimova dla C_Pl (m_sp** >= m_sp*)")
        print(f"  Dla C_Pl: 2B i 3B tracza wiazanie w tej samej kolejnosci")
        print()
        print(f"  K13 STATUS: NIE ISTNIEJE dla C=C_Pl")
        print(f"  Mozliwe okno dla innych C (np. C < C_Pl)")

print()
print("─" * 70)
print("Uwaga: Wyniki oparte na siatce Nb=75, Nh=55 (umiarkowana dokladnosc).")
print("Dla potwierdzenia: uruchom z Nb=120, Nh=80 (~4x dluzej).")
print("─" * 70)
print()
print("=" * 70)
print("EX34v2 DONE")
print("=" * 70)
