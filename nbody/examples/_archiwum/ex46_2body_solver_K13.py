"""
ex46_2body_solver_K13.py
========================
Dedykowany solver 2-ciałowy TGP — precyzyjne m_sp** dla K13

MOTYWACJA
---------
Ex34v2 znalazł m_sp* = 0.0831 l_Pl⁻¹ (próg 3B, z bisekcji E0_2D=0).
Dla m_sp** (próg 2B) ex34v2 użył solwera izosceles V3=OFF, który
jest METODOLOGICZNIE NIEADEKWATNY: geometria izosceles wymusza
d12=d13=a przy b→∞, co daje a→∞ — nie można wyizolować pary (1,2).

Ten skrypt implementuje POPRAWNY solver 2-ciałowy:
  -1/(2μ) ψ''(d) + V_2B(d; C, m_sp) ψ(d) = E ψ(d)
gdzie μ = 1/2 (masa zredukowana dla dwóch jednakowych ciał m=1)
i V_2B = -C²*exp(-m_sp*d)/d + β*C²*exp(-m_sp*d)/d²

Jest to 1D problem własny na siatce FD z warunkami Dirichleta.

METODA
------
1. Dla każdego m_sp w [0.070, 0.095]:
   a. Bisekcja w C: znajdź C_Q(2B) gdzie E0_2B(C) = 0
   b. Siatka FD: N=3000 punktów, d ∈ [0.3, 40] l_Pl
2. Interpolacja: znajdź m_sp** gdzie C_Q(2B)(m_sp**) = C_Pl
3. Porównanie z m_sp* = 0.0831 (z ex34v2): potwierdź okno Efimova

WYNIK
-----
  m_sp**:  próg 2B (dokładny)
  m_sp*:   próg 3B (z ex34v2, stały = 0.0831)
  Okno: (m_sp**, m_sp*) lub brak
  Status K13: POTWIERDZONY lub BRAK OKNA

Autor: TGP Analysis Session v27, 2026-03-22
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh
import warnings
warnings.filterwarnings('ignore')

C_PL = 1.0 / (2.0 * np.sqrt(np.pi))   # = 0.28209...
BETA = 1.0   # = gamma (warunek próżniowy N0-5)
MU   = 0.5   # masa zredukowana: μ = m₁m₂/(m₁+m₂) = 1/2 dla m₁=m₂=1

print("=" * 70)
print("EX46: Dedykowany solver 2-ciałowy TGP — precyzyjne m_sp** (K13)")
print(f"      C_Pl = {C_PL:.6f},  mu = {MU},  beta = {BETA}")
print("=" * 70)
print()

# -------------------------------------------------------------------
# 1D Schrödinger FD dla problemu 2-ciałowego
# -------------------------------------------------------------------

def V_2body(d, C, m_sp):
    """Potencjał 2-ciałowy TGP (bez bariery zerowego punktu — czysty 2B)."""
    if d < 1e-8:
        return 1e8
    em = np.exp(-m_sp * d)
    v_attr = -C**2 * em / d
    v_rep  = +C**2 * BETA * em / d**2
    return v_attr + v_rep

def V_2body_array(d_arr, C, m_sp):
    """Wektoryzowany potencjał 2-ciałowy."""
    em = np.exp(-m_sp * d_arr)
    v  = -C**2 * em / d_arr + C**2 * BETA * em / d_arr**2
    return v

def solve_2body_E0(C, m_sp, N=2000, d_lo=0.3, d_hi=40.0):
    """
    Rozwiązuje 1D Schrödingera dla pary TGP.
    Zwraca najniższy poziom energii E0 [E_Pl].
    E0 < 0: para ZWIĄZANA; E0 > 0: para NIEZWIĄZANA.
    """
    d_arr = np.linspace(d_lo, d_hi, N)
    dd = d_arr[1] - d_arr[0]

    V = V_2body_array(d_arr, C, m_sp)
    V = np.clip(V, -1e4, 1e4)

    # Hamiltonian FD
    diag  = 1.0 / (MU * dd**2) + V
    off   = -0.5 / (MU * dd**2) * np.ones(N - 1)

    from scipy.sparse import diags as sp_diags
    H = sp_diags([off, diag, off], [-1, 0, 1], format='csr')

    try:
        vals, _ = eigsh(H, k=1, which='SA', tol=1e-8, maxiter=10000)
        return float(vals[0])
    except Exception:
        return np.inf

def find_C_Q2B(m_sp, C_lo=0.10, C_hi=0.60, tol=1e-5, max_iter=50):
    """
    Bisekcja: znajdź C_Q(2B)(m_sp) = najniższe C gdzie E0_2B < 0.
    Zwraca C_Q(2B) lub None jeśli nie znaleziono w [C_lo, C_hi].
    """
    E_lo = solve_2body_E0(C_lo, m_sp)
    E_hi = solve_2body_E0(C_hi, m_sp)

    if E_lo < 0:
        return C_lo  # Para wiąże już przy C_lo — prawdopodobnie za małe C_lo
    if E_hi >= 0:
        return None  # Para nie wiąże nawet przy C_hi

    # Bisekcja
    for _ in range(max_iter):
        C_mid = 0.5 * (C_lo + C_hi)
        E_mid = solve_2body_E0(C_mid, m_sp)
        if E_mid < 0:
            C_hi = C_mid
        else:
            C_lo = C_mid
        if C_hi - C_lo < tol:
            break

    return 0.5 * (C_lo + C_hi)

# -------------------------------------------------------------------
# KROK 1: Skan C_Q(2B) dla m_sp ∈ [0.065, 0.100]
# -------------------------------------------------------------------

print("─" * 70)
print("KROK 1: Skan C_Q(2B)(m_sp) — próg wiązania pary")
print("─" * 70)
print()

m_scan = np.arange(0.065, 0.101, 0.005)
CQ2B = {}

print(f"  {'m_sp':>8}  {'C_Q(2B)':>10}  {'C_Pl w oknie 2B?':>20}")
print("  " + "─" * 44)

for m_sp in m_scan:
    cq = find_C_Q2B(round(m_sp, 4))
    CQ2B[round(m_sp, 4)] = cq
    if cq is not None:
        in_window = "C_Pl > C_Q(2B) [WIAŻE]" if C_PL > cq else "C_Pl < C_Q(2B) [nie wiaze]"
        print(f"  {m_sp:>8.4f}  {cq:>10.5f}  {in_window:>20}")
    else:
        print(f"  {m_sp:>8.4f}  {'N/A':>10}  {'nie znaleziono':>20}")

# -------------------------------------------------------------------
# KROK 2: Interpolacja m_sp** gdzie C_Q(2B) = C_Pl
# -------------------------------------------------------------------

print()
print("─" * 70)
print("KROK 2: Interpolacja m_sp** gdzie C_Q(2B)(m_sp**) = C_Pl")
print("─" * 70)
print()

# Zbierz pary (m_sp, C_Q(2B)) gdzie C_Q(2B) istnieje
m_arr  = sorted(m for m in CQ2B if CQ2B[m] is not None)
cq_arr = [CQ2B[m] for m in m_arr]

# Szukaj przejścia C_Q(2B) przez C_Pl
m_dstar = None
for i in range(len(m_arr) - 1):
    if cq_arr[i] is not None and cq_arr[i+1] is not None:
        # C_Q(2B) rośnie z m_sp (para trudniej się wiąże przy większym m_sp)
        # m_sp**: C_Q(2B)(m_sp**) = C_Pl
        # Czyli: C_Q(2B)(m_sp) - C_Pl zmienia znak
        f_i   = cq_arr[i]   - C_PL
        f_ip1 = cq_arr[i+1] - C_PL
        if f_i * f_ip1 < 0:
            frac = -f_i / (f_ip1 - f_i)
            m_dstar = m_arr[i] + frac * (m_arr[i+1] - m_arr[i])
            print(f"  Zero-crossing: C_Q(2B) przechodzi przez C_Pl między")
            print(f"    m_sp = {m_arr[i]:.4f} (C_Q2B={cq_arr[i]:.5f})")
            print(f"    m_sp = {m_arr[i+1]:.4f} (C_Q2B={cq_arr[i+1]:.5f})")
            break

# Jeśli nie znaleziono zero-crossing, spróbuj finer bisekcji w okolicach spodziewanego m_sp** ~ 0.08
if m_dstar is None:
    print("  Brak zero-crossing w grubym skanie. Próbuję finer scan wokół m_sp ~ 0.076-0.085...")
    m_fine = np.arange(0.072, 0.088, 0.002)
    cq_fine = []
    for m_sp in m_fine:
        cq = find_C_Q2B(round(m_sp, 4), C_lo=0.15, C_hi=0.50)
        cq_fine.append(cq)
        status = f"{cq:.5f}" if cq else "N/A"
        print(f"    m_sp={m_sp:.3f}: C_Q(2B)={status}")

    for i in range(len(m_fine) - 1):
        if cq_fine[i] and cq_fine[i+1]:
            f_i   = cq_fine[i]   - C_PL
            f_ip1 = cq_fine[i+1] - C_PL
            if f_i * f_ip1 < 0:
                frac = -f_i / (f_ip1 - f_i)
                m_dstar = float(m_fine[i]) + frac * float(m_fine[i+1] - m_fine[i])
                print(f"\n  Zero-crossing znalezione!")
                break

# -------------------------------------------------------------------
# KROK 3: Werdykt — okno Efimova 2D
# -------------------------------------------------------------------

print()
print("=" * 70)
print("WERDYKT: Okno Efimova 2D (solver 2-ciałowy)")
print("=" * 70)
print()

m_star = 0.0831  # z ex34v2 (bisekcja E0_3D=0)

print(f"  C_Pl     = {C_PL:.6f}")
print(f"  m_sp*    = {m_star:.4f} l_Pl^-1  [3B próg, ex34v2]")

if m_dstar is not None:
    print(f"  m_sp**   = {m_dstar:.4f} l_Pl^-1  [2B próg, ex46 DOKŁADNY]")
    print()

    if m_dstar < m_star:
        width = m_star - m_dstar
        print(f"  OKNO EFIMOVA ISTNIEJE: m_sp ∈ ({m_dstar:.4f}, {m_star:.4f}) l_Pl^-1")
        print(f"  Szerokość: Δm_sp = {width:.4f} l_Pl^-1")
        print()
        print(f"  K13 STATUS: POTWIERDZONY ✓")
        print(f"  Efimov window przy C_Pl: para niewiązana, trójka związana")
        print()
        print(f"  Porównanie z ex34 interpolacją:")
        print(f"    ex34 interp:  m_sp** ≈ 0.0806")
        print(f"    ex46 dokładny: m_sp** = {m_dstar:.4f}")
        diff = abs(m_dstar - 0.0806)
        print(f"    Różnica: {diff:.4f} l_Pl^-1 ({diff/0.0806*100:.1f}%)")
    else:
        print(f"  BRAK okna Efimova (m_sp** >= m_sp*)")
        print(f"  K13 STATUS: NIEZATWIERDZONY")
else:
    print(f"  m_sp** = NIE ZNALEZIONE w zakresie skanowania")
    print(f"  Możliwe: C_Q(2B) < C_Pl dla wszystkich m_sp w zakresie")
    print(f"           => para wiąże się zawsze przy C_Pl (brak okna Efimova)")
    print(f"  Lub:     C_Q(2B) > C_Pl dla wszystkich m_sp (para nigdy nie wiąże)")
    print()
    # Diagnoza z dostępnych danych
    vals_below = [(m, c) for m, c in CQ2B.items() if c is not None and c < C_PL]
    vals_above = [(m, c) for m, c in CQ2B.items() if c is not None and c > C_PL]
    if vals_below:
        print(f"  Para WIAŻE (C_Q2B < C_Pl) dla: {[m for m,c in vals_below]}")
    if vals_above:
        print(f"  Para NIE WIAŻE (C_Q2B > C_Pl) dla: {[m for m,c in vals_above]}")

print()
print("─" * 70)
print(f"Uwaga: N=2000 punktów siatki FD, d ∈ [0.3, 40] l_Pl.")
print(f"Dla wyższej dokładności: N=5000, d_hi=80.")
print("─" * 70)
print()
print("=" * 70)
print("EX46 DONE — 2-body solver for K13")
print("=" * 70)
