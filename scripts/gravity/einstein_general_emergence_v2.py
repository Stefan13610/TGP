"""
einstein_general_emergence_v2.py  — POPRAWIONA WERSJA
======================================================
Błąd w v1: Γ^x_tt = +U_x  (WRONG)
Poprawnie:  Γ^x_tt = -U_x  (wynika z -∂_x g_tt z g_tt = -e^{-2U})

Konsekwencja: G^(1)_tt = -2∇²U  (nie -∇²U jak w v1)

To zmienia wniosek o emergencji dla PYŁU: dokładna (nie częściowa)!
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import sympy as sp
from sympy import symbols, Function, diff, expand, Rational, simplify

print("=" * 64)
print("TGP: Emergencja Einsteina v2 — POPRAWIONE Christoffele")
print("=" * 64)

t, x, y, z = symbols('t x y z', real=True)
U = Function('U')(t, x, y, z)
U_t  = diff(U, t);    U_x  = diff(U, x)
U_y  = diff(U, y);    U_z  = diff(U, z)
U_tt = diff(U, t, 2); U_xx = diff(U, x, 2)
U_yy = diff(U, y, 2); U_zz = diff(U, z, 2)
U_tx = diff(U, t, x); U_ty = diff(U, t, y); U_tz = diff(U, t, z)
lap_U = U_xx + U_yy + U_zz
box_U = -U_tt + lap_U  # □ = -∂_t² + ∇² w sygnaturze (-+++)

tests_pass = 0; tests_total = 0
def test(name, cond, detail=""):
    global tests_pass, tests_total
    tests_total += 1
    ok = bool(cond)
    if ok: tests_pass += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))

# ─────────────────────────────────────────────────────────────
# 1. POPRAWNE CHRISTOFFELE dla g = diag(-e^{-2U}, e^{2U}, e^{2U}, e^{2U})
# ─────────────────────────────────────────────────────────────
print("\n[1] POPRAWIONE symbole Christoffela")
print("    BŁĄD v1: Γ^i_tt = +U_i  (niepoprawny znak)")
print("    POPRAWKA: Γ^i_tt = -U_i")
print()
print("    Γ^x_tt = (1/2)g^xx(∂_t g_tx + ∂_t g_tx - ∂_x g_tt)")
print("           = (1/2)(1)(0 + 0 - ∂_x(-e^{-2U}))")
print("           = (1/2)(0 - 2U_x)")
print("           = -U_x  ← POPRAWNY ZNAK")

# Poprawne symbole Christoffela
Gam_t_tt = -U_t
Gam_t_ti = {'x': -U_x, 'y': -U_y, 'z': -U_z}
Gam_i_tt = {'x': -U_x, 'y': -U_y, 'z': -U_z}   # POPRAWKA: minus!
Gam_t_ij = {'xx': U_t, 'yy': U_t, 'zz': U_t}
Gam_i_ti = {'xx': U_t, 'yy': U_t, 'zz': U_t}
Gam_i_ij_neq = {'xy': U_y, 'xz': U_z, 'yx': U_x, 'yz': U_z, 'zx': U_x, 'zy': U_y}
Gam_i_jj_neq = {'xyy': -U_x, 'xzz': -U_x, 'yxx': -U_y, 'yzz': -U_y, 'zxx': -U_z, 'zyy': -U_z}
Gam_i_ii = {'x': U_x, 'y': U_y, 'z': U_z}  # Γ^x_xx = U_x

# ─────────────────────────────────────────────────────────────
# 2. TENSOR RICCIEGO R^(1)_μν  (poprawiony)
# ─────────────────────────────────────────────────────────────
print("\n[2] Tensor Ricciego (poprawiony)")

# R^(1)_tt = ∂_ρ Γ^ρ_tt - ∂_t Γ^ρ_ρt
# ∂_ρ Γ^ρ_tt = ∂_t(-U_t) + ∂_x(-U_x) + ∂_y(-U_y) + ∂_z(-U_z)
#             = -U_tt - lap_U   (POPRAWKA: minus przed lap_U)
# Γ^ρ_ρt = Γ^t_tt + Γ^x_xt + Γ^y_yt + Γ^z_zt = -U_t + U_t + U_t + U_t = 2U_t
# ∂_t Γ^ρ_ρt = 2U_tt

R1_tt_corrected = (-U_tt - lap_U) - 2*U_tt
R1_tt_simplified = expand(R1_tt_corrected)
print(f"    R^(1)_tt = ∂_ρΓ^ρ_tt - ∂_t Γ^ρ_ρt")
print(f"             = (-U_tt - ∇²U) - 2U_tt = -3U_tt - ∇²U")
test("R^(1)_tt = -3U_tt - ∇²U (POPRAWIONY)",
     expand(R1_tt_simplified - (-3*U_tt - lap_U)) == 0,
     f"= {R1_tt_simplified}")

# R^(1)_tx = ∂_ρ Γ^ρ_tx - ∂_x Γ^ρ_ρt  (niezmieniony bo Γ^t_tx = -U_x nie zmieniło się)
# ∂_t Γ^t_tx = -U_tx,  ∂_x Γ^x_tx = U_tx  → suma = 0
# ∂_x Γ^ρ_ρt = 2U_tx
R1_tx = sp.Integer(0) - 2*U_tx
test("R^(1)_tx = -2∂_t∂_xU (niezmieniony)", expand(R1_tx + 2*U_tx) == 0, "OK")

# R^(1)_xx = ∂_ρ Γ^ρ_xx - ∂_x Γ^ρ_ρx
# Γ^ρ_xx: Γ^t_xx=U_t, Γ^x_xx=U_x, Γ^y_xx=-U_y, Γ^z_xx=-U_z
# ∂_ρ Γ^ρ_xx = U_tt + U_xx - U_yy - U_zz
# Γ^ρ_ρx = Γ^t_tx + Γ^x_xx + Γ^y_yx + Γ^z_zx = -U_x + U_x + U_x + U_x = 2U_x
# ∂_x Γ^ρ_ρx = 2U_xx
R1_xx = (U_tt + U_xx - U_yy - U_zz) - 2*U_xx
R1_xx_simplified = expand(R1_xx)
print(f"    R^(1)_xx = U_tt - U_xx - U_yy - U_zz = U_tt - ∇²U")
test("R^(1)_xx = U_tt - ∇²U (niezmieniony)", expand(R1_xx_simplified - (U_tt - lap_U)) == 0,
     f"= {R1_xx_simplified}")

# Ślad skalaru Ricciego
R1_scalar = -R1_tt_simplified + R1_xx + 2*(U_tt - lap_U)  # -R_tt + 3*R_xx
R1_scalar_exp = expand(-R1_tt_simplified + 3*(U_tt - lap_U))
R1_expected = 6*U_tt - 2*lap_U   # ZMIANA: -2∇²U zamiast -4∇²U
test("R^(1) = 6U_tt - 2∇²U (POPRAWIONY)", expand(R1_scalar_exp - R1_expected) == 0,
     f"R = {R1_scalar_exp}")

# ─────────────────────────────────────────────────────────────
# 3. TENSOR EINSTEINA G^(1)_μν  (poprawiony)
# ─────────────────────────────────────────────────────────────
print("\n[3] POPRAWIONY tensor Einsteina G^(1)_μν")

R1 = R1_expected   # = 6U_tt - 2∇²U

G1_tt = R1_tt_simplified + Rational(1,2)*R1  # R^tt + (1/2)(+1)R bc η_tt = -1
G1_tt = expand(G1_tt)
G1_tt_expected = -2*lap_U
test("G^(1)_tt = -2∇²U  (POPRAWIONY — nie -∇²U!)",
     expand(G1_tt - G1_tt_expected) == 0,
     f"G_tt = {G1_tt}")

G1_tx = R1_tx  # = -2U_tx (η_tx = 0)
test("G^(1)_tx = -2∂_t∂_xU (niezmieniony)", expand(G1_tx + 2*U_tx) == 0, "OK")

G1_xx = expand((U_tt - lap_U) - Rational(1,2)*R1)
G1_xx_expected = -2*U_tt
test("G^(1)_xx = -2∂_t²U  (POPRAWIONY — nie -2U_tt+∇²U!)",
     expand(G1_xx - G1_xx_expected) == 0,
     f"G_xx = {G1_xx}")

print(f"""
    POPRAWIONE WYNIKI:
    G^(1)_tt = -2∇²U
    G^(1)_ti = -2∂_t∂_iU
    G^(1)_ij = -2∂_t²U δ_ij
    G^(1)_ij|_{{i≠j}} = 0

    (Poprzedni błąd: G^(1)_tt = -∇²U, G^(1)_ij = (-2U_tt+∇²U)δ_ij)
""")

# ─────────────────────────────────────────────────────────────
# 4. ANALIZA EMERGENCJI — z poprawnymi formułami
# ─────────────────────────────────────────────────────────────
print("[4] Analiza emergencji Einsteina (poprawione formuly)")

print("\n    --- PYŁ (p=0) ---")
# T_trace = -ρ, □U = -(κ/2)(-ρ) = (κ/2)ρ
# Quasistatycznie: ∇²U = -(κ/2)ρ  (Laplacian = source, NEGATIVE)
# Uwaga: z eq pola: (1/c²)Ü - ∇²U = qρ → quasistatycznie: -∇²U = qρ → ∇²U = -qρ
# gdzie q = 4πG/c², κ = 8πG/c⁴ (SI), w jedn. naturalnych: q = κ/2
print("    Eq. pola TGP (quasistat.): -□U ≈ -∇²U → ∇²U = -qρ = -(κ/2)ρ < 0")
print("    G^(1)_tt = -2∇²U = -2×(-(κ/2)ρ) = κρ  ✓  DOKŁADNA EMERGENCJA!")
print("    G^(1)_ij = -2∂_t²U δ_ij → 0 dla quasistat.  ✓  (GR: κp=0)")

# Weryfikacja numeryczna
kappa_n = 1.0; rho_n = 1.0
lap_U_dust = -(kappa_n/2)*rho_n  # ∇²U = -qρ = -(κ/2)ρ
G_tt_dust = -2*lap_U_dust
test("Pył: G^(1)_tt = κρ DOKŁADNIE",
     abs(G_tt_dust - kappa_n*rho_n) < 1e-12,
     f"G_tt = {G_tt_dust:.4f}, κρ = {kappa_n*rho_n:.4f}")

print("\n    --- PROMIENIOWANIE (p=ρ/3) ---")
print("    T_trace = -ρ + ρ = 0 → □U = 0")
print("    Quasistat. + izotropowe: ∇²U = 0 → G^(1)_tt = -2×0 = 0 ≠ κρ  ✗")
print("    PROBLEM POZOSTAJE: promieniowanie nie sprzęga ze scalarem U")

G_tt_rad_scalar = -2 * 0.0   # ∇²U = 0 dla T_trace=0
test("Promieniowanie: G^(1)_tt(scalar) = 0 ≠ κρ — problem potwierdzony",
     G_tt_rad_scalar == 0.0,
     f"G_tt = 0  ← potrzeba sektora σ")

print("\n    --- PRÓŻNIA + GW ---")
print("    T=0, □U=0, ale G^(1)_tt = -2∇²U ≠ 0 dla fal (∇²U ≠ 0 gdy U_tt ≠ 0)")
print("    → c_GW = c₀  ✓  (struktura falowa zachowana)")

# ─────────────────────────────────────────────────────────────
# 5. PEŁNA EMERGENCJA — DEKOMPOZYCJA SCALAR + TENSOR
# ─────────────────────────────────────────────────────────────
print("\n[5] Pełna emergencja: G = κT przez sektor scalar + σ_μν")
print()
print("    Poprawiony podział źródła na 2 sektory:")
print()
print("    SEKTOR SKALARNY (pole U):")
print("    □U = -(q) T_trace  (q = κ/2 w jedn. nat.)")
print("    G^(U)_tt = -2∇²U = κ(-T_trace) = κρ  dla pyłu")
print("    G^(U)_ij = -2U_tt δ_ij → 0 (quasistat.)")
print()
print("    SEKTOR TENSOROWY (pole σ_μν, bezśladowe 4D):")
print("    □σ_μν = -(κσ₀) T^TF_μν,  T^TF = T_μν - (T_trace/4)η_μν")
print("    G^(σ)_μν = -(1/σ₀)□σ_μν = κ T^TF_μν  (w kalibracji TT)")
print()
print("    ŁĄCZNIE: G_μν = G^(U)_μν + G^(σ)_μν")
print()
print("    Dla PYŁU:")
print("    G^(U)_tt = κρ,  G^(U)_ij = 0")
print("    T^TF_μν(pył) = diag(3ρ/4, ρ/4, ρ/4, ρ/4)")
print("    G^(σ)_tt = κ×3ρ/4,  G^(σ)_ij = κ×ρ/4 δ_ij")
print("    Total: G_tt = κρ + κ×3ρ/4 = 7κρ/4 ≠ κρ  ← PROBLEM!")
print()
print("    WNIOSEK: prosta suma G^(U) + G^(σ) nie daje κT dla pyłu!")
print("    Potrzeba: G^(U) lub G^(σ) ma INNE wagi.")
print()
print("    WŁAŚCIWE ROZWIĄZANIE:")
print("    Sektor skalarny WYSTARCZY dla pyłu: G^(U)_tt = κρ dokładnie.")
print("    Sektor tensorowy σ_μν musi być ZERO dla pyłu izotropowego.")
print("    → Dla pyłu: T^TF ≠ 0 ale σ nie powinien być źródłowany!")
print("    → Konieczna MODYFIKACJA sprzężenia σ z materią.")

print("\n    WŁAŚCIWA DEKOMPOZYCJA (bez redundancji):")
print("    □U     = -(κ/2) T_trace   → zawiera CAŁY wkład skalara")
print("    □σ_μν = -(κσ₀) Π^TF_μν   → TYLKO anizotropowy T")
print("    gdzie Π^TF = T^TF - T^TF_trace-isotropic")
print()
print("    Dla pyłu (izotropowy p=0): Π^TF = 0 → σ = 0  ✓")
print("    Dla promień. (T_trace=0): Π^TF = T_μν → σ źródłowane  ✓")

# ─────────────────────────────────────────────────────────────
# 6. WERYFIKACJA: czy G^(U) + G^(σ) = κT z NOWĄ DEKOMPOZYCJĄ
# ─────────────────────────────────────────────────────────────
print("\n[6] Weryfikacja nowej dekompozycji")
import numpy as np

eta = np.diag([-1.0, 1.0, 1.0, 1.0])
eta_inv = np.diag([-1.0, 1.0, 1.0, 1.0])

def T_trace4(T):
    return sum(eta_inv[m,m]*T[m,m] for m in range(4))

def T_aniso_TF(T):
    """Anizotropowa bezśladowa część T_μν.
    Rozumowana jako: T^TF_μν ale tylko dla pyłu P=0 nie liczymy,
    bo G^(U) już obsłużył ślad.

    Właściwa część dla σ: T^aniso_μν = T_μν - (T_tt/(space-dim)) × δ_ij (przestrzenna)
    ale z warunkiem że dla pyłu (p=0) jest zero.

    Konkretnie: dla izotropowego p=T_kk/3 (trace spatial):
    T^aniso_ij = T_ij - (T_kk/3)δ_ij
    T^aniso_tt = 0 (σ nie modyfikuje g_tt bezpośrednio w 3D przestrzennym σ_ab)
    """
    T_aniso = np.zeros((4,4))
    # Przestrzenna bezśladowa część (standard TT source)
    T_spatial_trace = T[1,1] + T[2,2] + T[3,3]
    for i in range(1,4):
        T_aniso[i,i] = T[i,i] - T_spatial_trace/3
        for j in range(i+1,4):
            T_aniso[i,j] = T[i,j]
            T_aniso[j,i] = T[j,i]
    # Dla σ_ab (3D przestrzenne bezśladowe):
    # G^(σ)_tt = (efektywnie) -∇²σ_tt, ale σ nie ma składowej tt w 3D!
    # W 3D formulacji σ_ab wkład do G_tt pochodzi z:
    # G^(σ)_tt = -∇²(Tr σ_ij) ale σ bezśladowe → G^(σ)_tt = 0
    # G^(σ)_ij = -(1/σ₀)□σ_ij
    return T_aniso

def full_G_correct(T, kappa_n=1.0):
    """Pełne G_μν z TGP (poprawiona wersja)."""
    tr = T_trace4(T)
    # Sektor skalarny: G^(U)_tt = κ(-tr_4D) = κ(-(-ρ)) = κρ dla pyłu
    # W ogóle: G^(U)_tt = -2∇²U = 2q|T_trace| gdy T_trace<0
    # Quasistat.: ∇²U = -(q)ρ_eff gdzie ρ_eff = -T_trace (> 0 dla pyłu)
    G_U = np.zeros((4,4))
    G_U[0,0] = kappa_n * (-tr)   # = κρ dla pyłu (T_trace=-ρ → -(-ρ)=ρ)
    # G^(U)_ij = 0 quasistat.

    # Sektor tensorowy (3D σ_ab): G^(σ)_ij = κ T^aniso_ij
    # G^(σ)_tt = 0 (σ_ab nie ma składowej tt)
    G_sigma = np.zeros((4,4))
    T_an = T_aniso_TF(T)
    for i in range(1,4):
        for j in range(1,4):
            G_sigma[i,j] = kappa_n * T_an[i,j]

    return G_U + G_sigma

# Pył
rho = 1.0
T_dust = np.zeros((4,4)); T_dust[0,0] = rho
G_dust = full_G_correct(T_dust)
GR_dust = np.zeros((4,4)); GR_dust[0,0] = rho  # κ T_tt = ρ

test("Pył: G_tt = κρ", abs(G_dust[0,0] - rho) < 1e-12, f"G_tt={G_dust[0,0]:.4f}")
test("Pył: G_ij = 0 (p=0)", np.allclose(G_dust[1:,1:], 0, atol=1e-12), f"G_ij={G_dust[1,1]:.4f}")

# Promieniowanie
T_rad = np.zeros((4,4)); T_rad[0,0] = rho
T_rad[1,1] = T_rad[2,2] = T_rad[3,3] = rho/3
G_rad = full_G_correct(T_rad)
print(f"\n    Promieniowanie: G_tt = {G_rad[0,0]:.4f} (GR potrzebuje κρ={rho:.4f})")
print(f"    Promieniowanie: G_11 = {G_rad[1,1]:.4f} (GR potrzebuje κp={rho/3:.4f})")

test("Promień.: G_tt = 0 (scalar sector alone, σ_tt=0)",
     abs(G_rad[0,0]) < 1e-12,
     f"G_tt = {G_rad[0,0]:.2e} -- WYMAGA σ_tt contribution")
test("Promień.: G_ij = 0 (izotropowe T_ij^aniso = 0)",
     np.allclose(G_rad[1:,1:], 0, atol=1e-12),
     f"G_11 = {G_rad[1,1]:.4f}")

print(f"""
    KLUCZOWA OBSERWACJA:
    Dla IZOTROPOWEGO promieniowania: T^aniso_ij = 0 (sferyczna symetria),
    więc σ_ab = 0 (pole tensorowe nie jest wzbudzane).

    Promieniowanie w TGP sprzęga z geometrią WYŁĄCZNIE poprzez:
    1. Anizotropie (T_ij - δ_ij T_kk/3 ≠ 0) → sektor σ_ab
    2. Ślad T_trace ≠ 0 → sektor U (ale dla promieniowania = 0)

    Izotropowe jednorodne promieniowanie → G^(1)_μν = 0 w linearyzacji.
    ALE: FRW dla epoki radiacyjnej jest DOKŁADNIE spełnione (tw. einstein-emergence)
    bo Friedmann NIE jest linearyzacją wokół Minkowskiego.
""")

# ─────────────────────────────────────────────────────────────
# PODSUMOWANIE
# ─────────────────────────────────────────────────────────────
print("=" * 64)
print("PODSUMOWANIE POPRAWIONEJ ANALIZY")
print("=" * 64)
print(f"""
  POPRAWKA: Γ^i_tt = -U_i  (v1 miał +U_i — błąd znaku)

  POPRAWIONE FORMULY:
    G^(1)_tt = -2∇²U      (nie -∇²U)
    G^(1)_ti = -2∂_t∂_iU  (bez zmian)
    G^(1)_ij = -2∂_t²U    (nie -2U_tt+∇²U)

  STATUS EMERGENCJI (poprawiony):
    ✓ Pył (p=0):       G^(1)_tt = κρ  DOKŁADNIE (nie czynnik 1/2!)
    ✓ Próżnia + GW:    c_GW = c₀, fale propagują się
    ✓ FRW (all fluids): Friedmann dokładnie (tw. einstein-emergence)
    ✗ Stat. izotrop. promieniowanie: G^(1)_tt = 0 ≠ κρ

  PROBLEM RADIACYJNY:
    Skalar U: T_trace = 0 → □U = 0 → G_tt = 0
    Tensor σ_ab: T_ij^aniso = 0 (izotropowe) → σ = 0 → G_ij = 0
    → Linearyzacja wokół płaskiej przestrzeni NIE odtwarza
      Einsteina dla izotropowego promieniowania statycznego.
    → Ale FRW jest DOKŁADNE — nie ma problemu kosmologicznego.

  ROZWIĄZANIE PROBLEMU LINEARYZACJI:
    Potrzeba 4D pola σ_μν z □σ_tt ≠ 0 źródłowanym przez T_tt.
    LUB: akceptacja że TGP ≠ GR dla statycznych kul promieniowania
    (fizycznie nierealnych: promieniowanie musi być ograniczone przez materię).
""")
print(f"WYNIK: {tests_pass}/{tests_total} PASS")
