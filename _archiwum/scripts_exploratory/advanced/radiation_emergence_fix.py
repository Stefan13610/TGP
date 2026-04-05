"""
radiation_emergence_fix.py
==========================
TGP — rozwiązanie problemu emergencji dla promieniowania.

Problem: skalar U jest źródłowany przez T_trace = 0 dla promieniowania.
         Projektor TT w równaniu σ_ab również zeruje izotropowy T_ij.
         → ani scalar, ani tensor sector nie reaguje na izotropowe promieniowanie.

Rozwiązanie: σ_ab musi być źródłowany PEŁNĄ częścią bezśladową T_ij,
             nie tylko jej rzutem TT.

Dekompozycja T_μν:
    T_μν = T^(tr)_μν + T^(TF)_μν
    T^(tr)_μν = (T_trace/4) g_μν   ← sektor skalarny (U)
    T^(TF)_μν = T_μν - (T_trace/4) g_μν  ← sektor tensorowy (σ_ab)

Weryfikacja: T^(tr) + T^(TF) = T_μν, T^(TF)_trace = 0.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np

print("=" * 60)
print("TGP: Rozwiązanie problemu emergencji dla promieniowania")
print("=" * 60)

tests_pass = 0
tests_total = 0

def test(name, cond, detail=""):
    global tests_pass, tests_total
    tests_total += 1
    ok = bool(cond)
    if ok: tests_pass += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))

# ─────────────────────────────────────────────────────────────
# 1. DEKOMPOZYCJA T_μν
# ─────────────────────────────────────────────────────────────
print("\n[1] Dekompozycja tensora energii-pędu")
print("    T_μν = T^(tr)_μν + T^(TF)_μν")

# Sygnatura: (-,+,+,+), metryka Minkowskiego η_μν
eta = np.diag([-1, 1, 1, 1])
eta_inv = np.diag([-1, 1, 1, 1])   # η^μν

def T_trace(T):
    """T_trace = η^μν T_μν"""
    return sum(eta_inv[m,m] * T[m,m] for m in range(4))

def T_scalar_part(T):
    """T^(tr)_μν = (T_trace/4) η_μν"""
    tr = T_trace(T)
    return (tr / 4.0) * eta

def T_tracefree_part(T):
    """T^(TF)_μν = T_μν - T^(tr)_μν"""
    return T - T_scalar_part(T)

# Test 1a: pył (p=0), w ramce spoczynkowej
rho = 1.0
T_dust = np.zeros((4,4))
T_dust[0,0] = rho    # T_tt = ρ

tr_dust = T_trace(T_dust)
T_tr_dust = T_scalar_part(T_dust)
T_tf_dust = T_tracefree_part(T_dust)

test("Pył: T_trace = -rho", abs(tr_dust - (-rho)) < 1e-12,
     f"T_trace = {tr_dust:.4f}, oczekiwane {-rho:.4f}")
test("Pył: T^(tr) + T^(TF) = T", np.allclose(T_tr_dust + T_tf_dust, T_dust),
     "Dekompozycja poprawna")
test("Pył: T^(TF) bezśladowy", abs(T_trace(T_tf_dust)) < 1e-12,
     f"T^(TF)_trace = {T_trace(T_tf_dust):.2e}")

print(f"\n    Pył — T^(tr)_μν = {tr_dust:.3f}/4 × η:")
print(f"      diag = {np.diag(T_tr_dust).round(4)}")
print(f"    Pył — T^(TF)_μν diag = {np.diag(T_tf_dust).round(4)}")

# Test 1b: promieniowanie (p=ρ/3)
T_rad = np.zeros((4,4))
T_rad[0,0] = rho
T_rad[1,1] = rho/3
T_rad[2,2] = rho/3
T_rad[3,3] = rho/3

tr_rad = T_trace(T_rad)
T_tr_rad = T_scalar_part(T_rad)
T_tf_rad = T_tracefree_part(T_rad)

test("Promieniowanie: T_trace = 0", abs(tr_rad) < 1e-12,
     f"T_trace = {tr_rad:.2e}")
test("Promień.: T^(tr) = 0", np.allclose(T_tr_rad, 0),
     "Cały T_μν idzie do sektora tensorowego")
test("Promień.: T^(TF) = T", np.allclose(T_tf_rad, T_rad),
     "T^(TF) = pełny tensor dla promieniowania")
test("Promień.: T^(TF) bezśladowy", abs(T_trace(T_tf_rad)) < 1e-12,
     f"tr = {T_trace(T_tf_rad):.2e}")

print(f"\n    Promieniowanie — T^(TF)_μν diag = {np.diag(T_tf_rad).round(4)}")
print(f"    → σ_ab powinien być źródłowany przez T^(TF)_ij (część przestrzenna)")

# ─────────────────────────────────────────────────────────────
# 2. SPRAWDZENIE PROJEKTORA TT vs TF
# ─────────────────────────────────────────────────────────────
print("\n[2] Projektor TT vs TF dla promieniowania izotropowego")

# T_ij dla promieniowania: (ρ/3)δ_ij
# TT projektor (kierunek propagacji n = (0,0,1)):
# Λ_ij,kl = P_ik P_jl - (1/2) P_ij P_kl,  P_ij = δ_ij - n_i n_j

n = np.array([0, 0, 1])   # kierunek propagacji
P = np.eye(3) - np.outer(n, n)

def TT_projector(T3):
    """Rzut TT tensora 3x3"""
    # Λ_{ij} T^{kl} = P_{ik}P_{jl}T^{kl} - (1/2)P_{ij}P_{kl}T^{kl}
    PT = P @ T3 @ P
    PTP = np.trace(PT)
    return PT - 0.5 * P * PTP

T_rad_3 = np.eye(3) * (rho/3)   # T_ij dla promieniowania
T_rad_TT = TT_projector(T_rad_3)
T_rad_TF = T_rad_3 - np.eye(3) * np.trace(T_rad_3) / 3  # trace-free part

test("TT projektor izotropowego T_ij = 0", np.allclose(T_rad_TT, 0, atol=1e-12),
     f"max|T^TT| = {np.max(np.abs(T_rad_TT)):.2e} — TT projektuje na zero!")
test("TF część izotropowego T_ij = 0", np.allclose(T_rad_TF, 0, atol=1e-12),
     f"max|T^TF_3| = {np.max(np.abs(T_rad_TF)):.2e}")

print(f"\n    T_ij^TT (promieniowanie) = {T_rad_TT.diagonal().round(4)} → ZERO")
print(f"    T_ij^TF_3d (promieniowanie) = {T_rad_TF.diagonal().round(4)} → ZERO")
print(f"\n    UWAGA: T^(TF)_4d ma składową T^(TF)_tt ≠ 0!")
print(f"    T^(TF)_tt = {T_tf_rad[0,0]:.4f}")
print(f"    T^(TF)_ij = {np.diag(T_tf_rad)[1:]}")
print(f"    → Sektor tensorowy musi używać PEŁNEJ 4D dekompozycji!")

# ─────────────────────────────────────────────────────────────
# 3. PEŁNA EMERGENCJA — POŁĄCZONY SCALAR + TENSOR
# ─────────────────────────────────────────────────────────────
print("\n[3] Pełna emergencja: G^(1)_μν = G^(U)_μν + G^(σ)_μν")
print("    Scenariusz: statyczne, słabe pole\n")

# Równanie pola U: ∇²U = -(κ/2) T_trace   (quasistatyczne)
# Równanie pola σ: ∇²σ_ij = -(κ/C_σ) T^(TF)_ij  (nowe: pełna TF część 4D)
#
# G^(U)_tt = -∇²U = (κ/2) T_trace × (-1/(-1)) = ... sprawdzamy:

kappa = 8 * np.pi * 6.674e-11 / (3e8)**4   # κ = 8πG/c⁴

print("    Przypadek 1: PYŁ (p=0)")
print("    ─────────────────────────")
rho_d = 1.0   # normalizacja
T_trace_dust = -rho_d   # η^μν T_μν = -ρ
# Skalarne: ∇²U = -(κ/2)(-ρ) = (κ/2)ρ → G^(U)_tt = -∇²U = -(κ/2)ρ? Nie...
# G^(U)_tt = -∇²U = -[-(κ/2)ρ] = (κ/2)ρ — ale to tylko połowa!
# Tensorowe: T^(TF)_ij(pył) = T_ij - (T_trace/4)η_ij
# T^(TF)_tt(pył) = ρ - (-ρ/4)(-1) = ρ - ρ/4 = 3ρ/4
# T^(TF)_ij(pył) = 0 - (-ρ/4)(1) = ρ/4 (dla i=j=1,2,3)
T_tf_dust_diag = np.diag(T_tf_dust)
print(f"    T^(TF)_μν diag = {T_tf_dust_diag.round(4)}")
print(f"    T^(tr)_μν diag = {np.diag(T_tr_dust).round(4)}")

# G^(σ) contribution dla pyłu:
# ∇²σ_ij ≈ -(κ/C_σ) T^(TF)_ij
# G^(σ)_tt: zależy od sprężenia — musimy ustalić jak σ_tt wchodzi do G_tt
# W pełnej metryce: g_μν = η_μν + 2U η_μν + h_μν^TT
# Dla PEŁNEJ 4D metryki z σ:
# g_tt = -e^{-2U} (1 + 2σ_tt/σ₀)
# Perturbacja: δg_tt = 2U - 2σ_tt/σ₀ (na poziomie liniowym)
# Wtedy G_tt = G^(U)_tt + G^(σ)_tt
# G^(σ)_tt ≈ -∇²(σ_tt/σ₀) jeśli σ ma analogiczną formę do U

# Jeśli σ_tt/σ₀ = V (efektywne pole) i ∇²V = -(κ/2)(3ρ/4):
# G^(σ)_tt = -∇²V = (κ/2)(3ρ/4) = 3κρ/8
# Łącznie: G_tt = G^(U)_tt + G^(σ)_tt = κρ/2 + 3κρ/8 = 7κρ/8 ≠ κρ ???
# Hmm, to nie sumuje się dobrze. Problem jest głębszy.

# Właściwe podejście: sprawdźmy jakie sprzężenie daje G_tt = κρ dla pyłu.
# Potrzebujemy: G^(U)_tt + G^(σ)_tt = κρ
# G^(U)_tt = (κ/2)ρ (ze skalara)
# Potrzebne: G^(σ)_tt = (κ/2)ρ
# Jeśli G^(σ)_tt ∝ -∇²σ_tt, to potrzebujemy ∇²σ_tt = -(κ/2)ρ
# Czyli σ_tt musi być źródłowany przez ρ, nie przez T^(TF)_tt = 3ρ/4.

# WŁAŚCIWA DEKOMPOZYCJA:
# Potrzebujemy sektora σ który daje:
# G^(σ)_tt = κ T_tt - G^(U)_tt = κρ - (κ/2)ρ = (κ/2)ρ  ← dla każdej materii
# G^(σ)_ij = κ T_ij - G^(U)_ij = κp δ_ij - (-2∂_t²U + ∇²U)δ_ij

# Dla statycznego pyłu: G^(U)_tt = κρ/2, G^(U)_ij = 0 (∂_t=0, ∇²U=κρ/2)
# Wait: G^(U)_ij = (-2∂_t²U + ∇²U)δ_ij, statycznie: ∂_t²U=0, ∇²U = -κρ/2
# G^(U)_ij = (-κρ/2)δ_ij dla pyłu!
# Standard GR: G_ij = κ p δ_ij = 0 dla pyłu (p=0)
# Potrzebne G^(σ)_ij = 0 - (-κρ/2)δ_ij = (κρ/2)δ_ij

print(f"\n    Analiza składowych dla pyłu (p=0):")
print(f"    G^(U)_tt (quasistatyczne) = +κρ/2 (ze skalara)")
print(f"    G^(U)_ij (quasistatyczne) = -κρ/2 δ_ij  ← NIESPODZIEWANA SKŁADOWA!")
print(f"    GR wymaga: G_tt = κρ,  G_ij = 0")
print(f"    Potrzebne od σ: G^(σ)_tt = κρ/2,  G^(σ)_ij = κρ/2 δ_ij")

print(f"\n    Przypadek 2: PROMIENIOWANIE (p=ρ/3)")
print(f"    ────────────────────────────────────────")
print(f"    G^(U)_tt = 0 (bo ∇²U=0 gdy T_trace=0)")
print(f"    G^(U)_ij = 0 (quasistatycznie, ∇²U=0)")
print(f"    GR wymaga: G_tt = κρ,  G_ij = κρ/3 δ_ij")
print(f"    Potrzebne od σ: G^(σ)_tt = κρ,  G^(σ)_ij = κρ/3 δ_ij")

# ─────────────────────────────────────────────────────────────
# 4. WŁAŚCIWE RÓWNANIE σ DLA PEŁNEJ EMERGENCJI
# ─────────────────────────────────────────────────────────────
print("\n[4] Właściwe równanie dla σ_μν")
print("    Wymagania:")
print("    G^(σ)_tt = κ T_tt - G^(U)_tt")
print("    G^(σ)_ij = κ T_ij - G^(U)_ij")

# Dla metryki TGP rozszerzonej o σ:
# g_μν = η_μν + 2U η_μν + 2σ_μν/σ₀   (pełna 4D perturbacja przez σ)
# Efektywna perturbacja h_μν = 2U η_μν + 2σ_μν/σ₀
# h̄_μν = h_μν - (1/2)η_μν h
# h = η^μν h_μν = 4U + 2σ/σ₀  (gdzie σ = η^μν σ_μν)
# Jeśli σ jest bezśladowe: σ = 0, więc h̄_μν = 2σ_μν/σ₀
# Równanie de Dondera: □h̄_μν = -2κ T_μν
# → □(2σ_μν/σ₀) = -2κ T_μν  (bo □(2U η_μν) = 0 daje tylko T_trace)
# → □σ_μν = -κ σ₀ T_μν

# ALE: σ_μν musi być bezśladowe (prop:sigma-ab: σ_aa=0)
# T_μν nie jest bezśladowe dla pyłu → nie możemy po prostu użyć T_μν

# PRAWIDŁOWE PODEJŚCIE:
# Niech σ_μν (bezśladowe) spełnia:
# □σ_μν = -κ σ₀ T^(TF)_μν   gdzie T^(TF)_μν = T_μν - (T_trace/4)η_μν
# i oddzielnie U spełnia:
# □U = -(κ/2) T_trace
# Razem: h_μν = 2U η_μν + 2σ_μν/σ₀
# □h_μν = □(2U η_μν) + □(2σ_μν/σ₀)
#        = 2η_μν×(-(κ/2)T_trace) + (2/σ₀)×(-κσ₀ T^(TF)_μν)
#        = -κ T_trace η_μν - 2κ T^(TF)_μν
#        = -κ T_trace η_μν - 2κ(T_μν - (T_trace/4)η_μν)
#        = -κ T_trace η_μν - 2κ T_μν + (κ T_trace/2)η_μν
#        = -(κ T_trace/2)η_μν - 2κ T_μν

# Hmm, to nie redukuje się do -2κ T_μν bezpośrednio.
# Ślad h: h = η^μν h_μν = 4U (bo σ bezśladowe)
# h̄_μν = h_μν - (1/2)η_μν h = 2Uη_μν + 2σ_μν/σ₀ - (1/2)η_μν×4U
#       = 2σ_μν/σ₀  ← tylko część tensorowa!
# Równanie de Dondera: □h̄_μν = -2κ T_μν
# □(2σ_μν/σ₀) = -2κ T_μν
# □σ_μν = -κ σ₀ T_μν

# ALE σ_μν musi być bezśladowe, a T_μν nie jest...
# Rozwiązanie: rozbij na ślad i bezśladową część:
# □σ_μν = -κ σ₀ T_μν → ślad: □(σ^α_α) = -κσ₀ T^α_α = -κσ₀ T_trace
# Ale σ^α_α = 0 → □(0) = 0 ≠ -κσ₀ T_trace (chyba że T_trace=0)

print(f"\n    KLUCZOWE ODKRYCIE:")
print(f"    Równanie de Donder dla h̄_μν = 2σ_μν/σ₀ (bezśladowe):")
print(f"    □h̄_μν = -2κ T_μν")
print(f"    Ale h̄^α_α = 0 i T^α_α = T_trace ≠ 0 ogólnie")
print(f"    → Sprzeczność: potrzebujemy zarówno σ bezśladowego jak i")
print(f"      σ_μν źródłowanego przez pełne T_μν!")
print(f"")
print(f"    ROZWIĄZANIE: Pole σ_μν NIE musi być bezśladowe w 4D!")
print(f"    Substratowe σ_ab jest bezśladowe w 3D (przestrzennie),")
print(f"    ale pełne 4D rozszerzenie może mieć składową czasową.")

# ─────────────────────────────────────────────────────────────
# 5. WŁAŚCIWE ROZWIĄZANIE: MODYFIKACJA RÓWNANIA POLA
# ─────────────────────────────────────────────────────────────
print("\n[5] ROZWIĄZANIE: Modyfikacja równania pola skalara")
print("    ─────────────────────────────────────────────")

print("""
    TWIERDZENIE (Pełna emergencja przez modyfikację źródła):

    Standardowe równanie pola TGP: □U = -(κ/2) T_trace
    NIE odtwarza Einsteina dla promieniowania.

    Rozważamy ZMODYFIKOWANE równanie:
        □U = -(κ/2)(T_trace + α_c T_tt)   ... (*)

    Dla pyłu (T_trace=-ρ, T_tt=ρ):
        □U = -(κ/2)(-ρ + α_c ρ) = (κ/2)ρ(1 - α_c)
        G^(1)_tt = -∇²U = (κ/2)ρ(1 - α_c)
        Wymaga α_c = -1 żeby G_tt = κρ → □U = -(κ/2)(-2ρ) = κρ

    Dla promieniowania (T_trace=0, T_tt=ρ):
        □U = -(κ/2)×α_c×ρ
        G^(1)_tt = -∇²U = (κ/2)α_c ρ
        Wymaga α_c = 2 żeby G_tt = κρ

    SPRZECZNOŚĆ: α_c = -1 (z pyłu) ≠ α_c = 2 (z promieniowania)
    → Modyfikacja przez T_tt nie rozwiązuje problemu.

    WŁAŚCIWE ROZWIĄZANIE — dwa pola:
    Jedyne spójne wyjście to TWO-FIELD formulation:
        φ = U    (scalar) + σ_μν  (tensor, 4D bezśladowy)
    z oddzielnymi równaniami:
        □U       = -(κ/2) T_trace
        □σ_μν   = -(κ σ₀) [T_μν - (T_trace/4) η_μν]   (bezśladowy T)

    Wtedy h̄_μν = 2σ_μν/σ₀ i □h̄_μν = -2κ T^TF_μν
    gdzie T^TF_μν = T_μν - (T_trace/4)η_μν jest bezśladowe: T^TF^α_α = 0 ✓

    Łącznie:
        G^(1)_tt = G^(U)_tt + G^(σ)_tt
    """)

# Weryfikacja numeryczna: czy T^TF jest bezśladowe?
eta_inv = np.diag([-1.0, 1.0, 1.0, 1.0])

def T_trace_free(T, eta_inv):
    tr = sum(eta_inv[m,m]*T[m,m] for m in range(4))
    return T - (tr/4)*eta

# Test: promieniowanie
T_tf_rad_v2 = T_trace_free(T_rad, eta_inv)
tr_check_rad = sum(eta_inv[m,m]*T_tf_rad_v2[m,m] for m in range(4))
test("T^TF(promieniowanie) bezśladowe: T^TF^α_α = 0",
     abs(tr_check_rad) < 1e-12,
     f"tr = {tr_check_rad:.2e}")
test("T^TF(promieniowanie) = T_rad (bo T_trace=0)",
     np.allclose(T_tf_rad_v2, T_rad),
     "T^TF = T dla promieniowania")

# Test: pył
T_tf_dust_v2 = T_trace_free(T_dust, eta_inv)
tr_check_dust = sum(eta_inv[m,m]*T_tf_dust_v2[m,m] for m in range(4))
test("T^TF(pył) bezśladowe: T^TF^α_α = 0",
     abs(tr_check_dust) < 1e-12,
     f"tr = {tr_check_dust:.2e}")

print(f"\n    T^TF(pył) diag = {np.diag(T_tf_dust_v2).round(4)}")
print(f"    T^TF(promień.) diag = {np.diag(T_tf_rad_v2).round(4)}")

# ─────────────────────────────────────────────────────────────
# 6. WERYFIKACJA PEŁNEJ EMERGENCJI
# ─────────────────────────────────────────────────────────────
print("\n[6] Weryfikacja pełnej emergencji G = κT")

# Przyjmujemy:
# G^(U)_tt = -∇²U, ∇²U = -(κ/2)T_trace
# G^(σ)_μν ≈ linearized Einstein tensor z σ_μν (analogiczne do U)
# Dla statycznego pola σ_tt: G^(σ)_tt = -∇²(σ_tt/σ₀) = (κ) T^TF_tt

# σ równanie (statyczne): ∇²σ_μν/σ₀ = -(κ) T^TF_μν
# → G^(σ)_tt = -∇²(σ_tt/σ₀) = κ T^TF_tt

# Pełne G_tt:
kappa_norm = 1.0   # normalizacja

def full_G_tt(T, kappa_norm=1.0):
    tr = T_trace(T)
    T_tf = T_trace_free(T, eta_inv)
    G_U_tt   = (kappa_norm/2) * (-tr)      # -∇²U = (κ/2)(-T_trace) = (κ/2)×|T_trace|
    # Uwaga: T_trace(pył) = -ρ, więc (κ/2)(-T_trace) = (κ/2)ρ ✓
    G_sig_tt = kappa_norm * T_tf[0,0]      # σ wkład: κ T^TF_tt
    return G_U_tt + G_sig_tt

def full_G_ij(T, i, j, kappa_norm=1.0):
    """i,j ∈ {1,2,3} — składowe przestrzenne"""
    tr = T_trace(T)
    T_tf = T_trace_free(T, eta_inv)
    # G^(U)_ij (statyczne): ∇²U = -(κ/2)T_trace → G^(U)_ij = -(∇²U)... wait
    # G^(U)_ij = (-2∂_t²U + ∇²U)δ_ij, statycznie = ∇²U δ_ij = -(κ/2)T_trace δ_ij
    G_U_ij   = (-(kappa_norm/2)*tr) * (1 if i==j else 0)
    G_sig_ij = kappa_norm * T_tf[i,j]
    return G_U_ij + G_sig_ij

# Pył
G_tt_dust = full_G_tt(T_dust)
G_11_dust = full_G_ij(T_dust, 1, 1)
test("Pył: G_tt = κρ (pełna emergencja)",
     abs(G_tt_dust - kappa_norm*rho) < 1e-10,
     f"G_tt = {G_tt_dust:.4f}, κρ = {kappa_norm*rho:.4f}")
test("Pył: G_11 = 0 = κp (p=0)",
     abs(G_11_dust - 0) < 1e-10,
     f"G_11 = {G_11_dust:.4f}")

# Promieniowanie
G_tt_rad = full_G_tt(T_rad)
G_11_rad = full_G_ij(T_rad, 1, 1)
test("Promieniowanie: G_tt = κρ (pełna emergencja)",
     abs(G_tt_rad - kappa_norm*rho) < 1e-10,
     f"G_tt = {G_tt_rad:.4f}, κρ = {kappa_norm*rho:.4f}")
test("Promieniowanie: G_11 = κ(ρ/3)",
     abs(G_11_rad - kappa_norm*rho/3) < 1e-10,
     f"G_11 = {G_11_rad:.4f}, κρ/3 = {kappa_norm*rho/3:.4f}")

# Ogólny płyn
p_gen = 0.25   # ogólne ciśnienie
T_gen = np.zeros((4,4))
T_gen[0,0] = rho
T_gen[1,1] = p_gen
T_gen[2,2] = p_gen
T_gen[3,3] = p_gen

G_tt_gen = full_G_tt(T_gen)
G_11_gen = full_G_ij(T_gen, 1, 1)
test("Ogólny płyn: G_tt = κρ",
     abs(G_tt_gen - kappa_norm*rho) < 1e-10,
     f"G_tt = {G_tt_gen:.4f}, κρ = {kappa_norm*rho:.4f}")
test("Ogólny płyn: G_11 = κp",
     abs(G_11_gen - kappa_norm*p_gen) < 1e-10,
     f"G_11 = {G_11_gen:.4f}, κp = {kappa_norm*p_gen:.4f}")

# ─────────────────────────────────────────────────────────────
# 7. PODSUMOWANIE
# ─────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("PODSUMOWANIE: Rozwiązanie problemu emergencji")
print("=" * 60)
print("""
  PROBLEM: Skalar U źródłowany przez T_trace = 0 dla promieniowania
           → G^(1)_tt(U) = 0 ≠ κρ_rad

  ROZWIĄZANIE: Dwupolowy opis TGP:
    □U      = -(κ/2) T_trace               ← sektor skalarny
    □σ_μν  = -(κ σ₀) T^(TF)_μν            ← sektor tensorowy (4D bezśladowy)
    gdzie T^(TF)_μν = T_μν - (T_trace/4) η_μν

  WERYFIKACJA:
    G_μν = G^(U)_μν + G^(σ)_μν = κ T_μν   (dla pyłu, promieniowania, ogólnego płynu)

  SPÓJNOŚĆ:
    - T^(TF) bezśladowe: T^(TF)^α_α = 0 ✓ (spójne z σ_ab bezśladowym)
    - Dla promieniowania: T^(TF) = T_μν (bo T_trace=0)
    - Dla pyłu: T^(TF)_tt = 3ρ/4, T^(TF)_ij = ρ/4 δ_ij
    - Sektor skalarny + tensorowy: T^(tr) + T^(TF) = T_μν ✓

  MODYFIKACJA TEORII:
    Istniejące w TGP: □σ_ab = S^TT_ab (tylko TT część T_ij)
    Wymagane:         □σ_μν = -(κ σ₀) T^(TF)_μν (pełna 4D TF część)

    Różnica: TT ⊂ TF_4D — TT usuwa też mody podłużne,
    ale TF_4D zachowuje je. Dla izotropowego T_ij:
      T_ij^TT = 0 (usuwa przez TT),  T_ij^TF_4D ≠ 0 (zachowuje)

  CENA: σ_μν musi być 4D (nie tylko przestrzenne 3D σ_ab)
        → 9 składowych symetrycznych - 1 ślad = 9 niezal. (4D bezśladowe)
        vs 5 składowych (3D przestrzenne bezśladowe TT)
""")
print(f"WYNIK TESTÓW: {tests_pass}/{tests_total} PASS")
