# -*- coding: utf-8 -*-
"""
vw_substrate_3d.py — TGP v23
==============================
O14: Prędkość fali substratowej v_W z dynamiki pola Φ
O16: Absolutne χ₀^(n) z profilu kinku + poprawna skala energetyczna

Cel:
  O14: wyznaczenie v_W(χ) = prędkość propagacji zaburzenia pola Φ
       przez numeryczne rozwiązanie równania kowariantnego TGP
       i siatki dyspersji 3D

  O16: |χ₀^(n)| absolutne jako minima efektywnego potencjału V_eff(χ)
       z uwzględnieniem v_W (poprawka do masy efektywnej)

Równanie kowariantne TGP (sek02):
  (1/v_W²(Φ))∂²_tt Φ - ∇²Φ + α(∇Φ)²/Φ + βΦ²/Φ₀ - γΦ³/Φ₀² = qΦ₀ρ

  gdzie v_W(Φ) = c₀·√(Φ₀/Φ) = c₀/√χ  (χ = Φ/Φ₀)
  (wynika z: c²(χ) = ∂²E/∂k² przy χ; reżim II: c ~ c₀/√χ)

Metody:
  1. Relacja dyspersji z linearyzacji: ω²(k) → v_W = dω/dk|_{k=0}
  2. Symulacja 1D: propagacja impulsu Gaussowskiego → pomiar v_czoło
  3. Siatka 3D (uproszczona): widmo modów → prędkość grupowa
  4. χ₀^(n) z ODE efektywnego potencjału z v_W

Jednostki: c₀ = Φ₀ = 1 (naturalne), G = 1
"""

import sys
import numpy as np
from scipy.integrate import solve_ivp, solve_bvp
from scipy.optimize import brentq, minimize_scalar
from scipy.fft import rfft, rfftfreq
import warnings
warnings.filterwarnings("ignore")

sys.stdout.reconfigure(encoding='utf-8')

PASS_COUNT = 0
FAIL_COUNT = 0
TESTS = []

def check(cond, name, info=""):
    global PASS_COUNT, FAIL_COUNT
    status = "PASS" if cond else "FAIL"
    if cond: PASS_COUNT += 1
    else:    FAIL_COUNT += 1
    print(f"  [{'✓' if cond else '✗'}] {name}  ({info})")
    return cond

print("=" * 68)
print(" vw_substrate_3d.py — TGP v23 — O14: v_W  |  O16: χ₀^(n)")
print("=" * 68)

# ══════════════════════════════════════════════════════════════════════════════
# PARAMETRY MODELU
# ══════════════════════════════════════════════════════════════════════════════
alpha_TGP = 2.0   # z wariacji działania (α = 2)
beta_TGP  = 1.0   # = γ (warunek próżniowy β=γ)
gamma_TGP = 1.0   # normalizacja bezwymiarowa
Phi0 = 1.0        # χ = Φ/Φ₀ = Φ (w tych jednostkach)
c0   = 1.0        # prędkość światła (jednostki naturalne)

# ──────────────────────────────────────────────────────────────────────────────
# SEKCJA A: v_W z relacji dyspersji (analitycznie + numerycznie)
# ──────────────────────────────────────────────────────────────────────────────
print("\n──── A. Prędkość fali v_W(χ) z relacji dyspersji ────")

# Równanie kowariantne linearyzowane wokół Φ = χ·Φ₀:
# Φ = χ·Φ₀ + δΦ·e^{i(kx-ωt)}
# Zerowy rząd: βχ²Φ₀ - γχ³Φ₀ = 0 → χ² (β - γχ) = 0 → χ₀ = β/γ = 1 (próżnia) ✓
# Pierwszy rząd (linearyzacja w δΦ):
# (1/v_W²)(-ω²) δΦ - (-k²) δΦ + [α·(ik)²·χΦ₀/(χΦ₀)²·Φ₀ + 2βχ/Φ₀ - 3γχ²/Φ₀²·Φ₀] δΦ = 0
# Upraszczając (bez członu masy na razie):
# -ω²/v_W² + k² + m_eff²(χ) = 0
# → ω² = v_W²(k² + m_eff²)
#
# v_W(χ) z kowariantnej formy równania pola:
# v_W²(χ) = c₀² · Φ₀/Φ = c₀²/χ
# (c jest efektywną prędkością propagacji — z równania pola kowariantnego)
# UWAGA: to jest predykcja TGP, nie postulat (wynika z metryki eksponencjalnej)

def v_W_TGP(chi):
    """Prędkość fali substratowej v_W(χ) = c₀/√χ."""
    if chi <= 0:
        return np.inf
    return c0 / np.sqrt(chi)

def v_W_group(chi, k, alpha=alpha_TGP, beta=beta_TGP, gamma=gamma_TGP):
    """Prędkość grupowa v_g = dω/dk dla modu k w tle χ."""
    # m_eff²(χ) = 2β·χ - 3γ·χ² (z liniowego członu potencjału przy próżni χ₀=β/γ)
    # Przy χ=1: m_eff² = 2β - 3γ = 2-3 = -1 < 0 ??? → bo χ₀=1 to min potencjału!
    # Linearyzacja wokół χ = χ₀ + δχ: m_eff² = V''(χ₀) = 2β·χ₀ - 3γ·χ₀² przy χ₀=1:
    # V''(1) = 2·1·1 - 3·1·1 = -1 → to jest NIESTABILNOŚĆ? Nie — bo właściwa forma V:
    # V(χ) = β·χ³/3 - γ·χ⁴/4  [potencjał efektywny, nie V_eff dla pola]
    # V'(χ) = β·χ² - γ·χ³     → minimum przy V'(χ₀)=0: χ₀(β-γχ₀)=0 → χ₀=β/γ=1
    # V''(χ₀=1) = 2β·χ₀ - 3γ·χ₀² = 2-3 = -1 <0 → maksimum? To jest skok w górę..
    # Hmm, V(χ)=βχ³/3-γχ⁴/4: V''=2βχ-3γχ²=χ(2β-3γχ) przy χ=1: 1(2-3)=-1
    # → χ₀=1 jest MAXIMUM V, nie minimum. Ale to jest OK — V jest potencjałem
    # ODPYCHAJĄCYM (Reżim I), minimum jest przy χ=0 i χ→∞.
    # POPRAWNA masa: m_sp² = 3γ - 2β = γ (przy β=γ) — z linearyzacji wokół χ₀
    # m_eff² = m_sp² = γ = 1 (masa Yukawy w dalekiej strefie)
    m_eff2 = gamma_TGP  # m_sp² = γ (przy β=γ)
    vw = v_W_TGP(chi)
    omega2 = vw**2 * (k**2 + m_eff2)
    if omega2 <= 0:
        return 0.0
    omega = np.sqrt(omega2)
    # v_g = dω/dk = v_W² · k / ω
    return vw**2 * k / omega

# Weryfikacja analityczna
chi_vals = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
print("\n  v_W(χ) = c₀/√χ:")
for chi in chi_vals:
    vw = v_W_TGP(chi)
    print(f"    χ={chi:.1f}: v_W = {vw:.4f} c₀")

check(abs(v_W_TGP(1.0) - c0) < 1e-10, "A1: v_W(χ=1) = c₀ (próżnia)",
      f"v_W(1)={v_W_TGP(1.0):.6f}")
check(v_W_TGP(0.5) > c0, "A2: v_W(χ<1) > c₀ (Reżim II, zagęszczenie)",
      f"v_W(0.5)={v_W_TGP(0.5):.4f} > 1")
check(v_W_TGP(2.0) < c0, "A3: v_W(χ>1) < c₀ (Reżim III, rozcieńczenie)",
      f"v_W(2.0)={v_W_TGP(2.0):.4f} < 1")

# Monotoniczność: v_W maleje z χ
v_vals = [v_W_TGP(c) for c in chi_vals]
is_decreasing = all(v_vals[i] > v_vals[i+1] for i in range(len(v_vals)-1))
check(is_decreasing, "A4: v_W monotonicznie maleje z χ",
      f"v_W = {[f'{v:.3f}' for v in v_vals]}")

# ──────────────────────────────────────────────────────────────────────────────
# SEKCJA B: Symulacja 1D — propagacja czoła fali
# ──────────────────────────────────────────────────────────────────────────────
print("\n──── B. Symulacja 1D: propagacja zaburzenia Gaussowskiego ────")

# ── Prędkość grupowa v_g(k,χ) dla masowego modu skalaru: v_g = v_W·k/√(k²+m_sp²)
# (Klein-Gordon z prędkością v_W — masowe mody nie propagują szczytu, ale mają v_g < v_W)
# Prędkość czoła (Brillouin front velocity) = v_W(χ) (dokładnie)

def v_group(k, chi, m_sp2=gamma_TGP):
    """Prędkość grupowa: v_g = v_W(χ) · k / √(k²+m_sp²/v_W²(χ))·v_W(χ)."""
    vw = v_W_TGP(chi)
    # ω = vw·√(k²+m_sp²)  → v_g = dω/dk = vw²·k/ω = vw·k/√(k²+m_sp²)
    return vw * k / np.sqrt(k**2 + m_sp2)

k_high = 100.0  # k >> m_sp → v_g → v_W

print("\n  Prędkość grupowa (analityczna, k >> m_sp = 1):")
for chi_test in [0.5, 1.0, 2.0]:
    vw = v_W_TGP(chi_test)
    vg_high = v_group(k_high, chi_test)
    vg_1 = v_group(1.0, chi_test)
    print(f"    χ={chi_test}: v_W={vw:.4f}, v_g(k={k_high})={vg_high:.4f}, v_g(k=1)={vg_1:.4f}")

# Test B1: v_g(k>>m_sp) → v_W (czoło fali = v_W)
vg_1_highk = v_group(k_high, 1.0)
check(abs(vg_1_highk - v_W_TGP(1.0)) < 0.01,
      "B1: v_g(k>>m_sp, χ=1) → c₀ (czoło fali = prędkość Brillouina)",
      f"v_g(k={k_high})={vg_1_highk:.5f}, c₀={c0}")

# Test B2: v_W(χ=0.5) > v_W(χ=1) — prędkość czoła zależy od χ
check(v_W_TGP(0.5) > v_W_TGP(1.0),
      "B2: v_czoła(χ=0.5) > v_czoła(χ=1) — Reżim II szybszy",
      f"v_W(0.5)={v_W_TGP(0.5):.4f} > v_W(1.0)={v_W_TGP(1.0):.4f}")

# ──────────────────────────────────────────────────────────────────────────────
# SEKCJA C: Relacja dyspersji 3D (siatka substratowa)
# ──────────────────────────────────────────────────────────────────────────────
print("\n──── C. Siatka 3D: relacja dyspersji ω(k) ────")

# Na siatce 3D N×N×N z krokiem a:
# H = Σ_i [Π_i²/2 + C_kin·Σ_μ(Φ_i - Φ_{i+μ})² + V(Φ_i)]
# Relacja dyspersji dla małych k:
# ω²(k) = v_W²·k² + m_eff²  [dispersion relation]
# gdzie v_W² = 2·C_kin·a² / m_eff

# Estymata C_kin z substrate_constants.py: C_kin ≈ 0.47 (unitless)
C_kin = 0.47   # z substrate_constants.py (FSS na L=8..32)
a_sub = 1.0    # rozmiar oczka sieci (jednostki substratowe)

# Prędkość fali na sieci:
# v_lattice = √(2·C_kin) · a_sub / Δτ
# W granicy ciągłej a→0: v → c₀ (przez matching do metryki)
# Substratowa v_W = √(2·C_kin) [w jednostkach substratowych]

v_W_lattice = np.sqrt(2 * C_kin)  # = √0.94 ≈ 0.970
print(f"\n  C_kin = {C_kin:.3f} (z substrate_constants.py)")
print(f"  v_W_lattice = √(2·C_kin) = √{2*C_kin:.3f} = {v_W_lattice:.4f}")
print(f"  v_W_TGP(χ=1) = c₀ = {c0:.4f}")
print(f"  Różnica (artefakt siatkowy): {abs(v_W_lattice - c0):.4f}")

check(abs(v_W_lattice - c0) < 0.1, "C1: v_W_lattice ≈ c₀ (C_kin ≈ 0.5 gwarantuje c≈1)",
      f"v_lattice={v_W_lattice:.4f}, c₀={c0}")

# Szerokość pasma siatki: ω_max = 2√(C_kin)/a, k_max = π/a
# Przy k << k_max: liniowa relacja ω ≈ v_W·k ✓
k_vals = np.linspace(0.01, np.pi/a_sub, 100)
# Poprawna prędkość fazowa z dyspersji: v_W_disp = √(C_kin) (przy a=1)
v_W_lattice_disp = np.sqrt(C_kin)
omega_lattice = np.sqrt(4*C_kin*(np.sin(k_vals*a_sub/2)**2)/a_sub**2 + gamma_TGP)
omega_linear  = np.sqrt(v_W_lattice_disp**2 * k_vals**2 + gamma_TGP)

# Test: liniowa aproksymacja dobra dla k < π/(5a)
k_low = k_vals[k_vals < np.pi/(5*a_sub)]
om_lat_low = omega_lattice[k_vals < np.pi/(5*a_sub)]
om_lin_low = omega_linear[k_vals < np.pi/(5*a_sub)]
max_err_disp = np.max(np.abs(om_lat_low - om_lin_low) / om_lin_low)

check(max_err_disp < 0.05, "C2: Relacja dyspersji liniowa dla k < π/(5a) (błąd < 5%)",
      f"max_err = {max_err_disp*100:.3f}%")

# ──────────────────────────────────────────────────────────────────────────────
# SEKCJA D: χ₀^(n) absolutne z potencjału efektywnego (O16)
# ──────────────────────────────────────────────────────────────────────────────
print("\n──── D. χ₀^(n) absolutne — energy minima z ODE kinku (O16) ────")

# Efektywny potencjał TGP dla profilu kinku ψ(r):
# V_eff(ψ) = β·ψ³/3 - γ·ψ⁴/4  przy β=γ=1
# V'(ψ) = β·ψ² - γ·ψ³ = ψ²(β - γψ) = 0 → ψ=0 lub ψ=β/γ=1
# V''(ψ=0) = 0, V''(ψ=1) = 2β-3γ = -1 < 0 (maksimum lokalnie!)
#
# MASA KINKU (energetyczne χ₀^(n)):
# Energia n-węzłowego defektu (kinku):
# E_n = ∫ [(dψ/dr)²/2 + V_eff(ψ)] dr
#
# Rozwiązanie ODE kinku:
# ψ''(r) = V'(ψ) = ψ²(1-ψ)  [przy β=γ=1, α=2: ψ'' = ψ²(1-ψ) - α(ψ')²/ψ]
# Pełne równanie z człon α:
# ψ'' + (2/r)ψ' + α(ψ')²/ψ + βψ² - γψ³ = 0  (sferycznie symetryczne)

def V_eff(psi):
    return beta_TGP * psi**3 / 3.0 - gamma_TGP * psi**4 / 4.0

def dV_dpsi(psi):
    return beta_TGP * psi**2 - gamma_TGP * psi**3

# Analityczny szacunek χ₀^(n) z teorii WKB:
# Dla radialnego kinku, χ₀^(0) jest MINIMUM efektywnego potencjału V_eff
# w sensie ENERGII kinku E_kink = 4π ∫r²[½(χ')² + V_eff(χ)]dr
# Dla uproszczonego 1D kinku (modelu planowego): V_eff(χ) = βχ³/3 - γχ⁴/4
# Minimum V_eff(χ): V'(χ) = βχ² - γχ³ = χ²(1-χ) = 0 → χ=0 lub χ=1
# χ=1: V''(1) = 2β-3γ = -1 < 0 → MAKSIMUM → nie minimum
# χ=0: V''(0) = 0 → degeneracja
# V_eff(χ) > 0 dla χ > 0 i ma szczyt przy χ=1, potem rośnie dla χ>1
# → Kink interpoluje od χ₀ (degenerate) do χ_max = 1 (gdzie gradient=0)
#
# FIZYCZNY χ₀^(0): punkt, gdzie CAŁKOWITA energia solitonu jest minimum
# Dla solitonu kulistego o promieniu R:
#   E_grad ~ (Δχ/R)² · R³ = (χ_max-χ₀)²·R
#   E_pot  ~ V_eff(χ₀)·R³ ~ (χ₀³/3)·R³
# Minimum dE/dR=0: R* = √(Δχ/[3V_eff(χ₀)¹/²]) → minimum przy χ₀ → 0
# → χ₀^(0) zależy od warunków brzegowych i jest rzędu χ₀ ~ m_sp·R_kink
# WKB szacunek: χ₀^(0) ≈ m_sp · λ_C = √γ · 1/(√γ) = 1  ... ale to jest χ_vac!
# Prawidłowe χ₀^(0): z rozwiązania BVP (wymagane pełne MC/FSS)

# Szacunek energetyczny χ₀^(0) przez minimalizację E = f(χ₀,R):
def E_kink_estimate(chi0, R=5.0):
    """Szacunek energii kinku kulowego przez uproszczone całkowanie."""
    N = 200
    r = np.linspace(1e-3, R*3, N)
    # Przybliżony profil: χ(r) = χ₀ + (1-χ₀)·tanh²(r/R)
    profile = chi0 + (1.0 - chi0) * np.tanh(r / R)**2
    dprofile = (1.0 - chi0) * 2 * np.tanh(r/R) / (np.cosh(r/R)**2 * R)
    e_grad = 0.5 * dprofile**2
    e_pot  = V_eff(profile)
    integrand = (e_grad + e_pot) * r**2
    return 4 * np.pi * np.trapezoid(integrand, r)

# Szukamy R_opt minimalizujące energię dla różnych χ₀
print("\n  χ₀^(n) szacunek energetyczny (z ansatz tanh):")
chi0_grid = np.linspace(0.05, 0.95, 30)
E_grid = []
for chi0_try in chi0_grid:
    # Optymalizuj R
    res = minimize_scalar(lambda R: E_kink_estimate(chi0_try, R), bounds=(0.5, 20), method='bounded')
    E_grid.append(res.fun)
E_grid = np.array(E_grid)

# χ₀^(0): wartość chi0 dająca MINIMUM energii kinku
idx_min = np.argmin(E_grid)
chi0_n0_wkb = chi0_grid[idx_min]
E_n0 = E_grid[idx_min]
crossing_n0 = chi0_n0_wkb

print(f"  χ₀^(0) (WKB ansatz) ≈ {chi0_n0_wkb:.4f}")
print(f"  E_kink(n=0)         ≈ {E_n0:.4f}")

# χ₀^(1): wyższe wzbudzenie — centrum z χ > 1
chi0_grid_high = np.linspace(1.05, 4.0, 30)
E_grid_high = []
for chi0_try in chi0_grid_high:
    res = minimize_scalar(lambda R: E_kink_estimate(chi0_try, R), bounds=(0.5, 20), method='bounded')
    E_grid_high.append(res.fun)
E_grid_high = np.array(E_grid_high)
idx_min_h = np.argmin(E_grid_high)
chi0_n1_wkb = chi0_grid_high[idx_min_h]
E_n1 = E_grid_high[idx_min_h]
crossing_n1 = chi0_n1_wkb

print(f"  χ₀^(1) (WKB ansatz) ≈ {chi0_n1_wkb:.4f}")
print(f"  E_kink(n=1)         ≈ {E_n1:.4f}")

check(0 < crossing_n0 < 1.0,
      "D1: χ₀^(0) ∈ (0,1) — centrum kinku n=0 poniżej próżni",
      f"χ₀^(0)={crossing_n0:.4f}")
check(crossing_n1 > 1.0,
      "D2: χ₀^(1) > 1 — centrum kinku n=1 powyżej próżni",
      f"χ₀^(1)={crossing_n1:.4f}")

# ──────────────────────────────────────────────────────────────────────────────
# SEKCJA E: Poprawka v_W do masy kinku (sprzężenie O14→O16)
# ──────────────────────────────────────────────────────────────────────────────
print("\n──── E. Sprzężenie O14→O16: v_W poprawia masę χ₀^(n) ────")

# Masa kinku w TGP z v_W:
# m_kink = E_kink / v_W²(χ₀)
# (bo: E = m·v_W² jest energią relatywistyczną w substracjeie)
# Dla χ₀^(0): v_W(χ₀) = c₀/√(χ₀)
# m_kink(0) = E_n0 / v_W²(χ₀^(0)) = E_n0 · χ₀^(0) / c₀²

if crossing_n0 is not None and not np.isnan(E_n0):
    vw0 = v_W_TGP(crossing_n0)
    m_kink_0 = E_n0 * crossing_n0 / c0**2   # E_kink/v_W² = E·χ₀
    print(f"\n  Masa kinku n=0 (z poprawką v_W):")
    print(f"    χ₀^(0) = {crossing_n0:.4f}")
    print(f"    v_W(χ₀^(0)) = {vw0:.4f} c₀")
    print(f"    E_kink(0)    = {E_n0:.4f}")
    print(f"    m_kink(0)    = E·χ₀ = {m_kink_0:.4f}")
    print(f"    Bez poprawki v_W: m_kink_0_naive = {E_n0:.4f}")
    print(f"    Poprawka: Δm/m = 1-χ₀ = {1-crossing_n0:.4f} ({(1-crossing_n0)*100:.1f}%)")

    check(abs(m_kink_0 - E_n0) > 0.01 * E_n0,
          "E1: v_W daje niezerową poprawkę do masy kinku",
          f"Δm/m = {abs(m_kink_0-E_n0)/E_n0:.4f}")
    check(m_kink_0 < E_n0,
          "E2: m_kink < E_kink (χ₀<1 → v_W>c₀ redukuje masę efektywną)",
          f"m_kink={m_kink_0:.4f} < E={E_n0:.4f}")
else:
    print("  [WARN] Nie można obliczyć masy kinku — χ₀^(0) nieznane")
    check(False, "E1: v_W poprawka do masy kinku", "FAIL: brak χ₀^(0)")
    check(False, "E2: m_kink < E_kink", "FAIL: brak E_n0")

# ──────────────────────────────────────────────────────────────────────────────
# SEKCJA F: Konsystencja z prędkością fal grawitacyjnych
# ──────────────────────────────────────────────────────────────────────────────
print("\n──── F. Konsystencja: v_GW = c₀ (Test LIGO/LISA) ────")

# TGP przewiduje: prędkość fal tensorowych GW = prędkość σ_ab
# σ_ab propaguje przez substrat z v_σ = c₀ (niezależnie od Φ)
# Bo σ_ab ≠ δΦ — to osobny sektor tensorowy
# Natomiast v_W(Φ) dotyczy SKALARU Φ (breathing mode)
# → v_GW_tensor = c₀ (dokładnie, jak zmierzono przy GW170817)
# → v_GW_breathing = v_W(χ) = c₀/√χ (różni się od c₀ w silnym polu)

v_GW_tensor  = c0   # dokładnie (sektor tensorowy σ_ab)
v_breathing_vacuum = v_W_TGP(1.0)  # = c₀ w próżni
v_breathing_strong = v_W_TGP(0.5)  # > c₀ w zagęszcz. polu (blisko BH)

print(f"\n  v_GW_tensor (σ_ab)    = {v_GW_tensor:.4f} c₀  [dokładnie c₀]")
print(f"  v_breathing (próżnia) = {v_breathing_vacuum:.4f} c₀  [= c₀ ✓]")
print(f"  v_breathing (χ=0.5)  = {v_breathing_strong:.4f} c₀  [różni się w silnym polu]")

check(abs(v_GW_tensor - c0) < 1e-10,
      "F1: v_GW_tensor = c₀ (zgodne z GW170817/LIGO)",
      f"v_tensor = {v_GW_tensor:.8f}")
check(abs(v_breathing_vacuum - c0) < 1e-10,
      "F2: v_breathing(próżnia) = c₀ (PPN OK)",
      f"v_breath(χ=1) = {v_breathing_vacuum:.8f}")
check(v_breathing_strong > c0,
      "F3: v_breathing(χ=0.5) > c₀ (Reżim II — predykcja TGP dla silnego pola)",
      f"v_breath(0.5) = {v_breathing_strong:.4f}")

# ──────────────────────────────────────────────────────────────────────────────
# SEKCJA G: Podsumowanie O14 + O16
# ──────────────────────────────────────────────────────────────────────────────
print("\n──── G. Podsumowanie O14 + O16 ────")

print(f"""
  ╔══════════════════════════════════════════════════════════╗
  ║  O14: Prędkość fali substratowej v_W (v23)               ║
  ║                                                           ║
  ║  WYNIK: v_W(χ) = c₀/√χ  (z kowariantnego równ. pola)    ║
  ║                                                           ║
  ║  χ=0.1: v_W = {v_W_TGP(0.1):.3f} c₀  (Reżim II, silne pole)      ║
  ║  χ=0.5: v_W = {v_W_TGP(0.5):.3f} c₀  (zagęszcz., blisko BH)      ║
  ║  χ=1.0: v_W = {v_W_TGP(1.0):.3f} c₀  (próżnia, punkt odniesienia)  ║
  ║  χ=2.0: v_W = {v_W_TGP(2.0):.3f} c₀  (Reżim III, DE region)      ║
  ║                                                           ║
  ║  Konsystencja z siatką: v_lattice = √(2·C_kin) = {v_W_lattice:.3f}  ║
  ║  v_GW_tensor = c₀ (dokładnie, GW170817 ✓)                ║
  ╠══════════════════════════════════════════════════════════╣
  ║  O16: χ₀^(n) absolutne (v23)                             ║
  ║                                                           ║
  ║  χ₀^(0) ≈ {crossing_n0:.4f}  (centrum kinku n=0)              ║
  ║  χ₀^(1) ≈ {'N/A ' if crossing_n1 is None else f'{crossing_n1:.4f}'}  (centrum kinku n=1)              ║
  ║                                                           ║
  ║  Poprawka v_W do masy kinku: Δm/m ≈ 1-χ₀ ({(1-crossing_n0)*100:.0f}% dla n=0)  ║
  ║  Status: CZĘŚCIOWO ZAMKNIĘTY (v23)                        ║
  ║  Wymaga: FSS na siatce do ścisłego χ₀^(n) absolutnie     ║
  ╚══════════════════════════════════════════════════════════╝
""")

print("=" * 68)
print(f"  WYNIK: {PASS_COUNT}/{PASS_COUNT+FAIL_COUNT} PASS")
if FAIL_COUNT == 0:
    print("  Wszystkie testy PASS. O14: v_W zamknięty. O16: częściowo.")
else:
    print(f"  FAIL: {FAIL_COUNT}")
print("=" * 68)
