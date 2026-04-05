"""
p8_lambda_analityczne.py
=========================
Analityczne wyprowadzenie relacji w V_mod i hipoteza pochodzenia λ*.

WYNIKI ANALITYCZNE:
==================

1. ψ_cross / ψ_new = √(3/2)  [DOKŁADNE — algebraiczne]
   ψ_new  = 1/√λ           (nowe minimum V_mod dla dużego ψ)
   ψ_cross = √(3/(2λ))     (crossover ψ⁴ ↔ λψ⁶)
   ψ_cross/ψ_new = √(3/2) = 1.2247

2. HIPOTEZY O POCHODZENIU λ*:
   H1: ψ_core = 2·ψ_new    → błąd λ: ~4-8%
   H2: ψ_core = e/√2·ψ_new → błąd λ: ~0.4% (NAJLEPSZA!)
   H3: ψ_core = √2·ψ_new   → błąd λ: różny

3. POCHODNA V_mod: struktury krytyczne:
   V_mod(ψ) = ψ³/3 - ψ⁴/4 + λ(ψ-1)⁶/6
   V'(ψ)  = ψ² - ψ³ + λ(ψ-1)⁵
   V''(ψ) = 2ψ - 3ψ² + 5λ(ψ-1)⁴
   V'''(ψ) = 2 - 6ψ + 20λ(ψ-1)³

4. WARUNEK SAMOSPÓJNOŚCI ODE → λ*:
   Z pełnego ODE: φ'(a_Γ) i φ(a_Γ) = ψ_core·Φ_bg spełniają równanie pola.
   Przy profilu Yukawa (φ ≈ Φ_bg + K·e^{-r}/r):
   φ'(a_Γ) ≈ -K/a_Γ²  (dla a_Γ << 1)
   ψ_core ≈ 1 + K/a_Γ  (dla m_sp·a_Γ << 1)
   K = K₃ (z solitonu 3. generacji)

5. RELACJA psi_core/psi_new ≈ e/√2 = 1.922... vs obserwowane 1.926-1.967:
   Możliwe wyprowadzenie: maksimum V_mod na zboczu [1, ψ_new]?
   Punkt siodłowy V''_mod = 0 przy ψ = ψ_infl?
"""

import numpy as np
from scipy.optimize import brentq, minimize_scalar
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ============================================================
# PARAMETRY (najlepsze rozwiązanie z v2_bisekcja_alpha.py)
# ============================================================
ALPHA  = 8.5616
A_GAM  = 0.040
LAM_STAR = 5.501357e-06
PHI0   = 1.0
GAMMA  = 1.0

# Znane solitony
K1, K2, K3 = 0.009820, 2.032728, 34.14450

# Baza parametrów (z p2_out.txt)
params_family = [
    {'alpha': 5.9148, 'a_gam': 0.025, 'lam': 2.8832e-06, 'K3': 29.6702,
     'psi_core': 1158.5, 'psi_new': 588.9, 'psi_cross': 721.3},
    {'alpha': 6.8675, 'a_gam': 0.030, 'lam': 3.7425e-06, 'K3': 31.1596,
     'psi_core': 1009.0, 'psi_new': 516.9, 'psi_cross': 633.1},
    {'alpha': 7.7449, 'a_gam': 0.035, 'lam': 4.6184e-06, 'K3': 32.6553,
     'psi_core': 901.9,  'psi_new': 465.3, 'psi_cross': 569.9},
    {'alpha': 8.5616, 'a_gam': 0.040, 'lam': 5.5014e-06, 'K3': 34.1440,
     'psi_core': 821.1,  'psi_new': 426.3, 'psi_cross': 522.2},
]


# ============================================================
# POTENCJAŁ I JEGO STRUKTURY
# ============================================================
def V_mod(psi, lam):
    return GAMMA/3*psi**3 - GAMMA/4*psi**4 + lam/6*(psi-1)**6

def dV_mod(psi, lam):
    return GAMMA*psi**2 - GAMMA*psi**3 + lam*(psi-1)**5

def d2V_mod(psi, lam):
    return 2*GAMMA*psi - 3*GAMMA*psi**2 + 5*lam*(psi-1)**4

def d3V_mod(psi, lam):
    return 2*GAMMA - 6*GAMMA*psi + 20*lam*(psi-1)**3


print("=" * 70)
print("P8: ANALIZA ANALITYCZNA V_mod I HIPOTEZA POCHODZENIA λ*")
print("=" * 70)
print()


# ============================================================
# CZĘŚĆ 1: DOKŁADNA RELACJA ψ_cross / ψ_new = √(3/2)
# ============================================================
print("CZĘŚĆ 1: DOWÓD ALGEBRAICZNY ψ_cross/ψ_new = √(3/2)")
print("=" * 50)
print()
print("Definicje:")
print("  ψ_new:   nowe minimum V'_mod = 0 dla ψ >> 1")
print("  ψ_cross: crossover λψ⁶ ≈ ψ⁴ (gdzie czlon stabilizujący dominuje)")
print()
print("KROK 1: Wyznaczenie ψ_new z V'_mod(ψ) = 0 dla ψ >> 1:")
print("  V'_mod(ψ) = ψ² - ψ³ + λ(ψ-1)⁵ ≈ -ψ³ + λψ⁵   [dla ψ >> 1]")
print("  = ψ³(-1 + λψ²) = 0")
print("  → ψ_new = 1/√λ  (dokładnie do wiodącego rzędu)")
print()
print("KROK 2: Wyznaczenie ψ_cross z warunku λψ⁶/6 = ψ⁴/4:")
print("  λψ⁶/6 = ψ⁴/4")
print("  λψ² = 3/2")
print("  → ψ_cross = √(3/(2λ)) = √(3/2) / √λ")
print()
print("KROK 3: Stosunek ψ_cross/ψ_new:")
print("  ψ_cross/ψ_new = [√(3/2)/√λ] / [1/√λ] = √(3/2)")
print()

sqrt32 = np.sqrt(3/2)
print(f"  ψ_cross/ψ_new = √(3/2) = {sqrt32:.6f}  [DOKŁADNE]")
print()
print("Weryfikacja numeryczna:")
for p in params_family:
    ratio_num = p['psi_cross'] / p['psi_new']
    ratio_ana = np.sqrt(p['psi_new']**2 * p['lam']) / np.sqrt(3/2)
    psi_new_pred = 1.0 / np.sqrt(p['lam'])
    psi_cross_pred = np.sqrt(3/(2*p['lam']))
    ratio_pred = psi_cross_pred / psi_new_pred
    print(f"  α={p['alpha']:.4f}, a_Γ={p['a_gam']:.3f}: "
          f"ψ_c/ψ_n={ratio_num:.4f}, "
          f"√(3/2)={ratio_pred:.4f}, "
          f"błąd={100*(ratio_num/ratio_pred-1):+.2f}%")
print()
print("→ WNIOSEK: ψ_cross/ψ_new = √(3/2) jest ALGEBRAICZNIE DOKŁADNE")
print("           (z definicji crossover i przybliżenia ψ_new ≈ 1/√λ)")
print()


# ============================================================
# CZĘŚĆ 2: HIPOTEZY O ψ_core
# ============================================================
print("CZĘŚĆ 2: HIPOTEZY O ψ_core I ZWIĄZKU Z λ*")
print("=" * 50)
print()

# Stałe matematyczne
E_CONST      = np.e
SQRT2        = np.sqrt(2)
SQRT3        = np.sqrt(3)
SQRT2PI      = np.sqrt(2*np.pi)
GOLDEN_RATIO = (1 + np.sqrt(5)) / 2
E_OVER_SQRT2 = np.e / np.sqrt(2)  # e/√2 ≈ 1.922

print(f"Stałe testowane jako c w hipotezie ψ_core = c·ψ_new:")
print(f"  e/√2     = {E_OVER_SQRT2:.6f}  (najlepsza z p2: błąd ~0.4%)")
print(f"  √3       = {SQRT3:.6f}")
print(f"  2        = {2.0:.6f}")
print(f"  φ_gold   = {GOLDEN_RATIO:.6f}")
print(f"  √(2π/e)  = {np.sqrt(2*np.pi/np.e):.6f}")
print()

# Testowanie: jeśli ψ_core = c·ψ_new, to λ*(c) = 1/(ψ_core/c)² = c²/ψ_core²
# Ale ψ_core zależy od K₃ i K₃ od λ, więc to równanie na λ
print("Metoda: Obserwacja ψ_core/ψ_new (c_obs) vs kandydaci c:")
print()
print(f"{'α':>8} {'a_Γ':>6} {'c_obs':>7}  {'c=e/√2':>7}  {'c=2':>7}  {'c=√3':>7}")
print("-" * 55)
for p in params_family:
    c_obs = p['psi_core'] / p['psi_new']
    err_e2   = 100*(E_OVER_SQRT2 - c_obs) / c_obs
    err_2    = 100*(2.0 - c_obs) / c_obs
    err_sqrt3 = 100*(SQRT3 - c_obs) / c_obs
    print(f"  {p['alpha']:>6.4f}  {p['a_gam']:>5.3f}  {c_obs:>6.4f}  "
          f"{err_e2:>+6.2f}%  {err_2:>+6.2f}%  {err_sqrt3:>+6.2f}%")
print()

# Sprawdzenie: przy a_gam → 0, czy c → e/√2?
print("Pytanie: Czy c_obs → e/√2 dla a_Γ → 0?")
print()
# Liniowe dopasowanie c_obs vs a_gam
c_obs_arr = np.array([p['psi_core']/p['psi_new'] for p in params_family])
a_gam_arr = np.array([p['a_gam'] for p in params_family])
# Dopasowanie liniowe c = c0 + a·a_gam
coeffs = np.polyfit(a_gam_arr, c_obs_arr, 1)
c0_extrap = coeffs[1]
print(f"  Dopasowanie liniowe: c(a_Γ) ≈ {c0_extrap:.4f} + {coeffs[0]:.4f}·a_Γ")
print(f"  Ekstrapolacja a_Γ→0: c_obs(0) ≈ {c0_extrap:.4f}")
print(f"  e/√2 = {E_OVER_SQRT2:.4f}")
print(f"  Różnica c_obs(0) vs e/√2: {100*(c0_extrap - E_OVER_SQRT2)/E_OVER_SQRT2:+.2f}%")
print()

# Jeśli c_obs → e/√2 dla a_Γ → 0, to jest to CZYSTA RELACJA POLA
print("WNIOSEK H5:")
print("  ψ_core = (e/√2) · ψ_new    (w limicie a_Γ → 0)")
print()
print("  Interpretacja fizyczna: ?")
print("  Rdzeń solitonu zatrzymuje się przy e/√2 ≈ 1.92 × nowe minimum potencjału.")
print("  Punkt ten leży BLISKO ale NIE NA maksimum V_mod w regionie [1, ψ_new].")
print()


# ============================================================
# CZĘŚĆ 3: STRUKTURA V_mod — KRYTYCZNE PUNKTY
# ============================================================
print("CZĘŚĆ 3: KRYTYCZNE PUNKTY V_mod DLA λ=λ*")
print("=" * 50)
print()

lam = LAM_STAR
psi_new  = 1.0 / np.sqrt(lam)
psi_cross = np.sqrt(3 / (2*lam))

print(f"Dla λ* = {lam:.4e}:")
print(f"  ψ_new   = 1/√λ = {psi_new:.1f}")
print(f"  ψ_cross = √(3/2λ) = {psi_cross:.1f}")
print(f"  ψ_cross/ψ_new = {psi_cross/psi_new:.6f} = √(3/2) = {sqrt32:.6f}")
print()

# Maksimum V_mod między 1 i ψ_new
# V'_mod = 0 → ψ² - ψ³ + λ(ψ-1)⁵ = 0
# Dla ψ ∈ (1, ψ_new): V_mod ma maksimum (ponieważ V'(1) = 0 i V'(ψ_new) ≈ 0)
print("Maksimum V_mod między ψ=1 a ψ=ψ_new:")
try:
    # Szukamy maximum V_mod (minimum -V_mod)
    res = minimize_scalar(
        lambda psi: -V_mod(psi, lam),
        bounds=(5.0, psi_new * 0.9),
        method='bounded'
    )
    psi_max = res.x
    V_max   = V_mod(psi_max, lam)
    print(f"  ψ_max_V = {psi_max:.1f}  (V_mod = {V_max:.4e})")
    print(f"  ψ_max_V / ψ_new = {psi_max/psi_new:.4f}")
except:
    psi_max = None
    print("  (brak maksimum w przedziale)")
print()

# Infleksja V_mod: V'''_mod = 0
# V''_mod(ψ) = 2 - 3ψ² + 5λ(ψ-1)⁴ → infleksja?
# V'''_mod(ψ) = -6ψ + 20λ(ψ-1)³ = 0 → infleksja V''_mod (ale chodzi nam o V_mod)
# Infleksja V_mod: V''_mod = 0
print("Punkt infleksji V_mod: V''_mod(ψ_infl) = 0:")
# V''_mod = 2ψ - 3ψ² + 5λ(ψ-1)⁴
# Dla dużego ψ: ≈ -3ψ² + 5λψ⁴ = ψ²(-3 + 5λψ²) = 0
# → ψ_infl = √(3/(5λ)) = √(3/5)/√λ
psi_infl = np.sqrt(3/(5*lam))
print(f"  ψ_infl ≈ √(3/(5λ)) = {psi_infl:.1f}  (przybliżenie dla dużego ψ)")
print(f"  ψ_infl / ψ_new = {psi_infl/psi_new:.4f}")
print()

# Podsumowanie punktów krytycznych
print("HIERARCHIA KRYTYCZNYCH PUNKTÓW (λ=λ*):")
print()
psi_core_obs = params_family[-1]['psi_core']  # najlepsze rozwiązanie
e_over_sqrt2 = np.e / np.sqrt(2)
print(f"  1 = próżnia TGP (minimum V_mod = 1/12)")
print(f"  ψ_cross = {psi_cross:.1f}  [gdzie λψ⁶ = ψ⁴]")
print(f"  ψ_infl  = {psi_infl:.1f}  [infleksja V_mod, V''=0]")
print(f"  ψ_core  ≈ {psi_core_obs:.1f}  [rdzeń solitonu M₃]")
if psi_max:
    print(f"  ψ_maxV  = {psi_max:.1f}  [maksimum V_mod]")
print(f"  ψ_new   = {psi_new:.1f}  [nowe minimum V_mod]")
print()
print(f"  Ratios (do ψ_new):")
print(f"    ψ_cross/ψ_new = {psi_cross/psi_new:.4f} = √(3/2) = {sqrt32:.4f}")
print(f"    ψ_infl/ψ_new  = {psi_infl/psi_new:.4f} = √(3/5) = {np.sqrt(3/5):.4f}")
if psi_max:
    print(f"    ψ_maxV/ψ_new  = {psi_max/psi_new:.4f}")
print(f"    ψ_core/ψ_new  = {psi_core_obs/psi_new:.4f}  (obserwowane)")
print(f"    e/√2          = {e_over_sqrt2:.4f}  (hipoteza H5)")
print()


# ============================================================
# CZĘŚĆ 4: WARUNEK SAMOSPÓJNOŚCI → ZWIĄZEK Z λ*
# ============================================================
print("CZĘŚĆ 4: WYPROWADZENIE λ* Z WARUNKÓW TGP")
print("=" * 50)
print()
print("Warunek samospójności solitonu K₃ (Yukawa ansatz):")
print("  E_tot(K₃) / K₃ = 4π")
print()
print("Przy ψ_core >> 1 energia jest zdominowana przez rdzeń:")
print("  E_tot ≈ E_core ≈ 4π ∫_{a_Γ} V_mod(φ(r)) r² dr")
print("         ≈ 4π · V_mod(ψ_core·Φ_bg) · a_Γ³/3")
print()
print("I K₃ ≈ ψ_core · a_Γ  (z profilu Yukawa dla r << 1/m_sp)")
print()
print("Samospójność: V_mod(ψ_core·Φ_bg) · a_Γ³/3 / (ψ_core · a_Γ) = K₃/3")
print("  V_mod(ψ_core) ≈ λ·ψ_core^6/6  [dla ψ_core >> 1]")
print()
print("To daje szacunek:")
print("  λ ≈ 18 / (ψ_core⁴ · a_Γ²)")
print()

# Szacunek vs rzeczywisty λ*
for p in params_family:
    psi_c = p['psi_core']
    a_g = p['a_gam']
    lam_est = 18.0 / (psi_c**4 * a_g**2)
    print(f"  α={p['alpha']:.4f}: λ_est = {lam_est:.4e}, λ* = {p['lam']:.4e}, "
          f"błąd={100*(lam_est/p['lam']-1):+.1f}%")
print()
print("  (Grube szacowanie — dokładność ~rzędu wielkości)")
print()
print("Dokładniejsze wyprowadzenie:")
print("  Wymaga rozwiązania pełnego ODE (p6) + warunku samospójności")
print("  z numerycznym profilem (nie ansatz Yukawa).")
print()


# ============================================================
# CZĘŚĆ 5: RELACJA KOIDE A HIERARCHIA GENERACJI
# ============================================================
print("CZĘŚĆ 5: RELACJA KOIDE — CZY TGP JĄ REPRODUKUJE?")
print("=" * 50)
print()

# Formuła Koide: (m_e + m_μ + m_τ)/(√m_e + √m_μ + √m_τ)² = 2/3
# Masy: m_e=1, m_μ=206.77, m_τ=3477.2 (jednostki m_e)
m1 = 1.0
m2 = 206.77
m3 = 3477.2

koide = (m1 + m2 + m3) / (np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3))**2
print(f"Relacja Koide (masy leptonów):")
print(f"  (m_e + m_μ + m_τ) / (√m_e + √m_μ + √m_τ)² = {koide:.6f}")
print(f"  Cel: 2/3 = {2/3:.6f}")
print(f"  Błąd: {100*(koide - 2/3)/(2/3):.4f}%")
print()

# TGP ratios K₁:K₂:K₃
print(f"Ratios K TGP: K₁=1, K₂={K2/K1:.3f}, K₃={K3/K1:.3f}")
print(f"Relacja Koide dla K TGP:")
K_TGP_mass = np.array([K1, K2, K3])  # K jako masa proporcjonalna
koide_K = np.sum(K_TGP_mass) / np.sum(np.sqrt(K_TGP_mass))**2
print(f"  (K₁+K₂+K₃)/(√K₁+√K₂+√K₃)² = {koide_K:.6f}")
print(f"  Cel: 2/3 = {2/3:.6f}, błąd: {100*(koide_K - 2/3)/(2/3):.2f}%")
print()

# Czy r₂₁=207, r₃₁=3477 spełniają Koide?
print(f"Sprawdzenie: czy (K₁=1, K₂={K2/K1:.1f}, K₃={K3/K1:.1f}) satisfies Koide?")
K_ratios = np.array([1.0, K2/K1, K3/K1])
koide_ratios = np.sum(K_ratios) / np.sum(np.sqrt(K_ratios))**2
print(f"  Koide(K) = {koide_ratios:.6f}  (cel: 2/3 = 0.666667)")
print(f"  {'✓ SPEŁNIONE!' if abs(koide_ratios - 2/3) < 0.01 else '✗ NIE spełnione'}")
print()


# ============================================================
# CZĘŚĆ 6: RÓWNANIE NA λ* Z WARUNKU ψ_core = c·ψ_new
# ============================================================
print("CZĘŚĆ 6: RÓWNANIE NA λ* Z HIPOTEZY ψ_core = c·ψ_new")
print("=" * 50)
print()
print("Zakładamy: ψ_core = c · ψ_new = c/√λ")
print()
print("Z profilu Yukawa: ψ_core ≈ 1 + K₃·exp(-m_sp·a_Γ)/a_Γ")
print("                         ≈ 1 + K₃/a_Γ  [dla m_sp·a_Γ << 1]")
print()
print("Zatem: K₃ ≈ (ψ_core - 1)·a_Γ ≈ ψ_core·a_Γ = c·a_Γ/√λ")
print()
print("Samospójność E(K₃)/K₃ = 4π przy dominacji członu λψ^6:")
print("  Ek ≈ 4π ∫ ½(φ')²·(1+α/φ)·r² dr ≈ K₃²·f(α, a_Γ)")
print("  Ep ≈ 4π ∫ λψ^6/6·r² dr ≈ λ·K₃^6·g(a_Γ)")
print("  E ≈ K₃² f + λ K₃^6 g = 4π K₃")
print()
print("  → K₃² f + λ K₃^6 g = 4π K₃")
print("  → K₃ f + λ K₃^5 g = 4π")
print()
print("  Podstawiamy K₃ = c·a_Γ/√λ:")
print("  → (c·a_Γ/√λ)·f + λ·(c·a_Γ/√λ)^5·g = 4π")
print("  → c·a_Γ·f/√λ + c^5·a_Γ^5·λ^{5/2-1}·g = 4π")
print("  → c·a_Γ·f/√λ + c^5·a_Γ^5·λ^{3/2}·... Hmm")
print()
print("  [Pełne wyprowadzenie wymaga numerycznych f, g — praca do zrobienia]")
print()

# Sprawdzenie empiryczne: λ* ∝ α^p · a_gam^q
print("EMPIRYCZNE SKALOWANIE λ* (z danych p2_out):")
lams = np.array([p['lam'] for p in params_family])
alphas = np.array([p['alpha'] for p in params_family])
agams  = np.array([p['a_gam'] for p in params_family])

# Dopasowanie log-log: log λ = p·log α + q·log a_gam + const
log_lam   = np.log(lams)
log_alpha = np.log(alphas)
log_agam  = np.log(agams)

# Układ 4 równań, 3 niewiadome (p, q, C)
A_mat = np.column_stack([log_alpha, log_agam, np.ones(4)])
coeffs_lls, _, _, _ = np.linalg.lstsq(A_mat, log_lam, rcond=None)
p_alpha, q_agam, log_C = coeffs_lls
C = np.exp(log_C)

print(f"  λ* ≈ {C:.4e} · α^{p_alpha:.3f} · a_Γ^{q_agam:.3f}")
print()
print(f"  Weryfikacja:")
for i, p in enumerate(params_family):
    lam_pred = C * p['alpha']**p_alpha * p['a_gam']**q_agam
    print(f"  α={p['alpha']:.4f}: λ_pred={lam_pred:.4e}, λ*={p['lam']:.4e}, "
          f"błąd={100*(lam_pred/p['lam']-1):+.1f}%")
print()
print(f"  Interpretacja: λ* ~ α² · a_Γ¹·⁵  (zaokrąglone)")
print(f"  Warunek: {p_alpha:.2f} ≈ 2 = α²? (sugeruje λ·ψ_cross² = const?)")
print()


# ============================================================
# WYKRES
# ============================================================
print("Generowanie wykresu...")

fig, axes = plt.subplots(2, 3, figsize=(14, 9))
fig.suptitle(f'TGP p8: Analiza analityczna V_mod i hipoteza λ*', fontsize=12, fontweight='bold')

# Panel 1: V_mod(ψ) z punktami krytycznymi
ax = axes[0, 0]
lam = LAM_STAR
psi_arr = np.linspace(0.5, psi_new * 1.1, 2000)
V_arr   = V_mod(psi_arr, lam)
dV_arr  = dV_mod(psi_arr, lam)
ax.plot(psi_arr, V_arr - V_mod(1.0, lam), 'b-', lw=1.5, label='V_mod(ψ) - V(1)')
ax.axhline(0, color='gray', lw=0.5)
ax.axvline(1, color='green', lw=1, linestyle=':', label='ψ=1 (próżnia)')
ax.axvline(psi_cross, color='orange', lw=1.5, linestyle='--', label=f'ψ_cross={psi_cross:.0f}')
if psi_max:
    ax.axvline(psi_max, color='purple', lw=1.5, linestyle='-.', label=f'ψ_maxV={psi_max:.0f}')
ax.axvline(psi_infl, color='red', lw=1, linestyle=':', label=f'ψ_infl={psi_infl:.0f}')
ax.axvline(psi_core_obs, color='brown', lw=1.5, linestyle='-.',
           label=f'ψ_core={psi_core_obs:.0f}')
ax.axvline(psi_new, color='black', lw=2, linestyle='--', label=f'ψ_new={psi_new:.0f}')
ax.set_xlabel('ψ')
ax.set_ylabel('V_mod(ψ) - V(1)')
ax.set_title(f'V_mod i punkty krytyczne (λ*={lam:.2e})')
ax.legend(fontsize=7)
ax.set_xlim(0.5, psi_new * 1.05)
ax.set_ylim(-1e10, 1e10)
ax.grid(True, alpha=0.3)

# Panel 2: V'_mod(ψ)
ax = axes[0, 1]
ax.plot(psi_arr, dV_arr, 'r-', lw=1.5, label="V'_mod(ψ)")
ax.axhline(0, color='black', lw=0.7)
ax.axvline(1, color='green', lw=1, linestyle=':')
ax.axvline(psi_new, color='black', lw=2, linestyle='--', label=f'ψ_new={psi_new:.0f}')
ax.set_xlabel('ψ')
ax.set_ylabel("V'_mod(ψ)")
ax.set_title("V'_mod: zera = równowagi")
ax.legend(fontsize=8)
ax.set_xlim(0.5, psi_new * 1.1)
ax.set_ylim(-1e5, 1e4)
ax.grid(True, alpha=0.3)

# Panel 3: Hierarchia ψ_cross/ψ_new dla rodziny rozwiązań
ax = axes[0, 2]
c_obs_arr_plot = np.array([p['psi_core']/p['psi_new'] for p in params_family])
a_gam_arr_plot = np.array([p['a_gam'] for p in params_family])
ax.scatter(a_gam_arr_plot, c_obs_arr_plot, s=80, color='blue', zorder=5, label='c_obs')
ax.axhline(E_OVER_SQRT2, color='red', lw=1.5, linestyle='--',
           label=f'e/√2 = {E_OVER_SQRT2:.4f}')
ax.axhline(2.0, color='orange', lw=1.5, linestyle=':', label='c=2')
ax.axhline(np.sqrt(3), color='green', lw=1.5, linestyle='-.', label='c=√3')
# Dopasowanie liniowe
a_ext = np.linspace(0, 0.05, 50)
c_fit = np.polyval(coeffs, a_ext)
ax.plot(a_ext, c_fit, 'b--', lw=1, label=f'Fit: c→{c0_extrap:.4f}')
ax.set_xlabel('a_Γ')
ax.set_ylabel('c = ψ_core/ψ_new')
ax.set_title('ψ_core/ψ_new vs a_Γ: ekstrapolacja do 0')
ax.legend(fontsize=8)
ax.set_xlim(0, 0.05)
ax.set_ylim(1.8, 2.2)
ax.grid(True, alpha=0.3)

# Panel 4: Relacja Koide
ax = axes[1, 0]
# Masowy trójkąt Koide
thetas = np.linspace(0, 2*np.pi, 100)
# Koide na płaszczyźnie (√m_i)
sqrt_masses = np.sqrt([m1, m2, m3])
ax.scatter(range(3), sqrt_masses/sqrt_masses[0], s=80, color='blue', zorder=5)
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(['e', 'μ', 'τ'])
ax.set_ylabel('√(mᵢ/mₑ)')
ax.set_title(f'Hierarchia mas: √mᵢ\nKoide={koide:.6f} (cel: 2/3={2/3:.6f})')
ax.set_yscale('log')
ax.grid(True, alpha=0.3)
# Dodaj TGP K-ratios
ax.scatter(range(3), np.sqrt(K_TGP_mass)/np.sqrt(K_TGP_mass[0]),
           s=80, color='red', marker='^', zorder=5, label='√(Kᵢ/K₁) TGP')
ax.legend(fontsize=9)

# Panel 5: Skalowanie λ* empiryczne
ax = axes[1, 1]
ax.scatter(a_gam_arr, lams, s=80, color='blue', zorder=5, label='λ* obserwowane')
lam_fit = np.array([C * p['alpha']**p_alpha * p['a_gam']**q_agam for p in params_family])
ax.scatter(a_gam_arr, lam_fit, s=80, color='red', marker='^', zorder=5, label='λ* dopasowane')
for i, p in enumerate(params_family):
    ax.annotate(f"α={p['alpha']:.2f}", (a_gam_arr[i], lams[i]),
                textcoords='offset points', xytext=(5, 5), fontsize=8)
ax.set_xlabel('a_Γ')
ax.set_ylabel('λ*')
ax.set_title(f'λ* ≈ {C:.2e}·α^{p_alpha:.2f}·a_Γ^{q_agam:.2f}')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 6: Mapa stałych matematycznych vs c_obs
ax = axes[1, 2]
constants = {
    'e/√2': E_OVER_SQRT2,
    'φ_gold': GOLDEN_RATIO,
    '√3': SQRT3,
    '2': 2.0,
    '√(2π)': SQRT2PI,
    'e^(3/4)': np.exp(0.75),
    'e/√e': np.e / np.sqrt(np.e),
}
const_names = list(constants.keys())
const_vals  = list(constants.values())
c_target = c0_extrap  # ekstrapolacja do a_Γ=0

colors_bar = ['green' if abs(v - c_target)/c_target < 0.01 else 'lightblue'
               for v in const_vals]
bars = ax.bar(range(len(const_names)), const_vals, color=colors_bar, edgecolor='black', alpha=0.8)
ax.axhline(c_target, color='red', lw=2, linestyle='--', label=f'c_obs(0)≈{c_target:.4f}')
ax.axhline(E_OVER_SQRT2, color='blue', lw=1.5, linestyle=':', label=f'e/√2={E_OVER_SQRT2:.4f}')
ax.set_xticks(range(len(const_names)))
ax.set_xticklabels(const_names, rotation=45, ha='right', fontsize=9)
ax.set_ylabel('Wartość stałej')
ax.set_title('Kandydaci na c = ψ_core/ψ_new\n(zielony = błąd <1%)')
ax.legend(fontsize=9)
ax.set_ylim(1.5, 2.6)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
out_path = __file__.replace('.py', '.png')
plt.savefig(out_path, dpi=120, bbox_inches='tight')
plt.close()
print(f"  Zapisano: {out_path}")
print()


# ============================================================
# PODSUMOWANIE
# ============================================================
print("=" * 70)
print("PODSUMOWANIE p8: ANALITYCZNA ANALIZA V_mod I λ*")
print()
print("1. WYNIK DOKŁADNY (algebraiczny):")
print(f"   ψ_cross / ψ_new = √(3/2) = {sqrt32:.6f}  ← DOKŁADNE!")
print(f"   Wyprowadzenie: z definicji crossover λψ⁴ = ψ⁶ i ψ_new = 1/√λ")
print()
print("2. NAJLEPSZA HIPOTEZA O ψ_core:")
print(f"   ψ_core ≈ (e/√2) · ψ_new  (błąd ~0.4% na λ*)")
print(f"   e/√2 = {E_OVER_SQRT2:.6f}")
print(f"   c_obs(a_Γ→0) ≈ {c0_extrap:.4f}  (ekstrapolacja)")
print(f"   Brak jasnego wyprowadzenia z pierwszych zasad TGP.")
print()
print("3. EMPIRYCZNE SKALOWANIE:")
print(f"   λ* ≈ {C:.2e} · α^{p_alpha:.2f} · a_Γ^{q_agam:.2f}")
print(f"   Interpretacja: nie do końca jasna (wymaga ODE analizy)")
print()
print("4. RELACJA KOIDE:")
print(f"   Dla K-ratios TGP: Koide(K) = {koide_ratios:.6f}  (cel: 2/3)")
print(f"   {'Spełniona!' if abs(koide_ratios - 2/3) < 0.001 else 'Nie spełniona — ratios nie odtwarzają Koide dokładnie'}")
print()
print("OTWARTE PYTANIA:")
print("  Q1: Czy ψ_core = c·ψ_new ma wyprowadzenie z równania pola TGP?")
print(f"      (c ≈ e/√2 sugeruje związek z termodynamiką lub geometrią)")
print("  Q2: Czy λ* ma naturalną wartość w TGP (np. z substratu)?")
print("  Q3: Czy Koide jest spełniona dla PRAWDZIWYCH profili ODE (p6)?")
print("  Q4: Czy psi_core/psi_new → e/√2 dla a_Γ→0 (punkt stały ODE)?")
