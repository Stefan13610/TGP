"""
ex166_alpha1_vs_alpha2.py
==========================
Implikacje α=1 (substrate) vs α=2 (canonical) dla predykcji TGP.

KONTEKST (ex163):
  ODE substratowe (α=1, K∝g²): φ-FP daje r₂₁=206.77, m_τ z 83 ppm ✅
  ODE kanoniczne (α=2, K∝g⁴): bariera duchowa → φ-FP NIE DZIAŁA ❌

PYTANIE: Czy zmiana α=2→1 w solitonie psuje inne predykcje TGP?

PLAN:
  1. Operator kinetyczny: D_kin[Φ] z α=1 vs α=2
  2. Równanie pola TGP: co się zmienia?
  3. Grawitacja: PPN parametry (γ_PPN, β_PPN)
  4. Kosmologia: κ = 3/(4Φ₀), ciemna energia
  5. Masa bozonu przestrzennego m_sp
  6. Lensing: kąt odchylenia
  7. Synteza: co NAPRAWDĘ zależy od α, a co nie?
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np

print("=" * 72)
print("ex166: α=1 (substrate) vs α=2 (canonical) — implikacje dla TGP")
print("=" * 72)

# =====================================================================
# 1. Operator kinetyczny
# =====================================================================
print("\n--- 1. Operator kinetyczny D_kin ---")
print("""
  Ogólny operator:
    D_kin[φ] = ∇²φ + (α/φ)(∇φ)² = (1/(1+2α)) · φ^{-2α} ∇²(φ^{1+2α})

  α=2: D = ∇²φ + 2(∇φ)²/φ = (1/3φ²) ∇²(φ³)     [kanoniczny]
  α=1: D = ∇²φ + (1/φ)(∇φ)² = (1/φ²) ∇²(φ²)/2   [substratowy]
       = (1/(2φ²)) ∇²(φ²) — NIE! Policzmy:
       d/dr[φ²] = 2φφ', d²/dr²[φ²] = 2(φ')² + 2φφ''
       (1/(2φ²)) · ∇²(φ²) = (φ'' + (φ')²/φ) + (2/r)φ'  (w sferycznych)
       Hmm, nie upraszcza się tak prosto.

  Kluczowe pytanie: jakie φ^n czyni D_kin liniowym?
    D_kin[φ] = 0 ↔ ∇²(φ^{1+2α}) = 0
    α=2: ∇²(φ³) = 0 → φ³ jest harmoniczne
    α=1: ∇²(φ³) = 0? NIE. Dla α=1: φ^{1+2} = φ³ też!
    Hmm: 1+2α = 1+2·1 = 3 → ∇²(φ³) = 0 dla α=1 RÓWNIEŻ??
    Sprawdzenie: D = ∇²φ + (α/φ)(∇φ)²
    d²/dr² φ³ = d/dr(3φ²φ') = 6φ(φ')² + 3φ²φ''
    (1/3φ²)·[6φ(φ')² + 3φ²φ''] = 2(φ')²/φ + φ''
    → D_kin(α=2) = φ'' + 2(φ')²/φ = (1/3φ²)∇²(φ³) ✓

    Dla α=1: D(α=1) = φ'' + (φ')²/φ
    φ^{1+2·1} = φ³ → (1/3φ²)∇²(φ³) = φ'' + 2(φ')²/φ ≠ φ'' + (φ')²/φ
    → Różnica: (φ')²/φ !

  Poprawnie: D(α) = 0 ↔ ∇²(φ^{1+2α}) + rest = 0
  Ale to jest for radial 1D. W pełnym 3D:
  ∇²f + (2/r)f' → ∇²(φ^{1+2α}) = (1+2α)[(1+2α-1)φ^{2α-1}(φ')² + φ^{2α}∇²φ]
  = (1+2α)φ^{2α}[φ'' + 2α(φ')²/φ + (2/r)φ']
  = (1+2α)φ^{2α} · D_kin[φ]  (z α w D)
  → D_kin(α)[φ] = 0 ↔ ∇²(φ^{1+2α}) = 0 ✓ (dla każdego α)
""")

print("  Weryfikacja: ∇²(φ^n) = 0 z n = 1+2α")
print("  α=1: n=3 → φ³ harmoniczne")
print("  α=2: n=5 → φ⁵ harmoniczne")
print("  Oba mają tę samą strukturę — różnią się WYŁĄCZNIE w n.")
# Wait: let me recheck
# D(α) = ∇²φ + (α/φ)(∇φ)² + (2/r)φ' (spher.)
# ∇²(φ^n) = n(n-1)φ^{n-2}(φ')² + nφ^{n-1}(∇²φ + 2φ'/r)
# = nφ^{n-1}[∇²φ + (2/r)φ' + (n-1)(φ')²/φ]
# = nφ^{n-1} · [D(α=n-1)]
# So: ∇²(φ^n) = 0 ↔ D(α=n-1) = 0 → α = n-1
# α=1 → n=2, NOT n=3! Let me recompute:
# ∇²(φ²) = 2[(φ')² + φ∇²φ + 2φφ'/r]
# Hmm nah: ∇²(φ²) in 3D spherical:
# = 2φ∇²φ + 2(∇φ)² + 2·(2/r)·2φφ'... no.
# Flat 3D: ∇²(φ²) = 2φ∇²φ + 2|∇φ|²
# = 2φ[∇²φ + |∇φ|²/φ] = 2φ · D(α=1)
# So α=1 → ∇²(φ²) = 0 ↔ D(α=1) = 0 → n=2

print("""
  KOREKTA:
  ∇²(φⁿ) = nφ^{n-1}[∇²φ + (n-1)(∇φ)²/φ] (flat, without 2/r term)
  W sferycznym: ∇²(φⁿ) = nφ^{n-1}[∇²φ + (2/r)φ' + (n-1)(∇φ)²/φ]
                        = nφ^{n-1} · D_kin(α=n-1)[φ]

  → α=1: n=2 → ∇²(φ²) = 2φ · D(α=1) → φ² harmoniczne
  → α=2: n=3 → ∇²(φ³) = 3φ² · D(α=2) → φ³ harmoniczne ← to jest w TGP!

  WNIOSEK:
  α=2 → Φ³ harmoniczne → metryka g_µν = Φ^{2/3}η_µν → a = Φ^{1/3}
  α=1 → Φ² harmoniczne → metryka g_µν = Φ η_µν?     → a = Φ^{1/2}??

  Ale: emergentna metryka TGP to Φ^{2/3} (z α=2, sek03).
  Z α=1: metryka byłaby Φ^{p} z innym p!
""")

# =====================================================================
# 2. Metryka emergentna i kosmologia
# =====================================================================
print("=" * 72)
print("--- 2. Metryka emergentna: α=1 vs α=2 ---")
print("=" * 72)
print("""
  Z sek03 TGP: metryka emergentna g_µν = Φ^{2/(1+2α)} η_µν

  Wyprowadzenie: D_kin[Φ] = 0 → ∇²(Φ^{1+2α}) = 0
  Identyfikacja Φ^{1+2α} jako potencjału Newtona:
    Φ^{1+2α} = 1 - 2Φ_N/c²
  → g_00 = -(1 + 2Φ_N/c²) ≈ -(Φ/Φ₀)^{2(1+2α)}
  → g_ij = (Φ/Φ₀)^{2/(1+2α)} · δ_ij?

  Hmm, to wymaga dokładniejszego odczytania.
  Kluczowy krok: relacja Φ → g_µν.
  W TGP z α=2: g_µν ~ Φ^{2/3} η_µν (konforemna).
  Potęga 2/3 = 2/(1+2·2) = 2/5? NIE.

  Sprawdźmy z sek03 dokładnie.
  Relacja: ds² = Φ^{2p} dx² z p do ustalenia.
  Równanie pola: ∇²Φ + (α/Φ)(∇Φ)² = source
  W metryce Φ^{2p}: ∇²_{g}Φ = Φ^{-2p}[∇²Φ + ...correction terms]

  Alternatywne podejście: Φ³ harmoniczne → identyfikujemy z r:
  Φ³ ∝ 1/r (Newtonian) → Φ ∝ r^{-1/3} → g_µν ∝ Φ^{2/3} ∝ r^{-2/9}?
  Hmm, to nie jest standardowe.

  Dokładniej z TGP (sek03, konwencja):
  Stacjonarny rozkład daleko od źródła:
    ∇²(Φ^{1+2α}) = -q·ρ_mat
  Rozwiązanie: Φ^{1+2α} = Φ₀^{1+2α} + q·M/(4πr)
  → Φ ≈ Φ₀[1 + q·M/(4π(1+2α)Φ₀^{1+2α}·r)]

  PPN parametr γ_PPN: z porównania z GR.
  W GR: Φ_N = -GM/r
  W TGP: Φ^{1+2α} ≈ Φ₀^{1+2α}(1 + 2Φ_N/c²)
  gdzie Φ_N = -q·M/(4π·2·Φ₀^{1+2α}·r) · c²

  Lensing: γ_PPN = 1 w GR.
  W TGP: γ_PPN zależy od RELACJI między g_00 i g_ij.

  KLUCZOWY PUNKT: α wpływa na metrykę emergentną,
  a zatem na WSZYSTKIE predykcje grawitacyjne.
""")

# =====================================================================
# 3. Porównanie predykcji
# =====================================================================
print("=" * 72)
print("--- 3. Co zależy od α, a co nie? ---")
print("=" * 72)

Phi0 = 115.0

# κ = 3/(4Φ₀) — czy to zależy od α?
# κ pochodzi z porównania FRW TGP z FRW GR.
# W TGP: H² = (?)·ρ, gdzie (?) zależy od jak Φ wchodzi do metryki.
# Z α=2: a ∝ Φ^{1/3} → H = Φ̇/(3Φ) → H² ∝ Φ̇²/Φ² ∝ ρ/(Φ₀ · coś)
# → κ = 3/(4Φ₀) (z sek03)
#
# Z α=1: a ∝ Φ^{p} (inny p) → H = pΦ̇/Φ → H² ∝ p²Φ̇²/Φ² → κ(α=1) ≠ κ(α=2)?

print("""
  ┌────────────────────────────────────────────────────────────────────┐
  │  Predykcja          │  Zależy od α?  │  α=2 (kanon.) │  α=1 (sub.) │
  ├────────────────────────────────────────────────────────────────────┤
  │  φ-FP + r₂₁        │  TAK (ODE sol.)│  ❌ (ghost)   │  ✅ (207)   │
  │  Koide K=2/3        │  NIE (algeb.)  │  ✓ (0.667)   │  ✓ (0.667)  │
  │  K (kwarkowe)       │  NIE (algeb.)  │  ✓ (0.731)   │  ✓ (0.731)  │
  │  Metryka g_µν       │  TAK           │  Φ^{2/3}     │  Φ^{?}      │
  │  κ = 3/(4Φ_eff)     │  TAK           │  7/(2·115)≈0.030 │  ???     │
  │  N_e                │  TAK           │  (1/3)ln(1/ε₀)│ ???        │
  │  n_s                │  TAK           │  1-2/N_e     │  ???         │
  │  PPN (γ, β)         │  TAK           │  γ=1, β=1    │  ???         │
  │  m_sp               │  INDIRECTLY    │  √γ          │  √γ         │
  │  G_SM emergence     │  NIE           │  ✓            │  ✓           │
  └────────────────────────────────────────────────────────────────────┘
""")

# =====================================================================
# 4. Hipoteza dualizmu α
# =====================================================================
print("=" * 72)
print("--- 4. Hipoteza dualizmu: α_grav ≠ α_sol ---")
print("=" * 72)
print("""
  OBSERWACJA KLUCZOWA:
  • Soliton (cząstka) żyje wewnątrz ISTNIEJĄCEJ metryki (próżnia TGP).
  • Metryka (grawitacja) jest MAKROSKOPOWYM efektem pola Φ.

  → Mogą mieć RÓŻNE efektywne α!

  α_grav = 2: metryka emergentna g_µν = Φ^{2/3} η_µν
    (z prop:substrate-action, geometryczne sprzężenie K=φ⁴)
    → PPN, κ, n_s — poprawne (potwierdzone ex107, ex110)

  α_sol = 1: soliton w próżni (perturbacja wokół Φ₀)
    (z lem:K_phi2, perturbacyjne K_eff = Ja²φ̄²)
    → φ-FP, masy leptonów — poprawne (potwierdzone ex157)

  INTERPRETACJA:
  Soliton to perturbacja φ = φ̄ + δφ, z φ̄ ≈ 1 (próżnia).
  Efektywna K dla solitonu: K_eff(φ̄=1) z lem:K_phi2.
  Ale soliton exploruje g ∈ [g₀, 1], więc K(g) = g² jest traktowane
  jako K_sub dla PEŁNEGO profilu solitonu.

  Makroskopowo (kosmologia): Φ zmienia się od 0 do Φ₀.
  Tu pełne K(Φ) = K_geo·Φ⁴ jest poprawne → α=2.

  To jest jak RENORMALIZACJA:
  • UV (skala solitonu, r ~ 1/m_sp): efektywne α=1
  • IR (skala kosmologiczna, R ~ H⁻¹): efektywne α=2
""")

# =====================================================================
# 5. Test: czy α=1 dla solitonu jest spójne z α=2 dla grawitacji?
# =====================================================================
print("=" * 72)
print("--- 5. Spójność α_sol=1 z α_grav=2 ---")
print("=" * 72)
print("""
  Soliton „widzi" α_sol=1 bo operuje w reżimie:
    g(r) ∈ [g₀, 1], g₀ ∈ [0.8, 1.8]
    K(g) = g² (lem:K_phi2, poprawne perturbacyjnie wokół g=1)

  Grawitacja „widzi" α_grav=2 bo operuje w reżimie:
    Φ(x) od 0 do Φ₀ (pełen zakres)
    K(Φ/Φ₀) = (Φ/Φ₀)⁴ (prop:substrate-action, pełne nieliniowe)

  CZY TO SPÓJNE?
  TAK — jeśli K_sub(g) = g² jest PRZYBLIŻENIEM K_full = g⁴ wokół g=1:
    K_full(g) = g⁴ = [1 + (g-1)]⁴ ≈ 1 + 4(g-1) + 6(g-1)² + ...
    K_sub(g) = g² = [1 + (g-1)]² ≈ 1 + 2(g-1) + (g-1)² + ...

    Oba ≈ 1 + O(g-1) wokół g=1, ale współczynniki RÓŻNE.
    K_sub jest niższym przybliżeniem (mniej dokładny daleko od g=1).

  Ale: soliton elektronu ma g₀ = 0.87, mion g₀ = 1.41.
  Dla g₀ = 0.87: K_full = 0.87⁴ = 0.572, K_sub = 0.87² = 0.757
  Dla g₀ = 1.41: K_full = 1.41⁴ = 3.95,  K_sub = 1.41² = 1.99
  → Różnica ZNACZNA (nie perturbacyjna)!

  WNIOSEK: K_sub ≠ K_full w zakresie solitonu.
  To NIE jest przybliżenie — to są dwa RÓŻNE modele fizyczne.
""")

# Numeryczny test: wartości K(g) w zakresie solitonu
print("  Numeryczne porównanie K(g) w zakresie solitonu:")
print(f"  {'g':>6} {'K_sub=g²':>12} {'K_full=g⁴':>12} {'K_sub/K_full':>14}")
print("  " + "-" * 48)
for g in [0.5, 0.7, 0.87, 1.0, 1.2, 1.41, 1.73, 2.0]:
    Ks = g**2
    Kf = g**4
    ratio = Ks/Kf if Kf > 0 else 0
    print(f"  {g:6.2f} {Ks:12.4f} {Kf:12.4f} {ratio:14.4f}")

print("""
  K_sub/K_full = 1/g² → 1 tylko przy g=1.
  Dla g₀^e = 0.87: K_sub = 1.32·K_full (32% wyższy)
  Dla g₀^μ = 1.41: K_sub = 0.50·K_full (50% niższy)

  To oznacza: ODE substratowe i kanoniczne dają RÓŻNE A_tail(g₀)
  → RÓŻNE solitony → α=1 i α=2 to naprawdę DWIE TEORIE solitonów.
""")

# =====================================================================
# 6. Rozstrzygnięcie: co mówią dane?
# =====================================================================
print("=" * 72)
print("--- 6. Rozstrzygnięcie empiryczne ---")
print("=" * 72)

from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
R_MAX = 120.0
M_E, M_MU, M_TAU = 0.511, 105.658, 1776.86

def solve_ode(g0, alpha):
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-8)
        if alpha == 1:
            src = 1 - g
        elif alpha == 2:
            src = g**2 * (1 - g)  # canonical source
        else:
            src = 1 - g  # general substrate-like
        cross = (alpha / g) * gp**2
        if r < 1e-10:
            return [gp, (src - cross) / 3.0]
        return [gp, src - cross - 2*gp/r]
    try:
        sol = solve_ivp(rhs, (0, R_MAX), [g0, 0], rtol=1e-10, atol=1e-12, max_step=0.05)
        return sol.t, sol.y[0]
    except:
        return np.array([0]), np.array([g0])

def A_of_g0(g0, alpha):
    r, g = solve_ode(g0, alpha)
    mask = (r >= 25) & (r <= 90)
    if np.sum(mask) < 30: return 0.0
    rf = r[mask]; df = (g[mask]-1)*rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return np.sqrt(bc[0]**2 + bc[1]**2)

# Test: φ-FP dla α=1 z różnymi source terms
print("\n  Test: α=1 z source = 1-g (substratowy) vs α=1 z source = g²(1-g)")
print(f"  {'source':>20} {'g₀^e':>10} {'g₀^μ':>10} {'r₂₁':>10} {'status':>10}")
print("  " + "-" * 65)

for source_label in ["1-g (substrate)", "g²(1-g) (hybrid)"]:
    def solve_custom(g0, src_type=source_label):
        def rhs(r, y):
            g, gp = y
            g = max(g, 1e-8)
            if "substrate" in src_type:
                src = 1 - g
            else:
                src = g**2 * (1 - g)
            cross = (1.0 / g) * gp**2  # α=1
            if r < 1e-10:
                return [gp, (src - cross) / 3.0]
            return [gp, src - cross - 2*gp/r]
        sol = solve_ivp(rhs, (0, R_MAX), [g0, 0], rtol=1e-10, atol=1e-12, max_step=0.05)
        return sol.t, sol.y[0]

    def A_custom(g0, src_type=source_label):
        r, g = solve_custom(g0, src_type)
        mask = (r >= 25) & (r <= 90)
        if np.sum(mask) < 30: return 0.0
        rf = r[mask]; df = (g[mask]-1)*rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc = np.linalg.lstsq(M, df, rcond=None)[0]
        return np.sqrt(bc[0]**2 + bc[1]**2)

    r21_PDG = M_MU / M_E
    g0_scan = np.linspace(0.4, 1.1, 50)
    g0_e = None

    def r21_res(g0):
        a1 = A_custom(g0)
        a2 = A_custom(PHI*g0)
        if a1 < 1e-10: return 1e6
        return (a2/a1)**4 - r21_PDG

    for i in range(len(g0_scan)-1):
        if PHI*g0_scan[i+1] > 2.5: continue
        try:
            ra = r21_res(g0_scan[i])
            rb = r21_res(g0_scan[i+1])
            if ra * rb < 0:
                g0_e = brentq(r21_res, g0_scan[i], g0_scan[i+1], xtol=1e-8)
                break
        except: pass

    if g0_e is not None:
        A_e = A_custom(g0_e)
        A_mu = A_custom(PHI*g0_e)
        r21_calc = (A_mu/A_e)**4
        print(f"  {source_label:>20} {g0_e:10.6f} {PHI*g0_e:10.6f} {r21_calc:10.2f} {'✅' if abs(r21_calc - r21_PDG) < 1 else '❌':>10}")
    else:
        print(f"  {source_label:>20} {'---':>10} {'---':>10} {'---':>10} {'❌':>10}")

# =====================================================================
# 7. WNIOSKI
# =====================================================================
print(f"\n{'='*72}")
print("WNIOSKI ex166")
print(f"{'='*72}")
print(f"""
  1. α WPŁYWA NA:
     - Metrykę emergentną: g_µν = Φ^{{2/(1+2α)}} η_µν
       α=2: Φ^{{2/5}}?? → TRZEBA sprawdzić z sek03 dokładnie
       α=1: Φ^{{2/3}}?? → to zależy od wyprowadzenia w sek03
     - Relację D_kin = 0 → ∇²(Φ^{{1+2α}}) = 0
     - ODE solitonowe i A_tail(g₀)
     - Wszystkie predykcje kosmologiczne (κ, N_e, n_s)

  2. α NIE WPŁYWA NA:
     - Koide K (algebraiczne, nie zależy od ODE)
     - Emergencję G_SM (topologiczna)
     - m_sp = √γ (zależy od V, nie od K)

  3. HIPOTEZA DUALIZMU α:
     α_grav = 2 (kosmologia, PPN): z prop:substrate-action (K=φ⁴)
     α_sol = 1 (soliton, masy): z lem:K_phi2 (K=φ²)
     NIE SĄ przybliżeniem — to dwa RÓŻNE reżimy!
     K_sub/K_full = 1/g² → 32% różnica przy g₀^e.

  4. PROBLEM:
     Jeśli α_grav ≠ α_sol, TGP ma de facto DWA parametry α,
     co osłabia jedność teorii.
     ALBO: α jest EFEKTYWNE i interpoluje między 1 (UV) a 2 (IR).
     → Renormalizacja α(scale) — do zbadania.

  5. ALTERNATYWA (OCKHAM):
     Może α=1 jest fundamentalny WSZĘDZIE?
     → Trzeba przeliczyć PPN, κ, n_s z α=1.
     → Jeśli OK → teoria prostsza (K=φ², α=1 jedyny)
     → Jeśli NIE → dualizm α niezbędny.

  STATUS: OTWARTY — wymaga przeliczenia grawitacji z α=1.
  REKOMENDACJA: ex167 — PPN parametry i κ z α=1.
""")
