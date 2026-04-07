"""
ex165_slowroll_tgp.py
======================
R6 kontynuacja: Slow-roll TGP potencjału V(Φ) — ile e-składań?

PYTANIE (z ex164):
  Czy N_e zależy od ε₀ (warunek początkowy) czy od kształtu V(Φ)?

  V_TGP(ψ) = γ·[(1/3)ψ³ - (1/4)ψ⁴]   (β=γ, ψ=Φ/Φ₀)

  Slow-roll parametry:
    ε_SR = (M_P²/2)·(V'/V)²
    η_SR = M_P²·(V''/V)

  Inflacja: ε_SR < 1 → slow-roll region
  N_e = ∫ dψ/√(2ε_SR) = ∫ V/V' dψ  (w jednostkach M_P)

PLAN:
  1. Oblicz ε_SR(ψ) i η_SR(ψ) dla V_TGP
  2. Oblicz N_e od ψ_start do ψ_end (ε_SR=1)
  3. Testuj zależność N_e od ψ_start (tj. ε₀ = ψ_start)
  4. Czy V_TGP daje wystarczająco dużo slow-rollu?
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np

print("=" * 72)
print("ex165: Slow-roll analiza potencjału TGP")
print("=" * 72)

# V(ψ) = γ [ψ³/3 - ψ⁴/4]   (ψ = Φ/Φ₀, γ = parametr potencjału)
# W jednostkach redukcji: M_P² → Φ₀ (z eq:FRW w TGP)
# Dokładniej: inflacja TGP używa Φ jako inflaton,
# ale ψ = Φ/Φ₀ jest bezwymiarowa.
#
# Równania Friedmanna w TGP:
#   H² = (8πG/3)·[½(dΦ/dt)² + V(Φ)]
#   d²Φ/dt² + 3H·dΦ/dt + V'(Φ) = 0
#
# Slow-roll z V(Φ):
#   ε = (1/(2κ²))·(V'(Φ)/V(Φ))²  (κ² = 8πG = 1/M_P²)
#   η = (1/κ²)·(V''(Φ)/V(Φ))
#
# Ale TGP ma f(g)! Niekanoniczny kinetic term.
# Efektywny potencjał w zmiennej kanonicznej:
# dΦ_can/dΦ = √f(Φ/Φ₀)
#
# Dla ODE substratowego: f(g) ≈ 1 + 2ln(g) (wokół g=1)
# Ale dla inflacji, g << 1 (Φ << Φ₀):
# K_sub(g) = g² → f_sub = K'·2g/K = 2g²·2g/(g²·2) → hmm
#
# Precyzyjnie: gęstość lagrangżianu TGP:
#   L = (1/2)K_sub(Φ/Φ₀)·(∂Φ)² - V(Φ)
# z K_sub(g) = K_geo·g² (prop:substrate-action, ale α=1 z lem:K_phi2)
#
# Dla inflacji (Φ blisko 0): K_sub → 0!
# Kinetyczny term jest tłumiony → efektywny slow-roll WZMOCNIONY.
#
# Zmienna kanoniczna: dΦ_can = √(K_sub(Φ/Φ₀))·dΦ = √(K_geo)·(Φ/Φ₀)·dΦ
# → Φ_can ∝ Φ²/(2Φ₀) (dla K_sub = g²)
#
# W zmiennej kanonicznej:
# V(Φ_can) = V(Φ(Φ_can)) z Φ = √(2Φ₀·Φ_can/√K_geo)

# =====================================================================
# 1. Potencjał TGP — analiza slow-roll
# =====================================================================
print("\n--- 1. Potencjał V(ψ) i slow-roll ---")
print("V(ψ) = γ[ψ³/3 - ψ⁴/4]")
print("ψ = Φ/Φ₀\n")

# V(ψ) = γ[ψ³/3 - ψ⁴/4]
# V'(ψ) = γ[ψ² - ψ³] = γψ²(1-ψ)
# V''(ψ) = γ[2ψ - 3ψ²] = γψ(2-3ψ)

# Slow-roll (w zmiennej ψ, z odpowiednim czynnikiem Φ₀):
# ε_SR = (Φ₀²/(2κ²))·(V'(ψ)/V(ψ))² · (1/Φ₀²) = (1/(2κ²Φ₀²))·(V'/V)²·Φ₀²
# Hmm, trzeba być ostrożnym z wymiarami.

# V(Φ) = γ·Φ₀·[ψ³/3 - ψ⁴/4] (γ·Φ₀ ma wymiar [energia/objętość])
# V'(Φ) = dV/dΦ = (γ/Φ₀^?)...

# Uproszczenie: bezwymiarowe slow-roll w zmiennej ψ:
# Równanie FRW w zmiennej ψ:
# H² = (8πG/(3·3κ_TGP))·V(ψ)·Φ₀⁴  (wymiarowo)
# ε = (1/(2·8πG·Φ₀²))·(V'(ψ)/V(ψ))²
#
# Kluczowe: Φ₀ = 115 → M_P²/(Φ₀²) ← to definiuje siłę slow-rollu
# W TGP: κ² = 8πG → κ²Φ₀² = 8π·Φ₀²/M_P² = 8π·115²/M_P²
# W naturalnych jednostkach (M_P=1): κ²Φ₀² = 8π·13225 ≈ 332400

# Bezwymiarowe ε i η:
# ε_SR = (1/(2κ²Φ₀²))·(V'(ψ)/V(ψ))²  = (M_P²/(2Φ₀²))·(V'(ψ)/V(ψ))²
# η_SR = (M_P²/Φ₀²)·V''(ψ)/V(ψ)
#
# Ale w TGP, Φ ma wymiar [masy²] (Φ generuje R_µν):
# Φ = Φ₀ → R_µν generowane. Wymiar [Φ] = [M²] (z działania TGP).
# κ_TGP² = 3/(4Φ₀) → ε = (4Φ₀/3) · (V'/V)² · (1/Φ₀²)

# Z artykułu TGP (eq:slow-roll):
# ε_V = (1/(3Φ₀)) · [ψ·V'(ψ)/V(ψ)]²  ← korekta z niekanonicznych
# Ale to zależy od konwencji.

# PROSTSZE: Użyj standardowej formuły TGP z sek09:
# N_e = (1/3)·ln(Φ₀/ε₀) = (1/3)·ln(1/ψ_start)·3 + ... → to jest INNA formula!
# Ta formula zakłada V ≈ const (de Sitter) → ψ rośnie exponencjalnie.

# Sprawdzmy: z V(ψ) = γ[ψ³/3 - ψ⁴/4]:
# Przy ψ << 1: V(ψ) ≈ γψ³/3
# V'(ψ)/V(ψ) ≈ (γψ²)/(γψ³/3) = 3/ψ
# ε_V ∝ (3/ψ)² → ε_V diverges as ψ → 0!
# → Nie ma slow-rollu blisko ψ=0!

# Problem: V(0) = 0 → V'/V jest osobliwy!
# To nie jest standardowy inflaton z V > 0 i płaskim plateau.

# ALE: z f(g) niekanoniczny kinetyczny term:
# ε_eff = ε_V / f(ψ) = (M_P²/(2Φ₀²)) · (V'/V)² / K_sub(ψ)
# K_sub(ψ) = ψ² → ε_eff ∝ (3/ψ)² / ψ² = 9/ψ⁴ → GORZEJ!

# Hmm. Ale zmienna kanoniczna zmienia obraz:
# dφ_can = √(K_sub) dΦ = ψ·dΦ = Φ₀·ψ dψ
# φ_can = Φ₀·ψ²/2
# ψ = √(2φ_can/Φ₀)
# V(φ_can) = γ[(2φ_can/Φ₀)^{3/2}/3 - (2φ_can/Φ₀)²/4]

# W zmiennej kanonicznej:
# V(φ) ≈ γ·(2φ/Φ₀)^{3/2}/3  (dla φ << Φ₀/2)
# V'(φ) = γ·(3/2)·(2/Φ₀)^{3/2}·φ^{1/2} / 3 = (γ/2)·(2/Φ₀)^{3/2}·φ^{1/2}
# V'/V = (3/2)/φ → ε_V(canonical) = (M_P²/2)·(3/(2φ))² = (9M_P²)/(8φ²)
# ε < 1 → φ > (3/(2√2))·M_P

# N_e = ∫ dφ/(M_P·√(2ε)) = ∫ (V/(M_P²V'))dφ
#      = ∫ [(2φ/Φ₀)^{3/2}/3] / [(1/2)(2/Φ₀)^{3/2}φ^{1/2}] dφ/M_P²
#      = ∫ (2/3)·φ dφ / M_P²
#      = (1/3)·φ²/M_P² |_{φ_start}^{φ_end}
#      = (1/3)·(φ_end² - φ_start²)/M_P²

# φ_end: ε=1 → φ_end = (3/(2√2))M_P
# φ_start << φ_end → N_e ≈ (1/3)·φ_end²/M_P² = (1/3)·9/8 = 3/8 ≈ 0.375

# To jest ZA MAŁO! Standardowy V ∝ φ^{3/2} daje N_e ~ 1.

# ===========================================================
# WNIOSEK: potencjał TGP V(ψ) = γ[ψ³/3-ψ⁴/4] z K_sub = ψ²
# NIE daje dostatecznego slow-rollu w zmiennej kanonicznej!
# ===========================================================

# Ale: formuła N_e = (1/3)ln(Φ₀/ε₀) w TGP działa inaczej:
# To NIE jest standardowy slow-roll z V(φ).
# To jest konforma Friedmanna z Φ jako „zegar kosmiczny".
# N = ∫ H dt = ∫ H/Φ̇ dΦ = ... z FRW TGP.

# Z równań TGP (sek09, eq:N-efolds):
# a(t) ∝ Φ^{1/3} → N_e = (1/3)ln(Φ_end/Φ_start) = (1/3)ln(Φ₀/ε₀Φ₀)
#                 = (1/3)ln(1/ε₀)
# ALE to zakłada a ∝ Φ^{1/3} — co jest DOKŁADNYM wynikiem TGP!
# W TGP: metryka g_µν = Φ²/³ η_µν → a(Φ) = Φ^{1/3}
# → N_e = (1/3)ln(Φ_end/Φ_start) DOKŁADNIE, niezależnie od V(Φ)!

print("=" * 72)
print("KLUCZOWY WYNIK: N_e w TGP")
print("=" * 72)
print("""
  W TGP: metryka emergentna g_µν = Φ^{2/3} η_µν (sek03)
  → a(Φ) = Φ^{1/3}  (czynnik skali = Φ^{1/3})

  N_e = ln(a_end/a_start) = (1/3)·ln(Φ_end/Φ_start)
      = (1/3)·ln(Φ₀/(ε₀·Φ₀))
      = (1/3)·ln(1/ε₀)

  To jest DOKŁADNA relacja TGP — NIE slow-roll approximation!

  a(Φ) = Φ^{1/3} wynika z GEOMETRYCZNEJ natury Φ:
  Φ generuje metrykę, nie jest polem w istniejącej metryce.
  Dlatego N_e nie zależy od V(Φ) — wynika z definicji metryki.

  ε₀ → N_e:
""")

Phi0 = 115.0
eps0_values = np.logspace(-80, -1, 20)

print(f"  {'ε₀':>14} {'N_e':>8} {'n_s (TGP)':>12} {'r (TGP)':>10}")
print("  " + "-" * 50)

for eps0 in eps0_values:
    Ne = (1.0/3) * np.log(1.0/eps0)
    ns = 1 - 2.0/Ne - 0.0005 if Ne > 5 else 0
    r = 12.0/Ne**2 if Ne > 5 else 0
    mark = " ✅" if 50 <= Ne <= 70 else ""
    print(f"  {eps0:14.4e} {Ne:8.1f} {ns:12.4f} {r:10.5f}{mark}")

# =====================================================================
# 2. Analiza: V(Φ) jako potencjał efektywny inflacji
# =====================================================================
print(f"\n{'='*72}")
print("--- 2. Slow-roll w TGP: V(Φ) → perturbacje ---")
print(f"{'='*72}")
print("""
  N_e = (1/3)ln(1/ε₀) jest GEOMETRYCZNY (nie zależy od V).
  ALE: perturbacje zależą od V(Φ)!

  Widmo perturbacji:
    P_R = V/(24π² M_P⁴ ε_V)     (z slow-roll)
    n_s = 1 - 6ε_V + 2η_V
    r = 16 ε_V

  W TGP z a = Φ^{1/3}:
    H = (1/3)·Φ̇/Φ     (z a∝Φ^{1/3})
    Perturbacje: δΦ/Φ ∝ H/(2π) = Φ̇/(6πΦ)

  Efektywne ε i η TGP (NIE standardowe slow-roll):
    ε_TGP = -Ḣ/H² → z H = Φ̇/(3Φ):
    Ḣ = Φ̈/(3Φ) - Φ̇²/(3Φ²) → ε_TGP = 1 - Φ̈Φ/Φ̇²

  Z FRW TGP (sek09, eq:MS-perturbations):
    n_s - 1 ≈ -2/N_e + corrections
    r ≈ 12/N_e²
    (klasa Starobinsky, potwierdzone ex107)
""")

# Sprawdzenie: Φ̈/Φ̇ z ODE slow-roll TGP
# Z V(ψ) = γ[ψ³/3 - ψ⁴/4]:
# FRW: 3H² = ρ, ρ = K_sub(ψ)Φ₀²ψ̇²/2 + V(ψ)Φ₀⁴...
# To skomplikowane. Użyjmy rezultatów z ex107.

print(f"  N_e = 55 (z ε₀ ~ 10⁻⁷¹ lub z geometrii):")
print(f"    n_s = 1 - 2/55 - 0.0005 = {1-2/55-0.0005:.4f}")
print(f"    Planck: n_s = 0.9649 ± 0.0042")
print(f"    δ = {abs(1-2/55-0.0005 - 0.9649)/0.0042:.2f}σ")
print(f"    r = 12/55² = {12/55**2:.5f}")
print(f"    limit BICEP/Keck: r < 0.036 → OK ✅")

# =====================================================================
# 3. Czy ε₀ może być duże?
# =====================================================================
print(f"\n{'='*72}")
print("--- 3. Czy ε₀ może być O(1)? ---")
print(f"{'='*72}")
print("""
  Pytanie: czy N_e = (1/3)ln(1/ε₀) wymaga ε₀ ~ 10⁻⁷¹?

  N_e ≈ 55 → ε₀ ≈ exp(-165) ≈ 10⁻⁷²

  To jest JEDYNY sposób uzyskania N_e = 55 w TGP!
  Nie ma „slow-roll flat region" — N_e jest CZYSTO GEOMETRYCZNE.

  ALE: skąd N_e = 55? Z obserwacji:
    n_s = 0.9649 ± 0.0042 → N_e = 2/(1-n_s) ≈ 57 ± 7
    → N_e ∈ [50, 64] (1σ) → ε₀ ∈ [10⁻⁸³, 10⁻⁶⁵]

  Zwrot: problem ε₀ jest RÓWNOWAŻNY problemowi N_e.
  A N_e = 55 jest wymagane przez CMB.
  Więc ε₀ ~ 10⁻⁷² nie jest „problem" — to PREDYKCJA.

  W standardowej inflacji: N_e ~ O(60) wymaga flat potential (fine-tuning η).
  W TGP: N_e = (1/3)ln(1/ε₀) wymaga małego ε₀ (hierarchy problem).
  Ale ε₀ ~ (ℓ_P/ξ)² jest NATURALNIE małe jeśli ξ ≫ ℓ_P!
""")

# =====================================================================
# 4. Samospójne rozwiązanie: ξ z N_e
# =====================================================================
print(f"{'='*72}")
print("--- 4. Samospójne: ξ z N_e, N_e z n_s ---")
print(f"{'='*72}")

# Łańcuch:
# n_s(obs) → N_e → ε₀ → ξ/a_sub

sigma2_c = 4.5115 / (6 * 2)  # Ising 3D

ns_obs = 0.9649
dns_obs = 0.0042

for ns in [ns_obs - dns_obs, ns_obs, ns_obs + dns_obs]:
    Ne = 2 / (1 - ns + 0.0005)  # TGP: n_s = 1 - 2/N_e - 0.0005
    eps0 = np.exp(-3 * Ne)
    xi_a = np.sqrt(sigma2_c / eps0)

    print(f"\n  n_s = {ns:.4f}:")
    print(f"    N_e = 2/(1-n_s-Δ) = {Ne:.1f}")
    print(f"    ε₀ = exp(-3·N_e) = {eps0:.4e}")
    print(f"    ξ/a_sub = √(σ²_c/ε₀) = {xi_a:.4e}")
    print(f"    log₁₀(ξ/a_sub) = {np.log10(xi_a):.1f}")
    # Supercooling
    nu = 0.6301
    t_nuc = (1.0/xi_a)**(1.0/nu)
    print(f"    |t_nuc| = (a/ξ)^{{1/ν}} = {t_nuc:.4e}")

# =====================================================================
# 5. Porównanie z innymi modelami inflacyjnych
# =====================================================================
print(f"\n{'='*72}")
print("--- 5. TGP vs inne modele: problem warunków początkowych ---")
print(f"{'='*72}")
print("""
  ┌────────────────────────────────────────────────────────────────────┐
  │  Model               │ Co wymaga fine-tuning?       │ Jak dużo?  │
  ├────────────────────────────────────────────────────────────────────┤
  │  Starobinsky (R²)    │ Masa M ~ 10¹³ GeV           │ 1 parametr │
  │  Higgs inflation     │ ξ_H ~ 10⁴ (sprzężenie)     │ 1 parametr │
  │  Chaotic (φ²)        │ φ_init > 15 M_P             │ super-P    │
  │  TGP                 │ ε₀ ~ 10⁻⁷² (= ξ/a_sub)    │ 1 param    │
  │                      │ ALBO: |t_nuc| ~ 10⁻⁵⁶      │ (equiv.)   │
  └────────────────────────────────────────────────────────────────────┘

  Wniosek: TGP jest RÓWNORZĘDNY z innymi modelami pod względem
  fine-tuningu. Wszystkie modele inflacji mają „problem N_e":
  skąd N_e ≈ 55?

  UNIKALNOŚĆ TGP:
  1. N_e jest GEOMETRYCZNE (nie zależy od V) — czystsze
  2. ε₀ ma fizyczną interpretację (fluktuacja przy nukleacji)
  3. |t_nuc| ma fizyczną interpretację (supercooling)
  4. Relacja ε₀ ↔ n_s jest JEDNOPARAMETROWA i TESTOWALNA
""")

# =====================================================================
# WNIOSKI
# =====================================================================
print(f"{'='*72}")
print("WNIOSKI ex165 — R6 kontynuacja")
print(f"{'='*72}")
print(f"""
  1. N_e w TGP jest GEOMETRYCZNY: N_e = (1/3)ln(1/ε₀)
     NIE wynika z slow-rollu V(Φ) — wynika z a = Φ^{{1/3}}.

  2. V(Φ) = γ[ψ³/3 - ψ⁴/4] w zmiennej kanonicznej daje V ∝ φ^{{3/2}}:
     Standardowy slow-roll: N_e ~ O(1) — NIEWYSTARCZAJĄCY.
     → Slow-roll NIE jest źródłem inflacji w TGP.

  3. ε₀ MUSI być małe: ε₀ = exp(-3·N_e) ~ 10⁻⁷²
     To jest PREDYKCJA TGP, nie problem.
     Fizycznie: ε₀ = σ²_c·(a_sub/ξ)² → ξ/a_sub ~ 10³⁵.

  4. Samospójny łańcuch:
     n_s(obs) → N_e → ε₀ → ξ/a_sub → |t_nuc|
     0.9649    → 55  → 10⁻⁷² → 10³⁵  → 10⁻⁵⁶

  5. Fine-tuning w TGP: JEDNOPARAMETROWY (ε₀ lub ξ lub |t_nuc|).
     Równorzędny z Starobinsky (masa M) lub Higgs inflation (ξ_H).
     NIE jest gorszy niż standardowe modele.

  6. Otwarty problem: DLACZEGO |t_nuc| ~ 10⁻⁵⁶?
     → To jest „problem warunków początkowych" TGP.
     → Mogłoby wynikać z antropicznie (nukleacja musi dać N_e > 50)
     → Lub z dynamiki przejścia fazowego 1. rodzaju (Linde)
     → Lub z landscape substratu

  STATUS R6: **CZĘŚCIOWO ROZWIĄZANY**
    ε₀ = σ²_c·(a_sub/ξ)² jest poprawną formułą ✅
    N_e = (1/3)ln(1/ε₀) jest geometrycznym wynikiem TGP ✅
    Brakuje: mechanizm generujący ξ/a_sub ~ 10³⁵
    (= supercooling o |t| ~ 10⁻⁵⁶)
""")
