---
title: "UV.3 program — explicit Φ₀ wave-function renormalization Z_Φ = 14/3"
date: 2026-05-02
cycle: UV.3
status: PROPOSED
parent: "[[../op-uv2-mtgp-absolute-scale/CRITIQUE_repackaged_circularity_2026-05-02.md]]"
predecessors:
  - "[[../op-uv-as-ngfp/Phase3_results.md]]"            # UV.1 NGFP {g*, λ*, η_N*}
  - "[[../op-gamma1-phi-eff-anchor-resolution/README.md]]" # γ.1 Φ_eff = 8π identification
  - "[[../op-lambda1-e2-amplitude-emergence/]]"         # λ.1 (10/3)e² connection
  - "[[../../core/sek00_summary/sek00_summary.tex]]"    # canonical Φ₀^bare / Φ_eff table
  - "[[../../core/formalizm/dodatekN_erg_renormalizacja.tex]]" # ERG flow
tags:
  - TGP
  - UV3
  - wave-function-renormalization
  - Z_Phi
  - 14/3
  - 4.6667
  - bare-vs-effective
  - core-anchored
  - PROPOSED
---

# UV.3 — jawna nazwa renormalizacji Φ₀: **Z_Φ = 14/3 ≈ 4.667**

> **Cel:** zastąpić cyrkularny program UV.2 (`K_struct = N_A·2π² ≈ 173`,
> post-hoc fit do `M_TGP/M_GUT`) jawnym sformułowaniem renormalizacji UV
> która **już jest** w rdzeniu TGP, ale nie była dotąd nazwana. Sek00
> tabela (linie 385–388) podaje:
>
> $$\Phi_0^{\rm bare} = 168\,\Omega_\Lambda \approx 115, \quad
> \Phi_{\rm eff} = \Phi_0^{\rm bare}\cdot\tfrac{3}{14} \approx 24{,}66, \quad
> \boxed{\;Z_\Phi \equiv \frac{\Phi_0^{\rm bare}}{\Phi_{\rm eff}} = \frac{14}{3} \approx 4{,}6667\;}$$
>
> To jest **TGP wave-function/screening renormalization** (UV → IR). Nie
> wymaga `M_GUT` jako anchora i nie wprowadza nowych dim-less kombinacji.
> Wartość `14/3` wynika **algebraicznie** z definicji potencjału akcji:
> $P(1)/V(1) = (\gamma/56)/(\gamma/12) = 12/56 = 3/14$ (eq. 74 sek00).

## Hypothesis (UV.3-H)

W TGP istnieje jeden i tylko jeden czynnik renormalizacji UV→IR dla pola
substratu Φ. Nazywamy go `Z_Φ` (wave-function renormalization w
ekranowanym dielektryku) i jego wartość jest zafiksowana strukturalnie:

$$Z_\Phi = \frac{P(1)}{V(1)}\bigg|^{-1} = \frac{14}{3} = 4{,}6\overline{6}$$

Wszystkie pozorne "różne wartości Φ₀" w rdzeniu (≈ 115 vs ≈ 24.66 vs ≈ 25
vs 24.783) są dwiema warstwami: **Φ₀^bare** (poziom UV — kosmologiczna
kalibracja Warstwy II, `168·Ω_Λ`) i **Φ_eff** (poziom IR — ekranowany
dielektryk, używany w `α_em`, `α_s`, `κ`, `a_Γ`). Mapa między nimi to
**wyłącznie** stała `Z_Φ = 14/3`.

## Reference frame — różnica vs UV.2 (krytyka)

| element | UV.2 (BLOCKING, post-hoc fit) | UV.3 (anchored w rdzeniu) |
|---|---|---|
| obiekt renormalizacji | M_TGP (skala masowa) | Φ₀ (wartość pola substratu) |
| anchor zewnętrzny | M_GUT_2loop ≈ 2·10¹⁶ GeV (band ±20 %) | Ω_Λ (Planck, ±0.7 %) ×16 8 |
| dim-less czynnik | K_struct = N_A·2π² ≈ 173 (4-cand fit) | Z_Φ = 14/3 ≈ 4.667 (algebraic) |
| pochodzenie | "S³ vol × photon-ring" narracja | P(1)/V(1) = (γ/56)/(γ/12) **strict** |
| status falsyfikowalności | nieoperacyjny w paśmie M_GUT 10–30 % | operacyjny: drift Z_Φ vs 14/3 ≤ 0.1 % |
| niezależne predykcje | brak (wszystko propaguje 0.30 %) | Φ_eff·Z_Φ = `168·Ω_Λ` ZGADZA SIĘ z `≈115` (sek00:385) |

UV.2 fittował **liczbę 173** żeby pasować do `(M_Pl/M_GUT)·√(g*/N_A)`.
UV.3 **deklaruje 14/3** z definicji `P/V` i sprawdza czy istniejące
wartości w rdzeniu pasują pod **jeden** czynnik renormalizacji.

## Phase plan (5 + 6 + 5 = 16 sub-tests, anti-circularity gate Phase 0)

### Phase 0 — anti-circularity audit (1 obowiązkowa kontrola)

- **U0.1 Inputs–Outputs balance sheet**: każdy sub-test musi jawnie podać
  (a) external inputs (Planck Ω_Λ, PDG α_s/α_em, sek00 P/V definicje)
  (b) structural axioms (β = γ, K(g) = g⁴, P(g) = β g⁷/7 − γ g⁸/8)
  (c) derived outputs (Φ₀^bare, Φ_eff, Z_Φ, κ).
  
  Tautology test: **żaden sub-test nie może** zwracać wartości
  algebraicznie kasującej się z wejściem. Falsifiability test: dla
  każdego output istnieje wartość input która by go odrzuciła.

### Phase 1 — inwentaryzacja Φ₀ + nazwanie Z_Φ (5)

- **U1.1 Catalog wszystkich wartości Φ₀/Φ_eff w rdzeniu** (grep + ledger):
  - `Φ₀^bare ≈ 115` (sek00:385, "168·Ω_Λ", Warstwa II)
  - `Φ_eff = Φ₀·3/14 ≈ 24.66` (sek00:386, eq. 77)
  - `Φ_eff = 24.65` (sek08:445, "bare, z Λ_obs" — terminologia inna!)
  - `Φ_eff = 24.783` (dodatekV:130, Brannen α_s lock)
  - `Φ_eff = 8π ≈ 25.13` (γ.1 P3.5, T-Λ structural pure)
  - `Φ_eff = (10/3)e² ≈ 24.63` (γ.1 P3.5+, T-Λ corrected)
  - `Φ₀(Λ) ≈ 23.3` (dodatekQ Q.4, 96πG₀ρ_Λ/H₀²)
  - `Φ₀(r₂₁) ≈ 25.4` (dodatekQ Q.4, (α_K √r₂₁)^(3/5))
  - `Φ₀(κ) ≈ 24.7` (dodatekQ Q.4, 3/(4κ_obs))
  
  **Gate:** wszystkie wartości w jednym z **dwóch** klastrów
  (UV: ≈ 115; IR: ≈ 24.6–25.1). Pas IR ma rozrzut ≤ 2 %.

- **U1.2 Algebraiczna derywacja czynnika 3/14 z P/V**:
  Z sek00 eq. 74: `P(1) = γ/56`, `V(1) = γ/12`. Stąd
  `P(1)/V(1) = 12/56 = 3/14` (sympy exact, drift 0.0000 %).
  **Niezależne** od `β = γ` warunku w sensie wartości liczbowej (pojawia
  się wyłącznie z eksponentów 7 i 8 w `P` oraz 3 i 4 w `V`).

- **U1.3 Definicja Z_Φ ≡ Φ₀^bare/Φ_eff = 14/3**:
  Sympy LOCK. **Drift target:** 0.0 % (definicja).

- **U1.4 Cross-check Z_Φ vs istniejący ledger**:
  - `Φ₀^bare = 168·Ω_Λ_Planck = 168·0.6847 = 115.03`
  - `Φ_eff^cosmo = 36·Ω_Λ = 24.65`
  - `(Φ₀^bare)/(Φ_eff^cosmo) = 115.03/24.65 = 4.667` = `14/3`. **Drift target:** ≤ 0.1 %.

- **U1.5 Anti-circularity check (Phase 0 obowiązkowy raport)**:
  pokazać że Z_Φ = 14/3 NIE jest output sympy zwracającym własne wejście:
  - input: definicje `P(g), V(g)` (struktura akcji, sek08a)
  - output: liczba `14/3` z `P(1)/V(1)`
  - test: czy zmiana eksponentów `7 → 9` w P daje `Z_Φ = (1/9 − 1/10)·12 = 12/90 = 2/15`?
    TAK — Z_Φ jest **funkcją** struktury akcji, nie liczbą wpisaną.

**Score gate:** ≥ 4/5 PASS (z obowiązkowym U1.2) → Phase 2 forward.

### Phase 2 — UV→IR cascade derivation + ERG cross-check (6)

- **U2.1 P/V derivation z sympy** (β/7·g⁷ − γ/8·g⁸ → potencjał + projekcja):
  jawnie pokazać że `Z_Φ = 14/3` wynika **algebraicznie** z eksponentów,
  nie z dopasowania liczbowego. Independent inputs: `(7, 8, 3, 4)` —
  exponents w P i V. Drift: 0.0 % (exact).

- **U2.2 κ-niezmiennik pod Z_Φ** (sek00:387):
  `κ = 3/(4Φ_eff) = 7/(2Φ₀^bare)` musi być spełnione **algebraicznie**:
  `3/(4·Φ_eff) = 3/(4·Φ₀^bare/Z_Φ) = 3·Z_Φ/(4·Φ₀^bare) = 3·(14/3)/(4·Φ₀^bare) = 14/(4·Φ₀^bare) = 7/(2·Φ₀^bare)` ✓
  Drift: 0.0 % (exact). **Pierwszy realny test struktury** Z_Φ — pokazuje
  że wartość 14/3 jest **wymuszona** przez parametryzację κ na obu poziomach.

- **U2.3 a_Γ-niezmiennik (sek00:388, samospójność DESI DR2)**:
  `a_Γ ≈ 1/Φ_eff` z hipotezy `a_Γ·Φ_0 = 1`. Pod Z_Φ:
  `a_Γ·Φ₀^bare = a_Γ·Z_Φ·Φ_eff = (14/3)·(a_Γ·Φ_eff) = (14/3)·1 = 4.667`.
  Czyli "1" w `a_ΓΦ_0 = 1` (DESI DR2 1.005 ± 0.005) jest **w jednostkach Φ_eff**.
  Test: drift `a_Γ·Φ_eff` vs `1` (DESI DR2). **Drift target:** ≤ 1 %.

- **U2.4 ERG kontrast: Z_Φ vs K_IR/K_UV** (dodatekN twierdzenie 4):
  - `K_IR/K_UV = 1.13` (LPA' Wilsona–Fishera, 13 % renormalizacja kinetyczna)
  - `Z_Φ = 14/3 ≈ 4.667` (substrate dielectric screening, P(1)/V(1) inverse)
  - **Niezależne renormalizacje na różnych skalach**: K_IR/K_UV to
    wave-function renormalization w fazie krytycznej (T ≈ T_c); Z_Φ to
    UV→IR projekcja w fazie low-T (`T ≪ T_c`, sek08c). Brak konfliktu.
  - Drift gate: brak (różne wielkości; zapis interpretacyjny).

- **U2.5 γ.1 multi-anchor compatibility** (open problem γ.1 trade-off):
  Pod Z_Φ = 14/3:
  - `Φ_eff = 8π → Φ₀^bare = 14·8π/3 = 112π/3 ≈ 117.29`
  - `Φ_eff = (10/3)e² → Φ₀^bare = 14·(10/3)e²/3 = 140e²/9 ≈ 114.9`
  - `Φ_eff = 24.783 → Φ₀^bare = 14·24.783/3 = 115.65`
  Wszystkie w paśmie `Φ₀^bare ∈ [114.9, 117.3]` = ±1.0 % wokół `≈ 115`.
  **Z_Φ jest niezmiennikiem γ.1 trade-off** — γ.1 multi-anchor NIE rozprasza Z_Φ.

- **U2.6 Anti-tautology check**:
  Sprawdzić czy `Φ₀^bare` przewidziane z **kosmologii** (`168·Ω_Λ_Planck`)
  zgadza się z `Z_Φ·Φ_eff` używając **gauge-coupling** Φ_eff (Brannen 24.783).
  - Cosmo predicted: `168·0.6847 = 115.03` (z PDG/Planck Ω_Λ)
  - Gauge predicted: `(14/3)·24.783 = 115.65` (z PDG α_s)
  - Drift: `0.54 %` 
  - **Critical:** 0.54 % drift między dwoma niezależnymi anchorami
    (Ω_Λ Planck vs α_s PDG) **odzwierciedla** γ.1 Ω_Λ↔α_s trade-off.
    Nie jest cyrkularny — to jest realna **predykcja**: pod jednym Z_Φ,
    obie skale **muszą** się spotkać w paśmie γ.1 ±1 %.

**Score gate:** ≥ 5/6 PASS → Phase 3 forward.

### Phase 3 — predictions + cross-cycle integration (5)

- **U3.1 UV.2 K_struct ≈ 173 reinterpretacja**:
  pokazać że `K_struct = N_A·2π² = 173.15` z UV.2 to **propagatywny artefakt**
  cyrkularności χ.1 + UV.2, NIE prawdziwa renormalizacja UV. Realna
  renormalizacja UV w TGP to **Z_Φ na poziomie pola Φ**, nie K_struct na
  poziomie skali masowej M_TGP. Pokazać że pod Z_Φ pole `Φ` ma odpowiednio
  zinterpretowaną UV-completion **bez** żadnego dim-less factor `~173`.

- **U3.2 Falsyfikowalna predykcja**:
  Jeśli przyszły pomiar `Ω_Λ` (CMB-S4 2030+) zmieni się o > 1 %, to
  `Φ₀^bare = 168·Ω_Λ` musi się przesunąć synchronicznie z `Φ_eff·14/3`,
  inaczej Z_Φ = 14/3 jest **falsyfikowane**. Test: korelacja Δ Ω_Λ vs Δ α_s
  (PDG) pod hipotezą Z_Φ niezmiennego.

- **U3.3 Status promotion vs UV.2 BLOCKING verdict**:
  - UV.2 `M_TGP DERIVED FULL` rollback do `NUMEROLOGICALLY ANCHORED` (krytyka 2026-05-02)
  - UV.3 promotion: `Z_Φ = 14/3 STRUCTURAL DERIVED` (Phase 1 algebraic + Phase 2 cross-checks)
  - `Φ₀^bare = 168·Ω_Λ` zachowany jako **single dimensionful anchor**
    (Warstwa II Planck calibration), Φ_eff = `(3/14)·Φ₀^bare` derived

- **U3.4 Cross-cycle integration z γ.1**:
  γ.1 H5 dał `Φ_eff = 8π` jako pure structural; pod Z_Φ:
  `Φ₀^bare^pure = 14·8π/3 = 112π/3 ≈ 117.29`. Falsyfikuje się jeśli
  przyszłe Planck/CMB-S4 da `Φ₀^bare > 117.5` (out of band).

- **U3.5 4-channel UV.3 convergence**:
  | # | Channel | Form | Verdict |
  |---|---|---|---|
  | 1 | algebraic P/V | 12/56 = 3/14 (sek00 eq. 74) | exact |
  | 2 | κ parametrization | 3/(4Φ_eff) = 7/(2Φ_0^bare) (sek00:387) | exact |
  | 3 | Cosmological calibration | 168·Ω_Λ_Planck = 115.03 | ≤ 0.1 % vs (14/3)·Φ_eff |
  | 4 | Gauge-coupling cross-check | (14/3)·24.783 = 115.65 vs cosmo 115.03 | ≤ 0.6 % (γ.1 trade-off) |

**Score gate:** ≥ 4/5 PASS → UV.3 program END (renormalization explicit).

## Promotions post-UV.3 (jeśli pełna konwergencja)

- **Z_Φ = 14/3 NEW ANCHOR** w status_map (wave-function renormalization, EXACT)
- **Φ₀^bare ↔ Φ_eff** explicit dwuwarstwowy formalism (UV: ≈ 115, IR: ≈ 24.66)
- **Sek00 tabela 385–388** rozwinięta o jawną definicję `Z_Φ` z derivacją `P/V`
- **UV.2 K_struct = N_A·2π²** rollback do "numerologiczna obserwacja na poziomie skali masowej, propaguje cyrkularność χ.1"
- **γ.1 multi-anchor** podkonsumowany przez Z_Φ niezmiennik (Φ₀^bare niezmienne ±1 % pod γ.1 trade-off)

## Falsification path (UV.3 zostanie odrzucone jeśli)

1. Sympy P(1)/V(1) z sek00 eq. 74 nie daje `12/56 = 3/14` (definicyjna sprzeczność)
2. `(14/3)·Φ_eff` (gdziekolwiek wziąć Φ_eff w paśmie γ.1) odbiega od `168·Ω_Λ_Planck` o > 1 %
3. κ parametrization `3/(4Φ_eff) ≠ 7/(2Φ₀^bare)` algebraicznie (test U2.2 fails)
4. Przyszłe Planck/CMB-S4 da `Ω_Λ` poza zakresem `Φ₀^bare ∈ [114, 118]`/(168) = `[0.679, 0.702]`
5. ERG flow (dodatekN) wykryje **dodatkowy** czynnik renormalizacji UV→IR pomijający Z_Φ = 14/3

## Open frontiers post-UV.3 (jeśli zamknięte)

- **Mechanizm dynamiczny Z_Φ** (czemu eksponenty `(7,8,3,4)` w P i V, a nie inne)
- **Powiązanie z η_N* = -2** UV.1 NGFP (czy Z_Φ jest IR-side renormalizacja
  podczas gdy η_N* = -2 jest UV-side? Możliwa **single-scale** unifikacja)
- **Re-derywacja `α_em`, `α_s`, `κ`** explicit w terms of `Φ₀^bare` (nie Φ_eff)
  — sprawdzić czy formuły w sek09, dodatekO, dodatekV są równoważne pod
  podstawieniem `Φ_eff = (3/14)·Φ₀^bare`

## Cross-references

- [[../op-uv2-mtgp-absolute-scale/CRITIQUE_repackaged_circularity_2026-05-02.md]] — krytyka UV.2 fittingu
- [[../op-uv-as-ngfp/Phase3_results.md]] — UV.1 NGFP {g* = 0.71, λ* = 0.19, η_N* = -2}
- [[../op-gamma1-phi-eff-anchor-resolution/README.md]] — γ.1 H5 Φ_eff = 8π identification
- [[../../core/sek00_summary/sek00_summary.tex]] linie 60–124, 380–390 — kanoniczna tabela Φ₀ / Φ_eff
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] linia 446 — sek08 "Φ_eff = Φ_0(1+δ) z 1-loop" (terminology)
- [[../../core/formalizm/dodatekN_erg_renormalizacja.tex]] — ERG flow, K_IR/K_UV = 1.13
- [[../../core/formalizm/dodatekQ_coarse_graining_formal.tex]] Q.4 — `a_Γ·Φ₀ = 1` samospójność
- [[../../core/formalizm/dodatekV_su3_formalizacja.tex]] linie 130, 141, 151 — Φ_eff w α_s
- [[../../core/formalizm/dodatekO_u1_formalizacja.tex]] linia 437 — Φ_0 = 24.66 w α_em
