---
title: "Honest framing — cztery korzenie UV/CG TGP (α=2, c₀, 𝒜→α_s, K_geo): synteza statusu do submission"
date: 2026-06-26
type: meta-synthesis
status: 🟢 ACTIVE — 1-stronicowa nota syntetyczna (Limitations / honest framing)
purpose: "Pojedyncze źródło prawdy o statusie 4 stałych derywowalnych tylko przez UV/CG; gotowy blok do sekcji Limitations submission."
synthesizes: ["#37 op-c0-derivation-from-substrate", "#38/#39 op-Phi-field-identity / op-bond-order-RG-selection", "#43 op-A-derivation-from-CG", "#48 op-Kgeo-from-D-uniqueness", "#49 op-CG-alpha-eff-convergence"]
parent: "[[../STATE.md]]"
tags: [meta, honest-framing, uv-cg, limitations, anti-lakatos]
---

# Honest framing — cztery korzenie UV/CG TGP

## §0 — Teza (jedno zdanie)

Makro-fenomenologia TGP (sektor leptonowy 1→3 do 0,006%, PPN exact, OP-7 tensor 96,9% PASS,
kosmologia) jest domknięta, **ale cztery kluczowe stałe nie są derywowane z substratu** — są
**selekcjami aksjomatycznymi na gęstości**, których jedyna droga do derywacji prowadzi przez ten
sam, niedomknięty most coarse-grainingu Γ→Φ (NGFP / UV). To uczciwie obniża pozorną parsymonię,
**ale nie zmienia bilansu parametrów** (#42) ani nie falsyfikuje teorii.

## §1 — Mapa czterech korzeni (KOMPLETNA, 4/4)

| Korzeń | Czym jest | Status | Dlaczego warunkowy | Cykl(e) |
|---|---|---|---|---|
| **α = 2** | wykładnik kinetyczny `K(Φ)=K_geo·Φ⁴` (klasa konforemna C1–C3) | **CLOSED-NEGATIVE** (aksjomat na gęstości) | substrat NIE generuje α=2: kanoniczny `s=0` → `α_eff=−½` (chain-rule exact + FSS value-blind, werdykt (B)); α=2 wymaga `s=5` (`K∝ŝ¹⁰`) — niekanoniczne (#38), RG-irrelevant (#39), escape `η~O(5)` zamknięty (#49); **mapa obstrukcji KOMPLETNA, NON-DERIVABLE z mikro `H_Γ`** (#53) | #38/#39 + **#49 + #53** |
| **c₀** | normalizacja kinetyczna tensora `C_σ≡κ_E` (= sprzężenie metryczne `c₀=C(ψ=1)`) | **FREE-PARAMETER** (UV) | liniowa rozbieżność UV (−16/35), brak izolowanego bieguna spin-2 ⇒ brak scheme-independent continuum; „c₀≈4π" to *matching* + kalibracja GW150914, nie derywacja | #33/#34 + #37 |
| **𝒜 → α_s** | stała konfinementu `𝒜=a_Γ/φ=C_F²α_s²` | **POSTULATE-CONDITIONAL** | cały łańcuch wisi na 1 postulacie `K_geo·m_sp²=π·Φ₀²`, redukującym się do niedomkniętego CG-1/CG-3 (ex200 4/8); α_s 0,03σ = consistency-check warunkowy, NIE first-principles | #43 |
| **K_geo** | bezwymiarowy prefaktor geometryczny `K(φ)=K_geo·φ⁴` | **POSTULATE-CONFIRMED** | `thm:D-uniqueness` fiksuje FORMĘ (φ⁴) + α=2, ale `K_geo` to wolna stała całkowania (C3 ją definiuje, nie liczy); absorbowalna; π geometryczne, ale skala `K_geo·m_sp²` nie | #48 |

## §2 — Wspólny mianownik

Wszystkie cztery zbiegają do jednego: **pełnego domknięcia mostu Γ→Φ (continuum/NGFP)**. Stan tego mostu:
- **NGFP analitycznie** (`op-uv-as-ngfp`, AS Reuter/Litim): 7/7 PASS — `g*·λ*`, `η_N*=−2`, `N_A=500/57`.
- **Most substrat→makro** (CG-1/CG-3, `status_map`): **[SZKIC]/[OTWARTY]** — `ex200` α_eff niezbieżny (4/8),
  `ex202` σ_TGP nie zamknięte (7/8, T6 FAIL), CG-1 (Banach) [OTWARTY].
- Domknięty: CG-2 (LPA' 8/8), CG-3 (homogenizacja H¹, 5/5, #31), CG-5 (Φ₀ 8/8).

⟹ **Luka jest jedna, dobrze zlokalizowana, wieloletnia** (analityka funkcjonalna + ERG + homogenizacja
de Giorgi), o niskim priorytecie inżynieryjnym a wysokim fundamentalnym.

### §2.1 — Korzeń α=2: mapa obstrukcji KOMPLETNA (#53, 2026-06-27)

Cykl `op-CG-Kij-from-Hgamma` domyka korzeń α=2 **analitycznie** (werdykt **NON-DERIVABLE**,
value-blind). `α=2` (sprzężenie geometryczne `K_ij=J(φ_iφ_j)²`) **nie wyprowadza się** z mikro
`H_Γ` (eq:B-H) żadnym coarse-grainingiem — cztery niezależne ramiona obstrukcji:

1. **Gaussian:** bilinearny bond `−Jŝ_iŝ_j` (jedyny obecny w `H_Γ`) → kinetyka kanoniczna w `ŝ`
   → `α_eff=−½` w gęstości `Φ=⟨ŝ²⟩` (potwierdza #49 jako baseline/anchor).
2. **η-escape:** `η_3D-Ising≈0,036 ≪ 1` → wymiar anomalny nie zmostkuje luki (B-REFUTED).
3. **RG-relevance:** `Δ[Φⁿ(∇Φ)²]=(n+2)Δ_ε+2 > d=3` dla `n=−1,0,1,2` (`Δ_ε≈1,413`) — **cały**
   sektor kinetyczny kompozytu `Φ=ε` RG-irrelewantny ⟹ wykładnik nie pinowany przez FP =
   nieuniwersalny datum UV = efektywnie aksjomat (potwierdza #39).
4. **Stopień bondu:** `α=2` wymagałby szesciopolowego bondu `(ŝ_iŝ_j)³` — nieobecnego w `H_Γ`,
   o wolnym współczynniku, RG-irrelewantnego ⟹ nowy aksjomat, nie derywacja.

**Bonus (R5):** trzy wykładniki `{−½, 1, 2}` to nie jedna rodzina, lecz trzy konstrukcje
(`H_Γ≠F_kin`): bilinear→`−½`; `(φφ)²` jako energia bondu→`α=1`; `(φφ)²` jako współczynnik
sztywności→`α=2`. Tylko bilinear jest mikroskopowo obecny.

> **DOUBT W-CGK-1 (rama).** „Headline #49" `−1→+4` (Δe=5) zestawia `e` substratu w **gęstości**
> z `e` manuskryptu w **amplitudzie** (rama mieszana). W spójnej ramie luka = **Δe=4 (amplituda)**
> lub **Δe=2 (gęstość)**, nie 5. **Werdykt (B) #49 pozostaje robust** (luka ≠ 0 w obu ramach);
> korekta dotyczy wyłącznie magnitudy argumentu η, nie wniosku. #49 LOCKED, niezmienione.

⟹ α=2 jest **nieredukowalnym, uczciwym aksjomatem** (nie ukrytą porażką). To wzmacnia, nie
falsyfikuje, parsymonię ledgera (#42, `N_axiom=6` bez zmian).

## §3 — Dlaczego to NIE jest słabość falsyfikująca

1. **Bilans parametrów bez zmian (#42):** α=2 ∈ N_axiom=6 (już liczone jako aksjomat); c₀ ∈ FREE(8);
   α_s = TRADED(2); K_geo wbudowane. Uczciwy `N_free=10 ≪ 19` (SM) stoi; korona leptonowa (1→3) stoi.
2. **Makro-predykcje niezależne od UV:** struktura GW (2 TT + breathing, `c_GW=c`, `m_σ²=2m_s²`),
   PPN (γ=β=1), masy leptonów — wszystkie `c₀`/`κ_E`-independent.
3. **Dwustronna falsyfikowalność zachowana:** `PR-025 forward (b)` — jeśli domknięcie CG wymusi
   `K_geo·m_sp²≠π·Φ₀²`, most 𝒜=C_F²α_s² jest **sfalsyfikowany**.
4. **Werdykt (B) #49 ≠ falsyfikacja:** „substrat nie derywuje α=2" **ratyfikuje** istniejącą uczciwą
   pozycję (α=2 = selekcja, nie derywacja), domykając zawieszony residual numerycznie.

## §4 — Drop-in prose do submission (EN, Limitations)

> **Limitations — the UV/coarse-graining roots.** TGP's macroscopic phenomenology is self-contained,
> but four quantities are *axiomatic selections on the density field rather than substrate derivations*:
> the kinetic exponent α=2, the tensor normalization c₀ (≡ the free UV constant C_σ), the confinement
> constant 𝒜 (hence α_s), and the geometric prefactor K_geo. Each is derivable only by closing the same
> substrate-to-continuum (coarse-graining/NGFP) bridge, which remains open (CG-1 Banach existence; α_eff
> non-convergence on accessible lattices). We make this explicit rather than absorbing it into the headline:
> a value-blind, pre-registered finite-size-scaling study confirms that the canonical substrate yields
> K(Φ)∝Φ^{−1} (α_eff=−½), not the K(Φ)∝Φ⁴ that α=2 would require — so α=2 is carried as a conformal-class
> selection (C1–C3), not a derivation. This lowers the apparent parsimony but leaves the honest free-parameter
> count unchanged (N_free≈10 vs ~19 for the Standard Model) and does not affect the UV-independent predictions
> (lepton masses, PPN, GW mode content c_GW=c). The bridge closure is a falsifiable program, not a hidden
> assumption: a closed coarse-graining that forced K_geo·m_sp²≠π·Φ₀² would refute the α_s consistency bridge.

## §5 — Anti-Lakatos
✓ Bias dwustronny: korzenie jawne (nie chowane), aksjomaty nie liczone jako parametry. ✓ (B)/FREE/
CONDITIONAL/CONFIRMED — wszystkie statusy negatywne/warunkowe jawne. ✓ Falsyfikowalność (PR-025) zachowana.
✓ Zero nowych claimów — synteza zalockowanych werdyktów #37–#53 (α=2: mapa obstrukcji kompletna, #53).

## Cross-references
- [[../research/op-CG-alpha-eff-convergence-2026-06-26/Phase_FINAL_close.md]] (#49, α=2)
- [[../research/op-CG-Kij-from-Hgamma-2026-06-27/Phase_FINAL_close.md]] (#53, α=2 mapa obstrukcji KOMPLETNA — NON-DERIVABLE)
- [[../research/op-c0-derivation-from-substrate-2026-06-22/]] (#37, c₀)
- [[../research/op-A-derivation-from-CG-2026-06-25/Phase_FINAL_close.md]] (#43, 𝒜)
- [[../research/op-Kgeo-from-D-uniqueness-2026-06-26/Phase_FINAL_close.md]] (#48, K_geo)
- [[../research/op-parameter-counting-balance-sheet-2026-06-25/Phase_FINAL_close.md]] (#42, ledger)
- [[../research/op-uv-as-ngfp/]] (NGFP) · [[../core/_meta_latex/status_map.tex]] (CG-1..CG-5)
- [[PRE_REGISTERED_FALSIFIERS.md]] PR-025 (forward falsifier (b))
