---
title: "op-bh-alpha-threshold — Structural status of α(ψ) in Candidate D BH action"
date: 2026-04-28
cycle: BH.1
status: ACTIVE
predecessor: "[[../op-m92/OP_M92_P0plus_candD_multisource_results.md]] (multi-source ISSUE: α_SI spans 19 orders of magnitude)"
flask: "tgp-core (DOI 10.5281/zenodo.19670324)"
related:
  - "[[../../INDEX.md]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
  - "[[../op-m92/OP_M92_readiness_summary.md]]"
  - "[[../op-m92/OP_M92_Phase1_PLAN.md]]"
tags:
  - TGP
  - BH
  - alpha-psi
  - threshold
  - multi-source
  - candidate-D
---

# op-bh-alpha-threshold — Structural status of α(ψ) in Candidate D BH action

> **Cel:** Zamknąć **Phase 0+ multi-source ISSUE** — pokazać, czy "stała"
> α z Candidate D action `S_int ∝ α · T^μν J_μ J_ν` jest **pojedynczą
> uniwersalną stałą fizyczną** (single number in SI), czy też z konieczności
> **funkcją scalar-field α(ψ)** z progiem ψ_th ≈ 1.05. Cykl strukturalny
> mini-pattern (3 fazy) — analog SC.1, NIE duplikuje pełnego Phase 1 PLAN
> (15 miesięcy paper-rigor M9.2-D derivation).

---

## Tło

W OP-M92 Phase 0+ Candidate D structural sketch (2026-04-25) heurystycznie
wyznaczono:

```
α_geom ~ 0.1     w jednostkach geometrycznych (M_BH = 1)
                  przy photon ring r = 3.88 M, ψ = 1.168
                  dla scenariusza (e) target shift +14.56%
```

Sub-test multi-source (Phase 0+ tego samego dnia) pokazał **fundamentalny
problem strukturalny**:

| Source | M (M_⊙) | R_S (m) | α_SI (m²) |
|--------|---------|---------|-----------|
| Sgr A* | 4.3·10⁶ | 1.27·10¹⁰ | **1.61·10¹⁹** |
| M87* | 6.5·10⁹ | 1.92·10¹³ | **3.69·10²⁵** |
| GW150914 final | 65 | 1.92·10⁵ | **3.69·10⁹** |
| Neutron star | 1.4 | 4.13·10³ | **1.71·10⁶** |

α_SI skaluje się jak **M_BH²** (przez R_S² conversion z geom units),
co daje **19 rzędów wielkości** przy 9.5 rzędach w masie BH.

⇒ α **nie jest** pojedynczą uniwersalną stałą fizyczną pod naiwnym
Candidate D action z α dimensionful.

Sketch resolution paths z OP_M92_P0plus_candD_multisource_results.md:
- **Path A** (L_TGP intrinsic length scale): FAILS (1.5·10³ rozjazd)
- **Path B** (substrate density enhancement ρ_sub ~ M/R_S³): FAILS
- **Path C** (α/M_Pl² prefactor): FAILS unphysically (α ~ 10⁸⁹)
- **Path D** (scenario (e) target nie-uniwersalny): nie-valid (geometria photon ring uniwersalna)
- **Path E** (α(ψ) z progiem): **VIABLE** — preferowana
- **Path F** (composite Candidate D + B): możliwa hybrid

**Path E hipoteza:**
```
α(ψ) = α_0 · (ψ - ψ_th)ⁿ · Θ(ψ - ψ_th)
ψ_th ≈ 1.05    (lab ψ ≈ 1, photon ring ψ = 1.168)
α_0 dimensionless O(1)
```

Lab regime (ψ ≈ 1) — α się **nie aktywuje** → SI suppression naturalna.
Strong-field (BH photon ring ψ = 1.168, NS surface ψ ~ 1.4) — α aktywne.
ψ-zależność absorbuje M² scaling bo `ψ(r=3.88M)` jest **uniwersalne w
geom units** ale `ψ(lab Earth)` ≈ 1 niezależnie od pobliskich BH.

---

## 3-fazowy plan

### Phase 1 — Multi-source dimensional + mass-scaling audit

**Cel:** Sformalizować i zamknąć Phase 0+ ISSUE — pokazać że α z
Candidate D **nie może** być uniwersalną stałą pod **żadną** prostą
unit-bridge transformation (analog SC.1.Phase1: H₀ unit-cousin
α_PB ↔ α_0 → REJECTED).

**Hipoteza nullowa H₀:** α z Candidate D **da się** zinterpretować
jako pojedyncza fizyczna stała (po zastosowaniu jakiegoś czynnika
konwersji f między geom-units a SI).

**Hipoteza alternatywna H₁:** α **musi** być funkcją scalar-field α(ψ)
(albo równoważną nielokalną modyfikacją), aby M² scaling zostało
zaabsorbowane.

**Testy:**
- T1.1 — formalna analiza wymiarowa Candidate D action `S_int = α ∫ T^μν J_μ J_ν √-g d⁴x` (sympy `dimsys_SI`); identyfikacja dim[α].
- T1.2 — explicit M² scaling demonstration: α_SI = α_geom · R_S² dla 4 źródeł, sympy/numerical (replikacja istniejącego wyniku w czystej formie + tabelaryzacja).
- T1.3 — closure paths A/B/C explicit zawodzą (pokaż dla każdej path single dimensionless `f` nie istnieje).
- T1.4 — lab-suppression check: `ψ_lab ≈ 1` — pokaz że pod hipotezą Path E `α(ψ_lab) ≈ 0` exponentially/polynomially → naturalne tłumienie eksperymentów lab-scale.
- T1.5 — Path E vs Path D/F summary table: które rozwiązania pozostają strukturalnie viable po Phase 1.

**Wynik oczekiwany:** H₀ rejected w 1 cyklu (analog SC.1.Ph1).
Stwierdzenie: α(ψ) (lub strukturalnie ekwiwalentny mechanizm) jest
wymaganą modyfikacją Candidate D dla multi-source uniwersalności.

**Status:** zamknięcie negatywne H₀ — cel jest pokazać że Path E (nie konstanta)
to JEDYNA strukturalnie żywotna kontynuacja.

### Phase 2 — Path E α(ψ) substrate-physics derivation

**Cel:** Wyprowadzić formę funkcji `α(ψ) = α_0 (ψ - ψ_th)ⁿ Θ(ψ - ψ_th)`
z TGP core + EFT classification (analog SC.1.Phase2 Abrikosov–Gorkov
first-principles). Cel: określić **n**, **ψ_th**, oraz związek z
`α_0 ≈ 4.04` z closure_2026-04-26 T-α threshold.

**Strategia:**
- Sympy-symbolic EFT analysis: z lokalnej kowariantnej akcji ∫ f(ψ) T^μν J_μ J_ν √-g d⁴x na 4D Lorentzian, jakie f(ψ) są naturalne?
- Power-law z progiem (ψ - ψ_th)ⁿ vs analytical (1 - exp(-κ(ψ - ψ_th))) vs Heaviside threshold — klasyfikacja form, które redukują się do α_0 ≈ const w strong-field i znikają w lab.
- Konwergencja z α_0 z T-α threshold — sprawdzenie czy stała `α_0` w Path E to **ta sama** α z T-α (ψ-substrate threshold) czy strukturalnie różna.
- Predykcja `n` z TGP scaling: relacja do κ_TGP ≈ 2.012 (jeśli n = κ_TGP, czy działa multi-source?).
- ψ_th wyznaczone z multi-source consistency: dla danych 4 źródeł (Sgr A* photon ring ψ=1.168, M87* ψ=1.168 [photon ring uniwersalne!], GW150914 final ψ ~ 1.2, NS surface ψ ~ 1.4) jaki ψ_th + n minimalizuje multi-source RMS_log α_SI/α_universal?

**Falsyfikacja:** jeśli **żadne** (n, ψ_th, α_0) nie daje zgodności
multi-source ≤ 30% spread, Path E jest strukturalnie niewystarczający
i wymagana jest Path F (composite D+B) lub powrót do Phase 1 PLAN.

### Phase 3 — Multi-source falsification map

**Cel:** Zarejestrować **uniwersalne** Path E predictions α(ψ) jako
multi-source falsification map (analog SC.1.Phase3 multi-LnH9):

- ngEHT photon ring odlewki: Sgr A* (2030+), M87* (2030+), kandydaci pośredni (NGC 5128 Cen A, M84 — masa BH ~10⁸-10⁹ M_⊙).
- LIGO/LISA strong-field merger ringdown: GW150914-class (M~65 M_⊙) vs supermassive merger LISA (M~10⁶ M_⊙).
- NS surface measurement (NICER pulsars): czy α(ψ) aktywuje się przy ψ ~ 1.4?

**Predykcja kluczowa:** universal +14.56% deviation w b_crit photon ring
dla **wszystkich** Schwarzschild sources, jeśli Path E z (n, ψ_th, α_0)
z Phase 2 — ZGODNE z fenomenologią bez M² scaling problem.

**Sukces:** publikacja BH4-BH7 entries w PREDICTIONS_REGISTRY z konkretnymi
horyzontami eksperymentalnymi i clean discriminator factors.

**Falsyfikacja:** jeśli dwie ngEHT measurements (Sgr A* i M87*) dają
**różne** deviation factors poza 1σ ngEHT precision → Path E falsified,
multi-source Issue otwarte → eskalacja do pełnego Phase 1 PLAN
(15-month covariant derivation).

---

## Co robimy najpierw

**Phase 1 — TERAZ.** Mała deterministyczna analiza wymiarowa + mass-scaling
demonstration. 1-cyklowe zamknięcie negatywne H₀ (analog SC.1.Ph1).

Phase 2 i 3 queued — wymagają user-go-ahead lub zamknięcia Phase 1.

---

## Cross-references

- `research/op-m92/OP_M92_P0plus_candD_multisource_results.md` — multi-source ISSUE źródłowy, paths A-F sketch
- `research/op-m92/OP_M92_Phase1_PLAN.md` — pełny 9-15-month PLAN (świadomie NIE duplikujemy)
- `research/op-m92/OP_M92_readiness_summary.md` — Phase 0/0+ readiness gate
- `research/closure_2026-04-26/alpha_psi_threshold/results.md` — α_0 ≈ 4.04 T-α origin (relacja do α_0 w Path E?)
- `research/op-sc-alpha-origin/program.md` — wzorzec strukturalnego mini-cyklu 3-fazowego
- `INDEX.md` — flask deposits map
- `PREDICTIONS_REGISTRY.md` — BH1-BH3 sektor (photon ring), kandydaci na BH4-BH7

---

## Decyzja po Phase 1

- Jeśli Phase 1 = H₀ rejected (oczekiwane): negatywny audyt zamyka cykl BH.1.Ph1, Path E zostaje JEDYNĄ żywotną kontynuacją, propose Phase 2 do user-a.
- Jeśli Phase 1 = H₀ supported (mało prawdopodobne — wymagałoby ukrytej unit-bridge konwersji niewykrytej w Phase 0+): rozszerz audyt o explicit f-derivation, raport do user-a.

---

## Relacja do Phase 1 PLAN (op-m92/OP_M92_Phase1_PLAN.md)

**Ten cykl NIE zastępuje** pełnego Phase 1 PLAN (paper-rigor 9-15-month
covariant derivation z 4D Lagrangian symetrii, perturbative expansion,
loop checks). Jest **strukturalnym mini-audytem** o zakresie:

- formalizuje multi-source ISSUE (Phase 0+ był tylko sketch)
- klasyfikuje Path E jako jedyną strukturalnie żywotną drogę
- wyprowadza formę α(ψ) na poziomie EFT (nie pełna covariant action)
- rejestruje multi-source predictions dla ngEHT/LISA falsification

Po zamknięciu BH.1.Phase3 user może zdecydować:
- (a) zamknij sektor BH na poziomie EFT — pełny Phase 1 PLAN nie potrzebny
- (b) eskaluj do pełnego Phase 1 PLAN (15-month) jeśli BH.1 odsłoni nowe wymagania paper-rigor
- (c) hybrydowo: BH.1 + selected paper-rigor sub-tests z Phase 1 PLAN

Decyzja queued na po-Phase-3.
