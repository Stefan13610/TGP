---
title: "Decyzja konwencji PN: energy-PN vs phase-PN dla T01 / M911-P1"
date: 2026-05-07
parent: "[[README.md]]"
type: convention-decision
tgp_owner: audyt/T01_LIGO3G_falsifier
tags:
  - convention
  - decision
  - PN-counting
  - phase-PN
  - energy-PN
  - ppE
  - 3PN
  - 2PN
  - M911
  - T01
  - EXT-5
related:
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[PPN_TO_PPE_MAPPING.md]]"
  - "[[FALSIFIER_STATEMENT_DRAFT.md]]"
  - "[[../../research/op-newton-momentum/M9_1_pp_P1_results.md]]"
---

# CONVENTION_DECISION — PN counting dla T01

> **Cel pliku.** W [[PPN_TO_PPE_MAPPING.md]] §1 (KRYTYCZNA UWAGA)
> zidentyfikowano niejednoznaczność: **U³ deviation w `g_tt`** może
> być nazwana "3PN" lub "2PN" w zależności od konwencji. Ten plik
> rozstrzyga, **którą konwencję rekomenduję dla T01**, dlaczego
> i jakie są implikacje. Decyzja jest **rekomendacją audytu** —
> autorska *adopcja* tej konwencji w `PREDICTIONS_REGISTRY.md` jest
> ostatecznym krokiem.

## §1 — Problem

### 1.1 Dwie konwencje, ten sam fakt fizyczny

W M9_1_pp_P1_results.md §3.2 wyprowadzono explicit:

```
Δg_tt / (-c²) = − (5/6) U³ + (23/12) U⁴ − (19/6) U⁵ + (337/72) U⁶ + …
```

z `U = GM/(rc²)` (Newtonian potential).

Tę samą deviation (U³ w `g_tt`) literatura nazywa **dwoma różnymi
sposobami**:

#### Konwencja A — "Energy PN" (intuicyjna, używana w README/M9_1_pp_P1)

Każdą power U ponad leading liczy jako "1 PN order":
- 0PN = Newton (U^1 w g_tt)
- 1PN = O(U²) deviation w g_tt
- 2PN = O(U³) deviation w g_tt
- **(5/6) U³ → 3PN deviation w g_tt** ← README T01

To jest naturalne dla *statycznego* PPN expansion (Will 1971;
Will–Nordtvedt 1972). Mercury, Cassini, LLR są testowane na 1PN
β/γ_PPN, czyli O(U²) effects w `g_tt`. M9_1_pp_P1 używa tej konwencji.

#### Konwencja B — "Phase PN" (standardowa w GW literature)

Każdy "PN order" w fazie waveformu inspiralu odpowiada *connected*
power w SPA:
- 0PN phase = quadrupole leading: `Ψ ~ u⁻⁵`, gdzie u = (πMf)^(1/3)
- 1PN phase = `Ψ ~ u⁻³` (Wagoner–Will 1976)
- 2PN phase = `Ψ ~ u⁻¹` (Blanchet–Damour–Iyer 1995)
- 3PN phase = `Ψ ~ u^+1` (Blanchet–Faye–Iyer 2002, complete)
- **U³ deviation w `g_tt` → 2PN deviation w fazie waveformu**

Cutler–Flanagan 1994, Yunes–Pretorius 2009, Blanchet *Living Reviews*
2014 — wszystko używa konwencji B.

### 1.2 Mapowanie między konwencjami

Dla **statycznej, sferycznie symetrycznej** metryki:

```
power U^N w g_tt  ↔  energy-PN: N − 1     ↔  phase-PN: N − 2     ↔  ppE: b = 2N − 7
```

Sprawdzenie:

| Power U^N | Energy-PN | Phase-PN | b_ppE |
|-----------|-----------|----------|-------|
| U¹ (Newton) | 0PN     | (leading) | −5    |
| U² (1PN)   | 1PN      | 1PN      | −3    |
| U³ (2PN)   | 2PN      | **2PN**  | **−1** |
| U⁴ (3PN)   | 3PN      | 3PN      | +1    |
| U⁵ (4PN)   | 4PN      | 4PN      | +3    |

**KRYTYCZNA UWAGA do tabeli.** Mapping `Phase-PN: N − 2` zakłada
że *każda* power U w g_tt propaguje się 1:1 do fazy waveformu w
SPA. Dla **modyfikacji TYLKO g_tt** (jak w M9.1''), to jest *blisko*
prawdy — ale wymaga walidacji w cyklu `op-ppE-mapping/` (zob.
[[PPN_TO_PPE_MAPPING.md]] §2.1, założenia A1–A6).

### 1.3 Konsekwencje dla niejednoznaczności

W rozmowie konwersacyjnej powiedzenie "TGP ma 3PN deviation":
- *interpretowane przez fizyka grawitacyjnego* → "phase 3PN ≡
  Blanchet–Faye–Iyer kompletne 3PN coefficient → b_ppE = +1"
- *interpretowane przez audytora TGP* → "energy 3PN ≡ U³ w g_tt
  → b_ppE = -1"

**To są ZUPEŁNIE INNE rzeczy** w terms of:
- which detector channel constrains it (3PN vs 2PN phase Fisher bounds różnią się o ~10×),
- which other modified gravity theories share the b_ppE (b=−1: dCS, sGB, EÆ; b=+1: Brans-Dicke higher-PN mixing, MG quadratic).

## §2 — Rekomendacja: konwencja PHASE (B)

### 2.1 Argumenty za PHASE convention

1. **Standardowa w peer-review.** Wszystkie ppE papers (Yunes, Yagi,
   Chamberlain, Sampson, Will–Yunes), wszystkie LIGO Tests-of-GR
   papers (LSC GWTC-1/2/3 ToGR), wszystkie ET science cases
   (Maggiore 2020, ET Steering Cmt 2023) — wszystko jest w PHASE
   convention.

2. **Bezpośredni klucz do rejestru.** PREDICTIONS_REGISTRY ma już
   wpisy GW1–GW6, BH5, ε.1 z explicit "f-frequency dependent
   phase deviation" sformułowanymi w PHASE convention. M911-P1
   wpis musi być spójny.

3. **ppE Fisher bounds w literaturze.** Tabele Chamberlain–Yunes
   2017, Yagi–Yunes 2016 są w PHASE; jeśli T01 używa ENERGY
   convention, każdy peer-review reader musi przekonwertować →
   confusion + reduced credibility.

4. **Bezpośrednia eksperymentalna obserwabla.** GW detector
   bezpośrednio mierzy *fazę* sygnału, nie *metryczne power U*.
   Phase convention jest closer to data.

### 2.2 Argumenty przeciw

1. **M9_1_pp_P1_results.md jest w ENERGY.** Cykl historycznie
   używa energy convention. Switch wymaga przepisania komentarzy
   w README.md i naprawienia narrative w `TGP_FOUNDATIONS.md` §3.

2. **PPN tests (Mercury/Cassini) są naturalnie w ENERGY.** β_PPN,
   γ_PPN są 1PN w energy convention; 1PN w phase convention to
   *coś innego* (1PN phase). Mixing w jednym papierze jest
   niewygodne.

### 2.3 Propozycja kompromisu — DUAL nomenclature

Adoptować PHASE jako **primary** w wpisach do
`PREDICTIONS_REGISTRY.md` i peer-review documents, ALE zachować
ENERGY w wewnętrznych M9.x analitycznych dokumentach (gdzie
β_PPN, γ_PPN są 1PN zgodnie z PPN literature).

**Implementacja w rejestrze:**

```markdown
### M911-P1 — 3PN-energy / 2PN-phase deviation in inspiral phase

**Convention:** Below "3PN" refers to *energy/PPN* counting (Will 1971,
master formula sympy P23): U³ in g_tt expansion. In *phase* counting
(Cutler–Flanagan 1994; standard in waveform literature), this is the
**2PN phase correction**, ppE coefficient **b_ppE = −1**.

**Mapping table (resolves convention ambiguity):**

| Power in g_tt | Energy-PN (PPN) | Phase-PN (waveform) | ppE b |
|---------------|-----------------|----------------------|-------|
| U³            | 3PN-energy      | 2PN-phase           | −1    |

**Predykcja:** δβ_ppE^(2PN-phase) = κ · (5/6), gdzie κ is SPA
chain prefactor (~0.4–1.5, locked in `op-ppE-mapping/` Phase 1).
…
```

## §3 — Implikacje dla istniejących artifacts T01

### 3.1 [[PPN_TO_PPE_MAPPING.md]] §3 (dictionary table)

Aktualnie tabela w §3 ma kolumnę "Konwencja 'metric PN'" i "Konwencja
'phase PN'" osobno, z dual b_ppE wartościami. **TRZYMAJ AS IS** —
to jest dictionary, nie decyzja.

W §1 zdanie "**przyjmuję konwencję PHASE**" jest **już w zgodzie
z tym dokumentem**. Bez zmian.

### 3.2 [[FALSIFIER_STATEMENT_DRAFT.md]] §1 tabela

Tabela detekcji ma "ppE 3PN coefficient" — to jest *PHASE 3PN*,
czyli b_ppE = +1, nie b_ppE = −1!

**KOREKTA WYMAGANA.** (5/6) U³ w g_tt → **2PN-phase** → b = −1.
Wpis falsifier powinien mówić **"2PN-phase coefficient"**, nie
"3PN coefficient". Aktualnie tekst w FALSIFIER_STATEMENT_DRAFT.md §1
zawiera niejednoznaczność.

Korekta: zob. §4 niniejszego pliku — patch instructions do FALSIFIER_STATEMENT_DRAFT.md.

### 3.3 [[SENSITIVITY_BACK_OF_ENVELOPE.md]] §3 tabela bounds

Tabela §3.1 poprawnie podaje bounds dla **2PN (b=−1)** ~10⁻¹/10⁻³/3·10⁻⁴
i **3PN (b=+1)** ~3·10⁻¹/3·10⁻³/10⁻³ jako **dwie linie**. Decyzja
PHASE-convention oznacza, że **tylko linia 2PN (b=−1) jest
relewantna** dla M911-P1.

Tabela §4.2 (TGP vs bounds) musi używać linii 2PN.

**KOREKTA WYMAGANA.** Aktualne numbers w §4.2 ("β_ppE_5σ_bound LIGO-O3
~3·10⁻¹") są w istocie 3PN line; **dla 2PN (właściwa M911-P1)
powinny być ~10⁻¹** w O3 (lepszy bound, jaśniej "borderline"
pre-falsifier).

### 3.4 [[README.md]] sekcja "Diagnoza (z EXT-5)" i "Problemy jakie rodzi" #4

Aktualnie problem #4 mówi "U³ to człon 3PN" — to **konwencja
ENERGY**. Po decyzji PHASE convention, należy przeformułować:

> "U³ to member 2PN-phase (równoważnie 3PN w PPN/energy convention).
> LIGO obecnie ma 3PN-phase complete (Blanchet-Faye-Iyer 2002),
> 3.5PN-phase aktualnie używane (Mishra 2016). LIGO 3G celuje w
> 4PN-phase. Dla 2PN-phase coefficient (M911-P1 b_ppE=-1), bounds
> z LIGO-O3 są ~10⁻¹, ET-D ~10⁻³."

## §4 — Patch instructions

### 4.1 Patch do [[FALSIFIER_STATEMENT_DRAFT.md]]

W §1 tabela "Threshold detekcji":

**ZMIANA:** wszystkie "Δβ_ppE^(3PN)" → **"Δβ_ppE^(2PN-phase)"** lub
po prostu **"Δβ_ppE^(b=−1)"**. Wartości pozostają (są correct dla
b=-1):
- LIGO-O3: ~10⁻¹ (NIE ~3·10⁻¹ jak w aktualnym drafcie — to było 3PN-phase)
- LIGO-O5 single: ~10⁻¹
- LIGO-O5 stack 100 BBH: ~3·10⁻²
- ET-D single: ~10⁻³
- CE single: ~3·10⁻⁴

W §1 falsifier clause:

**ZMIANA:** "ppE 3PN coefficient β_ppE^(3PN)" → **"ppE 2PN-phase
coefficient β_ppE^(b=−1)"**.

### 4.2 Patch do [[SENSITIVITY_BACK_OF_ENVELOPE.md]] §4.2

**ZMIANA:** w tabeli "TGP vs bounds":
- Wszystkie "β_ppE_5σ_bound" odniesienia → "β_ppE_5σ_bound^(2PN-phase, b=−1)".
- Liczby pozostają (zostały już w §3.1 podane jako line "2PN").

### 4.3 Patch do [[README.md]] problem #4

Tekst aktualny:
> "U³ to człon 3PN (post-Newtonian rzędu 3). LIGO obecnie ma 3.5PN
> waveformy. LIGO 3G będzie miał ~4PN+."

**ZMIANA:**
> "U³ deviation w `g_tt` to **2PN-phase** correction w stationary phase
> approximation waveformu (równoważnie: 3PN-energy w konwencji
> PPN/Will). LIGO obecnie ma 3.5PN-phase complete (Mishra 2016),
> LIGO 3G celuje w 4PN-phase. Dla 2PN-phase coefficient (M911-P1
> b_ppE=−1), aktualne LIGO-O3 bounds są ~10⁻¹; ET-D + CE osiągną
> ~10⁻³ na single event, ~10⁻⁴ na stack. Wymaga spójnego rachunku
> 2PN-phase w TGP — co znaczy: pełny dE/dt + dE/dr + SPA chain
> w M9.1''. **Pełne 2PN-phase matching jest dużą pracą** ale
> wykonalne (cykl `op-ppE-mapping/` Phase 1)."

### 4.4 Patch do [[NEEDS.md]] N7 (regime ważności)

**DODAĆ** do N7: "z explicit konwencją PN counting (zob.
[[CONVENTION_DECISION.md]]: PHASE convention adoptowana)."

## §5 — Kiedy patche zaaplikować

**OPCJA A (zalecane):** zaaplikować patche **TERAZ**, w ramach
T01 v2.1 (tej iteracji audytu). To czyni T01 internally consistent
i unambiguous przed wklejeniem falsifier do rejestru.

**OPCJA B:** odroczyć patche do zamknięcia cyklu `op-ppE-mapping/`
(Phase 1 zwalida wybór). To zmniejsza ryzyko, że patche będą
kolejny raz aktualizowane.

**Rekomendacja audytu:** **OPCJA A** — patche są mechaniczne,
nie wprowadzają nowych założeń, tylko klarują nomenclature. Cykl
`op-ppE-mapping/` może odkryć że b_ppE jest *różne* od −1
(efektywnie: cosmological constant-like behavior, lub nowy
mode), ale wtedy korekta będzie *na danych*, nie na konwencji.

## §6 — Cytaty i źródła konwencji

### Konwencja phase-PN (rekomendowana)

- Cutler & Flanagan, Phys. Rev. D **49**, 2658 (1994) — kanonicalny
  paper SPA inspiral.
- Blanchet, *Living Rev. Relativ.* **17**, 2 (2014), §6 — phase-PN
  dictionary tables.
- Yunes & Pretorius, Phys. Rev. D **80**, 122003 (2009),
  §II — ppE definition w phase-PN.
- Buonanno et al., Phys. Rev. D **80**, 084043 (2009) — TaylorF2
  phase-PN reference.

### Konwencja energy-PN (PPN historical)

- Will, *Theory and Experiment in Gravitational Physics* (CUP 2nd ed. 2018),
  Ch. 4 — PPN formalism.
- Will & Nordtvedt, ApJ **177**, 757 (1972) — original PPN expansion.
- Damour & Esposito-Farèse, Class. Quantum Grav. **9**, 2093 (1992) —
  PPN scalar-tensor.

### Mapping między konwencjami

- Sampson, Yunes, Cornish, Phys. Rev. D **88**, 064056 (2013),
  §II.B — explicit dictionary energy-PN ↔ phase-PN ↔ ppE b.
- Mirshekari & Will, Phys. Rev. D **87**, 084070 (2013), Tabela I —
  PPN parameters → ppE mapping.

## Cross-references

- [[README.md]] — diagnoza T01 (zaktualizować problem #4 zgodnie z §4.3)
- [[NEEDS.md]] — N7 (zaktualizować zgodnie z §4.4)
- [[FALSIFIER_STATEMENT_DRAFT.md]] §1 — patch zgodnie z §4.1
- [[PPN_TO_PPE_MAPPING.md]] — dictionary; brak zmian (już dual-presentation)
- [[SENSITIVITY_BACK_OF_ENVELOPE.md]] §4.2 — patch zgodnie z §4.2
- [[../../research/op-newton-momentum/M9_1_pp_P1_results.md]] — pochodzenie
  liczb 5/6, 23/12, …; energy-convention używana, MOŻE zostać
  zachowana w research-folderze (peer-review T01 paper przepisuje
  do phase convention przy submisji)
- [[../../PREDICTIONS_REGISTRY.md]] — target dla wpisu (po patch FALSIFIER_STATEMENT_DRAFT)
