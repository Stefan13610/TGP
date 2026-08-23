---
title: "Phase 0 — LOCK: plombowanie reguły A/B/C (α_eff FSS) + read-lock + audyt ex200/ex202 vs CG34 + balance gate (op-CG-alpha-eff-convergence)"
date: 2026-06-26
cycle: op-CG-alpha-eff-convergence-2026-06-26
parent: "[[./README.md]]"
phase: 0
status: LOCKED
pre_registration_date: "2026-06-26"   # IMMUTABLE — plomba PRZED nowym rachunkiem FSS α_eff
anti_lakatos_lock: ARMED
---

# Phase 0 — LOCK + Audyt Fazy A (op-CG-alpha-eff-convergence)

> **Cel cyklu (README §0.1):** rozstrzygnąć value-blind, czy efektywny wykładnik kinetyczny
> `α_eff` blokowo-uśrednionego substratu zbiega do **2** (continuum/FSS), czy obstrukcja jest
> **strukturalna** (niespójność `α_eff = s−1` z CG34 #31). Native output = α_eff (bezwymiarowy);
> pochodny status = α=2 jako **substrate-DERIVED** — tylko przy werdykcie (A).

---

## §0 — Autoryzacja + WIP gate

| Warunek wejścia | Stan |
|---|---|
| Explicit autoryzacja użytkownika („działaj z Fazą A") | ✅ 2026-06-26 |
| Wolny WIP slot (po #48 op-Kgeo CLOSED, slot zwolniony) | ✅ STATE #48 §WIP |
| Werdykty #31/#43/#48 IMMUTABLE — zakaz re-litygacji | ✅ cykl NIE re-litiguje; rozstrzyga 1 zlokalizowaną sprzeczność (ex200 ≈2 vs CG34 s−1) |
| Reguła + progi PRZED rachunkiem | ✅ plomba §3, immutable timestamp |

**Transition:** `folder_status: (new) → active`.

---

## §1 — Read-lock źródeł (read-only; status z manuskryptu/skryptów, NIE życzenie)

### §1.1 — `status_map.tex` — obstrukcja CG (most Γ→Φ)
- **CG-3 / ex200** (l.1514-1517): blokowe uśrednianie Isinga; `Φ_B≥0` OK, ale **α_eff niedostatecznie
  zbieżny na małych siatkach; ex200 4/8 PASS** (T2,T3,T5,T7 FAIL — wymaga większych L).
- **ex202 / sektor kwarkowy** (l.1519-1523): `𝒜` uniwersalne (1%), ale **skalowanie `𝒜∼√σ/Φ₀` nie
  zamknięte; ex202 7/8 PASS (T6 FAIL: σ_TGP)**.
- **K_geo conformal** (l.72): `K(φ)=K_geo·φ⁴`, „selekcja w klasie konforemnej (na gęstości), **NIE
  derywacja z substratu**".

### §1.2 — `op-CG34-continuum-closure` (#31, 2026-06-14) — KLUCZOWY PREDECESSOR
- **CG-3 = ZAMKNIĘTY NUM** (5/5): homogenizacja `Φ_B→Φ` w `H¹`, normy Cauchy {1,34;0,38;0,11},
  rate `L_B^{−1,43}`. ⟹ **continuum limit pola istnieje** (homogenizacja działa).
- **CG-4 = PARTIAL:** ZAMKNIĘTE — `c*>0` stabilne (red-flag c*→0 = artefakt `⟨|∇Φ|²⟩`), `β=γ`,
  `K_hom-forma=K_IR`. **DO-PINNIĘCIA — `α=2 ↔ K(φ)`:** samodzielna algebra `α_eff = s−1`
  (`Z∝φ^{2s}`); **substrat `Z∝φ²` (s=1) ⇒ α_eff = 0**, nie 2 ⟹ **niespójność lematu A3**
  (premisa `s=1`, konkluzja `s=0`/`K_1∼1/Φ`), zgłoszona do rdzenia, **nie sfabrykowano α=2**.
- **Substrat `-J(φ_iφ_j)²` PATOLOGICZNY dla pure-MC:** λ małe → runaway `⟨ρ⟩~1600`; λ duże →
  frozen `ξ→0`; **brak okna scale-separated**. ⟹ Trigger B (tech-debt) armed: Phase 1 NIE dziedziczy
  tego substratu bezkrytycznie.

### §1.3 — `op-Kgeo-from-D-uniqueness` (#47-48, 2026-06-26) — bezpośredni poprzednik
- Werdykt **(C) POSTULATE-CONFIRMED** (9/9): `thm:D-uniqueness` fiksuje FORMĘ `K(φ)=K_geo·φ⁴`+α=2,
  ale **K_geo to wolny prefaktor** (C3 go *definiuje*, nie liczy). K_geo dołącza do α=2/c₀ jako
  **irreducibly conditional pending CG**. ⟹ ten cykl atakuje **dokładnie ten sam korzeń** od strony
  numerycznej (czy substrat *generuje* α=2/`Z∝φ⁶`, czy nie).

### §1.4 — Inspekcja kodu (read-only, baseline)
- **`ex200_cg3_continuum_verification.py`:** `run_tests()` (l.336-340) używa **`test_L = [16]`**
  (POJEDYNCZY L!), `test_b = [2,4]`, `n_therm=500`, `n_measure=100`. **T3** (`α_eff≈2`, l.379-389):
  estymator = **log-log fit gęstości gradientowej vs Φ** (`extract_kinetic_coefficient`, l.134-170),
  cut `r²>0,3`, tolerancja **`|α−2|<1,5`**. ⟹ „α_eff≈2 PASS" jest **przetolerowany** i oparty o
  estymator, który CG34 zdiagnozował jako **artefakt-prone** (ten sam `⟨|∇Φ|²⟩` co red-flag c*→0).
  **FSS niemożliwe przy single-L** — sprzeczne z deklaracją „continuum limit".
- **`ex202_quark_m0_derivation.py`:** baseline uruchomiony **2026-06-26** → **7/8 PASS**, **T6 FAIL**
  (`A~√σ/Φ₀`: A_pred=0,0176 vs A_obs≈2,47·10⁻⁵ → odchylenie ~712×). Potwierdza `status_map`:
  most `𝒜~√σ/Φ₀` **niezamknięty**. (T4/T5 shifted-Koide 2/3 = dopasowanie m₀, nie predykcja skali.)

---

## §2 — AUDYT REKONCYLIACYJNY (Faza A — główne ustalenie)

**Dwa sprzeczne sygnały dla α_eff:**

| Źródło | Wynik α_eff | Metoda | Wiarygodność |
|---|---|---|---|
| `ex200` T3 (single-L=16) | „≈2 PASS" | log-log gęstości gradientowej, tol 1,5, r²>0,3 | ⚠ NISKA — single-L (brak FSS), estymator artefakt-prone (CG34), tolerancja 1,5 |
| `op-CG34` algebra | `α_eff = s−1 = 0` (substrat s=1) | samodzielna algebra `Z∝φ^{2s}` | 🟢 WYŻSZA — analityczna, niezależna, zgłoszona uczciwie |

**Rozstrzygnięcie audytu (value-blind, NIE werdykt cyklu):** sygnały **NIE są równoważne** — `ex200`
mierzy zaszumiony estymator przy jednym L z tolerancją tak szeroką (1,5), że **nie odróżnia α=2 od
α=0,5..3,5**; algebra CG34 jest precyzyjna i wskazuje `s−1`. ⟹ **Obstrukcja jest strukturalna
(skala/wykładnik `Z(φ)`), NIE „mały L".** „Mały L" to objaw (FSS niezbieżne), nie przyczyna.

**Konsekwencja dla designu Phase 1:** naiwne „uruchom ex200 przy większym L" **nie wystarczy** i
ryzykuje powtórzenie artefaktu. Phase 1 MUSI: (a) zamienić estymator na structure-factor (CG34),
(b) zrobić **prawdziwe FSS** (≥3 L), (c) zmierzyć **wykładnik `s` z `Z(φ)`** bezpośrednio (nie tylko
α_eff z gradientów), (d) użyć niepatologicznego substratu.

> **Obserwacja a priori (value-blind, NIE przesądza):** CG34 `s=1` + ex200 przetolerowany
> ⟹ sygnał w stronę **(B) REFUTED-SUBSTRATE** (α_eff→0) **lub (C) INDETERMINATE** (patologia substratu).
> (A) wymagałoby emergentnego `Z∝φ⁶` — możliwe tylko jeśli coarse-graining generuje nietrywialny
> wykładnik anomalny. Reguła §3 rozstrzyga rachunkiem, nie a priori.

---

## §3 — 🔒 PLOMBA reguły falsyfikacji (value-blind, IMMUTABLE 2026-06-26)

**Metryka rozstrzygnięcia:** `ᾱ := α_eff⁽⁰⁾` — finite-size-ekstrapolowane efektywne α_eff z
blokowo-uśrednionego substratu, estymator **structure-factor** (NIE log-log ex200), **bez α_s i bez
mas kwarków** (T-anti-circ), przy **≥3 rozmiarach L** z trendem FSS; `R²_FSS` = jakość dopasowania
ekstrapolacji.

| Werdykt | Warunek | Konsekwencja |
|---|---|---|
| **(A) DERIVED-SUBSTRATE** | `\|ᾱ−2\| ≤ 0,3` **∧** `R²_FSS ≥ 0,7` (≥3 L) | substrat GENERUJE α=2 (emergentne `Z∝φ⁶`, s=3); niespójność A3 (#31) ZREKONCYLIOWANA; α=2: aksjomat→**substrate-DERIVED**; zasila K_geo (#47-48)/𝒜 (#43); zmienia ledger #42 |
| **(B) REFUTED-SUBSTRATE** | `\|ᾱ−2\| ≥ 1,0` **∧** `R²_FSS ≥ 0,7` | substrat NIE generuje α=2; niespójność A3 (#31) **POTWIERDZONA realna**; α=2 ściśle aksjomatyczne-na-gęstości; ścieżka substratowa **CLOSED-NEGATIVE** (genuine negatyw); rekomendacja: reframe dodatekQ2 A3 |
| **(C) INDETERMINATE** | `ᾱ ∈ [1,0;1,7)∪(2,3;3,0]` **LUB** `R²_FSS < 0,7` **LUB** patologia substratu (brak okna scale-separated) | α_eff⁽⁰⁾ nieoznaczalne ⇒ **ratyfikacja status quo** (α=2 aksjomatyczne; K_geo POSTULATE-CONDITIONAL #48). A priori najbardziej prawdopodobne (CG34) |

**Anti-moving-goalposts (BINDING):** progi `|ᾱ−2|`: 0,3 / 1,0 oraz `R²_FSS`: 0,7 są **immutable**.
Jakakolwiek zmiana progu po rozpoczęciu Phase 1 = **HALT-B**. Reguła trójwartościowa, dwustronna+
(anti-Lakatos Q7): pozytywny (A), negatywny (B), warunkowy (C) — żaden nie faworyzowany; (B) jest
realnie osiągalny i jawny.

**Sprzężenie zewnętrzne:** (A) → thm:alpha2 status + ledger #42 + most K_geo/𝒜; (B) → dodatekQ2 A3
reframe + CLOSED-NEGATIVE; (C) → status quo #48. Propagacja = **user-gated**, per A/B/C.

---

## §4 — Balance gate (budżet nowych stałych)

- **Budżet nowych stałych = 0.** α_eff jest **MIERZONE** z substratu, nie wprowadzane. `s` (wykładnik
  `Z(φ)`) to własność modelu substratu, nie nowa stała teorii.
- Jeśli Phase 1 wykaże, że α_eff⁽⁰⁾ nieoznaczalne (patologia/FSS) → to **wynik (C)**, nie nowa stała.
- **Zero hardcoded** w Phase 1; ≥1 niezależny solver; T-anti-circ (zero α_s/mas kwarków; zakaz
  „odwrotnego" dopasowania do α=2).

---

## §5 — Plan Phase 1 (do osobnego, user-gated uruchomienia)

**Cel:** zmierzyć `ᾱ = α_eff⁽⁰⁾` + `R²_FSS` i zastosować plombę §3.

1. **Estymator (KRYTYCZNA poprawka vs ex200):** structure-factor `Z(k) = ⟨|Φ̃_B(k)|²⟩`, wykładnik
   z niskoimpulsowego skalowania propagatora `G(k) ~ 1/(Z k² + m²)` → `K_eff`, **NIE** log-log
   gęstości gradientowej (CG34: artefakt-prone, źródło c*→0).
2. **Bezpośredni pomiar `s`:** `Z(φ)` jako funkcja lokalnego `φ` (binowanie po φ), fit `Z∝φ^{2s}`;
   wynik `α_eff = s−1`. To rozstrzyga niespójność A3 wprost (s=1 ⇒ B; s=3 ⇒ A).
3. **FSS:** L ∈ {16,24,32} (opcj. 48 jeśli zasoby), skan b ∈ {2,4,8}; ekstrapolacja `α_eff(L)→α_eff⁽⁰⁾`,
   pomiar `R²_FSS`.
4. **Substrat niepatologiczny:** wybrać model ze scale-separated window (CG34: `-J(φ_iφ_j)²` runaway/
   frozen) — np. standardowy `φ⁴` Z₂ blisko `κ_c≈0,189` (CG34 Phase 1), LUB czysto analityczna
   homogenizacja `Z(φ)` jeśli MC niestabilne.
5. **Wydajność:** pure-Python ex200 za wolny dla large-L (potrójna pętla); wymaga wektoryzacji
   (numpy checkerboard) lub `numba`; flaga zasobowa.
6. **T-anti-circ ENFORCED:** zero α_s/mas kwarków; zakaz wymuszania α=2; ≥1 niezależny solver.

Następnie: zastosować plombę §3 (value-blind) → werdykt A/B/C → propagacja user-gated.

---

## §6 — Anti-Lakatos LOCK (ARMED)
✓ Reguła + progi zaplombowane PRZED nowym rachunkiem (immutable 2026-06-26). ✓ Read-lock źródeł
read-only. ✓ Audyt rekoncyliacyjny zapisany value-blind (ex200 przetolerowany vs CG34 s−1), bez
przesądzania werdyktu. ✓ Circularity guard armed (zero α_s/mas kwarków; zakaz odwrotnego dopasowania).
✓ Reguła trójwartościowa A/B/C, (B) realnie osiągalny i jawny. ✓ Budżet nowych stałych = 0. ✓ #31/#43/#48
IMMUTABLE, zero re-litygacji. ✓ Trigger B (tech-debt substratu) armed — CG34 patologia jawna, nie obejście.
✓ Estymator poprawiony (structure-factor) z explicit uzasadnieniem (artefakt ex200), nie ukryty.

**Sign-off:** Claudian @ 2026-06-26 — Phase 0 LOCKED + Audyt Fazy A WYKONANY. Następny krok: Phase 1
FSS (osobne, user-gated uruchomienie; §5).

## Cross-references
- [[./README.md]] · [[../op-CG34-continuum-closure-2026-06-14/Phase_FINAL_close.md]] (#31 — niespójność A3)
- [[../op-Kgeo-from-D-uniqueness-2026-06-26/Phase_FINAL_close.md]] (#47-48 — C POSTULATE-CONFIRMED)
- [[../../core/_meta_latex/status_map.tex]] CG-3/CG-4 [SZKIC], ex200 4/8, ex202 7/8
- skrypty: `tooling/scripts/ex200_cg3_continuum_verification.py`, `tooling/scripts/ex202_quark_m0_derivation.py`
- [[../op-uv-as-ngfp/]] (NGFP — opcjonalna L2 reduction)
