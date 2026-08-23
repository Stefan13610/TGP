---
title: "Phase 0 — LOCK: plombowanie reguły A/B/C + read-lock źródeł + balance gate (op-Kgeo-from-D-uniqueness)"
date: 2026-06-26
cycle: op-Kgeo-from-D-uniqueness-2026-06-26
parent: "[[./README.md]]"
phase: 0
status: LOCKED
pre_registration_date: "2026-06-26"   # IMMUTABLE — plomba PRZED rachunkiem K_geo^(0)
anti_lakatos_lock: ARMED
---

# Phase 0 — LOCK (op-Kgeo-from-D-uniqueness)

> **Cel cyklu (z README §0.1):** wyznaczyć **K_geo niezależnie od poziomu-0**
> (`thm:D-uniqueness` + geometria całkowania po kącie rury kolorowej), **BEZ** wstrzykiwania
> α_s ani mas kwarków (circularity guard T-anti-circ), i zastosować value-blind regułę A/B/C
> do stosunku `R := K_geo^(0)·m_sp²/(π·Φ₀²)`. Native output = K_geo (bezwymiarowy); pochodny
> observable = α_s(M_Z) — **genuine predykcją staje się TYLKO przy werdykcie (A) DERIVED**.

---

## §0 — Autoryzacja + WIP gate

| Warunek wejścia (z README STATUS) | Stan |
|---|---|
| Explicit autoryzacja użytkownika („rozpocznij cykl") | ✅ 2026-06-26 |
| Wolny WIP slot (CYCLE_LIFECYCLE WIP-limit 5) | ✅ po sesji #46 propagacja #36-43 KOMPLETNA; kolejka user wyczerpana (STATE #43-#46); realny WIP wyczyszczony |
| Werdykt #43 IMMUTABLE — zakaz re-litygacji | ✅ cykl NIE re-litiguje #43; atakuje 1 zlokalizowaną lukę (K_geo niezależne) |
| Reguła + progi PRZED rachunkiem | ✅ plomba poniżej §3, immutable timestamp |

**Transition:** `folder_status: parking → active` (README YAML zaktualizowane).

---

## §1 — Read-lock źródeł (read-only; status z manuskryptu, NIE życzenie)

### §1.1 — `dodatekX_quark_sector.tex` — łańcuch i postulat
- **prop:X-A-from-tube-tension** (l.954-1048), status **AN+NUM [warunkowy]**.
- **eq:X-sigma-physical** (l.962-971): `σ_phys = (K_geo·m_sp²)·π(C_F α_s/(π Φ₀))² = [K_geo·m_sp²/(π Φ₀²)]·C_F²α_s²`.
- **eq:X-K-msp-hypothesis** (l.975-979): `K_geo·m_sp² ≝ π·Φ₀²` (oznaczone `=!`, hipoteza) ⟹
  kasacja Φ₀ ⟹ **eq:X-A-equals-CF2-alphas2** (l.981-989): `𝒜 = a_Γ/φ = C_F²α_s²(M_Z)`.
- **Weryfikacja numeryczna (l.992-1006):** `C_F²α_s² = 0,02518` vs `𝒜_TGP = a_Γ/φ = 0,02472`
  → ratio **1,019 (1,9%)**; pełne całkowanie poprawia do **1,4%**. ⚠ **To używa α_s jako INPUT**
  → dokładnie circularity, której Phase 1 MUSI uniknąć (R liczone bez α_s).
- **Warunek (l.1032-1039):** „K_geo jest parametrem geometrycznym substratu (tw. D-uniqueness),
  m_sp²=β∼Φ₀ na dużych skalach, czynnik π z całkowania po kącie rury. **Pełne wyprowadzenie
  wymaga zamknięcia coarse-grainingu (CG-1/CG-3).**"
- **Następne kroki (l.1042-1047):** (1) wyprowadzić eq:X-K-msp z aksjomatów poziomu-0;
  (2) bieganie α_s do ∼1 GeV; (3) profil rury poza gaussem (K₀-Bessel).

### §1.2 — `thm:D-uniqueness` (sek08_formalizm.tex l.962-1000) — KLUCZOWA OBSERWACJA SCOPE
Twierdzenie: wśród operatorów `L_kin = -½K(φ)(∇φ)²` warunki (C1) stałość α, (C2) K(0)=0,
(C3) `K(φ)=K_geo·φ⁴` wyznaczają **jednoznacznie α=2** oraz formę `K(φ)=K_geo·φ⁴`.

> **OBSERWACJA (value-blind, NIE werdykt):** Krok 2 dowodu (l.1023-1035) pokazuje
> `K(φ)=C·φ^{2α}, C>0` — twierdzenie ustala **FORMĘ** (φ⁴) i **α=2**, ale **`K_geo` to dowolny
> dodatni prefaktor C; jego WARTOŚĆ nie jest wyznaczona przez D-uniqueness**. Potwierdza to
> `status_map` l.72: `K(φ)=K_geo·φ⁴` ⇒ „selekcja w klasie konforemnej C1-C3 (na gęstości),
> **NIE derywacja z substratu**". ⟹ Ścieżka 2 (README §2) napotyka strukturalną przeszkodę
> już na wejściu: D-uniqueness sama nie pinuje K_geo. Czy K_geo^(0) da się wyznaczyć z
> *geometrii rury* (Ścieżka 1) lub *CG* (Ścieżka 3) niezależnie — to zadanie Phase 1.
> **Sygnał a priori w stronę (C), ale NIE przesądzony** (reguła value-blind rozstrzyga rachunkiem).

### §1.3 — `status_map.tex` — obstrukcja CG (most Γ→Φ)
- **CG-1/CG-3** (l.1323-1329): „Brak mostu Γ→Φ" — prop:Q-weak-continuum + prop:Q-alpha-from-phi-squared;
  status **[SZKIC], nie pełne domknięcie**.
- **CG-3 / ex200** (l.1514-1517): blokowe uśrednianie Isinga; `Φ_B≥0` OK, ale **α_eff
  niedostatecznie zbieżny na małych siatkach; ex200 4/8 PASS** (T2,T3,T5,T7 FAIL — wymaga większych L).
- **ex202 / sektor kwarkowy** (l.1519-1523): `𝒜` uniwersalne (1%), ale **skalowanie
  `𝒜∼√σ/Φ₀` nie zamknięte; ex202 7/8 PASS (T6 FAIL: σ_TGP)**.
- **K_geo conformal** (l.72, l.563, l.1597-1624): `K_geo·φ⁴`, `α(α-1)=2 ⇒ α=2`; kanoniczne
  pole `Φ̃=√K_geo·φ³/3`.

### §1.4 — Skrypty referencyjne (do Phase 1, NIE uruchamiane w Phase 0)
- `research/nbody/examples/ex200_lyapunov_beta_scan_yukawa.py` + `tooling/scripts/ex200_cg3_continuum_verification.py` (α_eff convergence, CG-3).
- `research/nbody/examples/ex202_multipole_vs_feynman_triple_overlap.py` + `tooling/scripts/ex202_quark_m0_derivation.py`.
- `tooling/scripts/color_tube_variational_tgp.py` (σ̂=πA² zgodne 0,5%; profil rury) + `color_tube_advanced_tgp.py` (Bessel K₀ vs gauss).

### §1.5 — Źródło luki (#43)
`op-A-derivation-from-CG-2026-06-25/Phase_FINAL_close.md`: werdykt **POSTULATE-CONDITIONAL**
(4/4). Łańcuch L1 σ̂=πA² (ansatz; Bessel K₀ daje ×0,3-0,6 → nie unikalny) · L2 A_color
derywowane · L3 m_sp²=γ · L4 σ_phys definicyjne · **L5 = eq:X-K-msp-hypothesis (jedyny
load-bearing postulat)**. L5 redukuje się do mostu Γ→Φ. „Następny krok (§5): wyprowadzenie
K_geo z poziomu-0 D-uniqueness" — **to jest dokładnie ten cykl**.

---

## §2 — TGP-native check (z README §0.3) — potwierdzenie pre-Phase-1
Q1-Q8 wypełnione w README §0.3. Kluczowe:
- **Q4 circularity guard (KRYTYCZNY):** K_geo^(0) liczone **bez α_s i bez mas kwarków**;
  T-anti-circ obowiązkowy w Phase 1 (zero PDG α_s w wejściu; ratio numeryczny §1.1 l.992-1006
  NIE może być użyty jako wejście — to byłaby cyrkularność).
- **Q8:** cykl NIE re-litiguje #43; wyznacza K_geo^(0) (czego #43 NIE robił).
- **ASK-RULE (NATIVE_PATTERNS §1):** Phase 1 dotknie Φ_eq/profil rury → Pattern 2.1 relevant;
  red flag do pilnowania: NIE linearyzować D_kin do Klein-Gordon bez weak-field justification
  (Warning A/B Pattern 2.1). Jeśli Phase 1 sięgnie po std-physics dla K_geo bez native ścieżki
  → Trigger A/C → STOP-and-ASK.

---

## §3 — 🔒 PLOMBA reguły falsyfikacji (value-blind, IMMUTABLE 2026-06-26)

**Metryka rozstrzygnięcia:** `R := K_geo^(0)·m_sp²/(π·Φ₀²)`, gdzie `K_geo^(0)` wyprowadzone
niezależnie od poziomu-0 (Ścieżka 1/2/3), bez α_s ani mas kwarków.

| Werdykt | Warunek na R | Konsekwencja |
|---|---|---|
| **(A) DERIVED** | `R ∈ [0,95; 1,05]` (≤5%) | most 𝒜=C_F²α_s² domknięty od poziomu-0 ⇒ α_s(M_Z) staje się **genuine first-principles** predykcją; α_s przechodzi TRADED(#42)→DERIVED (zmiana ledger parameter-counting) |
| **(B) REFUTED-BRIDGE** | `R ∉ [0,80; 1,25]` (>25% w którąkolwiek stronę) | postulat eq:X-K-msp-hypothesis **sfalsyfikowany** ⇒ uruchamia PR-025 forward falsifier (b); 𝒜 NIE jest mostem do α_s |
| **(C) POSTULATE-CONFIRMED** | `R ∈ [0,80;0,95)∪(1,05;1,25]` **LUB** K_geo^(0) nieoznaczalne bez domknięcia CG-1/CG-3 (ex200 α_eff niezbieżny przy dostępnym L) | status #43 POSTULATE-CONDITIONAL **RATYFIKOWANY** — K_geo irreducibly conditional |

**Anti-moving-goalposts (BINDING):** progi 5%/25% są immutable. Jakakolwiek zmiana progu
po rozpoczęciu Phase 1 = **HALT-B**. Reguła jest dwustronna (anti-Lakatos Q7): pozytywny (A),
negatywny (B), warunkowy (C) — żaden nie faworyzowany.

**Sprzężenie zewnętrzne:** (B) → PR-025 forward (b) `K_geo·m_sp²≠π·Φ₀² ⟹ 𝒜≠C_F²α_s²`.
Propagacja (dodatekX rem:X-rescope-43 / PREDICTIONS_REGISTRY / PR-025) = **user-gated**, per A/B/C.

---

## §4 — Balance gate (budżet nowych stałych)

- **Budżet nowych stałych = 0.** K_geo **NIE jest nową stałą** — ma być POCHODNĄ z poziomu-0
  (D-uniqueness / geometria rury / CG). Jeśli Phase 1 wykaże, że K_geo jest *wolnym parametrem*
  (prefaktor C niewyznaczony, §1.2) → to jest **wynik (C)**, nie wprowadzenie nowej stałej.
- m_sp², Φ₀, π, C_F: istniejące wielkości poziomu-0 / geometrii grupy; nie nowe.
- **Zero hardcoded** w Phase 1; ≥1 niezależny solver; T-anti-circ (zero α_s w wejściu).

---

## §5 — Plan Phase 1 (do osobnego uruchomienia)

Trzy ścieżki (README §2), w kolejności mocy:
1. **Geometryczna** — całkowanie po kącie rury: czy czynnik `π·Φ₀²` wychodzi z geometrii
   przekroju (Bessel K₀ vs gauss; #43 L1: Bessel ×0,3-0,6 ⟹ test unikalności ansatzu).
2. **D-uniqueness** — czy K_geo wynika z `thm:D-uniqueness` niezależnie od m_sp (§1.2: a priori
   przeszkoda — twierdzenie fixuje formę φ⁴, nie wartość prefaktora C).
3. **Lattice/CG** — ex200 z większym L: czy α_eff zbiega (obecnie 4/8). Niezbieżność przy
   dostępnym L → silny sygnał (C).

Następnie: wyliczyć `R`, zastosować plombę §3 (value-blind), i — jeśli (A) — policzyć
α_s(M_Z) first-principles vs PDG (teraz genuine, nie consistency).

---

## §6 — Anti-Lakatos LOCK (ARMED)
✓ Reguła + progi zaplombowane PRZED rachunkiem (immutable 2026-06-26). ✓ Read-lock źródeł
read-only (status z manuskryptu). ✓ Circularity guard armed (zero α_s/mas kwarków w K_geo^(0)).
✓ Reguła dwustronna (A/B/C). ✓ Budżet nowych stałych = 0. ✓ #43 IMMUTABLE, zero re-litygacji.
✓ Obserwacja §1.2 (D-uniqueness nie pinuje K_geo) zapisana value-blind, NIE jako przesądzony werdykt.

**Sign-off:** Claudian @ 2026-06-26 — Phase 0 LOCKED. Następny krok: Phase 1 (osobne uruchomienie).

## Cross-references
- [[./README.md]] · [[../op-A-derivation-from-CG-2026-06-25/Phase_FINAL_close.md]] (#43 — źródło luki)
- [[../../partial_proofs/quark_sector/dodatekX_quark_sector.tex]] prop:X-A-from-tube-tension, eq:X-K-msp-hypothesis
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] thm:D-uniqueness (l.962-1000)
- [[../../core/_meta_latex/status_map.tex]] CG-1/CG-3 [SZKIC], ex200 4/8, ex202 7/8
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] PR-025 (forward falsifier (b))
