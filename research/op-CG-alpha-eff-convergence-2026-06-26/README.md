---
title: "op-CG-alpha-eff-convergence — czy α_eff blokowo-uśrednionego substratu zbiega do 2 (continuum/FSS), czy obstrukcja jest strukturalna (niespójność α_eff=s−1 z CG34/#31)?"
date: 2026-06-26
type: research-cycle
folder_status: closed-resolved   # CLOSED 2026-06-26: (B) REFUTED-SUBSTRATE. Faza A LOCK + Phase 1 FSS + FINAL w 1 sesji — patrz Phase_FINAL_close.md
claim_status: "structural — (B) REFUTED-SUBSTRATE: substrat NIE generuje α=2 (α_eff=−1/2 exact; e_inf=−0.12 MC); α=2 aksjomatyczne-na-gęstości"
parent: "[[../../core/_meta_latex/status_map.tex]]"

# ============== KICKOFF CONTRACT (mandatory post-2026-05-10) ==============
contract:
  # --- L1: Native (MANDATORY) ---
  L1_native:
    output_observable: "α_eff (bezwymiarowy efektywny wykładnik kinetyczny substratu: K_eff(φ) ∝ φ^{2α_eff}, równoważnie α_eff = s−1 dla Z(φ) ∝ φ^{2s}). Brak bezpośredniego instrumentu (structural). Pochodny status: α=2 — czy 'aksjomatyczne na gęstości' (status quo) czy 'substrate-DERIVED'."
    measurement_instrument: "Lattice Monte Carlo substratu Z₂ (φ⁴/Ising 3D) + blokowe uśrednianie Φ_B = ⟨ŝ²⟩_blok; estymator structure-factor (CG34 Phase 1 — NIE log-log fit gęstości gradientowej z ex200, oznaczony jako artefakt-prone). Cross-check przez downstream α=2 (sektor leptonowy/metryka)."
    native_coefs_constrained:
      - "α (wykładnik kinetyczny TGP; manuskrypt: K(φ)=K_geo·φ⁴ ⇒ α=2)"
      - "s (wykładnik skalowania Z(φ)∝φ^{2s}; CG34: α_eff=s−1; substrat Φ=⟨ŝ²⟩ ⇒ s=1 ⇒ α_eff=0)"
    falsification_rule: >
      "Niech ᾱ := finite-size-extrapolowane α_eff⁽⁰⁾ z blokowo-uśrednionego substratu, mierzone
      estymatorem structure-factor (NIE artefakt-prone log-log ex200), BEZ α_s i BEZ mas kwarków
      (T-anti-circ), przy ≥3 rozmiarach siatki L z trendem FSS. Wtedy (progi value-blind, plombowane):
      (A) DERIVED-SUBSTRATE — jeśli |ᾱ−2| ≤ 0,3 ORAZ R²_FSS ≥ 0,7 przy ≥3 L: substrat GENERUJE α=2
          (emergentne Z∝φ⁶, s=3); niespójność lematu A3 (#31) ZREKONCYLIOWANA; α=2 przechodzi
          'aksjomat na gęstości' → 'substrate-DERIVED'; zasila K_geo (#47-48)/𝒜 (#43);
      (B) REFUTED-SUBSTRATE — jeśli |ᾱ−2| ≥ 1,0 (w szczególności ᾱ≈0 z algebry s−1) ORAZ R²_FSS ≥ 0,7:
          substrat NIE generuje α=2; niespójność A3 (#31) POTWIERDZONA jako realna; α=2 pozostaje
          ściśle aksjomatyczne-na-gęstości, ścieżka substratowa CLOSED-NEGATIVE (genuine wynik
          negatywny — publikowalny); rekomendacja rdzenia: reframe lematu A3 (dodatekQ2);
      (C) INDETERMINATE — jeśli ᾱ ∈ [1,0;1,7)∪(2,3;3,0] (niejednoznaczne) LUB R²_FSS < 0,7 (FSS
          niezbieżny przy dostępnym L) LUB patologia substratu (brak okna scale-separated:
          runaway/frozen, CG34 §2): α_eff⁽⁰⁾ nieoznaczalne ⇒ ratyfikacja status quo (α=2 aksjomatyczne;
          K_geo POSTULATE-CONDITIONAL). A priori NAJBARDZIEJ prawdopodobne (CG34: substrat patologiczny; ex200 przetolerowany)."
    pre_registration_date: "2026-06-26"   # LOCKED immutable @ Phase 0 (PRZED jakimkolwiek nowym rachunkiem FSS); progi 0,3/1,0 value-blind zaplombowane

  # --- L2: Cross-framework reduction (OPTIONAL — last stage) ---
  L2_framework_reduction:
    target_frameworks: ["Wilsonian ERG / FRG fixed-point (K_IR/K_UV)", "Asymptotic-Safety NGFP (op-uv-as-ngfp)"]
    reduction_type: "not-attempted"   # ocena w Phase FINAL; fail NIE unieważnia L1
    validation_transfer: ""
    failure_disposition: "L1-stands"

  # --- L3: Falsification map (consistency) ---
  L3_falsification_map:
    - { bound: "ex200 CG-3 α_eff convergence (status_map l.1514-1517)", constrains: "α_eff via T3 (K∝Φ²)", window: "4/8 PASS @ L=16 single-size (T2,T3,T5,T7 FAIL)", status: "pending-FSS" }
    - { bound: "ex202 σ_TGP scaling (status_map l.1519-1523)", constrains: "𝒜~√σ/Φ₀ (powiązany most)", window: "7/8 PASS, T6 FAIL (σ_TGP, baseline 2026-06-26: A_pred/A_obs ~712×)", status: "open" }
    - { bound: "CG34 lemat A3 (op-CG34 Phase 3)", constrains: "α_eff=s−1; substrat s=1⇒α=0", window: "niespójność zgłoszona, nie rozstrzygnięta", status: "pending-resolution" }

# ============== END KICKOFF CONTRACT ==============

tgp_status:
  level: T3                 # głęboki korzeń strukturalny (substrat → makro)
  kind: audit               # Faza A = audyt + pre-rejestracja; Phase 1 = derivation (user-gated)
  output_type: structural   # α_eff bezwymiarowy strukturalny; max claim C (KICKOFF §2.2)
  core_compatibility: review-only
  may_edit_core: false
  has_needs_file: false
  has_findings_file: false
  exports_findings: false
  open_bridges: ["op-CG34-continuum-closure-2026-06-14", "op-Csigma-coarse-graining-2026-06-14", "op-uv-as-ngfp"]
  depends_on: ["#31 op-CG34-continuum-closure (CG-3 ZAMKNIĘTY NUM; CG-4 PARTIAL; niespójność A3)", "#47-48 op-Kgeo-from-D-uniqueness (C POSTULATE-CONFIRMED)"]
  impacts: ["dodatekQ CG-3/CG-4 status", "dodatekQ2 lemat A3", "thm:alpha2 / α=2 status (aksjomat vs derywacja)", "K_geo (#47-48) / 𝒜 (#43) wspólny korzeń", "c₀ (#37) wspólny mianownik UV/CG"]
  source_of_status: ["[[../op-CG34-continuum-closure-2026-06-14/Phase_FINAL_close.md]]", "[[../op-Kgeo-from-D-uniqueness-2026-06-26/Phase_FINAL_close.md]]"]
---

# op-CG-alpha-eff-convergence — Faza A: audyt + pre-rejestracja FSS α_eff (most CG-3/CG-4)

> **STATUS: 🟡 ACTIVE — Faza A (Phase 0 LOCK + audyt ex200/ex202 vs CG34 + plan Phase 1 FSS).**
> Cykl wskazany przez #43/#48 jako **najpilniejszy pojedynczy filar numeryczny** wspólnego korzenia
> UV/CG (α=2 / K_geo / 𝒜 / c₀). **Reframe Fazy A (kluczowy):** obstrukcja α_eff to **NIE 'mały L'**,
> lecz **strukturalna niespójność** — CG34 (#31) ustalił uczciwą algebrę `α_eff = s−1`
> (substrat `Z∝φ²`, s=1 ⇒ α_eff=0, **nie** 2), podczas gdy `ex200.run_tests()` biega tylko przy
> `L=16, b∈{2,4}` z tolerancją T3 = 1,5 — jego „α_eff≈2 PASS" jest **przetolerowanym sygnałem
> sprzecznym** z algebrą s−1. Faza A pre-rejestruje value-blind regułę **A/B/C** rozstrzygającą,
> czy ścieżka substratowa do α=2 jest osiągalna (A), sfalsyfikowana (B), czy nieoznaczalna (C).

## §0 — Cel + native-first contract

[CITE: meta/CYCLE_KICKOFF_TEMPLATE.md §1-§2; meta/PPN_AS_PROJECTION.md §3 (three-layer L1/L2/L3); meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md §1 (ASK-RULE); meta/M9_RESTRUCTURE_NOTE.md §3]

### §0.1 — Native observable target

#43 (op-A-derivation-from-CG) i #47-48 (op-Kgeo) zlokalizowały, że **trzy** kluczowe stałe TGP
(α=2, c₀, 𝒜→α_s) są `POSTULATE-CONDITIONAL`/aksjomatyczne — wszystkie czekają na domknięcie
**mostu coarse-grainingu Γ→Φ (CG-3/CG-4)**. Najbardziej skonkretyzowany, mierzalny filar tego mostu
to **zbieżność efektywnego wykładnika kinetycznego α_eff** blokowo-uśrednionego substratu:

$$K_{\rm eff}(\phi) \propto \phi^{2\alpha_{\rm eff}}, \qquad \alpha_{\rm eff} = s-1 \ \text{dla}\ Z(\phi)\propto\phi^{2s}.$$

Manuskrypt twierdzi `K(φ)=K_geo·φ⁴ ⇒ α=2` (thm:D-uniqueness — FORMA; thm:alpha2). Ale `status_map`
l.72 oznacza to jako **„selekcja w klasie konforemnej (na gęstości), NIE derywacja z substratu"**, a
CG34 (#31) ujawnił, że substrat `Φ=⟨ŝ²⟩` ma `s=1` ⇒ `α_eff=0` — sprzeczność z `α=2` (wymagałoby `s=3`).

**Cel Fazy A:** (1) zaudytować i zrekoncyliować dwa sprzeczne sygnały (`ex200` przetolerowany ≈2 vs
CG34 algebra `s−1=0`); (2) **pre-zarejestrować value-blind regułę A/B/C** rozstrzygającą natywnie
(lattice MC + structure-factor estimator, BEZ α_s/mas kwarków) czy α_eff⁽⁰⁾ zbiega do 2. Native
output = α_eff (bezwymiarowy strukturalny); pochodny status = α=2 jako derywacja substratowa **tylko
przy werdykcie (A)**.

### §0.2 — Pre-registered falsification rule (🔒 PLOMBED @ Phase 0 LOCK 2026-06-26)

Reguła trójwartościowa **(A) DERIVED-SUBSTRATE / (B) REFUTED-SUBSTRATE / (C) INDETERMINATE** —
patrz `contract.L1_native.falsification_rule` powyżej i `Phase0_LOCK.md §3`. Progi (|ᾱ−2|: 0,3 / 1,0;
R²_FSS: 0,7) **value-blind**, 🔒 **zaplombowane 2026-06-26 PRZED jakimkolwiek nowym rachunkiem FSS**
(immutable `pre_registration_date`). **Anti-moving-goalposts:** po Phase 1 zmiana progu = HALT-B.

**Sprzężenie zewnętrzne:** (A) → α=2 'aksjomat na gęstości' → 'substrate-DERIVED' (zmienia ledger
#42 + odblokowuje most K_geo/𝒜); (B) → potwierdza niespójność A3 (#31), uruchamia reframe dodatekQ2 +
ścieżka substratowa CLOSED-NEGATIVE; (C) → ratyfikacja status quo (α=2 aksjomatyczne; #47-48
POSTULATE-CONDITIONAL niezmieniony). Reguła dwustronna+ (anti-Lakatos): żaden wynik nie faworyzowany.

### §0.3 — TGP-native check (mandatory, pre-Phase-1)

- **Q1 (native, nie projection?):** ✅ α_eff to natywny wykładnik kinetyczny substratu TGP, nie
  parametr obcego frameworku. NGFP/ERG (op-uv-as-ngfp) = opcjonalna L2 reduction, nie primary.
- **Q2 (observable units?):** α_eff bezwymiarowy (wykładnik); output_type strukturalny (max claim C).
- **Q3 (instrument?):** lattice MC substratu + blokowe uśrednianie; estymator structure-factor (CG34).
- **Q4 (circularity guard?):** ⚠ KRYTYCZNY — α_eff⁽⁰⁾ MUSI być liczone bez α_s i bez mas kwarków.
  T-anti-circ obowiązkowy w Phase 1 (zero PDG w wejściu; zakaz „odwrotnego" dopasowania do α=2).
- **Q5 (pre-registration?):** reguła + progi w Phase 0 LOCK, immutable timestamp 2026-06-26.
- **Q6 (forced-zero declare-not-verify?):** N/A (brak α_i/ζ_i w zakresie — to cykl substratowy CG).
- **Q7 (anti-Lakatos two-sided?):** ✅ reguła ma wynik POZYTYWNY (A), NEGATYWNY (B), WARUNKOWY (C);
  szczególnie: (B) jest **realnie osiągalny** (CG34 s−1=0) i NIE jest faworyzowany ani ukrywany.
- **Q8 (ASK-RULE NATIVE_PATTERNS §1)?:** Trigger B (predecessor tech-debt) **rozważony** — substrat
  `-J(φ_iφ_j)²` jest jawnie oznaczony przez CG34 jako **patologiczny dla pure-MC**; Phase 1 NIE
  dziedziczy go bezkrytycznie, lecz testuje alternatywny niepatologiczny model substratu (lub
  czysto analityczną homogenizację) — patrz `Phase0_LOCK.md §5`.

### §0.4 — Pre-flight methodology read confirmation (mandatory, KICKOFF §2.6)

- [x] Przeczytano `PPN_AS_PROJECTION.md` §0-§2 (+ §3 three-layer) — native obserwable PRIMARY;
      framework reduction (ERG/NGFP) = consistency map, ostatni etap; tu output strukturalny (nie PPN).
- [x] Przeczytano `TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §0-§1 — ASK-RULE Triggery A-D, tech-debt
      warning §0.3 (substrat-MC most = potencjalny tech-debt; Trigger B armed dla CG34 substratu).
- [x] Przeczytano `M9_RESTRUCTURE_NOTE.md` §1.4 + §3 — M9.1'' jako Path-2 anchor, NIE framework
      (nierelewantne dla cyklu substratowego, ale potwierdzone read).
- [x] Przeczytano `CYCLE_KICKOFF_TEMPLATE.md` §1-§2 — kickoff contract, ordering L1→L2, §2.2 status map.
- [x] Read-lock źródeł (read-only, `Phase0_LOCK.md §1`): `status_map` CG-3/CG-4 + ex200 4/8 + ex202 7/8;
      `op-CG34 Phase_FINAL_close` (niespójność A3, patologia substratu); `op-Kgeo Phase_FINAL` (C);
      skrypty `ex200_cg3_continuum_verification.py` + `ex202_quark_m0_derivation.py` (inspekcja kodu).

**Sign-off:** Claudian @ 2026-06-26 (Phase 0 LOCK). Pełny audit trail → `Phase0_LOCK.md`.

## §1 — Phase 0: scope mapping (✅ WYKONANE @ LOCK — patrz Phase0_LOCK.md)

1. **Read-lock źródeł** (read-only): `status_map` CG-3/CG-4/ex200/ex202; `op-CG34` werdykt + algebra
   `α_eff=s−1`; `op-Kgeo` werdykt (C); inspekcja kodu `ex200`/`ex202`.
2. **Audyt rekoncyliacyjny** (Phase0_LOCK §2): `ex200.run_tests()` = single-L (L=16, b∈{2,4}),
   T3 tol=1,5, r²>0,3 cut, estymator log-log gęstości gradientowej (CG34: artefakt-prone) ⟹ jego
   „α_eff≈2" jest słaby i sprzeczny z CG34 `s−1=0`. `ex202` baseline 7/8 (T6 FAIL: σ_TGP ~712×).
3. **Plombowanie reguły** (A/B/C + progi 0,3/1,0/0,7) + immutable timestamp.
4. **Balance gate:** budżet nowych stałych = 0 (α_eff MIERZONE, nie nowa stała).

## §2 — Phase 1: native derivation (α_eff⁽⁰⁾ z FSS) — USER-GATED (plan: Phase0_LOCK §5)

- **Estymator (poprawiony):** structure-factor / propagator (CG34 Phase 1), NIE log-log gęstości
  gradientowej ex200 (źródło artefaktu c*→0 i przetolerowanego α_eff).
- **FSS:** ≥3 rozmiary L (np. {16,24,32}, opcjonalnie 48 jeśli zasoby) × skan b; ekstrapolacja
  α_eff(L,b)→α_eff⁽⁰⁾; pomiar R²_FSS.
- **Substrat:** model niepatologiczny (CG34: `-J(φ_iφ_j)²` runaway/frozen) — Phase 1 §5 wybiera
  wariant ze scale-separated window LUB czysto analityczną homogenizację `Z(φ)`.
- **T-anti-circ:** zero α_s/mas kwarków; zakaz „odwrotnego" wymuszania α=2; ≥1 niezależny solver.

## §3 — Phase 2-N: native obserwable computation

- Wyliczyć α_eff⁽⁰⁾ + R²_FSS; zastosować plombę §3 Phase0_LOCK (value-blind).
- Rozstrzygnąć niespójność A3: czy emergentne `Z∝φ⁶ (s=3)` [⇒ A], czy `Z∝φ² (s=1)` [⇒ B], czy
  nieoznaczalne [⇒ C].

## §FINAL — Optional L2 framework reduction

[OPTIONAL] Mapowanie α_eff na ERG/FRG fixed-point (K_IR/K_UV ≈ 1,13) lub NGFP (op-uv-as-ngfp η_N*=−2).
Fail NIE unieważnia L1.

## §FINAL+1 — L3 falsification map + propagacja (user-gated, per A/B/C)

- (A): dodatekQ CG-4 [PARTIAL]→[ZAMKNIĘTY], thm:alpha2 status, ledger #42; (B): dodatekQ2 A3 reframe
  + ścieżka substratowa CLOSED-NEGATIVE; (C): ratyfikacja status quo. Wpis STATE.

## §estimate — Skala

Faza A: 1 sesja (ten deliverable). Phase 1 FSS: 1-2 sesje (large-L MC kosztowne; pure-Python ex200
za wolne — wymaga wektoryzacji/numba lub mniejszego N_sweep). Wynik (C) a priori najbardziej
prawdopodobny (CG34 patologia substratu + ex200 4/8); wartość = **przekształcenie sprzecznych
sygnałów w jeden rozstrzygnięty werdykt** (analog #37/#39/#48), ewentualnie zaskakujące (A)/(B).

## Cross-references
- [[../op-CG34-continuum-closure-2026-06-14/Phase_FINAL_close.md]] (#31 — CG-3 ZAMKNIĘTY NUM; niespójność A3)
- [[../op-Kgeo-from-D-uniqueness-2026-06-26/Phase_FINAL_close.md]] (#47-48 — C POSTULATE-CONFIRMED)
- [[../op-A-derivation-from-CG-2026-06-25/Phase_FINAL_close.md]] (#43 — most 𝒜 POSTULATE-CONDITIONAL)
- [[../../core/_meta_latex/status_map.tex]] CG-3/CG-4 [SZKIC], ex200 4/8, ex202 7/8
- [[../op-uv-as-ngfp/]] · [[../op-Csigma-coarse-graining-2026-06-14/]] (powiązane CG/UV)
- skrypty: `tooling/scripts/ex200_cg3_continuum_verification.py`, `tooling/scripts/ex202_quark_m0_derivation.py`
