---
title: "op-parameter-counting-balance-sheet — program-wide uczciwy bilans inputów: czy headline '40 predykcji z 3 inputów' jest UCZCIWY, OPTYMISTYCZNY, czy MYLĄCY? Inwentarz wszystkich wolnych parametrów (ciągłych) + aksjomatów selekcji (dyskretnych) vs predykcje; porównanie do SM (19+). Synteza #36-#41 + manuskrypt."
date: 2026-06-25
type: research-cycle
folder_status: parking
parent: "[[../../README.md]]"

# ============== KICKOFF CONTRACT (BINDING) ==============
contract:
  L1_native:
    output_observable: "Uczciwy licznik: N_free (ciągłe wolne parametry numeryczne) + N_axiom (dyskretne aksjomaty selekcji) vs N_pred (genuine predykcje, NIE fity). Werdykt o uczciwości headline'u '40 predykcji z 3 inputów'. Porównanie do SM (19+ wolnych parametrów)."
    measurement_instrument: "Manuskrypt LIVE (README highlights; sek00/sek04; TGP_FOUNDATIONS §3.5.3 Φ_0=EFT free param, §3.5.3.1 g̃≈0,98 postulate; dodatekJ/dodatekX masy) + ustalenia #36-#41 (α=2 #36/#38/#39 aksjomatyczne; c₀ #37 wolny UV; quark 4 inputy #41; Ω_Λ↔g̃ trade)."
    native_coefs_constrained:
      - "N_free: g₀^e, Φ_0 (EFT free, FOUNDATIONS §3.5.3), c₀/C_σ (wolny UV #37), quark anchors (4, #41), g̃/β=γ (postulate)"
      - "N_axiom (dyskretne): N=3, Z₂, α=2 (selekcja klasy C1-C3, #36/#38/#39), φ-FP (golden ladder), Koide-B=√2"
      - "N_pred genuine: r_21/r_31 leptony (φ+Koide), r_31 kwarki (#41), Koide K=2/3, m_W, n_s, r, α_s(z 𝒜 warunkowe)..."
      - "klasyfikacja każdej pozycji: FREE-NUMERICAL / STRUCTURAL-AXIOM / DERIVED / TRADED (input↔postulate)"
    falsification_rule: "§0.2 — HEADLINE-HONEST / HEADLINE-OPTIMISTIC / HEADLINE-MISLEADING, progi LOCKED przed tally."
    pre_registration_date: "2026-06-25"
  L2_framework_reduction:
    target_frameworks: ["Standard Model (19+ wolnych parametrów: 9 mas fermionów + 4 CKM + 3 cechowania + 2 Higgs + θ_QCD; +7 neutrina)"]
    reduction_type: "honest-comparison (NIE redukcja; bilans porównawczy)"
    failure_disposition: "L1-stands"
  L3_falsification_map:
    - { bound: "headline README '40 predykcji z 3 inputów (+structural selection axioms)'", constrains: "uczciwy N_free vs 3", window: "progi §0.2", status: "pending tally" }
    - { bound: "SM 19+ wolnych parametrów", constrains: "czy TGP istotnie ekonomiczniejszy", window: "N_free_TGP ≪ 19?", status: "pending" }
    - { bound: "Ω_Λ 'promoted to prediction' (T-Λ closure)", constrains: "czy genuine win czy input↔g̃ trade (FOUNDATIONS §3.5.3.1)", window: "klasyfikacja TRADED", status: "pending" }
    - { bound: "ustalenia #36-#41 IMMUTABLE", constrains: "α=2 aksjomat, c₀ wolny, quark 4 inputy — wejścia do tally, nie re-litygacja", window: "anti-Lakatos", status: "LOCKED" }

# ============== END KICKOFF CONTRACT ==============

tgp_status:
  level: L1
  kind: meta-audit (balance-sheet)
  output_type: epistemic (uczciwy licznik + werdykt o headline)
  core_compatibility: review-only
  may_edit_core: false
  exports_findings: true
  depends_on:
    - "#36 (α=2 = aksjomat selekcji na gęstości, P4: 'czy 3 inputy uczciwe')"
    - "#37 (c₀ = wolny parametr UV)"
    - "#38/#39 (α=2 nieredukowalnie aksjomatyczne)"
    - "#41 (quark sektor: 4 inputy → 2 predykcje; m_t scheme-zależny)"
    - "TGP_FOUNDATIONS §3.5.3 (Φ_0=EFT free param), §3.5.3.1 (g̃≈0,98 postulate; β=γ open)"
    - "README highlights (40 predykcji, 3 inputy + structural selection axioms — nota #36 P4)"
  impacts:
    - "README headline + Abstract (uczciwy framing inputów)"
    - "tgp_letter/tgp_companion (submission — headline)"
    - "STATE.md (pozycja epistemiczna programu)"

predecessors:
  - "[[../op-quark-mass-core-g0-rescue-test-2026-06-25/]] (#41 quark 4 inputy)"
  - "[[../op-alpha2-status-propagation-audit-2026-06-22/]] (#36 P4: pytanie '3 inputy uczciwe?')"
  - "[[../op-c0-derivation-from-substrate-2026-06-22/]] (#37 c₀ wolny)"
  - "[[../op-M03-balance-sheet-retrofit-2026-05-06/]] (precedens balance-sheet methodologia)"

related:
  - "[[../../TGP_FOUNDATIONS.md]] §3.5.3 (Φ_0 EFT free), §3.5.3.1 (g̃ postulate)"
  - "[[../../audyt/NUMERICAL_ANCHORS_REGISTRY.md]] (anchor vs derivation framework)"
  - "[[../../meta/PRE_REGISTERED_FALSIFIERS.md]]"

classification: META-AUDIT — program-wide uczciwy bilans inputów vs predykcji; werdykt o headline; SM comparison
priority: high (#36 P4 deliverable; bezpośrednio dotyka uczciwości głównego claimu TGP; po domknięciu korzeni α=2/c₀/quark to najwyższa dźwignia epistemiczna)
goal: "Value-blind: zbudować uczciwy ledger WSZYSTKICH wolnych parametrów (ciągłych: g₀^e, Φ_0, c₀, quark anchors, g̃/β=γ) + aksjomatów selekcji (dyskretnych: N=3, Z₂, α=2, φ-FP, Koide-B) vs genuine predykcji. Sklasyfikować Ω_Λ (TRADED czy WIN). Policzyć N_free uczciwie. Porównać do SM (19+). Werdykt: HEADLINE-HONEST / OPTIMISTIC / MISLEADING."
estimated_effort: "~1 sesja (Phase 0 + Phase 1 tally + FINAL)"

risk_flags:
  - "R1 (META/HIGH, DWUSTRONNE): bias może iść w OBIE strony — (a) inflacja N_free (liczenie symetrii Z₂ jak parametru — niesprawiedliwe vs SM, gdzie grupa cechowania nie jest 'parametrem'); (b) deflacja (ukrywanie Φ_0/c₀/g̃). Antidotum: jawne kryteria FREE vs AXIOM §0.2; porównanie z analogiczną konwencją SM."
  - "R2: 'predykcja' vs 'fit' — niektóre 'predykcje' to fity (1 param→1 obserwabla). Liczyć tylko genuine (N_input < N_output w bloku)."
  - "R3: Ω_Λ 'promotion' — pokusa policzyć jako win; faktycznie g̃≈0,98 postulate (FOUNDATIONS §3.5.3.1) ⇒ TRADED. Nie podwójnie liczyć."
  - "R4: zakaz re-litygacji #36-#41 (IMMUTABLE wejścia) i treści manuskryptu (read-only)."
  - "R5: porównanie SM musi być symetryczne (ta sama konwencja co dla TGP: grupa cechowania/symetria NIE liczona jako 'parametr' po obu stronach)."

phase_plan:
  Phase_0: "Balance + 8/8 gate + kryteria FREE/AXIOM/DERIVED/TRADED LOCKED + progi werdyktu + szkielet ledgera"
  Phase_1: "Pełny tally: enumeracja wszystkich pozycji z klasyfikacją + cytatem źródła; N_free, N_axiom, N_pred; SM comparison symetryczny; werdykt wyliczony"
  Phase_FINAL: "Werdykt HEADLINE-*; rekomendacja uczciwego framingu; propagacja (README/submission); PR brak (audyt)"

tags:
  - parameter-counting
  - balance-sheet
  - input-honesty
  - meta-audit
  - SM-comparison
  - free-vs-axiom
  - value-blind
  - anti-Lakatos-LOCKED
  - synthesis-36-41
  - cycle-scaffold-2026-06-25
---

# op-parameter-counting-balance-sheet-2026-06-25

> **Cel (#36 P4 deliverable):** zbudować **program-wide uczciwy bilans inputów** i
> rozstrzygnąć, czy headline TGP **„40 predykcji z 3 inputów"** jest **UCZCIWY**,
> **OPTYMISTYCZNY**, czy **MYLĄCY** — licząc WSZYSTKIE wolne parametry (ciągłe) +
> aksjomaty selekcji (dyskretne) vs genuine predykcje, z symetrycznym porównaniem do
> SM (19+ wolnych parametrów). Synteza ustaleń #36–#41 + manuskrypt.

## §0 — Cel + contract

### §0.0 — Tożsamość cyklu
1. **Meta-audyt syntetyzujący, NIE nowa fizyka.** Wejścia: ustalenia #36-#41 (IMMUTABLE) +
   manuskrypt (read-only). Zakaz re-litygacji (R4).
2. **Bias DWUSTRONNY (R1).** Pokusa NIE jest jednokierunkowa: równie niesprawiedliwe byłoby
   zawyżać N_free (licząc Z₂ jak parametr) jak zaniżać (ukrywając Φ_0/c₀/g̃). Kryteria
   FREE/AXIOM LOCKED §0.2; porównanie SM symetryczne (R5).
3. **Każdy werdykt = sukces.** HEADLINE-OPTIMISTIC (uczciwy licznik > 3, ale TGP nadal ≪ SM)
   jest najprawdopodobniejszy i w pełni wartościowy.
4. Każda faza = osobne „działaj".

### §0.1 — Native observable target
Uczciwy licznik **(N_free, N_axiom, N_pred)** + werdykt o headline + SM comparison.

### §0.2 — Pre-registered classification + verdict rule (LOCKED 2026-06-25)

**Kryteria klasyfikacji (każda pozycja → dokładnie jedna klasa):**
- **FREE-NUMERICAL** = wielkość ciągła strojona do danych bez first-principles derywacji
  (analog masy/sprzężenia SM). Liczona do N_free.
- **STRUCTURAL-AXIOM** = dyskretny wybór (symetria, liczba całkowita, selekcja klasy) —
  NIE strojony ciągle. Liczony do N_axiom (osobno; analog: grupa cechowania SM też NIE
  liczona jako parametr — R5 symetria konwencji).
- **DERIVED** = genuine predykcja: w bloku N_output > N_input (np. leptony 1→3).
- **TRADED** = pozorny win, faktycznie wymiana input↔postulate (np. Ω_Λ↔g̃). NIE liczona
  jako redukcja N_free.

**Reguła werdyktu (progi LOCKED; konwencja symetryczna z SM — symetrie NIE liczone):**
> - **HEADLINE-HONEST** ⟺ uczciwy `N_free ≤ 4` (≈ deklarowane „3", tolerancja +1 na zaokrąglenie/Φ_0).
> - **HEADLINE-OPTIMISTIC** ⟺ `5 ≤ N_free ≤ 12` ORAZ `N_free ≪ N_free^SM` (≈19) ORAZ blok
>   leptonowy genuine DERIVED (1→3). Tj. headline zaniża licznik, ale TGP realnie ekonomiczniejszy.
> - **HEADLINE-MISLEADING** ⟺ `N_free > 12` (zbliża się do SM) LUB większość „predykcji" to
>   fity (N_pred_genuine < N_free) LUB blok leptonowy NIE jest genuine DERIVED.

```yaml
pre_registration_date: 2026-06-25
recovery_scope:
  allowed: ["Enumeracja z cytatem źródła", "Klasyfikacja per kryteria §0.2", "SM comparison symetryczny"]
  forbidden:
    - "Re-litygacja #36-#41 lub manuskryptu (IMMUTABLE/read-only)"
    - "Liczenie Z₂/symetrii jako parametru (asymetria vs SM)"
    - "Liczenie Ω_Λ jako redukcji (TRADED, nie win)"
    - "Post-hoc zmiana progów 4/12"
    - "Inflacja LUB deflacja N_free pod z góry przyjęty werdykt"
```

### §0.3 — TGP-native check
- [x] Q1-Q8: meta-audyt; brak BD-drift; ASK-RULE: „input" rozróżnione FREE vs AXIOM jawnie.

### §0.4 — Pre-flight reads
- [x] #36-#41 cykle (wejścia IMMUTABLE)
- [x] TGP_FOUNDATIONS §3.5.3 (Φ_0 EFT free), §3.5.3.1 (g̃≈0,98 postulate, β=γ open)
- [x] README highlights + #36 P4 nota
- [x] NUMERICAL_ANCHORS_REGISTRY (anchor vs derivation framework)
- [ ] **Phase 1:** sek00/sek04 pełna lista predykcji + PREDICTIONS_REGISTRY status (genuine vs fit)

**Sign-off:** Claudian @ 2026-06-25.

---
## §1 — Phase 0 [Patrz Phase0_balance.md — LOCKED]
## §2 — Phase 1 tally [wymaga „działaj"]
## §FINAL — Werdykt HEADLINE-*

## Status
🟡 **PARKING — scaffold + Phase 0 LOCK 2026-06-25.** Autoryzacja: user „Rewizja parameter-countingu". **Phase 1 wymaga „działaj".**
