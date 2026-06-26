---
title: "op-L08-quark-g0-tail-vs-core-audit — audyt epistemiczny: czy zakres sek08b:529 g₀∈[0,817;0,891] to RDZENIOWE sprzężenie g₀ (jak założył sufit HALT-B 2,68×), czy WARTOŚĆ OGONOWA g_min<1 (mierzona z zewnątrz)? Path α z REALITY_CONTACT — rozstrzygnięcie czy HALT-B kwarkowy to artefakt normalizacji, czy realna niewydolność strukturalna."
date: 2026-06-25
type: research-cycle
folder_status: parking
parent: "[[../../audyt/L08_kink_fermion_closure/README.md]]"

# ============== KICKOFF CONTRACT (BINDING post-2026-05-10) ==============
contract:
  # --- L1: Native (MANDATORY) ---
  L1_native:
    output_observable: "Klasyfikacja epistemiczna (kategoria) wielkości g w zakresie sek08b:529 [0,817; 0,891] — pytanie binarne+: czy ta liczba to (A) RDZENIOWE sprzężenie g₀ definiujące soliton (domena, na której zbudowano sufit strukturalny T11=2,68× w cyklu op-L08-Phase6-quark-sector-mass-formula HALT-B), czy (B) WARTOŚĆ OGONOWA g_min<1 z oscylacyjnego ogona (reżim przypróżniowy, mierzony z zewnątrz, ≠ wartość w układzie własnym solitonu)? Wielkość pochodna: czy argument sufitu 2,68× pozostaje ważny po rozstrzygnięciu kategorii."
    measurement_instrument: "Tekst rdzenia LIVE (read-only): sek08b_ghost_resolution.tex (g₀ rdzeniowe leptonów {e≈1; μ≈2,0; τ≈2,34} l.346; ogon g_min<1 l.60/328-329/347-349; claim universalności l.528-529) + dodatek app:ogon-masy (definicja A_tail(g₀), hip:J-mass-Atail4, eq:J-tail, eq:J-ode) + why_n3 PHASE2 (uniwersalny wzór masowy + kalibracja leptonowa) + cykl-poprzednik op-L08-Phase6-quark-sector-mass-formula (faktyczne użycie g₀ w T11 ceiling)"
    native_coefs_constrained:
      - "Definicja g₀ (rdzeniowe sprzężenie wejściowe ODE solitonu) vs g(r) (profil) vs g_min (ekstremum ogona oscylacyjnego) — rozróżnienie z eq:J-ode / eq:J-tail"
      - "Domena g użyta w T11 ceiling poprzednika: max(m)/min(m) ≤ (A_max/A_min)²·(1+ε)^(e²/2) z ε wyliczonym z [0,817;0,891]"
      - "Empiryczny rozstęp rdzeniowego g₀ leptonów: {1,0; 2,0; 2,34} ⇒ rdzeniowy zakres szer. ≈ 2,34× (vs ogonowy [0,817;0,891] szer. ≈ 1,09×)"
    falsification_rule: "Patrz §0.2 — reguła trójdzielna NORM-OVERLOAD / NORM-COHERENT / NORM-INDETERMINATE, LOCKED PRZED jakimkolwiek rachunkiem ratunkowym."
    pre_registration_date: "2026-06-25"

  # --- L2: Cross-framework reduction (OPTIONAL) ---
  L2_framework_reduction:
    target_frameworks:
      - "Spójność wewnętrzna rdzenia TGP (sek08b ↔ dodatek app:ogon-masy ↔ why_n3 ↔ cykl HALT-B)"
    reduction_type: "internal-consistency-audit (NIE cross-framework; audyt notacyjno-pojęciowy)"
    failure_disposition: "L1-stands"

  # --- L3: Falsification map (consistency) ---
  L3_falsification_map:
    - { bound: "g₀ rdzeniowe leptonów {1,0; 2,0; 2,34} (sek08b:346)", constrains: "Czy [0,817;0,891] ⊂ zakres rdzeniowy {1..2,34}? — NIE (cały < 1)", window: "porównanie zbiorów", status: "pending Phase 1" }
    - { bound: "g_min < 1 ogona oscylacyjnego (sek08b:60,328-329)", constrains: "Czy [0,817;0,891] ⊂ region ogonowy g<1? — TAK prima facie", window: "porównanie zbiorów", status: "pending Phase 1" }
    - { bound: "Sufit T11 = 2,68× (cykl HALT-B, ε≈0,091 z [0,817;0,891])", constrains: "Ważność sufitu zależy od kategorii g; jeśli g=ogonowe ⇒ ε rdzeniowe ≫ 0,091 ⇒ sufit VOIDED", window: "logiczna ważność argumentu", status: "pending Phase 1" }
    - { bound: "PR-004/PR-010/PR-014-candidate + werdykt HALT-B IMMUTABLE", constrains: "Werdykt HALT-B cyklu-poprzednika NIETYKALNY (ten audyt NIE modyfikuje go retroaktywnie; bada jego ZAKRES WAŻNOŚCI, nie jego poprawność dla testowanej hipotezy)", window: "anti-Lakatos invariant", status: "LOCKED" }

# ============== END KICKOFF CONTRACT ==============

tgp_status:
  level: L1
  kind: audit
  output_type: epistemic-status (kategoria wielkości + ważność argumentu strukturalnego)
  core_compatibility: review-only
  may_edit_core: false
  has_needs_file: false
  has_findings_file: false
  exports_findings: true
  open_bridges:
    - "L08 (kink fermion closure) — problem #3 quark sub-component: status HALT-B re-examined"
    - "sek08b:529 universalność kwarkowa — kandydat korekty notacyjnej (g₀ przeciążone?)"
  depends_on:
    - "op-L08-Phase6-quark-sector-mass-formula-2026-05-16 (HALT-B; T11 ceiling 2,68×; R3 flag: 'audit range may correspond to different normalization')"
    - "op-FFS-quark-object-2026-05-20 (A− cond.: kwarki = obiekty FFS z ogonem/strukturą zewnętrzną; Φ_0_local relacyjne)"
    - "why_n3 Phase 1-5 (uniwersalny wzór m = c_M·A_tail²·g₀^(e²/2) + kalibracja leptonowa)"
    - "core/sek08b_ghost_resolution: g₀ rdzeniowe {1;2,0;2,34}; ogon g_min<1; claim l.529 g₀∈[0,817;0,891]"
    - "dodatek app:ogon-masy: A_tail(g₀), eq:J-ode, eq:J-tail, hip:J-mass-Atail4"
  impacts:
    - "audyt/L08_kink_fermion_closure problem #3 (quarks) — status HALT-B: utrzymany LUB zdjęty do reopened"
    - "core/sek08b_ghost_resolution l.528-529 — korekta notacyjna (jeśli NORM-OVERLOAD)"
    - "core/sek07_predykcje l.296 R12 OPEN (m_b/m_t recovery) — kontekst"
    - "PREDICTIONS_REGISTRY — PR-014 formalizacja (jeśli NORM-COHERENT) LUB licencja nowego cyklu rescue-test (jeśli NORM-OVERLOAD)"
    - "README/headline '40 predykcji z 3 inputów' — uczciwość liczenia sektora kwarkowego"
  source_of_status:
    - "STATE.md sesja #39 WIP: 'Następna najważniejsza rzecz: quark-mass HALT-B (sufit 2,68× vs 80000×)'"
    - "op-L08-Phase6-quark-sector-mass-formula §9.1 Path α (audit refinement) — REKOMENDOWANY next"

predecessors:
  - "[[../op-L08-Phase6-quark-sector-mass-formula-2026-05-16/]] (HALT-B; T11 ceiling 2,68× vs req. 80000×; R3 flag normalization)"
  - "[[../op-FFS-quark-object-2026-05-20/]] (A− cond.; kwark = FFS z ogonem; ładunki ułamkowe z windingu)"
  - "[[../why_n3/]] Phase 1-5 (uniwersalny wzór masowy + kalibracja leptonowa)"

related:
  - "[[../../audyt/L08_kink_fermion_closure/README.md]] (problem #3 quark sub-component)"
  - "[[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§3 (BINDING contract, L1 observable-first)"
  - "[[../../meta/CALIBRATION_PROTOCOL.md]] §3.6 (analytical pre-derivation + 8/8 gate)"
  - "[[../../meta/PRE_REGISTERED_FALSIFIERS.md]] (PR-014 candidate; formalizacja contingent na werdykt)"
  - "[[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]] l.346 (g₀ rdzeniowe), l.60/328-329/347-349 (ogon g_min<1), l.528-529 (claim universalności)"
  - "[[../../core/sek07_predykcje/sek07_predykcje.tex]] l.296 (R12 OPEN m_b/m_t)"

classification: AUDIT — epistemiczny status zakresu g w sek08b:529; rozstrzygnięcie kategorii (rdzeń vs ogon) determinujące ważność sufitu HALT-B
priority: high (REALITY_CONTACT: jedyny ciężki problem OTWARTY i badalny czystą analizą teraz; galaktyki wyczerpane, α=2/c₀ domknięte aksjomatycznie, GW strukturalne data-gated)
goal: "Value-blind audyt: ustalić Z TEKSTU RDZENIA + dodatku ogon-masy + cyklu-poprzednika, czy liczba g w sek08b:529 [0,817; 0,891] oznacza (A) RDZENIOWE sprzężenie g₀ (domena sufitu T11=2,68×) czy (B) WARTOŚĆ OGONOWĄ g_min<1 (reżim przypróżniowy mierzony z zewnątrz, ≠ rest-frame solitonu). Konsekwencja determinowana mechanicznie: jeśli (B), domena sufitu jest błędna kategorią ⇒ argument 2,68× VOIDED ⇒ HALT-B reopened + licencja osobnego cyklu rescue-test; jeśli (A), sufit stoi ⇒ HALT-B CONFIRMED-STRENGTHENED + PR-014 formalizacja + uczciwa demotacja predykcji. Werdykt trójdzielny: NORM-OVERLOAD / NORM-COHERENT / NORM-INDETERMINATE."
estimated_effort: "~1 sesja (Phase 0 + Phase 1 audit/sympy + Phase FINAL)"
target_window: "Phase 1: ekstrakcja definicji (eq:J-ode, eq:J-tail, A_tail(g₀)) + przyporządkowanie kategorii [0,817;0,891] + test ważności sufitu T11 pod obiema kategoriami; honest verdict per faktyczne definicje, NIE per preferencja ratunku."

six_requirements_target:
  - "P1: Ekstrakcja definicji g₀ vs g(r) vs g_min z eq:J-ode + eq:J-tail (co dokładnie jest wejściem ODE, a co ekstremum profilu) — FP"
  - "P2: Lokalizacja [0,817;0,891] względem (a) zbioru rdzeniowych g₀ leptonów {1;2,0;2,34}, (b) regionu ogonowego g_min<1 — FP, rozstrzygnięcie zbiorowo-przedziałowe"
  - "P3: Rekonstrukcja, którą wielkość T11 ceiling poprzednika wstawił do max(m)/min(m) ≤ (A_max/A_min)²·(1+ε)^(e²/2) — FP, audyt logiczny argumentu"
  - "P4 (CENTRALNY): Werdykt kategorii — NORM-OVERLOAD vs NORM-COHERENT vs NORM-INDETERMINATE per reguła §0.2, WYLICZONY z P1-P3, NIE wybrany"
  - "P5: Test ważności sufitu pod każdą kategorią: jeśli g=rdzeniowe ⇒ jaki realny ε rdzeniowy (z {1;2,0;2,34} szer. 2,34×)? jaki sufit wtedy? — FP, liczba bez fitu do PDG"
  - "P6: S05 single-Φ + werdykt HALT-B poprzednika NIETYKALNE (DEC; ten audyt nie modyfikuje retroaktywnie żadnego LOCKa)"

risk_flags:
  - "R1 (META/HIGH): pokusa ratunku — 'tail insight' usera czyni NORM-OVERLOAD atrakcyjnym; antidotum: werdykt WYLICZONY z definicji tekstu rdzenia PRZED jakąkolwiek rekalkulacją stosunków mas; rescue-test = OSOBNY cykl/PR, ZAKAZANY w tej sesji"
  - "R2: notacja g₀ w rdzeniu może być faktycznie przeciążona (to samo '$g_0$' dla wejścia ODE i dla wartości profilu) — wtedy werdykt NORM-OVERLOAD jest UCZCIWY, nie naciągany; ale wymaga jawnego cytatu eq:J-ode pokazującego co jest argumentem"
  - "R3: app:ogon-masy może definiować A_tail i g₀ tak, że [0,817;0,891] to JESZCZE INNA wielkość (np. g* próg duchowy, albo g_min konkretnych leptonów) — wtedy NORM-INDETERMINATE lub re-scope; zbiór kategorii musi pokryć tę możliwość"
  - "R4: nawet jeśli NORM-OVERLOAD (sufit VOIDED), to NIE jest dowód, że kwarki SĄ reprodukowalne — to tylko zdjęcie błędnego no-go; rescue-test (osobny PR) może dać HALT-B z INNEGO powodu. Zakaz nadinterpretacji 'sufit padł ⇒ kwarki działają'."
  - "R5: BD-drift — brak; pracujemy z klasyczną masą solitonu i definicjami profilu, bez propagatora Φ"

phase_plan:
  Phase_0: "Balance sheet + 8/8 ☑ gate + zbiór kategorii CLOSED + reguła werdyktu LOCKED + forbidden moves + pre-derywacja (co tekst MÓWI prima facie, przed pełną ekstrakcją)"
  Phase_1: "Sympy/audyt: P1 ekstrakcja definicji (eq:J-ode/eq:J-tail) → P2 przyporządkowanie [0,817;0,891] → P3 rekonstrukcja domeny T11 → P4 werdykt kategorii (wyliczony) → P5 sufit pod każdą kategorią; circularity guard: brak PDG quark mass w fazie kategoryzacji"
  Phase_FINAL: "Honest verdict NORM-*; dyspozycja HALT-B (utrzymany/reopened); korekta sek08b:529 (jeśli OVERLOAD); PR-014 formalizacja (jeśli COHERENT) lub licencja rescue-test cycle (jeśli OVERLOAD); propagacja"

tags:
  - L08
  - quark-sector
  - normalization-audit
  - tail-vs-core
  - g0-overload
  - HALT-B-reexamination
  - structural-ceiling-validity
  - value-blind
  - anti-Lakatos-LOCKED
  - decision-trichotomous
  - cycle-scaffold-2026-06-25
---

# op-L08-quark-g0-tail-vs-core-audit-2026-06-25

> **Cel (Path α z REALITY_CONTACT + rekomendacja §9.1 cyklu HALT-B):** value-blind
> audyt epistemiczny rozstrzygający, **czym jest liczba** `g ∈ [0,817; 0,891]` w
> sek08b:529 ("ten sam ODE działa na leptony i~kwarki"): **rdzeniowym sprzężeniem
> g₀** (na którym zbudowano sufit strukturalny `T11 = 2,68×` cyklu HALT-B) — czy
> **wartością ogonową `g_min < 1`** z oscylacyjnego ogona solitonu (reżim
> przypróżniowy mierzony z zewnątrz, **różny od wartości w układzie własnym kwarka**
> — wprost uwaga inicjująca użytkownika 2026-06-25). **Werdykt trójdzielny:
> NORM-OVERLOAD / NORM-COHERENT / NORM-INDETERMINATE.**

## §0 — Cel + native-first contract

### §0.0 — Tożsamość cyklu (nienegocjowalne)

1. **TO JEST AUDYT, NIE RATUNEK.** Werdykt HALT-B cyklu
   `op-L08-Phase6-quark-sector-mass-formula-2026-05-16` jest **IMMUTABLE i NIETYKALNY**.
   Ten cykl **nie podważa**, że uniwersalny wzór z `g₀ ∈ [0,817; 0,891]` daje 0/5
   stosunków — to zostało policzone poprawnie **dla testowanej tam hipotezy**. Pytanie
   tego cyklu jest węższe i czysto epistemiczne: **czy hipoteza testowana w HALT-B
   (że domeną g₀ kwarków jest [0,817; 0,891]) była właściwą kategorią fizyczną** — czy
   też cykl HALT-B przetestował **wartość ogonową g_min** pod nazwą "g₀ rdzeniowe",
   przez przeciążenie notacji w sek08b:529.
2. **KAŻDY z trzech werdyktów jest pełnoprawnym sukcesem metodologicznym.**
   NORM-COHERENT (sufit stoi, HALT-B potwierdzony, predykcja zdemotowana uczciwie)
   jest **równie wartościowy** jak NORM-OVERLOAD. Wynik NIE jest z góry znany; uwaga
   usera o ogonie czyni OVERLOAD *hipotezą do sprawdzenia*, NIE konkluzją.
3. **Rescue-test = osobny cykl + osobny PR.** Nawet przy NORM-OVERLOAD (sufit VOIDED),
   faktyczna rekalkulacja 5 stosunków mas kwarków z rdzeniowym g₀ jest **POZA zakresem
   tej sesji** (forbidden move #1). Ten audyt może najwyżej *zlicencjonować* taki cykl.
4. Każda faza = osobne „działaj"; po każdej fazie raport + wpis STATE.md.

### §0.1 — Native observable target

**Co rozstrzygamy (epistemiczna obserwabla, NIE liczbowa predykcja):**

Kategoria wielkości `g` w przedziale **[0,817; 0,891]** (sek08b:529), wybór ze zbioru
CLOSED:

- **(A) RDZENIOWE g₀** — bezpośredni parametr wejściowy ODE solitonu (eq:J-ode);
  wielkość, która dla leptonów przyjmuje `{e: ≈1,0; μ: ≈2,0; τ: ≈2,34}` (sek08b:346).
- **(B) WARTOŚĆ OGONOWA g_min** — ekstremum oscylacyjnego ogona profilu `g(r)`
  (eq:J-tail), `g_min < 1`, reżim przypróżniowy (sek08b:60, 328-329, 347-349);
  to jest wielkość **mierzona z zewnątrz**, różna od wartości rdzeniowej (rest-frame).
- **(C) INNA wielkość** — np. próg duchowy `g*`, albo zupełnie odrębna normalizacja
  (np. parametr podgrupy chiralnej) — wymaga jawnego cytatu z app:ogon-masy.

**Instrument:** wyłącznie tekst LIVE rdzenia + dodatku + cykli (read-only). **Żadnych
danych PDG kwarków w fazie kategoryzacji** (circularity guard: kategoria g nie może
być wybrana przez to, która kategoria "ratuje" stosunki mas).

### §0.2 — Pre-registered falsification rule (LOCKED 2026-06-25, PRZED rachunkiem)

> **Reguła werdyktu (trójdzielna; kryteria wyliczalne z tekstu, nie z preferencji):**
>
> Niech `D_core = {1,0; 2,0; 2,34}` (rdzeniowe g₀ leptonów, sek08b:346) i niech
> `R_tail = {g : g_min < 1}` (region ogonowy, sek08b:60/328-329). Niech `I = [0,817; 0,891]`.
>
> - **NORM-OVERLOAD** ⟺ ekstrakcja eq:J-ode/eq:J-tail wykazuje, że `I` indeksuje
>   **wartości profilu/ogona (kategoria B)**, a NIE parametr wejściowy ODE (kategoria A),
>   **ORAZ** sek08b:529 używa symbolu „$g_0$" dla tej wartości ogonowej (przeciążenie
>   notacji udokumentowane jawnym cytatem). Warunek wspierający (NIE wystarczający sam):
>   `I ⊄ conv(D_core)` (cały przedział < 1, leży poniżej najmniejszego rdzeniowego g₀=1).
>   **Konsekwencja LOCKED:** domena sufitu T11 (ε≈0,091 z `I`) jest **błędną kategorią**
>   ⇒ argument 2,68× **VOIDED** (nie sfalsyfikowany — *nieważny dla rdzeniowego g₀*) ⇒
>   HALT-B **reopened** (status: STRUCTURAL_INSUFFICIENCY → INDETERMINATE-PENDING-RESCUE);
>   licencja na osobny cykl `op-quark-mass-core-g0-rescue-test` (nowy PR, własny Phase 0).
>
> - **NORM-COHERENT** ⟺ ekstrakcja wykazuje, że `I` indeksuje **rdzeniowe g₀ (kategoria A)**
>   — tj. sek08b:529 twierdzi (poprawnie wg definicji ODE), że kwarki mają rdzeniowe
>   sprzężenie w wąskim oknie [0,817; 0,891], **różnym** od leptonowego {1..2,34}.
>   **Konsekwencja LOCKED:** sufit T11=2,68× jest **ważny** ⇒ HALT-B **CONFIRMED-STRENGTHENED**
>   ⇒ PR-014 formalizacja w rejestrze + uczciwa demotacja: nagłówek schodzi do
>   „lepton sektor derived; quark sektor (H) — universalność OBALONA empirycznie".
>
> - **NORM-INDETERMINATE** ⟺ tekst rdzenia + app:ogon-masy + why_n3 **nie wystarczają**
>   do jednoznacznego przyporządkowania (kategoria C bez rozstrzygnięcia, lub sprzeczne
>   sygnały między dokumentami). **Konsekwencja LOCKED:** GAP udokumentowany; HALT-B status
>   QUO (utrzymany jako testujący hipotezę A, z jawną adnotacją „domena niezweryfikowana");
>   rekomendacja: minimalne doprecyzowanie definicji w app:ogon-masy (osobny housekeeping).

```yaml
pre_registration_date: 2026-06-25
pre_registration_hash: <auto-set by git commit SHA>
recovery_scope:
  allowed_directions:
    - "Ekstrakcja definicji wyłącznie z LIVE tekstu (eq:J-ode, eq:J-tail, app:ogon-masy, why_n3)"
    - "Audyt logiczny argumentu T11 (która wielkość weszła do ε) — bez modyfikacji werdyktu HALT-B"
    - "Obliczenie hipotetycznego sufitu pod kategorią A vs B (liczba strukturalna, BEZ fitu do PDG)"
  forbidden_directions:
    - "Rekalkulacja 5 stosunków mas kwarków z rdzeniowym g₀ (= rescue-test; OSOBNY cykl/PR)"
    - "Modyfikacja/retraktacja werdyktu HALT-B poprzednika (IMMUTABLE)"
    - "Wybór kategorii g motywowany tym, która 'ratuje' stosunki mas (value-blind violation)"
    - "Post-hoc redefinicja przedziału [0,817;0,891] lub progu"
    - "Wprowadzenie nowego pola/stałej (S05; budżet nowych stałych 0)"
  if_recovery_exhausted:
    - "NORM-INDETERMINATE: GAP udokumentowany; housekeeping doprecyzowania definicji w app:ogon-masy"
```

### §0.3 — TGP-native check (mandatory, pre-Phase-1)

- [x] **Q1 (Pattern coverage):** Pattern 2.7 (asymptotic matching ogona A_tail); audyt definicyjny, nie nowy rachunek pola.
- [x] **Q2 (Red flags):** NONE — brak BD-form, brak m_Φ, brak Φ-quantum carrier; klasyczny profil solitonu.
- [x] **Q3 (Inherited LOCKs):** HALT-B IMMUTABLE; T11 ceiling = przedmiot audytu logicznego (NIE modyfikacji); kalibracja leptonowa why_n3 LIVE; FFS quark=obiekt z ogonem (A− cond.) LIVE.
- [x] **Q4 (Standard-physics tools):** ekstrakcja definicji + sympy do rekonstrukcji ε/sufitu pod obiema kategoriami; universal narzędzia.
- [x] **Q5 (m_Φ usage):** N/A — profil solitonu g(r), bez propagatora Φ.
- [x] **Q6 (GR limit):** N/A.
- [x] **Q7 (ASK-RULE Trigger A — form-meaning split):** „g₀" w TGP ma **dwa możliwe znaczenia** (wejście ODE rdzenia vs wartość profilu/ogona); rozróżnienie ICH jest CAŁYM przedmiotem cyklu — declarative, jawnie deklarowane PRZED rachunkiem.
- [x] **Q8 (BD-drift audit plan):** self-audit w Phase FINAL (niskie ryzyko; brak propagatora).

### §0.4 — Pre-flight methodology read confirmation (per KICKOFF §2.6)

- [x] Przeczytano [[../../STATE.md]] sesje #36-#39 (kontekst: α=2/c₀ domknięte aksjomatycznie; #39 WIP nominuje quark-mass HALT-B jako następną najważniejszą rzecz)
- [x] Przeczytano [[../op-L08-Phase6-quark-sector-mass-formula-2026-05-16/README.md]] + Phase_FINAL_close (BINDING contract + T11 ceiling 2,68× + R3 normalization flag + §9.1 Path α rekomendacja)
- [x] Przeczytano [[../op-FFS-quark-object-2026-05-20/Phase_FINAL_close.md]] (kwark = obiekt FFS z ogonem/strukturą zewnętrzną; Φ_0_local relacyjne — wsparcie ontologiczne dla rozróżnienia rdzeń/ogon)
- [x] Przeczytano [[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]] l.34-62 (prob:kinetic-singularity), l.326-374 (rem:heavy-solitons + prop:Atail-preserved: g₀ rdzeniowe {1;2,0;2,34}, ogon g_min<1, A_tail invariant), l.526-531 (claim universalności kwarkowej g₀∈[0,817;0,891])
- [x] Przeczytano [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] §0-§2 (format; PR-014 candidate; anti-moving-goalposts)
- [ ] **DO ZROBIENIA w Phase 1 (część merytoryczna):** pełna ekstrakcja `dodatek app:ogon-masy` (eq:J-ode, eq:J-tail, hip:J-mass-Atail4, definicja A_tail(g₀)) + `why_n3 PHASE2` (uniwersalny wzór + kalibracja) — to jest RACHUNEK Phase 1, nie pre-flight.

**Sign-off:** Claudian (theoretical physics agent) @ 2026-06-25.

### §0.5 — Sympy/audyt substance plan (Phase 1)

| Test | Klasa | Pytanie | Pre-decyzja PASS |
|---|---|---|---|
| T1 | **FIRST_PRINCIPLES** | Co jest argumentem wejściowym eq:J-ode (rdzeniowe g₀) vs co jest wartością profilu g(r)/g_min? | Jawne rozróżnienie z równania |
| T2 | **FIRST_PRINCIPLES** | Czy [0,817;0,891] ⊂ conv(D_core={1;2,0;2,34})? | Test przynależności: NIE (cały <1) |
| T3 | **FIRST_PRINCIPLES** | Czy [0,817;0,891] ⊂ R_tail (g_min<1 region z eq:J-tail)? | Test przynależności |
| T4 | **LITERATURE_ANCHORED** | Jaki symbol/słowo sek08b:529 wiąże z [0,817;0,891] — „$g_0$"? cytat verbatim | Cytat + klasyfikacja notacji |
| T5 | **FIRST_PRINCIPLES** | Którą wielkość T11 ceiling poprzednika wstawił do ε=(g_max−g_min)/g_ref? rekonstrukcja | Audyt logiczny: ε≈0,091 z `I` |
| T6 | **FIRST_PRINCIPLES** | Werdykt kategorii (A/B/C) WYLICZONY z T1-T5 per reguła §0.2 | NORM-* determinowany, nie wybrany |
| T7 | **FIRST_PRINCIPLES** | Sufit hipotetyczny pod kategorią A (ε z `I`) vs pod kategorią B (ε rdzeniowe z {1..2,34}, szer. 2,34×) | Dwie liczby strukturalne, bez fitu PDG |
| T8 | **FIRST_PRINCIPLES** | Circularity guard: czy werdykt T6 użył JAKIEJKOLWIEK danej PDG kwarka? | free-symbols audit: NIE |
| T9 | **DECLARATIVE** | S05 single-Φ + HALT-B IMMUTABLE preserved | separate from PASS count |

**Substance:** 7 FP + 1 LIT + 1 DEC (target ≥75% FP). **0 hardcoded T_pass.** Werdykt
empiryczny per faktyczne definicje.

---

## §1 — Phase 0: balance sheet + 8/8 gate
[Patrz `Phase0_balance.md` — LOCKED]

## §2 — Phase 1: audyt kategoryzacji
[Patrz `Phase1_audit.md` + `Phase1_sympy.py/.txt` — wymaga „działaj"]

## §FINAL — Honest verdict NORM-* + dyspozycja HALT-B
[Patrz `Phase_FINAL_close.md`]

---

## Status

🟡 **PARKING — scaffold + Phase 0 LOCK 2026-06-25**. Pre-flight §0.4 complete (pełna
ekstrakcja app:ogon-masy + why_n3 = zadanie Phase 1). Autoryzacja: user kickoff
2026-06-25 („działaj z tym, rozpisz nowy cykl" + uwaga o ogonie/rest-frame).

**Phase 1 FAST-AUDIT wymaga osobnego „działaj".**

**Cross-references:**
- [[../op-L08-Phase6-quark-sector-mass-formula-2026-05-16/Phase_FINAL_close.md]] (HALT-B, T11, §9.1 Path α)
- [[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]] (l.346 rdzeń, l.60/328-329 ogon, l.528-529 claim)
- [[../op-FFS-quark-object-2026-05-20/Phase_FINAL_close.md]] (kwark=FFS z ogonem)
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] (PR-014 candidate)
