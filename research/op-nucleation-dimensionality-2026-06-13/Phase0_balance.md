---
title: "Phase 0 — pre-registration: op-nucleation-dimensionality — value-blind audyt selekcji wymiaru D=3 z maszynerii TGP (topologia defektów + bg-stabilność CE-H + nukleacja + marginalność) oraz przeglądu ND sortowania; falsyfikatory F-ND-A/B/C/D/E LOCKED"
type: phase0_balance
status: LOCKED
locked_date: 2026-06-13
cycle: op-nucleation-dimensionality-2026-06-13
authorization: "User 2026-06-13: 'jesteś ekspertem fizyki teoretycznej; twoje zadanie zająć się cyklem op-nucleation-dimensionality' + 'kontynuuj' → fraza aktywacyjna. Per precedens domu #15/#16/#17 (STATE.md: fraza aktywacyjna pokrywa Phase 0 LOCK). Lektury obowiązkowe §0.4 wykonane PRZED niniejszym LOCK. Phase 1+ = osobne 'działaj' każdorazowo."
methodology_binding: "CYCLE_KICKOFF_TEMPLATE §1–§3 (L1 + §2.6 pre-flight); CYCLE_LIFECYCLE (claim_status taxonomy; output_type: structural ⇒ ceiling C); CALIBRATION_PROTOCOL §3.6 (analytical pre-derivation, sign conventions, assumptions); SCOPING_op-nucleation-dimensionality §2 (kalibracja epistemiczna user — H-SORT working-mechanism BINDING) + §4 (uczciwy test D>3; D_obs comparison-only)."
anti_lakatos_lock: PRESERVED
PR_candidate: "BRAK — cykl strukturalny (output_type: structural). Nowy PR tylko jeśli wyłoni się observable (np. ostro testowalna predykcja związana z D=3); wtedy następny wolny = PR-025, decyzja wyłącznie user w Phase FINAL."
object_under_audit: "core/sek07a_wymiar_wzmocniony/sek07a_wymiar_wzmocniony.tex (prop:wymiar-quantitative) — read-only; jego INTERPRETACJA SELEKCYJNA = hipoteza pod testem, NIE input."
queue_next_after_cycle: "decyzja user (kandydat: op-phi-radiative-dof-audit)"
---

# Phase 0 — pre-registration LOCKED (przed jakimkolwiek rachunkiem)

## §0 — Tożsamość cyklu (nienegocjowalne)

1. **Cykl AUDYTOWY, nie derywacyjny-od-zera.** Rdzeń TGP już zawiera preferencyjny argument
   za D=3 (`sek07a` prop:wymiar-quantitative: homotopia `N_sekt(d)` + potencjał trzech
   reżimów `Δ_d` + wskaźnik `Q(d)`, Q(3)=3 i 0 dla d≠3; konkluzja „d=3 jedynym realistycznym
   wyborem", a zarazem deklaracja „Argument pozostaje **preferencyjny**"). Cytowanie konkluzji
   sek07a jako potwierdzenia selekcji = **ruch zakazany** (założenie odpowiedzi). Zadaniem
   cyklu jest **value-blind rozstrzygnięcie**: czy selekcja D=3 jest DERIVED (z mechanizmu,
   niezależnie od d), czy ARTEFAKTEM konstrukcji Q(d) dopasowanej do D_obs=3.
2. **D_obs = 3 jest comparison-only anchorem, NIGDY inputem.** Każdy rachunek per-D
   (D = 1,2,3,4,5,6) wykonywany jest symetrycznie; porównanie z D_obs=3 wyłącznie PO locku
   wyników per-D, w osobnej sekcji (circularity guard FP w każdym sympy: free-symbols /
   control-flow audit na „3" jako uprzywilejowanym wejściu).
3. **HONEST_NEGATIVE i PARTIAL = pełnoprawne wyniki.** „Maszyneria nie selekcjonuje D
   jednoznacznie" lub „selekcja preferencyjna, nie derived" są poprawnymi werdyktami —
   ceiling claim_status = **C** (output_type: structural; D = liczba całkowita, brak
   observable z jednostkami).
4. Każda faza = osobne „działaj"; po każdej fazie raport + decision menu + wpis STATE.md.

### §0.4 — Pre-flight methodology read confirmation (per KICKOFF §2.6)

- [x] Przeczytano STATE.md (sesje #9–#18; kontekst cykli D/A/B/BA/PR-004/GST)
- [x] Przeczytano CAŁY `core/sek07a_wymiar_wzmocniony/sek07a_wymiar_wzmocniony.tex`
      (OBIEKT AUDYTU — Część I homotopia, Część II Δ_d, wskaźnik Q(d), wniosek)
- [x] Przeczytano `meta/SCOPING_op-nucleation-dimensionality_2026-06-12.md` (osie a–d;
      §2 kalibracja epistemiczna user BINDING; §4 uczciwy test D>3 + D_obs comparison-only)
- [x] Przeczytano cykl wzorcowy `research/op-galactic-substrate-tail-2026-06-13/` (Phase 0
      struktura: falsyfikatory CLOSED, forbidden moves, risk register, anticipated outcomes)
- [x] Przeczytano `research/op-CE-H-3D-native-interaction-retry-2026-05-23/` (γ-1 CLEAN PASS:
      V_int = −2π v² n₁n₂ log(L/r₀); R²=0.9998; klasa kanału log/1-r; G_χ=1/(4πr))
- [x] Przeczytano `meta/CYCLE_KICKOFF_TEMPLATE.md` §1–§3 + `meta/CYCLE_LIFECYCLE.md`
      (claim_status taxonomy; output_type ⇒ ceiling; §X TGP-native check)
- [x] Przeczytano `meta/PRE_REGISTERED_FALSIFIERS.md` (numeracja PR; ostatni PR-022 APPENDED,
      PR-023 reserved bariogeneza, PR-024 reserved GST — niezarejestrowany przy HONEST_NEGATIVE)
- [x] Przejrzano skrypty: `tooling/scripts/ls9_epsilon0_nucleation.py`,
      `tooling/scripts/stability/dimensional_analysis.py`,
      `research/galaxy_scaling/gs9a_dimensional_structure.py` (istniejące narzędzia D-zależne)

**Sign-off:** Claudian (agent wykonawczy cyklu ND) @ 2026-06-13

## §1 — Kontrakt kickoff L1 (structural-first)

```yaml
L1_native:
  output_observable: "D ∈ ℤ₊ (wymiar przestrzenny generowanej przestrzeni); wtórnie N_sekt(D), Δ_D, Q(D), warunek bg-stabilizacji(D)"
  measurement_instrument: "D_obs = 3 — comparison-only anchor PO locku per-D; NIGDY input"
  native_coefs_constrained: ["N_sekt(D) z π_{k}(M_ord)", "A_D,B_D,C_D z stałych TGP {β,γ,Φ₀,λ}", "ν_D⁻¹ (nietrywialne przejście fazowe; d_c^Ising=4)", "warunek równowagi L*>0 bg-CE-H vs D"]
  falsification_rule: "F-ND-A/B/C/D/E (klasy + progi LOCKED 2026-06-13)"
  pre_registration_date: "2026-06-13"
L2_framework_reduction:
  target_frameworks: ["klasyfikacja defektów Mermin/Kibble (π_k(M) ↔ kowymiar k+1)", "Derrick virial scaling w D"]
  reduction_type: "analytical-exact"
  failure_disposition: "L1-stands"
L3_falsification_map:
  - { bound: "D_obs = 3", constrains: "selekcja wymiaru", window: "comparison-only po locku", status: "pending" }
  - { bound: "3 generacje materii (hyp:generacje)", constrains: "N_sekt(D)=3", window: "konsystencja deklaratywna; comparison-only", status: "pending" }
```

**Uwaga statusowa (CYCLE_LIFECYCLE §2.2):** `output_type: structural` ⇒ max claim_status **C**.
Cykl świadomie NIE aspiruje do A/A−: jego wartość = uczciwe rozstrzygnięcie statusu
istniejącego argumentu rdzenia (preferencyjny vs derived), nie nowa observable.

## §2 — Pytania (zbiory CLOSED; kolejność WYMUSZONA)

**Q-D1 (selekcja wymiaru — rdzeń cyklu):** czy maszyneria TGP (LIVE: S05 + Z₂ + U(1) + RP²;
zero nowych pól) wyróżnia D=3 jako wymiar generowanej przestrzeni — i z jaką siłą logiczną
(DERIVED / PREFERENTIAL / NONE)? Hipotezy (zbiór CLOSED):

- **H-SELECT-3-DERIVED** — koniunkcja niezależnych osi (topologia ∧ stabilność ∧ marginalność)
  wycina D=3 jednoznacznie, a każdy czynnik wycinający jest wyprowadzony z mechanizmu
  niezależnie od d (brak parametru dobranego pod d=3, brak Θ-czynnika uzasadnionego ex post).
- **H-SELECT-3-PREFERENTIAL** — D=3 jest wyróżniony, ale wyłącznie jako koniunkcja warunków
  „preferencyjnych" (co najmniej jeden czynnik Q(d) — typowo Θ(ν_d⁻¹) lub „fizyka trywialna
  w d≥4" — jest argumentem miękkim/jakościowym, nie ostrym progiem derived). Zgodne z własną
  deklaracją sek07a („preferencyjny"); osłabia „jedyny realistyczny wybór" do „najmocniejszy
  kandydat".
- **H-SELECT-OTHER** — uczciwy przegląd D=1..6 ujawnia inny wymiar równie/silniej uprzywilejowany
  (np. D=4 spełnia Δ_D ostrzej; jakiś D>3 niesie pełną hierarchię defektów) ⇒ argument sek07a
  obalony lub niejednoznaczny.
- **H-NO-SELECTION** — żadna oś nie wycina pojedynczego D bez wstrzyknięcia D_obs=3 ⇒
  selekcja wymiaru NIE jest własnością maszynerii TGP (HONEST_NEGATIVE).
- **GAP** — rozstrzygnięcie wymaga elementu spoza LIVE (np. dynamiki nukleacji bąbla w D
  wymiarach niewyprowadzonej; rozmaitości porządku nie ustalonej jednoznacznie) — luka
  deklarowana, zero nowych pól.

**Q-D2 (przegląd asymetrii / sortowania ND — TYLKO po rozstrzygnięciu Q-D1):** mechanizm
sortowania klasy H-SORT per D — czy selektywność działa tylko/najlepiej w jakimś D?
**Wiążące ograniczenie (SCOPING §2, kalibracja epistemiczna user):** H-SORT = mechanizm
roboczy DOPUSZCZAJĄCY spójność, NIE ustalona bariogeneza; **zakaz cytowania H-SORT jako
ustalonego** (forbidden move #12). Hipotezy (CLOSED): **SORT-PEAKS-3 / SORT-PEAKS-OTHER /
SORT-MONOTONE (brak piku) / GAP.** Q-D2 ma charakter INFORMATIONAL-strukturalny: jego wynik
NIE może podnieść claim_status cyklu ani „uratować" Q-D1.

**Reguła kolejności (WYMUSZONA):** Q-D1 najpierw (fast audyt §5). Jeśli Q-D1 = H-NO-SELECTION
lub GAP na poziomie fundamentów — Q-D2 nie jest wykonywane (zamknięcie po fazie Q-D1).

## §3 — Falsyfikatory (zbiory CLOSED; klasy i progi LOCKED 2026-06-13)

### F-ND-A — oś topologiczna (homotopia defektów; werdykt Phase 1)

Klasy (CLOSED): **TOPO-SELECTS-3 / TOPO-SELECTS-OTHER / TOPO-NO-SELECTION / GAP.**

Rachunek: dla każdego D ∈ {1,2,3,4,5,6} policz `N_sekt(D) = #{k∈[1,D] : π_{k−1}(M_ord) ≠ 0}`
dla JAWNIE ustalonej rozmaitości porządku M_ord, symbolicznie (sympy/grupy homotopii znanych
przestrzeni). **Pre-zarejestrowany spór do rozstrzygnięcia (NIE werdykt):** sek07a używa
`M_ord = SO(3)/Z₂` i podaje `π₂(SO(3)/Z₂) = ℤ`, podczas gdy **π₂ dowolnej grupy Liego (i jej
ilorazu przez grupę dyskretną) znika** (π₂(SO(3))=0 ⇒ π₂(SO(3)/Z₂)=0); źródłem defektów
punktowych π₂=ℤ jest rozmaitość `RP² = S²/Z₂` (π₂(RP²)=π₂(S²)=ℤ), jak w SCOPING §3a. F-ND-A
MUSI rozstrzygnąć, która rozmaitość jest GENUINE rozmaitością porządku TGP (z aksjomatów,
nie z wygody), i policzyć N_sekt(D) dla TEJ rozmaitości — bez uprzywilejowania D=3.

- **TOPO-SELECTS-3:** dla genuine M_ord pełna hierarchia (ściana ∧ struna ∧ punkt stabilny)
  oraz N_sekt = 3 zachodzi dla D=3 i **dla żadnego innego D∈{1,2,4,5,6}** — przy uczciwym
  liczeniu (te same reguły π_k dla każdego D). Wymaga: (i) defekty punktowe stabilne ⟺
  π_{D−1}(M_ord)≠0; (ii) jawny audyt D=4,5,6 (czy π_{D−1}≠0 też tam daje punkty — wtedy
  NIE-selekcja).
- **TOPO-SELECTS-OTHER:** N_sekt=3 lub pełna hierarchia realizuje się też (lub wyłącznie)
  dla D≠3.
- **TOPO-NO-SELECTION:** liczba sektorów rośnie monotonicznie z D (np. każde D>3 ma ≥3
  sektory) ⇒ „dokładnie 3" wymaga dodatkowego kryterium spoza topologii (przejście do F-ND-B).
- **GAP:** rozmaitość porządku nieustalona jednoznacznie z aksjomatów (SO(3)/Z₂ vs RP² vs
  inna) — luka deklarowana; rozstrzygnięcie poza LIVE.

### F-ND-B — oś stabilności (Derrick / potencjał trzech reżimów; werdykt Phase 1/2)

Klasy (CLOSED): **STAB-SELECTS-3-DERIVED / STAB-SELECTS-3-FITTED / STAB-SELECTS-OTHER / STAB-NO-SELECTION / GAP.**

Dwa pod-rachunki (oba wymagane):

1. **Derrick / bg-stabilizacja CE-H vs D:** skalowanie energii statycznego solitonu
   E_grad ∝ L^(D−2), E_pot ∝ L^D, człon stabilizujący tłem (CE-H, m_σ²=2λΦ₀²) ∝ L^(?).
   Czy istnieje okno D z równowagą L*>0 bez fine-tuningu? (Klasyczny Derrick: brak statycznych
   solitonów skalarnych dla D>1 bez stabilizatora; TGP omija przez bg — pytanie: dla jakich D
   bg-stabilizacja daje L*>0.)
2. **Audyt dyskryminanty `Δ_D` sek07a (KLUCZOWY value-blind test):** sek07a twierdzi
   `Δ_3 = 4B₃²−12A₃C₃>0` spełnione „z marginesem" przez `B₃/√(A₃C₃) ≈ 3.4`, a `Δ_4` „ostrzejszy"
   i fizyka d=4 „trywialna" (Θ(ν₄⁻¹)=0). **F-ND-B MUSI ustalić, czy A_D,B_D,C_D są wyprowadzone
   ze stałych TGP {β,γ,Φ₀,λ} jako funkcje D NIEZALEŻNIE (wtedy wartość 3.4 jest predykcją), czy
   wprowadzone jako parametry, których stosunek dobrano tak, by Δ_3>0 a Δ_2<0** (wtedy = fitted).
   Test obu kierunków: policz B_D/√(A_D C_D) dla D=2,3,4,5 z TYCH SAMYCH reguł i sprawdź, czy
   próg √(...) jest przekraczany selektywnie dla D=3.

- **STAB-SELECTS-3-DERIVED:** A_D,B_D,C_D derived z {β,γ,Φ₀,λ} niezależnie od d; Δ_D>0 z trzema
  reżimami zachodzi dla D=3 i nie dla D∈{1,2,4,5}; marginesy są PREDYKCJĄ, nie dopasowaniem.
- **STAB-SELECTS-3-FITTED:** D=3 wychodzi, ALE co najmniej jeden z {A_D,B_D,C_D, próg √} jest
  dobrany/uzasadniony ex post pod d=3 (np. „naturalne β~10⁻²" wybrane tak, by d=2 wypadło) ⇒
  selekcja przez stabilność = artefakt (zasila H-SELECT-3-PREFERENTIAL).
- **STAB-SELECTS-OTHER / STAB-NO-SELECTION / GAP** — analogicznie.

### F-ND-C — oś nukleacji + marginalności grawitacyjnej (werdykt Phase 2)

Klasy (CLOSED): **NUCL-MARG-SELECTS-3 / NUCL-MARG-SELECTS-OTHER / NUCL-MARG-NO-SELECTION / GAP.**

- Nukleacja w D wymiarach: akcja bąbla S_E(D) (thin-wall: S_E ∝ σ^D/ε^(D−1)-klasa), bariera
  i R_c jako funkcje D (reuse maszynerii `ls9_epsilon0_nucleation.py`, uogólnionej symbolicznie
  na D). Pytanie: czy nukleacja faworyzuje jakieś D (np. przez tempo Γ ∝ exp(−S_E(D)))?
- Marginalność grawitacyjno-wzrostowa w D ((1/2)v_c² = G_D M/R^(D−2)-klasa; analog γ-3 R=ct;
  forma wzrostu indicial p(D); sferyczność a_l(D)). Czy naddeterminowany układ {U,ρ₀} domyka się
  EXACT tylko dla pewnych D?
- **NUCL-MARG-SELECTS-3** wymaga: selekcja D=3 z nukleacji LUB marginalności, wyprowadzona
  symbolicznie, bez wstrzyknięcia D_obs=3 i bez dobierania anchorów pod 3. **NUCL-MARG-NO-SELECTION**
  (najbardziej prawdopodobny per §8): te osie są D-zależne ilościowo, ale nie wycinają
  pojedynczego D — księgowość działa „dla każdego D z odpowiednią stałą".

### F-ND-D — oś sortowania ND (Q-D2; werdykt późnej fazy; TYLKO po Q-D1)

Klasy (CLOSED): **SORT-PEAKS-3 / SORT-PEAKS-OTHER / SORT-MONOTONE / GAP.**
Wiążące: H-SORT = working-mechanism (forbidden #12). Werdykt F-ND-D NIE podnosi claim_status
cyklu. Rachunek: wydajność sortowania (KB1 topologiczna selekcja ścian kowymiaru 1 + KB2
cross-term) jako funkcja D; hipoteza pre-derywacyjna §4.4: D≥2 otwiera „kanał boczny"
(partnerzy się omijają) ⇒ wydajność maleje z D, podczas gdy defekty punktowe wymagają D≥3 —
czy iloczyn wycina D=3? (INFORMATIONAL.)

### F-ND-E — agregat „okno życia" (Phase FINAL)

Mapowanie (CLOSED) — werdykt zależy od koniunkcji osi A–C (D oraz „derived vs fitted"):

| F-ND-A | F-ND-B | F-ND-C | F-ND-E (agregat) |
|---|---|---|---|
| SELECTS-3 | SELECTS-3-DERIVED | SELECTS-3 lub NO-SELECTION | **DIM-3-DERIVED** — maszyneria TGP selekcjonuje D=3 mechanizmem; sek07a status: wzmocniony do „derived" (rewizja = user) |
| SELECTS-3 | SELECTS-3-FITTED | * | **DIM-3-PREFERENTIAL** — D=3 wyróżniony, ale ≥1 oś = preferencyjna/fitted; sek07a „jedyny realistyczny wybór" osłabiony do „najmocniejszy kandydat" (zgodnie z własną deklaracją „preferencyjny") |
| SELECTS-OTHER | * | * | **SEK07A-CHALLENGED** — uczciwy test D>3 wskazuje inny/niejednoznaczny wymiar; rozbieżność raportowana (np. π₂(SO(3)/Z₂) vs RP²) |
| NO-SELECTION | NO-SELECTION | NO-SELECTION | **NO-DIM-SELECTION (HONEST_NEGATIVE)** — selekcja wymiaru nie jest własnością maszynerii bez D_obs=3 |
| GAP (dowolna oś) | — | — | **GAP_CLOSURE** — luka fundamentalna deklarowana (np. rozmaitość porządku nieustalona) |

Decyzja o JAKIEJKOLWIEK rewizji statusu sek07a w rdzeniu (preferencyjny↔derived↔challenged) =
WYŁĄCZNIE user (Phase FINAL). Cykl nie edytuje rdzenia (forbidden #2).

## §4 — Analytical pre-derivation (§3.6.1–§3.6.9; PRZED rachunkiem)

### §4.1 — Homotopia defektów (oczekiwane wartości; standardowa klasyfikacja)

Defekt kowymiaru (k+1) w D-wymiarowej przestrzeni jest klasyfikowany przez π_k(M_ord)
(Mermin/Kibble). Defekt punktowy (kowymiar D) ⟺ π_{D−1}(M_ord) ≠ 0. Fakty matematyczne
(LOCKED matematyką, nie projektem):
- π_k(grupa Liego) — π₂ = 0 zawsze; π₁(SO(3)) = ℤ₂; π₃(SO(3)) = ℤ.
- π_k(RP²): π₀=0 (spójna), π₁=ℤ₂, π₂=ℤ (bo nakrycie uniwersalne S²).
- π_k(S²): π₁=0, π₂=ℤ, π₃=ℤ.
**Konsekwencja pre-derywacyjna (NIE werdykt):** twierdzenie „π₂(SO(3)/Z₂)=ℤ" (sek07a) jest
matematycznie podejrzane; „π₂(RP²)=ℤ" (SCOPING) jest poprawne. Genuine rozmaitość porządku
TGP musi być ustalona z aksjomatów (S05+Z₂+U(1)+RP²; dodatek app:E-spin-half) — to rachunek
F-ND-A, wynik OTWARTY. **Zakaz** wyboru rozmaitości tak, by dać π₂≠ℤ w D=3 (selection-by-result).

### §4.2 — Liczenie sektorów symetrycznie w D (enumeracja założeń)

N_sekt(D) liczone TYMI SAMYMI regułami dla każdego D: defekty wymiarów D−1,…,0 odpowiadają
π₀,…,π_{D−1}(M_ord). „Pełna hierarchia" (ściana ∧ struna ∧ punkt) wymaga π₀≠0 ∧ π₁≠0 ∧ π₂≠0
ORAZ D≥3 (by punkt = kowymiar 3 istniał). Uczciwy audyt D>3: czy w D=4,5,6 pojawiają się
DODATKOWE sektory (π₃,π₄,…≠0) — jeśli tak, „dokładnie 3" jest własnością D=3 tylko przy
obcięciu liczenia (które trzeba uzasadnić bez D_obs).

### §4.3 — Skalowanie Derricka i bg-stabilizacja (oczekiwania)

Statyczny soliton skalarny: E[Φ_λ(x)=Φ(x/λ)] = λ^(D−2) E_grad + λ^D E_pot. dE/dλ|₁ = 0 ⇒
(D−2)E_grad + D·E_pot = 0 — dla E_grad,E_pot>0 brak rozwiązania przy D≥2 (Derrick). TGP
wprowadza człon tła (CE-H): oczekiwana klasa stabilizatora zmienia bilans; pytanie ilościowe
= przy jakim skalowaniu członu bg równowaga L*>0 istnieje i dla jakich D. **Oczekiwanie
pre-derywacyjne (do testu):** stabilizacja bg może działać w paśmie D, niekoniecznie tylko D=3
— a zatem stabilność sama prawdopodobnie NIE wycina D=3 jednoznacznie (zasila H-SELECT-3-PREFERENTIAL).

### §4.4 — Sortowanie ND (konwencja + oczekiwanie; LOCKED przed rachunkiem)

Ściana = kowymiar 1 zawsze (każde D). Wymuszony porządek pary jest 1D wzdłuż normalnej do
ściany; w D≥2 partnerzy mogą się omijać kanałem bocznym ⇒ oczekiwana wydajność sortowania
maleje z D. Defekty punktowe wymagają D≥3. **Oczekiwanie:** iloczyn (wydajność ↓ z D) × (punkty
istnieją od D=3) mógłby dawać szczyt przy D=3 — ale to dokładnie ten typ „iloczynu wskaźników",
który grozi reverse-engineeringiem; stąd F-ND-D = INFORMATIONAL i nie podnosi claim_status.

### §4.5 — Konwencje i guardy (§3.6.6–§3.6.8; LOCKED)

- D traktowane jako symboliczny parametr; wszystkie wielkości liczone jako funkcje D, a D=3
  podstawiane DOPIERO w sekcji porównawczej.
- **Circularity guard FP (każdy sympy):** audyt, że „3" nie wchodzi do gałęzi decyzyjnej
  (T_pass) jako literał uprzywilejowany; D_obs=3, „żyjemy w 3D", „3 generacje" nieobecne w
  wyprowadzeniu (wyłącznie w sekcji comparison-only).
- Stałe TGP {β,γ,Φ₀,λ} symboliczne; m_σ=√(2λ)Φ₀ pochodna LOCKED CE-H. Zero wartości liczbowych
  dobranych pod d=3 (w szczególności „naturalne β~10⁻²" musi być uzasadnione niezależnie od d).

### §4.6 — Precyzja (§3.6.9)

Fakty homotopii = EXACT (algebra grup). Wielkości numeryczne (Δ_D, L*, S_E(D)) raportowane
symbolicznie + ±5% przy ewentualnym dopasowaniu numerycznym profilu. Każdy „margines"
(np. 3.4 > √3) musi mieć jawne źródło stałych.

### §4.7 — Klasyfikacja stałych (§3.6.13)

| Stała | Klasa | Nota |
|---|---|---|
| β, γ, Φ₀, λ | (α) TGP_FUNDAMENTAL | symboliczne; zero wartości dobranych pod d=3 |
| m_σ = √(2λ)Φ₀ | (α) pochodna | LOCKED CE-H |
| A_D, B_D, C_D | pochodne (do audytu) | MUSZĄ być funkcjami {β,γ,Φ₀,λ,D}; jeśli wolne ⇒ STAB-FITTED |
| D | parametr całkowity | symboliczny; podstawienie 3 tylko w comparison-only |
| D_obs = 3 | comparison-only | NIGDY input; guard FP w każdym sympy |
| ν_D, η_D (wykładniki krytyczne) | LITERATURE (Ising/O(n)) | d_c=4 fakt RG; użycie deklaratywne |
| **Nowe pola / stałe fundamentalne** | **budżet: 0** | luka = GAP |

## §5 — Plan faz (audyt-first)

| Faza | Zakres | Werdykt | Gate |
|---|---|---|---|
| **Phase 1 — FAST AUDYT (Q-D1 topologia+stabilność)** | F-ND-A (homotopia N_sekt(D), D=1..6, genuine M_ord, audyt SO(3)/Z₂ vs RP²) → F-ND-B audyt Δ_D (derived vs fitted A,B,C; B/√(AC) dla D=2..5) + Derrick/bg-stabilizacja vs D | F-ND-A, F-ND-B | „działaj" |
| Phase 2 (jeśli Q-D1 nie = NO-SELECTION/GAP) | F-ND-C nukleacja S_E(D) + marginalność grawitacyjna w D (symbolicznie) | F-ND-C | „działaj" |
| Phase 3 (TYLKO jeśli Q-D1 rozstrzygnięte i user chce Q-D2) | F-ND-D sortowanie ND (INFORMATIONAL; H-SORT working-mechanism) | F-ND-D | „działaj" |
| Phase FINAL | agregat F-ND-E; DOUBTS register; dyspozycja statusu sek07a = WYŁĄCZNIE user | F-ND-E | „działaj" |

Każda faza: PhaseN_derivation.md + PhaseN_sympy.py/.txt (PASS/FAIL per FP; 0 hardcoded T_pass;
circularity-guard FP: audyt „3" jako uprzywilejowanego literału). NO-SELECTION/GAP w Phase 1
⇒ od razu Phase FINAL (zamknięcie), bez Q-D2.

## §6 — Forbidden moves (LOCKED; 14)

1. **Użycie D_obs=3 (lub „żyjemy w 3D", „3 generacje", „3 kolory") jako inputu** wyprowadzenia
   — wyłącznie comparison-only po locku per-D; guard FP w każdym sympy.
2. Edycja rdzenia (sek07a oraz S05/Z₂/U(1)/RP²/γ-*/CE-H — read-only); rewizja statusu sek07a
   = dyspozycja user w Phase FINAL, nie edycja przez cykl.
3. **Cytowanie konkluzji sek07a („d=3 jedyny realistyczny wybór") jako potwierdzenia** selekcji
   — sek07a = obiekt audytu, nie input.
4. Nowe pola spoza S05+Z₂+U(1)+RP² lub nowe stałe fundamentalne (luka = GAP).
5. Dobieranie A_D,B_D,C_D, progu √(...), lub „naturalnych" wartości stałych (β~10⁻²) tak, by
   wynik wypadł na D=3 (value-blind: stałe symboliczne; podstawienia uzasadnione niezależnie od d).
6. Liczenie N_sekt / Δ_D / S_E innymi regułami dla D=3 niż dla D≠3 (asymetria liczenia).
7. Obcięcie przeglądu do D≤4 bez uczciwego audytu D=5,6 (sek07a kończy na d=4 — cykl musi iść dalej).
8. Wybór rozmaitości porządku M_ord pod pożądane π₂ (selection-by-result); M_ord ustalone z aksjomatów.
9. Hardcoded T_pass / „3" w gałęzi decyzyjnej sympy.
10. Zmiana klas/progów falsyfikatorów ex post (zbiory CLOSED; IMMUTABLE).
11. Podniesienie claim_status powyżej C przez Q-D2/sortowanie lub przez „ładny" iloczyn wskaźników.
12. **Cytowanie H-SORT jako ustalonej bariogenezy** (kalibracja epistemiczna user, SCOPING §2).
13. Miękkie domknięcie / przeciąganie cyklu po NO-SELECTION (HONEST_NEGATIVE = sukces metodologiczny).
14. Scope-creep do op-phi-radiative-dof-audit lub innych cykli kolejki (osobna pre-rejestracja).

## §7 — Risk register

| ID | Ryzyko | Severity | Mitygacja |
|---|---|---|---|
| R-ND-1 | **Reverse-engineering Q(d) pod D_obs=3** — najsilniejsza pokusa (rdzeń już konkluduje D=3) | META/HIGH | D_obs comparison-only; klasy CLOSED; audyt „derived vs fitted" wbudowany w F-ND-B; guard FP |
| R-ND-2 | Niejednoznaczność rozmaitości porządku (SO(3)/Z₂ vs RP² vs S²) ⇒ π₂ różne | HIGH | F-ND-A ustala M_ord z aksjomatów; rozbieżność = GAP, nie wybór wygodny |
| R-ND-3 | Stabilność (Derrick/bg) prawdopodobnie NIE wycina pojedynczego D ⇒ pokusa dodania Θ-czynnika | HIGH | §4.3 oczekiwanie zapisane; każdy Θ-czynnik musi mieć ostry próg derived, inaczej = PREFERENTIAL |
| R-ND-4 | „Naturalne wartości stałych" (β~10⁻²) jako ukryty fit | MED-HIGH | forbidden #5; stałe symboliczne; podstawienia uzasadnione niezależnie od d |
| R-ND-5 | Asymetria audytu D<3 vs D>3 (łatwo wykluczyć D=1,2; trudniej D=4,5,6 uczciwie) | MED-HIGH | forbidden #7; obowiązkowy audyt D=5,6; te same reguły |
| R-ND-6 | Czynnik „fizyka trywialna w d≥4 / pole średnie" (Θ(ν⁻¹)) jest miękki (d_c=4 to fakt RG, ale „trywialna" ≠ „niemożliwa") | MED | klasyfikacja jako PREFERENTIAL, nie DERIVED, jeśli to jedyny nóż na d=4 |
| R-ND-7 | Nukleacja w D: brak gotowej symbolicznej maszynerii (ls9 nie liczy selekcji D) | MED | F-ND-C może = GAP (uczciwie); zero dorabiania |
| R-ND-8 | Pokusa „TGP wyprowadza dlaczego 3D" jako efektowny wynik | META/HIGH | progi LOCKED; anticipated outcomes §8; ceiling C; rewizja sek07a = user |
| R-ND-9 | Mieszanie Q-D2 (sortowanie/H-SORT) z Q-D1 by „dopełnić" selekcję | MED | kolejność wymuszona; F-ND-D INFORMATIONAL; forbidden #11/#12 |
| R-ND-10 | Hipoteza „3 generacje ↔ 3 sektory" cytowana jako dowód selekcji | MED | comparison-only (L3 map); konsystencja deklaratywna, nie input |

## §8 — Anticipated outcomes (INFORMATIONAL; kierunki pokusy zapisane PRZED rachunkiem)

1. **Najbardziej prawdopodobny wynik: DIM-3-PREFERENTIAL (H-SELECT-3-PREFERENTIAL).** Oczekiwanie:
   oś topologiczna może autentycznie wyróżniać D=3 (pełna hierarchia ściana/struna/punkt wymaga
   dokładnie D=3 dla rozmaitości z π₂≠0, np. RP²), ALE oś stabilności (Derrick/bg) i marginalności
   prawdopodobnie działają w paśmie D, a kluczowy nóż na d=4 w sek07a (Θ(ν₄⁻¹)=0, „fizyka trywialna")
   jest argumentem miękkim, nie ostrym progiem. ⇒ Uczciwy werdykt: D=3 jest najmocniejszym
   kandydatem (preferencyjnym), zgodnie z własną deklaracją sek07a, a „jedyny realistyczny wybór"
   jest sformułowaniem mocniejszym niż uzasadnia rachunek. To NIE jest porażka — to doprecyzowanie
   statusu (C, structural).
2. **Kierunek pokusy #1:** reverse-engineering Q(d) — przyjęcie czynników Θ i stałych A,B,C tak,
   by Q(3)=3 i Q(d≠3)=0. Antidotum: audyt derived-vs-fitted (F-ND-B), guard FP, klasy CLOSED.
3. **Kierunek pokusy #2:** rozstrzygnięcie sporu π₂(SO(3)/Z₂) vs π₂(RP²) „po stronie wygody"
   (wybór rozmaitości dającej punkty w D=3). Antidotum: M_ord z aksjomatów; rozbieżność = GAP.
4. **Kierunek pokusy #3:** zatrzymanie przeglądu na d=4 (jak sek07a). Antidotum: obowiązkowy D=5,6.
5. **Możliwy wynik mocny (DIM-3-DERIVED):** gdyby genuine M_ord = RP² (π₂=ℤ) dawało pełną
   hierarchię WYŁĄCZNIE w D=3 z uczciwym audytem D>3, a stabilność/marginalność były przynajmniej
   konsystentne — to byłby autentyczny strukturalny wynik (ceiling C). Wymagałby wzmożonego audytu
   (każdy czynnik derived), bo to „dlaczego 3D" — najbardziej kuszący rezultat.
6. **Możliwy wynik negatywny (SEK07A-CHALLENGED / NO-DIM-SELECTION):** gdyby D>3 niósł równie
   pełną hierarchię lub żadna oś nie wycinała pojedynczego D — argument rdzenia wymagałby rewizji.
   Raportowany uczciwie; rozbieżność z sek07a (np. błąd π₂) = wartościowy wkład, nie porażka.

## §9 — Outcome sets per FP (§3.6.11)

Każdy FP Phase 1–3 deklaruje zbiór {PASS / PARTIAL_compute / FAIL} (+ PARTIAL_concept_mismatch
tam, gdzie TGP-native nie ma odpowiednika — jawne uzasadnienie per §3.6.11(b)); budżet
PARTIAL_compute: max 1 na cykl. Werdykty F-ND-* wyłącznie w klasach CLOSED §3. Fakty homotopii =
EXACT (brak PARTIAL dla algebry grup).

## §10 — Anti-Lakatos compliance (Phase 0)

Zbiory hipotez/klas CLOSED przed rachunkiem ✓ · D_obs=3 comparison-only z guardem FP ✓ ·
**obiekt audytu (sek07a) jawnie zadeklarowany jako hipoteza pod testem, nie input** ✓ ·
audyt derived-vs-fitted wbudowany (F-ND-B) ✓ · uczciwy test obu kierunków D<3 i D>3 (D=5,6
obowiązkowe) ✓ · spór matematyczny π₂(SO(3)/Z₂) vs π₂(RP²) zarejestrowany jako OTWARTY (nie
rozstrzygnięty pod wygodę) ✓ · 0 nowych pól/stałych ✓ · 0 edycji rdzenia (rewizja statusu = user) ✓ ·
H-SORT working-mechanism (zakaz cytowania jako ustalonej bariogenezy) ✓ · ceiling claim_status C
zadeklarowany (brak inflacji statusu przez Q-D2) ✓ · anticipated outcomes z kierunkami pokusy ✓ ·
kolejka po cyklu = decyzja user (bez scope-creep) ✓.

**Phase 0 LOCKED 2026-06-13. Następny krok: Phase 1 FAST AUDYT (F-ND-A + F-ND-B) — wymaga user „działaj".**
