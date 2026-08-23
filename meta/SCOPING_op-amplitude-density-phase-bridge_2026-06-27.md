---
title: "SCOPING — op-amplitude-density-phase-bridge: czy K_eff przechodzi Φ^{-1} → Φ^{+4} (reżim amplitudowy σ → reżim gęstościowy/metryczny Φ)? Mapa terenu + rejestr ryzyk PRZED Phase 0."
date: 2026-06-27
type: meta-scoping
status: "PRE-PHASE-0 NOTE (notatka robocza; NIE pre-rejestracja; NIE cykl; zero werdyktów; wymaga własnego Phase 0 + 'działaj')"
origin: "User 2026-06-27: reinterpretacja werdyktu (B) #49 — α_eff=−1/2 jako przedmetryczny reżim amplitudowy (σ pierwotne, Φ=σ² kompozyt), α=2 jako pometryczny atraktor metryczny. Trzy osie do zbadania: (1) faza przejściowa + wykres ewolucji wskaźnika, (2) realność K_eff(Φ,ℓ): Φ^{-1}→Φ^4 i sens fizyczny (narodziny metryki?), (3) dowód stabilności (≈ oś 1 inaczej ubrana)."
parent_cycle: "[[../research/op-CG-alpha-eff-convergence-2026-06-26/Phase_FINAL_close.md]] CLOSED-RESOLVED (B) REFUTED-SUBSTRATE (LOCKED, IMMUTABLE)"
relates_to:
  - "[[HONEST_FRAMING_UV_CG_ROOTS.md]] (cztery korzenie UV/CG)"
  - "[[../research/op-emergent-metric-from-interaction-2026-05-09/Phase_FINAL_close.md]] (emergent metric, CLOSED)"
  - "[[../research/op-nucleation-dimensionality-2026-06-13/Phase_FINAL_close.md]] (nukleacja, DIM-3-PREFERENTIAL)"
  - "[[../core/_meta_latex/status_map.tex]] l.72 (selekcja na gęstości), l.489-501 (α=1 vs α=2)"
anti_lakatos_note: "Mapa terenu / pytanie strukturalne. ŻADNE twierdzenie nie nabywa statusu przez tę notę. Werdykt (B) #49 (α_eff=−1/2 z substratu, CLOSED-NEGATIVE), status_map l.72 (α=2 = selekcja na gęstości, NIE derywacja), #42 ledger (α=2 ∈ N_axiom=6), #48 K_geo, #37 c₀ — wszystkie LOCKED. Pre-derywacja §3 = OCZEKIWANIE, nie próg. Endpoint +4 NIE może być dostrajany pod wynik."
tags:
  - scoping
  - amplitude-density-bridge
  - kinetic-exponent-crossover
  - emergent-metric
  - pre-metric-phase
  - conformal-stability-selection
  - UV-CG-roots-followup
  - anti-Lakatos-LOCKED
---

# SCOPING — most amplituda→gęstość: czy K_eff biegnie Φ⁻¹ → Φ⁺⁴?

> **Notatka robocza (etap analizy).** Cel: ocenić, czy reinterpretacja werdyktu (B) #49
> jako **dwufazowego crossoveru** (reżim amplitudowy σ → reżim gęstościowy/metryczny Φ)
> jest (a) fizycznie nietrywialna, (b) wykonalna obliczeniowo, (c) zgodna z LOCKED rdzeniem.
> **To NIE jest start cyklu** — to rozpoznanie terenu + jawny rejestr największych problemów,
> które trzeba rozbroić PRZED Phase 0.

## §1 — Teza do zbadania (jedno zdanie)

`α_eff = −½` (`K∝Φ⁻¹`) opisuje **przedmetryczny reżim amplitudowy** substratu (pole pierwotne
`σ`, `Φ=σ²=⟨ŝ²⟩` to kompozyt); `α = 2` (`K∝Φ⁴`) opisuje **pometryczny reżim gęstościowy**,
w którym `Φ` stało się samodzielnym polem gęstości przestrzeni pod warunkami stabilności
konforemnej. Pytanie cyklu: **czy istnieje realny (nie-tautologiczny) operator `K_eff(Φ, ℓ)`,
który interpoluje `Φ⁻¹ → Φ⁴` wraz ze skalą/kondensacją — i czy ten crossover to fizycznie
„narodziny metryki", czy coś subtelniejszego?**

## §2 — Trzy osie badawcze (od usera, doprecyzowane)

| Oś | Pytanie | Pożądany deliverable | Trudność |
|---|---|---|---|
| **A — faza przejściowa** | Czy crossover `−½ → +2` jest gładki, czy dynamicznym skokiem? Jak wygląda ewolucja wskaźnika? | **Wykres `e_eff(ℓ)` (lub `e_eff(⟨Φ⟩)`)** — wykładnik efektywny vs skala/gęstość | WYSOKA |
| **B — realność `K_eff(Φ,ℓ)`** | Czy interpolacja `Φ⁻¹→Φ⁴` jest realnym obiektem, czy artefaktem zmiennej? Co znaczy fizycznie — narodziny metryki czy coś subtelniejszego? | Jawna konstrukcja `K_eff(Φ,ℓ)` + interpretacja (`√(−g)`-dressing vs wymiar anomalny vs zmiana DOF) | WYSOKA |
| **C — dowód stabilności** | Czy `α=2` jest wyróżnione jako atraktor/wartość samostabilizująca w rodzinie `α`? | Kryterium stabilności konforemnej selekcjonujące `α=2` (lub falsyfikacja: selekcjonuje inne α / nie selekcjonuje) | BARDZO WYSOKA |

**Uwaga usera (trafna):** oś C ≈ oś A „inaczej ubrana". Zgoda — w obu chodzi o to samo:
**co fizycznie ustawia wartość docelową wykładnika.** A patrzy na to dynamicznie (ewolucja),
C strukturalnie (atraktor). Warto trzymać je razem; rozdzielić tylko jeśli Phase 0 pokaże,
że wymagają różnych narzędzi.

## §3 — Pre-derywacja / oczekiwanie analityczne (OCZEKIWANIE, nie próg; zapisane PRZED Phase 0)

**Fakt zakotwiczony (LOCKED, #49):** dla `Φ=σ^{2p}`, kanoniczna kinetyka `(∇σ)²` daje
`K(Φ)∝Φ^{e(p)}`, `e(p)=1/p−2`. Kompozyt fizyczny `p=1` ⟹ `e=−1` (`α_eff=−½`). `α=2` (`e=+4`)
wymaga `p=1/6` — bez sensu substratowego. To jest **czysty jakobian** zamiany zmiennych.

**Konsekwencja krytyczna (sedno ryzyka R1):** jeśli `Φ=σ²` obowiązuje **globalnie i dokładnie**,
to `K(Φ)=1/(4Φ)` oraz kanoniczne `(∇σ)²` to **DOSŁOWNIE TA SAMA TEORIA** w dwóch zmiennych.
Wtedy **nie ma czego „flow-ować"** — wykres `e_eff(ℓ)` byłby stały (=const wyznaczona przez
wybór zmiennej), a crossover byłby tautologią. **Realny crossover wymaga, by relacja `σ↔Φ`
sama zależała od skali/fazy** (`Φ = F_ℓ(σ)`, `F_ℓ` nietrywialnie zależne od ℓ).

**Dwie kandydujące drogi do nietrywialności (do rozstrzygnięcia w Phase 1):**

- **Droga A (RG / wymiar anomalny) — OCZEKIWANA jako ZABLOKOWANA.** „Wykładnik biegnie"
  przez fluktuacje. Wymaga `Δe = 5`, `η ~ O(5)`; WF 3D daje `η≈0,036`. **#49 zamknął to
  negatywnie.** Jeśli oś A sprowadza się do tego — wynik z góry znany (NEGATIVE). Forbidden:
  re-litygacja #49.
- **Droga B (czynnik geometryczny/ramowy) — OCZEKIWANA jako jedyna żywa.** `+4` NIE jest
  wymiarem anomalnym, lecz **dressingiem przez emergentną metrykę** (`√(−g)=c₀ψ`, rama
  konforemna `g_{μν}=f(Φ)η_{μν}`, `f∝Φ`), który **fizycznie nie istnieje przed kondensacją**.
  To omija ścianę `η~O(5)`, bo nie żąda od fluktuacji wygenerowania `Δe=5`. **Hipoteza robocza
  (do testu, NIE claim):** `K_eff(Φ,ℓ) = [czynnik amplitudowy 1/Φ] × [czynnik metryczny narastający
  z ℓ]`, gdzie czynnik metryczny włącza się przy nukleacji i przesuwa efektywny wykładnik
  `−1 → +4`. **Rachunek nadrzędny nad tą hipotezą.**

**Strukturalne wsparcie Drogi B (cytaty, nie dowód):** `sek08c` — `φ⁴` współwystępuje z
`√(−g)=c₀ψ` i pojawia się przez interaction-driven emergent metric ([[../research/op-emergent-metric-from-interaction-2026-05-09/]] CLOSED). Cykl nukleacyjny
([[../research/op-nucleation-dimensionality-2026-06-13/]]) ma dynamikę „bąbla generującego
przestrzeń" (W-ND-4) — naturalny kandydat na „fazę przejściową" osi A.

## §4 — NAJWIĘKSZE PROBLEMY (rejestr ryzyk — to jest sedno tej notatki)

> Uporządkowane wg dotkliwości. **R1, R4, R5 są egzystencjalne** dla cyklu — jeśli się ich
> nie rozbroi w Phase 0, oś badawcza jest albo tautologią, albo sprzecznością z rdzeniem.

### R1 (KRYTYCZNE — tautologia zmiennej) — „nie ma czego flow-ować" — ✅ ROZBROJONE WARUNKOWO 2026-06-27
Jeśli `Φ=σ²` globalnie, `−½` i kanoniczne `σ` to ta sama teoria. Wykres ewolucji wskaźnika
(oś A) byłby wtedy **wykresem stałej**, a „crossover" — przepisaniem tego samego dwoma
literami. **Co rozbraja:** Phase 0 MUSI jawnie zdefiniować, **co czyni relację `σ↔Φ`
zależną od skali ℓ** (np. `Φ` jako parametr porządku z własną, generowaną dynamicznie
sztywnością, a nie sztywno `σ²`). Bez tego — cykl nie ma przedmiotu.
> **WYNIK mini-testu R1** ([[../research/op-amplitude-density-phase-bridge-2026-06-27/Phase0_R1_FINDINGS.md]]):
> goły bulk-crossover = **TAUTOLOGIA** potwierdzona (T1+T3': obie ramy → to samo pole kanoniczne
> `χ=√2√Φ`). Teza żywa TYLKO jako **niezgodność par `(K_sub,V_sub)` vs `(K_TGP,V_M911)`** w ramie
> kanonicznej, zakotwiczona w `Φ=0` (nukleacja). Nowe wejście egzystencjalne: **`V_sub`** (potencjał
> substratu) — bez niego test niewykonalny. Naiwny dressing metryczny też niewystarcza (T2: `m=5` vs TGP `m=1`).

### R2 (WYSOKIE — droga A jest martwa) — re-litygacja #49
Każda wersja „wykładnik biegnie przez wymiar anomalny" jest **zamknięta negatywnie (#49,
η~O(5))**. Ryzyko: niechcący odtworzyć #49 i „odkryć" znany NEGATIVE. **Co rozbraja:**
explicit commit na Drogę B (geometryczny dressing) + forbidden move: zakaz traktowania `Δe`
jako wymiaru anomalnego.

### R3 (WYSOKIE — cyrkularność emergent-metric) — skąd naprawdę bierze się +4
Droga B opiera `+4` na `√(−g)` i `f(Φ)∝Φ`. Ale **metryka w TGP była wyprowadzana przy
założeniach o strukturze pola.** Ryzyko: `+4` „wyemerguje" tylko dlatego, że metryka już
zakłada `α=2`. **Co rozbraja:** audyt [[../research/op-emergent-metric-from-interaction-2026-05-09/]]
— sprawdzić, czy derywacja metryki jest **niezależna od `α=2`** (czy nie zakłada tego, co
chcemy pokazać). Jeśli zakłada → Droga B jest cyrkularna i upada.

### R4 (KRYTYCZNE — α ma fizyczny ząb) — sektor mas nie jest frame-invariant
Obrona „to ta sama teoria w innych zmiennych" **zawodzi dokładnie tam, gdzie boli**:
`status_map` l.499 — PPN, `κ`, `N_e`, `n_s` są α-NIEZALEŻNE, **ale ODE solitonu (→ masy
leptonów) JEST α-zależne.** Czyli `−½` vs `+2` to realna różnica fizyczna w sektorze materii,
nie tylko wybór ramy. **Co rozbraja:** dwufazowy obraz MUSI pokazać, jak sektor mas (używający
efektywnie `+2`/`+1`) łączy się z `−½`-substratem bez rozjazdu predykcji mas. To najtrudniejszy
test spójności.

### R5 (KRYTYCZNE — wewnętrzny konflikt α=1 vs α=2) — w grze są TRZY wykładniki, nie dwa
`status_map` l.489-501: w sektorze solitonowym/materii rdzeń **już znalazł `α=1` (φ²)
preferowane** („bariera duchowa w α=2 wyklucza φ-FP"; `cor:alpha1-preferred`, ex166/167),
podczas gdy sektor kinetyczno-metryczny używa `α=2` (φ⁴). Substrat daje `α_eff=−½`.
**Czyli czysty obraz usera „−½ → +2" pomija trzecią wartość: `α=1`.** To jednocześnie:
- **szansa:** TGP już żyje z różnymi efektywnymi α w różnych sektorach — multi-reżim jest
  precedensem, nie herezją;
- **zagrożenie:** trójka {−½, 1, 2} burzy prostą narrację dwufazową; dowód stabilności (oś C)
  musi wyjaśnić, dlaczego **różne sektory siadają na różnych α**, a nie jeden wszechświat na 2.
**Co rozbraja:** Phase 0 musi zintegrować l.489-501 jako wejście, NIE pominąć. Jeśli się nie
da pogodzić — teza dwufazowa wymaga przeformułowania na „multi-sektorową".

### R6 (ŚREDNIE — definicja ℓ na przedmetrycznym substracie)
`K_eff(Φ,ℓ)` zakłada pojęcie skali/długości `ℓ`. Ale **przed emergencją metryki „długość"
i „pęd" nie są dobrze zdefiniowane** — są kompozytami grafu Γ. Ryzyko: `ℓ` przemyca metrykę,
której istnienie jest częścią tezy. **Co rozbraja:** zdefiniować `ℓ` czysto kombinatorycznie
(skala blokowania uśredniania na Γ, liczba węzłów/krok RG), bez odwołania do metryki.

### R7 (ŚREDNIE — anti-Lakatos: dostrajanie endpointu)
Pokusa: zbudować `K_eff` tak, by **z definicji** lądował na `+4`. To by było moving-goalposts.
**Co rozbraja:** endpoint `+4` musi być **wyliczony**, nie wbudowany; circularity guard
(zakaz odwrotnego dopasowania do α=2, jak w #49 §3).

### R8 (NISKIE — „stabilny" ≠ „jednoznacznie selekcjonowany")
Oś C może pokazać, że `α=2` jest stabilne — ale to NIE znaczy, że jest **jedyne** stabilne
(α=1 też może być). Argument `α=100`-glassy daje tylko górne odcięcie, nie unikalność.
**Co rozbraja:** sformułować oś C jako selekcję (najlepiej: jedyność w paśmie), nie samą
stabilność punktową; spójne z R5.

## §5 — Co realnie ten cykl MOŻE i CZEGO NIE MOŻE osiągnąć (kalibracja oczekiwań)

- **NIE domknie mostu CG** ani nie „wyprowadzi" α=2 z substratu. `N_free` bez zmian.
  Werdykt (B) #49 pozostaje LOCKED.
- **MOŻE (sukces):** przekształcić status α=2 z „aksjomat-łatka" w „**atraktor IR /
  selekcja przez stabilność konforemną**" — mocniejsza, uczciwsza pozycja do publikacji
  (analog: negatyw G_SPA=48 jako pełnoprawny wynik metodologiczny).
- **MOŻE (porażka = też wynik):** pokazać, że crossover jest tautologią zmiennej (R1) lub
  cyrkularny (R3) lub sprzeczny z sektorem mas (R4) — wtedy HONEST_NEGATIVE, a dwufazowa
  interpretacja zostaje odrzucona jawnie.

## §6 — Szkic formatu i kosztu (wzór, do LOCK w Phase 0)

- **Format:** 2-3 fazy merytoryczne. Phase 1 = oś B (konstrukcja `K_eff(Φ,ℓ)` + test
  tautologii R1, Droga B vs A); Phase 2 = oś A (numeryczny wykres `e_eff(ℓ)` na Γ, FSS jak #49);
  Phase 3 = oś C (kryterium stabilności/selekcji, integracja z R4/R5).
- **Narzędzia:** sympy (jakobian, `√(−g)`-dressing, znaki/stabilność konforemna), reuse
  estymatora FSS z #49 (`Phase1_fss.py`, IMMUTABLE — kopiować, nie modyfikować).
- **Budżet nowych stałych: 0.** Output type: **structural** (brak observable z jednostkami;
  ceiling claim_status = C / STRUCTURAL_VERIFIED).
- **Dane obserwacyjne: 0** (czysto strukturalny; brak PR-### chyba że pojawi się falsyfikowalny
  styk z sektorem mas — wtedy osobny wpis).

## §7 — Forbidden moves (szkic, do LOCK w Phase 0)

1. Zakaz re-litygacji #49 / status_map l.72 / #42 / #48 / #37 (LOCKED, IMMUTABLE).
2. Zakaz traktowania `Δe = −1→+4` jako wymiaru anomalnego (Droga A martwa; tylko Droga B).
3. Zakaz dostrajania `K_eff` tak, by z definicji dawało `+4` (R7; endpoint wyliczany).
4. Zakaz odwołania do metryki przy definicji skali `ℓ` (R6; `ℓ` kombinatoryczne na Γ).
5. Zakaz pominięcia l.489-501 (α=1 vs α=2) — MUSI być wejściem (R5).
6. Zafiksować konwencję `Φ=σ^{2p}` i znaczenie „skali ℓ" PRZED jakimkolwiek rachunkiem.
7. Zakaz claimowania „derywacji α=2" — maksimum to „selekcja przez stabilność" (anti-Lakatos).

## §8 — Rejestracja follow-upu

**`op-amplitude-density-phase-bridge`** — REGISTERED jako kandydat (NOT activated; wymaga
własnego Phase 0 + user „działaj"). §2-§7 = szkic zakresu. Rozwija reinterpretację werdyktu
(B) #49. **Egzystencjalne bramki do przejścia w Phase 0 (inaczej cykl nie rusza):** rozbrojenie
R1 (tautologia), R3 (cyrkularność emergent-metric), R4 (sektor mas), R5 (α=1 vs α=2).

> **Relacja do czterech korzeni UV/CG ([[HONEST_FRAMING_UV_CG_ROOTS.md]]):** ten cykl NIE
> domyka mostu Γ→Φ. Atakuje wyłącznie **interpretację** jednego korzenia (α=2): czy `−½`-substrat
> i `+2`-makro to dwie fazy jednego obiektu, czy nieusuwalny rozjazd. Pozostałe korzenie
> (c₀ #37, 𝒜 #43, K_geo #48) nietknięte.

---

## Status

🟡 **PRE-PHASE-0 NOTE — notatka robocza, utworzona 2026-06-27.** Zero werdyktów cyklu. Nie cykl.

**Postęp diagnostyczny 2026-06-27** (folder [[../research/op-amplitude-density-phase-bridge-2026-06-27/]]):
- ✅ **R1 (tautologia)** — rozbrojone: bulk-crossover = tautologia (T1+T3'); mapa `ŝ→⟨ŝ²⟩` to
  prawdziwy coarse-graining (T_D), więc teza żywa, ale to droga #49 (`−½`).
- ✅ **Bramka #0 V_sub** — CLEARED: `U(φ)=(β/3)φ³−(γ/4)φ⁴` + warstwy mikro/Landau (dodatekB).
- ✅ **Kluczowy wynik Phase 1 first-pass** ([[../research/op-amplitude-density-phase-bridge-2026-06-27/Phase1_Vsub_FINDINGS.md]]):
  `α=2` to **aksjomat v2** (sprzężenie geometryczne `K_ij=J(φ_iφ_j)²`), NIE derywacja z mikro
  `H_Γ` (bilinear bond → kanoniczna kinetyka → `K∝Φ⁻¹` = #49). Rdzeń sam to deklaruje
  (rem:B-v2-status). Teza dwufazowa usera = **poprawna interpretacja**; „most" = wyprowadzenie
  aksjomatu geometrycznego z `H_Γ` (#49 mówi: near-Gaussian tego nie robi).
- 🔜 **Pozostałe bramki:** R3 (emergent-metric), R4 (sektor mas), R5 (teraz **trzy** wykładniki:
  `−½` mikro / `α=1` H_Γ-energia / `α=2` aksjomat v2).

**Pytanie rdzenia PRAWDZIWEGO cyklu (zaostrzone):** czy istnieje JAKIKOLWIEK (non-Gaussian /
strong-coupling) coarse-graining `H_Γ → K_ij=J(φ_iφ_j)²`?

✅ **ROZSTRZYGNIĘTE 2026-06-27** — cykl [[../research/op-CG-Kij-from-Hgamma-2026-06-27/Phase_FINAL_close.md]]:
**F-CGK-D = NON-DERIVABLE.** `α=2` nie wyprowadza się z mikro `H_Γ` (kompletna mapa obstrukcji:
Gaussian→−½ / η→unitarność / RG→irrelewancja wszystkich `O_n` / α=2 wymaga szesciopolowego,
wolnego, irrelewantnego bondu). Ratyfikuje #49/#39/status_map l.72. **Korzeń α=2 domknięty
analitycznie.** Bonus: R5 (trzy wykładniki) wyjaśniony jako artefakt `H_Γ≠F_kin`. DOUBT W-CGK-1:
„headline #49" Δe=5 miesza ramy (spójnie 4/2); werdykt (B) robust.

**Cross-references:**
- [[../research/op-CG-alpha-eff-convergence-2026-06-26/Phase_FINAL_close.md]] (#49 — (B) REFUTED-SUBSTRATE, LOCKED)
- [[HONEST_FRAMING_UV_CG_ROOTS.md]] (cztery korzenie UV/CG)
- [[../research/op-emergent-metric-from-interaction-2026-05-09/Phase_FINAL_close.md]] (emergent metric — audyt cyrkularności R3)
- [[../research/op-nucleation-dimensionality-2026-06-13/Phase_FINAL_close.md]] (nukleacja — kandydat „fazy przejściowej")
- [[../core/_meta_latex/status_map.tex]] l.72 (selekcja na gęstości), l.489-501 (α=1 vs α=2, R5)
- [[CYCLE_KICKOFF_TEMPLATE.md]] · [[templates/op-cycle-kickoff-template-v2-2026-05-11.md]] (wzór README przy aktywacji)
