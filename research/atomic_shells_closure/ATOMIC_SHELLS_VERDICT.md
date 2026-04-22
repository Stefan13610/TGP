# Atomic Shells — Verdict 2026-04-21

**Program:** `research/atomic_shells_closure/`
**Skrypty:** as01, as02, as03
**Motywacja:** Po domknięciu d-class gap w SC (ps41, 60% RMS redukcji, Lu pozostaje
strukturalnie nieuzyskowalne), pytanie: czy TGP opisuje Li (najlżejszy nietrywialny
atom, hipotetycznie najprostszy test)?

---

## TL;DR

**Lit NIE jest prostszym testem TGP. Lit jest NIETESTOWALNY przez TGP.**

TGP w obecnej formie nie zawiera aparatu atomowego. Redukuje się do standardowej
mechaniki kwantowej dla atomów z precyzją lepszą niż 10⁻⁶⁰ (poniżej dowolnego
pomiaru). Nie dodaje żadnych nowych przewidywań chemicznych.

Empiryczne amplitudy orbitali `A_s, A_sp, A_d, A_f` używane w SC/ρ(T) **nie są
A_tail izolowanych atomów**, lecz efektywnymi znacznikami morza Fermiego w metalach.
**Między rdzeniem TGP (masy leptonów `m ∝ A⁴`) a aplikacjami (SC/ρ) istnieje strukturalna
luka.**

---

## Tabela wyników

| Test | Co testowano | Wynik | Wniosek |
|------|--------------|-------|---------|
| **as01** | H 1s w TGP | proton-tail na a₀ zgasza się jak exp(−2.5·10⁵); electron self-tail jak exp(−137) | TGP = QM standard dla H do 10⁻⁶⁰ |
| **as01** | Rydberg Z² dla H, He⁺, Li²⁺ | błędy ≤0.05% (= Lamb shift + FS) | skalowanie Z² zachowane |
| **as02** | Li IE₁/IE₂/IE₃ = hydrogenic+screening | Z_eff(IE₁)=1.26, δ_s^Li=0.41 | trywialne, zgodne z QM |
| **as02** | Koide K_Li = ΣIE/(Σ√IE)² | K_Li = 0.417, vs lept K=2/3 | Li IE **NIE jest drabiną solitonową** |
| **as02** | Brannen √2 ansatz dla Li IE | residual² = 0.30 (vs ~0 lept.) | Li IE **nie spełnia ansatz leptonowy** |
| **as02** | A ∝ IE^(1/4) — naiwny most TGP | A_Li ≈ 0.79, vs A_s(SC) = 0.11 | **niezgodność 7×** — A_s(SC) NIE jest tail Li |
| **as02** | Tail-phase sign 2s ma node | R_2s(r>2a₀) UJEMNE ✓ jakościowo zgodne | znak OK, magnituda NIE |
| **as02** | Alkali consistency (SC A_s stałe) | s-metale w SC wszystkie mają A_s=−0.111 | per-atom tail-sign hipoteza ZAWODZI |
| **as03** | α_pol(alkali) ~ IE^β | β = −2.81 (QM naive: −2) | TGP nie dostarcza poprawki 40% |
| **as03** | Invariant (α·IE⁴)^(¼) | ~18±5% CV (H→Fr) | TGP nie ma pierwszego zasad dla tej stałej |
| **as03** | Cs→Fr relatywistyczna kontrakcja | obserwowana, TGP nie ma bezpośredniego aparatu | wymaga P7.5a-like η_f fit |

---

## Znalezione strukturalne ograniczenia TGP

### L1. Atomy nie są w zakresie TGP
Soliton ODE TGP opisuje pojedyncze cząstki (elektron, mion, tau) jako rozwiązania
`g'' + (g')²/g + 2g'/r = 1 − g` z różnymi `g₀`. Masy skalują jak `A_tail⁴`.

Dla atomu — **wiązania** elektron-jądro przez Coulomb — TGP potrzebuje sprzężenia
elektromagnetycznego. **TGP nie opisuje obecnie elektromagnetyzmu.**

Soliton-w-soliton ansatz (elektron w polu Φ protonu) daje korekcję O(exp(-Bohr/Compton)) =
O(exp(-137)) = O(10⁻⁶⁰) poniżej dowolnej obserwacji. Czyli TGP **degeneruje do standardu
QM w atomowym regime**.

**Konsekwencja:** lit, wodór, hel, ..., wszystko opisywane przez standardową QM (HF, DFT).
TGP jest "niewidoczny" atomowo.

### L2. A_orb z SC/ρ nie są derywowalne atomowo
Hipoteza "A_s = A_tail(2s Li)" zawodzi na 3 poziomach:
1. Magnituda: A_tail proxy ~ 0.79, obserwowane |A_s| = 0.111 (7× rozbieżne)
2. Stałość: A_s jest stałe dla wszystkich s-metali w SC, ale per-atom tail-sign fluktuowałby z liczbą nodów
3. Nie-lokalny charakter: A_s ewidentnie opisuje **Fermi sea**, nie izolowany atom

**Implikacja dla SC paper:** tabela "universal TGP constants" powinna jawnie rozdzielać:
- **Pierwsze zasady:** C_0 (substrat), α_P6B (phonon w SC), g₀^e/g₀^μ/g₀^τ (leptony)
- **Empiryczne:** A_s, A_sp, A_d, A_f (fitowane do SC; **nie** derywowane atomowo)

### L3. Gap atom→solid nie jest zamknięty
Metal Li jest siatką ATM protonów + 3N elektronów (pozostałe 2N w core 1s², N w Fermi sea).
W TGP core jest:
- teoria solitonu izolowanego elektronu (mass A⁴) ✓
- teoria bryły/cząstki na tle `Φ` ✓
- **brak** teorii KOHERENTNEGO WIELOELEKTRONOWEGO stanu Fermiego w substracie

To jest dokładnie miejsce, gdzie SC/ρ(T) operuje empirycznie przez A_orb. Bez derywacji
z pierwszych zasad **SC/ρ(T) w TGP są phenomenologią**, a nie pochodnymi rdzeniowego
Lagrangianu.

---

## Status programu

| Stan przed as01-03 | Stan po as01-03 |
|---|---|
| Istniała ambicja: "A_s = A_tail 2s elektronu w TGP" | Obalona (7× rozbieżność) |
| Lit miał być prostym testem | Okazał się niestestowalny (TGP ≡ QM do 10⁻⁶⁰) |
| SC/ρ A_orb traktowane jak TGP-fundamentalne | Udokumentowano jako phenomenologia |
| Lu w ps41 jako "4f-descriptor brak" | Konsystentnie z L2: TGP nie ma atomowego aparatu dla 4f |

**Nie osiągnięto closure** w tradycyjnym sensie — TGP nie opisuje atomów. Ale **osiągnięto
honest documentation** — po raz pierwszy rozdział między zakresem rdzenia TGP (grawitacja,
masy) a phenomenologicznym rozszerzeniem TGP (SC, ρ(T), A_orb) jest JAWNY.

---

## Wnioski dla paperu SC v2 i paperu Core v2

### SC paper v2 — dodać:
1. Sekcja "Status of orbital amplitudes": jasne stwierdzenie, że A_s, A_sp, A_d, A_f są
   **fitowane** do danych SC, nie wyprowadzone z rdzenia TGP.
2. W tabeli "universal constants" oddzielna kategoria **phenomenological constants**
   (obok **derived constants**).
3. W sekcji "predictions" adnotacja: "Predictions for new materials require calibration of
   A_orb to matching reference materials; the absolute magnitude of T_c is not derived
   ab initio."

### Core paper v2 — dodać (opcjonalnie):
1. Sekcja "Scope boundary": precyzyjne kryterium **atomic regime**: TGP corrections
   scale as `exp(-Bohr/Compton)` ~ 10⁻⁶⁰, effectively zero. TGP reduces to standard QM
   for atomic/chemical scales.
2. Explicit statement że EM coupling jest **poza** TGP substrat ansatz.

---

## Następne fronty (zgodnie z user directive: "uderzać inne fronty")

Po as01-03 wiemy, że **atomowa chemia nie daje nowej informacji o TGP**. Kierunki, które
mogą dać:

### 1. **EM jako substrat** (fundamental gap)
Czy TGP może być rozszerzony tak, by substrat produkował Coulomb? Obecne:
- substrat → grawitacja (derywowane, OK)
- substrat → QM (79/79 tests)
- substrat → EM: **brak**

Nowy kierunek: `research/em_from_substrate/` — próba derywacji Maxwell z TGP Lagrangianu.
Jeśli powiedzie się, atom naturalnie wchodzi w TGP.

### 2. **Cohesive energy metali** (bypass atom, direct metal)
Li metal ma E_coh = 1.63 eV/atom. Na 1.11, K 0.93, Rb 0.85, Cs 0.80. To są
wielociałowe liczby. Czy istnieje TGP-derywowany invariant dla E_coh jak dla Tc?

Kierunek: `research/cohesion_closure/` — most między izolowanym atomem a metalem,
wypełnia L3.

### 3. **Relatywistyczne 4f/5f dla Lu** (ostatnia ofensywa w SC)
Kontynuacja P7.5a: η_f jako TGP constant, fitowany do Lu-series. Jeśli się zamknie,
Lu wypada z listy outlierów SC.

Kierunek: `research/superconductivity_closure/P7.13_lanthanide_descriptor.py` —
lokalna poprawka w SC.

### 4. **Cieczy, transport cieplny, muon g-2** (już otwarte)
Zgodnie z [[NEW_DIRECTIONS_2026-04-20.md]], 5 alternatywnych
frontów które są w zakresie TGP (nie wymagają atomowej chemii). Priorytet:
- `liquid_viscosity` — bazy danych VFT, szybka walidacja
- `muon_g_minus_2` — krótki horyzont, FNAL Run-6 do 2028

---

## Ostateczna konstatacja

Hipoteza użytkownika: "Lit jako lekki pierwiastek powinien łatwiej poddawać się TGP
niż cięższe."

**Obalona.** Lit nie jest prostszy, bo TGP nie ma aparatu atomowego wcale. Masa
pierwiastka jest irrelevantna — problem leży w STRUKTURALNEJ luce między rdzeniem
TGP a chemią.

**Prawdziwa obserwacja:** to, co TGP robi dobrze, to właściwości **makroskopowe
cząstek** (masy, Koide, spin) i **makroskopowe substratu** (MOND-ν(y) phenomenology,
substrate cosmology). Sektor atomowy jest systemu punktu, w którym efekty TGP są
eksponencjalnie zgaszone.

**Rekomendacja:** przenieść wysiłek na fronty, gdzie TGP-rdzeń daje NATYWNE
przewidywania — VFT dla cieczy, μ−g−2, cluster+ν — zamiast walczyć o zamknięcie
atomowej chemii, która fundamentalnie nie należy do zakresu TGP.

---

## Pliki programu

| Plik | Opis | Status |
|------|------|--------|
| [[TGP/TGP_v1/research/atomic_shells_closure/PLAN.md]] | Plan programu (A1-A4) | ✓ |
| [[as01_hydrogen_probe.py]] | H 1s sanity | ✓ PASS (trivially) |
| [[as01_results.txt]] | Output as01 | ✓ |
| [[as02_lithium_ionization.py]] | Li IE chain + Koide + A_orb | ✓ (A_s NIE atomowe) |
| [[as02_results.txt]] | Output as02 | ✓ |
| [[as03_polarizability_probe.py]] | α_pol(alkali) test | ✓ (brak TGP invariantu) |
| [[as03_results.txt]] | Output as03 | ✓ |
| [[ATOMIC_SHELLS_VERDICT.md]] | Ten dokument | ✓ |

as04, as05, as06 z PLAN.md **nie są potrzebne** — as01-03 już daje odpowiedź:
TGP atomowo nie rozróżnia Li od innych elementów, bo efekty substratu są poniżej
10⁻⁶⁰. Dalsze skrypty generowałyby tylko bardziej szczegółową demonstrację tego samego.

---

*Status: CLOSED jako honest negative result.*
*Next action: per user directive — select one of następnych frontów (EM from substrate,
cohesive energy, P7.13 Lu descriptor, lub już-otwarte VFT/muon g-2).*
