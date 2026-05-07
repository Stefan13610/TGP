---
title: "Recenzja zewnętrzna TGP_v1 — 5 strukturalnych długów (2026-05-06)"
date: 2026-05-06
parent: "[[README.md]]"
type: external-review
tgp_owner: audyt
tags:
  - audit
  - external-review
  - peer-review
  - structural-debt
  - foundations
related:
  - "[[README.md]]"
  - "[[SUMMARY_2026-05-04.md]]"
  - "[[CLOSURE_SUMMARY_2026-05-06.md]]"
  - "[[PRIORITY_MATRIX.md]]"
  - "[[L01_rho_operational/README.md]]"
  - "[[../TGP_FOUNDATIONS.md]]"
tgp_status:
  folder_status: audit
  level: mixed
  kind: external-review
  core_compatibility: review-only
  last_reviewed_against_core: 2026-05-06
  may_edit_core: false
  exports_findings: true
  has_needs_file: false
  has_findings_file: false
  open_bridges: []
  depends_on: []
  impacts: []
  source_of_status: []
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-06
revisions:
  - v1 (2026-05-06): Pierwsza wersja, 5 punktów EXT-1..5
  - v2 (2026-05-06): EXT-1 przeformułowany po korekcie autora dot. ax:c–ax:G
    (varying constants); status zmieniony z "wysokie ryzyko falsyfikacji"
    na "otwarte ryzyko, wymaga jawnego FRW calc"
---

# Recenzja zewnętrzna TGP_v1 — 5 strukturalnych długów

> ## ⚠ UPDATE 2026-05-06 v2 — korekta EXT-1
>
> **Autor wskazał** (sesja 2026-05-06), że recenzent w wersji v1 pominął
> rolę aksjomatów `ax:c`–`ax:G` (`core/sek04_stale/sek04_stale.tex` lin. 29,
> 250–254): w TGP **c₀, ℏ₀, G₀ nie są stałymi fundamentalnymi**, lecz
> wartościami referencyjnymi w obecnej epoce; same c, ℏ, G są funkcjami
> pola Φ:
>
> ```
> c(Φ) = c₀ √(Φ₀/Φ)        ℏ(Φ) = ℏ₀ √(Φ₀/Φ)        G(Φ) = G₀ Φ₀/Φ
> ```
>
> To otwiera ścieżkę domknięcia EXT-1, której v1 nie rozważał: **H_TGP(z)
> w erze radiacyjnej nie musi być H_GR(z)** — pochodzi z Φ-EOM w FRW
> z varying constants, gdzie foton jest "pasażerem geometrii" (porusza
> się po g_eff[Φ]), nie źródłem ρ przez L_mat. EXT-1 v2 reflektuje tę
> dwustronność: **T^μ_μ_EM = 0 trzyma się strukturalnie** (konformalna
> niezmienniczość Lagrange'a EM w 4D), ale **konsekwencja kosmologiczna
> "brak fenomenologii radiacji" jest podejrzana**, bo Φ-dynamika +
> varying c,ℏ,G ma własną strukturę H_TGP(z).
>
> Pełny tekst EXT-1 v2 poniżej zastępuje v1. Pozostałe punkty EXT-2..5
> bez zmian.

## Cel

Niniejsza recenzja jest **zewnętrznym, niezależnym przeglądem** TGP_v1
z perspektywy fizyka teoretycznego, który nie pracował w cyklu rozwojowym
TGP. W odróżnieniu od audytów S/L/D/M (`README.md` § Klasy), które
inwentaryzują *wewnętrzne* sprzeczności i luki, ten dokument identyfikuje
**5 długów strukturalnych**, które rozstrzygną — w opinii recenzenta —
czy TGP zostanie z czasem rozpoznana jako *teoria wyprowadzona*, czy
jako *konstrukcja kompatybilna z obserwacjami*. Różnica między tymi
dwoma statusami jest istotna: pierwszy daje TGP siłę falsyfikującą,
drugi czyni ją jednym z wielu modeli zgodnych z danymi w słabym polu.

## Podstawa recenzji

Materiał przeczytany w sesji 2026-05-06:

- `TGP_FOUNDATIONS.md` (binding ontologia, §1–§10)
- `core/sek01_ontologia/sek01_ontologia.tex` (lin. 1–1435)
- `core/sek02_pole/sek02_pole.tex` (lin. 1–262, full)
- `core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex`
  (lin. 1–300; G.0 closure ADDENDUM 2026-05-02; L01 formal definition
  2026-05-04 lin. 148–183)
- `core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex` (lin. 1–700)
- `core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex` (header)
- `core/sek05_ciemna_energia/sek05_ciemna_energia.tex` (lin. 1–300)
- `core/sek08_formalizm/sek08_formalizm.tex` (selektywnie:
  H_Γ lin. 736+; T_μν^(Φ) lin. 7989–8003; SEC violation lin. 8006–8047;
  G=κT identity lin. 10607–10615; prop:ghost-free-fundamental
  lin. 2588–2604)
- `audyt/README.md`, `audyt/PRIORITY_MATRIX.md`,
  `audyt/L01_rho_operational/README.md`, `audyt/L03_K_phi_stability/README.md`

## Zakres recenzji

Recenzja **zakłada poprawność** wewnętrznych zamknięć audytowych
(S01–S06 CLOSED, L01–L04 EXECUTED, L03 EXECUTED 2026-05-06, D01
EXECUTED 2026-05-06). Punkty poniżej dotyczą **długów wykraczających
poza obecny rejestr S/L/D/M** — problemów, które są skutkami strukturalnymi
samej *konstrukcji* TGP, a nie chaotycznego rytmu rozwoju.

---

## Punkt 1 (v2) — Kosmologia radiacyjna w TGP z varying c, ℏ, G

### Diagnoza

L01 formal definition (cykl 2026-05-04,
`research/op-L01-rho-stress-energy-bridge-2026-05-04/`) zamyka
operacyjną definicję `ρ` przez tożsamość:

```
ρ(x) ≡ −T^μ_μ(x) / c_0²
```

(`sek08a_akcja_zunifikowana.tex` lin. 148–183, mapping SM 5 sektorów).
W połączeniu z aksjomatami stałych (`sek04_stale.tex` lin. 29, 250–254):

```
c(Φ) = c₀ √(Φ₀/Φ)     ℏ(Φ) = ℏ₀ √(Φ₀/Φ)     G(Φ) = G₀ Φ₀/Φ
```

bezpośrednie konsekwencje dla sektorów SM:

| Sektor SM | T^μ_μ | ρ_TGP | Status w erze radiacyjnej |
|-----------|-------|-------|----------------------------|
| Fermion masywny (Dirac) | m·\|ψ\|² | m·\|ψ\|²/c_0² ≥ 0 | mała ilościowo (n_b/n_γ ≈ 10⁻⁹) |
| Foton (EM) | 0 (konformalny) | **0** | strukturalnie zerowy |
| Yang-Mills klasyczny | 0 | **0** | strukturalnie zerowy |
| QCD kwantowy (kondensat) | ~Λ_QCD⁴ | Λ_QCD⁴/c_0² (anomalia śladu) | znacząca, ale stała w czasie |
| Promieniowanie (p = ρ_e/3) | 0 | **0** | strukturalnie zerowy |
| Pył (p = 0) | ρ_e | ρ_e/c_0² | mała w erze radiacyjnej |
| Ciemna energia (p = −ρ_e) | 4ρ_e | 4ρ_e/c_0² | zaniedbywalna w wczesnej epoce |

**Strukturalnie:** promieniowanie nie generuje Φ przez L_mat. T^μ_μ_EM = 0
jest **twierdzeniem matematycznym** (Weyl-niezmienniczość Lagrange'a EM
−¼F_μν F^μν √(−g) w 4D), niezależnym od tego, czy w mianowniku stoi
c₀ czy c(Φ). To trzyma się.

**Fenomenologicznie:** ale H_TGP(z) w erze radiacyjnej **nie musi być
równe H_GR(z)** — pochodzi z Φ-EOM w FRW, gdzie dynamika Φ +
varying c(Φ), ℏ(Φ), G(Φ) tworzą własną strukturę kinematyczną
ekspansji. Foton w TGP jest **pasażerem geometrii** (porusza się po
g_eff[Φ]), ale geometria ewoluuje przez Φ-dynamikę, więc foton
*doświadcza* redshiftu, decouplingu, BBN-relevantnej kinematyki —
mimo że nie jest *źródłem* Φ.

### Pliki dotknięte

- `sek08a_akcja_zunifikowana.tex` lin. 148–183 (L01 formal definition,
  comment block z mapping)
- `sek04_stale/sek04_stale.tex` lin. 27–82 (ax:c jako pole),
  lin. 178–238 (prop:c-from-metric — wyprowadzenie z metryki),
  lin. 250–254 (ax:c–ax:G hierarchia), lin. 508–556 (tabela skalowań)
- `research/op-L01-rho-stress-energy-bridge-2026-05-04/`
  (formal_definition.md, SM_sector_mapping.md, photon_treatment.md)
- `audyt/L01_rho_operational/POST_ACTION_UPDATE_2026-05-04.md`
- `core/sek08a` lin. 119–127 (eq:L-mat-unified)
- **BRAKUJĄCY plik:** `research/op-FRW-radiation-era-varying-c/`
  (cykl jeszcze nie wykonany — patrz "Ścieżki domknięcia")

### Problemy jakie rodzi

1. **Brak istniejącego rachunku.** Nie zlokalizowałem w `research/`
   jawnego rachunku Φ-EOM w FRW dla z > z_rec z porównaniem do BBN
   i CMB. Cykl `research/op-newton-momentum/M9.x` zawiera FRW
   *zasadniczo*, ale przede wszystkim w kontekście weak-field PPN
   i obecnej ekspansji, nie ery radiacyjnej. To **otwarte zadanie**,
   nie zamknięty problem.

2. **Reżim ekstremalnie wczesny** (z > 10¹⁰): jeśli Φ_then << Φ₀, to
   c(Φ) >> c₀, ℏ >> ℏ₀, G >> G₀. M9.1'' (zaprojektowana dla weak field
   wokół ψ = 1) może łamać założenia perturbacyjne — wyższe rzędy w
   1/ψ mogą dominować. **Czy M9.1'' jest valid w erze radiacyjnej?**
   Otwarte. Może potrzebna oddzielna analiza w reżimie ψ << 1.

3. **Kinematyka rekombinacji** (z ≈ 1100, T ≈ 0.3 eV). Tempo
   decouplingu Thomson zależy od σ_T·n_e·c, gdzie:
   - σ_T = (8π/3) (e²/m_e c²)² zawiera c⁻⁴
   - energia Rydberga (Saha) zawiera ℏ⁻²
   - prędkość fotonu = c
   - w TGP wszystkie skalują z Φ jednocześnie
   
   Czy kombinacje się **znoszą** (zachowując standardową kinematykę),
   czy **łamią** (modyfikując sound horizon i pierwszy peak CMB)?
   Bez explicit rachunku nie wiadomo.

4. **Konkretne pytania liczbowe.** Aby TGP odzyskała fenomenologię BBN:
   - H_TGP(z=10⁹) ≈ H_GR(z=10⁹) z dokładnością ~5% (BBN tolerance dla ⁴He);
   - położenie pierwszego peaka CMB l_1 = 220.0 ± 0.5;
   - sound horizon r_s(z_drag) zgodny z Plancka 147 Mpc;
   - N_eff = 3.046 (relativistic species).
   
   Spójność wszystkich *jednocześnie* w modelu z 1–2 swobodnymi
   parametrami (γ, Φ₀ kalibrowane do innych obserwabli) jest *bardzo
   nietrywialna*. To nie jest "raczej wyjdzie".

5. **Dwustronność diagnozy.** Jeśli explicit FRW calc daje **zgodność**
   — to jest *fenomenalny* wynik: BBN+CMB z 1–2 wolnymi parametrami
   zamiast standardowych ~6 ΛCDM. Jeśli daje **niezgodność** — to
   silny sygnał, że L_mat wymaga rozszerzenia (powrót do ścieżki D
   poniżej, naruszającej S04). Obie możliwości są *równie ważne*
   i obie wymagają rachunku.

6. **Open NEEDS w L01** — `research/op-L01-rho-stress-energy-bridge-2026-05-04/`
   sygnalizuje N1 (quantum trace anomaly EM) i N2 (QCD vacuum)
   jako otwarte. W połączeniu z varying constants te otwarte punkty
   są tym, co rozstrzyga, czy fotony naprawdę są "tylko pasażerami",
   czy mają drugorzędny wkład do ρ przez anomalię.

### Potencjalne ścieżki domknięcia

**Ścieżka A — explicit Φ-EOM w FRW dla ery radiacyjnej (cykl
op-FRW-radiation-era-varying-c).** Wziąć Φ-EOM z `eq:field-eq-reproduced`
+ M9.1'' metric + ax:c–ax:G + ρ_matter (mała) + ρ_QCD anomalia (stała).
Policzyć H_TGP(z) dla z ∈ [10³, 10¹⁰]. Porównać z H_GR(z). Jeśli
różnica < 5% w obszarze BBN i < 0.5% w obszarze rekombinacji — TGP
odzyskuje standardową fenomenologię *bez* ρ_rad jako źródła. Jeśli nie
— jasny sygnał o niezgodności i konieczności pivotu L_mat.

**Ścieżka B — BBN abundance calc (op-BBN-TGP).** Z H_TGP(z) ze
ścieżki A wziąć Boltzmann codes (np. AlterBBN port pod TGP), policzyć
abundance ⁴He, D/H, ³He/H, ⁷Li/H. Porównać z PDG (Y_p = 0.245 ± 0.003,
D/H = 2.55 ± 0.05 · 10⁻⁵). Tu dochodzą efekty varying ℏ przez energie
wiązań nuklearnych — to nietrywialne, bo ℏ wpływa na Coulomb
barrier i kinetyki reakcji. Pełny rachunek wymaga rewizji
nuklearnych cross-sections w funkcji Φ.

**Ścieżka C — CMB peaks (op-CMB-TGP).** CAMB/CLASS port pod TGP
cosmology z varying c, ℏ, G. Position pierwszego peaka (l ≈ 220) jest
fenomenalnie wrażliwy na sound horizon r_s(z_drag), który zależy od
H(z), c_s(z), c(Φ_drag). Jeśli TGP daje l_1 zgodne z Plancka — DERIVED.
Jeśli nie — sygnał konieczności modyfikacji.

**Ścieżka D — alternatywne sprzężenie L_mat dla pól bezmasowych
(rezerwowa, jeśli A/B/C zawiodą).** L_mat = −(q/Φ_0) φ ρ + L_rad(φ, F_μν)
z dodatkowym sprzężeniem dilatonowym dla pól cechowania. To narusza
ax:metric-coupling (S04 zamknięty 2026-05-04 przez B9 MICROSCOPE 6/6
PASS), więc droga ostatnia. Ale jeśli ścieżki A/B/C niezgodne, jedyna
opcja zachowania predykcyjności TGP w erze radiacyjnej.

**Ścieżka E — przyznanie zakresu (najszczersza, jeśli rachunki
wskażą trudność).** TGP jako "GR w limicie post-recombination",
z explicit zaznaczeniem, że dla z > z_rec teoria nie pretenduje do
zastąpienia standardowej kosmologii. To drastycznie obniża rangę
TGP, ale jest *honest*. Powinno być rozważane tylko po nieudaniu się
ścieżek A/B/C.

### Rekomendowany priorytet

**P1 — krytyczny, otwarte ryzyko.** L01 jest oznaczony jako EXECUTED
2026-05-04 w `PRIORITY_MATRIX.md`, ale "executed" oznacza tu *formal
definition*, nie *fenomenologiczne sprawdzenie w erze radiacyjnej*.
Bez ścieżek A+B+C explicit policzonych, **nie wiadomo**, czy TGP
ratuje fenomenologię BBN/CMB przez varying constants, czy nie. Status
EXT-1: **otwarte**, nie *prawie na pewno przegrane* (jak w v1).

**Subiektywna ocena recenzenta:** prawdopodobieństwo, że ścieżka A da
H_TGP(z) zgodne z BBN/CMB w 5% / 0.5% tolerance, wynosi ~55–65% —
niepewność, nie a priori przegrane. Niepewność wynika z faktu, że
*z jednej strony* varying constants dają TGP pełną swobodę
re-strukturyzacji H(z), *z drugiej strony* jednoczesna spójność BBN
i CMB w 1–2 wolnych parametrach jest nietrywialnym constraint na
strukturę Φ-EOM.

### Powiązanie z istniejącym audytem

Rozszerza [[L01_rho_operational/README.md]] poza obecny zakres
(operacyjna definicja → kosmologia radiation-era z varying c, ℏ, G).
NEEDS N1, N2 w `research/op-L01-rho-stress-energy-bridge-2026-05-04/`
powinny zostać przepromowane do P1 i połączone z nowym cyklem
`op-FRW-radiation-era-varying-c/`.

### Uwaga merytoryczna

Punkt 1 v1 błędnie zakładał, że "brak ρ_rad jako źródła Φ → brak fazy
radiation-dominated w TGP". To było zaniedbanie roli `ax:c–ax:G` —
varying c, ℏ, G dają TGP możliwość re-strukturyzacji H(z) niezależnej
od standardowego mechanizmu ρ_rad ~ a⁻⁴. Punkt 1 v2 reflektuje tę
korektę: T^μ_μ_EM = 0 trzyma się strukturalnie, ale fenomenologiczne
konsekwencje wymagają jawnego rachunku FRW. Status: otwarte z wysokim
priorytetem.

---

## Punkt 2 — Wyprowadzenie warunku zerowej sumy z H_Γ

### Diagnoza

`sek05_ciemna_energia.tex` lin. 240–293 (`prop:Lambda-positive`)
opiera **całość mechanizmu Λ_eff** na warunku zerowej sumy:

```
∫_Σ φ √h d³x = 0     (eq:zero-sum-phi, lin. 250)
```

który zwany jest "zasadą zerowej sumy" (`ax:zero` w sek01).
Bez tego warunku ⟨φ²_min⟩ niekoniecznie > 0, a Λ_eff = (8πG_0/c_0⁴)
⟨U(φ_min)⟩ niekoniecznie > 0. Cała sek05 (ciemna energia z fluktuacji)
zawisa na tym jednym wymaganiu.

Warunek jest **globalny** — więz na całej hypersurface przestrzennej Σ.
Nie wynika z lokalnej dynamiki (Φ-EOM), nie wynika z H_Γ jawnie.

### Pliki dotknięte

- `sek05_ciemna_energia.tex` lin. 146–161 (subsection
  "Związek z zasadą zerowej sumy"), 240–293 (prop:Lambda-positive)
- `sek01_ontologia.tex` aksjomat ax:zero (numer linii do zlokalizowania
  w pełnym pliku 1435 lin.)
- `sek08_formalizm.tex` lin. 736+ (H_Γ definition v2 GL-bond)
- `axioms/substrat/` (substrat directory — referencje do roli Z₂)

### Problemy jakie rodzi

1. **Lokalność vs globalność.** Globalny więz bez Lagrange'a multiplier
   z dynamiki to klasyczna czerwona flaga w fizyce teoretycznej. Albo
   jest konsekwencją gauge/redundancji (jak constraint Hamiltoniana
   w GR z dyfeomorfizmu), albo wprowadza nielokalność w postaci
   "system jako całość zna globalny rozkład materii".

2. **Operacyjna nieobserwowalność.** Warunek dotyczy całej hypersurface
   Σ jednocześnie. Lokalny eksperymentator nie może go zweryfikować ani
   złamać. To czyni go **nie-falsyfikowalnym** — co jest niezgodne
   z duchem TGP (deklarowanym w `TGP_FOUNDATIONS.md` § 5.3:
   "GR jako numeryczny analog, nie izomorfizm; otwiera drzwi
   falsyfikacji").

3. **Mechaniczna ad-hoc-owość.** Wprowadzenie ax:zero jako *aksjomat*
   przyznaje, że jest on potrzebny, ale nie wynika z głębszej
   struktury. W ten sposób TGP ma 3 fundamenty:
   - jedno pole Φ z Z₂ (§1 FOUNDATIONS),
   - hamiltonian H_Γ z GL-bond,
   - **warunek zerowej sumy ∫ φ √h = 0**.
   Trzeci wygląda na *montowanie* zamiast *wyprowadzenia*.

4. **Kompatybilność z lokalnymi rozwiązaniami.** Jeśli rozważam
   izolowaną cząstkę w TGP (radialny kink Φ(r)), to lokalnie φ > 0
   (rośnie wokół źródła). ∫ φ √h > 0 dla każdej skończonej kuli. Więz
   ∫_Σ φ = 0 wymaga, by **gdzieś indziej** φ był ujemny — ale
   `sek02_pole` lin. 116–117 explicit zakazuje φ < 0 ("Φ(r) > 0 wszędzie,
   przestrzeń jest generowana, nie 'ujemna'"). Sprzeczność?
   Albo "deficyt φ < 0" oznacza coś innego niż dosłowna ujemność —
   ale wtedy potrzebne jest doprecyzowanie operacyjne.

5. **Cosmological constant problem zwija się tylko częściowo.**
   Λ_eff ~ γ/12 ~ Φ_0 H_0²/(12 c_0²). Wymaga γ w skali H_0² — ale
   *dlaczego* γ jest w skali H_0², a nie M_Pl²? Standardowy CC problem
   (10¹²⁰ rozbieżność QFT vacuum vs obs) został przesunięty z "dlaczego
   Λ jest takie małe" na "dlaczego γ jest takie małe". Bez wyprowadzenia
   γ z UV, problem wraca tylnymi drzwiami.

### Potencjalne ścieżki domknięcia

**Ścieżka A — derywacja zerowej sumy z Z₂-symetrii substratu.**
Hipoteza: ∫_Σ φ √h = 0 może być konsekwencją tożsamości typu Bianchi
dla pola φ na grafie Γ. Konkretnie: jeśli H_Γ = −J Σ (φ_i φ_j)² ma
chiralną Z₂ (`sek01_ontologia.tex:83`), to operator parzystości
P: φ_i ↔ −φ_i pozostawia H_Γ invariantnym. Suma ∫ φ √h nad
hypersurface zachowującą P powinna być zerem z konstrukcji.
**Cykl proponowany:** `op-zero-sum-derivation/` z explicit derywacją
∫ φ = 0 jako tożsamość Z₂ (nie aksjomat).

**Ścieżka B — Lagrange'a multiplier w akcji zunifikowanej.**
Dodać do S_TGP człon λ ∫ φ √h. Wariacja po λ daje warunek;
po φ daje EOM z λ jako efektywnym Λ. To by *zlikwidowało* ax:zero
jako aksjomat, redukując go do dynamicznego więzu z explicit
mnożnikiem. Cena: λ jest dodatkowym parametrem, którego skala wymaga
wyjaśnienia (ale to znów jest fine-tuning γ ~ H_0², więc nic nie tracimy
w stosunku do obecnego stanu).

**Ścieżka C — Wprowadzenie φ_eff = φ − ⟨φ⟩_Σ.** Reinterpretacja:
"zerowa suma" jest *definicją* φ_eff jako odchylenia od średniej.
Wówczas ∫ φ_eff = 0 jest tożsamością. Λ_eff przestaje wynikać
z warunku, a zaczyna wynikać z ⟨φ⟩_Σ — co jest mierzalne (tło
kosmologiczne). Wymaga przekształcenia całej sek05.

**Ścieżka D — domknięcie przez nielokalność.** Przyznanie, że TGP jest
nielokalna na skali kosmologicznej (nie sprzeczne z relativity jeśli
nielokalność jest *spacelike* na hypersurface), z explicit dyskusją
horyzontów i kompatybilności z FRW. Cykl `op-nonlocal-foundations/`.

### Rekomendowany priorytet

**P2 — wysoki.** Bez ścieżki A lub B, sek05 (główna predykcja
kosmologiczna TGP) opiera się na aksjomacie, którego status nie
jest ani derived, ani falsyfikowalny.

### Powiązanie z istniejącym audytem

Nie ma odpowiednika w obecnym S/L/D/M. Proponuję klasę **L07** —
"warunek zerowej sumy: aksjomat vs derywacja".

---

## Punkt 3 — Niezależne wyprowadzenie hiperbolicznej metryki M9.1''

### Diagnoza

Forma metryki efektywnej była iterowana:

| Iteracja | Forma g_tt | Status | Powód odrzucenia |
|----------|-----------|--------|------------------|
| Forma I (potęgowa) | −c²/ψ | DEPRECATED 2026-04-25 | β_PPN = 4 vs obs 1, 3·10⁴σ |
| Forma II (eksponencjalna) | −c² e^(−2U) | DEPRECATED | równoważna do O(U), różnica na 2PN+ |
| Forma III (M9.1'' hiperboliczna) | −c²(4−3ψ)/ψ | CANONICAL 2026-04-25 | γ_PPN = β_PPN = 1 *exact* |

(`sek08a_akcja_zunifikowana.tex` lin. 196–251, M9.1''
canonical replacement).

**Problem:** "exact by construction" jest *podejrzanie zręcznym
wyborem*. Master formula (P23 sympy LOCK 5/5 PASS) reprodukuje
γ = β = 1, ale **nie ma argumentu, który *przed* zobaczeniem PPN
forsowałby hiperboliczną formę**. To znaczy, że hipoteza została
*dobrana tak*, by dawała poprawny PPN.

### Pliki dotknięte

- `sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex`
  (G.0 preamble, eq:metric-M911-canonical)
- `sek08a_akcja_zunifikowana.tex` lin. 186–251, 240–290 (canonical
  M9.1'' replacement)
- `research/op-newton-momentum/M9.1*.md` (M9.1, M9.1', M9.1'' historia)
- `audyt/S01_metric_four_forms/` (cztery formy metryki — zamknięte
  przez G.0 closure, ale bez derywacji *z pierwszych zasad*)

### Problemy jakie rodzi

1. **Status epistemiczny.** TGP_v1 deklaruje hierarchię (W/E/P/H —
   wyprowadzenie/empiryczna/propozycja/hipoteza). Metryka M9.1'' ma
   status (E) w `TGP_FOUNDATIONS.md` § 2 — "metryka efektywna,
   emergencja w limicie GR, FRW, PPN". Ale **emergencja nie znaczy
   derywacja**: M9.1'' jest *postulatem*, którego konsekwencje są
   sprawdzalne. To bliżej (P) niż (E).

2. **Brak fundamentu z H_Γ.** `sek02:189–203` przyznaje wprost:
   "Dokładna postać tej relacji [g_eff = g_eff[Φ, ∂Φ, …]] jest
   **otwartym problemem TGP**." Hipoteza minimalna (lin. 210–221)
   ma konkurentów. M9.1'' pojawiła się jako **trzecia iteracja**
   po falsyfikacji form I i II — zatem jest *empirycznie wybrana*,
   a nie *strukturalnie wyprowadzona*.

3. **Kruchość liczbowa.** γ_PPN = β_PPN = 1 *exact* jest świetne, ale
   2PN współczynniki c_2 = −1, c_3 = +5/3, c_4 = −10/3 są specyficzne
   dla M9.1''. Oznacza to przewidywanie:
   |Δg_tt| = (5/6) U³ deviation od GR (testowalne LIGO 3G — sek08c
   M9.1'' P1). Ale *dlaczego ta wartość*? Bez fundamentu, jeśli LIGO 3G
   da inną wartość, M9.1'' upada — i TGP nie ma kolejnej iteracji
   *gotowej z konstrukcji*.

4. **Filozoficzny dług**. "Metryka emergentna" jest sztandarowym
   hasłem TGP. Ale jeśli emergentna metryka jest **postulowana
   funkcyjnie** (g_eff[ψ] = − c_0² (4−3ψ)/ψ jako ansatz), to różnica
   między TGP a teorią skalarno-tensorową z konkretnym ansatzem
   sprzężenia konformalnego jest *retoryczna*, nie strukturalna.
   `TGP_FOUNDATIONS.md` § 5.1 jawnie odrzuca "teorię skalarno-tensorową
   w stylu Brans-Dicke / Horndeski" — ale operacyjnie M9.1'' jest
   szczególnym przypadkiem skalarno-tensorowego ansatzu.

### Potencjalne ścieżki domknięcia

**Ścieżka A — derywacja M9.1'' z budżetu substratu.**
`sek08c_metryka_z_substratu` ma w nazwie "z substratu", co sugeruje
ambicję derywacji. Realizacja: explicit coarse-graining H_Γ na
graf regularny, identyfikacja efektywnej metryki przez Wilsonowską
RG flow. Wymaga cyklu `op-metric-from-substrate/` z wyprowadzeniem
w 4 krokach: H_Γ → continuum action → effective metric ansatz →
matching to M9.1''. Punkt 4 musi być *konsekwencją*, nie *wyborem*.

**Ścieżka B — wariacyjne kryterium minimalności.**
M9.1'' jako *jedyna* metryka spełniająca:
(i) γ_PPN = β_PPN = 1,
(ii) hiperboliczność horyzontu ψ = 4/3 jako natural BH cutoff,
(iii) konformalna struktura kompatybilna z α = 2 selection.
Jeśli (i)+(ii)+(iii) okażą się **wystarczającymi i konstruktywnymi**
warunkami, to M9.1'' staje się *wybrana z konstrukcji*, a nie
*dobrana*. Cykl `op-metric-uniqueness/`.

**Ścieżka C — przyznanie statusu (P) w hierarchii.** Najszczersze
rozwiązanie: zmienić status M9.1'' z (E) na (P) — *propozycja*,
nie *emergencja*. Przepisać `TGP_FOUNDATIONS.md` § 2, by hierarchia
była: (W) substrat, (W) Φ-EOM, **(P) M9.1''**, (P/H) materia.
Cena: TGP traci jeden z głównych argumentów marketingowych
("emergencja metryki"), zyskuje zaufanie naukowe.

**Ścieżka D — derywacja z budżetu Verlindego/Padmanabhana.**
W tradycji emergentnej grawitacji (sek08c referuje do Sakharov,
Verlinde) metryka wynika z entropic budgetu na horyzoncie. Spróbować
explicit oblicić, czy hiperboliczne g_tt = −(4−3ψ)/ψ wynika z
założenia: powierzchnia Hubble'a + entropia ∝ A_horizon / 4G + Z₂
substratowa. Cykl `op-entropic-metric/`.

### Rekomendowany priorytet

**P2 — wysoki.** Bez derywacji M9.1'' z fundamentu, TGP jest
operacyjnie skalarno-tensorową teorią ze szczególnym ansatzem.
To nie unieważnia jej, ale pozbawia jednego z głównych
roszczeń ontologicznych.

### Powiązanie z istniejącym audytem

[[S01_metric_four_forms]] zamknięty CLOSED-RESOLVED via G.0 2026-05-02.
Ale "RESOLVED" oznacza tu *spójność wewnętrzna* (4 formy → 1 forma),
nie *derywacja z pierwszych zasad*. Proponuję rozszerzenie: nowa
klasa **S07** — "metryka M9.1'' jako postulat vs derywacja".

---

## Punkt 4 — Domknięcie Phase 5/6 why_n3 (fermiony jako kinki)

### Diagnoza

`TGP_FOUNDATIONS.md` § 4 deklaruje warstwę 3c jako hipotezę:

> **3c — Kinki / defekty** (cząstka = radialny kink Φ + topologia
> chiralna). Hipoteza/roadmap alternatywnego opisu fermionów jako
> struktur w samym Φ; emergent Dirac propagator. **Strukturalny szkic
> 5-fazowy zamknięty 2026-05-01** w `research/why_n3/`. Otwarte (Phase 6+):
> X = e²/4 analytic derivation, A^(5−α) vs A²·g_0^(e²/2) reconciliation
> dla τ/e.

Phase 1–5 są zamknięte strukturalnie. **Phase 6+ pozostaje otwarta**,
i to ona zawiera kluczowe domknięcia analytyczne.

### Pliki dotknięte

- `research/why_n3/` (Phase 1–5 closed; Phase 6+ open)
- `TGP_FOUNDATIONS.md` § 4 (3c hipoteza)
- `core/sek08_formalizm.tex` lin. 9658+ (energie kinku K_n(α))
- `audyt/L05_mass_exponent_drift/` (k=4 LP-4 vs p=5−α R3 2026-05-01,
  P2 OPEN)

### Problemy jakie rodzi

1. **Spin-statistics theorem.** Poprawny opis fermionów wymaga, by
   stany dwucząstkowe były antysymetryczne (exchange phase = −1).
   Phase 3 `why_n3` używa "RP² + Berry phase π" — to znana ścieżka
   (Finkelstein-Rubinstein 1968), ale wymaga, by Berry phase
   *spójnie indukowała Lorentzowską strukturę spinową*. Bez explicit
   konstrukcji emergentnego propagatora Diraca z odpowiednimi
   antykomutacyjnymi własnościami, "kink jako fermion" pozostaje
   *roszczeniem strukturalnym*, nie *konstrukcją operacyjną*.

2. **Trzy generacje (e/μ/τ).** Phase 5 closed strukturalnie ma
   "uniwersalną formułę masy"
   m_obs(g_0, α) = c_M · A_tail² · g_0^(e²(1−α/4))
   (`research/why_n3/Phase5_*.md`). Dla α=2 daje g_0^(e²/2) ≈ g_0^3.69.
   Reprodukuje m_μ/m_e i m_τ/m_e z PDG <0.01%. **Ale e² w wykładniku
   jest empirycznym dopasowaniem** — bez derywacji wykładnika
   z głębszej struktury, formuła jest *spektakularnym numerologicznym
   sukcesem*, nie wyprowadzeniem. L05 (mass exponent drift) jasno to
   sygnalizuje (k=4 LP-4 vs p=5−α R3 — niezgodność wykładników między
   różnymi sektorami).

3. **Kwarki, neutrina, bozony cechowania.** TGP_v1 koncentruje się
   na leptonach. Kwarki (g_0 ∈ [0.817, 0.891]) są "uniwersalne" via
   ten sam ODE substratowy (`sek08b:529`), ale **explicit predykcje
   mas kwarków** nie są w PREDICTIONS_REGISTRY. Neutrina (Σm_ν) są
   w D01 jako anchor lock, nie jako derywacja. Bozony cechowania
   (W, Z, gluon) nie mają realizacji w warstwie 3c.

4. **Algebra Diraca.** Pola fermionowe w SM mają strukturę
   Cliffordowską {γ^μ, γ^ν} = 2g^μν. Z kinka skalarnego Φ z Z₂
   wyprowadzić tę algebrę jest nietrywialne. Skyrme model
   (1962) jest klasycznym precedensem (skyrmion = baryon w QCD-like
   chiralnej teorii), ale wymaga grupy chiralnej SU(N)_L × SU(N)_R,
   nie samej Z₂. TGP ma tylko Z₂ — to za mało dla pełnej algebry
   spinowej.

5. **Zamiast Diraca: czy emergentna SUSY?** Niektóre programy
   (np. Wen-Levin string-net) generują emergentne fermiony przez
   string-net condensation. To wymaga dyskretnej geometrii i
   deconfinement. TGP ma graf Γ, więc droga w zasadzie otwarta —
   ale niewystarczalność Z₂ pozostaje.

### Potencjalne ścieżki domknięcia

**Ścieżka A — explicit emergent Dirac propagator.**
Cykl `op-why_n3-Phase6-dirac/` z konstrukcją 2-cząstkowego stanu
fermionowego z dwóch kinków, weryfikacją antysymetrii pod exchange,
identyfikacją γ^μ z transportu Berry'ego po RP². Spin-1/2 ma być
*konsekwencją*, nie *postulatem*.

**Ścieżka B — derywacja e² w wykładniku masy.**
Cykl `op-why_n3-Phase6-mass-exponent/`. Hipoteza: e² w
g_0^(e²/2) wynika z ilości stopni swobody w ogonie oscylacyjnym
(2 polaryzacje × Berry phase 2π). Wymaga explicit obliczenia:
dlaczego e² (ładunek), a nie g (sprzężenie grawitacyjne)?

**Ścieżka C — przyznanie statusu hipotezy w prognozach.**
Najszczersze rozwiązanie: warstwa 3c pozostaje (H) *na zawsze*
w `TGP_FOUNDATIONS.md`, dopóki Phase 6+ nie zamknie się analitycznie.
TGP funkcjonuje jako teoria grawitacji + materia jako warstwa 3a/3b
(Dirac fields + g_eff coupling), bez roszczenia, że wszystko emerguje
z Φ. To **odróżnia** TGP od programu unifikacyjnego — czyniąc ją
**teorią grawitacji emergentnej**, jak Verlinde 2010.

**Ścieżka D — rozszerzenie symetrii substratu.**
Z₂ jest niewystarczająca dla algebry spinowej. Rozszerzenie do
SU(2) lub Spin(3) na poziomie substratu *narusza* §1 FOUNDATIONS
("dyskretna Z₂ chiralna na substracie") — to jest poważna decyzja
pivot substratu (zakazana bez rozmowy z autorem). Ale jeśli Phase 6+
nie zamknie się z samą Z₂, rozważenie tej ścieżki będzie konieczne.

### Rekomendowany priorytet

**P2 — wysoki.** Bez Phase 6+ TGP jest *teorią grawitacji emergentnej
+ container na SM*, nie *unifikacją*. Decyzja: czy TGP chce być
kompletną teorią materii (wymaga zamknięcia 3c), czy ograniczyć się
do grawitacji (3a/3b wystarczają)?

### Powiązanie z istniejącym audytem

Powiązane z [[L05_mass_exponent_drift]] (k=4 vs p=5−α drift między
sektorami). L05 jest P2 OPEN. Proponuję skonsolidować pod nową klasą
**L07** lub **L08** — "warstwa 3c (kinki jako fermiony) — domknięcie
analityczne".

---

## Punkt 5 — Test różnicujący: |Δg_tt| = (5/6) U³ vs LIGO 3G

### Diagnoza

`TGP_FOUNDATIONS.md` § 3 sygnalizuje konkretną predykcję różnicującą:

> M9.1'' przewiduje explicit |Δg_tt| = (5/6) U³ deviation od GR
> (testowalne LIGO 3G).

Ten jeden konkret odróżnia TGP od GR w obserwowalnym reżimie — i jako
taki jest **najmocniejszym istniejącym kandydatem** na test
falsyfikujący.

### Pliki dotknięte

- `TGP_FOUNDATIONS.md` § 3 (eksplicit deviation)
- `sek08c_metryka_z_substratu` (M9.1'' P1, derywacja deviation)
- `research/op-newton-momentum/M9.x` (PPN higher-order coefficients)
- `PREDICTIONS_REGISTRY.md` (status predykcji)

### Problemy jakie rodzi

1. **Brak explicit sensitivity analysis.** TGP_v1 deklaruje
   "testowalne LIGO 3G", ale nie ma w `research/` policzenia:
   - jaki SNR w Einstein Telescope / Cosmic Explorer wymagany,
   - przy jakich masach BBH/BNS deviation jest detektowalna,
   - jaki strain h(f) odpowiada (5/6) U³ deviation,
   - parameter degeneracy z innymi modyfikacjami GR (ppE framework).

2. **Brak wpisania do PREDICTIONS_REGISTRY z statusem
   FALSIFIABLE.** Predykcja (5/6) U³ powinna być (W)/(P) z explicit
   warunkiem falsyfikacji ("jeśli LIGO 3G nie zobaczy deviation
   o określonej amplitudzie, M9.1'' upada"). Bez tego jest *opisem
   konsekwencji*, nie *kontraktem z obserwacjami*.

3. **Konkurencja modeli.** ppE framework (parametrized post-Einstein,
   Yunes-Pretorius) ma 5 parametrów modyfikujących GR. TGP musi
   pokazać, że (5/6) U³ jest **specyficznym wyborem** w tej
   parametryzacji (które α_ppE, β_ppE odpowiada), żeby być
   rozróżnialna od innych modyfikacji.

4. **Reżim ważności.** U³ to człon 3PN (post-Newtonian rzędu 3).
   LIGO obecnie ma 3.5PN waveformy. LIGO 3G będzie miał ~4PN+.
   Wymaga spójnego rachunku 3PN w TGP (czy TGP daje *wszystkie*
   3PN współczynniki, czy tylko jeden? Pełne 3PN matching to dużo
   pracy).

5. **Dualizm modeli waveform.** Czy w TGP używamy GR waveform
   z deviation (małe perturbacja) czy *pełnego* TGP waveform
   (numerical relativity w TGP)? Brak wyboru jest blokujące.

### Potencjalne ścieżki domknięcia

**Ścieżka A — explicit calc Δh(f) dla GW150914-like event.**
Cykl `op-LIGO-3G-deviation/`. Użyć GW150914 jako szablon (M_1+M_2,
spinów, distance), policzyć GR waveform + (5/6) U³ deviation
w fazie inspiral. Określić, jaki SNR potrzebny do detekcji.
Wynik: konkretna liczba dla Einstein Telescope (ET) i Cosmic Explorer
(CE).

**Ścieżka B — mapowanie na ppE.** Cykl `op-ppE-mapping/`. Pokazać,
że (5/6) U³ deviation odpowiada konkretnym (α_ppE, β_ppE, b) z
Yunes-Pretorius parameterization. To umieszcza TGP w obecnym
benchmarkingu modyfikacji GR.

**Ścieżka C — falsifier statement w PREDICTIONS_REGISTRY.**
Najlżejsza akcja: dodać w `PREDICTIONS_REGISTRY.md` explicit wpis:
"M9.1'' P1 — (5/6) U³ deviation. Falsyfikacja: jeśli ET/CE z 5σ
nie widzi deviation > X·10^(−Y) w fazie inspiral dla M_BBH > Z M_⊙,
M9.1'' jest sfalsyfikowana. Odpowiedzialna grupa: research/op-LIGO-3G."

**Ścieżka D — peer-review submit.** Najmocniejszy ruch: napisać paper
"Strong-field test of M9.1'' metric: predictions for Einstein Telescope"
i wysłać do PRD lub Class. Quantum Grav. To wprowadza TGP do
mainstream falsification machinery.

### Rekomendowany priorytet

**P3 — średni, ale strategicznie ważny.** Pozostałe 4 punkty są
strukturalne (długi teoretyczne). Ten jest *strategiczny* — dotyczy
*sposobu pokazania światu*, że TGP jest falsyfikowalna.

### Powiązanie z istniejącym audytem

Nie ma odpowiednika w S/L/D/M (te są o spójności wewnętrznej, nie
o testach zewnętrznych). Proponuję klasę **T01** — "test różnicujący"
jako pierwsza pozycja w nowej klasie **T** (testy falsyfikujące).

---

## Priorytetyzacja końcowa (v2)

| ID | Skrót | Klasa | Priorytet | Charakter |
|----|-------|-------|-----------|-----------|
| **EXT-1** | Kosmologia radiacyjna z varying c, ℏ, G → BBN/CMB | rozszerzenie L01 + nowy cykl op-FRW | **P1** | **Otwarte ryzyko** (v2): T^μ_μ_EM=0 strukturalnie, ale H_TGP(z) z varying constants może (50–65%) odzyskać BBN/CMB. Wymaga jawnego rachunku — jeszcze nie wykonanego |
| **EXT-2** | Warunek zerowej sumy → derywacja | nowy L07 | **P2** | Aksjomat load-bearing dla całego sek05 — wymagana derywacja z H_Γ albo Lagrange'a multiplier |
| **EXT-3** | M9.1'' z pierwszych zasad | nowy S07 | **P2** | "Emergencja metryki" wymaga explicit derywacji, nie ansatzu po falsyfikacji |
| **EXT-4** | Phase 6+ why_n3 fermiony | rozszerzenie L05 → L08 | **P2** | TGP jako unifikacja vs TGP jako emergentna grawitacja — decyzja statusowa |
| **EXT-5** | LIGO 3G falsifier statement | nowa klasa T01 | **P3** | Strategiczny — peer-review path do mainstream |

## Diagnoza syntetyczna (v2)

TGP_v1 ma **dwa strukturalne długi krytyczne** (EXT-2, EXT-3),
**jeden krytyczny otwartego ryzyka** (EXT-1, status zmieniony w v2),
**jeden strategiczny** (EXT-4) i **jeden taktyczny** (EXT-5). Każdy
z nich jest tematem na cykl badawczy ~3–6 miesięcy. Domknięcie
wszystkich pięciu **zmieniłoby TGP z konstrukcji kompatybilnej
z obserwacjami w teorię wyprowadzoną** — to znaczy: każda predykcja
miałaby derywację z H_Γ + Z₂ + jednej hipotezy metrycznej, bez
trzeciej iteracji ansatzu po falsyfikacji.

**EXT-1 ma teraz status dwustronny** (v2): może być fenomenalnym
sukcesem TGP (BBN+CMB z 1–2 wolnych parametrów dzięki varying c, ℏ,
G), albo silnym sygnałem do pivotu L_mat. Bez explicit FRW calc
nie wiadomo, w którą stronę. Cykl `op-FRW-radiation-era-varying-c/`
jest najpilniejszym zadaniem TGP_v1 z perspektywy zewnętrznej.

Jeśli żaden z trzech P1/P2 strukturalnych długów nie zostanie
domknięty w ciągu 12 miesięcy, TGP pozostanie **dobrze prowadzonym
projektem badawczym** (z metodologią sympy LOCKs i cykli audytu),
ale **nie kandydatem na fundament**. To nie jest porażka — to jest
*honest status*. Najgorsza ścieżka to udawanie, że strukturalne długi
zostały zamknięte przez "annotations only".

## Cross-references

- [[README.md]] — indeks audytu (S/L/D/M klasy)
- [[PRIORITY_MATRIX.md]] — macierz P1/P2/P3/P4 (EXT-1..5 do dodania)
- [[SUMMARY_2026-05-04.md]] — pełny raport ekspercki S/L/D/M
- [[CLOSURE_SUMMARY_2026-05-06.md]] — closure summary L03+D01+M03
- [[L01_rho_operational/README.md]] — ρ operational (rozszerzane przez EXT-1)
- [[L05_mass_exponent_drift/README.md]] — mass exponent drift
  (powiązane z EXT-4)
- [[S01_metric_four_forms/README.md]] — 4 formy metryki
  (CLOSED-RESOLVED via G.0; EXT-3 idzie głębiej)
- [[../TGP_FOUNDATIONS.md]] — aksjomatyka §1, hierarchia §2,
  grawitacja §5, pęd §6
- [[../meta/AUDYT_TGP_2026-05-01.md]] — meta-audit 2026-05-01
- [[../meta/PLAN_DOMKNIECIA_MASTER.md]] — plan domknięcia luk historycznych

## Status recenzji

**Recenzja zewnętrzna z 2026-05-06.** Materiał dla autora TGP do
rozważenia. Nie wymaga akcji w core/ ani w research/ (review-only).
Decyzje dotyczące promocji długów EXT-1..5 do oficjalnych klas
S/L/D/M/T pozostają w gestii autora.

---

**Recenzent:** Claude Opus 4.7 (1M context), sesja 2026-05-06,
na zlecenie autora TGP (Mateusz Serafin).
**Materiał bazowy:** ~50k linii LaTeX/Markdown w `core/`, `audyt/`,
`research/`, `TGP_FOUNDATIONS.md`.
**Czas pracy:** 1 sesja, ~2.5h (v1 + v2 update).
**Status:** review-only, brak edycji w core/.

**Historia rewizji:**
- **v1 (2026-05-06)**: Pierwsza wersja. EXT-1 sklasyfikowany jako
  "WYSOKIE ryzyko falsyfikacji" — argumentacja zakładała implicite,
  że bez ρ_rad jako źródła Φ, H_TGP(z) w erze radiacyjnej nie ma
  fenomenologicznego odpowiednika H_GR(z) ~ a⁻² i BBN/CMB upadają.
- **v2 (2026-05-06)**: Korekta po wskazówce autora. Aksjomaty
  `ax:c–ax:G` w `sek04_stale.tex:29, 250–254` definiują c, ℏ, G jako
  **funkcje pola Φ**, nie stałe fundamentalne. To otwiera mechanizm
  re-strukturyzacji H_TGP(z) z varying constants, niezależny od
  standardowego ρ_rad ~ a⁻⁴. T^μ_μ_EM = 0 trzyma się strukturalnie
  (Weyl-niezmienniczość 4D), ale konsekwencje fenomenologiczne w erze
  radiacyjnej są **otwarte do rachunku**, nie *przegrane a priori*.
  Status EXT-1 zmieniony na "P1 otwarte ryzyko" z subiektywną oceną
  55–65% szansy na fenomenologiczne domknięcie ścieżką A. Pozostałe
  punkty EXT-2..5 bez zmian.
