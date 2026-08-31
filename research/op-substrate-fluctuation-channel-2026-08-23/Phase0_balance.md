# Phase 0 — LOCK: op-substrate-fluctuation-channel-2026-08-23

**Status: PHASE0-LOCKED 2026-08-23. Zero obliczeń przed zamknięciem tego pliku.**
**Autoryzacja:** user 2026-08-23: „Zapisz to i zajmij się w tej sesji 1+2 jak
proponowałeś" (burza mózgów: [[../../meta/BRAINSTORM_2026-08-23_brakujace-puzzle.md]],
punkty 1+2).
**Rejestr:** structural-emergence (inwentarz kanałów na poziomie 0), NIE
empirical-novelty. Zero claimów obserwacyjnych.

---

## 1. Kontekst i delimitacja

- **Poziom:** 0 (substrat, statystyka H_Γ / warstwa Gaussowska dodatekB),
  NIE poziom 1 (akcja continuum, solitony). Cykl jest KOMPLEMENTARNY do
  `op-bath-two-sectors-2026-08-23` (poziom 1, tła ψ_n, osobny agent) —
  **zakaz dublowania jego Q1/Q2**.
- **Luka adresowana:** bilans programu „most do grawitacji" (STATE #64+):
  przetestowane kanały oddziaływań = amplitudowy (Yukawa m~m₀, sektor
  solitonowy), Goldstone (ładunkowy — odpada), geometryczny (refrakcja —
  jedyny pozytywny). **Kanał fluktuacyjny (indukowany fluktuacjami
  substratu, typ Casimira) nigdy nie był liczony.** Jednocześnie dodatekB
  wymaga bliskości krytyczności dla granicy continuum
  (`prop:continuum-conditions`) i daje m_eff² = γ(1+T_Γ)
  (`cor:entropy-potential`) z γ ≈ 12Λ_eff/Φ₀ — czyli ξ=1/√γ jest skali
  kosmologicznej, a substrat wewnątrz horyzontu jest efektywnie bliski
  krytyczności. Te dwa fakty korpusu nigdy nie zostały połączone
  z inwentarzem kanałów.
- **Druga luka (idea 1, poziom 0):** czy statystyka substratu w kąpieli
  o gęstości Φ_b ≠ Φ_vac zmienia KRZYWIZNĘ potencjału efektywnego
  (możliwa emergencja znaku tachionowego z rozrzedzenia tła) — na
  poziomie MFT z wiązaniem gradientowym v2. To jest poziom-0 odpowiednik
  hipotezy „dwóch sektorów" (N2), ale w INNEJ maszynerii (statystyka,
  nie ODE solitonowe) — wynik dowolnego znaku nie przesądza Q2 tamtego
  cyklu i odwrotnie.

## 2. Pytania binarne

> **QF:** czy na poziomie Gaussowskim substratu istnieje kanał
> oddziaływania defekt–defekt o **uniwersalnym znaku** (przyciąganie
> niezależne od „ładunku"/znaku defektu), odróżnialny od kanału
> wymiany klasycznej sygnaturą zasięgu (2m vs m) i wykładnikiem
> potęgowym na krytyczności (−2 vs −1)?

> **QB:** czy w MFT substratu v2 (on-site + wiązanie gradientowe
> `rem:B-v2-status`) krzywizna potencjału efektywnego w punkcie
> jednorodnym może stać się **ujemna** (tachionowa) dla kąpieli
> rozrzedzonej Φ_b < Φ_vac — i czy wkład wiązania gradientowego ma
> określony znak (stabilizujący/destabilizujący) dla wszystkich Φ_b ≥ 0?

## 3. Modele (ZAMKNIĘTE — zero zmian po starcie obliczeń)

### 3.1 QF — warstwa Gaussowska substratu

Pole fluktuacji δφ wokół próżni jednorodnej na sieci kubicznej L³
(periodycznej), energia dyskretna
`H = Σ_x [ ½ Σ_a (δφ(x+ê_a)−δφ(x))² + ½ m² δφ(x)² ]`, β_T = 1
(normalizacja K(1)=1 — zgodna z propagatorem `eq:B-propagator` dodatekB;
skala T wchodzi liniowo w F i nie zmienia znaków/wykładników).
Funkcja Greena sieci: `G_m(r) = (1/L³) Σ_k e^{ikr} / (m² + Σ_a 2(1−cos k_a))`.

**Przypadek krytyczny:** m = 0 z jawnie usuniętym modem zerowym k=0
(projekcja; fizycznie: ustalona średnia globalna). Przypadki masywne:
pełna suma.

Dwa typy defektu punktowego (deklaracja modelowa):
- **M-S (źródłowy):** sprzężenie liniowe −q_i δφ(x_i). Oddziaływanie
  klasyczne (wymiana): `U_S(d) = −q₁q₂ G_m(d)` (do wyprowadzenia exact
  w Phase 1).
- **M-D (kotwiczący/pinning):** δφ(x_i) = v_i (zamrożenie wartości pola —
  model rdzenia solitonu zamrażającego substrat). Energia swobodna
  z rozkładu brzegowego Gaussa: `F(d) = ½ ln det Σ(d) + ½ vᵀ Σ(d)⁻¹ v`,
  `Σ(d) = [[G(0), G(d)],[G(d), G(0)]]`; oddziaływanie
  `F_int(d) = F(d) − F(∞)`. Część **fluktuacyjna** (v-niezależna):
  `F_fl(d) = ½ ln(1 − (G(d)/G(0))²)`.

### 3.2 QB — MFT v2 w kąpieli

Węzeł s w kąpieli z zamrożonych sąsiadów s_b (koordynacja z=6, jak
dodatekB), hamiltonian v2:
`V_eff(s; s_b) = (m₀²/2)s² + (λ₀/4)s⁴ + zJ s² s_b² (s²−s_b²)²`.
Zmienne zredukowane r=m₀²/J, u=λ₀/J; punkt odniesienia = punkt stały WF
dodatekB: r*=−2.251, u*=3.917 (`eq:B-WF`), plus skan (r,u) w otoczeniu.
Wielkość badana: `C(s_b) = ∂²V_eff/∂s²|_{s=s_b}` (krzywizna w punkcie
jednorodnym = lokalna stabilność konfiguracji jednorodnej o gęstości
Φ_b=s_b²), rozłożona na wkład on-site i wkład wiązania.

## 4. Kryteria (PRE-REJESTROWANE)

### QF (werdykt = koniunkcja QF-1..QF-4)

- **QF-1 (uniwersalność znaku):** `F_fl(d) < 0` i monotonicznie rosnące
  do 0 z d (= przyciąganie) dla WSZYSTKICH m ∈ {0, 0.05, 0.2, 0.5}
  i wszystkich d w oknach §5; jednocześnie część klasyczna M-D
  (v-zależna) i kanał M-S zmieniają znak z sgn(v₁v₂) / sgn(q₁q₂)
  (tabela znaków: (v₁,v₂) ∈ {(a,a),(a,−a),(a,0)}, a=0.3). Osiągalny
  FAIL: niemonotoniczność, znak zależny od v w części fluktuacyjnej,
  albo brak flipu znaku w kanale klasycznym.
- **QF-2 (sygnatura zasięgu):** dopasowanie
  `ln|F_fl| = a − κ_D d − 2 ln d` daje κ_D = 2μ(m) ± 10%,
  a `ln G_m = a − κ_S d − ln d` daje κ_S = μ(m) ± 10%, gdzie
  μ(m) = arccosh(1+m²/2) (dokładny zasięg sieciowy); dla m ∈ {0.2, 0.5}.
  Osiągalny FAIL: κ_D ≈ κ_S (brak struktury „podwójnej wymiany").
- **QF-3 (wykładniki krytyczne):** przy m=0: slope ln|U_S| vs ln d
  = −1.0 ± 0.1; slope ln|F_fl| vs ln d = −2.0 ± 0.2; dryf wykładnika
  między L=64 a L=128 < 0.1. Osiągalny FAIL: wykładniki równe, poza
  tolerancją, lub L-zależne.
- **QF-4 (kontrole negatywne):** (a) m=2.0 (ξ<1): |F_fl(d=6)|/|F_fl(d=2)|
  < 10⁻³ (brak dalekozasięgowego artefaktu FFT); (b) kontrola obrazów
  periodycznych: |G^{L=64}(d)−G^{L=128}(d)|/|G^{L=128}(d)| < 10⁻² dla
  d ≤ 16, m ≥ 0.2; (c) sanity: G_m(d) > 0 dla wszystkich d (funkcja
  Greena masywnego operatora dodatnia).

**QF-PASS** = QF-1 ∧ QF-2 ∧ QF-3 ∧ QF-4. **QF-FAIL** = ¬QF-1 lub ¬QF-2
lub ¬QF-3 przy przejściu QF-4. **QF-INCONCLUSIVE** = pad QF-4 (problem
maszynerii, nie fizyki) — wtedy żadnych wniosków.

### QB (werdykt = QB-1 i QB-2 rozstrzygane NIEZALEŻNIE)

- **QB-1 (znak wkładu wiązania, sympy exact):** wyznaczyć
  `ΔC_bond(s_b) = ∂²[zJ s²s_b²(s²−s_b²)²]/∂s²|_{s=s_b}` symbolicznie.
  Werdykt binarny: ΔC_bond ≥ 0 dla wszystkich s_b (wiązanie gradientowe
  STABILIZUJE punkt jednorodny — emergencja tachionu z wiązania
  niemożliwa na MFT) vs istnieje s_b z ΔC_bond < 0. Oba wyniki
  informatywne; zgłosić wprost.
- **QB-2 (próg rozrzedzenia):** czy istnieje Φ_c ∈ (0, Φ_vac) takie, że
  C(s_b) < 0 dla Φ_b < Φ_c przy (r*, u*, z=6)? Jeśli tak: wyznaczyć
  Φ_c/Φ_vac (dokładne równanie + wartość numeryczna) oraz mapę znaku
  C na skanie r ∈ [−3, −1], u ∈ [2, 6] (siatka 41×41, z=6).
  Osiągalny FAIL: C ≥ 0 wszędzie (brak progu) albo Φ_c ≥ Φ_vac
  (jednorodna próżnia sama niestabilna — patologia do zgłoszenia).
- **QB-3 (odczyt, deskryptywnie, ZERO claimów):** czy wzór
  {kąpiel pełna → C>0 stabilny; rozrzedzenie → C<0 tachionowy} jest
  strukturalnie równoległy do hipotezy dwóch sektorów N2 — wyłącznie
  komentarz porównawczy, bez przenoszenia na poziom 1.

## 5. Okna dopasowań (ZAMKNIĘTE przed obliczeniami)

| Przypadek | L | okno d (oś sieci) |
|---|---|---|
| m=0.5, fit κ | 64 | d ∈ [2, 8] |
| m=0.2, fit κ | 64 | d ∈ [5, 15] |
| m=0.05 (tylko QF-1) | 128 | d ∈ [10, 30] |
| m=0, fit potęgowy | 64 | d ∈ [6, 16] |
| m=0, fit potęgowy | 128 | d ∈ [8, 24] |

Zakaz przesuwania okien po obejrzeniu danych. Jeśli fit w oknie jest
jawnie zdominowany artefaktem (R² < 0.99), zgłosić jako QF-INCONCLUSIVE
dla tego punktu — NIE dobierać okna.

### Amendment A1 (PRZED jakimkolwiek obliczeniem, 2026-08-23) — estymator przypadku krytycznego

Projekcja modu zerowego (m=0, §3.1) wnosi do propagatora skończonego L
znaną, r-NIEZALEŻNĄ stałą addytywną (efekt skończonego rozmiaru rzędu
1/L, porównywalny z 1/(4πd) na końcu okien §5) — czysty fit potęgowy
byłby zdominowany tą stałą z powodów maszynerii, nie fizyki. Dlatego
QF-3 używa estymatorów:
- **S-typ:** fit nieliniowy `G(d) = A·d^p + B` (3 parametry) w oknie §5;
  wymóg: p = −1.0 ± 0.1; dryf p między L=64 a L=128 < 0.1.
- **D-typ:** propagator połączony `G_c(d) = G(d) − B̂`, `G_c(0) = G(0) − B̂`
  (B̂ z fitu S-typu); `F_fl` z g_c = G_c(d)/G_c(0); wymóg: slope
  ln|F_fl| vs ln d = −2.0 ± 0.2 w oknie §5.
Falsyfikowalność zachowana: gdyby propagator sieciowy nie był ~1/d
(np. ~1/d² albo wykładniczy), p wypada poza tolerancję. Przypadki
masywne: bez zmian (pełne sumy, bez stałej). Amendment zamknięty przed
napisaniem jakiegokolwiek skryptu; żadne dane nie były obejrzane.

## 6. Drzewo decyzyjne

- **QF-PASS:** poziom 0 ma czwarty kanał o sygnaturze grawitacyjnej
  (uniwersalne przyciąganie); do NEEDS: (i) skalowanie dla obiektów
  ROZCIĄGŁYCH (punktowy daje potencjał ~d⁻² na krytyczności — to NIE
  jest jeszcze Newton d⁻¹; jawnie odnotować), (ii) N-ciałowość/addytywność,
  (iii) związek z γ ≈ 12Λ_eff/Φ₀ (zasięg kosmologiczny) — user-gate.
- **QF-FAIL:** kanał fluktuacyjny nie ma deklarowanej sygnatury —
  inwentarz kanałów poziomu 0 zamknięty negatywnie; zgłosić wprost.
- **QB-2 z progiem Φ_c:** poziom 0 daje ilościowy mechanizm „tachion
  z rozrzedzenia" (spinodala substratu) — dopisek porównawczy do analizy
  N2 możliwy TYLKO jako user-gated NEED (rdzeń nietykany).
- **QB-1 ≥ 0 wszędzie:** znak tachionowy W maszynerii 2 NIE pochodzi
  z wiązania gradientowego substratu na MFT — zawęża pochodzenie
  (pcha na poziom 1 lub poza-MFT); wynik negatywny = informacja.

## 7. Forbidden moves

1. Zmiana kryteriów §4, okien §5, modeli §3 po starcie obliczeń.
2. Dublowanie Q1/Q2 cyklu `op-bath-two-sectors` (poziom 1) — zakaz.
3. Edycje rdzenia .tex — zakaz (wszystko przez NEEDS, user-gate).
4. Claimy grawitacyjne/obserwacyjne z wyników poziomu 0 — zakaz
   (dozwolony wyłącznie inwentarz kanałów + odczyt deskryptywny).
5. Każdy skrypt musi mieć test zdolny dać FAIL (lekcja A4
   op-lattice-bath-runaway); SUMMARY musi cytować faktyczne liczby.
6. Wynik negatywny → zgłoszenie wprost w Phase_FINAL_close.
7. Commit/push — tylko za jawną zgodą użytkownika.

## 8. Fazy i deliverables

- **Phase 1 (QF analitycznie, sympy):** wyprowadzenie U_S(d) i F(d)
  M-D z całek Gaussa (exact, macierz 2×2); rozwinięcie F_fl przy małym
  G(d)/G(0); formy kontinuum (Yukawa e^{−mr}/4πr; −½(G/G(0))²
  ⟹ e^{−2mr}/r²); → `Phase1_analytic.py` + `Phase1_output.txt`.
- **Phase 2 (QF numerycznie):** G_m przez FFT (L=64/128), tabela znaków,
  fity κ i wykładników wg §4–§5; → `Phase2_lattice.py` + `Phase2_output.txt`.
- **Phase 3 (QB):** sympy exact ΔC_bond i C(s_b); próg Φ_c przy (r*,u*);
  mapa znaku na skanie (r,u); → `Phase3_bath_curvature.py` + `Phase3_output.txt`.
- **Zamknięcie:** `Phase_FINAL_close.md` (werdykty QF/QB + drzewo §6),
  `NEEDS.md` (user-gated), wpis STATE.md.

*LOCK zamknięty 2026-08-23 przed napisaniem jakiegokolwiek kodu.*
