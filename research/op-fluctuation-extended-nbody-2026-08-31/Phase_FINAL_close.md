# Phase FINAL — zamknięcie: op-fluctuation-extended-nbody-2026-08-31

**Status: CLOSED-EXECUTED 2026-08-31 (LOCK → Phase 1–3 → close, jedna sesja).**
**Werdykty: QE (pytanie N1 rodzica): NIE — czysty wynik negatywny.
QN (pytanie N2 rodzica): TAK — nieaddytywność z uniwersalnym znakiem dodatnim.**

Rejestr: structural-emergence, poziom 0. Zero claimów obserwacyjnych.
LOCK: [[Phase0_balance.md]] — kryteria i okna nietknięte po starcie obliczeń.

---

## 1. QE — obiekty rozciągłe: reżim −1/d NIE istnieje (dwuciałowo)

Pytanie binarne N1 rodzica: czy dla kul zamrożonego pola promienia R istnieje
reżim (R, d), w którym kanał fluktuacyjny na krytyczności daje −1/d.

**Odpowiedź: NIE** (maszyneria 9/9 PASS; L=128, klastry R∈{1,2,3} =
7/33/123 węzłów, log-det pełnych macierzy kowariancji; krytycznie
propagator connected z B̂ = −1.564e−3 wg odziedziczonego Amendment A1):

- **Wykładnik dalekiego pola pozostaje −2:** slope = −2.070 (R=1),
  −2.091 (R=2), −2.147 (R=3); wszystkie R² ≥ 0.99986; dryf L=96↔128
  (R=2) = 0.054 < 0.1.
- **Przy kontakcie wykładnik lokalny STROMIEJE (do ~−3), nigdy nie
  łagodnieje do −1:** pełne profile p_loc(d) w [[Phase2_output.txt]];
  zakres od −3.05 (d=2R+3) monotonicznie do −2.01 (d=27); najdłuższy
  przebieg |p_loc+1|≤0.15 = **0** dla każdego R (kryterium QE-4:
  wymagane ≥3).
- **R wchodzi wyłącznie w amplitudę, zgodnie z obrazem pojemnościowym**
  (Phase 1 P1-4 exact: det(I−c²Σ_A⁻¹JΣ_B⁻¹J) = 1−c²C_AC_B):
  A(R)/[½C_R²(4π)⁻²] = 0.933 / 0.938 / 0.962 dla R=1/2/3 (zbiega do 1
  z R — poprawki wyższych multipoli maleją).
- Uniwersalność znaku (F<0, monotoniczne) i sygnatura zasięgu 2μ
  **przeżywają rozciągłość**: κ_D = 2μ ± 1.9%/3.7% (R=1/2, m=0.2).

Konsekwencja dla programu „most do grawitacji": Newtonowskiego −1/d nie
da się uzyskać z kanału fluktuacyjnego przez samą rozciągłość dwóch ciał;
pozostałe nieliczone ścieżki to efekty zbiorowe/ośrodkowe (N-ciałowość
w skończonej gęstości defektów, ekranowanie) lub poziom 1 — patrz NEEDS.

## 2. QN — N-ciałowość: nieaddytywny, znak dodatni (osłabia), kontrast z klasycznym

- **Phase 1 (sympy exact, 5/5):** ΔF₃ = g₁₂g₁₃g₂₃ − ½(g₁₂²g₁₃² + g₁₂²g₂₃²
  + g₁₃²g₂₃²) + O(g⁵); człon wiodący rzędu 3, **dodatni** dla g_ij>0.
- **Phase 3 (numerycznie, 3/3):** znak ΔF₃ **dodatni we wszystkich 20
  punktach** (m∈{0.2, krytyczny-connected} × {trójkąt T, kolinearna C} ×
  d∈{4..12}); wielkość na krytyczności |ΔF₃|/|ΣF_par| = 1.1–4.0%
  (malejąca ~1/d: slope ln|ΔF₃| = −3.03 vs analitycznie −3);
  zgodność z rozwinięciem rzędu ≤4: 1.6e−4 względnie (m=0.2, T, d=8).
- **Kontrast:** kanał klasyczny (źródłowy) addytywny DOKŁADNIE (<1e−12,
  superpozycja Gaussowska potwierdzona numerycznie).

Odczyt deskryptywny (zero claimów): kanał fluktuacyjny ma słabą,
uniwersalnie ODPYCHAJĄCĄ poprawkę trójciałową — dla „grawitacyjnej"
interpretacji kanału oznacza to odchył od superpozycji ~1% przy d~10
i zanikający z odległością; NIE przenosić na sektor solitonowy
(op-nbody-additivity liczył INNY kanał — klasyczny, addytywny z poprawką
~exp; poziomy rozdzielone).

## 3. Incydenty (wzorzec jawności A3/A4)

- **P1-5a (bug maszynerii, przed werdyktem):** pierwotna wersja
  Phase1_analytic.py wywołała `sp.diff(., sp.log(d))` — błąd składniowy
  sympy (crash, nie FAIL testu). Poprawka: pochodna logarytmiczna
  d·d/dd ln|F|; kryterium P1-5 bez zmian. Pierwotny traceback w historii
  sesji; skorygowany przebieg = [[Phase1_output.txt]].

## 4. Anti-Lakatos

✓ LOCK zamknięty przed kodem; zero zmian kryteriów/okien po starcie.
✓ Wynik negatywny QE zgłoszony wprost jako główny wynik cyklu.
✓ Kontrole negatywne wykonane: QE-0 (tożsamość R=0, <1e−16), Fischer
(F<0 wszędzie), dryf L=96↔128, QN-3 (addytywność klasyczna <1e−12).
✓ Rdzeń .tex NIETKNIĘTY (wnioski wyłącznie przez NEEDS, user-gate).
✓ Zakaz dublowania op-bloch-chain-stability (równoległa sesja) dotrzymany
— katalog nietykany, STATE.md odłożony do zwolnienia.
✓ Okna i profile raportowane w całości (bez selekcji punktów).

## 5. Deliverables

[[Phase0_balance.md]] · [[Phase1_analytic.py]] + [[Phase1_output.txt]] ·
[[Phase2_extended.py]] + [[Phase2_output.txt]] · [[Phase3_nbody.py]] +
[[Phase3_output.txt]] · [[NEEDS.md]] · [[README.md]]
