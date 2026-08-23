# Dynamika wstążek — Krok 1: test topologiczny stabilizatora (B)

**Data:** 2026-07-30
**Typ:** WYNIKI (domykanie strzałki substrat→wstążki, warstwa dynamiczna, krok 1/2)
**Weryfikacja:** `RIBBONS_stepB_topology.py` (inputy grupowo-topologiczne u źródła)
**Status:** krok 1 z 2. Rozstrzyga MECHANIZM ochrony; nie domyka jeszcze rozmiaru/energii.

---

## 0. Pytanie kroku 1

Czy stabilizator (B) — „framing wstążki = centralny −1 z 2T (ten sam element, co daje spin)"
— daje **realną** przeszkodę topologiczną przed usunięciem (zapadnięciem do próżni),
czy jest życzeniowy? (Gdyby życzeniowy → odpada, schodzimy do (A) Skyrme.)

## 1. Baza topologiczna (rozstrzygająca)

**M = SO(3)/T = S³/2T** — sferyczna forma przestrzenna. 2T ⊂ SU(2) działa **swobodnie**
przez lewe mnożenie (lewe mnożenie w grupie jest zawsze wolne) ⟹ nakrycie uniwersalne
**M̃ = S³**. Stąd, bez zgadywania:

| grupa homotopii | wartość | konsekwencja fizyczna |
|---|---|---|
| π₁(M) | **2T** (rząd 24) | typy dysklinacji (7 klas) — z poprzedniego wyniku |
| π₂(M) | **0** | **brak defektów punktowych, brak ucieczki monopolowej** |
| π₃(M) | **ℤ** | ładunek samo-splotu (Hopf/skyrmion) |

Zweryfikowane u źródła: |2T|=24, −1 ∈ 2T centralny rzędu 2; framing 2π → −1, framing 4π → +1.

## 2. Wynik (uczciwie, z rozróżnieniem)

**π₂ = 0 działa dwustronnie — i to jest sedno:**

- **Sama klasa π₁ NIE wystarcza** do ochrony zwartej pętli. Goła pętla klasy −1 może się
  zagoić przy kurczeniu (brak przeszkody π₂). Więc „bariera z gołej klasy dysklinacji"
  byłaby przeceniona — nie robię tego błędu.
- **Realna ochrona idzie przez π₃ = ℤ = samo-splot (self-linking).** Kluczowa tożsamość:
  framing rzutuje na π₁ przez **parzystość samo-splotu SL**:
  > SL **nieparzysty** ⟺ element **−1** (spin) ⟺ **Q_π₃ ≠ 0** ⟺ pętla **nieusuwalna**.
  To jest **dokładnie relacja Finkelsteina-Rubinsteina (−1)^B** — ten sam mechanizm, co filar
  spinu (`qm_spin`). Ładunek, który czyni wstążkę zachowaną cząstką, = jej samo-splot =
  spięty z jej spinem. **Jedna topologia** (duch §0).

## 3. Werdykt kroku 1

> **(B) jest REALNE — ale poprzez π₃-samo-splot, nie przez gołą klasę π₁.**
> Framing −1 (spin) ⟺ nieparzysty samo-splot ⟺ niezerowy ładunek π₃=ℤ ⟺ wstążka
> **nie może zostać usunięta** (ładunek zachowany). Metastabilność-przed-rozpadem jest
> **topologiczna, nie hack** (§5.1 ✅), i **spięta ze spinem** (§5.6, filar) — jedna struktura.

**CAVEAT (nie zamiatam):** ładunek π₃ jest skyrmionowy → **Derrick kurczy go do punktu**.
Czyli (B) daje „nieusuwalny", ale **nie** „skończony preferowany rozmiar". **Rozmiar wymaga
skali** (człon typu Skyrme). To jest **Krok 2** — i teraz jest on uzasadniony: jest realny,
topologicznie chroniony obiekt, któremu brakuje tylko skali (odwrotność błędu bounce-hierarchy,
gdzie „stabilizowano" artefakt).

## 4. Klasyfikacja audytowa

**DERIVED (topologia, zweryfikowane):**
- M = S³/2T; π₁=2T, π₂=0, π₃=ℤ.
- Ochrona wstążki = π₃ samo-splot; framing −1 ⟺ SL nieparzysty ⟺ Q≠0 ⟺ nieusuwalna.
- Spięcie metastabilności ze spinem (FR (−1)^B) — jedna topologia.

**OTWARTE → Krok 2:**
- Skończony rozmiar (Derrick) — potrzebny człon skali (Skyrme-Faddeev-podobny).
- Rachunek energii **z kontrolą pudła R i zbieżności siatki N od startu** (test artefaktu,
  którego brak zabił bounce-hierarchy).
- Explicite: skonstruować pole S³→M z Q=1 (nieparzysty SL) i pokazać metastabilne minimum
  energii przy skończonym rozmiarze; sprawdzić niezależność od R.

---

**Pliki:** `RIBBONS_stepB_topology.py`
**Powiązane:** `WYNIKI_substrat-wstazki_2T_2026-07-30.md`, `research/qm_spin/README.md` (FR/filar)
