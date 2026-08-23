# Dynamika wstążek — Krok 2: skończony rozmiar (skyrmion B=1)

**Data:** 2026-07-30
**Typ:** WYNIKI (domykanie strzałki substrat→wstążki, warstwa dynamiczna, krok 2/2)
**Model:** Skyrme na M (M̃=S³), hedgehog B=1 (decyzja użytkownika)
**Weryfikacja:** `RIBBONS_stepB2_skyrme_v3.py` (uruchomione, z kontrolami artefaktu)

---

## 0. Uczciwa historia rachunku (dyscyplina)

To NIE poszło od razu — i to jest ważne dla wiarygodności:

- **v1** (L-BFGS, init = liniowa rampa π(1−r/R), pudło R=30–60) **UTKNĄŁ** i dał
  **podręcznikowy artefakt pudła**: `<r>≈R/2` (6.3/9.0/14.3/19.5/29.6 dla R=15..60),
  E∝R, wirial rosnący 276→1239. Skrypt miał ZASZYTY fałszywy print „minimum przy λ=1".
  **Złapane i odrzucone** — to dokładnie ten błąd (narracja > dowód), który zabił bounce-hierarchy.
- **v2** (solve_bvp) **nie zbiegł** (osobliwość r=0 + błąd initu).
- **Przyczyna obu:** zły init + zbyt duże pudło (rozmiar skyrmionu ~1, nie 30).
- **v3** (relaksacja, init zlokalizowany 4·arctan(e⁻ʳ), pudło R~8–12): **znajduje prawdziwy
  skyrmion B=1.**

## 1. Wynik (v3, z kontrolami)

**Skan pudła R** (kontrola artefaktu — tu v1 poległo):

| R | E | E2 | E4 | wirial E2/E4 | ⟨r⟩ | B |
|---|---|---|---|---|---|---|
| 6 | 11.715 | 5.982 | 5.734 | 1.043 | 0.954 | 1.000 |
| 8 | 11.691 | 5.958 | 5.734 | 1.039 | 0.954 | 1.000 |
| 10 | 11.690 | 5.957 | 5.734 | 1.039 | 0.954 | 1.000 |
| 12 | 11.690 | 5.957 | 5.734 | 1.039 | 0.954 | 1.000 |
| 16 | 11.690 | 5.956 | 5.734 | 1.039 | 0.954 | 1.000 |

- **E, ⟨r⟩, wirial STAŁE względem R** → zlokalizowany soliton O(1), **NIE artefakt pudła**.
- **B = 1.000** (stopień π₃) → to jest ładunek Q=1 z kroku 1.
- **⟨r⟩ = 0.954** (rozmiar O(1)), nie ~R/2.

**Zbieżność siatki N** (R=10): E→11.690, ⟨r⟩→0.954 zbiegają (N=1600+); B=1 dokładnie.

**Derrick / balans:** E2=5.96 i E4=5.73 **porównywalne** (E4/E2≈0.96) → człon Skyrme'a
**realnie balansuje** σ (w v1 był znikomy 0.42 vs 117 = reżim zapadania). Krzywa
E(λ)=λ·E2+E4/λ ma minimum przy λ≈1; sam σ maleje monotonicznie z λ→0 (zapada). To Derrick w równowadze.

## 2. Klasyfikacja audytowa (co DERIVED, co WŁOŻONE)

**DERIVED:**
- Q=1 wstążka realizuje się jako **metastabilny, zlokalizowany soliton skończonego rozmiaru**
  (skyrmion B=1), niezależny od pudła i zbieżny z siatką.
- To **ten sam hedgehog B=1**, co filar spinu (`qm_spin`) → spin-½ (FR) i wstążka to **jeden obiekt**.

**WŁOŻONE / OGRANICZENIA (nie zamiatam):**
- **Człon Skyrme'a (4-pochodne) jest INPUTEM funkcjonału energii, nie wyprowadzony z reguł
  substratu.** Krok 1 (nieusuwalność via π₃) jest czysto topologiczny; ale skończony ROZMIAR
  wymaga tego założonego członu. To jest dokładnie otwarty §5.4 (minimalny background-free
  zestaw ruchów Φ) — czy człon Skyrme'a jest „dozwolony/minimalny" pozostaje OTWARTE.
- **Wirial = 1.039, nie 1.000** — ~4% resztka z ucięcia pudła (Dirichlet f(R)=0 tnie potęgowy
  ogon) i regularyzacji r=0. Zbieżny, ale nie zeruje się. **Nie twierdzę „wirial=1 dokładnie".**
- To **odtwarza standardowy, podręcznikowy skyrmion B=1** — nowość NIE jest w istnieniu solitonu
  (znane), lecz w tym, że wstążka mapuje się na niego, a target M=S³/2T niesie kolor (2T).

## 3. Stan strzałki substrat→wstążki po kroku 2

| Warstwa | Status |
|---|---|
| Szkielet topologiczny (target 2T wymuszony, spin+kolor z 1 grupy) | ✅ DERIVED |
| Metastabilność-przed-usunięciem (nieusuwalność) | ✅ DERIVED (π₃ samo-splot, spięte ze spinem) |
| Skończony rozmiar (dynamika) | ✅ pokazane — ale przez WŁOŻONY człon Skyrme'a |
| Background-free / minimalny zestaw ruchów Φ (§5.4) | ⚠ OTWARTE (i teraz spięte z „skąd człon Skyrme'a") |

## 4. Jednozdaniowo

> Wstążka Q=1 realizuje się jako metastabilny, zlokalizowany skyrmion B=1 (E,⟨r⟩,B niezależne
> od pudła, zbieżne) — ten sam hedgehog, co filar spinu; nieusuwalność jest czysto topologiczna
> (π₃, krok 1), ale skończony rozmiar wymaga **włożonego** członu Skyrme'a, co przenosi ciężar
> na otwarty §5.4 (minimalny, background-free zestaw ruchów Φ).

---

**Pliki:** `RIBBONS_stepB2_skyrme_v3.py` (kanoniczny; v1/v2 usunięte jako nieudane — lekcja w §0)
**Powiązane:** `WYNIKI_stepB_dynamika-topologia_2026-07-30.md` (krok 1),
`WYNIKI_substrat-wstazki_2T_2026-07-30.md` (szkielet), `research/qm_spin/README.md` (filar)
