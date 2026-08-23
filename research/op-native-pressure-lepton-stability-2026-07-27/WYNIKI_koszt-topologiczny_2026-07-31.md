# Domknięcie luki: Q = ładunek topologiczny; topologia kosztuje budżet

> # ⚠ CZĘŚCIOWO WYCOFANE 2026-07-31 — patrz `WYNIKI_psi-orientacja_KOREKTA_2026-07-31.md`
> **Bound N ≥ 5B (koszt symplicjalny) POZOSTAJE WAŻNY.** Wycofane jest natomiast
> **wnioskowanie R ∝ B^{1/3}**: bound objętościowy okazał się **niewiążący** — w ansatzu
> hedgehog nawinięcie B wymaga trawersu kąta **Bπ**, co daje **R ∝ B**, a dla B>1 sferyczny
> hedgehog i tak nie jest konfiguracją minimalną. **Skalowanie R(B) dla B>1 jest OTWARTE.**
> Zdanie „B=1 jest najtańszą topologią (5 sympleksów)" pozostaje w mocy.

**Data:** 2026-07-31
**Typ:** WYNIKI (domknięcie głównej luki z `WYNIKI_budzet-skala`)
**Weryfikacja:** `BUDGET_topological_cost.py` (uruchomione)
**Falsyfikator zadeklarowany przed liczeniem:** brak skończonej minimalnej triangulacji /
brak ograniczenia stopnia liczbą sympleksów ⟹ hipoteza pada. **Nie zadziałał.**

---

## 0. Luka, którą to zamyka

`WYNIKI_budzet-skala` dało **R_min ∝ Q^{1/3}**, ale z **założonym** Q = ∫(ψ−1)dV.
To była główna słabość: „zachowana ilość" wisiała w powietrzu.

**Teraz: Q = B (ładunek topologiczny/nawinięcie), a związek z budżetem jest kombinatoryczny.**

## 1. Mechanizm: reprezentacja topologii kosztuje węzły

Substrat jest **dyskretny** — sek08c: „N_B jest liczbą **topologiczną** (ilość węzłów
w kostce, nie pole fizyczne)". Więc nawinięcie musi być **reprezentowane na skończonej
liczbie węzłów**. To generuje twardy koszt:

> **Bound (rygorystyczny, nie heurystyka):** odwzorowanie symplicjalne stopnia B między
> triangulowanymi S³ — każdy sympleks dziedziny pokrywa **co najwyżej jeden** sympleks celu
> (ze znakiem). Aby pokryć cel **B razy**, potrzeba ≥ B·(liczba sympleksów celu).

**Minimalna triangulacja S³ = brzeg 4-sympleksu** — zweryfikowane:
V=5, E=10, F=10, T=**5**, χ = 5−10+10−5 = **0** ✓ (χ(S³)=0).

⟹ **N(B) ≥ 5B.** Koszt rośnie **liniowo** z ładunkiem topologicznym.

| B | N_min(B) |
|---|---|
| **1** | **5** ← fundamentalny fermion (filar spin-½, B=1) |
| 2 | 10 |
| 3 | 15 |
| 8 | 40 |
| 27 | 135 |

## 2. Stąd rozmiar

Gęstość węzłów n = const (budżet) ⟹ V ≥ N(B)/n = 5B/n ⟹

$$R_{min}(B)=\Big(\frac{3\cdot 5B}{4\pi n}\Big)^{1/3}\ \propto\ B^{1/3}$$

Zweryfikowane numerycznie: R_min(B)/R_min(1) = **B^{1/3} dokładnie** (1, 1.2599, 1.4422,
2.0000, 3.0000, 4.0000 dla B=1,2,3,8,27,64).

## 3. Co to znaczy (i dlaczego jest mocniejsze niż poprzednio)

- **„Zachowana ilość" nie jest już dowolna.** To **nawinięcie**, które (a) już mamy
  (π₃(M)=ℤ), (b) już jest **chronione topologicznie** (nie da się go usunąć), (c) jest
  tym samym B, z którego pochodzi **spin-½** (filar).
- **Topologia kosztuje budżet** — to jest szukany most między warstwą topologiczną
  a warstwą budżetową. Obie były osobno ustalone; teraz są sprzężone.
- **B=1 jest najtańszą możliwą topologią** (5 sympleksów) ⟹ fundamentalny fermion jest
  **obiektem minimalnego budżetu**. To natywnie tłumaczy, dlaczego B=1 jest poziomem
  podstawowym, bez postulowania.

## 4. Stan wieży po tym kroku

| Piętro | Mechanizm | Status |
|---|---|---|
| substrat | węzły + budżet informacyjny (N_B topologiczne, B=N_B·s₀) | rdzeń |
| wstążki | dysklinacje reperu, π₁=2T=Q₈⋊ℤ₃ → **spin + kolor** | ✅ (2T wymuszone) |
| domykanie | reguła kombinatoryczna: pojedyncza wstążka nie domyka się; trialność-0 konieczna | ✅ |
| lokalizacja | **kolaps** ⟂ cap budżetowy; **rozpad** ⟂ przeszkoda topologiczna | ✅ obustronnie |
| skala | **R ∝ B^{1/3}**, koszt N ≥ 5B | ✅ ten dokument |

**Żadne piętro nie używa: minimalizacji energii, członu Skyrme'a, płaskiego tła, fitu.**

## 5. Ograniczenia (uczciwie)

1. **Bound N ≥ 5B dotyczy odwzorowań symplicjalnych** — substrat niekoniecznie jest
   symplicjalny. To **idealizacja**. **Skalowanie (∝B) jest odporne; stała 5 zależy
   od dyskretyzacji** i nie należy jej traktować jako predykcji liczbowej.
2. Daje **dolne** ograniczenie rozmiaru (anty-kolaps), **nie pełny profil** obiektu.
3. **Gęstość węzłów przyjęta stała**; jej związek z ψ (n ∝ 1/h, z sek08c) **nie został
   tu użyty** — pełne sprzężenie budżet↔profil pozostaje **OTWARTE**.
4. R ∝ B^{1/3} jako *kształt* jest generyczne dla nasycenia; natywne jest **źródło**
   (koszt topologiczny + budżet), nie wykładnik.

## 6. Jednozdaniowo

> Zachowaną „ilością" jest **ładunek topologiczny B**, a most do budżetu jest
> kombinatoryczny: reprezentacja nawinięcia stopnia B na dyskretnym substracie wymaga
> **N ≥ 5B** sympleksów, więc przy stałej gęstości węzłów **R_min ∝ B^{1/3}** — przy czym
> **B=1 (spin-½) okazuje się najtańszą możliwą topologią**, czyli naturalnym poziomem
> podstawowym.

---

**Plik:** `BUDGET_topological_cost.py`
**Powiązane:** `WYNIKI_budzet-skala_2026-07-31.md` (cap, domknięcie obustronne),
`ONTOLOGIA_energia-relacyjna-budzet_2026-07-31.md`, `research/qm_spin/README.md` (B=1, filar)
