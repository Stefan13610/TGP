# Profil z budżetu: kształt, rozmiar i głębokość obiektu bez energii

> # ⛔ WYCOFANE 2026-07-31 — czytać `WYNIKI_psi-orientacja_KOREKTA_2026-07-31.md`
> Domknięcie sprzężenia ψ↔orientacja wykazało, że **metryka nie widzi orientacji**
> (g_ij = h·δ_ij niezależnie od R), więc **topologia siedzi w orientacji, NIE w ψ**.
> Wymóg „trawersu ψ", na którym oparty jest CAŁY ten dokument, **nie jest wymogiem
> topologicznym** — był ansatzem z `qm_spin`. Minimalizacja budżetu daje **h ≡ 1**
> (ψ płaskie, brak dołka). **Wyniki „liniowy profil / bag z ostrą krawędzią" oraz
> wyznaczenie ψ_core NIE OBOWIĄZUJĄ.** Zachowane wyłącznie jako zapis błędnej ścieżki.

**Data:** 2026-07-31
**Typ:** WYNIKI (domknięcie warstwy lokalizacji — od ograniczeń do KSZTAŁTU)
**Weryfikacja:** `BUDGET_profile.py` (uruchomione; obie formy metryki)
**Falsyfikator zadeklarowany przed liczeniem:** profil zdegenerowany (ψ≡1 / brak
skończonego rdzenia) ⟹ budżet nie wyznacza profilu. **Nie zadziałał.**

---

## 0. Zasada (jawne założenie, testowane)

> **Obiekt = konfiguracja MINIMALNEGO ZUŻYCIA BUDŻETU** spełniająca (a) dyskretność
> (|Δψ| ≤ δ_max na krok węzłowy) i (b) topologię (pełny trawers ψ_core→1).

To natywny odpowiednik dawnej „minimalizacji energii", ale na prymitywie, który
w TGP **faktycznie istnieje** (budżet/entropia). **To jest założenie, nie twierdzenie** —
i jest tym, co w tym rachunku można obalić.

## 1. Kształt: nie soliton, lecz „bag"

Minimalizacja daje jednoznacznie (ψ_core=0.7, δ_max=0.1):

```
i :  ψ_i     |Δψ|
0 : 0.7000   —
1 : 0.8000   0.1000  ← tempo NASYCONE
2 : 0.9000   0.1000  ← tempo NASYCONE
3 : 1.0000   0.1000  ← tempo NASYCONE
4+: 1.0000   0.0000  (płaska próżnia)
```

**Ten sam kształt w obu formach metryki** (forma I: h=ψ; M9.1'': h=ψ/(4−3ψ)) —
3 kroki nasycone w obu. Kontrola formy przeszła.

> **Profil = liniowy trawers z maksymalnym dozwolonym tempem + ostre przejście do próżni.**
> **NIE jest to gładki soliton z ogonem eksponencjalnym.** Obiekt typu **bag/kropla**,
> nie „chmura". To jest **predykcja strukturalna, która może zawieść** (ogony są mierzalne).

Test degeneracji: bez wymogu trawersu (ψ_core=1) obiekt **znika**, koszt = 0 — poprawnie.
Z trawersem: skończony rdzeń. Falsyfikator nie zadziałał.

## 2. Domknięcie głębokości (dotąd wolny parametr)

Dwa **przeciwstawne** więzy:
- **(a)** minimalizacja budżetu chce ψ_core → 1 (płytko = tanio),
- **(b)** koszt topologiczny wymaga N_rdzeń ≥ 5B węzłów (`WYNIKI_koszt-topologiczny`).

Przy tempie nasyconym r_c = (1−ψ_core)/δ_max, więc N = (4/3)πr_c³ ≥ 5B.
Minimalizacja **saturuje** ten wiąz (równość):

$$r_c(B)=\Big(\frac{15B}{4\pi}\Big)^{1/3},\qquad \psi_{core}(B)=1-\delta_{max}\,r_c(B)$$

| B | r_c [węzły] | ψ_core | N_rdzeń |
|---|---|---|---|
| 1 | 1.061 | 0.8939 | 5.0 |
| 2 | 1.337 | 0.8663 | 10.0 |
| 3 | 1.530 | 0.8470 | 15.0 |
| 8 | 2.122 | 0.7878 | 40.0 |
| 27 | 3.182 | 0.6818 | 135.0 |

> **Głębokość i rozmiar wyznaczone JEDNOCZEŚNIE, oba ∝ B^{1/3}.**
> Żaden nie jest wolnym parametrem — **jedyną wolną skalą pozostaje δ_max**
> (rozdzielczość substratu, czyli fundamentalna długość; por. aksjomat A7 ℓ_P=const).

## 3. Bilans warstwy lokalizacji (kompletny)

| Pytanie | Odpowiedź | Mechanizm |
|---|---|---|
| dlaczego nie kolapsuje? | cap budżetowy | `WYNIKI_budzet-skala` |
| dlaczego nie rozpływa się? | przeszkoda topologiczna | `WYNIKI_domykanie-solitony` |
| jak duży? | R ∝ B^{1/3} | koszt N≥5B + budżet |
| jak głęboki? | (1−ψ_core) ∝ B^{1/3} | saturacja (a)⟂(b) |
| jaki kształt? | liniowy trawers + ostra krawędź | min. budżetu + dyskretność |

**Nic z tego nie używa: energii, członu Skyrme'a, płaskiego tła, fitu.**

## 4. Predykcje falsyfikowalne (do przyszłej konfrontacji)

1. **Brak ogonów eksponencjalnych** — ostra krawędź, profil liniowy w rdzeniu.
2. **R ∝ B^{1/3}** oraz **głębokość ∝ B^{1/3}** — gęstość nasycenia.
3. **B=1 = obiekt minimalny** (5 węzłów) — fundamentalny fermion najtańszy.

## 5. Ograniczenia (uczciwie)

1. **Zasada minimalizacji budżetu jest ZAŁOŻENIEM** — natywnym, ale nie wyprowadzonym z rdzenia.
2. Liniowość profilu wynika dość **bezpośrednio** z (koszt wypukły + limit tempa) — to nie
   jest głęboka niespodzianka, lecz konsekwencja setupu. Nietrywialna jest **ostra krawędź**
   i brak ogona (to odróżnia od standardowych solitonów).
3. **δ_max wolne** — cała skala bezwzględna nieustalona; wszystkie wyniki są w jednostkach δ_max.
   **Nie fitować δ_max do niczego znanego.**
4. Stała „5" (z bounda symplicjalnego) **zależy od dyskretyzacji** — nie jest predykcją liczbową.
5. Idealizacja radialna/sferyczna; **sprzężenie z orientacją reperu (sektor 2T) nie policzone** —
   liczyłem sektor ψ, zakładając trawers z `qm_spin` (f slaved to ψ). Pełne sprzężenie
   ψ ↔ orientacja pozostaje OTWARTE.

## 6. Jednozdaniowo

> Minimalizacja **budżetu** (nie energii) przy dyskretności substratu i wymogu trawersu
> topologicznego daje **obiekt typu bag**: liniowy rdzeń, ostra krawędź, **bez ogonów**,
> o rozmiarze i głębokości **wyznaczonych jednocześnie i obu ∝ B^{1/3}**, z jedyną wolną
> skalą δ_max — co domyka warstwę lokalizacji od strony kształtu, nie tylko ograniczeń.

---

**Plik:** `BUDGET_profile.py`
**Powiązane:** `WYNIKI_koszt-topologiczny_2026-07-31.md`, `WYNIKI_budzet-skala_2026-07-31.md`,
`ONTOLOGIA_energia-relacyjna-budzet_2026-07-31.md`, `research/qm_spin/README.md`
