# Skala z budżetu: cap na Φ blokuje kolaps — lokalizacja bez energii

> # ⛔ WYCOFANE 2026-07-31 (audyt adwersarialny) — patrz `AUDYT_blok-budzetowy_2026-07-31.md`
> 1. **SAMOZNIŚCZENIE:** własna korekta autora (h≡1) daje ψ≡1 ⟹ Q=∫(ψ−1)dV = **0** ⟹ **R_min = 0**.
>    Mechanizm §3 jest po korekcie PUSTY, nie „zastąpiony".
> 2. **R_min ∝ Q^{1/3} to TAUTOLOGIA** — cap na gęstość + zachowana ilość dają to w KAŻDEJ teorii
>    (audytor: identycznie dla 8 losowych kształtów). Zero treści specyficznej dla TGP.
> 3. **„Próżnia = minimum budżetu" to ARTEFAKT NORMALIZACJI** n₀=n₀ₜ. Dla n₀=3, n₀ₜ=1 minimum
>    wypada przy h=0.577 — próżnia NIE jest minimum.
> 4. **„ψ_max<4/3" to algebra**, nie budżet: 4h/(1+3h)<4/3 dla każdego h>0.
> 5. **Falsyfikator był pseudo-falsyfikatorem** — „cap istnieje" wynika z AM-GM (h+1/h≥2), było
>    gwarantowane a priori. Nie mógł zadziałać.
> **Cały dokument nieważny.** Zachowany jako zapis błędnej ścieżki.

**Data:** 2026-07-31
**Typ:** WYNIKI (warstwa lokalizacji, ścieżka natywna: budżet zamiast energii)
**Weryfikacja:** `BUDGET_scale_cap.py` (uruchomione)
**Falsyfikator zadeklarowany PRZED liczeniem:** brak górnego ograniczenia zużycia budżetu
przy fh=1 ⟹ brak capu ⟹ hipoteza pada. **Nie zadziałał** — cap istnieje.

---

## 0. Podstawa (cytat z rdzenia, `sek08c` def:info-budget)

> **B = N_B · s₀**, s₀ = m_sp² = γ (entropia mikroskopowa na węzeł).
> „**B jest stałą** — nie zależy od Φ(x), ponieważ **N_B jest liczbą topologiczną**."
> Budżet dzielony: n_sp + n_czas = B; h = n_sp/n₀, f = n_czas/n₀ₜ;
> założenie multiplikatywne n_sp·n_czas = const ⟹ **fh = 1**.

**Zastrzeżenie uczciwości:** rdzeń sam zaznacza, że wiąz addytywny NIE daje fh=1 i że
multiplikatywność jest **dodatkowym** założeniem. Używam **obu naraz** (fh=1 oraz
n_sp+n_czas ≤ B). To jest **moja interpretacja**, nie cytat z twierdzenia.

## 1. Cap istnieje (wynik główny)

Zużycie budżetu przy fh=1: **g(h) = n₀h + n₀ₜ/h**.

- **Minimum w h = 1** (g_min = 2 w normalizacji n₀=n₀ₜ=1),
- **g(h) → ∞ gdy h→0 ORAZ h→∞** — zużycie rośnie w obie strony.

⟹ dla skończonego budżetu B istnieje **skończony przedział dozwolony [h_min, h_max]**:

| B | h_min | h_max | ψ_max (forma I) | ψ_max (M9.1'') |
|---|---|---|---|---|
| 2.5 | 0.500 | 2.000 | 2.000 | 1.1429 |
| 3.0 | 0.382 | 2.618 | 2.618 | 1.1827 |
| 4.0 | 0.268 | 3.732 | 3.732 | 1.2240 |
| 10 | 0.101 | 9.899 | 9.899 | 1.2899 |
| 100 | 0.010 | 99.99 | 99.99 | 1.3289 |

**Pole jest ograniczone OBUSTRONNIE** — nie tylko z góry. W M9.1'' **ψ_max → 4/3 przy B→∞**,
zawsze **ściśle poniżej** horyzontu; ψ=4/3 to granica absolutna (cały budżet na przestrzeń,
zegar staje, f=0). *(NIE używam tu obalonego mostu ψ↔g₀ — to fit 2-punktowy, inny argument.)*

## 2. Efekt uboczny o dużej wadze: próżnia = minimum budżetu

**Minimum zużycia budżetu wypada dokładnie przy h=1, czyli ψ=1 — czyli w próżni.**

To **natywnie tłumaczy** poprzedni wynik (`WYNIKI_cisnienie-lokalizacja`), gdzie wyszło,
że wszystko relaksuje do ψ≡1 — ale **bez odwoływania się do energii**:

> Układ dąży do ψ=1 nie dlatego, że to minimum energii, lecz dlatego, że to **stan
> najmniejszego zużycia budżetu**. „Relaksacja do próżni" = **ekonomia budżetu**, nie
> minimalizacja energii.

Zgodne z ontologią §1 (energia nie jest prymitywem): ten sam fakt wyprowadzony z prymitywu,
którym jest budżet/entropia.

## 3. Cap blokuje kolaps Derricka (mechanizm skali)

Obiekt: ψ(r) = 1 + A·χ(r/R). Zachowana „ilość substratu" Q ~ ∫(ψ−1)dV = c·A·R³.
Kolaps R→0 przy stałym Q wymaga **A → ∞**. Ale budżet narzuca **A ≤ A_max = ψ_max − 1**:

$$\frac{Q}{cR^3}\le A_{max}\;\Longrightarrow\; R \ge R_{min}=\Big(\frac{Q}{c\,A_{max}}\Big)^{1/3}>0$$

| ψ_max | A_max | R_min (Q=1) | R_min (Q=8) |
|---|---|---|---|
| 1.10 | 0.100 | 2.154 | 4.309 |
| 1.20 | 0.200 | 1.710 | 3.420 |
| 4/3 | 0.333 | 1.442 | 2.885 |

> **KOLAPS ZABLOKOWANY.** Istnieje skończony minimalny rozmiar R_min > 0.
> **Zero minimalizacji energii. Zero członu Skyrme'a.** Skala z więzu budżetowego.

To jest dokładnie to, czego brakowało po dyskwalifikacji płaskiego Skyrme'a (§5.4).

## 4. Predykcja strukturalna

**R_min ∝ Q^{1/3}** ⟺ **stała gęstość maksymalna = gęstość nasycenia**
(Q=1→R=1, Q=8→R=2, Q=27→R=3). Podpis materii nasyconej — jak R ∝ A^{1/3} dla jąder.

**Uczciwie:** samo skalowanie Q^{1/3} jest **generyczne dla dowolnego capu**; natywne
i nietrywialne jest tu **ŹRÓDŁO capu (budżet)**, nie sam wykładnik. Nie ogłaszam
„wyprowadziliśmy prawo jądrowe".

## 5. Domknięcie obustronne (synteza z sektorem topologicznym)

| Kierunek destrukcji | Co blokuje | Status |
|---|---|---|
| **kolaps** (R→0) | **cap budżetowy** na amplitudę | ✅ ten dokument |
| **rozpłynięcie** (R→∞, A→0, powrót do próżni) | **przeszkoda topologiczna** — wstążka nie domyka się sama | ✅ `WYNIKI_domykanie-solitony` |

⟹ **rozmiar obiektu jest ograniczony OBUSTRONNIE: budżet od dołu, topologia od góry.**
Oba więzy są natywne, relacyjne, i **żaden nie używa minimalizacji energii**.

## 6. Co pozostaje OTWARTE (nie zamiatam)

1. **Co dokładnie jest zachowane (Q)?** Przyjąłem ∫(ψ−1)dV. W sektorze wstążek powinien to być
   **ładunek topologiczny**; związek Q ↔ nawinięcie **nie jest ustalony**. To główna luka.
2. **Podwójne założenie budżetowe** (addytywne ≤B + multiplikatywne fh=1) — interpretacja, nie twierdzenie.
3. **Wartość B** (a więc ψ_max) nieustalona — nie fitować do niczego znanego.
4. Kształt χ przyjęty sztywny (skalowanie jednoparametrowe) — pełny profil nie liczony.

## 7. Jednozdaniowo

> Budżet informacyjny substratu daje **cap na Φ** (zużycie rośnie w obie strony od ψ=1,
> więc skończone B ⟹ skończony przedział), co **blokuje kolaps** i daje **R_min ∝ Q^{1/3}**
> — skala bez energii i bez Skyrme'a; w połączeniu z przeszkodą topologiczną (wstążka nie
> domyka się sama) rozmiar jest ograniczony obustronnie, a próżnia ψ=1 okazuje się
> **stanem minimalnego budżetu**, co natywnie tłumaczy relaksację bez odwołania do energii.

---

**Plik:** `BUDGET_scale_cap.py`
**Powiązane:** `ONTOLOGIA_energia-relacyjna-budzet_2026-07-31.md` (ustalenie ontologiczne),
`WYNIKI_domykanie-solitony_2026-07-31.md` (górne ograniczenie), `core/sek08c` (def:info-budget)
