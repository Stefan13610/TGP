# Sprzężenie ψ↔orientacja — domknięcie, które KORYGUJE dwa wcześniejsze wyniki

> # ⛔ SAMA KOREKTA WYCOFANA 2026-07-31 (audyt) — patrz `AUDYT_blok-budzetowy_2026-07-31.md`
> Audyt wykazał, że **korekta była zrobiona na pustej podstawie**:
> 1. **h≡1 to ARTEFAKT KONSTRUKCJI, nie wynik.** W funkcji celu h wchodzi WYŁĄCZNIE przez g(h)
>    (min w h=1) i **nie ma ANI JEDNEGO członu sprzęgającego h z orientacją**. Model nie mógł
>    dać nic innego. Audytor: 0/262 losowych profili f dało h≠1.
> 2. **BUG:** h[0] jest zmienną martwą (waga w₀=0²=0, więz tylko na h[1:]) — h₀ przyjmował
>    wartość startową (x0[0]=3.0 → h₀=3.0).
> 3. **Przy konsekwentnym sprzężeniu** (|df| mierzone w odległości własnej, bo g=h·δ) rdzeń
>    **dołuje W DRUGĄ STRONĘ**: h≈1.26–1.50 > 1.
> ⟹ Pytanie „czy rdzeń ψ dołuje" jest **OTWARTE**, a nie rozstrzygnięte na „nie".
> Wycofanie profilu „bag" było zrobione **z błędnego powodu** (choć sam profil i tak nie jest ustalony).
> 4. „δ_max wypada z budżetu" = **REPARAMETRYZACJA** (1 parametr → 2: B i c), a dodatkowo
>    n_orient **łamie przesłankę prop:antipodal** („dzielony WYŁĄCZNIE między h i f"), z której
>    pochodzi fh=1 → g(h)=h+1/h. **Sprzeczność wewnętrzna.**

**Data:** 2026-07-31
**Typ:** WYNIKI + **AUTOKOREKTA** (domknięcie ostatniej luki warstwy lokalizacji)
**Weryfikacja:** `BUDGET_psi_orientation.py` (uruchomione)

---

## 0. Ustalenie decydujące: metryka NIE WIDZI orientacji

Reper e^a_i = √h·R^a_i, R ortogonalna ⟹ g_ij = e^a_i e^a_j = h·(RᵀR)_ij = **h·δ_ij**.

Zweryfikowane numerycznie: max|g − h·I| ≈ 6·10⁻¹⁶ dla losowych R ∈ SO(3).

> **Orientacja R jest niewidoczna w metryce.** Zatem **topologia (nawinięcie) siedzi
> w ORIENTACJI, a nie w profilu ψ.**

## 1. ⚠ KOREKTA #1: wymóg „trawersu ψ" był NIEUPRAWNIONY

W `WYNIKI_profil-z-budzetu_2026-07-31.md` narzuciłem, że topologia wymaga trawersu
ψ_core → 1, powołując się na `qm_spin` (f = π(1−g)/(1−g₀)). **To była parametryzacja
(ansatz), nie wymóg topologiczny** — sam `qm_spin` stwierdza, że B=1 jest „niezależne
od g₀ i kształtu profilu".

**Test:** minimalizacja budżetu przy obecnym nawinięciu daje **h ≡ 1 na całym profilu**
(max|h−1| = 0.000000).

> **Rdzeń ψ NIE dołuje. Obiekt jest CZYSTĄ TEKSTURĄ ORIENTACJI o płaskiej metryce.**
> Wynik „liniowy profil ψ / bag z ostrą krawędzią" **NIE OBOWIĄZUJE** — opierał się
> na wymuszonym trawersie ψ, którego topologia nie wymaga.

## 2. ✅ ZYSK: limit tempa WYPADA z budżetu (był założeniem)

Rozszerzenie budżetu o orientację (naturalne: s₀ = entropia na węzeł = log liczby stanów,
a stany węzła obejmują orientację):

$$n_{sp} + n_{czas} + n_{orient} \le B,\qquad g(h)=h+1/h \ge 2$$

⟹ **|Δf|_max = (B − 2)/c** — maksymalna zmiana orientacji na krok węzłowy.

| B | \|Δf\|_max | kroki na trawers π |
|---|---|---|
| 2.2 | 0.200 | 15.7 |
| 2.5 | 0.500 | 6.3 |
| 3.0 | 1.000 | 3.1 |

> **δ_max przestaje być wolnym parametrem** — wypada z **nadwyżki budżetu (B−2)**.
> To usuwa jedno z założeń poprzedniego rachunku. Realny zysk.

## 3. Lokalizacja działa — ale innym mechanizmem

Trawers okazał się **SKUPIONY** (4 z 20 kroków), z tempem **nasyconym** na (B−2)/c.

**Mechanizm:** miara powłokowa (w ∝ i²) czyni nawinięcie **najtańszym w centrum**
(mało węzłów). Więc winding koncentruje się w rdzeniu i biegnie z maksymalnym tempem.

⟹ **rozmiar rdzenia = (trawers)/|Δf|_max** — skończony, wyznaczony przez budżet.

**Kolaps nadal zablokowany**, ale przez **limit tempa** (wyprowadzony), a nie przez
**cap amplitudy** (który przy h≡1 w ogóle się nie aktywuje). Mechanizm czystszy niż poprzednio.

## 4. ⚠ KOREKTA #2: skalowanie R(B) jest NIEPEWNE

`WYNIKI_koszt-topologiczny` dało **R ∝ B^{1/3}** z bounda objętościowego N ≥ 5B.
Ale w ansatzu hedgehog nawinięcie B wymaga trawersu **Bπ**, więc przy stałym tempie:
**R ∝ B (liniowo)**, nie B^{1/3}.

- Bound objętościowy N≥5B daje R ≳ B^{1/3} — **słabszy, niewiążący**.
- W ansatzu sferycznym wiążący jest trawers radialny ⟹ **R ∝ B**.
- **ALE:** dla B>1 sferyczny hedgehog **nie jest konfiguracją minimalną** (w modelach
  solitonowych stany B>1 są powłokowe/wielościenne, nie sferyczne).

> **Werdykt: B=1 jest solidne; skalowanie dla B>1 jest OTWARTE i zależne od ansatzu.**
> Wcześniejsze „R ∝ B^{1/3} = gęstość nasycenia" **wycofuję jako nieustalone.**

## 5. Co zostaje po korektach (stan uczciwy)

**PRZEŻYWA:**
- Topologia w orientacji; metryka nie widzi R (zweryfikowane).
- **Limit tempa z budżetu** — wyprowadzony, nie założony (zysk netto).
- **Lokalizacja**: winding skupia się w rdzeniu (miara powłokowa) i nasyca tempo ⟹ skończony rozmiar.
- Kolaps zablokowany (przez tempo, nie cap).
- Przeszkoda topologiczna (wstążka nie domyka się sama) — nietknięta.
- 2T / spin+kolor / reguły domykania — nietknięte (to inna warstwa).

**WYCOFANE / SKORYGOWANE:**
- „Profil ψ liniowy, bag z ostrą krawędzią" — **wycofane** (opierał się na złym wymogu).
- „R ∝ B^{1/3}, gęstość nasycenia" — **wycofane jako nieustalone** (ansatz-zależne).
- „Cap amplitudy blokuje kolaps" — **zastąpione** przez limit tempa (cap nieaktywny przy h≡1).

**OTWARTE:**
- Skalowanie R(B) dla B>1 (konfiguracje niesferyczne).
- Czy h≡1 przeżyje, gdy budżet orientacji sprzęgnie się z sektorem 2T (dysklinacje,
  nie tylko hedgehog) — liczyłem hedgehog radialny.
- Stała c (koszt orientacji na jednostkę kąta) — wolna.

## 6. Jednozdaniowo

> Domknięcie sprzężenia pokazało, że **metryka nie widzi orientacji**, więc topologia
> mieszka w R, a nie w ψ — co **unieważnia mój wymóg trawersu ψ i wycofuje wynik
> „bag" oraz skalowanie B^{1/3}** — ale jednocześnie **wyprowadza limit tempa z budżetu**
> (usuwając założenie) i pokazuje, że obiekt lokalizuje się jako **czysta tekstura
> orientacji o płaskiej metryce**, z rozmiarem wyznaczonym przez nadwyżkę budżetu.

---

**Plik:** `BUDGET_psi_orientation.py`
**Koryguje:** `WYNIKI_profil-z-budzetu_2026-07-31.md`, `WYNIKI_koszt-topologiczny_2026-07-31.md`
**Nietknięte:** `WYNIKI_domykanie-solitony`, `WYNIKI_s54_domkniecie_frame-native`
