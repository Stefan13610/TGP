# Audyt adwersarialny bloku budżetowego — WYNIK: blok się sypie

**Data:** 2026-07-31
**Typ:** AUDYT NIEZALEŻNY (dwóch audytorów, soczewki: fundamenty / rachunki)
**Metoda:** własne skrypty audytorów (`AUDIT4_*.py`), niezależne konstrukcje (SL(2,3)/GF(3)),
uruchamianie kodu autora u źródła
**Werdykt:** **warstwa budżetowa upada w całości.** Ocalała wyłącznie warstwa grupowa 2T.

---

## 0. Werdykt audytora (cytat)

> „Korekta autora poszła **za mało daleko i w złą stronę**: wycofał dwa wyniki na podstawie
> rachunku (h≡1), który jest tożsamością własnej konstrukcji — a nie wycofał
> `WYNIKI_budzet-skala`, którego centralny mechanizm jego własna korekta unieważnia
> doszczętnie (przy h≡1 zachowana ilość Q=∫(ψ−1)dV = 0, więc R_min = 0)."

## 1. Ustalenia krytyczne

| # | Ustalenie | Waga |
|---|---|---|
| 1 | **h≡1 to artefakt konstrukcji.** W funkcji celu h wchodzi tylko przez g(h) (min w h=1) i **nie ma członu sprzęgającego h z orientacją** — model nie mógł dać nic innego. 0/262 losowych profili f dało h≠1. | **KRYTYCZNE** |
| 2 | **BUG:** h[0] to zmienna martwa (waga w₀=0²=0, więz tylko na h[1:]); h₀ przyjmował wartość startową (x0[0]=3.0 → h₀=3.0). | **KRYTYCZNE** |
| 3 | **Przy konsekwentnym sprzężeniu** (|df| w odległości własnej, bo g=h·δ) rdzeń **dołuje w drugą stronę**: h≈1.26–1.50 **>** 1. ⟹ pytanie o dołek jest OTWARTE. | **KRYTYCZNE** |
| 4 | **Lokalizacja = artefakt UV.** Obiekt ma zawsze ~π/δ_max ≈ 3 komórki siatki ⟹ R_phys/L = 0.300, 0.150, …, **0.0047** dla M=10…640. Rozmiar ∝ stała sieci, **znika w granicy kontinuum**. **To ten sam wzorzec, co dawne „skalowanie z pudłem"**, przesunięty na koniec UV. | **KRYTYCZNE** |
| 5 | **Samozniszczenie `WYNIKI_budzet-skala`:** h≡1 ⟹ ψ≡1 ⟹ Q=0 ⟹ R_min=0. Autor tego nie zauważył i nie wycofał dokumentu. | **KRYTYCZNE** |
| 6 | **Nadokreślenie w RDZENIU:** n_sp+n_czas=B **oraz** n_sp·n_czas=const dają x²−2x+1 ⟹ **h≡1 identycznie** — Φ nie może się zmieniać. Oba wiązy naraz zamrażają pole. | **WYSOKIE (dotyczy rdzenia)** |
| 7 | **n_orient łamie przesłankę** prop:antipodal („budżet dzielony **wyłącznie** między h i f, brak ukrytych stopni swobody"), z której pochodzi fh=1 → g(h)=h+1/h. **Sprzeczność wewnętrzna.** | **WYSOKIE** |
| 8 | „δ_max wypada z budżetu" = **reparametryzacja** 1→2 parametry (B, c); (B,c)=(3,1/δ) odtwarza dowolne δ. **Nie jest zyskiem.** | **WYSOKIE** |
| 9 | **R_min ∝ Q^{1/3} = tautologia** (cap + zachowana ilość dają to w każdej teorii; identycznie dla 8 losowych kształtów). | **WYSOKIE** |
| 10 | **„Próżnia = minimum budżetu" = artefakt normalizacji** n₀=n₀ₜ. Dla n₀=3,n₀ₜ=1 minimum przy h=0.577. | **WYSOKIE** |
| 11 | **Bound symplicjalny NIE ogranicza węzłów.** T ≤ C(V,4) ⟹ V ~ B^{1/4}, R ~ B^{1/12}. Most „topologia kosztuje BUDŻET (węzły)" **nie przechodzi** — konflacja sympleksów z węzłami. | **WYSOKIE** |
| 12 | **„Weryfikacja χ(S³)=0" jest pusta** — χ=0 dla każdej triangulacji rozmaitości nieparzystowymiarowej. Nic nie zweryfikowała. | **ŚREDNIE** |
| 13 | **Oba „falsyfikatory zadeklarowane przed liczeniem" to pseudo-falsyfikatory** o zerowym prawdopodobieństwie zadziałania: „cap istnieje" ⟸ AM-GM (h+1/h≥2); „skończona triangulacja" ⟸ zwartość PL. Testowały fakty znane a priori, nie hipotezy o TGP. | **WYSOKIE (metodologiczne)** |
| 14 | **Domykanie: „trialność-0 nie wystarcza (21/30)" = artefakt** wliczenia **klas centralnych** ({e},{−1}) jako „wstążek". Po ich usunięciu: **13/13 domyka się**. Teza „mocniejsze niż SM" upada. | **WYSOKIE** |

## 1b. Ustalenia audytora #1 (soczewka: fundamenty — czytał rdzeń u źródła)

| # | Ustalenie | Waga |
|---|---|---|
| A | **GŁÓWNY PRZEMYT: B jest przez rdzeń USTALONE na 2, a autor je PODNIÓSŁ.** Wiąz addytywny obowiązuje jako **równość** także w Φ=Φ₀, gdzie n_sp=n₀, n_czas=n₀ₜ ⟹ **B = n₀+n₀ₜ**. Cała tabela §1 `WYNIKI_budzet-skala` (B=2.5…100) wymaga **nadwyżki B>2, której rdzeń zabrania**. Dyskryminant B²−4n₀n₀ₜ>0 jest niezerowy **wyłącznie dzięki wstrzykniętej nadwyżce**. Bez niej cap i limit tempa **znikają identycznie**. | **KRYTYCZNE** |
| B | Przy kanonicznym **B=2**: |Δf|_max=(B−2)/c = **0** ⟹ **żadne nawinięcie nie jest reprezentowalne, obiekt nie istnieje.** | **KRYTYCZNE** |
| C | **Reaktywowałem addytywność, którą rdzeń JAWNIE ODRZUCA** w tym samym akapicie: *„Sam wiąz addytywny nie prowadzi do fh=1. Wymagamy silniejszego warunku: **entropia substratowa jest nieaddytywna**"*. Suma n_sp+n_czas jako „zużycie" jest dokładnie tym, co rdzeń unieważnia. | **KRYTYCZNE** |
| D | **n_orient łamie przesłankę, z której bierze fh=1** — i to w jednym równaniu. A fh=1 jest założeniem (ii) thm:metric-from-budget-M911 ⟹ **rozszerzenie unieważnia metrykę, na której liczę.** | **KRYTYCZNE** |
| E | **Liniowość kosztu kąta jest NOŚNA.** Przy koszcie **kwadratowym** (co dałaby jakakolwiek akcja gradientowa) optimum to Δf ∝ 1/wᵢ — **rozmyte po całej domenie, bez nasycenia i bez skończonego rdzenia**. Lokalizacja wisi wyłącznie na arbitralnym wyborze liniowości. | **KRYTYCZNE** |
| F | „Trawers skupiony" to rozwiązanie **bang-bang programu liniowego** (rosnąca waga i² + wiąz pudełkowy) — **wymuszone algebraicznie, nie emergentne**. | **WYSOKIE** |
| G | **Bilans parametrów:** usunąłem 1 (δ_max), wprowadziłem 2 (B, c) ⟹ **netto STRATA sprzedana jako „realny zysk"**. | **WYSOKIE** |
| H | **Rdzeń NIGDZIE nie minimalizuje budżetu** — budżet występuje wyłącznie jako **wiąz zachowawczy** (B=const). „Minimalizacja budżetu" nie ma kotwicy w rdzeniu ⟹ **ten sam błąd co minimalizacja energii, w nowym przebraniu.** | **KRYTYCZNE (metodologiczne)** |
| I | **ψ_max→4/3 NIE jest fitem** — kod nie zawiera g₀, PDG ani skanu do celu; **zarzut cyrkularności OBALONY, tu jestem czysty**. Ale wynik jest **PUSTY**: ψ=4h/(1+3h)<4/3 dla każdego h>0 — własność parametryzacji M9.1'', budżet nie wnosi nic. | **ŚREDNIE** |
| J | **Osobny problem RDZENIA:** most g₀,crit=1,874 ≡ ψ_horizon=4/3 żyje w `sek08a` lin.994–996 jako „koincydencja na 4 cyfrach", mimo że **1,874 ≠ 1,333** (wymaga nieujawnionego w tym miejscu mapowania g₀↔ψ). | **WYSOKIE (rdzeń)** |

**Werdykt audytora #1:** *„Warstwa budżetowa sypie się… obie ilościowe konsekwencje (cap amplitudy, limit tempa) są proporcjonalne do B−(n₀+n₀ₜ), a rdzeń zeruje tę wielkość."*

**Uczciwość intencji (cytat):** *„Żadnego dopasowania do PDG, ℓ_P, gęstości nasycenia ani g₀. Falsyfikatory deklarowane przed liczeniem. Autokorekty prawdziwe i kosztowne dla autora. Pod względem intencji autor jest uczciwy — problemem jest, że oflagował **słabszą wersję każdego zarzutu**, niż wynosi wersja rzeczywista."*

**Droga ratunku wskazana przez audytora:** porzucić addytywną sumę i wyprowadzić miarę zużycia
budżetu z **rzeczywistej, nieaddytywnej struktury entropii substratowej**, na którą rdzeń się
powołuje, ale **której nie podaje**. Dopóki ta miara nie istnieje, każdy „cap z budżetu"
odtworzy ten sam błąd.

## 2. Co OCALAŁO (niezależnie odtworzone)

- **Struktura 2T** (odtworzona konstrukcją macierzową SL(2,3), niezależnie od kwaternionów autora):
  |2T|=24, [2T,2T]=Q₈, **−1∈Q₈ ⟹ spin bezbarwny**, trialność = homomorfizm na ℤ₃ (576/576 par),
  7 klas, 5 par domykających, **statystyka złożenia (e vs −1) NIEWYMUSZONA** (4/4 trójki).
  Audytor: *„Ta warstwa jest nietknięta i uczciwie zaraportowana"*.
- **Uniqueness 2T** (z wcześniejszego audytu, odporne na cutoff).
- **Bound symplicjalny T_dom ≥ 5·|B|** — poprawne twierdzenie **o odwzorowaniach symplicjalnych**
  („5" użyte w bezpieczną stronę), ale **NIE o węzłach ani o budżecie**.
- **g_ij = h·δ_ij niezależnie od R** — poprawne, ale trywialne (przepisanie RᵀR=I).
- „Pojedyncza wstążka nie anihiluje" i „trialność-0 konieczna" — prawdziwe, ale **tautologiczne**
  (g≠e; tri(e)=0). Działają jako przeszkoda topologiczna, ale jako „wyniki" mają zerową treść.

## 3. Wzorzec do zapamiętania (najważniejsza lekcja)

> Audytor: *„Głębszy wzorzec jest ten sam co w dawnym »skalowaniu z rozmiarem pudła«, tylko
> przesunięty na koniec UV: obiekt ma zawsze ~3 komórki siatki, więc jego rozmiar fizyczny jest
> proporcjonalny do stałej sieci i znika w granicy kontinuum. To nie jest skala dynamiczna —
> to obcięcie."*

**Dodatkowy wzorzec (mój, autora):** trzykrotnie w tej sesji zbudowałem „wynik" na kryterium
zaimportowanym z zewnątrz (płaskie tło → absolutna energia → minimalizacja budżetu). Za każdym
razem audyt to złapał. **Zamiana nazwy prymitywu (energia→budżet) nie czyni zasady wariacyjnej
natywną**, jeśli rdzeń jej nie uzasadnia.

## 4. Status dokumentów po audycie

| Dokument | Status |
|---|---|
| `WYNIKI_budzet-skala_2026-07-31.md` | ⛔ **WYCOFANY w całości** |
| `WYNIKI_profil-z-budzetu_2026-07-31.md` | ⛔ **WYCOFANY** (wycofany wcześniej, ale z błędnego powodu) |
| `WYNIKI_psi-orientacja_KOREKTA_2026-07-31.md` | ⛔ **sama korekta WYCOFANA** (pusta podstawa + bug) |
| `WYNIKI_koszt-topologiczny_2026-07-31.md` | ⚠ **zawężony do twierdzenia o sympleksach**; most do budżetu/węzłów OBALONY |
| `WYNIKI_domykanie-solitony_2026-07-31.md` | ⚠ §2(a) „21/30" **OBALONE** (klasy centralne); reszta = tautologie |
| `ONTOLOGIA_energia-relacyjna-budzet_2026-07-31.md` | ⚠ ontologia (energia relacyjna) **stoi**; ale §5 hipoteza budżetowa **upadła** |
| `WYNIKI_s54_domkniecie_frame-native`, `WYNIKI_substrat-wstazki_2T` | ✅ nietknięte |

## 5. Stan wieży po audycie (uczciwie)

```
substrat  → wstążki      → domykanie   → lokalizacja
  (?)        2T ✅          tautologie      ✗ UPADŁA
```

- **Ocalałe DERIVED:** uniqueness 2T; spin bezbarwny (−1∈Q₈); struktura Q₈⋊ℤ₃; bound symplicjalny
  (jako czyste twierdzenie).
- **Warstwa lokalizacji: BRAK.** Nie mamy mechanizmu skali. Wszystkie trzy próby (Skyrme, cap
  budżetowy, limit tempa) upadły — każda z innego powodu.
- **Otwarte i nierozwiązane:** co lokalizuje obiekt; czy rdzeń ψ dołuje; nadokreślenie wiązów
  w samym rdzeniu (pkt 1.6) — **to wymaga rozstrzygnięcia w rdzeniu, nie w tej gałęzi.**

---

**Skrypty audytorów:** `AUDIT4_h_decoupling.py`, `AUDIT4_localization.py`, `AUDIT4_simplicial.py`,
`AUDIT4_scale_cap.py`, `AUDIT4_core_consistency.py`, `AUDIT4_ribbons_independent.py`
