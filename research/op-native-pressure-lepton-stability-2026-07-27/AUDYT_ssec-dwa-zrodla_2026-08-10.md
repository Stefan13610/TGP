# Audyt `ssec:dwa-zrodla` (sek08_formalizm, lin. 2578–2916) — fundament trzech reżimów

**Data:** 2026-08-10
**Typ:** AUDYT ŹRÓDŁOWY + weryfikacja numeryczna całek nakładania
**Skrypty:** `REZIMY_audit_dwa-zrodla.py` (v1 — **skażony pudłem**), `REZIMY_audit_dwa-zrodla_v2.py` (poprawiony)
**Kontekst:** to jest sekcja, do której §3 odsyła jako do „zamknięcia luki formalnej".

---

## 0. Werdykt

> **Reżim I (grawitacja) — STOI.** **Reżim II (odpychanie) — znak stoi, wzór na skalowanie NIE.**
> **Reżim III (studnia/confinement) — NIE JEST USTANOWIONY.** Wyprowadzony z funkcjonału
> w obszarze, w którym **rdzeń sam deklaruje go za niewiarygodny**, przy pomocy dwóch
> **wzajemnie sprzecznych** wzorów na skale przejścia, z których jeden zależy od członu
> **nigdy nieobliczonego**.
> **Dobra wiadomość:** stabilne minimum (odległość równowagowa pary) **przeżywa usunięcie
> całej wątpliwej części** — patrz §5.

## 1. Najpierw: mój własny błąd w tym audycie

Pierwsza wersja skryptu liczyła ∫Φ₁Φ₂ dla Φ ~ A/(r+r₀) **bez ekranowania**. Wyszło
I₂ ≈ 828 niemal niezależnie od d — i o mało nie zaraportowałem tego jako „I₂ ~ d^(−0.045)".
**To był artefakt pudła**: całka jest IR-rozbieżna, więc mierzyłem objętość siatki, nie fizykę.
Złapałem to i przeliczyłem z ekranowaniem.

Drugi błąd: w wydruku v2 napisałem „I₃/I₂ **maleje** gdy d maleje". **To nieprawda** —
rośnie (0.093 przy d=8 → 0.427 przy d=0.5). Kierunek jest **zgodny** z tezą rdzenia.
Prostuję poniżej w §4.

## 2. Co stoi (bez zastrzeżeń)

**(a) `rem:Elin-yukawa` — wzorcowo uczciwy fragment.** Profil dalekiego pola ma czynnik
e^(−√γ d); z γ ~ 10⁻⁵² m⁻² i r₀ ~ l_P wychodzi γ̂ ~ 10⁻¹²², więc korekcja
> „e^(−√γ̂ d) > 1 − 10⁻⁵¹ i jest **pomijalna**. Przybliżenie 1/d jest dokładne
> z precyzją **ponad 50 rzędów wielkości**."

Sprawdzone — rachunek się zgadza. **Reżim I (grawitacja, E_lin ~ −1/d) jest solidny.**

**(b) `prop:ghost-free-fundamental`** — zerowanie f(Φ)=1+4ln(Φ/Φ₀) przy ψ≈0.78 jest
**artefaktem rozwinięcia**, a w działaniu fundamentalnym człon ψ⁵(∇ψ)² jest dodatnio
określony na całym basenie ψ∈(0,4/3). **To jest poprawne rozwiązanie problemu ducha** —
ale ma cenę, patrz §3(a).

**(c) `rem:Veff-Eint-bridge` — uczciwa higiena epistemiczna.** Rdzeń **sam przyznaje**,
że V_eff w §3 był wprowadzony „**fenomenologicznie**", a ta sekcja domyka lukę. Dobry zwyczaj.

## 3. Cztery problemy

### (a) ⛔ Rdzeń deklaruje zakres ważności — i reżim III leży poza nim

Lin. 2660–2662, **dosłownie**:
> „Funkcjonał E_int z f(Φ) jest **wiarygodny jedynie w reżimie Φ ≳ 0,8 Φ₀**
> (tj. **|δΦ/Φ₀| ≲ 0,2**)."

Lin. 2839, uzasadnienie reżimu III, **dosłownie**:
> „Dla małych d Φ₁,₂ są **duże** w regionie nakładania, więc człon Φ₁Φ₂(Φ₁+Φ₂) ~ Φ³
> daje wkład silniejszy…"

**To są dwa zdania w tej samej podsekcji, i one się wykluczają.** Reżim III z definicji
wymaga |δΦ/Φ₀| ~ 1 — pięciokrotnie poza zadeklarowanym zakresem.

Dodatkowo: f(Φ) zeruje się przy Φ/Φ₀ = e^(−1/4) = **0,7788** — dokładnie na granicy
deklarowanych 0,8. Poniżej **człon kinetyczny tego funkcjonału jest ujemny**.
Czyli studnia jest liczona funkcjonałem, który w tym obszarze ma **ujemną energię kinetyczną**.
`prop:ghost-free-fundamental` mówi, że to artefakt — zgoda, ale to **nie ratuje E_int**,
tylko potwierdza, że **E_int nie wolno tam używać**.

### (b) ⛔ Twierdzenie o rozbieżności E_β jest podwójnie fałszywe

Lin. 2813–2816, **dosłownie**:
> „Całka ∫Φ₁Φ₂d³x rośnie szybko, gdy d maleje… Dla Φ₁,₂ ~ 1/r całka **dywerguje
> logarytmicznie**, ale jest **regularyzowana przez nieliniową strukturę profilu
> blisko źródeł**."

Numerycznie (d=2, bez ekranowania):

| L | I₂ | I₂/L | I₂/ln L |
|---|---|---|---|
| 10 | 98.98 | 9.90 | 42.99 |
| 20 | 237.10 | 11.86 | 79.15 |
| 40 | 528.86 | 13.22 | 143.36 |
| 80 | 1128.78 | 14.11 | 257.59 |

**I₂/L jest w przybliżeniu stałe; I₂/ln L rośnie 6-krotnie.** Zatem:
1. rozbieżność jest **LINIOWA w L, nie logarytmiczna**;
2. jest w **IR (duże r)**, więc regularyzacja **„blisko źródeł" jej nie usuwa** —
   to naprawa w złym miejscu.

Fizycznie obcina ją dopiero **masa Yukawy**, ale przy λ ~ 10⁶¹ r₀ (z tego samego rdzenia).

### (c) ⛔ Człon E_α jest w dekompozycji, ale nigdy nieobliczony

`eq:Eint-decomp`: E_int = E_lin + E_β + E_γ + **E_α(d)**.
Paragrafy (i), (ii), (iii) liczą trzy człony. **Paragrafu (iv) nie ma.**

A `eq:scales` podaje **r_well ~ (α/β)·r₀** — czyli **skala studni zależy od α**,
które wchodzi wyłącznie przez nieobliczony E_α. `prop:trzy-rezimy-beta-gamma` z kolei
używa **tylko trzech członów, bez α**.

### (d) ⛔ Dwa sprzeczne wzory na te same skale

| | `thm:three-regimes` (eq:scales) | `prop:trzy-rezimy-beta-gamma` |
|---|---|---|
| r_rep | (β/γ)·(qM/Φ₀) | 2β − √(4β²−18βC) |
| r_well | (α/β)·r₀ | 2β + √(4β²−18βC) |
| porządek | **r_well < r_rep** | **d_well > d_rep** ⛔ |
| zależność od masy | r_well **niezależne** od M | d_well **maleje** z M |
| zależność od α, r₀ | **tak** | **nie** |

Numerycznie (β=1): C=0.02 → d_rep=0.092, d_well=3.908; C=0.20 → d_rep=1.368, d_well=2.632.
**Porządek zawsze odwrotny** względem tabeli reżimów. To nie jest literówka w etykietach —
to **dwa niezgodne wyprowadzenia tej samej wielkości**.

## 4. Ilościowo: czy E_γ w ogóle może zdominować?

Z ekranowanym profilem, przy β=γ (warunek próżniowy) i Φ₀=25 (`prop:Lambda-eff`):

|E_γ/E_β| = (3/2)(1/Φ₀)·(I₃/I₂)

| d | I₃/I₂ | \|E_γ/E_β\| |
|---|---|---|
| 8.0 | 0.0928 | 0.0056 |
| 4.0 | 0.1848 | 0.0111 |
| 2.0 | 0.2980 | 0.0179 |
| 1.0 | 0.3847 | 0.0231 |
| 0.5 | 0.4270 | **0.0256** |

**Stosunek ROŚNIE, gdy d maleje — kierunek zgodny z tezą rdzenia** (to jest ta poprawka
mojego błędu z §1). Ale przy najmniejszym d w zakresie E_γ jest wciąż **~40× za słaby**,
by zdominować. Wzrost jest powolny (4,6× na 16× zmiany d).

**Uczciwe zastrzeżenie:** to zależy od profilu. Przy profilu silnie spiętrzonym (Φ ≫ Φ₀)
stosunek może się odwrócić — **ale wtedy jesteśmy poza zakresem ważności z §3(a)**.
> **Albo E_γ dominuje, albo funkcjonał obowiązuje. Nie da się mieć obu naraz.**

## 5. Co z tego jest KONSTRUKTYWNE — i mocniejsze, niż sądziłem

Sprawdziłem, czy stabilne minimum (odległość równowagowa pary — „rozmiar balonu")
zależy od wątpliwego członu kwartycznego. **Nie zależy.**

Po **usunięciu E_γ** zostaje V = −4πC²/d + 8πβC²/d², czyli klasyczna forma −a/d + b/d²:
- V → 0⁻ przy d→∞ (przyciąganie), V → +∞ przy d→0 (odpychanie),
- **minimum przy d\* = 2b/a = 4β**, V(d\*) = −a²/(4b) < 0.

Porównanie: z pełnymi trzema członami (β=1, C=0.1) minimum wypadało przy **d₂ = 3.483**;
bez członu kwartycznego — przy **d\* = 4.000**. **Różnica ~13%.**

> **Stabilna odległość równowagowa pary jest ODPORNA na wywalenie całej zepsutej części.**
> Wynika z konkurencji **przyciągania (E_lin) i odpychania (E_β)** — dwóch członów,
> których znaki są pewne. To jest szczebel, który **stoi**.

I dokładnie to jest potrzebne w obrazie balonów: **naturalna skala upakowania**,
wyprowadzona z (β, C), nie dopasowana. Nie potrzebuje confinementu.

## 6. Czego NIE twierdzę

- **Nie twierdzę, że confinement w TGP nie istnieje.** Twierdzę, że `ssec:dwa-zrodla`
  **go nie ustanawia**, i że §3 i sek09 powołują się na niego jako na wynik.
- Nie sprawdzałem trzech istniejących skryptów rdzenia (`two_body_Veff.py`,
  `three_regimes_quantitative.py`, `two_source_potential.py`) — możliwe, że któryś
  już to łapie albo używa innej regularyzacji. **To osobny check.**
- Reżim II: **znak** (odpychanie, E_β>0) jest pewny z samego wzoru; **skalowanie**
  (1/d)·ln(d/r₀) nie ma pokrycia w rachunku.

## 7. Jednozdaniowo

> Fundament trzech reżimów niesie **solidną grawitację**, **odpychanie o pewnym znaku
> lecz błędnym skalowaniu**, i **studnię, która nie jest ustanowiona** — wyprowadzoną poza
> zadeklarowanym zakresem ważności funkcjonału, dwoma sprzecznymi wzorami, z których jeden
> zależy od nieobliczonego członu; **za to stabilna odległość równowagowa pary przeżywa
> usunięcie całej wątpliwej części i pozostaje wynikiem, na którym można budować.**

---

**Źródła (u źródła):** `core/sek08_formalizm/sek08_formalizm.tex` lin. 2578–2607 (E_int, eq:energy-corrected),
2610–2663 (rem:kinetic-origin + **zakres ważności lin. 2660–2662**), 2665–2691 (prop:ghost-free-fundamental),
2766–2772 (eq:Eint-decomp z E_α), 2776–2803 (E_lin + Yukawa), 2807–2827 (E_β + teza o rozbieżności),
2831–2843 (E_γ), 2847–2872 (thm:three-regimes + eq:scales), 2888–2916 (rem:Veff-Eint-bridge)
**Rachunek:** `REZIMY_audit_dwa-zrodla_v2.py`
