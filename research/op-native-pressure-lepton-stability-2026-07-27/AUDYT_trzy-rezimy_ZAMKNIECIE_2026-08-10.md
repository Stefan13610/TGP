# Zamknięcie audytu trzech reżimów — pięć skryptów, trzy niezgodne rachunki, jeden ocalały wynik

**Data:** 2026-08-10
**Typ:** ZAMKNIĘCIE AUDYTU (wszystkie skrypty rdzenia sprawdzone i uruchomione)
**Skrypty rdzenia:** `two_body_Veff.py`, `three_regimes_quantitative.py`,
`two_source_potential.py`, `effective_potential.py`
**Weryfikacja:** `REZIMY_audit_trzy-rezimy.py`, `REZIMY_audit_dwa-zrodla_v2.py`,
`REZIMY_audit_skrypty_rdzenia.py`, `REZIMY_trzy-niezgodne-rachunki.py`

---

## 0. Werdykt końcowy

> **Studnia (reżim III / confinement) nie jest ustanowiona przez żaden z pięciu rachunków
> w rdzeniu.** Trzy z nich liczą tę samą wielkość i dają **wzajemnie sprzeczne** wyniki —
> różny porządek skal, różną zależność od masy, **różny mechanizm**. Żaden nie jest
> oznaczony jako obalający pozostałe.
> **Ocalały wynik: d\* = 4β** — stabilna odległość równowagowa pary, odporna na usunięcie
> całej wątpliwej części, potwierdzona niezależnie przez skrypt rdzenia co do cyfry.

## 1. Nowe znalezisko: `two_source_potential.py` dostaje trzy reżimy BEZ członu kwartycznego

Skrypt raportuje **13/13 PASS**. Ale jego własny test T1a brzmi:
```
[PASS] T1a: Trzy rezimy bez E_gamma (E_lin+E_beta): 2 zmiany znaku
            -- zmiany znaku przy d ~ [0.322, 6.000]
```

**Trzy reżimy z dwóch członów — bez tego, który wg rdzenia „tworzy studnię".**
Sprawdziłem mechanizm (profil Yukawy, m=1):

| d | E_lin ~ −e^(−d)/d | I₂(d) | I₂/\|E_lin\| |
|---|---|---|---|
| 0.10 | −9.048 | 5.685 | **0.628** |
| 0.30 | −2.469 | 4.654 | 1.885 |
| 1.00 | −0.368 | 2.311 | 6.282 |
| 4.00 | −0.005 | 0.115 | 25.11 |
| 8.00 | −0.00004 | 0.002 | **50.19** |

Mechanizm jest banalny:
- **małe d**: E_lin ~ **−1/d dywerguje**, a I₂ jest **skończone** ⟹ przyciąganie
- **średnie d**: I₂ jeszcze duże, E_lin już stłumione ⟹ odpychanie
- **duże d**: I₂ ~ e^(−2md) znika **dwa razy szybciej** niż E_lin ~ e^(−md) ⟹ przyciąganie

> ⛔ **„Studnia" w tym skrypcie to po prostu osobliwość −1/d punktowego źródła.**
> To **ten sam człon**, który na dużych d nazywamy grawitacją. Nie ma w niej ani członu
> kwartycznego, ani confinementu — jest w niej brak regularyzacji źródła.

## 2. Test, który przechodzi, choć pokazuje załamanie rachunku

Ten sam skrypt, T4b — **rzeczywiste liczby**:
```
[PASS] T4b: E_gamma poglabia studnie przy d=1
            V_bez=0.9057,  V_z=-16.8694
```
To jest **zmiana o 19,6× wartości wyjściowej, ze zmianą znaku**.

**E_gamma nie jest poprawką — całkowicie dominuje.** Rozwinięcie perturbacyjne, na którym
stoi cały rachunek, **tam nie zbiega**. Skrypt zalicza to jako sukces („pogłębia studnię").

Dodatkowo: `Obliczam I3 dla kilku wartosci d (d >= 0.8)` — **małe d są pominięte**,
czyli dokładnie ten obszar, który miałby stanowić studnię.

## 3. `effective_potential.py` — brak wartości weryfikacyjnej

Cały output:
```
Vacuum case: beta=gamma=5.0, C=0.3 ... (x4)
Saved .../effective_potential.png
Saved .../regime_map.png
Done.
```
**Żadnej liczby, żadnego testu, żadnego kryterium.** Rysuje wykresy. Nie potwierdza niczego.

## 4. Trzy niezgodne rachunki tej samej wielkości

| źródło | r_well | r_rep | porządek | uwaga |
|---|---|---|---|---|
| `thm:three-regimes` (eq:scales) | (α/β)·r₀ | (β/γ)(qM/Φ₀) | **r_well < r_rep** | zależy od **α — nigdy nieobliczonego** |
| `prop:trzy-rezimy-beta-gamma` | 2β+√(4β²−18βC) | 2β−√(...) | **d_well > d_rep** ⛔ | brak α, brak r₀ |
| `two_source_potential.py` | ~0.32/m | ~6.0/m | **r_well < r_rep** | trzy reżimy **bez E_γ** |

Cztery niezgodności:
1. **dwa różne porządki skal** (r_well<r_rep vs d_well>d_rep),
2. **dwie różne zależności od masy** (niezależne od M vs malejące z M),
3. **dwa różne mechanizmy studni** (człon kwartyczny vs osobliwość −1/d),
4. **jeden zależy od członu E_α, który nigdy nie został policzony.**

**Żaden nie jest w rdzeniu oznaczony jako obalający pozostałe.** Współistnieją, każdy
z własnym „PASS".

## 5. Pełny bilans po audycie

| element | status |
|---|---|
| **Reżim I — grawitacja** (E_lin ~ −1/d) | ✅ **stoi**; Yukawa pomijalna do 50 rzędów (sprawdzone) |
| **Reżim II — odpychanie** (E_β > 0) | ⚠ **znak pewny**; skalowanie (1/d)ln(d/r₀) **fałszywe** (rozbieżność liniowa, IR, nie log/UV); zasięg absurdalny (10²⁶ m dla protonu) |
| **Reżim III — studnia / confinement** | ⛔ **NIE USTANOWIONA** — poza zakresem ważności funkcjonału (rdzeń sam deklaruje \|δΦ/Φ₀\|≲0,2), punkt **niestabilny** (V''<0 — bariera, nie studnia), 17 rzędów pod Planckiem, trzy sprzeczne wzory |
| **„makro ⟹ tylko grawitacja"** | ⛔ **fałszywe o 20+ rzędów** (M_crit = 8×10¹⁹ M_☉ — gromady galaktyk są poniżej) |
| **d\* = 4β** (odległość równowagowa) | ✅ **ODPORNE** — przeżywa usunięcie E_γ (4.000 vs 3.483, 13%), potwierdzone przez skrypt rdzenia co do cyfry |
| E_α w dekompozycji | ⛔ **nigdy nieobliczony**, a skale od niego zależą |
| eq:Eint = E[Φ₁+Φ₂]−E[Φ₁]−E[Φ₂] | ⛔ **niespójne** (Φ₁+Φ₂ → 2Φ₀); skrypty cicho używają δΦ i mają rację |

## 6. Wzorzec, który się powtarza — po obu stronach

W tej sesji złapałem u siebie: zaśmiecony pudłem test całek, odwrócony kierunek I₃/I₂,
przypisanie nośnika na podstawie nazwy obiektu. **W rdzeniu znalazłem ten sam wzorzec:**

- `two_body_Veff.py`: SUMMARY „reproduces the three-regime structure" przy **jednym**
  przejściu przez zero we własnym outpucie,
- `three_regimes_quantitative.py`: punkt niestabilny (V''=−5,23) oznaczony `(stable)`,
- `two_source_potential.py`: PASS na teście, którego liczby pokazują **rozpad rozwinięcia** (19,6×),
- `effective_potential.py`: brak jakiegokolwiek kryterium.

To nie jest zarzut o nieuczciwość — to **brak testu, który mógłby nie przejść**.
Cztery skrypty, ~17 „PASS", i **ani jeden FAIL** mimo trzech sprzecznych wyników.

## 7. Co realnie zostaje do budowania

**Jeden szczebel: d\* = 4β.** Wynika z konkurencji dwóch członów o pewnych znakach
(przyciąganie E_lin + odpychanie E_β) i **nie potrzebuje** ani studni, ani confinementu,
ani członu kwartycznego. To jest naturalna **odległość równowagowa pary** — czyli
w obrazie balonów **skala upakowania**.

**Jego jedyny problem jest ten sam, co wszędzie:** przy kalibracji Φ₀≈25 z ciemnej energii
wychodzi r₀ = 2,65×10²⁵ m, więc d\* = 1,06×10²⁶ m — **kosmologiczne, nie cząstkowe**.

> **Blokada całego programu jest teraz jedna i konkretna: co ustala Φ₀ poza domeną
> kosmologiczną.** Rdzeń deklaruje Φ₀ jako „EFT scale-dependent free parameter".
> Dopóki jest wolny, **każda skala wyjdzie kosmologiczna**, a d\* nie jest predykcją.

## 8. Czego NIE twierdzę

- Nie twierdzę, że confinement w TGP jest niemożliwy — twierdzę, że **żaden z pięciu
  istniejących rachunków go nie pokazuje**, a §3 i `sek09` powołują się na niego jako na wynik.
- Nie twierdzę, że skrypty są bezwartościowe — `two_source_potential.py` ma porządną
  weryfikację I₂ (analityczna vs numeryczna, błąd ~1e−2) i **poprawnie** identyfikuje,
  że trzy reżimy nie potrzebują E_γ. To użyteczny wynik, tylko źle opisany.
- Nie sprawdziłem, czy przy **innej** kalibracji Φ₀ skale trafiają w fizykę cząstek.
  To jest następne pytanie, nie wniosek.

---

**Uruchomione:** `tooling/scripts/gravity/two_body_Veff.py`, `tooling/scripts/gravity/effective_potential.py`,
`tooling/scripts/stability/three_regimes_quantitative.py`, `tooling/scripts/profiles/two_source_potential.py`
**Powiązane:** `ANALIZA_trzy-rezimy_intuicja-gradientowa_2026-08-10.md`,
`AUDYT_ssec-dwa-zrodla_2026-08-10.md`, `AUDYT_skrypty-trzy-rezimy_2026-08-10.md`
