# Główne założenia warstwy `nbody`

Skrót dla czytelnika kodu. Szczegóły i mapowanie na pliki: [ANALIZA_NBODY_INTEGRACJA.md](ANALIZA_NBODY_INTEGRACJA.md), glosariusz w [`__init__.py`](__init__.py).

## 1. Rdzeń pola TGP (N0)

- Przestrzeń opisana jest polem skalarnym Φ; próżnia Φ = Φ₀ jest stabilna przy **β = γ** (warunek próżni z potencjału efektywnego).
- Masa ekranowania liniaryzacji: **m_sp² = 3γ − 2β**; przy β = γ: **m_sp = √γ**.
- Źródła punktowe: siła „ładunku” **C** w profilu typu Yukawy jest powiązana z interpretacją bezwymiarową używaną w kodzie (patrz `tgp_field.py`, `pairwise.py`).

## 2. Droga B (źródła vs linearyzacja)

Linearyzacja przy Φ = Φ₀ daje oscylacyjny ogon **∝ sin(r)/r**, a nie czystą Yukawę. W pakiecie **profil Yukawy dla ciał punktowych jest traktowany jako warunek źródłowy / Droga B** (defekt z nałożonym zanikiem), a nie jako swobodne rozwiązanie pełnego równania w całej przestrzeni. Zob. §1.2 w analizie integracyjnej.

## 3. Konwencje dynamiki wielociałowej

- **m_i = C_i** — w module przyjęta jest ekwiwalencja masy inercyjnej i siły źródła (jak w `__init__.py`).
- **softening** w `dynamics_v2.py` to regulator numeryczny (√(r²+ε²)), nie twierdzenie fizyczne; wnioski wymagają sprawdzenia zbieżności w ε.
- Siły 2-ciałowe: zamknięta postać w sektorze parametryzowanym przez C, β, γ, m_sp.
- Siły 3-ciałowe z **Φ⁴**: człon niezerowy; w granicy Coulomba znane są formuły zamknięte; dla ogólnej geometrii Yukawy całka potrójna **nie ma prostej formy zamkniętej** — stąd rozróżnienie na `three_body_terms` (przybliżenia) i `three_body_force_exact` (całka numeryczna 2D Feynmana). Zob. [ANALIZA_3CIALA.md](ANALIZA_3CIALA.md).

## 4. Równowagi i Earnshaw

W TGP pojawia się **bariera repulsywna** ∝ 1/d², która zmienia jakość problemu w porównaniu z czystym Newtonem/Coulombem: **statyczne równowagi 3+ ciał w wolnej przestrzeni są dopuszczalne** w tym modelu (klasyczne wyniki w notatkach TeX i w analizie 3-ciał).

## 5. Skrypty `examples/ex*.py` — osobna warstwa założeń

Każdy plik `exNNN_*.py` może dodawać **własne** założenia fenomenologiczne (np. konkretne ODE na g(r), ghost boundary, dopasowanie ogona **A_tail**, mapa **m ∝ A_tail⁴**, warunki Koidego). Te założenia **nie wynikają automatycznie** z punktów 1–4; są hipotezami testowanymi numerycznie. Status zamknięć (PASS/FAIL, interpretacja) jest w docstringach skryptów i w dokumencie spójności w korzeniu repozytorium: **`ANALIZA_KRYTYCZNA_v6.md`** (starsze wersje w `_archiwum/stare_analizy/`).

## 6. Klasyfikacja ścisłości (przyjęta w pakiecie)

Jak w `__init__.py`: **EXACT** (analityka), **APPROXIMATE**, **NUMERICAL**, **HEURISTIC** — przy czytaniu wyników warto sprawdzać, która etykieta dotyczy danego modułu lub sekcji skryptu.
