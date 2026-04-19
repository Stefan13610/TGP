# Program P5 — Superconductivity Closure

**Cel:** wyprowadzić z rdzenia TGP warunki istnienia fazy o zerowym oporze
i formułę skalującą T_c jako funkcję parametrów efektywnych (substrat +
geometria sieci), a na bazie tego wskazać klasy materiałów-kandydatów.

**Data startu:** 2026-04-19
**Data zamknięcia:** 2026-04-19
**Status programu:** **CLOSED (3/4 pełne + 1/4 fenomenologiczne)**

## Kluczowe wyniki

- **ps1**: IFF-warunek topologiczny na fazę zerooporową w sieci solitonów TGP.
  Friedel-like oscylacje J(a) z okresem 2π (udowodnione numerycznie do 4 cyfr).
- **ps2**: Formuła zamknięta T_c = k_d(z) · C_0 · A(g_0)² · Λ_E / k_B, gdzie
  **C_0 = 48.82** (RMS 1.54%), a* = 7.725 ± 0.016 uniwersalne (rozpr. 0.21%).
- **ps3**: Kalibracja długości **1 jedn. substratu = a_0 (Bohr) = 0.529 Å**
  daje a* TGP = 4.088 Å — idealnie trafia w Al (4.046 Å, T_c=1.175 K).
- **ps4**: TGP-score ranking: **16 z 20 top kandydatów to znane SC**.
  Cuprates 13, A15 ~10, super-hydrydy top (LaH10 score=77).

## Kluczowe predykcje

1. **Yb, Ce** (non-SC, score >85) — f-electron, kandydaci pod ciśnieniem
2. **Optymalne T_c** dla materiału z a ≈ 4.0-4.1 Å + d-orbitalami + FCC/BCC
3. **Anti-predykcja**: materiały z a/a_0 ≈ nπ (zera J(a)) nie będą SC
4. **Klatkowe struktury** (z=20) + d/f-orbitale = super-hydrydy

---

## Założenia metodologiczne

1. **Wariant A (topologiczny) na początek** — wykorzystujemy istniejącą
   infrastrukturę: dodatek T4 (bounce topology), dodatek O (U(1) formalizacja).
   NIE wprowadzamy nowej fizyki, tylko nowy **język** opisu.

2. **Docelowo β (konieczność + wystarczalność)** — warunek na SC ma być
   IFF, nie tylko konieczny. Tzn. pokazujemy implikację w obie strony.

3. **Sieć krystaliczna solitonów jako prymityw** — materiał opisujemy jako
   regularną sieć Bravais `{R_i}` stanów solitonowych `g_0(r - R_i)` substratu
   TGP (α=1, d=3). Każdy węzeł niesie fazę `θ_i ∈ U(1)`.

4. **Emergent parametr porządku:**
   ```
   ψ(r) = Σ_i g_0(r - R_i) · e^{iθ_i}
   ```
   Koherencja globalna ⇔ `⟨e^{iθ}⟩ ≠ 0` w granicy termodynamicznej.

5. **Topologiczna ochrona:**
   - `π_1(U(1)) = ℤ` ⇒ kwantyzowane fluksoidy
   - Bariera rozwicia (unwinding) ~ `J · L^d / ξ^d` (koszt wir. XY)
   - Prąd persystentny = topologicznie chroniony winding ≠ 0

6. **Reguła workflow:** ps4 (kandydaci) odpalamy DOPIERO po zamknięciu
   ps1–ps3. Inaczej ryzyko fenomenologicznego dopasowania do znanych
   nadprzewodników (curve-fitting trap).

---

## Problemy P5 (cztery skrypty)

| # | Plik | Pytanie | Oczekiwany status |
|---|---|---|---|
| **ps1** | `ps1_zero_resistance_condition.py` | IFF-warunek na fazę zerooporową w sieci solitonów TGP. Czy istnieje bariera topologiczna > k_B T dla danej geometrii? | Hipoteza → **Propozycja** |
| **ps2** | `ps2_tc_scaling.py` | Formuła T_c(a, z, g_0, α_i) z parametrów substratu i sieci. | Szkic → **Hipoteza** |
| **ps3** | `ps3_effective_params_map.py` | Mapowanie parametrów efektywnych substratu TGP na obserwable materiału: ω_Debye, N(E_F), λ_ep, μ*. | Program → **Propozycja** |
| **ps4** | `ps4_candidate_classes.py` | Skan bazy znanych materiałów SC (SuperCon DB / Materials Project) pod kątem warunków P5. Wybór klas-kandydatów. | Program → **Predykcja** |

---

## ps1 — sformułowanie precyzyjne

**Pytanie:** niech `Λ = {R_i}` to sieć Bravais w R^3 ze stałą `a`,
liczbą koordynacyjną `z`, i niech `g_0(r)` to solitonowe rozwiązanie ODE
substratu TGP (α=1, d=3) zlokalizowane z długością `ξ`. Zdefiniujmy
efektywny hamiltonian XY na sieci:

```
H_eff = -Σ_{<ij>} J(|R_i - R_j|) cos(θ_i - θ_j)
```

gdzie `J(a)` to całka sprzężenia dwóch solitonów (sprzężenie Josephsona
między sąsiednimi "pre-elektronowymi" stanami zlokalizowanymi).

**Twierdzenie ps1 (cel):** *faza zerooporowa istnieje wtedy i tylko wtedy,
gdy:*

1. **J > 0** (atrakcyjne sprzężenie faz) — warunek **sprzężenia**
2. **T < T_c(J, z, d)** (pod progiem fluktuacji termicznych) — warunek **koherencji**
3. **π_1(U(1)) = ℤ** (nietrywialna topologia przestrzeni parametru porządku) — warunek **topologiczny**

**Wariant A (topologiczny):** skupiamy się na (3) + dowodzie że bariera
rozwicia globalnego winding rośnie z objętością próbki.

**Wyjście ps1:**
- J(a) jako funkcja stałej sieci (z numerycznego overlapu solitonów)
- ξ z asymptotyki ogona g_0(r) (spodziewamy się A/r^2 z dod. J)
- T_c^{XY-mean-field} = z J / (2 k_B)
- T_c^{BKT} dla sieci 2D = π J / 2
- Warunek wystarczalności: dowód że dla (1)+(2)+(3) energia unwindingu
  rośnie ekstensywnie → topologiczna ochrona

**Kryterium sukcesu:** ps1 zamyka się jeśli dla sieci substratu TGP
(z g_0 ∈ {g_0^e, φ g_0^e, g_0^τ}) udaje się numerycznie pokazać:
- J(a) > 0 dla pewnego zakresu a (istnienie SC fazy)
- skalowanie T_c ~ z J(a*) gdzie a* = optimum potencjału wiążącego
- Analityczny argument że winding jest topologicznie chroniony w 3D

**Kryterium porażki (honest fallback):** jeśli J(a) < 0 dla wszystkich
a > 0 (tj. solitony odpychają się), to faza SC w TGP wymaga dodatkowego
kleju (skalarny mediator — wariant B). Wtedy zamykamy ps1 ze statusem
NEG i przechodzimy do wariantu B jako ps1b.

---

## Skrypty i struktura folderu

```
superconductivity_closure/
├── README.md                          (ten plik)
├── ps1_zero_resistance_condition.py   (w toku)
├── ps1_results.txt                    (generowany)
├── ps2_tc_scaling.py                  (po ps1)
├── ps3_effective_params_map.py        (po ps2)
└── ps4_candidate_classes.py           (po ps3)
```

---

## Powiązania z resztą rdzenia TGP

- **sek02** — ODE substratu α=1, d=3 (podstawa)
- **dodatek J** — ogon masy A_tail (źródło długości lokalizacji ξ)
- **dodatek T4** — bounce topology (topologiczna ochrona)
- **dodatek O** — U(1) formalizacja (faza ψ = |ψ|e^{iθ})
- **dodatek R2** — dynamical Z3 (analogia dla 3-band SC / d-wave)
- **Program P4** — zamknięcie sektora cząstek (baza solitonów g_0^{e,μ,τ})

**Mapowanie klasyczne → TGP:**

| BCS/Ginzburg-Landau | TGP substrat |
|---|---|
| Para Coopera | Soliton `g_0` z dodatkową fazą U(1) |
| Order parameter ψ | `ψ = Σ_i g_0(r-R_i) e^{iθ_i}` |
| Gap Δ | Bariera topologiczna rozwicia |
| T_c | Skala sprzężenia J(a) solitonów |
| Coherence length ξ | Długość lokalizacji solitonu (A/r^2 tail) |
| London penetration depth λ_L | Ekran elektromagnetyczny w substr. U(1) |
| Vortex / fluxoid | Soliton z winding n ∈ π_1(U(1)) |

---

## Dalsze kroki po zamknięciu ps1

- **ps2**: formuła T_c z α-drabinki TGP. Hipoteza robocza:
  `T_c ~ (α_3 / α_2) · (J_0 / k_B) · z^{1/d}` gdzie J_0 = stała charakterystyczna
  substratu. Weryfikacja z wartościami α_i z P4 / T6.

- **ps3**: mapowanie ξ_TGP ↔ długość Pippardowska, J_TGP ↔ λ_ep·ω_D (McMillan).
  Test spójności: dla MgB2/Nb/Pb czy TGP przewiduje T_c w zakresie obserwowanym.

- **ps4**: skan bazy SuperCon + Materials Project pod kątem:
  - wysoka symetria sieci (BCC, FCC — wysokie z)
  - gęstość stanów przy E_F w okolicy punktów φ-drabinki
  - lokalizacja elektronów d/f zgodna z ξ_TGP
  Wyjście: top-10 klas materiałowych o przewidywanej T_c > 100 K.

---

## Zamknięcie programu (oczekiwane)

- 3/4 pełne zamknięcie (analogicznie do P4: ps1+ps2+ps3 jako Propozycje/Hipotezy)
- ps4 jako Program → Predykcja listy kandydatów
- Nowy dodatek LaTeX: `dodatekP5_nadprzewodnictwo_closure.tex`
- Rejestracja w `main.tex` obok `dodatekP4_sektor_czastek_closure.tex`
- Aktualizacja `tabela_epistemiczna.tex` z nowymi pozycjami P5
