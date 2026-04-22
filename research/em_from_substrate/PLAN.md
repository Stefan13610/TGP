# em_from_substrate — program badawczy

**Session:** 2026-04-21
**Cel:** Zderzyć formalizm U(1)/EM z dodatekO + sek09 z niezależną numeryką,
zidentyfikować RZECZYWISTE luki między teorią a weryfikacją, rozszerzyć
o dwuciałowy (Coulomb) test substratowy.

---

## Kontekst — co JUŻ jest w TGP

### Formalizm (zamknięty analitycznie)

1. **dodatekO_u1_formalizacja.tex** — luki O-1, O-2, O-3 zamknięte:
   - O-1: transwersalna dekompozycja (Helmholtz) na substracie
   - O-2: kwantowanie ładunku z winding n∈ℤ zwartej fazy θ
   - O-3: α_em z parametrów substratu (formuła)

2. **sek09_cechowanie.tex** (thm:photon-emergence, 5-krokowy dowód):
   - ψ_i = φ_i·e^(iθ_i) → H = −J·Re(ψ_i*·ψ_j)
   - L_phase = Jv²a²/2 · Σ_μ(∂_μθ)²
   - A_μ ≡ (ℏ₀/e)·∂_μθ
   - F_μν = ∂_μA_ν − ∂_νA_μ z plakietek
   - 1/μ₀ = 2Jv²a_sub²·e²/ℏ₀²
   - α_em = J²·a_sub⁴·v⁴·ε₀/(4π·ℏ·c₀)

3. **Hierarchia defektów** (sek09, ex207): U(1) → SU(2)×U(1) → SU(3) z substratu
   odpowiednio 1-/2-/3-składnikowego (ψ, Ψ∈ℂ², Ξ∈ℂ³).

### Weryfikacja numeryczna (już zrobiona)

| Skrypt | Status | Co testuje |
|---|---|---|
| `ex109_u1_gauge_emergence.py` | 12 testów T1–T12 | gradient fazy → A_μ, plakietki → F_μν, niezmienniczość cech., bezmasowość, kwantyzacja wiru, granica ciągła, relacja 1/μ₀ |
| `alpha_em_rg_flow.py` | 9/9 PASS | odwrócony przepływ RG: α(m_e)=1/137.036 → α(ℓ_P)≈1/94.09 |
| `alpha_em_substrate_v2.py` | hipoteza robocza | 3 ścieżki do α_em (η_WF, RG, topologiczna θ_Z₂) — rząd wielkości |
| `tgp_agamma_phi0_test.py` | n-1 PASS | a_Γ·Φ₀=1 z 3 ścieżek (S1, S2c, S3) + DESI DR2 |
| `ex207_homotopy_defect_chain.py` | 20/20 PASS | topologia π_k defektów U(1)/SU(2)/SU(3) |

### Co jest otwarte

Po analizie istniejącego kodu + dokumentacji, RZECZYWISTE luki to:

**G1. Bezpośrednie wyprowadzenie α_em z (J, v, a_sub)** — wersja z dodatekO
eq:alpha-em-substrate pokazuje formułę, ale `alpha_em_substrate_v2.py` nie
testuje JEJ — testuje alternatywne hipotezy (η_WF/(4π N_c), θ_Z²). Potrzebny
skrypt który BIERZE formułę z twierdzenia i liczy α_em dla zadanych (J, v, a_sub)
w jednostkach ℓ_P, pokazując czy daje 1/137.036 przy fizycznie sensownych
parametrach substratu.

**G2. Coulomb z dwóch wirów fazowych** — sek09 wyprowadza F_μν z plakietek, ale
NIE pokazano numerycznie że dwa wiry fazowe n=±1 oddziałują jak
V(r) = ±e²/(4πε₀r). To jest właściwy test „EM emerguje z substratu".
ex109 testuje pojedyncze pole — nie interakcję między wirami.

**G3. Relacja J_amp ↔ J_phase** — sek09 traktuje substrat jako JEDNO pole ψ.
Ale sek04/sek08 (grawitacja) dają J_amp (sprzężenie amplitudy φ) różne od
J_phase. Czy ma to sens? Czy J_amp i J_phase mogą być TĄ SAMĄ stałą
z RG flow? ex207 pokazuje hierarchię ale nie kontrastuje sprzężeń.

**G4. Ładunek elektronu z winding number** — sek09 pokazuje że ∮∂θ dl = 2πn
(kwantyzacja ℤ), ale nie wyprowadza DLACZEGO elektron ma n=+1 (a nie np. n=3).
Czy n dla leptonu jest związane z g₀^e ≈ 0.86941? To przynosi mechanizm
TGP dla generacji ładunku.

**G5. Propagator fotonu z MC substratu** — granica ciągła ex109 dostaje
ω² = k² na sieci, ale dla małego L. Full propagator D_μν(k) wymaga MC
większej objętości i/lub spektralnego analizera fluktuacji ψ_i. To wzmacnia
T10 (brak masy) do pełnej funkcji korelacyjnej.

**G6. Substrat-EM ↔ atom** — czy pole A_μ wygenerowane z substratu sprzęga
się do elektronu (np. bound state w Coulombie) tak samo jak fizyczny foton?
Test: równanie Kleina-Gordona w zewnętrznym A_μ zbudowanym z wiru → 13.6 eV
dla hydrogenic.

---

## Plan skryptów

| # | Skrypt | Cel | Kluczowy test |
|---|---|---|---|
| em00 | `em00_baseline_diagnostic.py` | Zebrać wyniki istniejących skryptów, identyfikować faktyczne luki | Matrix „co jest done / open" |
| em01 | `em01_alpha_em_direct_substrate.py` | **G1**: α_em = J²a⁴v⁴ε₀/(4π ℏ c₀) dla fizycznych (J, v, a_sub) | Czy α_em(J=0.366, v=1, a=ℓ_P) → 1/137? |
| em02 | `em02_two_vortex_coulomb.py` | **G2**: Dwa wiry fazowe n=±1 na 3D siatce → oddziaływanie V(r) | Fit V(r) ∝ 1/r, stała=e²/(4πε₀) |
| em03 | `em03_J_amp_vs_J_phase.py` | **G3**: RG flow obu sprzężeń z różnych skal | Czy J_amp(ℓ_P)=J_phase(ℓ_P)? |
| em04 | `em04_winding_to_charge.py` | **G4**: winding n ↔ ładunek elektronu | Czy n=1 jest najniższym stanem stabilnym? |
| em05 | `em05_photon_propagator_mc.py` | **G5**: propagator D_μν(k) z MC substratu | Pole k=0, residuum → 1/μ₀ |
| em06 | `em06_atom_in_vortex.py` | **G6**: Kleina-Gordona w A_μ(wir) → E_1s | Czy −13.6 eV emerguje? |
| VERDICT | `EM_VERDICT.md` | Podsumowanie — co zamknięte, co zostaje otwarte | — |

---

## Kryteria falsyfikacji (twarde)

- **S1**: em01 PASS jeśli |α_em(J,v,a)/α_obs − 1| < 10% dla (J, v, a_sub) z fizycznych
  rozważań (J~0.366 z prop:alpha-em, v~parametr porządku z sek09, a_sub~ℓ_P).
  Jeśli błąd > 10× → formuła eq:alpha-em-substrate jest FAŁSZYWA lub brakuje
  renormalizacji — trzeba otwierać nowy front.

- **S2**: em02 PASS jeśli V(r) wiru_+/wiru_− ma slope 1/r z dokładnością < 5%
  na długościach 5–20 a_sub. Jeśli slope ≠ 1 → substrat NIE daje Coulomba
  w 3+1D (poważny kryzys — być może potrzebny 4D czasoprzestrzenny MC).

- **S3**: em06 PASS jeśli E_1s ∈ [−14, −13] eV dla elektronu w A_μ wygenerowanym
  przez wir kwantowy n=1. Jeśli wartość nie zbiega do −13.6 eV → sprzężenie
  elektron–substrat jest inne niż e²/(4πε₀r).

---

## Zasady epistemiczne

1. **Nie powtarzamy weryfikacji** — jeśli ex109 już przetestował element,
   em0X nie duplikuje tego; cytujemy i idziemy dalej.
2. **Uczciwie dokumentujemy** — jeśli em01 fails dla formuły z dodatekO,
   piszemy to w VERDICT bez upiększeń. To jest właśnie ten moment kiedy
   atomowa porażka (as02/03) musi być zapłacona — EM musi ratować albo
   TGP ma drugą poważną lukę.
3. **Kolejność**: em00 → em01 → em02 są KRYTYCZNE (formuła + Coulomb).
   em03–em06 to "nice to have" jeśli czas pozwala.

---

## Po em_from_substrate

Zgodnie z planem użytkownika: NASTĘPNY program to **cohesion_closure**
(energia kohezji metali, wiązania chemiczne) — patrz CLOSURE_PLAN (do utworzenia).
