# Program: Atomic shells — TGP na poziomie atomu

**Data otwarcia:** 2026-04-21
**Motywacja:** Po częściowym zamknięciu d-class gap w SC (ps41 → 60% RMS redukcji, Lu pozostał 50× off), ujawnił się **strukturalny problem**: empiryczne amplitudy orbitali
`A_s = −0.111`, `A_sp = +0.207`, `A_d = +0.310`, `A_f = +2.034` używane w SC i ρ(T)
są **dopasowane**, a nie wyprowadzone. Równocześnie rdzeń TGP wyprowadza masy leptonów
z potęgi A_tail⁴ z precyzją 0.013%. Obie te struktury operują na "amplitudzie ogona",
ale nigdy nie zostały połączone.

**Pytanie badawcze:** Czy TGP opisuje **poprawnie atom**? Jeśli tak — z czego dokładnie
wynikają wartości A_orb? Jeśli nie — na czym się wykłada?

---

## 0. Dlaczego lit?

Lit jest **pierwszym nietrywialnym atomem**:

| Pierwiastek | Struktura | Status testu TGP |
|---|---|---|
| H (Z=1) | 1 e⁻ w 1s | test SANITY: czy TGP reprodukuje Rydberga? |
| He (Z=2) | 2 e⁻ w 1s², closed shell | test 2-elektronowej korelacji (singlet) |
| **Li (Z=3)** | **[He] 2s¹** | **pierwszy atom z core + valence — NATURALNY TEST A_s** |
| Na (Z=11) | [Ne] 3s¹ | poziom kontrolny (cięższy alkalik) |
| Cs (Z=55) | [Xe] 6s¹ | upper bound — relatywistyka zaczyna dominować |

Li ma **3 protony + 3 elektrony + 3 lub 4 neutrony (Li-6, Li-7)**. 2s elektron walencyjny jest w polu jednego protonu efektywnego (Z_eff ≈ 1.26) — idealny test "A_s"
w najprostszym możliwym kontekście atomowym.

**Cały ciąg jonizacyjny Li:**
- IE₁ = 5.3917 eV (odrywanie 2s)
- IE₂ = 75.640 eV (odrywanie z 1s² → 1s¹)
- IE₃ = 122.454 eV (hydrogenic Li²⁺)

IE₃ jest **dokładnie hydrogenic**: 13.6057·3² = 122.45 eV. Ten jeden numer daje nam
rygorystyczny test — jeśli TGP tego nie dostaje, cała reszta nie ma znaczenia.

---

## 1. Status: co TGP MA, co TGP musi dostarczyć

### 1.1 Co TGP już ma (rdzeń)

| Element | Ustalone | Referencja |
|---|---|---|
| Soliton ODE `g'' + (g')²/g + 2g'/r = 1 − g` | TAK | `mass_scaling_k4`, `why_n3` |
| Tail: `g(r) ~ 1 − A·e⁻ʳ/r` | TAK | R5, ODE numerics |
| Masy leptonów: `m ∝ A⁴` | TAK (0.013%) | `r5_k_squared_mechanism.py` |
| `g₀^e = 0.86941`, `g₀^μ = φ·g₀^e`, `g₀^τ = 1.73027` | TAK | ps2 (P4) |
| Koide K=2/3 dla leptonów naładowanych | WYPROWADZONE | ps2 (P4) |
| `ℏ = π·χ·A_tail(Φ)` | TAK (3 drogi) | `qm_foundations/q0_analytical.py` |
| `χ ≈ 0.86` uniwersalne (CV 7%) | TAK | q0_analytical |
| Metryka substratowa `ds² = −c₀²g·dt² + (1/g)·δᵢⱼdxⁱdxʲ` | TAK | sek08c |

### 1.2 Czego TGP NIE ma (gap diagnozowany)

| Element | Stan | Krytyczność |
|---|---|---|
| **Stan związany H w TGP** | brak skryptu | 🔥 sanity blocker |
| **Atom wieloelektronowy jako rozwiązanie ODE** | brak | 🔥 główny test |
| **Derywacja A_orb(s/sp/d/f) z atomu** | brak (empiryczne z SC) | 🔥 główny cel |
| **Energia Rydberga z ℏ = πχA_tail** | brak weryfikacji | 🔥 sanity |
| **Defekt kwantowy δ_s^Li z TGP** | brak | średnio |
| **EM w TGP** | brak dokumentacji | fundamentalnie otwarte |

**Główna zagadka:** SC/ρ(T) używają A_s = **−0.111** (ujemna!). W rdzeniu TGP A_tail
dla **swobodnego elektronu** to A_e ≈ +0.1246 (dodatnia, z A_tail⁴ = m_e). Dlaczego
A_orb dla **s-walencji w metalu** ma przeciwny znak? To jest sygnał że:
- albo A_orb w SC NIE jest tym samym co A_tail rdzenia,
- albo podłoże atomowe zmienia fazę,
- albo coś w SC/ρ(T) jest nadrozwinięte empirycznie.

Probe P3 odpowie, który scenariusz jest poprawny.

---

## 2. Cele (A1–A4)

### A1. Sanity: wodorowy 1s
Dla jednego elektronu w potencjale coulombowskim pojedynczego protonu-solitonu:
1. Policzyć E_1s z efektywnego równania (Schrödinger z `ℏ(Φ)` lokalnym).
2. Porównać z E_1s^exp = −13.6057 eV.
3. Weryfikacja skali `a_0^TGP = ℏ²/(m_e·e²)` z wyprowadzonymi ℏ, m_e.

**Falsyfikator A1:** jeśli TGP przewiduje |E_1s − 13.6 eV| > 1%, atomowe zastosowania
są od razu wykluczone.

### A2. Test Rydbergowej skali Z² (H, He⁺, Li²⁺)
Dla hydrogenicznych jonów (1 elektron, ładunek Z=1,2,3):
- E_1s ∝ Z²: H: −13.6, He⁺: −54.4, Li²⁺: −122.4
- Czy TGP odtwarza dokładnie skalowanie Z²?
- Czy relatywistyczna korekta `α² Z⁴` pojawia się naturalnie z lokalnego `ℏ(Φ)`?

**Falsyfikator A2:** jeśli skala Z² jest załamana powyżej 0.1%, TGP nie opisuje poprawnie
hydrogenicznego atomu.

### A3. Li 2s¹ walencja — pierwszy test wieloelektronowy
Zadanie: policzyć IE₁(Li) = 5.39 eV z TGP.
- Opcja (a): **Traktuj Li jako H-like ze screeningiem**: `E_2s = −13.6·Z_eff²/n²`. Dla n=2, IE₁=5.39 → Z_eff = 1.26. Czy TGP wyprowadza Z_eff = 1.26 z przesłaniania 1s² przez substrat?
- Opcja (b): **Pełne 3-solitonowe rozwiązanie** (jądro + 2 core e + 1 valence e) w TGP — czy ten jedyny walencyjny soliton ma A_tail powiązane z A_s = −0.111?

**Falsyfikator A3:** jeśli TGP nie daje ani Z_eff ≈ 1.26, ani |A_valence| ≈ 0.11, to
A_s z SC jest czystym fitowaniem, a Li — strukturalny brak.

### A4. Most do SC: A_orb(s/sp/d/f) z atomów
Dla alkali (Li, Na, K, Rb, Cs) i pełnych rzędów (s, sp, d, f):
- Wylicz A_valence(atom) z TGP z A3.
- Porównaj z empirycznymi A_s = −0.111, A_sp = +0.207, A_d = +0.310, A_f = +2.034.
- Test korelacji: czy TGP-derywowane A(atom) pokrywa się co do znaku i względnego uporządkowania?

**Falsyfikator A4:** brak korelacji między A(atom) a A_orb(SC) → A_orb są phenomenologiczne,
nie strukturalne. Kończymy SC w paperze v1 z honest note "A_orb fitted".

---

## 3. Plan numeryczny (as01–as06)

| # | Skrypt | Cel | Przewidywany wynik |
|---|--------|-----|--------------------|
| as01 | `as01_hydrogen_probe.py` | H 1s w TGP (A1) | E = −13.6 eV? |
| as02 | `as02_rydberg_Z2_scaling.py` | Hydrogenic H/He⁺/Li²⁺ (A2) | Z² exactness? |
| as03 | `as03_lithium_screening.py` | Li IE₁ via Z_eff (A3a) | Z_eff ≈ 1.26? |
| as04 | `as04_lithium_full_TGP.py` | Li jako 3-soliton konfiguracja (A3b) | A_valence ~ 0.11? |
| as05 | `as05_Aorb_derivation.py` | A_s/sp/d/f z atomów (A4) | korelacja z SC? |
| as06 | `as06_alkali_series.py` | Li→Cs, kontrola alkali-relatywizmu | Cs odchylenia? |

**Dependency graph:**
```
as01 → as02 (Rydberg extension)
as01 → as03 (Li screening uses H baseline)
as03 + as04 → as05 (A_orb synthesis)
as05 → as06 (alkali series validates A_s)
```

---

## 4. Hipotezy robocze (do sfalsyfikowania)

### H1. "TGP = standardowa QM w granicy wolnej"
Dla H w TGP `ℏ(Φ) → ℏ₀` w granicy `Φ → Φ₀`, więc 1s TGP = 1s Dirac. Przewidywanie:
E₁s ≈ −13.6 eV z relatywistycznymi poprawkami α²Z⁴·const.

### H2. "A_orb są efektywne, nie strukturalne"
A_s = −0.111 z SC jest dopasowane do **chemicznej hybrydyzacji w siatce** (nie
pojedynczego atomu). TGP wyprowadzi A_valence dla **izolowanego atomu Li** i wyjdzie
**inny znak i inna magnituda** niż A_s(SC). Wtedy A_orb są fit, nie pierwsze zasady.

### H3. "TGP ma derywowany defekt kwantowy"
TGP z lokalnym ℏ(Φ) przewidzi δ_s^Li = (n − n*)·n ≈ 0.4 z substratowej polaryzacji
core 1s². Jeśli tak — stara struktura Madelunga jest konsekwencją substratu.

### H4. "Strukturalny brak — atomy potrzebują nowego mechanizmu"
Jeśli as01–as04 zawiodą (|E_H − 13.6| > 1%), to TGP **nie opisuje atomów**. Program
P4 (leptony) był wynikiem SZCZEGÓLNYM — masy wolnych solitonów. Atomy wymagają
oddzielnego wyprowadzenia (np. multi-soliton bound states z nowym ansatz). To jest
**pierwszy uczciwie udokumentowany fundamentalny gap**.

---

## 5. Literatura startowa (Li)

| Pozycja | Wartość | Źródło |
|---|---|---|
| IE₁ (Li) | 5.39171 eV | NIST ASD |
| IE₂ (Li) | 75.6400 eV | NIST ASD |
| IE₃ (Li) | 122.4544 eV | NIST ASD |
| IE₁ (H) | 13.59844 eV | CODATA |
| m_e | 0.5109989 MeV | PDG |
| ℏ | 1.054572·10⁻³⁴ J·s | CODATA |
| a_0 | 0.52917721 Å | CODATA |
| Li polaryzowalność atomu | 164.19 a.u. | CRC |
| Li core polaryzowalność (Li⁺) | 0.1923 a.u. | Tang (2003) |
| δ_s^Li (quantum defect) | 0.4017 | Drake (1998) |
| E(Li ground, atomowa) | −7.47806 Ha = −203.487 eV | HF benchmark |
| E(Li⁺ 1s²) | −7.27991 Ha = −198.098 eV | HF benchmark |

---

## 6. Relacje z innymi sektorami TGP

- **mass_scaling_k4**: dostarcza `m ∝ A⁴` i mechanizm tail-core — będzie użyte w as03/as04.
- **qm_foundations**: `ℏ = πχA_tail` — użyjemy w as01 jako podstawowe equation.
- **particle_sector_closure (P4)**: g₀^e = 0.86941 — referencyjne g₀ dla elektronu.
- **superconductivity_closure**: A_s, A_sp, A_d, A_f empiryczne — do porównania w as05.
- **rho_normal_state_closure**: te same A_orb klasy — drugi test walidacyjny.

---

## 7. Falsyfikowalność programu

Na końcu as01–as06 będą cztery możliwe stany:

| Stan | Kryterium | Wniosek |
|---|---|---|
| **S1. CLOSURE** | as01 ✓, as05 korelacja r > 0.9 | TGP atomowe jest strukturalnie zamknięte. A_orb są derywowane. |
| **S2. PARTIAL** | as01 ✓, as05 korelacja 0.5<r<0.9 | TGP atomowe OK dla prostych atomów; złożone A_orb wymagają poprawek. |
| **S3. STRUCTURAL LIMIT** | as01 ✓, as05 r < 0.5 | TGP OK dla H i H-like; A_orb z SC są PHENOMENOLOGIAMI bez pierwszego zasady. |
| **S4. FUNDAMENTAL GAP** | as01 ✗ (|E_H − 13.6| > 1%) | TGP nie opisuje atomów jak jest. Potrzebny nowy mechanizm. |

**Każdy z tych stanów jest POSTĘPEM.** Program ma strukturę symetrycznego ostrzału —
albo się zamyka, albo daje **pierwszy uczciwie udokumentowany fundamentalny limit**.

---

## 8. Otwarte pytania (do zastanowienia podczas as01–as06)

1. Czy TGP zawiera elektromagnetyzm, czy tylko grawitację/substrat? Jeśli nie, to jest
   zasadnicze pytanie — atom Coulombowski wymaga siły 10³⁹ razy silniejszej niż TGP oferuje.
2. Czy ψ elektronu w atomie = `g(r)−1` dla tego samego solitonu, czy jest NOWYM polem?
   W rdzeniu TGP masa elektronu jest statyczną własnością solitonu, nie falową funkcją.
3. Czy g₀^e = 0.86941 obowiązuje TAKŻE gdy elektron jest związany w atomie? Czy g₀
   dryfuje z `Φ_local(atom)`? Jeśli tak — kolejne pytanie: ile?
4. Matematyczny status wielu solitonów w tym samym substracie: czy superpozycja jest
   liniowa (jak w qm_superposition), czy wymaga nieliniowych poprawek na odległości atomowe?

---

## 9. Link do rdzenia TGP

- [[TGP/TGP_v1/papers/core/tgp_core.tex]] — Paper 1 (fundament)
- [[TGP/TGP_v1/research/qm_foundations/README.md]] — emergentna QM
- [[TGP/TGP_v1/research/mass_scaling_k4/README.md]] — m∝A⁴
- [[TGP/TGP_v1/research/particle_sector_closure/README.md]] — sektor leptonowy
- [[TGP/TGP_v1/research/superconductivity_closure/D_CLASS_SC_CLOSURE.md]] — co domknął ps41, co pozostało
- [[TGP/TGP_v1/research/superconductivity_closure/RHO_SC_CROSSCHECK.md]] — A_orb jako wspólny invariant

---

*Status: OPEN. Next action: `as01_hydrogen_probe.py` — test Rydberga jako sanity blockera.*
