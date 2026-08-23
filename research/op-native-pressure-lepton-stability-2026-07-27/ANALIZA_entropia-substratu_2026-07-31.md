# Analiza: czym jest „entropia substratu" i dlaczego rdzeń się na nią powołuje

**Data:** 2026-07-31
**Typ:** ANALIZA ŹRÓDŁOWA (tylko analiza — bez liczenia, bez budowania)
**Metoda:** czytanie rdzenia u źródła (`axioms/substrat/dodatekB_substrat.tex`,
`axioms/notacja/dodatekA_notacja.tex`, `core/sek07_predykcje`, `core/sek08c`)
**Pytanie użytkownika:** substrat ma dostarczać tylko minimalne warunki dla złożonych obiektów;
entropia wygląda na próbę stabilizacji zaburzeń — czy to „błąd młodości programu"?

---

## 0. Odpowiedź krótka

**Tak — entropia substratu w TGP jest konstruktem SEKTORA KOSMOLOGICZNEGO, nie prymitywem
substratu.** Została wprowadzona, by zamrozić ewolucję pola tła (problem Ġ/G vs LLR) i dopasować
sektor ciemnej energii. Jej parametr był pierwotnie **skanowany przeciw danym DESI**, a jego
późniejsze „domknięcie" jest **tożsamością normalizacyjną**, nie wyprowadzeniem. Współczynnik
termiczny przy tym członie wynosi dziś **~10⁻³³ eV** — czyli przy skalach cząstek jest
**fizycznie nieobecny**.

## 1. Co dokładnie mówi rdzeń (cytaty)

**Definicja parametru** (`dodatekB`, thm:s0-from-GL, lin. 521–533):
> „przy warunku próżniowym β=γ i normalizacji Φ₀=1, **parametr entropii substratowej N0-7**
> wyraża się bezpośrednio przez masę pola: **s₀ = m_sp² = γ**, w jednostkach naturalnych
> **k_B T_sub|próżnia = 1**."

**Dowód** (lin. 535–564): z propagatora GL MFT ⟨δφ(k)²⟩ = k_BT/(KΦ₀²(k²+m_sp²)),
a następnie: „Zatem **(przy KΦ₀² = 1, k_B T_sub = 1)**: s₀ = m_sp² = γ."

**Funkcjonał entropii** (cor:entropy-potential, lin. 566–588):
> S_Γ(φ) = k_B N_B γ (φ − ln φ − 1),  V_eff^total = U(φ) − **T_Γ**·γ(φ − ln φ − 1),
> gdzie **T_Γ = T_H = H/(2π)** to temperatura Hawkinga-Gibbonsa **horyzontu Hubble'a**.
> „Dla T_Γ ≪ 1 (**dziś: T_H ~ 10⁻³³ eV**): m_eff ≈ m_sp."

**Skąd γ** (rem:O22-partial-closure, lin. 590–601):
> „s₀ **nie jest wolny** — jest zdeterminowany przez γ = m_sp², które jest już wyznaczone
> **z warunku ciemnej energii**: γ ≈ 12Λ_eff/Φ₀."

**Czym jest N0-7** (`dodatekA`, lin. 411–417):
> „N0-7 — **Zmienność G i ograniczenie kosmologiczne**: G(Φ)=G₀Φ₀/Φ (A8);
> Ġ/G = −ψ̇/ψ. Dla **zamrożonego pola** (|ψ̇| ≪ Hψ): |Ġ/G| ≲ H₀δ_ψ."

**Historia parametru** (`sek07_predykcje`, komentarz lin. 1051–1052):
> „O22 ZAMKNIĘTE v29: **skan s0 ∈ {0.001, 0.01, 0.05, 0.1}**, T0=0.1, OMEGA_M=0.315;
> s0=0.001: **w_de(z=0) = −0.82043**, Delta_w=+0.17957 **[POZA DESI]**."

## 2. Diagnoza — cztery ustalenia

### (a) Entropia narodziła się jako łatka na Ġ/G, nie jako prymityw
Nazwa parametru brzmi dosłownie „parametr entropii substratowej **N0-7**", a N0-7 to
**napięcie kosmologiczne**: zmienność G kontra ograniczenie obserwacyjne (LLR). Mechanizm
entropii daje polu **efektywną masę** m_eff² = γ(1+T_Γ), a więc **zamraża pole tła** —
czyli dokładnie to, czego wymaga N0-7 („dla zamrożonego pola |ψ̇| ≪ Hψ"). To jest
**stabilizacja zaburzeń**, dokładnie jak w Twojej intuicji.

### (b) s₀ był pierwotnie SKANOWANY przeciw danym kosmologicznym
Komentarz w `sek07_predykcje` pokazuje skan s₀ ∈ {0.001, 0.01, 0.05, 0.1} z oceną
w_de(z=0) względem **DESI**. To jest **dopasowanie parametru do obserwacji kosmologicznej** —
nie wyprowadzenie z reguł substratu.

### (c) „Domknięcie" s₀ = m_sp² = γ jest tożsamością normalizacyjną
Dowód jawnie **ustawia dwie stałe na 1** (k_B T_sub = 1 **oraz** KΦ₀² = 1) i dopiero wtedy
s₀ = 1/⟨δφ²⟩ = m_sp². Czyli treść to: *„odwrotność wariancji = masa² w jednostkach, w których
temperatura i sztywność wynoszą 1"*. To nie jest odkrycie związku entropia↔masa — to wybór
jednostek. **„s₀ nie jest wolny" jest prawdziwe tylko warunkowo**: s₀ zostało przywiązane
do γ, ale γ jest **wyznaczone przez ciemną energię** (γ ≈ 12Λ_eff/Φ₀). Parametr nie zniknął —
**przesunął się do sektora kosmologicznego**.

### (d) Kierunek wyprowadzenia jest ODWRÓCONY względem zamysłu programu
| Zamysł (Twój) | Stan faktyczny w rdzeniu |
|---|---|
| substrat → minimalne warunki → złożone obiekty | ciemna energia (Λ_eff) → γ → s₀ → „entropia substratu" |
| prymityw mikroskopowy | wielkość kalibrowana obserwacją kosmologiczną |
| — | temperatura substratu **T_Γ = temperatura horyzontu Hubble'a** |

Substrat nie dostarcza tu warunków — **substrat dostaje swoje parametry z kosmologii**.

## 3. Konsekwencja: dlaczego moja warstwa budżetowa musiała upaść

Budżet z `sek08c` to **B = N_B·s₀**, a s₀ pochodzi z konstrukcji opisanej wyżej. Czyli
budowałem lokalizację cząstek na wielkości, która:
- jest kalibrowana **ciemną energią**,
- ma współczynnik termiczny **T_H ~ 10⁻³³ eV** (przy skalach cząstek: **zero**),
- powstała, by **zamrozić kosmologiczną ewolucję G**.

To był **błąd kategorii**: import łatki kosmologicznej do mikrofizyki. Audytorzy złapali objawy
(nadokreślenie B=2, nieaddytywność, brak zasady wariacyjnej), ale **przyczyna źródłowa jest tutaj**.

## 4. Niuans, którego nie zamiatam (na korzyść rdzenia)

Faktyczny funkcjonał rdzenia **S_Γ ∝ (φ − ln φ − 1)** ma tę samą jakościową strukturę, którą
zgadywałem jako g(h)=h+1/h: wypukły, minimum w φ=1, rozbieżny w obie strony (φ→0: −ln φ→∞;
φ→∞: φ→∞). Więc **intuicja o „koszcie odchylenia od próżni" nie jest bezpodstawna** — jest
w rdzeniu. **Ale** w rdzeniu ten człon jest mnożony przez T_Γ ≈ 0, więc **nie ma mocy
lokalizującej**. Różnica między „kształt istnieje" a „kształt coś robi" jest tu decydująca.

## 5. Werdykt wobec hipotezy „błąd młodości"

**Potwierdzona, z zastrzeżeniem zakresu.**
- Jako **narzędzie sektora kosmologicznego** (zamrożenie pola, w_de, Ġ/G) entropia jest
  wewnętrznie spójna i ma jasną rolę — **to nie jest błąd w swoim sektorze**.
- Jako **prymityw substratu** — tak, to relikt: nazwa („entropia substratowa") sugeruje poziom
  fundamentalny, podczas gdy pochodzenie, kalibracja i skala są kosmologiczne.
- **Etykieta jest myląca i to ona wygenerowała mój błąd.** „Budżet informacyjny substratu"
  brzmi jak prymityw mikroskopowy; faktycznie jest parametrem tła kosmologicznego.

## 6. Co z tego wynika dla programu (bez rekomendacji działań — analiza)

1. **Nie ma dziś w TGP mikroskopowego prymitywu entropijnego.** Jest kosmologiczny człon
   termiczny o zerowym wpływie przy skalach cząstek.
2. Audytorska „droga ratunku" („wyprowadzić miarę zużycia z nieaddytywnej entropii, na którą
   rdzeń się powołuje") jest **iluzoryczna w tym sektorze** — ta nieaddytywna struktura nie
   istnieje jako obiekt mikroskopowy; rdzeń powołuje się na nią, by uzasadnić fh=1, ale jej nie podaje.
3. **Wciąż otwarte i nietknięte:** czy substrat ma jakikolwiek własny prymityw ilościowy
   (poza regułami kombinatorycznymi/topologicznymi, które JAKO JEDYNE przeżyły audyt).
4. Możliwe, że odpowiedź brzmi: substrat dostarcza **tylko** reguł kombinatorycznych
   (dopuszczalne kształty, dopuszczalne sklejenia), a wszystko ilościowe jest **relacyjne** —
   co byłoby spójne z Twoją ontologią energii relacyjnej i z tym, że **jedyne, co przetrwało
   audyt, to teoria grup i topologia**.

---

**Źródła (zweryfikowane u źródła):**
- `axioms/substrat/dodatekB_substrat.tex` lin. 518–601 (thm:s0-from-GL, cor:entropy-potential, rem:O22)
- `axioms/notacja/dodatekA_notacja.tex` lin. 411–417 (definicja N0-7)
- `core/sek07_predykcje/sek07_predykcje.tex` lin. 1051–1052 (skan s₀ vs DESI)
- `core/sek08c_metryka_z_substratu.tex` (def:info-budget — konsument s₀)
