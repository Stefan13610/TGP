# gs7 — Program badawczy: mikro-mechanizm a₀ w TGP

## Kontekst

Fenomenologia działa (gs1-gs2): a₀ = cH₀/(2π) daje BTFR, Freeman limit, phantom DM.
Wszystkie mechanizmy mikro z równania solitonu **OBALONO** (gs4-gs6):
- Ogony sin(r)/r: losowe fazy → kasują się (δ_rand/δ_N ~ 10⁻⁶ dla N=10¹¹)
- Self-terms |∇δ|²: już w zmierzonej masie
- Cross-terms: szum, nie sygnał
- Nieliniowość g'²/g: 12 rzędów za słaba (GM/Rc² ~ 5×10⁻⁷)

**Pytanie centralne**: jaki mechanizm TGP produkuje a₀ ≈ 1.2×10⁻¹⁰ m/s²?

## Plan badawczy — 5 opcji do zbadania

### Opcja A: Warunek brzegowy substratu

**Idea**: Substrat TGP jest skończonym, rozszerzającym się medium. Na granicy
przyczynowej (horyzont c/H₀) warunki brzegowe wpływają na propagację grawitacji.
To NIE jest efekt perturbacyjny — to globalna własność geometrii substratu.

**Co policzyć**:
1. Równanie pola TGP w rozszerzającym się substracie z warunkiem brzegowym na r = c/H₀
2. Rozwiązanie statyczne (potencjał galaktyczny) z tym warunkiem brzegowym
3. Porównanie z Newtonem: czy modyfikacja daje efektywne a₀?
4. Sprawdzenie czy a₀ = cH₀/(2π) wynika naturalnie

**Kluczowe pytanie**: Czy warunek brzegowy na horyzoncie kosmologicznym
"dociera" do skali galaktycznej? (W GR: tak — efekt de Sittera)

**Kryteria sukcesu**:
- [ ] Modyfikacja siły grawitacyjnej zaczyna się przy a ~ a₀
- [ ] a₀ ∝ cH₀ (proporcjonalność do tempa ekspansji)
- [ ] Poprawny profil krzywej rotacji (płaska, nie rosnąca)

---

### Opcja B: Grawitacja entropowa (Verlinde-like)

**Idea**: Grawitacja w TGP emerge z entropii substratu.
Miara entropii: S = ∫√g d³x (objętość substratu).
Siła entropowa: F = T · dS/dr.
"Temperatura" T ~ ℏH₀/(2πk_B) (temperatura horyzontu de Sittera).
→ Daje: a₀ ~ cH₀ (zgodne z obserwacją!)

**Co policzyć**:
1. Zdefiniować entropię substratu TGP: S[g] = ?
2. Obliczyć siłę entropową F = T·∂S/∂r dla profilu galaktycznego
3. Porównanie z MOND: czy daje tę samą interpolację?
4. Sprawdzenie Verlinde'owskiego wyprowadzenia: ΔS_elastic = (c³/6GℏH₀)·M·a₀·r²
5. Predykcje dla klastrów (Verlinde twierdzi, że rozwiązuje problem klastrów)

**Kluczowe pytanie**: Czy TGP ma naturalną definicję entropii substratu?
(g = metryka substratu → √g = element objętości → S ∝ V_substrate?)

**Kryteria sukcesu**:
- [ ] Entropia substratu dobrze zdefiniowana
- [ ] F_entropic dominuje nad F_Newton przy a < a₀
- [ ] Interpolacja MOND odtworzona (RAR)
- [ ] Predykcja dla klastrów lepsza niż MOND

---

### Opcja C: Zmodyfikowana dyspersja

**Idea**: Substrat TGP ma dyskretną strukturę na pewnej skali l.
Propagacja grawitacji ma zmodyfikowaną relację dyspersji:
ω² = c²k²·(1 - l²k² + ...)
Na długich falach (skala galaktyczna): dodatkowe termy modyfikują potencjał.
Jeśli l ~ c/H₀: modyfikacja zaczyna się przy a ~ cH₀ ~ a₀.

**Co policzyć**:
1. Potencjał z zmodyfikowaną dyspersją: Φ(r) = -GM/r · f(r/l)
2. Forma f(r/l) dla różnych dyspersji (Lorentz, sinus, lattice)
3. Wynikowa krzywa rotacji v²(r)
4. Porównanie z MOND i danymi SPARC

**Kluczowe pytanie**: Jaka skala l? Jeśli l ~ c/H₀ ~ 4 Gpc, to modyfikacja
na skali galaktycznej (kpc) jest perturbacyjna: r/l ~ 10⁻⁶.
Czy perturbacja rzędu (r/l)² ~ 10⁻¹² wystarczy?

**Ryzyko**: Podobne do wielu "emergent MOND" propozycji w literaturze
(Milgrom, Famaey, Blanchet). Trudno o oryginalność.

**Kryteria sukcesu**:
- [ ] f(r/l) daje płaską krzywą rotacji
- [ ] Skala l wynika z TGP (nie jest ad hoc)
- [ ] Predykcje różne od standardowego MOND (testowalność)

---

### Opcja D: Dwuskalowe TGP

**Idea**: Równanie solitonu ma DWIE naturalne skale:
1. Mikroskopowa: l_micro (rozmiar cząstki/solitonu)
2. Kosmologiczna: l_cosmo = c/H₀ (horyzont ekspansji)

Przyspieszenie de Sittera:
a_dS = c·√(Λ/3) = c·H₀·√(Ω_Λ) = 5.4×10⁻¹⁰ m/s²

**a₀ ≈ a_dS/(2π) ≈ 8.6×10⁻¹¹ m/s²** — bliskie obserwowanemu 1.2×10⁻¹⁰!

**Co policzyć**:
1. Skąd czynnik 2π? (okres oscylacji solitonu? geometria sfery?)
2. Dokładna relacja: a₀ = f(c, H₀, Ω_Λ) — jaka forma f?
3. Czy w TGP: Λ = μ²c² (masa efektywna pola) → a₀ = μc²/(2π)?
4. Ewolucja a₀(z) z redshiftem (a₀ zmienia się z H(z)?)
5. Porównanie z danymi: a₀(z=0) vs a₀(z=1) — czy obserwacje ograniczają?

**Kluczowe pytanie**: Dlaczego a₀ ≈ cH₀ a nie c²/l_Planck czy inna kombinacja?
W TGP: jeśli substrat rozszerza się z tempem H₀, to naturalna skala
przyspieszenia = c·H₀. Czynnik 1/(2π) z cykliczności?

**Kryteria sukcesu**:
- [ ] Czynnik 2π wyprowadzony (nie ad hoc)
- [ ] a₀(z) predykcja zgodna z obserwacjami (lub dająca się testować)
- [ ] Związek z Ω_Λ wyjaśniony (dlaczego a_dS a nie a_H = cH₀?)

---

### Opcja E: μ(r) zależne od skali

**Idea**: Równanie solitonu TGP na skalach galaktycznych ma INNY współczynnik:
```
g'' + g'²/g + 2g'/r + μ²(g - 1) = source
```
- μ² = 1 na skali mikroskopowej (cząstka/soliton)
- μ² = (a₀·r/(c²))² ~ 10⁻¹² na skali galaktycznej

Wtedy długość fali ogona: λ = 2π/μ ~ 10⁶ × r_gal
→ Ogon NIE oscyluje na skali galaktycznej → NIE kasuje się!
→ Koherentna modyfikacja potencjału grawitacyjnego

**Co policzyć**:
1. Rozwiązanie g'' + g'²/g + 2g'/r + μ²(g-1) = source dla μ ≪ 1
2. Forma potencjału: Φ(r) = -GM/r · (1 + f(μr))
3. Krzywa rotacji: v²(r) = GM/r · (1 + r·f'(μr)/f + ...)
4. Warunek na μ aby v²→const przy dużych r → daje a₀?
5. Skąd μ(r)? Renormalizacja grupowa? Efektywna teoria pola?

**Kluczowe pytanie**: Czy μ² ≪ 1 na dużych skalach wynika z TGP,
czy jest parametrem ad hoc? Jeśli μ² = (H₀/ω_soliton)², to:
μ = H₀/ω_sol ~ 10⁻¹⁸/10⁻⁴² = 10²⁴ (za duże!)
Potrzeba innej motywacji dla μ ≪ 1.

**Kryteria sukcesu**:
- [ ] Płaska krzywa rotacji z μ ≪ 1
- [ ] μ wynika z TGP (nie ad hoc)
- [ ] Reprodukuje BTFR (v⁴ = GMa₀)
- [ ] Zgodne z gs1/gs2 fenomenologią

---

## Kolejność badań

1. **Opcja D** (dwuskalowe TGP) — szybki test: sprawdzić skąd 2π, ile wynosi a_dS/(2π)
2. **Opcja E** (μ zależne od skali) — najbardziej konkretna, można od razu policzyć
3. **Opcja A** (warunek brzegowy) — wymaga rozwiązania równania z BC na horyzoncie
4. **Opcja B** (entropia) — wymaga zdefiniowania entropii substratu
5. **Opcja C** (dyspersja) — ryzyko duplikacji istniejącej literatury

## Kryteria odrzucenia

Mechanizm jest **odrzucony** jeśli:
- Daje a₀ różne od cH₀/(2π) o więcej niż rząd wielkości
- Nie daje wykładnika 4 w BTFR (v⁴ ∝ M)
- Wymaga fine-tuningu parametrów
- Jest wewnętrznie sprzeczny z innymi wynikami TGP (solitony, kosmologia)

## Powiązania

- [[TGP/TGP_v1/research/galaxy_scaling/gs6_assessment.py]] — ocena obalonych mechanizmów
- [[TGP/TGP_v1/research/galaxy_scaling/gs1_flat_well_model.py]] — fenomenologia (cel do reprodukcji)
- [[TGP/TGP_v1/research/cosmo_tensions/ct6_mechanism_diagnosis.py]] — analogiczny problem ze skalami
