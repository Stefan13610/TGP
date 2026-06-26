---
title: "op-PSR-spin-substrate-channel-2026-06-13 (EXPLORATORY)"
folder_status: exploratory-note
claim_status: "STRUCTURAL_CONDITIONAL — spin does NOT change the orbital verdict"
date: 2026-06-13
relates_to: "[[../op-PSR-Pdot-energy-balance-2026-06-13/]] (FALSIFIED — NOT reopened)"
authorization: "User 2026-06-13: prośba o policzenie wkładu spinu pulsara, nawet jeśli pogarsza obraz"
---

# Wkład spinu pulsara do bilansu energii (eksploracja warunkowa)

**Status metodologiczny:** to NIE jest recovery zamkniętego, sfalsyfikowanego cyklu
op-PSR-Pdot (zakazane, Phase 0 §6.6). To oddzielna, eksploracyjna odpowiedź na
pytanie autora, czy spin pulsara — w obrazie TGP, gdzie obrót to rezonans z
substratem Φ — otwiera kanał energetyczny zmieniający R = P_φ/P_GR. Liczby kanału
radiacyjnego są warunkowe (mechanizm rezonansu ruchowego = STRUCTURAL_CONDITIONAL).

## Pytanie autora
Pęd/obrót w TGP = rezonans z otaczającą przestrzenią; przy pulsarach nałożenie
efektu orbitalnego i obrotowego może być silne. Policzyć, nawet jeśli pogarsza.

## Wynik (skrypt `spin_energy_balance.py`)

| # | Test | Wynik |
|---|---|---|
| P1 | Spin osiowosymetryczny, stały: czy M⃛_ij ≠ 0? (warunek radiacji) | **M⃛_ij = 0** (symbolicznie) → zero radiacji do substratu |
| P1 | Spin z elipsoidalnością ε | M⃛_xx = 8 I_z Ω³ ε sin(2Ωt) → radiacja przy **2Ω_spin** |
| P2 | P_GR orbitalne / P_φ = 1/6 P_GR (J0737) | 2.36×10²⁵ W / 3.94×10²⁴ W |
| P2 | Luminancja spin-down A | 7.7×10²⁶ W (→ magnetosfera/EM, nie orbita) |
| P3b | Odcisk spin-orbita (1.5PN) na Ṗ_b | ~5×10⁻¹³ × kwadrupol → pomijalny |
| P3c | Rezerwuar spin-down vs offset potrzebny | 196× większy, ale **odłączony** od orbity |
| P4 | Nowy kanał TGP: skalar przy 2Ω_spin = 88 Hz | ε_sat ~ 7×10⁻⁷; osobna obserwabla (CW), nie wchodzi w Ṗ_b |

## Werdykt (warunkowy)
Włączenie spinu **nie ratuje** werdyktu orbitalnego: R pozostaje 1/6 (Gałąź A).
Dwie przyczyny strukturalne:
1. **Symetria (mechanism-independent):** stały, osiowosymetryczny obrót ma w
   układzie LAB stacjonarny rozkład masy → zerowy zmienny multipol → zero
   strumienia energii do nieskończoności, w OTW i w TGP jednakowo. Rezonans z
   substratem, jeśli realny, jest tu polem **bliskim/stacjonarnym**
   (analog gravitomagnetyzmu / frame-dragging), nie kanałem radiacyjnym.
2. **Rozdzielczość spektralna:** kanał orbitalny (2f_orb ≈ 0.23 mHz) i spinowy
   (2ν_spin ≈ 88 Hz) dzieli ~6 rzędów. W strumieniu (∝⟨ψ̇²⟩) człony krzyżowe
   między nimi uśredniają się do zera → **brak konstruktywnego nałożenia**,
   wbrew intuicji "silnego nałożenia".

**Produkt uboczny (pozytywny):** spin surfaces a genuine NEW TGP prediction —
skalarna składowa promieniowania spin-down przy 2Ω_spin (P4). To osobny przyszły
test (continuous-wave / budżet spin-down), nie ratunek Ṗ_b.

## Rozszerzenie 2026-06-13b: całkowanie po nasyceniu (skrypt `saturation_charge_integration.py`)

Pytanie autora: skoro w silnym źródle konfiguracja przestrzeni nasyca się (jak przy
horyzoncie ψ=4/3, gdzie brak rekonfiguracji i oddziaływań), to całkowanie po
„nakładających się warstwach materii" nie daje liniowego q·M — głębokie warstwy
ekranowane. Czy efektywny ładunek skalarny C_eff < q·M tłumi P_φ i ratuje werdykt?

Równania TGP użyte: horyzont ψ=4/3 (sek08c); profil złożony χ_glob i nasycenie
χ~a₀/r M-niezależne (sek08_formalizm 2466); usztywnienie f(Φ)=1+4ln χ (2521);
jedno sprzężenie q²=4πG wiąże ORBITĘ i promieniuje (Pdot T2a LOCKED).

| # | Test | Wynik |
|---|---|---|
| P1 | NS pod horyzontem? ψ_surf=1.082, ψ_core≤1.329 vs ψ_hor=4/3=1.333 | **SUB-horyzont** — brak strefy „zamrożonej"; pulsar oddziałuje w pełni |
| P2 | Siła nasycenia f(Φ)=1+4ln χ w rdzeniu | tylko **~1.3–1.9×** (O(1)), NIE reżim χ≫1 (M-niezależny) |
| P3 | **dR/dS dla wspólnego przeskalowania q→q(1−S)** | **R = 1/6 niezmiennie; dR/dS = 0** (symbolicznie) |
| P4 | Różnicowe ekranowanie S_A≠S_B (różna kompaktowość) → dipol −1PN | P_dip/P_quad ~ **11×** → **POGARSZA** werdykt |

**Werdykt (warunkowy):** całkowanie pulsara jest w TGP sensowne, ale nie ratuje
falsyfikacji — z powodu **jednopolowości**. To samo q²=4πG wiąże orbitę i
promieniuje; nasycenie przeskalowuje OBA sektory tym samym czynnikiem, a obserwabla
R = P_φ/P_GR (znormalizowana do OTW z mas wyznaczonych z orbity) **kasuje wspólne
przeskalowanie dokładnie** (dR/dS=0). Jedyny efekt różnicowy — dipol z S_A≠S_B —
**pogłębia** falsyfikację (−1PN, ~11× kwadrupol). Tłumienie radiacji BEZ tłumienia
wiązania to Gałąź D (screening) = samobójcza. Brak wyjścia w teorii jednopolowej.
Spójne z W-PDOT-4 (HIGH).

## Rozszerzenie 2026-06-13c: odpowiedź rezonansowa OTOCZENIA (skrypt `environment_response.py`)

Pytanie autora: poprzednie rachunki liczyły tylko własny układ źródła (naga emisja
dwuciałowa). Czy policzono odpowiedź **ośrodka** — substratu skonfigurowanego przez
masę (kaskada „masa→konfiguracja→kolejna konfiguracja")? Substrat jako medium:
`((1/c²)∂_tt − ∇² + γ)δΦ = źródło`, dyspersja `ω²=c²k²+c²γ`, jedyna własna częstość
= `ω_gap = c√γ = H₀` (m_sp ~ H₀/c).

| # | Test | Wynik |
|---|---|---|
| A | Przezroczystość: ω_GW vs ω_gap | ω_GW/ω_gap ≈ **6×10¹⁴** → v_g=c, brak absorpcji/dyspersji; korekta ~10⁻²⁹ |
| B | Wake (Czerenkow): Mach = v_orb/c_s | **2×10⁻³ ≪ 1** → podświetlnie, brak kanału wake/rezonansu ruchowego |
| C | Rezonans wsteczny otoczenie→orbita | detuning ω_orb/ω_gap ≈ **3×10¹⁴** → zero anti-dampingu |
| D | Co mierzy Ṗ_b? | energię **opuszczającą orbitę** (near-zone), nie docierającą do ∞ → absorpcja w dół strumienia nieistotna |
| E | Dressing tła na R | tło działa na skalar i tensor jednakowo (oba to δΦ) → kasuje się w R |

**Werdykt (warunkowy):** włączenie odpowiedzi otoczenia NIE zmienia werdyktu.
Ośrodek jest przy ω_orb przezroczysty, podświetlny i nierezonansowy. Decydujący
argument (D): Ṗ_b mierzy energię opuszczającą orbitę, więc cokolwiek substrat robi
w dół strumienia (pochłania, rekonfiguruje kaskadowo) nie rusza Ṗ_b; jedyne wyjście
to rezonansowe pompowanie energii Z POWROTEM do orbity, a substrat nie ma modu przy
ω_orb (jedyny mod = H₀, 14 rzędów dalej). Kaskada konfiguracji to statyczny profil
χ_glob(r) — stacjonarny, nie promieniuje; jako tło kasuje się w stosunku R.

## Rozszerzenie 2026-06-13d: ujęcie 3-polowe / n-body (skrypt `threefield_nbody.py`)

Pytanie autora: ująć to 3-polowo — P1=skalar ciała 1, P2=skalar ciała 2, P3=pole
otoczenia (tło wszechświata) — i policzyć złożenia P3-P1, P3-P2, P3-P1-P2.

**Odczyt A (w TGP, jedno pole; P1/P2/P3 = wkłady δΦ₁+δΦ₂+Φ̄):**

| # | Wynik |
|---|---|
| P1 sym. | **R = P_φ/P_GR = 1/6 dla DOWOLNEGO kwadrupola I⃛_ij** — stały czynnik geometryczny, niezależny od dekompozycji źródeł |
| P2 n-body | `<I⃛²>`: P3-P1 = 32d⁴m₁²m₂⁴ω⁶/M⁴; P3-P2 = 32d⁴m₁⁴m₂²ω⁶/M⁴; cross = 64d⁴m₁³m₂³ω⁶/M⁴; **full = 32μ²d⁴ω⁶** (tożsamość full=self₁+2cross+self₂ ✓) |
| | R = 1/6 dla P3-P1, P3-P2 ORAZ P3-P1-P2 — człon krzyżowy BUDUJE kwadrupol układu (już w sfalsyfikowanym rachunku), nie dodaje kanału |
| P3 | tło Φ̄ ewoluuje w skali 1/H₀; ∂_tΦ̄ przy ω_orb ~ 2×10⁻¹⁴ → radiacyjnie bezwładne |

**Odczyt B (trzy NIEZALEŻNE pola dynamiczne = nowa teoria):**
- 3 skalary (helicity-0): każdy promieniuje kwadrupol z tym samym 1/6 → R_total = n/6 ≥ 1/6, **pogarsza**.
- Ṗ_b zgadza się z kwadrupolem **helicity-2** do 6×10⁻⁵; suma pól helicity-0 NIE syntezuje helicity-2 (teoria reprezentacji Lorentza). Dopasowanie do OTW wymaga, by P3 było realnym polem **spin-2 = grawiton**.
- FOUNDATIONS reguła 6 **zakazuje grawitonu** (fundamentalnego i kompozytowego). Więc jedyna struktura wielopolowa, która ratuje, leży **poza aksjomatami TGP** = nowa teoria z własnym Phase 0 + re-derywacją Newton/1PN/GWTC-3/NICER.

**Werdykt:** w obu odczytach werdykt orbitalny stoi. Pomysł 3-polowy, doprowadzony
uczciwie do końca, **re-derywuje** dlaczego jedyne wyjście to porzucenie rdzenia
TGP (jedno pole skalarne, brak grawitonu). To pozytywny wynik: ostrzy, czym MUSI
być amendment — polem spin-2 — a nie rekombinacją n-body w obrębie TGP.

## Rozszerzenie 2026-06-13e: znak korelacji źródeł (skrypt `correlation_signs.py`)

Zarzut autora: czynniki nie powinny się jawnie DODAWAĆ — korelacja źródeł powinna
ODEJMOWAĆ. Test: jawny znak członu krzyżowego 2⟨X₁:X₂⟩ per multipol (J0737, CM).

| Multipol | człon krzyżowy 2⟨X₁:X₂⟩ | suma | parzystość w x |
|---|---|---|---|
| Monopol q·m | (Q̇=0, brak radiacji) | — | — |
| **Dipol** q·m·x | **−2 d²m₁²m₂²ω⁴q²/M²** (UJEMNY) | **= 0** | nieparzysty → anty-równoległe → **ODEJMUJE** |
| **Kwadrupol** m·xx | **+64 d⁴m₁³m₂³ω⁶/M⁴** (DODATNI); korelacja znormalizowana = **+1** | ≠ 0 | parzysty → równoległe → **DODAJE** |

**Rozstrzygnięcie:** autor ma rację — korelacja ODEJMUJE, ale **tylko dla momentów
nieparzystych w x** (dipol, pęd). Przeciwne strony CM dają przeciwne znaki →
anty-równoległe → kasacja. I TGP **już to bankuje**: skalarny dipol D = q·Σmᵢxᵢ = 0
(T0). Dla kwadrupola (∝ x·x, parzysty) znak x się podnosi do kwadratu → oba tensory
są **idealnie wyrównane** (korelacja +1) → dodawanie jest **wymuszone geometrią**,
nie wyborem modelu. To nie da się obejść dekompozycją P1/P2/P3. Aby kwadrupol
ODEJMOWAŁ, potrzeba drugiego źródła kwadrupolowego w przeciw-fazie tensorowej
niosącego ~P_GR = pole spin-2 = grawiton (zakazany). Werdykt stoi; R = 1/6.
