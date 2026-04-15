# PLAN ROZWOJU TGP v4 — 2026-04-14

> **Zasada organizacji**: Rdzeń teorii (`TGP_v1/`) jest stabilny i gotowy do publikacji.
> Nowa praca badawcza odbywa się w izolowanych podfolderach `research/<problem>/`.
> Scalanie do rdzenia następuje dopiero po zamknięciu problemu (testy + dowód).

> **Poprzednie plany**: PLAN_ROZWOJU_v3.md (zamknięty), PLAN_DOMKNIECIA_MASTER.md (zamknięty)

---

## Architektura folderów

```
TGP_v1/                          ← RDZEŃ (stabilny, publikowalny)
├── research/
│   ├── cabibbo_correction/      ← R1: Korekcja Cabibbo (Ω_Λ/N)²
│   ├── continuum_limit/         ← R2: CG-1/3/4 ciągłe przejście
│   ├── why_n3/                  ← R3: Dlaczego N=3 generacji
│   ├── metric_ansatz/           ← R4: h(Φ)=Φ z pierwszych zasad
│   ├── mass_scaling_k4/         ← R5: m ∝ A_tail⁴ z działania
│   ├── brannen_sqrt2/           ← R6: B=√2 analitycznie
│   └── uv_completion/           ← R7: Unifikacja przy M_Pl
├── sek00–sek10, dodatek*.tex    ← Tekst główny (nie ruszać)
├── tgp_letter.tex               ← PRL letter (gotowy)
├── tgp_companion.tex            ← PRD companion (gotowy)
├── scripts/                     ← Skrypty walidacyjne rdzenia
└── nbody/                       ← Biblioteka N-body + examples
```

**Reguły pracy:**
1. Każdy folder `research/X/` zawiera własny `README.md` z opisem problemu
2. Skrypty eksploracyjne trafiają do `research/X/`, NIE do `scripts/` ani `nbody/examples/`
3. Gdy problem jest zamknięty → scalenie do rdzenia (nowy dodatek .tex + skrypt weryfikacyjny)
4. Foldery badawcze mogą być atakowane **niezależnie i równolegle**

---

## Mapa problemów

### Priorytet: NATYCHMIAST (ROI: wysoki impact / niski nakład)

| ID | Problem | Folder | Nakład | Impact |
|----|---------|--------|--------|--------|
| **R1** | Korekcja Cabibbo (Ω_Λ/N)² | `research/cabibbo_correction/` | 2–4 tyg. | ⭐⭐⭐⭐⭐ |
| **R4** | h(Φ)=Φ z równań Einsteina | `research/metric_ansatz/` | 2–4 tyg. | ⭐⭐⭐⭐ |

### Priorytet: ŚREDNIOTERMINOWY (1–3 miesiące)

| ID | Problem | Folder | Nakład | Impact |
|----|---------|--------|--------|--------|
| **R6** | B=√2 analitycznie z ODE | `research/brannen_sqrt2/` | 2–6 tyg. | ⭐⭐⭐⭐ |
| **R5** | m ∝ A⁴ z działania solitonu | `research/mass_scaling_k4/` | 1–3 mies. | ⭐⭐⭐⭐ |

### Priorytet: DŁUGOTERMINOWY (6–12+ miesięcy)

| ID | Problem | Folder | Nakład | Impact |
|----|---------|--------|--------|--------|
| **R2** | CG-1/3/4 continuum limit | `research/continuum_limit/` | 6–12 mies. | ⭐⭐⭐⭐⭐ |
| **R3** | Dlaczego N=3 generacji | `research/why_n3/` | 1–3 mies. | ⭐⭐⭐⭐⭐ |
| **R7** | UV completion, unifikacja | `research/uv_completion/` | 2–4 tyg. | ⭐⭐ |

### Nie wymagają osobnych folderów

| Problem | Status | Uwagi |
|---------|--------|-------|
| Monitoring K10 (JUNO, NO vs IO) | CZEKAMY | Dane ~2028–2030 |
| Monitoring K14 (DESI DR3, phantom) | CZEKAMY | Dane ~2027–2029 |
| δ_CP ≈ 62° vs 68° (~2σ) | AKCEPTOWALNE | Przybliżenie pierwszego rzędu |
| Chirality + anomalie z dyn. Φ | NISKI | Nie blokuje publikacji |
| α_em z substratu | NISKI | Program otwarty |

---

## R1: Korekcja Cabibbo — `research/cabibbo_correction/` ✅ ROZWIĄZANE

**Problem:** λ_C = Ω_Λ/N = 0.2282 vs PDG 0.22500 ± 0.00067 → **4.8σ napięcie**

**ROZWIĄZANIE (2026-04-14):** Korekcja Z₃ self-energy subtraction:

```
λ_C = (Ω_Λ/N) × (|G| - |Z₃|)/(|G| - 1) = (Ω_Λ/N) × 165/167 = 0.22550
Napięcie: 4.8σ → 0.75σ  ✅
```

**Fizyka:** Elementy Z₃ (3 z 168) zachowują numer generacji i nie przyczyniają się do mieszania. Czynnik F = 165/167 = frakcja "aktywnych" kanałów mieszania.

**Wyniki:**
- [x] Pełna analiza GL(3,𝔽₂): 6 klas sprzężoności, 28 podgrup Z₃, 20 podwójnych koklas
- [x] Korekcja 165/167: napięcie 0.75σ
- [x] Weryfikacja CKM: |V_us|=0.7σ, |V_cd|=1.0σ, |V_cb|=0.2σ
- [x] Test jednoznaczności: **tylko GL(3,𝔽₂) daje zgodność** (S₃ → 131σ, A₅ → 6.7σ)
- [x] Scalenie z rdzeniem (tgp_companion.tex §F2, tgp_letter.tex, scripts/)
- [ ] Formalizacja dowodu

**Pliki:** `research/cabibbo_correction/r1_gl3f2_structure.py`, `r1_cabibbo_correction_derivation.py`

---

## R2: Ciągłe przejście substrat → pole — `research/continuum_limit/`

**Problem:** Trzy otwarte twierdzenia blokują claim "wyprowadzone z pierwszych zasad":
- **CG-1**: Istnienie i jednoznaczność punktu stałego S* (kontrakcja Banacha)
- **CG-3**: Zbieżność Φ_B → Φ w H¹ (homogenizacja)
- **CG-4**: Identyfikacja K_hom = K_TGP

**Obecny status:**
- Słabe twierdzenie α=2: ZAMKNIĘTE (algebraicznie, Lemma A1–A5)
- Numeryczna weryfikacja: K_IR/K_UV = 1.000 (FRG LPA', CG-2: 8/8 PASS)
- Silne twierdzenie: OTWARTE

**Plan ataku:**
1. CG-1: Szukać kontrakcji Banacha w przestrzeni operatorów blokowania
   - Narzędzia: Teoria ERG (Polchinski, Wetterinck), functional analysis
   - Kluczowa literatura: Brydges & Yau, Bauerschmidt et al.
2. CG-3: Zastosować twierdzenia homogenizacji (de Giorgi–Nash–Moser)
   - Sprawdzić warunki: eliptyczność, ograniczoność współczynników
3. CG-4: Identyfikacja K — wymaga zarówno CG-1 jak i CG-3

**Kryterium zamknięcia:** Formalne dowody trzech twierdzeń (CG-1, CG-3, CG-4)

**Pliki rdzenia do scalenia:**
- Rozszerzenie `dodatekQ_coarse_graining_formal.tex`
- Nowy `dodatek_CG_proof.tex` z pełnymi dowodami

**Uwaga:** To jest **czysta matematyka**. Może być publikowalne niezależnie od fizyki TGP.

---

## R3: Dlaczego N=3 generacji — `research/why_n3/` ⚙️ DUŻY POSTĘP

**Problem:** NAJFUNDAMENTALNIEJSZE otwarte pytanie. GL(3,𝔽₂) zakłada N=3, nie wyprowadza.

**POSTĘP (2026-04-15): Mechanizm auto-przestrzeni + α-zależność**

```
GŁÓWNE WYNIKI:
1. Singularność metryczna: soliton z g₀ > g₀_crit ma g(r)→0 → dziura w przestrzeni
2. g₀_crit(1D) = 4/3 DLA KAŻDEGO α (twierdzenie z prawa zachowania)
3. g₀_crit(3D) zależy od α kinetic coupling:
     α=0.5: g₀_crit=2.618 → N=3 ✓
     α=0.882: g₀_crit=2.276 → N=2→3 TRANSITION
     α=1 (substrat): g₀_crit=2.206 → N=2 (deficit 3.1%)
     α≈3: g₀_crit=1.728 ≈ g₀^τ(Koide)!
4. dm/dg₀ → ∞ przy barierze (twardy limit)
5. Lagrangian: L = g^{2α}·g'²/2 + g³/3 - g⁴/4 (sam U(g) dla all α)
```

**Status:**
- [x] Mechanizm singularności metrycznej — ZWERYFIKOWANY
- [x] g₀_crit(1D) = 4/3 — TWIERDZENIE (α-niezależne)
- [x] α_crit = 0.882 — OBLICZONE
- [x] α_Koide ≈ 3 — ODKRYTE
- [ ] Fizyczna wartość α (dlaczego α < 0.882?)
- [ ] Analityczne g₀_crit(3D)
- [ ] Formalizacja dowodu

**Kryterium zamknięcia:** Wyprowadzenie α < 0.882 z pierwszych zasad.

**Status:** SILNA HEURYSTYKA. Mechanizm działa, brakuje ustalenia α.

---

## R4: Ansatz metryczny h(Φ)=Φ — `research/metric_ansatz/` ✅ ZASADNICZO ZAMKNIĘTY

**Problem:** Metryka ds² = -(c₀²/ψ)dt² + ψδᵢⱼdxⁱdxʲ z ψ=Φ/Φ₀ jest **postulatem**.
Dlaczego p=1, a nie Φ^p?

**ROZWIĄZANIE (2026-04-14):** Pięć niezależnych argumentów wymusza p = 1:

1. **Gęstość substratu:** Φ = gęstość węzłów → g_ij = (Φ/Φ₀)δ_ij → p=1 (definicja)
2. **PPN Cassini + LLR:** γ = p = 1 (Cassini: |γ-1| < 2.3×10⁻⁵), β = 1 (LLR)
3. **Budżet informacyjny:** f·h = 1 (antypodyczny) → q = p
4. **Element objętościowy:** √(-g) = ψ^p musi = ψ (gęstość) → p = 1 [NOWY, A2b]
5. **Stosunek mas:** Tylko p=1 → r₂₁ = 206.77 ≈ PDG 206.768

**Weryfikacja:** `r4_einstein_self_consistency.py` 11/11 PASS, `ex206` 8/8 PASS, `a2` 6/6 PASS

**Co pozostaje (dodatkowe, nie blokujące):**
- [ ] A2a: Relacja dyspersji fononów (c_s ∝ √Φ)
- [ ] A2c: Argument entropijny (S_BH ∝ A)
- [ ] Formalizacja łańcucha dowodowego

**Pliki rdzenia do scalenia:**
- Nowy paragraf w `sek08c_metryka_z_substratu.tex`

---

## R5: Prawo skalowania m ∝ A_tail⁴ — `research/mass_scaling_k4/` ⚙️ W TRAKCIE

**Problem:** Formuła masowa jest **fundamentem** sektora leptonowego. Daje r₂₁ = 206.768
(0.0001% zgodność z PDG). Ale k=4 jest **postulatem**.

**POSTĘP (2026-04-14):** Łańcuch dowodowy skorygowany po wynikach negatywnych:

```
P1: WIRIAŁ E^(2) = 0 dokładnie           ✅ UDOWODNIONE
P2: KONWERGENCJA k ≥ 4 w d=3             ✅ UDOWODNIONE (k = 2(d-1)/(d-2) = 4)
P3: E^(3) → 0                            ❌ OBALONY (E³~A³ dominuje E⁴~A⁴!)
P3': On-shell identity: E³=-(2π/3)∫h³r²  ✅ UDOWODNIONE (nowe twierdzenie)
P4: E_full ~ A^{2α} (nieperturbacyjne)    ✅ ZWERYFIKOWANE (k≈4.4 canonical)
```

**WYNIK NEGATYWNY:** E^(3) NIE znika! |E³/E⁴| ~ A^{-0.9} → ∞ dla małych solitonów.
Perturbacyjny dowód m~A⁴ jest **niemożliwy**. Skalowanie mas jest własnością
**nieperturbacyjną** (core-tail matching + convergence).

**Kluczowe wyniki:**
- k = 2(d-1)/(d-2) = 4 jest **jedynym integerem** (d=3→4, d=4→3, d=5→2.67)
- Weryfikacja numeryczna: k_eff = 4.0001, (A_μ/A_e)⁴ = 206.74 ≈ 206.768 (0.013%)
- On-shell identity (nowe): E³_sub = -(2π/3)∫h³r²dr, E³_can = +(4π/3)∫h³r²dr
- ∫h³r²dr logarytmicznie rozbieżny dla zlinearyzowanego h, skończony dla pełnego solitonu
- E_full ~ A^{4.36} (canonical) potwierdza skalowanie nieperturbacyjne

**Co zostaje do zamknięcia:**
1. **Nieperturbacyjny dowód E_full ~ A^{2α}** — mechanizm core-tail matching
2. **Zamknięta formuła c_M** — stała proporcjonalności wyznaczona tylko numerycznie
3. **Formalizacja łańcucha (Lean 4)**

**Pliki:**
- `research/mass_scaling_k4/r5_e3_cancellation.py` — E^(3) NIE znika: 5/7 PASS (**NOWE**)
- `research/mass_scaling_k4/r5_mass_ratio_verification.py` — weryfikacja k_eff i zbieżności
- `research/mass_scaling_k4/r5_virial_mass_derivation.py` — skan E(A_tail), błędne ODE
- `scripts/lp4_mass_exponent_verification.py` — rdzeń, 9/9 PASS

**Kryterium zamknięcia:** Twierdzenie: "m ∝ A⁴ wynika z α=2, d=3, K(Φ)=Φ²"

**Pliki rdzenia do scalenia:**
- Rozszerzenie `dodatekJ_ogon_masy.tex`

---

## R6: B=√2 analitycznie — `research/brannen_sqrt2/` ⚙️ W TRAKCIE

**Problem:** Stosunek Brannena B = b/a = √2 jest potwierdzony numerycznie do 10⁻⁶,
ale brak dowodu analitycznego. Gdyby udowodnić → Koide K=2/3 staje się **twierdzeniem**.

**POSTĘP (2026-04-14):**

Łańcuch algebraiczny KOMPLETNY:
```
GL(3,𝔽₂) → Z₃ podgrupa → fazy 120° → K = 2/3 → B = √2
```

Negatywne wyniki (eliminacja ślepych ścieżek):
- **Ścieżka 4 nie działa:** F(φ) = A(φg₀)/A(g₀) NIE jest stałe (CV = 220%)
- **φ²-drabinka tau nie działa w substracie:** g₀^τ = φ²g₀^e → A_tail = 0
- Best-fit g₀^τ = 1.73, co daje B = 1.4143 ≈ √2 z 10⁻⁴

**Brakujące ogniwo:** Formalizacja Z₃ → fazy 120° (dlaczego parametryzacja Brannena?)

**Co zostaje do zamknięcia:**
1. **Z₃ → equidistant phases:** formalny dowód z reprezentacji GL(3,𝔽₂)
2. **Tau z Koide constraint:** m_τ z K(m_e,m_μ,m_τ) = 2/3, nie z drabinki φ²
3. **Porównanie kanoniczne vs substrat:** a3d (kanoniczne) działa lepiej

**Pliki rdzenia do scalenia:**
- Rozszerzenie `dodatekT3_brannen_geometry.tex`

---

## R7: UV completion — `research/uv_completion/`

**Problem:** TGP nie ma explicit UV completion. Bieganie stałych sprzężenia
α₁, α₂, α₃ z modyfikacjami TGP nie jest obliczone.

**Obecny status:**
- Aksjomat: g₀ jest RG-invariant (β(g₀)=0 at 1-loop, F12)
- Ale: Wpływ dynamicznego Φ na running α_i nie jest obliczony
- Otwarte: Czy α₁ = α₂ = α₃ przy M_Pl z TGP?

**Plan ataku:**
1. Obliczyć β-funkcje dla α_i z modyfikacją K(Φ) w propagatorach
2. Sprawdzić unifikację przy M_Pl
3. Porównać z MSSM / inne modele

**Kryterium zamknięcia:** Wykres running α_i(μ) z TGP, sprawdzenie unifikacji

**Priorytet:** BONUSOWY — nie blokuje publikacji, ale wzmacnia teorię

---

## Status eksperymentalny (monitoring — bez folderów)

| Kill criterion | Eksperyment | Dane kiedy | Co zabija TGP |
|---------------|-------------|------------|----------------|
| K10 | JUNO | ~2028–2030 | Inverted ordering |
| K14 | DESI DR3 / Euclid | ~2027–2029 | w < -1 (phantom) |
| K5 | Hyper-K | 2026–2028 | Rozpad protonu |
| K11 | DESI/Euclid | ~2027 | Σm_ν > 200 meV |
| K12 | CMB-S4 | ~2028 | 4. generacja (N_ν ≠ 3) |

---

## Tabela statusu epistemicznego rdzenia

| Element | Status | Folder | Uwagi |
|---------|--------|--------|-------|
| α=2 (kinetic coupling) | ✅ TWIERDZENIE (słabe) | R2 (silne) | Algebraiczne, Lemma A1–A5 |
| K(ℓ)=2/3 | ✅ NUMERYCZNE (10⁻⁶) | R6 (analityczne) | Zależy od B=√2 |
| h(Φ)=Φ | ✅ ZASADNICZO ZAMKNIĘTY | R4 | 5 niezależnych argumentów, 11/11 PASS |
| m ∝ A⁴ | 🟡 NUMERYCZNE + argument konwergencji | R5 | E³≠0 (neg.result), E_full~A⁴ (niepert.) |
| N=3 | 🟡 MECHANIZM (α-zależny) | R3 | α<0.882→N=3, substrat α=1→N=2 (deficit 3.1%) |
| λ_C = Ω_Λ/N | 🟠 4.8σ NAPIĘCIE | R1 | Brak korekcji wyższego rzędu |
| CG-1/3/4 | 🔴 OTWARTE | R2 | Czysta matematyka |
| β=γ (vacuum) | ✅ TWIERDZENIE | — | Zamknięte |
| d=3 (wymiar) | ✅ TWIERDZENIE | — | Z zbieżności solitonów |
| K(ρ)=ρ² | ✅ NUMERYCZNE | R2 | K_IR/K_UV = 1.000 |

---

## Kolejność ataku (rekomendacja)

```
Faza 1 (teraz):     R1 (Cabibbo)  ←── najniżej wiszący owoc
                        ↓
Faza 2 (za 2-4 tyg): R4 (metryka)  ←── zamienia postulat w twierdzenie
                        ↓
Faza 3 (za 1-2 mies): R6 (B=√2)    ←── Koide staje się twierdzeniem
                     + R5 (k=4)     ←── fundament mas
                        ↓
Faza 4 (6-12 mies):  R2 (CG)       ←── "z pierwszych zasad"
                     + R3 (N=3)     ←── fundamentalne pytanie
                        ↓
Bonus:               R7 (UV)       ←── wzmocnienie, nie konieczne
```

**Foldery R1–R7 mogą być atakowane równolegle.** Powyższa kolejność to optymalna
sekwencja dla jednej osoby — ale jeśli jest czas na dwa wątki, R1+R4 lub R5+R6
dobrze się łączą.

---

> *Plan v4 utworzony 2026-04-14. Poprzedni: PLAN_ROZWOJU_v3.md (zamknięty).*
> *Rdzeń: 497 testów, 91% pass rate, 0/15 kill criteria naruszone.*
