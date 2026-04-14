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
| **R3** | Dlaczego N=3 generacji | `research/why_n3/` | nieznany | ⭐⭐⭐⭐⭐ |
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
- [ ] Scalenie z rdzeniem (tgp_companion.tex §F2)
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

## R3: Dlaczego N=3 generacji — `research/why_n3/`

**Problem:** NAJFUNDAMENTALNIEJSZE otwarte pytanie. GL(3,𝔽₂) zakłada N=3, nie wyprowadza.

**Obecne heurystyki (żadna nie jest dowodem):**
- Bariera duchowa: d=3 → k=4 → WKB daje dokładnie 3 stany związane
- |GL(3,𝔽₂)| = 168 = (2N+1)·2^N·N! — ale to konsekwencja N=3, nie przyczyna
- Anomaly cancellation wymaga N_gen = N_color = 3 — powiązanie, nie dowód

**Ścieżki ataku:**
1. **Topologiczna**: Czy π₁(konfiguracji solitonów w d=3) wymusza 3 klasy?
2. **Dynamiczna**: Czy 4. generacja jest dynamicznie zakazana? (g₀^(4) > g₀_crit = 8/5)
   - Już sprawdzone numerycznie (H8: PASS) — ale to nie odpowiada "dlaczego"
3. **Informacyjna**: Czy entropia Shanona na substracie Z₂ w d=3 wymusza 3 sektory?
4. **Algebraiczna**: Czy GL(N,𝔽₂) dla N≠3 prowadzi do sprzeczności fizycznych?
   - GL(2,𝔽₂) = S₃: 6 elementów → za mało struktur
   - GL(4,𝔽₂): 20160 elementów → prowadzi do 4 generacji (zakazanych eksperymentalnie)

**Kryterium zamknięcia:** Wyprowadzenie N=3 z fizycznego argumentu (nie tautologia)

**Status:** Może wymagać nowej fizyki. Ryzyko: NIEZNANE.

---

## R4: Ansatz metryczny h(Φ)=Φ — `research/metric_ansatz/`

**Problem:** Metryka ds² = -(c₀²/ψ)dt² + ψδᵢⱼdxⁱdxʲ z ψ=Φ/Φ₀ jest **postulatem**.
Dlaczego p=1, a nie Φ^p?

**Obecne uzasadnienie:**
- Numerycznie: 2 kryteria wybierają p=1 (fiber bundle + r₂₁)
- Ale: Shapiro delay i PPN γ=1 są zdegenerowane (dowolne p daje γ=1)

**Plan ataku (trzy niezależne ścieżki):**

### A2a: Fonony na substracie
- Obliczyć relację dyspersji ω(k) fononów na Γ
- Sprawdzić czy c_sound ∝ √Φ wymusza liniowy coupling h=Φ
- Nakład: 2 tygodnie

### A2b: Test równań Einsteina (NAJBARDZIEJ OBIECUJĄCY)
- Wstawić g_ij = Φ^p·δ_ij do równań pola Einsteina
- Sprawdzić dla jakich p zachodzi:
  - (a) Brak duchów (ghost-free)
  - (b) Pozytywna energia
  - (c) Poprawny limit newtonowski
- Hipoteza: Tylko p=1 spełnia wszystkie warunki jednocześnie
- Nakład: 2–3 tygodnie

### A2c: Argument informacyjny
- "Bity przestrzenne" ∝ Φ w d=3
- Entropia Bekenstein-Hawking wymaga liniowy coupling
- Nakład: 1–2 tygodnie

**Kryterium zamknięcia:** Twierdzenie: "Spośród h(Φ)=Φ^p, tylko p=1 daje spójną teorię"

**Pliki rdzenia do scalenia:**
- Nowy paragraf w `sek08c_metryka_z_substratu.tex`

---

## R5: Prawo skalowania m ∝ A_tail⁴ — `research/mass_scaling_k4/`

**Problem:** Formuła masowa jest **fundamentem** sektora leptonowego. Daje r₂₁ = 206.768
(0.0001% zgodność z PDG). Ale k=4 jest **postulatem**.

**Obecne uzasadnienie:**
- Argument wymiarowy: w d=3, zbieżność ogona wymaga n > 2d/(d-1) = 3, więc k=4
  (ex188_A4_dimensional_argument.py)
- Alternatywy odrzucone: E²=0 (OK), E³ (FAIL, -0.647), prop:K-exponent (OBALONY, ex150)

**Plan ataku:**
1. **Z twierdzenia wiriałowego:**
   - Soliton w polu z K(Φ)=Φ^α — wyprowadzić relację między energią a A_tail
   - Sprawdzić czy E_bind ∝ A⁴ wynika z α=2 i d=3
2. **Z funkcjonału działania:**
   - S[g] = ∫[K(g)·(∇g)² + P(g)] — czy masa = ∂S/∂(parametr)?
   - Jakie potęgi A_tail pojawiają się naturalnie?
3. **Z analizy asymptotycznej ogona:**
   - g(r) ~ 1 - A·e^{-μr}/r^ν dla r→∞
   - Masa jako residue w transformacie Mellina

**Kryterium zamknięcia:** Twierdzenie: "m ∝ A⁴ wynika z α=2, d=3, K(Φ)=Φ²"

**Pliki rdzenia do scalenia:**
- Rozszerzenie `dodatekJ_ogon_masy.tex`

---

## R6: B=√2 analitycznie — `research/brannen_sqrt2/`

**Problem:** Stosunek Brannena B = b/a = √2 jest potwierdzony numerycznie do 10⁻⁶,
ale brak dowodu analitycznego. Gdyby udowodnić → Koide K=2/3 staje się **twierdzeniem**.

**Obecny status:**
- a3d_soliton_brannen_r.py: 5/6 PASS, T5 = EXPECTED FAIL
- B_num = 1.414212... ≈ √2 (|δ| < 10⁻⁶)
- Związek: K = (1 + B²/2)/N = (1 + 1)/3 = 2/3

**Plan ataku:**
1. **Analiza ODE solitonowego:**
   - g'' + (2/r)g' = dP/dg z K(g)=g^{2α}
   - Warunki brzegowe: g(0)=g₀, g'(0)=0, g(∞)=1
   - Szukać relacji między g₀ a parametrami Brannena (a, b, θ)
2. **Symetria ODE:**
   - Czy ODE ma ukrytą symetrię SO(2) lub dyskretną S₃ wymuszającą B=√2?
   - Analiza Lie punktowych symetrii ODE
3. **Perturbacyjnie wokół B=√2:**
   - Założyć B = √2 + ε, sprawdzić czy ε=0 jest punktem stacjonarnym
   - Minimalizacja jakiego funkcjonału daje B=√2?

**Kryterium zamknięcia:** Dowód analityczny B=√2 z ODE solitonowego

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
| h(Φ)=Φ | 🟡 POSTULAT | R4 | 2 kryteria numeryczne |
| m ∝ A⁴ | 🟡 POSTULAT + heurystyka | R5 | Argument wymiarowy |
| N=3 | 🟡 HEURYSTYKA | R3 | Bariera duchowa |
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
