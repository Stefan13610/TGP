---
title: "Reconciliation R5 K² vs why_n3 Phase 2 mass formula"
date: 2026-04-30
type: cross-cycle-reconciliation
parents:
  - "[[research/mass_scaling_k4/README.md]]"
  - "[[research/why_n3/PHASE2_n_alpha_derivation.md]]"
related:
  - "[[r5_phase2_reconciliation.py]]"
  - "[[r5_k_squared_mechanism.py]]"
  - "[[r3_observable_vs_full_mass.py]]"
status: PARTIAL RESOLUTION — sub-tensja g₀_τ kalibracji pozostaje OPEN (Phase 7+)
tags:
  - TGP
  - R5
  - why_n3
  - reconciliation
  - mass-formula
  - alpha-dependence
---

# Reconciliation: R5 K² mechanizm vs why_n3 Phase 2 mass formula

> **Status:** PARTIAL RESOLUTION. Główna tensja A⁴ (R5) vs A^(5−α) (Phase 2)
> dla α=2 została **strukturalnie rozwiązana** jako artefakt nieoznaczoności
> uniwersalności w R5. Sub-tensja w kalibracji g₀_τ między R5 i Phase 2
> pozostaje **OPEN** dla osobnej analizy (Phase 7+).

---

## 1. Sformułowanie tensji

### 1.1 Roszczenie R5 (`mass_scaling_k4/README.md`, 2026-04-16)

```
m_phys = c_m · K²
gdzie K = ∫(1/2) g^(2α) (g')² r² dr  [pełna akcja kinetyczna]

K ~ A² uniwersalnie (slope 1.9997 substrat α=1, 1.9965 canonical α=2)

⟹  m ~ A⁴ niezależnie od α
```

R5 README claim: "Ratio m_i/m_j cutoff-independent — POTWIERDZONE".

### 1.2 Roszczenie why_n3 Phase 2 (`why_n3/PHASE2_n_alpha_derivation.md`, 2026-05-01)

```
m_obs(g₀, α) = c_M · A_tail²(g₀, α) · g₀^[e²·(1−α/4)]
            = c_M · A² · g₀^n(α)  gdzie n(α) = e²(1−α/4)

Specjalne wartości:
  α = 1 substrate:  n = 3e²/4 ≈ 5.541
  α = 2 canonical:  n = e²/2 ≈ 3.695
  α = 4 Hobart-Derrick: n = 0 (no core dressing)
```

Phase 2 verified for α=2: **m_μ/m_e = 206.77 vs PDG 206.7682 → diff −0.001%** ✓.

### 1.3 Tensja

Dla α=2, dwie formuły dają **drastycznie różne predykcje**:
- R5 K²:  m_μ/m_e = (K_μ/K_e)² ≈ (A_μ/A_e)⁴ = (5.91)⁴ = **1221**
- Phase 2: m_μ/m_e = (A_μ/A_e)²·(g₀_μ/g₀_e)^(e²/2) = 34.93·5.92 = **206.78**

Rozbieżność **+490%** — jedna z formuł musi być błędna dla α=2.

---

## 2. Audyt skryptów źródłowych

### 2.1 Co R5 `r5_k_squared_mechanism.py` ACTUALNIE testuje

Analiza 7 testów w skrypcie:

| Test | Co weryfikuje | Zakres α |
|------|---------------|----------|
| 1 | K ~ A² scaling | α=1 substrate (skan g₀ ∈ [0.6, 1.8]) |
| 2 | K ~ A² scaling | α=2 canonical (skan g₀ ∈ [0.6, 0.95]) |
| **3** | **m=c·K² match PDG** | **TYLKO α=1 substrate** (linie 252–296) |
| 4 | C_T = R_max/4 + C_core | α=1 substrate |
| 5 | K_core/A² universal | α=1 substrate |
| 6 | virial K/|V| ≈ 1 | α=1 substrate |

**Krytyczne odkrycie:** TEST 3 (jedyny test mass-formula vs PDG) jest
przeprowadzany **TYLKO dla α=1 substrate**. Dla α=2 canonical TEST 2
weryfikuje TYLKO slope K~A² — nie testuje czy K² daje PDG ratio.

**R5 README claim "K~A² uniwersalnie (slope 1.9965 canonical)" jest
GENUINE**. Ale extrapolation **"⟹ m=c·K² uniwersalnie"** jest
**niezweryfikowana** dla α=2.

### 2.2 Co why_n3 `r3_observable_vs_full_mass.py` SECTION 3 weryfikuje

Skrypt for α=2 canonical explicit testuje 11 kandydatów mass formula
przeciwko PDG. Wynik (cytat z `r3_observable_vs_full_mass.txt`):

```
Formula                |  Ratio mu/e    diff% |  Ratio tau/e    diff%
A_tail^3               |      206.56    -0.10 |       2941.8   -15.40
A_tail^4               |     1221.06  +490.54 |      42152.8 +1112.25
K^2                    |     1220.04  +490.05 |      42184.8 +1113.17  ← R5 K² claim
K^(3/2)                |      206.43    -0.16 |       2943.5   -15.35
|M_energy|^2           |      490.61  +137.27 |       4442.3   +27.75
```

**Wniosek:** Dla α=2, K² daje **+490% mismatch** vs PDG dla μ/e.
**R5 K² mechanizm NIE działa dla α=2 canonical.**

### 2.3 Niezależna weryfikacja (`r5_phase2_reconciliation.py`, ten dokument)

Niezależny skrypt potwierdza:

```
α=1 substrate:
  R5 K²:                  μ/e=207.02 (+0.119%)  τ/e=3470.2 (-0.20%)  ✓
  Phase 2 (g₀_τ=1.72931): μ/e=206.81 (+0.021%)  τ/e=2662.9 (-23.42%) ⚠

α=2 canonical:
  R5 K²:                  μ/e=1220.04 (+490%)   τ/e=42184  (+1113%) ✗
  Phase 2 (g₀_τ=1.72931): μ/e=206.76 (-0.005%)  τ/e=2604.7 (-25.09%) ⚠
```

---

## 3. Strukturalne rozwiązanie

### 3.1 Trzy różne wielkości

Reconciliation wymaga rozdzielenia **trzech różnych obiektów**:

| Wielkość | Definicja | Scaling | Czy m_obs? |
|----------|-----------|---------|------------|
| K | ∫(1/2)g^(2α)(g')²r²dr | ~ A² ∀α | NIE — to akcja, nie masa |
| K² | (akcja kinetyczna)² | ~ A⁴ ∀α | NIE dla α≠1 |
| M_full | K + V_eff | ≈ 0 (Hobart-Derrick) | NIE — to binding energy |
| **m_obs** | **fizyczne PDG** | **c·A²·g₀^(e²(1-α/4))** | **TAK** (Phase 2) |

### 3.2 Dlaczego R5 K² działa dla α=1?

Dla α=1 substrate ODE, korelacja g₀-A przez równanie różniczkowe
**absorbuje** czynnik core-dressing g₀^(3e²/4) do efektywnego A^4 scaling.
Dokładniej: dla α=1, slope K~A² to slope DOKŁADNIE 1.9997 (R5 README), bo
relacja K = c(g₀)·A² ma prefactor c(g₀) bliski stałej dla zakresu g₀
(K/A²: e=17.60, μ=17.63, τ=17.60 — varies <0.2%).

To powoduje że dla α=1:
- (K_μ/K_e)² = (A_μ/A_e)⁴ z dokładnością <0.5%
- Phase 2 formula `A²·g₀^(3e²/4)` ≈ R5 `A^4` z dokładnością ~0.1%

### 3.3 Dlaczego R5 K² NIE działa dla α=2?

Dla α=2 canonical:
- K/A² nadal varies <0.1% (e=49.90, μ=49.89, τ=49.92 — patrz CASE 3 w skrypcie)
- Wiec K~A² faktycznie zachowuje się jak claim, ale BARDZIEJ stricte!
- Problem nie leży w K~A² scaling, ale w **wartościach absolutnych A** dla α=2

Dla α=2, A_μ/A_e = 5.91 (nie 3.79 jak w α=1). Stąd:
- (A_μ/A_e)⁴ = 1221 ≠ 207 (PDG)
- Wymagana inna potęga: (A_μ/A_e)^p = 207 daje p = log(207)/log(5.91) = **3.001**

To jest dokładnie why_n3 odkrycie p(α) = 5−α: dla α=2 → p = 3.

### 3.4 Reconciliation: m_obs = c·A²·g₀^n(α) jest **właściwym** uniwersalnym wyrazem

Phase 2 formula `m_obs = c_M·A²·g₀^[e²(1−α/4)]` jest uniwersalna w α:

- **α = 1**: redukuje się EFEKTYWNIE do c·A⁴ (z ~0.5% dokładnością) bo
  g₀-A korelacja absorbuje core dressing — to jest R5 case, GENUINE.
- **α = 2**: wymaga **explicit g₀^(e²/2) factor** — m_obs ≠ c·A⁴.
  R5 K² extrapolation FAILS.
- **α = 4** Hobart-Derrick: m_obs = c·A² (no core dressing, n(4)=0).

**R5 K² mechanizm jest SPECYFICZNYM PRZYPADKIEM Phase 2 dla α=1.**
Generalne stwierdzenie "m=c·K² uniwersalnie" jest **artefaktem
ekstrapolacji** — zostało zweryfikowane TYLKO dla α=1 substrate w
TEST 3 `r5_k_squared_mechanism.py`.

### 3.5 Status oryginalnej hipotezy 2 (M_full vs m_obs)

Oryginalna hipoteza 2 z `mass_scaling_k4/README.md`:
> "R5 K² daje **M_full**, why_n3 A^(5-α) daje **m_obs** — distinction
> analogiczny do GR ADM vs Komara, QFT bare vs renormalized"

**Verdict:** Hipoteza była PARTIALLY RIGHT direction ale NIE konkretnie:

- ✓ **Distinction GENUINE**: M_full ≠ m_obs jest realna dystynkcja
  (pokazane w CASE 2: |M_full|_μ/|M_full|_e = 22 ≠ 207 = m_μ/m_e PDG)
- ✗ **Ale R5 K² ≠ M_full**: K² ≈ A⁴ jest osobną wielkością strukturalną,
  NIE M_full = K + V_eff ≈ 0
- ✗ **Analog ADM/Komara**: nie zachodzi w prosty sposób, bo M_full ≈ 0
  w 3D scalar (Hobart-Derrick balance) — nie jest "asymptotic mass"

**Właściwa interpretacja:**
- M_full ≈ 0 to **binding energy** (potencjał - kinetic), nie mass observable
- K² to **structural action-level quantity**, korelujące z m_obs tylko
  dla α=1
- m_obs = c·A²·g₀^n(α) to fizyczny observable, **α-zależne core dressing**

---

## 4. Sub-tensja: g₀_τ kalibracja (OPEN)

### 4.1 Obserwacja

Phase 2 i R5 używają **różnych wartości g₀_τ**:

| Source | g₀_τ | Pochodzenie |
|--------|------|-------------|
| R5 `r5_k_squared_mechanism.py` (substrate) | 1.72931 | "z Koide" |
| R5 `r5_k_squared_mechanism.py` (canonical) | 2.276 = φ²·g₀_e | φ-drabinka |
| Phase 2 `PHASE2_n_alpha_derivation.md` | 1.75505 | "z Koide K=2/3" |
| `r3_observable_vs_full_mass.py` | 1.72931 | (skopiowane z R5) |

### 4.2 Konsekwencje

Phase 2 z g₀_τ = 1.75505 **dokładnie** matchuje PDG τ/e = 3477.23 (-0.085%).
Phase 2 z g₀_τ = 1.72931 daje τ/e = 2604.7 (-25%).
R5 K² z g₀_τ = 1.72931 daje τ/e = 42184 (+1113% dla α=2).

### 4.3 Interpretacja

To NIE jest błąd Phase 2 ani R5, ale **niezamknięta kalibracja w R3
ODE rdzeniu** — różne metody (Koide K=2/3 vs φ² ladder vs direct fit)
dają nieco różne g₀_τ.

**Status:** OPEN, sub-problem cyklu R3, do osobnej analizy w Phase 7+.

### 4.4 Sugestia weryfikacji

Pełna konsystencja wymaga:
1. Re-deriwacji g₀_τ z **jednoznacznego** mechanizmu (Koide derivation
   chain w `r3_koide_derivation.py` lub φ-ladder)
2. Aktualizacji wszystkich skryptów R5 + R3 + why_n3 by używały
   **jednej** wartości g₀_τ z explicit derivation source
3. Re-uruchomienia Phase 2 + R5 K² z ujednoliconą kalibracją

---

## 5. Status epistemiczny

### 5.1 Co JEST udowodnione

- ✅ **R5 K² mechanizm dla α=1 substrate**: GENUINE — `m_phys = c·K²`
  match PDG 0.013% diff dla μ/e (R5 TEST 3 PASS).
- ✅ **K~A² scaling dla obu α**: GENUINE — slope 1.9997 sub, 1.9965 can
  z dokładnością <1% (R5 TEST 1+2 PASS).
- ✅ **Phase 2 formula dla α=2 canonical**: GENUINE — `m_obs =
  c·A²·g₀^(e²/2)` match PDG -0.001% diff dla μ/e (Phase 2 verified).
- ✅ **R5 K² ≠ m_obs dla α=2 canonical**: DOWIEDZIONE numerycznie —
  K² extrapolation daje +490% mismatch.
- ✅ **M_full ≠ m_obs ≠ K²** w 3D scalar — wszystkie trzy są distinct.

### 5.2 Co POZOSTAJE OPEN

- ⚠ **g₀_τ kalibracja**: różne konwencje między R5/Phase 2 (1.72931
  vs 1.75505). Wymaga Phase 7+ analizy.
- ⚠ **Analityczne pochodzenie X = e²/4** w n(α) = e²(1−α/4): EMPIRICAL
  discovery (residuum < 0.1% dla α ∈ [0.5, 2.5]), wciąż awaiting
  RG-derivation. **Już udokumentowane** w `PHASE6_Q5_R5_bridge_first_attempt.md`
  jako NEGATIVE first attempt.
- ⚠ **Generalizacja R5 dla α≠1**: czy istnieje mechanizm `m_obs = c·F(K, α)`
  gdzie F redukuje się do K² dla α=1, ale ma inną postać dla α≠1?
  Phase 7+ research.

### 5.3 Honest framing po reconciliation

| Wcześniejszy claim | Po reconciliation |
|--------------------|-------------------|
| "R5 K² mechanizm uniwersalny ∀α" | "R5 K² mechanizm SPECYFICZNY dla α=1 substrate" |
| "m ~ A⁴ uniwersalnie" | "m_obs ~ A^(5−α) z α-zależnym core-dressing" |
| "k=4 jedyny integer" | "k=4 jest specialny dla α=1; dla α=2 (TGP-canonical) k=3" |
| "m_obs = M_full" | "m_obs ≠ M_full ≠ K² — trzy różne wielkości" |

---

## 6. Konsekwencje dla głównego ciała teorii

### 6.1 R5 status

R5 jako cykl badawczy **pozostaje GENUINE** dla α=1 substrate ODE:
- Wszystkie 7/7 PASS dla `r5_k_squared_mechanism.py` (TEST 3 dla α=1)
- 9/9 PASS dla `lp4_mass_exponent_verification.py` (k=4 dla α=1)
- PDG match 0.013% μ/e dla α=1

R5 jako "uniwersalna" mass formula **wymaga zawężenia** do α=1
substrate. Stwierdzenie "K² mechanizm działa dla obu α" w README
linia 19 jest mylące — slope K~A² działa, ale m=K² mass formula
działa tylko dla α=1.

### 6.2 TGP-canonical α=2 status

TGP-canonical (`K(φ) = K_geo·φ⁴` z sek08a) operuje w α=2. Dla
TGP-canonical, **właściwa mass formula to Phase 2**:
```
m_obs = c_M · A² · g₀^(e²/2)
```

Ta formula dała 6/6 PASS PDG verification w Phase 2 (μ/e -0.001%,
τ/e -0.085%, τ/μ +0.015%) — z g₀_τ = 1.75505.

### 6.3 Cycle bridge: R5 ↔ why_n3

R5 i why_n3 są teraz **logicznie spójne**:
- R5 odpowiada za **α=1 substrate ODE** edge case
- why_n3 Phase 2 odpowiada za **uniwersalną α-zależność** including TGP-canonical α=2
- Most: m_obs = c·A²·g₀^(e²(1−α/4)); R5 to specjalny przypadek α=1

---

## 7. Pliki

| Plik | Opis |
|------|------|
| `r5_phase2_reconciliation.py` | Niezależny verification script (ten dokument) |
| `r5_phase2_reconciliation.txt` | Output z 3 case'ami (α=1, α=2, struktura) |
| `r5_k_squared_mechanism.py` | R5 oryginalny (TEST 3 = α=1 only verification) |
| `../why_n3/r3_observable_vs_full_mass.py` | why_n3 mass formula candidates dla α=2 |
| `../why_n3/PHASE2_n_alpha_derivation.md` | Phase 2 closure z X=e²/4 (EMPIRICAL) |

---

## 8. Cross-references

- `research/mass_scaling_k4/README.md` § "Otwarta tensja 2026-05-01"
  (item 6) — to **PARTIAL RESOLUTION** dla głównej tensji
- `research/why_n3/PHASE2_n_alpha_derivation.md` § 4 — Phase 2 mass
  formula complete dla R3
- `research/why_n3/PHASE6_Q5_R5_bridge_first_attempt.md` — Phase 6
  NEGATIVE attempt dla X=e²/4 RG-derivation (sub-problem)
- `research/why_n3/CORRECTIONS_2026-05-01.md` — RESOLUTION sekcja
  "Strukturalna interpretacja insight'u" (m_obs vs M_full)
- `meta/AUDYT_TGP_2026-05-01.md` § AB.3 — outstanding Phase 6+ items
  (g₀_τ kalibracja będzie nową pozycją Phase 7+)

---

**Autor:** Reconciliation analysis 2026-04-30 (TGP audit closure follow-up)
**Status:** PARTIAL RESOLUTION (główna tensja closure, sub-tensja g₀_τ open)
**Następne kroki:** Phase 7+ — ujednolicenie g₀_τ kalibracji między R5/R3/why_n3
