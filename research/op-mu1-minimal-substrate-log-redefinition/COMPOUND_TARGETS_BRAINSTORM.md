---
title: "μ.1 follow-up — Brainstorm targets dla compound emergence form exp(Σε)"
date: 2026-05-02
type: brainstorm
parent: "[[README.md]]"
trigger: po μ.1 NO-GO closure user spytał gdzie compound exp(Σε) mogłoby być derivable
status: BRAINSTORM — analiza portfolio TGP
tags:
  - TGP
  - mu1
  - brainstorm
  - compound-search
  - portfolio-review
---

# Brainstorm — gdzie w TGP `exp(Σε)` z derivable Σε?

## TL;DR

**Compound exp(Σε) NIE jest natywną matematyką TGP-programu.** Po przeszukaniu
portfolio (12 cykli przegląnięte): TGP rozwiązuje swoje problemy przez
**algebraiczne counting** (Cabibbo: GL(3,𝔽₂) factor 165/167), **geometric
ratios** (atomic shells: 2n²), **soliton tail amplitudes** (mass formula),
**RG flow** (running α — testowany w M.4, NEG). **Żaden** główny TGP-success
nie używa `exp(Σε)` natively.

**Konsekwencja:** szukanie compound-target jest **strukturalnie naciągane**.
Najsilniejsza możliwość — PMNS/Majorana topology (op-zeta winding + double-cover) —
wymagałaby istotnego deep-dive bez gwarantowanego payoff.

**Rekomendacja:** zrezygnować z compound-search, wybrać inny target z
TGP-portfolio (G1 lub open cycle).

---

## Macierz kandydatów (post-search)

| # | Target | TGP-cykl | Form natywna | Σε derivable? | Score |
|---|--------|----------|--------------|---------------|-------|
| 1 | Cosmologic. e-foldings N_e≈60 | cosmo_tensions, hubble_tension, op-cosmology-closure | log a(t) z FRW (nie exp(Σε)) | NIE — counting H-crossings retrospekcyjne | **2/10** |
| 2 | Running α(μ) RG flow | op-alpha-fine-structure, op-omega1-substrate-em-coupling | exp[−β₀ ln μ] | nie z topologii (β empirical) | **3/10** |
| 3 | Cabibbo/Wolfenstein CKM | cabibbo_correction, op-eta-wolfenstein | macierz produktów rotacji | TAK — GL(3,𝔽₂) algebra | **7.5/10** ✓ |
| 4 | PMNS neutrino mixing | op-iota, op-zeta, op-mu, neutrino_msw | macierz produktów + Majorana | TAK potencjalnie — winding + 2π double-cover | **6-7/10** ◐ |
| 5 | Atomic shells 2n² | atom_from_soliton, atomic_shells_closure | power law (nie exp) | log compound, ale power-law base | 4/10 |
| 6 | g-2 anomaly | muon_g_minus_2 | seria α/π korekt | empirical perturbacyjnie | 4/10 |
| 7 | Yukawa hierarchy | (rozproszone) | m_t/m_e ~ exp(11.5) | NIE — to hierarchy problem | 3/10 |
| 8 | BH entropy S=A/4 | op-bh-alpha-threshold, op-eht | exp(microstates) | trudne dla makro-BH | 3/10 |
| 9 | Casimir vacuum | casimir_mof | ζ-regularizacja | nie compound exp | 3/10 |
| 10 | Neutrino mass hierarchy | neutrino_msw, op-omicron1 | m_3/m_1 ~ 100 = exp(4.6)? | uncertain | 4/10 |

---

## Top 2 kandydaci (najsilniejsze)

### **#3 Cabibbo / CKM** (Score 7.5/10) — ALE JUŻ CLOSED

**Stan obecny:**
- `cabibbo_correction/` — 4.8σ → **0.75σ** po Z₃ self-energy correction
- Wolfenstein phase 1: 5/5 PASS
- **Compound structure:** CKM = ∏ rotacji R₁·R₂·R₃ z phase factors

**GL(3,𝔽₂) Σε-like derivation:**
- |GL(3,𝔽₂)| = 168
- |Z₃ subgroups| = 28
- N_eff = 3.036 z (|G|−|Z₃|)/2 = 165/2 ≈ 82.5 lub 165/167 = 0.988
- To jest **algebraiczny argument**, nie compound exp

**Kompleks:** TGP-program już rozwiązał Cabibbo **bez** compound exp(Σε).
Struktura okazała się być algebraiczna, nie compound. **Reopening dla compound
re-interpretation — niska wartość**, bo działający mechanizm jest prostszy.

### **#4 PMNS neutrino + Majorana** (Score 6-7/10) — ACTIVE, NIEZBADANE

**Stan obecny:**
- `op-iota-charge-pmns-unification/` — phase 1-3 active
- `op-zeta/` — phase 2 z `phase2_pmns_derivation.txt`
- `op-mu-pmns-phase-hardening/` — phase z δ_CP
- `op-nu-majorana-phase-mbb/` — Majorana m_ββ

**Potencjalny compound source:**
- Majorana fermions: ψ → ψ + 2π pod rotation, **double-cover** SU(2)→SO(3)
- Σε = 2 mogłoby pochodzić z **2 Majorana phases** (per particle 4π, normalized to 2π = 2 units)
- **Winding argument** w op-zeta phase2 — wymaga deep-dive

**Open question:** czy Σ_winding charges = 2 dla Majorana neutrino dają compound
emergence form pod amplitude/phase split (μ.1 candidate d, score 5/10)?

**Realność:** umiarkowana. Wymaga 1-2 dni deep-dive op-zeta i op-mu phases bez
gwarantowanego success. Mogłoby być POWERFUL jeśli winding rzeczywiście wymusza
N_e = 2 dla Majorana.

---

## Wyklutowane (silne NEGATIVE)

### Cosmologic e-foldings (Score 2/10)
- TGP już TESTOWAŁ ten path: cosmo_tensions ct7 — "definityvny werdykt TGP poza
  zasięgiem"
- H-backreaction mechanizmy (tachyonic, soliton population, RG running) —
  **wszystkie 8-9 rzędów wielkości za małe**
- N_e nie ma natywnie formy exp(Σε); pochodzi z całki ∫H dt/dz

### Running α(μ) RG flow (Score 3/10)
- M.4 z λ.1 testował exactly to: γ_φ ≈ 0.03 vs target 3.69 → 125× miss
- Compound form RG jest naturalna ale Σε = β-coefficient jest empirical
- **Recurrent path** — nie warto re-testować w innym zakresie

---

## Verdict ostateczny dla compound search

**μ.1 + portfolio review → compound emergence form exp(Σε) z derivable Σε
NIE JEST natywna w TGP.**

**Main reason:**
- TGP success cases (Cabibbo, atomic shells, Wolfenstein, mass scaling) używają
  innych struktur: algebraic counting, geometric ratios, soliton geometries
- Tylko 1 kandydat (PMNS/Majorana) ma realistyczny szansa, ale wymaga zaangażowania
  bez jasnego payoff

**Strategiczna rekomendacja:**
1. **Zamknij compound-search** — dwa cykle (M.6, μ.1) wyczerpały podstawowe
   ścieżki; portfolio review nie pokazuje silniejszych celów
2. **Wróć do G1** (Φ_eff anchor inconsistency resolution) — pierwotny next-target
   po λ.1, niezablokowany przez compound NIL
3. **Alternatywa:** jeśli chcesz jeszcze jedną próbę compound, **op-zeta
   Majorana deep-dive** ma najsensowniejszą szansę (~50%) ale ROI niepewny

---

## Bonus: structural lesson z μ.1 + brainstorm

**TGP nie jest compound-program.** TGP-substrate produkuje:
- **Power-law tails** w solitonach (mass formula `m = c·A²·g₀^X`)
- **Algebraiczne counting** w mixing (Cabibbo GL(3,𝔽₂))
- **Geometric ratios** w shells (2n²)
- **Soliton-soliton interactions** w sub-tensji τ

`exp(Σε)` form jest **natural w QFT** (RG, partition function, path integral),
ale w TGP-substrate **nie znalazłem natywnego mechanizmu**. To może być
**fundamentalna cecha** TGP: jest to non-perturbative geometric framework,
a compound exp(Σε) jest perturbative mathematics.

**To jest valuable null result:** TGP-program nie potrzebuje compound exp(Σε)
dla swoich successes. λ.1 NEG + μ.1 NO-GO + portfolio NIL **łącznie pokazują**
że X = e²/2 mass formula jest bardziej prawdopodobnie:
- Numerical coincidence (0.0007% match z exp(2)/2 może być przypadek)
- Albo deep-structural fakt którego nasze obecne narzędzia nie ujawniają
- ALE NIE jest naturalnie compound-derivable z TGP-substrate

---

## Files

- `README.md` — μ.1 main verdict (NO-GO closure)
- `PLAN.md` — pre-implementation plan
- `phase1_psi_*.py` — reparametryzacja proof
- `phase2_compound_*.py` — compound + topology test
- `COMPOUND_TARGETS_BRAINSTORM.md` (this file) — portfolio review po μ.1

## Następne (proposed)

Wybór:
1. **G1 cycle** (Φ_eff anchor) — pierwotny next-target, fresh
2. **op-zeta deep-dive** — 50% szansa znalezienia Majorana Σ-source, niska gwarancja
3. **Inny open cycle z portfolio** — np. op-omega1 substrate-EM coupling lock,
   op-rho1 71Ge cross-section, op-omicron1 Σm_ν cosmo
