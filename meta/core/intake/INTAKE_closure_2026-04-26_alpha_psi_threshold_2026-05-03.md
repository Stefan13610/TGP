---
title: "INTAKE — closure_2026-04-26/alpha_psi_threshold → core (alternative track for OP-M92 α-universality)"
date: 2026-05-03
type: intake
status: PENDING-HUMAN-REVIEW
source_folder: research/closure_2026-04-26/alpha_psi_threshold
target_section: core/sek06_czarne_dziury (option: dodać α(ψ) threshold subsection) + WEP MICROSCOPE update
hotspot_addressed: H-B9 (alternative track; main closure: § W)
agent_requesting: Sesja 3.5 Q6 (agent)
parent: "[[meta/core/CORE_INTAKE.md]]"
related:
  - "[[research/closure_2026-04-26/alpha_psi_threshold/results.md]]"
  - "[[research/closure_2026-04-26/alpha_psi_threshold/setup.md]]"
  - "[[research/op-m92/OP_M92_P0plus_candD_multisource_results.md]]"
  - "[[meta/core/CORE_HOTSPOTS.md]]"
tags:
  - intake
  - core-promotion
  - cross-validation
  - alpha-psi
  - WEP-MICROSCOPE
  - T-alpha
---

# INTAKE — `closure_2026-04-26/alpha_psi_threshold` → core (alternative track for OP-M92)

## 1. Wynik proponowany do promocji

**Path E (T-α): `α(ψ) z ψ-threshold`** strukturalnie rozwiązuje
multi-source α-universality issue z OP-M92 Phase 0+:

```
α(ψ) = α₀ × (ψ - 1)² × Θ(ψ - 1)         [α₀ ≈ 4 dimensionless O(1)]
```

**Mechanizm:**
- `ψ_ph = 1.168` jest UNIVERSAL przy photon ring dla wszystkich
  Schwarzschild BH (geom-units) → `α(ψ_ph) = 0.0282 × α₀ ≈ 0.114`
  identyczna dla SgrA*, M87*, GW150914, NS
- `α(ψ_Earth)` tłumione faktorem `(G·M_Earth/c²R)² ≈ 5×10⁻¹⁹`
  → WEP MICROSCOPE margin **wzrasta 6.7× → 4×10¹⁶×**

**Konsekwencja dla M9.2-D:** lead candidate status RESTORED.

## 2. Source w research/

- `research/closure_2026-04-26/alpha_psi_threshold/results.md` (10.9 KB) —
  TL;DR + 5/5 PASS structural tests
- `research/closure_2026-04-26/alpha_psi_threshold/alpha_psi_threshold.py`
  (13.7 KB) + `.txt` (6.4 KB) — sympy + numerical multi-source verification
- `research/closure_2026-04-26/alpha_psi_threshold/setup.md` (7.7 KB) —
  audit design (Path E motywacja + 5 strukturalnych tests)

## 3. Cross-validation z audytem 2026-05-01

**KLUCZOWY KONTEKST:** Audyt 2026-05-01 § W **niezależnie** zaadresował
M9.2 WEP MICROSCOPE composition test (H-B9):

- `research/op-newton-momentum/B9_*.py` (6/6 PASS)
- M9.2 composition test sympy + numerical
- Path: composition test bez α-threshold

**Dwie różne strategie**:

| Ścieżka | Folder | PASS | Filozofia |
|---------|--------|-----:|-----------|
| Audyt main (§ W) | `op-newton-momentum/B9_wep_microscope_composition.py` | 6/6 | Composition test M9.2 (constant α implicitly) |
| Cross-validation | `closure_2026-04-26/alpha_psi_threshold/alpha_psi_threshold.py` | 5/5 | α(ψ) threshold zastępuje constant α; wzmacnia margin 6.7× → 4×10¹⁶× |

**To jest różnica jakościowa, nie ilościowa.** § W zamyka M9.2 WEP test
**przy zachowaniu** constant α; T-α **zastępuje** constant α funkcją
α(ψ), eliminując jednocześnie OP-M92 multi-source issue oraz
zwiększając margines WEP o 16 rzędów wielkości.

## 4. Proponowana lokacja w core

| Target | Co dokładnie | Czy istnieje | Diff |
|--------|--------------|--------------|------|
| `core/sek06_czarne_dziury.tex` | Nowa subsection: "α(ψ) threshold and OP-M92 multi-source α-universality resolution" | NO | (nowa subsection ~30 linii + tabela multi-source) |
| `core/sek06_czarne_dziury.tex` przy α dyskusji | Update: zastąp "constant α" notation z "`α(ψ) = α₀(ψ-1)²Θ(ψ-1)`" + cross-reference do nowej subsection | YES | dopisać 1-2 linie + cross-reference |
| `core/dodatekU` lub WEP-relevant section | Update WEP MICROSCOPE margin z 6.7× na 4×10¹⁶× po Path E | TBD (lokalizacja WEP discussion) | manualne sprawdzenie + 1-line update |
| `PREDICTIONS_REGISTRY.md` | Update OP-M92 entries z `multi-source caveat` na `RESOLVED via α(ψ) threshold` | YES (markery `M92`) | strikethrough caveats; positive verdict |

## 5. Hotspot, który ten wynik adresuje

**H-B9** w [[meta/core/CORE_HOTSPOTS.md]] §B (status: `RESOLVED-AUDIT-FULL`
przez § W). Ten INTAKE proponuje **upgrade z constant α do α(ψ) threshold**.

Plus: T-α **rozwiązuje OP-M92 multi-source α-universality** problem,
który nie był osobnym hotspotem audytu, ale jest paralelnym strukturalnym
luki — Phase 0+ multi-source check (2026-04-25) wykrył że naive Candidate D
z constant α implikuje `α_SI ~ M_BH² (19 OOM span NS↔M87*)`. T-α to
zamyka.

## 6. Co MUSI być sprawdzone

- [ ] Wynik kompiluje się (sympy 5/5 PASS confirmed w `results.md`)
- [ ] Brak konfliktu z aktualnym `main.tex` (read-only)
- [ ] `core_compatibility` w YAML folderu źródłowego (Sesja 4)
- [ ] **CRITICAL:** zgodność T-α z § W M9.2 composition test —
      czy T-α α(ψ_Earth) = α₀×5×10⁻¹⁹ przesuwa M9.2 wyniki w sposób
      wymagający re-runu § W przez Form-IV `√(-g)` z α(ψ) threshold?
      (W niewerbalnej wersji `α(ψ_Earth) → 0` więc m_field jest
      effectively M9.2-canonical, ale formalnie należy zweryfikować)
- [ ] PREDICTIONS_REGISTRY caveats `M92`/`OP-M92` należy zlokalizować
      i przygotować patch

## 7. Co się STANIE po akceptacji

1. **Człowiek** edytuje `core/sek06_czarne_dziury.tex` dodając α(ψ) subsection
2. **Człowiek** updateuje WEP MICROSCOPE numerical w odpowiedniej sekcji
3. **Człowiek** updateuje `PREDICTIONS_REGISTRY.md` OP-M92 entries
4. **Człowiek** dopisuje wpis do [[meta/core/CORE_INVENTORY.md]] §1.7
5. **Człowiek** updateuje status H-B9 w `CORE_HOTSPOTS.md`
6. **Człowiek** zamyka INTAKE

## 8. Decyzja człowieka

- [ ] **ACCEPT** — Path E α(ψ) threshold promowany do core (data: ____)
- [ ] **DEFER** — czekamy na osobny audit M9.2 z α(ψ) (powód: ____)
- [ ] **REJECT** — § W constant-α composition test wystarcza (powód: ____)
- [ ] **NEEDS-MORE-EVIDENCE** — (lista: ____)

---

## 9. Notatki

### 9.1 PRO

1. **WEP margin 4×10¹⁶× to JAKOŚCIOWY upgrade** — od "marginal pass"
   (6.7×) do "future-proof against any plausible MICROSCOPE-2 sensitivity"
2. **Rozwiązuje OP-M92 multi-source α-universality issue** — single
   structural fix dla 19 OOM α_SI spread NS↔M87*
3. **`α_0 ≈ 4` jest natural O(1) constant** — nie fitting

### 9.2 CONTRA / risks

1. **T-α modyfikuje aksjomat α** — to jest semi-architectural change;
   wymaga sprawdzenia konsystencji z `core/sek08*` formalizmem
2. **`α(ψ_Earth) = α₀×5×10⁻¹⁹`** — jeśli M9.2 m_field zależy od α
   (a powinno), trzeba re-run § W z α(ψ) zamiast constant
3. **`(ψ-1)² Θ(ψ-1)` jest non-analytic at ψ=1** — non-smooth
   transition wymaga sprawdzenia czy nie psuje variational principle
   `S_TGP` w słabym polu (sek08a)

### 9.3 Rekomendacja agent

**ACCEPT z warunkiem (NEEDS-MORE-EVIDENCE):** Path E to silne
strukturalne rozwiązanie, ale wymaga **explicit verification** że
M9.2 § W wyniki przeżywają tranzycję `constant α → α(ψ)`. Sugerowane:

1. (a) Człowiek akceptuje INTAKE i zleca dedicated cycle "M9.2-α(ψ)
   re-run" przed pełną promocją do core, **albo**
2. (b) Człowiek akceptuje promocję, ale ze "Verification needed" footnote
   w nowej subsection sek06, otwarcie zaznaczając że pełna konsystencja
   z M9.2 § W wymaga osobnego sprawdzenia.

Decyzja należy do człowieka.
