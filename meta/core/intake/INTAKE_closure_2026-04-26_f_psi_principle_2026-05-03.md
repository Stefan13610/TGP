---
title: "INTAKE — closure_2026-04-26/f_psi_principle → core (alternative track for f(ψ) deeper principle)"
date: 2026-05-03
type: intake
status: PENDING-HUMAN-REVIEW
source_folder: research/closure_2026-04-26/f_psi_principle
target_section: core/sek08c_metryka_z_substratu (option: dodać T-FP principle subsection)
hotspot_addressed: H-B8 (alternative track; main closure: § V)
agent_requesting: Sesja 3.5 Q6 (agent)
parent: "[[meta/core/CORE_INTAKE.md]]"
related:
  - "[[research/closure_2026-04-26/f_psi_principle/results.md]]"
  - "[[research/closure_2026-04-26/f_psi_principle/setup.md]]"
  - "[[research/closure_2026-04-26/CLOSURE_2026-04-26_SUMMARY.md]]"
  - "[[meta/core/CORE_HOTSPOTS.md]]"
tags:
  - intake
  - core-promotion
  - cross-validation
  - f-psi
  - T-FP-principle
---

# INTAKE — `closure_2026-04-26/f_psi_principle` → core (alternative track for f(ψ) deeper principle)

## 1. Wynik proponowany do promocji

**Zasada T-FP (Substrate Polynomial-Degree Normalization):** `f(ψ)` jest
*jednoznacznie wybranym* dimensionless rationalem typu `V(Φ)/Φⁿ` z
**`n = deg V = 4`**, znormalizowanym do próżni:

```
f(ψ) = [V(Φ) / Φᵈ] / [V(Φ_eq) / Φ_eqᵈ]    z d = 4 unique
```

Cztery warunki fizyczne wymuszające `n = 4` UNIQUE:
1. **Asymptotic finite nonzero** przy ψ → ∞ (eliminuje n ≥ 5)
2. **Singular at 0** (eliminuje n ≤ 2)
3. **Zero-inheritance** `f(ψ_0)=0 ⟺ V(ψ_0)=0` (eliminuje n ∉ {n: V's zero})
4. **No spurious zeros** (eliminuje n=3 — extra zero)

**Trzy "konwergentne motywacje" P2-C/D/E** (oryginalne triple convergence)
są **trzy konsekwencje** tej jednej zasady, nie trzy independent postulaty.

**Bonus:** automatycznie `f'(1) = -4`, `f''(1) = +8` → `β_PPN = γ_PPN = 1`.

## 2. Source w research/

- `research/closure_2026-04-26/f_psi_principle/results.md` (10.4 KB) —
  TL;DR + 12/12 PASS structural tests
- `research/closure_2026-04-26/f_psi_principle/f_psi_principle.py` (12.6 KB)
  + `.txt` (6.2 KB) — sympy + numerical scan over n ∈ {2,3,4,5,6}
- `research/closure_2026-04-26/f_psi_principle/setup.md` (6.8 KB) —
  audit design (4 warunki + 12 sub-testów)

## 3. Cross-validation z audytem 2026-05-01

**KLUCZOWY KONTEKST:** Audyt 2026-05-01 § V **niezależnie** zamknął ten
sam strukturalny problem (H-B8 — M9.1'' Lagrangean derivation status):

- `research/op-newton-momentum/B8_lagrangean_independence_check.py` (5/5 PASS)
- "Closure via honesty acknowledgment" — explicit mówi że niezależna
  derivation `c₂=-1` z first principles to long-term cycle

**Dwie niezależne ścieżki**:

| Ścieżka | Folder | PASS | Filozofia |
|---------|--------|-----:|-----------|
| Audyt main (§ V) | `op-newton-momentum/B8_lagrangean_independence_check.py` | 5/5 | Independence check Lagrangean (acknowledged calibration) |
| Cross-validation | `closure_2026-04-26/f_psi_principle/f_psi_principle.py` | 12/12 | T-FP principle z `n = deg V = 4` UNIQUE z 4 warunków |

T-FP daje **głębszy strukturalny argument** niż § V — § V mówi "calibration
acknowledged", T-FP daje pierwszą zasadę "skąd `n = 4` UNIQUE".

## 4. Proponowana lokacja w core

| Target | Co dokładnie | Czy istnieje | Diff |
|--------|--------------|--------------|------|
| `core/sek08c_metryka_z_substratu.tex` | Nowa subsection: "T-FP principle: substrate polynomial-degree normalization" — derywacja `f(ψ) = (4-3ψ)/ψ` z V(Φ) bez postulatu | NO (Form-IV ma markery STATUS) | (nowa subsection ~25 linii LaTeX + cross-reference do `eq:M911-form-IV`) |
| `core/sek08c_metryka_z_substratu.tex` przy Form-IV | Notatka "T-FP daje first-principles motivation" zamiast "triple convergence P2-C/D/E" | YES (M9.1'' description) | dopisać 2-3 linie + odsyłacz do nowej subsection |

## 5. Hotspot, który ten wynik adresuje

**H-B8** w [[meta/core/CORE_HOTSPOTS.md]] §B (status: `RESOLVED-AUDIT-FULL`
przez § V). Ten INTAKE proponuje **upgrade** z "honesty acknowledgment"
do "first-principles structural derivation".

Cytat hotspota (oryginalny audyt § B.8): "M9.1'' 'potrójna motywacja'
cyrkularna — 3 zasady (E1/E2/E3) sprowadzają się do jednego ad-hoc
postulatu f(4/3)=0".

T-FP **explicite rozwiązuje** ten zarzut — pokazuje że `f(ψ) = (4-3ψ)/ψ`
nie jest postulatem, lecz UNIQUE rationalem n=4 wymuszonym przez
4 warunki fizyczne (z których 3 mapują na P2-C/D/E).

## 6. Co MUSI być sprawdzone

- [ ] Wynik kompiluje się (sympy/numerical 12/12 PASS confirmed w `results.md`)
- [ ] Brak konfliktu z aktualnym `main.tex` (read-only)
- [ ] `core_compatibility: current` w YAML folderu źródłowego (Sesja 4)
- [ ] `source_of_status` z linkami do `f_psi_principle.py` (Sesja 4)
- [ ] Ścieżka § V (independence check) i T-FP dają zgodne `c₂=-1` derivation
      lub explicit acknowledgment różnicy filozoficznej

## 7. Co się STANIE po akceptacji

1. **Człowiek** edytuje `core/sek08c_metryka_z_substratu.tex` dodając
   T-FP subsection
2. **Człowiek** dopisuje wpis do [[meta/core/CORE_INVENTORY.md]] §1.13
3. **Człowiek** zmienia w `closure_2026-04-26` aggregator YAMLu
   `subfolder_summary` markery
4. **Człowiek** updateuje status H-B8 w `CORE_HOTSPOTS.md` na
   `RESOLVED-AUDIT-FULL + RESOLVED-FULL-PROMOTED`
5. **Człowiek** zamyka INTAKE z `status: ACCEPTED-AND-PROMOTED`

## 8. Decyzja człowieka

- [ ] **ACCEPT** — T-FP promowany jako first-principles motivation w core (data: ____)
- [ ] **DEFER** — wynik OK, promocja przesunięta (powód: ____)
- [ ] **REJECT** — § V acknowledgment wystarcza (powód: ____)
- [ ] **NEEDS-MORE-EVIDENCE** — (lista: ____)

---

## 9. Notatki

### 9.1 PRO

1. T-FP **podnosi epistemiczny status** M9.1'' z "triple convergence
   acknowledged calibration" do "n = 4 UNIQUE z 4 warunków fizycznych"
2. **Bonus** matching `β_PPN = γ_PPN = 1` z `f'(1)=-4, f''(1)=8` jest
   automatic, nie post-hoc fit
3. Rozwiązuje § B.8 kwestię cyrkularności w sposób strukturalny

### 9.2 CONTRA / risks

1. T-FP zakłada `V(Φ)` polynomial degree 4 z aksjomatów (single-Φ Z₂ +
   lokalność/renormalizowalność v2 GL-substrate) — czy ten argument jest
   sam pełnoprawny? Może wymaga `core/sek01_ontologia` konsystencji check
2. **`c₂=-1` nadal calibration** — T-FP daje `f(ψ)`, ale `c₂` to
   property nieliniowości α=2 K(φ)=φ⁴; czy jest to derived w T-FP czy
   tylko `f(ψ)` plus separate `c₂` derivation z α=2?

### 9.3 Rekomendacja agent

**ACCEPT** — T-FP jest najsilniejszą formalną podstawą dla M9.1''
hiperbolicznej formy, jakiej audyt brakuje. § V "honesty acknowledgment"
to słaba formuła zamknięcia; T-FP daje twardszy argument. Promocja
zaowocowałaby M9.1'' jako "M9.1'' z derived f(ψ) z T-FP", co jest
nettowym zwiększeniem siły teorii.

Decyzja należy do człowieka.
