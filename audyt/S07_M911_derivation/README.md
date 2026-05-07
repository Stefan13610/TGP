---
title: "S07 — Metryka M9.1'' jako postulat vs derywacja z pierwszych zasad"
date: 2026-05-06
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/S07_M911_derivation
tags:
  - audit
  - structural
  - metric
  - M911
  - derivation-vs-postulate
  - emergent-gravity
  - S07
  - EXT-3
related:
  - "[[../EXTERNAL_REVIEW_2026-05-06.md]]"
  - "[[../README.md]]"
  - "[[../PRIORITY_MATRIX.md]]"
  - "[[../S01_metric_four_forms/README.md]]"
  - "[[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
tgp_status:
  folder_status: audit
  level: L1
  kind: audit
  core_compatibility: partial
  last_reviewed_against_core: 2026-05-06
  may_edit_core: false
  exports_findings: false
  has_needs_file: false
  has_findings_file: false
  open_bridges: ["op-metric-from-substrate", "op-metric-uniqueness", "op-entropic-metric"]
  depends_on: ["S01 CLOSED via G.0 2026-05-02"]
  impacts: ["sek08c metryka", "sek08a M9.1'' canonical", "TGP_FOUNDATIONS § 2 status (E)"]
  source_of_status:
    - "[[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-3"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-06
---

# S07 — Metryka M9.1'' jako postulat vs derywacja z pierwszych zasad

## Klasa: STRUKTURALNA (zewnętrzna recenzja EXT-3) • Priorytet: **P2**

## Diagnoza (z EXT-3)

Forma metryki efektywnej była **iterowana**:

| Iteracja | Forma g_tt | Status | Powód odrzucenia |
|----------|-----------|--------|------------------|
| Forma I (potęgowa) | −c²/ψ | DEPRECATED 2026-04-25 | β_PPN = 4 vs obs 1, **3·10⁴σ** |
| Forma II (eksponencjalna) | −c² e^(−2U) | DEPRECATED | równoważna do O(U), różnica na 2PN+ |
| **Forma III (M9.1'' hiperboliczna)** | **−c²(4−3ψ)/ψ** | **CANONICAL 2026-04-25** | **γ_PPN = β_PPN = 1 *exact*** |

(`sek08a_akcja_zunifikowana.tex` lin. 196–251, M9.1'' canonical replacement).

**Problem:** "exact by construction" jest **podejrzanie zręcznym
wyborem**. Master formula (P23 sympy LOCK 5/5 PASS) reprodukuje
γ = β = 1, ale **nie ma argumentu, który *przed* zobaczeniem PPN
forsowałby hiperboliczną formę**. To znaczy, że hipoteza została
*dobrana tak*, by dawała poprawny PPN.

## Pliki dotknięte

- [[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]
  (G.0 preamble, eq:metric-M911-canonical)
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]
  lin. 186–251, 240–290 (canonical M9.1'' replacement)
- [[../../research/op-newton-momentum/]]
  (M9.1, M9.1', M9.1'' historia)
- [[../S01_metric_four_forms/README.md]] (cztery formy metryki —
  zamknięte przez G.0 closure, ale bez derywacji *z pierwszych zasad*)

## Problemy jakie rodzi

1. **Status epistemiczny.** TGP_v1 deklaruje hierarchię (W/E/P/H —
   wyprowadzenie/empiryczna/propozycja/hipoteza). Metryka M9.1'' ma
   status (E) w `TGP_FOUNDATIONS.md` § 2 — "metryka efektywna,
   emergencja w limicie GR, FRW, PPN". **Ale emergencja nie znaczy
   derywacja**: M9.1'' jest *postulatem*, którego konsekwencje są
   sprawdzalne. To bliżej (P) niż (E).

2. **Brak fundamentu z H_Γ.** `sek02:189–203` przyznaje wprost:
   > "Dokładna postać tej relacji [g_eff = g_eff[Φ, ∂Φ, …]] jest
   > **otwartym problemem TGP**."

   Hipoteza minimalna (lin. 210–221) ma konkurentów. M9.1'' pojawiła
   się jako **trzecia iteracja** po falsyfikacji form I i II — zatem
   jest *empirycznie wybrana*, a nie *strukturalnie wyprowadzona*.

3. **Kruchość liczbowa.** γ_PPN = β_PPN = 1 *exact* jest świetne, ale
   2PN współczynniki c_2 = −1, c_3 = +5/3, c_4 = −10/3 są specyficzne
   dla M9.1''. Oznacza to przewidywanie:
   |Δg_tt| = (5/6) U³ deviation od GR (testowalne LIGO 3G — **EXT-5/T01**).
   Ale *dlaczego ta wartość*? **Bez fundamentu, jeśli LIGO 3G da inną
   wartość, M9.1'' upada — i TGP nie ma kolejnej iteracji *gotowej
   z konstrukcji*.**

4. **Filozoficzny dług.** "Metryka emergentna" jest sztandarowym
   hasłem TGP. Ale jeśli emergentna metryka jest **postulowana
   funkcyjnie** (g_eff[ψ] = − c_0² (4−3ψ)/ψ jako ansatz), to różnica
   między TGP a teorią skalarno-tensorową z konkretnym ansatzem
   sprzężenia konformalnego jest *retoryczna*, nie strukturalna.
   `TGP_FOUNDATIONS.md` § 5.1 jawnie odrzuca "teorię skalarno-tensorową
   w stylu Brans-Dicke / Horndeski" — ale operacyjnie M9.1'' jest
   szczególnym przypadkiem skalarno-tensorowego ansatzu.

## Potencjalne ścieżki domknięcia

### Ścieżka A — derywacja M9.1'' z budżetu substratu

`sek08c_metryka_z_substratu` ma w nazwie "z substratu", co sugeruje
ambicję derywacji. **Realizacja:** explicit coarse-graining H_Γ na
graf regularny, identyfikacja efektywnej metryki przez Wilsonowską
RG flow.

**Cykl proponowany:** `op-metric-from-substrate/` z wyprowadzeniem
w 4 krokach:
1. H_Γ → continuum action
2. effective metric ansatz
3. matching to M9.1''
4. punkt 3 musi być *konsekwencją*, nie *wyborem*

### Ścieżka B — wariacyjne kryterium minimalności

M9.1'' jako **jedyna** metryka spełniająca:
- (i) γ_PPN = β_PPN = 1
- (ii) hiperboliczność horyzontu ψ = 4/3 jako natural BH cutoff
- (iii) konformalna struktura kompatybilna z α = 2 selection

Jeśli (i)+(ii)+(iii) okażą się **wystarczającymi i konstruktywnymi**
warunkami, to M9.1'' staje się *wybrana z konstrukcji*, a nie
*dobrana*.

**Cykl proponowany:** `op-metric-uniqueness/`

### Ścieżka C — przyznanie statusu (P) w hierarchii

**Najszczersze rozwiązanie:** zmienić status M9.1'' z (E) na (P) —
*propozycja*, nie *emergencja*. Przepisać `TGP_FOUNDATIONS.md` § 2,
by hierarchia była:
- (W) substrat
- (W) Φ-EOM
- **(P) M9.1''**
- (P/H) materia

**Cena:** TGP traci jeden z głównych argumentów marketingowych
("emergencja metryki"), zyskuje **zaufanie naukowe**.

### Ścieżka D — derywacja z budżetu Verlindego/Padmanabhana

W tradycji emergentnej grawitacji (sek08c referuje do Sakharov,
Verlinde) metryka wynika z entropic budgetu na horyzoncie. Spróbować
explicit oblicić, czy hiperboliczne g_tt = −(4−3ψ)/ψ wynika z
założenia:
- powierzchnia Hubble'a + entropia ∝ A_horizon / 4G + Z₂ substratowa

**Cykl proponowany:** `op-entropic-metric/`

## Rekomendowany priorytet

**P2 — wysoki.** Bez derywacji M9.1'' z fundamentu, TGP jest
**operacyjnie skalarno-tensorową teorią ze szczególnym ansatzem**.
To nie unieważnia jej, ale pozbawia jednego z głównych roszczeń
ontologicznych.

## Powiązanie z istniejącym audytem

[[../S01_metric_four_forms/README.md]] zamknięty CLOSED-RESOLVED via
G.0 2026-05-02. Ale "RESOLVED" oznacza tu *spójność wewnętrzna*
(4 formy → 1 forma), **nie *derywacja z pierwszych zasad***.
S07 idzie głębiej — wymaga *derywacji*, nie tylko *spójności*.

## Cross-references

- [[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-3 — recenzja źródłowa
- [[../README.md]] — indeks audytu
- [[../PRIORITY_MATRIX.md]] — do update z S07 P2
- [[../S01_metric_four_forms/README.md]] — closed via G.0 (spójność wewnętrzna)
- [[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]
- [[../../research/op-newton-momentum/]] — M9.1, M9.1', M9.1'' historia
- [[../../TGP_FOUNDATIONS.md]] § 2 hierarchia W/E/P/H
