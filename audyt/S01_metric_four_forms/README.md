---
title: "S01 — cztery sprzeczne formy metryki w sek08c"
date: 2026-05-04
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/S01_metric_four_forms
tags:
  - audit
  - metric
  - sek08c
  - critical
  - structural-contradiction
related:
  - "[[../SUMMARY_2026-05-04.md]]"
  - "[[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]"
  - "[[../../TGP_FOUNDATIONS.md]]"
  - "[[../../meta/AUDYT_TGP_2026-05-01.md]]"
tgp_status:
  folder_status: audit
  level: L2
  kind: audit
  core_compatibility: broken
  last_reviewed_against_core: 2026-05-04
  may_edit_core: false
  exports_findings: false
  has_needs_file: true
  has_findings_file: false
  open_bridges: ["sek08c-rewrite", "M9.1''-canonical-form"]
  depends_on: []
  impacts:
    - "[[../S02_volume_element_M9]]"
    - "[[../S03_beta_PPN_convention]]"
  source_of_status:
    - "[[../../meta/AUDYT_TGP_2026-05-01.md]] §A.1"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# S01 — cztery sprzeczne formy metryki w `sek08c`

## Klasa: STRUKTURALNA SPRZECZNOŚĆ • Priorytet: P1

## Diagnoza

W jednym pliku rdzenia
[[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]
**współistnieją cztery wzajemnie sprzeczne formy metryki**, każda z
własnym „twierdzeniem jedyności" lub „derywacją". Audit z 2026-05-01
oznaczył je inline komentarzami `STATUS`, ale **wyprowadzenia same
pozostały w pliku** — czytelnik dostaje muzeum sprzecznych derywacji
z naklejkami.

### Forma (I) — boxed power

```
g_tt = -c²/ψ,    g_rr = ψ
```

- Lokalizacja: `eq:metric-full-derived` lin. 154
- Status: **FALSIFIED 2026-04-25** (M9.1: β_PPN=4 vs Mercury 1±10⁻⁴,
  niezgodność 3·10⁴σ)
- Charakter: zaboxowane jako „twierdzenie derywacyjne"
- Obecność w pliku: nadal kanoniczna w sek08c lin. 131–200; używana w
  sek08a do wariacji + κ=3/(4Φ_0)

### Forma (II) — eksponencjalna

```
g_tt = -c²·e^{-2U}
```

- Lokalizacja: lin. 208–211
- Status: **OBSOLETE / PRZEJŚCIOWA**
- Charakter: była używana dla `β_PPN=1` *dokładnie*
- Brak: niezależnego wyprowadzenia z budżetu; sprzeczna z (I) i (III)

### Forma (III) — antypodalne twierdzenie jedyności

```
f = ψ^{-1/2},   h = ψ^{+1/2}
```

- Lokalizacja: `thm:antipodal-uniqueness` lin. 345–371
- Status: **NIESPÓJNA WEWNĘTRZNIE** (M9.1' §2.1, 2026-04-25)
- Charakter: „twierdzenie jedyności" — linearyzacja daje γ=1 ✓, ale
  pełne O(ε²) daje **β_metric=3** (nie 1 jak deklaruje
  `rem:antipodal-implications`); z α=2 dynamiką (c₂=-1) → β_PPN=7
- „Jedyność" nie obejmuje pełnej nieliniowości

### Forma (IV) — M9.1'' hiperboliczna

```
g_tt = -c²·(4-3ψ)/ψ,    g_rr = ψ/(4-3ψ)
```

- Lokalizacja: [[../../TGP_FOUNDATIONS.md]] lin. 64–69
- Status: **OPEN, KANONICZNA wg FOUNDATIONS** post-2026-04-25
- Charakter: daje β_PPN=γ_PPN=1 algebraically (master formula z c₂=-1)
- **Krytyczne:** ta forma **NIE jest fizycznie obecna w `sek08c.tex`**,
  tylko w nagłówkowym komentarzu i w `research/op-newton-momentum/M9_1_pp_*`

## Wpływ na rdzeń

| Element | Forma używana | Status post-pivot M9.1'' |
|---------|---------------|---------------------------|
| `sek08a` wariacja → κ=3/(4Φ_0) | (I) `√(-g)=c·ψ` | UNIEWAŻNIONE — patrz [[../S02_volume_element_M9]] |
| `Φ-EOM` (`prop:field-eq-from-action`) | (I) | UNIEWAŻNIONE |
| `M9.2` m_field | (I) | UNIEWAŻNIONE 5/5 PASS |
| `M9.3` Peters-Mathews | (I) | UNIEWAŻNIONE 5/5 PASS |
| `thm:antipodal-uniqueness` | (III) | TWIERDZENIE BŁĘDNE (β=3) |
| linearyzacja `sek08c` lin. 183–191 | (I) → β=2 | SPRZECZNE Z (IV) β=1 |
| TGP_FOUNDATIONS §3 deklaracja | (IV) | NIE WYPROWADZONE w sek08c |

## Status w `meta/AUDYT_TGP_2026-05-01.md`

§A.1 — **CRITICAL**, oznaczony „CLOSED markery inline w sek08c.tex"
(§ M.6). Faktyczny stan: **CLOSED-annotation-only**. Audit dodał
nagłówkowy komentarz lin. 64–129 i inline `STATUS` markery przy
formach (I)–(III), ale:

- Twierdzenia jedyności pozostały w pliku.
- Wyprowadzenia kroku 3 mostu pozostały w pliku.
- Forma (IV) M9.1'' **NIE została do pliku dodana** jako kanoniczna
  derywacja.
- Linearyzacja β=2 (lin. 183–191) pozostała sprzeczna z M9.1'' β=1.

Audit § M.4 explicit acknowledges: `√(-g) w M9.x derivations = OBSOLETE`.

## Rekomendacja

**Pełny rewrite `sek08c.tex` z M9.1'' jako jedyną kanoniczną formą.**

1. Sekcja I (1–80% pliku): pełna derywacja M9.1'' z budżetu informacyjnego
   substratu (forma (IV)) — wszystkie kroki, łącznie z `√(-g)=c·ψ/(4-3ψ)`,
   PPN expansion, master formula β=1 z explicit acknowledged kinetic
   correction.
2. Sekcja II (dodatek historyczny): formy (I), (II), (III) z explicit
   FALSIFIED/OBSOLETE/NIESPÓJNE — **bez** twierdzeń jedyności, jako
   „pivot record" 2026-03 → 2026-04-25.
3. Linearyzacja lin. 183–191 → usunąć lub przepisać dla M9.1''
   (powinna dać β=2 dla (I), β=1 dla (IV) — explicit różnicowanie).
4. `thm:antipodal-uniqueness` → przepisać jako *twierdzenie pod-warunkowe*
   („pod założeniami liniowości β-warunek nie obejmuje O(ε²)") albo
   wycofać.

**Estymata:** 1–2 tygodnie pracy LaTeX + 1 tydzień review.

**Zależności wychodzące:**
- [[../S02_volume_element_M9]] — re-run M9.x z poprawnym `√(-g)`
- [[../S03_beta_PPN_convention]] — formalna derywacja c₂=-1

## Pliki dotknięte

| Plik | Zakres edycji |
|------|---------------|
| `core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex` | full rewrite |
| `core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex` | sync `√(-g)` z (IV) |
| `core/sek08_formalizm/sek08_formalizm.tex` | sync `rem:hierarchia-sektorow` |
| `TGP_FOUNDATIONS.md` lin. 64–69 | (już canonical, OK) |

## Open NEEDS

Patrz [[NEEDS.md]].

## Cross-references

- [[../SUMMARY_2026-05-04.md]] §I.S1
- [[../PRIORITY_MATRIX.md]] klaster A
- [[../../meta/AUDYT_TGP_2026-05-01.md]] §A.1, §M.1–M.6
- [[../../research/op-newton-momentum/M9_1_pp_P3_results.md]]
