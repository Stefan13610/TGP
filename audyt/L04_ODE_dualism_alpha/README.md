---
title: "L04 — dualizm formulacji ODE: K=g² (α=1) vs K=g⁴ (α=2)"
date: 2026-05-04
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/L04_ODE_dualism_alpha
tags:
  - audit
  - ontology
  - ODE
  - alpha
  - kinetic-formulation
  - canonical-vs-substrate
related:
  - "[[../SUMMARY_2026-05-04.md]]"
  - "[[../L03_K_phi_stability]]"
  - "[[../L05_mass_exponent_drift]]"
  - "[[../../meta/PLAN_DOMKNIECIA_MASTER.md]]"
tgp_status:
  folder_status: audit
  level: L1
  kind: audit
  core_compatibility: partial
  last_reviewed_against_core: 2026-05-04
  may_edit_core: false
  exports_findings: false
  has_needs_file: true
  has_findings_file: false
  open_bridges: ["ODE-canonicalization-decision"]
  depends_on: []
  impacts:
    - "[[../L05_mass_exponent_drift]]"
    - "[[../L03_K_phi_stability]]"
    - "[[../D01_drifting_numbers]]"
  source_of_status:
    - "[[../../meta/PLAN_DOMKNIECIA_MASTER.md]] LP-6, N-1"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# L04 — dualizm formulacji ODE: dwie żywe teorie

## Klasa: LUKA ONTOLOGICZNA / DECYZYJNA • Priorytet: P2

## Diagnoza

TGP ma **dwie konkurujące formulacje ODE solitonu**, z których każda
ma własny `α` (kinetic exponent), własną K(φ), własne ghost properties
i własne `g₀^e`. Obie żyją w plikach jako *równolegle używane*.

LP-6 ([[../../meta/PLAN_DOMKNIECIA_MASTER.md]] § N-1) jawnie tabelaryzuje:

| Formulacja | `α` | K(φ) | g₀^e | r₂₁ | r₃₁ | ghost | stability |
|-----------|-----|------|------|-----|-----|-------|-----------|
| Uproszczone | — | — | 0.899 | 206.77 ✓ | 591 ✗ | — | — |
| Pełne f(g) | — | — | 0.834 | 206.77 ✓ | bariera | — | — |
| **Substratowa** | **1** | **K=g²** | **0.869** | **206.77 ✓** | **3477 ✓** | ghost-free | stabilna g₀>1.3 |
| **Kanoniczna** | **2** | **K=g⁴** | — | 206.77 ✓ | — | ghost przy g*=0.779 | niestabilna |

LP-6 (12/12 PASS) zamyka „dictionary" — ale rekomendacja jest jednoznaczna:

> **Wniosek**: Formulacje są RÓWNOWAŻNE w słabym polu (PPN, κ, n_s, α_s,
> Koide). Różnią się TYLKO w ODE solitonu. K=g² jest preferowana:
> ghost-free, stabilna, reprodukuje pełne spektrum.

I:

> **Rekomendacja**: W sek08 jasno zadeklarować ODE substratowe jako
> kanoniczną formę, inne jako przybliżenia historyczne.

**Ta rekomendacja nie została wcielona.**

## Sprzeczność z README

[[../../README.md]] explicit deklaruje:

> Kinetic coupling α = 2: **algebraic theorem, not a fit**

— to wspiera *kanoniczną* formulację (K=g⁴, α=2). Tymczasem PLAN_DOMKNIECIA
LP-6 oraz LP-4 (k=4 wzór masy) **jednoznacznie używają substratowej
(K=g², α=1)** dla obliczeń liczbowych:

- LP-4 9/9 PASS: substrate ODE K=g² daje r₂₁ = 206.74 (δ=0.013%) z k_eff = 4.000
- Kanoniczna K=g⁴ jest niestabilna dla g₀ > 1.3 ⇒ nie da się policzyć τ.
- Substratowa K=g² daje wszystkie trzy masy leptonów z 0.006%–0.013%.

## To nie jest „dictionary" — to dwie teorie

W zewnętrznym artykule cytującym TGP, czytelnik pyta: *„czy α=2 czy α=1?"*.
Odpowiedź „obie, w zależności co liczysz" to **strukturalna sprzeczność**,
nie elastyczność modelu. Standard QFT vs Standard QFT-z-różnym-K to dwie
różne teorie pól, nie jedna.

Audit wskazuje to wprost (audit § C.2 i inne):

> standard kinetic K = ½(∂φ)² + V wymaga V''(min) > 0. TGP używa
> K(φ) = K_geo·φ⁴ ... niestandardowy operator kinetyczny — wymaga
> osobnego dowodu pełnej linearyzacji EOM.

## Pivot R3 (why_n3 2026-05-01) — zaostrzenie problemu

[[../../research/why_n3/CORRECTIONS_2026-05-01.md]] wprowadza relację
**`m_obs = c · A_tail^(5−α)`** — masa zależy od `α`:

- α=1 (substratowa) ⇒ p=4 (LP-4 oryginalne)
- α=2 (kanoniczna) ⇒ p=3 (R3 nowy)

To znaczy: **fundamentalna predykcja TGP (mass formula) zmienia się
zależnie od wyboru formulacji**. To nie jest „equivalent in weak field"
— to *różne predykcje liczbowe* dla mas leptonowych.

W obecnym stanie:

- Manuskrypt deklaruje α=2 jako algebraic theorem
- Skrypty używają α=1 (K=g²) bo daje stabilne ODE
- Wzór masy `M ∝ A^k` jest k=4 (LP-4) lub k=3 (R3) zależnie od α

Fizyczny rozłam głęboki. Status quo: każdy paper / cykl wybiera tę
formulację, która działa dla danej predykcji. To **systemic problem**,
nie cosmetic.

## Status w audycie

[[../../meta/AUDYT_TGP_2026-05-01.md]] nie wymienia bezpośrednio jako
CRITICAL/HIGH (prawdopodobnie audit przyjął LP-6 dictionary jako
zamknięty). [[../../meta/PLAN_DOMKNIECIA_MASTER.md]] N-1 explicit
acknowledges: *„Nie jest to sprzeczność — r₂₁ jest strukturalny (φ-FP),
niezależny od formy ODE. Bezwzględne g₀ to artefakt normalizacji. Ale:
theory musi się zdecydować na JEDNĄ kanoniczną formę ODE."*

→ status: **decyzja autorska wciąż otwarta**.

## Rekomendacja

Otworzyć dedykowany cykl `op-ODE-canonicalization-decision/`:

### Decision Path

Krótka decyzja autorska: **K=g² (α=1) substratowa** czy **K=g⁴ (α=2)
kanoniczna**?

Argumenty za K=g² (preferowane wg LP-6):

- Ghost-free (g* nieosiągalne)
- Stabilna numerycznie dla g₀ > 1.3
- Reprodukuje wszystkie 3 masy leptonów (m_e, m_μ, m_τ) z <0.02% PDG
- Quark universality: ten sam ODE dla u, c, t i d, s, b
- Cross-sector Koide tryplety mieszane (b,c,t) z δ = 0.42%

Argumenty za K=g⁴ (preferowane wg README):

- α=2 jest „algebraic theorem" (nie fit)
- Standard QFT-style K(φ)=φ²·∂²φ jest bliższy
- Konformalne sprzęganie d=4

### Implementation phase (po decyzji)

1. **Globalny lock formulacji** w `core/sek08*` jako *jedyna kanoniczna*.
2. **Druga formulacja → dodatek historyczny** z explicit „equivalent in
   weak field, different in strong field" — bez pretensji do równoważności.
3. **Mass formula resolution**: jeśli K=g² → k=4 (LP-4 jak było). Jeśli
   K=g⁴ → p=3 (R3 nowy). Jednoznaczność.
4. **Update README**: skoro α=2 deklarowane → albo prawdziwie używać K=g⁴,
   albo poprawić README na K=g² substratową.

**Estymata:** 1 tydzień decyzji + 2 tygodnie propagacji.

## Pliki dotknięte (estymata)

- `core/sek08*` — wszystkie miejsca cytujące K(φ), α, ODE solitonu (~30
  lokalizacji)
- `dodatekA_notacja.tex` — sync notation
- `dodatekB_substrat.tex` — explicit relation z Hamiltonianem v2 GL
- `research/why_n3/` — CORRECTIONS_2026-05-01 wymaga reinterpretacji
- `research/mass_scaling_k4/` — sync z decyzją
- `scripts/lp4_*`, `scripts/lp6_*`, `scripts/a3d_soliton_brannen_r.py` — sync
- README.md — albo zmienić α=2 deklarację, albo prawdziwie używać K=g⁴

## Open NEEDS

Patrz [[NEEDS.md]].

## Cross-references

- [[../SUMMARY_2026-05-04.md]] §II.L4
- [[../PRIORITY_MATRIX.md]] klaster D
- [[../../meta/PLAN_DOMKNIECIA_MASTER.md]] LP-6, § N-1, LP-4
- [[../../research/why_n3/CORRECTIONS_2026-05-01.md]]
- [[../L05_mass_exponent_drift]] (zależność: wybór α determinuje k vs p)
- [[../L03_K_phi_stability]] (zależność: K(φ) określa spektralną analizę)
- [[../D01_drifting_numbers]] (g₀^e drift)
