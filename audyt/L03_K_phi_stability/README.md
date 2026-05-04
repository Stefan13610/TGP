---
title: "L03 — V''(1)<0 vs K(φ)=K_geo·φ⁴ stability"
date: 2026-05-04
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/L03_K_phi_stability
tags:
  - audit
  - ontology
  - stability
  - kinetic
  - V_eff
related:
  - "[[../SUMMARY_2026-05-04.md]]"
  - "[[../../axioms/notacja/dodatekA_notacja.tex]]"
  - "[[../../meta/AUDYT_TGP_2026-05-01.md]]"
tgp_status:
  folder_status: audit
  level: L1
  kind: audit
  core_compatibility: partial
  last_reviewed_against_core: 2026-05-04
  may_edit_core: false
  exports_findings: false
  has_needs_file: false
  has_findings_file: false
  open_bridges: ["full-spectral-analysis"]
  depends_on: []
  impacts: []
  source_of_status:
    - "[[../../meta/AUDYT_TGP_2026-05-01.md]] §C.2"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# L03 — V''(1)<0 vs K(φ)=K_geo·φ⁴ stability

## Klasa: LUKA ONTOLOGICZNA • Priorytet: P2

## Diagnoza

W kanonicznym QFT warunek stabilności wymaga:

```
V''(φ_min) > 0
```

co zapewnia `m_eff² > 0` przy linearyzacji wokół minimum, brak ghost
modes, dodatnio-określone propagatory.

TGP ma:

```
V(φ) = (β/3)·φ³ - (γ/4)·φ⁴, β=γ (warunek próżni)
V'(φ) = β·φ² - γ·φ³ = γφ²(1-φ)  ⇒ V'(1)=0 ✓
V''(φ) = 2β·φ - 3γ·φ² = γφ(2-3φ)  ⇒ V''(1) = -γ < 0
```

⇒ **V''(1) jest *ujemne*** — formalnie `φ=1` jest *maksimum* potencjału,
nie minimum! To wygląda jak naruszenie warunku stabilności.

## „Rozwiązanie": niestandardowy kinetic prefactor

TGP używa **niestandardowego operatora kinetycznego**:

```
K(φ) = K_geo · φ⁴
L_kin = (1/2) · K(φ) · g_eff^μν · ∂_μφ · ∂_νφ
```

Po linearyzacji `φ = 1 + δφ` wokół φ=1:

```
m_eff² = -V''(1) / K(1) = γ / K_geo > 0
```

⇒ **dodatnia masa**, jeśli K_geo > 0. Stability uratowana.

## Co brakuje

Audit § C.2 zamyka *annotation-only*:

> Stabilność potwierdzona przez pełną linearizację δφ EOM:
> m_eff² = -V''(1)/K(1) = γ/K_geo > 0 (po właściwym rescaling). Sympy
> LOCK linearizacji w sek08a krok 2 + dodatekA:404-412 jest wynaleziony
> nie z naive Taylor V, ale z covariant EOM expansion.

**Brakuje w `core/`:**

1. **Pełnej analizy spektralnej** — czy w okolicach `φ=1` istnieją tylko
   massive modes (m² > 0), czy także ghost modes (m² < 0 z odwrotnym
   znakiem K)? K(φ)=K_geo·φ⁴ ma `K(0) = 0` i `K' (0) = 0` — to *miękko
   znika* przy zerze, co może generować osobliwości w propagatorze przy
   `φ → 0`.

2. **Mode counting w obecności złamania Z₂** — Z₂ symmetry breaking daje
   1 Goldstone (dla ciągłego), ale Z₂ jest dyskretne ⇒ nie ma Goldstone
   bosona. Powinien być tylko 1 massive mode. Czy jest?

3. **Tachyonic projection check** — czy linearyzacja na `Φ_eq[ρ]` (z
   źródłem) zachowuje stabilność dla wszystkich ρ ≥ 0?

4. **Ghost wall przy `φ < g* ≈ 0.78`** wspomniany w LP-6 (PLAN_DOMKNIECIA
   N-1) — gdzie K(φ) traci pozytywność. Dla formy kanonicznej K=g⁴
   istnieje ghost przy g* = exp(-1/(2α)) = 0.779 (dla α=2). Skutki
   fizyczne tej ściany ghost'ów na predykcje?

## Wpływ na predykcje

- **Spektrum mas leptonowych**: linearyzacja `Φ-EOM` dla solitonu daje
  masę `m_e ∝ A_tail^k` (LP-4). Stabilność solitona zależy od pełnej
  analizy spektralnej.
- **GW propagation**: dyspersja fluktuacji `δΦ` propagujących na tle
  `Φ_eq` zależy od `K(φ)·∂² + V''(φ)`. Jeśli `K(φ)·V''(φ) < 0` w jakimś
  obszarze, mamy tachyonic propagation.
- **N-body symulacje**: stability solitonowych konfiguracji w obecności
  multiple sources wymaga gwarancji że nigdzie nie wpadamy w
  unstable/ghost regime.

## Status w audycie

§ C.2 — MEDIUM, „CLOSED annotation-only". Faktyczny stan: pełna analiza
spektralna nie jest wykonana. Cytat „covariant EOM expansion" w sek08a
krok 2 i dodatkuA:404–412 jest tylko *zarysem*, nie pełnym dowodem.

## Rekomendacja

Otworzyć dedykowany cykl `op-spectral-analysis-Phi/`:

### Phase 1 — pełna spektralna analiza linearyzacji

Dla `Φ-EOM` linearyzowanego wokół `Φ_eq[ρ]` (różne profile: vacuum,
Yukawa, soliton):

- Diagonalize `K(φ)·∂² + V''(φ)` operator
- Sprawdzić znak wszystkich eigenvalues
- Identify wszystkie modes (scalar massive, possible breathing, etc.)
- Confirm m² > 0 wszędzie gdzie φ > 0

### Phase 2 — ghost wall analysis

Dla `g* ≈ 0.78` (kanoniczna K=g⁴) lub wybranej formulacji (po L04
decyzji):

- Pokazać, że w fizycznych konfiguracjach (`Φ_eq` od materii) `φ` nigdy
  nie wpada w ghost regime
- Lub: pokazać, że ghost wall ma fizyczną interpretację (boundary między
  fazami)

### Phase 3 — tachyonic check

Linearyzacja na różnych tle (vacuum, Newtonian, strong field) — żaden nie
pojawia się tachyon mode.

**Estymata:** 4–6 tygodni (formalna analiza + numerical verification).

## Pliki dotknięte

| Plik | Zakres edycji |
|------|---------------|
| `sek08a_akcja_zunifikowana.tex` krok 2 | rozszerzyć linearization analysis |
| `dodatekA_notacja.tex` lin. 404–412 | sync z full spectral analysis |
| nowy: `research/op-spectral-analysis-Phi/` | nowy cykl |

## Cross-references

- [[../SUMMARY_2026-05-04.md]] §II.L3
- [[../PRIORITY_MATRIX.md]] klaster D
- [[../../meta/AUDYT_TGP_2026-05-01.md]] §C.2
- [[../L04_ODE_dualism_alpha]] (powiązany — różne K(φ) formy)
