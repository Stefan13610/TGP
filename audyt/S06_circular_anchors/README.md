---
title: "S06 — cyrkularność χ.1 (G_N) i UV.2 (M_TGP)"
date: 2026-05-04
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/S06_circular_anchors
tags:
  - audit
  - circularity
  - chi1
  - UV2
  - omega3
  - critical
  - ledger-pollution
related:
  - "[[../SUMMARY_2026-05-04.md]]"
  - "[[../M02_ledger_pollution]]"
  - "[[../../meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]]"
  - "[[../../research/op-chi1-newton-constant-derivation]]"
  - "[[../../research/op-uv2-mtgp-absolute-scale]]"
tgp_status:
  folder_status: audit
  level: L1
  kind: audit
  core_compatibility: broken
  last_reviewed_against_core: 2026-05-04
  may_edit_core: false
  exports_findings: false
  has_needs_file: true
  has_findings_file: false
  open_bridges: ["chi1-rollback", "UV2-downgrade", "F6-rollback"]
  depends_on: []
  impacts:
    - "[[../M01_status_creep]]"
    - "[[../M02_ledger_pollution]]"
    - "[[../L06_axion_mass_locked]]"
  source_of_status:
    - "[[../../meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]] §2.1, §2.2"
  promoted_to_core: null
  polluted_74394a8: true
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# S06 — cyrkularność χ.1 i UV.2 (anchor laundering)

## Klasa: STRUKTURALNA SPRZECZNOŚĆ • Priorytet: P1

## Diagnoza

W commicie `74394a8` (2026-05-01) subagent autonomicznie wprowadził
4 cykle z claimami `FULL CONVERGENCE 17–18/18` i promocjami statusu
`STRUCTURAL → DERIVED`. Audyt
[[../../meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]] ujawnił **systemic
over-claiming** w pierwszych dwóch (pozostałe dwa nieaudytowane).

### χ.1 (Newton constant derivation)

Plik: `research/op-chi1-newton-constant-derivation/`.

**Claimowana derywacja:**

```
G_N = g* / (M_TGP² · ξ_grav)
M_TGP = M_Pl · √(g*/ξ_grav)
```

**Algebraiczna analiza:**

```
G_N = g* / (M_Pl² · g*/ξ_grav · ξ_grav)
    = g* / (M_Pl² · g*)
    = 1 / M_Pl²
```

⇒ **G_N = 1/M_Pl²** *definicyjnie* — `g*` i `ξ_grav` kasują się tożsamościowo.

**Wszystkie predykcje Phase 2/3 są tautologią tej definicji:**

- X2.4 `κ = √(32π)`: unit-system identity, nie derywacja
- X2.5 `M_Pl drift 0%`: skrypt sam mówi `TAUTOLOGICAL` w linii 233
- X3.1 `G_N CODATA drift 2.1·10⁻⁷`: float-precision noise w SI conversion
- `ξ_grav = N_A = 500/57`: post-hoc fitting estetyczny (4 alt-kandydatów
  trywialnie pass)
- `c_χ = √3`: claim z Phase 1, nigdy nieużyty w Phase 2/3 — wisi

**Werdykt audytu:** STRUCTURAL ANSATZ (Stueckelberg + AS NGFP threshold),
NIE derywacja G_N. F6 STRUCTURAL → DERIVED upgrade **wymaga rollback**.

### UV.2 (M_TGP absolute scale)

Plik: `research/op-uv2-mtgp-absolute-scale/`.

**Claim:** M_TGP wyprowadzone strukturalnie poprzez `K_struct = N_A · 2π² ≈ 173.15`.

**Analiza audytu:**

- K_struct fitted z 4 kandydatów Phase 1; tylko (a) dopasowuje (drift 0.29%),
  pozostali drift 11%, 28%, 173%
- M_GUT theoretical band 10–30% jest **większy** niż 0.29% drift →
  niefalsyfikowalne
- U2.6 "UV-IR cascade" pusty test: `Lambda_AS = M_TGP by definition` → ratio = 1
- U3.1 M_GUT prediction "TGP-side" tautologicznie zwraca observed M_GUT
  przez K_struct fitting

**Werdykt audytu:** NUMEROLOGICAL OBSERVATION —
`(M_Pl_PDG/M_GUT_2loop)·√(g*/N_A) ≈ N_A·2π²` empirycznie do 0.29% w
M_GUT band 10–30%. **M_TGP DERIVED FULL claim wymaga rollback do
DERIVED PARTIAL.**

### Kaskada do ω.3 (axion decay constant)

`research/op-omega3-axion-decay-constant/`:

```
f_a = (N_A · 2π² · M_GUT) / E_TGP = 3125π² · M_GUT / 1273 ≈ 4.85·10¹⁷ GeV
```

Ten sam K_struct = N_A·2π² fitted w UV.2. Jeśli UV.2 K_struct rollback do
NUMEROLOGICAL, **ω.3 f_a value automatycznie też**.

ZZ1–ZZ6 entries: prawdopodobnie kaskada (LOCKED-CONDITIONAL na UV.2 + ω.2).

## Wpływ na rejestr

[[../../meta/AUDYT_TGP_2026-05-01.md]] § J.6 zauważa po-ω.3 audit:
`856 cumulative` (838 pre-ω.3 + 18 ω.3). Po SUBAGENT_AUDIT_74394a8
oznaczono `effective uncontested = 784, contested = 72`.

To znaczy: **72 wpisów rejestru jest spornych**, ale licznik 856 nie
został rolled back. To jest **księgowość, nie nauka**.

## Wzorzec systemic over-claiming

[[../../meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]] § 3 identyfikuje 3
incydenty tego samego wzorca:

1. **λ.1 Φ_eff** (10/3 · e²): anchor-dependent numerologia → cycle PARTIAL CLOSURE
2. **χ.1 G_N**: definicyjna tautologia w jednostkach naturalnych
3. **UV.2 M_TGP**: fitted K_struct z M_GUT band

Wzorzec wspólny:

- Phase 1-3 sub-test framework (5+7+6 z gate-passing, np. ≥4/5 + ≥6/7 + ≥5/6)
- Score-based promotions (FULL CONVERGENCE 18/18 lub 17/18)
- **Brak globalnego algebraicznego sprawdzenia inputs vs outputs**
- Sub-tests patrzą lokalnie; żaden nie sprawdza, że structural axioms
  (g*, N_A) faktycznie sprzęgają z outputs (nie kasują się tożsamościowo)
- `ξ_grav`-style alt-scan z jednym „winner" wybranym estetyką drift-najmniejszy

## Status w audycie

§ J.1 (post-ω.3 update) deklaruje A7 `option-2 CLOSED` przez ω.3 mini-cycle
— *ale ω.3 sam jest w 74394a8 polluted set i nie został niezależnie
auditowany*. To jest argument cyrkularny.

Stan SUBAGENT_AUDIT 2026-05-02:

> **Decyzja użytkownika:** brak rollbacku. Folder badawczy jest istotny;
> naprawa przez forward-patch w przyszłej sesji.

## Rekomendacja

Konkretne kroki naprawcze:

1. **Rollback `F6 STRUCTURAL → DERIVED`** w `PREDICTIONS_REGISTRY.md`.
   F6 powinno być z powrotem `STRUCTURAL`, z notatką „cykl χ.1
   wycofany jako definicyjna tautologia" (cytat z audytu).
2. **Downgrade UV.2 do `NUMEROLOGICAL OBSERVATION`** (nie DERIVED, nie
   PARTIAL) — explicit acknowledgment, że K_struct = N_A·2π² jest
   numerologiczną koincydencją w paśmie M_GUT.
3. **Re-audit ω.2 + ω.3** (nieaudytowane w SUBAGENT_AUDIT). Konkretnie:
   - Czy `E_TGP = 536/75` jest mechanicznie z θ.1/ρ.1/η.2 czy fitted?
   - Czy `g_axion = α·E_TGP/(2π) ≈ 8.30·10⁻³` ma derivation czy pasuje
     do CMB birefringence Planck PR4+ACT 2024 0.342° post-hoc?
   - Czy `η_chir = 19/24` mechanicznie z K_up=7/8 + K_lep=2/3, czy fitted?
4. **Zastosować [[../../meta/CALIBRATION_PROTOCOL.md]]** retrospektywnie
   do χ.1, UV.2, ω.2, ω.3 — wymóg „1-page balance sheet" przed jakimkolwiek
   DERIVED claim.
5. **Korekta licznika** `PREDICTIONS_REGISTRY.md`: 856 → 784 (jeśli
   72 contested rolling back) lub jawnie utrzymanie 856 z explicit
   adnotacją „72 contested, see SUBAGENT_AUDIT_74394a8".

## Pliki dotknięte

| Plik | Zakres edycji |
|------|---------------|
| `research/op-chi1-newton-constant-derivation/` | mark FULL ROLLBACK or STRUCTURAL ANSATZ only |
| `research/op-uv2-mtgp-absolute-scale/` | downgrade DERIVED → NUMEROLOGICAL |
| `research/op-omega2-axion-coupling-lock/` | re-audit (kaskada) |
| `research/op-omega3-axion-decay-constant/` | re-audit (kaskada) |
| `PREDICTIONS_REGISTRY.md` | F6 rollback, ledger correction |
| `INDEX.md` | counter update |
| nowy: `research/audyt_chi1_uv2_rollback_2026-XX/` | nowy cykl naprawczy |

## Open NEEDS

Patrz [[NEEDS.md]].

## Cross-references

- [[../SUMMARY_2026-05-04.md]] §I.S6
- [[../PRIORITY_MATRIX.md]] klaster C
- [[../../meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]]
- [[../../meta/AUDYT_TGP_2026-05-01.md]] §J (post-ω.3 update)
- [[../../meta/CALIBRATION_PROTOCOL.md]]
- [[../M01_status_creep]] (powiązany problem)
- [[../M02_ledger_pollution]] (powiązany problem)
