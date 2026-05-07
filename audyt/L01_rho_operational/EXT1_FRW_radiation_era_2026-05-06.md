---
title: "L01 EXT-1 — Kosmologia radiacyjna z varying c, ℏ, G (BBN/CMB)"
date: 2026-05-06
parent: "[[README.md]]"
type: audit-extension
tgp_owner: audyt/L01_rho_operational
tags:
  - audit
  - external-review
  - L01-extension
  - cosmology
  - radiation-era
  - varying-constants
  - BBN
  - CMB
  - P1
related:
  - "[[README.md]]"
  - "[[POST_ACTION_UPDATE_2026-05-04.md]]"
  - "[[../EXTERNAL_REVIEW_2026-05-06.md]]"
  - "[[../../research/op-L01-rho-stress-energy-bridge-2026-05-04/NEEDS.md]]"
tgp_status:
  folder_status: audit
  level: L1
  kind: audit-extension
  core_compatibility: review-only
  last_reviewed_against_core: 2026-05-06
  may_edit_core: false
  exports_findings: false
  has_needs_file: false
  has_findings_file: false
  open_bridges: ["op-FRW-radiation-era-varying-c", "op-BBN-TGP", "op-CMB-TGP"]
  depends_on: ["L01_rho_operational EXECUTED 2026-05-04"]
  impacts: ["L01 NEEDS N1 (quantum EM trace anomaly)", "L01 NEEDS N2 (QCD vacuum)"]
  source_of_status:
    - "[[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-1 v2"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-06
---

# L01 EXT-1 — Kosmologia radiacyjna z varying c, ℏ, G

## Klasa: ROZSZERZENIE L01 (zewnętrzna recenzja) • Priorytet: **P1 OTWARTE RYZYKO**

## Diagnoza (z EXT-1 v2)

L01 formal definition (cykl 2026-05-04) zamknął operacyjną definicję `ρ`:

```
ρ(x) ≡ −T^μ_μ(x) / c_0²
```

W połączeniu z aksjomatami `ax:c–ax:G` (`sek04_stale.tex` lin. 29, 250–254):

```
c(Φ) = c₀ √(Φ₀/Φ)     ℏ(Φ) = ℏ₀ √(Φ₀/Φ)     G(Φ) = G₀ Φ₀/Φ
```

→ w TGP **c₀, ℏ₀, G₀ NIE są stałymi fundamentalnymi**, lecz wartościami
referencyjnymi w obecnej epoce; same `c, ℏ, G` są funkcjami pola Φ.

**Strukturalnie:** promieniowanie nie generuje Φ przez L_mat. T^μ_μ_EM = 0
jest twierdzeniem matematycznym (Weyl-niezmienniczość Lagrange'a EM
−¼F_μν F^μν √(−g) w 4D), niezależnym od tego, czy w mianowniku stoi
c₀ czy c(Φ). **To trzyma się.**

**Fenomenologicznie:** ALE H_TGP(z) w erze radiacyjnej **NIE musi być
równe H_GR(z)** — pochodzi z Φ-EOM w FRW, gdzie dynamika Φ +
varying c(Φ), ℏ(Φ), G(Φ) tworzą własną strukturę kinematyczną
ekspansji. Foton w TGP jest **pasażerem geometrii** (porusza się po
g_eff[Φ]), ale geometria ewoluuje przez Φ-dynamikę.

## Status: P1 OTWARTE RYZYKO (nie *przegrane a priori*)

**Subiektywna ocena recenzenta:** prawdopodobieństwo, że ścieżka A da
H_TGP(z) zgodne z BBN/CMB w 5%/0.5% tolerance, wynosi **~55–65%**.

## Konkretne pytania liczbowe

Aby TGP odzyskała fenomenologię BBN+CMB:

1. H_TGP(z=10⁹) ≈ H_GR(z=10⁹) z dokładnością **~5%** (BBN tolerance ⁴He)
2. Położenie pierwszego peaka CMB l_1 = **220.0 ± 0.5**
3. Sound horizon r_s(z_drag) zgodny z Planckiem **147 Mpc**
4. N_eff = **3.046** (relativistic species)

Spójność wszystkich *jednocześnie* w modelu z 1–2 swobodnymi
parametrami (γ, Φ₀ kalibrowane do innych obserwabli) jest **bardzo
nietrywialna**. To nie jest "raczej wyjdzie".

## Pliki dotknięte

- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]
  lin. 148–183 (L01 formal definition, comment block z mapping)
- [[../../core/sek04_stale/sek04_stale.tex]] lin. 27–82 (ax:c jako pole),
  178–238 (prop:c-from-metric), 250–254 (ax:c–ax:G hierarchia),
  508–556 (tabela skalowań)
- [[../../research/op-L01-rho-stress-energy-bridge-2026-05-04/]]
  (formal_definition.md, SM_sector_mapping.md, photon_treatment.md)
- [[POST_ACTION_UPDATE_2026-05-04.md]]
- [[../../core/sek08a]] lin. 119–127 (eq:L-mat-unified)
- **BRAKUJĄCY plik:** `research/op-FRW-radiation-era-varying-c/`
  (cykl jeszcze nie wykonany — najpilniejsze zadanie)

## Problemy jakie rodzi

1. **Brak istniejącego rachunku.** Nie zlokalizowany jawny rachunek
   Φ-EOM w FRW dla z > z_rec z porównaniem do BBN i CMB.

2. **Reżim ekstremalnie wczesny** (z > 10¹⁰): jeśli Φ_then << Φ₀, to
   c(Φ) >> c₀, ℏ >> ℏ₀, G >> G₀. **Czy M9.1'' jest valid w erze
   radiacyjnej?** Otwarte. Może potrzebna oddzielna analiza ψ << 1.

3. **Kinematyka rekombinacji** (z ≈ 1100, T ≈ 0.3 eV):
   - σ_T = (8π/3)(e²/m_e c²)² zawiera c⁻⁴
   - energia Rydberga (Saha) zawiera ℏ⁻²
   - prędkość fotonu = c
   - w TGP wszystkie skalują z Φ jednocześnie

   Czy kombinacje się **znoszą** (zachowując standardową kinematykę),
   czy **łamią** (modyfikując sound horizon i pierwszy peak CMB)?

4. **Konkretne pytania liczbowe** (BBN ⁴He + CMB l_1 + r_s + N_eff)
   — spójność wszystkich *jednocześnie* w 1–2 wolnych parametrach.

5. **Dwustronność diagnozy:**
   - **Zgodność** → fenomenalny wynik (BBN+CMB z 1–2 wolnymi parametrami
     vs ~6 ΛCDM)
   - **Niezgodność** → silny sygnał konieczności rozszerzenia L_mat
     (powrót do ścieżki D, narusza S04)

6. **Open NEEDS w L01** — N1 (quantum trace anomaly EM) i N2 (QCD vacuum)
   pozostają otwarte. W połączeniu z varying constants te otwarte punkty
   są tym, co rozstrzyga, czy fotony są "tylko pasażerami".

## Potencjalne ścieżki domknięcia

### Ścieżka A — explicit Φ-EOM w FRW (cykl `op-FRW-radiation-era-varying-c`)

Wziąć Φ-EOM z `eq:field-eq-reproduced` + M9.1'' metric + ax:c–ax:G +
ρ_matter (mała) + ρ_QCD anomalia (stała). Policzyć **H_TGP(z)** dla
z ∈ [10³, 10¹⁰]. Porównać z H_GR(z).

- **Sukces:** różnica < 5% w BBN, < 0.5% w rekombinacji → TGP odzyskuje
  fenomenologię standardową *bez* ρ_rad jako źródła
- **Niepowodzenie:** jasny sygnał konieczności pivotu L_mat

### Ścieżka B — BBN abundance calc (`op-BBN-TGP`)

Z H_TGP(z) ze ścieżki A wziąć Boltzmann codes (np. AlterBBN port pod
TGP), policzyć abundance ⁴He, D/H, ³He/H, ⁷Li/H. Porównać z PDG:
- Y_p = 0.245 ± 0.003
- D/H = 2.55 ± 0.05 · 10⁻⁵

Tu dochodzą efekty varying ℏ przez energie wiązań nuklearnych — to
nietrywialne, bo ℏ wpływa na Coulomb barrier i kinetyki reakcji.

### Ścieżka C — CMB peaks (`op-CMB-TGP`)

CAMB/CLASS port pod TGP cosmology z varying c, ℏ, G. Position
pierwszego peaka (l ≈ 220) jest fenomenalnie wrażliwy na sound
horizon r_s(z_drag).

### Ścieżka D — alternatywne sprzężenie L_mat (rezerwowa)

L_mat = −(q/Φ_0) φ ρ + L_rad(φ, F_μν) z dodatkowym sprzężeniem
dilatonowym dla pól cechowania. **Narusza ax:metric-coupling** (S04
zamknięty 2026-05-04 przez B9 MICROSCOPE 6/6 PASS) — droga ostatnia.

### Ścieżka E — przyznanie zakresu (najszczersza)

TGP jako "GR w limicie post-recombination", z explicit zaznaczeniem,
że dla z > z_rec teoria nie pretenduje do zastąpienia standardowej
kosmologii. Drastycznie obniża rangę TGP, ale jest *honest*.

## Rekomendowany priorytet

**P1 — krytyczny, otwarte ryzyko.** L01 jest oznaczony jako EXECUTED
2026-05-04 w PRIORITY_MATRIX, ale "executed" oznacza tu *formal
definition*, nie *fenomenologiczne sprawdzenie w erze radiacyjnej*.
Bez ścieżek A+B+C explicit policzonych, **nie wiadomo**, czy TGP
ratuje fenomenologię BBN/CMB przez varying constants, czy nie.

## Powiązanie z istniejącym audytem

- Rozszerza [[README.md]] (L01_rho_operational) poza obecny zakres
  (operacyjna definicja → kosmologia radiation-era z varying c, ℏ, G)
- NEEDS N1, N2 w
  [[../../research/op-L01-rho-stress-energy-bridge-2026-05-04/NEEDS.md]]
  PROMOTED do P1
- Najpilniejszy nowy cykl: `research/op-FRW-radiation-era-varying-c/`

## Cross-references

- [[README.md]] — L01 audit issue (EXECUTED 2026-05-04)
- [[POST_ACTION_UPDATE_2026-05-04.md]] — L01 closure formal
- [[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-1 v2 — recenzja źródłowa
- [[../PRIORITY_MATRIX.md]] — do update z EXT-1 P1
- [[../../research/op-L01-rho-stress-energy-bridge-2026-05-04/NEEDS.md]]
- [[../../core/sek04_stale/sek04_stale.tex]] — ax:c–ax:G axiomy
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]
  L01 formal definition
