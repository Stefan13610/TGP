---
title: "Inwentaryzacja stanu rdzenia + opcje cleanup (2026-05-04)"
date: 2026-05-04
parent: "[[README.md]]"
type: cleanup-inventory
tgp_owner: audyt
tags:
  - cleanup
  - inventory
  - core-LaTeX
  - ADDENDUM-strategy
  - literal-removal
related:
  - "[[SESSION_REPORT_2026-05-04.md]]"
  - "[[../research/op-g0-r3-from-canonical-projection]]"
  - "[[../research/op-L01-rho-stress-energy-bridge-2026-05-04]]"
  - "[[../research/op-L04-ODE-canonicalization-2026-05-04]]"
---

# Inwentaryzacja stanu rdzenia + opcje cleanup

## 1. Cel inwentaryzacji

Po sesji 2026-05-04, w której zaadresowano **9 audit issues** (S01-S06,
L01, L02, L04 + częściowo L05) przez 6 cykli + 3 NON-BREAKING addytywne
edycje, **rdzeń LaTeX nadal stosuje strategię ADDENDUM** (G.0 świadomie
wybrała: nie literal removal, tylko addytywne v2.0 warstwy).

To inwentaryzuje co konkretnie pozostaje "muzeum pivotów" i przedstawia
trzy opcje cleanup z pros/cons.

## 2. Stan rdzenia per plik

### 2.1 `core/sek08a_akcja_zunifikowana.tex` (1021 linii)

**Struktura dwuwarstwowa (G.0 strategia, 2026-05-02):**

| Linie | Warstwa | Treść | Status |
|-------|---------|-------|--------|
| 1-805 | **v1.x** | V_orig=(β/3)φ³−(γ/4)φ⁴, √(-g)=c·φ, q·c²/Φ_0=2πG_0 | DEPRECATED ale aktywna |
| 806-1021 | **v2.0** | V_M911, √(-g)=c·ψ/(4-3ψ), q·c²/Φ_0=(4/5)πG_0 | KANONICZNA G.0-closed |

**Edycje sesji 2026-05-04 (NON-BREAKING):**

- L01 formal `ρ ≡ -T^μ_μ/c_0²` definition w komentarzu po `eq:L-mat-unified`
  (lin. 102-129+30 nowych)

**Co pozostaje "muzeum pivotów" w body:**

- `eq:V-selfinterference` (lin. 83-93): V_orig (β/3)φ³−(γ/4)φ⁴ jest fizycznie
  nadal w body — z footnote o canonical, ale equation pozostaje
- `eq:sqrt-g-eff` (lin. 141-145): √(-g)=c·φ (forma I) jest fizycznie w body
- `prop:kappa-corrected` (somewhere lin. 386+): κ=3/(4Φ_0) z falsified √(-g)
- `eq:newton-limit` (somewhere): q·c²/Φ_0=2πG_0 (vs v2.0: (4/5)πG_0)

**Cross-references HIGH-impact:** ~30 plików referuje deprecated
propositions. Literal removal wymaga grep-and-replace (P33 audit).

### 2.2 `core/sek08c_metryka_z_substratu.tex` (592 linii)

**Cztery formy metryki współistnieją:**

| Forma | Linie | g_tt | Status post-G.0 |
|-------|-------|------|------------------|
| (I) `eq:metric-full-derived` | ~154 (BOXED) | -c²/ψ | **FALSIFIED 2026-04-25** (β_PPN=4 vs Mercury 3·10⁴σ) |
| (II) Eksponencjalna | ~208-211 | -c²·e^{-2U} | **OBSOLETE / PRZEJŚCIOWA** |
| (III) `thm:antipodal-uniqueness` | ~345-371 | f=ψ⁻¹/², h=ψ¹/² | **NIESPÓJNA WEWNĘTRZNIE** (β_metric=3) |
| (IV) M9.1'' hiperboliczna | tylko nagłówek (lin. 87-104, komentarz) | -c²(4-3ψ)/ψ | **KANONICZNA wg FOUNDATIONS post-G.0** |

**Edycje G.0 (2026-05-02):**

- Preamble G.0 CLOSURE block (lin. 7-49): wyjaśnia że A1+A2+A3 są
  CLOSED-RESOLVED via G.0, że formy I/II/III są deprecated
- 3 inline `STATUS CLOSED` markery przy każdej z form (I), (II), (III)
- Tabela porównawcza M9.1'' (gdzieś w body)

**Co pozostaje "muzeum pivotów" w body:**

- Cała sekcja derywacji formy (I) z budżetu informacyjnego (lin. ~131-200)
- Linearyzacja β=2 dla potęgowej (lin. 183-191) — sprzeczne z M9.1''
- Forma (II) eksponencjalna (lin. ~208-211) bez niezależnej derywacji
- `thm:antipodal-uniqueness` (lin. 345-371) — twierdzenie *jedyności*
  wewnętrznie niespójne
- Linie ~131-590: ~460 linii z deprecated formulami

**Forma (IV) M9.1'' kanoniczna NIE istnieje fizycznie w sek08c body** —
tylko w preamble komentarzu i `TGP_FOUNDATIONS.md:64-69`.

### 2.3 `core/sek08_formalizm.tex` (11826 linii)

**Edycje G.0 (2026-05-02):**

- Intro G.0 closure block (~30 lin)
- 2 inline annotations (FRW prop, footnote)
- `rem:materia-hierarchia` (lin. 9199-9282) ZAKTUALIZOWANY 2026-05-01
  z Phase 2 universal mass formula reference + R5 ↔ Phase 2 bridge

**ax:metric-coupling (lin. 11257-11297):**

- audit A4 note z 2026-05-01 explicit zaktualizowany Option-2 decision
- Fully synced z L01 cycle

**`thm:D-uniqueness` (lin. 956-1048):** czysty, formal proof α=2.

**Co pozostaje:**

- Stare hipotezy/twierdzenia dotyczące pre-G.0 metryki — minor,
  najczęściej z explicit annotations
- Generic clean-up okazjonalnych pre-pivot remarks

### 2.4 `TGP_FOUNDATIONS.md` (272 linii)

**Status pre-update:**

- Data: 2026-04-25
- §3 hierarchia formalizmu już wspomina M9.1'' (lin. 39-40)
- §3 warstwa 0 już ma σ_ab composite (lin. 37) ← post Path B PRIMARY
- §7 M3-M8 archive — używa β/γ bez subskryptu (potencjalna kolizja
  semantyczna z (β/γ)_GL)

**Co potrzebne:**

- Update §7 z explicit reference do `app:A-beta-gamma-distinction` z L02
- Update §3 z reference do Phase 2 universal mass formula (m_obs vs M_full)
- Update §4 (warstwa 3b) z reference do L01 formal `ρ ≡ -T^μ_μ/c_0²`
- Update §10 (canonical sources) z reference do L01/L04 cykli

### 2.5 `README.md` (319 linii)

**Status:**

- `Closure 2026-04-26 update` jest aktualne (lin. 92-103)
- α_s lock 0.1184 (lin. 95) ← B3-v1 partial
- m_H 125.31 (lin. 97) ← C9 PDG-updated
- Σm_ν 59.6 (lin. 100) ← **NIE zaktualizowane** do anchor B4 59.01
- "Kinetic coupling α = 2: algebraic theorem, not a fit" (lin. 102) ← OK

**Co potrzebne:**

- Update Σm_ν → 59.01 (B4-anchored)
- Dodać reference do L04 RESOLVED (mass formula `m_obs = c·A²·g₀^[e²(1−α/4)]`)
- Dodać reference do L01 EXECUTED (formal ρ definition)
- Dodać reference do L02 EXECUTED (β/γ distinction)
- Update sekcja "Predictions" z explicit mass formula?

### 2.6 `INDEX.md` (371 linii)

**Status:** master ledger 856 cumulative (z forward-patch 72 contested).

**Co potrzebne:**

- Dodać entry dla L01, L02, L04 cykli (3 nowych research/op-* foldery)
- Update master verification ledger (jeśli L01/L02/L04 dodają wpisy)
- Sync z S04 N1+N5 closure via L01

## 3. Trzy opcje cleanup

### Opcja A: Status quo (ADDENDUM strategy continued)

**Akcja:** brak zmian fizycznych w body. Tylko meta-updates (FOUNDATIONS,
README, INDEX) z reference do nowych cykli.

**Pros:**

- ✓ 100% bezpieczne — zero risk cross-ref breakage
- ✓ Compile clean już osiągnięty (537 stron pdflatex)
- ✓ G.0 strategia explicit zalecała (Phase 4 § 7.4): „brak ryzyka
  rozsypania innych derivacji"
- ✓ Traceability: czytelnik widzi ŁAŃCUCH derivacji v1.x → v2.0
- ✓ ~1-2h pracy

**Cons:**

- ✗ "Muzeum pivotów" pozostaje w body — sek08c ma 4 formy metryki,
  sek08a ma deprecated propositions razem z v2.0
- ✗ Zewnętrzny referee czyta i widzi sprzeczne równania, mimo że są
  z explicit annotations
- ✗ Manuskrypt nie jest "publication-ready quality"

**Estymata:** 1-2h (głównie meta-updates).

### Opcja B: Hybrydowe (literal cleanup body + historical appendix)

**Akcja:** Usunąć deprecated equations/propositions z body, zachować je
w dedicated `dodatek_pivot_history.tex` z explicit FALSIFIED/OBSOLETE
markers.

**Pros:**

- ✓ Body czysty — czytelnik widzi tylko kanoniczne formuły
- ✓ Historical traceability zachowana w appendix
- ✓ Cross-references zachowane (deprecated → appendix)
- ✓ Manuskrypt staje się publication-ready
- ✓ Ryzyko cross-ref breakage minimalne (tylko relocation, nie deletion)

**Cons:**

- ✗ Wymaga ostrożnej refaktoryzacji ~30 plików HIGH-impact (P33 audit)
- ✗ pdflatex compile clean check po każdej zmianie
- ✗ ~1-2 tygodni pracy

**Konkretne kroki:**

1. **Krok B.1**: Utworzyć `dodatek_pivot_history.tex` z 3 sekcjami:
   - `\section{Pre-G.0 sek08a v1.x}` (V_orig, √(-g)=c·φ, κ z M9.1)
   - `\section{Pre-G.0 sek08c forms (I), (II), (III)}` (FALSIFIED)
   - `\section{Pre-Phase 2 mass formulas LP-4 / R5 K²}` (specjalne case α=1)
2. **Krok B.2**: Refactor sek08a body → tylko v2.0 (V_M911, √(-g)=c·ψ/(4-3ψ),
   etc.) z reference do `dodatek_pivot_history` przy każdej deprecated
3. **Krok B.3**: Refactor sek08c body → tylko M9.1'' kanoniczna z explicit
   derivation z budżetu (jak w preamble), reszta → appendix
4. **Krok B.4**: Cross-reference audit — wszystkie ~30 P33 plików sync
5. **Krok B.5**: pdflatex compile check (powinno być 537+ stron clean)
6. **Krok B.6**: Update FOUNDATIONS, README, INDEX

**Estymata:** 1-2 tygodni systematic work.

### Opcja C: Pełen literal cleanup (deprecated → trash)

**Akcja:** Całkowicie usunąć deprecated propositions z core. Bez appendix.

**Pros:**

- ✓ Najczystszy możliwy stan rdzenia
- ✓ Manuskrypt minimalna objętość

**Cons:**

- ✗ **WYSOKIE RYZYKO** cross-ref breakage (~30 HIGH-impact files)
- ✗ Utrata traceability historycznej (pivot rationale niewidoczna)
- ✗ Każda re-derivacja wymaga ponownego sprawdzenia konsystencji
- ✗ ~2-3 tygodni pracy + intensive testing
- ✗ G.0 strategy explicit ostrzega przed tym

**Estymata:** 2-3 tygodni + ryzyko cofnięcia.

## 4. Rekomendacja

**Opcja B (hybrydowe)** — *jeśli są zasoby na 1-2 tygodnie*. Daje:

- Publication-ready quality manuskryptu
- Traceability zachowana w appendix
- Cross-references intact (deprecated → appendix, nie usunięte)
- Compile clean

**Albo Opcja A natychmiast + Opcja B później** — sesja 2026-05-04 może
wykonać tylko meta-updates (FOUNDATIONS, README, INDEX z L01/L02/L04
references), a body cleanup zostawić dla dedicated cyklu cleanup.

**NIE rekomenduję Opcji C** — zbyt ryzykowne, G.0 strategy explicit
ostrzega.

## 5. Plan dla sesji 2026-05-04 — wybór strategiczny

### 5.A Wariant minimalny (1-2h, BEZPIECZNY)

Wykonać tylko **meta-updates**:

1. **TGP_FOUNDATIONS.md** §3, §4, §7, §10 — references do L01/L02/L04
2. **README.md** — Σm_ν fix do 59.01, reference do L01/L02/L04 EXECUTED
3. **INDEX.md** — entry dla 3 nowych cykli (op-L01, op-L02, op-L04)

To synchronizuje meta-poziom z post-sesja 2026-05-04 stanem, **bez
ryzyka** cross-ref breakage.

### 5.B Wariant pełen Opcja B (1-2 tygodnie, REKOMENDOWANY)

Po wariancie minimalnym, wykonać literal cleanup body:

4. Utworzyć `dodatek_pivot_history.tex`
5. Refactor sek08a body (deprecated → appendix)
6. Refactor sek08c body (forms I/II/III → appendix, body M9.1'')
7. Cross-ref audit (~30 P33 plików)
8. pdflatex compile clean check

**Decyzja autora wymagana.**

## 6. Co konkretnie sesja 2026-05-04 może zrobić **teraz**

Bezpieczna w jednej sesji (~1-2h):

- ✓ TGP_FOUNDATIONS.md sync z L01/L02/L04
- ✓ README.md Σm_ν fix + L01/L02/L04 references
- ✓ INDEX.md nowe entries

**Niebezpieczne (wymagają dedicated cyklu):**

- ✗ sek08a body refactor
- ✗ sek08c body refactor
- ✗ Cross-ref audit ~30 plików

## Cross-references

- [[SESSION_REPORT_2026-05-04.md]] — kontekst sesji
- [[../research/op-g0-r3-from-canonical-projection/README.md]] § 7.4 (G.0 ADDENDUM strategy explicit)
- [[../research/op-L01-rho-stress-energy-bridge-2026-05-04]]
- [[../research/op-L02-...]] (L02 was inline edit, not cycle)
- [[../research/op-L04-ODE-canonicalization-2026-05-04]]
- [[PRIORITY_MATRIX.md]]
