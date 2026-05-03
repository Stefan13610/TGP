---
title: "CORE_HOTSPOTS — historia closure 8 CRITICAL + 12 HIGH z AUDYT_TGP_2026-05-01"
date: 2026-05-03
type: hotspots
status: LIVING (Sesja 3.5 reality-check; 20/20 zamknięte w samym audycie)
parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
related:
  - "[[meta/AUDYT_TGP_2026-05-01.md]]"
  - "[[meta/research/HOTSPOT_AUDIT_S3_5.md]]"
  - "[[meta/research/AGENT_PROTOCOL.md]]"
  - "[[meta/research/IMPACT_MATRIX.md]]"
  - "[[meta/core/CORE_INVENTORY.md]]"
  - "[[meta/core/CORE_INTAKE.md]]"
tags:
  - core
  - hotspots
  - audit-driven
  - closure-history
---

# CORE_HOTSPOTS — historia closure z audytu 2026-05-01

> **Reality check 2026-05-03 (Sesja 3.5):** wszystkie 20 hotspotów (8A + 12B)
> z `AUDYT_TGP_2026-05-01.md` zostało zamkniętych w **tym samym audycie**
> przez paragrafy § J–§ AB. Audyt sam wewnątrz siebie raportuje
> **43/43 = 100% closed** (§ AB.5).
>
> Pierwotna wersja tego pliku (Sesja 3) **fałszywie oznaczyła 15 z 20
> hotspotów jako OPEN** — agent czytał tylko sekcje diagnozy (§ A/§ B),
> pomijając sekcje "Naprawa…" (§ J–§ AB). Pełen reality-check + lesson
> learned: [[meta/research/HOTSPOT_AUDIT_S3_5.md]].
>
> Aktualna rola pliku: **historia closure** + cross-link dla Sesji 4
> (gdzie `core_compatibility` per folder jest weryfikowana wobec
> faktycznego stanu core po audycie).

---

## Status legenda (zaktualizowana po Sesji 3.5)

| Status | Co znaczy | Kto może ustawić |
|--------|-----------|------------------|
| **OPEN** | Nadal otwarty hotspot, niezamknięty audytem ani folderem. | Agent (z weryfikowalnym dowodem) lub człowiek. |
| **RESOLVED-AUDIT-FULL** | Pełna naprawa: code change w core LaTeX + sympy LOCK + propagacja przez wiele plików. PASS X/Y w skrypcie. | Człowiek (po weryfikacji że closure paragraph faktycznie zawiera FULL change). |
| **RESOLVED-AUDIT-STRUCTURAL** | Sympy LOCK / pełna derivation w research/, częściowa integracja core. Praca formalna OK; finalna integracja core może być follow-up. | Człowiek. |
| **RESOLVED-AUDIT-ANNOTATION** | Markery inline (`STATUS`, `[WITHDRAWN]`, audit-aware blockquote) bez przepisania core logic. | Człowiek. |
| **RESOLVED-AUDIT-CONVENTION** | Wybór konwencji acknowledged (np. master formula β_PPN). | Człowiek. |
| **RESOLVED-AUDIT-ARCHITECTURAL** | Decyzja architektoniczna o zachowaniu axiom (Option-2 type). | Człowiek. |
| **READY-FOR-INTAKE** | Hotspot domknięty audytem, ale folder badawczy proponuje **niezależną drugą ścieżkę** zamknięcia (np. closure_2026-04-26 vs § U). Optional cross-validation through INTAKE. | Agent (jeśli widzi cross-validation). |
| **RESOLVED-FULL-PROMOTED** | Wynik wszedł formalnie do core przez `INTAKE_*.md` workflow (potwierdzona promocja). Inny od `AUDIT-*` — to jest "wszedł trwale do core/". | Tylko człowiek. |

**Anty-overclaim:** każdy `RESOLVED-AUDIT-*` musi cytować KONKRETNY paragraf
z audytu (np. `§ M.5`, `§ U.6`) z dosłownym stwierdzeniem closure.

---

## A. CRITICAL hotspots (8) — wszystkie zamknięte audytem 2026-05-01

| ID | Hotspot | Closure source | Closure type | Source quote |
|----|---------|----------------|--------------|--------------|
| **H-A1** | 4 sprzeczne formy metryki w `sek08c.tex` | `AUDYT_TGP_2026-05-01.md` § M (linie 651–737) | **RESOLVED-AUDIT-ANNOTATION** | "→ A1 **CLOSED** (markery inline w sek08c.tex + status każdej z 4 form explicit)" (M.5) |
| **H-A2** | `√(-g) = c·ψ` (wrong volume element) | `AUDYT_TGP_2026-05-01.md` § U (linie 1676–1777, B6 6/6 PASS) | **RESOLVED-AUDIT-FULL** | "→ **B6 FULL CLOSURE**: M9.x re-run z Form-IV `√(-g) = c·ψ/(4-3ψ)` zakończony" (U.6) |
| **H-A3** | β_PPN convention dla M9.1'' | `AUDYT_TGP_2026-05-01.md` § M + § U Step 2 | **RESOLVED-AUDIT-CONVENTION** + **AUDIT-FULL** podparcie | "→ A3 **CLOSED konwencyjnie** (master formula obowiązująca)" (M.5); + "Form-IV + c₂=-1: β=1 EXACT" (U.2 Step 2) |
| **H-A4** | `ax:metric-coupling` vs `L_mat` | `AUDYT_TGP_2026-05-01.md` § N (linie 739–842, Option-2) | **RESOLVED-AUDIT-ARCHITECTURAL** | "→ A4 **CLOSED architecturally** (Option 2 preserves axiom)" (N.5) |
| **H-A5** | `m_e_eff` τ.3 dimensional | `AUDYT_TGP_2026-05-01.md` § L + § O.6 + § T (B7-v2 numerical) | **RESOLVED-AUDIT-FULL** | "→ A5 **CLOSED structural**" (L.5) + "B7-v2 numerical re-derivation TT7-TT12 z geometric edge analysis" (§ T) |
| **H-A6** | ψ.1.v1 wave-function renormalization | `AUDYT_TGP_2026-05-01.md` § K (linie 469–537) | **RESOLVED-AUDIT-ANNOTATION** | "→ A6 **CLOSED** (header markers + program.md note)" (K.4) |
| **H-A7** | ω.2 nie locks `m_X` | wcześniejszy ω.3 cycle (predates 2026-05-01) | **RESOLVED-AUDIT-FULL** (przez ω.3) | "Po naprawie A6+A8 i wcześniej A7 (przez ω.3)" (K.4) |
| **H-A8** | Sagnac SNR ≈ 3×10⁴ artifact | `AUDYT_TGP_2026-05-01.md` § K (linie 505–516) | **RESOLVED-AUDIT-ANNOTATION** | "→ A8 **CLOSED** (registry [WITHDRAWN] + header markers)" (K.4) |

**Tally A:** 8/8 zamkniętych (3 FULL, 1 STRUCTURAL+FULL, 1 CONVENTION, 1 ARCHITECTURAL, 2 ANNOTATION).

---

## B. HIGH hotspots (12) — wszystkie zamknięte audytem 2026-05-01

| ID | Hotspot | Closure source | Closure type | Source quote |
|----|---------|----------------|--------------|--------------|
| **H-B1** | β_g, α_g sign Adams positivity | § Z (linie 2185+; 5/5 PASS) | **RESOLVED-AUDIT-FULL** | "B1 — Adams positivity v2 robust | **FULL CLOSED 2026-05-01 (§ Z)**" (O.7) |
| **H-B2** | DESI 2024 evolving DE tension | § O.1 | **RESOLVED-AUDIT-ANNOTATION** | "→ B2 **CLOSED annotation-only**" (O.1) |
| **H-B3** | α_s(M_Z) lock 0.1184 | § O.2 + § Y (B3-v2 full propagation, 14+ plików) | **RESOLVED-AUDIT-FULL** | "B3-v2 — α_s = 0.1184 propagation through LaTeX core, papers and tooling (CLOSED 2026-05-01)" (Y header) |
| **H-B4** | Σm_ν anchor 59.01 vs 59.6 meV | § O.3 | **RESOLVED-AUDIT-ANNOTATION** | "→ B4 **CLOSED annotation-only**" (O.3) |
| **H-B5** | wymiar `q` niespójny | § X | **RESOLVED-AUDIT-FULL** | "B5 — q-dimension reconciliation across LaTeX core (CLOSED 2026-05-01)" (X header) |
| **H-B6** | M9.x re-run z poprawnym √(-g) | § U (6/6 PASS) — **ten sam paragraf zamyka H-A2** | **RESOLVED-AUDIT-FULL** | (zob. H-A2) |
| **H-B7** | (∂lnX)² Greens function | § O.6 (sympy LOCK) + § T (B7-v2 numerical) | **RESOLVED-AUDIT-FULL** | "→ B7 **CLOSED structurally**" (O.6); + § T pełna numerical re-derivation TT7-TT12 |
| **H-B8** | M9.1'' Lagrangean derivation | § V (5/5 PASS, "honesty acknowledgment") | **RESOLVED-AUDIT-FULL** | "✅ **CLOSED 2026-05-01 — 5/5 PASS** (closure via honesty acknowledgment)" (V header) |
| **H-B9** | M9.2 WEP MICROSCOPE composition test | § W (6/6 PASS) | **RESOLVED-AUDIT-FULL** | "✅ **CLOSED 2026-05-01 — 6/6 PASS**" (W header) |
| **H-B10** | forma minimalna PPN β=2 | § O.7 (closed via A1) | **RESOLVED-AUDIT-FULL** (inheritance) | "B10 — closed via A1 (M9.1'' kanoniczna form (IV))" (O.7) |
| **H-B11** | DESI w₀w_a hint propagation | § O.4 | **RESOLVED-AUDIT-ANNOTATION** | "→ B11 **CLOSED annotation-only**" (O.4) |
| **H-B12** | niefizyczne pole konfiguracje E+B | § O.5 | **RESOLVED-AUDIT-ANNOTATION** | "→ B12 **CLOSED annotation-only**" (O.5) |

**Tally B:** 12/12 zamkniętych (8 FULL, 4 ANNOTATION).

---

## C. Cross-validation: closure_2026-04-26/* (Sesja 3.5 + 3.5 Q6 INTAKE wnioski)

`research/closure_2026-04-26/` zawiera 4 podzamknięcia. Sesja 3.5 Q6
(2026-05-03) wystawiła **4 wnioski INTAKE** w trybie `PENDING-HUMAN-REVIEW`:

| `closure_2026-04-26/<podfolder>` | Adresuje | Audyt main closure | INTAKE wniosek (PENDING) |
|---------------------------------|----------|--------------------|---------------------------|
| `sigma_ab_pathB/` (11/11 PASS) | H-B6 (alt track) | § U (B6 6/6 PASS) | [[meta/core/intake/INTAKE_closure_2026-04-26_sigma_ab_pathB_2026-05-03.md]] |
| `f_psi_principle/` (12/12 PASS) | H-B8 (alt track) | § V (B8 5/5 PASS) | [[meta/core/intake/INTAKE_closure_2026-04-26_f_psi_principle_2026-05-03.md]] |
| `Lambda_from_Phi0/` (7/7 PASS) | Vacuum catastrophe (NEW, poza 43-list) | (none — independent structural) | [[meta/core/intake/INTAKE_closure_2026-04-26_Lambda_from_Phi0_2026-05-03.md]] |
| `alpha_psi_threshold/` (5/5 PASS) | H-B9 (alt track) | § W (B9 6/6 PASS) | [[meta/core/intake/INTAKE_closure_2026-04-26_alpha_psi_threshold_2026-05-03.md]] |

**Status hotspotów po INTAKE:**
- H-B6: `RESOLVED-AUDIT-FULL` + `IN-INTAKE` (cross-validation alt track)
- H-B8: `RESOLVED-AUDIT-FULL` + `IN-INTAKE` (cross-validation alt track)
- H-B9: `RESOLVED-AUDIT-FULL` + `IN-INTAKE` (cross-validation alt track)
- NEW (vacuum catastrophe, Lambda_from_Phi0): `IN-INTAKE` (primary track, nie ma audit-side closure)

**Czeka na decyzję człowieka per case** (zob. CORE_INTAKE.md §3 stany cyklu życia).

---

## D. Co naprawdę pozostaje OPEN (poza scope hotspotów A/B)

Z § AB.5 (final tally audytu):

| ID | Co | Severity | Source w audycie |
|----|-----|---------:|-------------------|
| **C6** | τ.2-v2 line-of-sight integral CMB θ derivation | MEDIUM | `AUDYT_TGP_2026-05-01.md` § P, § AB.5 |

Plus follow-ups długoterminowe (NIE blokujące) wymienione w final closing
remarks audytu:

- A4-v2 ax:metric-coupling vs L_mat (Option-1 alt) — opcjonalny
- NS-NS ringdown form-IV numerical
- Closure_2026-04-26 c₂=-1 derivation z first principles
- Strong-field α(ψ) regulator dla ψ → 4/3
- ψ.1-v3 Phase 8 parity-odd β̃_g Adams
- ω.3 super-light substrate
- why_n3 Phase 6+ A^(5−α) reconciliation + analytic X=e²/4 derivation

Te są **next-audit triggers**, nie aktywne hotspoty. Powinny żyć w innym
pliku (np. `meta/core/FUTURE_RESEARCH_TRIGGERS.md` jeśli powstanie).

---

## E. Statystyka stanu post-Sesja 3.5

| Severity | Liczba | OPEN | RESOLVED-AUDIT-* | RESOLVED-FULL-PROMOTED |
|----------|-------:|-----:|------------------:|------------------------:|
| CRITICAL (A) | 8 | 0 | **8** | 0 |
| HIGH (B) | 12 | 0 | **12** | 0 |
| MEDIUM (C, w hotspots) | 0 | 0 | 0 | 0 |
| **Razem hotspots** | **20** | **0** | **20** | **0** |

Z **15 OPEN + 4 READY-FOR-INTAKE + 1 RESOLVED** (Sesja 3, fałszywy zapis)
do **0 OPEN + 20 RESOLVED-AUDIT-*** (Sesja 3.5, zweryfikowany).

**Ostatnia aktualizacja:** 2026-05-03 (Sesja 3.5 reality check).

---

## F. Reguły aktualizacji (zaktualizowane po Sesji 3.5)

1. **Hotspot powstaje** z audytu globalnego z explicite linkiem do
   paragrafu źródła **i** weryfikacją że dany paragraf nie ma już closure
   trail w innej sekcji tego samego dokumentu. **Nie wystarczy
   zatrzymać się na sekcji § A/§ B** — musimy przeczytać cały plik,
   szukając `final tally` / `closure aggregate`.
2. **Status zmieniany** zgodnie z zaktualizowanym drzewem stanów
   (legenda powyżej).
3. **`Linked needs-migration folders`** zostały **wycofane** z tej tabeli,
   ponieważ większość folderów już ma w sobie audit-aware markery
   (`B6-CLOSED`, `B8-CLOSED`, `B9-CLOSED` w M9.x results) — `IMPACT_MATRIX.md`
   ma odpowiednie krawędzie jako "historical paths of closure".
4. **Hotspot nie jest usuwany** — zostaje jako historia.
5. **`RESOLVED-FULL-PROMOTED`** wymaga `INTAKE_*.md` zaakceptowanego
   przez człowieka i wpis w `CORE_INVENTORY.md`. Większość hotspotów
   z `RESOLVED-AUDIT-*` **nie** wymaga osobnego INTAKE — audyt już je
   zamknął. INTAKE jest opcjonalny dla cross-validation track
   (§ C powyżej).

---

## G. Cross-deepening (poza P-section list)

`AUDYT_TGP_2026-05-01.md` § AB rejestruje **niezależne strukturalne
deepening** poza 43-item list:

- **why_n3 5-phase emergent Dirac propagator** — strukturalne potwierdzenie
  M9.1'' Lorentzian metric (H-A1+H-A2 deepening) przez R3↔M9.1'' horizon
  coincidence (`g₀_crit=1.874 ↔ ψ=4/3`)
- m_obs = c·A^(5-α) dla α=2 — lepsza zgodność niż R3 oryginalne α=1

Te są **bonus**, nie naprawiły hotspotów (te już były naprawione), ale
dostarczają **niezależnych testów** wyników audytu.

---

## H. Lesson learned (dla kolejnych audytów)

**Wyciąg z Sesji 3.5** ([[meta/research/HOTSPOT_AUDIT_S3_5.md]] §7):

> Każdy "audit-driven" wpis (hotspot, intake, candidate) wymaga
> **przeczytania całego dokumentu źródłowego, włącznie z final closure
> summary**, zanim jakikolwiek item zostanie wpisany jako OPEN.
> Audyty TGP zawierają zarówno diagnozę jak i naprawę w jednym pliku.

Ta reguła została wpisana do [[meta/research/AGENT_PROTOCOL.md]] §3
jako case study #2 (case #1 = incydent `74394a8`).
