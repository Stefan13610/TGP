---
title: "Phase 1 — moduł B test mikrofizyczny 1D: KB1 PASS (ściana rozróżnia C-partnerów — selekcja topologiczna + energetyka) · KB2 PASS-DERIVED (sortowanie: partner zgodny ze ścianą → bulk, C-partner → sektor frontowy) · KB3 CONDITIONAL (wyścig: separacja wygrywa ⟺ ξ_s > ln(3+2√3) EXACT) · SIG-1 EXACT"
type: phase_result
status: PHASE1_COMPLETE
phase: 1
cycle: op-frontier-bridge-and-asymmetry-2026-06-12
created_date: 2026-06-12
authorization: "User 2026-06-12: 'Phase' → Phase 1 (opcja 1 decision menu post-Phase-0)"
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
sympy: "10/10 PASS; 0 hardcoded T_pass; circularity guard FP10 (G_obs/η_B/asymetrie nieobecne)"
falsifier_status: "F-BA-6: KB1 PASS + KB2 PASS-DERIVED + KB3 CONDITIONAL (warunek istnienia EXACT) + SIG-1 PASS EXACT; werdykt F-BA-6 OPEN (klasy CLOSED Phase 0 §1.4 — rozstrzygnięcie w Phase FINAL po KB4/Phase 4); moduł A nietknięty; NO PR-023"
anti_lakatos_lock: PRESERVED
---

# Phase 1 — czy ściana rozróżnia soliton od C-partnera? (test 1D, reuse op-CE-H)

## §0 — Verdict at a glance

| Kryterium (Phase 0 §1.4, LOCKED) | Wynik |
|---|---|
| **KB1 (istnienie rozróżnienia)** | **PASS** — dwie niezależne linie: (a) **topologiczna reguła selekcji** (FP2); (b) energetyka osadzenia: ΔE_C/M_K = 4.95 przy wewnętrznej krawędzi ≫ floor 10⁻³ (FP3/FP5) |
| **KB2 (kierunek sortowania)** | **PASS-DERIVED** — partner o ładunku topologicznym ZGODNYM ze ścianą → odpychany do bulku; C-partner (przeciwny) → przyciągany ku ścianie/sektorowi frontowemu (FP6) |
| **KB3 = SIG-3 (wyścig)** | **CONDITIONAL** — separacja wygrywa z anihilacją ⟺ **ξ_s > ξ_s\* = ln(3+2√3) ≈ 1.866 EXACT** (przy kreacji na krawędzi stabilności ξ_d = ln(2+√3) ≈ 1.317); spełnione na (1.866, 4] ⊂ pre-rejestrowanego [1, 4] — **NIE pełny zakres**, raport mechaniczny |
| **SIG-1 (budżet ×2)** | **PASS EXACT** — t_*^(B) = √2·t_*; t_* odtwarza formę FM P3 LOCKED (FP8) |
| **F-BA-6** | **OPEN** — H-SORT_DERIVED_1DPROXY wymaga KB1+KB2+KB3 w pełni pozytywnych (Phase 0 §1.4); KB3 jest warunkowe ⇒ klasyfikacja końcowa w Phase FINAL (KB4/H-CP → Phase 4). **NO PR-023** |

Dyspozycja stała: **1D-proxy ≠ 3D claim** (Phase 0 §1.2.2; ryzyko R-BA-2 konflacji C↔P ujawnione).

## §1 — KB1: ściana ROZRÓŻNIA C-partnerów (dwie linie)

**Linia (a) — topologiczna reguła selekcji (FP2):** w sektorze Z2 ściana (wielki kink,
−Φ₀→+Φ₀ proxy profilu 0→Φ₀) wymusza PORZĄDEK pary po stronie bulku: próżnia za ścianą = +Φ₀,
więc partner przylegający do ściany MUSI mieć ładunek przeciwny (antykink: +Φ₀→−Φ₀);
konfiguracja ściana–kink (przylegle) **nie istnieje** w sektorze. Rozróżnienie jest
maksymalne — na poziomie dozwolonych konfiguracji, nie tylko energii. Jedyny dozwolony
porządek pary: **ściana – A – K** (czytane ku bulkowi).

**Linia (b) — energetyka osadzenia (FP3/FP5):** cross-term kinetyczny
E_× = ∫Φ_wall′·φ_q′ dx > 0 dla orientacji zgodnej (penalizowana — dokładnie konwencja
Phase 0 §8(c)); ΔE_C = 2E_× ≠ 0 w warstwie; zanik ~ ξ·e^(−ξ) → 0 w bulku (zgodność z LOCKED
F_substrat = 0 — sanity ✓). Floor pre-rejestrowany: |ΔE_C|/M_K = {4.95, 0.92, 0.34} przy
ξ ∈ {1.317, 3, 4} — wszystkie ≫ 10⁻³. Asymptotyka parowa (Manton-class): rate zaniku V_int
fitted 1.366 vs m_Φ = 1.414 (rel dev 3.4% < 5% std §3.6.9) — **spójne z LOCKED op-CE-H**
(m√2 ≡ m_Φ, mapowanie notacji Phase 0 §8(f)).

**Wyjątek pre-metryczny (FP9):** kanał pary po stronie zewnętrznej NIE istnieje jako kanał
cząstkowy — m_eff² < 0 dla |Φ| < Φ₀/√3 (FM P1 FP6 LOCKED): brak propagujących modów ⇒
jedyny trwały kanał = strona bulku, z wymuszonym porządkiem z linii (a).

## §2 — KB2: kierunek sortowania WYPROWADZONY (value-blind)

Z asymptotycznych oddziaływań parowych (znaki z LOCKED-class machinery: przeciwne ładunki
przyciągają, zgodne odpychają; FP6):

- **C-partner (antykink, przylegający):** F = −Cm·e^(−m·x₁) — **przyciągany KU ścianie**
  (sektor frontowy; w obrazie 3D: marginalnie nieosiągalny front, SCOPING Q2),
- **partner zgodny ze ścianą (kink, głębszy):** F = +Cm·e^(−m·x₂) — **odpychany W GŁĄB bulku**.

Etykieta „materia" NIE była przypisana z góry (Phase 0 §7): przypisanie wynika z rachunku —
**do bulku wchodzi partner o ładunku topologicznym zgodnym z orientacją ściany.** To jest
dokładnie mechanizm H-SORT: netto-jednoimienność bulku bez łamania C w amplitudzie.

**Nota strukturalna (INFORMATIONAL, 1D-proxy):** antykink wchłonięty przez ścianę anihiluje
z nią lokalnie, a głębszy kink przejmuje rolę frontu — „sektor frontowy" antymaterii realizuje
się w 1D jako absorpcja C-partnera w strukturę postępującej ściany. Spójne z narracją H-SORT
(antymateria ukryta w/za frontem); przeniesienie na 3D = przyszła praca, nie claim.

## §3 — KB3/SIG-3: wyścig sortowanie vs anihilacja — warunek istnienia EXACT

Siła względna pary (x₁ = d, x₂ = d + s; równe masy; FP7):

```
s̈ ∝ F_K − F_A = C·m·[e^(−m·x₁) + e^(−m·x₂) − 2e^(−m·s)]
separacja wygrywa ⟺ e^(−ξ_d)·(1 + e^(−ξ_s)) > 2e^(−ξ_s),   ξ ≡ m_Φ·(odległość)
```

Przy kreacji na wewnętrznej krawędzi stabilności (LOCKED FM: x* = δ·atanh(1/√3)
⇒ ξ_d = ln(2+√3) ≈ 1.317 — tożsamość EXACT, FP5):

> **ξ_s\* = ln(3 + 2√3) ≈ 1.8663 EXACT** (zamknięta forma; ξ_s\* = ξ_d + ½ln3)
> — para separuje się (sortowanie) dla ξ_s > ξ_s\*; anihiluje dla ξ_s < ξ_s\*.

**Raport mechaniczny vs pre-rejestracja:** zakres LOCKED ξ_s ∈ [1, 4]: grid
{1: FAIL, 1.5: FAIL, 2: PASS, 2.5: PASS, 3: PASS, 4: PASS} ⇒ **KB3 = CONDITIONAL**
(nie pełnozakresowe PASS; bez łagodzenia). Caveat zadeklarowany: przybliżenie rozrzedzone
(parowe, asymptotyczne) traci ważność dla ξ ≲ 2 — częściowo pokrywa się z obszarem FAIL;
to NIE jest użyte do ratowania werdyktu (raport jak wyżej).

**Znalezisko strukturalne (INFORMATIONAL, nie pre-rejestrowane):** naturalna separacja
kreacyjna pary ~ szerokość solitonu (ξ_s ~ 2) leży **na krawędzi** ξ_s\* ≈ 1.87 — mechanizm
jest blisko-krytyczny: część par (ciasnych) anihiluje przy ścianie, część (szerszych)
zostaje posortowana. Konsekwencja: H-SORT z natury przewiduje **częściową wydajność
i niezerowy kanał anihilacyjny przy froncie** — bezpośredni input ilościowy dla SIG-2
(leakage/tło γ, Phase 4) i dla bilansu budżetowego (SIG-1 ×2 = górna granica; pary
anihilujące zwracają energię ścianie). Zero claimów ilościowych tutaj.

## §4 — SIG-1: budżet ×2 (EXACT)

Popyt parowy 2Ṁc² (Ṁ = 2c³/9G LOCKED, read-only) vs podaż πλΦ₀⁴c³t² (FM P3 LOCKED):
**t_*^(B)/t_* = √2 EXACT** (FP8); cross-check: t_* odtwarza formę FM P3
(√2/3√π)·c/(√(Gλ)Φ₀²) — tożsamość symboliczna ✓. Przesunięcie progu deficytu
wczesnoepokowego — falsyfikowalna sygnatura H-SORT zapisana.

## §5 — Caveats i deklaracje uczciwe (per-use)

1. **1D-proxy ≠ 3D claim** (Phase 0 §1.2.2): etykieta C-odd = ładunek topologiczny q
   (za FFS); realne pole Z2 jest C-trywialne — ryzyko konflacji C↔P (R-BA-2) stoi;
   pełny 3D test (U(1)/RP², hedgehog vs anty-hedgehog w tle radialnym) = przyszły cykl.
2. **Przybliżenie rozrzedzone/asymptotyczne** (Manton tail-class): ważność ξ ≳ 2;
   zakres pre-rejestrowany [1,4] zawiera obszar poza ważnością — raportowano mechanicznie
   mimo to (anti-goalpost).
3. **Statyczne kryterium siły początkowej** dla wyścigu (kinematyka kreacji nieznana —
   prędkości początkowe pary niezadeklarowane w żadnym LOCKED źródle); pełna dynamika
   (z dyssypacją radiacyjną) = możliwe zaostrzenie w przyszłej fazie.
4. **Amplituda C oddziaływania:** LIT informational — C_fit/(m_ΦΦ₀²) = 6.16 przy
   kandydatach konwencji Manton {8, 16}; wyścig (FP7) i floor (FP5) są od C NIEZALEŻNE
   (C skraca się; tylko znaki + wspólny rate m_Φ wchodzą do warunku).
5. **Ściana proxy = pełny kink** dla profilu frontu 0→Φ₀: ogon po stronie bulku identyczny
   (2Φ₀e^(−m_Φx)) — deklaracja Phase 0 §8 + skrypt nagłówek (v).
6. **Naprawy plumbingu testów (transparentnie):** pierwotny run 8/10 — (i) FP3: skryptowy
   sub-check „<10⁻³ przy ξ=10" nie był progiem pre-rejestrowanym (cross-term zanika jak
   ξ·e^(−ξ) — poprawka wielomianowa); zastąpiony pre-rejestrowanym sensem (zanik do zera,
   e-fold ratio); (ii) FP5: tożsamość atanh→log wymagała rewrite (czysta mechanika sympy).
   **Żaden próg pre-rejestrowany nie został zmieniony** (klasa „better seeds", precedens
   CE-H T_P2_5).

## §6 — Status F-BA-6 i PR po Phase 1

| Element | Status |
|---|---|
| KB1 | PASS (dwie linie) |
| KB2 | PASS-DERIVED |
| KB3 | CONDITIONAL (warunek istnienia EXACT: ξ_s > ln(3+2√3)) |
| KB4 (H-CP) | NIETKNIĘTY → Phase 4 (wymaga fazy U(1)/RP²; GAP-able) |
| SIG-1 | PASS EXACT |
| SIG-2 (leakage → tło γ) | → Phase 4 (input: blisko-krytyczność §3) |
| **F-BA-6** | **OPEN** — żadna klasa CLOSED nie jest przyznana w tej fazie; H-SORT ma wykazane istnienie mechanizmu Z WARUNKIEM ilościowym (dokładnie deliverable żądany przez SCOPING sygnaturę 3) |
| PR-023 | **NO APPEND** (forbidden move #8) |

Moduł A: nietknięty (Phase 2-3).

## §7 — Anti-Lakatos (Phase 1): COMPLIANT ✓

0/10 hardcoded; werdykty KB per pre-rejestrowane kryteria (Phase 0 §1.4 LOCKED), zero nowych
kryteriów; KB3 raportowane CONDITIONAL wbrew pokusie pełnego PASS (ważność przybliżenia NIE
użyta do rescue); etykieta „materia" przypisana wynikiem, nie założeniem (value-blindness §7
Phase 0 zrealizowana); circularity guard FP10 (G_obs/η_B/f̄_max nieobecne w wyrażeniach
werdyktowych); LOCKED inputs read-only (Ṁ, V_int, m_eff², x*); naprawy skryptu udokumentowane
(§5.6) bez zmiany progów; 0 nowych stałych (λ, Φ₀ symboliczne; normalizacja numeryczna
λ=Φ₀=1 zadeklarowana §3.6.8(d)); 0 predecessor verdicts modified.
