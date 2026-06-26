---
title: "Phase 1 — F-FM-Q4: frontier definition — Q4 RESOLVED_STRUCTURAL (identification); property ledger EXACT; stability boundary Φ₀/√3 (NEW)"
type: phase_result
status: PHASE1_COMPLETE
phase: 1
cycle: op-frontier-microphysics-2026-06-11
created_date: 2026-06-11
authorization: "User 2026-06-11: 'zająć się cyklem op-frontier-microphysics' (activation auth covers Phase 0 + Phase 1, house precedent: FCR sesja #15)"
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
sympy: "8/8 PASS; 0 hardcoded T_pass"
falsifier_resolved: "F-FM-Q4 = RESOLVED_STRUCTURAL (Q4-C identification; conditional on LOCKED concept ontology — risk R-FM-5); F-FM-V NOT touched (Phase 2); F-FM-D unchanged"
anti_lakatos_lock: PRESERVED
---

# Phase 1 — czym jest frontier? (concept §10.6 Q4)

## §0 — Verdict at a glance

**F-FM-Q4 = RESOLVED_STRUCTURAL: Q4-C (IDENTIFICATION).** Pytanie Q4 („granica w Φ-space
**czy** fizyczna granica wszechświata?") jest przy ontologii TGP **źle postawioną dychotomią**:
przestrzeń jest generowana strukturą Φ (pozycja C, concept §1.1 LOCKED), więc **brzeg
przestrzeni generowanej i warstwa przejściowa pola Φ to jeden i ten sam obiekt**. Frontier
otrzymuje pierwszą precyzyjną definicję matematyczną + ledger własności EXACT.

## §1 — Wykluczenia (kryteria KQ1-KQ4, value-blind, pre-LOCKED Phase 0 §1.3)

| Kandydat | Kryterium | Werdykt |
|---|---|---|
| **Q4-A** (brzeg w pre-istniejącej przestrzeni) | KQ1: wymaga background manifold — wprost sprzeczne z deklaracją ontologiczną pozycji C (przestrzeń NIE jest tłem; concept §1.1 LOCKED) | **EXCLUDED** (semantyczne, z LOCKED źródła) |
| **Q4-B** (granica tylko w Φ-space, bez locus przestrzennego) | KQ2: γ-3 LOCKED daje skończony promień R = ct w każdej chwili; KQ3: EQ-2 ⟨Φ⟩_frontier(x,t) jest funkcją położenia, a kreacja jest boundary-localized (FCR A-i LOCKED-claim); sympy: gradient profilu ściany zlokalizowany (max na locus, → 0 poza) | **EXCLUDED** (FP1) |
| **Q4-C** (identyfikacja) | KQ1 ✓ (zero tła), KQ2 ✓ (locus = R(t) = ct), KQ3 ✓ (warstwa przejściowa = nośnik kreacji), KQ4 ✓ („zewnętrze" wchodzi wyłącznie jako wartość referencyjna Φ → 0 typu E1 — żadnego przedłużania metryki poza R) | **UNIQUE SURVIVOR** ⭐ |

## §2 — Definicja matematyczna frontu (D-Q4)

> **Frontier(t)** := warstwa przejściowa pola ⟨Φ⟩ interpolująca |Φ|: Φ₀ (bulk, E2-saturacja)
> → 0 (referencja E1-like), zlokalizowana na sferze **R(t) = ct** (γ-3 LOCKED), o szerokości
> **δ = 2/m_Φ = √2/(√λ·Φ₀)**, gdzie m_Φ = √(V″(Φ₀)) = √(2λ)·Φ₀.
> Locus tej warstwy **jest** brzegiem przestrzeni generowanej — czytania „Φ-space boundary"
> i „physical edge" to dwa opisy jednego obiektu.

Profil referencyjny: half-kink limitu zdegenerowanego (deklarowane przybliżenie thin-wall,
Coleman-style; Phase 0 §8(e)); thin-wall validity: δ/R(t) → 0 (FP5).

## §3 — Ledger własności (EXACT, symbolic; sympy 8/8)

| Własność | Wartość | FP |
|---|---|---|
| Ciśnienie napędowe | ΔV = V(0) − V(Φ₀) = **λΦ₀⁴/4 > 0** — kreacja na froncie energetycznie „z górki"; teza koncepcyjna §2.2 („gradient ⇒ driving force") zweryfikowana **energetycznie** | FP2 |
| Szerokość ściany | δ = 2/m_Φ EXACT (kink spełnia statyczne EOM — residual 0) | FP3 |
| Napięcie powierzchniowe | σ = ∫₀^{Φ₀}√(2V)dΦ = **(2/3)√(λ/2)·Φ₀³** EXACT | FP4 |
| Dynamika ściany | F = ΔV − 2σ/R > 0 dla R > R_c = 2σ/ΔV = (8√2/3)/(√λΦ₀); wypadkowe parcie na zewnątrz utrzymane ⇒ γv rośnie ⇒ **v → c**: front asymptotycznie zerowy (null) — **CONSISTENT z γ-3** (γ-3 pozostaje LOCKED źródłem R = ct; tu tylko check spójności) | FP5 |
| **Granica stabilności (NEW)** | m_eff²(Φ) = V″(Φ) = λ(3Φ²−Φ₀²) ⇒ **stabilna masywna materia istnieje tylko tam, gdzie \|Φ\| > Φ₀/√3** — ściśle WEWNĄTRZ ściany (głębokość x* = δ·atanh(1/√3) ≈ 0.659δ od centrum) | FP6 |

## §4 — Konsekwencja strukturalna (input do Phase 2, NIE rozstrzygnięcie)

FP6 daje pierwszy twardy mikrofizyczny fakt o kreacji: **materia nie może ustabilizować się
NA locusie frontu** (tam m_eff² < 0) — wchodzi do bulku przez **wewnętrzną krawędź ściany**
(|Φ| = Φ₀/√3). To zawęża semantykę v_c (Phase 0 §1.2: prędkość przy WEJŚCIU do source-free
bulk) do dobrze określonego zdarzenia: przejście przez wewnętrzną krawędź warstwy.
**Tiebreaker B-k3 vs B-k4 NIE jest tu rozstrzygany** — mechanizmy M1/M2/M3 i kryteria K1-K4
czekają na Phase 2 (forbidden move #3 respektowany: żadne kryterium nie zostało użyte
przedwcześnie).

## §5 — Walidacja pre-derivacji (§3.6.9)

p₊(ε) z równania indicjalnego p² + p/3 − ε = 0: p₊(2/3) = 2/3 EXACT (EdS), p₊(3/2) =
(√55−1)/6 EXACT; log₁₀G = {2.0252, 3.2485}, rel dev ≤ 4×10⁻⁵ względem wartości oczekiwanych
Phase 0 §8(b) (standard ±5%). Δ_bulk(ε) = |3ε−2|/4: **jedyne zero ε = 2/3**; Δ_bulk(3/2) = 5/8.

## §6 — Anti-Lakatos (Phase 1): COMPLIANT ✓

0/8 hardcoded; zbiór kandydatów Q4 CLOSED pre-declared, wykluczenia kryteriami semantycznymi
z LOCKED źródeł (nie wartością wzrostu; G_obs nieobecne — FP8 free-symbols audit); werdykt
RESOLVED_STRUCTURAL **warunkowy względem ontologii concept-paper** (R-FM-5 ujawnione, klasa
rygoru §3.6.12); przybliżenie thin-wall deklarowane per-use (§10.1 mitigation); 0 predecessor
verdicts modified; tiebreaker nietknięty; 0 nowych stałych fundamentalnych (λ, Φ₀ wyłącznie
symboliczne).
