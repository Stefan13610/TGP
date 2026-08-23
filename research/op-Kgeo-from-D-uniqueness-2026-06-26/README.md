---
title: "op-Kgeo-from-D-uniqueness — czy K_geo·m_sp²=π·Φ₀² jest derywowalne z poziomu-0 (D-uniqueness + geometria rury), czy nieredukowalnie postulowane?"
date: 2026-06-26
type: research-cycle
folder_status: closed-resolved   # CLOSED 2026-06-26: (C) POSTULATE-CONFIRMED (9/9). Phase 0 LOCK + Phase 1 + FINAL w 1 sesji — patrz Phase_FINAL_close.md
claim_status: "C — STRUCTURAL_VERIFIED (pending-bridge: CG-1/CG-3)"
parent: "[[../../partial_proofs/quark_sector/dodatekX_quark_sector.tex]]"

# ============== KICKOFF CONTRACT (mandatory post-2026-05-10) — DRAFT do LOCK ==============
contract:
  # --- L1: Native (MANDATORY) ---
  L1_native:
    output_observable: "K_geo (bezwymiarowy parametr geometryczny substratu) ORAZ — przez most — α_s(M_Z) (bezwymiarowe, PDG-mierzalne 0,1179±0,0009)"
    measurement_instrument: "α_s(M_Z): world-average PDG (DIS + e⁺e⁻ + lattice + τ-decays). K_geo: pochodna poziomu-0 (brak bezpośredniego instrumentu — testowany przez most do α_s)."
    native_coefs_constrained:
      - "K_geo (czynnik geometryczny rury kolorowej, prop:X-A-from-tube-tension)"
      - "𝒜 = a_Γ/φ (uniwersalna stała konfinementu; status #43 POSTULATE-CONDITIONAL)"
    falsification_rule: >
      "Niech K_geo^(0) = niezależnie wyprowadzona z poziomu-0 (thm:D-uniqueness + całkowanie po
      kącie rury), bez użycia α_s ani mas kwarków (circularity guard). Wtedy:
      (A) DERIVED — jeśli K_geo^(0)·m_sp²/(π·Φ₀²) ∈ [0,95; 1,05] (≤5%), most 𝒜=C_F²α_s² jest
          domknięty od poziomu-0 ⇒ α_s(M_Z) staje się genuine first-principles predykcją;
      (B) REFUTED-BRIDGE — jeśli K_geo^(0)·m_sp²/(π·Φ₀²) ∉ [0,80; 1,25] (>25% w którąkolwiek
          stronę), postulat eq:X-K-msp-hypothesis jest sfalsyfikowany (uruchamia PR-025
          forward falsifier (b)); 𝒜 NIE jest mostem do α_s;
      (C) POSTULATE-CONFIRMED — jeśli wynik w [0,80;0,95)∪(1,05;1,25] LUB K_geo^(0) okazuje się
          nieoznaczalny bez domknięcia CG-1/CG-3 (ex200 α_eff niezbieżny przy dostępnym L),
          status POSTULATE-CONDITIONAL (#43) RATYFIKOWANY — K_geo irreducibly conditional."
    pre_registration_date: "2026-06-26"   # LOCKED immutable @ Phase 0 (przed jakimkolwiek rachunkiem K_geo^(0)); progi 5%/25% value-blind zaplombowane

  # --- L2: Cross-framework reduction (OPTIONAL — last stage) ---
  L2_framework_reduction:
    target_frameworks: ["lattice-QCD string tension σ ↔ Λ_QCD", "Casimir C_F α_s/(π Φ_0) color-tube ansatz"]
    reduction_type: "not-attempted"   # do oceny w Phase FINAL; fail NIE unieważnia L1
    validation_transfer: ""
    failure_disposition: "L1-stands"

  # --- L3: Falsification map (consistency) ---
  L3_falsification_map:
    - { bound: "PDG α_s(M_Z)=0,1179±0,0009", constrains: "K_geo via most 𝒜=C_F²α_s²", window: "0,03σ (pole) — consistency only, NIE primary", status: "pending" }
    - { bound: "ex200 α_eff convergence (Γ→Φ coarse-graining)", constrains: "K_geo niezależność", window: "4/8 PASS @ current L", status: "pending" }

# ============== END KICKOFF CONTRACT ==============

tgp_status:
  level: T3                 # głęboki korzeń strukturalny (substrat → makro)
  kind: derivation          # próba derywacji load-bearing postulatu
  output_type: structural   # K_geo = parametr strukturalny; α_s observable WARUNKOWO (gdy DERIVED)
  core_compatibility: review-only
  may_edit_core: false
  has_needs_file: false
  has_findings_file: false
  exports_findings: false
  open_bridges: ["op-CG34-continuum-closure-2026-06-14", "op-Csigma-coarse-graining-2026-06-14", "op-uv-as-ngfp"]
  depends_on: ["#43 op-A-derivation-from-CG (POSTULATE-CONDITIONAL)", "#41 op-quark-mass-core-g0-rescue-test (𝒜 numeryczne)"]
  impacts: ["partial_proofs/quark_sector/dodatekX prop:X-A-from-tube-tension", "meta/PRE_REGISTERED_FALSIFIERS PR-025 forward (b)", "PREDICTIONS_REGISTRY α_s status", "wspólny korzeń z α=2 (#39) i c₀ (#37)"]
  source_of_status: ["[[../op-A-derivation-from-CG-2026-06-25/Phase_FINAL_close.md]]"]
---

# op-Kgeo-from-D-uniqueness — cykl INICJUJĄCY track UV/CG (most Γ→Φ)

> **STATUS: 🟢 CLOSED-RESOLVED — (C) POSTULATE-CONFIRMED (2026-06-26, sympy 9/9).**
> Phase 0 LOCK + Phase 1 + FINAL w 1 sesji. **Wynik:** R := K_geo⁽⁰⁾·m_sp²/(π·Φ₀²) jest
> **nieoznaczalne value-blind** — żadna z 3 ścieżek nie daje niecyrkularnej liczby K_geo⁽⁰⁾:
> (1) D-uniqueness fixuje FORMĘ `K(φ)=K_geo·φ⁴`+`α=2`, ale K_geo to wolna stała całkowania
> (C3 ją *definiuje*, nie liczy); (2) K_geo absorbowalne przez `Φ̃=√K_geo·φ³/3` (brak
> niezmiennika poziomu-0); (3) geometria rury daje π (z kąta, `σ̂=πA²`, potwierdzone bez α_s),
> ale NIE skalę `K_geo·m_sp²`; CG (ex200 4/8) niezbieżny. ⟹ gałąź (C) reguły (pre-zarejestrowana)
> **ratyfikuje #43 POSTULATE-CONDITIONAL**. α_s pozostaje consistency-check warunkowy. Audit
> trail: `Phase0_LOCK.md` + `Phase1_Kgeo.py`/`.txt` + `Phase_FINAL_close.md`.

## §0 — Cel + native-first contract

[CITE: meta/CYCLE_KICKOFF_TEMPLATE.md §1; meta/PPN_AS_PROJECTION.md §3.1; meta/M9_RESTRUCTURE_NOTE.md §2]

### §0.1 — Native observable target

#43 ustalił, że cały łańcuch 𝒜 = C_F²α_s² (α_s 0,03σ) wisi na **jednym load-bearing
postulacie** — identyfikacji `eq:X-K-msp-hypothesis`:
$$K_{\rm geo}\cdot m_{\rm sp}^2 = \pi\cdot\Phi_0^2.$$
Łańcuch: L1 σ̂=πA² (ansatz) · L2 A_color=C_F α_s/(πΦ_0) (derywowane) · L3 m_sp²=γ
(derywowane) · L4 σ_phys=K_geo·m_sp²·σ̂ (definicyjne) · **L5 = ten postulat** ⟹ 𝒜=C_F²α_s².

**Cel inicjujący:** wyznaczyć **K_geo niezależnie od poziomu-0** (z `thm:D-uniqueness` +
geometrii całkowania po kącie rury kolorowej), **BEZ** wstrzykiwania α_s lub mas kwarków
(circularity guard T-anti-circ). Native output = K_geo (bezwymiarowy); pochodny observable =
α_s(M_Z), który staje się **genuine** predykcją *tylko jeśli* werdykt = DERIVED.

To jest **wspólny mianownik** trzech korzeni warunkowych: 𝒜 (#43, ten cykl), α=2 (#39, NGFP),
c₀ (#37, Ward/UV) — wszystkie czekają na domknięcie UV/CG. Cykl inicjujący atakuje najbardziej
skonkretyzowany z nich (K_geo ma jawne równanie kandydujące).

### §0.2 — Pre-registered falsification rule (🔒 PLOMBED @ Phase 0 LOCK 2026-06-26)

Reguła trójwartościowa **(A) DERIVED / (B) REFUTED-BRIDGE / (C) POSTULATE-CONFIRMED** —
patrz `contract.L1_native.falsification_rule` powyżej. Progi (5% / 25%) **value-blind**,
🔒 **zaplombowane 2026-06-26 PRZED rachunkiem** (immutable `pre_registration_date`).
**Anti-moving-goalposts:** po Phase 1 zmiana progu = HALT-B. Plomba: metryka rozstrzygnięcia
= stosunek `R := K_geo^(0)·m_sp²/(π·Φ₀²)`; (A) `R∈[0,95;1,05]`; (B) `R∉[0,80;1,25]`;
(C) `R∈[0,80;0,95)∪(1,05;1,25]` LUB K_geo^(0) nieoznaczalne bez domknięcia CG-1/CG-3.

**Sprzężenie z PR-025:** wynik (B) REFUTED-BRIDGE = bezpośrednie uruchomienie forward
falsyfikatora (b) z PR-025 (`K_geo·m_sp²≠π·Φ₀² ⟹ 𝒜≠C_F²α_s²`). Wynik (A) DERIVED = α_s
przechodzi z TRADED (#42) do genuine DERIVED (zmienia ledger parameter-counting). Wynik (C) =
ratyfikacja status quo #43 (najbardziej prawdopodobny a priori, given ex200 4/8).

### §0.3 — TGP-native check (mandatory, pre-Phase-1)

- **Q1 (native, nie projection?):** ✅ K_geo to natywny parametr geometryczny substratu, nie
  parametr obcego frameworku. α_s observable jest pochodną, nie primary.
- **Q2 (observable units?):** K_geo bezwymiarowy; α_s(M_Z) bezwymiarowy, PDG-mierzalny.
- **Q3 (instrument?):** PDG world-average α_s; K_geo testowany przez most.
- **Q4 (circularity guard?):** ⚠ KRYTYCZNY — K_geo^(0) MUSI być liczone bez α_s/mas kwarków.
  T-anti-circ obowiązkowy w Phase 1 (zero PDG α_s w derywacji K_geo).
- **Q5 (pre-registration?):** reguła + progi w Phase 0 LOCK, immutable timestamp.
- **Q6 (forced-zero declare-not-verify?):** N/A (brak α_i/ζ_i w zakresie).
- **Q7 (anti-Lakatos two-sided?):** ✅ reguła ma wynik POZYTYWNY (A), NEGATYWNY (B), i
  WARUNKOWY (C) — żaden nie jest faworyzowany.
- **Q8 (nie re-litiguje #43?):** ✅ #43 werdykt IMMUTABLE; ten cykl wyznacza K_geo^(0)
  (czego #43 NIE robił — #43 tylko zlokalizował, że K_geo nie jest niezależnie ustalone).

### §0.4 — Pre-flight methodology read confirmation (mandatory, KICKOFF §2.6)

- [x] Przeczytano `PPN_AS_PROJECTION.md` §3.1 — three-layer L1/L2/L3 (native obserwable
      PRIMARY; framework reduction = consistency map, ostatni etap)
- [x] Przeczytano `TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §1-§4 — anti-BD-drift; ASK-RULE
      triggery A-D; Pattern 2.1 (Φ_eq, D_kin ≠ Klein-Gordon), §3 red flags
- [x] Przeczytano `M9_RESTRUCTURE_NOTE.md` §1.4 + §3 — M9.1'' jako Tier-2 anchor, NIE framework
- [x] Przeczytano `CYCLE_KICKOFF_TEMPLATE.md` §1-§2 — kickoff contract, ordering L1→L2
- [x] Read-lock źródeł (read-only, Phase0_LOCK.md §1): `dodatekX` prop:X-A-from-tube-tension
      (l.954-1048) + eq:X-K-msp-hypothesis (l.975) + warunek (l.1032-1039); `thm:D-uniqueness`
      (sek08 l.962-1000); `status_map` CG-1 (l.1323-1329) / CG-3 (l.1514-1517) + K_geo conformal
      (l.72); ex200 4/8 (l.1517) + ex202 7/8 T6 FAIL (l.1523); #43 Phase_FINAL_close

**Sign-off:** Claudian @ 2026-06-26 (Phase 0 LOCK). Pełny audit trail → `Phase0_LOCK.md`.

## §1 — Phase 0: scope mapping (✅ WYKONANE @ LOCK — patrz Phase0_LOCK.md)

1. **Read-lock źródeł** (read-only): `dodatekX prop:X-A-from-tube-tension` + `eq:X-K-msp-hypothesis`;
   `status_map.tex` CG-1/CG-3 [SZKIC]; `thm:D-uniqueness`; ex200 (α_eff convergence) + ex202.
2. **Pre-flight read confirmation** (§2.6 template): PPN_AS_PROJECTION §3.1, NATIVE_PATTERNS
   §1-§4, M9_RESTRUCTURE §1.4+§3, KICKOFF §1-§2.
3. **Plombowanie reguły** (A/B/C + progi 5%/25%) + immutable timestamp.
4. **Balance gate:** budżet nowych stałych = 0 (K_geo NIE jest nową stałą — ma być POCHODNĄ
   z D-uniqueness; jeśli okaże się wolnym parametrem → to wynik C).

## §2 — Phase 1: native derivation (K_geo z poziomu-0)

- **Ścieżka 1 (geometryczna):** całkowanie po kącie rury kolorowej — czy czynnik π·Φ₀²
  wychodzi z geometrii przekroju rury (Bessel K_0 vs gauss ansatz; #43 L1 zaznaczył, że
  Bessel daje ×0,3–0,6, więc ansatz nieunikalny — to kluczowy test).
- **Ścieżka 2 (D-uniqueness):** czy K_geo wynika z `thm:D-uniqueness` (poziom-0 selekcja
  klasy konforemnej) niezależnie od m_sp.
- **Ścieżka 3 (lattice/CG):** ex200 z większym L — czy α_eff zbiega (obecnie 4/8 PASS).
  Jeśli niezbieżny przy dostępnym L → silny sygnał (C).
- **Sympy/numeryka:** ≥1 niezależny solver, 0 hardcoded, T-anti-circ (zero α_s w wejściu).

## §3 — Phase 2-N: native obserwable computation

- Wyliczyć K_geo^(0) ze ścieżki najmocniejszej; obliczyć stosunek K_geo^(0)·m_sp²/(π·Φ₀²).
- Zastosować regułę A/B/C (value-blind, z plomby Phase 0).
- Jeśli (A): wyliczyć α_s(M_Z) first-principles i porównać z PDG (teraz genuine, nie consistency).

## §FINAL — Optional L2 framework reduction

[OPTIONAL] Mapowanie K_geo na lattice string-tension σ ↔ Λ_QCD (Casimir scaling). Fail NIE
unieważnia L1.

## §FINAL+1 — L3 falsification map + propagacja

- Aktualizacja `dodatekX rem:X-rescope-43` + `PREDICTIONS_REGISTRY` α_s status + `PR-025`
  (rezultat forward (b)) — **user-gated**, per wynik A/B/C.
- Wpis STATE.

## §estimate — Skala

1–2 sesje dla **rozstrzygnięcia A/B/C** (nie pełne domknięcie CG — to wieloletni track).
Wynik C (ratyfikacja POSTULATE-CONDITIONAL) jest a priori najbardziej prawdopodobny (ex200
4/8); wartość cyklu = **przekształcenie „nie wiemy" w precyzyjną mapę obstrukcji** (analog
#37 c₀, #39 α=2), ewentualnie zaskakujące (A) lub (B).

## Cross-references
- [[../op-A-derivation-from-CG-2026-06-25/Phase_FINAL_close.md]] (#43 — źródło luki)
- [[../op-quark-mass-core-g0-rescue-test-2026-06-25/Phase_FINAL_close.md]] (#41 — 𝒜 numeryczne)
- [[../../partial_proofs/quark_sector/dodatekX_quark_sector.tex]] prop:X-A-from-tube-tension, eq:X-K-msp-hypothesis
- [[../../core/_meta_latex/status_map.tex]] CG-1/CG-3 [SZKIC]
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] PR-025 (forward falsifier (b))
- [[../op-CG34-continuum-closure-2026-06-14/]] · [[../op-Csigma-coarse-graining-2026-06-14/]] · [[../op-uv-as-ngfp/]] (powiązane CG/UV)
