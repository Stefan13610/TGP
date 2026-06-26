---
title: "op-PR004-SPARC-fit-execution — wykonanie LOCKED falsyfikatora PR-004: χ²_red(TGP Newton-baryon) vs χ²_red(MOND simple) na SPARC 175"
type: research_cycle
status: CLOSED-RESOLVED
phase: FINAL
folder_status: closed-resolved
claim_status: "CLOSED-RESOLVED PR004-TRIGGERED-FALSIFIED-MECHANISM (LOCKED 2026-06-13) — TGP g_eff[Φ̄ ≈ Φ₀] (Newton + wyłącznie bariony, S05) przegrywa z MOND simple na SPARC 175 przy t = 5.4σ (próg 5σ, reguła IMMUTABLE PR-004); Q1+Q2: 5.5σ. χ²_red GLOBAL 578 vs 50; mediany 85 vs 10.5; TGP lepsze tylko w 25/175 (HSB barionowo zdominowane). Per kontrakt: 'framework needs structural amendment, NOT continued recovery'."
created_date: 2026-06-12
closed_date: 2026-06-13
authorization: "User 2026-06-12: 'spróbujmy' (po audycie REALITY_CONTACT: PR-004 = LOCKED-PENDING-FIT, dane dostępne)"
PR_executed: "PR-004 (LOCKED 2026-05-13) — decision rule IMMUTABLE, wykonana mechanicznie"
data: "SPARC Lelli+2016 (astroweb.cwru.edu/SPARC): MassModels_Lelli2016c.mrt (3391 pkt, 175 galaktyk) + SPARC_Lelli2016c.mrt (Q-flags) — LITERATURE_ANCHORED, pobrane 2026-06-12, kopie w ./data/"
anti_lakatos_lock: PRESERVED
---

# op-PR004-SPARC-fit-execution

## Wynik (mechaniczny, reguła IMMUTABLE)

**PR-004 TRIGGERED:** χ²_red(TGP) > χ²_red(MOND simple) przy **t = 5.4σ** (paired per-galaxy,
bootstrap frac = 1.0000; Q1+Q2: 5.5σ) ⇒ mechanizm rotacyjny TGP z g_eff[Φ̄] background
**INSUFFICIENT**; per kontrakt: (a) osobna ρ_DM = S05 violated / (b) dedykowany mechanizm;
if_recovery_exhausted: **„framework needs structural amendment, NOT continued recovery"**.

| Miara | TGP (Newton+bariony) | MOND simple |
|---|---|---|
| χ²_red GLOBAL (3391 pkt) | **578.1** | 50.0 |
| χ²_red MEDIAN per-galaxy | **85.2** | 10.5 |
| Galaktyki wygrane | 25/175 | 150/175 |

Zero dopasowywanych parametrów w obu modelach (Υ_d = 0.5, Υ_b = 0.7, a₀ = 1.2×10⁻¹⁰ FIXED
pre-LOCK; guard FP6: zero wywołań optymalizatora).

## Pliki

- [[Phase0_balance.md]] — pipeline LOCKED PRZED danymi (Υ fixed, χ² def, operacjonalizacja 5σ,
  filtry, forbidden moves)
- [[Phase1_fit.py]] / [[Phase1_fit.txt]] — wykonanie 6/6 PASS (FP2: audyt poprawności —
  tożsamość inwersji 8×10⁻¹⁶, deep-MOND asymptota ✓; dyspozycja konwencji benchmarku ~2.0)
- [[Phase_FINAL_close.md]] — zamknięcie + dyspozycja strukturalna
- ./data/ — kopie SPARC (LITERATURE_ANCHORED)
