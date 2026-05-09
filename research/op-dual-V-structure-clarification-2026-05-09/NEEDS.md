---
title: "NEEDS — op-dual-V-structure-clarification-2026-05-09"
date: 2026-05-09
type: needs-list
status: COMPLETE
parent: "[[./README.md]]"
---

# NEEDS — op-dual-V-structure-clarification

## CRITICAL needs (D1-D3)

### D1: Sympy formal verification dual-V mathematical consistency
**Priority:** CRITICAL
**Phase:** 1 (sympy)
**Status:** OPEN

**Description:** Sprawdzić czy V_M9.1'' (gravity) + V_orig (matter)
**simultaneously valid** w TGP framework jest mathematically consistent.

**Resolution path:** sympy:
- Verify V_M9.1'' EOM (R3 ODE) i V_orig EOM (Phi expansion) operate w
  **różnych regimach** (static gravity vs dynamic matter)
- Check że matter Lagrangian L_mat z V_orig nie modify gravity sektor
- Document kompatybilność

### D2: G.0 explicit scope verification (text reading)
**Priority:** CRITICAL
**Phase:** 1 (sympy + reading)
**Status:** ✅ PARTIAL — pre-cycle finding (A4 marker)

**Description:** Eksplicite cytowanie G.0 statements ograniczających V_M9.1''
do gravity sector.

**Resolution path:** Phase 1 sympy text comparison cytat → check.

### D3: Sek08a annotation update recommendation
**Priority:** HIGH (deliverable)
**Phase:** 1 (results)
**Status:** OPEN

**Description:** Sformułuj **precise recommended annotation update** dla
sek08a linie 95-110 (DEPRECATED 2026-05-02).

**Pre-cycle proposed text:**
```
V_orig = -β Φ³/(3 Φ_0) + γ Φ⁴/(4 Φ_0²)
[DEPRECATED FOR GRAVITATIONAL SECTOR 2026-05-02 via G.0 closure;
 see prop:V-M911-canonical for gravity replacement.
 Matter sector usage MAINTAINED — see op-dual-V-structure-clarification-2026-05-09
 (A4 marker realization). Used in Phase 5 Mach inertia, T-Λ ρ_vac coupling.]
```

## SUPPORTING needs (D4-D5)

### D4: Cross-cycle re-classification post-Path-C
**Priority:** HIGH (deliverable)
**Phase:** 1 (results)
**Status:** OPEN

**Description:** Updated table z V usage per cykl po Path C confirmation.

| Cykl | V usage | Sektor | Status post-Path C |
|---|---|---|---|
| G.0 closure | V_M9.1'' | gravity | ✅ canonical |
| T-Λ closure | V_orig | matter | ✅ valid (no longer "residual gap") |
| Phase 5 MAG | V_orig | matter | ✅ valid (no longer "ambiguous") |
| sek08a master | mixed (annotated) | both | ⚠️ annotation needs update |

### D5: Update parent cycles z dual-V resolution
**Priority:** MEDIUM (cross-cycle)
**Phase:** 1 (post-results)
**Status:** OPEN

**Description:** Update:
- [[../op-V-canonical-consistency-audit-2026-05-09/]] — residual gaps RESOLVED
- [[../op-MAG-Phase5-V-reference-clarification-2026-05-09/]] — Path C confirmed
- [[../op-Phi-vacuum-scale-2026-05-09/]] — P11 BLOCKER fully resolved

## Priority matrix

| ID | Title | Priority | Phase | Status |
|---|---|---|---|---|
| D1 | Dual-V mathematical consistency | CRITICAL | 1 | OPEN |
| D2 | G.0 explicit scope verification | CRITICAL | 1 | ✅ PARTIAL |
| D3 | Sek08a annotation recommendation | HIGH | 1 | OPEN |
| D4 | Cross-cycle re-classification | HIGH | 1 | OPEN |
| D5 | Parent cycles update | MEDIUM | 1 | OPEN |

## Cross-references

- [[./README.md]]
- [[./Phase0_balance.md]]
- [[../op-g0-r3-from-canonical-projection/Phase1_results.md#5.3]] — A4 marker
- [[../op-MAG-Phase5-V-reference-clarification-2026-05-09/Phase1_clarification_results.md]] — Path C origin
