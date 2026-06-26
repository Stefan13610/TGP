---
title: "op-disformal-hamiltonian-viability — formalne rozstrzygnięcie werdyktu sektora radiacyjnego via DISFORMAL VIABILITY (sygnatura/hiperboliczność g_eff + DOF count slaved-TT + niezależność od O12), zastępując nierobustny argument induced-TT z op-disformal-stability Phase 2"
type: research_cycle
status: 🟡 CLOSED-RESOLVED — KOREKTA: kanał skalarno-disformalny BROKEN (pod-wynik), sektor radiacyjny UNDERDETERMINED (2026-06-14, adwersaryjna kontrola)
verdict_correction: "[[./ADVERSARIAL_REVIEW_2026-06-14.md]] — F-VIA-E 'BROKEN-via-viability' NADMIERNIE ZAOSTRZONE (błąd zakresu: σ_ab=właściwy radiator, B-niezależny, κ_E niepinowane). Geometria CONFIRMED (screening skalara BROKEN), zakres REFUTED. Poprawny status sektora: UNDERDETERMINED (=cykl-rodzic)."
phase: FINAL
folder_status: closed
final: "[[./Phase_FINAL_close.md]] — F-VIA-E LOCKED; sektor radiacyjny TGP_v1 SFALSYFIKOWANY; propagacja FOUNDATIONS CL-6 + REALITY_CONTACT_AUDIT + STATE; domyka op-disformal-stability (poprawny argument)"
phase0_lock: "[[./Phase0_balance.md]] 🔒 LOCKED 2026-06-14"
phase1_result: "[[./Phase1_derivation.md]] (sympy 5/5) — F-VIA-A=SIGNATURE-FLIP, F-VIA-B=SLAVED, F-VIA-C=EMPTY, F-VIA-D=NEEDS-|u|>~1, F-VIA-E=BROKEN-via-viability"
implied_verdict: "sektor radiacyjny TGP_v1 → BROKEN-via-viability (g_eff flip sygnatury |u|=1=r_V; trylemat ∅ O12-niezależnie; induced-TT formalnie void). Terminalny — Phase FINAL = świadoma autoryzacja."
created_date: 2026-06-14
registered_by: "user 2026-06-14 (sesja #26): audyt op-disformal-stability → 'C' (napraw + audyt werdyktu + scoping)"
scoping: "[[../../meta/SCOPING_op-disformal-hamiltonian-viability_2026-06-14.md]]"
parent_cycle: "[[../op-disformal-stability-2026-06-14/]] (Phase FINAL WSTRZYMANA — argument induced-TT nierobustny)"
audit_basis: "[[../op-disformal-stability-2026-06-14/AUDIT_verdict_2026-06-14.md]] + [[../op-disformal-stability-2026-06-14/AUDIT_verdict_sympy.txt]]"
cycle_category: "VERDICT-RESOLUTION (solidny argument viability zastępuje induced-TT; werdykt sektora radiacyjnego)"
expected_duration: "1–2 fazy; BROKEN-via-viability / NOT-BROKEN oba pełnoprawne"
predecessor_verdicts_LOCKED:
  - "PR-025 TRIGGERED — LOCKED"
  - "op-gravitational-sector-survival INDETERMINATE — LOCKED"
  - "op-disformal-radiation-resolution UNDERDETERMINED — LOCKED"
  - "op-disformal-stability Phase 1 (B<0 zdrowy skalar) — POPRAWNE; Phase 2 argument induced-TT — NIEROBUSTNY (audyt)"
independent_of: "O12 (trylemat ma być O12-niezależny), op-nucleation-dimensionality"
anti_lakatos_lock: "INHERITED; aktywny od Phase 0"
---

# op-disformal-hamiltonian-viability (REGISTERED-QUEUED — Phase 0 pending)

> **Status:** zarejestrowany po audycie op-disformal-stability. Wymaga własnego Phase 0 + „działaj".

## Pytanie wiodące

Audyt ([[../op-disformal-stability-2026-06-14/AUDIT_verdict_2026-06-14.md]]) wykazał EXACT:
$g_{\rm eff}=\mathrm{diag}(-A,A+bG^2,A,A)$, wartość własna radialna $A(1+u)$ **flipuje sygnaturę
przy $|u|=1$ ($=r_V$) dla $B<0$**. Trylemat:

| reżim | $g_{\rm eff}$ | skalar | screening |
|---|---|---|---|
| B<0, $\|u\|>1$ | flip ✗ | zdrowy | silny |
| B>0, $\|u\|>2$ | OK | ghost ✗ | silny |
| B<0, $\|u\|<1$ | OK | zdrowy | brak ✗ (→PR-025) |

> **Czy {g_eff Lorentz} ∩ {skalar zdrowy} ∩ {silne ekranowanie} = ∅ dla każdego B (niezależnie
> od O12)?** TAK ⟹ **BROKEN-via-viability** (solidny dowód, zastępuje induced-TT). Okno ⟹ **NOT-BROKEN**.

## Czym różni się od op-disformal-stability

op-disformal-stability dało **poprawną Phase 1** (B<0 = zdrowy skalar), ale **nierobustną Phase 2**
(BROKEN via $c_T$ induced-TT — tryb, który rdzeń `rem:GW-scope-2026` oznacza jako niefizyczny).
Ten cykl rozstrzyga **tę samą konkluzję poprawnym narzędziem**: sygnatura/viability $g_{\rm eff}$
(własność metryki, nie slaved-TT). Próg $|u|=1$, który Phase 2 przypisała „niestabilności tensora",
to **naprawdę** degeneracja $g_{\rm eff}$.

## Twarde wymogi Phase 0 (pełne w [[../../meta/SCOPING_op-disformal-hamiltonian-viability_2026-06-14.md]])

- F-VIA-A (sygnatura g_eff) / F-VIA-B (DOF count slaved-TT) / F-VIA-C (niezależność od O12) /
  F-VIA-D (skaling screeningu, warunkowy) / F-VIA-E (agregat).
- Zakaz induced-TT jako dowodu (tylko jako „błędna ścieżka"); zakaz strojenia B/M_*; konwencja zafiksowana.
- Budżet nowych stałych 0.

## Aktywacja

Nowy agent: czytaj [[../../meta/SCOPING_op-disformal-hamiltonian-viability_2026-06-14.md]] +
[[../op-disformal-stability-2026-06-14/AUDIT_verdict_2026-06-14.md]] (+ AUDIT_verdict_sympy) +
rdzeń `hyp:disformal`/`prop:cT`/`rem:GW-scope-2026`, potem user „działaj" → Phase 0 LOCK.
