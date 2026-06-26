---
title: "Phase 0 — pre-registration LOCK: rozstrzygnąć werdykt sektora radiacyjnego via DISFORMAL VIABILITY (sygnatura g_eff + DOF count slaved-TT + niezależność od O12), zastępując nierobustny argument induced-TT"
date: 2026-06-14
type: phase-balance
phase: 0
status: 🔒 LOCKED 2026-06-14 (przed formalnym rachunkiem; AUDIT scouting NIE liczy się jako rachunek cyklu)
cycle: op-disformal-hamiltonian-viability-2026-06-14
scoping: "[[../../meta/SCOPING_op-disformal-hamiltonian-viability_2026-06-14.md]]"
authorization: "User 2026-06-14 (sesja #26): 'działaj'"
parent_verdicts_LOCKED:
  - "PR-025 TRIGGERED — LOCKED"
  - "op-gravitational-sector-survival INDETERMINATE — LOCKED"
  - "op-disformal-radiation-resolution UNDERDETERMINED — LOCKED"
  - "op-disformal-stability Phase 1 (B<0 zdrowy skalar) POPRAWNE; Phase 2 (induced-TT) NIEROBUSTNY — audyt"
tags: [pre-registration, disformal-viability, signature, DOF-count, O12-independence, anti-Lakatos-LOCKED]
---

# Phase 0 — pre-registration LOCK

## §1 — Pytanie
Czy intersekcja **{$g_{\rm eff}$ Lorentzowska} ∩ {skalar bez ghost/instability} ∩ {nietrywialne
ekranowanie $\dot P_b$}** jest PUSTA dla każdego $B(\Phi)$ (znak i magnituda, niezależnie od O12)?
Cykl strukturalny: 0 danych. Reuse operatorów $Z^{\mu\nu}$, $g_{\rm eff}$ (EXACT, LOCKED).

## §2 — Wejścia LOCKED (read-only)
| Input | Forma | Źródło |
|---|---|---|
| Operator skalara | $Z^{\mu\nu}=2(A-bX)\eta^{\mu\nu}-4b\,\partial\bar\phi\partial\bar\phi$; $b=B/M_*^4$, $u=bX/A$ | survival Filar I / op-disformal-radiation-resolution |
| Metryka efektywna | $g_{\rm eff}=A\eta+b\,\partial\bar\Phi\partial\bar\Phi$ (disformal, `hyp:disformal`) | sek08 |
| Zdrowie skalara | no-ghost $u<1$; gradient $u<1/3$ ∨ $u>1$; zdrowy+ekranujący ⟺ B<0 | op-disformal-stability Phase 1 (POPRAWNE) |
| Skaling ekranowania | $S=1/\|1-u\|$ (heurystyka, leading-order) | op-disformal-radiation-resolution Phase 1 |
| Metryka emergentna | $g_{\rm eff}=g_{\rm eff}(\Phi,\partial\Phi)$ ⇒ δg slaved do δΦ | `rem:GW-scope-2026` |
| Konwencja | mostly-plus; statyczne $X=(\nabla\bar\Phi)^2>0$ | LOCKED Phase 0 (forbidden #4) |

## §3 — Falsyfikatory LOCKED (reguły IMMUTABLE; werdykt z flag)
| ID | Test | Reguła |
|---|---|---|
| **F-VIA-A** | sygnatura $g_{\rm eff}$: wartości własne vs $u$ | wartość radialna $A(1+u)<0$ dla $\|u\|>1$ (B<0) ⟹ `SIGNATURE-FLIP`; B>0 zawsze Lorentz |
| **F-VIA-B** | DOF count: czy induced-TT ma własny człon kinetyczny | δg∝δΦ (slaved) ⟹ `SLAVED` (argument induced-TT Phase 2 formalnie void); własny kinetyk ⟹ `INDEPENDENT` |
| **F-VIA-C** | trylemat $\forall B$: {Lorentz}∩{skalar zdrowy}∩{screening $\|u\|\gtrsim1$} | = ∅ ⟹ `EMPTY`; istnieje okno ⟹ `WINDOW` |
| **F-VIA-D** | (warunkowy) czy nietrywialne ekranowanie wymaga $\|u\|\gtrsim1$ | tak ⟹ trylemat domyka (qualitatively robust); nie ⟹ przelicz |
| **F-VIA-E** | agregat z flag | (A=FLIP ∧ B=SLAVED ∧ C=EMPTY) ⟹ **BROKEN-via-viability**; (C=WINDOW) ⟹ **NOT-BROKEN** (sign+mag-pin) |

## §4 — Pre-derywacja (oczekiwanie: BROKEN-via-viability; rachunek nadrzędny)
B<0,|u|>1: flip sygnatury (g_eff niefizyczna) · B>0,|u|>2 (screening): ghost skalara (L'<0) ·
B<0,|u|<1: g_eff OK, ale brak ekranowania (S→1) → PR-025 konforemny stoi. Intersekcja oczekiwana ∅.
**Jedyny element nie-twardy: link „screening ⇒ |u|≳1"** (qualitatively robust: perturbacyjnie mały
człon disformalny nie da O(1) supresji; specyficzny skaling 1/|1−u| dziedziczony jako heurystyka).

## §5 — Forbidden moves (IMMUTABLE)
1. Liczby PR-025/survival/parent LOCKED.
2. Zakaz induced-TT jako DOWODU (tylko jako „błędna ścieżka"; F-VIA-B formalnie ją zamyka).
3. Zakaz strojenia B/M_* (O12) pod wynik; trylemat O12-niezależny (F-VIA-C).
4. Konwencja sygnatury+znak X zafiksowana (§2).
5. Skaling screeningu = heurystyka dziedziczona; jawnie oznaczony (F-VIA-D), nie udawany jako twardy.
6. Zakaz miękkiego domknięcia ORAZ przedwczesnego negatywu — werdykt wyłącznie z flag.
7. Budżet nowych stałych 0.

## §6 — Anti-Lakatos checklist
✓ Phase 0 LOCKED przed formalnym rachunkiem. ✓ Zbiór flag CLOSED. ✓ Werdykt z flag. ✓ BROKEN i
NOT-BROKEN oba pre-deklarowane. ✓ Element heurystyczny (skaling screeningu) jawnie oznaczony jako
jedyna nie-twarda przesłanka. ✓ Werdykty poprzedników nietknięte. ✓ Budżet stałych 0.
