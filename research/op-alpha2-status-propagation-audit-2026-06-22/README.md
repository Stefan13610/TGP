---
title: "op-alpha2-status-propagation-audit — audyt spójności propagacji statusu α=2 = selekcja na gęstości (po #31 op-A3 / #32)"
date: 2026-06-22
type: research-cycle
status: 🟢 CLOSED — DO-POPRAWY naprawione (P1+P2+P3+P4 zastosowane 2026-06-22; main.tex build exit 0, 553 str.)
session: "#36"
parent_cycles:
  - "[[../op-A3-alpha-resolution-2026-06-14/Phase_FINAL_close.md]]"
  - "[[../op-amplitude-density-global-audit-2026-06-16/Phase_FINAL_close.md]]"
template: "[[../op-sigma-status-propagation-audit-2026-06-20/README.md]]"
coordination: "[[../../STATE.md]]"
tags: [audit, status-propagation, alpha2, kinetic-coupling, selection-vs-derivation, parameter-counting, anti-Lakatos]
---

# op-alpha2-status-propagation-audit

## Cel
Audyt spójnościowy (value-blind, anti-Lakatos): znaleźć **wszystkie** miejsca w rdzeniu (.tex)
i dokumentach głównych/**submission**, gdzie kanoniczne sprzężenie kinetyczne **α = 2** (K(φ) = φ⁴)
jest opisane jako „wyprowadzone / derived / proved / algebraic theorem **z substratu**" w sposób
**niespójny** z zasądzonym w sesjach #31 (op-A3) i #32 (audit gęstości) statusem:

> **α = 2 NIE jest wyprowadzone z substratu pod ontologią Φ = ⟨ŝ²⟩ (gęstość).** Substrat dałby
> α = ½ (op-A3, sympy 5/5). α = 2 ma status **selekcji aksjomatycznej na gęstości** (jednoznaczność
> w klasie C1–C3), NIE derywacji. Relacja EL (α = p/2) i fenomenologia α = 2 (PPN, masy, Koide)
> pozostają niepodważone.

## Dlaczego ten cykl (priorytet)
Wskazany jako **najważniejsza rzecz do zbadania w TGP_v1**: α = 2 / K = φ⁴ jest logicznym
**korzeniem** całej teorii (metryka, solitony, masy leptonów, PPN). Rdzeń (sek08/sek10) został
naprawiony w #31/#32, ale naprawa **nie spropagowała** do warstwy publikacyjnej — co tworzy
żywą niespójność „manuskrypt vs streszczenie/list", widoczną dla recenzenta zewnętrznego.
Wzorzec identyczny z #35 (κ_E), stawka wyższa (dokumenty do submission).

## Pliki cyklu
- `Phase0_lock.md` — fakt referencyjny F + rozróżnienia R1–R4 + reguła klasyfikacji + zakres + falsyfikator + forbidden moves.
- `Phase0_balance.md` — gate (dotyka statusów epistemicznych; lean jawny).
- `Phase1_audit.md` — tabela 11 trafień + klasyfikacja SPÓJNE / DO-POPRAWY / NIEJASNE.
- `Phase_FINAL_close.md` — werdykt + rekomendacje P1–P4 + build verification (po autoryzacji).

## Werdykt (Phase FINAL, value-blind)
**DO-POPRAWY** — propagacja statusu z #31/#32 **NIEPEŁNA**. **8 trafień twardych:**
`README.md` (×2), `tgp_letter.tex`, `tgp_companion.tex`, `tgp_core.tex` (×2),
`status_map.tex`, `dodatekH`. **SPÓJNE (kotwice):** sek08 `rem:alpha2-pivot-status-pl`,
sek10:208–216, `TGP_FOUNDATIONS.md`. Wynik NIE jest negatywny.

## Deliverable wtórny (zgłoszony, NIE rozstrzygany)
Jeśli α = 2 to selekcja aksjomatyczna — nagłówek „**40 predykcji z 3 inputów**" zależy od
decyzji, czy warunki selekcji C1–C3 liczą się jako wejście. Audyt **ujawnia** tę zależność;
decyzja należy do rdzenia (osobna autoryzacja).

## Anti-Lakatos
Relacja EL (R1), fenomenologia α = 2 (R2), thm:alpha2-jako-selekcja-w-klasie (R3),
Φ = ⟨ŝ²⟩ kanoniczne (R4) — zalockowane jako nietykalne. Budżet nowych stałych = 0.
Rdzeń/papers NIE edytowane — rekomendacje zgłoszone, edycje czekają na autoryzację usera.

## Status
🟢 **CLOSED** — user autoryzował **Pełny P1+P2+P3 + P4** (2026-06-22). 11 edycji zastosowanych
(README ×3, letter ×2, companion ×2, core ×2, status_map ×2, dodatekH ×1). `main.tex` build
**exit 0, 553 strony**, brak nowych dangling refs. **Osobny finding:** 3 standalone papers mają
pre-existing fatal errors makr (`\ZZ`, `\gone^*`, `\Z`) niezwiązane z α=2 — do osobnej decyzji.
