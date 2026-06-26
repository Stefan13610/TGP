---
title: "op-amplitude-density-global-audit"
date: 2026-06-16
type: research-cycle
status: CLOSED
verdict: INCONSISTENT-LABELING
---

# op-amplitude-density-global-audit (2026-06-16)

**Cel:** audyt odpowiedzialnościowy po integracji Opcji B (#31). Zweryfikować value-blind,
czy rozróżnienie amplituda φ / gęstość Φ=⟨ŝ²⟩ jest globalnie spójne w sekcjach downstream
(sek08c metryka, sek07/sek06 predykcje, soliton, masy, sek02/sek09 równanie pola, dodatekQ/Q2).

**Werdykt: INCONSISTENT** (4× G-INCONSISTENT, reguła LOCKED).
- Nośna **fizyka rdzenia SPÓJNA**: pole kanoniczne = **gęstość Φ** (`def:Phi`), **α=2 postulowane
  na gęstości** (selekcja C1–C3, `rem:alpha2-pivot-status`); metryka/soliton/masy/PPN wszystkie
  w tej samej zmiennej (11× G-CONSISTENT, 0 sprzeczności wewnątrz).
- Dewiant = **moja edycja Opcji B #31** (3 uwagi: sek08 `rem:amplitude-vs-density-alpha`,
  sek10 `rem:K_to_f_amplitude`, dodatekQ2 `rem:A3-correction-alpha`) + 1 arytmetyczna zaległość
  (sek10:145 `g²`→`g⁴`). Wprowadziły błędne ramowanie „pole kanoniczne = amplituda, α=½ w gęstości",
  sprzeczne z nośnym „Φ=gęstość, α=2 postulowane".

**Cykl obalił własne założenie:** Opcja B #31 NIE jest spójna w dół — wymaga przeramowania.

## Fazy
- [[Phase0_lock.md]] — protokół value-blind, gate'y, reguła werdyktu (LOCKED).
- [[Phase2_classification.md]] — inwentaryzacja (3 Explore) + klasyfikacja firsthand; tabele A/B.
- [[Phase_FINAL_close.md]] — werdykt + lista 4 poprawek + rekomendacja anty-overload φ.

## ✅ Poprawki ZASTOSOWANE (user „Działaj — wszystkie 4", 2026-06-16)
1. ✅ sek08 `rem:amplitude-vs-density-alpha` — przeramowane na Φ=gęstość pole kanoniczne, α=2 postulat.
2. ✅ sek10 `rem:K_to_f_amplitude` — usunięte „g=amplituda" (g≡ψ=Φ/Φ₀ gęstość); amplituda = reprez. substratu (α=½).
3. ✅ dodatekQ2 `rem:A3-correction-alpha` — przeramowane (α=2 postulat na gęstości).
4. ✅ sek10:145 `K_sub(g)=g²` → `K(g)=g⁴`.
- **Build: main.pdf 552 s., pdflatex 0 błędów fatalnych; nowe cross-refy rozwiązane.**
+ rekom. (niski prio., otwarte): osobny symbol mikro-amplitudy w sek10 §K_to_f (anty-overload φ).
