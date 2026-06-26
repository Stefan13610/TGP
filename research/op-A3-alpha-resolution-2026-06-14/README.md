---
title: "op-A3-alpha-resolution — rozstrzygnięcie niespójności α=2 ↔ K(φ) (lemat A3 / thm:D-uniqueness)"
type: research_cycle
status: "🟢 CLOSED — DERIVED-INCONSISTENCY (α=2 nie wynika z substratu pod Φ=⟨ŝ²⟩) 2026-06-14"
phase: FINAL
folder_status: active
created_date: 2026-06-14
closed_date: 2026-06-14
authorization: "User 2026-06-14 (sesja #31): 'tak działaj z op-A3-alpha-resolution'"
origin: "Residual op-CG34-continuum-closure Phase 3 (niespójność lematu A3)"
verdict: "DERIVED-INCONSISTENCY: relacja EL poprawna; α=2 ⟺ Φ∝u (amplituda), sprzeczne z Φ=⟨ŝ²⟩∝u²; pod ontologią α=1/2 (K~u⁴) lub 0 (K~u²)"
anti_lakatos_lock: PRESERVED
---

# op-A3-alpha-resolution (🟢 CLOSED — DERIVED-INCONSISTENCY)

> **Pytanie:** Czy α=2 (człon (∇Φ)²/Φ w równaniu pola TGP) **wynika** z geometrycznego sprzężenia
> substratu K(φ)∝φ⁴ przez Φ=φ²? **Odpowiedź (sympy, value-blind, 5/5): NIE — to niespójność konwencji.**

## Werdykt
- **Relacja EL poprawna:** dla L=−½K(φ)(∇φ)², K=φ^p → współczynnik (∇φ)²/φ = p/2 = α. sek08 **niepodważone**.
- **Niespójność w thm:D-uniqueness C3:** wstawia substratowe K∝u⁴ (pole **mikro** u, spiny) do pola **makro** φ=Φ/Φ₀, pomijając Φ=⟨ŝ²⟩∝u². Po poprawnej zmianie zmiennych (sek10 eq:kinetic_macro) K_eff(ψ)∝ψ (liniowy) → **α=1/2, nie 2**.
- **Potwierdzone z 2 miejsc rdzenia:** thm:D-uniqueness (α=2) vs sek10 eq:kinetic_macro (liniowy K_eff → α=1/2) są wzajemnie sprzeczne.
- **α=2 ⟺ Φ∝u (amplituda)** — sprzeczne z ontologią Φ=⟨ŝ²⟩∝u².

## Znaczenie
- **Rozstrzyga** residual op-CG34 Phase 3 i niespójność lematu A3: nie był to błąd rachunku — realna niespójność micro/macro.
- **α=2 pozostaje fenomenologicznie wiarygodne** (PPN, masy, soliton), ale **nie jest wyprowadzone** z substratu pod Φ=⟨ŝ²⟩.

## Opcje naprawy (do decyzji rdzenia)
- **B (prawdopodobnie zamierzona):** kanoniczne pole kinetyczne = **amplituda** φ (K∝φ⁴, α=2); Φ=⟨ŝ²⟩ = osobna gęstość. Rozróżnić notacyjnie.
- **C:** α=2 jako aksjomat fenomenologiczny; thm:D-uniqueness = argument spójnościowy (selekcja), nie wyprowadzenie.

## Pliki
- [[./Phase0_lock.md]] — LOCK pytania, materiał dowodowy z rdzenia, reguła decyzyjna.
- [[./Phase1_alpha_sympy.py]] + [[./Phase1_output.txt]] — derywacja value-blind (5/5 PASS).
- [[./Phase_FINAL_close.md]] — werdykt + anatomia niespójności + rekomendacje rdzenia (§5–6).

## Rekomendacje rdzenia (zgłoszone, NIE wykonane)
1. sek08 thm:D-uniqueness C3 — doprecyzować pole (amplituda u vs gęstość Φ=⟨ŝ²⟩).
2. sek10 §K_to_f — usunąć sprzeczność eq:kinetic_macro↔claim α=2; poprawić e^{2ln g}→1+2ln g (nie 1+4ln g, linia 167).
3. dodatekQ2 lemat A3 — ujednolicić premisę/konkluzję.
4. Rozróżnić notacyjnie φ (amplituda) od Φ=⟨ŝ²⟩ (gęstość).
