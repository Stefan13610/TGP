---
title: "Phase FINAL — op-A3-alpha-resolution CLOSE. Werdykt value-blind: α=2 = DERIVED-INCONSISTENCY. Relacja EL (α=p/2 dla K=f^p) POPRAWNA. ALE: substrat K(u)∝u⁴ (prop:substrate-action) + ontologia Φ=⟨ŝ²⟩∝u² po POPRAWNEJ zmianie zmiennych (sek10 eq:kinetic_macro) daje K_eff(ψ)∝ψ (liniowy) → α=1/2, NIE 2. thm:D-uniqueness C3 dostaje α=2 przez POMINIĘCIE transformacji Φ∝u² (traktuje macro Φ/Φ₀ jakby było mikro u). Wewnętrzna sprzeczność też w sek10 (eq:kinetic_macro liniowy vs claim α=2). α=2 spójne tylko gdy Φ∝u (amplituda) — sprzeczne z Φ=⟨ŝ²⟩. Rekomendacja korekty rdzenia (bez edycji core)."
type: phase_close
status: 🟢 CLOSED — DERIVED-INCONSISTENCY (α=2 nie wynika z substratu pod Φ=⟨ŝ²⟩) 2026-06-14
phase: FINAL
cycle: op-A3-alpha-resolution-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14 (sesja #31): 'tak działaj z op-A3-alpha-resolution'"
script: "[[./Phase1_alpha_sympy.py]]"
output: "[[./Phase1_output.txt]]"
verdict: "DERIVED-INCONSISTENCY: relacja EL poprawna; α=2 wymaga Φ∝u (amplituda), sprzeczne z Φ=⟨ŝ²⟩∝u²; pod Φ=⟨ŝ²⟩ poprawne α=1/2 (K~u⁴) lub 0 (K~u²)"
anti_lakatos_lock: PRESERVED
core_recommendation: "thm:D-uniqueness C3 + sek10 §K_to_f + dodatekQ2 A3: ujednolicić; rozstrzygnąć Φ∝u vs Φ∝u²"
---

# Phase FINAL — CLOSE (op-A3-alpha-resolution)

> **Werdykt (value-blind, sympy 5/5):** $\alpha=2$ = **DERIVED-INCONSISTENCY**.
> Relacja EL ($\alpha=p/2$ dla $K=f^p$) jest **poprawna** (sek08 niepodważone). **Ale** geometryczne
> sprzężenie substratu $K(u)\propto u^4$ + ontologia $\Phi=\langle\hat s^2\rangle\propto u^2$ daje — po
> **poprawnej** zmianie zmiennych — $\alpha=\tfrac12$, **nie 2**. $\alpha=2$ pojawia się w
> thm:D-uniqueness C3 tylko przez **pominięcie** transformacji $\Phi\propto u^2$.

## §1 — Co dokładnie ustalono (sympy, value-blind)
| Test | Wynik |
|---|---|
| **F-α-A** | **PASS** — relacja EL $\alpha=K'f/(2K)=p/2$ dla $K=f^p$ poprawna (reprodukuje sek08 eq:K-ode). C3 arytmetyka OK: GDYBY $K(\Phi/\Phi_0)=(\Phi/\Phi_0)^4$, to $\alpha=2$. |
| **F-α-B** | **PASS** — substrat $K(u)\propto u^{2s}$ + $\Phi\propto u^2$ → $K_{\rm eff}(\psi)\propto\psi^{s-1}$ → $\alpha=(s-1)/2$. $s{=}1(K{\sim}u^2){:}~\alpha{=}0$; $s{=}2(K{\sim}u^4){:}~\alpha{=}\tfrac12$. |
| **F-α-C** | **PASS** — $\alpha=2$ **nie wynika** z substratu pod $\Phi\propto u^2$ → **INCONSISTENCY**. |
| **F-α-D** | **PASS** — $\alpha=2$ ⟺ $\Phi\propto u$ (amplituda) [sprzeczne z $\Phi=\langle\hat s^2\rangle$] lub $K\propto u^{10}$ [brak derywacji]. |

## §2 — Anatomia niespójności (precyzyjnie)
Łańcuch claimowany (thm:D-uniqueness, sek08 eq:967–1050):
1. $\mathcal L_{\rm kin}=-\tfrac12K(\varphi)(\nabla\varphi)^2$, $\varphi=\Phi/\Phi_0$ (macro). [eq:971]
2. C1+C2: $K(\varphi)=C\varphi^{2\alpha}$. [eq:K-ode]
3. C3: substrat → $K(\varphi)=K_{\rm geo}\varphi^4$ ⟹ $2\alpha=4$ ⟹ $\alpha=2$.

**Błąd w kroku 3:** prop:substrate-action daje $K\propto u^4$ w polu **mikroskopowym** $u$ (spiny, $H_\Gamma=-J\sum(u_iu_j)^2$). Krok 3 wstawia tę potęgę do pola **makroskopowego** $\varphi=\Phi/\Phi_0$, **jak gdyby** $\varphi=u$. Ale ontologia (sek08 poziom 0; sek10 eq:Phi_phi_id) mówi $\Phi=\langle\hat s^2\rangle\propto u^2$, więc $\varphi=\Phi/\Phi_0\propto u^2$. Poprawnie: $K\propto u^4=(\varphi)^2$, czyli $2\alpha=2$, $\alpha=1$ (proste liczenie potęg) — a z **pełną** transformacją (czynnik $1/\psi$ z $(\nabla u)^2$, sek10 eq:kinetic_macro) $K_{\rm eff}(\psi)\propto\psi^1$ → $\alpha=\tfrac12$.

**Sprzeczność jest też WEWNĄTRZ sek10:** eq:kinetic_macro (linie 114–119) jawnie transformuje $K_{\rm geo}u^4(\nabla u)^2=K_0\,\psi(\nabla\psi)^2$ — **liniowy** w $\psi$ (→ $\alpha=\tfrac12$) — a jednocześnie linia 137 twierdzi „$\alpha=2$". Liniowe $K_{\rm eff}$ i $\alpha=2$ są niezgodne.

## §3 — Obie strony (anti-Lakatos, uczciwie)
**PRO $\alpha=2$:**
- Relacja EL poprawna; arytmetyka C3 wewnętrznie spójna *gdy* $K(\Phi/\Phi_0)=(\Phi/\Phi_0)^4$.
- MC LK-1a-g raportuje „$\alpha=2$ w 1.2σ" — ale dodatekQ2 podaje $\alpha_{\rm eff}=6.48\pm3.82$, więc 1.2σ jest **bardzo słabe** (przedział obejmuje też $\alpha=\tfrac12$).
- $\alpha=2$ może być **fenomenologicznie wymagane** (PPN, profil solitonu, widmo mas) — wtedy jest **niezależnym wejściem**, nie wyprowadzonym z substratu.

**CONTRA $\alpha=2$ (z substratu):**
- Zmiana zmiennych $\Phi\propto u^2$ jednoznacznie (sympy + sek10 eq:kinetic_macro) daje $\alpha=\tfrac12$.
- $\alpha=2$ wymagałby $\Phi\propto u$ (amplituda), co łamie $\Phi=\langle\hat s^2\rangle$ — fundamentalną definicję.

## §4 — Werdykt
**$\alpha=2$ = DERIVED-INCONSISTENCY.** Twierdzenie „α=2 **wyprowadzone** z geometrycznego sprzężenia
substratu" (thm:D-uniqueness) jest **wewnętrznie niespójne** z ontologią $\Phi=\langle\hat s^2\rangle\propto u^2$.
- Pod $\Phi=\langle\hat s^2\rangle$: poprawne $\alpha=\tfrac12$ (dla $K\sim u^4$) lub $\alpha=0$ (dla $K\sim u^2$).
- $\alpha=2$ jest spójne **tylko** gdy kanoniczne pole kinetyczne to **amplituda** $u$ (Φ∝u), nie gęstość $\langle\hat s^2\rangle$.

**To rozstrzyga residual op-CG34 Phase 3** (i niespójność lematu A3): nie była to pomyłka w moim rachunku — jest realna niespójność konwencji micro/macro w rdzeniu.

## §5 — Opcje naprawy (konstruktywne; do decyzji rdzenia)
1. **Opcja B (najprawdopodobniej zamierzona):** zadeklarować, że **kanonicznym polem kinetycznym jest amplituda** $\varphi$ (z $K(\varphi)=K_{\rm geo}\varphi^4$, $\alpha=2$), a $\Phi=\langle\hat s^2\rangle$ to **osobna** wielkość (gęstość/faza metryczna) używana w ontologii i metryce. Wtedy trzeba **rozróżnić** dwa pola w notacji (obecnie zlane jako „Φ"/„φ"). To usuwa sprzeczność, ale wymaga przeglądu, gdzie używa się $\Phi\propto u^2$ vs $\varphi$.
2. **Opcja C:** przyjąć $\alpha=2$ jako **fenomenologiczny aksjomat** (wymagany przez PPN/masy), a thm:D-uniqueness przeklasyfikować z „wyprowadzenia" na „spójnościowy" (selekcja w klasie, jak sugeruje rem:alpha2-pivot-status).
3. **Odrzucić Opcję A** ($K\sim u^{10}$ — sztuczne, brak derywacji).

## §6 — Rekomendacje dla rdzenia (NIE wykonane — forbidden move #2; zgłoszone)
1. **sek08 thm:D-uniqueness C3:** doprecyzować, w jakim polu ($u$ amplituda czy $\Phi=\langle\hat s^2\rangle$) zachodzi $K\propto(\cdot)^4$. Obecnie eq:971 ($\varphi=\Phi/\Phi_0$) + C3 (substrat $u^4$) są niespójne pod $\Phi\propto u^2$.
2. **sek10 §K_to_f:** usunąć sprzeczność eq:kinetic_macro (liniowy $K_{\rm eff}$) ↔ claim α=2; poprawić rozwinięcie linia 167 ($e^{2\ln g}\to1+2\ln g$, nie $1+4\ln g$).
3. **dodatekQ2 lemat A3:** ujednolicić (premisa $Z\sim\phi^2$ ↔ konkluzja $K_1\sim1/\Phi$) zgodnie z wybraną opcją (§5).
4. **Rozróżnić notacyjnie** pole amplitudowe $\varphi$ od gęstości $\Phi=\langle\hat s^2\rangle$ w sek08/sek10.

## §7 — Anti-Lakatos (checklist)
✓ Werdykt **wyliczony** sympy (value-blind), nie wybrany; relacja EL **niepodważona** (uczciwie: błąd jest tylko w C3, nie w całym formalizmie).
✓ **Obie strony** przedstawione (PRO: MC 1.2σ słabe, fenomenologia; CONTRA: zmiana zmiennych).
✓ $\alpha=2$ **NIE odrzucone jako fenomenologia** — odrzucone tylko jako „wyprowadzone z substratu" pod $\Phi=\langle\hat s^2\rangle$.
✓ Rdzeń **NIE edytowany**; opcje naprawy (§5) + rekomendacje (§6) zgłoszone osobno.
✓ Niespójność potwierdzona z DWÓCH niezależnych miejsc rdzenia (thm:D-uniqueness vs sek10 eq:kinetic_macro) — nie pojedynczy odczyt.

## §8 — Status końcowy
**🟢 CLOSED — DERIVED-INCONSISTENCY.** Cykl rozstrzygnął residual op-CG34/lemat A3: niespójność α=2↔K(φ)
jest **realna** (konwencja micro/macro), potwierdzona sympy i z dwóch miejsc rdzenia. $\alpha=2$ pozostaje
**fenomenologicznie wiarygodne**, ale **nie jest wyprowadzone** z substratu pod $\Phi=\langle\hat s^2\rangle$ —
wymaga albo rozróżnienia pól (amplituda vs gęstość, Opcja B), albo statusu aksjomatu (Opcja C). Rekomendacje
zgłoszone; decyzja należy do rdzenia (osobna autoryzacja).
