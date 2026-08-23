---
title: "FINDINGS — mini-test R1 (DIAGNOSTYKA): goly crossover -1/2->+2 = tautologia; teza zywa tylko jako (K,V)-niezgodnosc zakotwiczona w Phi=0"
date: 2026-06-27
type: phase0-findings
status: "DIAGNOSTIC RESULT (nie werdykt cyklu; brama R1 rozbrojona z konkretnym przeformulowaniem)"
parent: "[[Phase0_R1_sketch.md]]"
scoping: "[[../../meta/SCOPING_op-amplitude-density-phase-bridge_2026-06-27.md]]"
locked_inputs: ["#49 (alpha_eff=-1/2, e(p)=1/p-2), IMMUTABLE"]
artifacts: ["[[R1_test.py]] / [[R1_test.txt]] (value-blind, nienaruszony)", "[[R1_test_v2.py]] / [[R1_test_v2.txt]] (korekta 2 artefaktow)"]
anti_lakatos_note: "v1 zachowany jako zapis value-blind; 2 artefakty (T3 punkt nie-stacjonarny V'(1)=gamma/3; T4 brzeg odciety przez positive) zdiagnozowane JAWNIE, nie zamiecione; poprawione w v2. 0 hardcoded. Endpoint +4 nie wbudowany."
---

# FINDINGS — mini-test R1

## §0 — Wynik w jednym zdaniu

**Gładki crossover wykładnika `−½ → +2` „w objętości" jest tautologią zamiany zmiennej
(`Φ=σ²` daje to samo pole kanoniczne `χ=√2·√Φ`); teza dwufazowa przeżywa wyłącznie w ostrej
postaci: „czy pary `(K_sub, V_sub)` i `(K_TGP, V_M911)` są NIE-równoważne przez redefinicję
pola" — co jest pytaniem o **niezgodność potencjałów w ramie kanonicznej**, fizycznie
zakotwiczonym jedynie w `Φ=0` (nukleacja).**

## §1 — Co test rozstrzygnął (DECYDUJĄCE)

| Wynik | Test | Status |
|---|---|---|
| `K(Φ)` pojedynczego pola = czysta konwencja; `χ(Φ)=σ` dla dowolnego `p` | T1 (v1) | ✅ decydujące |
| `e(p)=1/p−2` zreprodukowane EXACT; `p=1→−1`, `p=1/6→+4`, `p=2→−3/2` | T1 (v1) | ✅ |
| Obie ramy (gęstość `Φ`, amplituda `σ`) dają TO SAMO pole kanoniczne `χ=√2√Φ` | T3' (v2) | ✅ decydujące |
| ⟹ **goły crossover = TAUTOLOGIA** (redefinicja zachowuje całą fizykę) | T1+T3' | ✅ |
| Naiwny dressing metryczny `√−g∝Φ^m` NIE mostkuje: wymaga `m=5`; TGP ma `m=1` (net 0), `m=2` (net 1) | T2 (v1) | ✅ |
| Mapy kanoniczne par różnią się funkcyjnie: `χ_TGP∝Φ³` vs `χ_sub∝Φ^{1/2}` | T3'' (v2) | ✅ |
| Jedyny lokus nieodwracalności redefinicji = `Φ=0` | T4' (v2) | ✅ |

## §2 — Dwa artefakty v1 (zdiagnozowane jawnie, anti-Lakatos)

1. **T3 v1 = False (artefakt):** „masa" `V''/2K` policzona w `ψ=1`, które **NIE jest ekstremum**
   `V_M911` (bo `V'(1)=γ/3≠0` — cecha TGP: próżnia podtrzymywana źródłem, `W(1)=γ/3`).
   Formuła zakłada `V'=0`. Nierówność `γ` vs `2γ/3` ⟹ artefakt punktu, NIE dowód fizyki.
   Poprawka v2 (T3') liczy pole kanoniczne wprost ⟹ ramy identyczne.
2. **T4 v1 = False (artefakt):** `sigma` zadeklarowane `positive` ⟹ sympy odciął brzeg `σ=0`.
   Lokus `F'=0` JEST na `Φ=0`. Poprawka v2 (T4') z `sig` real potwierdza brzeg.

## §3 — Przeformułowanie R1 (ostra, poprawna postać)

Redefinicja pola `Φ=F(σ)` zachowuje **całą** fizykę (masy, S-macierz). Zatem:

> **Crossover jest REALNY wtedy i tylko wtedy, gdy substrat `(K_sub=Φ⁻¹, V_sub)` i makro-TGP
> `(K_TGP=Φ⁴, V_M911)` NIE są równoważne przez redefinicję pola.** Test: wyraź oba potencjały
> w ich polach kanonicznych `χ_sub∝Φ^{1/2}`, `χ_TGP∝Φ³` i sprawdź, czy `V_sub(χ)` = `V_M911(χ)`
> (mod stała/skala). **Równe ⟹ ta sama teoria (crossover fikcyjny).** **Różne ⟹ substrat i TGP
> to RÓŻNE teorie** — a wtedy `−½` z #49 znaczy „substrat to inna teoria niż TGP", nie „zła rama".

To przesuwa cel z „wykładnika" na **niezgodność potencjałów** — i czyni **`V_sub` (potencjał
substratu) jedynym brakującym wejściem**, którego nie ma na Phase 0.

## §4 — Konsekwencje dla trzech osi usera

- **Oś A (wykres `e_eff(ℓ)`):** w objętości rysowałby **tautologię** (stałą wyznaczoną przez
  wybór zmiennej). Sensowna „ewolucja" istnieje tylko wokół `Φ=0`. ⟹ przeskalować oś A z
  „bulk crossover" na „zachowanie w nukleacji `Φ→0`".
- **Oś B (realność `K_eff(Φ,ℓ)`):** odpowiedź — gołe `K` NIE jest fizyczne; realny obiekt to
  **para `(K,V)`** + zachowanie na brzegu `Φ=0`. „Narodziny metryki" = ewentualna niezgodność
  par, nie zmiana wykładnika.
- **Oś C (stabilność/selekcja):** przeformułować na „czy dwie teorie `(K,V)` są różne i czy
  `α=2` jest selekcjonowane na brzegu nukleacji" (spójne z R5: α=1 vs α=2 vs −½).

## §5 — Bramka R1: ROZBROJONA warunkowo

R1 NIE jest już „ślepą" bramką: wiemy, że (a) naiwna teza (gładki bulk crossover) jest martwa
(tautologia), (b) żywa teza ma jedno, konkretne wejście: **`V_sub` + test zgodności potencjałów
w ramie kanonicznej**, (c) arena fizyczna = `Φ=0` (nukleacja), nie objętość.

**Następna bramka przed Phase 0 (zaktualizowana kolejność):**
1. **Pozyskać `V_sub`** (potencjał substratu w zmiennej amplitudowej `σ`) — z `H_Γ` / aksjomatów.
   Bez niego test z §3 jest niewykonalny. To nowa bramka egzystencjalna #0.
2. Następnie R3 (czy emergent-metric niezależna od α=2) — bez zmian.
3. R4/R5 (sektor mas; α=1 vs α=2) — bez zmian.

## §6 — Czego NIE wolno wywnioskować

- NIE: „udowodniono, że substrat ≠ TGP". Test pokazał tylko, że to jest **właściwe pytanie**
  i że odpowiedź wymaga `V_sub` (Phase 1). #49 LOCKED bez zmian.
- NIE: „α=2 zderywowane / obalone". Status `selekcja na gęstości` (status_map l.72) nietknięty.
- NIE: redukcja `N_free`. Zero zmian w ledgerze #42.

## Cross-references
- [[Phase0_R1_sketch.md]] · [[R1_test.py]] · [[R1_test.txt]] · [[R1_test_v2.py]] · [[R1_test_v2.txt]]
- [[../../meta/SCOPING_op-amplitude-density-phase-bridge_2026-06-27.md]] (R1 §4)
- [[../op-CG-alpha-eff-convergence-2026-06-26/Phase_FINAL_close.md]] (#49, LOCKED)
