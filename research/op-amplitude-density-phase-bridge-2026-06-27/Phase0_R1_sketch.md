---
title: "Phase 0 (szkic) — mini-test R1: czy crossover Φ⁻¹→Φ⁺⁴ jest tautologią zmiennej, czy realną fizyką?"
date: 2026-06-27
type: phase0-sketch
status: "DIAGNOSTIC PRE-PHASE-0 (nie cykl; nie pre-rejestracja; brak werdyktu falsyfikującego; test rozpoznawczy bramki R1)"
parent: "[[../../meta/SCOPING_op-amplitude-density-phase-bridge_2026-06-27.md]]"
locked_inputs:
  - "#49 op-CG-alpha-eff-convergence: α_eff=−1/2 z substratu (e(p)=1/p−2, p=1 ⟹ e=−1), CLOSED-NEGATIVE, IMMUTABLE"
  - "status_map l.72 (α=2 = selekcja na gęstości), l.489-501 (α=1 vs α=2)"
anti_lakatos_note: "Test DIAGNOSTYCZNY (nie falsyfikator). Liczy 4 pytania symboliczne value-blind. Endpoint +4 NIE wbudowany. 0 hardcoded T_pass. Werdykt (B) #49 LOCKED."
---

# Phase 0 (szkic) — mini-test R1

> **Cel bramki R1:** rozstrzygnąć PRZED startem, czy „crossover wykładnika kinetycznego
> `−½ → +2`" niesie treść fizyczną, czy jest tylko przepisaniem tej samej teorii w innej
> zmiennej (tautologia). Jeśli tautologia bez ucieczki → cała teza dwufazowa pada (uczciwy
> NEGATIVE, oszczędza miesiące). Jeśli istnieje ucieczka → wiemy DOKŁADNIE, jakiego obiektu
> trzeba pilnować (nie samego wykładnika).

## §1 — Pytanie operacyjne

Substrat: kanoniczna kinetyka amplitudy `(∂σ)²`. Gęstość makro: `Φ = F(σ)`. Po zamianie
zmiennych `(∂σ)² = K(Φ)(∂Φ)²`, `K(Φ) = (dσ/dΦ)²`. Dla `Φ=σ²` ⟹ `K∝Φ⁻¹` (= wynik #49).

**Pytanie:** czy wykładnik `e` w `K∝Φ^e` jest obiektem fizycznym, czy konwencją zmiennej?
I jeśli konwencją — co (jaki dodatkowy obiekt) czyni go fizycznym?

## §2 — Cztery pod-testy symboliczne (sympy, value-blind)

| Test | Pyta | Metoda | Co rozstrzyga |
|---|---|---|---|
| **T1** | Czy gołe `K(Φ)` pojedynczego pola jest konwencją? | kanonizacja `χ=∫√K dΦ`; czy `χ ∝ σ` dla dowolnego `p` | jeśli tak → R1 realne (tautologia dla pola swobodnego) |
| **T2** | Czy „dressing metryczny" (`√−g∝Φ^m`) mostkuje `−1→+4`? | net `E=m+e`; solve `m` dla `E=+4`; porównaj z kandydatami `m∈{1,2}` | czy naiwna Droga B wystarcza |
| **T3** | Co przywraca treść fizyczną? | z potencjałem `V`: niezmiennik `m²_phys=V''/(2K)` w obu ramach | identyfikuje, że anchor = para `(K,V)`, nie `K` samo |
| **T4** | Czy skalowo-zależna mapa `F_ℓ` daje nietautologiczny crossover? | jakobian `F'`; gdzie traci odwracalność | gdzie crossover MOŻE być fizyczny (Φ=0?) |

## §3 — Reguła interpretacji (DIAGNOSTYCZNA, nie falsyfikująca)

- T1 PASS (χ∝σ) ⟹ **R1 potwierdzone realne**: gołego wykładnika nie wolno traktować jako fizyki.
- T3 pokazuje niezmiennik ⟹ **ucieczka z R1 istnieje**, ale przenosi cel z „wykładnika" na „parę (K,V)".
- T2: jeśli wymagane `m` ≠ kandydaci TGP ⟹ naiwna Droga B niewystarczająca (potrzeba dokładnej akcji sek08c — Phase 1).
- T4: lokus nieodwracalności wskazuje, gdzie jedynie crossover może być fizyczny (kandydat: Φ=0 / nukleacja).

**To NIE jest werdykt o teorii.** To mapa: które wersje tezy dwufazowej są martwe (tautologia),
a która jest żywa (i czego wymaga).

## §4 — Forbidden (mini-test)
1. Zakaz re-litygacji #49 (Droga A / wymiar anomalny poza zakresem).
2. Zakaz wbudowania `+4` w `K_eff` (endpoint wyliczany/porównywany, nie zakładany).
3. 0 nowych stałych; 0 danych obserwacyjnych; 0 hardcoded `T_pass`.

Wynik: [[R1_test.txt]] · skrypt: [[R1_test.py]]
