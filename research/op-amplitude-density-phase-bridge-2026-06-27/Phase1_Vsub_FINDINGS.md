---
title: "FINDINGS Phase 1 first-pass — V_sub ODNALEZIONY; alpha=2 = aksjomat v2 (sprzezenie geometryczne), NIE derywacja z mikro H_Gamma; teza dwufazowa usera RATYFIKOWANA jako poprawna interpretacja"
date: 2026-06-27
type: phase1-findings
status: "DIAGNOSTIC RESULT (first-pass; nie domkniecie cyklu; bramka V_sub CLEARED)"
parent: "[[Phase0_R1_FINDINGS.md]]"
scoping: "[[../../meta/SCOPING_op-amplitude-density-phase-bridge_2026-06-27.md]]"
artifacts: ["[[Phase1_Vsub_sharp.py]] / [[Phase1_Vsub_sharp.txt]]"]
core_sources:
  - "axioms/substrat/dodatekB_substrat.tex (eq:B-H, eq:K-geometric, eq:U-GL, rem:B-v2-status, prop:substrate-action)"
  - "core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex (K_sub(g)=g^2, H_Gamma != F_kin, g=phi^2=Phi/Phi0)"
locked_inputs: ["#49 (alpha_eff=-1/2), IMMUTABLE", "status_map l.72 (alpha=2 selekcja na gestosci)"]
anti_lakatos_note: "Wynik RATYFIKUJE szczery samoopis rdzenia (rem:B-v2-status: alpha=2 = 'bezposredni skutek aksjomatu, bez wyprowadzenia MK->GL'). Zero nowych claimow o teorii. 0 hardcoded. Endpoint +4 nie wbudowany (T_B wyliczony z K_ij)."
---

# FINDINGS Phase 1 first-pass — V_sub + status alpha=2

## §0 — Wynik w jednym zdaniu

**`V_sub` istnieje jawnie w rdzeniu** (potencjał GL substratu); ostry test pokazuje, że
**`α=2` NIE wynika z mikroskopowego `H_Γ` (bilinearny bond → kinetyka kanoniczna → `K∝Φ⁻¹`,
dokładnie #49), lecz z osobnego aksjomatu v2 — sprzężenia geometrycznego `K_ij=J(φ_iφ_j)²`** —
co rdzeń sam deklaruje jako „pivot aksjomatyczny". To **ratyfikuje #49 od strony analitycznej**
i ustawia tezę dwufazową usera jako **poprawną interpretację**: `−½` żyje na poziomie mikro
(`ŝ`, `H_Γ`), `+2` na poziomie aksjomatu v2 (`φ`); luka `−½→+2` = **status aksjomatu
geometrycznego (postulat vs derywacja)**.

## §1 — V_sub ZNALEZIONY (bramka egzystencjalna #0 CLEARED)

Potencjał substratu istnieje w rdzeniu w trzech, spójnych warstwach:

| Warstwa | Postać | Zmienna | Źródło |
|---|---|---|---|
| Mikro on-site | `V_ŝ = (m₀²/2)ŝ² + (λ₀/4)ŝ⁴` | `ŝ` (spin substratu, Z₂) | eq:B-H |
| Landau | `V_GL(v) = (r/2)v² + (u/4)v⁴`, `v=⟨ŝ⟩` | `v` | eq:VGL-full |
| GL amplitudowy | `U(φ) = (β/3)φ³ − (γ/4)φ⁴` `=` `−γφ³(3φ−4)/12` (przy β=γ) | `φ=(Φ/Φ₀)^{1/2}` | eq:U-GL |

Słownik: `φ² = Φ/Φ₀ = g`, `Φ=⟨ŝ²⟩` (parzyste w `ŝ`). To było **brakujące wejście** ostrego
testu R1 (Phase0 FINDINGS §3) — teraz dostępne.

## §2 — Ostry test (sympy, exact) — wyniki

| Test | Wynik | Znaczenie |
|---|---|---|
| **T_A** | mikro bond `−Jŝ_iŝ_j` → prefaktor kinetyczny `K_s=J/2` (stały, `e_s=0`) | **kanoniczna kinetyka w `ŝ`** → w `Φ=ŝ²`: `K∝Φ⁻¹` = **#49** |
| **T_B** | geom v2 `K_ij=J(φ_iφ_j)²` na `(φ_i−φ_j)²` → `K(φ)=Jφ⁴`, `e_φ=4` | `α=e/2=2` (manuskrypt) — wyliczone, nie wbudowane |
| **T_C** | `(φ_iφ_j)² ~ ŝ_i²ŝ_j²` ≠ `ŝ_iŝ_j` (bilinear) | **kwartyczny `K_ij` NIE jest w `H_Γ`** → `α=2` to osobny obiekt |
| **T_D** | `dΦ/dŝ=2ŝ` → 0 w `ŝ=0`, sign-blind (`ŝ`,`−ŝ`→ ta sama `Φ`) | mapa `ŝ→⟨ŝ²⟩` = **coarse-graining, NIE redefinicja** → R1 omijalne, ale to droga #49 (`−½`) |

## §3 — Kluczowe rozstrzygnięcie (i dlaczego NIE jest sprzecznością z rdzeniem)

Rdzeń (`rem:B-v2-status`, 2026-04-24) deklaruje wprost:

> „Wynik `α=2` jest w v2 **bezpośrednim skutkiem aksjomatu** poprzez prop:substrate-action,
> **bez konieczności wyprowadzenia MK→GL**."

Czyli `α=2` **jest** w rdzeniu aksjomatem (sprzężenie geometryczne v2 = „pivot aksjomatyczny"),
a NIE derywacją z mikro `H_Γ`. Mój T_C to **potwierdza analitycznie** (kwartyczny `K_ij` ≠
bilinearny bond), zgodnie z samoopisem rdzenia i z `status_map` l.72 („selekcja na gęstości,
NIE derywacja z substratu"). **Zero nowej sprzeczności; ratyfikacja istniejącego, uczciwego statusu.**

## §4 — Co to znaczy dla tezy dwufazowej usera

**Teza usera = POPRAWNA INTERPRETACJA luki, nie jej zamknięcie:**

- `α_eff=−½` ⟺ poziom **mikro** (`ŝ`, `H_Γ` bilinearny, near-Gaussian) — reżim „amplitudowy/przedaksjomatyczny".
- `α=2` ⟺ poziom **aksjomatu v2** (`φ`, sprzężenie geometryczne `K_ij=J(φ_iφ_j)²`) — reżim „gęstościowy/postulowany".
- Mapa między nimi (`ŝ→Φ=⟨ŝ²⟩`) to **prawdziwy coarse-graining** (T_D: nieodwracalny), więc
  crossover **NIE jest tautologią** (R1 ostatecznie omijalne tą drogą). **ALE** ten właśnie
  coarse-graining to droga #49 → daje `−½`, nie `+2`.

⟹ **„Most" którego szuka teza = wyprowadzenie aksjomatu geometrycznego v2 (`K_ij=J(φ_iφ_j)²`)
z mikro `H_Γ`.** #49 mówi: near-Gaussian unitarny coarse-graining tego NIE robi (daje `−½`,
`η~O(5)` zamyka ucieczkę). Pozostaje wąskie, fundamentalne pytanie na PRAWDZIWY cykl:
**czy istnieje JAKIKOLWIEK (non-Gaussian / strong-coupling) coarse-graining `H_Γ→F_kin`,
który generuje kwartyczne `K_ij`?** Jeśli nie — `α=2` jest trwałym aksjomatem (uczciwie).

## §5 — Krajobraz trzech wykładników (zaostrzony R5)

Test ujawnił, że w grze są **trzy** wartości, każda z innym źródłem:

| Wykładnik | `K∝` | Źródło | Status |
|---|---|---|---|
| `α_eff=−½` | `Φ⁻¹` | mikro `H_Γ` (ŝ-kanoniczny) + `Φ=⟨ŝ²⟩` | #49 (derywowane, LOCKED) |
| `α=1` | `g²` (=`φ²` w `g`) ↔ `φ²` | `H_Γ` jako energia bondu (`H_Γ≠F_kin`) | sek08b (sektor solitonu) |
| `α=2` | `φ⁴` | aksjomat v2 `K_ij=J(φ_iφ_j)²` | manuskrypt (selekcja) |

R5 (α=1 vs α=2) jest więc głębsze niż „dwie fazy" — to **trzy różne funkcjonały kinetyczne**.
Prawdziwy cykl musi wyjaśnić, czemu różne sektory siadają na różnych `α` (lub zredukować je).

## §6 — Czego NIE wolno wywnioskować
- NIE: „obalono α=2" / „znaleziono sprzeczność w rdzeniu". Rdzeń jest szczery (rem:B-v2-status).
- NIE: „most CG domknięty" / „α=2 zderywowane". `N_free` bez zmian; #49, status_map l.72 LOCKED.
- NIE: „teza usera udowodniona". Pokazano tylko, że jest **dobrze postawiona** i wskazano
  jej dokładny przedmiot (status aksjomatu geometrycznego v2).

## §7 — Zaktualizowane bramki przed PRAWDZIWYM Phase 0 cyklu
- ✅ **#0 V_sub** — CLEARED (istnieje; §1).
- ✅ **R1 tautologia** — rozbrojone (T_D: mapa to CG, nie redef; ale = droga #49).
- 🔜 **Pytanie rdzenia cyklu (NOWE, ostre):** czy `H_Γ → K_ij=J(φ_iφ_j)²` ma JAKĄKOLWIEK
  drogę coarse-grainingu (non-Gaussian)? — to jest właściwy, wąski temat cyklu.
- 🔜 **R3** (emergent-metric niezależna od α=2), **R4** (sektor mas), **R5** (teraz: 3 wykładniki).

## Cross-references
- [[Phase0_R1_FINDINGS.md]] · [[Phase1_Vsub_sharp.py]] · [[Phase1_Vsub_sharp.txt]]
- [[../op-CG-alpha-eff-convergence-2026-06-26/Phase_FINAL_close.md]] (#49, LOCKED)
- [[../../meta/HONEST_FRAMING_UV_CG_ROOTS.md]] (cztery korzenie)
- axioms/substrat/dodatekB_substrat.tex (eq:B-H, eq:K-geometric, eq:U-GL, rem:B-v2-status)
- core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex (H_Gamma != F_kin)
