---
title: "Phase 0 — op-Csigma-coarse-graining: LOCK falsyfikatorów F-CG-A..E + próg survival κ_E=5/6 (value-blind, PRZED rachunkiem). Obiekt: C_σ = sztywność kinetyczna σ_ab (eq:S-sigma); κ_E=C_σσ_0²/(GR-lock)."
date: 2026-06-14
type: phase-balance
phase: 0
status: 🔒 LOCKED 2026-06-14
cycle: op-Csigma-coarse-graining-2026-06-14
authorization: "User 2026-06-14 (sesja #31): autoryzacja fazy 1 = 'działaj' (Phase 0 LOCK → Phase 1)"
object_id: "C_σ = stała kinetyczna σ_ab (eq:S-sigma); σ_0 = normalizacja (eq:h-TT-from-sigma); κ_E = C_σσ_0²/(wartość GR-lock)"
scoping: "[[../../meta/SCOPING_op-Csigma-coarse-graining_2026-06-14.md]]"
parent_verdicts_LOCKED:
  - "op-sigma-kinetic-Csigma: sektor radiacyjny UNDERDETERMINED-fine-tuned; survival ⟺ κ_E=5/6 (miara zero); F-CS-A=GAP (C_σ niewyprowadzone)"
  - "PR-025: branże konforemne FALSIFIED (gałąź B = 7/6 → 2646σ); LOCKED"
  - "thm:amplitude-matching: ξ_eff=4πG_0σ_0Φ_0 (amplituda dopasowana do GR); m_σ→0 przyjęte"
  - "rem:sigma-params: C_σ (lub σ_0) = problem otwarty rdzenia ('wyznaczalny w zasadzie, obecnie niezobliczony')"
anti_lakatos_lock: "INHERITED + ACTIVE od tej fazy"
tags: [pre-registration, C-sigma, coarse-graining, kappa-E-pinning, value-blind-lock, anti-Lakatos-LOCKED]
---

# Phase 0 — LOCK (op-Csigma-coarse-graining)

> **Cel Phase 0:** zamrozić falsyfikatory i progi PRZED jakimkolwiek rachunkiem coarse-grainingu,
> by werdykt Phase FINAL był value-blind. Zero wartości liczbowych C_σ tutaj — tylko reguły decyzji.

## §1 — Pytanie (LOCKED)
Czy stała kinetyczna $C_\sigma$ w
$$S_\sigma=\int d^4x\sqrt{-g_{\rm eff}}\Big[\tfrac{C_\sigma}{2}\,\partial_\mu\sigma_{ab}\partial^\mu\sigma^{ab}-\tfrac{m_\sigma^2}{2}\sigma_{ab}\sigma^{ab}+\tfrac{\xi_{\rm eff}}{\Phi_0^2}\sigma_{ab}\partial^a\Phi\partial^b\Phi\Big]\quad(\text{eq:S-sigma})$$
jest **wyprowadzalna przez coarse-graining** korelacji kierunkowej
$\sigma_{ab}=K_{ab}-\tfrac13\mathrm{Tr}K\,\delta_{ab}$, $K_{ab}=\langle\hat s_i\hat s_{i+\hat a_b}\rangle$ (def:sigma-ab)
z hamiltonianu substratu $H_\Gamma$ (eq:B-H) — a jeśli tak, czy złożona w
$\kappa_E=C_\sigma\sigma_0^2/(\text{GR-lock})$ daje $5/6$ (SURVIVE) czy $\neq5/6$ (FALSIFIED)?

## §2 — Wejścia LOCKED (z rdzenia + rodzica; NIE rewidowane w tym cyklu)
| Input | Forma | Źródło |
|---|---|---|
| $H_\Gamma$ (substrat) | $\sum_i[\tfrac{m_0^2}{2}\hat s_i^2+\tfrac{\lambda_0}{4}\hat s_i^4]-J\sum_{\langle ij\rangle}\hat s_i\hat s_j$ (v1); GL-bond v2 $J\sum A_{ij}\hat s_i^2\hat s_j^2(\hat s_j^2-\hat s_i^2)^2$ | eq:B-H, eq:H-Gamma |
| Obiekt σ_ab | $\sigma_{ab}=K_{ab}-\tfrac13\mathrm{Tr}K\,\delta$; sym., bezśladowy, 5 komp. (spin-2) | def:sigma-ab |
| Akcja σ_ab | eq:S-sigma (powyżej); $C_\sigma>0$ z $J>0$ (ghost-free) | eq:S-sigma, prob:tensor-modes |
| Amplituda (pinned) | $\xi_{\rm eff}=4\pi G_0\sigma_0\Phi_0$; $m_\sigma\to0$ | thm:amplitude-matching, rem:sigma-params |
| Normalizacja | $h^{TT}_{ij}=2\sigma_{ij}/\sigma_0$ | eq:h-TT-from-sigma |
| Strumień/survival | $\dot P_b=\kappa_E P_{GR}+\tfrac16 P_{GR}$; survival ⟺ $\kappa_E=5/6$; $\kappa_E\approx1\to7/6$=gałąź B | op-sigma-kinetic-Csigma Phase 1 |
| Prędkość | σ_ab propaguje z $c_0$ (Box=$c_0^{-2}\partial_t^2-\nabla^2$) | rem:cGW-tensor |
| Parametry substratu | $\{J,\,a_{\rm sub},\,\mu,\,m_0^2,\,\lambda_0\}$ | eq:B-H, rem:param-counting |

## §3 — Falsyfikatory LOCKED (próg pre-rachunkowy)
| ID | Test | Reguła decyzji (LOCKED) |
|---|---|---|
| **F-CG-A** | Czy kinetyka σ_ab emerguje **Lorentz-inwariantnie** z $H_\Gamma$ ($C_\sigma>0$, propagacja $c_0$, człon czasowy $(\partial_t\sigma)^2$)? | TAK → przejdź do B; NIE (brak emergentnego członu czasowego/Lorentza) → **GAP strukturalny** |
| **F-CG-B** | Forma $C_\sigma(\{J,a_{\rm sub},\dots\})$ z gradient expansion | **DERIVED** (jawna formuła z prefaktorem) / **PARTIAL** (skaling, prefaktor O(1) otwarty) / **GAP** |
| **F-CG-C** | $\sigma_0$ i złożenie $\kappa_E=C_\sigma\sigma_0^2/(\text{GR-lock})$; czy $C_\sigma$ i $\sigma_0$ to 1 czy 2 parametry | wartość $\kappa_E$ lub przedział |
| **F-CG-D** | $\kappa_E$ vs survival **5/6** | $\kappa_E=5/6\pm$tol → **SURVIVE**; $\neq$ → **FALSIFIED**; tylko skaling → **PARTIAL-lean** |
| **F-CG-E** | Agregat | DERIVED ∧ κ_E≠5/6 → **FALSIFIED**; ∧ κ_E=5/6 → **SURVIVE**; GAP/PARTIAL → sektor pozostaje **UNDERDETERMINED-fine-tuned** (status węższy, lean wg skalowania) |

**Próg tolerancji (LOCKED):** survival wymaga $\kappa_E=5/6$ EXACT (miara zero). „±tol" oznacza: jeśli
gradient expansion da $\kappa_E$ z prefaktorem znanym co do rzędu, to wartość $5/6$ musi leżeć w
przedziale skalowania, by uznać PARTIAL-lean SURVIVE; jeśli skaling wyklucza $5/6$ o rzędy →
**de facto FALSIFIED**. Żaden swobodny parametr NIE może być dostrojony, by trafić w $5/6$ (§5 pkt 1).

## §4 — Pre-derywacja (oczekiwanie; INFORMATIONAL, rachunek nadrzędny)
Najbardziej prawdopodobne **PARTIAL/GAP** (R-tensor-hard HIGH): coarse-graining sektora **bezśladowego
tensorowego** jest trudniejszy niż skalarnego (op-gamma-RG-running domknął skalar tylko numerycznie
LPA', analitycznie OTWARTY — dodatekQ CG-1/CG-3/CG-4 OTWARTE). Oczekiwany skaling $C_\sigma\sim J\cdot f(a_{\rm sub})$
bez prefaktora O(1). Drugi scenariusz: **DERIVED → κ_E≠5/6 → FALSIFIED** (najmocniejszy, definitywne domknięcie).
Trzeci (mało prawdopodobny): κ_E=5/6 SURVIVE (wymaga „spisku" dwóch sektorów — R-conspiracy).

## §5 — Forbidden moves (LOCKED)
1. **Zakaz strojenia $C_\sigma/\sigma_0$ do $\kappa_E=5/6$** — survival musi WYNIKNĄĆ z coarse-grainingu, nie być dobrane (value-blind).
2. Zakaz mylenia obiektu: $C_\sigma$ = sztywność σ_ab (nośnik GW), NIE skalar / NIE g_eff / NIE induced-TT (lekcja 3 korekt sesji #30).
3. Zakaz rewizji PR-025 / parent / rem:cGW-tensor / thm:amplitude-matching (LOCKED).
4. **GAP/PARTIAL deklarowane JAWNIE** jeśli coarse-graining tensorowy przekracza obecną machinerię — zakaz fabrykowania prefaktora O(1) (anti-Lakatos, lekcja #30; rodzic zrobił to poprawnie = F-CS-A GAP).
5. Budżet nowych stałych **0** ($C_\sigma$ istnieje w eq:S-sigma; $\{J,a_{\rm sub},\mu,m_0^2,\lambda_0\}$ to parametry substratu).
6. **Wymóg Lorentza** ($c_0$) na kinetykę σ_ab — nie naruszać (rem:cGW-tensor).

## §6 — Anti-Lakatos (checklist Phase 0)
✓ Obiekt zidentyfikowany PRZED rachunkiem ($C_\sigma$ = sztywność σ_ab; §1–§2).
✓ Falsyfikatory + próg 5/6 LOCKED value-blind PRZED Phase 1 (§3).
✓ GAP/PARTIAL dopuszczalne i a priori oczekiwane (§4) — brak ochrony frameworku.
✓ Zakaz strojenia survival (§5 pkt 1).
✓ Liczby poprzedników LOCKED (§2); budżet stałych 0 (§5 pkt 5).
