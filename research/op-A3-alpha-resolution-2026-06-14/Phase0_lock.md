---
title: "Phase 0 — op-A3-alpha-resolution: LOCK pytania i reguły decyzyjnej. Rozstrzygnąć niespójność α=2↔K(φ) wykrytą w op-CG34-continuum-closure (lemat A3): czy α=2 (człon (∇Φ)²/Φ w równaniu pola TGP) WYNIKA z geometrycznego sprzężenia substratu przez Φ=φ², czy jest niespójnością konwencji micro/macro. Value-blind, anti-Lakatos."
type: phase_lock
status: PHASE0_LOCKED
phase: 0
cycle: op-A3-alpha-resolution-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14 (sesja #31): 'tak działaj z op-A3-alpha-resolution'"
origin: "Residual op-CG34-continuum-closure Phase 3: α_eff=s−1 (Z∝φ^{2s}); substrat s=1→α=0, NIE 2"
target: "sek08 thm:D-uniqueness (α=2), sek10 (K(φ)=K_geo·φ⁴, eq:kinetic_macro), dodatekQ2 lemat A3, ontologia Φ=⟨ŝ²⟩"
anti_lakatos_lock: ACTIVATED
---

# Phase 0 — LOCK (op-A3-alpha-resolution)

> **Pytanie wiodące:** Czy $\alpha=2$ (współczynnik członu $(\nabla\Phi)^2/\Phi$ w równaniu pola TGP)
> **wynika** z geometrycznego sprzężenia kinetycznego substratu ($K\propto\varphi^4$) przez identyfikację
> $\Phi=\langle\hat s^2\rangle\propto\varphi^2$ — czy jest **niespójnością konwencji** (micro $\varphi$ vs macro $\Phi/\Phi_0$)?

## §1 — Materiał dowodowy (z rdzenia, zlokalizowany)
| Źródło | Treść | Wariabl |
|---|---|---|
| **sek08 thm:D-uniqueness** (eq:EL-general, eq:K-ode) | $\mathcal L_{\rm kin}=-\tfrac12K(\varphi)(\nabla\varphi)^2$, $\varphi=\Phi/\Phi_0$; EL coeff $\alpha=\tfrac{K'\varphi}{2K}$; $K=C\varphi^{2\alpha}$. C3: $K=K_{\rm geo}\varphi^4\Rightarrow2\alpha=4\Rightarrow\alpha=2$ | **macro** $\varphi=\Phi/\Phi_0$ |
| **sek10 lemat K_φ2** (eq linia 82) | $H_\Gamma$ gradient expansion: $K(\varphi)=Ja^2\varphi^2$ | **micro** spiny $\varphi_i$ |
| **sek10 prop:substrate-action** | $K_{ij}=J(\varphi_i\varphi_j)^2$ → $K(\varphi)=K_{\rm geo}\varphi^4$ | **micro** spiny |
| **sek10 eq:Phi_phi_id** (linia 104) | $\Phi=(\varphi^2/\varphi_{\rm ref}^2)\Phi_0$ ⟹ $\Phi\propto\varphi^2$ (ontologia $\Phi=\langle\hat s^2\rangle$) | most micro↔macro |
| **sek10 eq:kinetic_macro** (linie 114–119) | transformacja $K_{\rm geo}\varphi^4(\nabla\varphi)^2 = K_0\,\psi(\nabla\psi)^2$ (**liniowy** w $\psi=\Phi/\Phi_0$) | poprawna zmiana zmiennych |
| **dodatekQ2 lemat A3** | premisa $Z\sim\varphi^2$, konkluzja $K_1\sim1/\Phi$ → „α=2" | niespójny (op-CG34 Phase 3) |

**Podejrzenie (do rozstrzygnięcia, NIE przesądzone):** thm:D-uniqueness C3 wstawia $K_{\rm geo}\varphi^4$ do
**macro** $\varphi=\Phi/\Phi_0$, podczas gdy substrat daje $K\propto\varphi_{\rm micro}^4$, a $\Phi\propto\varphi_{\rm micro}^2$ — więc po zmianie zmiennych (sek10 eq:kinetic_macro) $K_{\rm eff}$ jest **liniowy** w $\psi$, nie $\psi^4$.

## §2 — Obiekt i konwencja LOCKED
- Równanie pola TGP: $\nabla^2\Phi+\alpha\dfrac{(\nabla\Phi)^2}{\Phi}+\ldots=0$ (sek08 eq:EL-general, eq:Dkin-unique) — **α = bezpośredni współczynnik** członu $(\nabla\Phi)^2/\Phi$.
- Relacja LOCKED (sek08, poprawna): dla $\mathcal L=-\tfrac12 K(\varphi)(\nabla\varphi)^2$, $K(\varphi)=C\varphi^{2\alpha}$ ⟺ EL coeff $=\alpha$. (Tę relację weryfikujemy, NIE podważamy.)
- Ontologia LOCKED: $\Phi=\langle\hat s^2\rangle$ ⟹ $\Phi\propto\varphi_{\rm micro}^2$ (sek08 poziom 0, sek10 eq:Phi_phi_id).

## §3 — Falsyfikatory / reguła decyzyjna LOCKED
| ID | Test | Reguła |
|---|---|---|
| **F-α-A** | Relacja EL: $K\propto\varphi^{2\alpha}\Leftrightarrow$ coeff $=\alpha$ (macro var) | sympy — musi PASS (sanity; jeśli FAIL, błąd po mojej stronie) |
| **F-α-B** | Transformacja substratu: $K(\varphi_{\rm micro})\propto\varphi_{\rm micro}^{2s}$ + $\Phi\propto\varphi_{\rm micro}^2$ → $\alpha=?$ | sympy wylicza $\alpha$ value-blind |
| **F-α-C** | Czy $\alpha=2$ z substratu? | $\alpha(\text{substrat})=2$ → **DERIVED** (α=2 potwierdzone); $\neq2$ → **INCONSISTENCY** (chain złamany) |
| **F-α-D** | Co potrzeba dla α=2? | wyznaczyć wymaganą potęgę $K(\varphi_{\rm micro})$ lub relację $\Phi(\varphi)$ |

**Reguła agregatu:** F-α-C = DERIVED jeśli poprawnie stransformowany substrat daje α=2; INCONSISTENCY jeśli ≠2 (z jawnym wskazaniem gdzie chain się łamie i jaka korekta by go naprawiła).

## §4 — Forbidden moves (anti-Lakatos)
1. **Zakaz przesądzania** α=2 ani α≠2 — sympy liczy value-blind.
2. **Zakaz rewizji rdzenia** (thm:D-uniqueness, sek10) — werdykt = rekomendacja, nie edycja core.
3. **Obie strony** jawnie: argument za α=2 (MC LK-1 „1.2σ", uniqueness) i przeciw (zmiana zmiennych).
4. Relacja EL ($\alpha=K'\varphi/2K$) jest rdzeniowa i POPRAWNA — nie podważam jej; sprawdzam tylko C3 (czy $K(\Phi/\Phi_0)\propto(\Phi/\Phi_0)^4$).
5. Jeśli INCONSISTENCY — podać **konstruktywne opcje naprawy** (Φ∝φ amplitude? K∝φ^{10}? inna konwencja?).

## §5 — Oczekiwanie (INFORMATIONAL, nie wiąże)
Z op-CG34 Phase 3: prawdopodobnie **INCONSISTENCY** (substrat poprawnie stransformowany daje α=0 lub 1/2,
nie 2; thm:D-uniqueness C3 myli micro/macro). Ale: MC LK-1 podobno mierzy α=2 (1.2σ, słabo). Werdykt rozstrzygnie sympy.

## §6 — Status
**🔒 PHASE 0 LOCKED.** Anti-Lakatos aktywny. Przejście do derywacji sympy (value-blind).
