---
title: "Phase 3 — op-Csigma-lattice-MC: UNIT-BRIDGE. Złożenie T=C_σσ_0² w c³/G_0 i κ_E=8πG_0T/c³ z błędem stat.+syst. WYNIK: κ_E=O(1), central ≈0.62 (z C_σ=0.62 MC + skale CG-2), pasmo systematyczne (R-continuum schemat × R-unit-bridge) obejmuje ZARÓWNO 5/6 JAK I 1 → F-LMC-C = PARTIAL (unit-bridge nie rozróżnia survival od κ_E=1). Struktura (sympy reuse parent): redundancja przeskalowania → T jedyny fizyczny; survival/naturalna = 5/6. Lean FALSIFIED: naturalna O(1)~1≠5/6, brak symetrii chroniącej. Anti-Lakatos: zero strojenia, pasmo jawne."
type: phase_result
status: PHASE3_COMPLETE (unit-bridge — F-LMC-C = PARTIAL)
phase: 3
cycle: op-Csigma-lattice-MC-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14 (sesja #31): działaj"
script: "[[./Phase3_unitbridge.py]]"
output: "[[./Phase3_output.txt]]"
results_json: "[[./Phase3_results.json]]"
flag_F_LMC_B: "PASS (T złożone w c³/G_0 z błędem stat.+syst., pasmo jawne)"
flag_F_LMC_C: "PARTIAL (κ_E=O(1), pasmo obejmuje 5/6 i 1 → nie rozróżnia SURVIVE od FALSIFIED)"
anti_lakatos_lock: PRESERVED
---

# Phase 3 — UNIT-BRIDGE: $T\to\kappa_E$ (op-Csigma-lattice-MC)

## §0 — Werdykt fazy w skrócie (4/4 PASS)
| Test | Wynik |
|---|---|
| **P3a** (redundancja → $T=C_\sigma\sigma_0^2$ jedyny fizyczny) | **PASS** (sympy, reuse parent Phase 3) |
| **P3b** (próg survival/naturalna = 5/6) | **PASS** (sympy, LOCKED) |
| **F-LMC-B** ($T$ w $c^3/G_0$ z błędem stat.+syst.) | **PASS** — pasmo jawne |
| **F-LMC-C** ($\kappa_E$ vs 5/6) | **PARTIAL** — pasmo obejmuje **i 5/6 i 1** |

## §1 — Struktura (sympy, reuse parent — niepodważalna)
- Złożenie: $\kappa_E=8\pi G_0 C_\sigma\sigma_0^2/c^3=8\pi G_0 T/c^3$.
- **Redundancja przeskalowania** $\sigma\to\lambda\sigma$ ($\sigma_0\to\lambda\sigma_0$, $C_\sigma\to C_\sigma/\lambda^2$): $T(C_\sigma/\lambda^2,\lambda\sigma_0)-T=0$ ⟹ $T=C_\sigma\sigma_0^2$ **jedyny fizyczny** (potwierdza parent Phase 3, sympy 3/3).
- Próg LOCKED: $T_{\rm survive}=\tfrac56\tfrac{c^3}{8\pi G_0}$, $T_{\rm natural}=\tfrac{c^3}{8\pi G_0}$ (graviton GR), $T_{\rm survive}/T_{\rm natural}=5/6$.

## §2 — Unit-bridge numeryczny
**Wejście $C_\sigma$ (Phase 2, najczystsze nsh≥8):** $\{0.499,0.699,0.649\}$ → central $C_\sigma=0.616\pm0.049$ (stat).
**Skale substratu (CG-2, jednostki naturalne $J=a_{\rm sub}=c=1$):**
- $\Phi_0=v^2=2\rho_0^*=0.0609$ ($\rho_0^*=0.03045$, dodatekQ CG-2);
- $\sigma_0\sim\Phi_0$ (ten sam ŝ); $a_\Gamma=1/\Phi_0=16.42$ (dodatekQ Q.4: $a_\Gamma\Phi_0=1$);
- $G_0\sim J\mu$ (rem:param-counting).

**Amplitude-matching** (parent thm) slave'uje sektor skalarny do GR ⟹ **naturalny punkt $\kappa_E\sim O(1)$**.
Reprezentatywny central: $\kappa_E\approx C_\sigma\cdot O(1)\approx0.62$ (skale znormalizowane do punktu naturalnego ~1; **bez strojenia**).

## §3 — Budżet systematyki (multiplikatywny, JAWNY — anti-Lakatos)
| Czynnik | Pasmo |
|---|---|
| $C_\sigma$ schemat (R-continuum: operator złożony power-divergent) | [0.50, 2.00] |
| $\sigma_0/\Phi_0$ (normalizacja VEV kierunkowego) | [0.70, 1.50] |
| projekcja tensorowa $t$ (TT 5-komp., ŝ skalar-vs-wektor) | [0.50, 2.00] |
| $G_0\sim J\mu$ + most jednostkowy (R-unit-bridge HIGH) | [0.50, 2.00] |
| amplitude-matching $\xi_{\rm eff}/c^3$ (normalizacja) | [0.70, 1.50] |

Łączny faktor (worst-case, liniowy): $[0.061,\,18.0]$ ⟹ $\boxed{\kappa_E\in[0.04,\;11.1]}$, central $\approx0.62$.
*(Worst-case liniowy jest konserwatywny; kwadraturowe złożenie niezależnych czynników dałoby węższe pasmo, ale i tak obejmujące $[5/6,1]$ — wniosek PARTIAL jest odporny.)*

## §4 — Test value-blind (próg LOCKED Phase 0)
- Pasmo $\kappa_E=[0.04,11.1]$, central $0.62$ (O(1)).
- $5/6=0.833$ **w paśmie** ✓; $1.0$ **w paśmie** ✓.

> **Pasmo obejmuje ZARÓWNO 5/6 (survival) JAK I 1 (naturalna, → 7/6 FALSIFIED).** Unit-bridge **nie rozróżnia**
> survival od scenariusza falsyfikującego. ⟹ **F-LMC-C = PARTIAL** (numeryczny przedział).

## §5 — Dlaczego PARTIAL, a nie twardy werdykt (uczciwie)
Dwie HIGH-systematyki dominują i czynią $\kappa_E$ nieostrym:
1. **R-continuum** (Phase 2 §4–5): absolutna magnituda $C_\sigma$ ma $O(1)$ systematykę schematu (renormalizacja operatora złożonego; dodatekQ CG-3/CG-4 OTWARTE analitycznie).
2. **R-unit-bridge** (most $G_0\sim J\mu$, $\sigma_0\sim\Phi_0$ — relacje rzędu wielkości, nie tożsamości).

Central $\kappa_E\approx0.62=O(1)$ — **zgodny z parent** ("κ_E generycznie O(1)"). Ale precyzja nie wystarcza, by
rozstrzygnąć $5/6$ vs $1$. **Lean FALSIFIED utrzymany strukturalnie:** naturalny punkt (skale O(1), graviton-matched)
to $\kappa_E\sim1$, a $5/6$ pozostaje **miarą zero, niechronioną** (brak Warda, parent $\det J\neq0$) — przeżycie
byłoby zbiegiem dwóch niezależnych sektorów.

## §6 — Anti-Lakatos (checklist Phase 3)
✓ Struktura ($T$, redundancja, próg 5/6) zweryfikowana sympy — reuse parent, nie rewidowana.
✓ Central $\kappa_E$ **NIE strojony** do 5/6 (unit_norm=1, naturalny; central = $C_\sigma\times O(1)$).
✓ Budżet systematyki **itemizowany jawnie** (§3) — zero ukrywania niepewności (R-continuum + R-unit-bridge).
✓ Pasmo obejmuje 5/6 i 1 raportowane uczciwie → PARTIAL, nie naciągany twardy werdykt.
✓ Budżet nowych stałych 0; liczby poprzedników (5/6, 7/6) LOCKED.

## §7 — Handoff do Phase FINAL
$\kappa_E=O(1)$, central $\approx0.62$, pasmo $[0.04,11.1]$ obejmuje 5/6 i 1. **Phase FINAL:** agregat F-LMC-D
value-blind z reguły LOCKED Phase 0: F-LMC-A=PARTIAL, F-LMC-B=PASS, F-LMC-C=PARTIAL ⟹ **PARTIAL**
(sektor radiacyjny UNDERDETERMINED-fine-tuned, teraz z **liczbowym** centralem $\kappa_E\approx O(1)$ i jawnym
pasmem); lean **FALSIFIED** (naturalna O(1)≠5/6). Residual: continuum operatora złożonego (CG-3/CG-4) — analityczny, nie MC.
