---
title: "op-CG4-substrate-closure — domknięcie residuum CG-4: niepatologiczny model substratu dla scheme-independent C_σ (pinowanie κ_E → twardy werdykt sektora radiacyjnego GW)"
date: 2026-06-20
type: cycle_readme
status: 🟢 CLOSED (2026-06-20) — pełny cykl 0+1+2+3+FINAL; substrat RESOLVED (M0); C_σ dowiedzenie UV-czuły = wolny parametr (sektor radiacyjny nieusuwalny bąblem)
cycle: op-CG4-substrate-closure-2026-06-20
parent_cycles:
  - "[[../op-CG34-continuum-closure-2026-06-14/Phase_FINAL_close.md]] (CG-4 PARTIAL; residua: α=2↔K(φ) + N5 substrat)"
  - "[[../op-Csigma-lattice-MC-2026-06-14/Phase_FINAL_close.md]] (κ_E≈0.62, pasmo [0.04,11.1]; residual = scheme-indep. continuum operatora złożonego)"
  - "[[../op-Csigma-coarse-graining-2026-06-14/Phase_FINAL_close.md]] (C_σ>0, O(1); prefaktor = GAP)"
target: "dodatekQ CG-4 [CZĘŚCIOWY NUM] → [ZAMKNIĘTY NUM] lub GAP/FALSIFIED-hard; PLAN_NUMERYCZNY_CG3_CG4 N5"
anti_lakatos_lock: ACTIVATED
authorization: "User 2026-06-20: 'działaj z op-CG4-substrate-closure'"
tags: [CG-4, substrate-model, C-sigma, kappa-E-pinning, radiation-sector, R-continuum, composite-operator, anti-Lakatos]
---

# op-CG4-substrate-closure

> **Cel jednozdaniowy:** znaleźć i zwalidować model substratu $H_\Gamma$ (i) **stabilny**
> (bez runaway/frozen), (ii) z **czystym ciągłym punktem krytycznym Z₂** (kontrolowane $\xi\to\infty$),
> (iii) reprodukujący wymaganą **emergentną kinetykę** ($(\nabla\Phi)^2/\Phi$, $\Phi=\langle\hat s^2\rangle$),
> tak by **scheme-independent continuum** współczynnika $p^2$ operatora złożonego $O_{ab}=\partial\hat s\,\partial\hat s$
> dało $C_\sigma$ z dokładnością **< faktor 1.2** → zwężenie pasma $\kappa_E$ → **twardy werdykt** sektora
> radiacyjnego GW (SURVIVE / FALSIFIED-hard).

## §1 — Dlaczego ten cykl (ścieżka krytyczna)

Sektor grawitacyjny radiacyjny TGP jest **jedyną otwartą drogą ratunku** po falsyfikacji:
- metryka konforemna: FALSIFIED (LOCKED);
- M9.1'' $(4-3\psi)/\psi$ sygnatura GW: **FALSIFIED-OBSERVATIONAL 5.02σ** (GWTC-3, 2026-05-09);
- disformalny screening skalara: BROKEN (viability);
- **σ_ab (nośnik GW): UNDERDETERMINED-fine-tuned** — przeżycie ⟺ $\kappa_E=5/6$ dokładnie (miara zero).

Przeżycie ⟺ $\kappa_E=8\pi G_0 C_\sigma\sigma_0^2/c^3 = 5/6$ **dokładnie**; wartość naturalna $\kappa_E\sim1\Rightarrow$
suma 7/6 ⟹ FALSIFIED 2646σ. Pasmo lattice-MC $[0.04,11.1]$ (central 0.62) **nie rozróżnia** 5/6 od 1.

**Ostatnia brama** (lattice-MC FINAL §4 + CG-34 FINAL §4): scheme-independent continuum zrenormalizowanego
współczynnika $p^2$ operatora **złożonego**. Blokada strukturalna:
- (a) operator złożony jest **UV power-divergent** (R-continuum) → wymaga subtrakcji additywnej + mieszania operatorów;
- (b) testowany substrat $-J(\phi_i\phi_j)^2$ jest **patologiczny** (małe λ → runaway $\langle\rho\rangle\!\sim\!1600$; duże λ → frozen $\xi\to0$) ⟹ **brak czystego punktu krytycznego** na skończonej sieci ⟹ continuum niedostępne.

Residual jest **analityczny+numeryczny** (model substratu), nie statystyczny — większe MC nie pomoże bez lepszego modelu.

## §2 — Dwa residua CG-4 (z op-CG34 Phase_FINAL §4)

| # | Residual | Stan wejściowy | Rola w tym cyklu |
|---|---|---|---|
| **R-A** | **model substratu / N5** (niepatologiczny: stabilny + czysty punkt krytyczny + emergentne $(\nabla\Phi)^2/\Phi$) | zablokowany (patologia $-J(\phi_i\phi_j)^2$) | **GŁÓWNY** — warunek konieczny scheme-indep. $C_\sigma$ |
| **R-B** | **α=2 ↔ mikroskopowy $K(\phi)$** (niespójność lematu A3) | algebra $\alpha_{\rm eff}=s-1$; substrat $s\!=\!1\Rightarrow\alpha\!=\!0$ | **POBOCZNY** — rozstrzygnięty konceptualnie #32 (α=2 = **postulat na gęstości**, nie derywacja); tu: spójność wybranego modelu z `rem:alpha2-pivot-status` |

**Uwaga anti-Lakatos (z sesji #32):** substrat (amplituda $\sqrt\Phi$) daje α=½ (≠2). Cykl **NIE** próbuje
wyprowadzić α=2 z substratu — α=2 pozostaje **selekcją aksjomatyczną C1–C3 na gęstości**. R-B jest tu wyłącznie
testem **niesprzeczności** wybranego modelu (czy struktura $(\nabla\Phi)^2/\Phi$ jest generowana — nie czy jej
współczynnik = 2).

## §3 — Struktura faz (rekomendowana; każda = osobne „działaj")

| Faza | Treść | Metoda | Status |
|---|---|---|---|
| **0** | LOCK kryteriów C-A..C-D, falsyfikatory, reguła agregatu, forbidden moves; balance sheet | — | 🔒 **LOCKED** |
| **1** | Silnik analityczny: klasyfikacja kandydatów substratu (stabilność + punkt krytyczny), reprodukcja $\alpha_{\rm eff}=s-1$, znak $C_\sigma$ z bąbla 3D + power-counting operatora złożonego | sympy | ⏳ oczekuje |
| **2** | MC najlepszego kandydata: walidacja silnika (faza, $\xi$, pik $\chi$ w $\beta_c$) + pomiar $C_\sigma$ z structure factor operatora złożonego, skan rozmiaru → continuum | numeryka (3D MC) | — |
| **3** | Renormalizacja: subtrakcja additywna UV + mieszanie operatorów → **scheme-independent** $C_\sigma$; propagacja do $\kappa_E$ + pasmo | analityczny + num | — |
| **FINAL** | Agregat C-D (value-blind, reguła LOCKED) → CG-4 ZAMKNIĘTY NUM / GAP / FALSIFIED-hard; rekomendacje rdzenia | — | — |

## §4 — Pliki

- [[Phase0_lock.md]] — kryteria zalockowane (PASS/FAIL, falsyfikatory, agregat, forbidden moves)
- [[Phase0_balance.md]] — **obowiązkowy** balance sheet (CALIBRATION_GATE_ENFORCEMENT)
- (Phase1_*, Phase2_*, Phase3_*, Phase_FINAL_* — w miarę autoryzacji)

## §5 — Reuse (budżet nowych stałych = 0)

- CG-2 (`tgp_erg_lpa_prime`): $K_{\rm IR}(\rho)=\rho$, $\nu=0.749$, $\Phi_0=2\rho_0^*=0.0609$.
- CG-34 Phase 1: silnik φ⁴ Z₂/d3 zwalidowany, $\kappa_c\approx0.189$; metoda **structure factor** (nie $\langle|\nabla\Phi|^2\rangle$).
- op-Csigma-lattice-MC: bąbel 3D $\Pi(p)=\tfrac{1}{8\pi m}-\tfrac{p^2}{96\pi m^3}$ (znak $C_\sigma>0$ DERIVED); $\kappa_E$ central 0.62.
- dodatekQ2 lematy A1–A5; PLAN_NUMERYCZNY_CG3_CG4 (N5).
- sesja #32: `rem:alpha2-pivot-status` (α=2 = postulat na gęstości; substrat amplituda → α=½).
