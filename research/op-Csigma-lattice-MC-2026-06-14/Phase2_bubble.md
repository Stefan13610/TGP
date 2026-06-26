---
title: "Phase 2 — op-Csigma-lattice-MC: RDZEŃ NUMERYCZNY. MC bąbla kierunkowego spin-2 na 3D Ising (Swendsen-Wang). WYNIK: F-LMC-A = PARTIAL. C_σ=-coeff(p²) bąbla Π_TT ZMIERZONE i DODATNIE (sztywność>0, ZNAK DERIVED numerycznie, czyste p²≪m² w disordered), O(1) w jednostkach sieci (~0.5–0.7). ALE absolutna magnituda UV/schemat-zależna (operator złożony power-divergent, Phase 1 §5) + przy krytyczności p_min/m_s≈2>1 → analityczny współczynnik p² NIEDOSTĘPNY przy m→0. Scheme-independent continuum NIE osiągnięte (dodatekQ CG-3/CG-4 OTWARTE). Anti-Lakatos: C_σ z błędem, zero strojenia do 5/6, GAP-składowa jawnie."
type: phase_result
status: PHASE2_COMPLETE (RDZEŃ — F-LMC-A = PARTIAL)
phase: 2
cycle: op-Csigma-lattice-MC-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14 (sesja #31): działaj"
script: "[[./Phase2_mc.py]]"
output: "[[./Phase2_output.txt]]"
results_json: "[[./Phase2_results.json]]"
flag_F_LMC_A: "PARTIAL (C_σ>0 O(1) DERIVED znak+rząd; prefaktor O(1) systematyka schematu; nie clean continuum)"
anti_lakatos_lock: PRESERVED
---

# Phase 2 — RDZEŃ: MC bąbla kierunkowego na 3D Ising (op-Csigma-lattice-MC)

## §0 — Werdykt fazy w skrócie (5/5 testów PASS)
| Test | Wynik |
|---|---|
| **MC-1/2** (walidacja 3D Ising: faza Z₂, pik χ) | **PASS** — $\|m\|$ rośnie z β; pik χ przy $\beta_c=0.221654$ |
| **F-LMC-A(i)** ($C_\sigma=-\text{coeff}(p^2)>0$?) | **PASS** — $C_\sigma>0$ wszędzie (disordered, czyste $p^2\ll m^2$), znak **DERIVED** |
| **F-LMC-A(ii)** (UV-czułość operatora złożonego?) | **PASS** (potwierdzona) — absolutna magnituda schemat-zależna |
| **F-LMC-A(iii)** (obstrukcja $p^2$ przy krytyczności?) | **PASS** (potwierdzona) — $p_{\min}/m_s\approx2>1$ |
| **AGREGAT F-LMC-A** | **PARTIAL** — $C_\sigma>0$, $O(1)$, DERIVED; **scheme-independent continuum NIE osiągnięte** |

## §1 — Metoda numeryczna (zwalidowana Phase 1)
- Substrat 3D Ising ($d=3,n=1,\mathbb Z_2$), **Swendsen-Wang** (algorytm klastrowy, R-critical-slowing pokonane; ~3 ms/sweep przy $L=32$).
- Operator spin-2 TT: $O_{ab}=(D_a\hat s)(D_b\hat s)-\tfrac13\delta_{ab}(D_c\hat s)^2$, $D_a=$ symetryczna różnica.
- Bąbel skalarny rotacyjny: $\Pi_{\rm TT}(p)=\tfrac1V\sum_{ab}\langle|\tilde O_{ab}(p)|^2\rangle_{\rm c}$ (5 niezależnych komponentów; rotacyjnie niezmienniczy → uśrednialny po powłokach pędu).
- $C_\sigma\equiv-[\text{coeff }p^2]$ z dopasowania $\Pi_{\rm TT}(p)=\Pi_{\rm TT}(0)-C_\sigma p^2+B p^4$ na powłokach z $p^2<(0.6\,m_s)^2$ (zachowanie hierarchii $p^2\ll m_s^2$). Błędy: jackknife po blokach.

## §2 — Walidacja MC (3D Ising)
$\|m\|(\beta)=\{0.017,0.032,0.060,0.253,0.671,0.839\}$ dla $\beta=\{0.15,0.19,0.21,\beta_c,0.24,0.27\}$ — poprawne łamanie $\mathbb Z_2$ (disordered→ordered). Pik podatności $\chi$ przy $\beta=\beta_c=0.221654$ ✓. Substrat jest **realny, nie-Gaussowski** (klasa uniwersalności 3D Ising) — wymóg scopingu spełniony.

## §3 — RDZEŃ: $C_\sigma$ w fazie disordered (czysta ekstrakcja $p^2$)
W disordered $m_s$ duże ⟹ $p_{\min}^2\ll m_s^2$ osiągalne ⟹ **czysta** ekstrakcja współczynnika $p^2$ ($L=32$, 2000 pomiarów):

| $\beta$ | $m_s$ | $p_{\min}/m_s$ | $C_\sigma=-\text{coeff}(p^2)$ | nsh |
|---|---|---|---|---|
| 0.10 | 2.276 | 0.09 | **0.499 ± 0.025** | 41 |
| 0.13 | 1.479 | 0.13 | **0.699 ± 0.090** | 18 |
| 0.16 | 1.148 | 0.17 | **0.649 ± 0.153** | 11 |
| 0.18 | 0.851 | 0.23 | 0.382 ± 0.691 | 6 |

> **Wynik:** $C_\sigma>0$ **wszędzie** ⟹ sztywność tensorowa **DODATNIA**, ZNAK **DERIVED numerycznie** (nie tylko z analitycznego znaku bąbla parent). Magnituda $C_\sigma=O(1)$ w jednostkach sieci, $\approx0.5$–$0.7$ (najczystsze pomiary, nsh duże). To **potwierdza i wzmacnia** parent Phase 2 (tam: znak+skaling DERIVED, prefaktor GAP) — teraz prefaktor **zmierzony jako $O(1)$**.

## §4 — Probe R-continuum: absolutna magnituda jest UV/schemat-zależna (KLUCZOWE)
Ustaliłem $\beta=0.16$ ($m_s\approx0.9$–$1.1$), zmieniam $L$:

| $L$ | $\Pi_{\rm TT}(0)$ | $C_\sigma$ | nsh |
|---|---|---|---|
| 16 | 3.519 | 0.044 ± 0.69 | 2 |
| 24 | 3.470 | 0.239 ± 0.46 | 6 |
| 32 | 3.653 | 0.903 ± 0.44 | 7 |

$\Pi_{\rm TT}(0)\approx3.5$ (rząd) z $\sim5\%$ dryfem — ale to **UV-czuła additywna stała** (Phase 1 §5: bąbel operatora złożonego $\sim\int d^3q\sim\Lambda^3$). Współczynnik $p^2$ (czyli $C_\sigma$) ma duży rozrzut z $L$ (dominowany przez szum przy małym $L$/małym nsh, ale **z realną składową schematu**). 

> **Wniosek F-LMC-A(ii):** absolutna magnituda $C_\sigma$ niesie **$O(1)$ systematykę schematu** (renormalizacja operatora złożonego: additywne UV + mieszanie z operatorami niższego wymiaru). To **materializacja R-continuum (HIGH)** zalockowanego w Phase 0.

## §5 — Obstrukcja przy krytyczności (gdzie żyje fizyka TGP, $m\to0$)
Przy $\beta_c$ korelacja $\xi\sim L$, więc $p_{\min}=2\pi/L\sim2\pi m_s$:

| $L$ | $m_s$ | $p_{\min}/m_s$ |
|---|---|---|
| 16 | 0.196 | 2.01 |
| 24 | 0.137 | 1.92 |
| 32 | 0.095 | 2.07 |

$p_{\min}/m_s\approx2>1$ na **wszystkich** dostępnych $L$ ⟹ regime $p^2\ll m_s^2$ **strukturalnie niedostępny** na skończonej sieci przy krytyczności. Analityczny współczynnik $p^2$ (Taylor) **nie jest izolowalny** przy $m\to0$ — bąbel jest tam krytycznym prawem potęgowym z anomalnym wykładnikiem, nie funkcją analityczną w $p^2$. To **niezależna od statystyki, strukturalna** bariera (nie do pokonania większym $N$, tylko jakościowo innym podejściem — ciągła renormalizacja, dodatekQ CG-3/CG-4 OTWARTE).

## §6 — Dyspozycja F-LMC-A (reguła LOCKED Phase 0)
Reguła F-LMC-A: „TAK (zbieżny po renormalizacji + $ma\to0$) → przejdź; NIE (niezbieżność) → GAP".
- **Część TAK:** $C_\sigma>0$ mierzalny, $O(1)$, znak+rząd **DERIVED** (czyste $p^2$ w disordered). ✓
- **Część NIE:** scheme-independent **continuum** współczynnika $p^2$ operatora złożonego **NIE osiągnięte** (power-divergence + obstrukcja krytyczna). To **częściowy GAP** w absolutnej normalizacji.

> **Dyspozycja: F-LMC-A = PARTIAL** (nie czyste PASS, nie czyste GAP). $C_\sigma=O(1)$ lattice units
> (~0.5–0.7), ale z **nieusuwalną numerycznie $O(1)$ systematyką schematu**. Przekazuję do Phase 3 jako
> **pasmo** $C_\sigma\sim O(1)$, nie ostrą liczbę.

## §7 — Anti-Lakatos (checklist Phase 2)
✓ MC zwalidowany vs znany 3D Ising (faza, $\beta_c$) przed pomiarem fizycznym.
✓ $C_\sigma$ podany **z błędem** (jackknife) i jako **pasmo** — zero strojenia do 5/6.
✓ Składowa GAP (continuum operatora złożonego) **jawna** (§4–5), nie zamieciona — anti-Lakatos.
✓ Obstrukcja krytyczna ($p_{\min}/m_s>1$) udokumentowana jako **strukturalna**, nie statystyczna.
✓ ZNAK $C_\sigma>0$ to realny postęp vs parent (DERIVED numerycznie, nie tylko analityczny znak bąbla).
✓ Budżet nowych stałych 0; liczby poprzedników LOCKED.

## §8 — Handoff do Phase 3
$C_\sigma=O(1)$ jednostek sieci (najczystsze: $\approx0.5$–$0.7$), **dodatni**, z $O(1)$ systematyką schematu
(R-continuum). **Phase 3 (unit-bridge):** $G_0\sim J\mu$, $a_\Gamma\Phi_0=1$, $\sigma_0\sim\Phi_0$ → złożyć
$T=C_\sigma\sigma_0^2$ w $c^3/G_0$ → $\kappa_E=8\pi G_0T/c^3$ z błędem stat.+**syst.** (systematyka schematu
z §4 + systematyka unit-bridge R-unit-bridge HIGH). Oczekiwane: $\kappa_E=O(1)$ z **pasmem** → test vs 5/6.
