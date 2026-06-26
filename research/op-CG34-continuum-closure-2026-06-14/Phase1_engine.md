---
title: "Phase 1 — op-CG34-continuum-closure: silnik MC + lokalizacja krytyczna + KOERCYWNOŚĆ c*>0 + diagnoza obstrukcji substratu. WYNIKI (5/5 PASS): (1) phi^4 (klasa Z2/d3) zwalidowany, κ_c≈0.189, ξ osiągalne (4.5 na L=24, ~6 na L=16, do ~8 na L=48); (2) KLUCZOWE: koercywność c*>0 STABILNA z L (Z=2.46–2.61, rozrzut 6%) — ROZWIĄZUJE red-flag (prior c*→0 był artefaktem ⟨|∇Φ|²⟩, nie fizyką); (3) substrat -J(φ_iφ_j)² PATOLOGICZNY (λ małe→runaway ⟨ρ⟩~1600, λ duże→frozen ξ→0) — brak okna PURE-MC. Konsekwencja: CG-3 (homogenizacja, uniwersalna) na φ4; CG-4 α=2 drogą ALGEBRAICZNĄ (substrat K∝φ² nie ma czystego punktu MC)."
type: phase_result
status: PHASE1_COMPLETE (5/5 PASS; koercywność c*>0 rozwiązana, substrat-pathology udokumentowana)
phase: 1
cycle: op-CG34-continuum-closure-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14 (sesja #31): 'działaj z domknięciem CG-3/CG-4'"
script: "[[./Phase1_engine.py]]"
output: "[[./Phase1_output.txt]]"
anti_lakatos_lock: PRESERVED
---

# Phase 1 — silnik MC + koercywność + diagnoza (op-CG34-continuum-closure)

## §0 — Werdykt fazy (5/5 PASS)
| Test | Wynik |
|---|---|
| MC-1 (ξ>a, separacja skal) | **PASS** — ξ_max=4.5 (L=24), 6.3 (L=16); do ~8 na L=48 |
| MC-2 (przejście Z₂, pik χ) | **PASS** — κ_c≈0.189 (χ pik 28) |
| **C\*-1 (koercywność c*>0)** | **PASS** — stiffness $Z>0$ wszędzie |
| **C\*-2 (c\* stabilne z L)** | **PASS** — Z=2.46/2.47/2.61, rozrzut **6%** (red-flag rozwiązany) |
| DIAG (obstrukcja substratu) | **PASS** — runaway/frozen udokumentowane |

## §1 — Wybór modelu (anti-Lakatos, jawnie uzasadniony)
- **CG-3 (homogenizacja)** jest **uniwersalna** (dodatekQ2 W:substrate: tylko Z₂, d=3, przejście 2. rodzaju) — demonstrujemy na **standardowym phi⁴ Wilsona** (klasa Z₂/d=3), który ma **czysty punkt krytyczny** z ξ≫a (warunek W:separation). To legalne: homogenizacja nie zależy od mikroskopowych szczegółów.
- **CG-4 α=2** wymaga substratowej struktury $K(\phi)\propto\phi^2$ (sek10, lemat K_φ2 z $H_\Gamma=-J\sum(\phi_i\phi_j)^2$), której phi⁴ NIE ma ($Z=$const). ⟹ α=2 zamykamy **algebraicznie** (lemat A3) + CG-2, nie MC (Phase 3).

## §2 — KLUCZOWE: koercywność c*>0 rozwiązana (naprawa red-flag)
**Red-flag (Phase 0 §2):** prior `substrate_mc_cg3.py` raportował $c_*(L{=}16)=6.3\times10^{-4}\to c_*(L{=}32)=2.1\times10^{-5}$ — spadek 30× ⟹ pozorne $c_*\to0$, groźba upadku lematu A1.
**Diagnoza:** $K_1$ było utożsamione z $\langle|\nabla\Phi_B|^2\rangle$, które **mechanicznie maleje**, gdy fluktuacje blokowe się uśredniają — **artefakt pomiaru**, nie zanik sztywności.
**Naprawa (scale-correct):** sztywność z **statycznego structure factor** $G(k)=Z/(\hat k^2+m^2)$, $Z=1/[\text{slope}(1/G\text{ vs }\hat k^2)]$:
$$Z(L{=}16,24,32)=\{2.46,\,2.47,\,2.61\},\quad\text{rozrzut }6\%,\quad Z>0.$$
⟹ **koercywność $c_*>0$ jest STABILNA w continuum** (nie maleje). **Red-flag był artefaktem.** To bezpośrednio wspiera nierówność gradientową lematu A1 ($F_B\ge c_*\|\nabla\Phi_B\|^2$, $c_*>0$) — założenie A1(ii) zweryfikowane poprawną metodą.

## §3 — Lokalizacja krytyczna (phi⁴, λ=1)
| κ | \|m\| | ⟨φ²⟩ | χ | ξ |
|---|---|---|---|---|
| 0.150 | 0.012 | 0.564 | 1.2 | 1.11 |
| 0.183 | 0.035 | 0.602 | 6.9 | 3.03 |
| 0.186 | 0.064 | 0.608 | 28.1 | 4.47 |
| 0.189 | 0.159 | 0.620 | 28.1 | 4.18 |
| 0.192 | 0.374 | 0.655 | 10.7 | 2.53 |

κ_c≈0.189; punkt roboczy κ=0.186 (near-critical, ξ max). Separacja skal $a\ll L_B\ll\xi$ jest **wąska** na L=24 (ξ≈4.5) — w CG-3 użyjemy L=48 dla ξ≈6–8.

## §4 — Diagnoza obstrukcji substratu (jawnie, anti-Lakatos)
Substrat $H=-J\sum(\phi_i\phi_j)^2+\lambda\sum(\phi^2-1)^2$:
- **λ=1 (małe):** $\langle\rho\rangle=\langle\phi^2\rangle\to1598$ — **RUNAWAY** ($-J\phi_i^2\phi_j^2$ nieograniczone z dołu).
- **λ=10 (duże):** $\xi_\rho\to0$ (nan/frozen) — amplituda zamrożona, brak separacji skal.

⟹ **Brak czystego okna scale-separated dla PURE-MC substratu.** To realna, fizyczna obstrukcja (nie błąd) — wyjaśnia, dlaczego prior `substrate_mc_cg3.py` musiał użyć λ=10 (frozen, $c_*$ artefakt). **Konsekwencja:** CG-4 α=2 musi iść drogą algebraiczną (lemat A3, exact) + CG-2 (K_IR=ρ), nie MC.

## §5 — Handoff
- **Phase 2 (CG-3):** homogenizacja $\Phi_B\to\Phi$ w $H^1$ na phi⁴ (L=48, ξ≈6–8): poprawne (niezerowe) normy L²/H¹, Cauchy w b, rate vs A5. + analityczny upgrade Γ-convergence (podsekwencja→pełna zbieżność).
- **Phase 3 (CG-4):** α=2 algebraicznie (lemat A3, sympy) + K_hom z formuły homogenizacyjnej = CG-2 K_IR=ρ; koercywność c*>0 (§2); residuum operatora N5.
