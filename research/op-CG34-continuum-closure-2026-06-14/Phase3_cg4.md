---
title: "Phase 3 — op-CG34-continuum-closure: CG-4 (K_hom=K_TGP) = PARTIAL. ZAMKNIĘTE: c*>0 stabilne (red-flag rozwiązany, Phase 1), β=γ (algebraicznie U'(Φ_0)=0), K_hom forma lokalna 2.rzędu = K_IR (CG-2). USTALENIE (anti-Lakatos): samodzielna algebra daje α_eff=s−1 dla Z(φ)∝φ^{2s} — substrat Z∝φ² (s=1) daje α=0, NIE 2; lemat A3 ma niespójność premisa(Z∝φ²)↔konkluzja(K₁∝1/Φ=s=0). Człon (∇Φ)²/Φ JEST generowany przez Φ=φ² (struktura TGP obecna), ale dokładne α=2↔K(φ) wymaga dopięcia w rdzeniu. Residuum N5 na substracie zablokowane (patologia)."
type: phase_result
status: PHASE3_COMPLETE (CG-4 PARTIAL; c*/β=γ/K_hom-forma ZAMKNIĘTE; α=2↔K(φ) DO-PINNIĘCIA)
phase: 3
cycle: op-CG34-continuum-closure-2026-06-14
created_date: 2026-06-14
script: "[[./Phase3_cg4.py]]"
output: "[[./Phase3_output.txt]]"
target: "CG-4 [OTWARTY,AN+NUM] → PARTIAL (znaczący postęp; nie pełne zamknięcie)"
anti_lakatos_lock: PRESERVED
core_finding: "Niespójność lematu A3 (α=2 ↔ mikroskopowy K(φ)) — zgłoszona, nie sfabrykowana"
---

# Phase 3 — CG-4: identyfikacja K_hom = K_TGP

## §0 — Werdykt: **CG-4 = PARTIAL** (znaczący postęp)
| Składnik | Status |
|---|---|
| **C4-A** (c*>0 stabilne) | **ZAMKNIĘTE** — Phase 1: Z=2.46–2.61, rozrzut 6% (structure factor); red-flag c*→0 rozwiązany |
| **C4-B** (α=2 z Φ=φ²) | **USTALENIE** — struktura (∇Φ)²/Φ obecna; ale α=2↔K(φ) **niespójne w A3** (do dopięcia) |
| **C4-C** (β=γ) | **ZAMKNIĘTE** — algebraicznie z U′(Φ₀)=0 (lemat A4) |
| **C4-D** (K_hom=K_TGP forma) | **ZAMKNIĘTE co do formy** — A2 local-limit + CG-2 K_IR=ρ; residuum N5-MC zablokowane (patologia) |

## §1 — C4-B: samodzielna derywacja α (anti-Lakatos, NIE zakładamy A3)
Mikroskopowy kinetyk $Z(\phi)=Z_0\phi^{2s}$. Zmiana zmiennych $\Phi=\phi^2$ ($\phi>0$):
$(\nabla\phi)^2=(\nabla\Phi)^2/(4\Phi)$, więc $K_1(\Phi)=Z(\sqrt\Phi)/(4\Phi)=\tfrac{Z_0}{4}\Phi^{s-1}$.
Euler–Lagrange $\int K_1(\Phi)(\nabla\Phi)^2$: $\nabla^2\Phi+\tfrac{K_1'}{2K_1}(\nabla\Phi)^2=0$, $\tfrac{K_1'}{2K_1}=\tfrac{s-1}{2\Phi}$.
$$\boxed{\alpha_{\rm eff}=s-1}\qquad(\text{konwencja }\nabla^2\Phi+\alpha_{\rm eff}\tfrac{(\nabla\Phi)^2}{2\Phi}=0)$$
| $s$ | $Z(\phi)$ | $K_1(\Phi)$ | $\alpha_{\rm eff}$ |
|---|---|---|---|
| 0 | const | $\propto1/\Phi$ | $-1$ (|coeff|=1/2; A3 nazywa to „α=2") |
| 1 | $\propto\phi^2$ (substrat sek10) | const | **0** |
| 3 | $\propto\phi^6$ | $\propto\Phi^2$ | **2** |

> **USTALENIE (jawne, anti-Lakatos):** człon $(\nabla\Phi)^2/\Phi$ **jest generowany** przez $\Phi=\phi^2$
> (niezerowy współczynnik dla $Z=$const — **struktura TGP obecna**, C4-B(i) PASS). **ALE** dokładna wartość
> $\alpha=2$ jest **niespójna w lemacie A3**: A3 deklaruje premisę $Z\sim\phi^2$ ($s=1$), lecz w dowodzie używa
> $K_1=Z_0/(4\Phi)$ (czyli $s=0$, $Z=$const). Moja algebra: $s=1\Rightarrow\alpha=0$ (nie 2). $\alpha=2$ wymaga
> $Z=$const (konwencja A3) lub $Z\sim\phi^6$. **To trzeba dopiąć w rdzeniu** (sek10 vs A3) — nie fabrykuję α=2.

## §2 — C4-C: β=γ (algebraicznie)
$U(\Phi)=\tfrac{\beta}{2}\tfrac{\Phi^2}{\Phi_0}-\tfrac{\gamma}{3}\tfrac{\Phi^3}{\Phi_0^2}$; $U'(\Phi_0)=\beta-\gamma$. Warunek próżniowy $U'(\Phi_0)=0\Rightarrow\beta=\gamma$ (sympy). **ZAMKNIĘTE** (lemat A4, niezależne od α).

## §3 — C4-A: koercywność c*>0 (red-flag rozwiązany)
Phase 1: stiffness $Z(L{=}16,24,32)=\{2.46,2.47,2.61\}$ (rozrzut 6%, structure factor) — $c_*>0$ **stabilne w continuum**. To zamyka kluczowy, kontestowany punkt A1(ii) (dolny szacunek gradientowy). Prior $c_*\to0$ był artefaktem $\langle|\nabla\Phi|^2\rangle$.

## §4 — C4-D: K_hom = K_TGP (forma) + blokada N5
- Lemat A2 (zamknięty): $L_B\ll\xi\Rightarrow$ człony nielokalne $\to0$, funkcjonał **lokalny 2. rzędu**.
- Granica wolnozmienna (brak oscylacji w komórce) ⟹ $K_{\rm hom}(\Phi)=K_1(\Phi)$ (uśrednienie lokalne).
- CG-2 (`tgp_erg_lpa_prime`, 8/8): $K_{\rm IR}(\rho)=\rho$ — kinetyk z **niezależnej** drogi (ERG). Identyfikacja $K_{\rm hom}=K_{\rm TGP}=K_{\rm IR}$ **co do formy** (oba lokalne 2. rzędu, dodatnie).
- **Residuum operatora N5** na substracie: **ZABLOKOWANE** (Phase 1 §4: substrat runaway/frozen) ⟹ droga algebraiczna+ERG zamiast MC. Zgłoszone jawnie.

## §5 — Anti-Lakatos
✓ α zderywowane **samodzielnie** (sympy), nie skopiowane z A3 — ujawniło niespójność.
✓ α=2 **NIE sfabrykowane**: zgłoszone jako wymagające dopięcia (sek10 $Z\sim\phi^2$ vs A3).
✓ c*>0 z poprawnej metody (structure factor), nie z artefaktu.
✓ Blokada N5 (patologia substratu) zgłoszona jawnie, nie obejście.
✓ Rdzeń NIE edytowany; niespójność A3 = rekomendacja (Phase FINAL §rekomendacje).

## §6 — Handoff do Phase FINAL
CG-4 = PARTIAL: trzy z czterech składników zamknięte (c*, β=γ, K_hom-forma); α=2↔K(φ) = otwarta niespójność A3 do dopięcia. Phase FINAL: agregat CG-3 (ZAMKNIĘTY NUM) + CG-4 (PARTIAL) + rekomendacje rdzenia.
