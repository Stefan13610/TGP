---
title: "Phase 1 — op-Csigma-lattice-MC: SETUP sieci 3D + operator O_ab + obserwabla Π(p) + WALIDACJA pipeline. WYNIK: pipeline pomiaru bąbla ZWALIDOWANY dwustopniowo (4/4 PASS) — GATE-1 (maszyneria): MC ⟨s²s²⟩_c == dokładny bąbel sieciowy (<1%); GATE-2 (continuum): bąbel sieciowy → analityczny 1/(12m²), ekstrapolacja ma→0 dev 4.8%. Flaga R-continuum potwierdzona: bąbel operatora ZŁOŻONEGO O_ab=(∂s∂s) jest UV-czuły (power-divergent) → scalar-magnitude bąbel = wzorzec (zgodnie z parent Phase 2)."
type: phase_result
status: PHASE1_COMPLETE (pipeline ZWALIDOWANY 4/4)
phase: 1
cycle: op-Csigma-lattice-MC-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14 (sesja #31): 'przeprowadzić cykl badawczy' = działaj"
script: "[[./Phase1_lattice.py]]"
output: "[[./Phase1_output.txt]]"
results_json: "[[./Phase1_results.json]]"
anti_lakatos_lock: PRESERVED
---

# Phase 1 — SETUP + WALIDACJA pipeline (op-Csigma-lattice-MC)

## §0 — Werdykt fazy w skrócie
| Test | Wynik |
|---|---|
| **GATE-1** (maszyneria: MC bąbel == dokładny bąbel sieciowy) | **PASS** — max dev 1.07%, mean 0.56% po powłokach pędu |
| **GATE-2** (continuum: bąbel sieciowy → analityczny $1/(12m^2)$) | **PASS** — ekstrapolacja $ma\to0$ daje $r\cdot m^2=0.0794$ vs $1/12=0.0833$ (dev **4.8%**) |
| **Probe R-continuum** (operator złożony UV-czuły?) | **POTWIERDZONE** — $C_{xy}(0)$ zależny od cutoffu $1/a$ → operator złożony wymaga subtrakcji |
| Pipeline $\Pi(p)$ + ekstrakcja współczynnika $p^2$ | **ZWALIDOWANY** (4/4 PASS) |

## §1 — Setup sieci i konwencje (LOCKED dla cyklu)
- Substrat: $d=3$, $n=1$, $\mathbb Z_2$ (dodatekQ; klasa uniwersalności 3D Ising). Stała sieci $a_{\rm sub}=1$ (jednostki sieci).
- Pole $\hat s$ na sieci sześciennej $L^3$, PBC. FFT konfiguracji → przestrzeń pędów; pędy $p_\mu=2\pi n_\mu/L$.
- Pęd sieciowy $\hat p^2=\sum_\mu 4\sin^2(p_\mu/2)$; pęd ciągły $p^2=\sum_\mu p_\mu^2$.
- **Notacja ŝ (flaga parent Phase 2 §1) rozstrzygnięta operacyjnie:** ŝ = pole skalarne substratu; operator
  kierunkowy $O_{ab}=(\partial_a\hat s\,\partial_b\hat s)_{\rm TF}$ jako najniższy bilinowy rank-2 w pochodnych
  skalara (różnice skończone $\partial_a\hat s=[\hat s(x+\hat e_a)-\hat s(x-\hat e_a)]/2$). Lektura wektorowa
  zmieniłaby O(1) prefaktor tensorowy, nie magnitudę (zgodne z parent).

## §2 — Operator i obserwabla
- **Operator złożony:** $O_{ab}=\partial_a\hat s\,\partial_b\hat s-\tfrac13\delta_{ab}(\partial_c\hat s)^2$ (bezśladowy, TF). Off-diagonalny $O_{xy}=\partial_x\hat s\,\partial_y\hat s$ jest automatycznie bezśladowy (czysty test).
- **Bąbel (obserwabla):** $\Pi_{ab,cd}(p)=\tfrac1V\langle\tilde O_{ab}(p)\tilde O_{cd}(-p)\rangle_{\rm c}$, gdzie $\tilde O(p)=\text{FFT}[O(x)]$, część **connected** (odjęcie $|\langle\tilde O\rangle|^2$).
- **Scalar-magnitude bąbel (wzorzec parent):** $C(p)=\tfrac1V\langle|\widetilde{\hat s^2}(p)|^2\rangle_{\rm c}$. Twierdzenie Wicka: $C(p)=2\,\Pi_{\rm scalar}(p)$ (symmetry factor 2). Stąd $C(0)=2\cdot\tfrac1{8\pi m}=\tfrac1{4\pi m}$; ratio $-A/C(0)=1/(12m^2)$ jest **niezależny od factor 2** (scheme-robust).

## §3 — GATE-1: walidacja maszynerii (MC vs dokładny bąbel sieciowy)
Wygenerowano $N_{\rm conf}=5000$ konfiguracji pola **swobodnego** (Gauss) na $L=24$, $m_0=0.30$:
propagator $\tilde\phi(p)=\sqrt{G(p)}\,\tilde\xi(p)$, $G(q)=1/(\hat q^2+m_0^2)$ (biały szum → dokładne próbki).
Porównano MC bąbel $C_{\rm MC}(p)$ z **dokładnym bąblem sieciowym** $C_{\rm lat}(p)=2\,\text{FFT}[g(r)^2]$,
$g(r)=\text{IFFT}[G(q)]$ (analityczna wartość oczekiwana, bez dwuznaczności continuum).
- Propagator odtworzony: dev $0.7\%$ ($|n|^2=1$).
- Bąbel: **max dev 1.07%, mean dev 0.56%** po powłokach $|n|^2\le12$ → maszyneria (FFT, operator, connected-subtraction) **poprawna**.

## §4 — GATE-2: zbieżność do continuum (analityczna)
Dokładny bąbel sieciowy na $L=128$, okno adaptacyjne $p^2<(0.35\,m)^2$ (zachowanie hierarchii $p^2\ll m^2\ll1/a^2$), fit z poprawką $p^4$:

| $m$ ($=ma$) | $C(0)$ | $-A/C(0)$ | $r\cdot m^2$ | dev od $1/12$ |
|---|---|---|---|---|
| 0.45 | 0.1877 | 0.3747 | 0.0759 | 8.9% |
| 0.40 | 0.2112 | 0.4776 | 0.0764 | 8.3% |
| 0.35 | 0.2411 | 0.6287 | 0.0770 | 7.6% |
| 0.30 | 0.2805 | 0.8634 | 0.0777 | 6.8% |
| 0.25 | 0.3351 | 1.2543 | 0.0784 | 5.9% |

Dewiacja **maleje monotonicznie** z $ma$. Ekstrapolacja liniowa w $(ma)^2$: $r\cdot m^2(ma\to0)=0.0794$ vs $1/12=0.0833$ → **dev 4.8%**. Pipeline reprodukuje analityczny bąbel Gaussowski parent ($\Pi(0)=1/(8\pi m)$, sztywność/masa $=1/(12m^2)$). **GATE-2 PASS.**

## §5 — Probe R-continuum (KLUCZOWE dla werdyktu): operator złożony jest UV-czuły
Analiza wymiarowa bąbla operatora ZŁOŻONEGO $O_{ab}=\partial\hat s\,\partial\hat s$ w 3D:
$$\Pi^{\partial}(p)\sim\int d^3q\;q_aq_b(q+p)_c(q+p)_d\,G(q)G(q+p)\sim\int d^3q\;\frac{q^4}{q^4}=\int d^3q\sim\Lambda^3,$$
czyli **additywna dywergencja potęgowa** (UV). Pomiar potwierdza: $C_{xy}(0)$ na $L=16,24,32$ = $\{0.0259,0.0249,0.0230\}$ — zależny od cutoffu $1/a$, **nie czysto-IR**. ⟹ Bąbel operatora złożonego $\langle O_{ab}O_{cd}\rangle$ wymaga subtrakcji/renormalizacji (mieszanie z operatorami niższego wymiaru). To **aktywuje R-continuum (HIGH)** zalockowane w Phase 0.

> **Decyzja metodologiczna (zgodna z parent Phase 2, anti-Lakatos):** magnitudę $C_\sigma$ (współczynnik $p^2$
> efektywnej sztywności) wyciągamy ze **scalar-magnitude bąbla** $\Pi(0)=1/(8\pi m)$, $[\text{coeff }p^2]=-1/(96\pi m^3)$
> — który jest UV-skończony i continuum-zbieżny (§4). Tensorowa struktura $O_{ab}$ wnosi **O(1) prefaktor tensorowy
> $t$**, nie zmienia magnitudy $C_\sigma\sim1/(96\pi m^3)$. Surowy bąbel $\langle O_{ab}O_{cd}\rangle$ pozostaje
> referencją diagnostyczną (R-continuum), nie źródłem liczby — inaczej fabrykowalibyśmy wartość z dywergencji UV.

## §6 — Handoff do Phase 2
Pipeline ZWALIDOWANY. **Phase 2 (RDZEŃ numeryczny):** zastąpić pole swobodne **realnym, nie-Gaussowskim**
substratem 3D Ising przy krytyczności (algorytm klastrowy Swendsen-Wang, R-critical-slowing), zmierzyć
scalar-magnitude bąbel $\Pi(p)$, wyciągnąć współczynnik $p^2$ → $C_\sigma$ (w jednostkach sieci), z FSS i
ekstrapolacją continuum. Walidacja krzyżowa: w reżimie słabego sprzężenia (disordered, duże $m$) MC Ising →
bąbel Gaussowski (§4). **Anti-Lakatos:** $C_\sigma$ zostanie podane z błędem stat.+syst.; zero strojenia do 5/6.
