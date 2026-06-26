---
title: "Phase 0 — op-Csigma-lattice-MC: LOCK falsyfikatorów F-LMC-A..D, progu κ_E=5/6 (value-blind), reguły decyzyjnej i forbidden moves PRZED jakąkolwiek liczbą. Anti-Lakatos lock aktywowany. Obiekt LOCKED: T=C_σσ_0² (jeden parametr, redundancja przeskalowania, parent Phase 3). Próg survival 5/6 EXACT (miara zero, parent), naturalna κ_E=1→7/6 FALSIFIED."
type: phase_lock
status: PHASE0_LOCKED
phase: 0
cycle: op-Csigma-lattice-MC-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14 (sesja #31): 'przeprowadzić cykl badawczy op-Csigma-lattice-MC' = działaj"
parent_cycle: "[[../op-Csigma-coarse-graining-2026-06-14/Phase_FINAL_close.md]] (PARTIAL; residual GAP = liczbowe T)"
scoping: "[[../../meta/SCOPING_op-Csigma-lattice-MC_2026-06-14.md]]"
anti_lakatos_lock: ACTIVATED
---

# Phase 0 — LOCK (op-Csigma-lattice-MC)

> **Cel cyklu:** wyznaczyć **liczbowo** tensor stiffness $T=C_\sigma\sigma_0^2$ z kierunkowego bąbla
> $\Pi_{ab,cd}(p)=\langle O_{ab}(p)O_{cd}(-p)\rangle_{\rm c}$ na sieci 3D Ising, złożyć w
> $\kappa_E=8\pi G_0T/c^3$ i wydać **value-blind** werdykt vs $5/6$. **Wszystkie progi i reguły
> zalockowane TUTAJ, przed jakąkolwiek liczbą** (anti-Lakatos).

## §1 — Obiekt LOCKED (dziedziczony z parent, niepodważalny w tym cyklu)
- Fizyczny parametr radiacyjny sektora tensorowego = **JEDEN**: $T\equiv C_\sigma\sigma_0^2$ (tensor stiffness).
  Redundancja przeskalowania $\sigma\to\lambda\sigma$, $\sigma_0\to\lambda\sigma_0$, $C_\sigma\to C_\sigma/\lambda^2$
  (parent Phase 3, sympy 3/3). **Zakaz** liczenia $C_\sigma$ i $\sigma_0$ osobno jako fizycznych.
- Złożenie: $\boxed{\kappa_E=8\pi G_0 T/c^3}$ (parent Phase 2/3, sympy).
- $C_\sigma$ = współczynnik przy $p^2$ w propagatorze kompozytu (bąblu) $\Pi(p)$
  (parent Phase 2 §1: „współczynnik przy $p^2$ w $\Pi$ daje $C_\sigma$").

## §2 — Próg LOCKED (value-blind, EXACT, miara zero)
| Wielkość | Wartość | Źródło |
|---|---|---|
| **survival** | $\kappa_E=5/6$ EXACT | parent (skalar konforemny zabiera $1/6$; tensor musiałby nieść $5/6$) |
| **naturalna** | $\kappa_E=1$ → total $7/6$ | PR-025 gałąź B (2646σ **FALSIFIED**) |
| $T_{\rm survive}$ | $\tfrac{5}{6}\cdot\tfrac{c^3}{8\pi G_0}=\tfrac{5c^3}{48\pi G_0}$ | parent Phase 3 |
| $T_{\rm natural}$ | $\tfrac{c^3}{8\pi G_0}$ (= sztywność grawitonu GR) | parent Phase 3 |
| stosunek | $T_{\rm survive}/T_{\rm natural}=5/6$ | parent Phase 3 |

**Próg $5/6$ jest LOCKED TUTAJ.** Żadna wartość $T$ z MC nie może go zmienić.

## §3 — Falsyfikatory LOCKED

| ID | Test | Reguła decyzyjna (LOCKED) |
|---|---|---|
| **F-LMC-A** | Czy $\Pi(p)$ ma mierzalny, dodatni, **continuum-zbieżny** współczynnik $p^2$ (skończone $C_\sigma>0$)? | **TAK** (zbieżny po renormalizacji operatora złożonego + ekstrapolacji $m a\to0$) → przejdź do B. **NIE** (brak skończonej sztywności / niezbieżność / dywergentne mieszanie) → **GAP strukturalny** (koniec, werdykt = GAP). |
| **F-LMC-B** | Ekstrakcja $T=C_\sigma\sigma_0^2$ w jednostkach $c^3/G_0$ (unit-bridge: $G_0\sim J\mu$, $a_\Gamma\Phi_0=1$, $\sigma_0\sim\Phi_0$). | Wartość $\pm$ błąd (stat.+syst.) → C. Jeśli systematyka unit-bridge **dominuje** i daje tylko rząd wielkości → **PARTIAL** (przedział). |
| **F-LMC-C** | $\kappa_E=8\pi G_0T/c^3$ vs $5/6$ (próg §2). | $\kappa_E=5/6\pm\sigma_{\rm tot}$ (pasmo obejmuje 5/6, **nie** obejmuje 1) → **SURVIVE**. $\kappa_E\neq5/6$ poza $\sigma_{\rm tot}$ → **FALSIFIED**. Pasmo obejmujące **i 5/6 i 1** → **PARTIAL**. |
| **F-LMC-D** | **AGREGAT** | DERIVED (A,B PASS) ∧ $\kappa_E\neq5/6$ → **FALSIFIED hard**; ∧ $\kappa_E=5/6$ → **SURVIVE**; A=GAP → **GAP**; szerokie pasmo (5/6∧1) → **PARTIAL** (sektor UNDERDETERMINED-fine-tuned, liczbowy przedział). |

**Reguła agregatu jest funkcją wyniku, NIE wyborem.** Zostanie wyliczona w Phase FINAL (sympy), nie wybrana ręcznie.

## §4 — Forbidden moves (LOCKED, anti-Lakatos)
1. **Zakaz strojenia $T$ do κ_E=5/6.** Survival musi WYNIKNĄĆ z MC, nie być dobrane. Żadnego dopasowywania prefaktora/normalizacji „aż wyjdzie 5/6".
2. **Zakaz liczenia $C_\sigma$, $\sigma_0$ osobno** jako fizycznych — tylko niezmienniczy $T=C_\sigma\sigma_0^2$ (redundancja, parent Phase 3).
3. **Zakaz rewizji parent / PR-025 / rdzenia LOCKED** (liczby 7/6, 5/6, $\det J\neq0$, $\kappa_E=8\pi G_0T/c^3$).
4. **GAP/PARTIAL jawnie** jeśli continuum nie zbiega lub unit-bridge nieokreślony — **zakaz fabrykowania wartości**.
5. **Systematyki raportowane z błędami** (continuum, FSS, renormalizacja operatora złożonego, unit-bridge) — zakaz ukrywania niepewności.
6. **Budżet nowych stałych = 0.** $T$ istnieje; $\{J,a_{\rm sub},\mu,m_0^2,\lambda_0\}$ = parametry substratu (nie nowe).
7. **Walidacja pipeline OBOWIĄZKOWA** przed pomiarem fizycznym: pipeline musi reprodukować analityczny bąbel Gaussowski ($\Pi(0)=1/(8\pi m)$, ratio $-[\text{coeff }p^2]/\Pi(0)=1/(12m^2)$) na **kontrolowanym** polu swobodnym, inaczej cały pomiar = nieważny.

## §5 — Reuse machinery (LOCKED)
- dodatekQ CG-2 `tgp_erg_lpa_prime.py` (8/8 PASS; $\rho_0^*=0.03045$, $\nu=0.749$) — wzorzec numeryczny + parametry punktu stałego WF (d=3, n=1, Z₂).
- parent Phase 2 (bąbel analityczny — wzorzec walidacji) + Phase 3 (struktura $T=C_\sigma\sigma_0^2$, unit-bridge).
- dodatekQ Q.4 ($a_\Gamma\cdot\Phi_0=1$) — most jednostkowy.

## §6 — Ryzyka aktywne od teraz (z scopingu §6)
- **R-continuum (HIGH):** bąbel operatora ZŁOŻONEGO $O_{ab}=\partial\hat s\,\partial\hat s$ ma silniejszą dywergencję UV (additywna $\sim1/a$ + multiplikatywna renormalizacja). Jeśli zrenormalizowany współczynnik $p^2$ nie zbiega przy $ma\to0$ → **F-LMC-A = GAP**.
- **R-unit-bridge (HIGH):** $T$ w $c^3/G_0$ wymaga $G_0\sim J\mu$ + $a_\Gamma\Phi_0=1$; systematyka O(1) może zdominować → **PARTIAL**.
- **R-critical-slowing (MED):** blisko krytyczności → wymóg algorytmu klastrowego (Swendsen-Wang/Wolff).
- **R-tensor-projection (MED):** definicja sieciowa $O_{ab}$ bezśladowego + notacja ŝ skalar-vs-wektor.

## §7 — Oczekiwanie (INFORMATIONAL, NIE wiąże werdyktu)
Lean parent: **DERIVED → κ_E=O(1)≠5/6 → FALSIFIED hard** (brak symetrii chroniącej 5/6). Drugi: **PARTIAL**
(unit-bridge/continuum systematyka → przedział obejmujący 5/6 i 1). Trzeci (mało prawdopodobny): **GAP**
(continuum operatora złożonego niezbieżne). Czwarty (spisek): κ_E=5/6 SURVIVE. **Werdykt rozstrzygnie MC, nie to oczekiwanie.**

## §8 — Status
**🔒 PHASE 0 LOCKED.** Anti-Lakatos lock AKTYWNY. Próg 5/6, falsyfikatory F-LMC-A..D, reguła agregatu,
forbidden moves — wszystko zalockowane przed pierwszą liczbą. Przejście do Phase 1 (setup + walidacja pipeline).
