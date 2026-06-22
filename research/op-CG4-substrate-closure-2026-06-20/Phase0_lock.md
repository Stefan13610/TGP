---
title: "Phase 0 — op-CG4-substrate-closure: LOCK kryteriów. Cel: niepatologiczny model substratu (stabilny + czysty punkt krytyczny Z₂ + emergentne (∇Φ)²/Φ) umożliwiający scheme-independent C_σ < faktor 1.2 → zwężenie pasma κ_E → twardy werdykt sektora radiacyjnego GW. Wszystkie progi PASS/FAIL, falsyfikatory, reguła agregatu i forbidden moves zalockowane TUTAJ, przed pierwszą liczbą."
type: phase_lock
status: PHASE0_LOCKED
phase: 0
cycle: op-CG4-substrate-closure-2026-06-20
created_date: 2026-06-20
authorization: "User 2026-06-20: 'działaj z op-CG4-substrate-closure'"
target: "CG-4 [CZĘŚCIOWY NUM] → [ZAMKNIĘTY NUM] / GAP / FALSIFIED-hard; N5 residuum operatora"
anti_lakatos_lock: ACTIVATED
---

# Phase 0 — LOCK (op-CG4-substrate-closure)

> **Reguła nadrzędna (anti-Lakatos):** wszystkie kryteria PASS/FAIL, progi dyskryminacji, falsyfikatory,
> reguła agregatu i forbidden moves są zamknięte **w tym pliku, przed pierwszą liczbą Phase 1**. Werdykt
> Phase FINAL będzie **wyliczony** z reguły §4, nie wybrany.

## §1 — Stan wejściowy (co dokładnie jest otwarte)

Z [[../op-CG34-continuum-closure-2026-06-14/Phase_FINAL_close.md]] §4 i
[[../op-Csigma-lattice-MC-2026-06-14/Phase_FINAL_close.md]] §4:

| Element | Stan przed cyklem |
|---|---|
| $C_\sigma$ (sztywność tensora σ_ab) | znak DERIVED, $O(1)$ ZMIERZONE (~0.5–0.7), **prefaktor scheme-dependent** |
| $\kappa_E=8\pi G_0 C_\sigma\sigma_0^2/c^3$ | central ≈0.62, pasmo **[0.04, 11.1]** (obejmuje 5/6 ∧ 1) |
| substrat $-J(\phi_i\phi_j)^2$ | **PATOLOGICZNY** (runaway/frozen) — brak czystego punktu krytycznego |
| operator złożony $O_{ab}=\partial\hat s\,\partial\hat s$ | **UV power-divergent** (R-continuum) — wymaga subtrakcji + mieszania |
| α=2 ↔ $K(\phi)$ | $\alpha_{\rm eff}=s-1$; #32: α=2 = **postulat na gęstości** (nie derywacja); A3 do ujednolicenia |

## §2 — Definicje operacyjne (LOCKED)

- **Substrat niepatologiczny** = model spełniający jednocześnie C-A (stabilność) ∧ C-B (czysty punkt krytyczny).
- **Próg przeżycia** (z parent #30/#31, niezmienny): $\kappa_E=5/6$ **dokładnie** ⟹ SURVIVE; $\kappa_E\neq5/6$ ⟹ FALSIFIED-hard.
- **Dokładność dyskryminująca** (LOCKED): $C_\sigma$ scheme-independent z **względnym pasmem $<$ faktor 1.2**
  (tj. $C_\sigma^{\max}/C_\sigma^{\min}<1.2$), co przekłada się na pasmo $\kappa_E$ węższe niż odstęp 1 ↔ 5/6 (=1.2×).
- **Φ = gęstość** $\langle\hat s^2\rangle$ (pole kanoniczne, sek01 `def:Phi`); $\phi$ = amplituda; $\Phi=\phi^2$.

## §3 — Kryteria domknięcia LOCKED (PASS / PARTIAL / FAIL)

### C-A — Stabilność termodynamiczna
| Reguła PASS | Reguła FAIL |
|---|---|
| Potencjał efektywny **ograniczony z dołu**; $\langle\rho\rangle=\langle\hat s^2\rangle$ skończone, stabilne w rozmiarze $L$ (rozrzut <20%); brak runaway ($\langle\rho\rangle$ nie rośnie z $L$) i brak frozen ($\xi>2a$ w fazie krytycznej) | $\langle\rho\rangle$ rośnie z $L$ (runaway) **lub** $\xi\to0$ przy każdym dostępnym $T$ (frozen) |

### C-B — Czysty ciągły punkt krytyczny Z₂
| Reguła PASS | Reguła FAIL |
|---|---|
| Istnieje $T_c$ z **ciągłym** przejściem Z₂: pik podatności $\chi$ skalujący się z $L$; $\xi(T\!\to\!T_c)$ rosnące z $L$ (kontrolowane $\xi\to\infty$); brak histerezy/skoku $\langle\rho\rangle$ (nie 1. rzędu) | Przejście 1. rzędu (skok/histereza) **lub** brak rozbieżności $\xi$ |

### C-C — Emergentna struktura kinetyczna (spójność, NIE derywacja α)
| Reguła PASS | Reguła FAIL / USTALENIE |
|---|---|
| Człon $(\nabla\Phi)^2/\Phi$ **generowany** (niezerowy współczynnik przy $\Phi=\phi^2$); zmierzony wykładnik $s$ modelu zgodny z $\alpha_{\rm eff}=s-1$; **niesprzeczność** z `rem:alpha2-pivot-status` (α=2 = postulat na gęstości) | struktura $(\nabla\Phi)^2/\Phi$ nieobecna (sprzeczność z TGP) ⟹ FAIL; jeśli model wymusza $\alpha_{\rm eff}\neq$ wartości spójnej z #32 ⟹ USTALENIE jawne (nie fabrykować α=2) |

### C-D — Scheme-independent $C_\sigma$ → dyskryminacja $\kappa_E$ (BRAMA TWARDA)
| Reguła PASS | Reguła PARTIAL | Reguła FAIL/GAP |
|---|---|---|
| $C_\sigma$ po subtrakcji UV + mieszaniu **continuum-zbieżne**, pasmo $<$ faktor 1.2; pasmo $\kappa_E$ **wyklucza 5/6 albo wyklucza wszystko poza 5/6** | $C_\sigma$ zbieżne, ale pasmo $\kappa_E$ wciąż obejmuje 5/6 ∧ 1 (jak #31) | continuum nie istnieje (operator nie renormalizuje się scheme-independent) ⟹ GAP twardy |

## §4 — Reguła agregatu LOCKED (value-blind)

```
WARUNEK BRAMY: substrat niepatologiczny = C-A PASS ∧ C-B PASS.
  Jeśli ŻADEN kandydat nie spełnia C-A ∧ C-B → CG-4 = GAP (substrat-blocked; jawnie).

Jeśli BRAMA otwarta (∃ kandydat C-A ∧ C-B PASS):
  ∧ C-D PASS ∧ κ_E-band wyklucza 5/6        → CG-4 ZAMKNIĘTY ⟹ sektor radiacyjny = FALSIFIED-hard
  ∧ C-D PASS ∧ κ_E-band ⊂ {5/6 ± faktor1.2} → CG-4 ZAMKNIĘTY ⟹ sektor radiacyjny = SURVIVE (fine-tuned, ale chronione)
  ∧ C-D PARTIAL (band obejmuje 5/6 ∧ 1)      → CG-4 PARTIAL ⟹ sektor UNDERDETERMINED-fine-tuned (bez zmiany vs #31)
  ∧ C-D FAIL/GAP                              → CG-4 GAP ⟹ residual operatorowy nieusuwalny tą drogą
C-C jest wspierające (spójność), nie wchodzi do bramy werdyktu sektora.
```

- Werdykt **wyliczony** z tej reguły w Phase FINAL (sympy/skrypt agregatu), nie wybrany ręcznie.
- **Lean strukturalny (jawny, z #30/#31): FALSIFIED-hard** — $\kappa_E=5/6$ to miara zero; wartość naturalna $\sim1$.
  Reguła dopuszcza SURVIVE **tylko** gdy pomiar realnie wskaże 5/6 ± faktor 1.2. Zero strojenia do 5/6.

## §5 — Falsyfikatory (konkretne, dwustronne)

1. **FALSIFIED-hard sektora** ⟺ (∃ niepatologiczny substrat) ∧ (scheme-indep. $C_\sigma$, pasmo <1.2×) ∧ ($\kappa_E$ wyklucza 5/6). To jest **dopuszczalny i oczekiwany** wynik (lean).
2. **SURVIVE** ⟺ to samo, ale $\kappa_E$ pokrywa 5/6 i wyklucza 1. Wymaga **pozytywnej** zgodności, nie braku falsyfikacji.
3. **GAP** ⟺ żaden niepatologiczny substrat (C-A∧C-B FAIL dla wszystkich kandydatów) **lub** operator nie renormalizuje się scheme-independent (C-D GAP).
4. **PARTIAL** ⟺ continuum zbieżne, ale za szerokie (jak #31) — brak postępu werdyktowego, ale udokumentowany model.

## §6 — Forbidden moves (anti-Lakatos)

1. **Zakaz strojenia** parametrów substratu/schematu, by $\kappa_E$ „wyszło" 5/6. $\kappa_E$ mierzone, raportowane z jawnym pasmem.
2. **Zakaz fabrykowania α=2** z substratu (#32: α=2 = postulat na gęstości). C-C testuje tylko obecność struktury, nie wartość.
3. **Zakaz utożsamiania $C_\sigma$ z surową (nieodjętą) magnitudą operatora złożonego** — obowiązkowa subtrakcja additywna UV + analiza mieszania (R-continuum). Surowa wartość ≠ scheme-independent.
4. **Walidacja silnika MC OBOWIĄZKOWA** przed pomiarem fizycznym (faza Z₂, $\langle\phi^2\rangle$, $\xi$, pik $\chi$) — inaczej pomiary nieważne (forbidden move CG-34 #4).
5. **GAP/PARTIAL jawnie** — zakaz naciągania do „ZAMKNIĘTY", jeśli pasmo nie zwęża się <1.2× lub brak niepatologicznego substratu.
6. **Zakaz rewizji rdzenia** w tym cyklu — status update = rekomendacja, nie edycja core .tex (forbidden CG-34 #5). Budżet nowych stałych 0.
7. **Zakaz multi-candidate selection bez argumentu strukturalnego** (wzorzec κ.1/θ.1 z CALIBRATION_GATE) — wybór modelu substratu musi mieć uzasadnienie fizyczne (uniwersalność/stabilność), nie minimalizację dryfu.

## §7 — Kandydaci substratu do skanu Phase 1 (LOCKED lista wejściowa)

| ID | Model $H_\Gamma$ | Hipoteza wejściowa |
|---|---|---|
| **M0** | φ⁴ Z₂ Ginzburg–Landau: $\sum(\nabla s)^2 + r\,s^2 + u\,s^4$, $u>0$ | KNOWN-GOOD (klasa Isinga 3D); silnik CG-34 Phase 1; baseline stabilny + WF punkt krytyczny |
| **M1** | bond $-J\sum(s_i s_j)^2$ (obecny) | PATOLOGICZNY (runaway/frozen) — kontrola negatywna |
| **M2** | gradient-bond v2 (README): $+J\sum A_{ij}s_i^2 s_j^2(s_j^2-s_i^2)^2$ + on-site $V(s)$ | ≥0 zawsze, ale kierunek płaski (uniform) → ryzyko frozen; wymaga on-site bound |
| **M3** | φ⁶ tricritical: $r\,s^2+u\,s^4+w\,s^6$, $w>0$ | bounded-below nawet $u<0$; strojenie do trikrytyczności (czysty punkt) |

Skan Phase 1 = analityczna klasyfikacja (stabilność bounded-below + typ przejścia mean-field/Landau) → wybór
kandydata(ów) do MC Phase 2 **z argumentem strukturalnym** (forbidden #7).

## §8 — Reuse i budżet

- Reuse: CG-2 ($K_{\rm IR}=\rho$, $\nu=0.749$, $\Phi_0=0.0609$), CG-34 silnik+structure-factor, bąbel 3D op-Csigma-lattice-MC, A1–A5, `rem:alpha2-pivot-status` (#32). Budżet nowych stałych = **0**.
- Phase0_balance.md (obowiązkowy gate) — patrz [[Phase0_balance.md]].

## §9 — Status

**🔒 PHASE 0 LOCKED (2026-06-20).** Anti-Lakatos aktywny. Kryteria C-A..C-D, reguła agregatu §4, falsyfikatory
§5, forbidden moves §6, lista kandydatów §7 — zalockowane przed pierwszą liczbą. Lean jawny: FALSIFIED-hard.
Przejście do Phase 1 (silnik analityczny: klasyfikacja kandydatów + $\alpha_{\rm eff}=s-1$ + znak $C_\sigma$) — po „działaj Phase 1".
