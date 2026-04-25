# OP-EHT — strong-field photon ring TGP M9.1'' vs EHT

**Data otwarcia:** 2026-04-25
**Status:** **OPEN** — strukturalna ramka, sub-tests T1-T5 do uruchomienia
**Następnik:** OP-7 (CLOSED 2026-04-25, 94/97 = 96.9% PASS) — niezależny.
**Punkt wyjścia:** [[research/op7/EHT_quick_results.md]] (5/6 PASS,
INCONCLUSIVE-leaning-NEGATIVE: M9.1'' photon ring **+14.6% over GR uniformly**,
M87* marginal, Sgr A* outside 4.4% EHT 2022 envelope).
**Cross-references:** [[research/op-newton-momentum/M9_1_pp_P3_results.md]]
(strong-field test "open"); [[research/op7/OP7_T1_results.md]]
(σ_ab = 0 for static spherical — sigma sector cannot rescue);
[[paper/tgp_core.tex]] §applications BH shadows + falsifiability F4.

---

## 1. Cel OP-EHT

Domknąć (lub falsyfikować) **strong-field consistency** TGP M9.1''
hyperbolic metric w reżimie U_N ~ 0.2-0.3 (photon ring), gdzie:
- post-Newtonian expansion z PPN c_n nie zbiega się elementarnie,
- analitical 7-term tail daje robust +14.6% deviation w b_crit
  (impact parameter) i shadow diameter,
- M87* ledwie mieści się w EHT 2019 ~10% systematic,
- Sgr A* tension przy 19.7% > EHT 2022 4.4% bound.

**Pytanie centralne:**

> Czy deviację +14.6% w b_crit można wyjaśnić strukturalnie w obrębie
> M9.1'' (truncation, calibration, environment) lub czy jest to genuine
> TGP prediction wymagająca M9.2 pivot lub deferral do ngEHT 2030+?

Jeśli +14.6% jest **artefaktem** (truncation/calibration) — TGP zostaje
unfalsified, OP-EHT zamyka POSITIVE.
Jeśli +14.6% jest **genuine + akceptowalne** w EHT systematyce — OP-EHT
zamyka CONDITIONAL POSITIVE z ngEHT falsification window.
Jeśli +14.6% jest **genuine + niefalsyfikowane przez EHT 2024+** — OP-EHT
zamyka CONDITIONAL NEGATIVE; M9.1'' jako nieadekwatne dla strong-field,
wymaga M9.2 pivot.

## 2. Tło — co już wiemy

### 2.1 EHT-quick wynik (2026-04-25)

```
Test                              | M9.1''      | GR Schwarzschild | TGP/GR
----------------------------------|-------------|------------------|-------
r_ph (areal)                      | 3.880 r_g   | 3.000 r_g        | 1.293
b_crit                            | 5.952 r_g   | 5.196 r_g        | 1.146
Shadow diameter (uniform)         | +14.6%      | baseline         | 1.146
M87* shadow (45.48 / 39.70 µas)   | +8.3% vs EHT (42.0±3.0)  | within 10%
Sgr A* shadow (62.01 / 54.14 µas) | +19.7% vs EHT (51.8±2.3) | OUTSIDE 4.4%
```

Konwergencja numeryczna: spread < 1e-6 across 5 runs (r_out ∈ {500,
1000, 2000, 5000} r_g, n_terms ∈ {5, 6, 7}). **+14.6% jest robust**
przy 7-term PN tail; pytanie: czy 8-, 10-, 15-term tail ten samosć?

### 2.2 OP-7 T1 zamyka σ_ab path

Dla static spherical metric M9.1'':
- σ_ab = K_ab − (1/3)δ_ab Tr(K) z K_ab = ⟨(∂_a ŝ)(∂_b ŝ)⟩_B
- W próżni izotropowej (∂_a ŝ izotropowe) ⇒ σ_ab = 0
- ⇒ żaden tensor sector contribution w static BH

Konkluzja OP-7 T1 + EHT-quick: **żaden σ_ab dynamics nie może uratować
+14.6% deviation w izolowanym static spherical BH**. Pivot musi
przyjść z innego strukturalnego źródła.

### 2.3 Trzy strukturalne paths identifikowane w EHT_quick §6

1. **Pivot w strong-field** — np. matter-coupling correction
   (q-renormalization), non-vacuum BH environment, lub M9.2
   momentum back-reaction.
2. **Genuine TGP deviation undetectable obecnie** — mass uncertainties
   M87* ±10%, Sgr A* ±5%, ring modeling +/- ~5%. Combined envelope
   może rozszerzyć się do ~15-20%, gdzie TGP staje się marginally
   consistent. Falsifikacja z ngEHT 2030+ przy 1% precision.
3. **Recalibration `Phi_0`** — matching A = GM/(2c²) jest 1PN.
   Higher-PN matter coupling (G(Φ) = G_0/ψ w strong field) może
   przesunąć A i tym samym r_ph + b_crit.

Wszystkie trzy są **strukturalnie różne** od T1-σ_ab path; wymagają
oddzielnych sub-testów.

### 2.4 Co OP-EHT musi dostarczyć

1. **PN truncation robustness audit** (T1): czy +14.6% jest robust
   przy n=10, 12, 15 term tail? Czy dalsze c_n wnoszą >1% poprawkę
   w b_crit?
2. **Mass calibration MC** (T2): propagacja niepewności pomiaru
   masy BH, distance, ring modeling przez +14.6% prediction.
   Czy effective EHT envelope rozszerza się na tyle że TGP mieści się?
3. **Q-renormalization audit** (T3): czy G(Φ) = G_0/ψ matter-coupling
   w strong field przesuwa A = GM/(2c²) na tyle że b_crit obniża się
   do +5% lub mniej?
4. **Non-vacuum environment** (T4): accretion disk + plasma ray-tracing
   correction. M87* ma silny accretion flow, Sgr A* relatively quiet.
   Korekta `δb_crit / b_crit` od astrophysical environment.
5. **ngEHT falsification window** (T5): przy 1% precision obu masy
   BH, distance, i ring diameter, w którym roku TGP M9.1'' zostanie
   detected/falsified vs GR? Final verdict OP-EHT.

---

## 3. Plan testów OP-EHT

| # | Test | Cel | Metoda | Status |
|---|------|-----|--------|--------|
| **T1** | PN truncation robustness | Pokazać czy +14.6% jest robust przy n=10, 12, 15 PN terms | sympy: extend recurrence c_n; numerical: compare b_crit at n=7,10,12,15 | **OPEN** |
| **T2** | Mass calibration MC | Monte Carlo over (M_BH, D, ring modeling) → effective EHT envelope; czy TGP mieści się z 95% CL? | numpy MC, GRAVITY/EHT priors | **OPEN** |
| **T3** | Q-renormalization audit | Czy G(Φ) = G_0/ψ w strong-field przesuwa A? Recompute b_crit z full Φ-dependent matching | sympy: A_renorm(ψ); numerical re-run photon ring | **OPEN** |
| **T4** | Non-vacuum environment | Accretion disk + plasma RT correction; estymować δb_crit dla M87* (high accretion) i Sgr A* (low) | analitycznie z literature; ray-tracing tests | **OPEN** |
| **T5** | ngEHT 2030+ falsification + final verdict | Określić rok detekcji/falsyfikacji TGP M9.1'' przy ngEHT 1% precision; final OP-EHT verdict | analiza priors ngEHT | **OPEN** |

**Kryterium sukcesu OP-EHT POSITIVE:**
- T1 pokazuje +14.6% jest artefaktem (zbiega do <5% przy n>10) **OR**
- T2 pokazuje effective envelope absorbuje +14.6% (przy current uncertainties) **OR**
- T3 + T4 łącznie redukują b_crit do < 5% deviation (consistent z EHT M87* i Sgr A*)

**Kryterium sukcesu OP-EHT CONDITIONAL POSITIVE:**
- T1-T4 fail to reduce deviation, ale T5 pokazuje że current EHT
  systematyka nie odróżnia TGP od GR przy 90% CL
- TGP staje się falsifiable w ngEHT 2030+ at 1% level — explicit
  prediction window

**Kryterium falsyfikacji OP-EHT:**
- T1 robust to n=15 (deviation +14.6% genuine)
- T2 effective envelope ≤ 8% (Sgr A* tension stays)
- T3 q-renormalization daje < 2% poprawkę (insufficient)
- T4 environment correction nie przekracza 3%
- T5 EHT 2024 już discriminates → TGP M9.1'' falsified strong-field
- ⇒ M9.2 pivot jest mandatory

---

## 4. Pliki OP-EHT (planowane)

- `OP_EHT_setup.md` — ten plik.
- `op_eht_t1_pn_truncation.py` — T1 implementacja (sympy: extend
  PN recurrence c_n; numpy: re-run photon ring at n=7,10,12,15).
- `op_eht_t1_pn_truncation.txt` — raw output T1.
- `OP_EHT_T1_results.md` — werdykt T1.
- `op_eht_t2_mass_calibration_mc.py` — T2 (numpy MC: M_BH, D,
  ring modeling priors → effective envelope).
- `op_eht_t2_mass_calibration_mc.txt` — raw output T2.
- `OP_EHT_T2_results.md` — werdykt T2.
- `op_eht_t3_q_renormalization.py` — T3 (sympy: A_renorm(ψ) z
  G(Φ) = G_0/ψ; numerical re-run).
- `op_eht_t3_q_renormalization.txt` — raw output T3.
- `OP_EHT_T3_results.md` — werdykt T3.
- `op_eht_t4_environment.py` — T4 (accretion disk + plasma δb_crit).
- `op_eht_t4_environment.txt` — raw output T4.
- `OP_EHT_T4_results.md` — werdykt T4.
- `op_eht_t5_ngeht_window.py` — T5 (ngEHT priors → falsification year).
- `op_eht_t5_ngeht_window.txt` — raw output T5.
- `OP_EHT_T5_results.md` — werdykt T5 + final OP-EHT verdict.

---

## 5. Cross-references

- [[research/op7/EHT_quick_results.md]] — punkt wyjścia +14.6%.
- [[tooling/scripts/gravity/eht_photon_ring_m911.py]] — base script do reuse.
- [[research/op-newton-momentum/M9_1_pp_P3_results.md]] —
  M9.1'' P3 audit (3 PASS, 1 conditional tension RESOLVED via OP-7,
  1 OPEN = OP-EHT, 0 falsified).
- [[research/op7/OP7_T1_results.md]] — σ_ab = 0 dla static spherical.
- [[paper/tgp_core.tex]] §applications "Black-hole shadows" + F4 falsifiability.
- [[KNOWN_ISSUES.md]] M9.1'' P3 strong-field section.

---

## 6. Bottom line — OP-EHT priority po OP-7 closure

OP-EHT jest **najwyższym priorytetem** w post-OP-7 era:
- **Niezależny od OP-7** (T1 confirms σ_ab = 0 dla static spherical)
- **Już ma robust prediction** (+14.6% z 5/6 PASS EHT-quick)
- **Bezpośrednio testowalny** (M87* + Sgr A* dane już istnieją)
- **Najbliższy do falsyfikacji** lub akceptacji TGP w strong-field

### Path forward post-OP-EHT:

1. **T1-T4 batch run** (numerical, możliwe w 1-2 dni)
2. **T5 verdict** (ngEHT timeline analysis)
3. **Section 7 OP-EHT row update** w paper (paralelnie do OP-7 row)
4. **Falsifiability section update** (F4 → strong-field shadow,
   ngEHT 2030+)
5. **M9.2 conditional opening** — tylko jeśli OP-EHT zamknie NEGATIVE
