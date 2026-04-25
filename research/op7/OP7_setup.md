# OP-7 — sektor tensorowy TGP, mechanizm `σ_ab`

**Data otwarcia:** 2026-04-25
**Status:** OPEN (start)
**Następnik:** M9.1'' P3 (NOT-FALSIFIED, weak field) i P4 (paper integration).
**Punkt wyjścia ontologiczny:** `tgp-core-paper/paper/tgp_core.tex` §2,
`Remark: One substrate, two projections` (Φ ∝ ⟨ŝ²⟩, σ_ab ∝ ⟨ŝ·ŝ_{+â}⟩^TF).
**Aktualizacja po T2 (2026-04-25):** definicja σ_ab z H_Γ jest **gradientowym
naprężeniem** `σ_ab = K_ab − (1/3)δ_ab Tr(K)`, gdzie `K_ab = ⟨(∂_a ŝ)(∂_b ŝ)⟩_B`.
Postać `⟨ŝ·ŝ_{+â}⟩^TF` (nearest-neighbor) jest *równoważna do leading order
gradient expansion* na siatce, ale gradient-strain forma jest **kanoniczna**:
pochodzi z kinetycznego termu H_Γ, naturalnie traceless+symmetric, ma 5 d.o.f.
w 3D, znika w izotropowej próżni. Patrz [[OP7_T2_results.md]].
**Wcześniejsze pre-pivotowe rozważania:** `tooling/scripts/substrate/tensor_from_substrate.py`,
`tooling/scripts/gw/disformal_waveform.py`, `tooling/scripts/gw/gw_breathing_mode.py`.

---

## 1. Cel OP-7

Zamknąć **strukturalnie** sektor tensorowy GW w TGP — czyli udzielić
dowodu (lub falsyfikacji), że pojedyncze pole skalarne `Φ` widziane przez
hiperboliczną metrykę M9.1'' wraz z jego kompozytową projekcją tensorową
`σ_ab` z substratu reprodukuje **dwie polaryzacje TT** GW (h+, h×) z
amplitudą zgodną z GR i ograniczeniami LIGO/Virgo/KAGRA.

**Pytanie centralne:**

> Czy dla metryki ansatz `g_μν[ψ, σ_ab]` rozszerzającej M9.1'' o substratową
> projekcję tensorową, perturbacje na próżni dają **dwie** propagujące mody
> TT, z prędkością c₀ i amplitudą zgodną z formułą kwadrupolową?

Jeśli TAK — TGP zamyka GW prediction zgodnie z GR observation.
Jeśli NIE — TGP jest falsyfikowane na poziomie GW170817/GW150914.

## 2. Tło — co już wiemy

### 2.1 Twierdzenie no-tensor (M9.1'')

Dla M9.1'' hyperbolic ansatz:

```
ds² = -c₀²·(4-3ψ)/ψ·dt² + ψ/(4-3ψ)·δ_ij dx^i dx^j,    ψ = Φ/Φ₀
```

metryka jest **diagonalna i konformalna**, więc perturbacja `ψ → ψ + δψ`
daje:

```
δg_tt = (∂_ψ f)·δψ = (4c₀²/ψ²)·δψ        [skalar]
δg_ii = (∂_ψ h)·δψ = [4/(4-3ψ)²]·δψ      [skalar, izotropowy]
δg_ij^TT ≡ 0                              [brak TT]
```

Pod dekompozycją SVT (Scalar-Vector-Tensor) względem SO(3): **tylko
breathing mode**. Vector i Tensor sektory są tożsamościowo zerowe.

To jest stara `thm:no-tensor` przeniesiona do nowej formy hiperbolicznej.

### 2.2 σ_ab jako kompozytowa projekcja substratu

Już w `tgp_core.tex` §2 jest stwierdzone (jako postulat coarse-grainingu):

```
Φ      ∝ ⟨ŝ²⟩_B                          (skalar, izotropowy)
σ_ab   ∝ ⟨ŝ_i · ŝ_{i+ê_a}⟩_B^TF           (tensor, traceless+symmetric)
```

`σ_ab` jest **kompozytową** projekcją tego samego pola `ŝ`, więc TGP
**nie staje się** teorią scalar-tensor — pozostaje jednoskładnikową
teorią substratową, z dwiema projekcjami efektywnymi.

`tensor_from_substrate.py` pokazał numerycznie:

- σ_ab ma 5 niezależnych stopni swobody (5 składowych spin-2);
- σ_ab = 0 w próżni izotropowej (warunek brzegowy);
- σ_ab ≠ 0 w obecności anizotropowego źródła (kwadrupol masy).

### 2.3 Co OP-7 musi dostarczyć

1. **Dowód strukturalny** (T1): no-tensor dla M9.1'' (formalnie).
2. **Definicja σ_ab** (T2): z H_Γ, jako kompozytowy operator coarse-grainingu.
3. **EOM dla σ_ab** (T3): wariacyjnie z S_TGP[ŝ], z masą m_σ i sprzężeniem do Q_ab^TT.
4. **Sprzężenie do metryki** (T4): rozszerzony ansatz `g_μν[ψ, σ_ab]`,
   który przy σ=0 redukuje do M9.1'' i przy σ≠0 produkuje TT modes.
5. **Amplituda GW** (T5): formuła kwadrupolowa, dopasowanie do GW150914.
6. **Konsystencja** (T6): PPN niezmienione, c_GW = c₀, ghost-free, Z₂.

---

## 3. Plan testów OP-7

| # | Test | Cel | Metoda | Status |
|---|------|-----|--------|--------|
| **T1** | No-tensor dla M9.1'' | Dowód, że δg_eff[δψ] ma tylko mod S | sympy: SVT decomposition perturbacji | **POSITIVE 2026-04-25** ([[OP7_T1_results.md]]) |
| **T2** | Definicja σ_ab z H_Γ | Pokazać, że gradient strain `K_ab - (1/3)δ_ab Tr(K)` jest dobrze zdefiniowanym operatorem coarse-grainingu | analitycznie + sympy + lattice MC | **POSITIVE 2026-04-25** ([[OP7_T2_results.md]]) |
| **T3** | Dynamika σ_ab | Wyprowadzić □σ_ab + m_σ²σ_ab = -ξ T_ab^TT z S_TGP | wariacyjnie | **POSITIVE (44/47=94%) 2026-04-26**: T3.1-T3.4 ([[OP7_T3_results.md]]) + T3.5-T3.6 ([[OP7_T3_extended_results.md]]); Φ₀/m_σ tension RESOLVED przez decoupling |
| **T4** | Metryka rozszerzona | Postać `g_ij = h(ψ)δ_ij + Λ(ψ)·σ_ij`; zachowanie M9.1'' przy σ=0 | analitycznie | **POSITIVE (13/13=100%) 2026-04-25** ([[OP7_T4_results.md]]): Λ(ψ)=const=1 strukturalnie unikalne; scenario A (decoupling) RATIFIED |
| **T5** | Formuła kwadrupolowa | h_+, h_× ∝ Q̈_ij/r; dopasowanie do GW150914 | analitycznie + dane LIGO | **partial (T3.4) 2026-04-25**: ξ/G ≈ 1.06 |
| **T6** | Konsystencja | PPN niezmienione, c_GW = c₀, ghost-free, Z₂ | sympy + skrypty | **partial (T3.3) 2026-04-25**: ghost-free PASS; pełne PPN open |

**Kryterium sukcesu OP-7:**
- T1 ✓ (strukturalny no-tensor M9.1'')
- T2-T4 dają spójną definicję + dynamikę σ_ab
- T5 produkuje 2 mody TT z amplitudą GR (do współczynnika ξ_eff
  dopasowywanego do G)
- T6 wszystkie passy

**Kryterium falsyfikacji OP-7:**
- T2 lub T3 nie dają dobrze zdefiniowanego operatora — TGP musiałby
  rozszerzyć aksjomat o niezależne pole tensorowe (tj. zostać
  scalar-tensor, łamiąc TGP_FOUNDATIONS §1).
- T5 daje amplitudę różniącą się od GR > LIGO bound (typowo 5%).

---

## 4. Pliki OP-7

- `OP7_setup.md` — ten plik.
- `op7_t1_no_tensor.py` — T1 implementacja (sympy SVT decomp).
- `op7_t1_no_tensor.txt` — raw output T1.
- `OP7_T1_results.md` — werdykt T1.
- `op7_t2_sigma_from_HGamma.py` — T2 implementacja (sympy + lattice MC).
- `op7_t2_sigma_from_HGamma.txt` — raw output T2 (12/12 PASS).
- `OP7_T2_results.md` — werdykt T2.
- `op7_t3_sigma_dynamics.py` — T3.1 sympy structural EOM derivation (Path A + Path B).
- `op7_t3_sigma_dynamics.txt` — raw output T3.1 (11/11 PASS).
- `op7_t3_2_m_sigma_scale.py` — T3.2 m_σ hipotezy A/B/C + GW170817 bound.
- `op7_t3_2_m_sigma_scale.txt` — raw output T3.2 (4/7 PASS, Φ₀/m_σ tension).
- `op7_t3_3_ghost_analysis.py` — T3.3 ghost+dispersion+massless analysis.
- `op7_t3_3_ghost_analysis.txt` — raw output T3.3 (5/5 PASS).
- `op7_t3_4_xi_coupling.py` — T3.4 ξ coupling + GW150914 matching.
- `op7_t3_4_xi_coupling.txt` — raw output T3.4 (5/5 PASS, ξ/G ≈ 1.06).
- `OP7_T3_results.md` — synteza werdyktu T3 (mixed: 25/28 = 89% PASS, structural pos with Φ₀/m_σ tension).
- `op7_t3_5_bethe_salpeter.py` — T3.5 spectral + Bethe-Salpeter analysis.
- `op7_t3_5_bethe_salpeter.txt` — raw output T3.5 (10/10 PASS).
- `op7_t3_6_symmetry_protection.py` — T3.6 symmetry/decoupling/ULDM scenarios.
- `op7_t3_6_symmetry_protection.txt` — raw output T3.6 (9/9 PASS).
- `OP7_T3_extended_results.md` — synteza T3-extended (19/19 PASS; tension RESOLVED przez decoupling).
- `op7_t4_metric_coupling.py` — T4 implementacja (sympy ansatz analysis + 5 candidates Λ(ψ)).
- `op7_t4_metric_coupling.txt` — raw output T4 (13/13 PASS).
- `OP7_T4_results.md` — werdykt T4 (Λ=const=1 unique; scenario A RATIFIED).
- `eht_photon_ring_m911.py` (w tooling/scripts/gravity/) — EHT-quick analysis.
- `eht_photon_ring_m911.txt` — raw output EHT (5/6 PASS).
- `EHT_quick_results.md` — werdykt EHT (INCONCLUSIVE-leaning-NEGATIVE: +14.6% photon ring; Sgr A* tension).
- `TGP_CLOSURE_PLAN_2026-04-25.md` — master closure plan (OP-7 + parallel programs).
- (kolejne: T5, T6 jeszcze open; T4 zamknięte 2026-04-25).

---

## 5. Cross-references

- `tgp-core-paper/paper/tgp_core.tex` §2 (Remark: One substrate,
  two projections) — postulat σ_ab.
- `tgp-core-paper/paper/tgp_core.tex` §7 status table OP-7 row.
- `TGP/TGP_v1/research/op-newton-momentum/M9_1_pp_P3_results.md`
  — GW170817 conditional tension, EHT open, oba zależne od OP-7.
- `TGP/TGP_v1/tooling/scripts/substrate/tensor_from_substrate.py`
  — wcześniejsza heurystyka σ_ab (pre-M9.1'', wymaga rewrite).
- `TGP/TGP_v1/tooling/scripts/gw/disformal_waveform.py`
  — pre-pivotowy disformal mechanism (porażka na 18 rzędów).
- `TGP/tgp-core-paper/KNOWN_ISSUES.md` C4 (LIGO 5% bound, GW170817).

---

## 6. Bottom line

OP-7 jest **najbardziej krytycznym** otwartym programem operacyjnym TGP
po M9.1''. Jego wynik zamyka pytania:

- Czy TGP jest falsyfikowane przez **GW170817** (skalar-only nie pasuje
  do dwupolaryzacyjnego strain)?
- Czy TGP może mieć **emergentny spin-2** zachowując pojedyncze
  pole `Φ` na poziomie aksjomatu?
- Czy ax:metric-from-potential **wymaga uogólnienia** o σ_ab?

T1 zaczyna od najkrótszego, najpewniejszego ogniwa — strukturalnego
dowodu no-tensor dla M9.1''. Pozytywny wynik T1 nie jest "sukcesem" w
sensie GW — jest **potwierdzeniem motywacji** dla σ_ab.
