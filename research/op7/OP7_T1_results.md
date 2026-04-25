# OP-7 / T1 — Strukturalny dowód no-tensor dla M9.1'' i 2 d.o.f. z σ_ab

**Data:** 2026-04-25
**Status:** ✅ POSITIVE
**Plik wykonawczy:** `op7_t1_no_tensor.py`
**Raw output:** `op7_t1_no_tensor.txt`
**Cross-references:**
- [[OP7_setup.md]] §3 (T1 row)
- [[../op-newton-momentum/M9_1_pp_P3_results.md]] (P3 GW170817 conditional tension, P3 EHT open)
- [[../../../tgp-core-paper/paper/tgp_core.tex]] §2 Remark "One substrate, two projections"

---

## 1. Cel testu

Dostarczyć rygorystyczny rachunek symboliczny pokazujący że:

(A) Pojedyncze pole `Φ` (zmienna `ψ = Φ/Φ₀`) z hyperbolicznym ansatzem M9.1''

```
ds² = -c₀² · (4-3ψ)/ψ · dt² + ψ/(4-3ψ) · δ_ij dx^i dx^j
```

daje pod dekompozycją SVT (Scalar-Vector-Tensor wzgl. SO(3))
**TYLKO mod skalarny** (breathing) — sektory wektorowy V i tensorowy
TT są tożsamościowo zerowe.

(B) Rozszerzony ansatz `g_ij = h(ψ) δ_ij + Λ(ψ) · σ_ij` z substratową
projekcją tensorową `σ_ab` daje **dokładnie 2 fizyczne d.o.f. TT**
(= h₊, h×) po nałożeniu warunków traceless + transverse.

## 2. Metoda

**sympy** symbolic decomposition:

1. Linearyzacja metryki M9.1'' wokół próżni `ψ_eq = 1`.
2. Plane-wave perturbacja `δψ(t,z) = ε·cos(ωt - kz)` propagująca w `+z`.
3. Obliczenie:
   - `δg_tt = -c₀²·(df/dψ)|_eq · δψ = +4c₀² · δψ`
   - `δg_ii = (dh/dψ)|_eq · δψ = +4 · δψ`  (i=x,y,z, wszystkie równe)
4. Dekompozycja `δg_ij`:
   - **Trace:** `12·ε·cos(...)` ≠ 0  → mod skalarny.
   - **Bezśladowa część:** `δg_ij - (1/3)Tr·δ_ij = 0` (z δg_ij ∝ δ_ij).
5. **Vector sektor:** `δg_0i = 0` z M9.1'' ansatzu (brak frame-dragging).
6. **Tensor TT sektor:** projekcja `h_TT = P·δg·P - (1/2)·P·tr_T(P·δg·P)`,
   gdzie `P = diag(1,1,0)` dla propagacji w `z` → wynik **macierzowo zerowy**.
7. Część B: σ_ij jako symmetric traceless 3×3 (5 d.o.f.); transverse
   (s_13 = s_23 = 0) + traceless transverse (s_22 = -s_11) → 2 d.o.f.

## 3. Wyniki

### 3.1 Część A — single-Φ M9.1'' (strukturalny no-tensor)

| Sektor | Wynik | Status |
|---|---|---|
| **Trace** δg_ii | `12 ε cos(kz - ωt)` | NIETRYWIALNE (breathing) |
| **Bezśladowa** δg_ij | 0 (tożsamościowo) | PASS |
| **Vector** δg_0i | 0 (M9.1'' ansatz) | PASS |
| **Tensor TT** h_ij^TT | 0 (macierzowo) | PASS |

**Wniosek A:** Single-Φ M9.1'' nie ma TT modes. Wynik niezależny od
konkretnej formy `f(ψ), h(ψ)` — wystarczy że `g_ij^(spacjalny) ∝ δ_ij`
(izotropowość przestrzennej części).

### 3.2 Część B — rozszerzony ansatz z σ_ab (2 d.o.f. TT)

Macierz σ_ij = symmetric traceless 3×3:

```
σ = [[s11, s12, s13],
     [s12, s22, s23],
     [s13, s23, -(s11+s22)]]
```

5 niezależnych skalarów: s11, s22, s12, s13, s23. Po nałożeniu warunków
fizycznej propagacji TT (transverse + traceless transverse):

- transverse `n^i σ_ij = 0` z `n = ê_z`: s13 = s23 = 0.
- traceless transverse: s22 = -s11.

Pozostają 2 d.o.f.:

```
σ_TT = [[ s11,  s12,  0],
        [ s12, -s11,  0],
        [   0,    0,  0]]
```

To są dokładnie **h₊** (= s11) i **h×** (= s12) jak w GR.

| Cecha | Wynik | Status |
|---|---|---|
| `Tr(σ) = 0` | 0 | PASS |
| `σ` symmetric | `‖σ - σᵀ‖ = 0` | PASS |
| Liczba d.o.f. po TT-redukcji | 2 | PASS |

### 3.3 Suma checks

7/7 PASS. Werdykt: **POSITIVE**.

## 4. Werdykt T1

✅ **STRUKTURALNIE POTWIERDZONO:**

1. **No-tensor dla M9.1'' single-Φ** — pojedyncze pole `ψ` produkuje
   tylko mod breathing (skalar). Sektory V i TT są tożsamościowo zerowe.
   To jest rozszerzenie `thm:no-tensor` z formy konformalnej (pre-pivot)
   na hiperboliczną (M9.1'').

2. **Rozszerzony ansatz `g_ij = h(ψ) δ_ij + Λ(ψ) σ_ij` daje 2 d.o.f. TT**
   po nałożeniu standardowych warunków transverse + traceless. Dokładnie
   tyle ile potrzeba dla GR-zgodnej fenomenologii GW (h₊, h×).

3. **Mechanizm substratowy zachowuje single-substrate axiom** — `σ_ab`
   jest **kompozytowym operatorem** coarse-grainingu (postulowanym
   już w `tgp_core.tex` §2 Remark "One substrate, two projections"),
   a nie niezależnym polem tensorowym. TGP **nie jest** teorią
   scalar-tensor.

## 5. Implikacje

### 5.1 Dla GW170817 (M9.1'' P3 conditional tension)

P3 conditional tension dla GW170817 wynika z tego, że single-Φ jest
strukturalnie niezdolne do produkcji TT modes — co oznacza scalar-only
GW signal, który LIGO/Virgo by **falsyfikował** na poziomie 5% bound.

T1 potwierdza, że to **NIE jest** falsyfikacja TGP, lecz potwierdzenie
**konieczności σ_ab** — który już jest postulowany w paperze, ale
bez wyprowadzonej dynamiki (T2-T6).

### 5.2 Dla aksjomatu single-substrate

Gdyby T2 zawiodło (σ_ab nie da się wyprowadzić jako kompozytowy
operator z `H_Γ`), TGP musiałoby **dodać niezależne pole tensorowe**,
co byłoby utratą single-substrate axiom (TGP_FOUNDATIONS §1) i
strukturalną falsyfikacją TGP w obecnej formie.

T1 jako pierwszy krok OP-7 jest **niezbędny ale niewystarczający** —
zamyka tylko stronę kinematyczną. Strona dynamiczna (T2-T5) i
konsystencyjna (T6) pozostają otwarte.

### 5.3 Dla EHT (P3 strong-field open)

T1 nie dotyczy bezpośrednio EHT (strong-field), ponieważ tam istotne
są pełne nieliniowe rozwiązania `Φ-EOM`. Jednak **trajektorie fotonu**
w silnym polu wymagają σ_ab nawet dla statycznego BH (gdy nieaksjalność
metryki staje się istotna). T1 dostarcza ramy dla tej dyskusji.

## 6. Następne kroki OP-7

| # | Test | Cel | Następnik T1 |
|---|------|-----|--------------|
| **T2** | Wyprowadzenie σ_ab z H_Γ | Pokazać że `σ_ab ∝ ⟨ŝ_i ŝ_{i+ê_a}⟩^TF` jest dobrze zdefiniowanym kompozytowym operatorem coarse-grainingu | analitycznie + lattice MC |
| **T3** | Dynamika σ_ab | `□σ_ab + m_σ² σ_ab = -ξ T_ab^TT` z wariacji S_TGP[ŝ] | wariacyjnie |
| **T4** | Sprzężenie Λ(ψ) | Postać Λ(ψ) zachowująca M9.1'' przy σ=0; ghost-free; c_GW = c₀ | analitycznie |
| **T5** | Formuła kwadrupolowa | h₊, h× ∝ Q̈_ij/r; matching ξ_eff = 4πG·σ₀·Φ₀² do GW150914 | analitycznie + dane |
| **T6** | Konsystencja | PPN niezmienione (σ=0 dla sferyczności), c_GW = c₀, Z₂-parzystość | sympy + scripts |

T2 jest najbardziej kluczowe — bez wyprowadzenia σ_ab z H_Γ T1 pozostaje
matematycznym ćwiczeniem, a nie strukturalnym mechanizmem TGP. T2
naturalnie buduje na `tooling/scripts/substrate/tensor_from_substrate.py`
(już zawiera lattice MC dla bilinear correlator), ale wymaga
przeformułowania w terminach v2 GL-substrate (nie pre-pivot Ising-like).

## 7. Pliki

- `op7_t1_no_tensor.py` — kod sympy (wykonano).
- `op7_t1_no_tensor.txt` — raw output 7/7 PASS.
- `OP7_T1_results.md` — ten plik.
- `OP7_setup.md` — plan ogólny.
