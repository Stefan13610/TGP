---
title: "op-Csigma-lattice-MC — liczbowe wyznaczenie T=C_σσ_0² (tensor stiffness) z kierunkowego bąbla ⟨O_ab O_cd⟩ na sieci 3D Ising (MC/FRG), by przypiąć κ_E=8πG_0T/c³ liczbowo i wydać TWARDY werdykt sektora radiacyjnego (FALSIFIED / SURVIVE / GAP). Klasa dodatekQ CG-2, wielosesyjny."
type: research_cycle
status: "🟡 CLOSED-RESOLVED — PARTIAL (UNDERDETERMINED-fine-tuned, liczbowy κ_E≈0.62 O(1); lean FALSIFIED) 2026-06-14"
phase: FINAL
folder_status: active
created_date: 2026-06-14
closed_date: 2026-06-14
registered_by: "user 2026-06-14 (sesja #31): 'zarejestruj op-Csigma-lattice-MC'"
executed_by: "user 2026-06-14 (sesja #31): 'przeprowadzić cykl badawczy' = działaj (Phase 0→FINAL)"
verdict: "PARTIAL — κ_E≈0.62=O(1), pasmo [0.04,11.1] obejmuje 5/6 i 1; lean strukturalny FALSIFIED; residual GAP = scheme-independent continuum operatora złożonego (dodatekQ CG-3/CG-4, analityczny)"
phases_done: "Phase0_lock · Phase1_setup (pipeline 4/4 PASS) · Phase2_bubble (RDZEŃ, C_σ>0 O(1), 5/5 PASS) · Phase3_unitbridge (κ_E O(1), 4/4) · Phase_FINAL_close"
scoping: "[[../../meta/SCOPING_op-Csigma-lattice-MC_2026-06-14.md]]"
parent_cycle: "[[../op-Csigma-coarse-graining-2026-06-14/]] (PARTIAL; UNDERDETERMINED-fine-tuned węższy; residual GAP = wartość T)"
grandparent_cycle: "[[../op-sigma-kinetic-Csigma-2026-06-14/]] (UNDERDETERMINED-fine-tuned; C_σ=GAP)"
target_open_problem: "rem:sigma-params (wartość T niezobliczona) + rem:param-counting (3→2 LICZBOWE); dodatekQ CG-2/CG-3/CG-4"
cycle_category: "DERIVATION-NUMERICAL (lattice-MC kierunkowego bąbla; wielosesyjny)"
expected_duration: "wielosesyjny (klasa dodatekQ CG-2); DERIVED/PARTIAL/GAP wszystkie pełnoprawne"
object_id: "T = C_σσ_0² (JEDEN parametr, redundancja przeskalowania); κ_E = 8πG_0T/c³"
anti_lakatos_lock: "INHERITED; aktywny od Phase 0"
---

# op-Csigma-lattice-MC (🟡 CLOSED-RESOLVED — PARTIAL, lean FALSIFIED)

> **Status:** WYKONANY (Phase 0→FINAL, sesja #31, 2026-06-14). Numeryczny (3D Ising, Swendsen-Wang).
> **Werdykt: PARTIAL** — $\kappa_E\approx0.62=O(1)$ (pasmo $[0.04,11.1]$ obejmuje 5/6 i 1), sektor radiacyjny
> **UNDERDETERMINED-fine-tuned (teraz LICZBOWY)**, **lean strukturalny FALSIFIED**. $C_\sigma>0$, $O(1)$
> (~0.5–0.7) **ZMIERZONE**. Residual GAP **przesunięty/zidentyfikowany**: scheme-independent continuum
> operatora złożonego (R-continuum: power-divergence + obstrukcja krytyczna $p_{\min}/m_s>1$) → wymaga
> **analitycznej** renormalizacji (dodatekQ **CG-3/CG-4**), nie MC. Pełny werdykt: [[./Phase_FINAL_close.md]].

## Pliki cyklu (wykonane)
- [[./Phase0_lock.md]] — LOCK falsyfikatorów F-LMC-A..D, próg 5/6, forbidden moves (anti-Lakatos aktywny).
- [[./Phase1_setup.md]] + `Phase1_lattice.py` — setup + walidacja pipeline (**4/4 PASS**: maszyneria <1%, continuum dev 4.8%).
- [[./Phase2_bubble.md]] + `Phase2_mc.py` — RDZEŃ: MC bąbla, $C_\sigma>0$ $O(1)$ zmierzone (**5/5 PASS**); F-LMC-A=PARTIAL.
- [[./Phase3_unitbridge.md]] + `Phase3_unitbridge.py` — unit-bridge → $\kappa_E\approx0.62$ (**4/4**); F-LMC-C=PARTIAL.
- [[./Phase_FINAL_close.md]] + `Phase_FINAL_sympy.py` — agregat value-blind → **PARTIAL, lean FALSIFIED**.

## Pytanie wiodące
Czy $T=C_\sigma\sigma_0^2$ jest wyznaczalne **liczbowo** z kierunkowego bąbla
$\Pi_{ab,cd}(p)=\langle O_{ab}(p)O_{cd}(-p)\rangle_{\rm c}$, $O_{ab}=(\partial_a\hat s\,\partial_b\hat s)_{\rm TF}$,
na sieci 3D Ising — a złożone w $\kappa_E=8\pi G_0T/c^3$ daje $5/6$ (SURVIVE) czy $\neq$ (FALSIFIED)?
Survival = miara zero, **niechroniona** (parent op-Csigma-coarse-graining); naturalna κ_E=1 → 7/6 FALSIFIED.

## Trasa (fazy; pełny opis w scopingu §2)
1. **Phase 1** — setup sieci 3D Ising/GL (eq:B-H, parametry dodatekQ CG-2) + operator $O_{ab}$ + obserwable $\Pi(p)$; algorytm klastrowy.
2. **Phase 2** — pomiar bąbla, ekstrakcja $C_\sigma$ ze współczynnika $p^2$ (RDZEŃ numeryczny); ekstrapolacja continuum + FSS; walidacja vs bąbel analityczny (parent Phase 2).
3. **Phase 3** — unit-bridge ($G_0\sim J\mu$, $a_\Gamma\Phi_0=1$) → $T$ w jednostkach $c^3/G_0$ → $\kappa_E$ z błędem stat.+syst.
4. **Phase FINAL** — $\kappa_E$ vs 5/6 (próg LOCKED Phase 0) → **FALSIFIED hard** / SURVIVE / PARTIAL / GAP.

## Falsyfikatory (szkic; LOCK w Phase 0)
F-LMC-A (bąbel p²>0, continuum-zbieżny?) · F-LMC-B (T w c³/G_0, unit-bridge) · F-LMC-C (κ_E vs 5/6) · F-LMC-D (agregat).

## Reuse
dodatekQ CG-2 (`tgp_erg_lpa_prime.py` 8/8 PASS; ρ₀*=0.03045, ν=0.749) — wzorzec+parametry; op-Csigma-coarse-graining Phase 2 (bąbel analityczny — walidacja) + Phase 3 (T=C_σσ_0², unit-bridge); dodatekQ Q.4 (a_Γ·Φ_0=1).

## Twarde wymogi / ryzyka (scoping §5/§6)
- Zakaz strojenia T do 5/6 (value-blind); zakaz liczenia C_σ,σ_0 osobno (tylko niezmienniczy T); zakaz fabrykowania wartości (GAP/PARTIAL jawnie).
- **Systematyki raportowane z błędami** (continuum, FSS, renormalizacja operatora złożonego, unit-bridge).
- R-continuum (HIGH): renormalizacja bąbla operatora ZŁOŻONEGO (dodatekQ CG-3/CG-4 otwarte). R-unit-bridge (HIGH): systematyka może dać tylko PARTIAL. R-critical-slowing (MED). R-tensor-projection (MED).

## Oczekiwany wynik (INFORMATIONAL)
Najprawdopodobniej **DERIVED → κ_E=O(1)≠5/6 → FALSIFIED hard** (lean parent: brak symetrii chroniącej 5/6).
Możliwe PARTIAL (przedział obejmujący 5/6 i 1) lub GAP (continuum niezbieżne). Mało prawdopodobne κ_E=5/6 SURVIVE.

## Aktywacja
Nowy agent: czytaj [[../../meta/SCOPING_op-Csigma-lattice-MC_2026-06-14.md]] + parent
([[../op-Csigma-coarse-graining-2026-06-14/Phase_FINAL_close.md]], Phase 2 bąbel + Phase 3 T) +
dodatekQ CG-2 (`tgp_erg_lpa_prime.py`, parametry punktu stałego) + eq:B-H / def:sigma-ab.
Potem user „działaj" → Phase 0 LOCK.
