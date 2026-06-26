---
title: "op-Csigma-coarse-graining — wyprowadzenie C_σ (i σ_0) sektora σ_ab z H_Γ przez coarse-graining, by przypiąć κ_E i rozstrzygnąć sektor radiacyjny (FALSIFIED / SURVIVE / GAP). Klasa op-gamma-RG-running, wielosesyjny."
type: research_cycle
status: 🟡 CLOSED-RESOLVED — PARTIAL (UNDERDETERMINED-fine-tuned, węższy; lean FALSIFIED) 2026-06-14
phase: FINAL
folder_status: closed
close: "[[Phase_FINAL_close.md]] + [[Phase_FINAL_sympy.py]]/[[Phase_FINAL_sympy.txt]] — agregat F-CG-E value-blind (reguła LOCKED): PARTIAL; lean FALSIFIED; spawn op-Csigma-lattice-MC (liczbowe T)"
created_date: 2026-06-14
registered_by: "user 2026-06-14 (sesja #31): 'rozpisz ten cykl'"
activated_by: "user 2026-06-14 (sesja #31): 'działaj' (faza 1); 'działaj z fazą 2'; 'Działaj z fazą 3'"
phase0: "[[Phase0_balance.md]] 🔒 LOCKED — F-CG-A..E + próg 5/6 (value-blind)"
phase1: "[[Phase1_derivation.md]] + [[Phase1_sympy.py]]/[[Phase1_sympy.txt]] — σ_ab=KOMPOZYT bilinowy; F-CG-A=PARTIAL-pending; F-CG-B/C/D=OPEN; C_σ NIE wyprowadzone (anti-Lakatos)"
phase2: "[[Phase2_derivation.md]] + [[Phase2_sympy.py]]/[[Phase2_sympy.txt]] — metoda=propagator kompozytu (bąbel); C_σ>0+skaling DERIVED, prefaktor=GAP; κ_E=8πG_0C_σσ_0²/c³; BRAK symetrii lockującej→κ_E O(1), 5/6 miara zero. PARTIAL, lean FALSIFIED"
phase3: "[[Phase3_derivation.md]] + [[Phase3_sympy.py]]/[[Phase3_sympy.txt]] — REDUNDANCJA σ→λσ (sympy 3/3): C_σ i σ_0 = JEDEN parametr T=C_σσ_0² (uzasadnia rem:param-counting 3→2). κ_E=8πG_0T/c³; T niezlockowane (brak diff-inwariancji)→lean FALSIFIED. F-CG-C=RESOLVED+PARTIAL"
scoping: "[[../../meta/SCOPING_op-Csigma-coarse-graining_2026-06-14.md]]"
parent_cycle: "[[../op-sigma-kinetic-Csigma-2026-06-14/]] (UNDERDETERMINED-fine-tuned; brama = C_σ z H_Γ)"
target_open_problem: "sek08 rem:sigma-params / rem:param-counting — C_σ 'niezobliczony; problem otwarty'; redukcja 3→2 parametrów TGP"
cycle_category: "DERIVATION-FROM-SUBSTRATE (coarse-graining sektora tensorowego; wielosesyjny)"
expected_duration: "wielosesyjny (klasa op-gamma-RG-running ~10+ faz); GAP/PARTIAL/DERIVED wszystkie pełnoprawne"
object_id: "C_σ = stała kinetyczna σ_ab (eq:S-sigma); κ_E = C_σσ_0²/(GR-lock); NOŚNIK GW (poprawny obiekt)"
anti_lakatos_lock: "INHERITED; aktywny od Phase 0"
---

# op-Csigma-coarse-graining (🟡 CLOSED-RESOLVED — PARTIAL, lean FALSIFIED)

> **Werdykt końcowy (value-blind, reguła LOCKED Phase 0):** **PARTIAL** ⟹ sektor radiacyjny
> **UNDERDETERMINED-fine-tuned (STATUS WĘŻSZY)**, lean strukturalny **FALSIFIED**. Pełny close:
> [[Phase_FINAL_close.md]].
>
> **Co cykl ustalił:** (1) σ_ab = KOMPOZYT bilinowy ⟨ŝŝ⟩ (Phase 1); (2) kinetyka z **propagatora
> kompozytu (bąbel)** — $C_\sigma>0$, skaling DERIVED, prefaktor GAP; **brak Warda** ⟹ κ_E O(1)-bounded
> (Phase 2); (3) **redundancja przeskalowania** ⟹ $C_\sigma,\sigma_0$ = JEDEN parametr $T=C_\sigma\sigma_0^2$,
> uzasadnia `rem:param-counting` 3→2 (Phase 3). $\kappa_E=8\pi G_0T/c^3$; survival 5/6 = miara zero,
> niechroniona; naturalna κ_E=1 → 7/6 FALSIFIED.
>
> **Residual GAP (jedyny):** liczbowa wartość $T=C_\sigma\sigma_0^2$ → **lattice-MC kierunkowego bąbla**
> ⟨O_ab O_cd⟩ na sieci 3D Ising. **Spawn:** `op-Csigma-lattice-MC` (REGISTERED-QUEUED). To przeprowadza
> PARTIAL → DERIVED/FALSIFIED-hard. Anti-Lakatos: prefaktor NIE sfabrykowany; rdzeń nie edytowany.

## Pytanie wiodące
Czy $C_\sigma$ (sztywność kinetyczna σ_ab) jest wyprowadzalna z $H_\Gamma$ przez coarse-graining korelacji
kierunkowej $\sigma_{ab}=\langle\hat s_i\hat s_{i+\hat a_b}\rangle^{\rm TF}$ — a złożona w
$\kappa_E=C_\sigma\sigma_0^2/(\text{GR-lock})$ daje $5/6$ (SURVIVE) czy $\neq$ (FALSIFIED)?
Survival to **miara zero** (op-sigma-kinetic-Csigma); naturalna κ_E≈1 → gałąź B PR-025 (falsyfikacja).

## Trasa (fazy; pełny opis w scopingu §2)
1. **Phase 1** — formalne H_Γ dla σ_ab (człon kierunkowy eq:B-H; σ_ab jako kompozyt; inwentarz {J,a_sub,μ,...}).
2. **Phase 2** — coarse-graining → C_σ (RDZEŃ): gradient expansion, ekstrakcja współczynnika $(\partial\sigma)^2$; wymóg Lorentza (c_0).
3. **Phase 3** — σ_0 + złożenie κ_E (rozstrzygnąć C_σ vs σ_0 = jeden czy dwa parametry).
4. **Phase FINAL** — κ_E vs 5/6 (próg LOCKED Phase 0) → FALSIFIED / SURVIVE / GAP.

## Falsyfikatory (szkic; LOCK w Phase 0)
F-CG-A (emergentna kinetyka Lorentz?) · F-CG-B (forma C_σ) · F-CG-C (σ_0+κ_E) · F-CG-D (κ_E vs 5/6) · F-CG-E (agregat).

## Reuse
op-gamma-RG-running (wzorzec skalarny H_Γ→S[Φ]); dodatekQ_coarse_graining; GL-bond v2; def:sigma-ab; thm:amplitude-matching; op-sigma-kinetic-Csigma (wejście strukturalne).

## Twarde wymogi / ryzyka (scoping §5/§6)
- Zakaz strojenia C_σ/σ_0 do 5/6 (value-blind); zakaz fabrykowania prefaktora (GAP/PARTIAL jawnie — lekcja sesji #30).
- R-tensor-hard (HIGH): coarse-graining tensorowy trudniejszy niż skalarny — możliwy PARTIAL/GAP (zgodnie z rem:sigma-params „otwarte").
- Wymóg Lorentza (c_0) na kinetykę σ_ab.

## Oczekiwany wynik (INFORMATIONAL)
Najprawdopodobniej PARTIAL/GAP (skaling bez prefaktora) → sektor pozostaje UNDERDETERMINED-fine-tuned,
ale węższy (jeśli skaling wyklucza 5/6 o rzędy → de facto FALSIFIED). Możliwe DERIVED→κ_E≠5/6→FALSIFIED
(definitywne domknięcie). Mało prawdopodobne κ_E=5/6 SURVIVE (spisek dwóch sektorów).

## Aktywacja
Nowy agent: czytaj [[../../meta/SCOPING_op-Csigma-coarse-graining_2026-06-14.md]] + cykl-rodzic
([[../op-sigma-kinetic-Csigma-2026-06-14/Phase1_derivation.md]]) + rdzeń (def:sigma-ab, eq:S-sigma,
thm:amplitude-matching, eq:B-H l.6700, rem:sigma-params) + wzorzec op-gamma-RG-running-derivation.
Potem user „działaj" → Phase 0 LOCK.
