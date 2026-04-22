# `nbody/examples` — mapa statusów

Ten plik porządkuje aktywne skrypty w `nbody/examples/` według ich bieżącej roli
po synchronizacji `nbody` z rdzeniem teorii.

## Zasada ogólna

Numer `exNNN` **nie oznacza** sam z siebie statusu teoretycznego. Po synchronizacji
obowiązuje prostszy podział:

- `CANONICAL` — skrypty, na których opiera się aktualna warstwa `nbody`
- `ACTIVE-EXPLORATORY` — aktywne badania i rozwinięcia, ale nie punkt odniesienia
- `LEGACY-TRANSLATIONAL` — starszy język, ważny historycznie lub jako mapa przejścia
- `ARCHIVE` — materiały przeniesione do `_archiwum/`

Przy konflikcie pierwszeństwo mają:

1. `../THEORY_SYNC_NBODY.md`
2. `../ZALOZENIA_NBODY.md`
3. ten plik
4. dopiero potem pojedyncze docstringi `exNNN`

## 1. CANONICAL

To jest bieżący tor, który warto utrzymywać i na którym warto budować dalsze
analizy.

### A. Regresja i EOM N-ciał

- `verify_nbody_canonical_quick.py`
- `verify_nbody_bridge_extended.py`
- `verify_nbody_eom_quick.py`
- `ex138_eom_tgp_coulomb_leapfrog.py`
- `ex139_yukawa_feynman_leapfrog.py`
- `ex140_yukawa_feynman_N4_scaling_rk.py`
- `ex141_net_force_translational_invariance.py`
- `ex142_yukawa_overlap_quadrature_convergence.py`
- `ex143_net_torque_so3_invariance.py`
- `ex144_conserved_P_L_leapfrog.py`
- `ex145_minus_grad_V_matches_force.py`
- `ex146_yukawa_feynman_grad_matches_force.py`
- `ex147_com_acceleration_zero.py`

### B. Lyapunov / chaos / Hamiltonian consistency

- `verify_nbody_lyapunov_quick.py`
- `ex148_lyapunov_tgp_vs_newton_pairwise.py`
- `ex149_lyapunov_spectrum_newton.py`
- `ex150_lyapunov_pairwise_vs_coulomb3b.py`
- `ex151_lyapunov_scan_beta_pairwise.py`
- `ex152_lyapunov_yukawa_feynman_short.py`
- `ex153_lyapunov_spectrum_sum_newton.py`
- `ex154_lyapunov_scan_softening_csv.py`
- `ex155_lyapunov_grid_beta_softening_csv.py`
- `ex156_lyapunov_newton_rk4_vs_leapfrog.py`
- `ex157_lyapunov_tgp_pairwise_rk4_vs_leapfrog.py`
- `ex158_lyapunov_newton_seed_spread_leapfrog.py`
- `ex159_lyapunov_coulomb3b_rk4_vs_leapfrog.py`
- `ex160_lyapunov_yukawa_feynman_leapfrog_only.py`
- `ex161_lyapunov_pairwise_beta_d_rep_scan.py`
- `ex162_lyapunov_yukawa_feynman_fd_vs_split_jacobian.py`
- `ex163_lyapunov_yukawa_feynman_split_vs_analytic_jacobian.py`
- `ex164_lyapunov_yukawa_feynman_long_window_analytic.py`
- `ex165_lyapunov_yukawa_feynman_refine_dt_nquad.py`
- `ex166_lyapunov_yukawa_feynman_convergence_grid_csv.py`
- `ex167_leapfrog_energy_drift_yukawa_vs_coulomb3b_short.py`
- `ex168_plot_ex166_convergence_csv.py`
- `ex169_lyapunov_yukawa_feynman_leapfrog_fd_vs_analytic_jacobian.py`
- `ex170_lyapunov_yukawa_feynman_spectrum_sum_leapfrog.py`
- `ex171_lyapunov_coulomb3b_spectrum_sum_leapfrog.py`
- `ex172_lyapunov_pairwise_spectrum_sum_leapfrog.py`
- `ex173_lyapunov_pairwise_leapfrog_fd_vs_analytic_spectrum.py`
- `ex174_lyapunov_newton_vs_yukawa_feynman_matched_energy.py`
- `ex175_lyapunov_burrau_newton_yukawa_three_hamiltonians.py`
- `ex176_lyapunov_newton_vs_coulomb3b_matched_energy.py`
- `ex177_lyapunov_newton_vs_pairwise_matched_energy.py`
- `ex178_lyapunov_matched_energy_family_table.py`
- `ex179_lyapunov_burrau_newton_coulomb_three_hamiltonians.py`
- `ex180_lyapunov_burrau_newton_pairwise_three_hamiltonians.py`
- `ex181_lyapunov_matched_energy_seven_row_overview.py`
- `ex182_lyapunov_seven_row_overview_csv.py`
- `ex183_summarize_ex182_seven_row_csv.py`
- `ex184_lyapunov_seven_row_t_final_scan_csv.py`
- `ex185_summarize_ex184_t_scan_csv.py`
- `ex186_lyapunov_matched_h_family_csv.py`
- `ex187_summarize_ex186_matched_h_csv.py`
- `ex188_lyapunov_matched_h_family_t_final_scan_csv.py`
- `ex189_summarize_ex188_matched_h_t_scan_csv.py`
- `ex190_leapfrog_energy_drift_matched_h_family_short.py`
- `ex191_rk45_energy_diag_matched_h_family_short.py`
- `ex192_matched_h_leapfrog_vs_rk45_energy_table.py`
- `ex193_matched_h_lf_vs_rk45_energy_csv.py`
- `ex194_summarize_ex193_matched_h_lf_vs_rk_csv.py`
- `ex198_lyapunov_p1_closure.py`
- `ex199_lyapunov_beta_scan.py`
- `ex200_lyapunov_beta_scan_yukawa.py`
- `ex207_lyapunov_multi_ic_beta_scan_v3.py`
- `ex208_lyapunov_spectrum_pairing.py`

### C. Most klasyczny -> EFT -> N-body

- `ex195_soliton_ksub_full_vs_lpa.py`
- `ex196_phi_fp_full_vs_lpa.py`
- `ex197_optimal_g0_full_form.py`
- `ex205_path_c_yukawa_from_defect.py`
- `ex206_soliton_interaction_energy.py`

Te pliki są dziś kluczowe dla interpretacji:

- klasyczny defekt ma ogon `sin(r)/r`
- źródło Yukawy w `nbody` jest obiektem EFT
- FULL `K_sub = g^2` jest kanoniczne

`verify_nbody_bridge_extended.py` jest celowo węższy niż
`verify_nbody_lyapunov_quick.py`: ma spinać nową ścieżkę mostu jednym stabilnym
runnerem, a nie uruchamiać cały cięższy pakiet P1.

## 2. ACTIVE-EXPLORATORY

To są aktywne eksploracje lub rozwinięcia, ale nie należy ich traktować jako
jednej kanonicznej warstwy `nbody`.

### A. Ścieżka 9 / A_tail / Koide / topologia

- `ex104`–`ex124`
- `ex125`–`ex137`

Status:

- aktywne i wartościowe,
- część wyników jest nadal używana jako kontekst,
- ale nie są one już najlepszym miejscem do czytania podstaw `nbody`,
- niektóre pliki zachowują język przejściowy między starszą a nowszą formulacją.

### B. Rozszerzenia sektorowe i syntetyczne

Przykładowe grupy:

- `ex209` / `ex209b` — reżimy `V3/V2`, shape-space i asymptotyka dużego `t` dla `I_Y`
- `ex218+` — pochodne relacje `Phi0`, `Jc`, `R12`, `K`
- `ex224+` — łańcuchy syntezy, counting, flavor, CKM/PMNS, neutrina
- `ex242+` — EW / QCD / anomaly / baryogenesis / UV / phenomenology
- `ex271+` — aktualne summary / scorecard / lattice / signatures

Status:

- aktywne badania,
- dobre jako laboratorium liczbowe,
- nie traktować automatycznie jako "kanonicznego rdzenia `nbody`".

## 3. LEGACY-TRANSLATIONAL

To są aktywne pliki poza `_archiwum/`, ale o podwyższonym ryzyku dryfu
pojęciowego.

Typowe sygnały:

- język `Formulation B`
- ghost-point jako główna intuicja
- `g0 ~ 1.24` zamiast kanonicznego `g0^e ~ 0.869...`
- stare ścieżki selekcji `alpha*`, `alpha_K_old`
- narracja sprzed rozdzielenia `sin(r)/r` i EFT Yukawy

Najczęściej dotyczy to części skryptów z zakresu:

- `ex55`–`ex103`
- wybranych plików w `ex58`–`ex70`
- niektórych starszych syntez w późniejszych numerach

Szczególnie spójny blok legacy tworzą dziś `ex66`–`ex79`: stary program
`A_tail^4`, `g0* ~ 1.24`, `alpha*`, `delta(alpha)` i bifurkacje wokół dawnych
reguł selekcji.

Drugim zwartym blokiem legacy jest `ex80`–`ex103`: algebra `alpha*`, profile
`F(alpha)` / `G3(alpha)`, sektor `tau`, sum rules i mnożniki hierarchii w
starym języku selekcji.

Nie oznacza to, że te pliki są "złe". Oznacza tylko, że należy je czytać jako:

- mapę przejścia,
- eksplorację historycznych hipotez,
- źródło pomocnicze, nie źródło kanoniczne.

## 4. ARCHIVE

Wszystko w `examples/_archiwum/` traktujemy jako archiwum.

To obejmuje m.in.:

- dawny tor Efimov / bound-state / FDM / rotacje
- pierwsze wersje defect / Path B
- stare dynamiki 3-ciałowe
- wczesne hipotezy, które później zostały poprawione, zawężone albo odrzucone

Jeśli archiwum jest cytowane przez nowszy plik, to nowszy plik ma pierwszeństwo.

## 5. Szybka reguła użytkowa

Jeśli chcesz uruchomić tylko rzeczy naprawdę reprezentatywne dla obecnego
`nbody`, zacznij od:

1. `verify_nbody_eom_quick.py`
2. `verify_nbody_canonical_quick.py`
3. `verify_nbody_bridge_extended.py`
4. `verify_nbody_lyapunov_quick.py`
5. `ex195`–`ex197`
6. `ex205`

Jeśli chcesz czytać historię dojścia do obecnej formulacji, dopiero potem
wchodź w:

- `ex104`–`ex137`
- `ex55`–`ex103`
- `_archiwum/`
