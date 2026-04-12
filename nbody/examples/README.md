# `nbody/examples` — skrypty numeryczne

## Czytaj najpierw

Po synchronizacji `nbody` z rdzeniem teorii nie wszystkie `exNNN` mają ten sam
status. Kanoniczna mapa:

- [STATUS_MAP.md](STATUS_MAP.md)

Najkrótsza ścieżka do bieżącej warstwy `nbody`:

1. `verify_nbody_canonical_quick.py`
2. `verify_nbody_bridge_extended.py`
3. `verify_nbody_eom_quick.py`
4. `verify_nbody_lyapunov_quick.py`
5. `ex195_soliton_ksub_full_vs_lpa.py`
6. `ex196_phi_fp_full_vs_lpa.py`
7. `ex197_optimal_g0_full_form.py`
8. `ex205_path_c_yukawa_from_defect.py`

Przy konflikcie między starszym docstringiem a aktualnym opisem teorii
pierwszeństwo mają:

1. `../THEORY_SYNC_NBODY.md`
2. `STATUS_MAP.md`
3. ten README

## Uruchamianie

Z tego katalogu:

```text
python ex106_path9_formalization.py
python ex113_oj3_tau_koide.py
```

Albo z IDE z ustawionym katalogiem roboczym na `examples/`.

## Indeks (skrót)

### Regresja „starego” toru Efimov / okien kwantowych

Status: **legacy regression / archive bridge**, nie główny tor obecnego `nbody`.

| Plik | Opis |
|------|------|
| `verify_nbody_canonical_quick.py` | Krótki runner współczesnej ścieżki kanonicznej: `ex195` → `ex196` → `ex197` → `ex205` → `verify_nbody_eom_quick.py`. |
| `verify_nbody_bridge_extended.py` | Pełniejszy, stabilny runner mostu: `ex195` → `ex196` → `ex197` → `ex205` → `verify_nbody_eom_quick.py` → reprezentatywne checki Lyapunova (`ex148`, `ex170`, `ex172`, `ex198`). |
| `verify_all.py` | Porównuje zaimplementowane w pliku wzory (Feynman, ZP, progi C, FD) z wartościami referencyjnymi z sesji ex23–ex33. Skrypty ex23–ex33 są w `_archiwum/`; sam `verify_all.py` jest samowystarczalny. |

### Tor roboczy 2026 (ścieżka 9, kosmologia, solitony, O-J3)

Status: **active-exploratory**. Ta grupa zawiera ważne wyniki i aktywne skrypty,
ale nie całość należy traktować jako kanoniczny rdzeń obecnego `nbody`.

| Plik | Temat |
|------|--------|
| `ex104_kappa_from_action.py` | κ z poprawnej akcji (LLR / ψ) |
| `ex105_ns_full_pipeline.py` | Pełny pipeline n_s |
| `ex106_path9_formalization.py` | Ścieżka 9: ODE, A_tail, φ-FP, r₂₁ |
| `ex107_ns_mukhanov_sasaki.py` | Mukhanov–Sasaki / n_s |
| `ex108_chiral_mass_split.py` | Chiralność / masa L–R |
| `ex109_u1_gauge_emergence.py` | Emergencja U(1) |
| `ex110_agamma_phi0_desi_dr2.py` | a_Γ·Φ₀ vs DESI |
| `ex111_tau_mass_soliton_energy.py` | Scenariusz S3: energia solitonu vs masa (wynik negatywny) |
| `ex112_soliton_energy_ksub.py` | OP-G / OP-E, K_sub, E_core |
| `ex113_oj3_tau_koide.py` | O-J3: Koide + A_tail, g₀^τ |
| `ex114_tail_phase_map.py` | Mapa (B,C,A) vs g₀; triada; CSV + PNG baz/alt; `*_twowin.png`; testy P1–P9 |
| `ex114_tail_okna_podsumowanie.md` | Algebra fitu, heurystyka asymptotyki, **procedura domknięcia okien** (kryteria P4–P9) |
| `ex115_tail_window_scan.py` | Skan siatki \((r_L,r_R)\): `ex115_window_scan.csv`, `ex115_window_scan_admissible.png` |
| `ex116_tail_rL_from_epsilon.py` | Dynamiczne \(r_L^*(\varepsilon)\) przy stałym \(r_R=35\): `ex116_rL_vs_epsilon.csv` |
| `ex117_tail_linearization.md` | Wyprowadzenie \(h''+2h'/r+h\approx 0\) z ODE ex106/ex114 |
| `ex117_linear_operator_residual.py` | Triada: iloraz \(\|L[h]\|_{\rm RMS}/\|h\|_{\rm RMS}\) → `ex117_linear_residual.csv` |
| `ex118_rmse_vs_linear_residual.py` | Pearson(\(\zeta\), RMSE/A) triada × 2 pasma → `ex118_rmse_vs_zeta.csv` |
| `ex119_tail_pipeline_digest.py` | **Digest** ex114–ex118: `ex119_pipeline_digest.csv` + `.png` |
| `ex120_ex115_scan_overlay.py` | Wizualizacja skanu ex115: `ex120_ex115_overlay.png` + `*_meta.csv` |
| `ex121_tail_suite_runner.py` | **Runner** ex114→ex120 + ex122–ex124; log: `ex121_tail_suite_log.txt` (`--dry-run`; ex124 ~1 min) |
| `PLAN_TAIL_DOMKNIECIE_ROZWOJ.md` | Plan domknięcia / rozwoju (fazy A–D) |
| `ex122_cross_term_ratio.py` | \(|(\alpha/g)g'^2|/|V'|\) triada → `ex122_cross_term_ratios.csv` |
| `ex122_tail_cross_term_sketch.md` | Szkic ilorazu + uwaga o p99 vs max |
| `ex123_koide_epistemics.py` | Tabela \(r_{31}^{\varphi^2}\) vs Koide przy \(r_{21}\) ze Ścieżki 9 → `ex123_koide_epistemics.csv` |
| `ex123_koide_epistemics.md` | Uzasadnienie epistemiczne (dwa domknięcia) |
| `ex124_dense_g0_solver_compare.py` | Gęsta siatka g₀ + reg vs K_sub → `ex124_dense_solver_compare.csv`, `ex124_dense_A_compare.png` |
| `ex125_qk32_phase_derivation.py` … `ex137_top4_alpha_phi_ladder.py` | Eksploracje Z₃, Brannen, top4/top5, WKB — szczegóły w docstringach; CSV/PNG w `_outputs/` |
| `ex138_eom_tgp_coulomb_leapfrog.py` | `build_tgp_integration_pair('coulomb_3b')` + leapfrog (Coulomb na \(I\)) |
| `ex139_yukawa_feynman_leapfrog.py` | `build_tgp_integration_pair('yukawa_feynman')` + leapfrog (dokładne \(I_Y\), wolniejsze) |
| `ex140_yukawa_feynman_N4_scaling_rk.py` | \(N=4\): koszt acc vs \(N=3\), leapfrog + krótki RK45 (ten sam backend) |
| `ex141_net_force_translational_invariance.py` | \(\sum_i C_i a_i \approx 0\) — wszystkie `TGP_INTEGRATION_BACKENDS` |
| `ex142_yukawa_overlap_quadrature_convergence.py` | Zbieżność `n_quad` w `yukawa_overlap_exact` (koszt vs błąd) |
| `ex143_net_torque_so3_invariance.py` | \(\sum \mathbf{x}_i\times\mathbf{F}_i \approx 0\) — inwariancja \(SO(3)\) |
| `ex144_conserved_P_L_leapfrog.py` | Zachowanie \(\mathbf{P}\), \(\mathbf{L}\) na krótkim leapfrog (`coulomb_3b`) |
| `ex145_minus_grad_V_matches_force.py` | \(-\nabla V\) (FD) vs siły z `acc_fn` (`coulomb_3b`) |
| `ex146_yukawa_feynman_grad_matches_force.py` | Jak `ex145`, backend `yukawa_feynman`, \(N=3\) |
| `ex147_com_acceleration_zero.py` | \(\ddot{\mathbf{R}}_{\rm cm}\equiv \mathbf{a}_{\rm cm}\approx 0\) (wewnętrzne \(V\)) |
| `ex148_lyapunov_tgp_vs_newton_pairwise.py` | P1: \(\lambda_{\max}\) (Benettin) — Newton vs TGP pairwise, problem pitagorejski; `--quick` |
| `ex149_lyapunov_spectrum_newton.py` | P1: 6 wykładników (Newton); domyślnie leapfrog; ``--rk4`` |
| `ex150_lyapunov_pairwise_vs_coulomb3b.py` | P1: \(\lambda_{\max}\) — `pairwise` vs `coulomb_3b` (FD Jacobian) |
| `ex151_lyapunov_scan_beta_pairwise.py` | P1: skan `beta` (`gamma=beta`); opcjonalnie `--out` CSV |
| `ex152_lyapunov_yukawa_feynman_short.py` | P1: `yukawa_feynman` vs `coulomb_3b` (RK4+Benettin); po biegu: stdout RK45 `max|dE/E0|`; `--n-quad`, `--quick` |
| `ex153_lyapunov_spectrum_sum_newton.py` | P1: spektrum 18D (Newton); domyślnie leapfrog (`Σλ≈0`); ``--rk4`` |
| `ex154_lyapunov_scan_softening_csv.py` | P1: skan `softening` → `_outputs/ex154_softening_scan.csv` |
| `ex155_lyapunov_grid_beta_softening_csv.py` | P1: siatka `beta` × `softening` → `_outputs/ex155_beta_softening_grid.csv` |
| `ex156_lyapunov_newton_rk4_vs_leapfrog.py` | P1: `λ_max` Newton — RK4 vs leapfrog (styczna) |
| `ex157_lyapunov_tgp_pairwise_rk4_vs_leapfrog.py` | P1: `λ_max` TGP pairwise — RK4 vs leapfrog (FD `J`) |
| `ex158_lyapunov_newton_seed_spread_leapfrog.py` | P1: rozrzut `λ_max` vs seed styczny (Newton, leapfrog); `--out` |
| `ex159_lyapunov_coulomb3b_rk4_vs_leapfrog.py` | P1: `coulomb_3b`, FD `J`, RK4 vs leapfrog (kosztowne) |
| `ex160_lyapunov_yukawa_feynman_leapfrog_only.py` | P1: `yukawa_feynman`, leapfrog+tangent, FD `J` (krótko) |
| `ex161_lyapunov_pairwise_beta_d_rep_scan.py` | P1: skan β z kolumną `d_rep` (2-ciało, `C_ref=min C`), jawny `J` dla `V_2` |
| `ex162_lyapunov_yukawa_feynman_fd_vs_split_jacobian.py` | P1: `yukawa_feynman` — `λ_max` przy pełnym FD `J` vs split (`V_2` analitycznie) |
| `ex163_lyapunov_yukawa_feynman_split_vs_analytic_jacobian.py` | P1: split (FD `V_3`) vs pełny analityczny `J` (`three_body_force_jacobian_exact`) |
| `ex164_lyapunov_yukawa_feynman_long_window_analytic.py` | P1: dłuższe `t` + wyższe `n_quad`, leapfrog+analityczny `J`; stdout RK45 `max|dE/E0|` (referencja, nie orbita LF) |
| `ex165_lyapunov_yukawa_feynman_refine_dt_nquad.py` | P1: ten sam `t_final`+seed, coarse vs fine `dt`/`n_quad`, analityczny `J` (zbieżność względem siatki) |
| `ex166_lyapunov_yukawa_feynman_convergence_grid_csv.py` | P1: siatka `dt`×`n_quad` → CSV `_outputs/ex166_yukawa_convergence_grid.csv` (analityczny `J`) |
| `ex167_leapfrog_energy_drift_yukawa_vs_coulomb3b_short.py` | P1: leapfrog — krotki horyzont, max `|ΔE/E0|` dla `yukawa_feynman` vs `coulomb_3b` (Burrau) |
| `ex168_plot_ex166_convergence_csv.py` | P1: odczyt CSV `ex166`, statystyki + opcjonalny PNG `_outputs/ex168_ex166_lambda_grid.png` |
| `ex169_lyapunov_yukawa_feynman_leapfrog_fd_vs_analytic_jacobian.py` | P1: leapfrog Benettin — pelny FD `J` vs analityczny `J`, ten sam seed |
| `ex170_lyapunov_yukawa_feynman_spectrum_sum_leapfrog.py` | P1: spektrum `k=6N`, Yukawa Feynman, `|sum(lambda)|` ~ 0 (LF + analityczny `J`) |
| `ex171_lyapunov_coulomb3b_spectrum_sum_leapfrog.py` | P1: spektrum `k=6N`, `coulomb_3b`, `|sum(lambda)|` ~ 0 (LF + FD `J`) |
| `ex172_lyapunov_pairwise_spectrum_sum_leapfrog.py` | P1: spektrum `k=6N`, `pairwise` (`V_2`), `|sum(lambda)|` ~ 0 (LF + analityczny `J`) |
| `ex173_lyapunov_pairwise_leapfrog_fd_vs_analytic_spectrum.py` | P1: `pairwise`, spektrum `k=6N`, FD `J` vs analityczny `J`, ten sam seed |
| `ex174_lyapunov_newton_vs_yukawa_feynman_matched_energy.py` | P1: Burrau `x`; Newton (`v=0`) vs `yukawa_feynman` z dopasowaniem `H`; leapfrog+analityczne `J` |
| `ex175_lyapunov_burrau_newton_yukawa_three_hamiltonians.py` | P1: to samo `x` — `H_N`, Yukawa `v=0`, Yukawa `H_Y=H_N`; trzy `lambda_max` (kontekst fazowy) |
| `ex176_lyapunov_newton_vs_coulomb3b_matched_energy.py` | P1: jak `ex174`, lecz `coulomb_3b` + dopasowanie `H`; Newton `J` analityczny, Coulomb FD `J` |
| `ex177_lyapunov_newton_vs_pairwise_matched_energy.py` | P1: jak `ex148` (`V_2` + dopasowanie `H`), leapfrog + analityczne `J`; siatka jak `ex174` |
| `ex178_lyapunov_matched_energy_family_table.py` | P1: jedna tabela Newton + Yukawa + `coulomb_3b` + pairwise przy `H=H_N` |
| `ex179_lyapunov_burrau_newton_coulomb_three_hamiltonians.py` | P1: jak `ex175`, lecz `coulomb_3b` + FD `J` (trzy `H` na tym samym `x`) |
| `ex180_lyapunov_burrau_newton_pairwise_three_hamiltonians.py` | P1: jak `ex175`, lecz pairwise `V_2` (trzy `H`; RNG jak `ex177`) |
| `ex181_lyapunov_matched_energy_seven_row_overview.py` | P1: Newton + 6 wierszy TGP (`v=0` i `H=H_N` × 3 backendy) |
| `ex182_lyapunov_seven_row_overview_csv.py` | P1: jak `ex181` → CSV `_outputs/ex182_seven_row_overview.csv` |
| `ex183_summarize_ex182_seven_row_csv.py` | P1: odczyt / walidacja CSV `ex182`, opcjonalny PNG |
| `ex184_lyapunov_seven_row_t_final_scan_csv.py` | P1.C: skan `t_final` × siedem wierszy → `_outputs/ex184_seven_row_t_final_scan.csv` |
| `ex185_summarize_ex184_t_scan_csv.py` | P1.C: walidacja / podsumowanie CSV `ex184`, opcjonalny PNG |
| `ex186_lyapunov_matched_h_family_csv.py` | P1: jak `ex178` → CSV `_outputs/ex186_matched_h_family.csv` (układ kolumn jak `ex182`) |
| `ex187_summarize_ex186_matched_h_csv.py` | P1: odczyt / walidacja CSV `ex186`, opcjonalny PNG |
| `ex188_lyapunov_matched_h_family_t_final_scan_csv.py` | P1.C: skan `t_final` × cztery wiersze matched `H` → `_outputs/ex188_matched_h_t_final_scan.csv` |
| `ex189_summarize_ex188_matched_h_t_scan_csv.py` | P1.C: walidacja / podsumowanie CSV `ex188`, opcjonalny PNG |
| `ex190_leapfrog_energy_drift_matched_h_family_short.py` | P1.C: krótki leapfrog `|ΔE/E0|` dla rodziny matched-`H` jak `ex178` |
| `ex191_rk45_energy_diag_matched_h_family_short.py` | P1.C: krótki RK45 (DOP853) — diagnostyka energii matched-`H` (jak `ex178`) |
| `ex192_matched_h_leapfrog_vs_rk45_energy_table.py` | P1.C: tabela leapfrog vs RK45, wspólny `t_final` (`ex190`+`ex191`) |
| `ex193_matched_h_lf_vs_rk45_energy_csv.py` | P1.C: jak `ex192` → `_outputs/ex193_matched_h_lf_vs_rk_energy.csv` |
| `ex194_summarize_ex193_matched_h_lf_vs_rk_csv.py` | P1.C: odczyt / walidacja CSV `ex193` |
| `verify_nbody_lyapunov_quick.py` | P1: `ex148`–`ex194` z `--quick` |
| `verify_nbody_eom_quick.py` | Kolejno: `ex138`, `ex141`–`ex147` (bez ciężkich `ex139`/`ex140`, bez ex148+) |

Nowsze diagnostyki geometrii / asymptotyki `I_Y`:

| Plik | Temat |
|------|--------|
| `ex209_v3_v2_regime_msp_scan.py` | Reżimy `V3/V2` vs `m_sp d`, geometria i skala `C` |
| `ex209b_lambda_shape_space_map.py` | Mapa shape-space dla wykładnika tłumienia `lambda(q1,q2)` sterującego asymptotyką dużego `t` dokładnego `I_Y`; zapisuje `_outputs/ex209b_lambda_shape_space.csv` |

### Ścieżka 8–9 — eksploracja i pośrednie weryfikacje (ex55–ex103)

Status: **legacy-translational / exploratory**. Używać ostrożnie przy pracy z
aktualną formulacją; część skryptów zachowuje starszy język (`ghost`, dawne
`g0`, przejściowe hipotezy selekcji).

Skrypty `ex55`–`ex103` to łańcuch numeryczny wokół **A_tail**, **α\***, okien, stabilności i wzorów sumacyjnych. Część wcześniejszych hipotez została zamknięta negatywnie lub zastąpiona przez ex106+ (patrz archiwum planów w `_archiwum/` repozytorium). Przy konflikcie **nadrzędne są ex106 / ex112 / ex113** oraz aktualny opis w **`ANALIZA_KRYTYCZNA_v6.md`** (korzeń `TGP_v1/`, poza folderem `nbody`).

W praktyce szczególnie `ex66`–`ex79` należy dziś czytać jako spójny blok
`legacy-translational`: masa z `A_tail^4`, `g0* ~ 1.24`, `alpha*`, bifurkacje i
powiązane reguły selekcji.

Podobnie `ex80`–`ex103` to dalszy ciąg tego samego starego toru:
algebraiczne reguły dla `alpha*`, profile `F(α)` / `G3(α)`, sektor `tau`,
sum rules i mnożniki hierarchii. To nadal materiał przejściowy, a nie
współczesna ścieżka kanoniczna.

Ciąg dalszy (2026-04): **ex125–ex137** — eksploracja Z₃ / Brannen / top4–top5 / WKB (patrz nagłówki plików i `_outputs/` po uruchomieniu).

### Podkatalogi

| Katalog | Zawartość |
|---------|-----------|
| `_archiwum/` | Starsze eksperymenty (ex1–ex54, ex związane z Efimov, K13–K18, GW, FDM, …). Status: **ARCHIVE**. Zob. [_archiwum/README.md](_archiwum/README.md). |
| `_outputs/` | Zrzuty stdout z uruchomień (`ex88_output.txt`, …). Zob. [_outputs/README.md](_outputs/README.md). |

## Konwencja nazw

- `exNNN_opis.py` — aktywny eksperyment.
- `exNNNv2_*.py` — wariant tego samego numeru (np. szybsza lub alternatywna implementacja).
