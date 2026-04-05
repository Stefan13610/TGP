# Zrzuty wyjścia

Pliki `ex*_output.txt` to **zapis stdout** z wybranych uruchomień skryptów (porównania, audyty). Nie są wymagane do działania kodu; można je odtworzyć, ponownie uruchamiając odpowiedni `ex*.py` w katalogu `examples/`.

**ex114:** `ex114_tail_phase_map.csv` (z kolumnami `B_tail_altwin`, `C_tail_altwin`, …), `ex114_generations.csv`, `ex114_generations_altwin.csv`, `ex114_tail_phase_BC.png`, `ex114_tail_phase_generations.png`, `ex114_tail_phase_generations_twowin.png` — `ex114_tail_phase_map.py` (mapa (B,C), triada, porównanie okien tail [20,35] vs [22,36]).

**ex115:** `ex115_window_scan.csv`, `ex115_window_scan_admissible.png` — `ex115_tail_window_scan.py` (obszar dopuszczalnych okien tail na triadzie).

**ex116:** `ex116_rL_vs_epsilon.csv` — `ex116_tail_rL_from_epsilon.py` (minimalne \(r_L\) z progu średniej \(|g-1|\) na triadzie, \(r_R=35\)).

**ex117:** `ex117_linear_residual.csv` — `ex117_linear_operator_residual.py` (reszta operatora liniowego na \(h=g-1\)); teoria: `ex117_tail_linearization.md`.

**ex118:** `ex118_rmse_vs_zeta.csv` — `ex118_rmse_vs_linear_residual.py` (powiązanie RMSE fitu z \(\zeta\)).

**ex119:** `ex119_pipeline_digest.csv`, `ex119_pipeline_digest.png` — `ex119_tail_pipeline_digest.py` (zbiorczy raport tail).

**ex120:** `ex120_ex115_overlay.png`, `ex120_ex115_overlay_meta.csv` — `ex120_ex115_scan_overlay.py` (warstwa na CSV ex115).

**ex121:** `ex121_tail_suite_log.txt` — `ex121_tail_suite_runner.py` (pełny łańcuch tail + ex122–ex124, ~2–3 min).

**ex122:** `ex122_cross_term_ratios.csv` — `ex122_cross_term_ratio.py`.

**ex123:** `ex123_koide_epistemics.csv` — `ex123_koide_epistemics.py`.

**ex124:** `ex124_dense_solver_compare.csv`, `ex124_dense_A_compare.png` — `ex124_dense_g0_solver_compare.py`.

**ex154:** `ex154_softening_scan.csv` — `ex154_lyapunov_scan_softening_csv.py` (P1: `lambda_max` vs `softening`, TGP pairwise).

**ex151 (opcjonalnie):** `ex151_beta_scan.csv` — `ex151_lyapunov_scan_beta_pairwise.py --out ...` (P1: `lambda_max` vs `beta`).

**ex153 (opcjonalnie):** spektrum 18D — `ex153_lyapunov_spectrum_sum_newton.py --out ...` (kolumny `j`, `lambda` + wiersz `sum`).

**ex155:** `ex155_beta_softening_grid.csv` — `ex155_lyapunov_grid_beta_softening_csv.py`.

**ex166 / ex168:** `ex166_yukawa_convergence_grid.csv` — `ex166_lyapunov_yukawa_feynman_convergence_grid_csv.py`; opcjonalnie `ex168_ex166_lambda_grid.png` — `ex168_plot_ex166_convergence_csv.py` (po wygenerowaniu CSV).

**ex182:** `ex182_seven_row_overview.csv` — `ex182_lyapunov_seven_row_overview_csv.py` (P1: siedem wierszy jak `ex181`, eksport CSV).

**ex183:** opcjonalnie `ex183_ex182_lambda_bar.png` — `ex183_summarize_ex182_seven_row_csv.py` (bez `--quick`; wymaga `matplotlib`).

**ex184:** `ex184_seven_row_t_final_scan.csv` — `ex184_lyapunov_seven_row_t_final_scan_csv.py` (P1.C: kilka `t_final`, za każdym razem siedem wierszy jak `ex181`).

**ex185:** opcjonalnie `ex185_ex184_lambda_vs_t.png` — `ex185_summarize_ex184_t_scan_csv.py` (bez `--quick`; wymaga `matplotlib`).

**ex186:** `ex186_matched_h_family.csv` — `ex186_lyapunov_matched_h_family_csv.py` (P1: cztery wiersze jak `ex178`, kolumny jak `ex182`).

**ex187:** opcjonalnie `ex187_ex186_lambda_bar.png` — `ex187_summarize_ex186_matched_h_csv.py` (bez `--quick`; wymaga `matplotlib`).

**ex188:** `ex188_matched_h_t_final_scan.csv` — `ex188_lyapunov_matched_h_family_t_final_scan_csv.py` (P1.C: kilka `t_final`, za każdym razem cztery wiersze jak `ex178`).

**ex189:** opcjonalnie `ex189_ex188_lambda_vs_t.png` — `ex189_summarize_ex188_matched_h_t_scan_csv.py` (bez `--quick`; wymaga `matplotlib`).

**ex190:** brak pliku wyjściowego — `ex190_leapfrog_energy_drift_matched_h_family_short.py` (stdout: max `|ΔE/E0|` leapfroga; ta sama rodzina IC co `ex178`).

**ex191:** brak pliku wyjściowego — `ex191_rk45_energy_diag_matched_h_family_short.py` (stdout: RK45 `success`, `max|ΔE/E0|`; IC jak `ex178`).

**ex192:** brak pliku wyjściowego — `ex192_matched_h_leapfrog_vs_rk45_energy_table.py` (stdout: tabela LF vs RK przy wspólnym `t_final`).

**ex193:** `ex193_matched_h_lf_vs_rk_energy.csv` — `ex193_matched_h_lf_vs_rk45_energy_csv.py` (te same liczby co `ex192`, eksport CSV).

**ex194:** brak pliku wyjściowego — `ex194_summarize_ex193_matched_h_lf_vs_rk_csv.py` (odczyt / walidacja CSV `ex193`).

**ex158 (opcjonalnie):** `ex158_seed_lambda.csv` — `ex158_lyapunov_newton_seed_spread_leapfrog.py --out ...`.
