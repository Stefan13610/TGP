#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_nbody_lyapunov_quick.py
==============================
Regresja P1 (Lyapunov / chaos): ``ex148``--``ex194`` z ``--quick``.

Nie wlacza sie do ``verify_nbody_eom_quick.py`` (osobny koszt czasu).

Uruchomienie z korzenia ``TGP_v1/``:

    python nbody/examples/verify_nbody_lyapunov_quick.py
"""

from __future__ import annotations

import os
import subprocess
import sys

_SCRIPTS = (
    "ex148_lyapunov_tgp_vs_newton_pairwise.py",
    "ex149_lyapunov_spectrum_newton.py",
    "ex150_lyapunov_pairwise_vs_coulomb3b.py",
    "ex151_lyapunov_scan_beta_pairwise.py",
    "ex152_lyapunov_yukawa_feynman_short.py",
    "ex153_lyapunov_spectrum_sum_newton.py",
    "ex154_lyapunov_scan_softening_csv.py",
    "ex155_lyapunov_grid_beta_softening_csv.py",
    "ex156_lyapunov_newton_rk4_vs_leapfrog.py",
    "ex157_lyapunov_tgp_pairwise_rk4_vs_leapfrog.py",
    "ex158_lyapunov_newton_seed_spread_leapfrog.py",
    "ex159_lyapunov_coulomb3b_rk4_vs_leapfrog.py",
    "ex160_lyapunov_yukawa_feynman_leapfrog_only.py",
    "ex161_lyapunov_pairwise_beta_d_rep_scan.py",
    "ex162_lyapunov_yukawa_feynman_fd_vs_split_jacobian.py",
    "ex163_lyapunov_yukawa_feynman_split_vs_analytic_jacobian.py",
    "ex164_lyapunov_yukawa_feynman_long_window_analytic.py",
    "ex165_lyapunov_yukawa_feynman_refine_dt_nquad.py",
    "ex166_lyapunov_yukawa_feynman_convergence_grid_csv.py",
    "ex167_leapfrog_energy_drift_yukawa_vs_coulomb3b_short.py",
    "ex168_plot_ex166_convergence_csv.py",
    "ex169_lyapunov_yukawa_feynman_leapfrog_fd_vs_analytic_jacobian.py",
    "ex170_lyapunov_yukawa_feynman_spectrum_sum_leapfrog.py",
    "ex171_lyapunov_coulomb3b_spectrum_sum_leapfrog.py",
    "ex172_lyapunov_pairwise_spectrum_sum_leapfrog.py",
    "ex173_lyapunov_pairwise_leapfrog_fd_vs_analytic_spectrum.py",
    "ex174_lyapunov_newton_vs_yukawa_feynman_matched_energy.py",
    "ex175_lyapunov_burrau_newton_yukawa_three_hamiltonians.py",
    "ex176_lyapunov_newton_vs_coulomb3b_matched_energy.py",
    "ex177_lyapunov_newton_vs_pairwise_matched_energy.py",
    "ex178_lyapunov_matched_energy_family_table.py",
    "ex179_lyapunov_burrau_newton_coulomb_three_hamiltonians.py",
    "ex180_lyapunov_burrau_newton_pairwise_three_hamiltonians.py",
    "ex181_lyapunov_matched_energy_seven_row_overview.py",
    "ex182_lyapunov_seven_row_overview_csv.py",
    "ex183_summarize_ex182_seven_row_csv.py",
    "ex184_lyapunov_seven_row_t_final_scan_csv.py",
    "ex185_summarize_ex184_t_scan_csv.py",
    "ex186_lyapunov_matched_h_family_csv.py",
    "ex187_summarize_ex186_matched_h_csv.py",
    "ex188_lyapunov_matched_h_family_t_final_scan_csv.py",
    "ex189_summarize_ex188_matched_h_t_scan_csv.py",
    "ex190_leapfrog_energy_drift_matched_h_family_short.py",
    "ex191_rk45_energy_diag_matched_h_family_short.py",
    "ex192_matched_h_leapfrog_vs_rk45_energy_table.py",
    "ex193_matched_h_lf_vs_rk45_energy_csv.py",
    "ex194_summarize_ex193_matched_h_lf_vs_rk_csv.py",
)


def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.abspath(os.path.join(here, "..", ".."))
    py = sys.executable
    failed = []
    for name in _SCRIPTS:
        path = os.path.join(here, name)
        print(f"==> {name} --quick", flush=True)
        r = subprocess.run(
            [py, path, "--quick"],
            cwd=repo,
            env={**os.environ, "PYTHONUNBUFFERED": "1"},
        )
        if r.returncode != 0:
            failed.append(name)
    if failed:
        print("FAILED:", ", ".join(failed))
        raise SystemExit(1)
    print("verify_nbody_lyapunov_quick: all PASS")


if __name__ == "__main__":
    main()
