#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex114_tail_phase_map.py
=======================
Mapa fazowa asymptotyki ogona: (B_tail, C_tail, A_tail) vs g₀.

Wzorowane na ex106_path9_formalization.py (ten sam ODE, odbicia przy g*,
ten sam ansatz: (g-1)·r ≈ B·cos(r) + C·sin(r)).

Wyjścia:
  - examples/_outputs/ex114_tail_phase_map.csv
  - examples/_outputs/ex114_generations.csv  (triada, okno bazowe [20,35])
  - examples/_outputs/ex114_generations_altwin.csv  (ta sama triada, okno [22,36])
  - examples/_outputs/ex114_tail_phase_BC.png
  - examples/_outputs/ex114_tail_phase_generations.png  (tylko okno bazowe)
  - examples/_outputs/ex114_tail_phase_generations_twowin.png  (B–C: baz vs alt, obok siebie)

Testy:
  P1: CSV zapisany
  P2: A_tail > 0 na całej siatce
  P3: przy g₀* powtórny integrał daje te same A,B,C co wiersz siatki
  P4–P5: dwa okna tail → |ΔA|/A < 8%, |Δφ| < 0.35 rad
  P6: RMSE/A < 10% (jakość fitu cos/sin w [r_L, r_R])
  P7: CSV triady + wykres generations (gdy jest g₀* i matplotlib)
  P8: triada: |Δφ| między oknami < 0.35 rad (stabilność klasy fazowej)
  P9: triada: średnie |g-1| w oknie bazowym < próg (strefa asymptotyki)

Szczegóły algebraiczne / procedura okien: ex114_tail_okna_podsumowanie.md
"""

from __future__ import annotations

import io
import sys
import warnings
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

warnings.filterwarnings("ignore")

# --- parametry zsynchronizowane z ex106 ---------------------------------
ALPHA = 2.0
G_GHOST = np.exp(-1.0 / (2.0 * ALPHA))
PHI = (1.0 + np.sqrt(5.0)) / 2.0
R21_PDG = 206.768

R_MAX = 40.0
R_START = 1e-4
R_TAIL_L = 20.0
R_TAIL_R = 35.0
MAX_STEP = 0.02
RTOL = 1e-10
ATOL = 1e-13
G_BOUNCE = G_GHOST + 0.005

# Drugie okno (przesunięte) — test odporności
R_TAIL_L_ALT = 22.0
R_TAIL_R_ALT = 36.0

OUT_DIR = Path(__file__).resolve().parent / "_outputs"

RESULTS: list[tuple[str, str, str]] = []


def check(cond: bool, label: str, detail: str = "") -> bool:
    status = "PASS" if cond else "FAIL"
    RESULTS.append((label, status, detail))
    icon = "[PASS]" if cond else "[FAIL]"
    line = f"  {icon} {label}"
    if detail:
        line += f"\n         => {detail}"
    print(line)
    return cond


def Vprime(g):
    return g**2 * (1.0 - g)


def f_kin(g):
    return 1.0 + 2.0 * ALPHA * np.log(np.maximum(g, 1e-30))


def rhs_regularized(r, y):
    g, gp = y
    g = max(g, G_BOUNCE + 1e-7)
    fg = f_kin(g)
    if abs(fg) < 1e-10:
        return [gp, 0.0]
    driving = Vprime(g)
    cross = (ALPHA / g) * gp**2
    if r < 1e-10:
        return [gp, (driving - cross) / (3.0 * fg)]
    damp = fg * 2.0 * gp / r
    return [gp, (driving - cross - damp) / fg]


def event_hit_ghost(r, y):
    return y[0] - G_BOUNCE


event_hit_ghost.terminal = True
event_hit_ghost.direction = -1


def integrate_soliton(g0: float, r_max: float | None = None, max_bounces: int = 8):
    if r_max is None:
        r_max = max(R_MAX, 15.0 * g0)
    r0 = R_START
    y0 = [g0, 0.0]
    segs_r, segs_g, segs_gp = [], [], []

    for bounce_num in range(max_bounces + 1):
        sol = solve_ivp(
            rhs_regularized,
            [r0, r_max],
            y0,
            method="DOP853",
            max_step=MAX_STEP,
            rtol=RTOL,
            atol=ATOL,
            events=[event_hit_ghost],
            dense_output=False,
        )
        segs_r.append(sol.t)
        segs_g.append(sol.y[0])
        segs_gp.append(sol.y[1])

        if sol.t_events[0].size > 0 and bounce_num < max_bounces:
            r_b = float(sol.t_events[0][0])
            gp_b = float(sol.y_events[0][0, 1])
            r0 = r_b + 1e-6
            y0 = [G_BOUNCE + 1e-5, -gp_b]
        else:
            break

    r = np.concatenate(segs_r)
    g = np.concatenate(segs_g)
    gp = np.concatenate(segs_gp)
    idx = np.argsort(r)
    return r[idx], g[idx], gp[idx]


def fit_tail_detailed(
    r_arr: np.ndarray,
    g_arr: np.ndarray,
    r_L: float = R_TAIL_L,
    r_R: float = R_TAIL_R,
):
    """
    Zwraca (A, B, C, rmse) dla y = (g-1)*r ≈ B cos r + C sin r.
    rmse = sqrt(mean((y - y_hat)^2)).
    """
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    if np.sum(mask) < 10:
        return 0.0, 0.0, 0.0, float("nan")
    r_fit = r_arr[mask]
    delta_fit = (g_arr[mask] - 1.0) * r_fit
    cos_r = np.cos(r_fit)
    sin_r = np.sin(r_fit)
    X = np.column_stack([cos_r, sin_r])
    coefs, _, _, _ = np.linalg.lstsq(X, delta_fit, rcond=None)
    B, C = float(coefs[0]), float(coefs[1])
    y_hat = B * cos_r + C * sin_r
    rmse = float(np.sqrt(np.mean((delta_fit - y_hat) ** 2)))
    A = float(np.sqrt(B**2 + C**2))
    return A, B, C, rmse


def atail_bc(g0: float, r_L: float = R_TAIL_L, r_R: float = R_TAIL_R):
    r, g, _gp = integrate_soliton(g0)
    return fit_tail_detailed(r, g, r_L, r_R)


def ratio_func(g0: float) -> float:
    A1, *_ = atail_bc(g0)[:4]
    A2, *_ = atail_bc(PHI * g0)[:4]
    if A1 < 1e-10:
        return 1e10
    return (A2 / A1) ** 4 - R21_PDG


def find_g0_star() -> float | None:
    try:
        v_lo = ratio_func(1.15)
        v_hi = ratio_func(1.35)
        if v_lo * v_hi < 0:
            return float(brentq(ratio_func, 1.15, 1.35, xtol=1e-5))
        g0s = np.linspace(1.05, 1.45, 60)
        ratios = np.array([ratio_func(g) for g in g0s])
        sc = np.where(np.diff(np.sign(ratios)))[0]
        if len(sc) > 0:
            i = int(sc[0])
            return float(brentq(ratio_func, g0s[i], g0s[i + 1], xtol=1e-5))
    except Exception:
        pass
    return None


def wrap_angle_diff(d: float) -> float:
    """Różnica kątów w (-pi, pi]."""
    d = (d + np.pi) % (2 * np.pi) - np.pi
    return float(d)


def dist_B_eq_C(B: float, C: float) -> float:
    """Odległość punktu (B,C) od prostej B=C w metryce euklidesowej."""
    return abs(B - C) / np.sqrt(2.0)


def tail_window_mean_abs_g_minus_one(
    g0: float,
    r_L: float = R_TAIL_L,
    r_R: float = R_TAIL_R,
) -> float:
    """Średnia |g(r)-1| na [r_L, r_R] — diagnostyka „czy okno jest w strefie g≈1”."""
    r, g, _ = integrate_soliton(float(g0))
    mask = (r >= r_L) & (r <= r_R)
    if np.sum(mask) < 5:
        return float("nan")
    return float(np.mean(np.abs(g[mask] - 1.0)))


def generation_table(
    g0_star: float,
    r_L: float = R_TAIL_L,
    r_R: float = R_TAIL_R,
) -> list[dict]:
    """Metryki (B,C) dla triady jak w ex106: g₀^e, φ·g₀^e, φ²·g₀^e."""
    out = []
    triple = [
        ("e", g0_star),
        ("mu", PHI * g0_star),
        ("tau", PHI**2 * g0_star),
    ]
    for name, g0 in triple:
        A, B, C, rmse = atail_bc(g0, r_L, r_R)
        phi = float(np.arctan2(C, B))
        d_bc = dist_B_eq_C(B, C)
        out.append(
            {
                "gen": name,
                "g0": float(g0),
                "A_tail": A,
                "B_tail": B,
                "C_tail": C,
                "phase_rad": phi,
                "phase_deg": float(np.degrees(phi)),
                "dist_B_eq_C": float(d_bc),
                "dist_B_eq_C_over_A": float(d_bc / max(A, 1e-15)),
                "abs_C_over_A": float(abs(C) / max(A, 1e-15)),
                "abs_B_over_A": float(abs(B) / max(A, 1e-15)),
                "delta_from_pi4_deg": float(np.degrees(wrap_angle_diff(phi - np.pi / 4))),
                "rmse_tail": rmse,
            }
        )
    return out


def main():
    print("=" * 70)
    print("EX114: MAPA (B_tail, C_tail, A_tail) vs g₀  [ex106-equivalent ODE]")
    print("=" * 70)
    print(f"  Okno bazowe tail: r ∈ [{R_TAIL_L}, {R_TAIL_R}]")
    print(f"  Okno alternatywne: r ∈ [{R_TAIL_L_ALT}, {R_TAIL_R_ALT}]")
    print()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    g0_star = find_g0_star()
    print(f"  g₀* (φ-FP) = {g0_star:.6f}" if g0_star else "  g₀* (φ-FP): nie znaleziono")
    print()

    # Siatka: gęściej przy małym g₀, obejmuje e, μ, τ z ex106
    g0_lo, g0_hi = 1.10, 3.40
    n_main = 52
    g0_grid = np.linspace(g0_lo, g0_hi, n_main)
    if g0_star is not None:
        g0_grid = np.unique(np.sort(np.append(g0_grid, g0_star)))

    rows = []
    max_rel_A = 0.0
    max_dphi = 0.0
    max_rmse_ratio = 0.0

    for g0 in g0_grid:
        r, g, _ = integrate_soliton(float(g0))
        A, B, C, rmse = fit_tail_detailed(r, g, R_TAIL_L, R_TAIL_R)
        A2, B2, C2, rmse2 = fit_tail_detailed(r, g, R_TAIL_L_ALT, R_TAIL_R_ALT)
        phi = float(np.arctan2(C, B))
        phi2 = float(np.arctan2(C2, B2))
        dphi = abs(wrap_angle_diff(phi2 - phi))
        rel_A = abs(A - A2) / max(A, 1e-15)
        max_rel_A = max(max_rel_A, rel_A)
        max_dphi = max(max_dphi, dphi)
        rmse_ratio = rmse / max(A, 1e-15)
        max_rmse_ratio = max(max_rmse_ratio, rmse_ratio)

        rows.append(
            {
                "g0": float(g0),
                "A_tail": A,
                "B_tail": B,
                "C_tail": C,
                "phase_rad": phi,
                "phase_deg": np.degrees(phi),
                "rmse_tail": rmse,
                "A_tail_altwin": A2,
                "B_tail_altwin": B2,
                "C_tail_altwin": C2,
                "phase_rad_altwin": phi2,
                "phase_deg_altwin": np.degrees(phi2),
                "dphi_altwin_rad": dphi,
                "rel_delta_A_altwin": rel_A,
            }
        )

    csv_path = OUT_DIR / "ex114_tail_phase_map.csv"
    header = list(rows[0].keys())
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in header) + "\n")
    print(f"  Zapisano: {csv_path}")
    print()

    gen_rows: list[dict] = []
    gen_rows_alt: list[dict] = []
    gen_csv_path = OUT_DIR / "ex114_generations.csv"
    gen_alt_csv_path = OUT_DIR / "ex114_generations_altwin.csv"
    max_dphi_triad = 0.0
    max_mean_abs_gm1_triad = 0.0
    if g0_star is not None:
        print("  --- Średnie |g-1| w oknie bazowym (triada) ---")
        triad_pairs = [
            ("e", g0_star),
            ("mu", PHI * g0_star),
            ("tau", PHI**2 * g0_star),
        ]
        for name, g0t in triad_pairs:
            mdev = tail_window_mean_abs_g_minus_one(g0t, R_TAIL_L, R_TAIL_R)
            if not np.isnan(mdev):
                max_mean_abs_gm1_triad = max(max_mean_abs_gm1_triad, mdev)
            print(f"    {name:3s}  g₀={g0t:.5f}  śr.|g-1| = {mdev:.5f}")
        print()

        gen_rows = generation_table(g0_star, R_TAIL_L, R_TAIL_R)
        gen_rows_alt = generation_table(g0_star, R_TAIL_L_ALT, R_TAIL_R_ALT)
        gh = list(gen_rows[0].keys())
        with open(gen_csv_path, "w", encoding="utf-8") as f:
            f.write(",".join(gh) + "\n")
            for row in gen_rows:
                f.write(",".join(str(row[k]) for k in gh) + "\n")
        with open(gen_alt_csv_path, "w", encoding="utf-8") as f:
            f.write(",".join(gh) + "\n")
            for row in gen_rows_alt:
                f.write(",".join(str(row[k]) for k in gh) + "\n")
        print(f"  Zapisano: {gen_csv_path}")
        print(f"  Zapisano: {gen_alt_csv_path}")
        print()
        print(f"  --- Triada: okno r∈[{R_TAIL_L},{R_TAIL_R}] vs r∈[{R_TAIL_L_ALT},{R_TAIL_R_ALT}] ---")
        for row, rowa in zip(gen_rows, gen_rows_alt):
            dphi_t = abs(wrap_angle_diff(row["phase_rad"] - rowa["phase_rad"]))
            max_dphi_triad = max(max_dphi_triad, dphi_t)
            print(
                f"    {row['gen']:3s}  g₀={row['g0']:.5f}  "
                f"φ_base={row['phase_deg']:7.2f}°  φ_alt={rowa['phase_deg']:7.2f}°  "
                f"|Δφ|={np.degrees(dphi_t):5.2f}°"
            )
            print(
                f"         B_base={row['B_tail']:+.5f} C_base={row['C_tail']:+.5f}  "
                f"B_alt={rowa['B_tail']:+.5f} C_alt={rowa['C_tail']:+.5f}"
            )
        print()
        print("  --- Metryki (tylko okno bazowe) ---")
        for row in gen_rows:
            print(
                f"    {row['gen']:3s}  A={row['A_tail']:.5f}  "
                f"|B−C|/(√2·A)={row['dist_B_eq_C_over_A']:.4f}  "
                f"Δ(φ−π/4)={row['delta_from_pi4_deg']:+.2f}°"
            )
        print()

    # --- wykres ---
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        g0s = np.array([r["g0"] for r in rows])
        As = np.array([r["A_tail"] for r in rows])
        Bs = np.array([r["B_tail"] for r in rows])
        Cs = np.array([r["C_tail"] for r in rows])
        Bs_alt = np.array([r["B_tail_altwin"] for r in rows])
        Cs_alt = np.array([r["C_tail_altwin"] for r in rows])
        phs = np.array([r["phase_rad"] for r in rows])

        fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.5))
        ax = axes[0, 0]
        ax.plot(g0s, As, "b-", lw=1.2)
        ax.set_xlabel("g₀")
        ax.set_ylabel("A_tail")
        ax.set_title("Amplituda ogona √(B²+C²)")
        ax.grid(True, alpha=0.3)
        if g0_star:
            ax.axvline(g0_star, color="gray", ls="--", lw=0.8, label="g₀*")
            ax.axvline(PHI * g0_star, color="orange", ls=":", lw=0.8)
            ax.axvline(PHI**2 * g0_star, color="green", ls=":", lw=0.8)
            ax.legend(loc="best", fontsize=7)

        ax = axes[0, 1]
        ax.plot(g0s, Bs, label="B", lw=1.0)
        ax.plot(g0s, Cs, label="C", lw=1.0)
        ax.set_xlabel("g₀")
        ax.set_ylabel("współczynnik")
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title("B_tail, C_tail (fit cos/sin)")

        ax = axes[1, 0]
        ax.plot(g0s, np.degrees(phs), "purple", lw=1.2)
        ax.set_xlabel("g₀")
        ax.set_ylabel("φ_faza = atan2(C,B) [°]")
        ax.grid(True, alpha=0.3)
        ax.set_title("Faza ogona")

        ax = axes[1, 1]
        sc = ax.scatter(Bs, Cs, c=g0s, cmap="viridis", s=12, alpha=0.85)
        plt.colorbar(sc, ax=ax, label="g₀")
        ax.set_xlabel("B_tail")
        ax.set_ylabel("C_tail")
        ax.set_title("Parametrycznie (B,C), kolor = g₀")
        ax.grid(True, alpha=0.3)
        ax.set_aspect("equal", adjustable="box")

        plt.tight_layout()
        png_path = OUT_DIR / "ex114_tail_phase_BC.png"
        fig.savefig(png_path, dpi=140)
        plt.close(fig)
        print(f"  Wykres: {png_path}")

        # --- wykres: triada + proste B=C, B=0, C=0 ---
        if g0_star is not None and gen_rows:
            fig2, axg = plt.subplots(figsize=(7.2, 7.0))
            sc2 = axg.scatter(Bs, Cs, c=g0s, cmap="viridis", s=10, alpha=0.35, label="siatka g₀")
            plt.colorbar(sc2, ax=axg, label="g₀", fraction=0.046)

            b_all = list(Bs) + [r["B_tail"] for r in gen_rows]
            c_all = list(Cs) + [r["C_tail"] for r in gen_rows]
            pad = 0.08 * (max(max(abs(x) for x in b_all), max(abs(x) for x in c_all), 0.2))
            lim_lo = min(min(b_all), min(c_all)) - pad
            lim_hi = max(max(b_all), max(c_all)) + pad
            t = np.linspace(lim_lo, lim_hi, 120)

            axg.plot(t, t, "k--", lw=1.0, alpha=0.7, label="B = C")
            axg.axhline(0.0, color="gray", ls=":", lw=0.9, alpha=0.8, label="C = 0")
            axg.axvline(0.0, color="gray", ls="-.", lw=0.9, alpha=0.8, label="B = 0")

            colors = {"e": "tab:blue", "mu": "tab:orange", "tau": "tab:green"}
            for row in gen_rows:
                axg.scatter(
                    [row["B_tail"]],
                    [row["C_tail"]],
                    s=120,
                    marker="*",
                    c=colors.get(row["gen"], "red"),
                    edgecolors="black",
                    linewidths=0.6,
                    zorder=5,
                    label=f"{row['gen']} (g₀={row['g0']:.3f})",
                )
            axg.set_xlabel("B_tail")
            axg.set_ylabel("C_tail")
            axg.set_title("Klasy fazowe: triada vs proste referencyjne")
            axg.set_aspect("equal", adjustable="box")
            axg.grid(True, alpha=0.3)
            axg.legend(loc="upper left", fontsize=7, framealpha=0.9)
            axg.set_xlim(lim_lo, lim_hi)
            axg.set_ylim(lim_lo, lim_hi)

            png2 = OUT_DIR / "ex114_tail_phase_generations.png"
            fig2.savefig(png2, dpi=140, bbox_inches="tight")
            plt.close(fig2)
            print(f"  Wykres: {png2}")

            # --- te same proste + triada: okno baz vs okno alt (wspólna skala) ---
            if gen_rows_alt:
                import matplotlib.colors as mcolors

                fig3, (axL, axR) = plt.subplots(1, 2, figsize=(13.0, 6.4))
                vmin, vmax = float(np.min(g0s)), float(np.max(g0s))
                norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

                b_all = (
                    list(Bs)
                    + list(Bs_alt)
                    + [r["B_tail"] for r in gen_rows]
                    + [r["B_tail"] for r in gen_rows_alt]
                )
                c_all = (
                    list(Cs)
                    + list(Cs_alt)
                    + [r["C_tail"] for r in gen_rows]
                    + [r["C_tail"] for r in gen_rows_alt]
                )
                pad = 0.08 * max(max(abs(x) for x in b_all), max(abs(x) for x in c_all), 0.2)
                lim_lo = min(min(b_all), min(c_all)) - pad
                lim_hi = max(max(b_all), max(c_all)) + pad
                tlin = np.linspace(lim_lo, lim_hi, 120)
                colors = {"e": "tab:blue", "mu": "tab:orange", "tau": "tab:green"}

                def decorate_panel(ax, Bx, Cx, grows, subtitle):
                    ax.scatter(
                        Bx, Cx, c=g0s, cmap="viridis", s=10, alpha=0.35, norm=norm
                    )
                    ax.plot(tlin, tlin, "k--", lw=1.0, alpha=0.7)
                    ax.axhline(0.0, color="gray", ls=":", lw=0.9, alpha=0.8)
                    ax.axvline(0.0, color="gray", ls="-.", lw=0.9, alpha=0.8)
                    for row in grows:
                        ax.scatter(
                            [row["B_tail"]],
                            [row["C_tail"]],
                            s=120,
                            marker="*",
                            c=colors.get(row["gen"], "red"),
                            edgecolors="black",
                            linewidths=0.6,
                            zorder=5,
                            label=row["gen"],
                        )
                    ax.set_xlabel("B_tail")
                    ax.set_ylabel("C_tail")
                    ax.set_title(subtitle, fontsize=10)
                    ax.set_aspect("equal", adjustable="box")
                    ax.grid(True, alpha=0.3)
                    ax.set_xlim(lim_lo, lim_hi)
                    ax.set_ylim(lim_lo, lim_hi)

                decorate_panel(
                    axL,
                    Bs,
                    Cs,
                    gen_rows,
                    f"okno bazowe: r ∈ [{R_TAIL_L}, {R_TAIL_R}]",
                )
                decorate_panel(
                    axR,
                    Bs_alt,
                    Cs_alt,
                    gen_rows_alt,
                    f"okno alt: r ∈ [{R_TAIL_L_ALT}, {R_TAIL_R_ALT}]",
                )
                hL, labL = axL.get_legend_handles_labels()
                axL.legend(hL, labL, loc="upper left", fontsize=7, framealpha=0.9)
                hR, labR = axR.get_legend_handles_labels()
                axR.legend(hR, labR, loc="upper left", fontsize=7, framealpha=0.9)

                sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
                sm.set_array([])
                cbar = fig3.colorbar(sm, ax=[axL, axR], fraction=0.035, pad=0.02)
                cbar.set_label("g₀")

                fig3.suptitle(
                    "Płaszczyzna (B, C): porównanie okien fitu (wspólna skala osi)",
                    fontsize=11,
                    y=1.02,
                )
                plt.tight_layout()
                png3 = OUT_DIR / "ex114_tail_phase_generations_twowin.png"
                fig3.savefig(png3, dpi=140, bbox_inches="tight")
                plt.close(fig3)
                print(f"  Wykres: {png3}")
    except Exception as ex:
        print(f"  (matplotlib pominięty: {ex})")
    print()

    # --- testy ---
    As_arr = np.array([r["A_tail"] for r in rows])
    check(bool(csv_path.is_file() and csv_path.stat().st_size > 80), "P1: CSV zapisany", str(csv_path))
    check(bool(np.all(As_arr > 0)), "P2: A_tail > 0 na całej siatce", f"min(A)={np.min(As_arr):.6g}")

    if g0_star is not None:
        A_ref, B_ref, C_ref, _ = atail_bc(g0_star)
        row_star = next(r for r in rows if abs(r["g0"] - g0_star) < 1e-12)
        dA = abs(row_star["A_tail"] - A_ref) / A_ref
        dB = abs(row_star["B_tail"] - B_ref) / max(abs(B_ref), 1e-12)
        dC = abs(row_star["C_tail"] - C_ref) / max(abs(C_ref), 1e-12)
        check(
            dA < 1e-6 and dB < 1e-6 and dC < 1e-6,
            "P3: Ten sam g₀* w siatce = powtórny integrał (A,B,C)",
            f"rel|ΔA|={dA:.3e} rel|ΔB|={dB:.3e} rel|ΔC|={dC:.3e}",
        )
    else:
        check(False, "P3: g₀* — brak FP, test pominięty", "")

    check(max_rel_A < 0.08, "P4: Odporność okna |ΔA|/A < 8%", f"max={max_rel_A:.4f}")
    check(max_dphi < 0.35, "P5: Odporność okna |Δφ| < 0.35 rad", f"max={max_dphi:.4f}")
    check(max_rmse_ratio < 0.10, "P6: RMSE/A < 10% (jakość fitu w oknie)", f"max={max_rmse_ratio:.4f}")

    if g0_star is not None:
        ok_gen = (
            gen_csv_path.is_file()
            and gen_alt_csv_path.is_file()
            and gen_csv_path.stat().st_size > 50
            and gen_alt_csv_path.stat().st_size > 50
        )
        check(
            ok_gen,
            "P7: ex114_generations.csv + ex114_generations_altwin.csv",
            f"{gen_csv_path.name}, {gen_alt_csv_path.name}",
        )
        check(
            max_dphi_triad < 0.35,
            "P8: triada |Δφ| baz vs alt < 0.35 rad",
            f"max={max_dphi_triad:.4f} rad ({np.degrees(max_dphi_triad):.2f}°)",
        )
        # Próg: triada w strefie, gdzie tło jest wyraźnie bliskie próżni w oknie pomiarowym
        check(
            max_mean_abs_gm1_triad < 0.12,
            "P9: triada śr.|g-1| < 0.12 w oknie bazowym",
            f"max={max_mean_abs_gm1_triad:.4f}",
        )
    else:
        check(True, "P7: brak g₀* — pominięto CSV triady", "")
        check(True, "P8: brak g₀* — pominięto", "")
        check(True, "P9: brak g₀* — pominięto", "")

    n_fail = sum(1 for _l, s, _d in RESULTS if s == "FAIL")
    print()
    print("=" * 70)
    print(f"EX114: {len(RESULTS) - n_fail}/{len(RESULTS)} testów PASS")
    print("=" * 70)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
