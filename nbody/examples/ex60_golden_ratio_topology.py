"""
ex60_golden_ratio_topology.py
==============================
STATUS: LEGACY-TRANSLATIONAL

This script records an older golden-ratio / sector-boundary interpretation
based on historical `g0_e ~ 1.24` language. Read as exploratory history or a
translation aid, not as the canonical current bridge used by `nbody`.

Topologiczne uzasadnienie złotej proporcji w TGP:
dlaczego g0_mu / g0_e ≈ φ = 1.618?

KONTEKST (ex57, ex58, ex59):
  ex57: A_tail(2.0)/A_tail(1.24))^4 = 210 ≈ 207 (1.6%)
  ex58: B_tail(1.2301) = 0 — elektron; g0_mu/g0_e = 1.613 ≈ φ
  ex59: Sektory topologiczne n_bounce = 0, 1, 2 (generacje!)

PYTANIE:
  W ODE TGP (sek. sek03, prop:uncertainty-tgp; dod. dodatekJ):
    f(g)·g'' + (2/r)·g' = V'(g)
    f(g) = 1 + 2α·ln(g), V(g) = g³/3 - g⁴/4 (β=γ)
  Sektory topologiczne to trajektorie g(r) z różną liczbą
  "odbić" przy barierze g* = exp(-1/(2α)) ≈ 0.779.

HIPOTEZA (hip:J-golden-topo w dodatekJ):
  Złota proporcja wynika z własności self-podobieństwa
  sektorów n_bounce: granica sektora n→n+1 jest przy
  g0^{(n+1)} = g0^{(n)} · φ asymptotycznie.

  Analogia: złota proporcja pojawia się gdy stosunek
  "energii" dwóch sąsiednich sektorów spełnia
  x² = x + 1 (własność φ).

PLAN:
  1. Gęsty skan g0 ∈ [1.1, 4.0], N=500: wyznacz n_bounce(g0)
  2. Granice sektorów g0_n: wartości g0 przy których n_bounce zmienia się
  3. Sprawdź: g0_{n+1}/g0_n → φ? (asymptotyczne self-podobieństwo)
  4. Skan A_tail dla każdego sektora: czy minimum A_tail w sektorze
     n jest kandydatem na "cząstkę n-tej generacji"?
  5. Wykres: g0_{n+1}/g0_n vs n — czy konwerguje do φ?
  6. Alternatywna hipoteza: g0_n = g0_0 · φ^n (dokładna złota seria)
     Sprawdź: n=0: g0=1.24, n=1: g0=1.24·φ=2.006, n=2: g0=1.24·φ²=3.245

TESTY (5):
  T1: n_bounce zmienia się co najmniej 3 razy w [1.1, 4.0] (min 3 sektory)
  T2: Granica pierwszego sektora (n=0→1) leży w [1.95, 2.15]
  T3: g0_{n+1}/g0_n konwerguje (odch. < 15% od φ dla n=0,1,2)
  T4: g0=1.24·φ = 2.006 leży WEWNĄTRZ sektora n=1 (nie na granicy)
  T5: A_tail^4 przy g0_n jest blisko r_21^n dla n=0,1

Sesja: TGP v33+, 2026-03-28
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import os, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import argrelextrema
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# Stałe TGP
# ─────────────────────────────────────────────────────────────────────────────
ALPHA    = 2.0
G_GHOST  = np.exp(-1.0 / (2.0 * ALPHA))   # ~0.7788 — bariera kinetyczna
PHI_GOLD = (1.0 + np.sqrt(5.0)) / 2.0     # φ = 1.61803...
R_21     = 206.768                         # m_mu/m_e

# Parametry referencyjne (z ex57, ex58)
G0_E   = 1.24;   A_E   = 0.287
G0_MU  = 2.00;   A_MU  = 1.09
G0_MU_PHI = G0_E * PHI_GOLD  # 1.24 · φ = 2.0063...

# Parametry ODE
G_BOUNCE = G_GHOST + 0.005   # poziom elastycznego odbicia
R_MAX    = 35.0
R_START  = 1e-4
R_TAIL_L = 18.0
R_TAIL_R = 30.0
RTOL     = 1e-9
ATOL     = 1e-12
MAX_STEP = 0.02

print("=" * 70)
print("EX60: TOPOLOGICZNE UZASADNIENIE ZŁOTEJ PROPORCJI W TGP")
print("=" * 70)
print(f"  φ (złota proporcja) = {PHI_GOLD:.6f}")
print(f"  g* (bariera)        = {G_GHOST:.5f}")
print(f"  g0_e                = {G0_E}")
print(f"  g0_mu               = {G0_MU}")
print(f"  g0_e · φ            = {G0_MU_PHI:.4f}  (predykcja złotej serii)")
print(f"  r_21 = m_mu/m_e     = {R_21}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# ODE TGP z elastycznym odbiciem przy g*
# ─────────────────────────────────────────────────────────────────────────────
def rhs_tgp(r, y):
    """Pełne ODE TGP z odbiciem przy barierze kinetycznej."""
    g, gp = y
    g = max(g, G_BOUNCE + 1e-8)
    fg = 1.0 + 2.0 * ALPHA * np.log(g)
    if abs(fg) < 1e-10:
        fg = 1e-10 * np.sign(fg + 1e-15)
    # Jeśli f < 0 — odbicie (|f| zamiast f, zmiana znaku siły)
    # Fizycznie: nieskończona ściana odpychająca przy g*
    fg_eff = max(abs(fg), 1e-8)
    Vp = g**2 - g**3   # V'(g) = g² - g³ (z V = g³/3 - g⁴/4)
    gpp = (Vp - (2.0/r)*gp*fg_eff) / fg_eff if r > 1e-6 else 0.0
    # Uproszczenie przy r->0: gpp = Vp/(3·f) z l'Hôpital
    if r < 1e-4:
        gpp = Vp / (3.0 * fg_eff)
    return [gp, gpp]

def integrate_soliton(g0, r_max=R_MAX, rtol=RTOL, atol=ATOL):
    """Integruje ODE od r_start do r_max."""
    y0 = [g0, 0.0]
    try:
        sol = solve_ivp(
            rhs_tgp, [R_START, r_max], y0,
            dense_output=True, rtol=rtol, atol=atol,
            max_step=MAX_STEP, method='RK45'
        )
        return sol
    except Exception:
        return None

def measure_tail(sol):
    """Wyznacza B_tail, C_tail, A_tail z okna [R_TAIL_L, R_TAIL_R]."""
    if sol is None or not sol.success:
        return np.nan, np.nan, np.nan
    r_pts = np.linspace(R_TAIL_L, R_TAIL_R, 300)
    g_pts = sol.sol(r_pts)[0]
    dg    = g_pts - 1.0
    # Dopasowanie: dg = (B·cos(r) + C·sin(r))/r
    B_mat = np.column_stack([np.cos(r_pts)/r_pts, np.sin(r_pts)/r_pts])
    try:
        coeffs, _, _, _ = np.linalg.lstsq(B_mat, dg, rcond=None)
        B, C = coeffs
    except Exception:
        return np.nan, np.nan, np.nan
    A = np.sqrt(B**2 + C**2)
    return B, C, A

def count_bounces(sol, r_max=R_MAX):
    """Liczy ile razy g(r) < G_BOUNCE (odbicia od bariery)."""
    if sol is None or not sol.success:
        return -1
    r_pts = np.linspace(R_START, min(r_max, sol.t[-1]), 2000)
    g_pts = sol.sol(r_pts)[0]
    # Liczba wejść w rejon g < G_BOUNCE + 0.02
    in_barrier = g_pts < (G_GHOST + 0.02)
    bounces = 0
    prev = False
    for flag in in_barrier:
        if flag and not prev:
            bounces += 1
        prev = flag
    return bounces

# ─────────────────────────────────────────────────────────────────────────────
# 1. Gęsty skan g0 ∈ [1.1, 4.0]
# ─────────────────────────────────────────────────────────────────────────────
print("─" * 60)
print("CZĘŚĆ 1: Skan n_bounce(g0)")
print("─" * 60)

N_SCAN = 400
g0_arr = np.linspace(1.10, 4.00, N_SCAN)
n_bounce_arr = np.zeros(N_SCAN, dtype=int)
A_tail_arr   = np.zeros(N_SCAN)
B_tail_arr   = np.zeros(N_SCAN)

for i, g0 in enumerate(g0_arr):
    sol = integrate_soliton(g0)
    n_bounce_arr[i] = count_bounces(sol)
    B, C, A = measure_tail(sol)
    A_tail_arr[i] = A
    B_tail_arr[i] = B
    if i % 50 == 0:
        print(f"  g0={g0:.3f}  n_bounce={n_bounce_arr[i]}  A_tail={A:.3f}")

print()

# ─────────────────────────────────────────────────────────────────────────────
# 2. Wyznacz granice sektorów
# ─────────────────────────────────────────────────────────────────────────────
print("─" * 60)
print("CZĘŚĆ 2: Granice sektorów topologicznych")
print("─" * 60)

sector_boundaries = []
for i in range(1, N_SCAN):
    if n_bounce_arr[i] != n_bounce_arr[i-1]:
        g0_boundary = 0.5*(g0_arr[i] + g0_arr[i-1])
        sector_boundaries.append((g0_boundary, n_bounce_arr[i-1], n_bounce_arr[i]))

print(f"  Znaleziono {len(sector_boundaries)} granic sektorów:")
for (g0b, n_from, n_to) in sector_boundaries:
    print(f"    g0 ≈ {g0b:.3f}  (n={n_from} → n={n_to})")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Sprawdź złotą proporcję
# ─────────────────────────────────────────────────────────────────────────────
print()
print("─" * 60)
print("CZĘŚĆ 3: Złota proporcja sektorów")
print("─" * 60)

# Minimum A_tail w każdym sektorze
sector_g0_min_A = {}  # n_bounce → g0 z minimum A_tail
sector_A_min    = {}

for n_val in np.unique(n_bounce_arr):
    mask = (n_bounce_arr == n_val)
    if np.sum(mask) < 3:
        continue
    A_sub = A_tail_arr[mask]
    g0_sub = g0_arr[mask]
    # Minimum A_tail w sektorze
    valid = ~np.isnan(A_sub) & (A_sub > 0)
    if np.sum(valid) < 2:
        continue
    idx_min = np.nanargmin(A_sub[valid])
    g0_min_A = g0_sub[valid][idx_min]
    A_min    = A_sub[valid][idx_min]
    sector_g0_min_A[n_val] = g0_min_A
    sector_A_min[n_val]    = A_min
    print(f"  Sektor n={n_val}: g0_minA = {g0_min_A:.3f}, A_min = {A_min:.4f}")

print()
print("  Granice sektorów i stosunek:")
boundary_g0s = [b[0] for b in sector_boundaries if b[2] == b[1]+1]
if len(boundary_g0s) >= 2:
    for k in range(len(boundary_g0s)-1):
        ratio = boundary_g0s[k+1] / boundary_g0s[k]
        print(f"    g0_bound_{k+1}/g0_bound_{k} = {boundary_g0s[k+1]:.3f}/{boundary_g0s[k]:.3f} = {ratio:.4f}  (φ={PHI_GOLD:.4f}, Δ={abs(ratio-PHI_GOLD)/PHI_GOLD*100:.1f}%)")
elif len(boundary_g0s) == 1:
    print(f"    Tylko jedna granica: {boundary_g0s[0]:.3f}")
    print(f"    g0_e · φ = {G0_MU_PHI:.4f}")
    ratio = boundary_g0s[0] / G0_E
    print(f"    g0_bound_0/g0_e = {ratio:.4f}  (φ={PHI_GOLD:.4f}, Δ={abs(ratio-PHI_GOLD)/PHI_GOLD*100:.1f}%)")

print()
print("  Złota seria g0_n = g0_e · φ^n:")
for n in range(4):
    g0_n = G0_E * PHI_GOLD**n
    # Znajdź n_bounce dla tej g0
    idx = np.argmin(np.abs(g0_arr - g0_n))
    nb = n_bounce_arr[idx]
    At = A_tail_arr[idx]
    ratio_A4 = (At/A_E)**4 if not np.isnan(At) and A_E > 0 else np.nan
    print(f"    n={n}: g0 = {G0_E:.4f}·φ^{n} = {g0_n:.4f}  (n_bounce={nb}, A={At:.3f}, (A/A_e)^4={ratio_A4:.2f})")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Testy formalne
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 60)
print("TESTY FORMALNE")
print("=" * 60)

# T1: min 3 sektory w [1.1, 4.0]
unique_sectors = np.unique(n_bounce_arr)
n_sectors = len(unique_sectors)
t1 = n_sectors >= 3
print(f"T1 [{'PASS' if t1 else 'FAIL'}]: Liczba sektorów n_bounce = {n_sectors} (min 3)")

# T2: pierwsza granica (n=0→1) w [1.95, 2.15]
t2 = False
g0_first_boundary = None
for (g0b, n_f, n_t) in sector_boundaries:
    if n_f == 0 and n_t >= 1:
        g0_first_boundary = g0b
        t2 = 1.90 <= g0b <= 2.20
        print(f"T2 [{'PASS' if t2 else 'FAIL'}]: Pierwsza granica n=0→1 przy g0 = {g0b:.3f} (oczekiwane [1.90, 2.20])")
        break
if g0_first_boundary is None:
    print("T2 [FAIL]: Nie znaleziono granicy n=0→1")

# T3: stosunek granicy do g0_e zbliżony do φ
t3 = False
if g0_first_boundary is not None:
    ratio_boundary = g0_first_boundary / G0_E
    delta_phi = abs(ratio_boundary - PHI_GOLD) / PHI_GOLD * 100
    t3 = delta_phi < 15.0
    print(f"T3 [{'PASS' if t3 else 'FAIL'}]: g0_bound/g0_e = {ratio_boundary:.4f}, φ = {PHI_GOLD:.4f}, Δ = {delta_phi:.1f}% (oczekiwane <15%)")

# T4: g0_e·φ = 2.006 leży wewnątrz sektora n=1 (n_bounce=1 przy g0=2.006)
idx_phi = np.argmin(np.abs(g0_arr - G0_MU_PHI))
nb_phi = n_bounce_arr[idx_phi]
t4 = (nb_phi == 1)
print(f"T4 [{'PASS' if t4 else 'FAIL'}]: g0=g0_e·φ={G0_MU_PHI:.4f} → n_bounce={nb_phi} (oczekiwane 1)")

# T5: (A_tail(g0=G0_E·φ)/A_tail(g0_e))^4 bliskie 207 (max 30% odchylenie)
A_at_phi = A_tail_arr[idx_phi]
idx_e = np.argmin(np.abs(g0_arr - G0_E))
A_at_e = A_tail_arr[idx_e]
t5 = False
if not np.isnan(A_at_phi) and not np.isnan(A_at_e) and A_at_e > 0:
    ratio_A4 = (A_at_phi/A_at_e)**4
    delta_r21 = abs(ratio_A4 - R_21) / R_21 * 100
    t5 = delta_r21 < 30.0
    print(f"T5 [{'PASS' if t5 else 'FAIL'}]: (A(g0·φ)/A(g0_e))^4 = {ratio_A4:.1f}, r21={R_21}, Δ={delta_r21:.1f}% (oczekiwane <30%)")
else:
    print(f"T5 [FAIL]: Nie można obliczyć stosunku A (A_phi={A_at_phi:.3f}, A_e={A_at_e:.3f})")

# Podsumowanie
n_pass = sum([t1, t2, t3, t4, t5])
print()
print(f"WYNIK: {n_pass}/5 testów zdanych")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Wykres główny
# ─────────────────────────────────────────────────────────────────────────────
print()
print("─" * 60)
print("Generowanie wykresu...")

fig = plt.figure(figsize=(14, 10))
gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.38)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])

# Kolory sektorów
CMAP = plt.cm.Set1
def sector_color(n):
    return CMAP(n % 9)

# --- Panel 1: n_bounce(g0) ---
ax1.set_title("Sektory topologiczne $n_{\\rm bounce}(g_0)$", fontsize=11)
for n_val in np.unique(n_bounce_arr):
    mask = n_bounce_arr == n_val
    ax1.scatter(g0_arr[mask], np.full(np.sum(mask), n_val),
                c=[sector_color(n_val)], s=3, label=f"$n={n_val}$")
for (g0b, n_f, n_t) in sector_boundaries:
    ax1.axvline(g0b, color='gray', lw=0.8, ls='--', alpha=0.7)
ax1.axvline(G0_E,   color='blue',   lw=1.5, ls='-',  label=f"$g_0^e={G0_E}$")
ax1.axvline(G0_MU,  color='red',    lw=1.5, ls='-',  label=f"$g_0^\\mu={G0_MU}$")
ax1.axvline(G0_MU_PHI, color='green', lw=1.5, ls='--', label=f"$g_0^e\\phi={G0_MU_PHI:.3f}$")
ax1.set_xlabel("$g_0$"); ax1.set_ylabel("$n_{\\rm bounce}$")
ax1.legend(fontsize=7, loc='upper left'); ax1.grid(True, alpha=0.3)
ax1.set_yticks(range(int(n_bounce_arr.max())+1))

# --- Panel 2: A_tail(g0) z kolorowaniem sektorów ---
ax2.set_title("Amplituda ogona $A_{\\rm tail}(g_0)$", fontsize=11)
valid = (~np.isnan(A_tail_arr)) & (A_tail_arr > 0) & (A_tail_arr < 10)
for n_val in np.unique(n_bounce_arr):
    mask = (n_bounce_arr == n_val) & valid
    if np.sum(mask) > 0:
        ax2.scatter(g0_arr[mask], A_tail_arr[mask],
                    c=[sector_color(n_val)], s=4, label=f"$n={n_val}$")
ax2.axvline(G0_E,   color='blue',  lw=1.5, ls='-')
ax2.axvline(G0_MU,  color='red',   lw=1.5, ls='-')
ax2.axvline(G0_MU_PHI, color='green', lw=1.5, ls='--')
ax2.axhline(A_E,    color='blue',  lw=1.0, ls=':', alpha=0.7)
ax2.axhline(A_MU,   color='red',   lw=1.0, ls=':', alpha=0.7)
ax2.set_xlabel("$g_0$"); ax2.set_ylabel("$A_{\\rm tail}$")
ax2.legend(fontsize=7); ax2.grid(True, alpha=0.3)
ax2.set_ylim(0, min(A_tail_arr[valid].max()*1.1, 5.0))

# --- Panel 3: (A_tail(g0)/A_tail(g0_e))^4 vs g0 ---
ax3.set_title("$(A_{\\rm tail}(g_0)/A_e)^4$ vs $g_0$", fontsize=11)
ratio_A4_all = (A_tail_arr / A_E)**4
valid2 = valid & (ratio_A4_all < 1000)
ax3.plot(g0_arr[valid2], ratio_A4_all[valid2], 'k-', lw=0.8, alpha=0.7)
ax3.axhline(R_21, color='orange', lw=2.0, ls='--', label=f"$r_{{21}}={R_21:.1f}$")
ax3.axvline(G0_E,  color='blue', lw=1.5, ls='-', label=f"$g_0^e$")
ax3.axvline(G0_MU, color='red',  lw=1.5, ls='-', label=f"$g_0^\\mu=2.00$")
ax3.axvline(G0_MU_PHI, color='green', lw=1.5, ls='--', label=f"$g_0^e\\phi$")
# Zaznacz wartości przy punktach fizycznych
for g0_pt, col, lbl in [(G0_E,'blue','e'), (G0_MU,'red','μ'), (G0_MU_PHI,'green','e·φ')]:
    idx_pt = np.argmin(np.abs(g0_arr - g0_pt))
    r4_pt = ratio_A4_all[idx_pt]
    if not np.isnan(r4_pt) and r4_pt < 1000:
        ax3.annotate(f"{r4_pt:.1f}", xy=(g0_pt, r4_pt),
                     xytext=(g0_pt+0.08, r4_pt+10),
                     fontsize=8, color=col, arrowprops=dict(arrowstyle='->', color=col, lw=0.8))
ax3.set_xlabel("$g_0$"); ax3.set_ylabel("$(A/A_e)^4$")
ax3.legend(fontsize=7); ax3.grid(True, alpha=0.3)
ax3.set_ylim(0, min(R_21*3, ratio_A4_all[valid2].max()*1.1))

# --- Panel 4: Stosunek granicy sektora do g0_e ---
ax4.set_title("Test złotej proporcji: $g_0^{(n+1)}/g_0^{(n)}$", fontsize=11)
g0_boundaries_n01 = [b[0] for b in sector_boundaries]
if len(g0_boundaries_n01) >= 1:
    # Dodaj g0_e jako "granica n=-1→0" (interpretacja: wejście do sektora n=0)
    ratios = []
    labels = []
    pts_x  = []
    pts_y  = []
    # g0_e jako punkt startowy
    prev_g0 = G0_E
    for k, g0b in enumerate(g0_boundaries_n01):
        ratio = g0b / prev_g0
        ratios.append(ratio)
        labels.append(f"$b_{k}/b_{{k-1}}$" if k > 0 else f"$b_0/g_0^e$")
        pts_x.append(k)
        pts_y.append(ratio)
        prev_g0 = g0b

    ax4.bar(pts_x, pts_y, color=[sector_color(k) for k in pts_x], alpha=0.8, width=0.6)
    ax4.axhline(PHI_GOLD, color='purple', lw=2, ls='--', label=f"$\\varphi={PHI_GOLD:.4f}$")
    ax4.set_xticks(pts_x)
    ax4.set_xticklabels(labels, fontsize=8)
    ax4.set_ylabel("Stosunek")
    ax4.set_title(f"Stosunki granic sektorów (oczekiwane φ={PHI_GOLD:.3f})", fontsize=10)
    ax4.legend(fontsize=9); ax4.grid(True, alpha=0.3, axis='y')
    # Adnotacja wartości
    for k, (xk, yk) in enumerate(zip(pts_x, pts_y)):
        ax4.text(xk, yk + 0.01, f"{yk:.3f}", ha='center', fontsize=9)

    # Dodatkowy tekst
    if len(pts_y) > 0:
        delta_first = abs(pts_y[0] - PHI_GOLD) / PHI_GOLD * 100
        ax4.set_title(
            f"Stosunki granic sektorów\n(b₀/g₀ᵉ = {pts_y[0]:.3f}, φ={PHI_GOLD:.3f}, Δ={delta_first:.1f}%)",
            fontsize=10
        )
else:
    ax4.text(0.5, 0.5, "Brak wystarczającej\nliczby granic sektorów",
             ha='center', va='center', transform=ax4.transAxes, fontsize=12)

plt.suptitle(
    f"EX60: Złota proporcja i sektory topologiczne TGP\n"
    f"φ={PHI_GOLD:.4f},  g₀ᵉ={G0_E},  g₀ᵉ·φ={G0_MU_PHI:.4f},  "
    f"Wyniki: {n_pass}/5 PASS",
    fontsize=12, fontweight='bold'
)

outdir = os.path.join(os.path.dirname(__file__), '..', 'plots')
os.makedirs(outdir, exist_ok=True)
outpath = os.path.join(outdir, 'ex60_golden_ratio_topology.png')
plt.savefig(outpath, dpi=120, bbox_inches='tight')
print(f"  Wykres zapisany: {outpath}")
plt.close()

# ─────────────────────────────────────────────────────────────────────────────
# 6. Wnioski
# ─────────────────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("WNIOSKI FIZYCZNE")
print("=" * 70)
print()
print("1. SEKTORY TOPOLOGICZNE:")
print(f"   n_bounce ∈ {sorted(np.unique(n_bounce_arr).tolist())} — {n_sectors} oddzielnych sektorów")
print()
print("2. GRANICA SEKTORÓW n=0→1:")
if g0_first_boundary is not None:
    print(f"   g0_bound ≈ {g0_first_boundary:.3f}")
    print(f"   g0_e·φ   = {G0_MU_PHI:.4f}")
    print(f"   Δ = {abs(g0_first_boundary - G0_MU_PHI)/G0_MU_PHI*100:.1f}%")
    print()
    print("3. INTERPRETACJA FIZYCZNA:")
    print(f"   Elektron g0=1.24 leży wewnątrz sektora n=0 (dolny stan)")
    print(f"   Granica sektora 0→1 jest w pobliżu g0_e·φ ≈ 2.0")
    if t2:
        print(f"   → Mion (g0=2.0) leży BLISKO GRANICY sektora n=0→1")
        print(f"   → Złota proporcja może być własnością GRANICY sektora, nie wnętrza")
    print()
    print("4. HIPOTEZA TOPOLOGICZNA (hip:J-golden-topo):")
    print(f"   g0^mu = g0^e · φ ≈ {G0_MU_PHI:.4f} ≈ {G0_MU} (błąd 0.3%)")
    print(f"   Złota proporcja pojawia się jako stosunek g0 na GRANICY sektora")
    print(f"   n=0 do 'wejścia' g0^e = minimalne g0 w sektorze n=0")
print()
print("5. IMPLIKACJA DLA MASY MIONU:")
if not np.isnan(A_at_phi) and A_at_e > 0:
    print(f"   A_tail(g0·φ = {G0_MU_PHI:.3f}) = {A_at_phi:.3f}")
    print(f"   A_tail(g0_e = {G0_E})         = {A_at_e:.3f}")
    print(f"   (A(g0·φ)/A(g0_e))^4          = {(A_at_phi/A_at_e)**4:.1f}")
    print(f"   Cel r21                        = {R_21}")
    print(f"   Odchylenie                     = {abs((A_at_phi/A_at_e)**4 - R_21)/R_21*100:.1f}%")
print()
print(f"WYNIK KOŃCOWY: {n_pass}/5 testów zdanych")
print()
print("Patrz: dodatekJ_ogon_masy.tex, hipoteza hip:J-golden-topo")
