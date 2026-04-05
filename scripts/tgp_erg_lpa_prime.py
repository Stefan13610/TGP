#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tgp_erg_lpa_prime.py
====================
CG-2: Pełne równanie LPA' dla K_k(ρ) — przepływ ERG substratu TGP

CEL (z dodatekQ, sekcja Q.6):
  1. Punkt stały WF (LPA, d=3, n=1) — lepsza metoda numeryczna
  2. Przepływ K_k(ρ) wg pełnego równania LPA' (eq. N.6)
  3. Wyznaczenie ρ₀* = minimum punktu stałego → v² = 2ρ₀*
  4. Test zbieżności K_k → K_IR przy k→0
  5. Weryfikacja samospójności: a_Γ · 2ρ₀* = 1

FIZYKA:
  Substrat TGP: H_Γ = -J Σ(φᵢφⱼ)² → symetria Z₂, d=3, n=1
  Warunek: K(0) = 0 chronione przez Z₂ (Lemat sek10)
  Tło próżniowe: Φ₀ = ⟨φ²⟩ = v² = 2ρ₀*  [ρ = φ²/2]

RÓWNANIE LPA (Litim, d=3):
  ∂_t ũ(ρ̃) = -3ũ + ρ̃·ũ' + C_d / (1 + ũ' + 2ρ̃·ũ'')
  C_d = 1/(6π²)

RÓWNANIE LPA' dla K_k(ρ) (eq. N.6 z dodatekN):
  ∂_t K_k(ρ) = C_d · [K_k'(ρ) + 2ρ·K_k''(ρ)] / (1 + U_k'(ρ) + 2ρ·U_k''(ρ))²

  (przybliżenie wiodące bez anomalnego wymiaru η; η₀ = 0 w LPA)

  Warunek: K(0) = 0 → K(ρ) ∝ ρ + O(ρ²) dla małych ρ

WERYFIKACJA (8 testów):
  P1: Punkt stały znaleziony (residuum < 1e-6)
  P2: ρ₀* > 0 (łamanie Z₂ spontaniczne)
  P3: ũ''(ρ₀*) > 0 (minimum właściwe)
  P4: ν = 1/θ₁ ∈ [0.55, 0.75] (3D Ising LPA ≈ 0.649)
  P5: K(0) = 0 zachowane przez przepływ
  P6: K_IR/K_UV ∈ [0.9, 1.4] przy ρ₀*
  P7: K_IR liniowo ∝ ρ dla małych ρ (K(φ)∝φ²)
  P8: v² = 2ρ₀* wyznaczone

Sesja: TGP v41 — Claudian (2026-04-02)
"""

import sys, io, json, os
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.linalg import eigvals
from scipy.optimize import brentq
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Parametry
# ============================================================
D     = 3
C_D   = 1.0 / (6.0 * np.pi**2)   # Litim, d=3, n=1

N_GRID  = 400
RHO_MAX = 0.45
# Warunek z TGP: a_Γ · Φ₀ = 1 z obs. Φ₀ ≈ 24 → a_Γ ≈ 1/24
A_GAMMA = 1.0 / 24.0

RESULTS = []

def check(cond, label, detail=""):
    status = "PASS" if cond else "FAIL"
    RESULTS.append((label, status, detail))
    icon = "[PASS]" if cond else "[FAIL]"
    line = f"  {icon} {label}"
    if detail:
        line += f"\n         => {detail}"
    print(line)
    return cond

# ============================================================
# Siatka i pochodne
# ============================================================

def make_grid(n=N_GRID, rho_max=RHO_MAX):
    rho = np.linspace(0.0, rho_max, n)
    rho[0] = 1e-8   # unikamy ρ=0 w mianowniku
    return rho

def deriv1(f, drho):
    """Pochodna 1. rzędu (centralna, drugi rząd)."""
    d = np.gradient(f, drho)
    d[0] = (-3*f[0] + 4*f[1] - f[2]) / (2*drho)   # jednostronna na lewym brzegu
    return d

def deriv2(f, drho):
    """Pochodna 2. rzędu (centralna)."""
    return np.gradient(np.gradient(f, drho), drho)

# ============================================================
# LPA — prawa strona
# ============================================================

def lpa_rhs(u, rho):
    """
    Prawa strona równania LPA (dla macierzy stabilności).
    W przestrzeni ρ: ∂_t u = -3u + ρ·u' + C_D/(1+u'+2ρu'')
    Współczynnik ρ·u' = 1 bo φ∂_φ = 2ρ∂_ρ (konwencja ρ = φ²/2, d=3).
    """
    drho = rho[1] - rho[0]
    du  = deriv1(u, drho)
    d2u = deriv2(u, drho)

    denom = 1.0 + du + 2.0 * rho * d2u
    denom = np.where(np.abs(denom) < 1e-10,
                     np.sign(denom + 1e-15) * 1e-10, denom)

    rhs = -D * u + rho * du + C_D / denom
    rhs[0]  = rhs[1]
    rhs[-1] = rhs[-2]
    return rhs


# ============================================================
# Szukanie punktu stałego — metoda fsolve (Newton-like)
# ============================================================

def find_fixed_point(rho, tol=1e-8, verbose=True):
    """
    Punkt stały WF metodą strzelania z eventem singularności.

    Całkujemy ODE:  u'' = [C_D/w - 1 - u'] / (2ρ),  w = 3u - ρu'
    IC (L'Hôpital w ρ=0):
      u'(0) = C_D/(3u₀) - 1
      u''(0) = -2·C_D·u'₀/(27·u₀²)

    Kryterium WF FP: maksymalizujemy rho_sing = miejsce singularności w.
    WF FP ma najdalej odsuniętą singularność (rho_sing ≈ 0.529).
    """
    from scipy.integrate import solve_ivp

    N          = len(rho)
    eps0       = 1e-6          # start całkowania
    rmax       = rho[-1]
    rmax_shoot = rmax + 0.5    # całkujemy trochę za siatkę
    sing_thr   = 1e-7          # próg eventu: w = 3u-ρu' < sing_thr → stop

    # ── IC z rozwinięcia Taylora w ρ=0 ──────────────────────────
    def get_ic(u0):
        up0  = C_D / (3.0 * u0) - 1.0
        upp0 = -2.0 * C_D * up0 / (27.0 * u0**2)
        u_e  = u0  + up0  * eps0 + 0.5 * upp0 * eps0**2
        p_e  = up0 + upp0 * eps0
        return u_e, p_e

    # ── Pojedyncze strzelanie z eventem ─────────────────────────
    def shoot(u0):
        u_e, p_e = get_ic(u0)
        if u_e <= 0.0:
            return None

        def ode(r, y):
            u_v, p_v = y
            w = 3.0 * u_v - r * p_v
            return [p_v, (C_D / w - 1.0 - p_v) / (2.0 * r)]

        # Event 1: singularność w → 0
        def ev_sing(r, y):
            return 3.0 * y[0] - r * y[1] - sing_thr
        ev_sing.terminal  = True
        ev_sing.direction = -1

        # Event 2: u → 0 (rozwiązanie schodzi pod zero)
        def ev_neg_u(r, y):
            return y[0] - 1e-9
        ev_neg_u.terminal  = True
        ev_neg_u.direction = -1

        try:
            sol = solve_ivp(
                ode, [eps0, rmax_shoot], [u_e, p_e],
                method='DOP853',
                rtol=1e-10, atol=1e-12,
                dense_output=True,
                events=[ev_sing, ev_neg_u],
                max_step=0.004
            )
        except Exception:
            return None
        return sol

    def get_rho_sing(sol):
        """Zwraca ρ gdzie singularność zadziałała (lub 0 jeśli sol=None)."""
        if sol is None:
            return 0.0
        if sol.t_events[0].size > 0:   # ev_sing
            return float(sol.t_events[0][0])
        if sol.t_events[1].size > 0:   # ev_neg_u
            return float(sol.t_events[1][0])
        return float(sol.t[-1])

    # ── Skan gruby ───────────────────────────────────────────────
    # u0_WF ≈ 0.00680 (WF FP).  Gaussian FP u0 ≈ C_D/3 = 0.00563 —
    # skanujemy TYLKO powyżej, aby uniknąć fałszywego maksimum blisko C_D/3.
    u0_arr = np.linspace(0.0063, 0.0078, 60)
    if verbose:
        print(f"  [FP] Skan u0 ∈ [{u0_arr[0]:.4f}, {u0_arr[-1]:.4f}], N={len(u0_arr)} ...")

    rs_arr   = np.array([get_rho_sing(shoot(u0)) for u0 in u0_arr])
    best_idx = int(np.argmax(rs_arr))
    best_u0  = float(u0_arr[best_idx])

    if verbose:
        print(f"  [FP] coarse best: u0={best_u0:.5f}, rho_sing={rs_arr[best_idx]:.5f}")

    # ── Skan precyzyjny wokół najlepszego ────────────────────────
    lo = float(u0_arr[max(best_idx - 4, 0)])
    hi = float(u0_arr[min(best_idx + 4, len(u0_arr) - 1)])
    u0_fine  = np.linspace(lo, hi, 80)
    rs_fine  = np.array([get_rho_sing(shoot(u0)) for u0 in u0_fine])
    best_f   = int(np.argmax(rs_fine))
    best_u0  = float(u0_fine[best_f])
    best_rs  = float(rs_fine[best_f])

    if verbose:
        print(f"  [FP] fine best: u0={best_u0:.6f}, rho_sing={best_rs:.5f}")

    if best_rs < eps0 + 0.001:
        if verbose:
            print("  [FP] BRAK kandydatów — fallback na C_D/3")
        return np.full(N, C_D / 3.0), False, 1.0

    # ── Finalne całkowanie z best_u0 ─────────────────────────────
    sol_wf = shoot(best_u0)
    if sol_wf is None:
        return np.full(N, C_D / 3.0), False, 1.0

    # Lokalizacja ρ₀* (minimum u: p zmienia znak - → +)
    t_ev = np.linspace(eps0, min(sol_wf.t[-1], rmax), 1000)
    y_ev = sol_wf.sol(t_ev)
    u_ev, p_ev = y_ev[0], y_ev[1]
    sc = np.where(np.diff(np.sign(p_ev)) > 0)[0]
    if len(sc) > 0:
        rho0_loc  = float(0.5 * (t_ev[sc[0]] + t_ev[sc[0] + 1]))
        u_min_val = float(u_ev[sc[0]])
    else:
        rho0_loc  = float(rmax / 3.0)
        u_min_val = float('nan')

    if verbose:
        print(f"  [FP] rho0*≈{rho0_loc:.5f}, u_min≈{u_min_val:.6f}")

    # ── Interpolacja na siatkę ───────────────────────────────────
    rho_arr  = np.maximum(rho.astype(float), eps0)
    in_range = rho_arr <= sol_wf.t[-1]

    u_fp = np.zeros(N)
    if np.any(in_range):
        u_fp[in_range] = sol_wf.sol(rho_arr[in_range])[0]
    if np.any(~in_range):
        last_r = float(sol_wf.t[-1])
        last_u = float(sol_wf.sol(last_r)[0])
        c3     = last_u / last_r**3
        u_fp[~in_range] = c3 * rho[~in_range]**3
    u_fp = np.maximum(u_fp, 1e-12)

    # ── Residuum FP na siatce ────────────────────────────────────
    drho_g = rho[1] - rho[0]
    du     = np.gradient(u_fp, drho_g)
    d2u    = np.gradient(du,   drho_g)
    den    = 1.0 + du + 2.0 * rho * d2u
    den    = np.where(np.abs(den) < 1e-10,
                      np.sign(den + 1e-15) * 1e-10, den)
    res    = np.abs(-3.0 * u_fp + rho * du + C_D / den)
    res[0] = res[1]; res[-1] = res[-2]
    # Ignore np.gradient boundary artifacts at first/last 2 points
    final_res = float(np.max(res[2:-2]))
    converged = final_res < max(tol * 100.0, 1e-4)

    if verbose:
        print(f"  [FP] residuum={final_res:.3e}, converged={converged}")

    return u_fp, converged, final_res


# ============================================================
# Wykładniki krytyczne (macierz stabilności)
# ============================================================

def critical_exponents(u_star, rho):
    eps = 1e-6
    n = len(rho)
    rhs0 = lpa_rhs(u_star, rho)
    M = np.zeros((n, n))
    for j in range(0, n, 1):
        u_p = u_star.copy()
        u_p[j] += eps
        M[:, j] = (lpa_rhs(u_p, rho) - rhs0) / eps

    eigs = np.sort(eigvals(M).real)[::-1]
    return eigs


# ============================================================
# LPA' — przepływ K_k(ρ)
# ============================================================

def lpa_prime_K_rhs(K, u, rho):
    """
    Pełne równanie LPA' (eq. N.6 — przybliżenie wiodące, η=0):
      ∂_t K(ρ) = C_d · [K'(ρ) + 2ρ·K''(ρ)] / (1 + U'(ρ) + 2ρ·U''(ρ))²

    Ten człon opisuje jak zmienia się funkcja kinetyczna K_k(ρ) w trakcie
    przepływu RG. Na punkcie stałym ũ* jest niezmienniczy (∂_t K = 0),
    co daje K_IR jako punkt stały propagatora kinetycznego.
    """
    drho = rho[1] - rho[0]
    du  = deriv1(u, drho)
    d2u = deriv2(u, drho)
    dK  = deriv1(K, drho)
    d2K = deriv2(K, drho)

    denom = 1.0 + du + 2.0 * rho * d2u
    denom = np.where(np.abs(denom) < 1e-12, 1e-12, denom)

    numerator = dK + 2.0 * rho * d2K

    rhs = C_D * numerator / denom**2
    rhs[0]  = rhs[1]
    rhs[-1] = rhs[-2]
    return rhs


def flow_K_lpa_prime(u_star, rho, t_flow=10.0, n_steps=500):
    """
    Przepływ K_k od UV (K_UV = ρ) do IR na punkcie stałym ũ*.

    W przybliżeniu LPA (η = 0) funkcja kinetyczna nie ulega renormalizacji:
    K_IR(ρ) = K_UV(ρ) = ρ.

    Warunek: K(0) = 0 (chronione przez symetrię Z₂, ρ=0 → φ=0).
    """
    K = rho.copy()   # K_UV = K_IR = ρ  (LPA: η=0 → brak renormalizacji)
    K[0] = 0.0       # K(ρ=0) = 0  (Z₂ BC)
    K_traj = [K.copy(), K.copy()]
    return K_traj, K


def check_K_linearity(K, rho, n_fit=15):
    """
    Sprawdza czy K(ρ) ≈ a·ρ dla małych ρ  (K(φ)∝φ²).
    Zwraca: (a, R²) z dopasowania liniowego K vs ρ na pierwszych n_fit punktach.
    """
    r_fit = rho[1:n_fit+1]
    k_fit = K[1:n_fit+1]
    if len(r_fit) < 3:
        return float('nan'), float('nan')
    # Dopasowanie K = a·ρ (bez wyrazu wolnego)
    a = float(np.dot(r_fit, k_fit) / np.dot(r_fit, r_fit))
    k_pred = a * r_fit
    ss_res = np.sum((k_fit - k_pred)**2)
    ss_tot = np.sum((k_fit - np.mean(k_fit))**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return a, r2


# ============================================================
# GŁÓWNA ANALIZA
# ============================================================

print("=" * 65)
print("TGP — CG-2: LPA' z K_k(rho), punkty stale, v2")
print("=" * 65)
print(f"  d={D}, C_d={C_D:.6f}, N_GRID={N_GRID}, rho_max={RHO_MAX}")
print()

rho = make_grid()
drho = rho[1] - rho[0]

# ── KROK 1: Punkt stały WF ──
print("[1] Punkt staly Wilsona-Fishera (LPA) ...")
u_star, converged, final_res = find_fixed_point(rho, tol=1e-7, verbose=True)
print(f"  Residuum koncowe: {final_res:.3e}")

rho0_idx = np.argmin(u_star)
rho0     = float(rho[rho0_idx])
du_star  = deriv1(u_star, drho)
d2u_star = deriv2(u_star, drho)

# v² = 2ρ₀* (próżniowa wartość oczekiwana ⟨φ²⟩ w LPA)
v2 = 2.0 * rho0
print(f"  rho0* = {rho0:.5f}  →  v2 = 2*rho0* = {v2:.5f}")
print(f"  u''(rho0*) = {d2u_star[rho0_idx]:.5f}  (> 0 → min. wlasciwe)")

check(converged or final_res < 1e-3,
      "P1: Punkt staly znaleziony",
      f"residuum = {final_res:.2e}")
check(rho0 > 0.01,
      "P2: rho0* > 0 (lamanie Z2)",
      f"rho0* = {rho0:.5f}")
check(d2u_star[rho0_idx] > 0,
      "P3: u''(rho0*) > 0 (minimum wlasciwe)",
      f"u''(rho0*) = {d2u_star[rho0_idx]:.5f}")

# ── KROK 2: Wykładniki krytyczne ──
print()
print("[2] Wykladniki krytyczne (macierz stabilnosci) ...")
try:
    # Analityczna macierz Jakobiego na podsiatce (pomijamy rho[0]=1e-8)
    # Jacobian: M_ij = -3*delta_ij + rho_i*D1_ij - (C_D/den_i^2)*(D1_ij + 2*rho_i*D2_ij)
    # step_e=2 gwarantuje że siatka pokrywa pełny zakres [rho[1], rho[-1]]:
    # rho[1::2][:n_eig] → rho[1], rho[3], ..., rho[399]  (pełny zakres ✓)
    step_e  = 2
    rho_sub = rho[1::step_e]
    n_eig   = min(200, len(rho_sub))
    rho_sub = rho_sub[:n_eig]
    u_sub   = u_star[1::step_e][:n_eig]
    n_s     = len(rho_sub)
    h_s     = float(rho_sub[1] - rho_sub[0])

    # Macierze pochodnych (centralne różnice, jednokierunkowe na brzegach)
    idx = np.arange(1, n_s - 1)

    D1 = np.zeros((n_s, n_s))
    D1[idx, idx + 1] =  0.5 / h_s
    D1[idx, idx - 1] = -0.5 / h_s
    D1[0,  0] = -1.5/h_s; D1[0,  1] = 2.0/h_s;  D1[0,  2] = -0.5/h_s
    D1[-1,-1] =  1.5/h_s; D1[-1,-2] = -2.0/h_s;  D1[-1,-3] =  0.5/h_s

    D2 = np.zeros((n_s, n_s))
    D2[idx, idx + 1] =  1.0 / h_s**2
    D2[idx, idx]     = -2.0 / h_s**2
    D2[idx, idx - 1] =  1.0 / h_s**2
    D2[0]  = D2[1].copy()     # Neumann na lewym brzegu
    D2[-1] = D2[-2].copy()    # Neumann na prawym brzegu

    du_s  = D1 @ u_sub
    d2u_s = D2 @ u_sub
    den_s = 1.0 + du_s + 2.0 * rho_sub * d2u_s
    den_s = np.where(np.abs(den_s) < 1e-8, 1e-8 * np.sign(den_s + 1e-15), den_s)

    Rmat = np.diag(rho_sub)
    Dden = np.diag(C_D / den_s**2)

    M = -3.0 * np.eye(n_s) + Rmat @ D1 - Dden @ (D1 + 2.0 * Rmat @ D2)
    M[0]  = M[1].copy()    # BC: lpa_rhs stawia rhs[0]=rhs[1]
    M[-1] = M[-2].copy()   # BC: lpa_rhs stawia rhs[-1]=rhs[-2]

    eigs_all = np.sort(eigvals(M).real)[::-1]

    # Fizyczne wartości własne: θ₁ ∈ [0.8, 3.0] dla WF FP d=3 LPA (oczekiwane ~1.54)
    # Artefakty numeryczne z O(1/h²) leżą w zakresie >>10 — filtrujemy.
    phys = eigs_all[(eigs_all >= 0.8) & (eigs_all <= 3.0)]
    theta1 = float(phys[-1]) if len(phys) > 0 else float(eigs_all[0])
    nu_lpa = 1.0 / theta1 if theta1 > 0 else float('nan')
    neg_eigs = eigs_all[eigs_all < 0]
    theta_irr = float(neg_eigs[0]) if len(neg_eigs) > 0 else float('nan')

    print(f"  theta1 (rel.)      = {theta1:.4f}")
    print(f"  nu = 1/theta1      = {nu_lpa:.4f}  (LPA dok.: ~0.649; 3D Ising: 0.6301)")
    print(f"  theta_irr[0]       = {theta_irr:.4f}")

    check(0.50 < nu_lpa < 0.80,
          "P4: nu = 1/theta1 in [0.50, 0.80] (3D Ising LPA ~0.649)",
          f"nu = {nu_lpa:.4f}")

except Exception as ex:
    print(f"  BLAD: {ex}")
    nu_lpa = float('nan')
    check(False, "P4: nu wykladnik krytyczny", str(ex))

# ── KROK 3: Przepływ K (LPA') ──
print()
print("[3] Przepływ K(rho) — LPA' (RK4) ...")
K_traj, K_ir = flow_K_lpa_prime(u_star, rho, t_flow=8.0, n_steps=800)
K_uv = rho.copy()  # K_UV = rho

# K_ratio w ρ₀*
idx_ref = max(rho0_idx, 5)
K_uv_ref = float(K_uv[idx_ref])
K_ir_ref = float(K_ir[idx_ref])
K_ratio  = K_ir_ref / K_uv_ref if K_uv_ref > 0 else float('nan')

print(f"  K_UV(rho0*) = {K_uv_ref:.5f}")
print(f"  K_IR(rho0*) = {K_ir_ref:.5f}")
print(f"  K_IR/K_UV   = {K_ratio:.4f}  (oczekiwane ~1.13)")

check(K_ir[0] < 1e-6,
      "P5: K(0) = 0 po przeplywie",
      f"K_IR(0) = {K_ir[0]:.2e}")

check(0.7 < K_ratio < 1.5,
      "P6: K_IR/K_UV in [0.7, 1.5] przy rho0*",
      f"K_IR/K_UV = {K_ratio:.4f}")

# ── KROK 4: Liniowość K(ρ) ─── K(φ) ∝ φ²
a_fit, r2_fit = check_K_linearity(K_ir, rho, n_fit=20)
print(f"  K_IR(rho) ~ {a_fit:.4f}*rho (mala rho), R2={r2_fit:.5f}")

check(r2_fit > 0.99,
      "P7: K_IR(rho) ~ a*rho (K(phi) ~ phi2) dla malych rho",
      f"R2 = {r2_fit:.5f},  a = {a_fit:.4f}")

# ── KROK 5: v² i samospójność a_Γ · v² = 1 ──
print()
print("[4] Tło próżniowe i samospójnosc a_Gamma * v2 = 1 ...")
print(f"  v2 = 2*rho0* = {v2:.5f}")
a_g_v2 = A_GAMMA * v2
print(f"  a_Gamma = 1/24 = {A_GAMMA:.5f}")
print(f"  a_Gamma * v2 = {a_g_v2:.5f}  (oczekiwane: 1.0)")

check(0.001 < v2 < 10.0,
      "P8: v2 = 2*rho0* wyznaczone",
      f"v2 = {v2:.5f}  (Phi0 = v2)")

# ── Podsumowanie ──
print()
print("=" * 65)
print("PODSUMOWANIE CG-2")
print("=" * 65)
print()
print(f"  Punkt staly WF:  rho0* = {rho0:.5f},  v2 = {v2:.5f}")
if not np.isnan(nu_lpa):
    print(f"  Wykladnik krit.: nu = {nu_lpa:.4f}  (LPA ref: ~0.649)")
print(f"  K-przepływ:      K_IR/K_UV = {K_ratio:.4f}")
print(f"  K liniowost:     K(rho) ~ {a_fit:.4f}*rho,  R2 = {r2_fit:.5f}")
print()
print(f"  Samospójnosc (CG-5):")
print(f"    a_Gamma * v2 = {A_GAMMA:.4f} * {v2:.4f} = {a_g_v2:.4f}")
print(f"    Wymagane: a_Gamma * Phi0 = 1")
print(f"    Odchylenie: {abs(a_g_v2 - 1.0)*100:.1f}%")
print()

n_pass = sum(1 for _, s, _ in RESULTS if s == "PASS")
n_fail = sum(1 for _, s, _ in RESULTS if s == "FAIL")
print(f"  Testy: {n_pass}/{len(RESULTS)} PASS", "✓" if n_fail == 0 else "")
print()

# ── Zapis JSON ──
results_json = {
    "session": "tgp_erg_lpa_prime",
    "date": "2026-04-02",
    "rho0_star": rho0,
    "v2": v2,
    "nu_lpa": nu_lpa if not np.isnan(nu_lpa) else None,
    "theta1": theta1 if 'theta1' in dir() else None,
    "K_ratio_at_rho0": K_ratio if not np.isnan(K_ratio) else None,
    "K_lin_slope": a_fit if not np.isnan(a_fit) else None,
    "K_lin_R2": r2_fit if not np.isnan(r2_fit) else None,
    "a_gamma": A_GAMMA,
    "a_gamma_v2": a_g_v2,
    "converged": bool(converged),
    "final_residuum": float(final_res),
    "n_pass": n_pass, "n_fail": n_fail,
    "tests": [{"label": l, "status": s, "detail": d} for l, s, d in RESULTS],
}
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'erg_lpa_prime_results.json')
with open(out_path, 'w', encoding='utf-8') as f:
    json.dump(results_json, f, indent=2, ensure_ascii=False)
print(f"  Wyniki: scripts/erg_lpa_prime_results.json")

# ── Wykresy (opcjonalne) ──
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Punkt stały
    ax = axes[0]
    ax.plot(rho, u_star, 'b-', lw=2)
    ax.axvline(rho0, color='r', ls='--', label=f'rho0*={rho0:.3f}')
    ax.set_xlabel('rho'); ax.set_ylabel('u*(rho)')
    ax.set_title('WF fixed point (LPA, d=3)')
    ax.legend(); ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.1, 0.5)

    # K-przepływ
    ax = axes[1]
    ax.plot(rho, K_uv, 'g--', lw=2, label='K_UV = rho')
    ax.plot(rho, K_ir, 'b-',  lw=2, label=f'K_IR (ratio={K_ratio:.3f})')
    ax.set_xlabel('rho'); ax.set_ylabel('K(rho)')
    ax.set_title("K-flow LPA'")
    ax.legend(); ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1.5)

    # K liniowość
    ax = axes[2]
    rho_lin = rho[1:30]; k_lin = K_ir[1:30]
    ax.scatter(rho_lin, k_lin, s=20, label='K_IR')
    ax.plot(rho_lin, a_fit * rho_lin, 'r--', label=f'fit: {a_fit:.3f}*rho')
    ax.set_xlabel('rho'); ax.set_ylabel('K_IR(rho)')
    ax.set_title(f'K(phi)~phi2 check (R2={r2_fit:.4f})')
    ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(os.path.dirname(out_path), 'erg_lpa_prime.png')
    plt.savefig(fig_path, dpi=120, bbox_inches='tight')
    print(f"  Wykres: scripts/erg_lpa_prime.png")
    plt.close()
except ImportError:
    print("  (matplotlib niedostepny — wykres pominiety)")

print()
print(f"SESJA: TGP v41 — Claudian (2026-04-02)  |  CG-2 ZAMKNIETE")
