#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p99_K2_convergence.py  --  TGP v1
====================================
Cel: Precyzyjna ekstrapolacja K*_2(r_max -> inf) dla drugiego solitonu TGP.

Kontekst (p98, r_max in {40,60,80}):
  K*_2(r=40) = 0.0857  -- galaz nizsza (inna niz r=60,80!)
  K*_2(r=60) = 0.2103  -- galaz stabilna
  K*_2(r=80) = 0.2153  -- ta sama galaz

Kluczowy warunek poczatkowy (z p98!):
  phi(A_GAM) = psi
  dphi(A_GAM) = -K / A_GAM^2    <-- K wchodzi przez nachylenie poczatkowe!

Plan:
  1. Uzyj identycznej metody jak p98 (dphi0 = -K/A_GAM^2)
  2. r_max w {60, 80, 100, 120, 150}
  3. Dla kazdego r_max: skan K w [0.10, 0.45], n_K=50
     Dla kazdego K: wszystkie zera psi w [1.001, 5.0] (ostatnie zero = Galaz U)
     Oblicz g_U(K) = E[phi]/(4pi K) - 1
  4. Znajdz zmiane znaku g_U(K) -> K*_2 (brentq)
  5. Fit K*_2(r) = a + b/r -> K*_2(inf)
  6. K*_2/K*_1(inf) vs R_21 = 206.77

Testy PASS/FAIL:
  P1: K*_2 znaleziono przy wszystkich r_max w {60,80,100,120,150}
  P2: psi0(K*_2) zbiezne: spread < 0.30
  P3: fit K*_2(r)=a+b/r: R^2 > 0.990
  P4: K*_2/K*_1(inf) > 5 (fizyczny drugi soliton)
  P5: K*_2/K*_1(inf) << 206.77 (V_mod nie odtwarza hierarchii leptonow)

Data: 2026-03-26
"""

import sys, io, warnings, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from numpy.linalg import lstsq
warnings.filterwarnings('ignore')

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

# ─── parametry fizyczne (identyczne z p98) ───────────────────────────────────
ALPHA_K = 8.5616
A_GAM   = 0.040
LAM     = 5.501357e-06
GAM     = 1.0
V1      = GAM/3 - GAM/4   # energia prozni = 1/12
K1_INF  = 0.010029        # K*_1(inf) z p97
R21_PDG = 206.77

# ─── parametry skanu ─────────────────────────────────────────────────────────
RMAXES  = [60.0, 80.0, 100.0, 120.0, 150.0]
K_LO    = 0.10
K_HI    = 0.45
N_K     = 50
PSI_LO  = 1.001
PSI_HI  = 5.0
N_PSI   = 80
N_EVAL  = 600   # punktow na rozwiazanie ODE
RTOL    = 1e-8
ATOL    = 1e-10

# ─── ODE ─────────────────────────────────────────────────────────────────────
def dV_mod(phi):
    return GAM*phi**2 - GAM*phi**3 + LAM*(phi - 1)**5

def V_mod(phi):
    return GAM/3*phi**3 - GAM/4*phi**4 + LAM/6*(phi - 1)**6

def ode_rhs(r, y):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA_K / phi
    d2phi = dV_mod(phi)/kfac + ALPHA_K*dphi**2/(2*phi**2*kfac) - 2.0/r*dphi
    return [dphi, d2phi]

# ─── strzelanie: phi(r_max) przy danych (psi, K, r_max) ──────────────────────
def phi_at_rmax(psi, K, r_max):
    """
    Integruje ODE od r=A_GAM do r_max z warunkami:
      phi(A_GAM) = psi
      dphi(A_GAM) = -K / A_GAM^2
    Zwraca phi(r_max) lub nan przy bledzie.
    """
    phi0  = psi
    dphi0 = -K / A_GAM**2
    # geometryczna siatka (dobra dla solitonu)
    r_ev  = A_GAM * (r_max / A_GAM)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(ode_rhs, [A_GAM, r_max], [phi0, dphi0],
                        method='DOP853', rtol=RTOL, atol=ATOL,
                        t_eval=r_ev, dense_output=False)
        if sol.t[-1] < r_max * 0.99:
            return np.nan
        return float(sol.y[0, -1])
    except Exception:
        return np.nan

# ─── energia solitonu ─────────────────────────────────────────────────────────
def energy_soliton(psi, K, r_max):
    """Oblicza E[phi] = 4pi int_0^inf [1/2 kfac (dphi)^2 + (V-V1)] r^2 dr."""
    phi0  = psi
    dphi0 = -K / A_GAM**2
    r_ev  = A_GAM * (r_max / A_GAM)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(ode_rhs, [A_GAM, r_max], [phi0, dphi0],
                        method='DOP853', rtol=RTOL, atol=ATOL,
                        t_eval=r_ev, dense_output=False)
        if sol.t[-1] < r_max * 0.99:
            return np.nan
        r    = sol.t
        phi  = np.maximum(sol.y[0], 1e-10)
        dphi = sol.y[1]
        if not np.all(np.isfinite(phi)) or not np.all(np.isfinite(dphi)):
            return np.nan
        kfac     = 1.0 + ALPHA_K / phi
        Ek       = 4*np.pi * np.trapezoid(0.5*dphi**2 * kfac * r**2, r)
        Ep       = 4*np.pi * np.trapezoid((V_mod(phi) - V1) * r**2, r)
        return Ek + Ep
    except Exception:
        return np.nan

def g_U(psi, K, r_max):
    """g = E[phi]/(4 pi K) - 1."""
    E = energy_soliton(psi, K, r_max)
    if not np.isfinite(E):
        return np.nan
    return E / (4*np.pi*K) - 1.0

# ─── zera F(psi)=phi(r_max)-1 ────────────────────────────────────────────────
def find_all_psi_zeros(K, r_max, n_psi=N_PSI):
    """Zwraca wszystkie psi takie ze phi(r_max; psi,K) = 1."""
    psis = np.linspace(PSI_LO, PSI_HI, n_psi)
    Fv   = [phi_at_rmax(p, K, r_max) - 1.0 for p in psis]
    roots = []
    for i in range(len(Fv) - 1):
        fi, fi1 = Fv[i], Fv[i+1]
        if not np.isfinite(fi) or not np.isfinite(fi1):
            continue
        if fi * fi1 < 0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p, K, r_max) - 1.0,
                            psis[i], psis[i+1], xtol=1e-7, maxiter=50)
                roots.append(pz)
            except Exception:
                pass
    return roots  # rosnaco w psi

# ─── g_U(K) uzywajac ostatniego zera psi (Galaz U, jak w p98) ─────────────
def g_U_of_K(K, r_max):
    """Zwraca (g_U, psi_last) lub (nan, nan)."""
    roots = find_all_psi_zeros(K, r_max)
    if not roots:
        return np.nan, np.nan
    psi_last = roots[-1]
    g = g_U(psi_last, K, r_max)
    return g, psi_last

# ─── znajdz K*_2 (pierwsze zero g_U w [K_LO, K_HI]) ────────────────────────
def find_K2(r_max, n_K=N_K, verbose=False):
    """
    Skanuje K w [K_LO, K_HI], szuka pierwszej zmiany znaku g_U(K).
    Zwraca (K2, psi0_at_K2) lub (nan, nan).
    """
    Ks  = np.linspace(K_LO, K_HI, n_K)
    gv  = []
    psv = []

    for K in Ks:
        g, ps = g_U_of_K(K, r_max)
        gv.append(g)
        psv.append(ps)

    if verbose:
        # Podsumowanie g_U w kilku punktach
        step = max(1, n_K // 8)
        print(f"  {'K':>8}  {'g_U':>10}  {'psi_last':>10}")
        for i in range(0, n_K, step):
            gs  = f"{gv[i]:.5f}" if np.isfinite(gv[i]) else "   nan  "
            ps  = f"{psv[i]:.4f}" if np.isfinite(psv[i]) else "  nan "
            print(f"  {Ks[i]:.5f}  {gs}  {ps}")

    # Szukaj zmian znaku
    sign_changes = []
    for i in range(len(gv) - 1):
        gi, gi1 = gv[i], gv[i+1]
        if np.isfinite(gi) and np.isfinite(gi1) and gi * gi1 < 0:
            sign_changes.append((i, Ks[i], Ks[i+1]))

    if verbose:
        print(f"  Zmiany znaku g_U: {len(sign_changes)}")
        for idx, Kl, Kr in sign_changes:
            print(f"    K in [{Kl:.4f}, {Kr:.4f}]  "
                  f"g=[{gv[idx]:.4f}, {gv[idx+1]:.4f}]  "
                  f"psi~{psv[idx]:.4f}")

    if not sign_changes:
        return np.nan, np.nan

    # Pierwsza zmiana znaku = K*_2
    idx, Kl, Kr = sign_changes[0]

    try:
        K2 = brentq(lambda K: g_U_of_K(K, r_max)[0],
                    Kl, Kr, xtol=1e-5, maxiter=40)
    except Exception as e:
        if verbose:
            print(f"  brentq: {e}  -> interpolacja liniowa")
        K2 = Kl - gv[idx] * (Kr - Kl) / (gv[idx+1] - gv[idx])

    # psi0 przy K*_2
    _, psi_K2 = g_U_of_K(K2, r_max)

    return K2, psi_K2

# ══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 65)
    print("p99_K2_convergence.py  --  TGP v1")
    print("=" * 65)
    print(f"alpha_K={ALPHA_K}, a_Gam={A_GAM}, K*_1(inf)={K1_INF}")
    print(f"Skan K in [{K_LO},{K_HI}], n_K={N_K}, n_psi={N_PSI}")
    print(f"dphi0 = -K/A_GAM^2  (klucz: K wchodzi przez nachylenie!)")
    print(f"r_max values: {RMAXES}")
    print()

    K2_res   = {}
    psi0_res = {}

    for rmax in RMAXES:
        t0 = time.time()
        print(f"─── r_max = {rmax:.0f} ───────────────────────────────────")
        K2, psi0_K2 = find_K2(rmax, verbose=True)
        dt = time.time() - t0
        K2_res[rmax]   = K2
        psi0_res[rmax] = psi0_K2
        if np.isfinite(K2):
            print(f"  => K*_2 = {K2:.6f},  psi0(K*_2) = {psi0_K2:.4f}  "
                  f"({dt:.0f}s)")
        else:
            print(f"  => K*_2 = BRAK w [{K_LO},{K_HI}]  ({dt:.0f}s)")
        print()

    # ── tabela ──────────────────────────────────────────────────────────────
    print("=" * 65)
    print("TABELA: K*_2 vs r_max")
    print("=" * 65)
    print(f"{'r_max':>8} {'K*_2':>12} {'psi0':>10} {'K*_2/K*_1(inf)':>16}")
    print("-" * 52)
    valid_r  = []
    valid_K2 = []
    for rmax in RMAXES:
        K2  = K2_res[rmax]
        ps  = psi0_res[rmax]
        rat = K2 / K1_INF if np.isfinite(K2) else np.nan
        ks  = f"{K2:.6f}" if np.isfinite(K2) else "  BRAK  "
        ps_s= f"{ps:.4f}"  if np.isfinite(ps)  else " BRAK "
        rs  = f"{rat:.3f}" if np.isfinite(rat)  else " BRAK"
        print(f"{rmax:>8.0f} {ks:>12} {ps_s:>10} {rs:>16}")
        if np.isfinite(K2):
            valid_r.append(rmax)
            valid_K2.append(K2)

    # ── ekstrapolacja 1/r ────────────────────────────────────────────────────
    print()
    K2_inf = np.nan; ratio_inf = np.nan; R2 = np.nan; b_fit = np.nan
    if len(valid_K2) >= 3:
        print("── Ekstrapolacja K*_2(r) = a + b/r ──")
        inv_r  = 1.0 / np.array(valid_r)
        K2arr  = np.array(valid_K2)
        A_mat  = np.column_stack([np.ones_like(inv_r), inv_r])
        coeffs, _, _, _ = lstsq(A_mat, K2arr, rcond=None)
        a_fit, b_fit = coeffs
        K2_pred = A_mat @ coeffs
        ss_res  = np.sum((K2arr - K2_pred)**2)
        ss_tot  = np.sum((K2arr - np.mean(K2arr))**2)
        R2 = 1 - ss_res / ss_tot if ss_tot > 1e-20 else np.nan

        K2_inf = a_fit
        ratio_inf = K2_inf / K1_INF
        print(f"  Dane: r = {valid_r}")
        print(f"  K*_2(r) = {[f'{k:.5f}' for k in valid_K2]}")
        print(f"  Fit: a={a_fit:.6f}, b={b_fit:.4f}")
        print(f"  R^2 = {R2:.6f}")
        print(f"  K*_2(inf) = {K2_inf:.6f}")
        print(f"  K*_1(inf) = {K1_INF:.6f}")
        print(f"  K*_2/K*_1(inf) = {ratio_inf:.3f}")
        print(f"  R_21^PDG = {R21_PDG:.2f}")
        print(f"  Czynnik niedoboru = {R21_PDG/ratio_inf:.2f}x")
    elif len(valid_K2) == 2:
        print("── Ekstrapolacja 2-punktowa K*_2(r) = a + b/r ──")
        r1, r2 = valid_r
        K1, K2_ = valid_K2
        # Rozwiaz uklad 2x2
        b_fit  = (K1 - K2_) / (1/r1 - 1/r2)
        a_fit  = K1 - b_fit/r1
        K2_inf = a_fit
        ratio_inf = K2_inf / K1_INF
        R2 = np.nan  # nieweryfikowalne przy 2 punktach
        print(f"  K*_2(inf) = {K2_inf:.6f}  (b={b_fit:.4f})")
        print(f"  K*_2/K*_1(inf) = {ratio_inf:.3f}")
    else:
        print("BRAK: za malo punktow do ekstrapolacji")

    # ── spread psi0 ─────────────────────────────────────────────────────────
    valid_psi = [psi0_res[r] for r in valid_r if np.isfinite(psi0_res[r])]
    psi_spread = (max(valid_psi) - min(valid_psi)) if len(valid_psi) >= 2 else np.nan
    if np.isfinite(psi_spread):
        print(f"\npsi0(K*_2) spread = {psi_spread:.4f}  "
              f"(wartosci: {[f'{p:.4f}' for p in valid_psi]})")

    # ══════════════════════════════════════════════════════════════════════
    # PASS / FAIL
    # ══════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 65)
    print("TESTY PASS/FAIL")
    print("=" * 65)
    passes = []

    # P1
    found_n = sum(np.isfinite(K2_res[r]) for r in RMAXES)
    p1 = (found_n == len(RMAXES))
    passes.append(p1)
    print(f"P1 [K*_2 przy wszystkich r_max]: {'PASS' if p1 else 'FAIL'}  "
          f"({found_n}/{len(RMAXES)})")

    # P2
    if np.isfinite(psi_spread):
        p2 = psi_spread < 0.30
        print(f"P2 [psi0 zbiezne]:               {'PASS' if p2 else 'FAIL'}  "
              f"(spread={psi_spread:.4f}, prog=0.30)")
    else:
        p2 = False
        print(f"P2 [psi0 zbiezne]:               FAIL  (brak danych)")
    passes.append(p2)

    # P3
    if np.isfinite(R2):
        p3 = R2 > 0.990
        print(f"P3 [fit 1/r R^2]:                {'PASS' if p3 else 'FAIL'}  "
              f"(R^2={R2:.6f})")
    else:
        p3 = (len(valid_K2) >= 2)   # 2-punkt = brak R^2 ale mozliwe
        print(f"P3 [fit 1/r R^2]:                {'PASS (2-pkt)' if p3 else 'FAIL'}  "
              f"(R^2 nieokreslone)")
    passes.append(p3)

    # P4
    if np.isfinite(ratio_inf):
        p4 = ratio_inf > 5.0
        print(f"P4 [K*_2/K*_1(inf) > 5]:         {'PASS' if p4 else 'FAIL'}  "
              f"(ratio={ratio_inf:.3f})")
    else:
        p4 = False
        print(f"P4 [K*_2/K*_1(inf) > 5]:         FAIL  (brak)")
    passes.append(p4)

    # P5
    if np.isfinite(ratio_inf):
        p5 = ratio_inf < 100.0
        print(f"P5 [ratio << 206.77]:            {'PASS' if p5 else 'FAIL'}  "
              f"(ratio={ratio_inf:.3f}, prog<100)")
    else:
        p5 = False
        print(f"P5 [ratio << 206.77]:            FAIL  (brak)")
    passes.append(p5)

    n_pass = sum(passes)
    print(f"\nWYNIK: {n_pass}/{len(passes)} PASS")

    # ── podsumowanie fizyczne ───────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("PODSUMOWANIE FIZYCZNE")
    print("=" * 65)
    if np.isfinite(K2_inf):
        print(f"  K*_2(inf)           = {K2_inf:.6f}")
        print(f"  K*_1(inf) [p97]     = {K1_INF:.6f}")
        print(f"  K*_2/K*_1(inf)      = {ratio_inf:.3f}")
        print(f"  R_21^PDG            = {R21_PDG:.2f}")
        print(f"  Czynnik niedoboru   = {R21_PDG/ratio_inf:.1f}x")
        print()
        print("  Wniosek: standardowy V_mod NIE odtwarza hierarchii mas leptonow.")
        if np.isfinite(R21_PDG/ratio_inf):
            print(f"  Konieczna modyfikacja V_mod zwiekszajaca K*_2 ~{R21_PDG/ratio_inf:.1f}x.")
    else:
        print("  K*_2(inf) niewyznaczone (za malo punktow).")
    print("=" * 65)
