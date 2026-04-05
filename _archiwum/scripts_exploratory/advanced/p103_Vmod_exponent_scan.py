#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p103_Vmod_exponent_scan.py  --  TGP v1
======================================
Sciezka 1 Koide: skan wykladnika n w V_mod = gamma/3*phi^3 - gamma/4*phi^4 + lam/n*(phi-1)^n
Cel: znalezc n dajace K*_2/K*_1 = 206.77 (stosunek mas mion/elektron).

Wniosek z p102: zmiana lambda nie zmienia K*_2/K*_1 (= 16.66 = const).
=> Konieczna zmiana formy funkcjonalnej: testujemy rozne wykladniki n.

Fizyczna intuicja:
  - K*_1 (psi_0~1.24): (0.24)^n -> 0 dla n>>1  => K*_1 niezalezne od n
  - K*_2 (psi_0~4.1 ): (3.1 )^n -> inf dla n>>1 => K*_2 rosnie gwaltownie z n
  => K*_2/K*_1 powinno rosnac z n

Metoda (jak p102):
  - K*_2: dwufazowy skan K z PSI2=[3.50,7.00]
    Faza 1: K in [K2_P102*0.25, K2_P102*4.0] = [0.043, 0.685], N_K=50
            (dokladnie jak p102; krok~0.013 gwarantuje braktowanie K*_2~0.171)
    Faza 2: K in [0.685, 3.50], N_K=50
            (dla duzych n gdzie K*_2 moze byc >> K*_2(n=6))
  - Diagnostyka naprawa z p102: sledz psi_0 wewnatrz g_func (unika re-run po brentq)

Data: 2026-03-26
"""
import sys, io, warnings, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ── stale fizyczne ────────────────────────────────────────────────────────────
ALPHA_K    = 8.5616
GAM        = 1.0
V1         = GAM/3.0 - GAM/4.0   # = 1/12 (niezalezne od n)
LAM        = 5.501357e-06          # lambda_phys
AGAM       = 0.040
R_MAX      = 80.0
R21_TARGET = 206.77

# ── wartosci referencyjne z p102 ─────────────────────────────────────────────
K1_REF  = 0.010281
K2_P102 = 0.171272

# ── okna psi_0 ────────────────────────────────────────────────────────────────
PSI1_LO, PSI1_HI = 1.10, 1.80   # K*_1: wyklucza spurious psi~1.008
PSI2_LO, PSI2_HI = 3.50, 7.00   # K*_2: jak p102 (sprawdzone)

# ── zakresy K (dwufazowe dla K*_2) ────────────────────────────────────────────
K1_LO = K1_REF * 0.70            # ~0.0072
K1_HI = K1_REF * 1.50            # ~0.0154

# Faza 1: dok. jak p102 — krok ~0.013, gwarantuje bracket K*_2(n=6)=0.171
K2_LO_P1 = K2_P102 * 0.25        # ~0.043
K2_HI_P1 = K2_P102 * 4.0         # ~0.685
N_K2_P1  = 50

# Faza 2: rozszerzona dla duzych n
K2_LO_P2 = K2_P102 * 4.0         # ~0.685
K2_HI_P2 = K1_REF  * 350.0       # ~3.60  (pokrywa K*_2_cel ~2.13)
N_K2_P2  = 50

# ── parametry numeryczne ──────────────────────────────────────────────────────
N_PSI1 = 50
N_PSI2 = 60
N_K1   = 20
N_EVAL = 600
RTOL   = 1e-8
ATOL   = 1e-10

# ── wykladniki do skanu ───────────────────────────────────────────────────────
N_VALUES = [4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20]


# ── potencjal i ODE ───────────────────────────────────────────────────────────

def dV_mod_n(phi, n):
    """dV/dphi = gamma*phi^2 - gamma*phi^3 + lam*(phi-1)^(n-1)"""
    delta = phi - 1.0
    try:
        corr = LAM * float(delta)**(n - 1)
        if not np.isfinite(corr):
            corr = np.sign(delta)**(n-1) * 1e15
    except (OverflowError, ValueError):
        corr = np.sign(delta)**(n-1) * 1e15
    return GAM * phi**2 - GAM * phi**3 + corr


def V_mod_n(phi, n):
    """V_mod = gamma/3*phi^3 - gamma/4*phi^4 + lam/n*(phi-1)^n"""
    delta = phi - 1.0
    try:
        stab = (LAM / n) * float(delta)**n
        if not np.isfinite(stab):
            stab = 1e15
    except (OverflowError, ValueError):
        stab = 1e15
    return GAM/3.0 * phi**3 - GAM/4.0 * phi**4 + stab


def ode_rhs_n(r, y, n):
    phi, dphi = y
    phi  = max(phi, 1e-10)
    kfac = 1.0 + ALPHA_K / phi
    phi2 = phi**2
    return [dphi,
            dV_mod_n(phi, n) / kfac
            + ALPHA_K * dphi**2 / (2.0 * phi2 * kfac)
            - 2.0 / r * dphi]


def phi_at_rmax(psi, K, n, r_max=R_MAX):
    dphi0 = -K / AGAM**2
    r_ev  = AGAM * (r_max / AGAM)**np.linspace(0, 1, N_EVAL)
    for method in ('DOP853', 'Radau'):
        try:
            sol = solve_ivp(lambda r, y: ode_rhs_n(r, y, n),
                            [AGAM, r_max], [psi, dphi0],
                            method=method, rtol=RTOL, atol=ATOL,
                            t_eval=r_ev, dense_output=False)
            if sol.t[-1] < r_max * 0.99:
                if method == 'DOP853':
                    continue
                return np.nan
            return float(sol.y[0, -1])
        except Exception:
            if method == 'DOP853':
                continue
            return np.nan
    return np.nan


def energy_soliton(psi, K, n, r_max=R_MAX):
    dphi0 = -K / AGAM**2
    r_ev  = AGAM * (r_max / AGAM)**np.linspace(0, 1, N_EVAL)
    for method in ('DOP853', 'Radau'):
        try:
            sol = solve_ivp(lambda r, y: ode_rhs_n(r, y, n),
                            [AGAM, r_max], [psi, dphi0],
                            method=method, rtol=RTOL, atol=ATOL,
                            t_eval=r_ev, dense_output=False)
            if sol.t[-1] < r_max * 0.99:
                if method == 'DOP853':
                    continue
                return np.nan
            r_arr   = sol.t
            phi_arr = np.maximum(sol.y[0], 1e-10)
            dphi_arr = sol.y[1]
            if not (np.all(np.isfinite(phi_arr)) and np.all(np.isfinite(dphi_arr))):
                if method == 'DOP853':
                    continue
                return np.nan
            kfac = 1.0 + ALPHA_K / phi_arr
            Ek   = 4*np.pi * np.trapezoid(0.5 * dphi_arr**2 * kfac * r_arr**2, r_arr)
            Vpot = np.clip(np.array([V_mod_n(float(p), n) for p in phi_arr]),
                           -1e12, 1e12)
            Ep   = 4*np.pi * np.trapezoid((Vpot - V1) * r_arr**2, r_arr)
            E = Ek + Ep
            if not np.isfinite(E):
                if method == 'DOP853':
                    continue
                return np.nan
            return E
        except Exception:
            if method == 'DOP853':
                continue
            return np.nan
    return np.nan


def find_psi_zero(K, n, psi_lo, psi_hi, n_psi):
    """Pierwsze zero phi(r_max)=1 w [psi_lo, psi_hi]. Zwraca psi_zero lub nan."""
    psis = np.linspace(psi_lo, psi_hi, n_psi)
    Fv   = [phi_at_rmax(p, K, n) - 1.0 for p in psis]
    for i in range(len(Fv) - 1):
        fi, fi1 = Fv[i], Fv[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1)):
            continue
        if fi * fi1 < 0:
            try:
                return brentq(lambda p: phi_at_rmax(p, K, n) - 1.0,
                              psis[i], psis[i+1], xtol=1e-7, maxiter=50)
            except Exception:
                pass
    return np.nan


def g_at_psi(psi, K, n):
    E = energy_soliton(psi, K, n)
    return E / (4.0*np.pi*K) - 1.0 if (np.isfinite(E) and K > 0) else np.nan


def g_U_narrow(K, n, psi_lo, psi_hi, n_psi):
    """Oblicz g_U i psi_0 dla danego K w oknie psi."""
    pz = find_psi_zero(K, n, psi_lo, psi_hi, n_psi)
    if np.isnan(pz):
        return np.nan, np.nan
    return g_at_psi(pz, K, n), pz


def find_Kstar_range(n, K_lo, K_hi, n_K,
                     psi_lo=PSI2_LO, psi_hi=PSI2_HI, n_psi=N_PSI2,
                     verbose=False):
    """
    Znajdz K* (K*_2) w [K_lo, K_hi].
    Metodologia jak p102: skan g_U, szukaj zmiany znaku, brentq.
    Sledzi psi_0 wewnatrz g_func (fix: poprawne g_star bez re-run).
    """
    Ks  = np.linspace(K_lo, K_hi, n_K)
    gv  = []
    psv = []
    for K in Ks:
        g, pz = g_U_narrow(K, n, psi_lo, psi_hi, n_psi)
        gv.append(g)
        psv.append(pz)

    # Zbierz zmiany znaku
    sign_changes = []
    for i in range(len(gv) - 1):
        gi, gi1 = gv[i], gv[i+1]
        if np.isfinite(gi) and np.isfinite(gi1) and gi * gi1 < 0:
            sign_changes.append((Ks[i], Ks[i+1], gi, gi1, psv[i]))

    if verbose:
        n_valid = sum(np.isfinite(g) for g in gv)
        print(f"        g_U valid: {n_valid}/{len(gv)}, zmian znaku: {len(sign_changes)}")
        valid = [(Ks[i], gv[i], psv[i]) for i in range(len(gv)) if np.isfinite(gv[i])]
        for K_v, g_v, pz_v in valid[:6]:
            pz_s = f"{pz_v:.4f}" if np.isfinite(pz_v) else "nan"
            print(f"          K={K_v:.5f} g={g_v:.4f} psi={pz_s}")
        if len(valid) > 6:
            print(f"          ... ({len(valid)-6} wiecej)")
        if sign_changes:
            Kl, Kr, gl, gr, pz0 = sign_changes[0]
            print(f"        1. zmiana znaku: K=[{Kl:.5f},{Kr:.5f}], g=[{gl:.4f},{gr:.4f}]")

    if not sign_changes:
        return np.nan, np.nan, np.nan

    # Bierz pierwsza zmiane znaku
    Kl, Kr, gl, gr, _ = sign_changes[0]

    # Naprawa: sledz ostatnie psi_0 w g_func (poprawny g_star)
    _last_psi = [np.nan]

    def g_func(K):
        g, pz = g_U_narrow(K, n, psi_lo, psi_hi, n_psi)
        _last_psi[0] = pz
        return g if np.isfinite(g) else np.nan

    try:
        # Weryfikacja znaku na kranicach przed brentq
        ga = g_func(Kl); gb = g_func(Kr)
        if not (np.isfinite(ga) and np.isfinite(gb) and ga * gb < 0):
            # Uzyj oryginalnych wartosci z coarse scan
            Kstar = Kl - gl * (Kr - Kl) / (gr - gl)
            psi_star = find_psi_zero(Kstar, n, psi_lo, psi_hi, n_psi)
            E_star   = energy_soliton(psi_star, Kstar, n)
            g_star   = E_star / (4*np.pi*Kstar) - 1.0 if np.isfinite(E_star) else np.nan
            return Kstar, psi_star, g_star

        Kstar = brentq(g_func, Kl, Kr, xtol=1e-5, maxiter=40)
        psi_star = _last_psi[0]       # psi_0 z ostatniej ewaluacji g_func (przy Kstar)
        if not np.isfinite(psi_star):
            psi_star = find_psi_zero(Kstar, n, psi_lo, psi_hi, n_psi)
        E_star = energy_soliton(psi_star, Kstar, n)
        g_star = E_star / (4*np.pi*Kstar) - 1.0 if np.isfinite(E_star) else np.nan
    except Exception:
        Kstar = Kl - gl * (Kr - Kl) / (gr - gl)
        psi_star = find_psi_zero(Kstar, n, psi_lo, psi_hi, n_psi)
        E_star   = energy_soliton(psi_star, Kstar, n)
        g_star   = E_star / (4*np.pi*Kstar) - 1.0 if np.isfinite(E_star) else np.nan

    return Kstar, psi_star, g_star


def find_K2star(n, verbose=False):
    """
    Dwufazowe szukanie K*_2:
    Faza 1: K in [K2_LO_P1, K2_HI_P1] — dokladnie jak p102 (krok~0.013)
    Faza 2: K in [K2_LO_P2, K2_HI_P2] — rozszerzone dla duzych n
    """
    if verbose:
        print(f"      Faza 1: K=[{K2_LO_P1:.3f},{K2_HI_P1:.3f}], N_K={N_K2_P1}...", flush=True)
    K2, psi2, g2 = find_Kstar_range(n, K2_LO_P1, K2_HI_P1, N_K2_P1,
                                     PSI2_LO, PSI2_HI, N_PSI2, verbose=verbose)
    if np.isfinite(K2):
        return K2, psi2, g2, 1

    if verbose:
        print(f"      Faza 2: K=[{K2_LO_P2:.3f},{K2_HI_P2:.3f}], N_K={N_K2_P2}...", flush=True)
    K2, psi2, g2 = find_Kstar_range(n, K2_LO_P2, K2_HI_P2, N_K2_P2,
                                     PSI2_LO, PSI2_HI, N_PSI2, verbose=verbose)
    if np.isfinite(K2):
        return K2, psi2, g2, 2

    return np.nan, np.nan, np.nan, 0


# ── GLOWNA PETLA ──────────────────────────────────────────────────────────────
t0 = time.time()
print("=" * 72)
print("p103: Skan wykladnika n  --  V_mod = g/3*phi^3 - g/4*phi^4 + lam/n*(phi-1)^n")
print(f"lambda={LAM:.3e}, a_Gamma={AGAM}, r_max={R_MAX}")
print(f"Cel: K*_2/K*_1 = {R21_TARGET}  (K*_2_cel ~ {R21_TARGET*K1_REF:.4f})")
print(f"K*_2 Faza1: [{K2_LO_P1:.3f},{K2_HI_P1:.3f}] N={N_K2_P1} krok~{(K2_HI_P1-K2_LO_P1)/(N_K2_P1-1):.4f}")
print(f"K*_2 Faza2: [{K2_LO_P2:.3f},{K2_HI_P2:.3f}] N={N_K2_P2} krok~{(K2_HI_P2-K2_LO_P2)/(N_K2_P2-1):.4f}")
print("=" * 72)

results = []

for n in N_VALUES:
    tn = time.time()
    print(f"\n=== n = {n} ===", flush=True)

    # ── K*_1 ──────────────────────────────────────────────────────────────────
    print(f"  K*_1: K=[{K1_LO:.4f},{K1_HI:.4f}] x psi=[{PSI1_LO},{PSI1_HI}]...",
          end=' ', flush=True)
    K1, psi1, g1 = find_Kstar_range(n, K1_LO, K1_HI, N_K1,
                                     PSI1_LO, PSI1_HI, N_PSI1)
    if np.isfinite(K1):
        print(f"K*_1={K1:.6f}  psi_0={psi1:.4f}  g={g1:.2e}")
    else:
        print("BRAK")

    # ── K*_2 ──────────────────────────────────────────────────────────────────
    print(f"  K*_2: ...", end=' ', flush=True)
    K2, psi2, g2, phase = find_K2star(n, verbose=False)
    if np.isfinite(K2):
        print(f"K*_2={K2:.6f}  psi_0={psi2:.4f}  g={g2:.2e}  [Faza{phase}]")
    else:
        print("BRAK (obie fazy)")

    # ── stosunek ──────────────────────────────────────────────────────────────
    if np.isfinite(K1) and np.isfinite(K2):
        r21 = K2 / K1
        frac = r21 / R21_TARGET
        tag = "  <== CEL!" if r21 >= R21_TARGET else f"  ({frac:.1%} celu)"
        print(f"  K*_2/K*_1 = {r21:.4f}{tag}")
    else:
        r21 = np.nan

    dt = time.time() - tn
    print(f"  [czas n={n}: {dt:.0f}s]", flush=True)
    results.append((n, K1, psi1, g1, K2, psi2, g2, r21))


# ── TABELA ZBIORCZA ───────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("TABELA WYNIKOW — skan wykladnika n")
print(f"  {'n':>4}  {'K*_1':>10}  {'psi0_1':>7}  {'K*_2':>10}  {'psi0_2':>7}  {'r21':>8}  {'r21/cel':>8}")
print("  " + "-" * 63)
for n, K1, psi1, g1, K2, psi2, g2, r21 in results:
    K1s  = f"{K1:.6f}"   if np.isfinite(K1)   else "    --    "
    K2s  = f"{K2:.6f}"   if np.isfinite(K2)   else "    --    "
    p1s  = f"{psi1:.4f}" if np.isfinite(psi1) else "   --  "
    p2s  = f"{psi2:.4f}" if np.isfinite(psi2) else "   --  "
    r21s = f"{r21:.3f}"  if np.isfinite(r21)  else "   --  "
    frcs = f"{r21/R21_TARGET:.4f}" if np.isfinite(r21) else "   --  "
    print(f"  {n:>4}  {K1s:>10}  {p1s:>7}  {K2s:>10}  {p2s:>7}  {r21s:>8}  {frcs:>8}")


# ── PODSUMOWANIE FIZYCZNE ─────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("PODSUMOWANIE FIZYCZNE  --  Sciezka 1 Koide:")
print(f"lambda = {LAM:.3e} (stale).  Cel: K*_2/K*_1 = {R21_TARGET}")
print(f"Referencja n=6 (z p102): K*_2/K*_1 ~ 16.66")
print()

valid = [(n, r) for n, K1, p1, g1, K2, p2, g2, r in results if np.isfinite(r)]
if valid:
    r_vals = [r for _, r in valid]
    n_vals = [n for n, _ in valid]
    n_max, r21_max = max(valid, key=lambda x: x[1])
    n_min, r21_min = min(valid, key=lambda x: x[1])

    print(f"n z oboma K*: {n_vals}")
    print(f"K*_2/K*_1 zakres: [{r21_min:.2f}, {r21_max:.2f}]")
    print(f"Max K*_2/K*_1 = {r21_max:.2f} przy n = {n_max}")

    if len(r_vals) >= 2:
        diffs = [r_vals[i+1] - r_vals[i] for i in range(len(r_vals)-1)]
        mono  = all(d > -2.0 for d in diffs)
        print(f"Trend: {'monoton. rosn.' if mono else 'NIE monot.'}")

    if r21_max >= R21_TARGET:
        for n_i, r_i in valid:
            if r_i >= R21_TARGET:
                print(f"\n==> CEL OSIAGNIETY: n={n_i} daje K*_2/K*_1={r_i:.2f}")
                break
        # Interpolacja
        for i in range(len(valid)-1):
            if valid[i][1] < R21_TARGET <= valid[i+1][1]:
                n_lo, r_lo = valid[i]; n_hi, r_hi = valid[i+1]
                frac_n = (R21_TARGET - r_lo) / (r_hi - r_lo)
                n_interp = n_lo + frac_n * (n_hi - n_lo)
                print(f"    Interpolacja n*: {n_interp:.2f}  (miedzy n={n_lo} i n={n_hi})")
    else:
        deficit = R21_TARGET / r21_max
        print(f"\n==> Cel NIEOSIAGNIETY. Max r21={r21_max:.2f} (n={n_max}), brakuje x{deficit:.1f}")
        print(f"    Trend: {' '.join(f'n={n}:{r:.1f}' for n,r in valid)}")
else:
    print("BRAK poprawnych wynikow!")


# ── PASS / FAIL ───────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("PASS/FAIL:")
K1_n6  = next((K1 for n, K1, *_ in results if n == 6 and np.isfinite(K1)), np.nan)
r21_valid = [(n, r) for n, *_, r in results if np.isfinite(r)]
r_all  = [r for _, r in r21_valid]

P1 = np.isfinite(K1_n6)
P2 = P1 and abs(K1_n6 / K1_REF - 1.0) < 0.05
P3 = len(r21_valid) >= 5
P4 = len(r_all) >= 3 and all(r_all[i+1]-r_all[i] > -2.0 for i in range(len(r_all)-1))
P5 = any(r >= 50.0    for _, r in r21_valid)
P6 = any(r >= R21_TARGET for _, r in r21_valid)

for i, (P, desc) in enumerate([
    (P1, "K*_1 znaleziony dla n=6 (cross-check z p102)"),
    (P2, f"K*_1(n=6) stabilny: |K/K1_REF-1|<5%  (K1_REF={K1_REF})"),
    (P3, "K*_2 znaleziony dla >=5 wartosci n"),
    (P4, "K*_2/K*_1 monotonicznie rosnie z n (tolerancja 2)"),
    (P5, "max K*_2/K*_1 >= 50  (postep ku celowi)"),
    (P6, f"K*_2/K*_1 = {R21_TARGET} osiagniety (CEL KOIDE)"),
], 1):
    print(f"  P{i}: {'PASS' if P else 'FAIL'} -- {desc}")

n_pass = sum([P1,P2,P3,P4,P5,P6])
total  = time.time() - t0
print(f"\nWynik: {n_pass}/6 PASS  (czas calkowity: {total/60:.1f} min)")
