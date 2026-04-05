#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p106_alphaSN_vs_agamma.py  --  TGP v1
======================================
Test hipotezy: alpha_SN(a_Gamma) = A - B*a_Gamma
  => czy A = 9.0 dokładnie (w granicy a_Gamma -> 0)?

Z p105_v2: alpha_SN(a_Gamma = 0.040) = 8.9964 +/- 0.0004

Dla każdego a_Gamma w AGAM_VALS:
  - Bisekcja w alpha_K in [AK_LO, AK_HI]
  - Szukamy alpha_SN(a_Gamma) z precyzją BISECT_TOL
  - Fit liniowy: alpha_SN = A - B*a_Gamma -> ekstrapolacja do a_Gamma=0 -> A=?

Data: 2026-03-26
"""
import sys
import io
import warnings
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from concurrent.futures import ThreadPoolExecutor, as_completed

warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ── Stałe fizyczne ─────────────────────────────────────────────────────────────
GAM    = 1.0
V1     = GAM/3.0 - GAM/4.0
LAM    = 5.501357e-06
R_MAX  = 80.0
N_EVAL = 500
RTOL   = 1e-8
ATOL   = 1e-10

# ── Parametry skanowania ───────────────────────────────────────────────────────
AGAM_REF = 0.040
AK_REF   = 8.5616
AKSN_REF = 8.9964    # z p105_v2
K2_AKSN_REF = 0.012771   # K*₂ przy alpha_SN dla AGAM_REF
PSI2_AKSN_REF = 1.2549   # psi₂ przy alpha_SN dla AGAM_REF

# Wartości a_Gamma do sprawdzenia
AGAM_VALS = [0.020, 0.030, 0.040, 0.050, 0.060]

# Bisekcja: bracket [AK_LO, AK_HI] musi zawierać alpha_SN dla każdego a_Gamma
# Z hipotezy: alpha_SN ~ 9 - 0.09*agam:
#   agam=0.020: ~8.998   agam=0.060: ~8.995
# Bezpieczny bracket: [8.90, 9.10]
AK_LO = 8.90
AK_HI = 9.10
BISECT_TOL = 0.002   # wystarczająca precyzja do testu trendu

# Okno adaptacyjne wokół K*₂ i psi₂
K_FRAC   = 3.0    # okno K: [K_prev/K_FRAC, K_prev*K_FRAC]
PSI_DELT = 0.50   # okno psi: psi_prev +/- PSI_DELT
PSI_LO_ABS = 1.05  # absolutne minimum psi

N_K      = 25
N_PSI    = 80
N_WORKERS = 16


# ── Funkcje ODE ───────────────────────────────────────────────────────────────

def V_mod(phi):
    return GAM/3.0*phi**3 - GAM/4.0*phi**4 + LAM/6.0*(phi - 1.0)**6


def make_rhs(alpha_K):
    def rhs(r, y):
        phi, dphi = y
        phi  = max(phi, 1e-10)
        kfac = 1.0 + alpha_K / phi
        dV   = GAM*phi**2 - GAM*phi**3 + LAM*(phi - 1.0)**5
        return [dphi,
                dV/kfac + alpha_K*dphi**2 / (2.0*phi**2*kfac) - 2.0/r*dphi]
    return rhs


def solve_full(psi, K, rhs_fn, agam):
    dphi0 = -K / agam**2
    r_ev  = agam * (R_MAX/agam)**np.linspace(0, 1, N_EVAL)
    for method in ('DOP853', 'Radau'):
        try:
            sol = solve_ivp(rhs_fn, [agam, R_MAX], [psi, dphi0],
                            method=method, rtol=RTOL, atol=ATOL, t_eval=r_ev)
            if sol.t[-1] < R_MAX * 0.99:
                if method == 'DOP853': continue
                return None
            return sol
        except Exception:
            if method == 'DOP853': continue
            return None
    return None


def phi_end(psi, K, rhs_fn, agam):
    sol = solve_full(psi, K, rhs_fn, agam)
    return float(sol.y[0, -1]) if sol else np.nan


def is_mono(sol, tol=0.01):
    return bool(np.all(np.diff(sol.y[0]) <= tol))


def energy_sol(sol, alpha_K):
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    if not (np.all(np.isfinite(phi)) and np.all(np.isfinite(dphi))):
        return np.nan
    kfac = 1.0 + alpha_K / phi
    Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
    Vpot = np.clip([V_mod(float(p)) for p in phi], -1e12, 1e12)
    Ep   = 4*np.pi*np.trapezoid((np.array(Vpot) - V1)*r**2, r)
    E    = Ek + Ep
    return E if np.isfinite(E) else np.nan


# ── Równoległa g_U(K) ─────────────────────────────────────────────────────────

def compute_gU_for_K(args):
    K, alpha_K, psi_lo, psi_hi, n_psi, rhs_fn, agam = args
    warnings.filterwarnings('ignore')
    psis = np.linspace(psi_lo, psi_hi, n_psi)
    Fv   = [phi_end(p, K, rhs_fn, agam) - 1.0 for p in psis]
    for i in range(len(Fv) - 1):
        fi, fi1 = Fv[i], Fv[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1) and fi*fi1 < 0):
            continue
        try:
            pz  = brentq(lambda p: phi_end(p, K, rhs_fn, agam) - 1.0,
                         psis[i], psis[i+1], xtol=1e-8, maxiter=60)
            sol = solve_full(pz, K, rhs_fn, agam)
            if sol and is_mono(sol):
                E = energy_sol(sol, alpha_K)
                g = E/(4*np.pi*K) - 1.0 if np.isfinite(E) else np.nan
                return K, pz, g, E
        except Exception:
            pass
    return K, np.nan, np.nan, np.nan


def find_K2star(alpha_K, agam, K2_prev, psi2_prev, rhs_fn):
    """Szuka K*₂ (gałąź B) przy alpha_K i agam. Adaptacyjne okno."""
    K_lo   = max(K2_prev / K_FRAC, 1e-4)
    K_hi   = min(K2_prev * K_FRAC, 5.0)
    psi_lo = max(psi2_prev - PSI_DELT, PSI_LO_ABS)
    psi_hi = min(psi2_prev + PSI_DELT, 4.0)

    Ks   = np.logspace(np.log10(K_lo), np.log10(K_hi), N_K)
    args = [(K, alpha_K, psi_lo, psi_hi, N_PSI, rhs_fn, agam) for K in Ks]

    gv   = [np.nan] * len(Ks)
    psv  = [np.nan] * len(Ks)

    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = {ex.submit(compute_gU_for_K, a): i for i, a in enumerate(args)}
        for fut in as_completed(futures):
            i = futures[fut]
            try:
                K_out, pz, g, E = fut.result()
                gv[i]  = g
                psv[i] = pz
            except Exception:
                pass

    for i in range(len(gv) - 1):
        gi, gi1 = gv[i], gv[i+1]
        if not (np.isfinite(gi) and np.isfinite(gi1) and gi*gi1 < 0):
            continue
        _last = [np.nan, np.nan, np.nan]

        def gfunc(K, _l=_last, _pl=psi_lo, _ph=psi_hi,
                  _r=rhs_fn, _ak=alpha_K, _ag=agam):
            _, pz2, g2, E2 = compute_gU_for_K((K, _ak, _pl, _ph, N_PSI, _r, _ag))
            _l[0] = pz2; _l[1] = g2; _l[2] = E2
            return g2 if np.isfinite(g2) else np.nan

        try:
            ga = gfunc(Ks[i]); gb = gfunc(Ks[i+1])
            if np.isfinite(ga) and np.isfinite(gb) and ga*gb < 0:
                Kstar = brentq(gfunc, Ks[i], Ks[i+1], xtol=1e-6, maxiter=50)
                return Kstar, _last[0] if np.isfinite(_last[0]) else psi2_prev, _last[1]
        except Exception:
            Kstar = Ks[i] - gv[i]*(Ks[i+1]-Ks[i]) / (gv[i+1]-gv[i])
            _, ps, g2, _ = compute_gU_for_K((Kstar, alpha_K, psi_lo, psi_hi, N_PSI, rhs_fn, agam))
            return Kstar, ps if np.isfinite(ps) else psi2_prev, g2

    return np.nan, np.nan, np.nan


def test_branch_B(alpha_K, agam, K2_seed, psi2_seed):
    """Zwraca (alive, K2, psi2, g, dt) dla danego (alpha_K, agam)."""
    t0     = time.time()
    rhs_fn = make_rhs(alpha_K)
    K2_new, psi2_new, g2 = find_K2star(alpha_K, agam, K2_seed, psi2_seed, rhs_fn)
    dt     = time.time() - t0

    if np.isfinite(K2_new) and np.isfinite(g2) and abs(g2) < 0.15:
        return True, K2_new, psi2_new, g2, dt
    return False, K2_seed, psi2_seed, np.nan, dt


# ── Bisekcja alpha_SN dla jednego agam ────────────────────────────────────────

def find_alphaSN_for_agam(agam):
    """Bisekcja alpha_SN(agam). Zwraca (lo, hi, K2_lo, psi2_lo, historia)."""
    # Skalowanie nasion K*₂: K*₂ ∝ agam (z K*₁ ∝ agam z p101)
    scale   = agam / AGAM_REF
    K2_seed = K2_AKSN_REF * scale   # ~K*₂ przy bifurkacji
    psi2_seed = PSI2_AKSN_REF       # ~niezależne od agam (jak psi₁)

    print(f"\n{'='*60}")
    print(f"a_Gamma = {agam:.3f}  "
          f"(K2_seed={K2_seed:.6f}, psi2_seed={psi2_seed:.4f})")
    print(f"Bisekcja alpha_K in [{AK_LO}, {AK_HI}], tol={BISECT_TOL}")

    # Weryfikacja nasion: sprawdź czy w AK_LO jest ALIVE i w AK_HI jest DEAD
    print(f"  Weryfikacja nasienna AK_LO={AK_LO} ...", end='', flush=True)
    alive_lo, K2_lo, psi2_lo, g_lo, dt_lo = test_branch_B(AK_LO, agam, K2_seed, psi2_seed)
    print(f"  {'ALIVE' if alive_lo else 'DEAD'}  "
          f"K2={K2_lo:.6f}  [{dt_lo:.0f}s]")

    if not alive_lo:
        # AK_LO też jest DEAD — trzeba zejść niżej. Próbuj 8.80
        print(f"  UWAGA: AK_LO={AK_LO} jest DEAD! Próba 8.80 ...", end='', flush=True)
        alive_lo, K2_lo, psi2_lo, g_lo, dt_lo = test_branch_B(8.80, agam, K2_seed, psi2_seed)
        print(f"  {'ALIVE' if alive_lo else 'DEAD'}  [{dt_lo:.0f}s]")
        if alive_lo:
            lo_bound = 8.80
        else:
            print(f"  BŁĄD: nie znaleziono ALIVE dla agam={agam:.3f}")
            return np.nan, np.nan, np.nan, np.nan, []
    else:
        lo_bound = AK_LO

    print(f"  Weryfikacja AK_HI={AK_HI} ...", end='', flush=True)
    alive_hi, _, _, _, dt_hi = test_branch_B(AK_HI, agam, K2_lo * 0.1, psi2_lo)
    print(f"  {'ALIVE' if alive_hi else 'DEAD'}  [{dt_hi:.0f}s]")

    if alive_hi:
        print(f"  UWAGA: AK_HI={AK_HI} jest ALIVE! Próba 9.20 ...", end='', flush=True)
        alive_hi2, _, _, _, dt_hi2 = test_branch_B(9.20, agam, K2_lo * 0.01, psi2_lo)
        print(f"  {'ALIVE' if alive_hi2 else 'DEAD'}  [{dt_hi2:.0f}s]")
        if not alive_hi2:
            hi_bound = 9.20
        else:
            print(f"  BŁĄD: nie znaleziono DEAD dla agam={agam:.3f}")
            return np.nan, np.nan, np.nan, np.nan, []
    else:
        hi_bound = AK_HI

    lo = lo_bound
    hi = hi_bound
    history = []
    step = 0

    print(f"  Start bisekcji: [{lo:.4f}, {hi:.4f}]")

    while (hi - lo) > BISECT_TOL:
        step += 1
        mid = round((lo + hi) / 2.0, 6)
        # Użyj aktualnych najlepszych nasion dla ALIVE
        print(f"  Krok {step}: ak={mid:.6f}  (dak={hi-lo:.5f}) ...",
              end='', flush=True)

        alive, K2_mid, psi2_mid, g_mid, dt = test_branch_B(mid, agam, K2_lo, psi2_lo)

        if alive:
            lo = mid
            K2_lo   = K2_mid
            psi2_lo = psi2_mid
            status = f"ALIVE  K2={K2_mid:.6f}  psi2={psi2_mid:.4f}  g={g_mid:.2e}"
        else:
            hi = mid
            status = "DEAD"

        print(f"  {status}  [{dt:.0f}s]")
        history.append((mid, alive, K2_mid if alive else np.nan,
                        psi2_mid if alive else np.nan, g_mid, dt))

    ak_SN = (lo + hi) / 2.0
    print(f"  -> alpha_SN({agam:.3f}) = {ak_SN:.6f} +/- {(hi-lo)/2:.6f}")
    return lo, hi, K2_lo, psi2_lo, history


# ── Główna pętla ──────────────────────────────────────────────────────────────

if __name__ == '__main__':
    t0_total = time.time()
    print("="*72)
    print("p106: alpha_SN(a_Gamma) — test hipotezy A = 9.0")
    print(f"Hipoteza: alpha_SN = A - B*a_Gamma, czy A=9?")
    print(f"Znane: alpha_SN(0.040) = {AKSN_REF}  (p105_v2)")
    print(f"a_Gamma do sprawdzenia: {AGAM_VALS}")
    print(f"Bisekcja: tol={BISECT_TOL}, bracket=[{AK_LO},{AK_HI}]")
    print(f"Watki: {N_WORKERS}  N_K={N_K}  N_PSI={N_PSI}")
    print("="*72)

    results = {}

    for agam in AGAM_VALS:
        if abs(agam - AGAM_REF) < 1e-6:
            # Wartość referencyjna — już znana z p105_v2
            lo_r, hi_r = AKSN_REF - 0.0004, AKSN_REF + 0.0004
            results[agam] = {
                'lo': lo_r, 'hi': hi_r,
                'ak_SN': AKSN_REF,
                'err': 0.0004,
                'history': [],
                'ref': True
            }
            print(f"\na_Gamma = {agam:.3f}: REFERENCYJNA  "
                  f"alpha_SN = {AKSN_REF:.6f} +/- 0.000400  (p105_v2)")
            continue

        lo_r, hi_r, K2_final, psi2_final, hist = find_alphaSN_for_agam(agam)
        if np.isfinite(lo_r):
            ak_SN = (lo_r + hi_r) / 2.0
            err   = (hi_r - lo_r) / 2.0
            results[agam] = {
                'lo': lo_r, 'hi': hi_r,
                'ak_SN': ak_SN, 'err': err,
                'history': hist,
                'ref': False
            }
        else:
            results[agam] = {
                'lo': np.nan, 'hi': np.nan,
                'ak_SN': np.nan, 'err': np.nan,
                'history': hist,
                'ref': False
            }

    # ── Tabela wyników ────────────────────────────────────────────────────────
    print(f"\n{'='*72}")
    print("TABELA WYNIKÓW: alpha_SN(a_Gamma)")
    print(f"  {'a_Gamma':>8}  {'alpha_SN':>10}  {'err':>8}  {'przedział'}")
    print("  " + "-"*55)

    agams_ok = []
    aksns_ok = []
    for agam in AGAM_VALS:
        r = results[agam]
        if np.isfinite(r['ak_SN']):
            lo_s = f"({r['lo']:.4f}"
            hi_s = f"{r['hi']:.4f})"
            ref_s = " [ref]" if r['ref'] else ""
            print(f"  {agam:>8.3f}  {r['ak_SN']:>10.6f}  "
                  f"{r['err']:>8.6f}  {lo_s}, {hi_s}{ref_s}")
            agams_ok.append(agam)
            aksns_ok.append(r['ak_SN'])
        else:
            print(f"  {agam:>8.3f}       BRAK")

    # ── Fit liniowy alpha_SN(a_Gamma) = A - B*a_Gamma ────────────────────────
    print(f"\n{'='*72}")
    print("FIT: alpha_SN(a_Gamma) = A - B*a_Gamma")

    if len(agams_ok) >= 2:
        agams_arr = np.array(agams_ok)
        aksns_arr = np.array(aksns_ok)

        c = np.polyfit(agams_arr, aksns_arr, 1)  # c[0]=slope, c[1]=intercept
        A_fit = c[1]  # intercept = alpha_SN(a_Gamma=0)
        B_fit = -c[0]  # slope (negated since alpha_SN = A - B*agam)

        aksns_pred = np.polyval(c, agams_arr)
        residuals  = aksns_arr - aksns_pred
        rms        = np.sqrt(np.mean(residuals**2))
        R2         = 1.0 - np.var(residuals) / np.var(aksns_arr) if np.var(aksns_arr) > 0 else np.nan

        print(f"  alpha_SN = {A_fit:.6f} - {B_fit:.4f} * a_Gamma")
        print(f"  Interpolacja a_Gamma=0: A = alpha_SN(0) = {A_fit:.6f}")
        print(f"  RMS residuals = {rms:.6f}")
        print(f"  R^2 = {R2:.6f}" if np.isfinite(R2) else "")

        # Test hipotezy A = 9
        dev_from_9 = abs(A_fit - 9.0)
        dev_pct    = dev_from_9 / 9.0 * 100
        print(f"\n  Odchylenie A od 9.0: {dev_from_9:.6f}  ({dev_pct:.4f}%)")

        if dev_pct < 0.1:
            print(f"  HIPOTEZA A=9 POTWIERDZONA (delta < 0.1%)")
        elif dev_pct < 1.0:
            print(f"  HIPOTEZA A=9 CZESCIOWO POTWIERDZONA (delta < 1%)")
        else:
            print(f"  HIPOTEZA A=9 OBALONA (delta >= 1%)")

        # Residua dla każdego punktu
        print(f"\n  Residua fitu:")
        for ag, ak, pr, res in zip(agams_arr, aksns_arr, aksns_pred, residuals):
            print(f"    agam={ag:.3f}  alpha_SN={ak:.6f}  pred={pr:.6f}  res={res:+.6f}")

        # Dodatkowe: fit kwadratowy (ewentualnie)
        if len(agams_ok) >= 3:
            c2 = np.polyfit(agams_arr, aksns_arr, 2)
            A2 = c2[2]
            print(f"\n  Fit kwadratowy (sprawdzenie): A(0) = {A2:.6f}")

    else:
        print("  Za mało punktów do fitu (potrzeba >= 2).")

    # ── Podsumowanie ──────────────────────────────────────────────────────────
    print(f"\n{'='*72}")
    print("PODSUMOWANIE:")
    print(f"  alpha_SN(0.040) = {AKSN_REF:.6f}  [p105_v2, precyzja 0.0004]")
    for agam in AGAM_VALS:
        r = results[agam]
        if np.isfinite(r['ak_SN']) and not r['ref']:
            ratio_to_9 = r['ak_SN'] / 9.0
            print(f"  alpha_SN({agam:.3f}) = {r['ak_SN']:.6f} +/- {r['err']:.6f}  "
                  f"(vs 9.0: {(9.0-r['ak_SN']):.6f} = {(1-ratio_to_9)*100:.3f}%)")

    print(f"\n[Czas calkowity: {(time.time()-t0_total)/60:.1f} min]")
    print("="*72)
