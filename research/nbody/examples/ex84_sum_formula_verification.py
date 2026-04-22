"""
ex84_sum_formula_verification.py
==================================
STATUS: LEGACY-TRANSLATIONAL

This script belongs to the older `alpha*` / sum-formula verification program.
Treat it as historical exploratory context rather than a canonical synchronized
`nbody` path.

BEZPOSREDNIA WERYFIKACJA FORMULY SUMY S = 2*pi - 11/10

CEL:
  Sprawdzic, czy funkcja fz0(alpha) = (A(phi*z0(alpha)) / A(z0(alpha)))^4
  ma DOKLADNIE DWA ZERA rowne r21 = 206.768 w przedziale alpha in [2.3, 2.9],
  i czy ich suma S = alpha*_1 + alpha*_2 = 2*pi - 11/10.

ROZNICA od ex79/ex83:
  - ex79 hardkoduje alfa* z zewnetrznego zrodla; ex83 szuka przez g0star
  - ex84 liczy fz0(alpha) BEZPOSREDNIO na gestej siatce i brentq
  - Jeden tryb dzialania: wyznacz z0(alpha) -> oblicz f(z0, alpha) -> znajdz zera

WERYFIKACJA KRYTYCZNA (ex80 vs ex83):
  ex80 raportowal: alpha*_1 = 2.44143, S = 5.18318462 ~ 2pi-11/10 (0.13 ppm)
  ex83 raportowal: alpha*_1 = 2.44033, S = 5.18215 (rozbieznosc ~200 ppm)
  => Ktory wynik jest poprawny?
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')

# ========== STALE ==========
PHI       = (1.0 + np.sqrt(5.0)) / 2.0
R21       = 206.7682830   # CODATA
S_FORMULA = 2.0 * np.pi - 1.1
ALPHA_TGP = 2.0           # parametr ODE TGP (staly)
G_OFF     = 0.005
R_MAX     = 120.0         # dluzsze calkowanie dla lepszej amplitudy
WIN_COARSE = [(20, 36), (30, 46), (40, 56), (50, 66), (62, 78)]  # okna [rL, rR]

# ========== ODE ==========

def integrate_soliton(g0, alpha, r_max=R_MAX, rtol=1e-11, atol=1e-13):
    """
    Scalkuj ODE TGP: f(g)*g'' + (2/r)*f(g)*g' + (alpha/g)*(g')^2 = V'(g)
    z f(g) = 1 + 2*alpha*ln(g),  V'(g) = g^2*(1-g).
    Odbija sie od granicy duchowej gb = exp(-1/(2*alpha)) + G_OFF.
    Zwraca (r, g) po concatenacji wszystkich odcinkow.
    """
    gb = np.exp(-1.0 / (2.0 * alpha)) + G_OFF

    def rhs(r, y):
        g, gp = y
        g = max(g, gb + 1e-8)
        fg = 1.0 + 2.0 * alpha * np.log(g)
        if abs(fg) < 1e-10:
            return [gp, 0.0]
        Vprime = g * g * (1.0 - g)
        curl   = (alpha / g) * gp * gp
        if r < 1e-10:
            return [gp, (Vprime - curl) / (3.0 * fg)]
        return [gp, (Vprime - curl - fg * 2.0 * gp / r) / fg]

    def ev_bounce(r, y):
        return y[0] - gb
    ev_bounce.terminal  = True
    ev_bounce.direction = -1

    y0 = [g0, 0.0]; r0 = 1e-10
    ra, ga = [], []
    for _ in range(10):
        sol = solve_ivp(rhs, (r0, r_max), y0, events=ev_bounce,
                        dense_output=True, rtol=rtol, atol=atol,
                        max_step=0.03)
        ra.append(sol.t); ga.append(sol.y[0])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh  = sol.t_events[0][-1]
            st  = sol.sol(rh)
            y0  = [st[0], -st[1]]
            r0  = rh
        else:
            break
    r_all = np.concatenate(ra)
    g_all = np.concatenate(ga)
    idx   = np.argsort(r_all)
    return r_all[idx], g_all[idx]


def fit_amplitude_B(r, g, rL, rR):
    """
    Dopasuj (g(r)-1)*r = B*cos(r) + C*sin(r) w oknie [rL, rR].
    Zwraca A = sqrt(B^2 + C^2) i wspolczynnik B.
    """
    mask = (r >= rL) & (r <= rR)
    if mask.sum() < 30:
        return 0.0, 0.0, 0.0
    rf  = r[mask]
    df  = (g[mask] - 1.0) * rf
    M   = np.column_stack([np.cos(rf), np.sin(rf)])
    bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
    B, C = bc[0], bc[1]
    return float(np.sqrt(B*B + C*C)), float(B), float(C)


def A_inf_extrapolated(g0, alpha):
    """
    Oszacuj A_tail = lim_{r->inf} |g(r)-1|*r przez ekstrapolacje z okien.
    Dopasowanie: A(rL) = A_inf * (1 + a/rL + b/rL^2).
    """
    r, g = integrate_soliton(g0, alpha)
    Av, rv = [], []
    for rL, rR in WIN_COARSE:
        if rR > r[-1]:
            break
        A, _, _ = fit_amplitude_B(r, g, rL, rR)
        if A > 1e-4:
            Av.append(A)
            rv.append(float(rL))
    if len(Av) < 2:
        return Av[-1] if Av else 0.0
    rv = np.array(rv, dtype=float)
    Av = np.array(Av, dtype=float)
    # Dopasowanie z korekcja 1/r i 1/r^2
    try:
        p, _ = curve_fit(
            lambda x, ai, a, b: ai * (1.0 + a/x + b/x**2),
            rv, Av, p0=[Av[-1], 0., 0.], maxfev=5000
        )
        return float(p[0])
    except Exception:
        try:
            p, _ = curve_fit(
                lambda x, ai, a: ai * (1.0 + a/x),
                rv, Av, p0=[Av[-1], 0.], maxfev=2000
            )
            return float(p[0])
        except Exception:
            return float(Av[-1])


def B_coeff_window(g0, alpha, rL=28.0, rR=42.0):
    """Wspolczynnik B ogona oscylacyjnego w oknie [rL, rR]."""
    r, g = integrate_soliton(g0, alpha)
    _, B, _ = fit_amplitude_B(r, g, rL, rR)
    return B


def find_z0(alpha, lo=1.05, hi=2.5, n_scan=60):
    """
    Znajdz z0(alpha): wartosc g0 gdzie B_coeff(g0, alpha) = 0.
    Pierwsze zero (sektor n=0, 'elektron').
    """
    gs = np.linspace(lo, hi, n_scan)
    Bs = []
    for g in gs:
        try:
            Bs.append(B_coeff_window(g, alpha))
        except Exception:
            Bs.append(np.nan)
    Bs = np.array(Bs)
    # Znajdz pierwsze przejscie przez zero
    for i in range(len(Bs) - 1):
        if np.isnan(Bs[i]) or np.isnan(Bs[i+1]):
            continue
        if Bs[i] * Bs[i+1] < 0:
            try:
                return brentq(lambda x: B_coeff_window(x, alpha),
                              gs[i], gs[i+1], xtol=1e-10, rtol=1e-12,
                              maxiter=200)
            except Exception:
                pass
    return None


def fz0(alpha):
    """
    f(z0(alpha)) = (A(phi*z0) / A(z0))^4.
    Zwraca NaN jesli z0 nie znalezione.
    """
    z0 = find_z0(alpha)
    if z0 is None:
        return np.nan
    Ae = A_inf_extrapolated(z0, alpha)
    Am = A_inf_extrapolated(PHI * z0, alpha)
    if Ae < 1e-6:
        return np.nan
    return (Am / Ae) ** 4


def worker(alpha):
    return alpha, fz0(alpha)


# ========== MAIN ==========

if __name__ == '__main__':
    print("=" * 70)
    print("EX84: BEZPOSREDNIA WERYFIKACJA FORMULY SUMY S = 2pi - 11/10")
    print("=" * 70)
    print(f"  R21_CODATA  = {R21:.7f}")
    print(f"  S_formula   = 2*pi - 11/10 = {S_FORMULA:.10f}")
    print(f"  PHI (zlota) = {PHI:.10f}")
    print(f"  R_MAX       = {R_MAX}")
    print()

    # ----------------------------------------------------------------
    # KROK 1: Gesta siatka — znalezienie przedzialu z zerami
    # ----------------------------------------------------------------
    print("KROK 1: Skan fz0(alpha) na siatce [2.30, 2.85] (50 punktow)")
    print("-" * 70)
    N_CORES = min(16, cpu_count())
    print(f"  Uzywam {N_CORES} rdzeni")

    alpha_coarse = np.linspace(2.30, 2.85, 50)

    with Pool(N_CORES) as pool:
        results_coarse = pool.map(worker, alpha_coarse)

    avals_c = np.array([r[0] for r in results_coarse])
    fvals_c = np.array([r[1] for r in results_coarse])

    print(f"  {'alpha':>8}  {'f(z0)':>10}  {'f-R21':>9}")
    print(f"  {'-'*8}  {'-'*10}  {'-'*9}")
    crossings = []
    for i, (a, f) in enumerate(zip(avals_c, fvals_c)):
        if np.isnan(f):
            print(f"  {a:8.4f}  {'NaN':>10}  {'---':>9}")
            continue
        mark = ""
        if i > 0 and not np.isnan(fvals_c[i-1]):
            if (fvals_c[i-1] - R21) * (f - R21) < 0:
                crossings.append((avals_c[i-1], a))
                mark = "  <-- PRZEJSCIE"
        dist = f - R21
        print(f"  {a:8.4f}  {f:10.4f}  {dist:+9.4f}{mark}")

    print()
    print(f"  Znalezione przejscia przez R21: {len(crossings)}")
    for lo, hi in crossings:
        print(f"    [{lo:.5f}, {hi:.5f}]")
    print()

    # ----------------------------------------------------------------
    # KROK 2: Precyzyjne zera przez brentq w kazdym przejsciu
    # ----------------------------------------------------------------
    print("KROK 2: Precyzyjne zera fz0(alpha) = R21 przez brentq")
    print("-" * 70)

    alpha_stars = []
    for lo, hi in crossings:
        def f_minus_R(a):
            val = fz0(a)
            return (val - R21) if not np.isnan(val) else R21  # fallback

        try:
            a_star = brentq(f_minus_R, lo, hi, xtol=1e-10, rtol=1e-12, maxiter=200)
            f_star = fz0(a_star)
            alpha_stars.append(a_star)
            print(f"  alpha* = {a_star:.10f}   f(z0) = {f_star:.8f}   |f-R21| = {abs(f_star-R21):.4e}")
        except Exception as e:
            print(f"  BLAD brentq w [{lo:.4f}, {hi:.4f}]: {e}")

    print()

    # ----------------------------------------------------------------
    # KROK 3: Analiza sumy i porownanie z formula
    # ----------------------------------------------------------------
    print("=" * 70)
    print("KROK 3: ANALIZA SUMY I FORMULY")
    print("=" * 70)

    if len(alpha_stars) >= 2:
        a1, a2 = sorted(alpha_stars[:2])
        S = a1 + a2
        D = a2 - a1
        P = a1 * a2

        print(f"  alpha*_1  = {a1:.10f}")
        print(f"  alpha*_2  = {a2:.10f}")
        print()
        print(f"  S = alpha*_1 + alpha*_2 = {S:.10f}")
        print(f"  S_formula = 2*pi - 11/10 = {S_FORMULA:.10f}")
        dS = S - S_FORMULA
        print(f"  S - S_formula           = {dS:+.6e}  ({dS/S_FORMULA*1e6:+.3f} ppm)")
        print()
        print(f"  D = alpha*_2 - alpha*_1 = {D:.10f}")
        print(f"  3/10                    = {0.3:.10f}")
        print(f"  D - 3/10                = {D-0.3:+.6e}  ({(D-0.3)/0.3*1e6:+.0f} ppm)")
        print()
        print(f"  P = alpha*_1 * alpha*_2 = {P:.10f}")
        P_hyp = (np.pi - 0.7) * (np.pi - 0.4)
        print(f"  (pi-7/10)(pi-4/10)      = {P_hyp:.10f}")
        print(f"  P - P_hyp               = {P - P_hyp:+.6e}  ({(P-P_hyp)/P_hyp*1e6:+.0f} ppm)")
        print()
        print(f"  alpha*_1 ~ pi - 7/10 = {np.pi - 0.7:.10f}  err = {a1-(np.pi-0.7):+.6e}")
        print(f"  alpha*_2 ~ pi - 4/10 = {np.pi - 0.4:.10f}  err = {a2-(np.pi-0.4):+.6e}")
        print()
        print(f"  Weryfikacja ex80: alpha*_1 = 2.44143051, alpha*_2 = 2.74175411")
        print(f"  Ten skrypt:       alpha*_1 = {a1:.8f}, alpha*_2 = {a2:.8f}")
        print(f"  Roznica a1:       {a1 - 2.44143051:+.6e}")
        print(f"  Roznica a2:       {a2 - 2.74175411:+.6e}")
        print()

        # Dodatkowe kandydaty dla S
        print("  Kandydaci dla S:")
        cands = [
            ('2*pi - 11/10',  2*np.pi - 1.1),
            ('2*pi - ln(3)',   2*np.pi - np.log(3)),
            ('pi + e - 1/pi', np.pi + np.e - 1/np.pi),
            ('2*(pi - 11/20)', 2*(np.pi - 0.55)),
            ('3*pi - e',       3*np.pi - np.e),
        ]
        for name, val in sorted(cands, key=lambda c: abs(c[1]-S)):
            print(f"    {name:<30} = {val:.8f}  diff = {S-val:+.3e}  ({abs(S-val)/S*1e6:.2f} ppm)")

    elif len(alpha_stars) == 1:
        print(f"  Znaleziono tylko JEDNO zero: alpha* = {alpha_stars[0]:.8f}")
        print("  Mozliwe: drugie zero poza [2.30, 2.85] lub f_min > R21")
    else:
        print("  BRAK zer — f(z0, alpha) nie osiaga R21 w [2.30, 2.85]?")
        print("  Sprawdz wydruk KROKU 1.")

    # ----------------------------------------------------------------
    # KROK 4: z0(alpha*) — wartosci elektronowa i mionowa
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("KROK 4: WARTOSCI z0 PRZY ZERACH (kwantyzacja)")
    print("=" * 70)
    for i, a_star in enumerate(sorted(alpha_stars[:2])):
        z0 = find_z0(a_star)
        if z0 is not None:
            z0_phi = PHI * z0
            Ae = A_inf_extrapolated(z0, a_star)
            Am = A_inf_extrapolated(z0_phi, a_star)
            ratio = (Am/Ae)**4 if Ae > 1e-6 else np.nan
            print(f"  Galezia {i+1}: alpha* = {a_star:.8f}")
            print(f"    z0(alpha*)   = {z0:.8f}")
            print(f"    phi*z0       = {z0_phi:.8f}")
            print(f"    A(z0)        = {Ae:.8f}")
            print(f"    A(phi*z0)    = {Am:.8f}")
            print(f"    (Am/Ae)^4    = {ratio:.6f}  (cel: {R21:.4f})")
            print(f"    przy alpha_TGP=2: z0^e~1.24, z0^mu~2.00")
            print()

    # ----------------------------------------------------------------
    # TESTY
    # ----------------------------------------------------------------
    print("=" * 70)
    print("TESTY")
    print("=" * 70)

    n_pass = 0

    T1 = len(alpha_stars) >= 2
    print(f"  T1: f(z0,alpha) = R21 ma >= 2 rozwiazania: {'PASS' if T1 else 'FAIL'}")
    if T1: n_pass += 1

    if len(alpha_stars) >= 2:
        a1, a2 = sorted(alpha_stars[:2])
        S  = a1 + a2

        T2 = abs(S - S_FORMULA) / S_FORMULA < 10e-6  # < 10 ppm
        ppm_S = abs(S - S_FORMULA) / S_FORMULA * 1e6
        print(f"  T2: |S - (2pi-11/10)| < 10 ppm: {'PASS' if T2 else 'FAIL'}  ({ppm_S:.3f} ppm)")
        if T2: n_pass += 1

        T3 = abs(a1 - 2.44143051) < 0.002 and abs(a2 - 2.74175411) < 0.002
        print(f"  T3: alpha* blisko ex80-wartosci (0.002): {'PASS' if T3 else 'FAIL'}"
              f"  (da1={a1-2.44143051:+.4e}, da2={a2-2.74175411:+.4e})")
        if T3: n_pass += 1

        T4 = abs(D - 0.3) / 0.3 < 0.005  # < 0.5%
        print(f"  T4: D ~ 3/10 (<0.5%): {'PASS' if T4 else 'FAIL'}  ({abs(D-0.3)/0.3*100:.3f}%)")
        if T4: n_pass += 1

        T5 = abs(P - P_hyp) / P_hyp < 0.001  # < 0.1%
        print(f"  T5: P ~ (pi-7/10)(pi-4/10) (<0.1%): {'PASS' if T5 else 'FAIL'}  ({abs(P-P_hyp)/P_hyp*100:.4f}%)")
        if T5: n_pass += 1
    else:
        print(f"  T2-T5: FAIL (brak 2 zer)")

    print(f"\nWYNIK: {n_pass}/5")

    # ----------------------------------------------------------------
    # WNIOSEK
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("WNIOSEK")
    print("=" * 70)
    if len(alpha_stars) >= 2:
        a1, a2 = sorted(alpha_stars[:2])
        S = a1 + a2
        dS = S - S_FORMULA
        ppm = dS / S_FORMULA * 1e6
        print(f"  Zera f(z0, alpha) = R21:  alpha*_1 = {a1:.8f},  alpha*_2 = {a2:.8f}")
        print(f"  Suma:  S = {S:.10f}")
        print(f"  Wzor:  2*pi - 11/10 = {S_FORMULA:.10f}")
        print(f"  Odchylenie:  {dS:+.4e}  ({ppm:+.3f} ppm)")
        if abs(ppm) < 1.0:
            print()
            print("  *** FORMULA S = 2*pi - 11/10 POTWIERDZONA (< 1 ppm) ***")
            print("  *** r21 wynika z geometrii ODE TGP bez wolnych parametrow ***")
        elif abs(ppm) < 10.0:
            print()
            print("  ** Formula S = 2*pi - 11/10 potwierdzona (< 10 ppm) **")
        else:
            print()
            print(f"  ! Formula S = 2*pi - 11/10 NIEZGODNA ({ppm:.1f} ppm)")
            print("  Mozliwe: G_OFF, WIN_COARSE lub R_MAX wymagaja kalibracji")
    print("=" * 70)
