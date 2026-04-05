"""
ex157_one_param_prediction.py
==============================
R9: 1-parametrowa predykcja widma leptonowego

STRUKTURA:
  ODE substratowe: g'' + (1/g)(g')^2 + (2/r)g' = 1-g
  Dwa warunki:
    (1) phi-FP:  g0^mu = phi * g0^e
    (2) Koide:   K(A_e, A_mu, A_tau) = 2/3  =>  g0^tau = g0^tau(g0^e)

  Jeden wolny parametr: g0^e
  Funkcja A(g0) wyznaczona przez ODE.

  r21 = (A(phi*g0_e) / A(g0_e))^4 — funkcja g0^e
  r31 = (A(g0_tau(g0_e)) / A(g0_e))^4 — funkcja g0^e

PLAN:
  1. Znalezc g0^e takie ze r21 = 206.768 (PDG)
  2. Z tym g0^e wyznaczyc g0^tau z K=2/3
  3. Obliczyc r31 — to jest PREDYKCJA bez parametrow
  4. Porownac z PDG r31 = 3477.15
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
R_MAX = 150.0

# PDG values
R21_PDG = 206.768  # m_mu/m_e
R31_PDG = 3477.15  # m_tau/m_e

def solve_substrate(g0, r_max=R_MAX):
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-8)
        source = 1.0 - g
        cross = (1.0 / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / 3.0]
        return [gp, source - cross - 2.0 * gp / r]
    y0 = [g0, 0.0]
    sol = solve_ivp(rhs, (0.0, r_max), y0, rtol=1e-12, atol=1e-14,
                    max_step=0.02, dense_output=True)
    return sol.t, sol.y[0]

def extract_BC(r, g, rL=25, rR=100):
    mask = (r >= rL) & (r <= rR)
    if np.sum(mask) < 40:
        return 0.0, 0.0
    rf = r[mask]
    df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return bc[0], bc[1]

def A_of_g0(g0):
    r, g = solve_substrate(g0)
    B, C = extract_BC(r, g)
    return np.sqrt(B**2 + C**2)

def koide_K(A_e, A_mu, A_tau):
    m = np.array([A_e**4, A_mu**4, A_tau**4])
    sm = np.sqrt(m)
    return np.sum(m) / np.sum(sm)**2

# Cache A(g0) to avoid recomputation
_A_cache = {}
def A_cached(g0):
    g0r = round(g0, 12)
    if g0r not in _A_cache:
        _A_cache[g0r] = A_of_g0(g0)
    return _A_cache[g0r]

print("=" * 72)
print("ex157: 1-parametrowa predykcja widma leptonowego")
print("=" * 72)
print(f"  ODE: g'' + (1/g)(g')^2 + (2/r)g' = 1-g")
print(f"  phi-FP: g0^mu = phi * g0^e")
print(f"  Koide: K = 2/3")
print(f"  PDG: r21 = {R21_PDG}, r31 = {R31_PDG}")

# --- 1. r21(g0_e) jako funkcja g0_e ---
print(f"\n--- 1. r21(g0_e) ---")
g0e_scan = np.linspace(0.75, 0.95, 30)
r21_scan = []

for g0e in g0e_scan:
    Ae = A_of_g0(g0e)
    Am = A_of_g0(PHI * g0e)
    if Ae > 1e-10:
        r21 = (Am/Ae)**4
    else:
        r21 = np.nan
    r21_scan.append(r21)
    if abs(r21 - R21_PDG) < 50 or g0e in [g0e_scan[0], g0e_scan[-1]]:
        print(f"  g0_e={g0e:.4f}: r21={r21:.2f}")

# --- 2. Dokladna lokalizacja g0_e z r21 = 206.768 ---
print(f"\n--- 2. g0_e z r21 = {R21_PDG} ---")

def r21_residual(g0e):
    Ae = A_of_g0(g0e)
    Am = A_of_g0(PHI * g0e)
    if Ae < 1e-10:
        return 1e6
    r21 = (Am/Ae)**4
    return r21 - R21_PDG

# Szukamy przejscia przez zero
r21_arr = np.array(r21_scan)
g0e_opt = None
for i in range(len(r21_arr)-1):
    if np.isnan(r21_arr[i]) or np.isnan(r21_arr[i+1]):
        continue
    r1 = r21_arr[i] - R21_PDG
    r2 = r21_arr[i+1] - R21_PDG
    if r1 * r2 < 0:
        try:
            g0e_opt = brentq(r21_residual, g0e_scan[i], g0e_scan[i+1],
                             xtol=1e-12, rtol=1e-14)
            break
        except:
            pass

if g0e_opt is None:
    print("  NIE ZNALEZIONO — rozszerzam skan")
    g0e_scan2 = np.linspace(0.60, 0.98, 100)
    for g0e in g0e_scan2:
        Ae = A_of_g0(g0e)
        Am = A_of_g0(PHI * g0e)
        if Ae > 1e-10:
            r21 = (Am/Ae)**4
            if abs(r21 - R21_PDG) < 5:
                print(f"  Blisko: g0_e={g0e:.6f}, r21={r21:.2f}")

if g0e_opt is not None:
    Ae = A_of_g0(g0e_opt)
    Am = A_of_g0(PHI * g0e_opt)
    r21_check = (Am/Ae)**4

    print(f"\n  >>> g0_e(r21=PDG) = {g0e_opt:.12f} <<<")
    print(f"  g0_mu = phi*g0_e  = {PHI*g0e_opt:.12f}")
    print(f"  A_e               = {Ae:.10f}")
    print(f"  A_mu              = {Am:.10f}")
    print(f"  r21               = {r21_check:.6f}  (target: {R21_PDG})")

    # --- 3. g0_tau z K=2/3 ---
    print(f"\n--- 3. g0_tau z K = 2/3 ---")

    def koide_res_tau(g0t):
        At = A_of_g0(g0t)
        if At < 1e-10:
            return 1.0
        return koide_K(Ae, Am, At) - 2.0/3.0

    # Skan
    g0t_scan = np.linspace(1.4, 2.0, 60)
    kres = []
    for g0t in g0t_scan:
        try:
            kr = koide_res_tau(g0t)
            kres.append((g0t, kr))
        except:
            pass

    g0_tau_opt = None
    for i in range(len(kres)-1):
        g1, r1 = kres[i]
        g2, r2 = kres[i+1]
        if r1 * r2 < 0:
            try:
                g0_tau_opt = brentq(koide_res_tau, g1, g2, xtol=1e-12, rtol=1e-14)
                break
            except:
                pass

    if g0_tau_opt is not None:
        At = A_of_g0(g0_tau_opt)
        K_check = koide_K(Ae, Am, At)
        r31_pred = (At/Ae)**4

        print(f"  g0_tau(K=2/3)     = {g0_tau_opt:.12f}")
        print(f"  A_tau             = {At:.10f}")
        print(f"  K                 = {K_check:.12f}  (target: 0.666667)")
        print(f"")
        print(f"  ===================================")
        print(f"  PREDYKCJA: r31    = {r31_pred:.4f}")
        print(f"  PDG:       r31    = {R31_PDG}")
        print(f"  DELTA:            = {abs(r31_pred - R31_PDG)/R31_PDG*100:.4f}%")
        print(f"  ===================================")

        # Stosunki g0
        print(f"\n  g0_tau/g0_e  = {g0_tau_opt/g0e_opt:.8f}")
        print(f"  g0_tau/g0_mu = {g0_tau_opt/(PHI*g0e_opt):.8f}")
        print(f"  g0_mu/g0_e   = {PHI:.8f} (phi, exact)")

        # m_tau predykcja
        m_e_MeV = 0.511  # MeV
        m_mu_MeV = 105.658  # MeV
        m_tau_PDG_MeV = 1776.86  # MeV
        m_tau_pred_MeV = r31_pred * m_e_MeV

        print(f"\n  Predykcja m_tau:")
        print(f"    m_tau(pred) = r31 * m_e = {r31_pred:.2f} * {m_e_MeV} MeV")
        print(f"    m_tau(pred) = {m_tau_pred_MeV:.2f} MeV")
        print(f"    m_tau(PDG)  = {m_tau_PDG_MeV:.2f} MeV")
        print(f"    delta       = {abs(m_tau_pred_MeV - m_tau_PDG_MeV)/m_tau_PDG_MeV*100:.4f}%")

        # --- 4. Sprawdzenie: ile parametrow? ---
        print(f"\n--- 4. Podsumowanie parametrow ---")
        print(f"  INPUT:")
        print(f"    ODE substratowe (0 parametrow)")
        print(f"    phi-FP: g0_mu = phi*g0_e (0 parametrow)")
        print(f"    Koide: K = 2/3 (0 parametrow)")
        print(f"    r21 = {R21_PDG} (1 parametr — wyznacza g0_e)")
        print(f"  OUTPUT:")
        print(f"    r31 = {r31_pred:.4f} (PREDYKCJA)")
        print(f"    PDG = {R31_PDG}")

        # --- 5. Czulosc na bledy ---
        print(f"\n--- 5. Czulosc ---")
        # Jak zmienia sie r31 gdy r21 zmieni sie o 1%?
        for dr in [-0.01, -0.001, 0.001, 0.01]:
            r21_test = R21_PDG * (1 + dr)
            def r21_res_test(g0e):
                ae = A_of_g0(g0e)
                am = A_of_g0(PHI * g0e)
                if ae < 1e-10: return 1e6
                return (am/ae)**4 - r21_test

            try:
                g0e_t = brentq(r21_res_test, g0e_opt*0.99, g0e_opt*1.01,
                               xtol=1e-10, rtol=1e-12)
                ae_t = A_of_g0(g0e_t)
                am_t = A_of_g0(PHI * g0e_t)

                def kr_t(g0t):
                    at = A_of_g0(g0t)
                    if at < 1e-10: return 1.0
                    return koide_K(ae_t, am_t, at) - 2.0/3.0

                g0t_t = brentq(kr_t, g0_tau_opt*0.95, g0_tau_opt*1.05,
                               xtol=1e-10, rtol=1e-12)
                at_t = A_of_g0(g0t_t)
                r31_t = (at_t/ae_t)**4
                print(f"  r21*(1{dr:+.3f}) = {r21_test:.3f}: "
                      f"g0_e={g0e_t:.8f}, r31={r31_t:.2f} "
                      f"(delta_r31 = {(r31_t-r31_pred)/r31_pred*100:+.2f}%)")
            except:
                print(f"  r21*(1{dr:+.3f}): bisekcja nieudana")

    # --- 6. Alternatywne rozwiazanie: czy moga byc inne g0_tau z K=2/3? ---
    print(f"\n--- 6. Wszystkie g0_tau z K = 2/3 ---")
    g0t_full = np.linspace(0.5, 2.1, 200)
    kres_full = []
    for g0t in g0t_full:
        try:
            kr = koide_res_tau(g0t)
            kres_full.append((g0t, kr))
        except:
            pass

    n_solutions = 0
    for i in range(len(kres_full)-1):
        g1, r1 = kres_full[i]
        g2, r2 = kres_full[i+1]
        if r1 * r2 < 0:
            try:
                g0t_sol = brentq(koide_res_tau, g1, g2, xtol=1e-10)
                At_sol = A_of_g0(g0t_sol)
                r31_sol = (At_sol/Ae)**4
                n_solutions += 1
                mark = " <<< PDG" if abs(r31_sol - R31_PDG) / R31_PDG < 0.05 else ""
                print(f"  Rozwiazanie #{n_solutions}: g0_tau={g0t_sol:.8f}, "
                      f"r31={r31_sol:.2f}{mark}")
            except:
                pass
    print(f"  Razem: {n_solutions} rozwiazan z K=2/3")

# --- Wnioski ---
print(f"\n{'='*72}")
print("WNIOSKI ex157")
print(f"{'='*72}")
print(f"""
  STRUKTURA PREDYKCYJNA:
    ODE substratowe + phi-FP + Koide
    1 parametr wejsciowy (r21 lub g0_e)
    => r31 jest PREDYKCJA

  KLUCZOWE PYTANIE: delta(r31) < 1%?
    Jesli TAK => TGP PREDUKUJE mase tau
    Jesli NIE => cos jest nie tak z ODE lub warunkami

  Status epistemologiczny:
    ODE substratowe  [AX] — aksjomat (K_sub = g^2)
    phi-FP           [POST+NUM] — postulat numeryczny
    Koide K=2/3      [POST] — postulat (Z_3)
    m ~ A^4          [POST+NUM]
    g0_e             [OBS] — z dopasowania r21
""")
