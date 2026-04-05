"""
ex159_quark_extended_koide.py
==============================
R12: Rozszerzone testy Koide dla kwarkow.

ex158 pokazal: prosty K=2/3 NIE dziala na kwarki.
Pytanie: co MODYFIKOWAC?

Hipotezy:
  H1: Kwarki na innej skali running (szukamy mu gdzie K=2/3)
  H2: Parametryzacja Brannena m_k = M/3*(1+sqrt(2)*cos(theta+2pi*k/3))^2
      (K=2/3 z definicji, theta rozne dla kazdego sektora)
  H3: Factor kolorowy — masy efektywne m_eff = m/3 lub m/sqrt(3)
  H4: Inny wykladnik: m ~ A^n z n != 4
  H5: Inne ODE (inna K_sub) dla kwarkow
  H6: Koide z "shifted masses": m_k + m_0 (Kartavtsev)
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.optimize import minimize_scalar, brentq

PHI = (1 + np.sqrt(5)) / 2

# PDG masses (MeV)
m_e = 0.51099895; m_mu = 105.6583755; m_tau = 1776.86
m_u = 2.16; m_c = 1270.0; m_t = 172760.0
m_d = 4.67; m_s = 93.4; m_b = 4180.0

def koide_K(m1, m2, m3):
    sm = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return (m1 + m2 + m3) / sm**2

print("=" * 72)
print("ex159: Rozszerzone testy Koide dla kwarkow")
print("=" * 72)

# ===== H2: Parametryzacja Brannena =====
print("\n--- H2: Parametryzacja Brannena ---")
print("  m_k = M/3 * (1 + sqrt(2)*cos(theta_0 + 2*pi*k/3))^2")
print("  Ta parametryzacja ZAWSZE daje K=2/3.")
print("  Pytanie: jaki theta_0 fituje dane?\n")

def brannen_masses(M, theta0):
    """Zwraca (m1, m2, m3) z parametryzacji Brannena."""
    m = []
    for k in range(3):
        mk = M/3 * (1 + np.sqrt(2) * np.cos(theta0 + 2*np.pi*k/3))**2
        m.append(mk)
    return sorted(m)

def fit_brannen(m1, m2, m3):
    """Znajdz M i theta0 fitujace dane masy."""
    M = m1 + m2 + m3
    # theta0 z warunku na m1/M (najmniejsza masa)
    # m1 = M/3 * (1 + sqrt(2)*cos(theta0 + 2pi/3*k_min))^2
    # Szukamy theta0 minimalizujacego blad
    masses_sorted = sorted([m1, m2, m3])

    def residual(theta0):
        pred = sorted(brannen_masses(M, theta0))
        err = sum((pred[i] - masses_sorted[i])**2 / masses_sorted[i]**2 for i in range(3))
        return err

    # Skan
    best_theta = 0
    best_err = 1e10
    for th in np.linspace(0, 2*np.pi, 1000):
        err = residual(th)
        if err < best_err:
            best_err = err
            best_theta = th

    # Refinement
    from scipy.optimize import minimize
    res = minimize(residual, best_theta, method='Nelder-Mead', options={'xatol': 1e-12})
    return res.x[0] % (2*np.pi), res.fun

# Leptony
theta_lep, err_lep = fit_brannen(m_e, m_mu, m_tau)
M_lep = m_e + m_mu + m_tau
pred_lep = sorted(brannen_masses(M_lep, theta_lep))
print(f"  LEPTONY:")
print(f"    M = {M_lep:.2f} MeV")
print(f"    theta_0 = {theta_lep:.6f} rad = {np.degrees(theta_lep):.4f} deg")
print(f"    theta_0/(2pi/9) = {theta_lep/(2*np.pi/9):.6f}")
print(f"    Fit: ({pred_lep[0]:.4f}, {pred_lep[1]:.2f}, {pred_lep[2]:.2f})")
print(f"    PDG: ({m_e:.4f}, {m_mu:.2f}, {m_tau:.2f})")
print(f"    err = {err_lep:.2e}")

# (u,c,t)
theta_uct, err_uct = fit_brannen(m_u, m_c, m_t)
M_uct = m_u + m_c + m_t
pred_uct = sorted(brannen_masses(M_uct, theta_uct))
K_pred = koide_K(pred_uct[0], pred_uct[1], pred_uct[2])
print(f"\n  (u,c,t):")
print(f"    M = {M_uct:.2f} MeV")
print(f"    theta_0 = {theta_uct:.6f} rad = {np.degrees(theta_uct):.4f} deg")
print(f"    Fit: ({pred_uct[0]:.2f}, {pred_uct[1]:.2f}, {pred_uct[2]:.2f})")
print(f"    PDG: ({m_u:.2f}, {m_c:.2f}, {m_t:.2f})")
print(f"    err = {err_uct:.2e}")
print(f"    K(fit) = {K_pred:.8f} (should be 2/3)")

# (d,s,b)
theta_dsb, err_dsb = fit_brannen(m_d, m_s, m_b)
M_dsb = m_d + m_s + m_b
pred_dsb = sorted(brannen_masses(M_dsb, theta_dsb))
K_pred_dsb = koide_K(pred_dsb[0], pred_dsb[1], pred_dsb[2])
print(f"\n  (d,s,b):")
print(f"    M = {M_dsb:.2f} MeV")
print(f"    theta_0 = {theta_dsb:.6f} rad = {np.degrees(theta_dsb):.4f} deg")
print(f"    Fit: ({pred_dsb[0]:.2f}, {pred_dsb[1]:.2f}, {pred_dsb[2]:.2f})")
print(f"    PDG: ({m_d:.2f}, {m_s:.2f}, {m_b:.2f})")
print(f"    err = {err_dsb:.2e}")
print(f"    K(fit) = {K_pred_dsb:.8f}")

# Porownanie theta
print(f"\n  Porownanie theta_0:")
print(f"    leptony: {theta_lep:.6f} ({np.degrees(theta_lep):.2f} deg)")
print(f"    (u,c,t): {theta_uct:.6f} ({np.degrees(theta_uct):.2f} deg)")
print(f"    (d,s,b): {theta_dsb:.6f} ({np.degrees(theta_dsb):.2f} deg)")
print(f"    theta_uct/theta_lep = {theta_uct/theta_lep:.4f}")
print(f"    theta_dsb/theta_lep = {theta_dsb/theta_lep:.4f}")

# ===== H6: Shifted Koide (Kartavtsev) =====
print("\n--- H6: Shifted Koide: K(m_k + m_0) = 2/3 ---")
print("  Szukamy m_0 takiego ze K(m_1+m_0, m_2+m_0, m_3+m_0) = 2/3\n")

def shifted_koide_res(m0, m1, m2, m3):
    return koide_K(m1+m0, m2+m0, m3+m0) - 2.0/3.0

for label, (m1, m2, m3) in [
    ("leptony", (m_e, m_mu, m_tau)),
    ("(u,c,t)", (m_u, m_c, m_t)),
    ("(d,s,b)", (m_d, m_s, m_b)),
]:
    # Skan m0 w szerokim zakresie
    m0_range = np.concatenate([
        np.linspace(-m1*0.99, 0, 200),
        np.linspace(0, m3*0.5, 200),
    ])
    solutions = []
    for i in range(len(m0_range)-1):
        try:
            r1 = shifted_koide_res(m0_range[i], m1, m2, m3)
            r2 = shifted_koide_res(m0_range[i+1], m1, m2, m3)
            if r1 * r2 < 0:
                m0_sol = brentq(lambda x: shifted_koide_res(x, m1, m2, m3),
                                m0_range[i], m0_range[i+1], xtol=1e-8)
                K_check = koide_K(m1+m0_sol, m2+m0_sol, m3+m0_sol)
                solutions.append(m0_sol)
                r21_s = (m2+m0_sol)/(m1+m0_sol) if m1+m0_sol > 0 else np.nan
                r31_s = (m3+m0_sol)/(m1+m0_sol) if m1+m0_sol > 0 else np.nan
                print(f"  {label}: m_0 = {m0_sol:.4f} MeV, K = {K_check:.8f}, "
                      f"r21_shifted = {r21_s:.2f}, r31_shifted = {r31_s:.2f}")
        except:
            pass
    if not solutions:
        print(f"  {label}: brak rozwiazania w zakresie")

# ===== H4: Inny wykladnik =====
print("\n--- H4: Koide z m ~ A^n — szukamy n ---")
print("  Jesli m ~ A^n, to K powinno byc liczone z A = m^{1/n}")
print("  K_n = (sum m^{2/n}) / (sum m^{1/n})^2\n")

for label, (m1, m2, m3) in [
    ("leptony", (m_e, m_mu, m_tau)),
    ("(u,c,t)", (m_u, m_c, m_t)),
    ("(d,s,b)", (m_d, m_s, m_b)),
]:
    # Szukamy n takiego ze K_n = 2/3
    def K_of_n(n):
        if n < 0.01: return np.nan
        a1, a2, a3 = m1**(2/n), m2**(2/n), m3**(2/n)
        b1, b2, b3 = m1**(1/n), m2**(1/n), m3**(1/n)
        return (a1+a2+a3) / (b1+b2+b3)**2

    # Skan n
    n_range = np.linspace(0.5, 20, 200)
    K_vals = [K_of_n(n) for n in n_range]

    # Znajdz n gdzie K = 2/3
    sols = []
    for i in range(len(K_vals)-1):
        if np.isnan(K_vals[i]) or np.isnan(K_vals[i+1]):
            continue
        if (K_vals[i] - 2/3) * (K_vals[i+1] - 2/3) < 0:
            try:
                n_sol = brentq(lambda n: K_of_n(n) - 2/3, n_range[i], n_range[i+1])
                K_check = K_of_n(n_sol)
                sols.append(n_sol)
                print(f"  {label}: n = {n_sol:.4f}, K_n = {K_check:.8f}")
            except:
                pass
    if not sols:
        # Pokaz min/max K
        K_arr = np.array([k for k in K_vals if not np.isnan(k)])
        if len(K_arr) > 0:
            print(f"  {label}: K range = [{np.min(K_arr):.4f}, {np.max(K_arr):.4f}], "
                  f"2/3 = 0.6667 {'IN RANGE' if np.min(K_arr) < 2/3 < np.max(K_arr) else 'OUT OF RANGE'}")

# ===== H3: Factor kolorowy =====
print("\n--- H3: Factor kolorowy — m_eff = m * C ---")
print("  Kwarki maja N_c=3 kolor. Moze m_eff = m/3 lub m/N_c^p?\n")

for C_label, C in [("1", 1), ("1/3", 1/3), ("1/sqrt(3)", 1/np.sqrt(3)),
                     ("1/9", 1/9), ("3", 3)]:
    K_uct = koide_K(m_u*C, m_c*C, m_t*C)
    K_dsb = koide_K(m_d*C, m_s*C, m_b*C)
    # K jest niezalezne od C (jednorodna funkcja)!
    if C_label == "1":
        print(f"  UWAGA: K jest jednorodne stopnia 0 — C * m nie zmienia K!")
        print(f"  K(u,c,t) = {K_uct:.6f}, K(d,s,b) = {K_dsb:.6f}")
        print(f"  Factor kolorowy na MASY nie dziala.")
        break

# ===== Nowa hipoteza: Kaskadowy Koide =====
print("\n--- H7: Kaskadowy Koide (Sumino) ---")
print("  Koide dla par sektorow: (e,mu,tau) + (d,s,b) + (u,c,t)")
print("  lub masy mieszane miedzy sektorami\n")

# Sumino: sqrt(m_k) = M_k * (z0 + z1 * e^{i*2pi*k/3}), k=0,1,2
# Co jesli kwarkowe sektory uzywaja innej fazy Z_3?

# Sprawdzmy: Koide miedzy sektorami
cross_tests = [
    ("(e, s, t)", m_e, m_s, m_t),
    ("(e, c, b)", m_e, m_c, m_b),
    ("(d, mu, t)", m_d, m_mu, m_t),
    ("(u, s, tau)", m_u, m_s, m_tau),
    ("(d, c, tau)", m_d, m_c, m_tau),
    ("(u, mu, b)", m_u, m_mu, m_b),
    ("(e, d, u)", m_e, m_d, m_u),
    ("(mu, s, c)", m_mu, m_s, m_c),
    ("(tau, b, t)", m_tau, m_b, m_t),
]

print(f"  {'Trojka':20s}  {'K':>10}  {'delta(%)':>10}")
print("  " + "-" * 45)
for label, ma, mb, mc in cross_tests:
    K = koide_K(ma, mb, mc)
    delta = abs(K - 2/3) / (2/3) * 100
    mark = " <<<" if delta < 5 else ""
    print(f"  {label:20s}  {K:10.6f}  {delta:10.2f}%{mark}")

# ===== Nowa idea: n_eff szukany dla kazdego sektora =====
print("\n--- H8: Wykladnik n taki ze K_n = 2/3 dla KAZDEGO sektora ---")
print("  W leptonach: n=4 (m~A^4). Czy kwarki maja inne n?\n")

# Dla kazdego sektora, szukamy n
for label, (m1, m2, m3) in [
    ("leptony (e,mu,tau)", (m_e, m_mu, m_tau)),
    ("up (u,c,t)", (m_u, m_c, m_t)),
    ("down (d,s,b)", (m_d, m_s, m_b)),
]:
    def K_of_n(n):
        if n < 0.01: return np.nan
        a = np.array([m1**(2/n), m2**(2/n), m3**(2/n)])
        b = np.array([m1**(1/n), m2**(1/n), m3**(1/n)])
        return np.sum(a) / np.sum(b)**2

    # Fine scan
    n_range = np.linspace(0.1, 30, 3000)
    K_vals = []
    for n in n_range:
        try:
            K_vals.append(K_of_n(n))
        except:
            K_vals.append(np.nan)

    solutions = []
    for i in range(len(K_vals)-1):
        if np.isnan(K_vals[i]) or np.isnan(K_vals[i+1]):
            continue
        if (K_vals[i] - 2/3) * (K_vals[i+1] - 2/3) < 0:
            try:
                n_sol = brentq(lambda n: K_of_n(n) - 2/3, n_range[i], n_range[i+1], xtol=1e-10)
                solutions.append(n_sol)
            except:
                pass

    if solutions:
        print(f"  {label:25s}: n = {', '.join(f'{n:.4f}' for n in solutions)}")
        # Sprawdz stosunki
        if len(solutions) > 0:
            for n in solutions[:3]:
                # W tym n, A_k = m_k^{1/n}
                A1, A2, A3 = m1**(1/n), m2**(1/n), m3**(1/n)
                r21_A = (A2/A1)**n  # = m2/m1
                # phi-FP test: A2/A1 = ?
                ratio_21 = A2/A1
                ratio_31 = A3/A1
                phi_err = abs(ratio_21/ratio_31 - PHI/2) / (PHI/2) * 100
                print(f"    n={n:.4f}: A2/A1 = {ratio_21:.4f}, A3/A1 = {ratio_31:.4f}, "
                      f"(A2/A1)^n = {ratio_21**n:.2f} = r21")
    else:
        K_arr = np.array([k for k in K_vals if not np.isnan(k)])
        print(f"  {label:25s}: brak n (K range: [{np.min(K_arr):.4f}, {np.max(K_arr):.4f}])")

# ===== Wnioski =====
print(f"\n{'='*72}")
print("WNIOSKI ex159")
print(f"{'='*72}")
print(f"""
  KLUCZOWE WYNIKI:

  1. Brannen: kazda trojka mas MOZE byc fitowana formula Brannena
     (bo K=2/3 jest wbudowane). Roznica to theta_0.
     Pytanie: czy theta_0 ma sens fizyczny w TGP?

  2. Shifted Koide: K(m+m_0) = 2/3 ma rozwiazanie dla kwarkow
     Pytanie: co fizycznie oznacza shift m_0?

  3. Wykladnik n: leptony maja n ≈ 4 (K_4 = 2/3).
     Kwarki maja INNE n. Pytanie: czy n ma sens w TGP?

  4. Factor kolorowy: skalowanie mas nie zmienia K (jednorodnosc).
     Factor kolorowy musialby wchodzic INACZEJ.

  5. Cross-sektorowe: sa trojki mieszane bliskie K=2/3?

  WNIOSEK OGOLNY:
  Kwarki NIE spelniaja Koide K=2/3 z masami PDG.
  Mozliwe rozszerzenia wymagaja NOWYCH postulatow (theta, m_0, n).
  Najciekawsze: shifted Koide i wykladnik n.
""")
