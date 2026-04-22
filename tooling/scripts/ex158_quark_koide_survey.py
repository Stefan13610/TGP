"""
ex158_quark_koide_survey.py
============================
R12: Rozszerzenie schematu phi-FP + Koide na sektor kwarkowy.

KROK 1: Zbadaj masy kwarkow PDG i sprawdz:
  - Koide K = (sum m)/(sum sqrt(m))^2 dla (u,c,t) i (d,s,b)
  - Czy istnieje phi-FP w stosunkach mas?
  - Jakie g0 daja poprawne stosunki r21, r31?

Masy kwarkow PDG (running masses w schemacie MS-bar):
  Leptony (pole masses):
    m_e = 0.511 MeV, m_mu = 105.658 MeV, m_tau = 1776.86 MeV

  Kwarki (MS-bar at mu = 2 GeV dla lekkich, pole mass dla ciezkich):
    m_u = 2.16 MeV,  m_c = 1270 MeV,   m_t = 172760 MeV
    m_d = 4.67 MeV,  m_s = 93.4 MeV,   m_b = 4180 MeV

  UWAGA: masy kwarkowe sa scheme-dependent! Uzywamy:
    - Lekkie (u,d,s): MS-bar at 2 GeV (PDG 2024)
    - Ciezkie (c,b): MS-bar at m_c, m_b
    - Top: pole mass

  Alternatywne konwencje daja inne K!
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np

PHI = (1 + np.sqrt(5)) / 2

# --- Masy PDG ---
# Leptony (pole mass, MeV)
m_e   = 0.51099895
m_mu  = 105.6583755
m_tau = 1776.86

# Kwarki — kilka konwencji
# Konwencja A: PDG 2024, MS-bar at standard scales
m_u_A = 2.16    # MS-bar at 2 GeV, +0.49 -0.26
m_d_A = 4.67    # MS-bar at 2 GeV, +0.48 -0.17
m_s_A = 93.4    # MS-bar at 2 GeV, +8.6 -3.4
m_c_A = 1270.0  # MS-bar at m_c
m_b_A = 4180.0  # MS-bar at m_b
m_t_A = 172760.0  # pole mass (direct measurements)

# Konwencja B: running masses at mu = 2 GeV (Xing et al.)
m_u_B = 2.16
m_d_B = 4.67
m_s_B = 93.4
m_c_B = 1090.0   # running at 2 GeV
m_b_B = 3990.0   # running at 2 GeV
m_t_B = 168800.0  # running at 2 GeV (from pole via RGE)

# Konwencja C: running at mu = 1 GeV (Koide-friendly)
m_u_C = 2.9      # approximate
m_d_C = 6.4
m_s_C = 128.0
m_c_C = 1350.0
m_b_C = 4800.0
m_t_C = 172000.0

def koide_K(m1, m2, m3):
    """K = (m1+m2+m3) / (sqrt(m1)+sqrt(m2)+sqrt(m3))^2"""
    sm = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return (m1 + m2 + m3) / sm**2

def koide_3K(m1, m2, m3):
    """3K — format uzywany w ex153"""
    return 3 * koide_K(m1, m2, m3)

print("=" * 72)
print("ex158: Koide survey — leptony i kwarki")
print("=" * 72)

# --- 1. Leptony (weryfikacja) ---
print("\n--- 1. LEPTONY (weryfikacja) ---")
K_lep = koide_K(m_e, m_mu, m_tau)
r21_lep = m_mu / m_e
r31_lep = m_tau / m_e
print(f"  m_e = {m_e:.4f}, m_mu = {m_mu:.4f}, m_tau = {m_tau:.2f} MeV")
print(f"  K = {K_lep:.8f}  (target: 0.66666667)")
print(f"  3K = {3*K_lep:.8f}  (target: 2.00000000)")
print(f"  delta K = {abs(K_lep - 2/3)/(2/3)*100:.6f}%")
print(f"  r21 = {r21_lep:.3f}, r31 = {r31_lep:.2f}")
print(f"  sqrt(m_mu/m_e) = {np.sqrt(r21_lep):.4f}")
print(f"  sqrt(m_tau/m_e) = {np.sqrt(r31_lep):.4f}")

# --- 2. Kwarki — Koide K ---
print("\n--- 2. KWARKI — Koide K ---")

sectors = {
    "A: (u,c,t) standard": (m_u_A, m_c_A, m_t_A),
    "A: (d,s,b) standard": (m_d_A, m_s_A, m_b_A),
    "B: (u,c,t) at 2 GeV": (m_u_B, m_c_B, m_t_B),
    "B: (d,s,b) at 2 GeV": (m_d_B, m_s_B, m_b_B),
    "C: (u,c,t) at 1 GeV": (m_u_C, m_c_C, m_t_C),
    "C: (d,s,b) at 1 GeV": (m_d_C, m_s_C, m_b_C),
}

for label, (m1, m2, m3) in sectors.items():
    K = koide_K(m1, m2, m3)
    r21 = m2 / m1
    r31 = m3 / m1
    delta = abs(K - 2/3) / (2/3) * 100
    mark = " <<<" if delta < 5 else ""
    print(f"  {label:30s}: K={K:.6f} (delta={delta:.2f}%), "
          f"r21={r21:.1f}, r31={r31:.1f}{mark}")

# --- 3. Rozszerzone formuły Koide ---
print("\n--- 3. Rozszerzone formuly Koide ---")
print("  Generalized: K_n = (sum m^n) / (sum m^{n/2})^2 dla roznych n")

for n_label, n_val in [("n=1 (standard)", 1), ("n=1/2", 0.5), ("n=2", 2), ("n=1/3", 1/3)]:
    print(f"\n  {n_label}:")
    for label, (m1, m2, m3) in [
        ("leptony", (m_e, m_mu, m_tau)),
        ("(u,c,t)", (m_u_A, m_c_A, m_t_A)),
        ("(d,s,b)", (m_d_A, m_s_A, m_b_A)),
    ]:
        s_num = m1**n_val + m2**n_val + m3**n_val
        s_den = (m1**(n_val/2) + m2**(n_val/2) + m3**(n_val/2))**2
        K_n = s_num / s_den
        delta = abs(K_n - 2/3) / (2/3) * 100
        mark = " <<<" if delta < 5 else ""
        print(f"    {label:15s}: K_{n_label[:3]}={K_n:.6f} (delta={delta:.2f}%){mark}")

# --- 4. Stosunki phi w masach kwarkowych ---
print("\n--- 4. Stosunki phi w masach kwarkowych ---")
print(f"  phi = {PHI:.6f}")

quark_masses = {
    "u": m_u_A, "d": m_d_A, "s": m_s_A,
    "c": m_c_A, "b": m_b_A, "t": m_t_A
}

# Sektor up-type
print("\n  Sektor (u, c, t):")
print(f"  m_c/m_u = {m_c_A/m_u_A:.2f}")
print(f"  m_t/m_c = {m_t_A/m_c_A:.2f}")
print(f"  m_t/m_u = {m_t_A/m_u_A:.2f}")
print(f"  sqrt(m_c/m_u) = {np.sqrt(m_c_A/m_u_A):.4f}")
print(f"  sqrt(m_t/m_c) = {np.sqrt(m_t_A/m_c_A):.4f}")
print(f"  (m_c/m_u)^(1/4) = {(m_c_A/m_u_A)**0.25:.4f}")
print(f"  (m_t/m_c)^(1/4) = {(m_t_A/m_c_A)**0.25:.4f}")

# Sektor down-type
print("\n  Sektor (d, s, b):")
print(f"  m_s/m_d = {m_s_A/m_d_A:.2f}")
print(f"  m_b/m_s = {m_b_A/m_s_A:.2f}")
print(f"  m_b/m_d = {m_b_A/m_d_A:.2f}")
print(f"  sqrt(m_s/m_d) = {np.sqrt(m_s_A/m_d_A):.4f}")
print(f"  sqrt(m_b/m_s) = {np.sqrt(m_b_A/m_s_A):.4f}")
print(f"  (m_s/m_d)^(1/4) = {(m_s_A/m_d_A)**0.25:.4f}")
print(f"  (m_b/m_s)^(1/4) = {(m_b_A/m_s_A)**0.25:.4f}")

# --- 5. Test phi-FP: czy r21 = (A_mu/A_e)^4 z g0_mu = phi*g0_e? ---
# W leptonach: m ~ A^4, A ~ A_tail(g0)
# Jesli phi-FP: A_mu/A_e dziala z ODE substratowego
# To: jakie g0 daja r21 kwarkowe?
print("\n--- 5. Implikowane g0 z ODE substratowego ---")
print("  W leptonach: g0_e = 0.8694, r21 = 206.77")
print("  Pytanie: jakie g0^u i g0^d daloby kwarkowe r21?")

# Lepton reference
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

R_MAX = 150.0

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

def koide_K_from_A(Ae, Am, At):
    m = np.array([Ae**4, Am**4, At**4])
    sm = np.sqrt(m)
    return np.sum(m) / np.sum(sm)**2

# r21 jako funkcja g0_base (z phi-FP: g0_2 = phi*g0_base)
print("\n  r21(g0_base) z ODE substratowe + phi-FP:")

g0_bases = np.linspace(0.5, 1.05, 40)
r21_of_g0 = []
for g0b in g0_bases:
    g0m = PHI * g0b
    if g0m > 2.15:
        r21_of_g0.append(np.nan)
        continue
    Ae = A_of_g0(g0b)
    Am = A_of_g0(g0m)
    if Ae > 1e-10:
        r21_of_g0.append((Am/Ae)**4)
    else:
        r21_of_g0.append(np.nan)

# Wypisz kluczowe wartosci
targets = {
    "r21_lep (m_mu/m_e)": 206.768,
    "r21_up (m_c/m_u)": m_c_A / m_u_A,
    "r21_down (m_s/m_d)": m_s_A / m_d_A,
}

for target_label, target_r21 in targets.items():
    # Znajdz g0 dajace ten r21
    found = False
    for i in range(len(r21_of_g0)-1):
        r1, r2 = r21_of_g0[i], r21_of_g0[i+1]
        if np.isnan(r1) or np.isnan(r2):
            continue
        if (r1 - target_r21) * (r2 - target_r21) < 0:
            # Interpolacja
            def r21_res(g0b):
                Am = A_of_g0(PHI * g0b)
                Ae = A_of_g0(g0b)
                if Ae < 1e-10: return 1e6
                return (Am/Ae)**4 - target_r21
            try:
                g0_opt = brentq(r21_res, g0_bases[i], g0_bases[i+1],
                                xtol=1e-10, rtol=1e-12)
                Ae = A_of_g0(g0_opt)
                Am = A_of_g0(PHI * g0_opt)
                r21_check = (Am/Ae)**4

                print(f"\n  {target_label} = {target_r21:.2f}:")
                print(f"    g0_base = {g0_opt:.8f}")
                print(f"    g0_2    = phi*g0 = {PHI*g0_opt:.8f}")
                print(f"    r21     = {r21_check:.4f}")
                found = True

                # Teraz szukaj g0_tau z K=2/3
                def koide_res(g0t):
                    At = A_of_g0(g0t)
                    if At < 1e-10: return 1.0
                    return koide_K_from_A(Ae, Am, At) - 2.0/3.0

                # Skan
                g0t_range = np.linspace(max(0.3, PHI*g0_opt*0.5),
                                        min(2.15, PHI*g0_opt*2.0), 80)
                for j in range(len(g0t_range)-1):
                    try:
                        kr1 = koide_res(g0t_range[j])
                        kr2 = koide_res(g0t_range[j+1])
                        if kr1 * kr2 < 0:
                            g0t_sol = brentq(koide_res, g0t_range[j], g0t_range[j+1],
                                             xtol=1e-10)
                            At = A_of_g0(g0t_sol)
                            K = koide_K_from_A(Ae, Am, At)
                            r31 = (At/Ae)**4
                            print(f"    g0_tau(K=2/3) = {g0t_sol:.8f}, r31 = {r31:.2f}")
                    except:
                        pass

            except Exception as e:
                print(f"    Bisekcja nieudana: {e}")

    if not found:
        # r21 poza zakresem — sprawdz ekstrema
        r21_arr = np.array([x for x in r21_of_g0 if not np.isnan(x)])
        if len(r21_arr) > 0:
            print(f"\n  {target_label} = {target_r21:.2f}: POZA ZAKRESEM")
            print(f"    r21 range: [{np.min(r21_arr):.1f}, {np.max(r21_arr):.1f}]")

# --- 6. Koide bez phi-FP ---
print("\n--- 6. Test: Koide BEZ phi-FP (ogolne g0_1, g0_2, g0_3) ---")
print("  Pytanie: czy Koide K=2/3 moze byc spelniony dla kwarkowych r21")
print("  z DOWOLNYMI g0, nie tylko phi-FP?")

# Dla sektora (d,s,b): r21 = m_s/m_d ≈ 20
# Sprawdzmy: przy jakim r21 i Koide K=2/3, jaki jest r31?
# K = 2/3 + r21 znane => r31 wyznaczone (algebraicznie!)

print("\n  Algebraiczne rozwiazanie: K = (1 + r21 + r31) / (1 + sqrt(r21) + sqrt(r31))^2 = 2/3")
print("  Rozwiazujemy dla r31 przy danym r21:")

def solve_r31_from_koide(r21_target):
    """
    K = (1 + r21 + r31) / (1 + sqrt(r21) + sqrt(r31))^2 = 2/3
    Let x = sqrt(r31). Then:
    3(1 + r21 + x^2) = 2(1 + sqrt(r21) + x)^2
    3 + 3*r21 + 3*x^2 = 2(1 + sqrt(r21))^2 + 4(1+sqrt(r21))*x + 2*x^2
    x^2 - 4(1+sqrt(r21))*x + 3 + 3*r21 - 2(1+sqrt(r21))^2 = 0
    """
    s = np.sqrt(r21_target)
    # a*x^2 + b*x + c = 0
    a = 1.0
    b = -4.0 * (1 + s)
    c = 3.0 + 3.0*r21_target - 2.0*(1 + s)**2
    disc = b**2 - 4*a*c
    if disc < 0:
        return []
    x1 = (-b + np.sqrt(disc)) / (2*a)
    x2 = (-b - np.sqrt(disc)) / (2*a)
    results = []
    for x in [x1, x2]:
        if x > 0:
            r31 = x**2
            # Verify
            K = (1 + r21_target + r31) / (1 + np.sqrt(r21_target) + np.sqrt(r31))**2
            results.append((r31, K))
    return results

test_cases = [
    ("leptony (r21=206.77)", 206.768, 3477.15),
    ("(u,c,t) r21=587.96", m_c_A/m_u_A, m_t_A/m_u_A),
    ("(d,s,b) r21=20.00", m_s_A/m_d_A, m_b_A/m_d_A),
]

for label, r21, r31_pdg in test_cases:
    solutions = solve_r31_from_koide(r21)
    print(f"\n  {label}:")
    print(f"    r21 = {r21:.2f}, r31_PDG = {r31_pdg:.2f}")
    for r31_sol, K_check in solutions:
        delta = abs(r31_sol - r31_pdg) / r31_pdg * 100
        mark = " <<<" if delta < 10 else ""
        print(f"    r31(K=2/3) = {r31_sol:.2f}, K_check = {K_check:.8f}, "
              f"delta_PDG = {delta:.2f}%{mark}")

# --- 7. Wnioski ---
print(f"\n{'='*72}")
print("WNIOSKI ex158")
print(f"{'='*72}")
print(f"""
  1. Koide K dla kwarkow:
     - Leptony: K = 0.66666 (exact 2/3)
     - (u,c,t): K = ? (zalezy od konwencji mas)
     - (d,s,b): K = ? (zalezy od konwencji mas)

  2. Kluczowy problem: masy kwarkowe sa scheme-dependent
     (MS-bar at different scales daja rozne K)

  3. Algebraicznie: K=2/3 + r21 => r31 wyznaczone
     Test: czy r31(Koide) = r31(PDG)?

  4. phi-FP + ODE substratowe:
     - Dziala dla leptonow (r21 = 207, r31 = 3477)
     - Ograniczenie: ODE ma skonczony zakres r21
""")
