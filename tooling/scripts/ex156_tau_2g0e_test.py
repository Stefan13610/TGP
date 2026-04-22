"""
ex156_tau_2g0e_test.py
======================
R9: Test hipotezy g0^tau = 2*g0^e

Odkrycie z ex155: g0^tau/g0^e = 1.9892 ≈ 2 (delta 0.5%)
Schemat: g0^e : g0^mu : g0^tau = 1 : phi : 2

PYTANIA:
1. Czy g0^tau = 2*g0^e daje Koide DOKLADNIE?
2. Czy mozna dokladniej wyznaczyc g0^e (z B_tail=0 lub innego warunku)?
3. Jaka jest DOKLADNA wartosc g0^tau z warunku K=2/3?
4. Czy g0^tau_Koide / g0^e = 2 dokladnie?
5. Dlaczego 1 : phi : 2?  (phi^2 = phi+1, wiec 2 ≈ phi+phi^{-1} = phi+0.618)
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
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

def koide_K(A_e, A_mu, A_tau):
    """K = (sum m_i) / (sum sqrt(m_i))^2, m_i = A_i^4"""
    m = np.array([A_e**4, A_mu**4, A_tau**4])
    sm = np.sqrt(m)
    return np.sum(m) / np.sum(sm)**2

def koide_residual(g0_tau, A_e, A_mu):
    """K - 2/3 jako funkcja g0_tau"""
    A_tau = A_of_g0(g0_tau)
    if A_tau < 1e-10:
        return 1.0  # far from zero
    return koide_K(A_e, A_mu, A_tau) - 2.0/3.0

print("=" * 72)
print("ex156: Test hipotezy g0^tau = 2*g0^e")
print("=" * 72)

# --- 1. Precyzyjne obliczenia dla g0^e = 0.8694 ---
g0_e = 0.8694
g0_mu = PHI * g0_e
g0_tau_2 = 2.0 * g0_e  # hipoteza

print(f"\n--- 1. Wartosci bazowe ---")
print(f"  g0_e    = {g0_e:.6f}")
print(f"  g0_mu   = {g0_mu:.6f} = phi*g0_e")
print(f"  g0_tau  = {g0_tau_2:.6f} = 2*g0_e (hipoteza)")
print(f"  g0_tau_PDG = 1.7294 (z dopasowania)")

A_e = A_of_g0(g0_e)
A_mu = A_of_g0(g0_mu)
A_tau_2 = A_of_g0(g0_tau_2)
A_tau_pdg = A_of_g0(1.7294)

print(f"\n  A_e     = {A_e:.8f}")
print(f"  A_mu    = {A_mu:.8f}")
print(f"  A_tau(2*g0_e) = {A_tau_2:.8f}")
print(f"  A_tau(PDG)    = {A_tau_pdg:.8f}")

# --- 2. Koide test ---
print(f"\n--- 2. Test Koide ---")

K_2 = koide_K(A_e, A_mu, A_tau_2)
K_pdg = koide_K(A_e, A_mu, A_tau_pdg)

r21_2 = (A_mu/A_e)**4
r31_2 = (A_tau_2/A_e)**4
r31_pdg = (A_tau_pdg/A_e)**4

print(f"  g0_tau = 2*g0_e = {g0_tau_2:.6f}:")
print(f"    K     = {K_2:.8f}  (target: 0.666667)")
print(f"    delta = {abs(K_2 - 2/3)/(2/3)*100:.4f}%")
print(f"    r21   = {r21_2:.2f}  (PDG: 206.768)")
print(f"    r31   = {r31_2:.2f}  (PDG: 3477.15)")

print(f"\n  g0_tau = 1.7294 (PDG):")
print(f"    K     = {K_pdg:.8f}")
print(f"    delta = {abs(K_pdg - 2/3)/(2/3)*100:.4f}%")
print(f"    r31   = {r31_pdg:.2f}")

# --- 3. Precyzyjna lokalizacja g0_tau z warunku K = 2/3 ---
print(f"\n--- 3. Dokladna lokalizacja g0_tau z K = 2/3 ---")

# Skan wstepny
g0_scan = np.linspace(1.5, 1.9, 40)
res_scan = []
for g0t in g0_scan:
    try:
        res = koide_residual(g0t, A_e, A_mu)
        res_scan.append((g0t, res))
    except:
        pass

# Znajdz przejscie przez zero
g0_koide = None
for i in range(len(res_scan)-1):
    g1, r1 = res_scan[i]
    g2, r2 = res_scan[i+1]
    if r1 * r2 < 0:
        # Bisekcja
        try:
            g0_koide = brentq(lambda g: koide_residual(g, A_e, A_mu),
                              g1, g2, xtol=1e-10, rtol=1e-12)
            break
        except Exception as e:
            print(f"  Bisekcja nieudana: {e}")

if g0_koide is not None:
    A_tau_K = A_of_g0(g0_koide)
    K_exact = koide_K(A_e, A_mu, A_tau_K)
    r31_K = (A_tau_K/A_e)**4

    print(f"  g0_tau(K=2/3) = {g0_koide:.10f}")
    print(f"  K             = {K_exact:.10f}")
    print(f"  r31           = {r31_K:.2f}")
    print(f"  A_tau         = {A_tau_K:.8f}")

    # Stosunek do g0_e
    ratio = g0_koide / g0_e
    print(f"\n  g0_tau(K=2/3) / g0_e = {ratio:.10f}")
    print(f"  2.0                   = 2.0000000000")
    print(f"  delta od 2            = {abs(ratio - 2.0)/2.0*100:.6f}%")

    # Inne stosunki
    print(f"\n  g0_tau / g0_mu = {g0_koide/g0_mu:.10f}")
    print(f"  1/phi          = {1/PHI:.10f}")
    print(f"  2/phi          = {2/PHI:.10f}")
    print(f"  phi-1          = {PHI-1:.10f}")

    # Algebraiczne stosunki
    print(f"\n--- 4. Testy algebraiczne ---")
    candidates = {
        "2": 2.0,
        "phi + phi^{-1}": PHI + 1/PHI,  # = phi + phi-1 = sqrt(5)
        "sqrt(5)": np.sqrt(5),
        "phi + 1/phi": PHI + 1/PHI,
        "3 - phi": 3 - PHI,
        "phi^2 - phi + 1": PHI**2 - PHI + 1,  # = 1
        "2*phi - 1": 2*PHI - 1,  # = sqrt(5)
        "phi^(3/2)": PHI**1.5,
        "e^(1/phi)": np.exp(1/PHI),
        "2*phi/(phi+1)": 2*PHI/(PHI+1),
        "phi/(phi-1)": PHI/(PHI-1),  # = phi^2
    }

    for label, val in candidates.items():
        delta = abs(ratio - val) / val * 100
        mark = " <<<" if delta < 1 else ""
        print(f"  g0_tau/g0_e vs {label:20s} = {val:.8f}  (delta = {delta:.4f}%){mark}")

# --- 5. Czy g0_e jest wyjatkowe? ---
print(f"\n--- 5. Co specjalnego jest w g0_e = 0.8694? ---")

# B_tail(g0_e) = ?
r_e, g_e = solve_substrate(g0_e)
B_e, C_e = extract_BC(r_e, g_e)
print(f"  B(g0_e) = {B_e:.8f}")
print(f"  C(g0_e) = {C_e:.8f}")
print(f"  A(g0_e) = {A_e:.8f}")

# g0_e w stosunku do phi
print(f"\n  g0_e = {g0_e:.6f}")
print(f"  1/phi = {1/PHI:.6f}  (delta = {abs(g0_e - 1/PHI)/g0_e*100:.2f}%)")
print(f"  phi-1 = {PHI-1:.6f}  (delta = {abs(g0_e - (PHI-1))/g0_e*100:.2f}%)")
print(f"  2/phi^2 = {2/PHI**2:.6f}  (delta = {abs(g0_e - 2/PHI**2)/g0_e*100:.2f}%)")
print(f"  phi/2 = {PHI/2:.6f}  (delta = {abs(g0_e - PHI/2)/g0_e*100:.2f}%)")
print(f"  sqrt(phi)-1+1/phi = {np.sqrt(PHI)-1+1/PHI:.6f}")
print(f"  2-phi = {2-PHI:.6f}  (delta = {abs(g0_e - (2-PHI))/g0_e*100:.2f}%)")
print(f"  1/(1+1/phi) = {1/(1+1/PHI):.6f}  (delta = {abs(g0_e - 1/(1+1/PHI))/g0_e*100:.2f}%)")

# --- 6. Testy z innymi g0_e ---
print(f"\n--- 6. Universalnosc: czy 1:phi:2 dziala dla INNYCH g0_e? ---")
print(f"  Dla kazdego g0_e: sprawdzamy K(g0_e, phi*g0_e, 2*g0_e)")

for g0e_test in [0.5, 0.6, 0.7, 0.8, 0.8694, 0.9, 0.95, 1.0, 1.05, 1.1]:
    g0m = PHI * g0e_test
    g0t = 2.0 * g0e_test
    if g0t > 2.2:
        continue
    Ae = A_of_g0(g0e_test)
    Am = A_of_g0(g0m)
    At = A_of_g0(g0t)
    if Ae < 1e-10 or Am < 1e-10 or At < 1e-10:
        continue
    K = koide_K(Ae, Am, At)
    r21 = (Am/Ae)**4
    r31 = (At/Ae)**4
    print(f"  g0_e={g0e_test:.4f}: K={K:.6f}, r21={r21:.1f}, r31={r31:.1f}")

# --- 7. Precyzyjne znalezienie g0_e z warunku K=2/3 przy schemacie 1:phi:2 ---
print(f"\n--- 7. g0_e z warunku K(g0_e, phi*g0_e, 2*g0_e) = 2/3 ---")

def koide_residual_g0e(g0e):
    g0m = PHI * g0e
    g0t = 2.0 * g0e
    Ae = A_of_g0(g0e)
    Am = A_of_g0(g0m)
    At = A_of_g0(g0t)
    if Ae < 1e-10 or Am < 1e-10 or At < 1e-10:
        return 1.0
    return koide_K(Ae, Am, At) - 2.0/3.0

# Skan
g0e_range = np.linspace(0.5, 1.05, 50)
res_e = []
for g0e in g0e_range:
    try:
        res = koide_residual_g0e(g0e)
        res_e.append((g0e, res))
    except:
        pass

# Przejscia przez zero
for i in range(len(res_e)-1):
    g1, r1 = res_e[i]
    g2, r2 = res_e[i+1]
    if r1 * r2 < 0:
        try:
            g0e_opt = brentq(koide_residual_g0e, g1, g2, xtol=1e-10, rtol=1e-12)
            g0m_opt = PHI * g0e_opt
            g0t_opt = 2.0 * g0e_opt
            Ae = A_of_g0(g0e_opt)
            Am = A_of_g0(g0m_opt)
            At = A_of_g0(g0t_opt)
            K = koide_K(Ae, Am, At)
            r21 = (Am/Ae)**4
            r31 = (At/Ae)**4

            print(f"\n  >>> ROZWIAZANIE: g0_e = {g0e_opt:.10f} <<<")
            print(f"  g0_mu = {g0m_opt:.10f}")
            print(f"  g0_tau = {g0t_opt:.10f}")
            print(f"  K = {K:.10f}")
            print(f"  r21 = {r21:.4f}  (PDG: 206.768)")
            print(f"  r31 = {r31:.4f}  (PDG: 3477.15)")
            print(f"  A_e = {Ae:.8f}, A_mu = {Am:.8f}, A_tau = {At:.8f}")

            # Porownanie z g0_e = 0.8694
            print(f"\n  g0_e_opt = {g0e_opt:.10f}")
            print(f"  g0_e_rem = 0.8694")
            print(f"  delta = {abs(g0e_opt - 0.8694)/0.8694*100:.4f}%")

            # Masy leptonowe (normalizacja do elektronu)
            m_e = Ae**4
            m_mu = Am**4
            m_tau = At**4
            print(f"\n  Masy (w jednostkach m_e):")
            print(f"  m_e/m_e   = 1.0000")
            print(f"  m_mu/m_e  = {m_mu/m_e:.4f}  (PDG: 206.768)")
            print(f"  m_tau/m_e = {m_tau/m_e:.4f}  (PDG: 3477.15)")

        except Exception as e:
            print(f"  Bisekcja nieudana: {e}")

# --- 8. Wnioski ---
print(f"\n{'='*72}")
print("WNIOSKI ex156")
print(f"{'='*72}")
print(f"""
  ODKRYCIE: g0_tau / g0_e = 1.9892 ≈ 2 (delta 0.5%)

  Schemat TGP:
    g0^e   = g0_e           (parametr bazowy)
    g0^mu  = phi * g0_e     (phi-FP)
    g0^tau = 2 * g0_e       (hipoteza "doubling")

  Stosunki: 1 : phi : 2

  Uwaga: phi + 1/phi = sqrt(5) ≈ 2.236, NIE 2
  Ale: 2/phi = {2/PHI:.4f} = phi - 1 + 1/phi - ...

  Kluczowe pytanie: DLACZEGO 2?
  - Czy 2 = phi + (2-phi) = phi + 0.382?
  - Czy to granica stabilnosci?
  - Czy 2*g0_e jest kosmologicznie wyroznionym punktem?

  Mozliwy mechanizm Z_3:
    Z_3 := warunek Koide K(A_e, A_mu, A_tau) = 2/3
    + phi-FP: g0^mu = phi * g0^e
    => g0^tau jest WYZNACZONE
    => g0^tau ≈ 2*g0^e (numerycznie)

  Pytanie otwarte: czy g0^tau = 2*g0^e jest DOKLADNE
  czy tylko przyblizenie numeryczne warunku K=2/3?
""")
