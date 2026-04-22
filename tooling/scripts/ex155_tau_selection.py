"""
ex155_tau_selection.py
======================
R9: Szukanie zasady selekcji g0^tau w ODE substratowym.

ODE: g'' + (1/g)(g')^2 + (2/r)g' = 1-g  (K_sub = g^2)

Wiemy:
  g0_e  = 0.8694 (phi-FP)
  g0_mu = phi * g0_e = 1.4068 (phi-FP)
  g0_tau = 1.7294 (z dopasowania r31 do PDG)

r21 = 206.1, r31 = 3466, K = 2/3 (Koide) -- wszystko poprawne.

PYTANIE: Co SELEKCJONUJE g0_tau = 1.7294?
Kandydaci:
  (a) B_tail(g0_tau) = 0? -> SPRAWDZIC
  (b) g0_tau = phi * g0_mu / phi = g0_mu? -> NIE (1.73 != 1.41)
  (c) Granica sektora (n_cross)? -> SPRAWDZIC
  (d) g0_tau taki ze CV(A^2) = 1? -> SPRAWDZIC
  (e) Minimum jakiejs funkcji? -> SPRAWDZIC
  (f) g0_tau = phi * g0_e + g0_e = g0_e * (1+phi)? -> SPRAWDZIC
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp

PHI = (1 + np.sqrt(5)) / 2
R_MAX = 120.0

def solve_substrate(g0):
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-8)
        source = 1.0 - g
        cross = (1.0 / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / 3.0]
        return [gp, source - cross - 2.0 * gp / r]
    y0 = [g0, 0.0]
    sol = solve_ivp(rhs, (0.0, R_MAX), y0, rtol=1e-10, atol=1e-12,
                    max_step=0.04, dense_output=True)
    return sol.t, sol.y[0]

def extract_BC(r, g, rL=25, rR=80):
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

def count_crossings(r, g):
    """Policz ile razy g(r) przecina g=1."""
    u = g - 1.0
    signs = np.sign(u[1:]) * np.sign(u[:-1])
    return np.sum(signs < 0)

print("=" * 72)
print("ex155: Selekcja g0^tau w ODE substratowym")
print("=" * 72)

g0_e = 0.8694
g0_mu = PHI * g0_e
g0_tau_pdg = 1.7294

print(f"  g0_e   = {g0_e:.4f}")
print(f"  g0_mu  = {g0_mu:.4f} = phi*g0_e")
print(f"  g0_tau = {g0_tau_pdg:.4f} (PDG-dopasowane)")

# --- 1. Profil solitonu i n_cross ---
print("\n--- 1. Topologia solitonow (n_cross) ---")
g0_scan = np.linspace(0.50, 2.20, 80)
for g0 in g0_scan:
    r, g = solve_substrate(g0)
    nc = count_crossings(r, g)
    A = A_of_g0(g0)
    B, C = extract_BC(r, g)
    # Minimum g
    g_min = np.min(g)
    mark = ""
    if abs(g0 - g0_e) < 0.015: mark = " <- e"
    elif abs(g0 - g0_mu) < 0.015: mark = " <- mu"
    elif abs(g0 - g0_tau_pdg) < 0.015: mark = " <- tau"
    if nc > 0 or g0 > 1.3 or mark:
        print(f"  g0={g0:.4f}: n_cross={nc}, g_min={g_min:.4f}, "
              f"A={A:.4f}, B={B:.4f}{mark}")

# --- 2. Stosunki zlotej proporcji ---
print("\n--- 2. Stosunki zlotej proporcji ---")
ratios = {
    "phi*g0_e": PHI * g0_e,
    "phi^2*g0_e": PHI**2 * g0_e,
    "(1+phi)*g0_e": (1+PHI) * g0_e,
    "2*g0_e": 2 * g0_e,
    "phi*g0_mu": PHI * g0_mu,
    "g0_e + g0_mu": g0_e + g0_mu,
    "sqrt(g0_e*g0_mu)*phi": np.sqrt(g0_e*g0_mu)*PHI,
    "g0_mu + g0_e*(phi-1)": g0_mu + g0_e*(PHI-1),
    "g0_mu * sqrt(phi)": g0_mu * np.sqrt(PHI),
}
for label, val in ratios.items():
    delta = abs(val - g0_tau_pdg) / g0_tau_pdg * 100
    mark = " <<<" if delta < 2 else ""
    print(f"  {label:25s} = {val:.4f}  (delta = {delta:.1f}%){mark}")

# --- 3. CV(A^2) = 1 ---
print("\n--- 3. Skan: g0_tau dajace CV(A^2) = 1 ---")
A_e = A_of_g0(g0_e)
A_mu = A_of_g0(g0_mu)
print(f"  A_e  = {A_e:.6f}")
print(f"  A_mu = {A_mu:.6f}")

g0_tau_scan = np.linspace(1.0, 2.1, 100)
for g0t in g0_tau_scan:
    A_t = A_of_g0(g0t)
    if A_t < 1e-8:
        continue
    A2 = np.array([A_e**2, A_mu**2, A_t**2])
    cv = np.std(A2) / np.mean(A2)
    K = (A_e**4 + A_mu**4 + A_t**4) / (A_e**2 + A_mu**2 + A_t**2)**2
    # Koide: K = 2/3
    delta_K = abs(K - 2.0/3.0) / (2.0/3.0) * 100
    if abs(cv - 1.0) < 0.05 or delta_K < 2:
        r31 = (A_t/A_e)**4
        print(f"  g0_tau={g0t:.4f}: CV={cv:.4f}, K={K:.6f} "
              f"(delta_K={delta_K:.2f}%), r31={r31:.1f}")

# --- 4. Granica sektora n_cross 0->1 ---
print("\n--- 4. Granica sektora n_cross: 0 -> 1 ---")
prev_nc = 0
for g0 in np.linspace(0.5, 2.2, 200):
    r, g = solve_substrate(g0)
    nc = count_crossings(r, g)
    if nc != prev_nc:
        print(f"  Przejscie n_cross {prev_nc} -> {nc} przy g0 ~ {g0:.4f}")
        delta = abs(g0 - g0_tau_pdg) / g0_tau_pdg * 100
        print(f"    (g0_tau = {g0_tau_pdg:.4f}, delta = {delta:.1f}%)")
    prev_nc = nc

# --- 5. Wnioski ---
print(f"\n--- 5. Wnioski ---")
print(f"  g0_tau/g0_e  = {g0_tau_pdg/g0_e:.4f}")
print(f"  g0_tau/g0_mu = {g0_tau_pdg/g0_mu:.4f}")
print(f"  g0_mu/g0_e   = {g0_mu/g0_e:.4f} = phi = {PHI:.4f}")

if __name__ == "__main__":
    pass  # run inline
