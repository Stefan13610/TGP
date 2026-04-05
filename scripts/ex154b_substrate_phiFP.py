"""
ex154b: Test phi-FP na ODE substratowym z wartosciami z rem:K-substrate-tau.
g0^e = 0.8694, g0^mu = phi*g0^e = 1.4068, g0^tau = 1.7294
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp

PHI = (1 + np.sqrt(5)) / 2
R_MAX = 120.0

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
    sol = solve_ivp(rhs, (0.0, r_max), y0, rtol=1e-10, atol=1e-12,
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

def koide_QK(m1, m2, m3):
    s1 = np.sqrt(abs(m1)) + np.sqrt(abs(m2)) + np.sqrt(abs(m3))
    return (abs(m1)+abs(m2)+abs(m3)) / s1**2 * 3 if s1 > 1e-15 else np.nan

print("=" * 72)
print("ex154b: ODE substratowe — test phi-FP z wartosciami rem:K-substrate-tau")
print("=" * 72)

# Wartosci z rem:K-substrate-tau
g0_e = 0.8694
g0_mu = PHI * g0_e  # 1.4068
g0_tau_rem = 1.7294

print(f"\n  g0_e   = {g0_e:.4f}")
print(f"  g0_mu  = {g0_mu:.4f} (phi*g0_e)")
print(f"  g0_tau = {g0_tau_rem:.4f} (z rem:K-substrate-tau)")

# Oblicz solitony
results = {}
for label, g0 in [("e", g0_e), ("mu", g0_mu), ("tau", g0_tau_rem)]:
    r, g = solve_substrate(g0)
    B, C = extract_BC(r, g)
    A = np.sqrt(B**2 + C**2)
    results[label] = {"g0": g0, "A": A, "B": B, "C": C, "r": r, "g": g}
    print(f"\n  {label}: g0={g0:.4f}, A={A:.6f}, B={B:.6f}, C={C:.6f}")

# Stosunki mas
A_e = results["e"]["A"]
A_mu = results["mu"]["A"]
A_tau = results["tau"]["A"]

print(f"\n--- Stosunki mas ---")
r21 = (A_mu/A_e)**4 if A_e > 0 else 0
r31 = (A_tau/A_e)**4 if A_e > 0 else 0
QK = koide_QK(A_e**4, A_mu**4, A_tau**4)

print(f"  r_21 = (A_mu/A_e)^4 = {r21:.2f}  (PDG: 206.768)")
print(f"  r_31 = (A_tau/A_e)^4 = {r31:.2f}  (PDG: 3477.15)")
print(f"  Q_K  = {QK:.6f}  (PDG: 1.500000)")

A2 = np.array([A_e**2, A_mu**2, A_tau**2])
cv = np.std(A2)/np.mean(A2)
print(f"  CV(A^2) = {cv:.4f}  (Z_3: 1.0000)")

# Skan g0_tau szukajac najlepszego Q_K i r_31
print(f"\n--- Skan g0_tau ---")
print(f"  {'g0_tau':>8}  {'A_tau':>10}  {'r_21':>10}  {'r_31':>10}  {'Q_K':>8}")
print("  " + "-" * 55)

for g0t in np.linspace(1.0, 2.2, 50):
    r, g = solve_substrate(g0t)
    B, C = extract_BC(r, g)
    A_t = np.sqrt(B**2 + C**2)
    if A_t < 1e-8:
        continue
    r31_t = (A_t/A_e)**4
    QK_t = koide_QK(A_e**4, A_mu**4, A_t**4)
    mark = ""
    if abs(QK_t - 1.5) < 0.02: mark = " <K>"
    if abs(r31_t - 3477) < 200: mark += " <r31>"
    print(f"  {g0t:8.4f}  {A_t:10.6f}  {r21:10.2f}  {r31_t:10.1f}  {QK_t:8.4f}{mark}")

# Test phi^2
g0_phi2 = PHI**2 * g0_e
r_p2, g_p2 = solve_substrate(g0_phi2)
B_p2, C_p2 = extract_BC(r_p2, g_p2)
A_p2 = np.sqrt(B_p2**2 + C_p2**2)
r31_p2 = (A_p2/A_e)**4 if A_e > 0 else 0
QK_p2 = koide_QK(A_e**4, A_mu**4, A_p2**4)
print(f"\n  phi^2*g0_e = {g0_phi2:.4f}: A={A_p2:.6f}, r_31={r31_p2:.1f}, Q_K={QK_p2:.4f}")
