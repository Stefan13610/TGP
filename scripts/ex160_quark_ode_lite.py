"""
ex160 lite: ODE substratowe dla kwarkow — kluczowe testy.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
R_MAX = 150.0

M_E=0.511; M_MU=105.658; M_TAU=1776.86
M_U=2.16; M_C=1270.0; M_T=172760.0
M_D=4.67; M_S=93.4; M_B=4180.0

def solve_substrate(g0):
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-8)
        if r < 1e-10:
            return [gp, (1-g - gp**2/g) / 3.0]
        return [gp, 1-g - gp**2/g - 2*gp/r]
    sol = solve_ivp(rhs, (0, R_MAX), [g0, 0], rtol=1e-10, atol=1e-12, max_step=0.05)
    return sol.t, sol.y[0]

def A_of_g0(g0):
    r, g = solve_substrate(g0)
    mask = (r >= 25) & (r <= 100)
    if np.sum(mask) < 40: return 0.0
    rf = r[mask]; df = (g[mask]-1)*rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return np.sqrt(bc[0]**2 + bc[1]**2)

def koide_K_A(Ae, Am, At):
    m = np.array([Ae**4, Am**4, At**4])
    return np.sum(m) / np.sum(np.sqrt(m))**2

print("=" * 72)
print("ex160 lite: ODE substratowe dla kwarkow")
print("=" * 72)

# Precompute A(g0) grid
print("\nPrecomputing A(g0)...", flush=True)
g0_grid = np.linspace(0.40, 2.12, 80)
A_grid = np.array([A_of_g0(g0) for g0 in g0_grid])
print("Done.", flush=True)

# Inverse: find g0 for given A target
def find_g0_for_A(A_target):
    """Find g0 such that A(g0) = A_target using grid + refinement."""
    solutions = []
    for i in range(len(A_grid)-1):
        if (A_grid[i]-A_target)*(A_grid[i+1]-A_target) < 0:
            try:
                g0s = brentq(lambda g: A_of_g0(g)-A_target,
                             g0_grid[i], g0_grid[i+1], xtol=1e-8)
                solutions.append(g0s)
            except: pass
    return solutions

# For each sector: find g0_1 from phi-FP + r21, then g0_3 from r31
sectors = [
    ("leptony", M_MU/M_E, M_TAU/M_E),
    ("down(d,s,b)", M_S/M_D, M_B/M_D),
    ("up(u,c,t)", M_C/M_U, M_T/M_U),
]

print(f"\n{'Sektor':>15} {'r21':>8} {'r31':>10} {'g0_1':>10} {'g0_3':>10} {'K_ODE':>10} {'g0_3/g0_1':>10}")
print("-"*80)

for label, r21, r31 in sectors:
    # Find g0_1 from phi-FP: (A(phi*g0)/A(g0))^4 = r21
    def r21_res(g0):
        a1 = A_of_g0(g0)
        a2 = A_of_g0(PHI*g0)
        if a1 < 1e-10: return 1e6
        return (a2/a1)**4 - r21

    # Scan for sign change
    g0_1 = None
    for i in range(len(g0_grid)-1):
        g0a, g0b = g0_grid[i], g0_grid[i+1]
        if PHI*g0b > 2.12: continue
        try:
            ra = r21_res(g0a); rb = r21_res(g0b)
            if ra*rb < 0:
                g0_1 = brentq(r21_res, g0a, g0b, xtol=1e-10)
                break
        except: pass

    if g0_1 is None:
        print(f"  {label:>15} {r21:8.1f} {r31:10.1f} {'—':>10} {'—':>10} {'—':>10} {'—':>10}")
        continue

    A1 = A_of_g0(g0_1)
    A2 = A_of_g0(PHI*g0_1)

    # Find g0_3 from r31: (A(g0_3)/A(g0_1))^4 = r31
    A3_target = A1 * r31**0.25
    g0_3_list = find_g0_for_A(A3_target)

    if not g0_3_list:
        # Check if A3_target is in range
        print(f"  {label:>15} {r21:8.1f} {r31:10.1f} {g0_1:10.6f} {'OOR':>10} {'—':>10} {'—':>10}")
        print(f"    A3_target={A3_target:.4f}, A_max={np.max(A_grid):.4f}")
        continue

    for g0_3 in g0_3_list:
        A3 = A_of_g0(g0_3)
        K = koide_K_A(A1, A2, A3)
        ratio = g0_3/g0_1
        r31_c = (A3/A1)**4
        print(f"  {label:>15} {r21:8.1f} {r31_c:10.1f} {g0_1:10.6f} {g0_3:10.6f} {K:10.6f} {ratio:10.6f}")

# Detailed analysis
print(f"\n{'='*72}")
print("DETALE")
print(f"{'='*72}")

for label, r21, r31 in sectors:
    def r21_res(g0):
        a1 = A_of_g0(g0)
        a2 = A_of_g0(PHI*g0)
        if a1 < 1e-10: return 1e6
        return (a2/a1)**4 - r21

    g0_1 = None
    for i in range(len(g0_grid)-1):
        g0a, g0b = g0_grid[i], g0_grid[i+1]
        if PHI*g0b > 2.12: continue
        try:
            ra = r21_res(g0a); rb = r21_res(g0b)
            if ra*rb < 0:
                g0_1 = brentq(r21_res, g0a, g0b, xtol=1e-10)
                break
        except: pass

    if g0_1 is None: continue

    A1 = A_of_g0(g0_1)
    A2 = A_of_g0(PHI*g0_1)
    A3_target = A1 * r31**0.25
    g0_3_list = find_g0_for_A(A3_target)

    print(f"\n  {label}:")
    print(f"    g0_1 = {g0_1:.8f}, A_1 = {A1:.8f}")
    print(f"    g0_2 = {PHI*g0_1:.8f}, A_2 = {A2:.8f}")
    print(f"    r21 = {(A2/A1)**4:.4f} (target {r21:.4f})")

    if g0_3_list:
        g0_3 = g0_3_list[-1]  # highest
        A3 = A_of_g0(g0_3)
        K = koide_K_A(A1, A2, A3)
        print(f"    g0_3 = {g0_3:.8f}, A_3 = {A3:.8f}")
        print(f"    r31 = {(A3/A1)**4:.2f} (target {r31:.2f})")
        print(f"    K = {K:.8f} (2/3 = 0.66666667)")
        print(f"    delta(K) = {abs(K-2/3)/(2/3)*100:.4f}%")
        print(f"    g0_3/g0_1 = {g0_3/g0_1:.6f}")
        print(f"    g0_2/g0_1 = {PHI:.6f} (phi)")

        # What K would be if g0_3 were selected by Koide?
        # Find g0_3 from K=2/3
        def koide_res(g0t):
            at = A_of_g0(g0t)
            if at < 1e-10: return 1.0
            return koide_K_A(A1, A2, at) - 2/3

        for i in range(len(g0_grid)-1):
            try:
                ka = koide_res(g0_grid[i]); kb = koide_res(g0_grid[i+1])
                if ka*kb < 0:
                    g0_k = brentq(koide_res, g0_grid[i], g0_grid[i+1], xtol=1e-8)
                    Ak = A_of_g0(g0_k)
                    r31_k = (Ak/A1)**4
                    print(f"    --- Koide K=2/3 predykcja ---")
                    print(f"    g0_3(K) = {g0_k:.8f}, r31(K) = {r31_k:.2f}, "
                          f"r31/r31_PDG = {r31_k/r31:.4f}")
            except: pass
    else:
        print(f"    g0_3: OUT OF RANGE (A3_target={A3_target:.4f} > A_max={np.max(A_grid):.4f})")

print(f"\n{'='*72}")
print("WNIOSKI")
print(f"{'='*72}")
