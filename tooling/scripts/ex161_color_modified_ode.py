"""
ex161_color_modified_ode.py
============================
R12: Czy modyfikacja ODE z czynnikiem kolorowym daje Koide dla kwarkow?

IDEA: Kwarki maja kolor (N_c=3). W TGP, soliton kwarkowy moze miec
zmodyfikowane sprzezenie substratowe:
  K_sub(g) = g^(2*alpha_eff)  zamiast g^2

Albo zmodyfikowany potencjal:
  V'(g)/K(g) = (1-g) * C_color

Albo zmodyfikowany czlon kinetyczny:
  g'' + (alpha_eff/g)(g')^2 + (2/r)g' = 1-g

W ODE substratowym alpha_eff=1 (z K_sub=g^2).
Dla kwarkow moze alpha_eff != 1.

PLAN:
  1. Skan alpha_eff: dla kazdego alpha_eff, oblicz r21(g0,alpha) i K(alpha)
  2. Znajdz alpha_eff dajace K=2/3 dla down-type i up-type
  3. Czy ten alpha_eff jest prostym wyrazeniem z N_c?
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
R_MAX = 120.0

M_D=4.67; M_S=93.4; M_B=4180.0
M_U=2.16; M_C=1270.0; M_T=172760.0

def solve_ode_alpha(g0, alpha_eff, r_max=R_MAX):
    """ODE: g'' + (alpha_eff/g)(g')^2 + (2/r)g' = 1-g"""
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-8)
        src = 1.0 - g
        cross = (alpha_eff / g) * gp**2
        if r < 1e-10:
            return [gp, (src - cross) / 3.0]
        return [gp, src - cross - 2*gp/r]
    sol = solve_ivp(rhs, (0, r_max), [g0, 0], rtol=1e-10, atol=1e-12, max_step=0.05)
    return sol.t, sol.y[0]

def A_of_g0_alpha(g0, alpha_eff):
    r, g = solve_ode_alpha(g0, alpha_eff)
    mask = (r >= 25) & (r <= 90)
    if np.sum(mask) < 30: return 0.0
    rf = r[mask]; df = (g[mask]-1)*rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return np.sqrt(bc[0]**2 + bc[1]**2)

def koide_K_A(Ae, Am, At):
    m = np.array([Ae**4, Am**4, At**4])
    return np.sum(m) / np.sum(np.sqrt(m))**2

print("=" * 72)
print("ex161: ODE z czynnikiem kolorowym alpha_eff")
print("=" * 72)
print("ODE: g'' + (alpha/g)(g')^2 + (2/r)g' = 1-g")
print("Substratowe: alpha=1. Kanoniczne: alpha=2.\n")

# Skan alpha_eff dla down-type (r21=20, r31=895)
print("--- 1. Skan alpha_eff: sektor (d,s,b) ---")
print(f"  r21_target = {M_S/M_D:.2f}, r31_target = {M_B/M_D:.2f}\n")

r21_dsb = M_S/M_D
r31_dsb = M_B/M_D

alpha_range = [0.3, 0.5, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0]
print(f"  {'alpha':>6} {'g0_1':>10} {'r21':>10} {'g0_3':>10} {'r31':>10} {'K':>10} {'delta_K%':>10}")
print("  " + "-" * 70)

for alpha in alpha_range:
    # Find g0_1 from phi-FP + r21
    def r21_res(g0):
        a1 = A_of_g0_alpha(g0, alpha)
        a2 = A_of_g0_alpha(PHI*g0, alpha)
        if a1 < 1e-10: return 1e6
        return (a2/a1)**4 - r21_dsb

    g0_1 = None
    for g0_test in np.linspace(0.4, 1.1, 50):
        try:
            if PHI*g0_test > 3.0: continue
            r = r21_res(g0_test)
            if abs(r) < 50:
                # We're in the right ballpark, refine
                pass
        except: pass

    # Binary search
    g0_scan = np.linspace(0.4, 1.1, 40)
    for i in range(len(g0_scan)-1):
        if PHI*g0_scan[i+1] > 3.0: continue
        try:
            ra = r21_res(g0_scan[i])
            rb = r21_res(g0_scan[i+1])
            if ra * rb < 0:
                g0_1 = brentq(r21_res, g0_scan[i], g0_scan[i+1], xtol=1e-8)
                break
        except: pass

    if g0_1 is None:
        print(f"  {alpha:6.2f} {'—':>10}")
        continue

    A1 = A_of_g0_alpha(g0_1, alpha)
    A2 = A_of_g0_alpha(PHI*g0_1, alpha)
    r21_c = (A2/A1)**4

    # Find g0_3 from r31
    A3_target = A1 * r31_dsb**0.25

    g0_3 = None
    g0t_scan = np.linspace(max(0.3, g0_1*1.1), min(3.0, g0_1*3.5), 50)
    A_t_scan = [A_of_g0_alpha(g, alpha) for g in g0t_scan]

    for i in range(len(g0t_scan)-1):
        if (A_t_scan[i]-A3_target)*(A_t_scan[i+1]-A3_target) < 0:
            try:
                g0_3 = brentq(lambda g: A_of_g0_alpha(g, alpha)-A3_target,
                              g0t_scan[i], g0t_scan[i+1], xtol=1e-6)
                break
            except: pass

    if g0_3 is None:
        print(f"  {alpha:6.2f} {g0_1:10.6f} {r21_c:10.2f} {'OOR':>10}")
        continue

    A3 = A_of_g0_alpha(g0_3, alpha)
    r31_c = (A3/A1)**4
    K = koide_K_A(A1, A2, A3)
    delta_K = abs(K-2/3)/(2/3)*100

    mark = " <<<" if delta_K < 5 else ""
    print(f"  {alpha:6.2f} {g0_1:10.6f} {r21_c:10.2f} {g0_3:10.6f} {r31_c:10.1f} {K:10.6f} {delta_K:10.2f}%{mark}")

# Analogicznie dla up-type
print(f"\n--- 2. Skan alpha_eff: sektor (u,c,t) ---")
print(f"  r21_target = {M_C/M_U:.2f}, r31_target = {M_T/M_U:.2f}\n")

r21_uct = M_C/M_U
r31_uct = M_T/M_U

print(f"  {'alpha':>6} {'g0_1':>10} {'r21':>10} {'g0_3':>10} {'r31':>10} {'K':>10} {'delta_K%':>10}")
print("  " + "-" * 70)

for alpha in alpha_range:
    def r21_res_u(g0):
        a1 = A_of_g0_alpha(g0, alpha)
        a2 = A_of_g0_alpha(PHI*g0, alpha)
        if a1 < 1e-10: return 1e6
        return (a2/a1)**4 - r21_uct

    g0_1 = None
    g0_scan = np.linspace(0.4, 1.1, 40)
    for i in range(len(g0_scan)-1):
        if PHI*g0_scan[i+1] > 3.0: continue
        try:
            ra = r21_res_u(g0_scan[i])
            rb = r21_res_u(g0_scan[i+1])
            if ra * rb < 0:
                g0_1 = brentq(r21_res_u, g0_scan[i], g0_scan[i+1], xtol=1e-8)
                break
        except: pass

    if g0_1 is None:
        print(f"  {alpha:6.2f} {'—':>10}")
        continue

    A1 = A_of_g0_alpha(g0_1, alpha)
    A2 = A_of_g0_alpha(PHI*g0_1, alpha)
    r21_c = (A2/A1)**4

    A3_target = A1 * r31_uct**0.25
    g0_3 = None
    g0t_scan = np.linspace(max(0.3, g0_1*1.1), min(4.0, g0_1*4.0), 50)
    A_t_scan = [A_of_g0_alpha(g, alpha) for g in g0t_scan]

    for i in range(len(g0t_scan)-1):
        if (A_t_scan[i]-A3_target)*(A_t_scan[i+1]-A3_target) < 0:
            try:
                g0_3 = brentq(lambda g: A_of_g0_alpha(g, alpha)-A3_target,
                              g0t_scan[i], g0t_scan[i+1], xtol=1e-6)
                break
            except: pass

    if g0_3 is None:
        print(f"  {alpha:6.2f} {g0_1:10.6f} {r21_c:10.2f} {'OOR':>10}")
        continue

    A3 = A_of_g0_alpha(g0_3, alpha)
    r31_c = (A3/A1)**4
    K = koide_K_A(A1, A2, A3)
    delta_K = abs(K-2/3)/(2/3)*100

    mark = " <<<" if delta_K < 5 else ""
    print(f"  {alpha:6.2f} {g0_1:10.6f} {r21_c:10.2f} {g0_3:10.6f} {r31_c:10.1f} {K:10.6f} {delta_K:10.2f}%{mark}")

print(f"\n{'='*72}")
print("WNIOSKI ex161")
print(f"{'='*72}")
print("""
  Szukamy alpha_eff takiego ze Koide K=2/3 dla kwarkow.
  alpha_eff = 1: substratowe (leptony OK, kwarki NIE)
  alpha_eff = 2: kanoniczne z f(g) = 1+4ln(g)

  Jesli K(alpha) przechodzi przez 2/3 przy jakims alpha:
    => kwarki uzywaja INNEGO sprzezenia substratowego
    => N_c = 3 modyfikuje alpha_eff
""")
