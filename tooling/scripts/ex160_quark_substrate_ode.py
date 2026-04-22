"""
ex160_quark_substrate_ode.py
==============================
R12: ODE substratowe dla kwarkow — pelny test.

STRATEGIA:
  1. ODE substratowe jest UNIVERSALNE (nie zalezy od sektora)
  2. phi-FP: g0^(2) = phi * g0^(1) — universalny
  3. g0^(1) wyznaczone przez r21 = m_2/m_1
  4. g0^(3) wyznaczone przez r31 = m_3/m_1 (z PDG)
  5. Pytanie: czy K(g0_1, phi*g0_1, g0_3_PDG) ma jakis wzorzec?

SEKTORY:
  Leptony: r21 = 206.77, r31 = 3477.2 → K = 0.66667 (2/3!)
  (d,s,b): r21 = 20.0,   r31 = 895.1  → K = ?
  (u,c,t): r21 = 588.0,  r31 = 79982  → K = ?

HIPOTEZA ZEROWA:
  Moze K_ODE (z dopasowanego g0_3) NIE jest 2/3 dla kwarkow,
  ale jest jakims prostym wyrazeniem (np. 2/3 * f(N_c)).
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
R_MAX = 150.0

# PDG ratios
M_E = 0.511; M_MU = 105.658; M_TAU = 1776.86
M_U = 2.16; M_C = 1270.0; M_T = 172760.0
M_D = 4.67; M_S = 93.4; M_B = 4180.0

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

def koide_K_A(Ae, Am, At):
    m = np.array([Ae**4, Am**4, At**4])
    sm = np.sqrt(m)  # = A^2
    return np.sum(m) / np.sum(sm)**2

print("=" * 72)
print("ex160: ODE substratowe dla kwarkow — pelny test")
print("=" * 72)

# ===== 1. Dla kazdego sektora: znajdz g0_base z r21 =====
print("\n--- 1. g0_base z phi-FP + r21 ---")

sectors = {
    "leptony": {"m1": M_E, "m2": M_MU, "m3": M_TAU,
                "r21": M_MU/M_E, "r31": M_TAU/M_E},
    "down (d,s,b)": {"m1": M_D, "m2": M_S, "m3": M_B,
                      "r21": M_S/M_D, "r31": M_B/M_D},
    "up (u,c,t)": {"m1": M_U, "m2": M_C, "m3": M_T,
                    "r21": M_C/M_U, "r31": M_T/M_U},
}

results = {}

for label, sec in sectors.items():
    r21_target = sec["r21"]
    r31_target = sec["r31"]

    # Znajdz g0 dajace r21 z phi-FP
    def r21_res(g0b):
        Ae = A_of_g0(g0b)
        Am = A_of_g0(PHI * g0b)
        if Ae < 1e-10: return 1e6
        return (Am/Ae)**4 - r21_target

    # Skan
    g0_scan = np.linspace(0.50, 1.05, 80)
    r21_vals = []
    for g0 in g0_scan:
        try:
            Ae = A_of_g0(g0)
            Am = A_of_g0(PHI * g0)
            r21_vals.append((Am/Ae)**4 if Ae > 1e-10 else np.nan)
        except:
            r21_vals.append(np.nan)

    g0_base = None
    for i in range(len(r21_vals)-1):
        if np.isnan(r21_vals[i]) or np.isnan(r21_vals[i+1]):
            continue
        if (r21_vals[i] - r21_target) * (r21_vals[i+1] - r21_target) < 0:
            try:
                g0_base = brentq(r21_res, g0_scan[i], g0_scan[i+1],
                                 xtol=1e-12, rtol=1e-14)
                break
            except:
                pass

    if g0_base is None:
        print(f"\n  {label}: r21 = {r21_target:.2f} — NIE ZNALEZIONO g0_base!")
        print(f"    r21 range: [{min(x for x in r21_vals if not np.isnan(x)):.1f}, "
              f"{max(x for x in r21_vals if not np.isnan(x)):.1f}]")
        results[label] = None
        continue

    Ae = A_of_g0(g0_base)
    Am = A_of_g0(PHI * g0_base)
    r21_check = (Am/Ae)**4

    print(f"\n  {label}:")
    print(f"    r21_target = {r21_target:.4f}")
    print(f"    g0_base    = {g0_base:.10f}")
    print(f"    g0_2       = {PHI*g0_base:.10f}")
    print(f"    A_1        = {Ae:.8f}")
    print(f"    A_2        = {Am:.8f}")
    print(f"    r21_check  = {r21_check:.4f}")

    # Znajdz g0_3 dajace r31 z PDG
    def r31_res(g0t):
        At = A_of_g0(g0t)
        if At < 1e-10 or Ae < 1e-10: return 1e6
        return (At/Ae)**4 - r31_target

    g0t_scan = np.linspace(0.3, 2.15, 100)
    r31_vals = []
    for g0t in g0t_scan:
        try:
            At = A_of_g0(g0t)
            r31_vals.append((At/Ae)**4 if Ae > 1e-10 else np.nan)
        except:
            r31_vals.append(np.nan)

    g0_3_solutions = []
    for i in range(len(r31_vals)-1):
        if np.isnan(r31_vals[i]) or np.isnan(r31_vals[i+1]):
            continue
        if (r31_vals[i] - r31_target) * (r31_vals[i+1] - r31_target) < 0:
            try:
                g0_3 = brentq(r31_res, g0t_scan[i], g0t_scan[i+1],
                              xtol=1e-10)
                At = A_of_g0(g0_3)
                K = koide_K_A(Ae, Am, At)
                r31_check = (At/Ae)**4
                g0_3_solutions.append((g0_3, At, K, r31_check))
            except:
                pass

    if g0_3_solutions:
        for idx, (g0_3, At, K, r31c) in enumerate(g0_3_solutions):
            print(f"    g0_3 (sol #{idx+1}) = {g0_3:.10f}")
            print(f"    A_3            = {At:.8f}")
            print(f"    r31_check      = {r31c:.2f} (target: {r31_target:.2f})")
            print(f"    K(ODE)         = {K:.8f}")
            print(f"    delta(K-2/3)   = {abs(K-2/3)/(2/3)*100:.4f}%")
            print(f"    g0_3/g0_1      = {g0_3/g0_base:.6f}")
            print(f"    g0_3/g0_2      = {g0_3/(PHI*g0_base):.6f}")

        results[label] = {
            "g0_base": g0_base,
            "g0_2": PHI * g0_base,
            "Ae": Ae, "Am": Am,
            "solutions": g0_3_solutions
        }
    else:
        print(f"    g0_3: NIE ZNALEZIONO w zakresie!")
        r31_arr = np.array([x for x in r31_vals if not np.isnan(x)])
        print(f"    r31 range: [{np.min(r31_arr):.1f}, {np.max(r31_arr):.1f}]")
        print(f"    r31_target = {r31_target:.1f}")
        results[label] = {"g0_base": g0_base, "Ae": Ae, "Am": Am, "solutions": []}

# ===== 2. Porownanie K miedzy sektorami =====
print(f"\n{'='*72}")
print("--- 2. Porownanie K miedzy sektorami ---")
print(f"{'='*72}")

print(f"\n  {'Sektor':20s} {'g0_1':>10} {'g0_2':>10} {'g0_3':>10} {'r21':>10} {'r31':>10} {'K':>10}")
print("  " + "-" * 80)
for label, sec in sectors.items():
    r = results.get(label)
    if r is None:
        print(f"  {label:20s} — nie znaleziono —")
        continue
    if r["solutions"]:
        g0_3, At, K, r31c = r["solutions"][-1]  # najwyzsze g0_3
        print(f"  {label:20s} {r['g0_base']:10.6f} {r['g0_2']:10.6f} {g0_3:10.6f} "
              f"{sec['r21']:10.2f} {r31c:10.1f} {K:10.6f}")
    else:
        print(f"  {label:20s} {r['g0_base']:10.6f} {r['g0_2']:10.6f} {'—':>10} "
              f"{sec['r21']:10.2f} {'—':>10} {'—':>10}")

# ===== 3. A_tail(g0) — profil universalny =====
print(f"\n--- 3. A_tail(g0) — profil universalny ---")
print("  Wszystkie sektory uzywaja TEGO SAMEGO ODE i TEJ SAMEJ funkcji A(g0).")
print("  Roznica jest TYLKO w wyborze g0 dla kazdej generacji.\n")

print(f"  {'g0':>8} {'A_tail':>12} {'A^4':>12} {'log10(A^4)':>12}")
print("  " + "-" * 50)
for g0 in [0.5, 0.6, 0.7, 0.8, 0.85, 0.87, 0.90, 0.95, 1.0,
           1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1]:
    A = A_of_g0(g0)
    print(f"  {g0:8.4f} {A:12.6f} {A**4:12.4f} {np.log10(max(A**4, 1e-30)):12.4f}")

# ===== 4. Test: co jesli kwarki maja inny phi? =====
print(f"\n--- 4. Inne proporcje g0_2/g0_1 niz phi ---")

for label, sec in sectors.items():
    r = results.get(label)
    if r is None or not r["solutions"]:
        continue

    g0_3, At, K_ode, r31c = r["solutions"][-1]
    g0_1 = r["g0_base"]
    g0_2 = r["g0_2"]

    # Ratios
    ratio_21 = g0_2 / g0_1
    ratio_31 = g0_3 / g0_1
    ratio_32 = g0_3 / g0_2

    print(f"\n  {label}:")
    print(f"    g0_1 = {g0_1:.8f}")
    print(f"    g0_2 = {g0_2:.8f} (= phi*g0_1 by construction)")
    print(f"    g0_3 = {g0_3:.8f} (from r31 PDG)")
    print(f"    g0_2/g0_1 = {ratio_21:.6f} (phi = {PHI:.6f})")
    print(f"    g0_3/g0_1 = {ratio_31:.6f}")
    print(f"    g0_3/g0_2 = {ratio_32:.6f}")

    # Test czy g0_3/g0_1 jest prosta funkcja phi
    print(f"    Testy g0_3/g0_1:")
    tests = {
        "2": 2.0,
        "phi^2": PHI**2,
        "phi+1": PHI+1,
        "3": 3.0,
        "2*phi": 2*PHI,
        "phi^(3/2)": PHI**1.5,
        "e": np.e,
        "pi/phi": np.pi/PHI,
        "sqrt(2*phi)": np.sqrt(2*PHI),
    }
    for tl, tv in tests.items():
        delta = abs(ratio_31 - tv) / tv * 100
        mark = " <<<" if delta < 3 else ""
        print(f"      vs {tl:12s} = {tv:.6f}  (delta = {delta:.2f}%){mark}")

# ===== 5. Test BEZ phi-FP: wolne g0_1, g0_2, g0_3 =====
print(f"\n--- 5. Bez phi-FP: wolne g0 dla kazdego sektora ---")
print("  Szukamy g0_1, g0_2, g0_3 takich ze:")
print("  (A2/A1)^4 = r21, (A3/A1)^4 = r31, BEZ warunku g0_2 = phi*g0_1\n")

for label, sec in sectors.items():
    r21_t = sec["r21"]
    r31_t = sec["r31"]

    # g0_1 jest wolne, g0_2 z r21, g0_3 z r31
    # Dla danego g0_1, szukamy g0_2 z (A(g0_2)/A(g0_1))^4 = r21
    # To wymaga odwrocenia A(g0) — zlozony problem

    # Prostsza metoda: skanuj g0_1, dla kazdego znajdz g0_2 i g0_3
    print(f"  {label}: r21={r21_t:.2f}, r31={r31_t:.2f}")

    # Precompute A(g0) na siatce
    g0_grid = np.linspace(0.3, 2.15, 120)
    A_grid = np.array([A_of_g0(g0) for g0 in g0_grid])

    # Dla kazdego g0_1, znajdz g0_2 i g0_3
    best_K = None
    best_combo = None

    for i1, g0_1 in enumerate(g0_grid):
        A1 = A_grid[i1]
        if A1 < 1e-10: continue

        A2_target = A1 * r21_t**0.25
        A3_target = A1 * r31_t**0.25

        # Znajdz g0_2 taki ze A(g0_2) ≈ A2_target
        idx2 = np.argmin(np.abs(A_grid - A2_target))
        if abs(A_grid[idx2] - A2_target) / A2_target > 0.02:
            continue

        # Znajdz g0_3 taki ze A(g0_3) ≈ A3_target
        idx3 = np.argmin(np.abs(A_grid - A3_target))
        if abs(A_grid[idx3] - A3_target) / A3_target > 0.02:
            continue

        g0_2 = g0_grid[idx2]
        g0_3 = g0_grid[idx3]
        A2 = A_grid[idx2]
        A3 = A_grid[idx3]

        r21_c = (A2/A1)**4
        r31_c = (A3/A1)**4
        K = koide_K_A(A1, A2, A3)

        if best_K is None or abs(r21_c/r21_t - 1) + abs(r31_c/r31_t - 1) < \
           abs(best_combo[3]/r21_t - 1) + abs(best_combo[4]/r31_t - 1):
            best_K = K
            best_combo = (g0_1, g0_2, g0_3, r21_c, r31_c, K,
                          g0_2/g0_1, g0_3/g0_1)

    if best_combo:
        g1, g2, g3, r21c, r31c, K, rat21, rat31 = best_combo
        print(f"    g0 = ({g1:.4f}, {g2:.4f}, {g3:.4f})")
        print(f"    r21 = {r21c:.2f} (target {r21_t:.2f})")
        print(f"    r31 = {r31c:.2f} (target {r31_t:.2f})")
        print(f"    K   = {K:.6f} (2/3 = 0.666667)")
        print(f"    g0_2/g0_1 = {rat21:.4f} (phi = {PHI:.4f}, delta = {abs(rat21 - PHI)/PHI*100:.1f}%)")
        print(f"    g0_3/g0_1 = {rat31:.4f}")

# ===== Wnioski =====
print(f"\n{'='*72}")
print("WNIOSKI ex160")
print(f"{'='*72}")
print(f"""
  STRUKTURA:
  - ODE substratowe daje universalna funkcje A_tail(g0)
  - phi-FP (g0_2 = phi*g0_1) reprodukuje r21 KAZDEGO sektora
  - Problem: g0_3 wyznaczone przez r31(PDG) daje K != 2/3 dla kwarkow

  PYTANIE CENTRALNE: co SELEKCJONUJE g0_3 w kazdym sektorze?
  - Leptony: K=2/3 selekcjonuje g0_3 (Koide)
  - Kwarki: K != 2/3, wiec INNY mechanizm selekcji

  MOZLIWE KIERUNKI:
  1. Kwarki uzywaja innego ODE (zmodyfikowane K_sub z powodu koloru)
  2. Kwarki uzywaja innej proporcji niz phi-FP
  3. Kwarki uzywaja shifted Koide (m + m_0)
  4. Masy kwarkowe sa scheme-dependent — moze na odpowiedniej skali K=2/3?
""")
