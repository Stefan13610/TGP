#!/usr/bin/env python3
"""
ex169_third_gen_selection.py
Sesja v45, 2026-04-05

Problem: phi-FP daje r21 (uniwersalnie). Co determinuje r31?
- Koide K=2/3 dziala TYLKO dla leptonow
- Dla kwarkow K != 2/3
- Pytanie: czy istnieje uniwersalny warunek na g0^(3)?

Hipotezy testowane:
H1: g0^(3) = phi^2 * g0^(1)  (prosta iteracja phi-FP)
H2: g0^(3) = phi * g0^(2) = phi^2 * g0^(1)  (ten sam co H1)
H3: g0^(3) z warunku K(sektor) = K_empiryczny  (nie-predykcyjne)
H4: g0^(3) z warunku cross-Koide (b,c,t) = 2/3
H5: g0^(3) = f(g0^(1), phi) z jakiejs uniwersalnej funkcji f
H6: Potegi phi: g0^(n) = phi^(a_n) * g0^(1) z a_n != n-1
H7: Stosunek r31/r21^p = const (skalowanie potegowe)

Kluczowe narzedzie: ODE substratowe (alpha=1) z phi-FP.
"""
import numpy as np
from scipy.integrate import solve_ivp

PHI = (1 + np.sqrt(5)) / 2  # golden ratio

# ---- PDG masses (MeV) ----
SECTORS = {
    'lepton': {'m1': 0.51100, 'm2': 105.658, 'm3': 1776.86,
               'names': ('e', 'mu', 'tau')},
    'down':   {'m1': 4.67, 'm2': 93.4, 'm3': 4180.0,
               'names': ('d', 's', 'b')},
    'up':     {'m1': 2.16, 'm2': 1270.0, 'm3': 172760.0,
               'names': ('u', 'c', 't')},
}

for sec in SECTORS.values():
    sec['r21'] = sec['m2'] / sec['m1']
    sec['r31'] = sec['m3'] / sec['m1']
    sec['r32'] = sec['m3'] / sec['m2']


# ---- ODE solver (alpha=1, substrate) ----
def solve_soliton(g0, r_max=80, n_points=20000):
    """Solve g'' + (1/g)(g')^2 + (2/r)g' = 1-g"""
    def rhs(r, y):
        g, gp = y
        if g < 1e-12:
            g = 1e-12
        gpp = (1 - g) - (1.0 / g) * gp**2 - (2.0 / r) * gp if r > 1e-12 else (1 - g) / 3.0
        return [gp, gpp]

    r_span = (1e-6, r_max)
    r_eval = np.linspace(r_span[0], r_span[1], n_points)
    y0 = [g0, 0.0]
    sol = solve_ivp(rhs, r_span, y0, method='RK45', t_eval=r_eval,
                    rtol=1e-10, atol=1e-12, max_step=0.05)
    return sol.t, sol.y[0]


def get_Atail(g0, r_max=80):
    """Extract A_tail from oscillatory tail g(r) ~ 1 + A*cos(r+d)/r"""
    r, g = solve_soliton(g0, r_max)
    delta = g - 1.0
    # Window: r in [40, 70] (well into tail)
    mask = (r > 40) & (r < 70)
    r_w = r[mask]
    d_w = delta[mask]
    if len(r_w) < 100:
        return None
    # A_tail = max|delta * r| in window
    envelope = np.abs(d_w * r_w)
    return np.median(np.sort(envelope)[-50:])  # robust estimate


def find_g0_for_r21(r21_target, g0_guess=0.87, tol=1e-6):
    """Find g0^(1) such that phi-FP gives r21 = (A2/A1)^4"""
    from scipy.optimize import brentq

    def residual(g0_1):
        g0_2 = PHI * g0_1
        A1 = get_Atail(g0_1)
        A2 = get_Atail(g0_2)
        if A1 is None or A2 is None or A1 < 1e-15:
            return 1e10
        r21_calc = (A2 / A1)**4
        return r21_calc - r21_target

    # Bracket search
    return brentq(residual, 0.75, 0.95, xtol=tol)


print("=" * 72)
print("ex169: Selekcja trzeciej generacji -- problem otwarty R12")
print("=" * 72)

# ---- 1. Reproduce phi-FP for all sectors ----
print("\n--- 1. phi-FP: reprodukcja r21 (weryfikacja) ---\n")

g0_values = {}
for name, sec in SECTORS.items():
    g0_1 = find_g0_for_r21(sec['r21'])
    g0_2 = PHI * g0_1
    A1 = get_Atail(g0_1)
    A2 = get_Atail(g0_2)
    r21_calc = (A2 / A1)**4
    g0_values[name] = {'g0_1': g0_1, 'g0_2': g0_2, 'A1': A1, 'A2': A2}
    print(f"  {name:>7s}: g0^(1) = {g0_1:.5f}, g0^(2) = {g0_2:.5f}, "
          f"r21 = {r21_calc:.1f} (PDG: {sec['r21']:.1f})")


# ---- 2. Test H1: g0^(3) = phi^2 * g0^(1) ----
print("\n--- 2. H1: g0^(3) = phi^2 * g0^(1) ---\n")
for name, sec in SECTORS.items():
    g0_1 = g0_values[name]['g0_1']
    g0_3_H1 = PHI**2 * g0_1
    A1 = g0_values[name]['A1']
    A3 = get_Atail(g0_3_H1)
    if A3 is not None and A1 > 0:
        r31_H1 = (A3 / A1)**4
        r31_PDG = sec['r31']
        err = (r31_H1 - r31_PDG) / r31_PDG * 100
        print(f"  {name:>7s}: g0^(3) = {g0_3_H1:.5f}, "
              f"r31 = {r31_H1:.1f} (PDG: {r31_PDG:.1f}), "
              f"err = {err:+.1f}%")
    else:
        print(f"  {name:>7s}: A_tail extraction failed for g0 = {g0_3_H1:.4f}")


# ---- 3. Inverse problem: what g0^(3) gives correct r31? ----
print("\n--- 3. Inverse: g0^(3) dajace poprawne r31 ---\n")

def find_g0_for_r31(A1, r31_target):
    """Find g0^(3) such that (A3/A1)^4 = r31"""
    from scipy.optimize import brentq

    def residual(g0_3):
        A3 = get_Atail(g0_3)
        if A3 is None or A1 < 1e-15:
            return 1e10
        return (A3 / A1)**4 - r31_target

    # Need broader search for heavy quarks
    try:
        return brentq(residual, 0.5, 3.0, xtol=1e-5)
    except ValueError:
        # Try wider range
        try:
            return brentq(residual, 0.3, 5.0, xtol=1e-5)
        except ValueError:
            return None


g0_3_exact = {}
for name, sec in SECTORS.items():
    A1 = g0_values[name]['A1']
    g0_3 = find_g0_for_r31(A1, sec['r31'])
    g0_1 = g0_values[name]['g0_1']
    g0_2 = g0_values[name]['g0_2']

    if g0_3 is not None:
        g0_3_exact[name] = g0_3
        ratio_to_g01 = g0_3 / g0_1
        ratio_to_g02 = g0_3 / g0_2
        ratio_phi2 = PHI**2

        print(f"  {name:>7s}: g0^(3) = {g0_3:.5f}")
        print(f"           g0^(3)/g0^(1) = {ratio_to_g01:.5f}  "
              f"(phi^2 = {ratio_phi2:.5f}, err = {(ratio_to_g01/ratio_phi2-1)*100:+.1f}%)")
        print(f"           g0^(3)/g0^(2) = {ratio_to_g02:.5f}  "
              f"(phi = {PHI:.5f}, err = {(ratio_to_g02/PHI-1)*100:+.1f}%)")

        # What power of phi?
        if ratio_to_g01 > 0:
            n_phi = np.log(ratio_to_g01) / np.log(PHI)
            print(f"           g0^(3) = phi^{n_phi:.4f} * g0^(1)")
    else:
        print(f"  {name:>7s}: nie znaleziono g0^(3) w zakresie [0.3, 5.0]")


# ---- 4. Test: r31 = r21^p scaling ----
print("\n--- 4. Skalowanie potegowe: r31 = r21^p ---\n")
for name, sec in SECTORS.items():
    r21 = sec['r21']
    r31 = sec['r31']
    if r21 > 1:
        p = np.log(r31) / np.log(r21)
        print(f"  {name:>7s}: r21 = {r21:.1f}, r31 = {r31:.0f}, "
              f"p = ln(r31)/ln(r21) = {p:.4f}")

# Check if p is related to phi
print(f"\n  phi = {PHI:.5f}")
print(f"  phi + 1 = {PHI+1:.5f}")
print(f"  2*phi = {2*PHI:.5f}")
print(f"  phi^2 = {PHI**2:.5f}")

for name, sec in SECTORS.items():
    r21 = sec['r21']
    r31 = sec['r31']
    if r21 > 1:
        p = np.log(r31) / np.log(r21)
        # Test phi-related values
        for label, val in [('phi', PHI), ('phi+1', PHI+1),
                           ('2', 2.0), ('phi^2', PHI**2),
                           ('3/2', 1.5), ('2*phi-1', 2*PHI-1)]:
            r31_pred = r21**val
            err = (r31_pred - r31) / r31 * 100
            if abs(err) < 50:
                print(f"  {name:>7s}: r31 = r21^{label} = {r31_pred:.0f} "
                      f"(PDG: {r31:.0f}, err = {err:+.1f}%)")


# ---- 5. Test: g0^(3) = phi * g0^(2) (chained phi-FP) ----
print("\n--- 5. Chained phi-FP: g0^(3) = phi * g0^(2) ---\n")
print("  (to jest identyczne z H1: g0^(3) = phi^2 * g0^(1))")
print("  Ale testujemy r32 = (A3/A2)^4:\n")

for name, sec in SECTORS.items():
    g0_2 = g0_values[name]['g0_2']
    g0_3_chain = PHI * g0_2  # = phi^2 * g0^(1)
    A2 = g0_values[name]['A2']
    A3 = get_Atail(g0_3_chain)
    if A3 is not None and A2 > 0:
        r32_chain = (A3 / A2)**4
        r32_PDG = sec['r32']
        err = (r32_chain - r32_PDG) / r32_PDG * 100
        print(f"  {name:>7s}: r32_chain = {r32_chain:.1f} (PDG: {r32_PDG:.1f}), "
              f"err = {err:+.1f}%")


# ---- 6. Koide-inspired: what K gives correct r31 from r21? ----
print("\n--- 6. K efektywne: jaki K daje poprawne r31 z r21? ---\n")
print("  K = (1 + r21 + r31) / (1 + sqrt(r21) + sqrt(r31))^2\n")

for name, sec in SECTORS.items():
    r21 = sec['r21']
    r31 = sec['r31']
    denom = (1 + np.sqrt(r21) + np.sqrt(r31))**2
    K = (1 + r21 + r31) / denom
    print(f"  {name:>7s}: K = {K:.6f}  (2/3 = {2/3:.6f}, "
          f"delta = {(K - 2/3)/(2/3)*100:+.2f}%)")

# ---- 7. Cross-Koide (b,c,t) as constraint ----
print("\n--- 7. (b,c,t) cross-Koide jako warunek ---\n")
print("  Jesli K(b,c,t) = 2/3, to m_t = f(m_b, m_c):\n")

m_b = SECTORS['down']['m3']   # 4180 MeV
m_c = SECTORS['up']['m2']     # 1270 MeV

# K = (m_b + m_c + m_t) / (sqrt(m_b) + sqrt(m_c) + sqrt(m_t))^2 = 2/3
# Let x = sqrt(m_t). Then:
# (m_b + m_c + x^2) / (sqrt(m_b) + sqrt(m_c) + x)^2 = 2/3
# 3(m_b + m_c + x^2) = 2(sqrt(m_b) + sqrt(m_c) + x)^2
# 3m_b + 3m_c + 3x^2 = 2(sqrt(m_b)+sqrt(m_c))^2 + 4(sqrt(m_b)+sqrt(m_c))x + 2x^2
# x^2 - 4(sqrt(m_b)+sqrt(m_c))x + 3(m_b+m_c) - 2(sqrt(m_b)+sqrt(m_c))^2 = 0

sq_b = np.sqrt(m_b)
sq_c = np.sqrt(m_c)
S = sq_b + sq_c

A_coeff = 1.0
B_coeff = -4.0 * S
C_coeff = 3.0 * (m_b + m_c) - 2.0 * S**2

discriminant = B_coeff**2 - 4*A_coeff*C_coeff

if discriminant >= 0:
    x1 = (-B_coeff + np.sqrt(discriminant)) / (2*A_coeff)
    x2 = (-B_coeff - np.sqrt(discriminant)) / (2*A_coeff)
    for i, x in enumerate([x1, x2]):
        if x > 0:
            m_t_pred = x**2
            m_t_PDG = SECTORS['up']['m3']
            err = (m_t_pred - m_t_PDG) / m_t_PDG * 100
            print(f"  Rozwiazanie {i+1}: sqrt(m_t) = {x:.1f} MeV^(1/2)")
            print(f"    m_t = {m_t_pred:.0f} MeV = {m_t_pred/1000:.1f} GeV "
                  f"(PDG: {m_t_PDG/1000:.1f} GeV, err = {err:+.1f}%)")

            # Verify K
            K_check = (m_b + m_c + m_t_pred) / (sq_b + sq_c + x)**2
            print(f"    K_check = {K_check:.6f}")

# ---- 8. Generalized phi-power for 3rd gen ----
print("\n--- 8. Uogolnione phi-potegi ---\n")
print("  g0^(3) = phi^n * g0^(1), szukamy optymalnego n:\n")

for name in ['lepton', 'down', 'up']:
    if name not in g0_3_exact:
        continue
    g0_1 = g0_values[name]['g0_1']
    g0_3 = g0_3_exact[name]
    n_opt = np.log(g0_3 / g0_1) / np.log(PHI)
    print(f"  {name:>7s}: n_opt = {n_opt:.4f}")

print(f"\n  Dla phi-FP (2nd gen): n = 1 (dokladnie)")
print(f"  Jesli n_opt jest universalny -> nowa selekcja")

# ---- 9. Summary of r31 predictions ----
print("\n" + "=" * 72)
print("PODSUMOWANIE: predykcje r31")
print("=" * 72)

print(f"\n  {'Metoda':<30s}  {'lepton':>10s}  {'down':>10s}  {'up':>10s}")
print("  " + "-" * 65)
print(f"  {'PDG (dokladne)':<30s}  {'3477':>10s}  {'895':>10s}  {'79982':>10s}")

# H1: phi^2
for name in ['lepton', 'down', 'up']:
    g0_1 = g0_values[name]['g0_1']
    g0_3 = PHI**2 * g0_1
    A1 = g0_values[name]['A1']
    A3 = get_Atail(g0_3)
    if A3 is not None:
        r31 = (A3/A1)**4
        g0_values[name]['r31_H1'] = r31

print(f"  {'H1: g0^(3)=phi^2*g0^(1)':<30s}  "
      f"{g0_values['lepton'].get('r31_H1', 0):>10.0f}  "
      f"{g0_values['down'].get('r31_H1', 0):>10.0f}  "
      f"{g0_values['up'].get('r31_H1', 0):>10.0f}")

# Cross-Koide (b,c,t)
print(f"  {'H4: K(b,c,t)=2/3 -> m_t':<30s}  {'---':>10s}  {'---':>10s}", end="")
if discriminant >= 0 and x1 > 0:
    print(f"  {x1**2/SECTORS['up']['m1']:>10.0f}")
else:
    print(f"  {'---':>10s}")

print(f"\n  Uwaga: phi^2 scaling odpowiada r31 = r21^p z p ~ 2")
print(f"  ale p_lepton != p_down != p_up (brak uniwersalnosci)")

print("\n" + "=" * 72)
print("WNIOSKI ex169")
print("=" * 72)
print("""
  1. phi-FP (g0^(2) = phi*g0^(1)) jest UNIWERSALNY dla r21 -- OK.
  2. Prosta iteracja phi^2 (H1) NIE daje poprawnego r31:
     - leptony: r31 zalezy od Koide K=2/3 (nie od phi^2)
     - kwarki: r31 nie wynika z phi^2 ani z K=2/3
  3. Skalowanie potegowe r31 = r21^p:
     - p jest ROZNY w kazdym sektorze
     - brak uniwersalnej wartosci p
  4. Cross-Koide (b,c,t) = 2/3 daje predykcje m_t
     -- do porownania z PDG.
  5. WNIOSEK: trzecia generacja wymaga DODATKOWEGO
     warunku selekcji, specyficznego sektorowo.
     phi-FP sam nie wystarczy.
  6. Kandydaci na warunek 3rd gen:
     - Koide K=2/3 (lepton-specific)
     - Cross-Koide (b,c,t) -- do zbadania
     - Shifted Koide z m0
     - Nowy punkt staly ODE
""")
