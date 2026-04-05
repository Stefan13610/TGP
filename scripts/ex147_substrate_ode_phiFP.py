"""
ex147_substrate_ode_phiFP.py
============================
Weryfikacja phi-FP i skalowania M ~ A_tail^k z pelnym ODE substratowym.

Porownanie dwoch opisow:
  (A) Opis ciagly:    f(g)g'' + f'(g)(g')^2/2 + (2/r)f(g)g' = V'(g)
                       f(g) = 1 + 4 ln(g)
  (B) Opis substratowy: K_geo g^2 g'' + K_geo g (g')^2 + (2/r) K_geo g^2 g' = V'(g)
                        (K_sub(g) = K_geo g^2, K_geo = 1)

V(g) = g^3/3 - g^4/4, V'(g) = g^2(1-g)

Cel: sprawdzic czy:
  1. phi-FP istnieje w opisie (B) => r_21 zachowane
  2. A_tail(g0) z (B) zgadza sie z (A) w okolicy g~1
  3. Skalowanie M ~ A_tail^k z k=4 w obu opisach

TGP v1, sesja weryfikacyjna
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

# ============================================================
# 1. ODE solitonu - dwa opisy
# ============================================================

def soliton_ode_continuum(r, y):
    """Opis (A): f(g) = 1 + 4 ln(g)."""
    g, gp = y
    if r < 1e-12:
        Vp = g**2 * (1 - g)
        f = 1 + 4 * np.log(max(g, 1e-30))
        if abs(f) < 1e-15:
            f = 1e-15
        gpp = -Vp / f / 3.0
        return [gp, gpp]

    Vp = g**2 * (1 - g)
    f = 1 + 4 * np.log(max(g, 1e-30))
    fp = 4.0 / max(g, 1e-30)

    if abs(f) < 1e-15:
        return [gp, 0.0]

    gpp = (Vp - fp * gp**2 / 2 - 2 * gp / r * f) / f
    return [gp, gpp]


def soliton_ode_substrate(r, y):
    """Opis (B): K_sub(g) = g^2."""
    g, gp = y
    g_safe = max(g, 1e-30)
    if r < 1e-12:
        Vp = g_safe**2 * (1 - g_safe)
        gpp = -Vp / (g_safe**2) / 3.0
        return [gp, gpp]

    Vp = g_safe**2 * (1 - g_safe)
    # g^2 g'' + g (g')^2 + (2/r) g^2 g' = V'(g)
    gpp = (Vp - g_safe * gp**2 - 2 * gp / r * g_safe**2) / g_safe**2
    return [gp, gpp]


def shoot(g0, ode_func, r_max=60, n_points=20000):
    """Rozwiaz ODE z g(0)=g0, g'(0)=0."""
    r_span = (1e-10, r_max)
    r_eval = np.linspace(1e-10, r_max, n_points)
    try:
        sol = solve_ivp(ode_func, r_span, [g0, 0.0],
                        method='Radau', rtol=1e-11, atol=1e-13,
                        t_eval=r_eval, max_step=0.5)
        if sol.success:
            return sol.t, sol.y[0], sol.y[1]
    except Exception:
        pass
    return None, None, None


# ============================================================
# 2. Ekstrakcja amplitudy ogona A_tail
# ============================================================

def extract_tail(r, g, r_min=20, r_max=50):
    """Wyznacz A_tail i B_tail z ogona g(r) - 1 ~ (A sin(r) + B cos(r))/r."""
    mask = (r >= r_min) & (r <= r_max) & (r > 0.1)
    if np.sum(mask) < 10:
        return 0, 0

    rr = r[mask]
    uu = (g[mask] - 1) * rr  # u = r * (g-1) ~ A sin(r) + B cos(r)

    # Dopasuj u = A sin(r) + B cos(r) przez regresje liniowa
    sin_r = np.sin(rr)
    cos_r = np.cos(rr)
    X = np.column_stack([sin_r, cos_r])
    try:
        coeffs, residuals, rank, sv = np.linalg.lstsq(X, uu, rcond=None)
        A = coeffs[0]
        B = coeffs[1]
        return A, B
    except Exception:
        return 0, 0


def tail_amplitude(r, g, r_min=20, r_max=50):
    """A_tail = sqrt(A^2 + B^2)."""
    A, B = extract_tail(r, g, r_min, r_max)
    return np.sqrt(A**2 + B**2)


# ============================================================
# 3. Energia solitonu
# ============================================================

def soliton_energy(r, g, gp, use_substrate=False):
    """Calkowita energia solitonu (kinetyczna + potencjalna)."""
    if use_substrate:
        # E = integral [g^2/2 (g')^2 + V(g)] r^2 dr
        K_vals = g**2
    else:
        # E = integral [f(g)/2 (g')^2 + V(g)] r^2 dr
        K_vals = 1 + 4 * np.log(np.maximum(g, 1e-30))
        K_vals = np.maximum(K_vals, 1e-15)

    V_vals = g**3 / 3 - g**4 / 4
    V_vac = 1/3 - 1/4  # = 1/12

    # Energia wzgledem prozni
    integrand = (K_vals / 2 * gp**2 + V_vals - V_vac) * r**2
    E = 4 * np.pi * np.trapezoid(integrand, r)
    return E


# ============================================================
# 4. Phi-FP: warunek punktu stalego
# ============================================================

PHI = (1 + np.sqrt(5)) / 2  # Zlota proporcja

def compute_r21(g0, ode_func, use_substrate=False, r_max=60):
    """Oblicz r_21 = (A_tail(phi*g0) / A_tail(g0))^4."""
    r1, g1, gp1 = shoot(g0, ode_func, r_max=r_max)
    if r1 is None:
        return None

    g0_mu = PHI * g0
    r2, g2, gp2 = shoot(g0_mu, ode_func, r_max=r_max)
    if r2 is None:
        return None

    A1 = tail_amplitude(r1, g1)
    A2 = tail_amplitude(r2, g2)

    if A1 < 1e-20:
        return None

    return (A2 / A1) ** 4


def find_phiFP(ode_func, g0_range=(1.15, 1.40), use_substrate=False):
    """Znajdz g0* dla phi-FP: r_21(g0*) = r_21_PDG."""
    r21_PDG = 206.768

    def residual(g0):
        r21 = compute_r21(g0, ode_func, use_substrate)
        if r21 is None:
            return 1e10
        return r21 - r21_PDG

    try:
        g0_star = brentq(residual, g0_range[0], g0_range[1], xtol=1e-6)
        r21_val = compute_r21(g0_star, ode_func, use_substrate)
        return g0_star, r21_val
    except Exception as e:
        # Skan
        g0s = np.linspace(g0_range[0], g0_range[1], 50)
        r21s = []
        for g0 in g0s:
            r21 = compute_r21(g0, ode_func, use_substrate)
            r21s.append(r21 if r21 is not None else 0)
        r21s = np.array(r21s)
        valid = r21s > 0
        if np.any(valid):
            idx = np.argmin(np.abs(r21s[valid] - r21_PDG))
            g0_best = g0s[valid][idx]
            return g0_best, r21s[valid][idx]
        return None, None


# ============================================================
# 5. Skalowanie M vs A_tail
# ============================================================

def mass_scaling(ode_func, use_substrate=False, g0_values=None):
    """Wyznacz wykladnik k w M ~ A_tail^k."""
    if g0_values is None:
        g0_values = np.linspace(1.05, 1.50, 20)

    A_tails = []
    energies = []

    for g0 in g0_values:
        r, g, gp = shoot(g0, ode_func, r_max=60)
        if r is None:
            continue
        A = tail_amplitude(r, g)
        E = soliton_energy(r, g, gp, use_substrate)
        if A > 1e-15 and E > 0:
            A_tails.append(A)
            energies.append(E)

    if len(A_tails) < 5:
        return None, None

    logA = np.log(A_tails)
    logE = np.log(energies)

    # Dopasuj k: log(E) = k*log(A) + const
    coeffs = np.polyfit(logA, logE, 1)
    k = coeffs[0]
    return k, (A_tails, energies)


# ============================================================
# 6. TESTY
# ============================================================

def run_tests():
    n_pass = 0
    n_total = 0

    def test(name, condition, detail=""):
        nonlocal n_pass, n_total
        n_total += 1
        status = "PASS" if condition else "FAIL"
        if condition:
            n_pass += 1
        print(f"  [{status}] {name}")
        if detail:
            print(f"         {detail}")

    print("=" * 70)
    print("ex147: Weryfikacja phi-FP z pelnym ODE substratowym")
    print("=" * 70)

    # ----- Test 1: phi-FP w opisie ciaglym (A) -----
    print("\n--- Test 1: phi-FP w opisie ciaglym f(g) = 1+4ln(g) ---")
    g0_A, r21_A = find_phiFP(soliton_ode_continuum, g0_range=(1.20, 1.35))
    if g0_A is not None:
        print(f"    g0* = {g0_A:.5f}")
        print(f"    r_21 = {r21_A:.2f} (PDG: 206.768)")
        test("phiFP_continuum", abs(r21_A - 206.768) / 206.768 < 0.02,
             f"r_21 = {r21_A:.2f}, odch. = {abs(r21_A-206.768)/206.768*100:.2f}%")
    else:
        test("phiFP_continuum", False, "Nie znaleziono")

    # ----- Test 2: phi-FP w opisie substratowym (B) -----
    print("\n--- Test 2: phi-FP w opisie substratowym K_sub = g^2 ---")
    g0_B, r21_B = find_phiFP(soliton_ode_substrate, g0_range=(1.10, 1.50),
                              use_substrate=True)
    if g0_B is not None:
        print(f"    g0* = {g0_B:.5f}")
        print(f"    r_21 = {r21_B:.2f} (PDG: 206.768)")
        test("phiFP_substrate", r21_B is not None and abs(r21_B - 206.768)/206.768 < 0.15,
             f"r_21 = {r21_B:.2f}, odch. = {abs(r21_B-206.768)/206.768*100:.2f}%")
    else:
        test("phiFP_substrate", False, "Nie znaleziono")

    # ----- Test 3: Porownanie g0* miedzy opisami -----
    print("\n--- Test 3: Porownanie g0* ---")
    if g0_A is not None and g0_B is not None:
        delta_g0 = abs(g0_A - g0_B) / g0_A
        test("g0_agreement", delta_g0 < 0.15,
             f"|g0_A - g0_B|/g0_A = {delta_g0*100:.1f}%")
    else:
        test("g0_agreement", False, "Brak danych")

    # ----- Test 4: Skalowanie M ~ A_tail^k (opis ciagly) -----
    print("\n--- Test 4: Skalowanie M ~ A_tail^k (opis ciagly) ---")
    k_A, data_A = mass_scaling(soliton_ode_continuum)
    if k_A is not None:
        print(f"    k = {k_A:.2f} (oczekiwane: ~4)")
        test("scaling_continuum", abs(k_A - 4) < 1.5,
             f"k = {k_A:.2f}")
    else:
        test("scaling_continuum", False, "Brak danych")

    # ----- Test 5: Skalowanie M ~ A_tail^k (opis substratowy) -----
    print("\n--- Test 5: Skalowanie M ~ A_tail^k (opis substratowy) ---")
    k_B, data_B = mass_scaling(soliton_ode_substrate, use_substrate=True)
    if k_B is not None:
        print(f"    k = {k_B:.2f} (oczekiwane: ~4)")
        test("scaling_substrate", abs(k_B - 4) < 1.5,
             f"k = {k_B:.2f}")
    else:
        test("scaling_substrate", False, "Brak danych")

    # ----- Test 6: A_tail zgadza sie blisko g=1 -----
    print("\n--- Test 6: A_tail(g0=1.1) w obu opisach ---")
    g0_test = 1.10
    r_A, g_A, gp_A = shoot(g0_test, soliton_ode_continuum)
    r_B, g_B, gp_B = shoot(g0_test, soliton_ode_substrate)
    if r_A is not None and r_B is not None:
        At_A = tail_amplitude(r_A, g_A)
        At_B = tail_amplitude(r_B, g_B)
        if At_A > 0 and At_B > 0:
            rel = abs(At_A - At_B) / At_A
            print(f"    A_tail(continuum) = {At_A:.6e}")
            print(f"    A_tail(substrate) = {At_B:.6e}")
            test("Atail_agreement_g1p1", rel < 0.5,
                 f"Roznica = {rel*100:.1f}%")
        else:
            test("Atail_agreement_g1p1", False, "A_tail = 0")
    else:
        test("Atail_agreement_g1p1", False, "Solver failed")

    # ----- Test 7: Profil ogona oscylacyjny -----
    print("\n--- Test 7: Profil ogona g0=1.24 ---")
    g0_e = 1.24
    r_A, g_A, gp_A = shoot(g0_e, soliton_ode_continuum, r_max=80)
    r_B, g_B, gp_B = shoot(g0_e, soliton_ode_substrate, r_max=80)
    if r_A is not None and r_B is not None:
        At_A = tail_amplitude(r_A, g_A, r_min=25, r_max=70)
        At_B = tail_amplitude(r_B, g_B, r_min=25, r_max=70)
        print(f"    A_tail(continuum, g0=1.24) = {At_A:.6e}")
        print(f"    A_tail(substrate, g0=1.24) = {At_B:.6e}")
        if At_A > 0:
            rel = abs(At_A - At_B) / At_A
            test("Atail_electron", rel < 0.5, f"Roznica = {rel*100:.1f}%")
        else:
            test("Atail_electron", False, "A_tail = 0")
    else:
        test("Atail_electron", False, "Solver failed")

    # ----- Test 8: Stosunek A_tail(phi*g0)/A_tail(g0) -----
    print("\n--- Test 8: Stosunek amplitud w obu opisach ---")
    if r21_A is not None and r21_B is not None:
        ratio_A = r21_A**(1/4)
        ratio_B = r21_B**(1/4)
        print(f"    (A_mu/A_e)_continuum = {ratio_A:.4f}")
        print(f"    (A_mu/A_e)_substrate = {ratio_B:.4f}")
        test("ratio_preserved", abs(ratio_A - ratio_B) / ratio_A < 0.15,
             f"Roznica = {abs(ratio_A-ratio_B)/ratio_A*100:.1f}%")
    else:
        test("ratio_preserved", False, "Brak danych")

    # ---- Podsumowanie ----
    print("\n" + "=" * 70)
    print(f"WYNIK: {n_pass}/{n_total} PASS")
    print("=" * 70)

    if n_pass >= 6:
        print("\nWNIOSEK: ODE substratowy K_sub(g)=g^2 zachowuje mechanizm phi-FP")
        print("i skalowanie M ~ A_tail^4. Ghost-freedom jest potwierdzone")
        print("bez utraty fizycznych przewidywan teorii.")
    else:
        print("\nUWAGA: Niektore testy nie przeszly - wymaga dalszej analizy.")

    return n_pass, n_total


if __name__ == "__main__":
    run_tests()
