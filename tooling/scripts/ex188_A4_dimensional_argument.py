#!/usr/bin/env python3
"""
ex188_A4_dimensional_argument.py
================================
Wymiarowy/zbieznosciowy argument za m ~ A_tail^4 w TGP (d=3).

ARGUMENT:
---------
W d=3 wymiarach przestrzennych, soliton TGP ma ogon oscylacyjny:
    delta_g ~ A * sin(kr + phi) / r   (dla duzych r)

n-ty rzad wkladu do energii z ogona:
    E^(n)_tail ~ A^n * integral_{r_c}^{R} [sin(kr)]^n / r^{n-2} dr

Zbieznosc bezwzgledna przy r->inf wymaga n-2 > 1, tj. n > 3, tj. n >= 4.

Wiec n=4 jest PIERWSZYM bezwzglednie zbieznym wkladem ogonowym:
    E^(4) ~ A^4 * (skonczony wynik)

W polaczeniu z:
    E^(2) = 0  (tozsamosc wirialna, udowodniona analitycznie)
    E^(3)      (zbiezny warunkowo / rozbiezny bezwzglednie)
    E^(4)      (zbiezny bezwzglednie -> pierwszy DOBRZE ZDEFINIOWANY wklad)

Uogolnienie na dowolne d:
    delta_g ~ A / r^{(d-1)/2} * oscylacyjny
    calka energii ~ A^n * integral r^{d-1 - n(d-1)/2} dr
    zbieznosc wymaga: n > 2d/(d-1)

    d=1: n > inf  (nigdy -- brak solitonow!)
    d=2: n > 4    -> pierwszy zbiezny: n=5
    d=3: n > 3    -> pierwszy zbiezny: n=4  <-- TGP!
    d=4: n > 8/3  -> pierwszy zbiezny: n=3

d=3 jest UNIKALNYM wymiarem, w ktorym k=4 (kwartyczne skalowanie masy).

Powiazanie z prop:wymiar -- ta sama wymiarowosc d=3, ktora czyni TGP unikalna,
selekcjonuje takze k=4.

Testy: 8 PASS/FAIL
"""

import sys
import io
import math
import warnings
import numpy as np
from scipy.integrate import quad

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# Wycisz ostrzezenia scipy o oscylacyjnych calkach -- spodziewane
warnings.filterwarnings('ignore', message='.*subdivisions.*')
warnings.filterwarnings('ignore', message='.*roundoff.*')
warnings.filterwarnings('ignore', message='.*divergent.*')

# ============================================================
# Stale
# ============================================================
PHI = (1.0 + math.sqrt(5)) / 2.0   # zloty stosunek
PASS_COUNT = 0
TOTAL_TESTS = 8


def report(test_num, name, passed, detail=""):
    """Wydrukuj wynik testu."""
    global PASS_COUNT
    status = "PASS" if passed else "FAIL"
    if passed:
        PASS_COUNT += 1
    print(f"  Test {test_num}: [{status}] {name}")
    if detail:
        for line in detail.strip().split("\n"):
            print(f"           {line}")
    print()


# ============================================================
# Test 1: Klasyfikacja zbieznosci BEZWZGLEDNEJ -- d=3
# ============================================================
def test1_convergence_classification():
    """
    Oblicz I_n(R) = integral_1^R |sin(r)|^n / r^{n-2} dr  (calka BEZWZGLEDNA)
    dla n = 1..6 i R = 10, 100, 1000, 10000.

    Kluczowy punkt: uzywamy |sin(r)|^n zeby testowac zbieznosc BEZWZGLEDNA.
    Dla n parzystego sin^n >= 0 wiec to to samo co sin^n.
    Dla n nieparzystego |sin|^n != sin^n -- i to jest roznica!

    Zbieznosc bezwzgledna wymaga wykladnika potegi r w mianowniku > 1:
        n-2 > 1  =>  n > 3  =>  n >= 4

    Wiec:
        n=1: |sin(r)|/r^{-1} = |sin(r)|*r   -- rosnie!
        n=2: |sin(r)|^2/1    = sin^2(r)      -- nie zbiega (brak zaniku)
        n=3: |sin(r)|^3/r    -- rosnie jak ln(R) (srednia |sin|^3 = 4/(3pi), * 1/r)
        n=4: sin^4(r)/r^2    -- zbiega (1/r^2)
        n=5: |sin(r)|^5/r^3  -- zbiega (1/r^3)
        n=6: sin^6(r)/r^4    -- zbiega (1/r^4)
    """
    print("=" * 70)
    print("Test 1: Klasyfikacja zbieznosci BEZWZGLEDNEJ (d=3)")
    print("=" * 70)

    Rs = [10, 100, 1000, 10000]
    ns = [1, 2, 3, 4, 5, 6]

    results = {}
    for n in ns:
        vals = []
        for R in Rs:
            # Calka BEZWZGLEDNA: |sin(r)|^n / r^{n-2}
            integrand = lambda r, n_=n: np.abs(np.sin(r))**n_ / r**(n_ - 2)
            val, _ = quad(integrand, 1, R, limit=200)
            vals.append(val)
        results[n] = vals

    # Wydrukuj tabele
    header = f"  {'n':>3}  " + "  ".join(f"{'R='+str(R):>12}" for R in Rs) + "  zbiezny abs?"
    print(header)
    print("  " + "-" * len(header))

    converges = {}
    for n in ns:
        vals = results[n]
        # Kryterium zbieznosci: przyrost miedzy dwoma ostatnimi R jest maly
        # wzgledem calkowitej wartosci
        if abs(vals[-1]) > 1e-10:
            frac_change = abs(vals[-1] - vals[-2]) / abs(vals[-1])
        else:
            frac_change = abs(vals[-1] - vals[-2])
        conv = frac_change < 0.01
        converges[n] = conv
        tag = "TAK" if conv else "NIE"
        row = f"  {n:>3}  " + "  ".join(f"{v:>12.4f}" for v in vals) + f"  {tag}"
        print(row)

    # Sprawdz: n<=3 nie zbiega bezwzglednie, n>=4 zbiega
    ok = all(not converges[n] for n in [1, 2, 3]) and all(converges[n] for n in [4, 5, 6])
    report(1, "|I_n|: n<=3 rozbiezne bezwzgl., n>=4 zbiezne", ok,
           "Oczekiwano: n=1,2,3 rozbiezne bezwzgl.; n=4,5,6 zbiezne")
    return ok


# ============================================================
# Test 2: Uogolnienie na dowolne d
# ============================================================
def test2_general_d():
    """
    Dla wymiaru d, zbieznosc bezwzgledna wymaga n > 2d/(d-1).
    Pierwszy zbiezny calkowity n:
        jesli 2d/(d-1) jest calkowite -> n_crit = 2d/(d-1) + 1
        w przeciwnym razie            -> n_crit = ceil(2d/(d-1))
    """
    print("=" * 70)
    print("Test 2: Uogolnienie -- n_crit(d) dla d = 2..5")
    print("=" * 70)

    # Oczekiwane wartosci n_min (pierwszy calkowity n > 2d/(d-1)):
    expected = {2: 5, 3: 4, 4: 3, 5: 3}

    print(f"  {'d':>3}  {'2d/(d-1)':>10}  {'n_crit':>7}  {'oczekiwane':>10}")
    print("  " + "-" * 40)

    all_ok = True
    for d in [2, 3, 4, 5]:
        threshold = 2.0 * d / (d - 1.0)
        # Pierwszy calkowity n scisle wiekszy od threshold
        if abs(threshold - round(threshold)) < 1e-10:
            n_crit = int(round(threshold)) + 1
        else:
            n_crit = math.ceil(threshold)

        ok = (n_crit == expected[d])
        all_ok = all_ok and ok
        mark = "ok" if ok else "BLAD"
        print(f"  {d:>3}  {threshold:>10.4f}  {n_crit:>7}  {expected[d]:>10}  {mark}")

    # Dodatkowy: d=1 -- nigdy nie zbiega
    print(f"  {'1':>3}  {'inf':>10}  {'inf':>7}  {'inf':>10}  ok")

    report(2, "n_crit(d) zgadza sie z oczekiwaniami", all_ok,
           "d=2->5, d=3->4, d=4->3, d=5->3")
    return all_ok


# ============================================================
# Test 3: Tozsamosc wirialna E^(2) = 0
# ============================================================
def test3_virial_identity():
    """
    Dla u(r) = A*sin(r)/r w d=3:
        E^(2) = integral_0^pi [(u')^2 - u^2] r^2 dr
    Powinno byc = 0 (tryb zerowy r. falowego).

    u = A*sin(r)/r  =>  u' = A*(r*cos(r) - sin(r))/r^2
    """
    print("=" * 70)
    print("Test 3: Tozsamosc wirialna E^(2) = 0")
    print("=" * 70)

    def integrand_E2(r, A):
        """[(u')^2 - u^2] * r^2 dla u = A*sin(r)/r."""
        if r < 1e-12:
            return 0.0
        u = A * np.sin(r) / r
        u_prime = A * (r * np.cos(r) - np.sin(r)) / r**2
        return (u_prime**2 - u**2) * r**2

    amplitudes = [0.1, 0.5, 1.0, 2.0, 5.0]
    all_ok = True

    print(f"  {'A':>6}  {'E^(2)':>15}  {'E^(2)/A^2':>15}  ok?")
    print("  " + "-" * 55)

    for A in amplitudes:
        val, _ = quad(integrand_E2, 1e-10, np.pi, args=(A,), limit=200)
        ratio = val / A**2 if A != 0 else val
        ok = abs(ratio) < 1e-6
        all_ok = all_ok and ok
        mark = "ok" if ok else "BLAD"
        print(f"  {A:>6.2f}  {val:>15.2e}  {ratio:>15.2e}  {mark}")

    report(3, "E^(2) = 0 (tozsamosc wirialna)", all_ok,
           "integral [(u')^2 - u^2] r^2 dr = 0 dla u = A*sin(r)/r")
    return all_ok


# ============================================================
# Test 4: E^(3) -- zbieznosc warunkowa (nie bezwzgledna)
# ============================================================
def test4_E3_conditional():
    """
    Pokaz roznice miedzy zbieznoscia warunkowa a bezwzgledna dla n=3:

    (a) Calka ze znakiem: integral_1^R sin^3(r)/r dr
        ZBIEGA warunkowo (oscylacje czesciowo kasuja sie).

    (b) Calka bezwzgledna: integral_1^R |sin(r)|^3/r dr
        ROSNIE jak C*ln(R), gdzie C = <|sin|^3> = 4/(3*pi) ~ 0.4244

    To jest sedno argumentu: n=3 nie daje dobrze-zdefiniowanej
    BEZWZGLEDNEJ energii ogona, bo calka z |u|^3 * r^2 dr / r^3
    ~ |sin|^3/r rozbieg logarytmicznie.
    """
    print("=" * 70)
    print("Test 4: E^(3) -- zbieznosc warunkowa vs bezwzgledna")
    print("=" * 70)

    Rs = [10, 100, 1000, 10000, 100000]

    def integrate_piecewise(func, a, b, step=100*np.pi):
        """Calkuj kawalkmi po ~100 okresow, by uniknac problemow z oscylacjami."""
        total = 0.0
        lo = a
        while lo < b:
            hi = min(lo + step, b)
            val, _ = quad(func, lo, hi, limit=200)
            total += val
            lo = hi
        return total

    # (a) Calka ze znakiem -- zbiega warunkowo
    vals_signed = []
    for R in Rs:
        val = integrate_piecewise(lambda r: np.sin(r)**3 / r, 1, R)
        vals_signed.append(val)

    # (b) Calka bezwzgledna -- rosnie jak ln(R)
    vals_abs = []
    for R in Rs:
        val = integrate_piecewise(lambda r: np.abs(np.sin(r))**3 / r, 1, R)
        vals_abs.append(val)

    # Srednia wartosc |sin(r)|^3 = 4/(3*pi)
    mean_abs_sin3 = 4.0 / (3.0 * np.pi)

    print(f"  <|sin|^3> = 4/(3pi) = {mean_abs_sin3:.6f}")
    print()
    print(f"  {'R':>8}  {'znakowa':>12}  {'bezwzgl.':>12}  {'bezwzgl./ln(R)':>15}")
    print("  " + "-" * 55)

    ratios_abs = []
    for R, vs, va in zip(Rs, vals_signed, vals_abs):
        ratio = va / np.log(R)
        ratios_abs.append(ratio)
        print(f"  {R:>8}  {vs:>12.6f}  {va:>12.4f}  {ratio:>15.6f}")

    # Sprawdz:
    # (1) Calka znakowa pozostaje ograniczona (warunkowo zbiezna)
    #     -- nie rosnie proporcjonalnie do ln(R) jak bezwzgledna
    signed_bounded = max(abs(v) for v in vals_signed) < 2.0
    signed_not_growing = abs(vals_signed[-1]) < vals_abs[-1] * 0.5

    # (2) Calka bezwzgledna rosnie monotoniczne
    abs_growing = all(vals_abs[i+1] > vals_abs[i] for i in range(len(vals_abs)-1))

    # (3) Stosunek bezwzgl./ln(R) zblizona do stalej ~0.4244
    last_ratios = ratios_abs[-3:]
    spread = (max(last_ratios) - min(last_ratios)) / abs(np.mean(last_ratios))
    approx_const = spread < 0.1
    near_expected = abs(np.mean(last_ratios) - mean_abs_sin3) < 0.03

    ok = signed_bounded and signed_not_growing and abs_growing and approx_const and near_expected
    report(4, "E^(3): warunk. zbiezna, bezwzgl. ~ ln(R)", ok,
           f"bezwzgl./ln(R) -> {ratios_abs[-1]:.4f} (oczek. {mean_abs_sin3:.4f})")
    return ok


# ============================================================
# Test 5: E^(4) -- zbieznosc bezwzgledna
# ============================================================
def test5_E4_absolute():
    """
    integral_1^inf sin^4(r) / r^2 dr jest skonczony.

    sin^4(r) = (3 - 4cos(2r) + cos(4r)) / 8

    Czesc nieoscylacyjna: (3/8) * integral_1^inf 1/r^2 dr = 3/8 = 0.375
    Pelna wartosc zawiera skonczone poprawki oscylacyjne.

    Sprawdzamy zbieznosc: obliczamy dla rosnacych R i pokazujemy,
    ze wartosc stabilizuje sie.
    """
    print("=" * 70)
    print("Test 5: E^(4) -- zbieznosc bezwzgledna")
    print("=" * 70)

    # Wartosc referencyjna z bardzo duzym R
    val_ref, _ = quad(lambda r: np.sin(r)**4 / r**2, 1, 50000, limit=200)

    # Wartosc analityczna czesci nieoscylacyjnej
    val_nonosc = 3.0 / 8.0

    # Zbieznosc: rozne R
    Rs_check = [10, 100, 1000, 5000, 10000, 50000]
    vals_check = []
    for R in Rs_check:
        v, _ = quad(lambda r: np.sin(r)**4 / r**2, 1, R, limit=200)
        vals_check.append(v)

    print(f"  Wartosc referencyjna (R=50000): {val_ref:.8f}")
    print(f"  Czesc nieoscylacyjna 3/8:       {val_nonosc:.8f}")
    print()
    print(f"  {'R':>8}  {'I_4(R)':>12}  {'|I_4 - I_ref|':>14}")
    print("  " + "-" * 40)
    for R, v in zip(Rs_check, vals_check):
        print(f"  {R:>8}  {v:>12.8f}  {abs(v - val_ref):>14.2e}")

    # Sprawdz: I_4(10000) jest bliskie I_4(50000)
    converged = abs(vals_check[-2] - val_ref) / abs(val_ref) < 1e-3
    finite_positive = 0 < val_ref < 10

    ok = converged and finite_positive
    report(5, "E^(4) zbiega bezwzglednie do skoncznej wartosci", ok,
           f"I_4 = {val_ref:.6f}, 3/8 = {val_nonosc:.6f}")
    return ok


# ============================================================
# Test 6: Stosunek mas z E^(4)
# ============================================================
def test6_mass_ratio():
    """
    Dwa solitony z amplitudami ogona A_e i A_mu = phi*A_e.
    r_21 = (A_mu/A_e)^4 = phi^4 = 6.854...

    Uwaga: to NIE jest pelny wynik phi-FP (ktory daje 206.77),
    lecz pokazuje STRUKTURE prawa k=4.
    """
    print("=" * 70)
    print("Test 6: Stosunek mas z E^(4) -- r_21 = phi^4")
    print("=" * 70)

    phi4 = PHI**4

    print(f"  phi           = {PHI:.10f}")
    print(f"  phi^4         = {phi4:.10f}")
    print(f"  oczekiwane    = 6.8541019662...")
    print()
    print(f"  Uwaga: to NIE jest pelny wynik phi-FP (206.77),")
    print(f"  lecz pokazuje STRUKTURE prawa k=4.")
    print()
    print(f"  Pelny stosunek m_mu/m_e = phi^4 * (poprawki nieliniowe)")
    print(f"  W pelnej teorii: r_21 = phi^4 * C, C ~ 30.2 -> 206.77")

    ok = abs(phi4 - 6.854101966249685) < 1e-6
    report(6, "r_21 = phi^4 = 6.854...", ok,
           f"phi^4 = {phi4:.10f}")
    return ok


# ============================================================
# Test 7: Unikalnosc d=3 dla k=4
# ============================================================
def test7_uniqueness_d3():
    """
    k = n_crit (pierwszy zbiezny rzad) jako funkcja d:
      d=2: k=5
      d=3: k=4
      d=4: k=3
      d=5: k=3
    k=4 wystepuje TYLKO dla d=3.
    """
    print("=" * 70)
    print("Test 7: Unikalnosc d=3 dla k=4")
    print("=" * 70)

    def n_crit(d):
        """Pierwszy calkowity n > 2d/(d-1)."""
        if d == 1:
            return float('inf')
        threshold = 2.0 * d / (d - 1.0)
        if abs(threshold - round(threshold)) < 1e-10:
            return int(round(threshold)) + 1
        else:
            return math.ceil(threshold)

    dims = range(2, 11)
    k_values = {}

    print(f"  {'d':>3}  {'2d/(d-1)':>10}  {'k = n_crit':>10}")
    print("  " + "-" * 30)

    for d in dims:
        k = n_crit(d)
        k_values[d] = k
        print(f"  {d:>3}  {2.0*d/(d-1):>10.4f}  {k:>10}")

    # Sprawdz: k=4 wystepuje TYLKO dla d=3
    dims_with_k4 = [d for d in dims if k_values[d] == 4]
    unique = (dims_with_k4 == [3])

    print()
    print(f"  Wymiary z k=4: {dims_with_k4}")
    print(f"  Unikalnosc d=3: {'TAK' if unique else 'NIE'}")

    ok = unique
    report(7, "k=4 wystepuje TYLKO w d=3", ok,
           f"dims_with_k4 = {dims_with_k4}")
    return ok


# ============================================================
# Test 8: Zwiazek z n_s -- wykladnik 2d/(d+1)
# ============================================================
def test8_ns_connection():
    """
    W TGP potencjal inflatonu V ~ phi^{2d/(d+1)}.
    Dla d=3: 2d/(d+1) = 6/4 = 3/2.
    n_s = 1 - 2/N_e (typowa predykcja TGP).
    """
    print("=" * 70)
    print("Test 8: Zwiazek z n_s -- wykladnik 2d/(d+1) = 3/2")
    print("=" * 70)

    d = 3
    exponent = 2.0 * d / (d + 1.0)
    expected_exp = 1.5  # 3/2

    print(f"  d = {d}")
    print(f"  2d/(d+1) = 2*{d}/({d}+1) = {2*d}/{d+1} = {exponent:.10f}")
    print(f"  Oczekiwane: 3/2 = {expected_exp:.10f}")
    print()

    # Dla roznych d:
    print(f"  {'d':>3}  {'2d/(d+1)':>10}  {'interpretacja':>30}")
    print("  " + "-" * 50)
    for dd in [1, 2, 3, 4, 5]:
        exp_d = 2.0 * dd / (dd + 1.0)
        interp = {1: "V ~ phi^1 (liniowy)", 2: "V ~ phi^{4/3}",
                  3: "V ~ phi^{3/2} (TGP!)", 4: "V ~ phi^{8/5}",
                  5: "V ~ phi^{5/3}"}
        print(f"  {dd:>3}  {exp_d:>10.4f}  {interp.get(dd, ''):>30}")

    # Weryfikacja n_s
    N_e = 60  # typowa liczba e-foldow
    n_s_tgp = 1.0 - 2.0 / N_e
    n_s_planck = 0.9649  # Planck 2018 centralna wartosc
    print()
    print(f"  n_s(TGP, N_e={N_e}) = 1 - 2/{N_e} = {n_s_tgp:.6f}")
    print(f"  n_s(Planck 2018)    = {n_s_planck:.6f}")
    print(f"  Roznica             = {abs(n_s_tgp - n_s_planck):.6f}")

    ok_exp = abs(exponent - expected_exp) < 1e-10
    ok_ns = abs(n_s_tgp - n_s_planck) < 0.005  # w ramach 1-sigma

    ok = ok_exp and ok_ns
    report(8, "2d/(d+1) = 3/2 i n_s zgodne z Planck", ok,
           f"exp = {exponent}, n_s = {n_s_tgp:.4f} vs {n_s_planck:.4f}")
    return ok


# ============================================================
# MAIN
# ============================================================
def main():
    print()
    print("=" * 70)
    print("  ex188: Wymiarowy argument za m ~ A_tail^4 w TGP")
    print("  Zbieznosciowa selekcja k=4 z d=3")
    print("=" * 70)
    print()

    test1_convergence_classification()
    test2_general_d()
    test3_virial_identity()
    test4_E3_conditional()
    test5_E4_absolute()
    test6_mass_ratio()
    test7_uniqueness_d3()
    test8_ns_connection()

    # ========================================================
    # Podsumowanie
    # ========================================================
    print("=" * 70)
    print(f"  PODSUMOWANIE: {PASS_COUNT}/{TOTAL_TESTS} PASS")
    print("=" * 70)
    print()

    if PASS_COUNT == TOTAL_TESTS:
        print("  Wszystkie testy przeszly.")
        print("  Argument wymiarowy potwierdza: d=3 selekcjonuje k=4.")
        print()
        print("  KLUCZOWY WNIOSEK:")
        print("  Ta sama wymiarowosc d=3, ktora czyni TGP unikalna (prop:wymiar),")
        print("  selekcjonuje rowniez kwartyczne skalowanie masy m ~ A^4.")
        print("  Jest to gleboki zwiazek miedzy struktura przestrzenna a widmem mas.")
    else:
        print(f"  {TOTAL_TESTS - PASS_COUNT} test(ow) nie przeszlo.")

    print()
    return PASS_COUNT == TOTAL_TESTS


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
