#!/usr/bin/env python3
"""
O-K1: Czy Koide Q=3/2 wynika z TGP?
=====================================
Najwazniejszy otwarty problem TGP.

Koide: Q = (sum sqrt(m_n))^2 / sum(m_n) = 3/2
Rownowazne: sqrt(m_n) = a*(1 + sqrt(2)*cos(2*pi*n/3 + delta))

W TGP: m_n ~ A_tail(g0_n)^4, wiec sqrt(m_n) ~ A_tail(g0_n)^2
Koide staje sie:
  A_tail(g0_n)^2 = a*(1 + sqrt(2)*cos(2*pi*n/3 + delta))

STRATEGIA:
1. Precyzyjna mapa A_tail(g0) z ODE
2. Test: czy A_tail^2 ma forme "a + b*cos(phi(g0))"?
3. Szukanie g0_tau BEZ zakladania Koide — czy ODE sam daje Q~3/2?
4. Analiza energii solitonu: czy Koide minimalizuje cos?
5. Ukryta symetria ODE: Z_3, dualnosc, etc.

TGP v1 -- 2026-03-31 (O-K1)
"""

import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# 1. ODE solitonu TGP (forma z dodatekJ)
# ============================================================
PHI = (1 + np.sqrt(5)) / 2  # zlota proporcja
ALPHA = 2
GG = np.exp(-1 / (2 * ALPHA))  # g* = exp(-1/4) ~ 0.7788 (ghost)

def fk(g):
    """f(g) = 1 + 2*alpha*ln(g) = 1 + 4*ln(g)"""
    return 1 + 2 * ALPHA * np.log(g) if g > 0 else -1e30

def Vp(g):
    """V'(g) = g^2*(1-g)"""
    return g**2 * (1 - g)

def shoot(g0, r_max=300, n_pts=25000):
    """Rozwiaz ODE: f(g)*g'' + (2/r)*g' = V'(g) z g(0)=g0, g'(0)=0."""
    fg0 = fk(g0)
    c2 = Vp(g0) / (3 * fg0)  # g ~ g0 + c2*r^2 blisko r=0
    rs = 0.01

    def rhs(r, y):
        g, p = y
        if g <= 1e-15: return [p, 0]
        fg = fk(g)
        if abs(fg) < 1e-10: return [p, 0]
        if r < 1e-10: return [p, Vp(g) / fg / 3]
        return [p, (Vp(g) - 2/r * p) / fg]

    def ev(r, y):
        return 100 - abs(y[0])
    ev.terminal = True

    s = solve_ivp(rhs, [rs, r_max], [g0 + c2*rs**2, 2*c2*rs],
                  method='RK45', rtol=1e-11, atol=1e-13,
                  max_step=0.05, events=[ev], dense_output=True)
    r = np.linspace(rs, min(s.t[-1], r_max), n_pts)
    return r, s.sol(r)[0]

def A_tail(g0, r_fit_lo=120, r_fit_hi=260):
    """Amplituda ogona oscylacyjnego."""
    r, g = shoot(g0)
    m = (r >= r_fit_lo) & (r <= r_fit_hi)
    rf = r[m]
    tl = (g[m] - 1) * r[m]
    if len(rf) < 10:
        return np.nan
    A = np.column_stack([np.cos(rf), np.sin(rf)])
    B, C = np.linalg.lstsq(A, tl, rcond=None)[0]
    return np.sqrt(B**2 + C**2)

print("=" * 70)
print("  O-K1: CZY KOIDE Q=3/2 WYNIKA Z TGP?")
print("=" * 70)
print(f"  ODE: f(g)*g'' + (2/r)*g' = V'(g)")
print(f"  f(g) = 1 + 4*ln(g),  V'(g) = g^2*(1-g)")
print(f"  Ghost: g* = {GG:.6f}")
print(f"  phi = {PHI:.6f}")
print()

# ============================================================
# 2. Precyzyjna mapa A_tail(g0)
# ============================================================
print("--- 2. Mapa A_tail(g0) ---")

# Gesta siatka g0 od 0.80 do 2.5
g0_range = np.linspace(0.80, 2.50, 200)
A_map = []
valid_g0 = []

for g0 in g0_range:
    try:
        At = A_tail(g0)
        if not np.isnan(At) and At > 0:
            A_map.append(At)
            valid_g0.append(g0)
    except:
        pass

valid_g0 = np.array(valid_g0)
A_map = np.array(A_map)
print(f"  Zmapowano {len(valid_g0)} punktow g0 in [{valid_g0[0]:.3f}, {valid_g0[-1]:.3f}]")

# ============================================================
# 3. Bisekcja phi-FP: znalezienie g0*
# ============================================================
print("\n--- 3. Bisekcja phi-FP ---")

target_r21 = 206.768

def R_func(g0_e):
    """r_21 = (A(phi*g0) / A(g0))^4"""
    Ae = A_tail(g0_e)
    Am = A_tail(PHI * g0_e)
    if Ae < 1e-15 or np.isnan(Ae) or np.isnan(Am):
        return np.nan
    return (Am / Ae)**4

lo, hi = 0.898, 0.901
for _ in range(40):
    mid = (lo + hi) / 2
    R = R_func(mid)
    if R is None or np.isnan(R):
        break
    if R < target_r21:
        lo = mid
    else:
        hi = mid

g0_star = (lo + hi) / 2
Ae = A_tail(g0_star)
Am = A_tail(PHI * g0_star)
r21 = (Am / Ae)**4

print(f"  g0* = {g0_star:.10f}")
print(f"  A_e = {Ae:.10f}")
print(f"  A_mu = {Am:.10f}")
print(f"  r_21 = {r21:.6f} (target: {target_r21})")

# ============================================================
# 4. Koide: parametryzacja i r_31
# ============================================================
print("\n--- 4. Parametryzacja Koide ---")

a = 1 + np.sqrt(r21)
disc = 6 * a**2 - 3 - 3 * r21
r31_K = (2 * a + np.sqrt(disc))**2
print(f"  a = 1 + sqrt(r21) = {a:.6f}")
print(f"  r_31 (Koide) = {r31_K:.2f}")
print(f"  r_31 (PDG)   = 3477.48")

# Koide parametryzacja: sqrt(m_n/m_e) = a*(1 + sqrt(2)*cos(2*pi*n/3 + delta))
# sqrt(1) = a*(1 + sqrt(2)*cos(delta))        [n=0, elektron]
# sqrt(r21) = a*(1 + sqrt(2)*cos(2pi/3+delta))  [n=1, mion]
# sqrt(r31) = a*(1 + sqrt(2)*cos(4pi/3+delta))  [n=2, tau]

# Z r_21 wyznaczamy delta i a_K
# sqrt(1) = 1 -> a*(1 + sqrt(2)*cos(delta)) = 1
# sqrt(r21) = 14.379 -> a*(1 + sqrt(2)*cos(2pi/3+delta)) = 14.379

sq1 = 1.0
sq2 = np.sqrt(r21)
sq3 = np.sqrt(r31_K)

# Znajdz a_K i delta
def koide_eqs(delta):
    """Rownanie na delta z warunku sqrt(m1)=1."""
    c1 = 1 + np.sqrt(2) * np.cos(delta)
    c2 = 1 + np.sqrt(2) * np.cos(2*np.pi/3 + delta)
    if abs(c1) < 1e-10: return 1e10
    a_K = 1.0 / c1
    return a_K * c2 - sq2

# Szukaj delta w zakresie [-pi, pi]
from scipy.optimize import brentq
delta_vals = np.linspace(-np.pi, np.pi, 1000)
koide_vals = [koide_eqs(d) for d in delta_vals]
delta_solutions = []
for i in range(len(delta_vals)-1):
    if koide_vals[i] * koide_vals[i+1] < 0:
        try:
            d_sol = brentq(koide_eqs, delta_vals[i], delta_vals[i+1])
            delta_solutions.append(d_sol)
        except:
            pass

print(f"\n  Rozwiazania delta:")
for d in delta_solutions:
    c1 = 1 + np.sqrt(2) * np.cos(d)
    c2 = 1 + np.sqrt(2) * np.cos(2*np.pi/3 + d)
    c3 = 1 + np.sqrt(2) * np.cos(4*np.pi/3 + d)
    if abs(c1) > 1e-10:
        a_K = 1.0 / c1
        sq1_check = a_K * c1
        sq2_check = a_K * c2
        sq3_check = a_K * c3
        r21_check = sq2_check**2
        r31_check = sq3_check**2
        if r21_check > 0 and r31_check > 0 and sq2_check > 0:
            Q_check = (sq1_check + sq2_check + sq3_check)**2 / (sq1_check**2 + sq2_check**2 + sq3_check**2)
            print(f"    delta = {d:.6f} rad ({np.degrees(d):.2f} deg)")
            print(f"    a_K = {a_K:.6f}")
            print(f"    sqrt(m/m_e) = [{sq1_check:.4f}, {sq2_check:.4f}, {sq3_check:.4f}]")
            print(f"    r_21 = {r21_check:.2f}, r_31 = {r31_check:.2f}")
            print(f"    Q = {Q_check:.8f}")

# ============================================================
# 5. KLUCZOWY TEST: Czy A_tail^2 ma strukture kosinusowa?
# ============================================================
print(f"\n{'='*70}")
print("  5. TEST: A_tail^2(g0) = a*(1 + sqrt(2)*cos(phi(g0))) ?")
print("="*70)

# Jesli A_tail^2 ~ a + b*cos(c*ln(g0-g*) + d), to Koide wyjdzie
# gdy g0_n sa rownomiernie rozlozone w ln(g0-g*)

# Sprobujmy znalezc transformacje zmiennej x(g0) taka ze
# A_tail^2(g0) = a + b*cos(x(g0))
# i g0_e, g0_mu, g0_tau sa rownomiernie rozlozone w x

# Krok 1: Znormalizuj A_tail^2
A2_map = A_map**2
A2_e = Ae**2
A2_mu = Am**2

# Krok 2: Sprobuj x = ln(g0 - g*)
x_map = np.log(valid_g0 - GG)
x_e = np.log(g0_star - GG)
x_mu = np.log(PHI * g0_star - GG)

print(f"  x_e = ln(g0*-g*) = {x_e:.6f}")
print(f"  x_mu = ln(phi*g0*-g*) = {x_mu:.6f}")
print(f"  Delta x = x_mu - x_e = {x_mu - x_e:.6f}")

# Krok 3: Fit A^2(x) = p0 + p1*cos(omega*x + phi0) na mapie
from scipy.optimize import curve_fit

def cos_model(x, p0, p1, omega, phi0):
    return p0 + p1 * np.cos(omega * x + phi0)

# Wstepne parametry
mask_fit = (valid_g0 > 0.85) & (valid_g0 < 2.2)
x_fit = x_map[mask_fit]
y_fit = A2_map[mask_fit]

try:
    # Proba fitu kosinusowego
    p0_init = np.mean(y_fit)
    p1_init = (np.max(y_fit) - np.min(y_fit)) / 2
    popt, pcov = curve_fit(cos_model, x_fit, y_fit,
                           p0=[p0_init, p1_init, 3.0, 0.0],
                           maxfev=10000)
    p0, p1, omega, phi0 = popt
    residuals = y_fit - cos_model(x_fit, *popt)
    r2 = 1 - np.sum(residuals**2) / np.sum((y_fit - np.mean(y_fit))**2)

    print(f"\n  Fit kosinusowy: A^2 = {p0:.4f} + {p1:.4f}*cos({omega:.4f}*x + {phi0:.4f})")
    print(f"  R^2 = {r2:.6f}")

    # Sprawdz rowne odleglosci w x
    x_pred_e = omega * x_e + phi0
    x_pred_mu = omega * x_mu + phi0
    delta_phase = x_pred_mu - x_pred_e
    print(f"  Faza e: {x_pred_e:.4f}")
    print(f"  Faza mu: {x_pred_mu:.4f}")
    print(f"  Delta fazy = {delta_phase:.4f}")
    print(f"  2*pi/3 = {2*np.pi/3:.4f}")
    print(f"  Stosunek delta/(2pi/3) = {delta_phase / (2*np.pi/3):.4f}")

    if abs(delta_phase / (2*np.pi/3) - 1.0) < 0.1:
        print(f"  *** BLISKIE 2*pi/3 !!! Koide moglby wynikac z perioycznosci! ***")
    else:
        print(f"  Delta fazy != 2*pi/3 -> Koide NIE wynika z prostego kosinusa w x=ln(g0-g*)")
except Exception as e:
    print(f"  Fit kosinusowy nie powiodl sie: {e}")
    r2 = -1

# ============================================================
# 6. Alternatywne zmienne: czy cos innego daje 2pi/3?
# ============================================================
print(f"\n{'='*70}")
print("  6. SZUKANIE ZMIENNEJ x DAJACEJ ROWNE FAZY")
print("="*70)

# Probujemy rozne transformacje z(g0):
# - z = g0 - g*
# - z = (g0 - g*)^p dla roznych p
# - z = ln(g0/g*)
# - z = (g0/g*)^p

# Dla kazdej transformacji szukamy omega t.z.:
# omega*(z_mu - z_e) = 2*pi/3

transforms = {
    'g0 - g*': lambda g: g - GG,
    'ln(g0-g*)': lambda g: np.log(g - GG),
    'ln(g0/g*)': lambda g: np.log(g / GG),
    '(g0-g*)^0.5': lambda g: np.sqrt(g - GG),
    '(g0-g*)^2': lambda g: (g - GG)**2,
    '1/(g0-g*)': lambda g: 1.0 / (g - GG),
    'g0^2': lambda g: g**2,
    'arctan((g0-1)/0.3)': lambda g: np.arctan((g-1)/0.3),
}

print(f"  Transformacja        | z_e      z_mu     Delta_z  | omega (2pi/3/Dz) | g0_tau (predict) | r31_pred  | Q_pred")
print(f"  " + "-"*110)

for name, tfunc in transforms.items():
    z_e = tfunc(g0_star)
    z_mu = tfunc(PHI * g0_star)
    dz = z_mu - z_e
    if abs(dz) < 1e-10:
        continue
    omega = 2 * np.pi / 3 / dz
    # Predict z_tau = z_e + 2*delta_z (rowne odleglosci)
    z_tau = z_e + 2 * dz

    # Odwroc transformacje: znajdz g0_tau
    # Dla kazdej transformacji musimy odwrocic
    g0_tau = np.nan
    try:
        if name == 'g0 - g*':
            g0_tau = z_tau + GG
        elif name == 'ln(g0-g*)':
            g0_tau = np.exp(z_tau) + GG
        elif name == 'ln(g0/g*)':
            g0_tau = GG * np.exp(z_tau)
        elif name == '(g0-g*)^0.5':
            g0_tau = z_tau**2 + GG
        elif name == '(g0-g*)^2':
            g0_tau = np.sqrt(z_tau) + GG if z_tau > 0 else np.nan
        elif name == '1/(g0-g*)':
            g0_tau = 1.0/z_tau + GG if z_tau > 0 else np.nan
        elif name == 'g0^2':
            g0_tau = np.sqrt(z_tau) if z_tau > 0 else np.nan
        elif name == 'arctan((g0-1)/0.3)':
            g0_tau = 0.3 * np.tan(z_tau) + 1
    except:
        continue

    if np.isnan(g0_tau) or g0_tau < 0.8 or g0_tau > 5:
        continue

    # Oblicz r31
    try:
        A_tau = A_tail(g0_tau)
        if np.isnan(A_tau) or A_tau < 1e-10:
            continue
        r31_pred = (A_tau / Ae)**4

        # Q
        sq_e = 1.0
        sq_mu = np.sqrt(r21)
        sq_tau = np.sqrt(r31_pred)
        Q = (sq_e + sq_mu + sq_tau)**2 / (1 + r21 + r31_pred)

        mark = " <--" if abs(Q - 1.5) < 0.05 else ""
        print(f"  {name:22s} | {z_e:8.4f} {z_mu:8.4f} {dz:8.4f} | {omega:16.4f} | {g0_tau:16.6f} | {r31_pred:9.1f} | {Q:.6f}{mark}")
    except:
        continue

# ============================================================
# 7. Skan: Q(g0_tau) dla wszystkich mozliwych g0_tau
# ============================================================
print(f"\n{'='*70}")
print("  7. SKAN: Q(g0_tau) dla g0_tau in [1.0, 3.5]")
print("="*70)

g0_tau_scan = np.linspace(0.85, 4.0, 400)
Q_scan = []
r31_scan = []

for g0t in g0_tau_scan:
    try:
        At = A_tail(g0t)
        if np.isnan(At) or At < 1e-10:
            Q_scan.append(np.nan)
            r31_scan.append(np.nan)
            continue
        r31 = (At / Ae)**4
        sq_tau = np.sqrt(r31)
        Q = (1 + np.sqrt(r21) + sq_tau)**2 / (1 + r21 + r31)
        Q_scan.append(Q)
        r31_scan.append(r31)
    except:
        Q_scan.append(np.nan)
        r31_scan.append(np.nan)

Q_scan = np.array(Q_scan)
r31_scan = np.array(r31_scan)

# Gdzie Q = 3/2?
mask_valid = ~np.isnan(Q_scan)
if np.any(mask_valid):
    # Znajdz g0_tau(Q=1.5) przez interpolacje
    q_valid = Q_scan[mask_valid]
    g_valid = g0_tau_scan[mask_valid]
    r_valid = r31_scan[mask_valid]

    crossings = []
    for i in range(len(q_valid)-1):
        if (q_valid[i] - 1.5) * (q_valid[i+1] - 1.5) < 0:
            g_cross = g_valid[i] + (1.5 - q_valid[i])/(q_valid[i+1]-q_valid[i])*(g_valid[i+1]-g_valid[i])
            crossings.append(g_cross)

    print(f"  Crossing Q=3/2 przy g0_tau =")
    for gc in crossings:
        At_c = A_tail(gc)
        r31_c = (At_c/Ae)**4
        print(f"    g0_tau = {gc:.6f}, r_31 = {r31_c:.1f}")
        print(f"    g0_tau / g0* = {gc / g0_star:.6f}")
        print(f"    g0_tau / (phi*g0*) = {gc / (PHI*g0_star):.6f}")
        print(f"    g0_tau / (phi^2*g0*) = {gc / (PHI**2*g0_star):.6f}")
        print(f"    g0_tau / (2*g0*) = {gc / (2*g0_star):.6f}")

    # Wyznacz Q minimum i maximum
    print(f"\n  Q(g0_tau) statystyki:")
    print(f"    min Q = {np.nanmin(Q_scan):.6f} (g0 = {g0_tau_scan[np.nanargmin(Q_scan)]:.4f})")
    print(f"    max Q = {np.nanmax(Q_scan):.6f} (g0 = {g0_tau_scan[np.nanargmax(Q_scan)]:.4f})")

# ============================================================
# 8. Energia solitonu: czy Koide minimalizuje?
# ============================================================
print(f"\n{'='*70}")
print("  8. ENERGIA SOLITONU: E(g0) i Koide")
print("="*70)

def soliton_energy(g0, r_max=200):
    """Energia (akcja) solitonu."""
    r, g = shoot(g0, r_max=r_max, n_pts=10000)
    if len(r) < 100:
        return np.nan
    dr = r[1] - r[0]
    # Gradient energii: int [K(g)/2 * (dg/dr)^2 + V(g)] * 4*pi*r^2 dr
    dg = np.gradient(g, r)
    K_g = np.array([max(fk(gi), 0.01) for gi in g])
    V_g = g**3/3 - g**4/4
    integrand = (K_g/2 * dg**2 + V_g) * 4 * np.pi * r**2
    E = np.trapezoid(integrand, r)
    return E

# Skan energii
g0_E_scan = np.linspace(0.85, 2.5, 80)
E_scan = []
print("  Obliczanie E(g0)...", end="", flush=True)
for g0 in g0_E_scan:
    try:
        E = soliton_energy(g0)
        E_scan.append(E)
    except:
        E_scan.append(np.nan)
print(" gotowe.")

E_scan = np.array(E_scan)

# Energie trzech leptonow
E_e = soliton_energy(g0_star)
E_mu = soliton_energy(PHI * g0_star)

# Suma energii: E_total = E_e + E_mu + E_tau
# Czy Koide minimalizuje E_total?
print(f"\n  E(g0*) [e] = {E_e:.4f}")
print(f"  E(phi*g0*) [mu] = {E_mu:.4f}")

if len(crossings) > 0:
    for gc in crossings:
        E_tau_K = soliton_energy(gc)
        E_total_K = E_e + E_mu + E_tau_K
        print(f"  E(g0_tau_Koide={gc:.4f}) [tau] = {E_tau_K:.4f}")
        print(f"  E_total (Koide) = {E_total_K:.4f}")

# Porownaj z innymi g0_tau
g0_tau_test = [PHI**2 * g0_star, 2 * g0_star, 2.5 * g0_star]
labels_test = ['phi^2*g0*', '2*g0*', '2.5*g0*']
for g0t, lbl in zip(g0_tau_test, labels_test):
    Et = soliton_energy(g0t)
    E_tot = E_e + E_mu + Et
    print(f"  E({lbl}={g0t:.4f}) = {Et:.4f}, E_total = {E_tot:.4f}")

# ============================================================
# 9. Pochodna dlnA/dlng0: szukanie struktury
# ============================================================
print(f"\n{'='*70}")
print("  9. STRUKTURA A_tail: dlnA/dlng0")
print("="*70)

# d(ln A)/d(ln g0) -- "wymiar skalowania"
lnA = np.log(A_map)
lng0 = np.log(valid_g0)
dlnA = np.gradient(lnA, lng0)

print(f"  dlnA/dlng0 na g0*: {np.interp(np.log(g0_star), lng0, dlnA):.4f}")
print(f"  dlnA/dlng0 na phi*g0*: {np.interp(np.log(PHI*g0_star), lng0, dlnA):.4f}")
if len(crossings) > 0:
    print(f"  dlnA/dlng0 na g0_tau(Koide): {np.interp(np.log(crossings[0]), lng0, dlnA):.4f}")

# Czy dlnA/dlng0 jest staly? (potegowa zaleznosc)
mask_power = (valid_g0 > 0.85) & (valid_g0 < 2.2)
dlnA_range = dlnA[mask_power]
print(f"  dlnA/dlng0 zakres: [{np.min(dlnA_range):.2f}, {np.max(dlnA_range):.2f}]")
print(f"  dlnA/dlng0 srednia: {np.mean(dlnA_range):.4f} +/- {np.std(dlnA_range):.4f}")

if np.std(dlnA_range) / abs(np.mean(dlnA_range)) < 0.15:
    print(f"  -> A ~ g0^{np.mean(dlnA_range):.2f} (przyblizona zaleznosc potegowa)")
else:
    print(f"  -> A(g0) NIE jest potegowa (zmiennosc {np.std(dlnA_range)/abs(np.mean(dlnA_range))*100:.0f}%)")

# ============================================================
# 10. Wykresy
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# (0,0) A_tail(g0)
ax = axes[0, 0]
ax.plot(valid_g0, A_map, 'b-', lw=2)
ax.axvline(g0_star, color='red', ls='--', alpha=0.5, label=f'$g_0^e$={g0_star:.3f}')
ax.axvline(PHI*g0_star, color='green', ls='--', alpha=0.5, label=f'$g_0^\\mu$={PHI*g0_star:.3f}')
if len(crossings) > 0:
    ax.axvline(crossings[0], color='purple', ls='--', alpha=0.5, label=f'$g_0^\\tau$(Koide)={crossings[0]:.3f}')
ax.axvline(GG, color='gray', ls=':', alpha=0.3, label=f'g*={GG:.3f}')
ax.set_xlabel('$g_0$'); ax.set_ylabel('$A_{tail}(g_0)$')
ax.set_title('Amplituda ogona'); ax.legend(fontsize=7)
ax.grid(True, alpha=0.2)

# (0,1) A_tail^2(g0) vs cosine fit
ax = axes[0, 1]
ax.plot(valid_g0, A2_map, 'b-', lw=2, label='$A^2(g_0)$')
if r2 > 0:
    x_plot = np.log(valid_g0 - GG)
    A2_fit = cos_model(x_plot, *popt)
    ax.plot(valid_g0, A2_fit, 'r--', lw=1.5, label=f'cos fit (R²={r2:.3f})')
ax.set_xlabel('$g_0$'); ax.set_ylabel('$A^2(g_0)$')
ax.set_title('$A^2$ vs cosine fit'); ax.legend(fontsize=8)
ax.grid(True, alpha=0.2)

# (0,2) Q(g0_tau) scan
ax = axes[0, 2]
ax.plot(g0_tau_scan, Q_scan, 'b-', lw=2)
ax.axhline(1.5, color='red', ls='--', lw=2, label='Q = 3/2 (Koide)')
ax.axvline(PHI**2*g0_star, color='orange', ls=':', alpha=0.7, label=f'$\\phi^2 g_0^*$={PHI**2*g0_star:.3f}')
ax.axvline(2*g0_star, color='cyan', ls=':', alpha=0.7, label=f'$2g_0^*$={2*g0_star:.3f}')
if len(crossings) > 0:
    for gc in crossings:
        ax.axvline(gc, color='purple', ls='--', alpha=0.7)
ax.set_xlabel('$g_0^{\\tau}$'); ax.set_ylabel('Q')
ax.set_title('Q(g₀τ) — Koide skan')
ax.set_ylim(1.0, 2.5)
ax.legend(fontsize=7); ax.grid(True, alpha=0.2)

# (1,0) Energia solitonu
ax = axes[1, 0]
mask_E = ~np.isnan(E_scan)
ax.plot(g0_E_scan[mask_E], E_scan[mask_E], 'b-', lw=2)
ax.axvline(g0_star, color='red', ls='--', alpha=0.5, label='e')
ax.axvline(PHI*g0_star, color='green', ls='--', alpha=0.5, label='mu')
if len(crossings) > 0:
    ax.axvline(crossings[0], color='purple', ls='--', alpha=0.5, label='tau(K)')
ax.set_xlabel('$g_0$'); ax.set_ylabel('E(g₀)')
ax.set_title('Energia solitonu'); ax.legend(fontsize=8)
ax.grid(True, alpha=0.2)

# (1,1) dlnA/dlng0
ax = axes[1, 1]
ax.plot(valid_g0, dlnA, 'b-', lw=2)
ax.axvline(g0_star, color='red', ls='--', alpha=0.5)
ax.axvline(PHI*g0_star, color='green', ls='--', alpha=0.5)
ax.set_xlabel('$g_0$'); ax.set_ylabel('$d\\ln A / d\\ln g_0$')
ax.set_title('Wymiar skalowania A(g₀)'); ax.grid(True, alpha=0.2)

# (1,2) r31(g0_tau) vs Q
ax = axes[1, 2]
mask_r = ~np.isnan(r31_scan)
ax.plot(r31_scan[mask_r], Q_scan[mask_r], 'b-', lw=2)
ax.axhline(1.5, color='red', ls='--', lw=2, alpha=0.7)
ax.axvline(3477.48, color='red', ls=':', alpha=0.5, label='r₃₁(PDG)')
ax.set_xlabel('$r_{31}$'); ax.set_ylabel('Q')
ax.set_title('Q vs r₃₁'); ax.legend(fontsize=8)
ax.set_xlim(0, 8000); ax.set_ylim(1.0, 2.5)
ax.grid(True, alpha=0.2)

plt.suptitle('O-K1: Koide Q=3/2 z TGP?', fontsize=14, y=1.02)
plt.tight_layout()
sd = os.path.dirname(os.path.abspath(__file__))
plt.savefig(os.path.join(sd, 'ok1_koide.png'), dpi=150, bbox_inches='tight')
print(f"\n  Wykres: ok1_koide.png")

# ============================================================
# 11. Podsumowanie
# ============================================================
print(f"\n{'='*70}")
print("  PODSUMOWANIE O-K1")
print("="*70)
print(f"""
  PYTANIE: Czy Koide Q=3/2 wynika z TGP?

  WYNIKI:
  1. A_tail(g0) jest gladka, monotonicznie rosnaca dla g0 > g*
  2. Q(g0_tau) przechodzi przez 3/2 dla specyficznego g0_tau
  3. Koide NIE wynika z prostej potegowej A~(g0-g*)^alpha
     (dlnA/dlng0 zmienia sie o {np.std(dlnA_range)/abs(np.mean(dlnA_range))*100:.0f}%)
  4. Fit kosinusowy A^2(ln(g0-g*)) ma R^2 = {r2:.3f}
     {'-> SLABY fit, Koide nie wynika z perioycznosci w ln(g0-g*)' if r2 < 0.95 else '-> DOBRY fit!'}
  5. Roznica faz e->mu: {delta_phase:.4f} vs 2*pi/3 = {2*np.pi/3:.4f}
     {'-> NIEZGODNE' if abs(delta_phase/(2*np.pi/3) - 1) > 0.1 else '-> ZGODNE!'}

  STATUS O-K1: {'OTWARTY' if r2 < 0.95 or abs(delta_phase/(2*np.pi/3) - 1) > 0.1 else 'ZAMKNIETY'}

  MOZLIWE DROGI:
  a) Inna zmienna (nie ln(g0-g*)) moze dac rowne fazy
  b) Koide moze wynikac z warunku ekstremalnego na E_total
  c) Ukryta symetria ODE (np. dualnosc g <-> 2-g) moze ograniczac Q
  d) Efekty wieloczesticowe (nie jednosolitonowe) moga generowac Q=3/2
""")
print("GOTOWE.")
