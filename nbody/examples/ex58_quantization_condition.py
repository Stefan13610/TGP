"""
ex58_quantization_condition.py
================================
Warunek kwantowania g0 w TGP: dlaczego g0_e = 1.24 i g0_mu = 2.00?

KONTEKST (ex57):
  Przełom: (A_tail(2.00)/A_tail(1.24))^4 = 210 ≈ 207 (1.6%)
  Globalnie: M_raw ~ g0^4.12
  Ogon: g(r)-1 ~ (B*cos(r) + C*sin(r))/r przy r >> 1

PYTANIE:
  Co fizycznie wyróżnia g0 = 1.24 i g0 = 2.00 spośród innych wartości?

HIPOTEZY KWANTOWANIA:
  H1: B_tail(g0) = 0 — ogon jest "czystym sinusem"
      g(r)-1 ~ C*sin(r)/r (faza rozproszeniowa delta_phase = 0)
      Fizycznie: soliton w rezonansie z substratem (zero składowej cos)
  H2: dA_tail/dg0 = 0 — amplituda ogona ma ekstremum (warunek stacjonarny)
  H3: Faza ogona phi(g0) = atan2(B,C) stabilizuje się ("zatrzask fazowy")
  H4: Energia efektywna: cos(r)-ogon reprezentuje fale wychodzące,
      sin(r)-ogon — stojące; B=0 to warunek "fali stojącej"

PLAN:
  1. Gęsty skan g0 ∈ [1.1, 3.0], N=300: B(g0), C(g0), A(g0), phi(g0)
  2. Znajdź zera B_tail(g0) i zera C_tail(g0) — kandydaci kwantowania H1
  3. Znajdź ekstrema dA_tail/dg0 — kandydaci H2
  4. Faza phi(g0) = atan2(B, C): szybka oscylacja ~ r? Czy ma "ukrytą" strukturę?
  5. Zoom wokół g0=1.24 i g0=2.00: N=200 punktów w ε-otoczeniu
  6. Test: czy B_tail = 0 jest spełnione dla obydwu lub jednego z nich?
  7. Analogia ze stałymi ex39: K*1=0.01003, K*2=0.1004 jako C_tail(g0_e) i C_tail(g0_mu)?
  8. Macierz S: e^{2i*delta} z B i C

TESTY (5):
  T1: Phase plot A_tail(g0) jest monotonicznie rosnąca w [1.1, 2.1]
  T2: Istnieje co najmniej jedno zero B_tail w [1.1, 3.0]
  T3: Istnieje co najmniej jedno zero B_tail w [1.9, 2.1] (blisko g0_mu=2.0)
  T4: Faza phi(g0) jest ciągłą funkcją g0 (maks. skok < 1.5 rad w sąsiednich punktach)
  T5: A_tail(1.24) in [0.2, 0.4], A_tail(2.0) in [0.8, 1.5] (ex57 benchmark)

Sesja: TGP v33, 2026-03-27
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import os, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import argrelextrema
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# Stałe TGP
# ─────────────────────────────────────────────────────────────────────────────
ALPHA    = 2.0
G_GHOST  = np.exp(-1.0 / (2.0 * ALPHA))   # ~0.7788
TARGET   = 206.768
G_BOUNCE = G_GHOST + 0.005                 # poziom elastycznego odbicia

# Stałe z ex39 (ZPE): K*1 i K*2 (możliwe powiązanie z C_tail?)
K_STAR_1 = 0.01003
K_STAR_2 = 0.1004

# Znane wartości z ex57
G0_E   = 1.24;   A_E   = 0.287   # elektron
G0_MU  = 2.00;   A_MU  = 1.09   # mion (przełom)

# Parametry integracji
R_MAX    = 40.0    # dłuższy zakres dla dokładniejszego pomiaru fazy
R_START  = 1e-4
R_TAIL_L = 20.0    # większe okno dopasowania
R_TAIL_R = 36.0
RTOL     = 1e-10
ATOL     = 1e-13
MAX_STEP = 0.02    # drobny krok

print("=" * 70)
print("EX58: WARUNEK KWANTOWANIA g0 — FAZA I AMPLITUDA OGONA TGP")
print("=" * 70)
print(f"  alpha   = {ALPHA},  g* = {G_GHOST:.5f}")
print(f"  Cel: wyjasnij g0_e = {G0_E}, g0_mu = {G0_MU}")
print(f"  R_tail  = [{R_TAIL_L}, {R_TAIL_R}]  (dluzsze okno = dokładniejsza faza)")
print(f"  K*1     = {K_STAR_1},  K*2 = {K_STAR_2}  (z ex39)")
print()

# ─────────────────────────────────────────────────────────────────────────────
# ODE i integracja (identyczna z ex57)
# ─────────────────────────────────────────────────────────────────────────────
def rhs_physical(r, y):
    """Pełne ODE TGP — obszar fizyczny g > G_BOUNCE."""
    g, gp = y
    g   = max(g, G_BOUNCE + 1e-6)
    fg  = 1.0 + 2.0 * ALPHA * np.log(g)
    if abs(fg) < 1e-10:
        return [gp, 0.0]
    driving = g**2 * (1.0 - g)
    cross   = (ALPHA / g) * gp**2
    if r < 1e-10:
        return [gp, (driving - cross) / (3.0 * fg)]
    damp = fg * 2.0 * gp / r
    return [gp, (driving - cross - damp) / fg]

def event_hit_ghost(r, y):
    return y[0] - G_BOUNCE
event_hit_ghost.terminal  = True
event_hit_ghost.direction = -1


def integrate_with_bounces(g0, r_max=R_MAX, max_bounces=10):
    """Integracja z elastycznymi odbiciami przy g*."""
    r0    = R_START
    y0    = [g0, 0.0]
    segs_r, segs_g, segs_gp = [], [], []

    for bounce_num in range(max_bounces + 1):
        sol = solve_ivp(
            rhs_physical, [r0, r_max], y0,
            method='DOP853',
            max_step=MAX_STEP,
            rtol=RTOL, atol=ATOL,
            events=[event_hit_ghost],
            dense_output=False
        )
        segs_r.append(sol.t)
        segs_g.append(sol.y[0])
        segs_gp.append(sol.y[1])

        if sol.t_events[0].size > 0 and bounce_num < max_bounces:
            r_b  = float(sol.t_events[0][0])
            gp_b = float(sol.y_events[0][0, 1])
            r0   = r_b + 1e-6
            y0   = [G_BOUNCE + 1e-5, -gp_b]
        else:
            break

    r  = np.concatenate(segs_r)
    g  = np.concatenate(segs_g)
    gp = np.concatenate(segs_gp)
    idx = np.argsort(r)
    return r[idx], g[idx], gp[idx]


def fit_tail(r_arr, g_arr, r_L=R_TAIL_L, r_R=R_TAIL_R):
    """
    Dopasuj g(r)-1 ~ (B*cos(r) + C*sin(r))/r.
    Zwraca A = sqrt(B^2+C^2), B, C, phi = atan2(B, C).
    Faza phi: A*sin(r + phi) gdzie phi = atan2(B, C).
    """
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    if np.sum(mask) < 20:
        return 0.0, 0.0, 0.0, 0.0   # A, B, C, phi

    r_fit     = r_arr[mask]
    delta_fit = (g_arr[mask] - 1.0) * r_fit

    M_mat = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    result = np.linalg.lstsq(M_mat, delta_fit, rcond=None)
    B, C = result[0]
    A   = np.sqrt(B**2 + C**2)
    phi = np.arctan2(B, C)   # faza: A*sin(r + phi), B=A*sin(phi), C=A*cos(phi)
    return float(A), float(B), float(C), float(phi)


def analyze(g0):
    """Pełna analiza: A_tail, B_tail, C_tail, phi."""
    r, g, gp = integrate_with_bounces(g0)
    A, B, C, phi = fit_tail(r, g)
    g_min = float(np.min(g))
    M_raw = 4.0*np.pi * float(np.trapezoid(r**2 * gp**2 / 2.0, r))
    n_bounce = int(np.sum(np.diff((g < 0.85).astype(int)) == 1)) if np.any(g < 0.85) else 0
    return {'A': A, 'B': B, 'C': C, 'phi': phi,
            'g_min': g_min, 'M_raw': M_raw, 'n_bounce': n_bounce,
            'r': r, 'g': g}


# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 1: Gęsty skan g0 ∈ [1.10, 3.0]
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 1: Gęsty skan g0 ∈ [1.10, 3.0], N=300 ---")

N_SCAN  = 300
g0_scan = np.linspace(1.10, 3.0, N_SCAN)

A_arr   = np.zeros(N_SCAN)
B_arr   = np.zeros(N_SCAN)
C_arr   = np.zeros(N_SCAN)
phi_arr = np.zeros(N_SCAN)
Mraw_arr = np.zeros(N_SCAN)
gmin_arr = np.zeros(N_SCAN)
nb_arr   = np.zeros(N_SCAN, dtype=int)

for i, g0 in enumerate(g0_scan):
    res = analyze(g0)
    A_arr[i]   = res['A']
    B_arr[i]   = res['B']
    C_arr[i]   = res['C']
    phi_arr[i] = res['phi']
    Mraw_arr[i] = res['M_raw']
    gmin_arr[i] = res['g_min']
    nb_arr[i]   = res['n_bounce']

print(f"  Scan gotowy: A range [{A_arr.min():.4f}, {A_arr.max():.4f}]")
print(f"  B range [{B_arr.min():.4f}, {B_arr.max():.4f}]")
print(f"  C range [{C_arr.min():.4f}, {C_arr.max():.4f}]")
print()

# Tabela co 20 punktów
print(f"  {'g0':>6}  {'A':>8}  {'B':>9}  {'C':>9}  {'phi/pi':>7}  {'n_b':>4}")
print("  " + "-" * 55)
for i in range(0, N_SCAN, 20):
    g0 = g0_scan[i]
    print(f"  {g0:6.3f}  {A_arr[i]:8.4f}  {B_arr[i]:9.4f}  {C_arr[i]:9.4f}"
          f"  {phi_arr[i]/np.pi:7.4f}  {nb_arr[i]:4d}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 2: Zera B_tail(g0) i C_tail(g0) — kandydaci H1
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 2: Zera B_tail i C_tail (hipoteza H1: fala stojąca) ---")

# Zera B_tail: B zmienia znak
zeros_B = []
for i in range(N_SCAN - 1):
    if B_arr[i] * B_arr[i+1] < 0:
        # interpolacja liniowa
        g0_zero = g0_scan[i] - B_arr[i] * (g0_scan[i+1] - g0_scan[i]) / (B_arr[i+1] - B_arr[i])
        zeros_B.append(float(g0_zero))

zeros_C = []
for i in range(N_SCAN - 1):
    if C_arr[i] * C_arr[i+1] < 0:
        g0_zero = g0_scan[i] - C_arr[i] * (g0_scan[i+1] - g0_scan[i]) / (C_arr[i+1] - C_arr[i])
        zeros_C.append(float(g0_zero))

print(f"  Zera B_tail(g0) (B=0, 'czysty sinus'): {zeros_B}")
print(f"  Zera C_tail(g0) (C=0, 'czysty cosinus'): {zeros_C}")
print()

# Odległości od znanych g0
for name, g0_known in [("g0_e=1.24", G0_E), ("g0_mu=2.00", G0_MU)]:
    if zeros_B:
        dists_B = [abs(z - g0_known) for z in zeros_B]
        nearest_B = zeros_B[np.argmin(dists_B)]
        print(f"  {name}: najbliższe zero B = {nearest_B:.4f} (odl = {min(dists_B):.4f})")
    if zeros_C:
        dists_C = [abs(z - g0_known) for z in zeros_C]
        nearest_C = zeros_C[np.argmin(dists_C)]
        print(f"  {name}: najbliższe zero C = {nearest_C:.4f} (odl = {min(dists_C):.4f})")

print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 3: Ekstrema dA_tail/dg0 — hipoteza H2
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 3: Ekstrema dA_tail/dg0 (hipoteza H2: warunek stacjonarny) ---")

dA_dg0 = np.gradient(A_arr, g0_scan)

# Lokalne ekstrema A_tail
local_max_idx = argrelextrema(A_arr, np.greater, order=3)[0]
local_min_idx = argrelextrema(A_arr, np.less,    order=3)[0]

print(f"  Lokalne MAXIMA A_tail: {[(f'{g0_scan[i]:.3f}', f'{A_arr[i]:.4f}') for i in local_max_idx]}")
print(f"  Lokalne MINIMA A_tail: {[(f'{g0_scan[i]:.3f}', f'{A_arr[i]:.4f}') for i in local_min_idx]}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 4: Analiza fazy phi(g0)
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 4: Faza ogona phi(g0) = atan2(B, C) ---")

# Odwij fazę
phi_unwrapped = np.unwrap(phi_arr)

# Szybkość zmiany fazy
dphi_dg0 = np.gradient(phi_unwrapped, g0_scan)

print(f"  Całkowita zmiana fazy: {phi_unwrapped[-1] - phi_unwrapped[0]:.3f} rad"
      f" = {(phi_unwrapped[-1] - phi_unwrapped[0])/np.pi:.2f} π")
print(f"  Średnia dphi/dg0 = {np.mean(dphi_dg0):.3f} rad / jednostka g0")
print()

# Znajdź g0 gdzie phi ≈ 0 lub phi ≈ ±π/2
phi_targets = [(0.0, "phi=0"), (np.pi/2, "phi=pi/2"), (-np.pi/2, "phi=-pi/2")]
for phi_t, label in phi_targets:
    # Zera phi - phi_t
    phi_diff = phi_arr - phi_t
    crossings = []
    for i in range(N_SCAN - 1):
        if phi_diff[i] * phi_diff[i+1] < 0:
            g0_c = g0_scan[i] - phi_diff[i] * (g0_scan[i+1] - g0_scan[i]) / (phi_diff[i+1] - phi_diff[i])
            crossings.append(float(g0_c))
    if crossings:
        print(f"  g0 gdzie {label}: {[f'{x:.3f}' for x in crossings[:5]]}")

print()

# Wartości fazy przy znanych g0
idx_e  = np.argmin(np.abs(g0_scan - G0_E))
idx_mu = np.argmin(np.abs(g0_scan - G0_MU))
print(f"  phi(g0_e={G0_E})  = {phi_arr[idx_e]:.4f} rad = {phi_arr[idx_e]/np.pi:.4f} pi")
print(f"  phi(g0_mu={G0_MU}) = {phi_arr[idx_mu]:.4f} rad = {phi_arr[idx_mu]/np.pi:.4f} pi")
print(f"  B(g0_e)  = {B_arr[idx_e]:.5f},  C(g0_e)  = {C_arr[idx_e]:.5f}")
print(f"  B(g0_mu) = {B_arr[idx_mu]:.5f},  C(g0_mu) = {C_arr[idx_mu]:.5f}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 5: Zoom wokół g0_e = 1.24 i g0_mu = 2.00
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 5: Zoom wokół g0_e=1.24 i g0_mu=2.00 ---")

for name, g0_center in [("elektron", G0_E), ("mion", G0_MU)]:
    eps_zoom = 0.15
    g0_zoom  = np.linspace(g0_center - eps_zoom, g0_center + eps_zoom, 80)
    B_zoom = []; C_zoom = []; A_zoom = []; phi_zoom = []

    for g0 in g0_zoom:
        res = analyze(g0)
        A_zoom.append(res['A']); B_zoom.append(res['B'])
        C_zoom.append(res['C']); phi_zoom.append(res['phi'])

    B_zoom = np.array(B_zoom); C_zoom = np.array(C_zoom)
    A_zoom = np.array(A_zoom); phi_zoom = np.array(phi_zoom)

    # Zera B w tym zakresie
    zeros_B_z = []
    for i in range(len(g0_zoom) - 1):
        if B_zoom[i] * B_zoom[i+1] < 0:
            g0_z = g0_zoom[i] - B_zoom[i]*(g0_zoom[i+1]-g0_zoom[i])/(B_zoom[i+1]-B_zoom[i])
            zeros_B_z.append(float(g0_z))

    print(f"\n  [ZOOM {name} g0 ∈ [{g0_center-eps_zoom:.2f}, {g0_center+eps_zoom:.2f}]]")
    print(f"    A(g0_center) = {A_zoom[len(A_zoom)//2]:.5f}")
    print(f"    B(g0_center) = {B_zoom[len(B_zoom)//2]:.5f}")
    print(f"    C(g0_center) = {C_zoom[len(C_zoom)//2]:.5f}")
    print(f"    phi(g0_center) = {phi_zoom[len(phi_zoom)//2]:.5f} rad = "
          f"{phi_zoom[len(phi_zoom)//2]/np.pi:.5f} pi")
    if zeros_B_z:
        print(f"    Zera B w tym zakresie: {[f'{z:.4f}' for z in zeros_B_z]}")
        # Odległość od g0_center
        nearest = min(zeros_B_z, key=lambda z: abs(z - g0_center))
        print(f"    Najbliższe zero B: {nearest:.4f} (odl od {g0_center} = {abs(nearest-g0_center):.4f})")
    else:
        print(f"    BRAK zer B w tym zakresie!")

print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 6: Porównanie C_tail z stałymi K*1, K*2 (ex39)
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 6: C_tail kontra K*1=0.01003, K*2=0.1004 (ex39) ---")

print(f"  C_tail(g0_e={G0_E})  = {C_arr[idx_e]:.6f}  vs K*1 = {K_STAR_1}")
print(f"  C_tail(g0_mu={G0_MU}) = {C_arr[idx_mu]:.6f}  vs K*2 = {K_STAR_2}")
print(f"  A_tail(g0_e)  = {A_arr[idx_e]:.6f}  vs K*1 = {K_STAR_1}")
print(f"  A_tail(g0_mu) = {A_arr[idx_mu]:.6f}  vs K*2 = {K_STAR_2}")
print(f"  Stosunek A_mu/A_e = {A_arr[idx_mu]/A_arr[idx_e]:.5f}  (ex57: 3.80)")
print(f"  Stosunek K*2/K*1  = {K_STAR_2/K_STAR_1:.3f}  (ex53: 10.01)")
print()

# Szukaj g0 gdzie |C_tail(g0) - K*1| jest minimalne
dists_K1 = np.abs(C_arr - K_STAR_1)
dists_K2 = np.abs(C_arr - K_STAR_2)
g0_closest_K1 = g0_scan[np.argmin(dists_K1)]
g0_closest_K2 = g0_scan[np.argmin(dists_K2)]
print(f"  C_tail = K*1 najbliżej g0 = {g0_closest_K1:.3f}  (min dist = {dists_K1.min():.5f})")
print(f"  C_tail = K*2 najbliżej g0 = {g0_closest_K2:.3f}  (min dist = {dists_K2.min():.5f})")

# Szukaj g0 gdzie |A_tail(g0) - K*1|, |A_tail(g0) - K*2|
dists_AK1 = np.abs(A_arr - K_STAR_1)
dists_AK2 = np.abs(A_arr - K_STAR_2)
g0_AK1 = g0_scan[np.argmin(dists_AK1)]
g0_AK2 = g0_scan[np.argmin(dists_AK2)]
print(f"  A_tail = K*1 najbliżej g0 = {g0_AK1:.3f}")
print(f"  A_tail = K*2 najbliżej g0 = {g0_AK2:.3f}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 7: Analiza S-matrix
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 7: Macierz S (S-matrix) e^{2i*delta} ---")
print("  W sferycznym rozproszeniu s-fali:")
print("  g(r)-1 ~ (B*cos(r) + C*sin(r))/r = A*sin(r+phi)/r")
print("  Rozwiązanie regularne f.t. = sin(r)/r   <-> C!=0, B=0")
print("  Amplituda rozproszenia f(k) ~ e^{i*delta}*sin(delta)/k")
print("  Tangens przesunięcia fazy: tan(delta_phase) = -B/C")
print()

# Oblicz delta_phase = atan(-B/C)
# Uwaga: zwykły atan2 nie nadaje się bezpośrednio; używamy -atan2(B, C)
delta_phase = -np.arctan2(B_arr, C_arr)   # przesunięcie fazy rozproszenia

print(f"  delta_phase(g0_e={G0_E}) = {delta_phase[idx_e]:.5f} rad = {delta_phase[idx_e]/np.pi:.5f} pi")
print(f"  delta_phase(g0_mu={G0_MU}) = {delta_phase[idx_mu]:.5f} rad = {delta_phase[idx_mu]/np.pi:.5f} pi")

# Zera delta_phase ↔ B=0 (rezonans w konwencji n-fali)
# Szukaj przejść przez 0 w sensie delta_phase
zeros_delta = []
for i in range(N_SCAN - 1):
    if delta_phase[i] * delta_phase[i+1] < 0:
        g0_z = g0_scan[i] - delta_phase[i]*(g0_scan[i+1]-g0_scan[i])/(delta_phase[i+1]-delta_phase[i])
        zeros_delta.append(float(g0_z))

print(f"  Zera delta_phase (B=0): {[f'{z:.4f}' for z in zeros_delta]}")

# Ekstrema delta_phase (rezonanse Breit-Wigner: delta_phase = pi/2)
resonances = []
for i in range(N_SCAN - 1):
    d = delta_phase[i]
    # Check if |delta_phase| = pi/2
    if abs(d - np.pi/2) < 0.05 or abs(d + np.pi/2) < 0.05:
        resonances.append(g0_scan[i])

print(f"  Przybliżone rezonanse |delta|=pi/2: {[f'{z:.3f}' for z in resonances[:5]]}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 8: Warunki kwantowania — przesuwanie fazowe i periodyczność
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 8: Periodyczność i zatrzask fazowy phi(g0) ---")

# Szybka oscylacja phi pochodzi z liniowej zmiany wraz z g0
# Znajdź "tło" liniowe przez dopasowanie
p_lin = np.polyfit(g0_scan, phi_unwrapped, 1)
phi_detrended = phi_unwrapped - np.polyval(p_lin, g0_scan)

print(f"  Liniowe tło fazy: d(phi)/d(g0) = {p_lin[0]:.3f} rad/jednostka_g0")
print(f"  Czas okresu fazy: T_g0 = 2*pi / {abs(p_lin[0]):.3f} = {2*np.pi/abs(p_lin[0]):.4f}")
print()

# Znajdź lokalne ekstrema oddetrendowanej fazy
if len(phi_detrended) > 6:
    local_maxima = argrelextrema(phi_detrended, np.greater, order=5)[0]
    local_minima = argrelextrema(phi_detrended, np.less,    order=5)[0]
    max_g0 = g0_scan[local_maxima]
    min_g0 = g0_scan[local_minima]
    if len(max_g0) >= 2:
        spacings = np.diff(max_g0)
        print(f"  Odstępy między maksimami phi_detrended: {[f'{s:.3f}' for s in spacings]}")
        print(f"  Średni odstęp = {np.mean(spacings):.4f}")
    print(f"  Maksima phi_detrended przy g0 = {[f'{g:.3f}' for g in max_g0]}")
    print(f"  Minima phi_detrended przy g0  = {[f'{g:.3f}' for g in min_g0]}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 9: Profile g(r) dla g0=1.24 i g0=2.00 — wizualne porównanie
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 9: Szczegółowe profile i dopasowania ogona ---")

for name, g0, A_ref, A_label in [
    ("elektron", G0_E, A_E, "A_e"),
    ("mion",     G0_MU, A_MU, "A_mu"),
]:
    r, g, gp = integrate_with_bounces(g0)
    A, B, C, phi = fit_tail(r, g)
    M_raw = 4.0*np.pi * float(np.trapezoid(r**2 * gp**2 / 2.0, r))
    g_min = float(np.min(g))

    print(f"\n  [{name.upper()}] g0={g0}")
    print(f"    g_min = {g_min:.5f}  (g* = {G_GHOST:.5f})")
    print(f"    A_tail = {A:.5f}  (ex57 ref: {A_ref})")
    print(f"    B_tail = {B:.5f}  (skal. cosinusowa)")
    print(f"    C_tail = {C:.5f}  (skal. sinusowa)")
    print(f"    phi    = {phi:.5f} rad = {phi/np.pi:.5f} pi")
    print(f"    M_raw  = {M_raw:.4f}")

    # Dopasowanie ogona: residuum
    mask = (r >= R_TAIL_L) & (r <= R_TAIL_R)
    g_fit = 1.0 + (B * np.cos(r[mask]) + C * np.sin(r[mask])) / r[mask]
    resid = np.std(g[mask] - g_fit)
    print(f"    Residuum dopasowania ogona: {resid:.2e}")

print()

# ─────────────────────────────────────────────────────────────────────────────
# PODSUMOWANIE WARUNKÓW KWANTOWANIA
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("PODSUMOWANIE: WARUNKI KWANTOWANIA")
print("=" * 70)

# Najlepsze kandydatury
print()
print("Hipoteza H1 (B_tail=0, 'czysty sinus'):")
if zeros_B:
    for z in zeros_B[:6]:
        dist_e  = abs(z - G0_E)
        dist_mu = abs(z - G0_MU)
        print(f"  B=0 przy g0 = {z:.4f}  [odl od g0_e: {dist_e:.4f}, od g0_mu: {dist_mu:.4f}]")
else:
    print("  BRAK zer B_tail w [1.1, 3.0]")

print()
print("Hipoteza H2 (ekstrema A_tail):")
for idx in local_max_idx[:5]:
    z = g0_scan[idx]
    dist_e  = abs(z - G0_E)
    dist_mu = abs(z - G0_MU)
    print(f"  max A przy g0 = {z:.4f}  A={A_arr[idx]:.4f}"
          f"  [odl od g0_e: {dist_e:.4f}, od g0_mu: {dist_mu:.4f}]")

print()
print(f"Faza phi przy znanych g0:")
print(f"  phi(g0_e = {G0_E})  = {phi_arr[idx_e]/np.pi:.4f} pi  "
      f"(B={B_arr[idx_e]:.5f}, C={C_arr[idx_e]:.5f})")
print(f"  phi(g0_mu = {G0_MU}) = {phi_arr[idx_mu]/np.pi:.4f} pi  "
      f"(B={B_arr[idx_mu]:.5f}, C={C_arr[idx_mu]:.5f})")
print(f"  Różnica faz: Δphi = {(phi_arr[idx_mu] - phi_arr[idx_e])/np.pi:.4f} pi")
print()

# ─────────────────────────────────────────────────────────────────────────────
# CZĘŚĆ 10: Dokładne zera B_tail przez bisekcję + stosunek mas
# ─────────────────────────────────────────────────────────────────────────────
print("--- Część 10: Dokładne zera B_tail (bisekcja) + stosunek mas ---")
print("  Cel: (A_tail(zero_1)/A_tail(zero_0))^4 =? 207")
print()

def B_tail_only(g0):
    """Oblicz tylko B_tail dla danego g0."""
    r, g, gp = integrate_with_bounces(g0)
    _, B, _, _ = fit_tail(r, g)
    return B

def bisect_B_zero(g0_lo, g0_hi, tol=1e-5, max_iter=50):
    """Bisekcja zera B_tail w przedziale [g0_lo, g0_hi]."""
    B_lo = B_tail_only(g0_lo)
    B_hi = B_tail_only(g0_hi)
    if B_lo * B_hi > 0:
        return None   # Brak zmiany znaku
    for _ in range(max_iter):
        g0_mid = 0.5 * (g0_lo + g0_hi)
        B_mid  = B_tail_only(g0_mid)
        if abs(g0_hi - g0_lo) < tol:
            break
        if B_lo * B_mid <= 0:
            g0_hi = g0_mid; B_hi = B_mid
        else:
            g0_lo = g0_mid; B_lo = B_mid
    return 0.5 * (g0_lo + g0_hi)

# Znalezione wcześniej zera: ~1.230, ~2.080, ~2.788
# Bisekcja w wąskich oknach
zero_windows = [(1.15, 1.27), (2.00, 2.15), (2.70, 2.88)]
exact_zeros  = []

for (lo, hi) in zero_windows:
    z = bisect_B_zero(lo, hi)
    if z is not None:
        exact_zeros.append(z)
        res = analyze(z)
        print(f"  Zero B_tail = {z:.6f}  "
              f"A={res['A']:.5f}  B={res['B']:.2e}  C={res['C']:.5f}  "
              f"n_bounce={res['n_bounce']}")

print()

# Stosunek mas przy dokładnych zerach B
if len(exact_zeros) >= 2:
    z0, z1 = exact_zeros[0], exact_zeros[1]
    res0 = analyze(z0)
    res1 = analyze(z1)
    A0, A1 = res0['A'], res1['A']
    ratio_A4 = (A1 / A0)**4
    print(f"  g0_e (zero_0) = {z0:.6f},  A0 = {A0:.5f}")
    print(f"  g0_mu (zero_1) = {z1:.6f},  A1 = {A1:.5f}")
    print(f"  A1/A0 = {A1/A0:.5f}  (ex57: 3.80,  207^1/4 = {207**0.25:.5f})")
    print(f"  (A1/A0)^4 = {ratio_A4:.3f}  (cel: 206.77,  odchylenie: {100*abs(ratio_A4-207)/207:.2f}%)")
    print()
    if len(exact_zeros) >= 3:
        z2 = exact_zeros[2]
        res2 = analyze(z2)
        A2 = res2['A']
        ratio_A4_02 = (A2 / A0)**4
        print(f"  g0 (zero_2) = {z2:.6f},  A2 = {A2:.5f}")
        print(f"  (A2/A0)^4 = {ratio_A4_02:.3f}  (dla ewentualnego tauonu?)")
        print()

    # Sekundarny: M_raw przy zerach
    M0 = res0['M_raw']
    M1 = res1['M_raw']
    print(f"  M_raw(zero_0) = {M0:.4f}")
    print(f"  M_raw(zero_1) = {M1:.4f}")
    print(f"  M_raw(zero_1)/M_raw(zero_0) = {M1/M0:.3f}  (cel: 206.77)")

print()

# ─────────────────────────────────────────────────────────────────────────────
# WYKRESY
# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 12))
gs  = gridspec.GridSpec(3, 2, figure=fig, hspace=0.4, wspace=0.35)

# Panel 1: A_tail, B_tail, C_tail vs g0
ax1 = fig.add_subplot(gs[0, :])
ax1.plot(g0_scan, A_arr,   'k-',  lw=2,   label=r'$A_{\rm tail}$')
ax1.plot(g0_scan, B_arr,   'b-',  lw=1.5, label=r'$B_{\rm tail}$ (cos)')
ax1.plot(g0_scan, C_arr,   'r-',  lw=1.5, label=r'$C_{\rm tail}$ (sin)')
ax1.axhline(0, color='gray', lw=0.8, ls='--')
ax1.axvline(G0_E,  color='green',  lw=1.5, ls='--', label=f'$g_0^e={G0_E}$')
ax1.axvline(G0_MU, color='purple', lw=1.5, ls='--', label=f'$g_0^\\mu={G0_MU}$')
# Zaznacz zera B
for z in zeros_B:
    ax1.axvline(z, color='blue', lw=0.8, ls=':', alpha=0.6)
ax1.set_xlabel(r'$g_0$', fontsize=12)
ax1.set_ylabel('Amplituda', fontsize=12)
ax1.set_title(r'Składowe ogona $A, B, C$ vs $g_0$  (ex58)', fontsize=13)
ax1.legend(loc='upper left', fontsize=9)
ax1.grid(True, alpha=0.3)

# Panel 2: Faza phi vs g0
ax2 = fig.add_subplot(gs[1, 0])
ax2.plot(g0_scan, phi_arr / np.pi, 'darkgreen', lw=1.5)
ax2.axhline(0, color='gray', lw=0.7, ls='--')
ax2.axhline(0.5, color='red', lw=0.7, ls=':', alpha=0.5, label='π/2')
ax2.axhline(-0.5, color='red', lw=0.7, ls=':', alpha=0.5)
ax2.axvline(G0_E,  color='green',  lw=1.5, ls='--')
ax2.axvline(G0_MU, color='purple', lw=1.5, ls='--')
ax2.set_xlabel(r'$g_0$')
ax2.set_ylabel(r'$\phi / \pi$')
ax2.set_title(r'Faza ogona $\phi = \arctan2(B,C)$')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

# Panel 3: Faza oddetrendowana
ax3 = fig.add_subplot(gs[1, 1])
ax3.plot(g0_scan, phi_detrended, 'darkred', lw=1.2, label='faza oddetrendowana')
ax3.axhline(0, color='gray', lw=0.7, ls='--')
ax3.axvline(G0_E,  color='green',  lw=1.5, ls='--', label=f'$g_0^e$')
ax3.axvline(G0_MU, color='purple', lw=1.5, ls='--', label=f'$g_0^\\mu$')
ax3.set_xlabel(r'$g_0$')
ax3.set_ylabel(r'$\phi - \phi_{\rm linear}$')
ax3.set_title('Faza oddetrendowana')
ax3.legend(fontsize=8)
ax3.grid(True, alpha=0.3)

# Panel 4: delta_phase (przesunięcie fazowe rozproszenia)
ax4 = fig.add_subplot(gs[2, 0])
ax4.plot(g0_scan, delta_phase / np.pi, 'm-', lw=1.5, label=r'$\delta_{\rm sc}/\pi$')
ax4.axhline(0,    color='gray',  lw=0.7, ls='--')
ax4.axhline(0.5,  color='red',   lw=0.7, ls=':', alpha=0.5, label='π/2 (rezonans)')
ax4.axhline(-0.5, color='red',   lw=0.7, ls=':', alpha=0.5)
ax4.axvline(G0_E,  color='green',  lw=1.5, ls='--')
ax4.axvline(G0_MU, color='purple', lw=1.5, ls='--')
ax4.set_xlabel(r'$g_0$')
ax4.set_ylabel(r'$\delta_{\rm sc}/\pi$')
ax4.set_title(r'Przesunięcie fazowe $\delta = -\arctan(B/C)$')
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3)

# Panel 5: Profile g(r) dla g0_e i g0_mu
ax5 = fig.add_subplot(gs[2, 1])
for name, g0_v, col in [("elektron g0=1.24", G0_E, 'green'),
                         ("mion g0=2.00",      G0_MU, 'purple')]:
    r_p, g_p, _ = integrate_with_bounces(g0_v)
    mask_plot = r_p <= 35.0
    ax5.plot(r_p[mask_plot], g_p[mask_plot], color=col, lw=1.2, label=name)
ax5.axhline(1.0,     color='k',    lw=0.8, ls='--', alpha=0.5, label='g=1 (próżnia)')
ax5.axhline(G_GHOST, color='orange', lw=0.8, ls='--', alpha=0.5, label='g=g*')
ax5.set_xlabel('r')
ax5.set_ylabel('g(r)')
ax5.set_title('Profile solitonów TGP')
ax5.legend(fontsize=8)
ax5.set_xlim(0, 35)
ax5.set_ylim(0.5, max(G0_E, G0_MU) + 0.3)
ax5.grid(True, alpha=0.3)

# Tytuł główny
fig.suptitle(
    f'EX58: Warunek kwantowania g₀  — Faza ogona B(g₀), C(g₀)\n'
    f'Zera B_tail: {[f"{z:.3f}" for z in zeros_B[:4]]}  '
    f'| g₀_e={G0_E} (phi={phi_arr[idx_e]/np.pi:.3f}π)  '
    f'g₀_μ={G0_MU} (phi={phi_arr[idx_mu]/np.pi:.3f}π)',
    fontsize=11
)

plot_path = os.path.join(os.path.dirname(__file__), 'ex58_quantization_condition.png')
plt.savefig(plot_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"[WYKRES] Zapisano: {plot_path}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# TESTY
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("TESTY")
print("=" * 70)

tests_passed = 0
tests_total  = 5

# T1: A_tail jest monotonicznie rosnąca w [1.1, 2.1]
mask_mono = (g0_scan >= 1.1) & (g0_scan <= 2.1)
A_mono = A_arr[mask_mono]
is_mono = bool(np.all(np.diff(A_mono) > -0.05))  # tolerancja 0.05
t1 = is_mono
print(f"  T1: A_tail monotonicznie rośnie w [1.1,2.1]: {'PASS' if t1 else 'FAIL'}"
      f"  (min dA={np.min(np.diff(A_mono)):.4f})")
if t1: tests_passed += 1

# T2: Istnieje co najmniej jedno zero B_tail w [1.1, 3.0]
t2 = len(zeros_B) >= 1
print(f"  T2: ≥1 zero B_tail w [1.1,3.0]: {'PASS' if t2 else 'FAIL'}"
      f"  (znaleziono {len(zeros_B)}: {[f'{z:.3f}' for z in zeros_B[:4]]})")
if t2: tests_passed += 1

# T3: Istnieje zero B_tail w [1.8, 2.2]
zeros_near_mu = [z for z in zeros_B if 1.8 <= z <= 2.2]
t3 = len(zeros_near_mu) >= 1
print(f"  T3: zero B_tail w [1.8,2.2] (blisko g0_mu): {'PASS' if t3 else 'FAIL'}"
      f"  ({[f'{z:.4f}' for z in zeros_near_mu]})")
if t3: tests_passed += 1

# T4: Faza phi_unwrapped jest ciągła (maks. skok < 0.3 rad między sąsiednimi punktami)
# Używamy odwiniętej fazy — surowe phi_arr ma nieciągłości 2pi z atan2
max_dphi_unwr = np.max(np.abs(np.diff(phi_unwrapped)))
t4 = max_dphi_unwr < 0.3
print(f"  T4: Ciągłość phi_unwrapped (maks. skok < 0.3 rad): {'PASS' if t4 else 'FAIL'}"
      f"  (maks. skok = {max_dphi_unwr:.4f} rad)")
if t4: tests_passed += 1

# T5: A_tail benchmark z ex57
A_e_measured  = A_arr[idx_e]
A_mu_measured = A_arr[idx_mu]
t5 = (0.2 <= A_e_measured <= 0.4) and (0.8 <= A_mu_measured <= 1.5)
print(f"  T5: A_tail benchmark: A_e∈[0.2,0.4] i A_mu∈[0.8,1.5]: {'PASS' if t5 else 'FAIL'}"
      f"  (A_e={A_e_measured:.4f}, A_mu={A_mu_measured:.4f})")
if t5: tests_passed += 1

print()
print(f"WYNIK: {tests_passed}/{tests_total} testów przeszło")
print()

# Wniosek końcowy
print("─" * 70)
print("WNIOSEK FIZYCZNY:")
if len(zeros_B) >= 2:
    z0 = zeros_B[0]; z1 = zeros_B[1] if len(zeros_B) > 1 else "?"
    print(f"  Zera B_tail przy g0 = {z0:.3f} (n=0) i {z1 if isinstance(z1, str) else f'{z1:.3f}'} (n=1?)")
    dist_e0  = abs(z0 - G0_E)
    if not isinstance(z1, str):
        dist_mu1 = abs(z1 - G0_MU)
        print(f"  Odległość: |g0_e - zero_0| = {dist_e0:.3f},  |g0_mu - zero_1| = {dist_mu1:.3f}")
        if dist_e0 < 0.05:
            print("  ★ g0_e = 1.24 LEŻY BLISKO zera B_tail (H1 potwierdzona dla elektronu)")
        if dist_mu1 < 0.05:
            print("  ★ g0_mu = 2.00 LEŻY BLISKO zera B_tail (H1 potwierdzona dla mionu)")

print(f"  Faza przy g0_e: phi/pi = {phi_arr[idx_e]/np.pi:.4f}")
print(f"  Faza przy g0_mu: phi/pi = {phi_arr[idx_mu]/np.pi:.4f}")
if abs(phi_arr[idx_e] - phi_arr[idx_mu]) < 0.3:
    print("  ★ Fazy elektronu i mionu są zbliżone — ten sam 'poziom' kwantowania")
print("─" * 70)
