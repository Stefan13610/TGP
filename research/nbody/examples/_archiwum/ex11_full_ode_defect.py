"""
ex11_full_ode_defect.py
========================
Pelne nieliniowe ODE dla defektu topologicznego TGP.

POPRAWKA vs ex10: wlasciwy znak RHS.
  Z funktionalu energii E = int {1/2 f|grad Phi|^2 + V(Phi)} d^3x
  z V(Phi) = (beta/3 Phi_0)*Phi^3 - (gamma/4 Phi_0^2)*Phi^4:

  V'(g*Phi_0) / Phi_0 = beta*g^2 - gamma*g^3 = g^2*(1-g)  [beta=gamma=1]

  => RHS = +g^2*(1-g) > 0 dla g in (0,1)  [pole ROSNIE od g0 do 1]

PELNE ODE (sferyc. symetria, f(g) = 1 + 2*alpha*ln(g), alpha=2):
  f(g)*[g'' + 2*g'/r] + (alpha/g)*(g')^2 = g^2*(1-g)

  => g'' = [g^2*(1-g) - (alpha/g)*(g')^2 - f(g)*(2*g'/r)] / f(g)

GRANICA DUCHOWA: f(g*) = 0 dla g* = exp(-1/(2*alpha)) = exp(-1/4) ~ 0.7788
  Dla g < g*: f(g) < 0  (czlon kinetyczny jest "duchem" - unphysical)
  Obszar fizyczny: g > g* ~ 0.7788

PYTANIE 1: Czy pelne ODE daje ograniczone rozwiazania z g(0) >= g*?
PYTANIE 2: Jaki jest ogon dla duzych r? Yukawa czy oscylacyjny?
PYTANIE 3: Jaki jest C_eff i E_defect dla takich rozwiazan?
"""

import sys, os
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

ALPHA = 2.0
G_GHOST = np.exp(-1.0 / (2.0 * ALPHA))  # ~0.7788

print("=" * 65)
print("Pelne ODE TGP z f(g) != 1")
print("=" * 65)
print(f"alpha = {ALPHA}")
print(f"Granica duchowa: g* = exp(-1/4) = {G_GHOST:.6f}")
print(f"Obszar fizyczny: g > {G_GHOST:.4f}")
print()

# ===========================================================================
# 1. Uproszczone ODE (f=1) z POPRAWIONYM znakiem
# ===========================================================================
print("--- Czesc 1: f=1, RHS = g^2*(1-g)  [poprawiony znak] ---")
print()

def rhs_simple(r, y):
    g, gp = y
    if r < 1e-10:
        g0 = g
        gpp = g0**2 * (1.0 - g0) / 3.0
        return [gp, gpp]
    gpp = g**2 * (1.0 - g) - 2.0 * gp / r
    return [gp, gpp]

g0_vals = [0.80, 0.85, 0.90, 0.95, 0.99]
print(f"  {'g0':>6}  {'g(r=20)':>10}  {'zachowanie':>18}")
print("-" * 40)

results_simple = []
for g0 in g0_vals:
    sol = solve_ivp(rhs_simple, [1e-4, 40.0], [g0, 0.0],
                    method='DOP853', max_step=0.1,
                    rtol=1e-10, atol=1e-13, dense_output=True)
    r_arr = sol.t
    g_arr = sol.y[0]
    gp_arr = sol.y[1]

    idx20 = np.searchsorted(r_arr, 20.0)
    g_at_20 = g_arr[min(idx20, len(g_arr)-1)]

    # Sprawdz czy oscyluje czy monoton.
    crossings = np.sum(np.diff(np.sign(g_arr - 1.0)) != 0)
    if crossings > 2:
        behav = f"oscyluje ({crossings}x)"
    elif g_arr[-1] > 1.001:
        behav = "rosnie > 1"
    elif g_arr[-1] < 0.999 and g_arr[-1] > 0.001:
        behav = "monoton -> 1"
    else:
        behav = "inne"

    results_simple.append({'g0': g0, 'r': r_arr, 'g': g_arr, 'gp': gp_arr})
    print(f"  {g0:>6.3f}  {g_at_20:>10.5f}  {behav:>18}")

print()

# ===========================================================================
# 2. Pelne ODE z f(g) != 1
# ===========================================================================
print("--- Czesc 2: Pelne ODE z f(g) = 1 + 4*ln(g) ---")
print()

def rhs_full(r, y):
    g, gp = y
    fg = 1.0 + 2.0 * ALPHA * np.log(max(g, 1e-8))  # f(g)

    if r < 1e-10:
        if abs(fg) < 1e-10:
            gpp = 0.0
        else:
            gpp = g**2 * (1.0 - g) / (3.0 * fg)
        return [gp, gpp]

    driving = g**2 * (1.0 - g)
    cross   = (ALPHA / max(g, 1e-8)) * gp**2
    damp    = fg * 2.0 * gp / r

    if abs(fg) < 1e-10:
        # Na granicy duchowej: g' wyznaczone z warunku ciaglosci
        gpp = 0.0
    else:
        gpp = (driving - cross - damp) / fg

    return [gp, gpp]

def event_ghost(r, y):
    """Zatrzymaj gdy g osiaga granice duchowa g*."""
    return y[0] - G_GHOST - 0.001  # zatrzymaj przy g = g*+epsilon

event_ghost.terminal = True
event_ghost.direction = -1  # gdy g spada ponizej g*

print(f"  {'g0':>6}  {'g_min':>10}  {'g(r=20)':>10}  {'g<g*?':>8}")
print("-" * 42)

results_full = []
g0_vals_full = [0.85, 0.88, 0.90, 0.92, 0.95, 0.98, 0.99]

for g0 in g0_vals_full:
    fg0 = 1.0 + 2.0 * ALPHA * np.log(g0)
    if fg0 <= 0:
        print(f"  {g0:>6.3f}  g0 < g* (niefizykalny)")
        continue

    sol = solve_ivp(rhs_full, [1e-4, 40.0], [g0, 0.0],
                    method='DOP853', max_step=0.05,
                    rtol=1e-10, atol=1e-13,
                    events=event_ghost, dense_output=True)
    r_arr = sol.t
    g_arr = sol.y[0]

    g_min = g_arr.min()
    idx20 = np.searchsorted(r_arr, 20.0)
    g_at_20 = g_arr[min(idx20, len(g_arr)-1)] if idx20 < len(g_arr) else g_arr[-1]
    hit_ghost = g_min < G_GHOST + 0.002

    results_full.append({'g0': g0, 'r': r_arr, 'g': g_arr})
    print(f"  {g0:>6.3f}  {g_min:>10.5f}  {g_at_20:>10.5f}  {'TAK' if hit_ghost else 'nie':>8}")

print()

# ===========================================================================
# 3. Analiza ogona: Yukawa czy oscylacyjny?
# ===========================================================================
print("--- Czesc 3: Charakter ogona dla duzego r ---")
print()
print("Teoria: zliniearyzowane ODE kolo g=1 (delta = 1-g maly):")
print("  delta'' + (2/r)*delta' = -delta  =>  (Laplacian + 1)*delta = 0")
print("  Rozwiazania: A*sin(r)/r + B*cos(r)/r  (OSCYLACYJNE, nie Yukawa!)")
print()
print("Numeryczna weryfikacja ogona dla g0=0.90 (f=1 i pelne ODE):")
print()

def fit_tail(r_arr, g_arr, r_start=15.0, r_end=30.0):
    """Dopasuj ogon do A*sin(r+phi)/r."""
    mask = (r_arr >= r_start) & (r_arr <= r_end)
    r_fit = r_arr[mask]
    delta_fit = (1.0 - g_arr[mask]) * r_fit  # multiply by r to remove 1/r

    if len(r_fit) < 10:
        return None, None

    # delta * r = A*sin(r+phi) = A*sin(phi)*cos(r) + A*cos(phi)*sin(r)
    # fit: delta*r = a*cos(r) + b*sin(r)
    A_mat = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    try:
        coeffs, _, _, _ = np.linalg.lstsq(A_mat, delta_fit, rcond=None)
        a, b = coeffs
        A_amp = np.sqrt(a**2 + b**2)
        return A_amp, coeffs
    except:
        return None, None

for label, results in [("f=1 (uproszczone)", results_simple),
                        ("Pelne f(g)", results_full)]:
    print(f"  {label}:")
    for res in results:
        g0 = res['g0']
        if g0 not in [0.90, 0.95]:
            continue
        r = res['r']
        g = res['g']
        A_osc, _ = fit_tail(r, g)
        if A_osc is not None:
            # Sprawdz tez Yukawa
            mask = (r >= 15) & (r <= 30)
            r_fit = r[mask]
            delta_fit = (1.0 - g[mask])
            yuk_fit = np.exp(-r_fit) / r_fit
            if len(r_fit) > 5:
                c_yuk = np.dot(delta_fit, yuk_fit) / np.dot(yuk_fit, yuk_fit)
                resid_yuk = np.std(delta_fit - c_yuk * yuk_fit)
                resid_osc = A_osc / (r_fit.mean()) if r_fit.mean() > 0 else 0
                print(f"    g0={g0}: amp_osc={A_osc:.4f}")
                print(f"           C_yukawa_fit={c_yuk:.4f}, resid_yukawa={resid_yuk:.4f}")
    print()

# ===========================================================================
# 4. Droga B: Fenomenologiczny punkt zrodlowy
# ===========================================================================
print("=" * 65)
print("DROGA B: Yukawa jako punkt zrodlowy (fenomenologiczny)")
print("=" * 65)
print()
print("Wniosek z ODE:")
print("  Zliniearyzowane TGP daje ogon OSCYLACYJNY (sin(r)/r),")
print("  nie Yukawa (exp(-r)/r).")
print()
print("  Pelne f(g) nie zmienia tego - ogon jest wciaz oscylacyjny.")
print()
print("DROGA B mowi:")
print("  Akceptujemy ze profil Yukawa to ANSATZ (wejscie fenomenologiczne).")
print("  Uzasadnienie fizyczne: cial sa lokalizowane zaburzenia pola,")
print("  ktore w pewnym sensie 'emituja' Yukawa jako efektywny opis.")
print()
print("  Rownowalnik analogii: elektrodynamika klasyczna nie WYPROWADZA")
print("  ze elektron jest punktem - to jest zalozenie. Teoria jest spojna")
print("  i przewidywalna, mimo ze nie ma rownania dla struktury elektronu.")
print()
print("IMPLEMENTACJA DROGI B w TGP:")
print()
print("  Rownanie pola z zewnetrznym zrodlem:")
print("  (Laplacian - m_sp^2) * delta = -4*pi * C_i * delta^3(x - x_i)")
print()
print("  Rozwiazanie: delta_i(r) = C_i * exp(-m_sp * r) / r  (Yukawa) [ZDEFINIOWANE]")
print()
print("  m_sp = dopasowany do 2. obserwacji (zasieg sily)")
print("  C_i  = m_i * sqrt(G / 4*pi)  [z warunku Newtonowskiego, 1 pomiar]")
print()
print("  Przy tych 2 pomiarach: TGP jest W PELNI PREDYKTYWNA.")
print()

# Pokaz predykcje Drogi B
import numpy as np

G_Newton = 6.674e-11   # m^3 kg^-1 s^-2
m_Planck  = 2.176e-8   # kg
c_factor  = np.sqrt(G_Newton / (4 * np.pi))

print("  Predykcje Drogi B (jedyny wolny param: m_sp):")
print()

masses_kg = {
    'proton':     1.673e-27,
    'Ziemia':     5.972e24,
    'Slonce':     1.989e30,
}

# m_sp scenarios in SI
m_sp_si_vals = {
    '1/AU':   1.0 / 1.496e11,    # inverse AU
    '1/kpc':  1.0 / 3.086e19,    # inverse kpc
    '1/Mpc':  1.0 / 3.086e22,    # inverse Mpc
}

print(f"  {'Cialo':10s}  {'C [SI^1/2 m]':>15s}  C/m_proton")
print("-" * 45)
for name, m in masses_kg.items():
    C_SI = m * c_factor
    ratio = m / 1.673e-27
    print(f"  {name:10s}  {C_SI:>15.4e}  {ratio:.2e}")

print()
print("  Dla kazdego m_sp: efekty trojcialowe, N_crit, promieniowanie")
print("  sa OBLICZALNE DOKLADNIE (z calki Feynmana).")
print()
print("  Jedyna niepewnosc: wartosc m_sp.")
print("  Jedyna obserwacja potrzebna: czy sila TGP zanika na skali lambda=1/m_sp?")

# ===========================================================================
# Wykresy
# ===========================================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Wykres 1: Profile f=1 (poprawiony)
ax = axes[0]
for res in results_simple:
    g0 = res['g0']
    r, g = res['r'], res['g']
    mask = r < 25
    ax.plot(r[mask], g[mask], lw=1.5, label=f"g0={g0:.2f}")
ax.axhline(1.0, color='k', ls='--', lw=1, alpha=0.4)
ax.axhline(G_GHOST, color='r', ls=':', lw=1.5, label=f"g*={G_GHOST:.3f}")
ax.set_xlabel("r", fontsize=11)
ax.set_ylabel("g(r)", fontsize=11)
ax.set_title("ODE (f=1, poprawiony znak)", fontsize=11)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 25)

# Wykres 2: Profile pelne f(g)
ax = axes[1]
for res in results_full:
    g0 = res['g0']
    r, g = res['r'], res['g']
    mask = r < 25
    ax.plot(r[mask], g[mask], lw=1.5, label=f"g0={g0:.2f}")
ax.axhline(1.0, color='k', ls='--', lw=1, alpha=0.4)
ax.axhline(G_GHOST, color='r', ls=':', lw=1.5, label=f"g*={G_GHOST:.3f}")
ax.set_xlabel("r", fontsize=11)
ax.set_ylabel("g(r)", fontsize=11)
ax.set_title(f"Pelne ODE (f(g)=1+4ln(g))", fontsize=11)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 25)

# Wykres 3: Ogon - porownanie Yukawa vs sin(r)/r
ax = axes[2]
r_tail = np.linspace(5, 30, 300)
ax.plot(r_tail, np.sin(r_tail) / r_tail, 'b-', lw=2, label="sin(r)/r (oscylacyjny)")
ax.plot(r_tail, np.exp(-r_tail) / r_tail * 200, 'r--', lw=2, label="exp(-r)/r * 200 (Yukawa)")
if results_simple:
    res = [x for x in results_simple if x['g0'] == 0.95][0]
    r, g = res['r'], res['g']
    mask = r > 5
    delta = 1.0 - g[mask]
    # normalize
    i5 = np.argmin(np.abs(r[mask] - 5))
    norm = delta[i5] / (np.sin(5)/5) if abs(np.sin(5)/5) > 1e-8 else 1
    ax.plot(r[mask][:200], delta[:200]/max(abs(delta[:50]))*0.08,
            'g-', lw=2, label="TGP g0=0.95 (numeryczny)")
ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel("r", fontsize=11)
ax.set_ylabel("delta = 1-g(r)", fontsize=11)
ax.set_title("Ogon: oscylacyjny vs Yukawa", fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(5, 30)

plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "ex11_full_ode.png")
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print()
print(f"Wykres: {out}")
