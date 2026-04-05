"""
p6_ode_shooting_generacje.py
============================
POPRAWNE strzałkowanie ODE dla profili generacji w TGP.

PROBLEM z integracją inward (R_max → a_gam):
  Tryb rosnący e^{+m_eff r}/r jest NIESTABILNY inward —
  nawet drobny błąd IC eksploduje jak e^{m_eff*(R-r)}.
  Dla m_eff=0.323, R=30: czynnik wzmocnienia ≈ e^{9.98} ≈ 21 500.

POPRAWNA metoda (ten skrypt): integracja OUTWARD z rdzenia.

SFORMUŁOWANIE:
  Soliton TGP — profil φ(r) dla r > a_gam:
  (1 + α/φ)φ'' + (2/r)(1+α/φ)φ' - (α/2φ²)(φ')² = V'(φ)

  Warunek Neumanna na rdzeniu (źródło):
    φ'(a_gam) = -K / a_gam²
  (K = "ładunek" cząstki, stałe)

  Warunek asymptotyczny:
    φ(r→∞) → 1 (próżnia TGP)

  Samospójność:
    E[φ] = 4π·K  (energia pola = masa cząstki)

ALGORYTM:
  Dla danego K:
  1. Ustal IC: φ'(a_gam) = -K/a_gam²
  2. Strzelaj po φ(a_gam) = ψ_core (nieznane) — integracja outward
  3. Znajdź ψ_core* takie że φ(R_max) → 1  (warunek próżni)
  4. Oblicz E[φ_{ψ_core*}]

  Samospójność: g(K) = E(K) / (4π·K) - 1 = 0 — szukaj zer K.

WYNIKI:
  K₁ ≈ 0.0098 (ansatz): czy pełne ODE potwierdza?
  K₂ ≈ 2.033  (ansatz): czy profil φ → 1 asymptotycznie?
  K₃ ≈ 34.14  (ansatz): czy istnieje przy λ = λ*?
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ============================================================
# PARAMETRY
# ============================================================
ALPHA  = 8.5616
A_GAM  = 0.040
LAM    = 5.501357e-06
PHI0   = 1.0
GAMMA  = 1.0

M_EFF  = 1.0 / np.sqrt(1.0 + ALPHA)  # 0.323 — prawdziwa masa ekranowania
R_MAX  = 50.0  # dalekie tło — cel: φ(R_MAX) ≈ 1

print(f"Parametry: α={ALPHA}, a_Γ={A_GAM}, λ={LAM:.4e}")
print(f"Prawdziwa masa ekranowania:  m_eff = {M_EFF:.4f}  (Yukawa używa m_sp=1)")
print(f"R_max = {R_MAX},  exp(-m_eff·R_max) = {np.exp(-M_EFF*R_MAX):.2e}")
print()


def V_mod(phi):
    return GAMMA/3*phi**3 - GAMMA/4*phi**4 + LAM/6*(phi-1)**6

def dV_mod(phi):
    return GAMMA*phi**2 - GAMMA*phi**3 + LAM*(phi-1)**5

V1 = V_mod(1.0)


# ============================================================
# PEŁNE ODE (outward)
# ============================================================
def ode_rhs(r, y):
    """φ'' = V'(φ)/(1+α/φ) + α(φ')²/(2φ²(1+α/φ)) - (2/r)φ'"""
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA / phi
    ddphi = (dV_mod(phi) / kfac
             + ALPHA * dphi**2 / (2.0 * phi**2 * kfac)
             - (2.0 / r) * dphi)
    return [dphi, ddphi]


# ============================================================
# STRZAŁKOWANIE OUTWARD: znajdź ψ_core(K) takie że φ(R) → 1
# ============================================================
def phi_at_rmax(psi_core, K, r_max=R_MAX, n_eval=3000):
    """
    Integruje ODE outward od a_gam do r_max z IC:
      φ(a_gam) = psi_core
      φ'(a_gam) = -K / a_gam²

    Zwraca φ(r_max) — cel: 1.0.
    """
    phi0  = psi_core
    dphi0 = -K / A_GAM**2

    # Siatka logarytmiczna (gęściej przy rdzeniu)
    r_eval = A_GAM * (r_max / A_GAM)**np.linspace(0, 1, n_eval)

    # Event: zatrzymaj jeśli phi < 0 lub phi > 1e6
    def ev_lo(r, y): return y[0] - 1e-5
    ev_lo.terminal = True; ev_lo.direction = -1
    def ev_hi(r, y): return y[0] - 1e6
    ev_hi.terminal = True; ev_hi.direction = 1

    sol = solve_ivp(
        ode_rhs, [A_GAM, r_max], [phi0, dphi0],
        method='DOP853', rtol=1e-10, atol=1e-12,
        t_eval=r_eval, events=[ev_lo, ev_hi]
    )
    if len(sol.t) == 0:
        return np.nan

    phi_end = sol.y[0, -1]
    return phi_end


def find_psi_core(K, psi_lo=None, psi_hi=None):
    """
    Dla danego K: znajdź ψ_core takie że φ(R_MAX) = 1.

    Strategia:
    - Dużo ψ_core → φ rośnie i φ(R) >> 1
    - Małe ψ_core → φ maleje i φ(R) < 1 lub < 0

    Szukamy ψ_core* takie że F(ψ_core) ≡ φ(R_MAX) - 1 = 0.
    """
    if psi_lo is None:
        psi_lo = 1.001          # lekko powyżej próżni
    if psi_hi is None:
        psi_hi = 1.0 + K / A_GAM * 2  # górna granica (2× Yukawa core)

    # Sprawdź znaki na końcach
    F_lo = phi_at_rmax(psi_lo, K) - 1.0
    F_hi = phi_at_rmax(psi_hi, K) - 1.0

    if not (np.isfinite(F_lo) and np.isfinite(F_hi)):
        return None, None

    # Szukaj przedziału ze zmianą znaku
    if F_lo * F_hi > 0:
        # Skan wewnętrzny
        psi_vals = np.linspace(psi_lo, psi_hi, 30)
        F_vals   = np.array([phi_at_rmax(p, K) - 1.0 for p in psi_vals])
        for i in range(len(F_vals)-1):
            if np.isfinite(F_vals[i]) and np.isfinite(F_vals[i+1]):
                if F_vals[i] * F_vals[i+1] < 0:
                    psi_lo = psi_vals[i]
                    psi_hi = psi_vals[i+1]
                    F_lo   = F_vals[i]
                    F_hi   = F_vals[i+1]
                    break
        else:
            return None, None

    # Bisekcja
    try:
        psi_core = brentq(
            lambda p: phi_at_rmax(p, K) - 1.0,
            psi_lo, psi_hi,
            xtol=1e-6, rtol=1e-6, maxiter=50
        )
        return psi_core, phi_at_rmax(psi_core, K)
    except Exception:
        return None, None


# ============================================================
# ENERGIA profilu ze znalezionym ψ_core
# ============================================================
def profile_energy(psi_core, K, n_eval=5000):
    """Oblicza E[φ] i K_eff z profilu z IC (psi_core, -K/a_gam²)."""
    dphi0 = -K / A_GAM**2
    r_eval = A_GAM * (R_MAX / A_GAM)**np.linspace(0, 1, n_eval)

    sol = solve_ivp(
        ode_rhs, [A_GAM, R_MAX], [psi_core, dphi0],
        method='DOP853', rtol=1e-10, atol=1e-12,
        t_eval=r_eval
    )
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]

    K_eff = r[0]**2 * abs(dphi[0])  # = K (powinno być zgodne)

    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
    return Ek+Ep, Ek, Ep, K_eff, sol


# ============================================================
# SAMOSPÓJNOŚĆ g(K)
# ============================================================
def g_ode_full(K, verbose=False):
    """
    g(K) = E[φ_{ψ_core(K)}] / (4π·K) - 1

    Dla każdego K: wyznacza ψ_core tak by φ(R)→1, potem liczy E.
    """
    psi_core, phi_end = find_psi_core(K)
    if psi_core is None:
        if verbose: print(f"  K={K:.4e}: brak ψ_core (zmiana znaku nieznaleziona)")
        return np.nan, None

    if abs(phi_end - 1.0) > 0.01:
        if verbose: print(f"  K={K:.4e}: φ(R)={phi_end:.4f} (słaba zbieżność)")

    E, Ek, Ep, K_eff, sol = profile_energy(psi_core, K)
    g = E / (4*np.pi*K) - 1.0

    if verbose:
        print(f"  K={K:.4e}: ψ_core={psi_core:.4f}, φ(R)={phi_end:.5f}, "
              f"E={E:.4e}, g={g:.4e}")
    return g, {'psi_core': psi_core, 'E': E, 'Ek': Ek, 'Ep': Ep, 'K_eff': K_eff}


# ============================================================
# GŁÓWNA ANALIZA
# ============================================================
print("=" * 65)
print("P6: POPRAWNE STRZAŁKOWANIE OUTWARD — profile generacji TGP")
print("=" * 65)
print()

# --- Weryfikacja dla K_Yukawa ---
print("KROK 1: Weryfikacja dla znanych K_Yukawa:")
print()
yukawa_Ks = [0.009820, 2.032728, 34.14450]
yukawa_labels = ['K₁ (elektron)', 'K₂ (mion)', 'K₃ (taon?)']

for K_y, label in zip(yukawa_Ks, yukawa_labels):
    print(f"  {label}: K_Yukawa = {K_y:.5f}")
    g, info = g_ode_full(K_y, verbose=True)
    if info:
        print(f"    E/(4πK)-1 = {g:.4e}  {'✓ samospójny' if abs(g) < 0.1 else '✗ niespójny'}")
    print()

# --- Skan g_ODE(K) ---
print("KROK 2: Skan g_ODE(K):")
print()
K_scan = np.logspace(-3, 2, 50)
g_scan = []
info_scan = []

print(f"  {'K':>10}  {'ψ_core':>10}  {'φ(R_max)':>10}  {'g_ODE':>12}")
print("  " + "-"*48)
for K in K_scan:
    g, info = g_ode_full(K)
    g_scan.append(g)
    info_scan.append(info)
    psi_c   = info['psi_core'] if info else None
    phi_end = phi_at_rmax(psi_c, K) if psi_c else None
    if info:
        print(f"  {K:>10.4e}  {psi_c:>10.4f}  {phi_end:>10.5f}  {g:>12.4e}")
    else:
        print(f"  {K:>10.4e}  {'---':>10}  {'---':>10}  {'nan':>12}")

g_scan = np.array(g_scan)

# Znajdź zmiany znaku
zero_intervals = []
for i in range(len(g_scan) - 1):
    gi, gj = g_scan[i], g_scan[i+1]
    if np.isfinite(gi) and np.isfinite(gj) and gi * gj < 0:
        zero_intervals.append((K_scan[i], K_scan[i+1]))

print()
print(f"  Znaleziono {len(zero_intervals)} przedziałów ze zmianą znaku:")
for a_lo, a_hi in zero_intervals:
    print(f"    K ∈ [{a_lo:.4e}, {a_hi:.4e}]")

# Bisekcja
ode_roots = []
for K_lo, K_hi in zero_intervals:
    try:
        K_root = brentq(
            lambda K: g_ode_full(K)[0],
            K_lo, K_hi, xtol=1e-6, rtol=1e-6, maxiter=30
        )
        g_root, info_root = g_ode_full(K_root, verbose=False)
        ode_roots.append({'K': K_root, 'info': info_root})
    except Exception as e:
        print(f"    Błąd bisekcji: {e}")

print()
print(f"  WYNIKI: {len(ode_roots)} samospójnych solitonów ODE")
print()
print(f"  {'Nr':>3}  {'K_ODE':>12}  {'ψ_core':>10}  {'E':>12}  {'g':>10}")
print("  " + "-"*55)
for i, root in enumerate(ode_roots, 1):
    K     = root['K']
    info  = root['info']
    psi_c = info['psi_core'] if info else 0
    E     = info['E']        if info else 0
    g, _  = g_ode_full(K)
    print(f"  {i:>3}  {K:>12.5e}  {psi_c:>10.3f}  {E:>12.4e}  {g:>10.4e}")

print()

# Ratios
if len(ode_roots) >= 2:
    K1 = ode_roots[0]['K']
    K2 = ode_roots[1]['K']
    print(f"  ODE r₂₁ = K₂/K₁ = {K2/K1:.2f}  (cel: 207)")
if len(ode_roots) >= 3:
    K1 = ode_roots[0]['K']
    K3 = ode_roots[2]['K']
    print(f"  ODE r₃₁ = K₃/K₁ = {K3/K1:.2f}  (cel: 3477)")
print()


# ============================================================
# PROFILE WIZUALIZACJA
# ============================================================
print("Generowanie wykresu...")

fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.suptitle(f'TGP p6: Outward ODE vs Yukawa  (α={ALPHA}, a_Γ={A_GAM}, λ={LAM:.2e})',
             fontsize=11, fontweight='bold')

colors = ['blue', 'green', 'red']

# Panel 1: g_ODE(K) skan
ax = axes[0, 0]
mask = np.isfinite(g_scan)
ax.axhline(0, color='black', lw=0.7, linestyle='--')
if mask.any():
    ax.plot(K_scan[mask], g_scan[mask], 'b-', lw=1.5, label='g_ODE(K)')
    ax.scatter([r['K'] for r in ode_roots], [0]*len(ode_roots),
               s=100, color='red', zorder=5, label=f'Zera ODE ({len(ode_roots)})')
# Yukawa zeros dla porównania
ax.scatter(yukawa_Ks[:2], [0, 0], s=80, color='orange', marker='^',
           zorder=5, label='K_Yukawa (K₁,K₂)')
ax.set_xscale('log')
ax.set_xlabel('K')
ax.set_ylabel('g_ODE(K) = E/(4πK) − 1')
ax.set_title('Samospójność pełnego ODE')
ax.legend(fontsize=8)
ax.set_ylim(-2, 5)
ax.grid(True, alpha=0.3)

# Panel 2: Profile ODE vs Yukawa dla K₁
ax = axes[0, 1]
r_plot = A_GAM * (R_MAX/A_GAM)**np.linspace(0, 1, 1000)
for K_y, col, lab in zip(yukawa_Ks[:2], colors[:2], ['K₁', 'K₂']):
    phi_y = 1.0 + K_y * np.exp(-r_plot) / r_plot
    ax.plot(r_plot, phi_y, '--', color=col, lw=1, alpha=0.6, label=f'Yukawa {lab}')
for root, col, lab in zip(ode_roots[:2], colors[:2], ['ODE₁', 'ODE₂']):
    psi_c = root['info']['psi_core']
    K_r   = root['K']
    _, _, _, _, sol = profile_energy(psi_c, K_r, n_eval=1000)
    ax.plot(sol.t, sol.y[0], '-', color=col, lw=1.5, label=f'{lab}')
ax.axhline(1, color='black', lw=0.5, linestyle=':')
ax.set_xscale('log')
ax.set_xlabel('r')
ax.set_ylabel('φ(r)')
ax.set_title('Profile solitonów ODE (—) vs Yukawa (--)')
ax.legend(fontsize=8)
ax.set_ylim(0.9, 2.0)
ax.grid(True, alpha=0.3)

# Panel 3: ψ_core(K) — ODE vs Yukawa
ax = axes[1, 0]
K_arr_plot = K_scan[mask]
psi_c_ode  = np.array([info['psi_core'] if info else np.nan for info in info_scan])[mask]
psi_c_yuk  = np.array([1 + K / A_GAM for K in K_scan])[mask]  # Yukawa approx.
ax.loglog(K_arr_plot, psi_c_ode, 'b-', lw=1.5, label='ψ_core ODE')
ax.loglog(K_scan[mask], psi_c_yuk, 'g--', lw=1, alpha=0.7, label='ψ_core Yukawa ≈ 1+K/a_Γ')
ax.set_xlabel('K')
ax.set_ylabel('ψ_core')
ax.set_title('ψ_core(K): ODE vs Yukawa')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 4: Diagnoza — φ(R_max) - 1 dla K_Y Yukawa z różnymi ψ_core
ax = axes[1, 1]
K_test = yukawa_Ks[1]  # K₂ jako przykład
psi_range = np.linspace(1.001, 1 + 2*K_test/A_GAM, 50)
phi_ends  = np.array([phi_at_rmax(p, K_test) - 1.0 for p in psi_range])
ax.plot(psi_range, phi_ends, 'b-', lw=1.5)
ax.axhline(0, color='red', lw=1, linestyle='--', label='BC: φ(R)=1')
ax.axvline(1 + K_test/A_GAM, color='orange', lw=1, linestyle=':',
           label=f'ψ_core Yukawa = {1+K_test/A_GAM:.1f}')
ax.set_xlabel('ψ_core')
ax.set_ylabel('φ(R_max) - 1')
ax.set_title(f'Shooting dla K₂={K_test:.3f}: szukaj φ(R)=1')
ax.legend(fontsize=9)
ax.set_ylim(-5, 20)
ax.grid(True, alpha=0.3)

plt.tight_layout()
out_path = __file__.replace('.py', '.png')
plt.savefig(out_path, dpi=120, bbox_inches='tight')
plt.close()
print(f"  Zapisano: {out_path}")
print()


# ============================================================
# PODSUMOWANIE
# ============================================================
print("=" * 65)
print("PODSUMOWANIE p6: OUTWARD ODE SHOOTING")
print()
print(f"  Liczba samospójnych solitonów (pełne ODE): {len(ode_roots)}")
print()
if len(ode_roots) == 0:
    print("  g_ODE(K) nie ma zer — żaden soliton Yukawa nie przeżywa")
    print("  w pełnym ODE z poprawnymi BC. Możliwe przyczyny:")
    print("  a) g_ODE > 0 wszędzie: E > 4πK dla wszystkich K")
    print("     → brak równowagi pole↔cząstka")
    print("  b) ψ_core(K) nieznalezione: brak profilu spełniającego φ(∞)=1")
    print("     → solitony Yukawa są artefaktami ansatzu, nie rozwiązaniami ODE")
elif len(ode_roots) >= 1:
    for i, root in enumerate(ode_roots, 1):
        print(f"  Soliton {i}: K={root['K']:.5e}")
    if len(ode_roots) >= 2:
        K1, K2 = ode_roots[0]['K'], ode_roots[1]['K']
        print(f"  r₂₁(ODE) = {K2/K1:.2f}  (Yukawa: 207)")
    if len(ode_roots) >= 3:
        K1, K3 = ode_roots[0]['K'], ode_roots[2]['K']
        print(f"  r₃₁(ODE) = {K3/K1:.2f}  (Yukawa: 3477)")
print()
g_finite = g_scan[np.isfinite(g_scan)]
if len(g_finite) > 0:
    print(f"  g_ODE min = {g_finite.min():.4e}")
    print(f"  g_ODE max = {g_finite.max():.4e}")
    print()
    if g_finite.min() > 0:
        print("  g_ODE > 0 wszędzie: E > 4πK dla wszystkich K")
        print("  → Pełne ODE z niestandardowym term kinetycznym (1+α/φ)")
        print("    daje WIĘCEJ energii niż ansatz Yukawa przewiduje.")
        print("  → Brak samospójnych solitonów w pełnym ODE (przy λ=λ*)?")
        print()
        print("  Możliwa interpretacja:")
        print("  1. Cząstki TGP to quasi-cząstki w tle Φ₀, a nie izolowane solitony")
        print("  2. Samospójność wymaga Φ_bg ≠ Φ₀ (brane tło zależne od M)")
        print("  3. Właściwe równanie pola wymaga pełnego kowariantnego opisu")
        print("     (efekty geometryczne pola krzywizny nie są w uproszczonym ODE)")
