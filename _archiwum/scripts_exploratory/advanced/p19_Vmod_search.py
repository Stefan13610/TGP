"""
p19_Vmod_search.py
==================
Szukamy modyfikacji V_mod ktora daje r21=207 w pelnym ODE.

PROBLEM:
  V_mod(phi) = phi^3/3 - phi^4/4 + lam*(phi-1)^6/6
  Dla phi >> 1: V_mod ~ -phi^4/4 -> -inf
  To powoduje E_p << 0 dla duzych phi_core -> g = E/(4*pi*K) - 1 << 0 dla K>>K*1
  Wynik: tylko JEDEN soliton K*1=0.010414 (p14c).

POMYSL:
  Dodamy czlon PODNOSNIKOWY V_add(phi) ktory:
  a) V_add(1) = 0           (nie zmienia prozni)
  b) V'_add(1) = 0          (nie zmienia warunku prozni)
  c) V_add(phi) > 0 dla phi != 1  (stabilizuje)
  d) V_add(phi) >> phi^4/4 dla phi > phi_cross  (pokonuje ucieczke)

KANDYDAT:
  V_add(phi) = mu * (phi-1)^4  (czlon kwartyczny)

  V_add(1)  = 0              [OK]
  V'_add(1) = 0              [OK]
  V''_add(1) = 12*mu > 0     [zmienia m_eff!]

  Crossover z -phi^4/4:
  mu*(phi-1)^4 ~ phi^4/4 gdy mu ~ 1/4 * (phi/(phi-1))^4 ~ 1/4 dla phi >> 1
  => Potrzeba mu >= 1/4 by zatrzymac ucieczkę dla phi >> 1.

ALTERNATYWA:
  V_add(phi) = mu*(phi-1)^4 / (1 + nu*(phi-1)^2)^2  (nasycajace)
  albo: V_mod_new = phi^3/3 - mu*phi^4/4 + lam*(phi-1)^6/6  z mu < 1
  Dla mu=0: V_mod_new = phi^3/3 + lam*(phi-1)^6/6 -> STABILNY (brak phi^4 term)

  Ale czy mu < 1 zmienia V'(1)? V'_mod = phi^2 - phi^3 + lam*(phi-1)^5
  Z mu: V'_mod = phi^2 - mu*phi^3 + lam*(phi-1)^5
  V'(1) = 1 - mu + 0 = 1 - mu
  Dla mu=1: V'(1)=0 [OK, oryginalny]. Dla mu!=1: V'(1)!=0 [ZEPSUTE].

ROZWIAZANIE POPRAWNE:
  Parametryzacja: V = gamma/3*phi^3 - gamma/4*phi^4 + eta*(phi-1)^4/4 + lam*(phi-1)^6/6

  Sprawdzenie V'(1):
  V'(phi) = gamma*phi^2 - gamma*phi^3 + eta*(phi-1)^3 + lam*(phi-1)^5
  V'(1)   = gamma - gamma + 0 + 0 = 0  [ZAWSZE OK!]

  V''(phi) = 2*gamma*phi - 3*gamma*phi^2 + 3*eta*(phi-1)^2 + 5*lam*(phi-1)^4
  V''(1)   = 2*gamma - 3*gamma + 0 + 0 = -gamma < 0   <- PROBLEM!

  Hmm, phi=1 to maksimum nie minimum dla V_mod bazowego!

  Warunek Neumanna (a_Gam) stabilizuje, nie V''.
  Dla malych K linearyzacja daje oscylacje z m^2 = V''(1)/(1+alpha) < 0
  ale to jest PROBLEM OSOBNY od solitonu.

  Faktycznie V_mod ma phi=1 jako WROG (V'(1)=0, V''(1)<0) = punkt siodlowy,
  ale warunek brzegowy phi(inf)=1 jest STABILNY numerycznie bo
  oscylacje rozkładaja się eksponencjalnie przy m^2 < 0:
  m_eff = sqrt(-V''(1)/(1+alpha)) = sqrt(1/(1+alpha)) = M_EFF (masa tachyonowa???)

  To jest TACHYON? Nie - TGP to model z warunkiem brzegowym phi->1 jako
  minimum energii FREE (bez kinetyki). Ale V''(1) = -gamma < 0...

  Faktycznie: V_mod(psi) = psi^3/3 - psi^4/4, V'(1)=0 ale V''(1) = 2-3 = -1 < 0
  To oznacza ze psi=1 jest MAKSIMUM localnym V_mod, nie minimum!

  ALE: to jest V(psi) - V1 nie energia PELNA. Energia calkowita obejmuje
  tez kinetyczny i gradient terms. Minimum energii dla phi=const to phi=0 (V=0).
  Rozwiazanie phi=1 jest stabilne TYLKO dzieki warunkom brzegowym (phi->1 at inf).

  Dlatego solitony "siedzą" przy phi~1 pomimo ze V''(1)<0.

PLAN p19:
  Dodajemy V_add(phi) = eta*(phi-1)^4/4 z eta > 0.
  To modyfikuje tylko V dla phi != 1.

  Dla duzego phi:
  V_mod_new ~ -phi^4/4 + eta*phi^4/4 = (eta-1)/4 * phi^4

  Dla eta > 1: V_mod_new ~ (eta-1)/4 * phi^4 -> +inf (stabilne!)
  Dla eta = 1: V_mod_new ~ -phi^4/4 + phi^4/4 = 0 + phi^3/3 (kubiczne)
  Dla eta < 1: V_mod_new ~ -(1-eta)/4 * phi^4 -> -inf (wciaz niestabilne)

  => KLUCZOWY PROG: eta = 1.
  Dla eta >= 1 potencjal jest co najmniej kubiczny dla duzych phi.

  SKAN: eta od 0 do 2, szukamy kiedy pojawia sie drugi soliton w pelnym ODE.

UWAGA: zmiana eta modyfikuje m_eff!
  V_add''(phi) = 3*eta*(phi-1)^2
  V_add''(1) = 0  [OK - nie zmienia m_eff przy phi=1]

  => m_eff = sqrt(-V''(1)/(1+alpha)) = sqrt(1/(1+alpha)) NIEZMIENNE dla wszystkich eta!
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

ALPHA  = 8.5616
LAM    = 5.501357e-06
GAMMA  = 1.0
M_EFF  = 1.0 / np.sqrt(1.0 + ALPHA)
R_MAX_BASE = 80.0
V1     = GAMMA/3 - GAMMA/4   # V_mod(psi=1)

print("P19: Szukamy V_mod z r21=207 w pelnym ODE")
print("=" * 65)
print(f"  Bazowy ODE: alpha={ALPHA}, lam={LAM:.3e}")
print(f"  M_eff = {M_EFF:.4f}")
print()
print("  Modyfikacja: V_mod_new = V_mod + eta*(phi-1)^4/4")
print("  Dla phi>>1: V_mod_new ~ (eta-1)/4 * phi^4")
print("  Prog stabilnosci: eta = 1")
print()

# ============================================================
# ODE z V_mod_new(phi; eta)
# ============================================================

def make_ode(eta, lam=LAM, gamma=GAMMA, alpha=ALPHA):
    """Zwraca funkcje ODE dla danego eta."""
    def V_new(p):
        return gamma/3*p**3 - gamma/4*p**4 + lam/6*(p-1)**6 + eta/4*(p-1)**4
    def dV_new(p):
        return gamma*p**2 - gamma*p**3 + lam*(p-1)**5 + eta*(p-1)**3
    def ode(r, y):
        phi, dphi = y
        phi = max(phi, 1e-10)
        kfac = 1.0 + alpha / phi
        ddphi = (dV_new(phi)/kfac
                 + alpha*dphi**2/(2.0*phi**2*kfac)
                 - (2.0/r)*dphi)
        return [dphi, ddphi]
    return ode, V_new, dV_new


def phi_at_rmax(psi_core, K, a_gam, r_max, eta, n_eval=2000):
    """Strzal na zewnatrz: phi(r_max) dla danych ICs."""
    ode, _, _ = make_ode(eta)
    dphi0 = -K / a_gam**2
    def ev_lo(r, y): return y[0] - 1e-5
    ev_lo.terminal = True; ev_lo.direction = -1
    def ev_hi(r, y): return y[0] - 1e7
    ev_hi.terminal = True; ev_hi.direction = 1
    r_eval_raw = a_gam * (r_max/a_gam)**np.linspace(0, 1, n_eval)
    r_eval = np.clip(r_eval_raw, a_gam*(1+1e-12), r_max*(1-1e-12))
    sol = solve_ivp(ode, [a_gam, r_max], [psi_core, dphi0],
                    method='DOP853', rtol=1e-9, atol=1e-11,
                    t_eval=r_eval, events=[ev_lo, ev_hi])
    if len(sol.t) == 0 or sol.t[-1] < r_max * 0.99:
        return np.nan
    return sol.y[0, -1]


def compute_g(psi_core, K, a_gam, eta, n_eval=5000, lam=LAM):
    """Oblicz g = E/(4*pi*K) - 1 dla danego stanu."""
    ode, V_new, _ = make_ode(eta)
    dphi0 = -K / a_gam**2
    r_max = max(R_MAX_BASE, 20*a_gam, 30/M_EFF)
    r_eval_raw = a_gam * (r_max/a_gam)**np.linspace(0, 1, n_eval)
    r_eval = np.clip(r_eval_raw, a_gam*(1+1e-12), r_max*(1-1e-12))
    sol = solve_ivp(ode, [a_gam, r_max], [psi_core, dphi0],
                    method='DOP853', rtol=1e-9, atol=1e-11, t_eval=r_eval)
    r   = sol.t
    phi = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_new(phi)-V1)*r**2, r)
    return (Ek+Ep)/(4*np.pi*K) - 1.0


def find_psi_core_shoot(K, a_gam, eta, r_max=None, n_scan=50):
    """Znajdz psi_core z warunku phi(rmax)=1."""
    if r_max is None:
        r_max = max(R_MAX_BASE, 20*a_gam, 30/M_EFF)

    psi_guess = 1.0 + K/a_gam * np.exp(-M_EFF*a_gam)
    best_psi = np.nan

    for window in [0.5, 1.5, 3.0, 6.0]:
        psi_lo = max(1.001, psi_guess - window)
        psi_hi = psi_guess + window
        psi_scan = np.linspace(psi_lo, psi_hi, n_scan)
        F_scan = [phi_at_rmax(p, K, a_gam, r_max, eta) - 1.0
                  for p in psi_scan]
        for i in range(len(psi_scan)-1):
            if (np.isfinite(F_scan[i]) and np.isfinite(F_scan[i+1])
                    and F_scan[i]*F_scan[i+1] < 0):
                try:
                    best_psi = brentq(
                        lambda p: phi_at_rmax(p, K, a_gam, r_max, eta) - 1.0,
                        psi_scan[i], psi_scan[i+1],
                        xtol=1e-5, maxiter=20)
                except Exception:
                    pass
                if not np.isnan(best_psi): break
        if not np.isnan(best_psi): break
    return best_psi


def g_on_schw_line(K, C, eta):
    """g(K) na linii Schwarzschilda a_Gam = C*K."""
    a_gam = C * K
    psi_core = find_psi_core_shoot(K, a_gam, eta)
    if np.isnan(psi_core):
        return np.nan
    return compute_g(psi_core, K, a_gam, eta)


# ============================================================
# KROK 1: Jak eta zmienia g(K) na linii Schwarzschilda?
# ============================================================
C_ref = 3.841   # stala Schwarzschilda z K*1/a_Gam1
K_scan = np.logspace(-2, 1.0, 30)

eta_vals = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]

print("KROK 1: g(K) dla roznych eta na linii a_Gam = 3.841*K")
print("-" * 65)

g_results = {}
psi_results = {}

for eta in eta_vals:
    print(f"\n  eta = {eta:.1f}:")
    g_list  = []
    psi_list = []
    for K in K_scan:
        a_gam = C_ref * K
        psi_c = find_psi_core_shoot(K, a_gam, eta)
        if np.isnan(psi_c):
            g_list.append(np.nan); psi_list.append(np.nan)
            continue
        g_val = compute_g(psi_c, K, a_gam, eta)
        g_list.append(g_val); psi_list.append(psi_c)
        tag = ''
        if np.isfinite(g_val) and abs(g_val) < 0.05: tag = '  <-- ZERO!'
        print(f"    K={K:.4f}: psi={psi_c:.4f}, g={g_val:.4e}{tag}")

    g_arr = np.array(g_list)
    g_results[eta]   = g_arr
    psi_results[eta] = np.array(psi_list)

    # Znajdz zera
    zeros = []
    for i in range(len(K_scan)-1):
        gi, gj = g_arr[i], g_arr[i+1]
        psi_i = psi_results[eta][i]
        psi_j = psi_results[eta][i+1]
        if not (np.isfinite(gi) and np.isfinite(gj)): continue
        if gi*gj < 0:
            dpsi = abs(psi_j - psi_i) if (not np.isnan(psi_i) and not np.isnan(psi_j)) else 999
            if dpsi > 1.5: continue  # skok galezi
            K_zero = (K_scan[i]*K_scan[i+1])**0.5
            zeros.append(K_zero)
    print(f"    Zera g(K): {[f'{z:.4e}' for z in zeros]}")
    if len(zeros) >= 2:
        r21 = zeros[1]/zeros[0]
        print(f"    K*2/K*1 = {r21:.2f}")

# ============================================================
# KROK 2: Precyzyjne wyznaczenie eta* dajacego r21=207
# ============================================================
print()
print("=" * 65)
print("KROK 2: Szukamy eta* gdzie K*2/K*1 = 207")
print("=" * 65)

def compute_ratio_for_eta(eta, C=C_ref, K_scan_loc=None, verbose=False):
    """Oblicz stosunek K*2/K*1 dla danego eta."""
    if K_scan_loc is None:
        K_scan_loc = np.logspace(-2.5, 2.0, 45)

    g_vals = []
    psi_vals = []
    for K in K_scan_loc:
        a_gam = C * K
        psi_c = find_psi_core_shoot(K, a_gam, eta)
        if np.isnan(psi_c):
            g_vals.append(np.nan); psi_vals.append(np.nan)
            continue
        g_v = compute_g(psi_c, K, a_gam, eta)
        g_vals.append(g_v); psi_vals.append(psi_c)

    g_arr   = np.array(g_vals)
    psi_arr = np.array(psi_vals)

    zeros = []
    for i in range(len(K_scan_loc)-1):
        gi, gj = g_arr[i], g_arr[i+1]
        if not (np.isfinite(gi) and np.isfinite(gj)): continue
        if gi*gj < 0:
            dpsi = abs(psi_arr[i]-psi_arr[i+1]) if (not np.isnan(psi_arr[i]) and not np.isnan(psi_arr[i+1])) else 999
            if dpsi > 1.5: continue
            K_z = (K_scan_loc[i]*K_scan_loc[i+1])**0.5
            zeros.append(K_z)
            if verbose: print(f"      K*~{K_z:.4e}")

    if len(zeros) < 2:
        return np.nan, zeros, g_arr, K_scan_loc

    ratio = zeros[1] / zeros[0]
    if verbose:
        print(f"    K*1~{zeros[0]:.4e}, K*2~{zeros[1]:.4e}, r21={ratio:.2f}")
    return ratio, zeros, g_arr, K_scan_loc


# Skan eta granularnny
print("\nSkan eta [0.5, 5.0] -> r21:")
eta_fine = np.linspace(0.5, 5.0, 12)
ratio_fine = []

for eta_v in eta_fine:
    r, zeros, _, _ = compute_ratio_for_eta(eta_v, verbose=False)
    ratio_fine.append(r)
    tag = ''
    if not np.isnan(r):
        if abs(r - 207) / 207 < 0.15: tag = '  <-- BLISKI 207!'
        if abs(r - 207) / 207 < 0.05: tag = '  <-- TRAFIONY!'
    z_str = f"({len(zeros)} zer: {[f'{z:.3e}' for z in zeros]})" if zeros else "(brak zer)"
    r_str = f"{r:.2f}" if not np.isnan(r) else "nan"
    print(f"  eta={eta_v:.2f}: r21={r_str}  {z_str}{tag}")

# ============================================================
# KROK 3: Analiza V_mod dla roznych eta
# ============================================================
print()
print("=" * 65)
print("KROK 3: Ksztalt V_mod_new dla roznych eta")
print("=" * 65)
phi_plot = np.linspace(0.8, 5.0, 200)
V1_val   = V1

print(f"\n  phi    | {'V(eta=0)':>10} | {'V(eta=1)':>10} | {'V(eta=3)':>10}")
print("  " + "-"*46)
for phi in [1.0, 1.5, 2.0, 3.0, 4.0]:
    _, V0, _ = make_ode(0.0)
    _, V1e, _ = make_ode(1.0)
    _, V3e, _ = make_ode(3.0)
    print(f"  {phi:5.1f}  | {V0(phi)-V1_val:>10.4f} | {V1e(phi)-V1_val:>10.4f} | {V3e(phi)-V1_val:>10.4f}")

# ============================================================
# WYKRES
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('P19: Modyfikacja V_mod (eta*(phi-1)^4/4) - szukanie r21=207\n'
             f'(alpha={ALPHA}, a_Gam=3.841*K, lam={LAM:.2e})',
             fontsize=11, fontweight='bold')

# Panel 1: g(K) dla roznych eta
ax = axes[0]
for eta in [0.0, 1.0, 2.0, 3.0]:
    if eta not in g_results: continue
    g_arr = g_results[eta]
    mask = np.isfinite(g_arr)
    g_clip = np.clip(g_arr[mask], -10, 15)
    ax.semilogx(K_scan[mask], g_clip, '.-', lw=1.5, ms=4, label=f'eta={eta:.1f}')
ax.axhline(0, color='black', lw=1.2, linestyle='--')
ax.set_xlabel('K'); ax.set_ylabel('g')
ax.set_title('g(K) na linii a_Gam=3.841*K', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_ylim(-8, 12)

# Panel 2: r21 vs eta
ax = axes[1]
eta_valid = [eta_fine[i] for i in range(len(eta_fine)) if not np.isnan(ratio_fine[i])]
r_valid   = [ratio_fine[i] for i in range(len(eta_fine)) if not np.isnan(ratio_fine[i])]
if r_valid:
    ax.semilogy(eta_valid, r_valid, 'bo-', lw=1.5, ms=7)
ax.axhline(207.0,  color='red',    lw=1.5, linestyle='--', label='r21=207')
ax.axhline(28.7,   color='orange', lw=1.5, linestyle=':', label='r21=28.7 (eta=0)')
ax.set_xlabel('eta'); ax.set_ylabel('r21 = K*2/K*1')
ax.set_title('Stosunek mas vs eta', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Panel 3: V_mod dla roznych eta
ax = axes[2]
phi_v = np.linspace(0.9, 4.0, 200)
for eta in [0.0, 0.5, 1.0, 2.0, 3.0]:
    _, V_fn, _ = make_ode(eta)
    V_vals = np.array([V_fn(p) - V1 for p in phi_v])
    V_clip = np.clip(V_vals, -5, 20)
    ax.plot(phi_v, V_clip, lw=1.5, label=f'eta={eta:.1f}')
ax.axhline(0, color='black', lw=0.8, linestyle=':')
ax.axvline(1, color='gray', lw=0.8, linestyle=':')
ax.set_xlabel('phi'); ax.set_ylabel('V_mod(phi) - V1')
ax.set_title('Ksztalt V_mod_new(phi)', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_ylim(-5, 15)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"\nZapisano: {out}")

# ============================================================
# PODSUMOWANIE
# ============================================================
print()
print("=" * 65)
print("PODSUMOWANIE p19:")
print("=" * 65)
print()
print("  V_mod_new = V_mod + eta*(phi-1)^4/4")
print("  - V'(1) = 0 zachowane dla WSZYSTKICH eta   [OK]")
print("  - V''(1) = -gamma NIEZMIENNE (m_eff bez zmian) [OK]")
print("  - V_mod_new ~ (eta-1)/4 * phi^4 dla phi >> 1")
print(f"  - Prog stabilnosci: eta >= 1")
print()

# Wyniki ze skanu
valid_results = [(eta_fine[i], ratio_fine[i]) for i in range(len(eta_fine))
                 if not np.isnan(ratio_fine[i])]
if valid_results:
    print("  Wyniki skanu eta -> r21:")
    for eta_v, r_v in valid_results:
        print(f"    eta={eta_v:.2f}: r21={r_v:.1f}")
    print()

# Sprawdz czy gdzies mamy r21 ~ 207
near_207 = [(e, r) for e, r in valid_results if abs(r-207)/207 < 0.20]
if near_207:
    print(f"  BLISKIE 207: {near_207}")
else:
    max_r = max((r for _, r in valid_results), default=0)
    print(f"  Maks. r21 znalezione: {max_r:.1f}")
    if max_r < 207:
        print(f"  -> Za male. Trzeba wiekszego eta lub innej modyfikacji V_mod.")
    else:
        print(f"  -> MOZLIWE! Dalej sprawdzic dokladny eta*.")
print()
print("  NASTEPNE KROKI (jesli r21<207 dla wszystkich eta):")
print("  - Probowac modyfikacje V: V_6 = phi^3/3 - (1-eps)*phi^4/4 + lam*(phi-1)^6/6")
print("    (redukcja wspolczynnika kwartycznego o eps)")
print("  - Lub: V = phi^3/3 + lam_eff*(phi-1)^6/6 bez czlonu phi^4 (eps=1)")
print("  - Lub: dopasowaLin wejscie C jednoczesnie z eta")
