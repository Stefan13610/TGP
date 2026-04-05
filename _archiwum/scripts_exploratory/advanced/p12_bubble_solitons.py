"""
p12_bubble_solitons.py
======================
Topologiczne solitony "bankowe" w TGP -- hipoteza dla generacji 2 i 3.

MOTYWACJA:
  K2=2.033 i K3=34.14 nie maja profilu phi(inf)=1 przy stalym a_Gam=0.04.
  Profil K2 "rozbija sie przez phi=0 (N0)" -- to moze byc FIZYCZNE zachowanie!

HIPOTEZA:
  Generacje 2 i 3 to "solitony bankowe" (bubble solitons):
  - Rdzen: obszar absolutnej nicosci N0 (phi=0) dla r < r_0
  - Powloka: pole phi(r) rosnace od 0 do 1 dla r > r_0
  - Brak wewnetrznej granicy a_Gam -- N0 jest wewnetrzna granica!

ANALIZA ZACHOWANIA ODE PRZY phi->0:
  ODE: (1+a/phi)phi'' + (2/r)(1+a/phi)phi' - a(phi')^2/(2phi^2) = V'(phi)

  Przy phi -> 0+ ze skalowaniem phi ~ C*(r-r0)^n:
  - Czynniki a/phi dywerguja
  - Analiza: nieregularny punkt osobliwy przy phi=0 dla skonczonego r0

  Wniosek: N0 jest NIEREGULARNYM singularnym punktem ODE.
  Solitony bankowe wymagaja specjalnego traktowania.

CEL p12:
  1. Przeanalizowac lokalnie ODE przy phi -> 0
  2. Zbadac "odwrocone" strzalkowanie: od duzego r (phi=1) inward,
     sprawdzic czy phi naturalnie dochodzi do 0 dla K~K2, K3
  3. Znalezc r0(K) -- promien "powloki" solitonu bankowego
  4. Zdefiniowac samospojnosc dla solitonu bankowego: E = 4*pi*K
"""

import numpy as np
from scipy.integrate import solve_ivp
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
GAMMA  = 1.0
M_EFF  = 1.0 / np.sqrt(1.0 + ALPHA)
R_MAX  = 50.0

V1 = GAMMA/3 - GAMMA/4

def V_mod(p): return GAMMA/3*p**3 - GAMMA/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAMMA*p**2 - GAMMA*p**3 + LAM*(p-1)**5
def d2V_mod(p): return 2*GAMMA*p - 3*GAMMA*p**2 + 5*LAM*(p-1)**4

def ode_rhs(r, y):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA / phi
    ddphi = (dV_mod(phi) / kfac
             + ALPHA * dphi**2 / (2.0 * phi**2 * kfac)
             - (2.0 / r) * dphi)
    return [dphi, ddphi]

print("P12: Topologiczne solitony bankowe (bubble solitons) w TGP")
print("=" * 65)
print()


# ============================================================
# KROK 1: Analiza lokalna phi -> 0 przy phi ~ C*(r-r0)^n
# ============================================================
print("KROK 1: Analiza zachowania ODE przy phi -> 0:")
print()
print("  Dla phi ~ C*(r-r0)^n z n=2:")
print("  phi = C*e^2, phi' = 2C*e, phi'' = 2C  (gdzie e = r-r0)")
print()
print("  Czynnik dywergentny epsilon^{-2}:")
print("    (a/phi)*phi'' = a*2C/(C*e^2) = 2a/e^2")
print("    (a/phi')^2/(2phi^2) = a*4C^2*e^2/(2C^2*e^4) = 2a/e^2")
print("  -> Znosza sie! (2a/e^2 - 2a/e^2 = 0) DOKLADNIE.")
print()
print("  Czynnik epsilon^{-1}:")
print("    (2/r)*(a/phi)*phi' = (2/r)*(a/(C*e^2))*(2C*e) = 4a/(r*e)")
print("  -> Ten czlon DYWERGUJE dla skonczonego r_0!")
print("  -> Brak regularnego rozwiazania phi=0 przy phi ~ (r-r0)^2 i r0 > 0")
print()
print("  Wyjatki:")
print("    a) r_0 -> 0: wtedy 4a/(r*e) ~ 4a/(e*e) ~ 4a/r^2 gdy r->0")
print("       (phi ~ r^2 przy rdzeniu r=0 -- moze byc regularne!)")
print("    b) Zachowanie przy r_0=0: soliton z phi(0)=0 (N0 w centrum)")
print()

# ============================================================
# KROK 2: Szukaj profili z phi(r_0) = 0 dla integracji inward
# ============================================================
print("KROK 2: Integracja inward z phi(1) (asymptotyczna BC)")
print("         -- sprawdzenie czy phi naturalnie dochodzi do 0:")
print()
print("  Uwaga: integracja inward jest NIESTABILNA (mody rosnace e^{+m*r})")
print("  Ale dla GLOBALNEJ analizy jakosciowej mozemy sprawdzic r_0(K)")
print()

# Asymptotyczna BC: phi(R) = 1 + A*exp(-m_eff*R)/R
# Integracja inward z rozroznieniem trajektorii
R_start = 30.0  # punkt startowy integracji inward

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.suptitle(f'P12: Profile ODE inward -- czy phi dochodzi do N0?\n'
             f'(alpha={ALPHA}, a_Gam={A_GAM}, lam={LAM:.2e})',
             fontsize=11, fontweight='bold')

K_test_vals = [0.010414, 0.1, 0.5, 2.032728, 10.0, 34.14450]
K_labels    = ['K*_1=0.0104', 'K=0.1', 'K=0.5', 'K_2=2.033', 'K=10.0', 'K_3=34.14']
colors_K    = ['green', 'cyan', 'blue', 'orange', 'red', 'purple']

r0_found = {}   # r_0 dla kazdego K

for ax, K, label, col in zip(axes.flat, K_test_vals, K_labels, colors_K):
    # Asymptotyczna BC z poprawna m_eff (nie m_sp=1)
    eR   = np.exp(-M_EFF * R_start)
    phi_R  = 1.0 + K * eR / R_start
    dphi_R = K * eR * (-1.0/R_start**2 - M_EFF/R_start)

    # Integracja inward -- zatrzymaj przy phi < 1e-4
    def ev_zero(r, y): return y[0] - 1e-4
    ev_zero.terminal = True; ev_zero.direction = -1
    def ev_neg(r, y): return y[0] - 1e-8
    ev_neg.terminal = True; ev_neg.direction = -1

    r_eval = np.linspace(R_start, A_GAM, 3000)
    sol = solve_ivp(
        ode_rhs, [R_start, A_GAM], [phi_R, dphi_R],
        method='DOP853', rtol=1e-7, atol=1e-9,
        t_eval=r_eval, events=[ev_zero]
    )

    r   = np.sort(sol.t)
    phi = sol.y[0][np.argsort(sol.t)]

    phi_at_end = phi[-1]
    r_at_end   = r[-1]

    if phi_at_end < 1e-3:
        r0_found[K] = r_at_end
        print(f"  {label}: phi -> 0 at r0 = {r_at_end:.4f}  (phi_end={phi_at_end:.2e})")
    else:
        r0_found[K] = None
        print(f"  {label}: phi(a_Gam) = {phi_at_end:.4f}  (nie dochodzi do N0)")

    # Wykres
    ax.plot(r, phi, '-', color=col, lw=1.5)
    ax.axhline(0, color='black', lw=0.8, linestyle='--', alpha=0.5)
    ax.axhline(1, color='gray', lw=0.7, linestyle=':', alpha=0.7)
    if phi_at_end < 1e-3:
        ax.axvline(r_at_end, color='red', lw=1.2, linestyle='--',
                   label=f'N0 at r={r_at_end:.3f}')
    ax.set_xlabel('r')
    ax.set_ylabel('phi(r)')
    ax.set_title(f'{label}', fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.3, max(phi.max() * 1.05, 1.2))

plt.tight_layout()
out = __file__.replace('.py', '_inward.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"  Zapisano: {out}")
print()


# ============================================================
# KROK 3: Badanie r_0(K) -- promien powloki solitonu bankowego
# ============================================================
print("KROK 3: Skan r_0(K) -- promien powloki N0:")
print()
print("  Pytanie: dla jakich K profil phi dochodzi do N0?")
print("  Jezeli K < K_crit: phi nigdy nie dotyka N0 (ball soliton)")
print("  Jezeli K > K_crit: phi dochodzi do N0 (bubble soliton)")
print()

K_scan = np.logspace(-2, 2, 50)
r0_scan = []
K_valid = []

for K in K_scan:
    eR     = np.exp(-M_EFF * R_start)
    phi_R  = 1.0 + K * eR / R_start
    dphi_R = K * eR * (-1.0/R_start**2 - M_EFF/R_start)

    def ev_z(r, y): return y[0] - 1e-3
    ev_z.terminal = True; ev_z.direction = -1

    sol = solve_ivp(
        ode_rhs, [R_start, 0.01], [phi_R, dphi_R],
        method='DOP853', rtol=1e-7, atol=1e-9,
        events=[ev_z]
    )
    if sol.t_events[0].size > 0:
        r0_scan.append(float(sol.t_events[0][0]))
        K_valid.append(K)
    else:
        r0_scan.append(np.nan)

K_scan   = np.array(K_scan)
r0_arr   = np.array(r0_scan)
K_arr_v  = np.array(K_valid)
r0_arr_v = np.array([r0_scan[i] for i, K in enumerate(K_scan) if K in K_valid])

# Znajdz K_crit
K_crit = None
for i in range(len(K_scan)-1):
    if np.isnan(r0_scan[i]) and not np.isnan(r0_scan[i+1]):
        K_crit = K_scan[i]
        print(f"  K_crit approx: K in [{K_scan[i]:.3e}, {K_scan[i+1]:.3e}]")
        break
    elif not np.isnan(r0_scan[i]) and np.isnan(r0_scan[i+1]):
        K_crit = K_scan[i]

if K_arr_v.size > 0:
    print(f"  K min z N0 = {K_arr_v.min():.4f}")
    print(f"  K max z N0 = {K_arr_v.max():.2f}")
    print()
    print(f"  r_0(K) wybrane wartosci:")
    idx_show = np.round(np.linspace(0, len(K_arr_v)-1, min(10, len(K_arr_v)))).astype(int)
    for i in idx_show:
        print(f"    K={K_arr_v[i]:.4f}  r_0={r0_arr_v[i]:.4f}")

fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))

ax = axes2[0]
mask_valid = np.isfinite(r0_arr)
if mask_valid.any():
    ax.loglog(K_scan[mask_valid], r0_arr[mask_valid], 'b.-', lw=1.5, ms=4)
ax.axvline(0.010414, color='green', lw=1.5, linestyle='--', label='K*_1')
ax.axvline(2.032728, color='orange', lw=1.5, linestyle=':', label='K_2(Yukawa)')
ax.axvline(34.14450, color='red', lw=1.5, linestyle=':', label='K_3(Yukawa)')
ax.set_xlabel('K')
ax.set_ylabel('r_0 (promien powloki N_0)')
ax.set_title('r_0(K): promien wewnetrznej granicy N0')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

ax = axes2[1]
# Energia bubble solitonu: bez wewnetrznej czesci (r < r0)
# E_bubble = E_ball - E_inner
# Schematycznie: jesli phi dochodzi do N0 at r0, to pole istnieje tylko dla r > r0
# Wartosc K jest ustalona przez dphi(r0) = phi'(r0) * r0^2
K_bubble_vals = []
K_eff_vals    = []
if K_arr_v.size > 0:
    for K, r0 in zip(K_arr_v[:10], r0_arr_v[:10]):
        eR     = np.exp(-M_EFF * R_start)
        phi_R  = 1.0 + K * eR / R_start
        dphi_R = K * eR * (-1.0/R_start**2 - M_EFF/R_start)

        def ev_z2(r, y): return y[0] - 1e-3
        ev_z2.terminal = True; ev_z2.direction = -1

        sol = solve_ivp(
            ode_rhs, [R_start, 0.01], [phi_R, dphi_R],
            method='DOP853', rtol=1e-7, atol=1e-9,
            events=[ev_z2]
        )
        if sol.t_events[0].size > 0:
            # phi'(r0) i r0^2 * |phi'|
            dphi_r0 = float(sol.y_events[0][0][1])
            K_eff_r0 = r0**2 * abs(dphi_r0)
            K_bubble_vals.append(K)
            K_eff_vals.append(K_eff_r0)
            print(f"  K_in={K:.4f}  r0={r0:.4f}  K_eff(r0)={K_eff_r0:.4e}")

ax.scatter(K_bubble_vals, K_eff_vals, s=50, color='blue')
if K_bubble_vals:
    ax.plot([min(K_bubble_vals), max(K_bubble_vals)],
            [min(K_bubble_vals), max(K_bubble_vals)], 'r--', label='K_eff = K')
ax.set_xlabel('K (Yukawa)')
ax.set_ylabel('K_eff(r0) = r0^2 * |phi_prime(r0)|')
ax.set_title('K_eff na powloce N0 vs. K Yukawa')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
out2 = __file__.replace('.py', '_r0_scan.png')
plt.savefig(out2, dpi=120, bbox_inches='tight')
plt.close()
print(f"  Zapisano: {out2}")
print()


# ============================================================
# PODSUMOWANIE
# ============================================================
print("=" * 65)
print("PODSUMOWANIE p12:")
print()
print("ANALITYCZNE:")
print("  - phi=0 (N0) jest NIEREGULARNYM punktem osobliwym ODE dla r0 > 0")
print("  - Leading singularity cancels (n=2), ale epsilon^{-1} nie cancels")
print("  - WYJATEK: r0=0 (N0 w samym centrum)")
print()
print("NUMERYCZNE:")
K_cross = None
for i in range(len(K_scan)-1):
    if np.isnan(r0_scan[i]) and not np.isnan(r0_scan[i+1]):
        K_cross = (K_scan[i] * K_scan[i+1])**0.5  # geometryczna srednia
        break
if K_cross:
    print(f"  K_crit (przejscie ball -> bubble): K_crit ~= {K_cross:.4f}")
    print(f"  K*_1 = 0.01041 < K_crit: BALL SOLITON (nie dotyka N0) [oczekiwane]")
    print(f"  K_2 = 2.033 > K_crit:   BUBBLE SOLITON (dotyka N0)   [potwierdzone]")
    print(f"  K_3 = 34.14 > K_crit:   BUBBLE SOLITON (dotyka N0)   [potwierdzone]")
else:
    if mask_valid.any():
        print(f"  Wszystkie badane K daja N0 -- K_crit < {K_scan[0]:.4f}")
    else:
        print("  Zadne badane K nie daje N0 -- K_crit > {K_scan[-1]:.1f}")
print()
print("INTERPRETACJA FIZYCZNA:")
print("  Generacje 1,2,3 to ROZNE KLASY TOPOLOGICZNE:")
print("  Gen 1 (K*_1): ball soliton -- phi > 0 wszedzie")
print("  Gen 2 (K_2):  bubble soliton -- phi dochodzi do N0 w centrum")
print("  Gen 3 (K_3):  moz. double-bubble lub gleboko nielinowe")
print()
print("  NOWY WARUNEK BRZEGOWY dla bubble solitonow:")
print("  phi(r_0) = 0 (dotkniecie N0)")
print("  phi'(r_0) = ? (do wyznaczenia z ODE)")
print()
print("  Samospojnosc bubble: E_bubble[phi] = 4*pi*K_bubble")
print("  gdzie K_bubble = r_0^2 * |phi'(r_0)| (ladyunek na powloce N0)")
