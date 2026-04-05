"""
p11_stability_K1.py
===================
Analiza stabilnosci perturbacyjnej solitonu K*_1(ODE) w TGP.

TEZA:
  Soliton K*_1 = 0.010414 jest samospojny (g~0).
  Czy jest stabilny wobec perturbacji pola phi -> phi + delta*phi ?

OPERATOR TGP (z p7):
  H_TGP = -(1/(1+alpha/phi)) d^2/dr^2
        - (2/r) d/dr
        + V''_mod(phi)/(1+alpha/phi)
        - alpha*(phi')^2 / (2*phi^2*(1+alpha/phi))

  Zaburzenie: phi(r,t) = phi_0(r) + eps * u(r) * exp(i*omega*t)
  Stabilnosc: omega^2 > 0 dla wszystkich modow

METODA:
  1. Wygeneruj profil phi_0(r) przy K*_1
  2. Dyskretyzuj H_TGP na siatce logarytmicznej
  3. Oblicz widmo (eigenvalues) H_TGP
  4. Sprawdz czy wszystkie omega^2 > 0 (stabilny) czy sa ujemne (niestabilny)

DODATKOWO:
  Sprawdz warunek Derrricka-TGP:
    dE/dlambda|_{lambda=1} = 0 <=> (-Ek + 3*Ep + DB) = 0
  gdzie DB = granica termin z a_Gam.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ============================================================
# PARAMETRY  (z potwierdzonych wynikow p6/p9)
# ============================================================
ALPHA  = 8.5616
A_GAM  = 0.040
LAM    = 5.501357e-06
GAMMA  = 1.0
K_STAR = 0.010414   # K*_1(ODE) — potwierdzone
PSI_CORE_STAR = 1.2419  # psi_core* — potwierdzone
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

print(f"P11: Analiza stabilnosci solitonu K*_1 = {K_STAR}")
print(f"     psi_core* = {PSI_CORE_STAR},  a_Gam = {A_GAM},  R_max = {R_MAX}")
print()


# ============================================================
# KROK 1: Wygeneruj profil phi_0(r)
# ============================================================
N = 2000
r_eval = A_GAM * (R_MAX / A_GAM) ** np.linspace(0, 1, N)

dphi0 = -K_STAR / A_GAM**2
sol = solve_ivp(
    ode_rhs, [A_GAM, R_MAX], [PSI_CORE_STAR, dphi0],
    method='DOP853', rtol=1e-12, atol=1e-14,
    t_eval=r_eval
)

r    = sol.t
phi0 = sol.y[0]
dphi0_arr = sol.y[1]

print(f"Profil phi_0(r):")
print(f"  phi(a_Gam) = {phi0[0]:.6f}   (cel: {PSI_CORE_STAR:.4f})")
print(f"  phi(R_max) = {phi0[-1]:.6f}  (cel: 1.000)")
print(f"  phi_min    = {phi0.min():.6f}")
print(f"  phi_max    = {phi0.max():.6f}")

# Energia kontrolna
phi_pos = np.maximum(phi0, 1e-10)
Ek = 4*np.pi*np.trapezoid(0.5*dphi0_arr**2*(1+ALPHA/phi_pos)*r**2, r)
Ep = 4*np.pi*np.trapezoid((V_mod(phi_pos)-V1)*r**2, r)
E  = Ek + Ep
print(f"  E = {E:.6f}  (cel: {4*np.pi*K_STAR:.6f})")
print(f"  g = {E/(4*np.pi*K_STAR)-1:.2e}")
print(f"  Ek = {Ek:.5f},  Ep = {Ep:.5f},  Ek/|Ep| = {Ek/abs(Ep):.4f}")
print()


# ============================================================
# KROK 2: Derrick-TGP virial check
# ============================================================
# Skalowanie: phi_lambda(r) = phi_0(r/lambda)
# E(lambda) = lambda * Ek_std + lambda^{-1} * Ek_nl + lambda^3 * Ep + BC(lambda)
# gdzie Ek_std ~ standardowa kinetyczna, Ek_nl ~ czlon alpha/phi

# Standard: dE/dlambda|_1 = Ek_std - Ek_nl/1 + 3*Ep + dBC/dlambda|_1 = 0

# Numerycznie: E(lambda) przez skalowanie siatki
def energy_scaled(lam):
    r_lam  = r * lam   # nowe r: phi_0(r/lam) na starej siatce -> phi_0(r) na lam*r
    phi_s  = phi0      # phi_lambda(r) = phi_0(r/lam) — wartosc taka sama ale r-siatka przeskalowana
    dphi_s = dphi0_arr / lam  # dphi/dr -> dphi/(d(r/lam)*lam) = dphi_0/r * (1/lam)
    phi_s_pos = np.maximum(phi_s, 1e-10)
    Ek_s = 4*np.pi*np.trapezoid(0.5*(dphi_s)**2*(1+ALPHA/phi_s_pos)*r_lam**2, r_lam)
    Ep_s = 4*np.pi*np.trapezoid((V_mod(phi_s_pos)-V1)*r_lam**2, r_lam)
    return Ek_s + Ep_s

lam_arr = np.linspace(0.8, 1.2, 40)
E_lam = np.array([energy_scaled(l) for l in lam_arr])

# dE/dlam numerycznie w lam=1
dE_dlam = np.gradient(E_lam, lam_arr)
dE_at_1 = np.interp(1.0, lam_arr, dE_dlam)

print(f"Derrick-TGP virial:")
print(f"  dE/dlam |_{{lam=1}} = {dE_at_1:.4e}   (powinno byc 0 dla stabilnego)")
print(f"  E(lam=0.9)={energy_scaled(0.9):.5f}, E(lam=1.0)={energy_scaled(1.0):.5f}, E(lam=1.1)={energy_scaled(1.1):.5f}")
print()


# ============================================================
# KROK 3: Operator H_TGP — poprawna dyskretyzacja
# ============================================================
# Metoda: substytucja Schrodingera f(r) = r * u(r)
# -nabla^2 u = -f''/r  (w sferycznej symetrii)
# Operator staje sie 1D Schrodingera: -f''/(w(r)) + U_eff(r)*f = omega^2 * f
# SIATKA: rownoomierna w r (brak problemow z drobnym krokiem log)
#
# Warunki brzegowe: f(a_Gam)=0, f(R_max)=0 (Dirichlet)
#
# Poprawny U_eff z Lagrangianu TGP = +½(1+a/phi)(phi_t^2 - (nabla phi)^2) + V_mod:
#   H_TGP u = omega^2 u
#   omega^2*(1+a/phi_0)*u = nabla^2[(1+a/phi_0)*u]/(1+a/phi_0) - V''(phi_0)*u/(1+a/phi_0) + ...
# W przyblizenniu dominujacym:
#   U_eff(r) = -V''(phi_0)/(1+a/phi_0) + a*(phi_0')^2/(2*phi_0^2*(1+a/phi_0))
# (pochodna czlonu kinetycznego daje drugi wyraz)

N_s = 800
r_s = np.linspace(A_GAM, R_MAX, N_s + 2)  # rownoomierna siatka
dr  = r_s[1] - r_s[0]
r_in = r_s[1:-1]  # wewnetrzne (bez granic)

phi_s  = np.interp(r_in, r, phi0)
dphi_s = np.interp(r_in, r, dphi0_arr)
phi_s  = np.maximum(phi_s, 1e-10)

w_s   = 1.0 + ALPHA / phi_s
V2_s  = d2V_mod(phi_s)
# UWAGA: znak -V'' (TGP Lagrangian ma +V_mod => stabilnosc wymaga -V'')
U_eff = -V2_s / w_s + ALPHA * dphi_s**2 / (2.0 * phi_s**2 * w_s)

V2_inf = d2V_mod(1.0)
w_inf  = 1.0 + ALPHA
U_inf  = -V2_inf / w_inf   # = -(-1)/(1+alpha) = 1/(1+alpha) > 0

print(f"Operator H_TGP (1D Schrodinger, f=r*u, rownoomierna siatka N={N_s}):")
print(f"  dr = {dr:.4e}")
print(f"  U_eff min={U_eff.min():.4e},  max={U_eff.max():.4e}")
print(f"  U_eff(r_min)={U_eff[0]:.4e},  U_eff(r_max)={U_eff[-1]:.4e}")
print(f"  U_eff(phi=1): -V''(1)/(1+alpha) = {U_inf:.6f}  (= m_eff^2 = 1/(1+alpha))")
print()

# Macierz H: -d^2f/dr^2 / w(r) + U_eff(r)*f = omega^2*f
# -f'' / w approx [-f(i+1) + 2*f(i) - f(i-1)] / (dr^2 * w_i) + U_eff_i * f_i
N_in = len(r_in)
diag_main  = 2.0 / (dr**2 * w_s) + U_eff
diag_off   = -1.0 / (dr**2 * w_s[:-1])  # subdiagonal = superdiagonal (symmetric!)
# Note: for proper symmetric matrix, off-diag should use average w or take w at inner point
# Using w at inner point: diag_off[i] couples f[i] and f[i+1]
# diag_off for lower: -1/(dr^2 * w_s[i]) for coupling (i, i-1)
# diag_off for upper: -1/(dr^2 * w_s[i]) for coupling (i, i+1)
# These are asymmetric if w varies — symmetrize with w(i+1/2) = (w_i + w_{i+1})/2
w_half = (w_s[:-1] + w_s[1:]) / 2.0
diag_off_sym = -1.0 / (dr**2 * w_half)

H_sparse = diags([diag_off_sym, diag_main, diag_off_sym], [-1, 0, 1],
                 shape=(N_in, N_in), format='csr')

# Eigenvalues
print(f"Obliczanie eigenvalue H_TGP (macierz {N_in}x{N_in})...")
n_eigs = min(20, N_in - 2)
try:
    omega2, vecs = eigsh(H_sparse, k=n_eigs, which='SA')
    omega2 = np.sort(omega2)
    print(f"  Najnizsze eigenvalue (omega^2):")
    for i, o2 in enumerate(omega2[:10]):
        status = "STABILNY" if o2 > -1e-6 else "NIESTABILNY"
        print(f"    omega^2[{i}] = {o2:.4e}  {status}")
    print()
    n_neg = np.sum(omega2 < -1e-6)
    print(f"  Liczba ujemnych omega^2: {n_neg} / {n_eigs}")
    if n_neg == 0:
        print("  -> SOLITON K*_1 JEST STABILNY (spektralnie)")
    else:
        print(f"  -> SOLITON K*_1 MA {n_neg} MODY NIESTABILNE")
    print()
    # Continuum threshold: omega^2_cont = -V''(1)/(1+alpha) = 1/(1+alpha)
    m_eff_sq = U_inf   # poprawny: -V''(1)/(1+alpha) = 1/(1+alpha)
    n_below_cont = np.sum(omega2 < m_eff_sq - 1e-6)
    print(f"  Prog continuum: omega^2_cont = m_eff^2 = {m_eff_sq:.4f}")
    print(f"  Mody ponizej progu (zwiazane): {n_below_cont}")
    print(f"  Mody powyzej progu (rozpraszanie): {n_eigs - n_below_cont}")

except Exception as e:
    print(f"  Blad eigsh: {e}")
    omega2 = None
    vecs   = None


# ============================================================
# KROK 4: Wykresy
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle(f'p11: Stabilnosc solitonu K*_1={K_STAR}  '
             f'(alpha={ALPHA}, a_Gam={A_GAM})',
             fontsize=11, fontweight='bold')

# Panel 1: Profil phi_0(r)
ax = axes[0, 0]
ax.semilogx(r, phi0, 'b-', lw=2, label='phi_0(r)')
ax.axhline(1.0, color='gray', lw=0.8, linestyle=':')
ax.axvline(A_GAM, color='orange', lw=1, linestyle='--', label=f'a_Gam={A_GAM}')
ax.set_xlabel('r')
ax.set_ylabel('phi(r)')
ax.set_title('Profil solitonu K*_1')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(0.9, 1.4)

# Panel 2: E(lambda) — virial Derricka
ax = axes[0, 1]
ax.plot(lam_arr, E_lam, 'b-', lw=2)
ax.axvline(1.0, color='red', lw=1, linestyle='--', label='lambda=1')
ax.axhline(E, color='gray', lw=0.7, linestyle=':')
ax.set_xlabel('lambda (skala r)')
ax.set_ylabel('E(lambda)')
ax.set_title(f'Derrick: dE/dlam|_1={dE_at_1:.3e}')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 3: U_eff(r) — potencjal efektywny
ax = axes[1, 0]
ax.semilogx(r_in, U_eff, 'b-', lw=1.5, label='U_eff(r)')
ax.axhline(0, color='black', lw=0.7, linestyle=':')
ax.axhline(U_inf, color='red', lw=1, linestyle='--',
           label=f'U_inf={U_inf:.4f}')
ax.set_xlabel('r')
ax.set_ylabel('U_eff(r)')
ax.set_title('Potencjal efektywny H_TGP')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(-2, 2)

# Panel 4: Spectrum omega^2
ax = axes[1, 1]
if omega2 is not None:
    ax.scatter(range(len(omega2)), omega2, s=40, color='blue')
    ax.axhline(0, color='red', lw=1, linestyle='--', label='omega^2=0')
    ax.axhline(m_eff_sq, color='green', lw=1, linestyle=':',
               label=f'continuum: m_eff^2={m_eff_sq:.3f}')
    ax.set_xlabel('Nr modu')
    ax.set_ylabel('omega^2')
    ax.set_title('Widmo H_TGP (eigenvalues)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(min(-0.2, omega2[:5].min()-0.05), omega2[5]+0.1)
else:
    ax.text(0.5, 0.5, 'Blad eigenvalue', ha='center', va='center',
            transform=ax.transAxes)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")
print()

# ============================================================
# KROK 5: Mody wlasne — wizualizacja najnizszych
# ============================================================
if omega2 is not None and vecs is not None:
    print("Najnizsze mody perturbacyjne u_i(r):")
    fig2, axes2 = plt.subplots(1, min(4, n_eigs), figsize=(14, 4))
    if not hasattr(axes2, '__len__'):
        axes2 = [axes2]
    for idx, ax2 in enumerate(axes2):
        u = vecs[:, idx]
        u = u / (np.max(np.abs(u)) + 1e-30)
        ax2.plot(r_in, u, 'b-', lw=1.5)
        ax2.axhline(0, color='black', lw=0.5)
        ax2.set_title(f'u_{idx}: omega^2={omega2[idx]:.3e}', fontsize=9)
        ax2.set_xscale('log')
        ax2.set_xlabel('r')
        ax2.grid(True, alpha=0.3)
        print(f"  u_{idx}: omega^2={omega2[idx]:.4e}, wezly={np.sum(np.diff(np.sign(u)) != 0)}")
    plt.tight_layout()
    out2 = __file__.replace('.py', '_modes.png')
    plt.savefig(out2, dpi=110, bbox_inches='tight')
    plt.close()
    print(f"  Zapisano: {out2}")

print()
print("PODSUMOWANIE p11:")
print(f"  K*_1(ODE) = {K_STAR}")
print(f"  E[phi*] = {E:.5f} = 4*pi*K*_1 = {4*np.pi*K_STAR:.5f}  (blad: {abs(E-4*np.pi*K_STAR):.2e})")
print(f"  dE/dlam|_1 = {dE_at_1:.4e}  (Derrick-TGP virial)")
if omega2 is not None:
    print(f"  omega^2_min = {omega2[0]:.4e}")
    if omega2[0] > -1e-5:
        print("  STATUS: SOLITON STABILNY spektralnie")
    else:
        print(f"  STATUS: SOLITON NIESTABILNY - {np.sum(omega2<-1e-5)} ujemnych modow")
