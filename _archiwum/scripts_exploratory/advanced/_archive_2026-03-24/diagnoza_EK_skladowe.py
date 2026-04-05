"""
Diagnostyka skladowych E_kin/K i E_pot/K osobno
Szukamy: czy drugie minimum V tworzy drugie maximum E/K?
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PHI0 = 1.0; GAMMA = 1.0; M_SP = 1.0; Q = 1.0
Lambda = 4.0 * np.pi

def components_over_K(K, a_gam, alpha, delta):
    r    = np.linspace(a_gam, 30.0, 2000)
    phi  = PHI0 + K * np.exp(-M_SP*r) / r
    phi  = np.maximum(phi, 1e-8)
    dphi = K * np.exp(-M_SP*r) * (-M_SP*r - 1.0) / r**2

    E_kin_std = 4*np.pi * np.trapezoid(0.5*dphi**2 * r**2, r)
    E_kin_nl  = 4*np.pi * np.trapezoid(
        0.5*alpha*dphi**2/(PHI0*phi) * r**2, r)

    V  = GAMMA/3*phi**3 - GAMMA/4*phi**4 + delta/5*phi**5
    V0 = GAMMA/3 - GAMMA/4 + delta/5
    E_pot = 4*np.pi * np.trapezoid((V - V0) * r**2, r)

    psi_core = PHI0 + K*np.exp(-M_SP*a_gam)/a_gam   # psi w centrum

    return E_kin_std/K, E_kin_nl/K, E_pot/K, (E_kin_std+E_kin_nl+E_pot)/K, psi_core

# --- Parametry testowe ---
a_gam = 0.05; alpha = 10.0; delta = 0.20

K_vals = np.linspace(0.005, 15.0, 200)
Ek_std, Ek_nl, Ep, Etot, psi_core = [], [], [], [], []
for K in K_vals:
    a,b,c,d,e = components_over_K(K, a_gam, alpha, delta)
    Ek_std.append(a); Ek_nl.append(b)
    Ep.append(c);     Etot.append(d)
    psi_core.append(e)

Ek_std = np.array(Ek_std); Ek_nl = np.array(Ek_nl)
Ep = np.array(Ep); Etot = np.array(Etot)
psi_core = np.array(psi_core)

# Gdzie psi_core = psi_min (centrum pola trafia w minimum potencjalu)
psi_min = GAMMA / (2*delta)   # psi_min z V'=0: delta < gamma/4
print(f"psi_min (drugie minimum V) = {psi_min:.2f}")
print(f"delta = {delta}, a_gam = {a_gam}, alpha = {alpha}")

# K gdzie psi_core = psi_min
idx_min = np.argmin(np.abs(psi_core - psi_min))
K_at_psi_min = K_vals[idx_min]
print(f"K przy ktorym core osiaga psi_min: K = {K_at_psi_min:.4f}")
print(f"  E_tot/K tam: {Etot[idx_min]:.4f}")
print(f"  Lambda     : {Lambda:.4f}")
print()

# Lokalne ekstrema E_tot/K
dE = np.diff(Etot)
extrema_idx = [i for i in range(len(dE)-1) if dE[i]*dE[i+1] < 0]
print(f"Ekstrema E(K)/K:")
for idx in extrema_idx:
    typ = "MAX" if dE[idx] > 0 else "MIN"
    print(f"  {typ} w K={K_vals[idx]:.4f}, E/K={Etot[idx]:.4f}")

# Liczba przeciec z Lambda
n_cross = sum(
    (Etot[i]-Lambda)*(Etot[i+1]-Lambda) < 0
    for i in range(len(Etot)-1)
    if np.isfinite(Etot[i]) and np.isfinite(Etot[i+1])
)
print(f"Liczba przeciec z Lambda={Lambda:.2f}: {n_cross}")

# --- Wykres ---
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

ax = axes[0]
ax.plot(K_vals, Etot,   'b-',  lw=2.5, label=r'$E_{tot}/K$')
ax.plot(K_vals, Ek_std, 'g--', lw=1.5, label=r'$E_{kin,std}/K$')
ax.plot(K_vals, Ek_nl,  'c--', lw=1.5, label=r'$E_{kin,nl}/K$')
ax.plot(K_vals, Ep,     'r--', lw=1.5, label=r'$E_{pot}/K$')
ax.axhline(Lambda, color='k', ls=':', lw=1.5, label=fr'$\Lambda={Lambda:.1f}$')
ax.set_xlabel(r'$K$', fontsize=12)
ax.set_ylabel(r'$E_{...}(K)/K$', fontsize=12)
ax.set_title(f'Skladowe energii\na_G={a_gam}, alpha={alpha}, delta={delta}', fontsize=10)
ax.set_ylim(-200, 300)
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

ax2 = axes[1]
ax2.plot(K_vals, Etot, 'b-', lw=2)
ax2.axhline(Lambda, color='r', ls='--', lw=1.5, label=fr'$\Lambda$')
ax2.set_xlabel(r'$K$', fontsize=12)
ax2.set_ylabel(r'$E(K)/K$', fontsize=12)
ax2.set_title('E(K)/K — widok globalny', fontsize=10)
mask = np.isfinite(Etot)
ymax = min(200, np.max(Etot[mask]))
ax2.set_ylim(-50, ymax*1.1)
ax2.legend(); ax2.grid(True, alpha=0.3)

ax3 = axes[2]
ax3.plot(K_vals, psi_core, 'm-', lw=2, label=r'$\psi_{core}(K)$')
ax3.axhline(psi_min, color='r', ls='--', lw=1.5,
            label=fr'$\psi_{{min}}={psi_min:.2f}$')
ax3.axhline(1.0, color='gray', ls=':', lw=1)
ax3.set_xlabel(r'$K$', fontsize=12)
ax3.set_ylabel(r'$\psi(r=a_\Gamma)$', fontsize=12)
ax3.set_title('Wartosc pola w centrum vs K', fontsize=10)
ax3.legend(fontsize=10); ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('TGP/TGP_v1/scripts/advanced/diagnoza_EK.png',
            dpi=150, bbox_inches='tight')
print("\nWykres: diagnoza_EK.png")

# --- Rozszerzona diagnostyka: szukaj drugiego max ---
print()
print("Szukanie warunkow dla 2 maksimow E(K)/K:")
print(f"{'a_gam':>7} {'delta':>7} {'alpha':>7} {'n_ext':>6} {'n_cross':>8}")
for a_g in [0.02, 0.03, 0.05]:
    for dlt in [0.18, 0.20, 0.22, 0.24]:
        for alp in [5.0, 10.0, 20.0, 50.0]:
            K_v = np.linspace(0.003, 25.0, 300)
            Et = np.array([components_over_K(K, a_g, alp, dlt)[3] for K in K_v])
            dE = np.diff(Et)
            ext = [i for i in range(len(dE)-1)
                   if np.isfinite(dE[i]) and np.isfinite(dE[i+1])
                   and dE[i]*dE[i+1] < 0 and dE[i] > 0]  # tylko maxima
            nc = sum(
                (Et[i]-Lambda)*(Et[i+1]-Lambda) < 0
                for i in range(len(Et)-1)
                if np.isfinite(Et[i]) and np.isfinite(Et[i+1])
            )
            if nc >= 3 or len(ext) >= 2:
                print(f"{a_g:>7.3f} {dlt:>7.3f} {alp:>7.1f} "
                      f"{len(ext):>6d} {nc:>8d}  <-- {'3 przeciecia!' if nc>=3 else '2 maxima'}")
