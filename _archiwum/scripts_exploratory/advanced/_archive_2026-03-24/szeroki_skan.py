"""
Szeroki skan K do K_max=200: szukanie trzeciego przeciecia
przy bardzo duzych K (tam gdzie V5 = delta*psi^5 dominuje).

Fizyczny argument:
  E_kin ~ K^2/a_G  (maleje wzg K)
  E_pot ~ delta*K^5/(5*a_G^5) (rosnie dla duzych K)
=> E(K)/K ~ delta*K^4/a_G^5 ROSNIE dla K >> a_G
=> TRZECIE przeciecie z Lambda istnieje dla duzych K!
"""
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PHI0 = 1.0; GAMMA = 1.0; M_SP = 1.0; Q = 1.0
LAM = 5.514e-7  # V_mod: + LAM*(psi-1)^6/6 stabilizacja (efekt < 0.03% na r21)
Lambda = 4.0 * np.pi

def E_over_K(K, a_gam, alpha, delta, r_max=None, N=2000):
    if r_max is None:
        r_max = max(30.0, 10.0/M_SP)
    r    = np.linspace(a_gam, r_max, N)
    phi  = PHI0 + K * np.exp(-M_SP*r) / r
    phi  = np.maximum(phi, 1e-8)
    dphi = K * np.exp(-M_SP*r) * (-M_SP*r - 1.0) / r**2
    Ek   = 4*np.pi * np.trapezoid(
        0.5*dphi**2*(1+alpha/(PHI0*phi)) * r**2, r)
    V    = GAMMA/3*phi**3 - GAMMA/4*phi**4 + delta/5*phi**5 + LAM/6*(phi-1.0)**6
    V0   = GAMMA/3 - GAMMA/4 + delta/5
    Ep   = 4*np.pi * np.trapezoid((V-V0)*r**2, r)
    return (Ek + Ep) / K

def scan_wide(a_gam, alpha, delta, K_max=150, N_K=400):
    # Gestsza siatka przy malych K, rozrzedzona przy duzych
    K_low  = np.linspace(0.003, 5.0,   int(N_K*0.5))
    K_high = np.linspace(5.0,   K_max, int(N_K*0.5))
    K_vals = np.unique(np.concatenate([K_low, K_high]))
    EK     = np.array([E_over_K(K, a_gam, alpha, delta) for K in K_vals])
    f      = EK - Lambda
    crossings = []
    for i in range(len(f)-1):
        if np.isfinite(f[i]) and np.isfinite(f[i+1]) and f[i]*f[i+1] < 0:
            try:
                Kc = brentq(
                    lambda K: E_over_K(K, a_gam, alpha, delta) - Lambda,
                    K_vals[i], K_vals[i+1], xtol=1e-6)
                crossings.append(Kc)
            except Exception:
                pass
    return K_vals, EK, crossings

# ====================================================
# ANALITYCZNA PROGNOZA K3 z dominacji V5
# ====================================================
# Dla duzych K: E_pot ~ delta*K^5 * C_pot
# gdzie C_pot = 4*pi * integral (1/r^5)*r^2 dr od a_gam ~ 4*pi/(2*a_gam^2)
# E_pot/K ~ delta*K^4 * 4*pi/(2*a_gam^2) = Lambda
# K3 ~ (Lambda * 2*a_gam^2 / (4*pi*delta))^(1/4) / a_gam^(...)

def K3_estimate(a_gam, delta):
    """Bardzo gruba estymacja K3 z E_pot/K ~ Lambda."""
    # E_pot ~ 4pi * delta/5 * K^5 * integral r^-3 dr od a_gam
    # integral 1/r^3 * r^2 dr = integral 1/r dr od a_gam -> log-divergent
    # Uzywamy 4pi * delta/5 * K^5 / (2*a_gam^2) ~ Lambda * K
    # K^4 ~ Lambda * 5 * 2*a_gam^2 / (4*pi*delta)
    # K3 ~ (10*Lambda*a_gam^2 / (4*pi*delta))^(1/4)
    return (10*Lambda*a_gam**2 / (4*np.pi*delta))**0.25

print("="*65)
print("Analityczna estymacja K3 (trzecia generacja)")
print("="*65)
print(f"{'a_gam':>7} {'delta':>7} {'K3_est':>10} {'m3/m1_est':>12}")
print("-"*65)

a_gam_test = 0.05; alpha_test = 10.0
for delta in [0.10, 0.15, 0.20, 0.22, 0.24]:
    K3e = K3_estimate(a_gam_test, delta)
    K1e = 0.02   # z poprzednich wynikow (bardzo gruba aproksymacja)
    print(f"{a_gam_test:>7.3f} {delta:>7.3f} {K3e:>10.2f} {K3e/K1e:>12.1f}")

# ====================================================
# PELNY SKAN SZEROKI dla wybranych parametrow
# ====================================================
print()
print("="*65)
print("Pelny skan szeroki (K do 150)")
print("="*65)
print(f"{'a_gam':>6} {'alpha':>6} {'delta':>7} {'n':>4} "
      f"{'K1':>8} {'K2':>8} {'K3':>8} {'m2/m1':>8} {'m3/m1':>8}")
print("-"*65)

three_cases = []
for a_gam in [0.03, 0.05, 0.08]:
    for alpha in [5.0, 10.0, 20.0]:
        for delta in [0.15, 0.18, 0.20, 0.22]:
            K_v, EK, crossings = scan_wide(a_gam, alpha, delta, K_max=120, N_K=350)
            n = len(crossings)
            masses = [Lambda*K for K in crossings]
            K1 = crossings[0] if n>=1 else None
            K2 = crossings[1] if n>=2 else None
            K3 = crossings[2] if n>=3 else None
            m1 = masses[0] if n>=1 else None
            r21 = masses[1]/masses[0] if n>=2 else None
            r31 = masses[2]/masses[0] if n>=3 else None

            K1s  = f"{K1:.3f}"  if K1  else "-"
            K2s  = f"{K2:.3f}"  if K2  else "-"
            K3s  = f"{K3:.1f}"  if K3  else "-"
            r21s = f"{r21:.1f}" if r21 else "-"
            r31s = f"{r31:.1f}" if r31 else "-"

            if n >= 2:
                flag = "  <-- 3 GEN!" if n>=3 else ""
                print(f"{a_gam:>6.3f} {alpha:>6.1f} {delta:>7.3f} {n:>4d} "
                      f"{K1s:>8} {K2s:>8} {K3s:>8} {r21s:>8} {r31s:>8}{flag}")
                if n >= 3:
                    three_cases.append((a_gam, alpha, delta, crossings, masses))

# ====================================================
# WYKRES dla najlepszego przypadku (max crossings)
# ====================================================
if three_cases:
    best = three_cases[0]
    a_gam, alpha, delta, crossings, masses = best
else:
    # Brak 3, wez najlepszy z 2
    a_gam, alpha, delta = 0.05, 10.0, 0.20
    K_v, EK, crossings = scan_wide(a_gam, alpha, delta, K_max=120, N_K=350)
    masses = [Lambda*K for K in crossings]

K_v, EK, crossings2 = scan_wide(a_gam, alpha, delta, K_max=120, N_K=350)
masses2 = [Lambda*K for K in crossings2]

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
mask = np.isfinite(EK) & (EK > -500) & (EK < 500)
ax.plot(K_v[mask], EK[mask], 'b-', lw=2, label=r'$E(K)/K$')
ax.axhline(Lambda, color='r', ls='--', lw=1.5, label=fr'$\Lambda={Lambda:.2f}$')
colors = ['green', 'orange', 'purple']
gen_names = ['Gen 1', 'Gen 2', 'Gen 3']
for i, (Kc, mc) in enumerate(zip(crossings2[:3], masses2[:3])):
    ax.axvline(Kc, color=colors[i], ls=':', lw=1.5)
    ax.plot(Kc, Lambda, 'o', color=colors[i], ms=11,
            label=f'{gen_names[i]}: K={Kc:.3f}, m={mc:.2f}')
ax.set_xlabel(r'$K$', fontsize=12)
ax.set_ylabel(r'$E(K)/K$', fontsize=12)
ax.set_title(f'Szeroki skan: a_G={a_gam}, alpha={alpha}, delta={delta}', fontsize=10)
ax.set_ylim(-100, 200)
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

ax2 = axes[1]
# Zoom na male K
K_zoom = K_v[K_v < 5.0]
EK_zoom = EK[K_v < 5.0]
ax2.plot(K_zoom, EK_zoom, 'b-', lw=2)
ax2.axhline(Lambda, color='r', ls='--', lw=1.5)
for Kc, mc in zip(crossings2, masses2):
    if Kc < 5.0:
        ax2.axvline(Kc, color='g', ls=':', lw=1.5)
        ax2.plot(Kc, Lambda, 'go', ms=10)
ax2.set_xlabel(r'$K$', fontsize=12)
ax2.set_ylabel(r'$E(K)/K$', fontsize=12)
ax2.set_title('Zoom K < 5', fontsize=10)
ax2.set_ylim(-50, 150)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('TGP/TGP_v1/scripts/advanced/szeroki_skan.png',
            dpi=150, bbox_inches='tight')
print()
print(f"Wykres: szeroki_skan.png")
if masses2:
    print(f"Przeciecia: {[f'K={K:.3f} m={m:.2f}' for K,m in zip(crossings2,masses2)]}")
    if len(masses2) >= 2:
        print(f"m2/m1 = {masses2[1]/masses2[0]:.1f}  (cel: 207)")
    if len(masses2) >= 3:
        print(f"m3/m1 = {masses2[2]/masses2[0]:.1f}  (cel: 3477)")
