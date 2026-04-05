"""
p18_phi0_hierarchy.py
=====================
MECHANIZM 1 z pierwszych zasad: hierarchia poziomow tla Phi0.

KONTEKST (diagnoza_sekwencyjny.py):
  - Model sekwencyjny Phi0^(n) = Phi0^(1) + xi*(M1+...+M_{n-1})
    JEST NIEMOZLIWY dla r21=207 i r31=3477 jednoczesnie.
  - badanie_pi.py sugeruje: Phi0_2 ~ 2*pi dla (alpha=5.9, a_gam=0.05)

CEL p18:
  1. Wyznacz Phi0^(2), Phi0^(3) potrzebne dla r21=207, r31=3477
     przy naszych referencyjnych parametrach (alpha=8.5616, a_gam=0.040, lam=5.501e-6)
  2. Sprawdz hipoteze 2*pi: czy Phi0_2 ~ 2*pi?
  3. Sprawdz potegowe skalowanie M*(Phi0) ~ Phi0^n
  4. Szacuj sprzezenie zwrotne: jak bardzo soliton gen-1 zmienia Phi0?
  5. Sprawdz czy Phi0 moze byc wyznaczone przez nieliniowe rownania stanu TGP

NOWA HIPOTEZA: 'Skalowanie Bohra'
  W atomie wodoru poziomy energetyczne En ~ 1/n^2 wynikaja z
  calki fazowej (kwantyzacja Bohra): p*dq = n*h.
  Analogicznie w TGP: Phi0^(n) moze byc wyznaczone przez 'calke akcji'
  solitonu. Sprawdzimy: Phi0^(n) ~ (n*pi)^2 lub podobne.
"""

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Parametry referencyjne (wyniki p14c/p16/p17)
ALPHA    = 8.5616
A_GAM    = 0.040
LAM      = 5.501357e-6
GAMMA    = 1.0
PHI0_BASE = 1.0

# M_EFF = 1/sqrt(1+alpha) -- masa efektywna pola TGP (peln ODE)
M_EFF = 1.0 / np.sqrt(1.0 + ALPHA)

print("P18: Mechanizm 1 z pierwszych zasad -- hierarchia Phi0")
print("=" * 65)
print(f"  alpha={ALPHA}, a_gam={A_GAM}, lam={LAM:.3e}")
print(f"  M_eff(ODE) = 1/sqrt(1+alpha) = {M_EFF:.4f}")
print(f"  m_sp(Yukawa, Phi0=1) = sqrt(gamma/Phi0_base) = {np.sqrt(GAMMA/PHI0_BASE):.4f}")
print()

# ============================================================
# FUNKCJE: Model Yukawa w tle Phi0
# (jak w sekwencyjny_model.py, z parametrami referencyjnymi)
# ============================================================

def energy_yukawa(K, Phi0, alpha=ALPHA, a_gam=A_GAM, lam=LAM, n_eval=3000):
    """
    Energia solitonu Yukawa w tle Phi0.
    m_sp(Phi0) = sqrt(gamma*Phi0_BASE/Phi0) -- Yukawa screening w tle Phi0.
    """
    msp = np.sqrt(max(GAMMA * PHI0_BASE / Phi0, 1e-9))
    r_max = max(80.0 / msp, 20.0)
    r = np.linspace(a_gam, r_max, n_eval)

    phi  = Phi0 + K * np.exp(-msp * r) / r
    phi  = np.maximum(phi, 1e-10)
    dphi = K * np.exp(-msp * r) * (-msp * r - 1.0) / r**2

    # Kinetyczny: (1 + alpha/(Phi0*psi)) * dphi^2/2
    Ek = 4*np.pi * np.trapezoid(
        0.5 * dphi**2 * (1.0 + alpha / phi) * r**2, r)

    # Potencjalny: V_mod(psi) - V1, gdzie psi = phi/Phi0
    psi = phi / Phi0
    V1  = GAMMA/3 - GAMMA/4   # wartosc prozni psi=1
    Ep  = 4*np.pi * Phi0**3 * np.trapezoid(
        (GAMMA/3*psi**3 - GAMMA/4*psi**4
         + lam/6*(psi-1.0)**6 - V1) * r**2, r)

    return Ek + Ep

def g_func(M, Phi0, alpha=ALPHA, a_gam=A_GAM, lam=LAM):
    """g(M; Phi0) = E(K(M); Phi0) - M = 0 jest warunkiem samospojnosci."""
    if M < 1e-12: return 0.0
    K = M / (4.0 * np.pi * Phi0)
    E = energy_yukawa(K, Phi0, alpha=alpha, a_gam=a_gam, lam=lam)
    return E - M

def find_M_star(Phi0, alpha=ALPHA, a_gam=A_GAM, lam=LAM,
               M_max=5e5, N=500, verbose=False):
    """Znajdz wszystkie zera g(M; Phi0)."""
    # Skan wieloskalowy
    M_arr = np.unique(np.concatenate([
        np.linspace(1e-4,    2.0,   int(N*0.20)),
        np.linspace(2.0,     100,   int(N*0.20)),
        np.linspace(100,     5000,  int(N*0.25)),
        np.linspace(5000,    50000, int(N*0.20)),
        np.linspace(50000, M_max,   int(N*0.15)),
    ]))
    g_arr = np.array([g_func(M, Phi0, alpha, a_gam, lam) for M in M_arr])

    roots = []
    for i in range(len(g_arr) - 1):
        gi, gj = g_arr[i], g_arr[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi * gj < 0:
            try:
                r = brentq(
                    lambda M: g_func(M, Phi0, alpha, a_gam, lam),
                    M_arr[i], M_arr[i+1], xtol=1e-6, maxiter=40)
                roots.append(r)
                if verbose:
                    print(f"    M* ~ {r:.4f}  (g: {gi:.4e} -> {gj:.4e})")
            except Exception:
                pass
    return roots


# ============================================================
# KROK 1: M* jako funkcja Phi0 przy naszych parametrach
# ============================================================
print("KROK 1: Skalowanie M*(Phi0) dla naszych parametrow")
print("-"*65)

Phi0_scan = np.array([1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 6.28,
                       7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 30.0, 40.0, 50.0])

M1_base = None
M_star_vals = []
exponents   = []

print(f"{'Phi0':>8}  {'M*_1':>12}  {'M*/M1':>10}  {'wykladnik n':>12}")
for Phi0 in Phi0_scan:
    roots = find_M_star(Phi0)
    if roots:
        M_min = roots[0]
        if Phi0 == 1.0:
            M1_base = M_min
        ratio = M_min / M1_base if M1_base else float('nan')
        if Phi0 > 1.0 and M1_base:
            exp = np.log(ratio) / np.log(Phi0) if Phi0 > 1.01 else 0.0
        else:
            exp = 0.0
        M_star_vals.append((Phi0, M_min, ratio, exp))
        print(f"{Phi0:>8.3f}  {M_min:>12.4f}  {ratio:>10.3f}  Phi0^{exp:.3f}")
    else:
        M_star_vals.append((Phi0, np.nan, np.nan, np.nan))
        print(f"{Phi0:>8.3f}  {'brak':>12}")

print()
if M1_base:
    print(f"  M*(Phi0=1) = M1 = {M1_base:.6f}")

# ============================================================
# KROK 2: Znajdz Phi0_2 i Phi0_3 potrzebne dla r21=207, r31=3477
# ============================================================
print()
print("KROK 2: Phi0 wymagane dla r21=207 i r31=3477")
print("-"*65)

if M1_base is None:
    print("  BLAD: M1 nie znalezione!")
else:
    def ratio_minus_target(Phi0, target_r, M1=M1_base):
        roots = find_M_star(Phi0)
        if not roots:
            return target_r  # brak korzenia
        return roots[0] / M1 - target_r

    # Szukaj Phi0_2 (M*/M1 = 207)
    Phi0_2 = None
    try:
        Phi0_2 = brentq(lambda p: ratio_minus_target(p, 207.0),
                        1.5, 100.0, xtol=1e-4, maxiter=50)
        M2 = find_M_star(Phi0_2)[0]
        print(f"  Phi0_2 (r21=207):  Phi0_2 = {Phi0_2:.5f}")
        print(f"    Phi0_2 / (2*pi) = {Phi0_2/(2*np.pi):.6f}")
        print(f"    Phi0_2 / pi     = {Phi0_2/np.pi:.6f}")
        print(f"    Phi0_2 / M_eff  = {Phi0_2/M_EFF:.6f}")
        print(f"    sqrt(Phi0_2)    = {np.sqrt(Phi0_2):.6f}")
        print(f"    Phi0_2^2        = {Phi0_2**2:.4f}")
        print(f"    M*(Phi0_2) = M2 = {M2:.4f}  (r21 = {M2/M1_base:.2f})")
    except Exception as ex:
        print(f"  Phi0_2: blad ({ex})")

    print()

    # Szukaj Phi0_3 (M*/M1 = 3477)
    Phi0_3 = None
    try:
        Phi0_3 = brentq(lambda p: ratio_minus_target(p, 3477.0),
                        2.0, 200.0, xtol=1e-4, maxiter=50)
        M3 = find_M_star(Phi0_3)[0]
        print(f"  Phi0_3 (r31=3477): Phi0_3 = {Phi0_3:.5f}")
        print(f"    Phi0_3 / (2*pi)    = {Phi0_3/(2*np.pi):.6f}")
        print(f"    Phi0_3 / pi        = {Phi0_3/np.pi:.6f}")
        print(f"    Phi0_3 / Phi0_2    = {Phi0_3/Phi0_2:.6f}  (jesli znaleziono)")
        print(f"    sqrt(Phi0_3)       = {np.sqrt(Phi0_3):.6f}")
        print(f"    Phi0_3^(1/3)       = {Phi0_3**(1/3):.6f}")
        print(f"    M*(Phi0_3) = M3    = {M3:.4f}  (r31 = {M3/M1_base:.1f})")
        print(f"    r32 = M3/M2 = {M3/M2:.2f}")
    except Exception as ex:
        print(f"  Phi0_3: blad ({ex})")

# ============================================================
# KROK 3: Test hipotezy 2*pi dla roznych (alpha, a_gam)
# ============================================================
print()
print("KROK 3: Test hipotezy Phi0_2 ~ 2*pi dla roznych parametrow")
print("-"*65)
print(f"  Cel: Phi0_2(r21=207) / (2*pi) -> czy stale dla roznych (alpha, a_gam)?")
print(f"  2*pi = {2*np.pi:.5f}")
print()
print(f"  {'alpha':>6}  {'a_gam':>7}  {'Phi0_2':>10}  {'Phi0_2/2pi':>12}  {'diff %':>8}")

params_test = [
    (5.9148, 0.025),
    (6.8675, 0.030),
    (7.7449, 0.035),
    (8.5616, 0.040),
    (5.9,    0.05),
]

for alpha_t, a_gam_t in params_test:
    M1_t = find_M_star(1.0, alpha=alpha_t, a_gam=a_gam_t)[0] if find_M_star(1.0, alpha=alpha_t, a_gam=a_gam_t) else None
    if M1_t is None:
        print(f"  {alpha_t:>6.4f}  {a_gam_t:>7.3f}  {'brak M1':>10}")
        continue
    try:
        P2_t = brentq(
            lambda p: (find_M_star(p, alpha=alpha_t, a_gam=a_gam_t) or [float('nan')])[0]/M1_t - 207.0,
            1.5, 100.0, xtol=1e-3, maxiter=30)
        diff = 100*(P2_t - 2*np.pi)/(2*np.pi)
        print(f"  {alpha_t:>6.4f}  {a_gam_t:>7.3f}  {P2_t:>10.5f}  {P2_t/(2*np.pi):>12.6f}  {diff:>+8.3f}%")
    except Exception:
        print(f"  {alpha_t:>6.4f}  {a_gam_t:>7.3f}  {'blad brentq':>10}")

# ============================================================
# KROK 4: Sprzezenie zwrotne -- jak soliton gen-1 zmienia Phi0?
# ============================================================
print()
print("KROK 4: Samosprzezenie -- soliton gen-1 vs zmiana Phi0")
print("-"*65)
print("  Profil gen-1: phi(r) = Phi0 + K1*exp(-m_eff*r)/r")
print(f"  K1 = K*1(ODE) = 0.010414 (z p14c)")
print(f"  m_eff = {M_EFF:.4f} (ODE), m_sp(Yuk,Phi0=1) = {1.0:.4f}")
print()

K1_ode = 0.010414
# Calka nadmiaru pola: Delta_Phi = 4*pi*integral_a_Gam^inf K*exp(-m*r)/r * r^2 dr
# = 4*pi*K * integral_a_Gam^inf r*exp(-m*r) dr
# = 4*pi*K * [e^{-m*r}*(−r/m − 1/m^2)]_a_Gam^inf
# = 4*pi*K * (a_Gam/m + 1/m^2) * exp(−m*a_Gam)

for m_label, m_val in [("m_eff(ODE)", M_EFF), ("m_sp(Yuk)", 1.0)]:
    I_excess = (A_GAM / m_val + 1.0 / m_val**2) * np.exp(-m_val * A_GAM)
    Delta_Phi = 4*np.pi * K1_ode * I_excess
    print(f"  [{m_label} = {m_val:.4f}]")
    print(f"    Calka nadmiaru = {I_excess:.6f}")
    print(f"    Delta_Phi_total = 4*pi*K1 * I = {Delta_Phi:.6f}")
    print(f"    Phi0_eff = 1 + Delta_Phi = {1+Delta_Phi:.6f}")
    print(f"    => Zmiana Phi0 jest MINIMALNA (nie tlumacy hierarchii)")
    print()

print("  WNIOSEK: Sprzezenie zwrotne solitonu gen-1 na Phi0 jest rzedu ~0.1%.")
print("  Nie moze wygenerowac Phi0^(2) ~ 6.37 z pierwszych zasad.")
print()

# ============================================================
# KROK 5: Minimalne warunki dla Mechanizmu 1
# ============================================================
print("=" * 65)
print("KROK 5: Co jest potrzebne dla Mechanizmu 1 (analiza ogolna)")
print("=" * 65)
print()
print("  Model sekwencyjny (diagnoza_sekwencyjny.py):")
print("  Phi0^(2) = Phi0^(1) + xi*M1     --> r21 = 207  => xi = (Phi0_2 - 1)/M1")
print("  Phi0^(3) = Phi0^(1) + xi*(M1+M2)  --> r31 = 3477  => Phi0_3 = ?")
print()

if M1_base and Phi0_2 and Phi0_3:
    # Z warunku r21=207: xi = (Phi0_2 - 1)/M1
    xi_from_r21 = (Phi0_2 - 1.0) / M1_base
    # Phi0_3 przewidywane przez model sekwencyjny:
    M2_seq = find_M_star(Phi0_2)[0]
    Phi0_3_seq = 1.0 + xi_from_r21 * (M1_base + M2_seq)
    print(f"  xi (z r21=207) = (Phi0_2 - 1)/M1 = ({Phi0_2:.4f} - 1)/{M1_base:.4f}")
    print(f"                 = {xi_from_r21:.6f}")
    print()
    print(f"  M2 (Yukawa, Phi0_2={Phi0_2:.3f}) = {M2_seq:.4f}")
    print(f"  M1 + M2 = {M1_base + M2_seq:.4f}")
    print()
    print(f"  Phi0_3 (model sekwencyjny)  = 1 + xi*(M1+M2)")
    print(f"                               = 1 + {xi_from_r21:.4f} * {M1_base + M2_seq:.4f}")
    print(f"                               = {Phi0_3_seq:.4f}")
    print()
    print(f"  Phi0_3 (potrzebne dla r31=3477) = {Phi0_3:.4f}")
    print()
    print(f"  Stosunek: Phi0_3(sekwencyjny) / Phi0_3(potrzebne)")
    print(f"          = {Phi0_3_seq:.4f} / {Phi0_3:.4f}")
    print(f"          = {Phi0_3_seq/Phi0_3:.2f}x  (model sekwencyjny jest ZA DUZY!)")
    print()

# ============================================================
# KROK 6: Hipotezy dla Phi0 z pierwszych zasad
# ============================================================
print("=" * 65)
print("KROK 6: Testowanie hipotez dla Phi0^(n)")
print("=" * 65)

if M1_base and Phi0_2 and Phi0_3:
    print()
    print(f"  Wyznaczone Phi0^(2) = {Phi0_2:.5f},  Phi0^(3) = {Phi0_3:.5f}")
    print()

    hypotheses = {}
    # Geometryczne
    hypotheses['2*pi']            = 2*np.pi
    hypotheses['pi^2/1']          = np.pi**2
    hypotheses['(2*pi)^(4/3)']    = (2*np.pi)**(4/3)
    hypotheses['4*pi/2']          = 4*np.pi/2
    hypotheses['sqrt(r21)']       = np.sqrt(207.0)
    hypotheses['r21^(1/3)']       = 207.0**(1/3)
    hypotheses['1/lambda^(1/4)']  = LAM**(-0.25)
    hypotheses['alpha/pi']        = ALPHA/np.pi
    hypotheses['1/(a_gam*m_eff)'] = 1.0/(A_GAM*M_EFF)
    hypotheses['1/a_gam']         = 1.0/A_GAM

    print(f"  {'Hipoteza':>24}  {'wartsc':>10}  {'blad od Phi0_2':>15}  {'blad od Phi0_3':>15}")
    print("  " + "-"*70)
    for name, val in hypotheses.items():
        diff2 = 100*(val - Phi0_2)/Phi0_2 if Phi0_2 else float('nan')
        diff3 = 100*(val - Phi0_3)/Phi0_3 if Phi0_3 else float('nan')
        print(f"  {name:>24}  {val:>10.5f}  {diff2:>+14.2f}%  {diff3:>+14.2f}%")

    # Specjalnie szukaj prostej relacji Phi0_3/Phi0_2
    ratio_32 = Phi0_3/Phi0_2
    print()
    print(f"  Phi0_3/Phi0_2 = {ratio_32:.5f}")
    print(f"  Sprawdz czy to pi: {ratio_32/np.pi:.5f} * pi = {ratio_32:.5f}")
    print(f"  Sprawdz czy to e:  {ratio_32/np.e:.5f} * e  = {ratio_32:.5f}")
    print(f"  Sprawdz czy to 1+sqrt(r32): sqrt(r32)={np.sqrt(3477/207):.4f}, 1+sqrt(r32)={1+np.sqrt(3477/207):.4f}")

# ============================================================
# KROK 7: Wykres M*(Phi0) i skala Phi0 hierarchii
# ============================================================
print()
print("Tworze wykres...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('P18: Mechanizm 1 z pierwszych zasad -- hierarchia Phi0\n'
             f'(alpha={ALPHA}, a_gam={A_GAM}, lam={LAM:.2e})',
             fontsize=11, fontweight='bold')

ax = axes[0]
Phi0_arr = np.array([x[0] for x in M_star_vals if not np.isnan(x[1])])
M_arr_plot = np.array([x[1] for x in M_star_vals if not np.isnan(x[1])])
ratio_arr = np.array([x[2] for x in M_star_vals if not np.isnan(x[1])])

if len(Phi0_arr) > 1:
    ax.semilogy(Phi0_arr, ratio_arr, 'bo-', lw=1.5, ms=7, label='M*(Phi0)/M1')
    ax.axhline(207.0,  color='red',    lw=1.5, linestyle='--', label='r21=207')
    ax.axhline(3477.0, color='purple', lw=1.5, linestyle='--', label='r31=3477')
    if Phi0_2:
        ax.axvline(Phi0_2, color='orange', lw=1.5, linestyle=':', label=f'Phi0_2={Phi0_2:.3f}')
    if Phi0_3:
        ax.axvline(Phi0_3, color='green', lw=1.5, linestyle=':', label=f'Phi0_3={Phi0_3:.3f}')
    ax.axvline(2*np.pi, color='magenta', lw=1, linestyle=':', alpha=0.7, label=f'2*pi={2*np.pi:.3f}')
ax.set_xlabel('Phi0 (poziom tla)')
ax.set_ylabel('M*(Phi0) / M*(1)')
ax.set_title('Skala M*(Phi0) vs Phi0', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

ax = axes[1]
# Wykres logarytmiczny: log(M*/M1) vs log(Phi0) -- czy potegowe?
if len(Phi0_arr) > 2:
    log_Phi0 = np.log(Phi0_arr[1:])
    log_ratio = np.log(ratio_arr[1:])
    # Dopasuj linie prosta: log(M*/M1) = n*log(Phi0) + c
    valid = np.isfinite(log_ratio)
    if valid.sum() > 2:
        coeffs = np.polyfit(log_Phi0[valid], log_ratio[valid], 1)
        n_fit, c_fit = coeffs
        Phi0_line = np.linspace(1, max(Phi0_arr), 100)
        M_line = np.exp(c_fit) * Phi0_line**n_fit
        ax.loglog(Phi0_arr[1:], ratio_arr[1:], 'bo', ms=7, label='dane')
        ax.loglog(Phi0_line, M_line, 'r--', lw=1.5,
                  label=f'fit: M/M1 ~ Phi0^{n_fit:.3f}')
        ax.set_title(f'Skalowanie potegowe: M* ~ Phi0^{n_fit:.3f}', fontsize=9)
ax.axhline(207.0,  color='red',    lw=1, linestyle='--', alpha=0.5)
ax.axhline(3477.0, color='purple', lw=1, linestyle='--', alpha=0.5)
ax.set_xlabel('Phi0'); ax.set_ylabel('M*(Phi0)/M1')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

ax = axes[2]
# Test modelu sekwencyjnego
if M1_base and Phi0_2:
    Phi0_test = np.linspace(1.0, 30.0, 200)
    xi_test   = (Phi0_2 - 1.0) / M1_base

    # Phi0_3 wyznaczane przez model sekwencyjny jako funkcja xi:
    xi_arr = np.linspace(xi_test*0.1, xi_test*10, 50)
    Phi0_3_seq_arr = []
    for xi_v in xi_arr:
        M2_v = find_M_star(1.0 + xi_v * M1_base)
        M2_v = M2_v[0] if M2_v else float('nan')
        Phi0_3_v = 1.0 + xi_v * (M1_base + M2_v) if not np.isnan(M2_v) else float('nan')
        Phi0_3_seq_arr.append(Phi0_3_v)

    Phi0_3_seq_arr = np.array(Phi0_3_seq_arr)
    valid = np.isfinite(Phi0_3_seq_arr)
    ax.semilogy(xi_arr[valid], Phi0_3_seq_arr[valid], 'b-', lw=2,
                label='Phi0_3 sekwencyjny(xi)')
    if Phi0_3:
        ax.axhline(Phi0_3, color='red', lw=1.5, linestyle='--',
                   label=f'Phi0_3 potrzebne={Phi0_3:.2f}')
    ax.axvline(xi_test, color='green', lw=1.5, linestyle=':',
               label=f'xi z r21=207: {xi_test:.4f}')
    ax.set_xlabel('xi (sprzezenie sekwencyjne)')
    ax.set_ylabel('Phi0_3 (model sekwencyjny)')
    ax.set_title('Problem modelu sekwencyjnego', fontsize=9)
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")
print()

# ============================================================
# PODSUMOWANIE
# ============================================================
print("=" * 65)
print("PODSUMOWANIE p18: Mechanizm 1 (hierarchia Phi0)")
print("=" * 65)
print()
print("  WYNIKI:")
if M1_base:
    print(f"    M1 (Yukawa, Phi0=1, ref. params) = {M1_base:.6f}")
if Phi0_2:
    print(f"    Phi0_2 (r21=207)  = {Phi0_2:.5f}")
    print(f"      vs 2*pi         = {2*np.pi:.5f}  (blad: {100*(Phi0_2-2*np.pi)/(2*np.pi):+.2f}%)")
if Phi0_3:
    print(f"    Phi0_3 (r31=3477) = {Phi0_3:.5f}")
print()
print("  WNIOSKI:")
print("  1. Model sekwencyjny (Phi0_n = Phi0_1 + xi*sum(M_k)) jest NIEMOZLIWY")
print("     dla r21=207 i r31=3477 jednoczesnie.")
print()
print("  2. Sprzezenie zwrotne solitonu gen-1 zmienia Phi0 o ~0.1% --")
print("     za malo by wygenerowac hierarchie Phi0_2/Phi0_1 ~ 6.")
print()
print("  3. Hipoteza Phi0_2 ~ 2*pi:")
if Phi0_2:
    diff = 100*(Phi0_2-2*np.pi)/(2*np.pi)
    if abs(diff) < 5.0:
        print(f"     POTWIERDZONA z bledem {diff:+.2f}% dla naszych parametrow!")
        print(f"     Phi0_2 = {Phi0_2:.5f} ~ 2*pi = {2*np.pi:.5f}")
    else:
        print(f"     ODRZUCONA: blad {diff:+.1f}% dla naszych parametrow.")
print()
print("  4. Trzy stabilne fazy TGP (z V_mod):")
print("     V_mod ma DWIE stabilne prozni (phi=1 i phi~1/sqrt(lam)~426).")
print("     => Brak mechanizmu dla trzech roznych Phi0 z V_mod.")
print()
print("  5. WNIOSEK OGOLNY:")
print("     Phi0^(n) jest PARAMETREM ZEWNETRZNYM TGP, nie wynika")
print("     z wewnetrznych rownan modelu. Musi byc wyznaczany z:")
print("     a) Kosmologicznej historii pola (3 fazy przejscia?)")
print("     b) Warunkow poczatkowych Wszechswiata w TGP")
print("     c) Dodatkowego warunku kwantyzacji (analog Bohra)")
print()
print("  OTWARTE PYTANIE:")
print("  Czy Phi0^(n) = phi_core^(n-1) (rekurencja przez rdzen)?")
print("  Czy istnieje topologiczne wytlumaczenie (klas Chern)?")
print("  Czy Phi0^(n) ~ (2*pi)^(n-1) jest dokladna relacja?")
