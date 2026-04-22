"""
ex45_lyapunov_tgp_vs_newton.py
================================
Wykładnik Lyapunova: TGP vs Newton (3 ciała)

MOTYWACJA
---------
TGP zawiera barierę repulsywną V_β ~ 1/d² (reżim II) która może tłumić chaos
w układach 3-ciałowych. Pytanie: czy bariera repulsywna TGP redukuje wykładnik
Lyapunova λ w porównaniu do czystego Newtona?

METODA
------
1. Symulacja 3 ciał przez T = 5000 kroków (Runge-Kutta 4)
2. Perturbacja δ = 1e-8 w warunkach początkowych
3. Pomiar λ = lim(t→∞) ln|δ(t)/δ(0)| / t
4. Porównaj λ_TGP vs λ_Newton dla różnych C

PREDYKCJA
---------
λ_TGP < λ_Newton dla C bliskich progu wiązania (bariera tłumi zbliżenia)
Efekt wzmacnia się z C.

PARAMETRY
---------
- 3 ciała o równych masach w konfiguracji isosceles/random
- Masa: m = 1 [Planck]
- Siła źródłowa: C = 0.05, 0.10, 0.15, 0.20 [TGP] lub 0 [Newton]
- Czas: T = 5000 kroków, dt = 0.05
- Liczba perturbacji: Npert = 5 (uśrednianie po realizacjach)

Autor: TGP Analysis Session v27, 2026-03-22
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.special import k0 as K0, k1 as K1
import warnings
warnings.filterwarnings('ignore')

# -------------------------------------------------------------------
# Parametry globalne
# -------------------------------------------------------------------
BETA   = 1.0     # = gamma (warunek próżniowy N0-5)
GAMMA  = 1.0
M_SP   = 0.10    # [l_Pl^-1] masa pola — sektor Efimova
DT     = 0.05    # krok czasowy [t_Pl]
T_MAX  = 5000    # liczba kroków
DELTA  = 1e-8    # perturbacja początkowa
N_PERT = 5       # liczba perturbacji do uśrednienia

print("=" * 70)
print("EX45: Wykładnik Lyapunova — TGP vs Newton (3 ciała)")
print(f"      m_sp = {M_SP}, dt = {DT}, T = {T_MAX} kroków")
print("=" * 70)
print()

# -------------------------------------------------------------------
# Całka Feynmana dla V3 (uproszczona — przybliżenie siodłowe dla V3)
# -------------------------------------------------------------------
_pts, _wts = np.polynomial.legendre.leggauss(12)
_up = 0.5 * (1 + _pts); _uw = 0.5 * _wts
_UU, _VV = np.meshgrid(_up, _up, indexing='ij')
_WW = np.outer(_uw, _uw)
_A1 = _UU; _A2 = _VV * (1 - _UU); _A3 = (1 - _UU) * (1 - _VV)
_JAC = (1 - _UU)

def I_Y_scalar(d12, d13, d23, m):
    """Całka Feynmana I_triple dla trójki ciał."""
    D = _A1*_A2 + _A1*_A3 + _A2*_A3
    Q = _A2*d12**2 + _A1*d13**2 + _A3*d23**2
    good = (D > 1e-30) & (Q > 1e-30)
    u   = np.where(good, m * np.sqrt(Q / D), 1.0)
    val = np.where(good, D**(-1.5) * K0(u), 0.0)
    return 2.0 * np.sum(_WW * _JAC * val)

# -------------------------------------------------------------------
# Potencjały i siły
# -------------------------------------------------------------------

def V2_tgp(d, C, m):
    """Potencjał parowy TGP (Yukawa attracting + repulsive barrier)."""
    if d < 1e-6:
        return 1e6
    em  = np.exp(-m * d)
    # V_lin ~ -C^2/d * exp(-m*d) [przyciąganie]
    # V_beta ~ +C^2 * beta / d^2 [odpychanie]
    v_lin  = -C**2 * em / d
    v_beta = +C**2 * BETA * em / d**2
    v_gam  = -C**2 * GAMMA * em / d**3
    return v_lin + v_beta + v_gam

def F2_tgp(d, C, m):
    """Siła parowa TGP: F = -dV/dd (skalar, > 0 = odpychanie)."""
    if d < 1e-6:
        return 0.0
    em  = np.exp(-m * d)
    # dV_lin/dd:
    dV_lin  = C**2 * em * (1.0/d + m) / d
    dV_beta = -C**2 * BETA * em * (2.0/d + m) / d**2
    dV_gam  = C**2 * GAMMA * em * (3.0/d + m) / d**3
    return -(dV_lin + dV_beta + dV_gam)

def V3_approx(d12, d13, d23, C, m):
    """Przybliżenie siodłowe V3 (szybkie)."""
    P = d12 + d13 + d23
    if P < 1e-6:
        return 0.0
    return -6.0 * GAMMA * C**3 * 8.0 * np.pi**2 * np.exp(-m * P/2.0) / P

def newton_V2(d, C):
    """Potencjał newtonowski: V ~ -C^2/d (bez bariery)."""
    return -C**2 / max(d, 1e-4)

def newton_F2(d, C):
    """Siła newtonowska: F = -C^2/d^2 (< 0 = przyciąganie)."""
    return -C**2 / max(d, 1e-4)**2

# -------------------------------------------------------------------
# Symulator 3-ciałowy (2D dla prostoty: z = 0)
# -------------------------------------------------------------------

def compute_forces(pos, C, mode='tgp'):
    """Oblicz siły na 3 ciała."""
    forces = np.zeros_like(pos)
    for i in range(3):
        for j in range(i+1, 3):
            r = pos[j] - pos[i]
            d = np.linalg.norm(r)
            if d < 1e-6:
                continue
            rhat = r / d
            if mode == 'tgp':
                f_scalar = F2_tgp(d, C, M_SP)
            else:  # newton
                f_scalar = newton_F2(d, C)
            forces[i] -= f_scalar * rhat
            forces[j] += f_scalar * rhat
    return forces

def rk4_step(pos, vel, C, dt, mode='tgp'):
    """Jeden krok RK4."""
    m = 1.0  # jednostkowa masa

    def deriv(p, v):
        a = compute_forces(p, C, mode) / m
        return v, a

    v1, a1 = deriv(pos, vel)
    v2, a2 = deriv(pos + 0.5*dt*v1, vel + 0.5*dt*a1)
    v3, a3 = deriv(pos + 0.5*dt*v2, vel + 0.5*dt*a2)
    v4, a4 = deriv(pos + dt*v3, vel + dt*a3)

    pos_new = pos + dt/6.0 * (v1 + 2*v2 + 2*v3 + v4)
    vel_new = vel + dt/6.0 * (a1 + 2*a2 + 2*a3 + a4)
    return pos_new, vel_new

def run_trajectory(pos0, vel0, C, T, dt, mode='tgp'):
    """Uruchom symulację i zwróć trajektorię jako tablicę."""
    pos = pos0.copy()
    vel = vel0.copy()
    traj = np.zeros((T+1, 3, 2))
    traj[0] = pos
    for t in range(T):
        pos, vel = rk4_step(pos, vel, C, dt, mode)
        traj[t+1] = pos
    return traj

def compute_lyapunov(pos0, vel0, C, T, dt, delta, mode='tgp'):
    """
    Oblicz wykładnik Lyapunova metodą Benettin.
    Uruchamia referencyjną trajektorię i N_PERT perturbowanych,
    mierzy średni przyrost ln|δ(t)/δ(0)|/t.
    """
    lyap_estimates = []

    for _ in range(N_PERT):
        # Losowa perturbacja
        pert_pos = np.random.randn(3, 2)
        pert_pos /= np.linalg.norm(pert_pos)
        pert_pos *= delta

        pos_ref = pos0.copy()
        vel_ref = vel0.copy()
        pos_pert = pos0.copy() + pert_pos
        vel_pert = vel0.copy()

        log_sum = 0.0
        t_valid = 0

        for t in range(T):
            pos_ref, vel_ref   = rk4_step(pos_ref, vel_ref, C, dt, mode)
            pos_pert, vel_pert = rk4_step(pos_pert, vel_pert, C, dt, mode)

            # Dystans między trajektoriami
            diff = pos_pert - pos_ref
            d_cur = np.linalg.norm(diff)

            if d_cur < 1e-15 or not np.isfinite(d_cur):
                continue

            # Renormalizacja co 10 kroków
            if (t + 1) % 10 == 0:
                log_sum += np.log(d_cur / delta)
                t_valid += 1
                # Przeskaluj do delta
                pos_pert = pos_ref + diff / d_cur * delta
                vel_pert = vel_ref.copy()

        if t_valid > 0:
            lyap = log_sum / (t_valid * 10 * dt)
            lyap_estimates.append(lyap)

    if lyap_estimates:
        return np.mean(lyap_estimates), np.std(lyap_estimates)
    return np.nan, np.nan

# -------------------------------------------------------------------
# Konfiguracje startowe
# -------------------------------------------------------------------

def config_isosceles(d=5.0):
    """Konfiguracja izosceles: 1 na górze, 2 i 3 na dole."""
    pos = np.array([
        [0.0,  d * np.sqrt(3)/3],
        [-d/2, -d * np.sqrt(3)/6],
        [d/2,  -d * np.sqrt(3)/6],
    ])
    # Mała prędkość orbitalna (brak czystego spoczynku — stabilizacja)
    v0 = 0.05
    vel = np.array([
        [v0, 0.0],
        [-v0/2, -v0*np.sqrt(3)/2],
        [-v0/2,  v0*np.sqrt(3)/2],
    ])
    # Usuń środek masy
    vel -= vel.mean(axis=0)
    return pos, vel

def config_random(seed=42, d_mean=5.0):
    """Losowa konfiguracja 3 ciał w pudełku d_mean x d_mean."""
    rng = np.random.RandomState(seed)
    pos = rng.randn(3, 2) * d_mean / 2.0
    pos -= pos.mean(axis=0)
    vel = rng.randn(3, 2) * 0.1
    vel -= vel.mean(axis=0)
    return pos, vel

# -------------------------------------------------------------------
# GŁÓWNA ANALIZA: Lyapunov TGP vs Newton
# -------------------------------------------------------------------

C_values = [0.05, 0.10, 0.15, 0.20]
configs = [
    ("isosceles d=5", *config_isosceles(d=5.0)),
    ("isosceles d=8", *config_isosceles(d=8.0)),
    ("random seed=42", *config_random(seed=42)),
]

print("─" * 70)
print("WYNIKI: Wykładnik Lyapunova λ [t_Pl^-1]")
print("─" * 70)
print()

results = {}  # (C, mode, config_name) -> (lyap, std)

np.random.seed(0)

for cfg_name, pos0, vel0 in configs:
    print(f"Konfiguracja: {cfg_name}")
    print(f"  {'C':>6}  {'λ_TGP':>12}  {'σ_TGP':>8}  {'λ_Newton':>12}  {'σ_Newton':>8}  {'λ_TGP/λ_Newton':>14}")
    print(f"  {'─'*6}  {'─'*12}  {'─'*8}  {'─'*12}  {'─'*8}  {'─'*14}")

    for C in C_values:
        lam_tgp, std_tgp     = compute_lyapunov(pos0, vel0, C, T_MAX, DT, DELTA, mode='tgp')
        lam_newt, std_newt   = compute_lyapunov(pos0, vel0, C, T_MAX, DT, DELTA, mode='newton')

        results[(C, 'tgp', cfg_name)]    = (lam_tgp, std_tgp)
        results[(C, 'newton', cfg_name)] = (lam_newt, std_newt)

        if np.isfinite(lam_tgp) and np.isfinite(lam_newt) and lam_newt > 0:
            ratio = lam_tgp / lam_newt
            ratio_str = f"{ratio:.3f}"
        else:
            ratio_str = "N/A"

        print(f"  {C:>6.2f}  {lam_tgp:>12.5f}  {std_tgp:>8.5f}  {lam_newt:>12.5f}  {std_newt:>8.5f}  {ratio_str:>14}")

    print()

# -------------------------------------------------------------------
# SYNTEZA: Czy bariera TGP redukuje chaos?
# -------------------------------------------------------------------

print("─" * 70)
print("SYNTEZA: Bariera repulsywna TGP a chaos")
print("─" * 70)
print()

tgp_lt_newton_count = 0
total_valid = 0

ratios_by_C = {}
for C in C_values:
    ratios_C = []
    for cfg_name, _, _ in configs:
        lam_t = results.get((C, 'tgp', cfg_name), (np.nan,))[0]
        lam_n = results.get((C, 'newton', cfg_name), (np.nan,))[0]
        if np.isfinite(lam_t) and np.isfinite(lam_n) and lam_n > 1e-10:
            ratio = lam_t / lam_n
            ratios_C.append(ratio)
            total_valid += 1
            if ratio < 1.0:
                tgp_lt_newton_count += 1
    if ratios_C:
        ratios_by_C[C] = np.mean(ratios_C), np.std(ratios_C)
    else:
        ratios_by_C[C] = np.nan, np.nan

print(f"  {'C':>6}  {'<λ_TGP/λ_Newton>':>18}  {'σ':>8}  {'λ_TGP < λ_Newton?':>20}")
print(f"  {'─'*6}  {'─'*18}  {'─'*8}  {'─'*20}")
for C in C_values:
    mean_r, std_r = ratios_by_C[C]
    tgp_suppressed = mean_r < 1.0 if np.isfinite(mean_r) else False
    print(f"  {C:>6.2f}  {mean_r:>18.3f}  {std_r:>8.3f}  {'TAK' if tgp_suppressed else 'NIE':>20}")

print()
if total_valid > 0:
    frac = tgp_lt_newton_count / total_valid
    print(f"  Frakcja przypadków λ_TGP < λ_Newton: {frac:.1%} ({tgp_lt_newton_count}/{total_valid})")
    print()
    if frac > 0.6:
        print("  WNIOSEK: Bariera repulsywna TGP REDUKUJE chaos w większości konfiguracji.")
        print("  Predykcja λ_TGP < λ_Newton POTWIERDZONA (jakościowo).")
        verdict = "POTWIERDZONA"
    elif frac > 0.4:
        print("  WNIOSEK: Efekt mieszany — bariera TGP nie zawsze redukuje chaos.")
        print("  Zależy od konfiguracji startowej i wartości C.")
        verdict = "MIESZANA"
    else:
        print("  WNIOSEK: Bariera TGP nie redukuje chaosu w większości konfiguracji.")
        print("  Możliwe że zakres C lub konfiguracja jest zbyt daleka od progu.")
        verdict = "NIEUWIDOCZNIONA"

    print()
    print(f"  Status K_Lyap: {verdict}")

print()
print("─" * 70)
print("Uwaga: T=5000, dt=0.05, N_pert=5 — wyniki orientacyjne.")
print("Dla dokładności: T=50000, N_pert=20 (patrz ex45_highres.py).")
print("─" * 70)
print()
print("=" * 70)
print("EX45 DONE — Lyapunov TGP vs Newton")
print("=" * 70)
