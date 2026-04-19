"""
ps9_mgb2_two_gap.py  -  Program P6.A #1

Cel:
  MgB2 ma dwa pasma parujace:
    sigma-band (silne, 2Delta_sigma = 13 meV, T_c_intra ~50 K)
    pi-band   (slabsze, 2Delta_pi = 4 meV, T_c_intra ~10 K)
  Obserwowane T_c = 39 K jest mieszanka - BCS single-gap dalby 50K.

Hipoteza TGP:
  sigma-band: bond B-B covalent, sprzezenie A_sp (sp^2) + A_mu (sigma-like multi-e)
              J_sigma = C_0 * (A_sp^2 + A_mu^2) - wzmocnienie hybrydyzacji
  pi-band:   bond B-B pi (out-of-plane), sprzezenie A_sp slabsze
              J_pi = C_0 * (A_sp^2) * (1 - eta_pi)   eta_pi ~ 0.4
  T_c (observed) = f(T_c^sigma, T_c^pi)
    opcje: max(), harmonic mean, lub BCS-coupled = (T_sigma * T_pi)^(1/2)
"""

import numpy as np

# Stale
K_B = 8.617333e-5  # eV/K
A_BOHR_ANGSTROM = 0.52917721067
a_star_tgp_A = 7.725 * A_BOHR_ANGSTROM  # 4.088 A

# P4 amplitudy
A_e_P4  = 0.124587
A_mu_P4 = 0.472198
A_tau_P4 = 0.956027

# ps5 5c
Lambda_E_5c = 0.1309  # meV
A_map_5c = {"s": -0.1110, "sp": 0.2067, "d": 0.3096, "f": 2.0336}
sigma_a_5c = 2.5856
C_0 = 48.8222

def k_d(z):
    return {4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}.get(z, 2.936)

def M_gauss(a_A, sigma=sigma_a_5c):
    n = max(1, round(a_A / a_star_tgp_A))
    d = a_A - n * a_star_tgp_A
    return np.exp(-d**2 / sigma**2)

# =============================================================
# MgB2 geometria
# =============================================================
# struktura AlB2: heksagonalna, a=3.086 A (in-plane B-B), c=3.524 A
# Kazdy B ma 3 sasiadow sigma (w plaszczyznie), 1 sasiada pi (out-plane)
# z_sigma = 3, z_pi = 1 (lub 2 jesli liczymy gore i dol)

a_inplane = 3.086       # B-B in-plane
a_outplane = 3.524      # c-axis Mg-B

Tc_obs = 39.0

print("=" * 78)
print("  ps9_mgb2_two_gap.py")
print("=" * 78)
print()
print(f"  MgB2: a_in-plane = {a_inplane} A, a_c = {a_outplane} A, T_c_obs = {Tc_obs} K")
print()

# =============================================================
# Part A. Sigma-band (in-plane B-B)
# =============================================================
print("=" * 78)
print("  Part A. Sigma-band (B-B in-plane covalent)")
print("=" * 78)

# Efektywny z dla sigma: heksagonalna siatka 2D z 3 NN, k_d = ?
# Dla hexagonal 2D, XY 3D-like coupling dla warstw: uzywamy interpolacji
# z=3 nie jest w tabeli, extrapolation: k_d(3) ~ 0.5 (ponizej BKT progu)
# ALE w MgB2 jest slabe coupling c-axis -> quasi-3D z ~ 3+2 = 5 efektywne
z_sigma = 5  # hexagonal 2D (z=3) + weak c-axis (z~2)
k_sigma = np.interp(z_sigma, [4, 6, 8], [0.893, 2.202, 2.936])

# Wzmocnienie sigma z hybrydyzacji sp^2:
# hipoteza: A_sigma^2 = A_sp^2 + alpha * A_mu^2,  alpha < 1
# test rozny alpha
print(f"\n  {'alpha':>6}  {'A_sigma^2':>10}  {'J_sigma':>8}  {'T_sigma[K]':>11}")
print(f"  {'------':>6}  {'----------':>10}  {'--------':>8}  {'-----------':>11}")

M_in = M_gauss(a_inplane)
best_sigma = None
for alpha in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
    A_sig_sq = A_map_5c["sp"]**2 + alpha * A_map_5c["d"]**2
    J_sig = C_0 * A_sig_sq
    Tc_sub = k_sigma * J_sig * M_in
    Tc_sig = Tc_sub * (Lambda_E_5c * 1e-3) / K_B
    print(f"  {alpha:>6.2f}  {A_sig_sq:>10.4f}  {J_sig:>8.2f}  {Tc_sig:>11.2f}")
    if best_sigma is None or abs(Tc_sig - 50.0) < abs(best_sigma[1] - 50.0):
        best_sigma = (alpha, Tc_sig)

print(f"\n  Cel: T_sigma ~ 50 K (wartosc pure sigma-band pre-pi-mixing).")
print(f"  Najlepszy alpha = {best_sigma[0]:.2f} -> T_sigma = {best_sigma[1]:.2f} K")

# =============================================================
# Part B. Pi-band (out-of-plane)
# =============================================================
print()
print("=" * 78)
print("  Part B. Pi-band (B-B out-of-plane pi)")
print("=" * 78)

z_pi = 2       # out-of-plane tylko
k_pi = 0.30    # ekstrapolacja dla z=2 (1D chain)

# eta_pi redukuje A_sp dla pi
print(f"\n  {'eta_pi':>7}  {'A_pi^2':>8}  {'J_pi':>6}  {'T_pi[K]':>8}")
print(f"  {'-------':>7}  {'--------':>8}  {'------':>6}  {'--------':>8}")

M_out = M_gauss(a_outplane)
for eta in [0.0, 0.2, 0.4, 0.6, 0.8]:
    A_pi_sq = A_map_5c["sp"]**2 * (1 - eta)
    J_pi = C_0 * A_pi_sq
    Tc_sub = k_pi * J_pi * M_out
    Tc_pi = Tc_sub * (Lambda_E_5c * 1e-3) / K_B
    print(f"  {eta:>7.2f}  {A_pi_sq:>8.4f}  {J_pi:>6.2f}  {Tc_pi:>8.3f}")

# =============================================================
# Part C. Laczenie gap'ow - 4 formuly
# =============================================================
print()
print("=" * 78)
print("  Part C. Laczenie sigma + pi -> T_c(observed)")
print("=" * 78)
print()

# Wybierz konkretne T_sigma, T_pi
alpha_best = best_sigma[0]
A_sig_sq = A_map_5c["sp"]**2 + alpha_best * A_map_5c["d"]**2
J_sigma = C_0 * A_sig_sq
T_sigma = k_sigma * J_sigma * M_in * Lambda_E_5c * 1e-3 / K_B

eta_best = 0.4  # umiarkowana redukcja pi
A_pi_sq = A_map_5c["sp"]**2 * (1 - eta_best)
J_pi = C_0 * A_pi_sq
T_pi = k_pi * J_pi * M_out * Lambda_E_5c * 1e-3 / K_B

print(f"  T_sigma = {T_sigma:.2f} K (alpha={alpha_best})")
print(f"  T_pi    = {T_pi:.2f} K (eta={eta_best})")
print()
print("  Formula kombinowania:")
print(f"    max(T_s, T_pi)       = {max(T_sigma, T_pi):.2f} K")
print(f"    srednia              = {0.5*(T_sigma+T_pi):.2f} K")
print(f"    srednia harm         = {2*T_sigma*T_pi/(T_sigma+T_pi):.2f} K")
print(f"    geometric mean       = {np.sqrt(T_sigma*T_pi):.2f} K")
w_sigma = 0.8
print(f"    wazona (w_s=0.8)     = {w_sigma*T_sigma+(1-w_sigma)*T_pi:.2f} K")

# McMillan-like: T = T_s * (1 - (1 - T_pi/T_s)^gamma) gdy pi tlumi s
print()
print("  Hipoteza TGP: pi-band TLUMI sigma przez konkurencje o faze koherentna.")
print("  Formula tlumienia:")
for gamma in [0.5, 1.0, 1.5]:
    if T_pi < T_sigma:
        T_eff = T_sigma - (T_sigma - T_pi) * gamma / 2
        print(f"    T_s - (T_s - T_pi)*{gamma}/2 = {T_eff:.2f} K")

# =============================================================
# Part D. Werdykt
# =============================================================
print()
print("=" * 78)
print("  Part D. Werdykt ps9")
print("=" * 78)
print()
print(f"  Cel: T_c(MgB2) = {Tc_obs} K")
print()
print(f"  ps5 5c single-band: 8.16 K -> off 5x")
print(f"  ps9 pure sigma band:       {T_sigma:.2f} K")
print(f"  ps9 sigma-pi weighted avg: {w_sigma*T_sigma+(1-w_sigma)*T_pi:.2f} K")
print()
if abs(T_sigma - Tc_obs) / Tc_obs < 0.3:
    print("  Sigma-band ALONE wystarczy do wyjasnienia T_c MgB2 w TGP.")
    print("  Pi-band jest minorytetow redukujacym efekt w ~10% pozniej.")
elif abs(w_sigma*T_sigma+(1-w_sigma)*T_pi - Tc_obs) / Tc_obs < 0.3:
    print("  Wazona kombinacja sigma+pi pasuje z toleranсja 30%.")
else:
    print("  Model dwu-gap'owy ps9 tez nie trafia MgB2.")
    print("  Moze potrzebne inne A amplitudy (P4 nie mapuje B-sp3-orbital dokladnie).")

print()
print("=" * 78)
print("  ps9 complete.")
print("=" * 78)
