"""
ex7_nbody_accumulation.py
==========================
Dwa otwarte pytania TGP:

PYTANIE 1: Dla ilu cial N efekty trojcialowe staja sie dominujace?
  V3_total = C(N,3) * V3_triplet ~ N^3 * V3
  V2_total = C(N,2) * V2_pair   ~ N^2 * V2
  stosunek: (N-2)/3 * |V3/V2|_per_triplet -> 1  gdy N = N_crit

  N_crit = 3 / |V3_triplet/V2_pair|

PYTANIE 2: Promieniowanie skalarne TGP vs GR
  Pole TGP jest masywne (m_sp > 0) => istnieje odciecie czestotliwosci:
  - jesli czestotliwosc orbitalna omega < m_sp: BRAK promieniowania
  - jesli omega > m_sp: promieniowanie skalarne, moc:
      P = C^2 / (6*pi) * |a|^2 * sqrt(1 - m_sp^2/omega^2)

  Dla GR: P_GR ~ G*m^2*a^2 (zawsze, niezaleznie od okresu)

  Wniosek: krotkoookresowe orbity (omega >> m_sp) promieniuja tak samo
           jak GR (po przeskalowaniu), dlugoookresowe (omega < m_sp)
           w ogole nie promieniuja.
"""

import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from nbody.three_body_force_exact import yukawa_overlap_exact
from nbody.dynamics_v2 import potential_tgp
from nbody.tgp_field import screening_mass

beta, gamma = 1.0, 1.0
m = screening_mass(beta, gamma)

# ===========================================================================
# CZESC 1: Tabela N_crit
# ===========================================================================
print("=" * 65)
print("CZESC 1: Kiedy V3_total ~ V2_total? (N_crit)")
print("=" * 65)
print()
print("Schemat: N identycznych cial przy sredniej separacji d.")
print("V2_total = C(N,2)*V2_pair ~ N^2 * V2")
print("V3_total = C(N,3)*V3_triplet ~ N^3 * V3")
print("N_crit = 3 / |V3_triplet/V2_pair| + 2")
print()

C_list = [0.05, 0.10, 0.15, 0.20]
t_list = [1.0, 1.5, 2.0, 2.5, 3.0]

ratios_grid = np.zeros((len(t_list), len(C_list)))
Ncrit_grid  = np.zeros((len(t_list), len(C_list)))

header = "  t=m*d  " + "".join(f"  C={c:.2f}(N_crit)" for c in C_list)
print(header)
print("-" * len(header))

for it, t in enumerate(t_list):
    d = t / m
    row = f"  {t:.1f}    "
    for ic, C in enumerate(C_list):
        pos = np.array([[0,0,0],[d,0,0],[d/2,d*np.sqrt(3)/2,0]], dtype=float)
        Cv  = np.array([C, C, C])
        V2  = potential_tgp(pos, Cv, beta=beta, gamma=gamma)
        I_Y = yukawa_overlap_exact(d, d, d, m, n_quad=50)
        V3  = -6.0 * gamma * C**3 * I_Y
        r   = abs(V3 / V2)
        Nc  = int(3.0 / r + 2) if r > 1e-12 else 999999
        ratios_grid[it, ic] = r
        Ncrit_grid[it, ic]  = Nc
        row += f"  {r:.3%} ({Nc:5d})"
    print(row)

print()
print("Interpretacja:")
print("  Dla N > N_crit energia trojcialowa przekracza parami")
print("  => TGP jakosciowo rozni sie od sumy par (Grawitacja Newtonowska)")
print()

# ===========================================================================
# CZESC 2: Promieniowanie skalarne - odciecie czestotliwosci
# ===========================================================================
print("=" * 65)
print("CZESC 2: Promieniowanie skalarne TGP - odciecie m_sp")
print("=" * 65)
print()
print(f"Masa skalaru: m_sp = sqrt(beta) = {m:.4f}")
print(f"Dlugosc ekranowania: 1/m_sp = {1/m:.4f}")
print()

# Moc promieniowania dla oscylatora harmonicznego o amplitudzie A, omega
# x(t) = A*cos(omega*t), |a| = A*omega^2
# P_TGP = C^2/(6*pi) * A^2*omega^4 * sqrt(1 - m_sp^2/omega^2) dla omega > m_sp
# P_TGP = 0 dla omega <= m_sp

omega_vals = np.linspace(0.1, 5.0, 500)
C_rad = 1.0
A_rad = 1.0

P_TGP = np.where(
    omega_vals > m,
    C_rad**2 / (6*np.pi) * A_rad**2 * omega_vals**4 * np.sqrt(1 - m**2/omega_vals**2),
    0.0
)
# Dla porownania: masowe pole skalarne bez odciecia
P_massless = C_rad**2 / (6*np.pi) * A_rad**2 * omega_vals**4

print("  Moc promieniowania P(omega) dla C=1, A=1:")
print(f"  {'omega':>8}  {'P_TGP':>12}  {'P_massless':>12}  {'supresja':>10}")
print("-" * 50)
for om in [0.5, 0.8, 1.0, 1.1, 1.5, 2.0, 3.0, 5.0]:
    idx = np.argmin(np.abs(omega_vals - om))
    pt  = P_TGP[idx]
    pm  = P_massless[idx]
    sup = pt/pm if pm > 0 else 0.0
    flag = "  ZERO (omega<m_sp)" if om < m else ""
    print(f"  {om:>8.2f}  {pt:>12.4e}  {pm:>12.4e}  {sup:>10.4f}{flag}")

print()
print("  Wniosek TGP vs GR:")
print(f"  - Orbity z omega < {m:.3f} (okres T > {2*np.pi/m:.2f}) NIE promieniuja")
print("  - Dugie orbity sa dokladnie stabilne (brak zaniku)")
print("  - Krotkie orbity (omega >> m_sp) promieniuja jak masowe pole skalarne")
print()

# ===========================================================================
# CZESC 3: Konkretna predykcja - galaktyka
# ===========================================================================
print("=" * 65)
print("CZESC 3: Predykcja dla ukladu N-cialowego")
print("=" * 65)
print()
print("Pytanie: przy jakim N skupisko cial TGP 'czuje' swoja trojcialoscosc?")
print()

# Przy t=m*d=2, C=0.15: N_crit ~ 500
# Przy t=m*d=1, C=0.15: N_crit ~ ?
for C in [0.10, 0.15]:
    for t in [1.0, 2.0]:
        d = t / m
        pos = np.array([[0,0,0],[d,0,0],[d/2,d*np.sqrt(3)/2,0]], dtype=float)
        Cv  = np.array([C, C, C])
        V2  = potential_tgp(pos, Cv, beta=beta, gamma=gamma)
        I_Y = yukawa_overlap_exact(d, d, d, m, n_quad=50)
        V3  = -6.0 * gamma * C**3 * I_Y
        r   = abs(V3 / V2)
        Nc  = int(3.0 / r + 2)
        print(f"  C={C:.2f}, t=m*d={t:.1f}: N_crit = {Nc}")
        print(f"    (skupisko > {Nc} cial przy tej gestosci: energia trojakowa dominuje)")

# ===========================================================================
# WYKRESY
# ===========================================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Wykres 1: N_crit vs t dla roznych C
ax = axes[0]
t_fine = np.linspace(0.5, 4.0, 80)
for C in C_list:
    Nc_vals = []
    for t in t_fine:
        d = t / m
        pos = np.array([[0,0,0],[d,0,0],[d/2,d*np.sqrt(3)/2,0]], dtype=float)
        Cv  = np.array([C, C, C])
        V2  = potential_tgp(pos, Cv, beta=beta, gamma=gamma)
        I_Y = yukawa_overlap_exact(d, d, d, m, n_quad=40)
        V3  = -6.0 * gamma * C**3 * I_Y
        r   = abs(V3 / V2)
        Nc  = 3.0 / r + 2 if r > 1e-14 else 1e6
        Nc_vals.append(Nc)
    ax.semilogy(t_fine, Nc_vals, lw=2, label=f"C = {C:.2f}")

ax.axhline(100,  color='gray', ls=':', lw=1.2, label="N=100")
ax.axhline(1000, color='silver', ls=':', lw=1.2, label="N=1000")
ax.set_xlabel("t = m·d (bezwymiarowa separacja)", fontsize=12)
ax.set_ylabel("N_crit (liczba cial)", fontsize=12)
ax.set_title("Kiedy V3_total > V2_total?", fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim(0.5, 4.0)

# Wykres 2: Moc promieniowania TGP vs masless
ax = axes[1]
ax.plot(omega_vals, P_massless, 'b--', lw=1.5, label="Pole masowe (m=0)")
ax.plot(omega_vals, P_TGP,      'r-',  lw=2.5, label=f"TGP (m_sp={m:.2f})")
ax.axvline(m, color='orange', ls='--', lw=2, label=f"odciecie omega=m_sp={m:.2f}")
ax.fill_between(omega_vals, 0, P_TGP, where=(omega_vals < m),
                alpha=0.15, color='red', label="strefa bez promieniowania")
ax.set_xlabel("czestotliwosc orbitalna omega", fontsize=12)
ax.set_ylabel("Moc promieniowania P(omega)", fontsize=12)
ax.set_title("Promieniowanie skalarne TGP: odciecie m_sp", fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim(0.1, 5.0)
ax.set_ylim(-0.05, max(P_massless)*1.05)

plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "ex7_nbody_predictions.png")
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print()
print(f"Wykres zapisany: {out}")
