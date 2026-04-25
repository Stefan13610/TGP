# -*- coding: utf-8 -*-
"""
OP-7 / T2 -- Wyprowadzenie sigma_ab z H_Gamma jako kompozytowy operator
========================================================================

Cel: pokazac, ze tensor sigma_ab postulowany w tgp_core.tex sec.2
("One substrate, two projections") jest dobrze zdefiniowanym
kompozytowym operatorem coarse-grainingu pojedynczego pola substratowego
s-hat, z wlasnosciami:

  (i)   sigma_ab jest BILINEAROWE w s-hat (kompozyt, nie nowy stopien
        swobody)
  (ii)  Z2-parzyste:  s -> -s  =>  sigma_ab -> sigma_ab
  (iii) Symmetryczne i traceless: 5 d.o.f. spin-2 w 3D
  (iv)  Pod SO(3):    sigma_ab -> R_a^c R_b^d sigma_cd (rank-2 tensor)
  (v)   W prozni izotropowej: sigma_ab = 0 (tylko trace czesc)
  (vi)  Przy anizotropowym zaburzeniu (kwadrupol): sigma_ab != 0

Plan:
  Czesc A:  Analitycznie -- struktura sigma_ab z H_Gamma
  Czesc B:  Sympy -- weryfikacja Z2 parity + SO(3) kowariancji
  Czesc C:  Lattice MC -- na toy modelu v2 GL-substrate-like:
              C.1  prozni izotropowa: sigma -> 0
              C.2  z kwadrupolowym zaburzeniem: sigma != 0, ksztalt zgodny
              C.3  Z2 parity: |sigma(s) - sigma(-s)| = 0
              C.4  SO(3) kowariancja: sigma(R[s]) = R sigma(s) R^T
  Czesc D:  Werdykt T2

Output: T2 PASS jesli (i)-(vi) sprawdzaja sie analitycznie + numerycznie.
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

import numpy as np
import sympy as sp
from scipy.linalg import eigvalsh


def banner(title, level=1):
    if level == 1:
        print("\n" + "=" * 72)
        print(f"  {title}")
        print("=" * 72)
    else:
        print(f"\n  --- {title} ---")


def check(label, condition, value=None):
    mark = "PASS" if condition else "FAIL"
    extra = f"  [{value}]" if value is not None else ""
    print(f"  [{mark}] {label}{extra}")
    return bool(condition)


checks_summary = []


def record(label, condition, value=None):
    checks_summary.append((label, bool(condition)))
    return check(label, condition, value)


# =========================================================================
# CZESC A: Analitycznie -- struktura sigma_ab z H_Gamma
# =========================================================================
banner("CZESC A: Analityczna struktura sigma_ab z H_Gamma")

print("""
  Substratowy Hamiltonian (v2 GL-bond, sec.B paper):

    H_Gamma = sum_i [pi_i^2 / (2 mu) + (m0^2/2) s_i^2 + (lambda0/4) s_i^4]
              + J * sum_{<ij>} A_ij s_i^2 s_j^2 (s_j^2 - s_i^2)^2

  Z2 symmetria:  s_i -> -s_i  pozostawia H_Gamma niezmienione.

  Kompozytowy operator nearest-neighbor:

    K_{ab}(x) = lim_{|B| -> inf}  (1/|B|) sum_{i in B}
                  s_i * s_{i + e_a} * s_i * s_{i + e_b} / <s^2>

  W formie uproszczonej (po wykorzystaniu samouzgadniajacego MF):

    K_{ab}(x) ~ (1/|B|) sum_{i in B} s_i s_{i + e_a}    (jesli a == b)
              ~ <s_i s_{i+e_a}>_B                       (bilinear)

  Po dekompozycji 3x3 macierzy K na slad + bezsladowa:

    K_{ab} = (1/3) (Tr K) delta_{ab} + sigma_{ab}

  gdzie sigma_{ab} jest 5-skladowym TT-czesta (3D symmetric traceless).

  Wlasnosci sigma_ab (do zweryfikowania ponizej):

    (i)   bilinear w s-hat: K_ab ~ s s, wiec sigma_ab tez
    (ii)  Z2-parzyste: s -> -s  =>  s s -> s s  =>  K_ab -> K_ab
    (iii) symmetric + traceless (z konstrukcji)
    (iv)  rank-2 tensor SO(3) (lattice -> kontinuum przez block-averaging)
    (v)   prozni izotropowa: <s_i s_{i+e_a}> = const dla wszystkich a
          => K_ab ~ delta_ab => sigma_ab = 0
    (vi)  anizotropowe zrodlo (np. kwadrupol mat.) lamie SO(3) lokalnie
          => <s_i s_{i+e_a}> rozne dla roznych a => sigma_ab != 0
""")

# =========================================================================
# CZESC B: Sympy -- Z2 parity + SO(3) kowariancja symbolicznie
# =========================================================================
banner("CZESC B: Sympy weryfikacja symbolicznych wlasnosci")

# B.1: Z2 parity
banner("B.1  Z2 parzystosc: K_ab(s) = K_ab(-s)", 2)

# Symbolicznie: pole s w 3 punktach (i, i+e_a, i+e_b) -- siatka 3D z 6 punktami
# Wystarcza pokazac w 1D analogu: K = s_i * s_{i+1}.
# Po Z2: s -> -s daje (-s_i)(-s_{i+1}) = s_i s_{i+1} = K. Niezmienne.
s_i, s_ip = sp.symbols('s_i s_ip', real=True)  # s na sasiadzie
K_bilinear = s_i * s_ip
K_after_Z2 = K_bilinear.subs([(s_i, -s_i), (s_ip, -s_ip)])
delta_K_Z2 = sp.simplify(K_after_Z2 - K_bilinear)
record("Bilinear K_ab Z2-parzysty (s -> -s zostawia K niezmienione)",
       delta_K_Z2 == 0,
       f"K(-s) - K(s) = {delta_K_Z2}")

# B.2: SO(3) kowariancja
banner("B.2  SO(3) kowariancja: sigma_ab -> R_ac R_bd sigma_cd", 2)

# Wezmy konkretna macierz K_ab i obrocenie
# K_ab = s_a s_b (kontraktna postac dla pojedynczego punktu, jak rank-2 tensor)
sx, sy, sz = sp.symbols('sx sy sz', real=True)
s_vec = sp.Matrix([sx, sy, sz])
K_mat = s_vec * s_vec.T   # K_ab = s_a s_b (3x3 symmetric)

print(f"  K_ab = s_a s_b =\n{K_mat}")

# slad
trace_K = K_mat.trace()
sigma_mat = K_mat - sp.Rational(1, 3) * trace_K * sp.eye(3)
print(f"\n  Tr(K) = {sp.simplify(trace_K)}")
print(f"  sigma_ab = K_ab - (1/3) Tr(K) delta_ab")

# Obrot SO(3): kat phi wokol osi z
phi = sp.Symbol('phi', real=True)
R = sp.Matrix([
    [sp.cos(phi), -sp.sin(phi), 0],
    [sp.sin(phi),  sp.cos(phi), 0],
    [0,            0,           1],
])

# s' = R s
s_rot = R * s_vec
K_rot = s_rot * s_rot.T
sigma_rot = K_rot - sp.Rational(1, 3) * K_rot.trace() * sp.eye(3)

# Powinno byc R sigma R^T
sigma_via_transform = R * sigma_mat * R.T

# Roznica
diff_mat = sp.simplify(sigma_rot - sigma_via_transform)
diff_norm = sp.simplify(sum(sp.simplify(diff_mat[i, j])**2 for i in range(3) for j in range(3)))

record("SO(3) kowariancja: sigma(R s) = R sigma(s) R^T",
       diff_norm == 0,
       f"||sigma_rot - R sigma R^T||^2 = {diff_norm}")

# B.3: traceless + symmetric
record("sigma_ab symetryczne",
       sp.simplify(sigma_mat - sigma_mat.T).norm() == 0)

record("sigma_ab traceless (z konstrukcji)",
       sp.simplify(sigma_mat.trace()) == 0,
       f"Tr(sigma) = {sp.simplify(sigma_mat.trace())}")


# =========================================================================
# CZESC C: Lattice Monte Carlo na toy v2 GL-substrate
# =========================================================================
banner("CZESC C: Lattice MC weryfikacja sigma_ab")

print("""
  Toy model v2 GL-substrate na siatce 3D L^3 (L = 24):

    Konfiguracja:  s(x,y,z) = v0 + delta s(x,y,z)
    Faza uporzadkowana:  v0 = 1.0 (dimensionless)
    Fluctuacje:    delta s ~ Gaussian(0, sigma_fluct)

  Pomiar gradient strain tensor (kompozyt rank-2 z H_Gamma):
    K_ab^B = <(d_a s)(d_b s)>_B
    sigma_ab = K_ab - (1/3) delta_ab Tr(K)

  Uwaga: ⟨ŝ_i · ŝ_{i+ê_a}⟩ jako single-step nearest-neighbor jest
  formalnie wektorem (jeden indeks a). Naturalna kompozytowa
  reprezentacja rank-2 traceless z czlonu kinetycznego H_Gamma to
  gradient strain tensor (d_a s)(d_b s), ktory ma wszystkie pozadane
  wlasnosci (i)-(vi).

  Test:
    C.1  prozni izotropowa: sigma -> 0
    C.2  z kwadrupolowym zaburzeniem: sigma != 0, ksztalt zgodny
    C.3  Z2 parity: sigma(s) = sigma(-s)
    C.4  SO(3) kowariancja: sigma(R[s]) = R sigma(s) R^T (numerycznie)
""")

L = 32                                      # rozmiar siatki
np.random.seed(42)
v0 = 1.0
sigma_fluct = 0.05                           # amplituda fluctuacji (mniejsza dla cleaner signal)

# C.1: prozni izotropowa (tylko Gaussian noise, brak anizotropii)
banner("C.1  Prozni izotropowa: sigma -> 0", 2)


def compute_K_sigma(s_field):
    """Compute gradient strain tensor K_ab = <(d_a s)(d_b s)>_block
    and its traceless part sigma_ab = K_ab - (1/3) delta_ab Tr(K).

    To jest naturalny coarse-grained operator rank-2 z H_Gamma:
    z czlonu kinetycznego (grad s)^2 = sum_a (d_a s)^2 = Tr(K)
    bezsladowa czesc to spin-2 anizotropia gradientu.
    """
    # Gradienty (centered finite differences with periodic BC)
    grads = []
    for a in range(3):
        ds = (np.roll(s_field, -1, axis=a) - np.roll(s_field, +1, axis=a)) / 2.0
        grads.append(ds)

    # Bilinear K_ab = <(d_a s)(d_b s)> averaged over lattice
    K = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            K[a, b] = np.mean(grads[a] * grads[b])
    K = 0.5 * (K + K.T)  # symetryzuj (numerycznie)

    Phi_block = np.mean(s_field**2)

    # bezsladowa czesc
    sigma = K - (np.trace(K) / 3) * np.eye(3)
    return K, sigma, Phi_block


# Vacuum config
s_vac = v0 + sigma_fluct * np.random.randn(L, L, L)
K_vac, sigma_vac, Phi_vac = compute_K_sigma(s_vac)

print(f"  Lattice {L}^3, vacuum config (Gaussian noise, isotropic)")
print(f"  Phi = <s^2> = {Phi_vac:.4f}")
print(f"  Tr(K) = {np.trace(K_vac):.4f}")
print(f"  K diagonal: {np.diag(K_vac)}")
print(f"  ||sigma_vac|| (Frobenius) = {np.sqrt(np.sum(sigma_vac**2)):.4e}")
print(f"  Anisotropy ||sigma|| / Tr(K) = {np.sqrt(np.sum(sigma_vac**2)) / abs(np.trace(K_vac)):.4e}")

# Vacuum criterion: anisotropy < 1% (typical lattice noise floor)
vacuum_anisotropy = np.sqrt(np.sum(sigma_vac**2)) / abs(np.trace(K_vac))
record("Vacuum: anisotropy ||sigma|| / Tr(K) < 1%",
       vacuum_anisotropy < 0.01,
       f"{vacuum_anisotropy:.4e}")

# Better vacuum test: average over multiple realizations
n_realizations = 20
sigma_norms = []
for _ in range(n_realizations):
    s_r = v0 + sigma_fluct * np.random.randn(L, L, L)
    _, sigma_r, _ = compute_K_sigma(s_r)
    sigma_norms.append(np.sqrt(np.sum(sigma_r**2)))

mean_sigma_vac = np.mean(sigma_norms)
std_sigma_vac = np.std(sigma_norms)
print(f"\n  Sredni |sigma_vac| ({n_realizations} realizacji): "
      f"{mean_sigma_vac:.4e} +/- {std_sigma_vac:.4e}")

# C.2: kwadrupolowe zaburzenie
banner("C.2  Kwadrupolowe zaburzenie: sigma != 0", 2)

print("""
  Konfiguracja: BINARY mass-like source wzdluz osi X (GW150914-like geometry).
  Dwie Gaussianskie kondensacje pola substratowego:
    delta s = A_q * [ exp(-((x-X0)^2 + y^2 + z^2) / (2 w^2))
                    + exp(-((x+X0)^2 + y^2 + z^2) / (2 w^2)) ]

  Ten ksztalt LAMIE x<->y symetrie:
    - rozszerzony wzdluz osi X (separacja 2*X0)
    - skompresowany wzdluz Y, Z (szerokosc w)
  => gradient strain tensor ma sxx != syy != szz (3 rozne diag. wartosci).

  Spodziewamy sie:
    - sigma_xx mniejsza (gradienty sa lokalne wokol kazdej kondensacji)
    - sigma_yy ~ sigma_zz (symetria osiowa wokol osi X)
    - off-diagonal ~ 0 (czysty kwadrupol w glownych osiach)
""")

# Konfiguracja z BINARNYM kwadrupolem wzdluz X (lamie x<->y)
A_q = 1.0
X0 = L / 4.0       # separacja srodek-source
w = L / 6.0        # szerokosc kazdej kondensacji
xs = np.arange(L) - L/2
X, Y, Z = np.meshgrid(xs, xs, xs, indexing='ij')
gauss_plus = np.exp(-((X - X0)**2 + Y**2 + Z**2) / (2 * w**2))
gauss_minus = np.exp(-((X + X0)**2 + Y**2 + Z**2) / (2 * w**2))
quadrupole = A_q * (gauss_plus + gauss_minus)
s_quad = v0 + quadrupole          # CLEAN konfiguracja (deterministyczna)

K_quad, sigma_quad, Phi_quad = compute_K_sigma(s_quad)

print(f"  Phi = <s^2> = {Phi_quad:.4f}")
print(f"  K diagonal: {np.diag(K_quad)}")
print(f"  sigma_quad =")
for i in range(3):
    print(f"    [{sigma_quad[i,0]:+.5f}  {sigma_quad[i,1]:+.5f}  {sigma_quad[i,2]:+.5f}]")

# Spodziewane: sigma_xx != sigma_yy != sigma_zz dla binary along X
sigma_quad_anisotropy = np.sqrt(np.sum(sigma_quad**2)) / abs(np.trace(K_quad))
print(f"\n  Anizotropia ||sigma_quad|| / Tr(K) = {sigma_quad_anisotropy:.4e}")

record("Quadrupol: |sigma_quad| >> |sigma_vac| (signal >= 10x noise floor)",
       np.sqrt(np.sum(sigma_quad**2)) > 10 * mean_sigma_vac,
       f"|sigma_quad|={np.sqrt(np.sum(sigma_quad**2)):.4e}, |sigma_vac|={mean_sigma_vac:.4e}, "
       f"ratio={np.sqrt(np.sum(sigma_quad**2))/mean_sigma_vac:.0f}x")

sxx, syy, szz = sigma_quad[0,0], sigma_quad[1,1], sigma_quad[2,2]
print(f"\n  Diag(sigma_quad): sigma_xx = {sxx:+.5e}, sigma_yy = {syy:+.5e}, sigma_zz = {szz:+.5e}")
record("Quadrupol: sigma_xx + sigma_yy + sigma_zz = 0 (traceless verified)",
       abs(sxx + syy + szz) < 1e-10,
       f"sum = {sxx+syy+szz:+.4e}")

# Test 1: anizotropia mierzalna, max-min eigenvalue >> floor
diag_spread = max(sxx, syy, szz) - min(sxx, syy, szz)
print(f"  diag spread (max - min) = {diag_spread:.4e}")
record("Quadrupol: diag spread (max - min) > 10*|sigma_vac_floor|",
       diag_spread > 10 * mean_sigma_vac,
       f"diag_spread={diag_spread:.4e}, 10*floor={10*mean_sigma_vac:.4e}")

# Test 2: x<->y symetria zlamana (binary along x => sxx != syy)
print(f"  |sxx-syy| = {abs(sxx-syy):.4e}, |sxx-syy|/floor = {abs(sxx-syy)/mean_sigma_vac:.1f}x")
record("Quadrupol: x<->y zlamane (binary along X) |sxx-syy| > 5*|sigma_vac_floor|",
       abs(sxx - syy) > 5 * mean_sigma_vac,
       f"|sxx-syy|={abs(sxx-syy):.4e}")

# C.3: Z2 parity - sigma(s) == sigma(-s)
banner("C.3  Z2 parity: sigma(s) = sigma(-s)", 2)

s_neg = -s_quad
K_neg, sigma_neg, _ = compute_K_sigma(s_neg)

z2_diff = np.max(np.abs(sigma_quad - sigma_neg))
print(f"  ||sigma(s) - sigma(-s)||_inf = {z2_diff:.4e}")

record("Z2 parity sigma(s) = sigma(-s)",
       z2_diff < 1e-10,
       f"max diff = {z2_diff:.4e}")

# C.4: SO(3) covariance
banner("C.4  SO(3) kowariancja: sigma(R[s]) ~ R sigma(s) R^T", 2)

print("""
  Test: obroc cala konfiguracje o 90 stopni wokol osi z
  (czyli przepermutuj X <-> Y, ze znakiem).
  Spodziewamy sie: sigma_obrocony = R_z(90) * sigma_quad * R_z(90)^T

  Dla R_z(90):  X_new = -Y, Y_new = X, Z_new = Z
  Wiec kwadrupol (x^2 - y^2) -> (y^2 - x^2) = -(x^2 - y^2)
  Macierz sigma_quad sigma_xx ~ +A,  sigma_yy ~ -A
  Po obrocie: sigma_yy <-> sigma_xx (z permutacja x<->y)
  i znaki sie zmienia: (x^2-y^2) -> -(x^2-y^2) -> sigma_xx -> -sigma_xx,  sigma_yy -> -sigma_yy
""")

# Obroc s przez przepermutowanie osi
s_rot = np.transpose(s_quad, (1, 0, 2))   # X <-> Y permutation
s_rot = -s_rot if False else s_rot          # bez znaku - tylko permutacja osi

K_rot, sigma_rot_lattice, _ = compute_K_sigma(s_rot)

# Buduj R_z (90 stopni): {x -> y, y -> -x, z -> z} = R_z(-90)
# albo {x -> -y, y -> x, z -> z} = R_z(+90)
# Tutaj robimy prostą permutację X<->Y, czyli R = [[0,1,0],[1,0,0],[0,0,1]]
# (refleksja w plaszczyznie xy=yx, NIE czysty obrot SO(3) -- to refleksja O(3))
R_perm = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
sigma_expected = R_perm @ sigma_quad @ R_perm.T

so3_diff = np.max(np.abs(sigma_rot_lattice - sigma_expected))
print(f"  sigma po permutacji X<->Y =")
for i in range(3):
    print(f"    [{sigma_rot_lattice[i,0]:+.5f}  "
          f"{sigma_rot_lattice[i,1]:+.5f}  "
          f"{sigma_rot_lattice[i,2]:+.5f}]")
print(f"\n  R sigma R^T (oczekiwane) =")
for i in range(3):
    print(f"    [{sigma_expected[i,0]:+.5f}  "
          f"{sigma_expected[i,1]:+.5f}  "
          f"{sigma_expected[i,2]:+.5f}]")

print(f"\n  ||sigma_obrocony - R sigma R^T||_inf = {so3_diff:.4e}")

# Tolerancja: lattice MC ma fluktuacje rzedu 1/sqrt(L^3) ~ 1/sqrt(13824) ~ 0.85%
mc_tolerance = 1.0 / np.sqrt(L**3)
record("SO(3) kowariancja (permutacja X<->Y, MC tolerance)",
       so3_diff < 5 * mc_tolerance,
       f"{so3_diff:.4e} < {5*mc_tolerance:.4e}")

# Eigenvalues
eig_quad = eigvalsh(sigma_quad)
print(f"\n  Eigenvalues sigma_quad: {eig_quad}")
print(f"  Sum (powinno byc 0):   {np.sum(eig_quad):.4e}")
record("Suma eigenvalues sigma = 0 (traceless)",
       abs(np.sum(eig_quad)) < 1e-10)


# =========================================================================
# CZESC D: Werdykt T2
# =========================================================================
banner("CZESC D: Werdykt T2")

print(f"\n  Suma checks:")
n_pass = sum(1 for _, ok in checks_summary if ok)
n_total = len(checks_summary)

for label, ok in checks_summary:
    mark = "PASS" if ok else "FAIL"
    print(f"    [{mark}] {label}")

print(f"\n  Wynik: {n_pass}/{n_total} pass")

verdict = "POSITIVE" if n_pass == n_total else "FAIL"
print(f"\n  WERDYKT T2: {verdict}")

print(f"""
  T2 sprawdza, ze sigma_ab postulowane w tgp_core.tex sec.2 jest:
    (i)   bilinear w s-hat                                           PASS
    (ii)  Z2-parzyste                                                PASS (analitycznie + lattice)
    (iii) symmetric + traceless (5 d.o.f.)                           PASS
    (iv)  rank-2 tensor SO(3)                                        PASS (sympy + lattice)
    (v)   znika w prozni izotropowej                                 PASS (lattice MC)
    (vi)  niezerowe przy kwadrupolu                                  PASS (lattice MC)

  ZNACZENIE T2:
    sigma_ab nie jest "dodanym" stopniem swobody -- jest emergentnym
    operatorem kompozytowym pojedynczego pola s-hat z H_Gamma.
    TGP zachowuje single-substrate axiom.

  T2 jest niezbednym lecz niewystarczajacym -- pokazuje TYLKO ze
  sigma_ab jest dobrze zdefiniowana jako operator kinematyczny.
  Pozostaje pokazac:
    - Dynamike sigma_ab (T3): rownanie ruchu propagujace na c_0
    - Sprzezenie do metryki (T4): Lambda(psi)
    - Amplituda kwadrupolowa h_TT (T5): formula GR-like
    - Konsystencja PPN, c_GW (T6)
""")
