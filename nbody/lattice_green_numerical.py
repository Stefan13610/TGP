"""
Numerical comparison: continuum soliton vs lattice-regularized soliton.

Computes the 3D cubic lattice Green's function numerically and
shows that h(r) = r * psi(r) remains monotonically increasing.
"""
import numpy as np

# ============================================================
# 1. Lattice Green's function via Fourier integral
# ============================================================

def lattice_green_1d(n_max=200, N_grid=256):
    """
    Compute lattice Green's function G(n,0,0) for n = 0..n_max
    on a 3D simple cubic lattice with unit spacing.
    Uses numerical integration over the Brillouin zone with
    proper zero-mode regularization.
    """
    # Shift grid by half-step to avoid k=0 zero mode
    dk = 2 * np.pi / N_grid
    k = -np.pi + dk * (np.arange(N_grid) + 0.5)
    k1, k2, k3 = np.meshgrid(k, k, k, indexing='ij')

    # Lattice Laplacian eigenvalue: lambda(k) = 6 - 2(cos k1 + cos k2 + cos k3)
    lam = 6.0 - 2*np.cos(k1) - 2*np.cos(k2) - 2*np.cos(k3)

    G_values = np.zeros(n_max + 1)
    for n in range(n_max + 1):
        integrand = np.cos(n * k1) / lam
        G_values[n] = np.sum(integrand) * dk**3 / (2*np.pi)**3

    return G_values

print("Computing lattice Green's function (N=128)...")
n_max = 100
G = lattice_green_1d(n_max=n_max, N_grid=128)

print(f"G(0) = {G[0]:.6f}")
print(f"G(1) = {G[1]:.6f}")
print(f"G(5) = {G[5]:.6f}")
print(f"1/(4pi*5) = {1/(4*np.pi*5):.6f}")

# Normalization: for large n, G(n) ~ 1/(4*pi*n)
n_ref = 50
expected = 1.0 / (4*np.pi*n_ref)
C_norm = G[n_ref] / expected
print(f"\nNormalization: G({n_ref})/[1/(4pi*{n_ref})] = {C_norm:.6f}")

# Normalized: G_n such that G_n ~ 1/(4pi*r) at large r
# G_norm = G / C_norm  (but C_norm should be ~1 for correct integral)
# Let's just use raw G values which should already be ~ 1/(4pi*n)

# Watson's integral value for G(0):
# G_lattice(0,0,0) = W_3 / (6) where W_3 = 1.516... (Watson triple integral)
# Actually G(0) = (1/(2pi)^3) int dk / lambda(k)
# Known value: G(0) = 0.253... for 3D simple cubic
watson_G0 = 0.2527
print(f"\nWatson value G(0) = {watson_G0:.4f}")
print(f"Computed     G(0) = {G[0]:.4f}")
print(f"Agreement: {abs(G[0] - watson_G0)/watson_G0*100:.1f}%")

# ============================================================
# 2. Build soliton for r_s/a = 100 (physical-ish regime)
# ============================================================
print("\n" + "="*55)
print("SOLITON COMPARISON: r_s = 100a")
print("="*55)

r_s = 100.0  # r_s in lattice units
C_sol = 3 * r_s

n_arr = np.arange(1, n_max + 1, dtype=float)

# Continuum: u(n) = 1 + C_sol / n
u_cont = 1.0 + C_sol / n_arr
psi_cont = u_cont ** (1.0/3)
h_cont = n_arr * psi_cont

# Lattice: u(n) = 1 + C_sol * G(n) / G_asymptotic_norm
# G(n) ~ 1/(4pi*n) at large n, so we need u_latt ~ 1 + C_sol/n at large n
# => multiply G by 4*pi: 4pi*G(n) ~ 1/n
G_scaled = 4 * np.pi * G[1:n_max+1]
u_latt = 1.0 + C_sol * G_scaled
psi_latt = u_latt ** (1.0/3)
h_latt = n_arr * psi_latt

# Check monotonicity
dh_cont = np.diff(h_cont)
dh_latt = np.diff(h_latt)
print(f"h_continuum monotonic: {np.all(dh_cont > 0)}")
print(f"h_lattice   monotonic: {np.all(dh_latt > 0)}")

# Compare at key points
print(f"\n{'n':>5s} {'psi_cont':>10s} {'psi_latt':>10s} {'h_cont':>10s} {'h_latt':>10s} {'rel_diff':>10s}")
for n in [1, 2, 5, 10, 50, 100]:
    if n <= n_max:
        i = n - 1
        pc = psi_cont[i]
        pl = psi_latt[i]
        hc = h_cont[i]
        hl = h_latt[i]
        rd = abs(pc - pl)/pc * 100
        print(f"{n:5d} {pc:10.4f} {pl:10.4f} {hc:10.4f} {hl:10.4f} {rd:9.2f}%")

# ============================================================
# 3. r_s/a = 1 (extreme regime: BH size = lattice spacing)
# ============================================================
print("\n" + "="*55)
print("EXTREME REGIME: r_s = 1a (BH = lattice spacing)")
print("="*55)

r_s_ext = 1.0
C_ext = 3 * r_s_ext
u_latt_ext = 1.0 + C_ext * G_scaled
psi_latt_ext = u_latt_ext ** (1.0/3)
h_latt_ext = n_arr * psi_latt_ext

dh_ext = np.diff(h_latt_ext)
print(f"h_lattice monotonic: {np.all(dh_ext > 0)}")
if not np.all(dh_ext > 0):
    bad = np.where(dh_ext <= 0)[0]
    print(f"  Violations at n = {bad + 1}")
    for b in bad[:5]:
        print(f"    h({b+1}) = {h_latt_ext[b]:.6f}, h({b+2}) = {h_latt_ext[b+1]:.6f}")
else:
    print("  => NO PHOTON SPHERE even at r_s = a")

print(f"\npsi_latt(a) = {psi_latt_ext[0]:.6f}  (continuum: {(1+3)**(1/3):.6f})")
print(f"h_latt(a)   = {h_latt_ext[0]:.6f}  (continuum: {1*(1+3)**(1/3):.6f})")

# ============================================================
# 4. Analytical bound
# ============================================================
print("\n" + "="*55)
print("FUNDAMENTAL PROTECTION: the 1/3 factor")
print("="*55)
print("""
Photon sphere requires h'(r*) = 0, i.e. psi'/psi = -1/r.

Since psi = u^{1/3}:
  psi'/psi = (1/3) * u'/u

For u satisfying discrete Laplace eq, u is positive and
subharmonic: u'/u <= 1/r at the steepest point.

Therefore: |psi'/psi| = (1/3)|u'/u| <= 1/(3r) < 1/r

The factor 1/3 from the cube root PREVENTS photon sphere
for ANY solution of the (discrete or continuum) Laplace equation.

This is a STRUCTURAL result:
- Independent of lattice spacing a
- Independent of r_s/a ratio
- Independent of lattice geometry (cubic, FCC, etc.)
- Depends ONLY on the exponent 1/3 in psi = u^{1/3}

The exponent 1/3 comes from the nonlinear kinetic operator
with alpha = 2 (the TGP value). For alpha != 2, the exponent
would be different and a photon sphere might be possible.
""")

# What exponent would be needed?
print("What exponent p in psi = u^p would allow a photon sphere?")
print("  Need: p * |u'/u| >= 1/r near center")
print("  Since u'/u -> -1/r for harmonic u:")
print("  Need: p >= 1")
print(f"  TGP has: p = 1/3 (from alpha = 2)")
print(f"  GR-like: p = 1 would be borderline")
print(f"  p > 1 would guarantee photon sphere")
print()
print("The alpha-p relation: p = 1/(alpha+1)")
print(f"  alpha=2 (TGP):   p = 1/3  => NO photon sphere")
print(f"  alpha=0 (linear): p = 1    => BORDERLINE")
print(f"  alpha=-1/2:       p = 2    => photon sphere EXISTS")
