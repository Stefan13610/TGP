# Phase 0 pre-code verification of analytic claims locked in Phase0_balance.md
# NO experiment dynamics here - only checks of the locked analytics.
import numpy as np

a, b, c = 0.50, 1.60, 1.00
kappa   = 0.50
dx, dt  = 0.5, 0.02
w       = 1.50   # seed width

def V(s):
    s = np.asarray(s, dtype=float)
    return 0.5*a*s**2 - (b/3.0)*np.abs(s)**3 + 0.25*c*s**4

def Vp(s):
    s = np.asarray(s, dtype=float)
    return s*(a - b*np.abs(s) + c*s**2)

def Vpp(s):
    s = np.asarray(s, dtype=float)
    return a - 2.0*b*np.abs(s) + 3.0*c*s**2

print("=== Phase 0: weryfikacja rachunkow analitycznych (pre-code) ===")
print(f"potencjal: a={a}, b={b}, c={c}, kappa={kappa}, dx={dx}, dt={dt}")
print()

# 1. extrema of V
disc = b*b - 4*a*c
s_bar  = (b - np.sqrt(disc))/(2*c)
s_star = (b + np.sqrt(disc))/(2*c)
print("[1] Ekstrema V (pierwiastki a - b|s| + c s^2 = 0):")
print(f"    b^2-4ac = {disc:.4f}  (LOCK: 0.56)  -> {'OK' if abs(disc-0.56)<1e-12 else 'MISMATCH'}")
print(f"    bariera  |s_bar| = {s_bar:.4f}  Phi_bar = {s_bar**2:.4f}  (LOCK: ~0.426 / ~0.181)")
print(f"    prozniaN |s*|    = {s_star:.4f}  Phi*    = {s_star**2:.4f}  (LOCK: ~1.174 / ~1.378)")
print(f"    residua V'(s_bar), V'(s*): {Vp(s_bar):.2e}, {Vp(s_star):.2e}")
ok1 = abs(s_bar-0.426)<2e-3 and abs(s_star-1.174)<2e-3 and abs(Vp(s_bar))<1e-12 and abs(Vp(s_star))<1e-12
print(f"    -> {'PASS' if ok1 else 'FAIL'}")
print()

# 2. energies and curvatures
print("[2] Energie i krzywizny:")
print(f"    V(0)      = {V(0.0):.6f}")
print(f"    V(s_bar)  = {V(s_bar):.6f}   (wysokosc bariery)")
print(f"    V(s*)     = {V(s_star):.6f}  (LOCK: ~ -0.044, ma byc < 0)")
print(f"    V''(0)    = {Vpp(0.0):.4f}     (>0 => s=0 metastabilne)")
print(f"    V''(s*)   = {Vpp(s_star):.4f}   (>0 => s* stabilne)")
print(f"    V''(s_bar)= {Vpp(s_bar):.4f}   (<0 => bariera)")
s_dip = b/(3*c)
print(f"    min V''  przy |s|=b/3c={s_dip:.4f}: V''={Vpp(s_dip):.4f}")
ok2 = V(s_star) < 0 and Vpp(0.0)>0 and Vpp(s_star)>0 and Vpp(s_bar)<0 and abs(V(s_star)+0.044)<1e-3
print(f"    G6 pre-prediction: V''(s*)={Vpp(s_star):.3f} > V''(0)={Vpp(0.0):.3f}  "
      f"-> maly bump w CORE drozszy niz w FAR: {'potwierdzone' if Vpp(s_star)>Vpp(0.0) else 'NIE'}")
print(f"    -> {'PASS' if ok2 else 'FAIL'}")
print()

# 3. wall width vs dx (pinning risk locked in G5)
ell = np.sqrt(2*kappa/Vpp(s_star))
print("[3] Szerokosc scianki domenowej vs siatka:")
print(f"    ell ~ sqrt(2*kappa/V''(s*)) = {ell:.4f}  = {ell/dx:.2f} * dx")
print(f"    -> ryzyko pinningu realne (~2*dx), kontrola z G5 (dx=0.25 -> {ell/0.25:.2f} wezly) zasadna")
print()

# 4. explicit Euler stability
dt_diff = dx*dx/(4*kappa)
ss = np.linspace(-2.0, 2.0, 20001)
Vpp_max_scan = np.max(np.abs(Vpp(ss)))
dt_react = 2.0/Vpp_max_scan
print("[4] Stabilnosc jawnego Eulera:")
print(f"    dyfuzja:  dt < dx^2/(4*kappa) = {dt_diff:.4f}   (LOCK: 0.125), dt={dt} -> {'OK' if dt<dt_diff else 'FAIL'}")
print(f"    reakcja:  max|V''| na [-2,2] = {Vpp_max_scan:.3f} -> dt < 2/max|V''| = {dt_react:.4f}, dt={dt} -> {'OK' if dt<dt_react else 'FAIL'}")
dt_ctrl_diff = 0.25**2/(4*kappa)
print(f"    kontrola G5 (dx=0.25, dt=0.005): dt < {dt_ctrl_diff:.5f} -> {'OK' if 0.005<dt_ctrl_diff else 'FAIL'}")
print()

# 5. cluster superposition factor (basis for the superposition-trivial flag in G4)
print("[5] Wspolczynnik superpozycji klastra K=7 (heks pierscien + centrum), w=1.5:")
for D in (2, 3, 4):
    f_center = 1.0 + 6.0*np.exp(-D*D/(2*w*w))
    print(f"    D={D}: s_max(0) ~ {f_center:.3f} * A0   (6 sasiadow w odl. D)")
print("    -> D=2 silnie superpozycyjny (~3.5*A0), D=4 prawie addytywnie czysty (~1.17*A0);")
print("       flaga superposition-trivial z G4 ma operacyjny sens.")
print()

# 6. bump probe sanity (locked delta_bump: A=0.10, w=1.0) - potential-only per-site estimate
A_bump = 0.10
print("[6] Sanity delta_bump (A=0.10, w=1.0), oszacowanie czysto potencjalowe (bez czlonu grad):")
print(f"    FAR  (s=0):   dV ~ 0.5*V''(0)*A^2  = {0.5*Vpp(0.0)*A_bump**2:.6f} > 0")
print(f"    CORE (s=s*):  dV ~ 0.5*V''(s*)*A^2 = {0.5*Vpp(s_star)*A_bump**2:.6f} > 0  (drozszy niz FAR - zgodnie z predykcja G6)")
print(f"    FRONT: bump przesuwa interfejs w strone V(s*)={V(s_star):.4f}<0 -> jedyna lokalizacja z szansa na dE<=0")
print()

all_ok = ok1 and ok2 and dt < dt_diff and dt < dt_react
print(f"=== WYNIK FAZY 0: {'PASS - rachunki LOCK potwierdzone, model gotowy do Phase 1' if all_ok else 'FAIL - rozbieznosc z LOCK'} ===")
