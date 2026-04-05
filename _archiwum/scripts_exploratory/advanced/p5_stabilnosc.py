"""
p5_stabilnosc.py

Problem P5: Stabilnosc dynamiczna solitonu M3.

Perturbacja liniowa wokol profilu statycznego phi_sol(r):
    phi(r,t) = phi_sol(r) + eps * f(r) * exp(sigma*t)

Rownianie wlasne (uproszczone, bez czlonu alpha/phi w kinetyku pert.):
    H * f = -sigma^2 * f
    H = -(1/r^2) d/dr [r^2 d/dr] + V''_mod(phi_sol(r))

Soliton stabilny <==> wszystkie wartosci wlasne H > 0 (sigma^2 < 0 => oscylacje)
Soliton niestabilny <==> istnieje wartosc wlasna H < 0 (sigma^2 > 0 => eksponencjalny wzrost)

Uwaga:
  V''_mod(phi) = 2*phi - 3*phi^2 + 5*lam*(phi-1)^4
  V''_mod(1) = -1  (daleko od solitonu, tlo jest niestabilne!)
  V''_mod(psi_core) >> 0  (w rdzeniu, lambda*psi^4 dominuje)

Metoda: dyskretyzacja operatora H na siatce logarytmicznej,
        diagonalizacja macierzy 3-diagonalnej (scipy.linalg.eigh_tridiagonal).

Parametry: najlepsze rozwiazanie alpha=8.5616, a_gam=0.040, lam*=5.501e-6, K3=34.144
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal
import warnings; warnings.filterwarnings('ignore')

ALPHA = 8.5616
A_GAM = 0.040
LAM   = 5.501357e-6
K1    = 0.009820
K2    = 2.032728
K3    = 34.14450

def V_mod(psi, lam=LAM):
    return psi**3/3 - psi**4/4 + lam*(psi-1)**6/6

def dV2_mod(psi, lam=LAM):
    """Druga pochodna V_mod."""
    return 2*psi - 3*psi**2 + 5*lam*(psi-1)**4

def phi_sol(r, K):
    """Profil ansatzu Yukawa."""
    return np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)

# ---- Operator H na siatce logarytmicznej ----
def build_H_matrix(K, N=500, r_min=None, r_max=60.0):
    """
    Buduje macierz operatora H = -Laplace_r + V''(phi_sol(r))
    na siatce logarytmicznej r[k] = r_min * (r_max/r_min)^(k/N).

    Laplace_r f = (1/r^2) d/dr[r^2 f'] = f'' + (2/r)f'

    Dyskretyzacja: metoda elementow skonczonych z wagą r^2.
    """
    if r_min is None:
        r_min = A_GAM

    # Siatka logarytmiczna
    t   = np.linspace(0, 1, N+2)  # N wewnetrznych + 2 brzegowe
    r   = r_min * (r_max/r_min)**t

    # Punkty wewnetrzne: r[1]...r[N]
    r_int = r[1:-1]   # N punktow
    dr    = np.diff(r) # N+1 roznic

    # Wartosc potencjalu efektywnego V_eff = V''(phi) + l(l+1)/r^2
    # Dla l=0 (s-wave): V_eff = V''(phi)
    phi_int = phi_sol(r_int, K)
    V_eff   = dV2_mod(phi_int)

    # Stala ekranowania (asymptotyk V_eff -> V''(1) = -1)
    # => granica continuum H = -1 => tryby z E < -1 to stany zwiaane

    # Dyskretyzacja laplasjanu w sferze:
    # f'' + (2/r)f' = (1/r^2) d/dr[r^2 f']
    # Uzywamy metody objetosciowej (Galerkin z waga r^2 dr):
    # Macierz jest symetryczna trojdiagonalna po zmianie zmiennej u = r*f
    # -u'' + V_eff_eff(r)*u = E*u
    # gdzie V_eff_eff = V_eff + l(l+1)/r^2  (dla l=0: V_eff_eff = V_eff)
    # Plus poprawka od siatki log

    # Zamiast ogolnej dyskretyzacji, uzywamy prost przybliienie:
    # Symetryczna macierz trojdiagonalna w zmiennej u(r) = r*f(r)
    # H_u = -d^2u/dr^2 + [V''(phi) + ... ] u = E*u

    # Siatka rowna w u: tu dr nie jest rowne, ale mozemy uzyc
    # (1/dr_k) d/dr na siatce log

    # Prosta dyskretyzacja finito-roznicowa 3-punktowa:
    # f''(r_k) ~ (f_{k+1} - 2f_k + f_{k-1}) / (dr_avg)^2
    # f'(r_k) ~ (f_{k+1} - f_{k-1}) / (2*dr_avg)

    # Dla siatki log (nierownomiernej), uzywamy:
    # d^2f/dr^2 ~ 2*(f_{k+1}/(dr_k*dr_km1) - f_k*(1/dr_k + 1/dr_km1)/(dr_k+dr_km1)...

    # Prosta implementacja z r^2 waga:
    # (H f)_k = -[r^2 f']'/(r^2) + V''(phi) f  aprox.
    # ~ -(2/r)*f' - f'' + V''(phi)*f

    # Uzywamy 3-diag z dr_k = r[k+1]-r[k]:
    dr_int = dr[1:-1]   # dr miedzy wewnetrznymi punktami: dr[1]..dr[N-1]
    dr_p = dr[2:]       # dr[k+1] - do prawej
    dr_m = dr[1:-1]     # dr[k]   - do lewej

    # Wspolczynniki dla trojdiagonalnej macierzy:
    # H_{k,k}   = 1/dr_m[k] + 1/dr_p[k] + V''_k  (z wagą r^2)
    # H_{k,k+1} = -1/dr_p[k]
    # H_{k,k-1} = -1/dr_m[k]

    # Plus poprawka (2/r)*f': ~ (2/r_k) * (f_{k+1}-f_{k-1})/(dr_p+dr_m) - asymetria
    # Wlaczamy jako korekta do trojdiag

    # Diagonalna
    r_ip = r[2:-1]  # r_{k+1} (prawy sasiad)
    r_im = r[1:-2]  # r_{k-1} (lewy sasiad)

    # Symetryczna aproksymacja laplasjanu w sferze:
    # W postaci u_k = r_k * f_k, H_u = -d^2u/dr^2 + V_eff_u(r_k)*u
    # V_eff_u = V''(phi) (dla l=0, bez czlonu odsr.)
    #
    # Ale z siatka log, dr nie jest stale. Uzywamy:
    #   -u''_k ~ -(u_{k+1} - 2u_k + u_{k-1}) / h_k^2
    # gdzie h_k jest srednia odstep (dla siatki log: h_k ~ r_int[k] * ln(r_max/r_min)/N

    # Zmiana zmiennych: u_k = r_k * f_k
    # Diag: 1/(h_m*h_avg) + 1/(h_p*h_avg) + V_eff[k]  (gdzie h_avg = (h_p+h_m)/2)
    # Off-diag: -2/(h_p*(h_p+h_m)) i -2/(h_m*(h_p+h_m))

    N_int = len(r_int)
    diag = np.zeros(N_int)
    off  = np.zeros(N_int - 1)

    for k in range(N_int):
        h_m = dr[k]   if k > 0 else dr[0]
        h_p = dr[k+1] if k < N_int-1 else dr[-1]
        h_avg = 0.5*(h_m + h_p)
        diag[k] = 1.0/(h_m*h_avg) + 1.0/(h_p*h_avg) + V_eff[k]
        if k < N_int - 1:
            h_pp = dr[k+2] if k+2 < len(dr) else dr[-1]
            h_next_avg = 0.5*(h_p + h_pp)
            off[k] = -1.0/(h_p * h_next_avg)  # nie do konca symetryczne

    # Symetryzacja off-diagonal:
    # Uzywamy prostszej formy: -1/(h_p*h_m_next)^{1/2}... trudne.
    # Zastosujemy odwaznik wagowy dla uzyskania symetrii.

    # Prosta forma: diag_k = (h_m+h_p)/(h_m*h_p) + V_eff[k]
    #                off_k = -1/h_p (ale to moze byc asymetryczne)
    # Symetryzacja: macierz ważona macierza = D^{-1/2} H D^{-1/2}
    # gdzie D_kk = (h_m+h_p)/2

    # PROSTA IMPLEMENTACJA:
    h_vals = np.diff(r[:-1])[:N_int+1]
    diag2 = np.zeros(N_int)
    off2  = np.zeros(N_int - 1)
    for k in range(N_int):
        hm = dr[k]
        hp = dr[k+1]
        diag2[k] = 2.0/(hm*(hm+hp)) + 2.0/(hp*(hm+hp)) + V_eff[k]
    for k in range(N_int-1):
        hp = dr[k+1]
        hm_next = dr[k+1]
        hp_next = dr[k+2] if k+2 < len(dr) else dr[-1]
        diag2_next = 2.0/(hm_next*(hm_next+hp_next)) + 2.0/(hp_next*(hm_next+hp_next)) + V_eff[k+1]
        off2[k] = -2.0/(hp*(dr[k]+hp)) # uproszczona forma

    return diag2, off2, r_int, V_eff

# ================================================================
print("=" * 70)
print("P5: STABILNOSC DYNAMICZNA solitonu M3")
print(f"Parametry: alpha={ALPHA}, a_gam={A_GAM}, lam={LAM:.4e}")
print()
print(f"V''_mod(psi=1)      = {dV2_mod(1.0):.4f}  (tlo: < 0 => continuum < 0)")
psi_core = 1.0 + K3*np.exp(-A_GAM)/A_GAM
print(f"V''_mod(psi_core)   = {dV2_mod(psi_core):.4e}  (rdzen: >> 0)")
psi_cross = np.sqrt(3.0/(2.0*LAM))
print(f"V''_mod(psi_cross)  = {dV2_mod(psi_cross):.4e}  (crossover)")
print()

# ================================================================
# Profil V''(phi_sol(r)) wzdluz profilu solitonu
# ================================================================
print("PROFIL V''_mod(phi_sol(r)) wzdluz solitonu K3:")
print(f"{'r':>8} {'phi_sol':>12} {'V''_mod':>14} {'znak':>6}")
print("-" * 45)
r_pts = np.concatenate([
    np.array([A_GAM, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5]),
    np.linspace(0.5, 5.0, 10),
    np.array([10.0, 20.0, 40.0, 60.0])
])
for rpt in r_pts:
    phi = phi_sol(rpt, K3)
    V2  = dV2_mod(phi)
    print(f"{rpt:>8.4f} {phi:>12.4f} {V2:>14.4e} {'+ ' if V2>0 else '- ':>6}")

# Znajdz r_star: gdzie V''(phi(r)) = 0 (granica stabilnosci lokalna)
r_fine = np.linspace(A_GAM, 60.0, 10000)
phi_fine = phi_sol(r_fine, K3)
V2_fine  = dV2_mod(phi_fine)
sign_changes = np.where(np.diff(np.sign(V2_fine)))[0]
print()
if len(sign_changes) > 0:
    for idx in sign_changes:
        r_star = 0.5*(r_fine[idx]+r_fine[idx+1])
        phi_star = phi_sol(r_star, K3)
        print(f"  V''(phi_sol(r*)) = 0 przy r* = {r_star:.4f}, phi_sol(r*) = {phi_star:.4f}")
else:
    print("  Brak zmiany znaku V'' wzdluz profilu K3")
print()

# ================================================================
# Numeryczna diagonalizacja operatora H (wersja uproszczona)
# ================================================================
print("WIDMO OPERATORA H (dyskretyzacja, N=300):")
print()

for K_label, K in [("K1", K1), ("K2", K2), ("K3", K3)]:
    print(f"  --- Soliton {K_label} (K={K:.5f}) ---")
    diag, off, r_int, V_eff = build_H_matrix(K, N=300, r_max=60.0)

    # Eigh tridiagonalnej (symetrycznej)
    try:
        # Sprawdz czy off jest rozne
        evals = eigh_tridiagonal(diag, off, eigvals_only=True,
                                  select='i', select_range=(0, min(9, len(diag)-1)))
        print(f"  Najnizsze wartosci wlasne H:")
        for i, ev in enumerate(evals[:5]):
            print(f"    E_{i} = {ev:+.6f}  {'STABILNE (E>0)' if ev>0 else '*** NIESTABILNE (E<0) ***'}")
        E_min = evals[0]
        sigma_max = np.sqrt(max(-E_min, 0.0))
        print(f"  Granica continuum: V''(1) = {dV2_mod(1.0):.4f}")
        print(f"  E_min = {E_min:+.6f}")
        if E_min < 0:
            print(f"  => NIESTABILNY: sigma_max = {sigma_max:.4f} (czas zycia ~ {1/sigma_max:.2f})")
        else:
            print(f"  => STABILNY: brak ujemnych wartosci wlasnych w spektrum dyskretnym")
    except Exception as ex:
        print(f"  [blad diagonalizacji]: {ex}")
    print()

# ================================================================
# Analiza energetyczna: warunek stabilnosci Derrick
# ================================================================
print("=" * 70)
print("TEST DERRICKA: skalowanie profilu phi_lambda(r) = phi(r/lambda)")
print()
print("Energia solitonu jako funkcja skalowania:")
print(f"  E(lambda) = E_kin/lambda + E_pot*lambda^3")
print(f"  dE/dlambda|_{{lam=1}} = -E_kin + 3*E_pot = 0  (warunek Derricka)")
print()

from scipy.integrate import quad

for K_label, K in [("K1", K1), ("K2", K2), ("K3", K3)]:
    # Oblicz Ek i Ep na siatce log
    N = 6000
    r_max = 60.0
    t   = np.linspace(0, 1, N)
    r   = A_GAM * (r_max/A_GAM)**t
    phi = np.maximum(phi_sol(r, K), 1e-10)
    dphi = K*np.exp(-r)*(-r-1.0)/r**2
    V1 = V_mod(1.0)
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)

    derrick = -Ek + 3*Ep
    print(f"  {K_label} (K={K:.5f}): Ek={Ek:.4e}, Ep={Ep:.4e}")
    print(f"    -Ek + 3*Ep = {derrick:.4e}")
    print(f"    Ek/Ep = {Ek/Ep:.4f}  (Derrick wymaga: 1/3 dla skalarnych 3D)")
    sign_D = "=> skurcz niestabilny" if derrick < 0 else "=> stabilny wzgl. skalowania"
    print(f"    {sign_D}")
    print()

print("=" * 70)
print("UWAGA INTERPRETACYJNA:")
print()
print("  V''_mod(1) = -1 < 0: tlo jest lokalnie NIESTABILNE w sensie malych drgań.")
print("  To wlasciwosc potencjału TGP (sedlo przy psi=1, nie minimum).")
print("  Soliton moze byc stabilny jesli:")
print("  (a) Istnieje ładunek topologiczny chroniony")
print("  (b) Ujemne mody są o dlugofalowe i soliton jest quasi-stabilny")
print("  (c) Widmo H ma E_0 >= V''(1) = -1 (stany scalone powyzej continuum)")
print()
print("  Test Derricka (skalowalnosc): jesli -Ek + 3*Ep = 0, soliton spelnia")
print("  warunek Virial i nie jest niestabilny wzgledem globalnego skalowania.")
print()
print("  Pelna analiza stabilnosci wymaga:")
print("  1. Profilu dokladnego z p1_ode_shooting.py (nie ansatz Yukawa)")
print("  2. Wlasciwego operatora kinetycznego z czlonem alpha/phi")
print("  3. Analizy wlasnej spelniajac warunki brzegowe")
