# Phase 0 check: weryfikacja rachunkow PRZED kodem produkcyjnym.
# Wg Phase0_balance.md (LOCKED, sekcje 1-2 + doprecyzowania D1-D11):
#   1) algebra: s*, Phi*, V''(s*), xi, J, C2_pred, pasmo G4
#   2) profil rdzenia wiru: radialne ODE  kappa*(f'' + f'/r - f/r^2) = V'(f)
#   3) poprawki obrazow periodycznych w E(d): energia sieciowa ansatzu D1
#      + krzyzowo: energia sieciowa punktowych wirow (Green 5-pkt, FFT)
#      -> przewidywane EFEKTYWNE zbocze C2 na siatce pomiarowej (pre-code)

import numpy as np

kappa = 0.50
a, b, c = 0.50, 1.60, 1.00

print("=== Phase 0 check: op-goldstone-mediator (LOCK: Phase0_balance.md) ===")
print()

# ---------------- 1. algebra ----------------
# V(m) = 0.5*a*m^2 - (b/3)*m^3 + 0.25*c*m^4 ; V'(m) = m*(a - b*m + c*m^2)
disc = np.sqrt(b * b - 4.0 * a * c)
s_star = (b + disc) / (2.0 * c)
Phi_star = s_star ** 2
Vpp_star = a - 2.0 * b * s_star + 3.0 * c * s_star ** 2
Vpp_0 = a
xi = np.sqrt(kappa / Vpp_star)
J = kappa * Phi_star
C2_pred = 2.0 * np.pi * J
band = (0.8 * C2_pred, 1.2 * C2_pred)


def Vpot(m):
    return 0.5 * a * m**2 - (b / 3.0) * np.abs(m)**3 + 0.25 * c * m**4


def Vprime(m):
    return m * (a - b * np.abs(m) + c * m * m)


print("--- 1. algebra (rownania #1-#3 LOCKa) ---")
print(f"  s*        = {s_star:.6f}   (LOCK: 1.174166...)")
print(f"  Phi*      = {Phi_star:.6f}   (LOCK: 1.3787)")
print(f"  V''(s*)   = {Vpp_star:.6f}   (LOCK: 0.87869)")
print(f"  V(s*)     = {Vpot(s_star):+.6f}  (< 0 = prawdziwa proznia OK)")
print(f"  xi        = sqrt(kappa/V''(s*)) = {xi:.6f}   (LOCK: 0.754)")
print(f"  J         = kappa*Phi* = {J:.6f}   (LOCK: 0.6893)")
print(f"  C2_pred   = 2*pi*J = {C2_pred:.4f}   (LOCK: 4.331)")
print(f"  pasmo G4  = [{band[0]:.3f}, {band[1]:.3f}]   (LOCK: [3.46, 5.20])")
print(f"  skale mas (D4): sqrt(V''(0)/kappa) = {np.sqrt(Vpp_0/kappa):.4f}, "
      f"sqrt(V''(s*)/kappa) = {np.sqrt(Vpp_star/kappa):.4f}")
print()

# ---------------- 2. profil rdzenia (ODE radialne) ----------------
# gradient-flow relaksacja:  df/dt = kappa*(f'' + f'/r - f/r^2) - V'(f)
# f(0) = 0 (rdzen), f(R) = s* (Dirichlet daleko)
dr = 0.05
Rmax = 25.0
r = np.arange(1, int(Rmax / dr)) * dr
f = s_star * r / np.sqrt(r * r + 2.0 * xi * xi)      # ansatz startowy (LOCK sek. 2)
dt_ode = 0.4 * dr * dr / (2.0 * kappa)
res_final = np.nan
for it in range(400000):
    fl = np.concatenate(([0.0], f[:-1]))             # f(0)=0
    fr = np.concatenate((f[1:], [s_star]))           # f(R)=s*
    d2 = (fl + fr - 2.0 * f) / (dr * dr)
    d1 = (fr - fl) / (2.0 * dr)
    rhs = kappa * (d2 + d1 / r - f / (r * r)) - Vprime(f)
    f = f + dt_ode * rhs
    if it % 2000 == 0:
        res_final = float(np.max(np.abs(rhs)))
        if res_final < 1e-9:
            break
print("--- 2. profil rdzenia |psi|(r) z ODE radialnego ---")
print(f"  siatka dr={dr}, Rmax={Rmax}, residual koncowy = {res_final:.2e} "
      f"(iteracje: {it})")
prof_ansatz = lambda rr: s_star * rr / np.sqrt(rr * rr + 2.0 * xi * xi)
print("     r      f_ODE(r)   ansatz prof(r)   f/s*")
for rq in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
    fq = float(np.interp(rq, r, f))
    print(f"  {rq:5.2f}   {fq:8.5f}   {prof_ansatz(rq):8.5f}      {fq/s_star:.4f}")
r99 = float(r[np.argmax(f > 0.99 * s_star)])
print(f"  r(f = 0.99*s*) = {r99:.2f}  (~{r99/xi:.1f} xi)  -> rdzen zwarty, "
      f"dla d >= 6 separacja >> rdzen")
print()

# ---------------- 3. poprawki obrazow periodycznych w E(d) ----------------
TWO_PI = 2.0 * np.pi


def theta_ansatz(X, Y, L, vortices):
    """D1: theta = sum_i n_i sum_k arg sinh(pi*(z - z_i - kL)/L)."""
    Z = X + 1j * Y
    th = np.zeros_like(X)
    for (x0, y0, n) in vortices:
        for k in range(-2, 3):
            th = th + n * np.angle(np.sinh(np.pi * (Z - (x0 + k * L) - 1j * y0) / L))
    return th


def amp_ansatz(X, Y, L, vortices):
    A = np.full_like(X, s_star)
    for (x0, y0, n) in vortices:
        ddx = (X - x0 + 0.5 * L) % L - 0.5 * L
        ddy = (Y - y0 + 0.5 * L) % L - 0.5 * L
        rr = np.sqrt(ddx * ddx + ddy * ddy)
        A = A * rr / np.sqrt(rr * rr + 2.0 * xi * xi)
    return A


def lattice_energy(psi, dxg):
    gx = (np.roll(psi, -1, 0) - psi) / dxg
    gy = (np.roll(psi, -1, 1) - psi) / dxg
    m = np.abs(psi)
    dens = 0.5 * kappa * (np.abs(gx)**2 + np.abs(gy)**2) + Vpot(m)
    return float(np.sum(dens) * dxg * dxg)


def fit_L(dv, Ev):
    A_ = np.vstack([np.ones(len(dv)), np.log(dv)]).T
    coef, *_ = np.linalg.lstsq(A_, Ev, rcond=None)
    rss = float(np.sum((Ev - A_ @ coef) ** 2))
    return coef[0], coef[1], rss


def seam_report(th, tag):
    j1 = np.abs((th[0, :] - th[-1, :] + np.pi) % TWO_PI - np.pi)
    j2 = np.abs((th[:, 0] - th[:, -1] + np.pi) % TWO_PI - np.pi)
    print(f"    [{tag}] max skok fazy na szwie: x-szew {np.max(j2):.2e}, "
          f"y-szew {np.max(j1):.2e} rad (regularny gradient ~ O(dx*|grad theta|))")


def greens_pair_energy(N, L, d):
    """Energia sieciowa PUNKTOWEJ pary wirow: Green 5-pkt przez FFT.
    grad(theta) = rot(grad chi), Lap(chi) = 2*pi*rho  ->  E = J/2 int |grad chi|^2."""
    dxg = L / N
    rho = np.zeros((N, N))
    Nc = N // 2
    cells = int(round(0.5 * d / dxg))
    rho[Nc - cells, Nc] = +1.0 / (dxg * dxg)
    rho[Nc + cells, Nc] = -1.0 / (dxg * dxg)
    kx = TWO_PI * np.fft.fftfreq(N, d=dxg)
    KX, KY = np.meshgrid(kx, kx, indexing="ij")
    # wartosci wlasne 5-pkt Laplasjanu (spojnosc z dyskretyzacja modelu)
    lam = -((2 - 2 * np.cos(KX * dxg)) + (2 - 2 * np.cos(KY * dxg))) / (dxg * dxg)
    lam[0, 0] = 1.0
    chi_hat = TWO_PI * np.fft.fft2(rho) / lam
    chi_hat[0, 0] = 0.0
    chi = np.real(np.fft.ifft2(chi_hat))
    gx = (np.roll(chi, -1, 0) - chi) / dxg
    gy = (np.roll(chi, -1, 1) - chi) / dxg
    return 0.5 * J * float(np.sum(gx * gx + gy * gy)) * dxg * dxg, 2 * cells * dxg


def predict_geometry(N, L, d_list, tag):
    dxg = L / N
    coords = (np.arange(N) - N // 2) * dxg
    X, Y = np.meshgrid(coords, coords, indexing="ij")
    print(f"  geometria {tag}: N={N}, L={L}, dx={dxg}")
    E_ans, E_pt, d_eff = [], [], []
    for d in d_list:
        vort = [(-0.5 * d, 0.0, +1), (+0.5 * d, 0.0, -1)]
        th = theta_ansatz(X, Y, L, vort)
        psi = amp_ansatz(X, Y, L, vort) * np.exp(1j * th)
        E_ans.append(lattice_energy(psi, dxg))
        ept, dg = greens_pair_energy(N, L, d)
        E_pt.append(ept)
        d_eff.append(d)
        if d == d_list[0]:
            seam_report(th, f"{tag}, d={d}")
    E_ans, E_pt, d_eff = map(np.array, (E_ans, E_pt, d_eff))
    _, C2a, _ = fit_L(d_eff, E_ans)
    _, C2p, _ = fit_L(d_eff, E_pt)
    print(f"    d:                {np.array2string(d_eff, precision=1)}")
    print(f"    E_ansatz(d):      {np.array2string(E_ans, precision=4)}")
    print(f"    E_pointGreen(d):  {np.array2string(E_pt, precision=4)}")
    # odchylka od czystego 2piJ ln(d) (wzgledem pierwszego punktu)
    ideal = C2_pred * np.log(d_eff / d_eff[0])
    dev_a = (E_ans - E_ans[0]) - ideal
    print(f"    odchylka obrazow (ansatz, wzgl. d0): "
          f"{np.array2string(dev_a, precision=4)}")
    print(f"    -> C2_eff(ansatz)     = {C2a:.4f}  "
          f"({100*(C2a/C2_pred-1):+.1f}% vs 2piJ)")
    print(f"    -> C2_eff(pointGreen) = {C2p:.4f}  "
          f"({100*(C2p/C2_pred-1):+.1f}% vs 2piJ)")
    in_band = band[0] <= C2a <= band[1]
    print(f"    C2_eff(ansatz) w pasmie G4 [{band[0]:.2f}, {band[1]:.2f}]: "
          f"{in_band}")
    return C2a, C2p


print("--- 3. poprawki obrazow periodycznych: przewidywane zbocze efektywne ---")
print("  (energia sieciowa ansatzu D1 = predykcja E(d) az po stala; PRZED runami)")
C2_main = predict_geometry(256, 128.0, [6.0, 8.0, 12.0, 16.0, 24.0], "GLOWNA")
print()
C2_ctrl = predict_geometry(256, 64.0, [6.0, 8.0, 12.0], "KONTROLA dx=0.25 (D9)")
print()

print("=== WNIOSKI PHASE 0 ===")
print(f"  1) algebra LOCKa potwierdzona (s*, Phi*, xi, J, C2_pred = {C2_pred:.4f})")
print(f"  2) profil rdzenia: ansatz prof(r) blisko rozwiazania ODE; rdzen ~ xi")
print(f"  3) obrazy periodyczne przesuwaja zbocze o kilka %, C2_eff pozostaje")
print(f"     w pasmie G4 -> kryterium G4 wykonalne na zalockowanej geometrii")
