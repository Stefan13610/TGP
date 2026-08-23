# Phase 2: wir jako tlo stacjonarne (K2) — G3, G4 + G5-demo
# (LOCK: Phase0_balance.md, sekcje 2-4)

import numpy as np

kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt_flow = 0.02
dt_w = 0.05
N, L = 256, 128.0
dx = L / N
S_GUARD = 10.0
H_TOL = 1e-10
SEED = 20260704
TWO_PI = 2.0 * np.pi

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)
Phi_star = s_star ** 2
Vpp_star = a - 2.0 * b_par * s_star + 3.0 * c * s_star ** 2
xi2 = kappa / Vpp_star

coords = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
Z = X + 1j * Y
V1 = (-32.0, -32.0)
V2 = (+32.0, +32.0)


def gfun(m):
    return a - b_par * np.abs(m) + c * m * m       # V'(m)/m


def Vpp(m):
    return a - 2.0 * b_par * np.abs(m) + 3.0 * c * m * m


def Vpot(m):
    return 0.5 * a * m**2 - (b_par / 3.0) * m**3 + 0.25 * c * m**4


def lap(u):
    return (np.roll(u, 1, 0) + np.roll(u, -1, 0)
            + np.roll(u, 1, 1) + np.roll(u, -1, 1) - 4.0 * u) / (dx * dx)


def rhs(psi):
    return kappa * lap(psi) - psi * gfun(np.abs(psi))


def energy_grad(psi):
    gx = (np.roll(psi, -1, 0) - psi) / dx
    gy = (np.roll(psi, -1, 1) - psi) / dx
    return float(np.sum(0.5 * kappa * (np.abs(gx)**2 + np.abs(gy)**2)
                        + Vpot(np.abs(psi))) * dx * dx)


def theta_h_pair(z1, z2):
    """Para wspolliniowa w y (ta sama wsp. y): konstrukcja D1 (dokladna:
    asymptotyki x-szwu kasuja sie dla rownych y — dowod B1)."""
    th = np.zeros_like(X)
    for (x0, y0), n in ((z1, +1), (z2, -1)):
        for k in range(-2, 3):
            th = th + n * np.angle(np.sinh(np.pi * (Z - (x0 + k * L) - 1j * y0) / L))
    return th


def theta_v_pair(z1, z2):
    """Para wspolliniowa w x: konstrukcja transponowana (x <-> y);
    orientacja windingu odwrocona przez transpozycje -> znak -n."""
    Wt = Y + 1j * X
    th = np.zeros_like(X)
    for (x0, y0), n in ((z1, +1), (z2, -1)):
        for k in range(-2, 3):
            th = th - n * np.angle(np.sinh(np.pi * (Wt - (y0 + k * L) - 1j * x0) / L))
    return th


def make_pair():
    # POPRAWKA IMPLEMENTACYJNA (errata; intencja LOCKa D1): dla pary na
    # przekatnej goly iloczyn sinh zostawia skok fazy pi na szwie x
    # (zmierzone: 272 linki, max = pi -> relaksacja tworzy sztuczna ciemna
    # sciane). Kompozycja dokladna: (+1 A, -1 B) = (+1 A, -1 C) + (+1 C, -1 B),
    # C = (x_B, y_A): para pozioma (D1) + para pionowa (transponowana D1).
    A_, B_ = V1, V2
    C_ = (B_[0], A_[1])
    th = theta_h_pair(A_, C_) + theta_v_pair(C_, B_)
    Amp = np.full_like(X, s_star)
    for (x0, y0) in (V1, V2):
        ddx = (X - x0 + 0.5 * L) % L - 0.5 * L
        ddy = (Y - y0 + 0.5 * L) % L - 0.5 * L
        rr = np.sqrt(ddx * ddx + ddy * ddy)
        Amp = Amp * rr / np.sqrt(rr * rr + 2.0 * xi2)
    return Amp * np.exp(1j * th)


def wrap(x):
    return (x + np.pi) % TWO_PI - np.pi


def core_positions(psi):
    th = np.angle(psi)
    dthx = wrap(np.roll(th, -1, 0) - th)
    dthy = wrap(np.roll(th, -1, 1) - th)
    W = (dthx + np.roll(dthy, -1, 0) - np.roll(dthx, -1, 1) - dthy) / TWO_PI
    Phi = np.abs(psi) ** 2
    out = {}
    for sign in (+1, -1):
        cells = np.argwhere(W * sign > 0.5)
        if len(cells) == 0:
            return None
        ic, jc = cells[0]
        off = np.arange(-6, 7)
        ii, jj = (ic + off) % N, (jc + off) % N
        Pw = Phi[np.ix_(ii, jj)]
        wgt = np.clip(0.5 * Phi_star - Pw, 0.0, None)
        OX, OY = np.meshgrid(off * dx, off * dx, indexing="ij")
        sx = float(np.sum(wgt * OX) / wgt.sum())
        sy = float(np.sum(wgt * OY) / wgt.sum())
        out[sign] = (coords[ic] + 0.5 * dx + sx, coords[jc] + 0.5 * dx + sy)
    ntot = float(np.rint(W).sum())
    return out, ntot


print("=== Phase 2: wir jako tlo stacjonarne (LOCK: Phase0_balance.md) ===")
print(f"N={N} L={L} dx={dx}  para: {V1} (+1), {V2} (-1)")
print()

# ---------------- G3: relaksacja i stacjonarnosc ----------------
psi = make_pair()
# diagnostyka ansatzu (errata: skoki fazy poza rdzeniami musza zniknac)
th0 = np.angle(psi)
gxa = np.abs(wrap(np.roll(th0, -1, 0) - th0))
gya = np.abs(wrap(np.roll(th0, -1, 1) - th0))
far = (np.hypot((X - V1[0] + L/2) % L - L/2, (Y - V1[1] + L/2) % L - L/2) > 4) & \
      (np.hypot((X - V2[0] + L/2) % L - L/2, (Y - V2[1] + L/2) % L - L/2) > 4)
print(f"diagnostyka ansatzu: max |dtheta| na linkach poza rdzeniami = "
      f"{max(gxa[far].max(), gya[far].max()):.4f} (wadliwy ansatz mial pi)")
pos0 = core_positions(psi)
print(f"  kret ansatzu: n_tot={pos0[1]:.0f}; pozycje: "
      f"+1 {pos0[0][+1]}, -1 {pos0[0][-1]}")
print()
H_prev, H_inc = energy_grad(psi), 0.0
pos1000 = None
for stp in range(1, 2001):
    psi = psi + dt_flow * rhs(psi)
    if stp % 10 == 0:
        H = energy_grad(psi)
        if H > H_prev:
            H_inc = max(H_inc, H - H_prev)
        H_prev = H
    if stp == 1000:
        pos1000 = core_positions(psi)
psi_bg = psi
pos2000 = core_positions(psi_bg)
# determinizm
psi_d = make_pair()
for stp in range(1, 2001):
    psi_d = psi_d + dt_flow * rhs(psi_d)
det = bool(np.array_equal(psi_bg, psi_d))
res = float(np.max(np.abs(rhs(psi_bg))))
drift = max(np.hypot(pos2000[0][s][0] - pos1000[0][s][0],
                     pos2000[0][s][1] - pos1000[0][s][1]) for s in (+1, -1))
print("--- G3: stacjonarnosc po 2000 krokach relaksacji ---")
print(f"  pozycje rdzeni (krok 2000): +1 {tuple(round(v,3) for v in pos2000[0][+1])}, "
      f"-1 {tuple(round(v,3) for v in pos2000[0][-1])}")
print(f"  |psi| daleko od rdzeni: min = "
      f"{float(np.abs(psi_bg)[far].min()):.4f} (s*={s_star:.4f}; wadliwy "
      f"ansatz dawal 0.10)")
print(f"  residuum max|RHS| = {res:.3e}  (prog 1e-4)")
print(f"  dryf pozycji rdzeni (krok 1000 -> 2000) = {drift:.5f}  (prog 0.02)")
print(f"  kret: n_tot = {pos2000[1]:.0f};  determinizm bitowy: {det};  "
      f"H_mono: {H_inc <= H_TOL}")
g3 = res <= 1e-4 and drift <= 0.02 and pos2000[1] == 0
print(f"  -> {'PASS' if g3 else 'FAIL'}")
print()

f_bg = np.abs(psi_bg)
th_bg = np.angle(psi_bg)
eith = np.exp(1j * th_bg)
Vpp_bg = Vpp(f_bg)
g_bg = gfun(f_bg)

# ---------------- G4a: widma ----------------
SHIFT = 20.0
rng = np.random.default_rng(SEED)

# L_toy (rzeczywisty skalar)
v = rng.standard_normal((N, N))
v /= np.linalg.norm(v)
for it in range(4000):
    Lv = -kappa * lap(v) + Vpp_bg * v
    u_new = SHIFT * v - Lv
    v = u_new / np.linalg.norm(u_new)
lam_toy = float(np.sum(v * (-kappa * lap(v) + Vpp_bg * v)))
RRc1 = np.hypot((X - V1[0] + L/2) % L - L/2, (Y - V1[1] + L/2) % L - L/2)
RRc2 = np.hypot((X - V2[0] + L/2) % L - L/2, (Y - V2[1] + L/2) % L - L/2)
frac_core_toy = float(np.sum((v * v)[(RRc1 < 5) | (RRc2 < 5)]))

# L_full (pelna linearyzacja zespolona)
def A_full(u):
    loc = np.conj(eith) * u
    al, be = np.real(loc), np.imag(loc)
    return -kappa * lap(u) + eith * (Vpp_bg * al + 1j * g_bg * be)


u = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N)))
u /= np.sqrt(np.sum(np.abs(u) ** 2))
for it in range(4000):
    Au = A_full(u)
    u_new = SHIFT * u - Au
    u = u_new / np.sqrt(np.sum(np.abs(u_new) ** 2))
lam_full = float(np.real(np.sum(np.conj(u) * A_full(u))))
print("--- G4a: widma operatorow na tle wirowym ---")
print(f"  lambda_min(L_toy)  = {lam_toy:+.5f}  "
      f"(frakcja modu w rdzeniach: {frac_core_toy:.3f})  "
      f"[oczekiwanie pre-code: < -0.01 — pierscien V''<0 wiaze]")
print(f"  lambda_min(L_full) = {lam_full:+.6f}  "
      f"(rozdzielczosc ~1e-3; prog G4a: >= -1e-3)")
g4a = lam_full >= -1e-3
print(f"  -> G4a: {'PASS' if g4a else 'FAIL'}")
print()

# ---------------- G4b: noise-run (empiryczny test tachionu) ----------------
print("--- G4b: noise-run, pelna dynamika II rzedu, tau=150 ---")
noise = 1e-6 * (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N)))
psi_w = psi_bg + noise
psi_prev = psi_w.copy()
mask_core = (RRc1 < 5) | (RRc2 < 5)
tw, ln_core, ln_glob = [], [], []
E0 = None
fl = {"nonfinite": False, "runaway": False}
Es = []
for stp in range(1, 3001):
    psi_next = 2.0 * psi_w - psi_prev + dt_w * dt_w * rhs(psi_w)
    if stp % 10 == 0:
        if not np.all(np.isfinite(psi_next)):
            fl["nonfinite"] = True
            break
        if np.max(np.abs(psi_next)) > S_GUARD:
            fl["runaway"] = True
            break
        dpsi = psi_w - psi_bg
        tw.append(stp * dt_w)
        ln_core.append(np.log(np.sqrt(np.mean(np.abs(dpsi[mask_core]) ** 2)) + 1e-300))
        ln_glob.append(np.log(np.sqrt(np.mean(np.abs(dpsi) ** 2)) + 1e-300))
        ut = (psi_next - psi_prev) / (2 * dt_w)
        E = energy_grad(psi_w) + float(np.sum(0.5 * np.abs(ut) ** 2) * dx * dx)
        Es.append(E)
    psi_prev, psi_w = psi_w, psi_next
tw = np.array(tw)
sel = tw >= 50.0


def slope(yv):
    A_ = np.vstack([np.ones(sel.sum()), tw[sel]]).T
    coef, *_ = np.linalg.lstsq(A_, np.array(yv)[sel], rcond=None)
    return float(coef[1])


gam_core, gam_glob = slope(ln_core), slope(ln_glob)
q = len(Es) // 4
E_drift = abs(np.mean(Es[-q:]) - np.mean(Es[:q])) / abs(Es[0])
print(f"  gamma_fit: rdzenie = {gam_core:+.5f}, globalnie = {gam_glob:+.5f} "
      f"(prog 0.01; B2 mialo 0.125)")
print(f"  flagi: nf={fl['nonfinite']} run={fl['runaway']};  "
      f"dryft sekularny energii = {E_drift:.2e} (prog 1e-4)")
g4b = (gam_core <= 0.01 and gam_glob <= 0.01
       and not (fl["nonfinite"] or fl["runaway"]))
print(f"  -> G4b: {'PASS' if g4b else 'FAIL'}")
print()

# ---------------- G5: demonstracja pakietu (obserwacja) ----------------
print("--- G5 (demo, obserwacja): pakiet amplitudowy omega=1.1, b=6 ---")
om = 1.1
k0 = np.sqrt((om * om - Vpp_star) / kappa)
u0, sx_, sy_ = 1e-3, 8.0, 10.0
x0, y0 = -60.0, V1[1] + 6.0
env = u0 * np.exp(-((X - x0) ** 2) / (2 * sx_ * sx_)
                  - ((Y - y0) ** 2) / (2 * sy_ * sy_))
alpha0 = env * np.cos(k0 * (X - x0))
alpha_prev = env * np.cos(k0 * (X - x0) + om * dt_w)
psi_w = psi_bg + eith * alpha0
psi_prev = psi_bg + eith * alpha_prev
mask_tr = X > -10.0
alphas = []
next_s = None
for stp in range(1, 5001):
    psi_next = 2.0 * psi_w - psi_prev + dt_w * dt_w * rhs(psi_w)
    tau = stp * dt_w
    if stp % 10 == 0:
        if np.max(np.abs(psi_next)) > S_GUARD or not np.all(np.isfinite(psi_next)):
            print(f"  PRZERWANE (guard) przy tau={tau:.1f}")
            break
        dpsi = psi_w - psi_bg
        w2 = np.abs(dpsi) ** 2
        Wt = w2[mask_tr].sum()
        if Wt > 0.3 * w2.sum():
            xc = float((w2 * X)[mask_tr].sum() / Wt)
            if next_s is None and xc > 0.0:
                next_s = tau
        if next_s is not None and tau >= next_s and len(alphas) < 5:
            ut = (psi_next - psi_prev) / (2 * dt_w)
            gx = (np.roll(dpsi, -1, 0) - np.roll(dpsi, 1, 0)) / (2 * dx)
            gy = (np.roll(dpsi, -1, 1) - np.roll(dpsi, 1, 1)) / (2 * dx)
            Px = float(-np.real(np.sum(np.conj(ut)[mask_tr] * gx[mask_tr])) * dx * dx)
            Py = float(-np.real(np.sum(np.conj(ut)[mask_tr] * gy[mask_tr])) * dx * dx)
            alphas.append((Px, Py))
            next_s += 4.0
        if len(alphas) >= 5:
            break
    psi_prev, psi_w = psi_w, psi_next
if alphas:
    Pxm = float(np.mean([p[0] for p in alphas]))
    Pym = float(np.mean([p[1] for p in alphas]))
    al = float(np.degrees(np.arctan2(Pym, Pxm)))   # srednia WEKTOROWA
    al_s = [float(np.degrees(np.arctan2(p[1], p[0]))) for p in alphas]
    # wir przy y=-32, pakiet startuje y=-26: KU wirowi = P_y < 0 = alpha < 0
    kier = ("KU wirowi" if al < -0.5 else ("OD wiru" if al > 0.5
                                           else "nierozstrzygniete"))
    print(f"  alpha_demo (srednia wektorowa P) = {al:+.3f} deg; "
          f"probki katowe: {[f'{x:+.1f}' for x in al_s]}  ->  kierunek: {kier}")
    print(f"  (eikonal, Phase0: b=6, omega=1.1 -> +43 deg KU wirowi; "
          f"pomiar ILOSCIOWY = nastepny op)")
else:
    print("  brak probek (pakiet nie doszedl do okna) — raport")
print()

# ---------------- EWALUACJA ----------------
print("=== EWALUACJA (LOCKED) ===")
print(f"G1 (czesc Phase2): determinizm {det}; H_mono relaksacji "
      f"{H_inc <= H_TOL}; dryft energii dynamiki {E_drift:.2e} <= 1e-4: "
      f"{E_drift <= 1e-4}")
print(f"G3 (stacjonarnosc): {'PASS' if g3 else 'FAIL'}")
print(f"G4 (RDZEN, brak tachionu): a) lambda_min(L_full)={lam_full:+.6f} "
      f">= -1e-3: {g4a};  b) gamma_fit<=0.01: {g4b}  "
      f"-> {'PASS' if (g4a and g4b) else 'FAIL'}")
print(f"   diagnostyka: lambda_min(L_toy)={lam_toy:+.5f} "
      f"({'toy-sektor ma tachion na tle stacjonarnym — zdyskwalifikowany' if lam_toy < -0.01 else 'toy-sektor bez szybkiego tachionu'})")
print()
print(f"PODSUMOWANIE PHASE 2: G3={'PASS' if g3 else 'FAIL'} "
      f"G4={'PASS' if (g4a and g4b) else 'FAIL'}")
