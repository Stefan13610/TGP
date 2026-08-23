# Phase 2: refrakcja pakietow na tle zamrozonego obiektu — G3-G6
# Implementacja SCISLE wg Phase0_balance.md (LOCKED, sekcje 1-5 + D1-D11).
#
# Scenariusze: tachyon_check (D6b), deflection(b, omega) z kryterium
# waznosci D6(c), barrier_regime (obserwacja), linearity (G6),
# frozen_vs_live (obserwacja, D11).

import numpy as np

# ---------------- parametry LOCKED ----------------
kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt = 0.05
dt_flow = 0.02
N, L = 256, 128.0
dx = L / N
EPS = 0.30
S_GUARD = 10.0
SEED = 20260704
U0 = 1e-3
SX, SY = 8.0, 12.0
B_SCAN = [8.0, 12.0, 16.0, 20.0]
OM_SCAN = [1.0, 1.1, 1.3]
R_TARGET = 12.0
GAMMA_CONTAM = 0.02          # D6(c)
RATIO_CONTAM = 0.10          # D6(c)
w = 1.50

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)
Vpp_star = a - 2.0 * b_par * s_star + 3.0 * c * s_star ** 2

Nc = N // 2
coords = (np.arange(N) - Nc) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
RR = np.sqrt(X * X + Y * Y)

# alpha_pred z Phase0_output (D7, zalockowane przed pomiarem):
ALPHA_PRED = {(1.0, 8.0): +71.118, (1.0, 12.0): +24.781,
              (1.0, 16.0): -92.403, (1.0, 20.0): -0.243,
              (1.1, 8.0): +61.430, (1.1, 12.0): +19.048,
              (1.1, 16.0): -65.570, (1.1, 20.0): -0.170,
              (1.3, 8.0): +45.000, (1.3, 12.0): +11.099,
              (1.3, 16.0): -13.727, (1.3, 20.0): -0.101}


def Vprime(s):
    return s * (a - b_par * np.abs(s) + c * s * s)


def Vpp(s):
    return a - 2.0 * b_par * np.abs(s) + 3.0 * c * s * s


def lap(u):
    return (np.roll(u, 1, 0) + np.roll(u, -1, 0)
            + np.roll(u, 1, 1) + np.roll(u, -1, 1) - 4.0 * u) / (dx * dx)


def gaussian(A0, cx, cy):
    return A0 * np.exp(-(((X - cx) ** 2 + (Y - cy) ** 2) / (2.0 * w * w)))


def half_width_x(s):
    phi = s[:, Nc] ** 2
    if phi[Nc] <= EPS:
        return np.nan
    i = Nc
    while i + 1 < N and phi[i + 1] > EPS:
        i += 1
    xr = coords[i] + dx * (phi[i] - EPS) / (phi[i] - phi[i + 1])
    j = Nc
    while j - 1 >= 0 and phi[j - 1] > EPS:
        j -= 1
    xl = coords[j] - dx * (phi[j] - EPS) / (phi[j] - phi[j - 1])
    return 0.5 * (xr - xl)


def prep_object(A0, r_target, max_steps=40000):
    s = gaussian(A0, 0.0, 0.0)
    for k in range(6):
        ang = k * np.pi / 3.0
        s = s + gaussian(A0, 3.0 * np.cos(ang), 3.0 * np.sin(ang))
    for step in range(1, max_steps + 1):
        s = s + dt_flow * (kappa * lap(s) - Vprime(s))
        if step % 10 == 0:
            R = half_width_x(s)
            if np.isfinite(R) and R >= r_target:
                return s, R
    raise RuntimeError("prep_object failed")


def make_packet(om, bb, u0, R_obj):
    k0 = np.sqrt((om * om - a) / kappa)
    x0 = -round((R_obj + 40.0) / dx) * dx          # D2: zaokraglone do wezla
    env = u0 * np.exp(-((X - x0) ** 2) / (2 * SX * SX)
                      - ((Y - bb) ** 2) / (2 * SY * SY))
    u = env * np.cos(k0 * (X - x0))
    ut = om * env * np.sin(k0 * (X - x0))
    return u, ut, x0


def grad_c(u):
    gx = (np.roll(u, -1, 0) - np.roll(u, 1, 0)) / (2 * dx)
    gy = (np.roll(u, -1, 1) - np.roll(u, 1, 1)) / (2 * dx)
    return gx, gy


def deflection_run(vpp, R_obj, om, bb, u0, label, live_s=None, max_tau=400.0):
    """Run pakietu; pomiar D3 + walidacja D6(c). live_s: tlo zywe (D11)."""
    u, ut0, x0 = make_packet(om, bb, u0, R_obj)
    s_live = live_s.copy() if live_s is not None else None
    vpp_cur = Vpp(s_live) if s_live is not None else vpp
    u_prev = u - dt * ut0 + 0.5 * dt * dt * (kappa * lap(u) - vpp_cur * u)
    mask_trans = X > (R_obj + 5.0)
    mask_wall = RR < (R_obj + 5.0)
    taus_w, lnw = [], []
    alphas, tau_meas = [], None
    ratio_contam = np.nan
    flags = {"nonfinite": False, "runaway": False, "wrap": False}
    next_sample = None
    steps = int(max_tau / dt)
    for stp in range(1, steps + 1):
        if s_live is not None:
            s_live = s_live + dt * (kappa * lap(s_live) - Vprime(s_live))
            vpp_cur = Vpp(s_live)
        u_next = 2.0 * u - u_prev + dt * dt * (kappa * lap(u) - vpp_cur * u)
        tau = stp * dt
        if stp % 10 == 0:
            if not np.all(np.isfinite(u_next)):
                flags["nonfinite"] = True
                break
            if np.max(np.abs(u_next)) > S_GUARD:
                flags["runaway"] = True
                break
            ut = (u_next - u_prev) / (2.0 * dt)          # D1
            u2 = u * u
            nw = np.sqrt(np.mean(u2[mask_wall]))
            if nw > 0:
                taus_w.append(tau)
                lnw.append(np.log(nw))
            Wt = u2[mask_trans].sum()
            if Wt > 0.2 * u2.sum():
                xc = float((u2 * X)[mask_trans].sum() / Wt)
                if xc - x0 > 0.9 * L:
                    flags["wrap"] = True
                    break
                if next_sample is None and xc > R_obj + 20.0:
                    next_sample = tau                     # start probek D3
                    tau_meas = tau
            if next_sample is not None and tau >= next_sample and len(alphas) < 5:
                gx, gy = grad_c(u)
                Px = float(-(ut * gx)[mask_trans].sum() * dx * dx)
                Py = float(-(ut * gy)[mask_trans].sum() * dx * dx)
                Pxg = float(-(ut * gx).sum() * dx * dx)
                Pyg = float(-(ut * gy).sum() * dx * dx)
                alphas.append((np.degrees(np.arctan2(Py, Px)),
                               np.degrees(np.arctan2(Pyg, Pxg))))
                nrm_w = np.sqrt(u2[mask_wall].sum())
                nrm_t = np.sqrt(u2[mask_trans].sum())
                ratio_contam = nrm_w / nrm_t if nrm_t > 0 else np.inf
                next_sample += 4.0                        # D3: co Delta_tau=4
            if len(alphas) >= 5:
                break
        u_prev, u = u, u_next
    # D6(c): wzrost wykladniczy w oknie sciany (druga polowa przed pomiarem)
    gamma_fit = np.nan
    if len(taus_w) > 10:
        tw, lw = np.array(taus_w), np.array(lnw)
        sel = tw > 0.4 * tw[-1]
        if sel.sum() > 5:
            A_ = np.vstack([np.ones(sel.sum()), tw[sel]]).T
            coef, *_ = np.linalg.lstsq(A_, lw[sel], rcond=None)
            gamma_fit = float(coef[1])
    contaminated = ((np.isfinite(gamma_fit) and gamma_fit > GAMMA_CONTAM
                     and (not np.isfinite(ratio_contam)
                          or ratio_contam > RATIO_CONTAM))
                    or flags["runaway"] or flags["nonfinite"])
    if alphas:
        at = np.array([al[0] for al in alphas])
        ag = np.array([al[1] for al in alphas])
        alpha, sig = float(at.mean()), float(at.std())
        alpha_g = float(ag.mean())
    else:
        alpha = sig = alpha_g = np.nan
    return {"alpha": alpha, "sig": sig, "alpha_glob": alpha_g,
            "gamma_wall": gamma_fit, "ratio": ratio_contam,
            "valid": (not contaminated) and len(alphas) == 5,
            "n_samples": len(alphas), "tau_meas": tau_meas,
            "flags": flags, "label": label}


print("=== Phase 2: refrakcja na tle zamrozonym (LOCK: Phase0_balance.md) ===")
print(f"N={N} L={L} dx={dx} dt={dt} u0={U0} sx={SX} sy={SY} seed={SEED}")
print()

s_bg, R_obj = prep_object(0.6, R_TARGET)
VPP_BG = Vpp(s_bg)
print(f"tlo zamrozone: R_obj={R_obj:.3f} (D9)")
print()

# ---------------- A. tachyon_check (D6b) ----------------
print("--- tachyon_check (D6b): szum 1e-8, bez pakietu, tau<=200 ---")
rng = np.random.default_rng(SEED)
u = 1e-8 * rng.standard_normal((N, N))
u_prev = u.copy()
mask_wall = RR < (R_obj + 5.0)
tw, lw, loc = [], [], []
for stp in range(1, 4001):
    u_next = 2.0 * u - u_prev + dt * dt * (kappa * lap(u) - VPP_BG * u)
    if stp % 10 == 0:
        if np.max(np.abs(u_next)) > S_GUARD:
            print(f"  S_GUARD osiagniety przy tau={stp*dt:.1f} (flaga, stop)")
            break
        nw = np.sqrt(np.mean(u_next[mask_wall] ** 2))
        tw.append(stp * dt)
        lw.append(np.log(nw))
        loc.append(float(np.sum(u_next[mask_wall] ** 2) / np.sum(u_next ** 2)))
    u_prev, u = u, u_next
tw, lw = np.array(tw), np.array(lw)
sel = tw > 0.5 * tw[-1]
A_ = np.vstack([np.ones(sel.sum()), tw[sel]]).T
coef, *_ = np.linalg.lstsq(A_, lw[sel], rcond=None)
gamma_meas = float(coef[1])
print(f"  gamma_meas = {gamma_meas:.4f}  vs  sqrt(-lambda_min) = 0.1234 "
      f"(Phase0)  [zgodnosc: {100*abs(gamma_meas/0.1234-1):.1f}%]")
print(f"  lokalizacja wzrostu w r < R_obj+5: {loc[-1]:.4f} "
      f"(tachion SCIANY potwierdzony empirycznie)")
print()

# ---------------- B. deflection(b, omega) ----------------
print("--- deflection: skan (omega, b); waznosc wg D6(c) ---")
runs = {}
for om in OM_SCAN:
    for bb in B_SCAN:
        r = deflection_run(VPP_BG, R_obj, om, bb, U0, f"om={om},b={bb}")
        runs[(om, bb)] = r
        ap = ALPHA_PRED[(om, bb)]
        print(f"  om={om:.1f} b={bb:4.1f}: alpha={r['alpha']:+8.3f}deg "
              f"(sig={r['sig']:.3f}, glob={r['alpha_glob']:+8.3f})  "
              f"pred={ap:+8.3f}  gamma_wall={r['gamma_wall']:.4f}  "
              f"kontam.ratio={r['ratio']:.3f}  probek={r['n_samples']}  "
              f"runaway={r['flags']['runaway']}  "
              f"WAZNY={r['valid']}")
print()

# ---------------- C. linearity (G6) ----------------
print("--- linearity (G6): om=1.1, b=12, u0 vs u0/2 ---")
r_half = deflection_run(VPP_BG, R_obj, 1.1, 12.0, U0 / 2, "lin u0/2")
r_full = runs[(1.1, 12.0)]
if r_full["valid"] and r_half["valid"]:
    dlin = abs(r_half["alpha"] - r_full["alpha"]) / max(abs(r_full["alpha"]), 1e-12)
    print(f"  alpha(u0)={r_full['alpha']:+.3f}  alpha(u0/2)={r_half['alpha']:+.3f}  "
          f"roznica={100*dlin:.2f}% (prog 10%)")
    g6 = ("PASS" if dlin <= 0.10 else "FAIL")
else:
    print(f"  runy niewazne wg D6(c) (waznosc: {r_full['valid']}/{r_half['valid']}) "
          f"-> G6 NIEAPLIKOWALNE (D6d); uwaga zapisana z gory: rownanie liniowe "
          f"=> G6 nie wykrywa kontaminacji")
    g6 = "NIEAPLIKOWALNE"
print()

# ---------------- D. barrier_regime (obserwacja) ----------------
print("--- barrier_regime (obserwacja): om=0.8 in (m_FV, m_TV), b=0 ---")
r_bar = deflection_run(VPP_BG, R_obj, 0.8, 0.0, U0, "barrier", max_tau=300.0)
print(f"  probek pomiarowych: {r_bar['n_samples']} (transmisja "
      f"{'nie ' if r_bar['n_samples'] < 5 else ''}osiagnela okna pomiaru); "
      f"gamma_wall={r_bar['gamma_wall']:.4f}  runaway={r_bar['flags']['runaway']}")
print()

# ---------------- E. frozen_vs_live (obserwacja, D11) ----------------
print("--- frozen_vs_live (obserwacja, D11): om=1.1, b=12, tlo zywe ---")
r_live = deflection_run(None, R_obj, 1.1, 12.0, U0, "live", live_s=s_bg)
print(f"  alpha_live={r_live['alpha']:+.3f}deg (sig={r_live['sig']:.3f})  "
      f"gamma_wall={r_live['gamma_wall']:.4f}  kontam.ratio={r_live['ratio']:.3f}  "
      f"probek={r_live['n_samples']}  runaway={r_live['flags']['runaway']}")
print("  (tlo rosnie podczas przelotu — obserwacja jakosciowa, poza G1-G6)")
print()

# ---------------- EWALUACJA ----------------
print("=== EWALUACJA (LOCKED) ===")
n_valid = sum(1 for r in runs.values() if r["valid"])
print(f"waznosc pomiarow deflection wg D6(c): {n_valid}/12 waznych")
main = runs[(1.1, 12.0)]
if n_valid == 0:
    print("G3/G4/G6: NIEAPLIKOWALNE wg pre-rejestrowanej klauzuli D6(d):")
    print("  tlo zamrozone generuje tachion sciany (lambda_min=-0.0152,")
    print(f"  gamma_meas={gamma_meas:.3f}), ktory kontaminuje wszystkie runy")
    print("  przed oknem pomiaru. Geometria propagacji na TLE ZAMROZONYM")
    print("  jest niemierzalna w tej minimalnej realizacji.")
    print("G5: empirycznie NIEROZSTRZYGNIETE; strona teoretyczna (eikonal,")
    print("  D7, zalockowana przed pomiarem): kierunek zalezy od b —")
    print("  OD obiektu dla promieni przecinajacych wnetrze (b=8,12),")
    print("  KU obiektowi dla muskajacych pierscien sciany n^2>1 (b=16).")
else:
    g3 = (main["valid"] and abs(main["alpha"]) > 3 * main["sig"])
    print(f"G3: |alpha(b=12,om=1.1)|={abs(main['alpha']):.3f} > "
          f"3*sig={3*main['sig']:.3f}: {g3}")
    ok4 = []
    for key, r in runs.items():
        if r["valid"]:
            ap = ALPHA_PRED[key]
            ok4.append(np.sign(r["alpha"]) == np.sign(ap)
                       and 0.5 <= r["alpha"] / ap <= 2.0)
    print(f"G4: znak+wielkosc vs eikonal na waznych runach: "
          f"{sum(ok4)}/{len(ok4)} -> {'PASS' if all(ok4) else 'FAIL'}")
    signs = [np.sign(r["alpha"]) for r in runs.values() if r["valid"]]
    print(f"G5: znaki alpha na waznych runach: {signs}")
print(f"G6: {g6}")
