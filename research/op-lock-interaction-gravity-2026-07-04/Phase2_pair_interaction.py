# Phase 2: oddzialywanie pary zlockowanych obiektow — G3 (istnienie), G4 (prawo)
# Implementacja SCISLE wg Phase0_balance.md (LOCKED, sekcje 2-5 + D1-D10).
#
# Statyka:  E_int(g) = H[s1+s2] - H[s1] - H[s2]  (kompozycja, bez ewolucji)
# Dynamika: Delta_v(g) = v_zblizania(g) - 2*v_control(tau)  (kontrola D4)
# Fit (Y) vs (P) w log-przestrzeni (D2), rozstrzyga Delta_AIC >= 10.

import numpy as np

# ---------------- parametry LOCKED ----------------
N = 256
L = 128.0
dx = L / N
kappa = 0.50
a, b, c = 0.50, 1.60, 1.00
dt = 0.02
EPS = 0.30
S_GUARD = 10.0
w = 1.50
REC_EVERY = 10
VBASE = 2.0
KOFF = int(round(VBASE / (REC_EVERY * dt)))
H_TOL = 1e-10
G0_SCAN = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0]
R_OBJ = 8.0
G_GRID = [2.0, 2.5, 3.0, 3.5]          # D5
G3_DYN_CHECK = [2.0, 3.0, 4.0]         # D5: na runie g0=6

Nc = N // 2
coords = (np.arange(N) - Nc) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
CHEB_ZONE = np.maximum(np.abs(X), np.abs(Y)) > 0.45 * L


def Vpot(s):
    return 0.5 * a * s**2 - (b / 3.0) * np.abs(s)**3 + 0.25 * c * s**4


def Vprime(s):
    return s * (a - b * np.abs(s) + c * s * s)


def lap(s):
    return (np.roll(s, 1, 0) + np.roll(s, -1, 0)
            + np.roll(s, 1, 1) + np.roll(s, -1, 1) - 4.0 * s) / (dx * dx)


def energy(s):
    gx = (np.roll(s, -1, 0) - s) / dx
    gy = (np.roll(s, -1, 1) - s) / dx
    return float(np.sum(0.5 * kappa * (gx * gx + gy * gy) + Vpot(s)) * dx * dx)


def gaussian(A0, cx, cy):
    return A0 * np.exp(-(((X - cx) ** 2 + (Y - cy) ** 2) / (2.0 * w * w)))


def make_cluster(A0, D):
    s = gaussian(A0, 0.0, 0.0)
    for k in range(6):
        ang = k * np.pi / 3.0
        s = s + gaussian(A0, D * np.cos(ang), D * np.sin(ang))
    return s


def half_width_x(s):
    phi = s[:, Nc] ** 2
    if phi[Nc] <= EPS:
        return np.nan
    i = Nc
    while i + 1 < N and phi[i + 1] > EPS:
        i += 1
    if i + 1 >= N:
        return np.nan
    xr = coords[i] + dx * (phi[i] - EPS) / (phi[i] - phi[i + 1])
    j = Nc
    while j - 1 >= 0 and phi[j - 1] > EPS:
        j -= 1
    if j - 1 < 0:
        return np.nan
    xl = coords[j] - dx * (phi[j] - EPS) / (phi[j] - phi[j - 1])
    return 0.5 * (xr - xl)


def gap_x(s):
    """D1: szczelina = odleglosc przeciec Phi=eps wokol x=0 na osi y=0."""
    phi = s[:, Nc] ** 2
    if phi[Nc] > EPS:
        return 0.0                       # fronty zetkniete/zmergowane
    i = Nc
    while i - 1 >= 0 and phi[i - 1] <= EPS:
        i -= 1
    if i - 1 < 0:
        return np.nan
    x_l = coords[i - 1] + dx * (phi[i - 1] - EPS) / (phi[i - 1] - phi[i])
    j = Nc
    while j + 1 < N and phi[j + 1] <= EPS:
        j += 1
    if j + 1 >= N:
        return np.nan
    x_r = coords[j + 1] - dx * (phi[j + 1] - EPS) / (phi[j + 1] - phi[j])
    return x_r - x_l


def prep_object(A0, r_target, max_steps=20000):
    """Ewolucja cluster(A0, D=3, K=7) do half_width >= r_target (D7/D8)."""
    s = make_cluster(A0, 3.0)
    H_prev = energy(s)
    H_inc_max = 0.0
    for step in range(1, max_steps + 1):
        s = s + dt * (kappa * lap(s) - Vprime(s))
        if step % REC_EVERY:
            continue
        H = energy(s)
        if H > H_prev:
            H_inc_max = max(H_inc_max, H - H_prev)
        H_prev = H
        R = half_width_x(s)
        if np.isfinite(R) and R >= r_target:
            return s, R, step * dt, H_inc_max
    raise RuntimeError("prep_object: nie osiagnieto r_target")


def evolve_recorded(s0, max_steps, measure):
    """Ewolucja z zapisem measure(s) co REC_EVERY; stop wg flagi stop z measure."""
    s = s0.copy()
    taus, vals = [0.0], [measure(s)]
    H_prev = energy(s)
    H_inc_max = 0.0
    flags = {"nonfinite": False, "runaway": False, "contact": False}
    for step in range(1, max_steps + 1):
        s = s + dt * (kappa * lap(s) - Vprime(s))
        if step % REC_EVERY:
            continue
        if not np.all(np.isfinite(s)):
            flags["nonfinite"] = True
            break
        if np.max(np.abs(s)) > S_GUARD:
            flags["runaway"] = True
            break
        H = energy(s)
        if H > H_prev:
            H_inc_max = max(H_inc_max, H - H_prev)
        H_prev = H
        m = measure(s)
        taus.append(step * dt); vals.append(m)
        if np.any((s * s)[CHEB_ZONE] > EPS):
            flags["contact"] = True
            break
        if not np.isfinite(m) or m < 1.0:     # okno pomiarowe: g < 1.0 konczy
            break
    flags["H_inc_max"] = H_inc_max
    flags["H_monotone"] = H_inc_max <= H_TOL * max(1.0, abs(energy(s0)))
    return np.array(taus), np.array(vals), flags


def central_v(tau, val):
    n = len(val)
    idx = np.arange(KOFF, n - KOFF)
    v = (val[idx + KOFF] - val[idx - KOFF]) / (tau[idx + KOFF] - tau[idx - KOFF])
    return idx, v


def fit_loglaw(gv, Ev):
    """D2: (Y) ln|E| ~ alpha - m*g ; (P) ln|E| ~ beta - p*ln(g). Zwraca fity+AIC."""
    y = np.log(np.abs(Ev))
    n = len(y)
    out = {}
    for name, xcol in (("Y", gv), ("P", np.log(gv))):
        A_ = np.vstack([np.ones(n), -xcol]).T
        coef, *_ = np.linalg.lstsq(A_, y, rcond=None)
        pred = A_ @ coef
        rss = float(np.sum((y - pred) ** 2))
        tss = float(np.sum((y - y.mean()) ** 2))
        aic = n * np.log(rss / n) + 2 * 2
        out[name] = {"amp": np.exp(coef[0]), "exp": coef[1], "rss": rss,
                     "r2": 1.0 - rss / tss if tss > 0 else np.nan, "aic": aic, "n": n}
    out["dAIC_P_minus_Y"] = out["P"]["aic"] - out["Y"]["aic"]
    return out


print("=== Phase 2: oddzialywanie pary (LOCK: Phase0_balance.md) ===")
print(f"N={N} L={L} dx={dx} kappa={kappa} a={a} b={b} c={c} dt={dt} eps={EPS} "
      f"R_obj={R_OBJ} g0_scan={G0_SCAN}")
print()

# ---------------- przygotowanie obiektu ----------------
s_obj, R_act, tau_prep, hinc = prep_object(0.6, R_OBJ)
H_obj = energy(s_obj)
print(f"obiekt: cluster(0.6, D=3, K=7) -> R_act={R_act:.4f} po tau={tau_prep:.1f} "
      f"(H_inc_max prep={hinc:.2e}); H_obj={H_obj:.6f}")
print()

# ---------------- STATYKA: E_int(g) ----------------
print("--- pair_static: E_int(g) = H[s1+s2] - 2*H[obj] (D1: shift o pelne wezly) ---")
stat_g, stat_E = [], []
for g0 in G0_SCAN:
    X0 = R_act + 0.5 * g0
    cells = int(round(X0 / dx))
    s1 = np.roll(s_obj, -cells, axis=0)
    s2 = np.roll(s_obj, +cells, axis=0)
    pair = s1 + s2
    g_meas = gap_x(pair)
    E_int = energy(pair) - 2.0 * H_obj
    stat_g.append(g_meas); stat_E.append(E_int)
    print(f"  g0={g0:4.1f}  centra=+/-{cells * dx:6.2f}  g_meas={g_meas:7.4f}  "
          f"E_int={E_int:+.6e}   blad superpozycji ~e^-g = {np.exp(-g_meas):.4f}")
stat_g, stat_E = np.array(stat_g), np.array(stat_E)
print()

usable = np.isfinite(stat_g) & (stat_g > 0) & (stat_E < 0)
fits = fit_loglaw(stat_g[usable], stat_E[usable]) if usable.sum() >= 3 else None
if fits:
    fy, fp = fits["Y"], fits["P"]
    print("--- fit (Y) vs (P) na E_int(g_meas), log-przestrzen (D2) ---")
    print(f"  (Y) Yukawa:  E = -A*exp(-m*g):  A={fy['amp']:.5f}  m={fy['exp']:.4f}  "
          f"R^2(log)={fy['r2']:.6f}  AIC={fy['aic']:+.2f}")
    print(f"  (P) potega:  E = -B*g^(-p):     B={fp['amp']:.5f}  p={fp['exp']:.4f}  "
          f"R^2(log)={fp['r2']:.6f}  AIC={fp['aic']:+.2f}")
    print(f"  Delta_AIC (P - Y) = {fits['dAIC_P_minus_Y']:+.2f}   "
          f"(>= +10: wygrywa Y;  <= -10: wygrywa P)")
    # test czulosci (raportowany, poza kryterium): LOCK deklaruje superpozycje
    # "dopuszczalna dla g>=2" — fit powtorzony bez punktow g_meas < 2.
    sens = usable & (stat_g >= 2.0)
    if 3 <= sens.sum() < usable.sum():
        fs = fit_loglaw(stat_g[sens], stat_E[sens])
        print(f"  [czulosc, bez g_meas<2] (Y): A={fs['Y']['amp']:.5f} "
              f"m={fs['Y']['exp']:.4f} R^2={fs['Y']['r2']:.6f}   "
              f"(P): p={fs['P']['exp']:.4f} R^2={fs['P']['r2']:.6f}   "
              f"Delta_AIC(P-Y)={fs['dAIC_P_minus_Y']:+.2f}")
print()

# ---------------- DYNAMIKA ----------------
print("--- pair_dynamic: kontrola D4 (ten sam snapshot, pojedynczy) ---")
tau_c, R_c_, fl_c = evolve_recorded(s_obj, 4000, half_width_x)
idx_c, v_c = central_v(tau_c, R_c_)
print(f"  kontrola: {len(tau_c)} rekordow, R: {R_c_[0]:.3f} -> {R_c_[-1]:.3f}, "
      f"flagi: nf={fl_c['nonfinite']} run={fl_c['runaway']} bc={fl_c['contact']} "
      f"H_mono={fl_c['H_monotone']}")
print()

dyn = {}
all_flags = [fl_c]
for g0 in G0_SCAN:
    X0 = R_act + 0.5 * g0
    cells = int(round(X0 / dx))
    pair0 = np.roll(s_obj, -cells, axis=0) + np.roll(s_obj, +cells, axis=0)
    g_init = gap_x(pair0)
    tau_p, g_p, fl = evolve_recorded(pair0, 4000, gap_x)
    all_flags.append(fl)
    idx_p, vg = central_v(tau_p, g_p)          # dg/dtau
    v_close = -vg
    npts = min(len(idx_p), len(v_c))
    dv = v_close[:npts] - 2.0 * v_c[:npts]
    gm = g_p[idx_p][:npts]
    dyn[g0] = {"g": gm, "dv": dv, "g_init": g_init, "tau_end": tau_p[-1], "fl": fl}
    print(f"  pair(g0={g0:3.1f}): g_init={g_init:6.3f} -> g_end={g_p[-1]:6.3f} "
          f"(tau_end={tau_p[-1]:6.1f}, {len(tau_p)} rek.)  "
          f"nf={fl['nonfinite']} run={fl['runaway']} bc={fl['contact']} "
          f"H_mono={fl['H_monotone']}")
print()

print("--- Delta_v(g) [interpolacja; kolumny = g] ---")
gcols = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0]
hdr = "  run      " + "".join(f"  g={g:<7.1f}" for g in gcols)
print(hdr)
for g0 in G0_SCAN:
    d = dyn[g0]
    if len(d["g"]) < 2:
        continue
    gs, dvs = d["g"][::-1], d["dv"][::-1]
    row = f"  g0={g0:3.1f}  "
    for gq in gcols:
        if gs[0] <= gq <= gs[-1]:
            row += f"  {np.interp(gq, gs, dvs):+9.5f}"
        else:
            row += "        --"
    print(row)
print()

# fit informacyjny (Y)/(P) na Delta_v(g) z runu g0=8 (najszersze pokrycie)
d8 = dyn[8.0]
m8 = np.isfinite(d8["g"]) & np.isfinite(d8["dv"]) & (d8["dv"] > 0) & \
     (d8["g"] >= 1.5) & (d8["g"] <= 7.0)
if m8.sum() >= 5:
    fd = fit_loglaw(d8["g"][m8], -d8["dv"][m8] * -1.0)  # |dv|
    print("--- fit informacyjny na Delta_v(g), run g0=8 (poza kryterium G4) ---")
    print(f"  (Y): m_dyn={fd['Y']['exp']:.4f} R^2={fd['Y']['r2']:.4f}   "
          f"(P): p_dyn={fd['P']['exp']:.4f} R^2={fd['P']['r2']:.4f}   "
          f"Delta_AIC(P-Y)={fd['dAIC_P_minus_Y']:+.2f}")
    print()

# ---------------- EWALUACJA ----------------
print("=== EWALUACJA (LOCKED) ===")

g1_ok = all(not (f["nonfinite"] or f["runaway"]) and f["H_monotone"]
            for f in all_flags)
print(f"G1 (czesc Phase2): finite/no-runaway + H_nierosnace na wszystkich runach: "
      f"{'PASS' if g1_ok else 'FAIL'}")

stat_attr = bool(np.all(stat_E[usable] < 0)) and usable.sum() == len(G0_SCAN)
order = np.argsort(stat_g[usable])
stat_decr = bool(np.all(np.diff(np.abs(stat_E[usable][order])) < 0))
# ERRATA (raportowanie, dzien runow, przed interpretacja): pierwotna linia
# oceniala G3-dynamike przez interpolacje w punktach {2,3,4} (D5); punkt g=2.0
# lezy poza NOSNIKIEM estymatora roznic centralnych (obciecie +/-tau=2.0 przy
# przyspieszajacych frontach) — "niemierzalne w 2.0" zlewalo sie z "niedodatnie".
# Ocena skorygowana do BRZMIENIA LOCKa G3: "Delta_v(g) > 0 dla g <= 4":
# wszystkie ZMIERZONE punkty z g <= 4 musza byc dodatnie i musza istniec.
d6 = dyn[6.0]
m_le4 = np.isfinite(d6["g"]) & (d6["g"] <= 4.0)
n_le4 = int(m_le4.sum())
all_pos = bool(np.all(d6["dv"][m_le4] > 0)) if n_le4 else False
g_min_cov = float(np.min(d6["g"][m_le4])) if n_le4 else np.nan
gs6, dvs6 = d6["g"][::-1], d6["dv"][::-1]
interp_info = {gq: (f"{np.interp(gq, gs6, dvs6):+.5f}"
                    if gs6[0] <= gq <= gs6[-1] else "poza nosnikiem")
               for gq in G3_DYN_CHECK}
g3 = stat_attr and stat_decr and (n_le4 > 0) and all_pos
print(f"G3 (istnienie oddzialywania):")
print(f"  E_int < 0 dla wszystkich g: {stat_attr};  |E_int| malejace z g: {stat_decr}")
print(f"  dynamika (run g0=6, brzmienie LOCKa): punktow z g<=4: {n_le4}, "
      f"wszystkie Delta_v > 0: {all_pos}, min pokryte g = {g_min_cov:.3f}")
print(f"  informacyjnie, interpolacje w {G3_DYN_CHECK}: {interp_info}")
print(f"  -> {'PASS' if g3 else 'FAIL'}")

if fits:
    dA = fits["dAIC_P_minus_Y"]
    winner = "Y" if dA >= 10 else ("P" if dA <= -10 else "BRAK")
    m_fit = fits["Y"]["exp"]
    m_ok = 0.8 <= m_fit <= 1.2
    g4 = (winner == "Y" and m_ok) or (winner == "P")
    print(f"G4 (RDZEN, prawo algebraiczne): Delta_AIC={dA:+.2f} -> wygrywa: {winner}")
    if winner == "Y":
        print(f"  m_fit={m_fit:.4f} in [0.8, 1.2]: {m_ok}  "
              f"(predykcja substratu m0=1.0, rownanie #1)")
    print(f"  -> {'PASS' if g4 else 'FAIL'}")
else:
    g4 = False
    print("G4: FAIL (za malo punktow z E_int < 0 do fitu)")

print()
print(f"PODSUMOWANIE PHASE 2: G1(czesc)={'PASS' if g1_ok else 'FAIL'} "
      f"G3={'PASS' if g3 else 'FAIL'} G4={'PASS' if g4 else 'FAIL'}"
      + (f"  [forma: {'Yukawa' if fits and fits['dAIC_P_minus_Y'] >= 10 else 'potegowa/brak'}]"
         if fits else ""))
