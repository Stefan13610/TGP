# Phase 1: kalibracja sektora falowego — dyspersja omega(k) (G2 + czesc G1)
# Implementacja SCISLE wg Phase0_balance.md (LOCKED, sekcja 2-4 + D1, D4, D5).
#
# Sektor falowy: d^2u/dtau^2 = kappa*Lap(u) - V''(s_bg)*u, leapfrog (D1).
# Tla jednorodne: FV (s_bg = 0) i TV (s_bg = s*).
# Pomiar: fale plaskie k = 2*pi*m/L, omega ze zliczania przejsc przez zero.

import numpy as np

# ---------------- parametry LOCKED ----------------
kappa = 0.50
a, b_par, c = 0.50, 1.60, 1.00
dt = 0.05
N, L = 256, 128.0
dx = L / N
TAU_WIN = 100.0
STEPS = int(round(TAU_WIN / dt))
M_LIST = [10, 15, 20, 25, 30]          # D5
H_DRIFT_TOL = 1e-4

disc = np.sqrt(b_par * b_par - 4.0 * a * c)
s_star = (b_par + disc) / (2.0 * c)
Vpp_star = a - 2.0 * b_par * s_star + 3.0 * c * s_star ** 2

coords = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")


def lap(u):
    return (np.roll(u, 1, 0) + np.roll(u, -1, 0)
            + np.roll(u, 1, 1) + np.roll(u, -1, 1) - 4.0 * u) / (dx * dx)


def H_u(u, ut, vpp):
    gx = (np.roll(u, -1, 0) - u) / dx
    gy = (np.roll(u, -1, 1) - u) / dx
    return float(np.sum(0.5 * ut * ut + 0.5 * kappa * (gx * gx + gy * gy)
                        + 0.5 * vpp * u * u) * dx * dx)


def run_wave(u0, ut0, vpp, steps, probe=(10, 10)):
    """Leapfrog D1; zwraca probke u(punkt)(tau), serie H_u, flagi."""
    u_prev = u0 - dt * ut0 + 0.5 * dt * dt * (kappa * lap(u0) - vpp * u0)
    u = u0.copy()
    series, Hs, taus = [], [], []
    flags = {"nonfinite": False}
    for stp in range(1, steps + 1):
        u_next = 2.0 * u - u_prev + dt * dt * (kappa * lap(u) - vpp * u)
        if stp % 10 == 0:
            if not np.all(np.isfinite(u_next)):
                flags["nonfinite"] = True
                break
            ut = (u_next - u_prev) / (2.0 * dt)   # D1: roznica CENTRALNA w u(stp)
            series.append(float(u[probe]))
            Hs.append(H_u(u, ut, vpp))
            taus.append(stp * dt)
        u_prev, u = u, u_next
    return np.array(series), np.array(Hs), np.array(taus), flags, u


def freq_zero_cross(series, taus):
    """omega z przejsc przez zero: 2*pi * (n_okresow) / czas."""
    sgn = np.sign(series)
    idx = np.where(sgn[:-1] * sgn[1:] < 0)[0]
    if len(idx) < 4:
        return np.nan
    # czas pierwszego i ostatniego przejscia (interpolacja liniowa)
    def tcross(i):
        return taus[i] + (taus[i + 1] - taus[i]) * series[i] / (series[i] - series[i + 1])
    t0, t1 = tcross(idx[0]), tcross(idx[-1])
    n_half = len(idx) - 1
    return np.pi * n_half / (t1 - t0)


print("=== Phase 1: kalibracja dyspersji (LOCK: Phase0_balance.md, D5) ===")
print(f"N={N} L={L} dx={dx} kappa={kappa} dt={dt} okno tau={TAU_WIN}")
print()

results = {}
disp_ok = True
g1_ok = True
det_ok = True
for tag, vpp_val in (("FV", a), ("TV", Vpp_star)):
    print(f"--- tlo {tag}: V'' = {vpp_val:.5f} (m = {np.sqrt(vpp_val):.4f}) ---")
    print("    m     k       omega_meas   omega_teor   blad%   omega_siec   "
          "H_dryft(sek.)  H_osc")
    for m in M_LIST:
        k = 2.0 * np.pi * m / L
        u0 = 1e-3 * np.cos(k * X)
        ut0 = np.zeros_like(u0)
        s1, Hs, taus, fl, u_fin = run_wave(u0, ut0, vpp_val, STEPS)
        # determinizm bitowy: powtorka pierwszego przypadku
        if m == M_LIST[0]:
            s2, _, _, _, u_fin2 = run_wave(u0, ut0, vpp_val, STEPS)
            det_ok = det_ok and bool(np.array_equal(u_fin, u_fin2))
        om_meas = freq_zero_cross(s1, taus)
        om_th = np.sqrt(kappa * k * k + vpp_val)
        k2_lat = (2 - 2 * np.cos(k * dx)) / (dx * dx)
        om_lat = np.sqrt(kappa * k2_lat + vpp_val)
        err = 100 * (om_meas - om_th) / om_th
        q = len(Hs) // 4
        drift = abs(np.mean(Hs[-q:]) - np.mean(Hs[:q])) / abs(Hs[0])
        osc = (Hs.max() - Hs.min()) / abs(Hs[0])
        disp_ok = disp_ok and abs(err) <= 5.0
        g1_ok = g1_ok and drift <= H_DRIFT_TOL and not fl["nonfinite"]
        results[(tag, m)] = (om_meas, om_th, err, drift)
        print(f"   {m:3d}  {k:.4f}   {om_meas:9.5f}   {om_th:9.5f}   "
              f"{err:+5.2f}   {om_lat:9.5f}   {drift:9.2e}   {osc:.2e}")
    print()

# przypadek nieprzepuszczalny (kontrola): omega ponizej m_FV nie propaguje
# (pomiar jakosciowy: pakiet o k=0 -> oscylacja lokalna, bez transportu)
print("=== EWALUACJA (LOCKED) ===")
drifts = [results[key][3] for key in results]
g1_all = g1_ok and det_ok
print(f"G1 (czesc Phase1): determinizm bitowy: {det_ok}; brak NaN; "
      f"max dryft sekularny H_u = {max(drifts):.2e} <= {H_DRIFT_TOL}: "
      f"{g1_ok} (oscylacja leapfrog O((dt*omega)^2) raportowana, poza "
      f"kryterium; D4)  -> {'PASS' if g1_all else 'FAIL'}")
errs = [abs(results[key][2]) for key in results]
g2 = disp_ok
print(f"G2 (kalibracja dyspersji): max |blad| = {max(errs):.2f}% <= 5%: "
      f"{disp_ok}  -> {'PASS' if g2 else 'FAIL'}")
print()
print(f"PODSUMOWANIE PHASE 1: G1(czesc)={'PASS' if g1_all else 'FAIL'} "
      f"G2={'PASS' if g2 else 'FAIL'}")
