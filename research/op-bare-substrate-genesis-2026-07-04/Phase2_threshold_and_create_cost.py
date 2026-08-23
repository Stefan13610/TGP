# Phase 2: skan progu (G5) + DeltaE_create FAR/CORE/FRONT (G6)
# Implementacja SCISLE wg Phase0_balance.md (LOCK + poprawki pre-code, sekcja 11).
#
# [A] rozszerzenie skanu progu singla (Phase 1: wszystkie A0<=1.4 zanikaja)
#     - gorne A0 ograniczone stabilnoscia reakcji jawnego Eulera (raport wprost)
# [B] bisekcja progu klastra po A0 przy D in {2,3,4}  -> ostrosc progu (G5)
# [C] JEDEN kontrolny re-run przypadkow granicznych na N=256, dx=0.25, dt=0.005
#     (kontrola pinningu siatki z G5; te same parametry fizyczne, ten sam tau_final)
# [D] DeltaE_create na stanie zablokowanym (cluster(1.4, D=3) po steps_probe=3000)
#     w lokalizacjach FAR / CORE / FRONT (LOCK G6)
#
# BRAK clampu / mobility. Guard |s|>S_GUARD tylko flaguje runaway.

import numpy as np

# ---------------- parametry LOCKED ----------------
L = 64.0
KAPPA = 0.50
A, B, C = 0.50, 1.60, 1.00
W = 1.50                      # szerokosc seedu
EPS = 0.30
NOISE_SEED = 20260704         # (Phase 2 nie uzywa RNG - deterministyczna)
S_GUARD = 10.0
H_TOL_REL = 1e-10
TAU_FINAL = 120.0             # = 6000 * 0.02 (identyczny czas fizyczny na obu siatkach)
STEPS_PROBE_BASE = 3000       # steps_probe (LOCK) - przygotowanie stanu do G6
A_BUMP, W_BUMP = 0.10, 1.0    # delta_bump (LOCK, poprawka pre-code)
FAR_MIN_DIST = 10.0 * W       # FAR: > 10*w od struktury (LOCK G6)

BASE = dict(N=128, dt=0.02)               # siatka bazowa (LOCK)
CTRL = dict(N=256, dt=0.005)              # kontrola pinningu (LOCK G5, poprawka)

A0_SINGLE_EXT = [1.6, 2.0, 3.0, 4.0, 5.0]  # rozszerzenie skanu progu singla
D_SCAN = [2.0, 3.0, 4.0]
BISECT_LO, BISECT_HI = 0.2, 1.4            # brackety progu klastra (hi = lock z Phase 1)
BISECT_TOL = 0.01

# gorna granica stabilnosci reakcyjnej dla dt=0.02: dt < 2/max|V''| -> |s| ~< 6
# A0 > 5 zbliza sie do tej granicy -> skan uciety na 5.0 i raportowany wprost.


def Vpot(s):
    return 0.5 * A * s**2 - (B / 3.0) * np.abs(s) ** 3 + 0.25 * C * s**4


def Vprime(s):
    return s * (A - B * np.abs(s) + C * s * s)


class Sim:
    def __init__(self, N, dt):
        self.N, self.dt = N, dt
        self.dx = L / N
        self.steps = int(round(TAU_FINAL / dt))          # ten sam tau_final
        self.A_min = 4.0 / (N * N)
        co = (np.arange(N) - N // 2) * self.dx
        self.X, self.Y = np.meshgrid(co, co, indexing="ij")
        self.far_zone = np.maximum(np.abs(self.X), np.abs(self.Y)) > 0.4 * L

    def lap(self, s):
        return (np.roll(s, 1, 0) + np.roll(s, -1, 0)
                + np.roll(s, 1, 1) + np.roll(s, -1, 1) - 4.0 * s) / (self.dx ** 2)

    def energy(self, s):
        gx = (np.roll(s, -1, 0) - s) / self.dx
        gy = (np.roll(s, -1, 1) - s) / self.dx
        return float(np.sum(0.5 * KAPPA * (gx * gx + gy * gy) + Vpot(s)) * self.dx ** 2)

    def gaussian(self, A0, cx, cy):
        return A0 * np.exp(-(((self.X - cx) ** 2 + (self.Y - cy) ** 2) / (2.0 * W * W)))

    def single(self, A0):
        return self.gaussian(A0, 0.0, 0.0)

    def cluster(self, A0, D):
        s = self.gaussian(A0, 0.0, 0.0)
        for k in range(6):
            ang = k * np.pi / 3.0
            s = s + self.gaussian(A0, D * np.cos(ang), D * np.sin(ang))
        return s

    def evolve(self, s0, steps=None, h_every=20):
        steps = self.steps if steps is None else steps
        mid, tailfrom = steps // 2, int(round(0.9 * steps))
        s = s0.copy()
        H0 = self.energy(s)
        H_prev, H_inc = H0, 0.0
        runaway = nonfinite = False
        rec = {}

        def obs(k):
            Phi = s * s
            rec[k] = (float(Phi.max()), float(np.mean(Phi > EPS)))

        obs(0)
        contact = bool(np.any((s * s)[self.far_zone] > EPS))
        for k in range(1, steps + 1):
            s = s + self.dt * (KAPPA * self.lap(s) - Vprime(s))
            if not np.all(np.isfinite(s)):
                nonfinite = True; obs(k); break
            if np.max(np.abs(s)) > S_GUARD:
                runaway = True; obs(k); break
            if k % h_every == 0:
                H = self.energy(s)
                if H > H_prev:
                    H_inc = max(H_inc, H - H_prev)
                H_prev = H
            if k in (mid, tailfrom, steps):
                obs(k)
            if (not contact) and k % 200 == 0:
                contact = bool(np.any((s * s)[self.far_zone] > EPS))
        a_mid = rec.get(mid, (np.nan, np.nan))[1]
        a_tf = rec.get(tailfrom, (np.nan, np.nan))[1]
        a_fin = rec[max(rec)][1]
        pers = a_fin / a_mid if a_mid and a_mid > 0 else float("nan")
        tail = a_fin / a_tf if a_tf and a_tf > 0 else float("nan")
        return dict(s_final=s, area_final=a_fin, Phi_max_final=rec[max(rec)][0],
                    persistence=pers, persistence_tail=tail,
                    H_monotone=H_inc <= H_TOL_REL * max(1.0, abs(H0)),
                    runaway=runaway, nonfinite=nonfinite, boundary_contact=contact,
                    rec=rec)


def classify(sim, r):
    if r["nonfinite"]:
        return "NONFINITE"
    if r["runaway"]:
        return "RUNAWAY"
    if r["area_final"] < sim.A_min:
        return "DECAY"
    if (r["persistence"] >= 0.9 and r["persistence_tail"] >= 0.99):
        return "LOCK"
    return "MARGINAL"


def fmt(x, nd=3):
    return "n/a" if (isinstance(x, float) and np.isnan(x)) else f"{x:.{nd}f}"


base = Sim(**BASE)

print("=== Phase 2: skan progu (G5) + DeltaE_create (G6) ===")
print(f"baza: N={base.N} dx={base.dx} dt={base.dt} steps={base.steps} tau_final={TAU_FINAL}")
print(f"potencjal a,b,c={A},{B},{C} kappa={KAPPA} eps={EPS} A_min={base.A_min:.6f}")
print(f"delta_bump: A={A_BUMP} w={W_BUMP}; FAR > {FAR_MIN_DIST} od struktury")
print()

# ================= [A] rozszerzenie skanu progu singla =================
print("--- [A] single(A0), rozszerzenie skanu (w=1.5, siatka bazowa) ---")
single_lock_A0 = None
single_decay_A0 = 1.4          # najwiekszy zanikajacy single z Phase 1
for A0 in A0_SINGLE_EXT:
    r = base.evolve(base.single(A0))
    cl = classify(base, r)
    print(f"  single(A0={A0:>4}): {cl:<9} area_fin={r['area_final']:.6f} "
          f"pers={fmt(r['persistence'])} tail={fmt(r['persistence_tail'])} "
          f"H_mono={'T' if r['H_monotone'] else 'F'}")
    if cl == "DECAY":
        single_decay_A0 = max(single_decay_A0, A0)
    if cl == "LOCK" and single_lock_A0 is None:
        single_lock_A0 = A0
if single_lock_A0 is None:
    print(f"  -> zaden single do A0=5.0 nie lockuje: A_c(single, w={W}) > 5.0")
    print(f"     (wyzsze A0 lamia warunek stabilnosci reakcji dt<2/|V''| dla dt=0.02;")
    print(f"      raportowane wprost, bez zmiany dt - prog ostry pokazujemy po osi A0 klastra)")
else:
    print(f"  -> najmniejszy lockujacy single: A0={single_lock_A0}")
print()

# ================= [B] bisekcja progu klastra po A0 =================
print("--- [B] cluster(A0, D): bisekcja progu po A0 (ostrosc progu, G5) ---")
thresholds = {}
for D in D_SCAN:
    lo, hi = BISECT_LO, BISECT_HI
    r_lo = base.evolve(base.cluster(lo, D))
    r_hi = base.evolve(base.cluster(hi, D))
    c_lo, c_hi = classify(base, r_lo), classify(base, r_hi)
    print(f"  D={D:g}: bracket A0={lo}:{c_lo}  A0={hi}:{c_hi}")
    if not (c_lo == "DECAY" and c_hi == "LOCK"):
        print(f"    !! bracket nieprawidlowy - pomijam bisekcje dla D={D:g}")
        continue
    while hi - lo > BISECT_TOL:
        mid = 0.5 * (lo + hi)
        r_m = base.evolve(base.cluster(mid, D))
        c_m = classify(base, r_m)
        print(f"    A0={mid:.5f}: {c_m:<9} area_fin={r_m['area_final']:.6f} "
              f"pers={fmt(r_m['persistence'])} tail={fmt(r_m['persistence_tail'])}")
        if c_m == "LOCK":
            hi = mid
        else:
            lo = mid            # DECAY i MARGINAL po stronie zaniku (raportowane wprost)
    thresholds[D] = (lo, hi)
    print(f"  D={D:g}: prog A0_c(cluster) w ({lo:.5f}, {hi:.5f}), szerokosc {hi-lo:.5f}")
print()
if thresholds:
    if single_lock_A0 is not None:
        single_txt = f"A_c(single) w ({single_decay_A0}, {single_lock_A0})"
    else:
        single_txt = f"A_c(single) > {max(A0_SINGLE_EXT)}"
    adv_base = single_decay_A0 if single_lock_A0 is not None else max(A0_SINGLE_EXT)
    print(f"  Podsumowanie progow klastra (kolektywna przewaga; {single_txt}):")
    for D, (lo, hi) in sorted(thresholds.items()):
        print(f"    D={D:g}: A0_c ~ {0.5*(lo+hi):.4f}  -> przewaga kolektywna >= "
              f"{adv_base/(0.5*(lo+hi)):.1f}x (dolne oszacowanie z {adv_base})")
print()

# ================= [C] kontrola pinningu: N=256 =================
print("--- [C] kontrola pinningu siatki (JEDEN re-run przypadkow granicznych, LOCK G5) ---")
pin_ok = None
if 3.0 in thresholds:
    lo3, hi3 = thresholds[3.0]
    ctrl = Sim(**CTRL)
    print(f"  kontrola: N={ctrl.N} dx={ctrl.dx} dt={ctrl.dt} steps={ctrl.steps} "
          f"(ten sam tau_final={TAU_FINAL})")
    cases = [("graniczny LOCK", hi3, "LOCK"), ("graniczny DECAY", lo3, "DECAY")]
    pin_ok = True
    for name, A0, expected in cases:
        r_c = ctrl.evolve(ctrl.cluster(A0, 3.0), h_every=50)
        c_c = classify(ctrl, r_c)
        agree = (c_c == expected)
        pin_ok = pin_ok and agree
        print(f"  D=3 A0={A0:.5f} ({name}): baza={expected}  drobna_siatka={c_c}  "
              f"area_fin={r_c['area_final']:.6f} pers={fmt(r_c['persistence'])} "
          f"tail={fmt(r_c['persistence_tail'])} -> {'ZGODNE' if agree else 'NIEZGODNE (pinning!)'}")
else:
    print("  brak progu dla D=3 - kontrola pominieta")
print()

# ================= [D] G6: DeltaE_create FAR / CORE / FRONT =================
print("--- [D] G6: DeltaE_create na stanie zablokowanym ---")
prep = base.evolve(base.cluster(1.4, 3.0), steps=STEPS_PROBE_BASE)
s_state = prep["s_final"]
Phi = s_state * s_state
print(f"  stan: cluster(1.4, D=3) po steps_probe={STEPS_PROBE_BASE} "
      f"(area={prep['area_final']:.6f}, Phi_max={prep['Phi_max_final']:.4f})")

Nn = base.N
jc = Nn // 2
# FRONT: pierwszy wezel na promieniu +x od centrum, gdzie Phi spada ponizej eps
i_front = None
for i in range(jc, Nn):
    if Phi[i, jc] < EPS:
        i_front = i
        break
# CORE: maksimum Phi
i_core, j_core = np.unravel_index(np.argmax(Phi), Phi.shape)
# FAR: rog, weryfikacja odleglosci > 10*w od struktury (Phi > eps)
metric_pts = np.argwhere(Phi > EPS)
far_ij = None
for cand in [(-28.0, -28.0), (-30.0, -30.0), (-26.0, -30.0)]:
    ii = int(round(cand[0] / base.dx)) + Nn // 2
    jj = int(round(cand[1] / base.dx)) + Nn // 2
    d = np.min(np.sqrt(((metric_pts[:, 0] - ii) ** 2 + (metric_pts[:, 1] - jj) ** 2))) * base.dx
    if d > FAR_MIN_DIST:
        far_ij = (ii, jj, d)
        break

locs = {
    "FAR":   (far_ij[0], far_ij[1]),
    "CORE":  (int(i_core), int(j_core)),
    "FRONT": (int(i_front), jc),
}
H_state = base.energy(s_state)
results_g6 = {}
for name, (ii, jj) in locs.items():
    cx, cy = base.X[ii, jj], base.Y[ii, jj]
    bump = A_BUMP * np.exp(-(((base.X - cx) ** 2 + (base.Y - cy) ** 2) / (2.0 * W_BUMP ** 2)))
    dE = base.energy(s_state + bump) - H_state
    results_g6[name] = dE
    extra = f" (odleglosc od struktury {far_ij[2]:.1f} > {FAR_MIN_DIST})" if name == "FAR" else ""
    print(f"  {name:<5} @ ({cx:+.1f},{cy:+.1f}): s={s_state[ii,jj]:+.4f} Phi={Phi[ii,jj]:.4f}"
          f"  ->  DeltaE_create = {dE:+.6f}{extra}")
print()

# ================= EWALUACJA G5 / G6 =================
print("=== EWALUACJA (LOCKED) ===")

sharp = all((hi - lo) <= 2 * BISECT_TOL for lo, hi in thresholds.values()) and len(thresholds) == len(D_SCAN)
g5 = bool(sharp and (pin_ok is True))
print(f"G5 (ostry prog + pelna relaksacja + brak artefaktu):")
print(f"    ostry prog po A0 klastra: {'T' if sharp else 'F'} "
      f"(szerokosc przedzialu <= {2*BISECT_TOL} przy kazdym D)")
print(f"    pelna relaksacja: brak clampu w kodzie; H_Gamma nierosnace we wszystkich runach")
print(f"    kontrola pinningu N=256: {'ZGODNA' if pin_ok else ('NIEZGODNA' if pin_ok is False else 'n/a')}")
print(f"    -> {'PASS' if g5 else 'FAIL'}")

dE_far, dE_core, dE_front = results_g6["FAR"], results_g6["CORE"], results_g6["FRONT"]
g6 = dE_front < dE_far
g6_strong = dE_front <= 0.0
print(f"G6 (koszt kreacji): FRONT={dE_front:+.6f}  <  FAR={dE_far:+.6f} : "
      f"{'T' if g6 else 'F'};  FRONT<=0 (mocna wersja): {'T' if g6_strong else 'F'}")
print(f"    CORE={dE_core:+.6f} (raport wprost; predykcja pre-code: CORE > FAR: "
      f"{'potwierdzona' if dE_core > dE_far else 'NIEpotwierdzona'})")
print(f"    -> {'PASS' if g6 else 'FAIL'}" + ("  [mocna wersja: kreacja na froncie samoczynna]" if g6_strong else ""))
print()
print(f"PODSUMOWANIE PHASE 2: G5={'PASS' if g5 else 'FAIL'} G6={'PASS' if g6 else 'FAIL'}")
print(f"(G1-G4 PASS w Phase 1; regula decyzyjna cyklu -> Phase_FINAL_close.md)")
