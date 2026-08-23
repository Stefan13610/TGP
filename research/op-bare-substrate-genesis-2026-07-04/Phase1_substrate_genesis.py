# Phase 1: substrate genesis - scenariusze bare / single(A0) / cluster(A0,D)
# Implementacja SCISLE wg Phase0_balance.md (LOCK, wraz z poprawkami pre-code z sekcji 11).
#
# Dynamika: ds/dtau = kappa*Lap(s) - V'(s)   (przeplyw gradientowy, jawny Euler)
# V(s) = 0.5*a*s^2 - (b/3)*|s|^3 + 0.25*c*s^4,  Phi = s^2
# BRAK clampu / mobility / reguly amplitud. Guard |s|>S_GUARD tylko FLAGUJE runaway.

import numpy as np

# ---------------- parametry LOCKED (sekcja 2) ----------------
N = 128
L = 64.0
dx = L / N                     # 0.5
kappa = 0.50
a, b, c = 0.50, 1.60, 1.00
dt = 0.02
STEPS_DECAY = 6000
w = 1.50
A0_SCAN = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
D_SCAN = [2.0, 3.0, 4.0]
EPS = 0.30                     # Phi_metric_threshold
A_MIN = 4.0 / (N * N)          # minimalna metric_area (frakcja)
NOISE_AMP = 0.05               # bare: U(-0.05, +0.05)
RNG_SEED = 20260704
S_GUARD = 10.0                 # guard przeciw overflow (flaga runaway, NIE mechanizm locku)
FAR_CHEB = 0.4 * L             # strefa localized/boundary_contact (Chebyshev od centrum)
TAU_MID = 3000                 # mianownik persistence
TAIL_FROM = 5400               # persistence_tail = area(6000)/area(5400)
H_TOL_REL = 1e-10              # tolerancja numeryczna monotonicznosci H

coords = (np.arange(N) - N // 2) * dx
X, Y = np.meshgrid(coords, coords, indexing="ij")
FAR_ZONE = np.maximum(np.abs(X), np.abs(Y)) > FAR_CHEB

RECORD_AT = (0, 1500, 3000, 4500, 5400, 6000)


def Vpot(s):
    return 0.5 * a * s**2 - (b / 3.0) * np.abs(s) ** 3 + 0.25 * c * s**4


def Vprime(s):
    return s * (a - b * np.abs(s) + c * s * s)


def lap(s):  # 5-punktowy, periodyczny (LOCK sekcja 2)
    return (np.roll(s, 1, 0) + np.roll(s, -1, 0)
            + np.roll(s, 1, 1) + np.roll(s, -1, 1) - 4.0 * s) / (dx * dx)


def energy(s):
    gx = (np.roll(s, -1, 0) - s) / dx
    gy = (np.roll(s, -1, 1) - s) / dx
    return float(np.sum(0.5 * kappa * (gx * gx + gy * gy) + Vpot(s)) * dx * dx)


def gaussian(A0, cx, cy):
    return A0 * np.exp(-(((X - cx) ** 2 + (Y - cy) ** 2) / (2.0 * w * w)))


def make_bare():
    rng = np.random.default_rng(RNG_SEED)
    return rng.uniform(-NOISE_AMP, NOISE_AMP, size=(N, N))


def make_single(A0):
    return gaussian(A0, 0.0, 0.0)


def make_cluster(A0, D):
    s = gaussian(A0, 0.0, 0.0)
    for k in range(6):
        ang = k * np.pi / 3.0
        s = s + gaussian(A0, D * np.cos(ang), D * np.sin(ang))
    return s


def run(s0, steps, label):
    s = s0.copy()
    H0 = energy(s)
    H_prev = H0
    H_inc_max = 0.0
    runaway = False
    nonfinite = False
    contact_step = None
    rec = {}

    def observe(step):
        Phi = s * s
        rec[step] = (float(Phi.max()), float(np.mean(Phi > EPS)))

    observe(0)
    if np.any((s * s)[FAR_ZONE] > EPS):
        contact_step = 0

    for step in range(1, steps + 1):
        s = s + dt * (kappa * lap(s) - Vprime(s))
        if not np.all(np.isfinite(s)):
            nonfinite = True
            observe(step)
            break
        if np.max(np.abs(s)) > S_GUARD:
            runaway = True                      # forbidden move 3: dobicie do guardu = runaway
            observe(step)
            break
        H = energy(s)
        if H > H_prev:
            H_inc_max = max(H_inc_max, H - H_prev)
        H_prev = H
        if step in RECORD_AT:
            observe(step)
        if contact_step is None and step % 100 == 0:
            if np.any((s * s)[FAR_ZONE] > EPS):
                contact_step = step

    area_mid = rec.get(TAU_MID, (np.nan, np.nan))[1]
    area_tailfrom = rec.get(TAIL_FROM, (np.nan, np.nan))[1]
    area_fin = rec[max(rec)][1]
    persistence = (area_fin / area_mid) if area_mid and area_mid > 0 else float("nan")
    tail = (area_fin / area_tailfrom) if area_tailfrom and area_tailfrom > 0 else float("nan")
    Phi_fin = rec[max(rec)][0]
    return {
        "label": label,
        "s_final": s,
        "rec": rec,
        "H0": H0,
        "H_final": H_prev,
        "H_inc_max": H_inc_max,
        "H_monotone": H_inc_max <= H_TOL_REL * max(1.0, abs(H0)),
        "runaway": runaway,
        "nonfinite": nonfinite,
        "contact_step": contact_step,
        "boundary_contact": contact_step is not None,
        "area_final": area_fin,
        "Phi_max_final": Phi_fin,
        "persistence": persistence,
        "persistence_tail": tail,
        "decayed": area_fin < A_MIN,
        "localized_final": not bool(np.any((s * s)[FAR_ZONE] > EPS)),
        "s_max0": float(np.max(np.abs(s0))),
        "area0": rec[0][1],
    }


def fmt(x, nd=4):
    if isinstance(x, float) and np.isnan(x):
        return "n/a"
    return f"{x:.{nd}f}"


def row(r):
    return (f"{r['label']:<22} s_max(0)={fmt(r['s_max0'],3):>6}  "
            f"Phi_max_fin={fmt(r['Phi_max_final']):>8}  area_fin={fmt(r['area_final'],6):>9}  "
            f"pers={fmt(r['persistence'],3):>6}  tail={fmt(r['persistence_tail'],3):>6}  "
            f"{'DECAYED' if r['decayed'] else 'LOCKED ':<7}  "
            f"loc={'T' if r['localized_final'] else 'F'}  bc={'T' if r['boundary_contact'] else 'F'}  "
            f"H_mono={'T' if r['H_monotone'] else 'F'}  run={'T' if r['runaway'] else 'F'}")


def series(r):
    out = []
    for step in sorted(r["rec"]):
        pm, ar = r["rec"][step]
        out.append(f"      tau_step={step:>5}: Phi_max={pm:8.4f}  metric_area={ar:.6f}")
    return "\n".join(out)


print("=== Phase 1: substrate genesis (LOCK: Phase0_balance.md) ===")
print(f"N={N} L={L} dx={dx} kappa={kappa} a={a} b={b} c={c} dt={dt} "
      f"steps={STEPS_DECAY} w={w} eps={EPS} A_min={A_MIN:.6f} seed={RNG_SEED}")
print(f"brzegi: periodyczne, Laplasjan 5-pkt; guard |s|>{S_GUARD} (flaga, nie clamp)")
print()

results = {}

# ---------------- G1: determinizm (bare x2, bitowo) ----------------
print("--- runy ---")
r_bare_1 = run(make_bare(), STEPS_DECAY, "bare#1")
r_bare_2 = run(make_bare(), STEPS_DECAY, "bare#2")
deterministic = bool(np.array_equal(r_bare_1["s_final"], r_bare_2["s_final"]))
results["bare"] = r_bare_1
print(row(r_bare_1))
print(row(r_bare_2) + f"   [determinizm bitowy vs #1: {'T' if deterministic else 'F'}]")

# ---------------- single(A0) scan ----------------
singles = {}
for A0 in A0_SCAN:
    r = run(make_single(A0), STEPS_DECAY, f"single(A0={A0})")
    singles[A0] = r
    print(row(r))

# ---------------- wyznaczenie A_c (G3) ----------------
valid = {A0: r for A0, r in singles.items() if not (r["runaway"] or r["nonfinite"])}
decayed_A0 = [A0 for A0, r in valid.items() if r["decayed"]]
locked_A0 = [A0 for A0, r in valid.items() if not r["decayed"]]
A_sub = max(decayed_A0) if decayed_A0 else None      # najwiekszy podkrytyczny
A_sup = min(locked_A0) if locked_A0 else None        # najmniejszy nadkrytyczny

# ---------------- cluster(A_sub, D) ----------------
clusters = {}
if A_sub is not None:
    for D in D_SCAN:
        r = run(make_cluster(A_sub, D), STEPS_DECAY, f"cluster(A0={A_sub},D={D:g})")
        clusters[D] = r
        # flaga wg LOCK G4: s_max(0) klastra >= najmniejszy NADkrytyczny single.
        # Gdy zaden single nie lockuje (A_sup=None), komparator nie istnieje -> n/a.
        supflag = (A_sup is not None) and (r["s_max0"] >= A_sup)
        r["superposition_trivial"] = supflag
        sup_txt = ("n/a (brak nadkryt. singla)" if A_sup is None
                   else ("T" if supflag else "F"))
        r["sup_txt"] = sup_txt
        print(row(r) + f"   [superposition-trivial: {sup_txt}]")
else:
    print("BRAK podkrytycznego A0 w skanie -> cluster pominiety (G3 FAIL).")

print()
print("--- przebiegi (Phi_max, metric_area) ---")
for r in [r_bare_1] + [singles[A0] for A0 in A0_SCAN] + [clusters[D] for D in sorted(clusters)]:
    print(f"  {r['label']}:")
    print(series(r))
print()

# ---------------- ewaluacja G1-G4 (+G5 wstepnie) ----------------
print("=== EWALUACJA KRYTERIOW (LOCKED) ===")

all_runs = [r_bare_1, r_bare_2] + list(singles.values()) + list(clusters.values())
g1_finite = all(not (r["nonfinite"] or r["runaway"]) for r in all_runs)
g1_hmono = all(r["H_monotone"] for r in all_runs)
g1 = g1_finite and g1_hmono and deterministic
print(f"G1 (techniczne): finite/no-runaway={g1_finite}, H_nierosnace={g1_hmono}, "
      f"determinizm={deterministic}  -> {'PASS' if g1 else 'FAIL'}")

g2 = r_bare_1["area_final"] < A_MIN
print(f"G2 (bare pozostaje niemetryczne): metric_area_fin={r_bare_1['area_final']:.6f} "
      f"< A_min={A_MIN:.6f}: {'PASS' if g2 else 'FAIL'}")

# G3 wg LOCK: "Istnieje A0 (podkrytyczny), przy ktorym single(A0) po steps_decay
# ma metric_area < A_min". Falsyfikator: KAZDY A0 albo od razu lockuje, albo
# eksploduje. Brak nadkrytycznego singla w skanie NIE jest falsyfikatorem G3 -
# raportowany jako osobna obserwacja (A_c poza zakresem skanu).
g3 = len(decayed_A0) > 0
print(f"G3 (single decays): podkrytyczne A0={decayed_A0}, nadkrytyczne A0={locked_A0}")
if A_sup is not None:
    print(f"   A_c operacyjnie w przedziale ({A_sub}, {A_sup})  -> {'PASS' if g3 else 'FAIL'}")
else:
    print(f"   -> {'PASS' if g3 else 'FAIL (brak zanikajacego A0)'}")
    print(f"   OBSERWACJA (poza kryterium): zaden single z zalockowanego skanu nie lockuje;")
    print(f"   A_c > {max(A0_SCAN)} dla profilu w={w} -> doprecyzowanie progu w Phase 2.")

g4_pass_D = []
for D, r in sorted(clusters.items()):
    ok = ((not r["runaway"]) and (not r["nonfinite"])
          and r["area_final"] >= A_MIN
          and (not np.isnan(r["persistence"]) and r["persistence"] >= 0.9)
          and (not np.isnan(r["persistence_tail"]) and r["persistence_tail"] >= 0.99))
    if ok:
        g4_pass_D.append(D)
    print(f"G4 cluster(A0={A_sub}, D={D:g}): area>=A_min={r['area_final']>=A_MIN}, "
          f"pers={fmt(r['persistence'],3)}(>=0.9), tail={fmt(r['persistence_tail'],3)}(>=0.99), "
          f"superposition-trivial={r.get('sup_txt','F')}, "
          f"boundary_contact={'T' if r['boundary_contact'] else 'F'} -> {'LOCK' if ok else 'ZANIK/FAIL'}")
g4 = len(g4_pass_D) > 0
nontrivial_lock = [D for D in g4_pass_D if not clusters[D].get("superposition_trivial")]
print(f"G4 (kolektywny lock, RDZEN): {'PASS' if g4 else 'FAIL'}"
      + (f"  (lock przy D={g4_pass_D}; lock NIE-superpozycyjny przy D={nontrivial_lock if nontrivial_lock else 'BRAK'})"
         if g4 else "  -> klaster zanika jak single"))

print()
print("G5: wstepnie - ostrosc progu po A0 widoczna w skanie single; pelny skan progu,")
print("    przypadki graniczne i kontrola pinningu (N=256) -> Phase 2.")
print("G6: DeltaE_create (FAR/CORE/FRONT) -> Phase 2.")
print()
print(f"PODSUMOWANIE PHASE 1: G1={'PASS' if g1 else 'FAIL'} G2={'PASS' if g2 else 'FAIL'} "
      f"G3={'PASS' if g3 else 'FAIL'} G4={'PASS' if g4 else 'FAIL'} (G5/G6 -> Phase 2)")
