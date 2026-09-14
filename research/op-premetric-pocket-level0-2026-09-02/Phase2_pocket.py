#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-premetric-pocket-level0 (Phase 2) -- kieszen przedmetryczna na
poziomie 0. Model VERBATIM z op-bare-substrate-genesis (CLOSED, CYTAT
w LOCK sec. 1): s(x) 2D, Phi=s^2, ds/dtau = kappa Lap(s) - V'(s),
V = 0.5a s^2 - (b/3)|s|^3 + 0.25c s^4, a=0.5 b=1.6 c=1.0, kappa=0.5,
N=128 L=64 dx=0.5 dt=0.02, torus, Laplasjan 5-pkt, eps=0.30,
A_min = 4/128^2 (frakcja).

Scenariusze (LOCK sec. 1-2): A pocketA_R{4,8,16} (kieszen gola w +s*);
B dropB_R{4,8,16} (kropla antyfazowa -s*); C stripes (sciany Z2);
D pocketD_R{4,8} (kieszen na scianie). Kontrola pinningu: N=256
dx=0.25 dt=0.005 (C zawsze; pozytywy A/B/D).

REJESTR [INPUT]: seed=20260905; R in {4,8,16}; delta=1.07;
steps_max=30000 (tau_max=600); probkowanie co 50 krokow (dtau=1);
eps=0.30; s*=1.1741657 Phi*=1.3786652; s_bar=0.4258343 Phi_bar=0.1813348.

Uzycie: run | control <sid> | verdict
"""
import json
import os
import sys
import time
import numpy as np

A_P, B_P, C_P = 0.5, 1.6, 1.0
KAPPA = 0.5
EPS = 0.30                                     # dziedziczone
AMIN_FRAC = 4.0 / 128 ** 2                     # dziedziczone (frakcja)
S_STAR = (B_P + np.sqrt(B_P ** 2 - 4 * A_P * C_P)) / (2 * C_P)
S_BAR = (B_P - np.sqrt(B_P ** 2 - 4 * A_P * C_P)) / (2 * C_P)
PHI_STAR = S_STAR ** 2
DELTA = 1.07                                   # INPUT (szerokosc scianki)
SEED = 20260905                                # INPUT
L_BOX = 64.0
STEPS = 30000
SAMPLE = 50
RADII = (4.0, 8.0, 16.0)
BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-premetric-pocket-level0-2026-09-02/")
RESDIR = BASE + "Phase2_results/"
t0w = time.time()


def stamp(m):
    print("[t=%6.1fs] %s" % (time.time() - t0w, m), flush=True)


def banner(extra=""):
    print("REJESTR [INPUT]: model VERBATIM op-bare-substrate-genesis "
          "(a=0.5 b=1.6 c=1.0 kappa=0.5 eps=0.30 A_min=4/128^2); "
          "seed=%d; R=%s; delta=%.2f; s*=%.7f Phi*=%.7f s_bar=%.7f; "
          "steps=%d dtau_probe=1%s"
          % (SEED, RADII, DELTA, S_STAR, PHI_STAR, S_BAR, STEPS, extra),
          flush=True)


def Vpot(s):
    return 0.5 * A_P * s * s - (B_P / 3.0) * np.abs(s) ** 3 \
        + 0.25 * C_P * s ** 4


def Vprime(s):
    return s * (A_P - B_P * np.abs(s) + C_P * s * s)


class Sub2D:
    def __init__(self, N, dx, dt):
        self.N, self.dx, self.dt = N, dx, dt

    def lap(self, s):
        return (np.roll(s, 1, 0) + np.roll(s, -1, 0) + np.roll(s, 1, 1)
                + np.roll(s, -1, 1) - 4 * s) / self.dx ** 2

    def step(self, s):
        return s + self.dt * (KAPPA * self.lap(s) - Vprime(s))

    def H(self, s):
        e = np.sum(Vpot(s))
        for ax in (0, 1):
            d = (np.roll(s, -1, ax) - s) / self.dx
            e += 0.5 * KAPPA * np.sum(d * d)
        return float(e) * self.dx ** 2


def r_center(N, dx, cx=0.5, cy=0.5):
    x = (np.arange(N) + 0.5) * dx
    return np.sqrt((x[:, None] - cx * L_BOX) ** 2
                   + (x[None, :] - cy * L_BOX) ** 2)


def noise(N):
    rng = np.random.default_rng(SEED)
    f = rng.uniform(-0.05, 0.05, (128, 128))    # pole bazowe (dziedz.)
    if N == 128:
        return f
    return np.repeat(np.repeat(f, N // 128, 0), N // 128, 1)


def stripes_field(N, dx):
    x = (np.arange(N) + 0.5) * dx
    left = np.tanh((x - L_BOX / 4) / DELTA)
    right = np.tanh((3 * L_BOX / 4 - x) / DELTA)
    prof = np.where(x < L_BOX / 2, left, right)
    return S_STAR * prof[:, None] * np.ones((1, N))


def build(sid, N, dx):
    """Stan poczatkowy scenariusza (LOCK sec. 1)."""
    nz = noise(N)
    if sid.startswith("pocketA_R"):
        R = float(sid.split("R")[1])
        r = r_center(N, dx)
        w = 0.5 * (1 + np.tanh((r - R) / DELTA))
        return S_STAR * w + nz * (1 - w)
    if sid.startswith("dropB_R"):
        R = float(sid.split("R")[1])
        return S_STAR * np.tanh((r_center(N, dx) - R) / DELTA)
    if sid == "stripes":
        return stripes_field(N, dx)
    if sid.startswith("pocketD_R"):
        R = float(sid.split("R")[1])
        r = r_center(N, dx, cx=0.25, cy=0.5)
        w = 0.5 * (1 + np.tanh((r - R) / DELTA))
        return stripes_field(N, dx) * w + nz * (1 - w)
    raise SystemExit("nieznany scenariusz: " + sid)


SCEN = (["pocketA_R%g" % R for R in RADII]
        + ["dropB_R%g" % R for R in RADII]
        + ["stripes"] + ["pocketD_R%g" % R for R in (4.0, 8.0)])


def run_scen(sid, N=128, dx=0.5, dt=0.02, steps=STEPS, sample=SAMPLE,
             tag=""):
    eng = Sub2D(N, dx, dt)
    s = build(sid, N, dx)
    ser = []
    Hprev = None
    Hmono = True
    bcontact = False
    localized_obj = sid.startswith(("pocketA", "dropB"))
    for k in range(steps + 1):
        if k > 0:
            s = eng.step(s)
        if k % sample == 0:
            Phi = s * s
            bare = float(np.mean(Phi < EPS))
            Hn = eng.H(s)
            if Hprev is not None and Hn > Hprev + 1e-10:
                Hmono = False
            Hprev = Hn
            if localized_obj and bare > 0:
                mask = Phi < EPS
                idx = np.argwhere(mask)
                c = (N - 1) / 2.0
                cheb = float(np.max(np.max(np.abs(idx - c), axis=1))) * dx
                if cheb > 0.4 * L_BOX:
                    bcontact = True
            ser.append(dict(k=k, tau=k * dt, bare=bare,
                            metric=float(np.mean(Phi >= EPS)),
                            phimin=float(Phi.min()),
                            phimax=float(Phi.max()),
                            smin=float(s.min()), H=Hn))
    tail = ser[-1]["bare"] / ser[-4]["bare"] if ser[-4]["bare"] > 0 else \
        (1.0 if ser[-1]["bare"] == 0 else float('inf'))
    t_close = next((p["tau"] for p in ser if p["bare"] < AMIN_FRAC), None)
    t_life = next((p["tau"] for p in ser if -p["smin"] < S_BAR), None)
    res = dict(sid=sid, N=N, dt=dt, tag=tag, bare_end=ser[-1]["bare"],
               bare_end_nodes=ser[-1]["bare"] * N * N,
               metric_end=ser[-1]["metric"], phimin_end=ser[-1]["phimin"],
               smin_end=ser[-1]["smin"], tail=tail, t_close=t_close,
               t_life=t_life, H_mono=Hmono, boundary_contact=bcontact,
               series=ser)
    stamp("[%s%s N=%d] koniec: bare=%.3e (%.1f wezlow/128^2-ekw) "
          "tail=%.4f t_close=%s t_life=%s smin=%.3f Hmono=%s bc=%s"
          % (sid, tag, N, res["bare_end"],
             res["bare_end"] / AMIN_FRAC * 4, tail, t_close, t_life,
             res["smin_end"], Hmono, bcontact))
    return res, s


def save(res, s, name):
    os.makedirs(RESDIR, exist_ok=True)
    with open(RESDIR + name + ".json", "w") as f:
        json.dump(res, f)
    np.savez_compressed(RESDIR + name + ".npz", s=s)


def load(name):
    p = RESDIR + name + ".json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def run_all():
    banner("; Phase2 run")
    for sid in SCEN:
        res, s = run_scen(sid)
        save(res, s, sid)


def run_control(sid):
    banner("; kontrola pinningu N=256")
    res, s = run_scen(sid, N=256, dx=0.25, dt=0.005, steps=120000,
                      sample=200, tag="_ctrl")
    save(res, s, sid + "_ctrl")


def fit_scaling(Rs, taus):
    Rs, taus = np.array(Rs, float), np.array(taus, float)
    out = {}
    for name, X in (("lin", Rs), ("kwadrat", Rs ** 2)):
        A = np.vstack([X, np.ones_like(X)]).T
        coef, res_, _, _ = np.linalg.lstsq(A, taus, rcond=None)
        pred = A @ coef
        ss = 1 - np.sum((taus - pred) ** 2) / max(
            np.sum((taus - taus.mean()) ** 2), 1e-300)
        out[name] = (float(coef[0]), float(coef[1]), float(ss))
    return out


def verdict():
    out = []

    def em(x=""):
        print(x, flush=True)
        out.append(x)

    em("=" * 78)
    em("WERDYKTY op-premetric-pocket-level0 (litera LOCKa sec. 2)")
    banner("; verdict")
    em("=" * 78)
    # ---- A
    em("Q-PKT-A (kieszen gola):")
    heals, persist, inc = [], [], False
    for R in RADII:
        r = load("pocketA_R%g" % R)
        if r is None:
            inc = True
            continue
        em("  R=%-4g bare_end=%.3e t_close=%s tail=%.4f bc=%s Hmono=%s"
           % (R, r["bare_end"], r["t_close"], r["tail"],
              r["boundary_contact"], r["H_mono"]))
        if r["boundary_contact"] or not r["H_mono"]:
            inc = True
        elif r["t_close"] is not None and r["bare_end"] < AMIN_FRAC:
            heals.append((R, r["t_close"]))
        elif r["bare_end"] >= AMIN_FRAC and r["tail"] >= 0.99:
            persist.append(R)
        else:
            inc = True
    if persist:
        ctrls = [load("pocketA_R%g_ctrl" % R) for R in persist]
        ok = all(c is not None and c["bare_end"] >= AMIN_FRAC
                 and c["tail"] >= 0.99 for c in ctrls)
        qa = ("Q-PKT-A-PERSIST(%s)" % persist) if ok else \
            "Q-PKT-A-INCONCLUSIVE (kontrola N=256 niezgodna/brak)"
    elif len(heals) == len(RADII) and not inc:
        qa = "Q-PKT-A-HEAL-ALL"
        f = fit_scaling([h[0] for h in heals], [h[1] for h in heals])
        em("  tau_close(R): %s; fit lin: a=%.3f R2=%.5f; "
           "kwadrat: a=%.4f R2=%.5f"
           % (heals, f["lin"][0], f["lin"][2], f["kwadrat"][0],
              f["kwadrat"][2]))
    else:
        qa = "Q-PKT-A-INCONCLUSIVE"
    em("  => %s" % qa)
    # ---- B
    em("Q-PKT-B (kropla antyfazowa -- centralne):")
    lives, survive, incb = [], [], False
    for R in RADII:
        r = load("dropB_R%g" % R)
        if r is None:
            incb = True
            continue
        em("  R=%-4g t_life=%s bare_end=%.3e smin_end=%.3f tail=%.4f "
           "bc=%s Hmono=%s" % (R, r["t_life"], r["bare_end"],
                               r["smin_end"], r["tail"],
                               r["boundary_contact"], r["H_mono"]))
        if r["boundary_contact"] or not r["H_mono"]:
            incb = True
        elif r["t_life"] is not None:
            lives.append((R, r["t_life"]))
        elif -r["smin_end"] >= S_BAR and r["tail"] >= 0.99:
            survive.append(R)
        else:
            incb = True
    if survive:
        ctrls = [load("dropB_R%g_ctrl" % R) for R in survive]
        ok = all(c is not None and -c["smin_end"] >= S_BAR
                 and c["tail"] >= 0.99 for c in ctrls)
        qb = ("Q-PKT-B-TRWA(%s) -- PIERWSZY trwaly obiekt przedmetryczny"
              % survive) if ok else \
            "Q-PKT-B-INCONCLUSIVE (kontrola N=256 niezgodna/brak)"
    elif len(lives) == len(RADII) and not incb:
        qb = "Q-PKT-B-ZANIK"
        f = fit_scaling([x[0] for x in lives], [x[1] for x in lives])
        em("  tau_life(R): %s; fit lin: a=%.3f R2=%.5f; "
           "kwadrat: a=%.4f R2=%.5f (kandydat R^2 = curvature flow)"
           % (lives, f["lin"][0], f["lin"][2], f["kwadrat"][0],
              f["kwadrat"][2]))
    else:
        qb = "Q-PKT-B-INCONCLUSIVE"
    em("  => %s" % qb)
    # ---- C
    em("Q-PKT-C (sciany Z2 jako trwale struktury przedmetryczne):")
    rC = load("stripes")
    cC = load("stripes_ctrl")
    if rC is None or cC is None:
        qc = "Q-PKT-C-INCONCLUSIVE (brak biegu/kontroli)"
    else:
        em("  N=128: bare=%.4e (%.0f wezlow) tail=%.4f Phimin=%.3e "
           "Hmono=%s" % (rC["bare_end"], rC["bare_end"] * 128 * 128,
                         rC["tail"], rC["phimin_end"], rC["H_mono"]))
        em("  N=256: bare=%.4e tail=%.4f Phimin=%.3e"
           % (cC["bare_end"], cC["tail"], cC["phimin_end"]))
        ok = (rC["bare_end"] >= AMIN_FRAC and rC["tail"] >= 0.99
              and rC["H_mono"] and cC["bare_end"] >= AMIN_FRAC
              and cC["tail"] >= 0.99)
        qc = "Q-PKT-C-PASS-SHEET" if ok else "Q-PKT-C-FAIL"
    em("  => %s" % qc)
    # ---- D
    em("Q-PKT-D (pozostalosc po kieszeni na scianie):")
    qd_parts = []
    for R in (4.0, 8.0):
        r = load("pocketD_R%g" % R)
        if r is None or rC is None:
            qd_parts.append((R, "BRAK"))
            continue
        excess = r["bare_end"] - rC["bare_end"]
        em("  R=%-4g bare_end=%.4e baza(sciany)=%.4e excess=%.3e "
           "(prog A_min=%.3e) tail=%.4f"
           % (R, r["bare_end"], rC["bare_end"], excess, AMIN_FRAC,
              r["tail"]))
        if excess >= AMIN_FRAC and r["tail"] >= 0.99:
            c = load("pocketD_R%g_ctrl" % R)
            cok = (c is not None and cC is not None
                   and (c["bare_end"] - cC["bare_end"]) >= AMIN_FRAC
                   and c["tail"] >= 0.99)
            qd_parts.append((R, "RESIDUAL-OBJECT" if cok else
                             "INCONCLUSIVE (kontrola)"))
        else:
            qd_parts.append((R, "CLEAN-WALL"))
    em("  => Q-PKT-D: %s" % qd_parts)
    em()
    em("MAPOWANIE NA POZIOM 1 (deskryptywne, obowiazkowe): psi = Phi/Phi*"
       " (Phi*=%.4f); pas przedmetryczny Phi<%.2f <=> psi<%.3f"
       % (PHI_STAR, EPS, EPS / PHI_STAR))
    with open(BASE + "Phase2_output.txt", "w") as f:
        f.write("\n".join(out) + "\n")
    print("zapisano:", BASE + "Phase2_output.txt")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "run":
        run_all()
    elif mode == "control":
        run_control(sys.argv[2])
    elif mode == "verdict":
        verdict()
    else:
        raise SystemExit("nieznany tryb: " + mode)
