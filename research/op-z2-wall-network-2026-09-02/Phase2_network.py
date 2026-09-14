#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-z2-wall-network (Phase 2) -- siec scian Z2 z genezy wielodomenowej.
Model substratu VERBATIM (lancuch: op-bare-substrate-genesis LOCK
sec. 1-2): s(x) 2D, Phi=s^2, ds/dtau = kappa Lap(s) - V'(s),
V = 0.5a s^2 - (b/3)|s|^3 + 0.25c s^4, a=0.5 b=1.6 c=1.0, kappa=0.5,
N=128 L=64 dx=0.5 dt=0.02, torus, Laplasjan 5-pkt, eps=0.30.

LOCK sec. 1-2: start GRF pasmowy grid-niezalezny (mody |n|<=16,
obwiednia exp(-(|n|/8)^2), normalizacja std z budowy N=256),
A in {0.6, 1.0} x seedy {20260906..20260910}; steps=100000 (tau=2000),
probkowanie co 100 krokow. Estymatory: L_wall = dx * #krawedzi
s_i*s_j<0 (dlugosc manhattanska -- Amendment A1); N_dom (label +
sklejanie periodyczne); persistence_tail = X(100000)/X(90000).

Uzycie: run | control <seed> <A> | verdict
"""
import json
import os
import sys
import time
import numpy as np
from scipy import ndimage

A_P, B_P, C_P = 0.5, 1.6, 1.0
KAPPA = 0.5
EPS = 0.30
S_STAR = (B_P + np.sqrt(B_P ** 2 - 4 * A_P * C_P)) / (2 * C_P)
S_BAR = (B_P - np.sqrt(B_P ** 2 - 4 * A_P * C_P)) / (2 * C_P)
PHI_STAR = S_STAR ** 2
L_BOX = 64.0
STEPS = 100000
SAMPLE = 100
SEEDS = (20260906, 20260907, 20260908, 20260909, 20260910)
AMPS = (0.6, 1.0)
NMODE = 16
ENV_W = 8.0
FITWIN = (100.0, 1000.0)
BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-z2-wall-network-2026-09-02/")
RESDIR = BASE + "Phase2_results/"
t0w = time.time()


def stamp(m):
    print("[t=%7.1fs] %s" % (time.time() - t0w, m), flush=True)


def banner(extra=""):
    print("REJESTR [INPUT]: model VERBATIM op-bare-substrate-genesis "
          "(a=0.5 b=1.6 c=1.0 kappa=0.5 eps=0.30); GRF |n|<=16 "
          "env exp(-(|n|/8)^2) std z N=256; A=%s seedy=%s; steps=%d; "
          "L_wall=manhattan (A1); okno fitu tau=[100,1000]%s"
          % (AMPS, SEEDS, STEPS, extra), flush=True)


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


# ---------------------------------------------------- GRF (FROZEN, A1)
def grf_build(seed, N):
    rng = np.random.default_rng(seed)
    M = 2 * NMODE + 1
    C = rng.standard_normal((M, M, 2))
    Cc = C[..., 0] + 1j * C[..., 1]
    Csym = 0.5 * (Cc + np.conj(Cc[::-1, ::-1]))
    n = np.arange(-NMODE, NMODE + 1)
    env = np.exp(-(n[:, None] ** 2 + n[None, :] ** 2) / ENV_W ** 2)
    F = np.zeros((N, N), dtype=complex)
    for i in range(M):
        for j in range(M):
            F[(i - NMODE) % N, (j - NMODE) % N] = Csym[i, j] * env[i, j]
    return np.real(np.fft.ifft2(F)) * N ** 2


def grf(seed, N):
    f256 = grf_build(seed, 256)
    sd = float(np.std(f256))
    if N == 256:
        return f256 / sd
    return grf_build(seed, N) / sd


# ------------------------------------------------------- estymatory
def wall_length(s, dx):
    e = 0
    for ax in (0, 1):
        e += int(np.sum(s * np.roll(s, -1, ax) < 0))
    return e * dx


def label_periodic(mask):
    lab, n = ndimage.label(mask)
    if n == 0:
        return 0
    parent = np.arange(n + 1)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for ax in (0, 1):
        a = np.take(lab, 0, axis=ax).ravel()
        b = np.take(lab, -1, axis=ax).ravel()
        for x, y in set(zip(a.tolist(), b.tolist())):
            if x > 0 and y > 0:
                rx, ry = find(x), find(y)
                if rx != ry:
                    parent[max(rx, ry)] = min(rx, ry)
    return len(set(find(i) for i in range(1, n + 1)))


def observables(s, dx):
    Phi = s * s
    metric = Phi >= EPS
    plus = metric & (s > 0)
    minus = metric & (s < 0)
    return dict(plus=float(np.mean(plus)), minus=float(np.mean(minus)),
                bare=float(np.mean(~metric)),
                Lw=wall_length(s, dx),
                ndp=label_periodic(plus), ndm=label_periodic(minus))


# ------------------------------------------------------------- biegi
def run_gen(seed, A, N=128, dx=0.5, dt=0.02, steps=STEPS, sample=SAMPLE,
            tag=""):
    eng = Sub2D(N, dx, dt)
    s = A * grf(seed, N)
    ser = []
    Hp, mono = None, True
    for k in range(steps + 1):
        if k > 0:
            s = eng.step(s)
        if k % sample == 0:
            ob = observables(s, dx)
            Hn = eng.H(s)
            if Hp is not None and Hn > Hp + 1e-10:
                mono = False
            Hp = Hn
            ob.update(k=k, tau=k * dt, H=Hn)
            ser.append(ob)
    i90 = int(round(0.9 * steps / sample))
    last, s90 = ser[-1], ser[i90]

    def tail(key):
        return last[key] / s90[key] if s90[key] > 0 else \
            (1.0 if last[key] == 0 else float('inf'))
    res = dict(seed=seed, A=A, N=N, tag=tag, H_mono=mono,
               plus_end=last["plus"], minus_end=last["minus"],
               bare_end=last["bare"], Lw_end=last["Lw"],
               ndp_end=last["ndp"], ndm_end=last["ndm"],
               tail_Lw=tail("Lw"), tail_bare=tail("bare"), series=ser)
    stamp("[gen s=%d A=%.1f N=%d%s] koniec: +%.3f/-%.3f bare=%.4f "
          "Lw=%.1f (tailLw=%.4f tailbare=%.4f) Ndom=%d/%d Hmono=%s"
          % (seed, A, N, tag, last["plus"], last["minus"], last["bare"],
             last["Lw"], res["tail_Lw"], res["tail_bare"], last["ndp"],
             last["ndm"], mono))
    return res, s


def rid(seed, A, tag=""):
    return "gen_s%d_A%g%s" % (seed, A, tag)


def save(res, s, name):
    os.makedirs(RESDIR, exist_ok=True)
    ser = res.pop("series")
    with open(RESDIR + name + ".json", "w") as f:
        json.dump(dict(res, series=ser), f)
    res["series"] = ser
    np.savez_compressed(RESDIR + name + ".npz", s=s,
                        psi=s * s / PHI_STAR)


def load(name):
    p = RESDIR + name + ".json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def run_all():
    banner("; Phase2 run (10 biegow glownych)")
    for A in AMPS:
        for seed in SEEDS:
            res, s = run_gen(seed, A)
            save(res, s, rid(seed, A))


def run_control(seed, A):
    banner("; kontrola pinningu N=256")
    res, s = run_gen(seed, A, N=256, dx=0.25, dt=0.005,
                     steps=4 * STEPS, sample=4 * SAMPLE, tag="_ctrl")
    save(res, s, rid(seed, A, "_ctrl"))


# ------------------------------------------------------------- werdykt
def at_tau(ser, tau):
    return min(ser, key=lambda p: abs(p["tau"] - tau))


def classify(r):
    """Klasy Q-NET-C (LOCK sec. 2)."""
    mn = min(r["plus_end"], r["minus_end"])
    if mn < 0.02 and r["Lw_end"] < 0.05 * L_BOX:
        return "SINGLE-DOMAIN"
    if (mn >= 0.05 and 0.99 <= r["tail_Lw"] <= 1.01
            and r["tail_bare"] >= 0.99):
        return "WOUND-STRIPES"
    return "TRANSIENT"


def verdict():
    out = []

    def em(x=""):
        print(x, flush=True)
        out.append(x)

    em("=" * 78)
    em("WERDYKTY op-z2-wall-network (litera LOCKa sec. 2 + A1)")
    banner("; verdict")
    em("=" * 78)
    # ---- A
    em("Q-NET-A (mozaika; punkt kontrolny tau=100):")
    okA, failmode = True, []
    for A in AMPS:
        for seed in SEEDS:
            r = load(rid(seed, A))
            if r is None:
                em("  s=%d A=%.1f: BRAK" % (seed, A))
                if A == 1.0:
                    okA = False
                continue
            p100 = at_tau(r["series"], 100.0)
            met = p100["plus"] + p100["minus"]
            em("  s=%d A=%.1f: tau=100 +%.3f/-%.3f metric=%.3f "
               "bare=%.3f Lw=%.1f Hmono=%s"
               % (seed, A, p100["plus"], p100["minus"], met,
                  p100["bare"], p100["Lw"], r["H_mono"]))
            if A == 1.0:
                if met < 0.80:
                    okA = False
                    failmode.append((seed, "BARE"))
                elif min(p100["plus"], p100["minus"]) < 0.10:
                    okA = False
                    failmode.append((seed, "ONESIGN"))
    qa = "Q-NET-A-PASS-MOSAIC" if okA else \
        "Q-NET-A-FAIL (%s)" % failmode
    em("  => %s" % qa)
    # ---- B
    em("Q-NET-B (grubienie; okno tau=[100,1000], A=1.0):")
    slopes = []
    for seed in SEEDS:
        r = load(rid(seed, 1.0))
        if r is None:
            continue
        win = [p for p in r["series"]
               if FITWIN[0] <= p["tau"] <= FITWIN[1] and p["Lw"] > 0]
        full = [p for p in r["series"]
                if FITWIN[0] <= p["tau"] <= FITWIN[1]]
        if len(win) < len(full) or len(win) < 10:
            em("  s=%d: okno niewazne (siec znika przed tau=1000: "
               "Lw>0 w %d/%d probkach)" % (seed, len(win), len(full)))
            continue
        x = np.log([p["tau"] for p in win])
        y = np.log([p["Lw"] for p in win])
        Adm = np.vstack([x, np.ones_like(x)]).T
        coef, _, _, _ = np.linalg.lstsq(Adm, y, rcond=None)
        pred = Adm @ coef
        r2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
        slopes.append((seed, float(coef[0]), float(r2)))
        em("  s=%d: slope=%.3f R2=%.4f" % (seed, coef[0], r2))
    okcnt = sum(1 for _, sl, _ in slopes if -0.65 <= sl <= -0.35)
    if len(slopes) < 3:
        qb = ("Q-NET-B-INCONCLUSIVE (<3 waznych okien -- sieci znikaja"
              " przed tau=1000)")
    elif okcnt >= 3:
        qb = "Q-NET-B-PASS-COARSENING (slope -0.5+-0.15 w %d/%d)" \
            % (okcnt, len(slopes))
    else:
        qb = "Q-NET-B-FAIL (w tolerancji tylko %d/%d)" \
            % (okcnt, len(slopes))
    em("  => %s" % qb)
    # ---- C
    em("Q-NET-C (stan w tau_max=2000):")
    classes = {}
    for A in AMPS:
        for seed in SEEDS:
            r = load(rid(seed, A))
            if r is None:
                continue
            cl = classify(r)
            classes[(seed, A)] = cl
            em("  s=%d A=%.1f: %s (+%.3f/-%.3f Lw=%.1f tailLw=%.4f "
               "tailbare=%.4f Ndom=%d/%d)"
               % (seed, A, cl, r["plus_end"], r["minus_end"],
                  r["Lw_end"], r["tail_Lw"], r["tail_bare"],
                  r["ndp_end"], r["ndm_end"]))
    wound = [s for s in SEEDS if classes.get((s, 1.0)) == "WOUND-STRIPES"]
    single = [s for s in SEEDS if classes.get((s, 1.0)) == "SINGLE-DOMAIN"]
    if wound:
        ctr_ok = []
        for s in wound:
            c = load(rid(s, 1.0, "_ctrl"))
            ctr_ok.append(c is not None and classify(c) == "WOUND-STRIPES")
            em("  kontrola N=256 s=%d: %s" % (s, "ZGODNA (WOUND-STRIPES)"
               if ctr_ok[-1] else "BRAK/NIEZGODNA"))
        qc = ("Q-NET-C-PASS-TOPOLOGICAL (WOUND-STRIPES: %s)" % wound) \
            if any(ctr_ok) else \
            "Q-NET-C-INCONCLUSIVE (kandydaci %s bez zgodnej kontroli)" \
            % wound
    elif len(single) == len(SEEDS):
        qc = "Q-NET-C-FAIL-COARSEN-AWAY (5/5 SINGLE-DOMAIN)"
    else:
        qc = ("Q-NET-C-INCONCLUSIVE (TRANSIENT-y: %s)"
              % [s for s in SEEDS
                 if classes.get((s, 1.0)) == "TRANSIENT"])
    em("  => %s" % qc)
    # ---- widok poziomu 1 (deskryptywnie)
    em()
    em("WIDOK POZIOMU 1 (deskryptywnie, obowiazkowe; psi=Phi/Phi*,"
       " Phi*=%.4f):" % PHI_STAR)
    for A in AMPS:
        for seed in SEEDS:
            r = load(rid(seed, A))
            if r is None:
                continue
            met = r["plus_end"] + r["minus_end"]
            dsc = (met * L_BOX ** 2 / r["Lw_end"]) if r["Lw_end"] > 0 \
                else float('inf')
            em("  s=%d A=%.1f: frakcja przedmetryczna (ukryta dla psi) "
               "= %.4f; skala mozaiki L^2*metric/Lw = %.1f"
               % (seed, A, r["bare_end"], dsc))
    with open(BASE + "Phase2_output.txt", "w") as f:
        f.write("\n".join(out) + "\n")
    print("zapisano:", BASE + "Phase2_output.txt")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "run":
        run_all()
    elif mode == "control":
        run_control(int(sys.argv[2]), float(sys.argv[3]))
    elif mode == "verdict":
        verdict()
    else:
        raise SystemExit("nieznany tryb: " + mode)
