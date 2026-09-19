#!/usr/bin/env python3
# -*- coding: ascii -*-
"""op-matter-induced-creation -- wspolny driver ewolucji (MD sec. 6-8).

Zero decyzji merytorycznych: wszystkie progi/okna/detektory sa
argumentami wolanymi z Phase2/Phase3 zgodnie z MD (FROZEN).
Rejestrowane co dt_out=0.1: t, psi(0), min/max psi, E_core (gestosc
pelna), E_core^field, max_{r<=80}|psidot|, max_{r<=40}|psi-1|.
Profile pelne co prof_every=5 j.cz. Checkpoint npz co 100 j.cz.
"""
import os

import numpy as np

import engine_core as ec

RMAX_E = 80.0     # E_core r<=80 [LOCK]
RMAX_OBJ = 40.0   # max|psi-1| na r<=40 [LOCK]


def evolve(h, dt, lam, t0, t1, g, pi, sponge=True, ramp=True,
           dt_out=0.1, prof_every=5.0, ckpt_every=100.0,
           ckpt_path=None, keep_profiles=True):
    """Ewolucja od t0 do t1. ramp=True -> lam(t) wg LOCK sec.3
    (dla t<=600 zwraca dokladnie lam). Zwraca dict."""
    eng = ec.Engine(h, 200.0, sponge=sponge, lam=lam)
    n_out = int(round(dt_out/dt))
    n_prof = int(round(prof_every/dt))
    n_ckpt = int(round(ckpt_every/dt))
    nsteps = int(round((t1 - t0)/dt))
    n0 = int(round(t0/dt))

    ts, psi0s, mns, mxs, Es, Efs, Vs, Ds = [], [], [], [], [], [], [], []
    prof_t, prof = [], []

    def lam_t(t):
        return ec.lam_of_t(lam, t) if ramp else lam

    def record(n):
        t = n*dt
        eng.set_lam(lam_t(t))
        ts.append(t)
        psi0s.append(float(g[0]))
        mns.append(float(np.min(g)))
        mxs.append(float(np.max(g)))
        Es.append(eng.energy(g, pi, rmax=RMAX_E, matter=True))
        Efs.append(eng.energy(g, pi, rmax=RMAX_E, matter=False))
        pd = np.abs(pi/ec.Mfun(g))
        Vs.append(float(np.max(pd[eng.r <= RMAX_E])))
        Ds.append(float(np.max(np.abs(g[eng.r <= RMAX_OBJ] - 1.0))))

    record(n0)
    if keep_profiles:
        prof_t.append(t0)
        prof.append(g.copy())

    status = "OK"
    t_end = None
    for k in range(1, nsteps + 1):
        n = n0 + k
        t_cur = (n - 1)*dt
        t_new = n*dt
        eng.set_lam(lam_t(t_cur))
        try:
            g, pi = eng.step(g, pi, dt, lrho_new=eng.lrho_of(lam_t(t_new)))
        except ec.NonConvergence:
            status = "BREAKDOWN"
            t_end = t_new
            break
        st = eng.band_status(g)
        if st is not None:
            status = st
            t_end = t_new
            break
        if n % n_out == 0:
            record(n)
        if keep_profiles and n % n_prof == 0:
            prof_t.append(t_new)
            prof.append(g.copy())
        if ckpt_path is not None and n % n_ckpt == 0:
            np.savez_compressed(ckpt_path, g=g, pi=pi, t=t_new, h=h,
                                dt=dt, lam=lam)

    eng.set_lam(lam_t((n0 + nsteps)*dt if t_end is None else t_end))
    return dict(status=status, t_end=t_end, g=g, pi=pi, eng=eng,
                t=np.array(ts), psi0=np.array(psi0s),
                mn=np.array(mns), mx=np.array(mxs),
                E=np.array(Es), Ef=np.array(Efs),
                V=np.array(Vs), D=np.array(Ds),
                prof_t=np.array(prof_t),
                prof=(np.array(prof) if prof else np.zeros((0, 0))),
                r=eng.r, h=h, dt=dt, lam=lam)


def window_mask(arr, a, b):
    return (arr >= a - 1e-9) & (arr <= b + 1e-9)


def classify_on(res, wlo=500.0, whi=600.0):
    """Klasyfikacja fazy wlaczonej (MD sec.6): COLLAPSE / SETTLED-SUB /
    SETTLED-DEF / UNSETTLED."""
    out = {}
    if res["status"] != "OK":
        out["cls"] = "COLLAPSE"
        out["subtype"] = res["status"]
        out["t_end"] = res["t_end"]
        return out
    mp = window_mask(res["prof_t"], wlo, whi)
    psibar = res["prof"][mp].mean(axis=0)
    r = res["r"]
    D = float(np.max(np.abs(psibar[r <= RMAX_E] - 1.0)))
    ms = window_mask(res["t"], wlo, whi)
    V = float(np.max(res["V"][ms]))
    settled = V <= 0.01*max(D, 1e-12)
    psibar0 = float(psibar[0])
    if not settled:
        cls = "UNSETTLED"
    elif psibar0 < 5.0/6.0:
        cls = "SETTLED-SUB"
    else:
        cls = "SETTLED-DEF"
    out.update(cls=cls, subtype="", t_end=None, psibar0=psibar0,
               dpsi0=psibar0 - 1.0, V=V, D=D, settled=bool(settled),
               n_prof=int(np.sum(mp)), psibar=psibar)
    return out


def classify_off(res, T0=2.0*np.pi, t_off_end=700.0, t_max=1700.0,
                 wf_lo=1600.0, wf_hi=1700.0):
    """Klasyfikacja po wygaszeniu (MD sec.7): COLLAPSE /
    PERSISTENT-OBJECT / RETURN-TO-VACUUM / INCONCLUSIVE-RUN."""
    out = {}
    if res["status"] != "OK":
        out["cls"] = "COLLAPSE"
        out["subtype"] = res["status"]
        out["t_end"] = res["t_end"]
        return out
    t = res["t"]
    i700 = int(np.argmin(np.abs(t - t_off_end)))
    Eref = float(res["E"][i700])
    wp_hi = t_off_end + 100.0*T0
    mp = window_mask(t, t_off_end, wp_hi)
    Ewp = res["E"][mp]
    Dwp = res["D"][mp]
    mf = window_mask(t, wf_lo, wf_hi)
    E_end = float(res["E"][-1])
    D_end = float(res["D"][-1])
    cond_E = bool(np.all(Ewp >= 0.5*Eref))
    cond_D = bool(np.all(Dwp >= 0.02))
    # tau_obj: najdluzszy nieprzerwany prefiks spelniajacy OBA warunki
    mall = window_mask(t, t_off_end, t_max)
    tt = t[mall]
    ok = (res["E"][mall] >= 0.5*Eref) & (res["D"][mall] >= 0.02)
    tau = 0.0
    for i in range(len(ok)):
        if ok[i]:
            tau = float(tt[i] - t_off_end)
        else:
            break
    if Eref <= 0.0:
        cls = "INCONCLUSIVE-RUN"
    elif cond_E and cond_D:
        cls = "PERSISTENT-OBJECT"
    elif (E_end < 0.05*Eref) or bool(np.all(res["D"][mf] < 1e-3)):
        cls = "RETURN-TO-VACUUM"
    else:
        cls = "INCONCLUSIVE-RUN"
    out.update(cls=cls, subtype="", t_end=None, Eref700=Eref,
               E_end=E_end, E_ratio=(E_end/Eref if Eref != 0 else
                                     float("nan")),
               D_end=D_end, D_wf_max=float(np.max(res["D"][mf])),
               cond_E=cond_E, cond_D=cond_D, tau=tau,
               tau_censored=bool(tau >= t_max - t_off_end - 1e-6))
    return out


def save_series(path, res, extra=None):
    d = dict(t=res["t"], psi0=res["psi0"], mn=res["mn"], mx=res["mx"],
             E=res["E"], Ef=res["Ef"], V=res["V"], D=res["D"],
             prof_t=res["prof_t"], r=res["r"])
    if extra:
        d.update(extra)
    np.savez_compressed(path, **d)


def save_state(path, res, t):
    np.savez_compressed(path, g=res["g"], pi=res["pi"], t=t,
                        h=res["h"], dt=res["dt"], lam=res["lam"])


def jsonable(d):
    out = {}
    for k, v in d.items():
        if isinstance(v, (np.floating, np.integer)):
            out[k] = float(v)
        elif isinstance(v, np.ndarray):
            continue
        elif isinstance(v, (bool, int, float, str, type(None))):
            out[k] = v
    return out


def ensure_dir(p):
    if not os.path.isdir(p):
        os.makedirs(p)
    return p
