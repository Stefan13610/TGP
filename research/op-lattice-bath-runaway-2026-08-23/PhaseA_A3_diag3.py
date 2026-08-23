#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
Phase A / A3 -- DIAGNOSTYKA nr 3: zaleznosc (A, delta) od OKNA FITU.
Rdzenne skrypty fituja ogon w roznych oknach: a3d [20,35],
atail_* [50,200]. Jesli ogon ma wolno gasnacy transjent, (A, delta)
zaleza od okna -- wtedy wartosci p127/p131 moga byc artefaktem okna.
Testowane: S1 (alpha=2) i wariant biegnacy (eta=181/15) dla e i mu,
oraz S2 (alpha=1) dla mu.
"""
import numpy as np

PHI = (1.0 + 5.0 ** 0.5) / 2.0
ETA = 181.0 / 15.0


def rk4_generic(g0, accel, r_max=320.0, h=0.001, stride=10, g_floor=1e-6,
                d=3):
    r0 = 1e-4
    g = float(g0)
    a0 = accel(1.0, g, 0.0, True)
    gp = a0 * r0 / d
    g = g + 0.5 * a0 * r0 * r0 / d
    rs, gs = [], []
    r = r0
    nstep = int(np.ceil((r_max - r0) / h))
    for k in range(nstep):
        if k % stride == 0:
            rs.append(r)
            gs.append(g)
        k1g, k1p = gp, accel(r, g, gp, False)
        g2, p2 = g + h / 2 * k1g, gp + h / 2 * k1p
        k2g, k2p = p2, accel(r + h / 2, g2, p2, False)
        g3, p3 = g + h / 2 * k2g, gp + h / 2 * k2p
        k3g, k3p = p3, accel(r + h / 2, g3, p3, False)
        g4, p4 = g + h * k3g, gp + h * k3p
        k4g, k4p = p4, accel(r + h, g4, p4, False)
        g += h / 6 * (k1g + 2 * k2g + 2 * k3g + k4g)
        gp += h / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
        r += h
        if g <= g_floor or not np.isfinite(g):
            return np.array(rs), np.array(gs), True
    return np.array(rs), np.array(gs), False


def fit_win(r, g, w_lo, w_hi):
    y = (g - 1.0) * r
    m = (r >= w_lo) & (r <= w_hi)
    M = np.column_stack([np.cos(r[m]), np.sin(r[m])])
    (B, C), _, _, _ = np.linalg.lstsq(M, y[m], rcond=None)
    return float(np.hypot(B, C)), float(np.degrees(np.arctan2(C, B)))


def acc_S1(r, g, gp, at0):
    src = g * g * (1 - g) - (2.0 / g) * gp * gp
    return src if at0 else src - 2 * gp / r


def acc_S2(r, g, gp, at0):
    src = g * g * (1 - g) - (1.0 / g) * gp * gp
    return src if at0 else src - 2 * gp / r


def acc_run(r, g, gp, at0):
    ae = 2.0 / (1.0 + ETA * (g - 1.0) ** 2)
    src = g * g * (1 - g) - (ae / g) * gp * gp
    return src if at0 else src - 2 * gp / r


WINDOWS = [(20, 35), (30, 60), (45, 58), (50, 100), (50, 200),
           (100, 200), (150, 300)]

CASES = [
    ("S1 a=2  e  g0=0.8993 (ref A=0.1018 d=-81.1)", acc_S1, 0.8993),
    ("S1 a=2  mu g0=1.4550 (ref A=0.3861 d=+43.8)", acc_S1, 1.4550),
    ("S2 a=1  mu g0=1.4550 (ref A=0.3861 d=+43.8)", acc_S2, 1.4550),
    ("RUN eta e  g0=0.90548 (tgt -81.4)", acc_run, 0.90548),
    ("RUN eta mu g0=1.46510 (tgt +38.6)", acc_run, PHI * 0.90548),
]
for tag, acc, g0 in CASES:
    r, g, dead = rk4_generic(g0, acc)
    print("[%s]%s" % (tag, "  KOLAPS" if dead else ""))
    if dead:
        continue
    for (lo, hi) in WINDOWS:
        A, dd = fit_win(r, g, lo, hi)
        print("   okno [%3d,%3d]: A=%.4f  delta=%+7.2f" % (lo, hi, A, dd))
    print()
