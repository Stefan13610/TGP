#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
Phase A / A3 -- DIAGNOSTYKA nr 2: identyfikacja ukladu, z ktorego
pochodzi mapa fazowa dodatekH p127-p128 (delta_e=-81.1 przy g0=0.8993,
A=0.1018; delta_mu=+43.8 przy g0=1.4550, A=0.3861; siodlo ~1.99).

WSKAZOWKA STRUKTURALNA: 'siodlo' mapy przy g0 ~ 1.99 = prog kolapsu
2(a+2)/[2(a+2)-d] dla alpha=1 (= 2.0), NIE dla alpha=2 (= 1.6).
Hipoteza: mapa fazowa liczona przy alpha=1 (K=g^2), nie alpha=2.

Kandydaci (wszystkie z tym samym fitem (g-1)r ~ B cos r + C sin r,
okno [40,150], delta=atan2(C,B)):
  S1: alpha=2, zrodlo g^2(1-g)   [odrzucony w diag 1]
  S2: alpha=1, zrodlo g^2(1-g)   [prog 2.0 -- glowna hipoteza]
  S3: alpha=1, zrodlo (1-g)      [Form B z ode_koide_formA_exact]
  S4: F-S log (f_kin=1+4 ln g), zrodlo g^2(1-g), krzyz (2/g)g'^2
      [rownanie z a3d_soliton_brannen_r bez odbic]
Oraz: dla zwycieskiego kandydata -- wariant z biegnacym
alpha_eff(g)=2/(1+eta(g-1)^2), eta=181/15, g0_e=0.90548 (targety
LOCK A3: -81.4 / +38.6 / Delta=120.01).
"""
import numpy as np

PHI = (1.0 + 5.0 ** 0.5) / 2.0
ETA = 181.0 / 15.0
W_LO, W_HI = 40.0, 150.0


def rk4_generic(g0, accel, r_max=200.0, h=0.001, stride=20,
                g_floor=1e-6, d=3):
    """g'' = accel(r, g, g'); start szeregowy g''(0)=accel0/d."""
    r0 = 1e-4
    g = float(g0)
    a0 = accel(1.0, g, 0.0, True)
    gp = a0 * r0 / d
    g = g + 0.5 * a0 * r0 * r0 / d
    rs, gs = [], []
    r = r0
    dead = False
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
            dead = True
            break
    return np.array(rs), np.array(gs), dead


def fit_phase(r, g, w_lo=W_LO, w_hi=W_HI):
    y = (g - 1.0) * r
    m = (r >= w_lo) & (r <= w_hi)
    if m.sum() < 50:
        return np.nan, np.nan
    M = np.column_stack([np.cos(r[m]), np.sin(r[m])])
    (B, C), _, _, _ = np.linalg.lstsq(M, y[m], rcond=None)
    return float(np.hypot(B, C)), float(np.degrees(np.arctan2(C, B)))


def acc_S2(r, g, gp, at0):
    src = g * g * (1 - g) - (1.0 / g) * gp * gp
    return src if at0 else src - 2 * gp / r


def acc_S3(r, g, gp, at0):
    src = (1 - g) - (1.0 / g) * gp * gp
    return src if at0 else src - 2 * gp / r


def acc_S4(r, g, gp, at0):
    fk = 1.0 + 4.0 * np.log(max(g, 1e-6))
    if abs(fk) < 1e-8:
        fk = 1e-8
    drv = g * g * (1 - g) - (2.0 / g) * gp * gp
    if at0:
        return drv / fk
    return drv / fk - 2 * gp / r


REF = ((0.8993, 0.1018, -81.1), (1.4550, 0.3861, +43.8))
for tag, acc in (("S2 alpha=1, src g^2(1-g)", acc_S2),
                 ("S3 alpha=1, src (1-g)   ", acc_S3),
                 ("S4 F-S log (bez odbic)  ", acc_S4)):
    print("[%s]" % tag)
    for g0, A_ref, d_ref in REF:
        r, g, dead = rk4_generic(g0, acc)
        if dead:
            print("   g0=%.4f: KOLAPS" % g0)
            continue
        A, dd = fit_phase(r, g)
        print("   g0=%.4f: A=%.4f (ref %.4f)  delta=%+.2f (ref %+.1f)"
              % (g0, A, A_ref, dd, d_ref))
    print()

print("Siodlo / A_max -- skan g0 w [1.90, 1.999] dla S2 (prog 2.0):")
for g0 in (1.90, 1.95, 1.99):
    r, g, dead = rk4_generic(g0, acc_S2)
    if dead:
        print("   g0=%.3f: KOLAPS" % g0)
    else:
        A, dd = fit_phase(r, g)
        print("   g0=%.3f: A=%.4f  delta=%+.2f" % (g0, A, dd))

print()
print("Zwyciezca + biegnace alpha_eff (jesli S2 trafia):")


def make_acc_run(eta):
    def acc(r, g, gp, at0):
        ae = 2.0 / (1.0 + eta * (g - 1.0) ** 2)
        src = g * g * (1 - g) - (ae / g) * gp * gp
        return src if at0 else src - 2 * gp / r
    return acc


accR = make_acc_run(ETA)
for name, g0 in (("e", 0.90548), ("mu", PHI * 0.90548), ("tau", 4.0),
                 ("tau-", 3.999)):
    r, g, dead = rk4_generic(g0, accR, h=0.0005)
    if dead:
        print("   %-4s g0=%.5f: KOLAPS" % (name, g0))
        continue
    A, dd = fit_phase(r, g)
    print("   %-4s g0=%.5f: A=%.4f  delta=%+.2f  (targety: e -81.4, "
          "mu +38.6, tau -27.3)" % (name, g0, A, dd))
