#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
Addendum diagnostyczne do P1c (dokumentacja przyczyny FAIL; ZERO zmian
kryteriow). Pytania:
 (1) jaki jest pointwise przebieg Q_63 na profilu (gdzie jest kieszen
     ujemna, jaka glebokosc)?
 (2) co daje KOD RADIALNY (#63 weight-1) przy tym samym h co siatka 3D
     (h ~ 0.4, 0.3, ...)? Jesli radialny przy rownym h daje to samo co 3D,
     to most kartezjanski JEST spojny, a FAIL wynika z rozdzielczosci
     wymaganej przez model f_eps (sciana/rdzen), nie z dyskretyzacji 3D.
 (3) profil: g_min, g_max, F_eps na skrajach.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

PHI = (1 + 5 ** 0.5) / 2
G0 = PHI * 1.24915
AL = 2.0
EPS = 0.2
f_log = lambda g: 1 + 2 * AL * np.log(g)
fp_log = lambda g: 2 * AL / g
fpp_log = lambda g: -2 * AL / g ** 2
W1 = lambda g: g ** 2 * (1 - g)
W2 = lambda g: 2 * g - 3 * g ** 2


def mk(eps):
    s = lambda g: np.sqrt(f_log(g) ** 2 + eps ** 2)
    F = lambda g: 0.5 * (f_log(g) + s(g))
    Fp = lambda g: 0.5 * fp_log(g) * (1 + f_log(g) / s(g))
    Fpp = lambda g: 0.5 * (fpp_log(g) * (1 + f_log(g) / s(g))
                           + fp_log(g) ** 2 * eps ** 2 / s(g) ** 3)
    return F, Fp, Fpp


F, Fp, Fpp = mk(EPS)


def prof(R, N):
    def rhs(r_, y):
        g, gp = y
        g = max(g, 1e-3)
        drv = W1(g) - 0.5 * Fp(g) * gp ** 2
        if r_ < 1e-10:
            return [gp, drv / (3 * F(g))]
        return [gp, drv / F(g) - 2 * gp / r_]
    h = R / N
    r = (np.arange(N) + 0.5) * h
    sol = solve_ivp(rhs, [1e-6, R + h], [G0, 0.0], method='DOP853',
                    max_step=0.02, rtol=1e-10, atol=1e-13,
                    dense_output=True)
    r = r[r <= sol.t[-1]]
    v = sol.sol(r)
    g, gp = v[0], v[1]
    d2 = np.array([rhs(r[j], [g[j], gp[j]])[1] for j in range(len(r))])
    return r, g, gp, d2


def lam_min_radial(r, g, gp, d2):
    N = len(r)
    h = r[1] - r[0]
    Fm = F(0.5 * (g[:-1] + g[1:]))
    Q = W2(g) - 0.5 * Fpp(g) * gp ** 2 - Fp(g) * (d2 + 2 * gp / r)
    rm = 0.5 * (r[:-1] + r[1:])
    a = rm ** 2 * Fm
    d = np.empty(N)
    d[0] = a[0] / h ** 2
    d[1:-1] = (a[:-1] + a[1:]) / h ** 2
    d[-1] = (a[-1] + r[-1] ** 2 * F(g[-1])) / h ** 2
    d += r ** 2 * Q
    off = -a / h ** 2
    vals = eigh_tridiagonal(d / r ** 2, off / (r[:-1] * r[1:]),
                            select='i', select_range=(0, 2),
                            eigvals_only=True)
    return [float(x) for x in vals]


r4, g4, gp4, d24 = prof(60.0, 4000)
Q4 = W2(g4) - 0.5 * Fpp(g4) * gp4 ** 2 - Fp(g4) * (d24 + 2 * gp4 / r4)
i_min = int(np.argmin(Q4))
print("profil mu (fine): g in [%.4f, %.4f]; F_eps(g_min)=%.3e, "
      "F_eps(g_max)=%.3f" % (g4.min(), g4.max(), F(g4.min()), F(g4.max())))
print("Q_63 pointwise: min = %.4f przy r=%.3f (g=%.4f); Q(r_max)=%.4f"
      % (Q4[i_min], r4[i_min], g4[i_min], Q4[-1]))
print("Q_63 pierwszych 12 wezlow (r, g, Q):")
for j in range(0, 60, 5):
    print("  r=%.3f g=%.4f Q=%+.3f F=%.4f" % (r4[j], g4[j], Q4[j], F(g4[j])))

csg = CubicSpline(r4, g4)
csp = CubicSpline(r4, gp4)
csd = CubicSpline(r4, d24)
print("\nkod RADIALNY (#63 weight-1) przy roznych h (R=30, profil fine):")
for h in (0.3947, 0.30, 0.20, 0.10, 0.05, 0.025, 0.015):
    N = int(round(30.0 / h))
    rr = (np.arange(N) + 0.5) * h
    lam = lam_min_radial(rr, csg(rr), csp(rr), csd(rr))
    print("  h=%.4f N=%5d: lam_0..2 = %s"
          % (h, N, ["%+.4f" % x for x in lam]))
print("\n(3D N=76 dalo -8.8098, N=100 dalo -7.5371; radialny anchor "
      "h=0.015 R=60: -1.3896)")
