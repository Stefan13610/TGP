#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
Addendum diagnostyczne 2 do P1c (rozstrzygniecie hipotezy bledu konwencji
wagi; ZERO zmian kryteriow). Testy:
 (T1) radialne widmo F-WAZONE (#63: lam_min(F) = -7.8579 przy h=0.015,
      R=60) przy h in {0.3947, 0.30, 0.015} na R=30 -- czy pokrywa sie
      z 3D (-8.8098 / -7.5371)? Jesli TAK -> blad konwencji TESTU 3D.
 (T2) iloraz Rayleigha DOKLADNEGO modu w1 (u1/r z N=4000, ||v||_B=1,
      interpolowany na siatke 3D) w operatorze 3D weight-1:
      RQ = <v, A v> / <v, v>. Jesli RQ ~ -1.39 przy lam_min(eigsh)=-8.81,
      operator 3D reprezentuje kotwiczny mod POPRAWNIE, a nizsze wartosci
      to dodatkowe mody kieszeni Q (artefakt rozdzielczosci powloki).
 (T3) lokalizacja modu -8.81: rozklad promieniowy |v|^2 (srodek masy,
      frakcja w powloce |r-3.38|<0.3).
 (T4) min pointwise diag operatora 3D (Sum c/h^2 + Q) -- dolne
      oszacowanie skali modow pasozytniczych.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
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
s_ = lambda g: np.sqrt(f_log(g) ** 2 + EPS ** 2)
F = lambda g: 0.5 * (f_log(g) + s_(g))
Fp = lambda g: 0.5 * fp_log(g) * (1 + f_log(g) / s_(g))
Fpp = lambda g: 0.5 * (fpp_log(g) * (1 + f_log(g) / s_(g))
                       + fp_log(g) ** 2 * EPS ** 2 / s_(g) ** 3)


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


def tridiag(r, g, gp, d2):
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
    return d / r ** 2, off / (r[:-1] * r[1:])


r4, g4, gp4, d24 = prof(60.0, 4000)
csg = CubicSpline(r4, g4 - 1.0)
csp = CubicSpline(r4, gp4)
csd = CubicSpline(r4, d24)
d2_0 = float(csd(0.0))
R_END = float(r4[-1])

# referencja: waga-1 i F-wazone na siatce #63 (h=0.015, R=60)
d0, e0 = tridiag(r4, g4, gp4, d24)
lam1_ref = float(eigh_tridiagonal(d0, e0, select='i', select_range=(0, 0),
                                  eigvals_only=True)[0])
sF = np.sqrt(F(g4))
lamF_ref = float(eigh_tridiagonal(d0 / F(g4), e0 / (sF[:-1] * sF[1:]),
                                  select='i', select_range=(0, 0),
                                  eigvals_only=True)[0])
print("referencja radialna (h=0.015, R=60): lam_min(w1)=%.5f, "
      "lam_min(F)=%.5f (#63: -1.3896 / -7.8579)" % (lam1_ref, lamF_ref))

print("\n[T1] radialne widmo F-WAZONE i waga-1 przy grubych h (R=30):")
for h in (0.3947, 0.30, 0.10, 0.015):
    N = int(round(30.0 / h))
    rr = (np.arange(N) + 0.5) * h
    gg = 1.0 + csg(rr)
    pp = csp(rr)
    dd = csd(rr)
    d_, e_ = tridiag(rr, gg, pp, dd)
    l1 = float(eigh_tridiagonal(d_, e_, select='i', select_range=(0, 0),
                                eigvals_only=True)[0])
    sf = np.sqrt(F(gg))
    lF = float(eigh_tridiagonal(d_ / F(gg), e_ / (sf[:-1] * sf[1:]),
                                select='i', select_range=(0, 0),
                                eigvals_only=True)[0])
    print("  h=%.4f N=%5d: lam_min(w1)=%+.4f   lam_min(F)=%+.4f"
          % (h, N, l1, lF))
print("  (3D dalo: N=76/h=0.3947 -> -8.8098; N=100/h=0.30 -> -7.5371)")

# ---------------- operator 3D weight-1 (identyczny z Phase1_gate3d.py)
def minimg_r(N, h, L):
    x = np.arange(N) * h
    dx = np.mod(x - L / 2 + L / 2, L) - L / 2
    DX, DY, DZ = np.meshgrid(dx, dx, dx, indexing='ij')
    return np.sqrt(DX ** 2 + DY ** 2 + DZ ** 2)


def radial_fields(rr):
    inside = rr <= R_END
    g3 = np.where(inside, 1.0 + csg(np.minimum(rr, R_END)), 1.0)
    gp3 = np.where(inside, csp(np.minimum(rr, R_END)), 0.0)
    rr_safe = np.maximum(rr, 1e-12)
    lap3 = np.where(inside, csd(np.minimum(rr, R_END)) + 2 * gp3 / rr_safe,
                    0.0)
    lap3 = np.where(rr < 1e-9, 3 * d2_0, lap3)
    return g3, gp3, lap3


def build3d(g3, Q3, h):
    N = g3.shape[0]
    n = N ** 3
    idx3 = np.arange(n).reshape(N, N, N)
    diag = Q3.ravel().astype(float).copy()
    rows, cols, vals = [], [], []
    for ax in range(3):
        gn = np.roll(g3, -1, axis=ax)
        c = F(0.5 * (g3 + gn)) / h ** 2
        diag += (c + np.roll(c, 1, axis=ax)).ravel()
        r_ = idx3.ravel()
        c_ = np.roll(idx3, -1, axis=ax).ravel()
        v = (-c).ravel()
        rows += [r_, c_]
        cols += [c_, r_]
        vals += [v, v]
    rows.append(np.arange(n))
    cols.append(np.arange(n))
    vals.append(diag)
    return sp.coo_matrix((np.concatenate(vals),
                          (np.concatenate(rows), np.concatenate(cols))),
                         shape=(n, n)).tocsr(), diag

# mod w1 (kotwiczny) z N=4000: v = u/r
d0v, e0v = tridiag(r4, g4, gp4, d24)
vals0, vecs0 = eigh_tridiagonal(d0v, e0v, select='i', select_range=(0, 0))
u1 = vecs0[:, 0]
v1 = u1 / r4
csv1 = CubicSpline(r4, v1)
v1_0 = float(csv1(0.0))

L3 = 30.0
for N in (76, 100):
    h = L3 / N
    rr = minimg_r(N, h, L3)
    g3, gp3, lap3 = radial_fields(rr)
    Q3 = W2(g3) - 0.5 * Fpp(g3) * gp3 ** 2 - Fp(g3) * lap3
    A, diagA = build3d(g3, Q3, h)
    # [T4] min pointwise diag
    jmin = int(np.argmin(diagA))
    print("\n[T4] N=%d: min diag(A) = %+.3f przy r=%.3f (skala modow"
          " pasozytniczych)" % (N, float(diagA[jmin]),
                                float(rr.ravel()[jmin])))
    # [T2] Rayleigh dokladnego modu w1
    vd = np.where(rr <= R_END, csv1(np.minimum(rr, R_END)), 0.0)
    vd = np.where(rr < 1e-9, v1_0, vd)
    vv = vd.ravel()
    rq = float(vv @ (A @ vv)) / float(vv @ vv)
    print("[T2] N=%d: Rayleigh(mod w1 interpolowany) = %+.5f "
          "(kotwica w1: %.5f)" % (N, rq, lam1_ref))
    # [T3] lokalizacja modu najnizszego (tylko N=76 -- koszt)
    if N == 76:
        rng = np.random.default_rng(20260831)
        v0 = rng.standard_normal(N ** 3)
        v0 /= np.linalg.norm(v0)
        vals, V = spla.eigsh(A, k=3, which='SA', tol=1e-6, ncv=60,
                             maxiter=200000, v0=v0)
        o = np.argsort(vals)
        vmin = V[:, o[0]] ** 2
        vmin /= vmin.sum()
        rmean = float(np.sum(vmin * rr.ravel()))
        shell = float(np.sum(vmin[np.abs(rr.ravel() - 3.38) < 0.3]))
        core = float(np.sum(vmin[rr.ravel() < 2.0]))
        print("[T3] N=%d: lam=%s; mod najnizszy: <r>=%.3f, frakcja "
              "|r-3.38|<0.3: %.3f, frakcja r<2: %.3f"
              % (N, ["%+.4f" % x for x in np.sort(vals)], rmean, shell,
                 core))
