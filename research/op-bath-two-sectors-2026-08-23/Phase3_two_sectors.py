#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-bath-two-sectors (Phase 3, Q2) -- derywacyjny test dwoch sektorow:
czy znak tachionowy efektywnego potencjalu fluktuacji EMERGUJE z gestosci
zrodel w akcji STABILNEJ (M-S, eq:field-eq-reproduced, sek08a)?

LOCK: Phase0_balance.md sec. 2 (M-S), sec. 3 Phase 3; decyzje FROZEN:
Phase_method_decisions.md sec. Phase 3 (wariant wymiarowy: KOMORKA 1D
PERIODYCZNA, n ~ 1/d; sigma_s=0.5; q=1.0 PRIMARY + kontrola q=0.3 przy
d=4; operator fluktuacji z waga w=psi^4).

M-S:  psi'' + 2 psi'^2/psi + beta psi^2 - gamma psi^3 = -q rho~,
      beta = gamma = 1. Zrodla: siec gaussianow (okres d), norm. 1/okres.
P3a:  tla psi_n dla d in {inf, 8, 6, 4, 2} (Newton + kontynuacja w q).
      KONTROLA d=inf (osiagalny FAIL): zanik monotoniczny + m^2 =
      gamma +- 5% (fit ln(psi-1) na [8,20]); FAIL => STOP Q2.
P3b:  L phi = phi'' + (4 psi'/psi) phi' + [2 psi - 3 psi^2
      - 2 psi'^2/psi^2] phi; forma samosprzezona z w = psi^4:
      -(w phi')' - w Q phi = omega^2 w phi.  omega^2_min(d) na
      N in {400, 800}/okres. Odpowiedz statyczna: superkomorka
      10 okresow, L chi = -delta_h; klasyfikacja: zmiany znaku chi
      dla |x-x0| > 2 sigma_s (0 -> Yukawa; >=2 -> oscylacyjna).

Kryteria Q2 (LOCK, doslownie): PASS gdy istnieje d z omega^2_min < 0
LUB odpowiedzia oscylacyjna, zbieznie, przy czystej kontroli d=inf;
FAIL gdy omega^2_min > 0 i odpowiedz monotoniczna dla wszystkich d.

INPUT: sigma_s=0.5, q=1.0, wariant 1D. ZAKAZ przenoszenia czegokolwiek
z galezi tachionowej (Phase 1/2 nie zasilaja Q2 niczym).

Implementacja numeryczna: Jacobian/operator = macierze trojdiagonalne
(+ narozniki w periodycznych); solve_banded / eigh_tridiagonal / sparse.
"""
import numpy as np
from scipy.linalg import eigh, solve_banded, eigh_tridiagonal
import scipy.sparse as sp
import scipy.sparse.linalg as spla

BETA = 1.0
GAMMA = 1.0
SIGMA_S = 0.5     # INPUT (FROZEN)
Q_PRIM = 1.0      # INPUT (FROZEN)
D_LIST = (8.0, 6.0, 4.0, 2.0)
N_LIST = (400, 800)
M_SUPER = 10


def rho_lattice(x, d, sigma):
    r = np.zeros_like(x)
    for m in range(-6, 7):
        r += np.exp(-0.5 * ((x - d / 2 - m * d) / sigma) ** 2)
    return r / (np.sqrt(2 * np.pi) * sigma)


def neighbors(psi, periodic):
    if periodic:
        return np.roll(psi, -1), np.roll(psi, 1)
    ip = np.empty(len(psi))
    im = np.empty(len(psi))
    ip[:-1], ip[-1] = psi[1:], 1.0
    im[1:], im[0] = psi[:-1], 1.0
    return ip, im


def residual(psi, h, qrho, periodic):
    ip, im = neighbors(psi, periodic)
    d2 = (ip - 2 * psi + im) / h ** 2
    d1 = (ip - im) / (2 * h)
    return d2 + 2 * d1 ** 2 / psi + BETA * psi ** 2 - GAMMA * psi ** 3 \
        + qrho


def jac_bands(psi, h, periodic):
    """Trojdiagonalne pasma Jacobianu: (diag, upper, lower) + narozniki."""
    ip, im = neighbors(psi, periodic)
    d1 = (ip - im) / (2 * h)
    dg = -2 / h ** 2 - 2 * d1 ** 2 / psi ** 2 \
        + 2 * BETA * psi - 3 * GAMMA * psi ** 2
    up = 1 / h ** 2 + 4 * d1 / psi / (2 * h)    # sprzezenie do i+1
    lo = 1 / h ** 2 - 4 * d1 / psi / (2 * h)    # sprzezenie do i-1
    return dg, up, lo


def newton(N, h, qrho, periodic, tag=""):
    psi = np.ones(N)
    idx = np.arange(N)
    for frac in (0.25, 0.5, 0.75, 1.0):
        qr = frac * qrho
        ok = False
        for it in range(80):
            R = residual(psi, h, qr, periodic)
            nrm = float(np.max(np.abs(R)))
            # korekta 2 (Phase2_correction_note.md): prog >= podloga
            # kasowania w double dla czlonu (psi+ - 2 psi + psi-)/h^2
            tol = max(1e-10, 20 * np.finfo(float).eps
                      * float(np.max(psi)) / h ** 2)
            if nrm < tol:
                ok = True
                break
            dg, up, lo = jac_bands(psi, h, periodic)
            if periodic:
                J = sp.csr_matrix(
                    (np.concatenate([dg, up, lo]),
                     (np.concatenate([idx, idx, idx]),
                      np.concatenate([idx, (idx + 1) % N,
                                      (idx - 1) % N]))), shape=(N, N))
                dpsi = spla.spsolve(J, -R)
            else:
                ab = np.zeros((3, N))
                ab[0, 1:] = up[:-1]
                ab[1, :] = dg
                ab[2, :-1] = lo[1:]
                dpsi = solve_banded((1, 1), ab, -R)
            lam = 1.0
            for _ in range(40):
                trial = psi + lam * dpsi
                if float(np.min(trial)) > 1e-6 and float(np.max(np.abs(
                        residual(trial, h, qr, periodic)))) < nrm:
                    break
                lam *= 0.5
            psi = psi + lam * dpsi
        if not ok:
            print("  %s Newton NIE ZBIEGL (frac q=%.2f, ||R||=%.1e)"
                  % (tag, frac, nrm))
            return None
    return psi


def fluct_bands(psi, h, periodic):
    """A phi = omega^2 B phi; A tridiag (+narozniki), B=diag(w), w=psi^4.
    Zwraca (diagA, offA (i,i+1), w, cornerA)."""
    ip, im = neighbors(psi, periodic)
    d1 = (ip - im) / (2 * h)
    w = psi ** 4
    if periodic:
        wp = 0.5 * (w + np.roll(w, -1))
        wm = np.roll(wp, 1)
    else:
        wp = np.empty(len(w))
        wp[:-1] = 0.5 * (w[:-1] + w[1:])
        wp[-1] = 0.5 * (w[-1] + 1.0)
        wm = np.empty(len(w))
        wm[1:] = wp[:-1]
        wm[0] = 0.5 * (1.0 + w[0])
    Q = 2 * BETA * psi - 3 * GAMMA * psi ** 2 - 2 * d1 ** 2 / psi ** 2
    diagA = (wp + wm) / h ** 2 - w * Q
    offA = -wp / h ** 2          # (i, i+1); w periodycznym offA[-1]=naroznik
    return diagA, offA, w


def lowest_eigs(psi, h, periodic, k=4):
    diagA, offA, w = fluct_bands(psi, h, periodic)
    sw = np.sqrt(w)
    if not periodic:
        d = diagA / w
        e = offA[:-1] / (sw[:-1] * sw[1:])
        vals = eigh_tridiagonal(d, e, select='i',
                                select_range=(0, k - 1),
                                eigvals_only=True)
        return np.array(vals)
    N = len(psi)
    A = np.zeros((N, N))
    idx = np.arange(N)
    A[idx, idx] = diagA / w
    A[idx, (idx + 1) % N] = offA / (sw * sw[(idx + 1) % N])
    A[(idx + 1) % N, idx] = offA / (sw * sw[(idx + 1) % N])
    vals = eigh(A, eigvals_only=True, subset_by_index=(0, k - 1))
    return np.array(vals)


def sign_changes(y, thresh):
    s = np.sign(y[np.abs(y) > thresh])
    return int(np.sum(s[1:] != s[:-1]))


print("=" * 78)
print("Phase 3 (Q2): dwa sektory z jednej akcji stabilnej (M-S).")
print("Wariant FROZEN: komorka 1D periodyczna; sigma_s=0.5, q=1.0 [INPUT].")
print("Operator fluktuacji: waga w=psi^4; kotwica prozni: omega^2=gamma+k^2.")
print("=" * 78)

# ------------------------------------------------------- P3a: d=inf kontrola
print("\n[P3a] KONTROLA d=inf: pojedyncze zrodlo, [-40,40], Dirichlet psi=1")
stop_q2 = False
psi_inf = {}
for N_inf in (4000, 8000):
    L = 40.0
    h = 2 * L / N_inf
    x = np.arange(N_inf) * h - L + h / 2
    rho = np.exp(-0.5 * (x / SIGMA_S) ** 2) / (np.sqrt(2 * np.pi) * SIGMA_S)
    psi = newton(N_inf, h, Q_PRIM * rho, periodic=False,
                 tag="d=inf N=%d" % N_inf)
    if psi is None:
        stop_q2 = True
        continue
    psi_inf[N_inf] = (x, psi, h)
    mfit = (x >= 8) & (x <= 20) & (np.abs(psi - 1) > 1e-300)
    p = np.polyfit(x[mfit], np.log(np.abs(psi[mfit] - 1)), 1)
    m_meas = -p[0]
    m2 = m_meas ** 2
    dev = abs(m2 - GAMMA) / GAMMA * 100
    mm = (x >= 2) & (x <= 30)
    mono = bool(np.all(np.diff(np.abs(psi[mm] - 1)) <= 1e-14))
    ev = lowest_eigs(psi, h, periodic=False, k=3)
    print("  N=%d: max(psi)=%.6f; fit m=%.5f -> m^2=%.5f (gamma=1; "
          "odchyl %.2f%%, gate<=5%%); monotoniczny zanik [2,30]: %s; "
          "omega^2_min=%.5f (oczek. ~+gamma)"
          % (N_inf, float(np.max(psi)), m_meas, m2, dev, mono,
             float(ev[0])))
    if dev > 5.0 or not mono or ev[0] <= 0:
        stop_q2 = True

if stop_q2:
    print("\nKONTROLA d=inf FAIL -> STOP Q2 (procedura niewazna; "
          "osiagalny FAIL wystapil).")
    raise SystemExit(0)
print("  KONTROLA d=inf: PASS (Yukawa m^2=+gamma odtworzona).")

# ------------------------------------------------- P3a/P3b: siec periodyczna
print("\n[P3a/P3b] siec periodyczna: d in {8, 6, 4, 2}, q=1.0")
res_d = {}
for d in D_LIST:
    vals_grid = {}
    psi_keep = None
    for N in N_LIST:
        h = d / N
        x = np.arange(N) * h
        rho = rho_lattice(x, d, SIGMA_S)
        psi = newton(N, h, Q_PRIM * rho, periodic=True,
                     tag="d=%.0f N=%d" % (d, N))
        if psi is None:
            continue
        ev = lowest_eigs(psi, h, periodic=True, k=4)
        vals_grid[N] = float(ev[0])
        psi_keep = psi
        print("  d=%.0f N=%3d: psi in [%.6f, %.6f]; omega^2 najnizsze: %s"
              % (d, N, float(np.min(psi)), float(np.max(psi)),
                 ["%+.5f" % float(v) for v in ev]))
    if len(vals_grid) == 2:
        v1, v2 = vals_grid[N_LIST[0]], vals_grid[N_LIST[1]]
        mean = 0.5 * (v1 + v2)
        conv = abs(v1 - v2) <= 0.01 * max(abs(mean), 0.1)
        print("    zbieznosc siatek: |dv|=%.2e -> %s; "
              "omega^2_min(d=%.0f) = %+.5f"
              % (abs(v1 - v2), "ZBIEZNE" if conv else "NIEZBIEZNE",
                 d, mean))
        res_d[d] = dict(w2=mean, conv=conv, ok=True)
    else:
        res_d[d] = dict(w2=None, conv=False, ok=False)

# kontrola wrazliwosci q=0.3 przy d=4 (FROZEN w method_decisions)
h = 4.0 / N_LIST[1]
x = np.arange(N_LIST[1]) * h
rho = rho_lattice(x, 4.0, SIGMA_S)
psi_q = newton(N_LIST[1], h, 0.3 * rho, periodic=True, tag="d=4 q=0.3")
if psi_q is not None:
    ev = lowest_eigs(psi_q, h, periodic=True, k=2)
    print("  [wrazliwosc] d=4, q=0.3: omega^2_min=%+.5f" % float(ev[0]))

# ------------------------------------------- odpowiedz statyczna (superkom.)
print("\n[P3b] odpowiedz statyczna: superkomorka %d okresow, "
      "L chi = -delta_h w centrum zrodla" % M_SUPER)
resp = {}
N0 = N_LIST[0]
for d in D_LIST:
    if not res_d[d]["ok"]:
        continue
    h = d / N0
    x1 = np.arange(N0) * h
    rho1 = rho_lattice(x1, d, SIGMA_S)
    psi1 = newton(N0, h, Q_PRIM * rho1, periodic=True)
    if psi1 is None:
        continue
    psi_sc = np.tile(psi1, M_SUPER)
    Nsc = len(psi_sc)
    diagA, offA, w = fluct_bands(psi_sc, h, periodic=True)
    idx = np.arange(Nsc)
    A = sp.csr_matrix(
        (np.concatenate([diagA, offA, offA]),
         (np.concatenate([idx, idx, (idx + 1) % Nsc]),
          np.concatenate([idx, (idx + 1) % Nsc, idx]))),
        shape=(Nsc, Nsc))
    i0 = int((M_SUPER // 2) * N0 + round((d / 2) / h))
    rhs = np.zeros(Nsc)
    rhs[i0] = w[i0] * (1.0 / h)     # A chi = w delta_h  (L chi = -delta)
    chi = spla.spsolve(A, rhs)
    xs = np.arange(Nsc) * h
    x0 = xs[i0]
    mright = (xs > x0 + 2 * SIGMA_S) & (xs < x0 + (M_SUPER // 2 - 0.5) * d)
    y = chi[mright]
    nsc = sign_changes(y, 1e-10 * float(np.max(np.abs(chi))))
    cls = "OSCYLACYJNA" if nsc >= 2 else "MONOTONICZNA/YUKAWA"
    dec = float(np.abs(y[-1])) < float(np.abs(y[0]))
    print("  d=%.0f: zmiany znaku (|x-x0|>2 sigma_s): %d -> %s; "
          "|chi| zanika: %s" % (d, nsc, cls, dec))
    resp[d] = nsc

# kontrola odpowiedzi: proznia (czysta Yukawa referencyjna)
Nv = 4000
L = 40.0
h = 2 * L / Nv
diagA, offA, w = fluct_bands(np.ones(Nv), h, periodic=False)
ab = np.zeros((3, Nv))
ab[0, 1:] = offA[:-1]
ab[1, :] = diagA
ab[2, :-1] = offA[:-1]
rhs = np.zeros(Nv)
i0 = Nv // 2
rhs[i0] = w[i0] / h
chi = solve_banded((1, 1), ab, rhs)
x = np.arange(Nv) * h - L + h / 2
mright = (x > x[i0] + 2 * SIGMA_S) & (x < 30)
nsc_v = sign_changes(chi[mright], 1e-10 * float(np.max(np.abs(chi))))
print("  kontrola prozni (d=inf): zmiany znaku = %d (oczek. 0, Yukawa)"
      % nsc_v)

# ------------------------------------------------------------------ WERDYKT
print("\n" + "=" * 78)
print("PODSUMOWANIE Q2 (kryteria LOCK sec. 3 Phase 3, doslownie):")
any_tach = False
all_pos_mono = True
for d in D_LIST:
    r = res_d[d]
    if r["w2"] is None:
        all_pos_mono = False
        print("  d=%.0f: tlo niezbiegniete -> INCONCLUSIVE-pkt" % d)
        continue
    osc = resp.get(d, 0) >= 2
    tach = (r["w2"] < 0 and r["conv"]) or osc
    if tach:
        any_tach = True
    if not (r["w2"] > 0 and r["conv"] and not osc):
        all_pos_mono = False
    print("  d=%.0f: omega^2_min=%+.5f (conv=%s), odpowiedz: %s"
          % (d, r["w2"], r["conv"],
             "oscylacyjna" if osc else "monotoniczna"))
if any_tach:
    verdict = "Q2-PASS (znak tachionowy emerguje z gestosci)"
elif all_pos_mono:
    verdict = ("Q2-FAIL (gestosc NIE zmienia znaku; rozdwojenie W "
               "pozostaje problemem aksjomatycznym)")
else:
    verdict = "Q2-INCONCLUSIVE"
print("\nWERDYKT: %s" % verdict)
print("(kontrola d=inf: PASS; interpretacja -> Phase_FINAL_close.md)")
print("=" * 78)
