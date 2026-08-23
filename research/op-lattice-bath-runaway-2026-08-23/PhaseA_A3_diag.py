#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
Phase A / A3 -- DIAGNOSTYKA IMPLEMENTACJI (przed werdyktem A3).

Cel: pierwsze uruchomienie A3 (podstawienie bezposrednie alpha_eff
w czlonie krzyzowym) dalo delta_e=-75.5, delta_mu=+88.5, tau kolapsuje
-- niezgodne z dodatekH. LOCK dopuszcza korekte wylacznie dla
udokumentowanego bledu IMPLEMENTACJI. Ten skrypt:
 (a) WALIDUJE pipeline fitu na przypadku PLAIN alpha=2 przeciw
     niezaleznym liczbom rdzenia (p127-p128, mapa fazowa):
     g0_e=0.8993 -> A=0.1018, delta=-81.1 deg;
     g0_mu=1.4550 -> A=0.3861, delta=+43.8 deg; Delta=124.98 deg.
 (b) Testuje DWA warianty rownania z biegnacym alpha_eff:
     V1 (uzyte pierwotnie): g''+2g'/r+(ae(g)/g) g'^2 = g^2(1-g)
     V2 (EL z K(g)=g^(2 ae(g)), zrodlo g^2(1-g) po podzieleniu przez K):
        g''+2g'/r+[ae/g + ae'(g) ln g] g'^2 = g^2(1-g)
     dla eta_K=181/15, g0_e=0.90548, g0_mu=phi*g0_e, g0_tau=4.
Wynik: identyfikacja konwencji rdzenia LUB potwierdzenie
niereprodukowalnosci (wtedy A3 FAIL stoi).
"""
import numpy as np

PHI = (1.0 + 5.0 ** 0.5) / 2.0
ETA = 181.0 / 15.0
W_LO, W_HI = 40.0, 150.0


def rk4_g(g0, cross_fn, r_max=200.0, h=0.001, stride=20, g_floor=1e-6):
    r0 = 1e-4
    g = float(g0)
    src0 = g * g * (1 - g)
    gp = src0 * r0 / 3.0
    g = g + 0.5 * src0 * r0 * r0 / 3.0

    def rhs(r, gv, gpv):
        gc = max(gv, g_floor)
        return gpv, (gc * gc * (1 - gc) - cross_fn(gc) * gpv * gpv
                     - 2.0 * gpv / r)

    rs, gs = [], []
    r = r0
    dead = False
    nstep = int(np.ceil((r_max - r0) / h))
    for k in range(nstep):
        if k % stride == 0:
            rs.append(r)
            gs.append(g)
        k1g, k1p = rhs(r, g, gp)
        k2g, k2p = rhs(r + h / 2, g + h / 2 * k1g, gp + h / 2 * k1p)
        k3g, k3p = rhs(r + h / 2, g + h / 2 * k2g, gp + h / 2 * k2p)
        k4g, k4p = rhs(r + h, g + h * k3g, gp + h * k3p)
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
    M = np.column_stack([np.cos(r[m]), np.sin(r[m])])
    (B, C), _, _, _ = np.linalg.lstsq(M, y[m], rcond=None)
    A = float(np.hypot(B, C))
    return A, float(np.degrees(np.arctan2(C, B)))


print("(a) WALIDACJA pipeline: plain alpha=2 vs p127-p128 (mapa fazowa)")
cross_plain = lambda gc: 2.0 / gc
for g0, A_ref, d_ref in ((0.8993, 0.1018, -81.1), (1.4550, 0.3861, 43.8)):
    r, g, dead = rk4_g(g0, cross_plain)
    A, d = fit_phase(r, g)
    print("  g0=%.4f: A=%.4f (rdzen %.4f), delta=%+.2f (rdzen %+.1f)%s"
          % (g0, A, A_ref, d, d_ref, "  [KOLAPS]" if dead else ""))

print()
print("(b) warianty biegnacego alpha_eff (eta=181/15), g0: e=0.90548,")
print("    mu=phi*g0_e=%.5f, tau=4.0; targety: -81.4 / +38.6 / (tau -27.3)"
      % (PHI * 0.90548))


def ae(gc):
    return 2.0 / (1.0 + ETA * (gc - 1.0) ** 2)


def ae_prime(gc):
    return -4.0 * ETA * (gc - 1.0) / (1.0 + ETA * (gc - 1.0) ** 2) ** 2


variants = {
    "V1 cross=ae/g": lambda gc: ae(gc) / gc,
    "V2 cross=ae/g+ae'*ln g": lambda gc: ae(gc) / gc
        + ae_prime(gc) * np.log(gc),
}
for vname, cf in variants.items():
    print("  [%s]" % vname)
    for name, g0 in (("e", 0.90548), ("mu", PHI * 0.90548), ("tau", 4.0)):
        r, g, dead = rk4_g(g0, cf)
        if dead:
            print("    %-3s g0=%.5f: KOLAPS (g->0) -- brak ogona" % (name, g0))
            continue
        A, d = fit_phase(r, g)
        print("    %-3s g0=%.5f: A=%.4f  delta=%+.2f deg" % (name, g0, A, d))
