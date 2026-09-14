#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-premetric-pocket-level0 (Phase 1) -- bramka ciaglosci z zamknietym
cyklem op-bare-substrate-genesis (LOCK sec. 2 Phase 1). Silnik = import
z Phase2_pocket (ta sama maszyneria).
P1a: bare (szum U(-0.05,0.05), 6000 krokow) -> metric_area < A_min (G2).
P1b: single(A0=1.4, w=1.5) po 6000 krokach zanika (G3, najwiekszy
     podkrytyczny z tamtego skanu).
P1c: H_Gamma nierosnace we wszystkich biegach P1.
P1d: 1D profil sciany Z2: relaksacja do 1e-10; zbieznosc dx->dx/2
     (szerokosci pasow <1% wzgl.; Phi_min roznica <1% Phi_bar --
     odczyt: Phi_min ~ 0 z przeciecia, prog wzgledny zle okreslony).
Dowolny FAIL ==> STOP.
"""
import sys
import numpy as np

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-premetric-pocket-level0-2026-09-02/")
sys.path.insert(0, BASE)
from Phase2_pocket import (Sub2D, Vpot, Vprime, noise, banner,  # noqa: E402
                           r_center, KAPPA, EPS, AMIN_FRAC, S_STAR,
                           S_BAR)

PHI_BAR = S_BAR ** 2
lines = []


def em(x=""):
    print(x, flush=True)
    lines.append(x)


em("=" * 78)
em("Phase 1 -- bramka ciaglosci (silnik = import Phase2_pocket)")
banner("; Phase1 gate")
em("=" * 78)
allok = True

# ------------------------------------------------------------ P1a, P1b
for name, s0 in (
        ("P1a bare (szum)", noise(128).copy()),
        ("P1b single(A0=1.4,w=1.5)",
         1.4 * np.exp(-r_center(128, 0.5) ** 2 / (2 * 1.5 ** 2)))):
    eng = Sub2D(128, 0.5, 0.02)
    s = s0
    Hprev, mono = None, True
    for k in range(6000):
        s = eng.step(s)
        if k % 50 == 0:
            Hn = eng.H(s)
            if Hprev is not None and Hn > Hprev + 1e-10:
                mono = False
            Hprev = Hn
    marea = float(np.mean(s * s >= EPS))
    ok = marea < AMIN_FRAC and mono
    allok &= ok
    em("  %-28s metric_area=%.3e (<A_min=%.3e) Hmono=%s -> %s"
       % (name, marea, AMIN_FRAC, mono, "OK" if ok else "FAIL"))
em("  P1a/P1b/P1c: %s" % ("PASS" if allok else "FAIL"))
em()

# ------------------------------------------------------------ P1d: 1D
em("P1d: 1D profil sciany Z2 (Dirichlet -s*..+s*):")


def wall1d(Npts, dx, dt, itmax):
    x = (np.arange(Npts) + 0.5) * dx
    Ltot = Npts * dx
    s = S_STAR * np.tanh((x - Ltot / 2) / 1.07)
    s[0], s[-1] = -S_STAR, S_STAR
    res = None
    for k in range(itmax):
        lap = np.zeros_like(s)
        lap[1:-1] = (s[2:] - 2 * s[1:-1] + s[:-2]) / dx ** 2
        f = KAPPA * lap - Vprime(s)
        f[0] = f[-1] = 0.0
        s = s + dt * f
        if k % 2000 == 0:
            res = float(np.max(np.abs(f[1:-1])))
            if res <= 1e-10:
                break
    Phi = s * s
    dsx = np.diff(s) / dx
    H = float(np.sum(Vpot(s)) * dx + 0.5 * KAPPA * np.sum(dsx ** 2) * dx)
    sigma = H - Vpot(S_STAR) * Ltot

    def band(th):
        """Szerokosc pasa Phi<th z interpolacja liniowa OBU brzegow
        (korekta 1: Phase_correction_note_p1d_band.md)."""
        m = Phi < th
        if not m.any():
            return 0.0
        idx = np.where(m)[0]
        i0, i1 = int(idx[0]), int(idx[-1])
        xl = x[i0]
        if i0 > 0 and Phi[i0 - 1] != Phi[i0]:
            xl = x[i0 - 1] + dx * (Phi[i0 - 1] - th) / \
                (Phi[i0 - 1] - Phi[i0])
        xr = x[i1]
        if i1 < len(x) - 1 and Phi[i1 + 1] != Phi[i1]:
            xr = x[i1] + dx * (Phi[i1] - th) / (Phi[i1] - Phi[i1 + 1])
        return float(xr - xl)
    return dict(res=res, phimin=float(Phi.min()),
                w_eps=band(EPS), w_bar=band(PHI_BAR), sigma=sigma)


w1 = wall1d(4096, 0.05, 0.001, 400000)
w2 = wall1d(8192, 0.025, 0.00025, 1600000)
em("  dx=0.05 : res=%.2e Phi_min=%.3e szer(Phi<0.30)=%.4f "
   "szer(Phi<Phi_bar)=%.4f sigma_w=%.6f"
   % (w1["res"], w1["phimin"], w1["w_eps"], w1["w_bar"], w1["sigma"]))
em("  dx=0.025: res=%.2e Phi_min=%.3e szer(Phi<0.30)=%.4f "
   "szer(Phi<Phi_bar)=%.4f sigma_w=%.6f"
   % (w2["res"], w2["phimin"], w2["w_eps"], w2["w_bar"], w2["sigma"]))
rel_eps = abs(w1["w_eps"] - w2["w_eps"]) / w2["w_eps"]
rel_bar = abs(w1["w_bar"] - w2["w_bar"]) / w2["w_bar"]
dphi = abs(w1["phimin"] - w2["phimin"]) / PHI_BAR
ok1d = (w1["res"] <= 1e-10 and w2["res"] <= 1e-10 and rel_eps < 0.01
        and rel_bar < 0.01 and dphi < 0.01)
allok &= ok1d
em("  zbieznosc: d szer(eps)=%.3e%% d szer(bar)=%.3e%% "
   "dPhi_min/Phi_bar=%.3e (<1%%) -> %s"
   % (100 * rel_eps, 100 * rel_bar, dphi, "OK" if ok1d else "FAIL"))
em("  deskryptywnie: rdzen sciany Phi_min ~ %.1e -- sciana Z2 MA rdzen"
   " przedmetryczny (psi_min=%.2e); pas psi<0.218 szeroki %.2f j."
   % (w2["phimin"], w2["phimin"] / S_STAR ** 2, w2["w_eps"]))
em()
em("PODSUMOWANIE Phase 1: %s" % ("PASS" if allok else "FAIL -- STOP"))
with open(BASE + "Phase1_output.txt", "w") as f:
    f.write("\n".join(lines) + "\n")
print("zapisano:", BASE + "Phase1_output.txt")
if not allok:
    sys.exit(1)
