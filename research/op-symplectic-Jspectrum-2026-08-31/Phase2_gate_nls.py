#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-symplectic-Jspectrum (Phase 2) -- bramka maszynerii numerycznej
(osiagalne FAIL-e w OBIE strony). LOCK: Phase0_balance.md sec. 3 Phase 2;
decyzje FROZEN: Phase_method_decisions.md sec. 2, 4, 6.

Gate A (kotwica): reprodukcja operatora L_plus cyklu bloch -- problem
   wazony L_plus phi = lambda K phi na tle 3pi (npz), N=800, k=0:
   lambda_min = -1.222191 +/- 1e-4. FAIL => STOP.
C1 (znany STABILNY): soliton NLS kubicznego (i u_t + u_xx + |u|^2 u = 0,
   ramka omega=1 wchlonieta: K=2, U=g^2-g^4/2): max Re lambda <= tol
   mimo ujemnego kierunku L_plus. FAIL => STOP.
C2 (znany NIESTABILNY): NLS |u|^6 u (nadkrytyczny, K=2, U=g^2-g^8/4):
   MUSI wyjsc Re lambda > 0. FAIL => STOP.
C3 (proznia analitycznie): sigma(JL) wokol u=1 sektora tachionowego vs
   lambda^2(kappa) = -kappa^2(kappa^2-gamma)/4 (Phase 1), blad <= 1e-3
   w metryce |d(lambda^2)|/max(|lambda^2_ex|,1). FAIL => STOP.

[INPUT-ONTO] zlozenie zespolone u z g: u0 = g_d(x) (modul pola = g;
LOCK sec. 2 Rejestr WEJSC) -- flagowane w kazdym bloku wynikow.
"""
import numpy as np
from scipy.linalg import eig, eigh, solve_banded

BETA = 1.0
GAMMA = 1.0
BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-symplectic-Jspectrum-2026-08-31/")
NPZ = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
       "op-bloch-chain-stability-2026-08-31/Phase2_backgrounds.npz")
ANCHOR = -1.222191
ANCHOR_TOL = 1e-4
LAMBDA_BAND = 12.0   # ruling tol (Phase_method_decisions sec. 4)
GATE_C3 = 1e-3
LBOX, NBOX = (40.0, (1024, 2048))

FAILS = []


# ------------------------------------------------------------ maszyneria
def kin_matrix(K, h, d, k):
    """Macierz -d/dx(K d/dx) z faza Blocha e^{ikd} (stencil bloch)."""
    N = len(K)
    wmid = 0.5 * (K + np.roll(K, -1))
    dtype = float if k == 0.0 else complex
    A = np.zeros((N, N), dtype=dtype)
    idx = np.arange(N)
    A[idx, idx] = (wmid + np.roll(wmid, 1)) / h ** 2
    off = -wmid / h ** 2
    A[idx[:-1], idx[:-1] + 1] = off[:-1]
    A[idx[:-1] + 1, idx[:-1]] = off[:-1]
    ph = np.exp(1j * k * d) if k != 0.0 else 1.0
    A[N - 1, 0] = off[N - 1] * ph
    A[0, N - 1] = off[N - 1] * np.conj(ph) if k != 0.0 else off[N - 1]
    return A

def onshell_Qminus(g, K, h):
    """Q_minus = (K g')'/g dyskretnie: -(A_kin(k=0) g)/g -- L_minus g = 0
    dokladnie w fp (tozsamosc P1c; Phase_method_decisions sec. 2)."""
    A0 = kin_matrix(K, h, 1.0, 0.0)
    return -(A0 @ g) / g

def build_LpLm(g, h, d, k, model):
    """Operatory L_plus(k), L_minus(k) dla modelu 'canonical'/'nls3'/'nls7'."""
    if model == 'canonical':
        K = g ** 4
        gp = (np.roll(g, -1) - np.roll(g, 1)) / (2 * h)
        Qp = K * (2 * BETA * g - 3 * GAMMA * g ** 2 + 2 * gp ** 2 / g ** 2)
    elif model == 'nls3':
        K = 2.0 * np.ones_like(g)
        Qp = 2.0 - 6.0 * g ** 2          # U = g^2 - g^4/2
    elif model == 'nls7':
        K = 2.0 * np.ones_like(g)
        Qp = 2.0 - 14.0 * g ** 6         # U = g^2 - g^8/4
    else:
        raise ValueError(model)
    Qm = onshell_Qminus(g, K, h)
    Akin = kin_matrix(K, h, d, k)
    Lp = Akin + np.diag(Qp)
    Lm = Akin + np.diag(Qm)
    return Lp, Lm

def J_spectrum(Lp, Lm):
    """Pelne sigma(JL): M = [[0, Lm/2], [-Lp/2, 0]], eig gesty, sort Re."""
    N = Lp.shape[0]
    M = np.zeros((2 * N, 2 * N), dtype=Lp.dtype)
    M[:N, N:] = 0.5 * Lm
    M[N:, :N] = -0.5 * Lp
    lam = eig(M, right=False)
    return lam[np.argsort(-lam.real)]

def quad_symmetry(lams):
    """Diagnostyka symetrii lambda -> -lambda: max_i min_j |l_i + l_j|.
    [correction note 2026-08-31: pierwotny wzor sort(l)+sort(-l) byl
    tozsamosciowo 2*sort(l) i niczego nie testowal]"""
    out = 0.0
    for i0 in range(0, len(lams), 256):
        blk = lams[i0:i0 + 256]
        out = max(out, float(np.max(np.min(
            np.abs(blk[:, None] + lams[None, :]), axis=1))))
    return out

def tol_readings(lams):
    """(tol PRIMARY pasmowo-ograniczony |lambda|<=12, tol literalny)."""
    im_all = float(np.max(np.abs(lams.imag)))
    inband = lams[np.abs(lams) <= LAMBDA_BAND]
    im_band = float(np.max(np.abs(inband.imag))) if len(inband) else 0.0
    return (max(1e-6, 1e-3 * im_band), max(1e-6, 1e-3 * im_all))

def product_check(Lp, Lm, lams):
    """Cross-check LOCKa: lambda^2 = -nu/4, nu = spec(L_minus L_plus)."""
    nu = eig(Lm @ Lp, right=False)
    lam2_from_nu = np.sqrt((-nu / 4.0).astype(complex))
    maxre_nu = float(np.max(lam2_from_nu.real))
    maxre_direct = float(np.max(lams.real))
    return maxre_direct, maxre_nu, abs(maxre_direct - maxre_nu)

def newton_soliton(p, N, L):
    """Newton na pol-komorce parzystej dla D2 phi - phi + phi^p = 0."""
    h = L / N
    M = N // 2
    xc = np.arange(M + 1) * h - L / 2.0
    if p == 3:
        phi = np.sqrt(2.0) / np.cosh(xc)
    else:
        phi = 4.0 ** (1.0 / 6.0) * (1.0 / np.cosh(3.0 * xc)) ** (1.0 / 3.0)
    seed = phi.copy()
    for it in range(60):
        phim = np.empty(M + 3)
        phim[1:-1] = phi
        phim[0] = phi[1]        # lustro w j=0
        phim[-1] = phi[-2]      # lustro w j=M
        d2 = (phim[2:] - 2 * phi + phim[:-2]) / h ** 2
        R = d2 - phi + phi ** p
        if np.max(np.abs(R)) < 1e-13:
            break
        ab = np.zeros((3, M + 1))
        ab[1] = -2 / h ** 2 - 1 + p * phi ** (p - 1)
        ab[0, 1:] = 1 / h ** 2
        ab[2, :-1] = 1 / h ** 2
        ab[0, 1] = 2 / h ** 2   # podwojone sprzezenie lustrzane
        ab[2, M - 1] = 2 / h ** 2
        phi = phi + solve_banded((1, 1), ab, -R)
    full = np.empty(N)
    full[:M + 1] = phi
    full[M + 1:] = phi[M - 1:0:-1]
    d2f = (np.roll(full, -1) - 2 * full + np.roll(full, 1)) / h ** 2
    res = float(np.max(np.abs(d2f - full + full ** p)))
    seed_res = float(np.max(np.abs(
        (np.roll(np.concatenate([seed, seed[M - 1:0:-1]]), -1)
         - 2 * np.concatenate([seed, seed[M - 1:0:-1]])
         + np.roll(np.concatenate([seed, seed[M - 1:0:-1]]), 1)) / h ** 2
        - np.concatenate([seed, seed[M - 1:0:-1]])
        + np.concatenate([seed, seed[M - 1:0:-1]]) ** p)))
    return full, res, seed_res


print("=" * 78)
print("Phase 2 -- bramka maszynerii (LOCK sec. 3 Phase 2; ruling tol:")
print("Phase_method_decisions sec. 4 -- oba odczyty raportowane).")
print("=" * 78)

# ------------------------------------------------- Gate A: kotwica L_plus
print("\n[Gate A] reprodukcja L_plus cyklu bloch (problem wazony, tlo 3pi,")
print("  N=800, k=0): lambda_min vs kotwica %.6f +/- 1e-4" % ANCHOR)
print("  [INPUT-ONTO] u0 = g_d (modul pola = g; LOCK sec. 2)")
data = np.load(NPZ)
g3 = data["3pi__0.7__N800"]
d3 = float(data["3pi__0.7__d"][0])
h3 = d3 / len(g3)
Lp3, Lm3 = build_LpLm(g3, h3, d3, 0.0, 'canonical')
Kw = g3 ** 4
sw = np.sqrt(Kw)
H = Lp3 / sw[:, None] / sw[None, :]
lam_min = float(eigh(H, subset_by_index=(0, 0), eigvals_only=True)[0])
okA = abs(lam_min - ANCHOR) <= ANCHOR_TOL
print("  lambda_min(3pi, wazony) = %+.6f (|d|=%.2e) -> %s"
      % (lam_min, abs(lam_min - ANCHOR), "PASS" if okA else "FAIL"))
lmg = float(np.max(np.abs(Lm3 @ g3)))
lpgp = float(np.max(np.abs(
    Lp3 @ ((np.roll(g3, -1) - np.roll(g3, 1)) / (2 * h3)))))
print("  tozsamosci dyskretne: ||L_minus g||_inf = %.2e (mod fazowy, fp-0);"
      % lmg)
print("  ||L_plus g'||_inf = %.2e (translacja, oczekiwane O(h^2))" % lpgp)
if not okA:
    FAILS.append("GateA")

# --------------------------------------------------------------- C1 / C2
for tag, model, p in (("C1 (NLS kubiczny -- znany STABILNY)", 'nls3', 3),
                      ("C2 (NLS |u|^6 u -- znany NIESTABILNY)", 'nls7', 7)):
    print("\n[%s] pudelko L=%.0f, soliton Newton-polish (pol-komorka"
          % (tag, LBOX))
    print("  parzysta), operatory K=2, U wchloniete omega=1.")
    maxre = {}
    for N in NBOX:
        h = LBOX / N
        phi, res, seed_res = newton_soliton(p, N, LBOX)
        Lp, Lm = build_LpLm(phi, h, LBOX, 0.0, model)
        lams = J_spectrum(Lp, Lm)
        tolP, tolL = tol_readings(lams)
        mr = float(np.max(lams.real))
        maxre[N] = (mr, tolP, tolL)
        quad = quad_symmetry(lams)
        print("  N=%4d: seed-residuum=%.1e, Newton-residuum=%.1e;"
              % (N, seed_res, res))
        print("          max Re lambda = %+.6e; tol PRIMARY=%.3e,"
              " tol literalny=%.3e" % (mr, tolP, tolL))
        print("          symetria czworek +/-lambda: %.1e" % quad)
        if N == NBOX[0]:
            mrd, mrnu, dd = product_check(Lp, Lm, lams)
            print("          cross-check produktowy: max Re lambda(JL)="
                  "%.6e vs z L_-L_+: %.6e (|d|=%.1e)" % (mrd, mrnu, dd))
    mr0, tolP0, tolL0 = maxre[NBOX[0]]
    mr1 = maxre[NBOX[1]][0]
    dconv = abs(mr1 - mr0)
    thr = 0.01 * max(mr1, 0.01)
    if model == 'nls3':
        okP = mr0 <= tolP0
        okL = mr0 <= tolL0
        print("  C1 werdykt: max Re lambda(N=%d) = %.3e <= tol?"
              % (NBOX[0], mr0))
        print("    odczyt PRIMARY (pasmowy, tol=%.3e): %s;"
              " odczyt literalny (tol=%.3e): %s"
              % (tolP0, "PASS" if okP else "FAIL",
                 tolL0, "PASS" if okL else "FAIL"))
        print("    konsystencja siatek: |d max Re| = %.2e" % dconv)
        if not okP:
            FAILS.append("C1")
        c1_mr = mr1
    else:
        okC2 = mr0 > 0 and mr1 > 0 and dconv <= thr
        sep = mr1 / max(c1_mr, 1e-300)
        print("  C2 werdykt: Re lambda > 0? max Re lambda = %.6f (N=%d)"
              " / %.6f (N=%d), |d|=%.2e (prog zb. %.2e) -> %s"
              % (mr0, NBOX[0], mr1, NBOX[1], dconv, thr,
                 "PASS" if okC2 else "FAIL"))
        print("    separacja od C1 (deskryptywnie): %.1e x (wymog >= 10)"
              % sep)
        if not okC2:
            FAILS.append("C2")

# --------------------------------------------------------------------- C3
print("\n[C3] proznia u=1 sektora tachionowego, komorka d=2pi, 9 punktow k:")
print("  analityka (Phase 1): lambda^2(kappa) = -kappa^2(kappa^2-gamma)/4")
D_VAC = 2 * np.pi
maxerr_c3 = {}
for N in (400, 800):
    h = D_VAC / N
    g = np.ones(N)
    Lp0, Lm0 = None, None
    maxerr = 0.0
    for k in np.linspace(0.0, np.pi / D_VAC, 9):
        Lp, Lm = build_LpLm(g, h, D_VAC, float(k), 'canonical')
        P = (-0.25) * (Lm @ Lp)          # spec(P) = lambda^2
        lam2 = eig(P, right=False)
        kap = float(k) + 2 * np.pi * np.arange(-12, 13) / D_VAC
        lam2_ex = -0.25 * kap ** 2 * (kap ** 2 - GAMMA)
        order = np.argsort(np.abs(lam2_ex))
        ex8 = lam2_ex[order[:8]]
        num_sorted = lam2[np.argsort(np.abs(lam2))]
        # dopasowanie: dla kazdej z 8 dokladnych najblizsza numeryczna
        for ex in ex8:
            j = int(np.argmin(np.abs(num_sorted - ex)))
            err = abs(num_sorted[j] - ex) / max(abs(ex), 1.0)
            maxerr = max(maxerr, float(err))
    maxerr_c3[N] = maxerr
    print("  N=%3d: maxerr(8 galezi, wszystkie k) = %.3e -> %s"
          % (N, maxerr, "PASS" if maxerr <= GATE_C3 else "FAIL"))
okC3 = maxerr_c3[400] <= GATE_C3
ratio = maxerr_c3[400] / max(maxerr_c3[800], 1e-300)
print("  C3 werdykt (siatka bazowa N=400): %s; ratio N400/N800 = %.2f"
      " (oczekiwany ~4, rzad 2)" % ("PASS" if okC3 else "FAIL", ratio))
print("  Uwaga fizyczna (z Phase 1): 0<kappa<1 => Re lambda>0 -- proznia")
print("  tachionowa jest niestabilna TAKZE symplektycznie (max gamma/4).")
if not okC3:
    FAILS.append("C3")

# ----------------------------------------------------------------- werdykt
print("\n" + "=" * 78)
if FAILS:
    print("PHASE 2 GATE: FAIL (%s) => STOP -- kod niewazny (LOCK sec. 3)."
          % ", ".join(FAILS))
    raise SystemExit(1)
print("PHASE 2 GATE: PASS (Gate A kotwica + C1 + C2 + C3).")
print("Maszyneria wazna: przechodzi znany-stabilny (C1) MIMO ujemnego")
print("kierunku L_plus, wykrywa znany-niestabilny (C2), odtwarza")
print("analityke prozni (C3) i operator Hessianu cyklu bloch (Gate A).")
print("=" * 78)
