#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-bloch-chain-stability (Phase 3) -- RACHUNEK CENTRALNY: spektrum pasmowe
Blocha wokol samouzgodnionych tel lancucha (sektor tachionowy, akcja
kanoniczna K=g^4, U=g^7/7-g^8/8).

LOCK: Phase0_balance.md sec. 3 Phase 3; decyzje FROZEN:
Phase_method_decisions.md sec. 2, 4, 5.

P3a: omega^2_min(d) = min_k lambda_min(L_d(k)); operator (druga wariacja):
     -(K phi')' + Q phi = omega^2 K phi,  Q = K(2g - 3g^2 + 2g'^2/g^2);
     Bloch phi(x+d) = e^{ikd} phi(x); siatki k: 16 i 32 pkt w [0, pi/d];
     siatki przestrzenne N in {400, 800}. Wartosc glowna: (N=800, k=32).
     Zbieznosc: Delta_siatka i Delta_k <= 0.01 max(|omega^2_min|, 0.1).
     Mod translacyjny (Goldstone, brak zrodla) identyfikowany w k=0
     overlapem B-wazonym z g0' (prog 0.9; regula FROZEN przed startem);
     kryterium LOCKa stosowane LITERALNIE do pelnego minimum, deskryptywnie
     raportowane tez minimum bez modu translacyjnego.
     Cross-check kanoniczny: u=g^3/3 => -chi'' + (4g-5g^2) chi = omega^2 chi
     (waga 1) w k in {0, pi/d}, N=800 -- zgodnosc raportowana.

P3b (kontrola artefaktu, NIEUSUWALNA): proznia g=1 w superkomorce 4d dla
     kazdego zalockowanego d: 8 najnizszych galezi MUSI odtworzyc
     omega^2 = k^2 - 1 zwiniete do strefy [0, pi/(4d)]; blad wzgledny
     |num-ex|/max(|ex|,1) <= 1e-3.

P3c (kontrola sektora stabilnego, NIEUSUWALNA): pelna procedura
     (relaksacja Newtona + Bloch) w sektorze stabilnym m^2=+1 ze zrodlem
     gaussowskim (konstrukcja Q2 poprzednika: d_s=4, sigma_s=0.5, q=1):
     omega^2_min MUSI byc > 0. Kotwica: lambda_min(k=0) vs +1.88310.
"""
import json
import numpy as np
from scipy.linalg import eigh
import scipy.sparse as sp
import scipy.sparse.linalg as spla

BETA = 1.0
GAMMA = 1.0
BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-bloch-chain-stability-2026-08-31/")
NPZ = BASE + "Phase2_backgrounds.npz"
OUT_JSON = BASE + "Phase3_results.json"
NK = (16, 32)
GATE_VAC = 1e-3


def bloch_weighted(g, h, d, k, nev, sector='tach', vecs_out=False):
    """Problem uogolniony -(K phi')' + Q phi = omega^2 K phi (Bloch)."""
    N = len(g)
    gp = (np.roll(g, -1) - np.roll(g, 1)) / (2 * h)
    if sector == 'tach':
        w = g ** 4
        pot = w * (2 * BETA * g - 3 * GAMMA * g ** 2 + 2 * gp ** 2 / g ** 2)
    else:  # 'stab' -- operator Q2 poprzednika
        w = g ** 4
        pot = -w * (2 * BETA * g - 3 * GAMMA * g ** 2 - 2 * gp ** 2 / g ** 2)
    wmid = 0.5 * (w + np.roll(w, -1))
    diag = (wmid + np.roll(wmid, 1)) / h ** 2 + pot
    A = np.zeros((N, N), dtype=complex)
    idx = np.arange(N)
    A[idx, idx] = diag
    off = -wmid / h ** 2
    A[idx[:-1], idx[:-1] + 1] = off[:-1]
    A[idx[:-1] + 1, idx[:-1]] = off[:-1]
    ph = np.exp(1j * k * d)
    A[N - 1, 0] = off[N - 1] * ph
    A[0, N - 1] = off[N - 1] * np.conj(ph)
    sw = np.sqrt(w)
    H = A / sw[:, None] / sw[None, :]
    if vecs_out:
        vals, vecs = eigh(H, subset_by_index=(0, nev - 1))
        return vals, vecs / sw[:, None]
    return eigh(H, subset_by_index=(0, nev - 1), eigvals_only=True)


def bloch_uform(g, h, d, k, nev):
    """Cross-check kanoniczny: -chi'' + (4g-5g^2) chi = omega^2 chi."""
    N = len(g)
    pot = 4 * BETA * g - 5 * GAMMA * g ** 2
    A = np.zeros((N, N), dtype=complex)
    idx = np.arange(N)
    A[idx, idx] = 2 / h ** 2 + pot
    A[idx[:-1], idx[:-1] + 1] = -1 / h ** 2
    A[idx[:-1] + 1, idx[:-1]] = -1 / h ** 2
    ph = np.exp(1j * k * d)
    A[N - 1, 0] = -ph / h ** 2
    A[0, N - 1] = -np.conj(ph) / h ** 2
    return eigh(A, subset_by_index=(0, nev - 1), eigvals_only=True)


def scan_w2min(g, h, d, nk):
    ks = np.linspace(0.0, np.pi / d, nk)
    vmin, kmin = np.inf, 0.0
    for k in ks:
        v = bloch_weighted(g, h, d, float(k), 1)[0]
        if v < vmin:
            vmin, kmin = float(v), float(k)
    return vmin, kmin


print("=" * 78)
print("Phase 3 -- spektrum pasmowe Blocha (sektor tachionowy, akcja")
print("kanoniczna). LOCK sec. 3 Phase 3; operator i reguly FROZEN w")
print("Phase_method_decisions.md sec. 2/4/5.")
print("=" * 78)

data = np.load(NPZ)
bgs = []
for key in sorted(data.files):
    if key.endswith("__d"):
        dn, tag, _ = key.split("__")
        bgs.append((dn, tag, float(data[key][0])))

# --------------------------------------------------------- P3b (nieusuwalna)
print("\n[P3b] KONTROLA: proznia g=1 w superkomorce 4d -- zwiniete k^2-1")
D_ALL = {"pi": np.pi, "2pi": 2 * np.pi, "3pi": 3 * np.pi,
         "4pi": 4 * np.pi, "6pi": 6 * np.pi}
p3b_pass = True
for dn, d in D_ALL.items():
    dsc = 4 * d
    Nsc = 1600
    h = d / 400
    g = np.ones(Nsc)
    ks = np.linspace(0.0, np.pi / dsc, 8)
    maxerr = 0.0
    for k in ks:
        vals = bloch_weighted(g, h, dsc, float(k), 8)
        kap = k + 2 * np.pi * np.arange(-10, 11) / dsc
        ex = np.sort(kap ** 2 - 1.0)[:8]
        maxerr = max(maxerr, float(np.max(
            np.abs(vals - ex) / np.maximum(np.abs(ex), 1.0))))
    ok = maxerr <= GATE_VAC
    print("  d=%s (superkomorka %.4f, N=1600): maxerr(8 galezi) = %.3e -> %s"
          % (dn, dsc, maxerr, "PASS" if ok else "FAIL"))
    if not ok:
        p3b_pass = False

# --------------------------------------------------------- P3c (nieusuwalna)
print("\n[P3c] KONTROLA: sektor stabilny m^2=+1 ze zrodlem (konstrukcja Q2")
print("  poprzednika: d_s=4, sigma_s=0.5, q=1); wymog: omega^2_min > 0.")


def rho_lattice(x, d, sigma):
    r = np.zeros_like(x)
    for m in range(-6, 7):
        r += np.exp(-0.5 * ((x - d / 2 - m * d) / sigma) ** 2)
    return r / (np.sqrt(2 * np.pi) * sigma)


def newton_stab(N, h, qrho):
    """Newton periodyczny sektora stabilnego (verbatim wzorzec Q2)."""
    psi = np.ones(N)
    idx = np.arange(N)
    for frac in (0.25, 0.5, 0.75, 1.0):
        qr = frac * qrho
        ok = False
        for it in range(80):
            ip, im = np.roll(psi, -1), np.roll(psi, 1)
            d2 = (ip - 2 * psi + im) / h ** 2
            d1 = (ip - im) / (2 * h)
            R = d2 + 2 * d1 ** 2 / psi + BETA * psi ** 2 \
                - GAMMA * psi ** 3 + qr
            nrm = float(np.max(np.abs(R)))
            tol = max(1e-10, 20 * np.finfo(float).eps
                      * float(np.max(psi)) / h ** 2)
            if nrm < tol:
                ok = True
                break
            dg = -2 / h ** 2 - 2 * d1 ** 2 / psi ** 2 \
                + 2 * BETA * psi - 3 * GAMMA * psi ** 2
            up = 1 / h ** 2 + 4 * d1 / psi / (2 * h)
            lo = 1 / h ** 2 - 4 * d1 / psi / (2 * h)
            J = sp.csr_matrix(
                (np.concatenate([dg, up, lo]),
                 (np.concatenate([idx, idx, idx]),
                  np.concatenate([idx, (idx + 1) % N, (idx - 1) % N]))),
                shape=(N, N))
            dpsi = spla.spsolve(J, -R)
            lam = 1.0
            for _ in range(40):
                trial = psi + lam * dpsi
                ip, im = np.roll(trial, -1), np.roll(trial, 1)
                d2 = (ip - 2 * trial + im) / h ** 2
                d1t = (ip - im) / (2 * h)
                Rt = d2 + 2 * d1t ** 2 / trial + BETA * trial ** 2 \
                    - GAMMA * trial ** 3 + qr
                if float(np.min(trial)) > 1e-6 and \
                        float(np.max(np.abs(Rt))) < nrm:
                    break
                lam *= 0.5
            psi = psi + lam * dpsi
        if not ok:
            return None
    return psi


DS, SIG, QQ = 4.0, 0.5, 1.0
p3c_vals = {}
for N in (400, 800):
    h = DS / N
    x = np.arange(N) * h
    psi = newton_stab(N, h, QQ * rho_lattice(x, DS, SIG))
    if psi is None:
        print("  N=%d: Newton NIE ZBIEGL -> P3c FAIL" % N)
        p3c_vals[N] = None
        continue
    ks = np.linspace(0.0, np.pi / DS, 16)
    vmin, kmin = np.inf, 0.0
    for k in ks:
        v = bloch_weighted(psi, h, DS, float(k), 1, sector='stab')[0]
        if v < vmin:
            vmin, kmin = float(v), float(k)
    v_k0 = float(bloch_weighted(psi, h, DS, 0.0, 1, sector='stab')[0])
    p3c_vals[N] = vmin
    print("  N=%d: psi in [%.6f, %.6f]; omega^2_min = %+.5f (argmin k=%.4f);"
          " lambda_min(k=0) = %+.5f (kotwica poprzednika: +1.88310)"
          % (N, float(np.min(psi)), float(np.max(psi)), vmin, kmin, v_k0))
p3c_pass = all(v is not None and v > 0 for v in p3c_vals.values())
print("  P3c: %s (omega^2_min > 0 na obu siatkach: %s)"
      % ("PASS" if p3c_pass else "FAIL",
         ", ".join("%s" % v for v in p3c_vals.values())))

# ------------------------------------------------------------- P3a centralny
print("\n[P3a] spektra tel Phase 2 (tla: %s)"
      % ", ".join("%s/%s" % (dn, tag) for dn, tag, _ in bgs))
res_out = {}
for dn, tag, d in bgs:
    label = "%s/%s" % (dn, tag)
    print("\n  --- tlo %s (d=%.6f) ---" % (label, d))
    vals_cfg = {}
    argk = {}
    for N in (400, 800):
        g = data["%s__%s__N%d" % (dn, tag, N)]
        h = d / N
        for nk in NK:
            v, kmin = scan_w2min(g, h, d, nk)
            vals_cfg[(N, nk)] = v
            argk[(N, nk)] = kmin
            print("    N=%3d, k-siatka %2d: omega^2_min = %+.6f "
                  "(argmin k = %.6f)" % (N, nk, v, kmin))
    v_main = vals_cfg[(800, 32)]
    d_sp = abs(vals_cfg[(400, 32)] - v_main)
    d_k = abs(vals_cfg[(800, 16)] - v_main)
    tol = 0.01 * max(abs(v_main), 0.1)
    conv = d_sp <= tol and d_k <= tol
    print("    GLOWNA (N=800,k32): omega^2_min = %+.6f; Delta_siatka=%.2e,"
          " Delta_k=%.2e, prog=%.2e -> %s"
          % (v_main, d_sp, d_k, tol, "ZBIEZNE" if conv else "NIEZBIEZNE"))
    # k=0: spektrum + identyfikacja modu translacyjnego (regula FROZEN)
    g8 = data["%s__%s__N800" % (dn, tag)]
    h8 = d / 800
    vals0, vecs0 = bloch_weighted(g8, h8, d, 0.0, 4, vecs_out=True)
    gp = (np.roll(g8, -1) - np.roll(g8, 1)) / (2 * h8)
    w = g8 ** 4
    ovl = []
    for j in range(4):
        phi = np.real(vecs0[:, j])
        num = abs(float(np.sum(w * phi * gp) * h8))
        den = np.sqrt(float(np.sum(w * phi ** 2) * h8)
                      * float(np.sum(w * gp ** 2) * h8))
        ovl.append(num / den if den > 0 else 0.0)
    trans_idx = [j for j in range(4) if ovl[j] >= 0.9]
    nontrans = [float(vals0[j]) for j in range(4) if j not in trans_idx]
    print("    k=0 spektrum: %s" % ["%+.6f" % v for v in vals0])
    print("    overlap z g0' (B-wazony): %s -> mod translacyjny: %s"
          % (["%.3f" % o for o in ovl],
             ("indeks %s, lambda_trans = %s"
              % (trans_idx, ["%+.2e" % float(vals0[j]) for j in trans_idx]))
             if trans_idx else "NIEZIDENTYFIKOWANY w 4 najnizszych"))
    if nontrans:
        print("    deskryptywnie (nie-kryterialnie): min k=0 bez translacji"
              " = %+.6f" % nontrans[0])
    # cross-check u-forma (waga 1) w k in {0, pi/d}, N=800
    for k in (0.0, float(np.pi / d)):
        vw = float(bloch_weighted(g8, h8, d, k, 1)[0])
        vu = float(np.real(bloch_uform(g8, h8, d, k, 1)[0]))
        print("    cross-check u-forma k=%.4f: waga-K %+.6f vs u-forma"
              " %+.6f (|delta|=%.2e)" % (k, vw, vu, abs(vw - vu)))
    res_out[label] = dict(d=d, w2min=v_main, conv=bool(conv),
                          d_sp=d_sp, d_k=d_k, tol=tol,
                          argmin_k=argk[(800, 32)],
                          lam_trans=[float(vals0[j]) for j in trans_idx])

with open(OUT_JSON, "w") as f:
    json.dump(res_out, f, indent=1)
print("\nzapisano wyniki: %s" % OUT_JSON)

# ------------------------------------------------------------------ tabela
print("\n" + "=" * 78)
print("TABELA omega^2_min(d) (kryterium zbieznosci LOCK, doslownie):")
print("  P3b (proznia w superkomorce): %s" % ("PASS" if p3b_pass else "FAIL"))
print("  P3c (sektor stabilny): %s" % ("PASS" if p3c_pass else "FAIL"))
for label, r in res_out.items():
    print("  tlo %-14s d=%8.4f: omega^2_min = %+.6f (%s; argmin k=%.4f)"
          % (label, r["d"], r["w2min"],
             "ZBIEZNE" if r["conv"] else "NIEZBIEZNE", r["argmin_k"]))
any_pos = any(r["w2min"] > 0 and r["conv"] for r in res_out.values())
all_neg = all(r["w2min"] < 0 and r["conv"] for r in res_out.values())
print("\n  istnieje d z omega^2_min > 0 (zbiezne): %s" % any_pos)
print("  omega^2_min < 0 dla wszystkich policzalnych tel (zbieznie): %s"
      % all_neg)
print("  (werdykt Q wymaga Phase 4 -- kryteria LOCK sec. 3 Phase 4)")
print("=" * 78)
