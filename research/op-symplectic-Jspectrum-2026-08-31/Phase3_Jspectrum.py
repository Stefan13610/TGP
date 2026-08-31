#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-symplectic-Jspectrum (Phase 3) -- RACHUNEK CENTRALNY: sigma(JL) dla
4 tel lancucha sektora tachionowego (akcja kanoniczna K=g^4,
U=g^7/7-g^8/8; dynamika symplektyczna I rzedu i u_t = dE/du*).

LOCK: Phase0_balance.md sec. 2-3; decyzje FROZEN: Phase_method_decisions.md
sec. 1-5, 7. Operatory zweryfikowane w Phase 1 (sympy-exact) i Phase 2
(Gate A kotwica, C1/C2/C3 PASS).

Tla WYLACZNIE z Phase2_backgrounds.npz cyklu bloch (4: 3pi/0.7, 4pi/0.7,
6pi/0.7 2-garb, 6pi/S3single 1-garb). k=0 PRIMARY + skan 9 punktow
w [0, pi/d]; siatki N in {400, 800}; zbieznosc per k:
|d max Re lambda| <= 0.01 max(max Re lambda, 0.01).
tol: max(1e-6, 1e-3 max|Im lambda|) na siatce bazowej -- odczyt PRIMARY
pasmowy (|lambda| <= 12) ORAZ literalny-nieograniczony, oba raportowane
(ruling FROZEN: Phase_method_decisions sec. 4).

Werdykt per tlo i zbiorczy: Q-PASS / Q-FAIL / MIXED / Q-INCONCLUSIVE
(LOCK sec. 3 Phase 3, doslownie).

[INPUT-ONTO] zlozenie zespolone u z g: u0 = g_d(x) (modul pola = g;
LOCK sec. 2 Rejestr WEJSC) -- flaga w kazdym bloku wynikow.
"""
import json
import numpy as np
from scipy.linalg import eig

BETA = 1.0
GAMMA = 1.0
BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-symplectic-Jspectrum-2026-08-31/")
NPZ = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
       "op-bloch-chain-stability-2026-08-31/Phase2_backgrounds.npz")
OUT_JSON = BASE + "Phase3_results.json"
LAMBDA_BAND = 12.0
NKPTS = 9


def kin_matrix(K, h, d, k):
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


def build_LpLm(g, h, d, k):
    K = g ** 4
    gp = (np.roll(g, -1) - np.roll(g, 1)) / (2 * h)
    Qp = K * (2 * BETA * g - 3 * GAMMA * g ** 2 + 2 * gp ** 2 / g ** 2)
    A0 = kin_matrix(K, h, 1.0, 0.0)
    Qm = -(A0 @ g) / g          # on-shell: L_minus g = 0 exact (P1c)
    Akin = kin_matrix(K, h, d, k)
    return Akin + np.diag(Qp), Akin + np.diag(Qm)


def J_matrix(Lp, Lm):
    N = Lp.shape[0]
    M = np.zeros((2 * N, 2 * N), dtype=Lp.dtype)
    M[:N, N:] = 0.5 * Lm
    M[N:, :N] = -0.5 * Lp
    return M


def quad_symmetry(lams):
    out = 0.0
    for i0 in range(0, len(lams), 256):
        blk = lams[i0:i0 + 256]
        out = max(out, float(np.max(np.min(
            np.abs(blk[:, None] + lams[None, :]), axis=1))))
    return out


def tol_readings(lams_list):
    """Z listy widm (wszystkie k, siatka bazowa): (tol PRIMARY, literalny)."""
    im_band, im_all = 0.0, 0.0
    for lams in lams_list:
        im_all = max(im_all, float(np.max(np.abs(lams.imag))))
        inb = lams[np.abs(lams) <= LAMBDA_BAND]
        if len(inb):
            im_band = max(im_band, float(np.max(np.abs(inb.imag))))
    return max(1e-6, 1e-3 * im_band), max(1e-6, 1e-3 * im_all)


print("=" * 78)
print("Phase 3 -- RACHUNEK CENTRALNY: sigma(JL) na tlach lancucha")
print("(LOCK sec. 3 Phase 3; operatory FROZEN: Phase_method_decisions 1-2;")
print("bramki Phase 1/2: PASS). [INPUT-ONTO] u0 = g_d (modul pola = g).")
print("=" * 78)

data = np.load(NPZ)
BGS = [("3pi", "0.7"), ("4pi", "0.7"), ("6pi", "0.7"), ("6pi", "S3single")]

res_out = {}
verdicts_P = {}
verdicts_L = {}
for dn, tag in BGS:
    label = "%s/%s" % (dn, tag)
    d = float(data["%s__%s__d" % (dn, tag)][0])
    ks = np.linspace(0.0, np.pi / d, NKPTS)
    print("\n--- tlo %s (d=%.6f) --- [INPUT-ONTO] u0 = g_d" % (label, d))
    maxre = {}
    spectra_base = []
    diag_k0 = {}
    for N in (400, 800):
        g = data["%s__%s__N%d" % (dn, tag, N)]
        h = d / N
        for k in ks:
            Lp, Lm = build_LpLm(g, h, d, float(k))
            M = J_matrix(Lp, Lm)
            if k == 0.0:
                lams, vecs = eig(M, right=True)
            else:
                lams, vecs = eig(M, right=False), None
            maxre[(N, float(k))] = float(np.max(lams.real))
            if N == 400:
                spectra_base.append(lams)
            if k == 0.0:
                # diagnostyka k=0: tozsamosci, symetria, produkt, mody
                gp = (np.roll(g, -1) - np.roll(g, 1)) / (2 * h)
                lmg = float(np.max(np.abs(Lm @ g)))
                lpgp = float(np.max(np.abs(Lp @ gp)))
                qs = quad_symmetry(lams)
                P = -0.25 * (Lm @ Lp)
                lam2 = eig(P, right=False)
                mr_prod = float(np.max(np.sqrt(lam2.astype(complex)).real))
                # mody symetrii (regula FROZEN sec. 5): overlap >= 0.9
                idx_pos = np.where(lams.real > 1e-3)[0]
                excl = []
                nrm_g = np.sqrt(float(np.sum(g ** 2)))
                nrm_gp = np.sqrt(float(np.sum(gp ** 2)))
                for j in idx_pos:
                    a_c, b_c = vecs[:N, j], vecs[N:, j]
                    na = np.sqrt(float(np.sum(np.abs(a_c) ** 2)))
                    nb = np.sqrt(float(np.sum(np.abs(b_c) ** 2)))
                    o_tr = (abs(np.sum(np.conj(a_c) * gp)) / (na * nrm_gp)
                            if na > 0 else 0.0)
                    o_ph = (abs(np.sum(np.conj(b_c) * g)) / (nb * nrm_g)
                            if nb > 0 else 0.0)
                    if max(o_tr, o_ph) >= 0.9:
                        excl.append(j)
                re_sorted = np.sort(lams.real)[::-1]
                mask = np.ones(len(lams), dtype=bool)
                mask[excl] = False
                mr_excl = float(np.max(lams[mask].real)) if mask.any() \
                    else 0.0
                diag_k0[N] = dict(lmg=lmg, lpgp=lpgp, qs=qs,
                                  mr_prod=mr_prod, mr_excl=mr_excl,
                                  n_sym=len(excl),
                                  top4=[float(v) for v in re_sorted[:4]])
    tolP, tolL = tol_readings(spectra_base)
    # tabela per k + zbieznosc
    all_conv = True
    any_fail_P, any_fail_L = False, False
    all_stab_P, all_stab_L = True, True
    rows = []
    for k in ks:
        m4, m8 = maxre[(400, float(k))], maxre[(800, float(k))]
        dlt = abs(m8 - m4)
        thr = 0.01 * max(m8, 0.01)
        conv = dlt <= thr
        all_conv = all_conv and conv
        stab_P, stab_L = m8 <= tolP, m8 <= tolL
        if conv and not stab_P:
            any_fail_P = True
        if conv and not stab_L:
            any_fail_L = True
        all_stab_P = all_stab_P and stab_P and conv
        all_stab_L = all_stab_L and stab_L and conv
        rows.append(dict(k=float(k), maxre_N400=m4, maxre_N800=m8,
                         delta=dlt, thr=thr, conv=bool(conv)))
        print("  k=%.6f: max Re lambda = %+.6f (N400) / %+.6f (N800);"
              " |d|=%.2e (prog %.2e) -> %s"
              % (k, m4, m8, dlt, thr, "ZBIEZNE" if conv else "NIEZBIEZNE"))
    print("  tol (siatka bazowa N=400): PRIMARY(pasmowy |l|<=12) = %.4e,"
          " literalny = %.4e" % (tolP, tolL))
    for N in (400, 800):
        dg = diag_k0[N]
        print("  k=0, N=%d: ||L_-g||=%.1e (fp-0), ||L_+g'||=%.1e (O(h^2));"
              " symetria +/-lambda: %.1e" % (N, dg["lmg"], dg["lpgp"],
                                             dg["qs"]))
        print("    cross-check produktowy (lambda^2=-nu/4): max Re lambda ="
              " %.6f vs JL %.6f (|d|=%.1e)"
              % (dg["mr_prod"], maxre[(N, 0.0)],
                 abs(dg["mr_prod"] - maxre[(N, 0.0)])))
        print("    4 najwieksze Re lambda: %s; mody symetrii wsrod"
              " Re>1e-3: %d; max Re bez modow symetrii = %.6f"
              % (["%+.4f" % v for v in dg["top4"]], dg["n_sym"],
                 dg["mr_excl"]))
    if not all_conv:
        vP = vL = "INCONCLUSIVE"
    else:
        vP = "FAIL" if any_fail_P else ("PASS" if all_stab_P else "MIXED")
        vL = "FAIL" if any_fail_L else ("PASS" if all_stab_L else "MIXED")
    verdicts_P[label] = vP
    verdicts_L[label] = vL
    mr_main = maxre[(800, 0.0)]
    mr_glob = max(maxre[(800, float(k))] for k in ks)
    argk = [float(k) for k in ks
            if maxre[(800, float(k))] == mr_glob][0]
    print("  WERDYKT tla %s: PRIMARY %s (literalny: %s);"
          % (label, vP, vL))
    print("    max Re lambda(k=0, N800) = %+.6f; max po k = %+.6f"
          " (argmax k=%.4f)" % (mr_main, mr_glob, argk))
    res_out[label] = dict(
        d=d, tol_primary=tolP, tol_literal=tolL, rows=rows,
        maxre_k0_N800=mr_main, maxre_global_N800=mr_glob, argmax_k=argk,
        verdict_primary=vP, verdict_literal=vL,
        diag_k0={str(N): diag_k0[N] for N in diag_k0},
        input_onto="u0 = g_d(x) (modul pola = g; LOCK sec.2)")

with open(OUT_JSON, "w") as f:
    json.dump(res_out, f, indent=1)
print("\nzapisano wyniki: %s" % OUT_JSON)

# ------------------------------------------------------------------ werdykt
print("\n" + "=" * 78)
print("TABELA max Re lambda (N=800) per tlo (kryteria LOCK doslownie):")
for label, r in res_out.items():
    print("  tlo %-14s d=%8.4f: max Re lambda(k=0) = %+.6f, po k ="
          " %+.6f; tolP=%.1e -> %s (lit.: %s)"
          % (label, r["d"], r["maxre_k0_N800"], r["maxre_global_N800"],
             r["tol_primary"], r["verdict_primary"], r["verdict_literal"]))
vs = list(verdicts_P.values())
if all(v == "PASS" for v in vs):
    Q = "Q-PASS"
elif all(v == "FAIL" for v in vs):
    Q = "Q-FAIL"
elif any(v == "INCONCLUSIVE" for v in vs) and not any(
        v == "FAIL" for v in vs):
    Q = "Q-INCONCLUSIVE"
elif any(v == "FAIL" for v in vs) and any(v == "PASS" for v in vs):
    Q = "MIXED"
else:
    Q = "Q-FAIL(czesciowo; patrz per tlo)" if any(
        v == "FAIL" for v in vs) else "Q-INCONCLUSIVE"
vsL = list(verdicts_L.values())
if all(v == "PASS" for v in vsL):
    QL = "Q-PASS"
elif all(v == "FAIL" for v in vsL):
    QL = "Q-FAIL"
elif any(v == "FAIL" for v in vsL) and any(v == "PASS" for v in vsL):
    QL = "MIXED"
else:
    QL = "Q-INCONCLUSIVE"
print("\nWERDYKT Phase 3 (odczyt PRIMARY tol, ruling FROZEN): %s" % Q)
print("  (odczyt literalny tol -- raportowany, degeneracja udokumentowana"
      " w Phase_method_decisions sec. 4: %s)" % QL)
print("  kontrole Phase 2: PASS (Gate A/C1/C2/C3) -- warunek 'kontrole"
      " czyste' spelniony.")
print("=" * 78)
