#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-3d-lattice-bath-stability (Phase 3) -- RACHUNEK CENTRALNY: spektrum
Blocha wokol samouzgodnionych tel sieci sc 3D (sektor tachionowy, akcja
kanoniczna K=g^4, U=g^7/7-g^8/8, beta=gamma=1).

LOCK: Phase0_balance.md sec. 3 Phase 3; decyzje FROZEN:
Phase_method_decisions.md sec. 1, 2, 6; korekta 1 (eigsh tol=0 +
kontrola pokrycia translacji): Phase_correction_note_eigsh.md
(zapisana PRZED tym rachunkiem).

P3a: lambda_min(L_d(k)) w punktach Gamma, X, M, R (zalockowane);
     N in {32,48}; operator -div(K grad phi) + Q phi = w2 K phi,
     Q = K(2g-3g^2+2|grad g|^2/g^2); NIEUSUWALNE: identyfikacja 3 modow
     translacyjnych (overlap B-wazony z podprzestrzenia {d_i g} >= 0.9;
     lambda ~ O(h^2) w Gamma; ilorazy Rayleigha niezaleznie od ARPACK)
     PRZED interpretacja; w2_min = minimum PO ODJECIU modow
     translacyjnych (litera LOCKa); pelne minimum raportowane obok.
     Zbieznosc: |d w2_min(32->48)| <= 0.05 max(|w2_min|, 0.1).
     Cross-check u-formy (-lap + 4g-5g^2, waga 1) w {Gamma, R}, N=32.
P3b (NIEUSUWALNA): proznia g=1 w tej samej komorce, wszystkie zalockowane
     d, N in {32,48}, 4 punkty k, 10 galezi vs |k+G|^2-gamma;
     |d|/max(|ex|,1) <= 1e-2.
P3c (NIEUSUWALNA): sektor stabilny m^2=+gamma ze zrodlem gaussowskim 3D
     (d_s=4, sigma_s=0.5, q=1), relaksacja newton_krylov + Bloch;
     wymog: w2_min > 0 na obu siatkach.

INPUT: g0_mu=2.02117, d*1=3.0790, beta=gamma=1, eps (nieuzywane w
spektrach -- bez regularyzacji, LOCK sec. 2).
"""
import json
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import newton_krylov, NoConvergence

BETA = 1.0
GAMMA = 1.0
DSTAR1 = 3.0790                      # INPUT
BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-3d-lattice-bath-stability-2026-08-31/")
NPZ = BASE + "Phase2_backgrounds3d.npz"
OUT_JSON = BASE + "Phase3_results3d.json"
D_ALL = {"pi": np.pi, "2pi": 2 * np.pi, "dstar1": DSTAR1,
         "3pi": 3 * np.pi, "4pi": 4 * np.pi}
RNG_SEED = 20260831
K_fun = lambda g: g ** 4

t0_wall = time.time()


def stamp(msg):
    print("[t=%7.1fs] %s" % (time.time() - t0_wall, msg), flush=True)


def kpoints(d):
    q = np.pi / d
    return [("Gamma", (0.0, 0.0, 0.0), (1, 1, 1)),
            ("X", (q, 0.0, 0.0), (-1, 1, 1)),
            ("M", (q, q, 0.0), (-1, -1, 1)),
            ("R", (q, q, q), (-1, -1, -1))]


def build_op3d(g3, Phi, Q3, h, phases=(1, 1, 1), weight3=None):
    N = g3.shape[0]
    n = N ** 3
    idx3 = np.arange(n).reshape(N, N, N)
    diag = Q3.astype(float).ravel().copy()
    rows, cols, vals = [], [], []
    for ax in range(3):
        gn = np.roll(g3, -1, axis=ax)
        c = Phi(0.5 * (g3 + gn)) / h ** 2
        diag += (c + np.roll(c, 1, axis=ax)).ravel()
        ph = np.ones_like(g3)
        sl = [slice(None)] * 3
        sl[ax] = N - 1
        ph[tuple(sl)] = phases[ax]
        v = (-c * ph).ravel()
        r_ = idx3.ravel()
        c_ = np.roll(idx3, -1, axis=ax).ravel()
        rows += [r_, c_]
        cols += [c_, r_]
        vals += [v, v]
    rows.append(np.arange(n))
    cols.append(np.arange(n))
    vals.append(diag)
    A = sp.coo_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(n, n)).tocsr()
    if weight3 is not None:
        dinv = 1.0 / np.sqrt(weight3.ravel())
        D = sp.diags(dinv)
        A = D @ A @ D
    return A


_V0 = {}


def v0_det(n):
    if n not in _V0:
        v = np.random.default_rng(RNG_SEED).standard_normal(n)
        _V0[n] = v / np.linalg.norm(v)
    return _V0[n]


def lowest_eigs(A, nev=10, ncv=80, vecs=False):
    """eigsh 'SA'; tol=0 (korekta 1 -- pelna precyzja, krotnosci)."""
    n = A.shape[0]
    kw = dict(k=nev, which='SA', tol=0, ncv=min(ncv, n - 1),
              maxiter=200000, v0=v0_det(n))
    try:
        if vecs:
            vals, V = spla.eigsh(A, **kw)
            o = np.argsort(vals)
            return vals[o], V[:, o]
        return np.sort(spla.eigsh(A, return_eigenvectors=False, **kw)), None
    except spla.ArpackNoConvergence as e:
        return None, None


def d1_axes(g, h):
    return [(np.roll(g, -1, axis=ax) - np.roll(g, 1, axis=ax)) / (2 * h)
            for ax in range(3)]


def grad2_of(g, h):
    d1 = d1_axes(g, h)
    return d1[0] ** 2 + d1[1] ** 2 + d1[2] ** 2


def tach_Q(g, h):
    return K_fun(g) * (2 * BETA * g - 3 * GAMMA * g ** 2
                       + 2 * grad2_of(g, h) / g ** 2)


def stab_Q(g, h):
    return -K_fun(g) * (2 * BETA * g - 3 * GAMMA * g ** 2
                        - 2 * grad2_of(g, h) / g ** 2)


def exact_vacuum(kvec, d, m2, nb=10):
    G = 2 * np.pi / d * np.arange(-6, 7)
    KX, KY, KZ = np.meshgrid(kvec[0] + G, kvec[1] + G, kvec[2] + G,
                             indexing='ij')
    return np.sort((KX ** 2 + KY ** 2 + KZ ** 2).ravel())[:nb] + m2


def spectrum_bg(g, h, d, kname, kv, phs, nev=10, ncv=80):
    """Zwraca (vals, V_psi, sqw) dla tla g (sektor tachionowy, waga K)."""
    K3 = K_fun(g)
    A = build_op3d(g, K_fun, tach_Q(g, h), h, phases=phs, weight3=K3)
    vals, V = lowest_eigs(A, nev=nev, ncv=ncv, vecs=True)
    return vals, V, np.sqrt(K3.ravel()), A


def translation_analysis(g, h, vals, V, sqw, A):
    """Overlapy modow z podprzestrzenia translacyjna (B-wazone; FROZEN 0.9)
    + pokrycie c_i + ilorazy Rayleigha (korekta 1). Wektory w konwencji
    psi = B^(1/2) phi (eigsh na H); tau_i = B^(1/2) d_i g."""
    taus = []
    for ax in range(3):
        t = ((np.roll(g, -1, axis=ax) - np.roll(g, 1, axis=ax))
             / (2 * h)).ravel() * sqw
        taus.append(t)
    T = np.stack(taus, axis=1)                    # (n,3)
    # ortonormalizacja (Gram) w zwyklym ilorazie (psi-przestrzen = B-wazenie)
    Qt, _ = np.linalg.qr(T)
    ovl = []
    for j in range(V.shape[1]):
        psi = V[:, j]
        proj = Qt.T @ psi
        ovl.append(float(np.linalg.norm(proj) / np.linalg.norm(psi)))
    # pokrycie translacji przez zwrocone mody
    cov = []
    for i in range(3):
        t = T[:, i] / np.linalg.norm(T[:, i])
        cov.append(float(np.linalg.norm(V.T @ t)))
    # ilorazy Rayleigha (niezalezne od ARPACK)
    ray = []
    for i in range(3):
        t = T[:, i]
        ray.append(float((t @ (A @ t)) / (t @ t)))
    return ovl, cov, ray


# ===================================================================== RUN
print("=" * 78)
print("Phase 3 -- RACHUNEK CENTRALNY: Bloch 3D wokol tel sieci sc.")
print("LOCK sec. 3 Phase 3; FROZEN: method_decisions sec. 1/2/6;")
print("korekta 1 (eigsh tol=0): Phase_correction_note_eigsh.md.")
print("INPUT: g0_mu=2.02117, d*1=%.4f, beta=gamma=1." % DSTAR1)
print("=" * 78)

data = np.load(NPZ)
bgs = []
for key in sorted(data.files):
    if key.endswith("__d"):
        dn, tag, _ = key.split("__")
        bgs.append((dn, tag, float(data[key][0])))
print("tla z Phase 2: %s" % (", ".join("%s/%s" % (dn, tag)
                                       for dn, tag, _ in bgs)
                             if bgs else "BRAK"))

# --------------------------------------------------------- P3b (nieusuwalna)
stamp("[P3b] KONTROLA: proznia g=1 w komorce [0,d)^3 -- galezie |k+G|^2-1")
p3b_pass = True
p3b_tab = {}
for dn, d in D_ALL.items():
    for N in (32, 48):
        h = d / N
        g3 = np.ones((N, N, N))
        Q3 = tach_Q(g3, h)
        maxerr = 0.0
        for kname, kv, phs in kpoints(d):
            A = build_op3d(g3, K_fun, Q3, h, phases=phs, weight3=K_fun(g3))
            vals, _ = lowest_eigs(A, nev=10)
            if vals is None:
                maxerr = np.inf
                break
            ex = exact_vacuum(kv, d, -GAMMA, nb=10)
            maxerr = max(maxerr, float(np.max(
                np.abs(vals - ex) / np.maximum(np.abs(ex), 1.0))))
        ok = maxerr <= 1e-2
        p3b_tab["%s_N%d" % (dn, N)] = maxerr
        stamp("  d=%-6s N=%2d: maxerr(10 galezi, 4 pkt k) = %.3e -> %s"
              % (dn, N, maxerr, "PASS" if ok else "FAIL"))
        p3b_pass = p3b_pass and ok
print("  P3b: %s" % ("PASS" if p3b_pass else "FAIL"))

# --------------------------------------------------------- P3c (nieusuwalna)
stamp("[P3c] KONTROLA: sektor stabilny m^2=+1 ze zrodlem 3D "
      "(d_s=4, sigma=0.5, q=1)")
DS, SIG, QQ = 4.0, 0.5, 1.0


def rho_lattice3d(N, h, d, sigma):
    x = np.arange(N) * h
    c = d / 2
    rho = np.zeros((N, N, N))
    for mx in range(-2, 3):
        for my in range(-2, 3):
            for mz in range(-2, 3):
                DX = x - c - mx * d
                DY = x - c - my * d
                DZ = x - c - mz * d
                RR2 = (DX[:, None, None] ** 2 + DY[None, :, None] ** 2
                       + DZ[None, None, :] ** 2)
                rho += np.exp(-0.5 * RR2 / sigma ** 2)
    return rho / ((2 * np.pi) ** 1.5 * sigma ** 3)


def make_precond(N, h):
    n = N ** 3
    idx3 = np.arange(n).reshape(N, N, N)
    rows, cols, vals = [], [], []
    diag = np.full(n, -6.0 / h ** 2 - 1.0)
    for ax in range(3):
        r_ = idx3.ravel()
        c_ = np.roll(idx3, -1, axis=ax).ravel()
        v = np.full(n, 1.0 / h ** 2)
        rows += [r_, c_]
        cols += [c_, r_]
        vals += [v, v]
    rows.append(np.arange(n))
    cols.append(np.arange(n))
    vals.append(diag)
    A = sp.coo_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(n, n)).tocsc()
    ilu = spla.spilu(A, drop_tol=1e-5, fill_factor=15)
    return spla.LinearOperator((n, n), matvec=ilu.solve)


p3c_vals = {}
p3c_pass = True
for N in (32, 48):
    h = DS / N
    rho = rho_lattice3d(N, h, DS, SIG)
    M = make_precond(N, h)
    psi = np.ones(N ** 3)
    okrel = True
    for frac in (0.25, 0.5, 0.75, 1.0):
        def res_s(pv, _f=frac):
            p = pv.reshape(N, N, N)
            lap = np.zeros_like(p)
            grad2 = np.zeros_like(p)
            for ax in range(3):
                pp = np.roll(p, -1, axis=ax)
                pm = np.roll(p, 1, axis=ax)
                lap += (pp - 2 * p + pm) / h ** 2
                grad2 += ((pp - pm) / (2 * h)) ** 2
            return (lap + 2 * grad2 / p + BETA * p ** 2
                    - GAMMA * p ** 3 + _f * QQ * rho).ravel()
        try:
            psi = newton_krylov(res_s, psi, method='lgmres', inner_M=M,
                                f_tol=1e-9, maxiter=300, inner_maxiter=200)
        except NoConvergence as e:
            psi = np.asarray(e.args[0])
            okrel = False
    psi3 = psi.reshape(N, N, N)
    rmax = float(np.max(np.abs(res_s(psi, 1.0))))
    vmin, kargm = np.inf, "?"
    for kname, kv, phs in kpoints(DS):
        A = build_op3d(psi3, K_fun, stab_Q(psi3, h), h, phases=phs,
                       weight3=K_fun(psi3))
        vals, _ = lowest_eigs(A, nev=10)
        if vals is None:
            vmin = np.nan
            break
        if float(vals[0]) < vmin:
            vmin, kargm = float(vals[0]), kname
    p3c_vals[N] = vmin
    ok = okrel and np.isfinite(vmin) and vmin > 0
    p3c_pass = p3c_pass and ok
    stamp("  N=%2d: psi in [%.5f, %.5f], ||R||=%.1e; w2_min = %+.5f "
          "(argmin %s) -> %s" % (N, float(np.min(psi3)),
                                 float(np.max(psi3)), rmax, vmin, kargm,
                                 "PASS" if ok else "FAIL"))
print("  P3c: %s (wymog w2_min > 0 na obu siatkach)"
      % ("PASS" if p3c_pass else "FAIL"))

# ------------------------------------------------------------- P3a centralny
stamp("[P3a] spektra tel Phase 2 (sektor tachionowy)")
res_out = {"P3b": p3b_tab, "P3b_pass": bool(p3b_pass),
           "P3c": {str(k): (None if not np.isfinite(v) else v)
                   for k, v in p3c_vals.items()},
           "P3c_pass": bool(p3c_pass), "bgs": {}}
for dn, tag, d in bgs:
    label = "%s/%s" % (dn, tag)
    print("\n  --- tlo %s (d=%.6f) ---" % (label, d), flush=True)
    per_N = {}
    for N in (32, 48):
        g = data["%s__%s__N%d" % (dn, tag, N)]
        h = d / N
        kmins_clean, kmins_full = {}, {}
        trans_info = None
        incomplete = False
        for kname, kv, phs in kpoints(d):
            t1 = time.time()
            vals, V, sqw, A = spectrum_bg(g, h, d, kname, kv, phs)
            if vals is None:
                stamp("    N=%d k=%s: ARPACK NIEZBIEZNY -> INCOMPLETE"
                      % (N, kname))
                incomplete = True
                continue
            ovl, cov, ray = translation_analysis(g, h, vals, V, sqw, A)
            if kname == "Gamma" and min(cov) < 0.99:
                stamp("    N=%d Gamma: pokrycie translacji %s < 0.99 -> "
                      "eskalacja k=16, ncv=160"
                      % (N, ["%.3f" % c for c in cov]))
                K3 = K_fun(g)
                A2 = build_op3d(g, K_fun, tach_Q(g, h), h, phases=phs,
                                weight3=K3)
                vals, V = lowest_eigs(A2, nev=16, ncv=160, vecs=True)
                if vals is None:
                    incomplete = True
                    continue
                ovl, cov, ray = translation_analysis(g, h, vals, V, sqw, A2)
            tr_idx = [j for j in range(len(vals)) if ovl[j] >= 0.9]
            clean = [float(vals[j]) for j in range(len(vals))
                     if j not in tr_idx]
            kmins_full[kname] = float(vals[0])
            kmins_clean[kname] = clean[0] if clean else np.nan
            if kname == "Gamma":
                trans_info = dict(
                    n_trans=len(tr_idx),
                    lam_trans=[float(vals[j]) for j in tr_idx],
                    ovl=[round(o, 3) for o in ovl],
                    coverage=[round(c, 4) for c in cov],
                    rayleigh=[float(x) for x in ray])
            stamp("    N=%2d k=%-5s: 10 najn.: %s | trans idx %s [%.0fs]"
                  % (N, kname, ["%+.4f" % x for x in vals[:6]], tr_idx,
                     time.time() - t1))
        if incomplete or not kmins_clean:
            per_N[N] = None
            continue
        argk = min(kmins_clean, key=kmins_clean.get)
        per_N[N] = dict(w2min=float(kmins_clean[argk]), argk=argk,
                        w2min_full=float(min(kmins_full.values())),
                        argk_full=min(kmins_full, key=kmins_full.get),
                        kmins=dict(kmins_clean), kmins_full=dict(kmins_full),
                        trans=trans_info)
        print("    N=%2d: w2_min(clean po k) = %+.6f (argmin %s); "
              "pelne min = %+.6f (%s); Gamma: %d modow transl., "
              "lam_trans=%s, Rayleigh=%s"
              % (N, per_N[N]["w2min"], argk, per_N[N]["w2min_full"],
                 per_N[N]["argk_full"], trans_info["n_trans"],
                 ["%+.2e" % x for x in trans_info["lam_trans"]],
                 ["%+.2e" % x for x in trans_info["rayleigh"]]),
              flush=True)
    # zbieznosc siatkowa (LOCK)
    entry = dict(d=d)
    if per_N.get(32) and per_N.get(48):
        v32, v48 = per_N[32]["w2min"], per_N[48]["w2min"]
        dv = abs(v32 - v48)
        tol = 0.05 * max(abs(v48), 0.1)
        conv = dv <= tol
        print("    GLOWNA (N=48): w2_min = %+.6f; |delta(32->48)| = %.2e "
              "(prog %.2e) -> %s" % (v48, dv, tol,
                                     "ZBIEZNE" if conv else "NIEZBIEZNE"))
        entry.update(w2min=v48, w2min_N32=v32, conv=bool(conv), dv=dv,
                     tol=tol, argk=per_N[48]["argk"],
                     w2min_full=per_N[48]["w2min_full"],
                     kmins=per_N[48]["kmins"],
                     kmins_full=per_N[48]["kmins_full"],
                     trans=per_N[48]["trans"])
    else:
        entry.update(w2min=None, conv=False, incomplete=True)
        print("    INCOMPLETE (ARPACK/dane) -- raportowane wprost")
    # cross-check u-formy (waga 1) w {Gamma, R}, N=32
    g32 = data["%s__%s__N32" % (dn, tag)]
    h32 = d / 32
    for kname, kv, phs in [kp for kp in kpoints(d)
                           if kp[0] in ("Gamma", "R")]:
        Au = build_op3d(g32, lambda x: np.ones_like(x),
                        4 * BETA * g32 - 5 * GAMMA * g32 ** 2, h32,
                        phases=phs)
        vu, _ = lowest_eigs(Au, nev=1)
        Aw = build_op3d(g32, K_fun, tach_Q(g32, h32), h32, phases=phs,
                        weight3=K_fun(g32))
        vw, _ = lowest_eigs(Aw, nev=1)
        if vu is not None and vw is not None:
            print("    cross-check u-forma k=%-5s (N=32): waga-K %+.6f vs "
                  "u-forma %+.6f (|delta|=%.2e)"
                  % (kname, float(vw[0]), float(vu[0]),
                     abs(float(vw[0]) - float(vu[0]))), flush=True)
    res_out["bgs"][label] = entry

with open(OUT_JSON, "w") as f:
    json.dump(res_out, f, indent=1)
print("\nzapisano wyniki: %s" % OUT_JSON)

# ------------------------------------------------------------------ tabela
print("\n" + "=" * 78)
print("TABELA w2_min(d) (po odjeciu translacji; zbieznosc LOCK doslownie):")
print("  P3b: %s   P3c: %s" % ("PASS" if p3b_pass else "FAIL",
                               "PASS" if p3c_pass else "FAIL"))
for label, rr in res_out["bgs"].items():
    if rr.get("w2min") is None:
        print("  tlo %-14s d=%8.4f: INCOMPLETE" % (label, rr["d"]))
        continue
    print("  tlo %-14s d=%8.4f: w2_min = %+.6f (%s; argmin %s; "
          "pelne min %+.6f)"
          % (label, rr["d"], rr["w2min"],
             "ZBIEZNE" if rr["conv"] else "NIEZBIEZNE", rr["argk"],
             rr["w2min_full"]))
vals_ok = [rr for rr in res_out["bgs"].values() if rr.get("w2min") is not None]
any_pos = any(rr["w2min"] > 0 and rr["conv"] for rr in vals_ok)
all_neg = (len(vals_ok) > 0
           and all(rr["w2min"] < 0 and rr["conv"] for rr in vals_ok))
print("\n  istnieje d z w2_min > 0 (zbiezne): %s" % any_pos)
print("  w2_min < 0 dla wszystkich istniejacych tel (zbieznie): %s" % all_neg)
print("  (werdykt Q wymaga Phase 4; ruling kwantyfikatora PRZED Phase 4)")
print("  deskryptywnie vs 1D: bloch-chain w2_min = -1.222...-1.230")
print("=" * 78)
