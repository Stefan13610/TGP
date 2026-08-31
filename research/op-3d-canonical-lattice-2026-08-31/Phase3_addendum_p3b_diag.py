#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
Addendum diagnostyczne P3b (d=4pi, N=48; maxerr=2.482e-01 FAIL):
rozstrzygniecie, czy FAIL jest bledem operatora (dyskretyzacji), czy
artefaktem ARPACK (gubione krotnosci mimo tol=0 -- rodzina bledu
z Phase_correction_note_eigsh.md poprzednika).

Metoda: proznia g=1 => operator wagi-K jest DOKLADNIE FD-Laplasjanem
(K(1)=1) + diag(-1); jego pelne spektrum jest znane analitycznie
(symbol FD): lam(m) = sum_a (2-2cos(q_a h))/h^2 - 1,
q_a = (theta_a + 2pi m_a)/d, theta_a in {0,pi} wg punktu k, m_a=0..N-1.
Porownanie: (a) eigsh(k=10) vs FD-exact (czy eigsh zwraca podzbior
z przesunietymi pozycjami = zgubione kopie); (b) FD-exact vs continuum
|k+G|^2-1 (czy sama dyskretyzacja miesci sie w progu 1e-2);
(c) eigsh z eskalacja k=32, ncv=240 vs FD-exact (czy eskalacja odzyskuje
10 najnizszych galezi).

Zero zmian kryteriow; czysta diagnostyka przed correction note.
"""
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

BETA, GAMMA = 1.0, 1.0
RNG_SEED = 20260831
K_fun = lambda g: g ** 4

t0 = time.time()


def stamp(m):
    print("[t=%6.1fs] %s" % (time.time() - t0, m), flush=True)


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


def v0_det(n):
    v = np.random.default_rng(RNG_SEED).standard_normal(n)
    return v / np.linalg.norm(v)


def fd_exact(N, d, theta, nb):
    """Pelny multizbior FD-symbolu (proznia, waga 1), nb najnizszych."""
    h = d / N
    lam1 = []
    for th in theta:
        m = np.arange(N)
        q = (th + 2 * np.pi * m) / d
        lam1.append((2 - 2 * np.cos(q * h)) / h ** 2)
    L = (lam1[0][:, None, None] + lam1[1][None, :, None]
         + lam1[2][None, None, :]).ravel() - GAMMA
    return np.sort(L)[:nb]


def cont_exact(kvec, d, nb):
    G = 2 * np.pi / d * np.arange(-8, 9)
    KX, KY, KZ = np.meshgrid(kvec[0] + G, kvec[1] + G, kvec[2] + G,
                             indexing='ij')
    return np.sort((KX ** 2 + KY ** 2 + KZ ** 2).ravel())[:nb] - GAMMA


d = 4 * np.pi
N = 48
h = d / N
g3 = np.ones((N, N, N))
Q3 = K_fun(g3) * (2 * BETA * g3 - 3 * GAMMA * g3 ** 2)
print("P3b diag: d=4pi, N=48, h=%.4f" % h)
for kname, kv, phs in kpoints(d):
    theta = [0.0 if p == 1 else np.pi for p in phs]
    ex_fd = fd_exact(N, d, theta, 12)
    ex_ct = cont_exact(kv, d, 12)
    A = build_op3d(g3, K_fun, Q3, h, phases=phs, weight3=K_fun(g3))
    n = A.shape[0]
    stamp("k=%s ..." % kname)
    v10 = np.sort(spla.eigsh(A, k=10, which='SA', tol=0, ncv=80,
                             maxiter=200000, v0=v0_det(n),
                             return_eigenvectors=False))
    err_solver = np.max(np.abs(v10 - ex_fd[:10]))
    err_fd_ct = np.max(np.abs(ex_fd[:10] - ex_ct[:10])
                       / np.maximum(np.abs(ex_ct[:10]), 1.0))
    err_p3b = np.max(np.abs(v10 - ex_ct[:10])
                     / np.maximum(np.abs(ex_ct[:10]), 1.0))
    print("  eigsh(k=10,ncv=80): %s" % ["%+.4f" % x for x in v10])
    print("  FD-exact (12 najn.): %s" % ["%+.4f" % x for x in ex_fd])
    print("  continuum (12 najn.): %s" % ["%+.4f" % x for x in ex_ct])
    print("  |eigsh - FD-exact|_max(10) = %.3e  (0 => operator OK, blad"
          " = solver)" % err_solver)
    print("  |FD-exact - continuum|/max(|ex|,1) (10) = %.3e  (czysty blad"
          " dyskretyzacji)" % err_fd_ct)
    print("  metryka P3b (eigsh vs continuum) = %.3e" % err_p3b, flush=True)
    if err_p3b > 1e-2:
        stamp("  eskalacja: eigsh(k=32, ncv=240)...")
        v32 = np.sort(spla.eigsh(A, k=32, which='SA', tol=0, ncv=240,
                                 maxiter=200000, v0=v0_det(n),
                                 return_eigenvectors=False))
        err_s2 = np.max(np.abs(v32[:10] - ex_fd[:10]))
        err_p3b2 = np.max(np.abs(v32[:10] - ex_ct[:10])
                          / np.maximum(np.abs(ex_ct[:10]), 1.0))
        print("  eigsh(k=32): 10 najn.: %s" % ["%+.4f" % x for x in v32[:10]])
        print("  |eigsh32 - FD-exact|_max(10) = %.3e; metryka P3b po"
              " eskalacji = %.3e" % (err_s2, err_p3b2), flush=True)
print("\nWniosek drukowany przez liczby powyzej: jesli |eigsh-FDexact|>0"
      " przy k=10, a =0 przy k=32, FAIL byl artefaktem ARPACK (zgubione"
      " kopie degeneracji), nie bledem operatora/dyskretyzacji.")
