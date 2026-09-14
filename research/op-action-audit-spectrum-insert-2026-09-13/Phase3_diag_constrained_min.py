#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-action-audit-spectrum-insert (Phase 3, DIAGNOSTYKA DESKRYPTYWNA --
NIE werdyktotworcza, zero zmian kryteriow/progow/schematu FROZEN).

Cel: kontrola wiarygodnosci tabeli DeltaE_insert. Zamrozony schemat
"krok semi-implicit -> projekcja pinu" ma punkt staly spelniajacy
rhs_i=0 dla i>=2, ale rhs_1 = -d*A_t*a_0/m_1 != 0 (tozsamosc
zweryfikowana na biegach: pnorm_end ~ 76.6 przy d=4.56e-2, h=0.0125)
-- stad status TMAX zamiast STATIONARY. Tu liczymy DOKLADNE dyskretne
minimum warunkowe (Newton na wezlach swobodnych, psi_0 = A ustalone,
||dH||_inf <= 1e-13) startujac ze stanu koncowego flow i porownujemy
DeltaE_exact z DeltaE_flow. Potwierdzenie/sfalsyfikowanie skalowania
DeltaE ~ h (zerowa pojemnosc punktu w 3D).

Przypadki: A in {0.70, 1.25} x h in {0.025, 0.0125}, R=60
(przypadki STATIONARY-plateau; A=1.30 wylaczony -- BREAKDOWN).
"""
import json
import numpy as np
from scipy.linalg import solve_banded

BASE = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/"
        "op-action-audit-spectrum-insert-2026-09-13/")
RESDIR = BASE + "Phase3_results/"
OUT = BASE + "Phase3_diag_output.txt"

import importlib.util
spec = importlib.util.spec_from_file_location(
    "p3", BASE + "Phase3_insert_cost.py")
p3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p3)

lines = []


def em(s=""):
    print(s, flush=True)
    lines.append(s)


def dH_full(flow, g):
    """Gradient dyskretnej energii (bez podzialu przez masy)."""
    h = flow.h
    gm = 0.5 * (g[:-1] + g[1:])
    dg = np.diff(g) / h
    dH = flow.ctrap * h * flow.r ** 2 * p3.Ueffp(g)
    t_flux = flow.rm2 * p3.Keff(gm) * dg
    t_quad = 0.25 * h * flow.rm2 * p3.Keffp(gm) * dg ** 2
    dH[:-1] += -t_flux + t_quad
    dH[1:] += t_flux + t_quad
    return dH


def newton_constrained(flow, g0, A, tol=1e-13, itmax=200):
    """Newton na wezlach swobodnych (i>=1), psi_0=A; Jacobian
    trojdiagonalny przez FD z 3-kolorowaniem."""
    g = g0.copy()
    g[0] = A
    N1 = flow.N + 1
    for it in range(itmax):
        F = dH_full(flow, g)[1:]
        nrm = float(np.max(np.abs(F)))
        if nrm <= tol:
            return g, nrm, it
        eps = 1e-7
        Jl = np.zeros(N1 - 1)
        Jd = np.zeros(N1 - 1)
        Ju = np.zeros(N1 - 1)
        base = dH_full(flow, g)
        for c in range(3):
            gp = g.copy()
            idx = np.arange(1, N1)[(np.arange(1, N1) - 1) % 3 == c]
            gp[idx] += eps
            dcol = (dH_full(flow, gp) - base) / eps
            for i in idx:
                Jd[i - 1] = dcol[i]
                if i - 1 >= 1:
                    Ju[i - 1 - 1 + 1] = dcol[i - 1]   # element (i-1, i)
                if i + 1 <= N1 - 1:
                    Jl[i + 1 - 1 - 1] = dcol[i + 1]   # element (i+1, i)
        # solve_banded: ab[0]=naddiag (Ju przesuniete), ab[1]=diag,
        # ab[2]=poddiag
        ab = np.zeros((3, N1 - 1))
        ab[1] = Jd
        ab[0, 1:] = [Ju[j] for j in range(1, N1 - 1)]
        ab[2, :-1] = [Jl[j] for j in range(0, N1 - 2)]
        try:
            dg = solve_banded((1, 1), ab, -F)
        except Exception as e:
            return g, nrm, -1
        step = 1.0
        while step > 1e-6:
            gn = g.copy()
            gn[1:] += step * dg
            if np.all(np.isfinite(gn)) and np.max(
                    np.abs(dH_full(flow, gn)[1:])) < nrm:
                g = gn
                break
            step *= 0.5
        else:
            return g, nrm, -2
    return g, float(np.max(np.abs(dH_full(flow, g)[1:]))), itmax


em("=" * 78)
em("DIAGNOSTYKA DESKRYPTYWNA: dokladne dyskretne minimum warunkowe")
em("(Newton, ||dH_free||_inf <= 1e-13) vs punkt staly schematu flow")
em("(NIE werdyktotworcze; werdykt Q-D2 z Phase3_output.txt bez zmian)")
em("=" * 78)
rows = []
for lab, A in (("A070", 0.70), ("A125", 1.25)):
    for h in (0.025, 0.0125):
        jid = "%s_R60_h%g" % (lab, h)
        g_flow = np.load(RESDIR + jid + ".npz")["psi"]
        flow = p3.FlowRadialPin(60.0, h)
        with open(RESDIR + jid + ".json") as f:
            dE_flow = json.load(f)["dE_rel"]
        g_ex, nrm, it = newton_constrained(flow, g_flow, A)
        dE_ex = flow.energy_rel(g_ex)
        dev = float(np.max(np.abs(g_ex - g_flow)))
        em("%s: dE_flow=%+.8e dE_exact=%+.8e |roznica|=%.2e (%.3f%%)"
           % (jid, dE_flow, dE_ex, abs(dE_ex - dE_flow),
              100 * abs(dE_ex - dE_flow) / abs(dE_ex)))
        em("   Newton: ||dH||=%.2e po %d iter; max|psi_ex-psi_flow|=%.2e"
           % (nrm, it, dev))
        rows.append((lab, A, h, dE_ex))
em()
em("Skalowanie h (dokladne minima): ")
for lab, A in (("A070", 0.70), ("A125", 1.25)):
    d = {h: dE for l, a, h, dE in rows if l == lab}
    em("  %s: dE_exact(h=0.025)=%+.6e dE_exact(h=0.0125)=%+.6e "
       "stosunek=%.4f (0.5 = czyste skalowanie ~h, zerowa pojemnosc"
       " punktu w 3D)" % (lab, d[0.025], d[0.0125],
                          d[0.0125] / d[0.025]))
with open(OUT, "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("zapisano:", OUT)
