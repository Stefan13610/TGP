#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
op-lattice-bath-runaway (Phase A) -- AUDYT MASZYNERII 2 (BRAMKA CYKLU).
LOCK: Phase0_balance.md sec. 3 Phase A (A1-A6). Kryteria NIEZMIENNE.

NIEZALEZNA RE-IMPLEMENTACJA ODE (LOCK A1): rdzenne skrypty v47b uzywaja
scipy.solve_ivp (adaptacyjny RK45/DOP853). Tutaj: WLASNY stalokrokowy
klasyczny RK4 (bez scipy.integrate), wektoryzowany po warunkach
poczatkowych; prog kolapsu przez batch-refinement przedzialu (nie
bisekcja pojedynczych trajektorii). Zmienna calkowania dla Form A:
f = g^(alpha+1) (kolaps f->0 jest regularny; brak osobliwosci czlonu
(alpha/g)g'^2). ODE z biegnacym alpha_eff (A3) calkowane w g.

TESTY (wszystkie maja osiagalny FAIL -- LOCK sec. 4 pkt 2):
  A1: g0_crit(2,3) = 8/5 z |Delta| <= 1e-6; formula 2(a+2)/[2(a+2)-d]
      na >=3 parach (a,d) != (2,3). FAIL: |Delta| > 1e-3 lub brak
      generalizacji.
  A2: omega_tail = 1.000 +- 0.01 dla >=3 alpha; oscylacje w oknie
      |g-1| <= 0.2. Kontrola negatywna: syntetyczny zanik wykladniczy
      musi NIE przejsc procedury (pokazuje osiagalny FAIL).
  A3: ODE z alpha_eff(g)=2/(1+eta_K(g-1)^2), eta_K=181/15, g0_e=0.90548:
      delta_e=-81.4 +-2 deg, delta_mu=+38.6 +-2 deg, Delta(e->mu)=120 +-1;
      A_g ~ |g0-1| (1 +- 10%) przy prozni (plain alpha=2);
      okno fitu konczy sie <= 0.75*R; kontrola R w {R0, 1.5 R0}:
      parametry stabilne (FAIL: R-zaleznosc > 1%).
  A5: lokalizacja rownanie-po-rownaniu roznicy miedzy ODE maszynerii 2
      (dodatekH/O-L5/v47b) a "F-A kanoniczna" z AUDYT_KRYTYCZNY Dod. A
      + weryfikacja numeryczna obu zachowan.
  A6: blok uczciwosci p134e-g (powtorzony jawnie).

A4 (deskryptywny, bez PASS/FAIL): osobno -- PhaseA_A4_runner.py
przechwytuje faktyczne outputy 8 skryptow; tabela w PhaseA_report.md.
"""
import numpy as np

PHI = (1.0 + 5.0 ** 0.5) / 2.0
ETA_K = 181.0 / 15.0
G0_E_RUN = 0.90548              # dodatekH p130-p131 (running alpha_eff)
G0_MU_RUN = PHI * G0_E_RUN
G0_TAU_RUN = 4.0                # K-samodzielnosc (INPUT-adjacent, flagowane)

print("=" * 78)
print("Phase A: AUDYT MASZYNERII 2 (bramka cyklu op-lattice-bath-runaway)")
print("Niezalezna re-implementacja: wlasny stalokrokowy RK4 (bez")
print("scipy.integrate); rdzen v47b uzywa adaptacyjnego solve_ivp.")
print("=" * 78)

# =====================================================================
# Wspolny integrator: klasyczny RK4, staly krok, wektoryzowany po IC
# =====================================================================

def rk4_f_batch(alpha, d, f0_vec, r_max, h, sample_stride=0,
                f_floor=1e-9):
    """Form A w zmiennej kanonicznej f = g^(alpha+1):
       f'' + (d-1)/r f' = S(f),  S(f)=(a+1) f^p (1-f^q),
       p=(a+2)/(a+1), q=1/(a+1).
    Zwraca (collapsed, f_min, samples) -- samples=(r_s, f_s) gdy
    sample_stride>0 (tylko dla batcha dlugosci 1..3)."""
    a1 = alpha + 1.0
    p = (alpha + 2.0) / a1
    q = 1.0 / a1
    dm1 = float(d - 1)

    f = np.array(f0_vec, dtype=float).copy()
    n = f.size
    r0 = 1e-4
    S0 = a1 * f ** p * (1.0 - f ** q)
    fp = S0 * r0 / d
    f = f + 0.5 * S0 * r0 * r0 / d
    alive = np.ones(n, dtype=bool)
    f_min = f.copy()

    def rhs(r, fv, fpv):
        fc = np.clip(fv, 0.0, None)
        S = a1 * fc ** p * (1.0 - fc ** q)
        return fpv, S - dm1 * fpv / r

    nstep = int(np.ceil((r_max - r0) / h))
    rs, fs = [], []
    r = r0
    for k in range(nstep):
        if sample_stride and (k % sample_stride == 0):
            rs.append(r)
            fs.append(f.copy())
        k1f, k1p = rhs(r, f, fp)
        k2f, k2p = rhs(r + 0.5 * h, f + 0.5 * h * k1f, fp + 0.5 * h * k1p)
        k3f, k3p = rhs(r + 0.5 * h, f + 0.5 * h * k2f, fp + 0.5 * h * k2p)
        k4f, k4p = rhs(r + h, f + h * k3f, fp + h * k3p)
        df = h / 6.0 * (k1f + 2 * k2f + 2 * k3f + k4f)
        dp = h / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
        f = np.where(alive, f + df, f)
        fp = np.where(alive, fp + dp, fp)
        r += h
        f_min = np.minimum(f_min, f)
        newly_dead = alive & (f <= f_floor)
        if newly_dead.any():
            alive = alive & ~newly_dead
            if not alive.any() and not sample_stride:
                break
    collapsed = ~alive
    if sample_stride:
        return collapsed, f_min, (np.array(rs), np.array(fs))
    return collapsed, f_min, None


def find_gcrit(alpha, d, rounds=6, batch=41, r_max=400.0, h=0.005):
    """Prog kolapsu przez batch-refinement: w kazdej rundzie batch
    rownoodleglych g0, granica przezycie/kolaps zawezana do sasiedniej
    pary. Zwraca (gc, szerokosc_przedzialu, anomalia)."""
    a1 = alpha + 1.0
    gc_form = 2.0 * (alpha + 2.0) / (2.0 * (alpha + 2.0) - d)
    g_lo, g_hi = gc_form - 0.15, gc_form + 0.15
    anomaly = False
    for rnd in range(rounds):
        g0s = np.linspace(g_lo, g_hi, batch)
        collapsed, _, _ = rk4_f_batch(alpha, d, g0s ** a1, r_max, h)
        if collapsed[0] or not collapsed[-1]:
            # zly bracket -- poszerz (osiagalne tylko w 1. rundzie)
            g_lo -= 0.3
            g_hi += 0.3
            anomaly = True
            continue
        # pierwsza skolapsowana
        k = int(np.argmax(collapsed))
        if not np.all(collapsed[k:]) or np.any(collapsed[:k]):
            anomaly = True   # granica niemonotoniczna -- raportowac
        g_lo, g_hi = g0s[k - 1], g0s[k]
    return 0.5 * (g_lo + g_hi), (g_hi - g_lo), anomaly


# =====================================================================
# A1: prog kolapsu
# =====================================================================
print("\n" + "-" * 78)
print("[A1] prog kolapsu: g0_crit = 2(a+2)/[2(a+2)-d]")
print("-" * 78)

A1_CASES = [(2.0, 3), (1.0, 3), (3.0, 3), (2.0, 4), (1.5, 2)]
a1_pass = True
a1_rows = []
for alpha, d in A1_CASES:
    gc_form = 2.0 * (alpha + 2.0) / (2.0 * (alpha + 2.0) - d)
    gc_num, width, anom = find_gcrit(alpha, d)
    err = abs(gc_num - gc_form)
    tol = 1e-6
    ok = err <= tol and not anom
    if not ok:
        a1_pass = False
    a1_rows.append((alpha, d, gc_form, gc_num, err, width, anom))
    print("  alpha=%-4.1f d=%d  formula=%.10f  num=%.10f  |D|=%.2e"
          "  szer=%.1e  %s%s"
          % (alpha, d, gc_form, gc_num, err, width,
             "OK" if err <= tol else ("FAIL" if err > 1e-3 else "MARGINAL"),
             "  [ANOMALIA granicy]" if anom else ""))
print("  gate: (2,3) |D| <= 1e-6 oraz >=3 par (a,d)!=(2,3) zgodnych")
print("  A1: %s" % ("PASS" if a1_pass else "FAIL"))

# =====================================================================
# A2: ogon oscylacyjny, omega = 1 uniwersalne
# =====================================================================
print("\n" + "-" * 78)
print("[A2] ogon: omega_tail = 1.000 +- 0.01 dla alpha in {1,2,3} (d=3)")
print("-" * 78)


def tail_g_series(alpha, g0, r_max=200.0, h=0.004, stride=5):
    a1 = alpha + 1.0
    _, _, (rs, fs) = rk4_f_batch(alpha, 3, np.array([g0 ** a1]),
                                 r_max, h, sample_stride=stride)
    g = fs[:, 0] ** (1.0 / a1)
    return rs, g


def fit_omega(r, y, w_lo, w_hi):
    """Fit y ~ B cos(w r) + C sin(w r): siatka w + zawezenie."""
    m = (r >= w_lo) & (r <= w_hi)
    rr, yy = r[m], y[m]

    def resid(w):
        M = np.column_stack([np.cos(w * rr), np.sin(w * rr)])
        coef, res, _, _ = np.linalg.lstsq(M, yy, rcond=None)
        pred = M @ coef
        return float(np.sum((yy - pred) ** 2)), coef

    ws = np.arange(0.85, 1.1500001, 2e-4)
    rs_ = np.array([resid(w)[0] for w in ws])
    i = int(np.argmin(rs_))
    # zawezenie parabola
    if 0 < i < len(ws) - 1:
        y1, y2, y3 = rs_[i - 1], rs_[i], rs_[i + 1]
        denom = (y1 - 2 * y2 + y3)
        w_best = ws[i] + (0.5 * (y1 - y3) / denom) * 2e-4 if denom != 0 \
            else ws[i]
    else:
        w_best = ws[i]
    ss_res, coef = resid(w_best)
    ss_tot = float(np.sum((yy - yy.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return w_best, coef, r2


a2_pass = True
W_LO, W_HI = 40.0, 150.0
for alpha in (1.0, 2.0, 3.0):
    r, g = tail_g_series(alpha, 0.9)
    y = (g - 1.0) * r
    m = (r >= W_LO) & (r <= W_HI)
    ncross = int(np.sum(np.diff(np.sign(g[m] - 1.0)) != 0))
    in_win = float(np.max(np.abs(g[m] - 1.0)))
    w_fit, coef, r2 = fit_omega(r, y, W_LO, W_HI)
    ok = (abs(w_fit - 1.0) <= 0.01) and (ncross >= 10)
    if not ok:
        a2_pass = False
    print("  alpha=%.1f g0=0.9: omega_fit=%.5f  |w-1|=%.5f  "
          "zera(g-1)=%d  max|g-1| w oknie=%.3f (<=0.2: %s)  R2=%.6f  %s"
          % (alpha, w_fit, abs(w_fit - 1.0), ncross, in_win,
             "tak" if in_win <= 0.2 else "NIE", r2,
             "OK" if ok else "FAIL"))

# kontrola negatywna procedury: czysty zanik wykladniczy
r_syn = np.arange(0.02, 200.0, 0.02)
y_syn = 0.1 * np.exp(-r_syn)  # (g-1)*r dla ogona Yukawy ~ e^-r/r
m = (r_syn >= W_LO) & (r_syn <= W_HI)
ncross_syn = int(np.sum(np.diff(np.sign(y_syn[m])) != 0))
w_syn, _, r2_syn = fit_omega(r_syn, y_syn, W_LO, W_HI)
ctrl_fails = (ncross_syn < 10) or (r2_syn < 0.5)
print("  KONTROLA NEGATYWNA (syntetyczny e^-r, bez oscylacji): "
      "zera=%d, R2(fit osc)=%.3f -> procedura %s"
      % (ncross_syn, r2_syn,
         "poprawnie ODRZUCA (osiagalny FAIL potwierdzony)"
         if ctrl_fails else "NIE odrzuca -- BLAD PROCEDURY"))
if not ctrl_fails:
    a2_pass = False
print("  A2: %s" % ("PASS" if a2_pass else "FAIL"))

# =====================================================================
# A3: fazy i amplitudy (ODE z biegnacym alpha_eff -- dodatekH p130-131)
# =====================================================================
print("\n" + "-" * 78)
print("[A3] fazy ogona (alpha_eff(g)=2/(1+eta_K(g-1)^2), eta_K=181/15)")
print("    targets (LOCK): delta_e=-81.4+-2, delta_mu=+38.6+-2, "
      "Delta=120+-1; A_g~|g0-1| (1+-10%)")
print("-" * 78)


def rk4_g_running(g0_vec, r_max, h=0.001, stride=20, eta=ETA_K,
                  g_floor=1e-6):
    """ODE w g: g'' + 2 g'/r + (alpha_eff(g)/g) g'^2 = g^2 (1-g),
    alpha_eff(g) = 2/(1+eta (g-1)^2). Podstawienie bezposrednie
    (konwencja p129 'Rozwiazanie ODE z alpha_eff(g)')."""
    g = np.array(g0_vec, dtype=float).copy()
    n = g.size
    r0 = 1e-4
    src0 = g * g * (1.0 - g)
    gp = src0 * r0 / 3.0
    g = g + 0.5 * src0 * r0 * r0 / 3.0
    alive = np.ones(n, dtype=bool)

    def rhs(r, gv, gpv):
        gc = np.clip(gv, g_floor, None)
        aeff = 2.0 / (1.0 + eta * (gc - 1.0) ** 2)
        return gpv, (gc * gc * (1.0 - gc) - (aeff / gc) * gpv * gpv
                     - 2.0 * gpv / r)

    nstep = int(np.ceil((r_max - r0) / h))
    rs, gs = [], []
    r = r0
    for k in range(nstep):
        if k % stride == 0:
            rs.append(r)
            gs.append(g.copy())
        k1g, k1p = rhs(r, g, gp)
        k2g, k2p = rhs(r + 0.5 * h, g + 0.5 * h * k1g, gp + 0.5 * h * k1p)
        k3g, k3p = rhs(r + 0.5 * h, g + 0.5 * h * k2g, gp + 0.5 * h * k2p)
        k4g, k4p = rhs(r + h, g + h * k3g, gp + h * k3p)
        dg = h / 6.0 * (k1g + 2 * k2g + 2 * k3g + k4g)
        dp = h / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
        g = np.where(alive, g + dg, g)
        gp = np.where(alive, gp + dp, gp)
        r += h
        alive = alive & (g > g_floor) & np.isfinite(g)
    return np.array(rs), np.array(gs), ~alive


def fit_phase(r, y, w_lo, w_hi):
    """(g-1) r ~ B cos r + C sin r; A=sqrt(B^2+C^2);
    delta = atan2(C, B) [konwencja dodatekH lin. 1059]."""
    m = (r >= w_lo) & (r <= w_hi)
    M = np.column_stack([np.cos(r[m]), np.sin(r[m])])
    (B, C), _, _, _ = np.linalg.lstsq(M, y[m], rcond=None)
    A = float(np.hypot(B, C))
    d1 = float(np.degrees(np.arctan2(C, B)))     # atan2(C,B) -- primary
    d2 = float(np.degrees(np.arctan2(-C, B)))    # konwencja a3d
    return A, d1, d2, float(B), float(C)


R0_A3, R1_A3 = 200.0, 300.0
species = [("e", G0_E_RUN), ("mu", G0_MU_RUN), ("tau", G0_TAU_RUN)]
g0s = np.array([s[1] for s in species])
res_R = {}
for R in (R0_A3, R1_A3):
    rs, gs, dead = rk4_g_running(g0s, R)
    res_R[R] = (rs, gs, dead)

fits = {}
print("  okno fitu: [%.0f, %.0f]; 0.75*R0 = %.0f (koniec okna <= 0.75 R)"
      % (W_LO, W_HI, 0.75 * R0_A3))
for i, (name, g0) in enumerate(species):
    rs, gs, dead = res_R[R0_A3]
    y = (gs[:, i] - 1.0) * rs
    A, d1, d2, B, C = fit_phase(rs, y, W_LO, W_HI)
    fits[name] = (A, d1, d2)
    print("  %-3s g0=%.5f%s: A=%.6f  delta[atan2(C,B)]=%+8.2f deg  "
          "delta[atan2(-C,B)]=%+8.2f deg%s"
          % (name, g0, " (INPUT: K-samodz.)" if name == "tau" else "",
             A, d1, d2, "  [KOLAPS!]" if dead[i] else ""))

# ktora konwencja odtwarza targety? (obie sprawdzone, wybor raportowany)
tgt = {"e": -81.4, "mu": +38.6}
best_conv, best_err = None, 1e9
for conv, idx in (("atan2(C,B)", 1), ("atan2(-C,B)", 2)):
    err = max(abs(fits["e"][idx] - tgt["e"]),
              abs(fits["mu"][idx] - tgt["mu"]))
    print("  konwencja %-11s: |d_e-target|=%.2f, |d_mu-target|=%.2f"
          % (conv, abs(fits["e"][idx] - tgt["e"]),
             abs(fits["mu"][idx] - tgt["mu"])))
    if err < best_err:
        best_err, best_conv = err, conv
idx = 1 if best_conv == "atan2(C,B)" else 2
d_e, d_mu, d_tau = fits["e"][idx], fits["mu"][idx], fits["tau"][idx]
Delta_emu = d_mu - d_e
ph_ok = (abs(d_e - tgt["e"]) <= 2.0 and abs(d_mu - tgt["mu"]) <= 2.0
         and abs(Delta_emu - 120.0) <= 1.0)
print("  wybrana konwencja: %s (dodatekH podaje atan2(C_tail,B_tail))"
      % best_conv)
print("  delta_e=%+.2f (tgt -81.4+-2)  delta_mu=%+.2f (tgt +38.6+-2)  "
      "Delta(e->mu)=%.2f (tgt 120+-1)  -> %s"
      % (d_e, d_mu, Delta_emu, "OK" if ph_ok else "FAIL"))
print("  delta_tau(g0=4)=%+.2f deg (rdzen: -27.3; informacyjnie, "
      "poza gate A3)" % d_tau)

# amplituda przy prozni: plain alpha=2 (Form A, stala alpha)
print("  amplituda przy prozni (plain Form A, alpha=2):")
amp_ok = True
for g0 in (0.95, 0.98, 1.02, 1.05):
    r, g = tail_g_series(2.0, g0)
    A, _, _, _, _ = fit_phase(r, (g - 1.0) * r, W_LO, W_HI)
    ratio = A / abs(g0 - 1.0)
    ok = 0.90 <= ratio <= 1.10
    if not ok:
        amp_ok = False
    print("    g0=%.2f: A=%.6f  A/|g0-1|=%.4f  %s"
          % (g0, A, ratio, "OK" if ok else "FAIL"))

# kontrola R: R0 vs 1.5 R0 (to samo okno fitu)
print("  kontrola R (R0=%.0f vs 1.5R0=%.0f, okno bez zmian):" %
      (R0_A3, R1_A3))
r_ok = True
for i, (name, g0) in enumerate(species[:2]):   # gate: e, mu
    rs1, gs1, _ = res_R[R0_A3]
    rs2, gs2, _ = res_R[R1_A3]
    A1v, d1a, d1b, _, _ = fit_phase(rs1, (gs1[:, i] - 1) * rs1, W_LO, W_HI)
    A2v, d2a, d2b, _, _ = fit_phase(rs2, (gs2[:, i] - 1) * rs2, W_LO, W_HI)
    relA = abs(A2v - A1v) / A1v
    dd = abs((d2a if idx == 1 else d2b) - (d1a if idx == 1 else d1b))
    ok = relA <= 0.01 and dd <= 0.5
    if not ok:
        r_ok = False
    print("    %-3s: dA/A=%.2e  |d(delta)|=%.4f deg  %s"
          % (name, relA, dd, "OK" if ok else "FAIL"))
print("    UWAGA metodyczna: ODE calkowane jako IVP na zewnatrz (bez"
      " sciany pudla); kontrola R testuje dlugosc calkowania, nie"
      " odbicie Dirichleta. Dodatkowo czulosc na okno:")
for wlo, whi in ((40.0, 150.0), (70.0, 150.0), (40.0, 110.0)):
    rs1, gs1, _ = res_R[R0_A3]
    A_e2, de2a, de2b, _, _ = fit_phase(rs1, (gs1[:, 0] - 1) * rs1, wlo, whi)
    print("    okno [%3.0f,%3.0f]: A_e=%.6f delta_e=%+.2f deg"
          % (wlo, whi, A_e2, de2a if idx == 1 else de2b))

a3_pass = ph_ok and amp_ok and r_ok
print("  A3: %s" % ("PASS" if a3_pass else "FAIL"))

# =====================================================================
# A5: rozstrzygniecie niespojnosci istnienia (dodatekH vs AUDYT Dod. A)
# =====================================================================
print("\n" + "-" * 78)
print("[A5] niespojnosc istnienia: maszyneria 2 vs AUDYT_KRYTYCZNY Dod. A")
print("-" * 78)
print("""  ROWNANIE-PO-ROWNANIU (algebra):
  (M2)  maszyneria 2 (dodatekH/O-L5, wszystkie skrypty v47b):
        g'' + 2 g'/r + (alpha/g) g'^2 = + g^2 (1-g),  alpha=2
        = rownanie EL funkcjonalu E=int[ g^(2a) g'^2/2 + Wt(g) ] r^2 dr
          z Wt(g) = g^7/7 - g^8/8  (bo Wt'(g) = g^4 * g^2(1-g)),
          Wt''(1) = -1 < 0  -> proznia STATYCZNIE 'tachionowa'
          -> linearyzacja ogona: h'' + 2h'/r + h = 0 (OSCYLACJE, w=1).
  (AUD) AUDYT_KRYTYCZNY Dodatek A ('F-A kanoniczna'):
        u'' + 2 u'/r + 2 u'^2/u = u^3 - u^2 = - u^2 (1-u)
        = rownanie EL z W(u) = u^8/8 - u^7/7 (W' = u^7-u^6),
          W''(1) = +1 > 0  -> proznia stabilna
          -> h'' + 2h'/r - h = 0 (ZANIK wykladniczy; profil z g0>1
          nie zawraca -> runaway do niesk. -> brak solitonu).
  ROZNICA: wylacznie ZNAK czlonu zrodlowego (W -> -W). Zmienna,
  czlon kinetyczny (2u'^2/u), wymiar, BC -- identyczne.
  UWAGA 3. wariant: EL akcji F-S (W_FS = g^3/3-g^4/4) z K=g^4 daje
        u'' + 2u'/r + 2u'^2/u = (1-u)/u^2  -- INNE od obu powyzszych;
  'kanoniczna Form A' rdzenia NIE jest wiec EL(W_FS, K=g^4), tylko
  definiuje zrodlo g^2(1-g) na poziomie rownania po podzieleniu
  przez K. AUDYT wyprowadzil W'=u^7-u^6, tj. potencjal o znaku
  przeciwnym do implicit Wt maszynerii 2. Oba uklady sa wewnetrznie
  spojne; opisuja ROZNE teorie (roznica zlokalizowana).""")

print("  Weryfikacja numeryczna (wlasny RK4):")
print("  -- (M2) Form A alpha=2, d=3 (w zmiennej f):")
test_g0 = [1.2491, 1.5000, 2.0212, 3.1891]
coll, fmin, _ = rk4_f_batch(2.0, 3, np.array([g ** 3 for g in test_g0]),
                            300.0, 0.005)
m2_ok = True
for g0v, c, fm in zip(test_g0, coll, fmin):
    expect = "kolaps" if g0v > 1.6 else "ograniczone+oscylacje"
    got = "kolaps (f->0)" if c else \
        "ograniczone (min g=%.3f)" % fm ** (1 / 3.0)
    ok = c == (g0v > 1.6)
    if not ok:
        m2_ok = False
    print("     g0=%.4f: %s  [oczekiwane: %s]  %s"
          % (g0v, got, expect, "OK" if ok else "NIEZGODNE"))

print("  -- (AUD) u'' + 2u'/r = u^3 - u^2 - 2u'^2/u:")


def audyt_fa(u0, r_max=60.0, h=0.001, cap=50.0):
    r0 = 1e-4
    u, up = u0, 0.0
    src0 = u0 ** 3 - u0 ** 2
    up = src0 * r0 / 3.0
    u = u0 + 0.5 * src0 * r0 * r0 / 3.0

    def rhs(r, uv, upv):
        uc = max(uv, 1e-9)
        return upv, (uc ** 3 - uc ** 2 - 2.0 * upv * upv / uc
                     - 2.0 * upv / r)

    r = r0
    nstep = int(np.ceil((r_max - r0) / h))
    for _ in range(nstep):
        k1u, k1p = rhs(r, u, up)
        k2u, k2p = rhs(r + h / 2, u + h / 2 * k1u, up + h / 2 * k1p)
        k3u, k3p = rhs(r + h / 2, u + h / 2 * k2u, up + h / 2 * k2p)
        k4u, k4p = rhs(r + h, u + h * k3u, up + h * k3p)
        u += h / 6 * (k1u + 2 * k2u + 2 * k3u + k4u)
        up += h / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
        r += h
        if u > cap:
            return True, u, r      # runaway do gory
        if u < 1e-3:
            return True, u, r      # ucieczka do 0 (tez brak solitonu)
    return False, u, r


aud_ok = True
for g0v in test_g0:
    run, u_end, r_end = audyt_fa(g0v)
    if not run:
        aud_ok = False
    print("     u0=%.4f: %s (u=%.2f przy r=%.2f)  [AUDYT: runaway=TRUE, "
          "u_end~11.5]  %s"
          % (g0v, "RUNAWAY" if run else "ograniczone", u_end, r_end,
             "OK" if run else "NIEZGODNE"))

a5_resolved = m2_ok and aud_ok
print("""  ROZSTRZYGNIECIE: roznica ZLOKALIZOWANA (znak zrodla / W -> -W);
  obiekt, ktorego ogon idzie do Phase 1, to rownanie (M2)
  [g'' + 2g'/r + (alpha/g)g'^2 = g^2(1-g)] -- to z niego pochodza
  g0_crit=8/5, omega=1, fazy delta i A_g uzywane w dodatekH/O-L5.
  'F-A kanoniczna' z AUDYT Dod. A to INNY uklad (przeciwny znak
  potencjalu); jego runaway nie przeczy istnieniu solitonow (M2).
  Otwarta ryska (raportowana, nie gate): ktory znak W wynika z AKCJI
  TGP -- rdzen definiuje rownanie, nie wyprowadza W jawnie.""")
print("  A5: %s" % ("ROZSTRZYGNIETE" if a5_resolved
                    else "FAIL (zachowania niezgodne z dokumentacja)"))

# =====================================================================
# A6: korekta p134e-g (uczciwosc -- obowiazkowe powtorzenie)
# =====================================================================
print("\n" + "-" * 78)
print("[A6] korekta p134e-g (powtorzona jawnie, LOCK A6):")
print("-" * 78)
print("""  1. A_tau(g0^tau) jest MONOTONICZNIE ROSNACA -- nie istnieje
     fizyczne maksimum amplitudy.
  2. r_31 = 3477.5 z p131 to KRES SKANU (ograniczenie g0^tau <= 4.0),
     NIE fizyczne maksimum. ZAKAZ cytowania r_31 jako 'wyprowadzonego'.
  3. Delta(mu->tau) NIE osiaga +-120 deg dla zadnego g0^tau.
  4. Q_K = 3/2 jest PARAMETREM WEJSCIOWYM TGP (12 sciezek wyprowadzenia
     zawiodlo -- status_map O-L5); kazdy wynik zalezny od Q_K
     flagowany jako INPUT. g0_tau=4 (K-samodzielnosc) uzyte w A3/Phase1
     jest od tego lancucha zalezne -> flagowane.""")

# =====================================================================
# WERDYKT BRAMKI
# =====================================================================
print("\n" + "=" * 78)
print("WERDYKT BRAMKI (LOCK: A1-A3 PASS i A5 rozstrzygniete -> Phase 1):")
print("  A1: %s" % ("PASS" if a1_pass else "FAIL"))
print("  A2: %s" % ("PASS" if a2_pass else "FAIL"))
print("  A3: %s" % ("PASS" if a3_pass else "FAIL"))
print("  A5: %s" % ("ROZSTRZYGNIETE" if a5_resolved else "FAIL"))
gate = a1_pass and a2_pass and a3_pass and a5_resolved
print("  BRAMKA: %s" % ("OTWARTA -> Phase 1"
                        if gate else "ZAMKNIETA -> STOP CALEGO CYKLU"))
print("  (A4: tabela deskryptywna w PhaseA_report.md, outputy w"
      " PhaseA_A4_output_*.txt; A6: blok wyzej -- obowiazkowe"
      " niezaleznie od werdyktu.)")
print("=" * 78)
