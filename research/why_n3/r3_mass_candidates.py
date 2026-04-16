"""
R3: KANDYDACI NA FUNKCJONAL MASY -- ktory daje PDG ratio?

Kontext:
  - Empiryczne: m_phys/m_e = (A_tail_i/A_tail_e)^4  dziala do 0.24% (leptony)
  - Ale mechanizm jest nieznany: dlaczego A^4?
  - M_energy (pelna energia solitonu ~ A^2) NIE jest masa fizyczna.

PLAN: systematycznie testujemy kilku kandydatow funkcjonalu F[g]
  i patrzymy ktorego ratio F_mu/F_e i F_tau/F_e najlepiej pasuje do PDG.

KANDYDACI:
  H1: A_tail^4                                        [baseline, dziala 0.24%]
  H2: A_tail^p, p fitted                              [czy exponent = 4 exact?]
  H3: int (g-1)^4 r^2 dr                              [4th moment radial]
  H4: int (g-1)^2 (g')^2 r^2 dr                       [kinematics x amplitude]
  H5: int V_eff^2 r^2 dr                              [potentiel squared]
  H6: int [T^2 + V^2] r^2 dr                          [energy density squared]
  H7: int (1-g)^2 g^2 r^2 dr                          [potential (substrat) calka]
  H8: int [g''(r)]^2 r^2 dr                           [curvature squared]
  H9: (g0 - 1)^4                                      [depth only, no ODE]
  H10: A_tail^2 * <(g-1)^2>                           [hybrid]
  H11: int (g-1)^4 dr (without r^2 weight)            [linear 4th moment]
  H12: int (1-g)^4 r^4 dr                             [higher weight]
  H13: M_energy * A_tail^2                            [energy * amplitude^2]
  H14: (max|g-1|)^4                                   [amplitude only]

Dla kazdej hipotezy liczymy:
  - ratio mu/e, tau/e
  - diff% vs PDG
  - ile parametrow wolnych (powinno byc 0!)

CEL: znalezc kandydata zero-parameter, ktory jest dokladniejszy niz A^4
     (lub potwierdzic ze A^4 jest optimum).
"""

import numpy as np
from scipy.integrate import solve_ivp


def solve_soliton(g0, alpha=1.0, d=3, r_max=80.0):
    def ode(r, y):
        g, gp = y
        g = max(g, 1e-12)
        rhs = (1.0 - g) * g**(2.0 - 2.0*alpha)
        return [gp, rhs - (alpha/g)*gp*gp - ((d-1)/r)*gp]

    r0 = 1e-4
    c2 = (1.0 - g0) * g0**(2.0 - 2.0*alpha) / (2*d)
    y0 = [g0 + c2*r0*r0, 2*c2*r0]
    return solve_ivp(ode, (r0, r_max), y0, method='DOP853',
                      dense_output=True, rtol=1e-11, atol=1e-13, max_step=0.05)


def extract_A_tail(sol, r_fit_min=20.0):
    """Fit g-1 = A*sin(r+delta)/r for tail."""
    rs = np.linspace(r_fit_min, sol.t[-1] - 2, 2000)
    ys = sol.sol(rs)
    u = (ys[0] - 1) * rs
    S, C = np.sin(rs), np.cos(rs)
    M = np.column_stack([S, C])
    coefs, _, _, _ = np.linalg.lstsq(M, u, rcond=None)
    A_s, A_c = coefs
    return np.sqrt(A_s**2 + A_c**2)


def compute_functionals(sol, alpha=1.0, d=3, r_outer=70.0):
    """Compute all candidate functionals on a soliton solution."""
    rs = np.linspace(max(sol.t[0], 1e-4), r_outer, 8000)
    ys = sol.sol(rs)
    g = np.clip(ys[0], 1e-12, 10.0)
    gp = ys[1]
    r2 = rs**2

    # Higher derivative via ODE
    # g'' = (1-g)*g^(2-2a) - (a/g)*gp^2 - ((d-1)/r)*gp
    gpp = (1.0 - g) * g**(2.0 - 2.0*alpha) - (alpha/g)*gp**2 - ((d-1)/rs)*gp

    # Densities
    T = 0.5 * g**(2*alpha) * gp**2
    V = g**3/3.0 - g**4/4.0
    V1 = 1.0/12.0
    V_eff = V - V1
    # Substrat potential (raw)
    V_sub = (1.0 - g)**2 * g**2  # oryginalny potential (1-g)^2 g^2

    # g-1 minimum (peak below 1)
    g_dev = g - 1.0

    return {
        # --- funkcjonalne ---
        # Radialne z r^2 (3D standard)
        'I_gm1_2_r2':   np.trapezoid(g_dev**2 * r2, rs),            # int (g-1)^2 r^2
        'I_gm1_3_r2':   np.trapezoid(g_dev**3 * r2, rs),            # int (g-1)^3 r^2
        'I_gm1_4_r2':   np.trapezoid(g_dev**4 * r2, rs),            # int (g-1)^4 r^2
        'I_gm1_4_r4':   np.trapezoid(g_dev**4 * rs**4, rs),         # int (g-1)^4 r^4
        'I_gm1_4_r0':   np.trapezoid(g_dev**4, rs),                 # int (g-1)^4 (no weight)
        'I_gm1_2_gp2_r2': np.trapezoid(g_dev**2 * gp**2 * r2, rs),  # hybrid
        'I_Veff_r2':    np.trapezoid(V_eff * r2, rs),               # int V_eff r^2
        'I_Veff_sq_r2': np.trapezoid(V_eff**2 * r2, rs),            # int V_eff^2 r^2
        'I_Vsub_r2':    np.trapezoid(V_sub * r2, rs),               # int (1-g)^2 g^2 r^2
        'I_T_r2':       np.trapezoid(T * r2, rs),                   # int T r^2
        'I_T2_r2':      np.trapezoid(T**2 * r2, rs),                # int T^2 r^2
        'I_TV_r2':      np.trapezoid(T * V_eff * r2, rs),           # int T*V_eff r^2
        'I_gpp2_r2':    np.trapezoid(gpp**2 * r2, rs),              # int (g'')^2 r^2
        'M_energy':     4*np.pi*np.trapezoid((T + V_eff) * r2, rs), # full energy
        # --- szczytowe (peak) ---
        'abs_dev_max':  float(np.max(np.abs(g_dev))),                # max |g-1|
        'g_min':        float(np.min(g)),                            # min g
        'g_max':        float(np.max(g)),                            # max g
    }


def compute_all(g0, alpha=1.0, d=3, r_max=80.0, r_outer=70.0):
    sol = solve_soliton(g0, alpha, d, r_max)
    if not sol.success:
        return None
    A = extract_A_tail(sol)
    F = compute_functionals(sol, alpha, d, r_outer)
    F['A_tail'] = A
    F['g0'] = g0
    return F


def main():
    print("=" * 76)
    print("  R3 KANDYDACI NA FUNKCJONAL MASY: ktory daje PDG m_mu/m_e?")
    print("=" * 76)
    print()
    print("  PDG reference:")
    m_mu_pdg = 206.7682830
    m_tau_pdg = 3477.15
    print(f"    m_mu/m_e  = {m_mu_pdg:.6f}")
    print(f"    m_tau/m_e = {m_tau_pdg:.6f}")
    print()

    # Bridge calibration (g0 values from R5 electron/muon/tau substrate)
    g0_e, g0_mu, g0_tau = 0.86941, 1.40673, 1.72931
    print(f"  Substrat g0:")
    print(f"    g0_e   = {g0_e}")
    print(f"    g0_mu  = {g0_mu}")
    print(f"    g0_tau = {g0_tau}")
    print()

    print("  Solving solitons...")
    F_e = compute_all(g0_e)
    F_mu = compute_all(g0_mu)
    F_tau = compute_all(g0_tau)

    if F_e is None or F_mu is None or F_tau is None:
        print("  SOLVER FAILED")
        return

    print()
    print(f"  A_tail_e   = {F_e['A_tail']:.6f}")
    print(f"  A_tail_mu  = {F_mu['A_tail']:.6f}")
    print(f"  A_tail_tau = {F_tau['A_tail']:.6f}")
    print()

    # ======================================================================
    # Lista kandydatow do testu: (nazwa, funkcja F -> float, wykladnik p)
    # Testujemy F^p dla p = 1, 2, 4 jesli ma to sens
    # ======================================================================

    # Zbierzemy: f_e, f_mu, f_tau oraz ratios

    candidates = [
        # (label, key, exponent_p)  -> m = F^p
        ("A_tail^4",               'A_tail',           4),
        ("A_tail^2",               'A_tail',           2),
        ("A_tail^1",               'A_tail',           1),
        ("(g0-1)^4",               'g0-1',             4),
        ("I_gm1_4_r2^1",           'I_gm1_4_r2',       1),
        ("I_gm1_4_r4^1",           'I_gm1_4_r4',       1),
        ("I_gm1_4_r0^1",           'I_gm1_4_r0',       1),
        ("I_gm1_2_gp2_r2^1",       'I_gm1_2_gp2_r2',   1),
        ("(I_gm1_2_r2)^2",         'I_gm1_2_r2',       2),
        ("I_Vsub_r2^2",            'I_Vsub_r2',        2),
        ("|I_Veff_r2|^2",          'I_Veff_r2',        2),
        ("I_Veff_sq_r2^1",         'I_Veff_sq_r2',     1),
        ("I_T_r2^2",               'I_T_r2',           2),
        ("I_T2_r2^1",              'I_T2_r2',          1),
        ("|I_TV_r2|^1",            'I_TV_r2',          1),
        ("|I_TV_r2|^2",            'I_TV_r2',          2),
        ("I_gpp2_r2^1",            'I_gpp2_r2',        1),
        ("M_energy^2",             'M_energy',         2),
        ("abs_dev_max^4",          'abs_dev_max',      4),
    ]

    def get_F(dat, key):
        if key == 'g0-1':
            return abs(dat['g0'] - 1.0)
        return abs(dat[key])  # wez |.| gdzie potrzeba (np. V_eff moze byc ujemne)

    print("=" * 76)
    print("  KANDYDACI: ratio_mu = F_mu/F_e, ratio_tau = F_tau/F_e")
    print("  podniesione do wykladnika p")
    print("=" * 76)
    print()
    print(f"  {'Kandydat':<22} {'mu/e':>12} {'PDG diff%':>10}    {'tau/e':>12} {'PDG diff%':>10}")
    print(f"  {'-'*22} {'-'*12} {'-'*10}    {'-'*12} {'-'*10}")

    results = []
    for label, key, p in candidates:
        fe = get_F(F_e, key)
        fm = get_F(F_mu, key)
        ft = get_F(F_tau, key)
        if fe < 1e-20:
            continue
        # ratio podniesiony do p
        ratio_mu = (fm/fe)**p
        ratio_tau = (ft/fe)**p
        diff_mu = (ratio_mu - m_mu_pdg)/m_mu_pdg * 100
        diff_tau = (ratio_tau - m_tau_pdg)/m_tau_pdg * 100
        score = abs(diff_mu) + abs(diff_tau)
        results.append((label, ratio_mu, diff_mu, ratio_tau, diff_tau, score))
        print(f"  {label:<22} {ratio_mu:>12.4f} {diff_mu:>+9.3f}%    "
              f"{ratio_tau:>12.4f} {diff_tau:>+9.3f}%")

    print()
    print("=" * 76)
    print("  RANKING (posortowany po |diff_mu| + |diff_tau|)")
    print("=" * 76)
    print()
    print(f"  {'Kandydat':<22} {'mu/e diff%':>12} {'tau/e diff%':>14} {'score':>10}")
    print(f"  {'-'*22} {'-'*12} {'-'*14} {'-'*10}")
    results.sort(key=lambda x: x[5])
    for label, rmu, dmu, rtau, dtau, score in results:
        print(f"  {label:<22} {dmu:>+11.3f}% {dtau:>+13.3f}% {score:>10.3f}")

    # ======================================================================
    # FIT WYKLADNIKA: dla kazdej 'bazy' F, dopasujmy p tak aby
    #   (F_mu/F_e)^p = m_mu/m_e AND (F_tau/F_e)^p = m_tau/m_e
    # Jesli to jest spojne (jedno p dziala dla mu i tau), mamy hipoteze.
    # ======================================================================
    print()
    print("=" * 76)
    print("  FIT WYKLADNIKA p: dla kazdego F, znajdz p takie ze F^p daje PDG")
    print("  (jesli p_mu ~ p_tau, kandydat jest KONSYSTENTNY)")
    print("=" * 76)
    print()
    print(f"  {'Baza F':<22} {'p_mu':>10} {'p_tau':>10} {'|dp|':>10}  {'status':<15}")
    print(f"  {'-'*22} {'-'*10} {'-'*10} {'-'*10}  {'-'*15}")

    bases = [
        ('A_tail',          'A_tail'),
        ('g0-1',            'g0-1'),
        ('I_gm1_2_r2',      'I_gm1_2_r2'),
        ('I_gm1_4_r2',      'I_gm1_4_r2'),
        ('I_gm1_4_r4',      'I_gm1_4_r4'),
        ('I_gm1_2_gp2_r2',  'I_gm1_2_gp2_r2'),
        ('I_Vsub_r2',       'I_Vsub_r2'),
        ('|I_Veff_r2|',     'I_Veff_r2'),
        ('I_Veff_sq_r2',    'I_Veff_sq_r2'),
        ('I_T_r2',          'I_T_r2'),
        ('I_T2_r2',         'I_T2_r2'),
        ('|I_TV_r2|',       'I_TV_r2'),
        ('I_gpp2_r2',       'I_gpp2_r2'),
        ('M_energy',        'M_energy'),
        ('abs_dev_max',     'abs_dev_max'),
    ]

    consistent = []
    for label, key in bases:
        fe = get_F(F_e, key)
        fm = get_F(F_mu, key)
        ft = get_F(F_tau, key)
        if fe < 1e-20 or fm < 1e-20 or ft < 1e-20:
            continue
        # (F_mu/F_e)^p = m_mu/m_e -> p = log(m_mu/m_e) / log(F_mu/F_e)
        lrm = np.log(fm/fe)
        lrt = np.log(ft/fe)
        if abs(lrm) < 1e-6 or abs(lrt) < 1e-6:
            continue
        p_mu = np.log(m_mu_pdg) / lrm
        p_tau = np.log(m_tau_pdg) / lrt
        dp = abs(p_mu - p_tau)
        rel_dp = dp / (abs(p_mu) + 1e-15) * 100
        status = "KONSYST." if rel_dp < 0.5 else ("blisko" if rel_dp < 2.0 else "rozne")
        print(f"  {label:<22} {p_mu:>10.5f} {p_tau:>10.5f} {dp:>10.5f}  {status:<15}")
        if rel_dp < 2.0:
            consistent.append((label, p_mu, p_tau, dp))

    # ======================================================================
    # PODSUMOWANIE
    # ======================================================================
    print()
    print("=" * 76)
    print("  WNIOSKI")
    print("=" * 76)
    print()

    best = results[0]
    print(f"  1. NAJLEPSZY kandydat: {best[0]}")
    print(f"     mu/e diff: {best[2]:+.3f}%")
    print(f"     tau/e diff: {best[4]:+.3f}%")
    print(f"     score: {best[5]:.3f}")
    print()

    if consistent:
        print(f"  2. KONSYSTENTNE WYKLADNIKI (p_mu ~ p_tau, < 2% diff):")
        for label, pm, pt, dp in consistent:
            print(f"     {label:<22} p_mu={pm:.4f} p_tau={pt:.4f}")
    else:
        print(f"  2. Brak kandydatow z konsystentnym wykladnikiem.")
    print()

    print("  3. TEZA DO WERYFIKACJI:")
    print("     Jesli JEDYNIE A_tail^4 daje konsystentne p=4 dla mu i tau,")
    print("     to A_tail jest UNIVERSAL soliton invariant a p=4 jest SCALE DIMENSION")
    print("     tj. masa fizyczna transformuje sie jak 4-ty moment tailu.")
    print()
    print("  4. Znaczenie fizyczne: w teoriach efekt. (Skyrme, vortex) masa soliton")
    print("     jest 1/lambda^D (D=dim scaling), gdzie lambda = promien core.")
    print("     Tutaj: D=4 dla R3/R5. To 'brane-like' dim counting: soliton w 3D")
    print("     ktorego rozmiar core zmienia się jak 1/A, masa ~ 1/(core size)^4")
    print("     = A^4. => DIMENSION COUNTING D=4 mozliwa motywacja.")
    print()


if __name__ == "__main__":
    main()
