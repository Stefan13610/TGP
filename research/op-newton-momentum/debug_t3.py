"""Debug script: check numerical accuracy of T3 fit by running linearized
case with the SAME fit pipeline. Linearized should give a_2 ~ 0 exactly
(no nonlinear correction). Any nonzero a_2 in linearized = numerical
artifact, sets the noise floor for full nonlinear T3."""
import numpy as np
from m9_1_static import gaussian_source, solve_static
from scipy.optimize import curve_fit

sigma = 1.0
M = 1.0
print(f"{'q':>8} {'A_0':>10} {'A_eff':>10} {'a_2':>14} {'c_2_lin':>12} {'a_2 full':>14} {'c_2_full':>12}")
print('-' * 90)
for q in [0.5, 1.0, 2.0, 3.0]:
    rho = gaussian_source(M, sigma)
    A_0 = q*M/(4*np.pi)

    sol_lin = solve_static(beta=0.0, q=q, rho_func=rho, M_hint=M, sigma_hint=sigma,
                           r_max=200.0, n_pts=2000, linearized=True)
    sol_full = solve_static(beta=0.0, q=q, rho_func=rho, M_hint=M, sigma_hint=sigma,
                            r_max=200.0, n_pts=2000, linearized=False)

    r_test = np.linspace(8.0*sigma, 60.0*sigma, 500)
    eps_lin = sol_lin.sol(r_test)[0] / r_test
    eps_full = sol_full.sol(r_test)[0] / r_test

    def ans(r, a1, a2, a3, a4):
        return a1/r + a2/r**2 + a3/r**3 + a4/r**4

    p_lin, _ = curve_fit(ans, r_test, eps_lin, p0=[A_0, 0, 0, 0])
    p_full, _ = curve_fit(ans, r_test, eps_full, p0=[A_0, 0, 0, 0])

    a1_lin, a2_lin = p_lin[0], p_lin[1]
    a1_full, a2_full = p_full[0], p_full[1]
    c2_lin = a2_lin/a1_lin**2
    c2_full = a2_full/a1_full**2
    print(f"{q:8.3f} {A_0:10.5f} {a1_full:10.6f} {a2_lin:+14.4e} {c2_lin:+12.4f} {a2_full:+14.4e} {c2_full:+12.4f}")
