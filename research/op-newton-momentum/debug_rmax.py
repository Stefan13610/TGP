"""Test convergence of c_2 fit vs r_max."""
import numpy as np
from m9_1_static import gaussian_source, solve_static
from scipy.optimize import curve_fit

sigma = 1.0
M, q = 1.0, 1.0
rho = gaussian_source(M, sigma)
A_0 = q*M/(4*np.pi)

print('R_max  n_pts    A_eff           a_2             a_3            c_2')
print('-'*80)
for R, n in [(100,1500),(200,2000),(400,3000),(800,5000)]:
    sol = solve_static(beta=0.0, q=q, rho_func=rho, M_hint=M, sigma_hint=sigma,
                       r_max=R, n_pts=n, linearized=False)
    r_test = np.linspace(8.0, min(0.4*R,80.0), 500)
    eps_n = sol.sol(r_test)[0] / r_test
    def ans(r, a1, a2, a3, a4, a5):
        return a1/r + a2/r**2 + a3/r**3 + a4/r**4 + a5/r**5
    p, _ = curve_fit(ans, r_test, eps_n, p0=[A_0,0,0,0,0])
    a1, a2, a3 = p[0], p[1], p[2]
    c2 = a2/a1**2
    print(f'{R:5.0f}  {n:5d}  {a1:.7f}    {a2:+.5e}    {a3:+.5e}   {c2:+.5f}')
