#!/usr/bin/env python3
# -*- coding: ascii -*-
"""op-collapse-matter-source -- silnik dynamiki 2. rzedu.

KOPIA silnika `../op-r3-stationary-states-2026-09-14/engine_core.py`
(FROZEN, MD sec. 2-4 tamtego cyklu; cytat nizej). JEDYNA zmiana
merytoryczna (LOCK sec. 3 tego cyklu): czlon materii
  -dU_mat/dpsi w RHS (punktowo w F) oraz [U_mat(psi,r)-U_mat(1,r)]
  w ewaluatorze energii (tozsamosc sfaktoryzowana bez kancelacji --
  dziedziczona lekcja correction note 2 poprzednika).
Dla lam=0 silnik jest identyczny z poprzednikiem (dodawanie zer).

CYTAT (naglowek poprzednika):
  "Formy CYTAT (MD sec. 1; dziedziczone op-action-audit
   Phase1_output.txt):
     M(psi) = psi^6/(4-3psi)^2 ; M' = 12 psi^5 (2-psi)/(4-3psi)^3
     K(psi) = psi^4 ; K' = 4 psi^3
     U(psi) = psi^4/4 - psi^3/3 ; U' = psi^2(psi-1)
     K_geo = gamma = c0 = 1 [LOCK]. ZERO podlog/barier.
   pi-formulacja (MD sec. 2):
     psidot = pi/M(psi)
     pidot  = pi^2 M'/(2M^2) + F_sp(psi) - gamma_sp(r) pi
   Integrator (MD sec. 3): uogolniony Stormer-Verlet (symetryczny,
   2. rzedu) na parze (psi, pi); kroki implicit przez iteracje punktu
   stalego (stop: stagnacja maszynowa, correction note 1 (b)).
   Sponge (MD sec. 4): gamma_sp(r) = gamma0 * S((r-160)/40),
   S = smootherstep, gamma0 = 1.0 [INPUT-MD].
   Energia (MD sec. 2): prozniowo odjeta, czynnik 4pi; E_core: r<=80."

Czlon materii (WYNIK P1-H1 tego cyklu, Phase1_output.txt, gate PASS):
  U_mat(psi,r) = lam * rho_hat(r) * psi^2/(4-3psi)
  dU_mat/dpsi  = lam * rho_hat(r) * psi(8-3psi)/(4-3psi)^2
  U_mat-U_mat(1) = lam * rho_hat * (psi-1)(psi+4)/(4-3psi)
  rho_hat(r) = exp(-r^2/18) [LOCK, FROZEN]; lam >= 0.
"""
import numpy as np

GAMMA0_SP = 1.0          # INPUT-MD (dziedziczone)
R_SP0 = 160.0            # LOCK
R_SP1 = 200.0            # LOCK
PSI_LIMIT = 4.0/3.0
PSI_BAND_HI = PSI_LIMIT - 1e-6   # LOCK (pas klasyfikacyjny)
PSI_BAND_LO = 1e-6               # LOCK
FP_TOL = 1e-12           # awaryjny prog NonConvergence po cap
FP_MAXIT = 200           # stop = stagnacja maszynowa (dziedziczone)


def Mfun(p):
    return p**6/(4.0 - 3.0*p)**2


def Mpfun(p):
    return 12.0*p**5*(2.0 - p)/(4.0 - 3.0*p)**3


def Kfun(p):
    return p**4


def Kpfun(p):
    return 4.0*p**3


def Ufun(p):
    return p**4/4.0 - p**3/3.0


def Upfun(p):
    return p*p*(p - 1.0)


def dUfun(p):
    """U(psi) - U(1) = (psi-1)^2 (3 psi^2 + 2 psi + 1)/12 --
    tozsamosc DOKLADNA (correction note 2 poprzednika)."""
    d = p - 1.0
    return d*d*(3.0*p*p + 2.0*p + 1.0)/12.0


def rho_hat(r):
    """Zrodlo statyczne rho_hat(r) = exp(-r^2/18) [LOCK, FROZEN]."""
    return np.exp(-r*r/18.0)


def Umat_pfun(p, lr):
    """dU_mat/dpsi = (lam rho_hat) psi(8-3psi)/(4-3psi)^2 (P1-H1)."""
    return lr*p*(8.0 - 3.0*p)/(4.0 - 3.0*p)**2


def dUmatfun(p, lr):
    """U_mat(psi)-U_mat(1) = (lam rho_hat)(psi-1)(psi+4)/(4-3psi)
    (P1-H1; forma sfaktoryzowana bez kancelacji)."""
    return lr*(p - 1.0)*(p + 4.0)/(4.0 - 3.0*p)


class NonConvergence(Exception):
    pass


class Engine:
    """Radialna 3D siatka przesunieta r_i=(i+1/2)h; brzegi: zerowy
    strumien w r=0 i r=R (dziedziczone). lam = sila sprzezenia
    zrodla statycznego (lam=0 -> silnik poprzednika)."""

    def __init__(self, h, R, sponge=True, gamma0=GAMMA0_SP, lam=0.0):
        N = int(round(R/h))
        self.h = h
        self.R = R
        self.N = N
        self.lam = lam
        self.r = (np.arange(N) + 0.5)*h
        self.r2 = self.r**2
        self.rm = (np.arange(N - 1) + 1.0)*h      # r_{i+1/2}
        self.rm2 = self.rm**2
        self.lrho = lam*rho_hat(self.r)           # lam * rho_hat(r_i)
        if sponge:
            x = np.clip((self.r - R_SP0)/(R_SP1 - R_SP0), 0.0, 1.0)
            self.gsp = gamma0*(6*x**5 - 15*x**4 + 10*x**3)
        else:
            self.gsp = np.zeros(N)
        self.has_sponge = sponge and bool(np.any(self.gsp > 0))

    # ---- sila przestrzenna F_sp = -(1/(h r^2)) dE_sp/dpsi ---------
    def F(self, g):
        h = self.h
        gm = 0.5*(g[:-1] + g[1:])
        dg = np.diff(g)/h
        dH = h*self.r2*(Upfun(g) + Umat_pfun(g, self.lrho))
        t_flux = self.rm2*Kfun(gm)*dg
        t_quad = 0.25*h*self.rm2*Kpfun(gm)*dg**2
        dH[:-1] += -t_flux + t_quad
        dH[1:] += t_flux + t_quad
        return -dH/(h*self.r2)

    # ---- energia (prozniowo odjeta, 4pi; opcjonalnie r<=rmax) -----
    def energy(self, g, pi, rmax=None, matter=True):
        """matter=True: gestosc PELNA (zachowana; z czlonem U_mat
        prozniowo odjetym); matter=False: E^field deskryptywnie."""
        gm = 0.5*(g[:-1] + g[1:])
        dg = np.diff(g)/self.h
        dens = self.r2*(pi*pi/(2.0*Mfun(g)) + dUfun(g))
        if matter:
            dens = dens + self.r2*dUmatfun(g, self.lrho)
        grad = 0.5*self.rm2*Kfun(gm)*dg**2
        if rmax is None:
            return 4.0*np.pi*self.h*(float(np.sum(dens))
                                     + float(np.sum(grad)))
        mp = self.r <= rmax
        mg = self.rm <= rmax
        return 4.0*np.pi*self.h*(float(np.sum(dens[mp]))
                                 + float(np.sum(grad[mg])))

    # ---- krok: uogolniony Stormer-Verlet (dziedziczony) -----------
    def step(self, g, pi, dt):
        F1 = self.F(g)
        M1 = Mfun(g)
        Mp1 = Mpfun(g)
        c1 = Mp1/(2.0*M1*M1)
        hdt = 0.5*dt
        # (1) implicit pol-krok pedu -- stop: stagnacja maszynowa
        ph = pi
        d_prev = np.inf
        for _ in range(FP_MAXIT):
            new = pi + hdt*(F1 + ph*ph*c1 - self.gsp*ph)
            d = float(np.max(np.abs(new - ph)))
            ph = new
            if d == 0.0 or d >= d_prev:
                break
            d_prev = d
        else:
            if d > FP_TOL*max(1.0, float(np.max(np.abs(ph)))):
                raise NonConvergence("pi half-step")
        # (2) implicit krok pola -- stop: stagnacja maszynowa
        invM1 = 1.0/M1
        x = g
        d_prev = np.inf
        for _ in range(FP_MAXIT):
            new = g + hdt*ph*(invM1 + 1.0/Mfun(x))
            d = float(np.max(np.abs(new - x)))
            x = new
            if d == 0.0 or d >= d_prev:
                break
            d_prev = d
        else:
            if d > FP_TOL*max(1.0, float(np.max(np.abs(x)))):
                raise NonConvergence("psi step")
        g2 = x
        # (3) explicit pol-krok pedu
        F2 = self.F(g2)
        M2 = Mfun(g2)
        c2 = Mpfun(g2)/(2.0*M2*M2)
        pi2 = ph + hdt*(F2 + ph*ph*c2 - self.gsp*ph)
        return g2, pi2

    def band_status(self, g):
        """Pas klasyfikacyjny (LOCK): zwraca None lub etykiete."""
        if not np.all(np.isfinite(g)):
            return "BREAKDOWN"
        mx = float(np.max(g))
        mn = float(np.min(g))
        if mx >= PSI_BAND_HI:
            return "BREAKDOWN-BOUNDARY"
        if mn <= PSI_BAND_LO:
            return "BREAKDOWN-BOUNDARY-LOWER"
        return None


# ---- starty (LOCK; identyczne definicje jak u poprzednika) --------
def start_gauss(eng, a, sigma):
    return 1.0 + a*np.exp(-eng.r**2/(2.0*sigma**2))


def start_quasiR3(eng, a, sigma_w=15.0):
    return 1.0 + a*np.sinc(eng.r/np.pi)*np.exp(-eng.r**2
                                               / (2.0*sigma_w**2))


def start_vacuum(eng):
    return np.ones(eng.N)
