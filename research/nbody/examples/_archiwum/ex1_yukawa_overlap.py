"""
ex1_yukawa_overlap.py
=====================
Numeryczna analiza trójcentrowego overlapu Yukawy w TGP.

Implementuje DOKŁADNĄ jednowymiarową reprezentację:
    f_△(t) = 2 ∫₀^{1/3} p(Δ)·Δ^{-3/2}·K₀(t/√Δ) dΔ

gdzie dokładna gęstość p(Δ) na 2-simpleksie wynosi:
    p(Δ) = 2π/√3                                    dla Δ∈[1/4, 1/3]
    p(Δ) = (2/√3)[3·arcsin(1/(2√(1-3Δ))) - π/2]   dla Δ∈[0,   1/4]

Wyniki:
  - Wykres f_△(t) vs t
  - Wykres gęstości p(Δ)
  - Porównanie z asymptotypiami:
      t→∞: f_△ ~ C_inf · t^{-3/2} · exp(-√3·t)
      t→0: f_△ ~ -4π·ln(t) + C_0
  - Weryfikacja tożsamości ∫ Δ^{-3/2} dα = 2π
"""

import numpy as np
from scipy.integrate import quad
from scipy.special import k0 as K0_scipy
import matplotlib.pyplot as plt
import os

# ── Gęstość p(Δ) ──────────────────────────────────────────────────────

def p_density(Delta):
    """Dokładna gęstość miary uniformowej na 2-simpleksie dla Δ."""
    Delta = np.asarray(Delta, dtype=float)
    result = np.zeros_like(Delta)
    const = 2.0 * np.pi / np.sqrt(3.0)

    high = Delta >= 0.25
    result[high] = const

    low = (Delta > 0) & (Delta < 0.25)
    xi = 1.0 / (2.0 * np.sqrt(1.0 - 3.0 * Delta[low]))
    xi = np.clip(xi, -1.0, 1.0)
    result[low] = (2.0 / np.sqrt(3.0)) * (3.0 * np.arcsin(xi) - np.pi / 2.0)

    return result


# ── Całka f_△(t) — dokładna 1D ────────────────────────────────────────

def f_triangle_exact(t, n_quad=2000):
    """
    f_△(t) = 2 ∫₀^{1/3} p(Δ)·Δ^{-3/2}·K₀(t/√Δ) dΔ

    Całkowanie adaptacyjne z osobną obsługą punktów osobliwych.
    """
    def integrand(Delta):
        if Delta <= 0:
            return 0.0
        arg = t / np.sqrt(Delta)
        if arg > 600:
            return 0.0
        return 2.0 * p_density(np.array([Delta]))[0] * Delta**(-1.5) * K0_scipy(arg)

    # Osobliwość w Δ=0 (gęstość→0 szybciej niż Δ^{3/2}) — bezpieczna
    # Osobliwość w Δ=1/4 (nieciągłość pochodnej p) — podziel przedział
    eps = 1e-6
    val1, _ = quad(integrand, eps, 0.25 - eps, limit=200, epsrel=1e-8)
    val2, _ = quad(integrand, 0.25 + eps, 1.0/3.0 - eps, limit=200, epsrel=1e-8)
    return val1 + val2


# ── Asymptotyki ────────────────────────────────────────────────────────

def f_triangle_large_t(t):
    """Asymptotyka Yukawa t→∞: f ~ C_inf · t^{-3/2} · exp(-√3·t)"""
    C_inf = 8.0 * np.pi * np.sqrt(np.pi) / (np.sqrt(2.0) * 3.0**0.75)
    return C_inf * t**(-1.5) * np.exp(-np.sqrt(3.0) * t)


def f_triangle_small_t(t, C0=None):
    """Asymptotyka Coulomba t→0: f ~ -4π·ln(t) + C0"""
    if C0 is None:
        C0 = 13.65
    return -4.0 * np.pi * np.log(t) + C0


# ── Weryfikacja tożsamości C₁ = 2π ────────────────────────────────────

def verify_identity_C1():
    """Weryfikuje ∫_Δ₂ Δ^{-3/2} dα = 2π przez całkowanie 1D."""
    def integrand(Delta):
        if Delta <= 1e-8:
            return 0.0
        return p_density(np.array([Delta]))[0] * Delta**(-1.5)

    val, err = quad(integrand, 1e-6, 1.0/3.0 - 1e-8, limit=500, epsrel=1e-10)
    print(f"C₁ = ∫ p(Δ)·Δ^(-3/2) dΔ  =  {val:.8f}")
    print(f"Oczekiwane 2π             =  {2*np.pi:.8f}")
    print(f"Błąd względny             =  {abs(val - 2*np.pi)/(2*np.pi):.2e}")
    return val


# ── Plots ──────────────────────────────────────────────────────────────

def plot_density():
    """Wykres gęstości p(Δ) na [0, 1/3]."""
    Delta_vals = np.linspace(1e-4, 1.0/3.0 - 1e-4, 800)
    p_vals = p_density(Delta_vals)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(Delta_vals, p_vals, 'b-', linewidth=2, label=r'$p(\Delta)$')
    ax.axvline(0.25, color='r', linestyle='--', alpha=0.7, label=r'$\Delta=1/4$ (punkt krytyczny)')
    ax.axvline(1.0/3.0, color='g', linestyle='--', alpha=0.7, label=r'$\Delta=1/3$ (centrum, równoboczny)')
    ax.axhline(2*np.pi/np.sqrt(3), color='gray', linestyle=':', alpha=0.7,
               label=r'$2\pi/\sqrt{3}$ (wartość stała)')
    ax.set_xlabel(r'$\Delta = \alpha_1\alpha_2 + \alpha_1\alpha_3 + \alpha_2\alpha_3$', fontsize=12)
    ax.set_ylabel(r'$p(\Delta)$', fontsize=12)
    ax.set_title('Gęstość miary uniformowej na 2-simpleksie TGP', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plots_dir = os.path.join(os.path.dirname(__file__), '..', 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    out = os.path.join(plots_dir, 'ex1_density.png')
    plt.tight_layout()
    plt.savefig(out, dpi=120)
    plt.close()
    print(f"Zapisano: {os.path.abspath(out)}")


def plot_f_triangle():
    """Wykres f_△(t) z asymptotypiami."""
    t_vals = np.linspace(0.15, 6.0, 40)
    print("Obliczam f_triangle(t) - może chwilę potrwać...")
    f_vals = np.array([f_triangle_exact(t) for t in t_vals])

    t_large = np.linspace(2.0, 6.0, 100)
    t_small = np.linspace(0.05, 0.8, 60)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Lewy: pełny zakres
    ax = axes[0]
    ax.semilogy(t_vals, f_vals, 'b-o', markersize=4, label=r'$f_\triangle(t)$ dokładna')
    ax.semilogy(t_large, f_triangle_large_t(t_large), 'r--',
                label=r'Asymptotyka Yukawa: $C_\infty t^{-3/2}e^{-\sqrt{3}t}$')
    ax.set_xlabel(r'$t = m\cdot d$', fontsize=12)
    ax.set_ylabel(r'$f_\triangle(t)$', fontsize=12)
    ax.set_title('Overlap Yukawy — tr. równoboczny', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Prawy: małe t (reżim Coulomba)
    ax = axes[1]
    t_small_exact = t_vals[t_vals <= 1.5]
    f_small_exact = f_vals[t_vals <= 1.5]

    # Kalibracja C0 z jednego punktu
    if len(f_small_exact) > 0:
        t_ref = t_small_exact[0]
        C0_cal = f_small_exact[0] + 4.0 * np.pi * np.log(t_ref)
        ax.plot(t_small_exact, f_small_exact, 'b-o', markersize=4, label=r'$f_\triangle(t)$ dokładna')
        ax.plot(t_small, f_triangle_small_t(t_small, C0=C0_cal), 'g--',
                label=rf'Coulomb: $-4\pi\ln t + {C0_cal:.2f}$')
    ax.set_xlabel(r'$t = m\cdot d$', fontsize=12)
    ax.set_ylabel(r'$f_\triangle(t)$', fontsize=12)
    ax.set_title('Limit Coulomba ($t \\to 0$)', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plots_dir = os.path.join(os.path.dirname(__file__), '..', 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    out = os.path.join(plots_dir, 'ex1_f_triangle.png')
    plt.tight_layout()
    plt.savefig(out, dpi=120)
    plt.close()
    print(f"Zapisano: {os.path.abspath(out)}")


# ── Main ───────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("="*60)
    print("EX1: Dokładny overlap Yukawy f_△(t) w TGP")
    print("="*60)

    print("\n[1] Weryfikacja tożsamości C₁ = 2π:")
    C1 = verify_identity_C1()

    print("\n[2] Weryfikacja normalizacji ∫p(Δ)dΔ = 1/2:")
    norm, _ = quad(lambda D: p_density(np.array([D]))[0], 1e-6, 1/3 - 1e-6, limit=300)
    print(f"∫p(Δ)dΔ = {norm:.8f}  (oczekiwane 0.5)")

    print("\n[3] Kilka wartości f_△(t):")
    for t in [0.5, 1.0, 2.0, 3.0, 5.0]:
        fval = f_triangle_exact(t)
        fasy = f_triangle_large_t(t)
        ratio = fasy / fval if fval != 0 else float('nan')
        print(f"  t={t:.1f}: f={fval:.6f}, asympt={fasy:.6f}, stosunek={ratio:.4f}")

    print("\n[4] Rysowanie wykresów...")
    plot_density()
    plot_f_triangle()
    print("\nGotowe.")
