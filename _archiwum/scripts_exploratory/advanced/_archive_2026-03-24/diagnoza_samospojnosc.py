"""Diagnostyka: max(E(K)/K) vs Lambda dla roznych parametrow"""
import numpy as np

PHI0 = 1.0; GAMMA = 1.0; Q = 1.0

def E_over_K(K, a_gam=0.3, alpha=0.5):
    r = np.linspace(a_gam, 30.0, 2000)
    phi = PHI0 + K * np.exp(-r) / r
    phi = np.maximum(phi, 1e-8)
    dphi = K * np.exp(-r) * (-r - 1.0) / r**2
    E_k  = 4*np.pi * np.trapezoid(0.5 * dphi**2 * (1 + alpha/(PHI0*phi)) * r**2, r)
    V    = GAMMA/3*phi**3 - GAMMA/4*phi**4 - (GAMMA/3 - GAMMA/4)
    E_p  = 4*np.pi * np.trapezoid(V * r**2, r)
    return (E_k + E_p) / K

K_vals = np.linspace(0.05, 8.0, 60)

print("="*65)
print("Maksimum E(K)/K vs Lambda = 4*pi/q")
print("="*65)
print(f"{'a_Gam':>6} {'alpha':>6} {'max(E/K)':>10} {'K*':>6} {'Lambda':>8} {'n_cross':>8}")
print("-"*65)

for a_gam in [0.1, 0.3, 0.5, 1.0]:
    for alpha in [0.1, 1.0, 5.0]:
        EK = np.array([E_over_K(K, a_gam, alpha) for K in K_vals])
        valid = np.isfinite(EK)
        if valid.sum() == 0:
            continue
        max_EK  = np.max(EK[valid])
        K_star  = K_vals[valid][np.argmax(EK[valid])]
        Lambda  = 4*np.pi / Q
        # Liczba przeciec (E/K krzywa vs poziom Lambda)
        diff = EK[valid] - Lambda
        n_cross = sum(
            diff[i]*diff[i+1] < 0
            for i in range(len(diff)-1)
        )
        print(f"{a_gam:>6.1f} {alpha:>6.1f} {max_EK:>10.3f} {K_star:>6.3f} {Lambda:>8.4f} {n_cross:>8d}")

print()
print("Ile q potrzeba by max(E/K) > Lambda (warunek 2 przeciec):")
for a_gam in [0.1, 0.3, 0.5]:
    alpha = 1.0
    EK = np.array([E_over_K(K, a_gam, alpha) for K in K_vals])
    max_EK = np.max(EK[np.isfinite(EK)])
    q_min = 4*np.pi / max_EK
    print(f"  a_Gam={a_gam}: max(E/K)={max_EK:.3f} => q_min={q_min:.4f} (aktualnie q={Q})")

print()
print("WNIOSEK:")
print("  Przy q=1: Lambda=12.57 >> max(E/K)")
print("  Pojedynczy Yukawa nie moze dac 3 przeciec.")
print("  Potrzeba: kwantyzacji sieci LUB wielokrotnych rozw. Phi0")
