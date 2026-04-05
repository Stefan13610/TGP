"""
v3_diagnoza.py - Diagnoza podejscia v3 (regularyzowany kwartyczny)

Pytania:
1. Czy g(K) = E(K)/K - 4pi w ogole wraca do zera dla duzych K?
2. Jak zachowuje sie Ep(K) przy wzroscie K?
3. Porownanie: v2 (lambda*(psi-1)^6) vs v3 (regularyzacja)
"""
import numpy as np
import warnings; warnings.filterwarnings('ignore')

GAMMA = 1.0
ALPHA = 7.0
A_GAM = 0.030

def V_v2(psi, lam=5.514e-7):
    """V_mod z lambda*(psi-1)^6 (wersja 2)."""
    return GAMMA/3*psi**3 - GAMMA/4*psi**4 + lam/6*(psi-1)**6

def V_v3(psi, eps=1e-8, n=3):
    """V_mod = psi^3/3 - psi^4 / (4*(1+eps*(psi-1)^n)) (wersja 3)."""
    return GAMMA/3 * psi**3 - GAMMA/4 * psi**4 / (1 + eps*(psi-1)**n)

def energy_generic(K, alpha, a_gam, V_func, V1_val, N=3000):
    msp   = 1.0
    r_max = max(60.0/msp, 20.0)
    r     = np.linspace(a_gam, r_max, N)
    phi   = np.maximum(1.0 + K * np.exp(-msp*r)/r, 1e-10)
    dphi  = K * np.exp(-msp*r)*(-msp*r - 1.0)/r**2
    Ek    = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
    psi   = phi
    Ep    = 4*np.pi*np.trapezoid((V_func(psi) - V1_val)*r**2, r)
    return Ek, Ep, Ek+Ep

V1_v2 = V_v2(1.0)
V1_v3 = V_v3(1.0)

# K-siatka dla diagnostyki
K_arr = np.concatenate([
    np.linspace(0.01, 1.0,  50),
    np.linspace(1.0,  10.0, 50),
    np.linspace(10.0, 100.0, 50),
    np.linspace(100.0, 500.0, 30),
])

print("=" * 70)
print("DIAGNOSTYKA: g(K) = E(K)/K/4pi - 1 dla duzych K")
print(f"alpha={ALPHA}, a_gam={A_GAM}")
print("=" * 70)
print()

# ── Analiza behawioralna: Ep(K) przy duzych K ──────────────────────────────
print("Zachowanie Ep(K) dla duzych K:")
print(f"{'K':>8} {'Ep_v2':>14} {'Ep_v3(1e-8)':>14} {'Ep_v3(1e-6)':>14}")
print("-" * 55)
for K in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0]:
    _, Ep2, _ = energy_generic(K, ALPHA, A_GAM, V_v2, V1_v2)
    _, Ep3a,_ = energy_generic(K, ALPHA, A_GAM,
                               lambda p: V_v3(p, eps=1e-8), V1_v3)
    _, Ep3b,_ = energy_generic(K, ALPHA, A_GAM,
                               lambda p: V_v3(p, eps=1e-6), V1_v3)
    print(f"{K:>8.2f} {Ep2:>14.4e} {Ep3a:>14.4e} {Ep3b:>14.4e}")

print()
print("Wniosek: jesli Ep_v3 stale maleje (coraz bardziej ujemne) => brak 3. zera")
print("         jesli Ep_v2 osiaga minimum i wraca => mozliwe 3. zero")

print()
print("=" * 70)
print("g(K) = E(K)/K - 4*pi dla v2 (K do 200)")
print(f"{'K':>8} {'g_v2':>12} {'Ek_v2':>12} {'Ep_v2':>12}")
print("-" * 50)
for K in [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0]:
    Ek, Ep, E = energy_generic(K, ALPHA, A_GAM, V_v2, V1_v2)
    g = E/K - 4*np.pi
    print(f"{K:>8.2f} {g:>12.4e} {Ek:>12.4e} {Ep:>12.4e}")

print()
print("=" * 70)
print("g(K) = E(K)/K - 4*pi dla v3, eps=1e-8 (K do 500)")
print(f"{'K':>8} {'g_v3':>12} {'Ek_v3':>12} {'Ep_v3':>12} {'psi_core':>10}")
print("-" * 60)
for K in [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0,
          100.0, 150.0, 200.0, 300.0, 400.0, 500.0]:
    Ek, Ep, E = energy_generic(K, ALPHA, A_GAM,
                               lambda p: V_v3(p, eps=1e-8), V1_v3)
    g   = E/K - 4*np.pi
    psi_c = 1.0 + K * np.exp(-1.0*A_GAM) / A_GAM
    print(f"{K:>8.2f} {g:>12.4e} {Ek:>12.4e} {Ep:>12.4e} {psi_c:>10.1f}")

print()
print("=" * 70)
print("ANALIZA POTENCJALU: V_mod(psi) - V_mod(1) przy roznych psi")
print("Kluczowe: kiedy V_mod wraca powyzej V_mod(1)?")
print()
print(f"{'psi':>10} {'V_v2-V1':>14} {'V_v3(1e-8)-V1':>16} {'V_orig-V1':>14}")
print("-" * 60)
for psi in [1.0, 1.5, 2.0, 5.0, 10.0, 50.0, 100.0, 500.0,
            1000.0, 5000.0, 10000.0]:
    Vo = GAMMA/3*psi**3 - GAMMA/4*psi**4 - (GAMMA/3 - GAMMA/4)
    V2 = V_v2(psi) - V1_v2
    V3 = V_v3(psi, eps=1e-8) - V1_v3
    print(f"{psi:>10.1f} {V2:>14.4e} {V3:>16.4e} {Vo:>14.4e}")

print()
print("Dla v2: V_mod powraca powyzej V_mod(1) dla psi > psi_positive")
print("Dla v3: V_mod NIGDY nie wraca (rosnie jak psi^3/3 - psi/(4eps))")
print("  => v3 NIE tworzy 3. zera g(K)!")

print()
print("=" * 70)
print("WNIOSEK FUNDAMENTALNY:")
print()
print("  v2 (lambda*(psi-1)^6): Ep(K) ma lokalne minimum, potem ROSNIE")
print("  => E(K)/K wraca do 4*pi => 3. zero (ale artefakt numeryczny)")
print()
print("  v3 (regularyzacja kwartyczna): Ep(K) STALE maleje (coraz bardziej ujemne)")
print("  => E(K)/K nigdy nie wraca do 4*pi => brak 3. zera")
print()
print("  DLACZEGO?")
print("  v3: V_mod(psi) = psi^3/3 - psi/(4*eps) dla duzych psi")
print("  V_mod - V_mod(1) = psi^3/3 - psi/(4*eps) - 1/12")
print("  To jest UJEMNE dla psi < sqrt(3/(4*eps)) = psi_pos")
print("  i DODATNIE tylko dla psi > psi_pos (bardzo duze!)")
print()
print("  Calk Ep = [dodatnie od korzenia] + [ujemne od powloki]")
print("  Powloka jest o wiele wieksza objetosciowo => Ep zawsze ujemne")
print()
psi_pos = np.sqrt(3/(4*1e-8))
print(f"  Dla eps=1e-8: psi_pos = sqrt(3/4e-8) = {psi_pos:.0f}")
print(f"  K3 ~ psi_pos * a_gam = {psi_pos*A_GAM:.0f} (nieosiagalny zasiag!)")
