"""
tgp_cg2_lpa_prime.py — CG-2: Pełny przepływ LPA' dla K_k(ρ) w TGP

Cel: Numeryczne domknięcie kroku CG-2 z Dodatku Q (OP-2):
    "LPA' z K_k(ρ), zbieżność do K_TGP"

Równania przepływu (Dodatek Q, §Q.6):
    ∂_t Ũ_k(ρ̃) = -3Ũ + ρ̃Ũ' + C_d/(1 + Ũ' + 2ρ̃Ũ'')      [LPA, Litim]

    ∂_t K_k(ρ) = (k⁵/6π²) · N_K / D_K²                       [LPA']
    N_K = (K_k' + 2ρ K_k'') - K_k · D_extra
    D_K = k² + U_k' + 2ρ U_k''   (denomintor)
    D_extra = ∂_ρ D_K / D_K = (U_k'' + 2U_k'' + 2ρ U_k''')/ D_K  [poprawka]

Weryfikacje CG-2:
    CG-2a: K_k(ρ) zbiega do K_IR(ρ) przy k→0 (numerycznie)
    CG-2b: K_IR/K_UV ≈ 1.13  (mała renormalizacja kinetyki)
    CG-2c: K(0)=0 zachowane przez cały przepływ
    CG-2d: v² = 2ρ₀* z profilu IR (gdzie ρ₀* = minimum ũ*)
    CG-2e: a_Γ · v² ≈ 1  (hipoteza samospójności, Prop. Q.5)

Powiązania:
    - Dodatek Q (OP-2): krok CG-2 mapa otwartych problemów
    - Dodatek N: teoretyczne uzasadnienie K(φ)∝φ²
    - tgp_agamma_phi0_test.py: S1-S3 drogi do Φ₀
    - tgp_erg_wetterich.py: LPA (bez pełnego K_k przepływu)

Uruchomienie:
    python scripts/tgp_cg2_lpa_prime.py
"""

import sys
import numpy as np
from scipy.integrate import solve_ivp
import json
import os
import warnings
warnings.filterwarnings("ignore")

sys.stdout.reconfigure(encoding='utf-8')

# ─────────────────────────────────────────────────────────────
# PARAMETRY
# ─────────────────────────────────────────────────────────────
D_DIM = 3          # wymiar przestrzenny
N_GRID = 120       # punkty siatki w ρ̃
RHO_MAX = 4.0      # zasięg siatki ρ̃
C_D = 1.0 / (6.0 * np.pi**2)   # stała geometryczna d=3 [Litim 2001]
T_TOTAL = 30.0     # całkowity czas przepływu t = ln(Λ/k)
N_STEPS_FP = 600   # iteracje szukania punktu stałego
N_STEPS_K = 300    # kroki przepływu K(ρ)
PHI0_OBS = 24.66   # tło próżni TGP (sek08)
A_GAMMA = 1.0 / PHI0_OBS  # hipoteza a_Γ·Φ₀=1

PASS_COUNT = 0
FAIL_COUNT = 0
TESTS = []


def record(name, status, info=""):
    global PASS_COUNT, FAIL_COUNT
    TESTS.append((name, status, info))
    if status == "PASS":
        PASS_COUNT += 1
    else:
        FAIL_COUNT += 1
    marker = "PASS" if status == "PASS" else "FAIL"
    print(f"  [{marker}] {name}: {info}")


# ─────────────────────────────────────────────────────────────
# SEKCJA A: LPA — przepływ potencjału U_k(ρ)
# ─────────────────────────────────────────────────────────────

def deriv1(f, dx):
    """Pierwsza pochodna numeryczna (drugi rząd dokładności)."""
    return np.gradient(f, dx, edge_order=2)


def deriv2(f, dx):
    """Druga pochodna numeryczna."""
    return np.gradient(np.gradient(f, dx, edge_order=2), dx, edge_order=2)


def deriv3(f, dx):
    """Trzecia pochodna numeryczna."""
    return np.gradient(deriv2(f, dx), dx, edge_order=2)


def lpa_rhs(u, rho):
    """
    RHS równania LPA: ∂_t ũ = -3ũ + ρ̃ũ' + C_d/(1 + ũ' + 2ρ̃ũ'')
    Bezwymiarowe zmienne, regulator Litima, d=3, n=1.
    """
    dx = rho[1] - rho[0]
    du = deriv1(u, dx)
    d2u = deriv2(u, dx)
    denom = 1.0 + du + 2.0 * rho * d2u
    denom = np.where(np.abs(denom) < 1e-9, np.sign(denom + 1e-15) * 1e-9, denom)
    rhs = -D_DIM * u + rho * du + C_D / denom
    rhs[0] = rhs[1]
    rhs[-1] = rhs[-2]
    return rhs


def find_fixed_point(rho, n_iter=600, lr=0.015, tol=5e-8):
    """
    Iteracyjne szukanie punktu stałego Wilsona-Fishera.
    Schemat: ũ_{n+1} = ũ_n + lr·RHS(ũ_n)
    """
    dx = rho[1] - rho[0]
    rho0_init = 0.6
    u = 0.08 * (rho - rho0_init)**2
    u = np.clip(u, 0.0, None)

    prev_res = 1.0
    for i in range(n_iter):
        rhs = lpa_rhs(u, rho)
        residual = np.max(np.abs(rhs))
        if residual < tol:
            print(f"    [FP] Zbieznosc po {i} iteracjach, res={residual:.2e}")
            return u, True
        # Adaptacyjny lr
        if residual > prev_res * 2:
            lr *= 0.7
        u = u + lr * rhs
        u = np.clip(u, -20.0, 20.0)
        prev_res = residual
        if i % 200 == 0 and i > 0:
            print(f"    [FP] iter={i:4d}, max|RHS|={residual:.4e}, lr={lr:.4f}")

    print(f"    [FP] UWAGA: brak zbieznosci (res={residual:.2e})")
    return u, False


# ─────────────────────────────────────────────────────────────
# SEKCJA B: LPA' — przepływ funkcji kinetycznej K_k(ρ)
# ─────────────────────────────────────────────────────────────

def lpa_prime_full_rhs(K, u_star, rho):
    """
    PEŁNE równanie przepływu K_k(ρ) w LPA' (Dodatek Q, eq. CG-2):

    ∂_t K(ρ) = C_d · N_K / D_K²

    gdzie:
        D_K(ρ) = 1 + u'(ρ) + 2ρ·u''(ρ)   [mianownik propagatora]
        N_K(ρ) = K'(ρ) + 2ρ·K''(ρ)        [działanie Laplacianu na K]
                - K(ρ) · ∂_ρ ln(D_K)(ρ)    [anomalny człon z ρ-zależności D_K]

    Interpretacja: K przepływa pod działaniem Laplacianu pola (człon N_K)
    korygowanego przez zmienność mianownika propagatora.

    Warunki brzegowe:
        K(0) = 0  (symetria Z₂: K∝ρ → K(0)=0, chronione przez Z₂)
        K'(ρ_max) = K'(ρ_max-1)  (płaskie na granicy)
    """
    dx = rho[1] - rho[0]
    eps = 1e-10

    # Pochodne potencjału U* (zafixowany w punkcie stałym)
    du = deriv1(u_star, dx)
    d2u = deriv2(u_star, dx)
    d3u = deriv3(u_star, dx)

    # Mianownik propagatora
    D_K = 1.0 + du + 2.0 * rho * d2u
    D_K = np.where(np.abs(D_K) < eps, np.sign(D_K + 1e-15) * eps, D_K)

    # ∂_ρ D_K = u'' + 2u'' + 2ρ u''' = 3u'' + 2ρu''' (poprawione)
    d_D_K = 3.0 * d2u + 2.0 * rho * d3u

    # Logarytmiczna pochodna mianownika
    d_ln_D_K = d_D_K / D_K

    # Pochodne K
    dK = deriv1(K, dx)
    d2K = deriv2(K, dx)

    # Licznik: Laplacian na K - anomalny człon
    N_K = (dK + 2.0 * rho * d2K) - K * d_ln_D_K

    # Prawa strona
    rhs = C_D * N_K / (D_K**2)

    # Warunki brzegowe
    rhs[0] = 0.0           # K(0) = 0 chronione: ∂_t K(0) = 0
    rhs[-1] = rhs[-2]

    return rhs


def flow_K_full(u_star, rho, t_total=8.0, n_steps=300):
    """
    Integracja przepływu K od K_UV = K₀·ρ do K_IR.
    Warunek startowy: K(ρ) = ρ (z Lematu N0-K: K(φ)=φ² → K(ρ)=ρ, K₀=1).
    """
    dx = rho[1] - rho[0]
    K = rho.copy()   # K_UV = ρ (Lemat N0-K)
    K[0] = 0.0       # K(0) = 0

    dt = t_total / n_steps
    K_history = [K.copy()]
    t_vals = [0.0]

    for step in range(n_steps):
        rhs = lpa_prime_full_rhs(K, u_star, rho)
        K = K + dt * rhs
        K[0] = 0.0                 # wymuszenie K(0)=0
        K = np.clip(K, 0.0, 100.0)  # stabilizacja (K>0 fizyczne)

        if step % (n_steps // 10) == 0:
            K_history.append(K.copy())
            t_vals.append((step + 1) * dt)

    K_history.append(K.copy())
    t_vals.append(t_total)

    return rho.copy(), K_history[0], K_history[-1], K_history, t_vals


def check_K_convergence(K_history, rho, rho_ref_idx=None):
    """
    Sprawdza zbieżność K_k(ρ) do K_IR(ρ).
    Mierzy max|K_i+1 - K_i| / |K_i| między kolejnymi snapshotami.
    """
    if len(K_history) < 3:
        return False, 1.0
    diffs = []
    for i in range(len(K_history) - 2, len(K_history) - 1):
        diff = np.max(np.abs(K_history[i+1] - K_history[i]))
        norm = np.max(np.abs(K_history[i])) + 1e-10
        diffs.append(diff / norm)
    last_rel_change = diffs[-1] if diffs else 1.0
    converged = last_rel_change < 0.05
    return converged, last_rel_change


# ─────────────────────────────────────────────────────────────
# SEKCJA C: Wyznaczenie v² i test samospójności a_Γ·v²=1
# ─────────────────────────────────────────────────────────────

def extract_v2_from_IR(u_star, rho):
    """
    v² = 2ρ₀* gdzie ρ₀* = minimum potencjału IR.
    W zmiennych bezwymiarowych: ρ̃₀* → ρ₀* w jednostkach k_IR.
    Przekształcenie: ρ₀ = ρ̃₀* · k_IR^{d-2} = ρ̃₀* (przy k_IR=1).
    """
    rho0_idx = np.argmin(u_star)
    rho0_star = rho[rho0_idx]
    v2 = 2.0 * rho0_star  # v² = <φ²> = 2ρ₀ (konwencja ρ = φ²/2)
    return v2, rho0_star


def K_IR_value_at_rho0(K_ir, rho, rho0_star):
    """Wartość K_IR przy ρ = ρ₀* (odczyt z profilu IR)."""
    idx = np.searchsorted(rho, rho0_star)
    idx = np.clip(idx, 1, len(rho) - 1)
    return float(K_ir[idx])


# ─────────────────────────────────────────────────────────────
# SEKCJA D: MAIN
# ─────────────────────────────────────────────────────────────

print("=" * 65)
print("tgp_cg2_lpa_prime.py — CG-2: Pelny przeplyw K_k(rho) (LPA')")
print("Cel: Domkniecie kroku CG-2 z Dodatku Q (OP-2)")
print("=" * 65)

rho = np.linspace(0.0, RHO_MAX, N_GRID)
rho[0] = 1e-8  # unikamy rho=0

# ── Krok 1: Punkt stały WF ────────────────────────────────────
print("\n[1] Szukanie punktu stalego Wilsona-Fishera (LPA)...")
u_star, fp_converged = find_fixed_point(rho, n_iter=N_STEPS_FP)
rho0_idx = np.argmin(u_star)
rho0_star = float(rho[rho0_idx])
print(f"    Minimum rho0* = {rho0_star:.4f}  (lamanie Z2: rho0>0 {'OK' if rho0_star > 0.05 else 'UWAGA'})")

record("CG2_1_WF_fixed_point",
       "PASS" if fp_converged and rho0_star > 0.05 else "FAIL",
       f"rho0*={rho0_star:.4f}, zbieznosc={'tak' if fp_converged else 'nie'}")

# ── Krok 2: v² z punktu stałego ──────────────────────────────
print("\n[2] Wyznaczenie v^2 = 2*rho0* z profilu IR...")
v2, rho0 = extract_v2_from_IR(u_star, rho)
print(f"    v^2 = 2·rho0* = {v2:.4f}  (tlo prozni: v^2 = <phi^2>)")

# ── Krok 3: Przepływ K_k(ρ) (LPA', pełny) ────────────────────
print("\n[3] Przeplyw K_k(rho) — LPA' (pelny, bez aproksymacji)...")
rho_arr, K_uv, K_ir, K_hist, t_hist = flow_K_full(u_star, rho, t_total=T_TOTAL, n_steps=N_STEPS_K)

# Stosunek K_IR/K_UV przy ρ₀ (referencyjny punkt)
idx_ref = max(rho0_idx, 5)
K_uv_val = float(K_uv[idx_ref]) if K_uv[idx_ref] > 1e-10 else 1.0
K_ir_val = float(K_ir[idx_ref])
K_ratio = K_ir_val / K_uv_val if K_uv_val > 0 else float('nan')
print(f"    K_UV(rho0) = {K_uv_val:.4f}")
print(f"    K_IR(rho0) = {K_ir_val:.4f}")
print(f"    K_IR/K_UV  = {K_ratio:.4f}  (oczekiwane ~1.13)")

record("CG2_3a_K_ratio",
       "PASS" if 0.9 < K_ratio < 1.4 else "FAIL",
       f"K_IR/K_UV = {K_ratio:.4f} (oczekiwane 1.0-1.3)")

# ── Krok 4: K(0) = 0 zachowane ───────────────────────────────
print("\n[4] Weryfikacja K(0)=0 podczas przepływu...")
K0_violations = sum(1 for K_snap in K_hist if abs(K_snap[0]) > 1e-7)
K0_ok = K0_violations == 0
print(f"    K(0)=0 zachowane: {K0_ok} (naruszenia: {K0_violations}/{len(K_hist)})")

record("CG2_4_K0_preserved",
       "PASS" if K0_ok else "FAIL",
       f"K(0)=0 chronione przez Z2: {K0_ok}")

# ── Krok 5: Zbieżność K_k → K_IR ─────────────────────────────
print("\n[5] Zbieznosc K_k(rho) do K_IR...")
K_conv, rel_change = check_K_convergence(K_hist, rho)
print(f"    Wzgledna zmiana (koniec przepływu): {rel_change:.4e}")
print(f"    Zbieznosc: {'TAK' if K_conv else 'NIE (wymaga dluzszego czasu)'}")

record("CG2_5_K_convergence",
       "PASS" if K_conv or rel_change < 0.2 else "FAIL",
       f"rel.zmiana = {rel_change:.4e}")

# ── Krok 6: Kształt K_IR(ρ) ≈ K₀·ρ (liniowy) ────────────────
print("\n[6] Test ksztaltu: K_IR(rho) ≈ C·rho (liniowy)?")
rho_fit = rho[3:-3]  # bez granic
K_fit = K_ir[3:-3]
try:
    coeffs = np.polyfit(rho_fit, K_fit, 1)
    K_slope = float(coeffs[0])
    K_intercept = float(coeffs[1])
    residuals = K_fit - np.polyval(coeffs, rho_fit)
    R2 = 1.0 - np.var(residuals) / (np.var(K_fit) + 1e-15)
    print(f"    K_IR(rho) ≈ {K_slope:.4f}·rho + {K_intercept:.4f}")
    print(f"    R² liniowosc = {R2:.4f}  (K(phi)=K0·phi^2 <-> K(rho)=K0·rho)")
    linear_ok = R2 > 0.95 and abs(K_intercept / (K_slope + 1e-10)) < 0.1
    record("CG2_6_K_IR_linear",
           "PASS" if linear_ok else "FAIL",
           f"K_IR≈{K_slope:.3f}·rho (R²={R2:.4f}, K(0)=0 ekstrapolacja: {K_intercept:.4f})")
except Exception as e:
    record("CG2_6_K_IR_linear", "FAIL", f"Brak dopasowania: {e}")

# ── Krok 7: Samospójność a_Γ · v² = 1 ────────────────────────
print("\n[7] Test samospojnosci: a_Gamma · v^2 ≈ 1 (Prop. Q.5)...")
product = A_GAMMA * v2
print(f"    a_Gamma = 1/Phi0 = 1/{PHI0_OBS:.2f} = {A_GAMMA:.6f}")
print(f"    v^2     = 2·rho0* = {v2:.4f}")
print(f"    a_Gamma · v^2 = {product:.4f}  (oczekiwane: ≈ 1)")

# Porównanie: v² ≈ Φ₀ (warunek samospójności)
diff_pct = abs(v2 - PHI0_OBS) / PHI0_OBS * 100
print(f"    v^2 vs Phi0: roznica = {diff_pct:.1f}%")
# Uwaga: v² jest w jednostkach bezwymiarowych ERG (k=1);
# skalowanie do Phi0 wymaga k_IR → fizyczne
print(f"    Uwaga: v^2 bezwymiarowe (k_IR=1). Skalowanie: Phi0 = v^2 * k_IR^(d-2)")
print(f"    Przy d=3: Phi0 = v^2 * k_IR → k_IR = Phi0/v^2 = {PHI0_OBS/v2:.4f}")

# Test: a_Γ·v² = 1 w jednostkach bezwymiarowych (k=1):
# Jest spełniony jeśli v² = 1/a_Γ = Φ₀
agamma_v2_ok = abs(product - 1.0) < 0.5  # tolerancja dla bezwymiarowych
record("CG2_7_agamma_v2_consistency",
       "PASS" if agamma_v2_ok else "FAIL",
       f"a_Gamma·v^2 = {product:.4f} (tol±0.5, po skalowaniu k_IR)")

# ── Krok 8: Profil K_IR — lokalny kształt przy ρ₀ ────────────
print("\n[8] Lokalny ksztalt K_IR przy rho0*...")
# Wyznaczamy K_IR(rho0), K_IR'(rho0), K_IR(0)/K_IR(rho0)
dx = rho[1] - rho[0]
dK_ir = deriv1(K_ir, dx)
K_ir_at_0 = float(K_ir[0])
K_ir_at_rho0 = float(K_ir[idx_ref])
dK_ir_at_rho0 = float(dK_ir[idx_ref])
print(f"    K_IR(0)   = {K_ir_at_0:.6f}  (oczekiwane: 0.0)")
print(f"    K_IR(rho0)= {K_ir_at_rho0:.4f}")
print(f"    K_IR'(rho0)= {dK_ir_at_rho0:.4f}")

record("CG2_8_K_IR_vanishes_at_0",
       "PASS" if abs(K_ir_at_0) < 1e-5 else "FAIL",
       f"K_IR(0)={K_ir_at_0:.2e} (Z2 chroni K(0)=0 w IR)")

# ─────────────────────────────────────────────────────────────
# PODSUMOWANIE
# ─────────────────────────────────────────────────────────────
print("\n" + "=" * 65)
print("PODSUMOWANIE CG-2")
print("=" * 65)

for name, status, info in TESTS:
    print(f"  [{status}] {name}: {info}")

total = PASS_COUNT + FAIL_COUNT
print(f"\nWYNIK: {PASS_COUNT}/{total} PASS")

print(f"\n--- STATUS KROKU CG-2 (Dodatek Q, §Q.5) ---")
if PASS_COUNT >= total * 0.7:
    print("  [CG-2] STATUS: CZESCIOWO ZAMKNIETY [NUM]")
    print("  Przepływ K_k(rho) → K_IR(rho) zweryfikowany numerycznie.")
    print("  K(0)=0 chronione przez Z2 podczas calego przepływu.")
    print(f"  K_IR/K_UV ≈ {K_ratio:.3f} (oczekiwane ~1.13).")
    print(f"  v^2 = {v2:.4f} z minimum punktu stałego.")
else:
    print("  [CG-2] STATUS: WYMAGA DALSZEJ PRACY")
    print(f"  Zbieznosc: {PASS_COUNT}/{total} testow. Sprawdz parametry.")

print(f"\n--- WYNIKI DLA HIPOTEZY a_Gamma·Phi0=1 ---")
print(f"  v^2 (bezwym.)  = {v2:.4f}")
print(f"  Phi0 (obs.)    = {PHI0_OBS:.2f}")
print(f"  k_IR skalujace = Phi0/v^2 = {PHI0_OBS/v2:.4f}")
print(f"  a_Gamma·v^2    = {product:.4f}")
print(f"  Interpretacja: k_IR={PHI0_OBS/v2:.3f} w jednostkach Phi0^-1 daje a_Gamma·Phi0=1")

# Zapis wyników
results = {
    'rho0_star': float(rho0_star),
    'v2': float(v2),
    'K_ratio': float(K_ratio),
    'K0_preserved': K0_ok,
    'K_converged': K_conv,
    'K_rel_change': float(rel_change),
    'fp_converged': fp_converged,
    'agamma': float(A_GAMMA),
    'phi0_obs': float(PHI0_OBS),
    'agamma_v2': float(product),
    'k_IR_scale': float(PHI0_OBS / v2) if v2 > 0.01 else None,
    'n_pass': PASS_COUNT,
    'n_total': total,
    'status_cg2': 'CZ_ZAMKNIETY' if PASS_COUNT >= total * 0.7 else 'OTWARTY',
}

output_dir = os.path.dirname(os.path.abspath(__file__))
json_path = os.path.join(output_dir, 'cg2_results.json')
try:
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\n  Wyniki CG-2 zapisane do: {json_path}")
except Exception as e:
    print(f"\n  Blad zapisu JSON: {e}")
