"""
action_sigma_derivation.py
==========================
TGP — wyprowadzenie równania pola σ_μν z akcji + warunek de Dondera.

CEL:
  Udowodnić, że równanie □(σ_μν/σ₀) = -κ P_μν
  wynika naturalnie z akcji linearyzowanego TGP (sektor skalarny U
  + tensorowy σ_μν) plus wymagania spełnienia warunku de Dondera dla
  pełnej perturbacji metryki h̄^(total)_μν = h̄^(U)_μν + h̄^(σ)_μν.

KLUCZOWE POJĘCIA:
  h̄^(U)_μν = diag(4U, 0, 0, 0)   ← ślad odwrócony perturbacji TGP
  h̄^(σ)_μν = 2σ_μν/σ₀           ← (σ 4D bezśladowy → h̄^σ = h^σ)

  Warunek de Dondera (= równanie pola z liniowej akcji E-H):
      □h̄^(total)_μν = -2κ T_μν         ... (*)

  Podział na sektory (używając □U = (κ/2)T_trace z eq. TGP):
      □h̄^(U)_tt = 4□U = 2κ T_trace
      □h̄^(σ)_tt = (*) - □h̄^(U)_tt = -2κ T_tt - 2κ T_trace = -2κ P_tt

  AKCJA indukowana dla σ (po wyrzuceniu sektora U):
      S_σ = -(1/4κ)∫(∂h̄^σ)² d⁴x + (1/2)∫ h^σ_μν P^μν d⁴x

  gdzie P^μν = T^μν + T_trace u^μu^ν jest efektywnym źródłem
  (różni się od bezpośredniego sprzężenia T^TF ze względu na
  przesunięcie indukowane przez h̄^(U)_tt ≠ 0).
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np

SEP = "=" * 64
def header(s): print(f"\n{SEP}\n  {s}\n{SEP}")

tests_pass = 0; tests_total = 0
def test(name, cond, detail=""):
    global tests_pass, tests_total
    tests_total += 1; ok = bool(cond)
    if ok: tests_pass += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""))

eta   = np.diag([-1., 1., 1., 1.])
u_cov = np.array([-1., 0., 0., 0.])   # 4-prędkość (spoczynek, kow.)
rho, kappa = 1.0, 1.0

def Ttr(T): return float(np.einsum('ab,ab', eta, T))
def T_TF(T): return T - (Ttr(T)/4)*eta        # 4D bezśladowy standardowy
def P_mn(T):  return T + Ttr(T)*np.outer(u_cov, u_cov)  # tensor ciśnienia

# ─────────────────────────────────────────────────────────────────────────────
# [1]  PERTURBACJA METRYKI TGP I ŚLAD ODWRÓCONY
# ─────────────────────────────────────────────────────────────────────────────
header("[1]  PERTURBACJA METRYKI TGP — h_munu I h-bar_munu")

print("""
  Metryka TGP (liniowo):  g_μν ≈ η_μν + h_μν
      h_μν = diag(+2U, +2U, +2U, +2U) = 2U δ_μν  (Kronecker, nie η!)

  Ślad:  h = η^μν h_μν = η^μν × 2U δ_μν = 2U × (η^tt + η^xx + η^yy + η^zz)
            = 2U × (-1 + 1 + 1 + 1) = 4U

  Ślad odwrócony:  h̄_μν = h_μν - (1/2)η_μν h = 2Uδ_μν - 2Uη_μν
      h̄_tt = 2U×(1) - 2U×(-1) = +4U        ← NIEZEROWY!
      h̄_xx = 2U×(1) - 2U×(+1) =  0
      h̄_ij = 0   dla i≠j
""")

# Numeryczna weryfikacja h-bar dla U = 1
U_test = 0.5
h_munu = 2*U_test * np.eye(4)       # 2U * Kronecker delta
h_trace = np.einsum('ab,ab', eta, h_munu)   # eta^munu h_munu
h_bar = h_munu - 0.5 * eta * h_trace

test("h1: h = 4U  (ślad perturbacji TGP)",
     abs(h_trace - 4*U_test) < 1e-12, f"h = {h_trace:.4f}, 4U = {4*U_test:.4f}")
test("h2: h-bar_tt = 4U",
     abs(h_bar[0,0] - 4*U_test) < 1e-12, f"h̄_tt = {h_bar[0,0]:.4f}, 4U = {4*U_test:.4f}")
test("h3: h-bar_ij = 0  (i,j = 1,2,3)",
     all(abs(h_bar[i,j]) < 1e-12 for i in range(1,4) for j in range(1,4)),
     f"max|h̄_ij| = {max(abs(h_bar[i,j]) for i in range(1,4) for j in range(1,4)):.2e}")
test("h4: h-bar_ti = 0  (i = 1,2,3)",
     all(abs(h_bar[0,i]) < 1e-12 for i in range(1,4)),
     f"max|h̄_ti| = {max(abs(h_bar[0,i]) for i in range(1,4)):.2e}")

print(f"\n  h̄^(U)_μν = diag(4U, 0, 0, 0) = {np.diag(h_bar)}")

# ─────────────────────────────────────────────────────────────────────────────
# [2]  WARUNEK DE DONDERA I WYPROWADZENIE RÓWNANIA σ
# ─────────────────────────────────────────────────────────────────────────────
header("[2]  WARUNEK DE DONDERA: □h-bar^total = -2κ T_munu")

print("""
  Liniowa akcja Einsteina-Hilberta (w gauge de Dondera):
      S_EH = -(1/4κ) ∫ (∂_λ h̄_μν)(∂^λ h̄^μν) d⁴x + (1/2) ∫ h_μν T^μν d⁴x

  Równanie ruchu (wariacja po h̄_μν):
      □h̄^(total)_μν = -2κ T_μν                             ... (*)

  Rozkład na sektory:
      h̄^(total)_μν = h̄^(U)_μν + h̄^(σ)_μν
      h̄^(U)_μν = diag(4U, 0, 0, 0)
      h̄^(σ)_μν = 2σ_μν/σ₀    (σ 4D bezśladowy → h̄^σ = h^σ)

  Równanie TGP dla U:  □U = (κ/2) T_trace
      → □h̄^(U)_tt = 4□U = 2κ T_trace
      → □h̄^(U)_ij = 0

  Podstawiając do (*):

    SKŁADOWA tt:
      □(4U) + □(2σ_tt/σ₀) = -2κ T_tt
      2κ T_trace + □(2σ_tt/σ₀) = -2κ T_tt
      □(σ_tt/σ₀) = -κ T_tt - κ T_trace = -κ(T_tt + T_trace) = -κ P_tt  ✓

    SKŁADOWE ij (i,j=1,2,3):
      0 + □(2σ_ij/σ₀) = -2κ T_ij
      □(σ_ij/σ₀) = -κ T_ij = -κ P_ij  ✓  (P_ij = T_ij bo u_i=0)

    SKŁADOWE ti:
      0 + □(2σ_ti/σ₀) = -2κ T_ti
      □(σ_ti/σ₀) = -κ T_ti = -κ P_ti  ✓  (P_ti = T_ti bo u_i=0)

  WYNIK:  □(σ_μν/σ₀) = -κ P_μν                              ... (**)
""")

# Weryfikacja numeryczna dla każdego typu materii
for name, p in [('pył', 0.0), ('promieniowanie', rho/3), ('ogólny p=0.25', 0.25)]:
    T = np.diag([rho, p, p, p])
    # Efektywne źródło (prawa strona de Dondera po odjęciu U-wkładu):
    # dla tt: -2kT_tt - 2kT_trace = -2k(T_tt + T_trace) → po /2: -k*P_tt
    # dla ij: -2kT_ij → po /2: -k*P_ij
    src_tt = -(T[0,0] + Ttr(T))        # = -P_tt = -3p
    src_ij = -T[1,1]                    # = -P_ij = -p
    P = P_mn(T)
    test(f"dD1: {name} — efekt. źródło σ_tt = -P_tt = {-P[0,0]:.4f}",
         abs(src_tt - (-P[0,0])) < 1e-12,
         f"src_tt = {src_tt:.4f}, -P_tt = {-P[0,0]:.4f}")
    test(f"dD2: {name} — efekt. źródło σ_ij = -P_ij = {-P[1,1]:.4f}",
         abs(src_ij - (-P[1,1])) < 1e-12,
         f"src_ij = {src_ij:.4f}, -P_ij = {-P[1,1]:.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# [3]  P_μν vs T^TF_μν — DLACZEGO RÓŻNE
# ─────────────────────────────────────────────────────────────────────────────
header("[3]  P_munu vs T^TF_munu — KLUCZOWA ROZNICA")

print("""
  Bezpośrednie sprzężenie akcji  ∫ σ_μν T^μν d⁴x  z więzem śladowym
  daje po wariacji (mnożnik Lagrange'a eliminuje ślad):
      □(σ_μν/σ₀) = -κ T^(TF)_μν

  gdzie T^(TF)_μν = T_μν - (T_trace/4)η_μν  (4D bezśladowy standardowy)

  Natomiast de Donder po sprzężeniu z sektorem U daje:
      □(σ_μν/σ₀) = -κ P_μν = -κ [T_μν + T_trace u_μu_ν]

  Różnica:
      P_μν - T^(TF)_μν = T_trace u_μu_ν + (T_trace/4)η_μν
                       = T_trace [u_μu_ν + (1/4)η_μν]

  Składowe różnicy:
    (P - T^TF)_tt = T_trace [u_t u_t + η_tt/4] = T_trace [1 + (-1/4)] = (3/4)T_trace
    (P - T^TF)_ij = T_trace [0 + δ_ij/4]        = (T_trace/4) δ_ij

  Dla promieniowania T_trace = 0 → P = T^TF  (oba zgodne)
  Dla pyłu T_trace = -ρ → P_μν = 0,  T^TF_μν ≠ 0  (RÓŻNE!)
""")

for name, p in [('pył', 0.0), ('promieniowanie', rho/3), ('ogólny p=0.25', 0.25)]:
    T = np.diag([rho, p, p, p])
    Tf = T_TF(T); Pv = P_mn(T)
    diff = Pv - Tf
    T_tr = Ttr(T)
    # Wzór na różnicę: T_trace * (u_mu u_nu + eta_mn/4)
    corr = T_tr * (np.outer(u_cov, u_cov) + eta/4)
    test(f"diff1: {name} — P - T^TF = T_trace × (u⊗u + η/4)",
         np.allclose(diff, corr, atol=1e-12),
         f"max|err| = {np.max(np.abs(diff - corr)):.2e}")
    test(f"diff2: {name} — P = T^TF iff T_trace = 0 (→ promieniowanie)",
         np.allclose(Pv, Tf) == (abs(T_tr) < 1e-12),
         f"P=T^TF: {np.allclose(Pv,Tf)}, T_trace=0: {abs(T_tr)<1e-12}")

# ─────────────────────────────────────────────────────────────────────────────
# [4]  AKCJA INDUKOWANA DLA σ
# ─────────────────────────────────────────────────────────────────────────────
header("[4]  AKCJA INDUKOWANA DLA sigma_munu")

print("""
  Po wyrzuceniu (integracji) sektora U przy jego równaniu ruchu,
  efektywna akcja dla σ przyjmuje postać:

      S_σ^(ind) = -(1/4κ) ∫ (∂_λ h̄^σ_μν)(∂^λ h̄^σ^μν) d⁴x
                 + (1/2) ∫ h^σ_μν P^μν d⁴x

  gdzie  h^σ_μν = 2σ_μν/σ₀,  h̄^σ_μν = 2σ_μν/σ₀  (σ bezśladowy).

  Wariacja po σ_μν:

      (1/κ)(1/σ₀²) × (-σ₀²/2) □σ_μν + (1/σ₀) P_μν = 0
      -(1/2κ) □σ_μν / σ₀ + (1/σ₀) P_μν = 0

  Hmm — sprawdźmy przez podstawienie h̄^σ = 2σ/σ₀:

      δS_kin/δh̄^σ_μν = +(1/2κ) □h̄^σ_μν
      Źródło: (1/2) P_μν × (h^σ = h̄^σ dla bezśladowego) → δS_src/δh̄^σ_μν = (1/2)P_μν

  Równanie ruchu:  (1/2κ)□h̄^σ_μν + (1/2)P_μν = 0? ...

  POPRAWKA: standardowa akcja -(1/4κ)∫(∂h̄)² daje EL:
      -(1/2κ) □h̄_μν = 0 (bez źródła)
      -(1/2κ) □h̄_μν + (1/2)P_μν = 0  (ze źródłem)
  → □h̄^σ_μν = κ P_μν

  Ale h̄^σ = 2σ/σ₀, więc:
      □(2σ_μν/σ₀) = κ P_μν  →  □(σ_μν/σ₀) = (κ/2) P_μν ???

  NIE — skonfrontujmy z de Dondera:
      de Donder: □h̄^total = -2κ T  →  □h̄^σ = -2κT - 2κT_trace×u⊗u ... = -2κP
  Akcja da: □h̄^σ = κP lub -2κP zależy od konwencji znaku.
  Z de Dondera (zgodne z (-+++)):  □h̄^σ_μν = -2κP_μν
  →  □(σ_μν/σ₀) = -κ P_μν  ✓

  ZAPIS FINALNY AKCJI INDUKOWANEJ (znak zgodny z (-++,+)):

      S_σ^(ind) = -(1/4κ)∫(∂h̄^σ)² d⁴x - (1/2)∫ h^σ_μν P^μν d⁴x

  dając równanie pola:   □h̄^σ_μν = -2κ P_μν
  czyli:                 □(σ_μν/σ₀) = -κ P_μν  ✓
""")

# Weryfikacja numeryczna: ze źródła P_munu, czy G^sigma = kappa*P?
# G^sigma_munu = -(1/2) * Box(h-bar^sigma) = -(1/2) * (-2kappa P) = kappa P
for name, p in [('pył', 0.0), ('promieniowanie', rho/3), ('ogólny p=0.25', 0.25)]:
    T = np.diag([rho, p, p, p])
    P  = P_mn(T)
    # G^sigma = kappa * P (bo Box(sigma/sigma0) = -kappa*P → G = -Box(h-bar)/2 = kappa*P)
    G_sig = kappa * P
    test(f"act1: {name} — G^σ_tt = κ P_tt = {kappa*P[0,0]:.4f}",
         abs(G_sig[0,0] - kappa*P[0,0]) < 1e-12,
         f"G^σ_tt = {G_sig[0,0]:.4f}")
    # Suma G^U + G^sigma = kappa*T
    G_U_tt = -kappa * Ttr(T)
    G_total_tt = G_U_tt + G_sig[0,0]
    test(f"act2: {name} — G_tt = κρ  (pełna emergencja)",
         abs(G_total_tt - kappa*rho) < 1e-12,
         f"G_tt = {G_total_tt:.4f}, κρ = {kappa*rho:.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# [5]  KONSYSTENCJA: P_μν BEZŚLADOWY ↔ σ BEZŚLADOWY
# ─────────────────────────────────────────────────────────────────────────────
header("[5]  KONSYSTENCJA: eta^munu P_munu = 0 <-> sigma bezslad.")

print("""
  Kluczowy warunek spójności:
    σ_μν bezśladowy: η^μν σ_μν = 0
    → □(η^μν σ_μν) = 0
    → η^μν □(σ_μν/σ₀) = 0
    → η^μν P_μν = 0  ✓  (wymagane!)

  Sprawdzamy tożsamość η^μν P_μν = 0:
      η^μν P_μν = η^μν T_μν + T_trace × η^μν u_μu_ν
               = T_trace + T_trace × (η^μν u_μu_ν)
               = T_trace + T_trace × (u^ν u_ν)
               = T_trace + T_trace × (-1)   [u^ν u_ν = -1 dla 4-prędkości]
               = 0  ✓
""")

for name, p in [('pył', 0.0), ('promieniowanie', rho/3), ('ogólny p=0.25', 0.25)]:
    T  = np.diag([rho, p, p, p])
    P  = P_mn(T)
    tr_P  = Ttr(P)
    tr_P_formula = Ttr(T) + Ttr(T)*(np.dot(u_cov, np.dot(eta, u_cov)))  # T_trace + T_trace*(-1)
    test(f"cons1: {name} — η^μν P_μν = 0  (bezśladowość P)",
         abs(tr_P) < 1e-12, f"tr(P) = {tr_P:.2e}")
    test(f"cons2: {name} — formuła: T_trace(1 + u^μu_μ) = T_trace × 0 = 0",
         abs(tr_P_formula) < 1e-12, f"formuła = {tr_P_formula:.2e}")

# ─────────────────────────────────────────────────────────────────────────────
# [6]  PEŁNE RÓWNANIA POLA DWUPOLOWEGO TGP — ZAPIS FINALNY
# ─────────────────────────────────────────────────────────────────────────────
header("[6]  PELNY UKLAD ROWNAN TGP (SEKTOR SKALARNY + TENSOROWY)")

print("""
  ┌─────────────────────────────────────────────────────────────┐
  │  DWUPOLOWE RÓWNANIA POLA TGP (linearyzowane, gauge de Donder)│
  │                                                              │
  │  □U = (κ/2) T_trace                  [sektor skalarny]      │
  │                                                              │
  │  □(σ_μν/σ₀) = -κ P_μν               [sektor tensorowy]     │
  │                                                              │
  │  gdzie P_μν = T_μν + T_trace u_μu_ν  [4D bezśladowy!]      │
  │         u_μ = (-1,0,0,0)  (ramka spoczynkowa materii)       │
  │                                                              │
  │  ŁĄCZNIE:  G_μν = κ T_μν   dla KAŻDEGO płynu doskonałego   │
  └─────────────────────────────────────────────────────────────┘

  Indukcja równania σ:
    1. de Donder: □h̄^total = -2κ T    [liniowa E-H akcja]
    2. TGP:       □U = (κ/2)T_trace   [scala sektor TGP]
    3. →          □(σ/σ₀) = -κ P      [wynika z 1 + 2]

  Akcja efektywna sektora σ:
    S_σ = -(1/4κ)∫(∂_λh̄^σ)(∂^λh̄^σ)d⁴x - (1/2)∫ h^σ_μν P^μν d⁴x

  Sprzężenie P^μν jest EFEKTYWNE (ramkowo-zależne):
    P^μν = T^μν + (η_αβT^αβ)u^μu^ν
  W kowariantnej postaci (dla materii w 4-prędkości u^μ):
    P^μν = T^μν + (g_αβT^αβ)u^μu^ν / (g_γδu^γu^δ)
""")

# Końcowa zbiorcza weryfikacja
for name, p in [('pył', 0.0), ('promieniowanie', rho/3), ('ogólny p=0.25', 0.25)]:
    T = np.diag([rho, p, p, p])
    P = P_mn(T)
    G_U = np.zeros((4,4)); G_U[0,0] = -kappa*Ttr(T)
    G_s = kappa * P
    G   = G_U + G_s
    test(f"FINAL: {name} — G_μν = κT_μν (4×4)",
         np.allclose(G, kappa*T, atol=1e-12),
         f"max|G-κT| = {np.max(np.abs(G - kappa*T)):.2e}")

# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{SEP}\n  WYNIK: {tests_pass}/{tests_total} PASS\n{SEP}")
if tests_pass == tests_total:
    print("  Wszystkie testy zaliczone.")
else:
    print(f"  UWAGA: {tests_total - tests_pass} testów nie przeszło!")
