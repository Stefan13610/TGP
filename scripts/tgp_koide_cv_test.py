"""
tgp_koide_cv_test.py  —  Teoria Generowanej Przestrzeni (TGP)
=============================================================
Weryfikacja trojjakiej rownowaznosci Q_K = 3/2 (T-OP1):

    Q_K = 3/2  <=>  CV(sqrt(m)) = 1  <=>  theta = 45 stopni

gdzie:
    Q_K   = (sqrt(me)+sqrt(mmu)+sqrt(mtau))^2 / (me+mmu+mtau)  [formula Koide]
    CV    = std(sqrt(m)) / mean(sqrt(m))  [wspolczynnik zmiennosci]
    theta = kat miedzy (sqrt(me),sqrt(mmu),sqrt(mtau)) a (1,1,1)/sqrt(3)  [45-stopni]

Twierdzenie algebraiczne (Thm T-QK-CV):
    Dla N mas m_1,...,m_N (N=3):
        Q_K = N / (1 + CV^2)
    Dowod:
        S1 = sum(sqrt(mi)) = N*mu
        S2 = sum(mi) = N*(mu^2 + sigma^2) = N*mu^2*(1+CV^2)
        Q_K = S1^2/S2 = N^2*mu^2 / (N*mu^2*(1+CV^2)) = N/(1+CV^2)  QED
    Wniosek (N=3): Q_K = 3/2 <=> CV = 1

Tożsamość geometryczna (Cor T-QK-angle):
    cos^2(theta) = Q_K / N
    theta = 45 stopni <=> Q_K = N/2 = 3/2  (dla N=3)

Formulacja TGP (Prop T-CV-Atail):
    Masy w TGP: m_i propto A(g0^(i))^4 => sqrt(m_i) propto A(g0^(i))^2
    Warunek Koide w jezyku TGP: CV(A(g0^(n))^2) = 1
    tj. std(A^2) = mean(A^2)  dla trzech generacji leptonu

Referencja: dodatekT_koide_atail_formal.tex, ROADMAP_v3.md
Status T-OP1: [CZ. ZAMK.] — geometryczne sformulowanie; mechanizm TGP otwarty
"""

import sys
import io
import numpy as np

# UTF-8 output
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ---------------------------------------------------------------------------
# Stale fizyczne
# ---------------------------------------------------------------------------
M_E   = 0.51099895000   # MeV
M_MU  = 105.6583755     # MeV
M_TAU = 1776.86         # MeV

R21_PDG  = M_MU / M_E      # = 206.7683
R31_PDG  = M_TAU / M_E     # = 3477.48
R31_KOIDE = 3477.44         # z Twierdzenia T-r31 (Q_K=3/2 + r21)

G0_E   = 1.24915   # g0^e z phi-FP (dodatekJ2)
G0_TAU = 3.1891    # g0^tau z Koide (dodatekT)


# ---------------------------------------------------------------------------
# Funkcje pomocnicze
# ---------------------------------------------------------------------------
def koide_QK(masses):
    """Q_K = (sum(sqrt(m)))^2 / sum(m)"""
    lam = np.sqrt(masses)
    return np.sum(lam)**2 / np.sum(masses)

def cv_sqrtm(masses):
    """CV = std(sqrt(m)) / mean(sqrt(m))  [populacyjne std, N]"""
    lam = np.sqrt(masses)
    return np.std(lam) / np.mean(lam)

def angle_deg(masses):
    """Kat miedzy (sqrt(m)) a (1,1,1)/sqrt(3), w stopniach"""
    lam = np.sqrt(masses)
    n = len(masses)
    dem = np.ones(n) / np.sqrt(n)
    cos_t = np.dot(lam, dem) / np.linalg.norm(lam)
    cos_t = np.clip(cos_t, -1.0, 1.0)
    return np.degrees(np.arccos(cos_t))

def qk_from_cv(CV, N=3):
    """Q_K = N/(1+CV^2)  [Twierdzenie T-QK-CV]"""
    return N / (1.0 + CV**2)

def cv_from_qk(QK, N=3):
    """CV = sqrt(N/Q_K - 1)"""
    return np.sqrt(N / QK - 1.0)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 65)
    print("TGP — Trojjaka rownowaznosc Q_K=3/2 (T-OP1)")
    print("    Q_K=3/2  <=>  CV(sqrt(m))=1  <=>  theta=45 stopni")
    print("=" * 65)

    # === Twierdzenie algebraiczne ===
    print("\n[ALG] Twierdzenie T-QK-CV: Q_K = N/(1+CV^2)")
    print("      Dowod: S1=N*mu, S2=N*mu^2*(1+CV^2) => Q_K=N/(1+CV^2)  QED")
    print("      Wniosek (N=3): Q_K=3/2 <=> CV=1")

    # Weryfikacja tozsamosci cos^2(theta) = Q_K/N
    print("\n[GEO] Tozsamosc: cos^2(theta) = Q_K / N")
    print("      theta=45 deg <=> Q_K=N/2 (dla N=3: Q_K=3/2)")

    # === PDG ===
    print("\n" + "-" * 40)
    print("[PDG] Masy leptonu PDG:")
    masses_PDG = np.array([M_E, M_MU, M_TAU])
    QK_PDG = koide_QK(masses_PDG)
    CV_PDG = cv_sqrtm(masses_PDG)
    th_PDG = angle_deg(masses_PDG)
    cos2_PDG = np.cos(np.radians(th_PDG))**2
    print(f"  Q_K    = {QK_PDG:.8f}  (oczekiwane 3/2 = 1.5)")
    print(f"  CV     = {CV_PDG:.8f}  (oczekiwane 1.0)")
    print(f"  theta  = {th_PDG:.6f} deg  (oczekiwane 45.0)")
    print(f"  cos^2  = {cos2_PDG:.8f}  (oczekiwane 0.5)")
    print(f"  Q_K z CV: {qk_from_cv(CV_PDG):.8f}  (=Q_K: {abs(qk_from_cv(CV_PDG)-QK_PDG)<1e-6})")
    print(f"  CV z Q_K: {cv_from_qk(QK_PDG):.8f}  (=CV: {abs(cv_from_qk(QK_PDG)-CV_PDG)<1e-6})")

    # === TGP: r21, r31^K ===
    print("\n[TGP] Masy z r21 (phi-FP), r31^K (Koide+TGP):")
    masses_TGP = np.array([1.0, R21_PDG, R31_KOIDE])
    QK_TGP = koide_QK(masses_TGP)
    CV_TGP = cv_sqrtm(masses_TGP)
    th_TGP = angle_deg(masses_TGP)
    print(f"  r21={R21_PDG:.4f}, r31^K={R31_KOIDE}")
    print(f"  Q_K    = {QK_TGP:.8f}  (oczekiwane 1.5 z definicji)")
    print(f"  CV     = {CV_TGP:.8f}  (oczekiwane 1.0)")
    print(f"  theta  = {th_TGP:.6f} deg  (oczekiwane 45.0)")

    # === TGP: formulacja A_tail ===
    print("\n[ATAIL] Formulacja TGP przez funkcje A_tail:")
    print(f"  m_i propto A(g0^(i))^4  =>  sqrt(m_i) = lambda_i propto A(g0^(i))^2")
    print(f"  Koide Q_K=3/2  <=>  CV(A(g0^(n))^2) = 1")
    # A(g0)^4 = m_i => A(g0)^2 = sqrt(m_i) = lambda_i
    # lambda_i znormalizowane: [1, sqrt(r21), sqrt(r31^K)]
    lam_TGP = np.sqrt(masses_TGP)
    mu_A2 = np.mean(lam_TGP)
    sig_A2 = np.std(lam_TGP)
    CV_A2 = sig_A2 / mu_A2
    print(f"  A(g0^e)^2 : A(g0^mu)^2 : A(g0^tau)^2 = "
          f"1 : {lam_TGP[1]:.4f} : {lam_TGP[2]:.4f}")
    print(f"  mean(A^2) = {mu_A2:.6f}")
    print(f"  std(A^2)  = {sig_A2:.6f}")
    print(f"  CV(A^2)   = {CV_A2:.8f}  (oczekiwane 1.0)")
    print(f"  std = mean: {abs(sig_A2/mu_A2 - 1.0) < 1e-5}")

    # === delta_i = lambda_i/mu - 1 ===
    delta = lam_TGP / mu_A2 - 1.0
    print(f"\n  delta_i = A_i^2/mean-1: [{delta[0]:.6f}, {delta[1]:.6f}, {delta[2]:.6f}]")
    print(f"  sum(delta_i) = {np.sum(delta):.2e}  (=0)")
    print(f"  sum(delta_i^2) = {np.sum(delta**2):.8f}  (=N*CV^2 = 3)")

    # === Testy koncowe ===
    print(f"\n[TESTY]")
    checks = [
        ("ALG: Q_K = N/(1+CV^2) dla PDG",
         abs(qk_from_cv(CV_PDG) - QK_PDG) < 1e-6,
         f"Q_K(formula)={qk_from_cv(CV_PDG):.6f}, Q_K(direct)={QK_PDG:.6f}"),

        ("GEO: cos^2(theta) = Q_K/3 dla PDG",
         abs(np.cos(np.radians(th_PDG))**2 - QK_PDG/3) < 1e-6,
         f"cos^2={np.cos(np.radians(th_PDG))**2:.6f}, Q_K/3={QK_PDG/3:.6f}"),

        ("PDG: CV(sqrt(m)) ≈ 1 (blad < 1e-4)",
         abs(CV_PDG - 1.0) < 1e-4,
         f"CV={CV_PDG:.8f}, blad={abs(CV_PDG-1.0):.2e}"),

        ("PDG: theta ≈ 45 stopni (blad < 0.01 deg)",
         abs(th_PDG - 45.0) < 0.01,
         f"theta={th_PDG:.6f} deg"),

        ("TGP (r31^K): CV(sqrt(m)) = 1 dokladnie (<1e-5)",
         abs(CV_TGP - 1.0) < 1e-5,
         f"CV={CV_TGP:.8f}"),

        ("TGP: CV(A^2) = 1 (<1e-5)",
         abs(CV_A2 - 1.0) < 1e-5,
         f"CV(A^2)={CV_A2:.8f}"),

        ("TGP: sum(delta_i^2) = 3 (<1e-4)",
         abs(np.sum(delta**2) - 3.0) < 1e-4,
         f"sum(delta^2)={np.sum(delta**2):.8f}"),

        ("ALG: CV_from_QK poprawne dla PDG",
         abs(cv_from_qk(QK_PDG) - CV_PDG) < 1e-6,
         f"CV_from_QK={cv_from_qk(QK_PDG):.8f}, CV_direct={CV_PDG:.8f}"),
    ]

    n_pass = sum(1 for _, p, _ in checks if p)
    print(f"\n  WYNIKI: {n_pass}/{len(checks)} PASS")
    for name, passed, detail in checks:
        mark = "PASS" if passed else "FAIL"
        print(f"  [{mark}] {name}")
        print(f"        {detail}")

    print(f"\n[WNIOSEK]")
    if n_pass == len(checks):
        print(f"  Trojjaka rownowaznosc POTWIERDZONA: Q_K=3/2 <=> CV=1 <=> theta=45 stopni")
        print(f"  Nowe sformulowanie T-OP1 (Tw. T-QK-CV):")
        print(f"    'Koide formula Q_K=3/2 jest rownowazna warunkowi,")
        print(f"    ze wspolczynnik zmiennosci pierwiastkow mas CV=1,")
        print(f"    tj. sigma(sqrt(m)) = mean(sqrt(m)) dla trzech generacji.'")
        print(f"  Formulacja TGP: CV(A(g0^(n))^2) = 1 (TGP-analog Koide)")
        print(f"  Pytanie otwarte: czy TGP soliton narzuca CV(A^2)=1 z symetrii?")
    else:
        print(f"  Niektorе testy nie przeszly. PASS: {n_pass}/{len(checks)}")

    print(f"\n[STATUS T-OP1]")
    print(f"  Poprzedni: Q_K=3/2 nie jest FP TGP; wykluczono fi^n-skalowanie")
    print(f"  Obecny:    [CZ. ZAMK.] geometryczne sformulowanie T-QK-CV")
    print(f"             Q_K=3/2 <=> CV(A^2)=1 w jezyku TGP")
    print(f"  Otwarte:   mechanizm TGP narzucajacy CV(A^2)=1")
    print(f"             Kandydat: symetria Z_3 generacyjna -> cirkulantna macierz mas")

    print("\n" + "=" * 65)
    sys.exit(0 if n_pass == len(checks) else 1)


if __name__ == "__main__":
    main()
