"""
SZCZEBEL 0 WIEZY -- kontrola skladnikow PRZED budowaniem.

Pytanie: czy sigma_ab = <s_i s_{i+a}>^TF ma NIEZEROWA wartosc prozniowa?
  - Jesli TAK  -> istnieje rozmaitosc prozni (RP^2 / SO(3)/D_2) -> defekty -> nosnik spinu.
  - Jesli NIE  -> sigma_ab = 0 w prozni, brak rozmaitosci, RP^2 nie ma domu w substracie.

Test: 3D Ising (H_Gamma v1, eq:B-H z dodatekB) Metropolis.
Mierzymy C_mu = <s_i s_{i+mu}> osobno dla mu = x,y,z, budujemy
    M_ij = sum_mu C_mu * mu_i mu_j,   sigma = M - tr(M)/3 * I
i sprawdzamy SKALOWANIE |sigma| z rozmiarem ukladu L.

KRYTERIUM ROZSTRZYGAJACE (to jest sedno):
  - prawdziwy parametr porzadku:  |sigma| -> const > 0  gdy L -> inf
  - zwykla fluktuacja:            |sigma| ~ 1/sqrt(N) = L^{-3/2}  -> 0
Dopasowujemy wykladnik p w |sigma| ~ L^p i porownujemy z -1.5.
"""
import numpy as np

TC = 4.5115  # 3D Ising, J=1, k_B=1

def run(L, T, nsweep=4000, nburn=1500, seed=0):
    rng = np.random.default_rng(seed)
    s = rng.choice([-1.0, 1.0], size=(L, L, L))
    # szachownica
    idx = np.indices((L, L, L)).sum(axis=0) % 2
    A, B = (idx == 0), (idx == 1)
    acc = np.zeros(3); nacc = 0
    for sw in range(nsweep):
        for mask in (A, B):
            nb = (np.roll(s, 1, 0) + np.roll(s, -1, 0) +
                  np.roll(s, 1, 1) + np.roll(s, -1, 1) +
                  np.roll(s, 1, 2) + np.roll(s, -1, 2))
            dE = 2.0 * s * nb
            flip = (dE <= 0) | (rng.random(s.shape) < np.exp(-np.clip(dE, 0, 60) / T))
            s = np.where(mask & flip, -s, s)
        if sw >= nburn:
            c = np.array([np.mean(s * np.roll(s, -1, ax)) for ax in range(3)])
            acc += c; nacc += 1
    return acc / nacc  # [C_x, C_y, C_z]

def sigma_from_C(C):
    M = np.diag(C)                    # sum_mu C_mu * e_mu e_mu
    return M - np.trace(M) / 3 * np.eye(3)

print("=" * 78)
print("SZCZEBEL 0: czy sigma_ab ma niezerowa wartosc prozniowa?")
print("=" * 78)

for T, tag in [(0.7 * TC, "gleboko uporzadkowana (T=0.7Tc)"),
               (1.00 * TC, "krytyczna (T=Tc)")]:
    print(f"\n--- {tag} ---")
    print(f"{'L':>4}{'N':>8}{'C_x':>12}{'C_y':>12}{'C_z':>12}{'|sigma|':>14}{'1/sqrt(N)':>12}")
    Ls, sigs = [], []
    for L in (6, 8, 10, 12):
        C = run(L, T, seed=L)
        sg = np.abs(sigma_from_C(C)).max()
        Ls.append(L); sigs.append(sg)
        print(f"{L:>4}{L**3:>8}{C[0]:>12.6f}{C[1]:>12.6f}{C[2]:>12.6f}{sg:>14.3e}{L**-1.5:>12.3e}")
    p = np.polyfit(np.log(Ls), np.log(sigs), 1)[0]
    print(f"  dopasowany wykladnik: |sigma| ~ L^({p:+.2f})   [fluktuacja: -1.50 | porzadek: 0.00]")
    verdict = "FLUKTUACJA (sigma -> 0)" if p < -0.8 else ("PORZADEK (sigma != 0)" if p > -0.4 else "NIEROZSTRZYGNIETE")
    print(f"  ==> {verdict}")

print("\n" + "=" * 78)
print("KONTROLA POZYTYWNA: czy test w ogole wykrylby porzadek nematyczny?")
print("=" * 78)
print("Wymuszamy anizotropie: J_z = 1.30, J_x = J_y = 1.0 (sprzezenie kierunkowe).")

def run_aniso(L, T, Jz=1.30, nsweep=3000, nburn=1200, seed=0):
    rng = np.random.default_rng(seed)
    s = rng.choice([-1.0, 1.0], size=(L, L, L))
    idx = np.indices((L, L, L)).sum(axis=0) % 2
    A, B = (idx == 0), (idx == 1)
    acc = np.zeros(3); nacc = 0
    for sw in range(nsweep):
        for mask in (A, B):
            nb = (np.roll(s, 1, 0) + np.roll(s, -1, 0) +
                  np.roll(s, 1, 1) + np.roll(s, -1, 1) +
                  Jz * (np.roll(s, 1, 2) + np.roll(s, -1, 2)))
            dE = 2.0 * s * nb
            flip = (dE <= 0) | (rng.random(s.shape) < np.exp(-np.clip(dE, 0, 60) / T))
            s = np.where(mask & flip, -s, s)
        if sw >= nburn:
            acc += np.array([np.mean(s * np.roll(s, -1, ax)) for ax in range(3)]); nacc += 1
    return acc / nacc

print(f"{'L':>4}{'C_x':>12}{'C_y':>12}{'C_z':>12}{'|sigma|':>14}")
Ls, sigs = [], []
for L in (6, 8, 10, 12):
    C = run_aniso(L, 0.7 * TC, seed=100 + L)
    sg = np.abs(sigma_from_C(C)).max()
    Ls.append(L); sigs.append(sg)
    print(f"{L:>4}{C[0]:>12.6f}{C[1]:>12.6f}{C[2]:>12.6f}{sg:>14.3e}")
p = np.polyfit(np.log(Ls), np.log(sigs), 1)[0]
print(f"  wykladnik: |sigma| ~ L^({p:+.2f})  ==> {'PORZADEK wykryty' if p > -0.4 else 'brak'}")
print("\n  (Ta kontrola pokazuje, ze test DZIALA: gdy anizotropia jest WLOZONA, sigma != 0.")
print("   Pytanie brzmi, czy substrat wytwarza ja SPONTANICZNIE.)")
