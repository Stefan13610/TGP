# Program: analityczne EOM N ciał TGP + warstwa numeryczna

**Zakres:** ten dokument i powiązany kod w **`nbody/`** są utrzymywane jako spójna warstwa pakietu. Ewentualne różnice względem głównego tomu LaTeX w `TGP_v1/` nie blokują prac tutaj — integrujesz je osobno, gdy zechcesz.

**Cel nadrzędny:** mieć **jasno rozdzielone** (i udokumentowane) dwa problemy:

| Warstwa | Pytanie | Co musi być „zamknięte” |
|--------|---------|-------------------------|
| **1 — Analityka TGP** | Jakiego potencjału / jakich sił używamy i jak wygląda **wektorowa** postać \(F_i = -\nabla_i V\)? | Zamknięte wzory na \(V_2\), postać \(V_3\) przez \(I(d_{ij},d_{ik},d_{jk})\), reguła łańcuchowa na \(\nabla_i I\) **bez** wchodzenia w szczegóły kwadratury. |
| **2 — Przeliczenie na współrzędne** | Jak policzyć \(I\) i \(\partial I/\partial d_{ab}\) oraz jak z czasu dostać \(q(t)\)? | Backend: całka Feynmana 2D (`three_body_force_exact`), siatka 3D (`triple_overlap_numerical`), ewentualnie przyszłe zamknięcie \(I_Y\); integrator: `leapfrog_integrate` / RK. |

Implementacja warstwy 1 z callable \(\partial I/\partial d\): moduł **`eom_tgp.py`** (`irreducible_three_body_forces_from_I_derivatives`, `accelerations_tgp_nbody`).  
Domyślne zamknięcie Coulomba na \(I\) jest zgodne z `dynamics_v2.three_body_forces_approximate(..., use_yukawa=False)` (test w `if __name__ == "__main__"`).

**TeX do głównego tomu:** `\input{tgp_nbody_lagrangian_eom}` — plik [`tgp_nbody_lagrangian_eom.tex`](tgp_nbody_lagrangian_eom.tex) (ścieżka względem katalogu roboczego LaTeX; zwykle ustaw `\input{nbody/tgp_nbody_lagrangian_eom}` z korzenia `TGP_v1/`).  
Opcjonalnie chaos numeryczny: [`tgp_lyapunov_benettin.tex`](tgp_lyapunov_benettin.tex) (`\input{nbody/tgp_lyapunov_benettin}`).

**Integratory:** `dynamics_backends.build_tgp_integration_pair(backend, ...)` → `(acc_fn, pot_fn)` dla `leapfrog_integrate` / `rk45_integrate` (`pairwise`, `coulomb_3b`, `yukawa_feynman`, `yukawa_saddle_approx`).

---

## Najważniejsze domknięcia (priorytety w `nbody`)

1. **Formalny zapis EOM (LaTeX + kod)**  
   - Jedna sekcja w istniejącym pliku `.tex` w `nbody/` lub krótki dodatek: \(L = \sum_i \tfrac12 C_i \|\dot x_i\|^2 - V\), \(V\) jak wyżej, jawne \(F_i^{(ijk)}\).  
   - Status: częściowo pokryte przez `three_body_force_exact.py` (wstęp) i teraz `eom_tgp.py`.

2. **Backend \(I_Y\) dla dynamiki**  
   - Do symulacji z pełną Yukawą: **używać** `three_body_forces_exact` (już EXACT), a nie przybliżenia saddle-point z `dynamics_v2` (duży błąd przy małym \(m d\)).  
   - Otwarte badawczo: zamknięta forma \(I_Y\) (plan P4 w `PLAN_ROZWOJU_NBODY.md`).

3. **Spójny interfejs „analityka → numeryka”**  
   - Zrobione: `dynamics_backends.build_tgp_integration_pair` wybiera backend i zwraca `(acc_fn, pot_fn)` pod leapfrog / RK45 (`ex138`–`ex140`).

4. **Redukcja i stałe ruchu**  
   - Translacja: \(\sum_i \mathbf{F}_i = 0\) — **zweryfikowane** (`total_force_from_accelerations`, `ex141`).  
   - Obrót: jeśli \(V\) zależy wyłącznie od \(d_{ij}\) (w tym \(I_Y(d_{12},d_{13},d_{23})\)), \(V\) jest \(SO(3)\)-invariante \(\Rightarrow\) \(\sum_i \mathbf{x}_i\times\mathbf{F}_i=0\) — **zweryfikowane** (`total_torque_about_origin_from_accelerations`, `ex143`).  
   - Uwaga: moment pędu liczony względem innego punktu odniesienia przesuwa \(\mathbf{L}\) o człon od \(\sum \mathbf{F}_i\); przy \(\sum \mathbf{F}_i=0\) sam moment sił nie zależy od wyboru bieguna.

5. **Regularyzacja**  
   - `softening` wyłącznie jako parametr numeryczny; dokumentować zbieżność \(\varepsilon \to 0\) przy publikacji wyników dynamiki.

---

## Kroki „kontynuacji” (konkretne)

- [x] Moduł `eom_tgp.py` z rozdzieleniem warstw 1/2 i testem zgodności Coulomb.  
- [x] Fragment TeX `tgp_nbody_lagrangian_eom.tex` (Lagranżjan, \(\mathbf F_i\), mapowanie na moduły).  
- [x] `dynamics_backends.build_tgp_integration_pair` — backendy `pairwise`, `coulomb_3b`, `yukawa_feynman`, `yukawa_saddle_approx` (`TGP_INTEGRATION_BACKENDS`).  
- [ ] Wpięcie `\input{...}` w `main.tex` / dodatek (ręcznie, gdy zechcesz).  
- [x] `examples/ex139_yukawa_feynman_leapfrog.py` — krótki leapfrog z dokładnym $I_Y$ (`n_quad=22`, $N=3$).  
- [x] `examples/ex140_yukawa_feynman_N4_scaling_rk.py` — \(N=4\), pomiar czasu `acc_fn`, leapfrog + RK45.  
- [x] `examples/ex141_net_force_translational_invariance.py` — \(\sum_i C_i a_i \approx 0\) dla wszystkich backendów.  
- [x] `examples/ex142_yukawa_overlap_quadrature_convergence.py` — zbieżność `n_quad` dla `yukawa_overlap_exact`.  
- [x] `examples/ex143_net_torque_so3_invariance.py` — \(\sum_i \mathbf{x}_i\times(C_i\mathbf{a}_i)\approx 0\) (SO(3)).  
- [x] `examples/ex144_conserved_P_L_leapfrog.py` — dryf \(\mathbf{P}\), \(\mathbf{L}\) przy leapfrog + `coulomb_3b`.  
- [x] `examples/ex145_minus_grad_V_matches_force.py` — zgodność \(-\nabla V\) (różnice centralne) z `C a` dla `coulomb_3b`.  
- [x] `examples/ex146_yukawa_feynman_grad_matches_force.py` — to samo dla `yukawa_feynman` (\(N=3\)).  
- [x] `examples/verify_nbody_eom_quick.py` — szybka regresja (`ex138` + `ex141`–`ex147`).  
- [x] `tgp_nbody_lagrangian_eom.tex` — redukcja do środka mas + zapis Hamiltona (\(p_i=C_i\dot x_i\), \(H=\sum |p_i|^2/(2C_i)+V\)).  
- [x] `dynamics_backends`: `center_of_mass`, `center_of_mass_velocity`, `conjugate_momenta`, `center_of_mass_acceleration`.  
- [ ] P4 z planu: badanie zamknięcia \(I_Y\) (niezależne od integratora).  
- [ ] P1 z planu (chaos): faza „publication-ready” (jeszcze dłuższe `t`, analiza zbieżności) — częściowo: `lyapunov.py` (`acceleration_jacobian_yukawa_feynman_analytic`), `ex164`–`ex166` (okno / refinement / siatka CSV), `ex167` (krotki test `|ΔE/E0|` leapfroga), `ex168` (odczyt/wykres CSV `ex166`), `ex169` (LF FD `J` vs analityczny), `ex170` (suma spektrum `6N` Yukawa), `ex171` (suma spektrum `6N` `coulomb_3b`), `ex172` (suma spektrum `6N` `pairwise`), `ex173` (FD vs analityczny `J` na spektrum `pairwise`), `ex174` (Newton vs Yukawa Feynman przy dopasowanej energii, Burrau `x`), `ex175` (trzy konteksty hamiltonowskie na tym samym `x`), `ex176` (Newton vs `coulomb_3b` przy dopasowanej energii), `ex177` (Newton vs pairwise `V_2` przy dopasowanej energii, leapfrog), `ex178` (tabela Newton + 3× TGP matched `H`), `ex179` (trzy hamiltoniany Coulomb na Burrau `x`), `ex180` (trzy hamiltoniany pairwise na Burrau `x`), `ex181` (siedem wierszy `v=0` vs `H=H_N`), `ex182` (CSV jak `ex181`), `ex183` (odczyt CSV `ex182`), `ex184` (skan `t_final` × siedem wierszy, CSV), `ex185` (odczyt CSV `ex184`), `ex186` (CSV jak `ex178`), `ex187` (odczyt CSV `ex186`), `ex188` (skan `t_final` × cztery wiersze matched `H`, CSV), `ex189` (odczyt CSV `ex188`), `ex190` (krótki leapfrog `|ΔE/E0|` matched-`H` jak `ex178`), `ex191` (krótki RK45 / DOP853, ta sama rodzina IC), `ex192` (tabela LF vs RK, wspólny `t_final`), `ex193` (CSV jak `ex192`), `ex194` (odczyt CSV `ex193`), `tgp_lyapunov_benettin.tex`, `tgp_nbody_lagrangian_eom.tex` (akapit Lyapunov), `ex148`–`ex194`, `verify_nbody_lyapunov_quick.py`.

---

*Nawigacja:* [README.md](README.md) · [PLAN_ROZWOJU_NBODY.md](PLAN_ROZWOJU_NBODY.md) · [ZALOZENIA_NBODY.md](ZALOZENIA_NBODY.md)
