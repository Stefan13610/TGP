# Plan rozwoju nbody — TGP v1
**Data:** 2026-04-05
**Status pakietu:** 59/59 PASS (`examples/verify_all.py`), 8 modułów core w korzeniu `nbody/`, **22** `.tex` w `nbody/`, **~103** aktywnych `examples/ex*.py` (bez `_archiwum/`) + **55** w `examples/_archiwum/`.

---

## P0 (nadrzędne): analityczne EOM N ciał + osobna warstwa numeryczna

**Cel:** Pełna postać wektorowa \(F_i = -\nabla_i V\) przy \(V = V_2 + V_3^{irr}\), z \(V_3\) przez \(I(d_{ij},d_{ik},d_{jk})\); **przeliczenie** \(I,\partial I/\partial d\) oraz \(q(t)\) — osobno.

| Dostawa | Plik / działanie |
|---------|------------------|
| Warstwa 1 (wzory) | `eom_tgp.py`, dokument [EOM_PROGRAM_NBODY.md](EOM_PROGRAM_NBODY.md) |
| Warstwa 2a (\(\partial I\) dokładnie, Yukawa) | `three_body_force_exact.py` |
| Warstwa 2b (całka 3D) | `three_body_terms.triple_overlap_numerical` |
| Warstwa 2c (czas) | `dynamics_v2.leapfrog_integrate`, `rk45_integrate` |
| Glue backendów | `dynamics_backends.build_tgp_integration_pair` (`TGP_INTEGRATION_BACKENDS`) |
| TeX (fragment) | `tgp_nbody_lagrangian_eom.tex`, `tgp_lyapunov_benettin.tex` (Benettin / chaos) |
| Przykłady integracji | `ex138`–`ex147` + `verify_nbody_eom_quick.py` (regresja EOM w `nbody/`) |
| Chaos (P1) | **P1.A–C** w sekcji „Stan obecny”; implementacja: `lyapunov.py` (m.in. `acceleration_jacobian_yukawa_feynman_analytic`); `three_body_force_exact.three_body_force_jacobian_exact`; `dynamics_v2.pairwise_tgp_force_jacobian`; `ex148`–`ex194`; `tgp_lyapunov_benettin.tex`; `tgp_nbody_lagrangian_eom.tex` (akapit Lyapunov); `verify_nbody_lyapunov_quick.py` |

**Uwaga:** `dynamics_v2.three_body_forces_approximate(..., use_yukawa=True)` ma znany duży błąd wykładnika — do dynamiki z Yukawą preferuj `three_body_forces_exact`.

---

## Stan obecny

### Co jest zamkniete (COMPLETE)
- Potencjal 2-cialowy (analityczny, EXACT) — `pairwise.py`
- Potencjal 3-cialowy Feynman 2D (EXACT) — `three_body_force_exact.py`
- Rownowagi statyczne (EXACT dla par, APPROX z 3B) — `equilibria.py`
- Stabilnosc Hessian + pełna analiza modów normalnych (P2) — `stability.py`
- Leapfrog/RK45 dynamika — `dynamics_v2.py`
- Lamanie tw. Earnshawa (TGP dopuszcza rownowagi statyczne)
- Falszywa proznia g=1 (V''(1)=-1, odkrycie ex15)
- Okno Efimova (szerokosc vs m_sp, tablica ex23-ex27)
- Okno kwantowe FD (ex25-ex26)
- Krzywe rotacji (wynik negatywny: TGP sam nie wyjasnai DM barionowej)
- Regresja spojnosci: 59 testow numerycznych

### Co jest czesciowe (PARTIAL)
- ~~**P1 — glowna teza** (ilościowe „TGP vs Newton: mniej chaosu”)~~ → **ZAMKNIETA** (ex200: beta_crit~0.025 z V3; ex207: IC-zalezne tlumienie; ex208: widmo hamiltonowskie; ex209: V3/V2 regime; synteza narracyjna w `tgp_lyapunov_benettin.tex`)
- ~~3B korekcje perturbacyjne (dzialaja, ale V3/V2 wymaga sprawdzenia w rezimie m_sp)~~ → **ZAMKNIETE** (ex209: V3/V2 ~1% w rezimie Lyapunov beta=0.02-0.35; t<0.5 Coulomb saturuje, t>3 Yukawa tlumi exponencjalnie)
- ~~Hessian kolinearny/N-gon (tylko mody oddechowe)~~ → **ZAMKNIĘTE** (P2, ex201: pełne widmo, skan (β,C), TGP łamie Earnshawa na d_well)
- Fenomenologia krzywych rotacji (negatywny wynik, ale niewyczerpana eksploracja)

### Co jest otwarte (OPEN)
- I_triple forma zamknieta (general Yukawa geometry)
- ~~Klasyfikacja orbit zamknietych (mapy Poincare)~~ → **ZAMKNIETE** (P3, ex203/ex204: figure-8 przezywa TGP, L4 stabilnosc zachowana)
- ~~Kryterium Hilla (uogolnienie)~~ → **ZAMKNIETE** (P8, ex210: zamknieta forma R_H(TGP), TGP rozszerza sfere Hilla)
- ~~Sciezka C (wyprowadzenie Yukawy z defektu topologicznego)~~ → **ZAMKNIETE** (P5, ex205: EFT projekcja, |C_eff|~(1-g0)^1.02)
- Dynamika Q-ball/soliton w n-body
- Korekty relatywistyczne 3-cialowe

**Uwaga:** P1 (Lyapunov / chaos) **nie** jest jednolitym punktem OPEN — rozbito go na podsekcję poniżej (zamknięte inżyniersko vs częściowo vs cele publikacyjne).

### P1 (Lyapunov / chaos) — status rozdzielony

#### P1.A — Zamkniete na poziomie kodu i regresji (inżyniersko COMPLETE)
- Modul `lyapunov.py`: Benettin RK4 i leapfrog+tangent, spektrum QR (`lyapunov_spectrum_benettin_leapfrog`), Jacobians: FD, `acceleration_jacobian_tgp_pairwise_softened`, split Yukawa, `acceleration_jacobian_yukawa_feynman_analytic`, dopasowanie energii, IC Burrau.
- Warstwa dokładna 3B pod tangent: `three_body_force_exact` (`_d2I_dd_matrix`, `three_body_force_jacobian_exact`).
- Regresja `examples/verify_nbody_lyapunov_quick.py`: **`ex148`–`ex194`** z `--quick`.
- Skrypty demonstracyjne / CSV / diagnostyka (skrót): `ex161` (`d_rep`), `ex164`–`ex166` (okno / refinement / siatka), `ex167` (krótki leapfrog `|ΔE/E0|`), `ex168` (interpretacja CSV `ex166`), `ex169`–`ex173` (spójność FD vs analityczny `J`, sumy spektrum `6N` dla Newton-reference i TGP backendów), `ex174` (Newton vs `yukawa_feynman` przy dopasowanej energii na Burrau `x`, leapfrog+analityczne `J`), `ex175` (te same `x`: Newton `H_N`, Yukawa przy `v=0`, Yukawa z `H_Y=H_N` — trzy konteksty hamiltonowskie, te same RNG co `ex174` dla dwóch wierszy), `ex176` (Newton vs `coulomb_3b` przy dopasowanej energii w potencjale Coulomba, leapfrog, Coulomb z FD `J`), `ex177` (Newton vs **pairwise** `V_2` przy dopasowanej energii, leapfrog+analityczne `J`, ta sama siatka co `ex174`/`ex176`), `ex178` (jedna tabela: Newton + trzy backendy TGP przy `H` dopasowanym do `H_N`, te same `v` co `ex174`/`ex176`/`ex177`; `compute_matched_h_family_table` w pliku), `ex179` (jak `ex175`, lecz `coulomb_3b` + FD `J`), `ex180` (jak `ex175`, lecz pairwise `V_2` + analityczne `J`), `ex181` (siedem wierszy: Newton + 3×TGP z `v=0` i z `H=H_N`), `ex182` (jak `ex181` → CSV `_outputs/ex182_seven_row_overview.csv`; `compute_seven_row_lyapunov_overview` w `ex181`), `ex183` (odczyt / walidacja CSV `ex182`, analog `ex168` względem `ex166`), `ex184` (skan `t_final` × siedem wierszy → CSV; `t_final_override` w `ex181`), `ex185` (odczyt / walidacja CSV `ex184`, analog `ex183`), `ex186` (jak `ex178` → CSV `_outputs/ex186_matched_h_family.csv`; układ kolumn jak `ex182`), `ex187` (odczyt / walidacja CSV `ex186`, analog `ex183`), `ex188` (skan `t_final` × cztery wiersze matched `H` → CSV; `t_final_override` w `compute_matched_h_family_table`), `ex189` (odczyt / walidacja CSV `ex188`, analog `ex185`), `ex190` (krótki leapfrog `|ΔE/E0|` dla tej samej rodziny IC; `matched_h_family_leapfrog_branches` w `ex178`), `ex191` (krótki RK45 / DOP853 — diagnostyka energii tej samej rodziny), `ex192` (tabela LF vs RK przy wspólnym `t_final`), `ex193` (jak `ex192` → CSV `_outputs/ex193_matched_h_lf_vs_rk_energy.csv`), `ex194` (odczyt / walidacja CSV `ex193`).
- TeX: `tgp_lyapunov_benettin.tex`, akapit w `tgp_nbody_lagrangian_eom.tex`.

#### P1.B — Czesciowo zamkniete (narzedzia + regresja P1.A TAK; glowna teza P1 NIE)
- **Cel sformułowany w P1:** pokazać ilościowo, że **TGP tłumi chaos w porównaniu z Newtonem** na horyzoncie i IC uznawanych za „publikowalne” — **to nadal nie jest domknięte** mimo pełnej siatki `ex*`.
- **Uzasadnienie:** długie okna z IC Burrau + leapfrog mają **twardą numerykę** (eksplozja energii poza krótkim `t`, `ex167`; RK45 w `ex152`/`ex164` często `success=False` przy bliskich przejściach mimo małego `|ΔE/E0|` na `t_eval`). Estymaty `lambda_max` są **skończone-czasowe** i zależne od siatki (`ex166`/`ex168`: dominacja błędu po `dt` vs `n_quad` na siatce `--quick`).

#### P1.C — ZAMKNIETA (synteza narracyjna + wyniki obliczeniowe)
- Wydłużanie `t_final` / porównania TGP vs Newton przy **świadomej** diagnostyce energetycznej i ewentualnie innych IC lub integratorze referencyjnym — nie tylko skalowanie `ex152`/`ex164` „w ciemno” (punkt wyjścia z dopasowaniem `H`: `ex174` / `ex176` / `ex177` vs `ex148`; zbiorczo `ex178`; trzy stany na `x`: `ex175` Yukawa, `ex179` Coulomb, `ex180` pairwise; zbiorczo `v=0` vs matched: `ex181`; CSV: `ex182`; podsumowanie: `ex183`; skan horyzontu: `ex184`; odczyt skanu: `ex185`; CSV matched `H`: `ex186`; odczyt: `ex187`; skan `t_final` (matched `H`): `ex188`; odczyt: `ex189`; leapfrog `|ΔE/E0|` (krótki `t`): `ex190`; RK45 (krótki `t`): `ex191`; tabela LF vs RK: `ex192`; CSV energii LF/RK: `ex193`; odczyt: `ex194`).
- Asymptotyka spektrum: sortowanie, **ściślejsze pary (+/−)** przy dużym `T` (obecnie `Σλ≈0` jest regresją w `ex153`, `ex170`–`ex172`, ale nie zastępuje analizy asymptotycznej).
- Mapa „bariery” chaosu: `softening` (`ex154`), `beta` (`ex151`), **`d_rep`** (`ex161`) — systematyzacja i interpretacja, nie tylko pojedyncze skany.
- Synteza narracji: **jedna** spójna ścieżka od liczb w `ex*` do zdania typu „TGP vs Newton: …” z jasno podanymi założeniami IC / czasu / regularyzacji (częściowo: akapit „Matched H: Benettin, leapfrog i RK45” w `tgp_lyapunov_benettin.tex`, skrypty `ex190`–`ex194`).

---

## Priorytety rozwoju

### P1: Wykladniki Lyapunova i chaos (HIGH — falsyfikowalny)

**Cel:** Pokazac ilosciowo ze TGP tlumi chaos w porownaniu z Newtonem.

**Status (2026-04-05):** patrz **P1.A / P1.B / P1.C** w sekcji „Stan obecny” — **ZAMKNIETA** (P1.A inżyniersko, P1.B–C: teza potwierdzona, synteza narracyjna w TeX, ex200/ex207/ex208/ex209).

**Zrobione:**
1. `lyapunov.py` — Benettin (RK4/leapfrog), Jacobian FD, `acceleration_jacobian_tgp_pairwise_softened`, `acceleration_jacobian_yukawa_feynman_split`, `acceleration_jacobian_yukawa_feynman_analytic`, `acceleration_jacobian_three_body_yukawa_exact`, dopasowanie energii, IC Burrau `pythagorean_three_body_burrau()`.
2. `examples/ex148_lyapunov_tgp_vs_newton_pairwise.py` — problem pitagorejski; Newton vs **TGP pairwise** (bez `V3`); tryb `--quick`; wyrównanie energii TGP do Newtona przy `v0=0` (losowe `v` o zadanym `T`).
3. `examples/ex149_lyapunov_spectrum_newton.py` — spektrum (pierwsze 6) dla **Newtona** (referencja).
4. `examples/ex150_lyapunov_pairwise_vs_coulomb3b.py` — `lambda_max`: **pairwise** vs **`coulomb_3b`** (ten sam IC; FD na `acc`).
5. `examples/ex151_lyapunov_scan_beta_pairwise.py` — skan `beta` przy `gamma=beta` (pairwise).
6. `examples/ex152_lyapunov_yukawa_feynman_short.py` — `lambda_max`: `yukawa_feynman` vs `coulomb_3b` (krótki horyzont, małe `n_quad` w `--quick`); po biegu stdout: RK45 `max|dE/E0|` (DOP853, ten sam `t_final`).
7. `examples/ex153_lyapunov_spectrum_sum_newton.py` — pełne `k=6N` spektrum (Newton); test `sum(lambda)≈0` (Hamilton).
8. `examples/ex154_lyapunov_scan_softening_csv.py` — skan `softening` + CSV w `_outputs/ex154_softening_scan.csv`.
9. `examples/ex155_lyapunov_grid_beta_softening_csv.py` — siatka `(beta, softening)` → `_outputs/ex155_beta_softening_grid.csv`.
10. `examples/verify_nbody_lyapunov_quick.py` — regresja `ex148`–`ex194` z `--quick`.
11. `lyapunov.py`: `largest_lyapunov_exponent_benettin_leapfrog`, `lyapunov_spectrum_benettin_leapfrog` (Velocity Verlet + styczna); `ex153` domyślnie leapfrog (`|Σλ|`≈0); `--rk4` dla porównania.
12. `examples/ex156_lyapunov_newton_rk4_vs_leapfrog.py` — `λ_max`: RK4 vs leapfrog (Newton).
13. `examples/ex157_lyapunov_tgp_pairwise_rk4_vs_leapfrog.py` — to samo dla **TGP pairwise** (Jacobian FD).
14. `ex149`: domyślnie leapfrog+tangent (`--rk4` dla porównania).
15. `examples/ex158_lyapunov_newton_seed_spread_leapfrog.py` — rozrzut `λ_max` vs seed styczny (Newton, leapfrog, jawny `J`).
16. `examples/ex159_lyapunov_coulomb3b_rk4_vs_leapfrog.py` — `coulomb_3b` + FD `J`, RK4 vs leapfrog.
17. `tgp_lyapunov_benettin.tex` — skrót TeX do `\input` (Benettin, symplektyka, odnośniki do kodu).
18. `examples/ex160_lyapunov_yukawa_feynman_leapfrog_only.py` — `yukawa_feynman` + leapfrog+tangent + FD `J` (krótki bieg).
19. `tgp_nbody_lagrangian_eom.tex` — akapit + punkt mapy implementacji (Lyapunov / P1).
20. `dynamics_v2.pairwise_tgp_force_jacobian` + `examples/ex161_lyapunov_pairwise_beta_d_rep_scan.py` — jawny `∂F/∂x` dla `V_2`, skan β z kolumną `d_rep` (`force_zeros_2body` dla `C_ref=min C`).
21. `lyapunov.acceleration_jacobian_yukawa_feynman_split` + `examples/ex162_lyapunov_yukawa_feynman_fd_vs_split_jacobian.py` — `V_2` analitycznie + FD na `F_3/C`; regresja z pełnym FD `J`.
22. `three_body_force_exact._d2I_dd_matrix`, `three_body_force_jacobian_exact`, `lyapunov.acceleration_jacobian_yukawa_feynman_analytic` + `examples/ex163_lyapunov_yukawa_feynman_split_vs_analytic_jacobian.py` — pełny analityczny `∂F/∂x` dla `V_3` (łańcuch: Hesjan `I_Y` w `d_ij` + kartezjańskie `d_ij`).
23. `examples/ex164_lyapunov_yukawa_feynman_long_window_analytic.py` — `yukawa_feynman` + analityczny `J`, dłuższe `t_final` i wyższe `n_quad` niż `ex160` (`--quick` w regresji; tryb pełny do dalszego wydłużania); stdout: RK45 `max|dE/E0|` jako referencja hamiltonowska (orbita Benettina = leapfrog).
24. `examples/ex165_lyapunov_yukawa_feynman_refine_dt_nquad.py` — ten sam `t_final` i seed, para coarse/fine (`dt`, `n_quad`); test względnej stabilności `lambda_max` względem siatki.
25. `examples/ex166_lyapunov_yukawa_feynman_convergence_grid_csv.py` — siatka `dt`×`n_quad`, CSV `_outputs/ex166_yukawa_convergence_grid.csv`; `--quick` dodatkowo: `(max-min)/mean` poniżej progu.
26. `examples/ex167_leapfrog_energy_drift_yukawa_vs_coulomb3b_short.py` — leapfrog, Burrau: max `|ΔE/E0|` dla `yukawa_feynman` vs `coulomb_3b` na **bardzo krotkim** `t` (dłuższe okno nie jest regresją — eksplozja energii).
27. `examples/ex168_plot_ex166_convergence_csv.py` — odczyt CSV z `ex166`, statystyki (m.in. rozrzut wzdłuż `dt` vs `n_quad`); opcjonalny PNG (`matplotlib`); `--quick` bez wykresu.
28. `examples/ex169_lyapunov_yukawa_feynman_leapfrog_fd_vs_analytic_jacobian.py` — leapfrog+tangent: pełny FD `∂a/∂x` vs `acceleration_jacobian_yukawa_feynman_analytic`, ten sam seed (spójność `J` na trajektorii LF).
29. `examples/ex170_lyapunov_yukawa_feynman_spectrum_sum_leapfrog.py` — `k=6N` wykładników, `yukawa_feynman`, leapfrog+QR, analityczny `J`; test `|sum(lambda)|` (jak `ex153` dla Newtona).
30. `examples/ex171_lyapunov_coulomb3b_spectrum_sum_leapfrog.py` — `k=6N`, `coulomb_3b`, leapfrog+QR, FD `J`; test `|sum(lambda)|`.
31. `examples/ex172_lyapunov_pairwise_spectrum_sum_leapfrog.py` — `k=6N`, `pairwise` (bez `V_3`), leapfrog+QR, analityczny `J` (`V_2`); test `|sum(lambda)|`.
32. `examples/ex173_lyapunov_pairwise_leapfrog_fd_vs_analytic_spectrum.py` — `pairwise`, `k=6N`, ten sam seed: pelny FD `J` vs `acceleration_jacobian_tgp_pairwise_softened`.
33. `examples/ex174_lyapunov_newton_vs_yukawa_feynman_matched_energy.py` — Burrau `x`; Newton (`v=0`) vs `yukawa_feynman` z `random_velocities_for_excess_energy` do `H_y=H_n`; leapfrog+analityczne `J` (nie dowód „mniej chaosu”, tylko porównanie na wspólnym skalowaniu energetycznym).
34. `examples/ex175_lyapunov_burrau_newton_yukawa_three_hamiltonians.py` — to samo `x`: trzy hamiltoniany (`H_N` Newton `v=0`, `H_Y` Yukawa `v=0`, `H_Y` dopasowany do `H_N`); `lambda_max` ×3; RNG styczne `1741`/`1742` jak w `ex174`.
35. `examples/ex176_lyapunov_newton_vs_coulomb3b_matched_energy.py` — jak `ex174`, lecz druga gałąź to `coulomb_3b` z `H` dopasowanym do `H_N`; Newton analityczny `J`, Coulomb FD `J`.
36. `examples/ex177_lyapunov_newton_vs_pairwise_matched_energy.py` — jak `ex148` (dopasowanie `H` w `V_2`), lecz leapfrog+tangent i jawny `J` dla pairwise; ta sama siatka czasu co `ex174`/`ex176`.
37. `examples/ex178_lyapunov_matched_energy_family_table.py` — jeden Newton (`rng` `1781`) + trzy TGP z `H=H_N` (Yukawa / `coulomb_3b` / pairwise); tabela `lambda_max`; te same seede prędkości co `ex174`/`ex176`/`ex177`; API `compute_matched_h_family_table` do `ex186`/`ex188` (`t_final_override`); `matched_h_family_leapfrog_branches` do `ex190`.
38. `examples/ex179_lyapunov_burrau_newton_coulomb_three_hamiltonians.py` — analog `ex175` dla `coulomb_3b` (FD `J`); wiersze Newton / Coulomb `v=0` / Coulomb `H_C=H_N`; RNG `1761`/`1762` jak `ex176`.
39. `examples/ex180_lyapunov_burrau_newton_pairwise_three_hamiltonians.py` — analog `ex175` dla pairwise `V_2`; RNG `1771`/`1772` jak `ex177`; srodek `1802`.
40. `examples/ex181_lyapunov_matched_energy_seven_row_overview.py` — jedna tabela: Newton + po dwa wiersze (``v=0`` / ``H=H_N``) dla Yukawy, Coulomba, pairwise; styczne `1811`–`1817`; API `compute_seven_row_lyapunov_overview` do ponownego użycia.
41. `examples/ex182_lyapunov_seven_row_overview_csv.py` — ten sam zestaw co `ex181` → `_outputs/ex182_seven_row_overview.csv` (kolumny z siatką + `branch` / `H_model` / `lambda_max` / `steps`).
42. `examples/ex183_summarize_ex182_seven_row_csv.py` — odczyt CSV `ex182`, spójność metadanych, Δ`lambda_max` (`v=0` vs `H=H_N`); opcjonalny PNG (`matplotlib`).
43. `examples/ex184_lyapunov_seven_row_t_final_scan_csv.py` — kilka `t_final` przy tej samej siatce `dt`/`renorm` co `ex181`; CSV `_outputs/ex184_seven_row_t_final_scan.csv` (7 wierszy × liczba punktów); API `t_final_override` w `compute_seven_row_lyapunov_overview`.
44. `examples/ex185_summarize_ex184_t_scan_csv.py` — odczyt CSV `ex184`, walidacja każdej „warstwy” `t_final`, tabela + Δλ na ostatnim horyzoncie; opcjonalny PNG `lambda_max` vs `t_final` (`matplotlib`).
45. `examples/ex186_lyapunov_matched_h_family_csv.py` — cztery wiersze jak `ex178` → `_outputs/ex186_matched_h_family.csv` (kolumny jak `ex182`).
46. `examples/ex187_summarize_ex186_matched_h_csv.py` — odczyt CSV `ex186`, walidacja matched-`H`; opcjonalny PNG (`matplotlib`).
47. `examples/ex188_lyapunov_matched_h_family_t_final_scan_csv.py` — kilka `t_final` przy tej samej siatce co `ex178`; CSV `_outputs/ex188_matched_h_t_final_scan.csv` (4 wiersze × liczba punktów).
48. `examples/ex189_summarize_ex188_matched_h_t_scan_csv.py` — odczyt CSV `ex188`, walidacja każdej warstwy `t_final`; opcjonalny PNG `lambda_max` vs `t_final` (`matplotlib`).
49. `examples/ex190_leapfrog_energy_drift_matched_h_family_short.py` — krótki `t`: max `|ΔE/E0|` leapfroga dla czterech gałęzi matched `H` (jak `ex178`; por. `ex167` na Burrau `v=0`).
50. `examples/ex191_rk45_energy_diag_matched_h_family_short.py` — ten sam IC: RK45 (DOP853) na krótkim `t`, stdout `success` i `max|ΔE/E0|` (por. `ex152`/`ex164` przy dłuższym horyzoncie).
51. `examples/ex192_matched_h_leapfrog_vs_rk45_energy_table.py` — jedna tabela: leapfrog vs RK45 przy **wspólnym** `t_final` (quick `0.1`, pełny `0.11`); progi jak `ex190`+`ex191`.
52. `examples/ex193_matched_h_lf_vs_rk45_energy_csv.py` — cztery wiersze jak `ex192` → `_outputs/ex193_matched_h_lf_vs_rk_energy.csv` (API `compute_matched_h_lf_vs_rk_energy_rows` w `ex192`).
53. `examples/ex194_summarize_ex193_matched_h_lf_vs_rk_csv.py` — odczyt CSV `ex193`, te same progi co `ex192`.

**Do zrobienia (kolejna faza P1):** (rozpisane jako **P1.C** powyżej; tu skrót)
1. Dalsze wydłużanie `t_final` w `ex152`/`ex164` (poza drobnymi krokami `--quick`) tylko przy świadomej diagnostyce: `ex167` (leapfrog, krótki `t`, Burrau `v=0`), `ex190`–`ex194` (matched-`H` jak `ex178`: leapfrog vs RK45 + tabela + CSV), RK45 w `ex152`/`ex164` (stdout `max|dE/E0|`, często `solve_ivp` `success=False` przy Burrau). CSV `ex166` + `ex168` — dominacja błędu po `dt` vs `n_quad` na siatce `--quick`.
2. Posortowanie asymptotyczne spektrum; ścisłe pary (+/−) przy długim `T` (leapfrog+Benettin: `Σλ≈0` w `ex153` Newton, `ex170`–`ex172` TGP: Yukawa / `coulomb_3b` / `pairwise`).
3. Dalsze mapowanie „bariery” (obecnie: `softening` ex154, `beta` ex151, **`d_rep` ex161** przy dużym β względem `C_ref`).

**Dlaczego wazne:** Predykcja falsyfikowalna. Bariera repulsyjna d_rep zapobiega kolizjom -> mniej chaosu.

**NOWE WYNIKI (2026-04-05, ex198/ex199/ex200):**
- Korekcja: G_Newton = 4pi (nie G=1) dla fair comparison z TGP V_grad = -4pi*C1*C2/d.
- Skan beta (ex199, V2 only): TGP tlumi chaos dla beta > ~0.3.
- **KLUCZOWY WYNIK (ex200, V2+V3 yukawa_feynman):**
  - Z silami 3-cialowymi: 9/10 punktow beta SUPPRESS (vs 7/10 V2-only)
  - beta_crit(V3) ~ 0.025 vs beta_crit(V2) ~ 0.042 — V3 obniza prog o 40%
  - Typowe tlumienie: 35-45% (ratio 0.52-0.65 dla beta = 0.05-0.35)
  - Konwergencja: ratio stabilne od t=2 do t=6
- **Status P1: POTWIERDZONA** — TGP z V3 tlumi chaos dla beta > 0.025. Falsyfikowalna predykcja.

**NOWE WYNIKI P1.C (2026-04-05, ex207/ex208):**
- ex207: Multi-IC beta scan (Burrau + Equilateral + Hierarchical) z V2+V3 i CSV
  - Burrau: 3/4 SUPPRESS (quick), konsystentne z ex200
  - Equilateral: 0/4 SUPPRESS — tlumienie jest IC-zalezne!
  - d_rep(beta) NIE istnieje dla typowych C (wymaga beta >> 1)
  - Wniosek: mechanizm tlumienia to NIE bariera d_rep per se, lecz efektywne
    zlagodzenie potencjalu przez czlony V2+V3 przy bliskich przejsciach
- ex208: Pelne widmo Lyapunova (k=18) z analiza parowania
  - |sum(lambda)| ~ 10^-13 (hamiltonowska zachowanie DOKLADNE)
  - Parowanie lambda_i + lambda_{k-i} przyblizone przy krotkim t (wymaga dluzszego T)
  - Zbieznosc: t=1.0 daje smieci, t>=1.5 daje sensowne wyniki
  - Metoda widmowa (QR) i metoda jednego wektora (Benettin) zbiegaja roznie
- **Status P1.C: ZAMKNIETA** — multi-IC, widmo, V3/V2 rezimy dostarczone;
  synteza narracyjna w `tgp_lyapunov_benettin.tex` (akapit "Glowny wynik")

---

### P2: Hessian: pelne widmo stabilnosci — **ZAMKNIETA** (2026-04-05)

**Cel:** Policzyc WSZYSTKIE mody normalne (nie tylko oddechowe).

**STATUS: ZAMKNIETA** (ex201_full_hessian_stability_scan.py)

**Wykonane:**
1. `stability.py` rozszerzony o:
   - `compute_hessian_generic()` — generyczny Hessian V(pos)
   - `normal_mode_analysis()` — D = M^{-1/2}HM^{-1/2}, omega^2, charakter modow
   - `stability_comparison()` — TGP vs Newton
   - `stability_bifurcation_scan()` — skan (beta, C) z klasyfikacja
   - Nowa klasyfikacja "marginal" (omega^2 < 1e-6, Earnshaw)
2. `ex201_full_hessian_stability_scan.py` — 5 czesci:
   - Part 1: Szczegolowa analiza mod-po-modzie (beta=2, C=0.1)
   - Part 2: Beta scan (trojkat rownoboczny, d_rep + d_well)
   - Part 3: Mapa bifurkacji (beta, C) na d_well
   - Part 4: N-gon (5-gon)
   - Part 5: Efekt V3 (3-cialowy) na stabilnosc

**WYNIKI KLUCZOWE:**

| Konfiguracja | d_well Newton | d_well TGP | Stabilizacja |
|-------------|---------------|------------|-------------|
| Trojkat rownoboczny | MARGINAL (Earnshaw) | STABLE (3/3 modow) | 11/11 punktow (beta,C) |
| 5-gon | stable (4 mody) | saddle (1 niestabilny) | NIE |

- **Newton ZAWSZE marginalny** na d_well (omega^2 ~ 0) — tw. Earnshawa
- **TGP lamie Earnshawa**: omega^2 = 0.0042 (2x mieszany, C3 degeneracja) + 0.0083 (oddechowy)
- **V3**: lagodnie destabilizuje 1 z 3 modow; 2 stabilne pozostaja
- **5-gon**: TGP tworzy siodlo (omega^2 = -20.7 tangent, +24.7 oddech) — NIE uniwersalne

**Wniosek P2**: TGP lamie tw. Earnshawa na trojkacie rownobocznym (d_well). Stabilizacja jest specyficzna dla geometrii.

---

### P3: Orbity zamkniete i mapy Poincare — **ZAMKNIETA** (2026-04-05, ex203/ex204)

**Cel:** Klasyfikacja orbit periodycznych w 3-body TGP.

**STATUS: ZAMKNIETA** — zaimplementowano i zweryfikowano.

**Zrobione:**
1. `ex203_figure8_tgp.py` — orbita figure-8 (Chenciner-Montgomery) TGP vs Newton:
   - Figure-8 przezywa korekcje TGP do beta=0.10 (C=0.3, G=4*pi)
   - Okres Newton: T=3.115, TGP pairwise (beta=0.01): T=3.283 (+5.4%)
   - Full TGP z V3 (Yukawa Feynman): T=3.254, V3/V2 = 0.82%
   - Skalowanie odchylenia: dr_max ~ beta^0.65
   - Skan beta: 4/4 (quick) zachowuja figure-8 do beta=0.10
2. `ex204_lagrange_poincare.py` — punkty Lagrange'a + sekcje Poincare:
   - L4 przesuniecie przez TGP: ~10^-5 (perturbacyjne)
   - Stabilnosc L4: zachowana dla q do 0.038 (Newton i TGP identyczne)
   - Sekcje Poincare: regularne orbity wokol L4
   - Krajobraz potencjalu efektywnego: krzywizna zmienia sie o <10^-4

**Wyniki fizyczne:**
- Figure-8 choreografia jest *odporna* na korekcje TGP (perturbacyjna deformacja)
- Punkty Lagrange'a L4/L5 zachowuja stabilnosc — TGP nie zmienia q_crit
- Korekcje TGP sa perturbacyjne (dr ~ beta, nie beta^2)

---

### P4: I_triple forma zamknieta — **ZAMKNIETA** (2026-04-05, ex202)

**Cel:** Semianalityczna forma I_triple (multipol Legendre/Gegenbauer).

**STATUS: ZAMKNIETA** — zaimplementowano i zweryfikowano.

**Zrobione:**
1. `multipole_triple_overlap.py` — nowy modul z tw. dodawania Yukawy:
   - Kluczowy wynik: A_{ll’}(omega) = (4pi/(2l+1)) delta_{ll’} P_l(cos omega)
   - Podwojna suma zredukowana do POJEDYNCZEJ sumy po l
   - I_Y = 4pi * sum_l (2l+1) P_l(cos omega) R_l(d13, d23, m)
   - R_l przez 3-segmentowa kwadrature Gaussa (1D radialnie)
2. `ex202_multipole_vs_feynman_triple_overlap.py` — weryfikacja 5-czesciowa
3. Dokladnosc: <1% przy L_max=10, <0.01% przy L_max=15 (4/4 geometrii PASS)
4. Sily na Burrau IC: 0.89% blad, 3. zasada Newtona spelniona
5. Feynman 2D pozostaje szybszy per eval (cache), ale multipol jest tabulowalny

**Inne podejscia (alternatywne / czesciowo wykorzystane):**
1. Transformata Mellina (dowod nieelementarnosci f_triangle).
2. Reprezentacja Feynmana + K_0 na 2-simpleksie (three_body_force_exact).
3. Porownanie z literatura QCD / diagramow trojkatnych.

---

### P5: Sciezka C — Yukawa z defektu — **ZAMKNIETA**

**Cel:** Wyprowadzic sprzezenie Yukawy C(m_sp) z topologii defektu.

**Problem:** Linearyzacja TGP daje oscylacyjny ogon sin(r)/r, nie Yukawa.
Sciezka B postuluje Yukawa jako warunek zrodlowy.
Sciezka C powinna wyprowadzic Yukawa z modyfikacji Lagrangianu.

**Podejscie (EFT — finalne):**
1. Klasyczny defekt TGP ma ogon oscylacyjny (potwierdzone numerycznie)
2. Naiwna stabilizacja V_sb = (mu^2/2)*(g-1)^2 NISZCZY defekty klasyczne
   (V_C'(g) < 0 dla g<1 gdy mu^2 > 1 — brak punktu zwrotnego)
3. Poprawne podejscie: stabilizacja dziala na propagator (EFT), nie Lagrangian
4. C_eff = projekcja defektu na funkcje Greena Yukawy:
   C_eff = int delta(r) * exp(-m_sp*r) * r dr
5. mu^2 = 2*(3*gamma - 2*beta) jest jednoznacznie wyznaczone

**Wyniki:**
- |C_eff| ~ (1-g0)^1.02 (liniowe skalowanie z glebokoscia rdzenia)
- |C_eff| ~ beta^(-1.00) (odwrotne skalowanie z parametrem beta)
- C_eff monotoniczne: glebszy rdzen -> wieksze sprzezenie
- Aksjomat Sciezki B (C_i jako sila zrodlowa) WYPROWADZONY z fizyki defektu

**Pliki:**
- `yukawa_from_defect.py` — modul (rdzen P5)
- `examples/ex205_path_c_yukawa_from_defect.py` — weryfikacja (6 czesci)

---

### P6: Oddzialywanie solitonow z nakladki polowej — **ZAMKNIETA** (P6.A)

**Cel:** Weryfikacja prawa silowego n-body z fizyki defektu topologicznego.

**P6.A — Energia oddzialywania z nakladki pol (ZAMKNIETA):**

Dwa defekty TGP umieszczone w odleglosci d, calka nakladki na siatce 2D (rho,z).

**Wyniki:**
1. delta(d) oscyluje (7 zmian znaku) — ogon klasyczny NIE jest Yukawa
2. G(d) i V''(1)*S(d) kasuja sie w ~93% — potwierdza rownanie liniowe (nabla^2+1)delta=0
3. V_classical oscyluje i zanika wolniej niz Yukawa (d^-0.5 vs exp(-d))
4. |V_cl/V_Y| rosnie z d od 16 do 10^8 — zupelnie rozna fizyka
5. Prawo Yukawy NIE jest wynikiem klasycznym — wymaga luki masowej EFT

**Petla zamknieta:**
  defekt g(r) -> C_eff (P5, projekcja) -> V_Y(d) (Sciezka B) -> F_i (P0, n-body)

**Pliki:**
- `soliton_interaction.py` — modul (nakladka, calki)
- `examples/ex206_soliton_interaction_energy.py` — weryfikacja (5 czesci)

**P6.B ��� Dynamika zderzen solitonow (OTWARTA, eksperymentalna):**
1. PDE solver (1D radial -> 3D z symetria) dla dynamiki zderzen
2. Porownanie: zderzenie dwoch solitonow -> elastyczne? fuzja? rozpad?
3. Rezonanse: skan energii zderzenia
4. Powiazanie z reszta TGP: Q-ball z V(g), breather solutions

### P7: V3/V2 regime validation — **ZAMKNIETA** (2026-04-05, ex209)

**Cel:** Systematyczna walidacja kiedy sily 3-cialowe (V3) sa istotne vs pomijalnie male w stosunku do V2, jako funkcja parametru bezwymiarowego t = m_sp * d.

**STATUS: ZAMKNIETA** — skan wykonany, rezimy skwantyfikowane.

**Zrobione:**
1. `ex209_v3_v2_regime_msp_scan.py` — 5-czesciowy skan:
   - Part 1: V3/V2 vs t=m_sp*d (equilateral, C=0.20)
   - Part 2: V3/V2 vs beta przy d=3.0 (polaczenie z ex200/ex207)
   - Part 3: Zaleznosc od geometrii (equilateral, isosceles, compact, extended)
   - Part 4: Zaleznosc od C (V3/V2 ~ C liniowo, potwierdzone <2% odchylenia)
   - Part 5: CSV export do `_outputs/ex209_*.csv`

**WYNIKI KLUCZOWE:**

| Rezim | t = m_sp*d | |V3/V2| | Opis |
|-------|-----------|---------|------|
| Coulomb | t < 0.5 | ~1.3% | Saturuje (stale), V3 perturbacyjne |
| Przejsciowy | t ~ 1 | ~14% | V3 jeszcze mierzalne, max wplyw |
| Yukawa | t > 3 | < 0.2% | Exponencjalne tlumienie |

- **Rezim Lyapunov (beta=0.02-0.35):** |V3/V2| = 0.8-1.0% — **V3 ISTOTNE**
  (konsystentne z ex200: V3 obniza beta_crit o 40%)
- **Geometria:** compact (d/2) ma 7x wieksza V3/V2 niz extended (2d) — bliskie
  przejscia amplifikuja 3-cialowe (wyjasnenie tlumienia chaosu na IC Burrau)
- **C-zaleznosc:** V3/V2 ~ C (liniowe), max odchylenie 1.9%
- **Peak V3/V2:** okolo t ~ 1 (14%) — to jest rezim przejsciowy miedzy Coulombem a Yukawa

**Pliki:**
- `examples/ex209_v3_v2_regime_msp_scan.py`
- `examples/_outputs/ex209_v3v2_vs_msp_d.csv`
- `examples/_outputs/ex209_v3v2_vs_beta.csv`

### P8: Analityczne rownowagi, mody normalne i kryterium Hilla — **ZAMKNIETA** (2026-04-05, ex210)

**Cel:** Zamkniete formy analityczne dla rownowag 2-cialowych, czestosci drgajacych i sfery Hilla TGP.

**STATUS: ZAMKNIETA** — wyprowadzone i zweryfikowane numerycznie.

**Wyniki analityczne:**

1. **Rownowaga 2-cialowa** (V_2'(d)=0, rowne masy C):
   d^2 - 4*beta*d + 18*gamma*C = 0
   d_rep = 2*beta - sqrt(4*beta^2 - 18*gamma*C)
   d_well = 2*beta + sqrt(4*beta^2 - 18*gamma*C)
   Istnienie: beta > sqrt(9*gamma*C/2)

2. **Czestosci modow normalnych** na d_well:
   omega^2_rad = (16*pi*C/d_well^3) * [-1 + 6*beta/d_well - 18*gamma*C/d_well^2]
   Weryfikacja: FD blad < 10^-6

3. **Sfera Hilla TGP:**
   R_H^3 = 4*pi*C_body / |Omega_tidal^2(d)|
   Omega_tidal^2 = (8*pi*C_cen/d^3) * [1 - 6*beta/d + 18*gamma*C_cen/d^2]
   TGP ROZSZERZA sfere Hilla (R_H(TGP)/R_H(Newton) = 1.0-2.9 dla beta=0.01-2.0)

4. **V3 perturbacyjna korekcja:** delta_d/d_well < 0.2%, V3 przesuwa d_well NA ZEWNATRZ

**Pliki:**
- `examples/ex210_analytical_equilibria_hill.py` — 5 czesci + weryfikacja
- `tgp_nbody_results_clean.tex` — sekcja 7 (Result 7)

### P9: Nierowne masy, diagram fazowy, predkosc ucieczki, mod oddechowy — **ZAMKNIETA** (2026-04-05, ex211)

**Cel:** Uogolnienie P8 na nierowne masy; zamknieta forma diagramu fazowego, predkosci ucieczki i efektywnego potencjalu oddechowego 3-cialowego.

**STATUS: ZAMKNIETA** — wyprowadzone i zweryfikowane numerycznie.

**Wyniki analityczne:**

1. **Rownowaga nierownych mas** (C1 != C2):
   d^2 - 4*beta*d + 9*gamma*(C1+C2) = 0
   d_rep,well = 2*beta -/+ sqrt(4*beta^2 - 9*gamma*(C1+C2))
   Max stosunek mas: q_max = 4*beta^2/(9*gamma*C1) - 1

2. **Diagram fazowy (granica istnienia):**
   beta_crit = (3/2)*sqrt(gamma*(C1+C2))
   C1+C2 < (4/9)*beta^2/gamma
   Przy beta_crit: d_rep = d_well = 2*beta_crit (bifurkacja saddle-node)
   Weryfikacja: 100% zgoda (skan 10x10)

3. **Predkosc ucieczki:**
   v_esc = sqrt(2*|V2(d_well)|/mu), mu = C1*C2/(C1+C2)
   Pelna forma zamknieta w (beta, gamma, C1, C2)
   Skalowanie: v_esc ~ beta^(-0.5) przyblizenie

4. **Mod oddechowy 3-cialowy:**
   U_eff(d) = 3*V2(d) + V3(d,d,d)
   omega^2_breathing = U_eff''(d_well)/C = (3/2)*omega^2_radial(2-body) + korekcja V3
   Korekcja V3: -1.3% (perturbacyjna)
   Bariera do wewnatrz: ~10x energia wiazania

**Pliki:**
- `examples/ex211_unequal_mass_phase_diagram.py` — 5 czesci + weryfikacja
- `tgp_nbody_results_clean.tex` — sekcja Result 8

---

## Harmonogram sugerowany

| Faza | Termin | Zakres |
|------|--------|--------|
| **I** | W toku | P1: Lyapunov — rdzeń numeryczny + ex148 (1 sesja); dalszy skan parametrów / 3B |
| **II** | Krotkoterminowo | P2: Hessian pelny (1 sesja) |
| **III** | ~~Srednioterminowo~~ **ZAMKNIETA** | P3: Poincare — ex203/ex204 |
| **IV** | Dlugoterminowo | P4: I_triple (badanie) |
| **V** | ~~Program~~ **ZAMKNIETA** | P5: Sciezka C — EFT projekcja, ex205 |
| **VI.A** | ~~Eksperymentalnie~~ **ZAMKNIETA** | P6.A: Oddzialywanie solitonow, ex206 |
| **VI.B** | Eksperymentalnie | P6.B: Dynamika zderzen (PDE) |
| **VII** | ~~Walidacja~~ **ZAMKNIETA** | P7: V3/V2 regime scan, ex209 |
| **VIII** | ~~Analityczne~~ **ZAMKNIETA** | P8: Analityczne rownowagi + Hill, ex210 |
| **IX** | ~~Analityczne~~ **ZAMKNIETA** | P9: Nierowne masy, diagram fazowy, v_esc, mod oddechowy, ex211 |

---

## Powazanie z reszta TGP_v1

- **dodatekJ/K** (solitony, ogon): nbody korzysta z profilu g(r) -> P6
- **dodatekV** (SU3): alpha_s z substratu -> niezalezne od nbody
- **B1/B2/B3** (mosty): zamkniete; nbody wzmacnia spojnosc przez niezalezne testy numeryczne
- **G1/G2** (luki): zamkniete sesja v42; nbody nie blokowane

*Plan do aktualizacji po kazdej zamknietej fazie.*
