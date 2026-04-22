# Synchronizacja `nbody` z rdzeniem TGP

Ten dokument jest kanonicznym punktem odniesienia dla warstwy `nbody/`.
Porządkuje trzy rzeczy:

1. co w `nbody` jest klasyczne,
2. co jest mostem EFT,
3. co jest już efektywną dynamiką wielociałową.

## 1. Trzy warstwy

### Warstwa A — klasyczny defekt

Punkt wyjścia stanowi soliton / defekt pola:

`g^2 g'' + g (g')^2 + (2/r) g^2 g' = g^2 (1-g)`

W tej warstwie ogon przy `g -> 1` jest **oscylacyjny**, nie Yukawowy:

`delta(r) = 1 - g(r) ~ [A cos(r) + B sin(r)] / r`

Wniosek:

- klasyczny rdzeń nie daje sam z siebie `exp(-m r)/r`,
- `sin(r)/r` i Yukawa nie mogą być traktowane jako to samo przybliżenie.

### Warstwa B — most EFT

Warstwa `nbody` używa efektywnej masy ekranowania:

`m_sp^2 = 3 gamma - 2 beta`

oraz projekcji klasycznego defektu na źródło Yukawy:

`C_eff = integral_0^inf delta(r) exp(-m_sp r) r dr`

To jest roboczy most używany w pakiecie:

`defekt klasyczny -> (C_eff, m_sp) -> źródło Yukawy`

Most jest zaimplementowany operacyjnie w:

- `bridge_nbody.py`
- `yukawa_from_defect.py`
- `examples/ex205_path_c_yukawa_from_defect.py`

### Warstwa C — efektywne EOM N-ciał

Po przejściu przez most EFT każda cząstka jest opisana przez:

- `C_i` jako efektywną siłę źródła,
- `m_i = C_i`,
- wspólne `beta`, `gamma`,
- `m_sp`.

Wtedy:

`delta_i(r) = C_i exp(-m_sp r) / r`

i dynamika używa:

`V = V_2 + V_3`

gdzie:

`V_2(d) = -4 pi C_i C_j / d + 8 pi beta C_i C_j / d^2 - 12 pi gamma C_i C_j (C_i + C_j) / d^3`

oraz

`V_3 = (2 beta - 6 gamma) C_i C_j C_k I(d_ij, d_ik, d_jk)`

Równania ruchu:

`C_i x_i'' = -grad_i V`

## 2. Co jest kanoniczne w `nbody`

W `nbody` przyjmujemy następujący podział:

- `CLASSICAL`: soliton, ogon `sin(r)/r`
- `EFT-DERIVED`: `C_eff`, `m_sp`
- `N-BODY`: Yukawa-source `V_2`, `V_3`, EOM

To rozdzielenie jest ważniejsze niż starszy podział na "Droga B" i "Droga C",
bo eliminuje mieszanie języka fundamentalnego z efektywnym.

## 3. Jednoznaczne konwencje równań

### Screening

W `nbody` obowiązuje:

`m_sp^2 = 3 gamma - 2 beta`

Jeżeli pojawia się pomocnicze `mu^2`, to nie jest ono tym samym obiektem co
`m_sp^2`. W pomocniczym obrazie stabilizacyjnym może występować relacja:

`mu^2 = 2 m_sp^2`

ale fizyczna masa ekranowania używana przez `nbody` to zawsze `m_sp`.

### Exact Yukawa-overlap scaling

Dokładna całka Yukawowa `I_Y(d12,d13,d23;m_sp)` zależy tylko od
bezwymiarowych kombinacji `t_ij = m_sp d_ij`, więc:

`I_Y(d12,d13,d23;m_sp) = F(m_sp d12, m_sp d13, m_sp d23)`

Stąd wynika ścisła tożsamość:

`(d12 ∂/∂d12 + d13 ∂/∂d13 + d23 ∂/∂d23 - m_sp ∂/∂m_sp) I_Y = 0`

W `nbody` jest to teraz jawnie spięte w `three_body_force_exact.py` przez
dokładne pochodne po `d_ij` i po `m_sp`.

Wygodna parametryzacja shape-space używana roboczo to:

`d_min <= d_mid <= d_max`

`q1 = d_min / d_max`

`q2 = d_mid / d_max`

`t = m_sp d_max`

wtedy:

`I_Y(d12,d13,d23;m_sp) = F(t; q1, q2)`

z domeną:

`0 < q1 <= q2 <= 1,   q1 + q2 >= 1`

To jest dziś najwygodniejszy język do map `V_3/V_2` i porównania geometrii.

Dla geometrii ogólnej duże-`t` jest kontrolowane przez minimum funkcji fazowej
z dokładnej reprezentacji Feynmana:

`lambda(q1,q2) = min_{alpha in Delta_2} sqrt(Q/Delta)`

tak że wykładnicza część asymptotyki ma postać:

`I_Y(t;q1,q2) ~ exp(-lambda(q1,q2) t) * [prefaktor i poprawki]`

W `nbody` jest to teraz jawnie dostępne przez helpery
`yukawa_phase_argument`, `yukawa_overlap_shape_rate`,
`yukawa_overlap_geometry_rate` w `three_body_force_exact.py`.

Dla geometrii równobocznej (`d12=d13=d23=d`, `t=m_sp d`) duże-`t` daje
kontrolowaną asymptotykę:

`I_Y^eq(t) ~ A_eq t^(-3/2) exp(-sqrt(3) t)`

`A_eq = 4 sqrt(2) pi^(3/2) / 3^(3/4)`

Pierwsza poprawka poddominująca ma postać:

`I_Y^eq(t) ~ A_eq t^(-3/2) exp(-sqrt(3) t) [1 - 5/(8 sqrt(3) t) + O(t^(-2))]`

czyli tłumienie jest szybsze niż stare heurystyki typu `exp(-3t/2)`.

### Irreducible three-body term

W pełnym efektywnym potencjale punktowych źródeł wkład 3-ciałowy pochodzi z
obu nieliniowości:

- z `Phi^3`: `+ 2 beta`
- z `Phi^4`: `- 6 gamma`

czyli razem:

`V_3 = (2 beta - 6 gamma) C_i C_j C_k I`

Skrócona forma `V_3 = -6 gamma C_i C_j C_k I` opisuje tylko sam wkład
kwartykowy i nie powinna być używana jako pełne równanie EOM.

## 4. Most do innych warstw

`nbody` nie zamyka pełnego łańcucha `substrat -> metryka -> wszystkie observables`.
Zamyka jednak most roboczy wystarczający do obliczeń warstw wyższych:

`g0, beta, gamma`
`-> g(r)`
`-> C_eff, m_sp`
`-> V_2, V_3, EOM`
`-> energie, mody, Lyapunov, bound-state observables`

To jest obecnie najuczciwszy opis roli pakietu.

## 5. Pliki źródłowe

- `tgp_field.py`: profile i masa ekranowania
- `bridge_nbody.py`: kanoniczny interfejs mostu
- `eom_tgp.py`: postać analityczna EOM
- `three_body_terms.py`: energia 3-ciałowa
- `three_body_force_exact.py`: dokładne siły i Jacobiany 3-ciałowe
- `tgp_nbody_lagrangian_eom.tex`: zapis LaTeX EOM
- `tgp_nbody_results_clean.tex`: tekst publikacyjny
- `examples/STATUS_MAP.md`: mapa `canonical` / `active-exploratory` / `legacy`
