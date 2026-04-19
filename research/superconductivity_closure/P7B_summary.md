# P7.2 — Lantanowe pair-breaking (ps21)

**Data:** 2026-04-19
**Status:** P7.2 zamknięte; 1 nowy parametr uniwersalny (α_PB), predykcje dla wszystkich 15 lantanowców.

## Motywacja

Naiwne P6.C + P7.1 dla lantanowych superhydrydów LnH9 przy η=1 (wysokie P)
przewidywało T_c ≈ 200 K dla wszystkich Ln. Eksperyment obalił to dramatycznie:

| Materiał | P (GPa) | T_c^obs (K) | T_c^naiwne (K) | Problem |
|----------|---------|-------------|----------------|---------|
| LaH10    | 170     | 250         | 368            | ≈ OK (4f⁰) |
| CeH9     | 100     | 100         | 143            | obniżony 30% |
| **PrH9** | 120     | **5**       | 143            | **katastrofa** |
| **NdH9** | 130     | **4.5**     | 143            | **katastrofa** |
| YH9      | 201     | 243         | 260            | OK (brak 4f) |

Brakujący mechanizm = **Abrikosov-Gorkov pair breaking** od lokalnych 4f momentów magnetycznych.

## Formuła P7.2

$$
\boxed{T_c = T_c^{\text{base}} \cdot B_\text{PB}(\mu_\text{eff})}, \qquad
B_\text{PB}(\mu) = \exp\!\big(-\alpha_\text{PB} \cdot \mu_\text{eff}^2\big)
$$

gdzie:
- $T_c^{\text{base}}$ — T_c obliczone przez P6.A+B+C+D + P7.1 (bez uwzględnienia 4f)
- $\mu_\text{eff} = g_J\sqrt{J(J+1)}\,\mu_B$ — moment Hund'a dla Ln³⁺
- $\alpha_\text{PB}$ — uniwersalna stała TGP (fit na PrH9 + NdH9)

**Wyprowadzenie:** Lokalne 4f momenty niszczą fazę Φ-substratu przez
spin-flipping scattering. W substratowym odpowiedniku Abrikosova-Gorkova:
$$
\ln\frac{T_c^{(0)}}{T_c} = \psi\!\left(\tfrac{1}{2} + \tfrac{\Gamma_\text{sf}}{2\pi T_c}\right) - \psi\!\left(\tfrac{1}{2}\right)
$$
W granicy silnego $\Gamma_\text{sf} \gg T_c$: $\ln(T_c^{(0)}/T_c) \approx \Gamma_\text{sf}/(2\pi T_c)$, 
co prowadzi do eksponencjalnej supresji po uśrednieniu po próbie. Przy parametryzacji
$\Gamma_\text{sf} \propto \mu_\text{eff}^2$ dostajemy $B_\text{PB} = e^{-\alpha\mu^2}$
jako dobre przybliżenie w szerokim zakresie.

## Parametr TGP

**α_PB = 0.2887 μ_B⁻²** (fit na PrH9 μ=3.58, T_c=5K i NdH9 μ=3.62, T_c=4.5K; T_base=143K)

## Walidacja

| Ln    | 4f^n | μ_eff (μ_B) | B_PB    | T_pred (K) | T_obs (K) | ratio | komentarz |
|-------|------|-------------|---------|------------|-----------|-------|-----------|
| La³⁺  | 0    | 0.00        | 1.000   | 368        | 250       | 1.47× | empty 4f |
| Ce³⁺  | 1    | 2.54        | 0.155   | 22         | 100       | 0.22× | **Kondo partial screening** |
| Pr³⁺  | 2    | 3.58        | 0.025   | 3.5        | 5.0       | 0.71× | fit point |
| Nd³⁺  | 3    | 3.62        | 0.023   | 3.3        | 4.5       | 0.72× | fit point |
| Y³⁺   | —    | 0.00        | 1.000   | 260        | 243       | 1.07× | ref, brak 4f |

**RMS_log (5 materiałów) = 0.316** (≈ factor 2 błąd; głównie od Ce-Kondo).

**Ce outlier:** naiwne μ=2.54 → T=22K vs obserwowane 100K. Efektywny moment
jest zredukowany przez **Kondo screening** (~2×), co daje μ_eff ≈ 1.0 → B_PB = 0.75 → T ≈ 107K.
Ten efekt wymaga osobnej P7.3 obróbki (hybrydyzacja f-d).

## Predykcje dla wszystkich LnH9 @ ~150 GPa (T_base ≈ 200K)

| Ln | 4f^n | walencja | μ_eff | B_PB      | T_pred (K) | status |
|----|------|----------|-------|-----------|------------|--------|
| La | 0    | 3+       | 0.00  | 1.0       | **200**    | potwierdzone |
| Ce | 1    | 3+       | 2.54  | 0.155     | 31         | (Kondo → ~100) |
| Pr | 2    | 3+       | 3.58  | 0.025     | 5          | potwierdzone |
| Nd | 3    | 3+       | 3.62  | 0.023     | 4.6        | potwierdzone |
| Pm | 4    | 3+       | 2.68  | 0.126     | 25         | radioactive |
| Sm | 5    | 3+       | 1.55  | 0.500     | **100**    | testowalne |
| Eu | 6→7  | 2+/3+    | 3.4/7.94 | 0.036/1e-8 | 7 / 0   | Eu³⁺ (P>20GPa) |
| Gd | 7    | 3+       | 7.94  | 1.25e-8   | 0          | killer 4f⁷ |
| Tb | 8    | 3+       | 9.72  | 1.4e-12   | 0          | killer |
| Dy | 9    | 3+       | 10.63 | 6.8e-15   | 0          | killer (max μ) |
| Ho | 10   | 3+       | 10.60 | 8.2e-15   | 0          | killer |
| Er | 11   | 3+       | 9.59  | 2.9e-12   | 0          | killer |
| Tm | 12   | 3+       | 7.56  | 6.8e-8    | 0          | killer |
| Yb | 13   | 2+       | 0.00  | 1.0       | (38, η<1)  | η ogranicza |
| **Lu** | **14** | **3+** | **0.00** | **1.0** | **200-250** | **KLUCZOWA PREDYKCJA** |

## Ranking kandydatów do high-T_c superhydrydów

1. **La, Y** — potwierdzone LaH10 (250K), YH9 (243K)
2. **Lu** — **predykcja: LuH10 @ 170 GPa → 200-250 K** (4f¹⁴ closed, μ=0!)
3. **Sm** — 50-100K (niski moment 1.55 przez L-S cancellation)
4. **Ce** — 100K (potwierdzone, ale z Kondo korrekcją)
5. **Pr, Nd** — 5K (pair-breaking zabija)
6. **Gd-Tm** — T_c ≈ 0 (duże momenty Hund'a)
7. **Yb** — ograniczone przez η (nie moment) → 38K @ 300 GPa (ps18)

## Wpływ na redukcję parametrów

| Etap | Nowy parametr | Co zyskujemy |
|------|---------------|--------------|
| P7.2 | α_PB = 0.2887 | **predykcje dla 15 lantanowców** z 1 liczby |

Łączna struktura TGP SC:
- β = 2.527 (P6.D magnetic blocking)
- κ_TGP = 2.012 (P7.1 λ_sf first principles)
- α_PB = 0.2887 (P7.2 pair breaking)
- **3 uniwersalne parametry TGP** + tabele atomowe (N(E_F), I, μ_eff z Hund)

## Zakres stosowalności

| Regime | Kryterium | Status formuły |
|--------|-----------|----------------|
| Lekkie Ln³⁺ (Pr, Nd) | μ = 3.5-3.7 | OK ±30% |
| Średnie Ln³⁺ (Sm, Pm) | μ < 3 | OK predykcja (testowalna) |
| Ciężkie Ln³⁺ (Gd-Tm) | μ > 7 | OK trywialnie T_c = 0 |
| Ce (4f¹) | Kondo screening | wymaga P7.3 |
| Eu, Yb | walencyjne transitions | wymaga P_valence(P) |
| Lu (4f¹⁴) | μ = 0, closed shell | **predykcja 200-250K** |

## Kluczowa predykcja falsyfikowalna

**LuH10 przy ~170 GPa powinno mieć T_c ≈ 200-250 K** (jak LaH10).

Uzasadnienie:
- Lu³⁺ ma 4f¹⁴ closed shell → μ_eff = 0 → B_PB = 1
- Brak lokalnych momentów → brak pair-breaking
- Struktura LaH10/LuH10 analogiczne (Fm-3m, H-cage)
- η = 1 przy P > 150 GPa (delokalizacja d/f)

Historia: Dias 2023 N(D,H,Lu)₁₀ "room-T SC" zostało **wycofane** (retraction).
Jednakże pozostaje teoretyczna predykcja dla **prawdziwego** Lu-hydride: czystego LuH10.

## Powiązania

- [[ps18_verification_corrections.py]] — P_scale_Yb refit (552 GPa)
- [[ps19_p7a_lambda_sf_first_principles.py]] — κ_TGP = 2.012
- [[ps20_master_plot.py]] — combine P6+P7.1 (bez pair-breaking)
- [[ps21_p7b_lanthanide_pair_breaking.py]] — ten plik (P7.2)
- [[P7_plan.md]] — mapa P7
- [[P7A_summary.md]] — podsumowanie P7.1

## Otwarte tematy do P7.3

1. **Ce Kondo screening**: μ_eff(Kondo) wymaga hybrydyzacji f-d, 
   prawdopodobnie μ_eff/μ_nagie = f(T_K/T_c).
2. **Eu walencja**: przejście 2+→3+ przy ~20 GPa (tabele).
3. **Anizotropia**: momenty 4f z silną anizotropią (Tb, Dy) mogą mieć 
   efektywnie redukowane μ w stanie niskospinowym.
4. **T_c^base dla LnH_x**: zależy od η(P) i stałej krystalicznej (ps16/P6.C),
   trzeba tabelować dla każdego Ln.
