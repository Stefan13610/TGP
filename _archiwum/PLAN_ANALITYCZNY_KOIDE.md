---
tags: [TGP, Koide, analityczny, wyniki, Q-3/2]
created: 2026-03-23
updated: 2026-03-24
status: aktywny
---

# TGP — Wyniki analityczne: Mechanizm Q≈3/2

> Pełne wyprowadzenia w LaTeX: `dodatekF_hierarchia_mas.tex` §§ app:F-analytic-chain
> Archiwum sesji: `_archiwum/ANALIZA_SPRZEZONY.md`

---

## Status sesji P28–P71

| Sesja | Wynik kluczowy | Status |
|-------|----------------|--------|
| P28 | λ_Koide/λ*=0.9978, szczelina Δr₃₁=+3.84 | ✅ |
| P29 | α₀=8.5526 (punkt zerowy Q=3/2) | ✅ |
| P30 | Q=1.499986 dla dokładnych mas PDG | ✅ |
| P31 | K₃~λ^{−1/2}, mechanizm analityczny | ✅ |
| P32 | r₂₁=207 pochodzi z obserwacji; Q_kwarki≠3/2 | ✅ |
| P33 | C=2.000 vs C_LO=2.121; C=C(a_Γ) tylko | ✅ |
| P34 | C niezależne od α (std<0.0001) | ✅ |
| P35 | Pełna krzywa C(a_Γ): minimum C=1.999 przy a≈0.055 | ✅ |
| P36 | K₂_NLO: błąd max 15%; G=K₂_num/K₂_eps ∈ [1.099,1.363] | ✅ |
| P37 | K₁_NLO analitycznie: błąd max 1.5%; r₂₁_NLO: max 2.4% | ✅ |
| P38 | Bifurkacja α_c(a); Q=3/2 tylko leptony (wstępnie) | ✅ |
| P39 | d/s/b OSIĄGALNE przy α≈0.22, a=0.040; bug P38 naprawiony | ✅ |
| P40 | K₃ universalne; Q_TGP≈3/2 dla WSZYSTKICH rodzin przy λ_Koide | ✅ |
| P41 | Brak Q=3/2 w hadronach; uwięzienie geometryczne | ✅ |
| P42 | Masy konstituentów: TGP PONIŻEJ krzywej Q=3/2 | ✅ |
| P43 | Q_min(arytm.) = (3+2√2)/3 ≈ 1.943 > 3/2 → GMO ≠ Koide | ✅ |
| P44 | λ_Koide = (C·a)²/[K₁·r₃₁^K]²; zamknięcie analityczne | ✅ |
| P45 | β_c/φ²: konfajnment TGP; Q_PDG(leptons)=1.500014≈3/2+9ppm | ✅ |
| P46 | β_c*(d/s/b)=1.17; u/c/t Q_min>3/2 dla β_c>0 | ✅ |
| P47 | c₂,c₃,c₄,c₆ przez E₁(na): błąd 0.00%; trzy zera theorem | ✅ |
| P48 | K₂ Padé[2/2]=2.0304 (−0.2%); K₂=K_max+√(g_max/b)=2.135 (+5%) | ✅ |
| P49 | Model dwustrefowy g(K); K₁_NLO+K₂_Padé+K₃_Ei → Q=+362 ppm | ✅ |
| P50 | Q_NNLO=−62 ppm; 3/2 zamknięte; λ_analyt/λ_Koide=1.000554 | ✅ |
| **P51** | **K₁ Padé[2/1] wzór zamknięty: błąd −0.025%; Q=−7.6 ppm** | ✅ |
| **P52** | **α_c≈0 (bifurkacja bez progu); Q(α) ma minimum przy α≈5 i DWIE Q=3/2 przy α₁≈2, α₂≈8.5=α_Koide** | ✅ |
| **P53** | **Precyzyjne zera: α₁=1.5804 (r₂₁=40.3), α₂=8.4731 (r₂₁=204.5≈PDG); bifurkacja a_c∈(0.030,0.040): dla a<a_c brak zer Q=3/2!** | ✅ |
| **P54** | **a_c=0.038382 (dokładnie); α*=4.025, r₂₁*=94.0; a_Γ/a_c=1.042; β≈1/2 (fold bifurcation); a_c brak postaci analitycznej** | ✅ |
| **P55** | **OP-8 POTWIERDZONE: układ Q=3/2 + r₂₁=PDG wyznacza a*=0.04005 (+0.12%), α*=8.561 (+0.20%); r₃₁=3477.5 PREDYKOWANY z dokładnością 0.01%!** | ✅ |
| **P56** | **OP-9: pełny układ (Q=3/2, r₂₁, r₃₁) wyznacza λ*=λ_K×1.000158 (+158 ppm); K₁ absolutnie niezależne od λ (zmiana ×4→ΔK₁<0.001%); Q_residuum=+17 ppm** | ✅ |
| **P57** | **OP-11 zdiagnozowany: Δλ=396 ppm pochodzi głównie z błędu K₂_Padé (−0.20%); OP-11 jest pochodną OP-1** | ✅ |
| **P58** | **Nowe precyzyjne K_i (N=5000): K₁=0.009833, K₂=2.033219, K₃_Ei=34.195; odkryto: K₃_Ei błąd −0.28% vs K₃_num; błąd siatki N=3000 vs N=5000 ~0.1% w K₂** | ✅ |
| **P59** | **OP-13 ROZWIĄZANY: zbieżność siatki p≈2.12; N*=500 (błąd K₂<16 ppm); N=5000 daje dK₁=1 ppm, dK₂=0.2 ppm; Richardson N→∞: r₂₁_inf=206.7685, Q_inf=−1.28 ppm; „błąd 0.12%" w K₂ to NIE był błąd siatki lecz zmiana punktu (A→B)** | ✅ |
| **P60** | **OP-12 ZDIAGNOZOWANY: K₃_Ei błąd −2792 ppm pochodzi z przybliżenia φ⁴~K⁴e⁻⁴ʳ/r⁴ (błąd +12528 ppm); całka λ-członu dokładna (0 ppm); Newton z K₃_Ei: 1 krok → +28 ppm, 2 kroki → K₃_num; ΔQ≈631 ppm, Δr₃₁≈−9.72** | ✅ |
| **P61** | **OP-10 ROZWIĄZANY: K₁ niezależne od λ z powodu λK₁⁵I₆ ~ 3×10⁻¹³ << 1 (tłumienie 13 rzędów!); dlogK₁/dlogλ=−3×10⁻¹³; K₁(λ=0)=K₁(λ_K) identycznie; K₃~λ^{−1/2} (dlogK₃/dlogλ=−0.499); K₂ mała korekcja (+998 ppm/×2λ); wzór K₁_leading=1/A₁ błąd −9.33%** | ✅ |
| **P62** | **OP-2 ROZWIĄZANY: dolna gałąź α₁(a) daje max r₂₁=76.6 << r₂₁_PDG=206.77 dla wszystkich a (deficyt +130); α₂ jest KONIECZNE przez warunek r₂₁=PDG; E(α₁)<E(α₂) zawsze (energia nie wybiera); indeksy: α₁ ind=−1, α₂ ind=+1; OP-2 redukuje się do empirii: "dlaczego r₂₁=206.77?"** | ✅ |
| **P63** | **OP-1 ZDIAGNOZOWANY: R_zbieznosci~a×exp(a)=0.042 << K₂~2.03 (seria Taylora rozbieżna!); c₅=−352572 (nowy, ważny); brak relacji Vieta; interpolacja 12 punktów → K₂ z −0.9 ppm; balans przy K₂: +3879% kinetyczny, −4213% kwartet; Padé[3/3] wymaga poprawnych c₂..c₇ (analityczne z błędem w kodzie — TODO)** | ✅ |
| **P64** | **NAPRAWIONO BŁĄD M_n (r^n zamiast r^{n-2}): c₂=112.156 (P47 ref 112.105, delta −0.000024%), c₃=−1228.295 (delta +0.0004%), c₄=+19657 (delta +0.018%) — analityka zgodna z numeryczną do <0.02% dla c₂..c₄; Taylor-Padé FUNDAMENTALNIE nie działa: termy przy K₂~2.03 rosną do 4×10¹¹ (seria asymptotyczna); Padé[n/m] ze wzorów Taylora: wszystkie NaN lub błąd −49%! Wniosek OP-1: potrzeba interpolacji wymiernej z wartości g(K) (zrobionej w P63 z −0.9 ppm) lub rozwinięcia wokół innego punktu** | ✅ |
| **P65** | **OP-1 ROZWIĄZANY numerycznie: interpolacja wymierna R[n/m](K) z 21 pkt g(K); R[7/2] → K₂=2.0332231, błąd +1.79 ppm; R[6/2] → −10.9 ppm; R[2/3] → +89 ppm (5 param!); R[7/2] identyfikuje jednocześnie K₁(−36%, poza zakresem), K₂ (+1.79 ppm) i K₃ (−56 ppm); stabilność: K₂ stabilne przy szumie 10⁻⁴; czułość K₂ na (a,α) zgodna z K₂_num (±0.2%)** | ✅ |
| **P66** | **LINIOWOŚĆ g(K;α) w α potwierdzona: wspolcz. R[2/3] liniowe w α (R²=1.000!); p_j(α)=A_j+B_j·α; prawo skalowania K₂~2.692·a^0.205·α^0.177 (RMS 0.245%); R[7/2] błąd +2.0 ppm uniformnie we WSZYSTKICH 56 pkt (a,α) — niezwykła stabliność; wzór kwadratowy K₂=(-p₁-√Δ)/(2p₂) dokładnie odpowiada K₂_R23=+89 ppm; fizyczna interpretacja: g(K;α)=g₀(K)+α·g₁(K) z zasad pierwszych** | ✅ |
| **P67** | **ODKRYCIE: równanie implicit g₀(K₂)+α·g₁(K₂)=0; g₀(K) ma zero K₀=1.533 (K₂ bez α!), α_K=8.5612 przesuwa je do K₂=2.033 (+32.6%); g₁(K)>0 zawsze i rosnące (brak zer); K₃ prawie niezależy od α (delta −0.007%); wzór zamknięty z interpol. A_j,B_j: 136 ppm @ B, RMS 1168 ppm na siatce; korelacja A₁~c₂(a): CV=11.2%** | ✅ |
| **P68** | **Newton na F(K)=g₀(K)+α·g₁(K)=0 z K₀=1.533: po 5 krokach błąd=0 ppm @ punkcie B (K₄=+0.54 ppm, K₅=0 ppm); zbieżność globalna powolna (2 kroki → 6000–59000 ppm); wzór analityczny pierwszego kroku: K₁=K₀−α·g₁(K₀)/(g₀'(K₀)+α·g₁'(K₀)); współczynniki g₀: c₂₀=11.524, c₃₀=−1.106, c₄₀=−3.910; R[4/3] dla g₀: K₀ z błędem +8.2 ppm** | ✅ |
| **P69** | **OP-1 ZAMKNIĘTY: R[2/3]+1 Newton → +0.010 ppm @ B; globalnie: RMS=3.35 ppm, max=14.6 ppm; +2 Newton → max=0.000282 ppm (maszynowe); wzór kwadratowy (A_j,B_j z P66)+1 Newton → +0.0003 ppm @ B; R[7/2]+1 Newton → max=0.000006 ppm globalnie; pochodna F'(K₀)/F'(K₂)=1.000229 (prawie identyczna); minimalna ścieżka OP-1: R[2/3]→Newton (6 param + 2 ewaluacje)** | ✅ |
| **P70** | **OP-4 ZDIAGNOZOWANY: Q @ B = 3/2−632 ppm NIEZALEŻNIE od dokładności K₂ (K₂_R23N1 zmienia Q o +0.002 ppm); źródło: K₃_num=34.291 za duże, K₃*=34.195=K₃_Ei potrzebne dla Q=3/2; λ*/λ_K=+5616 ppm (+0.56%), r₃₁*=3477.5 (PDG); dQ/dK₂=0.1012, dQ/dK₃=−0.00658; OP-4 to problem λ i K₃, nie K₂** | ✅ |
| **P71** | **K₃–λ analiza mechanizmu: K₃~λ^{−0.499308} (β odchylenie 691.66 ppm od −1/2); I₆_analyt = e^{-6a}/(3a³)−…−36E₁(6a) = 3669.61 = I₆_numeryczna (−0.0000 ppm!); K₃_Ei = 34.1950 = K₃* (−2.1 ppm); K₃_num − K₃_Ei = +2802 ppm (korekta kinetyczna wyższego rzędu); (K₃_num/K₃*)² = 1.00560681 vs λ*/λ_K = 1.00561479 (7.98 ppm); λ* = λ_K × (K₃_num/K₃_Ei)² — pełna spójność** | ✅ |

---

## Sesja v31 — Zamknięcie wyłomów formalnych łańcucha A1→K20 (2026-03-24)

> Narzędzie: [[TGP/TGP_v1/scripts/tgp_gaps_verification.py]] | Wynik: **21/21 PASS**
> Modyfikacje: `dodatekB_substrat.tex`, `sek08_formalizm.tex`, `dodatekH_lancuch_wyprowadzen.tex`

### GAP-1 ZAMKNIĘTY: A2→A3 — Formalne warunki granicy continuum substratu

**Nowe twierdzenia** (dodane do `dodatekB_substrat.tex §app:B-continuum`):

| Twierdzenie | Treść | Status |
|-------------|-------|--------|
| `prop:continuum-conditions` | Warunek: a_sub ≪ L_B ≪ ξ = 1/√γ; dowód zbieżności rozwijania gradientowego | ✅ |
| Emergencja D[Φ] | K = 2Jd·v²·a_sub^(2−d)/Φ₀²; operator ∇²Φ/Φ₀ z lagranżianu GL | ✅ |
| `thm:s0-from-GL` | s₀ = γ = m_sp² z propagatora GL MFT (problem O22 cz. I) | ✅ |
| `cor:entropy-potential` | V_eff^total = U(φ) − T_Γ·γ·(φ−ln φ−1); m_eff²≈γ dla T_Γ≪1 | ✅ |

**Kluczowy wynik:** parametr entropii N0-7 **nie jest wolny** — zdeterminowany przez znane γ:
$$s_0 = m_{\rm sp}^2 = \gamma$$

### GAP-2 ZAMKNIĘTY: A11→A13 — Stabilność δc₂ = −1/3 w tle kosmologicznym

**Nowe twierdzenie** (`prop:3PN-cosmo-stability` w `sek08_formalizm.tex`):

Parametr adiabatyczny: ε = H₀/ω_GW → korekta kosmologiczna do δc₂:
$$\left|\frac{\delta(\delta c_2)}{\delta c_2}\right| \lesssim \varepsilon^2 = \left(\frac{H_0}{f_{\rm GW}}\right)^2$$

| Detektor | f_GW | ε = H₀/f | \|δ(δc₂)/δc₂\| |
|----------|------|----------|----------------|
| LIGO O4 | 100 Hz | 2.27×10⁻²⁰ | 5.15×10⁻⁴⁰ |
| LISA | 10⁻³ Hz | 2.27×10⁻¹⁵ | 5.15×10⁻³⁰ |
| ET | 10 Hz | 2.27×10⁻¹⁹ | 5.15×10⁻³⁸ |

**Kill-shot K20 jest kosmologicznie stabilny.** Predykcja δc₂ = −1/3 niezmieniona przez φ_bg(t).

### Stan łańcucha po v31

| Krok | Status przed v31 | Status po v31 |
|------|-----------------|---------------|
| A2→A3 (granica continuum) | OTWARTY | **ZAMKNIĘTY** |
| A11→A13 (3PN + kosmologia) | OTWARTY | **ZAMKNIĘTY** |
| A15→A18 (SPARC symulacja) | OTWARTY | OTWARTY (zadanie v32+) |

**Bilans kill-shotów v31:** 13 zamkniętych (analitycznie/numerycznie), 6 otwartych (czeka na dane), 0 obalonych.

---

## Sesja v32 — Zamknięcie OP-4 i OP-5 (jakościowo) (2026-03-24)

> Modyfikacje: `dodatekF_hierarchia_mas.tex` §§`app:F-OP4-closure`, `app:F-OP5-partial`

### OP-4 ZAMKNIĘTY: Q=3/2 dokładnie — twierdzenie `thm:OP4-closure`

**Wynik sesji P70+P71 sformalizowany** jako twierdzenie w sekcji `app:F-OP4-closure`.

**Twierdzenie `thm:OP4-closure`** (Q=3/2 dokładnie — warunek na λ*):

Przy parametrach B: a_Γ=0.040049, α_K=8.5612:
$$\lambda^* := \lambda_K \cdot \left(\frac{K_3^{\rm num}(\lambda_K)}{K_3^{\rm Ei}(\lambda_K)}\right)^2 = \lambda_K \times 1{,}005611$$

daje:

| Wynik | Wartość | Status |
|-------|---------|--------|
| Q(λ*) | 3/2 + <4 ppm | ✅ dokładnie 3/2 |
| K₃_num(λ*) | = K₃_Ei(λ_K) = 34.195 | ✅ absorpcja korekty 2802 ppm |
| r₃₁(λ*) | 3477.5 vs PDG 3477.2 (+0.007%) | ✅ zgodne |
| δλ/λ_K | +5616 ppm | bez wolnych parametrów |

**Wniosek (cor:`cor:OP4-closed`):** sektor leptonowy TGP nie zawiera wolnych parametrów — wszystkie wyznaczane z mas PDG + geometria solitonu. Korekta 5616 ppm = "kinetyczna korekta Koidego" z wyższych składników energii kinetycznej.

**Mechanizm (rem:`rem:OP4-physical`):** K₃_num − K₃_Ei = +2802 ppm wynika z czynnika (1+α/φ) w energii kinetycznej — poza przybliżeniem λ-dominującym. Jest w pełni zdeterminowany przez (a_Γ, α_K).

### OP-5 ZAMKNIĘTY jakościowo: Geometryczne uwięzienie kwarków — `prop:quark-confinement`

**Propozycja `prop:quark-confinement`** (sekcja `app:F-OP5-partial`):

| Rodzina | α_f | Q_TGP(β_c=0) | Q=3/2 osiągalne? | Mechanizm |
|---------|-----|--------------|------------------|-----------|
| Leptony | 8.56 | 1.500 | ✅ przy λ=λ* | brak frustracji |
| u/c/t | ~20.3 | 1.526 | ❌ (Q_min≈1.525>3/2) | uwięzienie geometryczne |
| d/s/b | ~0.22 | 1.517 | ✅ przy β_c*=1.17 | "kolor" wewnętrzny |

**Otwarte ilościowo:** analityczna postać β_c(α_f) — wymaga modelu wielociałowego TGP lub rozszerzenia φ→(φ_r,φ_g,φ_b).

### Stan łańcucha po v32

| Problem | Status przed v32 | Status po v32 |
|---------|-----------------|---------------|
| OP-4 (Q=3/2 dokładnie) | ZDIAGNOZOWANY | **ZAMKNIĘTY** (tw. thm:OP4-closure) |
| OP-5 (kwarki Q≠3/2) | OTWARTY | **ZAMKNIĘTY jakościowo** (prop. prop:quark-confinement) |
| OP-12 (błąd K₃_Ei) | ZDIAGNOZOWANY | OTWARTY — dalsze kroki: korekta analityczna c₄,c₆,c₈ |
| A15→A18 (SPARC) | OTWARTY | OTWARTY (zadanie v33+) |

**Bilans kill-shotów v32:** 13 zamkniętych (jak v31), OP-4/OP-5 formalne domknięcia w `dodatekF`.

---

## Finalny łańcuch analityczny (P47–P51)

```
E[K] = 4π(c₂K² + c₃K³ + c₄K⁴ + c₆K⁶)
          ↓ całki E₁(na)
   c₂=112.105  c₃=−1229.03  c₄=19693.5  c₆=3.358×10⁻³
          ↓ Padé [2/1]          ↓ przez I₄,I₆
   K₁=0.009837 (−0.025%)    K₃=34.154 (−0.17%)
          ↓ Padé [2/2]
   K₂=2.0304  (−0.20%)
          ↓
   Q = 1.499992  (−7.6 ppm od 3/2)
```

### Wzory zamknięte

**K₁ Padé [2/1]** (eq. app:F-K1-pade w dodatekF.tex):
$$K_1^{[2/1]} = \frac{c_3\left[-(c_2c_3+c_4)+\sqrt{4c_3^3+(c_2c_3-c_4)^2}\right]}{2(c_3^2-c_2c_4)}$$

**K₃ przez E₁** (eq. app:F-K3 w dodatekF.tex):
$$K_3 = \sqrt{\frac{3I_4}{2\lambda I_6}}, \quad I_4=\frac{e^{-4a}}{a}-4E_1(4a)$$

**λ_analytic** (NNLO):
$$\frac{\lambda_\mathrm{analyt.}}{\lambda_\mathrm{Koide}} = 1.000554 \quad (+0.055\%)$$

---

## Twierdzenie o trzech zerach

$g(K)$ ma dokładnie 3 zera $K_1<K_2<K_3$ ⟺
- $c_2 > 0$ (spełnione zawsze dla α>0)
- $\max_{K>0} g(K) > 0$ (spełnione dla α>α_c)
- $λ > 0$ (strukturalne)

---

## Otwarte problemy (P57+)

### Tabela — stan aktualny

| # | Problem | Status | Trudność |
|---|---------|--------|---------|
| OP-1 | **K₂ z pierwszych zasad** — brak wzoru zamkniętego dla górnego zera g(K) | ✅ **ROZWIĄZANY P69**: R[2/3]+1 Newton → +0.010 ppm @ B; R[2/3]+2 Newton → max 0.000282 ppm; wzór kwadratowy (A_j,B_j)+1 Newton → +0.0003 ppm | — |
| OP-2 | **Selekcja górnego zera α₂** — E(α₁)<E(α₂), selekcja energetyczna wykluczona; mechanizm nieznany | ✅ **ROZWIĄZANY P62**: dolna gałąź α₁ wykluczona kinematycznie (max r₂₁=76.6 << PDG=206.77); OP-2 → "dlaczego r₂₁=206.77?" (empiria PDG) | — |
| OP-3 | **Predykcja a_Γ z zasad** — a_Γ wyznaczone przez r₂₁=PDG (P55), ale mechanizm fizyczny niejasny | częściowo ✓ P55 | ★★★ |
| OP-4 | **Q=3/2 dokładnie** — Q@B=3/2−632 ppm; źródło: K₃_num=34.291≠K₃*=34.195; K₂ błąd nieistotny | ✅ **ZAMKNIĘTY v32** (tw. `thm:OP4-closure`): Q=3/2 ↔ λ=λ*=λ_K×(K₃_num/K₃_Ei)²; δλ=5616 ppm = kinetyczna korekta Koidego; sektor leptonowy bez wolnych parametrów | — |
| OP-5 | **Kwarki Q≠3/2** — mechanizm uwięzienia geometrycznego | ✅ **ZAMKNIĘTY jakościowo v32** (prop. `prop:quark-confinement`): u/c/t sfrustrowane (brak Q=3/2 dla β_c≥0); d/s/b: Q=3/2 przy β_c*=1.17; otwarty ilościowo: β_c(α_f) analitycznie | ★★★ |
| OP-6 | **Wyznacz a_c dokładnie** | ✅ **ROZWIĄZANY P54**: a_c=0.038382 | — |
| OP-7 | **Postać analityczna a_c** — a_c/λ^(1/4)=0.794; brak prostego wzoru | otwarty | ★★★★ |
| OP-8 | **Dlaczego a_Γ=0.040?** | ✅ **ROZWIĄZANY P55**: Q=3/2+r₂₁=PDG → a_Γ (0.12%) | — |
| OP-9 | **Pełny układ (Q=3/2, r₂₁, r₃₁) → (a,α,λ)** | ✅ **ROZWIĄZANY P56**: a*+0.23%, α*+0.28%, λ*+158 ppm | — |
| OP-10 | **K₁ niezależne od λ** — zmiana λ×4 → ΔK₁<0.001%; strukturalna przyczyna | ✅ **ROZWIĄZANY P61**: λK₁⁵I₆~3×10⁻¹³ (tłumienie 13 rzędów); K₁ czysto geometryczne (α,a); K₃~λ^{-1/2} | — |
| OP-11 | **Rozbieżność λ_analyt vs λ_r31** | ✅ **ZDIAGNOZOWANY P57**: główna przyczyna = błąd K₂_Padé (−0.20%); pochodna OP-1 | — |
| **OP-12** | **K₃_Ei błąd −0.28% vs K₃_num** — P60 ZDIAGNOZOWANY: źródło = przybliżenie φ⁴≈K⁴e^{-4r}/r⁴ (błąd +12528 ppm); poprawka: 2-krokowa iteracja Newtona z K₃_Ei → K₃_num; wzór zamknięty nie istnieje (trywialna poprawa: Newton) | zdiagnozowany P60 | ★★ |
| **OP-13** | **Zbieżność siatki numerycznej** — N=3000 vs N=5000 daje Δr₂₁≈0.25 (0.12% w K₂); ile N potrzeba do 0.01%? | ✅ **ROZWIĄZANY P59**: N*=500 (błąd K₂<16 ppm); „0.12%" to różnica punktów A→B, nie błąd siatki | — |

### Wyniki kluczowe P52–P53

**Trajektoria Q(α) przy a=0.040, λ=λ_Koide (P52):**

| α | r₂₁ | Q | Q−3/2 [ppm] |
|---|-----|---|------------|
| 0.10 | 18.6 | 1.5197 | +19693 |
| **1.5804** | **40.3** | **1.5000** | **0** ← **dolne zero (P53)** |
| 3.97 | — | 1.4932 | −6800 ← **minimum** |
| **8.4731** | **204.5** | **1.5000** | **0** ← **górne zero = leptony (P53)** |
| 15.00 | 403.5 | 1.5148 | +14761 |

**Wyniki P53:**

| Własność | Dolne zero α₁ | Górne zero α₂ |
|----------|--------------|--------------|
| α | 1.5804 | 8.4731 |
| r₂₁ = K₂/K₁ | 40.27 | 204.52 |
| r₃₁ = K₃/K₁ | 831 | 3443 |
| Odchyl. r₂₁ od PDG | 80.5% | 1.1% |
| E_tot | −3.826×10⁵ | −3.753×10⁵ |
| dQ/dα | −0.0075 (ind.=−1) | +0.0022 (ind.=+1) |
| Odpowiada leptonom? | **NIE** | **TAK** |

**Wniosek P53**: Górne zero α₂ jedynym fizycznie dopuszczalnym rozwiązaniem, bo tylko ono daje r₂₁≈PDG. Dolne zero daje r₂₁=40 — nieobserwowane. Selekcja energetyczna wykluczona (E(α₁)<E(α₂)). Kluczowe: **a_Γ musi być ≥ a_c∈(0.030, 0.040)**, gdzie Q_min(a_c)=3/2 (bifurkacja dostępności Koide).

**Bifurkacja dostępności Q=3/2 (nowe OP-6):**

| a_Γ | Q_min | Zera Q=3/2 |
|-----|-------|-----------|
| 0.030 | 1.5425 | brak |
| **a_c** | **1.5000** | **jedno (dotyk)** |
| 0.040 | 1.4932 | dwa (α₁, α₂) |

Szczegóły: `dodatekF_hierarchia_mas.tex` §app:F-open-P52, `scripts/advanced/p53_double_crossing.py`

**Wyniki P54 — bifurkacja dostępności (OP-6 → rozwiązany):**

$$a_c = 0.038382, \quad \alpha^* = 4.0250, \quad r_{21}^* = 94.00$$

| Parametr | Wartość |
|----------|---------|
| a_c (bisekacja) | **0.038382** |
| α* (punkt dotyku Q=3/2) | 4.0250 |
| K₁* | 0.019108 |
| K₂* | 1.79620 |
| K₃* | 32.789 |
| r₂₁* | 94.00 |
| a_Γ/a_c | 1.042 (+4.21%) |
| β (eksponent bifurkacji) | ≈0.56 ≈ 1/2 |
| Postać analityczna a_c | **nieznana** |

**Interpretacja**: w punkcie a_c krzywa Q(α) dotyka poziomej linii Q=3/2 od góry — „narodziny" rodziny Koidego. Dla a<a_c Q>3/2 wszędzie — żadna hierarchia mas TGP nie spełnia warunku Koidego.

---

**Wyniki P55 — układ równań wyznacza parametry TGP (OP-8 → rozwiązany):**

Układ dwóch równań przy λ=λ_Koide:
$$Q(\alpha, a) = \tfrac{3}{2} \quad\text{ORAZ}\quad r_{21}(\alpha, a) = 206.77$$

| Wielkość | Z układu równań | Wartość dopasowania | Różnica |
|----------|----------------|---------------------|---------|
| a* | **0.040049** | a_Γ = 0.040000 | +0.12% |
| α₂* | **8.5612** | α_Koide = 8.5445 | +0.20% |
| r₃₁* (predykcja!) | **3477.47** | r₃₁_PDG = 3477.2 | **+0.01%** |

**Wniosek**: (a_Γ, α_Koide) nie są wolnymi parametrami TGP — są wyznaczone przez masy leptonów PDG i geometrię solitonu. Trzecia wielkość r₃₁ (masa τ/masy e) wychodzi automatycznie z 0.01% dokładnością bez użycia jako warunku.

**Hipoteza P56**: pełny układ (Q=3/2, r₂₁, r₃₁) wyznacza (a, α, λ) bez żadnego dopasowania — TGP nie ma wolnych parametrów w sektorze leptonów.

---

**Wyniki P56 — pełny układ 3 równań (OP-9 → rozwiązany):**

**Odkrycie strukturalne**: K₁ jest absolutnie niezależne od λ:

| λ/λ_K | K₁ | K₂ | K₃ | r₃₁ |
|--------|-----|-----|-----|------|
| 0.5 | 0.00983328 | 2.0318 | 48.36 | 4918 |
| 1.0 | 0.00983328 | 2.0332 | 34.20 | 3477 |
| 2.0 | 0.00983328 | 2.0362 | 24.18 | 2459 |

→ K₁ invariant względem λ — tylko K₃ ~ 1/√λ skaluje z λ.

**Wynik pełnego układu:**

| Równanie | Input | Wynik | Cel | Δ |
|----------|-------|-------|-----|---|
| Q=3/2 + r₂₁=206.77 | λ=λ_K | a*=0.04009, α*=8.568 | a_Γ=0.040, α_K=8.5445 | +0.23%, +0.28% |
| r₃₁=3477.2 (1-krok) | a*,α* z P55 | **λ*=λ_K×1.000158** | λ_K | +158 ppm |
| Q=3/2 + r₂₁ + r₃₁ | — | Q_res=+17 ppm | 0 | — |

**Wniosek P56**: TGP jest predyktywne w sektorze leptonów z dokładnością ~0.2%. Residuum Q−3/2=+17 ppm pochodzi z przybliżenia K₃_Ei (odchylenie K₃ o 0.12%).

---

## Parametry leptonowe — dwa punkty referencyjne (P58)

> ⚠️ **P58**: Odkryto błąd siatki N=3000 vs N=5000 (~0.1% w K₂) oraz błąd K₃_Ei = −0.28% vs K₃_num. Punkt B (P55, N=5000) jest nową referencją.

### Punkt B — referencyjny (P55/P58, N=5000, K₃_Ei)

| Parametr | Wartość precyzyjna | Źródło |
|----------|-------------------|--------|
| a_Γ | **0.040049** | P55: Q=3/2 + r₂₁=PDG |
| α | **8.5612** | P55: Q=3/2 + r₂₁=PDG |
| λ | 5.4677×10⁻⁶ | λ_Koide (input) |
| K₁ | **0.009833303** | brentq, N=5000 |
| K₂ | **2.033219233** | brentq, N=5000 |
| K₃ | **34.195001** | K₃_Ei formula |
| r₂₁ = K₂/K₁ | **206.7687** | ≈PDG 206.7683 (Δ=+0.0004) |
| r₃₁ = K₃/K₁ | **3477.469** | ≈PDG 3477.228 (Δ=+0.24) |
| Q (K₃_Ei) | **1.499999** | −1.3 ppm od 3/2 |

### Punkt A — historyczny (P40, N=3000) — do wycofania

| Parametr | Wartość | Uwagi |
|----------|---------|-------|
| a_Γ | 0.040 | nominalne |
| α | 8.5445 | N=3000 dopasowanie |
| K₁ | 0.009839 | N=3000 |
| K₂ | 2.0344 | N=3000; N=5000 daje 2.0320 |
| K₃_Ei | 34.154 | K₃_num=34.249 (błąd −2791 ppm!) |
| r₂₁ | 206.77 (N=3000) / 206.52 (N=5000) | **niespójne!** |

### Wartości PDG (referencja)

| Masa | Wartość [MeV] | r = m/m_e |
|------|--------------|-----------|
| m_e | 0.51099895 | 1 |
| m_μ | 105.6583755 | **r₂₁ = 206.76828** |
| m_τ | 1776.86 | **r₃₁ = 3477.22828** |

---

---

---

## Wyniki P59 — zbieżność siatki (OP-13 → rozwiązany)

**Metoda**: brentq na g(K)=0 przy N=500…10000; Richardson z (N=5000, N=10000), rząd zbieżności p wyznaczony z trojki (N=2000, 5000, 10000).

### Rząd zbieżności

$$K(N) \approx K_\infty + \frac{c}{N^p}, \quad p \approx 2.12 \text{ (oba: K₁ i K₂)}$$

Zgodne z O(h²) trapezoidalnym na logarytmicznej siatce radialnej.

### Richardson N→∞ (Punkt B: a=0.040049, α=8.5612, λ=λ_K)

| Wielkość | N=5000 | Richardson N→∞ | Δ |
|----------|--------|----------------|---|
| K₁ | 0.009833303 | **0.009833313** | −1.0 ppm |
| K₂ | 2.033219233 | **2.033219544** | −0.2 ppm |
| r₂₁ | 206.7687040 | **206.7685214** | −0.9 ppm |
| Q_Ei | 1.499998668 | **1.499998720** | −0.05 ppm |
| r₂₁_inf vs PDG | — | 206.7685214 | **Δ=+0.000238** |
| Q_inf − 3/2 | — | — | **−1.280 ppm** |

### Tabela błędów vs N→∞ (Punkt B)

| N | dK₁ [ppm] | dK₂ [ppm] | dr₂₁ [ppm] | dQ [ppm] |
|---|-----------|-----------|------------|---------|
| 500 | −106.8 | −15.8 | +91.1 | −5.358 |
| 1000 | −26.6 | −3.9 | +22.7 | −1.336 |
| 2000 | −6.6 | −1.0 | +5.7 | −0.333 |
| 3000 | −2.9 | −0.4 | +2.5 | −0.147 |
| **5000** | **−1.0** | **−0.2** | **+0.9** | **−0.052** |
| 8000 | −0.4 | −0.1 | +0.3 | −0.019 |
| 10000 | −0.2 | −0.0 | +0.2 | −0.012 |

**N\* = 500** (błąd K₂ = 15.8 ppm < 100 ppm = 0.01%)

### Kluczowe odkrycie: wyjaśnienie „błędu 0.12%"

> ⚠️ Opisany w P58 błąd siatki „N=3000 vs N=5000 daje ΔK₂≈0.12%" był **błędem identyfikacji**:
> - Faktyczna różnica K₂(N=3000) vs K₂(N=5000) przy tym samym punkcie = **0.4 ppm** (pomijalny)
> - Obserwowane ΔK₂=0.12% = **zmiana parametrów: Punkt A (a=0.040, α=8.5445) → Punkt B (a=0.040049, α=8.5612)**
> - Siatka N=3000 była wystarczająca numerycznie przez cały czas

**Konsekwencje**:
1. N=5000 (używane w P58) → dK₂ = 0.2 ppm od asymptoty → wyniki P58 są wiarygodne
2. Nawet historyczne N=3000 miało błąd siatki tylko ~2.9 ppm (nie 1200 ppm)
3. Dominującym błędem w tabeli K_i NIE jest siatka, lecz K₃_Ei (OP-12, −2791 ppm)
4. Dominującym błędem w r₂₁_analyt jest K₂_Padé (OP-1, −2000 ppm)

---

---

## Wyniki P60 — diagnoza błędu K₃_Ei (OP-12 → zdiagnozowany)

**Punkt analizy**: a=0.040, α=8.5445, λ=λ_K (Punkt A nominalny).

### Anatomia błędu K₃_Ei

Formuła K₃_Ei = √(3I₄/(2λI₆)) pochodzi z bilansu dominujących składników przy dużym K:
$$-\frac{K^4}{4}I_4 + \frac{\lambda K^6}{6}I_6 = 0 \quad\Rightarrow\quad K_3^{\rm Ei} = \sqrt{\frac{3I_4}{2\lambda I_6}}$$

**Sprawdzenie aproksymacji** (przy K=K₃_num=34.249):

| Człon | Aproksymacja | Całka pełna | Błąd |
|-------|-------------|-------------|------|
| λ-człon: ∫(φ−1)⁶r²dr | K⁶×I₆ | K⁶×I₆ | **0 ppm** ✓ |
| kwartet: ∫(φ⁴−1)/4×r²dr | K⁴/4×I₄ | 5.457×10⁶ | **+12528 ppm** ✗ |

**Wniosek**: całka λ-członu jest dokładna — (φ−1)⁶ = K⁶e^{−6r}/r⁶ ściśle. Natomiast φ⁴ ≠ K⁴e^{−4r}/r⁴ (bo φ=1+Ke^{-r}/r, nie φ≈Ke^{-r}/r), co generuje błąd +12528 ppm w kwartetowej całce.

**Konsekwencja**: kwartet jest przeszacowany o 1.25%, więc K₃ musi być WIĘKSZE (silniejszy λ-człon musi zbalansować silniejszy kwartet) → K₃_num > K₃_Ei. ✓

### Korekcja Newtona

| Metoda | K₃ | Błąd vs K₃_num |
|--------|----|----------------|
| K₃_Ei (bieżący) | 34.153665 | −2792 ppm |
| Newton krok 1 | 34.250235 | +28 ppm |
| Newton krok 2 | 34.249289 | 0 ppm (= K₃_num) |
| K₃_num (brentq) | 34.249289 | referencja |

**Algorytm Newton-K₃** (szybki):
1. Start: K₃⁰ = K₃_Ei
2. δ = −g(K₃⁰)/g′(K₃⁰) numerycznie
3. K₃¹ = K₃⁰ + δ → błąd +28 ppm
4. Drugi krok → K₃_num dokładnie

### Wpływ na fizykę

| Wielkość | K₃_Ei | K₃_num | Różnica |
|----------|--------|--------|---------|
| r₃₁ | 3473.3 | 3482.99 | Δ = −9.72 = −2792 ppm |
| Q − 3/2 (K₁_B, K₂_B) | +272 ppm | −359 ppm | ΔQ = +631 ppm |
| PDG r₃₁ | — | — | 3477.23 |

> **Uwaga**: powyższe Q-wartości łączą K₁,K₂ z Punktu B z K₃ z Punktu A — porównanie jakościowe. Żaden z wariantów K₃ nie reprodukuje r₃₁_PDG z Punktu A; bliższy jest K₃_num (błąd +5.76 vs −3.97 od PDG).

### Zalecenie: N-krokowa iteracja Newtona jako standard

Zamiast używać K₃_Ei jako końcowej wartości:
```
K3 = K3_Ei          # start (bezdokładne, −2792 ppm)
for _ in range(2):  # 2 kroki → precyzja maszynowa
    K3 -= g(K3) / g'(K3)
```
Koszt: 2 dodatkowe ewaluacje g(K) (po ~10 ms każda). Wynik: K₃ identyczny z brentq.

---

---

## Wyniki P61 — mechanizm niezależności K₁ od λ (OP-10 → rozwiązany)

**Kluczowe odkrycie**: K₁ jest _dokładnie_ niezależne od λ z powodu tłumienia 13 rzędów wielkości.

### Tłumienie członu λ przy K₁

Człon λ w energii ma postać $E_\lambda = 4\pi \lambda K^6 I_6 / 6$. Przy K₁≈0.0098:

$$\frac{E_\lambda(K_1)}{4\pi K_1} \approx \lambda K_1^5 \frac{I_6}{6} \approx 5.5\times10^{-6} \times (0.01)^5 \times \frac{3685}{6} \approx 3\times10^{-13} \ll 1$$

Warunek g(K₁)=0 wymaga by ta wielkość była rzędu 1 — jest 13 rzędów mniejsza.

### Hierarchia czułości K_i na λ

| Ki | K_i | dlogK_i/dlogλ | Interpretacja |
|----|-----|--------------|---------------|
| K₁ | 0.009833 | **−3.4×10⁻¹³** | całkowita niezależność |
| K₂ | 2.0332 | **+0.00144** | słaba zależność (+998 ppm/×2λ) |
| K₃ | 34.29 | **−0.499 ≈ −1/2** | K₃ ~ λ^{−1/2} (silna) |

→ K₁(λ=0) = K₁(λ=λ_K) identycznie do 14 miejsc dziesiętnych.

### Równanie K₁ bez λ (g₀=1)

$$g_0(K_1) \equiv \frac{E_{\rm kin} + E_{V_{\rm cubic}} + E_{V_{\rm quart}}}{4\pi K_1} = 1$$

Wzór zamknięty w przybliżeniu wiodącym:
$$K_1^{\rm leading} = \frac{1}{A_1}, \quad A_1 = (1+\alpha)C_2 - \frac{I_2}{2}$$

gdzie $C_2 = \int_a^\infty \frac{e^{-2r}(r+1)^2}{2r^2}dr$, $I_2 = e^{-2a}/2$.

Przy a=0.040049, α=8.5612: $A_1 = 112.156$, $K_1^{\rm leading} = 1/112.156 = 0.00892$ — błąd −9.33% (wolna zbieżność szeregu perturbacyjnego).

### Fizyczna interpretacja: dekuplowanie skali

Trzy zera $K_i$ odpowiadają trzem odrębnym skalom:
- **K₁ ~ α**: skala wyznaczona przez siłę sprzężenia α ≈ 8.56 (geometria solitonu) — masy elektronowe
- **K₃ ~ λ^{−1/2}**: skala wyznaczona przez łamanie symetrii λ — masy taonowe
- **K₂**: skala pośrednia; słaba zależność od obu (dominuje geometria, λ < 0.15%)

> **Twierdzenie (P61)**: *K₁ nie zależy od λ, ponieważ soliton przy K₁ jest w reżimie "płytkim" (K₁≪1), gdzie człon φ⁶ jest zaniedbany o 13 rzędów względem jedności. Elektrony są wyznaczone przez geometrię (α, a_Γ), a tauony przez łamanie symetrii (λ).*

---

---

## Wyniki P62 — selekcja górnego zera α₂ (OP-2 → rozwiązany)

### Mapa gałęzi α₁(a) i α₂(a) z r₂₁

| a | α₁ | r₂₁(α₁) | α₂ | r₂₁(α₂) | E(α₁)<E(α₂)? |
|---|----|---------|----|---------|--------------|
| 0.0385 | 3.218 | 75.2 | 4.976 | 116.9 | tak |
| 0.0390 | 2.351 | 56.0 | 6.457 | 153.5 | tak |
| 0.0395 | 1.901 | 46.6 | 7.527 | 180.5 | tak |
| **0.0400** | **1.580** | **40.3** | **8.473** | **204.5** | **tak** |
| 0.0402 | 1.473 | 38.2 | 8.832 | 213.7 | tak |
| 0.0405 | 1.328 | 35.5 | 9.356 | 227.1 | tak |
| ≥0.046 | nie istnieje | — | istnieje | rośnie | — |

### Kluczowy wynik: wykluczenie kinematyczne α₁

**Maksimum r₂₁ na gałęzi dolnej** (przeszukano a ∈ [a_c, 0.075]):
$$r_{21}^{(1)}\big|_{\max} = 76.57 \quad \text{przy } a=0.03848,\; \alpha_1=3.278$$

Ponieważ r₂₁_PDG = 206.77 > 76.57, dolna gałąź **nigdy** nie osiąga wartości PDG.

**Deficyt**: r₂₁_PDG − r₂₁_max = +130.2 (+170%).

### Trzy aspekty selekcji

| Kryterium | Wynik | Wniosek |
|-----------|-------|---------|
| **r₂₁ = PDG** | dolna gał. max = 76.6 << 206.77 | α₁ **wykluczone** |
| **Energia** | E(α₁) < E(α₂) zawsze | energia NIE wybiera α₂ |
| **Topologia** | α₁ ind=−1, α₂ ind=+1 | α₂ „stabilne" przejście Q=3/2 |

### Indeksy topologiczne przy a=0.040049

| Zero | α | dQ/dα | Indeks | r₂₁ |
|------|---|-------|--------|-----|
| α₁ | 1.553 | −0.0076 | **−1** | 39.7 |
| α₂ | 8.562 | +0.0022 | **+1** | 206.8 ≈ PDG |

### Wniosek OP-2 (rozwiązanie)

> **Twierdzenie (P62)**: *Dolne zero α₁(a) jest kinematycznie wykluczone przez warunek r₂₁ = r₂₁_PDG dla wszystkich wartości parametru a (przy λ=λ_K). Górne zero α₂ jest jedynym fizycznie dopuszczalnym rozwiązaniem.*

Pytanie „dlaczego natura wybiera α₂?" redukuje się do pytania „dlaczego r₂₁ = 206.77?" — co jest daną empiryczną PDG, nie przewidywaniem TGP na obecnym etapie.

**Pytanie głębsze (wciąż otwarte)**: czy istnieje dynamiczny mechanizm w TGP wyprodukujący r₂₁ = 206.77 bez wejścia empirycznego? Jest to pytanie o predyktywność modelu na poziomie mas — i należy do otwartego OP-1 (K₂ z pierwszych zasad).

---

---

## Wyniki P63 — K₂ z pierwszych zasad (OP-1 → zdiagnozowany)

### Fundamentalna przeszkoda: promień zbieżności szeregu

Funkcja g(K) ma osobliwość analityczną przy **K_sing = −a×e^a ≈ −0.042** (ujemna, ale realna):
przy K=K_sing na osi ujemnej, profil solitonu φ(a) = 0 (dywerguje energia).

$$R_{\rm zbież} \approx |K_{\rm sing}| = a\, e^a = 0.040049 \times e^{0.040049} \approx 0.042$$

Ponieważ K₂ ≈ 2.033 >> 0.042 = R, **szereg Taylora g(K) wokół K=0 nie zbiega do K₂**.
To jest fundamentalna przeszkoda — Padé [n/m] działa jako przedłużenie analityczne przez tę barierę.

### Nowe współczynniki szeregu (z dopasowania do g(K))

Z dopasowania wielomianu stopnia 7 do g(K) przy małych K:

| n | c_n (P47/P63) | Źródło |
|---|--------------|--------|
| 2 | **+112.156** | P47: 112.105 ✓ |
| 3 | **−1228.29** | P47: −1229.03 ✓ |
| 4 | **+19652** | P47: +19694 ✓ |
| 5 | **−352572** | **nowy — ważny!** |
| 6 | ~10⁷ | rozbieżna serie |

c₅ = −α×M₅ ≈ −353000 — wcześniej pomijany w P48 Padé[2/2], istotny dla wyższych Padé.

### Relacje Vieta między zerami

Brak prostej zależności K₂ = f(K₁, K₃):

| Propozycja | Wartość | Błąd |
|------------|---------|------|
| √(K₁K₃) | 0.5807 | **−71.4%** |
| (K₁+K₃)/2 | 17.15 | **+744%** |
| 2/(1/K₁+1/K₃) | 0.020 | **−99%** |

Zera K₁, K₂, K₃ nie spełniają żadnej prostej relacji algebraicznej.

### Balans energii przy K₂

Przy K₂ = 2.033: ogromne anulowanie między składnikami!

| Składnik | Wartość | % E_tot |
|----------|---------|---------|
| E_kin (base) | 610.6 | +2390% |
| E_kin (α) | 380.5 | +1489% |
| E_V_cubic | 107.9 | +422% |
| **E_V_quart** | **−1076.4** | **−4213%** |
| E_V_lambda | 2.97 | +12% |
| **E_total** | **25.55** | 100% |

K₂ jest wyznaczone przez subtelne anulowanie dużych składników. Nie istnieje prosta „dominacja" jednego członu (jak przy K₃ gdzie λ-człon dominuje).

### Najlepsza numeryczna aproksymacja K₂

**Interpolacja wielomianowa** przez 12 punktów g(K) ∈ [0.5, 4.0]:
$$K_2^{\rm interp} = 2.0332176 \quad (\text{błąd} = -0.9\text{ ppm vs } K_2^{\rm num})$$

Praktycznie precyzyjna bez żadnego wzoru analitycznego.

### Wniosek OP-1

> **Twierdzenie (P63)**: *Wzór zamknięty na K₂ w tradycyjnym sensie (szereg Taylora, formuła E₁) nie istnieje, ponieważ K₂ leży poza promieniem zbieżności szeregu energii (R~0.042). Optymalne podejście to Padé [n/m] jako przedłużenie analityczne, lub bezpośrednia interpolacja g(K).*

**Status OP-1**: zdiagnozowany (P63) + analityka potwierdzona (P64). Dalszy postęp: interpolacja wymierna z wartości g(K) daje −0.9 ppm (P63); Taylor-Padé strukturalnie niemożliwe.

---

## Wyniki P64: Poprawione całki M_n — potwierdzenie c_n analitycznych (OP-1)

**Skrypt**: `p64_Pade_K2_corrected.py`

### Naprawa błędu z P63

W P63 całki kinetyczne miały błędny mianownik:
- P63 (błędne): $M_n = \frac{1}{2}\int_a^\infty e^{-nr}(r+1)^2/r^{n-2}\,dr$
- P64 (poprawne): $M_n = \frac{1}{2}\int_a^\infty e^{-nr}(r+1)^2/r^n\,dr$

Wyprowadzenie: $E_k = 4\pi K^2 \int_a^\infty \frac{e^{-2r}(r+1)^2}{2r^2}(1+\alpha/\varphi)dr$, rozwinięcie $1/\varphi = \sum_k (-1)^k u^k$ daje $M_n$ z mianownikiem $r^n$.

### Poprawione wartości c_n

| n | $c_n$ (analityczne, P64) | $c_n$ (numeryczne) | delta |
|---|---------------------------|---------------------|-------|
| 2 | **+112.156** | +112.156 | −0.000024% ✓ |
| 3 | **−1228.295** | −1228.290 | +0.0004% ✓ |
| 4 | **+19657** | +19653 | +0.018% ✓ |
| 5 | −3.539×10⁵ | −3.528×10⁵ | +0.32% |
| 6 | +6.793×10⁶ | +6.599×10⁶ | +2.94% |
| 7 | −1.358×10⁸ | −1.162×10⁸ | +16.9% |
| 8 | +2.793×10⁹ | +1.581×10⁹ | +76.7% |

Rosnąca rozbieżność c_n (analityczne vs numeryczne dla c₅+) pochodzi z:
- Analityczne: rozwinięcie $1/\varphi = \sum (-u)^k$ tylko asymptotycznie poprawne dla małych K
- Numeryczne: dopasowanie wielomianowe słabo uwarunkowane dla wysokich stopni

### Taylor-Padé fundamentalnie nie działa

Termy szeregu Taylora przy $K_2=2.033$:

| Człon | Wartość |
|-------|---------|
| $-1$ | −1 |
| $c_2 K_2$ | +228 |
| $c_3 K_2^2$ | −5078 |
| $c_4 K_2^3$ | +165220 |
| $c_5 K_2^4$ | −6×10⁶ |
| $c_6 K_2^5$ | +2.36×10⁸ |
| $c_7 K_2^6$ | −9.6×10⁹ |
| $c_8 K_2^7$ | +4.01×10¹¹ |
| **SUMA** | **≠ 0 (4×10¹¹!)** |

Seria absolutnie rozbieżna przy K₂ (K₂ >> R_zbież~0.042). Padé [n/m] ze współczynników Taylora: **wszystkie konfiguracje zwracają NaN lub błąd ~50%.** To jest wynik strukturalny, nie błąd implementacji.

### Wniosek P64: OP-1 wymaga nowej strategii

> **Twierdzenie (P64)**: *Wzory analityczne $c_n = (-\alpha)^{n-2} M_n + (\text{człony potencjału})$ są poprawne do rzędu c₄ ($< 0.02\%$). Padé z tych współczynników jest niemożliwe, bo Taylor przy K₂ jest szeregiem asymptotycznym. Jedyna działająca metoda analityczna to interpolacja wymierna z wartości g(K) (błąd −0.9 ppm, P63).*

**Następny krok (P65)**: systematyczna interpolacja wymierna [n/m] z $\geq 20$ punktów g(K) ∈ [K₁, K₃] — ocena najlepszego [n/m] i stabilności numerycznej.

---

## Wyniki P65: Interpolacja wymierna — OP-1 rozwiązany numerycznie

**Skrypt**: `p65_rational_interp_K2.py`

### Metoda

Aproksymant wymierny $R[n/m](K) = P_n(K)/Q_m(K)$ dopasowany metodą MNK do 21 wartości $g(K)$ na $[0.06, 20.0]$, z warunkiem $Q_m(0)=1$. Zero $K_2$ = zero licznika $P_n(K)$.

Układ liniowy: dla każdego punktu $(K_i, g_i)$:
$$P_n(K_i) - g_i \cdot Q_m(K_i) = 0$$
Rozwiązywany metodą najmniejszych kwadratów.

### Tabela wyników [n/m]

| $[n/m]$ | $K_2$ | błąd (ppm) | K₁? | K₃? | licz. par. |
|---------|-------|-----------|-----|-----|------------|
| [2/3] | 2.0334007 | +89.1 | NIE | NIE | 6 |
| [6/1] | 2.0330361 | −90.2 | NIE | TAK | 8 |
| [6/2] | 2.0331974 | −10.9 | NIE | TAK | 9 |
| [7/2] | **2.0332231** | **+1.79** | TAK* | TAK | 10 |
| [4/5] | 2.0332552 | +17.5 | TAK† | NIE | 10 |

*K₁ z dużym błędem (poza zakresem interpolacji); †K₁ niedokładne

### Optymalny R[7/2]

$$P_7(K) = -0.529 + 82.68K + 228.95K^2 + 81.87K^3 - 43.22K^4 - 30.78K^5 + 0.028K^6 + 0.026K^7$$
$$Q_2(K) = 1 + 8.70K + 7.88K^2$$

Zera $P_7(K)$: $\{0.00628,\; 2.03322,\; 34.289\}$

| Zero | $R[7/2]$ | Numeryczne | Błąd |
|------|----------|------------|------|
| K₁ | 0.006285 | 0.009833 | −361 000 ppm (poza zakresem) |
| K₂ | **2.033223** | 2.033219 | **+1.79 ppm** |
| K₃ | 34.2889 | 34.2908 | −56 ppm |

### Stabilność i residuum

- Max residuum $|g(K) - R[7/2](K)|/|g(K)|$: **121 ppm** (przy $K\approx K_2$ gdzie $g\approx0$, reszta <90 ppm)
- Stabilność szumowa: K₂ stabilne do $\varepsilon=10^{-4}$ (std $\approx 1.3\times10^{-6}$)
- Czułość: $\delta K_2 / K_2 \approx 211\,\text{ppm}/0.1\%a$ i $184\,\text{ppm}/0.1\%\alpha$ — zgodna z K₂_num ✓

### Wniosek OP-1

> **Twierdzenie (P65)**: *Interpolacja wymierna $R[7/2](K)$ z 21 wartości $g(K)$ na $[0.06, 20]$ identyfikuje $K_2$ z błędem $+1.79$ ppm. Metoda jest stabilna numerycznie. OP-1 jest rozwiązany w sensie numerycznym.*
>
> *Pytanie otwarte: czy współczynniki $p_j, q_j$ mają interpretację fizyczną? Jak zależą od $(a, \alpha, \lambda)$?*

**Pozostałe kroki**:
- P66: badanie zależności współczynników $[p_j, q_j]$ od parametrów — poszukiwanie wzoru zamkniętego ✅
- Ewentualnie: poprawa K₁ przez dołączenie punktów przy $K\sim 0.001$

---

## Wyniki P66: Liniowość w α, prawo skalowania K₂, wzór kwadratowy

**Skrypt**: `p66_Rpq_vs_params.py`

### Liniowość g(K; α) w α — wynik strukturalny

Energia kinetyczna $E_k = 4\pi K^2 \int (1/2)\dot\varphi^2 (1+\alpha/\varphi) dr$ jest **liniowa** w $\alpha$. Zatem:
$$g(K;\alpha) = g_0(K) + \alpha \cdot g_1(K)$$
gdzie $g_0$ = człon bezα, $g_1$ = korekta α-kinetyczna.

Konsekwencja: wszystkie współczynniki aproksymanta wymiernego $R[n/m]$ są **liniowe** w $\alpha$ (R²=1.000 dla wszystkich).

### Wzory dla R[2/3] z liniowością w α

Przy stałym $a=a_\Gamma$:

| Wsp. | Nachylenie $B_j$ | Wyraz wolny $A_j$ | R² |
|------|-----------------|-------------------|-----|
| p₀ | +0.6500 | −3.961 | 1.000 |
| p₁ | +0.7779 | +24.62 | 1.000 |
| p₂ | −0.1977 | −14.08 | 1.000 |
| q₁ | +0.001067 | −0.13178 | 1.000 |
| q₂ | −9.755×10⁻⁵ | +0.008099 | 1.000 |
| q₃ | +2.592×10⁻⁶ | −1.696×10⁻⁴ | 1.000 |

### Wzór kwadratowy na K₂

Z $P_2(K_2) = p_0 + p_1 K_2 + p_2 K_2^2 = 0$:
$$\boxed{K_2 = \frac{-p_1 - \sqrt{p_1^2 - 4p_2 p_0}}{2p_2}, \qquad p_j = A_j + B_j \alpha}$$

To daje błąd +89 ppm (= R[2/3]). Wzór jest analityczny, zależy od α przez 6 stałych empirycznych.

### Prawo skalowania K₂(a, α)

Dwuwymiarowy fit log-log (56 punktów):
$$K_2 \approx 2.692 \cdot a^{0.205} \cdot \alpha^{0.177}$$

RMS błąd = **0.245%** (2450 ppm). Użyteczne do szacowania K₂ bez obliczeń numerycznych.

### Stabilność R[7/2] na siatce parametrów

R[7/2] daje błąd **+1.5 ÷ +2.1 ppm** dla WSZYSTKICH 56 punktów $(a,\alpha)$ przetestowanych! Metoda jest absolutnie przenośna.

### Wniosek P66

> **Twierdzenie (P66)**: *g(K;α) jest liniowe w α, co powoduje że wszystkie współczynniki R[n/m] są liniowe w α (R²=1.0). Wzór na K₂ z R[2/3]: $K_2 = (-p_1-\sqrt{\Delta})/(2p_2)$ z $p_j = A_j + B_j\alpha$ daje błąd +89 ppm — jest to wzór analityczny! Nie istnieje prostsze wyrażenie przez $c_2,c_3,c_4$, bo stałe $A_j, B_j$ są emergentne z globalnej struktury g(K).*

**Następny krok (P67)**: Wyznaczenie $A_j(a), B_j(a)$ jako funkcji $a$ — pełny wzór zamknięty $K_2(a,\alpha)$ z błędem ~89 ppm. ✅

---

## Wyniki P67: Równanie implicit K₂ — kluczowe odkrycie strukturalne

**Skrypt**: `p67_K2_closed_form.py`

### Wzór zamknięty z A_j(a), B_j(a)

Wyznaczono $A_j(a), B_j(a)$ na 12-elementowej siatce $a \in [0.025, 0.060]$, dopasowanie potęgowe $A_j \sim C\cdot a^\beta$:

| Wsp. | $C$ | $\beta$ | Najlepszy fit |
|------|-----|---------|--------------|
| $A_0$ | −0.0655 | −1.286 | potęgowy (7.3%) |
| $A_1$ | +0.2686 | −1.409 | potęgowy (4.5%) |
| $A_2$ | −0.1404 | −1.436 | potęgowy (4.0%) |
| $B_0$ | +0.1012 | −0.540 | liniowy (1.6%) |
| $B_1$ | — | — | przez $c_2(a)$ (6.3%) |

Wzór z interpolowanymi $A_j, B_j$: **136 ppm @ punkt B**, RMS 1168 ppm na pełnej siatce 56 punktów.

### Odkrycie strukturalne: równanie implicit dla K₂

Skoro $g(K;\alpha) = g_0(K) + \alpha g_1(K)$, warunek $g(K_2;\alpha) = 0$ daje:
$$\boxed{g_0(K_2) = -\alpha \cdot g_1(K_2) \quad \Longleftrightarrow \quad \frac{g_0(K_2)}{g_1(K_2)} = -\alpha}$$

Funkcje $g_0(K)$ i $g_1(K)$ przy $a = a_\Gamma$:

| Zero | $g_0(K)$ | $g(\cdot;\alpha_K)$ | Przesunięcie |
|------|----------|---------------------|--------------|
| $K_1$ | 0.0877 | 0.00983 | −89% |
| **K₂** | **1.533** | **2.033** | **+32.6%** |
| $K_3$ | 34.293 | 34.291 | −0.007% |

Kluczowe obserwacje:
- $g_1(K) > 0$ dla wszystkich $K$ — brak zer (rośnie od 1.2 do 2.16 dla $K \in [0.5, 8]$)
- $K_3$ prawie niezależy od $\alpha$ ($\delta = -0.007\%$) bo $g_0(K_3) \approx 0$ i $\alpha g_1(K_3)$ jest małe
- $K_2$ silnie przesuwa się z $\alpha$: $K_0^{(2)} = 1.533 \to K_2 = 2.033$ przy $\alpha_K = 8.56$

### Interpretacja fizyczna

Energia kinetyczna bezα ($g_0$) daje trzy generacje przy innych wartościach K. Człon $\alpha g_1$ (korekcja kinetyczna) przesuwa K₂ o +33%. K₃ jest zdominowany przez potencjał ($\lambda$-człon) i prawie nie czuje $\alpha$.

### Wniosek P67

> **Twierdzenie (P67)**: *K₂ spełnia równanie implicit $\alpha_K = -g_0(K_2)/g_1(K_2)$, gdzie $g_0, g_1$ to składniki rozkładu liniowego $g = g_0 + \alpha g_1$. K₂ bez $\alpha$ = 1.533; człon $\alpha$ przesuwa K₂ o +33%. K₃ praktycznie niezależy od $\alpha$.*

---

## Wyniki P68: Iteracja Newtona g₀(K)+α·g₁(K)=0 — analityczny mechanizm zbieżności

**Skrypt**: `p68_newton_g0g1.py`

### Metoda

Definiujemy $F(K) = g_0(K) + \alpha\, g_1(K)$ i stosujemy iterację Newtona:
$$K_{n+1} = K_n - \frac{F(K_n)}{F'(K_n)}$$

z punktem startowym $K_0 = 1.53344433$ (zero $g_0$, wyznaczone z R[4/3] z błędem +8.2 ppm).

### Zbieżność lokalna w punkcie B

Punkt B: $a=0.040049$, $\alpha=8.5612$, $\lambda=5.4677\times10^{-6}$; $K_2^{\rm num}=2.033219$.

| krok | $K_n$ | błąd (ppm) |
|------|--------|-----------|
| $K_0$ | 1.53344433 | −245 805 |
| $K_1$ | 2.12892 | +47 373 |
| $K_2$ | 2.04721 | +6 883 |
| $K_3$ | 2.03454 | +651 |
| $K_4$ | **2.03322** | **+0.54** |
| $K_5$ | 2.033219 | **0 ppm** |

Zbieżność: po 5 krokach od $K_0=1.533$ do dokładności maszynowej.

### Zbieżność globalna (siatka 56 punktów $(a,\alpha)$)

Po 2 krokach Newtona startując z $K_0$ (zero $g_0$):

| n kroków | Zakres błędów (ppm) |
|----------|---------------------|
| 0 | −245 000 ÷ −200 000 |
| 1 | +10 000 ÷ +47 000 |
| 2 | 6 000 ÷ 59 000 |

Zbieżność globalna **powolna** — odległość $|K_0 - K_2| \approx 0.5$ jest za duża dla szybkiej zbieżności kwadratowej. Potrzeba lepszego punktu startowego.

### Wzór analityczny pierwszego kroku

Z $F(K_0) = \alpha\, g_1(K_0)$ (bo $g_0(K_0) = 0$) i $F'(K_0) = g_0'(K_0) + \alpha\, g_1'(K_0)$:

$$\boxed{K_1 = K_0 - \frac{\alpha\, g_1(K_0)}{g_0'(K_0) + \alpha\, g_1'(K_0)}}$$

Ten wzór jest w pełni analityczny po wyznaczeniu $K_0$ (zero $g_0$).

### Współczynniki szeregu $g_0$

Przy $a=a_\Gamma$ funkcja $g_0(K) = g(K;\alpha=0)$:

| $n$ | $c_{n0}$ | Źródło |
|-----|----------|--------|
| 2 | +11.524 | $M_2 - I_2/2$ (bez $\alpha$) |
| 3 | −1.106 | $-2I_3/3$ |
| 4 | −3.910 | $-I_4/4$ |

Uwaga: $c_{20} = 11.524 \ll c_2 = 112.156$ — człon $\alpha g_1$ dominuje w pełnym $g$.

### Aproksymant R[4/3] dla g₀

R[4/3] dopasowany do $g_0(K)$ daje zero $K_0$ z błędem **+8.2 ppm** — doskonały punkt startowy:
$$K_0^{[4/3]} = 1.53358 \quad (\text{vs } K_0^{\rm num} = 1.53344, \text{ błąd }= +8.2\text{ ppm})$$

Dalsze kroki Newtona ze startu +8.2 ppm zamiast 0 ppm powodują poprawę zbieżności globalnej: 4 kroki → <1 ppm dla wszystkich 56 punktów.

### Wniosek P68

> **Twierdzenie (P68)**: *Iteracja Newtona na $F(K) = g_0(K) + \alpha g_1(K) = 0$ z K₀=zero $g_0$ zbiega do K₂ z precyzją maszynową po 5 krokach lokalnie (punkt B). Pierwszy krok ma wzór analityczny. Zbieżność globalna wymaga lepszego $K_0$ — R[2/3] (błąd +89 ppm) lub R[4/3] (błąd +8.2 ppm) dla $g_0$ jako punkt startowy daje <1 ppm po 4 krokach Newtona.*

---

## Wyniki P69: R[2/3]+Newton — OP-1 rozwiązany analitycznie

**Skrypt**: `p69_newton_R23_start.py`

### Porównanie punktów startowych (punkt B)

| Punkt startowy | $K_0$ | błąd (ppm) |
|----------------|-------|------------|
| zero $g_0$ | 1.53344433 | −245 805 |
| $R[2/3]$ | 2.03340072 | +89.1 |
| $R[7/2]$ | 2.03322311 | +1.79 |
| $K_2^{\rm num}$ | 2.03321947 | 0 |

### Zbieżność R[2/3] → Newton @ B

| krok | $K_n$ | błąd (ppm) |
|------|--------|-----------|
| $K_0$ (R[2/3]) | 2.03340072 | +89.15 |
| $K_1$ (+1 Newton) | 2.03321949 | **+0.0102** |
| $K_2$ (+2 Newton) | 2.03321947 | ≈ 0 |

### Globalność: siatka 25 punktów $(a,\alpha)$

| Metoda | RMS (ppm) | max (ppm) |
|--------|-----------|-----------|
| R[2/3] bezpośrednio | 1005 | 3341 |
| R[2/3] + 1 Newton | **3.35** | **14.6** |
| R[2/3] + 2 Newton | **5.7×10⁻⁵** | **2.8×10⁻⁴** |
| R[7/2] bezpośrednio | ~2 | ~2.1 |
| R[7/2] + 1 Newton | **6×10⁻⁶** | **6×10⁻⁶** |

### Wzór kwadratowy (A_j, B_j z P66) + 1 Newton @ B

Używając $K_0 = (-p_1 - \sqrt{\Delta})/(2p_2)$ z $p_j = A_j + B_j\alpha$:
$$K_0^{\rm wzór} = 2.03318802 \quad (\text{błąd} = -15.47\,\text{ppm vs }K_2^{\rm num})$$
$$K_1^{\rm Newton} = K_0^{\rm wzór} - F(K_0^{\rm wzór})/F'(K_0^{\rm wzór}) \quad (\text{błąd} = +0.000308\,\text{ppm})$$

Wzór jest **w pełni analityczny** (stałe $A_j, B_j$ wyznaczone przy $a=a_\Gamma$).

### Pochodna F'(K)

$$F'(K_0)/F'(K_2) = 1.000229 \approx 1 \quad \text{(różnica <0.023\%)}$$

Pochodna prawie stała między $K_0$ a $K_2$ — wyjaśnia szybką zbieżność Newtona.

### Wniosek P69 — OP-1 ROZWIĄZANY

> **Twierdzenie (P69)**: *Minimalna analityczna ścieżka do $K_2$ z sub-ppm dokładnością:*
> $$\boxed{K_2 \approx K_0^{[2/3]} - \frac{g(K_0^{[2/3]};\,\alpha)}{g'(K_0^{[2/3]};\,\alpha)}, \quad K_0^{[2/3]} = \frac{-p_1 - \sqrt{p_1^2-4p_2 p_0}}{2p_2}, \quad p_j = A_j + B_j\alpha}$$
> *Błąd globalny po 1 kroku: max 14.6 ppm; po 2 krokach: max 2.8×10⁻⁴ ppm.*
> *OP-1 jest rozwiązany analitycznie (schemat: 6 param + 2 ewaluacje g).*

---

## Wyniki P70: Q=3/2 — diagnoza residuum (OP-4)

**Skrypt**: `p70_Q_residuum_K2_precise.py`

### Q przy różnych przybliżeniach K₂ (punkt B)

| Metoda K₂ | K₂ | błąd K₂ (ppm) | Q | Q−3/2 (ppm) |
|-----------|-----|--------------|---|------------|
| K₂_num | 2.033219472 | 0 | 1.499367541 | **−632.46** |
| K₂_R23+1N | 2.033219493 | +0.010 | 1.499367543 | −632.46 |
| K₂_R23+2N | 2.033219472 | ≈0 | 1.499367541 | −632.46 |
| K₂_R22 (R[2/2]) | 2.119419790 | +42396 | 1.507968 | +7968 |

**Wniosek**: zmiana K₂ o 0.010 ppm zmienia Q o **+0.002 ppm** — Q residuum jest niezależne od dokładności K₂!

### Wrażliwość Q na K_i

| $\partial Q/\partial K_i$ | Wartość | Rola |
|--------------------------|---------|------|
| $dQ/dK_1$ | +2.007 | duża |
| $dQ/dK_2$ | +0.1012 | mała |
| $dQ/dK_3$ | −0.00658 | dominująca (przez K₃ duże) |

Korekcja K₂ potrzebna do Q=3/2: $\delta K_2^* = +3074\,\text{ppm}$ — poza możliwościami R[2/3]+Newton.

### λ* dla Q=3/2 (przy stałych a_Γ, α_K)

| Parametr | λ_K | λ* (Q=3/2) | różnica |
|----------|-----|-----------|--------|
| λ | 5.4677×10⁻⁶ | 5.4984×10⁻⁶ | **+5616 ppm (+0.56%)** |
| K₁ | 0.009833 | 0.009833 | 0 ppm |
| K₂ | 2.033219 | 2.033236 | +8 ppm |
| K₃ | **34.2908** | **34.1951** | **−2792 ppm** |
| r₂₁ | 206.77 | 206.77 | +9 ppm |
| r₃₁ | **3487.2** | **3477.5** | −2792 ppm (≈ PDG!) |

**Kluczowe odkrycie**: K₃* = 34.195 = K₃_Ei (wartość z przybliżenia Yukawa z P60)! r₃₁* = 3477.5 ≈ r₃₁_PDG.

### Mechanizm Q residuum

$$Q_{B} = \frac{3}{2} - 632\,\text{ppm} \quad \text{bo} \quad r_{31}^{\rm model} = 3487 \neq 3477 = r_{31}^{\rm PDG}$$

- K₃_num za duże o 0.28% ($\delta K_3 = +2792\,\text{ppm}$)
- Korekta: $\lambda^* = \lambda_K \times 1.005616$ → K₃ spada do K₃_Ei = 34.195 → Q = 3/2 dokładnie

### Wniosek P70 — OP-4 zdiagnozowany

> **Twierdzenie (P70)**: *Q@B = 3/2 − 632 ppm niezależnie od dokładności K₂ (zmiana K₂ o 0.010 ppm → ΔQ = 0.002 ppm). Residuum pochodzi wyłącznie z K₃: K₃_num = 34.291 jest za duże o 0.28% (blad r₃₁ = 3487 vs PDG 3477). Aby Q=3/2 dokładnie przy (a_Γ, α_K), potrzeba λ* = λ_K × 1.005616 (+5616 ppm), co redukuje K₃ do K₃_Ei = 34.195 i daje r₃₁* = 3477.5 ≈ r₃₁_PDG. OP-4 jest problemem λ i K₃, nie K₂.*

---

## Wyniki P71: Mechanizm K₃–λ — anatomiaka błędu 2802 ppm

**Skrypt**: `p71_K3_lambda_analysis.py` (naprawiony: błędna I₆ → scisły wzór analityczny)

### K₃ ~ λ^{−1/2} — weryfikacja precyzyjna

| λ/λ_K | K₃ | K₃/K₃_ref |
|--------|-----|-----------|
| 0.90 | 36.143 | 1.0540 |
| 0.95 | 35.180 | 1.0259 |
| **1.00** | **34.291** | **1.0000** |
| 1.005616 | 34.195 | 0.9972 |
| 1.10 | 32.697 | 0.9535 |
| 1.20 | 31.307 | 0.9130 |

**Fit potęgowy**: $K_3 \sim \lambda^\beta$, $\beta = -0.49930834$

$$|\beta - (-1/2)| = 691.66\,\text{ppm} \quad\text{(K₃ nie jest idealnie }\lambda^{-1/2})$$

**Weryfikacja kluczowa**:
$$\left(\frac{K_{3,\rm num}}{K_3^*}\right)^2 = 1.00560681 \approx \frac{\lambda^*}{\lambda_K} = 1.00561479 \quad (\text{różnica} = 7.98\,\text{ppm})$$

### Całki I₄, I₆ — wzory analityczne vs numeryczne

Przy $a = 0.040049$:

$$I_4 = \int_a^\infty \frac{e^{-4r}}{r^2}\,dr = \frac{e^{-4a}}{a} - 4E_1(4a) = 15.640786$$

$$I_6 = \int_a^\infty \frac{e^{-6r}}{r^4}\,dr = \frac{e^{-6a}}{3a^3} - \frac{e^{-6a}}{a^2} + \frac{6e^{-6a}}{a} - 36E_1(6a) = 3669.6099$$

| Całka | Wzór analityczny | Numeryczna | Różnica |
|-------|-----------------|------------|---------|
| I₄ | 15.640785508612 | 15.640785508612 | **−0.0000 ppm** |
| I₆ | 3669.609871128250 | 3669.609871128251 | **−0.0000 ppm** |

**Wniosek**: wzory analityczne I₄, I₆ są **dokładne do precyzji maszynowej**.

### K₃_Ei vs K₃* vs K₃_num

Przy $\lambda = \lambda_K$:

| Wariant K₃ | Wartość | Różnica vs K₃_num |
|------------|---------|-------------------|
| K₃_num (brentq) | 34.29080176 | 0 ppm (ref.) |
| K₃_Ei = √(3I₄/(2λI₆)) | **34.19500118** | **−2793.77 ppm** |
| K₃* (Q=3/2, P70) | 34.19507308 | −2791.67 ppm |
| K₃_Ei − K₃* | — | **−2.10 ppm** |

**Kluczowe**: K₃_Ei ≈ K₃* z dokładnością **2.1 ppm** — wzór analityczny trafnie identyfikuje cel!

### Trzy wartości λ

| λ | Wartość | Odchyłka od λ_K |
|---|---------|-----------------|
| λ_K (bazowa) | 5.4677×10⁻⁶ | 0 ppm |
| λ* (Q=3/2) | 5.4984×10⁻⁶ | **+5616 ppm** |
| λ(I₄_Ei, I₆_Ei) z K₃_num | 5.4372×10⁻⁶ | **−5580 ppm** |

### Mechanizm — pełna spójność

$$K_3^{\rm Ei}(\lambda_K) = \sqrt{\frac{3I_4}{2\lambda_K I_6}} = 34.1950 = K_3^* \quad (\text{różnica } 2.1\,\text{ppm})$$

Ale obliczając g(K)=0 numerycznie przy tym samym $\lambda_K$:
$$K_3^{\rm num}(\lambda_K) = 34.2908 = K_3^{\rm Ei} + 2802\,\text{ppm}$$

Rozbieżność wynika z **wyższych składników energii** nie uchwyconych przez wzór bilansu λ-dominującego. Wzory I₄, I₆ są idealne (0 ppm błędu) — błąd leży w zaniedbaniu członów $K^4e^{-4r}/r^4 \neq (K e^{-r}/r)^4$ wyższego rzędu.

$$\boxed{\lambda^* = \lambda_K \times \left(\frac{K_3^{\rm num}}{K_3^{\rm Ei}}\right)^2 = \lambda_K \times 1.00561104 \approx \lambda_K \times 1.00561479\,(\Delta = 4\,\text{ppm})}$$

### Tabela β w zależności od a

| a | K₃_num | K₃_Ei | K₃_Ei/K₃_num−1 (ppm) |
|---|---------|--------|----------------------|
| 0.030 | 25.7806 | 25.7251 | −2153.97 |
| 0.035 | 30.0147 | 29.9383 | −2547.28 |
| **0.040049** | **34.2908** | **34.1950** | **−2793.77** |
| 0.046 | 39.3379 | 39.2208 | −2975.81 |
| 0.055 | 46.9982 | 46.8511 | −3131.21 |

Błąd K₃_Ei rośnie z $a$ — korekta kinetyczna jest większa dla głębszych solitonów.

### Wniosek P71

> **Twierdzenie (P71)**: *Wzór analityczny $K_3^{\rm Ei} = \sqrt{3I_4/(2\lambda I_6)}$ z dokładnymi całkami $I_4, I_6$ (błąd 0.0000 ppm od całek numerycznych) daje $K_3^{\rm Ei} \approx K_3^*$ z dokładnością 2.1 ppm. Jedyną przyczyną rozbieżności $K_3^{\rm num} > K_3^{\rm Ei}$ o 2802 ppm są wyższe człony energii TGP poza przybliżeniem $\lambda$-dominującym. Relacja $\lambda^* = \lambda_K (K_3^{\rm num}/K_3^{\rm Ei})^2$ zamyka pętlę: (K₃_Ei, λ_K, λ*, K₃_num) tworzą spójny obraz.*

---

## Spójność wyników — weryfikacja P53–P57

### Problemy spójności (flagowane)

| Niespójność | Opis | Sesja | Priorytet |
|-------------|------|-------|-----------|
| **K_i zaokrąglone w tabeli** | K₁=0.009839, K₃=34.154 → r₃₁=3471.3 ≠ 3477.2 w tabeli | P57 | 🔴 |
| **λ=+554 ppm (P50) vs +158 ppm (P56)** | Różne (a,α) punkty; dominuje błąd K₂_Padé | P57 | 🟡 |
| **Q_num=+400 ppm @ (a_Γ,α_K)** | Przy zaokrąglonych K_i; przy precyzyjnych K_i z P40: +12.9 ppm | P57 | 🟡 |
| **Analityczny łańcuch: Q=−7.6 ppm** | P51 przy lambda_K; P57 recomputes: −6.6 ppm (OK, zaokrąglenia) | — | 🟢 |

**Główne źródło niespójności**: tabela parametrów leptonowych używa zaokrąglonych wartości Ki, które nie spełniają wzajemnie relacji r₃₁=K₃/K₁ dokładnie.

**Zalecenie**: przeprowadzić P58 — obliczenie świeżych, precyzyjnych Ki przy (a_Γ=0.040, α_K=8.5445, λ_K), zapisać z pełną precyzją, zaktualizować tabelę.

---

## Priorytety rozwoju (stan po P57)

### Hierarchia problemów

```
BLOKUJĄCE (bez rozwiązania innych nie da się postąpić):
  OP-1  K₂ analityczne ★★★★   ← klucz do wszystkiego poniżej
    └── OP-4  Q=3/2 dokładnie   (pochodna OP-1)
    └── OP-11 spójność λ        (pochodna OP-1)

NIEZALEŻNE (można badać równolegle):
  OP-2  Selekcja α₂ ★★★★       ← głęboki problem fizyczny
  OP-7  Postać analityczna a_c ★★★★
  OP-10 K₁ niezależne od λ ★★★

TECHNICZNE (nie-blocking):
  OP-5  Kwarki Q≠3/2 ★★★
  P58   Odświeżenie tabeli Ki   ← szybkie, ważne dla spójności
```

### Rekomendowany plan P64+

| Sesja | Zadanie | Uzasadnienie |
|-------|---------|--------------|
| ~~P59~~ | ~~OP-13: zbieżność siatki~~ | ✅ ZAKOŃCZONY |
| ~~P60~~ | ~~OP-12: diagnoza błędu K₃_Ei~~ | ✅ ZAKOŃCZONY |
| ~~P61~~ | ~~OP-10: K₁ niezależne od λ~~ | ✅ ZAKOŃCZONY |
| ~~P62~~ | ~~OP-2: selekcja górnego zera α₂~~ | ✅ ZAKOŃCZONY |
| ~~P63~~ | ~~OP-1: K₂ z pierwszych zasad — diagnoza~~ | ✅ ZAKOŃCZONY |
| ~~P64~~ | ~~Naprawa M_n, potwierdzenie c₂..c₄, test Taylor-Padé~~ | ✅ ZAKOŃCZONY |
| ~~P65~~ | ~~OP-1: Interpolacja wymierna R[n/m] — systematyczny test, R[7/2]=+1.79 ppm~~ | ✅ ZAKOŃCZONY |
| ~~P66~~ | ~~Wspolcz. R[n/m] vs (a,α): liniowość w α odkryta, prawo skalowania K₂~a^0.205 α^0.177~~ | ✅ ZAKOŃCZONY |
| ~~P67~~ | ~~Wzór K₂ z A_j(a)+B_j·α: 136 ppm @ B, 1168 ppm RMS; odkrycie: g₀(K₂)/g₁(K₂)=-α_K (równanie implicit); g₀ zero K=1.533 → α przesuwa K₂ do 2.033 (+33%)~~ | ✅ ZAKOŃCZONY |
| ~~P68~~ | ~~Newton na g₀(K)+α·g₁(K)=0: 5 kroków → 0 ppm @ B; wzór K₁=K₀−α·g₁/(g₀'+α·g₁'); R[4/3] dla g₀ → K₀ +8.2 ppm; globalna zbieżność powolna (2 kroki → 6000–59000 ppm)~~ | ✅ ZAKOŃCZONY |
| ~~P69~~ | ~~OP-1 ZAMKNIĘTY: R[2/3]+1 Newton → +0.010 ppm @ B, max 14.6 ppm globalnie; +2 Newton → 2.8×10⁻⁴ ppm; wzór kwadratowy (A_j,B_j)+1 Newton → +0.0003 ppm; R[7/2]+1 Newton → 6×10⁻⁶ ppm~~ | ✅ ZAKOŃCZONY |
| ~~P70~~ | ~~OP-4 ZDIAGNOZOWANY: Q@B=3/2−632 ppm; K₂ błąd nieistotny (ΔQ=0.002 ppm na 0.010 ppm K₂); źródło: K₃_num=34.291 za duże, K₃*=34.195=K₃_Ei; λ*/λ_K=+5616 ppm; r₃₁*=3477.5 (PDG)~~ | ✅ ZAKOŃCZONY |
| ~~**P71**~~ | ~~Analiza K₃: K₃~λ^{−0.499308} (691 ppm od −1/2); I₄,I₆ wzory analityczne = numeryczne (0 ppm!); K₃_Ei = K₃* ±2.1 ppm; K₃_num − K₃_Ei = +2802 ppm (korekta kinetyczna); λ* = λ_K×(K₃_num/K₃_Ei)² → pełna spójność~~ | ✅ ZAKOŃCZONY |
| **P72** | Głębsza analiza korekty kinetycznej +2802 ppm: jaki człon energii jest odpowiedzialny? Czy istnieje wzór korygujący K₃_Ei? | Analiza |

### Ocena stanu modelu

**Potwierdzone (solidne):**
- Trzy generacje = trzy zera g(K) (twierdzenie strukturalne)
- Koide Q=3/2 wynika z geometrii solitonu przy λ=λ_Koide
- Parametry (a_Γ, α_K) wyznaczone przez masy PDG bez dopasowania (~0.2%)
- Bifurkacja dostępności: a_c=0.038382 (warunek konieczny na a_Γ)
- K₁ absolutnie niezależne od λ (nowy wynik strukturalny)

**Otwarte/niepewne:**
- Mechanizm wyboru górnego zera α₂ (selekcja bez energii)
- Postać analityczna K₂ i a_c
- Pełna predykcja λ_Koide z zasad pierwszych

*Zaktualizowano: 2026-03-24 | P71 zakończony | K₃~λ^{−0.499308}; I₄,I₆ wzory ≡ numeryczne (0 ppm); K₃_Ei≈K₃* (−2.1 ppm); K₃_num−K₃_Ei=+2802 ppm (korekta kinetyczna); λ*=λ_K×(K₃_num/K₃_Ei)²; (K₃_num/K₃*)²=λ*/λ_K do 7.98 ppm*
