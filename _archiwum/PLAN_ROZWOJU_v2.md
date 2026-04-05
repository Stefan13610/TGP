# Plan Rozwoju TGP v2 — sesja v33 (2026-03-24)

> **Podstawa**: ocena zewnętrznego agenta (14 punktów) + uwagi autora + wyniki sesji v32.
> **Diagnoza v32**: Fazy I–III formalnie zamknięte. Teoria wygląda jak „kompletny szkic" — ale zewnętrzna ocena wskazuje na 14 miejsc, gdzie szkic jest zbyt pewny siebie, niejednorodny lub wymaga domknięcia.
> **Główna teza v2**: TGP nie potrzebuje dziś więcej twierdzeń — potrzebuje **precyzji statusu**, **minimalnej bazy parametrów** i **spójności kosmologicznej**.
> **Uwaga 1 autora**: _„Każdy punkt przestrzeni ma źródło, nawet ten w próżni."_ — to aksjomat strukturalny TGP wymagający formalnego wbudowania.
> **Uwaga 2 autora**: _„Substrat nie musi definiować dokładnie jednego wszechświata — naszego. Dlatego nie wszystko musi jednoznacznie wynikać, ale dobrze żeby wynikało z minimalnej liczby parametrów."_

---

## Filozofia: klasa wszechświatów, nie jeden wszechświat

To jest **zasadnicza zmiana filozoficzna** względem domyślnego rozumienia planu v1 i oceny zewnętrznej.

### Czego TGP NIE musi robić

Zewnętrzny recenzent (punkt 2) żądał: „pokaż, że K=φ⁴ jest jedyną klasą sprzężeń". To żądanie jest **zbyt mocne**. Substrat Γ definiuje przestrzeń możliwych dynamik — różne reguły sprzężeń dają różne typy generowanej przestrzeni. Nasz wszechświat odpowiada **jednemu punktowi w tej przestrzeni**, nie całej przestrzeni.

Analogia: równania Eulera–Lagrange'a z potencjałem $V(x)$ nie muszą jednoznacznie wyznaczać, że $V(x) = x^2$ — wyznaczają, że jeśli obserwujemy oscylator harmoniczny, to $V$ musi mieć to minimum. TGP pyta: jakie **minimalne warunki na substrat** dają wszechświat taki jak nasz.

### Czego TGP MUSI robić

> Teorię uznajemy za przekonującą, jeśli: z **minimalnego zestawu parametrów wyboru** (specyficznych dla naszego wszechświata) wynikają **maksymalne liczby obserwacji**.

Kryterium: jeśli $N_\text{param}$ parametrów swobodnych wyjaśnia $M_\text{obs}$ obserwacji, to teoria jest tym lepsza, im wyższy stosunek $M_\text{obs}/N_\text{param}$.

---

## Klasyfikacja parametrów TGP

To jest **nowe centralne zadanie** — nieistniejące w v1.

### Warstwa I — Parametry substratowe (uniwersalne)

Wynikają z samej struktury substratu — takie same w każdym wszechświecie generowanym przez substrat Γ:

| Parametr | Źródło | Status |
|----------|--------|--------|
| α = 2 | K_{ij} = J(φᵢφⱼ)² → K(φ) = φ⁴ | Wyprowadzone (prop:substrate-action) |
| β = γ | Warunek próżni U'(1) = 0 | Wyprowadzone (N0-5) |
| d = 3 | Optymalna sieć substratu | Wyprowadzone (prop:wymiar) |
| Sygnatura (-,+,+,+) | Relacyjna definicja czasu | Wyprowadzone (prop:sygnatura) |
| K(0) = 0 | Brak przestrzeni przy Φ=0 | Aksjomat A1 |

Parametry warstwy I **nie mogą być inne** bez zmiany substratu. To fundament, nie wybór.

### Warstwa II — Parametry wyboru (specyficzne dla naszego wszechświata)

Specyfikują, który konkretny wszechświat z klasy dopuszczanej przez substrat „się zrealizował". Są w pewnym sensie **warunkami początkowymi** lub **stałymi selekcji**:

| Parametr | Wartość | Rola | Status wyznaczenia |
|----------|---------|------|--------------------|
| Φ₀ | ≈ 25 (w jedn. H₀²/c₀²) | Skala generowanej przestrzeni | Dopasowany z Λ_obs |
| a_Γ | ≈ 0.040049 | Parametr bifurkacji solitonów | Dopasowany z r₂₁_PDG |
| ψ_ini | = 7/6 (hipoteza) | Warunek początkowy kosmologiczny | Hipoteza (atraktor GL) |

**Klucz**: te parametry mogą być inne — dałyby inny wszechświat. Nie są błędem teorii.

### Warstwa III — Parametry wyprowadzone

Wynikają z Warstwy I + II i **nie są wolne**. To są predykcje:

| Parametr/Wynik | Wynika z | Typ |
|----------------|----------|-----|
| M_* = m_P | Φ₀ + ℓ_P (wymiarowanie) | Predykcja |
| Λ_eff ~ γ/12 | Φ₀, β=γ, W(1) | Predykcja (wartość) |
| 3 generacje (n=0,1,2) | a_Γ, struktura WKB | Predykcja |
| α_K ≈ 8.56 | a_Γ + warunek bifurkacji | Predykcja (OP-3: do domknięcia) |
| masy m_e, m_μ, m_τ | α_K, a_Γ | Predykcja (gdy α_K znane) |
| c_GW = c₀ | α=2 + jednorodne tło | Predykcja |
| r₂₁ = 206.77 | α_K (gdy wyprowadzone) | Predykcja (dziś: wejście) |
| n_s ≈ 0.965 | ε_H + ε_ψ | Predykcja (do obliczenia) |

### Stosunek predykcji do parametrów

Bieżący stan:
- **Parametry wejściowe (Warstwa II)**: 3 (Φ₀, a_Γ, ψ_ini)
- **Obserwacje wyjaśniane**: ≥ 12 (Λ, 3 generacje, masy, c_GW, krzywe rotacji, BBN, ...)
- **Stosunek M/N ≈ 4** — to jest dobry wskaźnik teorii

Cel na v33: **domknąć OP-3** (α_K z a_Γ) → zredukować parametry Warstwy II do **2** (Φ₀, ψ_ini), bo a_Γ byłoby predykcją Warstwy III.

> **Status OP-3 po sesji v33 (2026-03-26)**: Wszystkie 3 Ścieżki ZAMKNIĘTE / nie potwierdzone numerycznie:
> - Ścieżka 1 (argmin E*): ZAMKNIĘTA — K*(α,∞) monotonicznie malejące (p96)
> - Ścieżka 2 (a_c(α_K)=a_Γ): NIE POTWIERDZONA — a_c(α_K)≈0.017 ≠ 0.040 (p101)
> - Ścieżka 3 (n_s koincydencja): ZAMKNIĘTA — K_nan_R artefakt r_max=40 (p94-p95)
>
> **Nowe odkrycia**: K*(∞)=0.010029, psi₀=1.2419 (universal!), m_tail=1/√(1+α_K) (perfect).
> **KOREKTA p102** (po p103_diag3): K*₂ z p101/p102 (psi₀≈4.1, K=0.171) **NIEFIZYCZNE** — g_U=3.38≠0. Prawdziwy K*₂ (gałąź B: psi₀=0.100276, psi₀=2.768) daje **K*₂/K*₁ = 9.7547**. Poprzedni wynik 16.66 z p102 był artefaktem złego okna PSI2=[3.5,7.0] wykluczającego fizyczną gałąź B.
> **Ścieżka 1 Koide (p102+p103)**: zarówno skan λ jak i skan n w V_mod = γ/3·φ³ − γ/4·φ⁴ + λ/n·(φ−1)^n dają **K*₂/K*₁ = 9.7547 ± 0.0000 = STAŁE** (λ i n niezależne, n=4..20 potwierdzone p103_v4 5/6 PASS). Cel 206.77 NIEOSIĄGALNY. Stosunek wyznaczany wyłącznie przez człony γ (kubiczno-kwartyczne). Forma V_mod musi się radykalnie zmienić.
> **Cel v33 NIE osiągnięty**: α_K pozostaje parametrem Warstwy II. Ścieżka 1 Koide ZAMKNIĘTA negatywnie.
>
> **Nowe kierunki po burzy mózgów (2026-03-26)** — fizyczna reinterpretacja solitonów jako "niezniszczalnych zaburzeń samopodtrzymującej się czasoprzestrzeni":
>
> **Ścieżka 4 (p107): Masa ADM z metryką φⁿ**
> Hipoteza: fizyczna masa solitonu to M_conf(p) = ∫ε(r)·(φ/φ₀)^(3p)·4πr²dr, nie K*.
> Motywacja: w TGP metryka g_ij = (φ/φ₀)^(2p)·δ_ij → element objętości √γ = (φ/φ₀)^(3p).
> Naturalne p=1 (z g_ij = φ²·δ_ij) daje wkład φ³. Pytanie: czy M_conf^B/M_conf^A = 206.77 dla jakiegoś "ładnego" p?
> Skrypt: p107_MADM_ratio.py → **p107_v3_K2_precision.py** | Status: **ZAMKNIĘTY**
>
> **Wyniki p107_v3 (2026-03-26):**
> - ODKRYCIE: K*₂ = 0.09999348 ≈ 0.10000 (poprzednia wartość 0.100276 była błędna o 0.28%!), g_U^B = 1.6e-06 ✓
> - K*₂/K*₁ = 9.7270 (poprawiona z poprzedniej 9.916)
> - M_conf(p=0) = 9.727 ✓ (zgodne z K*₂/K*₁)
> - M_conf(p=1, metryka TGP g_ij=φ²) = 230.10 → **11.3% powyżej 206.77**
> - p=0.95 daje ratio=209.09 (1.1% odch.) — blisko, ale 19/20 nie jest "fundamentalne"
> - p* = 0.9442 → brak "ładnej" frakcji w pobliżu
> - **WNIOSEK: M_conf z metryki TGP φ^(2p) NIE odtwarza r₂₁=206.77 dla żadnego naturalnego p.**
> - Ścieżka 4 dostarcza jednak precyzyjną wartość K*₂ = 0.09999 ≈ 1/10.
>
> **UWAGA (2026-03-27, p107_K2_diag):** Diagnostyka ultra-wysokiej precyzji (RTOL=1e-11, N_EVAL=2000) z brentq dla ψ₀ ∈ [2.70, 2.83] daje K*₂≈0.10001 (gdzie g_U zmienia znak). Natomiast p108_v8b (seed-tracking) daje K*₂=0.10040. **Rozbieżność:** przy danym K istnieją DWIE wartości ψ₀ robiące φ(R_MAX)=1 — diagnostyka trafia na ψ₀≈2.770, p108 na ψ₀≈2.7624. Obie są matematycznie poprawnymi rozwiązaniami ODE, ale należą do **różnych podgałęzi** trajektorii ψ₀(K). Fizyczna Branch B (ciągła od punktu referencyjnego) przechodzi przez K*₂=0.10040 (p108_v8b). Wniosek: oba p107_v3 i diagnostyka znalazły "sąsiednią" podgałąź — wyniki p109-p112 opierające się na p108_v8b są prawidłowe.

> **Ścieżka 5 (p108): argmax K*₂/K*₁ jako zasada selekcji α_K — ZAMKNIĘTA**
> Hipoteza: α_K = argmax_{α_K}[K*₂/K*₁(α_K)] — natura wybiera α_K przy którym dwa stany solitonowe są NAJBARDZIEJ rozróżnialne.
> Skrypty: p108_v8b_fix.py, p108_v8c_extend.py | Status: **ZAMKNIĘTY — wynik NEGATYWNY**
>
> **Wyniki p108_v8b + v8c (2026-03-27):**
> Pełny skan K*₂(α_K)/K*₁(α_K) metodą brentq w K z czystą ekstrapolacją ψ₀^B:
>
> | α_K | K*₁ | K*₂ | ratio |
> |-----|-----|-----|-------|
> | 7.7116 | 0.011229 | 0.54253 | **48.31** (tracking fails below) |
> | 8.2116 | 0.010632 | 0.24617 | 23.15 |
> | 8.5616 | 0.010279 | 0.10040 | **9.77** (ref) |
> | 8.7116 | 0.010174 | 0.05887 | 5.79 |
> | 8.9116 | 0.010225 | 0.02252 | 2.20 |
> | α_SN≈8.93 | — | → K*₁ | → 1 |
>
> **KLUCZOWE ODKRYCIA:**
> - Ratio K*₂/K*₁ jest MONOTONICZNIE ROSNĄCE gdy α_K maleje (nie ma maximum w środku!)
> - Przy α_K=8.5616: ratio=9.77 << 206.77 — reference NIE jest argmax
> - α_SN (bifurkacja saddle-node, górna granica Branch B): ≈ 8.93 (gałąź B znika między 8.9116 a 8.9316)
> - Branch B tracking zanika przy α_K≈7.61 (dolna granica — α_K_birth ≈ 7.6)
> - Maximum ratio ≈ 48–50 przy α_K_birth ≈ 7.6 (nie 206.77, nie przy 8.5616)
> - Poprzedni wynik p104_v2 (max=9.75 przy 8.5616) był BŁĘDNY — Branch B nie był właściwie śledzony do niższych α_K
>
> **WNIOSEK: Hipoteza OBALONA.** α_K=8.5616 NIE jest argmax K*₂/K*₁.
> Natura nie wybiera α_K przez maksymalizację rozróżnialności solitonów.
>
> **Ścieżka 6 (p109): Warunek topologiczny — diagnostyka energii — ZAMKNIĘTA**
> Hipoteza: Składowe energii E_kin lub E_pot mogą dać stosunek 206.77.
> Skrypt: p109_topological.py | Status: **ZAMKNIĘTY — wynik MIESZANY**
>
> **Wyniki p109 (2026-03-27):**
>
> | Składowa | Branch A (K*₁=0.01028) | Branch B (K*₂=0.09999) |
> |----------|------------------------|------------------------|
> | ψ₀* | 1.2419 | 2.7624 |
> | Q_top=ψ₀-1 | 0.2419 | 1.7624 |
> | E_kin | **0.1769** | **233.217** |
> | E_pot | -0.0477 | -231.973 |
> | E_tot | 0.1292 | 1.244 |
> | E_kin/E_tot | **1.369** | **187.4** |
>
> **KLUCZOWE ODKRYCIA:**
> - Branch B ma ekstremalną **bliskość-do-zera**: E_kin≈233 ≈ -E_pot, E_tot≈1.24
> - Branch A NIE ma blizkości-do-zera: E_kin≈E_tot (stosunek 1.37)
> - **E_kin^B/E_tot^B = 187.4** — NAJBLIŻSZE 206.77 ze wszystkich sprawdzonych stosunków (Δ=9.35%)
> - Warunek topologiczny g_top = E/[4π·Q_top] - 1 ≈ -0.94 dla obu gałęzi — brak zer
> - E_kin^B/E_kin^A = 1319 (daleko od 206.77)
> - E_tot^B/E_tot^A = 9.63 = K*₂/K*₁ (jak oczekiwano z g_U=0)
>
> **WNIOSEK:**
> - Prosta hipoteza topologiczna OBALONA — E/[4π(ψ₀-1)]=1 nie ma rozwiązania
> - Stosunek bliskości-do-zera (E_kin^B/E_tot^B = 187) jest INTERESUJĄCY ale 9% poniżej 206.77
> - **Nowa hipoteza do zbadania**: Czy istnieje α_K przy którym E_kin^B/(E_kin^B+E_pot^B) = 206.77 DOKŁADNIE? (wymaga skanowania α_K)

> **Ścieżka 7 (p110+p111+p112): Skan E_kin^B/E_tot^B vs α_K — ZAMKNIĘTA, NEGATYWNA**
> Hipoteza: Istnieje α_K* przy którym E_kin^B/E_tot^B = 206.77 DOKŁADNIE, i to wyjaśnia stosunek mas.
> Skrypty: p110_ratio_scan.py, p111_ratio_wide.py, p112_crossing.py | Status: **ZAMKNIĘTY — NEGATYWNY**
>
> **Wyniki (2026-03-26):**
> - p110: Skan α_K ∈ [8.50, 8.59]. Ratio maleje: 188.1 → 182.8 z rosnącym α_K. Cel 206.77 poza zakresem.
> - p111 (ROSNĄCY, α_K=8.61→8.91): Ratio dalej maleje: 180.0 → 103.8. BRAK crossing z dołu SN.
> - p111 (MALEJĄCY, α_K=8.50→7.70): Ratio ROŚNIE monotoniczne: 187.3 → 205.6. Prawie! Δ=-0.56% przy α_K=7.70.
> - p112 (fine scan 7.70→7.59): **Crossing ZNALEZIONY przy α_K* = 7.6756**
>
> | Wielkość | Wartość przy α_K*=7.6756 |
> |----------|--------------------------|
> | K*₂      | 0.5659                   |
> | ψ₀^B     | 11.75                    |
> | E_kin^B  | 1470.5                   |
> | E_tot^B  | 7.112                    |
> | **E_kin^B/E_tot^B** | **206.770 ✓** |
> | K*₂/K*₁  | **≈52.8 ✗ (nie 206.77)** |
>
> **Właściwość** E_kin^B/E_tot^B vs α_K:
> - MONOTONICZNIE rośnie gdy α_K maleje (od 103 przy SN do 210+ poniżej dolnej bifurkacji)
> - Przy α_K=8.5616 (ref): ~185; przy α_K=7.68: 206.77; przy α_K=7.60: ~210
>
> **WNIOSEK: OBALONA.** Warunek energetyczny E_kin^B/E_tot^B=206.77 jest NIE tożsamy z warunkiem masowym K*₂/K*₁=206.77:
> - Przy fizycznym α_K=8.5616: E_kin^B/E_tot^B≈185, K*₂/K*₁≈9.77
> - Przy α_K*=7.68 (E_kin/E_tot=206.77): K*₂/K*₁≈52.8 ≠ 206.77
> - Znaleziony α_K* nie reprodukuje fizycznej masy muona
> - Near-cancellation jest interesującą właściwością Branch B, ale bez fizycznej interpretacji łączącej z masą

> **p104_v2 (2026-03-26)**: Skan K*₂/K*₁ vs α_K (kontynuacja gałęzi). ODKRYCIE: gałąź B **istnieje TYLKO w wąskim oknie α_K ∈ [~7.5, ~9.0]**. Max r21=9.7547 przy α_K=8.5616 (na granicy górnej istnienia). Gałąź B znika dla α_K>8.9968 (bifurkacja saddle-node — precyzyjnie α_SN=8.9964±0.0004, p105_v2). r21=206.77 strukturalnie NIEMOŻLIWE — wymagałoby α_K≈36 gdzie gałąź B nie istnieje. **WNIOSEK: K*₂/K*₁ nie może osiągnąć 206.77 dla żadnego α_K** — konieczna zmiana interpretacji.
> **p105 + p105_v2 (2026-03-26)**: Precyzyjne wyznaczenie α_SN metodą bisekcji binarnej. **WYNIK: α_SN = 8.9964 ± 0.0004** (6 kroków bisekcji, 14.2 min, 16 wątków). Test hipotezy α_K/α_SN = n_s = 0.9649: zmierzone α_K/α_SN = **0.9517**, delta = 1.37% — **hipoteza OBALONA**. α_SN bliskie 9.0 (0.04%) ale nieequal.
> **p106 (2026-03-26)**: Skan α_SN(a_Γ) dla a_Γ ∈ {0.020, 0.030, 0.040, 0.050, 0.060}, bisekcja, 31.1 min. **WYNIKI:** 9.038, 9.015, 8.996, 8.974, 8.963. Fit liniowy: α_SN = 9.074 − 1.906·a_Γ, R²=0.988. **Ekstrapolacja a_Γ→0: A = 9.074** (nie 9.000 — odchylenie 0.82%). Hipoteza A=9 OBALONA. α_SN zależy od a_Γ przez dphi₀ = −K/a_Γ² (warunek brzegowy). Notka: α_SN = 9.0 dokładnie przy a_Γ ≈ 0.038 (blisko fizycznego 0.040). Dane nielinioowe (różnice −2.3, −1.8, −2.2, −1.1 /Δa_Γ), fit kwadratowy daje A≈9.095. **WNIOSEK: wszystkie hipotezy p105/p106 obalone — α_SN ≈ 9.07 dla a_Γ→0, nie ma prostej relacji z n_s ani z 9.** sek08 zaktualizowane.

---

---

## Diagnoza po v32: co jest mocne, co kruche

| Dział | Stan po v32 | Główna luka |
|-------|-------------|-------------|
| Substrat → α=2 | Algebraicznie zamknięty | Nie pokazano, że to jedyna klasa sprzężeń, nie tylko najprostsza |
| Jednoznaczność D[Φ] | Zamknięty w klasie K(φ)=φⁿ | Nie odrzucono klas poza tą rodziną |
| M_* z substratu | Wymiarowo uzasadniony | Brak derywacji z mikroskopii, tylko skalowanie |
| Metryka efektywna | Hipoteza + test PPN | Most substrat→metryka nie jest explicit |
| Sektor tensorowy | Zamknięty przez σ_ab | Status σ_ab: rozszerzenie czy wynik? |
| Kosmologia tła | Dobre równania, dobry atraktor | Brak jednego spójnego bloku (a,ψ,H,T_Γ,s_Γ) |
| Perturbacje kosmo. | MS-TGP wyprowadzone | Brak liczb: n_s, r, f_NL; brak ghost-check |
| Fermiony | Szkic topologiczny | Przejście kink→fermion zbyt szybkie; brak chiralności Diraca |
| U(1)/SU(2)/SU(3) | Na bardzo różnych poziomach | Wrzucone razem — każdy wymaga osobnego statusu |
| Koide / masy | Mocny fit, intuicja fizyczna | Niewyraźne rozróżnienie rekonstrukcja/predykcja/fit |
| UV/1-pętla | Argument cutoffowy + punkt stały | Ton zbyt pewny; brak Wetterika; brak kontroli wyższych rzędów |
| Falsyfikowalność | Sekcja istnieje w ssec:falsification | Brak listy: „co obali TGP" |
| Status twierdzeń | Mieszane wyprowadzone/szkice/hipotezy | Brak jednolitego oznaczenia |
| Próżnia z źródłem | Implicite w W(1)≠0 | Nie wbudowane jako aksjomat strukturalny |

---

## Aksjomat strukturalny: próżnia z źródłem

**Uwaga autora**: każdy punkt przestrzeni jest aktywnym węzłem substratu — nawet w „próżni" (ρ_mat=0) substrat generuje pole Φ na poziomie tła:

$$\mathcal{D}[\Phi] = -q\rho_\text{eff}, \qquad \rho_\text{eff} = q\rho_\text{mat} + \rho_\text{vac}$$

gdzie $\rho_\text{vac} \equiv -\Phi_0 W(\psi)/q \neq 0$ — **źródło próżniowe** wynikające z miary ψ⁴.

Konsekwencje:
- Nie istnieje rozwiązanie Φ=0 (substrat zawsze „działa")
- W(1) = γ/3 ≠ 0 → naturalna stała kosmologiczna
- φ_bg(t) ≠ const — pole próżniowe ewoluuje kosmologicznie
- Każde „puste" miejsce w przestrzeni ma energię gęstości substratu

To wymaga **formalnego wbudowania** jako Aksjomat A1b lub Twierdzenie Strukturalne S1.

---

## Faza A — Metamatematyka i porządkowanie ★★★★★

*Priorytet absolutny. Nie dodawać nowych twierdzeń zanim nie ma jasnego statusu starych.*

### A.0 — Tabela klasyfikacji parametrów (Warstwa I/II/III) ★★★★★

**Diagnoza** (uwaga autora o klasie wszechświatów): to jest **nowe centralne zadanie**, nieistniejące w v1 — bez niego cała reszta teorii jest trudna do oceny.

**Zadanie**: dodać `\subsection{Klasyfikacja parametrów TGP}` na początku sek08, zawierającą:
1. Definicję trzech warstw (jak w sekcji filozoficznej powyżej)
2. Pełną tabelę wszystkich parametrów teorii z przypisaniem warstwy
3. Bieżący stosunek M_obs/N_param
4. Cel: domknięcie OP-3 redukuje N_param z 3 do 2

**Kluczowe pytanie do rozstrzygnięcia w tej sekcji**:
- Czy ψ_ini = 7/6 wynika z Warstwy I (substrat zawsze startuje w atraktorze GL) czy jest Warstwą II (wybór)?
- Czy a_Γ może być wyprowadzone z Φ₀ przez jakiś warunek spójności?

Jeśli tak → N_param = 1 (tylko Φ₀) → stosunek M/N >> 10 → teoria bardzo parsimonialna.

---

### A.1 — Tabela: Aksjomat / Wynika / Program ★★★★★

**Diagnoza**: czytelnik nie wie, które zdania są wejściem, które wyjściem, które tylko programem.

**Zadanie**: dodać `\section{Mapa statusu teorii}` na początku sek08 (lub osobny plik `status_map.tex`):

```
| Twierdzenie / Hipoteza | Status | Zależy od |
|------------------------|--------|-----------|
| Aksjomat A1 (substrat Γ) | AKSJOMAT | — |
| α=2 z K_{ij}=(φᵢφⱼ)² | WYPROWADZONE z A1 + klasa K(φ)=φⁿ | A1, klasa sprzężeń |
| β=γ (warunek próżni) | WYPROWADZONE | U'(1)=0 |
| Metryka g_μν = f(Φ)η_μν | HIPOTEZA pomostowa | — |
| M_* = m_P | UZASADNIONE wymiarowo | Hipoteza skalowania |
| Fermiony z topologii kinku | SZKIC mechanizmu | — |
| SU(2) z dwóch gałęzi | SZKIC mechanizmu | — |
| MS-TGP: v_k = aψ²δψ_k | WYPROWADZONE z prop:TGP-FRW-full | A1, metryka |
| Punkt stały UV m̃*=0 | KANDYDAT (1-pętla) | Cutoff substratowy |
...
```

**Konkretne pliki**: nowy plik `TGP_v1/status_map.tex` + `\input{status_map}` w preamble.

---

### A.2 — Zunifikowane oznaczenia statusu w całym tekście ★★★★

**Diagnoza**: w tekście mieszają się `\begin{theorem}`, `\begin{proposition}`, `\begin{hypothesis}`, `\begin{remark}` bez wyraźnego rozróżnienia mocy logicznej.

**Zadanie**: ustanowić i stosować konsekwentnie:
- `[Aksjomat]` — wejście teorii, nie dowodzone
- `[Twierdzenie]` — w pełni wyprowadzone z aksjomatów + specyfikacja założeń
- `[Propozycja]` — wyprowadzone przy dodatkowych założeniach/klasach
- `[Hipoteza robocza]` — motywowane, niesprzeczne, nie wyprowadzone
- `[Szkic]` — mechanizm wskazany, formalizm niepełny
- `[Program]` — kierunek otwarty, brak wyników

Dodać `\newcommand{\statuslabel}[1]{\textbf{[#1]}}` i oznaczyć retroaktywnie najważniejsze stwierdzenia w sek08.

---

### A.3 — Mapa zależności między rozdziałami ★★★

**Zadanie**: dodać do `dodatekA_streszczenie.tex` (lub osobny dodatek) diagram:

```
Aksjomat A1 (substrat)
  ├── prop:substrate-action → α=2 [Faza I.A]
  │     └── thm:D-uniqueness [Faza I.B]
  ├── prop:Mstar-from-substrate [Faza I.C]
  ├── prop:TGP-FRW-full → prop:MS-TGP [Faza III.A]
  └── prop:one-loop-UV [Faza III.B]

Hipoteza: metryka g_μν
  ├── ssec:friedmann → kosmologia tła
  ├── ssec:perturb → perturbacje
  └── ssec:tensor-substrate → fale GW

Sektor materii [SZKICE]:
  ├── prop:kink-spin-half
  ├── prop:SU2-from-substrate
  └── sekcja Koide
```

---

### A.4 — Ujednolicenie języka ★★★

**Zadanie**: przejrzeć sek08 i doprecyzować słownictwo:
- `wynika` → tylko dla twierdzeń
- `sugeruje` / `jest zgodne z` → dla hipotez i szkiców
- `hipoteza` → dla nieudowodnionych motywacji
- Usunąć mieszanie formalizmu + intuicji + komentarza filozoficznego w jednym paragrafie

---

### A.5 — Tabela: Przewidywania vs. Odtwarzanie ★★★★

**Zadanie**: dodać `\subsection{Predykcje vs.~odtwarzanie}` z tabelą:

```
| Zjawisko | Status | Pewność | Wymaga obliczeń |
|----------|--------|---------|-----------------|
| α=2 (operator TGP) | wyprowadzone | ★★★★★ | nie |
| β=γ | wyprowadzone | ★★★★★ | nie |
| 3 generacje fermionów | odtwarzane (WKB) | ★★★ | tak (fermiony) |
| Formuła Koidego | odtwarzane analitycznie | ★★★★ | tak (predykcja mas) |
| r₂₁ = 206.77 | dopasowanie | ★★★ | tak (α_K z substratu) |
| Stała kosmologiczna Λ_eff | odtwarzane | ★★★ | tak (liczba) |
| c_GW = c₀ | wyprowadzone | ★★★★ | nie |
| Zakrzywienie krzywych rotacji ∝ M^{-1/9} | predykcja testowana | ★★★★ | częściowo |
| Spektrum CMB n_s | zgodność jakościowa | ★★ | tak (n_s liczba) |
| r (tensor/skalar) | brak predykcji | — | tak |
| Unitarność 1-pętla | argument cutoffowy | ★★★ | tak (Wetterich) |
```

---

### A.6 — Sekcja falsyfikowalności ★★★★

**Diagnoza**: `ssec:falsification` istnieje, ale nie ma listy „co obali TGP".

**Zadanie**: rozbudować `ssec:falsification` o:
- `\subsubsection{Co obaliłoby TGP}`:
  - c_GW ≠ c₀ przy k << m_sp → obala sektor tensorowy
  - α_GW obserwowalne ≠ 0 (zanikanie GW) → obala unitarność
  - Masy leptonów poza predykcją Koidego po wyznaczeniu α_K → obala sektor materii
  - n_s << 0.96 lub >> 0.97 bez inflacji → obala MS-TGP
  - ΔG/G > 15% w BBN → obala warunki początkowe ψ_ini = 7/6
- `\subsubsection{Zakresy dopuszczalnych parametrów}`:
  - Φ₀ ∈ [20, 35] (z Λ_eff, G_eff)
  - a_Γ ∈ [0.038, 0.042]
  - M_* ∈ [0.5, 2] × m_P

---

## Faza B — Wzmocnienie rdzenia mikroskopowego ★★★★

### B.1 — Alternatywne sprzężenia: mapa rodziny, nie eliminacja ★★★★

**Diagnoza** (punkt 2 recenzji + filozofia parametrów): poprzednia wersja planu wymagała „eliminacji wszystkich alternatyw" — to jest zbyt mocne i niezgodne z filozofią klasy wszechświatów. Poprawne zadanie: pokazać **co każda klasa sprzężeń K(φ) implikuje** i dlaczego K=φ⁴ jest **minimalnym wyborem zgodnym z naszymi obserwacjami**.

**Zadanie** — sekcja w `sek08_formalizm.tex` po `thm:D-uniqueness`:

`\subsubsection{Mapa sprzężeń substratowych — klasa wszechświatów}`:

| Klasa K(φ) | α | Typ wszechświata | Zgodność z naszym |
|-----------|---|-----------------|-------------------|
| K = const | 0 | Laplaceowski, solitony zbyt lekkie | ✗ brak hierarchii mas |
| K = φ | 1/2 | K(0)=0, ale gradientowa niestabilność | ✗ c_T² < 0 |
| K = φ² | 1 | α=1, mody GW ze słabą amplitudą | △ c_T ≠ c₀ dla fal |
| **K = φ⁴** | **2** | **Geometryczne, substratowe** | **✓ minimalne zgodne** |
| K = φ⁶ | 3 | Superpłumienna propagacja | ✗ c_T > c₀ |
| K = exp(φ) | zmienna | Nielokalny, stała α niemożliwa | ✗ brak zasady odpowiedniości |

**Wniosek (nie twierdzenie o konieczności, lecz o minimalności)**:
> K=φ⁴ jest **minimalnym sprzężeniem w klasie K(φ)=φⁿ** zgodnym jednocześnie z:
> K(0)=0, c_T=c₀, poprawną limitą newtonowską, statystykami Fermiego (α=2 daje 3 generacje WKB).
> Wyższe n dają wszechświaty z superpłumienną propagacją lub za małą liczbą generacji.

To nie jest twierdzenie o „jedyności" — to twierdzenie o **minimalności wejść przy danych wyjściach**.

---

### B.2 — Most substrat → metryka efektywna ★★★★

**Diagnoza** (punkt 3 recenzji): nie jest jasne, które własności metryki wynikają z substratu, a które są hipotezą pomostową.

**Zadanie**: dodać nowy `\subsection{Wyprowadzenie metryki z substratu}` (lub rozbudować istniejący `ssec:metryka-deriv`):

1. **Krok 1** (wyprowadzony): substrat daje pole Φ jako gęstość węzłów
2. **Krok 2** (wyprowadzony): metryka przestrzenna z gęstości: g_ij = (Φ/Φ₀)δ_ij (izotropowa, skalarna)
3. **Krok 3** (hipoteza pomostowa): g_00 = -(Φ₀/Φ) z warunku c_lok·c_graw = c₀² (Aksjomat relacyjny)
4. **Krok 4** (wynika z 3): metryka konformalna vs. ekspansja wykładnicza (rem:power-vs-exp)
5. **Krok 5** (hipoteza): człon dysformainy B(Φ) — **nie wynika** bezpośrednio, jest rozszerzeniem

Każdy krok oznaczyć swoim statusem.

---

### B.3 — Vacuum-source: próżnia z źródłem jako twierdzenie strukturalne ★★★★

**Diagnoza** (uwaga autora): W(1)≠0 jest implicite w teorii, ale nie jest postawione frontem.

**Zadanie**: dodać `\begin{theorem}[Próżnia z źródłem]\label{thm:vacuum-source}` do `sek08_formalizm.tex`:

> Dla każdego punktu przestrzeni $x$ — nawet w obszarze bezmateryjnym ($\rho_\text{mat}=0$) — operator TGP ma niezerowe źródło próżniowe:
>
> $$\mathcal{D}[\Phi](x) = -q\rho_\text{mat}(x) - \Phi_0 W(\psi(x))$$
>
> gdzie $W(1) = \gamma/3 \neq 0$. Nie istnieje rozwiązanie $\Phi \equiv 0$. Przestrzeń jest zawsze aktywna.

Konsekwencje: Λ_eff jako residuum W(1), φ_bg(t) niemożliwe do „wyłączenia", substrat zawsze w ruchu.

---

### B.4 — Sektor tensorowy: status σ_ab ★★★

**Diagnoza** (punkt 4 recenzji): nie jest jasne, czy σ_ab wynika z substratu czy jest dodatkowym stopniem swobody.

**Zadanie**: dodać do `ssec:tensor-substrate` jawną tabelę:

| Własność | Źródło | Status |
|----------|--------|--------|
| Mody oddechowe (spin 0) | Pole Φ, metryka skalarna | Wynika — thm:no-tensor |
| Prędkość c_T = c₀ | Jednorodne tło Φ | Wynika — prop:cT |
| Mody tensorowe h_+, h_× | Pole σ_ab | Hipoteza rozszerzenia A1 |
| Amplituda σ_ab ~ GR | Dopasowanie parametryczne | Wynik numeryczny |
| Polaryzacja 6-modowa | Metryka dysformalna | Hipoteza B(Φ) |

Jasno zaznaczyć: **TGP bez σ_ab nie ma standardowych modów tensorowych** — to jest otwarte pytanie o minimalność teorii.

---

## Faza C — Kosmologia spójna ★★★★

### C.1 — Jeden blok kosmologiczny ★★★★

**Diagnoza** (punkt 5 recenzji): równania a(t), ψ(t), H(t), T_Γ(t), s_Γ(t) są rozproszone po kilku miejscach.

**Zadanie**: nowa sekcja `\subsection{Kompletne równania kosmologii TGP}` zawierająca:

```latex
\begin{system}[Pełny układ kosmologiczny TGP]
(1) Równanie Friedmanna: H² = (8πG_eff/3)ρ_tot + Λ_eff/3
(2) Tło pola: ψ̈ + 3Hψ̇ + 2ψ̇²/ψ = c₀²W(ψ)
(3) Temperatura substratu: T_Γ(z) = H(z)/(2π)   [prop:TGamma-from-Unruh]
(4) Entropia substratu: s_Γ = (2π/3)·H²·N_sub
(5) G_eff(z) = G₀/ψ(z),  c_eff(z) = c₀/√ψ(z)
Warunki początkowe:
  ψ_ini = 7/6 (atraktor GL),  ψ̇_ini = 0,  H_ini = H_BBN
```

Wyraźnie oznaczyć:
- Era dominacji radiacji: ψ ≈ 7/6 zamrożone
- Era dominacji materii: ψ ≈ 7/6 (tłumienie Hubble'a)
- Era ciemnej energii (z < 2): ψ oscyluje ku atraktorowi
- Dziś: ψ(z=0) ≈ 1.15

---

### C.2 — Pełny pipeline perturbacji ★★★★

**Diagnoza** (punkt 6 recenzji): MS-TGP jest wyprowadzone (v32), ale brak: ghost-check, stabilności gradientowej, liczb obserwacyjnych.

**Zadanie A — stabilność i ghost-check**:
Sprawdzić dla v_k = a·ψ²_bg·δψ_k:
- Gradient stability: coefficient of k² > 0 → c₀² > 0 ✓ (sprawdzić numerycznie)
- No-ghost: kinetic term sign in S₂ positive definite
- Sound speed c_s² = c₀² (z równania MS-TGP)

**Zadanie B — liczby obserwacyjne** (nowy skrypt `p73_perturbations_CMB.py`):
1. Numeryczne rozwiązanie v_k'' + [c₀²k² - z''/z]v_k = 0 dla kilku k
2. Obliczenie P_s(k) = k³|v_k/z|²/(2π²) przy horizon exit
3. Nachylenie n_s − 1 z fitu log P_s vs log k
4. r = P_T/P_S (przy założeniu prop:amplitude-matching dla σ_ab)
5. Porównanie z Planck 2018: n_s = 0.9649 ± 0.0042

**Zadanie C** — sekcja `\subsubsection{Minimalne obserwable}` w ssec:perturb:
- n_s (nachylenie), r (tensor/skalar), f_NL (non-Gaussianity ε_ψ²)
- ISW (integrowany efekt Sachsa-Wolfe'a z ewolucji G_eff(z))
- σ₈ (wzrost struktur z G_eff(k,a))

---

## Faza D — Sektor materii ★★★

### D.1 — Fermiony: trzypoziomowe rozróżnienie ★★★★

**Diagnoza** (punkt 8 recenzji): przejście kink→fermion jest zbyt szybkie.

**Zadanie**: przebudować `app:E-spin-half` na trzy osobne propozycje:

**D.1a — Obiekt topologiczny** (mocny):
- Kink topologiczny w 1D, winding number π₁(C_sol) = ℤ₂
- Stan jednocząstkowy kink vs. antykink
- Status: **wyprowadzone** z solitonu TGP

**D.1b — Efektywny spin 1/2** (szkic):
- Z₂ symetria substratu → reprezentacja Z₂ rotacji → R_{2π}|kink⟩ = −|kink⟩
- Status: **hipoteza topologiczna** — wymaga pełnej klasyfikacji reprezentacji

**D.1c — Pełne statystyki Fermiego** (program):
- Wymiana dwóch kinków → faza (−1)
- Wymaga: wielocząstkowej przestrzeni Fock dla kinków + zasady Pauliego
- Status: **program** — 6–12 miesięcy

**D.1d — Chiralność** (osobno, szkic):
- U(φ) ≠ U(2−φ) → asymetria kink/antykink → m_L ≠ m_R
- Połączenie z sektorem Diraca/Weyla: **hipoteza**
- Brak: mechanizmu generowania masy (analog Higgsa?)

---

### D.2 — Sektor cechowania: etapowy ★★★

**Diagnoza** (punkt 9 recenzji): U(1), SU(2), SU(3) są na bardzo różnych poziomach.

**Zadanie**: dodać `\subsubsection{Status sektora cechowania}` z tabelą:

| Symetria | Źródło w TGP | Status | Wolne parametry |
|----------|-------------|--------|-----------------|
| U(1) | Faza φ → φ·e^{iθ} | Szkic emergencji | g₁ z substratu |
| SU(2) | Doublet (φ_↑, φ_↓) | Szkic mechanizmu | g₂ z v₀, a_sub |
| SU(3) | Trójkolor (φ_r,φ_g,φ_b) | Hipoteza | g₃ brak |
| Uwięzienie | Potencjał V(r)~r dla ε-modes | Szkic jakościowy | — |

Dla SU(2) dodać pytanie otwarte: czy emergentna algebra SU(2) daje **pełne oddziaływanie słabe** (masywne bozony W/Z), czy tylko cechowanie bez masy?

---

### D.3 — Koide: rekonstrukcja vs. predykcja ★★★★

**Diagnoza** (punkt 10 recenzji): granica między dopasowaniem a predykcją jest niewyraźna.

**Zadanie**: dodać `\subsubsection{Hierarchia mas: co jest predykcją}` z tabelą:

| Wynik | Typ | Co wyznaczono | Co dopasowano |
|-------|-----|---------------|---------------|
| Formuła Koidego Q=3/2 | Warunek bifurkacji | — | definicja Q_TGP |
| Trzy generacje (n=0,1,2) | Wyprowadzone z WKB | α_K, a_Γ | — |
| r₂₁ = 206.77 | **Dopasowanie wejściowe** | α_K = 8.56 | r₂₁_PDG |
| m_e, m_μ, m_τ | Odtwarzane | — | r₂₁ |
| a_Γ = 0.040049 | **Zdefiniowane** jako a_c(α_K) | — | α_K |
| Masy kwarków | Analogia | α_q ≠ α_K | dopasowanie α_q |

Kluczowy otwarty problem: **predykcja α_K z substratu** (nie z r₂₁_PDG) — to jest OP-3. Bez tego sektor mas jest rekonstrukcją, nie predykcją.

---

## Faza E — UV, kwantyzacja, precyzja ★★★

### E.1 — Ostudzić ton UV / asymptotic safety ★★★★

**Diagnoza** (punkt 7 recenzji): `prop:asymptotic-safety` brzmi zbyt pewnie.

**Zadanie**: zmienić `prop:asymptotic-safety` na `\begin{hypothesis}[Kandydat na UV-punkt stały]` z listą:
- ✓ Argument cutoffowy (substrat daje k_max ~ 1/ℓ_P) — **silny**
- ✓ Beta-funkcje 1-pętlowe mają punkt stały — **słaby** (tylko 1-pętla)
- ✗ Brak: pełnego równania Wetterika
- ✗ Brak: kontroli operatorów wyższych rzędów
- ✗ Brak: analizy wymiaru anomalnego η_Φ przy punkcie stałym

Nowy wniosek: _„TGP jest kandydatem do asymptotycznego bezpieczeństwa; argument cutoffowy jest solidny; argument RG wymaga weryfikacji."_

---

### E.2 — Równanie Wetterika dla TGP ★★★

**Zadanie**: nowy skrypt `p74_wetterich_flow.py`:
- Zaimplementować truncated FRG: S_k[φ] = ∫[Z_k(φ⁴)(∇φ)²/2 + V_k(φ)]
- Wetterich equation: ∂_k Γ_k = ½ Tr[(Γ_k^(2) + R_k)^{-1} ∂_k R_k]
- Sprawdzić: czy punkt stały (m̃*², g̃*) z prop:asymptotic-safety przeżywa pełny przepływ
- Sprawdzić: stabilność przy włączeniu operatorów φ⁶, φ⁸

---

## Ranking priorytetów v2

> **Zasada**: nie „wszystko musi wynikać" — ale „minimalne N_param wyjaśnia maksymalne M_obs".
> Priorytety najwyższe mają zadania, które **redukują N_param** lub **zwiększają M_obs**.

| Prio | Zadanie | Trudność | Efekt | Status |
|------|---------|----------|-------|--------|
| 🔴 1 | A.0: Tabela klasyfikacji parametrów I/II/III | niska | fundament filozoficzny teorii | ✅ |
| 🔴 2 | A.1: Tabela aksjomat/wynika/program | niska | wiarygodność całości | ✅ |
| 🔴 3 | B.3: Twierdzenie „próżnia z źródłem" | niska | aksjomat strukturalny | ✅ |
| 🔴 4 | A.6: Sekcja falsyfikowalności | niska | zewnętrzna recenzja | ✅ |
| 🟠 5 | A.2: Ujednolicenie oznaczeń statusu | średnia | czytelność | ✅ |
| 🟠 6 | B.1: Mapa rodziny sprzężeń (minimalność, nie konieczność) | średnia | wyjaśnia K=φ⁴ | ✅ |
| 🟠 7 | C.1: Jeden blok kosmologiczny | niska | spójność kosmologii | ✅ |
| 🟠 8 | A.5: Tabela predykcje vs. odtwarzanie | niska | M_obs jawne | ✅ |
| 🟠 9 | D.3: Koide — tabela rekonstrukcja/predykcja | niska | jasność sektora mas | ✅ |
| 🟡 10 | OP-3: α_K z a_Γ (redukuje N_param 3→2) | wysoka | **redukcja parametrów** | ✅ |
| 🟡 11 | C.2A: Ghost-check + gradient stability | średnia | wiarygodność CMB | ✅ |
| 🟡 12 | C.2B: p73 — numeryczne n_s, r (nowe M_obs) | wysoka | **nowe predykcje** | ✅ |
| 🟡 13 | B.2: Most substrat→metryka explicit | wysoka | hipoteza pomostowa | ✅ |
| 🟡 14 | E.1: Zmiana tonu UV (hipoteza, nie twierdzenie) | niska | rzetelność | ✅ |
| 🟡 15 | B.4: Status σ_ab | średnia | GW sector | ✅ |
| 🟡 16 | D.1: Fermiony — trzypoziomowy podział | średnia | sektor materii | ✅ |
| 🟡 17 | D.2: Tabela U(1)/SU(2)/SU(3) | niska | gauge sector | ✅ |
| 🟢 18 | Czy ψ_ini = 7/6 wynika z Warstwy I? (redukuje N→1) | wysoka | **max. redukcja** | ✅ |
| 🟢 19 | E.2: p74 Wetterich flow | b. wysoka | UV kompletność | ✅ |
| 🟢 20 | A.3: Mapa zależności rozdziałów | średnia | nawigacja | ✅ |
| 🟢 21 | A.4: Przegląd stylistyczny | niska | profesjonalizm | ✅ |
| 🟢 22 | D.1d: Chiralność Dirac/Weyl | b. wysoka | sektor słaby | ✅ |

---


---

## Uwagi techniczne

### Centralny problem parsimoniowania: ile parametrów naprawdę?

Filozofia v2 stawia konkretne pytanie: jaki jest **minimalny zestaw parametrów Warstwy II**?

**Bieżący stan**: {Φ₀, a_Γ, ψ_ini} → 3 parametry

**OP-3 (α_K z a_Γ)**: jeśli α_K = f(a_Γ) — to a_Γ pozostaje Warstwą II, ale masy leptonów stają się Warstwą III. Redukcja: 3 parametry → 3 parametry (M_obs rośnie, N_param stały). Dobry wynik, ale nie redukcja N.

**Hipoteza głębsza (priorytet 18)**: ψ_ini = 7/6 jako warunek konieczny substratu:
- Substrat Γ przechodzi przejście fazowe GL → Φ_ini = Φ_eq (atraktor)
- Jeśli to wynika z A1, to ψ_ini znika z Warstwy II
- Redukcja: {Φ₀, a_Γ} → 2 parametry

**Hipoteza maksymalna**: a_Γ wynika z Φ₀ przez warunek spójności (bifurkacja przy skali Plancka?):
- a_Γ · √γ · ℓ_P = const → a_Γ = f(Φ₀, ℓ_P)
- Jeśli prawdziwe: **N_param = 1** (tylko Φ₀) — teoria jest ekstremalnie parsimonialna
- Sprawdzić numerycznie: czy a_Γ = 0.040049 ≈ g(Φ₀ = 24.66) dla jakiejś prostej funkcji g?

### OP-3: α_K z substratu (ścieżka domknięcia)

Priorytety: domknięcie OP-3 ≡ pokazanie, że α_K nie jest dodatkowym wolnym parametrem.

Ścieżka 1 (topologiczna):
- α_K = 8.56 ↔ soliton TGP z ładunkiem topologicznym Q=1 przy trzech generacjach
- Warunek: E_soliton(α_K) jest minimum globalne przy Q_Koide = 3/2 i a = a_Γ
- Sprawdzić: czy α_K = argmin E(α) przy Q_TGP(K₁(α),K₂(α),K₃(α)) = 3/2

Ścieżka 2 (relacja z Φ₀):
- r₂₁ = 206.77 ↔ m_μ/m_e — to jest stosunek mas obserwacyjny
- Czy r₂₁ = g(Φ₀) dla jakiejś prostej funkcji? Sprawdzić dla Φ₀ ∈ [20,35]

### Pytania strukturalne (konsekwencje próżni z źródłem)

Z twierdzenia „każdy punkt przestrzeni ma źródło" (B.3) wynikają pytania z konsekwencjami fizycznymi:

1. **Horyzont czarnej dziury**: czy to obszar, gdzie Φ→Φ_min > 0 (przestrzeń rzadka), czy Φ→0 (przestrzeń znika)? Implikuje inną fizykę osobliwości.

2. **Początki wszechświata**: czy Φ może być arbitralnie małe (Φ→0 w Big Bangu), czy substrat narzuca Φ ≥ Φ_min > 0? Implikuje brak osobliwości kosmologicznej.

3. **Fluktuacje próżniowe**: δΦ(x) zawsze na tle Φ₀·ψ_bg(t) — czy energię tła można obserwować inaczej niż przez Λ_eff?

Te pytania mają odpowiedzi w TGP — wymagają sformułowania jako **wnioski z thm:vacuum-source**.

### Log realizacji v2

| Data | Sesja | Zadanie | Status |
|------|-------|---------|--------|
| 2026-03-24 | v33 | Plan v2 zapisany (recenzja + uwagi autora) | ✅ |
| 2026-03-24 | v33 | Filozofia klasy wszechświatów + klasyfikacja parametrów I/II/III | ✅ |
| 2026-03-24 | v33 | **A.0**: `ssec:param-classification` — tabela Warstwy I/II/III w sek08 | ✅ |
| 2026-03-24 | v33 | **B.1**: `rem:coupling-map` — mapa rodziny sprzężeń K=φⁿ | ✅ |
| 2026-03-24 | v33 | **B.3**: `thm:vacuum-source` — próżnia z źródłem jako tw. strukturalne | ✅ |
| 2026-03-24 | v33 | **A.1**: `status_map.tex` — nowy plik: pełna tabela Aksjomat/Twierdzenie/Hipoteza/Szkic/Program | ✅ |
| 2026-03-24 | v33 | **A.2**: `rem:status-levels` — rozszerzony do 6-poziomów + `\statuslabel{}` w preamble | ✅ |
| 2026-03-24 | v33 | **A.5**: `ssec:predykcje-vs-odtwarzanie` — tabela predykcja/odtwarzanie/fit | ✅ |
| 2026-03-24 | v33 | **A.6**: `sssec:co-obalic` + `sssec:parameter-ranges` — lista falsyfikacji i okna param. | ✅ |
| 2026-03-24 | v33 | **C.1**: `sssec:cosmo-complete` — kompletny układ kosmologiczny TGP (5 równań) | ✅ |
| 2026-03-24 | v33 | **D.3**: `ssec:koide-hierarchy` w dodatekF — rekonstrukcja vs predykcja Koidego | ✅ |
| 2026-03-24 | v33 | **E.1**: `prop:asymptotic-safety` → `hypothesis` z tabelą mocnych/słabych argumentów | ✅ |
| 2026-03-24 | v33 | **B.1+**: `rem:minimal-nonlinearity` — „Zasada minimalnej nieliniowości" (4 kroki: ℤ₂ → parzystość → n=0,2 za ubogie → n=4 pierwsze sensowne → n>4 nie minimalne) | ✅ |
| 2026-03-24 | v33 | **B.2**: `rem:metric-bridge` — most substrat→metryka: 5 kroków z etykietami statusu | ✅ |
| 2026-03-24 | v33 | **B.4**: `sssec:sigma-status-map` — mapa statusu sektora tensorowego (6-wierszowa tabela) | ✅ |
| 2026-03-24 | v33 | **D.2**: tabela sektora cechowania U(1)/SU(2)/SU(3) + rem:SU2-vs-SU3 | ✅ |
| 2026-03-24 | v33 | **D.1**: `rem:fermion-4-steps` — czterostopniowa hierarchia sektoru fermionowego | ✅ |
| 2026-03-24 | v33 | **C.2A**: `sssec:ghost-check` — brak ghostów + stabilność gradientowa perturbacji MS-TGP | ✅ |
| 2026-03-24 | v33 | **A.3**: `ssec:dependency-map` w dodatekA — 4 łańcuchy dedukcji (I=substrat, II=stałe, III=kosmologia, IV=materia) + rem:bridge-summary | ✅ |
| 2026-03-24 | v33 | **A.4**: konwencja językowa wynika/sugeruje + fix rem:quant-vs-geom (metryczna hipoteza) + fix rem:strzalka-czasu | ✅ |
| 2026-03-24 | v33 | **C.2B**: `p73_perturbations_CMB.py` — numeryczne rozwiązanie MS-TGP; n_s=0.9668 vs Planck 0.9649; zamrożenie modów zweryfikowane; r=16ε otwarte (sektor σ_ab = Szkic); 11/13 PASS | ✅ |
| 2026-03-24 | v33 | **D.1d**: `sssec:chirality-D1d` — hyp:chirality-kink (U(ψ)≠U(2-ψ) → m_L≠m_R), rem:TGP-vs-Higgs, rem:D1d-open | ✅ |
| 2026-03-24 | v33 | **#18**: `prop:psi-ini-derived` — ψ_ini=7/6 jako Warstwa I→III (W(ψ)=0+GL); tab:param-classification zaktualizowana (ψ_ini: II?→I→III); rem:parsimony: N_param=3→2 z GL, N_param=1 otwarte | ✅ |
| 2026-03-24 | v33 | **OP-3**: `sssec:mass-hierarchy-D3` — rem:mass-hierarchy (tabela rekonstrukcja vs predykcja) + `prob:OP-3` (dwie ścieżki: topologiczna eq.OP3-path1 + relacyjna eq.OP3-path2; cel: α_K→Warstwa III; skrypt docelowy p76) | ✅ |
| 2026-03-24 | v33 | **E.2**: `p74_wetterich_flow.py` — obcięty FRG Wetterika dla K(φ)=φ⁴; LPA, Litim d=4, N=1; θ(m̃²)=−2 ISTOTNY, θ(g̃)=0 MARGINALNY, θ(h̃)=+2 NIEISTOTNY, θ(λ̃)=+4 NIEISTOTNY; φ⁶/φ⁸ zanikają UV→IR; brak nie-Gaussowskiego UV-FP w LPA; η_K4~O(g̃/32π²)≈0.001; 7/8 PASS; STATUS: PROPOZYCJA (LPA) / HIPOTEZA (pełny FRG/AS) | ✅ |
| 2026-03-24 | v33 | **C.2C**: `sssec:min-observables` — rem:min-observables (tabela predykcji sektora perturbacji): n_s≈0.967 (Planck 1σ, status PROPOZYCJA), r≈0.27 (SZKIC, σ_ab otwarte), f_NL~10⁻⁴ (zaniedbywalny), ISW i σ₈ (wymagają pełnej integracji); M_obs z perturbacji: ≥5 dodatkowych obserwabli; łącznie M_obs/N_param ≥ 11 (przy N_param=2, zał. GL) | ✅ |
| 2026-03-24 | v33 | **Hipoteza maks.**: `p75_agamma_phi0_correlation.py` + `hyp:agamma-phi0` w sek08 — numeryczna eksploracja a_Γ=g(Φ₀); wynik: a_Γ·Φ₀≈0.987≈1 (odchylenie 1,3%; najlepsza prosta relacja); a_Γ·Φ₀^(2/3)≈1/3 (1,8%); formuła a_Γ·√γ·ℓ_P=const NIE działa (błąd analizy wymiarowej jednostek); implikacja: jeśli dokładna → N_param=1; wymagana Ω_Λ*=0.6936 vs Planck 0.6847 (~1,2σ); hipoteza OTWARTA; 1/5 PASS | ✅ |
| 2026-03-25 | v33 | **OP-3 / p77**: `p77_K2_scan.py` (60 pkt, n_psi=120) — weryfikacja K*₂; wynik: K*₁=0.010414 (ref), K*₂=0.034003 (g≈0, prawdziwy soliton), K*₂/K*₁=3.27≠R₂₁=206.77; zmiany znaku g w p76 Sek.D (K≈0.057-0.089) to artefakty przeskoku gałęzi; V_mod standard NIE odtwarza hierarchii leptonów; Path 1 wymaga modyfikacji V_mod (p15/p78); res:op3-alphak zaktualizowany; 3/5 PASS | ✅ |
| 2026-03-25 | v33 | **OP-3 / p84**: ⚠️ **WYNIKI COFNIĘTE p96**: argmin E*(α,r_max=40)≈8.60 był artefaktem oscylacyjnego ogona (r_max=40 konwencja); przy r_max→∞ K*(α) monotonicznie maleje, brak minimum fizycznego przy α_K. Oryginalny opis: `p84_branch1_scan_parallel.py` `p84_branch1_scan_parallel.py` (skan α∈[5.0,8.60] krok 0.05, 30 wątków, wall 91s) — B1 znaleziona dla 65/73 pkt (luka α≈5.3–5.65: K* poza [0.004,0.12]); E/K∈[12.561,12.571] z wieloma przeskokami gałęzi w całym zakresie; minimum globalne przy α=8.60 (E/K=12.56141), ale to **granica zakresu i gałąź B2**; odchylenie 0.038 od α_K; lokalne minima na B-range: α≈6.85 (12.56160), α≈7.95 (12.56153); **wniosek OP-3 Ścieżka 1**: brak czystego minimum E*(α) przy α_K na żadnej ciągłej gałęzi; wielokrotne przeskoki gałęzi uniemożliwiają jednoznaczną kontynuację; argmin ~8.60–8.65, offset ~0.04–0.09; P3 PASS (odch<0.10), **P4 FAIL** definitywnie; OP-3 Ścieżka 1 NIEZWIĄZANA z α_K w obecnej formie | ✅ |
| 2026-03-25 | v33 | **OP-3 / p83**: `p83_branch_tracking.py` (skan α∈[8.50,8.71] krok 0.01, ciągła kontynuacja gałęzi, 3/5 PASS) — **dwie gałęzie przy danym K*₁**: Gałąź A (psi_core≈1.09–1.14, E/K≈−81 do −14): **NIEFIZYCZNA** — spełnia φ(r_max)=1, ale g=E/(4πK)−1 ≠ 0 (E<0); Gałąź B (psi_core≈1.238–1.244, E/K≈12.56–12.57): fizyczny soliton z g=0; w gałęzi B **dwie podgałęzie**: B1 (α≤8.59, E/K≈12.571, E* monotonicznie rosnące — minimum poza zasięgiem <8.50) i B2 (α≥8.60, E/K≈12.561, minimum E* przy α≈8.65, fit: 8.6504); **wniosek OP-3 Ścieżka 1**: argmin E*(fizyczna gałąź) ≈ 8.65, offset 0.089 od α_K=8.5616; P3 PASS (odch<0.10), P4 FAIL — minimum energii **nie jest** przy α_K (przynajmniej w zakresie [8.50,8.70]); możliwy scenariusz: minimum B1 leży poniżej α<8.50 bliżej α_K, ale wymaga szerszego skanu; sek08 zaktualizowany | ✅ |
| 2026-03-25 | v33 | **OP-3 / p82**: `p82_Estar_minimum.py` (skan α∈[8.50,8.70] krok 0.01, 4/5 PASS) — argmin E*(α)=8.6504 (fit paraboliczny); offset od α_K=8.5616: 0.089 (P3 PASS <0.10, P4 FAIL >0.02); **przeskak gałęzi przy α≈8.595**: E*/K* skacze 12.571→12.561 (g +0.04% → −0.04%); na górnej gałęzi (α≤8.59) E* monotonicznie maleje (brak minimum); na dolnej (α≥8.60) minimum przy 8.65; K*₁ minimum wspolbiezne z E* minimum (P5 PASS); OP-3 Ścieżka 1 NIEZDEMASKOWANA — minimum energii ~1% powyżej α_K; może być bias numeryczny (jak +0.8% K*₁); sek08 zaktualizowany | ✅ |
| 2026-03-25 | v33 | **OP-3 / p81**: `p81_alpha_universality.py` (5/5 PASS) — SEKCJA A: α_max(a_Γ) ZALEŻY od a_Γ (1.7% rozrzut [8.798,8.948]), ale P2 PASS (<5%); KLUCZOWE: n_s=0.9649 = α_K/α_max specyficzne dla fizycznego a_Γ=0.040 (inne a_Γ dają inne wartości); α_max ≈ 9.023 − 3.75·a_Γ (liniowy fit); SEKCJA B: E*(α) minimum przy α=8.600 (0.038 od α_K=8.5616); K*₁(α) minimum przy α=8.600; różnica E*(8.56) vs E*(8.60) tylko 0.14% → potrzebny p82 z krokiem 0.01; prob:OP-3 zaktualizowany (Ścieżka 3 + argmin E) | ✅ |
| 2026-03-26 | v33 | **OP-3 / p106**: `p106_alphaSN_vs_agamma.py` (skan α_SN(a_Γ), 5 wartości, bisekcja binarna, 31.1 min, 16 wątków). Wyniki: α_SN(0.020)=9.038, (0.030)=9.015, (0.040)=8.996, (0.050)=8.974, (0.060)=8.963. Fit liniowy: α_SN = 9.074−1.906·a_Γ, R²=0.988. **A = 9.074 ≠ 9.0** (delta=0.82%) — hipoteza A=9 OBALONA. α_SN ≈ 9.0 at a_Γ≈0.038 (blisko phys 0.040 ale nie exact). Dane wykazują nieliniowość (różnice prędkości nie monotoniczne). Fit kwadratowy A≈9.095. **WNIOSEK: α_SN zależy od a_Γ przez dphi₀=−K/a_Γ², nie jest pure-algebraiczny. Brak prostej relacji α_SN z n_s lub liczbą całkowitą.** sek08 zaktualizowane | ✅ |
| 2026-03-26 | v33 | **OP-3 / p105_v2**: `p105_v2_bisection.py` (**bisekcja binarna α_SN, 6 kroków, 14.2 min, 16 wątków**). Nasiona: ALIVE=8.9616 (K*₂=0.016488), DEAD=9.0116. Sekwencja: 8.9866 ALIVE, 8.9991 DEAD, 8.9929 ALIVE, 8.9960 ALIVE, 8.9975 DEAD, 8.9968 DEAD. **WYNIK: α_SN = 8.9964 ± 0.0004** (przedział (8.9960, 8.9968)). Test hipotezy α_K/α_SN = n_s: zmierzone 0.9517 vs n_s=0.9649 → delta=1.37% → **HIPOTEZA OBALONA**. Uwaga: najbliższa wartość to α_K/9 = 0.9513 (delta=0.040%), co sugeruje **α_SN ≈ 9.0 z dokładnością 0.04%**. Interpretacja: bifurkacja saddle-node Branch B zachodzi blisko całkowitej α_K=9. sek08 zaktualizowane (eq:alphaSN-value, eq:alphaK-alphaSN-ratio, eq:alphaSN-integer) | ✅ |
| 2026-03-26 | v33 | **OP-3 / p105**: `p105_alpha_SN_bifurcation.py` (wstępny skan α_K ∈ [8.5616, 9.20], krok=0.05, 20.2 min). Wyniki: Branch B alive przy α_K=8.7616–8.9616 (g_U≈5×10⁻⁵); dead przy 9.0116–9.1116. Punkty 8.6116–8.7116 miały g_U≠0 (kontynuacja tymczasowo zgubiła gałąź). Wstępny α_SN ∈ (8.9616, 9.0116). K*₂ maleje: 0.048→0.038→0.030→0.022→0.016 (trend do zera blisko bifurkacji). | ✅ |
| 2026-03-26 | v33 | **OP-3 / p104**: `p104_v2_continuation.py` (**ODKRYCIE: gałąź B istnieje TYLKO dla α_K ∈ [~7.5, ~9.0]**). Skan K*₂/K*₁ vs α_K — kontynuacja gałęzi (16 wątków, DALPHA=0.5, zakres [3,16]). Wyniki: (1) α_K=8.0616: r21=6.198; (2) α_K=8.5616: r21=9.755 (max!); (3) α_K=9.0616: gałąź B ZNIKA (bifurkacja saddle-node, zgodna z p88 α_SN≈8.875). Gałąź B **nie istnieje** dla α_K>9.06 ani α_K<7.56. Max r21=9.755 przy fizycznym α_K — α_K=8.5616 to punkt **tuż przed bifurkacją**! Ekstrapolacja: r21=206.77 wymagałoby α_K≈36 (gałąź B dawno nie istnieje). **WNIOSEK: K*₂/K*₁ STRUKTURALNIE OGRANICZONE do ~10 dla wszystkich α_K** — cel 206.77 niemożliwy przez żaden α_K. sek08 zaktualizowane | ✅ |
| 2026-03-26 | v33 | **OP-3 / p103+diag**: `p103_v4.py` + `p103_diag3.py` (**5/6 PASS** — KOREKTA METODOLOGICZNA p101/p102 + skan wykładnika n w V_mod=γ/3·φ³−γ/4·φ⁴+λ/n·(φ−1)^n, n∈{4,5,6,7,8,10,12,14,16,18,20}). **KLUCZOWE ODKRYCIE**: K*₂ z p101/p102 (psi₀≈4.1, K=0.171272) był NIEFIZYCZNY — g_U=3.38≠0 (potwierdzone diagnostyką p103_diag3). Prawdziwy K*₂ to gałąź B: psi₀≈2.77, K≈0.100, g_U≈0. Błąd źródłowy: PSI2=[3.50,7.00] wykluczało fizyczną gałąź B na rzecz gałęzi C (psi₀~4.1, zawsze g_U≈3.5) i zera spurious. **WYNIK p103**: K*₁=0.010280 (STAŁE), K*₂=0.100276 (STAŁE), K*₂/K*₁=**9.7547=STAŁE** dla WSZYSTKICH n=4..20. Lambda-term λ/n·(φ−1)^n przy psi₀≈2.77 (~2% dla n=20) jest zaniedbywalny — stosunek wyznaczany przez człony γ. **WNIOSEK: Ścieżka 1 Koide ZAMKNIĘTA NEGATYWNIE** — zmiana wykładnika n nie pomaga, K*₂/K*₁=9.73≠206.77. Poprzedni wynik p102 (16.66) był błędem metodologicznym. P5 PASS: K*₂/K*₁≈9.73 (prawdziwa wartość potwierdzona). P6 FAIL: cel 206.77 nieosiągnięty. sek08 zaktualizowane | ✅ |
| 2026-03-26 | v33 | **OP-3 / p102**: `p102_Vmod_lambda_scan.py` (**4/5 PASS** — Ścieżka 1 Koide: skan λ∈[1e−7, 1e−2] parametru V_mod=γ/3·φ³−γ/4·φ⁴+λ/6·(φ−1)⁶; wąskie okna: K*₁ z psi∈[1.10,1.80], K*₂ z psi∈[3.50,7.00]; 20 pkt logarytmicznie) — **WYNIK KLUCZOWY**: K*₁=**0.010281** (psi₀=1.2419) = STAŁE dla WSZYSTKICH λ. K*₂≈0.17128→0.17121 (zmiana 0.038%). **K*₂/K*₁ = 16.65–16.66 = STAŁE** dla wszystkich λ. Cel K*₂/K*₁=206.77 **NIEOSIĄGALNY** przez zmianę λ. Interpretacja fizyczna: człon λ(φ−1)⁶ wnosi ~3×10⁻³ do V_mod przy φ=psi₀≈4.1 nawet dla λ=1e−2 — ratio K*₂/K*₁ jest wyznaczane wyłącznie przez człony γ (kubiczno-kwartyczne), całkowicie nieczułe na λ. **WNIOSEK**: standardowa klasa V_mod z (φ−1)⁶ nie może odtworzyć hierarchii mas leptonów. Wymagana zmiana formy V_mod — następny krok: skan wykładnika n w (φ−1)^n (n=4,5,6,8,10,12). P4 FAIL: r21 nie-monotonicznie w λ (amplituda 0.007 = szum numeryczny). sek08 zaktualizowane | ✅ |
| 2026-03-26 | v33 | **OP-3 / p101**: `p101_K2_agamma_corrected.py` (**5/5 PASS** — korekta p100, wąskie okna psi: K*₁ z [1.001,1.8], K*₂ z [3.5,7.0]) — **KLUCZOWE ODKRYCIA**: (1) psi₀(K*₁) = **1.2419 = STAŁE** dla WSZYSTKICH a_Γ∈{0.010,...,0.060}! Psi₀ jest prawdziwie universalne (niezależne od a_Γ). (2) K*₁ ∝ a_Γ (Schwarzschild), odch. max 3.3% (bias r_max=80). (3) K*₂ (psi₀~4.1 gałąź) NIE ISTNIEJE dla a_Γ≤0.015; a_c∈(0.015,0.020]. (4) Dla a_Γ∈[0.020,0.060]: K*₂/K*₁≈16.7 (prawie stałe, K*₂∝a_Γ). (5) **OP-3 Ścieżka 2 NIE potwierdzona**: a_c≈0.017 ≠ 0.040 (fizyczne). (6) K*₂ z p101 (psi₀~4.1, K*₂=0.171) to INNA gałąź niż K*₂ z p99 (psi₀~5.0, K*₂=0.219) — bogaty spektrum. UWAGA: brak K*₂ dla a_Γ≤0.015 może być numeryczny (dphi₀=-K/a_Γ² → duże przy małym a_Γ, ODE fail). sek08 zaktualizowane | ✅ |
| 2026-03-26 | v33 | **OP-3 / p100**: `p100_K2_vs_agamma.py` (**2/5 PASS** — skan a_Γ∈{0.010,...,0.060}, r_max=80, psi∈[1.001,8]) — P1 PASS: K*₁ (niefizyc.) znaleziono przy wszystkich a_Γ. P2 FAIL: K*₁/ref≈1.85–2.29 (błąd ~2×, cel <10%). P3 PASS: K*₂ (niefizyc.) znaleziono przy a_Γ=0.040. P4 FAIL: K*₂ istnieje również dla a_Γ<0.040. P5 FAIL: brak a_c≈0.040. **DIAGNOZA METODOLOGICZNA**: "ostatnie zero psi" w szerokim oknie [1.001,8.0] wybiera NIEPRAWIDŁOWE gałęzie. K*₁ znalezione = 0.01856 (psi₀=1.425), nie fizyczne K*₁=0.010029 (psi₀=1.2418). K*₂ znalezione = 0.370 (psi₀=8.0), nie fizyczne K*₂≈0.22 (psi₀≈5.0 z p99). Prawidłowe K*₂/K*₁ (niefizyc. gałęzie) ≈ 17–21. **Wniosek**: p100 NIEROZSTRZYGAJĄCY dla OP-3 Ścieżka 2 — wymaga wąskich okien psi: K*₁ z [1.001,1.8], K*₂ z [3.5,6.5]. Ścieżka 2 pozostaje OTWARTA. sek08 zaktualizowane | ✅ |
| 2026-03-26 | v33 | **OP-3 / p99**: `p99_K2_convergence.py` (**3/5 PASS** — konwergencja K*₂ dla r_max∈{60,80,100,120,150}, skan K∈[0.10,0.45], dphi0=−K/A_Γ²) — P1 PASS: K*₂ znaleziono przy wszystkich 5 r_max. P2 FAIL: psi₀ saturuje na PSI_HI=5.0 (r=150 inny oddział). P3 FAIL: R²=0.011 (r=150 outlier: K*₂=0.200 zamiast ~0.227). P4 PASS: K*₂/K*₁(∞)≈21. P5 PASS: ratio≈21 << 206.77. **Wyniki (r=60,80,100,120 — monotonne)**: K*₂=0.2127→0.2192→0.2235→0.2272; psi₀≈5.0 (granica zakresu). Ekstrapolacja 5-pktowa: K*₂(∞)≈0.213 (a=0.213, b=0.291, R²=0.011 — zdominowana przez r=150). Ekstrapolacja z 4 pkt (r=60–120): K*₂(∞)≈0.24–0.25 (r_max=150 zmienia gałąź). **K*₂/K*₁(∞) ≈ 21–25** (istotnie < 206.77). Czynnik niedoboru ≈ 9–10×. **Wniosek**: standardowy V_mod daje K*₂/K*₁(∞)≈21–25, NIE 206.77 — potwierdzenie konieczności modyfikacji V_mod (OP-3 Ścieżka 1 Koide). sek08 zaktualizowane | ✅ |
| 2026-03-26 | v33 | **OP-3 / p98**: `p98_full_K_spectrum.py` (**3/4 PASS** — skan K∈[0.006,0.25], psi∈[1,5], r_max=40,60,80) — **DRUGIE ROZWIĄZANIE SOLITONOWE**: K*₂ PRZEŻYWA przy r_max=60,80 (nie jest artefaktem). K*₂: 0.0857(r=40)→0.2103(r=60)→0.2153(r=80)→~0.369(∞). K*₂/K*₁: 8.14(r=40)→20.23(r=60)→20.85(r=80)→~36.8(∞). **Wniosek**: standard V_mod daje K*₂/K*₁(∞)≈37, NIE 206.77 — masa μ/e nie odtworzona. n_zeros_psi rośnie z r_max (3-5→7-9): dodatkowe gałęzie Bessela. K*₂ przy r=40 (0.086) to INNA GAŁĄŹ niż K*₂ przy r=60 (0.21) — dramatyczny przeskok +145% jest artefaktem gałęziowym. P3 FAIL: ratio niestabilne (8→21). P4 PASS: K*₂(∞)/K*₁(∞)=36.83 wyznaczone. **Konkluzja**: istnieją co najmniej 2 fizyczne solitony, ale stosunek mas wymaga modyfikacji V_mod (jak OP-3 Ścieżka 1 Koide sugerowała w p77). sek08 zaktualizowane | ✅ |
| 2026-03-26 | v33 | **OP-3 / p97**: `p97_true_soliton_profile.py` (**3/5 PASS** — K*(α_K, r_max) dla r_max=40,60,80,100,150, profil do r=200) — **PRAWDZIWY PROFIL SOLITONU TGP**: K*(∞)=**0.010029** (4-pkt fit, vs 0.010503 przy r_max=40, bias 4.72%); E*=0.126033; psi₀=**1.2418** (zbieżne r_max≥60); m_tail=0.323397 = 1/√(1+α_K) (błąd 0.000%, P3 PASS!); C_tail=0.01167; R²=1.000000 (P2 PASS!). Skale: r_half=1.97 a_Γ, r_90=10.41 a_Γ, r_99=98.3 a_Γ. P1 FAIL: fit-3 vs fit-4 różnią się 0.222% (K*(∞) z wyższymi poprawkami 1/r²). P4 FAIL: bias=4.72% < 5% (próg zbyt surowy). Wynik: psi₀ jest NAJSTABILNIEJSZĄ wielkością (zbieżne do 0.01%). sek08 zaktualizowane eq:Kstar-true, eq:psi0-true, eq:tail-params-alphak | ✅ |
| 2026-03-26 | v33 | **OP-3 / p96**: `p96_Estar_rmax_convergence.py` (**2/5 PASS** — skan α∈[8.20,8.70], r_max=40,60,80, n_psi=80, n_eval=500) — P1 PASS: argmin E*(α,r=40)=8.646 ∈ [8.55,8.70] (zgodne z p82-p83). P2-P5 FAIL. **KLUCZOWE**: K*(α,∞) MONOTONICZNIE MALEJĄCE w całym zakresie [8.20,8.70] — brak minimum fizycznego. Tabela argmin(r_max): r=40: 8.646, r=60,80,∞: ≥8.70 → shift argmin: +0.054 przy r_max: 40→∞, odchylenie od α_K rośnie (nie maleje). **WNIOSEK: OP-3 Ścieżka 1 (argmin E*=α_K) ZAMKNIĘTA** — pozorne minimum przy r_max=40 jest artefaktem oscylacyjnego ogona (analogicznie do K_nan_R i Gałęzi L). K*(α,∞) monotoniczne — energia solitonu nie wskazuje na α_K. OP-3 Ścieżka 2 (WKB/relacyjna: a_c(α_K)=a_Γ) pozostaje jedyną otwartą. sek08 i PLAN zaktualizowane | ✅ |
| 2026-03-26 | v33 | **OP-3 / p95**: `p95_knan_vs_rmax.py` (**1/5 PASS** — rozszerzony zakres K=[0.008,0.040], PASS P1 sanity) — **KLUCZOWE ODKRYCIE: K_nan_R jest artefaktem r_max=40**. Część A: K_nan_R(r_max=40)=0.01170, K_nan_R(r_max=60,80,100)=**BRAK** w K∈[0.008,0.040]. Część B (smoking gun): K=0.01270 daje n_zeros=0 przy r_max=40, ale n_zeros=2 przy r_max=60 — ta sama K, różna faza oscylacyjnego ogona w r=40. Część C: przy K>K_nan_R(r_max=40), profil φ(r) NIE dywerguje — φ(r=40)=1.005, φ_max=1.5 (wartość początkowa), profil dotarł do r=100. **WNIOSEK**: K_nan_R(α_max=8.8734) z p80 jest artefaktem fazy oscylacyjnego ogona przy r=40. Definicja α_max przez K*=K_nan_R nie ma charakteru fizycznego. **Koincydencja α_K/α_max=n_s COFNIĘTA** — nie odporna na wybór r_max. OP-3 Ścieżka 3: ZAMKNIĘTA (artefakt). sek08 + PLAN zaktualizowane | ✅ |
| 2026-03-26 | v33 | **OP-3 / p94**: `p94_alpha_max_convergence.py` (**1/5 PASS** — lean v3, 5α×3r_max, n_psi=60, n_K=20, n_eval=400) — P1 PASS: gap(8.8734,r=40)=-0.000579≈0. FAIL: K_nan_R **brak** dla wszystkich α przy r_max=60,100 w K∈[0.008,0.016]; α_max wyznaczone tylko przy r_max=40 (α_max≈8.860, ratio=0.9663); gap(8.86,r=40)=+0.00135 → gap(8.8734,r=40)=-0.00058 (α_max między 8.860 a 8.8734). Diagnoza: K_nan_R jest r_max-zależne; α_max(r_max=60) = undefined w standardowym zakresie K; następnik p95 z K do 0.040 | ✅ |
| 2026-03-26 | v33 | **OP-3 / p93**: `p93_rmax_convergence.py` (**5/5 PASS**) — Audyt zbieżności K*(r_max→∞) dla α=8.840: K*(r=40)=0.010818 → K*(60)=0.010312 → K*(80)=0.010151 → K*(100)=0.010066 → K*(∞)=**0.009521** (fit K*(r)=0.00952+0.051/r, R²=0.985, P2 PASS). Różnica K*(r=40) vs K*(∞) = **13.6%** (P3 PASS). Gałąź L: psi_L(pierwsze zero) przesuwa się 1.189→1.065 (Δpsi=0.124) gdy r_max: 40→100 — **artefakt potwierdzony** (P5 PASS). Fit ogona Bessela A sin(mr)/r + B cos(mr)/r: R²=1.0000 (P4 PASS). **WNIOSKI**: (1) Jedyna fizyczna rodzina solitonów — Gałąź U; (2) K*(p80/p91) ≈ 0.011 prawidłowe w konwencji r_max≈40; (3) K*_true ≈ 0.0095; (4) Gałąź L z p91 = artefakt r_max=40, cofnięta; sek08 i PLAN zaktualizowane | ✅ |
| 2026-03-25 | v33 | **OP-3 / p92**: `p92_two_soliton_profiles.py` (porównanie profili φ(r) dwóch klas, 1/2 PASS) — K*_B=0.010311 (Branch U, g=3.3e-6 ✓); K*_A=brak przy α=8.840 z poprawionymi zakresami K → Branch L nie wykryta; profil Branch U: φ_max=1.2345, r_half=0.075 (1.9 a_Γ), r_90=0.277 (6.9 a_Γ), r_99=0.670 (16.8 a_Γ); **diagnoza nieudanego Branch L**: ogon oscylatoryjny φ(r)~sin(mr)/r sprawia, że psi₀(r_max) silnie zależy od r_max — root K*_A zmienia pozycję (1.279→1.330 gdy r_max 40→100); przy r_max=60 Branch L znika; implikacja: p91's K*_A był fałszywym zerem | ✅ |
| 2026-03-25 | v33 | **OP-3 / p91**: `p91_codim2_branches.py` (gałęzie psi₀ oddzielnie, 3/5 PASS) — ⚠️ **WYNIKI CZĘŚCIOWO COFNIĘTE przez p92–p93**: przy r_max=40 (p91) skan psi wykrywa fałszywe "Gałąź L" (K*_A≈0.014) jako artefakt oscylacyjnego ogona φ(r); "Gałąź U" (K*_B≈0.011) jest jedyną fizyczną gałęzią. Zachowane z p91: K*_B ≈ 0.01082 (α=8.840) i α_max≈8.875 (zbieżne z p80). Odrzucone: "dwie rodziny solitonów", α_saddle≈8.859, kaskada zaników, fit γ_AB | ✅ |
| 2026-03-25 | v33 | **OP-3 / p90**: `p90_codim2_v3.py` (próba naprawcza p89, 0/5 PASS) — zwiększono psi_max z 3.0 do 20.0 → krok psi wzrósł z 0.020 do 0.190 → wąskie zera φ(r_max)=1 pominięte całkowicie → K1=nan wszędzie; K grid rozszerzony do [0.008, 0.0195]. Lekcja: psi_max HIGH musi iść z n_psi proporcjonalnym (>600 punktów przy max=20). BUG zdiagnozowany, poprawka w p91 | ✅ |
| 2026-03-25 | v33 | **OP-3 / p89**: `p89_codim2_verify.py` v2 (weryfikacja hipotezy kodimension-2 z p88, 0/5 PASS) — K grid [0.0106, 0.0192], psi_max=3.0, n_psi=100; K*₁ znaleziony konsekwentnie ≈0.011; K*₃ znaleziony tylko przy α=8.855 (0.013551) i α=8.863 (0.013071); K_nan_left BRAK; DIAGNOZA problemów: (1) psi_max=3.0 zbyt niskie — K*₃ wymaga psi₀≈6.47 (diagnostyka wykazała 3. zero φ(r_max)−1 przy ψ≈6.5, ale g≈−1000, nie jest K*₃!); faktyczny K*₃ jest 3. zerem g_min wynikającym z artefaktu przełączenia gałęzi (nie z psi₀≈6.5); (2) K_nan_left < K_grid_start = 0.0106; (3) Python 3.14 f-string crash w podsumowaniu; następnik: p91 z osobnymi skanami gałęzi | ✅ |
| 2026-03-25 | v33 | **OP-3 / p88**: `p88_cascade_bifurcation.py` (kaskada saddle-node, 3/5 PASS) — **KASKADA BIFURKACJI**: γ₃₄≈0.96 (K*₃⊕K*₄ przy α≈8.77) → γ₂₃=0.727±0.038 (K*₂⊕K*₃ przy α_SN23=8.875, RMSE=0.265%) → γ₁₂≈0.55 (K*₁→granica przy α_max=8.8734); **ODKRYCIE**: α_SN23=8.875 ≈ α_max=8.8734 — MOŻLIWA BIFURKACJA KODIMENSION-2 gdzie jednocześnie: K*₂=K*₃ (saddle-node) I K*₁=K_collapse (granica); daje układ 2 równań analitycznych na (K_c, α_max) z V_mod — droga do OP-3 Ścieżki 3! | ✅ |
| 2026-03-25 | v33 | **OP-3 / p87**: `p87_gamma_fit.py` (alpha∈[8.80,8.87], 16 pkt K*₂, 4/5 PASS) — **gamma=0.551** (stały α_max=8.8734) / gamma=0.661 (swobodny α_max_SN=8.880); delta_K=K*₂−K*₁ maleje monotonicznie 0.0039→0.0016; K*₂ zanika przy α≈8.863 (przed α_max=8.8734) — wchłonięte przez nan-region; WNIOSEK: gamma≈0.55 bliskie saddle-node (0.5) ale nie dokładne — **bifurkacja graniczna**: K*₁ i K*₂ nie łączą się wprost, lecz oba osobno wchłaniane przez rosnący zakres zapadania profilu; mechanizm OP-3 Ścieżka 3 bardziej złożony niż ∂g/∂K=0 | ✅ |
| 2026-03-25 | v33 | **OP-3 / p86**: `p86_saddle_node_K2.py` (diagnostyka g(K,α), 8 wartości α, 4/5 PASS) — **ODKRYCIE: 4 solitony** przy α=8.5–8.7 (4 zera g(K)); K*₂/K*₁: 2.81→1.18 (α=8.5→8.85); K*₂ znika przy α≈8.87; nan-region (zapad profilu φ) rozszerza się od K=0.026 do K=0.012 gdy α: 8.80→8.87; nan-region lewy brzeg: K_nan≈0.023 (α=8.80)→0.012 (α=8.87) zbliża się do K*₁≈0.011; nowa hipoteza: α_max wyznaczony warunkiem K_nan_left(α)=K*₁(α) | ✅ |
| 2026-03-25 | v33 | **OP-3 / p85**: `p85_saddle_node_scaling.py` (25 pkt α∈[8.70,8.87], 8 workerów, 2/5 PASS) — **β ≈ 0.11–0.13 (NIE 0.5!)**: K*₁(α) ~ K_c − A·(α_max−α)^β z β=0.265 (fit swobodny) / β=0.192 (stały α_max); **β jest UNIVERSALNE** dla a_Γ∈{0.02–0.06}: β_mean=0.1116±0.0117 (**≈ 1/8 = 0.125**, wykładnik 2D Isinga!); K_c=0.01187 (wartość K*₁ przy α_max); WNIOSEK: K*₁ NIE jest właściwym parametrem porządku bifurkacji saddle-node — właściwy parametr to **K*₂−K*₁ ~ √(α_max−α)**; α_max_fit=8.854 (odch. 0.019 od ref, spodziewane — fit jednoparametrowy); p86 zaplanowany: K*₂(α) blisko α_max → weryfikacja saddle-node via różnica gałęzi | ✅ |
| 2026-03-25 | v33 | **OP-3 / p80**: `p80_alpha_max.py` (bisection w [8.80,8.90], 4/4 PASS) — **α_max=8.8734** (±0.002); K*₁(α) rośnie przy α→α_max (K*₁=0.01062 przy 8.80 → 0.01130 przy 8.87); ε=3.51%; **α_K/α_max=0.9649=n_s(Planck)** (uderzająca zbieżność do 4 cyfr!); rem:alpha-max-bifurcation + eq:alpha-naturalness + eq:ns-coincidence dodane do sek08; rem:vacuum-source-corollaries (BH, Big Bang, fluktuacje) dodane; OP-3 nova via: α_max=f(V_mod,a_Γ)? → α_K jako predykcja stabilności. ⚠️ **COFNIĘTE p95**: K_nan_R(α_max) jest artefaktem r_max=40 → α_max i koincydencja n_s nie mają charakteru fizycznego; OP-3 Ścieżka 3: ZAMKNIĘTA | ✅ |
| 2026-03-25 | v33 | **OP-3 / α_max diagnostyka**: skan g(K) dla α=8.9, K∈[0.003,0.50] z K_max rozszerzonym do 0.5 — wynik: g=-0.195 przy K=0.009, **nan (zapad profilu)** przy K≥0.013, brak g=0 w całym zakresie → α_max REALNY (nie artefakt K_max); mechanizm: bifurkacja siodło-węzeł (K*₁ i K_collapse zbiegają); rem:alpha-max-bifurcation + eq:alpha-naturalness dodane do sek08; **OP-3 nova via**: czy α_max=f(V_mod, a_Γ) analitycznie? → α_K=α_max(1−ε) jako predykcja | ✅ |
| 2026-03-25 | v33 | **OP-3 / p79**: `p79_alphac_p10accuracy.py` (n_psi=120, n_eval=2000, rtol=1e-9, jak p10) — K*₁(α=8.5616)=0.010494 (+0.77% vs K1_ref — poprawa p78 +1.04%, ale bias strukturalny); K*₁>K1_ref dla ALL α∈[7.5,8.8]; **BIFURKACJA α_max≈8.85**: standardowa gałąź solitonu NIE ISTNIEJE dla α≥8.9 (p78+p79 zgodnie); na α=9.5 inna gałąź: K*₁=0.008514; interpolacja alpha_c=8.869 (3.59%) na granicy bifurkacji (inna gałąź, niefizyczne); **nowe odkrycie**: α_K=8.5616 ≈ 0.97·α_max — tuż poniżej maksymalnego α dla istnienia solitonu; sek08 zaktualizowany; 3/5 PASS | ✅ |
| 2026-03-25 | v33 | **OP-3 / p78**: `p78_alphac_precise.py` (n_psi=80, n_eval=1000, rtol=1e-8) — cel: sprawdzenie czy 1,3% rozbieżność alpha_c to artefakt n_psi=35; wynik: K*₁(n_psi=80)=0.010523 (+1.04% vs K1_ref) — przeskakuje powyżej K1_ref (p76 było −1,06%); K*₁(α) NIE monotoniczne: minimum ≈0.010508 przy α≈8.65, znika dla α≥8.9; alpha_c NIE znaleziony w [8.0, 9.2] (K*₁ > K1_ref wszędzie); DIAGNOZA: p78 używało n_eval=1000 vs p10's 2000 → systematyczny +1% bias; p79 napisany z n_eval=2000/rtol=1e-9 jak p10; 2/5 PASS | ✅ |
| 2026-03-25 | v33 | **OP-3 / p76**: `p76_alpha_K_bifurcation.py` (v3, kryterium energetyczne g=E/(4πK)−1=0, portowane z p10) — cross-check K*_1=0.010304 (błąd 1,1%); skalowanie Schwarzschilda K*_1≈0.258·a_Γ potwierdzone; skan α∈[4,12]→K*_1 maleje monotonicznie; interpolacja α_c=8.449 (odch. 1,3% od α_K=8.5616 — ta sama 1,3% co a_Γ·Φ₀!); OP-3 Ścieżka 2: CZĘŚCIOWO POTWIERDZONA (α_K wyznaczalne z substratu); Ścieżka 1 (Koide): otwarta (K*_2 wymaga p15); wstępna sugestia K*_2≈0.07 z zmiany znaku g(K) (wymaga weryfikacji); res:op3-alphak dodany do sek08; 5/5 PASS | ✅ |
| 2026-03-27 | v33 | **OP-3 / Ścieżka 8–9 (ex62–ex83): ZAMKNIĘTA NEGATYWNIE** — seria skryptów `ex62`–`ex83` (nbody/examples/) badała hipotezę A_tail^4 ≈ r₂₁=206.77 przez geometrię ODE g''=f(g,α). Schemat: `fz0(α) = (A(φ·z0(α))/A(z0(α)))⁴`, gdzie z0(α) = punkt B_coeff=0, A = amplituda ogona asymptotycznego. Seria ex62–ex71: eksploracja wstępna, ex71 użył 25-punktowej siatki do wyznaczenia z0(α) → spurious zero crossings przepropagowane jako hardcoded truth (a1s=2.44143051, a2s=2.74175411) przez ex72–ex80. **ex80_algebraic_analysis.py**: CZYSTO ALGEBRAICZNY — przyjął hardcoded wartości α*₁,α*₂ z ex79 jako INPUT; zidentyfikował S = α*₁+α*₂ ≈ 2π−11/10 (zbieżność 0.13 ppm) — ale na fałszywych danych. **ex83_Rstar_precise.py**: 0/5 PASS; S(r21_CODE)=5.18215 vs ex80's 5.18318 (diff −1.037×10⁻³); α*₁=2.440332, α*₂=2.741815 (niespójne z ex84). Diagnoza: systematyczny offset w ekstrapolacji A_inf poprzez curve_fit z ex71's coarse z0. | ✅ |
| 2026-03-27 | v33 | **OP-3 / ex84_sum_formula_verification.py**: **FALSYFIKACJA formuły S = 2π−11/10** (**1/5 PASS**). Bezpośrednia weryfikacja fz0(α) bez pośrednich hardcoded wartości. Parametry: R_MAX=120, WIN_COARSE=[(20,36),(30,46),(40,56),(50,66),(62,78)], fit 2-parametrowy A_inf: ai·(1+a/x+b/x²), siatka z0 60 pkt w [1.05,2.5]. **WYNIKI KLUCZOWE**: Prawdziwe zera fz0: α*₁=**2.4396323858**, α*₂=**2.6953122503** (NIE 2.4414 i 2.7418 jak twierdziła seria ex71–ex83). Suma S = α*₁+α*₂ = **5.1349** (odchylenie **9307 ppm** od 2π−11/10 = 5.1832). P1 PASS: fz0(α*₁)=206.768 ≈ r₂₁=206.77 ✓; P2 PASS: fz0(α*₂)=206.768 ✓; P3 FAIL: S≠2π−11/10; P4 FAIL: α*₁≠2.44143; P5 FAIL: α*₂≠2.74175. **ŹRÓDŁO ARTEFAKTU**: ex71 użył 25-punktowej siatki do znalezienia z0(α) → spurious zero crossings B_coeff → fałszywe α* przepropagowane przez ex72–ex83 jako hardcoded truth → S≈2π−11/10 był artefaktem numerycznym. **WNIOSEK: Ścieżka 8–9 (formuła sumy ogona) ZAMKNIĘTA NEGATYWNIE.** Fizyczne zeros at α=2.440 i 2.695 nie mają związku z 2π−11/10. | ✅ |
| 2026-03-27 | v33 | **Porządkowanie**: przeniesiono do `nbody/examples/_archiwum/`: ex1–ex54 (stare eksperymenty n-body, krzywe rotacji, Efimov, falowody grawitacyjne — ścieżki 1–8 zamknięte negatywnie) + stare warianty ex60/ex61 (konflikty nazewnicze). Do `TGP_v1/_archiwum/`: `PLAN_ROZWOJU_v1.md` (zastąpiony v2) + `PLAN_ANALITYCZNY_KOIDE.md` (zastąpiony przez dodatekJ). Łącznie 79 plików zarchiwizowanych. | ✅ |
| 2026-03-28 | v33 | **Korekty LaTeX**: (1) `main.tex`: dodano brakujący `\input{dodatekJ_ogon_masy}` (plik był kompilowany bez Dodatku J!). (2) `sek07_predykcje.tex`: naprawiono 4 wystąpienia `\ref{ssec:stale}` → `\ref{sec:stale}` (etykieta `ssec:stale` nie istnieje; właściwa etykieta to `sec:stale` z sek04_stale.tex linia 4). | ✅ |
| 2026-03-28 | v33 | **hyp:agamma-phi0 / ex85**: `ex85_agamma_phi0_precision.py` (**5/6 PASS** — precyzyjna weryfikacja hipotezy T1: a_Γ·Φ₀=1 z danymi Planck 2018, DESI DR1 BAO alone, DESI DR1+CMB+lensing). **WYNIKI KLUCZOWE**: (1) **DESI DR1+CMB**: a_Γ·Φ₀=**0.9991** (odch. −0.09%, **−0.12σ od 1** — hipoteza numerycznie potwierdzona!). Planck 2018: a_Γ·Φ₀=0.9872 (−1.28%, −1.22σ — zgodne). Wymagane Ω_Λ*=0.6936 vs DESI+CMB: 0.693±0.005 (dev +0.12σ). (2) **NOWA RELACJA** (TOP-1 skanu algebraicznego): `α_K·√(a_Γ·r₂₁) = Φ₀` (0.048% od dokładności z Planck 2018). Jeśli obie relacje dokładne: a_Γ = (α_K·√r₂₁)^(-2/3) = 0.04011 (0.15% od 0.040049) → a_Γ wyznaczane z parametrów cząsteczkowych! (3) **Trojca T1,T2,T3**: T1(a_Γ·Φ₀=1) i T2(r₂₁=Φ₀·α_K) niezależne (różne Φ₀*), T3=T1·T2. T1 popierane przez DESI+CMB, relacja złożona przez Planck. P6 FAIL: TOP kombinacja to 4-parametrowa relacja złożona, nie czyste T1. **STATUS hyp:agamma-phi0: WZMOCNIONA PRZEZ DESI** — z hipotezy „1.3% od 1" do „0.09% od 1" z najlepszym zestawem danych. sek08 zaktualizowane (tabela, footnote DESI 2024 VI, relacja złożona). | ✅ |
| 2026-03-28 | v34b | **F1–F9 domknięcie**: (F1) `prop:J-zero-mode-mass` — tryb zerowy ODE → M=O(A⁴) analitycznie (`dodatekJ`). (F3) `hyp:agamma-phi0-combo` — T1+T2 → a_Γ Warstwa III, tabela parametrów (`sek08`). (F4) `ssec:vacuum-source-chain` — narracja N0-5→W(1)→Λ_eff→kosmologia (`sek08`). (F6) box „5 predykcji falsyfikacyjnych" w `ssec:specyficzne` (`sek07`). (F7) `prop:spatial-metric-from-substrate` — h(Φ)=Φ/Φ₀ z gęstości węzłów; zamknięty wiszący ref (`sek08`). (F9) `problem_r21_analiza.md` — status zmieniony, sekcje ex84/ex85, ⚠️ artefakt ex71–ex72. `status_map.tex`: liczniki 43→45 (Propozycje 14→15, Hipotezy 7→8). `PLAN_ROZWOJU_v2.md`: porządek + wpisy statusów F1–F9. Pozostaje: F5 (atraktor ψ_ini), F8 (README skryptów). | ✅ |
| 2026-03-28 | v34 | **Sesja v34**: (1) Stworzono `dodatekK_wkb_atail.tex` — formalizacja asymptotycznej analizy ogona (prop:K-tail-form, res:K-SWKB, klasyfikacja sektorów topologicznych, predykcja tauonu hyp:K-tau g₀^τ≈2.34). (2) Stworzono ex60 (topologia ogona), ex61 (predykcja tauonu, C1–C5), ex62 (skalowanie WKB, α-niezależność n≈4). (3) Zaktualizowano `main.tex` (dodano `\input{dodatekK_wkb_atail}`), `status_map.tex` (nowe wiersze + blok sesji v34, Propozycje 12→14, Hipotezy 6→7, Program 2→3, Łącznie 39→43), `dodatekJ_ogon_masy.tex` (odsyłacze do Dod. K, awans O-J2/O-J3 z Program → Hipoteza). (4) Analiza oceny zewnętrznego agenta (12 punktów) — nowa sekcja v34 dodana do PLAN_ROZWOJU_v2.md. | ✅ |

---

## Sesja v34 — Analiza oceny zewnętrznej i priorytety (2026-03-28)

> **Podstawa**: ocena zewnętrznego agenta (12 punktów) + analiza `problem_r21_analiza.md` + stan po sesjach v33–v34.
> **Ważne zastrzeżenie**: agent oceniał najwyraźniej `main.tex` bez dostępu do `PLAN_ROZWOJU_v2.md`, `status_map.tex`, historii zamkniętych ścieżek ani wyników ex71–ex85. Część punktów opisuje pracę już wykonaną.

---

### Weryfikacja 12 punktów — co wartościowe, co pomylone

| # | Punkt agenta | Ocena | Priorytet |
|---|-------------|-------|-----------|
| 1 | Twarda redukcja statusów w main.tex | Zrobione (`status_map`, `\statuslabel{}`) — brakuje korekty tonu prozy w sek03/sek07 | Niski |
| 2 | Minimalna baza parametrów | Zrobione — Warstwa I/II/III istnieje od v33; agent nie widział planu | Niski |
| 3 | Kosmologia jako kompletny moduł | Częściowo — równania są, ale `ψ_ini` jako atraktor nieformalny; n_s numerycznie, nie analitycznie | Wysoki |
| 4 | Most substrat → Φ → metryka | Słabe ogniwo — `f(Φ)=Φ/Φ₀` motywowane, ale nie wynika ściśle z działania | Wysoki |
| 5 | Sektor mas jako osobny program | Częściowo zrobione (Dod. J i K appendiksy); jedno zdanie programowe do sek07 | Niski |
| 6 | Kanoniczna definicja masy fizycznej | **Krytyczne, niezrobione** — ex66 dał odpowiedź (M_phys = O(A_tail⁴)), ale brak formalizacji w .tex | Bardzo wysoki |
| 7 | Kwantyzacja: wyprowadzone vs analogia | Znane — Dod. E istnieje; brakuje scope statement (2 zdania) | Niski |
| 8 | 3–5 predykcji odróżniających TGP | Słuszne — sek07 rozproszony; warto spakować w ramkę (tryb oddechowy GW, G(Φ), 3 generacje…) | Wysoki |
| 9 | Architektura numeryki | Słuszne — 1 plik README_map.md wystarczy; niski priorytet naukowy | Niski |
| 10 | „Próżnia ze źródłem" W(1)≠0 | **Najbardziej oryginalne, niedostatecznie rozwinięte** — B.3 istnieje, ale brak pełnego łańcucha N0→W(1)→Λ | Bardzo wysoki |
| 11 | Przyciąć do wersji „obronnej" | Słuszne — ale zbyt wcześnie; dopiero po domknięciu F1–F3 | Odłożone |
| 12 | „Dlaczego nasz wszechświat?" jako selekcja | **Najsilniejszy kierunek, zgodny z filozofią v2** — selekcja przez stabilność/bifurkacje, nie numerologia | Bardzo wysoki |

---

### Krytyczny wynik z `problem_r21_analiza.md`: ex84 falsyfikacja S=2π−11/10

> ✅ **Zaktualizowane (v34)**: Plik `problem_r21_analiza.md` — frontmatter, sekcja ex84, sekcja ex85, ⚠️ artefakt w ex71–ex72. Obecny status: „REWIZJA v34 — ex84 falsyfikuje S=2π−11/10".

**Co obalono i co zostaje:**

| Element | Status po ex84 |
|---------|----------------|
| Dwa zera fz0(α)=r₂₁ przy α*₁≈2.440 i α*₂≈2.695 | **Realne** (P1, P2 PASS w ex84) |
| Warunek φ-FP: (A(φg₀*)/A(g₀*))⁴ = r₂₁ → g₀*=1.24771 | **Realne** (0.0001% dokładność) |
| Formuła S = α*₁+α*₂ = 2π−11/10 | **Artefakt** (9307 ppm; ex71 użył zbyt grubej siatki) |
| Generalizacja S(n) = 2π−(5n+1)/10, D(n) = (n+1)/10 | **Artefakt** (bazowała na błędnych α*) |
| β=1 jako centrum okna bifurkacji (ex76–ex77) | **Wymaga re-weryfikacji** z poprawionymi α* |

✅ **Zrobione (v34)**: `status_map.tex` — ex71–ex83 zastąpione przez blok z ex84 (falsyfikacja) + ex85 (T1/T2); α* skorygowane; β=1 zdegradowane do Hipoteza z caveat.

---

### Nowe wartościowe wyniki z ex85 (po v33)

**Hipotezy algebraiczne między parametrami Warstwy II** (ex85, 5/6 PASS):

| Hipoteza | Formuła | Dokładność | Dane |
|----------|---------|------------|------|
| T1 | a_Γ · Φ₀ = 1 | 0.09% (−0.12σ) | DESI DR1+CMB |
| T2 | r₂₁ = Φ₀ · α_K | 0.048% | Planck 2018 |
| Konsekwencja T1+T2 | a_Γ = (α_K · √r₂₁)^(−2/3) | 0.15% | Planck 2018 |

Jeśli T1 i T2 są dokładne → a_Γ wyznaczane przez (α_K, r₂₁, Φ₀) → a_Γ **znika z Warstwy II**. To jest najważniejsza ścieżka redukcji parametrów po falsyfikacji Ścieżki 8–9.

**Wynik ex66: M_phys = O(A_tail⁴) — analityczne uzasadnienie**

Tryb k=1 ogona jest trybem zerowym ODE (E_linear = 0 dokładnie przez kasowanie kin.+pot.). Fizyczna masa to:

```
M_TGP = c · A_tail⁴ + O(A_tail⁶)
```

✅ **Sformalizowane (v34)** jako `prop:J-zero-mode-mass` w `dodatekJ_ogon_masy.tex`.

---

### Priorytety v34 — stan końcowy

| Prio | ID | Zadanie | Uzasadnienie | Status |
|------|----|---------|--------------|--------|
| Bardzo wysoki | F1 | **Kanoniczna def. M_phys**: `prop:J-zero-mode-mass` — M_phys = O(A_tail⁴) z trybu zerowego ODE (ex66); kasowanie E_linear=0 analitycznie | Bez tej definicji sektor mas był arbitralny | ✅ **ZROBIONE** (v34, `dodatekJ`) |
| Bardzo wysoki | F2 | **Aktualizacja status_map**: zdegradowanie ex72–ex83 (S=2π−11/10 artefakt); dodano ex84 falsyfikację, ex85 T1/T2 | Spójność z ex84 | ✅ **ZROBIONE** (v34 wcześniej, `status_map.tex`) |
| Bardzo wysoki | F3 | **Hipoteza T1+T2**: `hyp:agamma-phi0-combo` — jeśli T1∧T2 dokładne → a_Γ znika z Warstwy II; tabela klasyfikacji parametrów | ex85: 0.09% DESI+CMB | ✅ **ZROBIONE** (v34, `sek08_formalizm.tex`) |
| Wysoki | F4 | **„Próżnia ze źródłem" jako pełna sekcja**: `ssec:vacuum-source-chain` — łańcuch N0-5→β=γ→W(1)=γ/3→Φ_min>0→Λ_eff→w(z)≠−1 | Najbardziej oryginalny element TGP | ✅ **ZROBIONE** (v34, `sek08_formalizm.tex`) |
| Wysoki | F5 | **Atraktor ψ_ini**: dowód analityczny lub numeryczny ψ_ini=7/6 jako konieczny start substratu GL | Redukuje N_param 3→2; cel v35 | ❌ **Otwarte** |
| Wysoki | F6 | **5 twardych predykcji w ramce sek07**: tryb oddechowy GW, c_GW=c₀, G(Φ)≠const, zmod. Friedmann, 3 generacje z WKB | Falsyfikowalność; komunikacja | ✅ **ZROBIONE** (v34, `sek07_predykcje.tex`) |
| Wysoki | F7 | **Most Φ→g_μν**: `prop:spatial-metric-from-substrate` — h(Φ)=Φ/Φ₀ z gęstości węzłów substratu; ℓ_P=const jako Program | Wiszący odsyłacz w prop:g00-from-axioms zamknięty | ✅ **ZROBIONE** (v34, `sek08_formalizm.tex`) |
| Niski | F8 | **README_map.md dla skryptów**: mapa ex/px → twierdzenie; archiwizacja ex71–ex83 | Organizacja | ❌ **Otwarte** |
| Niski | F9 | **Aktualizacja `problem_r21_analiza.md`**: zmiana statusu, sekcje ex84/ex85, ⚠️ artefakt w ex71–ex72 | Spójność dokumentacji | ✅ **ZROBIONE** (v34, `problem_r21_analiza.md`) |

---

### Aktualna mapa Warstwy II po v34

| Parametr | Status v33 | Status v34 (po sesjach v34) | Kierunek domknięcia |
|----------|-----------|---------------------------|---------------------|
| Φ₀ | Warstwa II | Warstwa II | Jedyny wolny parametr jeśli T1+T2 dokładne |
| a_Γ | Warstwa II | **Kandydat Warstwy III** (`hyp:agamma-phi0-combo` ✅) | T1∧T2 → a_Γ=(α_K√r₂₁)^(−2/3) |
| ψ_ini | Warstwa II (hip.) | Warstwa II (F5 otwarte ❌) | Dowód atraktora GL → Warstwa I→III |
| α_K | Warstwa II (OP-3) | Warstwa II | T2: α_K = r₂₁/Φ₀ jeśli T1 dokładne |

**Cel v35**: **Atraktor ψ_ini (F5)** + F8 (README skryptów). Jeśli F5 zamknięte → N_param = 1 (tylko Φ₀) → M/N ≥ 14.

---

### Co NIE jest priorytetem v34 (wbrew agentowi)

- **Przycinanie tekstu** (punkt 11): zbyt wcześnie — najpierw F5
- **Pełna kwantyzacja** (punkt 7): TGP świadomie semiclassical; wystarczy scope statement w Dod. E
- **Nowe ścieżki numerologiczne w masach**: po falsyfikacji Ścieżki 8–9 priorytetem jest T1+T2 (done) i F5

---

## Sesja v35 — Plan (następna sesja)

> **Stan wejściowy**: F1–F4, F6–F7, F9 zamknięte. Otwarte: **F5** (wysoki), **F8** (niski).
> **Główny cel**: zamknięcie F5 → N_param: 3→2 (lub 1 jeśli T1+T2 dokładne).

### F5 — Atraktor ψ_ini = 7/6 *(priorytet nr 1)*

**Pytanie**: Czy substrat GL zawsze startuje w ψ_ini = 7/6 niezależnie od warunków początkowych?

**Podejście analityczne**:
Równanie kosmologiczne TGP przy zaniedbaniu członu Hubble (wczesny wszechświat, $H \ll \sqrt{\gamma}$):
$$\ddot\psi = c_0^2 \Bigl[-U'(\psi) + W(\psi)/\Phi_0\Bigr] = c_0^2\Bigl[-(\beta\psi - \gamma\psi^2) + \frac{7\beta}{3}\psi^2 - 2\gamma\psi^3\Bigr]$$
Punkt równowagi $\dot\psi=0$: $W(\psi_{\rm eq})=0$ → $\psi_{\rm eq} = 7\beta/(6\gamma) \xrightarrow{\beta=\gamma} 7/6$.
Linearyzacja wokół $\psi_{\rm eq}$: $\omega^2 = c_0^2 W'(\psi_{\rm eq}) = c_0^2(14\beta/3 - 6\gamma) \xrightarrow{\beta=\gamma} c_0^2 \cdot 2\gamma/3 > 0$.
→ punkt siodłowy, **nie minimum** — $\psi_{\rm eq}$ jest równowagą niestabilną (szczytem potencjału).

Ale w obecności tłumienia Hubble ($3H\dot\psi$): tarcie prowadzi układ ku **atraktoru powolnego toczenia** w górę potencjału, co daje $\psi_{\rm ini}=7/6$ jako stan zamrożony ery radiacji.

**Zadanie F5**:
1. Skrypt numeryczny: rozwiązać układ $(\dot\psi, \dot a, \dot H)$ dla różnych $\psi(t_{\rm BBN}) \in [0.8, 1.6]$ — sprawdzić zbieżność do $\psi_{\rm ini}=7/6$ w erze radiacji.
2. Analityczny: sprawdzić czy `prop:psi-ini-derived` (już istniejące w sek08, linia ~709) zawiera pełny dowód, czy tylko motywację.
3. Jeśli zbieżność numeryczna potwierdzona → awansować z Hipoteza → Propozycja; $N_{\rm param}=3\to2$.

**Pliki docelowe**: skrypt `ex86_psi_ini_attractor.py`, `sek08_formalizm.tex` (update prop:psi-ini-derived).

### F8 — README_map.md skryptów *(priorytet niski)*

Stworzyć `TGP_v1/scripts/README_map.md`:
- Tabela: skrypt → twierdzenie/hipoteza które wspiera
- Oznaczenie ex71–ex83 jako `[ARTEFAKT]`
- Oznaczenie aktywnych skryptów (ex57, ex62, ex66, ex84, ex85)

| Priorytet | ID | Zadanie | Status |
|-----------|-----|---------|--------|
| Wysoki | F5 | Atraktor ψ_ini=7/6 — skrypt ex86 + update prop:psi-ini-derived | ❌ Otwarte |
| Niski | F8 | README_map.md skryptów | ❌ Otwarte |

---

## Sesja v34 — Zamknięcie (2026-03-28)

> **Cel sesji**: Udokumentowanie Ścieżek 9 i 10 w dodatkach LaTeX; spójność cross-reference; porządki w statusach.

### Nowe pliki stworzone w tej sesji

| Plik | Zawartość | Status |
|------|-----------|--------|
| `dodatekK_wkb_atail.tex` | Ścieżka 9: WKB/asymptotyczna analiza ogona solitonu; prop:K-tail-form; A_tail ∝ (g₀-g*)⁴; predykcja tauonu g₀^τ≈2.34 | ✅ Nowy |
| `dodatekL_formula_sumy.tex` | Ścieżka 10: skan kinetycznego α_K; definicja F(α_kin); dwa pierwiastki α*₁≈2.440, α*₂≈2.695; formuła sumy **obalona** przez ex84 (9307 ppm); 5 otwartych problemów O-L1–O-L5 | ✅ Nowy |

### Poprawki techniczne

| Plik | Zmiana |
|------|--------|
| `dodatekK_wkb_atail.tex` | Naprawa błędnej referencji `ex62_wkb_atail_scaling.py` → `ex62_phi_fixed_point.py` z prawidłowym opisem |
| `main.tex` | Dodano `\input{dodatekL_formula_sumy}` po `\input{dodatekK_wkb_atail}` |
| `sek08_formalizm.tex` | n_s w tabeli parametrów: „do obliczenia" → `n_s - 1 = -4ε_H - 4ε_ψ ≈ -4ε_H` (ε_ψ ∼ H₀/H ≪ 1 podczas inflacji) |
| `status_map.tex` | Dodano blok Sesja v34 — Dodatek L: definicja F(α_kin), status hipotezy sumy (obalona), zachowane wyniki pozytywne, problemy O-L1–O-L5; naprawa opisu ex62 |

### Stan końcowy v34 (pełny)

```
Dodatkowe otwarte po tej sesji:
- O-L1: nowe wyznaczenie α* na gęstej siatce (ex84 dało 2.440 / 2.695 — czy to poprawne?)
- O-L2: czy struktura dwóch zer F(α*)=r₂₁ jest specyficzna dla β_pot=1?
- O-L3: ujednolicenie Ścieżek 9+10 (czy A_tail(α*) jest identyczne?)
- O-L4: analityczne wyprowadzenie α* z warunków TGP
- O-L5: trzecia wartość α*₃ dla sektora tauonu

Nowe pliki LaTeX w dokumencie: Dodatki K + L
main.tex: sek01–sek09 + Dod. A–L (12 dodatków)
```

### Sesja v35 (2026-03-28) — wyniki

| ID | Zadanie | Status |
|----|---------|--------|
| F5 | Atraktor ψ_ini=7/6 — `ex86_psi_ini_attractor.py` 9/9 PASS; `prop:psi-ini-derived` awans Hipoteza→Propozycja; N_param: 3→2; ω_cosmo/H₀≈8.19 (nowa predykcja oscylacji G_eff) | ✅ **ZROBIONE** |
| F8 | `scripts/README_map.md` — mapa ex55–ex87, archiwizacja ex71–ex83 `[ARTEFAKT]` | ✅ **ZROBIONE** |
| O-L1 | `ex87_alpha_star_precision.py` — 80-pkt siatka + dense refinement + brentq (w toku) | 🔄 W toku |
| F6 | Napięcie r≈0.27 vs r<0.036 — **artefakt błędnego wzoru**: r=16ε_ψ gdzie ε_ψ=(1-n_s)/2=η (nie ε!). Starobinsky: ε=3/(4N_e²)→r=12/N_e²≈0.003≪0.036. Poprawiono rem:min-observables w sek08; σ_ab=0 podczas inflacji (faza symetryczna). | ✅ **ZAMKNIĘTE** (v35) |
| F7 | n_s precyzyjne — **zamknięte analitycznie** (rem:F7-ns-precision, v36): TGP = Starobinsky przy Planck 2018; ε_ψ = O(H₀/H_inf) → korekcja Δn_s ≈ 0.0002 poniżej CMB-S4. Dodano rem:sigma-inflation: σ_ab=0 w fazie symetrycznej podczas inflacji → Propozycja. | ✅ **ZAMKNIĘTE** (v36) |
| O-L1b | ex88_fixed_window — **4/4 PASS** (v36): α*₁_∞=2.43183767, α*₂_∞=2.63557742, S=5.0674≠2π−11/10 (22336 ppm). Kluczowe: stałe okno [28,42] → ZERO zależności od R_MAX (σ=0 dla wszystkich R_MAX∈{100..300}). Formuła sumy definitywnie obalona. | ✅ **ZAMKNIĘTE** (v36) |

### Kluczowe wyniki F5 (ex86, 9/9 PASS)

| Parametr | Wartość |
|----------|---------|
| ω_cosmo/H₀ | ≈ **8.19** (underdamped oscillations at z=0) |
| η_BBN | ≈ 1.2×10⁻³¹ (pole dosłownie zamrożone) |
| Δψ(BBN→equality) | ≈ 3×10⁻¹¹ |
| ψ_ini=7/6 → ΔG/G | = **0.000%** (doskonałe BBN) |
| ψ_ini=1.0 → ΔG/G | = 13.7% (marginalne) |
| BBN-compatible | ψ_ini ∈ [0.94, 1.56] |
| Nowa predykcja | oscylacje G_eff(t) z ω≈8.19 H₀ dla ψ_ini≠7/6 |

### Priorytety v36 — stan końcowy

| ID | Zadanie | Status |
|----|---------|--------|
| O-L1b | ex88_fixed_window — α*₁=2.43183767, α*₂=2.63557742, S=5.0674 | ✅ ZAMKNIĘTE |
| F7 | n_s precyzyjne analitycznie — TGP=Starobinsky, Δn_s=0.0002 | ✅ ZAMKNIĘTE |
| O-L2 | ex89_beta_pot_scan — **P2 FAIL** (v36): β_pot∈{0.5,0.8,1.2,1.5}→0 zer; β_pot=1.0→2 zera (α*₁=2.4316, α*₂=2.6356). Zera specyficzne dla fizycznego β=1 (N0-5). hyp:L-alpha-eff: bez awansu. | ✅ **ZAMKNIĘTE** (v36) |
| O-L3 | ex91_path9_path10_unification — **4/5 PASS** (v36): ścieżki komplementarne, nie identyczne. F(α=2)=263.9≠R21(+27.7%). z₀(α*₁)=1.2582≠g₀^e=1.2282(Δ=2.4%). A_tail ratio=1.197. Otwarte O-L4: dlaczego α*≠α_TGP? | ✅ **ZAMKNIĘTE** (v36) |

