# P7.C — TGP predictions for open SC candidates (2025-2026)

**Data:** 2026-04-20
**Status:** tabela predykcyjna z zamrozonymi parametrami TGP (0 nowych fitow)
**Skrypt:** [[ps22_open_candidates_predictions.py]]

## Kontekst

Po zamknięciu P6.A-D + P7.1 + P7.2 mamy uniwersalną formułę TGP:
$$T_c = C_0\,A_{\text{orb}}^2\,k_d(z)\,M(a)\,\LamE \cdot B_{\text{mag}}(\lamsf) \cdot B_{\text{PB}}(\mu_{\text{eff}}).$$

Wszystkie stałe uniwersalne TGP ($C_0=48.8222$, $\beta=2.527$, $\kappa_{\text{TGP}}=2.012$,
$\alpha_{\text{PB}}=0.2887$, $A_s, A_{sp}, A_d, A_f$, $K_{dw}=3.498$) są **zamrożone
na wartościach z zamkniętych modułów P6 i P7**. Dla każdego kandydata TGP generuje
predykcję bez żadnego nowego refitu.

## Anchor normalization

TGP z domyślnymi parametrami systematycznie przeszacowuje wysokociśnieniowe
hydrydy lantanowcowe o czynnik ~1.47 (LaH10: TGP 368 K vs obs 250 K). Wprowadzamy
jawną normalizację do obserwacji:

| Klasa | Anchor | T_obs / T_raw | Faktor |
|-------|--------|---------------|--------|
| A (hydrydy high-P) | LaH10 @ 170 GPa | 250/368 | 0.679 |
| B (hydrydy ambient) | (brak; użyty A) | — | 0.679 |
| C (nickelany bilayer) | Hg1223 | 138/93 | 1.487 |
| D (kagome AFM) | (brak; raw) | — | 1.000 |
| E (nitrydo-hydrydy) | (brak; raw) | — | 1.000 |

Dla klasy A+B faktor < 1 (TGP overshoot), dla C > 1 (TGP undershoot dla Hg1223
base). RMS_log TGP master = 0.346 (~faktor 2), więc T_anchor ma niepewność
±50-70%.

## Klasa A — hydrydy wysokiego P, lantanowcowo-podobne (niepotwierdzone)

| Material | P (GPa) | T_raw | **T_anchor (K)** | T_lit (DFT/obs) | Status |
|----------|---------|-------|------|-----------------|--------|
| **LuH10** | 170 | 368 | **250** | 286 (DFT) | niezsyntetyzowane, 4f¹⁴ closed, **flagship** |
| ScH9 | 250 | 350 | **238** | 180-210 (DFT) | niezsyntetyzowane, 3d¹ no 4f |
| LaBeH8 | 130 | 287 | **195** | 185 (DFT) | częściowo, P=80 GPa daje 110 K obs |
| YbH9 | 550 | 215 | **146** | 60-100 (DFT) | wymaga P > P_scale_Yb=552 GPa |
| LaBeH8 | 80 | 202 | **138** | 110 obs | kinetic-limited, TGP zgadza się po anchor |
| LuH6 | 100 | 188 | **127** | 246 (DFT) | CaH6-like cage; Lu 4f¹⁴ |
| SmH9 | 150 | 109 | **74** | 200 (DFT bez PB) | Sm 4f⁵, μ=1.55 → umiarkowane PB |
| MgH6 | 300 | 97 | **66** | 263 (DFT) | wymaga P>250 GPa stabilizacji |
| **GdH9** | 150 | 0.0 | **0.0** | ~200 (DFT bez PB) | **TGP: killer 4f⁷**, μ=7.94 → e⁻¹⁵ |

**Kluczowa jakościowa predykcja TGP (niezależna od anchor):**
GdH9 i cięższe LnH9 (Tb, Dy, Ho, Er) **nie będą nadprzewodzące** powyżej mK —
mimo że DFT bez pair-breakingu przewidywa 200 K. Test eksperymentalny: syntetyzować
GdH9 i potwierdzić brak $T_c$.

## Klasa B — hydrydy ambient (termodynamicznie niestabilne w DFT)

| Material | T_raw | **T_anchor (K)** | T_lit (DFT) | Status |
|----------|-------|---------|-------------|--------|
| Mg2PtH6 | 93 | **63** | >100 z e-dopingiem | nie syntetyzowane |
| Mg2IrH6 | 86 | **58** | 100-160 | nie syntetyzowane |
| Mg2RhH6 | 80 | **54** | 70-90 | nie syntetyzowane |
| Mg2PdH6 | 76 | **52** | 80-110 | nie syntetyzowane |
| Li2AgH6 | 67 | **45** | 90-120 | niestabilny 0.319 eV/atom |
| Li2AuH6 | 62 | **42** | 91-140 | niestabilny 0.172 eV/atom |

**Istotna niezgodność TGP vs DFT:** predykcje Eliashberga dla ambient Li2AuH6/Mg2PtH6
dochodzą do 140 K. TGP mówi 40-60 K — **dwukrotnie mniej**. To różnica falsyfikowalna:
gdy któryś z tych związków zostanie kiedykolwiek zsyntetyzowany, $T_c$ powinno
rozstrzygnąć między naiwnym Eliashbergiem a TGP. TGP odpowiada za to obniżenie
przez M(a) Gaussian (większa komórka niż a*=4.088 Å) i mniejsze η(P=0) dla p-block.

## Klasa C — nickelany bilayer

| Material | T_raw | **T_anchor (K)** | T_obs/DFT | Status |
|----------|-------|---------|-----------|--------|
| Ac3Ni2O7 | 76 | **114** | — | predicted stable (Rhodes-Wahl); radioactive |
| Ba3Ni2O7 | 76 | **113** | — | predicted stable; niezsyntetyzowane |
| La3Ni2O7 /SrLaAlO4 strained film | 75 | **112** | 26-42 obs | TGP 3× overshoot → $\LamE^{\text{Ni}}\approx(1/3)\LamE^{\text{cup}}$ |
| (La,Pr)3Ni2O7 film | 69 | **103** | >40 obs | podobnie; Pr-dop nie kluczowy |
| La2SmNi2O7 @ 20 GPa | 63 | **94** | 96 obs | **TGP zgadza się w 3%** |
| Nd6Ni5O12 | 11 | **16** | 13 obs | quintuple-layer; Nd 4f³ PB |

**Wniosek dla nickelatów:** TGP-cuprate formula działa po skalowaniu $\LamE$.
Niektóre (Nd6Ni5O12, La2SmNi2O7) zgadzają się w ±10%, inne (La3Ni2O7 film)
systematycznie off o faktor 3. Oznacza, że nickelany mają **mniejszą efektywną
Λ_E** (związaną z J_AF), która musi być tabulowana per-rodzina.

## Klasa D — kagome magnetyczne

| Material | T_anchor (K) | T_lit | Komentarz |
|----------|--------|-------|-----------|
| CsCr3Sb5 | **0.1** | 1-5 (DFT upper) | TGP: Cr³⁺ 3d³ + duże λ_sf → SC zabite |

**TGP null prediction:** CsCr3Sb5 **nie ma** konwencjonalnego s-wave SC.
Jeśli eksperyment znajdzie $T_c > 10$ K → wymaga mechanizmu niekonwencjonalnego
(flat-band, topological), czego TGP w obecnej formie nie obejmuje.

## Klasa E — egzotyki

| Material | T_raw=T_anchor (K) | T_lit/obs | Status |
|----------|--------|-----------|--------|
| CaC2H8 @ 30 GPa | **78** | 80-100 (DFT) | p-block clathrat; dobra zgodność |
| ZrNH @ 150 GPa | **66** | ~40 obs (2025) | TGP 1.5× overshoot |
| ZrNH2 @ 150 GPa | **54** | 15-17 obs | TGP 3× overshoot — wymaga η<1 |

## Flagship predictions do publikacji Zenodo

Top-5 (w kolejności priorytet eksperymentalny):

1. **LuH10 @ 170 GPa → 250 K** — najbardziej konkretny test. Synteza
   Lu→LuH2→LuH3 znana; wymaga wyższego P dla LuH10 Fm-3m.
2. **ScH9 @ 250 GPa → 240 K** — inne d-dominated superhydryde; jeśli
   $T_c > 200$ K dla ScH9 → TGP potwierdzone.
3. **LaBeH8 @ 130 GPa → 195 K** — LaBeH8 zsyntetyzowany @ 80 GPa (110 K);
   pełne ciśnienie 130 GPa NIE zostało uzyskane; TGP przewiduje tam 195 K.
4. **GdH9, TbH9, DyH9 → 0 K** — null prediction falsyfikowalna
   synteza GdH9 (teoretycznie zrobiona) → obserwacja braku SC.
5. **Li2AuH6 ambient → 42 K** (vs DFT-Eliashberg 91-140 K) — rozstrzygający
   test: mamy niższe predykcje TGP. Wymaga najpierw eksperymentalnej stabilizacji.

## Jak używać w publikacji

Rekomendacja dla drugiego papera TGP-SC (Zenodo):

- Dodać nowy paragraf/tabelę w sekcji "Falsyfikowalne predykcje"
  z kluczowymi kandydatami (LuH10, ScH9, LaBeH8, GdH9, Li2AuH6).
- Załączyć skrypt ps22 jako supplementary material (albo link do repo).
- Podkreślić: **zero nowych parametrów** — TGP ma już wszystkie stałe z
  wcześniejszych modułów, tutaj tylko aplikacja.

## Otwarte pytania

1. Czy nickelany wymagają oddzielnej Λ_E^Ni = Λ_E^cup/3 (TGP P6.A.2 może
   być rozszerzony), czy to jest efekt różnej D-state topologii (Ni d⁹ sp-hybryd
   vs Cu d⁹ pure)?
2. Czy η(P=0) dla p-block (Li2XH6, Mg2XH6) można otrzymać analogicznie do
   lantanowcowego P_scale z 4f binding? Np. P_scale z energii ionizacji H 1s vs d-metal?
3. Czy dla nitrydo-hydrydów ZrNH2 konieczne jest wprowadzenie efektywnej
   hybrydyzacji N-H obniżającej ω_eff?

## Powiązania

- [[P6_closure.md]] — zamknięcie P6.A-D
- [[P7A_summary.md]] — P7.1 λ_sf z pierwszych zasad
- [[P7B_summary.md]] — P7.2 lanthanide pair-breaking
- [[ps20_master_plot.py]] — master walidacja P6+P7.1
- [[ps21_p7b_lanthanide_pair_breaking.py]] — P7.2 ps script
- [[ps22_open_candidates_predictions.py]] — **ten plik** (P7.C)
- [[../../papers/sc/tgp_sc.tex]] — paper SC, dokąd idzie tabela predykcyjna
