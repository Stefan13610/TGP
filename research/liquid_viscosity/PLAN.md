# Lepkość cieczy w stanach ekstremalnych (szkła, przechłodzone ciecze)

**Data startu:** 2026-04-20
**Status:** plan / open
**Kategoria:** TGP applications → miękka materia → transport pędu

## 1. Problem fizyczny

Dynamiczna lepkość $\eta$ jest jedną z najgorzej policzalnych wielkości
transportowych **z pierwszych zasad** (ab initio). Klasyczne relacje:

- **Stokes–Einstein**: $D\eta/(k_B T) = 1/(6\pi r)$ — łamana w
  przechłodzonych cieczach (Bhattacharyya 2008).
- **Vogel–Fulcher–Tammann (VFT)**: $\eta(T) = \eta_0\exp[B/(T-T_0)]$ —
  fenomenologiczny, $T_0$ nie ma fizycznego wyjaśnienia.
- **Mode-coupling theory (MCT)**: przewiduje idealne przejście szkliste
  w $T_c \gg T_g$ — błąd zwykle rzędu 20-50%.
- **Adam–Gibbs**: $\eta = \eta_0\exp[C/(T S_c(T))]$ — wymaga entropii konfiguracyjnej,
  trudno dostępnej eksperymentalnie.

**Uniwersalność Trachenko–Brazhkin** (2021): istnieje "fundamental lower bound"
$\eta_{\min}/\rho \approx \hbar/2m$ — skala kwantowa, niezależna od materiału.
Skąd ona się bierze?

## 2. Dlaczego TGP

W TGP **substrat ma sztywność** — pole $\Phi$ odczuwa wokół węzłów graf-substratu
silne gradientowe "naprężenia". W reżimie dobrze relaksowanym (wysokie $T$, mała
lepkość) węzły podsubstratowe poruszają się nieskoordynowanie. W granicy
$T\to T_g$ pojawia się **korelowana lokalna sztywność**:
$$\Phi(x)\to\Phi_\text{crit}\ \Rightarrow\ \text{local stiffness}\to\infty.$$

**Hipoteza:** lepkość dynamiczna w TGP wyraża się przez
$$\eta = \eta_0 \cdot \exp\!\left(\frac{\Delta[\Phi/\PhiZero]}{k_B T / E_\Phi}\right)$$
gdzie $\Delta[\Phi/\PhiZero]$ jest barierą konfiguracyjną substratu (wyprowadzalną
z $\beta, \gamma$ w potencjale), a $E_\Phi = \Lambda_E$ to naturalna skala energetyczna.

**Łączy to VFT, Adam–Gibbs i MCT w jedną formułę** z parametrami z rdzenia TGP,
jeśli hipoteza się potwierdzi.

Bonus: **dolny bound Trachenko** $\eta/\rho \sim \hbar/m$ powinien wychodzić z
zasady Heisenberga w kwantowaniu substratu — analogicznie jak "minimum thermal
conductivity" Cahill bierze się z kwantyzacji fononów.

## 3. Cele badawcze

### V1 — Lepkość jako funkcja $\Phi$-gradientu

Wyprowadzić z równania pola TGP (w stanie stacjonarnym z prędkością ścinania
$\dot\gamma$) efektywną naprężenie–odkształcenie i uzyskać formułę na
$\eta(T, P)$ z parametrów substratu.

### V2 — Przejście szkliste jako niestabilność substratu

Testować hipotezę: temperatura $T_g$ odpowiada warunkowi
$\langle\Phi\rangle = \Phi_\text{crit}$, gdzie $\Phi_\text{crit}$ jest rozwiązaniem
$\partial V/\partial\Phi = 0$ w metastabilnej gałęzi (mniejsza niż wartość
próżni).

### V3 — Fit VFT/MCT dla wzorcowych szkieł

Dopasować formułę TGP do danych dla:
- szkieł metalicznych: Vitreloy-1 (Zr41.2Ti13.8Cu12.5Ni10Be22.5), Pd43Cu27Ni10P20
- szkieł tlenkowych: SiO₂, B₂O₃, PMMA
- cieczy przechłodzonych: orto-terfenyl (OTP), glicerol, propylene carbonate
Oczekiwany wynik: 2-3 uniwersalne parametry TGP + tabelowane $E_\Phi^\text{mat}$.

### V4 — Kwantowy dolny bound

Wyprowadzić $(\eta/\rho)_\text{min} = f(\hbar, m, ...)$ z kwantyzacji
modów podstawowych substratu (per Trachenko-Brazhkin formula).

### V5 — Stokes–Einstein breakdown

Pokazać analitycznie kiedy $D \eta = \text{const} \cdot T$ się łamie, mapując
to na warunek $\nabla\Phi$-heterogeneity > scale $\ell_T$.

## 4. Plan numeryczny

- **ps01_viscosity_substrate_stress.py** — stress-strain z równania pola.
- **ps02_VFT_TGP_fit.py** — fit do 10 szkieł metalicznych/tlenkowych.
- **ps03_Tg_from_Phi_critical.py** — predykcja $T_g$ z $\rho, \Lambda_E$.
- **ps04_Trachenko_bound.py** — wyprowadzenie $\eta_\text{min}/\rho$.
- **ps05_SE_breakdown_map.py** — kryterium rozpadu Stokes-Einstein.

## 5. Literatura startowa

- Angell, *Formation of glasses from liquids and biopolymers* (1995) — fragility
- Trachenko-Brazhkin, *Minimal quantum viscosity from fundamental physical constants*,
  Sci. Adv. 6, eaba3747 (2020)
- Götze *Complex Dynamics of Glass-Forming Liquids* (2009) — MCT review
- Bhattacharyya-Bagchi, J. Chem. Phys. 119, 6667 (2008) — SE breakdown
- Biroli-Garrahan, J. Chem. Phys. 138, 12A301 (2013) — glass theory review

## 6. Relacje z innymi sektorami TGP

- **continuum_limit**: graf-substrat → ciecz ciągła → lepkość.
- **superconductivity_closure**: analogia "rigidity onset" między λ_sf a lepkość
  wysokotemperaturowa. Oba oparte na "substrate stiffening".
- **mass_scaling_k4**: bez k=4 kinetyki nie dostaniemy poprawnego skalowania
  kwantowego bound Trachenko.

## 7. Falsyfikowalność

- TGP przewiduje, że **wszystkie szkła** (metaliczne, tlenkowe, organiczne)
  mieszczą się w jednej formule z $\leq 3$ uniwersalnymi parametrami +
  jeden materiałowy. Obecnie VFT/Adam-Gibbs wymagają 3 materiałowych na każde.
- Konkretna predykcja: stosunek $T_g/T_m$ powinien być funkcją tylko
  $\Lambda_E^\text{mat}/E_\text{bind}$, nie klasy materiału.

## 8. Otwarte pytania

1. Czy $\Phi_\text{crit}$ glass-like jest TYM SAMYM obiektem co $\Phi_\text{crit}$
   z confinement-regime TGP (Three regimes — core paper)?
2. Jak włączyć wymiary - szkła 2D (Langmuir films) mają różne $T_g$?
3. Czy fragile-vs-strong (Angell) odpowiada różnym gałęziom $V(\Phi)$?

## 9. Link do rdzenia TGP

Core paper [[tgp_core.pdf]] § *Three regimes from non-monotonic
$\Phi$* — bezpośrednia baza (confinement regime = glass-like).

Potencjalnie sector mający największe połączenia z istniejącą fenomenologią
(VFT, Adam-Gibbs, Trachenko bound), więc dobra okazja do szerokiej walidacji
TGP na materiałach.
