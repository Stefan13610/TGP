# Anomalne przewodnictwo cieplne w kryształach molekularnych

**Data startu:** 2026-04-20
**Status:** plan / open
**Kategoria:** TGP applications → kondensowana materia → transport

## 1. Problem fizyczny

Standardowe modele transportu ciepła (Boltzmann + fonon dressed quasi-particles)
zawodzą w:

- **kryształach molekularnych** (naftalen, antracen, rubren): λ_th(T) ma plateau
  zamiast T⁻¹ w reżimie 50-300 K,
- **perowskitach hybrydowych** (MAPbI3, CsSnBr3): "ultraniskie" λ_th ≈ 0.3 W/m·K
  (ponieważ gęstość) i słaba zależność od masy kationu,
- **materii miękkiej** (soft crystals, plastic crystals): transport "glass-like" w krysztalicznej fazie.

Przy T = 300 K obliczone czasy relaksacji fononów z BTE są o rząd wielkości
krótsze niż obserwowane — więc transport nie jest fononowy w klasycznym sensie.
Pojawiają się hipotezy "fonon-roton soft modes", "diffuson-like channels",
"wavelike-particlelike duality" (Isaeva-Simoncelli 2019), żadna nie jest
uniwersalna.

## 2. Dlaczego TGP

W TGP pole $\Phi$ koduje lokalną gęstość substratu wygenerowaną przez masy:
$$\nabla^2\Phi + \frac{2(\nabla\Phi)^2}{\Phi} + \beta\frac{\Phi^2}{\PhiZero}
  - \gamma\frac{\Phi^3}{\PhiZero^2} = -q\PhiZero\rho.$$

Dla transportu ciepła kluczowe są:

- **Lokalne fluktuacje $\Phi$** odpowiadają "paczkom" substratu, które w SC-sektorze
  identyfikujemy z paramagnonami (λ_sf w P7.1). Tutaj mogą grać rolę nośników energii.
- **Metryka efektywna** $g_{ij}=e^{+2U/c_0^2}\delta_{ij}$: ciepło rozchodzi się
  po geodezjach substratu, nie po trywialnej sieci. Przy silnym gradiencie $\Phi$
  (molekularne "islands") propagacja fal staje się anizotropowa i podczas zderzeń
  z gradientem jest rozpraszana — nowy kanał rozpraszania, nieobecny w BTE.
- **N_0 / próg instabilności substratu**: w "miękkich" kryształach $\Phi$ może lokalnie
  spaść blisko wartości krytycznej, co efektywnie "wyłącza" niektóre mody.

## 3. Cele badawcze

### T1 — Dyfuzja energii w gradiencie $\Phi$

Wyprowadzić równanie dyfuzji ciepła w metryce $e^{+2U/c_0^2}\delta_{ij}$.
Oczekiwany wzór:
$$\partial_t u = \nabla\cdot\big(D(\Phi)\,\nabla u\big) + \xi(\Phi,\nabla\Phi)$$
gdzie $\xi$ to nowy czynnik rozpraszania TGP, zależny od $\nabla\Phi$ (źródło
zakrzywienia substratu).

### T2 — Plateau w kryształach molekularnych

Testować hipotezę, że plateau $\lambda_{\text{th}}(T)$ bierze się z faktu,
że powyżej pewnej $T$ typowa amplituda fluktuacji $\delta\Phi/\Phi$ staje się
niezależna od $T$ (saturacja), więc efektywna "phonon mean free path" przestaje
maleć z $T$.

### T3 — Skalowanie z gęstością / hybrydowych perowskitach

Wywieść skalowanie
$$\lambda_{\text{th}} \propto \rho^{p}\,M^{q}\,\LamE^{r}$$
i zmieścić znane dane (MAPbI3, CsPbBr3, Cs2AgBiBr6, PEDOT:PSS) w jednym wzorze.

### T4 — Uniwersalny Wiedemann-Franz z TGP

Sprawdzić, czy stosunek Lorentza $L = \kappa/\sigma T$ w metalach molekularnych
(e.g. Sr2RuO4, charge-transfer salts) otrzymuje poprawkę TGP modulującą klasyczne
$L_0 = \pi^2 k_B^2/3e^2$ przez czynnik związany z $\Lambda_E$ substratu.

## 4. Plan numeryczny (ps-skrypty)

- **ps01_thermal_substrate_diffusion.py** — numeryczne rozwiązanie równania
  dyfuzji ciepła w metryce TGP dla 1D kryształu; wyprowadzenie efektywnej $D(\Phi)$.
- **ps02_plateau_molecular_crystals.py** — fit $\lambda_{\text{th}}(T)$ dla
  naftalenu / antracenu / rubrenu z formułą TGP.
- **ps03_perovskite_hybrid_scan.py** — skaning 20+ perowskitów hybrydowych,
  predykcja $\lambda_{\text{th}}$ z $\rho, M, a$.
- **ps04_WF_ratio_TGP.py** — Lorenz number w solach przeniesienia ładunku.

## 5. Literatura startowa

- Simoncelli, Marzari, Mauri, *Unified theory of thermal transport* (2019)
- Isaeva et al., *Modeling heat transport in crystals and glasses* (2019)
- Ashcroft & Mermin rozdz. 25 (classic fononowy BTE) — baseline
- MAPbI3: Pisoni et al., J. Phys. Chem. Lett. 5, 2488 (2014), λ=0.3 W/m·K
- Rubren: Zeng et al., *Ultrafast thermal transport in organic semiconductors* (2020)

## 6. Relacje z innymi sektorami TGP

- **P6.A cuprates**: gęstość nośników z ZR singlet pool — tam też "phonon-less"
  kanał transportu (MIR band, Drude deficit) może być tym samym mechanizmem.
- **P7.1 λ_sf**: paramagnon-like fluktuacje $\Phi$ są wspólnym mianownikiem
  z naszym "phonon-like TGP energy carrier".
- **Continuum limit** (`research/continuum_limit/`): ciągłe granice substratu
  relacyjnego dla gradientu temperatury.

## 7. Falsyfikowalność

- Predykcja plateau $\lambda_{\text{th}}$ w anthracenie z $\beta, \gamma$ (TGP core)
  i jednego nowego parametru transportowego (np. $\eta_{\text{diff}}$).
- Jeśli dopasowanie do 4+ kryształów molekularnych wymaga różnego $\eta_{\text{diff}}$
  dla każdego → TGP transport-layer jest ad hoc (tak jak klasyczny BTE).

## 8. Otwarte pytania

1. Czy "gradient $\Phi$" rozpraszający ciepło jest tym samym obiektem co
   $\Bmag = 1/(1+\beta\lamsf)$ w P6.D, czy wymaga nowego operatora?
2. Jak parametryzować anizotropię molekularną (dipolowe momenty, π-staking)?
3. Czy pojawia się naturalny "minimum thermal conductivity" (analog Cahill) z
   granicy instabilności substratu ($\Phi \to \Phi_{\text{crit}}$)?

## 9. Link do rdzenia TGP

Core paper [[papers/core/tgp_core.pdf]] + Zenodo DOI.
SC paper [[papers/sc/tgp_sc.pdf]] jako wzór struktury aplikacyjnej.
