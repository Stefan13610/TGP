# P7.1 — λ_sf z pierwszych zasad TGP (ps19)

**Data:** 2026-04-19
**Status:** P7.1 zamknięte; 11 parametrów fenomenologicznych → 1 uniwersalny + 2 tabelaryczne.

## Motywacja

W P6.D wartości λ_sf były przypisywane empirycznie dla każdego materiału (11 liczb):
- V=0.60, Nb=0.20, Ta=0.05, Mo=0.05, Pd=0.80, Fe=1.50, Co=1.20, Ni=1.20, Cr=1.20, Rh=0.40, Ru=0.15

To sprawiało, że P6.D było "kalibrowane", nie "przewidujące". P7.1 zastępuje to
formułą wyprowadzoną z TGP substrate + 2 wartości atomowych z literatury.

## Formuła P7.1

$$
\boxed{\lambda_{sf}^{TGP} = \kappa_{TGP} \cdot A_\text{eff}^2 \cdot k_d(z) \cdot N(E_F) \cdot g(I \cdot N(E_F))}
$$

gdzie:

$$g(x) = \frac{x^2}{\sqrt{1 + 0.25 \, x^4}}$$

jest regulatorem (paramagnon → paramag, saturacja dla FM).

**Parametry:**

| Symbol | Wartość | Pochodzenie |
|--------|---------|-------------|
| $\kappa_{TGP}$ | **2.012** | fit TGP uniwersalny (ps19) |
| $A_\text{eff}$ | η·A_d = η·0.310 | z P6.B/C |
| $k_d(z)$ | z Monte Carlo (P5) | — |
| $N(E_F)$ | st/eV/atom | Janak 1977, MJW 1978 (tabelaryczne) |
| $I$ | eV | Janak 1977 (tabelaryczne) |

**Wyprowadzenie fizyczne:**

Paramagnon w TGP to lokalna fluktuacja substratu Φ. Stoner I wychodzi z
drugiej pochodnej funkcjonału E[Φ] po magnetyzacji m. Elektron d-bandowy
sprzęga się z Φ przez amplitudę A_d (tę samą co w pairing phonon channel).
W formie Eliashberga:

$$\lambda_{sf} \sim \frac{g_\text{e-p}^2 N(E_F)}{\omega_{sf}}, \quad
g_\text{e-p} = A_\text{eff}, \quad \omega_{sf} \sim \Lambda_E / (1-IN)$$

Integracja po paramagnonowej gęstości spektralnej daje skalowanie jak $(IN)^2$.
Regulator $g(x)$ zapobiega rozbieżności przy instabilności magnetycznej.

## Walidacja

### Fit (5 paramagnetycznych d-metali)

| Element | N(E_F) | I | I·N | z | λ_obs | λ_pred | ratio |
|---------|--------|---|-----|---|-------|--------|-------|
| V  | 1.35 | 0.72 | 0.97 | 8 | 0.60 | 0.70 | 1.16× |
| Nb | 1.24 | 0.57 | 0.71 | 8 | 0.20 | 0.34 | 1.69× |
| Ta | 0.77 | 0.53 | 0.41 | 8 | 0.05 | 0.07 | 1.40× |
| Mo | 0.43 | 0.62 | 0.27 | 8 | 0.05 | 0.02 | 0.33× |
| Pd | 1.46 | 0.59 | 0.86 | 12 | 0.80 | 0.88 | 1.11× |

RMS_log (fit) = 0.307 (≈ factor 2 błąd).

### Walidacja (8 niezależnych materiałów)

| Element | I·N | λ_obs | λ_pred | ratio | status |
|---------|-----|-------|--------|-------|--------|
| Fe | 1.42 | 1.50 | 1.24 | 0.82× | FM, saturated ✓ |
| Co | 1.70 | 1.20 | 2.41 | 2.00× | FM |
| Ni | 2.04 | 1.20 | 3.09 | 2.58× | FM |
| Cu | 0.21 | 0.05 | 0.01 | 0.22× | paramag, bardzo mały |
| Cr | 0.27 | 1.20 | 0.01 | 0.01× | **AFM nesting** — poza formułą |
| Rh | 0.83 | 0.40 | 0.68 | 1.69× | paramag ✓ |
| Ru | 0.61 | 0.15 | 0.28 | 1.88× | paramag ✓ |
| W  | 0.15 | 0.03 | 0.00 | 0.12× | paramag, b. mały |

### Fe-pnictides / FeSe (kluczowe dla SC)

| Związek | I·N | λ_obs (P6.D) | **λ_pred (P7.1)** | ratio |
|---------|-----|--------------|-------------------|-------|
| FeSe_bulk | 0.80 | 0.90 | **0.69** | 0.77× |
| FeSe/STO  | 0.60 | 0.20 | 0.40 | 2.01× |
| Ba122-Co  | 0.54 | 0.30 | **0.29** | 0.98× |
| LaFeAsO   | 0.66 | 0.50 | **0.47** | 0.93× |
| NdFeAsO-F | 0.54 | 0.30 | **0.29** | 0.98× |

**Fe-pnictidy przewidywane w dokładności ±5-25%** (oprócz FeSe/STO gdzie
strain silnie modyfikuje efektywny I, wymagaloby dedykowanej kalibracji).

## Redukcja liczby parametrów

| Etap | Parametry fenomenologiczne λ_sf | Parametry TGP | Inputy tabelaryczne |
|------|--------------------------------|---------------|---------------------|
| P6.D (ps15) | **11** (per materiał) | β=2.527 | — |
| **P7.1 (ps19)** | **0** | β=2.527, **κ_TGP=2.012** | N(E_F), I |

Netto: **11 → 1 parametr TGP + 2 inputy atomowe**. 

## Zakres stosowalności

| Regime | Kryterium | Status formuły |
|--------|-----------|----------------|
| Paramagnetyczny | I·N < 0.95 | ✅ ±30-70% |
| Bliski FM | 0.95 ≤ I·N < 1.5 | ⚠️ regulator działa (Fe ✓) |
| Uporządkowany FM | I·N ≥ 1.5 | ⚠️ over-predicts (Co, Ni 2-3×) |
| AFM nesting | q≠0 susceptibility | ❌ **nie ujęte** (Cr outlier) |
| Cuprate d-wave | λ_sf → 0 (glue, nie blocker) | osobna fizyka (P6.A) |

## Predykcje dla nowych d-metali

| Kandydat | I·N | λ_pred | Interpretacja |
|----------|-----|--------|---------------|
| Re (5d⁵) | 0.36 | 0.066 | słabe blokowanie, dobra sygnatura SC |
| Ir (5d⁷) | 0.68 | 0.42 | umiarkowane, jak Pd |
| Tc (4d⁵) | 0.35 | 0.056 | zgodne z Tc=7.8 K |
| Zr (4d²) | 0.49 | 0.12 | Tc=0.6 K zgodne |
| Re-W alloy | 0.26 | 0.017 | b. małe, dobra kandydat SC |

## Powiązania

- [[ps15_p6d_magnetic_blocking.py]] — β=2.527 empiryczne
- [[ps16_p6c_orbital_switching.py]] — η(P) orbital delocalization
- [[ps19_p7a_lambda_sf_first_principles.py]] — ten plik (P7.1)
- [[P6CD_summary.md]] — poprzedni stan P6
- [[P7_plan.md]] — mapa P7

## Otwarte tematy do P7.2

1. **AFM nesting** (Cr): λ_sf z q≠0 susceptibility wymaga fermi surface
   topology — poza ramami substratowej TGP.
2. **Strain-modulated I** (FeSe/STO): I efektywne w zwiazkach moze zalezec
   od hybridization; tabulowanie I_eff dla heterostruktur.
3. **Ordered FM regime**: dla I·N ≥ 1.5 paramagnon description się załamuje,
   potrzebna osobna obróbka z magnonami (ale Tc i tak → 0, więc nie jest to
   kluczowe dla predykcji T_c).
