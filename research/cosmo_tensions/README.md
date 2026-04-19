# Kosmologiczne napięcia jako jednorodny efekt substratu TGP

## Teza centralna

Trzy główne napięcia współczesnej kosmologii:
- **H₀ tension** (73 vs 67 km/s/Mpc)
- **S₈ tension** (0.76 vs 0.83)
- **DESI w(z) ≠ -1** (ewolucja ciemnej energii)

mają **wspólne źródło**: sprzężenie między lokalnym skupianiem materii
a globalną ekspansją substratu TGP.

### Mechanizm

W TGP przestrzeń **jest** substratem z metryką g_ij. Materia to solitony
(lokalne deformacje g ≠ 1). Ekspansja to globalna właściwość g = 1.

**Grawitacja dominuje lokalnie, przestrzeń dominuje globalnie.**

Gdy materia się skupia (formowanie struktur po rekombinacji):
1. Lokalnie g coraz bardziej odbiega od 1 (silniejsze solitony)
2. Globalnie pustki (voids) stają się czystszym g ≈ 1 → silniejsza ekspansja
3. Efektywne H_local > H_global (→ H₀ tension)
4. Substrat ma limit deformacji g₀_crit = 2.206 → nasycenie wzrostu (→ S₈)
5. <g> zmienia się w czasie → Λ_eff(z) ewoluuje (→ DESI w(z) ≠ -1)

## Struktura folderów

```
cosmo_tensions/          ← TEORIA (niniejszy folder)
├── README.md            ← ten plik
├── ct1_framework.tex    ← formalizm: zmodyfikowane Friedmanna z TGP
├── ct2_rho_tgp.py       ← obliczenia: ρ_TGP(a), w_eff(z)
├── ct3_hubble_fit.py    ← fit do H₀ data
├── ct4_s8_growth.py     ← growth factor z TGP tłumieniem
└── ct5_unified_fit.py   ← jeden model TGP → trzy napięcia

hubble_tension/          ← DANE H₀ (Cepheidy, TRGB, CMB, BAO)
s8_tension/              ← DANE S₈ (weak lensing, clusters, CMB lensing)
desi_dark_energy/        ← DANE DESI (BAO, w₀wₐ constraints)
```

## Kluczowe równania (rozwój w pliku ct1_framework.tex)

### 1. Energia substratu TGP

Gęstość energii substratu w stanie g(r):
$$\mathcal{E}[g] = \frac{1}{2}(\nabla g)^2 + V(g)$$

Stan próżni g = 1: $\mathcal{E}_{\rm vac} = V(1) \equiv \Lambda_{\rm TGP}$
Soliton: $\mathcal{E}_{\rm sol} = \mathcal{E}_{\rm vac} + \Delta\mathcal{E}$

### 2. Kosmologiczne uśrednienie

Średnia po objętości Hubble'a:
$$\langle g \rangle_H = \frac{1}{V_H} \int_{V_H} g(\mathbf{x}) \, d^3x$$

Rozpad: $g(\mathbf{x}) = 1 + \sum_i \delta g_i(\mathbf{x} - \mathbf{x}_i)$

gdzie $\delta g_i$ to solitony (cząstki) w pozycjach $\mathbf{x}_i$.

### 3. Zmodyfikowane równanie Friedmanna

$$H^2 = \frac{8\pi G}{3}\left(\rho_m + \rho_r + \rho_\Lambda + \rho_{\rm TGP}(a)\right)$$

gdzie $\rho_{\rm TGP}(a)$ to backreaction substratu:
$$\rho_{\rm TGP}(a) \propto f_{\rm clust}(a) \cdot \left(\frac{\delta_{\rm crit} - \langle\delta\rangle_{\rm sol}}{\delta_{\rm crit}}\right)$$

$f_{\rm clust}(a)$ = frakcja objętości w solitonach (rośnie w czasie).

### 4. Efektywne równanie stanu

$$w_{\rm TGP}(a) = -1 + \frac{a}{3\rho_{\rm TGP}} \frac{d\rho_{\rm TGP}}{da}$$

Ponieważ $f_{\rm clust}(a)$ rośnie → $\rho_{\rm TGP}$ zmienia się → $w \neq -1$.

## Hipoteza: TACHYONIC AMPLIFICATION — OBALONA (ct5/ct6)

Standard GR backreaction ~ (Φ/c²)² ~ 10⁻¹⁰ → za mały o 10 rzędów wielkości.

TGP ma trzy potencjalne mechanizmy amplifikacji:
1. Kinetic coupling K = ψ⁴ (czynnik ~4×)
2. Volume element √ψ (substrate budget)
3. ~~**TACHYONIC INSTABILITY V''(1) = -γ < 0**~~ ← OBALONA

### Dlaczego mechanizm tachyoniczny NIE DZIAŁA:

**Rozbieżność skali (scale mismatch)**:
- λ_tach = 2π/μ = 5632 Mpc (**super-Hubble!**)
- Formowanie struktur: 1-100 Mpc (sub-Hubble)
- Na skalach sub-Hubble: gradient przestrzenny ∇²ψ >> γψ → tachyon nieistotny
- Na skalach super-Hubble: brak istotnych perturbacji materii

| Skala | k²/(γ) | Tachyoniczny? |
|-------|---------|---------------|
| 8 Mpc | 12 553 | NIE (gradient 12553× silniejszy) |
| 100 Mpc | 80 | NIE |
| 500 Mpc | 3.2 | NIE |
| 1000 Mpc | 0.8 | TAK (ale brak źródła) |

**Wynik ilościowy (ct5)**: B_ψ/H₀² = 1.04×10⁻⁹ (wymagane: 0.174)
**Luka**: 8 rzędów wielkości — mechanizm jest ZANIEDBYWALNY.

ct3 dawał optymistyczną estymację 0.03-0.3 bo **pomijał gradienty przestrzenne**.
Gdy uwzględniamy pełną analizę k-space (ct5), efekt znika.

### Dodatkowa obserwacja: problem skali solitonów

Jeśli γ_micro = γ_cosmo, to solitony TGP mają rozmiar R ~ 900 Mpc (kosmologiczne!).
→ TGP **wymaga dwóch skal**: ℓ_micro (cząstki) i ℓ_cosmo (ekspansja).
To fundamentalna kwestia architektoniczna teorii.

### Otwarte kierunki

1. **Nie-perturbacyjne**: populacja solitonów zmienia efektywną ekspansję
2. **Running Λ**: γ biegnie ze skalą → Λ_eff(z) ≠ const (wymaga analizy RG)
3. **Osobne skale**: γ_micro ≠ γ_cosmo → efekty kolektywne
4. **Zmodyfikowana dyspersja**: struktura sieci TGP → wpływ na propagację fotonów
5. **Coupled dark sector**: masa solitonów zależy od ψ_background
6. **Inny K(ψ)**: K=ψ⁴ może nie być jedynym wyborem

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| `ct1_framework.tex` | Formalizm: backreaction z <3ψ̇²/ψ> | ✅ |
| `ct2_rho_tgp.py` | Model fenomenologiczny: B₀·f_struct(a), H(z), w(z), growth | ✅ |
| `ct3_dark_matter_backreaction.py` | Identyfikacja mechanizmów amplifikacji (tachyonic) | ⚠️ optymistyczny |
| `ct4_perturbation_growth.py` | Linearizowany wzrost perturbacji ψ, skala Jeansa | ✅ |
| `ct5_nonlinear_saturation.py` | Pełna analiza nasycenia: B_ψ/H₀² = 10⁻⁹ | ❌ **za mały** |
| `ct6_mechanism_diagnosis.py` | Diagnoza: dlaczego mechanizm zawodzi + alternatywy | 🔍 **KEY** |
| `ct7_soliton_cosmology.py` | **Ostatnie mechanizmy + HONEST VERDICT: TGP nie wyjaśnia napięć** | ❌ **FINAL** |

### Wyniki ct7: HONEST VERDICT (2026-04-18)

**Zbadane ostatnie 4 mechanizmy:**

| Mechanizm | B/H₀² | Potrzebne | Wynik |
|---|---|---|---|
| Populacja solitonów | ~0 | 0.174 | FAIL (cząstki za rzadkie, d/λ_C ~ 10¹⁸) |
| Running γ z RG | 2×10⁻³ | 0.174 | FAIL (η=0.044, only 1% zmiana) |
| Przejście fazowe substratu | 0 | 0.174 | FAIL (substrat zamrożony od z>>10¹⁰) |
| Architektura dwuskalowa | N/A | 0.174 | NO MECHANISM (γ_micro/γ_cosmo ~ 10⁸⁰) |

**WERDYKT:**
- **H₀ tension**: POZA ZASIĘGIEM TGP (8-9 rzędów za mały efekt)
- **S₈ tension**: NIEWYSTARCZAJĄCY (~2% zamiast ~8.5%)
- **DESI w(z)**: NIEKOMPATYBILNY (TGP daje w≥-1, DESI sugeruje phantom crossing)

**Konkluzja:** TGP = teoria GALAKTYCZNA. Ten sam mechanizm chameleonowy exp(-y⁰·⁸) który zapewnia bezpieczeństwo CMB i Układu Słonecznego eliminuje wpływ na kosmologiczne napięcia.
