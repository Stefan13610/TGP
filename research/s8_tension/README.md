# S₈ Tension w kontekście TGP

## Problem

Parametr S₈ = σ₈·√(Ωₘ/0.3) mierzy amplitudę fluktuacji gęstości materii
na skali 8 h⁻¹ Mpc. Dwa typy pomiarów dają systematycznie różne wartości:

| Metoda | S₈ | Źródło |
|--------|-----|--------|
| **CMB** (ekstrapolacja ΛCDM) | ~0.832 ± 0.013 | Planck 2018 |
| **Weak lensing** (późny Wszechświat) | ~0.76 ± 0.02 | KiDS-1000, DES Y3 |
| **Cluster counts** | ~0.77 ± 0.03 | eROSITA, SPT |
| **CMB lensing** | ~0.81 ± 0.02 | ACT, Planck lensing |

Rozbieżność: 2-3σ. Późne sondy konsekwentnie dają NIŻSZE S₈ niż
ekstrapolacja CMB. To sugeruje, że struktury rosną WOLNIEJ niż
przewiduje ΛCDM.

Świeży przegląd 2026 potwierdza, że problem nie zniknął —
jest traktowany jako drugie poważne napięcie (po H₀ tension).

## Interpretacja fizyczna

S₈ tension mówi: **struktury wielkoskalowe są gładsze niż powinny być**.
Albo:
1. Materia ciemna jest mniej „klumpowata" niż myśleliśmy
2. Grawitacja działa inaczej na dużych skalach
3. Coś tłumi wzrost struktur po rekombinacji

## Potencjalne połączenie z TGP

### Hipoteza 1: Tłumienie wzrostu przez substrat
TGP substrat z dynamiczną metryką g_ij może modyfikować
równanie wzrostu perturbacji δ̈ + 2Hδ̇ = 4πGρδ przez:
- Dodatkowe tarcie (viscosity substratu)
- Modyfikację efektywnego G na dużych skalach
- Korekty do propagatora grawitacyjnego

### Hipoteza 2: Skala odcięcia z continuum limit
Dyskretność substratu TGP wprowadza naturalną skalę odcięcia.
Jeśli ta skala wpływa na wzrost perturbacji, może tłumić
formowanie struktur na skalach ~8 Mpc.

### Hipoteza 3: Modyfikacja σ₈ przez solitonowe mody
Solitony TGP jako fundamental degrees of freedom mogą
modyfikować power spectrum na konkretnych skalach.

## Plan badawczy

1. **Przegląd stanu** (2026): zebrać najnowsze pomiary S₈
   (KiDS, DES, eROSITA, Euclid, DESI)
2. **Wzrost perturbacji w TGP**: jak substrat modyfikuje growth factor?
3. **Predykcja f·σ₈(z)**: TGP vs ΛCDM vs dane RSD
4. **Powiązanie z H₀ tension**: czy oba napięcia mają wspólne źródło w TGP?

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| `s8_tgp_analysis.py` | **10 pomiarów S₈, growth factor D(a) z TGP, porównanie z feedbackiem** | ⚠️ NEUTRALNY |

### Wyniki s8 (2026-04-18)

**Tension: 3.5σ** (Planck 0.832 vs lensing 0.771)
- Potrzebna supresja: 7.4%
- TGP dostarcza: ~0.001% (nu~1 przy y>>1 na skali 8 Mpc)
- Astrofizyka (AGN + neutrinos + systematics): 5-12% — **wystarczające**

**TGP nie rozwiązuje S₈, ale nie jest zagrożony:**
- S₈ prawdopodobnie tłumaczy się standardową astrofizyką
- TGP predykcja: ENVIRONMENT-DEPENDENT growth (void vs cluster)

## Powiązanie z ujednoliconym frameworkiem TGP

→ Patrz: `../cosmo_tensions/ct7_results.txt` — **definitywny werdykt**.

~~Mechanizm TGP dla S₈:~~ **NIEWYSTARCZAJĄCY**
- nu(y) ~ 1.0000 na skali 8 Mpc (y >> 1)
- Feedback baryonowy (3-5%) + neutrinos (1.5-3%) = wystarczające bez TGP

## Kluczowe referencje

- Heymans et al. 2021 — KiDS-1000
- DES Y3 Collaboration 2022
- Abdalla et al. — przegląd S₈ tension
- eROSITA cluster counts
- Euclid early results (2024/2025)
- Przegląd 2026 (do znalezienia)
