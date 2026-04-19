# Hubble Tension w kontekście TGP

## Problem

Dwa niezależne pomiary stałej Hubble'a H₀ dają systematycznie niespójne wyniki:

| Metoda | H₀ [km/s/Mpc] | Źródło |
|--------|----------------|--------|
| **Lokalne** (Cepheidy + SN Ia) | ~73.0–73.5 | SH0ES, Freedman et al. |
| **Wczesny Wszechświat** (CMB + ΛCDM) | ~67.4 ± 0.5 | Planck 2018/2020 |
| TRGB | ~69–70 | Freedman 2024/2025 |
| BAO + BBN | ~67–68 | DESI DR1/DR2 |

Rozbieżność: ~5σ (SH0ES vs Planck). To nie jest detal jednego eksperymentu —
to fundamentalny zgrzyt między **dwoma warstwami opisu kosmosu**:
fizyką wczesną (z=1100, CMB) vs późną (z<1, odległości lokalne).

## Potencjalne połączenie z TGP

### Hipoteza 1: Ewolucja metryki substratowej
TGP opisuje przestrzeń jako substrat z dynamiczną metryką g_ij.
Jeśli substrat ewoluuje kosmologicznie (np. g₀(z) zależy od redshiftu),
efektywna stała Hubble'a może różnić się między epokami:
```
H_eff(z) = H_ΛCDM(z) · (1 + ε(z))
```
gdzie ε(z) pochodzi z ewolucji substratu TGP.

### Hipoteza 2: Skalowanie stałej kosmologicznej
Jeśli Λ_eff wynika z energii zero-point substratu TGP,
to Λ może mieć słabą zależność od skali → modyfikacja H(z).

### Hipoteza 3: Dyskretność substratu
Continuum limit TGP może wprowadzać korekty skali Plancka
do propagacji fotonów na kosmologicznych dystansach → wpływ na distance ladder.

## Plan badawczy

1. **Przegląd stanu** (2026): zebrać najnowsze kompilacje H₀,
   w tym DESI DR2, SH0ES R22/R24, TRGB Freedman
2. **Model TGP dla H(z)**: jak ewolucja substratu modyfikuje Friedmanna?
3. **Predykcja**: czy TGP daje naturalną wartość H₀ między 67 a 73?
4. **Test**: porównanie z danymi BAO, SN Ia, CMB

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| `h0_tgp_analysis.py` | **Kompilacja 14 pomiarów H₀ + TGP backreaction budget + a₀(z) predykcja** | ❌ POZA ZASIĘGIEM |

### Wyniki h0 (2026-04-18)

**14 pomiarów H₀ skompilowanych:**
- Early Universe (weighted): H₀ = 67.5 ± 0.3
- Late Universe (weighted): H₀ = 72.1 ± 0.4
- Tension: ~8.5σ (ale JAGB daje 68.0 — konsystentne z Planck!)

**TGP backreaction budget:**
- Wszystkie mechanizmy (ct2-ct7): dH₀ ~ 0.07 km/s/Mpc
- Potrzebne (dla 73): dH₀ ~ 5.68 km/s/Mpc
- GAP: ~80x — STRUKTURALNY, nie parametryczny

**UNIKALNA PREDYKCJA TGP:**
- a₀(z) = c·H(z)/(2π) ewoluuje z redshiftem!
- Przy z=1: a₀ ~ 1.5 × a₀(z=0) = 1.8×10⁻¹⁰ m/s²
- TESTOWALNY z Euclid deep-field rotation curves
- ODRÓŻNIA TGP od MOND (gdzie a₀ = const)

## Powiązanie z ujednoliconym frameworkiem TGP

→ Patrz: `../cosmo_tensions/ct7_results.txt` — **definitywny werdykt**.

~~Mechanizm TGP dla H₀ tension:~~ **OBALONY (ct5/ct6/ct7)**
- ~~Amplifikacja tachyoniczna~~ → scale mismatch, B/H₀² = 10⁻⁹
- ~~Populacja solitonów~~ → d/λ_C ~ 10¹⁸, za rzadkie
- ~~Running gamma~~ → η=0.044, tylko ~1% zmiana

## Kluczowe referencje (do pobrania/przeczytania)

- Riess et al. (SH0ES) — najnowsze H₀ lokalne
- Planck 2020 — CMB H₀
- Freedman 2024/2025 — TRGB jako trzecia droga
- Di Valentino et al. — przegląd Hubble tension
- DESI DR2 (2025/2026) — BAO H₀
