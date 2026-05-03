---
title: "Hubble Tension w kontekście TGP"
date: 2026-05-03
tgp_status:
  folder_status: active
  level: L1
  kind: phenomenology
  core_compatibility: "unknown"
  last_reviewed_against_core: "unknown"
  may_edit_core: false
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges: []
  depends_on: []
  impacts: []
  source_of_status:
    - "README.md H1: 'Hubble Tension w kontekście TGP'"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: "2026-05-03"
---

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
- **DESI DR2 (arXiv:2503.14738, 2025-03)** — BAO+CMB+SN: H₀ = 68.17 ± 0.28
- **DESI Y3 (2026-04)** — 47 mln galaktyk, 2.8-4.2σ evolving DE

---

## Status post-cascade 2026-05-02

Wnioski tego folderu są **strukturalnie zachowane** po kaskadzie zamknięć
2026-04-26 → 2026-05-02 (closure_2026-04-26 + Phase 3.E + UV.3 + γ.1/δ.1/δ.2):

- `B_ψ/H_0² ~ 10⁻⁸` (gap 7.2 orders below required 0.17 dla H₀ tension)
  — **strukturalnie niezmienione** przez UV.3 (orthogonal: UV.3 = `Φ_0`
  absolute scale, hubble_tension = perturbacja backreaction)
- "TGP NIE rozwiązuje H₀ tension" — **wzmocnione** przez algebraiczne
  predykcje `Ω_Λ` z γ.1 + δ.1 + δ.2; cosmologia TGP ma teraz `Ω_Λ^pure = 2π/9`
  i `Ω_Λ^corr = 5e²/54` jako fixed prediction (bez room for H₀ shift)
- DESI DR2 (2025-03) actual `H₀ = 68.17 ± 0.28` — środek między Planck (67.36)
  i SH0ES (73.04), nadal ~5σ tension nierozwiązane
- TGP scope statement: GALAXY-SCALE GRAVITY + structural DE + CMB safety;
  **NOT** H_0 / S_8 tension solver

**Cross-link:**
- [[../op-cosmology-closure/M10_5_results.md]] — H₀/S₈ tensions audit (6/6 PASS, ct7 GREEN)
- [[../op-cosmology-closure/M10_R_results.md#125-post-m10-addenda-2026-04-28--2026-05-02]]
  — pełna updated falsification matrix
- [[../audyt_cosmology_drift_2026-05-03/README.md]] — audyt drift remediation 2026-05-03
