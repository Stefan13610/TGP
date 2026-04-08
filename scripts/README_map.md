# Mapa skryptów TGP v1

> Ostatnia aktualizacja: sesja v45+ (2026-04-08). Baseline: v35 (2026-03-28); nowe: ex188–ex194.
>
> Konwencja statusu: ✅ Aktywny | ⚠️ Wymaga weryfikacji | ❌ Artefakt/superseded | 📦 Archiwum

---

## Seria ex — Numeryczne eksperymenty (nbody/examples/)

### Ścieżka 9: Ogon solitonu i hierarchia mas (ex55–ex61)

| Skrypt | Status | Wspiera | Sesja |
|--------|--------|---------|-------|
| `ex55_nodal_spectrum_tgp.py` | ✅ | spektrum nodalne solitonów | v33+ |
| `ex56_ghost_boundary_amplification.py` | ✅ | brak duchów w sektorze ψ⁴ | v33+ |
| `ex57_mass_amplitude_scan.py` | ✅ | `prop:J-scaling`: A_tail ∝ (g₀−g*)^4 | v33+ |
| `ex58_quantization_condition.py` | ✅ | `prop:J-quantization`: B_tail(g₀)=0 | v33+ |
| `ex59_mass_stability_scan.py` | ✅ | stabilność solitonów | v33+ |
| `ex60_golden_ratio_topology.py` | ✅ | `hyp:J-golden-topo`: g₀^μ/g₀^e≈φ | v33+ |
| `ex61_atail_substrate_coupling.py` | ✅ | `hyp:K-tau`: g₀^τ≈2.34, m_τ≈1777 MeV | v33+ |

### Ścieżka 10: Skan kinetycznego α_K (ex62–ex85)

| Skrypt | Status | Wspiera | Wynik |
|--------|--------|---------|-------|
| `ex62_phi_fixed_point.py` | ✅ AKTYWNY | def. F(α_kin)=(A(φz₀)/A(z₀))⁴; dwa zera α* | Dod. L |
| `ex62b_precise.py` | ✅ | precyzyjna wersja ex62 | — |
| `ex63_scattering_phase.py` | ✅ | faza scatteringowa δ(g₀*) | v33+ |
| `ex64_alpha_star_mraw.py` | ✅ | M_raw ∝ A^2.013 (ogon rozbiega) | v33+ |
| `ex65_mtotal_verification.py` | ✅ | kasowanie E_linear=0 | v33+ |
| `ex66_mcore_physical_mass.py` | ✅ | `prop:J-zero-mode-mass`: M_phys=O(A⁴) | v33+ |
| `ex67_zero_mode_mass.py` | ✅ | tryb zerowy ODE, E_linear=0 dokładnie | v33+ |
| `ex68_mren_subtraction.py` | 📦 | renormalizacja masy (wynik: ogon dominuje) | v33+ |
| `ex69_g0star_selection.py` | ✅ | δ(α)=g₀*−z₀ liniowe w α; δ=0 przy α≈2.44 | v33+ |
| `ex70_alpha_star_precise.py` | ⚠️ | wymaga re-weryfikacji po ex84 | v33+ |
| **ex71_delta_map.py** | ❌ **ARTEFAKT** | 25-pkt siatka → spurious zeros propagowane do ex83 | v33+ |
| **ex72_alpha_product.py** | ❌ **ARTEFAKT** | bazuje na błędnych α* z ex71 | v33+ |
| **ex73_sum_formula.py** | ❌ **ARTEFAKT** | S=2π−11/10 oparte na błędnych α* | v33+ |
| **ex74_beta_bifurcation.py** | ❌ **ARTEFAKT** | β=1 jako centrum — wymaga re-weryfikacji | v33+ |
| **ex74b_beta_fine.py** | ❌ **ARTEFAKT** | jw. | v33+ |
| **ex75_beta_critical.py** | ❌ **ARTEFAKT** | jw. | v33+ |
| **ex76_beta_precise.py** | ❌ **ARTEFAKT** | jw. | v33+ |
| **ex77_window_symmetry.py** | ❌ **ARTEFAKT** | jw. | v33+ |
| **ex78_alpha_window.py** | ❌ **ARTEFAKT** | jw. | v33+ |
| **ex79_fz0_minimum.py** | ❌ **ARTEFAKT** | jw. | v33+ |
| **ex80_algebraic_analysis.py** | ❌ **ARTEFAKT** | S=2π−11/10 bazuje na błędnych α*₁=2.441,α*₂=2.742 | v33+ |
| **ex81_sum_vs_R.py** | ❌ **ARTEFAKT** | jw. | v33+ |
| **ex82_dSdR.py** | ❌ **ARTEFAKT** | jw. | v33+ |
| **ex83_Rstar_precise.py** | ❌ **ARTEFAKT** | jw. | v33+ |
| `ex84_sum_formula_verification.py` | ✅ KLUCZOWY | **FALSYFIKACJA S=2π−11/10**; prawidłowe: α*₁≈2.440, α*₂≈2.695 (9307 ppm od formuły) | v34 |
| `ex85_agamma_phi0_precision.py` | ✅ KLUCZOWY | `hyp:agamma-phi0`: a_Γ·Φ₀≈1 (5/6 PASS, 0.12σ); T1+T2 kombinacja | v33+ |

> **⚠️ Uwaga dla ex71–ex83**: Wyniki z tej serii (α*₁=2.4414, α*₂=2.7418, S=2π−11/10)
> są **artefaktami** zbyt grubej siatki w ex71 (25 punktów).
> Prawidłowe wartości z ex84 (gęsta siatka + brentq): α*₁≈2.440, α*₂≈2.695.
> Formuła S=2π−11/10 jest obalona (9307 ppm odchylenia).

### F5: Atraktor ψ_ini = 7/6 (ex86)

| Skrypt | Status | Wspiera | Wynik |
|--------|--------|---------|-------|
| `ex86_psi_ini_attractor.py` | ✅ NOWY | `prop:psi-ini-derived` → awans Hipoteza→Propozycja | 9/9 PASS (v35) |

---

## Skrypty kosmologiczne (scripts/cosmology/)

| Skrypt | Status | Cel |
|--------|--------|-----|
| `cosmological_evolution.py` | ✅ | ewolucja ψ(z) od BBN do dziś (ψ_ini=1.0) |
| `bbn_attractor_resolution.py` | ✅ | BBN constraint vs Φ₀; ΔG/G |
| `friedmann_derivation.py` | ✅ | wyprowadzenie zmod. równań Friedmanna |
| `p73_perturbations_CMB.py` | ✅ | MS-TGP, n_s, r (prop:TGP-power-spectrum) |
| `bbn_consistency.py` | ✅ | τ₀·H_BBN ≫ 1 (prop:BBN-constraint) |
| `lambda_eff_estimation.py` | ✅ | Λ_eff ≈ γ/12 (prop:Lambda-eff) |
| `cosmological_chain.py` | 📦 | łańcuch N0→Λ_eff (starsza wersja) |
| `tgp_cosmo.py` | 📦 | ogólna biblioteka kosmologiczna |

## Skrypty grawitacyjne (scripts/gravity/)

| Skrypt | Status | Cel |
|--------|--------|-----|
| `einstein_emergence_proof.py` | ✅ | wyłonienie równań Einsteina z TGP |
| `covariant_field_equation.py` | ✅ | kowariantna postać D[Φ] |
| `two_body_three_regimes.py` | ✅ | trzy reżimy TGP (prop:trzy-rezimy) |

## Skrypty zaawansowane (scripts/advanced/)

| Skrypt | Status | Cel |
|--------|--------|-----|
| `p74_wetterich_flow.py` | ✅ | FRG dla sektora kinetycznego Z_k(φ)=φ⁴ |
| `falsification_map.py` | ✅ | mapa falsyfikowalności (ssec:falsification) |
| `tensor_sector_emergence.py` | ✅ | sektor tensorowy σ_ab |

## Skrypty stabilności (scripts/stability/)

| Skrypt | Status | Cel |
|--------|--------|-----|
| `gl_phase_transition.py` | ✅ | przejście GL substratu → ψ_ini=7/6 |
| `three_regimes_quantitative.py` | ✅ | trzy reżimy ilościowo |
| `vacuum_selection.py` | ✅ | selekcja próżni U'(1)=0 |

---

## Mapa ex → LaTeX (kluczowe powiązania)

| Skrypt | Twierdzenie/Hipoteza | Status |
|--------|---------------------|--------|
| ex57 | `prop:J-scaling`: A_tail ∝ g₀^4.12 | Propozycja |
| ex58 | `prop:J-quantization`: B_tail=0 | Propozycja |
| ex60 | `hyp:J-golden-topo`: g₀^μ/g₀^e≈φ | Hipoteza |
| ex61 | `hyp:K-tau`: m_τ≈1777 MeV | Hipoteza |
| ex62 | `app:formula-sumy` def. F(α_kin) | Propozycja |
| ex66 | `prop:J-zero-mode-mass`: M=O(A_tail⁴) | Propozycja |
| ex84 | obalenie S=2π−11/10 | Twierdzenie (neg.) |
| ex85 | `hyp:agamma-phi0`: a_Γ·Φ₀≈1 | Hipoteza (5/6 PASS) |
| ex86 | `prop:psi-ini-derived`: ψ_ini=7/6 | Propozycja (9/9 PASS) |
| p73 | `cor:TGP-power-spectrum`: P_s(k) | Propozycja |
| p74 | UV completeness FRG | Szkic |

---

## Priorytety v35+ (otwarte)

| ID | Skrypt do napisania | Twierdzenie docelowe |
|----|---------------------|---------------------|
| O-L1 | re-run ex84 gęsta siatka (>200 pkt) | potwierdzenie α*₁≈2.440, α*₂≈2.695 |
| O-L2 | ex87 (β_pot scan) | niezależność od β_pot |
| O-L3 | ex88 (unifikacja Ścieżek 9+10) | czy A_tail(α*) = A_tail(g₀*) |
