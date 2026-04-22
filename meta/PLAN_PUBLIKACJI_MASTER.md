# Plan publikacji TGP — master roadmap

**Data:** 2026-04-05
**Źródło:** Krytyka zewnętrzna (agent recenzent) + analiza wewnętrzna drzewa źródeł
**Cel:** Dwie ścieżki publikacji: (A) N-body paper (szybka), (B) Full TGP paper (docelowa)

---

## Ścieżka A: N-body paper (cel: 2-3 tygodnie)

Publikacja `tgp_nbody_results_clean.tex` jako samodzielny paper w CQG/PRD/EPJC.
Stan: **13 stron**, 8 Results, Introduction + Bibliography, 10/10 reproduce PASS.

| ID | Zadanie | Status | Wysiłek | Deliverable |
|----|---------|--------|---------|-------------|
| A1 | Dodać .bib (bibliografia) | **ZAMKNIETE** | done | `tgp_nbody.bib` (35 entries) |
| A2 | Tło literaturowe (modified gravity N-body, Efimov, chaos) | **ZAMKNIETE** | done | Sekcja Introduction + in-text citations |
| A3 | Paczka reprodukowalności (= H7) | **ZAMKNIETE** | done | `MANIFEST_PUBLICATION.md` + reproduce_all (10/10 PASS) |
| A4 | Proofreading + journal format | **ZAMKNIETE** | done | Submission-ready .tex (natbib, 0 errors) |
| A5 | Supplementary material (ex scripts + CSV) | **ZAMKNIETE** | done | `supplementary_material.tex` + manifest |

**Ścieżka A: GOTOWA do submission na arxiv + journal.**

---

## Ścieżka B: Full TGP paper — priorytety z recenzji zewnętrznej

### Faza I: Szybkie domknięcia (dni, nie tygodnie)

| ID | Zadanie | Status | Istniejący materiał | Brakuje | Wysiłek |
|----|---------|--------|---------------------|---------|---------|
| H7 | Paczka reprodukowalności | **ZAMKNIETE** | README.md, requirements.txt, reproduce_all.py (10/10 PASS, 101s) | -- | done |
| H4 | Tabela GW: polaryzacje/prędkości/ograniczenia | **ZAMKNIETE** | `tgp_gw_summary_table.tex` | -- | done |
| H2 | Podsumowanie mostu metrycznego | **ZAMKNIETE** | `tgp_metric_bridge_table.tex` | -- | done |

### Faza II: Konceptualne domknięcia (1-2 tygodnie)

| ID | Zadanie | Status | Istniejący materiał | Brakuje | Wysiłek |
|----|---------|--------|---------------------|---------|---------|
| H1 | Słownik czasu substratu + analiza EEP | **ZAMKNIETE** | `tgp_time_glossary.tex` | -- | done |
| H5 | Metrologia: bezwymiarowe obserwable + CODATA | **ZAMKNIETE** | `tgp_metrology_observables.tex` | -- | done |

### Faza III: Obliczenia analityczne (2-6 tygodni)

| ID | Zadanie | Status | Istniejący materiał | Brakuje | Wysiłek |
|----|---------|--------|---------------------|---------|---------|
| H3 | Pełny PPN (Will 2014 standard) | **ZAMKNIETE** | `tgp_ppn_full.tex`: 10/10 param = GR; metric theory → ξ=α₃=ζᵢ=0; ∂γ/∂φ=0 → α₁=α₂=0 | Caveat: strong-field (NS) left for future | done |
| M1 | Integracja N-body z main manuscript | **ZAMKNIETE** | `dodatekY_nbody.tex` → Appendix Y w main.tex | -- | done |
| M3 | Soczewkowanie w g_eff | **ZAMKNIETE** | `tgp_lensing_formal.py` (8/8 PASS), `tgp_lensing_geff.tex` | -- | done |

### Faza IV: Pipeline kosmologiczny (4-8 tygodni)

| ID | Zadanie | Status | Istniejący materiał | Brakuje | Wysiłek |
|----|---------|--------|---------------------|---------|---------|
| H6 | Formalne likelihood Planck+DESI | **ZAMKNIETE** | `tgp_formal_likelihood.py`, `tgp_cosmo_likelihood.tex` | MCMC opcjonalny (emcee) | done |
| M2 | Perturbacje kosmologiczne + LSS | **ZAMKNIETE** | `tgp_perturbations_formal.py` (18/18 PASS), `tgp_perturbations_lss.tex` | -- | done |

---

## Kolejność wykonania (rekomendowana)

```
Tydzień 1-2:   H7 (reproduce) → H4 (GW table) → H2 (metric summary)
Tydzień 2-3:   H1 (czas) → H5 (metrologia)
Tydzień 3-6:   H3 (PPN pełny) ← KRYTYCZNE
Tydzień 4-8:   A1-A5 (N-body paper submission)
Tydzień 6-14:  H6 (kosmologia formalna) + M2 (perturbacje)
```

---

## Ryzyka

| Ryzyko | Prawdopodobieństwo | Wpływ | Mitygacja |
|--------|-------------------|-------|-----------|
| PPN: α₁,α₂ ≠ 0 (preferred frame) | ŚREDNIE | KRYTYCZNY | Jawny rachunek w H3; jeśli ≠ 0, potrzebny screening |
| Screening: m_sp ~ √γ ~ H₀ → dalekozasięgowe | NISKIE | WYSOKI | c_GW = c₀ jest tożsamością; γ_PPN = 1 z metryki wykładniczej; sprawdzić Cassini bound |
| Kosmologia: TGP w_de ≈ -0.82 vs DESI w₀ ≈ -0.45 | **ZAMKNIETE** | ŚREDNI | H6 wynik: Δχ²=+1.8, ΔAIC=+3.8 — TGP ≈ ΛCDM, Φ₀→0 preferowane |
| N-body paper: brak literatury porównawczej | NISKIE | NISKI | Przegląd: Heggie & Hut 2003, Valtonen & Karttunen, chaos 3-body |

---

## Metryki sukcesu

- **Ścieżka A (N-body):** Paper submitted → arxiv + journal
- **Ścieżka B (Full TGP):** Wszystkie H1-H7 zamknięte → submission-ready manuscript
- **Reprodukowalność:** `python verify_all.py && python verify_nbody_lyapunov_quick.py` → 0 exit code
- **Kosmologia:** χ²(TGP) vs χ²(ΛCDM) z publicznymi danymi Planck+DESI
  - ✅ **DONE**: χ²(TGP)=42.1 vs χ²(ΛCDM)=40.3, Δχ²=+1.8, ΔAIC=+3.8
  - TGP zawiera ΛCDM jako gładki limit (Φ₀→0)
  - Dane DESI DR1 BAO + Planck 2018 CMB + CC H(z), N=47

---

## Status H1-H7 (2026-04-06)

Wszystkie zadania H1-H7 **ZAMKNIĘTE**:
- H1 (czas substratowy): `tgp_time_glossary.tex` ✅
- H2 (most metryczny): `tgp_metric_bridge_table.tex` ✅
- H3 (PPN pełny): `tgp_ppn_full.tex` — 10/10 param = GR ✅
- H4 (GW tabela): `tgp_gw_summary_table.tex` ✅
- H5 (metrologia): `tgp_metrology_observables.tex` ✅
- H6 (kosmologia): `tgp_formal_likelihood.py` + `tgp_cosmo_likelihood.tex` ✅
- H7 (reprodukowalność): `reproduce_all.py` (10/10 PASS) ✅

**Ścieżka B (Full TGP paper): GOTOWA do kompilacji manuskryptu.**

### Status M1-M3 (2026-04-06)

- M1 (N-body integration): `dodatekY_nbody.tex` → Appendix Y ✅
- M2 (perturbacje + LSS): `tgp_perturbations_formal.py` (18/18 PASS) + `tgp_perturbations_lss.tex` ✅
  - G_eff/G₀ - 1 ~ 10⁻⁵¹ → TGP ≡ ΛCDM w sektorze LSS
  - Stabilność: no-ghost (Q_s=ψ⁴>0), gradient (c_s²=c₀²>0), tachyonic (m_sp²=γ>0)
  - α_eff = 3.7×10⁻²⁶, bound obserwacyjny: α_eff < 0.17 (2σ)
  - n_s = 1 - 2/N_e ∈ Planck 1σ dla N_e ~ 57-65
- M3 (soczewkowanie): `tgp_lensing_formal.py` (8/8 PASS) + `tgp_lensing_geff.tex` ✅
  - α̂_TGP = 4GM/(c₀²b) = α̂_GR (γ_PPN = 1, dokładnie)
  - Shapiro delay = GR (Cassini bound spełniony: γ-1 = 0)
  - Korekta Yukawy: Δα̂/α̂ ~ α_eff² ~ 10⁻⁵¹ (niewykrywalna)
  - ppN korekta: O(U²) ~ 10⁻¹² (piko-arcsec)
  - Silne pole: cień BH mniejszy o ~37% vs Schwarzschild (otwarte: wymaga pełnego rozwiązania)

---

*Plan zakończony: wszystkie fazy H1–H7 i M1–M3 zamknięte per 2026-04-06. Weryfikacja: ex194 → 45/45 PASS (2026-04-08).*
