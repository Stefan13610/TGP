---
title: "Phase 1 setup — α_PB ↔ α_0 unit-bridge audit"
date: 2026-04-28
cycle: SC.1.Phase1
status: PRE-EXECUTION
predecessor: "[[program.md]]"
tags:
  - TGP
  - SC
  - alpha-PB
  - dimensional-analysis
  - sympy
---

# Phase 1 — Setup: unit-bridge audit α_PB ↔ α_0

> **Cel:** Sprawdzić jednoznacznie czy α_PB ≈ 0.2887 μ_B⁻² (TGP SC v2) i α_0 ≈ 4.04 (T-α threshold) to **ta sama** stała w różnych systemach jednostek.

---

## H₀ vs H₁ recap

| Hipoteza | Twierdzenie | Falsyfikacja |
|---------|-------------|--------------|
| **H₀** (unit cousin) | α_PB = f · α_0 dla pewnej f z czystym czynnikiem jednostek | różne wymiary po pełnej redukcji do bazowych SI |
| **H₁** (strukturalnie różne) | α_PB i α_0 są niezależnymi stałymi z różnej fizyki | jakakolwiek reprezentacja `f(α_0)` powiela strukturalnie pasującą wymiarową α_PB |

---

## Procedura testów

### T1.1 — Pełna analiza wymiarowa (sympy `dimsys_SI`)

**Cel:** Wyprowadzić wymiar SI każdej α; sprawdzić czy są równe lub czy jakaś naturalna konwersja je łączy.

**Plan:**
1. `α_PB` ma wymiar `[μ_B⁻²]`. W SI: `μ_B = 9.274×10⁻²⁴ J/T`, więc `μ_B² = J²/T² = (kg·m²/s²)²/(kg/(A·s²))² = A²·s⁶·m⁴ / s⁴ · (1/kg²) = A²·m⁴·s²/kg²` (po uproszczeniu). Krótko: μ_B² = J²/T².
2. `α_0` ma wymiar `[1]` (dimensionless) — z T-α: α(ψ) = α_0 · (ψ-1)² Θ, gdzie ψ jest bezwymiarowe (substrate field normalized).
3. Aby `α_PB = f · α_0`, czynnik f musiałby mieć wymiar `[μ_B⁻²]` = `T²/J²`. Ale w teorii TGP nie ma fundamentalnej skali magnetycznej (nie ma ani μ_B w core ani β/κ_TGP nie ma wymiaru T). Więc f nie może być **pure structural constant** TGP.

**Wniosek oczekiwany:** **H₀ rejected na poziomie wymiarowym**. α_PB i α_0 nie są unit cousin.

### T1.2 — Geometric-unit equivalence check

**Cel:** Sprawdzić czy w geometric units (G = c = ℏ = k_B = 1) jakakolwiek naturalna kombinacja stałych daje μ_B⁻² jako "structural" output.

W geom units:
- μ_B = e·ℏ/(2m_e) [w SI]
- W geom: μ_B → α_em / (2 m_e^geom), gdzie α_em jest fine-structure, m_e^geom = m_e/M_Pl
- α_PB · μ_B² → α_PB · α_em² / (4 m_e²) [w geom]

Jeśli α_0 ↔ α_PB byłoby "ten sam" w geom units, oznaczałoby `α_0 · 1 = α_PB · μ_B² → α_PB = α_0 · 4 m_e² / α_em² ≈ 4.04 · 4 · (5.4×10⁻⁴⁴)² / (1/137)² ≈ ... numerycznie ekstremalnie małe`.

**Plan:** policzyć eksplicytnie i sprawdzić czy result ≈ 0.2887 (niemożliwe per szacunkowi), zalogować mismatch.

### T1.3 — Czy T-α residue ma reprezentację μ_B⁻²?

T-α: α(ψ_ph) = α_0 · 0.0282 ≈ 0.114 — bezwymiarowy.
α_PB · ⟨μ_eff²⟩ z LnH₉ (μ_eff ~ 3.6 μ_B): 0.2887 · 12.96 = 3.74. **Też O(1)!**

To ciekawe: α_PB · ⟨μ_eff²⟩ ≈ 3.74 jest blisko α_0 ≈ 4.04 numerycznie. Zbieżność czy zbieg okoliczności?

**T1.3a:** Sprawdzić eksplicytnie czy istnieje calibration choice μ_eff = `c · μ_B` (gdzie c jakiś TGP scaling) pozwalający wyprowadzić α_PB z α_0. Jeśli c jest **niezależnie znane** z TGP core (np. relacja do α_em, m_e/M_Pl), to mamy bridge.

**T1.3b:** Jeśli c wymaga *fitu* na ⟨μ_eff⟩ z LnH₉ (czyli c jest input, nie output), to **nie ma** strukturalnego bridge — tylko numeryczne zbieżność.

### T1.4 — Numerical coincidence vs structural identity

**Cel:** Distinguish mere numerical coincidence (within ~10% bracket) from structural identity (exact relation derivable z TGP core constants).

**Test:** sprawdzić czy ratio `α_PB · ⟨μ_eff²⟩ / α_0 = 3.74 / 4.04 = 0.926` można wyrazić jako **exact closed-form** w TGP core constants (κ_TGP, β, α_em, ...). Jeśli tak — coincidence jest strukturalna. Jeśli pasuje tylko numerycznie z błędem ~10% — to numerical coincidence (tj. brak bridge).

**Threshold akceptacji structural identity:** odchyłka < 0.1% i closed-form representation in ≤ 3 core constants.

---

## Materiał wykonawczy

- **Skrypt:** `phase1_unit_bridge_audit.py` (sympy dimensional analysis + numerical checks)
- **Output:** `phase1_unit_bridge_audit.txt` (full execution log)
- **Audit memo:** `Phase1_results.md` (post-execution conclusion)

## Środowisko

Python 3.x + sympy + numpy. Brak zewnętrznych zależności.

```bash
cd TGP/TGP_v1/research/op-sc-alpha-origin
PYTHONIOENCODING=utf-8 python -X utf8 phase1_unit_bridge_audit.py 2>&1 | tee phase1_unit_bridge_audit.txt
```

## Co znajdzie się w wynikach

- **PASS** dla każdego sub-testu T1.1–T1.4 jeśli H₁ jest poparte (oczekiwane).
- **FAIL** dla T1.1 oznaczałoby naprawdę nietrywialny strukturalny bridge — zaskoczenie, do raportu.

---

## Cross-references

- `program.md` — overall 3-phase plan
- `../closure_2026-04-26/alpha_psi_threshold/results.md` — α_0 origin
- `../../../tgp-sc-paper/paper/tgp_sc.tex` (lines 425-466) — α_PB definition
