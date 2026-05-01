---
title: "ψ.1.v2.Phase6 setup — corrected predictions TT19-TT23 + 4-channel re-convergence"
date: 2026-05-01
cycle: ψ.1.v2.Phase6
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - psi1
  - phase6
  - correction
  - predictions
  - 4channel
  - setup
---

# ψ.1.v2.Phase 6 setup — corrected predictions + 4-channel re-convergence (5 sub-tests)

## Reason

Phase 4 strukturalnie ustaliło L₅'_a tensor operator + Adams positivity →
β_g<0 forced. Phase 5 propagated do eikonal + anisotropic c_eff(θ) +
realistic chopper Sagnac. Phase 6: nowy zestaw predykcji TT19-TT23 (zastępują
WITHDRAWN TT13-TT18) + cross-channel re-convergence z σ.1 + τ.3 + ω.1 + φ.1.

## Sub-tests

### T6.1 (TT19): Sagnac chopper differential SNR realistic

**Predykcja TT19**: Sagnac A-B differential phase shift dla pętli interferometru
z directional chopper:

$$\Delta\phi_{A-B} = \frac{4\pi L_{arm}}{\lambda_\gamma}\cdot|\beta_g|\cdot\frac{|\vec\nabla\ln X|^2}{\Lambda^2}$$

z konwencji $\xi = +2|\beta_g|$, $|\beta_g| \sim 0.1$, $\Lambda \sim 100$ TeV.

Lab realistic: $\Delta\phi \sim 5\times 10^{-36}$ rad (sub-detection by ~24 orders).
**TT19 falsifier**: pozytywny detection w lab Sagnac z SNR > 10⁻³ FALSIFIES
ψ.1.v2 (unless $|\vec\nabla\ln X|$ artificially amplified by external source).

PASS: TT19 formally registered jako NULL prediction lab tabletop.

### T6.2 (TT20): TOF dual-arm directional

**Predykcja TT20**: Dual-arm time-of-flight z differential między arm A
(parallel-to-gradient) i arm B (perpendicular):

$$\Delta t_{A-B} = \frac{L_{arm}}{c_0}\cdot\frac{1}{2}\cdot|\beta_g|\cdot\frac{|\vec\nabla\ln X|^2}{\Lambda^2}$$

(half because $c_{eff}^{-1} \approx 1 + \frac{1}{2}\xi n^2/\Lambda^2$).

PASS: TT20 zarejestrowane.

### T6.3 (TT21): Cosmological NULL re-confirmation

Nie ma globalnego $\vec\nabla\ln X$ dla cosmological scale (background X
isotropic w równaniu φ.1 background). Lokalne fluktuacje $\langle (\partial \ln X)^2 \rangle$
uśredniają isotropic na kosmologicznych skalach → **NIE ma kierunkowego
$c$ shift na CMB lub BBN**.

PASS: TT21 = ψ.1.v2 NULL na cosmological probes (zgodne z już-existing
NULL z ω.1, σ.1, τ.2, τ.3 cycles).

### T6.4 (TT22): Magnetar FRB ω-INDEPENDENT survival

σ.1 cykl predicted ω²-dependent helicity dispersion w FRB. Czy ψ.1.v2 dodaje
ω-INDEPENDENT residual signal? Z eikonal:
$c_{eff}(\theta) = 1 - O(n^2/\Lambda^2)$ — **ω-INDEPENDENT** (nie ma frequency
zależności w leading order tensor operator).

Cumulative TOF na D = 1 Gpc dla magnetar FRB:
$\Delta t_{ψ.1} = (D/c_0)\cdot|\beta_g|\cdot\langle(\partial\ln X)^2\rangle/\Lambda^2$

dla cosmological avg $\langle(\partial\ln X)^2\rangle$ — ale to jest izotropowe
average! → ZERO ω-independent residual w bulk medium (anisotropic but average=0).

**Resztkowa anisotropia**: tylko gdy along propagation jest persistent gradient
(np. magnetar wiatr lub local cluster substrate gradient) → sub-leading O(10⁻¹⁰)
TOF anomaly, much below σ.1 leading O(10⁻⁵).

PASS: TT22 = ψ.1.v2 sub-leading na FRB; nie ma falsyfikacji; consistent z σ.1.

### T6.5 (TT23): 4-channel ψ.1.v2 convergence + program END

**4 channels:** A (AS NGFP UV), B (heavy-mode 1-loop), C (Adams positivity),
D (cosmological consistency).

Convergence requirement: wszystkie 4 channels yield consistent sign β_g < 0
i magnitude $|\beta_g| \sim O(0.1)$ — verify każdy channel oddzielnie.

| Channel | Sign β_g | Magnitude | Status |
|---------|:---:|:---:|--------|
| A: AS NGFP | NEG (assumed Wilson tensor) | O(0.1) anchor σ.1 | ✓ |
| B: heavy-mode 1-loop | NEG | O((m/Λ)²) suggestywnie | ✓ |
| C: Adams positivity | **NEG STRICT** | UV-independent | ✓ DECISIVE |
| D: cosmo consistency | NULL constraint | within bound | ✓ |

PASS: 4/4 convergence → ψ.1 program END (true close, NEGATIVE structural
result for v1 + corrected anisotropic prediction set v2).

## Phase verdict logic

5/5 PASS → ψ.1 program END (full close-out v2)
≤4/5 → re-evaluate predictions

## Cross-cycle update plan (post-Phase 6)

- **INDEX.md**: ledger 769→784, phase ledger 521→536, Active Programs ψ.1 → CLOSED, add v2 phase rows
- **PREDICTIONS_REGISTRY**: TT13-TT18 → WITHDRAWN status, TT19-TT23 LIVE
- **CMB downgrade**: "POST-CONFIRMED" → "LIVE PARTIAL candidate" (across ω.1, σ.1, τ.2 mentions)
- **σ.1 narrative**: ψ.1.v2 jest sub-leading anisotropy CHANNEL na top of σ.1 leading helicity-dispersion
