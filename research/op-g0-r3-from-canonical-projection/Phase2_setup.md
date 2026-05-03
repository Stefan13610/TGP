---
title: "G.0 Phase 2 — Verification anchor predictions na nowej akcji S_TGP[V_M911, M9.1'']"
date: 2026-05-02
phase: 2
parent: "[[README.md]]"
predecessor: "[[Phase1_results.md]]"
status: ACTIVE — start 2026-05-03
score_gate: "≥3/4 PASS dla Phase 3 forward"
tags:
  - TGP
  - G0
  - phase2
  - V-M911
  - mass-spectrum
  - PPN
  - FRW
  - vacuum-stability
---

# G.0 Phase 2 — Verification anchor predictions

## Cel Phase 2

Po Phase 1 CLOSED-POSITIVE z odkryciem `V_M911(ψ) = -ψ²(4-3ψ)²/12`, Phase 2
weryfikuje, że **wszystkie kluczowe predykcje TGP są zachowane** na nowej
akcji `S_TGP[K=ψ⁴, V=V_M911, √(-g)=c₀·ψ/(4-3ψ)]`:

1. **Vacuum + spectrum** — ψ=1 stable, m_sp² physical
2. **Mass spectrum lepton** — m_μ/m_e=206.77, m_τ/m_e=3477 (PDG)
3. **PPN parameters** — γ=β=1 (Solar System tests)
4. **FRW cosmology** — κ=3/(4Φ₀), |dG/G|/H₀<0.02 (LLR)

Wynik PASS w ≥3/4 = **Phase 3 forward** (sek08a integration audit).

---

## Sub-tasks (4)

### P21 — Sympy LOCK V_M911 + vacuum spectrum

**Cel:** Formalna sympy weryfikacja, że V_M911 jest **unique** pod
kanonicznymi constraint'ami, plus pełna analiza linearyzacji wokół vacuum
ψ=1 (m_sp², stability, eigenmody).

**Metoda:**
1. Sympy proof uniqueness V_M911 z wymagań:
   - K(ψ) = ψ⁴ (T-D-uniqueness)
   - Static EOM = R3 ODE
   - U_eff(ψ) = ψV/(4-3ψ) ma vacuum w ψ=1
2. Linearization U_eff wokół ψ=1: `δU = ½·U''(1)·(δψ)² + O(δψ³)`
3. Mass eigenvalue: m_sp² = U''(1)/K(1) = γ
4. Stability: U''(1) > 0 ✓
5. Higher-order corrections (anharmonicity ψ³, ψ⁴ terms)
6. Comparison z hypothesis hyp:vacuum-mass z sek08a

**PASS criterion:** Sympy LOCK wszystkich 5 punktów + zgodność z sek08a
m_sp² = γ.

**Plik:** `phase2_P21_vacuum_uniqueness.py`

---

### P22 — Mass spectrum lepton verification

**Cel:** Weryfikacja że R3 mass formula
`m_obs(g₀, α=2) = c_M · A_tail²(g₀) · g₀^(e²/2)` produkuje
m_μ/m_e=206.77 i m_τ/m_e=3477 (Koide K=2/3) na nowej akcji V_M911.

**Metoda:**
1. Re-solve R3 ODE dla g₀ ∈ {g_e, g_μ, g_τ} z **identyczne EOM** jak w
   why_n3 (G.0 Phase 1 G0a potwierdziło numerical match `0.000000`)
2. Compute A_tail(g₀) = lim r→∞ r·(g(r) - 1)
3. Apply mass formula z α=2: m_obs ∝ A_tail² · g₀^(e²/2)
4. Compute ratios m_μ/m_e, m_τ/m_e
5. Compare z PDG (206.7682830, 3477.23)

**PASS criterion:**
- m_μ/m_e: |error| < 0.01% PDG (current why_n3: <0.001%)
- m_τ/m_e: |error| < 0.5% PDG (current why_n3: ~0.006%)

**Plik:** `phase2_P22_mass_spectrum_verification.py`

---

### P23 — PPN γ=β=1 re-derivation

**Cel:** Weryfikacja, że PPN parameters γ_PPN=β_PPN=1 są zachowane na
nowej akcji V_M911 + M9.1''.

**Metoda:**
1. Linearize M9.1'' metric wokół ψ=1: ψ = 1 + δψ
   - g_tt = -c²(1 - 4δψ + 3δψ² + O(δψ³))
   - g_rr = (1 + δψ - 3δψ² + O(δψ³))/(...)
2. Identify Newtonian potential U: -2U = first-order coefficient
3. Compute γ_metric (z g_rr expansion) i β_metric (z g_tt O(U²))
4. Apply "master formula z kinetic correction" (sek08c convention):
   - β_PPN = β_metric + 2c₂/f'(1) (gdzie c₂ z α=2 dynamics)
5. Verify: γ_PPN = 1, β_PPN = 1
6. Cross-check z V_M911 source term in linearized field eq.

**PASS criterion:**
- γ_PPN = 1 (Cassini: |γ-1| < 2.3·10⁻⁵)
- β_PPN = 1 (Mercury: |β-1| < 10⁻⁴)

**Plik:** `phase2_P23_PPN_verification.py`

---

### P24 — FRW cosmology + κ derivation

**Cel:** Weryfikacja, że FRW κ=3/(4Φ₀) (sek08a claim) re-derives z
**poprawnego** √(-g) = c₀·ψ/(4-3ψ) i V=V_M911.

**Metoda:**
1. FRW background: ψ(t), spatial uniform
2. Reduced FRW action z V_M911 + √(-g)=c₀·ψ/(4-3ψ):
   - L_FRW = a³·ψ/(4-3ψ) × [-K(ψ)·ψ̇²/(2c²)·(4-3ψ)/ψ - V_M911 - matter]
3. Variation δS_FRW/δψ = 0 → ψ-EOM in cosmology
4. Linearize wokół ψ=1: δψ ≈ small
5. Identify κ z source coupling
6. Check |dG/G|/H₀ z evolution
7. Compare with sek08a target κ=3/(4Φ₀)

**PASS criterion:**
- κ = 3/(4Φ₀) (re-derived analytically)
- |dG/G|/H₀ ≤ 0.02 (LLR Williams+ 2012)
- Optional: BBN |ΔG/G| < 0.15 (Cyburt+ 2015)

**Plik:** `phase2_P24_FRW_cosmology.py`

---

## Score gate Phase 2

```
Score = sum (P2x PASS), x ∈ {1, 2, 3, 4}

≥ 3/4 PASS  →  Phase 3 forward (integracja z sek08a + audit revision)
2/4 PASS    →  partial; analiza FAILed sub-tasks; ewentualnie pivot V_M911
< 2/4 PASS  →  G.0 hypothesis weakened; review fundamental
```

**Kill criterion:** Jeśli P22 (mass spectrum) FAIL → V_M911 nie reprodukuje
empirii lepton; G.0 closure questionable. Najprawdopodobniej PASS bo G0a
profile match był exact.

---

## Sub-task ordering

**Order:** P21 → P22 → P23 → P24

| Sub-task | Trudność | Czas | Wymaga |
|---|---|---|---|
| P21 | Średnia | ~2 dni | sympy LOCK + linearization |
| P22 | Niska | ~1 dzień | re-run R3 solver z mass formula |
| P23 | Wysoka | ~3 dni | PPN expansion + master formula |
| P24 | Bardzo wysoka | ~3 dni | FRW reduction + linearization |

**Total Phase 2: ~9 dni** (~2 tygodnie kalendarzowe).

---

## Hard anchors po G.0 closure

| Anchor | Wartość | Gdzie testowane |
|---|---|---|
| ψ_vacuum | 1 (stable) | P21 |
| m_sp² | γ | P21 |
| g₀_crit | 1.874 | P22 (z R3 solver, niezmienione) |
| m_μ/m_e | 206.7683 | P22 |
| m_τ/m_e | 3477.23 | P22 |
| γ_PPN | 1.000 | P23 |
| β_PPN | 1.000 | P23 |
| κ | 3/(4Φ₀) | P24 |
| |dG/G|/H₀ | < 0.02 | P24 |

---

## Plik scaffolding Phase 2

```
research/op-g0-r3-from-canonical-projection/
├── README.md                                  ✓ updated z Phase 1 status
├── Phase1_setup.md                            ✓
├── Phase1_results.md                          ✓ CLOSED-POSITIVE
├── phase1_G0a_volume_integration.py           ✓ PASS 4/4
├── phase1_G0b_field_redefinition.py           ✓ NEGATIVE-INFORMATIVE
├── phase1_G0c_einstein_frame_projection.py    ✓ PASS 3/3
├── Phase2_setup.md                            ✓ ten dokument
├── phase2_P21_vacuum_uniqueness.py            ← następne
├── phase2_P22_mass_spectrum_verification.py
├── phase2_P23_PPN_verification.py
├── phase2_P24_FRW_cosmology.py
└── Phase2_results.md
```

---

## Open questions z Phase 1 (resolution targets w Phase 2)

1. **Skąd fizycznie pochodzi V_M911?** — P21 zacznie odpowiadać przez
   identyfikację (4-3ψ)² jako "M9.1'' geometric weight"

2. **Czy mass formula α=2 + e²/2 zachowane?** — P22 weryfikacja

3. **Czy PPN master formula z sek08c (`β = 1/2 + 1/2`) działa z V_M911?**
   — P23

4. **Czy κ=3/(4Φ₀) jest niezmienne, czy zmienia się przez √(-g) update?**
   — P24

---

**Status:** Phase 2 ACTIVE — start 2026-05-03 (jutro nominalnie, ale
kontynuujemy w tej sesji jeśli czas pozwoli).

**Order startu:** P21 → P22 → P23 → P24 (od najprostszego do najtrudniejszego).
