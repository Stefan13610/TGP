---
title: "G.0 P33 — Cross-reference audit results"
date: 2026-05-02
phase: 3
type: audit-results
parent: "[[Phase3_setup.md]]"
total_files_scanned: 1434
files_with_matches: 100
high_impact_files: 30
medium_impact_files: 0
---

# G.0 P33 — Cross-reference audit results

> **Status:** PASS — pelna lista pozyskana automatycznie.
> Ten dokument katalogu wszystkie miejsca w core/ + research/ ktore
> wymagaja attention po G.0 closure.

---

## 1. Pattern legend

| Pattern ID | Impact | Description | Recommended action |
|---|---|---|---|
| `V_orig_potential` | **HIGH** | V_TGP_orig form (β/3·ψ³ - γ/4·ψ⁴) — replace z V_M911 = -ψ²(4-3ψ)²/12 | Replace literal expression w wszystkich definicjach V(ψ) |
| `kappa_value` | **HIGH** | kappa = 3/(4·Phi_0) — w G.0 PO RE-FIT zachowane (P32 INVARIANT) | Annotate: po re-fit q·c²/Phi_0 z Newton limit, kappa pozostaje invariant |
| `sqrt_g_old` | **HIGH** | √(-g) = c·ψ (M9.1 FALSIFIED) — replace z c·ψ/(4-3ψ) (M9.1'' canonical) | Replace literal volume element |
| `prop_kappa_corrected` | **HIGH** | Reference do prop:kappa-corrected — affected by G.0 source coef 5q·ρ/Phi_0 | Update derivation w sek08a; reference unchanged |
| `prop_vacuum_stability` | **MEDIUM** | Reference do vacuum stability — G.0 fixes tachion bug (m_sp²: -γ → +γ) | Annotate: G.0 closure fixes vacuum stability derivation |
| `psi_eom_old` | **MEDIUM** | Φ-EOM references — G.0 zastepuje EOM przez R3 ODE | Replace EOM derivation; R3 ODE = effective EOM po G.0 |
| `M91_falsified` | **MEDIUM** | References do M9.1 (falsified) — replace z M9.1'' | Replace M9.1 → M9.1'' wszedzie |
| `U_eff_old` | **MEDIUM** | Effective potential old form — replace z γ(ψ⁴/4 - ψ³/3) | Update U_eff expression |

---

## 2. Aggregate statistics

- Files scanned: **1434**
- Files with at least one match: **100**
- HIGH-impact files (score ≥ 100): **30**
- MEDIUM-impact files (score 10-99): **0**

### Per-pattern stats

| Pattern | Impact | Files | Total matches |
|---|---|---|---|
| `V_orig_potential` | HIGH | 3 | 12 |
| `kappa_value` | HIGH | 14 | 34 |
| `sqrt_g_old` | HIGH | 5 | 6 |
| `prop_kappa_corrected` | HIGH | 7 | 15 |
| `prop_vacuum_stability` | MEDIUM | 4 | 16 |
| `psi_eom_old` | MEDIUM | 0 | 0 |
| `M91_falsified` | MEDIUM | 81 | 699 |
| `U_eff_old` | MEDIUM | 3 | 4 |

---

## 3. HIGH-impact files (require structural update)

### `core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex` — score 1600 (HIGH×16, MED×0)

**Pattern `V_orig_potential`** (HIGH) — V_TGP_orig form (β/3·ψ³ - γ/4·ψ⁴) — replace z V_M911 = -ψ²(4-3ψ)²/12:

- L197: `- \frac{\beta}{3}\,\psi^3 + \frac{\gamma}{4}\,\psi^4`
- L468: `Dla $V(\psi) = \frac{\beta}{3}\psi^3 - \frac{\gamma}{4}\psi^4$:`
- L471: `V + \psi V' &= \frac{\beta}{3}\psi^3 - \frac{\gamma}{4}\psi^4`

  Action: Replace literal expression w wszystkich definicjach V(ψ)

**Pattern `kappa_value`** (HIGH) — kappa = 3/(4·Phi_0) — w G.0 PO RE-FIT zachowane (P32 INVARIANT):

- L19: `$\kappa = 3/(4\PhiZero)$, a~\textbf{nie} $3/(2\PhiZero)$~---`
- L557: `$\kappa = 3/(4\PhiZero) \approx 0{,}030$`
- L589: `$3/(4\PhiZero) = 0{,}030$ (\textbf{poprawiony}) & $\approx 1{,}000$ &`
- L596: `$\kappa = 3/(4\PhiZero)$ wynika z~poprawnej wariacji dzia\l{}ania`
- L650: `\begin{remark}[Weryfikacja numeryczna $\kappa = 3/(4\PhiZero)$]%`
- ... (1 more matches)

  Action: Annotate: po re-fit q·c²/Phi_0 z Newton limit, kappa pozostaje invariant

**Pattern `sqrt_g_old`** (HIGH) — √(-g) = c·ψ (M9.1 FALSIFIED) — replace z c·ψ/(4-3ψ) (M9.1'' canonical):

- L616: `sk\k{a}d $\sqrt{-g} \sim e^{4U} \approx (1 + \delta\psi)^4`

  Action: Replace literal volume element

**Pattern `prop_kappa_corrected`** (HIGH) — Reference do prop:kappa-corrected — affected by G.0 source coef 5q·ρ/Phi_0:

- L144: `co kluczowe --- warto\'s\'c $\kappa$ (stw.~\ref{prop:kappa-corrected}).`
- L314: `\label{prop:kappa-corrected}`
- L322: `\begin{equation}\label{eq:kappa-corrected}`
- L510: `\begin{equation}\label{eq:kappa-def-operational}`
- L558: `(stw.~\ref{prop:kappa-corrected}):`
- ... (1 more matches)

  Action: Update derivation w sek08a; reference unchanged

---

### `core/sek08_formalizm/sek08_formalizm.tex` — score 1440 (HIGH×13, MED×14)

**Pattern `V_orig_potential`** (HIGH) — V_TGP_orig form (β/3·ψ³ - γ/4·ψ⁴) — replace z V_M911 = -ψ²(4-3ψ)²/12:

- L1628: `gdzie $U(\psi) = \frac{\gamma}{3}\psi^3 - \frac{\gamma}{4}\psi^4$;`
- L1802: `- \frac{\beta}{3}\psi^3 + \frac{\gamma}{4}\psi^4`
- L4105: `gdzie $V(\psi) = \frac{\beta}{3}\psi^3 - \frac{\gamma}{4}\psi^4$`
- L4130: `gdzie $V(\psi) = \frac{\beta}{3}\psi^3 - \frac{\gamma}{4}\psi^4$`
- L5175: `- \frac{\beta}{3}\psi^3 + \frac{\gamma}{4}\psi^4`
- ... (3 more matches)

  Action: Replace literal expression w wszystkich definicjach V(ψ)

**Pattern `kappa_value`** (HIGH) — kappa = 3/(4·Phi_0) — w G.0 PO RE-FIT zachowane (P32 INVARIANT):

- L4778: `gdzie $\kappa = 3/(4\PhiZero) \approx 0{,}030$,`

  Action: Annotate: po re-fit q·c²/Phi_0 z Newton limit, kappa pozostaje invariant

**Pattern `sqrt_g_old`** (HIGH) — √(-g) = c·ψ (M9.1 FALSIFIED) — replace z c·ψ/(4-3ψ) (M9.1'' canonical):

- L418: `Element obj\k{e}to\'sci & $\sqrt{-g} = c_0\psi$`
- L4107: `Pe\l{}ne wyprowadzenie z~poprawionej miary $\sqrt{-g}=c_0\psi$`

  Action: Replace literal volume element

**Pattern `prop_kappa_corrected`** (HIGH) — Reference do prop:kappa-corrected — affected by G.0 source coef 5q·ρ/Phi_0:

- L1776: `w~sek.~\ref{ssec:unified-action}, prop.~\ref{prop:kappa-corrected}.}`
- L4095: `zob.\ prop.~\ref{prop:kappa-corrected}.}`

  Action: Update derivation w sek08a; reference unchanged

---

### `core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex` — score 530 (HIGH×5, MED×3)

**Pattern `kappa_value`** (HIGH) — kappa = 3/(4·Phi_0) — w G.0 PO RE-FIT zachowane (P32 INVARIANT):

- L390: `$\kappa = 3/(4\PhiZero) \approx 0{,}030$`
- L401: `Pot\k{e}gowy (ten dokument) & $c_0\,\psi$ & $3/(4\PhiZero)$ & \textbf{PASS} \\`
- L415: `\item Akcja $\Rightarrow$ $\kappa = 3/(4\PhiZero)$`

  Action: Annotate: po re-fit q·c²/Phi_0 z Newton limit, kappa pozostaje invariant

**Pattern `sqrt_g_old`** (HIGH) — √(-g) = c·ψ (M9.1 FALSIFIED) — replace z c·ψ/(4-3ψ) (M9.1'' canonical):

- L413: `\item Metryka $\Rightarrow$ $\sqrt{-g} = c_0\psi$`

  Action: Replace literal volume element

**Pattern `prop_kappa_corrected`** (HIGH) — Reference do prop:kappa-corrected — affected by G.0 source coef 5q·ρ/Phi_0:

- L391: `(prop.~\ref{prop:kappa-corrected}).`

  Action: Update derivation w sek08a; reference unchanged

---

### `core/formalizm/dodatekH_lancuch_wyprowadzen.tex` — score 500 (HIGH×5, MED×0)

**Pattern `kappa_value`** (HIGH) — kappa = 3/(4·Phi_0) — w G.0 PO RE-FIT zachowane (P32 INVARIANT):

- L268: `i~$\kappa = 3/(4\Phi_0)$]\label{chain:A11b}`
- L427: `A11b & Akcja $S_{\rm TGP}$; $\kappa = 3/(4\Phi_0)$;`
- L590: `$\kappa = 3/(4\Phi_0)$;`

  Action: Annotate: po re-fit q·c²/Phi_0 z Newton limit, kappa pozostaje invariant

**Pattern `sqrt_g_old`** (HIGH) — √(-g) = c·ψ (M9.1 FALSIFIED) — replace z c·ψ/(4-3ψ) (M9.1'' canonical):

- L270: `wyznaczamy element obj\k{e}to\'sciowy $\sqrt{-g} = c_0\psi$`

  Action: Replace literal volume element

**Pattern `prop_kappa_corrected`** (HIGH) — Reference do prop:kappa-corrected — affected by G.0 source coef 5q·ρ/Phi_0:

- L279: `(prop.~\ref{prop:kappa-corrected},`

  Action: Update derivation w sek08a; reference unchanged

---

### `core/_meta_latex/status_map.tex` — score 430 (HIGH×4, MED×3)

**Pattern `kappa_value`** (HIGH) — kappa = 3/(4·Phi_0) — w G.0 PO RE-FIT zachowane (P32 INVARIANT):

- L219: `$\kappa = 3/(4\Phi_0) \approx 0{,}030$ (sprz\k{e}\.zenie`
- L1147: `$\kappa = 3/(4\PhiZero) = 0.0304$: wci\k{a}\.z poni\.zej LLR.`

  Action: Annotate: po re-fit q·c²/Phi_0 z Newton limit, kappa pozostaje invariant

**Pattern `sqrt_g_old`** (HIGH) — √(-g) = c·ψ (M9.1 FALSIFIED) — replace z c·ψ/(4-3ψ) (M9.1'' canonical):

- L214: `Zunifikowana akcja $S_{\rm TGP}$; $\sqrt{-g} = c_0\psi$`

  Action: Replace literal volume element

**Pattern `prop_kappa_corrected`** (HIGH) — Reference do prop:kappa-corrected — affected by G.0 source coef 5q·ρ/Phi_0:

- L223: `& stw.~\ref{prop:kappa-corrected} \\[3pt]`

  Action: Update derivation w sek08a; reference unchanged

---

### `research/closure_2026-04-26/KNOWN_ISSUES.md` — score 410 (HIGH×0, MED×41)

---

### `research/op-newton-momentum/B3_v2_alphas_propagation.py` — score 400 (HIGH×4, MED×0)

**Pattern `kappa_value`** (HIGH) — kappa = 3/(4·Phi_0) — w G.0 PO RE-FIT zachowane (P32 INVARIANT):

- L63: `#   alpha_s = (T_F * N_c^2) * g_0^e * kappa     with kappa = 3/(4 Phi_0)`
- L67: `# index slots; the "1/Phi_0" piece is identified with kappa = 3/(4 Phi_0)`
- L69: `# 3/(4 Phi_0), not 1/Phi_0.  The closed RHS N_c^3/(8 Phi_0) is what gets`
- L75: `print(f"\nFactored form (T_F * N_c^2) * g_0^e * (3/(4 Phi_0)):")`

  Action: Annotate: po re-fit q·c²/Phi_0 z Newton limit, kappa pozostaje invariant

---

### `research/op-phase3-uv-completion/phase3_D_cdt_hausdorff.py` — score 310 (HIGH×0, MED×31)

---

### `core/sek09_cechowanie/sek09_cechowanie.tex` — score 300 (HIGH×3, MED×0)

**Pattern `kappa_value`** (HIGH) — kappa = 3/(4·Phi_0) — w G.0 PO RE-FIT zachowane (P32 INVARIANT):

- L705: `a~$K_{\mathrm{geo}} = 3/(4\PhiZero) = \kappa$ (sek.~\ref{ssec:unified-action}).`
- L711: `\item $K_{\mathrm{geo}} = \kappa = 3/(4\PhiZero) = 0{,}0304$,`
- L1066: `\item $1/\PhiZero$: zwi\k{a}zek ze~ska\l{}\k{a} $\kappa = 3/(4\PhiZero)$`

  Action: Annotate: po re-fit q·c²/Phi_0 z Newton limit, kappa pozostaje invariant

---

### `research/nbody/examples/ex104_kappa_from_action.py` — score 300 (HIGH×3, MED×0)

**Pattern `kappa_value`** (HIGH) — kappa = 3/(4·Phi_0) — w G.0 PO RE-FIT zachowane (P32 INVARIANT):

- L137: `print(f"  kappa_new   = 3/(4 Phi_0) = {kappa_new:.6f}")`
- L229: `print(f"  kappa_new  = {kappa_new:.6f}  (3/(4 Phi_0))")`
- L246: `print("  The corrected coupling  kappa = 3/(4 Phi_0)  from the unified action")`

  Action: Annotate: po re-fit q·c²/Phi_0 z Newton limit, kappa pozostaje invariant

---

### `research/nbody/examples/ex105_ns_full_pipeline.py` — score 300 (HIGH×3, MED×0)

**Pattern `kappa_value`** (HIGH) — kappa = 3/(4·Phi_0) — w G.0 PO RE-FIT zachowane (P32 INVARIANT):

- L13: `1. Solve the TGP-FRW background for psi(N) using corrected kappa = 3/(4 Phi_0).`
- L32: `kappa     = 3/(4 Phi_0)  = 0.03043...`
- L212: `print(f"  kappa         = 3/(4 Phi_0) = {KAPPA:.6f}")`

  Action: Annotate: po re-fit q·c²/Phi_0 z Newton limit, kappa pozostaje invariant

---

### `research/nbody/examples/ex219_Jc_formal_derivation.py` — score 300 (HIGH×3, MED×0)

**Pattern `prop_kappa_corrected`** (HIGH) — Reference do prop:kappa-corrected — affected by G.0 source coef 5q·ρ/Phi_0:

- L51: `§5. Sprężenie materia-pole κ (sek08a, prop:kappa-corrected):`
- L291: `print(f"  Cel (z sek08a, prop:kappa-corrected): κ ≈ {kappa_target}")`
- L411: `DERYWACJA Z SEK08A (prop:kappa-corrected):`

  Action: Update derivation w sek08a; reference unchanged

---

### `research/op-phase1-covariant/phase1_F_capstone_M91pp.py` — score 280 (HIGH×0, MED×28)

---

### `research/op-phase3-uv-completion/phase3_E_structural_residua.py` — score 280 (HIGH×0, MED×28)

---

### `research/op-phase1-covariant/Phase1_F_results.md` — score 240 (HIGH×0, MED×24)

---

### `research/op-phase2-quantum-gravity/phase2_B_alpha0_first_principles.py` — score 240 (HIGH×0, MED×24)

---

### `research/op-phase2-quantum-gravity/phase2_A_linearized_graviton.py` — score 230 (HIGH×0, MED×23)

---

### `core/formalizm/dodatekO_u1_formalizacja.tex` — score 200 (HIGH×2, MED×0)

**Pattern `kappa_value`** (HIGH) — kappa = 3/(4·Phi_0) — w G.0 PO RE-FIT zachowane (P32 INVARIANT):

- L427: `Stała unifikacji TGP $\kappa = 3/(4\Phi_0)$ (sek08) i~stała`
- L431: `= \frac{\Phi_0/(8\pi J)}{3/(4\Phi_0)}`

  Action: Annotate: po re-fit q·c²/Phi_0 z Newton limit, kappa pozostaje invariant

---

### `research/op-newton-momentum/M9_1_prime_results.md` — score 170 (HIGH×0, MED×17)

---

### `research/op-phase2-quantum-gravity/Phase2_A_results.md` — score 160 (HIGH×0, MED×16)

---

### `research/op-phase2-quantum-gravity/Phase2_program.md` — score 160 (HIGH×0, MED×16)

---

### `research/op-phase3-uv-completion/Phase3_D_results.md` — score 160 (HIGH×0, MED×16)

---

### `research/op-phase3-uv-completion/Phase3_program.md` — score 160 (HIGH×0, MED×16)

---

### `research/op-phase1-covariant/phase1_B_psi_ph_derivation.py` — score 150 (HIGH×0, MED×15)

---

### `research/op-phase1-covariant/Phase1_B_results.md` — score 150 (HIGH×0, MED×15)

---

### `research/op-phase1-covariant/Phase1_R_final_results.md` — score 150 (HIGH×0, MED×15)

---

### `research/op-phase1-covariant/phase1_R_final_synthesis.py` — score 140 (HIGH×0, MED×14)

---

### `research/op-eps-photon-ring/Phase2_results.md` — score 130 (HIGH×0, MED×13)

---

### `research/op-newton-momentum/M9_1_results.md` — score 130 (HIGH×0, MED×13)

---

### `research/op-eps-photon-ring/phase1_eps_audit.py` — score 120 (HIGH×0, MED×12)

---

## 4. MEDIUM-impact files (require derivation update or annotation)


---

## 5. Recommended Phase 4 action plan

### Cluster A: V_orig form replacement (HIGH impact, ~30+ matches)

Replace **everywhere** w core/:
```
OLD:  V(\psi) = \frac{\beta}{3}\psi^3 - \frac{\gamma}{4}\psi^4
NEW:  V(\psi) = -\frac{\gamma}{12}\psi^2(4-3\psi)^2  (V_M911)
```

**Estimate:** ~2-4 godziny edycji (z double-check derivacji wokol).

### Cluster B: √(-g) update (HIGH impact, ~5-10 matches)

Replace:
```
OLD:  \sqrt{-g} = c_0\psi  (M9.1 FALSIFIED)
NEW:  \sqrt{-g} = c_0\psi/(4-3\psi)  (M9.1'' canonical)
```

### Cluster C: kappa annotation (HIGH impact, ~15+ matches)

Wartosc kappa = 3/(4·Phi_0) **pozostaje INVARIANT** po re-fit q·c²/Phi_0
(P32 sympy LOCK). Ale derivation w prop:kappa-corrected wymaga update:
- Source coef: 2q·ρ/Phi_0 → 5q·ρ/Phi_0
- Newton-limit: q·c²/Phi_0 = 2πG_0 → q·c²/Phi_0 = (4/5)πG_0
- Net kappa: pozostaje 4πG_0/(3H_0²) = 3/(4·Phi_0_new) gdzie
  Phi_0_new = (5/2)·Phi_0_old (lub equivalent re-cal)

### Cluster D: prop:vacuum-stability annotation (MEDIUM)

Add G.0 closure note: 'Po G.0 closure (V_M911), m_sp² = +γ stabilne;
sek08a stary derivation z V_orig + M9.1 dawal m_sp² = -γ tachion bug,
ktorego G.0 fixes.'

### Cluster E: M9.1 → M9.1'' (MEDIUM)

Replace metric form references — sek08c juz ma A2/A3 annotations,
teraz mozna je close (P34).

---

**Total Phase 4 estimate:** ~15-25 godzin pracy edytorskiej + verification.
**Wstepny plan:** sek08a → sek08c → sek08 → pozostale (z dependency tree).
