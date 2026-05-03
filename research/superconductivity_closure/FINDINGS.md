---
title: "FINDINGS — Program P5 — Superconductivity Closure"
date: 2026-05-03
parent: "[[README.md]]"
type: findings
tgp_owner: research/superconductivity_closure
source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)
tags:
  - findings
---

# FINDINGS — Program P5 — Superconductivity Closure

> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących
> plików folderu. **Żadne treści nie są wymyślone** — każdy item ma `source:`
> cytujący plik. Manual review w Sesji 6 (RESEARCH_BUS broadcasts).

## Numerical PASS counts (cited from .txt outputs)

- `ps32_results.txt` — 0/5 PASS — context: …======================================   Szlachetne (Cu,Ag,Au,Pt,Pd)            0/5 PASS, 0 BORDER, 5 FAIL  [ISSUE]   Alkaliczne (Li,Na,K)                   0/
- `ps32_results.txt` — 0/3 PASS — context: …0/5 PASS, 0 BORDER, 5 FAIL  [ISSUE]   Alkaliczne (Li,Na,K)                   0/3 PASS, 1 BORDER, 2 FAIL  [ISSUE]   Ferromagnetyki (Fe,Co,Ni)              0/3 P
- `ps32_results.txt` — 0/3 PASS — context: …0/3 PASS, 1 BORDER, 2 FAIL  [ISSUE]   Ferromagnetyki (Fe,Co,Ni)              0/3 PASS, 0 BORDER, 3 FAIL  [ISSUE]   Ziemie alk. (Ca,Sr,Ba,Mg)              0/4 P
- `ps32_results.txt` — 0/4 PASS — context: …0/3 PASS, 0 BORDER, 3 FAIL  [ISSUE]   Ziemie alk. (Ca,Sr,Ba,Mg)              0/4 PASS, 0 BORDER, 4 FAIL  [ISSUE]   Marginal SC (Zn,Cd)                    0/2 P
- `ps32_results.txt` — 0/2 PASS — context: …0/4 PASS, 0 BORDER, 4 FAIL  [ISSUE]   Marginal SC (Zn,Cd)                    0/2 PASS, 2 BORDER, 0 FAIL  [OK]   Insulatory (diamond,Si,Ge,NaCl,LiF)    0/5 PASS
- `ps32_results.txt` — 0/5 PASS — context: …0/2 PASS, 2 BORDER, 0 FAIL  [OK]   Insulatory (diamond,Si,Ge,NaCl,LiF)    0/5 PASS, 0 BORDER, 5 FAIL  [ISSUE]   Graphite                               0/1 PASS
- `ps32_results.txt` — 0/1 PASS — context: …0/5 PASS, 0 BORDER, 5 FAIL  [ISSUE]   Graphite                               0/1 PASS, 0 BORDER, 1 FAIL  [ISSUE]  =============================================

## TL;DR / Wynik / Verdict sections (cytaty)

### `README.md` — Kluczowe predykcje

> ## Kluczowe predykcje
> 
> 1. **Yb, Ce** (non-SC, score >85) — f-electron, kandydaci pod ciśnieniem
> 2. **Optymalne T_c** dla materiału z a ≈ 4.0-4.1 Å + d-orbitalami + FCC/BCC
> 3. **Anti-predykcja**: materiały z a/a_0 ≈ nπ (zera J(a)) nie będą SC
> 4. **Klatkowe struktury** (z=20) + d/f-orbitale = super-hydrydy
> 
> ---

## Eksportowalne formuły (boxed)

**`OSTRZAL_ROWNANIE_Tc_2026-04-20.md`:**

$$
T_c \;=\; \underbrace{\,C_0 \cdot A_{\mathrm{orb}}^2 \cdot k_d(z) \cdot M(a) \cdot \frac{\Lambda_E}{k_B}\,}_{T_c^{\mathrm{core}}} \;\cdot\; B_{\mathrm{orb}}(P)\,\cdot\,B_{\mathrm{mag}}(\lambda_{\mathr
$$

**`OSTRZAL_ROWNANIE_Tc_2026-04-20.md`:**

$$
\Lambda_E^{\mathrm{eff}}(\omega_{\mathrm{ph}}) \;=\; \Lambda_0 \,\left(\frac{\omega_{\mathrm{ph}}}{\omega_0}\right)^{\alpha_B},\qquad \alpha_B = 1.04
$$

**`OSTRZAL_ROWNANIE_Tc_2026-04-20.md`:**

$$
T_c \;\propto\; \omega_{\mathrm{ph}}^{\alpha_B} \;\propto\; M^{-\alpha_B/2} \;\approx\; M^{-0.52}
$$

**`OSTRZAL_ROWNANIE_Tc_2026-04-20.md`:**

$$
[T_c] \;=\; \underbrace{1}_{\mathrm{bezwym.}} \cdot \frac{[\mathrm{energia}]}{[\mathrm{energia/temperatura}]} \;=\; [\mathrm{temperatura}] \;\checkmark
$$

**`OSTRZAL_ROWNANIE_Tc_SERIA2_2026-04-20.md`:**

$$
r(|dlog|, \log T_c^{\mathrm{obs}}) = -0.13,\quad p = 0.54 \quad\text{(nieistotna)}
$$

---

## Cross-references

- [[README.md]] — opis folderu + YAML status
- [[NEEDS.md]] — otwarte luki tego folderu
- [[meta/research/RESEARCH_BUS.md]] — broadcast tych findings
- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa