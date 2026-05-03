---
title: "FINDINGS — op-quantum-closure"
date: 2026-05-03
parent: "[[README.md]]"
type: findings
tgp_owner: research/op-quantum-closure
source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)
tags:
  - findings
---

# FINDINGS — op-quantum-closure

> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących
> plików folderu. **Żadne treści nie są wymyślone** — każdy item ma `source:`
> cytujący plik. Manual review w Sesji 6 (RESEARCH_BUS broadcasts).

## Numerical PASS counts (cited from .txt outputs)

- `m11_1_audit.txt` — 6/6 PASS — context: …] M11.1.4: PASS     [v] M11.1.5: PASS     [v] M11.1.6: PASS    >> M11.1 CLOSED (6/6 PASS) — 1-loop V_eff audit complete:   >> structural β=γ preservation (M2b 
- `m11_2_betafn.txt` — 6/6 PASS — context: …=====================================================================   Result: 6/6 PASS     [PASS] M11.2.1 LPA FP reproducibility     [PASS] M11.2.2 LPA' η ex
- `m11_2_betafn.txt` — 6/6 PASS — context: …2.5 WF eigenspectrum     [PASS] M11.2.6 Branch I/II η reconciliation    ✓ M11.2 6/6 PASS — Branch II level 2 (β-functions) CLOSED.   ==========================
- `m11_3_gamma_k.txt` — 6/6 PASS — context: …=====================================================================   Result: 6/6 PASS     [PASS] M11.3.1 M10.5.4 reproduction (1.244 entry)     [PASS] M11.3
- `m11_3_gamma_k.txt` — 6/6 PASS — context: …act (sub-percent)     [PASS] M11.3.6 H_0 tension structural non-help    ✓ M11.3 6/6 PASS — Branch II level 3 (γ(k) RG running) CLOSED.  =======================
- `m11_4_structural.txt` — 6/6 PASS — context: …=====================================================================   Result: 6/6 PASS     [PASS] M11.4.1 β=γ vacuum condition at FRG WF FP     [PASS] M11.4.
- `m11_4_structural.txt` — 6/6 PASS — context: …5 n=2 minimality (B.2)     [PASS] M11.4.6 Branch I/II RG convergence    ✓ M11.4 6/6 PASS — Branch II level 4 (KNOWN_ISSUES structural closure) CLOSED.  =======
- `m11_E_emergent_source.txt` — 6/6 PASS — context: …= ρ_loc(Φ) restore translation invariance + Goldstone? # Closure-grade target: ≥6/6 PASS ######################################################################

## Eksportowalne formuły (boxed)

**`M11_1_audit_results.md`:**

$$
S[\Phi] = \int d^3x\Bigl[ V_\text{onsite}(\Phi) + \tfrac12 J\,\Phi^2(\nabla\Phi)^2 \Bigr]
$$

**`M11_1_audit_results.md`:**

$$
V_\text{eff}^\text{1-loop}(\bar\Phi) = V_\text{onsite}(\bar\Phi) + \tfrac{T}{2}\!\!\int_\text{BZ}\frac{d^3q}{(2\pi)^3}\,\ln\bigl[K_\Phi(q;\bar\Phi) + A(\bar\Phi)\bigr]
$$

**`M11_1_audit_results.md`:**

$$
\beta_\text{eff} = 3 c_3,\qquad \gamma_\text{eff} = -4 c_4
$$

**`M11_1_audit_results.md`:**

$$
\ln K_\Phi(q;\Phi) = \ln[2J\Phi^2 (3 - \Sigma_\mu \cos q_\mu)] = 2\ln\Phi + \ln[2J(3 - \Sigma_\mu \cos q_\mu)]
$$

**`M11_1_audit_results.md`:**

$$
\bigl\langle \ln K_\Phi(q;\Phi)\bigr\rangle_\text{BZ} = 2\ln\Phi + C_\text{BZ},\qquad C_\text{BZ} \equiv \int_\text{BZ}\frac{d^3q}{(2\pi)^3}\ln\bigl[2J(3-\Sigma\cos q_\mu)\bigr]
$$

---

## Cross-references

- [[README.md]] — opis folderu + YAML status
- [[NEEDS.md]] — otwarte luki tego folderu
- [[meta/research/RESEARCH_BUS.md]] — broadcast tych findings
- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa