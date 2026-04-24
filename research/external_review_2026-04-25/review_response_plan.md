# Response plan — external review 2026-04-25

**Trigger:** external agent review of TGP v1 paper delivered 2026-04-25.
**Scope:** concrete action items structured by priority, with acceptance
criteria and estimated cost. This is the master plan; per-priority
deliverables live in sibling files / other research folders.

---

## Review summary (pro / contra)

### Pro
- Explicit split `twierdzenie / numeryka / hipoteza` gives paper
  shape of "claim-plus-caveat" rather than "theory of everything".
- Kinetic operator `D_kin[φ] = ∇²φ + 2(∇φ)²/φ = (1/3φ²)∇²(φ³)`
  best piece: conditional theorem, α=2 selected within the class.
- `research/` folder well-organized; **negative checks** pattern
  (muon g-2, Wiedemann–Franz, Casimir) is methodologically healthy.

### Contra (ordered by severity)

| # | Issue | Status |
|---|-------|--------|
| C1 | Substrate looks stable at `s=0` → SSB seems impossible | **FALSE** (reviewer misread; axiom explicitly states `m₀² < 0` in ordered phase). But a **deeper variant is real**: in v2 the GL bond is purely gradient, so SSB depends ENTIRELY on `m₀² < 0` — no J-driven alternative. Paper doesn't flag this. |
| C2 | `U(φ) = β/3 φ³ − γ/4 φ⁴` with γ>0 → `−∞` at large φ | Real presentation issue. Full `V_eff(Φ)` is bounded below by bare `λ₀/4 Φ²`; the cubic+quartic form is a Taylor truncation valid near vacuum. Paper must state this explicitly. |
| C3 | "All 10 PPN = GR" overreach | Real. Static/isotropic ansatz gives `γ_PPN = β_PPN = 1`; the other 8 require moving matter + tensor sector (OP-7). |
| C4 | `c_GW = c₀` overreach | Real and physically serious. TGP has 1 scalar φ ⇒ at most scalar GW. GR's 2 tensor polarisations require OP-7 closure. Potential LIGO/Virgo tension (scalar component bounded at few %). |
| C5 | α=2 pivot — "moving the goalposts" risk | Partly trafne. KNOWN_ISSUES already records the pivot honestly, but core paper could be more radical: say "α=2 is selected within GL-substrate ansatz; minimal-substrate derivation falsified and retracted". |
| C6 | `p=1` defence mixes external and internal arguments | Minor. Paper should separate falsifiable evidence from TGP-internal consistency checks. |

### Cross-connection to on-going OP-2b investigation

Reviewer's C5 ("pivot moved physics from result to assumption") is
the **same class of issue** as the user's conjecture about MK-RG in
ŝ-variables losing the composite-field Jacobian:

    "we made a variable choice (ŝ) at a place where the result (β=γ)
     depends critically on the opposite choice (Φ, with explicit
     Jacobian)".

Both are cases where a **variable/axiom choice silently absorbs
physical content that used to be dynamical**. Fixing C5 and running
Test A (MK-RG in Φ-variables) are logically coupled.

---

## Priority 1 — paper-level clarifications (days)

All are edits to existing tex/md files. No new physics; just
correctly labelling what is assumption vs. result and bounding the
validity domain of current theorems.

### P1.1 — Substrate axiom patch (v2 SSB mechanism)

**File:** `axioms/substrat/dodatekB_substrat.tex` (+ mirror in
`core/sek01_ontologia/sek01_ontologia.tex` and `KNOWN_ISSUES.md`).

**Patch:** Add a `\begin{remark}[SSB mechanism in v2]` block after
the H_Gamma axiom, stating:

- In v1 the bond `−J Σ ŝ_i ŝ_j` was ferromagnetic and could drive
  SSB even at `m₀² > 0` via Jz > m₀² condition.
- In v2 the GL bond `J A_ij Φ_i Φ_j (Φ_j − Φ_i)²` vanishes on any
  uniform configuration, hence contributes NOTHING to mean-field v².
- Therefore v² = |m₀²|/λ₀ (not `(Jz − m₀²)/λ₀`) and SSB depends
  entirely on `m₀²(T_sub) < 0` below T_c.
- The existence of `T_c` with a sign-change of `m₀²` is itself
  axiomatic input, not a derived consequence of the dynamics.

Also update `v² = (Jz − m₀²)/λ₀` → `v² = |m₀²|/λ₀` everywhere
(currently wrong-for-v2 at `dodatekB_substrat.tex:300` and wherever
cross-referenced).

**Acceptance:** reviewer cannot claim "substrate dodatnio określony
wokół s=0 → brak SSB" without explicitly contradicting this remark.

### P1.2 — `U(φ)` validity domain

**File:** `core/sek08_formalizm.tex` (or wherever `U(φ)` is first
stated) + `M1_potential_inventory.md` if cross-reference needed.

**Patch:** Add a remark stating:

- `U(φ) = β/3 φ³ − γ/4 φ⁴` is a **Taylor truncation of V_eff(Φ)**
  around `Φ₀`, valid in the domain `|φ − 1| ≪ 1`.
- Full stability is on the level of `V_eff(Φ)`, which is bounded
  below by the bare `λ₀/4 Φ² > 0` at large Φ.
- The apparent unboundedness of U for large φ is a truncation
  artefact, not a feature of the theory.
- Cite M2a eq. (2.5) for the full V_eff form including `(T/2) ln Φ`.

**Acceptance:** paper no longer looks "globally unstable" at
face-value inspection.

### P1.3 — PPN claim scope

**File:** `core/` wherever PPN is claimed (TBD — find in P1.3
execution).

**Patch:** Restrict claim from "all 10 PPN = GR" to:
- Static, isotropic, weak-field ansatz: `γ_PPN = β_PPN = 1`.
- Remaining 8 parameters (ξ, α_{1,2,3}, ζ_{1,2,3,4}) require
  analysis of moving matter and tensor sector (OP-7).
- Flag this as OPEN in `KNOWN_ISSUES.md` under a new entry
  `2026-04-25: PPN-10 scope narrowed pending OP-7`.

**Acceptance:** no over-claim about Cassini/LLR/binary-pulsar
constraints beyond what γ,β alone cover.

### P1.4 — Gravitational waves

**File:** wherever `c_GW = c₀` is claimed (likely `core/sek08` or
a dedicated GW subsection).

**Patch:** Restrict to:
- Perturbations on `g_eff` (scalar φ fluctuations) propagate
  luminally (`c = c₀`).
- **Two-polarisation tensor GW of GR requires OP-7 closure**; TGP
  with a single scalar φ gives at most a scalar-type GW, and GR's
  transverse-traceless modes are not currently represented.
- Consequently GW150914 and similar LIGO/Virgo detections are
  NOT yet a closed prediction of TGP; flag as `OP-7 critical`.
- If OP-7 cannot produce two tensor polarisations, TGP is
  falsified by existing LIGO/Virgo scalar-mode bounds
  (`< few %` of tensor amplitude).

**Acceptance:** honest positioning that OP-7 is the GW-make-or-break
milestone.

### P1.5 — α=2 as selection

**File:** `KNOWN_ISSUES.md` (where pivot is already documented) +
front-matter of `core/sek08` or wherever α=2 is claimed.

**Patch:** One clear sentence at every α=2 assertion:

> "In the v2 axiomatic (GL-substrate), α=2 is selected within the
> class of Φ-covariant local second-order operators via conditions
> C1–C3. The earlier v1 derivation from a minimal bilinear bond
> substrate has been falsified (KNOWN_ISSUES 2026-04-24) and is
> retracted."

**Acceptance:** no reader can reasonably believe that α=2 is
"derived from first principles" in the minimal-substrate sense.

### P1.6 — Commit & push P1 batch

Single commit with all Priority-1 patches: message
`docs(core,axioms) review-2026-04-25: narrow claims + mark pivot impact`.
Push.

---

## Priority 2 — Test A: MK-RG in Φ-variables (days)

**Goal:** Test the user's conjecture that MK-RG in ŝ-variables
loses the composite-field Jacobian, which is exactly what drove
tree-level `β = γ`. If Φ-variable MK-RG preserves the identity, OP-2b
is closed and `thm:beta-eq-gamma-triple` is restored (possibly in
modified form).

**Script:** `research/op1-op2-op4/mk_rg_phi.py`.

Differences from `mk_rg_bgamma.py`:

1. Track operator basis
       V(Φ) = Σ_k c_{2k}/(2k) · Φ^k + μ ln Φ
   including explicit `μ ln Φ` term (tree: μ = T/2).
2. Measure `Φ^{-1/2}` becomes `Φ^{-1/2} e^{-V(Φ)/T}` (the full
   marginal).
3. MK-RG recursion: bond-move + decimation in Φ-space. Requires
   deriving the Φ-variable analogues of
       c_{2k}' = c_{2k} − 2 F_{2k}/(2k−1)!
   with proper treatment of the `ln Φ` and measure terms.
4. Flow μ(ℓ) and check whether it's relevant, marginal, or
   irrelevant at WF.

**Deliverables:**
- `mk_rg_phi.py` — Φ-variable MK-RG implementation.
- `mk_rg_phi_results.txt` — raw output.
- `M4_results.md` — interpretation, comparison with M3 (ŝ-version).
- If Test A succeeds (μ marginal/relevant, `β/γ → 1` at WF):
  update `thm:beta-eq-gamma-triple` to a restricted form;
  document closure of OP-2b.
- If Test A fails (μ irrelevant, `β/γ → −0.57` recovered):
  confirms OP-2b is fundamentally open; proceed to P3.

**Estimated cost:** 1-2 days (analytical derivation + implementation
+ validation against Gaussian sublimit).

---

## Priority 3 — Deferred technical extensions (weeks)

### P3.1 — Z_Φ wave-function renormalisation in MK-RG
Track the composite-field 2-point function at each RG step; multiply
`Φ`-coefficients by `Z_Φ^{appropriate power}`. Known-gap #2 in
M3_results.md §7.3.

### P3.2 — GL-bond operator in MK-RG
Extend single-site MK-RG to include explicit bond-operator
`G Φ_i Φ_j (ΔΦ)²` in the flow. Known-gap #1 in M3_results.md §7.3.

### P3.3 — OP-7 closure (tensor sector)
Derive whether TGP's single scalar φ can produce 2 tensor
polarisations via emergent mechanism (e.g., spin-2 bound state,
auxiliary σ_ab field, geometric projection). This is THE GW
make-or-break milestone after P1.4 repositioning.

### P3.4 — NPRG (nonperturbative RG)
Last-resort cross-check if P2 + P3.1 + P3.2 all fail to restore
`β = γ` structurally.

**Order:** P3.1 and P3.2 can run in parallel (independent tech
debts). P3.3 is orthogonal and can be opened anytime. P3.4 only
if P2+P3.1+P3.2 all negative.

---

## Acceptance criteria for overall response closure

External review is considered "addressed" when:

- [ ] C1–C6 all have either a paper-level patch (P1) or a flagged
      open problem with clear scope (P3.3, P3.4).
- [ ] Paper no longer contains any sentence a reviewer can quote
      as "over-claim" without directly contradicting a stated
      validity domain.
- [ ] KNOWN_ISSUES.md has a 2026-04-25 entry listing all six
      critiques and their disposition (patched / open / rejected
      with reasoning).
- [ ] OP-2b either closed (Test A succeeds) or confirmed open with
      a specific path forward (P3.1, P3.2, or P3.4).

---

## Files in this folder

- `review_response_plan.md` — this file.
- `review_audit_m0sq_sign.md` — factcheck of C1 (reviewer misread).
  (To be written during P1.1 execution.)
- Per-priority deliverables land in appropriate existing folders:
  - P1 patches → `core/`, `axioms/`, `KNOWN_ISSUES.md`.
  - P2 (Test A) → `research/op1-op2-op4/mk_rg_phi.py`, `M4_results.md`.
  - P3 → new research folders as needed.
