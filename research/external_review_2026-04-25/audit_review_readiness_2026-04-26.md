# Review-readiness consistency audit (A1)

**Date**: 2026-04-26
**Scope**: Pre-commit consistency check across `tgp-core-paper/paper/tgp_core.tex`,
`tgp-core-paper/KNOWN_ISSUES.md`, and supporting OP closure dossiers.
**Trigger**: Multiple Priority-1 patches (P1.2 RESOLVED, P1.4 RESOLVED) and
Phase 0+ OP-M92 commits accumulated 2026-04-25; audit before consolidating
into the P1.6 single commit.

---

## Verdict: 6 findings — 1 CRITICAL, 2 STALE, 2 WORDING, 1 SCOPE
## (+ Post-fix verification pass 2026-04-26 surfaced 1 new finding — also resolved)

> **Update 2026-04-26**: **All 6 findings RESOLVED**.
> Findings #1, #2 (KNOWN_ISSUES.md C1/C3 row refresh), #3 (orphan
> `prop:substrate-action` + bond exponent typo), #6 (Phase 0+
> paragraph), #4, #5 (path-citation hygiene). See per-finding
> resolution headers below for detail. P1.6 single-commit
> consolidation is the remaining milestone.

> **Update 2026-04-26 (post-fix verification pass)**: re-running
> the auditors after fixing Findings #1–#6 surfaced **NEW FINDING #7
> — this audit synthesis file was created at the wrong nested path
> during the prior session** (was at
> `TGP/tgp-core-paper/TGP/TGP_v1/research/external_review_2026-04-25/audit_review_readiness_2026-04-26.md`
> instead of the canonical
> `TGP/TGP_v1/research/external_review_2026-04-25/audit_review_readiness_2026-04-26.md`).
> The wrong location created a stray untracked `TGP/` subtree inside
> the `tgp-core-paper` git repo. **RESOLVED** by `mv` to the
> canonical path + `rmdir` of the empty parent chain inside
> tgp-core-paper. After the move, `git status` in tgp-core-paper
> shows only the 3 expected modifications (KNOWN_ISSUES.md, paper/
> tgp_core.pdf, paper/tgp_core.tex) and no stray untracked files.
> See section "Post-fix verification pass" at the bottom for the
> full re-audit.

The paper's *physics* is consistent with KNOWN_ISSUES.md and the OP closures.
However the **anchor surface** has drifted: claims about specific labels,
appendices, and review dispositions reference structures that either don't
exist in the tex (Finding #1, #3) or reflect outdated review state
(Finding #2, #6). Two further findings (#4, #5) are mechanical (path /
file-citation hygiene). A pre-publication reviewer following the breadcrumbs
in KNOWN_ISSUES.md → tex would hit a broken `\ref` and a missing appendix.

**Recommendation**: fix Findings #1, #3 BEFORE the P1.6 commit (these are
factually wrong claims). Findings #2, #6 should be folded into P1.6 as
disposition refreshes. Findings #4, #5 can be either deferred or batched.

---

## Audit checklist

| # | Task | Status |
|---|------|--------|
| 1 | C1–C6 dispositions (KNOWN_ISSUES vs paper text) | ✅ done |
| 2 | Cross-reference integrity (`\ref`/`\label`, file paths) | ✅ done |
| 3 | Wikilink integrity in KNOWN_ISSUES.md | ✅ done |
| 4 | Paper claims vs OP closures (OP-7, OP-EHT, OP-EHT-A, OP-M92) | ✅ done |
| 5 | Git status / uncommitted-changes reconciliation | ✅ done |
| 6 | Synthesis (this file) | ✅ done |

---

## Findings

### 🔴 FINDING #1 — CRITICAL: C1 disposition references non-existent labels

**Where**: `KNOWN_ISSUES.md:1585` (C1 row in C1-C6 disposition table)

**Claim** (verbatim):
> "Patched via new remarks `rem:B-ssb-v2` (dodatek B) and
> `rem:GL-breaking-axiomatic` (sek01)."

**Reality** (verified by grepping `paper/tgp_core.tex`):

- `rem:B-ssb-v2` — **NOT present**.
- `rem:GL-breaking-axiomatic` — **NOT present**.
- The paper does have `prop:Z2`, `ax:P3`, `prop:instability` (lines
  229–266) which functionally cover the SSB story, but they are not the
  remarks that the C1 row claims to exist.

**Impact**: an external reviewer cross-referencing C1 → tex anchor will
not find the cited labels. KNOWN_ISSUES.md is making a false claim
about paper state.

**Recommended fix**: rewrite the C1 disposition to reference what
*actually* exists in tex (`prop:Z2` + `ax:P3` + `prop:instability` for
the SSB story; `rem:alpha2-closure` for the v1→v2 GL bond pivot
discussion). Drop the references to `rem:B-ssb-v2` and
`rem:GL-breaking-axiomatic` unless those remarks are added in a
follow-up edit.

---

### 🔴 FINDING #3 — CRITICAL: orphan `\ref{prop:substrate-action}` and missing Appendix B  ✅ RESOLVED 2026-04-26

**Status**: RESOLVED via Option A (write Appendix B + fix bond exponent typo).

**Resolution summary**:
- `paper/tgp_core.tex` — bond in `eq:H-Gamma` (lines 213–215) corrected:
  `ŝ²ŝ²` → `ŝ⁴ŝ⁴` so the bond Taylor-expands to $K(\varphi)=\varphi^{4}$
  (matches `thm:alpha2`, `rem:alpha2-closure`, line 412 `K_ij = J(φ_iφ_j)²`).
- `paper/tgp_core.tex` — `\appendix` + `\setcounter{section}{1}` +
  Appendix B with `\label{app:substrate-action}` and the formal
  `\begin{proposition}\label{prop:substrate-action}` plus full
  block-averaging / gradient-expansion proof.
- Verification: `_audit_xref.py` re-run reports **0 orphan refs**.
  `pdflatex` builds 17-page PDF, no undefined references.
- Documentation: top-level entry added to
  `tgp-core-paper/KNOWN_ISSUES.md` (2026-04-26).

**Original investigation** (preserved below for record):

**Where**: 6 references in `tgp_core.tex` — lines 67, 162, 222, 304,
441, 917 — all of the form `App.~B, Prop.~\ref{prop:substrate-action}`.

**Reality**:
- `\label{prop:substrate-action}` — **NOT defined anywhere** in the tex.
- The document has **no `\appendix` command and no Appendix B section**
  (sections in the paper are `\section{Introduction}`,
  `\section{Axioms}`, `\section{From substrate to field}`,
  `\section{Effective field theory}`, `\section{Emergent metric}`,
  `\section{Three physical regimes}`,
  `\section{Status of theorems and open problems}`,
  `\section{Applications and falsifiability}`,
  `\section{Conclusion}`).

**Impact**: 6 broken cross-references. `\ref{prop:substrate-action}`
will render as `??` in compiled PDF. Critical because the proposition
is invoked to back-stop the GL gradient functional
$(K_{\mathrm{geo}}/2)\int\varphi^4(\nabla\varphi)^2 d^d x$ — i.e. the
v2 axiomatic substrate-to-continuum bridge.

**Recommended fix** (two options):
1. **Add Appendix B with `prop:substrate-action`** (the formal
   continuum-limit proposition that derives the GL gradient functional
   from the lattice bond — text is presumably ready in
   `TGP/TGP_v1/core/formalizm/dodatekQ_coarse_graining_formal.tex` or
   similar partial-proof file).
2. **Inline the proposition** at first invocation (line 222, in
   `ax:substrate`) and remove the App.~B references.

Option 1 is preferred for a peer-review-ready preprint: every claim
the body relies on should have a stated proof or reference.

---

### 🟡 FINDING #2 — STALE: C3 disposition row in KNOWN_ISSUES.md outdated

**Where**: `KNOWN_ISSUES.md:1587` (C3 row).

**Claim**:
> "**To be patched (P1.3)**: restrict to `γ_PPN = β_PPN = 1` from
> static/isotropic/weak-field ansatz. Full 10-PPN requires OP-7
> (moving matter + tensor sector)."

**Reality**: both halves of P1.3 are *already done*:

1. **Narrowing of `cor:ppn`**: paper line 692 explicitly scopes the
   corollary as "PPN parameters at 1PN, **within the ansatz of
   Theorem~\ref{thm:metric}**". The corollary itself states only
   `γ_PPN = β_PPN = 1` — no overreach.
2. **Full 10-PPN closure**: integrated into `rem:two-projections`
   (line 355–358):
   > "All 10 PPN parameters lie in experimental bounds (γ = β = 1
   > EXACT from M9.1″ P1; α₁₋₃ = ζ₁₋₄ = ξ_PPN = 0 from Z₂ + general
   > covariance, T6.1)."
   This is precisely what P1.3 was supposed to deliver, and it is
   delivered via OP-7 T6.1 (`OP7_T6_results.md`).

**Impact**: KNOWN_ISSUES.md C3 row signals outstanding work that has
already been completed.

**Recommended fix**: rewrite the C3 disposition in the same RESOLVED
form as C2/C4 — e.g.:
> "**RESOLVED 2026-04-26 (via OP-7 T6.1; see C4)**: `cor:ppn` narrowed
> to γ_PPN = β_PPN = 1 at 1PN within the Theorem~`thm:metric` ansatz
> (paper line 692); full 10-PPN closure delivered via OP-7 T6.1
> (`OP7_T6_results.md`), σ_ab corrections O(σ²)~10⁻⁶⁰; α₁₋₃ = ζ₁₋₄
> = ξ_PPN = 0 from Z₂ + general covariance. Integrated into
> `rem:two-projections` line 355–358."

---

### 🟡 FINDING #6 — STALE: paper OP-M92 row covers Phase 0 only; Phase 0+ work not integrated  ✅ RESOLVED 2026-04-26

**Status**: RESOLVED via three-edit refresh of `tgp_core.tex`.

**Resolution summary**:
- `tgp_core.tex:1012-1013` — status string upgraded to
  "OP-M92 (OPEN — Phase~0 + Phase~0+ shipped 2026-04-25;
  multi-source α-universality ISSUE FLAGGED)".
- `tgp_core.tex:~1117` — new ~30-line Phase 0+ summary paragraph
  appended to the OP-M92 description, covering all four cross-checks
  with quantitative results: variational scaffolding (5/5 POSITIVE),
  DESI cosmology (POSITIVE), WEP / 5th-force (POSITIVE TIGHT,
  MICROSCOPE η≈1.6e-16 vs 1.1e-15 bound, margin 6.7×), and
  multi-source α-universality (ISSUE FLAGGED — α_SI ∝ M_BH² across
  9 decades, candidate-independent, gates Phase 1).
- `tgp_core.tex:~1142-1144` — falsification-target footnote extended
  with "Phase~0 + Phase~0+ readiness package 2026-04-25, with a
  multi-source α-universality issue flagged for Phase~1".
- Verification: `pdflatex` builds cleanly (18 pages, was 17),
  `_audit_xref.py` reports 0 orphan refs.
- Documentation: top-level entry added to
  `tgp-core-paper/KNOWN_ISSUES.md` (2026-04-26).

**Original investigation** (preserved below for record):

**Where**: `tgp_core.tex:1012`, `1103`, `1142` — the paper consistently
says "Phase~0 readiness package shipped 2026-04-25".

**Reality**: KNOWN_ISSUES.md OP-M92 entry now reads "OPEN — PHASE 0
+ PHASE 0+ FULL (incl. multi-source consistency ISSUE FLAGGED)". Four
new commits between Phase 0 and HEAD:

| Commit | Phase 0+ artefact |
|---|---|
| `9b01bf3` | Candidate D structural sketch 5/5 POSITIVE |
| `e97cad0` | Candidate D vs DESI w(z) ≥ −1 cross-check (PASS) |
| `5068bfc` | Candidate D vs MICROSCOPE / LLR cross-check (PASS) |
| `19b6b38` | Multi-source α-universality **ISSUE FLAGGED** |

**Impact**: paper presents Candidate D as merely "leading the ranking"
(line 1108) without acknowledging that:
- Candidate D has now been numerically sketched (not just analytically
  promising);
- It passes a non-trivial cosmological cross-check (DESI w(z) ≥ −1);
- It passes WEP / 5th-force bounds;
- It has **a flagged self-consistency issue** in the multi-source
  α-universality sub-test.

The last item is the most important: a flagged issue in Candidate D's
multi-source behaviour deserves explicit mention so that an external
reviewer knows the M9.2 pivot is not yet free of internal tensions.

**Recommended fix**: extend paper's OP-M92 description (line 1103–1117)
with a one-paragraph Phase 0+ summary, ending in a forward reference to
`KNOWN_ISSUES.md` for the multi-source α-universality issue. This is a
~10-line edit, no figures or new equations.

---

### 🟠 FINDING #4 — WORDING: ambiguous path for `TGP_CLOSURE_PLAN_2026-04-25.md`  ✅ RESOLVED 2026-04-26

**Status**: RESOLVED. Line 989 now reads
`\texttt{research/op7/TGP\_CLOSURE\_PLAN\_2026-04-25.md}`,
matching the sibling citation on line 988.

**Original investigation** (preserved):

**Where**: `tgp_core.tex:989` — `\texttt{TGP\_CLOSURE\_PLAN\_2026-04-25.md}`.

**Reality**: file lives at
`TGP/TGP_v1/research/op7/TGP_CLOSURE_PLAN_2026-04-25.md`. The bare
basename in the texttt makes the citation ambiguous (is this in
`/`, `tgp-core-paper/`, `TGP_v1/research/op7/`, or somewhere else?).

**Recommended fix**: use the full path
`research/op7/TGP\_CLOSURE\_PLAN\_2026-04-25.md` to match the sibling
citation on line 988.

---

### 🟠 FINDING #5 — WORDING: path-prefix convention inconsistent in tex  ✅ RESOLVED 2026-04-26

**Status**: RESOLVED. Three outliers (lines 445, 921, 924) harmonized
to the dominant `research/...` convention. **Departure from this
audit's original recommendation**: I chose the dominant convention
(12 callsites already use `research/...`) over the audit-recommended
`TGP_v1/research/...` (only 3 callsites). Rationale: minimize diff
(3-line change vs 12-line change) and match the implicit
companion-repo mental cd. Documented in
`tgp-core-paper/KNOWN_ISSUES.md` (2026-04-26 entry on
"path-citation hygiene").

After fix: all 14 `\texttt{...research/...\}` callsites uniform;
zero `TGP_v1/` strings remain in `tgp_core.tex`. PDF builds clean,
xref auditor reports 0 orphan refs.

**Original investigation** (preserved):

**Where**: `tgp_core.tex` mixes two path-prefix conventions for
research files:
- **Without `TGP_v1/` prefix** (lines 988, 1022, 1023, 1024 etc.):
  `research/op7/...`, `research/op-eht/...`, `research/op-m92/...`,
  `research/casimir_mof/`, `research/muon_g_minus_2/`,
  `research/thermal_transport_molecular/`,
  `research/continuum_limit/cg_strong_numerical.py`.
- **With `TGP_v1/` prefix** (lines 445, 921, 924):
  `TGP_v1/research/op6/`, `TGP_v1/research/op6/m3a_block_rg_1d.py`.

Both refer to real files (the `op6/` directory does live at
`TGP_v1/research/op6/`). The reader's cd assumption — `TGP_v1/` as
project root or `TGP/` as project root — disambiguates which form
works, and it is currently undocumented.

**Recommended fix**: pick one convention (recommendation:
`TGP_v1/research/...` since op6 already uses it explicitly) and apply
it uniformly. Affects ~10 lines of tex, no semantic change.

---

### 🟢 RESOLVED items — no action needed

These were checked and found to be in good order:

- **C2 anchors** (line 479 `rem:U-truncation`, line 506
  `thm:field-eq`): match `KNOWN_ISSUES.md` C2 row exactly.
- **C4 anchors** (line 78–89 abstract footnote, line 308
  `rem:two-projections`, line 318 `eq:sigma-def`, line 338
  `eq:sigma-eom`, line 345 `eq:metric-extended`, line 736
  `cor:cGW`, line 1145 F3, line 1149 F4): all present, line
  numbers match `KNOWN_ISSUES.md` C4 row within ±1 line.
- **C5 / C6** (α=2 selection wording; p=1 defence): handled via
  `rem:alpha2-closure` (line 410). Body-text α=2 assertions
  (lines 64, 165, 370, 695, 1129) all flow through
  Theorem~`thm:alpha2` which is properly scoped to (C1)–(C3).
- **OP-7 status table** (line 970–989): row matches
  `OP7_T6_results.md` (94/97 = 96.9%, T6.6 ξ/G=1 EXACT).
- **OP-EHT status** (line 999): 13/18 = 72% PASS — matches.
- **OP-EHT-A status** (line 1000–1001): 8/12 = 67% PASS, NEGATIVE on
  naive proper-time — matches.
- **Wikilinks in `KNOWN_ISSUES.md`**: 19/19 resolve (smart resolver
  with index over `TGP/`).
- **Backtick file paths in `KNOWN_ISSUES.md`**: 113 references; all
  resolve (the one apparent miss `_results.md` is a regex artefact
  from the C2 row prose `+ .txt + _results.md`).

---

## Punch list (prioritised, for P1.6 commit)

| # | Severity | Action | Lines | ETA |
|---|----------|--------|-------|-----|
| 1 | 🔴 critical | Rewrite C1 disposition row to reference real labels (`prop:Z2`, `ax:P3`, `prop:instability`, `rem:alpha2-closure`); drop fictional `rem:B-ssb-v2`, `rem:GL-breaking-axiomatic` | KNOWN_ISSUES.md:1585 | 5 min |
| 2 | 🔴 critical | Either add Appendix B with `prop:substrate-action`, or inline the proposition and rewrite all 6 `App.~B, Prop.~\ref{prop:substrate-action}` callsites | tex:67, 162, 222, 304, 441, 917 + new appendix | 30–60 min if inlining; 1–2 h if writing full appendix |
| 3 | 🟡 stale | Refresh C3 row: "RESOLVED 2026-04-26 via OP-7 T6.1" with reference to `rem:two-projections` line 355 | KNOWN_ISSUES.md:1587 | 5 min |
| 4 | 🟡 stale | Add Phase 0+ summary paragraph to paper OP-M92 description (ends with multi-source α-universality ISSUE FLAGGED forward-ref to KNOWN_ISSUES.md) | tex:1103–1117 | 15 min |
| 5 | 🟠 wording | Add `research/op7/` prefix to `TGP_CLOSURE_PLAN_2026-04-25.md` citation | tex:989 | 1 min |
| 6 | 🟠 wording | Unify path-prefix convention (recommend `TGP_v1/research/...` everywhere) | tex: ~10 callsites | 10 min |

**Total estimated effort**: 1–3 hours depending on appendix scope.

---

## Implications for P1.6 single-commit strategy

The original P1.6 plan was: bundle P1.1–P1.5 patches into one commit
documenting the external review response. Audit findings extend this:

- **Add to P1.6**: Findings #1, #3 (factually wrong cross-references).
  These are not new physics; they are cross-reference hygiene that
  should not survive the same commit that claims "review C1–C4
  RESOLVED".
- **Fold into P1.6 if time permits**: Findings #2, #6 (disposition
  freshness). Otherwise carry into a follow-up `P1.7 audit-driven
  cleanup` commit.
- **Defer**: Findings #4, #5 (path-prefix hygiene). These don't
  affect review readability and can wait for a paper-formatting
  pass.

The audit confirms the paper's physics is sound; the gaps are
*paratext* — labels, appendices, paths. None of the six findings
question OP-7 closure, U(φ) bounded-below, σ_ab construction, or any
substantive claim.

---

## Post-fix verification pass (2026-04-26)

After applying fixes for all six original findings, a verification
pass was run to confirm no regressions and to surface any newly
introduced issues.

### Auditor re-runs

- **`_audit_xref.py`**: 59 labels, 33 refs, **0 orphan refs** (was 1
  before Finding #3 fix). 26 unused labels remain (all are internal
  appendix anchors or section/equation anchors not cross-referenced
  from elsewhere — expected and harmless).
- **`_audit_links.py`**: 41 bare references (basename-only filename
  citations resolving to a unique file via filesystem index). All
  are inside narrative prose in KNOWN_ISSUES.md changelog entries,
  not LaTeX citations. The 5 reported "missing" hits are regex
  artefacts from the auditor not unescaping LaTeX `\_` inside
  markdown backticks (e.g., reading
  `` `\texttt{research/op7/TGP\_CLOSURE\_PLAN\_2026-04-25.md}` `` as
  a literal path with backslash — the underlying file *does* exist,
  just at the un-escaped name). Auditor bug, not a content issue.
- **`pdflatex`**: clean compile, 18 pages, 661 KB. No undefined /
  multiply-defined reference warnings. Two-pass ref resolution
  succeeds.

### Regression spot-checks

| Finding | Expected anchor | Verified |
|---|---|---|
| #1 (C1 labels) | 10× `prop:Z2`/`ax:P3`/`prop:instability`/`rem:alpha2-closure` in KNOWN_ISSUES.md | ✅ 10 hits |
| #2 (C3 RESOLVED) | 17× `cor:ppn`/`rem:two-projections`/`OP-7 T6` | ✅ 17 hits |
| #3 (bond) | `\hat{s}_{i}^{4}\,\hat{s}_{j}^{4}` in tex | ✅ 3 hits (line 214 main + 1278, 1306 in Appendix B) |
| #3 (label) | `\label{prop:substrate-action}` | ✅ line 1271 |
| #4 (path prefix) | `research/op7/TGP\_CLOSURE\_PLAN\_2026-04-25.md` | ✅ line 989 |
| #5 (no `TGP_v1/` prefixes) | 0 `TGP_v1` matches in tex | ✅ 0 hits |
| #6 (status string) | "+ Phase~0+ shipped...ISSUE FLAGGED" in tex | ✅ line 1013–1014 |
| #6 (Phase 0+ paragraph) | new paragraph with all 4 cross-checks | ✅ lines 1117–1153 |
| #6 (falsification footnote) | "Phase~0 + Phase~0+ readiness package" | ✅ line 1181 |

### Git status review

`git -C TGP/tgp-core-paper status --porcelain`:
```
 M KNOWN_ISSUES.md
 M paper/tgp_core.pdf
 M paper/tgp_core.tex
```

Diff stat: `KNOWN_ISSUES.md +709 lines, tgp_core.tex +203 lines,
tgp_core.pdf regenerated (638→662 KB)`. No stray untracked files
(after FINDING #7 fix). All three modified files are expected.

### NEW FINDING #7 (RESOLVED)

**🔴 NEW FINDING #7 — STRAY: audit synthesis at wrong nested path** ✅ RESOLVED 2026-04-26

**Where**: `TGP/tgp-core-paper/TGP/TGP_v1/research/external_review_2026-04-25/audit_review_readiness_2026-04-26.md`
(doubly-prefixed path inside the `tgp-core-paper` git repo).

**Reality**: the file was created during the prior session by an
agent invocation that used a relative path like
`TGP/TGP_v1/research/external_review_2026-04-25/...` while its cwd
was already inside `TGP/tgp-core-paper/`. The result was a nested
copy at the wrong location, plus a stray untracked `TGP/` subtree
in the tgp-core-paper repo.

**Impact**: (i) the canonical path I had been citing in
`KNOWN_ISSUES.md` cross-references didn't actually point to a real
file; (ii) tgp-core-paper repo had a phantom `TGP/` directory
showing up in `git status`.

**Resolution**: `mv` to canonical path + `rmdir` of the empty parent
chain inside tgp-core-paper. Verified: file at
`TGP/TGP_v1/research/external_review_2026-04-25/audit_review_readiness_2026-04-26.md`
(16 KB), tgp-core-paper `git status` clean (only the 3 expected
modifications).

**Lesson**: when an agent creates new files using relative paths
that were originally drafted from the vault root, double-check
they don't end up nested inside the cwd of whatever tool issued the
mkdir/write. For this audit-related work specifically, future audit
artefacts should be created via absolute path or with explicit
verification of cwd.

### Bottom line

After all fixes (Findings #1–#7), the paper is in a clean
review-readiness state:
- 6 critical/stale/wording findings cleared.
- 1 newly-discovered nested-path artefact cleared.
- 0 orphan refs, clean PDF build, sane git diff (3 expected files).
- Audit synthesis file lives at canonical location.

**Ready for P1.6 single-commit consolidation.**
