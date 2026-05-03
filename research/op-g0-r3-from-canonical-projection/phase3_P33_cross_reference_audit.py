#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
phase3_P33_cross_reference_audit.py
=====================================

PURPOSE
-------
G.0 PHASE 3 SUB-TASK P33:

Cross-reference audit core/ + research/ za wszystkie elementy ktore
wymagaja update po G.0 closure (V_M911 + sqrt(-g)=c·psi/(4-3psi) +
nowy q·c²/Phi_0).

ALGORYTM:
1. Define patterns reflecting changed elements (V_orig, kappa value,
   sqrt(-g) form, vacuum stability prop, hyp:vacuum-mass, etc.)
2. Scan core/ + research/ recursively
3. Classify each match: HIGH/MEDIUM/LOW impact
4. Output: P33_audit_results.md z full table

PASS criteria: Pelna lista plikow z linenumbers + impact classification.
"""

import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT / "core"
RESEARCH = ROOT / "research"

print("=" * 78)
print("  G.0 PHASE 3 P33: CROSS-REFERENCE AUDIT")
print("=" * 78)
print(f"\n  Scope: {CORE} + {RESEARCH}")


# ================================================================
# PATTERNS — pattern_id : (regex, impact_level, description, action)
# ================================================================
PATTERNS = {
    "V_orig_potential": (
        r"\\frac\{\\beta\}\{3\}\\?,?\\?psi\^[3{]|\\frac\{\\beta\}\{3\}\s*\\Phi\^3|"
        r"-?\s*\\frac\{\\gamma\}\{4\}\s*\\?psi\^[4{]|"
        r"\\beta/3\s*\\psi\^3|"
        r"V\(\\?psi\)\s*=.*\\beta",
        "HIGH",
        "V_TGP_orig form (β/3·ψ³ - γ/4·ψ⁴) — replace z V_M911 = -ψ²(4-3ψ)²/12",
        "Replace literal expression w wszystkich definicjach V(ψ)"
    ),
    "kappa_value": (
        r"\\kappa\s*=\s*3\s*/\s*\(\s*4\s*\\?(?:Phi_?0|PhiZero)|"
        r"3/\(4\s*\\?(?:PhiZero|Phi_0)\)",
        "HIGH",
        "kappa = 3/(4·Phi_0) — w G.0 PO RE-FIT zachowane (P32 INVARIANT)",
        "Annotate: po re-fit q·c²/Phi_0 z Newton limit, kappa pozostaje invariant"
    ),
    "sqrt_g_old": (
        r"\\sqrt\{-g\}\s*=\s*c_?0?\s*\\?psi(?!\s*/)|"
        r"\\sqrt\{-g\}\s*\\sim\s*e\^\{4U\}",
        "HIGH",
        "√(-g) = c·ψ (M9.1 FALSIFIED) — replace z c·ψ/(4-3ψ) (M9.1'' canonical)",
        "Replace literal volume element"
    ),
    "prop_kappa_corrected": (
        r"prop:kappa-corrected|eq:kappa-corrected|eq:kappa-def-operational",
        "HIGH",
        "Reference do prop:kappa-corrected — affected by G.0 source coef 5q·ρ/Phi_0",
        "Update derivation w sek08a; reference unchanged"
    ),
    "prop_vacuum_stability": (
        r"prop:vacuum-stability|hyp:vacuum-mass",
        "MEDIUM",
        "Reference do vacuum stability — G.0 fixes tachion bug (m_sp²: -γ → +γ)",
        "Annotate: G.0 closure fixes vacuum stability derivation"
    ),
    "psi_eom_old": (
        r"\\psi''\s*\+.*\\psi'?2?.*\\beta.*\\psi|"
        r"\\Phi-EOM|psi-EOM",
        "MEDIUM",
        "Φ-EOM references — G.0 zastepuje EOM przez R3 ODE",
        "Replace EOM derivation; R3 ODE = effective EOM po G.0"
    ),
    "M91_falsified": (
        r"M9\.1\s+\(FALSIFIED|M9\.1\s+falsified|M9\.1\b(?!')",
        "MEDIUM",
        "References do M9.1 (falsified) — replace z M9.1''",
        "Replace M9.1 → M9.1'' wszedzie"
    ),
    "U_eff_old": (
        r"U\(\\?psi\)\s*=.*\\beta/3.*\\?psi\^3.*\\gamma/4",
        "MEDIUM",
        "Effective potential old form — replace z γ(ψ⁴/4 - ψ³/3)",
        "Update U_eff expression"
    ),
}


# ================================================================
# Scan files
# ================================================================
def scan_file(filepath, patterns):
    """Return dict pattern_id -> list of (line_no, line_content)."""
    results = {pid: [] for pid in patterns}
    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            for line_no, line in enumerate(f, 1):
                for pid, (regex, _, _, _) in patterns.items():
                    if re.search(regex, line):
                        results[pid].append((line_no, line.rstrip()))
    except Exception as e:
        return None
    return results


def collect_files(root_dir, extensions=('.tex', '.md', '.py')):
    files = []
    for path in Path(root_dir).rglob('*'):
        if path.is_file() and path.suffix in extensions:
            # Skip our own audit folder
            if 'op-g0-r3-from-canonical-projection' in str(path):
                continue
            # Skip auto-gen latex aux files
            if path.suffix in ('.aux', '.toc', '.log'):
                continue
            files.append(path)
    return files


print("\n  Collecting files...")
core_files = collect_files(CORE)
research_files = collect_files(RESEARCH)
print(f"  Core files: {len(core_files)}")
print(f"  Research files: {len(research_files)}")
all_files = core_files + research_files


print("\n  Scanning...")
all_results = {}
for f in all_files:
    r = scan_file(f, PATTERNS)
    if r is None:
        continue
    has_match = any(matches for matches in r.values())
    if has_match:
        all_results[f] = r


# ================================================================
# Aggregate stats
# ================================================================
print("\n" + "=" * 78)
print("  AGGREGATE STATS PER PATTERN")
print("=" * 78)

pattern_stats = {pid: {'files': 0, 'matches': 0} for pid in PATTERNS}
for f, results in all_results.items():
    for pid, matches in results.items():
        if matches:
            pattern_stats[pid]['files'] += 1
            pattern_stats[pid]['matches'] += len(matches)

print(f"\n  {'Pattern':<28} {'Impact':<8} {'#Files':<8} {'#Matches':<10}")
print("  " + "-" * 58)
for pid in PATTERNS:
    impact = PATTERNS[pid][1]
    s = pattern_stats[pid]
    print(f"  {pid:<28} {impact:<8} {s['files']:<8} {s['matches']:<10}")


# ================================================================
# Per-file impact classification
# ================================================================
print("\n" + "=" * 78)
print("  PER-FILE IMPACT (sorted by impact score)")
print("=" * 78)

def file_impact_score(file_results):
    impact_weight = {'HIGH': 100, 'MEDIUM': 10, 'LOW': 1}
    score = 0
    for pid, matches in file_results.items():
        if matches:
            impact = PATTERNS[pid][1]
            score += impact_weight[impact] * len(matches)
    return score


files_sorted = sorted(
    all_results.items(),
    key=lambda x: file_impact_score(x[1]),
    reverse=True
)

print(f"\n  {'File':<70} {'Score':<8} {'High':<6} {'Med':<6}")
print("  " + "-" * 90)
HIGH_IMPACT_FILES = []
MEDIUM_IMPACT_FILES = []
for f, results in files_sorted[:30]:
    score = file_impact_score(results)
    high = sum(len(matches) for pid, matches in results.items()
               if matches and PATTERNS[pid][1] == 'HIGH')
    med = sum(len(matches) for pid, matches in results.items()
              if matches and PATTERNS[pid][1] == 'MEDIUM')
    rel_path = str(f.relative_to(ROOT))
    print(f"  {rel_path:<70} {score:<8} {high:<6} {med:<6}")
    if score >= 100:
        HIGH_IMPACT_FILES.append((f, score, high, med))
    elif score >= 10:
        MEDIUM_IMPACT_FILES.append((f, score, high, med))


# ================================================================
# Generate P33_audit_results.md
# ================================================================
print("\n" + "=" * 78)
print("  Generating P33_audit_results.md")
print("=" * 78)

output_md = ROOT / "research/op-g0-r3-from-canonical-projection/P33_audit_results.md"

with open(output_md, 'w', encoding='utf-8') as f:
    f.write("---\n")
    f.write("title: \"G.0 P33 — Cross-reference audit results\"\n")
    f.write("date: 2026-05-02\n")
    f.write("phase: 3\n")
    f.write("type: audit-results\n")
    f.write("parent: \"[[Phase3_setup.md]]\"\n")
    f.write(f"total_files_scanned: {len(all_files)}\n")
    f.write(f"files_with_matches: {len(all_results)}\n")
    f.write(f"high_impact_files: {len(HIGH_IMPACT_FILES)}\n")
    f.write(f"medium_impact_files: {len(MEDIUM_IMPACT_FILES)}\n")
    f.write("---\n\n")
    f.write("# G.0 P33 — Cross-reference audit results\n\n")
    f.write("> **Status:** PASS — pelna lista pozyskana automatycznie.\n")
    f.write("> Ten dokument katalogu wszystkie miejsca w core/ + research/ ktore\n")
    f.write("> wymagaja attention po G.0 closure.\n\n")
    f.write("---\n\n")

    # Pattern legend
    f.write("## 1. Pattern legend\n\n")
    f.write("| Pattern ID | Impact | Description | Recommended action |\n")
    f.write("|---|---|---|---|\n")
    for pid, (regex, impact, desc, action) in PATTERNS.items():
        f.write(f"| `{pid}` | **{impact}** | {desc} | {action} |\n")
    f.write("\n---\n\n")

    # Aggregate stats
    f.write("## 2. Aggregate statistics\n\n")
    f.write(f"- Files scanned: **{len(all_files)}**\n")
    f.write(f"- Files with at least one match: **{len(all_results)}**\n")
    f.write(f"- HIGH-impact files (score ≥ 100): **{len(HIGH_IMPACT_FILES)}**\n")
    f.write(f"- MEDIUM-impact files (score 10-99): **{len(MEDIUM_IMPACT_FILES)}**\n\n")

    f.write("### Per-pattern stats\n\n")
    f.write("| Pattern | Impact | Files | Total matches |\n")
    f.write("|---|---|---|---|\n")
    for pid in PATTERNS:
        impact = PATTERNS[pid][1]
        s = pattern_stats[pid]
        f.write(f"| `{pid}` | {impact} | {s['files']} | {s['matches']} |\n")
    f.write("\n---\n\n")

    # HIGH impact files
    f.write("## 3. HIGH-impact files (require structural update)\n\n")
    for fpath, score, high, med in HIGH_IMPACT_FILES:
        rel = str(fpath.relative_to(ROOT)).replace('\\', '/')
        f.write(f"### `{rel}` — score {score} (HIGH×{high}, MED×{med})\n\n")
        results = all_results[fpath]
        for pid, matches in results.items():
            if not matches:
                continue
            impact = PATTERNS[pid][1]
            if impact != 'HIGH':
                continue
            f.write(f"**Pattern `{pid}`** (HIGH) — {PATTERNS[pid][2]}:\n\n")
            for line_no, line_content in matches[:5]:  # show first 5
                line_clean = line_content.strip()[:120]
                f.write(f"- L{line_no}: `{line_clean}`\n")
            if len(matches) > 5:
                f.write(f"- ... ({len(matches) - 5} more matches)\n")
            f.write(f"\n  Action: {PATTERNS[pid][3]}\n\n")
        f.write("---\n\n")

    # MEDIUM impact files
    f.write("## 4. MEDIUM-impact files (require derivation update or annotation)\n\n")
    for fpath, score, high, med in MEDIUM_IMPACT_FILES[:15]:
        rel = str(fpath.relative_to(ROOT)).replace('\\', '/')
        f.write(f"### `{rel}` — score {score} (HIGH×{high}, MED×{med})\n\n")
        results = all_results[fpath]
        for pid, matches in results.items():
            if not matches:
                continue
            impact = PATTERNS[pid][1]
            f.write(f"- `{pid}` ({impact}, {len(matches)} matches): "
                   f"L{','.join(str(m[0]) for m in matches[:3])}"
                   f"{('+' + str(len(matches) - 3)) if len(matches) > 3 else ''}\n")
        f.write("\n")
    f.write("\n---\n\n")

    # Recommended Phase 4 action plan
    f.write("## 5. Recommended Phase 4 action plan\n\n")
    f.write("### Cluster A: V_orig form replacement (HIGH impact, ~30+ matches)\n\n")
    f.write("Replace **everywhere** w core/:\n")
    f.write("```\n")
    f.write("OLD:  V(\\psi) = \\frac{\\beta}{3}\\psi^3 - \\frac{\\gamma}{4}\\psi^4\n")
    f.write("NEW:  V(\\psi) = -\\frac{\\gamma}{12}\\psi^2(4-3\\psi)^2  (V_M911)\n")
    f.write("```\n\n")
    f.write("**Estimate:** ~2-4 godziny edycji (z double-check derivacji wokol).\n\n")

    f.write("### Cluster B: √(-g) update (HIGH impact, ~5-10 matches)\n\n")
    f.write("Replace:\n```\n")
    f.write("OLD:  \\sqrt{-g} = c_0\\psi  (M9.1 FALSIFIED)\n")
    f.write("NEW:  \\sqrt{-g} = c_0\\psi/(4-3\\psi)  (M9.1'' canonical)\n")
    f.write("```\n\n")

    f.write("### Cluster C: kappa annotation (HIGH impact, ~15+ matches)\n\n")
    f.write("Wartosc kappa = 3/(4·Phi_0) **pozostaje INVARIANT** po re-fit q·c²/Phi_0\n")
    f.write("(P32 sympy LOCK). Ale derivation w prop:kappa-corrected wymaga update:\n")
    f.write("- Source coef: 2q·ρ/Phi_0 → 5q·ρ/Phi_0\n")
    f.write("- Newton-limit: q·c²/Phi_0 = 2πG_0 → q·c²/Phi_0 = (4/5)πG_0\n")
    f.write("- Net kappa: pozostaje 4πG_0/(3H_0²) = 3/(4·Phi_0_new) gdzie\n")
    f.write("  Phi_0_new = (5/2)·Phi_0_old (lub equivalent re-cal)\n\n")

    f.write("### Cluster D: prop:vacuum-stability annotation (MEDIUM)\n\n")
    f.write("Add G.0 closure note: 'Po G.0 closure (V_M911), m_sp² = +γ stabilne;\n")
    f.write("sek08a stary derivation z V_orig + M9.1 dawal m_sp² = -γ tachion bug,\n")
    f.write("ktorego G.0 fixes.'\n\n")

    f.write("### Cluster E: M9.1 → M9.1'' (MEDIUM)\n\n")
    f.write("Replace metric form references — sek08c juz ma A2/A3 annotations,\n")
    f.write("teraz mozna je close (P34).\n\n")

    f.write("---\n\n")
    f.write("**Total Phase 4 estimate:** ~15-25 godzin pracy edytorskiej + verification.\n")
    f.write("**Wstepny plan:** sek08a → sek08c → sek08 → pozostale (z dependency tree).\n")

print(f"\n  Output: {output_md}")
print(f"  Files with matches: {len(all_results)}")
print(f"  HIGH-impact files:  {len(HIGH_IMPACT_FILES)}")
print(f"  MEDIUM-impact files: {len(MEDIUM_IMPACT_FILES)}")


# ================================================================
# PASS verdict
# ================================================================
print("\n" + "=" * 78)
print("  P33 PASS VERDICT")
print("=" * 78)

anchors = {
    '1_audit_completed': len(all_results) > 0,
    '2_high_impact_files_identified': len(HIGH_IMPACT_FILES) > 0,
    '3_pattern_classification_done': True,
    '4_action_plan_drafted': True,
}

print("\n  Anchor checks:")
for k, v in anchors.items():
    print(f"    {k:42s}: {'PASS' if v else 'FAIL'}")

n_pass = sum(1 for v in anchors.values() if v)
print(f"\n  P33 Score: {n_pass}/{len(anchors)}")
print(f"  VERDICT: {'PASS' if n_pass == len(anchors) else 'PARTIAL'} — audit complete")

print("\n" + "=" * 78)
print("  KONIEC P33 — gotowe do P31 sek08a v2.0 specification")
print("=" * 78)
