#!/usr/bin/env python3
"""
identify_projection_cycles.py — Read-only triage scan dla retrofit metodologicznego 2026-05-10.

Scans research/op-*/ folders for "L2-only output" drift markers (per autor diagnosis 2026-05-10):
  - L2 markery (parameters w obcych frameworkach): β_ppE, β_PPN, γ_PPN, ξ_2, ppE, PPN, α_ppE
  - L1 native markery (observables w jednostkach fizycznych): arcsec, Hz, ms, strain, deflection,
    Shapiro, perihelion, φ(f), h(f), BBN fraction, σ_8

Klasyfikacja per cykl:
  - NATIVE_CLEAN          : L1 markers present (observable target jasny)
  - PROJECTION_SUSPECTED  : L2 markers present, brak L1 markers (drift suspect)
  - INTENTIONAL_PROJECTION: in known whitelist (np. op-GWTC3-reanalysis musiał użyć β-prior)
  - STRUCTURAL_OR_OTHER   : żadne markers (algebra, axiom verification, etc.)

Output: 3 CSV files w meta/triage/ + markdown summary do stdout.

Usage:
    python tooling/identify_projection_cycles.py
    python tooling/identify_projection_cycles.py --output-dir meta/triage_2026-05-10

Read-only — NIE modyfikuje żadnego YAML ani README.
"""

import argparse
import csv
import io
import re
import sys
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

VAULT_ROOT = Path(__file__).resolve().parent.parent
RESEARCH_DIR = VAULT_ROOT / "research"
SKIP_FOLDERS = {"_archive", "_sandbox", "_deprecated"}

# L2 markery — paramters w obcych frameworkach (drift indicators jeśli primary)
L2_PATTERNS = [
    r"\bbeta_ppE\b", r"\bβ_ppE\b", r"\bbeta_\{ppE\}\b",
    r"\balpha_ppE\b", r"\bα_ppE\b",
    r"\bbeta_PPN\b", r"\bβ_PPN\b", r"\bbeta_\{PPN\}\b",
    r"\bgamma_PPN\b", r"\bγ_PPN\b", r"\bgamma_\{PPN\}\b",
    r"\bxi_PPN\b", r"\bξ_PPN\b",
    r"\bb_ppE\b",
    r"\bppE\b", r"\bPPN\b",
    r"\bYunes-?Pretorius\b",
    r"\bppE basis\b", r"\bPPN basis\b",
    r"\bppE chart\b", r"\bPPN chart\b",
    r"\bWill parameters\b", r"\bWill basis\b",
]

# L1 native markery — observables w jednostkach fizycznych (native primary indicators)
L1_PATTERNS = [
    r"\barcsec\b", r"\barc-?sec\b", r"\barcsecond\b",
    r"\bmilliarcsec\b", r"\bμas\b", r"\bmas\b",
    r"\bHz\b(?!.{0,20}rate)",  # Hz but not "rate Hz"
    r"\bkHz\b", r"\bMHz\b",
    r"\bmicrosecond\b", r"\bμs\b", r"\bms\b",
    r"\bdimensionless strain\b", r"\bstrain ratio\b", r"\bstrain h\(",
    r"\bdeflection angle\b", r"\bdeflection of light\b",
    r"\bShapiro delay\b", r"\bShapiro time\b",
    r"\bperihelion shift\b", r"\bperihelion advance\b", r"\bperihelion precession\b",
    r"\bNordtvedt parameter\b",
    r"\bφ\(f\)\b", r"\bphi\(f\)\b", r"\bh\(f\)\b",
    r"\bphase deviation\b", r"\bphase shift\b",
    r"\bBBN.*fraction\b", r"\bD/H\b", r"\bHe-?4\b", r"\bLi-?7\b",
    r"\bsigma_8\b", r"\bσ_8\b",
    r"\bw_eff\b", r"\bw\(z\)\b",
    r"\bH\(z\)\b", r"\bH_0\b",
    r"\borbital decay\b",
    r"\beclipse timing\b",
    r"\bredshift z\b",
    r"\bspectral line shift\b",
]

# Whitelist intencjonalne projekcje (cykle gdzie L2 jest celem, NIE drift)
INTENTIONAL_PROJECTION_CYCLES = {
    "op-GWTC3-reanalysis",  # MUST use β-prior — LIGO data filtered tym priorem
    "op-ppE-mapping",  # explicit translation cycle (po 2026-05-10 reklasyfikacja)
    "op-PPN-from-derivation",  # nazwa wskazuje translation cycle
}

L2_RE = re.compile("|".join(L2_PATTERNS), re.IGNORECASE)
L1_RE = re.compile("|".join(L1_PATTERNS))  # case-sensitive (units like Hz, ms matter)

YAML_FOLDER_STATUS_RE = re.compile(
    r"^\s*folder_status:\s*[\"']?(?P<value>[\w\-]+)[\"']?\s*$",
    re.MULTILINE,
)
YAML_DATE_RE = re.compile(
    r"^\s*date:\s*[\"']?(?P<value>[\d\-]+)[\"']?\s*$",
    re.MULTILINE,
)
YAML_VERDICT_RE = re.compile(
    r"^\s*verdict:\s*[\"']?(?P<value>[A-Z_]+)[\"']?\s*$",
    re.MULTILINE,
)


def scan_cycle(cycle_dir: Path) -> dict:
    """Zwróć dict z metrykami cyklu."""
    name = cycle_dir.name
    readme = cycle_dir / "README.md"
    if not readme.exists():
        return {
            "cycle": name,
            "category": "NO_README",
            "l1_count": 0,
            "l2_count": 0,
            "folder_status": "?",
            "date": "?",
            "verdict": "?",
            "evidence_l1": "",
            "evidence_l2": "",
        }

    text = readme.read_text(encoding="utf-8", errors="replace")

    # Extract YAML metadata
    folder_status = "?"
    date = "?"
    verdict = "?"
    if m := YAML_FOLDER_STATUS_RE.search(text):
        folder_status = m.group("value")
    if m := YAML_DATE_RE.search(text):
        date = m.group("value")
    if m := YAML_VERDICT_RE.search(text):
        verdict = m.group("value")

    # Count markers
    l1_matches = L1_RE.findall(text)
    l2_matches = L2_RE.findall(text)
    l1_count = len(l1_matches)
    l2_count = len(l2_matches)

    # Sample evidence (unique, up to 5)
    l1_unique = list(dict.fromkeys(l1_matches))[:5]
    l2_unique = list(dict.fromkeys(l2_matches))[:5]

    # Categorize
    if name in INTENTIONAL_PROJECTION_CYCLES:
        category = "INTENTIONAL_PROJECTION"
    elif l1_count == 0 and l2_count == 0:
        category = "STRUCTURAL_OR_OTHER"
    elif l1_count >= 2:  # solid L1 evidence (≥2 different observable markers)
        if l2_count > l1_count * 2:
            category = "MIXED_L1_L2_HEAVY_PROJECTION"
        else:
            category = "NATIVE_CLEAN"
    elif l2_count >= 2 and l1_count <= 1:
        category = "PROJECTION_SUSPECTED"
    else:
        category = "STRUCTURAL_OR_OTHER"

    return {
        "cycle": name,
        "category": category,
        "l1_count": l1_count,
        "l2_count": l2_count,
        "folder_status": folder_status,
        "date": date,
        "verdict": verdict,
        "evidence_l1": "|".join(l1_unique),
        "evidence_l2": "|".join(l2_unique),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--output-dir",
        default="meta/triage_2026-05-10",
        help="Output directory for CSV files (default: meta/triage_2026-05-10)",
    )
    ap.add_argument("--quiet", action="store_true", help="Suppress per-cycle output")
    args = ap.parse_args()

    output_dir = VAULT_ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    if not RESEARCH_DIR.exists():
        print(f"❌ Research dir not found: {RESEARCH_DIR}", file=sys.stderr)
        return 2

    results = []
    for cycle_dir in sorted(RESEARCH_DIR.iterdir()):
        if not cycle_dir.is_dir() or cycle_dir.name in SKIP_FOLDERS:
            continue
        if cycle_dir.name.startswith("_"):
            continue
        results.append(scan_cycle(cycle_dir))

    # Bucket by category
    buckets = {
        "NATIVE_CLEAN": [],
        "PROJECTION_SUSPECTED": [],
        "INTENTIONAL_PROJECTION": [],
        "MIXED_L1_L2_HEAVY_PROJECTION": [],
        "STRUCTURAL_OR_OTHER": [],
        "NO_README": [],
    }
    for r in results:
        buckets[r["category"]].append(r)

    # Write CSVs
    fieldnames = [
        "cycle",
        "category",
        "l1_count",
        "l2_count",
        "folder_status",
        "date",
        "verdict",
        "evidence_l1",
        "evidence_l2",
    ]
    for category, rows in buckets.items():
        if not rows:
            continue
        path = output_dir / f"{category}.csv"
        with path.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)

    # Summary to stdout
    print("# Triage scan results — identify_projection_cycles.py")
    print(f"\nScanned: `{RESEARCH_DIR.relative_to(VAULT_ROOT)}` ({len(results)} cycles)")
    print(f"Output dir: `{output_dir.relative_to(VAULT_ROOT)}/`\n")
    print("| Category | Count | CSV file |")
    print("|---|---|---|")
    for category in [
        "NATIVE_CLEAN",
        "PROJECTION_SUSPECTED",
        "MIXED_L1_L2_HEAVY_PROJECTION",
        "INTENTIONAL_PROJECTION",
        "STRUCTURAL_OR_OTHER",
        "NO_README",
    ]:
        cnt = len(buckets[category])
        csv_link = f"`{category}.csv`" if cnt else "—"
        print(f"| {category} | {cnt} | {csv_link} |")

    # Top-priority items: PROJECTION_SUSPECTED
    if buckets["PROJECTION_SUSPECTED"]:
        print("\n## PROJECTION_SUSPECTED (drift candidates — manual triage required)\n")
        print("| Cycle | folder_status | verdict | L2 markers (sample) |")
        print("|---|---|---|---|")
        for r in sorted(buckets["PROJECTION_SUSPECTED"], key=lambda x: x["cycle"]):
            print(
                f"| {r['cycle']} | {r['folder_status']} | {r['verdict']} | "
                f"{r['evidence_l2'][:80]} |"
            )

    if buckets["MIXED_L1_L2_HEAVY_PROJECTION"]:
        print("\n## MIXED_L1_L2 (z heavy projection bias — review carefully)\n")
        print("| Cycle | L1/L2 | folder_status | verdict |")
        print("|---|---|---|---|")
        for r in sorted(buckets["MIXED_L1_L2_HEAVY_PROJECTION"], key=lambda x: x["cycle"]):
            print(
                f"| {r['cycle']} | {r['l1_count']}/{r['l2_count']} | "
                f"{r['folder_status']} | {r['verdict']} |"
            )

    print("\n**Next step:** manual triage decisions per cycle, output to "
          "`meta/PROJECTION_TRIAGE_2026-05-10.md`")
    return 0


if __name__ == "__main__":
    sys.exit(main())
