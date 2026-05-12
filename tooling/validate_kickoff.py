#!/usr/bin/env python3
"""
validate_kickoff.py — Technical enforcement gate for new research cycles (post-2026-05-11)

Implementation of Rec 4 z external review 2026-05-11. Checks that new
`research/op-*-YYYY-MM-DD/README.md` files contain the mandatory `contract::` block
i ## §0.4 Pre-flight methodology read confirmation section per
`meta/CYCLE_KICKOFF_TEMPLATE.md` (BINDING dla cykli post-2026-05-10).

Rationale:
  External review 2026-05-11 wykazała że 0/7 cykli 2026-05-11 cohort miało
  mandatory `contract::` block, mimo że CYCLE_KICKOFF_TEMPLATE jest BINDING
  od 2026-05-10. Brak technical enforcement → procedural drift.

  Per `meta/AUDIT_2026-05-11_sympy_substance.md` §4.3:
  - 24/104 testów = literal `T_pass = True` w 2026-05-11 cohort sympy
  - 0/112 testów = FIRST_PRINCIPLES derivation z TGP axioms
  - 0/7 cykli ma `contract::` blok

Usage:
  python tooling/validate_kickoff.py research/op-NAME-YYYY-MM-DD/README.md
  python tooling/validate_kickoff.py --all-new   # validate wszystkie cycles z date >= 2026-05-10
  python tooling/validate_kickoff.py --pre-commit  # used by git hook

Exit codes:
  0 = PASS (all required fields present)
  1 = FAIL (one or more required fields missing)
  2 = FILE_ERROR (file not found or unreadable)

Pure-stdlib implementation (no PyYAML dependency) for minimal install footprint.
"""

import sys
import re
import os
import argparse
from datetime import date, datetime
from pathlib import Path


# Cutoff date: cycles dated on or after this need BINDING contract::
CUTOFF_DATE = date(2026, 5, 10)

# Required contract:: subfields (per meta/CYCLE_KICKOFF_TEMPLATE.md §1)
REQUIRED_L1_NATIVE_FIELDS = [
    "output_observable",      # MUST: jednostki fizyczne, NIE β_ppE
    "measurement_instrument", # MUST: realny pomiar
    "falsification_rule",     # MUST: decision rule pre-data
    "pre_registration_date",  # MUST: immutable timestamp
]

REQUIRED_L2_FIELDS = [
    "target_frameworks",
    "reduction_type",
    "failure_disposition",
]

REQUIRED_YAML_FIELDS = [
    "output_type",  # observable | projection | structural — per BINDING §2.2
]

REQUIRED_BODY_MARKERS = [
    "§0.4",  # Pre-flight methodology read confirmation section
    "PPN_AS_PROJECTION",  # checkbox item from §0.4
    "TGP_NATIVE_COMPUTATIONAL_PATTERNS",
    "M9_RESTRUCTURE_NOTE",
    "CYCLE_KICKOFF_TEMPLATE",
]


def extract_yaml_block(content: str) -> str | None:
    """Extract front matter YAML between leading --- delimiters."""
    match = re.match(r"^---\s*\n(.*?)\n---\s*\n", content, re.DOTALL)
    if match:
        return match.group(1)
    return None


def extract_field_from_yaml(yaml_text: str, field: str) -> str | None:
    """Extract top-level field value from YAML (no nesting support — uses indent heuristic)."""
    pattern = rf"^{re.escape(field)}\s*:\s*(.*?)$"
    match = re.search(pattern, yaml_text, re.MULTILINE)
    if match:
        value = match.group(1).strip()
        return value if value else None
    return None


def has_contract_block(yaml_text: str) -> bool:
    """Check if `contract:` block (with colon, top-level) is present."""
    return bool(re.search(r"^contract\s*:\s*$", yaml_text, re.MULTILINE))


def extract_contract_subfield(yaml_text: str, parent: str, subfield: str) -> str | None:
    """
    Extract nested subfield z contract: > L1_native: > <subfield>.
    Indent-aware (looks for indented field within parent block).
    """
    # Find parent block
    parent_pattern = rf"^(\s*){re.escape(parent)}\s*:\s*$"
    parent_match = re.search(parent_pattern, yaml_text, re.MULTILINE)
    if not parent_match:
        return None

    parent_indent_len = len(parent_match.group(1))
    parent_end = parent_match.end()
    rest = yaml_text[parent_end:]

    # Read lines until indent drops back to parent level or lower
    for line in rest.split("\n"):
        if not line.strip():
            continue
        line_indent_len = len(line) - len(line.lstrip())
        # Stripped first chars must be deeper than parent_indent_len
        if line_indent_len <= parent_indent_len and line.strip():
            # Left parent block
            break
        # Check if this line is our subfield
        stripped = line.lstrip()
        sub_match = re.match(rf"^{re.escape(subfield)}\s*:\s*(.*?)$", stripped)
        if sub_match:
            value = sub_match.group(1).strip()
            # Strip quotes if quoted
            if value.startswith('"') and value.endswith('"'):
                value = value[1:-1]
            if value.startswith("'") and value.endswith("'"):
                value = value[1:-1]
            return value if value else None
    return None


def parse_cycle_date(readme_path: Path, yaml_text: str | None) -> date | None:
    """
    Determine cycle's date from:
      1. YAML `date:` field
      2. Date suffix in folder name (op-NAME-YYYY-MM-DD)
      3. None if neither available
    """
    if yaml_text:
        date_str = extract_field_from_yaml(yaml_text, "date")
        if date_str:
            try:
                return datetime.strptime(date_str.strip(), "%Y-%m-%d").date()
            except ValueError:
                pass

    folder_name = readme_path.parent.name
    match = re.search(r"(\d{4}-\d{2}-\d{2})$", folder_name)
    if match:
        try:
            return datetime.strptime(match.group(1), "%Y-%m-%d").date()
        except ValueError:
            pass

    return None


def validate_readme(readme_path: Path, verbose: bool = True) -> tuple[bool, list[str]]:
    """
    Validate README.md against post-2026-05-10 BINDING kickoff template.

    Returns:
        (pass: bool, missing_fields: list[str])
    """
    try:
        content = readme_path.read_text(encoding="utf-8")
    except (FileNotFoundError, PermissionError) as e:
        return False, [f"FILE_ERROR: {e}"]

    yaml_text = extract_yaml_block(content)
    if yaml_text is None:
        return False, ["NO_YAML_FRONT_MATTER"]

    cycle_date = parse_cycle_date(readme_path, yaml_text)
    if cycle_date and cycle_date < CUTOFF_DATE:
        if verbose:
            print(f"  SKIP: cycle date {cycle_date} < cutoff {CUTOFF_DATE} (pre-BINDING)")
        return True, []

    missing = []

    # Check 1: contract:: block exists
    if not has_contract_block(yaml_text):
        missing.append("YAML.contract: (BINDING block per CYCLE_KICKOFF_TEMPLATE §1)")
        # If contract block missing, all subfields also missing — don't proliferate
    else:
        # Check 2: L1_native subfields
        for field in REQUIRED_L1_NATIVE_FIELDS:
            value = extract_contract_subfield(yaml_text, "L1_native", field)
            if not value or value in ('""', "''", "[]", "{}"):
                missing.append(f"YAML.contract.L1_native.{field}")

        # Check 3: L2_framework_reduction subfields (less strict — at least reduction_type)
        l2_reduction_type = extract_contract_subfield(
            yaml_text, "L2_framework_reduction", "reduction_type"
        )
        if not l2_reduction_type:
            missing.append("YAML.contract.L2_framework_reduction.reduction_type")

    # Check 4: top-level YAML fields
    for field in REQUIRED_YAML_FIELDS:
        # output_type can be nested under tgp_status:
        if not _has_field_anywhere(yaml_text, field):
            missing.append(f"YAML.{field} (or YAML.tgp_status.{field})")

    # Check 5: §0.4 Pre-flight methodology read confirmation section in body
    body = content[content.find("\n---\n", 4) + 5:] if "\n---\n" in content[4:] else content
    section_04_present = re.search(r"##\s*§0\.4\s*[—\-]", body)
    if not section_04_present:
        missing.append("BODY: ## §0.4 — Pre-flight methodology read confirmation (BINDING §2.6)")
    else:
        # Check that all 4 methodology refs are mentioned in §0.4 vicinity
        section_04_block = _extract_section(body, "§0.4")
        for ref in ["PPN_AS_PROJECTION", "TGP_NATIVE_COMPUTATIONAL_PATTERNS",
                    "M9_RESTRUCTURE_NOTE", "CYCLE_KICKOFF_TEMPLATE"]:
            if ref not in section_04_block:
                missing.append(f"BODY: §0.4 missing reference do meta/{ref}.md")

    return (len(missing) == 0, missing)


def _has_field_anywhere(yaml_text: str, field: str) -> bool:
    """Check if field appears anywhere in YAML (top-level or nested under tgp_status)."""
    return bool(re.search(rf"^\s*{re.escape(field)}\s*:\s*\S", yaml_text, re.MULTILINE))


def _extract_section(body: str, marker: str) -> str:
    """Extract block of text from marker line until next ## heading."""
    start = body.find(marker)
    if start == -1:
        return ""
    # Find next ## heading
    next_heading = body.find("\n## ", start + len(marker))
    if next_heading == -1:
        return body[start:]
    return body[start:next_heading]


def find_new_cycles(vault_root: Path) -> list[Path]:
    """Find all research/op-*-YYYY-MM-DD/README.md with date >= CUTOFF."""
    research_dir = vault_root / "research"
    if not research_dir.is_dir():
        return []
    candidates = []
    for cycle_dir in research_dir.iterdir():
        if not cycle_dir.is_dir():
            continue
        if not cycle_dir.name.startswith("op-"):
            continue
        match = re.search(r"(\d{4}-\d{2}-\d{2})$", cycle_dir.name)
        if not match:
            continue
        try:
            cycle_date = datetime.strptime(match.group(1), "%Y-%m-%d").date()
        except ValueError:
            continue
        if cycle_date >= CUTOFF_DATE:
            readme = cycle_dir / "README.md"
            if readme.is_file():
                candidates.append(readme)
    return sorted(candidates)


def main():
    parser = argparse.ArgumentParser(
        description="Validate cycle README.md against BINDING CYCLE_KICKOFF_TEMPLATE."
    )
    parser.add_argument("readme", nargs="?", help="Path to README.md to validate")
    parser.add_argument(
        "--all-new", action="store_true",
        help="Validate all research/op-*-YYYY-MM-DD/README.md with date >= 2026-05-10"
    )
    parser.add_argument(
        "--pre-commit", action="store_true",
        help="Pre-commit hook mode — read staged file paths from stdin"
    )
    parser.add_argument(
        "--quiet", action="store_true", help="Only print FAIL details + summary"
    )
    args = parser.parse_args()

    verbose = not args.quiet
    vault_root = Path(__file__).parent.parent  # tooling/../  = TGP_v1/

    targets = []

    if args.readme:
        targets.append(Path(args.readme).resolve())
    elif args.all_new:
        targets = find_new_cycles(vault_root)
        if verbose:
            print(f"Scanning {len(targets)} post-{CUTOFF_DATE} cycles...")
    elif args.pre_commit:
        # Read staged file paths from stdin (git pre-commit hook)
        for line in sys.stdin:
            path = line.strip()
            if not path:
                continue
            if not path.endswith("/README.md"):
                continue
            if "research/op-" not in path:
                continue
            full_path = vault_root / path
            if full_path.is_file():
                targets.append(full_path)
    else:
        parser.print_help()
        return 0

    if not targets:
        print("No targets to validate.")
        return 0

    fail_count = 0
    pass_count = 0
    skip_count = 0

    for readme_path in targets:
        rel_path = readme_path.relative_to(vault_root) if vault_root in readme_path.parents else readme_path
        if verbose:
            print(f"\n[{rel_path}]")

        passed, missing = validate_readme(readme_path, verbose=verbose)

        if not missing and passed:
            # Maybe skipped due to date
            try:
                content = readme_path.read_text(encoding="utf-8")
                yaml_text = extract_yaml_block(content)
                cycle_date = parse_cycle_date(readme_path, yaml_text)
                if cycle_date and cycle_date < CUTOFF_DATE:
                    skip_count += 1
                    continue
            except Exception:
                pass
            pass_count += 1
            if verbose:
                print("  PASS — all BINDING fields present")
        elif passed:
            skip_count += 1
        else:
            fail_count += 1
            print(f"  FAIL — {len(missing)} missing field(s):")
            for field in missing:
                print(f"    - {field}")

    if verbose or fail_count > 0:
        total = pass_count + fail_count + skip_count
        print(f"\nSummary: {pass_count} PASS / {fail_count} FAIL / {skip_count} SKIP "
              f"(pre-cutoff) of {total} total")

    return 1 if fail_count > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
