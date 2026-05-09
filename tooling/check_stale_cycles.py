#!/usr/bin/env python3
"""
check_stale_cycles.py - Read-only weekly report of stale-active cycles.

Per meta/CYCLE_LIFECYCLE.md (BINDING 2026-05-09): cykl z folder_status=active
ALE bez commita >30 dni jest kandydatem do paused. Skrypt LISTUJE kandydatow,
NIE modyfikuje YAML (decyzja user'a).

Usage:
    python tooling/check_stale_cycles.py             # report
    python tooling/check_stale_cycles.py --threshold 14   # custom days
    python tooling/check_stale_cycles.py --strict        # >14 days = stale

Exit codes:
    0 = no stale-active cycles found (clean)
    1 = stale-active cycles found (review needed)

Output: Markdown table sortowana po days_since_commit DESC.
"""

import argparse
import io
import re
import subprocess
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

# Force UTF-8 output (Windows console default cp1250 chokes on emoji/Polish chars).
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

VAULT_ROOT = Path(__file__).resolve().parent.parent
RESEARCH_DIR = VAULT_ROOT / "research"

DEFAULT_THRESHOLD_DAYS = 30

FOLDER_STATUS_RE = re.compile(
    r"^\s*folder_status:\s*[\"']?(?P<value>[\w\-]+)[\"']?\s*$",
    re.MULTILINE,
)

SKIP_FOLDERS = {"_archive", "_sandbox"}


def get_last_commit_date(folder: Path) -> "datetime | None":
    rel = folder.relative_to(VAULT_ROOT).as_posix()
    try:
        result = subprocess.run(
            ["git", "log", "-1", "--format=%cI", "--", rel],
            cwd=str(VAULT_ROOT),
            capture_output=True,
            text=True,
            timeout=10,
        )
        out = result.stdout.strip()
        if not out:
            return None
        return datetime.fromisoformat(out)
    except (subprocess.TimeoutExpired, ValueError, FileNotFoundError):
        return None


def parse_folder_status(readme: Path) -> "str | None":
    try:
        text = readme.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return None
    match = FOLDER_STATUS_RE.search(text)
    return match.group("value") if match else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--threshold",
        type=int,
        default=DEFAULT_THRESHOLD_DAYS,
        help=f"Days without commit to consider stale (default: {DEFAULT_THRESHOLD_DAYS})",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Use 14-day threshold (CYCLE_LIFECYCLE.md stale_threshold_days)",
    )
    args = parser.parse_args()

    threshold = 14 if args.strict else args.threshold
    now = datetime.now(timezone.utc)

    stale: list[tuple[str, int, str]] = []  # (folder, days_since, last_commit_iso)
    no_history: list[str] = []

    for child in sorted(RESEARCH_DIR.iterdir()):
        if not child.is_dir() or child.name in SKIP_FOLDERS:
            continue
        readme = child / "README.md"
        if not readme.exists():
            continue
        status = parse_folder_status(readme)
        if status != "active":
            continue
        last_commit = get_last_commit_date(child)
        if last_commit is None:
            no_history.append(child.name)
            continue
        age = now - last_commit.astimezone(timezone.utc)
        if age >= timedelta(days=threshold):
            stale.append((child.name, age.days, last_commit.isoformat()))

    # Output
    print(f"# Stale-active cycles report\n")
    print(f"- Date: {now.isoformat()}")
    print(f"- Threshold: >={threshold} days without commit")
    print(f"- Vault: {VAULT_ROOT}\n")

    if not stale and not no_history:
        print("**Clean.** No stale-active cycles found.")
        return 0

    if stale:
        print(f"## Stale-active ({len(stale)})\n")
        print("| Cykl | Dni od commita | Ostatni commit |")
        print("|---|---|---|")
        for folder, days, iso in sorted(stale, key=lambda x: -x[1]):
            print(f"| `{folder}` | {days} | {iso[:10]} |")
        print()

    if no_history:
        print(f"## Active bez historii git ({len(no_history)})\n")
        for folder in no_history:
            print(f"- `{folder}` — nigdy nie commitowane")
        print()

    print("## Rekomendacja")
    print()
    print("Per [[meta/CYCLE_LIFECYCLE.md]] §Stale-detection: powyższe kandydaty")
    print("powinny być zreklasyfikowane na `paused` (z spisanym blocker'em w README")
    print("§Status), chyba że są aktywnie używane.")
    print()
    print("Mass-reclassify: `python tooling/reclassify_cycles_2026-05-09.py --apply`")
    print("(ale skrypt operuje na całym Bucket A+B+C — sprawdź dry-run najpierw).")

    return 1 if (stale or no_history) else 0


if __name__ == "__main__":
    sys.exit(main())
