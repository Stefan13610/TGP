#!/usr/bin/env python3
"""
check_status_drift.py - Read-only detection of status <-> folder_status mismatch.

Wiele cykli ma DWA pola w YAML README:
- status: <text>          (human-readable, np. "PHASE0_PHASE1_IN_PROGRESS",
                           "EARLY_HALT_2026-05-09", "🔒 CLOSED")
- folder_status: <enum>   (operacyjny, np. "active", "paused", "closed-resolved")

Diagnoza 2026-05-09: po cascade closure i mass-triage tekstowy 'status' field
nie zostal zaktualizowany w wielu README -> mowi co innego niz folder_status.

Skrypt LISTUJE niespójności (NIE modyfikuje, bo text status jest semantically
varied — automatyczna konwersja jest niebezpieczna).

Heurystyki dopasowania (status text -> oczekiwany folder_status):
- contains 'CLOSED', '🔒', 'EXECUTED', 'COMPLETE', 'FINAL'  -> closed-*
- contains 'EARLY_HALT', 'NULL'                            -> closed-NULL
- contains 'ACTIVE', 'OPEN', 'IN_PROGRESS', 'SETUP'        -> active
- contains 'PLACEHOLDER', 'INFORMAL', 'CONCEPTUAL'         -> parking
- contains 'PAUSED', 'BLOCKED'                             -> paused
- contains 'BRIDGE', 'PENDING', 'AWAITS'                   -> needs-bridge
- inne                                                     -> unknown (skip)

Jezeli text status sugeruje X ale folder_status to Y -> mismatch raport.

Usage:
    python tooling/check_status_drift.py             # report all mismatches
    python tooling/check_status_drift.py --only-paused   # filter
    python tooling/check_status_drift.py --strict       # exit 1 if any drift

Output: Markdown table.
"""

import argparse
import io
import re
import sys
from pathlib import Path

# Force UTF-8 output (Windows console default cp1250 chokes on emoji/Polish chars).
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
else:  # Python <3.7 fallback
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

VAULT_ROOT = Path(__file__).resolve().parent.parent
RESEARCH_DIR = VAULT_ROOT / "research"

SKIP_FOLDERS = {"_archive", "_sandbox"}

# Top-level i nested 'status: <value>' field.
STATUS_TEXT_RE = re.compile(
    r"^status:\s*(?P<value>.+?)\s*$",
    re.MULTILINE,
)
FOLDER_STATUS_RE = re.compile(
    r"^\s*folder_status:\s*[\"']?(?P<value>[\w\-]+)[\"']?\s*$",
    re.MULTILINE,
)


def classify_text_status(text_status: str) -> "str | None":
    """Heurystyczne mapowanie text status -> oczekiwany folder_status."""
    if not text_status:
        return None
    upper = text_status.upper()

    # Sprawdz w okreslonej kolejnosci (pierwszy match wygrywa).
    if any(k in upper for k in ("EARLY_HALT", "NULL_CLOSED", "STAGE_1_NULL")):
        return "closed-NULL"
    if any(k in upper for k in ("CLOSED", "🔒", "EXECUTED", "COMPLETE", "FINAL")):
        return "closed-resolved"
    if "FALSIFIED" in upper:
        return "closed-FALSIFIED"
    if any(k in upper for k in ("PLACEHOLDER", "INFORMAL", "CONCEPTUAL")):
        return "parking"
    if any(k in upper for k in ("PAUSED", "BLOCKED")):
        return "paused"
    if any(k in upper for k in ("BRIDGE", "PENDING", "AWAITS")):
        return "needs-bridge"
    if any(k in upper for k in ("ACTIVE", "OPEN", "IN_PROGRESS", "SETUP", "PHASE")):
        return "active"
    return None  # nieznany pattern -> nie raportuj


def parse_readme(readme: Path) -> "tuple[str | None, str | None]":
    """Returns (text_status, folder_status)."""
    try:
        text = readme.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return None, None

    # Tylko z YAML frontmatter (pomiedzy --- i ---).
    if not text.startswith("---"):
        return None, None
    end = text.find("\n---\n", 4)
    if end < 0:
        return None, None
    yaml_block = text[4:end]

    text_status = None
    m1 = STATUS_TEXT_RE.search(yaml_block)
    if m1:
        text_status = m1.group("value").strip().strip("\"'")

    folder_status = None
    m2 = FOLDER_STATUS_RE.search(yaml_block)
    if m2:
        folder_status = m2.group("value")

    return text_status, folder_status


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only-paused", action="store_true",
                        help="Filter to cycles with folder_status=paused only")
    parser.add_argument("--strict", action="store_true",
                        help="Exit 1 if any drift detected")
    args = parser.parse_args()

    drifts: list[tuple[str, str, str, str]] = []  # (folder, text, folder_st, expected)
    no_text: list[str] = []

    for child in sorted(RESEARCH_DIR.iterdir()):
        if not child.is_dir() or child.name in SKIP_FOLDERS:
            continue
        readme = child / "README.md"
        if not readme.exists():
            continue
        text_st, folder_st = parse_readme(readme)
        if folder_st is None:
            continue
        if args.only_paused and folder_st != "paused":
            continue
        if not text_st:
            no_text.append(child.name)
            continue

        expected = classify_text_status(text_st)
        if expected is None:
            continue  # text status zbyt niejasny do klasyfikacji
        if expected != folder_st:
            drifts.append((child.name, text_st[:50], folder_st, expected))

    # Output
    print(f"# Status drift report\n")
    print(f"- Vault: {VAULT_ROOT}")
    print(f"- Filter: {'paused only' if args.only_paused else 'all'}\n")

    if not drifts and not no_text:
        print("**Clean.** No drift detected.")
        return 0

    if drifts:
        print(f"## Drift ({len(drifts)})\n")
        print("| Cykl | Text status (50ch) | folder_status | Sugerowane (z text) |")
        print("|---|---|---|---|")
        for folder, text, fs, expected in sorted(drifts):
            print(f"| `{folder}` | `{text}` | `{fs}` | `{expected}` |")
        print()

    if no_text:
        print(f"## Brak text status field ({len(no_text)})\n")
        for f in no_text:
            print(f"- `{f}`")
        print()

    print("## Notatka\n")
    print("To raport heurystyczny — nie wszystkie 'driftów' wymagaj fixa.")
    print("Sugerowana akcja: dla każdego rzędu zdecydować czy tekst status")
    print("warto zaktualizować, czy folder_status wymaga zmiany. Skrypt nie")
    print("modyfikuje plików — manual review.")

    return 1 if (args.strict and (drifts or no_text)) else 0


if __name__ == "__main__":
    sys.exit(main())
