#!/usr/bin/env python3
"""
reclassify_cycles_2026-05-09.py — mass-triage cykli research/op-* + research/{non-op}/

Diagnoza 2026-05-09: 80+ cykli z folder_status: active, realnie WIP ~5.
Skrypt aplikuje reguly z meta/CYCLE_LIFECYCLE.md:

  - WIP-5 (hardcoded list ponizej) -> active
  - Closure marker (Phase_FINAL_close.md / Phase6_absolute_binding.md /
    Phase_FINAL.md / CLOSED.md) -> closed-resolved
  - folder_status: active z commitem >=14 dni i bez closure -> paused
  - folder_status: needs-bridge / parking / closed-* / null -> NIE RUSZAJ
  - folder_status: research, research-active -> active jesli w WIP-5,
    inaczej paused

Dwa tryby:
  python reclassify_cycles_2026-05-09.py             # dry-run (domyslnie)
  python reclassify_cycles_2026-05-09.py --apply     # aplikuje zmiany

Output:
  - Tabelka per-cykl: path, old_status, new_status, reason
  - Summary counts per nowy status
  - W trybie --apply: dodatkowo licznik plikow zmodyfikowanych

Wymagania: Python 3.8+, brak external deps (re + subprocess + pathlib).
"""

import re
import subprocess
import sys
from datetime import datetime, timedelta
from pathlib import Path

# === KONFIGURACJA (per STATE.md 2026-05-09) ===

VAULT_ROOT = Path(__file__).resolve().parent.parent
RESEARCH_DIR = VAULT_ROOT / "research"

# WIP-5 (active po reclassification). Critical path (S07) + 4 inne.
WIP_ACTIVE = {
    "op-S07-alternative-f-psi-derivation-2026-05-09",  # critical path
    "op-FRW-radiation-era-varying-c-2026-05-06",
    "op-emergent-metric-from-interaction-2026-05-09",
    "op-MAG-anomalous-moment-2026-05-09",
    "op-Phi-decomposition-photon-2026-05-07",
}

# Closure markers — obecnosc dowolnego z tych plikow w folderze cyklu
# oznacza ze cykl jest zamkniety (Phase 6 ABSOLUTE BINDING gate PASS).
CLOSURE_MARKERS = {
    "Phase_FINAL_close.md",
    "Phase_FINAL.md",
    "Phase6_absolute_binding.md",
    "CLOSED.md",
}

# Statusy ktorych NIE ruszamy (terminal lub explicit user-set).
PRESERVE_STATUSES = {
    "paused",  # post-2026-05-09: jawnie zamrozone
    "needs-bridge",
    "parking",
    "closed-resolved",
    "closed-NULL",
    "closed-superseded",
    "closed-FALSIFIED",
    "archive",
    "sandbox",
    "review",
    "audit",  # cykle wewnatrz research/ oznaczone jako audit
}

# System foldery do pominiecia (nie sa cyklami badawczymi).
SKIP_FOLDERS = {"_archive", "_sandbox", "external_review_2026-04-25"}

STALE_THRESHOLD_DAYS = 14  # >=14 dni bez commita + brak closure -> paused

# Regex do matchowania linii folder_status: <value> w YAML frontmatter.
# Obsluguje:
#   - top-level "folder_status: active"
#   - nested "  folder_status: active" (pod tgp_status)
FOLDER_STATUS_RE = re.compile(
    r"^(?P<indent>\s*)folder_status:\s*[\"']?(?P<value>[\w\-]+)[\"']?\s*$",
    re.MULTILINE,
)


# === FUNKCJE ===

def get_last_commit_date(folder: Path) -> "datetime | None":
    """Zwraca date ostatniego commita w folderze, lub None jesli brak."""
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
        # %cI = ISO 8601 (np. 2026-05-09T17:34:21+02:00)
        return datetime.fromisoformat(out)
    except (subprocess.TimeoutExpired, ValueError, FileNotFoundError):
        return None


def has_closure_marker(folder: Path) -> "str | None":
    """Zwraca nazwe marker'a zamkniecia jesli obecny, inaczej None."""
    for marker in CLOSURE_MARKERS:
        if (folder / marker).exists():
            return marker
    return None


def parse_frontmatter(readme: Path) -> "tuple[str, str | None]":
    """Czyta README i zwraca (full_text, current_folder_status)."""
    try:
        text = readme.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as e:
        return "", None
    match = FOLDER_STATUS_RE.search(text)
    if match is None:
        return text, None
    return text, match.group("value")


def decide_new_status(
    folder_name: str,
    current: "str | None",
    closure: "str | None",
    last_commit: "datetime | None",
    now: datetime,
) -> "tuple[str | None, str]":
    """Zwraca (new_status, reason) lub (None, reason) jesli nie ruszamy."""
    # 1. Brak README albo brak folder_status -> nie ruszaj (manual review).
    if current is None:
        return None, "no folder_status field (manual review)"

    # 2. Status terminal / explicit-user -> nie ruszaj.
    if current in PRESERVE_STATUSES:
        return None, f"preserved status '{current}'"

    # 3. WIP-5 -> active (no-op jesli juz active).
    if folder_name in WIP_ACTIVE:
        if current == "active":
            return None, "in WIP-5 already active (no change)"
        return "active", "in WIP-5 (promote)"

    # 4. Closure marker -> closed-resolved.
    if closure:
        return "closed-resolved", f"closure marker '{closure}' present"

    # 5. folder_status active/research/research-active + brak closure -> paused.
    if current in ("active", "research", "research-active"):
        if last_commit is None:
            return "paused", "no git history + not in WIP-5"
        age = now - last_commit
        if age >= timedelta(days=STALE_THRESHOLD_DAYS):
            return "paused", f"stale {age.days}d, not in WIP-5"
        # Ostatni commit <14 dni ale nie w WIP-5 -> paused (decyzja user'a).
        return "paused", f"recent ({age.days}d) but not in WIP-5"

    # 6. Inny status -> nie ruszaj.
    return None, f"unknown status '{current}' (manual review)"


def apply_change(readme: Path, new_status: str) -> bool:
    """Modyfikuje folder_status: <value> w pliku README (in-place).
    Zachowuje indent i otaczajacy YAML."""
    text, current = parse_frontmatter(readme)
    if current is None:
        return False
    new_text = FOLDER_STATUS_RE.sub(
        lambda m: f"{m.group('indent')}folder_status: {new_status}",
        text,
        count=1,
    )
    if new_text == text:
        return False
    readme.write_text(new_text, encoding="utf-8")
    return True


def main() -> int:
    apply_mode = "--apply" in sys.argv
    now = datetime.now().astimezone()

    # Znajdz wszystkie folderowe README cykli badawczych.
    candidates = []
    for child in sorted(RESEARCH_DIR.iterdir()):
        if not child.is_dir():
            continue
        if child.name in SKIP_FOLDERS:
            continue
        readme = child / "README.md"
        if not readme.exists():
            candidates.append((child, None))  # bez README
            continue
        candidates.append((child, readme))

    rows = []  # (folder_name, old, new, reason, action)
    for folder, readme in candidates:
        if readme is None:
            rows.append((folder.name, "(no README)", "—", "skip — no README", "skip"))
            continue
        text, current = parse_frontmatter(readme)
        closure = has_closure_marker(folder)
        last_commit = get_last_commit_date(folder)
        new_status, reason = decide_new_status(
            folder.name, current, closure, last_commit, now
        )
        if new_status is None:
            rows.append((folder.name, current or "(none)", "—", reason, "skip"))
            continue
        if new_status == current:
            rows.append((folder.name, current, new_status, reason, "no-op"))
            continue
        # Realna zmiana.
        if apply_mode:
            ok = apply_change(readme, new_status)
            action = "APPLIED" if ok else "FAIL"
        else:
            action = "DRY-RUN"
        rows.append((folder.name, current or "(none)", new_status, reason, action))

    # === Output ===
    print(f"\n{'=' * 100}")
    print(f"reclassify_cycles_2026-05-09.py — {'APPLY MODE' if apply_mode else 'DRY-RUN'}")
    print(f"vault: {VAULT_ROOT}")
    print(f"now: {now.isoformat()}")
    print(f"{'=' * 100}\n")

    # Pelna tabela.
    fmt = "{:<60} {:<20} {:<20} {:<8} {}"
    print(fmt.format("FOLDER", "OLD", "NEW", "ACTION", "REASON"))
    print("-" * 140)
    for folder, old, new, reason, action in rows:
        print(fmt.format(folder[:60], old[:20], new[:20], action, reason))

    # Summary.
    print("\n" + "=" * 100)
    print("SUMMARY:")
    by_action: dict[str, int] = {}
    by_new_status: dict[str, int] = {}
    for _, _, new, _, action in rows:
        by_action[action] = by_action.get(action, 0) + 1
        if action in ("APPLIED", "DRY-RUN"):
            by_new_status[new] = by_new_status.get(new, 0) + 1

    print("\nPer action:")
    for action, n in sorted(by_action.items()):
        print(f"  {action:<10} {n}")

    print("\nPer new status (only changes):")
    for status, n in sorted(by_new_status.items()):
        print(f"  {status:<25} {n}")

    print(f"\nTotal candidates: {len(rows)}")
    print(f"Mode: {'APPLY' if apply_mode else 'DRY-RUN (use --apply to commit)'}")
    print("=" * 100 + "\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
