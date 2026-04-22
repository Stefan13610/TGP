#!/usr/bin/env python3
r"""Fix stale `scripts/plots/` references in tooling/scripts/*.py docstrings/prints.

Context: during the reorg, `scripts/` was moved to `tooling/scripts/`. The
plots/ subfolder followed, so the real location is `tooling/scripts/plots/`.
The actual savefig() calls resolve paths via __file__, so they keep working.
This script only cleans up the stale literal references in docstrings and
print statements so documentation matches reality.

Scope:
  - All .py files under tooling/scripts/ (NOT under _archiwum/).
  - All .md files under tooling/scripts/ (if any).

Skipped by design:
  - _archiwum/** (historical snapshot; do not touch)

Idempotent: running twice is a no-op.
"""
from pathlib import Path

ROOT = Path('.')
SEARCH_DIRS = [ROOT / 'tooling' / 'scripts']
OLD = 'scripts/plots/'
NEW = 'tooling/scripts/plots/'


def main():
    touched = 0
    hits = 0
    skipped = 0
    errors = []
    for base in SEARCH_DIRS:
        if not base.exists():
            errors.append(f'MISSING: {base}')
            continue
        for p in sorted(base.rglob('*')):
            if not p.is_file():
                continue
            if p.suffix not in {'.py', '.md'}:
                continue
            # Never touch archived content.
            if '_archiwum' in p.parts:
                continue
            try:
                text = p.read_text(encoding='utf-8')
            except UnicodeDecodeError:
                skipped += 1
                continue
            count = text.count(OLD)
            if count == 0:
                continue
            # Be conservative: don't match the NEW string recursively.
            # Since NEW contains OLD as substring, we need to be careful.
            # Strategy: temporarily protect NEW with a sentinel, then
            # replace OLD, then restore.
            SENTINEL = '\x00TGPNEW\x00'
            tmp = text.replace(NEW, SENTINEL)
            tmp2 = tmp.replace(OLD, NEW)
            new_text = tmp2.replace(SENTINEL, NEW)
            actual_changes = sum(
                1 for _ in range(1)
                if new_text != text
            )
            if new_text == text:
                continue
            # Count genuine replacements (OLD not preceded by 'tooling/').
            # This is approximate; the sentinel logic above guarantees correctness.
            p.write_text(new_text, encoding='utf-8')
            delta = text.count(OLD) - new_text.count(OLD)
            hits += delta
            touched += 1
            print(f'[ok] {p}: {delta} replacement(s)')

    print()
    print(f'=== SUMMARY ===')
    print(f'Files touched: {touched}')
    print(f'Replacements: {hits}')
    print(f'Binary/unreadable skipped: {skipped}')
    if errors:
        for e in errors:
            print(f'ERROR: {e}')


if __name__ == '__main__':
    main()
