#!/usr/bin/env python3
r"""Analyze \ref orphans from DEPENDENCIES.md.

For each orphan label, search the corpus (.tex files) for \label{that-label}
to see if it exists but wasn't indexed (e.g., in _archiwum/) or truly missing.
Also reports similar-label matches (prefix suggestions) for fixable candidates.
"""
import re
from pathlib import Path
from collections import defaultdict

ROOT = Path('.')

# Parse orphan list from DEPENDENCIES.md
deps = (ROOT / 'DEPENDENCIES.md').read_text(encoding='utf-8')
in_ref = False
orphans = []  # list of (source_file, label)
for line in deps.splitlines():
    if line.startswith('### \\ref'):
        in_ref = True
        continue
    if line.startswith('### wikilinks'):
        in_ref = False
        continue
    if in_ref and line.startswith('- '):
        # Format: - `source` -> `label`
        m = re.match(r'- `([^`]+)` -> `([^`]+)`', line)
        if m:
            orphans.append((m.group(1), m.group(2)))

# Build complete \label{} index across ALL .tex files including _archiwum
label_index = defaultdict(list)  # label -> list of files
LABEL_RE = re.compile(r'\\label\s*\{([^}]+)\}')
for p in ROOT.rglob('*.tex'):
    text = p.read_text(encoding='utf-8', errors='ignore')
    for m in LABEL_RE.finditer(text):
        label_index[m.group(1)].append(str(p).replace('\\', '/'))

# Classify orphans
print(f'Orphan count: {len(orphans)}')
print(f'Total labels in index: {len(label_index)}')
print()

unique_labels = sorted(set(lbl for _, lbl in orphans))
print(f'Unique orphan labels: {len(unique_labels)}')
print()

print('=== Classification ===')
for lbl in unique_labels:
    count = sum(1 for _, l in orphans if l == lbl)
    sources = sorted(set(s for s, l in orphans if l == lbl))
    if lbl in label_index:
        # Label exists — must have failed to index for some reason
        files = label_index[lbl]
        print(f'\n[FOUND_IN_INDEX] `{lbl}` (used {count}x)')
        for f in files:
            print(f'  defined_in: {f}')
        for s in sources:
            print(f'  referenced_in: {s}')
    else:
        # Look for prefix or suffix matches
        candidates = [k for k in label_index if lbl in k or k in lbl]
        # Substring loose match
        loose = [k for k in label_index
                 if (lbl.lower() in k.lower() or k.lower() in lbl.lower())
                 and k not in candidates]
        print(f'\n[MISSING] `{lbl}` (used {count}x)')
        for s in sources:
            print(f'  from: {s}')
        if candidates:
            print(f'  similar_in_index: {candidates[:5]}')
        elif loose:
            print(f'  loose_match: {loose[:5]}')
