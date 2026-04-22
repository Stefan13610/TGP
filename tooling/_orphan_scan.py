#!/usr/bin/env python3
r"""One-shot orphan scan: .tex files not referenced by any \input{} in-tree.

Excludes standalone papers (papers_external/, tgp_core/letter/companion),
main.tex, and _archiwum/. research/nbody/ handled separately.
"""
import re
import os
from pathlib import Path

ROOT = Path('.')
INPUT_RE = re.compile(r'\\input\s*\{([^}]+)\}')

all_tex = set()
for p in ROOT.rglob('*.tex'):
    s = str(p).replace(os.sep, '/')
    if '_archiwum' in s:
        continue
    all_tex.add(s)

referenced = set()
for p in ROOT.rglob('*.tex'):
    try:
        text = p.read_text(encoding='utf-8', errors='ignore')
    except Exception:
        continue
    for m in INPUT_RE.finditer(text):
        tgt = m.group(1).strip()
        if not tgt.endswith('.tex'):
            tgt += '.tex'
        referenced.add(tgt.replace(os.sep, '/'))


def basename(p):
    return p.rsplit('/', 1)[-1]


# Skip-list: standalone documents
SKIP_BASENAMES = {'main.tex', 'tgp_core.tex', 'tgp_letter.tex', 'tgp_companion.tex'}

unref = []
for t in sorted(all_tex):
    bn = basename(t)
    if bn in SKIP_BASENAMES:
        continue
    # papers_external/* are standalone documents with their own main files
    if t.startswith('papers_external/'):
        continue
    # research/nbody — many files input each other; skip in this scan
    if t.startswith('research/nbody/'):
        continue
    is_ref = any(
        ref == t or ref.endswith('/' + bn) or ref == bn for ref in referenced
    )
    if not is_ref:
        unref.append(t)

print(f'Total in-scope .tex (excl. _archiwum/papers_external/nbody/masters): '
      f'{sum(1 for t in all_tex if basename(t) not in SKIP_BASENAMES and not t.startswith(("papers_external/", "research/nbody/")))}')
print(f'Unreferenced candidates: {len(unref)}')
for u in unref:
    print(f'  {u}')
