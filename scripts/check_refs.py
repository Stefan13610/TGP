"""Check for unresolved LaTeX cross-references in sek08_formalizm.tex."""
import re, glob, sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

labels = set()
for f in glob.glob('*.tex'):
    with open(f, encoding='utf-8') as fh:
        text = fh.read()
        pat = re.compile(r'\\label\{([^}]+)\}')
        for m in pat.finditer(text):
            labels.add(m.group(1))

refs = set()
with open('sek08_formalizm.tex', encoding='utf-8') as fh:
    content = fh.read()
    pat2 = re.compile(r'\\(?:eq)?ref\{([^}]+)\}')
    for m in pat2.finditer(content):
        refs.add(m.group(1))

unresolved = sorted(refs - labels)
print(f'Total labels across all .tex: {len(labels)}')
print(f'Total refs in sek08: {len(refs)}')
print(f'Unresolved: {len(unresolved)}')
for r in unresolved:
    print(f'  MISSING: {r}')

if not unresolved:
    print('\n  ALL REFERENCES RESOLVED')

# Check duplicate labels
from collections import Counter
all_labels = []
for f in glob.glob('*.tex'):
    with open(f, encoding='utf-8') as fh:
        pat3 = re.compile(r'\\label\{([^}]+)\}')
        for m in pat3.finditer(fh.read()):
            all_labels.append((m.group(1), f))

c = Counter(l for l, _ in all_labels)
dups = {k: v for k, v in c.items() if v > 1}
print(f'\n--- DUPLICATE LABELS ---')
if dups:
    print(f'Found {len(dups)} duplicates:')
    for k, v in sorted(dups.items()):
        files = [f for l, f in all_labels if l == k]
        print(f'  {k} ({v}x): {files}')
else:
    print('  No duplicate labels found')
