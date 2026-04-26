import re
import os

tex_path = os.path.join('..','..','..','tgp-core-paper','paper','tgp_core.tex')
tex_path = os.path.abspath(tex_path)
with open(tex_path,'r',encoding='utf-8') as f:
    s = f.read()

labels = sorted(set(re.findall(r'\\label\{([^}]+)\}', s)))
ref_pairs = re.findall(r'\\ref\{([^}]+)\}|\\eqref\{([^}]+)\}', s)
refs = sorted(set(a or b for a,b in ref_pairs))

orphan_refs = [r for r in refs if r not in labels]
unused_labels = [l for l in labels if l not in refs]

# File path citations
paths = sorted(set(re.findall(r'\\texttt\{([^}]+)\}', s)))
external_paths = [p for p in paths if '/' in p or '\\' in p]

print('=== Total labels:', len(labels))
print('=== Total refs:  ', len(refs))
print()
print('=== Orphan refs (referenced, no label) ===')
for r in orphan_refs:
    print('   ', r)
print()
print('=== Unused labels ===')
for l in unused_labels:
    print('   ', l)
print()
print('=== File-path citations in \\texttt{} ===')
for p in external_paths:
    print('   ', p)
