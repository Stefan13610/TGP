import re, os
os.chdir(r'C:/Users/Mateusz/Documents/ObsydnianMain')
with open('TGP/tgp-core-paper/KNOWN_ISSUES.md','r',encoding='utf-8') as f:
    s = f.read()

# Build index: filename -> list of full paths, plus a flat set for suffix matching.
index = {}
all_filenames = set()
for root, dirs, files in os.walk('TGP'):
    # skip heavy dirs
    if '.git' in root or '__pycache__' in root: continue
    for fn in files:
        index.setdefault(fn, []).append(os.path.join(root, fn).replace('\\','/'))
        all_filenames.add(fn)

raw_bt = sorted(set(re.findall(r'`([^`]+\.(?:md|py|tex|txt|json))`', s)))

# Normalize LaTeX-escape leaks: KNOWN_ISSUES.md quotes the .tex source
# verbatim (e.g. "research/op6/m3a\_block\_rg\_1d.py" inside backticks),
# so the captured string contains LaTeX backslash-underscore. The actual
# filesystem name uses plain underscore. Same for \# (rare but cheap).
def normalize(L):
    return L.replace(r'\_', '_').replace(r'\#', '#')

bt = sorted(set(normalize(L) for L in raw_bt))

def resolve(L):
    L = L.strip().lstrip('/')
    if any(c in L for c in ['{','}','*']):
        return True, 'glob-skipped'
    # Direct candidates first
    cands = [L, 'TGP/'+L, 'TGP/TGP_v1/'+L, 'TGP/tgp-core-paper/'+L]
    for c in cands:
        if os.path.exists(c):
            return True, c
    # Try by basename via index
    base = os.path.basename(L)
    if base in index:
        return True, 'BARE: '+'; '.join(index[base])
    # Suffix-fragment fallback: catches continuation-backtick artefacts
    # such as "`pre_results.md` + `.txt` + `_results.md`" where the third
    # backtick captures only the suffix `_results.md`. If any indexed
    # filename ends with the captured fragment, count as resolved-fragment.
    if '/' not in L:
        suffix_matches = [fn for fn in all_filenames if fn.endswith(L) and fn != L]
        if suffix_matches:
            return True, 'SUFFIX-FRAGMENT: ' + str(len(suffix_matches)) + ' file(s) end with this fragment'
    return False, None

really_missing = []
ambiguous_bare = []
suffix_fragments = []
for L in bt:
    found, where = resolve(L)
    if not found:
        really_missing.append(L)
    elif isinstance(where, str):
        if where.startswith('BARE'):
            ambiguous_bare.append((L, where))
        elif where.startswith('SUFFIX-FRAGMENT'):
            suffix_fragments.append((L, where))

print('== Truly missing files referenced in KNOWN_ISSUES.md ==')
for m in really_missing:
    print('   MISSING:', m)
print()
print('== Suffix-fragment artefacts (continuation-backtick captures; not real refs) ==')
for L, where in suffix_fragments:
    print('   ', L, '->', where)
print()
print('== Bare filename references (resolves but no path prefix) ==')
for L, where in ambiguous_bare[:20]:
    print('   ', L, '->', where[:120])
print(f'   ... (total {len(ambiguous_bare)} bare references)')
