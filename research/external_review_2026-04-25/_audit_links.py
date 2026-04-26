import re, os
os.chdir(r'C:/Users/Mateusz/Documents/ObsydnianMain')
with open('TGP/tgp-core-paper/KNOWN_ISSUES.md','r',encoding='utf-8') as f:
    s = f.read()

# Build index: filename -> list of full paths
index = {}
for root, dirs, files in os.walk('TGP'):
    # skip heavy dirs
    if '.git' in root or '__pycache__' in root: continue
    for fn in files:
        index.setdefault(fn, []).append(os.path.join(root, fn).replace('\\','/'))

bt = sorted(set(re.findall(r'`([^`]+\.(?:md|py|tex|txt|json))`', s)))

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
    return False, None

really_missing = []
ambiguous_bare = []
for L in bt:
    found, where = resolve(L)
    if not found:
        really_missing.append(L)
    elif isinstance(where, str) and where.startswith('BARE'):
        ambiguous_bare.append((L, where))

print('== Truly missing files referenced in KNOWN_ISSUES.md ==')
for m in really_missing:
    print('   MISSING:', m)
print()
print('== Bare filename references (resolves but no path prefix) ==')
for L, where in ambiguous_bare[:20]:
    print('   ', L, '->', where[:120])
print(f'   ... (total {len(ambiguous_bare)} bare references)')
