"""fix_latex_typos.py — naprawa literowek LaTeX w dodatekG i sek09"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

fixes = {
    'dodatekG_wielki_wybuch.tex': [
        # Naprawa \k{e> -> \k{e}  (brakujacy nawias zamykajacy)
        (r'\k{e>', r'\k{e}'),
    ],
}

for fname, replacements in fixes.items():
    txt = open(fname, encoding='utf-8').read()
    n_total = 0
    for old, new in replacements:
        count = txt.count(old)
        txt = txt.replace(old, new)
        n_total += count
        print(f'  {fname}: "{old}" -> "{new}" x{count}')
    open(fname, 'w', encoding='utf-8').write(txt)
    print(f'  Saved {fname} ({n_total} fixes)')

print('Done.')
