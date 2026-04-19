import os, re
files = sorted([f for f in os.listdir('.') if f.startswith('dodatek') and f.endswith('.tex')])
print(f'Total appendix files: {len(files)}')
total_section = 0
total_section_star = 0
for f in files:
    with open(f, 'r', encoding='utf-8') as fh:
        c = fh.read()
    sec = len(re.findall(r'^\\section\{', c, flags=re.M))
    sstar = len(re.findall(r'^\\section\*', c, flags=re.M))
    total_section += sec
    total_section_star += sstar
    if sec > 0 or sstar > 0:
        marker = ' <-- STARRED' if sstar > 0 and sec == 0 else ''
        print(f'{f:50s}  section={sec}  section*={sstar}{marker}')
print()
print(f'SUM: section (counted) = {total_section}, section* (not counted) = {total_section_star}')
