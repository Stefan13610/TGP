#!/usr/bin/env python3
r"""Fix 37 \ref orphans detected by DEPENDENCIES.md analysis.

Applies a carefully curated label remap. Also injects one missing
\label{ssec:T-Brannen} near the Brannen subsection heading.

Idempotent: running twice is a no-op (no more matches to replace).
"""
import re
from pathlib import Path

ROOT = Path('.')

# Each entry: (file, old_label, new_label, expected_count)
# `expected_count` = 0 means don't verify count (best effort).
FIXES = [
    # axioms/roznica_N0/dodatek0_aksjomatyka_roznicy.tex
    ('axioms/roznica_N0/dodatek0_aksjomatyka_roznicy.tex', 'app:R2', 'app:R2-qk-z3', 5),
    ('axioms/roznica_N0/dodatek0_aksjomatyka_roznicy.tex', 'ssec:aksjomaty', 'app:0-aksjomat', 1),
    ('axioms/roznica_N0/dodatek0_aksjomatyka_roznicy.tex', 'app:N', 'app:notacja', 1),
    ('axioms/roznica_N0/dodatek0_aksjomatyka_roznicy.tex', 'app:E-pi1', 'app:pi1-formal-proof', 1),

    # core/_meta_latex/status_map.tex
    ('core/_meta_latex/status_map.tex', 'sec:akcja-zunifikowana', 'ssec:unified-action', 1),
    ('core/_meta_latex/status_map.tex', 'prop:anomaly-Phi-robust', 'thm:anomaly-Phi-robust', 1),

    # core/formalizm/dodatekH_lancuch_wyprowadzen.tex
    ('core/formalizm/dodatekH_lancuch_wyprowadzen.tex', 'app:R-zero-mode', 'app:zero-mode-A4', 1),

    # core/sek02_pole/sek02_pole.tex
    ('core/sek02_pole/sek02_pole.tex', 'ssec:ppn', 'ssec:ppn-full', 1),
    ('core/sek02_pole/sek02_pole.tex', 'sec:ciemna-energia', 'sec:ciemna', 1),
    ('core/sek02_pole/sek02_pole.tex', 'ssec:metryka-z-substratu', 'ssec:metric-from-budget', 1),

    # core/sek04_stale/sek04_stale.tex
    ('core/sek04_stale/sek04_stale.tex', 'thm:metric-from-budget', 'thm:metric-from-substrate-full', 1),
    ('core/sek04_stale/sek04_stale.tex', 'sec:metryka-z-substratu', 'ssec:metric-from-budget', 1),

    # core/sek08_formalizm/sek08_formalizm.tex (8 occurrences)
    ('core/sek08_formalizm/sek08_formalizm.tex', 'sec:akcja-zunifikowana', 'ssec:unified-action', 8),

    # core/sek09_cechowanie/sek09_cechowanie.tex
    ('core/sek09_cechowanie/sek09_cechowanie.tex', 'app:T-koide-formal', 'app:T-koide', 1),

    # partial_proofs/defect_hierarchy/dodatekD2_defect_hierarchy_proof.tex
    ('partial_proofs/defect_hierarchy/dodatekD2_defect_hierarchy_proof.tex', 'app:O-u1', 'app:u1_formal', 1),

    # partial_proofs/koide_fp/dodatekT6_alpha3_subbreakthrough.tex
    ('partial_proofs/koide_fp/dodatekT6_alpha3_subbreakthrough.tex', 'app:T5', 'app:T5-c1', 2),
    ('partial_proofs/koide_fp/dodatekT6_alpha3_subbreakthrough.tex', 'eq:substrate-ODE', 'eq:P-substrate-ODE', 1),

    # partial_proofs/koide_fp/dodatekT_koide_atail_formal.tex
    # ssec:T-Brannen -> we'll ADD missing label, keep ref as-is (handled separately below)
    ('partial_proofs/koide_fp/dodatekT_koide_atail_formal.tex', 'eq:T-r-from-r21', 'eq:T-r-from-QK', 1),
    ('partial_proofs/koide_fp/dodatekT_koide_atail_formal.tex', 'rem:T-equipartition-status', 'rem:T-equipartition-motivation', 1),

    # partial_proofs/particle_sector/dodatekP4_sektor_czastek_closure.tex
    ('partial_proofs/particle_sector/dodatekP4_sektor_czastek_closure.tex', 'app:T3', 'app:T3-brannen', 1),

    # partial_proofs/quark_sector/dodatekX_quark_sector.tex
    ('partial_proofs/quark_sector/dodatekX_quark_sector.tex', 'ssec:X-confinement', 'app:X-confinement-bc', 1),

    # partial_proofs/superconductivity/dodatekP5_nadprzewodnictwo_closure.tex
    ('partial_proofs/superconductivity/dodatekP5_nadprzewodnictwo_closure.tex', 'app:B', 'app:substrat', 1),

    # partial_proofs/zero_mode/dodatekR_zero_mode_A4.tex
    # thm:R-virial-E2 (2x) -> point to existing prop:R-Elin-zero (virial identity in zero_mode)
    ('partial_proofs/zero_mode/dodatekR_zero_mode_A4.tex', 'thm:R-virial-E2', 'prop:R-Elin-zero', 2),
]

# Separate: add missing \label{ssec:T-Brannen} to subsection heading
# in dodatekT_koide_atail_formal.tex, right after the subsection line.
ADD_LABEL = {
    'partial_proofs/koide_fp/dodatekT_koide_atail_formal.tex': (
        # marker: the line RIGHT BEFORE which we insert the label
        r'\subsection{Parametryzacja Brannena --- mechanizm $\mathbb{Z}_3$}%',
        r'\label{ssec:T-Brannen}',
    ),
}


def replace_ref(text, old, new):
    r"""Replace \ref{old}, \eqref{old}, \pageref{old}, etc. with new label.

    Matches any \<refcmd>{old} where refcmd ends in 'ref' OR is exactly 'ref'.
    More precisely: matches the whole-brace form {old} only.
    """
    # Escape regex-significant chars in old
    pattern = re.compile(r'\{' + re.escape(old) + r'\}')
    new_text, n = pattern.subn('{' + new + '}', text)
    return new_text, n


def main():
    applied = 0
    skipped = 0
    errors = []
    for rel_path, old, new, expected in FIXES:
        p = ROOT / rel_path
        if not p.exists():
            errors.append(f'MISSING FILE: {rel_path}')
            continue
        text = p.read_text(encoding='utf-8')
        new_text, n = replace_ref(text, old, new)
        if n == 0:
            skipped += 1
            print(f'[skip] {rel_path}: `{old}` (0 occurrences — already fixed?)')
            continue
        if expected and n != expected:
            print(f'[warn] {rel_path}: `{old}` -> `{new}` expected {expected}, got {n}')
        else:
            print(f'[ok]   {rel_path}: `{old}` -> `{new}` ({n})')
        p.write_text(new_text, encoding='utf-8')
        applied += n

    # Add missing \label{ssec:T-Brannen}
    for rel_path, (marker, label) in ADD_LABEL.items():
        p = ROOT / rel_path
        text = p.read_text(encoding='utf-8')
        # Check if label already present
        if f'\\label{{{label.split("{")[1].rstrip("}")}}}' in text:
            # Simpler: check if "ssec:T-Brannen" is already a \label target
            pass
        if '\\label{ssec:T-Brannen}' in text:
            print(f'[skip] {rel_path}: label ssec:T-Brannen already present')
            continue
        if marker not in text:
            errors.append(f'MARKER NOT FOUND: {rel_path}: {marker}')
            continue
        new_text = text.replace(
            marker,
            marker + '\n' + label,
            1,
        )
        p.write_text(new_text, encoding='utf-8')
        print(f'[ok]   {rel_path}: inserted {label}')
        applied += 1

    print()
    print(f'=== SUMMARY ===')
    print(f'Applied edits: {applied}')
    print(f'Skipped (no match): {skipped}')
    if errors:
        print(f'ERRORS: {len(errors)}')
        for e in errors:
            print(f'  {e}')


if __name__ == '__main__':
    main()
