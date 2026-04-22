#!/usr/bin/env python3
r"""Fix 10 remaining wikilink orphans detected by DEPENDENCIES.md analysis.

Strategy:
  - Semantic remap where a reasonable existing target exists.
  - Convert aspirational wikilinks to plain-text code spans where no
    real target exists (so they no longer register as dead links but
    preserve authorial intent).

Idempotent: running twice is a no-op.
"""
from pathlib import Path

ROOT = Path('.')

# Each entry: (file, old_wikilink, new_wikilink_or_text)
# If new starts with "[[", it's still a wikilink; otherwise plain text.
FIXES = [
    # NUCLEAR_FROM_SOLITON_VERDICT.md: 4 aspirational refs to non-existent
    # section files -> map to actual existing appendices covering same topic.
    (
        'partial_proofs/nuclear_from_soliton/NUCLEAR_FROM_SOLITON_VERDICT.md',
        '[[sekcje_tex/sek05_leptony.tex]]',
        '[[dodatekF_hierarchia_mas.tex]]',
    ),
    (
        'partial_proofs/nuclear_from_soliton/NUCLEAR_FROM_SOLITON_VERDICT.md',
        '[[sekcje_tex/sek07_dyrak.tex]]',
        '[[dodatekE_kwantyzacja.tex]]',
    ),
    (
        'partial_proofs/nuclear_from_soliton/NUCLEAR_FROM_SOLITON_VERDICT.md',
        '[[sekcje_tex/sek09_topologia.tex]]',
        '[[dodatekU_su2_formalizacja.tex]]',
    ),
    (
        'partial_proofs/nuclear_from_soliton/NUCLEAR_FROM_SOLITON_VERDICT.md',
        '[[dodatekF_lepton_masses.tex]]',
        '[[dodatekF_hierarchia_mas.tex]]',
    ),

    # P6_closure.md (4 refs) -> P6_plan.md (the master plan for P6 closure).
    # Relative-path variant (3 files).
    (
        'research/casimir_mof/SUMMARY.md',
        '[[../superconductivity_closure/P6_closure.md]]',
        '[[../superconductivity_closure/P6_plan.md]]',
    ),
    (
        'research/liquid_viscosity/SUMMARY.md',
        '[[../superconductivity_closure/P6_closure.md]]',
        '[[../superconductivity_closure/P6_plan.md]]',
    ),
    (
        'research/muon_g_minus_2/SUMMARY.md',
        '[[../superconductivity_closure/P6_closure.md]]',
        '[[../superconductivity_closure/P6_plan.md]]',
    ),
    # Same-folder variant (1 file).
    (
        'research/superconductivity_closure/P7C_open_candidates.md',
        '[[P6_closure.md]]',
        '[[P6_plan.md]]',
    ),

    # desi_dark_energy/README.md: 2 aspirational refs -> convert to plain
    # inline text so they no longer render as broken wikilinks.
    (
        'research/desi_dark_energy/README.md',
        '[[TGP/TGP_v1/research/cluster_anomaly/]]',
        '`research/cluster_anomaly/` (planowane)',
    ),
    (
        'research/desi_dark_energy/README.md',
        '[[TGP/TGP_v1/research/galaxy_scaling/gs08b_ghost_free_verification.py]]',
        '`research/galaxy_scaling/gs08b_ghost_free_verification.py` (planowane)',
    ),
]


def main():
    applied = 0
    skipped = 0
    errors = []
    for rel_path, old, new in FIXES:
        p = ROOT / rel_path
        if not p.exists():
            errors.append(f'MISSING FILE: {rel_path}')
            continue
        text = p.read_text(encoding='utf-8')
        count = text.count(old)
        if count == 0:
            print(f'[skip] {rel_path}: `{old}` (0 occurrences)')
            skipped += 1
            continue
        new_text = text.replace(old, new)
        p.write_text(new_text, encoding='utf-8')
        print(f'[ok]   {rel_path}: `{old}` -> `{new}` ({count})')
        applied += count

    print()
    print(f'=== SUMMARY ===')
    print(f'Applied edits: {applied}')
    print(f'Skipped: {skipped}')
    if errors:
        for e in errors:
            print(f'ERROR: {e}')


if __name__ == '__main__':
    main()
