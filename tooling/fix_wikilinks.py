#!/usr/bin/env python3
"""Bulk fix wikilinks after TGP_v1 reorganization.

For each wikilink [[path/to/file.ext]] (or [[path/to/file.ext|alias]]):
  - If basename matches a known moved file, replace full path with just basename.
  - Preserve alias if present.
  - Skip wikilinks without '/'.
  - Skip files inside _archiwum/ and research/nbody/_archiwum_docs/.
"""

import os
import re
import sys
from pathlib import Path

ROOT = Path(r"C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1")

# Build set of known basenames from the move map (all filenames that were moved)
MOVED_BASENAMES = set()

# Sekcje glowne
for name in [
    "sek00_summary.tex", "sek01_ontologia.tex", "sek02_pole.tex",
    "sek03_rezimy.tex", "sek04_stale.tex", "sek05_ciemna_energia.tex",
    "sek06_czarne_dziury.tex", "sek07_predykcje.tex",
    "sek07a_wymiar_wzmocniony.tex", "sek08_formalizm.tex",
    "sek08a_akcja_zunifikowana.tex", "sek08b_ghost_resolution.tex",
    "sek08c_metryka_z_substratu.tex", "sek09_cechowanie.tex",
    "sek10_N0_wyprowadzenie.tex",
]:
    MOVED_BASENAMES.add(name)

# Aksjomaty
for name in [
    "dodatek0_aksjomatyka_roznicy.tex", "dodatekA_notacja.tex",
    "slownik_formalizmu.tex", "dodatekB_substrat.tex",
]:
    MOVED_BASENAMES.add(name)

# Core formalizm
for name in [
    "dodatekE_kwantyzacja.tex", "dodatekE_pi1_formal.tex",
    "dodatekH_lancuch_wyprowadzen.tex", "dodatekL_formula_sumy.tex",
    "dodatekM_erg_stabilizacja.tex", "dodatekN_erg_renormalizacja.tex",
    "dodatekO_u1_formalizacja.tex", "dodatekQ_coarse_graining_formal.tex",
    "dodatekU_su2_formalizacja.tex", "dodatekV_su3_formalizacja.tex",
]:
    MOVED_BASENAMES.add(name)

# Meta LaTeX
for name in [
    "status_map.tex", "tabela_epistemiczna.tex",
    "nota_c_phi_reconciliation.tex",
]:
    MOVED_BASENAMES.add(name)

# Partial proofs
for name in [
    "dodatekC_ringdown.tex", "dodatekK_wkb_atail.tex",
    "dodatekD2_defect_hierarchy_proof.tex", "dodatekD_trojcialowe.tex",
    "dodatekY_nbody.tex", "dodatekF_hierarchia_mas.tex",
    "dodatekF_v2_wkb_numerics.tex", "dodatekJ_ogon_masy.tex",
    "dodatekJ2_sciezka9_formalizacja.tex", "dodatekG_wielki_wybuch.tex",
    "dodatekI_v2_potencjal.tex", "dodatekK2_chiralnosc.tex",
    "dodatekT_koide_atail_formal.tex", "dodatekT2_koide_fp_algebra.tex",
    "dodatekT3_brannen_geometry.tex", "dodatekT5_c1_deficit_excess.tex",
    "dodatekT6_alpha3_subbreakthrough.tex", "dodatekT4_bounce_topology.tex",
    "dodatekP_quantum_emergence.tex", "dodatekP4_sektor_czastek_closure.tex",
    "dodatekP5_nadprzewodnictwo_closure.tex",
    "dodatekP6_broader_sc_mechanisms.tex",
    "dodatekQ2_most_gamma_phi_lematy.tex",
    "dodatekW_agamma_phi0_update.tex", "dodatekR_zero_mode_A4.tex",
    "dodatekR2_qk_z3_dynamical.tex", "dodatekX_quark_sector.tex",
    "dodatekZ_alphaK_status.tex",
]:
    MOVED_BASENAMES.add(name)

# Additionally: enumerate all .md and .tex files currently in TGP_v1 (excluding archive)
# so we also simplify paths to .md files that were moved (FERMION_*, PLAN.md, etc.)
# Note: basename must be globally unique for Obsidian default resolver. We'll still
# replace if basename is unique within the vault OR if the path looks moved.

def build_file_index():
    """Map basename -> list of paths (relative to ROOT) for all relevant files,
    excluding _archiwum."""
    idx = {}
    exts = (".md", ".tex", ".py", ".txt", ".pdf", ".bib")
    for dirpath, dirnames, filenames in os.walk(ROOT):
        # skip archive dirs
        dirnames[:] = [d for d in dirnames if d != "_archiwum" and d != "_archiwum_docs"]
        for fn in filenames:
            if fn.endswith(exts):
                rel = Path(dirpath, fn).relative_to(ROOT)
                idx.setdefault(fn, []).append(rel)
    return idx

FILE_INDEX = build_file_index()

# Exclude dirs from walking for processing
EXCLUDE_DIRS = {"_archiwum", "_archiwum_docs"}

WIKILINK_RE = re.compile(r"\[\[([^\[\]]+?)\]\]")


def should_skip_file(path: Path) -> bool:
    parts = path.relative_to(ROOT).parts
    for p in parts:
        if p in EXCLUDE_DIRS:
            return True
    return False


def process_wikilink(inner: str) -> str:
    """Given content inside [[ ... ]], return new content or original if no change."""
    # Split alias
    if "|" in inner:
        target, alias = inner.split("|", 1)
    else:
        target, alias = inner, None

    # Preserve anchor/heading: [[file#anchor]]
    if "#" in target:
        path_part, anchor = target.split("#", 1)
        anchor = "#" + anchor
    else:
        path_part, anchor = target, ""

    # Skip if no slash -> already basename-only
    if "/" not in path_part and "\\" not in path_part:
        return inner

    # Normalize separators
    norm = path_part.replace("\\", "/")
    basename = norm.rsplit("/", 1)[-1]

    # Skip if points into _archiwum
    if "_archiwum" in norm:
        return inner

    # Skip if points outside TGP_v1 (e.g. ../../../something outside)
    # Heuristic: if the basename is not in FILE_INDEX, leave alone
    # But also include the original MOVED_BASENAMES set as "known moved"
    # Check uniqueness in file index
    locations = FILE_INDEX.get(basename, [])
    is_known_moved = basename in MOVED_BASENAMES
    # Only simplify if basename is globally unique (single location) OR it's a known moved .tex
    # For .tex moved files: MOVED_BASENAMES are all globally unique by design
    if is_known_moved and len(locations) <= 1:
        new_target = basename + anchor
        if alias is not None:
            return f"{new_target}|{alias}"
        return new_target
    if basename in FILE_INDEX and len(locations) == 1:
        new_target = basename + anchor
        if alias is not None:
            return f"{new_target}|{alias}"
        return new_target

    # Otherwise leave alone (ambiguous basename OR not in vault)
    return inner


DRY_RUN = "--dry-run" in sys.argv


def process_file(path: Path):
    try:
        content = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        try:
            content = path.read_text(encoding="utf-8-sig")
        except Exception as e:
            return 0, str(e), []

    changes = 0
    unchanged_with_slash = []  # wikilinks with / that were NOT changed

    def repl(m):
        nonlocal changes
        inner = m.group(1)
        new_inner = process_wikilink(inner)
        if new_inner != inner:
            changes += 1
            return f"[[{new_inner}]]"
        else:
            # Record unchanged wikilinks that still contain /
            target = inner.split("|", 1)[0]
            if "/" in target and "_archiwum" not in target:
                unchanged_with_slash.append(inner)
            return m.group(0)

    new_content = WIKILINK_RE.sub(repl, content)

    if changes > 0 and not DRY_RUN:
        path.write_text(new_content, encoding="utf-8")
    return changes, None, unchanged_with_slash


def main():
    total_changes = 0
    per_file_changes = {}
    errors = []
    unchanged_by_file = {}  # file -> list of unchanged-with-slash wikilinks

    skipped_archive = 0

    for dirpath, dirnames, filenames in os.walk(ROOT):
        # Track skipped archive dirs
        skipped_dirs = [d for d in dirnames if d in EXCLUDE_DIRS]
        for d in skipped_dirs:
            # count md files in there
            for sd, _, sfs in os.walk(Path(dirpath, d)):
                for f in sfs:
                    if f.endswith(".md"):
                        skipped_archive += 1
        dirnames[:] = [d for d in dirnames if d not in EXCLUDE_DIRS]
        for fn in filenames:
            if not fn.endswith(".md"):
                continue
            path = Path(dirpath, fn)
            if should_skip_file(path):
                continue
            changes, err, unchanged = process_file(path)
            if err:
                errors.append((str(path), err))
                continue
            if changes > 0:
                rel = path.relative_to(ROOT)
                per_file_changes[str(rel)] = changes
                total_changes += changes
            if unchanged:
                rel = path.relative_to(ROOT)
                unchanged_by_file[str(rel)] = unchanged

    # Print summary
    print(f"DRY_RUN: {DRY_RUN}")
    print(f"TOTAL_CHANGES: {total_changes}")
    print(f"FILES_EDITED: {len(per_file_changes)}")
    print(f"SKIPPED_ARCHIVE_MD: {skipped_archive}")
    print("\n=== EDITED FILES ===")
    for f, c in sorted(per_file_changes.items(), key=lambda x: -x[1]):
        print(f"  {c:4d}  {f}")

    if unchanged_by_file:
        print("\n=== UNCHANGED WIKILINKS WITH '/' (suspects / untouched) ===")
        for f, links in sorted(unchanged_by_file.items()):
            print(f"  {f}:")
            # Dedup
            seen = set()
            for l in links:
                if l in seen:
                    continue
                seen.add(l)
                print(f"      [[{l}]]")

    if errors:
        print("\nERRORS:")
        for f, e in errors:
            print(f"  {f}: {e}")


if __name__ == "__main__":
    main()
