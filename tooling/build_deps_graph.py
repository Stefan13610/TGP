#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
build_deps_graph.py - Build folder-level dependency graph for TGP/TGP_v1.

Extracts four kinds of cross-folder dependencies:

    1. \input{path}        - LaTeX file inclusion (strong)
    2. \ref/\eqref/\autoref/\cref{label} - LaTeX cross-references
       (resolved via a label->file map built from all \label{X})
    3. \cite{key}          - bibliographic citations (counted per folder only)
    4. [[wikilink]]        - Obsidian-style wikilinks in .md files
                             (resolved via basename lookup)

Produces two files in the project root (TGP/TGP_v1/):

    DEPENDENCIES.md          - forward graph ("what does this folder load?")
    DEPENDENCIES_REVERSE.md  - reverse graph ("who loads this folder?")

Run from TGP/TGP_v1/ (or anywhere - script autodetects its root).
"""
from __future__ import annotations

import os
import re
import sys
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

# Folders that are completely ignored when collecting files.
# Note: "scripts" and "plots" were previously skipped at the top level; now
# they can exist legitimately (e.g. tooling/scripts/) so we no longer drop
# them by name. Generated / archive dirs stay excluded.
SKIP_DIR_NAMES = {"_archiwum", "_archiwum_docs", ".git", "__pycache__"}

# \input{preamble} is infrastructure noise -> drop it.
SKIP_INPUT_TARGETS = {"preamble", "tgp_main", "tgp_companionNotes",
                      "tgp_letterNotes"}

# Top-level folders the user wants surfaced in the "coarse" view.
COARSE_TOP_LEVEL = [
    "axioms", "core", "core/_meta_latex", "core/formalizm",
    "partial_proofs", "research", "papers_external", "meta", "tooling",
]

# Regexes -------------------------------------------------------------------

RE_INPUT       = re.compile(r"\\input\s*\{([^}]+)\}")
RE_INCLUDE     = re.compile(r"\\include\s*\{([^}]+)\}")
RE_LABEL       = re.compile(r"\\label\s*\{([^}]+)\}")
RE_REF         = re.compile(r"\\(?:ref|eqref|autoref|cref|Cref)\s*\{([^}]+)\}")
RE_CITE        = re.compile(r"\\cite[a-zA-Z]*\s*(?:\[[^\]]*\])?\s*\{([^}]+)\}")
# wikilinks: [[target]] or [[target|alias]] or [[target#heading]]
RE_WIKILINK    = re.compile(r"\[\[([^\]\|#]+)(?:#[^\]\|]+)?(?:\|[^\]]+)?\]\]")
# %-comments (rough; strips from first unescaped % to EOL)
RE_COMMENT     = re.compile(r"(?<!\\)%.*$", re.MULTILINE)

# ----------------------------------------------------------------------------
# Path / folder helpers
# ----------------------------------------------------------------------------

def find_root() -> Path:
    """Locate TGP/TGP_v1 root from script location or cwd."""
    here = Path(__file__).resolve()
    for p in [here.parent.parent, Path.cwd(), Path.cwd() / "TGP" / "TGP_v1"]:
        if (p / "main.tex").exists() and (p / "tgp_main.bib").exists():
            return p.resolve()
    raise SystemExit("Cannot locate TGP/TGP_v1 root (main.tex + tgp_main.bib).")


def is_skipped(path: Path, root: Path) -> bool:
    """True if path is inside any SKIP_DIR_NAMES directory under root."""
    try:
        rel = path.relative_to(root)
    except ValueError:
        return True
    parts = set(rel.parts)
    return bool(parts & SKIP_DIR_NAMES)


def fine_folder_of(file_rel: Path) -> str:
    """Return the immediate parent folder (relative, POSIX) for a file.

    Root-level files return '<root>'.
    """
    parent = file_rel.parent
    if str(parent) in ("", "."):
        return "<root>"
    return parent.as_posix()


def coarse_folder_of(fine: str) -> str:
    """Map a fine folder to one of COARSE_TOP_LEVEL (or '<root>')."""
    if fine == "<root>":
        return "<root>"
    # Check exact matches first (core/_meta_latex, core/formalizm).
    for top in COARSE_TOP_LEVEL:
        if fine == top or fine.startswith(top + "/"):
            # We want longest-prefix match; COARSE_TOP_LEVEL is ordered
            # with two-segment entries AFTER 'core', so keep scanning.
            pass
    # Longest-prefix match against COARSE_TOP_LEVEL.
    best = ""
    for top in COARSE_TOP_LEVEL:
        if (fine == top or fine.startswith(top + "/")) and len(top) > len(best):
            best = top
    if best:
        return best
    # First path segment fallback.
    return fine.split("/", 1)[0]


# ----------------------------------------------------------------------------
# File collection
# ----------------------------------------------------------------------------

SKIP_FILES = {"DEPENDENCIES.md", "DEPENDENCIES_REVERSE.md"}


# Extensions of files that may legitimately appear as wikilink targets
# in Obsidian. Indexed for resolution but not parsed for dependencies.
ASSET_EXTS = {".py", ".json", ".txt", ".csv", ".yaml", ".yml",
              ".pdf", ".png", ".jpg", ".jpeg", ".svg", ".gif", ".webp",
              ".ipynb", ".c", ".h", ".cpp", ".bib"}


def collect_files(root: Path):
    """Walk root collecting .tex, .md and asset files, honoring SKIP_DIR_NAMES.

    Returns (tex_files, md_files, asset_files).
    """
    tex_files: list[Path] = []
    md_files: list[Path] = []
    asset_files: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root):
        # Prune in-place.
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIR_NAMES]
        p_dir = Path(dirpath)
        for fn in filenames:
            if fn in SKIP_FILES:
                continue
            fp = p_dir / fn
            ext = Path(fn).suffix.lower()
            if ext == ".tex":
                tex_files.append(fp)
            elif ext == ".md":
                md_files.append(fp)
            elif ext in ASSET_EXTS:
                asset_files.append(fp)
    return sorted(tex_files), sorted(md_files), sorted(asset_files)


def read_text(p: Path) -> str:
    for enc in ("utf-8", "utf-8-sig", "cp1250", "latin-1"):
        try:
            return p.read_text(encoding=enc)
        except UnicodeDecodeError:
            continue
    return p.read_bytes().decode("utf-8", errors="replace")


# ----------------------------------------------------------------------------
# LaTeX parsing
# ----------------------------------------------------------------------------

def strip_comments(tex: str) -> str:
    return RE_COMMENT.sub("", tex)


def resolve_input_target(target: str, from_file: Path, root: Path) -> Path | None:
    r"""Resolve a \input{target} to an actual file under root, if possible.

    \input accepts paths relative to the master document. TGP uses paths
    relative to the TGP/TGP_v1 root (e.g. 'core/sek02_pole/sek02_pole').
    It may also be relative to the including file's folder.
    """
    target = target.strip()
    if not target:
        return None
    candidates = []
    # Strip any trailing .tex the author may have typed.
    bare = target[:-4] if target.endswith(".tex") else target

    # Root-relative.
    candidates.append(root / (bare + ".tex"))
    candidates.append(root / bare)
    # Include-file relative.
    candidates.append(from_file.parent / (bare + ".tex"))
    candidates.append(from_file.parent / bare)

    for c in candidates:
        if c.is_file():
            return c.resolve()
    return None


def build_label_map(tex_files: list[Path], root: Path) -> dict[str, Path]:
    r"""Map \label{X} -> owning file."""
    label_to_file: dict[str, Path] = {}
    duplicates: list[tuple[str, Path, Path]] = []
    for fp in tex_files:
        txt = strip_comments(read_text(fp))
        for m in RE_LABEL.finditer(txt):
            label = m.group(1).strip()
            if label in label_to_file:
                duplicates.append((label, label_to_file[label], fp))
            else:
                label_to_file[label] = fp.resolve()
    return label_to_file


# ----------------------------------------------------------------------------
# Wikilink resolver
# ----------------------------------------------------------------------------

def build_basename_index(root: Path, tex_files, md_files,
                         asset_files=None) -> dict[str, list[Path]]:
    """Index files by basename (case-insensitive) for wikilink lookup.

    Both extensionless stems (e.g. "sek02_pole") and full filenames
    (e.g. "sek02_pole.tex", "pairwise.py") are indexed — Obsidian
    wikilinks use whichever form the author wrote.
    """
    idx: dict[str, list[Path]] = defaultdict(list)
    files = list(tex_files) + list(md_files) + list(asset_files or [])
    for fp in files:
        resolved = fp.resolve()
        idx[fp.stem.lower()].append(resolved)
        idx[fp.name.lower()].append(resolved)
    return idx


def resolve_wikilink(target: str, basename_idx, root: Path | None = None) -> Path | None:
    target = target.strip()
    if not target:
        return None
    # Folder-like wikilinks (trailing /) resolve to the folder itself
    # if it exists. Obsidian treats [[folder/]] as a link to the folder.
    if root is not None and target.endswith("/"):
        folder_rel = target.rstrip("/")
        # Handle absolute-from-vault paths like "TGP/TGP_v1/research/foo"
        for prefix in ("TGP/TGP_v1/", ""):
            if folder_rel.startswith(prefix):
                cand = root / folder_rel[len(prefix):]
                if cand.is_dir():
                    return cand.resolve()
    # Obsidian supports both relative and bare wikilinks. Try the exact
    # key first, then progressively simpler forms.
    key = target.lower()
    hits = basename_idx.get(key, [])
    if hits:
        return hits[0]
    # Try the basename (last path component) both with and without ext.
    base = key.rsplit("/", 1)[-1]
    if base != key:
        hits = basename_idx.get(base, [])
        if hits:
            return hits[0]
        # Strip extension from base.
        if "." in base:
            stem = base.rsplit(".", 1)[0]
            hits = basename_idx.get(stem, [])
            if hits:
                return hits[0]
    # If target was given with an extension, try stripping.
    if "." in key:
        stem = key.rsplit(".", 1)[0]
        hits = basename_idx.get(stem, [])
        if hits:
            return hits[0]
    return None


# ----------------------------------------------------------------------------
# Dependency extraction
# ----------------------------------------------------------------------------

class DepCollector:
    """Accumulate per-folder dependencies at fine + coarse granularity."""

    def __init__(self):
        # fwd[src_folder][tgt_folder][kind] = count
        self.fine_fwd: dict[str, dict[str, Counter]] = defaultdict(
            lambda: defaultdict(Counter))
        self.coarse_fwd: dict[str, dict[str, Counter]] = defaultdict(
            lambda: defaultdict(Counter))
        self.folder_cite_count: Counter = Counter()     # fine folder
        self.folder_cite_count_coarse: Counter = Counter()
        self.all_fine_folders: set[str] = set()
        self.all_coarse_folders: set[str] = set()
        self.orphan_inputs: list[tuple[Path, str]] = []
        self.orphan_refs: list[tuple[Path, str]] = []
        self.orphan_wikilinks: list[tuple[Path, str]] = []
        self.totals = Counter()

    def register_folder(self, fine: str):
        self.all_fine_folders.add(fine)
        self.all_coarse_folders.add(coarse_folder_of(fine))

    def add(self, src_fp: Path, tgt_fp: Path, kind: str, root: Path):
        src_rel = src_fp.relative_to(root)
        tgt_rel = tgt_fp.relative_to(root)
        src_fine = fine_folder_of(src_rel)
        tgt_fine = fine_folder_of(tgt_rel)
        self.register_folder(src_fine)
        self.register_folder(tgt_fine)
        if src_fine == tgt_fine:
            return  # self-folder edges don't count at folder-level
        self.fine_fwd[src_fine][tgt_fine][kind] += 1
        src_c = coarse_folder_of(src_fine)
        tgt_c = coarse_folder_of(tgt_fine)
        if src_c != tgt_c:
            self.coarse_fwd[src_c][tgt_c][kind] += 1
        self.totals[kind] += 1


def process_tex_files(tex_files, label_map, bib_keys, dc: DepCollector,
                      root: Path):
    for fp in tex_files:
        if is_skipped(fp, root):
            continue
        txt = strip_comments(read_text(fp))

        # \input / \include ----------------------------------------------
        for m in list(RE_INPUT.finditer(txt)) + list(RE_INCLUDE.finditer(txt)):
            raw = m.group(1).strip()
            bare = raw[:-4] if raw.endswith(".tex") else raw
            if bare.split("/")[-1] in SKIP_INPUT_TARGETS:
                continue
            resolved = resolve_input_target(raw, fp, root)
            if resolved is None:
                dc.orphan_inputs.append((fp, raw))
                continue
            dc.add(fp, resolved, "input", root)

        # \ref etc. -----------------------------------------------------
        for m in RE_REF.finditer(txt):
            for label in m.group(1).split(","):
                label = label.strip()
                if not label:
                    continue
                tgt = label_map.get(label)
                if tgt is None:
                    dc.orphan_refs.append((fp, label))
                    continue
                dc.add(fp, tgt, "ref", root)

        # \cite --- counted per folder only -----------------------------
        cite_hits = 0
        for m in RE_CITE.finditer(txt):
            for key in m.group(1).split(","):
                key = key.strip()
                if not key:
                    continue
                if bib_keys and key not in bib_keys:
                    # Still count if key is unknown? We count known-in-bib.
                    # But user said: "count how many times bib keys are cited"
                    # -> require membership in bib.
                    continue
                cite_hits += 1
        if cite_hits:
            fine = fine_folder_of(fp.resolve().relative_to(root))
            dc.folder_cite_count[fine] += cite_hits
            dc.folder_cite_count_coarse[coarse_folder_of(fine)] += cite_hits
            dc.totals["cite"] += cite_hits


def process_md_files(md_files, basename_idx, dc: DepCollector, root: Path):
    for fp in md_files:
        if is_skipped(fp, root):
            continue
        # Register even if no deps.
        dc.register_folder(fine_folder_of(fp.resolve().relative_to(root)))
        txt = read_text(fp)
        for m in RE_WIKILINK.finditer(txt):
            target = m.group(1).strip()
            if not target:
                continue
            resolved = resolve_wikilink(target, basename_idx, root)
            if resolved is None:
                dc.orphan_wikilinks.append((fp, target))
                continue
            dc.add(fp, resolved, "wikilink", root)


# ----------------------------------------------------------------------------
# Bibliography key extraction
# ----------------------------------------------------------------------------

RE_BIB_ENTRY = re.compile(r"@\w+\s*\{\s*([^,\s]+)\s*,", re.IGNORECASE)


def load_bib_keys(bib_path: Path) -> set[str]:
    if not bib_path.exists():
        return set()
    txt = read_text(bib_path)
    return {m.group(1).strip() for m in RE_BIB_ENTRY.finditer(txt)}


# ----------------------------------------------------------------------------
# Rendering
# ----------------------------------------------------------------------------

def format_edge_list(edges: dict[str, Counter]) -> list[str]:
    """Render a dict {target_folder: Counter(kind->count)} as bullet lines."""
    if not edges:
        return ["  - --"]
    out: list[str] = []
    for tgt in sorted(edges.keys()):
        parts = []
        c = edges[tgt]
        total = sum(c.values())
        kind_bits = ", ".join(f"{k}:{c[k]}" for k in sorted(c) if c[k])
        out.append(f"  - {tgt} ({total}x | {kind_bits})")
    return out


def split_by_kind(edges: dict[str, Counter]):
    """Return three dicts: inputs, refs, wikilinks, each target->count."""
    by_input: Counter = Counter()
    by_ref:   Counter = Counter()
    by_wiki:  Counter = Counter()
    for tgt, c in edges.items():
        if c["input"]:
            by_input[tgt] += c["input"]
        if c["ref"]:
            by_ref[tgt] += c["ref"]
        if c["wikilink"]:
            by_wiki[tgt] += c["wikilink"]
    return by_input, by_ref, by_wiki


def render_forward(dc: DepCollector, bib_keys: set[str],
                   root: Path) -> str:
    L: list[str] = []
    today = date.today().isoformat()
    L.append("# Dependencies - TGP/TGP_v1 (forward graph)")
    L.append("")
    L.append(f"> Generated: {today} by `tooling/build_deps_graph.py`")
    L.append("> Sources: \\input{} + \\ref{} + \\cite{} + [[wikilink]]")
    L.append("")
    L.append("## Summary")
    L.append("")
    L.append(f"- Folders analyzed (fine): {len(dc.all_fine_folders)}")
    L.append(f"- Folders analyzed (coarse): {len(dc.all_coarse_folders)}")
    L.append(f"- Total dependencies found: "
             f"{dc.totals['input'] + dc.totals['ref'] + dc.totals['wikilink']}")
    L.append(f"  - `\\input`  edges: {dc.totals['input']}")
    L.append(f"  - `\\ref`    edges: {dc.totals['ref']}")
    L.append(f"  - `[[wiki]]` edges: {dc.totals['wikilink']}")
    L.append(f"- `\\cite{{}}` usages counted (bib keys): {dc.totals['cite']}")
    L.append(f"- Bibliography keys in tgp_main.bib: {len(bib_keys)}")
    L.append("")
    L.append(f"- Orphans (unresolved): "
             f"\\input={len(dc.orphan_inputs)}, "
             f"\\ref={len(dc.orphan_refs)}, "
             f"wikilink={len(dc.orphan_wikilinks)}")
    L.append("")
    # Coarse view -------------------------------------------------------
    L.append("## Coarse view (top-level folders)")
    L.append("")
    coarse_sorted = sorted(dc.all_coarse_folders)
    for src in coarse_sorted:
        targets = dc.coarse_fwd.get(src, {})
        # Only UNIQUE top-folders as deps, no counters (per spec).
        uniq = sorted(t for t in targets.keys() if t != src)
        deps_str = ", ".join(uniq) if uniq else "--"
        L.append(f"### {src}")
        L.append(f"- depends on: {deps_str}")
        L.append("")
    # Fine view ---------------------------------------------------------
    L.append("## Fine view (per subfolder)")
    L.append("")
    fine_sorted = sorted(dc.all_fine_folders)
    for src in fine_sorted:
        targets = dc.fine_fwd.get(src, {})
        by_input, by_ref, by_wiki = split_by_kind(targets)
        L.append(f"### {src}")

        def _fmt(counter: Counter) -> list[str]:
            if not counter:
                return ["  - --"]
            out = []
            for tgt in sorted(counter):
                out.append(f"  - {tgt} ({counter[tgt]}x)")
            return out

        L.append("- `\\input`:")
        L.extend(_fmt(by_input))
        L.append("- `\\ref` to other folders:")
        L.extend(_fmt(by_ref))
        L.append("- Wikilinks to:")
        L.extend(_fmt(by_wiki))
        L.append(f"- Cite count: {dc.folder_cite_count.get(src, 0)}")
        L.append("")
    # Orphans -----------------------------------------------------------
    if dc.orphan_inputs or dc.orphan_refs or dc.orphan_wikilinks:
        L.append("## Unresolved references (orphans)")
        L.append("")
        if dc.orphan_inputs:
            L.append(f"### \\input ({len(dc.orphan_inputs)})")
            for src, tgt in dc.orphan_inputs[:50]:
                rel = src.resolve().relative_to(root)
                L.append(f"- `{rel.as_posix()}` -> `{tgt}`")
            if len(dc.orphan_inputs) > 50:
                L.append(f"- ... +{len(dc.orphan_inputs) - 50} more")
            L.append("")
        if dc.orphan_refs:
            L.append(f"### \\ref ({len(dc.orphan_refs)})")
            shown = dc.orphan_refs[:50]
            for src, label in shown:
                rel = src.resolve().relative_to(root)
                L.append(f"- `{rel.as_posix()}` -> `{label}`")
            if len(dc.orphan_refs) > 50:
                L.append(f"- ... +{len(dc.orphan_refs) - 50} more")
            L.append("")
        if dc.orphan_wikilinks:
            L.append(f"### wikilinks ({len(dc.orphan_wikilinks)})")
            shown = dc.orphan_wikilinks[:50]
            for src, target in shown:
                rel = src.resolve().relative_to(root)
                L.append(f"- `{rel.as_posix()}` -> `[[{target}]]`")
            if len(dc.orphan_wikilinks) > 50:
                L.append(f"- ... +{len(dc.orphan_wikilinks) - 50} more")
            L.append("")
    return "\n".join(L) + "\n"


def render_reverse(dc: DepCollector, root: Path) -> str:
    # Build reverse maps.
    fine_rev: dict[str, dict[str, Counter]] = defaultdict(
        lambda: defaultdict(Counter))
    for src, tmap in dc.fine_fwd.items():
        for tgt, c in tmap.items():
            for kind, count in c.items():
                fine_rev[tgt][src][kind] += count

    coarse_rev: dict[str, dict[str, Counter]] = defaultdict(
        lambda: defaultdict(Counter))
    for src, tmap in dc.coarse_fwd.items():
        for tgt, c in tmap.items():
            for kind, count in c.items():
                coarse_rev[tgt][src][kind] += count

    L: list[str] = []
    today = date.today().isoformat()
    L.append("# Reverse Dependencies - TGP/TGP_v1")
    L.append("")
    L.append(f"> Generated: {today} by `tooling/build_deps_graph.py`")
    L.append("> Sources: \\input{} + \\ref{} + \\cite{} + [[wikilink]]")
    L.append("")
    L.append("## Coarse view")
    L.append("")
    for tgt in sorted(dc.all_coarse_folders):
        srcs = coarse_rev.get(tgt, {})
        uniq = sorted(s for s in srcs.keys() if s != tgt)
        dep_str = ", ".join(uniq) if uniq else "--"
        L.append(f"### {tgt}")
        L.append(f"- depended on by: {dep_str}")
        L.append("")
    L.append("## Fine view (per subfolder)")
    L.append("")
    for tgt in sorted(dc.all_fine_folders):
        srcs = fine_rev.get(tgt, {})
        L.append(f"### {tgt}")
        if not srcs:
            L.append("- depended on by: --")
            L.append("")
            continue
        L.append("- depended on by:")
        for src in sorted(srcs.keys()):
            c = srcs[src]
            total = sum(c.values())
            bits = ", ".join(f"{k}:{c[k]}" for k in sorted(c) if c[k])
            L.append(f"  - {src} ({total}x | {bits})")
        L.append("")
    return "\n".join(L) + "\n"


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main() -> int:
    root = find_root()
    print(f"[info] root = {root}", file=sys.stderr)

    bib_keys = load_bib_keys(root / "tgp_main.bib")
    print(f"[info] bib keys: {len(bib_keys)}", file=sys.stderr)

    tex_files, md_files, asset_files = collect_files(root)
    print(f"[info] tex files: {len(tex_files)}, md files: {len(md_files)}, "
          f"asset files indexed for wikilinks: {len(asset_files)}",
          file=sys.stderr)

    label_map = build_label_map(tex_files, root)
    print(f"[info] labels indexed: {len(label_map)}", file=sys.stderr)

    basename_idx = build_basename_index(root, tex_files, md_files, asset_files)

    dc = DepCollector()

    # Register every folder up-front (so empties appear).
    for fp in tex_files + md_files:
        try:
            dc.register_folder(fine_folder_of(fp.resolve().relative_to(root)))
        except ValueError:
            pass

    process_tex_files(tex_files, label_map, bib_keys, dc, root)
    process_md_files(md_files, basename_idx, dc, root)

    fwd_out = render_forward(dc, bib_keys, root)
    rev_out = render_reverse(dc, root)

    (root / "DEPENDENCIES.md").write_text(fwd_out, encoding="utf-8")
    (root / "DEPENDENCIES_REVERSE.md").write_text(rev_out, encoding="utf-8")
    print(f"[info] wrote {root / 'DEPENDENCIES.md'}", file=sys.stderr)
    print(f"[info] wrote {root / 'DEPENDENCIES_REVERSE.md'}", file=sys.stderr)

    # Echo totals for convenience.
    print(
        "[totals] inputs={inputs} refs={refs} wikilinks={wl} cites={cites} "
        "orphans(inputs/refs/wiki)={oi}/{orr}/{ow}".format(
            inputs=dc.totals["input"],
            refs=dc.totals["ref"],
            wl=dc.totals["wikilink"],
            cites=dc.totals["cite"],
            oi=len(dc.orphan_inputs),
            orr=len(dc.orphan_refs),
            ow=len(dc.orphan_wikilinks),
        ),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
