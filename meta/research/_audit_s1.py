"""
Sesja 1 — audyt struktury research/.

Read-only. Generuje JSON z metrykami per folder + link-count w INDEX/DEPENDENCIES.
Skrypt jest narzędziem audytowym, nie częścią teorii.
Output: stdout (JSON), zapisywany przez wywołującego do pliku.
"""

from __future__ import annotations
import json
import os
import re
import sys
from datetime import datetime
from pathlib import Path

# --- pozycjonowanie ---
HERE = Path(__file__).resolve().parent  # .../TGP/TGP_v1/meta/research
TGP_V1 = HERE.parent.parent              # .../TGP/TGP_v1
VAULT_ROOT = TGP_V1.parent.parent        # .../ObsydnianMain
CANONICAL = TGP_V1 / "research"
VAULT_RESEARCH = VAULT_ROOT / "research"

INDEX_MD = TGP_V1 / "INDEX.md"
DEPENDENCIES_MD = TGP_V1 / "DEPENDENCIES.md"
DEPENDENCIES_REVERSE_MD = TGP_V1 / "DEPENDENCIES_REVERSE.md"

# Foldery zatrute incydentem 74394a8 (z meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md)
POLLUTED_74394a8 = {
    "op-chi1-newton-constant-derivation",
    "op-uv2-mtgp-absolute-scale",
    "op-omega2-axion-coupling-lock",
    "op-omega3-axion-decay-constant",
}

# Specjalne foldery (po nazwie)
LEGACY_STAY = {"op6", "op7", "op1-op2-op4"}

KEYWORDS = {
    "PASS": r"\bPASS\b",
    "FAIL": r"\bFAIL\b",
    "CLOSED": r"\bCLOSED\b",
    "DERIVED": r"\bDERIVED\b",
    "LOCKED": r"\bLOCKED\b",
    "FALSIFIED": r"\bFALSIFIED\b",
    "WITHDRAWN": r"\bWITHDRAWN\b",
    "OBSOLETE": r"\bOBSOLETE\b",
    "SUPERSEDED": r"\bSUPERSEDED\b",
    "OPEN": r"\bOPEN\b",
    "BLOCKED": r"\bBLOCKED\b",
    "DEFERRED": r"\bDEFERRED\b",
    "STRUCTURAL": r"\bSTRUCTURAL\b",
    "PARTIALLY": r"\bPARTIALLY\b",
    "FULL_CONVERGENCE": r"\bFULL CONVERGENCE\b",
}


def scan_folder(folder: Path) -> dict:
    m = {
        "name": folder.name,
        "rel_path": folder.relative_to(TGP_V1).as_posix(),
        "is_polluted_74394a8": folder.name in POLLUTED_74394a8,
        "is_legacy_stay": folder.name in LEGACY_STAY,
        "has_program_md": (folder / "program.md").exists(),
        "has_readme": (folder / "README.md").exists(),
        "has_findings": (folder / "FINDINGS.md").exists(),
        "has_needs": (folder / "NEEDS.md").exists(),
        "phase_files": [],
        "readme_has_frontmatter": False,
        "readme_has_tgp_status": False,
        "readme_first_h1": "",
        "kw_counts": {k: 0 for k in KEYWORDS},
        "mtime_newest": None,
        "mtime_newest_iso": "",
        "file_count_total": 0,
        "md_count": 0,
        "py_count": 0,
        "txt_count": 0,
        "tex_count": 0,
        "subfolder_count": 0,
        "subfolder_names": [],
        "top_level_files": [],
    }

    # Phase files
    for p in folder.glob("Phase*.md"):
        m["phase_files"].append(p.name)
    m["phase_files"].sort()

    # README inspection
    readme = folder / "README.md"
    if readme.exists():
        try:
            content = readme.read_text(encoding="utf-8", errors="ignore")
            if content.lstrip().startswith("---"):
                m["readme_has_frontmatter"] = True
                if "tgp_status:" in content:
                    m["readme_has_tgp_status"] = True
            for line in content.splitlines():
                line = line.strip()
                if line.startswith("# "):
                    m["readme_first_h1"] = line[2:].strip()[:120]
                    break
        except Exception:
            pass

    # Top-level + subfolders
    try:
        for entry in folder.iterdir():
            if entry.is_dir():
                m["subfolder_count"] += 1
                m["subfolder_names"].append(entry.name)
            else:
                m["top_level_files"].append(entry.name)
        m["subfolder_names"].sort()
        m["top_level_files"].sort()
    except Exception:
        pass

    # Walk for counters + mtime + keywords
    newest = 0.0
    text_blob_size = 0
    MAX_TEXT_PER_FILE = 500_000  # bezpiecznik na ogromne .txt
    for root, dirs, files in os.walk(folder):
        for f in files:
            full = Path(root) / f
            m["file_count_total"] += 1
            ext = f.lower().rsplit(".", 1)[-1] if "." in f else ""
            if ext == "md":
                m["md_count"] += 1
            elif ext == "py":
                m["py_count"] += 1
            elif ext == "txt":
                m["txt_count"] += 1
            elif ext == "tex":
                m["tex_count"] += 1
            try:
                st = full.stat()
                if st.st_mtime > newest:
                    newest = st.st_mtime
                if ext in ("md", "txt", "py", "tex"):
                    if st.st_size > MAX_TEXT_PER_FILE:
                        # czytaj tylko 500KB head
                        with full.open("rb") as fh:
                            buf = fh.read(MAX_TEXT_PER_FILE)
                        text = buf.decode("utf-8", errors="ignore")
                    else:
                        text = full.read_text(encoding="utf-8", errors="ignore")
                    text_blob_size += len(text)
                    for k, pat in KEYWORDS.items():
                        m["kw_counts"][k] += len(re.findall(pat, text))
            except Exception:
                pass

    if newest > 0:
        m["mtime_newest"] = newest
        m["mtime_newest_iso"] = datetime.fromtimestamp(newest).strftime("%Y-%m-%d")
    m["text_blob_size_bytes"] = text_blob_size

    return m


def classify_heuristic(m: dict) -> str:
    """Heurystyka WSTĘPNA — nie nadaje statusu, tylko sugeruje."""
    name = m["name"]
    if name.startswith("closure_"):
        return "closure-aggregator"
    if name.startswith("external_review"):
        return "audit/review"
    if name in LEGACY_STAY:
        return "legacy-stay-in-place"
    if m["is_polluted_74394a8"]:
        return "needs-bridge (polluted-74394a8)"

    kw = m["kw_counts"]
    n_pass = kw["PASS"]
    n_fail = kw["FAIL"]
    n_falsified = kw["FALSIFIED"]
    n_withdrawn = kw["WITHDRAWN"]
    n_obsolete = kw["OBSOLETE"]
    n_superseded = kw["SUPERSEDED"]
    n_open = kw["OPEN"]
    n_phase = len(m["phase_files"])
    has_program = m["has_program_md"]

    # Archive heuristic — silne sygnały (WITHDRAWN/OBSOLETE/SUPERSEDED rzadko
    # występują w aktywnym tekście; FALSIFIED — często, jako reguła falsyfikacji
    # predykcji, NIE status folderu).
    strong_archive = n_withdrawn + n_obsolete + n_superseded
    if strong_archive >= 3:
        return "candidate-archive (heuristic, verify)"
    # FALSIFIED dominujące + niska aktywność:
    if n_falsified >= 10 and n_pass < 3:
        return "candidate-archive (heuristic, verify)"

    if name == "uv_completion":
        return "candidate-promoted-or-active (legacy uv)"

    if has_program and n_phase >= 4:
        if n_pass >= 20 and n_fail < 5:
            return "candidate-core-ready (heuristic, verify)"
        if n_pass >= 5:
            return "candidate-active (3-phase op-*)"
        return "candidate-active (3-phase, partial)"

    # Topic + minimal scaffolding
    if not has_program and n_phase == 0 and m["md_count"] <= 2 and m["py_count"] <= 3:
        return "candidate-needs-bridge-or-sandbox (heuristic)"

    if name.startswith("qm_"):
        return "candidate-needs-bridge (qm cluster)"

    return "candidate-active (heuristic)"


def compare_vault_vs_canonical():
    """Diff folderów wspólnych dla ./research vs TGP/TGP_v1/research."""
    if not VAULT_RESEARCH.exists():
        return {"vault_research_exists": False}
    canonical_set = {p.name for p in CANONICAL.iterdir() if p.is_dir()}
    vault_set = {p.name for p in VAULT_RESEARCH.iterdir() if p.is_dir()}
    only_in_vault = sorted(vault_set - canonical_set)
    only_in_canonical = sorted(canonical_set - vault_set)
    common = sorted(vault_set & canonical_set)

    file_level = []
    for name in common:
        v_dir = VAULT_RESEARCH / name
        c_dir = CANONICAL / name
        v_files = {f.name: f for f in v_dir.iterdir() if f.is_file()}
        c_files = {f.name: f for f in c_dir.iterdir() if f.is_file()}
        only_in_v = sorted(set(v_files) - set(c_files))
        only_in_c = sorted(set(c_files) - set(v_files))
        differ = []
        for fn in sorted(set(v_files) & set(c_files)):
            try:
                a = v_files[fn].read_bytes()
                b = c_files[fn].read_bytes()
                if a != b:
                    differ.append({
                        "file": fn,
                        "size_vault": len(a),
                        "size_canonical": len(b),
                    })
            except Exception:
                pass
        file_level.append({
            "folder": name,
            "only_in_vault": only_in_v,
            "only_in_canonical": only_in_c,
            "content_differs": differ,
            "vault_file_count": len(v_files),
            "canonical_file_count": len(c_files),
        })
    return {
        "vault_research_exists": True,
        "only_in_vault_root": only_in_vault,
        "only_in_canonical": only_in_canonical,
        "common_count": len(common),
        "per_folder_diff": file_level,
    }


def link_audit(folder_names: list[str]) -> dict:
    """
    Liczy, ile referencji do każdego folderu istnieje w INDEX.md / DEPENDENCIES*.md.
    To jest sanity-check dla wariantu minimalnego (oczekiwane: nie ruszamy
    fizycznie, więc 0 broken, ale chcemy znać GĘSTOŚĆ linkowania, żeby
    wiedzieć, które foldery są kruche).
    """
    results = {}
    sources = {
        "INDEX.md": INDEX_MD,
        "DEPENDENCIES.md": DEPENDENCIES_MD,
        "DEPENDENCIES_REVERSE.md": DEPENDENCIES_REVERSE_MD,
    }
    blobs = {}
    for src_name, src_path in sources.items():
        try:
            blobs[src_name] = src_path.read_text(encoding="utf-8", errors="ignore") if src_path.exists() else ""
        except Exception:
            blobs[src_name] = ""

    for name in folder_names:
        per_src = {}
        for src_name, blob in blobs.items():
            # Liczymy: literal "research/<name>" + literal "research\\<name>"
            count = blob.count(f"research/{name}") + blob.count(f"research\\{name}")
            per_src[src_name] = count
        results[name] = per_src
    return results


def main():
    folders = sorted(
        [p for p in CANONICAL.iterdir() if p.is_dir()],
        key=lambda p: p.name.lower()
    )
    metrics = []
    for f in folders:
        m = scan_folder(f)
        m["heuristic_class"] = classify_heuristic(m)
        metrics.append(m)

    folder_names = [m["name"] for m in metrics]
    link_counts = link_audit(folder_names)
    for m in metrics:
        m["link_count"] = link_counts.get(m["name"], {})
        m["link_total"] = sum(m["link_count"].values())

    vault_diff = compare_vault_vs_canonical()

    # Top-level files w research/ (programowe dokumenty, nie foldery)
    top_level_files = sorted([f.name for f in CANONICAL.iterdir() if f.is_file()])

    out = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "canonical_path": CANONICAL.relative_to(TGP_V1).as_posix(),
        "folder_count": len(metrics),
        "top_level_files_in_research": top_level_files,
        "polluted_74394a8": sorted(POLLUTED_74394a8),
        "vault_diff": vault_diff,
        "folders": metrics,
    }

    out_path = HERE / "_audit_s1_raw.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"OK wrote {out_path} ({out_path.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
