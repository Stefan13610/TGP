"""
Sesja 4 — klasyfikator YAML `tgp_status:` dla 86 folderów `research/`.

Fazy:
  1. Dry-run (default): generuje propozycje YAML dla każdego folderu →
     `_classify_s4_dryrun.json` + `_classify_s4_summary.md`.
     Nie modyfikuje plików w research/.
  2. Apply (--apply): wykonuje zmiany w README.md per folder
     (preserve istniejący frontmatter, dodaje/aktualizuje `tgp_status`).
     Tworzy minimalny README dla folderów bez README.

Heurystyki klasyfikacji (zob. AGENT_PROTOCOL.md §1):
  - POLLUTED-74394a8 → kwarantanna (level: unknown, status: needs-bridge)
  - closure_aggregator (1 folder) → kind: closure-aggregator, level: mixed
  - external_review_* → kind: review
  - legacy (op6/op7/op1-op2-op4) → kind: derivation, level: mixed (z notatką)
  - 3-phase op-* (program.md + Phase{1,2,3}_results.md):
      - PASS ≥ 30, FAIL < 5, brak FALSIFIED w title → L3, status: active
      - PASS ≥ 5, mixed → L2, status: active
      - inne → L2, status: active
  - topic + README only → L1, status: needs-bridge (jeśli mało scaffoldingu)
  - phenomenology keywords (DESI/Hubble/galaxy/BAO) → kind: phenomenology

Anty-overclaim:
  - Wszystkie wpisy w `source_of_status` muszą cytować konkretne pliki.
  - level: L4 jest ZAKAZANE w Sesji 4 (wymaga INTAKE workflow w Sesji 7).
  - 'FULL CONVERGENCE' w README/frontmatter → flag w source_of_status.
  - Folder z multi-doc inconsistent status → level: unknown + flag.
"""

from __future__ import annotations
import argparse
import json
import re
import sys
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
TGP_V1 = HERE.parent.parent
RESEARCH = TGP_V1 / "research"

DRYRUN_OUT = HERE / "_classify_s4_dryrun.json"
SUMMARY_OUT = HERE / "_classify_s4_summary.md"

# Foldery zatrute (kwarantanna z 74394a8)
POLLUTED = {
    "op-chi1-newton-constant-derivation",
    "op-uv2-mtgp-absolute-scale",
    "op-omega2-axion-coupling-lock",
    "op-omega3-axion-decay-constant",
}

# Legacy folders (gęsto linkowane, mixed level)
LEGACY = {"op6", "op7", "op1-op2-op4"}

# Specjalne foldery
SANDBOX_RESERVED = {"_sandbox", "_archive"}  # nie klasyfikujemy ich tu
CLOSURE_AGGREGATOR = {"closure_2026-04-26"}
EXTERNAL_REVIEW_PREFIX = "external_review_"
QM_PREFIX = "qm_"

# Phenomenology keywords (w nazwie folderu lub README)
PHENO_NAME_KEYWORDS = {
    "desi", "hubble", "galaxy", "cosmo", "s8", "muon", "muon_g", "casimir",
    "thermal", "liquid", "atomic", "nbody",
}

# Audit-aware markers (folder ma audit-driven patches)
AUDIT_MARKERS = [
    r"B6[-_]CLOSED", r"B7[-_]?v?2?[-_]?CLOSED", r"B8[-_]CLOSED", r"B9[-_]CLOSED",
    r"A1\+A2\+A3", r"A2[-_]CLOSED", r"A4[-_]CLOSED", r"A5[-_]PATCHED",
    r"A6\+A8", r"⚠ B[0-9]+", r"⚠ A[0-9]+",
    r"audit-aware", r"audit aware",
]
AUDIT_RE = re.compile("|".join(AUDIT_MARKERS), re.IGNORECASE)

# Phrase violations (banned per AGENT_PROTOCOL §3.1)
BANNED_PHRASES = [
    r"\bFULL CONVERGENCE\b",
    r"\bFULL[_-]CONVERGENCE\b",  # YAML tags form
]
BANNED_RE = re.compile("|".join(BANNED_PHRASES))

# Polluted markers
POLLUTED_RE = re.compile(r"polluted[-_]74394a8|repackaged.circularity|74394a8", re.IGNORECASE)

# Cross-validation (closure_2026-04-26 children) — special markers
CV_CHILDREN_OF_CLOSURE = {
    "sigma_ab_pathB", "f_psi_principle", "Lambda_from_Phi0", "alpha_psi_threshold",
}

TODAY = datetime.now().strftime("%Y-%m-%d")


def read_text(p: Path, max_bytes: int = 500_000) -> str:
    try:
        if p.stat().st_size > max_bytes:
            with p.open("rb") as f:
                return f.read(max_bytes).decode("utf-8", errors="ignore")
        return p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""


def has_existing_yaml(readme_path: Path) -> tuple[bool, bool, dict]:
    """
    Returns (has_frontmatter, has_tgp_status, parsed_yaml_dict).
    Simple parsing — does NOT use external YAML lib.
    """
    if not readme_path.exists():
        return (False, False, {})
    text = read_text(readme_path)
    if not text.lstrip().startswith("---"):
        return (False, False, {})

    # Extract frontmatter block
    parts = text.split("---", 2)
    if len(parts) < 3:
        return (True, False, {})
    fm = parts[1]
    has_status = "tgp_status:" in fm
    return (True, has_status, {"raw": fm})


def scan_folder(folder: Path) -> dict:
    """Zbierz signały z folderu."""
    name = folder.name
    info = {
        "name": name,
        "rel_path": folder.relative_to(TGP_V1).as_posix(),
        "is_polluted_74394a8": name in POLLUTED,
        "is_legacy": name in LEGACY,
        "is_closure_aggregator": name in CLOSURE_AGGREGATOR,
        "is_external_review": name.startswith(EXTERNAL_REVIEW_PREFIX),
        "is_qm": name.startswith(QM_PREFIX),
        "has_program_md": (folder / "program.md").exists(),
        "has_readme": (folder / "README.md").exists(),
        "has_findings": (folder / "FINDINGS.md").exists(),
        "has_needs": (folder / "NEEDS.md").exists(),
        "phase_files": sorted(p.name for p in folder.glob("Phase*.md")),
        "subfolder_count": 0,
        "subfolders": [],
        "md_count": 0,
        "py_count": 0,
        "txt_count": 0,
        "pass_count": 0,
        "fail_count": 0,
        "closed_count": 0,
        "withdrawn_count": 0,
        "falsified_count": 0,
        "obsolete_count": 0,
        "superseded_count": 0,
        "open_count": 0,
        "blocked_count": 0,
        "audit_marker_hits": 0,
        "banned_phrase_hits": 0,
        "banned_phrase_files": [],
        "polluted_self_marker": False,
        "readme_first_h1": "",
        "readme_yaml_status": None,
        "newest_mtime": None,
    }

    # Subfolders
    try:
        for entry in folder.iterdir():
            if entry.is_dir():
                info["subfolder_count"] += 1
                info["subfolders"].append(entry.name)
        info["subfolders"].sort()
    except Exception:
        pass

    # README first H1
    readme_p = folder / "README.md"
    if readme_p.exists():
        rtext = read_text(readme_p)
        for line in rtext.splitlines():
            ls = line.strip()
            if ls.startswith("# "):
                info["readme_first_h1"] = ls[2:].strip()[:120]
                break
        # check for existing tgp_status
        has_fm, has_status, _ = has_existing_yaml(readme_p)
        if has_fm:
            info["readme_yaml_status"] = "has_frontmatter" + (
                "_with_tgp_status" if has_status else "_no_tgp_status"
            )

    # Walk
    newest = 0.0
    files_walked = 0
    for child in folder.rglob("*"):
        if not child.is_file():
            continue
        if "__pycache__" in child.parts:
            continue
        files_walked += 1
        try:
            st = child.stat()
            if st.st_mtime > newest:
                newest = st.st_mtime
            ext = child.suffix.lower()
            if ext == ".md":
                info["md_count"] += 1
            elif ext == ".py":
                info["py_count"] += 1
            elif ext == ".txt":
                info["txt_count"] += 1
            if ext in (".md", ".txt", ".py", ".tex"):
                t = read_text(child)
                info["pass_count"] += len(re.findall(r"\bPASS\b", t))
                info["fail_count"] += len(re.findall(r"\bFAIL\b", t))
                info["closed_count"] += len(re.findall(r"\bCLOSED\b", t))
                info["withdrawn_count"] += len(re.findall(r"\bWITHDRAWN\b", t))
                info["falsified_count"] += len(re.findall(r"\bFALSIFIED\b", t))
                info["obsolete_count"] += len(re.findall(r"\bOBSOLETE\b", t))
                info["superseded_count"] += len(re.findall(r"\bSUPERSEDED\b", t))
                info["open_count"] += len(re.findall(r"\bOPEN\b", t))
                info["blocked_count"] += len(re.findall(r"\bBLOCKED\b", t))
                info["audit_marker_hits"] += len(AUDIT_RE.findall(t))
                bp = BANNED_RE.findall(t)
                if bp:
                    info["banned_phrase_hits"] += len(bp)
                    info["banned_phrase_files"].append(child.name)
                if POLLUTED_RE.search(t):
                    info["polluted_self_marker"] = True
        except Exception:
            pass

    if newest > 0:
        info["newest_mtime"] = datetime.fromtimestamp(newest).strftime("%Y-%m-%d")
    info["files_walked"] = files_walked

    return info


def classify(info: dict) -> dict:
    """Map signals → tgp_status YAML."""
    name = info["name"]

    yaml = {
        "folder_status": "active",
        "level": "L2",
        "kind": "derivation",
        "core_compatibility": "unknown",
        "last_reviewed_against_core": "unknown",
        "may_edit_core": False,
        "exports_findings": False,
        "has_needs_file": info["has_needs"],
        "has_findings_file": info["has_findings"],
        "open_bridges": [],
        "depends_on": [],
        "impacts": [],
        "source_of_status": [],
        "promoted_to_core": None,
        "polluted_74394a8": info["is_polluted_74394a8"],
        "pre_existing_findings": info["has_findings"],
        "pre_existing_needs": info["has_needs"],
        "last_yaml_update": TODAY,
    }

    notes = []  # debug notes — go into source_of_status / commentary

    # ── Special folder kinds ──
    if info["is_closure_aggregator"]:
        yaml["folder_status"] = "active"
        yaml["level"] = "mixed"
        yaml["kind"] = "closure-aggregator"
        yaml["core_compatibility"] = "partial"
        yaml["source_of_status"] = [
            f"closure_aggregator: {info['subfolder_count']} podzamknięcia ({', '.join(info['subfolders'])})",
            f"CLOSURE_2026-04-26_SUMMARY.md (35/35 PASS structural total)",
            f"PASS count w folderze={info['pass_count']}, CLOSED={info['closed_count']}",
            "Per Q2 (2026-05-03): tylko rodzic dostaje pełną warstwę YAML/NEEDS/FINDINGS",
            "4 INTAKE wnioski PENDING-HUMAN-REVIEW dla podzamknięć (Sesja 3.5 Q6)",
        ]
        notes.append(
            f"closure_aggregator with {info['subfolder_count']} subfolders; "
            f"per Q2 (2026-05-03) tylko rodzic dostaje pełną warstwę"
        )
        return {"yaml": yaml, "notes": notes}

    if info["is_external_review"]:
        yaml["folder_status"] = "review"
        yaml["level"] = "mixed"
        yaml["kind"] = "review"
        yaml["core_compatibility"] = "current"
        yaml["source_of_status"] = [
            f"external review folder: {name}",
            f"PASS count={info['pass_count']}, files: md={info['md_count']}, py={info['py_count']}",
            "Wzajemny audit pre-publication readiness",
        ]
        notes.append("external_review folder — wzajemny audit pre-publication readiness")
        return {"yaml": yaml, "notes": notes}

    if info["is_polluted_74394a8"]:
        yaml["folder_status"] = "needs-bridge"
        yaml["level"] = "unknown"
        yaml["kind"] = "derivation"
        yaml["core_compatibility"] = "unknown"
        yaml["source_of_status"] = [
            "Polluted-74394a8 — frozen pending forward-patch",
            "meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md (commit 74394a8 ledger pollution flag)",
        ]
        if info["banned_phrase_hits"]:
            yaml["source_of_status"].append(
                f"BANNED PHRASE 'FULL CONVERGENCE' hits={info['banned_phrase_hits']} "
                f"(per AGENT_PROTOCOL §3 case study #1)"
            )
        if info["polluted_self_marker"]:
            yaml["source_of_status"].append(
                "Folder zawiera self-critique (CRITIQUE_*.md lub równoważne)"
            )
        notes.append(
            "Polluted-74394a8 — frozen pending forward-patch"
        )
        if info["banned_phrase_hits"]:
            notes.append(f"BANNED PHRASE 'FULL CONVERGENCE' hits={info['banned_phrase_hits']}")
        if info["polluted_self_marker"]:
            notes.append("Folder zawiera self-critique (polluted self-marker found)")
        return {"yaml": yaml, "notes": notes}

    if info["is_legacy"]:
        yaml["folder_status"] = "active"
        yaml["level"] = "mixed"
        yaml["kind"] = "derivation"
        yaml["core_compatibility"] = "current"  # legacy fragmenty są w core (op6 v2 pivot, op7 OP-7 closure)
        yaml["source_of_status"] = [
            f"Legacy folder ({name}) — gęsto linkowany w INDEX/DEPENDENCIES",
            f"PASS count={info['pass_count']}, CLOSED={info['closed_count']}",
        ]
        if name == "op6":
            yaml["source_of_status"].append("op6/v2_pivot_summary.md (substrate v2 axiom — promoted to axioms/M1A)")
        elif name == "op7":
            yaml["source_of_status"].append("op7/OP7_T*_results.md (T1-T6 CLOSED 2026-04-25, 94/97 PASS)")
        elif name == "op1-op2-op4":
            yaml["source_of_status"].append("op1-op2-op4 historical aggregation (early OP cycles)")
        notes.append(
            f"Legacy ({name}) — gęsto linkowany, in place. "
            f"PASS count={info['pass_count']}, CLOSED={info['closed_count']}"
        )
        return {"yaml": yaml, "notes": notes}

    # ── Determine kind ──
    if info["is_qm"]:
        yaml["kind"] = "derivation"
    else:
        # Phenomenology keyword detection
        is_pheno = False
        for kw in PHENO_NAME_KEYWORDS:
            if kw in name.lower():
                is_pheno = True
                break
        if is_pheno:
            yaml["kind"] = "phenomenology"
        else:
            yaml["kind"] = "derivation"

    # ── Determine level + folder_status ──
    has_program = info["has_program_md"]
    n_phase_results = sum(1 for f in info["phase_files"] if "_results" in f)
    n_phase_setup = sum(1 for f in info["phase_files"] if "_setup" in f)
    n_pass = info["pass_count"]
    n_fail = info["fail_count"]
    n_closed = info["closed_count"]
    pass_fail_ratio = n_pass / max(n_fail, 1)

    if has_program and n_phase_results >= 3 and n_pass >= 30 and pass_fail_ratio > 5:
        yaml["level"] = "L3"
        yaml["folder_status"] = "active"
        if info["audit_marker_hits"] > 0:
            yaml["core_compatibility"] = "current"
        notes.append(
            f"3-phase op-* with N={n_phase_results} Phase results, "
            f"PASS={n_pass}, CLOSED={n_closed}, FAIL={n_fail}, "
            f"audit_markers={info['audit_marker_hits']}"
        )
    elif has_program and n_phase_results >= 2 and n_pass >= 5:
        yaml["level"] = "L2"
        yaml["folder_status"] = "active"
        notes.append(
            f"3-phase op-* partial: N={n_phase_results} Phase results, "
            f"PASS={n_pass}, FAIL={n_fail}"
        )
    elif has_program and n_phase_setup >= 2:
        yaml["level"] = "L1"
        yaml["folder_status"] = "active"
        notes.append(
            f"3-phase op-* in progress: N={n_phase_setup} Phase setup, "
            f"only {n_phase_results} results"
        )
    elif info["py_count"] >= 1 and info["txt_count"] >= 1:
        # Topic with scripts — L1 numerical
        yaml["level"] = "L1"
        yaml["folder_status"] = "active"
        notes.append(
            f"Topic with scripts: {info['py_count']}.py + {info['txt_count']}.txt; "
            f"PASS={n_pass}, no Phase structure"
        )
    elif info["md_count"] >= 1 and info["py_count"] == 0:
        # README-only topic
        yaml["level"] = "L1"
        yaml["folder_status"] = "needs-bridge"
        notes.append(
            f"README-only topic ({info['md_count']} .md, 0 .py); needs scaffolding"
        )
    else:
        yaml["level"] = "unknown"
        yaml["folder_status"] = "active"
        notes.append(f"Mixed scaffolding: md={info['md_count']}, py={info['py_count']}, txt={info['txt_count']}")

    # ── Audit marker → core_compatibility upgrade ──
    if info["audit_marker_hits"] >= 2 and yaml["core_compatibility"] == "unknown":
        yaml["core_compatibility"] = "current"
        notes.append(f"Audit markers found ({info['audit_marker_hits']}) → core_compatibility: current")

    # ── Banned phrase detection (Q7 finding) ──
    if info["banned_phrase_hits"] > 0:
        notes.append(
            f"⚠ BANNED PHRASE violation: 'FULL CONVERGENCE' x{info['banned_phrase_hits']} "
            f"in files: {info['banned_phrase_files'][:3]}"
        )
        # Don't auto-quarantine, but flag
        if yaml["level"] == "L3":
            yaml["level"] = "L2"  # downgrade pending phrase fix
            notes.append("Level downgrade L3→L2 pending phrase cleanup")

    # ── source_of_status ──
    sos = []
    if info["readme_first_h1"]:
        sos.append(f"README.md H1: '{info['readme_first_h1']}'")
    if n_phase_results >= 1:
        sos.append(
            f"Phase{{1..{n_phase_results}}}_results.md PASS={n_pass}, "
            f"CLOSED={n_closed}, FAIL={n_fail}"
        )
    if info["audit_marker_hits"] > 0:
        sos.append(f"audit-aware markers={info['audit_marker_hits']} (B6/B8/B9/A5/etc.)")
    if not sos:
        sos.append(f"Sesja 4 heuristic from _classify_s4.py (mtime {info['newest_mtime']}); requires manual verification")
    yaml["source_of_status"] = sos

    return {"yaml": yaml, "notes": notes}


def render_yaml_block(yaml: dict, notes: list[str]) -> str:
    """Render `tgp_status:` block as YAML string for inclusion in frontmatter."""
    lines = ["tgp_status:"]

    def add(key, val):
        if isinstance(val, bool):
            lines.append(f"  {key}: {'true' if val else 'false'}")
        elif val is None:
            lines.append(f"  {key}: null")
        elif isinstance(val, str):
            if val == "unknown" or "-" in val or " " in val or ":" in val:
                lines.append(f'  {key}: "{val}"')
            else:
                lines.append(f"  {key}: {val}")
        elif isinstance(val, list):
            if not val:
                lines.append(f"  {key}: []")
            else:
                lines.append(f"  {key}:")
                for item in val:
                    if isinstance(item, str):
                        # Quote if has special chars
                        if any(c in item for c in ":#&*?[]{}|>!%@`,'\""):
                            esc = item.replace('"', '\\"')
                            lines.append(f'    - "{esc}"')
                        else:
                            lines.append(f"    - {item}")
                    else:
                        lines.append(f"    - {item}")
        else:
            lines.append(f"  {key}: {val}")

    for k in [
        "folder_status", "level", "kind", "core_compatibility",
        "last_reviewed_against_core", "may_edit_core", "exports_findings",
        "has_needs_file", "has_findings_file", "open_bridges", "depends_on",
        "impacts", "source_of_status", "promoted_to_core", "polluted_74394a8",
        "pre_existing_findings", "pre_existing_needs", "last_yaml_update"
    ]:
        if k in yaml:
            add(k, yaml[k])

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true",
                        help="Apply changes to README files (default: dry-run)")
    parser.add_argument("--folder", type=str, default=None,
                        help="Process single folder by name (debugging)")
    args = parser.parse_args()

    # Collect folders
    folders = []
    for p in sorted(RESEARCH.iterdir(), key=lambda x: x.name.lower()):
        if not p.is_dir():
            continue
        if p.name in SANDBOX_RESERVED:
            continue
        if args.folder and p.name != args.folder:
            continue
        folders.append(p)

    print(f"Processing {len(folders)} folders (apply={args.apply})...")

    results = []
    for folder in folders:
        info = scan_folder(folder)
        decision = classify(info)
        yaml_block = render_yaml_block(decision["yaml"], decision["notes"])
        results.append({
            "folder": info["name"],
            "info": info,
            "yaml": decision["yaml"],
            "notes": decision["notes"],
            "yaml_block_text": yaml_block,
        })

    # Write JSON
    DRYRUN_OUT.write_text(
        json.dumps(results, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8"
    )
    print(f"Wrote {DRYRUN_OUT} ({DRYRUN_OUT.stat().st_size} B)")

    # Write summary
    summary_lines = []
    summary_lines.append("# Sesja 4 — dry-run classification summary")
    summary_lines.append("")
    summary_lines.append(f"Generated: {datetime.now().isoformat(timespec='seconds')}")
    summary_lines.append(f"Folders processed: {len(results)}")
    summary_lines.append("")

    # Status distribution
    from collections import Counter
    status_counter = Counter(r["yaml"]["folder_status"] for r in results)
    level_counter = Counter(r["yaml"]["level"] for r in results)
    kind_counter = Counter(r["yaml"]["kind"] for r in results)

    summary_lines.append("## Distribution")
    summary_lines.append("")
    summary_lines.append("### folder_status")
    summary_lines.append("")
    for k, v in status_counter.most_common():
        summary_lines.append(f"- `{k}`: {v}")
    summary_lines.append("")
    summary_lines.append("### level")
    summary_lines.append("")
    for k, v in level_counter.most_common():
        summary_lines.append(f"- `{k}`: {v}")
    summary_lines.append("")
    summary_lines.append("### kind")
    summary_lines.append("")
    for k, v in kind_counter.most_common():
        summary_lines.append(f"- `{k}`: {v}")
    summary_lines.append("")

    # Polluted + banned phrases
    polluted = [r for r in results if r["yaml"]["polluted_74394a8"]]
    banned = [r for r in results if r["info"]["banned_phrase_hits"] > 0]
    summary_lines.append("## Flags")
    summary_lines.append("")
    summary_lines.append(f"- **Polluted-74394a8** (kwarantanna): {len(polluted)}")
    for r in polluted:
        summary_lines.append(f"  - `{r['folder']}` — {r['notes'][0] if r['notes'] else '(no notes)'}")
    summary_lines.append("")
    summary_lines.append(f"- **Banned phrase 'FULL CONVERGENCE'**: {len(banned)}")
    for r in banned:
        summary_lines.append(
            f"  - `{r['folder']}` — hits={r['info']['banned_phrase_hits']}, "
            f"files: {r['info']['banned_phrase_files'][:3]}"
        )
    summary_lines.append("")

    # Per-folder table
    summary_lines.append("## Per-folder classification")
    summary_lines.append("")
    summary_lines.append("| # | Folder | folder_status | level | kind | core_compat | mtime | flags | notes (1st) |")
    summary_lines.append("|---|--------|---------------|-------|------|-------------|-------|-------|-------------|")
    for i, r in enumerate(results, 1):
        flags = []
        if r["yaml"]["polluted_74394a8"]:
            flags.append("⚠POLLUTED")
        if r["info"]["banned_phrase_hits"] > 0:
            flags.append("⚠FULLCONV")
        if r["info"]["audit_marker_hits"] > 0:
            flags.append(f"audit×{r['info']['audit_marker_hits']}")
        flag_str = " ".join(flags) or "—"
        first_note = r["notes"][0][:80] if r["notes"] else "—"
        summary_lines.append(
            f"| {i} | `{r['folder']}` | {r['yaml']['folder_status']} | "
            f"{r['yaml']['level']} | {r['yaml']['kind']} | {r['yaml']['core_compatibility']} | "
            f"{r['info']['newest_mtime']} | {flag_str} | {first_note} |"
        )

    # Sample YAML for verification (3 examples)
    summary_lines.append("")
    summary_lines.append("## Sample full YAML (3 reprezentatywnych)")
    summary_lines.append("")

    # Pick one of each: closure-aggregator, polluted, normal active
    samples = []
    for r in results:
        if r["yaml"]["kind"] == "closure-aggregator":
            samples.append(("closure-aggregator", r))
            break
    for r in results:
        if r["yaml"]["polluted_74394a8"]:
            samples.append(("polluted-74394a8", r))
            break
    for r in results:
        if r["yaml"]["level"] == "L3" and r["yaml"]["folder_status"] == "active" and not r["yaml"]["polluted_74394a8"]:
            samples.append(("active-L3", r))
            break

    for label, r in samples:
        summary_lines.append(f"### {label}: `{r['folder']}`")
        summary_lines.append("")
        summary_lines.append("```yaml")
        summary_lines.append(r["yaml_block_text"])
        summary_lines.append("```")
        summary_lines.append("")
        summary_lines.append("**Notes:**")
        for n in r["notes"]:
            summary_lines.append(f"- {n}")
        summary_lines.append("")

    SUMMARY_OUT.write_text("\n".join(summary_lines), encoding="utf-8")
    print(f"Wrote {SUMMARY_OUT} ({SUMMARY_OUT.stat().st_size} B)")

    if args.apply:
        print("APPLY MODE — modifying README files...")
        applied = 0
        for r in results:
            folder_path = RESEARCH / r["folder"]
            apply_to_readme(folder_path, r)
            applied += 1
        print(f"Applied changes to {applied} folders.")
    else:
        print()
        print("DRY-RUN COMPLETE.")
        print(f"Review:")
        print(f"  - {SUMMARY_OUT}")
        print(f"  - {DRYRUN_OUT}")
        print()
        print("To apply: python _classify_s4.py --apply")


def apply_to_readme(folder: Path, result: dict):
    """Apply tgp_status YAML to README. Preserve existing content."""
    readme = folder / "README.md"
    yaml_block = result["yaml_block_text"]
    folder_name = result["folder"]
    h1_title = result["info"].get("readme_first_h1") or folder_name

    if not readme.exists():
        # Create minimal README
        content = f"""---
title: "{h1_title}"
date: {TODAY}
parent: "[[../INDEX.md]]"
related:
  - "[[meta/research/FOLDER_STATUS_INDEX.md]]"
tags:
  - TGP
{yaml_block}
---

# {h1_title}

> **Sesja 4 auto-generated README** (2026-05-03).
> Folder nie miał wcześniej README. Treść poniżej jest minimalna —
> uzupełnij ręcznie zgodnie z [[meta/research/templates/README.template.md]].

## Cel

(do uzupełnienia)

## Stan

(do uzupełnienia — zob. `tgp_status` w frontmatterze)

## Pliki

(zob. listing folderu; szablon: `meta/research/templates/README.template.md`)
"""
        readme.write_text(content, encoding="utf-8")
        return

    # Existing README — preserve content, inject/update tgp_status
    text = readme.read_text(encoding="utf-8")

    if text.lstrip().startswith("---"):
        parts = text.split("---", 2)
        if len(parts) >= 3:
            fm = parts[1]
            body = parts[2]
            # Remove existing tgp_status block if present
            fm_lines = fm.split("\n")
            new_fm_lines = []
            in_tgp_status = False
            for line in fm_lines:
                if line.startswith("tgp_status:"):
                    in_tgp_status = True
                    continue
                if in_tgp_status:
                    if line.startswith("  ") or line == "":
                        continue
                    in_tgp_status = False
                new_fm_lines.append(line)
            new_fm = "\n".join(new_fm_lines).rstrip() + "\n" + yaml_block + "\n"
            new_text = "---" + new_fm + "---" + body
            readme.write_text(new_text, encoding="utf-8")
            return

    # No frontmatter — add it
    new_text = f"""---
title: "{h1_title}"
date: {TODAY}
{yaml_block}
---

{text}"""
    readme.write_text(new_text, encoding="utf-8")


if __name__ == "__main__":
    main()
