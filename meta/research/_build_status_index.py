"""
Build FOLDER_STATUS_INDEX.md from applied YAML statuses (Sesja 4).

Uses _classify_s4_dryrun.json as source of truth (matches what was applied).
Output: FOLDER_STATUS_INDEX.md (overwrites placeholder from Sesja 3).
"""

from __future__ import annotations
import json
from collections import Counter
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
DRYRUN = HERE / "_classify_s4_dryrun.json"
OUT = HERE / "FOLDER_STATUS_INDEX.md"


def main():
    if not DRYRUN.exists():
        print(f"ERROR: {DRYRUN} not found. Run _classify_s4.py first.")
        return

    data = json.loads(DRYRUN.read_text(encoding="utf-8"))
    n = len(data)
    today = datetime.now().strftime("%Y-%m-%d")

    # Counters
    folder_status = Counter(r["yaml"]["folder_status"] for r in data)
    level = Counter(r["yaml"]["level"] for r in data)
    kind = Counter(r["yaml"]["kind"] for r in data)
    core_compat = Counter(r["yaml"]["core_compatibility"] for r in data)
    polluted = sum(1 for r in data if r["yaml"]["polluted_74394a8"])
    banned = sum(1 for r in data if r["info"]["banned_phrase_hits"] > 0)

    lines = []
    lines.append("---")
    lines.append('title: "FOLDER_STATUS_INDEX — globalna mapa wszystkich folderów research/"')
    lines.append(f"date: {today}")
    lines.append("type: index")
    lines.append("status: GENERATED (Sesja 4) — auto-build z YAML w research/<X>/README.md")
    lines.append('parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"')
    lines.append("related:")
    lines.append('  - "[[meta/research/AGENT_PROTOCOL.md]]"')
    lines.append('  - "[[meta/research/AUDIT_RESEARCH_S1.md]]"')
    lines.append('  - "[[meta/research/templates/STATUS_BLOCK.yaml]]"')
    lines.append('  - "[[meta/research/CORE_CANDIDATES.md]]"')
    lines.append('  - "[[meta/research/_classify_s4.py]]"')
    lines.append("tags:")
    lines.append("  - index")
    lines.append("  - folder-status")
    lines.append("  - generated")
    lines.append("---")
    lines.append("")
    lines.append("# FOLDER_STATUS_INDEX — globalna mapa wszystkich folderów `research/`")
    lines.append("")
    lines.append("> **Auto-generated z `_classify_s4.py` + `_build_status_index.py`** —")
    lines.append("> Sesja 4 (2026-05-03). Source of truth: YAML `tgp_status:` w")
    lines.append("> `research/<X>/README.md` per folder (po `--apply`).")
    lines.append(">")
    lines.append("> **Regenerate:** `python meta/research/_build_status_index.py`")
    lines.append(">")
    lines.append(f"> **Liczba folderów:** {n} (po wyłączeniu `_sandbox/` i `_archive/`).")
    lines.append("")

    # ── Distribution ──
    lines.append("## 1. Rozkład statusów")
    lines.append("")
    lines.append("### 1.1 `folder_status`")
    lines.append("")
    lines.append("| Status | Liczba | % |")
    lines.append("|--------|-------:|---:|")
    for k, v in folder_status.most_common():
        lines.append(f"| `{k}` | {v} | {v*100//n}% |")
    lines.append("")
    lines.append("### 1.2 `level`")
    lines.append("")
    lines.append("| Level | Liczba | % |")
    lines.append("|-------|-------:|---:|")
    for k, v in sorted(level.items(), key=lambda x: x[0]):
        lines.append(f"| `{k}` | {v} | {v*100//n}% |")
    lines.append("")
    lines.append("### 1.3 `kind`")
    lines.append("")
    lines.append("| Kind | Liczba | % |")
    lines.append("|------|-------:|---:|")
    for k, v in kind.most_common():
        lines.append(f"| `{k}` | {v} | {v*100//n}% |")
    lines.append("")
    lines.append("### 1.4 `core_compatibility`")
    lines.append("")
    lines.append("| Compatibility | Liczba | % |")
    lines.append("|---------------|-------:|---:|")
    for k, v in core_compat.most_common():
        lines.append(f"| `{k}` | {v} | {v*100//n}% |")
    lines.append("")

    # ── Flags ──
    lines.append("## 2. Flagi")
    lines.append("")
    lines.append(f"- **Polluted-74394a8** (kwarantanna): **{polluted}** folderów")
    lines.append(f"- **Banned phrase 'FULL CONVERGENCE'**: **{banned}** folderów (zob. § 5)")
    lines.append("")

    # ── Per-folder full table ──
    lines.append("## 3. Pełna mapa (sortowana alfabetycznie)")
    lines.append("")
    lines.append("| # | Folder | folder_status | level | kind | core_compat | mtime | flags |")
    lines.append("|---|--------|---------------|-------|------|-------------|-------|-------|")
    for i, r in enumerate(data, 1):
        flags = []
        if r["yaml"]["polluted_74394a8"]:
            flags.append("⚠POLLUTED")
        if r["info"]["banned_phrase_hits"] > 0:
            flags.append(f"⚠FULLCONV×{r['info']['banned_phrase_hits']}")
        if r["info"]["audit_marker_hits"] > 0:
            flags.append(f"audit×{r['info']['audit_marker_hits']}")
        flag_str = " ".join(flags) or "—"
        mtime = r["info"]["newest_mtime"] or "?"
        lines.append(
            f"| {i} | `{r['folder']}` | {r['yaml']['folder_status']} | "
            f"{r['yaml']['level']} | {r['yaml']['kind']} | "
            f"{r['yaml']['core_compatibility']} | {mtime} | {flag_str} |"
        )
    lines.append("")

    # ── Per folder_status grouping ──
    lines.append("## 4. Per `folder_status`")
    lines.append("")
    for status_key in sorted(folder_status.keys()):
        members = sorted([r["folder"] for r in data if r["yaml"]["folder_status"] == status_key])
        lines.append(f"### 4.x. `{status_key}` ({len(members)})")
        lines.append("")
        for f in members:
            lines.append(f"- `research/{f}`")
        lines.append("")

    # ── Per level grouping (only L0-L4 + mixed/unknown) ──
    lines.append("## 5. Per `level`")
    lines.append("")
    for lvl_key in sorted(level.keys()):
        members = sorted([r["folder"] for r in data if r["yaml"]["level"] == lvl_key])
        lines.append(f"### 5.x. `{lvl_key}` ({len(members)})")
        lines.append("")
        for f in members:
            lines.append(f"- `research/{f}`")
        lines.append("")

    # ── Banned phrase violators ──
    lines.append("## 6. Banned phrase 'FULL CONVERGENCE' violators (Q7 finding)")
    lines.append("")
    lines.append("> Per [[meta/research/AGENT_PROTOCOL.md]] §3.0/§3.1: 'FULL CONVERGENCE' jest")
    lines.append("> zarezerwowane. Ta lista wymaga **phrase cleanup w future maintenance pass**")
    lines.append("> (zalecenie z [[meta/research/SANITY_PASS_S3_5.md]] §4.1).")
    lines.append("")
    lines.append("| Folder | Hits | Files (top 3) | Status |")
    lines.append("|--------|-----:|---------------|--------|")
    for r in data:
        if r["info"]["banned_phrase_hits"] == 0:
            continue
        files_short = ", ".join(f"`{f}`" for f in r["info"]["banned_phrase_files"][:3])
        polluted_str = "POLLUTED" if r["yaml"]["polluted_74394a8"] else r["yaml"]["folder_status"]
        lines.append(
            f"| `{r['folder']}` | {r['info']['banned_phrase_hits']} | {files_short} | {polluted_str} |"
        )
    lines.append("")

    # ── Audit-aware folders ──
    lines.append("## 7. Foldery z audit-aware markers (B6/B7/B8/B9/A5 patches)")
    lines.append("")
    lines.append("> Foldery które zawierają patche z audytu 2026-05-01 (B6-CLOSED, B8-CLOSED,")
    lines.append("> B9-CLOSED, A5-PATCHED, etc.). To są foldery o **wzmocnionym** zaufaniu —")
    lines.append("> ich `core_compatibility: current` jest udokumentowany.")
    lines.append("")
    lines.append("| Folder | Audit hits | folder_status | level |")
    lines.append("|--------|-----------:|---------------|-------|")
    audit_folders = [r for r in data if r["info"]["audit_marker_hits"] > 0]
    audit_folders.sort(key=lambda r: -r["info"]["audit_marker_hits"])
    for r in audit_folders:
        lines.append(
            f"| `{r['folder']}` | {r['info']['audit_marker_hits']} | "
            f"{r['yaml']['folder_status']} | {r['yaml']['level']} |"
        )
    lines.append("")

    # ── Anty-overclaim audit ──
    lines.append("## 8. Anti-overclaim audit (CRITICAL flagi)")
    lines.append("")
    issues = []
    for r in data:
        # No L4 in Sesja 4 (correct)
        if r["yaml"]["level"] == "L4":
            issues.append(f"- ❌ `{r['folder']}` ma `level: L4` w Sesji 4 (zakazane bez INTAKE)")
        # Polluted z 'FULL CONVERGENCE' jest oczekiwane (oryginalna kwestia)
        # Empty source_of_status z level >= L2 jest podejrzane
        if r["yaml"]["level"] in ("L2", "L3") and not r["yaml"]["source_of_status"]:
            issues.append(f"- ❌ `{r['folder']}` ma `level: {r['yaml']['level']}` ale puste `source_of_status`")
        # promoted_to_core != null (powinno być null w Sesji 4)
        if r["yaml"]["promoted_to_core"] is not None:
            issues.append(f"- ❌ `{r['folder']}` ma niezerowe `promoted_to_core` w Sesji 4 (zakazane)")

    if issues:
        for issue in issues:
            lines.append(issue)
    else:
        lines.append("- ✅ **0 CRITICAL findings** — wszystkie foldery przechodzą basic anti-overclaim audit.")
        lines.append("")
        lines.append("Twarde reguły potwierdzone:")
        lines.append("- 0 folderów z `level: L4` (zakazane bez INTAKE workflow w Sesji 7)")
        lines.append("- 0 folderów z niepustym `promoted_to_core` w Sesji 4")
        lines.append("- Każdy folder z `level >= L2` ma niepusty `source_of_status`")
    lines.append("")

    # ── Footer ──
    lines.append("## 9. Generation metadata")
    lines.append("")
    lines.append(f"- Generated: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"- Source: `_classify_s4_dryrun.json` (matches applied YAMLs)")
    lines.append(f"- Builder: `_build_status_index.py`")
    lines.append(f"- Folders: {n}")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"OK wrote {OUT} ({OUT.stat().st_size} B, {len(lines)} lines)")


if __name__ == "__main__":
    main()
