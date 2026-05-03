"""
Sesja 5 — generator FINDINGS.md + NEEDS.md per folder.

Reguły anty-overclaim (zob. AGENT_PROTOCOL.md §3):
  - WYŁĄCZNIE cytaty z istniejących plików.
  - ZERO nowych derivacji / claimów.
  - Każdy item ma `source:` cytujący plik + sekcję / linię.
  - Pusty extract → flag `has_findings_file: false` / `has_needs_file: false`.

Heurystyki ekstrakcji:

FINDINGS:
  - Phase{N}_results.md frontmatter (status, verdict, score) → numeric verdicts
  - "TL;DR" / "Wynik" / "Verdict" sections z Phase results → eksportable summary
  - "X/Y PASS" z .txt files → numerical confirmation
  - Boxed formulas ($$...$$, \[...\]) → exportable derivations
  - Tables z "Result" / "Wynik" / "Verdict" columns → structured findings

NEEDS:
  - "Open issues" / "Pytania otwarte" / "TODO" / "BLOCKED" sections
  - phaseN_score: pending implementation patterns
  - OPEN markers w tabelach
  - "wymaga" / "requires" / "Phase X pending" w tekście
  - For polluted: NEEDS reflects critique items
  - For closure_aggregator: NEEDS reflects KNOWN_ISSUES.md (jeśli istnieje)
"""

from __future__ import annotations
import argparse
import json
import re
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
TGP_V1 = HERE.parent.parent
RESEARCH = TGP_V1 / "research"
DRYRUN_OUT = HERE / "_extract_s5_dryrun.json"
SUMMARY_OUT = HERE / "_extract_s5_summary.md"

POLLUTED = {
    "op-chi1-newton-constant-derivation",
    "op-uv2-mtgp-absolute-scale",
    "op-omega2-axion-coupling-lock",
    "op-omega3-axion-decay-constant",
}
LEGACY = {"op6", "op7", "op1-op2-op4"}
SANDBOX_RESERVED = {"_sandbox", "_archive"}
TODAY = datetime.now().strftime("%Y-%m-%d")

# Markers to look for in section headers
FINDINGS_SECTION_RE = re.compile(
    r"^#+\s*(TL;?DR|TLDR|Wynik|Wyniki|Verdict|Final|Status|Score|Highlights|Key Result|Główny wynik|Główne rezultaty|Kluczowe (rezultaty|predykcje|wnioski))",
    re.IGNORECASE | re.MULTILINE,
)
NEEDS_SECTION_RE = re.compile(
    r"^#+\s*(Open|Otwarte|TODO|Pending|Pytania|Outstanding|Blocked|Brakuje|Limits?|Limitations|Caveats|Future)",
    re.IGNORECASE | re.MULTILINE,
)

# Patterns for inline data
PASS_RATIO_RE = re.compile(r"\b(\d+)\s*/\s*(\d+)\s*PASS\b")
SCORE_PASS_RE = re.compile(r"\b(\d+\.?\d*)\s*/\s*(\d+\.?\d*)\s*(PASS|GATE)\b")
BOXED_RE = re.compile(r"\$\$([^$]{10,300}?)\$\$", re.DOTALL)
PHASE_PENDING_RE = re.compile(
    r"phase\d_score:\s*pending|phaseN_score:\s*pending|Phase\s*\d\s*pending|Phase\s*\d\s*requires|out of scope|TODO|BLOCKED|wymaga\s+\w",
    re.IGNORECASE,
)
OPEN_MARKER_RE = re.compile(r"^\s*\|\s*[^|]*\|\s*\*?\*?OPEN", re.MULTILINE)


def read_text(p: Path, max_bytes: int = 800_000) -> str:
    try:
        if p.stat().st_size > max_bytes:
            with p.open("rb") as f:
                return f.read(max_bytes).decode("utf-8", errors="ignore")
        return p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""


def extract_section(text: str, section_re, max_chars: int = 800) -> list[tuple[str, str]]:
    """
    Find sections matching section_re in text, return list of (heading, body[:max_chars]).
    """
    results = []
    matches = list(section_re.finditer(text))
    for i, m in enumerate(matches):
        start = m.start()
        end = matches[i + 1].start() if i + 1 < len(matches) else min(len(text), start + max_chars * 3)
        # cut to next ## of same or higher level
        heading_match = re.match(r"^(#+)\s*(.*?)$", text[start:start + 200], re.MULTILINE)
        if heading_match:
            level = len(heading_match.group(1))
            heading_title = heading_match.group(2).strip()
            # find next heading of same or higher
            after_heading = start + heading_match.end()
            search_text = text[after_heading:end]
            same_or_higher_re = re.compile(rf"^(#{{{1},{level}}})\s+", re.MULTILINE)
            next_heading = same_or_higher_re.search(search_text)
            if next_heading:
                end = after_heading + next_heading.start()
            body = text[start:end]
            if len(body) > max_chars:
                body = body[:max_chars] + "...(truncated)"
            results.append((heading_title, body.strip()))
    return results


def extract_pass_counts(text: str) -> list[str]:
    results = []
    for m in PASS_RATIO_RE.finditer(text):
        ctx_start = max(0, m.start() - 80)
        ctx_end = min(len(text), m.end() + 80)
        ctx = text[ctx_start:ctx_end].replace("\n", " ").strip()
        results.append(f"{m.group(0)} — context: …{ctx}…")
    return results[:10]


def extract_boxed_formulas(text: str) -> list[str]:
    return [m.group(1).strip()[:200] for m in BOXED_RE.finditer(text)][:5]


def yaml_frontmatter_field(text: str, field: str) -> str | None:
    """Extract a single field value from YAML frontmatter."""
    if not text.lstrip().startswith("---"):
        return None
    parts = text.split("---", 2)
    if len(parts) < 3:
        return None
    fm = parts[1]
    m = re.search(rf"^{re.escape(field)}:\s*(.+?)$", fm, re.MULTILINE)
    if m:
        return m.group(1).strip().strip('"').strip("'")
    return None


def scan_for_findings(folder: Path) -> dict:
    """Extract findings (cited verdicts/results) from existing files."""
    findings = {
        "phase_results": [],     # list[(phase_name, status, verdict, score, file)]
        "tldr_sections": [],     # list[(file, heading, body)]
        "pass_counts": [],       # list[(file, "X/Y PASS context")]
        "boxed_formulas": [],    # list[(file, formula)]
        "readme_status_section": None,  # ("Stan"/"Status" body z README)
        "files_scanned": [],
    }

    # Check Phase*_results.md
    for f in folder.glob("Phase*_results.md"):
        text = read_text(f)
        findings["files_scanned"].append(f.name)
        status = yaml_frontmatter_field(text, "status")
        verdict = yaml_frontmatter_field(text, "verdict")
        score = yaml_frontmatter_field(text, "overall_score") or yaml_frontmatter_field(text, "score")
        program_status = yaml_frontmatter_field(text, "program_status")
        cycle = yaml_frontmatter_field(text, "cycle")
        findings["phase_results"].append({
            "file": f.name,
            "cycle": cycle,
            "status": status,
            "verdict": verdict,
            "score": score,
            "program_status": program_status,
        })
        # TL;DR section
        for heading, body in extract_section(text, FINDINGS_SECTION_RE, max_chars=600):
            findings["tldr_sections"].append({
                "file": f.name,
                "heading": heading,
                "body": body,
            })

    # Check results.md (closure_aggregator children, but Q2: parent only — so only top-level results.md)
    for fname in ("results.md", "RESULTS.md", "VERDICT.md"):
        f = folder / fname
        if f.exists():
            text = read_text(f)
            findings["files_scanned"].append(f.name)
            for heading, body in extract_section(text, FINDINGS_SECTION_RE, max_chars=600):
                findings["tldr_sections"].append({
                    "file": f.name,
                    "heading": heading,
                    "body": body,
                })

    # Check CLOSURE_*_SUMMARY.md (closure_aggregator)
    for f in folder.glob("CLOSURE_*_SUMMARY.md"):
        text = read_text(f)
        findings["files_scanned"].append(f.name)
        for heading, body in extract_section(text, FINDINGS_SECTION_RE, max_chars=800):
            findings["tldr_sections"].append({
                "file": f.name,
                "heading": heading,
                "body": body,
            })

    # README "Status" / "Stan" section
    readme = folder / "README.md"
    if readme.exists():
        text = read_text(readme)
        for heading, body in extract_section(text, FINDINGS_SECTION_RE, max_chars=400):
            findings["tldr_sections"].append({
                "file": "README.md",
                "heading": heading,
                "body": body,
            })

    # PASS counts from .txt files
    for f in folder.glob("*.txt"):
        text = read_text(f, max_bytes=200_000)
        for ctx in extract_pass_counts(text):
            findings["pass_counts"].append({"file": f.name, "context": ctx})

    # Boxed formulas (top-level .md)
    for f in folder.glob("*.md"):
        if f.name in ("FINDINGS.md", "NEEDS.md"):
            continue
        text = read_text(f)
        for formula in extract_boxed_formulas(text):
            findings["boxed_formulas"].append({"file": f.name, "formula": formula})

    # Limit volumes
    findings["pass_counts"] = findings["pass_counts"][:8]
    findings["boxed_formulas"] = findings["boxed_formulas"][:5]
    findings["tldr_sections"] = findings["tldr_sections"][:6]
    return findings


def scan_for_needs(folder: Path) -> dict:
    """Extract needs/open-issues (cited gaps) from existing files."""
    needs = {
        "open_sections": [],      # list[(file, heading, body)]
        "pending_phases": [],     # list[(file, "phaseN_score: pending"-context)]
        "open_table_rows": [],    # list[(file, row)]
        "wymaga_sentences": [],   # list[(file, sentence)]
        "known_issues": [],       # KNOWN_ISSUES.md OPEN items (closure_aggregator)
        "files_scanned": [],
    }

    files_to_scan = []
    files_to_scan.extend(folder.glob("Phase*_results.md"))
    files_to_scan.extend(folder.glob("Phase*_setup.md"))
    files_to_scan.extend(folder.glob("program.md"))
    files_to_scan.extend(folder.glob("PLAN*.md"))
    files_to_scan.extend(folder.glob("README.md"))
    files_to_scan.extend(folder.glob("results.md"))
    files_to_scan.extend(folder.glob("KNOWN_ISSUES.md"))
    files_to_scan.extend(folder.glob("SUMMARY.md"))

    for f in files_to_scan:
        if not f.exists():
            continue
        text = read_text(f)
        needs["files_scanned"].append(f.name)

        # Open / TODO / Pending sections
        for heading, body in extract_section(text, NEEDS_SECTION_RE, max_chars=500):
            needs["open_sections"].append({
                "file": f.name,
                "heading": heading,
                "body": body,
            })

        # Phase N pending markers
        for m in PHASE_PENDING_RE.finditer(text):
            ctx_start = max(0, m.start() - 60)
            ctx_end = min(len(text), m.end() + 60)
            ctx = text[ctx_start:ctx_end].replace("\n", " ").strip()
            needs["pending_phases"].append({"file": f.name, "context": ctx})

        # OPEN table rows
        for m in OPEN_MARKER_RE.finditer(text):
            row_start = m.start()
            row_end = text.find("\n", row_start)
            if row_end == -1:
                row_end = min(len(text), row_start + 200)
            row = text[row_start:row_end].strip()
            if row:
                needs["open_table_rows"].append({"file": f.name, "row": row[:200]})

    # Limit volumes
    needs["open_sections"] = needs["open_sections"][:6]
    needs["pending_phases"] = needs["pending_phases"][:8]
    needs["open_table_rows"] = needs["open_table_rows"][:8]
    return needs


def sanitize_banned_phrases(text: str) -> str:
    """
    Filter banned phrases (per AGENT_PROTOCOL §3) from cite-only extracted content.

    Strategy: replace "FULL CONVERGENCE" / "FULL_CONVERGENCE" with cite-quote
    indicating phrase redaction. PRESERVE adjacent PASS-count if present
    (so numerical info isn't lost).

    Examples:
      "6/6 PASS — FULL CONVERGENCE" → "6/6 PASS [phrase redacted per AGENT_PROTOCOL §3]"
      "FULL CONVERGENCE 16/16:" → "16/16 PASS [phrase redacted per AGENT_PROTOCOL §3]:"
      "tags: [FULL_CONVERGENCE]" → "tags: [REDACTED-PER-PROTOCOL]"
    """
    # Pattern 1: "X/Y PASS — FULL CONVERGENCE" or similar (PASS-count before phrase)
    text = re.sub(
        r"(\d+/\d+\s*PASS)\s*[—–-]\s*FULL CONVERGENCE\b\.?",
        r"\1 [phrase redacted per AGENT_PROTOCOL §3]",
        text,
    )
    # Pattern 2: "FULL CONVERGENCE X/Y" (phrase before PASS-count)
    text = re.sub(
        r"\bFULL CONVERGENCE\s+(\d+/\d+)\b",
        r"\1 PASS [phrase redacted per AGENT_PROTOCOL §3]",
        text,
    )
    # Pattern 3: bare "FULL CONVERGENCE" without nearby PASS count
    text = re.sub(
        r"\bFULL CONVERGENCE\b",
        "[phrase redacted per AGENT_PROTOCOL §3]",
        text,
    )
    # Pattern 4: YAML tag "FULL_CONVERGENCE"
    text = re.sub(
        r"\bFULL_CONVERGENCE\b",
        "REDACTED-PER-PROTOCOL",
        text,
    )
    return text


def build_findings_md(folder: Path, info: dict, status_yaml: dict, findings_data: dict) -> str:
    """Build FINDINGS.md content."""
    name = folder.name
    title = info.get("readme_h1") or name

    lines = []
    lines.append("---")
    lines.append(f'title: "FINDINGS — {title}"')
    lines.append(f"date: {TODAY}")
    lines.append('parent: "[[README.md]]"')
    lines.append("type: findings")
    lines.append(f"tgp_owner: research/{name}")
    lines.append("source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)")
    lines.append("tags:")
    lines.append("  - findings")
    lines.append("---")
    lines.append("")
    lines.append(f"# FINDINGS — {title}")
    lines.append("")
    lines.append("> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących")
    lines.append("> plików folderu. **Żadne treści nie są wymyślone** — każdy item ma `source:`")
    lines.append("> cytujący plik. Manual review w Sesji 6 (RESEARCH_BUS broadcasts).")
    lines.append("")

    # POLLUTED case — no findings, just quarantine note
    if name in POLLUTED:
        lines.append("## Status folderu: POLLUTED-74394a8 (kwarantanna)")
        lines.append("")
        lines.append("Folder w kwarantannie. **Findings są zamrożone** — claimy w plikach")
        lines.append("źródłowych są niezweryfikowane (zob. `CRITIQUE_*.md` i")
        lines.append("[[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]]).")
        lines.append("")
        lines.append("**`exports_findings: false`** — folder nie eksportuje wyników do")
        lines.append("`RESEARCH_BUS.md` do czasu forward-patch.")
        lines.append("")
        # Cite the critique file if exists
        for f in folder.glob("CRITIQUE_*.md"):
            lines.append(f"## Self-critique zachowane")
            lines.append("")
            lines.append(f"- Plik: `{f.name}` ({f.stat().st_size} B)")
            text = read_text(f)
            verdict = yaml_frontmatter_field(text, "Verdict") or yaml_frontmatter_field(text, "verdict")
            if verdict:
                lines.append(f"- Verdict (z frontmatter): {verdict}")
            # First sentence after the first ## heading
            m = re.search(r"^##\s+.*?\n\n([^\n]{50,400})", text, re.MULTILINE)
            if m:
                lines.append(f"- Cytat: {m.group(1).strip()}")
            lines.append("")
        lines.append("---")
        lines.append("")
        lines.append("*Folder **nie eksportuje findings** do RESEARCH_BUS. Zob. AGENT_PROTOCOL §3.2 quarantine clause.*")
        return "\n".join(lines)

    # Phase results verdicts (top section for op-* folders)
    if findings_data["phase_results"]:
        lines.append("## Phase results — frontmatter verdicts")
        lines.append("")
        lines.append("| Plik | Cycle | Status | Verdict | Score | Program |")
        lines.append("|------|-------|--------|---------|-------|---------|")
        for pr in findings_data["phase_results"]:
            cells = [pr["file"], pr["cycle"], pr["status"], pr["verdict"], pr["score"], pr["program_status"]]
            cells = [str(c) if c is not None else "—" for c in cells]
            lines.append("| " + " | ".join(f"`{c}`" if c != "—" else c for c in cells) + " |")
        lines.append("")

    # Numerical PASS counts
    if findings_data["pass_counts"]:
        lines.append("## Numerical PASS counts (cited from .txt outputs)")
        lines.append("")
        for pc in findings_data["pass_counts"]:
            ctx = pc["context"][:180].replace("|", "\\|")
            lines.append(f"- `{pc['file']}` — {ctx}")
        lines.append("")

    # TL;DR / Wynik sections
    if findings_data["tldr_sections"]:
        lines.append("## TL;DR / Wynik / Verdict sections (cytaty)")
        lines.append("")
        for tldr in findings_data["tldr_sections"]:
            lines.append(f"### `{tldr['file']}` — {tldr['heading']}")
            lines.append("")
            # Quote the body, escaping markdown breakage
            body = tldr["body"]
            # Indent as blockquote
            for line in body.splitlines():
                lines.append(f"> {line}")
            lines.append("")

    # Boxed formulas
    if findings_data["boxed_formulas"]:
        lines.append("## Eksportowalne formuły (boxed)")
        lines.append("")
        for bf in findings_data["boxed_formulas"]:
            lines.append(f"**`{bf['file']}`:**")
            lines.append("")
            lines.append("$$")
            lines.append(bf["formula"])
            lines.append("$$")
            lines.append("")

    # If empty
    if not (findings_data["phase_results"] or findings_data["pass_counts"] or
            findings_data["tldr_sections"] or findings_data["boxed_formulas"]):
        lines.append("## (Brak ekstrahowalnych findings)")
        lines.append("")
        lines.append("Folder nie ma plików z weryfikowalnymi PASS-ami / formułami / verdict-section.")
        lines.append("Status: `exports_findings: false`. Sesja 6 może wykryć to jako orphan; Sesja 8")
        lines.append("regression flag.")
        lines.append("")
        if findings_data["files_scanned"]:
            lines.append(f"Pliki przeszukane: {', '.join(f'`{f}`' for f in findings_data['files_scanned'])}")
        lines.append("")

    # Footer
    lines.append("---")
    lines.append("")
    lines.append("## Cross-references")
    lines.append("")
    lines.append("- [[README.md]] — opis folderu + YAML status")
    lines.append("- [[NEEDS.md]] — otwarte luki tego folderu")
    lines.append("- [[meta/research/RESEARCH_BUS.md]] — broadcast tych findings")
    lines.append("- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa")

    return "\n".join(lines)


def build_needs_md(folder: Path, info: dict, status_yaml: dict, needs_data: dict) -> str:
    """Build NEEDS.md content."""
    name = folder.name
    title = info.get("readme_h1") or name

    lines = []
    lines.append("---")
    lines.append(f'title: "NEEDS — {title}"')
    lines.append(f"date: {TODAY}")
    lines.append('parent: "[[README.md]]"')
    lines.append("type: needs")
    lines.append(f"tgp_owner: research/{name}")
    lines.append("source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)")
    lines.append("tags:")
    lines.append("  - needs")
    lines.append("---")
    lines.append("")
    lines.append(f"# NEEDS — {title}")
    lines.append("")
    lines.append("> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących")
    lines.append("> plików folderu. Każdy item ma `source:` cytujący plik.")
    lines.append("")

    # POLLUTED case
    if name in POLLUTED:
        lines.append("## Open: kwarantanna 74394a8")
        lines.append("")
        lines.append("| ID | Luka | Source | Status |")
        lines.append("|----|------|--------|--------|")
        lines.append(f"| N1 | Forward-patch decyzja człowieka po `{name}` polluted commit | [[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]] | OPEN |")
        # Cite critique
        for f in folder.glob("CRITIQUE_*.md"):
            text = read_text(f)
            verdict = yaml_frontmatter_field(text, "Verdict") or yaml_frontmatter_field(text, "verdict") or "BLOCKING"
            lines.append(f"| N2 | Self-critique unresolved | `{f.name}` (Verdict: {verdict}) | OPEN |")
        lines.append("")
        lines.append("## Reguła kwarantanny")
        lines.append("")
        lines.append("Per AGENT_PROTOCOL §3.2: agent **NIE WOLNO**:")
        lines.append("- awansować ten folder do `core-ready` lub `core-promoted`")
        lines.append("- dodawać `INTAKE_*.md` dla tego folderu")
        lines.append("- zmieniać statusu na `archive` (zatuszowałoby ślad)")
        lines.append("")
        lines.append("---")
        lines.append("")
        lines.append("## Cross-references")
        lines.append("")
        lines.append("- [[README.md]]")
        lines.append("- [[FINDINGS.md]] — kwarantannowane (frozen)")
        lines.append("- [[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]]")
        return "\n".join(lines)

    # Open sections
    if needs_data["open_sections"]:
        lines.append("## Otwarte sekcje (cytaty z plików)")
        lines.append("")
        for sec in needs_data["open_sections"]:
            lines.append(f"### `{sec['file']}` — {sec['heading']}")
            lines.append("")
            for line in sec["body"].splitlines():
                lines.append(f"> {line}")
            lines.append("")

    # Pending phases
    if needs_data["pending_phases"]:
        lines.append("## Phase pending / requirements (cytaty)")
        lines.append("")
        for pp in needs_data["pending_phases"][:6]:
            ctx = pp["context"][:200].replace("|", "\\|")
            lines.append(f"- `{pp['file']}` — …{ctx}…")
        lines.append("")

    # OPEN table rows
    if needs_data["open_table_rows"]:
        lines.append("## OPEN markery w tabelach")
        lines.append("")
        for row in needs_data["open_table_rows"][:6]:
            cleaned_row = row['row'].replace("|", " | ")
            lines.append(f"- `{row['file']}`: `{cleaned_row}`")
        lines.append("")

    # If empty
    if not (needs_data["open_sections"] or needs_data["pending_phases"] or needs_data["open_table_rows"]):
        lines.append("## (Brak ekstrahowalnych otwartych luk)")
        lines.append("")
        lines.append("Folder nie ma w plikach explicite oznaczonych otwartych luk")
        lines.append("(brak sekcji 'Open' / 'Pending' / 'TODO' / OPEN table rows).")
        lines.append("To może oznaczać: (a) folder jest faktycznie domknięty,")
        lines.append("(b) luki są implicite, niezformalizowane, (c) wymaga manualnego review.")
        lines.append("")
        if needs_data["files_scanned"]:
            lines.append(f"Pliki przeszukane: {', '.join(f'`{f}`' for f in needs_data['files_scanned'])}")
        lines.append("")

    # Footer
    lines.append("---")
    lines.append("")
    lines.append("## Cross-references")
    lines.append("")
    lines.append("- [[README.md]] — opis folderu + YAML status")
    lines.append("- [[FINDINGS.md]] — eksportowalne wyniki tego folderu")
    lines.append("- [[meta/research/CANDIDATE_BRIDGES.md]] — match z `FINDINGS.md` innych folderów")
    lines.append("- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa")

    return "\n".join(lines)


def yaml_frontmatter_field_from_readme(folder: Path, field: str) -> str | None:
    readme = folder / "README.md"
    if not readme.exists():
        return None
    text = read_text(readme)
    return yaml_frontmatter_field(text, field)


def collect_folder_info(folder: Path) -> dict:
    readme = folder / "README.md"
    info = {"name": folder.name, "readme_h1": None}
    if readme.exists():
        text = read_text(readme)
        for line in text.splitlines():
            ls = line.strip()
            if ls.startswith("# "):
                info["readme_h1"] = ls[2:].strip()
                break
    return info


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true",
                        help="Apply: write FINDINGS.md and NEEDS.md (default: dry-run)")
    parser.add_argument("--folder", type=str, default=None,
                        help="Single folder (debugging)")
    args = parser.parse_args()

    folders = []
    for p in sorted(RESEARCH.iterdir(), key=lambda x: x.name.lower()):
        if not p.is_dir() or p.name in SANDBOX_RESERVED:
            continue
        if args.folder and p.name != args.folder:
            continue
        folders.append(p)

    print(f"Processing {len(folders)} folders (apply={args.apply})...")

    results = []
    for folder in folders:
        info = collect_folder_info(folder)
        findings_data = scan_for_findings(folder)
        needs_data = scan_for_needs(folder)

        findings_md = build_findings_md(folder, info, {}, findings_data)
        needs_md = build_needs_md(folder, info, {}, needs_data)

        # Sanitize banned phrases (post-extraction filter, per AGENT_PROTOCOL §3)
        findings_md = sanitize_banned_phrases(findings_md)
        needs_md = sanitize_banned_phrases(needs_md)

        # "Empty" detection
        findings_empty = ("(Brak ekstrahowalnych findings)" in findings_md)
        needs_empty = ("(Brak ekstrahowalnych otwartych luk)" in needs_md)

        results.append({
            "folder": folder.name,
            "findings_data": findings_data,
            "needs_data": needs_data,
            "findings_size": len(findings_md),
            "needs_size": len(needs_md),
            "findings_empty": findings_empty,
            "needs_empty": needs_empty,
            "findings_preview": findings_md[:500],
            "needs_preview": needs_md[:500],
        })

        if args.apply:
            (folder / "FINDINGS.md").write_text(findings_md, encoding="utf-8")
            (folder / "NEEDS.md").write_text(needs_md, encoding="utf-8")

    # JSON output (debugging)
    summary_data = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "n_folders": len(results),
        "applied": args.apply,
        "n_findings_empty": sum(1 for r in results if r["findings_empty"]),
        "n_needs_empty": sum(1 for r in results if r["needs_empty"]),
        "results": [{k: v for k, v in r.items() if k != "findings_data" and k != "needs_data"} for r in results],
    }
    DRYRUN_OUT.write_text(json.dumps(summary_data, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print(f"Wrote {DRYRUN_OUT} ({DRYRUN_OUT.stat().st_size} B)")

    # Markdown summary
    lines = []
    lines.append(f"# Sesja 5 — extract dry-run summary")
    lines.append("")
    lines.append(f"Generated: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Folders: {len(results)} | Applied: {args.apply}")
    lines.append("")
    lines.append(f"- FINDINGS empty (no extractable): **{summary_data['n_findings_empty']}**")
    lines.append(f"- NEEDS empty (no extractable): **{summary_data['n_needs_empty']}**")
    lines.append("")
    lines.append("## Per-folder size table")
    lines.append("")
    lines.append("| # | Folder | FINDINGS B | NEEDS B | FINDINGS empty | NEEDS empty |")
    lines.append("|---|--------|-----------:|--------:|:--------------:|:-----------:|")
    for i, r in enumerate(results, 1):
        f_emp = "⚠" if r["findings_empty"] else "✓"
        n_emp = "⚠" if r["needs_empty"] else "✓"
        lines.append(f"| {i} | `{r['folder']}` | {r['findings_size']} | {r['needs_size']} | {f_emp} | {n_emp} |")
    lines.append("")
    SUMMARY_OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {SUMMARY_OUT} ({SUMMARY_OUT.stat().st_size} B)")


if __name__ == "__main__":
    main()
