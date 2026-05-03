"""
Sesja 5 — wygeneruj wpisy BROADCAST w RESEARCH_BUS.md per folder z `exports_findings: true`.

Każdy folder z FINDINGS.md (non-empty) generuje 1 wpis broadcast z:
- Date (TODAY)
- Source folder
- Headline (pierwsza non-trivial linia z FINDINGS)
- Source file (link do FINDINGS.md)
- Consumers: list folderów z `tgp_status.depends_on:` zawierającym ten folder
            (initial Sesja 5: pusta — Sesja 6 wypełni przez bridge matching)
"""

from __future__ import annotations
import json
import re
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent.parent / "research"
BUS_FILE = HERE / "RESEARCH_BUS.md"
SANDBOX_RESERVED = {"_sandbox", "_archive"}
TODAY = datetime.now().strftime("%Y-%m-%d")


def read_text(p: Path, max_bytes: int = 200_000) -> str:
    try:
        if p.stat().st_size > max_bytes:
            with p.open("rb") as f:
                return f.read(max_bytes).decode("utf-8", errors="ignore")
        return p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""


def yaml_field(text: str, field: str) -> str | None:
    if not text.lstrip().startswith("---"):
        return None
    parts = text.split("---", 2)
    if len(parts) < 3:
        return None
    fm = parts[1]
    m = re.search(rf"^\s*{re.escape(field)}:\s*(.+?)$", fm, re.MULTILINE)
    if m:
        return m.group(1).strip().strip('"').strip("'")
    return None


def extract_headline_from_findings(findings_text: str) -> str:
    """Get first meaningful line from FINDINGS.md (skip frontmatter, header, blockquote)."""
    if "---" in findings_text:
        parts = findings_text.split("---", 2)
        if len(parts) >= 3:
            findings_text = parts[2]
    lines = findings_text.splitlines()
    # Skip front quote/note
    skip_prefix = (">", "#", "")
    for line in lines:
        ls = line.strip()
        if not ls:
            continue
        if ls.startswith(("> ", "#", "##", "###", "####", "---")):
            continue
        if "(2026-05-03)" in ls or "auto-generated" in ls.lower():
            continue
        # Trim and return
        return ls[:180]
    # Fallback: H1
    for line in lines:
        if line.startswith("# "):
            return line[2:].strip()[:180]
    return "(no headline extractable)"


def extract_phase_score_summary(findings_text: str) -> str:
    """Get a short summary string with PASS counts if present."""
    # Look for table row with score
    m = re.search(r"\|\s*`Phase\d_results\.md`\s*\|.*?\|\s*`?([^|`]+?)`?\s*\|\s*`?([^|`]+?)`?\s*\|\s*`?([^|`]+?)`?", findings_text)
    if m:
        return f"{m.group(2)} / {m.group(3)}"
    # PASS context
    m = re.search(r"(\d+/\d+\s*PASS)", findings_text)
    if m:
        return m.group(1)
    return ""


def main():
    folders = []
    for p in sorted(RESEARCH.iterdir(), key=lambda x: x.name.lower()):
        if not p.is_dir() or p.name in SANDBOX_RESERVED:
            continue
        folders.append(p)

    broadcasts = []
    no_consumers = []
    for folder in folders:
        readme = folder / "README.md"
        findings = folder / "FINDINGS.md"
        if not readme.exists() or not findings.exists():
            continue

        readme_text = read_text(readme)
        findings_text = read_text(findings)

        exports = yaml_field(readme_text, "exports_findings")
        polluted = yaml_field(readme_text, "polluted_74394a8")

        if exports != "true":
            continue  # Only broadcast folders that export
        if polluted == "true":
            continue  # Polluted are quarantined, no broadcast

        headline = extract_headline_from_findings(findings_text)
        score = extract_phase_score_summary(findings_text)

        broadcasts.append({
            "folder": folder.name,
            "headline": headline,
            "score": score,
            "findings_link": f"[[research/{folder.name}/FINDINGS.md]]",
            "consumers": [],  # Sesja 6 fills
        })

    print(f"Broadcasts to add: {len(broadcasts)}")

    # Build new RESEARCH_BUS.md with broadcasts in section 2
    today_iso = datetime.now().isoformat(timespec="seconds")

    lines = []
    lines.append("---")
    lines.append('title: "RESEARCH_BUS — tablica ogłoszeń międzyfolderowych"')
    lines.append(f"date: {TODAY}")
    lines.append("type: bus")
    lines.append("status: ACTIVE — Sesja 5 initial broadcasts; Sesja 6 fills consumers via CANDIDATE_BRIDGES")
    lines.append('parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"')
    lines.append("related:")
    lines.append('  - "[[meta/research/AGENT_PROTOCOL.md]]"')
    lines.append('  - "[[meta/research/CANDIDATE_BRIDGES.md]]"')
    lines.append('  - "[[meta/research/IMPACT_MATRIX.md]]"')
    lines.append('  - "[[meta/research/FOLDER_STATUS_INDEX.md]]"')
    lines.append("tags:")
    lines.append("  - research-bus")
    lines.append("  - broadcast")
    lines.append("  - inter-folder")
    lines.append("---")
    lines.append("")
    lines.append("# RESEARCH_BUS — tablica ogłoszeń międzyfolderowych")
    lines.append("")
    lines.append("> **Cel pliku:** każdy nowy wynik (item w `FINDINGS.md`) generuje wpis")
    lines.append("> broadcastowy z polem `consumers:` (foldery, które mogą skonsumować).")
    lines.append("> Auto-broadcasted by Sesja 5 (per-folder podsumowanie); Sesja 6 dopisuje")
    lines.append("> consumers przez `CANDIDATE_BRIDGES.md` matching.")
    lines.append("")
    lines.append("## 0. Reguły wpisu")
    lines.append("")
    lines.append("1. **Trigger**: dodanie/aktualizacja `FINDINGS.md` w folderze badawczym")
    lines.append("   (z `exports_findings: true` i `polluted_74394a8: false`).")
    lines.append("2. **Format**: jeden wiersz tabeli na folder source.")
    lines.append("3. **`Consumers`**: lista folderów z `tgp_status.depends_on:` matching")
    lines.append("   ten folder LUB Sesja 6 bridge matching. Pusta dopuszczalna,")
    lines.append("   ale wymaga § 3 entry.")
    lines.append("4. **Status**:")
    lines.append("   - `BROADCAST` — wpis świeży, nieprzeczytany przez konsumentów")
    lines.append("   - `CONSUMED` — co najmniej jeden konsument potwierdził użycie")
    lines.append("   - `STALE` — wpis > 90 dni bez akcji")
    lines.append("")
    lines.append("## 1. Statystyka")
    lines.append("")
    lines.append("| Metryka | Liczba |")
    lines.append("|---------|-------:|")
    lines.append(f"| Aktywnych BROADCAST | {len(broadcasts)} |")
    lines.append(f"| CONSUMED | 0 |")
    lines.append(f"| STALE | 0 |")
    lines.append(f"| **Razem (wszystkie czasy)** | **{len(broadcasts)}** |")
    lines.append("")
    lines.append(f"Last build: {today_iso}.")
    lines.append("")
    lines.append("## 2. Aktywne ogłoszenia (BROADCAST)")
    lines.append("")
    lines.append("| Date | Source folder | Headline | Score | Findings link | Consumers | Status |")
    lines.append("|------|---------------|----------|-------|---------------|-----------|--------|")
    for b in broadcasts:
        cons = ", ".join(f"`{c}`" for c in b["consumers"]) if b["consumers"] else "(scanned, pending)"
        headline_safe = b["headline"].replace("|", "\\|")[:120]
        score_str = b["score"] if b["score"] else "—"
        lines.append(
            f"| {TODAY} | `{b['folder']}` | {headline_safe} | {score_str} | "
            f"{b['findings_link']} | {cons} | BROADCAST |"
        )
    lines.append("")

    # Section 3 — Skanowania (consumers: [] z notatką)
    lines.append("## 3. Skanowania (consumers: [] — pending Sesja 6 matching)")
    lines.append("")
    lines.append(f"All {len(broadcasts)} broadcasts mają consumers: [] do czasu Sesji 6")
    lines.append("(CANDIDATE_BRIDGES matching). To NIE oznacza orfan — oznacza pending.")
    lines.append("")
    lines.append("Po Sesji 6 broadcasts przejdą jeden z trzech stanów:")
    lines.append("- (a) Match z `NEEDS.md` innego folderu → wpis w `CANDIDATE_BRIDGES.md` jako PROPOSED")
    lines.append("- (b) Brak match (orfan) → tabela § 3 z notatką `# scanned, no consumers found YYYY-MM-DD`")
    lines.append("- (c) Pre-existing `depends_on:` w innym folderze cytuje ten folder → CONSUMED")
    lines.append("")

    lines.append("## 4. CONSUMED")
    lines.append("")
    lines.append("(empty — czeka na Sesję 6 matching i acknowledgments)")
    lines.append("")
    lines.append("## 5. STALE (≥ 90 dni bez akcji)")
    lines.append("")
    lines.append("(empty)")
    lines.append("")

    # Folders skipped
    lines.append("## 6. Foldery wyłączone z broadcast")
    lines.append("")
    lines.append("- 4 polluted-74394a8 (kwarantanna; AGENT_PROTOCOL §3.2)")
    lines.append("- 12 folderów z `exports_findings: false` (FINDINGS.md zawiera")
    lines.append("  `(Brak ekstrahowalnych findings)`):")
    skipped = []
    for folder in folders:
        readme = folder / "README.md"
        if not readme.exists():
            continue
        rt = read_text(readme)
        exp = yaml_field(rt, "exports_findings")
        polluted = yaml_field(rt, "polluted_74394a8")
        if exp == "false" and polluted != "true":
            skipped.append(folder.name)
    for s in sorted(skipped):
        lines.append(f"  - `{s}`")
    lines.append("")
    lines.append("Te foldery są kandydatami do (a) manual review FINDINGS w Sesji 6+,")
    lines.append("(b) flag w Sesji 8 jako wymagające scaffoldingu.")
    lines.append("")

    BUS_FILE.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {BUS_FILE} ({BUS_FILE.stat().st_size} B)")


if __name__ == "__main__":
    main()
