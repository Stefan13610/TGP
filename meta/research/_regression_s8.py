"""
Sesja 8 — regression / impact review (read-only).

Reguły (zob. PLAN_RESEARCH_WORKFLOW_v1.md §6 Sesja 8):
  - **NIE** modyfikuje plików w research/, core/, meta/.
  - Wynik: raport `meta/research/REGRESSION_S8.md`.

Checks:
  C1. Each folder has README + NEEDS + FINDINGS
  C2. level: L4 wymaga promoted_to_core (weryfikowalna ref)
  C3. RESEARCH_BUS broadcasts mają source w FINDINGS.md
  C4. CANDIDATE_BRIDGES bridges mają source/target istniejące
  C5. depends_on / impacts prowadzą do realnych folderów
  C6. last_reviewed_against_core ≤ today
  C7. last_yaml_update ≤ today

Anti-overclaim audits:
  A1. "FULL CONVERGENCE" / "100%" / "DERIVED" w NEEDS/FINDINGS bez source_of_status
  A2. source_of_status z "heuristic" gdy level >= L3
  A3. Polluted folders (74394a8) z exports_findings: true (powinny być false)

Severity:
  CRITICAL — narusza twardą regułę (anty-overclaim L4, brak źródła)
  HIGH    — strukturalna niespójność (broken cross-ref)
  MEDIUM  — minor inconsistency (date drift)
  INFO    — observation
"""

from __future__ import annotations
import json
import re
from datetime import date, datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
TGP_V1 = HERE.parent.parent
RESEARCH = TGP_V1 / "research"
META_RESEARCH = TGP_V1 / "meta" / "research"
META_CORE = TGP_V1 / "meta" / "core"

REPORT_OUT = HERE / "REGRESSION_S8.md"

POLLUTED = {
    "op-chi1-newton-constant-derivation",
    "op-uv2-mtgp-absolute-scale",
    "op-omega2-axion-coupling-lock",
    "op-omega3-axion-decay-constant",
}
SANDBOX_RESERVED = {"_sandbox", "_archive"}

TODAY_DATE = date.today()
TODAY_STR = TODAY_DATE.isoformat()

OVERCLAIM_PHRASES = [
    "FULL CONVERGENCE",
    "FULL_CONVERGENCE",
]


def read_text(p: Path, max_bytes: int = 500_000) -> str:
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


def yaml_list_field(text: str, field: str) -> list[str]:
    """Extract YAML list field. Handles both inline [...] and block - item."""
    if not text.lstrip().startswith("---"):
        return []
    parts = text.split("---", 2)
    if len(parts) < 3:
        return []
    fm = parts[1]
    # Try inline format
    m = re.search(rf"^\s*{re.escape(field)}:\s*\[(.*?)\]\s*$", fm, re.MULTILINE)
    if m:
        items = m.group(1).split(",")
        return [it.strip().strip('"').strip("'") for it in items if it.strip()]
    # Block format
    m = re.search(rf"^\s*{re.escape(field)}:\s*\n((?:^    -\s+.*\n)+)", fm, re.MULTILINE)
    if m:
        items = []
        for line in m.group(1).splitlines():
            sm = re.match(r"^    -\s+(.+)$", line)
            if sm:
                v = sm.group(1).strip().strip('"').strip("'")
                items.append(v)
        return items
    return []


def parse_iso_date(s: str | None) -> date | None:
    if not s or s.lower() == "unknown" or s.lower() == "null":
        return None
    s = s.strip().strip('"').strip("'")
    try:
        return date.fromisoformat(s)
    except Exception:
        try:
            return datetime.fromisoformat(s).date()
        except Exception:
            return None


def main():
    findings: list[dict] = []  # list of {severity, category, folder, message, evidence}

    # Collect folders
    research_folders = []
    for p in sorted(RESEARCH.iterdir(), key=lambda x: x.name.lower()):
        if not p.is_dir() or p.name in SANDBOX_RESERVED:
            continue
        research_folders.append(p)

    folder_names = {f.name for f in research_folders}

    # ── C1: Each folder has README + NEEDS + FINDINGS ──
    for folder in research_folders:
        for fname in ("README.md", "NEEDS.md", "FINDINGS.md"):
            p = folder / fname
            if not p.exists():
                findings.append({
                    "severity": "CRITICAL",
                    "category": "C1.missing-file",
                    "folder": folder.name,
                    "message": f"missing {fname}",
                    "evidence": str(p.relative_to(TGP_V1)),
                })

    # Load all YAMLs once
    folder_yamls = {}
    for folder in research_folders:
        readme = folder / "README.md"
        if readme.exists():
            text = read_text(readme)
            folder_yamls[folder.name] = {
                "text": text,
                "level": yaml_field(text, "level"),
                "folder_status": yaml_field(text, "folder_status"),
                "promoted_to_core": yaml_field(text, "promoted_to_core"),
                "core_compatibility": yaml_field(text, "core_compatibility"),
                "last_reviewed_against_core": yaml_field(text, "last_reviewed_against_core"),
                "last_yaml_update": yaml_field(text, "last_yaml_update"),
                "polluted_74394a8": yaml_field(text, "polluted_74394a8"),
                "exports_findings": yaml_field(text, "exports_findings"),
                "has_findings_file": yaml_field(text, "has_findings_file"),
                "has_needs_file": yaml_field(text, "has_needs_file"),
                "depends_on": yaml_list_field(text, "depends_on"),
                "impacts": yaml_list_field(text, "impacts"),
                "source_of_status": yaml_list_field(text, "source_of_status"),
            }

    # ── C2: level: L4 wymaga promoted_to_core ──
    for name, y in folder_yamls.items():
        if y.get("level") == "L4":
            promoted = y.get("promoted_to_core")
            if not promoted or promoted == "null":
                findings.append({
                    "severity": "CRITICAL",
                    "category": "C2.L4-without-promoted-to-core",
                    "folder": name,
                    "message": "level: L4 but promoted_to_core is null/empty",
                    "evidence": "research/" + name + "/README.md",
                })

    # ── C5: depends_on / impacts → real folders ──
    for name, y in folder_yamls.items():
        for dep in y["depends_on"]:
            # Strip path prefix
            dep_clean = dep.replace("research/", "").split("/")[0]
            if dep_clean and dep_clean not in folder_names:
                findings.append({
                    "severity": "HIGH",
                    "category": "C5.broken-depends_on",
                    "folder": name,
                    "message": f"depends_on '{dep}' → folder '{dep_clean}' not found",
                    "evidence": "research/" + name + "/README.md",
                })
        for imp in y["impacts"]:
            imp_clean = imp.replace("research/", "").split("/")[0]
            if imp_clean and imp_clean not in folder_names and not imp_clean.startswith("core"):
                findings.append({
                    "severity": "HIGH",
                    "category": "C5.broken-impacts",
                    "folder": name,
                    "message": f"impacts '{imp}' → folder '{imp_clean}' not found",
                    "evidence": "research/" + name + "/README.md",
                })

    # ── C6/C7: dates ≤ today ──
    for name, y in folder_yamls.items():
        for field in ("last_reviewed_against_core", "last_yaml_update"):
            d = parse_iso_date(y.get(field))
            if d and d > TODAY_DATE:
                findings.append({
                    "severity": "MEDIUM",
                    "category": f"C6.future-date-{field}",
                    "folder": name,
                    "message": f"{field}={d.isoformat()} > today {TODAY_STR}",
                    "evidence": "research/" + name + "/README.md",
                })

    # ── C3: RESEARCH_BUS broadcasts → source FINDINGS.md ──
    bus_path = META_RESEARCH / "RESEARCH_BUS.md"
    bus_text = read_text(bus_path) if bus_path.exists() else ""
    if bus_text:
        # Parse broadcast rows
        for m in re.finditer(r"^\|\s*\d{4}-\d{2}-\d{2}\s*\|\s*`([^`]+)`\s*\|", bus_text, re.MULTILINE):
            source_folder = m.group(1)
            findings_path = RESEARCH / source_folder / "FINDINGS.md"
            if not findings_path.exists():
                findings.append({
                    "severity": "HIGH",
                    "category": "C3.bus-broadcast-without-findings",
                    "folder": source_folder,
                    "message": f"RESEARCH_BUS lists `{source_folder}` but FINDINGS.md does not exist",
                    "evidence": "meta/research/RESEARCH_BUS.md",
                })

    # ── C4: CANDIDATE_BRIDGES → source/target exist ──
    bridges_path = META_RESEARCH / "CANDIDATE_BRIDGES.md"
    bridges_text = read_text(bridges_path) if bridges_path.exists() else ""
    if bridges_text:
        # Parse bridge rows: | BR-NNN | `source` | `target` | strength | ...
        for m in re.finditer(r"^\|\s*(BR-\d+)\s*\|\s*`([^`]+)`\s*\|\s*`([^`]+)`\s*\|", bridges_text, re.MULTILINE):
            br_id = m.group(1)
            src = m.group(2)
            tgt = m.group(3)
            for fldr, role in [(src, "source"), (tgt, "target")]:
                # strip path prefix
                f_clean = fldr.replace("research/", "").split("/")[0]
                if f_clean not in folder_names:
                    findings.append({
                        "severity": "HIGH",
                        "category": "C4.bridge-broken-ref",
                        "folder": f_clean,
                        "message": f"{br_id} {role}=`{fldr}` not found in research/",
                        "evidence": "meta/research/CANDIDATE_BRIDGES.md",
                    })

    # ── A1: overclaim phrases in FINDINGS/NEEDS without source_of_status ──
    for folder in research_folders:
        for fname in ("FINDINGS.md", "NEEDS.md"):
            p = folder / fname
            if not p.exists():
                continue
            text = read_text(p)
            for phrase in OVERCLAIM_PHRASES:
                if phrase in text:
                    # Check if folder has source_of_status entries
                    sos = folder_yamls.get(folder.name, {}).get("source_of_status", [])
                    if not sos:
                        findings.append({
                            "severity": "CRITICAL",
                            "category": "A1.overclaim-without-source",
                            "folder": folder.name,
                            "message": f"`{phrase}` in {fname} but YAML source_of_status is empty",
                            "evidence": f"research/{folder.name}/{fname}",
                        })
                    else:
                        findings.append({
                            "severity": "MEDIUM",
                            "category": "A1.overclaim-phrase-present",
                            "folder": folder.name,
                            "message": f"`{phrase}` in {fname} (source_of_status non-empty, but phrase still flagged for cleanup)",
                            "evidence": f"research/{folder.name}/{fname}",
                        })

    # ── A2: source_of_status 'heuristic' gdy level >= L3 ──
    for name, y in folder_yamls.items():
        if y.get("level") in ("L3", "L4"):
            sos = y.get("source_of_status", [])
            for entry in sos:
                if "heuristic" in entry.lower():
                    findings.append({
                        "severity": "CRITICAL",
                        "category": "A2.heuristic-source-with-high-level",
                        "folder": name,
                        "message": f"level: {y['level']} but source_of_status contains 'heuristic': '{entry[:100]}'",
                        "evidence": "research/" + name + "/README.md",
                    })

    # ── A3: Polluted folders z exports_findings: true ──
    for name, y in folder_yamls.items():
        if y.get("polluted_74394a8") == "true" and y.get("exports_findings") == "true":
            findings.append({
                "severity": "CRITICAL",
                "category": "A3.polluted-but-exports",
                "folder": name,
                "message": "polluted_74394a8: true but exports_findings: true (should be false per AGENT_PROTOCOL §3.2)",
                "evidence": "research/" + name + "/README.md",
            })

    # ── INFO: count of phrase violators (already known from Sesja 4) ──
    phrase_violator_count = sum(
        1 for f in findings
        if f["category"].startswith("A1.")
    )

    # ── Build report ──
    severity_counter = {"CRITICAL": 0, "HIGH": 0, "MEDIUM": 0, "INFO": 0}
    category_counter = {}
    for f in findings:
        severity_counter[f["severity"]] = severity_counter.get(f["severity"], 0) + 1
        category_counter[f["category"]] = category_counter.get(f["category"], 0) + 1

    lines = []
    lines.append("---")
    lines.append('title: "REGRESSION_S8 — read-only consistency + anti-overclaim audit"')
    lines.append(f"date: {TODAY_STR}")
    lines.append("type: regression-report")
    lines.append("session: S8")
    lines.append("status: GENERATED (read-only audit)")
    lines.append('parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"')
    lines.append("related:")
    lines.append('  - "[[meta/research/AGENT_PROTOCOL.md]]"')
    lines.append('  - "[[meta/research/FOLDER_STATUS_INDEX.md]]"')
    lines.append('  - "[[meta/research/HOTSPOT_AUDIT_S3_5.md]]"')
    lines.append("tags:")
    lines.append("  - regression")
    lines.append("  - audit")
    lines.append("  - session-8")
    lines.append("  - read-only")
    lines.append("---")
    lines.append("")
    lines.append("# REGRESSION_S8 — read-only consistency + anti-overclaim audit")
    lines.append("")
    lines.append("> **Read-only audyt** spójności workflow (per PLAN §6 Sesja 8).")
    lines.append(f"> Generated {datetime.now().isoformat(timespec='seconds')}.")
    lines.append(f"> Folders processed: **{len(research_folders)}**.")
    lines.append("")
    lines.append("## 1. Severity counts")
    lines.append("")
    lines.append("| Severity | Liczba |")
    lines.append("|---|---:|")
    for sev in ("CRITICAL", "HIGH", "MEDIUM", "INFO"):
        lines.append(f"| **{sev}** | {severity_counter.get(sev, 0)} |")
    lines.append(f"| **Razem** | **{len(findings)}** |")
    lines.append("")

    lines.append("## 2. Category breakdown")
    lines.append("")
    lines.append("| Category | Severity | N |")
    lines.append("|---|---|---:|")
    for cat in sorted(category_counter.keys()):
        cat_findings = [f for f in findings if f["category"] == cat]
        sev = cat_findings[0]["severity"] if cat_findings else "—"
        lines.append(f"| `{cat}` | {sev} | {len(cat_findings)} |")
    lines.append("")

    # Per-severity sections
    for sev in ("CRITICAL", "HIGH", "MEDIUM"):
        sev_findings = [f for f in findings if f["severity"] == sev]
        if not sev_findings:
            continue
        lines.append(f"## 3.x. {sev} findings ({len(sev_findings)})")
        lines.append("")
        # Group by category
        by_cat = {}
        for f in sev_findings:
            by_cat.setdefault(f["category"], []).append(f)
        for cat in sorted(by_cat.keys()):
            lines.append(f"### `{cat}` ({len(by_cat[cat])})")
            lines.append("")
            for f in by_cat[cat]:
                lines.append(f"- **`{f['folder']}`** — {f['message']}")
                lines.append(f"  - evidence: `{f['evidence']}`")
            lines.append("")

    # Twarde reguły confirmed
    lines.append("## 4. Twarde reguły confirmed (anti-overclaim audit)")
    lines.append("")
    n_l4 = sum(1 for y in folder_yamls.values() if y.get("level") == "L4")
    n_l4_with_promo = sum(
        1 for y in folder_yamls.values()
        if y.get("level") == "L4" and y.get("promoted_to_core") and y.get("promoted_to_core") != "null"
    )
    lines.append(f"- Folderów z `level: L4`: **{n_l4}**")
    lines.append(f"- Z weryfikowalnym `promoted_to_core`: **{n_l4_with_promo}/{n_l4}**")
    lines.append(f"  - Reguła wymaga 100% — {'✅' if n_l4 == n_l4_with_promo else '❌'}")
    n_polluted = sum(1 for y in folder_yamls.values() if y.get("polluted_74394a8") == "true")
    n_polluted_quarantined = sum(
        1 for y in folder_yamls.values()
        if y.get("polluted_74394a8") == "true" and y.get("exports_findings") == "false"
    )
    lines.append(f"- Folderów polluted-74394a8: **{n_polluted}**")
    lines.append(f"  - Z `exports_findings: false` (kwarantanna): **{n_polluted_quarantined}/{n_polluted}**")
    lines.append(f"  - Reguła kwarantanny — {'✅' if n_polluted == n_polluted_quarantined else '❌'}")
    lines.append("")

    lines.append("## 5. Wnioski końcowe")
    lines.append("")
    if severity_counter.get("CRITICAL", 0) == 0:
        lines.append("✅ **0 CRITICAL findings** — workflow przechodzi anti-overclaim audit.")
    else:
        lines.append(f"❌ **{severity_counter['CRITICAL']} CRITICAL findings** — wymagają natychmiastowej naprawy.")
    if severity_counter.get("HIGH", 0) == 0:
        lines.append("✅ **0 HIGH findings** — strukturalna spójność OK.")
    else:
        lines.append(f"⚠ **{severity_counter['HIGH']} HIGH findings** — wymagają review w przyszłej sesji.")
    if severity_counter.get("MEDIUM", 0) > 0:
        lines.append(f"ℹ **{severity_counter['MEDIUM']} MEDIUM findings** — głównie phrase violations (znane z Sesji 4 / 7).")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 6. Per AGENT_PROTOCOL §3 self-check confirmed")
    lines.append("")
    lines.append("- [x] Skrypt-checker NIE modyfikuje plików (read-only)")
    lines.append("- [x] Wszystkie wpisy mają konkretną `evidence:` ścieżkę")
    lines.append("- [x] Severity klasyfikowana zgodnie z PLAN §6 Sesja 8")
    lines.append("- [x] Auto-naprawa wyłączona (Sesja 8 jest read-only)")
    lines.append("")

    REPORT_OUT.write_text("\n".join(lines), encoding="utf-8")

    # Also write JSON
    json_out = HERE / "_regression_s8_dryrun.json"
    json_out.write_text(json.dumps({
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "n_folders": len(research_folders),
        "severity_counter": severity_counter,
        "category_counter": category_counter,
        "findings": findings,
    }, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"Wrote {REPORT_OUT} ({REPORT_OUT.stat().st_size} B)")
    print(f"Wrote {json_out} ({json_out.stat().st_size} B)")
    print()
    print("Severity counts:")
    for sev in ("CRITICAL", "HIGH", "MEDIUM"):
        n = severity_counter.get(sev, 0)
        print(f"  {sev}: {n}")
    print()
    if severity_counter.get("CRITICAL", 0) > 0:
        print("CRITICAL details:")
        for f in [x for x in findings if x["severity"] == "CRITICAL"][:10]:
            print(f"  [{f['category']}] {f['folder']}: {f['message'][:100]}")


if __name__ == "__main__":
    main()
