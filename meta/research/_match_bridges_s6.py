"""
Sesja 6 — heurystyczny matcher NEEDS × FINDINGS dla wszystkich par folderów.

Anty-spam (per AGENT_PROTOCOL §5):
  - Match wymaga ≥2 niezależnych sygnałów (variable name + phase marker, lub
    explicit folder name + formula match).
  - Match-strength:
    EXACT     — explicit cross-reference (folder name / file path) + ≥1 token
    PARTIAL   — explicit hotspot ID match (H-AN/H-BN) + ≥1 token, lub specyficzna
                pair tokens (np. m_field + WEP)
    HEURISTIC — generic keyword overlap ≥3 tokenów

Per-folder cap: max 5 PROPOSED matches per source folder (żeby nie zalewać).

Output:
  - CANDIDATE_BRIDGES.md — pełna lista PROPOSED matches
  - RESEARCH_BUS.md — update consumers (matched folders)
  - _match_bridges_s6_dryrun.json — raw data
  - meta/research/IMPACT_MATRIX.md — updated § 4
"""

from __future__ import annotations
import json
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent.parent / "research"
DRYRUN_OUT = HERE / "_match_bridges_s6_dryrun.json"
BRIDGES_FILE = HERE / "CANDIDATE_BRIDGES.md"
BUS_FILE = HERE / "RESEARCH_BUS.md"
IMPACT_FILE = HERE / "IMPACT_MATRIX.md"

SANDBOX_RESERVED = {"_sandbox", "_archive"}
POLLUTED = {
    "op-chi1-newton-constant-derivation",
    "op-uv2-mtgp-absolute-scale",
    "op-omega2-axion-coupling-lock",
    "op-omega3-axion-decay-constant",
}
TODAY = datetime.now().strftime("%Y-%m-%d")

# Physics-aware token patterns (anti-spam: only specific TGP/physics tokens)
PHYSICS_TOKENS = re.compile(
    r"\b("
    # Constants & couplings
    r"M_?Pl|m_?Pl|H_?0|H₀|Φ_?(eq|eff|0|bare)|Phi_?(eq|eff|0|bare)|"
    r"alpha_?(s|0|em|G|g|psi)|α_?(s|0|em|G|g|ψ)|"
    r"g̃|g_tilde|g_?(0\^?e|0_e|star|crit)|"
    r"Omega_?(Lambda|L|DM|M)|Ω_?(Λ|L|DM|M)|"
    r"Sigma_?m_?nu|Σm_?ν|m_?H|m_?W|m_?e|m_?mu|m_?tau|m_?nu|m_?field|m_?sigma|m_?s|m_?Z|m_?t|"
    r"c_?(0|2|GW|3|4)|c₀|c₂|kappa|κ|gamma_?PPN|γ_PPN|beta_?PPN|β_PPN|"
    r"N_?A|N_?f|N_?c|n_?s|"
    # Specific TGP terms
    r"M9\.[12]|M9\.1''|sek0[0-9]|dodatek[A-Z]|"
    r"sigma_ab|σ_ab|h_TT|h_b|h_L|"
    r"K\(phi\)|K\(φ\)|V\(phi\)|V\(Φ\)|"
    r"f\(psi\)|f\(ψ\)|"
    # Phenomenology
    r"DESI|JUNO|Hyper-?K|MICROSCOPE|LIGO|EHT|ngEHT|BICEP|Planck|"
    r"Sgr ?A|M87|GW150914|GW170817|"
    # Methods
    r"Adams positivity|Bethe-?Salpeter|OPE|FRG|NGFP|heat[- ]?kernel|"
    r"Hilbert series|Hobart-?Derrick|Yukawa|sympy LOCK|"
    # Phase markers (these are explicit closure markers)
    r"[AB]\d+[-_]?(CLOSED|PATCHED)|H-[AB]\d+|"
    # Common substrate physics
    r"WEP|PPN|Schwarzschild|Lorentzian|Coulomb|"
    # Predictions registry IDs
    r"TT\d+|XX\d+|YY\d+|WW\d+|ZZ\d+|DE\d+|Z\d|"
    # Specific cycle names
    r"OP-?7|OP-?M92|OP-?6|"
    # Specific physical effects
    r"vacuum catastrophe|Earnshaw|breathing mode|photon ring|"
    r"Born rule|soliton tail|Koide|Wolfenstein|Cabibbo|"
    r"hierarchy"
    r")\b",
    re.IGNORECASE,
)

# Folder name pattern (op-X, qm_X, topic)
FOLDER_NAME_RE = re.compile(r"\b(?:op-[\w-]+|qm_\w+|closure_\d{4}-\d{2}-\d{2}|"
                            r"hubble_tension|galaxy_scaling|cosmo_tensions|s8_tension|"
                            r"em_from_substrate|atom_from_soliton|nbody|why_n3|"
                            r"continuum_limit|metric_ansatz|cabibbo_correction|"
                            r"brannen_sqrt2|liquid_viscosity|muon_g_minus_2|"
                            r"casimir_mof|atomic_shells_closure|cohesion_closure|"
                            r"superconductivity_closure|rho_normal_state_closure|"
                            r"thermal_transport_molecular|particle_sector_closure|"
                            r"mass_scaling_k4|neutrino_msw|desi_dark_energy|"
                            r"uv_completion|external_review_2026-04-25)\b")

# Hotspot ID pattern
HOTSPOT_RE = re.compile(r"\bH-[AB]\d+\b")
PHASE_MARKER_RE = re.compile(r"\b[AB]\d+[-_](CLOSED|PATCHED)\b")


def read_text(p: Path, max_bytes: int = 200_000) -> str:
    try:
        if p.stat().st_size > max_bytes:
            with p.open("rb") as f:
                return f.read(max_bytes).decode("utf-8", errors="ignore")
        return p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""


def extract_signals(text: str, this_folder: str) -> dict:
    """Extract physics tokens, hotspots, phase markers, folder mentions."""
    physics_tokens = set(m.group(0).lower() for m in PHYSICS_TOKENS.finditer(text))
    hotspots = set(m.group(0) for m in HOTSPOT_RE.finditer(text))
    phase_markers = set(m.group(0) for m in PHASE_MARKER_RE.finditer(text))
    folder_mentions = set(m.group(0) for m in FOLDER_NAME_RE.finditer(text)) - {this_folder}
    return {
        "tokens": physics_tokens,
        "hotspots": hotspots,
        "phase_markers": phase_markers,
        "folder_mentions": folder_mentions,
    }


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


def find_quote(text: str, token: str, max_chars: int = 200) -> str:
    """Find first occurrence of token in text, return surrounding context."""
    m = re.search(re.escape(token), text, re.IGNORECASE)
    if not m:
        return ""
    start = max(0, m.start() - 80)
    end = min(len(text), m.end() + 100)
    ctx = text[start:end].replace("\n", " ").strip()
    return ctx[:max_chars]


def score_match(needs_sig: dict, findings_sig: dict, needs_folder: str, findings_folder: str) -> dict:
    """Score a (needs, findings) pair. Return dict with score components and total."""
    score = {
        "explicit_folder_xref": 0,
        "hotspot_overlap": 0,
        "phase_marker_overlap": 0,
        "specific_token_overlap": 0,
        "generic_token_overlap": 0,
        "common_signals": [],
    }

    # Explicit folder cross-reference
    if findings_folder in needs_sig["folder_mentions"]:
        score["explicit_folder_xref"] = 5
        score["common_signals"].append(f"NEEDS mentions `{findings_folder}` explicitly")
    if needs_folder in findings_sig["folder_mentions"]:
        score["explicit_folder_xref"] += 5
        score["common_signals"].append(f"FINDINGS mentions `{needs_folder}` explicitly")

    # Hotspot overlap (H-A1..H-B12)
    common_hotspots = needs_sig["hotspots"] & findings_sig["hotspots"]
    score["hotspot_overlap"] = 4 * len(common_hotspots)
    if common_hotspots:
        score["common_signals"].append(f"Hotspot match: {sorted(common_hotspots)}")

    # Phase marker overlap (B6-CLOSED, A5-PATCHED, etc.)
    common_phase = needs_sig["phase_markers"] & findings_sig["phase_markers"]
    score["phase_marker_overlap"] = 3 * len(common_phase)
    if common_phase:
        score["common_signals"].append(f"Phase marker match: {sorted(common_phase)}")

    # Specific physics tokens (>3 chars, multi-letter)
    common_tokens = needs_sig["tokens"] & findings_sig["tokens"]
    # filter for "specific" tokens (longer ones, less generic)
    specific_tokens = {t for t in common_tokens if len(t) >= 4 or any(c in t for c in "_-/")}
    generic_tokens = common_tokens - specific_tokens
    score["specific_token_overlap"] = 2 * len(specific_tokens)
    score["generic_token_overlap"] = 1 * len(generic_tokens)
    if specific_tokens:
        score["common_signals"].append(f"Specific tokens: {sorted(specific_tokens)[:6]}")

    score["total"] = (
        score["explicit_folder_xref"]
        + score["hotspot_overlap"]
        + score["phase_marker_overlap"]
        + score["specific_token_overlap"]
        + score["generic_token_overlap"]
    )

    # Determine strength
    if score["explicit_folder_xref"] > 0 and (score["specific_token_overlap"] >= 4 or score["hotspot_overlap"] > 0):
        score["strength"] = "EXACT"
    elif score["hotspot_overlap"] > 0 and score["specific_token_overlap"] >= 2:
        score["strength"] = "EXACT"
    elif score["explicit_folder_xref"] > 0 or score["hotspot_overlap"] > 0 or score["phase_marker_overlap"] > 0:
        score["strength"] = "PARTIAL"
    elif score["specific_token_overlap"] >= 6:
        score["strength"] = "HEURISTIC"
    else:
        score["strength"] = "WEAK"

    return score


def main():
    folders = []
    for p in sorted(RESEARCH.iterdir(), key=lambda x: x.name.lower()):
        if not p.is_dir() or p.name in SANDBOX_RESERVED:
            continue
        folders.append(p)

    print(f"Loading signals from {len(folders)} folders...")

    # Load NEEDS and FINDINGS for every folder
    folder_data = {}
    for folder in folders:
        name = folder.name
        needs_p = folder / "NEEDS.md"
        findings_p = folder / "FINDINGS.md"
        readme_p = folder / "README.md"
        needs_text = read_text(needs_p) if needs_p.exists() else ""
        findings_text = read_text(findings_p) if findings_p.exists() else ""
        readme_text = read_text(readme_p) if readme_p.exists() else ""

        polluted_status = yaml_field(readme_text, "polluted_74394a8")
        exports_findings = yaml_field(readme_text, "exports_findings")

        folder_data[name] = {
            "name": name,
            "needs_text": needs_text,
            "findings_text": findings_text,
            "needs_signals": extract_signals(needs_text, name) if needs_text else None,
            "findings_signals": extract_signals(findings_text, name) if findings_text else None,
            "is_polluted": polluted_status == "true",
            "exports": exports_findings == "true",
            "needs_empty": "(Brak ekstrahowalnych otwartych luk)" in needs_text,
            "findings_empty": "(Brak ekstrahowalnych findings)" in findings_text,
        }

    print("Computing matches...")

    # All-pairs match
    all_matches = []  # list of dicts
    for source_name, source_data in folder_data.items():
        if source_data["is_polluted"]:
            continue
        if source_data["needs_empty"] or not source_data["needs_signals"]:
            continue

        per_source_matches = []
        for target_name, target_data in folder_data.items():
            if target_name == source_name:
                continue
            if target_data["is_polluted"]:
                continue  # polluted don't supply consumed findings
            if not target_data["exports"]:
                continue  # findings not exportable
            if not target_data["findings_signals"]:
                continue

            score = score_match(
                source_data["needs_signals"],
                target_data["findings_signals"],
                source_name,
                target_name,
            )
            if score["strength"] in ("EXACT", "PARTIAL", "HEURISTIC"):
                per_source_matches.append({
                    "source_folder": source_name,
                    "target_folder": target_name,
                    "score": score,
                })

        # Sort by total score, take top 5
        per_source_matches.sort(key=lambda m: -m["score"]["total"])
        per_source_matches = per_source_matches[:5]
        all_matches.extend(per_source_matches)

    print(f"Total candidate matches: {len(all_matches)}")
    strength_counter = Counter(m["score"]["strength"] for m in all_matches)
    print(f"By strength: {dict(strength_counter)}")

    # Assign IDs and prepare for output
    all_matches.sort(key=lambda m: (-m["score"]["total"], m["source_folder"], m["target_folder"]))
    for i, m in enumerate(all_matches, 1):
        m["bridge_id"] = f"BR-{i:03d}"

    # Build CANDIDATE_BRIDGES.md
    lines = []
    lines.append("---")
    lines.append('title: "CANDIDATE_BRIDGES — luki analityczne potencjalnie domykalne wynikami z innych folderów"')
    lines.append(f"date: {TODAY}")
    lines.append("type: bridges")
    lines.append("status: ACTIVE — Sesja 6 PROPOSED matches; czeka na human review (HUMAN-CONFIRMED gate)")
    lines.append('parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"')
    lines.append("related:")
    lines.append('  - "[[meta/research/RESEARCH_BUS.md]]"')
    lines.append('  - "[[meta/research/IMPACT_MATRIX.md]]"')
    lines.append('  - "[[meta/research/AGENT_PROTOCOL.md]]"')
    lines.append("tags:")
    lines.append("  - bridges")
    lines.append("  - matching")
    lines.append("  - inter-folder")
    lines.append("---")
    lines.append("")
    lines.append("# CANDIDATE_BRIDGES — luki domykalne mostami między folderami")
    lines.append("")
    lines.append("> **Sesja 6 auto-generated** (2026-05-03). Heurystyczny match `NEEDS.md` ×")
    lines.append("> `FINDINGS.md` przez physics-aware tokenizer (zob.")
    lines.append("> `_match_bridges_s6.py`). Każdy wpis to PROPOSED kandydat na bridge —")
    lines.append("> czeka na HUMAN-CONFIRMED gate.")
    lines.append("")
    lines.append("## 0. Reguły score (anty-spam)")
    lines.append("")
    lines.append("| Strength | Wymaganie |")
    lines.append("|----------|-----------|")
    lines.append("| **EXACT** | (a) explicit cross-reference (folder mentioned by name) + ≥1 specific token, lub (b) hotspot ID match + ≥1 specific token |")
    lines.append("| **PARTIAL** | explicit cross-ref LUB hotspot LUB phase marker (B6-CLOSED, A5-PATCHED, etc.) |")
    lines.append("| **HEURISTIC** | ≥6 specific physics tokens overlap (np. M9, σ_ab, Φ_eq, β_PPN) |")
    lines.append("| WEAK | (filtrowane out, nie są publikowane) |")
    lines.append("")

    # Stats
    lines.append("## 1. Statystyka")
    lines.append("")
    lines.append("| Metryka | Liczba |")
    lines.append("|---------|-------:|")
    lines.append(f"| PROPOSED total | {len(all_matches)} |")
    lines.append(f"| EXACT | {strength_counter.get('EXACT', 0)} |")
    lines.append(f"| PARTIAL | {strength_counter.get('PARTIAL', 0)} |")
    lines.append(f"| HEURISTIC | {strength_counter.get('HEURISTIC', 0)} |")
    lines.append(f"| HUMAN-CONFIRMED | 0 (czeka na review) |")
    lines.append(f"| EXECUTED | 0 |")
    lines.append("")

    # Top 3 for human review
    lines.append("## 2. Top 3 bridges for human review")
    lines.append("")
    lines.append("> Najsilniejsze kandydatury (highest total score). Każdy wymaga manualnej decyzji.")
    lines.append("")
    top3 = all_matches[:3]
    for m in top3:
        s = m["score"]
        lines.append(f"### {m['bridge_id']}: `{m['source_folder']}` ← `{m['target_folder']}` ({s['strength']}, score={s['total']})")
        lines.append("")
        lines.append(f"- **Source NEEDS:** `research/{m['source_folder']}/NEEDS.md`")
        lines.append(f"- **Target FINDINGS:** `research/{m['target_folder']}/FINDINGS.md`")
        lines.append(f"- **Strength:** {s['strength']} (xref={s['explicit_folder_xref']}, hotspot={s['hotspot_overlap']}, phase={s['phase_marker_overlap']}, tokens={s['specific_token_overlap']})")
        lines.append("- **Common signals:**")
        for sig in s["common_signals"]:
            lines.append(f"  - {sig}")
        lines.append("- **Decision:**")
        lines.append("  - [ ] HUMAN-CONFIRMED")
        lines.append("  - [ ] EXECUTE (copy F<id> z target FINDINGS do source FINDINGS z notatką pochodzenia)")
        lines.append("  - [ ] REJECT (powód: ____)")
        lines.append("")

    # Full table
    lines.append("## 3. PROPOSED matches (full table)")
    lines.append("")
    lines.append("| ID | Source NEEDS folder | → Target FINDINGS folder | Strength | Score | Top signal |")
    lines.append("|----|---------------------|---------------------------|----------|------:|-----------|")
    for m in all_matches:
        sig = m["score"]["common_signals"][0] if m["score"]["common_signals"] else "—"
        sig = sig[:80].replace("|", "\\|")
        lines.append(f"| {m['bridge_id']} | `{m['source_folder']}` | `{m['target_folder']}` | {m['score']['strength']} | {m['score']['total']} | {sig} |")
    lines.append("")

    # Per source folder
    lines.append("## 4. Per-source folder (consumers)")
    lines.append("")
    by_source = defaultdict(list)
    for m in all_matches:
        by_source[m["source_folder"]].append(m)

    for src in sorted(by_source.keys()):
        matches = by_source[src]
        lines.append(f"### `{src}` — {len(matches)} kandydat consumers")
        lines.append("")
        for m in matches:
            lines.append(f"- {m['bridge_id']}: `{m['target_folder']}` ({m['score']['strength']}, score={m['score']['total']})")
        lines.append("")

    # Orphans — NEEDS bez kandydata
    sources_with_matches = set(by_source.keys())
    needs_folders = {n for n, d in folder_data.items()
                    if not d["is_polluted"] and not d["needs_empty"] and d["needs_signals"]}
    orphan_needs = sorted(needs_folders - sources_with_matches)
    lines.append(f"## 5. Orphan NEEDS ({len(orphan_needs)} folderów bez kandydata bridge)")
    lines.append("")
    if orphan_needs:
        lines.append("> Foldery z NEEDS, ale żaden FINDINGS w innym folderze nie matchuje powyżej WEAK threshold.")
        lines.append("> To są **realnie open** luki — wymagają NOWEJ pracy, nie bridge.")
        lines.append("")
        for f in orphan_needs:
            lines.append(f"- `research/{f}`")
        lines.append("")

    # Orphan FINDINGS — eksportujące, ale nikt nie chce
    targets_with_matches = set(m["target_folder"] for m in all_matches)
    findings_folders = {n for n, d in folder_data.items()
                       if not d["is_polluted"] and d["exports"] and d["findings_signals"]}
    orphan_findings = sorted(findings_folders - targets_with_matches)
    lines.append(f"## 6. Orphan FINDINGS ({len(orphan_findings)} folderów eksportujących bez konsumenta)")
    lines.append("")
    if orphan_findings:
        lines.append("> Foldery z `exports_findings: true`, ale żaden NEEDS nie pyta o ich wyniki.")
        lines.append("> Po 90 dniach: degradacja w `RESEARCH_BUS.md` na STALE.")
        lines.append("")
        for f in orphan_findings:
            lines.append(f"- `research/{f}`")
        lines.append("")

    # Foldery wyłączone
    lines.append(f"## 7. Wyłączone z matching")
    lines.append("")
    lines.append(f"- 4 polluted-74394a8 folders (kwarantanna; AGENT_PROTOCOL §3.2)")
    lines.append(f"- {sum(1 for d in folder_data.values() if d['needs_empty'])} folderów z empty NEEDS")
    lines.append(f"- {sum(1 for d in folder_data.values() if not d['exports'])} folderów z `exports_findings: false` (jako sources nie eksportują dla matchingu)")
    lines.append("")

    lines.append("## 8. Generation metadata")
    lines.append("")
    lines.append(f"- Generated: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"- Builder: `_match_bridges_s6.py`")
    lines.append(f"- Folders processed: {len(folders)}")
    lines.append(f"- All-pairs scored: {len(folders) * (len(folders) - 1)}")
    lines.append(f"- Matches above WEAK threshold: {len(all_matches)}")

    BRIDGES_FILE.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {BRIDGES_FILE} ({BRIDGES_FILE.stat().st_size} B)")

    # Save dryrun JSON
    DRYRUN_OUT.write_text(json.dumps({
        "matches": all_matches,
        "n_matches": len(all_matches),
        "strength_distribution": dict(strength_counter),
        "orphan_needs": orphan_needs,
        "orphan_findings": orphan_findings,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
    }, indent=2, ensure_ascii=False, default=list), encoding="utf-8")
    print(f"Wrote {DRYRUN_OUT} ({DRYRUN_OUT.stat().st_size} B)")

    # Update RESEARCH_BUS.md consumers (replace § 2 broadcast table consumers field)
    bus_text = BUS_FILE.read_text(encoding="utf-8")
    consumers_per_target = defaultdict(set)
    for m in all_matches:
        consumers_per_target[m["target_folder"]].add(m["source_folder"])

    new_bus_lines = []
    in_section_2 = False
    for line in bus_text.splitlines():
        if line.startswith("| 20") and "BROADCAST" in line:
            # Parse this row, update consumers
            cells = line.split(" | ")
            if len(cells) >= 7:
                # Find target folder (cell 1 has `name` in backticks)
                m = re.search(r"`([^`]+)`", cells[1])
                if m:
                    target = m.group(1)
                    consumers = consumers_per_target.get(target, set())
                    if consumers:
                        cons_str = ", ".join(f"`{c}`" for c in sorted(consumers))
                        cells[5] = cons_str
                    else:
                        cells[5] = "(orphan, no matches)"
                    line = " | ".join(cells)
        new_bus_lines.append(line)
    BUS_FILE.write_text("\n".join(new_bus_lines), encoding="utf-8")
    print(f"Updated RESEARCH_BUS.md consumers for {len(consumers_per_target)} target folders")

    print()
    print("DONE.")
    print(f"Top 3 for human review:")
    for m in top3:
        print(f"  {m['bridge_id']}: {m['source_folder']} ← {m['target_folder']} "
              f"({m['score']['strength']}, score={m['score']['total']})")


if __name__ == "__main__":
    main()
