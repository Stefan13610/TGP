"""
Sesja 5 post-step — update YAML flags has_findings_file / has_needs_file
po wygenerowaniu plików FINDINGS.md i NEEDS.md.

Plus update exports_findings: jeśli FINDINGS.md zawiera "(Brak ekstrahowalnych findings)"
→ exports_findings: false; w przeciwnym razie true.

Plus update last_yaml_update na dziś.
"""

from __future__ import annotations
import re
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent.parent / "research"
SANDBOX_RESERVED = {"_sandbox", "_archive"}
TODAY = datetime.now().strftime("%Y-%m-%d")


def update_readme_flags(folder: Path) -> dict:
    readme = folder / "README.md"
    findings = folder / "FINDINGS.md"
    needs = folder / "NEEDS.md"

    has_findings = findings.exists()
    has_needs = needs.exists()

    # Determine exports_findings: read FINDINGS.md, check for empty marker
    exports_findings = False
    if has_findings:
        ftext = findings.read_text(encoding="utf-8", errors="ignore")
        if "(Brak ekstrahowalnych findings)" not in ftext:
            exports_findings = True

    # Polluted folders never export findings
    if "polluted_74394a8: true" in (readme.read_text(encoding="utf-8", errors="ignore") if readme.exists() else ""):
        exports_findings = False

    if not readme.exists():
        return {"folder": folder.name, "skipped": "no README"}

    text = readme.read_text(encoding="utf-8")
    if not text.lstrip().startswith("---"):
        return {"folder": folder.name, "skipped": "no frontmatter"}

    parts = text.split("---", 2)
    if len(parts) < 3:
        return {"folder": folder.name, "skipped": "malformed frontmatter"}

    fm = parts[1]
    body = parts[2]

    # Replace flag values
    def replace_flag(fm_text, key, new_val):
        # Match "  key: value" within tgp_status block
        pattern = re.compile(rf"(\n  {re.escape(key)}:\s*)([^\n]*)")
        repl = f"\\g<1>{new_val}"
        new_fm, n = pattern.subn(repl, fm_text)
        return new_fm, n

    new_fm = fm
    changes = {}

    new_fm, n = replace_flag(new_fm, "has_findings_file", "true" if has_findings else "false")
    if n: changes["has_findings_file"] = "true" if has_findings else "false"
    new_fm, n = replace_flag(new_fm, "has_needs_file", "true" if has_needs else "false")
    if n: changes["has_needs_file"] = "true" if has_needs else "false"
    new_fm, n = replace_flag(new_fm, "exports_findings", "true" if exports_findings else "false")
    if n: changes["exports_findings"] = "true" if exports_findings else "false"
    new_fm, n = replace_flag(new_fm, "last_yaml_update", f'"{TODAY}"')
    if n: changes["last_yaml_update"] = TODAY

    new_text = "---" + new_fm + "---" + body
    readme.write_text(new_text, encoding="utf-8")
    return {"folder": folder.name, "changes": changes}


def main():
    folders = []
    for p in sorted(RESEARCH.iterdir(), key=lambda x: x.name.lower()):
        if not p.is_dir() or p.name in SANDBOX_RESERVED:
            continue
        folders.append(p)

    print(f"Updating YAML flags for {len(folders)} folders...")
    n_updated = 0
    n_skipped = 0
    exports_count = 0
    for folder in folders:
        result = update_readme_flags(folder)
        if "skipped" in result:
            n_skipped += 1
            continue
        n_updated += 1
        if result.get("changes", {}).get("exports_findings") == "true":
            exports_count += 1

    print(f"Updated: {n_updated}; Skipped: {n_skipped}")
    print(f"Folders with exports_findings: true: {exports_count}")


if __name__ == "__main__":
    main()
