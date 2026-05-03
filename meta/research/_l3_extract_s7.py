"""Extract Phase3 verdicts from 10 L3 folders for Sesja 7 INTAKE generation."""
from __future__ import annotations
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent.parent / "research"
OUT = HERE / "_l3_extract_s7.json"

L3 = [
    "op-alpha-fine-structure",
    "op-cross-sector-charge",
    "op-eps-photon-ring",
    "op-eta-wolfenstein",
    "op-eta2-denom-derivation",
    "op-sc-alpha-origin",
    "op-theta-quark-koide",
    "op-uv-as-ngfp",
    "op-xi-photon-ring",
    "op-zeta-mass-spectrum",
]


def yaml_field(text: str, field: str) -> str | None:
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


def first_paragraph_after_h1(text: str) -> str:
    """Get first paragraph after H1 (skipping frontmatter and headers)."""
    if "---" in text:
        parts = text.split("---", 2)
        if len(parts) >= 3:
            text = parts[2]
    lines = text.splitlines()
    in_h1 = False
    para = []
    for line in lines:
        if line.startswith("# "):
            in_h1 = True
            continue
        if in_h1:
            ls = line.strip()
            if ls.startswith(">") or ls.startswith("##"):
                if para:
                    break
                continue
            if not ls:
                if para:
                    break
                continue
            para.append(ls)
            if len(para) >= 3:
                break
    return " ".join(para)[:400]


def find_boxed_formula(text: str) -> list[str]:
    return re.findall(r"\$\$([^$]{10,200}?)\$\$", text, re.DOTALL)[:3]


def main():
    results = []
    for name in L3:
        folder = RESEARCH / name
        out = {"name": name}

        # Phase 3 results
        p3 = folder / "Phase3_results.md"
        if p3.exists():
            text = p3.read_text(encoding="utf-8", errors="ignore")
            out["cycle"] = yaml_field(text, "cycle")
            out["status"] = yaml_field(text, "status")
            out["verdict"] = yaml_field(text, "verdict")
            out["score"] = yaml_field(text, "overall_score") or yaml_field(text, "score")
            out["program_status"] = yaml_field(text, "program_status")
            out["new_prediction"] = yaml_field(text, "new_falsifiable_prediction")
            out["title"] = yaml_field(text, "title")
            out["phase3_first_para"] = first_paragraph_after_h1(text)
            out["phase3_formulas"] = find_boxed_formula(text)
        else:
            out["error"] = "no Phase3_results.md"

        # README title
        rm = folder / "README.md"
        if rm.exists():
            text = rm.read_text(encoding="utf-8", errors="ignore")
            out["readme_title"] = yaml_field(text, "title")
            for line in text.splitlines():
                if line.startswith("# "):
                    out["readme_h1"] = line[2:].strip()
                    break

        # Program.md cel
        prog = folder / "program.md"
        if prog.exists():
            text = prog.read_text(encoding="utf-8", errors="ignore")
            out["program_title"] = yaml_field(text, "title")
            # Try to find "## Cel" or "## Goal" or "Motywacja" section
            m = re.search(r"^##\s+(?:Cel|Goal|Motywacja|Motivation|TL;?DR|Summary)\b.*?\n+([^\n#]{50,400})",
                         text, re.MULTILINE | re.IGNORECASE)
            if m:
                out["program_cel"] = m.group(1).strip()[:400]

        results.append(out)

    OUT.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote {OUT} ({OUT.stat().st_size} B)")
    for r in results:
        print(f"--- {r['name']} ---")
        print(f"  cycle: {r.get('cycle', '?')}")
        print(f"  verdict: {r.get('verdict', '?')}")
        print(f"  score: {r.get('score', '?')}")
        print(f"  program_status: {r.get('program_status', '?')}")
        print(f"  new_prediction: {r.get('new_prediction', '?')}")
        print()


if __name__ == "__main__":
    main()
