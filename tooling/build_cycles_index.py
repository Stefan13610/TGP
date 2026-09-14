# -*- coding: utf-8 -*-
"""build_cycles_index.py — generuje CYCLES.tsv: jedna linia na cykl w research/.

Po co: INDEX.md ma 336 KB — agent nie moze go przeczytac, zeby zobaczyc mape cykli.
CYCLES.tsv to ta sama mapa w ~40 KB: greppowalna, sortowalna, tania w tokenach.

Zrodlo prawdy = frontmatter README.md kazdego folderu (pole `folder_status`
znormalizowane przez tooling/normalize_folder_status.py).

UZYCIE:
    python tooling/build_cycles_index.py          # zapisuje CYCLES.tsv w roocie
    python tooling/build_cycles_index.py --check  # tylko sprawdza aktualnosc (exit 1 = nieaktualny)
"""
import io, os, re, sys, datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "research")
OUT = os.path.join(ROOT, "CYCLES.tsv")

COLS = ["slug", "status", "data", "faza", "werdykt", "tytul"]
# "active" = wartosc sprzed normalizacji (folder wylaczony z normalize_folder_status.py) — sortuj na gorze
ORDER = {"wip": 0, "active": 0, "locked": 1, "paused": 2, "parking": 3,
         "closed": 4, "legacy": 5, "?": 6}


def read_fm(path):
    with io.open(path, encoding="utf-8") as f:
        txt = f.read()
    if not txt.startswith("---"):
        return {}
    m = re.search(r"\n---[ \t]*(\r?\n|$)", txt[3:])
    if not m:
        return {}
    fm = {}
    for line in txt[3:3 + m.start()].split("\n"):
        mm = re.match(r"^([a-zA-Z_][\w-]*):\s*(.*)$", line)
        if mm:
            v = re.sub(r"\s+#.*$", "", mm.group(2)).strip().strip('"').strip("'")
            fm[mm.group(1)] = v
    return fm


def clean(s, n=140):
    s = re.sub(r"\s+", " ", (s or "")).replace("\t", " ").strip()
    return (s[:n - 1] + "…") if len(s) > n else s


def phase_of(slug, fm):
    if fm.get("phase"):
        return clean(fm["phase"], 40)
    try:
        files = os.listdir(os.path.join(RES, slug))
    except OSError:
        return ""
    if any(f.startswith("Phase_FINAL") for f in files):
        return "FINAL"
    nums = [int(m.group(1)) for f in files for m in [re.match(r"^Phase(\d)", f)] if m]
    if nums:
        return "Phase %d" % max(nums)
    if any(f.lower().startswith("phase0") for f in files):
        return "Phase 0 (LOCK)"
    return ""


def build():
    rows = []
    for slug in sorted(d for d in os.listdir(RES) if os.path.isdir(os.path.join(RES, d))):
        p = os.path.join(RES, slug, "README.md")
        fm = read_fm(p) if os.path.exists(p) else {}
        status = fm.get("folder_status") or "?"
        # UWAGA: legacy_status to STATUS, nie werdykt — nie uzywac jako fallback.
        verdict = (fm.get("verdict") or fm.get("claim_status")
                   or fm.get("classification") or "")
        rows.append([
            slug,
            status,
            fm.get("date", ""),
            phase_of(slug, fm),
            clean(verdict, 100),
            clean(fm.get("title", ""), 110),
        ])
    rows.sort(key=lambda r: (ORDER.get(r[1], 9), r[2] or "", r[0]), reverse=False)
    rows.sort(key=lambda r: (ORDER.get(r[1], 9), ), reverse=False)

    head = [
        "# CYCLES.tsv — mapa cykli research/ (generowane: tooling/build_cycles_index.py)",
        "# wygenerowano: %s | cykli: %d" % (datetime.date.today().isoformat(), len(rows)),
        "# status: wip=w ruchu | locked=LOCK bez realizacji | paused=zamrozony | parking=pomysl | closed=zamkniety | legacy=generacja sprzed 2026-05",
        "# realny WIP jest w STATE.md §Active WIP — to jest indeks, nie zrodlo prawdy o priorytetach",
        "\t".join(COLS),
    ]
    return "\n".join(head + ["\t".join(r) for r in rows]) + "\n"


def main():
    txt = build()
    if "--check" in sys.argv:
        old = io.open(OUT, encoding="utf-8").read() if os.path.exists(OUT) else ""
        old_b = "\n".join(l for l in old.split("\n") if not l.startswith("# wygenerowano"))
        new_b = "\n".join(l for l in txt.split("\n") if not l.startswith("# wygenerowano"))
        if old_b != new_b:
            print("CYCLES.tsv NIEAKTUALNY — uruchom: python tooling/build_cycles_index.py")
            sys.exit(1)
        print("CYCLES.tsv aktualny")
        return
    with io.open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write(txt)
    print("zapisano %s (%d linii, %.1f KB)" % (OUT, txt.count("\n"), len(txt.encode("utf-8")) / 1024.0))


if __name__ == "__main__":
    main()
