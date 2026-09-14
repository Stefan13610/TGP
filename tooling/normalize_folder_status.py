# -*- coding: utf-8 -*-
"""normalize_folder_status.py — ujednolicenie `folder_status` we frontmatterach README cykli.

PROBLEM (stan 2026-09-14): 263 foldery w research/, z czego 117 ma `folder_status`
w 15 roznych wartosciach, 63 maja tylko wolnotekstowy `status:`, a 80 nie ma nic.
`folder_status: active` mial 27 wystapien przy realnym WIP = 1 (patrz STATE.md).
Agent nie jest w stanie odroznic cyklu zywego od porzuconego.

ROZWIAZANIE: `folder_status` staje sie polem MASZYNOWYM o 6 wartosciach:

    wip      — realnie w ruchu w tej/nastepnej sesji (zrodlo prawdy: STATE.md §Active WIP)
    locked   — Phase0_balance.md zapisany (LOCK), realizacja jeszcze nie ruszyla
    paused   — zaczety, swiadomie zamrozony
    parking  — pomysl zarejestrowany, Phase 0 nie istnieje
    closed   — cykl zamkniety (dowolny werdykt)
    legacy   — folder generacji sprzed 2026-05, status nigdy nie zapisany; referencja, nie WIP

Tresc opisowa NIE GINIE:
  - oryginalna wartosc -> `legacy_status:`
  - komentarz inline (`folder_status: closed-resolved  # CLOSED 2026-06-27: ...`) -> `verdict:`
    (tylko jesli `verdict:` jeszcze nie istnieje)

UZYCIE:
    python tooling/normalize_folder_status.py            # dry-run, raport
    python tooling/normalize_folder_status.py --apply    # zapis
"""
import io, os, re, sys, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "research")

# Cykle realnie w ruchu — synchronizowac recznie ze STATE.md §Active WIP.
WIP = {
    "op-r3-stationary-states-2026-09-14",
}
# Foldery wylaczone z operacji (np. pracuje na nich inny agent).
SKIP = {
    "op-r3-stationary-states-2026-09-14",
}

CLOSED_MARKERS = ("closed", "superseded", "null", "falsified", "negative closure",
                  "positive closure", "no-go closure", "amended-conditional",
                  "analytical decision documented", "framework established",
                  "partial positive", "stage_1_null", "phase 4 closed")
OPEN_MARKERS = ("active", "open", "in_progress", "in progress")


def split_fm(txt):
    """Zwraca (przed, linie_frontmattera, po) albo None gdy brak frontmattera."""
    if not txt.startswith("---"):
        return None
    m = re.search(r"\n---[ \t]*(\r?\n|$)", txt[3:])
    if not m:
        return None
    start = txt.index("\n", 0) + 1          # pierwsza linia po "---"
    end = 3 + m.start() + 1                 # poczatek linii zamykajacej "---"
    return txt[:start], txt[start:end].split("\n"), txt[end:]


def get(fm_lines, key):
    for i, line in enumerate(fm_lines):
        m = re.match(r"^%s:\s*(.*)$" % re.escape(key), line)
        if m:
            raw = m.group(1)
            comment = ""
            cm = re.search(r"\s+#\s*(.*)$", raw)
            if cm:
                comment = cm.group(1).strip()
                raw = raw[:cm.start()]
            return i, raw.strip().strip('"').strip("'"), comment
    return None, None, None


def is_locked_only(slug):
    """Phase0_balance.md istnieje, ale zadna faza realizacyjna nie ruszyla."""
    d = os.path.join(RES, slug)
    try:
        files = os.listdir(d)
    except OSError:
        return False
    has_lock = any(f.lower().startswith("phase0") for f in files)
    has_work = any(re.match(r"(?i)^phase([1-9]|_final|_method|_correction)", f) for f in files)
    return has_lock and not has_work


def classify(slug, fs, status):
    """Zwraca (nowa_wartosc, powod)."""
    if slug in WIP:
        return "wip", "STATE.md Active WIP"
    if fs:
        low = fs.lower()
        if low.startswith("closed") or low in ("amended-conditional",):
            return "closed", "folder_status '%s'" % fs
        if low == "parking":
            return "parking", "folder_status 'parking'"
        if low in ("deferred", "exploratory-note", "open-mixed-verdict", "paused",
                   "needs-bridge"):
            return "paused", "folder_status '%s'" % fs
        if low == "active":
            return "paused", "folder_status 'active', ale brak w STATE.md Active WIP"
        if low == "archive":
            return "legacy", "folder_status 'archive'"
        return "paused", "folder_status '%s' (nierozpoznane)" % fs
    if status:
        low = status.lower()
        if any(k in low for k in CLOSED_MARKERS):
            return "closed", "status: '%s'" % status[:60]
        if any(low.startswith(k) or ("— " + k) in low or ("- " + k) in low
               for k in OPEN_MARKERS) or any(k in low for k in OPEN_MARKERS):
            return "paused", "status: '%s' (otwarty, brak w STATE.md WIP)" % status[:60]
        return "paused", "status: '%s' (nierozpoznane)" % status[:60]
    return "legacy", "brak jakiegokolwiek pola statusu"


def main():
    apply = "--apply" in sys.argv
    dirs = sorted(d for d in os.listdir(RES) if os.path.isdir(os.path.join(RES, d)))
    stats = collections.Counter()
    changed = skipped = 0
    rows = []

    for slug in dirs:
        path = os.path.join(RES, slug, "README.md")
        if not os.path.exists(path):
            stats["BRAK README"] += 1
            rows.append((slug, "-", "-", "brak README.md"))
            continue
        if slug in SKIP:
            skipped += 1
            rows.append((slug, "-", "SKIP", "wylaczony z operacji (SKIP)"))
            continue

        with io.open(path, encoding="utf-8", newline="") as f:
            txt = f.read()
        parts = split_fm(txt)
        if parts is None:
            stats["BRAK FRONTMATTERA"] += 1
            rows.append((slug, "-", "-", "brak frontmattera YAML"))
            continue
        head, fm, tail = parts

        i_fs, fs, fs_comment = get(fm, "folder_status")
        _, status, _ = get(fm, "status")
        new, reason = classify(slug, fs, status)
        if new == "paused" and is_locked_only(slug):
            new, reason = "locked", "Phase0_balance.md bez fazy realizacyjnej (%s)" % reason
        stats[new] += 1
        rows.append((slug, fs or "(brak)", new, reason))

        if fs == new and not fs_comment:
            continue

        line_fs = "folder_status: %s" % new
        extra = []
        if fs and fs != new:
            extra.append('legacy_status: "%s"' % fs)
        if fs_comment:
            i_v, v, _ = get(fm, "verdict")
            if i_v is None:
                extra.append('verdict: "%s"' % fs_comment.replace('"', "'"))
            else:
                extra.append(None)  # verdict juz jest — komentarz zachowany w legacy_comment
                extra[-1] = 'legacy_comment: "%s"' % fs_comment.replace('"', "'")

        if i_fs is not None:
            fm[i_fs:i_fs + 1] = [line_fs] + extra
        else:
            # wstaw po `title:` jesli jest, inaczej na poczatku
            i_t, _, _ = get(fm, "title")
            pos = (i_t + 1) if i_t is not None else 0
            fm[pos:pos] = [line_fs] + extra

        changed += 1
        if apply:
            with io.open(path, "w", encoding="utf-8", newline="") as f:
                f.write(head + "\n".join(fm) + tail)

    print("=== %s ===" % ("ZAPISANO" if apply else "DRY-RUN (dodaj --apply zeby zapisac)"))
    print("folderow: %d | do zmiany: %d | pominietych (SKIP): %d\n" % (len(dirs), changed, skipped))
    print("--- rozklad docelowy ---")
    for k, c in stats.most_common():
        print("  %4d  %s" % (c, k))
    print("\n--- pelna mapa (slug | bylo -> bedzie | powod) ---")
    for slug, old, new, reason in rows:
        if old != new:
            print("  %-52s %-34s -> %-8s %s" % (slug[:52], old[:34], new, reason[:70]))


if __name__ == "__main__":
    main()
