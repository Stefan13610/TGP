# -*- coding: utf-8 -*-
"""similar.py — „czy my to juz liczylismy?" — wyszukiwarka podobienstwa po cyklach research/.

PO CO: grep odpowiada na pytanie „gdzie wystepuje slowo X". Nie odpowiada na pytanie
„czy ktos juz badal to zagadnienie", gdy nie znasz slowa, ktorego uzyl agent pol roku temu.
Przy 263 cyklach to pytanie wraca stale.

CO INDEKSUJE: per cykl skleja README.md + Phase0_balance.md + Phase_FINAL_close.md + NEEDS.md
(pytanie cyklu + werdykt + co zostalo otwarte). Jeden dokument = jeden cykl, wiec wynikiem
sa CYKLE, nie pliki.

JAK: TF-IDF + cosinus, czysta biblioteka standardowa (zero zaleznosci, zero pobierania modeli).
Normalizacja pod polski: lowercase, usuniecie diakrytykow, przyciecie koncowek fleksyjnych.
Tokeny techniczne (z cyfra lub wielka litera w srodku: psi, M9.1, Q-D1, g^tt) zostaja nietkniete.

ZWRACA SCIEZKI, NIE ZDANIA — nie ma jak zahalucynowac werdyktu.

UZYCIE:
    python tooling/similar.py "oscylony jako stany stacjonarne, masy jako czestosci"
    python tooling/similar.py "kreacja solitonu z prozni" --n 15
    python tooling/similar.py "sciany domenowe Z2" --status closed
    python tooling/similar.py --like op-metric-pair-M911-2026-09-02   # cykle podobne do danego
"""
import io, os, re, sys, math, unicodedata, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "research")
DOCS = ("README.md", "Phase0_balance.md", "Phase_FINAL_close.md", "NEEDS.md")

STOP = set("""
i w z na do nie to sie ze jest a o po za jak od tego ale czy dla przy juz tylko przez
lub oraz jego jej ich nas nam byl byla bylo sa byc moze tez tym tej ten ta te tych tym
co gdy gdyby jesli wiec bez pod nad miedzy wobec wedlug wg poza mimo dzieki
the of and to in is are was were for with that this it as be by on at from or an not
no yes we our can will has have had been which what when where
patrz oraz itd itp np tzn czyli sam sama samo bardzo caly cala cale
""".split())

TECH = re.compile(r"[0-9^_']")          # token techniczny — nie ruszamy
WORD = re.compile(r"[A-Za-zÀ-ſ][A-Za-z0-9À-ſ^_'.\-]*")


def deaccent(s):
    return "".join(c for c in unicodedata.normalize("NFD", s)
                   if unicodedata.category(c) != "Mn")


def norm(tok):
    """Lekka normalizacja fleksyjna: techniczne zostawiamy, slowa przycinamy do rdzenia."""
    if TECH.search(tok):
        return deaccent(tok.lower())
    t = deaccent(tok.lower())
    if len(t) < 3 or t in STOP:
        return None
    return t[:6]                        # crude stemming — dziala dla polskiej fleksji


def tokenize(text):
    out = []
    for m in WORD.finditer(text):
        t = norm(m.group(0))
        if t:
            out.append(t)
    return out


def strip_md(text):
    text = re.sub(r"```.*?```", " ", text, flags=re.S)      # bloki kodu
    text = re.sub(r"\[\[([^\]|]+)(\|[^\]]*)?\]\]", r"\1", text)  # wikilinki
    text = re.sub(r"[|>#*`_~\-]{1,}", " ", text)
    return text


def load_cycles():
    """Zwraca liste (slug, meta, tokeny, ktore_pliki)."""
    cycles = []
    for slug in sorted(d for d in os.listdir(RES) if os.path.isdir(os.path.join(RES, d))):
        d = os.path.join(RES, slug)
        parts, used = [], []
        for name in DOCS:
            p = os.path.join(d, name)
            if os.path.exists(p):
                try:
                    with io.open(p, encoding="utf-8", errors="replace") as f:
                        parts.append(f.read())
                    used.append(name)
                except OSError:
                    pass
        if not parts:
            continue
        raw = "\n".join(parts)
        meta = {}
        m = re.match(r"^---\n(.*?)\n---", parts[0], flags=re.S)
        if m:
            for line in m.group(1).split("\n"):
                mm = re.match(r"^([a-zA-Z_][\w-]*):\s*(.*)$", line)
                if mm:
                    meta[mm.group(1)] = re.sub(r"\s+#.*$", "", mm.group(2)).strip().strip('"').strip("'")
        cycles.append((slug, meta, tokenize(strip_md(raw)), used))
    return cycles


def build_tfidf(cycles):
    N = len(cycles)
    df = collections.Counter()
    tfs = []
    for _, _, toks, _ in cycles:
        tf = collections.Counter(toks)
        tfs.append(tf)
        for t in tf:
            df[t] += 1
    idf = {t: math.log((N + 1.0) / (c + 1.0)) + 1.0 for t, c in df.items()}
    vecs = []
    for tf in tfs:
        v = {t: (1.0 + math.log(c)) * idf[t] for t, c in tf.items()}
        nrm = math.sqrt(sum(x * x for x in v.values())) or 1.0
        vecs.append({t: x / nrm for t, x in v.items()})
    return vecs, idf


def vectorize(tokens, idf):
    tf = collections.Counter(tokens)
    v = {t: (1.0 + math.log(c)) * idf[t] for t, c in tf.items() if t in idf}
    nrm = math.sqrt(sum(x * x for x in v.values())) or 1.0
    return {t: x / nrm for t, x in v.items()}


def cosine(a, b):
    if len(a) > len(b):
        a, b = b, a
    return sum(x * b.get(t, 0.0) for t, x in a.items())


def clean(s, n):
    s = re.sub(r"\s+", " ", s or "").strip()
    return (s[:n - 1] + "…") if len(s) > n else s


def main():
    args = [a for a in sys.argv[1:]]
    n = 10
    status_filter = None
    like = None
    query_parts = []
    i = 0
    while i < len(args):
        if args[i] == "--n" and i + 1 < len(args):
            n = int(args[i + 1]); i += 2
        elif args[i] == "--status" and i + 1 < len(args):
            status_filter = args[i + 1]; i += 2
        elif args[i] == "--like" and i + 1 < len(args):
            like = args[i + 1]; i += 2
        else:
            query_parts.append(args[i]); i += 1

    if not query_parts and not like:
        print(__doc__)
        sys.exit(2)

    cycles = load_cycles()
    vecs, idf = build_tfidf(cycles)
    index = {slug: k for k, (slug, _, _, _) in enumerate(cycles)}

    if like:
        if like not in index:
            cand = [s for s in index if like.lower() in s.lower()]
            if len(cand) == 1:
                like = cand[0]
            else:
                print("Nie znam cyklu '%s'." % like)
                if cand:
                    print("Czy chodzilo o: " + ", ".join(cand[:8]))
                sys.exit(2)
        qv = vecs[index[like]]
        header = "Cykle podobne do: %s" % like
    else:
        qv = vectorize(tokenize(" ".join(query_parts)), idf)
        header = "Zapytanie: %s" % " ".join(query_parts)
        if not qv:
            print("Zadne slowo z zapytania nie wystepuje w korpusie.")
            sys.exit(1)

    scored = []
    for k, (slug, meta, _, used) in enumerate(cycles):
        if like and slug == like:
            continue
        if status_filter and meta.get("folder_status") != status_filter:
            continue
        s = cosine(qv, vecs[k])
        if s > 0:
            scored.append((s, slug, meta, used))
    scored.sort(reverse=True, key=lambda r: r[0])

    print(header)
    print("przeszukano %d cykli (README + Phase0_balance + Phase_FINAL_close + NEEDS)\n" % len(cycles))
    if not scored:
        print("Brak trafien.")
        return
    for rank, (s, slug, meta, used) in enumerate(scored[:n], 1):
        st = meta.get("folder_status", "?")
        verdict = meta.get("verdict") or meta.get("claim_status") or ""
        print("%2d. [%.3f] %-52s %s" % (rank, s, slug[:52], st))
        title = clean(meta.get("title", ""), 118)
        if title:
            print("      %s" % title)
        if verdict:
            print("      werdykt: %s" % clean(verdict, 118))
        print("      research/%s/  (%s)" % (slug, ", ".join(used)))
        print("")


if __name__ == "__main__":
    main()
