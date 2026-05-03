"""
Sesja 1 — generator raportu Markdown z _audit_s1_raw.json.
Read-only. Output: AUDIT_RESEARCH_S1.md (raport audytu).
"""

from __future__ import annotations
import json
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "_audit_s1_raw.json"
OUT = HERE / "AUDIT_RESEARCH_S1.md"


def fmt_int(n: int) -> str:
    return f"{n:,}".replace(",", " ")


def main():
    data = json.loads(RAW.read_text(encoding="utf-8"))
    folders = data["folders"]
    n = len(folders)

    # Sekcja 1 — agregaty
    counters = {
        "has_program_md": 0,
        "has_readme": 0,
        "has_findings": 0,
        "has_needs": 0,
        "readme_has_frontmatter": 0,
        "readme_has_tgp_status": 0,
        "is_polluted_74394a8": 0,
        "is_legacy_stay": 0,
    }
    for m in folders:
        for k in counters:
            if m.get(k):
                counters[k] += 1

    heur_counter = Counter(m["heuristic_class"] for m in folders)
    full_conv_count = sum(m["kw_counts"]["FULL_CONVERGENCE"] for m in folders)

    # Top 10 najgęściej zalinkowanych folderów (krytyczność na ewentualne mv)
    by_links = sorted(folders, key=lambda m: m["link_total"], reverse=True)[:15]
    # Foldery 0 linków (bezpieczne fizyczne ruchy)
    zero_link_folders = [m for m in folders if m["link_total"] == 0]

    # Foldery najnowsze + najstarsze
    by_mtime_desc = sorted(folders, key=lambda m: m["mtime_newest"] or 0, reverse=True)
    by_mtime_asc = sorted(folders, key=lambda m: m["mtime_newest"] or 0)

    # Foldery z silnym sygnałem PASS/CLOSED (potencjalni kandydaci core-ready)
    score = lambda m: m["kw_counts"]["PASS"] + m["kw_counts"]["CLOSED"] + m["kw_counts"]["DERIVED"]
    by_pass = sorted(folders, key=score, reverse=True)[:15]

    # Foldery z FALSIFIED/WITHDRAWN (potencjalne archive)
    archive_score = lambda m: m["kw_counts"]["FALSIFIED"] + m["kw_counts"]["WITHDRAWN"] + m["kw_counts"]["OBSOLETE"] + m["kw_counts"]["SUPERSEDED"]
    by_archive_signal = sorted(folders, key=archive_score, reverse=True)[:10]

    # Foldery zatrute 74394a8
    polluted_folders = [m for m in folders if m["is_polluted_74394a8"]]

    # Foldery z subfolders (nested) — closure-style
    with_subfolders = [m for m in folders if m["subfolder_count"] > 0]

    vault_diff = data["vault_diff"]

    # ----- buduj Markdown -----
    lines = []
    lines.append("---")
    lines.append('title: "Sesja 1 — Audyt struktury research/"')
    lines.append('date: 2026-05-03')
    lines.append('type: audit')
    lines.append('status: COMPLETED')
    lines.append('session: S1')
    lines.append(f'generated_at: "{data["generated_at"]}"')
    lines.append('parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"')
    lines.append('related:')
    lines.append('  - "[[meta/AUDYT_TGP_2026-05-01.md]]"')
    lines.append('  - "[[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]]"')
    lines.append('  - "[[INDEX.md]]"')
    lines.append('  - "[[DEPENDENCIES.md]]"')
    lines.append('tags:')
    lines.append('  - audit')
    lines.append('  - research-workflow')
    lines.append('  - session-1')
    lines.append('---')
    lines.append("")
    lines.append("# Sesja 1 — Audyt struktury `research/`")
    lines.append("")
    lines.append("> **Read-only audyt.** Nie nadaje statusów (to Sesja 4), nie modyfikuje folderów,")
    lines.append("> nie przenosi plików. Wynik jest wejściem dla Sesji 2 (szablony) i Sesji 4 (klasyfikacja).")
    lines.append("")
    lines.append(f"- Canonical: `{data['canonical_path']}`")
    lines.append(f"- Liczba folderów badawczych: **{n}**")
    lines.append(f"- Generated: `{data['generated_at']}`")
    lines.append(f"- Skrypty audytu: `meta/research/_audit_s1.py` + `meta/research/_audit_s1_summarize.py`")
    lines.append(f"- Surowe dane: `meta/research/_audit_s1_raw.json` ({fmt_int(RAW.stat().st_size)} B)")
    lines.append("")

    # --- 1. Agregaty ---
    lines.append("## 1. Agregaty")
    lines.append("")
    lines.append(f"| Metryka | Liczba | % z {n} |")
    lines.append("|---|---:|---:|")
    for k, v in counters.items():
        pct = (v * 100 // n) if n else 0
        lines.append(f"| `{k}` | {v} | {pct}% |")
    lines.append(f"| folders with subfolders | {len(with_subfolders)} | {len(with_subfolders)*100//n}% |")
    lines.append(f"| `FULL CONVERGENCE` occurrences (anti-overclaim flag) | {full_conv_count} | — |")
    lines.append("")
    lines.append("**Kluczowe obserwacje agregatów:**")
    lines.append("")
    lines.append(f"- {counters['has_findings']}/{n} folderów ma już `FINDINGS.md` (cel Sesji 5: 100%).")
    lines.append(f"- {counters['has_needs']}/{n} folderów ma już `NEEDS.md` (cel Sesji 5: 100%).")
    lines.append(f"- {counters['readme_has_tgp_status']}/{n} folderów ma już `tgp_status:` w README (cel Sesji 4: 100%).")
    lines.append(f"- {counters['has_program_md']}/{n} folderów ma `program.md` (kanoniczna konwencja 3-fazowa).")
    lines.append(f"- `FULL CONVERGENCE`: {full_conv_count} wystąpień — kandydaty do weryfikacji anty-overclaim w Sesji 8.")
    lines.append("")
    lines.append("> **Nota o liczbie folderów:** plan szacował \"~95\". Audyt liczy")
    lines.append(f"> faktycznie **{n}** podfolderów (różnica = top-level pliki w `research/`,")
    lines.append(f"> sekcja 10: programowe dokumenty + grafy). Liczba {n} jest precyzyjna.")
    lines.append("")

    # --- 2. Heurystyczna pre-klasyfikacja ---
    lines.append("## 2. Heurystyczna pre-klasyfikacja (NIE jest werdyktem)")
    lines.append("")
    lines.append("> **WAŻNE:** Heurystyka skanuje keywordy (`PASS`, `CLOSED`, `FALSIFIED`, …) i nazwy.")
    lines.append("> Nie jest decyzją statusową. Sesja 4 nadaje status z weryfikowalnym źródłem.")
    lines.append(">")
    lines.append("> **Znane ograniczenia heurystyki (kalibracja 2026-05-03):**")
    lines.append(">")
    lines.append("> - Słowo `FALSIFIED` w aktywnym tekście oznacza zazwyczaj **regułę**")
    lines.append(">   falsyfikacji predykcji (np. \"if observed → TGP falsified\"), NIE status")
    lines.append(">   archiwalny folderu. Dlatego archive-flag wymaga `WITHDRAWN`/`OBSOLETE`/")
    lines.append(">   `SUPERSEDED` (te są używane chirurgicznie) lub `FALSIFIED` przy bardzo")
    lines.append(">   niskim PASS.")
    lines.append("> - `candidate-active (heuristic)` jest klasą \"szumową\" — pokrywa wszystko,")
    lines.append(">   co nie pasuje do bardziej specyficznych klas. Sesja 4 musi je faktycznie")
    lines.append(">   sklasyfikować.")
    lines.append("> - `candidate-needs-bridge-or-sandbox` to topic-folders z minimalnym")
    lines.append(">   scaffoldingiem (1–2 .md, ≤3 .py). Część z nich jest realnie active,")
    lines.append(">   tylko jeszcze nie ma `program.md`.")
    lines.append("")
    lines.append("| Wstępna klasa heurystyczna | Liczba |")
    lines.append("|---|---:|")
    for cls, cnt in heur_counter.most_common():
        lines.append(f"| {cls} | {cnt} |")
    lines.append("")

    # Wypisz po klasach
    for cls in sorted(heur_counter.keys()):
        members = sorted([m["name"] for m in folders if m["heuristic_class"] == cls])
        lines.append(f"### 2.x. `{cls}` ({len(members)})")
        lines.append("")
        for name in members:
            lines.append(f"- `research/{name}`")
        lines.append("")

    # --- 3. Top 15 najgęściej linkowanych ---
    lines.append("## 3. Top 15 najgęściej zalinkowanych folderów")
    lines.append("")
    lines.append("> Liczba referencji do ścieżki `research/<folder>` w `INDEX.md` + `DEPENDENCIES.md` +")
    lines.append("> `DEPENDENCIES_REVERSE.md`. Wysoka liczba = drogie do fizycznego przeniesienia.")
    lines.append("> W wariancie minimalnym **niczego nie ruszamy**, ale to lista folderów, których w")
    lines.append("> Sesji 8.5 NIE ARCHIWIZUJEMY (zbyt drogie).")
    lines.append("")
    lines.append("| Folder | INDEX.md | DEPENDENCIES.md | DEP_REV.md | Total |")
    lines.append("|---|---:|---:|---:|---:|")
    for m in by_links:
        lc = m["link_count"]
        lines.append(f"| `{m['name']}` | {lc.get('INDEX.md', 0)} | {lc.get('DEPENDENCIES.md', 0)} | {lc.get('DEPENDENCIES_REVERSE.md', 0)} | **{m['link_total']}** |")
    lines.append("")
    lines.append(f"**Foldery z 0 linkami w INDEX/DEPENDENCIES:** {len(zero_link_folders)} sztuk")
    lines.append("(bezpieczne kandydaty do Sesji 8.5 archiwizacji, jeśli spełnią kryteria 3.4).")
    lines.append("")
    lines.append("<details><summary>Lista folderów z 0 linkami</summary>")
    lines.append("")
    for m in zero_link_folders:
        lines.append(f"- `research/{m['name']}` — {m['heuristic_class']}")
    lines.append("")
    lines.append("</details>")
    lines.append("")

    # --- 4. Najnowsze i najstarsze ---
    lines.append("## 4. Aktywność (mtime ostatniego pliku w folderze)")
    lines.append("")
    lines.append("### 4.1 Top 10 najnowszych (najprawdopodobniej `active`)")
    lines.append("")
    lines.append("| Folder | mtime | klasa | PASS | CLOSED |")
    lines.append("|---|---|---|---:|---:|")
    for m in by_mtime_desc[:10]:
        lines.append(f"| `{m['name']}` | {m['mtime_newest_iso']} | {m['heuristic_class']} | {m['kw_counts']['PASS']} | {m['kw_counts']['CLOSED']} |")
    lines.append("")
    lines.append("### 4.2 Top 10 najstarszych (kandydaci do `archive` / `core-promoted` / legacy)")
    lines.append("")
    lines.append("| Folder | mtime | klasa | PASS | FALSIFIED |")
    lines.append("|---|---|---|---:|---:|")
    for m in by_mtime_asc[:10]:
        lines.append(f"| `{m['name']}` | {m['mtime_newest_iso']} | {m['heuristic_class']} | {m['kw_counts']['PASS']} | {m['kw_counts']['FALSIFIED']} |")
    lines.append("")

    # --- 5. Top PASS/CLOSED ---
    lines.append("## 5. Top 15 folderów z najwyższym sygnałem PASS+CLOSED+DERIVED")
    lines.append("")
    lines.append("> Heurystyka kandydatów do `core-ready` lub już `core-promoted`.")
    lines.append("> **Anty-overclaim**: te liczby same nie wystarczą do nadania statusu — Sesja 4")
    lines.append("> wymaga `source_of_status` cytującego konkretny plik. Polluted folders z 74394a8")
    lines.append("> są wyłączone z dalszej promocji do czasu forward-patch (zob. Sesja 7).")
    lines.append("")
    lines.append("| Folder | PASS | CLOSED | DERIVED | LOCKED | mtime | links | flagi |")
    lines.append("|---|---:|---:|---:|---:|---|---:|---|")
    for m in by_pass:
        flags = []
        if m["is_polluted_74394a8"]:
            flags.append("⚠ POLLUTED-74394a8")
        if m["kw_counts"]["FULL_CONVERGENCE"] > 0:
            flags.append("⚠ FULL_CONVERGENCE")
        if m["is_legacy_stay"]:
            flags.append("legacy-stay")
        flag_str = ", ".join(flags) or "—"
        lines.append(f"| `{m['name']}` | {m['kw_counts']['PASS']} | {m['kw_counts']['CLOSED']} | {m['kw_counts']['DERIVED']} | {m['kw_counts']['LOCKED']} | {m['mtime_newest_iso']} | {m['link_total']} | {flag_str} |")
    lines.append("")

    # --- 6. Archive signal ---
    lines.append("## 6. Sygnał archiwalny (FALSIFIED/WITHDRAWN/OBSOLETE/SUPERSEDED)")
    lines.append("")
    lines.append("> Top 10 folderów z najwyższym sygnałem archiwalnym.")
    lines.append("> **NIE** kwalifikuje automatycznie do `_archive/` — Sesja 8.5 wymaga 2 źródeł")
    lines.append("> per folder + akceptację człowieka per folder.")
    lines.append("")
    lines.append("| Folder | FALSIFIED | WITHDRAWN | OBSOLETE | SUPERSEDED | links | mtime |")
    lines.append("|---|---:|---:|---:|---:|---:|---|")
    for m in by_archive_signal:
        kw = m["kw_counts"]
        if kw["FALSIFIED"] + kw["WITHDRAWN"] + kw["OBSOLETE"] + kw["SUPERSEDED"] == 0:
            continue
        lines.append(f"| `{m['name']}` | {kw['FALSIFIED']} | {kw['WITHDRAWN']} | {kw['OBSOLETE']} | {kw['SUPERSEDED']} | {m['link_total']} | {m['mtime_newest_iso']} |")
    lines.append("")

    # --- 7. Polluted folders ---
    lines.append("## 7. Foldery zatrute incydentem 74394a8")
    lines.append("")
    lines.append("> Z [[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]]. **Plan**: w Sesji 4 te foldery")
    lines.append("> dostają `level: unknown`, `core_compatibility: unknown`,")
    lines.append("> `folder_status: needs-bridge` z notatką. **Nie wolno** awansować ich w Sesji 7")
    lines.append("> bez forward-patch.")
    lines.append("")
    lines.append("| Folder | PASS | FULL_CONVERGENCE | mtime | links |")
    lines.append("|---|---:|---:|---|---:|")
    for m in polluted_folders:
        lines.append(f"| `{m['name']}` | {m['kw_counts']['PASS']} | {m['kw_counts']['FULL_CONVERGENCE']} | {m['mtime_newest_iso']} | {m['link_total']} |")
    lines.append("")

    # --- 8. Closure aggregators (folders with subfolders) ---
    lines.append("## 8. Closure-aggregator i foldery z podfolderami")
    lines.append("")
    lines.append("| Folder | subfolders | nazwy podfolderów |")
    lines.append("|---|---:|---|")
    for m in with_subfolders:
        lines.append(f"| `{m['name']}` | {m['subfolder_count']} | {', '.join('`' + s + '`' for s in m['subfolder_names'])} |")
    lines.append("")

    # --- 9. Vault-root vs canonical diff ---
    lines.append("## 9. Diff vault-rootowego `./research/` vs canonical")
    lines.append("")
    if not vault_diff.get("vault_research_exists"):
        lines.append("Vault-rootowy `./research/` nie istnieje.")
    else:
        only_vault = vault_diff["only_in_vault_root"]
        only_canon = vault_diff["only_in_canonical"]
        common = vault_diff["common_count"]
        lines.append(f"- Foldery wspólne: **{common}**")
        lines.append(f"- Tylko w vault root (`./research/`): **{len(only_vault)}**")
        lines.append(f"- Tylko w canonical (`TGP/TGP_v1/research/`): **{len(only_canon)}**")
        lines.append("")
        if only_vault:
            lines.append("**Tylko w vault root:**")
            for n_ in only_vault:
                lines.append(f"- `{n_}`")
            lines.append("")
        lines.append("**Foldery wspólne — file-level diff:**")
        lines.append("")
        lines.append("| Folder | only_in_vault | only_in_canonical | content_differs |")
        lines.append("|---|---:|---:|---:|")
        for entry in vault_diff["per_folder_diff"]:
            lines.append(f"| `{entry['folder']}` | {len(entry['only_in_vault'])} | {len(entry['only_in_canonical'])} | {len(entry['content_differs'])} |")
        lines.append("")
        lines.append("**Wniosek:** vault-root `./research/` jest **out-of-scope** tego planu (decyzja człowieka 2026-05-02).")
        lines.append("Sesja 1 raportuje stan, ale nie modyfikuje. W przyszłości człowiek może zdecydować:")
        lines.append("(a) usunąć vault-root duplikat, (b) zsynchronizować, (c) zostawić jako staging.")
        lines.append("")

    # --- 10. Top-level pliki w research/ (nie-foldery) ---
    lines.append("## 10. Top-level pliki w `research/` (programowe dokumenty, NIE foldery)")
    lines.append("")
    lines.append("Te pliki **zostają in place** — to programowe dokumenty (status / redirect / new-directions),")
    lines.append("nie foldery badawcze. Plan ich nie dotyka.")
    lines.append("")
    for f in data["top_level_files_in_research"]:
        lines.append(f"- `{f}`")
    lines.append("")

    # --- 11. Pełna lista folderów (kompaktowa tabela) ---
    lines.append("## 11. Pełna lista folderów (kompaktowa tabela)")
    lines.append("")
    lines.append("| # | Folder | klasa heurystyczna | program | phases | YAML | PASS | CLOSED | FALS | mtime | links |")
    lines.append("|---:|---|---|:-:|---:|:-:|---:|---:|---:|---|---:|")
    for i, m in enumerate(folders, 1):
        prog = "✓" if m["has_program_md"] else "—"
        n_phase = len(m["phase_files"])
        yaml = "✓" if m["readme_has_frontmatter"] else "—"
        kw = m["kw_counts"]
        lines.append(
            f"| {i} | `{m['name']}` | {m['heuristic_class']} | {prog} | {n_phase} | {yaml} | "
            f"{kw['PASS']} | {kw['CLOSED']} | {kw['FALSIFIED']} | {m['mtime_newest_iso']} | {m['link_total']} |"
        )
    lines.append("")

    # --- 12. Pytania do człowieka (po Sesji 1) ---
    lines.append("## 12. Pytania pozostające otwarte po Sesji 1")
    lines.append("")
    lines.append("Decyzje 1–6 z [[meta/PLAN_RESEARCH_WORKFLOW_v1.md#9-decyzje-człowieka-zatwierdzone-2026-05-02]]")
    lines.append("już są zatwierdzone. Pytania, które ten audyt **dodatkowo** odsłonił:")
    lines.append("")
    lines.append("**Q1.** Vault-root `./research/` (11 folderów) ma 100% pokrycie z canonical (zob. sekcja 9).")
    lines.append("Czy człowiek chce:")
    lines.append("- (a) zachować bez zmian (out-of-scope na zawsze),")
    lines.append("- (b) zaplanować kasowanie po zakończeniu workflow (Sesja 9.x dodatkowo),")
    lines.append("- (c) zsynchronizować z canonical jako mirror (np. via skrypt po commitach)?")
    lines.append("")
    lines.append("**Q2.** Foldery z subfolderami (sekcja 8) — czy każde podzamknięcie")
    lines.append("(np. `closure_2026-04-26/sigma_ab_pathB/`) traktujemy jako osobny **byt**")
    lines.append("z własnym YAMLem, NEEDS, FINDINGS, czy tylko rodzic dostaje pełną warstwę,")
    lines.append("a dzieci są lekkie?")
    lines.append("")
    lines.append("**Q3.** Top-level pliki w `research/` (sekcja 10): czy program-doc dostaje też")
    lines.append("YAML `tgp_status` (`folder_status: program-doc`, `kind: program-doc`),")
    lines.append("czy zostaje bez frontmattera?")
    lines.append("")
    if counters['has_findings'] + counters['has_needs'] > 0:
        lines.append("**Q4.** Foldery z `has_findings: true` lub `has_needs: true` już istnieją w niektórych")
        lines.append(f"miejscach ({counters['has_findings']} z FINDINGS, {counters['has_needs']} z NEEDS).")
        lines.append("Czy w Sesji 5 agent **nadpisuje** istniejące, czy tylko uzupełnia brakujące?")
        lines.append("Domyślna sugestia: **uzupełnia brakujące, istniejące zostawia + flag `pre_existing: true`**.")
        lines.append("")
    else:
        lines.append("**Q4.** *(nieaktualne)* — żaden folder nie ma jeszcze `FINDINGS.md`/`NEEDS.md`, więc")
        lines.append("Sesja 5 zaczyna od zera. Brak konfliktu z istniejącymi plikami.")
        lines.append("")

    # --- 13. Co dalej ---
    lines.append("## 13. Następne kroki")
    lines.append("")
    lines.append("Po akceptacji odpowiedzi Q1–Q4:")
    lines.append("")
    lines.append("- **Sesja 2**: szablony README/NEEDS/FINDINGS + 3 case-study previews.")
    lines.append("- **Sesja 3**: warstwa `meta/research/` + `meta/core/` + fizyczne `_sandbox`/`_archive`.")
    lines.append(f"- **Sesja 4**: klasyfikacja YAML w {n} folderach (batch po 10 → ~{(n + 9) // 10} batchy).")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("*Raport wygenerowany automatycznie. Surowe dane: [[_audit_s1_raw.json]].*")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"OK wrote {OUT} ({OUT.stat().st_size} bytes, {len(lines)} lines)")


if __name__ == "__main__":
    main()
