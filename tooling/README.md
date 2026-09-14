# Narzedzia i toolchain TGP

Narzedzia i skrypty pomocnicze projektu: toolchain numeryczny, skrypty weryfikacyjne, zrzuty wykresow, skrypty pomocnicze LaTeX i plik blokady zaleznosci.

## Zawartosc
- `NARZEDZIA_I_TOOLCHAIN.md` — opis toolchainu (kompilatory, zaleznosci, workflow)
- `requirements.lock` — zablokowane wersje zaleznosci Pythona
- `count_sections.py` — skrypt liczacy sekcje w zrodlach LaTeX
- `fix_wikilinks.py` — skrypt naprawiajacy wikilinki (Obsidian)
- `build_cycles_index.py` — generuje `CYCLES.tsv` w roocie (mapa 263 cykli, 1 linia = 1 cykl; tanszy zamiennik INDEX.md dla agentow); `--check` = exit 1 gdy nieaktualny
- `similar.py` — „czy my to juz liczylismy?": wyszukiwarka podobienstwa po cyklach (TF-IDF + cosinus, zero zaleznosci). Indeksuje README + Phase0_balance + Phase_FINAL_close + NEEDS kazdego cyklu; zwraca sciezki i werdykty, nie zdania. Tryby: zapytanie tekstowe, `--like <slug>`, `--status <wartosc>`, `--n <ile>`
- `normalize_folder_status.py` — ujednolica `folder_status` we frontmatterach README cykli do 6 wartosci (wip/locked/paused/parking/closed/legacy); dry-run domyslnie, `--apply` zapisuje; oryginal zachowywany w `legacy_status`
- `scripts/` — podfolder z pelnym zestawem skryptow numerycznych (F_alpha, Koide, substrate, ERG, BH shadow, cosmology, gw, gauge, particles, profiles, stability, ...); ma wlasny `README.md`
- `plots/` — wygenerowane wykresy PNG (einstein_emergence, alpha_eff_chi2, gl_phase_transition, growth_factor_tgp, phi0_cross_verification, O16_kink_mass_hierarchy, O22_entropy_cosmology_scan)
