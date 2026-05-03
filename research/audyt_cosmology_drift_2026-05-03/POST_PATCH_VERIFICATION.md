---
title: "POST_PATCH_VERIFICATION — checklist akceptacyjny"
date: 2026-05-03
parent: "[[README.md]]"
type: verification
---

# POST_PATCH_VERIFICATION — checklist akceptacyjny

## Smoke tests (Phase 5.1)

| # | Test | Expected |
|---|---|:---:|
| 1 | `grep "parametrem dopasowania, nie predykcj" TGP/TGP_v1/core/sek05*/*.tex` | TRUE (preserved) |
| 2 | `grep "rem:Lambda-post-uv3" TGP/TGP_v1/core/sek05*/*.tex` | TRUE (new remark) |
| 3 | `grep "Status update 2026-05-03\\|Status UPDATE 2026-05-03" TGP/TGP_v1/research/desi_dark_energy/README.md` | TRUE |
| 4 | `grep "12.5 Post-M10 addenda" TGP/TGP_v1/research/op-cosmology-closure/M10_R_results.md` | TRUE |
| 5 | `grep "F1.5\|F11\|F12" TGP/TGP_v1/research/op-cosmology-closure/M10_R_results.md` | TRUE (3 matches) |
| 6 | `grep "DESI DR2.*2503.14738" TGP/TGP_v1/research/desi_dark_energy/README.md` | TRUE |
| 7 | `grep "A.13.*Cosmology drift remediation" TGP/TGP_v1/research/closure_2026-04-26/KNOWN_ISSUES.md` | TRUE |

## Manual review checklist (Phase 5.2)

### sek05.tex
- [ ] `\rem:Lambda-post-uv3` wstawione przed `\rem:Lambda-quantitative`
- [ ] `prob:Lambda` flag zmieniony na "ROZWIĄZANY (post-2026-05-02)"
- [ ] Footnote G.0 v2.0 dodany do `\eq:U-phi-def`
- [ ] Stale-DR1 uwaga dodana do `\rem:wz-quantitative`
- [ ] Stare zdanie "Φ_0 jest parametrem dopasowania" PRESERVED (non-breaking)
- [ ] Diff ≤ 80 dodanych linii, 0 usuniętych linii

### desi_dark_energy/README.md
- [ ] Sekcja "Status update 2026-05-03" wstawiona po YAML, przed H1
- [ ] Sekcja 9 (post-cascade falsifikatory) dodana
- [ ] Tabele 1+4 mają kolumnę DR2
- [ ] Wszystkie cytowania `[[..]]` rozwiązują się
- [ ] Diff ≤ 60 dodanych linii, 0 usuniętych

### M10_R_results.md
- [ ] Sekcja 12.5 (Post-M10 addenda) wstawiona po sekcji 12, przed sekcją 13
- [ ] F1.5 / F11 / F12 wpisy dodane
- [ ] Cross-link do `[[../op-uv3-phi0-renormalization/]]`, `[[../op-gamma1-phi-eff-anchor-resolution/]]`, etc.
- [ ] YAML `last_yaml_update: "2026-05-03"`

### M10_5_results.md
- [ ] M10.5.5 ma uwagę o DR1 vs DR2 numbers
- [ ] Tabela sekcji 3 ma 2 nowe wiersze (Phase 3.E, UV.3)
- [ ] YAML `last_yaml_update: "2026-05-03"`

### KNOWN_ISSUES.md
- [ ] Entry A.13 dodany, oznakowany ✅ CLOSED 2026-05-03

### Auxiliary
- [ ] cosmo_tensions/README.md ma sekcję "Status post-cascade 2026-05-02"
- [ ] hubble_tension/README.md ma sekcję "Status post-cascade 2026-05-02"

### main.pdf rebuild
- [ ] `pdflatex main` exit code 0
- [ ] PDF diff ≤ 200 linii vs pre-patch (zmiany ograniczone do sek05)
- [ ] Brak `??` w PDF (wszystkie nowe `\ref{}` rozwiązane)

## Color-graded final status

Po wszystkich check:

| Plik | Pre-patch | Post-patch |
|---|:---:|:---:|
| sek05.tex | 🔴 RED | 🟢 GREEN |
| desi_README.md | 🔴 RED | 🟢 GREEN |
| M10_R_results.md | 🟡 YELLOW | 🟢 GREEN |
| M10_5_results.md | 🟡 YELLOW | 🟢 GREEN |

## Final verdict

✅ Audit closed po wszystkich smoke testach + manual review + main.pdf rebuild
   bez nowych errors.

⚠️ Audit hold-open jeśli którykolwiek check fails — log w [[patch_log.md]]
   i decyzja: rollback (selective git revert) lub rework patcha.
