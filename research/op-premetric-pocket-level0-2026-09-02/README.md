---
title: "op-premetric-pocket-level0 — fizyka kieszeni przedmetrycznej ψ→0 na poziomie 0 (N1 po op-dynamics-class-M911): los kieszeni fazy gołej, kropla antyfazowa Z2, ściany domenowe jako trwałe struktury przedmetryczne, pozostałość po kolapsie"
folder_status: closed
date: 2026-09-02
type: research-cycle
status: CLOSED
verdict: "P1 PASS (po korekcie 1 pomiaru); Q-PKT-A-HEAL-ALL (τ_close≈4.63·R, front v≈0.216); Q-PKT-B-ZANIK (τ_life≈0.98·R², R²=0.9997 — curvature flow; powłoka przedmetryczna długowieczna ∝R²); Q-PKT-C-PASS-SHEET (ściany Z2 trwałe: rdzeń ψ~2e−5, pas ψ<0.218 szer. 2.43, tail=1.0 obie siatki); Q-PKT-D CLEAN-WALL (excess=0 — kolaps nie zostawia obiektu). ODPOWIEDŹ N1: bounce / powłoka z przestrzału znaku / trwałość przedmetryczna = topologia Z2 niewidoczna dla poziomu 1."
tgp_owner: research/op-premetric-pocket-level0-2026-09-02
authorization: "User 2026-09-06: „ok, działaj z N1""
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-dynamics-class-M911-2026-09-02/NEEDS.md]]"
  - "[[../op-bare-substrate-genesis-2026-07-04/Phase_FINAL_close.md]]"
---

# op-premetric-pocket-level0 (2026-09-02)

Realizacja N1 z op-dynamics-class-M911: kolaps dipa słabopolowego
(poziom 1, dynamika bezwładna M9.1'') kończy się w ψ→0, gdzie model
efektywny traci ważność (B→0, w→0). Ten cykl pyta, co robi z taką
kieszenią SUBSTRAT (poziom 0, model dziedziczony VERBATIM
z op-bare-substrate-genesis: Φ=s², faza goła Φ=0 metastabilna,
próżnia s*=1.174, gradient flow w τ).

**Klucz strukturalny:** substrat ma Z2 (s→−s), poziom 1 (ψ=Φ/Φ*) jest
ślepy na znak; każda ścieżka −s*→+s* przechodzi przez s=0 ⟹ ściany
domenowe mają z konieczności rdzeń przedmetryczny Φ=0. Kolaps
z przestrzałem = lokalne przejście do przeciwnej studni.

**Pytania:** A — czy kieszeń goła zawsze się zasklepia (τ_close(R))?
B (centralne) — los kropli antyfazowej −s* (powłoka przedmetryczna):
zanik τ_life(R) / trwanie / fragmentacja? C — czy płaska ściana Z2
jest trwałą strukturą przedmetryczną (pas Φ<0.30 stabilny)?
D — czy kieszeń osadzona na ścianie zostawia po zasklepieniu
zlokalizowany nadmiar („koralik" — pre-rejestrowany pozytyw autora)?

Kryteria: [[Phase0_balance.md]].

## Log faz

- 2026-09-02: Phase 0 LOCK zapisany (autoryzacja „ok, działaj z N1").
- 2026-09-02: Phase 1 — pierwotnie P1d FAIL (błąd implementacji POMIARU
  szerokości pasa: bez interpolacji podsiatkowej, błąd O(dx)≈2% > próg
  1%); korekta 1 udokumentowana PRZED użyciem
  (`Phase_correction_note_p1d_band.md`, pierwotny output zachowany);
  po korekcie **PASS** (bare/single reprodukują G2/G3; ściana 1D:
  Φ_min=2.7e−5, pas Φ<0.30 szer. 2.4300, σ_w=0.429771, zbieżność
  0.01–0.03%).
- 2026-09-02: Phase 2 — **A-HEAL-ALL** (τ_close=8/24/63, liniowo 4.63·R);
  **B-ZANIK** (τ_life=9/52/243, kwadratowo 0.98·R², R²=0.9997);
  **C-PASS-SHEET** (tail=1.0000, N=128 i N=256; kontrola pinningu
  zgodna); **D-CLEAN-WALL** (excess=0.000, R=4 i 8). `Phase2_output.txt`,
  `Phase2_results/`.
- 2026-09-02: Zamknięcie: `Phase_FINAL_close.md` (odpowiedź N1 wprost)
  + `NEEDS.md` (user-gated: sieć ścian z genezy = ukryta struktura
  wielkoskalowa; dopisek core o ślepocie ψ na Z2; wersja 3D; most
  powłoka→poziom 1). **CLOSED.**
