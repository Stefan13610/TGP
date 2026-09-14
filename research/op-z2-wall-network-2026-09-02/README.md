---
title: "op-z2-wall-network — sieć ścian Z2 z genezy wielodomenowej: mozaika domen ±s*, prawo grubienia, konfiguracje topologicznie trwałe (ukryta dla ψ struktura wielkoskalowa substratu)"
folder_status: closed
date: 2026-09-02
type: research-cycle
status: CLOSED
verdict: "P1 PASS (po korekcie 1 komparatora GRF); Q-NET-A-PASS-MOSAIC (10/10 genez lockuje w mozaikę obu znaków); Q-NET-B-INCONCLUSIVE (sieci grubieją przez skalę pudła L=64 przed τ=1000 — wykładnik wymaga większego L); Q-NET-C-PASS-TOPOLOGICAL (s=20260908: trwałe paski nawinięte na torus, L_wall=128.0, tails=1.0000, kontrola N=256 niemal identyczna; łącznie 3/10 genez → stan nawinięty ~1/3). Trwała sieć przedmetryczna z losowej genezy ISTNIEJE; widok ψ: 7.81% objętości w pasach ψ<0.218, etykieta znaku niewidoczna."
tgp_owner: research/op-z2-wall-network-2026-09-02
authorization: "User 2026-09-06: „ok działaj z N1" (N1 z NEEDS op-premetric-pocket-level0)"
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-premetric-pocket-level0-2026-09-02/Phase_FINAL_close.md]]"
  - "[[../op-bare-substrate-genesis-2026-07-04/Phase_FINAL_close.md]]"
---

# op-z2-wall-network (2026-09-02)

Realizacja N1 z op-premetric-pocket-level0: skoro ściany Z2 są trwałymi
strukturami przedmetrycznymi, a geneza substratu nie preferuje znaku,
to losowa geneza powinna produkować mozaikę domen ±s* z siecią arkuszy
ψ≈0 — niewidoczną w opisie efektywnym. Model VERBATIM
z op-bare-substrate-genesis; starty: GRF pasmowy grid-niezależny,
A∈{0.6,1.0} × 5 seedów; τ_max=2000.

**Pytania:** A — czy powstaje mozaika obu znaków (vs jeden znak /
faza goła)? B — czy L_wall(τ) grubieje z wykładnikiem −1/2 (curvature
flow)? C — czy część genez kończy w trwałej konfiguracji nawiniętej
na torus (PASS-TOPOLOGICAL), czy wszystko grubieje do jednej domeny?

Kryteria: [[Phase0_balance.md]].

## Log faz

- 2026-09-02: Phase 0 LOCK zapisany (autoryzacja „ok działaj z N1")
  + Amendment A1 pre-code (czynnik stereologiczny 4/π estymatora
  manhattańskiego).
- 2026-09-02: Phase 1 — pierwotnie P1c FAIL (komparator GRF uśredniał
  2×2 przy rejestrze węzłowym ifft2 — błąd testu); korekta 1
  udokumentowana przed użyciem (`Phase_correction_note_p1c_colocation.md`,
  pierwotny output zachowany); po korekcie **PASS** (paski L_wall=128.0
  exact; kropla 64.0=8R; współpołożenie 0.0000; H monotone).
- 2026-09-02: Phase 2 — 10 biegów głównych + kontrola N=256:
  **A-PASS-MOSAIC** (10/10); **B-INCONCLUSIVE** (4/5 okien nieważnych —
  skala pudła; ważne: slope −0.249 z plateau pasków);
  **C-PASS-TOPOLOGICAL** (s=20260908 A=1.0: WOUND-STRIPES, kontrola
  N=256 zgodna: +0.516/−0.414, Lw=128.0, tails 1.0000; deskryptywnie
  3/10 nawiniętych). `Phase2_output.txt`, `Phase2_results/`.
- 2026-09-02: Zamknięcie: `Phase_FINAL_close.md` (odpowiedź N1 wprost)
  + `NEEDS.md` (user-gated: skalowanie pudła; 3D arkusze; dopisek core;
  obserwable poziomu 1 na sieci). **CLOSED.**
