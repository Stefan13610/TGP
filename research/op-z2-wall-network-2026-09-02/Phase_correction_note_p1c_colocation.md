---
title: "Correction note (PRZED użyciem skorygowanego wyniku): P1c — komparator GRF uśredniał 2×2 (rejestr komórkowy), a pola ifft2 są rejestrowane WĘZŁOWO: f256[2j]=f128[j] dokładnie dla pola pasmowego |n|≤16 — 0.27 to błąd wygładzenia komparatora, nie grid-zależność (std identyczne 0.0000)"
date: 2026-09-02
type: correction-note
tgp_owner: research/op-z2-wall-network-2026-09-02
status: DOCUMENTED-BEFORE-USE
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase1_output_pre_correction.txt]]"
---

# Korekta 1 — komparator współpołożenia GRF (P1c)

**Pierwotny wynik (zachowany):** `Phase1_output_pre_correction.txt` —
P1c FAIL: max|f128−avg2×2(f256)| = 0.2703 > 0.05, przy IDENTYCZNYM
std (|Δstd| = 0.0000).

**Diagnoza (błąd implementacji TESTU):** konstrukcja GRF przez
`np.fft.ifft2` próbkuje pole ciągłe w punktach węzłowych x_j = j·L/N;
dla pola pasmowego (|n|≤16 < N/2) próbki N=256 i N=128 są WSPÓŁPOŁOŻONE
co drugi węzeł: f256[2j,2k] = f128[j,k] dokładnie (te same
współczynniki, brak aliasingu). Komparator błędnie założył rejestr
komórkowy i uśredniał 4 sąsiednie węzły f256 — różnica 0.27 jest
błędem WYGŁADZENIA komparatora (≈ maks. krzywizna pola × dx²-składki),
nie właściwością konstrukcji. Dowód wewnętrzny w pierwotnym output:
std obu budów identyczne do 4 miejsc (to samo pole).

**Korekta (wyłącznie komparator testu):** porównanie
f128 vs f256[0::2, 0::2] (węzły współpołożone); progi LOCKa (0.05,
2%) NIETKNIĘTE. Zapisano PRZED ponownym uruchomieniem bramki i przed
jakimkolwiek biegiem Phase 2.
