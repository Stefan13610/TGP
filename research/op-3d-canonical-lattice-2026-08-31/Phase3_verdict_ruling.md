---
title: "Phase3_verdict_ruling — ruling kwantyfikatora d przy istnieniu częściowym (zapisany PO Phase 3, PRZED Phase 4)"
date: 2026-08-31
type: verdict-ruling
tgp_owner: research/op-3d-canonical-lattice-2026-08-31
status: RECORDED-BEFORE-PHASE4
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase3_output.txt]]"
---

# Ruling kwantyfikatora d (istnienie częściowe) — PRZED Phase 4

Stan po Phase 2/3 (fakty, bez interpretacji):

- Istnienie tła sieci sc: **wyłącznie d=2π** (oba starty NIESTAŁE-ZBIEŻNE;
  po dedup 2 tła numerycznie rozróżnialne >1e−8, spektralnie identyczne
  do 6 miejsc — efektywnie jedno tło fizyczne). d∈{π, 3.0790}:
  UCIECZKA-g→0; d=3π: NIEZBIEŻNE-siatkowo (A1.0) / KOLAPS (A0.7);
  d=4π: KOLAPS-DO-PRÓŻNI.
- Phase 3 (po odjęciu 3 modów translacyjnych, zidentyfikowanych
  w Γ: overlap=1.0, λ_trans=−1.84e−3=O(h²), Rayleigh zgodny):
  **ω²_min(2π) = −1.674350 (N=48; ZBIEŻNE: |Δ(32→48)|=1.47e−2 ≤
  8.37e−2; argmin Γ)** — dla OBU teł.

**Ruling (wzorzec bloch-chain, pre-rejestrowany w MD §5 i u poprzednika
§7 — stosowany bez zmian):**

1. **PRIMARY:** kwantyfikator „∃d: ω²_min(d)>0" ewaluowany **po tłach
   ISTNIEJĄCYCH** (tu: d=2π). Q-FAIL wymaga: ω²_min<0 zbieżnie dla
   WSZYSTKICH istniejących teł (spełnione już po Phase 3) ORAZ ucieczki
   ≤ 2·t*_izo(3D) w Phase 4 przy obu znakach perturbacji.
   Q-PASS wymagałby tła z ω²_min>0 (zbieżnie) + braku ucieczki do
   3·t*_izo — po Phase 3 wykluczone (żadne tło nie ma ω²>0), więc
   osiągalne werdykty PRIMARY to Q-FAIL albo Q-INCONCLUSIVE.
2. **Strict (raportowany obok, nie zmienia PRIMARY):** kwantyfikator po
   WSZYSTKICH zalockowanych d — punkty bez istniejącego tła nie niosą
   ω² (pytanie o stabilność sieci nieistniejącej jest puste); strict
   werdykt = „negatyw ograniczony do klasy istniejących teł; dla
   d∈{π, d*₁, 3π, 4π} sieć sc kanoniczna nie istnieje w zbadanej klasie
   startów" (negatyw istnienia, nie stabilności).
3. Phase 4 wykonywana dla KAŻDEGO tła ze zbieżnym Phase 3 (tu: oba tła
   d=2π, pełna macierz biegów dla każdego z osobna — mimo że są
   spektralnie identyczne do 6 miejsc; identyczność przebiegów będzie
   raportowana jako diagnostyka dedup).

Zapisano PRZED uruchomieniem Phase 4.
