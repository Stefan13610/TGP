---
title: "Phase_FINAL_close — zamknięcie: ODPOWIEDŹ N1 — kieszeń przedmetryczna re-metryzuje się w skończonym τ (Q-PKT-A-HEAL-ALL, τ_close≈4.6·R, front); przestrzał znaku daje długowieczną POWŁOKĘ przedmetryczną (Q-PKT-B-ZANIK wg curvature flow: τ_life≈0.98·R², R²=0.9997); trwałość przedmetryczna istnieje i jest TOPOLOGICZNA — ściany Z2 to stabilne struktury ψ≈0 (Q-PKT-C-PASS-SHEET; rdzeń ψ~2e−5, pas ψ<0.218 szeroki 2.43) niewidoczne dla poziomu 1; kolaps kieszeni NIE zostawia zlokalizowanego obiektu (Q-PKT-D CLEAN-WALL — pre-rejestrowany pozytyw autora nie zaszedł)"
date: 2026-09-02
type: phase-final-close
tgp_owner: research/op-premetric-pocket-level0-2026-09-02
status: CLOSED
verdict: "Phase 1: PASS (P1a bare 0 węzłów metrycznych; P1b single(1.4) zanika — ciągłość z op-bare-substrate-genesis; P1c H nierosnące; P1d profil 1D ściany zbieżny po korekcie 1 pomiaru: Δszer 0.012–0.027%, ΔΦ_min/Φ_bar=4.5e−4, σ_w=0.42977). Q-PKT-A: HEAL-ALL — kieszeń goła (R=4/8/16) zasklepia się frontem próżni w τ_close=8/24/63; skalowanie LINIOWE τ≈4.63·R (R²=0.998 > kwadratowe 0.991) ⟹ prędkość frontu v≈0.216. Q-PKT-B (centralne): ZANIK — kropla antyfazowa (−s*, powłoka przedmetryczna Z2) żyje τ_life=9/52/243, skalowanie KWADRATOWE τ≈0.98·R² (R²=0.9997 — curvature flow ścian) ⟹ powłoka jest DŁUGOWIECZNA (τ_life/τ_close ∝ R: dla R=16 już 3.9×), ale znika; brak stabilizacji rozmiaru. Q-PKT-C: PASS-SHEET — proste ściany Z2 trwałe na obu siatkach (tail=1.0000, H monotone; pas przedmetryczny Φ<0.30: 1024 węzłów N=128 / 7.81% N=256, szerokość na ścianę ~2.0/2.5 vs 1D 2.43 — rozdzielczość progu, klasyfikacja zgodna); rdzeń Φ_min→2.7e−5 (1D, ψ~2e−5). Q-PKT-D: CLEAN-WALL (R=4 i 8) — kieszeń osadzona na ścianie zasklepia się DOKŁADNIE do ściany bazowej (excess=0.000 vs próg 2.4e−4); pre-rejestrowany pozytyw autora (koralik) NIE zaszedł. ŁĄCZNA ODPOWIEDŹ N1 (drzewo §5, gałąź spójna): kolaps ψ→0 poziomu 1 kończy się na poziomie 0 re-metryzacją (bounce) w τ∝R; jeśli kolaps przestrzeliwuje znak (obserwowany przestrzał poziomu 1), pozostawia powłokę przedmetryczną o życiu ∝R²; JEDYNA trwała forma przedmetryczna to topologia Z2 (ściany) — strukturalnie NIEWIDOCZNA w opisie efektywnym ψ=Φ/Φ* (ślepym na znak s)."
anti_lakatos_lock: PRESERVED
tags: [premetric-pocket, level0, z2-walls, curvature-flow, bounce, topological-premetric-sheet, n1-closure, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_correction_note_p1d_band.md]]"
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-dynamics-class-M911-2026-09-02/Phase_FINAL_close.md]]"
  - "[[../op-bare-substrate-genesis-2026-07-04/Phase_FINAL_close.md]]"
---

# Phase FINAL — zamknięcie cyklu op-premetric-pocket-level0

**Status: CLOSED-EXECUTED (jedna sesja: LOCK → Phase 1 (+korekta 1
testu) → Phase 2 → zamknięcie).** Kryteria LOCKa stosowane DOSŁOWNIE;
zero zmian progów/modelu po starcie; jedna korekta implementacji
POMIARU (nie fizyki), udokumentowana przed użyciem.

## 0. Werdykty

| Pytanie | Werdykt | Jedno zdanie |
|---|---|---|
| P1 (ciągłość + 1D) | **PASS** (po korekcie 1) | bare/single reprodukują G2/G3 poprzednika; ściana 1D zbieżna (σ_w=0.42977, Δ<0.03%) |
| **Q-PKT-A** (kieszeń goła) | **HEAL-ALL** | τ_close = 8/24/63 (R=4/8/16), liniowo τ≈4.63·R (front, v≈0.216) |
| **Q-PKT-B** (kropla antyfazowa) | **ZANIK** (curvature flow) | τ_life = 9/52/243, τ≈0.98·R² (R²=0.9997) — powłoka przedmetryczna długowieczna, nie wieczna |
| **Q-PKT-C** (ściany Z2) | **PASS-SHEET** | proste ściany trwałe (tail=1.0000, obie siatki) — stabilne struktury ψ≈0 |
| **Q-PKT-D** (pozostałość) | **CLEAN-WALL** | excess = 0.000 (próg 2.4e−4) — kolaps nie zostawia obiektu; pozytyw autora NIE zaszedł |

## 1. Model i wejścia

Model VERBATIM z `op-bare-substrate-genesis` (CLOSED; cytat w LOCKu §1):
s(x) 2D N=128 L=64 dx=0.5 dt=0.02, Φ=s², ds/dτ=κ∇²s−V′(s),
V=0.5as²−(b/3)|s|³+0.25cs⁴, a=0.5 b=1.6 c=1.0, κ=0.5; ε=0.30,
A_min=4/128²; s*=1.1741657 (Φ*=1.3787), s_bar=0.4258 (Φ_bar=0.1813).
Wejścia cyklu: seed=20260905; R∈{4,8,16}; δ=1.07; steps=30000 (τ=600);
kontrola pinningu N=256 dx=0.25 dt=0.005. Mapowanie: ψ=Φ/Φ*;
pas przedmetryczny Φ<0.30 ⟺ ψ<0.218.

## 2. Wyniki liczbowe

- **A (kieszeń goła w +s*):** wszystkie zasklepione do zera węzłów
  gołych; τ_close(R): 8.0 / 24.0 / 63.0; fit liniowy a=4.625
  (R²=0.99781) > kwadratowy (R²=0.99070) ⟹ zasklepianie FRONTEM
  (prędkość ~0.216 j./j.τ), nie dyfuzją.
- **B (kropla −s* w +s* = zamknięta powłoka przedmetryczna):**
  τ_life(R): 9 / 52 / 243; fit kwadratowy a=0.9807, R²=0.99972
  (liniowy 0.97462) ⟹ klasyczny curvature flow ściany Z2
  (dR²/dτ≈const). Powłoka NIE stabilizuje rozmiaru; żyje ∝R²
  (dla R=16: τ_life=243 vs τ_close=63 kieszeni gołej — 3.9× dłużej,
  stosunek rośnie ∝R).
- **C (proste ściany Z2, tor):** trwałe: bare_area stała
  (tail=1.0000; N=128: 1024 węzły = 2 pasy szer. ~2.0; N=256: 7.81%
  = szer. ~2.5/ścianę vs 1D 2.43 — różnica to rozdzielczość progu ε
  na siatce, klasyfikacja zgodna); H monotone; s_min=−s* (rdzenie
  antyfazowe nienaruszone). Profil 1D: Φ_min=2.7e−5 (ψ≈2e−5!),
  pas Φ<ε szeroki 2.4300, pas Φ<Φ_bar 1.9103, napięcie σ_w=0.429771.
- **D (kieszeń na ścianie):** R=4 i R=8 — stan końcowy IDENTYCZNY
  z bazową ścianą (excess=0.000e+00; próg A_min=2.441e−4);
  tail=1.0000. Brak koralika.

## 3. ODPOWIEDŹ N1 (wprost; klasa zbadana: 2D, relaksacja w τ)

Kolaps dipa słabopolowego poziomu 1 (ψ→0, op-dynamics-class-M911) ma
na poziomie 0 następujące rozstrzygnięcie:

1. **Bounce:** kieszeń przedmetryczna bez zmiany znaku re-metryzuje
   się w skończonym τ ∝ R — front próżni metrycznej zamyka ją;
   nic nie zostaje.
2. **Przestrzał znaku (naturalny odpowiednik obserwowanego przestrzału
   poziomu 1) ⟹ powłoka przedmetryczna:** zamknięta ściana Z2
   o życiu τ ≈ 0.98·R² — obiekt DŁUGOWIECZNY (kwadratowo w rozmiarze),
   lecz kurczący się; to jest poziom-0 realizacja „długowiecznej
   kieszeni przedmetrycznej" z NEEDS N1.
3. **Trwałość przedmetryczna = topologia Z2:** jedyne wieczne struktury
   ψ≈0 to ściany domenowe (PASS-SHEET) — wymagają nietrywialnej
   topologii znaku s na dużą skalę, a opis efektywny ψ jest na znak
   ŚLEPY: sieć ścian jest niewidoczna dla poziomu 1. Wniosek
   strukturalny: „próżnia" ψ≡1 poziomu 1 może maskować sieć
   przedmetrycznych arkuszy substratu.
4. **Kreacja przez kolaps:** NIE zaszła (D: CLEAN-WALL) — kolaps
   kieszeni w tej klasie nie produkuje zlokalizowanych obiektów;
   dla hipotezy usera nośnikiem pozostaje długowieczność powłoki (2),
   nie trwała kreacja.

## 4. Korekty / incydenty / higiena (anti-Lakatos)

- ✓ **Korekta 1** (`Phase_correction_note_p1d_band.md`, zapisana PRZED
  ponownym biegiem bramki i przed Phase 2): pomiar szerokości pasa bez
  interpolacji podsiatkowej (błąd testu O(dx)≈2% > próg 1%);
  pierwotny FAIL zachowany (`Phase1_output_pre_correction.txt`);
  kryterium LOCKa nietknięte; wzorzec incydentu QF-4c poprzednika.
- ✓ Model/progi/promienie/seed nietknięte po starcie; kontrola N=256
  wykonana (C — obowiązkowa; A/B/D bez pozytywów wymagających);
  H_Γ nierosnące we wszystkich biegach; boundary_contact: 0.
- ✓ Katalogi innych cykli tylko odczyt; rdzeń .tex/STATE/git
  nietykane; pełne ścieżki bez `cd`; `ls` po zapisach. τ ≠ czas
  fizyczny; zero claimów grawitacyjnych/obserwacyjnych.
- Środowisko: CPython 3.14.2, numpy 2.4.3.

## 5. Ograniczenia klasy (pre-rejestrowane w LOCKu)

2D (poziom-0 model dziedziczony); relaksacja w τ (dynamika selekcji
substratu — korpusowa dla poziomu 0); ±s* dokładnie degenerowane
(brak tiltu Z2). W 3D powłoka B = zamknięta membrana (też curvature
flow, τ∝R² oczekiwane strukturalnie — do weryfikacji, NEEDS).

## 6. Pliki cyklu

`Phase0_balance.md` (LOCK) · `Phase1_gate.py` → `Phase1_output.txt`
(+ `Phase1_output_pre_correction.txt`) ·
`Phase_correction_note_p1d_band.md` · `Phase2_pocket.py` →
`Phase2_output.txt` + `Phase2_results/` (json/npz per scenariusz
+ kontrola) + `Phase2_run.log` + `Phase2_ctrl.log` · `NEEDS.md`
(user-gated) · `README.md`.

## 7. Mapowanie na drzewo LOCKa §5

Gałąź spójna: **A-HEAL-ALL ∧ B-ZANIK ∧ C-PASS-SHEET (+ D-CLEAN-WALL)**
⟹ „kieszeń re-metryzuje się w skończonym τ (bounce), przestrzał znaku
daje powłokę przedmetryczną o życiu τ_life(R), trwałość przedmetryczna
= topologia Z2 niewidoczna dla poziomu 1; NEEDS: dopisek core
(user-gate) + pytanie o 3D/sieć ścian" — dosłownie; [[NEEDS.md]].
