---
title: "Phase0_balance — LOCK: sieć ścian Z2 z genezy wielodomenowej (N1 op-premetric-pocket-level0) — czy losowa geneza substratu produkuje trwałą mozaikę domen ±s* z siecią arkuszy przedmetrycznych; prawo grubienia; konfiguracje topologicznie trwałe na torze"
date: 2026-09-02
type: phase0-lock
tgp_owner: research/op-z2-wall-network-2026-09-02
status: PHASE0-LOCKED
computations_performed: ZERO
authorization: "User 2026-09-06: „ok działaj z N1" (N1 z NEEDS op-premetric-pocket-level0: geneza wielodomenowa → statystyka sieci ścian — gęstość, grubienie, konfiguracje stabilne na torze)"
anti_lakatos_lock: ACTIVE
related:
  - "[[README.md]]"
  - "[[../op-premetric-pocket-level0-2026-09-02/Phase_FINAL_close.md]]"
  - "[[../op-premetric-pocket-level0-2026-09-02/NEEDS.md]]"
  - "[[../op-bare-substrate-genesis-2026-07-04/Phase0_balance.md]]"
---

# Phase 0 — LOCK cyklu `op-z2-wall-network`

**ZERO obliczeń wykonanych przed zapisaniem tego dokumentu.**

## 0. Pytanie i dziedziczona diagnoza

Poprzednik (op-premetric-pocket-level0): ściany Z2 substratu są TRWAŁYMI
strukturami przedmetrycznymi (rdzeń ψ~2e−5, pas ψ<0.218 szer. 2.43,
PASS-SHEET), zamknięte powłoki giną przez curvature flow (τ≈0.98·R²),
a poziom 1 (ψ) jest ślepy na znak s. Geneza substratu
(op-bare-substrate-genesis) lockuje kolektywnie i NIE preferuje znaku.
Wniosek do przetestowania: losowa geneza powinna produkować MOZAIKĘ
domen ±s* przedzieloną trwałą, grubiejącą siecią arkuszy
przedmetrycznych — pierwszy kandydat programu na ukrytą (dla ψ)
wielkoskalową strukturę substratu.

**Pytania binarne:**
- **Q-NET-A (mozaika):** czy zalockowany start losowy (GRF pasmowy,
  amplituda nadbarierowa) lockuje w mozaikę OBU znaków (nie: jeden znak
  wygrywa od razu / wszystko do fazy gołej)?
- **Q-NET-B (prawo grubienia):** czy długość sieci ścian L_wall(τ)
  grubieje potęgowo z wykładnikiem curvature-flow −1/2 (2D,
  parametr niezachowany)?
- **Q-NET-C (trwałość topologiczna):** czy część genez kończy
  w konfiguracji TRWAŁEJ (paski nawinięte na torus — sieć
  przedmetryczna, która NIE znika), czy każda geneza grubieje do
  jednej domeny (sieć zawsze przejściowa)?

Uwaga zalockowana: FAIL-y są pełnoprawnymi wynikami (jeden znak wygrywa
⟹ mozaika nie powstaje; wszystkie geneza → jedna domena ⟹ trwałość
przedmetryczna wymaga topologii wpisanej w warunki, nie genezy).
ZERO claimów kosmologicznych/obserwacyjnych; τ ≠ czas fizyczny.

## 1. Model ZAMKNIĘTY

**Model substratu VERBATIM** (łańcuch dziedziczenia: op-bare-substrate-
genesis LOCK §1–2, cytowany w LOCKu poprzednika §1): s(x) 2D, Φ=s²,
ds/dτ = κ∇²s − V′(s), V=0.5as²−(b/3)|s|³+0.25cs⁴, a=0.5 b=1.6 c=1.0,
κ=0.5; N=128, L=64, dx=0.5, dt=0.02, torus, Laplasjan 5-pkt; ε=0.30,
A_min=4/128² (frakcja); s*=1.1741657, s_bar=0.4258343. ZAKAZ zmian.

**Start (geneza losowa, FROZEN):** s₀ = A·f, gdzie f = GRF pasmowy
zbudowany GRID-NIEZALEŻNIE: współczynniki zespolone na modach
|n_x|,|n_y| ≤ 16 losowane rng(seed) standard_normal, hermityzowane,
obwiednia spektralna exp(−(|n|/8)²), wbudowane w siatkę N (jak
konstrukcja pasmowa linii poziomu 1), normalizacja std(f)=1 liczona
na budowie N=256 (wspólna dla obu siatek). **A ∈ {0.6, 1.0}**
(nadbarierowe frakcje pola: |s₀|>s_bar na ~46%/67% węzłów — obie
strony testują reżim locku). **Seedy: {20260906, 20260907, 20260908,
20260909, 20260910}** (5 genez × 2 amplitudy = 10 biegów głównych).

**Przebieg:** steps_max = 100000 (τ_max = 2000); próbkowanie co 100
kroków (Δτ=2); H_Γ nierosnące (guard).

**Obserwable (FROZEN):**
- maski domen: plus = (Φ≥ε)∧(s>0), minus = (Φ≥ε)∧(s<0),
  bare = (Φ<ε); frakcje plus_area, minus_area, bare_area;
- **L_wall(τ)** = dx · liczba krawędzi (4-sąsiedztwo, periodycznie)
  z s_i·s_j < 0 (standardowy estymator długości interfejsu);
- N_dom± = liczba spójnych komponentów masek plus/minus
  (sklejanie periodyczne union-find — maszyneria linii);
- persistence_tail(X) = X(krok 100000)/X(krok 90000) (ostatnie 10%).

**Kontrola pinningu (dziedziczona reguła):** dla biegów POZYTYWNYCH
Q-NET-C (trwała sieć) re-run N=256, dx=0.25, dt=0.005, ta sama
konstrukcja GRF (grid-niezależna) i τ_max — klasyfikacja musi się
jakościowo zgadzać.

## 2. Fazy i kryteria (zalockowane)

### Phase 1 — bramka estymatorów (FAIL ⟹ STOP)
- P1a: IC „paski" poprzednika (2 proste ściany): L_wall = 128.0
  DOKŁADNIE (2·L; 256 krawędzi zmiany znaku × dx) oraz N_dom+=1,
  N_dom−=1 (sklejanie periodyczne) — osiągalny FAIL estymatorów;
- P1b: IC kropla antyfazowa R=8: L_wall(τ=0) = 8R ± 10% (estymator
  krawędziowy mierzy długość MANHATTAŃSKĄ; dla izotropowego okręgu
  znany czynnik stereologiczny 4/π: (4/π)·2πR = 8R — Amendment A1
  pre-code, patrz §6); N_dom−=1; po τ=60 (> τ_life=52 poprzednika):
  L_wall=0 — ciągłość dynamiki;
- P1c: GRF grid-niezależny: f zbudowane na N=128 vs interpolacja
  biliniowa f z N=256 do węzłów N=128: max|Δ| ≤ 0.05 oraz
  |std₁₂₈−std₂₅₆|/std₂₅₆ ≤ 0.02;
- P1d: H_Γ nierosnące w krótkim biegu genezy (seed 20260906, A=1.0,
  2000 kroków).

### Phase 2 — geneza sieci (10 biegów głównych)
- **Q-NET-A (mozaika; punkt kontrolny τ=100):** PASS-MOSAIC ⟺
  we WSZYSTKICH 5 biegach A=1.0: plus_area ≥ 0.10 ∧ minus_area ≥ 0.10
  ∧ metric_area = plus+minus ≥ 0.80 w τ=100. Inaczej: FAIL-ONESIGN
  (któryś bieg: jeden znak <10% przy metric≥0.8) / FAIL-BARE
  (metric<0.8 — lock nie zaszedł) — raport per bieg. A=0.6
  deskryptywnie (wrażliwość na amplitudę, bez wpływu na werdykt).
- **Q-NET-B (grubienie):** dla biegów A=1.0 z L_wall(τ=1000)>0:
  slope regresji ln L_wall vs ln τ w oknie **τ∈[100,1000]** (FROZEN,
  zakaz przesuwania). PASS-COARSENING ⟺ slope = −0.5 ± 0.15
  w ≥3 z 5 biegów (ważne okna). FAIL ⟺ większość poza tolerancją.
  INCONCLUSIVE ⟺ <3 ważnych okien (sieci znikają przed τ=1000 —
  raport wprost, to też informacja).
- **Q-NET-C (trwałość topologiczna; stan w τ_max=2000):** klasyfikacja
  per bieg: SINGLE-DOMAIN (min(plus,minus) < 0.02 ∧ L_wall < 0.05·L)
  / WOUND-STRIPES (min(plus,minus) ≥ 0.05 ∧ persistence_tail(L_wall)
  ∈ [0.99, 1.01] ∧ persistence_tail(bare_area) ≥ 0.99)
  / TRANSIENT (inne — sieć wciąż grubieje w τ_max).
  **PASS-TOPOLOGICAL ⟺ ≥1 bieg A=1.0 w klasie WOUND-STRIPES
  potwierdzony kontrolą N=256 (ta sama klasa).**
  FAIL-COARSEN-AWAY ⟺ wszystkie biegi A=1.0 SINGLE-DOMAIN.
  Inaczej INCONCLUSIVE (TRANSIENT-y: raport L_wall(τ_max) i trend).
- **Deskryptywnie OBOWIĄZKOWO (widok poziomu 1):** dla stanu w τ_max
  każdego biegu: bare_area (frakcja objętości przedmetrycznej ukrytej
  w „próżni" ψ), rozkład L_wall, mapa ψ=Φ/Φ* zapisana do npz;
  średni rozmiar domeny ⟨d⟩ = 2·metric_area·L²/(L_wall·…) — raport
  L²·metric/L_wall (skala mozaiki).

### Phase 3 — zamknięcie
`Phase_FINAL_close.md`, `NEEDS.md` (user-gated), `README.md`.

## 3. Forbidden moves
1. Zmiana modelu (a,b,c,κ,ε,dx,dt), estymatorów, okna fitu, progów
   klasyfikacji, seedów, amplitud po starcie obliczeń — zakaz.
2. Zakaz dobierania seedów/okien po obejrzeniu danych; wszystkie
   10 biegów raportowane.
3. Los struktur wyłącznie z pola s; zero reguł pomocniczych.
4. Rdzeń .tex/STATE/git nietykane; katalogi innych cykli tylko odczyt.
5. Zero claimów kosmologicznych/obserwacyjnych; τ ≠ czas fizyczny.
6. Wyniki negatywne raportowane wprost; INCONCLUSIVE ≠ pozytyw.

## 4. Deliverables
`Phase1_gate.py`+output, `Phase2_network.py`+output+npz stanów,
`Phase_FINAL_close.md`, `NEEDS.md`, `README.md`.

## 5. Drzewo decyzyjne
```text
P1 FAIL → STOP (estymatory/maszyneria)
A-PASS ∧ B-PASS ∧ C-PASS-TOPOLOGICAL → geneza produkuje trwałą sieć
   przedmetryczną (ukrytą dla ψ): NEEDS PILNE: wersja 3D (arkusze),
   statystyka konfiguracji nawiniętych, user-gate dopisek core
   (struktura wielkoskalowa substratu)
A-PASS ∧ C-FAIL-COARSEN-AWAY → sieć zawsze przejściowa (żyje ~τ^1/2);
   trwałość przedmetryczna wymaga topologii narzuconej, nie genezy —
   NEEDS: warunki brzegowe/rozmiar vs czas grubienia (skala!)
A-FAIL-ONESIGN → geneza łamie Z2 efektywnie globalnie (raport mechanizmu)
A-FAIL-BARE → amplitudy nadbarierowe nie lockują (napięcie z genezą
   poprzednika — raport wprost)
B poza −0.5±0.15 → prawo grubienia niestandardowe (wpływ fazy gołej
   trzeciej — deskryptywnie mechanizm)
INCONCLUSIVE → raport + ewentualny re-lock τ_max
```

## 6. Log poprawek LOCK (wszystkie PRE-CODE)

1. **Amendment A1 (PRZED napisaniem jakiegokolwiek kodu):** pierwotny
   zapis P1b oczekiwał L_wall(τ=0)=2πR±10% dla kropli — niespójny
   z zalockowanym estymatorem (liczba krawędzi zmiany znaku = długość
   MANHATTAŃSKA interfejsu; dla izotropowej krzywej czynnik
   stereologiczny 4/π ⟹ okrąg daje 8R, ściana osiowa dokładnie 2L).
   Skorygowano oczekiwanie do 8R±10%. Kryteria Q-NET-B (slope log-log)
   i Q-NET-C (ogony ilorazowe) są niewrażliwe na stały czynnik
   estymatora — bez zmian. Żadne dane nie były obejrzane.

---

**LOCK ZAMKNIĘTY. Zmiany poniżej tej linii po starcie obliczeń
= forbidden move.**
