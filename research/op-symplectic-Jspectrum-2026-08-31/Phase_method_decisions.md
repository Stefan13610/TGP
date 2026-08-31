---
title: "Phase_method_decisions — decyzje implementacyjne ZAMROŻONE przed startem obliczeń (forma L₊/L₋ i czynnik ½ z E, dyskretyzacja on-shell L₋, odczyt tol/max|Im λ|, siatka k, konstrukcja bramek C1/C2/C3, diagnostyka modów symetrii)"
date: 2026-08-31
type: method-decisions
tgp_owner: research/op-symplectic-Jspectrum-2026-08-31
status: FROZEN-PRE-COMPUTE
computations_performed: ZERO
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/Phase_method_decisions.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/Phase3_bloch_spectrum.py]]"
---

# Decyzje metodyczne (zamrożone PRZED jakimkolwiek obliczeniem)

Wszystkie decyzje zapisane przed uruchomieniem pierwszego skryptu cyklu.
Kryteria/progi/tła LOCKa (Phase0_balance.md) NIETKNIĘTE — poniżej wyłącznie
doprecyzowania implementacyjne, które LOCK delegował lub pozostawił
niedomknięte (każde z osiągalnym FAIL i z uzasadnieniem).

## 1. Forma operatorów i czynnik ½ (rachunek ręczny; weryfikacja = P1c sympy)

Model LOCK §2 dosłownie: E[u] = ∫[½K(|u|)|∂_x u|² + U(|u|)]dx, K=g⁴,
U=g⁷/7−g⁸/8, i∂_t u = δE/δu*, tło rzeczywiste u₀=g_d, u = g_d + a + ib.

Rachunek ręczny (do weryfikacji sympy w P1c; FAIL sympy ⟹ STOP wg LOCKa):

```
E₂ = ½⟨a, L₊ a⟩ + ½⟨b, L₋ b⟩   (brak członu mieszanego: g rzeczywiste)
L₊ = −∂_x(K ∂_x) + Q₊,  Q₊ = U″(g) − K′(g)g″ − ½K″(g)g′²
L₋ = −∂_x(K ∂_x) + Q₋,  Q₋ = K′(g)g′²/(2g) + U′(g)/g
tożsamości on-shell (EL: (Kg′)′ = ½K′g′² + U′):
  Q₊ = K·(2βg − 3γg² + 2g′²/g²)   [= dokładnie Q z Phase 3 cyklu bloch]
  Q₋ = (Kg′)′/g                    [⟹ L₋ g = 0 exact — mod fazowy]
```

**Czynnik ½ w dynamice:** z konwencji Wirtingera δE/δu* = ½(δE/δa + iδE/δb)
i z normalizacji LOCKa ½K|u_x|² wynika ∂_t(a,b) = (½L₋ b, −½L₊ a), czyli
JL̂ = [[0, ½L₋], [−½L₊, 0]] i λ² ∈ spec(−¼ L₋L₊). Schematyczny cross-check
LOCKa „λ² = −ν, ν = spec(L₋L₊)" implementujemy DOKŁADNIE jako λ² = −ν/4
(czynnik ¼ jest konsekwencją normalizacji E z LOCKa §2, nie zmianą kryterium;
kryterium max Re λ ≤ tol jest niezależne od dodatniego przeskalowania czasu).
Wyprowadzenie czynnika = część P1c (sympy).

Bramki C1/C2 przechodzą przez TĘ SAMĄ maszynerię: NLS kubiczny/nadkrytyczny
zapisujemy w rodzinie kanonicznej (K, U) — K ≡ 2 (stała), U_C1 = g² − g⁴/2,
U_C2 = g² − g⁸/4 (ramka rotująca ω=1 wchłonięta do U; wtedy ½L₊ = −∂² +1−3φ²
itd. — formy podręcznikowe). Weryfikacja tego wchłonięcia = P1 (sympy).

## 2. Dyskretyzacja (FD2, zgodna z cyklem bloch — spójność z tłami i kotwicą)

- Pełna komórka, węzły x_j = j·h, h = d/N; N ∈ {400, 800} — te same siatki
  co tła w `Phase2_backgrounds.npz` (tła używane WPROST, bez interpolacji).
- Stencil −∂(K∂) verbatim z bloch: K_{j+½} = ½(K_j + K_{j+1}); Bloch-faza
  e^{ikd} na sprzęgle narożnym. L₊(k), L₋(k) Hermitowskie dla każdego k.
- **Q₊ dyskretne:** forma po podstawieniu EL, VERBATIM bloch:
  Q₊ = K·(2βg − 3γg² + 2·D1(g)²/g²), D1 centralne periodyczne. Powód:
  gate kotwicy λ_min(3π) = −1.222191 ± 1e−4 wymaga IDENTYCZNEGO operatora.
- **Q₋ dyskretne (decyzja):** forma on-shell Q₋ := [D(K·D g)]_j / g_j
  z TYM SAMYM stencilem półwęzłowym co część kinetyczna L₋. Skutek:
  dyskretny mod fazowy L₋·g = 0 DOKŁADNIE w arytmetyce fp (tożsamość
  analityczna on-shell z §1; kontrola: zgodność z formą off-shell
  K′g′²/(2g)+U′/g raportowana — oczekiwane O(h²); rozjazd O(1) = czerwona
  flaga implementacji).
- **JL̂:** macierz 2N×2N zespolona gęsta M = [[0, ½L₋],[−½L₊, 0]];
  solver `scipy.linalg.eig` (pełny, gęsty), sortowanie po Re λ (LOCK).
  Symetria czwórek ±λ, ±λ* raportowana jako sanity.
- **Cross-check produktowy:** eig(L₋·L₊) → ν; porównanie
  (max Re λ)² vs max(−Re ν)/4 w k=0 dla każdego tła (i w bramkach).

## 3. Siatka k i wartość główna

- LOCK: „≥8 punktów w [0, π/d]". Zamrażam: **9 punktów równomiernie
  z KRAŃCAMI** (linspace(0, π/d, 9)) — zawiera k=0 (PRIMARY) i krawędź.
- Wartość główna max Re λ: siatka N=800 (2N); zbieżność wg LOCKa:
  |Δ max Re λ| ≤ 0.01·max(max Re λ, 0.01), Δ między N=400 a N=800,
  per tło i per k (raport per tło; zbieżność wymagana we wszystkich k
  dla werdyktu tła).

## 4. Odczyt tol (ruling zapisany PRZED obliczeniami; oba odczyty raportowane)

LOCK: tol := max(1e−6, 1e−3·max|Im λ|) „na siatce bazowej" — bez
specyfikacji zbioru, po którym biegnie max|Im λ|. Dwa odczyty:

- **Literalny-nieograniczony:** max po CAŁYM policzonym σ(JL̂). Wada
  strukturalna: gałąź UV skaluje się jak h⁻² (|Im λ| ~ ½κ_max²K),
  więc tol ∝ N² — kryterium grid-rozbieżne i przy N=400 (tol ~ O(10))
  klasyfikowałoby jako „stabilne" każde tempo wzrostu O(1); to czyni
  FAIL bramki C1 praktycznie nieosiągalnym i przeczy LOCK §4 p.2
  („każdy test z osiągalnym FAIL") oraz wymogowi zbieżności.
- **PRIMARY (zamrożony): pasmowo-ograniczony:** max|Im λ| po
  σ(JL̂) ∩ {|λ| ≤ Λ_band}, **Λ_band := 12** (skala: analityczna gałąź
  próżni |λ| = ½κ√(κ²−γ) osiąga ≈12.2 przy κ=5, co obejmuje wszystkie
  fizycznie istotne pasma dla każdego d oraz bramek NLS; wybór zapisany
  PRZED obliczeniami, nie strojony). Fallback (pusty przekrój —
  nieoczekiwane): tol = 1e−6.

Oba tol (PRIMARY i literalny) liczone i raportowane w KAŻDYM punkcie
werdyktu; jeżeli werdykt zależy od odczytu — oba werdykty raportowane
jawnie (wzorzec Phase3_verdict_ruling cyklu bloch), bez zmiany progów.

## 5. Mody symetrii (diagnostyka deskryptywna, nie-kryterialna)

Dyskretyzacja rozszczepia bloki Jordana λ=0 (translacja: (g′,0);
faza: (0,g)) na pary λ ~ ±√(O(h²)) — możliwe pozorne Re λ ~ O(h).
Reguła pre-rejestrowana (lustro reguły Goldstone'a z bloch §2):
w k=0 identyfikujemy mody symetrii overlapem ≥ 0.9 składowej a z g′
lub składowej b z g (normy L²); ich λ raportujemy OSOBNO. Kryterium
LOCKa stosujemy LITERALNIE do pełnego max Re λ; deskryptywnie (bez
wpływu na werdykt) raportujemy też max Re λ z wyłączeniem modów
symetrii. Ochrona przed fałszywym Q-FAIL: pozorne λ ~ O(h) nie
przechodzi testu zbieżności (halving z N), więc scenariusz
prawdziwie-stabilny daje co najwyżej Q-INCONCLUSIVE, nie Q-FAIL.

## 6. Konstrukcja bramek Phase 2

- **Gate kotwicy L₊:** problem uogólniony ważony VERBATIM bloch
  (L₊φ = λ K φ po symetryzacji B^{−½}L₊B^{−½}, B = diag(K)), tło
  3π/0.7 z npz, N=800, k=0: λ_min musi być −1.222191 ± 1e−4.
  Macierz L₊ budowana przez NOWY kod tego cyklu (to jest gate
  reprodukcji operatora, nie kopiowanie liczby).
- **C1/C2:** pudełko periodyczne L_box = 40, siatki N ∈ {1024, 2048}
  (bazowa 1024), soliton wycentrowany; seedy analityczne
  φ_C1 = √2·sech(x−20), φ_C2 = 4^{1/6}·sech^{1/3}(3(x−20)); polish
  Newtona dyskretnego równania statyki D2φ − φ + φ^p = 0 (p=3, 7)
  na pół-komórce w klasie parzystej (pinning translacji; wzorzec
  relaksacji bloch §3), residuum cel ≤ 1e−12, lustrzane rozwinięcie.
  Weryfikacja seedu: ‖residuum analityczne‖∞ raportowane (test formy
  zamkniętej). k = 0 (pudełko izolowane, nie łańcuch).
  C1: max Re λ ≤ tol (odczyt §4; oba raportowane).
  C2: musi istnieć λ z Re λ > 0 (LOCK, literalnie); deskryptywnie
  wymagamy też separacji max Re λ(C2) ≥ 10× max Re λ(C1).
- **C3:** próżnia g≡1 na komórce d = 2π, N ∈ {400, 800}, 9 punktów k;
  porównanie 8 najniższych gałęzi λ² z analitycznym
  λ²(κ) = −¼κ²(κ²−γ), κ = k + 2πn/d (wyprowadzenie: Phase 1).
  **Metryka (frozen):** |λ²_num − λ²_ex| / max(|λ²_ex|, 1) ≤ 1e−3
  (podłoga 1 konieczna — gałąź przechodzi przez λ²=0; lustro metryki
  P1a cyklu bloch). FAIL osiągalny: zły czynnik ½/znak/faza ⟹ O(1).

## 7. Phase 3 — protokół

4 tła z npz: 3π/0.7, 4π/0.7, 6π/0.7 (2-garb), 6π/S3single (1-garb) —
klucze wprost z pliku; ZERO relaksacji. Dla każdego tła: N ∈ {400,800}
× 9 punktów k: pełne σ(JL̂), max Re λ, tol (oba odczyty, siatka bazowa
N=400), zbieżność per k; w k=0 dodatkowo: diagnostyka modów symetrii
(§5), cross-check produktowy (§2), tożsamość dyskretna ‖L₋g‖∞ (musi
być ~fp-0) i ‖L₊g′‖∞ (oczekiwane O(h²), raportowane). Werdykt per tło
i zbiorczy wg LOCKa §3. JSON: Phase3_results.json. Flaga [INPUT-ONTO]
(u₀ = g_d; złożenie u↔g z LOCKa §2) w każdym bloku wyników.

## 8. Phase 4 (tylko Q-PASS) — doprecyzowania

Superkomórka 4d (N_sc = 1600), perturbacja 0.01·‖tło‖∞ wzdłuż dawnego
modu runaway (mod amplitudowy k=0 z bloch Phase 3; re-diagonalizacja
L₊ ważonego w superkomórce), t_end = 3·3.62 = 10.86. Integrator:
split-step Stranga — krok liniowy exact w Fourierze dla −½∂² (K(1)=1)
+ reszta (−½∂((K−1)∂u) + (u/2ρ)(½K′|u_x|²+U′)) krokiem RK4; dt = 0.002,
kontrola dt/2. Gate LOCKa: |ΔQ|/Q i |ΔE|/E ≤ 1e−8. Brak ucieczki:
‖|u|−g_d‖∞ ≤ 3× amplituda początkowa perturbacji przez cały bieg.

## 9. Higiena wykonania

- `python` = CPython 3.14.2; numpy 2.4.3, scipy 1.17.1, sympy 1.14.0 —
  SPRAWDZONE importem przed startem (jedyne wykonane polecenie: wersje
  + odczyt kluczy npz; zero obliczeń merytorycznych).
- Zero `cd`; pełne ścieżki; po każdym zapisie weryfikacja `ls`.
  Outputy: `Phase*_output.txt` w katalogu cyklu (redirect stdout).
- Rdzeń `.tex`, STATE.md, git, katalogi innych cykli — NIETYKANE.

**FROZEN 2026-08-31, przed uruchomieniem jakiegokolwiek skryptu cyklu.**
