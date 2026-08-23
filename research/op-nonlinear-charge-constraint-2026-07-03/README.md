---
title: "op-nonlinear-charge-constraint — stabilność μ/τ przy zachowanym ładunku (właściwy Q-ball: Noether + Vakhitov–Kolokolov)"
type: research_cycle
status: CLOSED-EXECUTED
phase: 3
folder_status: closed
claim_status: "CLOSED-EXECUTED 2026-07-04 (sesja #63, osobny agent per autoryzacja #62). WERDYKTY: V1 NEGATYWNY (M0 nie ma budżetopodobnego ładunku zachowanego — pełna tabela C1–C5, sympy exact; per LOCK dalej wyłącznie gałąź M1 jako model-extension); V2 NEGATYWNY dla μ i τ (VK slope-positive na całych gałęziach ω≤0,25, deflacja ładunkowa nie usuwa modów głębokich, kontinuum tachioniczne przy każdym ω, próżnia ω duchowiona od ω_gh=0,2935 — klauzula tła z LOCK-a aktywna); V3 kierunek (i): niestabilność μ potwierdzona NIELINIOWO w M0-f_ε (wzrost wykładniczy zgodny z teorią liniową, zero saturacji, wyjście pola z dziedziny modelu g→0 w t*≈1,7–3,6). Hipoteza budżetowa autora obalona także w wersji nieliniowej/ładunkowej — w całej przetestowanej klasie."
created_date: 2026-07-03
closed_date: 2026-07-04
authorization: "User 2026-07-03 (sesja #62): 'ok zrób pożądki a potem rozpisz cykl badawczy N4 dla nowego agenta'; realizacja: sesja #63"
anti_lakatos_lock: PRESERVED
related:
  - "[[Phase0_balance.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-wall-dynamics-2026-07-03/README.md]]"
  - "[[../op-spectral-analysis-Phi-2026-07-03/README.md]]"
  - "[[../../audyt/L03_K_phi_stability/POST_ACTION_UPDATE_2026-07-04.md]]"
---

# op-nonlinear-charge-constraint — CLOSED-EXECUTED (2026-07-04)

**Pytanie cyklu:** czy istnieje ZACHOWANY ładunek (inwentarz Noether
zanurzenia M0 lub minimalne rozszerzenie U(1) — M1), przy którego
ustaleniu solitony μ/τ korony są orbitalnie stabilne (VK/GSS) — tj. czy
hipoteza budżetowa autora jest prawdziwa na poziomie nieliniowym, po
obaleniu wersji liniowej w #62?

**Odpowiedź: NIE — w całej klasie zamkniętej LOCK-iem.** Wyniki
negatywne zgłoszone wprost per [[Phase0_balance.md]] (kryteria V1–V3
niezmienione; jedyna dopuszczona korekta — konwencja znaku VK +
renormalizacja pudła — udokumentowana w Phase 1 PRZED P2b, zgodnie z
LOCK §8).

## V1 — inwentarz ładunków w M0 (Phase 1, sympy exact) → NEGATYWNY

[[Phase1_charge_inventory.py]] · [[Phase1_output.txt]]

- **P1a PASS:** EOM M0 wyprowadzone; energia zachowana exact
  (gęstość+strumień); statyczne EL(M0) = ODE korony a3d/CP-7 (exact).
- **P1b (kryterium V1):** pełna tabela — **0/9 kandydatów zachowanych**
  (C1: ∫(g−1)r²; C2: pęd kanoniczny ∫fġr²; C3: ∫(g−1)²r²; C4: h∈{1,g²};
  C5: pęd radialny + dylatacje c∈{0,1,3/2}; metoda: exactness
  operatorem Eulera na jetach — dywergencja zupełna on-shell ⟺
  zachowanie; residua potwierdzone sondami numerycznymi). Energia:
  zachowana (kontrola dodatnia; nie jest ładunkiem budżetowym).
  **Per LOCK: „hipoteza wymaga rozszerzenia M1"** — Phase 2 wyłącznie
  w gałęzi M1, jawnie model-extension (nie wchodzi do core).
- **P1c (M1):** ładunek Noether Q zachowany EXACT (θ-EL ≡ prawo
  lokalne dQ/dt = strumień brzegowy); redukcja ψ=φe^{−iωt} →
  W_eff = W − (ω²/2)fφ² (przy ω→0 dokładnie ODE korony — ciągłość
  gałęzi); rozkład GSS: L₊ (forma CP-7 z W→W_eff), L₋ z **L₋φ_ω = 0
  exact** (mod fazowy). **Konwencja VK zalockowana przed P2b:**
  Q := ω∫fφ²r²dr > 0; gałąź stabilna ⟺ dQ_sol/dω < 0; wielkości
  odjęte od ω-próżni (renormalizacja pudła R=60).
- **P1d (pytanie zalockowane):** człon ω² **NIE podnosi** krawędzi —
  **OBNIŻA ją**: c(ω) = −1 − 7ω² − 117ω⁴ + O(ω⁶) (seria exact z
  przesunięciem próżni φ_∞(ω) = 1 − 3ω² + 3ω⁴ + …); krawędź L₋ = 0
  exact dla każdego ω. Próżnia przecina ścianę duchową g* = e^{−1/4}
  przy **ω_gh = 0,2935** (dla ω > ω_gh tło kinetycznie zduchowione).
  **Nie istnieje ω_min z σ_ess ≥ 0.** Klauzula tła z V2 aktywna.

## V2 — rodzina Q-ball + VK w M1 (Phase 2) → NEGATYWNY dla μ i τ

[[Phase2_qball_family.py]] · [[Phase2_output.txt]] ·
[[Phase2b_Rcontrol.py]] · [[Phase2b_Rcontrol_output.txt]]

- **P2a:** gałęzie φ_ω (hard-wall, DOP853 rtol 1e−10; skan ω∈[0,1] krok
  0,05 W CAŁOŚCI): istnieją i są ciągłe z profilami CP-7 dla
  **ω ≤ 0,25** (e/μ/τ); dla ω ≥ 0,30 profil kolabuje na ścianę
  (29 odbić, nie wypełnia pudła) — zbieżnie z duchowieniem tła
  (ω_gh = 0,2935).
- **P2b (VK):** **dQ_sol/dω > 0 na całych gałęziach μ i τ**
  (slope-positive = niestabilna w konwencji zalockowanej w P1c);
  jedyny punkt slope-negative w skanie: e przy ω = 0,05 (elektron i
  tak stabilny bez ładunku — nie jest przedmiotem hipotezy).
- **P2c:** spektra L₊/L₋: deflacja kierunku ładunkowego
  (dQ/dφ ~ f′φ²+2fφ) i rodziny (dφ/dω) **nie usuwa modów głębokich**
  (N_c = N_loc wszędzie); μ: N_loc = 2,2,2 dla ω = 0–0,10 (przy
  ω = 0,10 mod pogłębia się do −3,35!), N_loc = 0 dla ω = 0,15–0,25
  wyłącznie dlatego, że krawędź nurkuje (−1,24…−3,12) pod mody —
  kontinuum pozostaje tachioniczne, to nie jest stabilizacja; τ:
  N_loc = 3→1–2, nigdy 0, mod −4,2 zawsze obecny. Zbieżność
  N∈{2000,4000,8000} zgodna; R-kontrola {40,80} dla ω=0,05–0,25:
  patrz [[Phase2b_Rcontrol_output.txt]]. L₋: λ_min ≈ 0 (mod fazowy;
  drobne odchylenia dyskretyzacyjne raportowane).
- **P2d (masy):** gate ω→0: r₂₁ drift 0,0005%, r₃₁ drift 0,0012%
  (< 0,1% **PASS** — baseline #62 odtworzony); przy ω > 0 dryf
  raportowany w całości bez progu (per LOCK): 5–25% (ω=0,05–0,10) do
  50–100% (ω=0,15–0,25) — „podkręcenie" Q-ballowe niszczy dopasowanie
  mas.
- **Cross-check f_ε (wyłącznie ε=0,2 + kontrola 0,1, per LOCK):**
  μ kolabuje przy ω ≥ 0,15 (ε=0,2) i ω ≥ 0,10 (ε=0,1); τ kolabuje
  zawsze (jak #62); e przetrwa do ω=0,30 — spójne z obrazem hard-wall.
- **Werdykt V2 (koniunkcja i–iii):** dla μ i τ **żadne ω** nie spełnia
  koniunkcji (VK fail wszędzie; (iii) fail lub bezprzedmiotowe).
  Dodatkowo klauzula tła z LOCK-a: kontinuum tachioniczne przy każdym
  ω ⇒ w zakresie, w którym kryterium w ogóle jest rozstrzygalne,
  wynik **NEGATYWNY zgłoszony wprost**; poza nim — nierozstrzygalne
  z powodu tła (też wynik, per LOCK).

## V3 — nieliniowy test dynamiczny w M0-f_ε (Phase 3) → kierunek (i): NIESTABILNOŚĆ POTWIERDZONA NIELINIOWO

[[Phase3_nonlinear_evolution.py]] · [[Phase3_output.txt]]

- Metoda: dokładny układ hamiltonowski semi-dyskretny (energia
  dyskretna zachowana exact w ODE; gate mierzy wyłącznie błąd RK4);
  g(r,0) = g_eq + a·v_deep, ġ(r,0) = 0, a ∈ {±0,01; ±0,03};
  **gate |ΔE|/E ≤ 2,4e−8 PASS** we wszystkich runach; zbieżność dt
  (0,004 vs 0,002): wyniki identyczne (t*, σ co do 4 cyfr).
- **Wynik (ε=0,2, główny):** wzrost wykładniczy a(t) = ⟨v_deep, δg⟩ z
  σ_fit = 0,97–1,74 vs √1,389 = 1,18 (3/4 runów w ±20%; odchylenie
  +47,8% w a=+0,01 = kontaminacja okna przez szybszy kierunek
  F-ważony, patrz niżej); **zero saturacji** (‖δg‖_∞ rośnie
  monotonicznie do 80–136% tła); **pole opuszcza dziedzinę modelu
  (g → 0, log niezdefiniowany) w skończonym czasie t* = 3,62** przy
  każdej amplitudzie i każdym znaku. Kontrola ε=0,1: to samo, szybciej
  (t* = 1,7–3,3). Rekurencji brak.
- Udokumentowana subtelność normalizacyjna (zapisana przed
  interpretacją, nie zmiana kryterium): λ_min = −1,3896 potwierdza
  #62 (waga 1), ale DOKŁADNA dynamika liniowa rządzi się problemem
  F-ważonym: λ_F = −7,86 (ε=0,2) / −52,4 (ε=0,1) ⇒ σ_F = 2,80/7,24 —
  miękka ściana (f_ε → ε²/4 przy g*) czyni region ścienny skrajnie
  szybkim. Zmierzone σ mieszają oba kierunki; kierunek werdyktu
  niezależny od normalizacji.
- **Werdykt V3:** kryterium (i) spełnione co do struktury (wzrost
  wykładniczy + brak saturacji poniżej 10% tła), z zastrzeżeniem
  normalizacyjnym jw.: **„niestabilność potwierdzona nieliniowo w
  M0-f_ε"** — i silniej: nieliniowość nie stabilizuje, lecz prowadzi
  do ucieczki pola z dziedziny modelu (dynamiczny odpowiednik
  statycznego kolapsu τ z #62 W2b). τ: poza zakresem P3 (brak
  reprezentanta EL w f_ε — kolaps, #62) — odnotowane.

## Synteza (dla STATE/#63)

Po CP-7 (#60: siodła + tachioniczne kontinuum), op-wall-dynamics
(#62: deflacja liniowa niewystarczająca; brak gładkiego zamiennika
ściany) i tym cyklu (#63): **wszystkie trzy ścieżki stabilizacji μ/τ
wewnątrz gładkiej teorii pola z odbiciem ad-hoc są zamknięte
negatywnie.** Hipoteza budżetowa autora w wersjach: liniowej (#62),
ładunkowej M0/M1 i nieliniowo-dynamicznej (#63) — obalona/negatywna.
Pozostałe kierunki (user-gated, [[NEEDS.md]] N4): dyskretność
substratu, inna symetria nośna budżetu, sektor sprzężony F-A (L04),
reinterpretacja metastabilnościowa.

## Limitations

1. Ładunki Q/E w pudle R=60 wymagają renormalizacji od ω-próżni
   (konwencja zalockowana w P1c); ogony odbite wnoszą człony brzegowe
   (diagnostyka dE = ωdQ: ratio ~1,6 zamiast 1 — raportowane,
   nie-kryterium).
2. Odbicie hard-wall pozostaje regularyzacją nie-EL (jak CP-7/#62);
   wyniki M1 dotyczą tej konwencji ściany (zamkniętej w LOCK §2).
3. σ_fit w P3 mierzy projekcję na v_deep (waga 1) w dynamice rządzonej
   operatorem F-ważonym — porównanie ±20% z √1,389 obarczone
   mieszaniem kierunków (udokumentowane w outputach).
4. M1 = rozszerzenie ontologii (U(1) na fazie zanurzenia) — użyte
   wyłącznie jako model-extension; NIE weszło do core.

## Deliverables

- [[Phase0_balance.md]] (LOCK, nietknięty) · [[Phase1_charge_inventory.py]] + [[Phase1_output.txt]]
- [[Phase2_qball_family.py]] + [[Phase2_output.txt]] · [[Phase2b_Rcontrol.py]] + [[Phase2b_Rcontrol_output.txt]]
- [[Phase3_nonlinear_evolution.py]] + [[Phase3_output.txt]]
- [[NEEDS.md]] (N1–N4, user-gated; zero edycji core wykonanych)
- [[../../audyt/L03_K_phi_stability/POST_ACTION_UPDATE_2026-07-04.md]] · wpis STATE.md #63

## Anti-Lakatos

✓ Kryteria/listy/siatki/konwencja ściany bez zmian po starcie obliczeń.
✓ Jedyna dopuszczona korekta (LOCK §8: konwencja znaku VK +
renormalizacja pudła) udokumentowana w Phase 1 PRZED P2b. ✓ Skany
raportowane w całości (w tym punkty INCOMPLETE/GHOSTED). ✓ Zero
re-fitowania g₀ i stałych korony (gate mas ω→0 PASS: 0,0005%).
✓ Wyniki negatywne zgłoszone wprost (V1, V2, V3). ✓ Niezbieżności i
kolapsy raportowane jako takie. ✓ M1 nie wszedł do core (NEEDS
user-gated). ✓ Rdzeń .tex NIETKNIĘTY. ✓ NIE commitowano.
