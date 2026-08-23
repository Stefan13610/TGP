# Wynik: rozpraszanie e–p w arenie relatywistycznego baby-Skyrme (falsyfikator §4)

**Data:** 2026-07-28
**Typ:** dokument-wynik (falsyfikator §4 z RAMY), weryfikacja u źródła — kod uruchomiony, liczby z faktycznych przebiegów
**Werdykt jednozdaniowy:** W relatywistycznej arenie baby-Skyrme operacyjny „ładunek" (przekaz pędu) **NIE jest Coulombowski** — znak siły ustala względna orientacja wewnętrzna χ, nie iloczyn ładunków topologicznych; oddziaływanie jest ekranowane wykładniczo (Yukawa), nie długozasięgowe. **Negatyw typu F-A: mocny, bo zamyka trop.**

---

## 0. Co testowaliśmy i dlaczego to nie jest toy

Falsyfikator §4 RAMY: *czy selektywność i znak przekazu pędu (rozpraszanie Coulomba)
WYPADAJĄ z dynamiki kształtów, bez postulowania EM.* Wybrana arena (bramka projektowa,
zatwierdzona przez użytkownika): **relatywistyczny baby-Skyrme 2+1D** — najmniejsza,
w której teza charge≠color da się testować, a sprzeczność spin-vs-Coulomb staje się liczbą.

**Reguła nie-cyrkularności (wymuszona konstrukcją):** do lagranżianu NIE wpisano żadnego
prawa parami q_iq_j, pola cechowania ani propagatora Coulomba. Jedyne składniki:
energia O(3)+Skyrme+potencjał na polu n∈S². „Ładunek" mierzony **operacyjnie** jako
przekaz pędu / znak deflekcji. Jeśli Coulomb miał się pojawić, musiał WYPAŚĆ.

Model (mostly-minus):
```
E = INT [ ½|∇n|² + ½κ²|∂₁n×∂₂n|² + μ²(1−n₃) ] d²x,   κ²=1, μ²=0.1
Q = (1/4π) INT n·(∂₁n×∂₂n) d²x ∈ ℤ = π₂(S²)
Dynamika: ∂ₜn=π, ∂ₜπ = −P_n[δU/δn] − (π·π)n   (hamiltonowska, bezstratna, 2. rz. w czasie)
```
„e" = antyskyrmion Q=−1, „p" = skyrmion Q=+1.

---

## 1. Walidacja areny (Etapy 1–2) — DERIVED (solidne narzędzie)

**Etap 1 — statyczny soliton (`stage1_static.py`):**
- Gradient-check analitycznego δE (nie zaufano wzorowi): błąd rel. **2.98e-10**.
- Virial Derricka E_skyrme = E_pot: **1.0001** (dokładny; E=19.66).
- Zbieżność siatki: |Q| = 0.991 → 0.996 → 0.998 (N=128→192→256), E_2D = 19.51 → 19.63 → 19.69 → E_1D.
- Sprzężenie ładunku: E(Q=+1) = E(Q=−1) do bitu; Q(+)=−Q(−).

**Etap 2 — dynamika bezstratna (`stage2_boost.py`):** boost Lorentza v=0.3.
- Dryf energii: **max 4.2e-4** (bezstratność ZMIERZONA, nie założona).
- Q zachowane (1.2e-4); prędkość COM = 0.298 ≈ v; pęd pola Px stały do 5 cyfr.
- Soliton leci sztywno, bez rozpadu/promieniowania.

> Arena i integrator są wiarygodne. To certyfikuje NARZĘDZIE, nie fizykę.

---

## 2. Wynik główny (Etap 3) — trzy niezależne testy, ten sam werdykt

### 2.1 Rozpraszanie dynamiczne, silne nakładanie (`stage3_scatter.py`, b=3, v=0.4)
`dvAy` projektyla: <0 = odchylenie ku partnerowi (przyciąganie), >0 (odpychanie).

| | χ=0 | χ=π/2 | χ=π | χ=3π/2 |
|---|---|---|---|---|
| UNLIKE (−,+) | −0.133 | −0.083 | **+0.279** | +0.081 |
| LIKE (−,−)   | +0.010 | −0.143 | −0.107 | −0.157 |

→ Znak **przeskakuje z χ**; LIKE głównie przyciąga; UNLIKE nie jest jednoznacznie
przyciągające. To nie jest struktura Coulomba.

### 2.2 Rozpraszanie dynamiczne, duże b (`stage3b_largeb.py`) — decydujący dyskryminator
Coulomb jest DŁUGOZASIĘGOWY → przy słabym nakładaniu by zdominował i był znakowo-określony.

| | χ=0 | χ=π/2 | χ=π | χ=3π/2 | znaki |
|---|---|---|---|---|---|
| b=5 UNLIKE | −0.040 | −0.105 | −0.202 | −0.226 | −−−− |
| b=5 LIKE   | +0.023 | −0.176 | −0.229 | −0.191 | +−−− |
| b=8 UNLIKE | −0.008 | +0.024 | −0.088 | −0.124 | 0+−− |
| b=8 LIKE   | +0.023 | −0.016 | −0.052 | −0.019 | +−−− |

→ **LIKE przyciąga dla większości χ przy każdym b** (Coulomb wymaga: LIKE zawsze odpycha
→ sfalsyfikowane). „UNLIKE wszystko−" przy b=5 jest pozorne: przy czystszym b=8 UNLIKE
sam przeskakuje znak, a LIKE przyciąga przy tych samych χ.

### 2.3 Statyczna energia oddziaływania (`stage3c_static.py`) — bez artefaktów pędu
E_int(d) = E_pary − 2·E_pojedynczy.

| UNLIKE | χ=0 | χ=π/2 | χ=π | χ=3π/2 | | LIKE | χ=0 | χ=π/2 | χ=π | χ=3π/2 |
|---|---|---|---|---|---|---|---|---|---|---|
| d=6 | +3.79 | −0.92 | −2.10 | −0.92 | | d=6 | +2.65 | +0.29 | −0.80 | +0.29 |
| d=8 | +1.18 | −0.20 | −0.98 | −0.20 | | d=8 | +0.66 | +0.01 | −0.48 | +0.01 |
| d=12| +0.21 | −0.01 | −0.20 | −0.01 | | d=12| +0.13 | −0.003| −0.13 | −0.003|

→ Znak E_int przeskakuje z χ. **Przy χ=π OBA (unlike i like) przyciągają** — Coulomb
tego zabrania. Iloczyn Q nie ustala znaku.

### 2.4 Zanik = Yukawa, nie Coulomb (sprawdzenie predykcji, nie fit)
E_int spada ~×10 na podwojeniu d (6→12). Potencjał nadaje polu masę:
lin. wokół próżni n₃=1 → m² = μ² = 0.1 → **m = 0.316**. 2D Yukawa `~K₀(mr)~e^{−mr}/√r`:
`e^{0.316·6}·√2 ≈ 9.5` — zgadza się z zaobserwowanym ×10. **Pole masywne → oddziaływanie
ekranowane, krótkozasięgowe.** Coulomb (potęgowy/log, długozasięgowy) NIE występuje.

---

## 3. Klasyfikacja (uczciwa)

| Element | Status |
|---|---|
| Solver baby-Skyrme, integrator, boost, zachowania E/p/Q | **DERIVED** (zwalidowane, gradient-check + zbieżność + zachowania) |
| „Ładunek operacyjny NIE jest Coulombowski w tej arenie" | **DERIVED (negatyw)** — 3 niezależne testy + zgodność z zanikiem Yukawy |
| Znak siły ustala orientacja χ (dipol), nie Q | **DERIVED** |
| „Emergentny Coulomb wymaga dodatkowej struktury (pole cechowania / Berry)" | **KONCEPT** (wniosek, nie policzone tu wprost) |

Nic tu nie jest FITTED ani CYRKULARNE: nie wpisano prawa ładunku; testowano, czy wypada.
Nie wypadło.

---

## 4. Dlaczego (fizyka, teraz empirycznie ugruntowana)

Pojedynczy baby-skyrmion **nie ma źródła monopolowego** dla pola długozasięgowego, a w tej
arenie **nie ma dynamicznego pola cechowania**, które prąd topologiczny mógłby wzbudzić.
Potencjał (konieczny do lokalizacji solitonu w 2D) czyni pole **masywnym** → wiodące
oddziaływanie to **ekranowany dipol**, którego orientacja = wewnętrzna faza χ. Stąd:
znak zależny od χ, zanik Yukawy, brak selektywności Q. To znane zachowanie baby-skyrmionów
(oddziaływanie dipolowe zależne od orientacji), tu potwierdzone od zera i bezcyrkularnie.

---

## 5. Zakres — czego to NIE przekreśla i NIE dowodzi

- **NIE testowano trasy Landau–Lifszyc / Berry** (dynamika 1. rzędu z tłem), gdzie emergentne
  EM istnieje z definicji. Na bramce projektowej wybrano wariant *relatywistyczny* (dom spin-½).
  Ten wynik mówi: relatywistyczna arena daje spin, ale **nie** Coulomb — nie mówi nic
  o wariancie LL.
- **NIE testowano jawnego emergentnego pola cechowania** sprzężonego z prądem topologicznym
  (musiałoby samo emergować, nie być wstawione) — to otwarta furtka.
- **„0" (neutralne) transparentne:** nie testowano osobno — w tym modelu **nie istnieje
  stabilny soliton Q=0** (para skyrmion–antyskyrmion anihiluje). To samo w sobie ustalenie.
- Init przez surową superpozycję + silne nakładanie przy małym b dodają szum; ale testy
  duże-b i statyczny są czyste i zgodne → werdykt odporny.
- Wynik dotyczy baby-Skyrme (2D, π₂) z potencjałem μ²(1−n₃). Inny potencjał/masa=0 nie
  ratuje długiego zasięgu (bezmasowe lumps O(3) nie są zlokalizowanymi solitonami).

---

## 6. Znaczenie dla programu (oś: charge↔color z jednej topologii)

Wynik **zaostrza sprzeczność zdiagnozowaną na bramce projektowej**: dwa filary żyją w różnych
ramach dynamicznych i **nie przychodzą za darmo z jednej relatywistycznej areny Skyrme**:
- spin-½ ← relatywistyczna topologia nawinięcia (jest),
- Coulomb/EM ← wymaga dynamicznego emergentnego pola (nie ma go w czystym σ-polu).

**Konsekwencja dla decyzji B (ruchy substratu Φ):** aby EM emergowało z Coulombowską
selektywnością, sama energia kształtu σ+Skyrme NIE wystarcza. Trzeba albo (i) dynamiki
LL/Berry (tło magnetyzacji — napięcie z background-free), albo (ii) **emergentnego pola
cechowania wzbudzanego przez prąd topologiczny**, które samo musi wypaść z reguł substratu,
nie być wstawione. To jest konkretne, sfalsyfikowane pytanie do B — nie ogólnik.

---

## 7. Pliki (ta sesja, folder `op-ep-scattering-babyskyrmion-2026-07-28`)
- `bs2d.py` — moduł: model, energia, wariacja (gradient-checked), Q, dynamika RK4, boost.
- `stage1_static.py` → `soliton_profile.npz` (zwalidowany soliton).
- `stage2_boost.py` — test bezstratności.
- `stage3_scatter.py`, `stage3b_largeb.py`, `stage3c_static.py` — falsyfikator §4 (3 testy).
- `debug_energy.py`, `debug_profile.py` — walidacja kodu energii i profilu (punkty odniesienia).

**Powiązane:** `RAMY_ladunek-kolor-nie-niezalezne_2026-07-28.md` (§4, §7 — zrealizowany falsyfikator).
