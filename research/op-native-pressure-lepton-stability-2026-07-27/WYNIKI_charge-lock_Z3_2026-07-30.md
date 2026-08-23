# Lock kolor↔ładunek: co realnie wypada, a co jest tautologią (mechanizm α)

**Data:** 2026-07-30 (skorygowane po audycie adwersarialnym 2026-07-31)
**Typ:** WYNIKI (strzałka wstążki→kwarki; oś §0)
**Weryfikacja:** `RIBBONS_charge_lock.py` + niezależny audyt (`AUDIT2_*.py`, SL(2,3)/GF(3))
**Status po audycie:** matematyka POTWIERDZONA (bez fitu, bez artefaktu). Ale nagłówek
„lock WYPADA strukturalnie" był PRZECENIONY. Realnie nowe: (i) uniqueness 2T, (ii) spin
bezbarwny. Sam „lock w trzecich" jest BLISKI TAUTOLOGII (patrz §AUDYT).

> ## ⚠ KOREKTA PO AUDYCIE (2026-07-31) — czytać najpierw
> Niezależny audytor odtworzył całą matematykę (osobno, SL(2,3)) — jest poprawna, bez fitu,
> bez artefaktu. Obalił jednak RETORYKĘ:
> 1. **„Lock" jest automatem, nie odkryciem:** Hom(G,U(1))=Hom(G^ab,U(1)) dla KAŻDEJ grupy —
>    każdy płaski ładunek U(1) MUSI faktoryzować przez abelianizację. „Ładunek=funkcja koloru"
>    to jedyne, co ładunek może robić, nie wyprowadzony lock.
> 2. **„3" to WEJŚCIE, nie wyjście:** wymaga DWÓCH wejść (kolor=ℤ₃ ORAZ ładunek=Hom(π₁,U(1))),
>    nie jednego. „Ładunek w trzecich" = ta sama ℤ₃ przez dualność Pontriagina.
> 3. **holonomia≠ładunek elektryczny:** Hom(π₁,U(1)) to fazy AB, nie ładunek cechowania (Gauss).
> **Realnie DERIVED (nietrywialne):** tylko (i) uniqueness 2T, (ii) −1∈[2T,2T] ⟹ spin bezbarwny.
> **§0 (lock kolor↔ładunek z jednej topologii) POZOSTAJE OTWARTE** — nie zamknięte tym wynikiem.

---

## 0. Nota o bugu (uczciwość)

Pierwszy przebieg dał „NIE Z3" (tylko 1 homomorfizm) i „spin −1 → trialność 1" — **sprzeczne
z ustalonym 2T^ab=ℤ₃**. To był **bug w kodzie**, nie w matematyce: funkcja `tri` etykietowała
warstwy przez arbitralny porządek zbioru, łamiąc addytywność ℤ₃. Poprawka: trialność jako
**prawdziwy homomorfizm 2T→ℤ₃** (log dyskretny wzgl. generatora), z asercją
`tri(gh)=tri(g)+tri(h) mod 3`. Po poprawce wynik zgodny z teorią. Zapisuję to, bo pokazuje,
że sprzeczny wynik → diagnoza kodu, nie naginanie narracji.

## 1. Test centralny

Ładunki U(1)_EM rozróżniające sektory dysklinacji (holonomie Aharonova-Bohma wokół wstążki)
= **Hom(π₁(M), U(1)) = Hom(2T, U(1)) = Hom(2T^ab, U(1)) = Hom(ℤ₃, U(1)) = ℤ₃.**

Zweryfikowane (`RIBBONS_charge_lock.py`):
- **3 charaktery** χ₀,χ₁,χ₂ (homomorfizmy 2T→U(1)) → **|Hom|=3 = ℤ₃**.
- Ładunki ułamkowe sektorów: **{0, 1/3, 2/3}** = 3 kolory.
- **spin −1: trialność 0 → ładunek ułamkowy 0** (bezbarwny — jak należy).
- elementy rzędu 3 (kolor/antykolor): trialności {1,2}.
- **Singlet kolorowy** (t=0+1+2=3≡0 mod 3) → ładunek ułamkowy **0 mod 1 = całkowity**.
  Pojedyncza kolorowa wstążka → ułamkowy. To sygnatura lock + konfinementu ułamka.

## 2. Co WYPADA (DERIVED, bez fitu)

- **Ładunek EM ułamkowy w TRZECICH** — z Hom(2T,U(1))=ℤ₃.
- **fractional(Q) = trialność/3** ⟹ **LOCK kolor↔ładunek STRUKTURALNY** (nie postulat grupy
  cechowania jak w SM). To jest oś §0 zrealizowana w części ℤ₃.
- **„3" pochodzi WYŁĄCZNIE z |2T/Q₈|=|ℤ₃|**, które było wymuszone trialnością (uniqueness 2T).
  Zero dopasowania do znanych ładunków.
- **Spin bezbarwny** (−1 → ładunek 0): nieabelowy rdzeń Q₈ nie niesie ładunku ułamkowego —
  checkowalne, nietrywialne, spełnione.
- **Singlety = ładunek całkowity**: 3 kolory sumują trialność do 0 mod 3.

## 3. Test na cyrkularność (jawnie)

- **Czy to tautologia „3 kolory → ładunek w 3"?** Nie do końca: nietrywialna treść to (a) ładunek
  = DOKŁADNIE abelianizacja ℤ₃, a nieabelowy Q₈ (spin) jest bezbarwny (charge 0) — checkowalne;
  (b) lock = fractional(Q) jest TĄ SAMĄ ℤ₃ co kolor, nie niezależny; (c) singlety→całkowity.
- **Czy lock jest wsadzony?** Nie: nie zakładam relacji Q↔kolor; wypada z kohomologii grupy
  Hom(π₁,U(1)). Wejściem jest tylko „kolor=ℤ₃".
- **Werdykt:** DERIVED-given-input, nie tautologia, nie fit. ALE to część ℤ₃, nie pełne §0.

## 4. Co NIE WYPADA (OTWARTE — nie fitować)

- **Dokładne wartości (+2/3 vs −1/3):** wymagają offsetu od nawinięcia B/2
  (Gell-Mann–Nishijima Q = B/2 + I₃-analog) + normalizacji U(1). Osobne wyprowadzenie.
  **Nie wolno tego dopasować do PDG** — to byłby powrót do fitu.
- **Emergencja EM z prądu topologicznego** (RAMY §3): identyfikacja „ℤ₃-ładunek topologiczny =
  ładunek elektryczny" jest KONCEPTEM, nie wyprowadzeniem. Wymaga własnej derywacji
  (skyrmion emergent-EM: ruch tekstury→E, nawinięcie→B).
- Uwaga: Hom(π₁,U(1)) klasyfikuje płaskie holonomie U(1); nazwanie tego „ładunkiem elektrycznym"
  czeka na §3-emergencję.

## 5. Jednozdaniowo (po audycie)

> Matematyka jest poprawna i bez fitu, ale realnie nowe/wyprowadzone jest tylko: uniqueness 2T
> i bezbarwność spinu (−1∈komutatorze); sam „lock ładunek↔kolor w trzecich" jest automatyczną
> konsekwencją Hom(G,U(1))=Hom(G^ab,U(1)) przy DWÓCH wejściach (kolor=ℤ₃, ładunek=holonomia),
> więc bliski tautologii — a identyfikacja holonomia=ładunek elektryczny jest nieudowodniona;
> **oś §0 pozostaje OTWARTA.**

## 6. AUDYT ADWERSARIALNY (2026-07-31) — pełny werdykt

Niezależny audytor (własna konstrukcja SL(2,3)/GF(3), pliki `AUDIT2_group_verify.py`,
`AUDIT2_uniqueness_probe.py`):

**POTWIERDZONE (bez zastrzeżeń):** |2T|=24, centrum ℤ₂, komutator Q₈, 2T^ab=ℤ₃,
Hom(2T,U(1))=ℤ₃, spin −1 bezbarwny, uniqueness 2T (odporne na cutoff n=1..15), korekta buga
legalna, brak fitu do PDG, brak artefaktu numerycznego.

**OBALONE/PRZECENIONE (retoryka):**
- „lock = fractional(Q)=trialność/3 to nietrywialny wynik" → **WĄTPLIWE**: automat z faktoryzacji
  charakterów przez abelianizację.
- „3 wyłącznie z |2T/Q₈|" → **przemilcza**, że |ℤ₃| było WŁOŻONE jako warunek wyboru 2T.
- „Hom(π₁,U(1)) = ładunek elektryczny" → **skok** (fazy AB ≠ ładunek Gaussa); oznaczone otwarte ✅.
- „§5.4 background-free reper" → **przemianowanie**; symetria T uzasadniona cyrkularnie (T bo kolor=ℤ₃).

**Klasyfikacja audytora (niezależna):**
- DERIVED solidne: uniqueness 2T; −1∈[2T,2T] ⟹ spin bezbarwny (jedyny nietrywialny fakt fizyczny).
- CYRKULARNE/TAUTOLOGICZNE: sam lock „w trzecich"; symetria T z żądanego wyniku.
- KONCEPT: holonomia→ładunek; reper→background-free.
- OTWARTE (uczciwie oznaczone): ±2/3,−1/3; emergencja EM; lokalizacja.

**Zgoda z „DERIVED-given-input, nie tautologia, nie fit":** „nie fit" TAK; „DERIVED" częściowo
(wymaga DWÓCH wejść); „nie tautologia" NIE dla samego locka.

---

**Pliki:** `RIBBONS_charge_lock.py`, `RIBBONS_s54_frame_structure.py`
**Powiązane:** `WYNIKI_s54_domkniecie_frame-native_2026-07-30.md` (2T=Q₈⋊ℤ₃),
`RAMY_ladunek-kolor-nie-niezalezne_2026-07-28.md` (§0, §3 EM emergentne, §4 falsyfikator)
