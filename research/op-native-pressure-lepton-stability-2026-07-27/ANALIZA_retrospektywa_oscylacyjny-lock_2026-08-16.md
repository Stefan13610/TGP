# Retrospekcja całego procesu (2026-08-16) — co umykało: dwie maszynerie stabilności i oscylacyjny lock

**Cel główny (przypomnienie):** wyprowadzić stabilne obiekty kwarko/solitono-podobne z substratu.
**Pytanie sesji:** przeanalizować całość procesu — czy coś umyka.
**Odpowiedź:** tak, dwie rzeczy strukturalne. Obie sprawdzone u źródła, jedna zweryfikowana rachunkiem (4/4 PASS z kontrolą negatywną).

---

## 1. ŚLEPA PLAMKA: wszystkie testy stabilności testowały konfigurację, którą ontologia TGP uznaje za niefizyczną

Inwentarz testów stabilności w całym korpusie:

| Test | Konfiguracja | Warunek brzegowy | Wynik |
|---|---|---|---|
| AUDYT_KRYTYCZNY dod. A (F-A kanoniczne) | 1 soliton | u(∞)→1 (próżnia) | runaway |
| L03 #60 (CP-7) | 1 soliton | próżnia | siodła μ:2/τ:3 + kontinuum tachioniczne |
| L03 #62 (W1/W2) | 1 soliton | próżnia | deflacja niewystarczająca; kolaps τ |
| L03 #63 (V1–V3) | 1 soliton | próżnia | runaway g→0, t*≈1.7–3.6 |
| sek03/sek08 dwa-zrodla + 4 skrypty | 2 źródła | Φ→Φ₀ w ∞ | studnia = artefakt (audyt 2026-08-10) |
| moje σ_ab (szczebel 0) | próżnia jednorodna | — | σ≡0 (pytanie źle postawione w ramie TGP) |

**Każdy** test: pojedynczy obiekt lub para w pustym tle z próżnią w nieskończoności. Ontologia TGP mówi: *pojedynczy obiekt nie istnieje; stabilizuje go ciśnienie sąsiadów*. Czyli wszystkie „negatywne wyniki" są **zgodne z teorią** — a konfiguracja, o której teoria twierdzi, że jest stabilna (skończona gęstość źródeł, wzajemne blokowanie), **nie została policzona nigdy** — ani w rdzeniu, ani w audytach, ani przeze mnie.

Dowód ślepej plamki wprost w L03 POST_ACTION 2026-07-04, lin. 34–37: lista dróg wyjścia = „dyskretność substratu, inna symetria niż U(1), sektor sprzężony F-A, reinterpretacja metastabilnościowa". **Na liście nie ma „skończonej gęstości / sąsiadów"** — centralnego mechanizmu ontologii autora.

## 2. DWIE MASZYNERIE STABILNOŚCI w rdzeniu — nigdy niepołączone

- **Maszyneria 1 — EFT Φ** (sek03_rezimy, sek08 ssec:dwa-zrodla, 4 skrypty): profile monotoniczne (Yukawa/screened), trzy reżimy. **Padła w audycie 2026-08-10** (studnia = ekstrapolacja rozbieżnego szeregu / osobliwość −1/d; skale absurdalne; 17 PASS bez możliwości FAIL).
- **Maszyneria 2 — ODE** (O-L5/why_n3, status_map lin. 1283–1309, dodatekH p131–p146): solitony radialne F-A kanonicznego, generacje = liczba węzłów, **analityczny próg kolapsu** g₀,crit = 2(α+2)/[2(α+2)−d] = 8/5 (potwierdzony 1e−10), N_gen=3 z kolapsu (τ przy g₀=1.5696, 1.9% pod progiem; 4. generacja 58.7% nad progiem — zabroniona), oraz **ogon oscylacyjny: linearyzacja w f=g^(α+1) daje r. Bessela sferycznego z ω=1, uniwersalne (niezależne od α)** — „oscylacje, nie zanik wykładniczy". DodatekH lin. 1149: „Stabilność 3D (Derrick): **nie stosuje się** (oscylujący ogon)".

**Nigdzie w korpusie te maszynerie się nie spotykają:** wszystkie rachunki dwuźródłowe (maszyneria 1) używały profili monotonicznych, podczas gdy maszyneria 2 wyprowadza ogon **oscylujący**. Nikt nie policzył oddziaływania dwóch źródeł z ogonem, który rdzeń faktycznie wyprowadza.

## 3. NOWY WYNIK: oscylujący ogon ⟹ dyskretna drabina odległości blokowania

Teoria liniowa (ważna w ogonie |δΦ|≪Φ₀ — czyli **dokładnie w reżimie, w którym funkcjonał wg audytu dwa-zrodla JEST wiarygodny**, Φ≳0.8Φ₀):

E_int(d) ∝ −e^(−κd)·cos(d+φ)/d

Skrypt `RETRO_oscillating_tail_lock.py`, **4/4 PASS**:
- **[T1]** dla κ∈{0, 0.1, 0.3, 0.5}: ≥4 stabilne minima (V″>0) w d∈(2,30); pierwsze d*≈6.0–6.1
- **[T2] kontrola negatywna:** Yukawa (bez cos) — **0 minimów** przy każdym κ ⟹ mechanizm pochodzi z oscylacji, nie z procedury
- **[T3]** odstęp minimów = **2π** z dokł. ≤1.3% (korekta 1/d maleje z d). *Korekta uczciwości: pierwsza wersja testu oczekiwała π — błąd autora testu (minima −cos/d są co 2π, między nimi maksima); poprawione przed użyciem wyniku, pierwotny FAIL udokumentowany w skrypcie.*
- **[T4]** łańcuch 3 źródeł: stabilna równowaga równoodległa (a=b=6.083, Hessian dodatnio określony: 0.089, 0.139)

**Co to daje:**
1. **„Lock" bez bezdennej studni** — dyskretne, stabilne odległości równowagi; dokładnie „wzajemne blokowanie konfiguracji" z ontologii balonów.
2. **Odwrócenie problemu ważności:** stara studnia żyła w obszarze, gdzie funkcjonał nie obowiązuje; lock oscylacyjny żyje w ogonie, gdzie obowiązuje.
3. **Skala nie jest kosmologiczna:** ω=1 w jednostkach rdzenia solitonu ODE ⟹ d* ≈ 6 r_core i drabina co 2π·r_core. To potencjalnie **rozpuszcza blokadę Φ₀** (d*=4β=1e26 m było artefaktem maszynerii 1 kalibrowanej ciemną energią).
4. Sieć wielu źródeł = **kryształ z ciśnieniem sąsiadów** — konfiguracja, której brak w inwentarzu z §1.

## 4. Reinterpretacja „tła tachionicznego" (KONCEPT, nie derywacja)

Oscylujący ogon statyczny ⟺ jednorodne tło przy/za progiem niestabilności (L03: „kontinuum tachioniczne", „tło duchowione" przy ω_gh=0.2935). W standardowej QFT — wada dyskwalifikująca. W ramie TGP („nie ma próżni") — możliwa **treść teorii**: jednorodna pusta przestrzeń NIE jest stanem stabilnym; stabilne są konfiguracje wypełnione źródłami (krystalizacja / pattern formation, analogia: punkt Lifshitza, RKKY/Friedel). Flaga: to reinterpretacja; wymaga pokazania, że niestabilność jednorodnego tła saturuje się na sieci źródeł, a nie ucieka dalej.

## 5. Zastrzeżenia (uczciwie)

- Teoria liniowa + suma parowa; κ i φ **nie wzięte z faktycznego ODE** (fazy istnieją w rdzeniu: δ_e=−81.4°, δ_μ=+38.6°, δ_τ=−27.3° — dodatekH lin. 1126–1129; A_g≈|g₀−1|). Różne fazy gatunków ⟹ potencjalnie różne d* dla par ee/eμ/μτ — niezbadane.
- Ogon jest ogonem solitonu, którego **dynamiczna** niestabilność (runaway g→0) jest potwierdzona w izolacji. Lock nie dowodzi jeszcze, że sieć podnosi mod niestabilny do ω²>0 — to jest **TEN brakujący rachunek**.
- Maszyneria 2 nie przeszła u mnie audytu takiego jak maszyneria 1 (a ma znane rysy: Q_K=3/2 jako wejście, korekta r₃₁ w dodatekH p134e–g). Nie budować na niej bez audytu.

## 6. Rekomendacja (kolejność)

1. **Rachunek centralny:** powtórzyć test V3 (runaway) dla solitonu **w kąpieli sąsiadów** — sieć/klatka periodyczna z zadaną gęstością n, ogony oscylacyjne z faktycznego ODE. Pytanie binarne: czy mod runaway dostaje ω²>0 przy jakimś n. To jest bezpośredni test głównego celu i ontologii jednocześnie.
2. Wyciągnąć (κ, φ, A) ogona z faktycznego rozwiązania ODE rdzenia (nie z modelu jak tutaj) i przeliczyć d* dla rzeczywistych par.
3. Audyt maszynerii 2 (jak zrobiony dla maszynerii 1) — zanim cokolwiek na niej stanie.

**Skrypt:** `RETRO_oscillating_tail_lock.py` (4/4 PASS, kontrola negatywna T2, korekta T3 udokumentowana).
**Źródła rdzenia:** `status_map.tex` 1283–1316; `dodatekH_lancuch_wyprowadzen.tex` 1120–1161; `audyt/L03_K_phi_stability/POST_ACTION_UPDATE_2026-07-04.md`; `dodatekE_pi1_formal.tex` 55–115.
