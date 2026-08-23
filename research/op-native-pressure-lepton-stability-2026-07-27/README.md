---
title: "op-native-pressure-lepton-stability — od ciśnienia Goldstone'a (N4d) przez łańcuch substrat→wstążki→kwarki po oscylacyjny lock"
type: research_cycle
status: OPEN-ACTIVE
phase: "po retrospektywie 2026-08-16; rachunek centralny (V3 w kąpieli sąsiadów) NIEWYKONANY"
folder_status: open-mixed-verdict
claim_status: "OPEN 2026-08-16 — cykl wieloetapowy z serią twardych audytów. OBALONE: N4d 'native pressure' w izolacji (E[u]≥0 z równością iff u≡1 w obu sektorach kanonicznych — strukturalnie niewykonalne), 'pressure+loops=111%' (overfitting: jawny skan scale/λ_loop), bounce-hierarchy (N_neg=artefakt pudła 12/19/25=floor(R/π); F-A kanoniczna: runaway dla wszystkich g₀), cała warstwa budżetowa (h≡1 artefakt+bug; lokalizacja=artefakt UV; B=2 z rdzenia ⟹ |Δf|_max=0), kolor ℤ₃ z substratu Isinga (rank-3 znika tożsamościowo; GL(3,𝔽₂) perfekcyjna), σ_ab bez próżni (|σ|~L^−2.03), reżim III/studnia (trzy niezgodne rachunki rdzenia; ~17 PASS bez możliwego FAIL; studnia=osobliwość −1/d). PRZEŻYŁO: uniqueness 2T (jedyna skończona podgrupa SU(2) nieabelowa z ℤ₃ w abelianizacji), spin bezbarwny (−1∈Q₈ ⟹ χ(−1)=1), E[u]≥0 (twierdzenie), d*=4β (odporne na usunięcie E_γ, 13%), oscylacyjny lock (4/4 PASS z kontrolą negatywną: dyskretna drabina minimów co 2π·r_core, d*≈6.0–6.1, stabilny łańcuch 3 źródeł). ZIDENTYFIKOWANA ŚLEPA PLAMKA: żaden test stabilności w całym korpusie nie policzył konfiguracji o skończonej gęstości źródeł — tej, o której ontologia TGP twierdzi, że jest stabilna."
created_date: 2026-07-27
closed_date: null
authorization: "kontynuacja po #63 (wybór ścieżki N4d z NEEDS N4 cyklu op-nonlinear-charge-constraint); README zrekonstruowane retrospektywnie 2026-08-23 z artefaktów cyklu (cykl nie miał Phase0_balance ani README od startu)"
anti_lakatos_lock: "BRAK klasycznego LOCKa Phase 0 — cykl eksploracyjny; dyscyplinę zapewniały audyty adwersarialne post-hoc (AUDYT_KRYTYCZNY 2026-07-28, AUDYT_blok-budzetowy 2026-07-31, audyty trzech reżimów 2026-08-10/15) i jawne WYCOFANIA w nagłówkach dokumentów"
related:
  - "[[../op-nonlinear-charge-constraint-2026-07-03/NEEDS.md]]"
  - "[[../op-spectral-analysis-Phi-2026-07-03/README.md]]"
  - "[[../op-ep-scattering-babyskyrmion-2026-07-28/WYNIK_ep-scattering-babyskyrmion_2026-07-28.md]]"
  - "[[../../audyt/L03_K_phi_stability/POST_ACTION_UPDATE_2026-07-04.md]]"
  - "[[ANALIZA_retrospektywa_oscylacyjny-lock_2026-08-16.md]]"
---

# op-native-pressure-lepton-stability (2026-07-27 → 2026-08-16, OPEN)

> **README zrekonstruowane retrospektywnie 2026-08-23.** Cykl powstał bez README
> i bez LOCKa Phase 0; ewoluował przez trzy przeorientowania celu. Poniższa mapa
> odtwarza chronologię z dokumentów źródłowych. Dokumenty oznaczone ⛔ są
> WYCOFANE (obalone własnym audytem) — zachowane jako rejestr ślepych zaułków.

## Cel

- **Pierwotny (2026-07-27, ścieżka N4d z NEEDS #63):** czy naturalne sprzężenie
  Goldstone'a między solitonami e/μ/τ tworzy potencjał ciśnienia V_bg zdolny
  przesunąć ujemne wartości własne μ (2) i τ (3) z CP-7 ku dodatnim.
- **Po resecie (2026-07-28, [[RAMY_ladunek-kolor-nie-niezalezne_2026-07-28.md]]):**
  ładunek i kolor nie są niezależne — muszą wypływać z JEDNEJ topologii; łańcuch
  substrat → wstążki → kwarki → solitony wg checklisty §5.
- **Finalny (retrospektywa 2026-08-16):** wyprowadzić stabilne obiekty
  kwarko/solitono-podobne z substratu; kluczowy brakujący rachunek —
  stabilność w konfiguracji o **skończonej gęstości sąsiadów**.

## Chronologia i werdykty

### Etap 0 (2026-07-27) — native Goldstone pressure, Phases 1–4 — ⛔ OBALONE 07-28
Ekstrakcja ładunków z ogona (q_e=1.200, q_μ=1.107, q_τ=1.049), równowaga
3-ciałowa (E_press=−5.046), hybryda N4d+N4c „111.7% celu, λ_τ→+0.48 STABLE".
Dokumenty: [[README_WZNOWIENIE.md]], [[CHECKPOINT_PHASE4_EXTENDED_2026-07-27.md]],
[[TIER2_COMPLETE_MECHANISM_EXPOSITION.md]], [[PHASE4_EXTENDED_QUICKSTART.md]],
[[PARALLEL_WORK_MANIFEST_2026-07-27.md]].

### Etap 1 (2026-07-28) — bounce-hierarchy → AUDYT KRYTYCZNY → reset programu
Rano: hipoteza bounce-hierarchy (korelacja r=−0.9065, „95% publication-ready")
— [[BOUNCE_HIERARCHY_COMPLETE_EXPOSITION.md]] ⛔. Po południu **AUDYT KRYTYCZNY**
(`TGP/TGP_v1/research/.../AUDYT_KRYTYCZNY_2026-07-28.md` — uwaga: w zagnieżdżonej
ścieżce, patrz §Uwagi): N_neg = artefakt pudła (12/19/25 dla R=40/60/80,
identycznie dla próżni); realny sygnał to mody poniżej krawędzi −1 (e:0, μ:2, τ:3);
„111%" = overfitting (jawny skan scale∈{1,…,10}); **F-A kanoniczna nie ma
solitonów crown w ogóle** (W″(1)=+1>0, runaway dla wszystkich g₀). Reset:
[[RAMY_ladunek-kolor-nie-niezalezne_2026-07-28.md]] (rejestr 6 obalonych zaułków,
filar spin-½, checklista §5).

### Etap 2 (2026-07-30) — substrat → wstążki: warstwa grupowa 2T — ✅ PRZEŻYŁA AUDYTY
- **Uniqueness 2T** ([[WYNIKI_substrat-wstazki_2T_2026-07-30.md]]): wśród
  wszystkich skończonych podgrup SU(2) tylko binarna tetraedralna 2T=Q₈⋊ℤ₃
  (rz. 24) jest nieabelowa z ℤ₃ w abelianizacji.
- Ochrona przez π₃, nie π₂ ([[WYNIKI_stepB_dynamika-topologia_2026-07-30.md]]);
  M=S³/2T: π₁=2T, π₂=0, π₃=ℤ.
- Skyrmion płaski ⛔ zdyskwalifikowany (łamie §5.4 background-free) —
  [[WYNIKI_stepB2_dynamika-rozmiar_2026-07-30.md]],
  [[WYNIKI_s54_background-free_reconciliation_2026-07-30.md]]; domknięcie
  frame-native: [[WYNIKI_s54_domkniecie_frame-native_2026-07-30.md]].
- Lock kolor↔ładunek ℤ₃ {0,1/3,2/3} ([[WYNIKI_charge-lock_Z3_2026-07-30.md]]) —
  ⚠ matematyka OK, ale bliski tautologii (Hom(G,U(1))=Hom(G^ab,U(1)) zawsze);
  realnie DERIVED tylko uniqueness 2T + **spin bezbarwny**.

### Etap 3 (2026-07-31, uruchomienia 08-09) — łańcuch do kwarków + upadek warstwy budżetowej
- [[STATE_lancuch-substrat-kwarki_2026-07-31.md]] — mapa strzałek, punkt wznowienia.
- [[WYNIKI_domykanie-solitony_2026-07-31.md]] — konfinement topologiczny (0/6 klas
  domyka się sama); „21/30" później obalone (po usunięciu klas centralnych 13/13).
- **[[WYNIKI_cisnienie-lokalizacja_2026-07-31.md]] — wynik negatywny MOCNY:**
  E[u] ≥ 0 z równością iff u≡1 w obu sektorach kanonicznych ⟹ ciśnienie
  sąsiadów w tym sektorze NIE lokalizuje; „pressure+loops" musiało być fitem.
- [[ONTOLOGIA_energia-relacyjna-budzet_2026-07-31.md]] — energia relacyjna,
  prymitywem budżet; ontologia stoi, ale hipoteza budżetowa z niej wywiedziona ⛔:
  [[WYNIKI_budzet-skala_2026-07-31.md]] ⛔, [[WYNIKI_profil-z-budzetu_2026-07-31.md]] ⛔,
  [[WYNIKI_psi-orientacja_KOREKTA_2026-07-31.md]] ⛔,
  [[WYNIKI_koszt-topologiczny_2026-07-31.md]] ⚠ (bound symplicjalny T≥5B zostaje
  jako czyste twierdzenie).
- **[[AUDYT_blok-budzetowy_2026-07-31.md]] — blok budżetowy UPADA W CAŁOŚCI**
  (24 ustalenia dwóch audytorów: h≡1 artefakt konstrukcji + bug martwej zmiennej;
  lokalizacja = artefakt UV, R∝stała sieci; rdzeń ustala B=2 ⟹ |Δf|_max=0;
  R_min∝Q^(1/3) tautologia).
- [[ANALIZA_entropia-substratu_2026-07-31.md]] — s₀ to konstrukt kosmologiczny
  (skanowany vs DESI), nie prymityw. [[ANALIZA_miara-ze-skali_wieza-v2_2026-07-31.md]]
  — Φ₀ jest jawnie „scale-dependent free parameter"; miary ze skali obecnie
  wyprowadzić się NIE da. [[ANALIZA_czy-TGP-ma-defekty_2026-07-31.md]] ⛔ WYCOFANY;
  [[ANALIZA_status-sigma_ROZSTRZYGNIETE_2026-07-31.md]] ⚠ — σ_ab vs σ^ij to dwa
  obiekty; ℤ₃ nie wypada z tensora rangi 2.

### Etap 4 (2026-08-10) — kolor umiera; audyt trzech reżimów
- [[WYNIKI_trzy-tropy-koloru_2026-08-10.md]]: trop 1 martwy (rank-3 znika
  tożsamościowo, max|T_ijk|≤1e−14), trop 2 martwy (GL(3,𝔽₂) perfekcyjna ⟹ zero
  charakterów ℤ₃), trop 3 potwierdzony (ℤ₂ twarda). **Znalezisko: rdzeń ma DWA
  niezgodne substraty** (dodatekB: ŝ∈ℝ/ℤ₂ vs sek09: Ξ∈ℂ³/SU(3)_c) — kolor w TGP
  jest POSTULOWANY.
- [[WYNIKI_szczebel0_sigma-nie-ma-prozni_2026-08-10.md]]: σ_ab≡0 w próżni Isinga
  (|σ|~L^−2.03) ⟹ RP² bez domu w substracie; filar spinu OTWARTY, nie obalony.
- [[ANALIZA_trzy-rezimy_intuicja-gradientowa_2026-08-10.md]],
  [[AUDYT_ssec-dwa-zrodla_2026-08-10.md]], [[AUDYT_skrypty-trzy-rezimy_2026-08-10.md]]:
  studnia bezdenna i poza zakresem ważności; punkt „stable" w skrypcie rdzenia ma
  V″=−5.23 (maksimum); skale absurdalne (d_well protonu 17 rzędów pod Planckiem,
  M_crit=8×10¹⁹ M_☉); ocalałe: reżim I (grawitacja) i **d\*=4β**.

### Etap 5 (2026-08-15) — [[AUDYT_trzy-rezimy_ZAMKNIECIE_2026-08-10.md]]
**Studnia (reżim III/confinement) nie jest ustanowiona przez żaden z pięciu
rachunków rdzenia**; trzy liczą tę samą wielkość i są wzajemnie sprzeczne
(inny rząd skal, inna zależność od masy, inny mechanizm); ~17 „PASS" bez testu
zdolnego dać FAIL. `two_source_potential.py` dostaje trzy reżimy BEZ członu
kwartycznego — studnia = osobliwość −1/d. Blokada programu: **co ustala Φ₀
poza domeną kosmologiczną**.

### Etap 6 (2026-08-16) — [[ANALIZA_retrospektywa_oscylacyjny-lock_2026-08-16.md]] — NOWE OTWARCIE
1. **Ślepa plamka:** wszystkie testy stabilności w korpusie (#60–#63, audyty,
   ten cykl) badały pojedynczy obiekt/parę w próżni; konfiguracja o skończonej
   gęstości źródeł — ta, o której ontologia twierdzi, że jest stabilna —
   nigdy nie policzona.
2. **Dwie maszynerie stabilności nigdy niepołączone:** EFT Φ (profile
   monotoniczne; padła w audycie 08-10) vs ODE (ogon OSCYLACYJNY, próg kolapsu
   g₀,crit=8/5, N_gen=3).
3. **Nowy wynik (4/4 PASS, kontrola negatywna):** oscylujący ogon
   E_int(d)∝−e^(−κd)·cos(d+φ)/d ⟹ dyskretna drabina stabilnych minimów co
   2π·r_core, pierwsze d\*≈6.0–6.1; stabilny łańcuch 3 źródeł (Hessian dodatni).
   Skala NIE kosmologiczna — potencjalnie rozpuszcza blokadę Φ₀.
   Skrypt: `RETRO_oscillating_tail_lock.py`.

## Stan otwarty (krytyczna ścieżka wg retrospektywy §6)

1. **RACHUNEK CENTRALNY (niewykonany):** powtórzyć test V3 (runaway) dla
   solitonu **w kąpieli sąsiadów** (sieć periodyczna, gęstość n, ogony
   oscylacyjne z faktycznego ODE). Pytanie binarne: czy mod runaway dostaje
   ω²>0 przy jakimś n.
2. Wyciągnąć (κ, φ, A) ogona z faktycznego ODE rdzenia (fazy istnieją:
   δ_e=−81.4°, δ_μ=+38.6°, δ_τ=−27.3°); przeliczyć d\* dla par ee/eμ/μτ.
3. **Audyt maszynerii 2** (ODE/O-L5/why_n3) — znane rysy: Q_K=3/2 jako wejście,
   korekta r₃₁ w dodatekH. Nie budować na niej bez audytu.

Dodatkowo otwarte: co ustala Φ₀ poza kosmologią; nośnik filaru spin-½
(why_n3/Δ_a, qm_spin/π₃ — niezbadane); trzy rozwidlenia dla koloru; dwa
niezgodne substraty w rdzeniu; nadokreślenie więzów (n_sp+n_czas=B ∧
n_sp·n_czas=const ⟹ h≡1).

## Uwagi porządkowe

- **Zagnieżdżona ścieżka:** część plików z 2026-07-27/28 (m.in. kluczowy
  `AUDYT_KRYTYCZNY_2026-07-28.md`, `TIER2_PHASE_SUMMARY_2026-07-27.md`,
  `TIER2_NATIVE_MECHANISM_SUMMARY.md`, `meta/TIER2_SESSION65_PHASE4_RESULTS_2026-07-27.md`)
  leży w artefaktowej ścieżce `<cykl>/TGP/TGP_v1/…` — ten sam wzorzec, który
  sprzątano w sesji #62. Do przeniesienia na właściwe miejsca (housekeeping).
- Cykl siostrzany: [[../op-ep-scattering-babyskyrmion-2026-07-28/WYNIK_ep-scattering-babyskyrmion_2026-07-28.md]]
  (falsyfikator §4 RAMY: „ładunek" w baby-Skyrme NIE jest Coulombowski — negatyw
  mocny, zamyka trop).
- Pełny indeks plików .md z jednozdaniowymi opisami: patrz
  [[STATE_lancuch-substrat-kwarki_2026-07-31.md]] (do 07-31) oraz sekcja
  Chronologia powyżej (08-10 → 08-16).
