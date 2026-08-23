# §5.4 background-free: uzgodnienie wstążki z kanonicznym substratem TGP

**Data:** 2026-07-30
**Typ:** USTALENIE FUNDAMENTALNE (research u źródła: sek08a, sek08c) + auto-audyt
**Wniosek:** moja konstrukcja wstążki (krok 1–2) ŁAMIE §5.4. Wymaga re-gruntowania.

---

## 0. Co ustaliłem z rdzenia (u źródła)

Przeczytane: `core/sek08a_akcja_zunifikowana.tex`, `core/sek08c_metryka_z_substratu.tex`.

- **Φ = pole skalarne** (coarse-grained z sieci węzłów substratu Γ; def:info-budget, ψ=Φ/Φ₀).
  Jedyna zmienna dynamiczna (rem:not-scalar-tensor pkt 4).
- **Metryka g_eff jest ALGEBRAICZNIE wyznaczona przez Φ** — NIE niezależna zmienna
  (rem:not-scalar-tensor pkt 1). g_ij=h(ψ)δ_ij, g_tt=−c₀²f(ψ), warunek budżetu **fh=1**
  (prop:antipodal-from-budget). **To jest maszyneria background-free (§5.4) i ISTNIEJE w rdzeniu.**
- **Akcja na samo-generowanej metryce**: S=∫√(−g_eff)[½K(ψ)(∇ψ)²−V(ψ)+L_mat], K=ψ⁴ (α=2).
  Miara √(−g_eff)=c₀ψ/(4−3ψ) (M9.1'') zależy od Φ. Nie ma płaskiego tła.
- **Filar spinu**: S³ z **vielbeinu** metryki g_ij=ψδ_ij (hedgehog e^a_i=√g R^a_i), B=1;
  czysto topologiczne (zależne tylko od warunków brzegowych, nie od profilu/energii).

## 1. Auto-audyt: moja wstążka łamie §5.4

| Element mojej konstrukcji | Problem wobec kanonicznego TGP |
|---|---|
| Φ: ℝ³→M=S³/2T (target order-parameter) | Kanoniczne Φ jest SKALARNE; S³ ma pochodzić z **reperu** metryki, nie z zewnętrznego targetu |
| Płaskie ℝ³ (miara r²dr, płaski Laplasjan) | Łamie §5.4: TGP nie ma tła; miara/metryka generowane przez Φ |
| Ręcznie dodany człon Skyrme'a (krok 2) | Niemininalny (§5.3); wprowadzony, by pokonać Derricka — a **Derrick to artefakt płaskiego tła** |

> **Werdykt: krok 2 (skyrmion na płaskim ℝ³ + Skyrme) NIE jest TGP-natywny.**
> Metastabilność z topologii (krok 1, π₃/framing) pozostaje sensowna jako WZORZEC, ale
> osadzenie energetyczne było background-zależne. Nie buduję dalej na kroku 2 w tej formie.

## 2. Co przeżywa, co pada

**Przeżywa:**
- Szkielet grupowy (2T wymuszone przez warunek trialności) — **czysta teoria grup, bez tła.** ✅
- Filar spinu jako WZOREC: liczba kwantowa z reperu samo-generowanej metryki, bez energii. ✅
- Krok 1 jako WZORZEC topologiczny (nieusuwalność via niezmiennik), ale musi być przeniesiony
  z płaskiego π₃ na reper metryki.

**Pada / do przerobienia:**
- Krok 2 w obecnej formie (płaski skyrmion + Skyrme). Człon Skyrme'a jako input.
- Założenie, że Φ mapuje do zewnętrznego targetu (zamiast: 2T z reperu).

## 3. Twardo-otwarty problem (uczciwie)

- Kanoniczny skalarny TGP (F-A, α=2, W''(1)=+1) **nie ma zlokalizowanego solitonu** (runaway,
  AUDYT DODATEK A). Więc topologiczna liczba kwantowa w reperze **nie daje automatycznie
  zlokalizowanej cząstki** — liczbę może nieść rozmyta konfiguracja.
- „Co lokalizuje wstążkę bez łamania §5.4" jest **głównym otwartym pytaniem** domknięcia §5.4.

## 4. Fork domknięcia §5.4 (decyzja teoretyczna użytkownika)

- **(A) Reper samo-generowanej metryki (TGP-natywne).** Kolor (2T) jako topologia wiązki
  reperów g_eff[Φ], jak spin (S³, B=1). Zero energii, zero Skyrme, zero Derricka →
  §5.4+§5.1 natywnie. Musi zmierzyć się z „co lokalizuje cząstkę".
- **(B) Minimalne rozszerzenie ruchów Φ.** Uzasadnić człon stabilizujący (Skyrme-podobny) jako
  minimalny, background-free składnik reguł substratu (miara na g_eff, nie płaska). Daje
  lokalizację, kosztem minimalności.
- **(C) Inne** (np. lokalizacja z samo-generowanej metryki: g_eff modyfikuje skalowanie tak,
  że Derrick nie wymusza zapadania — do sprawdzenia).

---

**Pliki źródłowe:** `core/sek08a_akcja_zunifikowana.tex`, `core/sek08c_metryka_z_substratu.tex`
**Powiązane:** `WYNIKI_stepB2_*` (krok 2, teraz zdyskwalifikowany jako nie-natywny),
`WYNIKI_stepB_*` (krok 1, wzorzec), `research/qm_spin/README.md` (filar: reper→liczba kwantowa)
