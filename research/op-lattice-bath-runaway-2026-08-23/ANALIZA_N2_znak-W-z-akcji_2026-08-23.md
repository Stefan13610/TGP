---
title: "ANALIZA N2 — znak W z akcji TGP: rdzeń zawiera OBA znaki (sek08a: wyprowadzenie z akcji vs Nota kanoniczna 2026-04-07); hipoteza autora 'pojedynczy obiekt vs kąpiel' ma dokładny odpowiednik strukturalny, ale bez derywacji"
date: 2026-08-23
type: analiza-post-close
tgp_owner: research/op-lattice-bath-runaway-2026-08-23
status: ROZSTRZYGNIĘTE-DOKUMENTACYJNIE (derywacja jednego znaku z akcji: OTWARTA)
verdict: "Źródło rozdwojenia znaku W zlokalizowane W SAMYM sek08a: (1) prop:field-eq-from-action (z akcji zunifikowanej, waga √−g_eff=c₀ψ) ⟹ próżnia STABILNA (m²=+γ, ogon Yukawa) — gałąź AUD/CP-7, brak solitonów; (2) inline 'Nota kanoniczna 2026-04-07' (lin. ~496–497): S[g]=∫[½g⁴(∇g)²+(β/7)g⁷−(γ/8)g⁸]d³x ⟹ EL = DOKŁADNIE równanie maszynerii 2 (źródło g²(1−g)), próżnia TACHIONOWA (W″(1)=−γ), ogon oscylacyjny, solitony z progiem 8/5. Nota kanoniczna paruje kanoniczną kinetykę K=g⁴ (L04) ze znakiem potencjału F-S — L04 canonicalization rozstrzygnął TYLKO kinetykę (α=2, K=g⁴), nigdy znak potencjału. Która akcja jest 'tą' akcją — NIEROZSTRZYGNIĘTE w rdzeniu."
related:
  - "[[PhaseA_report.md]]"
  - "[[ANALIZA_N1_pochodzenie-faz_2026-08-23.md]]"
  - "[[NEEDS.md]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
  - "[[../op-L04-ODE-canonicalization-2026-05-04/canonical_form_evidence.md]]"
  - "[[../op-native-pressure-lepton-stability-2026-07-27/ANALIZA_retrospektywa_oscylacyjny-lock_2026-08-16.md]]"
---

# N2 — znak W z akcji TGP

**Charakter dokumentu:** POST-CLOSE (śledztwo NEEDS N2 autoryzowane 2026-08-23).
Metoda: archeologia źródeł + algebra elementarna (jeden krok, jawnie niżej);
parowanie „akcja↔równanie" dla gałęzi (1) było już sympy-exact w CP-7 (#60).

## 1. Dwie akcje w samym sek08a

### Gałąź (1): akcja zunifikowana → próżnia STABILNA

`prop:field-eq-from-action` (sek08a, lin. ~471–520): wariacja akcji
z wagą √−g_eff = c₀ψ (potencjał wewnętrzny −(β/3)ψ³+(γ/4)ψ⁴) odtwarza
`eq:field-eq-reproduced`:

```
∇²Φ + 2(∇Φ)²/Φ + βΦ²/Φ₀ − γΦ³/Φ₀² = −qΦ₀ρ
```

Linearyzacja wokół ψ=1 (β=γ): ∂_ψ(βψ²−γψ³)|₁ = 2β−3γ = −γ ⟹
**∇²h = +γh ⟹ m² = +γ > 0: próżnia stabilna, ogon Yukawa e^(−√γ·r)/r.**
Zgodne z CP-7 #60: C1 (m_sp²=γ exact), C2 (próżnia F-A: N_neg=0).
W tej gałęzi statyczny profil z g₀>1 nie zawraca — **runaway, brak
solitonów korony** (AUDYT_KRYTYCZNY 2026-07-28 Dodatek A; potwierdzone
niezależnie w bramce A5 tego cyklu, 4/4).

### Gałąź (2): Nota kanoniczna 2026-04-07 → próżnia TACHIONOWA

Inline komentarz w sek08a (lin. ~496–497, wewnątrz dowodu
prop:field-eq-from-action!):

```
% kanoniczna: S[g] = ∫[½g⁴(∇g)² + (β/7)g⁷ − (γ/8)g⁸] d³x
```

EL tego funkcjonału (K=g⁴; jeden krok algebry):
U′ = βg⁶−γg⁷ ⟹ U′/K = βg²−γg³ = g²(1−g) przy β=γ=1 ⟹

```
g″ + (2/r)g′ + (2/g)g′² = g²(1−g)
```

— **DOKŁADNIE równanie maszynerii 2** (audytowane jako M2 w bramce A;
A1/A2 PASS: próg 8/5, ω=1). U″(1) = 6β−7γ = **−γ < 0: próżnia tachionowa,
ogon oscylacyjny, solitony istnieją dla g₀<8/5.**

### Sedno

Źródła obu gałęzi różnią się **względnym znakiem członu potencjałowego**
(βψ²−γψ³ po przeciwnych stronach równania). Obie żyją w TYM SAMYM
podrozdziale sek08a — gałąź (2) jako komentarz „kanoniczna" wewnątrz dowodu
gałęzi (1), bez zaznaczenia, że zmienia znak i fizykę próżni.

## 2. Skąd znak gałęzi (2): to znak potencjału F-S

CP-7 #60 (sympy-exact): „funkcjonał F-S (f=1+4ln g, **W′=g²(1−g)**) ⇔ ODE
korony a3d/ls10". Źródło g²(1−g) jest historycznie potencjałem **Formulacji
B/F-S** (W_FS=g³/3−g⁴/4, W_FS″(1)=−1 — stąd „kontinuum tachioniczne od −1"
w CP-7 C4/C5). Nota kanoniczna sparowała kanoniczną kinetykę K=g⁴
z tym właśnie znakiem źródła: **hybryda K(F-A) + znak W(F-S)**.

**L04 canonicalization (2026-05-04) rozstrzygnął wyłącznie kinetykę**
(thm:D-uniqueness: C1–C3 ⟹ α=2, K=g⁴ — [[../op-L04-ODE-canonicalization-2026-05-04/canonical_form_evidence.md]]).
**Znak potencjału nie był przedmiotem żadnego z trzech dowodów.** Dualizm
L04 nie został więc domknięty na poziomie potencjału — odrodził się jako
rozdwojenie znaku W.

## 3. Hipoteza autora: „pojedynczy obiekt vs kąpiel sąsiadów"

Hipoteza (2026-08-23): różnica znaku może odpowiadać temu, czy mierzymy
pojedynczy obiekt, czy obiekt w kąpieli.

**Ocena: korespondencja strukturalna JEST w korpusie, derywacji NIE MA.**

- Gałąź (1) (W″(1)>0): izolowany obiekt w **stabilnej próżni** — sektor
  grawitacyjny/słabopolowy (Yukawa m²=γ; CP-7: sektor grawitacyjny czysty).
- Gałąź (2) (W″(1)<0): **jednorodne tło niestabilne** — dokładnie obraz
  „nie ma próżni; stabilne są konfiguracje wypełnione źródłami"
  z retrospektywy 2026-08-16 §4 (krystalizacja, RKKY/Lifshitz); ogon
  oscylacyjny tej gałęzi jest statycznym odpowiednikiem tła
  tachionicznego i nośnikiem locka oscylacyjnego (drabina 2π).

Czyli: mapowanie {izolowany obiekt ↔ znak (1)}, {kąpiel/skończona gęstość ↔
znak (2)} jest spójne z całą strukturą wyników — **ale nigdzie w rdzeniu
nie ma rachunku, który wyprowadzałby zmianę znaku efektywnego W z przejścia
próżnia→skończona gęstość**. Taki rachunek (efektywny potencjał fluktuacji
wokół tła o gęstości n zamiast wokół ψ=1) byłby NATURALNYM elementem
rachunku centralnego N3 — i mógłby zdegradować problem „którą akcję wybrać"
do „które tło rozwijamy" (jedna akcja, dwa sektory).

## 4. Status i konsekwencje

1. **Rozstrzygnięte dokumentacyjnie:** źródło rozdwojenia zlokalizowane
   (sek08a: prop:field-eq-from-action vs Nota kanoniczna); maszyneria 2 MA
   kotwicę w rdzeniu (nota), ale ta kotwica jest w konflikcie ze znakiem
   wyprowadzenia z akcji w tym samym podrozdziale.
2. **Otwarte (istota N2):** derywacja JEDNEGO znaku z akcji + reguł
   kanonizacji, ALBO derywacja obu znaków jako dwóch sektorów jednej akcji
   (hipoteza kąpieli, §3). Do rozstrzygnięcia zanim maszyneria 2 będzie
   cytowana jako „kanoniczna".
3. Konsekwencja dla N1: fazy z p131 żyją w gałęzi (2) w wariancie
   logarytmicznym (eq:J-ode) — ich status dziedziczy status gałęzi (2).

## 5. Proponowane korekty rdzenia (USER-GATED, nie wykonane)

- sek08a: przekształcić inline „Notę kanoniczną" w jawny remark
  z flagą: EL noty = równanie solitonowe maszynerii 2; znak potencjału
  PRZECIWNY niż w eq:field-eq-reproduced; wybór znaku OTWARTY
  (odsyłacz do tej analizy).
- status_map O-L5: flaga przy „kanoniczna TGP ODE": źródło g²(1−g)
  odpowiada nocie kanonicznej, nie wyprowadzeniu z akcji zunifikowanej.
- Ewentualny nowy op (analytical decision, wzorzec L04): „znak W
  z akcji + hipoteza dwóch sektorów (próżnia vs kąpiel)".

## 6. POST-SCRIPTUM 2026-08-31 — rozstrzygnięcia po tej analizie

- **Hipoteza dwóch sektorów: OBALONA w klasie zbadanej**
  (op-bath-two-sectors, Phase 3, Q2-FAIL): w akcji gałęzi STABILNEJ
  samouzgodnione tła sieci źródeł d∈{8,6,4,2} dają ω²_min(d)
  = +1.34/+1.57/+1.88/+2.47 — dodatnie i ROSNĄCE z gęstością
  (kontrola d=∞ czysta: m²=γ±0.00%). Gęstość USZTYWNIA potencjał
  fluktuacji; znak tachionowy nie emerguje. **Wybór znaku W = otwarty
  problem AKSJOMATYCZNY (decyzja ontologiczna autora).**
- **Nota porównawcza (N5 op-substrate-fluctuation-channel, poziom 0 MFT):**
  na poziomie substratu znak tachionowy pochodzi z ROZRZEDZENIA
  (spinodala on-site: próg Φ_c/Φ_vac≈0.298 przy WF, skan 0.197–0.331),
  a wiązanie gradientowe ZAWSZE stabilizuje (ΔC_bond=+8zJs_b⁶≥0) —
  odwrotna rola gęstości niż w hipotezie dwóch sektorów; deskryptywnie
  spójne z Q2-FAIL poziomu 1 (oba poziomy: zagęszczanie stabilizuje).
- **Korekta rdzenia z §5: WYKONANA 2026-08-31** (user-gate): jawny
  remark `rem:W-sign-axiomatic` w sek08a (przy prop:field-eq-from-action),
  scalający lokalizację z tej analizy + Q2-FAIL + notę poziomu 0.
