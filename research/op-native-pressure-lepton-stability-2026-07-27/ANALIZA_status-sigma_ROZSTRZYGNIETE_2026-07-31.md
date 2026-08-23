# Status σ — ROZSTRZYGNIĘTE: to dwa różne obiekty, oba dają kierunek

> # ⚠ CZĘŚCIOWO WYCOFANE 2026-08-10 — patrz `WYNIKI_szczebel0_sigma-nie-ma-prozni_2026-08-10.md`
> **§2 i §5 („filar spinu ma nośnik w σ_ab", „RP² = rozmaitość próżni σ_ab") są BŁĘDNE.**
> Ising łamie ℤ₂ **wewnętrzną**, nie łamie **obrotów** ⟹ C_x=C_y=C_z exact ⟹ **σ_ab ≡ 0
> w jednorodnej próżni**. Obiekt bez wartości próżniowej **nie ma rozmaitości próżni**,
> więc nie klasyfikuje defektów. σ_ab żyje **tylko na gradientach** — tak jak σ^ij.
> **Zostaje w mocy:** §1 (σ_ab to obserwabla poziomu 0, dwa różne „σ") i §3 (ℤ₃ nie z rangi 2).
> **Status filaru spinu: OTWARTY** — padła jedna droga, nie filar.

**Data:** 2026-07-31
**Typ:** ROZSTRZYGNIĘCIE + **definitywne wycofanie mojego wyniku „TGP nie ma defektów"**
**Metoda:** `core/sek08_formalizm` (tabela poziomów hierarchii, thm:metric-form-uniqueness),
`core/_meta_latex/status_map.tex`, `core/sek08c` (RECOVERY 2026-05-09)
**Cel programu (ustalony przez użytkownika):** najpierw **struktury z KOLOREM**, potem
**ładunek z ruchu względnego** — ładunek jako pochodna koloru.

---

## 1. Rozstrzygnięcie: są DWA różne obiekty „σ"

| | **σ_ab** (poziom 0) | **σ^ij** (gradient composite) |
|---|---|---|
| Definicja | **⟨ŝ·ŝ₊â⟩^TF** — bezśladowy korelator wiązań substratu | (∂^iΦ)(∂^jΦ) − ⅓δ^ij(∇Φ)² |
| Pochodzenie | **obserwabla substratu Γ, poziom 0 — ta sama warstwa co Φ** | zbudowany z gradientów Φ |
| Czy nowe pole? | **NIE** — inny korelator TEGO SAMEGO substratu | NIE — funkcja Φ |
| Czy wynika z Φ? | **NIE** (Φ=⟨ŝ²⟩ to korelator 1-punktowy; σ_ab jest 2-punktowy) | TAK |

**Cytat rozstrzygający** (`sek08_formalizm`, tabela poziomów):
> **Poziom 0 — Substrat Γ=(V,E):** „Hamilton H_Γ, symetria ℤ₂, coarse-graining;
> **Φ = ⟨ŝ²⟩, σ_ab = ⟨ŝ ŝ₊â⟩^TF**"
>
> **Poziom 2 — Metryka efektywna g_μν(Φ, **σ_ab**)**

⟹ **σ_ab NIE jest doklejonym polem.** Jest **drugą obserwablą tego samego substratu**:
Φ to lokalna amplituda (korelator 1-punktowy), σ_ab to **korelator kierunkowy (wiązaniowy)**.
Zapis w `status_map` („Hipoteza; nie wynika z samego Φ") jest **formalnie poprawny** —
σ_ab faktycznie nie wynika z Φ — ale **mylący**, bo sugeruje rozszerzenie spoza substratu.
Napięcie w rdzeniu: **pozorne**, wynika z dwóch znaczeń litery σ + niefortunnego opisu statusu.

## 2. Definitywne wycofanie mojego wyniku negatywnego

**Cytat, który przesądza** (`sek08_formalizm`, thm:metric-form-uniqueness):
> „W **minimalnym sektorze skalarnym** TGP (bez rozszerzenia fazowego θ_i, tj. **bez pola σ_ab**)
> każda metryka efektywna… musi być **diagonalna i izotropowa**."

Czyli metryka konforemnie płaska g_ij = h·δ_ij, na której oparłem całe rozumowanie
„reper = cechowanie ⟹ brak defektów", obowiązuje **wyłącznie w sektorze BEZ σ_ab**.
**Poza tym sektorem metryka nie jest izotropowa i widzi kierunek.**

> ⛔ **`ANALIZA_czy-TGP-ma-defekty_2026-07-31.md` — WYCOFANA.** Jej teza („TGP nie ma pola
> zdolnego nieść defekty") jest **fałszywa**: dotyczyła okrojonego sektora, nie TGP.

**Użytkownik miał rację na obu poziomach:** orientacja jest relacyjna (z gradientu / z korelatora
wiązań), i nie da się jej rozpatrywać w izolacji.

## 3. Co to znaczy dla KOLORU (cel programu)

**σ_ab jest bezśladowym, symetrycznym tensorem rangi 2 — czyli dokładnie NEMATYCZNYM
parametrem porządku.** To ma twarde konsekwencje dla dostępnych rozmaitości próżni:

| Typ porządku (ranga 2) | Rozmaitość próżni | π₁ | Czy niesie ℤ₃ (kolor)? |
|---|---|---|---|
| jednoosiowy (uniaxial) | **RP² = S²/ℤ₂** | ℤ₂ | **NIE** |
| dwuosiowy (biaxial) | SO(3)/D₂ | **Q₈** | **NIE** (Q₈ nie ma elementu rzędu 3 — sprawdzone na początku sesji) |

**Kluczowa obserwacja pozytywna:** rozmaitość **RP²** to dokładnie to, czego używa
**why_n3 Phase 3** do wyprowadzenia **spinu-½** (π₂(RP²), faza Berry'ego γ=π).
⟹ **Filar spinu ma nośnik w σ_ab.** To **odwołuje mój alarm** z poprzedniej analizy:
spin-½ nie wisi na cechowaniu — wisi na nematycznym parametrze porządku poziomu 0.

**Kluczowa obserwacja negatywna (dla koloru):** tensor **rangi 2 nie może dać porządku
tetraedrycznego**. Porządek T (tetraedryczny), którego π₁ = 2T zawiera ℤ₃, wymaga
**tensora rangi 3** (faza tetraedratyczna w ciekłych kryształach jest opisywana tensorem
trzeciego rzędu). σ_ab = ⟨ŝ ŝ₊â⟩ jest **dwupunktowy**, więc rangi 2.

> **Wniosek dla celu „struktury z kolorem": ℤ₃ NIE wypada z σ_ab.**
> Potrzebny jest **korelator 3-punktowy** ⟨ŝ ŝ₊â ŝ₊b̂⟩ (ranga 3) — albo inne źródło ℤ₃.

## 4. Konkretne, dobrze postawione pytanie na dalej

> **[2026-08-10] WSZYSTKIE TRZY TROPY SPRAWDZONE — patrz `WYNIKI_trzy-tropy-koloru_2026-08-10.md`.**
> Trop 1 **MARTWY** (ℤ₂ twarda + inwersja sieci zabijają rangę 3, dwa niezależne powody).
> Trop 2 **MARTWY** — **BRAK POWIĄZANIA** (GL(3,𝔽₂) jest perfekcyjna ⟹ brak charakterów ℤ₃).
> Trop 3 **POTWIERDZONY**: ℤ₂ jest twarda (nośna dla Φ≥0 i klasy 3D Ising).


Cel: struktura z kolorem. Trzy tropy, wszystkie sprawdzalne u źródła:

1. **Czy substrat Γ dopuszcza korelator 3-punktowy** jako obserwablę poziomu 0
   (analogicznie do Φ i σ_ab)? Jeśli tak — ranga 3 → porządek tetraedryczny → π₁ ⊃ ℤ₃ → **kolor
   z tej samej maszynerii, co spin**. To byłaby realizacja osi §0 (spin i kolor z jednego substratu).
2. **Czy ℤ₃ z sektora cechowania** (`sek09`: ℤ₃ ⊂ GL(3,𝔽₂), „ℤ₃ 't Hooft ⟹ N=3 minimalne")
   jest tym samym ℤ₃, czy niezależnym? Jeśli tym samym — most już istnieje.
3. **Czy symetria ℤ₂ substratu** (H_Γ ma ℤ₂) da się rozszerzyć/złamać tak, by dać ℤ₃ —
   czy ℤ₂ jest twarda? To ograniczenie może wykluczyć trop 1.

## 5. Stan po rozstrzygnięciu

**Wycofane:** „TGP nie ma defektów" (fałszywe — dotyczyło sektora bez σ_ab).
**Odwołane:** alarm o filarze spinu (RP² = rozmaitość uniaxial σ_ab; nośnik istnieje).
**Ustalone:** σ_ab to obserwabla poziomu 0, nie doklejone pole; kierunek jest w substracie natywnie.
**Nowe ograniczenie (twarde):** kolor ℤ₃ **nie wypada z tensora rangi 2** — potrzebna ranga 3
lub inne źródło. To jest **konkretna, falsyfikowalna przeszkoda**, a nie mgła.

---

**Źródła (u źródła):** `core/sek08_formalizm` (tabela poziomów: Φ=⟨ŝ²⟩, σ_ab=⟨ŝŝ₊â⟩^TF;
thm:metric-form-uniqueness: izotropia TYLKO bez σ_ab); `core/_meta_latex/status_map.tex`;
`core/sek08c` RECOVERY 2026-05-09 (σ^ij gradient composite); `research/why_n3` (RP² → spin-½)
