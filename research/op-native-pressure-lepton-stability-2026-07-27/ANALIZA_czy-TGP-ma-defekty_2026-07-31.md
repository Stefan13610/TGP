# Czy TGP w ogóle ma pole zdolne nieść defekty? — kontrola składników wieży v2

> # ⛔ WYCOFANE 2026-07-31 — patrz `ANALIZA_status-sigma_ROZSTRZYGNIETE_2026-07-31.md`
> **Teza tego dokumentu jest FAŁSZYWA.** Oparłem ją na metryce izotropowej g_ij=h·δ_ij,
> która wg `thm:metric-form-uniqueness` obowiązuje **wyłącznie w „minimalnym sektorze
> skalarnym… BEZ pola σ_ab"**. TGP ma na **poziomie 0** drugą obserwablę substratu:
> **σ_ab = ⟨ŝ ŝ₊â⟩^TF** (bezśladowy korelator wiązań) — nie jest to doklejone pole, tylko
> inny korelator tego samego substratu Γ. Niesie kierunek natywnie ⟹ **defekty są możliwe**,
> a metryka nie jest konforemnie płaska. **Alarm o filarze spinu również ODWOŁANY**:
> RP² (rozmaitość nematyka jednoosiowego σ_ab) jest dokładnie tym, czego używa why_n3.

**Data:** 2026-07-31
**Typ:** ANALIZA STRUKTURALNA (kontrola „czy mechanizm ma składniki") — **wynik NEGATYWNY**
**Metoda:** czytanie rdzenia u źródła (`sek08a` rem:not-scalar-tensor, `_meta_latex/status_map.tex`,
`tabela_epistemiczna.tex`)
**Powód:** przed budowaniem czegokolwiek na wieży v2 (Kibble-Zurek) trzeba sprawdzić,
czy TGP ma parametr porządku, który może mieć nietrywialną rozmaitość próżni.

---

## 0. Wynik w jednym zdaniu

> **W kanonicznym TGP jedyną zmienną dynamiczną jest RZECZYWISTY SKALAR Φ. Rozmaitość próżni
> takiego pola jest PUNKTEM, więc wszystkie grupy homotopii są trywialne — TGP kanoniczne
> NIE MA defektów topologicznych żadnego rodzaju.** Reper R jest czystym cechowaniem.
> To podcina warstwę wstążek/2T — jedyną, która przetrwała oba audyty.

## 1. Dowód z rdzenia (cytaty)

**(a) Jedno pole dynamiczne.** `sek08a`, rem:not-scalar-tensor, pkt 1 i 4:
> „metryka g_μν^eff jest **wynikowa**, zdefiniowana **algebraicznie** przez Φ
> (**nie jest niezależną zmienną dynamiczną**)."
> „**Jedyną zmienną dynamiczną jest Φ** (nie para (g_μν, Φ))."

**(b) Pole tensorowe σ_ab jest tylko HIPOTEZĄ.** `_meta_latex/status_map.tex` lin. 219–223:
> „σ_ab: mody tensorowe z tensora naprężeń substratu — **Hipoteza** —
> **rozszerzenie A1; NIE WYNIKA Z SAMEGO Φ**"

oraz `tabela_epistemiczna.tex` lin. 25–27: sektor tensorowy σ_ab „dodaje **+1 nieredukowalny
parametr swobodny** C_σ".

## 2. Konsekwencja topologiczna (elementarna, nie wymaga rachunku)

Dla rzeczywistego skalara Φ z potencjałem mającym minimum przy φ=1:

- **rozmaitość próżni M_vac = {φ = 1} = PUNKT**,
- ⟹ π₀ = π₁ = π₂ = π₃ = **trywialne**,
- ⟹ **brak ścian domenowych, brak strun/dysklinacji, brak monopoli, brak tekstur.**

**Przejście g: 0→1 (inflacja TGP) nie jest łamaniem symetrii** — to zwykłe staczanie się pola
do jedynego minimum. **Kibble-Zurek nie ma tu czego produkować:** mechanizm KZ wymaga
**zdegenerowanej** rozmaitości próżni (różne obszary wybierają **różne** próżnie). Tu próżnia
jest jedna, więc nie ma czego „wybierać" i nie ma na czym powstać defektowi.

## 3. Reper R jest cechowaniem, nie stopniem swobody

Metryka TGP: **g_ij = h(ψ)·δ_ij** — konforemnie płaska, wyznaczona przez skalar.
Reper e^a_i = √h·R^a_i daje **tę samą metrykę dla KAŻDEGO R ∈ SO(3)** (bo RᵀR = I —
zweryfikowaliśmy to numerycznie w `BUDGET_psi_orientation.py`).

Jeżeli teoria zależy **wyłącznie** od metryki i Φ (a wg §1 tak jest), to **R jest lokalną
redundancją cechowania**, a nie polem fizycznym. Wtedy:
- „dysklinacje w R" są **artefaktami cechowania**, nie obiektami,
- nawinięcie B jest **zależne od cechowania**, więc nie jest obserwablą.

**To samo ustalenie, które wcześniej uznałem za wygodne** („metryka nie widzi orientacji"),
**obraca się teraz przeciw całej konstrukcji**: skoro nic nie widzi R, to R nic nie znaczy.

## 4. Co to robi z dotychczasowymi wynikami

| Wynik | Nowy status |
|---|---|
| Uniqueness 2T (jedyna podgrupa SU(2) z ℤ₃ w abelianizacji) | **matematyka stoi**, ale jako twierdzenie **warunkowe**: „JEŚLI parametr porządku ma strukturę SO(3)/T, TO…". **TGP kanoniczne nie ma takiego parametru.** |
| Spin bezbarwny (−1 ∈ Q₈) | j.w. — warunkowe |
| Reguły domykania wstążek | j.w. — warunkowe, dodatkowo tautologiczne (audyt) |
| Wieża v2 / Kibble-Zurek | **brak składników** w kanonicznym TGP (próżnia niezdegenerowana) |
| **Filar spin-½ z nawinięcia vielbeinu** (`qm_spin`) | ⚠ **POWAŻNE PYTANIE** — jeśli reper jest cechowaniem na konforemnie płaskiej metryce, to nawinięcie B jest zależne od cechowania. **Wymaga osobnego rozstrzygnięcia** (druga ścieżka, RP²+Berry z why_n3, używa innego pola Δ_a i może nie być tym objęta). |

**Nie rozstrzygam tu losu filaru** — sygnalizuję, że argument, który podciął wstążki, **dotyka
też jego**, i że to wymaga własnego, ostrożnego audytu. To nie jest teza, to alarm.

## 5. Jedyne wyjścia (uczciwie wyliczone, bez rekomendacji)

1. **Przyjąć σ_ab jako pole fizyczne** — wtedy istnieje tensor z orientacją, możliwa
   nietrywialna rozmaitość próżni i defekty. **Koszt:** status „Hipoteza", „nie wynika
   z samego Φ", **+1 wolny parametr** (C_σ). Czyli defekty **kupione**, nie wyprowadzone —
   łamie §5.3 (minimalność) i §5.4 (background-free).
2. **Rozszerzyć Φ o strukturę wewnętrzną** (Φ wielokomponentowe / z symetrią G łamaną do H).
   Wtedy M_vac = G/H może mieć π₁ = 2T. **Koszt:** to jest **nowy postulat o polu**,
   czyli zmiana aksjomatu A1 — poważna ingerencja w rdzeń.
3. **Uznać, że cząstki nie są defektami Φ** i szukać ich gdzie indziej (np. w sektorze
   cechowania sek09, gdzie ℤ₃ ⊂ GL(3,𝔽₂) już występuje).
4. **Przyjąć, że reper JEST fizyczny** — wymaga mechanizmu, który go usztywnia (torsja,
   sprzężenie do materii fermionowej). W TGP fermionów się dopiero szuka ⟹ **ryzyko cyrkularności**.

## 6. Dlaczego to jest dobra wiadomość (mimo że negatywna)

- To jest **kontrola składników przed budową** — dokładnie to, czego zabrakło trzy razy
  wcześniej w tej sesji (płaski Skyrme, cap budżetowy, limit tempa). Tym razem sprawdziłem
  **zanim** policzyłem.
- Wynik jest **elementarny i odporny**: nie zależy od numeryki, siatki, pudła ani normalizacji.
  Nie da się go obalić przez „inny wybór wagi".
- **Zawęża pole poszukiwań** do czterech konkretnych wyjść (§5), zamiast pozwalać budować
  kolejne piętra na strukturze, której nie ma.

## 7. Jednozdaniowo

> Kanoniczne TGP ma jedno pole dynamiczne — rzeczywisty skalar Φ — więc jego rozmaitość próżni
> jest punktem, wszystkie grupy homotopii są trywialne, reper jest cechowaniem, a mechanizm
> Kibble-Zurka nie ma czego produkować; **warstwa wstążek/2T pozostaje poprawną matematyką,
> ale bez nośnika w kanonicznym TGP**, a pytanie o filar spinu wymaga osobnego audytu.

---

**Źródła (u źródła):** `core/sek08a_akcja_zunifikowana.tex` (rem:not-scalar-tensor pkt 1,4);
`core/_meta_latex/status_map.tex` lin. 219–223 (σ_ab = Hipoteza, nie wynika z Φ);
`core/_meta_latex/tabela_epistemiczna.tex` lin. 22–27 (bilans parametrów, C_σ wolny);
`BUDGET_psi_orientation.py` [A] (g = h·δ niezależnie od R)
