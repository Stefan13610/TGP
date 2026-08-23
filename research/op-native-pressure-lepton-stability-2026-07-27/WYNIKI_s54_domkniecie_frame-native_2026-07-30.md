# §5.4 domknięte (ścieżka A) + karta wyników strzałki substrat→wstążki

**Data:** 2026-07-30
**Typ:** DOMKNIĘCIE §5.4 (background-free) + zamknięcie strzałki na poziomie strukturalnym
**Weryfikacja:** `RIBBONS_s54_frame_structure.py` + wcześniejsze skrypty tej sesji
**Decyzja:** ścieżka (A) — reper samo-generowanej metryki (użytkownik)

---

## 0. Rozwiązanie napięcia §5.4 (ontologiczne, nie nowa matematyka)

Napięcie: moja wstążka używała M=S³/2T jako **zewnętrznego** targetu na **płaskim** tle + Skyrme.
Rozwiązanie ścieżki (A): **to samo M=S³/2T, ale zinterpretowane jako przestrzeń konfiguracji
REPERU (frame) samo-generowanej metryki** g_ij=ψδ_ij (z rdzenia, background-free), z tetraedryczną
symetrią substratu T.

- **Matematyka identyczna** (π₁=2T, π₂=0, π₃=ℤ — zweryfikowane w krokach 1–2).
- **Ontologia teraz zgodna z §5.4:** frame bundle jest intrinsyczny dla metryki generowanej
  przez Φ; żadnego zewnętrznego pola/tła na wejściu. To DOKŁADNIE wzorzec filaru spinu
  (S³ z vielbeinu), rozszerzony o symetrię T.
- **Człon Skyrme'a ZNIKA:** liczby kwantowe (spin, kolor) to topologia reperu, nie minima
  energii → **brak Derricka, brak płaskiego tła, brak wkładanego członu.**

## 1. Spin i kolor z JEDNEGO reperu (zweryfikowane)

M = S³/2T (przestrzeń reperu):
- **π₃(M)=ℤ** → nawinięcie B = **spin** (jak filar; centralne −1 = generator FR).
- **π₁(M)=2T** → dysklinacje reperu = **kolor**.
- Oba z jednego M → oś §0 zrealizowana na poziomie substratu.

Struktura 2T (`RIBBONS_s54_frame_structure.py`): **2T = Q₈ ⋊ ℤ₃**
- komutator Q₈ (rząd 8) = bezbarwny rdzeń; **−1 (spin) ∈ Q₈, kolor(−1)=trywialny** → spin bezbarwny;
- abelianizacja ℤ₃ = kolor; **elementy rzędu 3 (kolor) NIE-centralne**, warstwy [1,2] = kolor/antykolor.
- ⟹ spin i kolor to **różne składniki jednej grupy, powiązane półproduktem** → **kandydat na lock §0**
  (rozstrzygnięcie należy do następnej strzałki, nie tej — nie przesądzam).

## 2. Uzasadnienie symetrii T (nie postulat)

T ⊂ SO(3) nie jest wybrane dowolnie: **wymuszone przez trialność-uniqueness** — spośród WSZYSTKICH
skończonych podgrup SO(3) tylko T daje π₁ z ℤ₃ w abelianizacji (2T), reszta gubi kolor (jak no-go
writhe→ℤ). Równoległość: filar why_n3 użył ℤ₂ (RP²=S²/ℤ₂) z symetrii substratu; tu T z tego samego
typu argumentu, tylko wymuszone przez wymóg koloru.

## 3. Karta wyników strzałki substrat→wstążki (wg RAMY §5)

| Reguła | Status |
|---|---|
| §5.1 metastabilny z topologii (nie hack) | ✅ klasa dysklinacji reperu = niezmiennik topologiczny |
| §5.2 składalny | ✅ iloczyn grupowy; nieabelowość → braiding fizyczny |
| §5.3 minimalny | ✅ **uniqueness** 2T; człon Skyrme'a USUNIĘTY |
| §5.4 background-free / samoreferencyjny | ✅ **DOMKNIĘTE** — reper samo-generowanej metryki, zero tła |
| §5.5 kolor+ładunek z jednej topologii | ⚠ kolor+spin z jednego reperu ✅; **ładunek** = treść NASTĘPNEJ strzałki; lock §0 OTWARTY |
| §5.6 spin z nawinięcia | ✅ π₃, centralne −1 (filar) |

## 4. Co pozostaje OTWARTE (uczciwie, poza zakresem tej strzałki)

- **Ładunek** i **lock kolor↔ładunek (§0)** — należą do strzałki wstążki→kwarki. 2T=Q₈⋊ℤ₃
  daje hook, ale lock nierozstrzygnięty.
- **Lokalizacja cząstki** (skończony rozmiar): pozostaje otwarta — ALE to **ten sam status co
  filar spinu** (B=1 to liczba kwantowa, nie z automatu stabilny soliton skończony w kanonicznym
  F-A, który jest runaway). Wstążka jako **nośnik liczb kwantowych** = background-free ✅; wstążka
  jako **zlokalizowana cząstka** = problem warstwy obiektowej, nie strukturalnej. Nie zamiatam:
  krok 2 (płaski Skyrme) był błędną próbą tego, zdyskwalifikowaną.

## 5. Jednozdaniowo

> Strzałka substrat→wstążki jest domknięta na poziomie STRUKTURALNYM i background-free:
> wstążki = dysklinacje reperu samo-generowanej metryki, target 2T wymuszony trialnością,
> niosące spin (π₃) i kolor (2T) z jednego reperu, bez Skyrme/Derricka/tła; ładunek, lock §0
> i lokalizacja cząstki pozostają OTWARTE jako treść wyższych warstw.

---

**Pliki:** `RIBBONS_s54_frame_structure.py`, `RIBBONS_group_compare.py`,
`RIBBONS_uniqueness_2T.py`, `RIBBONS_stepB_topology.py`
**Źródła rdzenia:** `core/sek08a_akcja_zunifikowana.tex`, `core/sek08c_metryka_z_substratu.tex`,
`research/qm_spin/README.md`
**Zdyskwalifikowane:** `RIBBONS_stepB2_skyrme_v3.py` (płaski Skyrme — nie-natywny, §5.4)
