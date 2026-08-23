# Szczebel 0 wieży: σ_ab NIE MA niezerowej próżni ⟹ RP² nie ma domu w substracie

**Data:** 2026-08-10
**Typ:** KONTROLA SKŁADNIKÓW przed budową — **wynik NEGATYWNY**
**Skrypt:** `TOWER_r0_sigma_orders.py` (3D Ising Metropolis, L=6..12, + kontrola pozytywna)
**Powód:** przed stawianiem wieży trzeba wiedzieć, czy σ_ab w ogóle ma rozmaitość próżni.

---

## 0. Wynik

> **σ_ab = ⟨ŝ ŝ₊â⟩^TF znika w jednorodnej próżni substratu Isinga — nie numerycznie,
> lecz TOŻSAMOŚCIOWO z symetrii.** Nie jest więc niezależnym parametrem porządku
> i **nie ma własnej rozmaitości próżni**. ⟹ **RP² nie jest rozmaitością próżni σ_ab.**
> Moje twierdzenie z `ANALIZA_status-sigma_ROZSTRZYGNIETE` („filar spinu ma nośnik w σ_ab")
> jest **BŁĘDNE i WYCOFANE**.

## 1. Argument ścisły (nie wymaga numeryki)

Grupa symetrii sieci kubicznej działa **przechodnio na trzy osie**. H_Γ jest wobec niej
niezmienniczy. Czyste stany Gibbsa fazy uporządkowanej Isinga łamią **wyłącznie ℤ₂ wewnętrzną**
(ŝ→−ŝ); **oba** stany (±) pozostają **kubicznie symetryczne**. Zatem exact:

```
C_x = ⟨ŝ_i ŝ_{i+x̂}⟩ = C_y = C_z  ⟹  M = diag(C,C,C)  ⟹  σ = M − (tr M/3)·I ≡ 0
```

**Ising łamie symetrię WEWNĘTRZNĄ, nie łamie OBROTÓW.** Nematyk wymaga złamania obrotów.
Substrat tego nie robi.

## 2. Potwierdzenie numeryczne

**T = T_c** (najostrzejszy test — tam fluktuacje są największe):

| L | N | C_x | C_y | C_z | \|σ\| |
|---|---|---|---|---|---|
| 6 | 216 | 0.378170 | 0.380170 | 0.379111 | 1.02e−03 |
| 8 | 512 | 0.367491 | 0.370522 | 0.369434 | 1.66e−03 |
| 10 | 1000 | 0.355733 | 0.355091 | 0.356312 | 6.21e−04 |
| 12 | 1728 | 0.349981 | 0.349932 | 0.350351 | 2.63e−04 |

**|σ| ~ L^(−2.03)** — szybciej niż L^(−1.5) oczekiwane dla czystej fluktuacji.
**Werdykt: FLUKTUACJA, σ → 0.**

**T = 0.7 T_c:** wartości 3.0e−04, 1.7e−05, 5.2e−05, 2.0e−04 — **niemonotoniczne, czysty szum**
przy C ≈ 0.875. Anizotropia względna ~10⁻⁴. Dopasowany wykładnik (−0.45) **nie ma sensu**
— dopasowuję szum. **Nie raportuję tego jako „nierozstrzygnięte" na korzyść porządku.**

## 3. Uczciwe ograniczenia tego testu

- **Kontrola pozytywna jest SŁABA.** Przy wymuszonym J_z=1.30 dostałem |σ| ≈ 1.5e−03 —
  tylko ~10× powyżej szumu, bo głęboko w fazie uporządkowanej **wszystkie korelatory
  nasycają się** przy ~0.92. Test wykrywa porządek, ale z małym marginesem.
- **Testowałem H_Γ w postaci v1** (bilinearnej, eq:B-H). Wersja **v2 (gradientowa)**
  sprawy nie ratuje — wręcz przeciwnie: wg `rem:B-ssb-v2` człon wiązaniowy v2
  „**z definicji znika na każdej konfiguracji jednorodnej**". Czyli w v2 próżnia jest
  **jeszcze bardziej** izotropowa. Wniosek zachodzi *a fortiori*.
- Argument §1 jest **ścisły**; numeryka jest tylko potwierdzeniem.

## 4. Co to znaczy — dokładnie, bez rozszerzania

σ_ab **istnieje** jako obserwabla poziomu 0 (to zostaje), ale jest **niezerowa wyłącznie
tam, gdzie jest niejednorodność** — czyli jest **niewolnikiem gradientów**, dokładnie jak
σ^ij = (∂Φ)(∂Φ)−⅓δ(∇Φ)². Moje rozróżnienie „σ_ab niezależny vs σ^ij pochodny" **traci ostrze**:
oba znikają w próżni, oba żyją tylko na gradientach.

**Konsekwencja:** obiekt bez własnej próżni **nie ma rozmaitości próżni**, więc
nie klasyfikuje defektów przez π_n. Ani RP², ani SO(3)/D₂.

## 5. Status filaru spinu — OTWARTY, NIE OBALONY

**Nie powtarzam błędu z `ANALIZA_czy-TGP-ma-defekty`.** Padła **jedna droga** do RP²
(przez σ_ab), a nie filar.

| Droga do spinu-½ | Status |
|---|---|
| RP² jako rozmaitość próżni σ_ab | ⛔ **martwa** (ten dokument) |
| `why_n3` Phase 3: RP² + faza Berry'ego π, pole **Δ_a** | ❓ **NIEZBADANE** — nie wiem, czym jest Δ_a i czy ma nośnik |
| `qm_spin`: B=1 z π₃(S³) + Finkelstein-Rubinstein na vielbeinie | ❓ **NIEZBADANE** po tym wyniku |

**Nie twierdzę, że filar upadł.** Twierdzę, że **nośnik, który mu przypisałem, nie istnieje**,
a prawdziwy nośnik jest **niesprawdzony**. To musi być następny check — u źródła, przed czymkolwiek.

## 6. Rachunek sumienia

To **piąty** raz w tym programie, gdy zbudowałem warstwę na składniku,
którego istnienia nie sprawdziłem (płaski Skyrme → cap budżetowy → limit tempa →
„brak defektów" → **nośnik σ_ab**). Różnica: tym razem złapałem to **sam, przed budową wieży**,
a nie po audycie. Ale wzorzec się nie zmienił i **nie zamierzam go zamiatać**:
**przypisałem nośnik na podstawie nazwy obiektu („korelator kierunkowy"), nie na podstawie
jego wartości próżniowej.**

---

**Rachunek:** `TOWER_r0_sigma_orders.py`
**Źródła:** `axioms/substrat/dodatekB_substrat.tex` eq:B-H (lin. 51–58), rem:B-ssb-v2 (lin. 83–88)
