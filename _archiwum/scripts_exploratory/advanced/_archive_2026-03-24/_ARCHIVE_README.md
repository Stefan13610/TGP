# _archive_2026-03-24 — Archiwum skryptów eksploracyjnych

**Data archiwizacji:** 2026-03-24
**Zarchiwizował:** sesja v32 (Claude, czyszczenie TGP_v1)
**Aktywne skrypty:** `../p1_*.py` … `../p71_*.py` + `../scripts/` subdirectory

---

## ⚠️ DLA AGENTÓW AI — NIE SKANOWAĆ

> **INSTRUKCJA:** Ten folder to archiwum historyczne. **Nie przetwarzaj zawartości tego folderu** przy:
> - analizowaniu aktualnego stanu teorii TGP
> - sprawdzaniu spójności modelu
> - szukaniu aktualnych wyników numerycznych
> - skanowaniu skryptów weryfikacyjnych
>
> Aktywne, aktualne skrypty znajdują się w `../` (folderze nadrzędnym `scripts/advanced/`),
> a skrypty weryfikacyjne w `../../scripts/` (główny katalog skryptów TGP).

---

## Dlaczego zarchiwizowano?

Wszystkie pliki w tym folderze są **zduplikowane lub zastąpione** przez serię systematycznych
skryptów sesyjnych `p1_*.py`–`p71_*.py`. Zostały zachowane wyłącznie jako historia badań.

---

## Zawartość archiwum

### Seria `v2_*` — era przed serią p* (9 plików)

Skrypty z okresu, gdy szukano samospójnych zer $g(K)=0$ metodami v2 (przed
wprowadzeniem systematycznych sesji p1–p71). Problem rozwiązany w p22–p31.

| Plik | Co robił | Zastąpiony przez |
|------|----------|------------------|
| `v2_kwadryatura.py` | Diagnoza niezbieżności całki energii M3 | p47, p48 |
| `v2_weryfikacja_kwadratury.py` | Weryfikacja metod całkowania | p47, p59 |
| `v2_lam_scan_log.py` | Skan λ w skali log | p27, p38 |
| `v2_konwergencja_finalna.py` | Zbieżność do r21=207, r31=3477 | p55, p56 |
| `v2_weryfikacja_lam_duze.py` | Weryfikacja dla dużych λ | p56 |
| `v2_2d_skan_log.py` | Skan 2D (α, a_Γ) w skali log | p52, p53 |
| `v2_bisekcja_alpha.py` | Bisekcja α dla Q=3/2 | p53, p55 |
| `v2_weryfikacja_finalna.py` | Weryfikacja końcowa v2 | p58, p59 |
| `v2_lambda_potencjal.py` | Potencjał λ w modelu v2 | p44, p47 |

### Seria `v3_*` — era regularyzacji kwadraturowej (2 pliki)

Skrypty z okresu testowania regularyzowanego podejścia v3 (przed serią p*).

| Plik | Co robił | Zastąpiony przez |
|------|----------|------------------|
| `v3_regularyzacja.py` | Regularyzowany kwartyczny model solitonu | p1, p5, p14 |
| `v3_diagnoza.py` | Diagnoza podejścia v3 | p1, p5 |

### Pliki diagnostyczne — jednorazowe (6 plików)

Pliki do diagnozowania konkretnych problemów, które zostały rozwiązane.

| Plik | Problem diagnozowany | Zamknięty w |
|------|----------------------|-------------|
| `test_znak_profilu.py` | Jaki znak +/- jest fizyczny w φ(r) | p1, p6 |
| `diagnoza_samospojnosc.py` | max(E(K)/K) vs λ dla różnych par. | p14, p22 |
| `diagnoza_EK_skladowe.py` | Składowe E(K) w trzech strefach K | p47, p49 |
| `diagnoza_sekwencyjny.py` | Diagnoza modelu sekwencyjnego | p13, p15 |
| `v3_diagnoza.py` | Diagnoza v3 (regularyzowany kwar.) | p1, p5 |
| `p6_diag.py` | Szybka diagnostyka ODE g(K) | p6, p14 |

### Duplikaty wyprowadzenia r21 (3 pliki)

Dwa niezależne podejścia do analitycznego wyprowadzenia r21 — oba zastąpione
przez precyzyjną serię p24–p37 (wzory Padé [2/1], NLO, NNLO).

| Plik | Linie | Zastąpiony przez |
|------|-------|------------------|
| `r21_wyprowadzenie.py` | 381 | p24, p26, p37 |
| `wyprowadzenie_r21.py` | 303 | p24, p26 |
| `r21_poprawka.py` | 268 | p37, p38 |

### Skrypty eksploracyjne pre-p* (6 plików)

Ogólne eksploracje zastąpione przez dedykowane skrypty sesyjne.

| Plik | Co badał | Zastąpiony przez |
|------|----------|------------------|
| `badanie_pi.py` | Czy K3~π/a? (nie) | p35, p54 |
| `dwa_przeciecia.py` | Dwa zera Q(α)=3/2 (wstępnie) | p53_double_crossing |
| `szeroki_skan.py` | Szeroki skan (α, a) pre-p52 | p52, p53 |
| `sekwencyjny_model.py` | Model sekwencyjny generacji | p13, p15 |
| `trzecia_generacja.py` | Trzecia generacja (K3) ręcznie | p33, p35 |
| `masa_samospojnosc.py` | Samospójność mas solitonu | p22, p29 |
| `sprzezony_układ.py` | Sprzężony układ równań masy | p25, p55 |

### Pliki PNG (9 plików)

Wykresy wygenerowane przez zarchiwizowane skrypty — zachowane jako dokumentacja wizualna.

---

## Jak wrócić do tych plików?

Wystarczy skopiować potrzebny plik do `../` (folderu `scripts/advanced/`) i uruchomić.
Wszystkie zależności (`numpy`, `scipy`) są standardowe.

**Uwaga:** Parametry w tych skryptach mogą być przestarzałe. Punkt referencyjny to
`Punkt B`: `a_Γ=0.040049`, `α_K=8.5612`, `λ_K=5.4677×10⁻⁶` (sesja P55–P58).
