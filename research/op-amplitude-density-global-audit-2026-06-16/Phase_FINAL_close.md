---
title: "op-amplitude-density-global-audit — Phase FINAL (werdykt)"
date: 2026-06-16
type: research-cycle
status: CLOSED
cycle: op-amplitude-density-global-audit
phase: FINAL
verdict: INCONSISTENT-LABELING (nośna fizyka spójna; dewiant = moja edycja Opcji B #31)
---

# Phase FINAL — werdykt

## Werdykt (wyliczony z reguły LOCKED, NIE wybrany)

**INCONSISTENT** — reguła Phase 0: „≥1× G-INCONSISTENT ⟹ INCONSISTENT".
Znaleziono **4× G-INCONSISTENT**.

**Doprecyzowanie charakteru (kluczowe):**
- Nośna FIZYKA rdzenia jest **SPÓJNA** (11× G-CONSISTENT, 0 sprzeczności wewnątrz):
  jedno pole **Φ=gęstość** (sek01 `def:Phi`), **α=2 postulowane** na gęstości
  (selekcja C1–C3, `rem:alpha2-pivot-status`); metryka, soliton, masy, PPN, T_Hawking
  — wszystkie w tej samej zmiennej i wzajemnie spójne.
- Sprzeczność jest **ramowania/etykietowania** i jest **zlokalizowana w mojej edycji Opcji B
  z sesji #31** (3 uwagi) + 1 arytmetyczna zaległość po moim niedokończonym fixie g²→g⁴.

**To znaczy: cykl wykrył, że MOJA WŁASNA edycja #31 jest niespójna z nośnym rdzeniem.**
Wstępne założenie cyklu („zweryfikować, że Opcja B jest spójna w dół") zostało **obalone**
przez value-blind audyt. To jest dokładnie cel audytu odpowiedzialnościowego.

## Na czym polega błąd (precyzyjnie)

| | Co mówi nośny rdzeń | Co mówi moja Opcja B #31 |
|---|---|---|
| Pole kanoniczne | **gęstość Φ** (`def:Phi`) | amplituda φ |
| α=2 | **na gęstości Φ** (postulat C1–C3) | na amplitudzie; w gęstości α=½ |
| Status substratu | substrat NIE derywuje α=2 (`rem:alpha2-pivot-status`) | substrat→α=½ w gęstości ⟹ amplituda kanoniczna |

Technicznie poprawny fakt op-A3 — „substrat daje α=½ w gęstości; α nie jest niezmiennikiem
$\varphi\to\Phi\propto\varphi^2$" — został przeze mnie **błędnie przeramowany** na „⟹ kanoniczne
pole = amplituda". Poprawny wniosek: **substrat (α=½) ≠ postulat (α=2) na tej samej gęstości Φ
⟹ α=2 jest selekcją, nie derywacją** (= dokładnie `rem:alpha2-pivot-status`, już w rdzeniu).
Gęstość Φ pozostaje polem kanonicznym; amplituda ŝ=√Φ jest reprezentacją mikroskopową.

## Lista G-INCONSISTENT (do naprawy)

1. **sek08 `rem:amplitude-vs-density-alpha`** (l.1076–1099, MOJA #31) — przeramować:
   pole kanoniczne = **gęstość Φ**; α=2 = postulat **na gęstości**; substrat (amplituda ŝ=√Φ)
   daje α=½ ⟹ NIE derywuje α=2 (stąd selekcja). NIE „amplituda jest polem kanonicznym".
2. **sek10 `rem:K_to_f_amplitude`** l.205 (MOJA #31) — usunąć kolizję: w tej podsekcji
   `g≡ψ=Φ/Φ₀` (gęstość, l.142). Substratowa transformacja daje α=½ (eq:kinetic_macro);
   kanoniczny α=2 to **postulat `K(g)=g⁴` na gęstości g** (`rem:canonical_g4`), nie „g=amplituda".
3. **dodatekQ2 `rem:A3-correction-alpha`** (l.342–363, MOJA #31) — to samo przeramowanie.
4. **sek10:145** `K_sub(g)=K_geo g²` → `g^4` (dokończyć korektę #31; spójnie z box α=2
   i eq:Ksub_expansion_check, które już są g⁴).

## G-AMBIGUOUS (rekomendacja, nie bloker)
- **Overload symbolu φ**: sek08a `φ=Φ/Φ₀` (gęstość bezw.) vs sek10:104 `φ=φ_ref√(Φ/Φ₀)`
  (mikro-amplituda). Rekomendacja: w sek10 §K_to_f użyć osobnego symbolu dla mikro-amplitudy
  (np. $a$ lub $\hat s$) i zarezerwować φ/g dla gęstości bezwymiarowej. Niski priorytet
  (nie psuje fizyki), ale eliminuje główne źródło pomyłki, która zrodziła mój błąd #31.
- sek07 kosmologia: doprecyzować φ≡Φ/Φ₀ w pasażach predykcyjnych.

## Anti-Lakatos
- ✓ Werdykt **wyliczony** z reguły LOCKED (4× G-INCONSISTENT), NIE wybrany pod „spójne".
- ✓ **Zgłosiłem, że to MOJA edycja #31 jest niespójna** — nie zatajone, nie przeramowane na
  „rdzeń był zły". Nośny rdzeń jest poprawny; dewiantem jest mój dodatek.
- ✓ Każdy SUSPECT agentów **zweryfikowany firsthand** w źródłach (def:Phi, sek02, sek08a/b, sek10).
- ✓ NIE „naprawiłem" przez ciche przemianowanie — różnica fizyczna (α=2 vs ½ na gęstości)
  nazwana wprost i przypisana do postulatu vs substratu.
- ✓ Rdzeń NIE edytowany w tym cyklu — lista poprawek czeka na autoryzację (przeramowanie moich
  uwag to zmiana znacząca, nie czysta kosmetyka, więc pytam mimo że to moje własne edycje).

## Postęp vs stan przed cyklem
- Przed: Opcja B #31 zintegrowana, zweryfikowana TYLKO „kompiluje się" (552 s., 0 błędów).
- Po: wykryto, że Opcja B #31 jest **niespójna z nośną konwencją Φ=gęstość/α=2-postulat**;
  precyzyjna lista 4 poprawek + 1 rekomendacja anty-overload. Nośna fizyka potwierdzona spójna.

## Residual
- 4 poprawki rdzenia (przeramowanie 3 moich uwag + 1 arytmetyka) — **wymaga „działaj"**.
- Po naprawie: ponowny build (pdflatex) dla potwierdzenia 0 błędów.
