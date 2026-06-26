---
title: "op-amplitude-density-global-audit — Phase 2 (klasyfikacja value-blind)"
date: 2026-06-16
type: research-cycle
status: DONE
cycle: op-amplitude-density-global-audit
phase: 2
---

# Phase 2 — klasyfikacja użyć wg gate'ów (value-blind)

Metoda: 3 równoległe Explore-agenty zinwentaryzowały użycia $\Phi/\varphi/\langle\hat s^2\rangle$;
**każdy SUSPECT zweryfikowałem firsthand w źródłach** (agenty są zwiadem, nie orzecznikiem).
Klasyfikacja wg gate'ów z `Phase0_lock.md`.

## Ustalenie fundamentalne (linchpin)

**Definicja źródłowa Φ (sek01 `def:Phi`, sek01:89, sek01:259):**
> „Pole $\Phi$ ... mierzy **gęstość wygenerowanej przestrzeni**"; $\Phi\propto\langle\hat s^2\rangle_{\mathcal B(x)}$;
> $\Phi=\langle\hat s^2\rangle_{\rm cg}$.

⟹ **Kanonicznym polem TGP jest GĘSTOŚĆ Φ** (nie amplituda). Cały nośny aparat makroskopowy
stosuje **α=2 do tej gęstości**, jako **postulat/selekcję** (C1–C3, `rem:alpha2-pivot-status`),
NIE jako derywację z substratu.

To **odwraca** wstępne założenie cyklu („zweryfikować, że Opcja B jest spójna w dół").
Audyt pokazał: nośny rdzeń jest spójny w konwencji **Φ=gęstość, α=2 postulowane**;
to **moje uwagi Opcji B (#31)** wprowadziły konkurencyjne ramowanie („pole kanoniczne = amplituda").

## Tabela A — NOŚNY RDZEŃ (konwencja: Φ=gęstość, α=2 postulowane)

| Lokalizacja | Kontekst | Zmienna | Fragment | Gate |
|---|---|---|---|---|
| sek01 `def:Phi` (l.8–25) | definicja pola | Φ = gęstość | „Φ mierzy gęstość wygenerowanej przestrzeni" | **G-CONSISTENT** (fundament) |
| sek01:89, :259 | mikro-most | Φ∝⟨ŝ²⟩ | `\Phi=\langle\hat s^2\rangle_{\rm cg}` | **G-CONSISTENT** |
| sek01:1261, :1316 | równanie pola | Φ (gęstość) | `\nabla^2\Phi + 2(\nabla\Phi)^2/\Phi` (α=2) | **G-CONSISTENT** |
| sek02:88 `eq:N-preview` | samointerferencja | Φ (gęstość) | `\mathcal N[\Phi]=\frac{\alpha}{\Phi_0}\frac{(\nabla\Phi)^2}{\Phi}+\dots,\ \alpha=2` | **G-CONSISTENT** |
| sek08a:167–171 | człon kinetyczny | φ≡Φ/Φ₀ (gęstość bezw.) | `K(\varphi)=K_{\rm geo}\varphi^4\ (\alpha=2)` | **G-CONSISTENT** |
| sek08b:206,240,244 | ODE solitonu | g≡Φ/Φ₀ (gęstość) | `K_{\rm geo}[g^2g''+2g(g')^2+\tfrac2r g^2 g']=\nabla\!\cdot\!(g^2\nabla g)` → α=2 | **G-CONSISTENT** |
| sek08c:408–432 | metryka g_eff | ψ≡Φ/Φ₀ (gęstość) | `h(\psi)=\psi/(4-3\psi)`, `\sqrt{-g_{\rm eff}}=c_0\psi/(4-3\psi)` | **G-CONSISTENT** (metryka←gęstość; gęstość = pole kanoniczne) |
| sek00:390 (kol. prawa) | tab. ODE | g (gęstość) | `g''+\tfrac2r g'+\tfrac{2}{g}(g')^2=V'` (α=2) | **G-CONSISTENT** |
| sek06:209,430 | T_Hawking, r_H | Φ_H (gęstość lokalna) | `T_H=\hbar(\Phi_H)\kappa/(2\pi c_0 k_B)` | **G-CONSISTENT** |
| sek09:172,388 | gęstość z amplitudy | Φ=⟨φ_i²⟩ | `\Phi=\langle\hat\phi_i^2\rangle` | **G-CONSISTENT** |
| dodatekQ:37 | coarse-graining | Φ=⟨φ_i²⟩ | `\Phi(x)=\frac1{N_B}\sum_{i\in B}\phi_i^2` | **G-CONSISTENT** |

**Wniosek tabeli A:** nośny rdzeń **wzajemnie spójny** — jedno pole Φ=gęstość, α=2 postulowane,
metryka/soliton/masy/PPN wszystkie w tej samej zmiennej. **0 sprzeczności wewnątrz.**

## Tabela B — DEWIANT (moje edycje Opcji B #31 + zaległości)

| Lokalizacja | Pochodzenie | Problem | Gate |
|---|---|---|---|
| **sek08 `rem:amplitude-vs-density-alpha`** (l.1076–1099) | **MOJA EDYCJA #31** | Twierdzi: „kanoniczne pole kinetyczne = **amplituda** φ; α=2 w amplitudzie; gęstość Φ = **osobna, NIE pole kinetyczne**, α=½". SPRZECZNE z nośnym rdzeniem, gdzie α=2 JEST na gęstości Φ, a Φ JEST polem kanonicznym. | **G-INCONSISTENT** |
| **sek10 `rem:K_to_f_amplitude`** l.205 | **MOJA EDYCJA #31** | „g w tej podsekcji to **amplituda**" — ale l.142 definiuje `g≡ψ=Φ/Φ₀` (**gęstość**). Kolizja symbolu g WEWNĄTRZ jednej podsekcji. | **G-INCONSISTENT** |
| **dodatekQ2 `rem:A3-correction-alpha`** (l.342–363) | **MOJA EDYCJA #31** | To samo dewiantowe ramowanie („α=2 tylko w amplitudzie; gęstość→α=½"). | **G-INCONSISTENT** (ramowanie) |
| **sek10:145** `K_sub(g)=K_geo g²` | pre-existing + **niedokończona moja korekta #31** g²→g⁴ | Box (l.149) twierdzi α=2 (⟹ f=1+4ln g ⟹ K=g⁴), eq:Ksub_expansion_check skorygowane do g⁴, ale l.145 wciąż `g²` (⟹ α=1). Arytmetyczna zaległość. | **G-INCONSISTENT** (drobne) |
| **symbol φ overload** | pre-existing | sek08a: `φ=Φ/Φ₀` (gęstość bezwymiarowa). sek10:104 `eq:Phi_phi_id`: `φ=φ_ref√(Φ/Φ₀)` (mikro-amplituda ∝√Φ). To samo „K(φ)=φ⁴" ⟹ α=2 (sek08a) vs α=1 (sek10). | **G-AMBIGUOUS** (przeciążony symbol) |
| sek07:103–135 | pre-existing | kosmologia miesza φ/Φ bez deklaracji roli | **G-AMBIGUOUS** (niski wpływ) |
| sek06:19 `Φ→∞` | pre-existing | nieznormalizowane Φ (vs ψ=Φ/Φ₀ gdzie indziej) | **G-NEUTRAL** (kosmetyka) |
| mass_scaling_k4 R5 vs Phase2 | pre-existing | dot. α=1 (substrat) vs α=2 (kanon) w teście mas — **inny problem**, nie amplituda/gęstość | **G-NEUTRAL** (poza zakresem) |

## Sedno sprzeczności (precyzyjnie)

Substrat (op-A3, sympy): mikroskopowe `K∝(amplituda ŝ)⁴`. Pod identyfikacją
`Φ=⟨ŝ²⟩∝ŝ²` pełny człon kinetyczny przepisany w gęstości (sek10 `eq:kinetic_macro`,
z jakobianem `(∇ŝ)²=(∇Φ)²/(4Φ)`) daje `K_eff(Φ)∝Φ`, czyli **α=½ w gęstości**.
Rdzeń **postuluje** zaś `K(Φ)=Φ⁴` ⟹ **α=2 w gęstości** (selekcja C1–C3).
**To dwie różne rzeczy i rdzeń o tym wie** (`rem:alpha2-pivot-status`: α=2 NIE jest derywacją z substratu).

**Mój błąd #31:** z technicznie poprawnego faktu („substrat→α=½ w gęstości; α nie jest
niezmiennikiem φ→√Φ") wyciągnąłem błędne ramowanie: „⟹ kanoniczne pole = amplituda, α=2 żyje
w amplitudzie". Tymczasem kanoniczne pole TGP **jest gęstością Φ**, a α=2 jest **postulatem na
gęstości** — substrat go nie wyprowadza, i to właśnie jest treść `rem:alpha2-pivot-status`.
Poprawne stwierdzenie: substrat (α=½) ≠ postulat (α=2) na gęstości ⟹ α=2 to selekcja, NIE derywacja.

## Agregat gate'ów
- G-CONSISTENT: 11 (cały nośny rdzeń, wzajemnie spójny)
- **G-INCONSISTENT: 4** (3× moje uwagi Opcji B #31 [sek08, sek10, dodatekQ2] + 1× arytmetyczna zaległość sek10:145)
- G-AMBIGUOUS: 2 (overload φ; sek07 kosmologia)
- G-NEUTRAL: 2

**≥1 G-INCONSISTENT ⟹ reguła LOCKED daje werdykt INCONSISTENT** (zob. Phase_FINAL).
