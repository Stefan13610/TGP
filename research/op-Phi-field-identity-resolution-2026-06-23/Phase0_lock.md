---
title: "Phase 0 — op-Phi-field-identity-resolution: LOCK pytania i reguły decyzyjnej. ZAKRES SKORYGOWANY 2026-06-23: pytanie 'amplituda vs gęstość' JUŻ rozstrzygnięte (gęstość kanoniczna; op-amplitude-density-global-audit-2026-06-16 + #36). Cykl przekierowany na realnie OTWARTE pytanie: substrate-realizability α=2 na gęstości — czy istnieje dopuszczalny (skalarny, Z₂) bond, którego gęstość Φ=⟨ŝ²⟩ daje α=2, czy to no-go? Value-blind, anti-Lakatos."
type: phase_lock
status: PHASE0_LOCKED
phase: 0
cycle: op-Phi-field-identity-resolution-2026-06-23
created_date: 2026-06-23
authorization: "User 2026-06-23: 'tak działaj z op-Phi-field-identity-resolution'"
origin: "Rekomendacja eksperta TGP (sesja #38): najważniejsza rzecz do zbadania = fundament α=2 / tożsamość pola"
scope_correction: "Pierwotny zakres (Opcja B vs C, amplituda vs gęstość) JUŻ rozstrzygnięty: gęstość kanoniczna (op-amplitude-density-global-audit-2026-06-16, 11× G-CONSISTENT) + α=2 = aksjomatyczna selekcja na gęstości (#36 2026-06-22). Cykl przekierowany na otwarty residual: substrate-realizability α=2."
target: "sek10 eq:kinetic_macro / lem:K_phi2 / prop:substrate-action; sek08 rem:amplitude-vs-density-alpha, rem:alpha2-pivot-status-pl; FOUNDATIONS §1 (fundament nieruchomy)"
anti_lakatos_lock: ACTIVATED
---

# Phase 0 — LOCK (op-Phi-field-identity-resolution)

## §0 — Korekta zakresu (anti-Lakatos, jawna)

> **PIERWOTNE pytanie cyklu** (z rekomendacji #38): rozstrzygnąć op-A3 §5 Opcję B vs C —
> czy kanonicznym polem TGP jest amplituda û (Φ∝u) czy gęstość Φ=⟨ŝ²⟩ (Φ∝u²).
>
> **USTALENIE (firsthand, rdzeń):** to pytanie jest **JUŻ ROZSTRZYGNIĘTE** i zalockowane:
> - `op-amplitude-density-global-audit-2026-06-16` (value-blind, 11× G-CONSISTENT, werdykt
>   INCONSISTENT-LABELING): **polem kanonicznym jest gęstość Φ** (`def:Phi`); „Opcja B =
>   amplituda kanoniczna" była błędną edycją #31 i została **obalona**.
> - `op-alpha2-status-propagation-audit-2026-06-22` (#36): **α=2 = aksjomatyczna selekcja na
>   gęstości** (C1–C3), spójnie w README/papers/meta.
> - `sek08 rem:amplitude-vs-density-alpha` (l.1076–1094) + `sek10 rem:K_to_f_amplitude`
>   (l.203–220) potwierdzają stan LOCKED.
>
> **Re-litygowanie tego = forbidden move (degeneracja Lakatosa).** Cykl NIE powtarza
> rozstrzygnięcia. Zamiast tego adresuje residual, którego **żaden** z trzech cykli nie zamknął.

## §1 — Pytanie wiodące (skorygowane)

> **Czy istnieje DOPUSZCZALNY substrat (skalarny, Z₂-parzysty — per FOUNDATIONS §1 fundament
> nieruchomy), którego gęstość Φ=⟨ŝ²⟩ daje gradientowy współczynnik α=2 na gęstości — czy
> jest to NO-GO?** Jeśli realizowalny: czy bond jest selekcjonowany jakąś zasadą (⟹ α=2
> derywowalne), czy nie-kanoniczny/dostrojony (⟹ α=2 pozostaje aksjomatem, ale z jawnym
> wskazaniem JAKI substrat by go dał)?

**Dlaczego to jest realnie otwarte:** op-A3 (sympy 5/5) policzył α_density tylko dla dwóch
bondów: s=1 (K∝ŝ², α=0) i s=2 (K∝ŝ⁴, α=½). **Nie** przeszukał klasy bondów — nie zapytał,
czy *jakiś* admissible bond daje α=2. To decyduje, czy fundament „emergencja z jednego pola"
przeżywa w sektorze kinetycznym, czy α=2 jest nieredukowalnie doszyte ręcznie.

## §2 — Materiał dowodowy (z rdzenia, zlokalizowany firsthand)

| Źródło | Treść | Pole |
|---|---|---|
| **sek10 lem:K_phi2** (l.39–56) | bond $H_\Gamma=-J\sum(\hat s_i\hat s_j)^2$ → ekstrakcja gradientu → $K(\hat s)=K_{\rm geo}\hat s^2$ | mikro $\hat s$ |
| **sek10 prop:substrate-action** (l.130–134) | konwencja wagowa $K_{ij}=J(\hat s_i\hat s_j)^2$ na $\sum K_{ij}(\hat s_i-\hat s_j)^2$ → $K(\hat s)=K_{\rm geo}\hat s^4$ | mikro $\hat s$ |
| **sek10 eq:Phi_phi_id** (l.103–108) | $\Phi=(\hat s^2/\hat s_{\rm ref}^2)\Phi_0$ ⟹ $\Phi\propto\hat s^2$ (ontologia $\Phi=\langle\hat s^2\rangle$) | most micro↔macro |
| **sek10 eq:kinetic_macro** (l.114–125) | $K_{\rm geo}\hat s^4(\nabla\hat s)^2 = K_0\,\psi(\nabla\psi)^2$ (**liniowy** w ψ → α=½) | poprawna transformacja |
| **sek08 rem:amplitude-vs-density-alpha** (l.1076–1094) | LOCKED: pole kanon. = gęstość Φ; α=2 postulat na gęstości; substrat (K∝ŝ⁴) → α=½; α NIE niezmiennik $\hat s\to\Phi\propto\hat s^2$ | stan rdzenia |
| **FOUNDATIONS §1** | fundament nieruchomy: jedno pole **skalarne** Φ z symetrią **Z₂**; pivot bondu *w obrębie skalarnego Z₂* dozwolony | ograniczenie dopuszczalności |

## §3 — Konwencja i relacje LOCKED (weryfikujemy, NIE podważamy)

- **Relacja EL** (op-A3 F-α-A, POPRAWNA): dla $\mathcal L=-\tfrac12 K(\varphi)(\nabla\varphi)^2$,
  $K=C\varphi^{p}$ ⟹ współczynnik członu $(\nabla\varphi)^2/\varphi$ = $p/2$ = α.
- **Ontologia** (LOCKED): $\Phi=\langle\hat s^2\rangle$ ⟹ $\Phi\propto\hat s^2$; $\hat s=\hat s_{\rm ref}\sqrt\psi$.
- **Rodzina bondów** (parametryzacja): admissible Z₂ skalarny bond daje substratowe sprzężenie
  kinetyczne $K(\hat s)\propto\hat s^{2s}$, $s\in\mathbb{Z}_{\ge 1}$ (Z₂ ⟹ parzyste w $\hat s$).
  Mapowanie do rzędu bondu: ekstrakcja gradientu z $(\hat s_i\hat s_j)^n$ → $s=n-1$;
  konwencja wagowa $(\hat s_i\hat s_j)^m$ → $s=m$.

## §4 — Falsyfikatory / reguła decyzyjna LOCKED (value-blind)

| ID | Test | Reguła |
|---|---|---|
| **F1** | SANITY: $\alpha_{\rm density}(s)=(s-1)/2$ z poprawnej transformacji; reprodukuje op-A3 (s=1→0, s=2→½) | musi PASS; jeśli FAIL — błąd po mojej stronie |
| **F2** | Rozwiązać $\alpha_{\rm density}=2$ dla $s$ | sympy zwraca $s^\star$ value-blind |
| **F3** | Czy $s^\star$ odpowiada **dopuszczalnemu** bondowi (skalarny, Z₂-parzysty, całkowity rząd)? | TAK → REALIZOWALNE; NIE (np. $s$ niecałkowite/ujemne) → NO-GO |
| **F4** | Czy $s^\star$ = bond **kanoniczny v2** ($(\hat s_i\hat s_j)^2$, s=2)? | TAK → α=2 DERYWOWALNE z kanonu; NIE → nie-kanoniczny (patrz F5) |
| **F5** | Cross-consistency: czy bond dający α=2 współistnieje z kanonicznym potencjałem $V\sim\Phi^3,\Phi^4$ (β=γ)? | bada, czy α=2 i V wymagają TEGO SAMEGO rzędu bondu, czy różnych członów |

**Reguła agregatu (werdykt WYLICZONY, nie wybrany):**
- **NO-GO** — jeśli żadne dopuszczalne (całkowite $s\ge1$) nie daje α=2 ⟹ α=2 nieredukowalnie aksjomatyczne na gęstości; „emergencja z substratu" dowodliwie niepełna w sektorze kinetycznym.
- **REALIZABLE-CANONICAL** — jeśli kanoniczny bond v2 (s=2) daje α=2 ⟹ α=2 derywowalne (sprzeczne z op-A3; oznaczałoby błąd — wymaga rewizji).
- **REALIZABLE-NONCANONICAL** — jeśli α=2 wymaga dopuszczalnego, ale **nie-kanonicznego** bondu (s≠2), nieselekcjonowanego obecną zasadą ⟹ α=2 pozostaje aksjomatem na gęstości, ALE konstruktywnie: wskazany konkretny substrat, który by go dał (+ ocena dostrojenia / współistnienia z V).

## §5 — Forbidden moves (anti-Lakatos)

1. **Zakaz re-litygowania** „amplituda vs gęstość" (rozstrzygnięte — §0). Gęstość kanoniczna LOCKED.
2. **Zakaz przesądzania** wyniku F2–F5 — sympy liczy value-blind.
3. **Zakaz rewizji rdzenia** — werdykt = rekomendacja, NIE edycja core (osobna autoryzacja).
4. **Relacja EL** ($\alpha=p/2$) rdzeniowa i POPRAWNA — nie podważam.
5. **Obie strony** jawnie: za realizowalnością (admissible high-order bond istnieje) i przeciw
   (nie-kanoniczny / dostrojony / niezgodny z V).
6. Budżet nowych stałych = 0. Jeśli REALIZABLE-NONCANONICAL — podać konstruktywnie, który bond.

## §6 — Oczekiwanie (INFORMATIONAL, nie wiąże)

Z $\alpha_{\rm density}=(s-1)/2$: α=2 ⟹ $s=5$ ⟹ $K(\hat s)\propto\hat s^{10}$ (bond rzędu n=6
w ekstrakcji / m=5 w konwencji wagowej). To **dopuszczalne** (skalar + Z₂-parzysty), ale **nie**
kanoniczny bond v2 (s=2). Prawdopodobny werdykt: **REALIZABLE-NONCANONICAL** — α=2 nie jest
no-go, lecz wymaga wysokorzędowego bondu nieselekcjonowanego obecną zasadą. Werdykt rozstrzygnie sympy.

## §7 — Status
**🔒 PHASE 0 LOCKED.** Zakres skorygowany (§0). Anti-Lakatos aktywny. Przejście do derywacji sympy (value-blind).
