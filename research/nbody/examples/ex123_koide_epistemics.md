# ex123 — Epistemika: Ścieżka 9 vs relacja Koide

**Cel:** rozdzielić jawnie dwa **różne** sposoby domknięcia trzeciej generacji przy ustalonym \(r_{21}\):

1. **Reguła geometryczna na \(g_0\)** (jak w ex106 T10):  
   \(r_{31}^{(\varphi^2)} = \bigl(A_{\mathrm{tail}}(\varphi^2 g_0^*)/A_{\mathrm{tail}}(g_0^*)\bigr)^4\).

2. **Relacja Koide** (algebra mas): przy tym samym \(r_{21}\) z Ścieżki 9,  
   \(r_{31}^{(\mathrm{K})}\) jest **jednoznacznie** wyznaczone przez warunek \(Q_K=3/2\) (funkcja `koide_r31_from_r21` w kodzie).

Te dwie liczby **nie muszą** być równe: zgodność z PDG w jednym podejściu nie implikuje zgodności w drugim bez dodatkowego twierdzenia łączącego \(\varphi^2\) z Koide.

**Wyjście numeryczne:** `ex123_koide_epistemics.csv` + konsola.

**Powiązanie z planem:** `PLAN_TAIL_DOMKNIECIE_ROZWOJ.md` — Faza C, punkt o niezależności od Koida.
