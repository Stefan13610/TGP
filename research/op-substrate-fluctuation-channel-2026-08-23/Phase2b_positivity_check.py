# -*- coding: utf-8 -*-
"""
Phase 2b (addendum diagnostyczny) — op-substrate-fluctuation-channel
Kontekst: Phase 2 QF-4c FAIL-as-implemented: G(d=22; L=64, m=2.0) = -2.33e-18
z FFT double. Diagnoza: prawdziwa wartosc ~ e^{-mu*22} = e^{-38.8} ~ 1e-17 —
na/pod progiem zaokraglen FFT (eps * G(0) ~ 2.4e-17). Test porownywal szum
maszynowy do zera — blad implementacji TESTU, nie maszynerii.

Weryfikacja niezalezna (ten skrypt):
  (a) prog zaokraglen: |G_fft(22)| < eps * G(0)  [wartosc = szum]
  (b) suma dokladna mpmath (dps=40) tego samego wyrazenia spektralnego:
      G(22,0,0) = (1/L^3) sum_k cos(k1*22) / (m^2 + sum 2(1-cos k_a)) > 0
  (c) sanity metody: mpmath vs FFT dla d=4 (obie wartosci >> progu) — zgodnosc
Argument analityczny (do Phase_FINAL): reprezentacja Neumanna/bladzenia
losowego z zabijaniem — wszystkie wyrazy dodatnie => G_m > 0 to twierdzenie;
numerycznie weryfikowalne tylko do precyzji arytmetyki.
"""
import numpy as np
from mpmath import mp, mpf, cos, pi

mp.dps = 40
L, m, d_bad, d_ctrl = 64, 2.0, 22, 4

# --- FFT (jak w Phase 2) ---
k = 2.0 * np.pi * np.fft.fftfreq(L)
cx = np.cos(k)
om = 2.0 * ((1 - cx)[:, None, None] + (1 - cx)[None, :, None]
            + (1 - cx)[None, None, :])
G = np.real(np.fft.ifftn(1.0 / (m * m + om)))
G0, Gbad_fft, Gctrl_fft = G[0, 0, 0], G[d_bad, 0, 0], G[d_ctrl, 0, 0]
floor = np.finfo(float).eps * G0
print(f"G_fft(0) = {G0:.6f}; prog zaokraglen eps*G(0) = {floor:.3e}")
print(f"G_fft({d_bad}) = {Gbad_fft:.3e}  |.| < prog: {abs(Gbad_fft) < floor}")

# --- mpmath: suma spektralna dokladna ---
c = [cos(2 * pi * n / L) for n in range(L)]
w1 = [2 * (1 - cv) for cv in c]
m2 = mpf(m) ** 2
# lista w23 dla wszystkich (n2,n3):
w23 = [w1[n2] + w1[n3] for n2 in range(L) for n3 in range(L)]

def G_axis_mp(d):
    tot = mpf(0)
    for n1 in range(L):
        num = cos(2 * pi * n1 * d / L)
        base = m2 + w1[n1]
        s = mpf(0)
        for w in w23:
            s += 1 / (base + w)
        tot += num * s
    return tot / L**3

Gbad_mp = G_axis_mp(d_bad)
Gctrl_mp = G_axis_mp(d_ctrl)
print(f"G_mp({d_bad})  = {mp.nstr(Gbad_mp, 12)}   (> 0: {Gbad_mp > 0})")
print(f"G_mp({d_ctrl})   = {mp.nstr(Gctrl_mp, 12)} vs FFT {Gctrl_fft:.12e}")
rel = abs(float(Gctrl_mp) - Gctrl_fft) / Gctrl_fft
print(f"sanity d={d_ctrl}: |mp-fft|/fft = {rel:.3e}")

ok_a = abs(Gbad_fft) < floor
ok_b = Gbad_mp > 0
ok_c = rel < 1e-12
print("-" * 60)
print(f"[{'PASS' if ok_a else 'FAIL'}] (a) wartosc FFT pod progiem zaokraglen")
print(f"[{'PASS' if ok_b else 'FAIL'}] (b) G({d_bad}) > 0 w dps=40: "
      f"{mp.nstr(Gbad_mp, 8)}")
print(f"[{'PASS' if ok_c else 'FAIL'}] (c) sanity metody mpmath (d={d_ctrl})")
if ok_a and ok_b and ok_c:
    print("WNIOSEK 2b: QF-4c pad = artefakt precyzji testu; pozytywnosc G_m")
    print("potwierdzona niezaleznie (dps=40). Kryterium LOCK (G_m>0) STOI;")
    print("werdykt kontroli QF-4c: PASS po udokumentowanej korekcie")
    print("implementacyjnej testu (wzorzec A3: blad implementacji ZNALEZIONY).")
else:
    print("WNIOSEK 2b: diagnoza NIEROZSTRZYGNIETA — QF pozostaje INCONCLUSIVE.")
    raise SystemExit(1)
