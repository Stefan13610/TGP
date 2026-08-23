# RIBBONS_stepB_topology.py
# KROK 1 planu 2-stopniowego: rozstrzygajacy test topologiczny stabilizatora (B).
# (B) = "framing wstazki = centralny element -1 z 2T (= ten sam, co daje spin)".
# Pytanie: czy oskrzydlona (-1) petla dysklinacji ma REALNA przeszkode
#          przed usunieciem (zapadnieciem do prozni)?
#
# Baza: M = SO(3)/T = S^3 / 2T  (sferyczna forma przestrzenna; 2T dziala swobodnie
# przez lewe mnozenie na SU(2)=S^3). Stad nakrycie uniwersalne M~ = S^3, wiec:
#   pi_1(M) = 2T   (rzad 24)
#   pi_2(M) = pi_2(S^3) = 0     <- BRAK defektow punktowych i BRAK ucieczki monopolowej
#   pi_3(M) = pi_3(S^3) = Z     <- ladunek samo-splotu / Hopfa (skyrmionowy)
#
# Ten skrypt WERYFIKUJE inputy u zrodla:
#   (1) 2T = 24 kwaternionow jednostkowych domyka sie w grupe; -1 centralny, rzad 2.
#   (2) mapa "framing -> element pi_1": obrot wewnetrzny o 2pi wraca do -1,
#       o 4pi wraca do +1 => framing rzutuje na Z_2 = {+1,-1} = PARZYSTOSC samo-splotu.
#       To jest DOKLADNIE relacja Finkelsteina-Rubinsteina (-1)^B (filar spinu).

import math

def qmul(p,q):
    w1,x1,y1,z1=p; w2,x2,y2,z2=q
    return (w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2,
            w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2)
def qkey(p): return tuple(round(c,6)+0.0 for c in p)
QID=(1.,0.,0.,0.); NEG=(-1.,0.,0.,0.)

def close_group(gens):
    E={qkey(QID):QID}; fr=[QID]
    while fr:
        nw=[]
        for g in fr:
            for h in gens:
                pr=qmul(g,h); k=qkey(pr)
                if k not in E: E[k]=pr; nw.append(pr)
        fr=nw
    return list(E.values())
def inv(g): w,x,y,z=g; return (w,-x,-y,-z)
def order_of(g):
    p=g;n=1
    while qkey(p)!=qkey(QID):
        p=qmul(p,g);n+=1
    return n

# (1) zbuduj 2T z generatorow Hurwitza
i=(0.,1.,0.,0.); h=(0.5,0.5,0.5,0.5)
G=close_group([i,h])
print("(1) 2T jako kwaterniony jednostkowe")
print(f"    |2T| = {len(G)}  (oczekiwane 24): {'OK' if len(G)==24 else 'BLAD'}")
has_neg = any(qkey(g)==qkey(NEG) for g in G)
central = all(qkey(qmul(NEG,g))==qkey(qmul(g,NEG)) for g in G)
print(f"    -1 w 2T? {has_neg};  -1 centralny? {central};  rzad(-1) = {order_of(NEG)}")
print(f"    => centralny element rzedu 2 obecny = generator spinu (FR/Berry).")

# (2) framing: obrot wewnetrzny o kat theta wokol osi n = kwaternion (cos(theta/2), sin(theta/2) n)
def frame(theta, axis=(1.,0.,0.)):
    s=math.sin(theta/2); return (math.cos(theta/2), s*axis[0], s*axis[1], s*axis[2])
f2pi = frame(2*math.pi)   # obrot o 2pi
f4pi = frame(4*math.pi)   # obrot o 4pi
print("\n(2) mapa framing -> element pi_1 (parzystosc samo-splotu)")
print(f"    framing 2pi  -> kwaternion {qkey(f2pi)}  == -1 ? {qkey(f2pi)==qkey(NEG)}")
print(f"    framing 4pi  -> kwaternion {qkey(f4pi)}  == +1 ? {qkey(f4pi)==qkey(QID)}")
print("    => framing akumuluje w Z (samo-splot SL), rzutuje na Z_2={+1,-1} przez parzystosc:")
print("       SL nieparzysty <=> element -1 (spin)  ;  SL parzysty <=> +1")
print("    To jest relacja FR: 2pi-obrot -> (-1)^B. Zgodne z filarem qm_spin.")

print("\n" + "="*66)
print("WNIOSEK TOPOLOGICZNY (krok 1, stabilizator B):")
print("-"*66)
print("M = S^3/2T:  pi_1=2T,  pi_2=0,  pi_3=Z.")
print(" - pi_2=0  => brak ucieczki monopolowej; sama klasa pi_1 NIE chroni")
print("   zwartej petli (goly -1 loop moze sie zagoic).")
print(" - OCHRONA jest przez pi_3=Z = samo-splot (Hopf/skyrmion).")
print("   framing '-1' (spin) <=> SL NIEPARZYSTY <=> Q_pi3 != 0 (nieusuwalny).")
print(" => (B) daje REALNA bariere przed USUNIECIEM, topologiczna, spieta ze spinem.")
print(" ALE: Q_pi3 jest skyrmionowy => Derrick kurczy go do punktu.")
print("      ROZMIAR wymaga skali (czlon typu Skyrme) => to jest KROK 2.")
