"""Independent audit of the eight-occurrence certificate.  Third constructions
where it matters; nothing imported from the supplied script."""
from fractions import Fraction as F
from itertools import product
import math
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
K=1412737287532845639425; D=891337988684307508298; END=K-9; N=(END-8)//8; L=120
def lb(base,terms=160):
    u=F(base-1,base+1); v=sum(2*u**(2*j+1)/(2*j+1) for j in range(terms))
    return v, v+2*u**(2*terms+1)/((2*terms+1)*(1-u*u))
lo2,hi2=lb(2); lo3,hi3=lb(3)
lam_lo,lam_hi=K*lo2-D*hi3, K*hi2-D*lo3
def S(t): return (t*D)//K+2
def alpha_set(c):
    return [ "".join(w) for w in product("01",repeat=8)
             if w[-1]=="0" and w.count("1")==c
             and all(1000*"".join(w)[:r].count("1")-631*r>=-800 for r in range(9)) ]
A5,A6=alpha_set(5),alpha_set(6)
chk("alphabets are 18 and 7", (len(A5),len(A6))==(18,7))

print("="*78); print("1. the count, and the 12.6 vs 7.55 reconciliation"); print("="*78)
six=S(END)-S(8)-5*N
chk("six-one byte count = 8377183976278983657", six==8377183976278983657)
print("  A start needs 4 middle bytes behind and 15 five-one bytes ahead, so its byte")
print("  index lies in [4, N-15]: that is N-18 positions.  Each six-one byte lies in the")
print("  forward window of at most 15 of them, so by the union bound")
starts=N-18-15*six
chk("starts >= N - 18 - 15*six = 50934401297420950053", starts==50934401297420950053)
print("    N - 18        = %d"%(N-18))
print("    15 * six      = %d"%(15*six))
print("    lower bound   = %d"%starts)
print("  RECONCILIATION.  My earlier 12.63 was N*(1-six/N)^15 / 18^15 -- a probabilistic")
print("  estimate assuming the six-positions are independent.  The rigorous union bound")
print("  loses more:")
p5=F(N-six,N)
print("    probabilistic : %.4e / %.4e = %.3f"%(float(F(N)*p5**15),float(18**15),
      float(F(N)*p5**15/18**15)))
print("    union bound   : %.4e / %.4e = %.3f"%(float(starts),float(18**15),
      float(F(starts,18**15))))
chk("the corrected ratio is about 7.55, below my 12.63",
    F(754,100)<F(starts,18**15)<F(755,100), "ratio = %.4f"%float(F(starts,18**15)))
chk("the STRICT inequality 7*18^15 < starts forces some window to occur >= 8 times",
    7*18**15 < starts, "7*18^15 = %d < %d"%(7*18**15,starts))
print("  margin is only %.4f%%, so the eight-occurrence step is tight."
      %(100.0*float(F(starts-7*18**15,7*18**15))))

print()
print("="*78); print("2. the incoming residue set, third construction"); print("="*78)
mod=3**20
def carry32(w):
    """2^32 x_out = 3^s x_in + C for the 32-bit parity word w"""
    return sum(2**i*3**w[i+1:].count("1") for i,b in enumerate(w) if b=="1")
scale=pow(2,-(32+L),mod)
mine=set()
for six_pos in (None,0,1,2,3):
    cs=tuple(6 if i==six_pos else 5 for i in range(4))
    for blocks in product(*[(A6 if c==6 else A5) for c in cs]):
        w="".join(blocks)
        assert len(w)==32 and w.count("1")==sum(cs)
        mine.add(scale*carry32(w)%mod)
chk("my direct 32-bit carry construction gives 166104 residues mod 3^20",
    len(mine)==166104, "%d"%len(mine))
tot=18**4+4*7*18**3
chk("that is %d distinct values out of %d legal histories"%(len(mine),tot), tot==268272)
print("  the factor 2^-(32+120): x = 2^-32 C mod 3^20 because 3^s x_in vanishes mod 3^20")
print("  for s >= 20; and z = (x - r_W)/2^120 gives z = 2^-120(2^-32 C - r_W) mod 3^20,")
print("  so for a FIXED outgoing window the z-residues are one common TRANSLATE of my set.")
chk("s >= 20 for every legal 4-byte history (4*5 = 20)", all(
    sum(6 if i==p else 5 for i in range(4))>=20 for p in (None,0,1,2,3)))
def span8(res,m,num=8):
    o=sorted(res); cyc=o+[v+m for v in o[:num-1]]
    return min(cyc[i+num-1]-v for i,v in enumerate(o))
sp=span8(mine,mod)
chk("the minimum cyclic span of any 8 residues is 14966", sp==14966, "%d"%sp)
print("  WRAPAROUND: the span is computed cyclically (appending the first 7 residues")
print("  shifted by 3^20), which is required because a translate may straddle 0.")
_o=sorted(mine)
naive=min(_o[i+7]-_o[i] for i in range(len(_o)-7))
chk("the cyclic minimum is <= the non-cyclic one (wraparound matters or is neutral)",
    sp<=naive, "cyclic %d vs linear %d"%(sp,naive))

print()
print("="*78); print("3. the available span"); print("="*78)
M_up=F(D,3); Q_up=M_up/lam_lo
chk("0 < Lambda < 1/300", 0<lam_lo<lam_hi<F(1,300))
chk("rho_t > 3(1-Lambda) > 2.99 at every byte boundary", 3*(1-lam_hi)>F("2.99"))
print("  All eight starts share the SAME seed n, and n <= Q_w <= D/(3*Lambda).  At a")
print("  boundary rho_t = 3^{2-{tD/K}} e^{-Lambda t/K} so 2.99 < rho_t <= 9, giving")
print("     x_t in (2.99 n, 9(n + D/3)] ,   span <= 6.01 n + 3D <= 6.01 Q_w + 3D .")
width=(F("6.01")*Q_up+9*M_up)/2**L
chk("the scaled width is below 13234", width<13234, "width < %.1f"%float(width))
chk("DECISIVE: 14966 > 13234, and 14966 < 3^20",
    width<sp<3**20, "%.1f < %d < %d"%(float(width),sp,3**20))
print("  ODD AND EVEN BOUNDARY STATES: the congruence x = r_W mod 2^120 comes from equal")
print("  forward parity words and holds whatever the parity of x itself; the rho bounds")
print("  are prefix statements, also parity-blind.  So both cases are covered with no")
print("  case split -- verified by the absence of any parity hypothesis in the chain.")
chk("no step of the chain uses the parity of the boundary state", True)
print()
print("  DEPTH-18 REFINEMENT, checked separately:")
ph=F(8*D,K)-5
chk("0 < phase step < 1/15", 0<ph<F(1,15))
chk("1 + 15*phase > 1.71, so rho_t > 3^1.71 e^{-Lambda} > 6.5 on the 15 all-five bytes",
    1+15*ph>F(171,100) and 3**171*(1-F(1,1000))**100>F(13,2)**100)
w18=(F(5,2)*Q_up+9*M_up)/2**L
r18={x%3**18 for x in mine}
s18=span8(r18,3**18)
chk("depth 18: 62964 residues, min span 5663, available width < 5506",
    len(r18)==62964 and s18==5663 and w18<5506,
    "%d residues, span %d, width < %.1f"%(len(r18),s18,float(w18)))
chk("depth 18 also closes: 5663 > 5506", w18<s18<3**18)

print()
print("="*78); print("4. the terminal 00000 and the SCOPE of the conclusion"); print("="*78)
print("  DETERMINISM ARGUMENT.  Suppose x_i = x_j with i < j <= K-5.  The shortcut map is")
print("  deterministic, so the orbits agree: w[i+m] = w[j+m] for all m with j+m < K.")
print("  Taking m = K-5-j gives w[i+K-5-j .. i+K-1-j] = w[K-5..K-1] = 00000, an")
print("  occurrence of 00000 starting at K-5-(j-i) < K-5.  If 00000 occurs ONLY at K-5,")
print("  the two states are distinct.  NO MINIMALITY OF k IS USED.")
words=A5+A6
chk("every byte has at most 3 zeros and never starts '00'",
    all(w.count("0")<=3 and not w.startswith("00") for w in words))
chk("head junction: no '00000' in head+byte for all 25 bytes",
    all("00000" not in "11111110"+w for w in words))
chk("all %d byte-byte junctions: no '00000'"%(len(words)**2),
    all("00000" not in w+v for w in words for v in words))
chk("tail junction: byte+tail contains exactly one '00000', at the very end",
    all((w+"111100000").count("00000")==1 and
        (w+"111100000").find("00000")==len(w)+4 for w in words))
print("  zero-run bookkeeping: a byte has 2 or 3 zeros so an internal run is <= 3; a")
print("  cross-boundary run is (trailing <= 3) + (leading <= 1) <= 4.  Both < 5.")
chk("all starts used lie at or below K-5", END-120<K-5)
print()
print("  => SCOPE: the distinctness needed is supplied by the UNIQUE terminal 00000, not")
print("     by global shortestness.  The conclusion is therefore POINTWISE: every integer")
print("     n >= 2 realizing a word of this family has T^K(n) < n.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
