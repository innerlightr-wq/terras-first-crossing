"""Independent audit of the 21/7 first-crossing certificate."""
from fractions import Fraction as F
from itertools import product
from math import gcd
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
K=1412737287532845639425; D=891337988684307508298
HEAD,TAIL="11111110","111100000"; MARKER="0111100000"
SPAN=3611678508; L=104
def lg(base,terms=240):
    x=F(base-1,base+1); lo=sum(2*x**i/i for i in range(1,terms,2))
    return lo, lo+2*x**(terms+1)/((terms+1)*(1-x*x))
lo2,hi2=lg(2); lo3,hi3=lg(3)
lam_lo,lam_hi=K*lo2-D*hi3, K*hi2-D*lo3
eta=F(8*D,K)-5
def a_of(t): return t*D//K+2
nb=(K-17)//8; end=K-len(TAIL)

print("="*78); print("0. the phase bound, which I nearly mis-flagged"); print("="*78)
print("  The certificate claims rho_t <= 3^(1+eta) < 256/81 before a six-one byte.")
print("  A decimal check looks like it FAILS (3^1.0474378 = 3.1605209 against")
print("  256/81 = 3.1604938).  It does not: exactly,")
print("     3^(1+eta) < 256/81  <=>  (8D/K - 4) ln3 < 8 ln2 - 4 ln3  <=>  D ln3 < K ln2,")
print("  i.e. precisely Lambda > 0.  So the bound is equivalent to the shell's own")
print("  defining inequality and holds with no slack to spare.")
chk("3^(1+eta) < 256/81 is EQUIVALENT to Lambda > 0", lam_lo>0 and (8*D-4*K)*hi3 < 8*K*lo2)
chk("and 256/81 - 2.999 < 0.162", F(256,81)-F("2.999")<F(".162"),
    "%.8f < 0.162"%float(F(256,81)-F("2.999")))
print("  LESSON: this one must be checked rationally; 7 decimal places is not enough.")

print(); print("="*78); print("1. the forced multiplicity"); print("="*78)
chk("0 < Lambda < 1e-18 and 0 < 20*eta < 1", 0<lam_lo<lam_hi<F(1,10**18) and 0<20*eta<1,
    "20*eta = %.6f"%float(20*eta))
chk("8 + 8*nbytes = end, a(8)=7, a(end)=D-4", 8+8*nb==end and a_of(8)==7 and a_of(end)==D-4)
sixes=D-11-5*nb
chk("six-one count = 8377183976278983657", sixes==8377183976278983657)
labels=7*21**12
chk("162*labels < sixes-19 <= 163*labels, so 163 starts share a 104-bit window",
    162*labels < sixes-19 <= 163*labels,
    "ratio = %.4f"%float(F(sixes-19,labels)))
print("  seven bytes behind and twelve ahead need 20*eta < 1 -- verified above.")

print(); print("="*78); print("2. the normalization, checked against my own derivation"); print("="*78)
print("  Certificate: z-residues are a translate of 2^-160 * sum_i 256^i 243^(6-i) C(b_i)")
print("  mod 3^35.  My independent derivation (concatenation C_uv = C_u 3^{d_v} + 2^{|u|} C_v)")
print("  gives per-byte scale 2^{8j - 8m - L} * 3^{5(m-1-j)} with m=7, L=104:")
print("     2^{-8*7-104} = 2^-160  and  2^{8j} = 256^j,  3^{5(6-j)} = 243^{6-j}.  IDENTICAL.")
chk("the two scale formulas agree symbolically for m=7, L=104", 8*7+104==160)
def carry(w): return sum(2**j*3**w[j+1:].count("1") for j in range(len(w)) if w[j]=="1")
A5=["".join(w)+"0" for w in product("01",repeat=7) if w.count("1")==5]
chk("21 five-one letters", len(A5)==21)
def cpp_carry(mask):
    c=0
    for j in range(7):
        if (mask>>j)&1: c=3*c+(1<<j)
    return c
chk("span.cpp's carry set equals mine",
    sorted(carry(w) for w in A5)==sorted(cpp_carry(m) for m in range(128) if bin(m).count("1")==5))
for m_,Lx in ((3,104),(4,104)):
    mod=3**(5*m_)
    mine={pow(2,-(8*m_+Lx),mod)*sum(pow(2,8*j,mod)*pow(3,5*(m_-1-j),mod)*carry(b)
          for j,b in enumerate(bs))%mod for bs in product(A5,repeat=m_)}
    readme={pow(2,-(8*m_+Lx),mod)*sum(pow(256,j,mod)*pow(243,m_-1-j,mod)*carry(b)
            for j,b in enumerate(bs))%mod for bs in product(A5,repeat=m_)}
    chk("depth %d: my formula and the README formula give the same residue set (%d)"
        %(m_,len(mine)), mine==readme)
print("  The depth-7 count 1372257936 and span 3611678508 were reproduced by compiling")
print("  and running span.cpp here (exit 0).  1.8e9 products cannot be held in a Python")
print("  set, so that pair is REPRODUCED rather than independently recomputed; the")
print("  normalization and carry convention above are what I verified independently.")

print(); print("="*78); print("3. the seed floor"); print("="*78)
floor=(2**L*SPAN - F(256,81)*F(D,3))/F(".162")
chk("floor > 4.5218e41", floor>F("4.5218")*10**41, "floor = %.6e"%float(floor))
chk("it matches the certificate's 452182363308404545602040205298328478303582",
    int(floor)==452182363308404545602040205298328478303582)
print("  derivation: the 163 anchor states lie in (2.999 n, (256/81)(n + D/3)], so their")
print("  span is < 0.162 n + (256/81) D/3; it is also >= 2^104 * sigma.  Solve for n.")
chk("floor * Lambda > D/3, which settles the ENDPOINT case with no failure hypothesis",
    floor*lam_lo>F(D,3), "floor*Lambda = %.4e vs D/3 = %.4e"%(float(floor*lam_lo),float(F(D,3))))
print("  (a failure at K needs n*Lambda < n(e^Lambda - 1) <= M_K <= D/3.)")

print(); print("="*78); print("4. earlier crossings: the h-range and the unique site"); print("="*78)
chk("-K*D/(3*floor*ln3) > -1 and K*Lambda/ln3 < 131, so -1 < h < 131",
    -K*F(D,3)/(floor*lo3)>-1 and K*lam_hi/lo3<131,
    "lower %.4f, upper %.4f"%(float(-K*F(D,3)/(floor*lo3)),float(K*lam_hi/lo3)))
chk("gcd(K,D)=1 so j = -h D^{-1} mod K is unique per h", gcd(K,D)==1)
inv=pow(D,-1,K); poss=[]
for h in range(131):
    j=(-h*inv)%K; t=j-2
    if not (8<=t<end and t%8==0): continue
    aa=a_of(t)
    if K*aa-D*j!=h or a_of(t+8)!=aa+5: continue
    elo,ehi=j*lo2-aa*hi3, j*hi2-aa*lo3
    if ehi<0: continue
    if floor*elo>F(aa,3): continue
    poss.append((h,j,t,aa,elo))
chk("the coarse memory test leaves exactly one site, h=60", len(poss)==1 and poss[0][0]==60,
    "%s"%[(p[0],p[1]) for p in poss])
chk("that site is j = 651838711708843936730", poss[0][1]==651838711708843936730)

print(); print("="*78); print("5. the refined memory bound and the closing comparison"); print("="*78)
def mcost(b): return sum(F(2**i,3**(b[:i].count("1")+1)) for i,bit in enumerate(b) if bit=="1")
A6=["".join(w)+"0" for w in product("01",repeat=7) if w.count("1")==6]
reg=[w for w in A5 if w!="00111110"]; bad="00111110"
chk("max regular cost 682/243, max six cost 1330/729, bad cost 844/243",
    max(map(mcost,reg))==F(682,243) and max(map(mcost,A6))==F(1330,729)
    and mcost(bad)==F(844,243))
chk("only 00111110 can cross early, and only at its SECOND zero",
    all(F(3**w[:j].count("1"),2**j)*F("2.99")>1 for w in reg+A6 for j in range(1,9))
    and all(F(3**bad[:j].count("1"),2**j)*F("2.99")>1 for j in (1,3,4,5,6,7,8))
    and F(3**bad[:2].count("1"),2**2)==F(1,4))
print("  so the exact safety condition for 00111110 is rho_start >= 4 (factor 1/4), which")
print("  is the condition I derived last round -- unchanged.")
chk("an EARLIER completed bad byte must start at rho>=4, giving cost/4 < 0.94",
    mcost(bad)/4<F(".94"), "%.6f < 0.94"%float(mcost(bad)/4))
chk("a regular byte gives cost/2.999 < 0.94 and a six-one byte likewise",
    max(map(mcost,reg))/F("2.999")<F(".94") and max(map(mcost,A6))/F("2.999")<F(".94"))
chk("the head contributes < 7/3", mcost(HEAD)<F(7,3))
h0,j0,t0,a0,elo0=poss[0]
cap=F(7,3)+F(".94")*((t0-8)//8)
chk("at the surviving site M_j/E < 0.952 * floor, hence strictly below the floor",
    cap/(elo0*floor)<F(".952"), "ratio = %.6f"%float(cap/(elo0*floor)))
chk("the two bounds are incompatible: 4.30250e41 < 4.52182e41",
    cap/elo0 < floor, "cap/E = %.6e < floor = %.6e"%(float(cap/elo0),float(floor)))

print(); print("="*78); print("6. marker uniqueness and scope"); print("="*78)
AL=A5+A6
chk("28 letters; MARKER absent from all %d three-byte products"%(len(AL)**3),
    len(AL)==28 and all(MARKER not in x+y+z for x,y,z in product(AL,repeat=3)))
chk("MARKER absent at every head junction",
    all(MARKER not in HEAD+x+y for x,y in product(AL,repeat=2)))
chk("MARKER occurs exactly once at every tail junction, at the very end",
    all((x+y+TAIL).find(MARKER)==len(x+y+TAIL)-len(MARKER) for x,y in product(AL,repeat=2)))
chk("the short marker 00000 is NOT unique here, so the long one is required",
    "00000" in "11111000"+"00111110")
print("  SCOPE.  The seed floor uses the WHOLE word (163 anchors spread over the middle")
print("  region), so the quantifier is over integers realizing a COMPLETE word.  It uses")
print("  no failure hypothesis and no minimality, so the conclusion is pointwise.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
