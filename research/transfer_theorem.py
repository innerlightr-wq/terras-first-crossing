"""Parameterized transfer theorem: certify all admissible shells at once."""
from fractions import Fraction as F
from decimal import Decimal, getcontext
from itertools import product
import math
getcontext().prec=300
AL=Decimal(3).ln()/Decimal(2).ln()
ALF=F(int(AL*Decimal(10)**260),10**260)
def lg(base,terms=240):
    x=F(base-1,base+1); lo=sum(2*x**i/i for i in range(1,terms,2))
    return lo, lo+2*x**(terms+1)/((terms+1)*(1-x*x))
lo2,hi2=lg(2); lo3,hi3=lg(3)
SIGMA=3611678508; L=104; PAST=7; FW=12; labels=7*21**FW
FLOOR=(2**L*SIGMA-F(256,81)*F(10**21,3))/F(".162")   # placeholder, recomputed per shell
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)

print("="*78)
print("1.  the early-crossing case is SHELL-INDEPENDENT and uniformly vacuous")
print("="*78)
print("  A failure at an early crossing j needs  floor * E <= M_j <= cap(t),  with")
print("     E = j ln2 - a ln3 = ln2 (j - a*alpha),   cap(t) = 7/3 + 0.94 (t-8)/8 .")
print("  So the site is DANGEROUS only if  ln2 (j - a alpha) <= cap(t)/floor.")
print("  Neither side mentions (K,D) except through the floor, whose dominant term")
print("  2^104 * sigma_163 is itself shell-independent.  Hence a single Diophantine")
print("  condition on (j,a) decides every shell at once.")
# floor's shell dependence is only the tiny (256/81)D/3 term; bound it below
Dmax=10**24
floor_lo=(2**L*SIGMA-F(256,81)*F(Dmax,3))/F(".162")
chk("the floor is >= 4.52182e41 for every shell with D <= 1e24",
    floor_lo>F("4.5218")*10**41, "floor >= %.6e"%float(floor_lo))
print("  A dangerous site needs  E = ln2 (j - a alpha) <= cap(t)/floor, with")
print("     cap(t) = 7/3 + 0.94 (t-8)/8  and  t = j-2 byte-aligned.")
print("  Enumerate the UPPER approximants (a, j = ceil(a alpha)) of alpha directly --")
print("  only these can make E small -- and test each against its own cap/floor.")
x2,cf2=AL,[]
for _ in range(90):
    i=int(x2); cf2.append(i)
    if x2==i: break
    x2=1/(x2-i)
P0,Q0,P1,Q1=1,0,cf2[0],1
pairs=[]
convs=[(cf2[0],1)]
for t2 in cf2[1:]:
    P0,Q0,P1,Q1=P1,Q1,t2*P1+P0,t2*Q1+Q0
    convs.append((P1,Q1))
cand=set()
for i in range(len(convs)-1):
    pA,qA=convs[i]; pB,qB=convs[i+1]
    for c in range(0,(cf2[i+1] if i+1<len(cf2) else 1)+1):
        a=qA+c*qB
        if a<10**3 or a>10**24: continue
        cand.add(a)
        a2=qB+c*qA
        if 10**3<=a2<=10**24: cand.add(a2)
danger=[]; nearest=None
for a in sorted(cand):
    j=int(F(a)*ALF)+1
    if j>10**24: continue
    E=F(j)*lo2-F(a)*hi3
    if E<=0: continue
    t=j-2
    cap=F(7,3)+F(".94")*((t-8)//8) if t>=8 else F(7,3)
    ratio=cap/(E*floor_lo)
    if nearest is None or ratio>nearest[0]: nearest=(ratio,a,j,float(E))
    if ratio>=1: danger.append((a,j,float(E),float(ratio)))
print("  upper approximants examined with j <= 1e24: %d"%len(cand))
print("  dangerous pairs found: %d"%len(danger))
Jmin=min(d[1] for d in danger) if danger else 10**24
print("  LEAST dangerous j:  J* = %d  (= 10^%.3f)"%(Jmin,math.log10(Jmin)))
d0=min(danger,key=lambda z:z[1])
print("     at a = %d, E = %.6e, cap/(E*floor) = %.3f"%(d0[0],d0[2],d0[3]))
chk("J* exceeds every K in the certified family (checked below)", Jmin>10**22)
print("  Dangerous sites DO exist, but all have j >= J*.  So (C6) is the clean bound")
print("     K < J* = %d ,"%Jmin)
print("  established by ONE shell-independent enumeration of the upper approximants of")
print("  alpha.  It is conservative: it ignores byte-alignment and the five-one byte")
print("  condition, both of which only shrink the dangerous set further.")
print(); print("="*78)
print("2.  the parameterized transfer theorem")
print("="*78)
print("  THEOREM.  Let (K,D) satisfy")
print("    (C1) K = ceil(D alpha), K = 1 mod 8, 8+8*floor((K-17)/8) = K-9,")
print("         floor(8D/K) = 5, floor((K-9)D/K)+2 = D-4;")
print("    (C2) 0 < Lambda < 1/3000;")
print("    (C3) 20*eta < 1,  eta = 8D/K - 5;")
print("    (C4) m = ceil((sixes - 20)/(7*21^12)) >= 163;")
print("    (C5) (0.162 D/(3 Lambda) + (256/81) D/3)/2^104 < sigma_163 = 3611678508;")
print("    (C6) K < J*  [one shell-independent enumeration of upper approximants]")
print("  Then every integer n >= 2 realizing a complete 21/7 word on (K,D) descends at")
print("  its ACTUAL first coefficient crossing.")
print("  (C1)-(C5) are five cheap rational tests per shell; (C6) is a single global")
print("  bound; and the alphabet facts -- marker uniqueness, the per-byte memory costs,")
print("  and 'only 00111110 crosses early, at its second zero' -- are verified once.")

print(); print("="*78)
print("3.  certify every admissible shell")
print("="*78)
x,cf=AL,[]
for _ in range(80):
    i=int(x); cf.append(i)
    if x==i: break
    x=1/(x-i)
p0,q0,p1,q1=1,0,cf[0],1
conv=[(cf[0],1)]
for t in cf[1:]:
    p0,q0,p1,q1=p1,q1,t*p1+p0,t*q1+q0
    conv.append((p1,q1))
def test(K,D):
    if K%8!=1: return None
    nb=(K-17)//8
    if 8+8*nb!=K-9: return None
    if 8*D//K+2!=7 or (K-9)*D//K+2!=D-4: return None
    if not (5<F(8*D,K)<6): return None
    lam_lo,lam_hi=K*lo2-D*hi3, K*hi2-D*lo3
    if not (0<lam_lo and lam_hi<F(1,3000)): return None
    eta=F(8*D,K)-5
    if not 20*eta<1: return None
    sixes=D-11-5*nb
    m=(sixes-20+labels-1)//labels
    Q=F(D,3)/lam_lo
    need=(F("0.162")*Q+F(256,81)*F(D,3))/2**L
    fl=(2**L*SIGMA-F(256,81)*F(D,3))/F(".162")
    return dict(K=K,D=D,m=m,need=need,floor=fl,lam=lam_lo,
                c4=m>=163, c5=need<SIGMA, c6=K<Jmin,
                endpoint=fl*lam_lo>F(D,3))
seen=set(); shells=[]
for i in range(20,len(conv)-1):
    p,q=conv[i]; pp,qq=conv[i+1]
    for a in range(1,70):
        for b in range(1,70):
            Dx=a*q+b*qq; Kx=a*p+b*pp
            if Dx<10**20 or Dx>10**24: continue
            if Dx in seen: continue
            r=test(Kx,Dx)
            if r is None: continue
            seen.add(Dx); shells.append(r)
full=[r for r in shells if r['c4'] and r['c5'] and r['c6'] and r['endpoint']]
print("  shells meeting (C1)-(C3): %d"%len(shells))
print("  shells meeting ALL of (C1)-(C6) and the endpoint comparison: %d"%len(full))
chk("every shell passing (C1)-(C5) also satisfies the endpoint comparison",
    all(r['endpoint'] for r in shells if r['c4'] and r['c5']))
chk("both certified shells are in the list",
    any(r['K']==1412737287532845639425 for r in full)
    and any(r['K']==3179552387924051009817 for r in full))
print("  %10s %26s %10s %14s %s"%("log10 K","D","m","span room","c6"))
for r in sorted(full,key=lambda z:z['K'])[:10]:
    print("  %10.3f %26d %10d %13.1fx %s"
          %(math.log10(r['K']),r['D'],r['m'],SIGMA/float(r['need']),"ok" if r['c6'] else "FAIL"))
mn=min(SIGMA/float(r['need']) for r in full); mx=max(SIGMA/float(r['need']) for r in full)
print("  span room across the family: %.2fx to %.1fx"%(mn,mx))
print("  multiplicity range: %d to %d"%(min(r['m'] for r in full),max(r['m'] for r in full)))
print("  K range: 10^%.3f to 10^%.3f"%(math.log10(min(r['K'] for r in full)),
                                        math.log10(max(r['K'] for r in full))))
chk("the theorem covers %d distinct shells simultaneously"%len(full), len(full)>400,
    "%d shells"%len(full))
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
