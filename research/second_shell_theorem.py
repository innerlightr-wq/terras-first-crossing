"""A complete first-crossing certificate on a SECOND shell."""
from fractions import Fraction as F
from decimal import Decimal, getcontext
from itertools import product
import math
getcontext().prec=260
def lg(base,terms=240):
    x=F(base-1,base+1); lo=sum(2*x**i/i for i in range(1,terms,2))
    return lo, lo+2*x**(terms+1)/((terms+1)*(1-x*x))
lo2,hi2=lg(2); lo3,hi3=lg(3)
SIGMA=3611678508; L=104; PAST=7; FW=12; labels=7*21**FW
HEAD,TAIL="11111110","111100000"; MARKER="0111100000"
A5=["".join(w)+"0" for w in product("01",repeat=7) if w.count("1")==5]
A6=["".join(w)+"0" for w in product("01",repeat=7) if w.count("1")==6]
ALL=A5+A6; BAD="00111110"
def mcost(b): return sum(F(2**i,3**(b[:i].count("1")+1)) for i,bit in enumerate(b) if bit=="1")
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)

SHELLS=[("21/7 (regression)",1412737287532845639425,891337988684307508298),
        ("second shell",     3179552387924051009817,2006074204580460545110)]
print("="*78); print("alphabet-only facts (shell independent)"); print("="*78)
chk("marker unique over all %d three-byte products"%len(ALL)**3,
    all(MARKER not in a+b+c for a,b,c in product(ALL,repeat=3)))
chk("marker absent at head junctions and unique at tail junctions",
    all(MARKER not in HEAD+a+b for a,b in product(ALL,repeat=2))
    and all((a+b+TAIL).find(MARKER)==len(a+b+TAIL)-10 for a,b in product(ALL,repeat=2)))
chk("only 00111110 crosses early, at its second zero",
    all(F(3**w[:j].count("1"),2**j)*F("2.999")>1 for w in ALL if w!=BAD for j in range(1,9))
    and F(3**BAD[:2].count("1"),2**2)==F(1,4))
chk("per-byte memory: earlier bad byte at rho>=4 and all others at rho>2.999 cost < 0.94",
    mcost(BAD)/4<F(".94") and max(mcost(w) for w in ALL if w!=BAD)/F("2.999")<F(".94")
    and mcost(HEAD)<F(7,3))

for name,K,D in SHELLS:
    print(); print("="*78); print("%s:  K = %d"%(name,K)); print("="*78)
    lam_lo,lam_hi=K*lo2-D*hi3, K*hi2-D*lo3
    eta=F(8*D,K)-5; nb=(K-17)//8; end=K-9; sixes=D-11-5*nb
    def a_of(t): return t*D//K+2
    chk("K = 1 mod 8, 8+8nb = K-9, a(8)=7, a(K-9)=D-4, 5<8D/K<6",
        K%8==1 and 8+8*nb==end and a_of(8)==7 and a_of(end)==D-4 and 5<F(8*D,K)<6)
    chk("0 < Lambda < 1/3000", 0<lam_lo and lam_hi<F(1,3000),
        "Lambda in (%.6e, %.6e)"%(float(lam_lo),float(lam_hi)))
    chk("20*eta < 1, so 7 behind and 12 ahead are all five-one", 20*eta<1,
        "20*eta = %.6f"%float(20*eta))
    chk("boundary rho > 2.999", 3*(1-lam_hi)>F("2.999"))
    chk("rho < 3^(1+eta) < 256/81 before a six-one byte (equivalent to Lambda>0)",
        (8*D-4*K)*hi3 < 8*K*lo2 and F(256,81)-F("2.999")<F(".162"))
    chk("all 8 head prefixes above the line",
        all(F(HEAD[:j].count("1"))-F(j)*(hi2/lo3) > 0 for j in range(1,9)))
    chk("all 9 tail starts above the line",
        all(F(len(TAIL)-j)*lo2-F(TAIL[j:].count("1"))*hi3 > lam_hi for j in range(len(TAIL))))
    m=(sixes-PAST-FW-1+labels-1)//labels
    Q=F(D,3)/lam_lo
    need=(F("0.162")*Q+F(256,81)*F(D,3))/2**L
    floor=(2**L*SIGMA-F(256,81)*F(D,3))/F(".162")
    chk("m = %d >= 163, so sigma_163 is reusable"%m, m>=163)
    chk("V/2^104 < sigma_163  (%.2fx room)"%(SIGMA/float(need)), need<SIGMA<3**(5*PAST))
    chk("ENDPOINT: floor*Lambda > D/3", floor*lam_lo>F(D,3),
        "%.4e > %.4e"%(float(floor*lam_lo),float(F(D,3))))
    inv=pow(D,-1,K); hb=int(K*lam_hi/lo3)+1
    sites=[];coarse=0
    for h in range(hb+1):
        j=(-h*inv)%K; t=j-2
        if not (8<=t<end and t%8==0): continue
        a=a_of(t)
        if K*a-D*j!=h or a_of(t+8)!=a+5: continue
        elo,ehi=j*lo2-a*hi3, j*hi2-a*lo3
        if ehi<0 or elo<=0: continue
        if floor*elo>F(a,3): coarse+=1; continue
        cap=F(7,3)+F(".94")*((t-8)//8)
        sites.append((h,j,t,a,cap/elo))
    print("    h-candidates %d -> byte-aligned sites %d (coarse test removed %d)"
          %(hb+1,len(sites),coarse))
    chk("EVERY early-crossing site excluded: cap/E < floor",
        all(s[4]<floor for s in sites),
        "worst ratio %.6f"%(max(float(s[4]/floor) for s in sites) if sites else 0.0))
    if sites:
        w=max(sites,key=lambda s:s[4])
        print("    tightest: h = %d, j = %d, cap/E = %.6e, floor = %.6e"
              %(w[0],w[1],float(w[4]),float(floor)))
print()
print("="*78); print("why the sites transfer"); print("="*78)
print("  A site is a pair (j, a) with t = j-2 byte-aligned and a = floor(tD/K)+2.")
print("  Both shells have D/K within 1e-39 of alpha, so floor(tD/K) agrees for every")
print("  t ~ 1e21: the SAME (j,a) pairs are sites on both.  Its ceiling cap/E depends")
print("  only on (j,a), and the floor is dominated by 2^104 sigma, which is")
print("  shell-independent.  Hence the comparison cap/E < floor transfers verbatim;")
print("  only the ADMISSIBLE RANGE of h, namely K*Lambda/ln3, is shell-dependent.")
print("  That is why a shell with more h-candidates is not harder to finish.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
