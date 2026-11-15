"""Phase 2 + 4: the reusable criterion, and a semiconvergent search for a second shell."""
from fractions import Fraction as F
from decimal import Decimal, getcontext
import math
getcontext().prec=260
AL=Decimal(3).ln()/Decimal(2).ln()
ALF=F(int(AL*Decimal(10)**220),10**220)
LN2=F(int(Decimal(2).ln()*Decimal(10)**220),10**220)
LN3=F(int(Decimal(3).ln()*Decimal(10)**220),10**220)
K0,D0=1412737287532845639425,891337988684307508298
SIGMA=3611678508; L=104; PAST=7; FW=12; S=21
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)

print("="*78); print("PHASE 4.  the reusable criterion"); print("="*78)
print("  sigma_m(R, 3^t) depends ONLY on the block alphabet, the incoming depth and m.")
print("  It is INDEPENDENT of the shell.  And m' >= m implies sigma_{m'} >= sigma_m.")
print("  So any shell forcing m >= 163 may reuse sigma = 3611678508 with NO new")
print("  enumeration.  The criterion becomes, with Q = D/(3 Lambda) and M = D/3:")
print("     (i)   5 < 8D/K < 6                                  [legal schedule]")
print("     (ii)  past + F <= floor(1/eta) - 2 = 19              [six-one separation]")
print("     (iii) m = ceil((sixes - past - F - 1)/(7 s^F)) >= 163")
print("     (iv)  (0.162 Q + (256/81) M) / 2^L  <  sigma_163")
print("     (v)   the terminal marker is unique                  [alphabet only]")
print("     (vi)  K*Lambda/ln3 candidate sites can be enumerated [early crossings]")
def shell_data(Dx):
    Kx=int(F(Dx)*ALF)+1
    lam=F(Kx)*LN2-F(Dx)*LN3
    if lam<=0: return None
    nb=(Kx-17)//8
    sixes=Dx-11-5*nb
    if sixes<=0: return None
    legal = 5 < F(8*Dx,Kx) < 6
    labels=7*S**FW
    m=(sixes-PAST-FW-1+labels-1)//labels
    Q=F(Dx,3)/lam; M=F(Dx,3)
    V=F("0.162")*Q+F(256,81)*M
    need=V/2**L
    return dict(K=Kx,D=Dx,lam=lam,sixes=sixes,legal=legal,m=m,Q=Q,V=V,need=need,
                hcount=float(Kx*lam/LN3))
d0=shell_data(D0)
chk("the 21/7 shell satisfies (i),(iii),(iv) with the reused sigma",
    d0['legal'] and d0['m']>=163 and d0['need']<SIGMA,
    "m = %d, need = %.4e < sigma = %.4e (%.1fx room)"
    %(d0['m'],float(d0['need']),SIGMA,SIGMA/float(d0['need'])))
chk("and (vi): only %d candidate early-crossing sites"%int(d0['hcount']), d0['hcount']<200)
print("  REGRESSION: the original result is recovered by this criterion.")

print(); print("="*78); print("PHASE 2.  semiconvergent search"); print("="*78)
print("  The best convergents have Lambda FAR too small (D/Lambda explodes, so (iv)")
print("  fails).  The 21/7 D is NOT a convergent denominator -- it is an intermediate")
print("  fraction chosen so Lambda lands in the window where (iv) and (vi) both hold.")
x,cf=AL,[]
for _ in range(70):
    i=int(x); cf.append(i)
    if x==i: break
    x=1/(x-i)
p0,q0,p1,q1=1,0,cf[0],1
conv=[(cf[0],1)]
for t in cf[1:]:
    p0,q0,p1,q1=p1,q1,t*p1+p0,t*q1+q0
    conv.append((p1,q1))
chk("D0 is not a convergent denominator", all(q!=D0 for _,q in conv))
# semiconvergents q = q_{m-1} + c*q_m
cands=[]
for i in range(2,len(conv)-1):
    pm1,qm1=conv[i-1]; pm,qm=conv[i]
    amax=cf[i+1] if i+1<len(cf) else 1
    for c in range(1,amax+1):
        q=qm1+c*qm
        if q<10**18 or q>10**26: continue
        d=shell_data(q)
        if d is None or not d['legal']: continue
        cands.append(d)
chk("D0 appears among the semiconvergents", any(d['D']==D0 for d in cands),
    "%d semiconvergent shells in 1e18..1e26"%len(cands))
ok=[d for d in cands if d['m']>=163 and d['need']<SIGMA and d['hcount']<10**6]
print("  %24s %12s %10s %14s %12s %s"%("D","log10 K","m","need/sigma","K*Lam/ln3","usable"))
for d in sorted(cands,key=lambda z:z['D'])[:40]:
    u = d['m']>=163 and d['need']<SIGMA and d['hcount']<10**6
    print("  %24d %12.2f %10d %14.4e %12.4e %s"
          %(d['D'],math.log10(d['K']),min(d['m'],10**9),float(d['need'])/SIGMA,
            d['hcount'],"YES" if u else "no"))
chk("at least one shell OTHER than the 21/7 one passes every cheap test",
    any(d['D']!=D0 for d in ok), "%d usable, %d of them new"%(len(ok),len([d for d in ok if d['D']!=D0])))
if ok:
    newo=[d for d in ok if d['D']!=D0]
    if newo:
        b=max(newo,key=lambda z:SIGMA/float(z['need']))
        print()
        print("  STRONGEST NEW CANDIDATE (largest anticipated margin):")
        print("    D = %d"%b['D']); print("    K = %d"%b['K'])
        print("    forced multiplicity m = %d  (>= 163 required)"%b['m'])
        print("    required span V/2^104 = %.6e"%float(b['need']))
        print("    reusable sigma_163     = %d"%SIGMA)
        print("    anticipated margin     = %.2fx"%(SIGMA/float(b['need'])))
        print("    early-crossing sites   = %.4e"%b['hcount'])
        print("    six-one bytes          = %.4e"%float(b['sixes']))
        print("  Its certificate still needs: marker uniqueness (alphabet only, already")
        print("  verified), the exact early-crossing filter at %d sites, and the per-byte"%int(b['hcount']))
        print("  memory bound -- but NO new residue enumeration.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
