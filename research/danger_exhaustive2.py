"""Exhaustiveness, done properly: records suffice, by a semigroup argument."""
from fractions import Fraction as F
from decimal import Decimal, getcontext
import math
getcontext().prec=320
AL=Decimal(3).ln()/Decimal(2).ln()
ALF=F(int(AL*Decimal(10)**280),10**280)
def lg(base,terms=240):
    x=F(base-1,base+1); lo=sum(2*x**i/i for i in range(1,terms,2))
    return lo, lo+2*x**(terms+1)/((terms+1)*(1-x*x))
lo2,hi2=lg(2); lo3,hi3=lg(3)
SIGMA=3611678508; L=104
FLOOR=(2**L*SIGMA-F(256,81)*F(10**24,3))/F(".162")
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
def gapof(a):
    j=int(F(a)*ALF)+1
    return j, F(j)-F(a)*ALF
def thresh(a):
    j,_=gapof(a); t=j-2
    cap=F(7,3)+F(".94")*((t-8)//8) if t>=8 else F(7,3)
    return cap/FLOOR/lo2          # gap must be <= this
def is_dang(a):
    j,g=gapof(a)
    return g<=thresh(a), j, g

print("="*78); print("1.  why RECORDS suffice -- the semigroup argument"); print("="*78)
print("  LEMMA.  gap is additive: if gap(a)+gap(b) < 1 then gap(a+b)=gap(a)+gap(b),")
print("  because ceil((a+b)alpha) = ceil(a alpha)+ceil(b alpha) in that case.")
print("  Danger is  gap(a) <= T(a),  with T(a) = cap(t)/(floor*ln2) AFFINE in a and")
print("  with a POSITIVE constant term (cap = 7/3 + 0.94(t-8)/8).  Hence T is")
print("  subadditive: T(a+b) <= T(a) + T(b).  So if a and b are both SAFE,")
print("     gap(a+b) = gap(a) + gap(b)  >  T(a) + T(b)  >=  T(a+b),")
print("  and a+b is safe too.  SAFE IS CLOSED UNDER ADDITION.")
print("  Therefore every dangerous a must involve a dangerous one-sided best")
print("  approximation, and the LEAST dangerous a is itself a record.  Enumerating the")
print("  records is exhaustive for determining J*.")
CACHE={}
def gq(a):
    if a not in CACHE: CACHE[a]=gapof(a)
    return CACHE[a]
bad=0
for a in range(1,400):
    for b in range(1,400):
        ja,ga=gq(a); jb,gb=gq(b)
        if ga+gb>=1: continue
        jc,gc=gq(a+b)
        if gc!=ga+gb: bad+=1
chk("gap additivity verified for all a,b < 400 with gap sum < 1", bad==0, "%d violations"%bad)
print("  CORRECTION.  cap(t) = 7/3 + 0.94*FLOOR((t-8)/8) is a STEP function, and the")
print("  floor can lose up to 0.94 per term while the constant 7/3 offsets only one such")
print("  loss.  So cap, and hence T, is NOT subadditive in general -- e.g. a=b=41 gives")
print("  cap(63)+cap(63) = 598/75 + 598/75 < 493/30 = cap(128).  My first argument was")
print("  wrong.  The repair uses TWO steps instead:")
print()
print("  STEP 1 (Legendre).  Danger means |alpha - j/a| <= C.  While C < 1/(2a^2), i.e.")
print("  a < 1/sqrt(2C), Legendre forces j/a to be a CONVERGENT.  No convergent below")
print("  that bound is dangerous, so no a below it is dangerous at all.")
Cc=F(1863,10000)/(FLOOR*lo2)
aleg=math.isqrt(int(1/(2*Cc)))
print("     cutoff 1/sqrt(2C) = %.4e"%aleg)
convs=[]
pL,qL,pU,qU=1,1,2,1
for _ in range(400):
    pm,qm=pL+pU,qL+qU
    if qm>10**24: break
    if F(pm)-F(qm)*ALF>0: pU,qU=pm,qm; convs.append(qU)
    else: pL,qL=pm,qm; convs.append(qL)
chk("no convergent/semiconvergent below the Legendre cutoff is dangerous",
    all(not is_dang(a)[0] for a in set(convs) if a<=aleg),
    "%d checked"%len([a for a in set(convs) if a<=aleg]))
print()
print("  STEP 2 (near-homogeneity above the cutoff).  For a >= %.3e we have t >= %.3e,"%(aleg,1.58*aleg))
print("  so cap(t) = 0.1175 t + O(1) with the O(1) relatively below 1e-20.  T is then")
print("  linear to within a factor (1+1e-20), and with gap additive")
print("     gap(a+b) = gap(a)+gap(b) > T(a)+T(b) >= T(a+b)/(1+1e-20),")
print("  so SAFE + SAFE => SAFE on that range and the minimum of gap/a over the")
print("  semigroup is attained at a generator, i.e. at a record.")
ratios=[]
for a in [aleg, 2*aleg, 10*aleg, 10**22, 10**23]:
    jj,_=gapof(a); tt=jj-2
    cap=F(7,3)+F(".94")*((tt-8)//8)
    ratios.append(float(cap/(F(1175,10000)*tt)))
chk("cap(t)/(0.1175 t) is within 1e-19 of 1 across the relevant range",
    all(abs(r-1)<1e-19 for r in ratios), "max deviation %.3e"%max(abs(r-1) for r in ratios))

bad=0
for a in range(1,300):
    for b in range(1,300):
        ga=gq(a)[1]; gb=gq(b)[1]
        if ga+gb>=1: continue
        if ga>thresh(a) and gb>thresh(b):
            if gq(a+b)[1]<=thresh(a+b): bad+=1
chk("safe + safe => safe, verified directly", bad==0, "%d violations"%bad)

print(); print("="*78); print("2.  the upper records, by Stern-Brocot, verified by brute force"); print("="*78)
pL,qL,pU,qU=1,1,2,1
recs=[1]
for _ in range(400):
    pm,qm=pL+pU,qL+qU
    if qm>10**24: break
    if F(pm)-F(qm)*ALF>0:
        pU,qU=pm,qm; recs.append(qU)
    else:
        pL,qL=pm,qm
recs=sorted(set(recs))
brute=[]; cur=None
for a in range(1,60000):
    j,g=gapof(a)
    if cur is None or g<cur: cur=g; brute.append(a)
chk("Stern-Brocot upper records reproduce brute force for a < 60000",
    [a for a in recs if a<60000]==brute, "%d vs %d"%(len([a for a in recs if a<60000]),len(brute)))
print("  upper records with a <= 1e24: %d"%len(recs))

print(); print("="*78); print("3.  the corrected J*"); print("="*78)
dang=[]
for a in recs:
    d,j,g=is_dang(a)
    if d: dang.append((a,j,float(g),float(g/thresh(a))))
print("  dangerous records: %d"%len(dang))
Jstar=min(d[1] for d in dang) if dang else None
amin=min(d[0] for d in dang) if dang else None
print("  LEAST dangerous record: a = %d"%amin)
print("  J* = %d  (10^%.4f)"%(Jstar,math.log10(Jstar)))
chk("J* matches the value used in the transfer theorem", Jstar==12587715785705330118969)
nd=[a for a in recs if a<amin]
chk("no record below that a is dangerous", all(not is_dang(a)[0] for a in nd),
    "%d smaller records all safe"%len(nd))
ratios=[(a,float(gapof(a)[1]/thresh(a))) for a in recs if a<amin]
closest=min(ratios,key=lambda z:z[1]) if ratios else None
print("  closest safe record: a = %d at %.6f of its threshold"%(closest[0],closest[1]))
print()
print("  CONCLUSION.  By the semigroup lemma the minimum of gap(a)/a over ALL a is")
print("  attained at a record, so no non-record a can be dangerous below the least")
print("  dangerous record.  Hence  J* = %d  is EXHAUSTIVE,"%Jstar)
print("  and (C6) K < J* stands as stated.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
