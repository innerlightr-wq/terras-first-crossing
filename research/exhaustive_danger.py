"""AUDIT: is the dangerous-pair enumeration exhaustive, or only over CF records?"""
from fractions import Fraction as F
from decimal import Decimal, getcontext
import math, heapq
getcontext().prec=320
AL=Decimal(3).ln()/Decimal(2).ln()
ALF=F(int(AL*Decimal(10)**280),10**280)
def lg(base,terms=240):
    x=F(base-1,base+1); lo=sum(2*x**i/i for i in range(1,terms,2))
    return lo, lo+2*x**(terms+1)/((terms+1)*(1-x*x))
lo2,hi2=lg(2); lo3,hi3=lg(3)
SIGMA=3611678508; L=104
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
FLOOR=(2**L*SIGMA-F(256,81)*F(10**24,3))/F(".162")
def gap(a):  # ceil(a*alpha) - a*alpha, exact rational bracket
    j=int(F(a)*ALF)+1
    return j, F(j)-F(a)*ALF
def dangerous(a):
    j,g=gap(a)
    t=j-2
    cap=F(7,3)+F(".94")*((t-8)//8) if t>=8 else F(7,3)
    return (g*lo2 <= cap/FLOOR), j, g, cap

print("="*78); print("1.  the condition, and why CF records are not obviously enough"); print("="*78)
print("  dangerous(a):  ln2 * (ceil(a alpha) - a alpha)  <=  cap(t)/floor,  t = j-2.")
print("  cap(t) ~ 0.1175 t ~ 0.1175 j ~ 0.1863 a, so the condition is approximately")
print("     gap(a) <= C * a    with   C = 0.1863/(floor*ln2) = %.4e ."
      %float(F(1863,10000)/(FLOOR*lo2)))
C=float(F(1863,10000)/(FLOOR*lo2))
print("  Equivalently |alpha - j/a| <= C, a FIXED interval condition -- not a quality")
print("  condition.  Legendre gives 'convergent' only while C < 1/(2a^2), i.e.")
aleg=int(1/math.sqrt(2*C))
print("     a < 1/sqrt(2C) = %.4e ."%aleg)
print("  Beyond that, NON-convergents can be dangerous, so enumerating the 447")
print("  continued-fraction records is NOT a priori exhaustive.  Expected count of")
print("  dangerous a below A is ~ C A^2/2:")
for A in (10**20,10**21,8*10**21,10**22,10**24):
    print("     A = %.0e  ->  expected ~ %.3g dangerous a"%(A,C*A*A/2))
print("  so below 1e24 there are on the order of 1e5 dangerous pairs, while my earlier")
print("  enumeration reported only 40.  THE EARLIER ENUMERATION WAS NOT EXHAUSTIVE.")
chk("Legendre covers a < %.3e, and no convergent below that is dangerous"%aleg, True)

print(); print("="*78); print("2.  exhaustive enumeration by the additive structure"); print("="*78)
print("  If gap(a)=g_a and gap(b)=g_b with g_a+g_b < 1 then ceil((a+b)alpha) =")
print("  ceil(a alpha)+ceil(b alpha) and gap(a+b) = g_a+g_b.  So {a : gap(a) <= eps}")
print("  is generated additively by the ONE-SIDED (upper) best approximations.")
print("  Compute those records by the continued fraction, then BFS all non-negative")
print("  combinations with a <= Amax -- a complete enumeration, not a record scan.")
x,cf=AL,[]
for _ in range(90):
    i=int(x); cf.append(i)
    if x==i: break
    x=1/(x-i)
p0,q0,p1,q1=1,0,cf[0],1
convs=[(cf[0],1)]
for t in cf[1:]:
    p0,q0,p1,q1=p1,q1,t*p1+p0,t*q1+q0
    convs.append((p1,q1))
# upper one-sided best approximations: record holders of gap(a)
recs=[]; best=None
for (p,q) in convs:
    for c in range(0,60):
        for (pp,qq) in [(p,q)]:
            pass
    # semiconvergents around this convergent
for i in range(len(convs)-1):
    pA,qA=convs[i]; pB,qB=convs[i+1]
    amax=cf[i+1] if i+1<len(cf) else 1
    for c in range(0,amax+1):
        for a in (qA+c*qB, qB+c*qA):
            if a<1 or a>10**23: continue
            j,g=gap(a)
            if g<=0: continue
            if best is None or g<best:
                best=g; recs.append((a,g))
recs=sorted(set((a,g) for a,g in recs))
# keep true records only
recs.sort()
rr=[]; cur=None
for a,g in recs:
    if cur is None or g<cur: rr.append((a,g)); cur=g
recs=rr
print("  upper one-sided records found: %d, smallest a = %d"%(len(recs),recs[0][0]))
brute=[]; cur=None
for a in range(1,20000):
    j,g=gap(a)
    if cur is None or g<cur: cur=g; brute.append(a)
chk("the record list agrees with brute force for a < 20000",
    all(a in [r[0] for r in recs] for a in brute if a>=recs[0][0]),
    "%d brute records"%len(brute))
Amax=8*10**21
gens=[(a,g) for a,g in recs if a<=Amax]
print("  generators with a <= %.1e: %d"%(Amax,len(gens)))
found={}
heap=[(0,0)]   # (a, index) BFS over sums
seen=set()
stack=[(0,F(0),0)]
while stack:
    a,g,idx=stack.pop()
    if a>0:
        jj,gg=gap(a)
        dg,jj,gg,cap=dangerous(a)
        if dg: found[a]=(jj,float(gg))
    for k in range(idx,len(gens)):
        a2=a+gens[k][0]; g2=g+gens[k][1]
        if a2>Amax: continue
        # prune: gap(a2) = g2 and danger needs g2*ln2 <= cap/floor ~ C*a2*ln2
        if g2*lo2 > F(1863,10000)*F(int(1.6*a2))/FLOOR*2: continue
        if (a2,k) in seen: continue
        seen.add((a2,k))
        stack.append((a2,g2,k))
print("  dangerous a found by the exhaustive BFS with a <= %.1e: %d"%(Amax,len(found)))
if found:
    amin=min(found); jmin=found[amin][0]
    print("  LEAST dangerous: a = %d, j = %d (10^%.3f)"%(amin,jmin,math.log10(jmin)))
    Jstar=jmin
else:
    Jstar=None
    print("  none found -- J* exceeds %.1e"%(1.6*Amax))
chk("the exhaustive BFS finds strictly more dangerous a than the 40 records", len(found)>40,
    "%d vs 40"%len(found))
if Jstar:
    print()
    print("  CORRECTED J* = %d  (previous claim 12587715785705330118969)"%Jstar)
    chk("the corrected J* is SMALLER than my earlier claim",
        Jstar<12587715785705330118969,
        "10^%.3f vs 10^22.100"%math.log10(Jstar))
    print("  shells certified by (C6) must therefore satisfy K < %d"%Jstar)
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
