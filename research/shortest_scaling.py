"""Part 4: how far does the interior comparison scale?  And a correction to the
previous round's report."""
from fractions import Fraction as F
from decimal import Decimal, getcontext
import math
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
def affine(bits):
    C=d=0
    for j,b in enumerate(bits):
        if b=='1': C=3*C+(1<<j); d+=1
    return C,d
def first_cross_at_end(bits):
    oc=0
    for t,b in enumerate(bits,1):
        oc+=int(b)
        if t<len(bits) and 3**oc<(1<<t): return False,t
    C,d=affine(bits)
    return (3**d<(1<<len(bits))), None

print("="*78)
print("4a. CORRECTION to research/JUNE_CORRIDOR_AUDIT.md")
print("="*78)
d,k=306,485
pos={(3**i).bit_length()-1 for i in range(d)}
X8=''.join('1' if j in pos else '0' for j in range(k))
for b in (1,2,3):
    ok,badt=first_cross_at_end(X8*b)
    print("  X8^%d : length %4d, first crossing at the end? %-5s %s"
          %(b,485*b,ok,"" if ok else "-- it already crossed at t = %d"%badt))
chk("X8 alone is a first-crossing word", first_cross_at_end(X8)[0])
chk("X8^b for b >= 2 is NOT a first-crossing word (it crosses at t = 485)",
    not first_cross_at_end(X8*2)[0] and first_cross_at_end(X8*2)[1]==485)
print("  My previous report treated X8^b as a WORD family when computing the vacuity")
print("  threshold.  It is only a SHELL family: 2^485 > 3^306, so the coefficient has")
print("  already fallen below one at t = 485 and X8^b never reaches t = 485b without")
print("  crossing.  The conclusion of that section (a fixed factor cannot follow an")
print("  unbounded family, because the admissible window grows while F stays fixed)")
print("  is unaffected -- the window width is M_w, which grows with d on any")
print("  line-hugging family -- but the vehicle was wrong and is withdrawn.")
print("  Valid unbounded shell families with 929 | 2^k - 3^d do exist: the pairs (k,d)")
print("  with 2^k = 3^d mod 929 form a lattice coset, and infinitely many of them are")
print("  critical (k = ceil(d alpha)).  Exhibit the first few:")
cnt=0
getcontext().prec=120
AL=Decimal(3).ln()/Decimal(2).ln()
found=[]
for dd in range(306, 400000):
    kk=int(AL*dd)+1
    if (pow(2,kk,929)-pow(3,dd,929))%929==0:
        found.append((kk,dd)); cnt+=1
        if cnt>=6: break
print("    critical shells (k,d) with 929 | 2^k - 3^d:", found)
chk("infinitely many critical shells carry the factor 929 (6 exhibited)", len(found)>=6)

print()
print("="*78); print("4b. the gain factor of the interior comparison"); print("="*78)
print("  U = (delta_j rho_k M_w + eps rho_j M_j)/(delta_j + eps), so")
print("     U ~ M_j + (M_w - M_j) * delta_j/(delta_j + eps).")
print("  The gain is real only when delta_j << eps, i.e. the INTERIOR drift is much")
print("  smaller than the FINAL deficit.  Both are values of the linear form:")
print("     delta_j ~ (alpha d_j - j) ln2 ,   eps ~ (k - alpha d) ln2 = Lambda .")
B=28*10**18
print("  But a failure forces Lambda < M_w/n < M_w/B.  With M_w ~ 74 that is")
print("     eps < %.4e ." % (74.0/B))
print("  And alpha d_j - j > 0 with d_j <= d, so by best approximation delta_j >= the")
print("  least positive value of alpha m - integer over m <= d_j.  Find the smallest")
print("  denominator that can reach eps:")
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
target=F(74,B)
first=None
for (pp,qq) in conv:
    val=F(pp)-F(qq)*F(int(AL*10**60),10**60)   # p - q*alpha ; lower convergents give <0
    g=abs(val)
    if g < target and first is None and qq>1:
        first=(qq,float(g))
print("    smallest convergent denominator q with |q*alpha - p| < %.3e :  q = %d (gap %.3e)"
      %(float(target),first[0],first[1]))
chk("reaching delta_j <= eps needs d_j of order 1e17 or more", first[0] > 10**16,
    "q = %d"%first[0])
print("  so BOTH d and d_j must be high-quality denominators, with d_j <= d yet a")
print("  STRICTLY BETTER approximation.  That is governed by the continued fraction of")
print("  alpha: if d_j = q_m and d is a neighbouring intermediate denominator, then")
print("     eps/delta_j  ~  a_{m+1} ,   the next partial quotient.")
print("  partial quotients of log2 3 up to the relevant range:")
qs=[q for (_,q) in conv]
rel=[(i,cf[i],qs[i-1] if i>0 else 1) for i in range(1,len(cf)) if qs[i-1]<10**20]
print("    max partial quotient with denominator below 1e20: %d"
      %max(a for _,a,_ in rel))
print("    the ten largest:", sorted((a for _,a,_ in rel), reverse=True)[:10])
chk("the partial quotients of log2 3 are small, bounding the gain factor",
    max(a for _,a,_ in rel) < 1000, "max = %d"%max(a for _,a,_ in rel))

print()
print("="*78); print("4c. what the gain can and cannot do"); print("="*78)
print("  To make the factor sieve bite we need U < F for a KNOWN factor F of 2^k - 3^d,")
print("  so F is small, hence we need U = O(1).  Two regimes:")
print()
print("  (i) M_w already small.  Then U <= M_w is already O(1) and the interior")
print("      comparison adds nothing -- and result 2 of this round shows M_w > 4.2472")
print("      is forced for d < 4.64e18, so the window is at least {0,...,4}.")
print("  (ii) M_w large (a line-hugging word, M_w ~ d/6 by result 1).  Then")
print("      U ~ M_j + M_w/a with a a partial quotient of log2 3, which is at most")
print("      %d in range.  With d >= 6.59e9 that leaves U >= %.3e -- still vastly"
      %(max(a for _,a,_ in rel), 6.59e9/6/max(a for _,a,_ in rel)))
print("      larger than any factor of 2^k - 3^d one could find.")
print("  So the interior comparison cannot convert a large-memory word into a")
print("  small-window one.  It is a per-word settler, not a family tool.  [PROVED]")
print()
print("  EMPIRICAL best case: minimise U over all interior j for real first-crossing")
print("  words, and report U/M_w.")
def best_U(bits):
    C,dd=affine(bits); PP,QQ=3**dd,1<<len(bits); AA=QQ-PP
    best=None
    for j in range(1,len(bits)):
        Cj,dj=affine(bits[:j]); pp,qq=3**dj,1<<j
        den=pp*QQ-qq*PP
        if den<=0: continue
        U=F((pp-qq)*C+AA*Cj,den)
        if best is None or U<best: best=U
    return best, F(C,PP)
rows=[]
for nm,w in (("X8 (critical line)",X8),
             ("X8 with 2 swaps", X8[:101]+'1101'+X8[105:410]+'1101'+X8[414:])):
    ok,_=first_cross_at_end(w)
    if not ok: continue
    U,Mw=best_U(w); rows.append((nm,float(U),float(Mw)))
    print("    %-22s min_j U = %10.4f   M_w = %10.4f   ratio %.4f"
          %(nm,float(U),float(Mw),float(U/Mw)))
chk("the best interior comparison still leaves U of the same order as M_w",
    all(r[1]/r[2] > 0.2 for r in rows), "ratios %s"%[round(r[1]/r[2],3) for r in rows])
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
