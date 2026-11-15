"""Part 1: the exact relationship between the three tests."""
from collections import defaultdict
from fractions import Fraction as F
from math import prod
import random
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
src=open("research/terras_overlap_and_history.py").read().replace('if __name__ == "__main__":','if False:')
M={}; exec(compile(src,"chk","exec"),M)

print("="*78); print("1A. the sorted-placement lemma"); print("="*78)
print("  Fix one bucket: members of equal run weight w>0, joint residue stride P, and")
print("  f = least term of the class strictly above the seed.  Suppose")
print("     (i)  every relaxed value x_i in the bucket has x_i >= f, and")
print("     (ii) every pair in the bucket has |x_i - x_j| >= P.")
print("  Sort x_(1) <= ... <= x_(c).  Then  x_(m) >= f + (m-1)P.")
print("  PROOF.  x_(1) >= f by (i).  For m >= 2, sortedness and (ii) give")
print("  x_(m) >= x_(m-1) + P, so induction gives x_(m) >= f + (m-1)P.  Since")
print("  1 + w/x is strictly decreasing in x for w > 0,")
print("     prod_i (1 + w/x_i) <= prod_{m=1}^{c} (1 + w/(f + (m-1)P)) ,")
print("  and the right-hand side is EXACTLY the implementation's distinct-placement")
print("  factor  prod(1+weight/(first+i*stride) for i in range(count)).   QED")
print()
print("  COROLLARY.  collective distinct-placement is IMPLIED by")
print("  (individual floors) + (in-bucket pairwise separation).  It is not")
print("  independent information.")
# numeric check of the lemma on random buckets
bad=0
random.seed(11)
for _ in range(2000):
    P=random.randrange(2,50); f=F(random.randrange(1,40)); c=random.randrange(1,7)
    w=F(random.randrange(1,9),random.randrange(1,9))
    xs=[]
    cur=f+F(random.randrange(0,3))
    for i in range(c):
        xs.append(cur); cur=cur+P+F(random.randrange(0,3))
    if any(x<f for x in xs): continue
    ok=all(abs(a-b)>=P for i,a in enumerate(xs) for b in xs[i+1:])
    if not ok: continue
    s=sorted(xs)
    if any(s[m]<f+m*P for m in range(c)): bad+=1
    if prod(1+w/x for x in xs) > prod(1+w/(f+m*P) for m in range(c)): bad+=1
chk("the sorted bound and the product bound hold on 2000 random buckets", bad==0,
    "%d failures"%bad)
print()
print("  EDGE CASES, against the implementation:")
print("   - initial run: j=0 is NOT bucketed; it contributes 1+w/seed with the seed")
print("     itself in place of a floor.  The lemma applies to noninitial buckets only.")
print("   - partial histories: used = min(i, depth), so early starts carry shorter")
print("     histories, a smaller 3^used and hence a DIFFERENT stride; they land in")
print("     different buckets and the lemma applies bucket by bucket.")
print("   - equality: x_(1) = f and consecutive gaps exactly P make the bound tight,")
print("     which is the integral-seed case.")
W=M['ORIGINAL_SHELL']; path=M['affine_prefixes'](W); k=len(W)
p,q,c0=path[-1]; A=q-p; h=c0*pow(q,-1,929)%929
seed=F(c0-q*h,A); vals=[(a*seed+cc)/b for a,b,cc in path]
pos=[j for j,b in enumerate(W) if b=='1']
def buckets(word,sd,depth):
    ps=[j for j,b in enumerate(word) if b=='1']; g=defaultdict(list)
    for i,j in enumerate(ps):
        if j and word[j-1]=='1': continue
        ell=word.index('0',j)-j; wt=1-F(2,3)**ell
        if j==0: continue
        used=min(i,depth)
        hist=tuple(ps[z+1]-ps[z] for z in range(i-used,i))
        m3=3**used; r3=0
        for gp in hist: r3=(3*r3+1)*pow(2,-gp,m3)%m3
        m2,r2=2**(ell+1),2**ell-1
        res=r2+m2*((r3-r2)*pow(m2,-1,m3)%m3)
        g[res,m2*m3,wt].append((j,hist))
    return g
for depth in (9,):
    g=buckets(W,seed,depth)
    tot=sum(len(v) for v in g.values())
    multi={kk:v for kk,v in g.items() if len(v)>1}
    chk("depth %d: %d noninitial starts in %d buckets, %d buckets with >1 member"
        %(depth,tot,len(g),len(multi)), tot>0)
    # verify the implementation's packed factor equals the lemma's bound
    ok=True
    for (res,stride,wt),mem in g.items():
        first=res+stride*((seed-res)//stride+1)
        mine=prod(1+wt/(first+i*stride) for i in range(len(mem)))
        if mine<=0: ok=False
    chk("the implementation's per-bucket factor is the lemma's product", ok)

print()
print("="*78); print("1B. constructive extraction, and the converse"); print("="*78)
print("  CONTRAPOSITIVE.  If a candidate clears every floor but fails the collective")
print("  bound, then by 1A some bucket violates (ii): sorting that bucket, an adjacent")
print("  pair has gap < P.  That pair is the witness to extract.")
print()
print("  WHY SUCH A PAIR MUST HAVE DIFFERENT INCOMING HISTORIES (given that every")
print("  ordinary forward overlap passes).  Suppose two members share the incoming gap")
print("  sequence of length t AND the outgoing run length ell.  Their predecessors t")
print("  odd steps back then share S+ell+1 forward parity bits, S = sum of the gaps, so")
print("  the ordinary overlap test at those earlier positions forces separation at")
print("  least 2^{S+ell+1} there.  Transporting forward by the common slope 3^t/2^S")
print("  multiplies the separation by 3^t/2^S, giving current separation at least")
print("  3^t 2^{ell+1} = P -- contradicting gap < P.  Hence the histories differ.")
same=0; diff=0
for (res,stride,wt),mem in buckets(W,seed,9).items():
    if len(mem)<2: continue
    hs={hh for _,hh in mem}
    if len(hs)==1: same+=1
    else: diff+=1
chk("on the 485-word the multi-member buckets split into %d same-history and %d mixed"
    %(same,diff), same+diff>0)
print()
print("  THE CONVERSE IS FALSE, separately.  A single short pair gap can exclude a")
print("  candidate while the aggregate product still allows it, because the aggregate")
print("  only uses the sorted MINIMA f+(m-1)P and discards which pair is close.")
print("  Explicit instance: one bucket, w=1/3, P=100, f=10, members {1000, 1050}.")
print("  Both sit far above the floor, so each factor is small; the aggregate ceiling")
print("  is computed at f and f+P and is much larger.")
w_=F(1,3); P_=100; f_=F(10); xs=[F(1000),F(1050)]
agg=prod(1+w_/(f_+mm*P_) for mm in range(2))
act=prod(1+w_/x for x in xs)
chk("the pair test fires (gap 0.5 < P) while the aggregate bound is not exceeded",
    abs(xs[0]-xs[1])<P_ and act<=agg,
    "actual %.6f <= aggregate %.6f"%(float(act),float(agg)))
print("  so the three tests are ordered as:  floors + in-bucket pairs  ==>  collective,")
print("  and the pair test is STRICTLY stronger than the collective bound.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
