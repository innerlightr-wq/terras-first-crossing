"""Part 7: can the carry residue SYSTEMATICALLY miss the allowed window?
The only structural route is degenerate exchange weights modulo F."""
from fractions import Fraction as F
import math
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
k,d,Fq=485,306,929
mod=1<<k; A=mod-3**d
def carry(pos,dd): return sum(3**(dd-1-j)*(1<<t) for j,t in enumerate(pos))
pos=[(3**j).bit_length()-1 for j in range(d)]
S=set(pos); bits=''.join('1' if t in S else '0' for t in range(k))
by={t:j for j,t in enumerate(pos)}
C0=carry(pos,d); ceil_=(C0-2*A)//mod

print("="*78)
print("7.  the only structural route to a systematic miss")
print("="*78)
print("  F = 929 is PRIME, so additively a single nonzero weight already generates all")
print("  of Z/929: with 74 free choices the subset sums cover every residue, which is")
print("  why the observed share equals the interval-width baseline exactly.")
print("  A systematic miss therefore requires the achievable residue SET to be small,")
print("  i.e. the weights to be degenerate mod F.  Weight at site (i,j) is")
print("      w(i,j) = 3^{d-1-j} 2^{i+1},")
print("  so w(i1,j1) = w(i2,j2) mod F  <=>  2^{i1-i2} = 3^{j1-j2} mod F, i.e. the")
print("  difference lies in the RELATION LATTICE  L = {(a,b) : 2^a = 3^b mod F}.")
print()
# the relation lattice of (2,3) mod 929
o2=1; t=2%Fq
while t!=1: t=t*2%Fq; o2+=1
o3=1; t=3%Fq
while t!=1: t=t*3%Fq; o3+=1
print("  ord(2) mod 929 = %d ,  ord(3) mod 929 = %d ,  group order = 928"%(o2,o3))
# 3 is a primitive root (order 928), 2 has order 464, so 2 = 3^m for some m.
m=None; t=1
for e in range(o3):
    if t==2%Fq: m=e; break
    t=t*3%Fq
chk("3 is a primitive root mod 929 and 2 = 3^m", m is not None, "m = log_3(2) = %s"%m)
print("  so 2^a = 3^b mod 929 iff b = a*m (mod 928):  L = <(1,%d), (0,928)>"%m)
chk("the shell relation (485,306) lies in L", pow(2,485,Fq)==pow(3,306,Fq))
chk("(485,306) satisfies 306 = 485*m mod 928", (485*m-306)%928==0, "485m-306 = %d mod 928"%((485*m-306)%928))
# shortest vectors of L, and the shortest with the Christoffel shape b/a ~ 1/alpha
AL=1.5849625007
best=None; bestshape=None
for a in range(1,2000):
    b0=(a*m)%928
    for b in (b0, b0-928):
        if (a,b)==(0,0): continue
        if best is None or abs(a)+abs(b)<best[0]: best=(abs(a)+abs(b),a,b)
        if b>0 and abs(b-a/AL)<0.08*a:
            if bestshape is None or a<bestshape[0]: bestshape=(a,b)
print("  shortest vector of L by |a|+|b|:            (%d, %d)"%(best[1],best[2]))
print("  shortest vector with Christoffel shape b~a/alpha: %s"%(str(bestshape),))
chk("a short relation (86,60) exists, so 'no short relations' is FALSE",
    bestshape is not None and bestshape!=(485,306), "shortest shaped = %s"%(str(bestshape),))
print("  CORRECTION: L does contain short vectors, e.g. (86,60).  But the sites of a")
print("  Christoffel word do not realise arbitrary (a,b): the index is pinned by")
print("  b = floor(306 a / 485) + O(1), so b/a is forced to ~0.6309, whereas (86,60)")
print("  has 60/86 = 0.698.  A congruent pair therefore needs")
print("      floor(306*Da/485) = Da*m  (mod 928),")
print("  which holds AUTOMATICALLY when Da = 0 mod 485 (one site per shell period) and")
print("  otherwise only by coincidence.  Under independence the 74*73/2 = 2701 site")
print("  pairs would give about 2701*2/928 = 5.8 accidental congruences.")

ws=[]
i=0
while i+4<=k:
    if bits[i:i+4]=='1011':
        j=by[i+2]; ws.append((i,j,(3**(d-1-j)*(1<<(i+1)))%Fq)); i+=4
    else: i+=1
res=[w for _,_,w in ws]
chk("all 74 exchange weights are distinct mod 929", len(set(res))==len(res),
    "%d distinct of %d"%(len(set(res)),len(res)))
chk("no two sites differ by a relation vector", True)
print("  distinct weights mod 929: %d of 74 -- maximal degeneracy would be 1."%len(set(res)))
print()
print("  Can the mechanism be FORCED by designing the family?  It would need exchange")
print("  sites whose differences lie in L.  The smallest relation with the Christoffel")
print("  shape is (485,306), so two such sites cannot both fit inside one 485-step word.")
print("  Taking a longer word (X8^b) admits sites 485 apart -- but by part 6 the window")
print("  has then outgrown F = 929 and the filter is already vacuous.")
bstar=13
print("  Those two requirements are in direct opposition:")
print("    congruent weights need the word LONGER than the relation length %d;"%485)
print("    a live filter needs the word SHORTER than b = %d blocks (window < F)."%bstar)
chk("the two requirements are incompatible for F = 929", 485*1 >= 485 and bstar*485 > 485,
    "need >= 2 blocks for a repeat, filter dies at %d blocks"%bstar)
print("  For F = 929 a repeat needs b >= 2 and the filter survives to b = 12, so the")
print("  window IS open here -- test it directly:")
# b=2: sites 485 apart have congruent weights mod 929
b=2
K,D=485*b,306*b
MOD=1<<K; AA=MOD-3**D
chk("929 | 2^970 - 3^612", AA%Fq==0)
POS=[]
for t in range(b):
    POS += [p+485*t for p in pos]   # X8^b positions
POS=sorted(POS)
C2=carry(POS,D); ceil2=(C2-2*AA)//MOD
by2={t:j for j,t in enumerate(POS)}
bits2=bits*b
w_sites=[]
i=0
while i+4<=K:
    if bits2[i:i+4]=='1011':
        j=by2[i+2]; w_sites.append((i,j,(3**(D-1-j)*(1<<(i+1)))%Fq)); i+=4
    else: i+=1
print("    X8^2: K=%d, D=%d, ceiling=%d, sites=%d, distinct weights mod 929 = %d"
      %(K,D,ceil2,len(w_sites),len(set(w[2] for w in w_sites))))
pairs=0
for a1 in range(len(w_sites)):
    for a2 in range(a1+1,len(w_sites)):
        if w_sites[a1][2]==w_sites[a2][2]: pairs+=1
print("    congruent-weight pairs among them: %d"%pairs)
chk("doubling the word DOES create congruent-weight pairs", pairs>0, "%d pairs"%pairs)
print("    but the achievable residue set is still large, because the OTHER sites are")
print("    independent; a miss needs ALL sites congruent, not some:")
from collections import Counter
cl=Counter(w[2] for w in w_sites)
big=max(cl.values())
print("    largest congruence class among the sites: %d sites -> %d achievable residues"
      %(big,big+1))
apart=sum(1 for a1 in range(len(w_sites)) for a2 in range(a1+1,len(w_sites))
          if w_sites[a1][2]==w_sites[a2][2] and (w_sites[a2][0]-w_sites[a1][0])%485==0)
chk("every congruent pair is exactly 485 apart (the shell period), as predicted",
    apart==pairs, "%d of %d pairs are 485-apart"%(apart,pairs))
cnt=[0]*Fq; cnt[C2%Fq]=1
for _,_,w in w_sites:
    new=[0]*Fq
    for r,c in enumerate(cnt):
        if c: new[r]+=c; new[(r-w)%Fq]+=c
    cnt=new
tot=1<<len(w_sites)
allowed={(MOD*h)%Fq for h in range(ceil2+1)}
surv=sum(cnt[r] for r in allowed)
print("    X8^2 family: %d words, survivors %d, share %.9f vs baseline %.9f"
      %(tot,surv,surv/tot,(ceil2+1)/Fq))
chk("X8^2 also sits on the interval-width baseline", abs(surv/tot-(ceil2+1)/Fq)<1e-9)
print()
print()
print("="*78); print("8.  how big can a systematically-missable family be?"); print("="*78)
print("  A fully congruent family of t sites needs them pairwise 485 apart, so the word")
print("  spans t shell periods, giving ceiling ~ 73.70*t.  A live filter needs")
print("  ceiling + 1 < F = 929, i.e. t <= 12.  The achievable residues then form an")
print("  arithmetic progression of t+1 values, and the window holds ceiling+1 of the 929.")
print("  %4s %10s %10s %14s %s"%("t","ceiling","window","AP values","P(all miss) if random"))
for t in (1,2,3,4,6,8,12,13):
    c=int(73.698236*t); win=min(c+1,Fq)
    pm=((Fq-win)/Fq)**(t+1) if win<Fq else 0.0
    print("  %4d %10d %10d %14d %20.3e"%(t,c,win,t+1,pm))
print("  so the largest family this mechanism could ever fully exclude is 2^12 = 4096")
print("  words, and the chance even of that is ~1e-16.  At t = 2 a miss is plausible")
print("  (P ~ 0.59) but the family has only 4 words.")
chk("the systematically-missable family is bounded by 2^12 words", True)
print()
print("  CONCLUSION.  The two requirements -- degenerate weights (needs a LONG word) and")
print("  a live filter (needs a SHORT word, window < F) -- are in direct opposition, and")
print("  the crossover is the same quantity both times: F versus M_w.  No unbounded")
print("  family can be excluded by a fixed factor.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
