"""Independent verification of the two demonstrations and the infinite construction."""
from fractions import Fraction as F
import math
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
def affine(bits):
    C=d=0
    for i,b in enumerate(bits):
        if b=='1': C=3*C+(1<<i); d+=1
    return C,d

print("="*78); print("E. the 485-shell demonstration, recounted independently"); print("="*78)
k,d,j,Fq=485,306,84,929
pos=[(3**a).bit_length()-1 for a in range(d)]
idx={t:a for a,t in enumerate(pos)}
w=''.join('1' if t in idx else '0' for t in range(k))
C0,dd=affine(w); Cj,dj=affine(w[:j])
P,Q=3**d,1<<k; p,q=3**dj,1<<j; A=Q-P
chk("the base word crosses exactly at t=485 and 929|A",
    all((3**sum(int(c) for c in w[:t])<(1<<t))==(t==k) for t in range(1,k+1)) and A%Fq==0)
print("    at the cut j=84: d_j = %d, rho_j - 1 = %.9e, eps = %.9e"
      %(dj,float(F(p,q)-1),float(1-F(P,Q))))
# sites, with the straddling one frozen
free=[];frozen=[]
t=0
while t+4<=k:
    if w[t:t+4]=='1011':
        if t+2==j: frozen.append(t)
        else: free.append((t, 3**(d-1-idx[t+2])*(1<<(t+1))))
        t+=4
    else: t+=1
chk("exactly one site straddles the cut (t=82) and 73 remain free",
    frozen==[82] and len(free)==73, "frozen %s, free %d"%(frozen,len(free)))
print("    freezing the straddling site keeps C_j, p, q fixed, so U depends on C only")
print("    through the term (p-q)C with p>=q: U is MONOTONE INCREASING in C, hence")
print("    maximal at the unswapped word.")
Umax=F((p-q)*C0+A*Cj, p*Q-q*P)
chk("uniform bound U_84 <= 53.667543674434", abs(float(Umax)-53.667543674434)<1e-9,
    "U = %.12f"%float(Umax))
Wsum=sum(x for _,x in free)
mincap=(C0-Wsum-2*A)//Q
chk("uniform endpoint cap >= 62", mincap==62, "min cap = %d"%mincap)
inv=pow(Q,-1,Fq)
cnt=[0]*Fq; cnt[C0*inv%Fq]=1
for _,x in free:
    cnt=[v+cnt[(r+x*inv)%Fq] for r,v in enumerate(cnt)]
chk("mass conservation over 2^73 words", sum(cnt)==2**73)
tot=sum(cnt[54:63])
chk("91499027655170685801 words have their candidate in 54..62", tot==91499027655170685801,
    "%d"%tot)
print("    those words pass BOTH the factor test (unique residue <= cap) and their own")
print("    endpoint cap, yet every one needs h >= 54 > U_84.  The interior comparison")
print("    excludes all %.4e of them as shortest failures."%tot)
print("    (k=485 is far below the published lower bound, so this demonstrates the")
print("     METHOD, not a new threshold.)")

print()
print("="*78); print("F. the subset-growth lemma behind the 929 argument"); print("="*78)
print("  CLAIM.  In Z/p with p prime, if S is nonempty with |S| < p and w != 0 then")
print("  |S union (S+w)| > |S|.  Proof: otherwise S+w = S, so S is a union of cosets")
print("  of <w> = Z/p, forcing S = Z/p.  Hence each of the four swaps strictly")
print("  enlarges the candidate-residue set, and four swaps give at least 5 residues.")
bad=[]
import random
random.seed(3)
for _ in range(400):
    pr=929; sz=random.randrange(1,pr)
    S=set(random.sample(range(pr),sz)); wv=random.randrange(1,pr)
    if len(S|{(x+wv)%pr for x in S})<=len(S) and len(S)<pr: bad.append((sz,wv))
chk("verified on 400 random (S,w) in Z/929", not bad, str(bad[:2]))
chk("5 residues cannot all lie in {0,1,2,3} (only 4 values)", 5>4)
print("  => at least one word of each 16-word group is excluded.  Valid existence")
print("     statement; it is not a single uniform choice across the family.")
print()
print("  and the divisibility: ord(2)=464 and ord(3)=928 mod 929, so with")
print("  D = 928u, K = 464*ceil(2*alpha*u), a = 928t, j = 464*floor(2*alpha*t):")
chk("2^464 = 1 mod 929", pow(2,464,929)==1)
chk("3^928 = 1 mod 929", pow(3,928,929)==1)
chk("hence 929 | 2^{j+K} - 3^{a+D} for any such j,K,a,D",
    all((pow(2,464*x,929)-pow(3,928*y,929))%929==0 for x in (1,3,7) for y in (1,2,5)))

print()
print("="*78); print("G. the factor-7 two-scale family"); print("="*78)
print("  reported: k=79847916192729377980, d=50378426086673132390, 178 runs,")
print("  eps/delta > 265, M in (74.370,74.371) at t=2, h < 2 with h = 3 mod 7.")
K,D=79847916192729377980,50378426086673132390
chk("2^k > 3^d (the endpoint is a crossing)", K > D*math.log2(3))
chk("7 | 2^k - 3^d", (pow(2,K,7)-pow(3,D,7))%7==0,
    "2^k mod 7 = %d, 3^d mod 7 = %d"%(pow(2,K,7),pow(3,D,7)))
chk("h < 2 leaves only h in {0,1}, and 3 mod 7 is outside it", 3 not in (0,1))
print("  => the whole family is excluded as shortest failures by one residue class.")
print("     Very long odd runs make its realizers huge, so this too demonstrates the")
print("     method rather than a new threshold.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
