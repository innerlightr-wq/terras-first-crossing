"""Independent audit of the June-to-Terras carry-residue route.
Fresh implementation; exact integer / rational arithmetic throughout."""
from fractions import Fraction as F
import math
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s" % ("PASS" if c else "FAIL", n, e))
    if not c: FAIL.append(n)

def crossing(n, kmax=200000):
    """first coefficient crossing of n: returns k,d,runs,endpoint,odd positions"""
    x,k,d,runs,prev = n,0,0,0,False
    pos=[]
    while k < kmax:
        odd = x%2==1
        if odd:
            pos.append(k); d+=1
            if not prev: runs+=1
            x=(3*x+1)//2
        else: x//=2
        k+=1; prev=odd
        if 3**d < (1<<k): return k,d,runs,x,pos
    return None
def carry(pos,d): return sum(3**(d-1-j)*(1<<t) for j,t in enumerate(pos))

print("="*78); print("1.  the June record seed, recomputed from scratch"); print("="*78)
n=63728127
k,d,runs,x,pos=crossing(n)
C=carry(pos,d); A=(1<<k)-3**d
chk("k=376, d=237, 97 odd runs, endpoint 49521089",
    (k,d,runs,x)==(376,237,97,49521089), "(%d,%d,%d,%d)"%(k,d,runs,x))
chk("affine identity 2^k x = 3^d n + C_w", (1<<k)*x == 3**d*n + C)
rho=F(3**d,1<<k)
chk("rho_k = 0.777067953...", abs(float(rho)-0.777067953)<1e-9, "rho = %.9f"%float(rho))
chk("it DESCENDS: endpoint < seed", x < n, "%d < %d"%(x,n))
chk("it passes the 82-run requirement (97 >= 82)", runs >= 82)
print("    So many runs alone does NOT make a counterexample: this word clears the")
print("    82-run bar but its crossing is deep (rho = 0.777, not ~1), and the descent")
print("    margin is huge.  M_w = C/3^d = %.4f while n(1-rho)/rho = %.4e."
      % (float(F(C,3**d)), float(n*(1-rho)/rho)))
chk("failure would need M_w > n(1-rho)/rho, which is false here",
    F(C,3**d) < n*(1-rho)/rho)

print()
print("="*78); print("2.  the endpoint-displacement identity and the residue test"); print("="*78)
print("  IDENTITY.  h := T^k(n) - n satisfies  2^k h = C_w - A n,  A = 2^k - 3^d > 0.")
print("  A is ODD, so for any divisor F of A:  h = C_w 2^{-k}  (mod F).")
print("  BOUNDS.  h >= 0 at a failure, and n >= 2 gives h <= floor((C_w - 2A)/2^k).")
print("  GCD TEST.  g := gcd(C_w, A) is odd and divides 2^k h, hence g | h.  If g >= s")
print("  (s = number of odd runs, and h < s at a failure) the only candidate is h = 0,")
print("  i.e. n = C_w/A; so if C_w/A is not an integer, every eligible seed descends.")
bad=[]
for m in range(3,20001,2):
    r=crossing(m)
    if not r: continue
    kk,dd,rr,xx,pp=r
    CC=carry(pp,dd); AA=(1<<kk)-3**dd
    hh=xx-m
    if (1<<kk)*hh != CC-AA*m: bad.append((m,"identity"))
    if AA%2==0: bad.append((m,"A even"))
    g=math.gcd(CC,AA)
    if hh>=0 and hh%g!=0: bad.append((m,"g|h"))
    if hh>=0 and hh > (CC-2*AA)//(1<<kk): bad.append((m,"ceiling"))
chk("identity, A odd, g|h and the n>=2 ceiling hold for every odd n < 20000", not bad,
    str(bad[:3]))

print()
print("="*78); print("3.  THE STRUCTURAL POINT: at F = A the test is COMPLETE"); print("="*78)
print("  Take F = A itself (no factoring needed).  Then")
print("     a failing eligible seed exists  <=>  h0 := C_w 2^{-k} mod A  <=  ceiling,")
print("  where ceiling = floor((C_w - 2A)/2^k).  Proof: n = (C_w - 2^k h)/A is forced,")
print("  it is automatically = Q_w mod 2^k because A is invertible mod 2^k (so it really")
print("  does realise w), h >= 0 and n >= 2 are exactly 0 <= h <= ceiling, and h must lie")
print("  in the class h0 mod A.  Writing C_w = A m + 2^k h0 gives m = n, so")
print("     ceiling >= h0  <=>  n >= 2 ,")
print("  i.e. the criterion is r(w) <= Q_w restated.  THE F = A TEST IS THE DESCENT")
print("  CRITERION IN NEW COORDINATES, not a new restriction.")
bad=[];tested=0
for m in range(3,40001,2):
    r=crossing(m)
    if not r: continue
    kk,dd,rr,xx,pp=r
    CC=carry(pp,dd); AA=(1<<kk)-3**dd
    h0=(CC*pow(1<<kk,-1,AA))%AA if AA>1 else 0
    ceil_=(CC-2*AA)//(1<<kk)
    # direct criterion: least eligible seed r(w) vs Q_w
    rw=None
    cls=CC*pow(AA,-1,1<<kk)%(1<<kk)   # = Q_w mod 2^k
    rw = cls if cls>=2 else cls+(1<<kk)
    direct = (F(rw) <= F(CC,AA))
    resid  = (h0 <= ceil_)
    if direct != resid: bad.append((m,direct,resid,h0,ceil_))
    tested+=1
chk("the F=A residue test agrees with 'r(w) <= Q_w' on every odd n < 40000",
    not bad, "%d words tested; %s"%(tested,str(bad[:2])))
print("  => the arithmetic lead is a REFORMULATION.  Its value is computational: a small")
print("     divisor F gives a cheap NECESSARY condition without computing r(w) or C_w/A.")

print()
print("="*78); print("4.  the gcd test's practical reach"); print("="*78)
hits=0; tot=0; gs=[]
for m in range(3,60001,2):
    r=crossing(m)
    if not r: continue
    kk,dd,rr,xx,pp=r
    CC=carry(pp,dd); AA=(1<<kk)-3**dd
    g=math.gcd(CC,AA); gs.append(g); tot+=1
    if g>=rr: hits+=1
print("  over %d first-crossing words with odd n < 60000:" % tot)
print("     gcd(C_w,A) = 1 in %.2f%% of words; max gcd = %d; g >= s in %d words (%.3f%%)"
      % (100.0*gs.count(1)/tot, max(gs), hits, 100.0*hits/tot))
chk("the gcd test is almost always vacuous (g is typically 1)", gs.count(1)/tot > 0.5)
print("  so the useful form is the PARTIAL filter h = C_w 2^{-k} mod F for a known")
print("  small factor F of A, not the gcd test.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
