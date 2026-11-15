"""Part 3 core: recover the 13 cuts exactly and examine their structure."""
from fractions import Fraction as F
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
def logb(base,terms=200):
    x=F(base-1,base+1)
    v=sum(2*x**(2*i+1)/(2*i+1) for i in range(terms))
    return v, v+2*x**(2*terms+1)/((2*terms+1)*(1-x*x))
lo2,hi2=logb(2); lo3,hi3=logb(3)
k=1412737287532845639425
d=891337988684307508298
Fq=929
print("="*78); print("3A. the record shell and its 13 cuts, recovered exactly"); print("="*78)
chk("2^k = 3^d = 608 mod 929", pow(2,k,Fq)==pow(3,d,Fq)==608)
chk("gcd(k,d)=1", __import__('math').gcd(k,d)==1)
lowL=k*lo2-d*hi3; highL=k*hi2-d*lo3
chk("0 < Lambda, and d*Lambda/ln2 in (130.5,130.6)",
    lowL>0 and F('130.5')<d*lowL/hi2 and d*highL/lo2<F('130.6'))
B=28*10**18
chk("B*Lambda in (2.84228, 2.84229)", F('2.84228')<B*lowL and B*highL/(1-highL)<F('2.84229'))
inv=pow(k,-1,d)
DETS=[16,23,39,46,53,69,76,83,99,106,113,129,136]
cuts=[]
for D in DETS:
    a=D*inv%d
    j=(k*a-D)//d
    assert 0<a<d and 0<j<k and k*a-d*j==D
    cuts.append((D,j,a))
print("  %5s %26s %26s %10s"%("D","j_D","a_D","j_D - a_D"))
for D,j,a in cuts:
    print("  %5d %26d %26d %10d"%(D,j,a,j-a))
chk("all 13 cuts recovered with k*a_D - d*j_D = D exactly",
    all(k*a-d*j==D for D,j,a in cuts))
print()
print("  A word AVOIDS cut D when its odd count at shortcut step j_D differs from a_D.")
print("  So a family avoids all 13 iff at each of these 13 positions the running odd")
print("  count misses one prescribed value.")
print()
print("="*78); print("3B. can a block family avoid all 13 cuts?"); print("="*78)
print("  Take interchangeable blocks of common length L and common odd count C.  Then")
print("  the odd count at any BLOCK BOUNDARY is C*(boundary index), independent of the")
print("  order; inside a block it varies with the order.  Write j_D = L*b_D + r_D.")
print("  If r_D = 0 the count is forced to C*b_D and avoidance needs C*b_D != a_D.")
print("  If r_D > 0 the count lies in [C*b_D, C*b_D + min(C,r_D)] and the order can")
print("  usually be chosen to miss a_D -- unless the interval is the single point a_D.")
print("  Test candidate block sizes:")
print("  %8s %8s %10s %10s %s"%("L","C","r_D=0 hits","forced hits","verdict"))
best=None
for L in (2,3,4,5,6,7,8,9,10,11,12,13,16,24,32,48,64,97,128):
    C=None
    # choose C = round(L/alpha) so blocks are near-critical
    for cc in range(1,L):
        if abs(F(cc)*hi3-F(L)*lo2)<F(1,2) or abs(F(cc)*lo3-F(L)*hi2)<F(1,2):
            C=cc; break
    if C is None: continue
    zero=0; forced=0
    for D,j,a in cuts:
        b,r=divmod(j,L)
        if r==0:
            zero+=1
            if C*b==a: forced+=1
        else:
            lo_=C*b; hi_=C*b+min(C,r)
            if lo_==a==hi_: forced+=1
    v="AVOIDABLE" if forced==0 else "%d forced collisions"%forced
    print("  %8d %8d %10d %10d %s"%(L,C,zero,forced,v))
    if forced==0 and best is None: best=(L,C)
chk("a block length with no forced collision at any of the 13 cuts exists",
    best is not None, "first is (L,C) = %s"%(str(best),))
if best:
    L,C=best
    print("  => with blocks of length %d and odd count %d, none of the 13 cuts is forced;"%(L,C))
    print("     every cut either misses a_D at the boundary or has order freedom inside.")
    ok=all((lambda b,r: (r==0 and C*b!=a) or (r>0 and not (C*b==a==C*b+min(C,r))))(*divmod(j,L))
           for D,j,a in cuts)
    chk("re-checked: all 13 cuts avoidable at (L,C) = (%d,%d)"%(L,C), ok)
    print("  block count m = k/L = %.4e, odd total C*m = %.4e vs d = %.4e"
          %(k/L, C*(k//L), d))
    print()
    print("  BUT a UNIFORM FULL TILING IS IMPOSSIBLE on this shell.  m blocks of (L,C)")
    print("  covering all k steps need C/L = d/k exactly, and gcd(k,d) = 1 forces L = k,")
    print("  C = d -- a single block, with no interchange freedom at all.")
    chk("gcd(k,d)=1 rules out any uniform full tiling with L < k",
        __import__('math').gcd(k,d)==1 and k%L!=0 or C*(k//L)!=d)
    print("  So the blocks must cover a PREFIX only: m blocks of (L,C) over the first")
    print("  m*L steps, then a FIXED tail of k - m*L steps carrying d - m*C ones.")
    print("  Cuts inside the blocked prefix keep order freedom; cuts beyond it have")
    print("  counts fixed by the tail and must miss a_D by the tail's own structure.")
    print("  %8s %8s %14s %14s %s"%("L","C","m_max","cuts inside","cuts in tail"))
    for (L2,C2) in ((2,1),(5,3),(8,5),(13,8),(97,61)):
        m_max=min(k//L2,d//C2)
        inside=sum(1 for D,jj,aa in cuts if jj<m_max*L2)
        print("  %8d %8d %14.4e %14d %14d"%(L2,C2,m_max,inside,13-inside))
    chk("every cut lies inside the maximal blocked prefix for each tested (L,C)",
        all(all(jj<min(k//L2,d//C2)*L2 for D,jj,aa in cuts) for (L2,C2) in ((2,1),(5,3),(8,5))))
    print("  => all 13 cuts sit inside the blocked region for these block sizes, so the")
    print("     order freedom is available at every cut.  What is NOT established here:")
    print("     (i) a tail realising d - m*C ones with the first crossing exactly at k;")
    print("     (ii) >= 82 maximal odd runs together with the avoidance choices;")
    print("     (iii) a nonempty candidate range 2.8e19 < n <= Q_w after the baseline")
    print("           tests;  (iv) a uniform signed-pair cover.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
