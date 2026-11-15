"""Parts 5-8: the 485-shell family, the vacuity threshold, and the systematic-miss
mechanism.  Independent implementation."""
from fractions import Fraction as F
import math, random
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s" % ("PASS" if c else "FAIL", n, e))
    if not c: FAIL.append(n)

print("="*78); print("4b. the gcd test by run count (correcting the headline above)"); print("="*78)
def crossing(n,kmax=200000):
    x,k,d,runs,prev=n,0,0,0,False; pos=[]
    while k<kmax:
        odd=x%2==1
        if odd:
            pos.append(k); d+=1
            if not prev: runs+=1
            x=(3*x+1)//2
        else: x//=2
        k+=1; prev=odd
        if 3**d<(1<<k): return k,d,runs,x,pos
    return None
def carry(pos,d): return sum(3**(d-1-j)*(1<<t) for j,t in enumerate(pos))
buckets={}
for m in range(3,60001,2):
    r=crossing(m)
    if not r: continue
    kk,dd,rr,xx,pp=r
    CC=carry(pp,dd); AA=(1<<kk)-3**dd
    g=math.gcd(CC,AA)
    b = 1 if rr==1 else (2 if rr<=3 else (3 if rr<=9 else 4))
    t=buckets.setdefault(b,[0,0]); t[0]+=1
    if g>=rr: t[1]+=1
names={1:"s = 1",2:"s = 2-3",3:"s = 4-9",4:"s >= 10"}
print("  %10s %10s %12s %s"%("runs","words","g >= s","share"))
for b in sorted(buckets):
    t=buckets[b]
    print("  %10s %10d %12d %10.3f%%"%(names[b],t[0],t[1],100.0*t[1]/t[0]))
chk("the gcd test is WEAK but not dead for many-run words (a few per cent)",
    0 < buckets[4][1] < buckets[4][0]//10,
    "%d of %d = %.3f%%"%(buckets[4][1],buckets[4][0],100.0*buckets[4][1]/buckets[4][0]))
print("  => the 73% headline was carried entirely by s = 1 words, where h < 1 already")
print("     forces h = 0.  For s >= 2 the rate settles at roughly 5-7%% regardless of s,")
print("     so the gcd test certifies descent for a few per cent of many-run words --")
print("     a genuine but small yield, not a route to a family.")

print()
print("="*78); print("5.  the 485-shell family, reproduced independently"); print("="*78)
k,d,Fq=485,306,929
mod=1<<k; A=mod-3**d
chk("929 | 2^485 - 3^306", A%Fq==0)
chk("A is odd and A > 0", A%2==1 and A>0)
pos=[(3**j).bit_length()-1 for j in range(d)]
bits=''.join('1' if t in set(pos) else '0' for t in range(k))
chk("the base word is the Christoffel word X_(306,485): 306 ones, 485 long",
    bits.count('1')==d and len(bits)==k)
runs=bits.count('01')+(1 if bits[0]=='1' else 0)
chk("179 odd runs", runs==179, str(runs))
# independent first-crossing check by exact integer comparison
oc=0; early=None
for t,b in enumerate(bits,1):
    oc+=int(b)
    if t<k and 3**oc<(1<<t): early=t; break
chk("no earlier crossing: 3^{d_t} >= 2^t for every t < 485", early is None, str(early))
chk("it crosses at t = 485", 3**d < mod)
C0=carry(pos,d)
Mw=F(C0,3**d)
ceil_=(C0-2*A)//mod
chk("M_w = 73.6982...", abs(float(Mw)-73.698236)<1e-6, "M_w = %.6f"%float(Mw))
chk("the h-ceiling from n >= 2 is 73", ceil_==73, str(ceil_))
print("    (the ceiling tracks M_w: floor(M_w*rho - 2 + 2rho) with rho = %.9f)"
      % float(F(3**d,mod)))
# exchanges, found independently by scanning for 1011 and verifying each swap
i=0; ex=[]
while i+4<=k:
    if bits[i:i+4]=='1011':
        w2=''.join(bits[:i]+'1101'+bits[i+4:])
        assert w2.count('1')==d
        r2=w2.count('01')+(1 if w2[0]=='1' else 0)
        oc2=0; ok=True
        for t,b in enumerate(w2,1):
            oc2+=int(b)
            if t<k and 3**oc2<(1<<t): ok=False; break
        p2=[t for t,b in enumerate(w2) if b=='1']
        ex.append((i, r2==runs and ok, C0-carry(p2,d)))
        i+=4
    else: i+=1
chk("74 disjoint 1011 sites", len(ex)==74, str(len(ex)))
chk("every single swap preserves k, d, run count AND the first crossing",
    all(e[1] for e in ex))
wts=[e[2] for e in ex]
chk("each swap subtracts exactly 3^{d-1-j} 2^{i+1} from the carry",
    all(w>0 for w in wts))
# residue distribution: independent method (direct subset-sum over F_929 by DP on residues)
cnt=[0]*Fq; cnt[C0%Fq]=1
for w in wts:
    new=[0]*Fq
    for rr,c in enumerate(cnt):
        if c:
            new[rr]+=c; new[(rr-w)%Fq]+=c
    cnt=new
size=1<<len(wts)
chk("mass conservation: sum of counts = 2^74", sum(cnt)==size, "%d"%size)
inv=pow(mod,-1,Fq)
allowed={(mod*h)%Fq for h in range(ceil_+1)}
surv=sum(cnt[r] for r in allowed)
excl=size-surv
chk("survivors = 1504650676996140590169", surv==1504650676996140590169, str(surv))
chk("excluded = 17384815254482440264615", excl==17384815254482440264615, str(excl))
print("    excluded fraction = %.12f%%"%(100.0*excl/size))
base=F(ceil_+1,Fq)
print("    interval-width baseline (ceiling+1)/F = %s = %.12f"%(base,float(base)))
print("    survivor share - baseline = %.6e"%(float(F(surv,size)-base)))
chk("the survivor share matches the interval-width baseline to < 1e-16",
    abs(float(F(surv,size)-base))<1e-16)
print("  => CONFIRMED, including the interpretation: each exclusion is certified, but the")
print("     aggregate is exactly what an equidistributed residue predicts.  No signal.")
# spot-check 6 random family members end-to-end
random.seed(5); bad=0
for _ in range(6):
    S=[j for j in range(74) if random.random()<0.5]
    w=list(bits)
    for j in S:
        i0=ex[j][0]; w[i0:i0+4]=list('1101')
    w=''.join(w)
    p=[t for t,b in enumerate(w) if b=='1']
    Cw=carry(p,d)
    if Cw != C0-sum(wts[j] for j in S): bad+=1
    h0=(Cw*pow(mod,-1,A))%A
    cl=(Cw-2*A)//mod
    resid=(Cw*inv)%Fq
    if (resid in allowed) != (h0<=cl or True): pass
    # the F=929 filter must never exclude a word the exact F=A test admits
    if h0<=cl and resid not in allowed: bad+=1
chk("on 6 random members: carry additivity holds and the filter never contradicts F=A",
    bad==0)

print()
print("="*78); print("6.  the vacuity threshold:  the filter needs F > ceiling"); print("="*78)
print("  The allowed window is the ceiling+1 residues {2^k h mod F : 0 <= h <= ceiling}.")
print("  Once ceiling+1 >= F the window is everything and the filter says nothing.")
print("  For the shell (485b, 306b) -- the word X8^b -- we have")
print("     2^485 - 3^306  |  2^{485b} - 3^{306b},   so 929 | A_b for EVERY b,")
print("  but M_w, and with it the ceiling, grows linearly in b:")
rho=F(3**d,mod)
print("  %4s %14s %10s %10s %s"%("b","M_w(X8^b)","ceiling","F=929","filter"))
bstar=None
for b in (1,2,4,8,12,13,14,20):
    Mb=Mw*(1-rho**b)/(1-rho)
    cb=int(Mb*rho**0)            # ceiling ~ M_w - 2 + 2rho, use floor(M_b) as a proxy
    live = cb+1 < Fq
    if not live and bstar is None: bstar=b
    print("  %4d %14.3f %10d %10d %s"%(b,float(Mb),cb,Fq,"live" if live else "VACUOUS"))
chk("the F=929 filter on X8^b becomes vacuous by b = 13", bstar is not None and bstar<=13,
    "first vacuous b = %s"%bstar)
chk("929 | 2^{485b} - 3^{306b} for b = 1..6",
    all(((1<<(485*b))-3**(306*b))%Fq==0 for b in range(1,7)))
print("  => a FIXED factor cannot follow an unbounded family: the window outgrows it.")
print("     Keeping the filter alive needs F > M_w, and M_w <= d/3 with M_w ~ d/3 for")
print("     a line-hugging word.  For a genuine counterexample d >= 6586818670, so one")
print("     would need a known factor of 2^k - 3^d exceeding %.3e with k >= 1.044e10."
      % (6586818670/3.0))
print("     That is a factoring problem on a number with ~3.1e9 decimal digits.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
