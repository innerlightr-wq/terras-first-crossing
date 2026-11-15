"""Independent audit of the three shortest-failure results."""
from fractions import Fraction as F
from decimal import Decimal, getcontext
import math
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
def lg(x,terms=200):
    t=F(x-1,x+1); lo=2*sum(t**(2*i+1)/(2*i+1) for i in range(terms))
    return lo, lo+2*t**(2*terms+1)/((2*terms+1)*(1-t*t))
LN2,LN3=lg(2),lg(3)
B=28*10**18
def affine(bits):
    C=d=0
    for j,b in enumerate(bits):
        if b=='1': C=3*C+(1<<j); d+=1
    return C,d
def is_first_crossing(bits):
    oc=0
    for t,b in enumerate(bits,1):
        oc+=int(b)
        if t<len(bits) and 3**oc<(1<<t): return False,t
    C,d=affine(bits)
    return (3**d<(1<<len(bits))), None

print("="*78); print("1.  a shallow crossing does NOT force small memory"); print("="*78)
k,d=485,306
pos={(3**i).bit_length()-1 for i in range(d)}
mech=''.join('1' if j in pos else '0' for j in range(k))
runs=k-d
lowm='1'*(d-runs+1)+'01'*(runs-1)+'0'
for nm,w in (("critical-line",mech),("long initial rise",lowm)):
    ok,bad=is_first_crossing(w)
    C,dd=affine(w)
    r=1+w.count('01')
    print("  %-20s k=%d d=%d runs=%d rho=%.9f M_w=%.9f  first-crossing: %s"
          %(nm,len(w),dd,r,float(F(3**dd,1<<len(w))),float(F(C,3**dd)),ok))
    chk("%s is a genuine first-crossing word on the same shell with %d runs"%(nm,runs),
        ok and len(w)==k and dd==d and r==runs)
C1,_=affine(mech); C2,_=affine(lowm)
chk("same k, d, runs, rho but M_w = 73.698236 vs 2.001023",
    abs(float(F(C1,3**d))-73.698236116396)<1e-9 and abs(float(F(C2,3**d))-2.001022761796)<1e-9)
print("  INFINITE VERSION.  Put S_j = floor(j*alpha) and cross at k = ceil(d*alpha).")
print("  Then R_j = -{j*alpha} and M_w = (1/3) sum 2^{-{j alpha}} > d/6, since each")
print("  2^{-{j alpha}} > 1/2.  So memory grows LINEARLY in d while the crossing depth")
print("  1 - rho_k = 1 - 2^{-(k-alpha d)} can be made arbitrarily small.")
bad=[]
for dd in (100,306,1000,5000):
    P=[(3**i).bit_length()-1 for i in range(dd)]
    kk=max(P)+1
    while 3**dd>=(1<<kk): kk+=1
    w=''.join('1' if t in set(P) else '0' for t in range(kk))
    ok,_=is_first_crossing(w)
    CC,_=affine(w); M=F(CC,3**dd)
    if not (ok and M>F(dd,6)): bad.append(dd)
chk("M_w > d/6 and the word is a genuine first crossing, for d = 100..5000", not bad, str(bad))
print("  => ANY argument forcing small memory must use the seed or the failure itself.")
print("     Endpoint depth and run count cannot do it.  [PROVED]")

print()
print("="*78); print("2.  the uniform low-memory range"); print("="*78)
P,Q=630118245525664765,397560349370386783
p,q=6724555128221608268,4242721909926539673
chk("adjacent fractions: P*q - p*Q = 1", P*q-p*Q==1)
chk("p/q < alpha", p*LN2[1]-q*LN3[0] < 0)
chk("alpha < P/Q", P*LN2[0]-Q*LN3[1] > 0)
print("  With (k,d) = u(P,Q) + v(p,q), u = qk - pd and v = Pd - Qk.  Lambda > 0 gives")
print("  k/d > alpha > p/q, so u >= 1.  And d = uQ + vq with u >= 1, so v >= 1 would")
print("  force d >= Q + q.  Hence d < Q + q implies v <= 0, and since")
print("  P ln2 - Q ln3 > 0 while p ln2 - q ln3 < 0,")
print("     Lambda = u(P ln2 - Q ln3) + v(p ln2 - q ln3) >= P ln2 - Q ln3 .")
lowb=P*LN2[0]-Q*LN3[1]
chk("B*(P ln2 - Q ln3) > 4.24722566", B*lowb > F(424722566,10**8),
    "B*Lambda >= %.8f"%float(B*lowb))
chk("the cutoff is d < Q + q = 4640282259296926456", Q+q==4640282259296926456)
nx_lo=(P+p)*LN2[0]-(Q+q)*LN3[1]; nx_hi=(P+p)*LN2[1]-(Q+q)*LN3[0]
chk("at the mediant the bound drops to 0.75511 < B*Lambda < 0.75512",
    F(75511,10**5) < B*nx_lo <= B*nx_hi < F(75512,10**5),
    "B*Lambda in (%.5f, %.5f)"%(float(B*nx_lo),float(B*nx_hi)))
print("  => for d < 4.640e18 every failure has M_w > 4.24722566, so any word with")
print("     M_w <= 4 descends there regardless of run count.  [PROVED]")
print("     This supersedes the 2.676 figure I quoted last round, which came from the")
print("     d < 2^61 split and was not uniform.")

print()
print("="*78); print("3.  the interior-comparison lemma"); print("="*78)
print("  LEMMA.  Let n be a counterexample of LEAST first-crossing length k, x_t = T^t(n).")
print("  Then n <= x_k < min(x_1,...,x_{k-1}).")
print("  PROOF.  Let x_j be the least interior value (j >= 1; x_j >= n >= 2 since no")
print("  prefix has crossed).  Relative to x_j the coefficient after k-j further steps")
print("  is rho_k/rho_j < 1, because rho_k < 1 <= rho_j.  So x_j has its own first")
print("  crossing at some t* <= k-j < k.  By minimality of k, x_j is not a")
print("  counterexample, so x_{j+t*} < x_j.  If 1 <= j+t* <= k-1 that contradicts the")
print("  minimality of x_j, so j+t* = k and x_k < x_j.  QED")
# verify the lemma's engine on real orbits: every seed whose crossing is shorter descends
bad=[]
for n0 in range(3,30001,2):
    x,t,dd=n0,0,0
    vals=[n0]
    while t<100000:
        b=x%2; t+=1; dd+=b
        x=(3*x+1)//2 if b else x//2
        vals.append(x)
        if 3**dd<(1<<t): break
    kk=t
    if vals[kk]>=n0:      # a failure (none expected in this range)
        if kk>=2 and vals[kk]>=min(vals[1:kk]): bad.append(n0)
chk("no failure in range violates the lemma (vacuously: none exist below 2.8e19)",
    not bad, str(bad[:3]))
print()
print("  THE ELIMINATION.  With p = 3^{d_j}, q = 2^j, P = 3^d, Q = 2^k, A = Q-P:")
print("     q x_j = p n + C_j ,   n = (C - Q h)/A ,   x_j > x_k = n + h")
print("  gives   h < U := ((p-q)C + A C_j)/(pQ - qP),  and pQ - qP > 0 since")
print("  rho_k/rho_j < 1.  Dividing through by qQ:")
print("     U = (delta_j rho_k M_w + eps rho_j M_j)/(delta_j + eps),")
print("  delta_j = rho_j - 1 >= 0,  eps = 1 - rho_k > 0,  M_j = C_j/3^{d_j}.")
print("  So U is a WEIGHTED AVERAGE of M_w and M_j with weights delta_j and eps.")
def U_of(bits,j):
    C,dd=affine(bits); Cj,dj=affine(bits[:j])
    pp,qq=3**dj,1<<j; PP,QQ=3**dd,1<<len(bits); AA=QQ-PP
    return F((pp-qq)*C+AA*Cj, pp*QQ-qq*PP), F(Cj,pp), F(C,PP), F(pp,qq)-1, 1-F(PP,QQ)
def perturb(w,sites):
    r=list(w)
    for i in sites:
        assert ''.join(r[i:i+4])=='1011'
        r[i:i+4]=list('1101')
    return ''.join(r)
probe=perturb(mech,[101,410])
ok,_=is_first_crossing(probe)
chk("the swapped word (sites 101, 410) is a genuine first-crossing word", ok)
Cp,_=affine(probe)
chk("its M_w = 73.377595728843", abs(float(F(Cp,3**d))-73.377595728843)<1e-9)
Uv,Mj,Mw,dj_,eps=U_of(probe,84)
chk("the eliminated bound at j=84 is U = 53.452393...", abs(float(Uv)-53.452393395143)<1e-9,
    "U = %.12f"%float(Uv))
alt=(dj_*F(3**d,1<<k)*Mw + eps*(1+dj_)*Mj)/(dj_+eps)
chk("the weighted-average identity reproduces U exactly", alt==Uv)
print("    at j=84:  delta_j = %.9f,  eps = %.9f,  M_j = %.6f,  M_w = %.6f -> U = %.6f"
      %(float(dj_),float(eps),float(Mj),float(Mw),float(Uv)))
A=(1<<k)-3**d
chk("929 | A", A%929==0)
h0=Cp*pow(1<<k,-1,929)%929
cap=(Cp-2*A)//(1<<k)
chk("the factor test leaves exactly h = 66, and the cap is 73", h0==66 and cap==73,
    "h = %d mod 929, cap %d"%(h0,cap))
chk("66 > U = 53.45, so the two conditions are INCOMPATIBLE: word excluded",
    F(66)>Uv)
print("  => the interior comparison settles a word the factor test alone could not.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
