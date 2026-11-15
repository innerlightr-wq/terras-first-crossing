"""Audit: the generalized restart lemma, the 1455-step 1/6 cutoff, the
no-low-return family, and the overlap-vs-packing question."""
from fractions import Fraction as F
from decimal import Decimal, getcontext
import math
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
getcontext().prec=200
AL=Decimal(3).ln()/Decimal(2).ln()

print("="*78); print("0. CORRECTION: my six-shell list was mis-transcribed"); print("="*78)
print("  My loop produced (485b, 306b) for b = 1..6.  The report printed an arithmetic")
print("  progression (3038,1917), (5591,3528), ... which is NOT what the loop returned.")
bad=[(3038,1917),(5591,3528),(8144,5139),(10697,6750),(13250,8361)]
chk("all five mis-transcribed pairs have 2^k < 3^d (not crossings)",
    all(2**k<3**d for k,d in bad))
chk("(3038,1917) is 282 mod 929, not 0", (pow(2,3038,929)-pow(3,1917,929))%929==282)
chk("the loop's actual output (485b,306b), b=1..6, is valid and divisible",
    all(2**(485*b)>3**(306*b) and (pow(2,485*b,929)-pow(3,306*b,929))%929==0
        for b in range(1,7)))
print("  => a transcription error in the prose, not in the code.  The valid shells are")
print("     exactly the reviewer's replacements.")

print()
print("="*78); print("1. the generalized interior-restart lemma"); print("="*78)
print("  LEMMA.  Let n be a counterexample of least first-crossing length k, x_i = T^i(n),")
print("  rho_i = 3^{d_i}/2^i.  For an interior i, let e(i) be the FIRST later index with")
print("  rho_{e(i)} < rho_i; it exists since rho_k < 1 <= rho_i.  Then x_{e(i)} < x_i.")
print("  PROOF.  Measured from x_i the coefficient after t steps is rho_{i+t}/rho_i, so")
print("  x_i's own first crossing is at t = e(i) - i <= k - i < k.  By minimality of k,")
print("  x_i is not a counterexample, hence it descends there: x_{e(i)} < x_i.  QED")
print("  This is STRICTLY stronger than the earlier version, which only compared the")
print("  global minimum interior value with x_k: here e(i) may be < k, so it compares")
print("  two INTERIOR values that no endpoint test ever relates.")
def affine(bits):
    C=d=0
    for i,b in enumerate(bits):
        if b=='1': C=3*C+(1<<i); d+=1
    return C,d
def prefixes(bits):
    """(d_i, C_i) for every i"""
    out=[(0,0)]; C=d=0
    for i,b in enumerate(bits):
        if b=='1': C=3*C+(1<<i); d+=1
        out.append((d,C))
    return out
def restart_bound(bits,i,e,C,Q,A):
    pre=prefixes(bits)
    di,Ci=pre[i]; de,Ce=pre[e]
    pi,qi=3**di,1<<i; pe,qe=3**de,1<<e
    den=qi*pe-qe*pi                       # < 0 since rho_e < rho_i
    assert den<0
    V=F(qe*Ci-qi*Ce, den)                 # n > V
    return F(C-A*V, Q)                    # h < this

print()
print("="*78); print("2. the fresh 1455-step mechanical word"); print("="*78)
k,d=1455,918
pos=set((3**a).bit_length()-1 for a in range(d))
w=''.join('1' if t in pos else '0' for t in range(k))
C,dd=affine(w); Q=1<<k; P=3**d; A=Q-P
chk("(1455,918) is a critical shell with 929 | A", dd==d and A>0 and A%929==0)
oc=0; ok=True
for t,b in enumerate(w,1):
    oc+=int(b)
    if t<k and 3**oc<(1<<t): ok=False; break
chk("it is a genuine first-crossing word", ok and 3**d<Q)
cap=(C-2*A)//Q
h0=C*pow(Q,-1,929)%929
chk("the factor 929 leaves h = 13 as the sole candidate under the cap",
    h0==13 and h0<=cap, "h = %d mod 929, cap = %d"%(h0,cap))
pre=prefixes(w)
rho=[F(3**pre[t][0],1<<t) for t in range(k+1)]
# every comparison with the FINAL endpoint
allow=0; tested=0
for i in range(1,k):
    di,Ci=pre[i]; pi,qi=3**di,1<<i
    den=qi*P-Q*pi
    if den>=0: continue
    V=F(Q*Ci-qi*C, den)
    tested+=1
    if F(13) < F(C-A*V,Q): allow+=1
chk("all %d comparisons with the final endpoint ALLOW h = 13"%tested, allow==tested,
    "%d of %d allow"%(allow,tested))
# the restart from i=84
i=84
e=next(t for t in range(i+1,k+1) if rho[t]<rho[i])
chk("the first relative crossing from step 84 is at step 569", e==569, "e(84) = %d"%e)
U=restart_bound(w,i,e,C,Q,A)
chk("it gives exactly h < 1/6", U==F(1,6), "bound = %s"%U)
chk("so h = 13 is EXCLUDED by this single interior restart", F(13)>=U)
print("    (k=1455 is below the published lower bound and seed integrality is relaxed,")
print("     so this calibrates the METHOD.)")

print()
print("="*78); print("3. the family that avoids low returns"); print("="*78)
print("  S_i = max{i, floor(i*alpha) - 1}, k = ceil(d*alpha).")
print("  %6s %7s %9s %16s %14s %s"%("d","k","runs","min proper rho","M_w","M_w > d/12"))
good=True
for dq in (306,612,918,1530):
    S=[max(i,int(AL*i)-1) for i in range(dq)]
    if len(set(S))!=dq: good=False; continue
    kk=int(AL*dq)+1
    if max(S)>=kk: good=False; continue
    ww=''.join('1' if t in set(S) else '0' for t in range(kk))
    CC,_=affine(ww)
    pr=prefixes(ww)
    mn=min(F(3**pr[t][0],1<<t) for t in range(1,kk))
    Mw=F(CC,3**dq)
    print("  %6d %7d %9d %16s %14.4f %s"
          %(dq,kk,1+ww.count('01'),mn,float(Mw),Mw>F(dq,12)))
    if mn<F(3,2) or Mw<=F(dq,12): good=False
chk("every proper nonempty prefix coefficient is >= 3/2 and M_w > d/12", good)
print("  => a word can keep EVERY interior coefficient far above 1 while the final")
print("     crossing is arbitrarily shallow and memory grows linearly.  So 'failure")
print("     forces a low return' needs an arithmetic exclusion of this whole class;")
print("     quantifying over a class is not by itself an escape from the realizer")
print("     inequality.  [the reviewer's point, CONFIRMED]")
print("  Note the restart lemma still applies here: it uses a drop from one HIGH")
print("  coefficient to another, both far above 1 -- which the endpoint test cannot see.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
