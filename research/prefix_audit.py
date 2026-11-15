"""Independent audit: three numerical corrections, and the refutation of my
partial-quotient cap."""
from fractions import Fraction as F
from decimal import Decimal, getcontext
import math
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
def lg(x,terms=260):
    t=F(x-1,x+1); lo=2*sum(t**(2*i+1)/(2*i+1) for i in range(terms))
    return lo, lo+2*t**(2*terms+1)/((2*terms+1)*(1-t*t))
LN2,LN3=lg(2),lg(3); B=28*10**18
def affine(bits):
    C=d=0
    for i,b in enumerate(bits):
        if b=='1': C=3*C+(1<<i); d+=1
    return C,d

print("="*78); print("A. the certified constant: my 4.24722567 was NOT a lower bound"); print("="*78)
P,Q=630118245525664765,397560349370386783
lo=B*(P*LN2[0]-Q*LN3[1]); hi=B*(P*LN2[1]-Q*LN3[0])
print("  certified enclosure of B*Lambda_min:")
print("    lower  %s" % str(Decimal(lo.numerator)/Decimal(lo.denominator))[:24])
print("    upper  %s" % str(Decimal(hi.numerator)/Decimal(hi.denominator))[:24])
chk("B*Lambda > 4.24722566", lo > F(424722566,10**8))
chk("B*Lambda > 4.24722567 is FALSE -- 4.24722567 exceeds the true value",
    not (lo > F(424722567,10**8)) and hi < F(424722567,10**8))
chk("the true value begins 4.24722566869589", str(Decimal(lo.numerator)/Decimal(lo.denominator))[:16]=="4.24722566869589")
print("  => my report printed a float rounded to NEAREST and then quoted it as a bound.")
print("     The usable statement is B*Lambda > 4.24722566.  (Same failure mode as the")
print("     Arb .str() lesson: a rounded display is not a certificate.)")

print()
print("="*78); print("B. the 27% figure is example-specific"); print("="*78)
k,d=485,306
pos={(3**i).bit_length()-1 for i in range(d)}
mech=''.join('1' if j in pos else '0' for j in range(k))
runs=k-d
lowm='1'*(d-runs+1)+'01'*(runs-1)+'0'
def best_ratio(bits):
    C,dd=affine(bits); PP,QQ=3**dd,1<<len(bits); AA=QQ-PP
    best=None
    for j in range(1,len(bits)):
        Cj,dj=affine(bits[:j]); pp,qq=3**dj,1<<j
        den=pp*QQ-qq*PP
        if den<=0: continue
        U=F((pp-qq)*C+AA*Cj,den)
        if best is None or U<best: best=U
    return best, F(C,PP)
for nm,w in (("mechanical",mech),
             ("two-swap",mech[:101]+'1101'+mech[105:410]+'1101'+mech[414:]),
             ("low-memory 1^128(01)^178 0",lowm)):
    U,Mw=best_ratio(w)
    print("  %-28s min_j U = %10.4f  M_w = %10.4f  ratio %.6f"%(nm,float(U),float(Mw),float(U/Mw)))
U3,M3=best_ratio(lowm)
chk("the low-memory word's ratio is 0.997451, not ~0.728",
    abs(float(U3/M3)-0.997451)<1e-6, "ratio = %.6f"%float(U3/M3))
print("  => the 27% reduction is a property of the HIGH-memory examples only.  For a")
print("     low-memory word the interior comparison buys essentially nothing, which is")
print("     consistent: U is an average of M_w and M_j, and there M_j ~ M_w.")

print()
print("="*78); print("C. the approximation must be taken from the correct SIDE"); print("="*78)
getcontext().prec=160
AL=Decimal(3).ln()/Decimal(2).ln()
ALF=F(int(AL*Decimal(10)**140), 10**140)
x,cf=AL,[]
for _ in range(90):
    i=int(x); cf.append(i)
    if x==i: break
    x=1/(x-i)
p0,q0,p1,q1=1,0,cf[0],1
conv=[(cf[0],1)]
for t in cf[1:]:
    p0,q0,p1,q1=p1,q1,t*p1+p0,t*q1+q0
    conv.append((p1,q1))
tol=F(74,B)
print("  we need the INTERIOR drift alpha*d_j - j > 0, i.e. j/d_j < alpha (a LOWER")
print("  approximant).  Check the side of each candidate denominator:")
for qq in (130441933147714940, 267118416222671843):
    pp=[a for a,b in conv if b==qq]
    if not pp: continue
    pp=pp[0]
    drift=F(qq)*ALF-pp
    print("    q = %-21d : alpha*q - p = %+.4e  -> %s"
          %(qq,float(drift),"LOWER (usable)" if drift>0 else "UPPER (wrong side)"))
d1=[(a,b) for a,b in conv if b==130441933147714940][0]
chk("130441933147714940 is on the WRONG side (alpha*q - p < 0)",
    F(d1[1])*ALF-d1[0] < 0)
first=None
for (pp,qq) in conv:
    dr=F(qq)*ALF-pp
    if dr>0 and dr<tol and first is None: first=(qq,float(dr))
print("  earliest LOWER convergent denominator with 0 < alpha*q - p < %.3e :"%float(tol))
print("    q = %d  (drift %.4e)"%(first[0],first[1]))
chk("that denominator is 267118416222671843", first[0]==267118416222671843,
    "got %d"%first[0])
print("  => my quoted 130441933147714940 was the right magnitude but the wrong side;")
print("     it has rho_j < 1, which a proper prefix of a first-crossing word cannot have.")

print()
print("="*78); print("D. MY PARTIAL-QUOTIENT CAP IS REFUTED"); print("="*78)
print("  I argued: if d_j = q_m is a convergent and d a NEIGHBOURING intermediate")
print("  denominator then eps/delta_j ~ a_{m+1}, so the gain is capped by the partial")
print("  quotients of log2 3 (<= 55 below 1e20).")
print("  The error is the coupling.  delta_j and eps are values of the linear form at")
print("  INDEPENDENT denominators d_j <= d; nothing forces them to be consecutive")
print("  convergents.  Taking d_j excellent and d merely good makes the ratio large.")
print()
print("  Demonstration, purely Diophantine.  Pick a LOWER convergent q_m with tiny")
print("  drift g = alpha*q_m - p_m, and then any d >= q_m whose own deficit")
print("  1 - {alpha*d} lands in a window far above g:")
print("  CONSTRUCTION.  Let q_m be a LOWER convergent with drift g = alpha*q_m - p_m > 0.")
print("  Put N = floor(1/g) - R and d = N*q_m, so alpha*d = N*p_m + N*g with")
print("  N*g = 1 - (frac(1/g) + R)*g; hence k = ceil(alpha*d) = N*p_m + 1 and")
print("     eps = k - alpha*d = (frac(1/g) + R)*g ~ R*g ,   delta_j = g ,   d_j = q_m <= d.")
print("  The ratio is R, chosen freely, subject only to eps < tol.")
print("  %6s %26s %12s %12s %12s"%("R","d_j = q_m","g","eps","eps/g"))
rows=[]
for R in (10**3, 10**5, 10**7):
    for (pp,qq) in conv:
        g=F(qq)*ALF-pp
        if g<=0 or qq<10**20: continue
        if R*g >= tol/2: continue
        N=int(1/g)-R
        if N<=0: continue
        dd=N*qq
        ad=F(dd)*ALF
        kk=N*pp+1
        epsv=F(kk)-ad
        if epsv<=0 or epsv>=tol: continue
        rows.append((R,qq,float(g),float(epsv),float(epsv/g)))
        print("  %6.0e %26d %12.3e %12.3e %12.3e"%(R,qq,float(g),float(epsv),float(epsv/g)))
        break
chk("explicit pairs with eps/delta_j far above 55 exist (ratios 1e3, 1e5, 1e7)",
    bool(rows) and max(r[4] for r in rows)>10**6,
    "max ratio achieved = %.3e"%(max(r[4] for r in rows) if rows else 0))
chk("every constructed eps stays inside the counterexample tolerance",
    all(r[3]<float(tol) for r in rows))
chk("every constructed interior drift is POSITIVE (correct side)",
    all(r[2]>0 for r in rows))
print("  and the ratio is UNBOUNDED: for d_j of quality g, the deficits 1-{alpha*d}")
print("  are equidistributed over d >= d_j, so a window (R*g, tol) is hit for any R as")
print("  soon as R*g < tol.  Since g can be made arbitrarily small by taking a deeper")
print("  convergent, R is unbounded.  MY CAP CLAIM WAS WRONG.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
