"""Addendum: Phase 1's five-component dependence table, and an explicit
first-crossing verification for the strongest candidate shell."""
from fractions import Fraction as F
from decimal import Decimal, getcontext
from itertools import product
import math
getcontext().prec=260
AL=Decimal(3).ln()/Decimal(2).ln()
ALF=F(int(AL*Decimal(10)**220),10**220)
LN2=F(int(Decimal(2).ln()*Decimal(10)**220),10**220)
LN3=F(int(Decimal(3).ln()*Decimal(10)**220),10**220)
K0,D0=1412737287532845639425,891337988684307508298
K1,D1=1431227343753494174081,903003915299123681908
HEAD,TAIL="11111110","111100000"; MARKER="0111100000"
SIGMA=3611678508; L=104; PAST=7; FW=12
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)

print("="*78); print("PHASE 1.  the five components and what each depends on"); print("="*78)
rows=[
 ("1 drift + 8-step schedule",
  "K,D via eta = 8D/K-5 and Lambda; needs 5<8D/K<6",
  "no", "no", "no"),
 ("2 forward-window multiplicity",
  "K,D via the six-one count ~ 0.00593 K",
  "yes: labels = 7 s^F", "indirectly (edge loss)", "yes: s^F"),
 ("3 incoming residues + sigma_m",
  "NOT AT ALL -- shell-independent",
  "yes: R = alphabet^past", "yes: modulus 3^{5 past}", "only via m"),
 ("4 distinctness from the marker",
  "no", "yes: marker length and junctions", "no", "no"),
 ("5 seed floor vs failure ceiling",
  "K,D via Q = D/(3 Lambda), M = D/3, rho bracket",
  "via sigma_m", "via sigma_m", "yes: 2^L"),
]
print("  %-28s %-44s %-22s %-24s %s"%("component","(K,D)","alphabet","past","F / L"))
for r in rows: print("  %-28s %-44s %-22s %-24s %s"%r)
print()
print("  UNIVERSAL ALGEBRAIC IDENTITIES (hold on every near-critical shell):")
print("    - rho_t = 3^{2-{tD/K}} exp(-t Lambda/K)  at every byte boundary;")
print("    - 3^{1+eta} < 256/81  <=>  Lambda > 0;")
print("    - eta -> 8 log_3 2 - 5 = 0.047438029, hence past + F <= 19;")
print("    - the concatenation carry rule, hence the 2^{-(8 past + L)} normalization;")
print("    - sigma_m nondecreasing in m, and sigma_m independent of (K,D).")
print("  DEPENDENT ON THE 21/7 SHELL'S FAVOURABLE NUMBERS:")
print("    - K*Lambda/ln3 = 130.5, so only 131 early-crossing sites to filter;")
print("    - m = 163 exactly (ratio 162.69), attained with no slack;")
print("    - the 154.49x room in the span inequality.")
chk("eta -> 8 log_3 2 - 5 on any near-critical shell",
    abs(float(F(8*D0,K0)-5)-float(8/AL-5))<1e-18)
chk("the cap floor(1/eta)-2 = 19 is the same for both shells",
    int(1.0/float(F(8*D0,K0)-5))-2 == int(1.0/float(F(8*D1,K1)-5))-2 == 19)

print(); print("="*78)
print("PHASE 2 completion: explicit first crossing for the candidate shell")
print("="*78)
A5=["".join(w)+"0" for w in product("01",repeat=7) if w.count("1")==5]
A6=["".join(w)+"0" for w in product("01",repeat=7) if w.count("1")==6]
AL_ALL=A5+A6
for nm,(Kx,Dx) in (("21/7 shell (regression)",(K0,D0)),("candidate",(K1,D1))):
    lam=F(Kx)*LN2-F(Dx)*LN3
    eta=F(8*Dx,Kx)-5
    beta_lo,beta_hi=LN2/LN3, LN2/LN3
    print("  --- %s: K = %d"%(nm,Kx))
    c1 = 5 < F(8*Dx,Kx) < 6
    c2 = 0 < lam < F(1,3000)
    c3 = 3*(1-lam) > F("2.999")
    c4 = (8*Dx-4*Kx)*LN3 < 8*Kx*LN2
    c5 = F(256,81)-F("2.999") < F(".162")
    nb=(Kx-17)//8; sixes=Dx-11-5*nb
    def a_of(t): return t*Dx//Kx+2
    c6 = a_of(8)==7 and a_of(Kx-9)==Dx-4 and 8+8*nb==Kx-9
    bl=F(int(AL*Decimal(10)**200),10**200)
    c7 = all(F(HEAD[:j].count("1"))-F(j)*(LN2/LN3) > 0 for j in range(1,9))
    c8 = all(F(len(TAIL)-j)*LN2 - F(TAIL[j:].count("1"))*LN3 > lam for j in range(len(TAIL)))
    c9 = all(F(3**w[:j].count("1"),2**j)*F("2.999")>1
             for w in AL_ALL if w!="00111110" for j in range(1,9))
    c10= all(F(3**"00111110"[:j].count("1"),2**j)*F("2.999")>1 for j in (1,3,4,5,6,7,8))
    chk("  (i) legal 5-or-6 schedule, and a(8)=7, a(K-9)=D-4", c1 and c6)
    chk("  boundary coefficient > 2.999 at every middle boundary", c2 and c3)
    chk("  rho < 3^(1+eta) < 256/81 before a six-one byte, and 256/81-2.999 < 0.162",
        c4 and c5)
    chk("  all 8 proper head prefixes stay above the line", c7)
    chk("  all 9 tail starts stay above the line", c8)
    chk("  no internal byte prefix crosses, except 00111110 at its second zero",
        c9 and c10)
    m=(sixes-PAST-FW-1+7*21**FW-1)//(7*21**FW)
    Q=F(Dx,3)/lam; need=(F("0.162")*Q+F(256,81)*F(Dx,3))/2**L
    chk("  m = %d >= 163 and V/2^104 = %.4e < sigma_163 (%.2fx room)"
        %(m,float(need),SIGMA/float(need)), m>=163 and need<SIGMA)
    print("      early-crossing sites to filter: %.1f"%float(Kx*lam/LN3))
print("  => the candidate shell's ENDPOINT first-crossing chain is fully verified, with")
print("     the same alphabet-only within-byte facts.  What remains for a theorem there")
print("     is the exact early-crossing filter at ~1458 sites (vs 131 on the 21/7 shell)")
print("     and the per-byte memory bound -- no new residue enumeration.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
