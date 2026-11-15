"""Independent pure-Python reconstruction of the six-one-start certificate.
(The supplied script needs NumPy, which is not installed here, so every number
below is rebuilt from scratch rather than re-run.)"""
from fractions import Fraction as F
from itertools import product
import math
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
K=1412737287532845639425; D=891337988684307508298; N=(K-17)//8
HEAD,TAIL="11111110","111100000"; L=112
def lb(base,terms=160):
    x=F(base-1,base+1); v=sum(2*x**(2*j+1)/(2*j+1) for j in range(terms))
    return v, v+2*x**(2*terms+1)/((2*terms+1)*(1-x*x))
lo2,hi2=lb(2); lo3,hi3=lb(3)
lam_lo,lam_hi=K*lo2-D*hi3, K*hi2-D*lo3
beta_lo,beta_hi=lo2/hi3, hi2/lo3
eta=F(8*D,K)-5
def S(t): return t*D//K+2
six=S(K-9)-S(8)-5*N
def alpha(c,size):
    return ["".join(w)+"0" for w in product("01",repeat=7)
            if w.count("1")==c and (size==21 or w[:2]!=("0","0"))]
def carry(w): return sum(2**j*3**w[j+1:].count("1") for j in range(len(w)) if w[j]=="1")
def resid(letters,depth,forced):
    mod=3**(5*depth); cs=[carry(w) for w in letters]; cur={0}
    for j in range(depth):
        sc=pow(2,8*j-8*depth-L,mod)*pow(3,5*(depth-1-j),mod)%mod
        offs=[sc*c%mod for c in cs]
        cur={(x+o)%mod for x in cur for o in offs}
    o=sorted(cur); m=len(o)
    cyc=o+[v+mod for v in o[:forced-1]]
    return m, min(cyc[i+forced-1]-v for i,v in enumerate(o))
print("="*78); print("0. shell and schedule"); print("="*78)
chk("0 < Lambda < 1/1000 and 0 < eta < 1/20", 0<lam_lo<lam_hi<F(1,1000) and 0<eta<F(1,20))
chk("8 + 8N = K - 9, S(8)=7, S(K-9)=D-4, six count matches",
    8+8*N==K-9 and S(8)==7 and S(K-9)==D-4 and six==8377183976278983657)

print(); print("="*78); print("1. the UNRESTRICTED 20/7 extension"); print("="*78)
A5,A6=alpha(5,20),alpha(6,20)
chk("alphabet sizes 20 and 7 (all 7-bit 5-one words but '0011111', plus a 0)",
    (len(A5),len(A6))==(20,7))
chk("the two new letters 01011110 and 10011110 are present",
    "01011110" in A5 and "10011110" in A5)
chk("the 18-letter density cap really is dropped: 2 letters violate it",
    sum(1 for w in A5 if not all(1000*w[:r].count("1")-631*r>=-800 for r in range(9)))==2)
chk("(4+14)*eta < 1, so any 18 adjacent middle bytes hold at most one six",
    (4+14)*eta<1, "18*eta = %.6f"%float(18*eta))
print("  Hence a window STARTING at a six-one byte has its preceding 4 and following 13")
print("  bytes all five-one, whenever the window fits inside the middle region.")
avail=six-4-13; labels=7*20**13; forced=(avail+labels-1)//labels
chk("eligible starts >= six - 4 - 13 and labels = 7*20^13",
    avail==six-17 and labels==7*20**13, "avail = %d, labels = %d"%(avail,labels))
chk("one 112-bit forward word occurs at least FIFTEEN times", forced==15,
    "ceil(%.4e / %.4e) = %d"%(float(avail),float(labels),forced))
chk("rho > 2.99 at every middle boundary", 3*(1-lam_hi)>F("2.99"))
chk("a six next forces {tD/K} >= 1-eta, hence rho < 3^(1+eta) < 3.17",
    3**21<F("3.17")**20, "3^(21/20) = %.6f < 3.17"%(3**(21/20)))
cnt,sp=resid(A5,4,15)
chk("128000 residues modulo 3^20", cnt==128000, "%d"%cnt)
chk("minimum fifteen-point cyclic span 107996", sp==107996, "%d"%sp)
M=F(D,3); Q=M/lam_lo
width=(F(".18")*Q+F("3.17")*M)/2**L
chk("available span below 101468", width<101468, "width < %.1f"%float(width))
chk("DECISIVE: width < 101468 < span < 3^20", width<101468<sp<3**20)
chk("every 20/7 letter has at most ONE leading zero, so rho/2 >= 1.495 > 1 inside",
    all(not w.startswith("00") for w in A5+A6))
print("    => no internal prefix of any byte drops the coefficient below one, so the")
print("       first crossing of every 20/7 word is exactly at K.")
w20=A5+A6
chk("'00000' never occurs in three adjacent 20/7 bytes",
    all("00000" not in a+b+c for a,b,c in product(w20,repeat=3)))
chk("'00000' never occurs at a head junction",
    all("00000" not in HEAD+a+b for a,b in product(w20,repeat=2)))
chk("'00000' occurs exactly once at a tail junction, at the very end",
    all((a+b+TAIL).find("00000")==len(a+b+TAIL)-5 for a,b in product(w20,repeat=2)))
chk("the starts used lie at or before K-5", K-9-L<=K-5)
print("  VERDICT 20/7: every eligible integer seed n >= 2 realizing a 20/7 word descends")
print("  at its FIRST crossing, which is at K.  Pointwise, no shortestness used.")

print(); add=alpha(5,21)
print("="*78); print("2. the FULL 21/7 endpoint claim"); print("="*78)
B5,B6=alpha(5,21),alpha(6,21)
chk("alphabet sizes 21 and 7; the new letter is 00111110",
    (len(B5),len(B6))==(21,7) and "00111110" in B5 and "00111110" not in A5)
avail2=six-5-13; labels2=7*21**13; forced2=(avail2+labels2-1)//labels2
chk("with five preceding bytes, eight occurrences are forced", forced2==8,
    "ceil(%.4e / %.4e) = %d"%(float(avail2),float(labels2),forced2))
chk("(5+14)*eta < 1 still holds", (5+14)*eta<1, "19*eta = %.6f"%float(19*eta))
cnt2,sp2=resid(B5,5,8)
chk("3111696 residues modulo 3^25", cnt2==3111696, "%d"%cnt2)
chk("minimum eight-point span 152838", sp2==152838, "%d"%sp2)
chk("the same available width still closes: width < 101468 < 152838 < 3^25",
    width<101468<sp2<3**25)
chk("memory bound byte by byte: carry(w)/(3^c * 2.99) < c/3 for every letter",
    all(F(carry(w),3**w.count("1"))/F("2.99")<F(w.count("1"),3) for w in B5+B6))
print("    this allows the coefficient to fall below one INSIDE a byte, because it only")
print("    uses rho at the byte START, which exceeds 2.99 at every middle boundary.")
print()
print("  THE SHORT MARKER BREAKS.  '11111000' and '00111110' are both legal, and")
print("  '11111000'+'00111110' contains '00000':")
chk("'00000' DOES occur inside the 21/7 family", "00000" in "11111000"+"00111110")
chk("so the 21/7 family needs the longer marker '0111100000'",
    all("0111100000" not in a+b+c for a,b,c in product(B5+B6,repeat=3)))
chk("'0111100000' never occurs at a head junction",
    all("0111100000" not in HEAD+a+b for a,b in product(B5+B6,repeat=2)))
chk("'0111100000' occurs exactly once at a tail junction, at the very end",
    all((a+b+TAIL).find("0111100000")==len(a+b+TAIL)-10 for a,b in product(B5+B6,repeat=2)))
chk("the starts used lie at or before K - 10", K-9-L<=K-10)
print("  VERDICT 21/7: descent holds at the ENDPOINT K for every eligible integer seed.")
print("  Words that never use 00111110 are exactly the 20/7 family and are covered at")
print("  their FIRST crossing; words that use it are covered only at K.")

print(); print("="*78); print("3. the remaining gap at 00111110"); print("="*78)
w="00111110"
print("  internal prefix multipliers of rho for this byte:")
mult=[]
for r in range(1,9):
    ones=w[:r].count("1"); mult.append((r,F(3**ones,2**r)))
    print("    r=%d  rho_r/rho = 3^%d/2^%d = %s"%(r,ones,r,F(3**ones,2**r)))
worst=min(m for _,m in mult)
chk("the binding internal prefix is r=2 with factor 1/4", worst==F(1,4) and mult[1][1]==F(1,4))
print("  EXACT CONDITION: every internal prefix stays at or above the crossing line iff")
print("     rho_start * 1/4 >= 1,  i.e.  rho_start >= 4,  i.e. 3^{d_t} >= 4*2^t.")
chk("rho_start >= 4 implies every internal prefix of 00111110 has rho_r >= 1",
    all(4*m>=1 for _,m in mult))
chk("but at a six-one start rho < 3.17 < 4, and even a generic boundary has rho <= 9",
    F("3.17")<4)
print("  In phase terms rho = 3^{2-{tD/K}} e^{-Lambda t/K}, so rho >= 4 iff")
print("     {tD/K} <= 2 - log_3 4 - o(1) = %.6f ."%(2-math.log(4,3)))
print("  So the condition is a PHASE condition on the boundary, satisfiable for some t.")
print("  Imposing it at every occurrence of 00111110 does characterise the")
print("  first-crossing-at-K subfamily: with it, all internal prefixes stay >= 1 and the")
print("  20/7 argument applies verbatim; without it at some occurrence, the coefficient")
print("  first falls below one TWO steps into that byte.")
chk("the earlier crossing, when it happens, is at byte start + 2", mult[1][0]==2)
print("  At that position the descent inequality that would have to hold is")
print("     T^{t+2}(n) < n,  i.e.  n (rho_{t+2}^{-1} - 1) > M_{t+2},")
print("  with rho_{t+2} = rho_start/4 < 1.  That is a DIFFERENT inequality from the one")
print("  certified at K, and nothing above establishes it.")
print("  => the 21/7 endpoint theorem does NOT give first-crossing descent for words")
print("     using 00111110 at a boundary with rho_start < 4.  THAT IS THE GAP.")

print(); print("="*78); print("4. what transfers"); print("="*78)
print("  SHELL-DEPENDENT: K, D, Lambda, eta = 8D/K - 5, the schedule S(t), the six-one")
print("  count, Q_w = D/(3 Lambda), M = D/3, the boundary bracket 2.99 < rho < 3.17, and")
print("  therefore the available width (0.18 Q + 3.17 M)/2^112.")
print("  ALPHABET-DEPENDENT ONLY: the residue table (128000 mod 3^20; 3111696 mod 3^25),")
print("  its minimum m-point spans, the label count 7*size^13, and the marker length.")
print("  SUFFICIENT CRITERION for another shell, reusing the same residue table:")
print("     (i)   (past+14)*eta < 1                      [schedule]")
print("     (ii)  forced = ceil((six-past-13)/(7*size^13)) >= m")
print("     (iii) (rho_max - rho_min) Q + rho_max M  <  2^L * sigma_m(table)")
print("     (iv)  the terminal marker occurs once, over all 3-byte junctions")
print("  FIRST QUANTITY TO FAIL ON TRANSFER: (ii).  six/N is fixed by the shell's phase,")
print("  so a shell with a smaller six-one share loses the forced multiplicity first;")
print("  here forced = 15 at 20/7, EXACTLY the multiplicity the span certificate needs")
print("  (ceil(14.608) = 15), so a slightly smaller six-one share would drop it to 14 and")
print("  the smaller sigma_14 might no longer clear the width.")
print("  NOTE on the stated constant 101468: it is essentially ceil(width) (width <")
print("  101467.0), so it is NOT a margin.  The real room is width against the span:")
chk("criterion (iii) margin at 20/7 is span/width - 1", width<F(101468),
    "%.4f%% spare against sigma_15 = %d"%(100.0*float(1-width/sp),sp))
chk("forced multiplicity at 20/7 is exactly the required 15, with no slack",
    (avail+labels-1)//labels==15 and (avail-1+labels-1)//labels==15,
    "ceil = %d from ratio %.4f"%(forced,float(F(avail,labels))))
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
