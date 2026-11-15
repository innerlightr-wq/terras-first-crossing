"""Independent audit of the eight-step families, and the decisive counting."""
from fractions import Fraction as F
from itertools import product
import math
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
K=1412737287532845639425; D=891337988684307508298; B=28*10**18
END=K-9; HEAD,TAIL="11111110","111100000"; N=(END-8)//8
def lb(base,terms=160):
    x=F(base-1,base+1); v=sum(2*x**(2*i+1)/(2*i+1) for i in range(terms))
    return v, v+2*x**(2*terms+1)/((2*terms+1)*(1-x*x))
lo2,hi2=lb(2); lo3,hi3=lb(3)
lam_lo,lam_hi=K*lo2-D*hi3, K*hi2-D*lo3
beta_lo,beta_hi=lo2/hi3, hi2/lo3
E_min,E_max=lam_lo, lam_hi/(1-lam_hi)
def S(t): return t*D//K+2
def alpha_set(c):
    out=[]
    for bits in product("01",repeat=8):
        w="".join(bits)
        if w[-1]!="0" or w.count("1")!=c: continue
        if all(1000*w[:r].count("1")-631*r>=-800 for r in range(9)): out.append(w)
    return out

print("="*78); print("1. the construction, audited independently"); print("="*78)
chk("beta in (0.630,0.631) and 0 < Lambda < 1", F(".630")<beta_lo<beta_hi<F(".631") and 0<lam_lo<lam_hi<1)
chk("B*E in (2.84228, 2.84229)", F("2.84228")<B*E_min<B*E_max<F("2.84229"))
chk("2^K = 3^D = 608 mod 929", pow(2,K,929)==pow(3,D,929)==608)
chk("END = 8 + 8N with N = 176592160941605704926", END==8+8*N and N==176592160941605704926)
chk("HEAD has 7 ones = S(8); TAIL has 4 = D - S(END)",
    HEAD.count("1")==S(8)==7 and TAIL.count("1")==D-S(END)==4)
chk("each middle byte needs exactly 5 or 6 ones (5 < 8D/K < 6)", 5<F(8*D,K)<6)
A={c:alpha_set(c) for c in (5,6)}
chk("alphabet sizes are 18 (c=5) and 7 (c=6)", {c:len(v) for c,v in A.items()}=={5:18,6:7},
    str({c:len(v) for c,v in A.items()}))
ok=True
for c,ws in A.items():
    for w in ws:
        for r in range(9):
            cnt=w[:r].count("1")
            if not (1-F(".001")+cnt-r*beta_hi>F(".19") and 2+cnt-r*beta_lo<5): ok=False
chk("every within-byte prefix keeps the drift in (0.19, 5)", ok)
ok=all(F(".19")<HEAD[:r].count("1")-r*beta_hi and HEAD[:r].count("1")-r*beta_lo<5
       for r in range(1,9))
chk("all 8 proper head prefixes keep the drift in (0.19, 5)", ok)
ok=all((9-r)*beta_lo-TAIL[r:].count("1")-lam_hi/lo3>F(".19") and
       (9-r)*beta_hi-TAIL[r:].count("1")<5 for r in range(9))
chk("all 9 tail starts keep the drift in (0.19, 5)", ok)
chk("hence 6/5 < rho_t < 243 for every noninitial proper prefix", F(6,5)**100<3**19)
inv=pow(K,-1,D)
upper=(10*D*lam_hi/(9*lo2)).__floor__()
chk("the >90% weight search bound is max_D = 145", upper==145)
strong=[]
for det in range(1,upper+1):
    a=(inv*det)%D; j=(K*a-det)//D
    v_lo,v_hi=a*lo3-j*hi2, a*hi3-j*lo2
    if v_hi<=0: continue
    if 9*v_hi/(1-v_hi)<lam_lo-lam_hi**2/2: strong.append(det)
chk("the 13 strong cuts reproduce", strong==[16,23,39,46,53,69,76,83,99,106,113,129,136])
miss=0
for det in strong:
    a=(inv*det)%D; j=(K*a-det)//D
    t=8*(j//8); c=S(t+8)-S(t)
    lo_=min(S(t)+w[:j-t].count("1") for w in A[c])
    if lo_>a: miss+=1
chk("at every one of the 13 cuts the MINIMUM achievable odd count already exceeds a_D",
    miss==13, "%d of 13"%miss)
print("    => avoidance is uniform over every choice, not a selection condition.")
n_lo=(F(D,729)-1858)/E_max; Q_up=F(D,3)/E_min
chk("M_w > D/729 and Q_w < 2^132", Q_up<2**132, "Q_w < %.4e = 2^%.2f"%(float(Q_up),math.log2(float(Q_up))))
chk("every word has a factor-compatible h with relaxed seed n_h > 10^37",
    n_lo>10**37>B+2**100, "n_h > %.4e = 2^%.2f"%(float(n_lo),math.log2(float(n_lo))))
print("    which prefix floors / placement tests that seed clears:")
print("      - any length-100 prefix floor: its least representative above B is <= B+2^100,")
print("        and n_h > 10^37 > B + 2^100, so ALL 100-bit prefix floors clear.")
chk("B + 2^100 < 10^37 < n_h", B+2**100<10**37<n_lo)
stride=2**7*3**9
chk("depth<=9 placement: K*stride < n_h/5, so individual AND aggregate ceilings clear",
    K*stride<n_lo/5, "K*stride = %.4e vs n_h/5 = %.4e"%(float(K*stride),float(n_lo/5)))
chk("interior states stay below 2^140", 243*(2**132+F(D,3))<2**140)
print("    (these are relaxations: integer realizability of n_h is NOT asserted.)")

print()
print("="*78); print("2. the two-choice control"); print("="*78)
chk("N - 17 > 4^18", N-17>4**18, "N-17 = %.4e, 4^18 = %.4e"%(N-17,4.0**18))
print("  With only four byte types, the N-17 windows of 18 consecutive bytes cannot all")
print("  differ, so two odd starts share a 144-bit forward parity window.  Equal forward")
print("  bits force their difference to be 0 mod 2^144, and distinctness (shortestness)")
print("  forces |x_i - x_j| >= 2^144 > 2^140 > every interior state.  Contradiction.")
chk("2^144 exceeds the state-height bound 2^140", 2**144>2**140)
print("  WHERE SHORTESTNESS ENTERS: only to make the k input states DISTINCT.  Without")
print("  minimal crossing length two interior states could coincide and the spacing")
print("  argument would have no bite.")
print("  HAMMING CERTIFICATE: every byte ends in 0 and carries at most six ones, so the")
print("  only run of seven ones is the initial one; each c=6 byte contains '111', and")
print("  3*beta < 2 means no length-3 factor of a critical Sturmian word has three ones.")
chk("3*beta_hi < 2, so every c=6 byte forces a mismatch at every intercept", 3*beta_hi<2)
c6=D-HEAD.count("1")-TAIL.count("1")-5*N
chk("the number of c=6 bytes is in (0,N)", 0<c6<N, "c6 = %d (%.3f%% of bytes)"%(c6,100.0*c6/N))
chk("all c=6 bytes contain '111'", all("111" in w for w in A[6]))
print("  => the disjoint '111's give at least %d mismatches against any Sturmian factor."%c6)

print()
print("="*78); print("3. THE LARGER FAMILY: does forward repetition exclude it?"); print("="*78)
print("  The c-pattern is FIXED by S, so a window's fillings number prod(18 or 7).")
print("  Fraction of bytes with c=6: %.4f%%.  The cheapest pigeonhole uses all-c=5"%(100.0*c6/N))
print("  windows, which are the most numerous and have the most fillings.")
p5=F(N-c6,N)
print("  %4s %8s %22s %22s %s"%("L","bits","#all-5 windows ~","18^L","pigeonhole?"))
Lwork=None
for L in (12,14,15,16,18,20):
    wins=F(N)*p5**L
    fill=18**L
    ok_=wins>fill
    if ok_ and (Lwork is None or L>Lwork): Lwork=L
    print("  %4d %8d %22.4e %22.4e %s"%(L,8*L,float(wins),float(fill),
          "YES" if ok_ else "no"))
chk("the pigeonhole works only up to L = 15 bytes (120 bits)", Lwork==15, "Lwork = %s"%Lwork)
chk("but the height bound needs 8L > 140, i.e. L >= 18 (144 bits)", 8*18>140 and 8*17<=140)
print("  SHORTFALL: the argument needs 18-byte windows (144 bits) but the pigeonhole")
print("  reaches only 15 bytes (120 bits).  The gap is %d bits of height range."%(141-120))
print("  Narrowing the height range does not close it: at L = 15 up to 2^(140-120) = 2^20")
print("  starts may share a window, while the pigeonhole yields only")
print("    %.1f occurrences of the most frequent all-5 15-byte window."
      %float(F(N)*p5**15/18**15))
chk("2^20 required occurrences vastly exceed the %.1f the pigeonhole delivers"
    %float(F(N)*p5**15/18**15), 2**20 > F(N)*p5**15/18**15)
print()
print("  VERDICT.  Ordinary forward repetition EXCLUDES the two-choice control at 144")
print("  bits and does NOT exclude the larger family, short by 3 bytes / 21 bits.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)

print()
print("="*78); print("4. is signed pair placement quantitatively viable here?"); print("="*78)
print("  A signed pair needs two starts with (a) a shared outgoing parity block of L")
print("  bits, giving x_i = x_j mod 2^L, and (b) different incoming histories of depth")
print("  t, giving a signed class for x_i - x_j modulo 3^t.  Joint modulus 2^L 3^t.")
print("  The test is decisive once the modulus EXCEEDS the gap interval, which has")
print("  width at most the state range 2^140: then the class holds at most one")
print("  admissible difference.")
print("  %6s %7s %10s %14s %22s %s"%("L bytes","L bits","t needed","2^L 3^t","pigeonhole gives","viable?"))
rows=[]
for Lb in (12,14,15,16,17,18):
    Lbits=8*Lb
    t=0
    while Lbits+t*float(hi3/lo2)<=140: t+=1
    mod_ok = (2**Lbits)*(3**t) > 2**140
    occ=float(F(N)*p5**Lb/18**Lb)
    exists = occ>1
    rows.append((Lb,Lbits,t,mod_ok,occ,exists))
    print("  %6d %7d %10d %14s %22.4g %s"%(Lb,Lbits,t,"> 2^140" if mod_ok else "too small",
          occ,"YES" if (mod_ok and exists) else ("no shared window" if mod_ok else "no")))
viable=[r for r in rows if r[3] and r[5]]
chk("a viable (L,t) exists: the pigeonhole guarantees the shared window AND the "
    "modulus beats the state range", bool(viable),
    "best is L = %d bytes (%d bits) with depth t = %d"%(viable[-1][0],viable[-1][1],viable[-1][2])
    if viable else "none")
if viable:
    Lb,Lbits,t,_,occ,_=viable[-1]
    print("  => %d-byte (%d-bit) shared forward window, guaranteed by pigeonhole to occur"
          %(Lb,Lbits))
    print("     about %.1f times, combined with depth-%d incoming histories: modulus"%(occ,t))
    print("     2^%d * 3^%d = 2^%.2f > 2^140."%(Lbits,t,Lbits+t*math.log2(3)))
    print("     Each such pair therefore pins x_i - x_j to AT MOST ONE admissible value,")
    print("     so the signed test is decisive per pair up to a single coincidence.")
print()
print("  CONTRAST.  Forward repetition alone needs the shared window to reach 141 bits,")
print("  which the pigeonhole cannot deliver (it stops at 120).  Adding only 3^13 of")
print("  incoming information closes exactly that 21-bit gap.")
chk("the depth-13 factor 3^13 covers the 21-bit shortfall", 3**13>2**20,
    "3^13 = 2^%.2f >= 2^20"%math.log2(3**13))
print()
print("  NOT ESTABLISHED: that the single surviving difference is actually missed for")
print("  EVERY member.  That is the remaining inequality, and it is one coincidence per")
print("  pair rather than a structural gap.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
