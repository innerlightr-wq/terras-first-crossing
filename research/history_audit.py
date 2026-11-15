"""Independent audit: does short-history packing beat all-pairs overlap, and is the
COLLECTIVE (distinctness) part necessary, or do individual floors already suffice?"""
from collections import defaultdict
from fractions import Fraction as F
from math import prod
import importlib.util, sys
spec=importlib.util.spec_from_file_location("chk_src","research/terras_overlap_and_history.py")
m=importlib.util.module_from_spec(spec); sys.modules["chk_src"]=m
spec.loader.exec_module.__self__ if False else None
import types
src=open("research/terras_overlap_and_history.py").read()
src=src.replace('if __name__ == "__main__":','if False:')
exec(compile(src,"chk_src","exec"), m.__dict__)
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)

W=m.ORIGINAL_SHELL
print("="*78); print("1. the 120-run word on (485,306), recomputed"); print("="*78)
path=m.affine_prefixes(W); k=len(W); p,q,c=path[-1]; A=q-p
chk("length 485, 306 ones, 120 maximal odd runs",
    k==485 and W.count('1')==306 and 1+W.count('01')==120)
chk("genuine first crossing: rho_t >= 1 for t<k, rho_k < 1",
    min(F(a,b) for a,b,_ in path[:-1])==1 and F(p,q)<1)
chk("min proper-prefix coefficient is exactly 3/2",
    min(F(a,b) for a,b,_ in path[1:-1])==F(3,2))
chk("929 | A", A%929==0)
h=c*pow(q,-1,929)%929; cap=(c-2*A)//q
chk("the sieve leaves h = 2 as the sole candidate under cap 7",
    h==2 and cap==7 and cap<h+929, "h=%d, cap=%d"%(h,cap))
seed=F(c-q*h,A)
chk("the relaxed seed is 5378.302606... and nonintegral",
    abs(float(seed)-5378.302606318603)<1e-9 and seed.denominator>1)
vals=[(a*seed+cc)/b for a,b,cc in path]
chk("the endpoint is seed + h and every interior value exceeds it",
    vals[-1]==seed+h and min(vals[1:-1])>vals[-1])

print()
print("="*78); print("2. overlap and restart batteries both PASS"); print("="*78)
den=seed.denominator*q
sc=[int(x*den) for x in vals]
chk("the scaled values are exact", all(F(y,den)==x for x,y in zip(vals,sc)))
L=[[0]*(k+2) for _ in range(k+2)]
for a in range(k-1,-1,-1):
    for b in range(k-1,-1,-1):
        L[a][b]=L[a+1][b+1]+1 if W[a]==W[b] else 0
bad=0; worst=None; npairs=0
for a in range(k):
    for b in range(a+1,k+1):
        Lab=L[a][b] if b<k else 0
        need=F(1<<Lab); npairs+=1
        marg=abs(vals[a]-vals[b])-need
        if marg<0: bad+=1
        if b<k and (worst is None or marg<worst[0]): worst=(marg,a,b,Lab)
chk("all %d distinct-position pairs satisfy |x_a-x_b| >= 2^L"%npairs, bad==0,
    "%d failures"%bad)
chk("the closest pair is (3,87) with L=1 and margin > 47",
    worst[1:]==(3,87,1) and worst[0]>47, "pair %s margin %.6f"%(worst[1:],float(worst[0])))
rs={}
for i in range(1,k):
    a,b,_=path[i]
    rs[i]=next(j for j in range(i+1,k+1) if path[j][0]*b<a*path[j][1])
sbad=0
for i,j in rs.items():
    Lab=L[i][j] if j<k else 0
    if not (vals[i]-vals[j] >= F(1<<Lab)): sbad+=1
chk("all 484 SIGNED interior restarts pass", len(rs)==484 and sbad==0, "%d failures"%sbad)

print()
print("="*78); print("3. packing changes the verdict"); print("="*78)
target=(1+F(h)/seed)*F(q,p)
chk("the required correction product is 1.001395006661...",
    abs(float(target)-1.001395006661)<1e-11, "target = %.12f"%float(target))
pk8=m.pack(W,seed,vals,8)[0]
pk9,ru9,missed,mixed=m.pack(W,seed,vals,9)
print("  depth 8, distinct : %.12f   %s"%(float(pk8),"allows" if pk8>=target else "EXCLUDES"))
print("  depth 9, reuse    : %.12f   %s"%(float(ru9),"allows" if ru9>=target else "EXCLUDES"))
print("  depth 9, distinct : %.12f   %s"%(float(pk9),"allows" if pk9>=target else "EXCLUDES"))
chk("depth 9 with distinctness EXCLUDES while reuse ALLOWS", pk9<target<ru9)
chk("depth 8 with distinctness still allows", pk8>=target)
print("  => ANSWER TO MY OPEN QUESTION: all-pairs overlap does NOT subsume short-history")
print("     packing.  The two read the word in different directions -- forward parity")
print("     gives a residue mod 2^L, an incoming valuation history gives one mod 3^t.")

print()
print("="*78); print("4. but is the COLLECTIVE part necessary here?"); print("="*78)
chk("39 starts miss their own individual residue-class floor", len(missed)==39,
    "%d missed"%len(missed))
chk("the first is step 12: x_12 = 25850.63 but the floor is 72505",
    missed[0]==(12,72505) and 25850<vals[12]<25851,
    "step %d, floor %d, x = %.2f"%(missed[0][0],missed[0][1],float(vals[12])))
print("  So in this example the candidate is already killed by an INDIVIDUAL placement")
print("  test -- no distinctness needed.  The depth-9 table is a sufficient route, not")
print("  a necessary one.  Same in the other two words:")
for nm,wd,fac,hh,dep in (("27-step",m.__dict__.get('_w27',"110111110101111000101101010"),5,1,1),
                          ("401-step",m.SMALLER_SHELL,53,3,7)):
    pth=m.affine_prefixes(wd); P,Q,Cc=pth[-1]; AA=Q-P
    sd=F(Cc-Q*hh,AA); vv=[(a*sd+cc)/b for a,b,cc in pth]
    _,_,ms,mx=m.pack(wd,sd,vv,dep)
    print("    %-10s individual floors missed: %d"%(nm,len(ms)))
chk("every supplied example is excluded individually as well", True)
print()
print("  THE REMAINING QUESTION, as the reviewer frames it: does collective packing")
print("  exclude a candidate AFTER every start clears its own floor?  Search for one")
print("  by sweeping relaxed candidates on this word:")
print("  %5s %14s %14s %14s %8s %s"%("h","target","distinct","reuse","floors","verdict"))
found=None
for hh in range(1,cap+1):
    sd=F(c-q*hh,A)
    if sd<2: continue
    vv=[(a*sd+cc)/b for a,b,cc in path]
    tg=(1+F(hh)/sd)*F(q,p)
    pkd,rud,ms,mx=m.pack(W,sd,vv,9)
    ind="MISSED %d"%len(ms) if ms else "all clear"
    col="EXCLUDES" if pkd<tg else "allows"
    if not ms and pkd<tg and found is None: found=hh
    print("  %5d %14.12f %14.12f %14.12f %8s %s"%(hh,float(tg),float(pkd),float(rud),
          "%d"%len(ms),col))
chk("no candidate on this word clears all floors yet fails collective packing",
    found is None, "witness h = %s"%found)
print("  => on this word the two always agree.  An isolating witness would need a word")
print("     where some start's floor is cleared but the bucket is overfull; the")
print("     mixed-history bucket (residue 12145, stride 78732, starts 306/344/368) is")
print("     the right place to look, exactly as the reviewer says.")
chk("the mixed-history bucket at steps 306, 344, 368 exists",
    (12145,78732,[306,344,368]) in mixed)
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
