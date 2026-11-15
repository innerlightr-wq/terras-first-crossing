"""Second-shell first-crossing theorem: exclude every early-crossing site."""
from fractions import Fraction as F
from decimal import Decimal, getcontext
from itertools import product
from math import gcd
import math
getcontext().prec=260
AL=Decimal(3).ln()/Decimal(2).ln()
ALF=F(int(AL*Decimal(10)**220),10**220)
def lg(base,terms=240):
    x=F(base-1,base+1); lo=sum(2*x**i/i for i in range(1,terms,2))
    return lo, lo+2*x**(terms+1)/((terms+1)*(1-x*x))
lo2,hi2=lg(2); lo3,hi3=lg(3)
SIGMA=3611678508; L=104; PAST=7; FW=12; S=21; labels=7*S**FW
HEAD,TAIL="11111110","111100000"
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
A5=["".join(w)+"0" for w in product("01",repeat=7) if w.count("1")==5]
A6=["".join(w)+"0" for w in product("01",repeat=7) if w.count("1")==6]
def mcost(b): return sum(F(2**i,3**(b[:i].count("1")+1)) for i,bit in enumerate(b) if bit=="1")
BAD="00111110"

def analyse(K,D,label,verbose=True):
    lam_lo,lam_hi=K*lo2-D*hi3, K*hi2-D*lo3
    if lam_lo<=0: return None
    eta=F(8*D,K)-5
    nb=(K-17)//8; end=K-9
    sixes=D-11-5*nb
    def a_of(t): return t*D//K+2
    ok_sched = (K%8==1 and 8+8*nb==end and a_of(8)==7 and a_of(end)==D-4
                and 5<F(8*D,K)<6)
    m=(sixes-PAST-FW-1+labels-1)//labels
    Q=F(D,3)/lam_lo
    need=(F("0.162")*Q+F(256,81)*F(D,3))/2**L
    floor=(2**L*SIGMA-F(256,81)*F(D,3))/F(".162")
    res=dict(K=K,D=D,sched=ok_sched,m=m,need=need,floor=floor,
             sep=20*eta<1, lam_lo=lam_lo, lam_hi=lam_hi,
             hmax=float(K*lam_hi/lo3), endpoint=floor*lam_lo>F(D,3))
    if not (ok_sched and m>=163 and need<SIGMA and res['sep']): return res
    # enumerate early-crossing sites
    inv=pow(D,-1,K)
    hbound=int(K*lam_hi/lo3)+1
    sites=[]; coarse=0
    for h in range(hbound+1):
        j=(-h*inv)%K; t=j-2
        if not (8<=t<end and t%8==0): continue
        a=a_of(t)
        if K*a-D*j!=h or a_of(t+8)!=a+5: continue
        elo,ehi=j*lo2-a*hi3, j*hi2-a*lo3
        if ehi<0: continue
        if elo<=0: continue
        if floor*elo>F(a,3): coarse+=1; continue
        cap=F(7,3)+F(".94")*((t-8)//8)
        sites.append((h,j,t,a,elo,cap,cap/elo))
    res['hbound']=hbound; res['coarse']=coarse; res['sites']=sites
    res['survivors']=[s for s in sites if s[6]>=floor]
    return res

print("="*78); print("alphabet-only facts (shell independent), reused"); print("="*78)
chk("an earlier completed 00111110 starts at rho>=4, cost/4 < 0.94", mcost(BAD)/4<F(".94"))
chk("regular and six-one bytes cost/2.999 < 0.94",
    max(map(mcost,[w for w in A5 if w!=BAD]))/F("2.999")<F(".94")
    and max(map(mcost,A6))/F("2.999")<F(".94"))
chk("the head costs < 7/3", mcost(HEAD)<F(7,3))
chk("only 00111110 crosses early, at its second zero",
    all(F(3**w[:j].count("1"),2**j)*F("2.999")>1
        for w in A5+A6 if w!=BAD for j in range(1,9))
    and F(3**BAD[:2].count("1"),2**2)==F(1,4))

print(); print("="*78); print("Q1.  the 1388 sites on the corrected strongest shell"); print("="*78)
K1,D1=1431227343753494174081,903003915299123681908
r=analyse(K1,D1,"candidate")
chk("schedule, separation, multiplicity and span all hold",
    r['sched'] and r['sep'] and r['m']>=163 and r['need']<SIGMA,
    "m = %d, span room %.2fx"%(r['m'],SIGMA/float(r['need'])))
chk("Lambda ~ 1.06557e-18 as independently spot-checked",
    F("1.0655")*10**-18<r['lam_lo'] and r['lam_hi']<F("1.0656")*10**-18,
    "Lambda in (%.6e, %.6e)"%(float(r['lam_lo']),float(r['lam_hi'])))
chk("the seed floor settles the ENDPOINT case", r['endpoint'],
    "floor*Lambda = %.4e > D/3 = %.4e"%(float(r['floor']*r['lam_lo']),float(F(D1,3))))
print("  h-range: -1 < h < %.1f  ->  %d integer candidates"%(r['hmax'],r['hbound']+1))
print("  of those, %d are genuine byte-aligned sites with E>0"%len(r['sites']))
print("  (%d were removed by the coarse test floor*E > a/3)"%r['coarse'])
chk("EVERY early-crossing site is excluded by the refined bound cap/E < floor",
    len(r['survivors'])==0, "%d survivors of %d sites"%(len(r['survivors']),len(r['sites'])))
if r['sites']:
    worst=max(r['sites'],key=lambda s:s[6])
    print("  tightest site: h = %d, j = %d"%(worst[0],worst[1]))
    print("     cap/E = %.6e   vs floor = %.6e   ratio %.6f"
          %(float(worst[6]),float(r['floor']),float(worst[6]/r['floor'])))
    print("  margin over all sites: %.4f%% of the floor"
          %(100.0*float(max(s[6] for s in r['sites'])/r['floor'])))

print(); print("="*78); print("Q2.  do survivors need phase classes?"); print("="*78)
if not r['survivors']:
    print("  No survivors, so no grouping is required.  The coarse test alone removed")
    print("  %d of %d candidates and the refined memory bound removed the rest."
          %(r['coarse'],r['coarse']+len(r['sites'])))
    chk("question Q2 is vacuous on this shell", True)
else:
    print("  %d survivors; grouping would be needed."%len(r['survivors']))
    chk("survivors exist -> grouping needed", False)

print(); print("="*78); print("Q3.  is another of the 474 shells easier?"); print("="*78)
x,cf=AL,[]
for _ in range(70):
    i=int(x); cf.append(i)
    if x==i: break
    x=1/(x-i)
p0,q0,p1,q1=1,0,cf[0],1
conv=[(cf[0],1)]
for t in cf[1:]:
    p0,q0,p1,q1=p1,q1,t*p1+p0,t*q1+q0
    conv.append((p1,q1))
cands={}
for i in range(20,len(conv)-1):
    p,q=conv[i]; pp,qq=conv[i+1]
    for a in range(1,60):
        for b in range(1,60):
            Dx=a*q+b*qq; Kx=a*p+b*pp
            if Dx<10**20 or Dx>10**24: continue
            if Kx%8!=1: continue
            nb=(Kx-17)//8
            if 8+8*nb!=Kx-9: continue
            if not (5<F(8*Dx,Kx)<6): continue
            if 8*Dx//Kx+2!=7 or (Kx-9)*Dx//Kx+2!=Dx-4: continue
            lam=Kx*lo2-Dx*hi3
            if lam<=0: continue
            sixes=Dx-11-5*nb
            m=(sixes-PAST-FW-1+labels-1)//labels
            need=(F("0.162")*F(Dx,3)/lam+F(256,81)*F(Dx,3))/2**L
            if m<163 or need>=SIGMA: continue
            hb=float(Kx*(Kx*hi2-Dx*lo3)/lo3)
            if hb>10**6: continue
            cands[Dx]=(Kx,m,float(need),hb)
print("  shells passing the full cheap screen: %d"%len(cands))
best=sorted(cands.items(),key=lambda kv:kv[1][3])[:8]
print("  %24s %10s %14s %12s %s"%("D","m","margin","h-sites","complete?"))
done=[]
for Dx,(Kx,m,nd,hb) in best:
    rr=analyse(Kx,Dx,"scan")
    comp = rr is not None and rr.get('survivors') is not None and len(rr['survivors'])==0
    if comp: done.append((Dx,Kx,m,SIGMA/nd,len(rr['sites'])))
    print("  %24d %10d %14.1fx %12.1f %s"%(Dx,m,SIGMA/nd,hb,"YES" if comp else "no"))
chk("at least one shell admits a COMPLETE first-crossing certificate", bool(done),
    "%d complete of %d examined"%(len(done),len(best)))
if done:
    e=min(done,key=lambda z:z[4])
    print()
    print("  EASIEST COMPLETE SHELL (fewest sites, all excluded):")
    print("    D = %d"%e[0]); print("    K = %d"%e[1])
    print("    m = %d, span margin %.1fx, early-crossing sites %d -- ALL excluded"
          %(e[2],e[3],e[4]))
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
