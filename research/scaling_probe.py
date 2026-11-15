"""Phases 1/2/5: what scales, candidate shells, and the first obstruction."""
from fractions import Fraction as F
from decimal import Decimal, getcontext
import math
getcontext().prec=240
AL=Decimal(3).ln()/Decimal(2).ln()
ALF=F(int(AL*Decimal(10)**200),10**200)
LN2=F(int(Decimal(2).ln()*Decimal(10)**200),10**200)
LN3=F(int(Decimal(3).ln()*Decimal(10)**200),10**200)
K0=1412737287532845639425; D0=891337988684307508298

print("="*78); print("PHASE 1.  eta is a property of log_2 3, not of the shell"); print("="*78)
eta0=F(8*D0,K0)-5
print("  eta = 8D/K - 5, and D/K = log_3 2 - Lambda/(K ln3), so")
print("     eta = 8 log_3 2 - 5 - 8 Lambda/(K ln3)  ->  8/alpha - 5 = %.9f"%float(8/AL-5))
print("  on the 21/7 shell eta = %.9f  (difference %.3e)"%(float(eta0),float(8/AL-5)-float(eta0)))
cap=int(1.0/float(8/AL-5))-2
print("  The six-one separation needs (past + F + 2) * eta < 1, i.e.")
print("     past + F  <=  floor(1/eta) - 2  =  %d .   UNIVERSAL across shells."%cap)
print("  The 21/7 certificate uses past=7, F=12: %d = %d, EXACTLY at the cap."%(7+12,cap))

print(); print("="*78); print("PHASE 5.  the criterion with K cancelled"); print("="*78)
print("  Contradiction needs floor > Q, i.e. (2^L sigma - 1.054 D)/0.162 > D/(3 Lambda).")
print("  With sigma_m ~ m 3^t/|R|, m ~ 0.00593 K/(7 s^F), 2^L = 256^(F+1),")
print("  3^t/|R| ~ (3^5/s)^past, and D ~ 0.631 K, the K cancels and the criterion is")
print("     (256/s)^F * (243/s)^past * Lambda  >  const .")
s=21
print("  with s=21: (256/21)^F = %.4f^F and (243/21)^past = %.4f^past"%(256/s,243/s))
best=None
for Fw in range(1,cap+1):
    pa=cap-Fw
    if pa<1: continue
    v=(256/s)**Fw*(243/s)**pa
    if best is None or v>best[0]: best=(v,Fw,pa)
print("  maximised under past + F <= %d at F=%d, past=%d, value %.4e"%(cap,best[1],best[2],best[0]))
print("  the 21/7 choice F=12, past=7 gives %.4e (%.1f%% of the optimum)"
      %((256/s)**12*(243/s)**7, 100*((256/s)**12*(243/s)**7)/best[0]))
print()
print("  THE EARLY-CROSSING ENUMERATION is the quantity that grows.  The h-range is")
print("     -1 < h < K*Lambda/ln3 ,  so the number of candidate sites is ~ K*Lambda/ln3.")
L0=K0*LN2-D0*LN3
print("  on the 21/7 shell K*Lambda/ln3 = %.4f -- only 131 candidates."%float(K0*L0/LN3))
print("  Lambda cannot be made small for free: the seed criterion needs")
print("     Lambda > const / ((256/s)^F (243/s)^past)  ~  %.3e ,"%(0.157/best[0]))
print("  so  K*Lambda  >=  K * %.3e  and the enumeration grows LINEARLY in K."%(0.157/best[0]))
print("  %12s %16s %20s"%("K","min K*Lambda/ln3","h-enumeration"))
for Kx in (1e19,1e21,1e24,1e27,1e30,1e33,1e36):
    hc=Kx*0.157/best[0]/float(LN3)
    print("  %12.0e %16.3e %20s"%(Kx,hc,"feasible" if hc<1e9 else "INFEASIBLE"))
print("  => FIRST UNAVOIDABLE OBSTRUCTION: past + F <= %d caps the gain at %.3e, so the"%(cap,best[0]))
print("     h-enumeration cannot be traded away and becomes infeasible around K ~ 1e30.")
print("     The cap comes from eta = 8 log_3 2 - 5, a property of log_2 3 alone.")

print(); print("="*78); print("PHASE 2.  candidate near-critical shells"); print("="*78)
print("  Need K = ceil(D alpha), 2^K > 3^D (automatic), a legal 5-or-6 schedule, enough")
print("  six-one bytes, and K*Lambda small enough to enumerate early crossings.")
print("  Lambda = (ceil(D alpha) - D alpha) ln2, so near-critical shells are the UPPER")
print("  approximants of alpha.  Build them from the continued fraction:")
x,cf=AL,[]
for _ in range(60):
    i=int(x); cf.append(i)
    if x==i: break
    x=1/(x-i)
p0,q0,p1,q1=1,0,cf[0],1
conv=[(cf[0],1)]
for t in cf[1:]:
    p0,q0,p1,q1=p1,q1,t*p1+p0,t*q1+q0
    conv.append((p1,q1))
print("  %22s %22s %12s %14s %14s %s"%("D","K=ceil(D alpha)","log10 K","K*Lam/ln3","six-one bytes","usable")) 
cands=[]
for (p,q) in conv:
    if q<10**15 or q>10**32: continue
    Kx=int(F(q)*ALF)+1
    lam=F(Kx)*LN2-F(q)*LN3
    if lam<=0: continue
    hcount=float(Kx*lam/LN3)
    nb=(Kx-17)//8
    sixes=q-11-5*nb
    ok = (5*Kx < 8*q*1 and True)
    legal = 5 < 8*q/Kx < 6
    enough = sixes > 2*7*21**12
    usable = legal and enough and hcount<1e9 and sixes>0
    cands.append((q,Kx,hcount,sixes,usable))
    print("  %22d %22d %12.2f %14.4e %14.4e %s"
          %(q,Kx,math.log10(Kx),hcount,float(sixes),"YES" if usable else "no"))
print("  (the 21/7 shell is D = %d)"%D0)
print("  its D a convergent denominator? %s"%any(q==D0 for _,q in conv))
good=[c for c in cands if c[4]]
print("  usable candidates: %d"%len(good))
if good:
    g=max(good,key=lambda c:c[1])
    print("  strongest (largest K with a feasible enumeration): D = %d, K = %d"%(g[0],g[1]))
    print("     K*Lambda/ln3 = %.4e, six-one bytes = %.4e"%(g[2],g[3]))
