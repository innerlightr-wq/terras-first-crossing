"""Verifier for the 81-run theorem.  Self-contained; uses only the certified
logarithm enclosures of logbounds.py.  Leaves verify_crossing.py untouched.

    python3 scripts/verify_81_runs.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from logbounds import L2, U2, L3, U3
from decimal import Decimal, getcontext

DEN = 10**80
LN2_LO = F(int(L2*DEN), DEN);   LN2_HI = F(int(U2*DEN)+1, DEN)
LN3_LO = F(int(L3*DEN), DEN);   LN3_HI = F(int(U3*DEN)+1, DEN)
assert LN2_LO <= L2 and U2 <= LN2_HI and LN3_LO <= L3 and U3 <= LN3_HI
AL_HI = F(317, 200)
assert LN3_HI/LN2_LO < AL_HI, "alpha < 317/200"
B, S = 28*10**18, 81
nfail = npass = 0
def chk(name, cond, extra=""):
    global nfail, npass
    if cond: npass += 1
    else: nfail += 1
    print("  [%s] %s %s" % ("PASS" if cond else "FAIL", name, extra))
def G(j):
    tot, p = F(0), F(1)
    for _ in range(j): tot += p; p *= AL_HI
    return tot
def clog2(x):
    h = 0
    while F(2)**h < x: h += 1
    return h
def e_br(p, q): return q*LN3_LO - p*LN2_HI, q*LN3_HI - p*LN2_LO

print("=== convergents of alpha = log2 3, with certified e = q ln3 - p ln2 > 0 ===")
getcontext().prec = 200
aa = Decimal(3).ln()/Decimal(2).ln()
x, cf = aa, []
for _ in range(90):
    i = int(x); cf.append(i)
    if x == i: break
    x = 1/(x - i)
p0,q0,p1,q1 = 1,0,cf[0],1
conv = [(cf[0],1)]
for t in cf[1:]:
    p0,q0,p1,q1 = p1,q1,t*p1+p0,t*q1+q0
    conv.append((p1,q1))
below = [(p,q) for (p,q) in conv if q > 0 and e_br(p,q)[0] > 0 and q < 10**40]
chk("alpha enclosed, %d admissible convergents p/q < alpha" % len(below), len(below) >= 20)

print("=== lower bound on d (pair I of the manuscript) ===")
PI = (9881527843552324, 6234549927241963)
chk("pair I is a convergent of alpha", PI in conv)
e1lo, e1hi = e_br(*PI)
chk("e_I in (5.22,5.24)e-18", F(522,100)*10**-18 < e1lo and e1hi < F(524,100)*10**-18,
    "(%.6e, %.6e)" % (float(e1lo), float(e1hi)))
dmin = (LN2_LO - F(PI[1]*S, B))/e1hi
chk("d > 1.28e17", dmin > F(128,100)*10**17, "d >= %.6e" % float(dmin))

print("=== the tail d >= 2^65 ===")
G81 = G(81)
chk("log2(81) <= 7", F(2)**7 >= 81)
chk("G_81 <= 2.7241e16", G81 < F(27241,10000)*10**16, "G_81 <= %.6e" % float(G81))
bad = [j for j in range(65, 500) if not G81*(7 + F(133,10)*(j+2)) <= F(2)**j]
chk("G_81[7+13.3(j+2)] <= 2^j for all 65 <= j <= 499", not bad, str(bad[:3]))
chk("j = 64 fails, so 2^65 is the least admissible split",
    not (G81*(7+F(133,10)*66) <= F(2)**64))

print("=== the head d < 2^65: five-row bootstrap ===")
def sel(lam, r):
    den = lam - F(r, B)
    if den <= 0: return None
    return 1 + F(S-r)*(1 + F(81,B))/den
def rowbound(d_prev, r, pq):
    p, q = pq
    elo, ehi = e_br(p, q)
    if elo <= 0: return None
    num = LN2_LO - d_prev*ehi
    if num <= 0: return None
    lam = num/q
    if lam <= F(r, B): return None
    cand = [1 + F(81)*(1+F(81,B))/lam]
    sb = sel(lam, r)
    if sb is not None: cand.append(sb)
    H = clog2(max(cand))
    return lam, H, F(H)*(G(S-r)+r)
ROWS = [(0, 75, 2044*10**15), (2, 72, 781*10**15), (3, 72, 493*10**15),
        (4, 74, 320*10**15), (7, 68, 738*10**14)]
d_prev = F(2)**65
print("  %3s %3s %4s %15s %17s %17s" % ("row","r","H","lambda","d bound","target"))
for i, (r, Ht, dt) in enumerate(ROWS):
    best = None
    for pq in below:
        res = rowbound(d_prev, r, pq)
        if res is None: continue
        if best is None or res[2] < best[2]: best = res + (pq,)
    chk("row %d (r=%d) has an admissible convergent" % (i, r), best is not None)
    if best is None: break
    lam, H, dnew, pq = best
    print("  %3d %3d %4d %15.6e %17.6e %17.6e   (p,q)=(%d,%d)"
          % (i, r, H, float(lam), float(dnew), float(dt), pq[0], pq[1]))
    chk("row %d: H = %d <= target %d and bound <= target" % (i, H, Ht),
        H <= Ht and dnew <= F(dt))
    d_prev = min(d_prev, dnew)
chk("CONTRADICTION: final d bound < d lower bound", d_prev < dmin,
    "%.6e < %.6e" % (float(d_prev), float(dmin)))
print()
print("%d checks, %d passed, %d failed" % (npass+nfail, npass, nfail))
sys.exit(1 if nfail else 0)
