"""Bottleneck analysis and parameter optimisation for the 75-run theorem.

The published proof fixes s=75, split J=61, and two Diophantine pairs.  Here we treat
(s, J, pair_I, pair_II) as free and ask for the largest s the SAME lemmas can reach.

Structure of the argument, for a hypothetical first-crossing failure with s odd runs:
  seed bound           n > B = 2.8e19                            (Rozier-Terracol Cor 5.4)
  comparator           0 < Lambda = k ln2 - d ln3 < s/n < s/B
  run-length lemma     d <= G_s * log2(n+1),  G_s = sum_{j<s} alpha^j,  alpha = log2 3
  pair (p,q), p/q<a    e := q ln3 - p ln2 > 0  and  q*Lambda + d*e >= ln2

Case d >= 2^J : Rhin gives Lambda >= k^-13.3, so n < s(2d)^13.3 and
                log2(n+1) < log2 s + 13.3(j+2) for 2^j <= d < 2^{j+1}.
                Need G_s*(log2 s + 13.3(j+2)) <= 2^j for every j >= J.
Case d <  2^J : pair II gives Lambda >= (ln2 - 2^J e2)/q2 =: Lmin, hence n < s/Lmin;
                pair I  gives d >= (ln2 - q1 s/B)/e1 =: dmin.
                Need G_s * ceil(log2(s/Lmin + 1)) < dmin.

Every inequality below is decided in exact rational arithmetic, with ln2 and ln3 replaced
by the certified enclosures of scripts/logbounds.py and used in the conservative direction.
The continued fraction only SUPPLIES candidate pairs; each candidate's properties are then
verified exactly, so the search needs no trust.
"""
import sys
sys.path.insert(0, 'scripts')
from fractions import Fraction as F
from logbounds import L2, U2, L3, U3
from decimal import Decimal, getcontext

B = 28 * 10**18
# The series enclosures have ~200-digit denominators; raising them to the 75th power is
# ruinous.  Round them OUTWARD to 20-decimal rationals, which stays certified because
# l2 <= L2 < ln2 < U2 <= u2, and verify that rounding here.
D = 10**60
l2 = F(int(L2 * D), D);        u2 = F(int(U2 * D) + 1, D)
l3 = F(int(L3 * D), D);        u3 = F(int(U3 * D) + 1, D)
assert l2 <= L2 and U2 <= u2 and l3 <= L3 and U3 <= u3
L2, U2, L3, U3 = l2, u2, l3, u3
AL_HI = F(317, 200)                       # alpha < 317/200, verified against the enclosure
assert U3 / L2 < AL_HI

def G(s):                                 # certified UPPER bound for G_s
    tot, p = F(0), F(1)
    for _ in range(s):
        tot += p; p *= AL_HI
    return tot

def e_bounds(p, q):                       # certified bracket for e = q ln3 - p ln2
    return q * L3 - p * U2, q * U3 - p * L2

def ceil_log2(x):                         # least integer L with 2^L >= x, exact
    L = 0
    while F(2) ** L < x: L += 1
    return L

def case2_ok(s, J, jmax=400):
    """for all j >= J: G_s*(log2 s + 13.3(j+2)) <= 2^j   (13.3 and log2 s bounded above)"""
    gs = G(s); ls = ceil_log2(F(s))       # log2 s <= ls
    for j in range(J, jmax):
        if not gs * (ls + F(133, 10) * (j + 2)) <= F(2) ** j:
            return False, j
    return True, None

def case1_ok(s, J, pI, pII):
    p1, q1 = pI; p2, q2 = pII
    e1_lo, e1_hi = e_bounds(p1, q1)
    e2_lo, e2_hi = e_bounds(p2, q2)
    if e1_lo <= 0 or e2_lo <= 0: return False, "pair not below alpha"
    num2 = L2 - F(2) ** J * e2_hi                       # ln2 - 2^J e2, lower bound
    if num2 <= 0: return False, "pair II too weak for this split"
    Lmin = num2 / q2                                    # lower bound on Lambda
    Lmax = ceil_log2(F(s) / Lmin + 1)                   # upper bound on log2(n+1)
    num1 = L2 - F(q1 * s, B)                            # ln2 - q1 s/B, lower bound
    if num1 <= 0: return False, "pair I too large for this s"
    dmin = num1 / e1_hi                                 # lower bound on d
    return (G(s) * Lmax < dmin), dict(Lmax=Lmax, dmin=float(dmin), lhs=float(G(s) * Lmax))

# ---- convergents of alpha, as candidate pairs (verified exactly afterwards) ----
getcontext().prec = 400
a = Decimal(3).ln() / Decimal(2).ln()
x, cf = a, []
for _ in range(120):
    i = int(x); cf.append(i)
    if x == i: break
    x = 1 / (x - i)
p0, q0, p1_, q1_ = 1, 0, cf[0], 1
conv = [(cf[0], 1)]
for t in cf[1:]:
    p0, q0, p1_, q1_ = p1_, q1_, t * p1_ + p0, t * q1_ + q0
    conv.append((p1_, q1_))
below = [(p, q) for (p, q) in conv if q > 0 and q * L3 - p * U2 > 0 and q < 10**20]

print("=== published configuration reproduced ===")
PI_PUB = (9881527843552324, 6234549927241963)
PII_PUB = (6724555128221608268, 4242721909926539673)
for s in (75, 76):
    ok1, info1 = case1_ok(s, 61, PI_PUB, PII_PUB)
    ok2, jbad = case2_ok(s, 61)
    print("  s=%d  split 2^61   case1=%-5s  case2=%-5s  %s"
          % (s, ok1, ok2, info1 if isinstance(info1, dict) else info1))

print()
print("=== bottleneck: how each side moves with s ===")
print("  %4s %14s %8s %16s %16s" % ("s", "G_s", "L_max", "G_s*L_max", "d_min(pair I)"))
for s in (70, 75, 76, 80, 85, 91):
    ok1, info = case1_ok(s, 61, PI_PUB, PII_PUB)
    if isinstance(info, dict):
        print("  %4d %14.4e %8d %16.4e %16.4e" % (s, float(G(s)), info['Lmax'], info['lhs'], info['dmin']))
    else:
        print("  %4d %14.4e   -- %s" % (s, float(G(s)), info))

print()
print("=== optimisation over (J, pair I, pair II), maximising s ===", flush=True)
# The two pairs enter through disjoint quantities, so they optimise INDEPENDENTLY:
#   pair I  maximises dmin = (ln2 - q1 s/B)/e1      (depends on s, not J)
#   pair II maximises Lmin = (ln2 - 2^J e2)/q2      (depends on J, not s)
cands = [(p, q) for (p, q) in below if 10**3 < q < 10**19]
ehi = {pq: e_bounds(*pq)[1] for pq in cands}

def best_pairI(s):
    bst = (F(0), None)
    for (p, q) in cands:
        num = L2 - F(q * s, B)
        if num <= 0: continue
        v = num / ehi[(p, q)]
        if v > bst[0]: bst = (v, (p, q))
    return bst

def best_pairII(J):
    bst, tw = (F(0), None), F(2) ** J
    for (p, q) in cands:
        num = L2 - tw * ehi[(p, q)]
        if num <= 0: continue
        v = num / q
        if v > bst[0]: bst = (v, (p, q))
    return bst

best = (0, None)
for J in range(50, 80):
    Lmin, pII = best_pairII(J)
    if pII is None: continue
    for s in range(70, 140):
        if not case2_ok(s, J)[0]: continue
        dmin, pI = best_pairI(s)
        if pI is None: continue
        Lmax = ceil_log2(F(s) / Lmin + 1)
        if G(s) * Lmax < dmin and s > best[0]:
            best = (s, (J, pI, pII, Lmax, float(dmin), float(G(s) * Lmax)))
    print("   J=%2d done, best so far s=%d" % (J, best[0]), flush=True)

print()
if best[1]:
    s, (J, pI, pII, Lmax, dmin, lhs) = best
    print("  LARGEST s for which BOTH cases close: %d  (published: 75)" % s)
    print("  split 2^%d;  pair I = (%d, %d);  pair II = (%d, %d)" % (J, pI[0], pI[1], pII[0], pII[1]))
    print("  case 1:  G_s*L_max = %.5e  <  d_min = %.5e   (L_max = %d)" % (lhs, dmin, Lmax))
    ok1, info = case1_ok(s, J, pI, pII); ok2, _ = case2_ok(s, J)
    print("  independent re-check: case1=%s case2=%s" % (ok1, ok2))
    ok1n, _ = case1_ok(s + 1, J, pI, pII); ok2n, _ = case2_ok(s + 1, J)
    print("  s+1 = %d fails: case1=%s case2=%s" % (s + 1, ok1n, ok2n))
else:
    print("  no configuration closed")
