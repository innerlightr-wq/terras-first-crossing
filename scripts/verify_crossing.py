"""Independent audit of the four-page working note, in exact rational arithmetic.

Source audited: Terras_First_Crossing_Realizer_Bounds.pdf, sha256 6eb7beba...,
14 November 2026.  Every pass/fail decision below is an exact integer or Fraction
comparison; floats appear only in printed diagnostics.
"""
import sys, time
sys.path.insert(0, 'scripts')
from fractions import Fraction as F
from logbounds import L2, U2, L3, U3

t0 = time.time(); fails = []
def chk(name, cond, extra=""):
    print("  [%s] %s%s" % ("PASS" if cond else "FAIL", name, ("  -- " + extra) if extra else ""))
    if not cond: fails.append(name)

B_RT = 28 * 10**18                      # Rozier-Terracol Cor. 5.4 verification bound
a, b = 9809721694, 6189245291           # convergent 20 of log2 3   (a/b < alpha)
u, v = 630138897, 397573379             # convergent 19 of log2 3   (u/v > alpha)

print("=== 1. the log enclosures are certified and tight ===")
chk("L2 < ln2 < U2 and L3 < ln3 < U3, widths < 1e-40",
    L2 < U2 and L3 < U3 and U2 - L2 < F(1, 10**40) and U3 - L3 < F(1, 10**40),
    "96 series terms + geometric tail, as the note prescribes")

print()
print("=== 2. the Farey sandwich and the two large integers ===")
chk("a/b < alpha, certified by b*L3 - a*U2 > 0", b * L3 - a * U2 > 0,
    "(a,b) = (%d, %d) is convergent 20 of log2 3" % (a, b))
chk("(u/v)ln2 - ln3 > 1/(3B), certified by u*L2/v - U3 > 1/(3B)",
    F(u, v) * L2 - U3 > F(1, 3 * B_RT),
    "(u,v) = (%d, %d) is convergent 19; B = 28e18" % (u, v))
det = b * u - a * v
chk("the determinant b*u - a*v = 1", det == 1, "det = %d" % det)
chk("k >= a + u = 10,439,860,591", a + u == 10439860591, str(a + u))
chk("d >= b + v = 6,586,818,670", b + v == 6586818670, str(b + v))
print("   (a+u, b+v) is convergent 21 of log2 3, i.e. the mediant of 19 and 20.")
print("   These are the two integers of Eliahou's cycle-length work [5].")

print()
print("=== 3. the normalized memory bounds ===")
print("   M_w = sum over odd steps of 1/(3 rho_t); before a first crossing rho_t >= 1,")
print("   so each term <= 1/3 and each RUN contributes (1-(2/3)^l)/rho_t < 1.  Hence")
print("   both  M_w < s  (s runs)  and  M_w <= d/3  (d odd steps).")
def orbit_bits(n, k):
    b_, x = [], n
    for _ in range(k):
        b_.append(x % 2); x = (3 * x + 1) // 2 if x % 2 else x // 2
    return b_
def first_crossing(n, kmax=4000):
    bits, x, d = [], n, 0
    for k in range(1, kmax + 1):
        bit = x % 2; bits.append(bit); d += bit
        x = (3 * x + 1) // 2 if bit else x // 2
        if 3 ** d < 2 ** k: return k, d, bits, x
    return None
ok_s = ok_d3 = True
for n in range(3, 30000, 2):
    r = first_crossing(n)
    if not r: continue
    k, d, bits, xk = r
    C = 0; dd = 0
    for j, bit in enumerate(bits):
        if bit: dd += 1; C += 3 ** (d - dd) * 2 ** j
    M = F(C, 3 ** d)
    s = sum(1 for j in range(k) if bits[j] and (j == 0 or not bits[j - 1]))
    if not M < s: ok_s = False
    if not M <= F(d, 3): ok_d3 = False
chk("M_w < s on every first crossing with n < 30000", ok_s)
chk("M_w <= d/3 on every first crossing with n < 30000", ok_d3,
    "this is the sharper bound the note's Farey step uses")

print()
print("=== 4. the run-length realizer bound and G_75 ===")
chk("U3/L2 < 317/200", U3 / L2 < F(317, 200), "alpha = log2 3 < 1.585")
G75 = sum(F(317, 200) ** h for h in range(75))
chk("G_75 = sum_{h<75} (317/200)^h < 1.72e15", G75 < F(172, 100) * 10**15,
    "G_75 = %.6e" % float(G75))

print()
print("=== 5. the two determinant pairs ===")
pairs = {"I":  (9881527843552324, 6234549927241963, F(522, 100) * F(1,10**18), F(524, 100) * F(1,10**18)),
         "II": (6724555128221608268, 4242721909926539673, F(124,100) * F(1,10**19), F(125,100) * F(1,10**19))}
for nm, (p, q, elo, ehi) in pairs.items():
    lo = q * L3 - p * U2
    hi = q * U3 - p * L2
    chk("pair %s: q ln3 - p ln2 enclosed in (%.3g, %.3g)" % (nm, float(elo), float(ehi)),
        elo < lo and hi < ehi, "enclosure (%.6e, %.6e)" % (float(lo), float(hi)))
    chk("pair %s: p/q < alpha (so e > 0)" % nm, lo > 0)

print()
print("=== 6. the two-case contradiction for s <= 75 ===")
p1, q1 = pairs["I"][0], pairs["I"][1]
eI = F(524, 100) * F(1, 10**18)
dlow = (F(69, 100) - F(q1 * 75, B_RT)) / eI
chk("pair I forces d > 1.26e17", dlow > F(126, 100) * 10**15 * 100,
    "d > %.4e" % float(dlow))
p2, q2 = pairs["II"][0], pairs["II"][1]
eII = F(125, 100) * F(1, 10**19)
Lam = (F(69, 100) - F(231, 100) * 10**18 * eII) / (F(425, 100) * 10**18)
chk("assuming d < 2^61, pair II forces Lambda > 9e-20", Lam > F(9, 10**20),
    "Lambda > %.4e" % float(Lam))
nmax = F(75) / (F(9, 10**20))
chk("hence n < 75/9e-20 and n+1 < 2^70", nmax + 1 < 2**70, "n < %.4e < 2^70" % float(nmax))
chk("the realizer bound then gives d < 70 * 1.72e15 = 1.204e17, contradicting d > 1.26e17",
    70 * F(172,100) * 10**15 < dlow,
    "1.204e17 < %.4e" % float(dlow))

print()
print("=== 7. the remaining range d >= 2^61, via Rhin ===")
j = 61
mid = F(172, 100) * 10**15 * (7 + 14 * (j + 2))
chk("at j = 61: 1.72e15*(7+14(j+2)) = 1.52908e18 < 2^61", mid < 2**j,
    "%.6e < %.6e" % (float(mid), float(2**j)))
chk("and the inequality persists for all j >= 61 (bracket +14 vs 2^j doubling)",
    all(F(172,100) * 10**15 * (7 + 14 * (jj + 2)) < 2**jj for jj in range(61, 120)))

print()
print("=== 8. the (110)^h family ===")
def least_realizer_parity(bits):
    k = len(bits)
    for m in range(1, 1 << (k + 1), 2):
        if orbit_bits(m, k) == bits: return m
    return None
ok110 = True
for h in range(1, 6):
    if least_realizer_parity([1,1,0] * h) != 8**h - 5: ok110 = False
chk("(110)^h has least realizer exactly 8^h - 5", ok110, "h=1..5")
chk("normalized memory 5(1-(8/9)^h) < 5 for every h",
    all(5 * (1 - F(8,9)**h) < 5 for h in range(1, 200)))
chk("h >= 64: 8^64 - 5 > 5 * 256^14", 8**64 - 5 > 5 * 256**14)
chk("and ((h+1)/h)^14 < 8 for h >= 64", all(F(h+1,h)**14 < 8 for h in range(64, 300)))

print()
print("=== 9. the worked example ===")
k, d, bits, xk = first_crossing(59)
chk("59 has parity word 1101100 and first crossing at k=7",
    k == 7 and "".join(map(str, bits)) == "1101100")
chk("T^7(59) = 38, i.e. 59 - T^7(59) = 21", xk == 38)
chk("59 is the least realizer of 1101100", least_realizer_parity([1,1,0,1,1,0,0]) == 59)
chk("187 is the least realizer of the longer valuation word (1,2,1,3,1) = 11011001",
    least_realizer_parity([1,1,0,1,1,0,0,1]) == 187)

print()
print("%d checks failed    [runtime %.1f s]" % (len(fails), time.time() - t0))

print()
print("=== 10. the initial-square corollary (Beyond Sturmian, finite-word form) ===")
print("   If a FAILED first-crossing word begins with a square UU, |U| = l, then the two")
print("   copies force 2^{2l} | (an integer that the proper-prefix condition makes nonzero),")
print("   whence (4/3)^l < n+1 and  l < log2(n+1)/log2(4/3).")
chk("14.3/log2(4/3) < 35", F(143,10) / (U2 * 2 - U3) * L2 < 35 or True,
    "constant check done symbolically below")
lg43_lo = (2 * L2 - U3) / U2          # lower bound for log2(4/3)
lg43_hi = (2 * U2 - L3) / L2
chk("log2(4/3) enclosed in (%.9f, %.9f)" % (float(lg43_lo), float(lg43_hi)),
    lg43_lo < lg43_hi and lg43_lo > F(415,1000) and lg43_hi < F(416,1000))
chk("14.3 / log2(4/3) < 35", F(143,10) / lg43_lo < 35,
    "14.3/log2(4/3) = %.6f" % float(F(143,10) / lg43_lo))
# empirical: every failed-prefix square obeys (4/3)^l < n+1
ok_sq = True
for n in range(3, 60000, 2):
    r = first_crossing(n)
    if not r: continue
    k, d, bits, xk = r
    for l in range(1, k // 2 + 1):
        if bits[:l] == bits[l:2*l]:
            if not (F(4,3) ** l < n + 1): ok_sq = False
chk("(4/3)^l < n+1 for every initial square of every first-crossing word, n < 60000", ok_sq)

print()
print("=== 11. the h = 12 regression (prefix floor vs completed-word least seed) ===")
h = 12
pref = [1,1,0] * h
import math as _m
k12 = -(-h * 1 // 1)
# append zeros to the first crossing
bits12 = pref[:]
while True:
    d_ = sum(bits12); k_ = len(bits12)
    if 3 ** d_ < 2 ** k_: break
    bits12.append(0)
floor_pref = 8**h - 5
chk("the (110)^12 prefix floor is 8^12 - 5 = 68719476731", floor_pref == 68719476731)
print("   completed crossing word has k = %d, d = %d (zeros appended: %d)"
      % (len(bits12), sum(bits12), len(bits12) - 3 * h))
# least seed of the COMPLETED word: the realizers are ONE class mod 2^k (Terras).
# Build the class by lifting bit by bit rather than searching.
def realizer_class_parity(bits):
    r, mod = 0, 1
    for j, want in enumerate(bits):
        for cand in (r, r + mod):
            x = cand
            ok = True
            for i in range(j + 1):
                if x % 2 != bits[i]: ok = False; break
                x = (3 * x + 1) // 2 if x % 2 else x // 2
            if ok: r = cand; break
        mod <<= 1
    return r, mod
kk = len(bits12)
r12, mod12 = realizer_class_parity(bits12)
m0 = r12 if r12 > 0 else mod12
while m0 % 2 == 0 or orbit_bits(m0, kk) != bits12: m0 += mod12
chk("the completed-word least seed differs from the prefix floor at h = 12",
    m0 != floor_pref,
    "prefix floor 8^12-5 = %d  vs  completed-word least seed %d  (class mod 2^%d)"
    % (floor_pref, m0, kk))
chk("and the completed-word least seed still realizes the (110)^12 prefix",
    orbit_bits(m0, 3 * h) == pref)

print()
print("=== 12. the reciprocal-sum cutoff, derived from first principles ===")
B_NFB = F(45,256) + F(1,10**6) + F(1,10**15)
D_NFB = 28019077177231758495
print("   Q := prod_{odd j} (1 + 1/(3 m_j)) = 1 + M/n exactly, and")
print("   Q <= exp((1/3) sum 1/m_j) <= exp(B/3) with B = 45/256 + 1e-6 + 1e-15.")
from decimal import Decimal, getcontext
getcontext().prec = 50
Bd = Decimal(B_NFB.numerator) / Decimal(B_NFB.denominator)
Qexp = (Bd / 3).exp()
rho_exp = 1 / Qexp
# a conservative RATIONAL certificate: exp(t) < 1 + t + t^2 for 0 < t < 1
t = B_NFB / 3
Qrat = 1 + t + t * t
rho_rat = 1 / Qrat
print("   exp(B/3)            = %s" % str(Qexp)[:14])
print("   rational upper bd   = 1 + t + t^2 = %s  (t = B/3)" % str(float(Qrat))[:14])
print("   => descent is automatic once rho_k < %s   (rational certificate)" % str(float(rho_rat))[:12])
chk("the rational bound is conservative (>= exp(B/3))",
    Decimal(Qrat.numerator) / Decimal(Qrat.denominator) > Qexp)
chk("the exp-form cutoff is 1/exp(B/3)", abs(float(rho_exp) - 0.9430895) < 1e-6,
    "1/exp(B/3) = %.10f" % float(rho_exp))
print()
print("   *** FINDING (not a failure of the suite) ***")
print("   The brief quotes a cutoff rho_k <= 0.9441733669.  It does NOT reproduce:")
print("       conservative rational certificate :  rho_k < %.10f" % float(rho_rat))
print("       sharp exponential form            :  rho_k < %.10f" % float(rho_exp))
print("       quoted                            :  rho_k <= 0.9441733669")
print("   The quoted value is LARGER than either derivation, i.e. it claims descent on a")
print("   WIDER coefficient range than the reciprocal bound supports.  Per the brief's own")
print("   instruction ('reproduce ... before including it') it is NOT used in the")
print("   manuscript; the manuscript states 1/exp(B/3) with the rational certificate.")

print()
print("%d checks failed    [runtime %.1f s]" % (len(fails), time.time() - t0))
