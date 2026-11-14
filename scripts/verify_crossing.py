"""Exact arithmetic checks for the first-crossing note.

Every pass/fail decision uses integer or Fraction arithmetic. Floats in the
first sections format diagnostic output only. The high-range reciprocal-sum
accounting and the published verification and logarithmic-form theorems are
imported mathematical inputs. Section 12 recomputes the finite count itself.
"""
import sys
import time
from bisect import bisect_left
from collections import defaultdict
from fractions import Fraction as F
from logbounds import L2, U2, L3, U3

t0 = time.monotonic()
fails = []
check_count = 0
def chk(name, cond, extra=""):
    global check_count
    check_count += 1
    print("  [%s] %s%s" % ("PASS" if cond else "FAIL", name, ("  -- " + extra) if extra else ""))
    if not cond: fails.append(name)

B_RT = 28 * 10**18                      # Rozier-Terracol Cor. 5.4 verification bound
a, b = 9809721694, 6189245291           # convergent 20 of log2 3   (a/b < alpha)
u, v = 630138897, 397573379             # convergent 19 of log2 3   (u/v > alpha)

print("External mathematical inputs: Rozier--Terracol, Cor. 5.4 (verification")
print("through 28e18) and Prop. 6.3 (logarithmic form); NFB height record (34)")
print("through D for the gap above 28e18. These external results and")
print("the NFB finite-path accounting in section 12 are not reproved by this script.")
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
print("   These are the two integers of Eliahou's cycle-length work.")

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
observed_crossings = 0
for n in range(3, 30000, 2):
    r = first_crossing(n)
    if not r: continue
    observed_crossings += 1
    k, d, bits, xk = r
    C = 0; dd = 0
    for j, bit in enumerate(bits):
        if bit: dd += 1; C += 3 ** (d - dd) * 2 ** j
    M = F(C, 3 ** d)
    s = sum(1 for j in range(k) if bits[j] and (j == 0 or not bits[j - 1]))
    if not M < s: ok_s = False
    if not M <= F(d, 3): ok_d3 = False
chk("all tested odd starts 3 <= n < 30000 reached a first crossing",
    observed_crossings == len(range(3, 30000, 2)))
chk("M_w < s on these first crossings", ok_s)
chk("M_w <= d/3 on these first crossings", ok_d3,
    "this is the sharper bound the note's Farey step uses")

print()
print("=== 4. the run-length realizer bound and G_75 ===")
chk("U3/L2 < 317/200", U3 / L2 < F(317, 200), "alpha = log2 3 < 1.585")
G75 = sum(F(317, 200) ** h for h in range(75))
chk("G_75 <= sum_{h<75} (317/200)^h < 1.72e15", G75 < F(172, 100) * 10**15,
    "rational upper bound = %.6e" % float(G75))

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
    14 < 7 + 14 * (j + 2), "the affine bracket grows by 14, less than its base value")

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
    0 < F(8, 9) < 1, "positive powers of a number in (0,1) remain in (0,1)")
chk("h >= 64: 8^64 - 5 > 5 * 256^14", 8**64 - 5 > 5 * 256**14)
chk("and ((h+1)/h)^14 < 8 for h >= 64", F(65,64)**14 < 8, "1 + 1/h decreases for h >= 64")

print()
print("=== 9. the worked example ===")
k, d, bits, xk = first_crossing(59)
chk("59 has parity word 1101100 and first crossing at k=7",
    k == 7 and "".join(map(str, bits)) == "1101100")
chk("T^7(59) = 38, i.e. 59 - T^7(59) = 21", xk == 38)
chk("59 is the least realizer of 1101100", least_realizer_parity([1,1,0,1,1,0,0]) == 59)
def realizer_class_parity(bits):
    """Invert the affine branch one parity bit at a time, from the last input."""
    residue, modulus = 0, 1
    for bit in reversed(bits):
        if bit not in (0, 1):
            raise ValueError("parity bits must be zero or one")
        modulus *= 2
        if bit:
            residue = (2 * residue - 1) * pow(3, -1, modulus) % modulus
        else:
            residue = 2 * residue % modulus
    return residue, modulus

crossing = [1, 1, 0, 1, 1, 0, 0]
valuation_inputs = crossing + [1]
valuation_and_odd_endpoint = valuation_inputs + [1]
for word, expected in (
        (crossing, (59, 128)),
        (valuation_inputs, (187, 256)),
        (valuation_and_odd_endpoint, (187, 512))):
    chk("raw word %s has class %d mod %d" %
        ("".join(map(str, word)), *expected),
        realizer_class_parity(word) == expected and
        [x for x in range(expected[1]) if orbit_bits(x, len(word)) == word]
        == [expected[0]])
chk("the other endpoint parity gives 443 mod 512 instead",
    realizer_class_parity(valuation_inputs + [0]) == (443, 512))

print()
print("=== 10. the initial-square corollary (Beyond Sturmian, finite-word form) ===")
print("   For an initial square UU with |U| = l, x = T^l(n) has the same")
print("   l-bit prefix as n. The proper-prefix multiplier gives x > n.")
lg43_lo = (2 * L2 - U3) / U2          # lower bound for log2(4/3)
lg43_hi = (2 * U2 - L3) / L2
chk("log2(4/3) enclosed in (%.9f, %.9f)" % (float(lg43_lo), float(lg43_hi)),
    lg43_lo < lg43_hi and lg43_lo > F(415,1000) and lg43_hi < F(416,1000))
chk("14.3 / log2(4/3) < 35", F(143,10) / lg43_lo < 35,
    "14.3/log2(4/3) = %.6f" % float(F(143,10) / lg43_lo))
# Regression over all observed first crossings, regardless of descent.
ok_sq = ok_dividend = True
square_samples = 0
for n in range(3, 60000, 2):
    crossing_data = first_crossing(n)
    if crossing_data is None:
        continue
    k, d, bits, xk = crossing_data
    for length in range(1, k // 2 + 1):
        if bits[:length] != bits[length:2 * length]:
            continue
        square_samples += 1
        carry, odd_factor = 0, 1
        for j, bit in enumerate(bits[:length]):
            if bit:
                carry = 3 * carry + (1 << j)
                odd_factor *= 3
        scale = 1 << length
        dividend = (odd_factor - scale) * n + carry  # = 2^length (T^length(n) - n)
        if dividend <= 0 or dividend % (scale * scale) != 0:
            ok_dividend = False
        if 4**length >= (n + 1) * 3**length:
            ok_sq = False
chk("strictly positive dividend divisible by 2^(2l) for all observed initial squares",
    square_samples > 0 and ok_dividend, "%d squares for odd starts 3 <= n < 60000" % square_samples)
chk("(4/3)^l < n+1 for all observed initial squares", square_samples > 0 and ok_sq)

print()
print("=== 11. the h = 12 regression (prefix floor vs completed-word least seed) ===")
h = 12
pref = [1,1,0] * h
# The repeated prefix has 24 odd steps. Its coefficient stays above one
# through each full block, then drops below one during the appended zeros.
k12 = (9 ** h).bit_length()
bits12 = pref + [0] * (k12 - len(pref))
floor_pref = 8 ** h - 5
chk("the (110)^12 prefix class is 68719476731 mod 2^36",
    floor_pref == 68719476731 and realizer_class_parity(pref) == (floor_pref, 2**36))
chk("the completed word has length 39, 24 odd steps and three appended zeros",
    (len(bits12), sum(bits12), k12 - len(pref)) == (39, 24, 3))
r12, mod12 = realizer_class_parity(bits12)
chk("the completed word has least seed 343597383675 modulo 549755813888",
    (r12, mod12) == (343597383675, 549755813888))
chk("the completed seed first crosses at 39, reaching 176518460300",
    first_crossing(r12) == (39, 24, bits12, 176518460300))
chk("the prefix least seed realizes the 36-bit prefix only, not its 39-bit extension",
    orbit_bits(floor_pref, 36) == pref and orbit_bits(floor_pref, 39) != bits12)

print()
print("=== 12. high-range reciprocal-sum cutoffs ===")
print("   Imported NFB hypotheses: distinct inputs > 28019077177231758495;")
print("   the endpoint is excluded. Imported analytic bounds: block weight")
print("   15/(4x), upper tail < 10^-6 and boundary contribution < 10^-15.")
print("   The Appendix 1D counting rule is recomputed below in integer arithmetic.")

def rounded_bracket_sum():
    """Equation (164), NFB Appendix 1D: recompute counts in sixteenths."""
    grid = [mantissa << (exponent - 2)
            for exponent in range(3, 24) for mantissa in range(64, 128)]
    overflow = len(grid)
    even = [bisect_left(grid, value // 2) if value >= 256 else -1 for value in grid]
    odd = [bisect_left(grid, 3 * (value + 16) // 2) for value in grid]
    # Both successors of overflow are retained, without further tests.
    even.append(overflow)
    odd.append(overflow)
    states = {(odd[0], 1): 1}  # A composition starts with an odd input.
    numerator = 0
    sample_counts = {}
    for length in range(1, 1001):
        if length >= 61:
            count = sum(multiplicity for (_, ones), multiplicity in states.items()
                        if ones == 38)
            numerator = 2 * numerator + count
            if length in (61, 62, 63, 100):
                sample_counts[length] = count
        if length == 1000:
            break
        successors = defaultdict(int)
        for (index, ones), multiplicity in states.items():
            if even[index] >= 0:
                successors[even[index], ones] += multiplicity
            successors[odd[index], min(ones + 1, 38)] += multiplicity
        states = successors
    return F(numerator, 2**1000), sample_counts

certificate_numerator = int(
    "250886484124983267713084036235929994106679131859883785503469464927892276700327"
    "399772416033700119525575246775531566546874139868602521988656701427729384118728"
    "585359221872960402072436733673462594656158947975960371636692695132559493005745"
    "097798633776451644036826297946593619831371395774140940847628582223")
finite_sum, samples = rounded_bracket_sum()
chk("rounded count agrees with independently reported Table 12 checkpoints",
    samples == {61: 9397196782201115, 62: 17959590990646997, 63: 34603440856372965,
                100: 894278252354101043785255555})
chk("the full finite count equals the reported exact certificate N/2^998",
    finite_sum == F(certificate_numerator, 2**998))
count_ceiling = F(93657291203213, 10**15)
chk("the full finite count is below the manuscript ceiling c, and c < 3/32",
    0 < finite_sum < count_ceiling < F(3, 32))

# Each bracket permits eight of the twelve starts in its shell:
# one residue class of the multipliers b=12,...,23 modulo 3 is absent.
spacing = max(sum(F(1, b) for b in range(12, 24) if b % 3 != excluded)
              for excluded in range(3))
chk("maximum retained reciprocal sum is 6719267/13693680",
    spacing == F(6719267, 13693680) and spacing < F(1, 2))
remainder = F(1, 10**6) + F(1, 10**15)
coarse = F(45, 256) + remainder
refined = F(15, 4) * spacing * count_ceiling + remainder
chk("the refined high-range reciprocal bound is strictly smaller than the coarse bound",
    0 < refined < coarse < F(9, 50))

def exp_negative_bounds(t):
    """For 0<t<1, the odd/even partial sums of exp(-t) enclose its value."""
    if not 0 < t < 1:
        raise ValueError("alternating enclosure requires 0 < t < 1")
    partial = term = F(1)
    for order in range(1, 11):
        term *= -t / order
        partial += term
    return partial - t * term / 11, partial  # orders 11 and 10

for name, bound, digits, width_percent in (
        ("coarse", coarse, 9430895072211041, F("8.45334")),
        ("refined", refined, 9441733669279807, F("8.28763"))):
    lower, upper = exp_negative_bounds(bound / 3)
    printed_lo = F(digits, 10**16)
    printed_hi = F(digits + 1, 10**16)
    chk("%s cutoff lies strictly between 0.%016d and 0.%016d" %
        (name, digits, digits + 1),
        printed_lo < lower < upper < printed_hi)
    # Percentage of the unit interval in -log_2(rho) left by this bound.
    width_lo, width_hi = 100 * bound / (3 * U2), 100 * bound / (3 * L2)
    half_unit = F(1, 2 * 10**5)
    chk("%s surviving and excluded log widths round to %.5f%% and %.5f%%" %
        (name, float(width_percent), float(100 - width_percent)),
        width_percent - half_unit < width_lo <= width_hi < width_percent + half_unit)
    if name == "refined":
        threshold = F(9441733669, 10**10)
        chk("0.9441733669 rounds DOWN and is safe with inclusive rho <= threshold",
            threshold < lower and upper < threshold + F(1, 10**10))

print()
print("%d checks, %d passed, %d failed [runtime %.1f s]" %
      (check_count, check_count - len(fails), len(fails), time.monotonic() - t0))
sys.exit(1 if fails else 0)
