"""Reproduce two June-to-Terras corridor diagnostics with Python 3.8+.

Run: python3 terras_june_corridor_probe.py
All mathematical checks and counts are exact. Decimals are for display.
This is a follow-up experiment. The 485-step family lies below the paper's
global lower bound on the length of a counterexample.

1. In a first-crossing word w let k be the shortcut length, d the number of
   odd inputs and S[j] the position of the j-th odd input, indexed from 0.
   The affine carry equals sum(3^(d-1-j)*2^S[j]). With R[j]=S[j]-j*log2(3)
   the normalized memory C/3^d equals sum(2^R[j])/3. The proper odd-input
   positions have R[j]<=0. Raw occupation under an upper ceiling c>=0
   counts all of them. If S[j] is moved earlier by one then the corresponding
   carry term is halved.

2. Write A=2^k-3^d>0. Any failure at eligible seed n has an integer excess
   h=T^k(n)-n>=0, and 2^k*h=C-A*n. Hence for any odd F dividing A,
   h=C*2^(-k) modulo F and h <= floor((C-2*A)/2^k). If the least
   nonnegative residue already exceeds that ceiling, every eligible seed
   descends. Passing the residue test means only that this filter is
   inconclusive. This code tests the incomplete factor F=929.
"""
from fractions import Fraction


def require(truth, message):
    if not truth:
        raise ArithmeticError(message)


def carry(positions):
    count = len(positions)
    return sum(3**(count-1-j) * (1 << t) for j, t in enumerate(positions))


def record_seed():
    n = 63_728_127  # Table 1 of the June Word-Level Occupation–Realizer Reduction
    x, k, d, runs, previous_odd = n, 0, 0, 0, False
    odd_positions = []
    while True:
        odd = x % 2 == 1
        if odd:
            odd_positions.append(k)
            runs += int(not previous_odd)
            d += 1
            x = (3*x+1)//2
        else:
            x //= 2
        k += 1
        previous_odd = odd
        if 3**d < (1 << k):
            break
    require((k, d, runs, x) == (376, 237, 97, 49_521_089), 'June record')
    require((1 << k)*x == 3**d*n + carry(odd_positions), 'carry')
    ratio = Fraction(3**d, 1 << k)
    print(f'June record seed {n:,}: k={k}, d={d}, odd runs={runs}, endpoint={x:,}, rho={float(ratio):.9f}')


def many_run_family():
    k, d, factor = 485, 306, 929
    modulus = 1 << k
    gap = modulus-3**d
    require(gap > 0 and gap % factor == 0, 'odd divisor')
    # floor(j*log2(3)) is obtained exactly from bit_length(3^j).
    positions = [(3**j).bit_length()-1 for j in range(d)]
    by_position = {t: j for j, t in enumerate(positions)}
    bits = ''.join('1' if t in by_position else '0' for t in range(k))
    odd_count = 0
    for t, bit in enumerate(bits, 1):
        odd_count += int(bit)
        if t < k:
            require(3**odd_count >= (1 << t), 'early crossing')
    runs = bits.count('01') + int(bits[0] == '1')
    require(odd_count == d and runs == 179, 'base word')
    total_carry = carry(positions)
    inverse = pow(modulus, -1, factor)
    ceiling = (total_carry-2*gap)//modulus
    require(ceiling == 73, 'common upper bound on a failing excess')

    # On each greedily chosen disjoint 1011 factor, allow a choice 1011/1101.
    # The second pattern moves one odd input one position earlier. The last
    # bit and first bit match, both patterns have the same two odd runs, and
    # proper prefix odd counts can only increase. Hence every choice has
    # the same k,d,runs and first crosses at k. At overlapping occurrences
    # the greedy scan uses only the first occurrence.
    i, weights = 0, []
    while i+4 <= k:
        if bits[i:i+4] == '1011':
            j = by_position[i+2]
            weights.append(3**(d-1-j) * (1 << (i+1)))
            i += 4
        else:
            i += 1
    require(len(weights) == 74, 'number of disjoint choices')
    # Each chosen exchange subtracts its independent weight from C.
    counts = [0]*factor
    counts[total_carry % factor] = 1
    for weight in weights:
        counts = [count + counts[(residue+weight) % factor]
                  for residue, count in enumerate(counts)]
    size = 1 << len(weights)
    require(sum(counts) == size, 'mass conservation')
    # Residues of C which permit 0<=h<=73, using the common bound for the
    # entire family. Individual words may have a lower upper bound.
    possible = sum(counts[(modulus*h) % factor] for h in range(ceiling+1))
    excluded = size-possible
    require(possible == 1_504_650_676_996_140_590_169, 'survivor count')
    require(excluded == 17_384_815_254_482_440_264_615, 'exclusion count')
    expected = Fraction(ceiling+1, factor)
    difference = Fraction(possible, size) - expected
    print(f'Family: k={k}, d={d}, odd runs={runs}, rho={float(Fraction(3**d, modulus)):.9f}')
    print(f'Paths={size:,}; excluded={excluded:,}; unresolved by this filter={possible:,}')
    print(f'Excluded fraction={float(Fraction(excluded, size)):.12%}')
    print(f'Allowed residue fraction={expected}; actual survivor share minus this fraction={float(difference):.6e}')
    print(f'Base path: required h mod {factor} = {total_carry*inverse%factor}, h ceiling = {ceiling}')


if __name__ == '__main__':
    record_seed()
    many_run_family()
