"""Finite Terras comparisons with exact rational certificates (Python 3.9+).

Run python3 terras_overlap_and_history.py. The test for two distinct input
values is the SIZE inequality |x_i-x_j| >= 2^L where L is their entire
common forward parity prefix within the word. We extend this comparison to
the endpoint here since h>0 and every interior value exceeds the endpoint.
A real orbit would obey the full integer congruence. Each listed candidate
has a rational nonintegral seed. Testing the size inequality separately
shows exactly which relaxed filters imply one another.

Minimality in crossing length would force distinct inputs and signed
descent at the first later coefficient below that of each interior start.
At every noninitial maximal odd run start x the last t incoming odd-step
gaps fix x mod 3^t. The outgoing run length ell fixes
x = 2^ell - 1 mod 2^(ell+1). In a joint residue class, distinct integers
above the seed occupy different terms of its arithmetic progression.
The run contributes 1+(1-(2/3)^ell)/x to the correction product, decreasing
with x. Reuse means allowing every member to occupy the first term.

If two starts in a bucket have identical incoming gap sequences, their
predecessors t odd steps back share S+ell+1 parity bits (S=sum gaps).
Their separation by at least 2^(S+ell+1), transported by the common slope
3^t/2^S, forces current separation by the full stride 3^t*2^(ell+1).
If every bucket has a single gap pattern and every member also clears the
first representative, all-pairs spacing already implies the packed ceiling.
"""
from collections import defaultdict
from fractions import Fraction as F
from math import prod


ORIGINAL_SHELL = (
    "1111011110101011011100011111011010101011100111011001111010111110110001"
    "0111110101001110011101011011110010110110111001110011010111100111011100"
    "0011111001111110001111111111101101101111101111101111110111111111011111"
    "1111111011110101111111110111011111101101110010111100110100111101001100"
    "1111110101101101101011011010101110110100111110011011011010110110100111"
    "1011110010101101101011011011101010110111010110101101001111100110001101"
    "10110000000100000000010001100010000110100010100100010101101100100"
)
SMALLER_SHELL = (
    "1111111001110100101011011111010101100111101001110101111110110011101101"
    "1010011111111111001011111001101010111111111111101011100001111111101101"
    "1100101100111010011001111101101110100101111010011001111100111001101011"
    "1101110100110111011001001011111110101001101100111110010011111100010001"
    "0110011101101011111111111100110011000010101001111101011111001011110010"
    "001110001001010011110111001010010110110101011110000"
)


def affine_prefixes(word):
    """Return (3^d_i, 2^i, C_i) at every prefix, including both endpoints."""
    p, q, c = 1, 1, 0
    path = [(p, q, c)]
    for bit in word:
        if bit == "1":
            p, c = 3*p, 3*c+q
        q *= 2
        path.append((p, q, c))
    return path


def pack(word, seed, values, depth):
    positions = [j for j, bit in enumerate(word) if bit == "1"]
    groups = defaultdict(list)
    for i, j in enumerate(positions):
        if j and word[j-1] == "1":
            continue
        ell = word.index("0", j)-j
        weight = 1-F(2, 3)**ell
        if j == 0:
            packed = reused = 1+weight/seed
            continue
        used = min(i, depth)
        history = tuple(positions[z+1]-positions[z] for z in range(i-used, i))
        mod3 = 3**used
        res3 = 0
        for gap in history:
            res3 = (3*res3+1)*pow(2, -gap, mod3) % mod3
        mod2, res2 = 2**(ell+1), 2**ell-1
        residue = res2 + mod2*((res3-res2)*pow(mod2, -1, mod3) % mod3)
        groups[residue, mod2*mod3, weight].append((j, history))
    floors_missed = []
    mixed = []
    for (residue, stride, weight), members in groups.items():
        first = residue + stride*((seed-residue)//stride + 1)
        count = len(members)
        packed *= prod(1+weight/(first+i*stride) for i in range(count))
        reused *= (1+weight/first)**count
        floors_missed.extend((j, first) for j, _ in members if values[j]<first)
        if len({history for _, history in members})>1:
            mixed.append((residue, stride, [j for j, _ in members]))
    return packed, reused, sorted(floors_missed), mixed


def certify(word, factor, expected_h, expected_cap, expected_runs, depth):
    path = affine_prefixes(word)
    k = len(word)
    p, q, c = path[-1]
    gap = q-p
    assert gap>factor>1 and gap%factor==0 and min(F(a,b) for a,b,_ in path[:-1])==1
    assert word.startswith("1") and word.endswith("0") and 1+word.count("01")==expected_runs
    h, cap = c*pow(q,-1,factor)%factor, (c-2*gap)//q
    assert (h,cap)==(expected_h,expected_cap) and 0<h<=cap<h+factor
    seed = F(c-q*h, gap)
    assert seed>=2 and seed.denominator>1
    values = [(a*seed+cc)/b for a,b,cc in path]
    assert values[-1]==seed+h and min(values[1:-1])>values[-1]
    # Express every rational path value with a single denominator to compare
    # pairwise distances by integer cross multiplication.
    denominator = seed.denominator*q
    scaled = [int(x*denominator) for x in values]
    assert all(F(y,denominator)==x for x,y in zip(values,scaled))
    restart = {}
    for i in range(1,k):
        a,b,_ = path[i]
        restart[i] = next(j for j in range(i+1,k+1) if path[j][0]*b < a*path[j][1])
    closest = None
    num_pairs = num_signed = 0
    for shift in range(1,k+1):
        shared = 0
        for i in range(k-shift,-1,-1):
            j = i+shift
            shared = shared+1 if j<k and word[i]==word[j] else 0
            required = denominator << shared
            delta = scaled[i]-scaled[j]
            margin = abs(delta)-required
            assert margin>=0, ("failed overlap", i,j,shared)
            if j<k and (closest is None or margin<closest[0]):
                closest = (margin,i,j,shared)
            if restart.get(i)==j:
                assert delta>=required
                num_signed+=1
            num_pairs+=1
    assert num_pairs == k*(k+1)//2 and num_signed==k-1
    target = (1+F(h)/seed)*F(q,p)
    for t in range(depth):
        assert pack(word,seed,values,t)[0]>=target
    previous = pack(word,seed,values,depth-1)[0]
    packed,reused,missed,mixed = pack(word,seed,values,depth)
    assert packed<target<reused
    assert missed and mixed
    if k==485:
        assert min(F(a,b) for a,b,_ in path[1:-1])==F(3,2)
        assert len(missed)==39 and missed[0]==(12,72505) and 25850<values[12]<25851
        assert closest[1:]==(3,87,1) and closest[0]>47*denominator
        assert (12145,78732,[306,344,368]) in mixed
    print(f"k={k}, d={word.count('1')}, runs={expected_runs}, h={h}, cap={cap}, n={float(seed):.12f}")
    print(f"input pairs={k*(k-1)//2}, min input margin={float(F(closest[0],denominator)):.12f}, restarts={num_signed}")
    print(f"target={float(target):.12f}, depth {depth-1}={float(previous):.12f}")
    print(f"depth {depth}: reuse={float(reused):.12f}, distinct={float(packed):.12f}, individual floors missed={len(missed)}")


if __name__ == "__main__":
    certify("110111110101111000101101010",5,1,2,8,1)
    certify(SMALLER_SHELL,53,3,5,97,7)
    certify(ORIGINAL_SHELL,929,2,7,120,9)
