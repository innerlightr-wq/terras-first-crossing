"""Three finite tests for the Terras first-crossing problem (Python 3.9+).

Run: python3 terras_records_and_sturmian.py

Conventions: a 1 denotes the shortcut step (3x+1)/2. A first-crossing
word has rho_t=3^d_t/2^t >= 1 for t<k, and rho_k<1. A failure has x_k>=n.
The record test locates exponent pairs, not actual visits by a word.
The Hamming test applies to a shortest failure and uses the published
n<k^14.3/3 and distinctness of the k input states. The pair example is
a calibration already excluded by a 14-bit prefix. All tests below use
integer or Fraction arithmetic; decimal approximations are not assertions.
"""
from fractions import Fraction as F
from math import gcd, isqrt


def log_bounds(base, terms=120):
    """Positive-series bounds for log(base) from 2*atanh((base-1)/(base+1))."""
    x = F(base-1, base+1)
    value = sum(2*x**(2*i+1)/(2*i+1) for i in range(terms))
    tail = 2*x**(2*terms+1)/((2*terms+1)*(1-x*x))
    return value, value+tail


def record_and_cuts():
    k = 1_412_737_287_532_845_639_425
    d = 891_337_988_684_307_508_298
    f = 929
    assert all(f % x for x in range(2, isqrt(f)+1))
    # 928=2^5*29. Therefore 3 has exact order 928 and 2=3^810 mod 929.
    assert pow(3, 928, f) == 1
    assert all(pow(3, 928//p, f) != 1 for p in (2, 29))
    assert pow(3, 810, f) == 2
    assert gcd(k, d) == 1 and (d-810*k) % 928 == 0
    assert pow(2, k, f) == pow(3, d, f) == 608

    lo2, hi2 = log_bounds(2)
    lo3, hi3 = log_bounds(3)
    low = k*lo2-d*hi3
    high = k*hi2-d*lo3
    dg_low, dg_high = d*low/hi2, d*high/lo2
    assert F('130.5') < dg_low < dg_high < F('130.6')
    B = 28*10**18
    assert F('2.84228') < B*low < B*high/(1-high) < F('2.84229')
    # If 0<a<d and 0<j-alpha*a<=g=k-alpha*d, the determinant
    # d*j-k*a = d*(j-alpha*a)-a*g lies in (-928,928).
    # It is 0 mod 928 on this lattice, and coprimality precludes zero.
    # Therefore (k,d) is a strict one-sided record.
    #
    # At a proper prefix put s=alpha*a-j>0. Weight above 9/10 means
    # 9*(2^s-1)<1-2^(-g), hence s<g/9 and D=k*a-d*j < 10*d*g/9.
    max_D = (10*dg_high/9).__floor__()
    assert max_D == 145
    inverse = pow(k, -1, d)
    cut_ids = []
    # 1-exp(-u) lies between u-u^2/2 and u for u>0.
    epsilon_low, epsilon_high = low-high*high/2, high
    for det in range(1, max_D+1):
        a = det*inverse % d
        j = (k*a-det)//d
        assert 0<a<d and 0<j<k and k*a-d*j==det
        v_low, v_high = a*lo3-j*hi2, a*hi3-j*lo2
        if v_high < 0:
            continue
        assert 0 < v_low <= v_high < 1
        # v <= exp(v)-1 <= v/(1-v) for 0<v<1.
        delta_low, delta_high = v_low, v_high/(1-v_high)
        if 9*delta_high < epsilon_low:
            cut_ids.append(det)
        else:
            assert 9*delta_low > epsilon_high
    assert cut_ids == [16,23,39,46,53,69,76,83,99,106,113,129,136]
    print('Factor-929 shell (k, d):', k, d)
    print('Possible interior cuts with weight > 90% (determinants):', cut_ids)
    a = 23*inverse % d
    print('Cut D=23 (j, odd count):', (k*a-23)//d, a)


def sturmian_separation():
    """Necessary Hamming separation for all k at or above the length floor.

    Relative to any length-k factor v of a critical Sturmian word, at most
    e/2 mismatches add ones: w has floor(beta*k) ones, v has that or one
    more. Thus rho_t < 3^(1+e/2). Since M_t <= k/3, each input x_t is
    below 2 * 3^(e/2) * k^14.3. At L=ceil(log2 of that bound), distinct
    inputs have different length-L windows. At most L+1 untouched windows
    can occur, while each mismatch spoils at most L. Thus k<=(e+2)L
    (also immediate when L>k).
    """
    K, E = 10_439_860_591, 114_474
    target = 2**10 * K**143 * 3**(5*E)
    L = ((target-1).bit_length()+9)//10
    assert L == 91_196 and (E+2)*L < K
    # For k>=K and e<=E, L < 91196.5 + 28.6*(k/K-1).
    # This follows from the half-bit bound at K and log2(x)<=2*(x-1).
    assert target < 2**911_955
    assert (E+2)*F('91196.5') < K
    assert (E+2)*F('28.6') < K
    print('Minimum mismatches certified for every shortest failure:', E+1)


def signed_pair():
    # The compared depths need the first 60 bits. Compute the full word's
    # Q to exhibit the interval covered on this existing calibration family.
    word = FULL_WORD
    path = [(1, 1, 0)]
    for bit in word:
        p, q, c = path[-1]
        if bit == '1':
            p, c = 3*p, 3*c+q
        path.append((p, 2*q, c))
    p, q, c = path[-1]
    assert len(word)==485 and word.count('1')==306 and 1+word.count('01')==120
    assert all(a>=b for a,b,_ in path[:-1]) and q>p
    Q = F(c, q-p)
    assert 7335 < Q < 7336
    def gap(n):
        a,b,c0 = path[5]
        p0,q0,c1 = path[51]
        return (a*n+c0)/b - (p0*n+c1)/q0
    assert word[5:13] == word[51:59] and word[13] != word[59]
    assert [i-word[:i].rfind('1') for i in (5, 51)] == [2, 3]
    # Equal forward bits imply difference=0 mod 256. Incoming odd gaps
    # imply the signed difference = 2^(-2)-2^(-3) = 2 mod 3.
    assert 512 % 256 == 0 and 512 % 3 == (pow(2,-2,3)-pow(2,-3,3)) % 3
    assert -6 < gap(F(2)) < -5 and 453 < gap(Q) < 454
    n = F(c-2*q, q-p)
    assert 330 < gap(n) < 331 and n.denominator>1
    p14,q14,c14=path[14]
    assert (-c14*pow(p14,-1,q14))%q14 == 15087 > Q
    print('Signed pair (5,51): -6 < gap < 454 on [2,Q], gap = 512 mod 768 required.')


FULL_WORD = (
    "1111011110101011011100011111011010101011100111011001111010111110110001"
    "0111110101001110011101011011110010110110111001110011010111100111011100"
    "0011111001111110001111111111101101101111101111101111110111111111011111"
    "1111111011110101111111110111011111101101110010111100110100111101001100"
    "1111110101101101101011011010101110110100111110011011011010110110100111"
    "1011110010101101101011011011101010110111010110101101001111100110001101"
    "10110000000100000000010001100010000110100010100100010101101100100"
)

if __name__ == '__main__':
    record_and_cuts()
    sturmian_separation()
    signed_pair()
