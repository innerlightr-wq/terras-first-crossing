"""Finite and infinite tests of prefix memory at a first coefficient crossing.
Run python3 terras_prefix_families.py (Python 3.9+, standard library).
All inequalities used by the verification use integers or rationals.

Infinite construction. Write alpha=log_2(3). For a bit word w use (rho,M)
defined by T^k(x)=rho*(x+M); M(vw)=M(v)+M(w)/rho(v). For a shortest
hypothetical failure the endpoint lies strictly below each interior value.
Indeed the smallest interior value has a shorter coefficient crossing by the
endpoint; it must descend there since no interior value is smaller. Writing
epsilon=1-rho(w), delta=rho(prefix)-1, the integer overshoot h satisfies

  h < U=(delta*rho(w)*M(w)+epsilon*rho(prefix)*M(prefix))/(delta+epsilon).

1. Choose an unbounded sequence of integers u by {2*alpha*u}->1 from below,
possible by irrational rotation. D=928u, K=464 ceil(2*alpha*u). Then
eta=K-alpha*D ->0+, so eventually K=ceil(alpha*D). Use the length-K suffix
with odd sites floor(i*alpha), i=0,...,D-1. Its proper-prefix coefficients
are at least 1 (up to the next odd input, t<=floor(i*alpha); after the last
odd input, t<=floor(D*alpha)). Its endpoint rho_s=2^-eta is below one.
The four zero-based sites 9,28,50,66 each allow a disjoint 1011->1101 swap.
Each independently moves one odd input one place earlier; all proper-prefix
coefficients remain >=1. M_s> D/6-2/3 since every original normalized
summand is 2^-{i*alpha}/3 >1/6, and at most four are halved.

2. Choose t as large as necessary with {2*alpha*t}->0+, a=928t and
j=464 floor(2*alpha*t), so 0<delta=3^a/2^j-1<(1-rho_s)/(100(D+1)).
Take the prefix 11110 followed by H=2a-j-3 copies of 110 and
L=2j-3a+2 copies of 10. Both counts are positive eventually. The counts
of odd inputs and steps in the prefix are a and j. An adjacent pair of high
blocks may instead be 111010, and a pair of low blocks may be 1100; lengths
and odd counts are preserved. The high blocks have rising coefficients
from boundary to boundary. In the low phase each block reaches its minimum
at its endpoint, and the last boundary still has coefficient 1+delta.
The combined word therefore first crosses at j+K: throughout the proper
suffix prefixes its coefficient is >=1+delta, and (1+delta)*rho_s<1.

3. The marker's M=65/81, rho=81/32. For any high block or high pair,
M<=5*(1-1/rho), with equality for 110 and its ordinary pair, so summing
over their starts gives prefix memory <65/81+(32/81)*5=25/9. In the low
phase each block has M<=1/rho-1, so telescoping gives an additional
amount <=1/(1+delta)-1/rho_high<1. Thus M(prefix)<34/9. Because
rho(prefix)>1 and M_s<=D/3, M(total)<34/9+D/3. With the above choice of
delta, epsilon=(1-rho_s)-delta*rho_s > 0.99*(1-rho_s).
It follows that U < (34/9+D/3)/(99*(D+1))+(1+delta)*34/9 < 4 for D>=1.
For example the first term is <=(34/9)/99, delta<1/100. Meanwhile
M(total)>=M_s/(1+delta) tends to infinity and rho(total) tends to 1.

4. Modulo 929, 2^464=3^928=1. Therefore 929 divides 2^(j+K)-3^(a+D).
For total carry C, a failure must satisfy h=C/2^(j+K) modulo 929.
Each of the four swaps changes C by a nonzero product of powers of 2 and 3
independently of the other swaps. Each nonzero translate enlarges the union
of a nonempty proper subset and its translate in the prime group Z/929.
Four switches therefore give at least five residues. At most four can be
0,1,2,3, so at least one of the sixteen words excludes every candidate h.
This is an existence statement for each such group, not a common choice
of switches for different prefixes. It excludes occurrences as SHORTEST
failures. The initial run of four ones is unique: each other odd run has
length at most three. In the suffix the original ones runs have length at
most two, since 2*alpha>3; the four isolated swaps increase this locally
to at most three. Thus no initial square can have period >=4. This does
not test all other realizer obstructions; no priority claim.

This program checks the elementary building blocks, a finite illustration
of the infinite proof and the separate finite family on the 485-step shell.
Its finite tests alone do not prove the asymptotic statement.
"""
from fractions import Fraction
from itertools import product
from math import isqrt



def affine(bits):
    c, d = 0, 0
    for j, bit in enumerate(bits):
        if bit == "1":
            c = 3*c + 2**j
            d += 1
    return c, d


def verify():
    k, d, j, factor = 485, 306, 84, 929
    positions = [(3**a).bit_length()-1 for a in range(d)]
    index = {pos: a for a, pos in enumerate(positions)}
    word = "".join("1" if t in index else "0" for t in range(k))
    C, odd = affine(word)
    Cj, prefix_odd = affine(word[:j])
    P, Q = 3**d, 2**k
    p, q = 3**prefix_odd, 2**j
    A = Q-P
    assert odd == d and A % factor == 0
    running = 0
    for t, symbol in enumerate(word, 1):
        running += int(symbol)
        assert (3**running < 2**t) == (t == k)
    weights = []
    fixed = []
    t = 0
    while t+4 <= k:
        if word[t:t+4] != "1011":
            t += 1
            continue
        # The only mixed cut is after the second symbol. Freezing that site
        # fixes p at the cut. Moving any other chosen odd input earlier never
        # increases either carry or decreases a proper-prefix coefficient.
        if t+2 == j:
            fixed.append(t)
        else:
            a = index[t+2]
            weights.append(3**(d-1-a) * 2**(t+1))
        t += 4
    assert fixed == [82] and len(weights) == 73
    upper = Fraction((p-q)*C + A*Cj, p*Q-q*P)
    min_cap = (C-sum(weights)-2*A)//Q
    assert Fraction(536675,10000)<upper<Fraction(536676,10000) and min_cap == 62
    inverse = pow(Q, -1, factor)
    counts = [0]*factor
    counts[C*inverse % factor] = 1
    for weight in weights:
        counts = [value + counts[(r+weight*inverse) % factor] for r,value in enumerate(counts)]
    assert sum(counts) == 2**73
    total = sum(counts[54:63])
    assert total == 91499027655170685801
    print(f"{total:,} words pass the factor and their endpoint cap but fail the j={j} comparison.")
    print(f"Uniform U_84 <= {float(upper):.12f}; uniform endpoint cap >= {min_cap}")



def log_interval(value):
    """Positive series for ln(value), with a geometric tail bound."""
    z = Fraction(value - 1, value + 1)
    terms = 130
    low = 2 * sum(z**(2*i+1)/(2*i+1) for i in range(terms))
    err = 2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
    return low, low+err


def block(bits, prime):
    carry, d = affine(bits)
    return pow(3, d, prime), carry % prime, pow(2, len(bits), prime)


def join(v, w, p):
    """For concatenation v followed by w return (3^d,C,2^k) mod p."""
    return w[0]*v[0] % p, (w[0]*v[1]+v[2]*w[1]) % p, v[2]*w[2] % p


def repeat(v, count, p):
    value = (1, 0, 1)
    while count:
        if count & 1:
            value = join(value, v, p)
        v = join(v, v, p)
        count //= 2
    return value


def verify_two_scales():
    for bits, rho, memory in (
        ('11110', Fraction(81,32), Fraction(65,81)),
        ('110', Fraction(9,8), Fraction(5,9)),
        ('110110', Fraction(81,64), Fraction(85,81)),
        ('111010', Fraction(81,64), Fraction(73,81)),
        ('10', Fraction(3,4), Fraction(1,3)),
        ('1010', Fraction(9,16), Fraction(7,9)),
        ('1100', Fraction(9,16), Fraction(5,9))):
        carry, odd = affine(bits)
        assert Fraction(3**odd,2**len(bits)) == rho
        assert Fraction(carry,3**odd) == memory

    p = 929
    assert all(p%i for i in range(2,isqrt(p)+1))
    assert pow(2,464,p)==pow(3,928,p)==1
    ln2, ln3 = log_interval(2), log_interval(3)
    lower, upper = ln3[0]/ln2[1], ln3[1]/ln2[0]
    check_value = 28*10**18
    P, Q = 630118245525664765, 397560349370386783
    assert check_value*(P*ln2[0]-Q*ln3[1]) > Fraction('4.24722566869')
    assert check_value*(P*ln2[1]-Q*ln3[0]) < Fraction('4.24722566870')
    # First strict positive prefix error at the tolerance 3.82e-18.
    p0,q0 = 216627100404256471,136676483074956903
    P0,Q0 = 206745572560704147,130441933147714940
    tol = Fraction('3.82e-18')
    assert q0*P0-Q0*p0==1 and q0*lower-p0>tol and Q0*upper-P0<0
    assert 0<(q0+Q0)*lower-(p0+P0)<(q0+Q0)*upper-(p0+P0)<tol
    # In this determinant-one basis an integer point of positive denominator
    # and positive error <tol must have BOTH coordinates >=1. Hence its
    # denominator >= q0+Q0, and equality is attained.

    # A finite instance of the proof. The calibration deliberately uses the
    # partial factor test alone; the seed floor can exclude this instance.
    D, K = 141984, 225040
    a, j = 698886859409568, 1107709464410944
    assert K%464==j%464==a%928==D%928==0
    assert 0 < K-D*upper < K-D*lower < 1
    gamma = a*ln3[1]-j*ln2[0]
    lam = K*ln2[0]-D*ln3[1]
    assert 0 < a*ln3[0]-j*ln2[1] <= gamma < 1 and lam>0
    # exp(gamma)-1 <=gamma/(1-gamma); 1-exp(-lam)>=lam/(1+lam).
    assert gamma/(1-gamma) < lam/(1+lam)/(100*(D+1))
    H, L = 2*a-j-3, 2*j-3*a+2
    assert min(H,L)>0 and 5+3*H+2*L==j and 4+2*H+L==a
    prefix=block('11110',p)
    for bits,copies in (('110',H),('10',L)):
        prefix=join(prefix,repeat(block(bits,p),copies,p),p)
    # Obtain every suffix odd position as a certified floor of i*alpha.
    scale=10**38
    lo=lower.numerator*scale//lower.denominator
    hi=-((-upper.numerator*scale)//upper.denominator)
    positions=[i*lo//scale for i in range(D)]
    assert all(i*hi//scale==position for i,position in enumerate(positions))
    idx={pos:i for i,pos in enumerate(positions)}
    w=''.join('1' if i in idx else '0' for i in range(K))
    sites=(9,28,50,66)
    assert all(w[i:i+4]=='1011' for i in sites)
    assert all('1111' not in w[i-3:i]+'1101'+w[i+4:i+7] for i in sites)
    carry=0
    for pos in positions:
        carry=(3*carry+pow(2,pos,p))%p
    total=join(prefix,(pow(3,D,p),carry,pow(2,K,p)),p)
    assert total[0]==total[2]==1
    changes=[-prefix[2]*pow(3,D-1-idx[z+2],p)*pow(2,z+1,p)%p for z in sites]
    assert all(changes)
    residues={(total[1]+sum(v*w for v,w in zip(choice,changes)))%p
              for choice in product((0,1),repeat=4)}
    assert len(residues)>=5 and any(h>=4 for h in residues)
    print(f'Finite two-scale illustration: {len(residues)} residues, '
          f'{sum(h>=4 for h in residues)} distinct residues outside [0,3].')


def require(test, message):
    if not test:
        raise AssertionError(message)


def block_carry_mod7(carry, offset, length):
    """Append length successive odd steps at the given shortcut offset."""
    return (pow(3,length,7)*carry + pow(2,offset,7)*
            (pow(3,length,7)-pow(2,length,7))) % 7

def verify_large_shell():
    p, q = 79_641_170_620_168_673_833, 50_247_984_153_525_417_450
    P, Q = 206_745_572_560_704_147, 130_441_933_147_714_940
    k, d = p + P, q + Q
    two, three = log_interval(2), log_interval(3)
    # Verify the smaller reported absolute-error denominator is on the upper side.
    # At tolerance T=3.82e-18 bits, the first lower denominator hitting T
    # is 267118416222671843. To verify minimality, write every integer
    # lattice point as u*(p0,q0)+v*(P,Q) in this unimodular basis.
    # g0=q0*alpha-p0>0, G=Q*alpha-P<0. If 0 < u*g0+v*G < T < g0,
    # and u*q0+v*Q>0, then u>=1 and v>=1; thus denominator>=q0+Q.
    # For u<=0, v>=1 is forced by denominator positive, contradicting
    # error positive. For v<=0 and u>=1 the error is at least g0.
    p0, q0 = 216627100404256471, 136676483074956903
    A_lo, A_hi = three[0]/two[1], three[1]/two[0]
    T = Fraction(382, 10**20)
    require(q0*P-Q*p0==1, "one-sided unimodular basis")
    require(q0*A_lo-p0>T and Q*A_hi-P<0, "neighbor signs")
    require(0<(q0+Q)*A_lo-(p0+P)<(q0+Q)*A_hi-(p0+P)<T,
            "earliest admissible lower approximation")
    # gamma = log(rho_p); Lambda = -log(rho_k).
    gamma_lo, gamma_hi = q*three[0]-p*two[1], q*three[1]-p*two[0]
    lam_lo, lam_hi = k*two[0]-d*three[1], k*two[1]-d*three[0]
    require(0 < gamma_lo < gamma_hi < Fraction(7,10**21), "proper near return")
    require(Fraction(1786,10**21) < lam_lo < lam_hi < Fraction(1787,10**21), "crossing gap")
    require(50 < 28*10**18*lam_lo < 28*10**18*lam_hi < 51,"scale of candidate correction")
    # For 0<x<1, e^x <= 1/(1-x), and for positive y, 1-e^-y >= y/(1+y).
    # Thus E/D >= [Lambda/(1+Lambda)] [(1-gamma)/gamma] where
    # E=1-rho_k and D=rho_p-1.
    ratio_lo = lam_lo/(1+lam_lo) * (1-gamma_hi)/gamma_hi
    require(ratio_lo > 265,"return weight")
    # Prefix 1^t 0 1^(q-t) 0^(p-q-1), t=2 mod6, 2<=t<q.
    # Its early even takes the coefficient to at least (3/2)^2/2=9/8.
    # After its second odd run the evens end at rho_p>1.
    # On the suffix of length P, place the odd steps at floor(i log_2(3))
    # for i<300 and then Q-300 consecutive odd steps starting at floor(300 log_2 3).
    # All remaining steps are even. Before the long run, at the next
    # odd position the coefficient is 3^i / 2^floor(i log_2(3)) >=1.
    # End of the entire word has rho_k>1/2 so the last even crosses first.
    N = 300
    positions = [(3**i).bit_length()-1 for i in range(N+1)]
    require(P > positions[-1]+(Q-N), "final even tail")
    require(lam_hi < two[0], "penultimate coefficient exceeds 1")
    # M_prefix=1+(2/3)^t-2(2/3)^q <=13/9.
    # 300 singleton contributions <=1/3 and the final run contributes <=1
    # in suffix coordinates. The prefix coefficient >1.
    U = Fraction(103,266) + Fraction(13,9)/(1-gamma_hi)
    require(U<Fraction(184,100),"uniform interior cutoff U<1.84")
    # Independent check of an elementary reason these words already descend:
    # if t>=100, 2^t|(n+1). If t<100, the second run starts at
    # y=[(3/2)^t(n+1)-1]/2, giving 2^(q-t)<=y+1 <=(3/2)^t(n+1).
    # Thus n+1 >= 2^(q-t)/(3/2)^t >= 2^(q-2t) >= 2^100.
    require(q-2*100>100 and (2**100-1)*lam_lo>103,
            "long run already exceeds the ordinary descent threshold")
    require(pow(2,3,7)==pow(3,6,7)==1,"orders divide 6")
    require(pow(2,k,7)==pow(3,d,7),"7 divides total gap")
    t=2
    # t=2 mod 6 gives the same answers for all lengths/offsets below.
    carry=block_carry_mod7(0,0,t)
    carry=block_carry_mod7(carry,t+1,q-t)
    for offset in positions[:-1]:
        carry=block_carry_mod7(carry, p+offset, 1)
    carry=block_carry_mod7(carry,p+positions[-1],Q-N)
    h=carry*pow(pow(2,k,7),-1,7)%7
    require(h==3,"h mod 7 must be 3")
    total=(q-3)//6+1
    require(total == 8_374_664_025_587_569_575,"distinct t-values")
    odd_runs_in_suffix=1+sum(b>a+1 for a,b in zip(positions,positions[1:]))
    require(odd_runs_in_suffix+2==178,"total odd runs")
    # Whole-word memory interval (t=2). Isolate the known finite terms;
    # the final run stops within (2/3)^(Q−N) of the infinite geometric tail.
    local=sum(Fraction(2**offset,3**i*3) for i,offset in enumerate(positions[:-1]))
    final_weight=Fraction(2**positions[-1],3**N)
    tiny=Fraction(2,3)**150
    suffix_lo=local+final_weight*(1-tiny)
    suffix_hi=local+final_weight
    Mlo=Fraction(13,9)-2*tiny+(1-gamma_hi)*suffix_lo
    Mhi=Fraction(13,9)+suffix_hi
    require(Fraction(74370,1000)<Mlo<Mhi<Fraction(74371,1000),"74.370 < M < 74.371")
    print(f'{total:,} words; k={k:,}; d={d:,}; 178 runs; epsilon/delta >265.')
    print('At t=2, 74.370 < M < 74.371. Across the family h<2 and h=3 (mod 7).')


if __name__ == '__main__':
    verify()
    verify_two_scales()
    verify_large_shell()
