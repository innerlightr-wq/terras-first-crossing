"""Finite-prefix packing on an explicit family without low interior returns.

Python 3.8+; standard library only. Run: python3 terras_local_packing.py

Let alpha=log_2(3), d>=2, k=ceil(alpha*d), e=k-alpha*d. Put odd inputs
at zero-based positions S_i=max(i,floor(i*alpha)-1), 0<=i<d; fill the
remaining positions through k-1 with even inputs. The S_i increase
strictly. At a nonempty proper prefix of length j carrying r<d odd steps,
1<=r and j<=S_r. Thus alpha*r-j >= min((alpha-1)*r,1+{alpha*r})
>=alpha-1. When r=d the smallest possible drift is alpha*d-(k-1)=1-e.
Consequently when e<=2-alpha every interior coefficient is >=3/2 (exactly
3/2 at the first step), and the first drop below 1 occurs at k. Each
memory summand 2^(S_i-alpha*i)/3 exceeds 1/12, hence M>d/12. Along
upper convergents e->0 these words have arbitrarily shallow endpoints.
Their number of odd runs tends to infinity: apart from a bounded initial
piece, successive gaps are 1 or 2, and their mean tends to alpha>1.
Moving one odd input earlier by an allowed 1011->1101 swap cannot reduce an
interior coefficient; an unchanged memory term gives the original strict
bound, while a changed term loses at most 1/6.

This checker instantiates the family on k=485, d=306, with one exchange at
zero-based position 312. It omits seed integrality and the known global
seed/length floors deliberately, to compare necessary conditions at the
same rational relaxed candidate. It is a calibration, not a Terras case
left open by current bounds.

Exact finite calibrations for shortest-failure restart inequalities.

Run with standard Python 3.8+. All comparisons use integers or Fractions.
The words here are below the published length floor; this compares tests.

For interior input i define e(i) as the first later shortcut index
j with rho_j < rho_i. This exists by k. Shortestness of the original
hypothetical failure implies x_i > x_{e(i)}. For L equal initial parities
at i and e(i), the parity-word congruence implies
x_i-x_{e(i)} >= 2^L. One can also sum these lower bounds along a chain of
restart edges. Writing x_i=rho_i*n+b_i, this gives a seed lower bound
(b_{e(i)}-b_i+2^L)/(rho_i-rho_{e(i)}).

"""

from collections import Counter
from fractions import Fraction as F
from math import prod


def require(condition, description):
    if not condition:
        raise ArithmeticError(description)


def run_product_upper(word, n, depth):
    """Group maximal odd runs by length and short history; pack strictly >n.

    In a shortest failure the odd inputs before the first crossing are
    distinct (an earlier repeat supplies a shorter failure at the cycle
    minimum). All noninitial odd inputs exceed n. For a run of length L
    from odd x, multiplying its L factors yields 1+(1-(2/3)^L)/x and
    x+1 has v_2=L. The preceding t=min(depth,i) gaps a between the odd
    positions give x mod 3^t from 2^a y_next=3*y+1: the 3^t coefficient
    of the starting y vanishes. Distinct starts in a fixed joint residue
    class are bounded below by its successive representatives above n.
    Different buckets can overlap, which only enlarges the upper bound.
    """
    pos = [i for i, b in enumerate(word) if b == '1']
    gaps = [b-a for a,b in zip(pos,pos[1:])]
    groups = Counter()
    for i,start in enumerate(pos):
        if start and word[start-1] == '1':
            continue
        length = word.index('0',start)-start
        weight = 1-F(2,3)**length
        if not start:
            upper = 1+weight/n
            continue
        used = min(depth,i)
        modulus3 = 3**used
        residue3 = 0
        for gap in gaps[i-used:i]:
            residue3 = (3*residue3+1)*pow(2**gap,-1,modulus3) % modulus3
        modulus2 = 2**(length+1)
        residue2 = 2**length-1
        offset = (residue3-residue2)*pow(modulus2,-1,modulus3) % modulus3
        groups[residue2+offset*modulus2,modulus2*modulus3,weight] += 1
    for (residue,modulus,weight),size in groups.items():
        first = residue + modulus*((n-residue)//modulus+1)
        upper *= prod(1+weight/(first+t*modulus) for t in range(size))
    return upper


def check():
    k,d,site = 485,306,312
    positions = {max(i,(3**i).bit_length()-2) for i in range(d)}
    base = ''.join('1' if i in positions else '0' for i in range(k))
    require(base[site:site+4] == '1011', 'swap site')
    word = base[:site]+'1101'+base[site+4:]
    require(1+word.count('01')==178 and word.count('1')==306, 'run/odd counts')
    p,q,c = 1,1,0
    prefixes = []
    for bit in word:
        if bit == '1':
            p,c = 3*p,3*c+q
        q *= 2
        prefixes.append((p,q,c))
    P,Q,C = prefixes.pop()
    gap = Q-P
    require(min(F(p,q) for p,q,c in prefixes)==F(3,2) and gap>0, 'first crossing')
    modulus = 929
    require(gap%modulus==0,'fixed factor')
    h0 = C*pow(Q,-1,modulus)%modulus
    require((C-2*gap)//Q == 37 and h0==21, 'sole eligible residue')
    n = F(C-Q*h0,gap)
    require(n>2 and n.denominator>1, 'positive relaxed seed is nonintegral')
    cuts = [F((p-q)*C+gap*c,p*Q-q*P) for p,q,c in prefixes]
    require(len(cuts)==484 and min(cuts)>F('36.9517144') and
            min(cuts)<F('36.9517145') and cuts.index(min(cuts))==0, 'all interior cuts')
    needed = F(Q,P)*(1+h0/n)
    bounds = [run_product_upper(word,n,depth) for depth in range(1,6)]
    require(all(y>needed for y in bounds[:-1]) and bounds[-1]<needed, 'depth 5 isolates')
    require(F('1.0023629420')<needed<F('1.0023629421'), 'product required')
    require(F('1.0034051858')<bounds[3]<F('1.0034051859'), 'depth 4')
    require(F('1.0022071267')<bounds[4]<F('1.0022071268'), 'depth 5')
    print(f'k={k}, d={d}, 178 runs, min proper rho=3/2, M={float(F(C,P)):.12f}')
    print(f'h={h0}, min U={float(min(cuts)):.12f}, required={float(needed):.12f}')
    print(f'At floor n: depth 4 <= {float(bounds[3]):.12f}, depth 5 <= {float(bounds[4]):.12f}')
    print('All assertions use exact rational arithmetic. Decimals are for display.')



def check_translated_crossings():
    """An irregular calibration passing all 484 translated first-crossing cuts.

    If r is the first smaller coefficient after interior time i, shortestness
    forces x_i>x_r. If the starts share the next L parity bits, their
    difference is a multiple of 2^L; combined their difference is >=2^L.
    At r=k the common observed suffix is empty. The stronger unsigned test
    on arbitrary pairs is checked separately and already excludes this word.
    """
    word = (
        '111110010111101101010110101111010101010111011010101101111101111100111011010'
        '011011011011100111000011011111001001101101011100111011111001111101101101101'
        '010011011011110110010110111010001111101000101110011011101110100110111010101'
        '110111011111011010001101101100111000111110010111011111100011111011011001001'
        '110001101111011001110011111110100101001011011010111011101110100111111101011'
        '001111000111110010011011110000111011110010101011110111101001110101111110000'
        '01111100011111110110001011011011000'
    )
    k = len(word)
    ps, qs, carries = [], [], []
    p = q = 1
    c = 0
    for bit in word:
        ps.append(p)
        qs.append(q)
        carries.append(c)
        if bit == '1':
            p, c = 3*p, 3*c+q
        q *= 2
    ps.append(p)
    qs.append(q)
    carries.append(c)
    rho = [F(pi,qi) for pi,qi in zip(ps,qs)]
    require(1+word.count('01')==123 and word.count('1')==306 and k==485,
            'irregular shell')
    require(min(rho[1:k])==F(729,512) and p<q, 'all proper drifts')
    gap = q-p
    h = c*pow(q,-1,929)%929
    require(gap%929==0 and h==20 and (c-2*gap)//q==22, 'factor and cap')
    n = F(c-q*h,gap)
    xs = [(pi*n+ci)/qi for pi,qi,ci in zip(ps,qs,carries)]
    require(xs[-1]==n+h and min(xs[1:k])>xs[-1] and n.denominator>1,
            'relaxed endpoint and all interior minima')
    cuts = []
    for i in range(1,k):
        r = next(r for r in range(i+1,k+1) if rho[r] < rho[i])
        length = 0
        while r+length < k and word[i+length] == word[r+length]:
            length += 1
        cuts.append((xs[i]-xs[r]-2**length,i,r,length))
    margin,i,r,length = min(cuts)
    require(len(cuts)==484 and (i,r,length)==(53,118,3) and margin>37,
            'all translated first crossing/LCP tests')
    target = (1+h/n)*q/p
    upper = [run_product_upper(word,n,depth) for depth in range(1,8)]
    require(all(y>target for y in upper[:-1]) and upper[-1]<target,
            'first packing exclusion at depth 7')
    require(word[:12]==word[261:273] and abs(xs[0]-xs[261]) < 2**12,
            'stronger arbitrary-pair test already excludes this example')
    print(f'Irregular example: 123 runs, h=20, required={float(target):.12f}, '
          f'depth 6={float(upper[5]):.12f}, depth 7={float(upper[6]):.12f}.')
    print('All translated first-crossing/LCP tests allow this relaxed candidate; '
          'the arbitrary-pair LCP test at indices 0 and 261 excludes it.')

def word(k, d, shift=0, swap=None):
    sites = {max(i, (3**i).bit_length() - 1 - shift) for i in range(d)}
    w = ''.join('1' if i in sites else '0' for i in range(k))
    if swap is not None:
        assert w[swap:swap+4]=='1011'
        w = w[:swap]+'1101'+w[swap+4:]
    return w


def affine_prefixes(w):
    rho, b = F(1), F(0)
    path=[(rho,b)]
    for bit in w:
        if bit == '1':
            rho, b = F(3,2)*rho, (3*b+1)/2
        else:
            rho, b = rho/2, b/2
        path.append((rho,b))
    assert all(r>=1 for r,_ in path[:-1]) and rho<1
    return path


def first_local_crossing(path,i):
    return next(j for j in range(i+1,len(path)) if path[j][0]<path[i][0])


def bounds(w):
    path=affine_prefixes(w)
    rho,b=path[-1]
    k=len(w)
    denom=2**k
    gap=denom-int(rho*denom)
    assert gap%929==0 and (b*denom).denominator==1
    h=(int(b*denom)*pow(denom,-1,929))%929
    seed=(b-h)/(1-rho)
    values=[r*seed+z for r,z in path]
    assert seed.denominator>1  # Deliberately relax the full seed congruence.
    cap=(int(b*denom)-2*gap)//denom
    assert h<=cap<h+929
    return path,h,seed,values,cap


def finite_restart():
    w=word(1455,918)
    path,h,seed,y,cap=bounds(w)
    assert (h,cap)==(13,220)
    endpoint_gap=min(v-y[-1] for v in y[1:-1])
    assert 141<endpoint_gap<142  # All interior-vs-endpoint cuts are satisfied.
    i,j=84,569
    assert first_local_crossing(path,i)==j and -5<y[i]-y[j]<-4
    ri,bi=path[i]; rj,bj=path[j]; r,b=path[-1]
    h_cut=b-(1-r)*(bj-bi)/(ri-rj)
    assert h_cut==F(1,6)
    print('Fresh 1455-step mechanical word, h=13: all endpoint tests pass;')
    print('restart 84->569 gives upper cutoff', float(h_cut),
          'and x84-x569=',float(y[i]-y[j]))


def parity_overlap():
    w=word(485,306,shift=1,swap=312)
    path,h,seed,y,cap=bounds(w)
    assert (h,cap)==(21,37)
    assert min(v-y[-1] for v in y[1:-1]) > 0
    assert [i for i in range(1,len(w)//2+1) if w[:i]==w[i:2*i]]==[1,2]
    i,j=354,419
    assert first_local_crossing(path,i)==j and w[i:i+65]==w[j:j+65]
    ri,bi=path[i]; rj,bj=path[j]; r,b=path[-1]
    Q=b/(1-r)
    low=(ri-rj)*2+bi-bj
    high=(ri-rj)*Q+bi-bj
    assert -10<low<high<836 and Q<36240
    exact_root=(bj-bi)/(ri-rj)
    assert 396<exact_root<397
    # 2^65 divides x_i-x_j, so no integer seed in [2,Q] can realize w.
    n_floor=(bj-bi+2**65)/(ri-rj)
    assert n_floor>F(15834*10**17)
    print('Shifted 485-step word, h=21: 65-bit overlap; seed equality would')
    print('require n=',float(exact_root),'; under shortestness n>=',float(n_floor))

def check_shell_arithmetic():
    """Previous proposed 929 calibrations used incorrect shortcut-step counts."""
    shells = [(485, 306), (3038, 1917), (5591, 3528), (8144, 5139),
              (10697, 6750), (13250, 8361)]
    for i, (k, d) in enumerate(shells):
        assert ((1 << k) > 3**d) == (i == 0)
        assert (pow(2, k, 929) == pow(3, d, 929)) == (i == 0)
    for b in range(1, 7):
        k, d = 485*b, 306*b
        assert (1 << (k-1)) <= 3**d < (1 << k)
        assert pow(2,k,929) == pow(3,d,929)


if __name__ == '__main__':
    check()
    check_translated_crossings()
    finite_restart()
    parity_overlap()
    check_shell_arithmetic()
