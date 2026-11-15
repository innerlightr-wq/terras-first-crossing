"""A compressed first-crossing family on the factor-929 record shell.
Run with Python 3, using standard library only. No giant words or powers built.

At every multiple of eight t between 8 and K-9 prescribe S(t)=floor(tD/K)+2.
The first word is 11111110, the last 111100000. Within a middle byte, take
c=S(t+8)-S(t) ones. The control has two choices; the larger family has all
bytes ending in zero with partial discrepancy >= -0.8 against slope .631.
The bounds in this file prove statements uniformly over every independent
choice. Passing an endpoint budget or a lower bound on a prefix realizer
is a relaxation of integer realizability. Exclusion by distinct states
uses minimal crossing length among all failures.
"""
from fractions import Fraction as F
from itertools import product

K = 1_412_737_287_532_845_639_425
D = 891_337_988_684_307_508_298
B = 28 * 10**18
END = K - 9
HEAD, TAIL = "11111110", "111100000"
N = (END-8)//8


def log_bounds(base):
    x = F(base-1, base+1)
    lo = sum(2*x**(2*i+1)/(2*i+1) for i in range(120))
    return lo, lo+2*x**241/(241*(1-x*x))


def boundary_count(t):
    return t*D//K + 2


def alphabet(c):
    choices = []
    for bits in product("01", repeat=8):
        w = "".join(bits)
        if w[-1] != "0" or w.count("1") != c:
            continue
        if all(1000*w[:r].count("1")-631*r >= -800 for r in range(9)):
            choices.append(w)
    return choices


def main():
    lo2, hi2 = log_bounds(2)
    lo3, hi3 = log_bounds(3)
    lam_lo, lam_hi = K*lo2-D*hi3, K*hi2-D*lo3
    beta_lo, beta_hi = lo2/hi3, hi2/lo3
    assert F(".630") < beta_lo < beta_hi < F(".631")
    assert 0 < lam_lo < lam_hi < 1
    assert lam_hi/lo3 < F(".001")
    # E = exp(Lambda)-1. Lambda <= E <= Lambda/(1-Lambda).
    E_min, E_max = lam_lo, lam_hi/(1-lam_hi)
    assert F("2.84228") < B*E_min < B*E_max < F("2.84229")
    assert pow(2,K,929) == pow(3,D,929) == 608

    assert END == 8+8*N and N == 176_592_160_941_605_704_926
    assert HEAD.count("1") == boundary_count(8) == 7
    assert TAIL.count("1") == D-boundary_count(END) == 4
    assert 5 < F(8*D,K) < 6  # therefore each byte has exactly 5 or 6 ones
    words = {c:alphabet(c) for c in (5,6)}
    assert {c:len(ws) for c,ws in words.items()} == {5:18,6:7}
    for c,ws in words.items():
        assert "1"*c+"0"*(8-c) in ws
        assert "1"*(c-1)+"01"+"0"*(7-c) in ws
        for w in ws:
            for r in range(9):
                count = w[:r].count("1")
                # At boundaries 1-.001 < S(t)-beta*t <= 2.
                assert 1-F(".001")+count-r*beta_hi > F(".19")
                assert 2+count-r*beta_lo < 5
    # Head has positive drift after the first step. Tail first rises and
    # then falls, ending just below zero. Check all nine proper tail starts.
    for r in range(1,9):
        assert F(".19") < HEAD[:r].count("1")-r*beta_hi
        assert HEAD[:r].count("1")-r*beta_lo < 5
    for r in range(9):
        # d - beta*K = -Lambda/log(3)
        remaining_ones = TAIL[r:].count("1")
        assert (9-r)*beta_lo-remaining_ones-lam_hi/lo3 > F(".19")
        assert (9-r)*beta_hi-remaining_ones < 5
    assert F(6,5)**100 < 3**19
    # Thus 6/5 < rho_t < 243 for every noninitial proper prefix.
    # There is no interior run of seven ones: each middle byte ends in 0,
    # contains at most six ones, and the tail starts with four ones. The
    # initial 1111111 occurs once. No initial square has period >=7.
    assert all(w.endswith("0") and w.count("1")<=6
               for ws in words.values() for w in ws)
    c6 = D-HEAD.count("1")-TAIL.count("1")-5*N
    assert 0<c6<N and all("111" in w for w in words[6])
    assert 3*beta_hi < 2  # every critical Sturmian length-3 word has <=2 ones
    # Every byte contains at least one run and ends in 0. Disjoint 111s
    # in the c=6 bytes each force a Hamming mismatch at every intercept.

    # Independent reconstruction of >90% candidate positions. If tau>9/10,
    # s=log_2(rho_j) < g/9 with g=Lambda/log(2).
    # For a odds at j, det=K*a-D*j = a*g + D*s < 10*D*g/9.
    inv = pow(K,-1,D)
    upper = (10*D*lam_hi/(9*lo2)).__floor__()
    assert upper == 145
    strong = []
    for det in range(1,upper+1):
        a = (inv*det)%D
        j = (K*a-det)//D
        v_lo, v_hi = a*lo3-j*hi2, a*hi3-j*lo2
        if v_hi <= 0:
            continue
        assert 0<v_lo<=v_hi<1
        if 9*v_hi/(1-v_hi) < lam_lo-lam_hi**2/2:
            strong.append(det)
            t = 8*(j//8)
            assert 8 <= t < END
            c = boundary_count(t+8)-boundary_count(t)
            assert min(boundary_count(t)+w[:j-t].count("1") for w in words[c]) > a
        else:
            assert 9*v_lo > lam_hi
    assert strong == [16,23,39,46,53,69,76,83,99,106,113,129,136]

    # Each odd step adds 1/(3*rho_t). Thus d/729 < M <= d/3.
    # Take h in 1,...,929 satisfying 2^K h = C_w (mod 929).
    # n_h = [M-(1+E)h]/E is the corresponding relaxed endpoint seed.
    n_lower = (F(D,729)-1858)/E_max
    Q_upper = F(D,3)/E_min
    assert E_max<1 and n_lower>10**37>B+2**100
    assert Q_upper<2**132
    # The first integer representative above B of any length-100 prefix
    # is <= B+2^100. This is a floor test, not an integrality test on n_h.
    stride = 2**7*3**9  # every noninitial run has length <=6
    assert K*stride < n_lower/5
    # Interior values are >6*n_h/5. Even the last of the first k permitted
    # slots in a depth<=9 equal-weight bucket is <=n_h+k*stride. Thus these
    # candidates also clear the individual and aggregate placement ceilings.
    assert 243*(2**132+F(D,3)) < 2**140
    # Control: >4^18 starts but only four possible byte types. For a
    # shortest failure the input states are distinct. A repeated 18-byte
    # window would give spacing >=2^144 among positive states <2^140.
    assert N-17>4**18
    print("Middle bytes:",N,"of which have six ones:",c6)
    print("Strong cuts avoided by both alphabets:",strong)
    print("Every word has a relaxed endpoint seed > 10^37; 100-bit prefix floors clear.")
    print("Control family: repeated 144-bit windows exclude shortest failure.")
    print("Larger family: 18 choices or 7; no full overlap/restart survival claimed.")


if __name__ == "__main__":
    main()
